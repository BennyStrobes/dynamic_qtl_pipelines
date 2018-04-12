args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(reshape)

get_allelic_df <- function(data, index) {
    if (data[index] == "NA") {
        return("NA")
    } else {
        a1 <- as.numeric(strsplit(data[index],';')[[1]])
        a2 <- as.numeric(strsplit(data[(index+1)],';')[[1]])
        environ <- as.numeric(strsplit(data[(index+2)],';')[[1]])
        df <- data.frame(a1=a1, a2=a2, environmental_var=environ, tot = (a1+a2), frac =(a1/(a1+a2)))
        return(df)
    }
}

make_as_plot <- function(allelic_df, site_number) {
    max_val <-max(max(allelic_df$a1), max(allelic_df$a2))

    p <- ggplot(allelic_df, aes(x=a1, y=a2, color=environmental_var)) +  geom_point() +
        theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        labs(x = "allele 1 counts", y = "allele 2 counts", color = "time step",title=paste0("site",site_number)) +
        geom_abline() +
        theme(legend.position="bottom") +
        scale_color_gradient(low="pink",high="blue") +
        scale_x_continuous(limits = c(0, max_val + 1), breaks = round(seq(0, max_val, by = 3),1)) + scale_y_continuous(limits = c(0,max_val+1), breaks = round(seq(0, max_val, by = 3),1)) 
    return(p)
}

make_as_plot_v2 <- function(allelic_df, site_number) {
    p <- ggplot(allelic_df, aes(x=environmental_var, y=frac, size=tot)) +  geom_point() +
        theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        labs(x = "Time step", y = "allele fraction", size = "total counts")
    return(p)
}


make_te_plot <- function(ensamble_id, rs_id, pvalue, beta, te_df, nb_conc) {
    # Change box plot line colors by groups
    p <- ggplot(te_df, aes(x=time_step, y=gene_counts, fill=genotype)) + geom_boxplot() +
       theme(text = element_text(size=7),plot.title = element_text(size=3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
       labs(x = "Time Step", y = "log(counts)", title=paste0(ensamble_id, " / ", rs_id, " / p=", pvalue, " \n beta_full=", strsplit(beta,"/")[[1]][1], "\n beta_null=",strsplit(beta,"/")[[1]][2] , "\n nb_conc_full=",strsplit(nb_conc,"/")[[1]][1], "\n nb_conc_null=",strsplit(nb_conc,"/")[[1]][2])) +
       theme(legend.position="bottom")
    return(p)
}

parameter_string = args[1]
visualization_dir = args[2]
output_dir = args[3]


input_file <- paste0(visualization_dir, parameter_string, "_dynamic_qtl_hits_summary.txt")
print(input_file)
stop = FALSE
f = file(input_file, "r")
while(!stop) {
  next_line = readLines(f, n = 1)
  data = strsplit(next_line,'\t')[[1]]
  ensamble_id <- data[1]
  rs_id <- data[2]
  pvalue <- data[3]

  beta <- data[4]
  environmental_vars <- as.numeric(strsplit(data[5],';')[[1]])
  gene_counts <- as.numeric(strsplit(data[6],';')[[1]])
  genotype <- as.numeric(strsplit(data[7],';')[[1]])
  nb_conc <- data[8]
  # allelic_counts <- get_allelic_counts(data)
  allelic_df_1 <- get_allelic_df(data,9)
  allelic_df_2 <- get_allelic_df(data,12)
  allelic_df_3 <- get_allelic_df(data,15)
  allelic_df_4 <- get_allelic_df(data,18)
  te_df <- data.frame(gene_counts=log(gene_counts), time_step=factor(environmental_vars), genotype=factor(genotype))
  if (allelic_df_1 == "NA") {
    num_sites = 0
  } else if (allelic_df_2 == "NA") {
    num_sites = 1
  } else if (allelic_df_3 == "NA") {
    num_sites = 2
  } else if (allelic_df_4 == "NA") {
    num_sites =3
  } else {
    num_sites = 4
  }
  num_sites=0
  output_file <- paste0(output_dir, parameter_string, "_", ensamble_id, "_", rs_id, "_dynamic_qtl_hit.png")

  te_plot <- make_te_plot(ensamble_id, rs_id, pvalue, beta, te_df, nb_conc)
  # ggsave(te_plot, file=output_file,width = 20,height=10.5,units="cm")
  if (num_sites > 0) {

    if (num_sites == 1) {
        as_plot1 <- make_as_plot(allelic_df_1, 1)
        as_plot_full <- plot_grid(as_plot1,NULL,NULL,NULL,nrow=1,ncol=4,label_size=10)
    } else if (num_sites == 2) {
        as_plot1 <- make_as_plot(allelic_df_1, 1)
        as_plot2 <- make_as_plot(allelic_df_2, 2)
        as_plot_full <- plot_grid(as_plot1,as_plot2,NULL,NULL,nrow=1,ncol=4,label_size=10)
    } else if (num_sites == 3) {
        as_plot1 <- make_as_plot(allelic_df_1, 1)
        as_plot2 <- make_as_plot(allelic_df_2, 2)
        as_plot3 <- make_as_plot(allelic_df_3, 3)
        as_plot_full <- plot_grid(as_plot1,as_plot2,as_plot3,NULL,nrow=1,ncol=4,label_size=10)
    } else if (num_sites == 4) {
        as_plot1 <- make_as_plot(allelic_df_1, 1)
        as_plot2 <- make_as_plot(allelic_df_2, 2)
        as_plot3 <- make_as_plot(allelic_df_3, 3)
        as_plot4 <- make_as_plot(allelic_df_4, 4)
        as_plot_full <- plot_grid(as_plot1, as_plot2, as_plot3, as_plot4,nrow=1,ncol=4,label_size=10)
    }
    merged <- plot_grid(te_plot, as_plot_full, nrow=2,ncol=1, label_size=10)
    #as_plot <- make_as_plot(allelic_df_1, 1)
    # as_plot <- make_as_plot_v2(allelic_df_1, 1)
    save_plot(output_file, merged, base_width=16, base_height=8)
  } else {
    ggsave(te_plot, file=output_file, width=20, height=10.5, units="cm")
  }

  ## Insert some if statement logic here
  if(length(next_line) == 0) {
    stop = TRUE
    close(f)
  }
}