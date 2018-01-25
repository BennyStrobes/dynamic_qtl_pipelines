source("betabinomial_glmm.R")
args = commandArgs(trailingOnly=TRUE)




plot_elbo <- function(full_elbo, null_elbo, gene_id, elbo_plot_file, pvalue) {
    df <- data.frame(elbo = c(full_elbo,null_elbo),iteration=c(1:length(full_elbo),1:length(null_elbo)), type = factor(c(rep("full",length(full_elbo)),rep("null",length(null_elbo)))))

    #PLOT!
    scatter <- ggplot(df, aes(x = iteration, y = elbo, colour = type)) + geom_point(size=.3) 
    scatter <- scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(colour="Model type",x = "iteration", y = "elbo", title = paste0(gene_id, " / ", pvalue))
    ggsave(scatter, file=elbo_plot_file,width = 15,height=10.5,units="cm")

}

dynamic_qtl_shell <- function(data, gene_num, model_version) {
    time_steps <- as.numeric(strsplit(as.character(data[gene_num,7]),",")[[1]])
    z <- as.integer(strsplit(as.character(data[gene_num,9]),",")[[1]])

    ys <- as.integer(strsplit(as.character(data[gene_num,10]),",")[[1]])
    ns <- as.integer(strsplit(as.character(data[gene_num,11]),",")[[1]])


    xNull <- cbind(numeric(length(time_steps))+1)
    xFull <- cbind(xNull,time_steps)

    if (model_version == "glm") {
        test_info <- betaBinomialGLM(ys,ns,xFull,xNull)

        pvalue <- test_info$lrtp
        loglr <- test_info$loglr
        beta <- test_info$fit_full$par$beta[2]

        results <- list(pvalue=pvalue, loglr=loglr, beta=beta)
    } else if (model_version == "glmm") {
        test_info <- betaBinomialGLMM(ys,ns,xFull,xNull,z)


        full_elbo <- test_info$v_full$elbo_progress[70:length(test_info$v_full$elbo_progress)]
        null_elbo <- test_info$v_null$elbo_progress[70:length(test_info$v_null$elbo_progress)]

        pvalue <- test_info$lrtp
        loglr <- test_info$loglr
        beta <- test_info$v_full$m[4]

        results <- list(pvalue=pvalue, loglr=loglr, beta=beta,full_elbo=full_elbo, null_elbo=null_elbo)
    }

    return(results)
}










#####################################
# Command Line Arguments
#####################################
qtl_results_dir = args[1]  # Directory to save results to
qtl_visualization_dir = args[2]  # Directory to save plots to
parameter_string = args[3] # String summarizing parameters. (to be used for file names)
model_version = args[4]  # Name of EAGLE model
elbo_diagnostic_dir = args[5]  # Directory to save elbo plots

# Input file and output file
input_file <- paste0(qtl_results_dir, parameter_string, "_dynamic_qtl_input_file.txt")
output_file <- paste0(qtl_results_dir, parameter_string, "_dynamic_qtl_results.txt")

# Load in input data
data <- read.table(input_file,header=TRUE)
num_genes <- dim(data)[1]

# Make output data frame
output_data <- cbind(data,pvalue=2,loglr=-1,beta=0)

for (gene_num in 1:num_genes){
    print(gene_num)
    gene_id <- data[gene_num, 6]
    results <- dynamic_qtl_shell(data, gene_num, model_version)
    output_data[gene_num,12] <- results$pvalue
    output_data[gene_num,13] <- results$loglr
    output_data[gene_num,14] <- results$beta

    if (model_version == "glmm") {
        elbo_plot_file <- paste0(elbo_diagnostic_dir, parameter_string,"_elbo_plot_", gene_id,".png")
        plot_elbo(results$full_elbo, results$null_elbo, gene_id, elbo_plot_file, results$pvalue)
    }

}



write.table(output_data,output_file, quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)