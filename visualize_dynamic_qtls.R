args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(reshape)
library(cowplot)




pvalue_histogram <- function(pvalues, output_file) {
    df <- data.frame(pvalue=pvalues)

    histo <- ggplot(df, aes(x = pvalues)) +
        geom_histogram(aes(y = ..density..), breaks = seq(0,1,by=.01)) +
        theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        labs(x = "pvalue", y = "density", title="nominal pvalues")
    ggsave(histo, file=output_file,width = 20,height=10.5,units="cm")
}


qq_plot <- function(pvalues, output_file) {
    real_pvalues <- sort(pvalues)

    uniform_pvalues <- sort(runif(length(real_pvalues)))


    df <- data.frame(pvalues=-log10(real_pvalues + .000000000001), expected_pvalues=-log10(uniform_pvalues + .000000000001))
    
    # PLOT!
    max_val <-max(max(-log10(uniform_pvalues + .000000000001)), max(-log10(real_pvalues + .000000000001)))

    #PLOT!
    scatter <- ggplot(df, aes(x = expected_pvalues, y = pvalues)) + geom_point() 
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(x = "Uniform", y = "Real")
    scatter <- scatter + geom_abline() +  theme(legend.position="none")
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))

    ggsave(scatter, file=output_file, width=20, height=10.5, units="cm")
}





# Command Line Args
real_merged_results = args[1]
parameter_string = args[2]
qtl_visualization_dir = args[3]

# Load in real results
real_data <- read.table(real_merged_results,header=TRUE)





#############################
# Make pvalue histogram
#############################
output_file <- paste0(qtl_visualization_dir, parameter_string, "_pvalue_histogram.png")
pvalue_histogram(real_data$pvalue, output_file)


#############################
# Make qqplot
#############################
output_file <- paste0(qtl_visualization_dir, parameter_string, "_qq_plot.png")
qq_plot(real_data$pvalue, output_file)