args = commandArgs(trailingOnly=TRUE)
library(ggplot2)




pvalue_histogram_plot <- function(dynamic_qtl_results, output_file) {
    histogram <- ggplot(data=dynamic_qtl_results, aes(dynamic_qtl_results$pvalue)) + 
        geom_histogram(breaks=seq(0, 1, by=.05)) + 
        labs(x="pvalue", y="Count")

    ggsave(histogram, file=output_file,width = 20,height=10.5,units="cm")

}

qq_plot <- function(dynamic_qtl_results, output_file) {
    # Sort real pvalues
    pvalues <- sort(dynamic_qtl_results$pvalue)
    # Simulate uniform distribution
    uniform <- sort(runif(length(pvalues)))

    df <- data.frame(pvalues=-log10(pvalues + .000000000001), expected_pvalues=-log10(uniform + .000000000001))

    # PLOT!
    max_val <-max(max(-log10(uniform + .000000000001)), max(-log10(pvalues + .000000000001)))
    #PLOT!
    scatter <- ggplot(df, aes(x = expected_pvalues, y = pvalues)) + geom_point() 
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(x = "Uniform", y = "Real")
    scatter <- scatter + geom_abline() +  theme(legend.position="none")
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))
    
    ggsave(scatter, file=output_file, width=20, height=10.5,units="cm")

}

qtl_results_dir = args[1]
qtl_visualization_dir = args[2]
parameter_string = args[3]


dynamic_qtl_results_file = paste0(qtl_results_dir, parameter_string, "_dynamic_qtl_results.txt")
dynamic_qtl_results <- read.table(dynamic_qtl_results_file, header=TRUE)


####################################
# Make qqplot
####################################
output_file <- paste0(qtl_visualization_dir, parameter_string, "_qq_plot.png")
qq_plot(dynamic_qtl_results, output_file)

####################################
# Make Pvalue histogram
####################################
output_file <- paste0(qtl_visualization_dir, parameter_string, "_pvalue_histogram.png")
pvalue_histogram_plot(dynamic_qtl_results, output_file)