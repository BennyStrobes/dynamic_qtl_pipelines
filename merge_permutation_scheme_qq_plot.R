args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(reshape)


sample_non_significant_hits <- function(pvalues, fraction_kept=.01,fraction_sampled=.001) {
    index <- floor(length(pvalues)*fraction_kept)
    to_keep <- pvalues[1:index]
    to_filter <- pvalues[(index+1):length(pvalues)]
    filtered <- sort(sample(to_filter,floor(length(to_filter)*fraction_sampled)))
    return(c(to_keep,filtered))
}



extract_pvalz <- function(file_name) {
    data <- read.table(file_name, header=TRUE)
    pvalues <- sort(data$pvalue)
    sampled_pvalues <- sample_non_significant_hits(pvalues)
    return(sampled_pvalues)
}

extract_uniform_pvalz <- function(file_name) {
    data <- read.table(file_name, header=TRUE)
    pvalues <- sort(data$pvalue)
    num_samples <- length(pvalues)
    uniform_pvalues <- sort(runif(length(pvalues)))
    sampled_uniform_pvalues <- sample_non_significant_hits(uniform_pvalues)
    return(sampled_uniform_pvalues)
}



# Command Line Args
parameter_string = args[1]
results_dir = args[2]
qtl_visualization_dir = args[3]


# Load in real results
real_results <- paste0(results_dir, parameter_string, "_permute_False_permutation_scheme_none_merged_dynamic_qtl_results.txt")
# Load in permuted results
perm_sample_null_results <- paste0(results_dir, parameter_string, "_permute_True_permutation_scheme_sample_null_merged_dynamic_qtl_results.txt")
perm_shuffle_all_results <- paste0(results_dir, parameter_string, "_permute_True_permutation_scheme_shuffle_all_merged_dynamic_qtl_results.txt")
perm_shuffle_hets_results <- paste0(results_dir, parameter_string, "_permute_True_permutation_scheme_shuffle_hets_merged_dynamic_qtl_results.txt")
perm_shuffle_lines_results <- paste0(results_dir, parameter_string, "_permute_True_permutation_scheme_shuffle_lines_merged_dynamic_qtl_results.txt")


real_pvalz <- extract_pvalz(real_results)
uniform_pvalz <- extract_uniform_pvalz(real_results)

perm_sample_null_pvalz <- extract_pvalz(perm_sample_null_results)
perm_shuffle_all_pvalz <- extract_pvalz(perm_shuffle_all_results)
perm_shuffle_hets_pvalz <- extract_pvalz(perm_shuffle_hets_results)
perm_shuffle_lines_pvalz <- extract_pvalz(perm_shuffle_lines_results)


all_pvalues <- c(real_pvalz, perm_sample_null_pvalz, perm_shuffle_all_pvalz, perm_shuffle_hets_pvalz, perm_shuffle_lines_pvalz)
type <- c(rep("real data",length(real_pvalz)), rep("perm: sample null",length(perm_sample_null_pvalz)), rep("perm: shuffle all",length(perm_shuffle_all_pvalz)), rep("perm: shuffle hets",length(perm_shuffle_hets_pvalz)), rep("perm: shuffle lines",length(perm_shuffle_lines_pvalz)))
uniform <- c(uniform_pvalz, uniform_pvalz, uniform_pvalz, uniform_pvalz, uniform_pvalz)

df <- data.frame(pvalues=-log10(all_pvalues + .000000000001), expected_pvalues=-log10(uniform + .000000000001), type=factor(type))

# PLOT!
max_val <-max(max(-log10(uniform + .000000000001)), max(-log10(all_pvalues + .000000000001)))
#PLOT!
scatter <- ggplot(df, aes(x = expected_pvalues, y = pvalues, colour = type)) + geom_point() 
scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
scatter <- scatter + labs(colour="Permutation Scheme",x = "Uniform", y = "Real")
scatter <- scatter + geom_abline()
scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))

output_file <- paste0(qtl_visualization_dir, parameter_string, "_qq_plot_compare_permutation_schemes.png")
ggsave(scatter, file=output_file,width = 20,height=10.5,units="cm")





