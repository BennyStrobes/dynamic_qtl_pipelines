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
results_dir = args[1]
qtl_visualization_dir = args[2]


# Load in results
results_1 <- paste0(results_dir, "te_log_linear_quadratic_basis_environmental_variable_pseudotime_predictions_4_biallelic_lines_2_biallelic_samples_15_biallelic_test_het_samples_6_optimizer_LBFGS_permutation_scheme_none_permute_False_merged_dynamic_qtl_results.txt")
results_2 <- paste0(results_dir, "te_log_linear_quadratic_basis_environmental_variable_pseudotime_predictions_4_biallelic_lines_2_biallelic_samples_15_biallelic_test_het_samples_6_optimizer_LBFGS_permutation_scheme_sample_null_permute_True_merged_dynamic_qtl_results.txt")
results_3 <- paste0(results_dir, "te_log_linear_quadratic_basis_environmental_variable_pseudotime_predictions_4_biallelic_lines_2_biallelic_samples_15_biallelic_test_het_samples_6_optimizer_LBFGS_permutation_scheme_shuffle_all_permute_True_merged_dynamic_qtl_results.txt")
results_4 <- paste0(results_dir, "te_log_linear_quadratic_basis_environmental_variable_pseudotime_predictions_4_biallelic_lines_2_biallelic_samples_15_biallelic_test_het_samples_6_optimizer_LBFGS_permutation_scheme_shuffle_hets_permute_True_merged_dynamic_qtl_results.txt")
results_5 <- paste0(results_dir, "te_log_linear_quadratic_basis_environmental_variable_pseudotime_predictions_4_biallelic_lines_2_biallelic_samples_15_biallelic_test_het_samples_6_optimizer_LBFGS_permutation_scheme_shuffle_lines_permute_True_merged_dynamic_qtl_results.txt")


pvalz_1 <- extract_pvalz(results_1)
pvalz_2 <- extract_pvalz(results_2)
pvalz_3 <- extract_pvalz(results_3)
pvalz_4 <- extract_pvalz(results_4)
pvalz_5 <- extract_pvalz(results_5)

uniform_pvalz <- extract_uniform_pvalz(results_1)



all_pvalues <- c(pvalz_1, pvalz_2, pvalz_3, pvalz_4, pvalz_5)
type <- c(rep("real",length(pvalz_1)), rep("perm: sample null",length(pvalz_2)), rep("perm: shuffle all",length(pvalz_3)), rep("perm: shuffle hets",length(pvalz_4)), rep("perm: shuffle lines",length(pvalz_5)))
uniform <- c(uniform_pvalz, uniform_pvalz, uniform_pvalz, uniform_pvalz,uniform_pvalz)

df <- data.frame(pvalues=-log10(all_pvalues + .000000000001), expected_pvalues=-log10(uniform + .000000000001), type=factor(type))

# PLOT!
max_val <-max(max(-log10(uniform + .000000000001)), max(-log10(all_pvalues + .000000000001)))
#PLOT!
scatter <- ggplot(df, aes(x = expected_pvalues, y = pvalues, colour = type)) + geom_point() 
scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
scatter <- scatter + labs(colour="Type",x = "Uniform", y = "Real")
scatter <- scatter + geom_abline()
scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))

output_file <- paste0(qtl_visualization_dir, "qq_plot_quadratic_basis_pseudotime4_compare_perm_schemes.png")
ggsave(scatter, file=output_file,width = 20,height=10.5,units="cm")





