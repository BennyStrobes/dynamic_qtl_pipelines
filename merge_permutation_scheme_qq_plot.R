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

results_1 <- paste0(results_dir, "te_log_linear_quadratic_basis_environmental_variable_time_steps_optimizer_Newton_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_none_permute_False_merged_dynamic_qtl_results.txt")
results_2 <- paste0(results_dir, "te_log_linear_quadratic_basis_environmental_variable_time_steps_optimizer_Newton_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_shuffle_all_permute_True_merged_dynamic_qtl_results.txt")
results_3 <- paste0(results_dir, "te_log_linear_quadratic_basis_environmental_variable_time_steps_optimizer_Newton_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_shuffle_lines_permute_True_merged_dynamic_qtl_results.txt")
results_4 <- paste0(results_dir, "te_log_linear_quadratic_basis_environmental_variable_time_steps_optimizer_Newton_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_sample_null_permute_True_merged_dynamic_qtl_results.txt")


pvalz_1 <- extract_pvalz(results_1)
pvalz_2 <- extract_pvalz(results_2)
pvalz_3 <- extract_pvalz(results_3)
pvalz_4 <- extract_pvalz(results_4)


uniform_pvalz1 <- extract_uniform_pvalz(results_1)
uniform_pvalz2 <- extract_uniform_pvalz(results_2)
uniform_pvalz3 <- extract_uniform_pvalz(results_3)
uniform_pvalz4 <- extract_uniform_pvalz(results_4)





all_pvalues <- c(pvalz_1, pvalz_2, pvalz_3, pvalz_4)
type <- c(rep("real",length(pvalz_1)), rep("shuffle_all",length(pvalz_2)), rep("shuffle_lines",length(pvalz_3)), rep("sample_null",length(pvalz_4)))
uniform <- c(uniform_pvalz1, uniform_pvalz2, uniform_pvalz3, uniform_pvalz4)

df <- data.frame(pvalues=-log10(all_pvalues + .000000000001), expected_pvalues=-log10(uniform + .000000000001), type=factor(type))

# PLOT!
max_val <-max(max(-log10(uniform + .000000000001)), max(-log10(all_pvalues + .000000000001)))
#PLOT!
scatter <- ggplot(df, aes(x = expected_pvalues, y = pvalues, colour = type)) + geom_point() 
scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
scatter <- scatter + labs(colour="Type",x = "Uniform", y = "Real")
scatter <- scatter + geom_abline()
scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))

output_file <- paste0(qtl_visualization_dir, "qq_plot_te_log_linear_quadratic_basis_environmental_variable_time_steps_genotype_round_covariate_method_comparison.png")
ggsave(scatter, file=output_file,width = 20,height=10.5,units="cm")





