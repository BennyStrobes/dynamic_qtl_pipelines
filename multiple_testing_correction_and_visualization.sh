#!/bin/bash
#SBATCH --time=3:00:00 --partition=broadwl --mem=5GB


parameter_string="$1"
qtl_results_dir="$2"
target_region_input_file="$3"
qtl_visualization_dir="$4"
total_jobs="$5"
gencode_file="$6"
joint_test_input_file="$7"
correction_factor_file="$8"
permutation_scheme="${9}"
min_num_biallelic_lines="${10}"
min_num_biallelic_samples="${11}"
min_num_het_test_variant_biallelic_samples="${12}"


#######################################################
# Part 1: Merge results from dynamic qtl parrallel runs
########################################################
# For real data:

real_parameter_string=$qtl_results_dir$parameter_string"_permutation_scheme_none_permute_False"

real_merged_results=$real_parameter_string"_merged_dynamic_qtl_results.txt"
# python merge_dynamic_qtl_runs.py $real_parameter_string $real_merged_results $target_region_input_file $total_jobs



# For permuted data
perm_parameter_string=$qtl_results_dir$parameter_string"_permutation_scheme_"$permutation_scheme"_permute_True"
perm_merged_results=$perm_parameter_string"_merged_dynamic_qtl_results.txt"
# python merge_dynamic_qtl_runs.py $perm_parameter_string $perm_merged_results $target_region_input_file $total_jobs

mod_parameter_string=$parameter_string"_permutation_scheme_"$permutation_scheme


# Rscript visualize_dynamic_qtls.R $real_merged_results $perm_merged_results $mod_parameter_string $qtl_visualization_dir 

#######################################################
# Part 2: Run FDR correction using real and permuted
########################################################
fdr_thresh=".05"
echo "starting"

######FDR METHOD 2
# output file for eFDR analysis
efdr_file=$qtl_results_dir$mod_parameter_string"_eFDR_results.txt"
# Run eFDR correction
# Rscript eFDR_correction.R $real_merged_results $perm_merged_results $efdr_file


# Assess genome wide significance of actual data based on the eFDR approach with FDR <= $fdr_thresh
# Output file for all significant variant gene pairs
significant_efdr_results=$qtl_results_dir$mod_parameter_string"_efdr_"$fdr_thresh"_significant.txt"
# Output file for significant egenes and their strongest associated variant
significant_efdr_gene_results=$qtl_results_dir$mod_parameter_string"_efdr_"$fdr_thresh"_significant_egenes.txt"
# python assess_significance_efdr_approach.py $efdr_file $real_merged_results $significant_efdr_results $significant_efdr_gene_results $fdr_thresh



echo "Start"
#python organize_significant_hits.py $mod_parameter_string $qtl_results_dir $joint_test_input_file $correction_factor_file $qtl_visualization_dir $parameter_string $target_region_input_file



Rscript visualize_significant_hits.R $mod_parameter_string $qtl_visualization_dir "/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te/temper_debug/te_round_genotype/"



if false; then














#######################################################
# Part 3: # Assess genome wide significance of actual data based on the qvalue threshold (FDR <= thresh) in null data
########################################################
# FDR threshold
qval_thresh=".05"
# Output file for all significant variant gene pairs
significant_results=$qtl_results_dir$mod_parameter_string"_qval_"$qval_thresh"_significant.txt"
# Output file for significant egenes and their strongest associated variant
significant_gene_results=$qtl_results_dir$mod_parameter_string"_qval_"$qval_thresh"_significant_egenes.txt"
python assess_significance.py $perm_merged_results $real_merged_results $qvalue_file $significant_results $significant_gene_results $qval_thresh

python organize_significant_hits.py $mod_parameter_string $qtl_results_dir $joint_test_input_file $correction_factor_file $qtl_visualization_dir $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples $parameter_string $target_region_input_file

Rscript visualize_significant_hits.R $mod_parameter_string $qtl_visualization_dir

python run_gene_set_enrichment.py $real_merged_results $significant_gene_results $qtl_visualization_dir $parameter_string $gencode_file


Rscript visualize_dynamic_qtls.R $real_merged_results $perm_merged_results $mod_parameter_string $qtl_visualization_dir 
fi