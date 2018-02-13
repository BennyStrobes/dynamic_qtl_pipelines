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
max_sites="$9"
permutation_scheme="${10}"



#######################################################
# Part 1: Merge results from dynamic qtl parrallel runs
########################################################
# For real data:
real_parameter_string=$qtl_results_dir$parameter_string"_permute_False_permutation_scheme_none"
real_merged_results=$real_parameter_string"_merged_dynamic_qtl_results.txt"
python merge_dynamic_qtl_runs.py $real_parameter_string $real_merged_results $target_region_input_file $total_jobs



# For permuted data
perm_parameter_string=$qtl_results_dir$parameter_string"_permute_True_permutation_scheme_"$permutation_scheme
perm_merged_results=$perm_parameter_string"_merged_dynamic_qtl_results.txt"
python merge_dynamic_qtl_runs.py $perm_parameter_string $perm_merged_results $target_region_input_file $total_jobs



#######################################################
# Part 2: # Compute qvalues on null (permuted) data
########################################################
mod_parameter_string=$parameter_string"_permutation_scheme_"$permutation_scheme
# output file for qvalue results
qvalue_file=$qtl_results_dir$mod_parameter_string"_qvalue_results.txt"
# Run qvalue analysis
Rscript compute_qvalues.R $perm_merged_results $qvalue_file




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


if false; then

python visualize_significant_hits.py $parameter_string $qtl_results_dir $joint_test_input_file $correction_factor_file $max_sites $qtl_visualization_dir


python run_gene_set_enrichment.py $real_merged_results $significant_gene_results $qtl_visualization_dir $parameter_string $gencode_file
fi

Rscript visualize_dynamic_qtls.R $real_merged_results $perm_merged_results $mod_parameter_string $qtl_visualization_dir 
