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



#######################################################
# Part 1: Merge results from dynamic qtl parrallel runs
########################################################
# For real data:
real_parameter_string=$qtl_results_dir$parameter_string"_permute_False"
real_merged_results=$real_parameter_string"_merged_dynamic_qtl_results.txt"
if false; then
python merge_dynamic_qtl_runs.py $real_parameter_string $real_merged_results $target_region_input_file $total_jobs
fi
# For permuted data
perm_parameter_string=$qtl_results_dir$parameter_string"_permute_True"
perm_merged_results=$perm_parameter_string"_merged_dynamic_qtl_results.txt"
if false; then
python merge_dynamic_qtl_runs.py $perm_parameter_string $perm_merged_results $target_region_input_file $total_jobs
fi



#######################################################
# Part 2: # Compute qvalues on null (permuted) data
########################################################
# output file for qvalue results
qvalue_file=$qtl_results_dir$parameter_string"_qvalue_results.txt"
# Run qvalue analysis
if false; then
Rscript compute_qvalues.R $perm_merged_results $qvalue_file
fi



#######################################################
# Part 3: # Assess genome wide significance of actual data based on the qvalue threshold (FDR <= thresh) in null data
########################################################
# FDR threshold
qval_thresh=".05"
# Output file for all significant variant gene pairs
significant_results=$qtl_results_dir$parameter_string"_qval_"$qval_thresh"_significant.txt"
# Output file for significant egenes and their strongest associated variant
significant_gene_results=$qtl_results_dir$parameter_string"_qval_"$qval_thresh"_significant_egenes.txt"
if false; then
python assess_significance.py $perm_merged_results $real_merged_results $qvalue_file $significant_results $significant_gene_results $qval_thresh
fi




python visualize_significant_hits.py $parameter_string $qtl_results_dir $joint_test_input_file $correction_factor_file $max_sites $qtl_visualization_dir




if false; then
python run_gene_set_enrichment.py $real_merged_results $qtl_visualization_dir $parameter_string $gencode_file
fi

if false; then
Rscript visualize_dynamic_qtls.R $real_merged_results $perm_merged_results $parameter_string $qtl_visualization_dir 
fi