#!/bin/bash
#SBATCH --time=3:00:00 --partition=broadwl --mem=5GB


parameter_string="$1"
qtl_results_dir="$2"
target_region_input_file="$3"
qtl_visualization_dir="$4"
total_jobs="$5"
gencode_file="$6"



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





python run_gene_set_enrichment.py $real_merged_results $qtl_visualization_dir $parameter_string $gencode_file


if false; then
Rscript visualize_dynamic_qtls.R $real_merged_results $parameter_string $qtl_visualization_dir 
fi