#!/bin/bash
#SBATCH --time=5:00:00 --partition=broadwl --mem=15GB
module unload python
module load Anaconda3

joint_test_input_file="$1"
model_version="$2"
output_stem="$3"
permute="$4"
job_number="$5"
total_jobs="$6"
optimization_method="$7"
permutation_scheme="${8}"
min_num_biallelic_lines="${9}"
min_num_biallelic_samples="${10}" 
min_num_het_test_variant_biallelic_samples="${11}"
covariate_method="${12}"

python dynamic_qtl_shell.py $joint_test_input_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples $covariate_method
