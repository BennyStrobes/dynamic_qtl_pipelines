#!/bin/bash
#SBATCH --time=20:00:00 --partition=broadwl --mem=5GB


joint_test_input_file="$1"
correction_factor_file="$2"
model_version="$3"
output_stem="$4"
permute="$5"
max_sites="$6"
job_number="$7"
total_jobs="$8"

python dynamic_qtl_shell.py $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $max_sites $job_number $total_jobs