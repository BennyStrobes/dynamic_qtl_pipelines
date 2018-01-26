#!/bin/bash
#SBATCH --time=2:00:00 --partition=broadwl --mem=10GB


joint_test_input_file="$1"
correction_factor_file="$2"
model_version="$3"
output_stem="$4"
job_number="$5"
total_jobs="$6"

python dynamic_qtl_shell.py $joint_test_input_file $correction_factor_file $model_version $output_stem $job_number $total_jobs