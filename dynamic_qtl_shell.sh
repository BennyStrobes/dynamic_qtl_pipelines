#!/bin/bash
#SBATCH --time=22:00:00 --partition=broadwl --mem=15GB
#SBATCH -o /project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te_quadratic/log_files/slurm.%N.%j.out # STDOUT
#SBATCH -e /project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te_quadratic/log_files/slurm.%N.%j.err # STDERR


module unload python
module load Anaconda3

joint_test_input_file="$1"
correction_factor_file="$2"
model_version="$3"
output_stem="$4"
permute="$5"
job_number="$6"
total_jobs="$7"
optimization_method="$8"
permutation_scheme="${9}"
covariate_method="${10}"
genotype_version="${11}"

python dynamic_qtl_shell.py $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $covariate_method $genotype_version
