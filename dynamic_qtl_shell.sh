#!/bin/bash
#SBATCH --time=5:00:00 --partition=broadwl --mem=15GB
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
min_num_biallelic_lines="${10}"
min_num_biallelic_samples="${11}" 
min_num_het_test_variant_biallelic_samples="${12}"
as_overdispersion_parameter_file="${13}"
as_overdispersion_parameter_sample_specific_file="${14}"
covariate_method="${15}"
te_nb_time_step_od_parameter_file="${16}"

python dynamic_qtl_shell.py $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $max_sites $job_number $total_jobs $optimization_method $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples $as_overdispersion_parameter_file $as_overdispersion_parameter_sample_specific_file $covariate_method $te_nb_time_step_od_parameter_file
