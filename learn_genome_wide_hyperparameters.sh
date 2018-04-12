#!/bin/bash
#SBATCH --time=20:00:00 --partition=gilad --mem=25GB

module unload python
module load Anaconda3

joint_test_input_file="$1"
correction_factor_file="$2"
as_overdispersion_parameter_file="$3"
as_overdispersion_parameter_sample_specific_file="$4"
min_num_biallelic_lines="$5"
min_num_biallelic_samples="$6"
min_num_het_test_variant_biallelic_samples="$7"
te_nb_time_step_od_parameter_file="$8"
variance_output_file="$9"


date
python learn_library_size_correction_factor.py $joint_test_input_file $correction_factor_file
if false; then

date
python learn_as_overdispersion_parameter.py $joint_test_input_file $as_overdispersion_parameter_file $as_overdispersion_parameter_sample_specific_file $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples



python learn_time_step_nb_overdispersion_parameter.py $joint_test_input_file $te_nb_time_step_od_parameter_file $correction_factor_file $variance_output_file
fi