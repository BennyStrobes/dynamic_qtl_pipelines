#######################################################################################################
# This is the driver script that runs EAGLE model on ipsc data
#######################################################################################################



###############################################################################
# Input Data
###############################################################################

# Directory created by "time_step_independent_qtl_pipelines" scripts
# Contains 1 file per sample with information on each test (variant, target region)
# Each file (sample) has the same number of lines (tests)
cht_input_file_dir="/project2/gilad/bstrober/ipsc_differentiation/time_step_independent_qtl_pipelines/wasp/cht_input_files/"

# File containing all of the target regions we are using. We are using this file to convert from gene positions to ensable id
target_region_input_file="/project2/gilad/bstrober/ipsc_differentiation/time_step_independent_qtl_pipelines/wasp/target_regions/target_regions_cis_distance_50000_maf_cutoff_0.1_min_reads_90_min_as_reads_20_min_het_counts_4_merged.txt"

# File containing conversions from ensamble ids to gene symbols
gencode_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/gencode.v19.annotation.gtf.gz"

# Files containing mapping from sample id to Nirmal's pseudotime predictions
# 3 state HMM
pseudotime_predictions_3_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/pseudotime_predictions_14_lines/3_state_output_14.csv"
# 4 state HMM
pseudotime_predictions_4_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/pseudotime_predictions_14_lines/4_state_output_14.csv"
# 5 state HMM
pseudotime_predictions_5_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/pseudotime_predictions_14_lines/5_state_output_14.csv"


total_expression_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess/processed_total_expression/quantile_normalized.txt"

###############################################################################
# Output directories (aasume all of these exist prior to starting analysis)
###############################################################################

# Root directory for this of all ipsc data based results
output_root="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_as/"

# Directory containing necessary input files to qtl tests
input_data_dir=$output_root"input_data/"


# Directory containing text files with results from dynamic qtl analysis
qtl_results_dir=$output_root"qtl_results/"

# Directory containing visualization of results found in qtl_results_dir
qtl_visualization_dir=$output_root"qtl_visualization/"






###############################################################################
# Dynamic QTL Calling
###############################################################################


##########################################
# Step 1: Create joint_test_input_file
##########################################
# joint_test_input_file is a table with 1 row per sample
# each row has 3 columns:
#### 1. Sample id
#### 2. environmental variable
#### 3. Absolute directory to CHT input file for that sample
# NOTE: This script is very specific to our data
# Takes less than a minute to run
# $environmental_variable is a parameter that describes how we parameterize the environmental variable. So far, this is done with:
### 1. 'time_steps': raw time format
environmental_variable_form="time_steps"
joint_test_input_file=$input_data_dir"joint_test_input_file_"$environmental_variable_form".txt"
if false; then
python create_joint_test_input_file.py $cht_input_file_dir $joint_test_input_file $environmental_variable_form $pseudotime_predictions_3_file $pseudotime_predictions_4_file $pseudotime_predictions_5_file $total_expression_file
fi






##########################################
# Step 3: Run dynamic QTL Model
##########################################
# Minimum number of cell lines that have at least one biallelic sample required for a site to be used
min_num_biallelic_lines="5"
# Minimum number of biallelic samples required for a site to be used
min_num_biallelic_samples="25"
# Minimum number of biallelic samples that have a heterozygous test variant required for a site to be used
min_num_het_test_variant_biallelic_samples="10"



#################
# Parameters
#################
# Name of model. Current options are:
### 1. 'joint_log_linear'
### 2. 'te_log_linear'
### 3. 'as_log_linear'
### 4. 'as_log_linear_fixed_sample_overdispersion'
### 5. 'as_log_linear_fixed_overdispersion'
### 6. 'te_log_linear_quadratic_basis'
### 7. 'te_gaussian_process'
### 8. 'te_log_linear_cubic_control'
### 9. 'te_log_linear_time_od'
### 10. 'te_log_linear_time_od_offline'
model_version="te_log_linear"

# Optimization method
optimization_method="LBFGS"

# covariates to use
### 1. "t15_troponin"
### 2. "none"
covariate_method="none"

# String used in output files to keep track of parameters used
parameter_string=$model_version"_environmental_variable_"$environmental_variable_form"_biallelic_lines_"$min_num_biallelic_lines"_biallelic_samples_"$min_num_biallelic_samples"_biallelic_test_het_samples_"$min_num_het_test_variant_biallelic_samples"_optimizer_"$optimization_method"_covariate_method_"$covariate_method
#parameter_string=$model_version"_environmental_variable_"$environmental_variable_form"_biallelic_lines_"$min_num_biallelic_lines"_biallelic_samples_"$min_num_biallelic_samples"_biallelic_test_het_samples_"$min_num_het_test_variant_biallelic_samples"_optimizer_"$optimization_method



# How many parallel nodes at once
total_jobs="50"

##################
# Run analysis
##################

##########################
# Run on Real data
permute="False"
permutation_scheme="none"
#########################
if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
    # Stem of all output files
    output_stem=$qtl_results_dir$parameter_string"_permutation_scheme_"$permutation_scheme"_permute_"$permute"_"$job_number"_"
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples $as_overdispersion_parameter_file $as_overdispersion_parameter_sample_specific_file $covariate_method $te_nb_time_step_od_parameter_file
done
fi





##########################
# Run permutation for all samples
permute="True"
permutation_scheme="shuffle_all"
##########################
if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
    # Stem of all output files
    output_stem=$qtl_results_dir$parameter_string"_permutation_scheme_"$permutation_scheme"_permute_"$permute"_"$job_number"_"
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples $as_overdispersion_parameter_file $as_overdispersion_parameter_sample_specific_file $covariate_method $te_nb_time_step_od_parameter_file
done
fi

 













##########################
# Run permutation for all samples
permute="True"
permutation_scheme="sample_null"

##########################
if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
    # Stem of all output files
    output_stem=$qtl_results_dir$parameter_string"_permutation_scheme_"$permutation_scheme"_permute_"$permute"_"$job_number"_"
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples $as_overdispersion_parameter_file $as_overdispersion_parameter_sample_specific_file
done
fi







##########################
# Run permutation independently in each cell line
permute="True"
permutation_scheme="shuffle_lines"
##########################
if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
    # Stem of all output files
    output_stem=$qtl_results_dir$parameter_string"_permutation_scheme_"$permutation_scheme"_permute_"$permute"_"$job_number"_"
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples $as_overdispersion_parameter_file $as_overdispersion_parameter_sample_specific_file
done
fi





##########################
# Run permutation independently in heterozygotes and homozygotes
permute="True"
permutation_scheme="shuffle_lines_ordered"
##########################
if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
    # Stem of all output files
    output_stem=$qtl_results_dir$parameter_string"_permutation_scheme_"$permutation_scheme"_permute_"$permute"_"$job_number"_"
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples $as_overdispersion_parameter_file $as_overdispersion_parameter_sample_specific_file $covariate_method
done
fi



##########################
# Run permutation independently in heterozygotes and homozygotes
permute="True"
permutation_scheme="shuffle_hets"
##########################
if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
    # Stem of all output files
    output_stem=$qtl_results_dir$parameter_string"_permutation_scheme_"$permutation_scheme"_permute_"$permute"_"$job_number"_"
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples $as_overdispersion_parameter_file $as_overdispersion_parameter_sample_specific_file
done
fi







min_num_biallelic_lines="2"
# Minimum number of biallelic samples required for a site to be used
min_num_biallelic_samples="15"
# Minimum number of biallelic samples that have a heterozygous test variant required for a site to be used
min_num_het_test_variant_biallelic_samples="6"


if false; then


permutation_scheme="shuffle_all"
sh multiple_testing_correction_and_visualization.sh $parameter_string $qtl_results_dir $target_region_input_file $qtl_visualization_dir $total_jobs $gencode_file $joint_test_input_file $correction_factor_file $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples



permutation_scheme="shuffle_lines_ordered"
sh multiple_testing_correction_and_visualization.sh $parameter_string $qtl_results_dir $target_region_input_file $qtl_visualization_dir $total_jobs $gencode_file $joint_test_input_file $correction_factor_file $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples






permutation_scheme="sample_null"
sbatch multiple_testing_correction_and_visualization.sh $parameter_string $qtl_results_dir $target_region_input_file $qtl_visualization_dir $total_jobs $gencode_file $joint_test_input_file $correction_factor_file $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples



permutation_scheme="shuffle_lines"
sbatch multiple_testing_correction_and_visualization.sh $parameter_string $qtl_results_dir $target_region_input_file $qtl_visualization_dir $total_jobs $gencode_file $joint_test_input_file $correction_factor_file $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples


permutation_scheme="shuffle_hets"
sbatch multiple_testing_correction_and_visualization.sh $parameter_string $qtl_results_dir $target_region_input_file $qtl_visualization_dir $total_jobs $gencode_file $joint_test_input_file $correction_factor_file $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples


permutation_scheme="shuffle_lines_ordered"
sh multiple_testing_correction_and_visualization.sh $parameter_string $qtl_results_dir $target_region_input_file $qtl_visualization_dir $total_jobs $gencode_file $joint_test_input_file $correction_factor_file $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples


Rscript merge_permutation_scheme_qq_plot.R $qtl_results_dir $qtl_visualization_dir
fi