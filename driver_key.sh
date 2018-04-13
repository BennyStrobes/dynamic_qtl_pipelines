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
output_root="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te/"

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
### 2. 'pseudotime_predictions_3': pseudotime predictions (from Nirmal) using hmm with 3 latent variables
### 3. 'pseudotime_predictions_4': pseudotime predictions (from Nirmal) using hmm with 4 latent variables
### 4. 'pseudotime_predictions_5': pseudotime predictions (from Nirmal) using hmm with 5 latent variables
### 5. 'uniform_4': 
### 6. 'median_pseudotime_4': 1 sample from each cell_line pseudotime state
### 7. 'time_steps_max_8': raw time format, but do not include any samples where time is greater than 8
### 8. 'time_steps_max_9': raw time format, but do not include any samples where time is greater than 9
environmental_variable_form="time_steps"
joint_test_input_file=$input_data_dir"joint_test_input_file_"$environmental_variable_form".txt"

if false; then
python create_joint_test_input_file.py $cht_input_file_dir $joint_test_input_file $environmental_variable_form $pseudotime_predictions_3_file $pseudotime_predictions_4_file $pseudotime_predictions_5_file $total_expression_file
fi





##########################################
# Step 2: Learn genome wide hyperparameters
##########################################

# File to contain learned library size correction factors for each sample (samples ordered same as $joint_test_input_file)
correction_factor_file=$input_data_dir"library_size_correction_factor_"$environmental_variable_form".txt"


# takes around 30 minutes to run
if false; then
sbatch learn_genome_wide_hyperparameters.sh $joint_test_input_file $correction_factor_file
fi



##########################################
# Step 3: Run dynamic QTL Model
##########################################

#################
# Parameters
#################
# Name of model. Current options are:
### 2. 'te_log_linear'
### 6. 'te_log_linear_quadratic_basis'
### 7. 'te_gaussian_process'
### 8. 'te_log_linear_cubic_control'
model_version="te_log_linear"

# Optimization method
optimization_method="LBFGS"

# covariates to use
### 1. "t15_troponin"
### 2. "none"
covariate_method="none"

# Genotype_version
### 1. "dosage"
### 2. "round"
genotype_version="round"

# String used in output files to keep track of parameters used
parameter_string=$model_version"_environmental_variable_"$environmental_variable_form"_optimizer_"$optimization_method"_genotype_"$genotype_version"_covariate_method_"$covariate_method



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
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $covariate_method $genotype_version
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
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $covariate_method $genotype_version
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
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $covariate_method $genotype_version
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
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $job_number $total_jobs $optimization_method $permutation_scheme $covariate_method $genotype_version
done
fi





permutation_scheme="shuffle_all"

sh multiple_testing_correction_and_visualization.sh $parameter_string $qtl_results_dir $target_region_input_file $qtl_visualization_dir $total_jobs $gencode_file $joint_test_input_file $correction_factor_file $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples






if false; then


permutation_scheme="sample_null"
sbatch multiple_testing_correction_and_visualization.sh $parameter_string $qtl_results_dir $target_region_input_file $qtl_visualization_dir $total_jobs $gencode_file $joint_test_input_file $correction_factor_file $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples



permutation_scheme="shuffle_lines"
sbatch multiple_testing_correction_and_visualization.sh $parameter_string $qtl_results_dir $target_region_input_file $qtl_visualization_dir $total_jobs $gencode_file $joint_test_input_file $correction_factor_file $permutation_scheme $min_num_biallelic_lines $min_num_biallelic_samples $min_num_het_test_variant_biallelic_samples




Rscript merge_permutation_scheme_qq_plot.R $qtl_results_dir $qtl_visualization_dir



fi







