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
target_region_input_file="/project2/gilad/bstrober/ipsc_differentiation/time_step_independent_qtl_pipelines/wasp/target_regions/target_regions_cis_distance_50000_maf_cutoff_0.1_min_reads_60_min_as_reads_10_min_het_counts_2_merged.txt"

# File containing conversions from ensamble ids to gene symbols
gencode_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/gencode.v19.annotation.gtf.gz"


###############################################################################
# Output directories (aasume all of these exist prior to starting analysis)
###############################################################################

# Root directory for this of all ipsc data based results
output_root="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data/"

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
python create_joint_test_input_file.py $cht_input_file_dir $joint_test_input_file $environmental_variable_form
fi








##########################################
# Step 2: Learn library size correction factor
##########################################
# File to contain learned library size correction factors for each sample (samples ordered same as $joint_test_input_file)
correction_factor_file=$input_data_dir"library_size_correction_factor.txt"
# takes around 30 minutes to run
if false; then
sbatch learn_library_size_correction_factor.sh $joint_test_input_file $correction_factor_file
fi




##########################################
# Step 3: Run dynamic QTL Model
##########################################

#################
# Parameters
#################
# Name of model. Current options are:
### 1. 'joint_log_linear'
model_version="joint_log_linear"


# Maximum number of exonic sites to allow per gene. If there are more than max_sites, take the top $max_sites highest expressing sites
max_sites="30"

# String used in output files to keep track of parameters used
parameter_string=$model_version"_environmental_variable_"$environmental_variable_form"_max_sites_"$max_sites

# How many parallel nodes at once
total_jobs="300"

##################
# Run analysis
##################
permute="False"
for job_number in $(seq 0 $(($total_jobs-1))); do 
    # Stem of all output files
    output_stem=$qtl_results_dir$parameter_string"_permute_"$permute"_"$job_number"_"
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $max_sites $job_number $total_jobs
done


permute="True"
if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
    # Stem of all output files
    output_stem=$qtl_results_dir$parameter_string"_permute_"$permute"_"$job_number"_"
    sbatch dynamic_qtl_shell.sh $joint_test_input_file $correction_factor_file $model_version $output_stem $permute $max_sites $job_number $total_jobs
done
fi


if false; then
# NEED TO CHANGE
parameter_string=$model_version"_environmental_variable_"$environmental_variable_form

sh multiple_testing_correction_and_visualization.sh $parameter_string $qtl_results_dir $target_region_input_file $qtl_visualization_dir $total_jobs $gencode_file

fi

