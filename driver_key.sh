#######################################################################################################
# This is the driver script that runs EAGLE models on ipsc data
#######################################################################################################



###############################################################################
# Input Data
###############################################################################

# Table of dimension (number of allelic sites) X (number of rna-seq samples)
# Each element is of form $numRefCounts"_"$numTotalCounts for that site in that sample
# It was created by 
allelic_counts_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_old_1_19_18/processed_allelic_counts/allelic_counts_gene_mapped_het_prob_999.txt"



###############################################################################
# Output directories (aasume all of these exist prior to starting analysis)
###############################################################################

# Root directory for this of all ipsc data based results
output_root="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data/"

# Directory containing text files with results from dynamic qtl analysis
qtl_results_dir=$output_root"qtl_results/"

# Directory containing visualization of results found in qtl_results_dir
qtl_visualization_dir=$output_root"qtl_visualization/"

# Directory containing gene set enrichment results
qtl_gene_set_enrichment_dir=$output_root"qtl_gene_set_enrichment/"

# Directory containing elbo diagnostic plots
elbo_diagnostic_dir=$output_root"elbo_diagnostic/"




###############################################################################
# Dynamic QTL Calling
###############################################################################

####################################
# Parameters
####################################
# Name of model. Current options are:
### 1. 'glm'
### 2. 'glmm'
model_version="glmm"
# Loci needs to have at least $min_het_lines heterozygous lines
min_het_lines="7"
# Restricting to only the heterozygous lines, at least $percent_biallelic % of the samples must be biallelic
percent_biallelic=".8"
# A sample must have at least $min_reads reads for it to be considered biallelic
min_reads="3"
# Boolean whether to take min transformation (either 'True' or 'False')
min_transformation="False"
# How many unique time steps to include in analysis
num_time_steps="16"

# String used in output files to keep track of parameters used
parameter_string=$model_version"_min_het_lines_"$min_het_lines"_percent_biallelic_"$percent_biallelic"_min_reads_"$min_reads"_min_transform_"$min_transformation"_num_time_steps_"$num_time_steps


####################################
# Run analysis
####################################

sh dynamic_qtl_run_and_visualize.sh $allelic_counts_file $qtl_results_dir $qtl_visualization_dir $parameter_string $model_version $min_het_lines $percent_biallelic $min_reads $min_transformation $num_time_steps $qtl_gene_set_enrichment_dir $elbo_diagnostic_dir


