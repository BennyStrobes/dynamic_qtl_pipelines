#!/bin/bash
#SBATCH --time=20:00:00 --mem=10G --partition=broadwl

allelic_counts_file="$1"
qtl_results_dir="$2"
qtl_visualization_dir="$3"
parameter_string="$4"
model_version="$5"
min_het_lines="$6"
percent_biallelic="$7"
min_reads="$8"
min_transformation="$9"
num_time_steps="${10}"
qtl_gene_set_enrichment_dir="${11}"
elbo_diagnostic_dir="${12}"


if false; then
# Prepare dynamic qtl input file
python prepare_dynamic_qtl_input_file.py $allelic_counts_file $qtl_results_dir $parameter_string $min_het_lines $percent_biallelic $min_reads $min_transformation $num_time_steps
fi

# Run dyanmic qtl analysis
Rscript run_dynamic_qtl.R $qtl_results_dir $qtl_visualization_dir $parameter_string $model_version $elbo_diagnostic_dir


if false; then
python multiple_testing_correction.py $qtl_results_dir $parameter_string


Rscript visualize_dynamic_qtls.R $qtl_results_dir $qtl_visualization_dir $parameter_string



python run_gene_set_enrichment.py $qtl_results_dir $qtl_gene_set_enrichment_dir $parameter_string
fi