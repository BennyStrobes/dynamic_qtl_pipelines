import numpy as np 
import os
import sys
import pdb

def get_cell_lines_from_sample_names(sample_names):
    cell_lines = []
    for sample_name in sample_names:
        cell_line = sample_name.split('_')[0]
        cell_lines.append(cell_line)
    return np.asarray(cell_lines)

# Compute how many unique cell lines are heterozygous at this locus
def get_number_of_heterozygous_lines(counts, cell_lines):
    het_cell_lines = []
    for i, count in enumerate(counts):
        if count == 'Nan':
            continue
        het_cell_lines.append(cell_lines[i])
    return len(np.unique(het_cell_lines))


def calculate_allelic_imbalence(ref_count, total_count):
    if total_count == 0:
        return .5
    else:
        return np.abs((float(ref_count)/total_count) - .5)

# Compute percent of samples that are heterozygous that show biallelic expression at this loci
def compute_percent_biallelic(counts, min_reads):
    # Keep track of which samples show biallelic expression and the total number of samples that are heterozygous
    biallelic_samples = 0
    total_samples = 0
    for count in counts:
        if count == 'Nan':
            continue
        total_samples = total_samples + 1
        ref_count = int(count.split('_')[0])
        total_count = int(count.split('_')[1])
        allelic_imbalence = calculate_allelic_imbalence(ref_count, total_count)
        if allelic_imbalence > .5:
            print('ASSUMPTION ERROROR')
        if total_count >= min_reads and allelic_imbalence < .49:
            biallelic_samples = biallelic_samples + 1
    return float(biallelic_samples)/total_samples


def transform_time_steps(time_steps, num_time_steps):
    bin_size = 16.0/num_time_steps
    time_steps = np.asarray(time_steps).astype(float)
    transformed_time_steps = np.floor(time_steps/np.ceil(bin_size)).astype(int)
    if max(transformed_time_steps) + 1 > num_time_steps:
        print('time step transformation error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        pdb.set_trace()
    return transformed_time_steps.astype(str)
# Now, this is a site that passes our filters
# We want to extract ordered list of:
### time steps
### cell lines
### alt counts
### total counts
def extract_ordered_list_of_samples_that_pass_filter(counts, min_reads, sample_names, min_transformation, num_time_steps):
    time_steps = []
    cell_lines = []
    alt_counts = []
    total_counts = []
    for i, count in enumerate(counts):
        if count == 'Nan':
            continue
        cell_line = sample_names[i].split('_')[0]
        time_step = sample_names[i].split('_')[1]
        ref_count = int(count.split('_')[0])
        total_count = int(count.split('_')[1])
        allelic_imbalence = calculate_allelic_imbalence(ref_count, total_count)
        if total_count >= min_reads and allelic_imbalence < .49:
            time_steps.append(time_step)
            cell_lines.append(cell_line)
            total_counts.append(str(total_count))
            if min_transformation == True:
                alt_count = min(ref_count, total_count - ref_count)
            elif min_transformation == False:
                alt_count = total_count - ref_count
            alt_counts.append(str(alt_count))
    transformed_time_steps = transform_time_steps(time_steps, num_time_steps)
    return np.asarray(transformed_time_steps), np.asarray(cell_lines), np.asarray(alt_counts), np.asarray(total_counts)


def convert_from_string_cell_lines_to_integer(line_cell_lines):
    converter = {}
    counter = 1
    for cell_line in line_cell_lines:
        if cell_line not in converter:
            converter[cell_line] = counter
            counter = counter + 1
    integer_cell_lines = []
    for cell_line in line_cell_lines:
        integer_cell_line = str(converter[cell_line])
        integer_cell_lines.append(integer_cell_line)
    return np.asarray(integer_cell_lines)

def prepare_input_file(allelic_counts_file, qtl_input_file, min_het_lines, percent_biallelic, min_reads, min_transformation, num_time_steps):
    # Open filehandle for input file
    f = open(allelic_counts_file)
    # Open filehandle for output file
    t = open(qtl_input_file, 'w')

    # Print header
    t.write('chrom_num\thet_position\tref_allele\talt_allele\trs_id\tgene_id\tordered_time_steps\tordered_cell_lines\tordered_cell_line_groups\tordered_alt_counts\tordered_total_counts\n')

    head_count = 0 # binary variable used to identify header
    # Loop through lines (het sites)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # This line is the header
            head_count = head_count + 1
            sample_names = np.asarray(data[1:])
            cell_lines = get_cell_lines_from_sample_names(sample_names)
            continue
        # Extract relevent info from line
        site_id = data[0]  # Name of heterozygous site
        counts = data[1:]  # Vector of allelic counts for each sample

        # Compute how many unique cell lines are heterozygous at this locus
        num_heterozygous_lines = get_number_of_heterozygous_lines(counts, cell_lines)
        # Ignore sites that do not have at least $min_het_lines unique cell lines
        if num_heterozygous_lines < min_het_lines:
            continue

        # Compute percent of samples that are heterozygous that show biallelic expression at this loci
        line_percent_biallelic = compute_percent_biallelic(counts, min_reads)
        # Ignore sites that do not have at least $percent biallelic % of samples
        if line_percent_biallelic < percent_biallelic:
            continue

        # Now, this is a site that passes our filters
        # We want to extract ordered list of:
        ### time steps
        ### cell lines
        ### alt counts
        ### total counts
        line_time_steps, line_cell_lines, line_alt_counts, line_total_counts = extract_ordered_list_of_samples_that_pass_filter(counts, min_reads, sample_names, min_transformation, num_time_steps)

        # Create ordered list of cell lines, but represent cell line as unique integers (base 1)
        # This is useful because this format is necessary for running EAGLE
        line_integer_cell_lines = convert_from_string_cell_lines_to_integer(line_cell_lines)


        # Extract relevent header info
        chrom_num = site_id.split('_')[0]
        position = site_id.split('_')[1]
        rs_id = site_id.split('_')[2]
        ref_allele = site_id.split('_')[3]
        alt_allele = site_id.split('_')[4]
        gene_id = site_id.split('_')[5]

        # Print to output file
        t.write(chrom_num + '\t' + position + '\t' + ref_allele + '\t' + alt_allele + '\t' + rs_id + '\t' + gene_id + '\t')
        t.write(','.join(line_time_steps) + '\t' + ','.join(line_cell_lines) + '\t' + ','.join(line_integer_cell_lines) + '\t' + ','.join(line_alt_counts) + '\t' + ','.join(line_total_counts) + '\n')
    t.close()



########################################################
# The goal of this script is to produce an input file for dynamic qtl calling
# containing only sites that pass the required filters
########################################################

########################################################
# Command line args
#########################################################
allelic_counts_file = sys.argv[1]
qtl_results_dir = sys.argv[2]
parameter_string = sys.argv[3]
min_het_lines = int(sys.argv[4])
percent_biallelic = float(sys.argv[5])
min_reads = int(sys.argv[6])
if sys.argv[7] == 'True':
    min_transformation = True
elif sys.argv[7] == 'False':
    min_transformation = False
else:
    print('min_transormation input error!')

num_time_steps = int(sys.argv[8])


# Dynamic qtl input file name (really the output file of this script)
qtl_input_file = qtl_results_dir + parameter_string + '_dynamic_qtl_input_file.txt'



prepare_input_file(allelic_counts_file, qtl_input_file, min_het_lines, percent_biallelic, min_reads, min_transformation, num_time_steps)