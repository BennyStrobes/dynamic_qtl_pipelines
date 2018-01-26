import numpy as np 
import os
import pdb
import sys
import gzip
from dynamic_qtl_driver import dynamic_qtl



# Parse through joint_test_input_file to extract ordered list of:
## 1. Sample Ids
## 2. Filehandles opening input cht files
## 3. Environmental variables
def parse_joint_test_input_file(joint_test_input_file):
    # Initialize output vectors
    sample_ids = []
    filehandles = []
    environmental_vars = []
    # Open input filehandle
    f = open(joint_test_input_file)
    head_count = 0  # skip header
    # Loop through each sample
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        # Add sample id and environmental variable to array
        sample_id = data[0]
        environmental_var = float(data[1])
        sample_ids.append(sample_id)
        environmental_vars.append(environmental_var)
        # Add filehandle
        file_name = data[2]
        f = gzip.open(file_name)
        header = f.readline()  # Skip header
        filehandles.append(f)
    return sample_ids, filehandles, environmental_vars

# Parse through correction_factor_file to extract ordered list of:
## 1. library_size_correction_factors
def parse_correction_factor_file(correction_factor_file):
    # Initialize output array
    correction_factors = []
    head_count = 0  # Used to skip header
    f = open(correction_factor_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        correction_factor = float(data[1])
        correction_factors.append(correction_factor)
    return np.asarray(correction_factors)

def extract_number_of_tests(joint_test_input_file):
    # Extract the cht input filename for one sample
    f = open(joint_test_input_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        file_name = data[2]
    f.close()
    # Count how many lines are in this cht sample
    f = gzip.open(file_name)
    counter = 0
    for line in f:
        counter = counter + 1
    f.close()
    return counter - 1


# Extract one line from each file handle and put into ordered list
def extract_test_info_from_filehandles(filehandles):
    test_infos = []  # Initialize output list
    # Loop through filehandles (samples)
    for filehandle in filehandles:
        line = filehandle.readline().decode('utf-8').rstrip()
        test_infos.append(line.split())
    return test_infos

# Convert full test line to an identifier that is specific to the variant, gene pair and not the sample
def test_info_to_test_id(info):
    return info[0] + '_' + info[1] + '_' + info[2] + '_' + info[3] + '_' + info[4] + '_' + info[7] + '_' + info[8]

# Check to make sure the lines from each sample are all testing the same test ((variant, gene) pair)
def check_to_make_sure_tests_are_all_the_same(test_infos):
    # Convert full test line to an identifier that is specific to the variant, gene pair and not the sample
    test_id = test_info_to_test_id(test_infos[0])
    for test_info in test_infos:
        if test_id != test_info_to_test_id(test_info):
            print('ASSUMTPTION ERROR, must stop')
            pdb.set_trace()
    return

# First extract gene_counts from test_infos
def extract_gene_counts(test_infos):
    # Initialize output vector
    gene_counts = []
    # Loop through samples
    for test_info in test_infos:
        gene_count = int(test_info[15])
        gene_counts.append(gene_count)
    return np.asarray(gene_counts)

# Extract genotype on allele_1 and allele_2 (h_1 and h_2)
def extract_test_snp_haplotype(test_infos):
    # Initialize output vector
    h_1 = []
    h_2 = []
    # Loop through samples
    for test_info in test_infos:
        haplotype = test_info[6]
        hap_allele_1 = int(haplotype.split('|')[0])
        hap_allele_2 = int(haplotype.split('|')[1])
        h_1.append(hap_allele_1)
        h_2.append(hap_allele_2)
    return np.asarray(h_1), np.asarray(h_2)

# Now extract allele specific read counts
def extract_allele_specific_read_counts(h_1, h_2, test_infos):
    # Create mapping from het_site_id to column in ys and ns
    mapping = {}
    count = 0
    for i,test_info in enumerate(test_infos):
        if h_1[i] == h_2[i]:  # test_snp is homozygous in this sample, so do not consider for allelic imbalence
            continue
        if test_info[10] == 'NA':  # Sample has no linked snps
            continue
        het_site_locations = test_info[9].split(';')
        het_probs = np.asarray(test_info[10].split(';')).astype(float)
        total_het_counts = np.asarray(test_info[12].split(';')).astype(int) + np.asarray(test_info[13].split(';')).astype(int)
        for het_site_index, het_prob in enumerate(het_probs):
            if het_prob > .9 and het_site_locations[het_site_index] not in mapping and total_het_counts[het_site_index] > 2:
                mapping[het_site_locations[het_site_index]] = count
                count = count + 1

    # Now, using the mapping fill in ys and ns
    K = len(mapping)  # Number of het sites
    N = len(h_1)  # Number of samples
    # Initialize allelic count matrices
    ys = np.zeros((N,K))
    ns = np.zeros((N,K))
    for i,test_info in enumerate(test_infos):
        if h_1[i] == h_2[i]:  # test_snp is homozygous in this sample, so do not consider for allelic imbalence
            continue
        if test_info[10] == 'NA':  # Sample has no linked snps
            continue
        het_site_locations = test_info[9].split(';')
        het_probs = np.asarray(test_info[10].split(';')).astype(float)
        ref_counts = np.asarray(test_info[12].split(';')).astype(int)
        alt_counts = np.asarray(test_info[13].split(';')).astype(int)
        for het_site_index, het_prob in enumerate(het_probs):
            if het_prob > .9 and het_site_locations[het_site_index] in mapping:
                column_number = mapping[het_site_locations[het_site_index]]
                ns[i, column_number] = ref_counts[het_site_index] + alt_counts[het_site_index]
                if h_1[i] == 0:
                    ys[i, column_number] = ref_counts[het_site_index]
                else:
                    ys[i, column_number] = alt_counts[het_site_index]
    return ys.astype(int), ns.astype(int)


# Re-organize test_infos data so that it is in a better format
# Specifically, we will extract a vector of length == Number_of_samples of:
### 1. gene_counts
### 2. allele_1 read counts (ys)
### 3. Sum of allele_1 and allele_2 read counts (ns)
### 4. Genotype on allele_1 (h_1)
### 5. Genotype on allele_2 (h_2)
def reorganize_data(test_infos):
    # First extract gene_counts from test_infos
    gene_counts = extract_gene_counts(test_infos)
    # Extract genotype on allele_1 and allele_2 (h_1 and h_2)
    h_1, h_2 = extract_test_snp_haplotype(test_infos)
    # Now extract allele specific read counts
    ys, ns = extract_allele_specific_read_counts(h_1, h_2, test_infos)
    return gene_counts, ys, ns, h_1, h_2

# For parrallelization purposes
# start_test is test_number we start at. end_test is test_number we end at
def determine_start_test_and_end_test(num_tests, total_jobs, job_number):
    bin_size = np.ceil(num_tests/total_jobs)
    start = bin_size*job_number
    end = bin_size*(job_number+1)
    return start, end

# Main Driver
def run_analysis(joint_test_input_file, correction_factor_file, model_version, output_stem, job_number, total_jobs):
    # Determine the number of tests
    num_tests = extract_number_of_tests(joint_test_input_file)
    # For parrallelization purposes
    # start_test is test_number we start at. end_test is test_number we end at
    start_test, end_test = determine_start_test_and_end_test(num_tests, total_jobs, job_number)
    # First extract ordered list of:
    ## 1. Sample Ids
    ## 2. Filehandles opening input cht files
    ## 3. Environmental variables
    ## 4. library size correction factors
    sample_ids, filehandles, environmental_vars = parse_joint_test_input_file(joint_test_input_file)
    library_size_correction_factors = parse_correction_factor_file(correction_factor_file)

    t = open(output_stem + 'dynamic_qtl_results.txt', 'w')
    # Loop through each of the tests
    for test_number in range(num_tests):
        if np.mod(test_number, 10000) == 0:
            print(test_number)
        # Extract one line from each file handle and put into ordered list
        test_infos = extract_test_info_from_filehandles(filehandles)
        # Only perform tests corresponding to this job number (for parrallelization purposes)
        if test_number < start_test or test_number >= end_test:
            continue
        # Check to make sure the lines from each sample are all testing the same test ((variant, gene) pair)
        check_to_make_sure_tests_are_all_the_same(test_infos)
        # Re-organize test_infos data so that it is in a better format
        # Specifically, we will extract a vector of length == Number_of_samples of:
        ### 1. gene_counts
        ### 2. allele_1 read counts (ys)
        ### 3. Sum of allele_1 and allele_2 read counts (ns)
        ### 4. Genotype on allele_1 (h_1)
        ### 5. Genotype on allele_2 (h_2)
        gene_counts, ys, ns, h_1, h_2 = reorganize_data(test_infos)

        # Run dynamic qtl test
        result = dynamic_qtl(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors, model_version)
        # Print results to output file
        t.write('\t'.join(test_infos[0][0:5]) + '\t' + test_infos[0][7] + '\t' + test_infos[0][8] + '\t' + str(result['loglike_null']) + '\t' + str(result['loglike_full']) + '\t')
        t.write(str(result['loglr']) + '\t' + str(result['pvalue']) + '\t' + str(result['fit_full']['par']['beta'][-1]) + '\n')
        if np.mod(test_number, 100) == 0:
            t.flush()
    t.close()


#######################
# Command line args
#######################
joint_test_input_file = sys.argv[1]  # Ordered list of samples and their corresponding cht input files
correction_factor_file = sys.argv[2]  # Ordered list of samples (same order as joint_test_input_file) that contains the library_size_correction_factor for each model_version
model_version = sys.argv[3]  # Name of glm that we are going to use
output_stem = sys.argv[4]  # Stem of output file
job_number = int(sys.argv[5])
total_jobs = int(sys.argv[6])


# Main Driver
run_analysis(joint_test_input_file, correction_factor_file, model_version, output_stem, job_number, total_jobs)