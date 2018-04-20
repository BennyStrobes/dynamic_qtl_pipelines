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
    cov = []
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
        cov.append(float(data[3]))
        # Add filehandle
        file_name = data[2]

        f = gzip.open(file_name)
        header = f.readline()  # Skip header
        filehandles.append(f)
    return sample_ids, filehandles, np.asarray(environmental_vars), np.asarray(cov)


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
        if test_info[10] == 'NA':  # Sample has no linked snps
            continue
        het_site_locations = test_info[9].split(';')
        linking_prob = np.asarray(test_info[11].split(';')).astype(float)
        het_probs = np.asarray(test_info[10].split(';')).astype(float)
        total_het_counts = np.asarray(test_info[12].split(';')).astype(int) + np.asarray(test_info[13].split(';')).astype(int)
        for het_site_index, het_prob in enumerate(het_probs):
            if het_prob > .9 and het_site_locations[het_site_index] not in mapping and total_het_counts[het_site_index] > 2 and linking_prob[het_site_index] > .9:
                mapping[het_site_locations[het_site_index]] = count
                count = count + 1

    # Now, using the mapping fill in ys and ns
    K = len(mapping)  # Number of het sites
    N = len(h_1)  # Number of samples
    # Initialize allelic count matrices
    ys = np.zeros((N,K))
    ns = np.zeros((N,K))
    for i,test_info in enumerate(test_infos):
        if test_info[10] == 'NA':  # Sample has no linked snps
            continue
        het_site_locations = test_info[9].split(';')
        het_probs = np.asarray(test_info[10].split(';')).astype(float)
        ref_counts = np.asarray(test_info[12].split(';')).astype(int)
        alt_counts = np.asarray(test_info[13].split(';')).astype(int)
        linking_prob = np.asarray(test_info[11].split(';')).astype(float)
        for het_site_index, het_prob in enumerate(het_probs):
            if het_prob > .9 and het_site_locations[het_site_index] in mapping and linking_prob[het_site_index] > .9:
                column_number = mapping[het_site_locations[het_site_index]]
                ns[i, column_number] = ref_counts[het_site_index] + alt_counts[het_site_index]
                if h_1[i] == 0:
                    ys[i, column_number] = ref_counts[het_site_index]
                else:
                    ys[i, column_number] = alt_counts[het_site_index]
    return ys.astype(int), ns.astype(int)


def extract_test_snp_dosage(test_infos):
    dosages = []
    # Loop through samples
    for test_info in test_infos:
        dosage = float(test_info[5])
        dosages.append(dosage)

        # Some basic error checking
        if dosage < 0.0 or dosage > 2.0:
            print('dosage assumption error!!')
            pdb.set_trace()
    return np.asarray(dosages)

# Re-organize test_infos data so that it is in a better format
# Specifically, we will extract a vector of length == Number_of_samples of:
### 1. gene_counts
### 2. allele_1 read counts (ys)
### 3. Sum of allele_1 and allele_2 read counts (ns)
### 4. Genotype on allele_1 (h_1)
### 5. Genotype on allele_2 (h_2)
def reorganize_data(test_infos):
    # Extract imputed dosage for this snp
    #dosage = extract_test_snp_dosage(test_infos)
    # Extract genotype on allele_1 and allele_2 (h_1 and h_2)
    h_1, h_2 = extract_test_snp_haplotype(test_infos)
    # Now extract allele specific read counts
    ys, ns = extract_allele_specific_read_counts(h_1, h_2, test_infos)
    return ys, ns, h_1, h_2

# For parrallelization purposes
# start_test is test_number we start at. end_test is test_number we end at
def determine_start_test_and_end_test(num_tests, total_jobs, job_number):
    bin_size = np.ceil(num_tests/total_jobs)
    start = bin_size*job_number
    end = bin_size*(job_number+1)
    return start, end

# Create dictionary the maps cell_line ids to an integer array. Where the array contains indices of each of the samples in that line
def get_cell_line_indices(joint_test_input_file):
    f = open(joint_test_input_file)
    dicti = {} # Initialize dictionary
    head_count = 0  # used to skip header
    index_counter = 0  # used to keep track of sample index
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # skip header
            head_count = head_count + 1
            continue
        sample_id = data[0]
        # extract cell line from sample id
        cell_line = sample_id.split('_')[0]
        # If cell line has never been seen before, add to dictionary
        if cell_line not in dicti:
            dicti[cell_line] = []
        # Add index to cell line specific array
        dicti[cell_line].append(index_counter)
        index_counter = index_counter + 1
    for key in dicti:
        dicti[key] = np.asarray(dicti[key])
    return dicti

# Make binary matrix where each element is True if site, sample pair is biallelic, and false if not
def construct_binary_biallelic_matrix(ys, ns, min_total_reads, biallelic_thresh):
    num_samples = ys.shape[0]
    num_sites = ys.shape[1]
    # Initialize output matrix
    binary_mat = np.ones((num_samples, num_sites), dtype=bool)
    # Loop through samples
    for n in range(num_samples):
        # Loop through allelic sites
        for k in range(num_sites):
            if ns[n, k] < min_total_reads:
                binary_mat[n, k] = False
            else:
                biallelic_fraction = np.abs((ys[n, k]/(ns[n, k])) - 0.5)
                if biallelic_fraction >= biallelic_thresh:
                    binary_mat[n, k] = False
    return binary_mat

# Extract number of samples that are biallelic and have a heterozygous test_variant
def get_num_het_test_variant_biallelic_samples(hets, biallelic_samples):
    count = 0
    num_samples = len(hets)
    for n in range(num_samples):
        if hets[n] == True and biallelic_samples[n] == True:
            count = count + 1
    return count


def get_num_biallelic_lines(biallelic_samples, cell_line_indices):
    biallelic_lines = {}
    for cell_line in cell_line_indices.keys():
        indices = cell_line_indices[cell_line]
        if np.sum(biallelic_samples[indices]) > 0:
            biallelic_lines[cell_line] = 1
    return len(biallelic_lines)


# If the dimension of ys and ns (Ie. the number of heterozygous exonic sites) is too large, the algorithm will become prohibitively slow
# So we subset to the top $max_sites with the largest number of minor allele read counts
def subset_allelic_count_matrices(ys, ns, cell_line_indices, min_num_biallelic_lines, min_num_biallelic_samples, min_num_het_test_variant_biallelic_samples, h_1, h_2):
    num_samples = ys.shape[0]
    num_sites = ys.shape[1]
    # Make binary matrix where each element is True if site, sample pair is biallelic, and false if not
    binary_biallelic_matrix = construct_binary_biallelic_matrix(ys, ns, 3, .49)
    # Initialize vector length of number of sites. Take on value true if site passes filters
    keep_sites = []
    for k in range(num_sites):
        # Number of biallic samples for this site
        num_biallelic_samples = np.sum(binary_biallelic_matrix[:, k])
        # Extract number of samples that are biallelic and have a heterozygous test_variant
        num_het_test_variant_biallelic_samples = get_num_het_test_variant_biallelic_samples(h_1 != h_2, binary_biallelic_matrix[:, k])
        # Extract number of lines that have a biallelic samples
        num_biallelic_lines = get_num_biallelic_lines(binary_biallelic_matrix[:, k], cell_line_indices)

        # Check if index passes filters
        if num_biallelic_samples >= min_num_biallelic_samples and num_het_test_variant_biallelic_samples >= min_num_het_test_variant_biallelic_samples and num_biallelic_lines >= min_num_biallelic_lines:
            keep_sites.append(True)
        else:
            keep_sites.append(False)
    ys_filtered = ys[:, keep_sites]
    ns_filtered = ns[:, keep_sites]

    return ys_filtered, ns_filtered



# Main Driver
def run_analysis(joint_test_input_file, model_version, output_stem, job_number, total_jobs, permute, optimization_method, permutation_scheme, min_num_biallelic_lines, min_num_biallelic_samples, min_num_het_test_variant_biallelic_samples, covariate_method):
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
    sample_ids, filehandles, environmental_vars, covs = parse_joint_test_input_file(joint_test_input_file)

    # Create dictionary the maps cell_line ids to an integer array. Where the array contains indices of each of the samples in that cell line
    cell_line_indices = get_cell_line_indices(joint_test_input_file)

    t = open(output_stem + 'dynamic_qtl_results.txt', 'w')
    # Loop through each of the tests
    for test_number in range(num_tests):
        # Extract one line from each file handle and put into ordered list
        test_infos = extract_test_info_from_filehandles(filehandles)

        # Only perform tests corresponding to this job number (for parrallelization purposes)
        if test_number < start_test or test_number >= end_test:
            continue
        print(test_number - start_test)
        # Check to make sure the lines from each sample are all testing the same test ((variant, gene) pair)
        check_to_make_sure_tests_are_all_the_same(test_infos)
        # Re-organize test_infos data so that it is in a better format
        # Specifically, we will extract a vector of length == Number_of_samples of:
        ### 1. gene_counts
        ### 2. allele_1 read counts (ys)
        ### 3. Sum of allele_1 and allele_2 read counts (ns)
        ### 4. Genotype on allele_1 (h_1)
        ### 5. Genotype on allele_2 (h_2)
        ys, ns, h_1, h_2 = reorganize_data(test_infos)

        # If the dimension of ys and ns (Ie. the number of heterozygous exonic sites) is too large, the algorithm will become prohibitively slow
        # So we subset to the top $max_sites with the largest number of minor allele read counts
        ys, ns = subset_allelic_count_matrices(ys, ns, cell_line_indices, min_num_biallelic_lines, min_num_biallelic_samples, min_num_het_test_variant_biallelic_samples, h_1, h_2)

        # Throw out some sites in allele specific only model
        if ys.shape[1] == 0 and model_version == 'as_log_linear':
            print('throw out :' + str(test_number))
            continue
        elif ys.shape[1] == 0 and model_version == 'as_log_linear_fixed_overdispersion':
            print('throw out :' + str(test_number))
            continue
        print(test_infos[0][2] + '\t' + test_infos[0][7])
        # Run dynamic qtl test
        result = dynamic_qtl(ys, ns, h_1, h_2, environmental_vars, model_version, permute, optimization_method, cell_line_indices, permutation_scheme, covs, covariate_method)
        # Print results to output file
        t.write('\t'.join(test_infos[0][0:5]) + '\t' + test_infos[0][7] + '\t' + test_infos[0][8] + '\t' + str(result['loglike_null']) + '\t' + str(result['loglike_full']) + '\t')
        # t.write(test_infos[0][0] + '\t' +  test_infos[0][1] + '\t' + test_infos[0][2].split("'")[1] + '\t' + test_infos[0][3].split("'")[1] + '\t' + test_infos[0][4].split("'")[1] + '\t' + test_infos[0][7] + '\t' + test_infos[0][8] + '\t' + str(result['loglike_null']) + '\t' + str(result['loglike_full']) + '\t')
        if model_version == 'as_log_linear':
            t.write(str(result['loglr']) + '\t' + str(result['pvalue']) + '\t' + str(result['fit_full']['par']['beta'][-1]) + '\t' + ','.join(np.atleast_1d(result['fit_full']['par']['conc']).astype(str)) + '\n')
        if np.mod(test_number, 10) == 0:
            t.flush()
    t.close()


#######################
# Command line args
#######################
joint_test_input_file = sys.argv[1]  # Ordered list of samples and their corresponding cht input files
model_version = sys.argv[2]  # Name of glm that we are going to use
output_stem = sys.argv[3]  # Stem of output file
permute = sys.argv[4]  # Binary (True/False) variable on whehter to permute the data or not
job_number = int(sys.argv[5])
total_jobs = int(sys.argv[6])
optimization_method = sys.argv[7]  # Technique to fit the GLM
permutation_scheme = sys.argv[8]  # How to permute NUll data (only applies if permute == 'True')
min_num_biallelic_lines = int(sys.argv[9])  # Minimum number of cell lines that have at least one biallelic sample required for a site to be used
min_num_biallelic_samples = int(sys.argv[10])  #  Minimum number of biallelic samples required for a site to be used
min_num_het_test_variant_biallelic_samples = int(sys.argv[11])  # Minimum number of biallelic samples that have a heterozygous test variant required for a site to be used
covariate_method = sys.argv[12]


# Main Driver
run_analysis(joint_test_input_file, model_version, output_stem, job_number, total_jobs, permute, optimization_method, permutation_scheme, min_num_biallelic_lines, min_num_biallelic_samples, min_num_het_test_variant_biallelic_samples, covariate_method)
