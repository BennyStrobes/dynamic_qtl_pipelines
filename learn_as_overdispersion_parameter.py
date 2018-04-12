import numpy as np
import os
import sys
import pdb
import pystan
import gzip


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
    return sample_ids, filehandles, np.asarray(environmental_vars)

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

# Re-organize test_infos data so that it is in a better format
# Specifically, we will extract a vector of length == Number_of_samples of:
### 1. allele_1 read counts (ys)
### 2. Sum of allele_1 and allele_2 read counts (ns)
### 3. Genotype on allele_1 (h_1)
### 4. Genotype on allele_2 (h_2)
def reorganize_data(test_infos):
    h_1, h_2 = extract_test_snp_haplotype(test_infos)
    # Now extract allele specific read counts
    ys, ns = extract_allele_specific_read_counts(h_1, h_2, test_infos)
    return ys, ns, h_1, h_2

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


def add_test_to_aggregrated_count_list(ys, ns, total_counts, allelic_counts, sample_names):
    num_sites = ys.shape[1]
    num_samples = ys.shape[0]
    for sample_num in range(num_samples):
        for site_num in range(num_sites):
            # Don't include samples with zero counts
            if ns[sample_num, site_num] == 0:
                continue
            total_counts.append(ns[sample_num, site_num])
            allelic_counts.append(ys[sample_num, site_num])
            sample_names.append(sample_num)
    return total_counts, allelic_counts, sample_names


def run_analysis(joint_test_input_file, as_overdispersion_parameter_file, as_overdispersion_parameter_sample_specific_file, min_num_biallelic_lines, min_num_biallelic_samples, min_num_het_test_variant_biallelic_samples):
    # Determine the number of tests
    num_tests = extract_number_of_tests(joint_test_input_file)
    # First extract ordered list of:
    ## 1. Sample Ids
    ## 2. Filehandles opening input cht files
    ## 3. Environmental variables
    ## 4. library size correction factors
    sample_ids, filehandles, environmental_vars = parse_joint_test_input_file(joint_test_input_file)
    # Loop through each of the tests


    # Create dictionary the maps cell_line ids to an integer array. Where the array contains indices of each of the samples in that cell line
    cell_line_indices = get_cell_line_indices(joint_test_input_file)

    # Initialize vectors to keep track of allelic counts and sample ids
    total_counts = []
    allelic_counts = []
    sample_names = []

    for test_number in range(num_tests):
        # Extract one line from each file handle and put into ordered list
        test_infos = extract_test_info_from_filehandles(filehandles)
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
        
        # Skip tests that don't have any allelic sites
        if ys.shape[1] == 0:
            continue

        # Add allelic counts from this test to our genome wide list of counts
        total_counts, allelic_counts, sample_names = add_test_to_aggregrated_count_list(ys, ns, total_counts, allelic_counts, sample_names)
    
    # Convert genome wide arrays to numpy vectors
    total_counts = np.transpose(np.asmatrix(np.asarray(total_counts)))
    allelic_counts = np.transpose(np.asmatrix(np.asarray(allelic_counts)))
    sample_names = np.asarray(sample_names) + 1
    num_samples = len(np.unique(sample_names))

    # randomly subset samples (select 1/70th of the original samples) for computational traction 
    # This seems like we are throwing away a lot of samples, but really does not affect our parameter estimation
    index_subset = np.random.choice(len(sample_names), int(np.floor(len(sample_names)/70)), replace=False)
    total_counts = total_counts[index_subset,:]
    allelic_counts = allelic_counts[index_subset, :]
    sample_names = sample_names[index_subset]


    # ONE OVERDISPERSION PARAMETER SHARED ACROSS ALL SAMPLES
    # Load in pystan object
    sm = pystan.StanModel(file='as_overdispersion_parameter.stan')

    # Get data into corect format for pystan
    data = dict(N=len(total_counts), ns=total_counts, ys=allelic_counts, concShape=1.01, concRate=0.01)
    # Run pystan optimization
    op = sm.optimizing(data=data)
    conc = np.atleast_1d(op['conc'])[0]

    t = open(as_overdispersion_parameter_file, 'w')
    t.write(str(conc))
    t.flush()
    t.close()

    # ONE OVERDISPERSION PARAMETER per sample
    # Load in pystan object
    sm = pystan.StanModel(file='as_overdispersion_parameter_sample_specific.stan')

    # Get data into corect format for pystan
    data = dict(N=len(total_counts), K=num_samples, ns=total_counts, ys=allelic_counts, sample_names=sample_names, concShape=1.01, concRate=0.01)
    # Run pystan optimization
    op = sm.optimizing(data=data)

    t = open(as_overdispersion_parameter_sample_specific_file, 'w')
    t.write('\t'.join(op['conc'].astype(str)))
    t.flush()
    t.close()




joint_test_input_file = sys.argv[1]
as_overdispersion_parameter_file = sys.argv[2]
as_overdispersion_parameter_sample_specific_file = sys.argv[3]
min_num_biallelic_lines = int(sys.argv[4])
min_num_biallelic_samples = int(sys.argv[5])
min_num_het_test_variant_biallelic_samples = int(sys.argv[6])


run_analysis(joint_test_input_file, as_overdispersion_parameter_file, as_overdispersion_parameter_sample_specific_file, min_num_biallelic_lines, min_num_biallelic_samples, min_num_het_test_variant_biallelic_samples)