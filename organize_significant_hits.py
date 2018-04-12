import numpy as np 
import os
import pdb
import sys
import gzip



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


# Parse through joint_test_input_file to extract ordered list of:
## 1. Sample Ids
## 2. Filehandles opening input cht files
## 3. Environmental variables
def parse_joint_test_input_file(joint_test_input_file):
    # Initialize output vectors
    sample_ids = []
    filehandles = []
    environmental_vars = []
    cell_lines = []
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
        cell_lines.append(sample_id.split('_')[0])
        # Add filehandle
        file_name = data[2]
        f = gzip.open(file_name)
        header = f.readline()  # Skip header
        filehandles.append(f)
    return sample_ids, filehandles, np.asarray(environmental_vars), np.asarray(cell_lines)

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


def extract_significant_variant_gene_pairs(egenes_file):
    f = open(egenes_file)
    head_count = 0
    pairs = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        rs_id = data[2]
        gene_id = data[5]
        pairs[rs_id + '_' + gene_id] = 1
    return pairs

def extract_null_variant_gene_pairs(file_name):
    f = open(file_name)
    pairs = {}
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        pvalue = float(data[-3])
        if pvalue == 0:
            rs_id = data[2]
            ensamble_id = data[5]
            pairs[rs_id + '_' + ensamble_id] = 1
    return pairs


############################################
# Extract vector of length number of tests where each element is a tuple
# The first element in the tuple is a binary variable that is a 1 if that variant-gene pair is a highest
# The second element of the tuple is a list of information describing the variant gene pair
def extract_hits_vector(egenes_file, hits_file):
    significant_variant_gene_pairs = extract_significant_variant_gene_pairs(egenes_file)
    hits_vector = []  # Initialize output vector
    head_count = 0  # Skip header
    used_genes = {}
    counter = 0
    hits = {}
    f = open(hits_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        ensamble_id = data[5]
        rs_id = data[2]
        # Extract relevent information from the current line (current line is 1 variant gene pair)
        pvalue = float(data[-3])
        if rs_id + '_' + ensamble_id in significant_variant_gene_pairs:  # Hit
            counter = counter + 1
            # hit information
            dicti = {}
            dicti['chrom_num'] = data[0]
            dicti['variant_pos'] = data[1]
            dicti['rs_id'] = data[2]
            dicti['ref_allele'] = data[3]
            dicti['alt_allele'] = data[4]
            dicti['ensamble_id'] = data[5]
            dicti['pvalue'] = pvalue
            dicti['beta'] = float(data[-2])
            dicti['conc'] = np.asarray(data[-1].split(','))
            hits[rs_id + '_' + ensamble_id] = dicti
    print(len(hits))
    return hits


def extract_hits_vector_v2(all_hits_file):
    significant_variant_gene_pairs = {}
    head_count = 0
    f = open(all_hits_file)
    hits_vector = []  # Initialize output vector
    used_genes = {}
    hits = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        pvalue = float(data[-3])
        ensamble_id = data[5]
        rs_id = data[2]
        if pvalue > 1e-7 or ensamble_id in used_genes:
            continue
        # hit information
        dicti = {}
        dicti['chrom_num'] = data[0]
        dicti['variant_pos'] = data[1]
        dicti['rs_id'] = data[2]
        dicti['ref_allele'] = data[3]
        dicti['alt_allele'] = data[4]
        dicti['ensamble_id'] = data[5]
        dicti['pvalue'] = pvalue
        dicti['beta'] = data[-2]
        dicti['conc'] = data[-1]
        hits[rs_id + '_' + ensamble_id] = dicti
        used_genes[ensamble_id] = 1
    print(len(hits))
    return hits




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



# Convert from dosage vector to vector of genotyeps
def dosage_to_genotype(dosage, ref_allele, alt_allele):
    converter = {}
    #converter[0] = ref_allele + ref_allele
    #converter[1] = ref_allele + alt_allele
    #converter[2] = alt_allele + alt_allele
    converter[0] = "0"
    converter[1] = "1"
    converter[2] = "2"


    genotype = []
    for dos in dosage:
        genotype.append(converter[dos])
    return np.asarray(genotype)

def visualize_hit(ys, ns, gene_counts, h_1, h_2, environmental_vars, library_size_correction_factors, test_dicti, output_file):
    plt.clf()
    # Plot to visualize total expression changes over time as a function of genotype
    fig = plt.figure(figsize=(40, 20))

    # call regplot on each axes
    ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=4)
    ax2 = plt.subplot2grid((2, 4), (1, 0))
    ax3 = plt.subplot2grid((2, 4), (1, 1))
    ax4 = plt.subplot2grid((2, 4), (1, 2))
    ax5 = plt.subplot2grid((2, 4), (1, 3))
    #sns.regplot(x=idx, y=df['x'], ax=ax1)
    #sns.regplot(x=idx, y=df['y'], ax=ax2)
    num_sites = ys.shape[1]
    if num_sites > 4:
        num_sites = 4

    gene_total_plot = gene_total_plotter(gene_counts, h_1+h_2, environmental_vars, test_dicti['rs_id'], test_dicti['ensamble_id'], test_dicti['ref_allele'], test_dicti['alt_allele'], test_dicti['pvalue'], test_dicti['beta'], library_size_correction_factors,ax1)
    #fig.savefig(output_file)
    num_sites = ys.shape[1]
    if num_sites > 4:
        num_sites = 4
    for exonic_site_num in np.arange(num_sites):
        if exonic_site_num == 0:
            axy = ax2
        elif exonic_site_num == 1:
            axy = ax3
        elif exonic_site_num == 2:
            axy = ax4
        elif exonic_site_num == 3:
            axy = ax5
        allelic_imbalence_plot = allelic_imbalence_plotter(h_1, h_2, environmental_vars, ys, ns, test_dicti['rs_id'], test_dicti['ensamble_id'], test_dicti['pvalue'], test_dicti['beta'], test_dicti['conc'].astype(float),axy, exonic_site_num)
    plt.tight_layout()
    fig.savefig(output_file)


# Convert from dosage vector to vector of genotyeps
def dosage_to_genotype(dosage, ref_allele, alt_allele):
    converter = {}
    #converter[0] = ref_allele + ref_allele
    #converter[1] = ref_allele + alt_allele
    #converter[2] = alt_allele + alt_allele
    converter[0] = "0"
    converter[1] = "1"
    converter[2] = "2"


    genotype = []
    for dos in dosage:
        genotype.append(converter[dos])
    return np.asarray(genotype)

# Plot to visualize total expression changes over time as a function of genotype
def gene_total_plotter(gene_counts, dosage, environmental_vars, rs_id, ensamble_id, ref_allele, alt_allele, pvalue, beta, library_size_correction_factors,ax1):
    # Convert from dosage vector to vector of genotyeps
    genotype = dosage_to_genotype(dosage, ref_allele, alt_allele)

    # gene_counts = (gene_counts/library_size_correction_factors)*np.mean(library_size_correction_factors)

    df = pd.DataFrame({rs_id: genotype, 'time_step': environmental_vars.astype(int), 'gene_counts': np.log(gene_counts)})
    ax = sns.boxplot(x="time_step", y="gene_counts", hue=rs_id, data=df, palette="Set3",width=.7,ax=ax1)
    plt.xlabel('Time Step')
    plt.ylabel('log(counts)')
    ax1.set_title(ensamble_id + ' / pvalue = ' + str(pvalue) + ' / beta = ' + str(beta))
    #sns.despine(offset=1, trim=True)
    return ax

def get_allelic_counts(ys, ns, h_1, h_2, site_num, environmental_vars):
    num_samples = ys.shape[0]
    a1_counts = []
    a2_counts = []
    environmental_vars_subset = []
    for sample_num in np.arange(num_samples):
        if h_1[sample_num] == h_2[sample_num]:  # homozygous variant
            continue
        if ns[sample_num, site_num] < 1:
            continue
        if h_1[sample_num] == 0:
            a1_counts.append(ys[sample_num, site_num])
            a2_counts.append(ns[sample_num, site_num] - ys[sample_num, site_num])
            environmental_vars_subset.append(environmental_vars[sample_num])
        elif h_2[sample_num] == 0:
            a2_counts.append(ys[sample_num, site_num])
            a1_counts.append(ns[sample_num, site_num] - ys[sample_num, site_num])
            environmental_vars_subset.append(environmental_vars[sample_num])
    return np.asarray(a1_counts), np.asarray(a2_counts), np.asarray(environmental_vars_subset)


def allelic_imbalence_plotter(h_1, h_2, environmental_vars, ys, ns, rs_id, ensamble_id, pvalue, beta, conc,axy, exonic_site_num):
    if ys.shape[1] != len(conc):
        print('fatal errooror')
        pdb.set_trace()
    # Order sites by smallest variance to largest
    ys = ys[:, conc.argsort()[::-1]]
    ns = ns[:, conc.argsort()[::-1]]
    # Extract allelic fractions at each site
    time_steps = []
    allelic_fractions = []
    cell_lines_arr = []
    identifiers = []
    exonic_sites = []
    depths = []
    num_exonic_sites = ys.shape[1]
    num_samples = len(environmental_vars)
    for sample_num in np.arange(num_samples):
        if h_1[sample_num] == h_2[sample_num]:  # homozygous variant
            continue
        if ns[sample_num, exonic_site_num] <= 2:
            continue
        if h_1[sample_num] == 0:
            allelic_fraction = 1.0 - float(ys[sample_num, exonic_site_num])/ns[sample_num, exonic_site_num]
        elif h_2[sample_num] == 0:
            allelic_fraction = (float(ys[sample_num, exonic_site_num])/ns[sample_num, exonic_site_num])
        else:
            print('eroroororo')
        depths.append(ns[sample_num,exonic_site_num])
        allelic_fractions.append(allelic_fraction)
        #allelic_fractions.append(abs(allelic_fraction-.5))
        time_steps.append(float(environmental_vars[sample_num]))

    # PLOT!
    df = pd.DataFrame({'time_step': np.asarray(time_steps).astype(int),'read_depth':depths,  'allelic_fraction': allelic_fractions})

    #ax = sns.pointplot(x="time_step", y="allelic_fraction", hue="identifiers", data=df)
    #ax = sns.regplot(x="time_step", y="allelic_fraction", data=df)
    ax = sns.regplot(data=df,x="time_step", y="allelic_fraction",ci=None, ax=axy)
    ax.set_title("Exonic site = " + str(exonic_site_num) + " / conc = " + str(conc[conc.argsort()[::-1]][exonic_site_num]))
    plt.ylim(ymax=1) # adjust the max leaving min unchanged
    plt.ylim(ymin=0)
    #ax = sns.boxplot(x="time_step", y="allelic_fraction", hue=exonic_sites,data=df, palette="Set3",width=.7)
    #plt.xlabel('Time Step')
    #plt.ylabel('Allelic fraction')
    #plt.title(ensamble_id + ' / pvalue = ' + str(pvalue) + ' / beta = ' + str(beta))
    #sns.despine(offset=1, trim=True)
    #iris = sns.load_dataset("iris")
    # Plot tip as a function of toal bill across days
    # = sns.regplot(x="sepal_length", y="sepal_width", data=iris)
    return ax



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


# Make mapping from gene identifier to ensamble id
def make_gene_mapping(target_region_input_file):
    # Get input for only 1 sample 
    # Only need one because genes are the same in all sample input files
    f = open(target_region_input_file)
    gene_mapping = {}  # Initialize gene mapping
    for line in f:
        line = line.rstrip()
        data = line.split()
        ensamble_id = data[6].split('_')[-1]
        gene_identifier = data[0] + '_' + data[7] + '_' + data[8]
        gene_mapping[gene_identifier] = ensamble_id
    f.close()
    return gene_mapping


###########################
# Command Line args
###########################
parameter_string = sys.argv[1]
qtl_results_dir = sys.argv[2]
joint_test_input_file = sys.argv[3]
correction_factor_file = sys.argv[4]
qtl_visualization_dir = sys.argv[5]
min_num_biallelic_lines = int(sys.argv[6])
min_num_biallelic_samples = int(sys.argv[7])
min_num_het_test_variant_biallelic_samples = int(sys.argv[8])
short_parameter_string = sys.argv[9]
target_region_input_file = sys.argv[10]

############################################
# Extract vector of length number of tests where each element is a tuple
# The first element in the tuple is a binary variable that is a 1 if that variant-gene pair is a highest
# The second element of the tuple is a list of information describing the variant gene pair
all_hits_file = qtl_results_dir + short_parameter_string + '_permutation_scheme_none_permute_False_merged_dynamic_qtl_results.txt'
egenes_file = qtl_results_dir + parameter_string + '_efdr_.05_significant_egenes.txt'
hits = extract_hits_vector(egenes_file, all_hits_file)
#hits = extract_hits_vector_v2(all_hits_file)

gene_mapping = make_gene_mapping(target_region_input_file)


# Determine the number of tests
num_tests = extract_number_of_tests(joint_test_input_file)


# First extract ordered list of:
## 1. Sample Ids
## 2. Filehandles opening input cht files
## 3. Environmental variables
## 4. library size correction factors
sample_ids, filehandles, environmental_vars, cell_lines = parse_joint_test_input_file(joint_test_input_file)
library_size_correction_factors = parse_correction_factor_file(correction_factor_file)

# Create dictionary the maps cell_line ids to an integer array. Where the array contains indices of each of the samples in that cell line
cell_line_indices = get_cell_line_indices(joint_test_input_file)

output_file = qtl_visualization_dir + parameter_string + '_dynamic_qtl_hits_summary.txt'

t = open(output_file, 'w')

max_sites = 4
counter = 0
# Loop through each of the tests
for test_number in range(num_tests):
    # Extract one line from each file handle and put into ordered list
    test_infos = extract_test_info_from_filehandles(filehandles)
    rs_id = test_infos[0][2]
    gene_identifier = test_infos[0][0] + '_' + test_infos[0][7] + '_' + test_infos[0][8]

    ensamble_id = gene_mapping[gene_identifier]
    # Only perform analysis on hits
    if rs_id + '_' + ensamble_id not in hits:
        continue
    hit_dicti = hits[rs_id + '_' + ensamble_id]
    # Re-organize test_infos data so that it is in a better format
    # Specifically, we will extract a vector of length == Number_of_samples of:
    ### 1. gene_counts
    ### 2. allele_1 read counts (ys)
    ### 3. Sum of allele_1 and allele_2 read counts (ns)
    ### 4. Genotype on allele_1 (h_1)
    ### 5. Genotype on allele_2 (h_2)


    gene_counts, ys, ns, h_1, h_2 = reorganize_data(test_infos)

    # If the dimension of ys and ns (Ie. the number of heterozygous exonic sites) is too large, the algorithm will become prohibitively slow
    # So we subset to the top $max_sites with the largest number of minor allele read counts
    ys, ns = subset_allelic_count_matrices(ys, ns, cell_line_indices, min_num_biallelic_lines, min_num_biallelic_samples, min_num_het_test_variant_biallelic_samples, h_1, h_2)
    

    genotype = (h_1 + h_2).astype(str)
    gene_counts = gene_counts/library_size_correction_factors
    num_sites = ys.shape[1]
    if num_sites > max_sites:
        num_sites = max_sites
    t.write(ensamble_id + '\t' + rs_id + '\t' + str(hit_dicti['pvalue']) + '\t' + str(hit_dicti['beta']) + '\t' + ';'.join(np.asarray(environmental_vars).astype(str)) + '\t')
    t.write(';'.join(np.asarray(gene_counts).astype(str)) + '\t' + ';'.join(genotype) + '\t' + hit_dicti['conc'][0])
    if num_sites > 10:
        conc = hit_dicti['conc'].astype(float)
        ys = ys[:, conc.argsort()[::-1]]
        ns = ns[:, conc.argsort()[::-1]]

    for site_num in range(max_sites):
        if site_num + 1 > num_sites:
            t.write('\tNA\tNA\tNA')
        else:
            a1_counts, a2_counts, environmental_vars_subset = get_allelic_counts(ys, ns, h_1, h_2, site_num, environmental_vars)
            if len(a1_counts) == 0:
                t.write('\tNA\tNA\tNA')
            else:
                t.write('\t' + ';'.join(a1_counts.astype(str)) + '\t' + ';'.join(a2_counts.astype(str)) + '\t' + ';'.join(environmental_vars_subset.astype(str)))
    t.write('\n')
    t.flush()

