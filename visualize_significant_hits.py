import numpy as np 
import os
import pdb
import sys
import gzip
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import seaborn as sns 
sns.set(font_scale=2.0)
sns.set_style("white")


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

############################################
# Extract vector of length number of tests where each element is a tuple
# The first element in the tuple is a binary variable that is a 1 if that variant-gene pair is a highest
# The second element of the tuple is a list of information describing the variant gene pair
def extract_hits_vector(egenes_file, hits_file):
    significant_variant_gene_pairs = extract_significant_variant_gene_pairs(egenes_file)
    hits_vector = []  # Initialize output vector
    head_count = 0  # Skip header
    used_genes = {}
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
        pvalue = float(data[-2])
        if rs_id + '_' + ensamble_id in significant_variant_gene_pairs:  # Hit
            # hit information
            dicti = {}
            dicti['chrom_num'] = data[0]
            dicti['variant_pos'] = data[1]
            dicti['rs_id'] = data[2]
            dicti['ref_allele'] = data[3]
            dicti['alt_allele'] = data[4]
            dicti['ensamble_id'] = data[5]
            dicti['pvalue'] = pvalue
            dicti['beta'] = float(data[-1])
            hits_vector.append((1,dicti))
            used_genes[ensamble_id] = 1
        else:  # Not a hit
            hits_vector.append((0,'Null'))
    return hits_vector



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

# If the dimension of ys and ns (Ie. the number of heterozygous exonic sites) is too large, the algorithm will become prohibitively slow
# So we subset to the top $max_sites with the largest number of minor allele read counts
def subset_allelic_count_matrices(ys, ns, max_sites):
    # The number of exonic sites
    K = ys.shape[1]
    if K <= max_sites:
        return ys, ns
    # Calculate the the sum number of reads per exonic site
    total_counts = np.sum(ns, axis=0)
    # Extract indices of top $max_sites exonic sites
    indices_to_keep = total_counts.argsort()[-max_sites:][::-1]
    return ys[:, indices_to_keep], ns[:, indices_to_keep]

def filter_allelic_count_matrices(ys, ns, min_reads):
    K = ys.shape[1]
    # Calculate the the sum number of reads per exonic site
    total_counts = np.sum(ns, axis=0)

    indexes = total_counts.argsort()[::-1]

    ys = ys[:, indexes]
    ns = ns[:, indexes]

    total_counts = np.sum(ns, axis=0)
    
    read_indexes = np.where(total_counts > min_reads)[0]

    return ys[:, read_indexes], ns[:, read_indexes]


# Convert from dosage vector to vector of genotyeps
def dosage_to_genotype(dosage, ref_allele, alt_allele):
    converter = {}
    converter[0] = ref_allele + ref_allele
    converter[1] = ref_allele + alt_allele
    converter[2] = alt_allele + alt_allele
    genotype = []
    for dos in dosage:
        genotype.append(converter[dos])
    return np.asarray(genotype)

# Plot to visualize total expression changes over time as a function of genotype
def gene_total_plotter(gene_counts, dosage, environmental_vars, cell_lines, rs_id, ensamble_id, ref_allele, alt_allele, pvalue, beta, library_size_correction_factors):
    # Convert from dosage vector to vector of genotyeps
    genotype = dosage_to_genotype(dosage, ref_allele, alt_allele)

    gene_counts = (gene_counts/library_size_correction_factors)*np.mean(library_size_correction_factors)

    df = pd.DataFrame({rs_id: genotype, 'time_step': environmental_vars.astype(int), 'cell_lines': cell_lines, 'gene_counts': np.log(gene_counts)})
    ax = sns.boxplot(x="time_step", y="gene_counts", hue=rs_id, data=df, palette="Set3",width=.7)
    plt.xlabel('Time Step')
    plt.ylabel('log(counts)')
    plt.title(ensamble_id + ' / pvalue = ' + str(pvalue) + ' / beta = ' + str(beta))
    sns.despine(offset=1, trim=True)
    return ax

def allelic_imbalence_plotter(h_1, h_2, environmental_vars, ys, ns, cell_lines, rs_id, ensamble_id, pvalue, beta):
    ys, ns = filter_allelic_count_matrices(ys, ns, 0)

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
        for exonic_site_num in np.arange(num_exonic_sites):
            if ns[sample_num, exonic_site_num] <= 2:
                continue
            if h_1[sample_num] == 0:
                allelic_fraction = float(ys[sample_num, exonic_site_num])/ns[sample_num, exonic_site_num]
            elif h_2[sample_num] == 0:
                allelic_fraction = 1 - (float(ys[sample_num, exonic_site_num])/ns[sample_num, exonic_site_num])
            else:
                print('eroroororo')
            depths.append(ns[sample_num,exonic_site_num])
            allelic_fractions.append(allelic_fraction)
            #allelic_fractions.append(abs(allelic_fraction-.5))
            time_steps.append(float(environmental_vars[sample_num]))
            cell_lines_arr.append(cell_lines[sample_num])
            exonic_sites.append(exonic_site_num)
            identifiers.append('site_' + str(exonic_site_num + 1) + '_' + cell_lines[sample_num])

    # PLOT!
    df = pd.DataFrame({'time_step': time_steps,'read_depth':depths, 'cell_lines': cell_lines_arr, 'exonic_sites': exonic_sites, 'allelic_fraction': allelic_fractions, 'identifiers': identifiers})

    #ax = sns.pointplot(x="time_step", y="allelic_fraction", hue="identifiers", data=df)
    #ax = sns.regplot(x="time_step", y="allelic_fraction", data=df)
    ax = sns.lmplot(data=df,x="time_step", y="allelic_fraction",col="exonic_sites", hue="exonic_sites",col_wrap=3,ci=None)
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

# PLOT!
def visualize_hit(ys, ns, gene_counts, h_1, h_2, environmental_vars, cell_lines, library_size_correction_factors, test_dicti, output_file):
    # Plot to visualize total expression changes over time as a function of genotype
    fig = plt.figure(figsize=(20,10))

    gene_total_plot = gene_total_plotter(gene_counts, h_1+h_2, environmental_vars, cell_lines, test_dicti['rs_id'], test_dicti['ensamble_id'], test_dicti['ref_allele'], test_dicti['alt_allele'], test_dicti['pvalue'], test_dicti['beta'], library_size_correction_factors)
    fig.savefig(output_file)
    #fig = plt.figure(figsize=(20,10))
    #allelic_imbalence_plot = allelic_imbalence_plotter(h_1, h_2, environmental_vars, ys, ns, cell_lines, test_dicti['rs_id'], test_dicti['ensamble_id'], test_dicti['pvalue'], test_dicti['beta'])
    #allelic_imbalence_plot.savefig(output_file)
    print('hit')

###########################
# Command Line args
###########################
parameter_string = sys.argv[1]
qtl_results_dir = sys.argv[2]
joint_test_input_file = sys.argv[3]
correction_factor_file = sys.argv[4]
max_sites = int(sys.argv[5])
qtl_visualization_dir = sys.argv[6]


############################################
# Extract vector of length number of tests where each element is a tuple
# The first element in the tuple is a binary variable that is a 1 if that variant-gene pair is a highest
# The second element of the tuple is a list of information describing the variant gene pair
all_hits_file = qtl_results_dir + parameter_string + '_permute_False_merged_dynamic_qtl_results.txt'
egenes_file = qtl_results_dir + parameter_string + '_qval_.05_significant_egenes.txt'
hits_vector = extract_hits_vector(egenes_file, all_hits_file)




# Determine the number of tests
num_tests = extract_number_of_tests(joint_test_input_file)

if num_tests != len(hits_vector):
    print('eroorrooror')

# First extract ordered list of:
## 1. Sample Ids
## 2. Filehandles opening input cht files
## 3. Environmental variables
## 4. library size correction factors
sample_ids, filehandles, environmental_vars, cell_lines = parse_joint_test_input_file(joint_test_input_file)
library_size_correction_factors = parse_correction_factor_file(correction_factor_file)

counter = 0
# Loop through each of the tests
for test_number in range(num_tests):
    # Extract one line from each file handle and put into ordered list
    test_infos = extract_test_info_from_filehandles(filehandles)
    # Only perform analysis on hits
    if hits_vector[test_number][0] == 0:  # Skip non-hits
        continue
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
    ys, ns = subset_allelic_count_matrices(ys, ns, max_sites)
    
    # Name of output file to save plot to
    output_file = qtl_visualization_dir + parameter_string + '_' + hits_vector[test_number][1]['ensamble_id'] + '_' + hits_vector[test_number][1]['rs_id'] + '_dynamic_qtl_summary.png'
   
    # PLOT! 
    #if ys.shape[1] == 0:
     #   print('miss')
     #   continue
    visualize_hit(ys, ns, gene_counts, h_1, h_2, environmental_vars, cell_lines, library_size_correction_factors, hits_vector[test_number][1], output_file)

