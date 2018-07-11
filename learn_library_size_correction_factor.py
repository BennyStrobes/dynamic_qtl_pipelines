import numpy as np 
import os
import sys
import pdb
import gzip
import pystan



# Create dictionary to convert from gene id to read count for this sample
def get_gene_id_to_read_count(cht_input_file):
    f = gzip.open(cht_input_file)
    head_count = 0 # used to skip header
    dicti = {} # initialize dictionary
    library_size = -1
    # stream lines of cht test input file
    for line in f:
        line = line.decode('utf-8').rstrip()
        data = line.split()
        if head_count == 0: # Skip header
            head_count = head_count + 1
            continue
        # Use exon locations to define gene_id
        gene_id = data[0] + '_' + data[7] + '_' + data[8]
        count = int(data[15])
        line_lib_size = int(data[16])
        if gene_id in dicti: # seen this gene before
            if count != dicti[gene_id]:
                print('ASSUMPTION ERROR')
                pdb.set_trace()
        dicti[gene_id] = count
        # Get library size for this sample
        if library_size == -1:
            library_size = line_lib_size
        else:
            if library_size != line_lib_size:
                print('erororoororor!!')
                pdb.set_trace()
    return dicti, library_size

# Convert list of dictionaries to a gene count matrix of dimeinsions N X P where N is number of samples and P is genes
def convert_data_to_matrix(conversions):
    N = len(conversions) # Num samples
    example_converter = conversions[0]
    P = len(example_converter) # Num genes
    # Initialize gene count matrix
    gene_counts = np.zeros((N,P))

    # Loop through genes
    for gene_index, gene_id in enumerate(sorted(example_converter.keys())):
        # Loop through samples
        for sample_index in range(N):
            gene_counts[sample_index, gene_index] = conversions[sample_index][gene_id]
    return gene_counts.astype(int)

# Create list of dictionaries. List is of length number of samples
# Each dictionary converts from gene_id to read count for this sample
def extract_read_counts(joint_test_input_file):
    # Open filehandle for input file
    f = open(joint_test_input_file)
    # Loop through samples
    head_count = 0  # used to skip header
    conversions = []  # Keep track of dictionaries for each sample
    library_sizes = [] # Keep track of true library size for each sample
    samples = []  # Keep track of sample ids
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # skip header
            head_count = head_count + 1
            continue
        sample_id = data[0]
        print(sample_id)
        cht_input_file = data[2]
        # Create dictionary to convert from gene id to read count for this sample
        converter, library_size = get_gene_id_to_read_count(cht_input_file)
        conversions.append(converter)
        library_sizes.append(library_size)
        samples.append(sample_id)
    return conversions, library_sizes, samples

# Command Line Args
joint_test_input_file = sys.argv[1]  # Input file containing one line per sample with location of cht input file
correction_factor_file = sys.argv[2]  # Output file containing one line per sample (in same order as joint_test_input_file) with correction_factors


# Load in pystan object
sm = pystan.StanModel(file='library_size.stan')


# Create list of dictionaries. List is of length number of samples
# Each dictionary converts from gene_id to read count for this sample
conversions, library_sizes, samples = extract_read_counts(joint_test_input_file)

# Convert list of dictionaries to a gene count matrix of dimeinsions N X P where N is number of samples and P is genes
gene_counts = convert_data_to_matrix(conversions)

# Get data into corect format for pystan
data = dict(N=gene_counts.shape[0], P=gene_counts.shape[1], gene_counts=gene_counts, concShape=1.01, concRate=0.01)

# Initialize correction factors to observed library sizes
init_param = dict(library_size=library_sizes)

# Run pystan optimization
op = sm.optimizing(data=data, init=init_param)

# Print library sizes
t = open(correction_factor_file, 'w')
t.write('sample_id\tlibrary_size_correction_factor\ttrue_library_size\n')
for i, sample_id in enumerate(samples):
    t.write(sample_id + '\t' + str(op['library_size'][i]) + '\t' + str(library_sizes[i]) + '\n')
t.close()
