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
    gene_names = []
    for gene_index, gene_id in enumerate(sorted(example_converter.keys())):
        # Loop through samples
        gene_names.append(gene_id)
        for sample_index in range(N):
            gene_counts[sample_index, gene_index] = conversions[sample_index][gene_id]
    return gene_counts.astype(int), gene_names

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
    time_steps = []
    count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # skip header
            head_count = head_count + 1
            continue
        count = count + 1
        sample_id = data[0]
        time_steps.append(int(sample_id.split('_')[1]))
        print(sample_id)
        cht_input_file = data[2]
        # Create dictionary to convert from gene id to read count for this sample
        converter, library_size = get_gene_id_to_read_count(cht_input_file)
        conversions.append(converter)
        library_sizes.append(library_size)
        samples.append(sample_id)
    return conversions, library_sizes, samples, np.asarray(time_steps) + 1

# Command Line Args
joint_test_input_file = sys.argv[1]  # Input file containing one line per sample with location of cht input file
output_file = sys.argv[2]  # Output file containing one line per sample (in same order as joint_test_input_file) with correction_factors
correction_factor_file = sys.argv[3]
variance_output_file = sys.argv[4]

# Load in pystan object
sm = pystan.StanModel(file='time_step_nb_overdispersion.stan')

library_size_correction_factors = np.loadtxt(correction_factor_file,dtype=str)[1:,1].astype(float)

# Create list of dictionaries. List is of length number of samples
# Each dictionary converts from gene_id to read count for this sample
conversions, library_sizes, samples, time_steps = extract_read_counts(joint_test_input_file)

# Convert list of dictionaries to a gene count matrix of dimeinsions N X P where N is number of samples and P is genes
gene_counts, gene_names = convert_data_to_matrix(conversions)

t = open('/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data/input_data/gene_counts.txt','w')
t.write('gene_name\t' + '\t'.join(np.asarray(samples)) + '\n')
for i, gene_name in enumerate(gene_names):
    t.write(gene_name + '\t')
    gene_vec = gene_counts[:,i].astype(str)
    t.write('\t'.join(gene_vec) + '\n')
t.close()

# Extract number of time steps
T = int(np.max(time_steps))

# Get data into corect format for pystan
data = dict(N=gene_counts.shape[0], library_size=library_size_correction_factors, P=gene_counts.shape[1], T=T, gene_counts=gene_counts, time_step=time_steps, concShape=1.01, concRate=0.01)



# Run pystan optimization
op = sm.optimizing(data=data)

# Print overdispersion parameters
t_handle = open(output_file, 'w')
for i, gene_name in enumerate(gene_names):
    t_handle.write(gene_name)
    for t in range(T):
        t_handle.write('\t' + str(op['conc'][i,t]))
    t_handle.write('\n')
t_handle.close()

variance_matrix = np.zeros((len(gene_names), T))
for i in range(len(gene_names)):
    for t in range(T):
        time_step_indices = np.where(np.asarray(time_steps)==(t+1))[0]
        variance_matrix[i,t] = np.var(gene_counts[time_step_indices,i])

t_handle = open(variance_output_file, 'w')
for i, gene_name in enumerate(gene_names):
    t_handle.write(gene_name)
    for t in range(T):
        t_handle.write('\t' + str(variance_matrix[i,t]))
    t_handle.write('\n')
t_handle.close()
    
