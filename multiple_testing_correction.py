import numpy as np 
import os
import sys
import pdb


# Load in list of tuples (pvalue, qtl_line)
# Also extract header
def get_pvalue_line_pairs(input_file):
    pairs = []  # Initialize output array
    head_count = 0 # Skip header
    # Stream input file
    f = open(input_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            header = line
            continue
        pvalue = float(data[-3])
        pair = (pvalue, line)
        pairs.append(pair)
    return pairs, header


# Perform benjamini hochberg for multiple testing correction
def bh_correction(input_file, output_file):
    # Output file handle
    t = open(output_file, 'w')

    # Load in list of tuples (pvalue, qtl_line)
    # Also extract header
    pvalue_line_pairs, header = get_pvalue_line_pairs(input_file)

    sorted_pvalue_line_pairs = sorted(pvalue_line_pairs, key=lambda x: x[0])

    t.write(header + '\n')

    alpha = .1  # FDR threshold
    m = len(sorted_pvalue_line_pairs)  # total number of tests
    valid = True

    # Loop through ordered list of qtl pairs
    for i, tupler in enumerate(sorted_pvalue_line_pairs):
        k = i + 1
        pvalue = tupler[0]
        line = tupler[1]
        if pvalue <= (float(k)/m)*alpha and valid == True:
            t.write(line + '\n')
        else:  # Once a test does not pass, all subsequent tests do not pass
            valid = False
    t.close()




# Command line args!

qtl_results_dir = sys.argv[1]
parameter_string = sys.argv[2]


# The raw results of dynamic qtl calling on ipsc data
input_file = qtl_results_dir + parameter_string + '_dynamic_qtl_results.txt'

# Subset input file to contain only genes that pass multiple testing correction
output_file = qtl_results_dir + parameter_string + '_dynamic_qtl_multiple_testing_results.txt'

# Perform benjamini hochberg for multiple testing correction
bh_correction(input_file, output_file)