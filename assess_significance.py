import numpy as np
import os
import sys
import pdb

# First, get index of sample that has the largest qvalue less than qval_thresh
def get_sample_index(qvalue_file, qval_thresh):
    f = open(qvalue_file)
    min_diff = 1000000
    best_index = -5
    for counter, line in enumerate(f):
        line = line.rstrip()
        data = line.split()
        if counter == 0:
            continue
        qval = float(data[1])
        diff = abs(qval - qval_thresh)
        if qval < qval_thresh and diff < min_diff:
            min_diff = diff
            best_index = counter
    return best_index

        


# Get pvalue threshold corresponding to qval_thresh in null
def get_pval_thresh(null_file, qvalue_file, qval_thresh):
    # First, get index of sample that has the largest qvalue less than qval_thresh
    index = get_sample_index(qvalue_file, qval_thresh)
    f = open(null_file)
    pval_thresh = 0
    for i, line in enumerate(f):
        line = line.rstrip()
        data = line.split()
        if i == index:
            pvaller = float(data[9])
    return pvaller

def filter_to_significant_results(real_file, significant_results_file, pval_thresh):
    f = open(real_file)
    t = open(significant_results_file, 'w')
    head_count = 0 
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            t.write(line + '\n')
            continue
        pval = float(data[9])
        if pval > pval_thresh:  # Not significant
            continue
        t.write(line + '\n')
    t.close()

def filter_to_egenes(significant_results_file, significant_gene_results_file):
    f = open(significant_results_file)
    t = open(significant_gene_results_file,'w')
    genes = {}
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            t.write(line + '\n')
            continue
        pval = float(data[9])
        gene = data[5]
        if gene not in genes:
            genes[gene] = (pval, line)
        else:
            old_pval = genes[gene][0]
            old_line = genes[gene][1]
            if pval < old_pval:
                genes[gene] = (pval, line)
            else:
                genes[gene] = (old_pval,old_line)
    f.close()
    for gene in genes.keys():
        liner = genes[gene][1]
        t.write(liner + '\n')
    t.close()

null_file = sys.argv[1]
real_file = sys.argv[2]
qvalue_file = sys.argv[3]
significant_results_file = sys.argv[4]  #  Output file
significant_gene_results_file = sys.argv[5]  # Ouptut file
qval_thresh = float(sys.argv[6])


# Get pvalue threshold corresponding to qval_thresh in null
pval_thresh = get_pval_thresh(null_file, qvalue_file, qval_thresh)

# Create list of variant-gene pairs that have pvalue less than pvalue threshold
filter_to_significant_results(real_file, significant_results_file, pval_thresh)

# Create list of significant genes and their highest associated variant
filter_to_egenes(significant_results_file, significant_gene_results_file)