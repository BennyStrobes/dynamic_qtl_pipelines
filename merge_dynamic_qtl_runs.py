import numpy as np 
import os
import sys
import pdb
import gzip

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



# Now print all of the parrallel output files to one merged output file
# Also, each line will contain the gene id
def merge_parrallel_files(real_parameter_string, output_file, gene_mapping, total_jobs):
    # open up output file handle
    t = open(output_file, 'w')
    # Write header
    t.write('chrom_num\tvariant_position\trs_id\tref_allele\talt_allele\tensamble_id\tnull_log_likelihood\tfull_log_likelihood\tloglr\tpvalue\tbeta\tconc_factors\n')
    # Loop through each parrallel job
    for job_number in range(total_jobs):
        print(job_number)
        input_file = real_parameter_string + '_' + str(job_number) + '_dynamic_qtl_results.txt'
        f = open(input_file)
        for line in f:
            line = line.rstrip()
            data = line.split()
            # Convert from gene identifier to ensamble id
            gene_identifier = data[0] + '_' + data[5] + '_' + data[6]
            ensamble_id = gene_mapping[gene_identifier]
            # print line to output file
            if len(data) == 13: # have  allelic reads
                t.write('\t'.join(data[0:5]) + '\t' + ensamble_id + '\t' + '\t'.join(data[7:]) + '\n')
            elif len(data) == 12: # no allelic reads
                t.write('\t'.join(data[0:5]) + '\t' + ensamble_id + '\t' + '\t'.join(data[7:]) + '\tNA\n')
            else:  # Just to make sure there was nothing we didnt account for
                print('erroro in assumptions')


    t.close()


###########################
# Command Line Args
###########################

real_parameter_string = sys.argv[1]
output_file = sys.argv[2]
target_region_input_file = sys.argv[3]
total_jobs = int(sys.argv[4])


# Make mapping from gene identifier to ensamble id
gene_mapping = make_gene_mapping(target_region_input_file)

# Now print all of the parrallel output files to one merged output file
# Also, each line will contain the gene id
merge_parrallel_files(real_parameter_string, output_file, gene_mapping, total_jobs)