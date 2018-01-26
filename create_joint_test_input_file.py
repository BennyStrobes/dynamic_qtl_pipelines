import numpy as np 
import os
import sys
import pdb

# Function to extract environmental_variable when environmental_variable_form=='time_steps'
def extract_environmental_variable_time_step_form(sample_name):
    return sample_name.split('_')[1]


input_directory = sys.argv[1]  # Directory containing all CHT input files
output_file = sys.argv[2]  # Output file
environmental_variable_form = sys.argv[3]  # String option describing how to parameterize the environmetal variable

t = open(output_file, 'w')  # Open output file handle
t.write('sample_id\tenvironmental_variable\tcht_input_file\n')

for file_name in sorted(os.listdir(input_directory)):
    if file_name.startswith('haplotype_read_') == False:
        continue
    if file_name.endswith('.txt.gz') == False:
        continue
    # Extract sample id from filename
    sample_name = file_name.split('.')[2]
    # Extract environmental variable (depends on environmental_variable_form)
    if environmental_variable_form == 'time_steps':
        environmental_variable = extract_environmental_variable_time_step_form(sample_name)
    
    # Print information to output file
    t.write(sample_name + '\t' + environmental_variable + '\t' + input_directory + file_name + '\n')
t.close()