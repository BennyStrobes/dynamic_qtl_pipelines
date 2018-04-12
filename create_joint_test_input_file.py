import numpy as np 
import os
import sys
import pdb

# Function to extract environmental_variable when environmental_variable_form=='time_steps'
def extract_environmental_variable_time_step_form(sample_name):
    return sample_name.split('_')[1]

# Function to extract environmental_variable when environmental_variable_form=='time_steps'
def extract_environmental_variable_time_step_form_max(sample_name, max_int):
    temp_time_step = int(sample_name.split('_')[1])
    if temp_time_step > max_int:
        time_step = 'NA'
    else:
        time_step = str(temp_time_step)
    return time_step


def extract_environmental_variable_pseudotime_form(sample_name, mapping_file):
    f = open(mapping_file)
    for line in f:
        line = line.rstrip()
        data = line.split(',')
        line_sample = data[0]
        pseudotime = data[1]
        if line_sample == sample_name:
            return pseudotime
    f.close()
    # ERROR, sample name not in nirmals file
    print('ASSUMPTION ERROROR in nirmals pseudotime file')
    pdb.set_trace()
    return

def extract_environmental_variable_uniform_form(sample_name, num_states):
    true_time_step = int(sample_name.split('_')[1])
    return str(int(np.floor(true_time_step/num_states)))

def get_troponin_expression(sample_name, total_expression_file):
    t15_sample = sample_name.split('_')[0] + '_15'
    f = open(total_expression_file)
    desired_gene = 'ENSG00000118194'
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            for i, ele in enumerate(data):
                if ele == t15_sample:
                    index = i
            continue
        if data[0] == desired_gene:
            return data[index]
    f.close()
    print('assumption errroro!!')
    pdb.set_trace()
    return

def extract_variable_pseudotime_median(sample_name, pseudotime_file, num_states):
    f = open(pseudotime_file)
    for line in f:
        line = line.rstrip()
        data = line.split(',')
        if sample_name == data[0]:
            cell_line = data[0].split('_')[0]
            time_step = int(data[0].split('_')[1])
            pseudotime = int(data[1])
    f.close()
    avail = []
    f = open(pseudotime_file)
    for line in f:
        line = line.rstrip()
        data = line.split(',')
        curr_cell_line = data[0].split('_')[0]
        curr_time_step = int(data[0].split('_')[1])
        curr_pseudotime = int(data[1])
        if curr_cell_line == cell_line and curr_pseudotime == pseudotime:
            tupler = (curr_time_step, curr_pseudotime, curr_cell_line)
            avail.append(tupler)
    avail_sorted = sorted(avail, key=lambda tup: tup[0])
    if len(avail_sorted) % 2 == 0:  # even
        lower_median_index = len(avail_sorted)/2 - 1
    else: # odd
        lower_median_index = len(avail_sorted)/2
    lower_median_tuple = avail_sorted[lower_median_index]
    if lower_median_tuple[0] == time_step:
        return str(lower_median_tuple[1])
    else:
        return 'NA'


input_directory = sys.argv[1]  # Directory containing all CHT input files
output_file = sys.argv[2]  # Output file
environmental_variable_form = sys.argv[3]  # String option describing how to parameterize the environmetal variable
pseudotime_predictions_3_file = sys.argv[4]
pseudotime_predictions_4_file = sys.argv[5]
pseudotime_predictions_5_file = sys.argv[6]
total_expression_file = sys.argv[7]

t = open(output_file, 'w')  # Open output file handle
t.write('sample_id\tenvironmental_variable\tcht_input_file\ttroponin_t15\n')

for file_name in sorted(os.listdir(input_directory)):
    if file_name.startswith('haplotype_read_counts_rand_hap') == False:
        continue
    if file_name.endswith('.txt.gz') == False:
        continue
    # Extract sample id from filename
    sample_name = file_name.split('.')[2]
    # Extract environmental variable (depends on environmental_variable_form)
    if environmental_variable_form == 'time_steps':
        environmental_variable = extract_environmental_variable_time_step_form(sample_name)
    elif environmental_variable_form == 'pseudotime_predictions_3':
        environmental_variable = extract_environmental_variable_pseudotime_form(sample_name, pseudotime_predictions_3_file)
    elif environmental_variable_form == 'pseudotime_predictions_4':
        environmental_variable = extract_environmental_variable_pseudotime_form(sample_name, pseudotime_predictions_4_file)
    elif environmental_variable_form == 'pseudotime_predictions_5':
        environmental_variable = extract_environmental_variable_pseudotime_form(sample_name, pseudotime_predictions_5_file)
    elif environmental_variable_form == 'uniform_4':
        environmental_variable = extract_environmental_variable_uniform_form(sample_name, 4)
    elif environmental_variable_form == 'median_pseudotime_4':
        environmental_variable = extract_variable_pseudotime_median(sample_name, pseudotime_predictions_4_file, 4)
    elif environmental_variable_form == 'median_pseudotime_5':
        environmental_variable = extract_variable_pseudotime_median(sample_name, pseudotime_predictions_5_file, 5)
    elif environmental_variable_form == 'median_pseudotime_3':
        environmental_variable = extract_variable_pseudotime_median(sample_name, pseudotime_predictions_3_file, 3)
    elif environmental_variable_form == 'time_steps_max_8':
        environmental_variable = extract_environmental_variable_time_step_form_max(sample_name, 8)
    elif environmental_variable_form == 'time_steps_max_9':
        environmental_variable = extract_environmental_variable_time_step_form_max(sample_name, 9)

    if environmental_variable == 'NA':
        continue
    troponin_expression = get_troponin_expression(sample_name, total_expression_file)
    # Print information to output file
    t.write(sample_name + '\t' + environmental_variable + '\t' + input_directory + file_name + '\t' + troponin_expression + '\n')
t.close()