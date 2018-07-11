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

def extract_cell_lines_largest_pseudotime_variable(sample_name, pseudotime_predictions_5_file):
    cell_line = sample_name.split('_')[0]
    f = open(pseudotime_predictions_5_file)
    maxy = -1
    for line in f:
        line = line.rstrip()
        data = line.split(',')
        line_cell_line = data[0].split('_')[0]
        if line_cell_line != cell_line:
            continue
        pseudo = int(data[1])
        if pseudo > maxy:
            maxy = pseudo

    f.close()
    if maxy == -1:
        print('fatal assumption errorr')
        pdb.set_trace()
    return str(maxy)

def extract_cell_line_pc(cell_line_pc_file, sample_name, pc_num):
    f = open(cell_line_pc_file)
    head_count = 0
    cell_liner = sample_name.split('_')[0]
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        if data[0] == cell_liner:
            return data[pc_num]

def get_hmm_cell_line_grouping(hmm_cell_line_groupings_file, sample_name):
    head_count = 0  # Skip header
    cell_line = sample_name.split('_')[0]
    f = open(hmm_cell_line_groupings_file)
    for line in f:
        line = line.rstrip()
        data = line.split(',')
        if head_count == 0:
            head_count = head_count + 1
            continue
        cell_line_cur = data[0]
        state = data[1]
        if cell_line_cur == cell_line:
            return state
    f.close()
    print('fundamentla assumption error')
    pdb.set_trace()

input_directory = sys.argv[1]  # Directory containing all CHT input files
output_file = sys.argv[2]  # Output file
environmental_variable_form = sys.argv[3]  # String option describing how to parameterize the environmetal variable
total_expression_file = sys.argv[4]
cell_line_pc_file = sys.argv[5]
hmm_cell_line_groupings_8_state_file = sys.argv[6]
hmm_cell_line_groupings_12_state_file = sys.argv[7]
hmm_cell_line_groupings_16_state_file = sys.argv[8]

t = open(output_file, 'w')  # Open output file handle
t.write('sample_id\tenvironmental_variable\tcht_input_file\ttroponin_t15\tcell_line_pc1\tcell_line_pc2\tcell_line_pc3\tcell_line_pc4\tcell_line_pc5\thmm_cell_line_grouping_8_state\thmm_cell_line_grouping_12_state\thmm_cell_line_grouping_16_state\n')

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

    cell_line_pc1 = extract_cell_line_pc(cell_line_pc_file, sample_name, 1)
    cell_line_pc2 = extract_cell_line_pc(cell_line_pc_file, sample_name, 2)
    cell_line_pc3 = extract_cell_line_pc(cell_line_pc_file, sample_name, 3)
    cell_line_pc4 = extract_cell_line_pc(cell_line_pc_file, sample_name, 4)
    cell_line_pc5 = extract_cell_line_pc(cell_line_pc_file, sample_name, 5)


    hmm_cell_line_grouping_8_state = get_hmm_cell_line_grouping(hmm_cell_line_groupings_8_state_file, sample_name)
    hmm_cell_line_grouping_12_state = get_hmm_cell_line_grouping(hmm_cell_line_groupings_12_state_file, sample_name)
    hmm_cell_line_grouping_16_state = get_hmm_cell_line_grouping(hmm_cell_line_groupings_16_state_file, sample_name)

    # Print information to output file
    t.write(sample_name + '\t' + environmental_variable + '\t' + input_directory + file_name + '\t' + troponin_expression + '\t' + cell_line_pc1 + '\t' + cell_line_pc2 + '\t' + cell_line_pc3 + '\t' + cell_line_pc4 + '\t' + cell_line_pc5 + '\t' + hmm_cell_line_grouping_8_state + '\t' + hmm_cell_line_grouping_12_state + '\t' + hmm_cell_line_grouping_16_state + '\n')
t.close()