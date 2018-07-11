#!/bin/bash
#SBATCH --time=20:00:00 --partition=broadwl --mem=25GB

module unload python
module load Anaconda3

joint_test_input_file="$1"
correction_factor_file="$2"



date
python learn_library_size_correction_factor.py $joint_test_input_file $correction_factor_file
