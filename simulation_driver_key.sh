#######################################################################################################
# This is the driver script that calls several simulation/debugging tests on various EAGLE type models
#######################################################################################################


##########################
# Output directories (aasume all of these exist prior to starting analysis)
##########################

## Root directory of all simulation based analysis
simulation_root="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/simulation/"

# Directory that contains data for part 1 analysis
initial_exploration_dir=$simulation_root"initial_exploration/"




#############################################################
# Part 1: Initial Exploration
# Explore betabinomial_glmm and betabinomial_glm 
#############################################################

Rscript test_bb_glmm.R