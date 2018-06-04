import numpy as np
import pystan
import pdb
from scipy import stats
import pickle
import pandas as pd



# Load in data into pystan format for null model (not including interaction term)
def load_in_null_data(gene_counts, dosage, environmental_vars, library_size_correction_factors, model_version, covs, covariate_method):
    N = len(gene_counts)
    T = int(np.max(environmental_vars) + 1)
    intercept = np.ones((N, 1))
    env_mat = np.transpose(np.asmatrix(environmental_vars))
    dosage_mat = np.transpose(np.asmatrix(dosage))
    if covariate_method == 'none':
        x_1 = np.hstack((intercept, env_mat, np.square(env_mat), dosage_mat))
    elif covariate_method == 'cell_line_pc1Xtime':
        covariate_interaction_arr = covs['cell_line_pc1']*environmental_vars
        covariate_interaction_mat = np.transpose(np.asmatrix(covariate_interaction_arr))
        covariate_interaction_arr_squared = covs['cell_line_pc1']*np.square(environmental_vars)
        covariate_interaction_mat_squared = np.transpose(np.asmatrix(covariate_interaction_arr_squared))
        x_1 = np.hstack((intercept, env_mat, np.square(env_mat), dosage_mat, covariate_interaction_mat, covariate_interaction_mat_squared))
    data = dict(N=N, P=x_1.shape[1], library_size=library_size_correction_factors, x_1=x_1, gene_counts=gene_counts, concShape=1.001, concRate=0.001)
    return data


def permute_environmental_vars_within_cell_line_ordered(environmental_vars, cell_line_indices):
    # initialize output vector
    environmental_vars_perm = np.zeros(len(environmental_vars))

    # Get largest and smallestenvironmental variable
    maxy = max(environmental_vars)
    miny = min(environmental_vars)

    # Create mapping from time step to permuted time step
    varz = np.arange(miny,maxy+1)
    perm_varz = np.random.permutation(varz)
    print(perm_varz)
    mapping = {}
    for i, var in enumerate(varz):
        mapping[var] = perm_varz[i]

    for i, environmental_var in enumerate(environmental_vars):
        environmental_vars_perm[i] = mapping[environmental_var]

    return environmental_vars_perm

def permute_environmental_vars_within_cell_line(environmental_vars, cell_line_indices):
    environmental_vars_perm = np.zeros(len(environmental_vars))
    # Seperately permute all samples within each cell line
    for cell_line in cell_line_indices:
        # Get indices of samples that belong to this cell line
        index = cell_line_indices[cell_line]
        environmental_vars_perm[index] = np.random.permutation(environmental_vars[index])
    return environmental_vars_perm

def permute_environmental_vars_all(environmental_vars):
    environmental_vars_perm = np.zeros(len(environmental_vars))
    # Seperately permute all samples within each cell line
    environmental_vars_perm = np.random.permutation(environmental_vars)
    return environmental_vars_perm


def permute_environmental_vars_within_hets(environmental_vars, h_1, h_2):
    # Seperately permute homozygous and heterozygous individuals
    hets = np.where(h_1 != h_2)[0]  # heterozygous individuals
    homos = np.where(h_1 == h_2)[0]  # homozygous individuals
    # Initialize permuted environmental vars vectors
    environmental_vars_perm = np.zeros(len(environmental_vars))
    # Permute heterozygous individuals
    environmental_vars_perm[hets] = np.random.permutation(environmental_vars[hets])
    # Permute homozygous individuals
    environmental_vars_perm[homos] = np.random.permutation(environmental_vars[homos])
    return environmental_vars_perm

# Optimize the null model
def null_model_optimization_shell(null_data, sm, optimization_method):
    # Run test by placing in try catch loop
    # If doesn't converge, then try it again with a different seed
    working = True
    iteration = 1
    while working:
        try:
            # Use same seed for null and alternate models
            seed = np.random.randint(10000000) + 1
            # Run dynamic qtls                                                                                   
            op_null = sm.optimizing(data=null_data, as_vector=False, seed=seed, algorithm=optimization_method, tol_obj=1e-17, tol_rel_obj=1e0, tol_grad=1e-13, tol_rel_grad=1e2, tol_param=1e-14)
            working = False
            # Make sure log likelihood is not nan
            if np.isnan(op_null['value']):
                working = True
                # Force to throw an exception
                exception_thrower = 5.0/0.0
        except:
            print('Starting over for ' + str(iteration) + ' time')
            iteration = iteration + 1
            # If you've tried 40 different start seeds and none of them converged, then do LBFGS
            if iteration > 10:
                print('LBFGS ran')
                seed = np.random.randint(10000000) + 1
                op_null = sm.optimizing(data=null_data, as_vector=False, seed=seed, algorithm='LBFGS', tol_obj=1e-17, tol_rel_obj=1e0, tol_grad=1e-13, tol_rel_grad=1e2, tol_param=1e-14)
    return op_null

# Sample total reads from negative binomial distribution
def sample_gene_counts_from_negative_binomial_distribution(y_1, nb_conc, library_size_correction_factors):
    # Initialize output vector
    sample_gene_counts = []
    
    # Number of samples we need to sample
    num_samples = len(library_size_correction_factors)
    # Loop through each sample
    for n in np.arange(num_samples):
        # mean of negative binomial
        mu_n = (y_1[n, 0])*library_size_correction_factors[n]
        # Compute parameters of alternative parameterization of negative binomial
        beta_n = nb_conc/mu_n
        alpha_n = mu_n*beta_n
        p_n = beta_n/(beta_n + 1.0)

        # Take random sample from parameterized negative binomial
        try:
            sample = np.random.negative_binomial(alpha_n, p_n)
        except:
            print('miss ' + str(n))
            print(alpha_n)
            print(p_n)
            print(beta_n)
            print(mu_n)
            sample = 0
        sample_gene_counts.append(sample)

    return sample_gene_counts



# Sample read counts from fitted null model
def draw_samples_from_fitted_null_te_model(nb_conc, beta, null_data):
    beta_mat = np.transpose(np.asmatrix(beta))
    # Predict relative log number of reads 
    xb_1 = np.dot(null_data['x_1'], beta_mat)
    # Predict relative number of reads 
    y_1 = np.exp(xb_1)

    # Sample total reads from negative binomial distribution
    sample_gene_counts = sample_gene_counts_from_negative_binomial_distribution(y_1, nb_conc, null_data['library_size'])

    return np.asarray(sample_gene_counts).astype(int)

def permute_dosage_vars_between_cell_line_blocks(dosage, cell_line_indices):
    num_cell_lines = len(cell_line_indices)
    dosage_perm = np.zeros(len(dosage))

    permuted_line_ordering = np.random.permutation(np.arange(num_cell_lines))

    line_names = []
    for cell_line in cell_line_indices.keys():
        line_names.append(cell_line)

    for cell_line_number, cell_line in enumerate(line_names):
        replacement_line = line_names[permuted_line_ordering[cell_line_number]]
        dosage_perm[cell_line_indices[cell_line]] = np.zeros(len(cell_line_indices[cell_line])) + dosage[cell_line_indices[replacement_line]][0]
    return dosage_perm




# Load in data into pystan format for full model (including interaction term)
def load_in_full_data(gene_counts, dosage, environmental_vars, library_size_correction_factors, permute, permutation_scheme, null_data, sm, cell_line_indices, optimization_method, model_version, covs, covariate_method):
    N = len(gene_counts)
    T = int(np.max(environmental_vars) + 1)
    intercept = np.ones((N, 1))
    env_mat = np.transpose(np.asmatrix(environmental_vars))
    dosage_mat = np.transpose(np.asmatrix(dosage))
    # Permute the data
    if permute == 'True':
        # Run permutation independently in each cell line
        if permutation_scheme == 'shuffle_lines':
            environmental_vars_perm = permute_environmental_vars_within_cell_line(environmental_vars, cell_line_indices)
            interaction_mat = np.transpose(np.asmatrix(dosage*environmental_vars_perm))
            interaction_mat_squared = np.transpose(np.asmatrix(dosage*np.square(environmental_vars_perm)))
        # Run permutation for all samples
        elif permutation_scheme == 'shuffle_all':
            environmental_vars_perm = permute_environmental_vars_all(environmental_vars)
            interaction_mat = np.transpose(np.asmatrix(dosage*environmental_vars_perm))
            interaction_mat_squared = np.transpose(np.asmatrix(dosage*np.square(environmental_vars_perm)))
        elif permutation_scheme == 'shuffle_all_include_covs':
            environmental_vars_perm = permute_environmental_vars_all(environmental_vars)
            environmental_vars = environmental_vars_perm
            interaction_mat = np.transpose(np.asmatrix(dosage*environmental_vars_perm))
            interaction_mat_squared = np.transpose(np.asmatrix(dosage*np.square(environmental_vars_perm)))
        elif permutation_scheme == 'shuffle_all_time':
            environmental_vars_perm = permute_environmental_vars_all(environmental_vars)
            environmental_vars = environmental_vars_perm
            env_mat = np.transpose(np.asmatrix(environmental_vars))
            interaction_mat = np.transpose(np.asmatrix(dosage*environmental_vars_perm))
            interaction_mat_squared = np.transpose(np.asmatrix(dosage*np.square(environmental_vars_perm)))
        #  Fit null model. Use parameters from fitted null models to draw samples (gene counts).
        #  Then run LRT on sampled data
        elif permutation_scheme == 'sample_null':
            interaction_mat = np.transpose(np.asmatrix(dosage*environmental_vars))
            interaction_mat_squared = np.transpose(np.asmatrix(dosage*np.square(environmental_vars)))
            environmental_vars_perm = environmental_vars
            # Optimize the null model
            op_null = null_model_optimization_shell(null_data, sm, optimization_method)
            if model_version == 'te_log_linear_quadratic_basis':
                gene_counts = draw_samples_from_fitted_null_te_model(np.atleast_1d(op_null['par']['nb_conc'])[0], op_null['par']['beta'], null_data)
                null_data['gene_counts'] = gene_counts
        else:
            print('error: permutation scheme ' + permutation_scheme + ' currently not implemented')
    elif permute == 'False':  # Do not permute the data
        environmental_vars_perm = environmental_vars
        interaction_mat = np.transpose(np.asmatrix(dosage*environmental_vars))
        interaction_mat_squared = np.transpose(np.asmatrix(dosage*np.square(environmental_vars)))
    if covariate_method == 'none':
        x_1 = np.hstack((intercept, env_mat, np.square(env_mat), dosage_mat, interaction_mat, interaction_mat_squared))
        null_data['x_1'] = x_1[:,:-2]
    elif covariate_method == 'cell_line_pc1Xtime':
        covariate_interaction_arr = covs['cell_line_pc1']*environmental_vars
        covariate_interaction_mat = np.transpose(np.asmatrix(covariate_interaction_arr))
        covariate_interaction_arr_squared = covs['cell_line_pc1']*np.square(environmental_vars)
        covariate_interaction_mat_squared = np.transpose(np.asmatrix(covariate_interaction_arr_squared))      
        x_1 = np.hstack((intercept, env_mat, np.square(env_mat), dosage_mat, covariate_interaction_mat, covariate_interaction_mat_squared, interaction_mat, interaction_mat_squared))
        null_data['x_1'] = x_1[:,:-2]


    full_data = dict(N=N, P=x_1.shape[1], library_size=library_size_correction_factors, x_1=x_1, gene_counts=gene_counts, concShape=1.001, concRate=0.001)
    return full_data, null_data, environmental_vars_perm


# Run model with only expr ~ intercept + time + time^2 + genotype.
# Use results for initialization of real model
def run_simple_initialization_for_quadratic(null_data_temp,full_data, seed,algorithm, sm):
    temp_data = {}
    if null_data_temp['P'] == 6:
        temp_data['x_1'] = null_data_temp['x_1'][:,:-2]
        temp_data['P'] = null_data_temp['P'] - 2
    else:
        temp_data['x_1'] = null_data_temp['x_1']
        temp_data['P'] = null_data_temp['P']
    temp_data['library_size'] = null_data_temp['library_size']
    temp_data['gene_counts'] = null_data_temp['gene_counts']
    temp_data['concShape'] = null_data_temp['concShape']
    temp_data['concRate'] = null_data_temp['concRate']
    temp_data['N'] = null_data_temp['N']
    op_init = sm.optimizing(data=temp_data, as_vector=False, seed=seed, algorithm=algorithm, tol_obj=1e-17, tol_rel_obj=1e0, tol_grad=1e-13, tol_rel_grad=1e2, tol_param=1e-14)
    beta_init = np.zeros(full_data['P'])
    beta_init[:temp_data['P']] = op_init['par']['beta']
    init_param = dict(nb_conc=op_init['par']['nb_conc'], beta=beta_init)
    return init_param




def run_dynamic_qtl(sm, null_data, full_data, dof, algorithm, iteration, model_version):
    # Use same seed for null and alternate models
    seed = np.random.randint(10000000) + 1

    # Run model with only expr ~ intercept + time + time^2 + genotype.
    # Use results for initialization of real model
    simple_init = run_simple_initialization_for_quadratic(null_data, full_data, seed, algorithm, sm)

    # Run pystan gradient based optimization on full model
    op_full = sm.optimizing(data=full_data, as_vector=False, init=simple_init, algorithm=algorithm, tol_obj=1e-17, tol_rel_obj=1e0, tol_grad=1e-13, tol_rel_grad=1e2, tol_param=1e-14)

    # Initialize null model with parameters defining the full model
    # initialization for joint model
    init_null = dict(nb_conc=op_full['par']['nb_conc'], beta=op_full['par']['beta'][:(len(op_full['par']['beta'])-dof)])


    # Run pystan gradient based optimization on null model
    if iteration == 1:
        #sm.optimizing(data=null_data, as_vector=False, seed=seed, algorithm=algorithm, tol_obj=1e-15, tol_rel_obj=1e1, tol_grad=1e-11, tol_rel_grad=1e4, tol_param=1e-12)
        op_null = sm.optimizing(data=null_data, as_vector=False, init=init_null, seed=seed, algorithm=algorithm, tol_obj=1e-17, tol_rel_obj=1e0, tol_grad=1e-13, tol_rel_grad=1e2, tol_param=1e-14)
    else:  # Don't initialize with full model on later iterations
        op_null = sm.optimizing(data=null_data, as_vector=False, seed=seed, algorithm=algorithm, tol_obj=1e-17, tol_rel_obj=1e0, tol_grad=1e-13, tol_rel_grad=1e2, tol_param=1e-14)

    # Compute chi-squared test statistic
    loglr = op_full['value'] - op_null['value']

    # Consider possibility that null did not fully converge
    if (loglr > 3):
        refit_init = simple_init
        refit_init['beta'] = refit_init['beta'][:-2]
        # Refit the null with new random init
        refit_null = sm.optimizing(data=null_data, as_vector=False, init=refit_init, algorithm=algorithm, tol_obj=1e-17, tol_rel_obj=1e0, tol_grad=1e-13, tol_rel_grad=1e2, tol_param=1e-14)
        # Use the null model that has the higher likelihood
        if refit_null['value'] > op_null['value']:
            op_null = refit_null
            loglr = op_full['value'] - op_null['value']
    return op_null, op_full, loglr










def dynamic_qtl(gene_counts, dosage, environmental_vars, library_size_correction_factors, model_version, permute, optimization_method, cell_line_indices, permutation_scheme, covs, covariate_method):
    # Load in correct model
    if model_version == 'te_log_linear_quadratic_basis':
        #sm = pystan.StanModel(file='te_log_linear.stan')
        #f = open('te_log_linear.pkl','wb')
        #pickle.dump(sm, f)
        #f.close()
        sm = pickle.load(open('te_log_linear.pkl', 'rb'))

    # Load in data into pystan format
    null_data = load_in_null_data(gene_counts, dosage, environmental_vars, library_size_correction_factors, model_version, covs, covariate_method)
    full_data, null_data, perm_env_vars = load_in_full_data(gene_counts, dosage, environmental_vars, library_size_correction_factors, permute, permutation_scheme, null_data, sm, cell_line_indices, optimization_method, model_version, covs, covariate_method)
    # Calculate the degrees of freedom of LRT
    dof = full_data['P'] - null_data['P']
    # Run test by placing in try catch loop
    # If doesn't converge, then try it again with a different seed
    working = True
    iteration = 1
    while working:
        try:
            # Run dynamic qtls
            op_null, op_full, loglr = run_dynamic_qtl(sm, null_data, full_data, dof, optimization_method, iteration, model_version)
            working = False
            # Make sure log likelihood is not nan
            if np.isnan(op_null['value']) or np.isnan(op_full['value']) or np.isnan(loglr):
                working = True
                # Force to throw an exception
                exception_thrower = 5.0/0.0
        except:
            print('Starting over for ' + str(iteration) + ' time')
            iteration = iteration + 1
            # If you've tried 10 different start seeds and none of them converged, then do LBFGS
            if iteration > 10:
                print('FAILURE')
                print(null_data)
                print(full_data)
                op_null, op_full, loglr = run_dynamic_qtl(sm, null_data, full_data, dof, "LBFGS", iteration, model_version)

    # Compute pvalue from chi-squared test statistic
    pvalue = 1.0 - stats.chi2.cdf(2.0*loglr, dof)
    
    return dict(pvalue=pvalue, loglr=loglr, fit_full=op_full, fit_null=op_null, loglike_null=op_null['value'], loglike_full=op_full['value'], dosage=np.squeeze(np.asarray(null_data['x_1'][:,-1])), test_time_steps=perm_env_vars)
