import GPy
import numpy as np
import pystan
import pdb
from scipy import stats
import pickle
import pandas as pd



# Load in data into pystan format for null model (not including interaction term)
def load_in_null_data(ys, ns, h_1, h_2, environmental_vars, model_version, covs, covariate_method):
    N = ys.shape[0]
    K = ys.shape[1]
    T = int(np.max(environmental_vars) + 1)
    intercept = np.ones((N, 1))
    env_mat = np.transpose(np.asmatrix(environmental_vars))
    h_1_mat = np.transpose(np.asmatrix(h_1))
    h_2_mat = np.transpose(np.asmatrix(h_2))
    if covariate_method == 'none':
        x_1 = np.hstack((intercept, env_mat, h_1_mat))
        x_2 = np.hstack((intercept, env_mat, h_2_mat))
    elif covariate_method == 't15_troponin':
        cov_mat = np.transpose(np.asmatrix(covs*environmental_vars))
        x_1 = np.hstack((intercept, env_mat, h_1_mat, cov_mat))
        x_2 = np.hstack((intercept, env_mat, h_2_mat, cov_mat))
    data = dict(N=N, K=K, T=T, time_step=(environmental_vars.astype(int) + 1), P=x_1.shape[1], x_1=x_1, x_2=x_2, ys=ys, ns=ns, concShape=1.001, concRate=0.001)
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
def sample_gene_counts_from_negative_binomial_distribution(y_1, y_2, nb_conc, library_size_correction_factors):
    # Initialize output vector
    sample_gene_counts = []
    
    # Number of samples we need to sample
    num_samples = len(library_size_correction_factors)
    # Loop through each sample
    for n in np.arange(num_samples):
        # mean of negative binomial
        mu_n = (y_1[n, 0] + y_2[n, 0])*library_size_correction_factors[n]
        # Compute parameters of alternative parameterization of negative binomial
        beta_n = nb_conc/mu_n
        alpha_n = mu_n*beta_n
        p_n = beta_n/(beta_n + 1.0)

        # Take random sample from parameterized negative binomial
        try:
            sample = np.random.negative_binomial(alpha_n, p_n)
        except:
            print(alpha_n)
            print(p_n)
            print(beta_n)
            print(mu_n)
        sample_gene_counts.append(sample)

    return sample_gene_counts


# Sample allelic counts from beta-binomial distribution
def sample_allelic_counts_from_beta_binomial_distribution(allelic_p, conc, ys, ns):
    # Extract number of samples and number of sites
    num_samples = ys.shape[0]
    num_exonic_sites = ys.shape[1]
    if num_exonic_sites == 0:
        print('zero sites')
        pdb.set_trace()
    # Initialize sampled allelic count matrices 
    sample_ys = np.zeros((num_samples, num_exonic_sites))
    sample_ns = np.zeros((num_samples, num_exonic_sites))

    # Loop through samples
    for n in np.arange(num_samples):
        # Loop through exonic sites
        for site_num in np.arange(num_exonic_sites):
            # If sites have zero reads in true data, sites have zero reads in simulated data
            if ns[n, site_num] == 0.0:
                continue
            # Predicted allelic fraction for this sample
            p_n = allelic_p[n, 0]
            # Convert to alternative parameterization of beta binomial
            alpha_n = p_n*conc[site_num]
            beta_n = (1.0-p_n)*conc[site_num]

            if alpha_n <= 0:
                print('n is less than or equal to zero')
                print(alpha_n)
                print(p_n)
                print(conc[site_num])

            # Sample p from beta
            sample_p = np.random.beta(alpha_n, beta_n)
            if np.isnan(sample_p):
                print(alpha_n)
                print(beta_n)
                print(p_n)
                print(conc[site_num])
            # Sample allelic counts from binomial using sampled p
            sample_y = np.random.binomial(ns[n, site_num], sample_p)
            # Fill in allelic count matrices
            sample_ns[n, site_num] = ns[n, site_num]
            sample_ys[n, site_num] = sample_y

    return sample_ys, sample_ns 

def draw_samples_from_fitted_null_as_model(conc, beta, null_data):
    beta_mat = np.transpose(np.asmatrix(beta))
    # Predict relative log number of reads on allele 1
    xb_1 = np.dot(null_data['x_1'], beta_mat)
    # Predict relative log number of reads on allele 2
    xb_2 = np.dot(null_data['x_2'], beta_mat)
    # Predict relative number of reads on allele 1
    y_1 = np.exp(xb_1)
    # Predict relative number of reads on allele 2
    y_2 = np.exp(xb_2)
    # Predicted allelic fraction for each sample
    allelic_p = y_1/(y_1 + y_2)

    # ERROR CHECKING!!
    if np.sum(np.isnan(allelic_p)) > 0:
        print('nan p')
        print(allelic_p)
        print(xb_1)
        print(xb_2)
        print(beta)
        print(null_data['x_1'])
        print(null_data['x_2'])

    sample_ys, sample_ns = sample_allelic_counts_from_beta_binomial_distribution(allelic_p, conc, null_data['ys'], null_data['ns'])
    return sample_ys.astype(int), sample_ns.astype(int)




# Load in data into pystan format for full model (including interaction term)
def load_in_full_data(ys, ns, h_1, h_2, environmental_vars, permute, permutation_scheme, null_data, sm, cell_line_indices,optimization_method, model_version, covs, covariate_method):
    N = ys.shape[0]
    K = ys.shape[1]
    T = int(np.max(environmental_vars) + 1)
    intercept = np.ones((N, 1))
    env_mat = np.transpose(np.asmatrix(environmental_vars))
    h_1_mat = np.transpose(np.asmatrix(h_1))
    h_2_mat = np.transpose(np.asmatrix(h_2))
    # Permute the data
    if permute == 'True':
        # Run permutation independently in each cell line
        if permutation_scheme == 'shuffle_lines':
            environmental_vars_perm = permute_environmental_vars_within_cell_line(environmental_vars, cell_line_indices)
            h_1_interaction_mat = np.transpose(np.asmatrix(h_1*environmental_vars_perm))
            h_2_interaction_mat = np.transpose(np.asmatrix(h_2*environmental_vars_perm))
        elif permutation_scheme == 'shuffle_lines_ordered':
            environmental_vars_perm = permute_environmental_vars_within_cell_line_ordered(environmental_vars, cell_line_indices)
            h_1_interaction_mat = np.transpose(np.asmatrix(h_1*environmental_vars_perm))
            h_2_interaction_mat = np.transpose(np.asmatrix(h_2*environmental_vars_perm))
        # Run permutation independently in heterozygotes and homozygotes
        elif permutation_scheme == 'shuffle_hets':
            environmental_vars_perm = permute_environmental_vars_within_hets(environmental_vars, h_1, h_2)
            h_1_interaction_mat = np.transpose(np.asmatrix(h_1*environmental_vars_perm))
            h_2_interaction_mat = np.transpose(np.asmatrix(h_2*environmental_vars_perm))
        # Run permutation for all samples
        elif permutation_scheme == 'shuffle_all':
            environmental_vars_perm = permute_environmental_vars_all(environmental_vars)
            h_1_interaction_mat = np.transpose(np.asmatrix(h_1*environmental_vars_perm))
            h_2_interaction_mat = np.transpose(np.asmatrix(h_2*environmental_vars_perm))
        #  Fit null model. Use parameters from fitted null models to draw samples (gene counts).
        #  Then run LRT on sampled data
        elif permutation_scheme == 'sample_null':
            h_1_interaction_mat = np.transpose(np.asmatrix(h_1*environmental_vars))
            h_2_interaction_mat = np.transpose(np.asmatrix(h_2*environmental_vars))
            # Optimize the null model
            op_null = null_model_optimization_shell(null_data, sm, optimization_method)
            if model_version == 'as_log_linear':
                # Sample read counts from fitted null model
                ys, ns = draw_samples_from_fitted_null_as_model(np.atleast_1d(op_null['par']['conc']), op_null['par']['beta'], null_data)
                # Also update null data to have same allelic counts
                null_data['ys'] = ys
        else:
            print('error: permutation scheme ' + permutation_scheme + ' currently not implemented')
    elif permute == 'False':  # Do not permute the data
        h_1_interaction_mat = np.transpose(np.asmatrix(h_1*environmental_vars))
        h_2_interaction_mat = np.transpose(np.asmatrix(h_2*environmental_vars))
    if covariate_method == 'none':
        x_1 = np.hstack((intercept, env_mat, h_1_mat, h_1_interaction_mat))
        x_2 = np.hstack((intercept, env_mat, h_2_mat, h_2_interaction_mat))
    elif covariate_method == 't15_troponin':
        cov_mat = np.transpose(np.asmatrix(covs*environmental_vars))
        x_1 = np.hstack((intercept, env_mat, h_1_mat, cov_mat, h_1_interaction_mat))
        x_2 = np.hstack((intercept, env_mat, h_2_mat, cov_mat, h_2_interaction_mat))
    # Change of basis functions
    full_data = dict(N=N, K=K, T=T, P=x_1.shape[1], time_step=(environmental_vars.astype(int) + 1), x_1=x_1, x_2=x_2, ys=ys, ns=ns, concShape=1.001, concRate=0.001)
    return full_data, null_data




def run_dynamic_qtl(sm, null_data, full_data, dof, algorithm, iteration, model_version):
    # Use same seed for null and alternate models
    seed = np.random.randint(10000000) + 1

    # Run pystan gradient based optimization on full model
    op_full = sm.optimizing(data=full_data, as_vector=False, seed=seed, algorithm=algorithm, tol_obj=1e-17, tol_rel_obj=1e0, tol_grad=1e-13, tol_rel_grad=1e2, tol_param=1e-14)

    # Initialize null model with parameters defining the full model
    # initialization for joint model
    if 'conc' in op_full['par'] and 'nb_conc' in op_full['par']:
        init_null = dict(conc=np.atleast_1d(op_full['par']['conc']), nb_conc=op_full['par']['nb_conc'], beta=op_full['par']['beta'][:(len(op_full['par']['beta'])-dof)])
    # Initialization for te only model
    elif 'conc' not in op_full and 'nb_conc' in op_full['par']:
        init_null = dict(nb_conc=op_full['par']['nb_conc'], beta=op_full['par']['beta'][:(len(op_full['par']['beta'])-dof)])
    # Initialization for as only model
    elif 'conc' in op_full['par'] and 'nb_conc' not in op_full['par']:
        init_null = dict(conc=np.atleast_1d(op_full['par']['conc']), beta=op_full['par']['beta'][:(len(op_full['par']['beta'])-dof)])
    elif 'conc' not in op_full['par'] and 'nb_conc' not in op_full['par']:
        init_null = dict(beta=op_full['par']['beta'][:(len(op_full['par']['beta'])-dof)])

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
        # Refit the null with new random init
        refit_null = sm.optimizing(data=null_data, as_vector=False, seed=seed, algorithm=algorithm, tol_obj=1e-17, tol_rel_obj=1e0, tol_grad=1e-13, tol_rel_grad=1e2, tol_param=1e-14)
        # Use the null model that has the higher likelihood
        if refit_null['value'] > op_null['value']:
            op_null = refit_null
            loglr = op_full['value'] - op_null['value']
    return op_null, op_full, loglr



def dynamic_qtl(ys, ns, h_1, h_2, environmental_vars, model_version, permute, optimization_method, cell_line_indices, permutation_scheme, covs, covariate_method):
    # Load in correct model
    if model_version == 'as_log_linear':
        #sm = pystan.StanModel(file='as_log_linear.stan')
        #f = open('as_log_linear.pkl','wb')
        #pickle.dump(sm, f)
        #f.close()
        sm = pickle.load(open('as_log_linear.pkl', 'rb'))


    # Load in data into pystan format
    null_data = load_in_null_data(ys, ns, h_1, h_2, environmental_vars, model_version, covs, covariate_method)
    full_data, null_data = load_in_full_data(ys, ns, h_1, h_2, environmental_vars, permute, permutation_scheme, null_data, sm, cell_line_indices, optimization_method, model_version, covs, covariate_method)
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
    
    #if pvalue <= .00001:
    #    test_dicti = {}
    #    test_dicti['rs_id'] = 'rsid'
    #    test_dicti['ensamble_id'] = 'ensamble'
    #    test_dicti['ref_allele'] = 'A'
    #    test_dicti['alt_allele'] = 'T'
    #    test_dicti['pvalue'] = pvalue
    #    test_dicti['beta'] = op_full['par']['beta'][-1]
    #    test_dicti['conc'] = np.atleast_1d(op_full['par']['conc'])
    #    output_file = '/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data/temper_debug/' + model_version + '_'  + permutation_scheme + '_' + str(np.random.randint(10000000)) + '.png'
    #    visualize_hit(full_data['ys'], full_data['ns'], full_data['gene_counts'], h_1, h_2, environmental_vars, library_size_correction_factors, test_dicti, output_file)

    return dict(pvalue=pvalue, loglr=loglr, fit_full=op_full, fit_null=op_null, loglike_null=op_null['value'], loglike_full=op_full['value'])
