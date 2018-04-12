import GPy
import numpy as np
import pystan
import pdb
from scipy import stats
import pickle
import pandas as pd



# Load in data into pystan format for null model (not including interaction term)
def load_in_null_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors, as_overdispersion_parameter, as_overdispersion_parameter_sample_specific, model_version, covs, covariate_method, te_nb_conc):
    N = len(gene_counts)
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
    # Change of basis functions
    if model_version == 'te_log_linear_quadratic_basis':
        x_1 = np.hstack((x_1, np.square(x_1[:,1])))
        x_2 = np.hstack((x_2, np.square(x_2[:,1])))
    elif model_version == 'te_log_linear_cubic_control':
        x_1 = np.hstack((x_1, np.power(x_1[:,1], 2)))
        x_2 = np.hstack((x_2, np.power(x_2[:,1], 2)))
    data = dict(N=N, K=K, T=T, nb_conc=te_nb_conc, time_step=(environmental_vars.astype(int) + 1), P=x_1.shape[1], library_size=library_size_correction_factors, x_1=x_1, x_2=x_2, ys=ys, ns=ns, gene_counts=gene_counts, concShape=1.001, concRate=0.001, as_overdispersion_parameter=as_overdispersion_parameter, as_overdispersion_parameter_sample_specific=as_overdispersion_parameter_sample_specific)
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
            op_null = sm.optimizing(data=null_data, as_vector=False, seed=seed, algorithm=optimization_method, tol_obj=1e-15, tol_rel_obj=1e1, tol_grad=1e-11, tol_rel_grad=1e4, tol_param=1e-12)
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
                op_null = sm.optimizing(data=null_data, as_vector=False, seed=seed, algorithm='LBFGS', tol_obj=1e-15, tol_rel_obj=1e1, tol_grad=1e-11, tol_rel_grad=1e4, tol_param=1e-12)
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




# Sample read counts from fitted null model
def draw_samples_from_fitted_null_joint_model(conc, nb_conc, beta, null_data):
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

    if np.sum(np.isnan(allelic_p)) > 0:
        print('nan p')
        print(allelic_p)
        print(xb_1)
        print(xb_2)
        print(beta)
        print(null_data['x_1'])
        print(null_data['x_2'])

    # Sample total reads from negative binomial distribution
    sample_gene_counts = sample_gene_counts_from_negative_binomial_distribution(y_1, y_2, nb_conc, null_data['library_size'])

    # No heterozygous sites exist
    if null_data['ys'].shape[1] == 0:
        sample_ys = null_data['ys']
        sample_ns = null_data['ns']
    else:  # heterozygous sites exist
        # Sample allelic counts from beta-binomial distribution
        sample_ys, sample_ns = sample_allelic_counts_from_beta_binomial_distribution(allelic_p, conc, null_data['ys'], null_data['ns'])

    return np.asarray(sample_gene_counts).astype(int), sample_ys.astype(int), sample_ns.astype(int)


# Sample read counts from fitted null model
def draw_samples_from_fitted_null_te_model(nb_conc, beta, null_data):
    beta_mat = np.transpose(np.asmatrix(beta))
    # Predict relative log number of reads on allele 1
    xb_1 = np.dot(null_data['x_1'], beta_mat)
    # Predict relative log number of reads on allele 2
    xb_2 = np.dot(null_data['x_2'], beta_mat)
    # Predict relative number of reads on allele 1
    y_1 = np.exp(xb_1)
    # Predict relative number of reads on allele 2
    y_2 = np.exp(xb_2)

    # Sample total reads from negative binomial distribution
    sample_gene_counts = sample_gene_counts_from_negative_binomial_distribution(y_1, y_2, nb_conc, null_data['library_size'])

    return np.asarray(sample_gene_counts).astype(int)

# Load in data into pystan format for full model (including interaction term)
def load_in_full_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors, permute, permutation_scheme, null_data, sm, cell_line_indices,optimization_method, model_version, as_overdispersion_parameter,as_overdispersion_parameter_sample_specific, covs, covariate_method, te_nb_conc):
    N = len(gene_counts)
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
            if model_version == 'joint_log_linear':
                # Sample read counts from fitted null model
                gene_counts, ys, ns = draw_samples_from_fitted_null_joint_model(np.atleast_1d(op_null['par']['conc']), np.atleast_1d(op_null['par']['nb_conc'])[0], op_null['par']['beta'], null_data)
                # Also update null data to have same allelic counts
                null_data['ys'] = ys
                null_data['gene_counts'] = gene_counts
            elif model_version == 'te_log_linear' or model_version == 'te_log_linear_quadratic_basis' or model_version == 'te_log_linear_cubic_control':
                gene_counts = draw_samples_from_fitted_null_te_model(np.atleast_1d(op_null['par']['nb_conc'])[0], op_null['par']['beta'], null_data)
                null_data['gene_counts'] = gene_counts
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
    if model_version == 'te_log_linear_quadratic_basis':
        x_1 = np.hstack((x_1[:, 0:3], np.square(x_1[:, 1]), x_1[:, -1], np.square(x_1[:, -1])))
        x_2 = np.hstack((x_2[:, 0:3], np.square(x_2[:, 1]), x_2[:, -1], np.square(x_2[:, -1])))
    elif model_version == 'te_log_linear_cubic_control':
        x_1 = np.hstack((null_data['x_1'], h_1_interaction_mat))
        x_2 = np.hstack((null_data['x_2'], h_2_interaction_mat))
    full_data = dict(N=N, nb_conc=te_nb_conc, K=K, T=T, P=x_1.shape[1], time_step=(environmental_vars.astype(int) + 1), library_size=library_size_correction_factors, x_1=x_1, x_2=x_2, ys=ys, ns=ns, gene_counts=gene_counts, concShape=1.001, concRate=0.001, as_overdispersion_parameter=as_overdispersion_parameter, as_overdispersion_parameter_sample_specific=as_overdispersion_parameter_sample_specific)
    return full_data, null_data


def quadratic_basis_initialize_optimization(full_data, seed, algorithm, sm):
    temp_x_1 = np.hstack((full_data['x_1'][:, 0:3], full_data['x_1'][:, 4]))
    temp_x_2 = np.hstack((full_data['x_2'][:, 0:3], full_data['x_2'][:, 4]))
    temp_data = dict(x_1=temp_x_1, x_2=temp_x_2, P=temp_x_1.shape[1], N=full_data['N'], K=full_data['K'], library_size=full_data['library_size'], ys=full_data['ys'], ns=full_data['ns'], gene_counts=full_data['gene_counts'], concShape=full_data['concShape'], concRate=full_data['concRate'], as_overdispersion_parameter=full_data['as_overdispersion_parameter'], as_overdispersion_parameter_sample_specific=full_data['as_overdispersion_parameter_sample_specific'])
    op_temp = sm.optimizing(data=temp_data, as_vector=False, seed=seed, algorithm=algorithm, tol_obj=1e-15, tol_rel_obj=1e1, tol_grad=1e-11, tol_rel_grad=1e4, tol_param=1e-12)
    return op_temp

def cubic_control_initialize_optimization(full_data, seed, algorithm, sm):
    temp_x_1 = full_data['x_1'][:,0:3]
    temp_x_2 = full_data['x_2'][:,0:3]
    temp_data = dict(x_1=temp_x_1, x_2=temp_x_2, P=temp_x_1.shape[1], N=full_data['N'], K=full_data['K'], library_size=full_data['library_size'], ys=full_data['ys'], ns=full_data['ns'], gene_counts=full_data['gene_counts'], concShape=full_data['concShape'], concRate=full_data['concRate'], as_overdispersion_parameter=full_data['as_overdispersion_parameter'], as_overdispersion_parameter_sample_specific=full_data['as_overdispersion_parameter_sample_specific'])
    op_temp = sm.optimizing(data=temp_data, as_vector=False, seed=seed, algorithm=algorithm, tol_obj=1e-15, tol_rel_obj=1e1, tol_grad=1e-11, tol_rel_grad=1e4, tol_param=1e-12)
    return op_temp


def run_dynamic_qtl(sm, null_data, full_data, dof, algorithm, iteration, model_version):
    # Use same seed for null and alternate models
    seed = np.random.randint(10000000) + 1

    if model_version == 'te_log_linear_quadratic_basis':
        # First optimize model on simpler version first. And use simpler model as initialization
        op_temp = quadratic_basis_initialize_optimization(full_data, seed, algorithm, sm)
        init_full = dict(nb_conc=np.atleast_1d(op_temp['par']['nb_conc'])[0], beta = np.asarray([op_temp['par']['beta'][0], op_temp['par']['beta'][1], op_temp['par']['beta'][2], 0.0, op_temp['par']['beta'][3], 0]))
        op_full = sm.optimizing(data=full_data, as_vector=False, init=init_full, seed=seed, algorithm=algorithm, tol_obj=1e-15, tol_rel_obj=1e1, tol_grad=1e-11, tol_rel_grad=1e4, tol_param=1e-12)
    elif model_version == 'te_log_linear_cubic_control':
        op_temp = cubic_control_initialize_optimization(full_data, seed, algorithm, sm)
        init_full = dict(nb_conc=np.atleast_1d(op_temp['par']['nb_conc'])[0], beta = np.asarray([op_temp['par']['beta'][0], op_temp['par']['beta'][1], op_temp['par']['beta'][2], 0.0, 0.0]))
        op_full = sm.optimizing(data=full_data, as_vector=False, init=init_full, seed=seed, algorithm=algorithm, tol_obj=1e-15, tol_rel_obj=1e1, tol_grad=1e-11, tol_rel_grad=1e4, tol_param=1e-12)
    else:
        # Run pystan gradient based optimization on full model
        op_full = sm.optimizing(data=full_data, as_vector=False, seed=seed, algorithm=algorithm, tol_obj=1e-15, tol_rel_obj=1e1, tol_grad=1e-11, tol_rel_grad=1e4, tol_param=1e-12)

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

# PLOT!
def visualize_hit(ys, ns, gene_counts, h_1, h_2, environmental_vars, library_size_correction_factors, test_dicti, output_file):
    plt.clf()
    # Plot to visualize total expression changes over time as a function of genotype
    fig = plt.figure(figsize=(40, 20))

    # call regplot on each axes
    ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=4)
    ax2 = plt.subplot2grid((2, 4), (1, 0))
    ax3 = plt.subplot2grid((2, 4), (1, 1))
    ax4 = plt.subplot2grid((2, 4), (1, 2))
    ax5 = plt.subplot2grid((2, 4), (1, 3))
    #sns.regplot(x=idx, y=df['x'], ax=ax1)
    #sns.regplot(x=idx, y=df['y'], ax=ax2)
    num_sites = ys.shape[1]
    if num_sites > 4:
        num_sites = 4

    gene_total_plot = gene_total_plotter(gene_counts, h_1+h_2, environmental_vars, test_dicti['rs_id'], test_dicti['ensamble_id'], test_dicti['ref_allele'], test_dicti['alt_allele'], test_dicti['pvalue'], test_dicti['beta'], library_size_correction_factors,ax1)
    #fig.savefig(output_file)
    num_sites = ys.shape[1]
    if num_sites > 4:
        num_sites = 4
    for exonic_site_num in np.arange(num_sites):
        if exonic_site_num == 0:
            axy = ax2
        elif exonic_site_num == 1:
            axy = ax3
        elif exonic_site_num == 2:
            axy = ax4
        elif exonic_site_num == 3:
            axy = ax5
        allelic_imbalence_plot = allelic_imbalence_plotter(h_1, h_2, environmental_vars, ys, ns, test_dicti['rs_id'], test_dicti['ensamble_id'], test_dicti['pvalue'], test_dicti['beta'], test_dicti['conc'].astype(float),axy, exonic_site_num)
    plt.tight_layout()
    fig.savefig(output_file)


# Convert from dosage vector to vector of genotyeps
def dosage_to_genotype(dosage, ref_allele, alt_allele):
    converter = {}
    #converter[0] = ref_allele + ref_allele
    #converter[1] = ref_allele + alt_allele
    #converter[2] = alt_allele + alt_allele
    converter[0] = "0"
    converter[1] = "1"
    converter[2] = "2"


    genotype = []
    for dos in dosage:
        genotype.append(converter[dos])
    return np.asarray(genotype)

# Plot to visualize total expression changes over time as a function of genotype
def gene_total_plotter(gene_counts, dosage, environmental_vars, rs_id, ensamble_id, ref_allele, alt_allele, pvalue, beta, library_size_correction_factors,ax1):
    # Convert from dosage vector to vector of genotyeps
    genotype = dosage_to_genotype(dosage, ref_allele, alt_allele)

    gene_counts = (gene_counts/library_size_correction_factors)*np.mean(library_size_correction_factors)

    df = pd.DataFrame({rs_id: genotype, 'time_step': environmental_vars.astype(int), 'gene_counts': np.log(gene_counts)})
    ax = sns.boxplot(x="time_step", y="gene_counts", hue=rs_id, data=df, palette="Set3",width=.7,ax=ax1)
    plt.xlabel('Time Step')
    plt.ylabel('log(counts)')
    ax1.set_title(ensamble_id + ' / pvalue = ' + str(pvalue) + ' / beta = ' + str(beta))
    #sns.despine(offset=1, trim=True)
    return ax

def allelic_imbalence_plotter(h_1, h_2, environmental_vars, ys, ns, rs_id, ensamble_id, pvalue, beta, conc,axy, exonic_site_num):
    if ys.shape[1] != len(conc):
        print('fatal errooror')
        pdb.set_trace()
    # Order sites by smallest variance to largest
    ys = ys[:, conc.argsort()[::-1]]
    ns = ns[:, conc.argsort()[::-1]]
    # Extract allelic fractions at each site
    time_steps = []
    allelic_fractions = []
    cell_lines_arr = []
    identifiers = []
    exonic_sites = []
    depths = []
    num_exonic_sites = ys.shape[1]
    num_samples = len(environmental_vars)
    for sample_num in np.arange(num_samples):
        if h_1[sample_num] == h_2[sample_num]:  # homozygous variant
            continue
        if ns[sample_num, exonic_site_num] <= 2:
            continue
        if h_1[sample_num] == 0:
            allelic_fraction = 1.0 - float(ys[sample_num, exonic_site_num])/ns[sample_num, exonic_site_num]
        elif h_2[sample_num] == 0:
            allelic_fraction = (float(ys[sample_num, exonic_site_num])/ns[sample_num, exonic_site_num])
        else:
            print('eroroororo')
        depths.append(ns[sample_num,exonic_site_num])
        allelic_fractions.append(allelic_fraction)
        #allelic_fractions.append(abs(allelic_fraction-.5))
        time_steps.append(float(environmental_vars[sample_num]))

    # PLOT!
    df = pd.DataFrame({'time_step': np.asarray(time_steps).astype(int),'read_depth':depths,  'allelic_fraction': allelic_fractions})

    #ax = sns.pointplot(x="time_step", y="allelic_fraction", hue="identifiers", data=df)
    #ax = sns.regplot(x="time_step", y="allelic_fraction", data=df)
    ax = sns.regplot(data=df,x="time_step", y="allelic_fraction",ci=None, ax=axy)
    ax.set_title("Exonic site = " + str(exonic_site_num) + " / conc = " + str(conc[conc.argsort()[::-1]][exonic_site_num]))
    plt.ylim(ymax=1) # adjust the max leaving min unchanged
    plt.ylim(ymin=0)
    #ax = sns.boxplot(x="time_step", y="allelic_fraction", hue=exonic_sites,data=df, palette="Set3",width=.7)
    #plt.xlabel('Time Step')
    #plt.ylabel('Allelic fraction')
    #plt.title(ensamble_id + ' / pvalue = ' + str(pvalue) + ' / beta = ' + str(beta))
    #sns.despine(offset=1, trim=True)
    #iris = sns.load_dataset("iris")
    # Plot tip as a function of toal bill across days
    # = sns.regplot(x="sepal_length", y="sepal_width", data=iris)
    return ax

def run_gp_regression(y, x):
    np.random.seed(2)
    dim = x.shape[1]
    kernel = GPy.kern.RBF(1, active_dims=[0])
    for i in range(1, dim):
        kernel = kernel + GPy.kern.RBF(1,active_dims=[i])
    kernel = kernel + GPy.kern.White(dim) + GPy.kern.Bias(input_dim=dim)
    #kernel = GPy.kern.RBF(dim)
    #kernel = GPy.kern.RBF(input_dim=dim) #+ GPy.kern.White(dim) #+ GPy.kern.Bias(input_dim=dim)
    poisson_likelihood = GPy.likelihoods.Poisson()
    laplace_inf = GPy.inference.latent_function_inference.Laplace()
    m = GPy.core.GP(X=x, Y=y, likelihood=poisson_likelihood, inference_method=laplace_inf, kernel=kernel)
    m.randomize()
    m.optimize()
    return m


def dynamic_qtl(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors, model_version, permute, optimization_method, cell_line_indices, permutation_scheme, as_overdispersion_parameter, as_overdispersion_parameter_sample_specific, covs, covariate_method, te_nb_conc):
    # Load in correct model
    if model_version == 'joint_log_linear':
        #sm = pystan.StanModel(file='joint_log_linear.stan')
        #f = open('joint_log_linear.pkl','wb')
        #pickle.dump(sm, f)
        sm = pickle.load(open('joint_log_linear.pkl', 'rb'))
    elif model_version == 'te_log_linear':
        sm = pickle.load(open('te_log_linear.pkl', 'rb'))
    elif model_version == 'as_log_linear':
        sm = pickle.load(open('as_log_linear.pkl', 'rb'))
    elif model_version == 'as_log_linear_fixed_overdispersion':
        sm = pickle.load(open('as_log_linear_fixed_overdispersion.pkl', 'rb'))
    elif model_version == 'as_log_linear_fixed_sample_overdispersion':
        sm = pickle.load(open('as_log_linear_fixed_sample_overdispersion.pkl', 'rb'))
    elif model_version == 'te_log_linear_quadratic_basis':
        sm = pickle.load(open('te_log_linear.pkl', 'rb'))
    elif model_version == 'te_gaussian_process':
        sm = 'fake'
    elif model_version == 'te_log_linear_cubic_control':
        sm = pickle.load(open('te_log_linear.pkl', 'rb'))
    elif model_version == 'te_log_linear_time_od':
        sm = pickle.load(open('te_log_linear_time_od.pkl', 'rb'))
    elif model_version == 'te_log_linear_time_od_offline':
        sm = pickle.load(open('te_log_linear_time_od_offline.pkl', 'rb'))

    # Load in data into pystan format
    null_data = load_in_null_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors, as_overdispersion_parameter, as_overdispersion_parameter_sample_specific, model_version, covs, covariate_method, te_nb_conc)
    full_data, null_data = load_in_full_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors, permute, permutation_scheme, null_data, sm, cell_line_indices, optimization_method, model_version, as_overdispersion_parameter, as_overdispersion_parameter_sample_specific, covs, covariate_method, te_nb_conc)
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
