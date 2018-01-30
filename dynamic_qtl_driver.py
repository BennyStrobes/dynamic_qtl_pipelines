import numpy as np
import pystan
import pdb
from scipy import stats
import pickle


# Load in data into pystan format for null model (not including interaction term)
def load_in_null_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors):
    N = len(gene_counts)
    K = ys.shape[1]
    intercept = np.ones((N, 1))
    env_mat = np.transpose(np.asmatrix(environmental_vars))
    h_1_mat = np.transpose(np.asmatrix(h_1))
    h_2_mat = np.transpose(np.asmatrix(h_2))
    x_1 = np.hstack((intercept, env_mat, h_1_mat))
    x_2 = np.hstack((intercept, env_mat, h_2_mat))
    data = dict(N=N, K=K, P=x_1.shape[1], library_size=library_size_correction_factors, x_1=x_1, x_2=x_2, ys=ys, ns=ns, gene_counts=gene_counts, concShape=1.001, concRate=0.001)
    return data


def permute_environmental_vars_within_cell_line(environmental_vars, cell_line_indices):
    # Seperately permute all samples within each cell line
    for cell_line in cell_line_indices:
        # Get indices of samples that belong to this cell line
        index = cell_line_indices[cell_line]
        environmental_vars[index] = np.random.permutation(environmental_vars[index])
    return environmental_vars


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

# Load in data into pystan format for full model (including interaction term)
def load_in_full_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors, permute):
    N = len(gene_counts)
    K = ys.shape[1]
    intercept = np.ones((N,1))
    env_mat = np.transpose(np.asmatrix(environmental_vars))
    h_1_mat = np.transpose(np.asmatrix(h_1))
    h_2_mat = np.transpose(np.asmatrix(h_2))
    # Permute the data
    if permute == 'True':
        #environmental_vars_perm = permute_environmental_vars_within_cell_line(environmental_vars, cell_line_indices)
        environmental_vars_perm = permute_environmental_vars_within_hets(environmental_vars, h_1, h_2)
        h_1_interaction_mat = np.transpose(np.asmatrix(h_1*environmental_vars_perm))
        h_2_interaction_mat = np.transpose(np.asmatrix(h_2*environmental_vars_perm))
    elif permute == 'False':  # Do not permute the data
        h_1_interaction_mat = np.transpose(np.asmatrix(h_1*environmental_vars))
        h_2_interaction_mat = np.transpose(np.asmatrix(h_2*environmental_vars))
    x_1 = np.hstack((intercept, env_mat, h_1_mat, h_1_interaction_mat))
    x_2 = np.hstack((intercept, env_mat, h_2_mat, h_2_interaction_mat))
    data = dict(N=N, K=K, P=x_1.shape[1], library_size=library_size_correction_factors, x_1=x_1, x_2=x_2, ys=ys, ns=ns, gene_counts=gene_counts, concShape=1.001, concRate=0.001)
    return data


def run_dynamic_qtl(sm, null_data, full_data, dof, algorithm, iteration):
    # Use same seed for null and alternate models
    seed = np.random.randint(10000000) + 1

    # Run pystan gradient based optimization on null model
    op_full = sm.optimizing(data=full_data, as_vector=False, seed=seed, algorithm=algorithm, iter=300)

    # Initialize null model with parameters defining the full model
    init_null = dict(conc=np.atleast_1d(op_full['par']['conc']), nb_conc=op_full['par']['nb_conc'], beta=op_full['par']['beta'][:(len(op_full['par']['beta'])-dof)])

    # Run pystan gradient based optimization on null model
    if iteration == 1:
        op_null = sm.optimizing(data=null_data, as_vector=False, init=init_null, seed=seed, algorithm=algorithm, iter=300)
    else:  # Don't initialize with full model on later iterations
        op_null = sm.optimizing(data=null_data, as_vector=False, seed=seed, algorithm=algorithm, iter=300)

    # Compute chi-squared test statistic
    loglr = op_full['value'] - op_null['value']

    # Consider possibility that null did not fully converge
    if (loglr > 3):
        # Refit the null with new random init
        refit_null = sm.optimizing(data=null_data, as_vector=False, seed=seed, algorithm=algorithm, iter=300)
        # Use the null model that has the higher likelihood
        if refit_null['value'] > op_null['value']:
            op_null = refit_null
            loglr = op_full['value'] - op_null['value']
    return op_null, op_full, loglr


def dynamic_qtl(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors, model_version, permute):
    # Load in data into pystan format
    null_data = load_in_null_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors)
    full_data = load_in_full_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors, permute)
    # Calculate the degrees of freedom of LRT
    dof = full_data['P'] - null_data['P']

    # Load in correct model
    if model_version == 'joint_log_linear':
        # sm = pystan.StanModel(file='joint_log_linear.stan')
        # f = open('joint_log_linear.pkl','wb')
        # pickle.dump(sm, f)
        sm = pickle.load(open('joint_log_linear.pkl', 'rb'))

    # Run test by placing in try catch loop
    # If doesn't converge, then try it again with a different seed
    working = True
    iteration = 1
    while working:
        try:
            # Run dynamic qtls
            op_null, op_full, loglr = run_dynamic_qtl(sm, null_data, full_data, dof, 'Newton', iteration)
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
                print('LBFGS ran')
                op_null, op_full, loglr = run_dynamic_qtl(sm, null_data, full_data, dof, 'LBFGS', iteration)

    # Compute pvalue from chi-squared test statistic
    pvalue = 1.0 - stats.chi2.cdf(2.0*loglr, dof)

    return dict(pvalue=pvalue, loglr=loglr, fit_full=op_full, fit_null=op_null, loglike_null=op_null['value'], loglike_full=op_full['value'])
