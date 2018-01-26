import numpy as np
import pystan
import pdb
from scipy import stats
import pickle



# Load in data into pystan format for null model (not including interaction term)
def load_in_null_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors):
    N = len(gene_counts)
    K = ys.shape[1]
    intercept = np.ones((N,1))
    env_mat = np.transpose(np.asmatrix(environmental_vars))
    h_1_mat = np.transpose(np.asmatrix(h_1))
    h_2_mat = np.transpose(np.asmatrix(h_2))
    x_1 = np.hstack((intercept, env_mat, h_1_mat))
    x_2 = np.hstack((intercept, env_mat, h_2_mat))
    data = dict(N=N, K=K, P=x_1.shape[1], library_size=library_size_correction_factors, x_1=x_1, x_2=x_2, ys=ys, ns=ns, gene_counts=gene_counts, concShape=1.01, concRate=0.01)
    return data

# Load in data into pystan format for full model (including interaction term)
def load_in_full_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors):
    N = len(gene_counts)
    K = ys.shape[1]
    intercept = np.ones((N,1))
    env_mat = np.transpose(np.asmatrix(environmental_vars))
    h_1_mat = np.transpose(np.asmatrix(h_1))
    h_2_mat = np.transpose(np.asmatrix(h_2))
    h_1_interaction_mat = np.transpose(np.asmatrix(h_1*environmental_vars))
    h_2_interaction_mat = np.transpose(np.asmatrix(h_2*environmental_vars))
    x_1 = np.hstack((intercept, env_mat, h_1_mat, h_1_interaction_mat))
    x_2 = np.hstack((intercept, env_mat, h_2_mat, h_2_interaction_mat))
    data = dict(N=N, K=K, P=x_1.shape[1], library_size=library_size_correction_factors, x_1=x_1, x_2=x_2, ys=ys, ns=ns, gene_counts=gene_counts, concShape=1.01, concRate=0.01)
    return data

def run_dynamic_qtl(sm, null_data, full_data, dof, algorithm):
    # Use same seed for null and alternate models
    seed = np.random.randint(10000000) + 1

    # Run pystan gradient based optimization on null model
    op_null = sm.optimizing(data=null_data, as_vector=False, seed=seed, algorithm=algorithm)

    # Use null betas as initialization for alternate/full betas
    beta_init = np.hstack((op_null['par']['beta'], np.zeros(dof)))

    # init = dict(beta=beta_init, conc=op_null['par']['conc'], nb_conc=op_null['par']['nb_conc'])
    init = dict(beta=beta_init)

    # Run pystan gradient based optimization on null model
    op_full = sm.optimizing(data=full_data, as_vector=False, seed=seed, algorithm=algorithm)
    return op_null, op_full


def dynamic_qtl(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors, model_version):
    # Load in data into pystan format
    null_data = load_in_null_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors)
    full_data = load_in_full_data(gene_counts, ys, ns, h_1, h_2, environmental_vars, library_size_correction_factors)
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
            op_null, op_full = run_dynamic_qtl(sm, null_data, full_data, dof, 'Newton')
            working = False
        except:
            print('Starting over for ' + str(iteration) + ' time')
            iteration = iteration + 1
            # If you've tried 50 different start seeds and none of them converged, then do LBFGS
            if iteration > 50:
                print('LBFGS ran')
                op_null, op_full = run_dynamic_qtl(sm, null_data, full_data, dof, 'LBFGS')

    # Compute chi-squared test statistic
    loglr = op_full['value'] - op_null['value']

    # Compute pvalue from chi-squared test statistic
    pvalue = 1.0 - stats.chi2.cdf(2.0*loglr, dof)

    return dict(pvalue=pvalue, loglr=loglr, fit_full=op_full, fit_null=op_null, loglike_null=op_null['value'], loglike_full=op_full['value'])
