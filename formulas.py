import sys, math
from numpy import multiply

sys.path.append('../lattice-estimator')

try:
    from estimator import *
except ImportError:
    print("Warning: Failed to import lattice_estimator, some options will not work")

try:
    import numpy as np
except ImportError:
    print("Warning: Failed to import numpy.")
    exit(0)

try:
    from scipy.special import lambertw
except ImportError:
    print("Warning: Failed to import scipy.")
    exit(0)



# Define global constants
PI = 3.14159265358979
e = 2.71828182845905
const = 2 * PI * e
ln2 = 0.693147180559945

# Auxiliary functions

# Low beta Gaussian Heuristic constant for use in NTRU Dense sublattice estimation.
gh_constant = {1: 0.00000, 2: -0.50511, 3: -0.46488, 4: -0.39100, 5: -0.29759, 6: -0.24880, 7: -0.21970, 8: -0.15748,
               9: -0.14673, 10: -0.07541, 11: -0.04870, 12: -0.01045, 13: 0.02298, 14: 0.04212, 15: 0.07014,
               16: 0.09205, 17: 0.12004, 18: 0.14988, 19: 0.17351, 20: 0.18659, 21: 0.20971, 22: 0.22728, 23: 0.24951,
               24: 0.26313, 25: 0.27662, 26: 0.29430, 27: 0.31399, 28: 0.32494, 29: 0.34796, 30: 0.36118, 31: 0.37531,
               32: 0.39056, 33: 0.39958, 34: 0.41473, 35: 0.42560, 36: 0.44222, 37: 0.45396, 38: 0.46275, 39: 0.47550,
               40: 0.48889, 41: 0.50009, 42: 0.51312, 43: 0.52463, 44: 0.52903, 45: 0.53930, 46: 0.55289, 47: 0.56343,
               48: 0.57204, 49: 0.58184, 50: 0.58852}

# Low beta \alpha_\beta quantity as defined in [AC:DucWoe21] for use in NTRU Dense subblattice estimation.
small_slope_t8 = {2: 0.04473, 3: 0.04472, 4: 0.04402, 5: 0.04407, 6: 0.04334, 7: 0.04326, 8: 0.04218, 9: 0.04237,
                  10: 0.04144, 11: 0.04054, 12: 0.03961, 13: 0.03862, 14: 0.03745, 15: 0.03673, 16: 0.03585,
                  17: 0.03477, 18: 0.03378, 19: 0.03298, 20: 0.03222, 21: 0.03155, 22: 0.03088, 23: 0.03029,
                  24: 0.02999, 25: 0.02954, 26: 0.02922, 27: 0.02891, 28: 0.02878, 29: 0.02850, 30: 0.02827,
                  31: 0.02801, 32: 0.02786, 33: 0.02761, 34: 0.02768, 35: 0.02744, 36: 0.02728, 37: 0.02713,
                  38: 0.02689, 39: 0.02678, 40: 0.02671, 41: 0.02647, 42: 0.02634, 43: 0.02614, 44: 0.02595,
                  45: 0.02583, 46: 0.02559, 47: 0.02534, 48: 0.02514, 49: 0.02506, 50: 0.02493, 51: 0.02475,
                  52: 0.02454, 53: 0.02441, 54: 0.02427, 55: 0.02407, 56: 0.02393, 57: 0.02371, 58: 0.02366,
                  59: 0.02341, 60: 0.02332}

# @cached_function
def ball_log_vol(n):
    return ((n/2.) * math.log(PI) - math.lgamma(n/2. + 1))

def log_gh(d, logvol=0):
    if d < 49:
        return (gh_constant[d] + logvol/d)

    return (1./d * (logvol - ball_log_vol(d)))

def delta(k):
    assert k >= 60
    delta = math.exp(log_gh(k)/(k-1))
    return (delta)

def check_overstreched(params):
    n    = params['n']
    lgq  = params['lgq']
    stds = params['std_s']
    stde = params['std_e']

    lnq = multiply(lgq, ln2)

    dense_det_log = math.log( math.sqrt(stds**2*n)+math.sqrt(stde**2*n))

    for beta in range(50, 1000):
        alpha_beta = delta(beta)**2
        alpha_beta_log = math.log(alpha_beta)
        m = math.round(0.5+lnq/(2*alpha_beta_log))

        rhs_log = 0.5*(m-1)*lnq-0.5*(m-1)**2*alpha_beta_log

        if dense_det_log<rhs_log:
            return beta

    return -1

def predicted_beta_bdd(params):

    n = params.n
    sigma = params.Xe.stddev
    q = params.q
    zeta = params.Xe.stddev / params.Xs.stddev

    lnq = math.log(q)

    beta_approx1 = 2*n/lnq*(math.log(n/(lnq))) #very rough approximation

    A = 2*n*lnq
    B = 2*(lnq - math.log(sigma*math.sqrt(const)))+math.log(beta_approx1/const)
    C = n/(2*lnq)
    D = lnq - math.log(zeta)

    Z  =  ((2*D*math.sqrt(C)+math.sqrt(A))/B)**2   #approximates beta/(ln(beta/2*const))
    ln_beta_const = -lambertw(-const/Z,-1) #approximates ln(beta/const)

    num   = 2*n*lnq*ln_beta_const
    denom = ln_beta_const+2*(lnq-math.log(sigma*math.sqrt(const))) - 2*math.sqrt(n/(2*lnq*Z))*(lnq-math.log(zeta))

    res = num/denom**2
    return res

def predicted_beta_usvp(d, lnq, sig, chi):
    x = np.divide(d, lnq)
    f1 = np.divide(np.multiply(d, np.log(x)), np.multiply(const, lnq - np.log(sig)))
    f2 = np.multiply(x, np.log(np.divide(d, lnq - np.log(sig))))

    # Beta calculation
    beta = np.divide(np.multiply(2 * d, np.multiply(lnq - np.log(chi), np.log(f1))), 
                     np.power(lnq + 0.5 * np.log(f2) - np.log(const * sig), 2))
    return beta

# Main functions

# Eq. (14)
def model_lambda_usvp(d, logq, std_s, std_e, params):
    
    sig = std_e
    chi = std_e/std_s

    lnq = np.multiply(logq, ln2)

    beta = predicted_beta_usvp(d, lnq, sig, chi)

    m2 = np.multiply(2 * d, beta * np.divide(lnq - np.log(chi), np.log(np.divide(beta, const))))

    return np.multiply(params[0], beta) + np.multiply(params[1], np.log(m2)) + params[2]

# Eq. (16)
def model_lambda_usvp_s(d, logq, params):
    lnq = np.multiply(logq, ln2)

    return np.multiply(np.multiply(params[0], np.divide(d, lnq)), 
                       np.log(np.divide(params[1] * d, lnq))) + np.multiply(params[2], np.log(d)) + params[3]

# Eq. (17)
def model_lambda_bdd(d, logq, std_s, std_e, std_s_num, params):
    sig = std_e 
    chi = std_e/std_s
    
    lnq = np.multiply(logq, ln2)

    beta = []

    if isinstance(logq, list):
        for lq in logq:
            params_est = LWE.Parameters(d, 2 ** lq, ND.UniformMod(std_s_num), ND.DiscreteGaussian(std_e))
            beta.append(predicted_beta_bdd(params_est))
    else:
        params_est = LWE.Parameters(d, 2 ** logq, ND.UniformMod(std_s_num), ND.DiscreteGaussian(std_e))
        beta.append(predicted_beta_bdd(params_est))

    # Intermediate calculations
    log_delta = np.log(beta) / (np.multiply(2, beta))
    m = np.sqrt(d * np.divide(lnq, log_delta))
    m2 = d * np.divide((lnq - np.log(chi)), log_delta)

    return np.multiply(params[0], beta) + params[1] * np.log(m2) + params[2]

# Eq. (20)
def model_lambda_bdd_s(d, logq, params):
    lnq = np.multiply(logq, ln2)

    return np.multiply(np.divide(params[0] * d, lnq), 
                       np.log(params[1] * d / lnq)) + np.multiply(params[2], np.log(d)) + params[3]

def model_n_bdd(l, logq, std_s, std_e, params):
    sigma = std_e
    zeta = std_e / std_s
    beta_approx = (l - np.log(100*l) - 16.4)/0.292 #approximate beta from lambda
    lnq = np.multiply(logq, ln2)
    n = l

    A = 2*lnq
    B = beta_approx
    C = np.log(beta_approx/const)
    D = lnq - np.log(zeta)
    E = 2*np.log(sigma*np.sqrt(const))

    denom = C*(A+2*D)**2
    nom = A*B*(A+C-E)**2
    return nom/denom

# Eq. (21)
def model_n_usvp(l, logq, params):
    lnq = np.multiply(logq, ln2)
    return np.multiply(np.divide(l + params[0] * np.log(lnq), params[1] * np.log(l) + params[2]) + params[3], lnq)

def model_n_usvp_s(l, logq, params):
    return "Not defined"

# Eq. (22)
def model_n_bdd_s(l, logq, std_s, std_e, params):
    sigma = std_e

    chi = std_e/std_s

    lnq = np.multiply(logq, ln2)

    # Rough approximation for beta
    beta_approx = l / 0.292
    num = (np.log(beta_approx / const) + lnq - np.log(const * sigma ** 2 / chi)) ** 2
    denom = 2 * (np.log(beta_approx / const) * lnq)
    leading_order = beta_approx * num / denom
    non_lead_order = 1 / (0.292) * np.log(8 * np.sqrt(leading_order * lnq * beta_approx / (beta_approx / const))) * num / denom

    return np.divide(np.multiply(l, params[0] * lnq), np.log(l)) + params[1] * np.log(lnq) * lnq + params[2] * lnq + params[3]