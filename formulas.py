import sys, math

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

def predicted_beta_(params):

    PI = 3.14159265358979
    e = 2.71828182845905
    const = 2*PI*e

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

# Main functions

# Eq. (14)
# TODO: add arbitrary error distribution
def model_lambda_usvp(d, logq, secret_dist, params):
    sig = 3.19 

    # Determine chi based on secret distribution
    if secret_dist == 'binary':
        chi = 2 * sig
    else:
        chi = np.sqrt(3 / 2) * sig

    lnq = np.multiply(logq, ln2)

    # Intermediate calculations
    x = np.divide(d, lnq)
    f1 = np.divide(np.multiply(d, np.log(x)), np.multiply(2 * PI * e, lnq - np.log(sig)))
    f2 = np.multiply(x, np.log(np.divide(d, lnq - np.log(sig))))

    # Beta calculation
    beta = np.divide(np.multiply(2 * d, np.multiply(lnq - np.log(chi), np.log(f1))), 
                     np.power(lnq + 0.5 * np.log(f2) - np.log(const * sig), 2))
    m2 = np.multiply(2 * d, beta * np.divide(lnq - np.log(chi), np.log(np.divide(beta, const))))

    return np.multiply(params[0], beta) + np.multiply(params[1], np.log(m2)) + params[2]

# Eq. (16)
def model_lambda_usvp_s(d, logq, params):
    lnq = np.multiply(logq, ln2)

    return np.multiply(np.multiply(params[0], np.divide(d, lnq)), 
                       np.log(np.divide(params[1] * d, lnq))) + np.multiply(params[2], np.log(d)) + params[3]

# Eq. (17)
# TODO: add arbitrary error distribution
def model_lambda_bdd(d, logq, secret_dist, params):
    sig = 3.19 

    if secret_dist == 'binary':
        chi = 2 * sig
    else:
        chi = np.sqrt(3 / 2) * sig

    lnq = np.multiply(logq, ln2)

    beta = []

    if isinstance(logq, list):
        for lq in logq:
            params_est = LWE.Parameters(d, 2 ** lq, ND.UniformMod(3), ND.DiscreteGaussian(3.19))
            beta.append(predicted_beta_(params_est))
    else:
        params_est = LWE.Parameters(d, 2 ** logq, ND.UniformMod(3), ND.DiscreteGaussian(3.19))
        beta.append(predicted_beta_(params_est))

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

# Eq. (21)
def model_n_usvp(l, logq, params):
    lnq = np.multiply(logq, ln2)
    return np.multiply(np.divide(l + params[0] * np.log(lnq), params[1] * np.log(l) + params[2]) + params[3], lnq)

# Eq. (22)
# TODO: add arbitrary error distribution
def model_n_bdd(l, logq, secret_dist, params):
    sigma = 3.19
    if secret_dist == 'binary':
        chi = 2 * sigma
    else:
        chi = np.sqrt(3 / 2) * sigma

    lnq = np.multiply(logq, ln2)

    # Rough approximation for beta
    beta_approx = l / 0.292
    num = (np.log(beta_approx / const) + lnq - np.log(const * sigma ** 2 / chi)) ** 2
    denom = 2 * (np.log(beta_approx / const) * lnq)
    leading_order = beta_approx * num / denom
    non_lead_order = 1 / (0.292) * np.log(8 * np.sqrt(leading_order * lnq * beta_approx / (beta_approx / const))) * num / denom

    return np.divide(np.multiply(l, params[0] * lnq), np.log(l)) + params[1] * np.log(lnq) * lnq + params[2] * lnq + params[3]
