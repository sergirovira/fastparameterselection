from config import *
from scipy.special import lambertw

import sys
sys.path.append('../lattice-estimator')

from estimator import *

# Define global constants
PI = 3.14159265358979
e = 2.71828182845905
const = 2 * PI * e
ln2 = 0.693147180559945

# Auxiliary functions

def predicted_beta_(params):
    # Extract parameters
    n = params.n
    sigma = params.Xe.stddev
    q = params.q
    zeta = params.Xe.stddev / params.Xs.stddev
    lnq = log(q)

    # Rough approximation for beta
    beta_approx1 = n / lnq * (log(n / lnq))

    # Intermediate calculations
    A = 2 * n * lnq
    B = 2 * log(q / (sigma * sqrt(const))) + log(beta_approx1 / const)
    C = n / (2 * lnq)
    D = lnq - ln(zeta)
    Z = (((2 * D * sqrt(C) + sqrt(A)) / B) ** 2).n()

    # Exact value of ln(beta_const) using Lambert W function
    ln_beta_const = -lambert_w(-1, -const / Z).n()

    num = 2 * n * lnq * ln_beta_const
    denom = ln_beta_const + 2 * log(q / (sigma * sqrt(const))) - 2 * sqrt(n / (2 * lnq * Z)) * log(q / zeta)

    res = num / denom ** 2

    return CC(res)  # Return the result

# Main functions

# Eq. (14)
def model_lambda_usvp(d, logq, secret_dist, params):
    sig = 3.19 

    # Determine chi based on secret distribution
    if secret_dist == 'binary':
        chi = 2 * sig
    else:
        chi = np.sqrt(3 / 2) * sig

    # Calculate lnq
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
def model_lambda_bdd(d, logq, secret_dist, params):

    lnq = np.multiply(logq, ln2)

    beta = []

    if isinstance(logq, list):
        for lq in logq:
            params_est = LWE.Parameters(d, 2 ** lq, ND.UniformMod(2), ND.DiscreteGaussian(secret_dist[0]))
            chi = params_est.Xe.stddev / params_est.Xs.stddev
            print("chi", chi)
            beta.append(predicted_beta_(params_est))
    else:
        params_est = LWE.Parameters(d, 2 ** logq, ND.UniformMod(2), ND.DiscreteGaussian(secret_dist[0]))
        beta.append(predicted_beta_(params_est))

    print("Xe: ", params_est.Xe.stddev, "Xs", params_est.Xs.stddev, "chi: ", chi)

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
    sigma = 3.19
    eta = np.sqrt(3 / 2) * 3.19

    lnq = np.divide(logq, 1.44269504088896340735992468100)

    return np.multiply(np.divide(l + params[0] * np.log(lnq), params[1] * np.log(l) + params[2]) + params[3], lnq)

# Eq. (22)
def model_n_bdd(l, logq, params):
    sigma = 3.19
    chi = np.sqrt(3 / 2) * 3.19
    lnq = np.divide(logq, 1.44269504088896340735992468100)

    # Rough approximation for beta
    beta_approx = l / 0.292
    num = (np.log(beta_approx / const) + lnq - np.log(const * sigma ** 2 / chi)) ** 2
    denom = 2 * (np.log(beta_approx / const) * lnq)
    leading_order = beta_approx * num / denom
    non_lead_order = 1 / (0.292) * np.log(8 * np.sqrt(leading_order * lnq * beta_approx / (beta_approx / const))) * num / denom

    return np.divide(np.multiply(l, params[0] * lnq), np.log(l)) + params[1] * np.log(lnq) * lnq + params[2] * lnq + params[3]

def model_n_bdd_s(l, logq, params):
    lnq = np.divide(logq,1.44269504088896340735992468100)
    return np.multiply(np.divide(params[0]*l + params[1], np.log(l)+ params[2]* np.log(np.log(l))),lnq) 

def model_n_usvp_s(l, logq, params):
    lnq = np.divide(logq,1.44269504088896340735992468100)
    chi = np.sqrt(3 / 2) * 3.19

    logdelta = np.log(l/0.292/const)/(2*l/0.292)
    D = (np.log(np.sqrt(l/0.292))+lnq-np.log(const*sigma))/(2*logdelta)
    beta_approx = (l/0.292- params[0]*np.log2(8*D)-(1/0.292)*16.4)+params[1] 
    num   = (0.5*np.log(beta_approx)+lnq-np.log(const*sigma))**2
    denom = 2*(np.log(beta_approx/const)*(lnq-np.log(chi)))
    leading_order = beta_approx*num / denom
    return leading_order #+params[2] 
