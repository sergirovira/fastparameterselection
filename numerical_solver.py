from numpy import pi, exp, log, log2, sqrt, multiply, divide
from scipy.optimize import fsolve, brentq, brenth, bisect, newton

import numpy as np

try:
    from estimator import *
except ImportError:
    print("Warning: Failed to import lattice_estimator, some options will not work")
    estimator_installed = 0

PI = 3.14159265358979
e = 2.71828182845905
const = 2 * PI * e
ln2 = 0.693147180559945

def numerical_beta_bdd(n, logq, std_s, std_e):
    beta_initial_guess = n / 4

    lnq = multiply(logq, ln2)

    zeta = std_s/std_e

    nom = lambda beta : 2 * n * lnq * log(beta/const)
    denom = lambda beta : log(beta/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
    eq = lambda beta : beta - nom(beta) / (denom(beta)**2)

    beta_solution = fsolve(eq, beta_initial_guess, full_output = False) # xtol = 10**(-15), factor = 1

    return beta_solution

def numerical_beta_usvp(n, logq, std_s, std_e):
    beta_initial_guess = n / 4

    zeta = std_s/std_e

    lnq = multiply(logq, ln2)

    nom = lambda beta : 2 * n * (lnq - log(zeta)) * log(beta/const)
    denom = lambda beta : lnq + log(sqrt(beta)/(const*std_e))
    eq = lambda beta : beta - nom(beta) / (denom(beta)**2)

    beta_solution = fsolve(eq, beta_initial_guess, full_output = False)

    return beta_solution

#This function takes as input the n computed by model_lambda_bdd
def numerical_n_bdd(n, logq,  std_s, std_e):
    # Use the fsolver to find n
    n_initial_guess = 1

    lnq = multiply(logq, ln2)

    zeta = std_s/std_e

    beta = numerical_beta_bdd(n, logq, std_s, std_e)

    nom = lambda n : 2 * n * lnq * log(beta/const)
    denom = lambda n : log(beta/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
    eq = lambda n : beta - nom(n) / (denom(n)**2)

    n_solution = fsolve(eq, n_initial_guess, full_output = False)
    
    return n_solution[0]

#This function takes as input the n computed by model_lambda_usvp
def numerical_n_usvp(n, logq,  std_s, std_e):
    n_initial_guess = 10

    lnq = multiply(logq, ln2)

    zeta = std_s/std_e

    beta = numerical_beta_usvp(n, logq, std_s, std_e)

    nom = lambda n : 2 * n * (lnq - log(zeta)) * log(beta/const)
    denom = lambda n : lnq + log(sqrt(beta)/(const*std_e))
    eq = lambda n : beta - nom(n) / (denom(n)**2)

    n_solution = fsolve(eq, n_initial_guess, full_output = False)
    
    return n_solution[0]

def numerical_logq_bdd(l, n, std_s, std_e):
    lnq_initial_guess = 3 * log(n)

    zeta = std_s/std_e

    beta = l/0.292 - 16.4

    nom = lambda lnq : 2 * n * lnq * log(beta/const)
    denom = lambda lnq : log(beta/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
    eq = lambda lnq : beta - nom(lnq) / (denom(lnq)**2)

    lnq_solution = fsolve(eq, lnq_initial_guess, full_output = False)
    
    return divide(lnq_solution[0], ln2)

def numerical_logq_usvp(l, n, std_s, std_e):
    lnq_initial_guess = 3 * log(n)

    zeta = std_s/std_e

    beta = l/0.292 - 16.4

    nom = lambda lnq : 2 * n * (lnq - log(zeta)) * log(beta/const)
    denom = lambda lnq : lnq + log(sqrt(beta)/(const*std_e))
    eq = lambda lnq : beta - nom(lnq) / (denom(lnq)**2)

    lnq_solution = fsolve(eq, lnq_initial_guess, full_output = False)
    
    return divide(lnq_solution[0], ln2)


def numerical_std_e_bdd(l, n, logq, std_s):
    lnq = logq * log(2)

    # find beta numerically
    beta_initial_guess = (l - 16.4) / 0.292
    d_optimal = lambda beta : sqrt(2 * n * lnq * beta / log(beta / const))
    eq8 = lambda beta : l - (0.292 * beta + log2(8 * d_optimal(beta)) + 16.4)

    beta_solution = fsolve(eq8, beta_initial_guess, full_output = False)
    #print("Beta initial guess: ", beta_initial_guess)
    #print("Beta: ", beta_solution)

    # find std_e numerically
    zeta_initial_guess = 10
    std_e_initial_guess = zeta_initial_guess * std_s

    std_e_initial_guess = 2

    nom = 2 * n * lnq * log(beta_solution/const)
    denom = lambda std_e : log(beta_solution/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(max(1,np.round(std_e/std_s)))) * sqrt(n * log(beta_solution/const) / (2 * lnq * beta_solution))
    eq6 = lambda std_e : beta_solution - nom / (denom(std_e)**2)

    std_e_solution = fsolve(eq6, std_e_initial_guess, full_output = False)

    #print("sol: ", std_e_solution)

    #parameters = LWE.Parameters(n, 2 ** logq, ND.UniformMod(2), ND.DiscreteGaussian(std_e_solution[0]))
    #estimator_output = LWE.primal_bdd(parameters)
    #print("Verification: ", estimator_output)

    return std_e_solution[0]

def numerical_std_e_usvp(l, n, logq, std_s):
    lnq = logq * log(2)

    # find std_e and beta numerically
    beta_initial_guess = (l - 16.4) / 0.292
    zeta_initial_guess = 10
    std_e_initial_guess = zeta_initial_guess * std_s
    
    nom = lambda std_e, beta : 2 * n * (lnq - log(max(1,round(std_e/std_s)))) * log(beta/const)
    denom = lambda std_e, beta : lnq + log(sqrt(beta)/(const*std_e))
    eq11 = lambda std_e, beta : beta - nom(std_e, beta) / (denom(std_e, beta)**2)
    
    d_optimal = lambda std_e, beta : sqrt(2 * n * (lnq - log(max(1,round(std_e / std_s)))) * beta / log(beta / const))
    eq12 = lambda std_e, beta : l - (0.292 * beta + log2(8 * d_optimal(std_e, beta)) + 16.4)

    def system_usvp_std_e(x): # x[0] = std_e, x[1] = beta
        f1 = eq11(x[0], x[1])
        f2 = eq12(x[0], x[1])
        return f1, f2
    
    std_e_solution, beta_solution = fsolve(system_usvp_std_e, [std_e_initial_guess, beta_initial_guess], full_output = False)
    return std_e_solution