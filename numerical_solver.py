from numpy import pi, exp, log, log2, sqrt, multiply, divide
from scipy.optimize import fsolve

import numpy as np

# import sys
# sys.path.append('../lattice-estimator')
# estimator_installed = 1

import matplotlib.pyplot as plt

# try:
#     from estimator import *
# except ImportError:
#     print("Warning: Failed to import lattice_estimator, some options will not work")
#     estimator_installed = 0


const = 2 * pi * exp(1)
ln2 = log(2)

def numerical_lambda_bdd(n, logq, std_s, std_e):
    lnq = logq * ln2
    zeta = max(1,round(std_e/std_s))

    # find beta numerically
    beta_initial_guess = n / 4

    nom = lambda beta : 2 * n * lnq * log(beta/const)
    denom = lambda beta : log(beta/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
    eq6 = lambda beta : beta - nom(beta) / (denom(beta)**2)

    beta_solution = fsolve(eq6, beta_initial_guess, full_output = False)
    
    # compute d (from 'FHE Formulas.pdf')
    d_optimal = sqrt(2 * n * lnq * beta_solution[0] / log(beta_solution[0] / const))

    # compute l
    l_solution = 0.292 * beta_solution[0] + log2(8*d_optimal) + 16.4

    return l_solution


def numerical_lambda_usvp(n, logq, std_s, std_e):
    lnq = logq * ln2
    zeta = max(1,round(std_e/std_s))

    # find beta numerically
    beta_initial_guess = n / 4

    nom = lambda beta : 2 * n * (lnq - log(zeta)) * log(beta/const)
    denom = lambda beta : lnq + log(sqrt(beta)/(const*std_e))
    eq11 = lambda beta : beta - nom(beta) / (denom(beta)**2)

    beta_solution = fsolve(eq11, beta_initial_guess, full_output = False)
    
    # compute d (as sobstitute in eq12)
    d_optimal = sqrt(2 * n * (lnq - log(zeta)) * beta_solution[0] / log(beta_solution[0] / const))

    # compute l
    l_solution = 0.292 * beta_solution[0] + log2(8 * d_optimal) + 16.4

    return l_solution


def numerical_n_bdd(l, logq,  std_s, std_e):
    lnq = logq * ln2
    zeta = max(1,round(std_e/std_s))

    # find n and beta numerically
    n_initial_guess = 100
    beta_initial_guess = (l - 16.4) / 0.292

    nom = lambda n, beta : 2 * n * lnq * log(beta/const)
    denom = lambda n, beta : log(beta/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
    eq6 = lambda n, beta : beta - nom(n, beta) / (denom(n, beta)**2)

    d_optimal = lambda n, beta : sqrt(2 * n * lnq * beta / log(beta / const))
    eq8 = lambda n, beta : l - (0.292 * beta + log2(8 * d_optimal(n, beta)) + 16.4)

    def system_bdd_l(x): # x[0] = n, x[1] = beta
        f1 = eq6(x[0], x[1])
        f2 = eq8(x[0], x[1])
        return f1, f2 

    n_solution, beta_solution = fsolve(system_bdd_l, [n_initial_guess, beta_initial_guess], full_output = False)
    return n_solution


def numerical_n_usvp(l, logq, std_s, std_e):
    lnq = logq * ln2
    zeta = max(1,round(std_e/std_s))

    # find n and beta numerically
    n_initial_guess = 100
    beta_initial_guess = (l - 16.4) / 0.292

    nom = lambda n, beta : 2 * n * (lnq - log(zeta)) * log(beta/const)
    denom = lambda beta : lnq + log(sqrt(beta)/(const*std_e))
    eq11 = lambda n, beta : beta - nom(n, beta) / (denom(beta)**2)

    d_optimal = lambda n, beta : sqrt(2 * n * lnq * beta / log(beta / const))
    eq12 = lambda n, beta : l - (0.292 * beta + log2(8 * d_optimal(n, beta)) + 16.4)

    def system_usvp_l(x): # x[0] = n, x[1] = beta
        f1 = eq11(x[0], x[1])
        f2 = eq12(x[0], x[1])
        return f1, f2 

    n_solution, beta_solution = fsolve(system_usvp_l, [n_initial_guess, beta_initial_guess], full_output = False)
    return n_solution


def numerical_logq_bdd(l, n, std_s, std_e):
    zeta = max(1,round(std_e/std_s))

    # find lnq and beta numerically
    lnq_initial_guess = 3 * log(n)
    beta_initial_guess = (l - 16.4) / 0.292
    
    nom = lambda lnq, beta : 2 * n * lnq * log(beta/const)
    denom = lambda lnq, beta : log(beta/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
    eq6 = lambda lnq, beta : beta - nom(lnq, beta) / (denom(lnq, beta)**2)

    d_optimal = lambda lnq, beta : sqrt(2 * n * lnq * beta / log(beta / const))
    eq8 = lambda lnq, beta : l - (0.292 * beta + log2(8 * d_optimal(lnq, beta)) + 16.4)

    def system_bdd_lnq(x): # x[0] = lnq, x[1] = beta
        f1 = eq6(x[0], x[1])
        f2 = eq8(x[0], x[1])
        return f1, f2
    
    lnq_solution, beta_solution = fsolve(system_bdd_lnq, [lnq_initial_guess, beta_initial_guess], full_output = False)
    logq_solution = lnq_solution / ln2
    return logq_solution


def numerical_logq_usvp(l, n, std_s, std_e):
    zeta = max(1,round(std_e/std_s))

    # find lnq and beta numerically
    lnq_initial_guess = 3 * log(n)
    beta_initial_guess = (l - 16.4) / 0.292
    
    nom = lambda lnq, beta : 2 * n * lnq * log(beta/const)
    denom = lambda lnq, beta : log(beta/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
    eq11 = lambda lnq, beta : beta - nom(lnq, beta) / (denom(lnq, beta)**2)

    d_optimal = lambda lnq, beta : sqrt(2 * n * (lnq - log(zeta)) * beta / log(beta / const))
    eq12 = lambda lnq, beta : l - (0.292 * beta + log2(8 * d_optimal(lnq, beta)) + 16.4)

    def system_usvp_lnq(x): # x[0] = lnq, x[1] = beta
        f1 = eq11(x[0], x[1])
        f2 = eq12(x[0], x[1])
        return f1, f2
    
    lnq_solution, beta_solution = fsolve(system_usvp_lnq, [lnq_initial_guess, beta_initial_guess], full_output = False)
    logq_solution = lnq_solution / ln2
    return logq_solution


def numerical_std_e_bdd(l, n, logq, std_s):
    lnq = logq * ln2

    # find beta numerically
    beta_initial_guess = (l - 16.4) / 0.292
    d_optimal = lambda beta : sqrt(2 * n * lnq * beta / log(beta / const))
    eq8 = lambda beta : l - (0.292 * beta + log2(8 * d_optimal(beta)) + 16.4)

    beta_solution = fsolve(eq8, beta_initial_guess, full_output = False)
   
    # find std_e numerically
    std_e_initial_guess = 3.19

    zeta = lambda std_e : max(1, np.round(std_e/std_s))
    nom = 2 * n * lnq * log(beta_solution/const)
    denom = lambda std_e : log(beta_solution/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta(std_e))) * sqrt(n * log(beta_solution/const) / (2 * lnq * beta_solution))
    eq6 = lambda std_e : beta_solution - nom / (denom(std_e)**2)

    std_e_solution = fsolve(eq6, std_e_initial_guess, full_output = False)

    # d_lwe = LWE.primal_bdd(LWE.Parameters(n, 2 ** logq, ND.UniformMod(2), ND.DiscreteGaussian(std_e_solution[0])))["d"]
    # print("d lwe", d_lwe)
    # print("d eq6", d_optimal(beta_solution))
    # print(" ")

    return std_e_solution[0]

def numerical_std_e_bdd_eq5(l, n, logq, std_s): #how to the two functions compare? keep only one function.
    lnq = logq * ln2

    # find std_e and beta numerically
    beta_initial_guess = (l - 16.4) / 0.292
    std_e_initial_guess = 3.19

    zeta = lambda std_e : max(1, np.round(std_e/std_s))

    d_optimal = lambda std_e, beta : sqrt(2 * n * (lnq - log(zeta(std_e))) * beta / log(beta / const))

    nom = lambda std_e, beta : d_optimal(std_e, beta) * log(beta/const)
    denom = lambda std_e, beta : log(beta/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta(std_e))) * n / d_optimal(std_e, beta)
    eq5 = lambda std_e, beta : beta - nom(std_e, beta) / denom(std_e, beta)

    eq8 = lambda std_e, beta : l - (0.292 * beta + log2(8 * d_optimal(std_e, beta)) + 16.4)

    def system_bdd_std_e(x): # x[0] = std_e, x[1] = beta
        f1 = eq5(x[0], x[1])
        f2 = eq8(x[0], x[1])
        return f1, f2

    std_e_solution, beta_solution = fsolve(system_bdd_std_e, [std_e_initial_guess, beta_initial_guess], full_output = False)
    # d_lwe = LWE.primal_bdd(LWE.Parameters(n, 2 ** logq, ND.UniformMod(2), ND.DiscreteGaussian(std_e_solution)))["d"]
    # print("d lwe", d_lwe)
    # print("d eq5", d_optimal(std_e_solution, beta_solution))
    # print(" ")
    return std_e_solution


def numerical_std_e_usvp(l, n, logq, std_s):
    lnq = logq * ln2

    # find std_e and beta numerically
    beta_initial_guess = (l - 16.4) / 0.292
    std_e_initial_guess = 3.19

    zeta = lambda std_e : max(1, np.round(std_e/std_s))
    
    nom = lambda std_e, beta : 2 * n * (lnq - log(zeta(std_e))) * log(beta/const)
    denom = lambda std_e, beta : lnq + log(sqrt(beta)/(const*std_e))
    eq11 = lambda std_e, beta : beta - nom(std_e, beta) / (denom(std_e, beta)**2)
    
    d_optimal = lambda std_e, beta : sqrt(2 * n * (lnq - log(zeta(std_e))) * beta / log(beta / const))
    eq12 = lambda std_e, beta : l - (0.292 * beta + log2(8 * d_optimal(std_e, beta)) + 16.4)

    def system_usvp_std_e(x): # x[0] = std_e, x[1] = beta
        f1 = eq11(x[0], x[1])
        f2 = eq12(x[0], x[1])
        return f1, f2
    
    std_e_solution, beta_solution = fsolve(system_usvp_std_e, [std_e_initial_guess, beta_initial_guess], full_output = False)
    return std_e_solution