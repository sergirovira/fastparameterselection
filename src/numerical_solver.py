from numpy import pi, exp, log, log2, sqrt, divide
from scipy.optimize import fsolve

import numpy as np


const = 2 * pi * exp(1)
ln2 = log(2)
e = exp(1)

def _delta(beta):
    return (beta / (2 * pi * e) * (pi * beta) ** (1 / beta)) ** (1 / (2 * (beta - 1)))


def numerical_lambda_bdd(n, logq, std_s, std_e):
    lnq = logq * ln2
    zeta = max(1,round(std_e/std_s))

    # find beta numerically
    beta_initial_guess = n / 4

    nom = lambda beta : 2 * n * lnq * log(beta/const)
    denom = lambda beta : log(beta/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
    eq6 = lambda beta : beta - nom(beta) / (denom(beta)**2)

    beta_solution = fsolve(eq6, beta_initial_guess, full_output = False)

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


def numerical_n_bdd(l, logq, std_s, std_e):

    eta_initial_guess = (l - 16.4) / 0.292
    eta_eq = lambda eta : l - (0.292*eta+log2(8*eta)+16.4)
    eta_solution = fsolve(eta_eq, eta_initial_guess, full_output = False)
    eta = eta_solution[0]

    lnq = logq * ln2
    d_optimal = lambda n, beta : sqrt(n * lnq / log(_delta(beta)))
    eq8 = lambda n, beta : l - (0.292 * beta + log2(8 * d_optimal(n, beta)) + 16.4)

    ln_std_e = log(std_e)
    zeta = max(0, ln_std_e - log(std_s))
    d = lambda n, beta : max(d_optimal(n, beta), n)
    eq_for_n_and_beta = lambda n, beta: eta - d(n, beta) + 1/log(_delta(beta)) * (lnq - ln_std_e - 0.5*log(const) - n / d(n, beta) * (lnq - zeta ) )

    def system_bdd_eta_n(x): # x[0] = n, x[1] = beta
        f1 = eq_for_n_and_beta(x[0], x[1])
        f2 = eq8(x[0], x[1])
        return f1, f2

    n_initial_guess = 100
    solutions_n_and_beta = fsolve(system_bdd_eta_n, [n_initial_guess, eta], full_output = True)
    n_solution = solutions_n_and_beta[0][0]
    if not solutions_n_and_beta[2]==1:
        print(solutions_n_and_beta[3])
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
    lnq_initial_guess = n/100

    zeta = std_e/std_s

    beta = (l-16.4)/0.292

    nom = lambda lnq : 2 * n * lnq * log(beta/const)
    denom = lambda lnq : log(beta/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
    eq = lambda lnq : beta - nom(lnq) / (denom(lnq)**2)

    lnq_solution = fsolve(eq, lnq_initial_guess, full_output = False)
    
    return divide(lnq_solution[0], ln2)

def numerical_logq_usvp(l, n, std_s, std_e):
    lnq_initial_guess = n/100

    zeta = std_e/std_s

    beta = (l-16.4)/0.292

    nom = lambda lnq : 2 * n * (lnq - log(zeta)) * log(beta/const)
    denom = lambda lnq : lnq + log(sqrt(beta)/(const*std_e))
    eq = lambda lnq : beta - nom(lnq) / (denom(lnq)**2)

    lnq_solution = fsolve(eq, lnq_initial_guess, full_output = False)

    return divide(lnq_solution[0], ln2)

# def numerical_logq_bdd(l, n, std_s, std_e):
#     eta_initial_guess = (l - 16.4) / 0.292
#     eta_eq = lambda eta : l - (0.292*eta+log2(8*eta)+16.4)
#     eta_solution = fsolve(eta_eq, eta_initial_guess, full_output = False)
#     eta = eta_solution[0]

#     d_optimal = lambda lnq, beta : sqrt(n * lnq / log(_delta(beta)))
#     eq8 = lambda lnq, beta : l - (0.292 * beta + log2(8 * d_optimal(lnq, beta)) + 16.4)

#     ln_std_e = log(std_e)
#     zeta = max(0, ln_std_e - log(std_s))
#     d = lambda lnq, beta : max(d_optimal(lnq, beta), n)
#     eq_for_lnq_and_beta = lambda lnq, beta: eta - d(lnq, beta) + 1/log(_delta(beta)) * (lnq - ln_std_e - 0.5*log(const) - n / d(lnq, beta) * (lnq - zeta ) )

#     def system_bdd_eta_lnq(x): # x[0] = lnq, x[1] = beta
#         f1 = eq_for_lnq_and_beta(x[0], x[1])
#         f2 = eq8(x[0], x[1])
#         return f1, f2

#     lnq_initial_guess = 3 * log(n)
#     solutions_lnq_and_beta = fsolve(system_bdd_eta_lnq, [lnq_initial_guess, eta], full_output = True)
#     lnq_solution = solutions_lnq_and_beta[0][0]

#     if not solutions_lnq_and_beta[2]==1:
#         print(solutions_lnq_and_beta[3])
#     logq_solution = lnq_solution / ln2
#     return logq_solution


# def numerical_logq_usvp(l, n, std_s, std_e):
#     zeta = max(1,round(std_e/std_s))

#     # find lnq and beta numerically
#     lnq_initial_guess = 3 * log(n)
#     beta_initial_guess = (l - 16.4) / 0.292
    
#     nom = lambda lnq, beta : 2 * n * lnq * log(beta/const)
#     denom = lambda lnq, beta : log(beta/const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
#     eq11 = lambda lnq, beta : beta - nom(lnq, beta) / (denom(lnq, beta)**2)

#     d_optimal = lambda lnq, beta : sqrt(2 * n * (lnq - log(zeta)) * beta / log(beta / const))
#     eq12 = lambda lnq, beta : l - (0.292 * beta + log2(8 * d_optimal(lnq, beta)) + 16.4)

#     def system_usvp_lnq(x): # x[0] = lnq, x[1] = beta
#         f1 = eq11(x[0], x[1])
#         f2 = eq12(x[0], x[1])
#         return f1, f2
    
#     lnq_solution, beta_solution = fsolve(system_usvp_lnq, [lnq_initial_guess, beta_initial_guess], full_output = False)
#     logq_solution = lnq_solution / ln2
#     return logq_solution


def numerical_std_e_bdd(l, n, logq, std_s):
    eta_initial_guess = (l - 16.4) / 0.292
    eta_eq = lambda eta : l - (0.292*eta+log2(8*eta)+16.4)
    eta_solution = fsolve(eta_eq, eta_initial_guess, full_output = False)

    lnq = logq * ln2
    d_optimal = lambda beta : sqrt(n * lnq / log(_delta(beta)))
    eq8 = lambda beta : l - (0.292 * beta + log2(8 * d_optimal(beta)) + 16.4)

    beta_solution = fsolve(eq8, eta_solution, full_output = False)
    d = max(d_optimal(beta_solution[0]), n)
    eta = eta_solution[0]
    beta = beta_solution[0]

    std_e_initial_guess = 5.502177429822036#0.001
    zeta = lambda ln_std_e : max(0, ln_std_e-log(std_s))
    eq_for_sigmae = lambda ln_std_e: eta - d + 1/log(_delta(beta)) * (lnq - ln_std_e - 0.5*log(const) -n/d*(lnq - zeta(ln_std_e) ) )
    std_e_solution = fsolve(eq_for_sigmae, std_e_initial_guess, full_output = True)

    #sigma_e_explicit = (d - eta - 1/log(_delta(beta))*(lnq - 0.5*log(const) - n/d*(lnq+log(std_s))))/(1/log(_delta(beta))*(n/d-1))
    #sigma_e_explicit2 = (eta + 1/log(_delta(beta))*(lnq - 0.5*log(const) - n*lnq/d) - d)/(1/log(_delta(beta)))

    if not std_e_solution[2]==1:
        print(std_e_solution[3])

    return exp(std_e_solution[0][0])


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