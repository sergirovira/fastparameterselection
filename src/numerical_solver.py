from numpy import pi, exp, log, sqrt, multiply, divide
from scipy.optimize import fsolve, brentq, brenth, bisect, newton

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

def numerical_std_e_bdd2(n, logq, std_s):
    zeta_initial_guess = 10

    std_e_guess = zeta_initial_guess*std_s

    beta = numerical_beta_bdd(n, logq, std_s, std_e_guess)

    print("beta: ", beta)
    print("sec: ", 0.292*beta + 16.4)

    beta = 80.0/0.292 - 16.4

    lnq = multiply(logq, ln2)

    nom = lambda zeta : 2 * n * lnq * log(beta/const)
    denom = lambda zeta : log(beta/const) + 2 * lnq - 2 * log(std_e_guess) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
    eq = lambda zeta : beta - nom(zeta) / (denom(zeta)**2)

    zeta_solution = fsolve(eq, zeta_initial_guess, full_output = False)
    sigma_solution = std_s*zeta_solution[0]
    
    return sigma_solution

def numerical_std_e_bdd(n, logq, std_s):
    std_e_guess = 3.19

    beta = numerical_beta_bdd(n, logq, std_s, std_e_guess)

    print("beta: ", beta)
    print("sec: ", 0.292 * beta + 16.4)

    #beta = (80.0 - 16.4) / 0.292

    lnq = multiply(logq, ln2)

    def nom():
        return 2 * n * lnq * log(beta / const)

    def denom(std_e):
        #print("zeta: ", zeta)
        return (
            log(beta / const)
            + 2 * lnq
            - 2 * log(std_e_guess)
            - log(const)
            - 2 * (lnq - log(std_s) - log(std_e)) * sqrt(n * log(beta / const) / (2 * lnq * beta))
        )

    def eq(std_e):
        return beta - nom() / (denom(std_e) ** 2)

    std_e = fsolve(eq, std_e_guess, full_output=False)

    return std_e[0]

def numerical_std_e_usvp(n, logq, std_s):
    zeta_initial_guess = 10

    std_e_guess = zeta_initial_guess*std_s

    beta = numerical_beta_usvp(n, logq, std_s, std_e_guess)

    lnq = multiply(logq, ln2)

    nom = lambda zeta : 2 * n * (lnq - log(zeta)) * log(beta/const)
    denom = lambda zeta : lnq + log(sqrt(beta)/(const*std_e_guess))
    eq = lambda zeta : beta - nom(zeta) / (denom(zeta)**2)

    zeta_solution = fsolve(eq, zeta_initial_guess, full_output = False)

    sigma_solution = std_s*zeta_solution[0]
    
    return sigma_solution