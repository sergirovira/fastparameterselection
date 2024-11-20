from math import sqrt, log, log2
from numerical_solver import *
from sage.all import RR

from latticeestimator.estimator import *
from latticeestimator.estimator.nd import NoiseDistribution
from latticeestimator.estimator.lwe_parameters import LWEParameters

PI = 3.14159265358979
e = 2.71828182845905
const = 2*PI*e
ln2 = log(2)

def _delta(beta):
    small = (
        (2, 1.02190),
        (5, 1.01862),
        (10, 1.01616),
        (15, 1.01485),
        (20, 1.01420),
        (25, 1.01342),
        (28, 1.01331),
        (40, 1.01295),
    )

    if beta <= 2:
        return 1.0219
    elif beta < 40:
        for i in range(1, len(small)):
            if small[i][0] > beta:
                return RR(small[i - 1][1])
    elif beta == 40:
        return small[-1][1]
    else:
        return beta / (2 * PI * e) * (PI * beta) ** (1 / beta) ** (1 / (2 * (beta - 1)))

def compute_d(beta, n, lnq):
    num   = 2*n*lnq
    denom = log(_delta(beta))
    return min(2*n, sqrt(num/denom))


def eq5(n, lnq, std_e, std_s, beta, d_opt):

    #zeta = max(1,round(float(std_e/std_s)))
    zeta = RR(1)
    if std_s < std_e:
        zeta = std_e / std_s
    d = d_opt #compute_d(beta, n, logq)
    nom = d*log(beta/const)
    #denom = log(beta/const)+2*log(q/(std_e*sqrt(const)))-2*(n/d)*log(q/zeta)
    denom = log(beta/const)+2*(lnq - log(std_e*sqrt(const)))-2*(n/d)*(lnq - log(zeta))
    return int(d), float(nom/denom)


Xs = NoiseDistribution.UniformMod(2)
std_s = Xs.stddev

# res = 2**56.3 #113649.86013751852*2**44.3
# FHEParam = LWEParameters(
#     n=2**7,
#     q= 2**64,#q=2**64,
#     Xs=NoiseDistribution.UniformMod(2),
#     Xe=NoiseDistribution.DiscreteGaussian(stddev=res),
# )
# primal_bdd_cost = LWE.primal_bdd(FHEParam, red_cost_model=RC.BDGL16) #RC.BDGL16, RC.MATZOV
# print(primal_bdd_cost)
# d    = primal_bdd_cost['d']
# beta = primal_bdd_cost['beta']
# eta  = primal_bdd_cost['eta']
# ll = log2(primal_bdd_cost['rop'])
# check_d, check_beta = eq5(2**6,2**64*ln2, res ,std_s, beta, d)
#
# print(f" eta: {eta} beta:{beta} lambda:{ll} check_beta: {check_beta} ")


lambdas = [80, 128, 192, 256]
nlogs = [10, 11, 12, 13]
#lambdas = [80]
#nlogs = [10]

lambda_tfhe = [112]
ns = [291, 325, 351, 418]

logqs = [64]
for l in lambda_tfhe:
    for n in ns:#nlog in nlogs:
        for logq in logqs:
            #res = numerical_std_e_bdd_eq5(l, 2**nlog, logq, std_s)
            #res = numerical_std_e_bdd_ln(l, 2**nlog, logq, std_s)
            res = numerical_std_e_bdd_eta(l, n, logq, std_s)
            if res==-1:
                print("------------------------------")
                continue
            #print(f"lambda:{l} n:{2**nlog} q:{logq } std_e: {res}")
            #print(res, 10**(-5), res<10**(-5))
            if res<1e-19:
                print(f"{l, n, logq} too small st. deviation, continue...")
                print("------------------------------")
                continue

            FHEParam = LWEParameters(
                n= n,#2**nlog,
                q= 2**logq,#q=2**64,
                Xs=NoiseDistribution.UniformMod(2),
                Xe=NoiseDistribution.DiscreteGaussian(stddev=res),
            )

            try:
                primal_bdd_cost = LWE.primal_bdd(FHEParam, red_cost_model=RC.BDGL16) #RC.BDGL16, RC.MATZOV
                d    = primal_bdd_cost['d']
                beta = primal_bdd_cost['beta']
                eta  = primal_bdd_cost['eta']
                ll = log2(primal_bdd_cost['rop'])
            except:
                print("Error in the estimator")
                assert(False)
            lnq = logq * ln2
            #for std_es in range(2**8, 2**13, 500):
                #checking eq(5) for all given values
            #check_d, check_beta = eq5(2**nlog,lnq, res ,std_s, beta, d)
            #check_d, check_beta = eq5(2**nlog,lnq, res ,std_s, 262, 1890)
            print(f"lambda:{l} n:{n} logq:{logq } std_e: {res} eta: {eta} beta:{beta} lambda:{ll} d:{d} ")
            print("------------------------------")
            #assert(False)






