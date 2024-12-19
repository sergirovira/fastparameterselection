from latticeestimator.estimator import *
from latticeestimator.estimator.nd import NoiseDistribution
from latticeestimator.estimator.lwe_parameters import LWEParameters

from math import sqrt, log, log2
from numerical_hybrid import *

def test_numerics_lambda():
    nlogs = [15]
    hs = [64]
    logqs = [800]
    for nlog in nlogs:
        for h in hs:
            for logq in logqs:
                res = numerical_lambda_hybrid(2**nlog, logq,NoiseDistribution.SparseTernary(2**nlog, p = h/2, m=h/2).stddev, 3.19, h)
                #print(f"n:{nlog} logq:{logq } h:{h}, ng: {res[0]}, beta: {res[1]}")


def run_estimator():
    nlogs = [15]
    hs = [64, 128, 192, 256, 512]
    logqs = [700,750,800,850,900,950,1000]
    #logqs = [800]
    hs = [256]
    for nlog in nlogs:
        for h in hs:
            for logq in logqs:
                FHEParam = LWEParameters(
                    n=2**nlog,
                    q= 2**logq,
                    Xs=NoiseDistribution.SparseTernary(2**nlog, p = h/2, m=h/2),
                    Xe=NoiseDistribution.DiscreteGaussian(stddev=3.19),
                )
                try:
                    primal_hybrid_cost = LWE.primal_hybrid(FHEParam, red_cost_model=RC.BDGL16)
                    ll = log2(primal_hybrid_cost['rop'])
                    beta = primal_hybrid_cost['beta']
                    eta  = primal_hybrid_cost['eta']
                    zeta = primal_hybrid_cost['zeta']
                    d = primal_hybrid_cost['d']
                    svp_cost = log2(primal_hybrid_cost["svp"])
                    red_cost = log2(primal_hybrid_cost["red"])
                    prob =  primal_hybrid_cost["prob"]
                    wt = primal_hybrid_cost["wt"]
                    print(f"n:{nlog} logq:{logq } h:{h} eta: {eta} beta:{beta} lambda:{ll} zeta:{zeta} d:{d} red:{red_cost} svp:{svp_cost} prob:{prob} wt:{wt} ")
                    print("------------------------------")
                except Exception as e:
                     print('error in the estimator:', e)
                     print("------------------------------")


if __name__ == "__main__":
    run_estimator()
