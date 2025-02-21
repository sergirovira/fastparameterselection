from sage.all import RR, binomial
from numpy import pi, exp, log, log2, sqrt, multiply, divide, ceil
import scipy.special 
from scipy.special import binom, comb

import numpy as np
import random

from scipy.optimize import fsolve, least_squares, minimize, root, differential_evolution, basinhopping

verbose = 0

#[logq, h, beta, ng, d, wg, lambda]
all_params14 = [
    [250.0, 64, 311, 7315, 17222, 8, 152.10, 150.67429846599785, 151.42513857065637, -24.937919046466178],
    [300.0, 64, 243, 7105, 17690, 6, 135.44, 132.29737467394443, 135.26930465118946, -28.599403752415867],
    [400.0, 64, 197, 5496, 21313, 5, 113.43, 110.83224522566272, 113.16482721001513, -19.658127291147853],
    [450.0, 64, 222, 3181, 26328, 6, 106.69, 103.31953631152784, 106.54493493166798, -5.448729013560462],
    [500.0, 64, 160, 4346, 23821, 4, 97.83, 93.90837846904454, 97.73283261200493, -14.86965355319587],
    [550.0, 64, 175, 2416, 27915, 5, 91.67, 89.77235594718668, 91.22294924605522, -3.800166486227074],
    [250.0, 128, 417, 5299, 21819, 11, 190.64, 185.6603026797267, 190.58843271982445, -32.80788088550894],
    [300.0, 128, 386, 3750, 25173, 11, 165.33, 163.6394943048666, 164.8020302323806, -15.867193415880255],
    [350.0, 128, 363, 2236, 28353, 11, 146.87, 144.1366533603707, 146.63015726863057, -4.276461942199443],
    [400.0, 128, 305, 1984, 28852, 9, 130.17, 125.56136420444274, 130.10927372263944, -4.654373106319603],
    [450.0, 128, 262, 1714, 29399, 8, 117.08, 114.400206128055, 116.83494230376697, -3.9298118320351234],
    [500.0, 128, 231, 1339, 30159, 7, 106.55, 101.90599172999536, 106.48718847481027, -2.67015608865621],
    [550.0, 128, 202, 1153, 30526, 6, 97.89, 91.90064673754325, 97.86671343744594, -2.516120241328436],
    [250.0, 192, 532, 3231, 26322, 16, 208.80, 206.60382390664884, 208.43981749595886, -16.808587423948467],
    [300.0, 192, 469, 1920, 29035, 15, 178.02, 174.63125886789211, 177.87394513173112, -4.52831817677164],
    [350.0, 192, 388, 1578, 29718, 12, 153.90, 146.65411634187433, 153.89489940924022, -4.171863993063626],
    [400.0, 192, 330, 1217, 30443, 11, 135.22, 132.54052725221905, 134.97575101040903, -2.278473729472701],
    [450.0, 192, 280, 1073, 30722, 9, 120.75, 115.55317547577633, 120.7108407077843, -2.5664218061196418],
    [500.0, 192, 255, 391, 32118, 10, 108.88, 104.1129576141822, 108.83109443803613, -0.014336789249684585],
    [550.0, 192, 219, 406, 32078, 8, 99.33, 92.8764752275092, 99.31729657495731, -0.07732276095073985],
    [250.0, 256, 417, 5299, 21819, 11, 190.64, 185.6603026797267, 190.58843271982445, -32.80788088550894],
    [300.0, 256, 386, 3750, 25173, 11, 165.33, 163.6394943048666, 164.8020302323806, -15.867193415880255],
    [350.0, 256, 363, 2236, 28353, 11, 146.87, 144.1366533603707, 146.63015726863057, -4.276461942199443],
    [400.0, 256, 305, 1984, 28852, 9, 130.17, 125.56136420444274, 130.10927372263944, -4.654373106319603],
    [450.0, 256, 262, 1714, 29399, 8, 117.08, 114.400206128055, 116.83494230376697, -3.9298118320351234],
    [500.0, 256, 231, 1339, 30159, 7, 106.55, 101.90599172999536, 106.48718847481027, -2.67015608865621],
    [550.0, 256, 202, 1153, 30526, 6, 97.89, 91.90064673754325, 97.86671343744594, -2.516120241328436]]


all_params15 = [ [700, 64, 186, 14350, 34934, 5, 126], [750, 64, 228, 10330, 43999, 6, 121],
                 [800, 64, 218, 9493, 45832, 6, 116], [850, 64, 208, 8745, 47454, 6, 111],
                 [900, 64, 208, 7327, 50501, 6, 107], [950, 64, 182, 8014, 49024, 5, 103],
                 [1000, 64, 200, 5162, 55061, 6, 100],

                 [700, 128, 366, 4166, 57182, 11, 149], [750, 128, 339, 3677, 58191, 10, 140],
                 [800, 128, 284, 5189, 55043, 8, 132], [850, 128, 289, 3112, 59331, 9, 125],
                 [900, 128, 264, 3175, 59203, 8, 119], [950, 128, 255, 2207, 61164, 8, 113],
                 [1000, 128, 229, 2734, 60086, 7, 108],

                 [700, 192, 390, 2897, 59804, 12, 155], [750, 192, 357, 2642, 60322, 11, 145],
                 [800, 192, 328, 2416, 60773, 10, 138], [850, 192, 311, 1652, 62330, 10, 129],
                 [900, 192, 283, 1806, 62005, 9, 122], [950, 192, 272, 917, 63796, 10, 116],
                 [1000, 192, 253, 803, 64049, 9, 110]
                 ]

all_params16 = [[700, 64, 502, 29972, 66900, 12, 204.3086856519731], [750, 64, 496, 27694, 72150, 12, 199.37220382253997],
                [800, 64, 299, 37512, 49128, 7, 174.33098274157615], [850, 64, 292, 36075, 52534, 7, 168.79487539660414],
                [900, 64, 251, 37402, 49384, 6, 163.78559399692077], [950, 64, 303, 31484, 63348, 7, 160.90719271637388],

                [1000, 64, 247, 34312, 56697, 6, 154.64053954098628], [700, 128, 539, 28087, 71276, 13, 249.88276051600195],
                [750, 128, 501, 27433, 72768, 12, 238.6364372746093], [800, 128, 566, 21107, 86998, 14, 231.7410914426918],
                [850, 128, 492, 22833, 83154, 12, 219.74323550817536], [900, 128, 491, 20366, 88615, 12, 211.576871641367],
                [950, 128, 458, 20137, 89112, 11, 203.43093137911046], [1000, 128, 475, 16489, 97048, 12, 196.5746228989299],

                [700, 192, 693, 20518, 88335, 17, 286.296735151414], [750, 192, 644, 19810, 89879, 16, 271.35354500033776],
                [800, 192, 674, 15094, 100094, 17, 258.8319896335122], [850, 192, 598, 16406, 97262, 15, 245.46310719147453],
                [900, 192, 588, 14138, 102114, 15, 233.9091332247279], [950, 192, 575, 12149, 106327, 15, 223.40413201505788],
                [1000, 192, 541, 11737, 107182, 14, 213.62315678892801],

                [700, 256, 840, 13621, 103275, 22, 308.4712511484923], [750, 256, 762, 13788, 102902, 20, 290.042959469996],
                [800, 256, 718, 12699, 105197, 19, 273.9726320045856], [850, 256, 679, 11646, 107409, 18, 259.33443985905603],
                [900, 256, 707, 6747, 117592, 20, 248.70589123795554], [950, 256, 664, 6265, 118573, 19, 235.78921886091717],
                [1000, 256, 623, 5967, 119170, 18, 223.9543039378368],

                [700, 512, 1030, 5050, 121177, 32, 342.8540670705324], [750, 512, 947, 4672, 121897, 29, 318.9133585374426],
                [800, 512, 877, 4260, 122739, 27, 297.85611374146856], [850, 512, 813, 3997, 123255, 25, 279.268050680993],
                [900, 512, 760, 3539, 124163, 24, 262.81633359803743], [950, 512, 707, 3498, 124248, 22, 248.03356251829015],
                [1000, 512, 664, 3162, 124919, 21, 234.78705890453406]]

all_params17 = [[800, 128, 699, 81027, 84065, 15, 350.09968802677224], [850, 128, 696, 77483, 92421, 15, 339.0132890998411],
                [900, 128, 649, 77111, 93297, 14, 327.8948900547392], [950, 128, 609, 76624, 94446, 13, 318.0851735688354],
                [1000, 128, 688, 67530, 115975, 15, 311.64359617713654],

                [700, 192, 882, 78641, 89721, 19, 456.5103580338542], [750, 192, 1036, 66339, 118846, 23, 441.2132564899414],
                [800, 192, 830, 73100, 102840, 18, 421.8560209553864], [850, 192, 1027, 57837, 138655, 23, 412.354268400881],
                [900, 192, 862, 63281, 125989, 19, 392.2665644085177], [950, 192, 1009, 50151, 156190, 23, 387.4267214325475],
                [1000, 192, 854, 56030, 142783, 19, 367.54015963955896], [700, 256, 1203, 63227, 126183, 27, 518.5415474215844],

                [750, 256, 1040, 66150, 119306, 23, 496.875156416131], [800, 256, 1072, 59937, 133821, 24, 475.4616704949971],
                [850, 256, 1064, 55805, 143348, 24, 456.56651201034595], [900, 256, 1016, 54080, 147284, 23, 439.07310472202175],
                [950, 256, 1048, 47774, 161558, 24, 422.94916480509175], [1000, 256, 1034, 44243, 169433, 24, 407.94532844240103],

                [700, 512, 1848, 35889, 187954, 46, 662.928570615052], [750, 512, 1756, 33078, 194012, 44, 625.5100119944597],
                [800, 512, 1667, 30700, 199099, 42, 591.5489890228017], [850, 512, 1831, 16256, 229510, 51, 579.4601853420198],
                [900, 512, 1504, 26837, 207301, 38, 532.2872140012998], [950, 512, 1498, 21375, 218804, 39, 506.1205897513737],
                [1000, 512, 1394, 21859, 217769, 36, 481.9016547943525]]

const = 2 * pi * exp(1)
ln2 = log(2)
e = exp(1)

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
    return (beta / (2 * pi * e) * (pi * beta) ** (1 / beta)) ** (1 / (2 * (beta - 1)))


def _delta_approx(beta):
    return (beta / (2 * pi * e))** (1 / (2 * (beta - 1)))

def log_delta(beta):
    return( (1 / (2 * (beta - 1))) *log2(beta / (2 * pi * e)) + (1 / (2 * (beta - 1))) * (1 / beta)*log2(pi * beta)  )




def entropy(x):
    return -x*log2(x) -(1-x)*log2(1-x)

def approx_binom(n, k):
    return n*entropy(k/n)+0.5*log2(n/(8*k*(n-k)))


def approx_startpoint(n, logq, h):
    sigma_s = sqrt(h/n)
    sigma_e = 3.19
    #eq1 = lambda beta, d: d*log2(beta)/beta - (1-n/d)*logq - n/d*log2(sigma_e/sigma_s) + log2(sigma_e)

    lnq = logq * ln2
    zeta = round(sigma_e/sigma_s)

    # find beta numerically
    beta_initial_guess = n / 4

    nom = lambda beta : 2 * n * lnq * log(beta/const)
    denom = lambda beta : log(beta/const) + 2 * lnq - 2 * log(sigma_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta/const) / (2 * lnq * beta))
    eq6 = lambda beta : beta - nom(beta) / (denom(beta)**2)

    beta_solution = fsolve(eq6, beta_initial_guess, full_output = False)

    # compute d (from 'FHE Formulas.pdf')
    d_optimal = sqrt(2 * n * lnq * beta_solution[0] / log(beta_solution[0] / const))
    return [beta_solution[0], d_optimal]

def probability_enum(n, h, ng, wg):
    try:
        if verbose: print(f"🔎 probability_enum inputs: n={n}, h={h}, ng={ng}, wg={wg}")

        denom = comb(int(n), int(ng), exact=True)
        if denom == 0:
            if verbose: print("⚠️ ERROR: binomial(n, ng) is 0, returning -inf")
            return -100  # Use a large negative number instead of -inf

        prob = comb(int(n-h), int(ng-wg), exact=True) * comb(int(h), int(wg), exact=True) / denom
        return log2(prob) if prob > 0 else -100
    except Exception as e:
        if verbose: print(f"❌ EXCEPTION in probability_enum: {e}")
        return -100  # Prevent crash


def ss_enum(ng, w):
    if not np.isnan(ng) and not np.isnan(wg):
        ss = sum(comb(int(ng), int(i), exact=True) * 2**i for i in range(0, w))
    else:
        ss = float('-inf')

    return log2(float(ss)) if ss > 0 else float('-inf')

epsilon = 1e-10  

def safe_log2(x):
    return log2(max(epsilon, x))

# def approx_binom(n, k):
#     return 0.5 * (k * safe_log2(n/k) + k)

def _delta(beta_):
    return max(epsilon, beta_)

def numerical_lambda_hybrid(n, logq, sigma_e, h):
    lnq = logq * ln2
    sigma_s = sqrt(h/n)
    xi = sigma_e/ sigma_s
    rt_min = float('inf')
    sol_tolerance_min = 10
    for wg in range(2,20):
        eq1 = lambda ng_, beta_, d_: (approx_binom(ng_,wg)+wg)+log2(d_) - (0.292*beta_+16.4+3)+2  #0.5*(ng_*entropy(wg/ng_)+wg) #- 0.5*log2(8*n*wg/ng_*(1-wg/ng_)) #0.5*(log2(math.comb(ng_, wg),2)+wg)
        eq2 = lambda ng_, beta_, d_: d_ - ceil(sqrt(n * lnq / log(_delta(beta_)))) + ng_ - 1 #adding ceil over sqrt() makes the eq. precise
        eq3 = lambda ng_, beta_, d_: (-d_-1)*log2(_delta(beta_))+ ((d_-n+ng_-1)*logq+ (n-ng_)*log2(xi))/d_ - 4 #success probability of Babai=1, i.e. ||b_d*|| = sigma_e \approx 4

        def system(x):
            f1 = eq1(x[0], x[1], x[2])
            f2 = eq2(x[0], x[1], x[2])
            f3 = eq3(x[0], x[1], x[2])
            return f1, f2, f3

        #print('approx start:', approx_startpoint(n,param[0], param[1]))
        initial_guess = approx_startpoint(n,logq, h)
        res = fsolve(system, [n/16, initial_guess[0]-10, initial_guess[1]-200], maxfev = 2**21, full_output=False)
        rt = 0.292*res[1]+log2(8*res[2])+16.4 - probability_enum(n, h, res[0], wg)
        sol_tolerance = eq1(res[0], res[1], res[2])+eq2(res[0], res[1], res[2])+eq3(res[0], res[1], res[2])
        print(wg, eq1(res[0], res[1], res[2]),eq2(res[0], res[1], res[2]), eq3(res[0], res[1], res[2]), rt, sol_tolerance)

        if(rt<rt_min and abs(sol_tolerance)<1e6*sol_tolerance_min):
            rt_min = rt
            ng_min = res[0]
            beta_min =res[1]
            d_min = res[2]
            wg_min = wg
            sol_tolerance_min = abs(sol_tolerance)
            print(wg, sol_tolerance)


    print(param[0],param[1], param[2], param[3], param[4], ": ", wg_min, beta_min,ng_min, d_min)
    print("------------------------------")
    return rt_min

def numerical_logq_hybrid(n, l, sigma_e, h):
    sigma_s = sqrt(h/n)
    xi = sigma_e/ sigma_s
    rt_min = float('inf')
    ng_min, beta_min, d_min, logq_min, wg_min = None, None, None, None, None
    sol_tolerance_min = 10

    beta_initial_guess = (l - log2(8*4*n) -  16.4) / 0.292
    logq_initial_guess = 700
    d_initial_guess = ceil(sqrt(n * logq_initial_guess * ln2 / log(_delta(beta_initial_guess))))

    for wg in range(2,30):
        # EXACT EQUATIONS:
        # eq1a = lambda ng_, beta_, d_: l - (0.292*beta_+16.4+3+safe_log2(d_)) + probability_enum(n, h, ng_, wg)-2
        # eq1b = lambda ng_, beta_, d_: l - ss_enum(ng_,wg)-2*safe_log2(d_) + probability_enum(n, h, ng_, wg)-5
        # eq2 = lambda ng_, beta_, d_, logq: d_ - ceil(sqrt(n * logq * ln2 / log(_delta(beta_)))) + ng_ - 1
        # eq3 = lambda ng_, beta_, d_, logq: (-d_ - 1) * safe_log2(_delta(beta_)) + ((d_ - n + ng_ - 1) * logq + (n - ng_) * safe_log2(xi)) / d_ - 4
        #
        #NON-EXACT EQUATIONS: 
        #eq1a = lambda ng_, beta_, d_: l - (0.292 * beta_ + 16.4 + 3 + safe_log2(d_)) + safe_log2(binom(n - h, ng_ - wg) * binom(h, wg) / max(epsilon, binom(n, ng_))) - 2
        #eq1b = lambda ng_, beta_, d_: l - safe_log2(binom(ng_, wg - 1)) - wg + 1 - 2 * safe_log2(d_) + safe_log2(binom(n - h, ng_ - wg + 1) * binom(h, wg - 1) / max(epsilon, binom(n, ng_))) - 4
        #eq2 = lambda ng_, beta_, d_, logq: d_ - ceil(sqrt(n * logq * ln2 / log(_delta(beta_)))) + ng_ - 1
        #eq3 = lambda ng_, beta_, d_, logq: (-d_ - 1) * safe_log2(_delta(beta_)) + ((d_ - n + ng_ - 1) * logq + (n - ng_) * safe_log2(xi)) / d_ - 4

        #EVEN MORE NON-EXACT EQUATIONS:
        eq1a = lambda ng_, beta_, d_: l - (0.292*beta_+16.4+3+log2(d_)) + approx_binom(n-h,ng_-wg)+approx_binom(h,wg) - approx_binom(n, ng_) -2
        eq1b = lambda ng_, beta_, d_:  l - approx_binom(ng_, wg)- wg - 2*log2(d_) + approx_binom(n-h,ng_-wg)+approx_binom(h,wg) - approx_binom(n, ng_) - 5
        eq2 = lambda ng_, beta_, d_, logq: d_ - sqrt(n * logq * ln2 / log(_delta(beta_))) + ng_ - 1 #adding ceil over sqrt() makes the eq. precise
        eq3 = lambda ng_, beta_, d_, logq: (-d_-1)*log_delta(beta_)+ ((d_-n+ng_-1)*logq+ (n-ng_)*log2(xi))/d_ - 4

        def system(x):
            f1a = eq1a(x[0], x[1], x[2])
            f1b = eq1b(x[0], x[1], x[2])
            f2 = eq2(x[0], x[1], x[2], x[3])
            f3 = eq3(x[0], x[1], x[2], x[3])
            return f1a, f1b, f2, f3

        ##[logq, h, beta, ng, d, wg, lambda]
        #res = fsolve(system, [n/16, param[3], param[4], logq_initial_guess], maxfev = 2**21, full_output=False)
        initial_guess = [param[3]-10, param[2]-5, param[4]+20, logq_initial_guess]
        bounds = ([n/10, n/200, n, 500], [n/2, n/100, 4*n, 1500])
        print(initial_guess)
        # print(bounds)
        # assert (False)
        #Determines the relative step size for the finite difference approximation of the Jacobian. The actual step is computed as x * diff_step.
        step_size = [5,2,10,20]
        res = least_squares(lambda x: system(x), initial_guess, diff_step=step_size, bounds=bounds, max_nfev=2**21).x
        #res = root(system, [param[3], param[2], param[4], logq_initial_guess], method='hybr').x
        #res = basinhopping(lambda x: sum(abs(y) for y in system(x)), [param[3], param[2], param[4], logq_initial_guess]).x


        rt = 0.292*res[1]+log2(8*res[2])+16.4 - probability_enum(n, h, res[0], wg)
        sol_tolerance = eq1a(res[0], res[1], res[2])+eq1b(res[0], res[1], res[2])+eq2(res[0], res[1], res[2], res[3])+eq3(res[0], res[1], res[2], res[3])
        print(f"wg: {wg} | rt: {rt:.5f} | tol: {sol_tolerance:.5f} | "
                          f"Eq1a: {eq1a(res[0], res[1], res[2]):.5f} | Eq1b: {eq1b(res[0], res[1], res[2]):.5f} | "
                          f"Eq2: {eq2(res[0], res[1], res[2], res[3]):.5f} | Eq3: {eq3(res[0], res[1], res[2], res[3]):.5f}")
        print(f"Solved values: ng={res[0]:.5f}, beta={res[1]:.5f}, d={res[2]:.5f}, logq={res[3]:.5f}")


        if(abs(sol_tolerance)<1e-5 and rt<rt_min): #if(rt<rt_min and abs(sol_tolerance)<1e6*sol_tolerance_min):
            rt_min = rt
            ng_min = res[0]
            beta_min =res[1]
            d_min = res[2]
            logq_min = res[3]
            wg_min = wg
            sol_tolerance_min = abs(sol_tolerance)
            print(wg, sol_tolerance)


    print(param[0], param[1], param[2], param[3], param[4], ": ", wg_min, logq_min, beta_min,ng_min, d_min)
    #print("------------------------------")

    return logq_min


#[logq, h, beta, ng, d, wg, lambda]
n = 2**14
for param in all_params14:

    if param[5]==0: continue
    wg = param[5]
    h = param[1]
    l = param[6]
    sigma_e = 3.19
    sigma_s = sqrt(h/n)
    xi = sigma_e/ sigma_s

    print(param)
    res = numerical_logq_hybrid(n, param[6], 3.19, param[1])
    print("------------------------------")
    print(n, param[6], param[1], wg)


# EXACT and NON-EXACT EQUATIONS:
#n=14: correcting constants: [-2, -5]
#n=15: correcting constants: [-2, -5]
#n=16: correcting constants: [-2, -5]
#n=17: correcting constants: [-2, -5]









