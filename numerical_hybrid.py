from sage.all import RR, binomial
from numpy import pi, exp, log, log2, sqrt, multiply, divide, ceil
from scipy.optimize import fsolve
import scipy.special

#[logq, h, beta, ng, d, wg, lambda]

all_params14 = [[700, 64, 116, 2253, 28240, 4, 76.15782878506808], [750, 64, 106, 1920, 28924, 4, 72.00766336007604], [800, 64, 108, 798, 31207, 4, 68.46461135817178], [850, 64, 104, 164, 32515, 5, 64.77966727705204], [900, 64, 93, 106, 32623, 5, 61.58602469289659], [950, 64, 91, 0, 33456, 0, 61.00200074471519], [1000, 64, 91, 0, 34325, 0, 61.038996729309524], [700, 128, 148, 265, 32324, 7, 77.77841686845393], [750, 128, 133, 125, 32616, 7, 73.27134882222275], [800, 128, 120, 0, 32890, 0, 69.44536145158476], [850, 128, 105, 88, 32670, 6, 65.10664522401244], [900, 128, 94, 43, 32769, 7, 61.99789614522741], [950, 128, 91, 0, 33456, 0, 61.00200074471519], [1000, 128, 91, 0, 34325, 0, 61.038996729309524], [700, 192, 153, 0, 32900, 0, 79.08179996359026], [750, 192, 133, 125, 32616, 7, 73.27134882222275], [800, 192, 119, 85, 32732, 7, 69.21360645769066], [850, 192, 106, 59, 32776, 7, 65.44942196868422], [900, 192, 94, 62, 32750, 6, 61.925730501613074], [950, 192, 91, 0, 33456, 0, 61.00200074471519], [1000, 192, 91, 0, 34325, 0, 61.038996729309524], [700, 256, 150, 151, 32563, 8, 78.35999644908343], [750, 256, 133, 125, 32616, 7, 73.27134882222275], [800, 256, 119, 59, 32758, 8, 69.25413802826124], [850, 256, 107, 0, 32913, 0, 65.65037080108321], [900, 256, 94, 62, 32750, 6, 61.925730501613074], [950, 256, 91, 0, 33456, 0, 61.00200074471519], [1000, 256, 91, 0, 34325, 0, 61.038996729309524], [700, 512, 152, 58, 32780, 11, 78.95141727596588], [750, 512, 134, 63, 32745, 9, 73.64342664307411], [800, 512, 119, 59, 32758, 8, 69.25413802826124], [850, 512, 106, 59, 32776, 7, 65.44942196868422], [900, 512, 95, 0, 32894, 0, 62.145546900238294], [950, 512, 91, 10, 33446, 10, 61.0016100106636], [1000, 512, 91, 10, 34315, 10, 61.038617969222045]]

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

all_params16 = [[700, 64, 502, 29972, 66900, 12, 204.3086856519731], [750, 64, 496, 27694, 72150, 12, 199.37220382253997], [800, 64, 299, 37512, 49128, 7, 174.33098274157615], [850, 64, 292, 36075, 52534, 7, 168.79487539660414], [900, 64, 251, 37402, 49384, 6, 163.78559399692077], [950, 64, 303, 31484, 63348, 7, 160.90719271637388], [1000, 64, 247, 34312, 56697, 6, 154.64053954098628], [700, 128, 539, 28087, 71276, 13, 249.88276051600195], [750, 128, 501, 27433, 72768, 12, 238.6364372746093], [800, 128, 566, 21107, 86998, 14, 231.7410914426918], [850, 128, 492, 22833, 83154, 12, 219.74323550817536], [900, 128, 491, 20366, 88615, 12, 211.576871641367], [950, 128, 458, 20137, 89112, 11, 203.43093137911046], [1000, 128, 475, 16489, 97048, 12, 196.5746228989299], [700, 192, 693, 20518, 88335, 17, 286.296735151414], [750, 192, 644, 19810, 89879, 16, 271.35354500033776], [800, 192, 674, 15094, 100094, 17, 258.8319896335122], [850, 192, 598, 16406, 97262, 15, 245.46310719147453], [900, 192, 588, 14138, 102114, 15, 233.9091332247279], [950, 192, 575, 12149, 106327, 15, 223.40413201505788], [1000, 192, 541, 11737, 107182, 14, 213.62315678892801], [700, 256, 840, 13621, 103275, 22, 308.4712511484923], [750, 256, 762, 13788, 102902, 20, 290.042959469996], [800, 256, 718, 12699, 105197, 19, 273.9726320045856], [850, 256, 679, 11646, 107409, 18, 259.33443985905603], [900, 256, 707, 6747, 117592, 20, 248.70589123795554], [950, 256, 664, 6265, 118573, 19, 235.78921886091717], [1000, 256, 623, 5967, 119170, 18, 223.9543039378368], [700, 512, 1030, 5050, 121177, 32, 342.8540670705324], [750, 512, 947, 4672, 121897, 29, 318.9133585374426], [800, 512, 877, 4260, 122739, 27, 297.85611374146856], [850, 512, 813, 3997, 123255, 25, 279.268050680993], [900, 512, 760, 3539, 124163, 24, 262.81633359803743], [950, 512, 707, 3498, 124248, 22, 248.03356251829015], [1000, 512, 664, 3162, 124919, 21, 234.78705890453406]]

all_params17 = [[800, 128, 699, 81027, 84065, 15, 350.09968802677224], [850, 128, 696, 77483, 92421, 15, 339.0132890998411], [900, 128, 649, 77111, 93297, 14, 327.8948900547392], [950, 128, 609, 76624, 94446, 13, 318.0851735688354], [1000, 128, 688, 67530, 115975, 15, 311.64359617713654], [700, 192, 882, 78641, 89721, 19, 456.5103580338542], [750, 192, 1036, 66339, 118846, 23, 441.2132564899414], [800, 192, 830, 73100, 102840, 18, 421.8560209553864], [850, 192, 1027, 57837, 138655, 23, 412.354268400881], [900, 192, 862, 63281, 125989, 19, 392.2665644085177], [950, 192, 1009, 50151, 156190, 23, 387.4267214325475], [1000, 192, 854, 56030, 142783, 19, 367.54015963955896], [700, 256, 1203, 63227, 126183, 27, 518.5415474215844], [750, 256, 1040, 66150, 119306, 23, 496.875156416131], [800, 256, 1072, 59937, 133821, 24, 475.4616704949971], [850, 256, 1064, 55805, 143348, 24, 456.56651201034595], [900, 256, 1016, 54080, 147284, 23, 439.07310472202175], [950, 256, 1048, 47774, 161558, 24, 422.94916480509175], [1000, 256, 1034, 44243, 169433, 24, 407.94532844240103], [700, 512, 1848, 35889, 187954, 46, 662.928570615052], [750, 512, 1756, 33078, 194012, 44, 625.5100119944597], [800, 512, 1667, 30700, 199099, 42, 591.5489890228017], [850, 512, 1831, 16256, 229510, 51, 579.4601853420198], [900, 512, 1504, 26837, 207301, 38, 532.2872140012998], [950, 512, 1498, 21375, 218804, 39, 506.1205897513737], [1000, 512, 1394, 21859, 217769, 36, 481.9016547943525]]

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

def probability_enum(n, h, ng, w):
    prob = 0
    ng = int(ng)
    for i in range(0,w):
        prob+=RR(binomial(n-h,ng-i)*binomial(h,i)/binomial(n, ng))
    return log2(prob)

def ss_enum(ng, w):
    ss = 0
    for i in range(0,w):
        ss+=RR(binomial(ng, i) * 2**i)
    return log2(ss)



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
    logq_initial_guess = 800
    d_initial_guess = ceil(sqrt(n * logq_initial_guess * ln2 / log(_delta(beta_initial_guess))))

    for wg in range(2,20):
        # EXACT EQUATIONS: ERROR: NotImplementedError: The Function binomial does not support numpy arrays as arguments
        # eq1a = lambda ng_, beta_, d_: l - (0.292*beta_+16.4+3+log2(d_)) + probability_enum(n, h, ng_, wg)-2 # log2(binomial(n-h,ng_-wg)*binomial(h,wg)/binomial(n, ng_) + binomial(n-h,ng_-wg+1)*binomial(h,wg-1)/binomial(n, ng_)) #approx_binom(n-h, ng_-wg)+approx_binom(h, wg)-approx_binom(n, ng_) - 2   #log2(binomial(n-h,ng_-wg)*binomial(h,wg)/binomial(n, ng_) + binomial(n-h,ng_-wg+1)*binomial(h,wg-1)/binomial(n, ng_)) #(approx_binom(ng_,wg)+wg)+log2(d_) - (0.292*beta_+16.4+3)+2  #0.5*(ng_*entropy(wg/ng_)+wg) #- 0.5*log2(8*n*wg/ng_*(1-wg/ng_)) #0.5*(log2(math.comb(ng_, wg),2)+wg)
        # eq1b = lambda ng_, beta_, d_:  l - ss_enum(ng_,wg)-2*log2(d_) + probability_enum(n, h, ng_, wg)-5 #approx_binom(n-h, ng_-wg)+approx_binom(h, wg)-approx_binom(n, ng_)
        # eq2 = lambda ng_, beta_, d_, logq: d_ - ceil(sqrt(n * logq * ln2 / log(_delta(beta_)))) + ng_ - 1 #adding ceil over sqrt() makes the eq. precise
        # eq3 = lambda ng_, beta_, d_, logq: (-d_-1)*log2(_delta(beta_))+ ((d_-n+ng_-1)*logq+ (n-ng_)*log2(xi))/d_ - 4 #success proba

        def system(x):
            f1a = eq1a(x[0], x[1], x[2])
            f1b = eq1b(x[0], x[1], x[2])
            f2 = eq2(x[0], x[1], x[2], x[3])
            f3 = eq3(x[0], x[1], x[2], x[3])
            return f1a, f1b, f2, f3

        #print('approx start:', approx_startpoint(n,param[0], param[1]))
        #initial_guess = approx_startpoint(n,logq, h)
        res = fsolve(system, [n/16, param[3], param[4], logq_initial_guess], maxfev = 2**21, full_output=False)
        rt = 0.292*res[1]+log2(8*res[2])+16.4 - probability_enum(n, h, res[0], wg)
        sol_tolerance = eq1a(res[0], res[1], res[2])+eq1b(res[0], res[1], res[2])+eq2(res[0], res[1], res[2], res[3])+eq3(res[0], res[1], res[2], res[3])
        print(wg, eq1a(res[0], res[1], res[2]),eq1b(res[0], res[1], res[2]), eq2(res[0], res[1], res[2], res[3]), eq3(res[0], res[1], res[2], res[3]), rt, sol_tolerance)

        if(abs(sol_tolerance)<1e-5 and rt<rt_min): #if(rt<rt_min and abs(sol_tolerance)<1e6*sol_tolerance_min):
            rt_min = rt
            ng_min = res[0]
            beta_min =res[1]
            d_min = res[2]
            logq_min = res[3]
            wg_min = wg
            sol_tolerance_min = abs(sol_tolerance)
            print(wg, sol_tolerance)


    print(param[0],param[1], param[2], param[3], param[4], ": ", wg_min, logq_min, beta_min,ng_min, d_min)
    #print("------------------------------")

    return 0 #TODO: return logq

#[logq, h, beta, ng, d, wg, lambda]
n = 2**16
for param in all_params16:

    if param[5]==0: continue
    wg = param[5]
    h = param[1]
    l = param[6]
    sigma_e = 3.19
    sigma_s = sqrt(h/n)
    xi = sigma_e/ sigma_s
    #res = numerical_lambda_hybrid(n, param[0], 3.19, param[1])

    # EXACT EQUATIONS:
    # eq1a = lambda ng_, beta_, d_: l - (0.292*beta_+16.4+3+log2(d_)) + probability_enum(n, h, ng_, wg)-2 # log2(binomial(n-h,ng_-wg)*binomial(h,wg)/binomial(n, ng_) + binomial(n-h,ng_-wg+1)*binomial(h,wg-1)/binomial(n, ng_)) #approx_binom(n-h, ng_-wg)+approx_binom(h, wg)-approx_binom(n, ng_) - 2   #log2(binomial(n-h,ng_-wg)*binomial(h,wg)/binomial(n, ng_) + binomial(n-h,ng_-wg+1)*binomial(h,wg-1)/binomial(n, ng_)) #(approx_binom(ng_,wg)+wg)+log2(d_) - (0.292*beta_+16.4+3)+2  #0.5*(ng_*entropy(wg/ng_)+wg) #- 0.5*log2(8*n*wg/ng_*(1-wg/ng_)) #0.5*(log2(math.comb(ng_, wg),2)+wg)
    # eq1b = lambda ng_, beta_, d_:  l - ss_enum(ng_,wg)-2*log2(d_) + probability_enum(n, h, ng_, wg)-5 #approx_binom(n-h, ng_-wg)+approx_binom(h, wg)-approx_binom(n, ng_)
    # eq2 = lambda ng_, beta_, d_, logq: d_ - ceil(sqrt(n * logq * ln2 / log(_delta(beta_)))) + ng_ - 1 #adding ceil over sqrt() makes the eq. precise
    # eq3 = lambda ng_, beta_, d_, logq: (-d_-1)*log2(_delta(beta_))+ ((d_-n+ng_-1)*logq+ (n-ng_)*log2(xi))/d_ - 4 #success probability of Babai=1, i.e. ||b_d*|| = sigma_e \approx 4
    #print(param[0],param[1], param[2], param[5], ": ", eq1a(param[3],param[2],param[4]), eq1b(param[3],param[2],param[4]), eq2(param[3],param[2],param[4], param[0]),eq3(param[3],param[2],param[4], param[0]))
    #----------------------

    #NON-EXACT EQUATIONS:
    eq1a = lambda ng_, beta_, d_: l - (0.292*beta_+16.4+3+log2(d_)) + log2(binomial(n-h,ng_-wg)*binomial(h,wg)/binomial(n, ng_)) -2
    eq1b = lambda ng_, beta_, d_:  l - log2(RR(binomial(ng_, wg)))- wg - 2*log2(d_) + log2(binomial(n-h,ng_-wg)*binomial(h,wg)/binomial(n, ng_)) - 5
    eq2 = lambda ng_, beta_, d_, logq: d_ - ceil(sqrt(n * logq * ln2 / log(_delta(beta_)))) + ng_ - 1 #adding ceil over sqrt() makes the eq. precise
    eq3 = lambda ng_, beta_, d_, logq: (-d_-1)*log2(_delta(beta_))+ ((d_-n+ng_-1)*logq+ (n-ng_)*log2(xi))/d_ - 4
    print(param[0],param[1], param[2], param[5], ": ", eq1a(param[3],param[2],param[4]), eq1b(param[3],param[2],param[4]), eq2(param[3],param[2],param[4], param[0]),eq3(param[3],param[2],param[4], param[0]))

    #res = numerical_logq_hybrid(n, param[6], 3.19, param[1])
    #print("------------------------------")

# EXACT EQUATIONS:
#n=14: correcting constants: [-2, -5]
#n=15: correcting constants: [-2, -5]
#n=16: correcting constants: [-2, -5]
#n=17: correcting constants: [-2, -5]









