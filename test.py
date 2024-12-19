from scipy.optimize import fsolve
from numpy import pi, exp, log, log2, sqrt, multiply, divide, ceil, max
import scipy.special
import math

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


# n:15 logq:700 h:128 eta: 2 beta:299 lambda:136.93750293515308 zeta:8209 d:49099 red:136.4268748378904 svp:135.19128999384753 prob:0.000255857412915181 wt:17
# n:15 logq:750 h:128 eta: 2 beta:282 lambda:129.9599337180137 zeta:7425 d:50759 red:129.69504814750738 svp:127.38419157814215 prob:0.000900574382913918 wt:16
# n:15 logq:800 h:128 eta: 2 beta:263 lambda:123.68384342306358 zeta:6911 d:51827 red:123.43861684737557 svp:121.00638738252408 prob:0.00150195804601612 wt:15
# n:15 logq:850 h:128 eta: 2 beta:255 lambda:117.89515759071558 zeta:5814 d:54129 red:117.52662862028342 svp:115.7458969930963 prob:0.0185552242809837 wt:15
# n:15 logq:900 h:128 eta: 2 beta:247 lambda:112.99590841004758 zeta:4798 d:56253 red:112.51309578319183 svp:111.18199718742733 prob:0.117085330639929 wt:15
# n:15 logq:950 h:128 eta: 2 beta:231 lambda:108.0767507652688 zeta:4416 d:56983 red:107.77952697412651 svp:105.65154451688593 prob:0.123360163161690 wt:14
# n:15 logq:1000 h:128 eta: 2 beta:216 lambda:103.54478126305857 zeta:4153 d:57524 red:103.37133921958673 svp:100.40268813476786 prob:0.126838697663571 wt:13

#[logq, h, beta, ng, d, wg]
all_params = [ [700, 64, 187, 14594, 34770, 10], [750, 64, 202, 12310, 40004, 11],
               [800, 64, 198, 11182, 42517, 11], [850, 64, 181, 11026, 42838, 10],
               [900, 64, 183, 9568, 46041, 10], [950, 64, 169, 9396, 46397, 9],
               [1000, 64, 159, 9000, 47230, 9],

               [700, 128, 299, 8209, 49099, 17], [750, 128, 282, 7425, 50759, 16],
               [800, 128, 263, 6911, 51827, 15], [850, 128, 255, 5814, 54129, 15],
               [900, 128, 247, 4798, 56253, 15], [950, 128, 231, 4416, 56983, 14],
               [1000, 128, 216, 4153, 57524, 13],

               [700, 256, 380, 3899, 58244, 24], [750, 256, 356, 3133, 59771, 23],
               [800, 256, 326, 2951, 60109, 21], [850, 256, 302, 2684, 60676, 20],
               [900, 256, 292, 1696, 62777, 21], [950, 256, 269, 1582, 62898, 19],
               [1000, 256, 251, 1423, 63264, 18] ]
ln2 = log(2)
e = exp(1)
# h = 128
# n = 2**15
# logq = 700
# lnq = logq * ln2
#
# wg = 17
# ng = 8209
# d = 49099
# beta = 299

def entropy(x):
    return -x*log2(x) -(1-x)*log2(1-x)

#print((0.5*(ng*entropy(wg/ng)+wg)+log(d,2)).n(), 0.292*beta+16.4+3)
#print((0.5*(log(binomial(ng, wg),2)+wg) + log(d,2)).n(), 0.292*beta+16.4+3 )
#print((log((n-ng)/ng,2) - log( (n-h+wg-ng)/(ng-wg) ,2) + (1-0.5/0.292)*log(1-wg/ng,2)).n())
#print((log((n-ng)/ng,2)).n(), log( (n-h+wg-ng)/(ng-wg) ,2).n(), (1-0.5/0.292)*log(1-wg/ng,2).n())
#print((log((n-ng)/ng,2) - log( (n-h+wg-ng/(ng-wg)) ,2) + 0.5*log(1-wg/ng,2)).n())

# eq1 = lambda ng_, beta_, d_: 0.5*(ng_*entropy(wg/ng_)+wg)+log(d_,2) - (0.292*beta_+16.4+3)
# eq2 = lambda ng_, beta_: log((n-ng_)/ng_,2) - log( (n-h+wg-ng_)/(ng_-wg) ,2) + (1-0.5/0.292)*log(1-wg/ng_,2)
# eq3 = lambda ng_, beta_, d_: d_ - sqrt(n * lnq / log(_delta(beta_))) + ng_ - 1
#
# def system(x):
#     f2 = eq2(x[0], x[1])
#     f1 = eq1(x[0], x[1], x[2])
#     f3 = eq3(x[0], x[1], x[2])
#     return f1, f2, f3

#res = fsolve(system, [4000, 150, 10000], full_output=True)
#print(res)

def approx_binom(n, k):
    return n*entropy(k/n)+0.5*log2(n/(8*k*(n-k)))

def prob_sum(n,ng, w):
    s = 0
    for i in range(1, w+1):
        s+=approx_binom(n-h, ng-i)+approx_binom(h, i)-approx_binom(n, ng)
    return s

#[logq, h, beta, ng, d, wg]
for param in all_params:
    wg = param[5]
    h = param[1]
    lnq = param[0] * ln2
    n = 2**15
    eq1 = lambda ng_, beta_, d_: 0.5*(approx_binom(ng_,wg)+wg)+log2(d_) - (0.292*beta_+16.4+3) -1 #0.5*(ng_*entropy(wg/ng_)+wg) #- 0.5*log2(8*n*wg/ng_*(1-wg/ng_)) #0.5*(log2(math.comb(ng_, wg),2)+wg)
    eq2 = lambda ng_: log2((n-ng_)/ng_)+log2( (n-h+wg-ng_)/(ng_-wg)) -0.5*log2(1-wg/ng_)
    eq3 = lambda ng_, beta_, d_: d_ - ceil(sqrt(n * lnq / log(_delta(beta_)))) + ng_ - 1 #adding ceil over sqrt() makes the eq. precise

    prob = lambda ng_ : approx_binom(n-h, ng_-wg)+approx_binom(h, wg)-approx_binom(n, ng_) #(n-h)*entropy((ng_-wg)/(n-h)) + h*entropy(wg/h) - n*entropy(ng_/n)
    eq4 = lambda ng_, d_ : 0.5*(approx_binom(ng_,wg)+wg)+2*log2(d_) -1 -prob(ng_)
    eq5 = lambda ng_, beta_, d_ :  (0.292*beta_+16.4+3 + log2(d_)) - prob(ng_)
    print(param[0],param[1], param[2], param[3], param[4], ": ", eq1(param[3], param[2], param[4]),eq2(param[3]), eq3(param[3], param[2], param[4]) )
    print(param[0],param[1], param[2], param[3], param[4], ": ", eq4(param[3],param[4]), eq5(param[3], param[2], param[4]), prob(param[3]))

    def system(x):
        f1 = eq1(x[0], x[1], x[2])
        f2 = eq2(x[0])
        f3 = eq3(x[0], x[1], x[2])
        return f1, f2, f3

    res = fsolve(system, [400, param[3], 2*n], maxfev = 2**21, full_output=False)
    #print(eq1(res[0], res[1], res[2]),eq2(res[0]), eq3(res[0], res[1], res[2]))
    #print(res)
    print("--------------")