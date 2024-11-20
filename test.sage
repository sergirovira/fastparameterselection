h = 64
n = 2**15

wg = 17
ng = 1833
d = 62437
beta = 246
# for ng in range(14594-10050, 14594+350, 100):
#     print(ng, (log(n-ng,2)-log(ng,2)-log(n-h+wg-ng,2)+log(ng-wg,2)+log(1-wg/ng,2)+(0.5*wg/ng)).n())


def entropy(x):
    return -x*log(x,2) -(1-x)*log(1-x,2)

print((0.5*(ng*entropy(wg/ng)+wg)+log(d,2)).n(), 0.292*beta+16.4+3)

