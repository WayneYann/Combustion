import t
from math import *
import matplotlib.pyplot as plt

a = 1.0
d = -100.5
r = 10.0

misdc_iters = [1, 5, 10]

dts = [0.25]
for j in range(10):
    dts.append(dts[j]/2)

E = []
for dl in [True, False]:
    for n in misdc_iters:
        e = []
        for dt in dts:
            e.append(t.solve_it(a, d, r, dt, max_iter=n, do_linear=dl))
        E.append(e)
        
        for i in range(len(e)-1):
            print 'iters: ', n, 'order: ', log(e[i+1]/e[i])/log(dts[i+1]/dts[i])
        
        plt.loglog(dts, e, label=str(n)+' linear? '+str(dl), marker='.')

plt.xlabel('dt')
plt.ylabel('L^1 error')
plt.legend()
plt.legend(loc=2)
plt.show()
