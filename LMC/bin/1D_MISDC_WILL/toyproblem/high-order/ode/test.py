import t
from math import *
import matplotlib.pyplot as plt

a = 0.0
d = -100.0
r = 10.0

misdc_iters = [1, 10, 100, 500, 5000]

dts = [0.25]
#for j in range(6):
#    dts.append(dts[j]/2)

E = []


for exact in [False]:
    for n in misdc_iters:
        e = []
        for dt in dts:
            e.append(t.solve_it(a, d, r, dt, max_iter=n, exact=exact))
        E.append(e)
        
        for i in range(len(e)-1):
            print 'iters: ', n, 'order: ', log(e[i+1]/e[i])/log(dts[i+1]/dts[i])
        
        plt.loglog(dts, e, label=str(n) + ' exact: ' + str(exact), marker='.')

plt.xlabel('dt')
plt.ylabel('L^1 error')
plt.legend()
plt.legend(loc=2)
#plt.show()
