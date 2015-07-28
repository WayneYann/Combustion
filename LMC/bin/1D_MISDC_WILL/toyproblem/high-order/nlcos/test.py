import t
from math import *
import matplotlib.pyplot as plt

a = 7.0
d = -3.0
r = 4.0

misdc_iters = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

dts = [0.125]
for j in range(5):
    dts.append(dts[j]/2)

E = []


for exact in [False]:
    for n in misdc_iters:
        e = []
        for dt in dts:
            e.append(t.solve_it(a, d, r, dt, max_iter=n))
        E.append(e)
        
        for i in range(len(e)-1):
            print 'iterations: ', n, '  order: ', log(e[i+1]/e[i])/log(dts[i+1]/dts[i])
        
        plt.loglog(dts, e, label=str(n), marker='.')

plt.xlabel('dt')
plt.ylabel('L^1 error')
plt.legend()
plt.legend(loc=2)
plt.show()
