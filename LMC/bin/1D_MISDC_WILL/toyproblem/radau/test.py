import t
from math import *
import matplotlib.pyplot as plt

a = 0.0
d = -1000.0
r = -30.0

def method_str():
    return 'Lobatto'

misdc_iters = [1, 2, 4, 8, 16, 1000, 2000]

# it should be that if piecewise constant w/ODE is stable for given dt
# then three-node Gauss Lobatto is stable for dt/5
dts = [0.125/5.0]

#for j in range(5):
#    dts.append(dts[j]/2)

fig = plt.figure(1, figsize=(10,6))

E = []
print method_str()
for n in misdc_iters:
    e = []
    for dt in dts:
        er = t.solve_it(d, r, dt, max_iter=n)
        if er == 0:
            er = 1e-50
        e.append(er)
    E.append(e)
    
    for i in range(len(e)-1):
        print 'iters: ', n, 'order: ', log(e[i+1]/e[i])/log(dts[i+1]/dts[i])
#    
#    plt.loglog(dts, e, label=str(n))

#fig.subplots_adjust(right=0.65)

#plt.xlabel('dt')
#plt.ylabel('L^1 error')
#plt.legend(loc=2, bbox_to_anchor=(1.0,0.8))
#plt.show()
