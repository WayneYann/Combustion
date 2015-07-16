import t
from math import *
import matplotlib.pyplot as plt

a = 0.0
d = -3.0
r = 1.0

# d = -4.0
# r = 3.0

def method_str():
    return 'Radau'

#misdc_iters = [1, 2, 4, 8, 16, 1000]
misdc_iters = [2, 3]

dts = [0.0125]
for j in range(2):
    dts.append(dts[j]/2)

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
    
    plt.loglog(dts, e, label=str(n))

fig.subplots_adjust(right=0.65)

plt.xlabel('dt')
plt.ylabel('L^1 error')
plt.legend(loc=2, bbox_to_anchor=(1.0,0.8))
plt.show()
