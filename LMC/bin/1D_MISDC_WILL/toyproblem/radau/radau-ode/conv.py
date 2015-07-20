import t
from math import *
import matplotlib.pyplot as plt
import itertools
from textwrap import wrap
marker = itertools.cycle(('x', '.', '*', '+')) 

a = 0.0
d = -3.0
r = 1.0

def method_str():
    return 'Gauss-Lobatto'

misdc_iters = [1, 2, 3, 4]

#dts = [0.125]
dts = [0.0625]
for j in range(7):
    dts.append(dts[j]/2)

fig = plt.figure(1, figsize=(10,6))

E = []
print method_str()
for n in misdc_iters:
    e = []
    for dt in dts:
        e.append(t.solve_it(d, r, dt, max_iter=n))
    E.append(e)
    
    for i in range(len(e)-1):
        if e[i+1] > 0 and e[i] > 0:
            print 'iters: ', n, 'order: ', log(e[i+1]/e[i])/log(dts[i+1]/dts[i])
        else:
            print 'iters: ', n, 'order:  ---'
    
    plt.loglog(dts, e, label=str(n), marker=marker.next())

fig.subplots_adjust(right=0.6)
fig.subplots_adjust(top=0.8)

plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$L^1$ error')
plt.legend(loc=2, bbox_to_anchor=(1.0,0.8), title='MISDC iterations')

title = "\n".join(wrap("$L^1$ error for MISDC using correction ODE on 3 "+
                       "Gauss-Lobatto nodes with constant forcing", 40))
plt.title(title)

plt.savefig('orders.pdf')

plt.show()

print E


