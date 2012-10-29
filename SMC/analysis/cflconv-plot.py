
import pickle
import matplotlib.pylab as plt
import numpy as np

from pyboxlib.utils import to_xy

with open('comps.pkl', 'r') as f:
    comps = pickle.load(f)


###############################################################################
# relative errors

plt.figure()

# sdc variants
for sdc in [ 3, 5, 9 ]:
    x, y = to_xy([ (x['cflfac'], x['rel']) for x in comps['sdc'] 
                   if x['sdc_nnodes'] == sdc ])

    plt.loglog(x, y, label='SDC' + str(sdc), marker='o')

    

# rk variant
x, y = to_xy([ (x['cflfac'], x['rel']) 
               for x in comps['rk'] ])

plt.loglog(x, y, label='RK3', marker='o')

plt.xlabel('cfl factor')
plt.ylabel('relative error')
plt.legend(loc='best')
plt.savefig('cflconv-convergence.png')


###############################################################################
# runtime

plt.figure()

# sdc variants
for sdc in [ 3, 5, 9 ]:
    x, y = to_xy([ (x['run'], x['rel']) for x in comps['sdc'] 
                   if x['sdc_nnodes'] == sdc ])

    plt.loglog(x, y, label='SDC' + str(sdc), marker='o')

# rk variant
x, y = to_xy([ (x['run'], x['rel']) 
               for x in comps['rk'] ])

plt.loglog(x, y, label='RK3', marker='o')

plt.xlabel('runtime')
plt.ylabel('relative error')
plt.legend(loc='best')
plt.savefig('cflconv-runtime.png')


plt.show()
    



