"""SDC_S3D run script."""

import re
import subprocess

from collections import defaultdict
from runutils import autorun


defaults = {
  'tfinal':  1.0e-6, 
  'dt':      2.5e-8,
  'cfl':    -1.0,
  'method':  2,
  'nnodes':  3,
  }


# runs

runs = [ ('1000', 1.0e-7, 'plt00010'), 
         ('0500', 5.0e-8, 'plt00021'), 
         ('0250', 2.5e-8, 'plt00041'), 
         ('0125', 1.25e-8, 'plt00081') ]

# autorun('conv/sdc13_dt0125', dt=1.25e-8,
#         nnodes=13, method=2, defaults=defaults)


# for name, dt, plt in runs:

#     autorun('conv/rk3_dt%s' % name, dt=dt, 
#             method=1, defaults=defaults)

#     autorun('conv/sdc3_dt%s' % name, dt=dt, 
#             nnodes=3, method=2, defaults=defaults)

#     autorun('conv/sdc5_dt%s' % name, dt=dt, 
#             nnodes=5, method=2, defaults=defaults)

# compute errors


rho_errors = defaultdict(list)
energy_errors = defaultdict(list)

for method in [ 'rk3', 'sdc3', 'sdc5' ]:
    for name, dt, plt in runs:
        p = subprocess.Popen(
            ['./fcompare.Linux.gfortran.mpi.omp.exe',
             '--infile1', 'conv/%s_dt%s/%s' % (method, name, plt),
             '--infile2', 'conv/sdc13_dt0125/plt00081'],
            stdout=subprocess.PIPE
            )

        out, err = p.communicate()

        m   = re.search('Variable_01\s+(\S+)', out)
        err = m.group(1)
        rho_errors[method].append(float(err))

        m   = re.search('Variable_05\s+(\S+)', out)
        err = m.group(1)
        energy_errors[method].append(float(err))

dt = [ x[1] for x in runs ]

print dt
print rho_errors['rk3']
print rho_errors['sdc3']
print rho_errors['sdc5']

print energy_errors['rk3']
print energy_errors['sdc3']
print energy_errors['sdc5']

# import matplotlib.pylab as plt
# plt.plot(dt, errors['rk3'], '-ok', label='RK3')
# plt.plot(dt, errors['sdc3'], '-sb', label='SDC3')
# plt.plot(dt, errors['sdc5'], '-^r', label='SDC5')
# plt.legend()
# plt.xlabel('dt')
# plt.ylabel('max abs error (density)')
# plt.show()




