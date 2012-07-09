"""SDC_S3D run script."""

import re
import subprocess

from collections import defaultdict
from runutils import autorun

#### config

run = True

probin = {
  'tfinal':  1.0e-6, 
  'dt':      2.5e-8,
  'cfl':    -1.0,
  'method':  2,
  }

sdc = {
  'qtype': '"Gauss-Lobatto"',
  'tol_residual': 1e-12,
  }

runs = [ ('1000', 1.0e-7, 'plt00010'), 
         ('0500', 5.0e-8, 'plt00021'), 
         ('0250', 2.5e-8, 'plt00041'), 
         ('0125', 1.25e-8, 'plt00081'),
         ]

#### run

if run:

    autorun('conv/sdc13_dt0125', dt=1.25e-8,
            nnodes=13, method=2, defaults=defaults)

    for name, dt, plt in runs:

        autorun('conv/rk3_dt%s' % name, dt=dt, 
                method=1, defaults=defaults)

        autorun('conv/sdc3_dt%s' % name, 
                probin, { 'dt': dt }, 
                sdc,    { 'nnodes': 3 })

        # autorun('conv/sdc5_dt%s' % name, 
        #         probin, { 'dt': dt }, 
        #         sdc,    { 'nnodes': 5, 'tol_residual': 1e-4 })

        autorun('conv/sdc5_dt%s' % name, 
                probin, { 'dt': dt }, 
                sdc,    { 'nnodes': 5 })

#### compute errors

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

#### plot!

import matplotlib.pylab as plt

dt = [ x[1] for x in runs ]

plt.figure()
plt.loglog(dt, rho_errors['rk3'], '-ok', label='RK3')
plt.loglog(dt, rho_errors['sdc3'], '-sb', label='SDC3')
plt.loglog(dt, rho_errors['sdc5'], '-^r', label='SDC5')
plt.legend(loc='best')
plt.xlabel('dt')
plt.ylabel('max abs error')
plt.title('density')

plt.figure()
plt.loglog(dt, energy_errors['rk3'], '-ok', label='RK3')
plt.loglog(dt, energy_errors['sdc3'], '-sb', label='SDC3')
plt.loglog(dt, energy_errors['sdc5'], '-^r', label='SDC5')
plt.xlabel('dt')
plt.ylabel('max abs error')
plt.title('energy')


e0 = energy_errors['rk3'][0]
x = [ dt[0], dt[-1] ]
y = [ e0, e0*(x[1]/x[0])**3 ]
plt.loglog(x, y, '--k', label='3rd order')

e0 = energy_errors['sdc3'][0]
x = [ dt[0], dt[-1] ]
y = [ e0, e0*(x[1]/x[0])**4 ]
plt.loglog(x, y, '--b', label='4th order')

e0 = energy_errors['sdc5'][0]
x = [ dt[0], dt[-1] ]
y = [ e0, e0*(x[1]/x[0])**8 ]
plt.loglog(x, y, '--r', label='8th order')


plt.legend(loc='best')

plt.show()
