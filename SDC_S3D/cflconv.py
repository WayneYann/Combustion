"""SDC_S3D convergence test (error vs CFL)."""

import glob
import re
import subprocess

from collections import defaultdict
from runutils import autorun

#### config

run = True

probin = {
  'tfinal':  1.0e-6, 
  }

sdc = {
  'qtype': '"Gauss-Lobatto"',
  'tol_residual': -1.0,
  'iters': 0,
  }

runs = [ ('050', 0.5),
         ('100', 1.0),
         ('150', 1.5),
         ('200', 2.0),
         ('250', 2.5),
         ]

#### run

if run:

    # autorun('cflconv/sdc13_cfl050', 
    #         probin, { 'cfl': 0.5, 'method': 2 },
    #         sdc,    { 'nnodes': 13, 'tol_residual': 1e-8 })

    for name, cfl in runs:

        # autorun('cflconv/rk3_cfl%s' % name,
        #         probin, { 'cfl': cfl, 'method': 1 },
        #         sdc,    {})

        # autorun('cflconv/sdc3_cfl%s' % name, 
        #         probin, { 'cfl': cfl, 'method': 2 },
        #         sdc,    { 'nnodes': 3 })

        # autorun('cflconv/sdc5_cfl%s' % name, 
        #         probin, { 'cfl': cfl, 'method': 2 },
        #         sdc,    { 'nnodes': 5, 'tol_residual': 1e-4, 'iters': 100 })

        autorun('cflconv/sdc9_cfl%s' % name, 
                probin, { 'cfl': cfl, 'method': 2 },
                sdc,    { 'nnodes': 9, 'tol_residual': 1e-4, 'iters': 100 })


#### compute errors

rho_errors = defaultdict(list)
energy_errors = defaultdict(list)

for method in [ 'rk3', 'sdc3', 'sdc5', 'sdc9' ]:
    for name, cfl in runs:

        # get last plot number
        plts = glob.glob('cflconv/%s_cfl%s/plt*' % (method, name))
        if plts:
            plt  = sorted(plts)[-1]
        else:
            continue

        # compute error
        p = subprocess.Popen(
            ['./fcompare.Linux.gfortran.mpi.omp.exe',
             '--infile1', plt,
             '--infile2', 'cflconv/sdc13_cfl050/plt00017'],
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

cfl = [ x[1] for x in runs ]

plt.figure()
plt.loglog(cfl, rho_errors['rk3'], '-ok', label='RK3')
plt.loglog(cfl, rho_errors['sdc3'], '-sb', label='SDC3')
plt.loglog(cfl, rho_errors['sdc5'], '-^r', label='SDC5')
plt.loglog(cfl, rho_errors['sdc9'], '-vm', label='SDC9')
plt.legend(loc='best')
plt.xlabel('cfl')
plt.ylabel('max abs error')
plt.title('density')

plt.figure()
plt.loglog(cfl, energy_errors['rk3'], '-ok', label='RK3')
plt.loglog(cfl, energy_errors['sdc3'], '-sb', label='SDC3')
plt.loglog(cfl, energy_errors['sdc5'], '-^r', label='SDC5')
plt.loglog(cfl, energy_errors['sdc9'], '-vm', label='SDC9')
plt.xlabel('cfl')
plt.ylabel('max abs error')
plt.title('energy')

e0 = energy_errors['rk3'][0]
x = [ cfl[0], cfl[-1] ]
y = [ e0, e0*(x[1]/x[0])**3 ]
plt.loglog(x, y, '--k', label='3rd order')

e0 = energy_errors['sdc3'][0]
x = [ cfl[0], cfl[-1] ]
y = [ e0, e0*(x[1]/x[0])**4 ]
plt.loglog(x, y, '--b', label='4th order')

e0 = energy_errors['sdc5'][0]
x = [ cfl[0], cfl[-1] ]
y = [ e0, e0*(x[1]/x[0])**8 ]
plt.loglog(x, y, '--r', label='8th order')

e0 = energy_errors['sdc9'][0]
x = [ cfl[0], cfl[-1] ]
y = [ e0, e0*(x[1]/x[0])**8 ]
plt.loglog(x, y, '--m', label='16th order')

plt.legend(loc='best')

plt.show()
