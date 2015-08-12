import numpy as np
import matplotlib.pyplot as plt
import pdb

import finalpde as pde

from math import *
from scipy.interpolate import interp1d, griddata
from scipy.integrate import cumtrapz, quad
import itertools
from textwrap import wrap

marker = itertools.cycle(('x', '.', '*', '+')) 

ratio = 2

a = -0.1
eps = 1.0
r = -10.0

Nx = 50
dt = (20.0/Nx)/2.0

FT = dt*1

N_A = 2

hs = []

def solve(k, i=2):
   global hs
   hs = []
   solns = []
   
   for j in range(k):
      dtl = dt/(ratio**j)
      
      (x, y) = pde.solve_pde(a, eps, r, Nx=Nx*(ratio**j), dt=dtl,
                             FT=FT, max_iter=i)
      print
      solns.append((x,y))
      
      hs.append(20.0/(Nx*(2**j)))
   return solns

def int_diff(x1, y1, x2, y2):
    y2i = griddata(x2, y2, x1)
    
    err = np.abs(y1 - y2i)
    return cumtrapz(err, x1)[-1]

def errors(solns):
   e = []
   for j in range(len(solns)-1):
      e.append(int_diff(solns[j][0],solns[j][1],solns[j+1][0],solns[j+1][1]))
   return e

def orders(e):
   o = []
   for j in range(len(e)-1):
      o.append(log(e[j]/e[j+1])/log(ratio))
   return o

def do_it(i):
   solns = solve(6, i=i)
   e = errors(solns)
   o = orders(e)
   print e
   print o
   plt.loglog(hs[:-1], e, marker=marker.next(), label=str(i))
   return e

fig = plt.figure(1, figsize=(10,6))

E = []

E.append(do_it(1))
E.append(do_it(2))
E.append(do_it(3))
E.append(do_it(4))

E.append(hs[:-1])

np.savetxt('error-output', np.array(E))

fig.subplots_adjust(right=0.6)
fig.subplots_adjust(top=0.8)

plt.legend(loc=2, bbox_to_anchor=(1.0,0.8), title='MISDC iterations')
plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$L^1$ error')
title = "\n".join(wrap("$L^1$ error for advection, diffusion, reaction equation"
                      +" using MISDC on 3 Gauss-Lobatto nodes", 40))
plt.title(title)
plt.savefig('orders.pdf')

plt.show(block=True)

#fig = plt.figure()
#for (x,y) in solns:
#   plt.plot(x,y)
#plt.ylim((-0.1, 1.1))
#plt.show(block=True)

