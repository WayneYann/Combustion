import pde3
import pde
import numpy as np
import matplotlib.pyplot as plt
import pdb

from math import *
from scipy.interpolate import interp1d, griddata
from scipy.integrate import cumtrapz, quad
import itertools
from textwrap import wrap

marker = itertools.cycle(('x', '.', '*', '+')) 

a = -1.0
eps = 10.0
r = -70.0
#eps = 0.0
#r = 0.0

Nx = 100
dt = (20.0/Nx)/5.0

FT = dt*1

N_A = 2

hs = []

def solve(k, i=2):
   global hs
   hs = []
   solns = []
   
   for j in range(k):
      dtl = dt/(2**j)
      (x, y) = pde3.solve_pde(a, eps, r, Nx=Nx*(2**j), dt=dtl,
                             FT=FT, max_iter=i)
      print
      solns.append((x,y))
      
      hs.append(20.0/(Nx*(2**j)))
   
   return solns

def solve_old(k, i=2):
   global hs
   hs = []
   solns = []
   
   for j in range(k):
      dtl = dt/(2**j)
      (x, y) = pde.solve_pde(a, eps, r, Nx=Nx*(2**j), dt=dtl,
                             FT=FT, max_iter=i, N_A=4, do_linear=False, do_modified=True)
      print
      solns.append((x,y))
      
      hs.append(20.0/(Nx*(2**j)))
   
   return solns

def do_advection_only3(k):
   solns = []
   e = []
   
   for j in range(k):
      dtl = dt/(2**j)
      (x, y) = pde3.solve_pde(a, 0.0, 0.0, Nx=Nx*2**(j), dt=dtl,
                             FT=FT, max_iter=1)
      
      exacty = 0.5*(np.tanh(10-2*(a*FT + x))+1)
      
      solns.append((x,y))
      e.append(cumtrapz(np.abs(y - exacty), x)[-1])
   return (solns, e)

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
      o.append(log(e[j]/e[j+1])/log(2))
   return o

def do_it(i):
   solns = solve_old(5, i=i)
   e = errors(solns)
   o = orders(e)
   print e
   print o
   plt.loglog(hs[:-1], e, marker=marker.next(), label=str(i))   

fig = plt.figure(1, figsize=(10,6))

do_it(1)
do_it(2)
do_it(3)
do_it(4)

fig.subplots_adjust(right=0.6)
fig.subplots_adjust(top=0.8)

plt.legend(loc=2, bbox_to_anchor=(1.0,0.8), title='MISDC iterations')
plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$L^1$ error')
#title = "\n".join(wrap("$L^1$ error for advection, diffusion, reaction equation"
#                      +" using MISDC on 3 Gauss-Lobatto nodes", 40))
#plt.title(title)
#plt.savefig('orders.pdf')

plt.show(block=True)

#fig = plt.figure()
#for (x,y) in solns:
#   plt.plot(x,y)
#plt.ylim((-0.1, 1.1))
#plt.show(block=True)

