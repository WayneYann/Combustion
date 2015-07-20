import pde
import numpy as np
import matplotlib.pyplot as plt
import pdb

from math import *
from scipy.interpolate import interp1d, griddata
from scipy.integrate import cumtrapz, quad

def int_diff(x1, y1, x2, y2):
    y2i = griddata(x2, y2, x1)
    err = np.abs(y1 - y2i)
    return cumtrapz(err, x1)[-1]

# coefficients
#a = -1.0
#eps = 30.0
#r = -400.0

a = -1.0
eps = 0.4
r = -7.0

Nx = 50
dt = (20.0/Nx)/5.0

FT = dt*10

N_A = 2

def solve(k):
   solns = []
   for j in range(k):
      (x, y) = pde.solve_pde(a, eps, r, Nx=Nx*2**j, dt=dt/(2**j),
                             FT=FT, N_A=4, do_linear=False, do_modified=False)
      print 'finished solve ', j+1
      solns.append((x,y))
   
   return solns

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

solns = solve(4)

e = errors(solns)
o = orders(e)

print e
print o
