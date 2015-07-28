from math import *
import numpy as np
from scipy.interpolate import interp1d, griddata
from scipy.integrate import cumtrapz, quad
import matplotlib.pyplot as plt
import pdb

import pde
import ml

FT = 1.0
Nx = 40
dt = 0.1

a = -1.0
eps = 20.0
r = -50.0

data = np.loadtxt('fine')
fine_x = data[:,0]
fine_y = data[:,1]

def int_diff(x1, y1, x2, y2):
    y2i = griddata(x2, y2, x1)
    err = np.abs(y1 - y2i)
    return cumtrapz(err, x1)[-1]

def int_error(mod, Nxx, dtt, max_iter=2):
    print "Nx: ", Nxx
    print "dt: ", dtt
    soln = mod.solve_pde(a, eps, r, Nx=Nxx, dt=dtt, FT=FT, max_iter=max_iter)
    
    return int_diff(soln[0], soln[1], fine_x, fine_y)

def compute_order(mod, Nxr, dtr, max_iter=2):
    e1 = int_error(mod, Nx,         dt,       max_iter=max_iter)
    e2 = int_error(mod, Nx*Nxr,     dt*dtr,   max_iter=max_iter)
    e3 = int_error(mod, Nx*Nxr**2, dt*dtr**2, max_iter=max_iter)
    e4 = int_error(mod, Nx*Nxr**3, dt*dtr**3, max_iter=max_iter)

    print 'error1: ', e1
    print 'error2: ', e2
    print 'error3: ', e3
    print 'error4: ', e4
    print 'order: ', log(e1/e2)/log(2)
    print 'order: ', log(e2/e3)/log(2)
    print 'order: ', log(e3/e4)/log(2)
    
    return (e1, e2, e3, e4)

def compare(soln1, soln2, soln3):
    return log(int_diff(soln1[0], soln1[1], soln2[0], soln2[1])/
               int_diff(soln2[0], soln2[1], soln3[0], soln3[1]))/log(2)

soln1 = pde.solve_pde(a, eps, r, Nx*4, dt, FT=FT)
soln2 = pde.solve_pde(a, eps, r, Nx*4, dt/2.0, FT=FT)
soln3 = pde.solve_pde(a, eps, r, Nx*4, dt/4.0, FT=FT)
print 'order: ', compare(soln1, soln2, soln3)
soln4 = pde.solve_pde(a, eps, r, Nx*4, dt/8.0, FT=FT)
print 'order: ', compare(soln2, soln3, soln4)
soln5 = pde.solve_pde(a, eps, r, Nx*16, dt/16.0, FT=FT)
print 'order: ', compare(soln3, soln4, soln5)
soln6 = pde.solve_pde(a, eps, r, Nx*32, dt/32.0, FT=FT)
print 'order: ', compare(soln4, soln5, soln6)
soln7 = pde.solve_pde(a, eps, r, Nx*64, dt/64.0, FT=FT)
print 'order: ', compare(soln5, soln6, soln7)

soln1 = ml.solve_pde(a, eps, r, Nx, dt, FT=FT)
soln2 = ml.solve_pde(a, eps, r, Nx*2, dt/2.0, FT=FT)
soln3 = ml.solve_pde(a, eps, r, Nx*4, dt/4.0, FT=FT)
print 'order: ', compare(soln1, soln2, soln3)
soln4 = ml.solve_pde(a, eps, r, Nx*8, dt/8.0, FT=FT)
print 'order: ', compare(soln2, soln3, soln4)
soln5 = ml.solve_pde(a, eps, r, Nx*16, dt/16.0, FT=FT)
print 'order: ', compare(soln3, soln4, soln5)
soln6 = ml.solve_pde(a, eps, r, Nx*32, dt/32.0, FT=FT)
print 'order: ', compare(soln4, soln5, soln6)
soln7 = ml.solve_pde(a, eps, r, Nx*64, dt/64.0, FT=FT)
print 'order: ', compare(soln5, soln6, soln7)
##Nx = 800
##pde_e = compute_order(pde, 1, 0.5)
##dt = 0.3
##ml1_e = compute_order(ml, 1, 0.5, max_iter=1)
##ml2_e = compute_order(ml, 1, 0.5, max_iter=2)

##gp = (Nx, Nx*2, Nx*4, Nx*8)

##pdb.set_trace()

#plt.clf()
#plt.loglog(gp, pde_e, label="PDE")
#plt.loglog(gp, ml1_e, label="GL 1")
#plt.loglog(gp, ml2_e, label="GL 2")
#plt.legend()
#plt.show(block=True)
