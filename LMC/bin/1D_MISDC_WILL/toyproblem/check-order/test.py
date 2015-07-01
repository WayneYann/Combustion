from math import *
import numpy as np
from scipy.interpolate import interp1d, griddata
from scipy.integrate import cumtrapz, quad
import pdb

import pde
import ml

FT = 1.0
Nx = 200
dt = 0.04

a = -1.0
eps = 20.0
r = -50.0

data = np.loadtxt('fine')
fine_x = data[:,0]
fine_y = data[:,1]

soln1 = pde.solve_pde(a, eps, r, Nx=Nx,   dt=dt,      FT=FT)
soln2 = pde.solve_pde(a, eps, r, Nx=Nx*2, dt=dt*0.5,  FT=FT)
soln3 = pde.solve_pde(a, eps, r, Nx=Nx*4, dt=dt*0.25, FT=FT)

fine_i1 = griddata(fine_x, fine_y, soln1[0])
fine_i2 = griddata(fine_x, fine_y, soln2[0])
fine_i3 = griddata(fine_x, fine_y, soln3[0])

err1 = np.abs(soln1[1] - fine_i1)
err2 = np.abs(soln2[1] - fine_i2)
err3 = np.abs(soln3[1] - fine_i3)

numer = cumtrapz(err1, soln1[0])[-1]
denom = cumtrapz(err2, soln2[0])[-1]

print 'error: ', numer
print 'error: ', denom
print 'order: ', log(numer/denom)/log(2)

numer = denom
denom = denom = cumtrapz(err3, soln3[0])[-1]
print 'error: ', numer
print 'error: ', denom
print 'order: ', log(numer/denom)/log(2)

soln1 = ml.solve_pde(a, eps, r, Nx=Nx,   dt=dt,     FT=FT)
soln2 = ml.solve_pde(a, eps, r, Nx=Nx*2, dt=dt*0.5, FT=FT)
soln3 = ml.solve_pde(a, eps, r, Nx=Nx*4, dt=dt*0.25, FT=FT)

fine_i1 = griddata(fine_x, fine_y, soln1[0])
fine_i2 = griddata(fine_x, fine_y, soln2[0])
fine_i3 = griddata(fine_x, fine_y, soln3[0])

err1 = np.abs(soln1[1] - fine_i1)
err2 = np.abs(soln2[1] - fine_i2)
err3 = np.abs(soln3[1] - fine_i3)

numer = cumtrapz(err1, soln1[0])[-1]
denom = cumtrapz(err2, soln2[0])[-1]

print 'error: ', numer
print 'error: ', denom
print 'order: ', log(numer/denom)/log(2)

numer = denom
denom = cumtrapz(err3, soln3[0])[-1]

print 'error: ', numer
print 'error: ', denom
print 'order: ', log(numer/denom)/log(2)

def do_errors(soln1, soln2, soln3):
    
    soln1i = griddata(soln1[0], soln1[1], soln2[0])
    soln2i = griddata(soln2[0], soln2[1], soln3[0])     
    
    err1 = np.abs(soln2[1] - soln1i)
    err2 = np.abs(soln3[1] - soln2i)
    
    numer = cumtrapz(err1, soln2[0])[-1]
    denom = cumtrapz(err2, soln3[0])[-1]

    print 'error: ', numer
    print 'error: ', denom

    print 'order: ', log(numer/denom)/log(2)
    print

def do_solns(Nxratio, dtratio):
    print 'ODE with constant  coefficients'
    soln1 = pde.solve_pde(a, eps, r, Nx=Nx,            dt=dt, FT=FT)
    soln2 = pde.solve_pde(a, eps, r, Nx=Nx*Nxratio,    dt=dt*dtratio,  FT=FT)
    soln3 = pde.solve_pde(a, eps, r, Nx=Nx*Nxratio**2, dt=dt*dtratio**2,  FT=FT)
    do_errors(soln1, soln2, soln3)
    
    print 'MISDC with Gauss-Lobatto. 1 iteration'
    soln1 = ml.solve_pde(a, eps, r, Nx=Nx, dt=dt, FT=FT,
                         max_iter=1)
    soln2 = ml.solve_pde(a, eps, r, Nx=Nx*Nxratio,    dt=dt*dtratio,  FT=FT,
                         max_iter=1)
    soln3 = ml.solve_pde(a, eps, r, Nx=Nx*Nxratio**2, dt=dt*dtratio**2,  FT=FT,
                         max_iter=1)
    do_errors(soln1, soln2, soln3)
    
    print 'MISDC with Gauss-Lobatto. 2 iterations'
    soln1 = ml.solve_pde(a, eps, r, Nx=Nx, dt=dt, FT=FT,
                         max_iter=2)
    soln2 = ml.solve_pde(a, eps, r, Nx=Nx*Nxratio,    dt=dt*dtratio,  FT=FT,
                         max_iter=2)
    soln3 = ml.solve_pde(a, eps, r, Nx=Nx*Nxratio**2, dt=dt*dtratio**2,  FT=FT,
                         max_iter=2)
    do_errors(soln1, soln2, soln3)


#print ' --- Refining in time --- '
#print
#do_solns(1,0.5)
#print ' --- Refining in space and time --- '
#print
#do_solns(2,0.5)
