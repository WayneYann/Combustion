from __future__ import division
import numpy as np
import sys
import pdb
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
from scipy.interpolate import interp1d
from scipy.optimize import newton
from math import *

###########################
# here are our parameters #
###########################

def solve_pde(a, eps, r, Nx=300, dt=None, FT=1.0, max_iter=1):
    # coefficients
    #a = -1.0
    #eps = 600.0
    #r = -1100.0

    #a = -1.0
    #eps = 0.4
    #r = -7.0

    # gridsize spacing
    endpt = 20.0
    h = float(endpt/Nx)
    # default timestep
    if not dt:
        dt = h/25.0
    else:
        dt = dt/5.0
    
    dt1 = dt/3.0
    dt2 = 2*dt/3.0
    # final time FT
    # number of timesteps
    Nt = int(ceil(FT/dt))

    # perform our Gaussian sub-interval integration
    def integ(f0, f1, f2, m):
        if m==1:
            return  (1/27.0)*(7*f0 + 15*f1 -  4*f2)*dt
        elif m==2:
            return (1/108.0)*( -f0 + 21*f1 + 16*f2)*dt

    def FR(z):
        return r*z*(z-1)*(z-0.5)

    def FRprime(z):
        return r*(3*z**2 - 3*z + 0.5)

    def FRprime2(z):
        return r*(6*z - 3)

    def FA(z):
        return np.dot(A, z) - a*bc/(2*h)

    def FD(z):
        return np.dot(D, z) + eps*bc/h**2

    def AD_RHS(z):
        return np.apply_along_axis(lambda z_m: FA(z_m) + FD(z_m), 1, z)

    def I_AD(z, m):    
        return integ(AD_RHS(z), m, dt)

    def I_R(z, m):
        Rz = np.apply_along_axis(lambda z_m: FR(z_m), 1, z)
        
        return integ(Rz, m, dt)

    def advance(n):
        # boundary values
        y[n+1][0] = y[n][0]
        y[n+1][-1] = y[n][-1]
            
        # stupid predictor
        y_prev1 = np.array(y[n][1:-1])
        y_prev2 = np.copy(y_prev1)
        y_curr1 = np.copy(y_prev1)
        y_curr2 = np.copy(y_prev1)
        yn = np.copy(y_prev1)
        
        
        # corrector iterations
        for k in range(max_iter):
            I_1 = integ(FA(yn)      + FD(yn)      + FR(yn),
                        FA(y_prev1) + FD(y_prev1) + FR(y_prev1),
                        FA(y_prev2) + FD(y_prev2) + FR(y_prev2),
                        1)
            y_AD1 = spsolve(Dinv1, yn + dt1*(eps*bc/h**2 - FD(y_prev1)) + I_1)
            
            rhs = (yn + dt1*(FD(y_AD1) - FD(y_prev1) - FR(y_prev1)) + I_1)
            
            for j in range(Nx - 2):
                soln = newton(lambda x: x - dt1*FR(x) - rhs[j],
                          y_curr1[j],
                          fprime = lambda x: 1 - dt1*FRprime(x),
                          fprime2 = lambda x: -dt1*FRprime2(x),
                          maxiter = 500)
                y_curr1[j] = soln
            
            I_2 = integ(FA(yn) + FD(yn) + FR(yn),
                        FA(y_prev1) + FD(y_prev1) + FR(y_prev1),
                        FA(y_prev2) + FD(y_prev2) + FR(y_prev2),
                        2)
            y_AD2 = spsolve(Dinv2, y_AD1 + dt2*(FA(y_curr1) - FA(y_prev1)
                                              + eps*bc/h**2 - FD(y_prev2)) + I_2)
            
            rhs = (y_curr1 + dt2*(FA(y_curr1) - FA(y_prev1)
                                + FD(y_AD2)   - FD(y_prev2) - FR(y_prev2)) + I_2)
            
            for j in range(Nx - 2):
                soln = newton(lambda x: x - dt2*FR(x) - rhs[j],
                          y_curr2[j],
                          fprime = lambda x: 1 - dt2*FRprime(x),
                          fprime2 = lambda x: -dt2*FRprime2(x),
                          maxiter = 500)
                y_curr2[j] = soln
            
            # move on to the next MISDC iteration...
            y_prev1 = np.copy(y_curr1)
            y_prev2 = np.copy(y_curr2)
        
        y[n+1][1:-1] = y_curr2

    x = np.linspace(0, endpt, num=Nx)
    y = np.zeros((Nt+2, Nx))

    y[0] = 0.5*(np.tanh(10-2*x)+1)

    # identity
    I = np.eye(Nx-2)
    # create the differentiation matrices
    Dx = (np.diag([1]*(Nx-3),1) - np.diag([1]*(Nx-3),-1))/(2*h)
    # laplacian operator
    Dxx = (-I*2 + np.diag([1]*(Nx-3),1) + np.diag([1]*(Nx-3),-1))/h**2

    # create the advection and diffusion operators
    # diffusion is implicit, so this is the (I - D)
    D = eps*Dxx
    A = a*Dx

    # we treat diffusion implicitly
    # for efficiency, we create the sparse matrices now
    # which we later use in the sparse solve
    # (I - dtp D) y = RHS
    Dinv1 = csc_matrix(I - dt1*D)
    Dinv2 = csc_matrix(I - dt2*D)

    # create the boundary condition
    bc = np.zeros(Nx-2)
    bc[0] = 1
    bc[-1] = 0

    plt.ion()
    skip = 1000
    for n in range(Nt):
        advance(n)
        if (n+1)%skip == 0 or n == Nt-1:
#            print 't = ', dt*(n+1)
            plt.cla()
            plt.plot(x,y[n+1])
            plt.ylim((0,1.1))
            plt.draw()
    
    soln = np.array([x, y[Nt]])
    np.savetxt('ft', soln.T)
    
    return soln

    #plt.show(block=True)
