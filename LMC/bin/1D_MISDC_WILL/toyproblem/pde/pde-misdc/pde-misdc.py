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
# number of advection substeps
N_A = 2
# number of reaction substeps
N_R = 4
# piecewise linear or constants
do_linear = True

endpt = 20.0
max_iter = 2
# coefficients
#a = -1.0
#eps = 10.0
##r = -400.0
#r = -200.0

a = -1.0
eps = 20.0
r = -300.0

# gridsize
Nx = 400
# spacing
h = float(endpt/Nx)
# timestep
dt = h/5
# final time
FT = 2.0
# number of timesteps
Nt = int(ceil(FT/dt))

a_quad_pts = [0, 1]

if N_R == 4:
    # Gauss-Lobatto quadrature points
    r_quad_pts = [0, (1-sqrt(1/5.0))*0.5, (1+sqrt(1/5.0))*0.5, 1]
    # quadrature weights computed using Mathematica
    weights = [[0.2206011330,   0.3793988670, -0.06781472846, 0.02060113296],
               [-0.07453559925, 0.5217491947,  0.5217491947, -0.07453559925],
               [0.02060113296, -0.06781472846, 0.3793988670,  0.2206011330],
               [1.0/6.0, 5.0/6.0, 5.0/6.0, 1.0/6.0]]
    weights_A = [(1-sqrt(1/5.0))*0.5, sqrt(1.0/5.0), (1-sqrt(1/5.0))*0.5, 1.0]
else:
    r_quad_pts = [0, 1]

# perform our Gaussian sub-interval integration
def int_4(f, m, l):
    return np.apply_along_axis(lambda f_j: np.dot(f_j, weights[m])*0.5*l, 0, f)
def int_2(f, m, l):
    return 0.5*(f[0]+f[1])*l
def integ(f, m, l):
    if N_R==4:
        return int_4(f, m, l)
    elif N_R==2:
        return int_2(f, m, l)

# calculate the value of the interpolating polynomial at an arbitrary time t
def interp_4(f, t):
    return (2*pow(dt,3)*f[0] + 10*pow(t,3)*(-f[0] + sqrt(5)*f[1] - sqrt(5)*f[2] 
            + f[3]) + pow(dt,2)*t*(-12*f[0] + 5*(1 + sqrt(5))*f[1] + 5*f[2]
            - 5*sqrt(5)*f[2] + 2*f[3]) - 5*dt*pow(t,2)*(-4*f[0] + f[1] 
            + 3*sqrt(5)*f[1] + f[2] - 3*sqrt(5)*f[2] + 2*f[3]))/(2.*pow(dt,3))

def interp_const(f, t):
    return 0.5*(f[0]+f[1])

def interp_linear(f, t):
    return (1-t/dt)*f[0] + (t/dt)*f[1]

def interp(f, t):
    if N_R == 4:
        return interp_4(f, t)
    elif N_R == 2:
        if do_linear:
            return interp_linear(f, t)
        else:
            return interp_const(f, t)

def rhs(f, z, t):
    return interp(f,t) + FR(z)

def rk4(f, t_n, y_n, dtp):
    M = 10
    h_rk = dtp/float(M)
    
    z = y_n
    
    for t in np.linspace(t_n, t_n+dtp-h_rk, M):
        k1 = h_rk*(rhs(f, z, t))
        k2 = h_rk*(rhs(f, z+0.5*k1, t+0.5*h_rk))
        k3 = h_rk*(rhs(f, z+0.5*k2, t+0.5*h_rk))
        k4 = h_rk*(rhs(f, z+k3, t+h_rk))
        
        z  = z + (k1 + 2*k2 + 2*k3 + k4)/6.0;
    return z

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
    ADz = AD_RHS(z)
    
    return 0.5*(ADz[0]+ADz[3])*dt*weights_A[m]

    return integ(ADz, m, dt)

def I_R(z, m):
    Rz = np.apply_along_axis(lambda z_m: FR(z_m), 1, z)
    
    #return 0.5*(Rz[0]+Rz[3])*dt
    
    return integ(Rz, m, dt)

def advance(n):
    # boundary values
    y[n+1][0] = y[n][0]
    y[n+1][-1] = y[n][-1]
        
    # stupid predictor
    y_prev = np.array([y[n][1:-1]]*N_R)
    y_curr = np.copy(y_prev)
    
    # corrector iterations
    for k in range(max_iter):
        y_AD = y[n][1:-1]
        
        for m in range(N_A-1):
            dtp = dt*(a_quad_pts[m+1] - a_quad_pts[m])
            
            y_AD = spsolve(Dinv[m],
                        y_AD + dtp*(FA(y_curr[m*N_R]) - FA(y_prev[m*N_R])
                                  + eps*bc/h**2       - FD(y_prev[m*N_R+N_R-1]))
                             + I_AD(y_prev, -1)
                             + I_R(y_prev, -1))
            
            y_q = np.copy(y_curr[m*N_R])
            for q in range(N_R-1):
                dtq = dtp*(r_quad_pts[q+1] - r_quad_pts[q])
                
                b = (y_q + dtq*(FA(y_curr[m*N_R]) - FA(y_prev[m*N_R])
                              + FD(y_AD)          - FD(y_prev[m*N_R + N_R - 1])
                              - FR(y_prev[m*N_R + q+1]))
                         + I_AD(y_prev, q)
                         + I_R(y_prev, q))
                
                for j in range(Nx - 2):
                    soln = newton(lambda x: x - dtq*FR(x) - b[j],
                                  y_prev[m*N_R + q + 1][j],
                                  fprime = lambda x: 1 - dtq*FRprime(x),
                                  fprime2 = lambda x: dtq*FRprime2(x),
                                  maxiter = 1000)
                    y_q[j] = soln
                
                y_curr[m*N_R + q+1] = y_q
        
        # move on to the next MISDC iteration...
        y_prev = np.copy(y_curr)
    y[n+1][1:-1] = y_curr[N_R - 1]

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
Dinv = np.empty(N_A - 1, dtype=object)
for m in range(N_A - 1):
    dtp = dt*(a_quad_pts[m+1] - a_quad_pts[m])
    Dinv[m] = csc_matrix(I - dtp*D)

# create the boundary condition
bc = np.zeros(Nx-2)
bc[0] = 1
bc[-1] = 0

plt.ion()
step = 10
for n in range(Nt):
    advance(n)
    if (n+1)%step == 0 or n == Nt-1:
        print 't = ', dt*(n+1)
        plt.cla()
        plt.plot(x,y[n+1])
        plt.ylim((0,1.1))
        plt.draw()

np.savetxt('ft', np.array([x, y[Nt]]).T)

plt.show(block=True)
