import numpy as np
import matplotlib.pyplot as plt
from math import *
import pdb

# this is unfortunately hard coded according to the number of substeps...
# Gauss-Lobatto quadrature points
quad_pts = [0, (1-sqrt(1/5.0))*0.5, (1+sqrt(1/5.0))*0.5, 1]

# number of advection substeps
N_A = len(quad_pts)
# quadrature weights computed using Mathematica
weights = [[0.2206011330,   0.3793988670, -0.06781472846, 0.02060113296],
           [-0.07453559925, 0.5217491947,  0.5217491947, -0.07453559925],
           [0.02060113296, -0.06781472846, 0.3793988670,  0.2206011330]]
# perform our Gaussian sub-interval integration
def int_4(f, m, l):
    return np.dot(f, weights[m])*0.5*l

# calculate the value of the interpolating polynomial at an arbitrary time t
def interp(f, dt, t):
    return (2*pow(dt,3)*f[0] + 10*pow(t,3)*(-f[0] + sqrt(5)*f[1] - sqrt(5)*f[2] + f[3]) + 
    pow(dt,2)*t*(-12*f[0] + 5*(1 + sqrt(5))*f[1] + 5*f[2] - 5*sqrt(5)*f[2] + 2*f[3]) - 
    5*dt*pow(t,2)*(-4*f[0] + f[1] + 3*sqrt(5)*f[1] + f[2] - 3*sqrt(5)*f[2] + 2*f[3]))/(2.*pow(dt,3))

def rhs(f, y, r, dt, n, t):
    return interp(f, dt, t) + FR(r, y, n*dt + t)#r*y*(y-1)*(y-0.5)#pow(y,r)

def rk4(f, r, t_n, y_n, n, dt, dtp):
    M = 10
    h = dtp/float(M)
    
    y = y_n
    
    for t in np.linspace(t_n, t_n+dtp-h, M):
        k1 = h*(rhs(f, y, r, dt, n, t))
        k2 = h*(rhs(f, y+0.5*k1, r, dt, n, t+0.5*h))
        k3 = h*(rhs(f, y+0.5*k2, r, dt, n, t+0.5*h))
        k4 = h*(rhs(f, y+k3, r, dt, n, t+h))
        
        y_n = y
        y   = y + (k1 + 2*k2 + 2*k3 + k4)/6.0;
    
    return y

def FA(a, y, t):
    #return a*(cos(t) - y)
    return a*y

def FD(d, y, t):
    #return d*(cos(t) - y)
    return d*y

def FR(r, y, t):
    #return -sin(t) + r*(cos(t) - y)
    return r*y*(y-1)*(y-0.5)

def AD_RHS(a, d, y, t, dt):
    return map(lambda j:   FA(a, y[j], t+quad_pts[j]*dt) 
                         + FD(d, y[j], t+quad_pts[j]*dt),
               range(4))

def I_R(r, y, t, m, dt):
    f = map(lambda j: FR(r, y[j], t+quad_pts[j]*dt), range(4))
    return int_4(f, m, dt)

def I_AD(a, d, y, m, t, dt):    
    return int_4(AD_RHS(a, d, y, t, dt), m, dt)

def solve_it(a, d, r, dt, max_iter=10, plot_it=False):
    #number of timesteps to run
    #N = 1
    N = int(1/dt)
    #N = 1
    
    C = a + d + r
    
    y = np.zeros(N+1)
    T = np.linspace(0,N*dt,num=N+1)
    
    # initial condition
    y[0] = 1.0
    
    for n in range(N):
        # stupid predictor
        y_prev = np.array([y[n]]*N_A)
        y_curr = np.array([y[n]]*N_A)
        
        # corrector iterations
        for k in range(max_iter):
            y_AD = y[n]
            for m in range(0,N_A-1):
                dtp = dt*(quad_pts[m+1] - quad_pts[m])
                
                t = dt*(n + quad_pts[m])
                tt = t+dtp
                
                y_AD = (y_AD + dtp*(
                          FA(a,y_curr[m],t) - FA(a,y_prev[m],t)
#                        + d*cos(tt)
                        - FD(d,y_prev[m+1],tt))
                        + I_AD(a, d, y_prev, m, dt*n, dt)
                        + I_R(r, y_prev, dt*n, m ,dt))/(1 - dtp*d)
                
                f_const = (FD(d, y_AD, tt) - FD(d, y_prev[m+1], tt)
                         + FA(a, y_curr[m], t) - FA(a, y_prev[m], t))
                
                f = f_const + AD_RHS(a, d, y_prev, dt*n,dt)
                
                soln = rk4(f, r, quad_pts[m]*dt, y_curr[m], n, dt, dtp)
                
                y_curr[m+1] = soln
            
            # move on to the next MISDC iteration...
            y_prev = y_curr
        
        y[n+1] = y_prev[N_A - 1]
    
    #exact = np.exp((a-d+r)*T)
    #exact = np.cos(T)
    #l1exact = np.sum(abs(exact))
    #l1err = np.sum(np.abs(y-exact))
    
    #plt.plot(T, y, T, exact)
    plt.plot(T, y)
    plt.show()
    
    #return l1err/l1exact
