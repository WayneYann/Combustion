import numpy as np
import matplotlib.pyplot as plt
from math import *
import pdb

def rk_rhs(y, c, r, t):
   return (c+r*y)/(1-r*t)

def rk4(c, r, y_n, dt):
    M = 1
    h = dt/float(M)
    
    y = y_n
    i = 0
    
    for t in np.linspace(0, dt-h, M):
        k1 = h*rk_rhs(y,        c, r, t)
        k2 = h*rk_rhs(y+0.5*k1, c, r, t+0.5*h)
        k3 = h*rk_rhs(y+0.5*k2, c, r, t+0.5*h)
        k4 = h*rk_rhs(y+k3,     c, r, t+h)
        
        y   = y + (k1 + 2*k2 + 2*k3 + k4)/6.0;
    
    return y

def solve_it(d, r, dt, max_iter=10):
    # final time
    #T = 3
    # number of timesteps to run
    #N = int(T/dt)
    
    N = int(1.0/dt)
    #N = 1
    
    C = d + r

    y = np.zeros(N+1)
    t = np.linspace(0,N*dt,num=N+1)

    residuals = np.ones(max_iter+1)

    # initial condition
    y[0] = 1.0

    #dt1 = dt/3.0
    #dt2 = 2.0*dt/3.0
    
    dt1 = dt/2.0
    dt2 = dt/2.0
    
    for n in range(N):
        # stupid predictor
        y_prev1 = y[n]
        y_prev2 = y[n]
        
        I_R1 = dt1*r*y[n]
        I_R2 = dt2*r*y[n]
        # corrector iterations
        for k in range(max_iter):
            #I_1 =  (1/27.0)*(7*y[n] + 15*y_prev1 -  4*y_prev2)*dt
            #I_2 = (1/108.0)*( -y[n] + 21*y_prev1 + 16*y_prev2)*dt
            
            I_1 = (5*y[n]/24.0 + y_prev1/3.0 - y_prev2/24.0)*dt
            I_2 = (-y[n]/24.0 + y_prev1/3.0 + 5*y_prev2/24.0)*dt
            
            y_AD1 = (y[n] + dt1*(-d*y_prev1) + (d+r)*I_1)/(1-d*dt1)
            c = d*y_AD1 - d*y_prev1 - r*y_prev1 + (d+r)*I_1/dt1
            #y_next1 = (-c*dt1 - y[n])/(r*dt1 - 1)
            y_next1 = rk4(c, r, y[n], dt1)
            
            y_AD2 = (y_AD1 + dt2*(-d*y_prev2) + (d+r)*I_2)/(1-d*dt2)
            c = d*y_AD2 - d*y_prev2 - r*y_prev2 + (d+r)*I_2/dt2
            #y_next2 = (-c*dt2 - y_next1)/(r*dt2 - 1)
            y_next2 = rk4(c, r, y_next1, dt2)
            
            y_prev1 = y_next1
            y_prev2 = y_next2
        
        y[n+1] = y_next2
    
    print 'soln:  ', y[N]
    print 'exact: ', np.exp(C*dt)
    
    exact = np.exp(C*t)
    l1exact = np.sum(exact)*dt
    l1err = np.sum(np.abs(y-exact))*dt
        
    print 'iters: ', max_iter, ' relative L^1 error: ', l1err/l1exact
        
    return l1err/l1exact
