import numpy as np
import matplotlib.pyplot as plt
from math import *
import pdb

def solve_it(d, r, dt, max_iter=10):
    # final time
    #T = 3
    # number of timesteps to run
    #N = int(T/dt)
    
    N = 1
    
    C = d + r

    y = np.zeros(N+1)
    t = np.linspace(0,N*dt,num=N+1)

    residuals = np.ones(max_iter+1)

    # initial condition
    y[0] = 1.0
    
    for n in range(N):
        # stupid predictor
        y_prev = y[n]
        
        I_R = dt*r*y[n]
        
        # corrector iterations
        for k in range(max_iter):
            I = 0.5*dt*(y[n]+y_prev)
            
            y_AD = (y[n] + dt*(-d*y_prev) + d*I + I_R)/(1-d*dt)
            #y_AD1 = (y[n] + dt1*(-d*y_prev1) + (d+r)*I_1)/(1-d*dt1)
            
            #c = d*y_AD1 - d*y_prev1 + d*I_1/dt1
            c = d*y_AD - d*y_prev + d*I/dt
            y_next = (-c + c*exp(r*dt))/r + exp(r*dt)*y[n]
            
            I_R = y_next - y[n] - c*dt
            
            y_prev = y_next
        
        y[n+1] = y_next
    
    print 'soln:  ', y[N]
    print 'exact: ', np.exp(C*dt)
    
    exact = np.exp(C*t)
    l1exact = np.sum(exact)*dt
    l1err = np.sum(np.abs(y-exact))*dt
        
    print 'iters: ', max_iter, ' relative L^1 error: ', l1err/l1exact
        
    return l1err/l1exact
