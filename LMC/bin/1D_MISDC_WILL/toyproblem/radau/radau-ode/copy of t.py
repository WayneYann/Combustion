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

    dt1 = dt/3.0
    dt2 = 2.0*dt/3.0
    
    for n in range(N):
        # stupid predictor
        y_prev1 = y[n]
        y_prev2 = y[n]
        
        I_R1 = dt1*r*y[n]
        I_R2 = dt2*r*y[n]
        # corrector iterations
        for k in range(max_iter):
            #I_1 = (1.0/9.0)*(5*y_prev1 - 2*y_prev2)*dt
            #I_2 = (2.0/9.0)*(2*y_prev1 +   y_prev2)*dt
            
            #I_1 = 0.0625*(9*y_prev1 - y_prev2)*dt
            #I_2 = 0.0625*(3*y_prev1 + 5*y_prev2)*dt
            
            I_1 =  (1/27.0)*(7*y[n] + 15*y_prev1 -  4*y_prev2)*dt
            I_2 = (1/108.0)*( -y[n] + 21*y_prev1 + 16*y_prev2)*dt
            
            #I_1 = 0.5*dt1*(y[n]+y_prev1)
            #I_2 = 0.5*dt2*(y_prev1 + y_prev2)
            
            #print 'diff1: ', abs(I_1 - I_1_o)
            
            #y_AD1 = (y[n] + dt1*(-d*y_prev1) + d*I_1 + I_R1)/(1-d*dt1)
            y_AD1 = (y[n] + dt1*(-d*y_prev1) + (d+r)*I_1)/(1-d*dt1)
            c = d*y_AD1 - d*y_prev1 - r*y_prev1 + (d+r)*I_1/dt1
            y_next1 = (-c*dt1 - y[n])/(r*dt1 - 1)
            
            #y_next1 = (-c + c*exp(r*dt1))/r + exp(r*dt1)*y[n]
            
            #I_R1 = y_next1 - y[n] - c*dt1
            
            #y_next1_old = (y[n] + dt1*(d*y_AD1 - d*y_prev1 - r*y_prev1)
            #               + (d+r)*I_1)/(1-r*dt1)
            
            #y_AD2 = (y_AD1 + dt2*(-d*y_prev2) + d*I_2 + I_R2)/(1-d*dt2)
            
            y_AD2 = (y_AD1 + dt2*(-d*y_prev2) + (d+r)*I_2)/(1-d*dt2)
            c = d*y_AD2 - d*y_prev2 - r*y_prev1 + (d+r)*I_2/dt2
            y_next2 = (-c*dt2 - y_next1)/(r*dt2 - 1)
            
            #y_next2 = (-c + c*exp(r*dt2))/r + exp(r*dt2)*y_next1
            #I_R2 = y_next2 - y_next1 - c*dt2
            #y_next2_old = (y_next1 + dt2*(d*y_AD2 - d*y_prev2 - r*y_prev2)
            #                   + (d+r)*I_2)/(1-r*dt2)
            
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
