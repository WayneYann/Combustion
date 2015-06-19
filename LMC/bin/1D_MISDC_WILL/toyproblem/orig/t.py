import numpy as np
import matplotlib.pyplot as plt
from math import *

def solve_it(a, d, r, dt, max_iter=10, do_linear=True, plot_it=False):
    # final time
    #T = 3
    # number of timesteps to run
    #N = int(T/dt)
    
    N = 1

    C = a + d + r

    y = np.zeros(N+1)
    t = np.linspace(0,N*dt,num=N+1)

    residuals = np.ones(max_iter+1)

    # initial condition
    y[0] = 1.0

    for n in range(N):
        A_n = a*y[n]
        D_n = d*y[n]
        I_R = r*y[n]*dt
        
        # stupid predictor
        y_prev = y[n]
        residuals[0] = dt*C*y[n]
        
        # corrector iterations
        for i in range(max_iter):
            A = y_prev*a
            D = y_prev*d
            
            y_AD = (y[n] + dt*0.5*(A_n + A + D_n - D) + I_R)/(1 - dt*d)
            
            
            # compute the forcing term
            f_const = d*y_AD - d*y_prev
            if do_linear:
                f_lin_prev = A_n + D_n
                f_lin_next = A + D
            else:
                avg = 0.5*(A_n + D_n + A + D)
                f_lin_prev = avg
                f_lin_next = avg
            
            # compute the slope and constant term
            m = (f_lin_next - f_lin_prev)/dt
            c = f_const + f_lin_prev
            
            # solve the ODE here
            # solving y_next' = r*y_next + c + mt
            y_next = ((exp(r*dt)*(m+c*r) - m - c*r - m*r*dt)/r**2 + exp(r*dt)*y[n])
            
            # compute the effect of the reaction term by integrating
            #I_R = (y_next - y[n])/dt - f_const - 0.5*(f_lin_prev + f_lin_next)
            I_R = (y_next - y[n]) - f_const*dt - 0.5*dt*(f_lin_prev + f_lin_next)
            y_prev = y_next
            residuals[i+1] = y[n] - y_next + dt*I_R + 0.5*dt*(a+d)*(y[n] + y_next)
            
            if False:
                print 'i=',i,'  y_next: ', y_next
                print 'i=',i,'  I_R:    ', I_R
                print 'i=',i,'  y_AD:   ', y_AD
                print 'i=',i,'  resid:  ', residuals[i+1]
                print
        
        y[n+1] = y_next
    
    
    exact = np.exp(C*t)
    l1exact = np.sum(exact)*dt
    l1err = np.sum(np.abs(y-exact))*dt
    
    
    if plot_it:
        print 'relative L^1 error: ', l1err/l1exact
        plt.plot(t,y,t,exact)
        plt.show()
        
    print 'relative L^1 error: ', l1err/l1exact
        
    return l1err/l1exact
    
    #print 'dt^2: ', dt**2
    #print 'errors: ', errors
    #print residuals

    #plt.plot(t,y,t,exact)
    #plt.plot(residuals)
    #plt.plot(errors)
    #plt.show()
