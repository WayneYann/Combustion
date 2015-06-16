import numpy as np
import matplotlib.pyplot as plt
from math import *

# this is unfortunately hard coded according to the number of substeps...
# Gauss-Lobatto quadrature points
quad_pts = [0, (1-sqrt(1/5.0))*0.5, (1+sqrt(1/5.0))*0.5, 1]
# number of substeps
M = len(quad_pts)
# quadrature weights computed using Mathematica
weights = [[0.2206011330,   0.3793988670, -0.06781472846, 0.02060113296],
           [-0.07453559925, 0.5217491947,  0.5217491947, -0.07453559925],
           [0.02060113296, -0.06781472846, 0.3793988670,  0.2206011330]]
# perform our Gaussian sub-interval integration
def int_4(f, m, l):
    return np.dot(f, weights[m])*0.5*l

def interp(f, dt, t):
    return (2*pow(dt,3)*f[0] + 10*pow(t,3)*(-f[0] + sqrt(5)*f[1] - sqrt(5)*f[2] + f[3]) + 
    pow(dt,2)*t*(-12*f[0] + 5*(1 + sqrt(5))*f[1] + 5*f[2] - 5*sqrt(5)*f[2] + 2*f[3]) - 
    5*dt*pow(t,2)*(-4*f[0] + f[1] + 3*sqrt(5)*f[1] + f[2] - 3*sqrt(5)*f[2] + 2*f[3]))/(2.*pow(dt,3))

def rk4(f, r, y_n, a, dt, dtp):
    M = 10#int(1/dt)
    h = dtp/float(M)
    
    y = y_n
    i = 0
    
    for t in np.linspace(a, a+dtp, M):
        k1 = h*(interp(f, dt, t)       + r*y)
        k2 = h*(interp(f, dt, t+0.5*h) + r*(y + 0.5*k1))
        k3 = h*(interp(f, dt, t+0.5*h) + r*(y + 0.5*k2))
        k4 = h*(interp(f, dt, t+h)     + r*(y + k3))
        
        y_n = y
        y = y + (k1 + 2*k2 + 2*k3 + k4)/6.0;

        i += 0.5*(y+y_n)*h
    
    return y
    #return (y, i)

def solve_it(C, dt, max_iter=10, plot_it=False):
    # final time
    #T = 3
    # number of timesteps to run
    #N = int(T/dt)
    
    N = 1

    y = np.zeros(N+1)
    t = np.linspace(0,N*dt,num=N+1)

    residuals = np.ones(max_iter+1)

    # initial condition
    y[0] = 1.0

    for n in range(N):
        # stupid predictor
        y_prev = np.array([y[n]]*M)
        y_next = y_prev
        
        delta = 0
        
        # corrector iterations
        for i in range(max_iter):
            for m in range(M-1):
                dtp = dt*(quad_pts[m+1] - quad_pts[m])
                
                f = C*y_prev - C*y_prev[m]
                f = np.zeros(4)
                
                y_next[m+1] = rk4(f, C, y_prev[m], quad_pts[m], dt, dtp)
                
                #y_next[m+1] = (y_next[m] - dtp*C*y_prev[m+1] + C*int_4(y_prev, m, dt))/(1 - dtp*C)
                
            y_prev = y_next
        
        y[n+1] = y_next[M-1]
    
    exact = np.exp(C*t)
    
    #print exact
    #print y
    
    l1exact = np.sum(exact)*dt
    l1err = np.sum(np.abs(y-exact))*dt
    
    print 'max_iter: ', max_iter, ' error: ', l1err
    
    if plot_it:
        print 'relative L^1 error: ', l1err/l1exact
        plt.plot(t,y,t,exact)
        plt.show()
    return l1err/l1exact

    #print 'relative L^1 error: ', l1err/l1exact
    #print 'dt^2: ', dt**2
    #print 'errors: ', errors
    #print residuals

    #plt.plot(t,y,t,exact)
    #plt.plot(residuals)
    #plt.plot(errors)
    #plt.show()
