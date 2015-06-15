import numpy as np
import matplotlib.pyplot as plt
from math import *

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

def solve_it(a, d, r, dt, max_iter=10, do_linear=True, plot_it=False):
    #number of timesteps to run
    #N = int(2/dt)
    N = 1
    
    C = a + d + r
    
    y = np.zeros(N+1)
    t = np.linspace(0,N*dt,num=N+1)
    
    # initial condition
    y[0] = 1.0
    
    for n in range(N):
        # stupid predictor
        y_prev = [y[n]]*N_A
        y_curr = [y[n]]*N_A
        
        # less stupid predictor
        for m in range(0, N_A-1):
            dtp = dt*(quad_pts[m+1] - quad_pts[m])
            y_AD = (y_prev[m] + dtp*a*y[n])/(1 - d*dtp)
            y_prev[m+1] = (y_prev[m] + dtp*(a*y_prev[m] + d*y_AD))/(1 - r*dtp)
        
        # corrector iterations
        for k in range(max_iter):
            y_AD = y[n]
            # run the advection substeps
            for m in range(0,N_A-1):
                dtp = dt*(quad_pts[m+1] - quad_pts[m])
                
                y_AD = (y_AD + dtp*(
                          a*y_curr[m] - a*y_prev[m] 
                        - d*y_prev[m+1])
                        + (a+d+r)*int_4(y_prev, m, dt))/(1 - dtp*d)
                
                y_curr[m+1] = (y_curr[m] + dtp*(a*y_curr[m] - a*y_prev[m]
                               + d*y_AD
                               -  d*y_prev[m+1] - r*y_prev[m+1])
                               + (a+d+r)*int_4(y_prev, m, dt))/(1 - r*dtp)
            # move on to the next MISDC iteration...
            y_prev = y_curr
        
        y[n+1] = y_prev[N_A - 1]
    
    exact = np.exp(C*t)
    l1exact = np.sum(exact)
    l1err = np.sum(np.abs(y-exact))
    
    if plot_it:
        print 'relative L^1 error: ', l1err/l1exact
        plt.plot(t,y,t,exact)
        plt.show()
    
    print 'L^1 error: ', l1err/l1exact
    return l1err/l1exact
