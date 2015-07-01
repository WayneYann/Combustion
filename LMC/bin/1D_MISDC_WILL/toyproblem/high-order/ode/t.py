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

def ode_solve_pt(f, r, t_n, y_n, dt, t):
    return (-60*sqrt(5)*f[2]*exp(r*t) + 60*f[3]*exp(r*t) - 10*f[2]*dt*exp(r*t)*r + 
    30*sqrt(5)*f[2]*dt*exp(r*t)*r - 20*f[3]*dt*exp(r*t)*r + 
    5*f[2]*pow(dt,2)*exp(r*t)*pow(r,2) - 
    5*sqrt(5)*f[2]*pow(dt,2)*exp(r*t)*pow(r,2) + 
    2*f[3]*pow(dt,2)*exp(r*t)*pow(r,2) + 
    2*pow(dt,3)*exp(r*t)*y_n*pow(r,4) - 
    5*f[1]*exp(r*t_n)*((1 + sqrt(5))*pow(dt,2)*pow(r,2)*(1 + r*t) - 
    (1 + 3*sqrt(5))*dt*r*(2 + r*t*(2 + r*t)) + 
    2*sqrt(5)*(6 + r*t*(6 + r*t*(3 + r*t)))) + 
    exp(r*t_n)*(-2*f[3]*(30 + r*(pow(dt,2)*r*(1 + r*t) - 
    5*dt*(2 + r*t*(2 + r*t)) + 5*t*(6 + r*t*(3 + r*t)))) + 
    5*f[2]*((-1 + sqrt(5))*pow(dt,2)*pow(r,2)*(1 + r*t) - 
    (-1 + 3*sqrt(5))*dt*r*(2 + r*t*(2 + r*t)) + 
    2*sqrt(5)*(6 + r*t*(6 + r*t*(3 + r*t))))) - 
    60*sqrt(5)*f[2]*exp(r*t)*r*t_n + 60*f[3]*exp(r*t)*r*t_n - 
    10*f[2]*dt*exp(r*t)*pow(r,2)*t_n + 
    30*sqrt(5)*f[2]*dt*exp(r*t)*pow(r,2)*t_n - 
    20*f[3]*dt*exp(r*t)*pow(r,2)*t_n + 
    5*f[2]*pow(dt,2)*exp(r*t)*pow(r,3)*t_n - 
    5*sqrt(5)*f[2]*pow(dt,2)*exp(r*t)*pow(r,3)*t_n + 
    2*f[3]*pow(dt,2)*exp(r*t)*pow(r,3)*t_n - 
    30*sqrt(5)*f[2]*exp(r*t)*pow(r,2)*pow(t_n,2) + 
    30*f[3]*exp(r*t)*pow(r,2)*pow(t_n,2) - 
    5*f[2]*dt*exp(r*t)*pow(r,3)*pow(t_n,2) + 
    15*sqrt(5)*f[2]*dt*exp(r*t)*pow(r,3)*pow(t_n,2) - 
    10*f[3]*dt*exp(r*t)*pow(r,3)*pow(t_n,2) - 
    10*sqrt(5)*f[2]*exp(r*t)*pow(r,3)*pow(t_n,3) + 
    10*f[3]*exp(r*t)*pow(r,3)*pow(t_n,3) + 
    5*f[1]*exp(r*t)*((1 + sqrt(5))*pow(dt,2)*pow(r,2)*(1 + r*t_n) - 
    (1 + 3*sqrt(5))*dt*r*(2 + r*t_n*(2 + r*t_n)) + 
    2*sqrt(5)*(6 + r*t_n*(6 + r*t_n*(3 + r*t_n)))) + 
    2*f[0]*(exp(r*t_n)*(30 + r*(-(pow(dt,3)*pow(r,2)) + 
    6*pow(dt,2)*r*(1 + r*t) - 10*dt*(2 + r*t*(2 + r*t)) + 
    5*t*(6 + r*t*(3 + r*t)))) + 
    exp(r*t)*(-30 + r*(pow(dt,3)*pow(r,2) - 6*pow(dt,2)*r*(1 + r*t_n) + 
    10*dt*(2 + r*t_n*(2 + r*t_n)) - 5*t_n*(6 + r*t_n*(3 + r*t_n))))))/(2.0*pow(dt,3)*exp(r*t_n)*pow(r,4))

def interp(f, dt, t):
    return (2*pow(dt,3)*f[0] + 10*pow(t,3)*(-f[0] + sqrt(5)*f[1] - sqrt(5)*f[2] + f[3]) + 
    pow(dt,2)*t*(-12*f[0] + 5*(1 + sqrt(5))*f[1] + 5*f[2] - 5*sqrt(5)*f[2] + 2*f[3]) - 
    5*dt*pow(t,2)*(-4*f[0] + f[1] + 3*sqrt(5)*f[1] + f[2] - 3*sqrt(5)*f[2] + 2*f[3]))/(2.*pow(dt,3))

def rk4(f, r, t_n, y_n, dt, dtp):
    M = 10
    h = dtp/float(M)
    
    y = y_n
    i = 0
    
    for t in np.linspace(t_n, t_n+dtp-h, M):
        k1 = h*(interp(f, dt, t)       + r*y)
        k2 = h*(interp(f, dt, t+0.5*h) + r*(y + 0.5*k1))
        k3 = h*(interp(f, dt, t+0.5*h) + r*(y + 0.5*k2))
        k4 = h*(interp(f, dt, t+h)     + r*(y + k3))
        
        y_n = y
        y   = y + (k1 + 2*k2 + 2*k3 + k4)/6.0;
        
        i += 0.5*(y+y_n)*h
    
    return (y, i)

def solve_it(a, d, r, dt, max_iter=10, plot_it=False, exact=False):
    #number of timesteps to run
    #N = int(2/dt)
    N = 1
    
    C = a + d + r
    
    y = np.zeros(N+1)
    t = np.linspace(0,N*dt,num=N+1)
    
    # initial condition
    y[0] = 1.0
    
    for n in range(N):
        I_R = r*y[n]*dt*np.diff(quad_pts)
        
        # stupid predictor
        y_prev = np.array([y[n]]*N_A)
        y_curr = np.array([y[n]]*N_A)
        
        # less stupid predictor
        
#        for m in range(0, N_A-1):
#            dtp = dt*(quad_pts[m+1] - quad_pts[m])
#            y_AD = (y_prev[m] + dtp*a*y[n])/(1 - d*dtp)
#            y_prev[m+1] = (y_prev[m] + dtp*(a*y_prev[m] + d*y_AD))/(1 - r*dtp)
#        
#        I_R[0] = r*int_4(y_prev, 0, dt)
#        I_R[1] = r*int_4(y_prev, 1, dt)
#        I_R[2] = r*int_4(y_prev, 2, dt)
        
        # corrector iterations
        for k in range(max_iter):
            y_AD = y[n]
            for m in range(0,N_A-1):
                dtp = dt*(quad_pts[m+1] - quad_pts[m])
                
                I_R[m] = r*int_4(y_prev, m, dt)
                
                y_AD = (y_AD + dtp*(
                          a*y_curr[m] - a*y_prev[m] 
                        - d*y_prev[m+1])
                        + (a+d)*int_4(y_prev, m, dt)
                        + I_R[m])/(1 - dtp*d)
                
                f_const = d*y_AD - d*y_prev[m+1] + a*y_curr[m] - a*y_prev[m]
                f = f_const + (a+d)*y_prev
                
                if exact:
                    y_curr[m+1] = ode_solve_pt(f, r,
                                               quad_pts[m]*dt, y_curr[m],
                                               dt, quad_pts[m+1]*dt)
                    I_R[m] = y_curr[m+1] - y_curr[m] - int_4(f, m, dt)
                                               
                    #I_R[m] = ((y_curr[m+1] - y_curr[m])
                    #         - f_const*dtp
                    #         - (a+d)*int_4(y_prev, m, dt))
                else:
                    soln = rk4(f, r, quad_pts[m]*dt, y_curr[m], dt, dtp)
                    y_curr[m+1] = soln[0]
                    I_R[m] = r*soln[1]
            
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
    
    print 'L^1 error: ', l1err/l1exact, ' dt: ', dt
    
    return l1err/l1exact
    
    #return l1err/l1exact
