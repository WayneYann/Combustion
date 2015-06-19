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

def ode_solve(f, r, c, y_n, dt, t):
    return (60*f[0] - 60*sqrt(5)*f[1] + 60*sqrt(5)*f[2] - 60*f[3] - 60*f[0]*exp(r*t) + 
     60*sqrt(5)*f[1]*exp(r*t) - 60*sqrt(5)*f[2]*exp(r*t) + 60*f[3]*exp(r*t) - 
     40*f[0]*dt*r + 10*f[1]*dt*r + 30*sqrt(5)*f[1]*dt*r + 10*f[2]*dt*r - 
     30*sqrt(5)*f[2]*dt*r + 20*f[3]*dt*r + 40*f[0]*dt*exp(r*t)*r - 
     10*f[1]*dt*exp(r*t)*r - 30*sqrt(5)*f[1]*dt*exp(r*t)*r - 
     10*f[2]*dt*exp(r*t)*r + 30*sqrt(5)*f[2]*dt*exp(r*t)*r - 
     20*f[3]*dt*exp(r*t)*r + 12*f[0]*pow(dt,2)*pow(r,2) - 
     5*f[1]*pow(dt,2)*pow(r,2) - 5*sqrt(5)*f[1]*pow(dt,2)*pow(r,2) - 
     5*f[2]*pow(dt,2)*pow(r,2) + 5*sqrt(5)*f[2]*pow(dt,2)*pow(r,2) - 
     2*f[3]*pow(dt,2)*pow(r,2) - 12*f[0]*pow(dt,2)*exp(r*t)*pow(r,2) + 
     5*f[1]*pow(dt,2)*exp(r*t)*pow(r,2) + 
     5*sqrt(5)*f[1]*pow(dt,2)*exp(r*t)*pow(r,2) + 
     5*f[2]*pow(dt,2)*exp(r*t)*pow(r,2) - 
     5*sqrt(5)*f[2]*pow(dt,2)*exp(r*t)*pow(r,2) + 
     2*f[3]*pow(dt,2)*exp(r*t)*pow(r,2) - 2*f[0]*pow(dt,3)*pow(r,3) - 
     2*c*pow(dt,3)*pow(r,3) + 2*f[0]*pow(dt,3)*exp(r*t)*pow(r,3) + 
     2*c*pow(dt,3)*exp(r*t)*pow(r,3) + 
     2*pow(dt,3)*exp(r*t)*y_n*pow(r,4) + 60*f[0]*r*t - 60*sqrt(5)*f[1]*r*t + 
     60*sqrt(5)*f[2]*r*t - 60*f[3]*r*t - 40*f[0]*dt*pow(r,2)*t + 10*f[1]*dt*pow(r,2)*t + 
     30*sqrt(5)*f[1]*dt*pow(r,2)*t + 10*f[2]*dt*pow(r,2)*t - 
     30*sqrt(5)*f[2]*dt*pow(r,2)*t + 20*f[3]*dt*pow(r,2)*t + 
     12*f[0]*pow(dt,2)*pow(r,3)*t - 5*f[1]*pow(dt,2)*pow(r,3)*t - 
     5*sqrt(5)*f[1]*pow(dt,2)*pow(r,3)*t - 5*f[2]*pow(dt,2)*pow(r,3)*t + 
     5*sqrt(5)*f[2]*pow(dt,2)*pow(r,3)*t - 2*f[3]*pow(dt,2)*pow(r,3)*t + 
     30*f[0]*pow(r,2)*pow(t,2) - 30*sqrt(5)*f[1]*pow(r,2)*pow(t,2) + 
     30*sqrt(5)*f[2]*pow(r,2)*pow(t,2) - 30*f[3]*pow(r,2)*pow(t,2) - 
     20*f[0]*dt*pow(r,3)*pow(t,2) + 5*f[1]*dt*pow(r,3)*pow(t,2) + 
     15*sqrt(5)*f[1]*dt*pow(r,3)*pow(t,2) + 5*f[2]*dt*pow(r,3)*pow(t,2) - 
     15*sqrt(5)*f[2]*dt*pow(r,3)*pow(t,2) + 10*f[3]*dt*pow(r,3)*pow(t,2) + 
     10*f[0]*pow(r,3)*pow(t,3) - 10*sqrt(5)*f[1]*pow(r,3)*pow(t,3) + 
     10*sqrt(5)*f[2]*pow(r,3)*pow(t,3) - 10*f[3]*pow(r,3)*pow(t,3))/(2.*pow(dt,3)*pow(r,4))

def fixed_pt(dt):
    return (25*(-45 + 3*(-25 + 7*sqrt(5))*dt + 16*(3 + 5*sqrt(5))*pow(dt,2) + 
    64*sqrt(5)*pow(dt,3) + 256*(-3 + sqrt(5))*pow(dt,4))*exp(4*dt) + 
    (15 + 40*dt + 32*pow(dt,2))*
    (-15*sqrt(5) + 6*(5 + 2*sqrt(5))*dt + 16*sqrt(5)*pow(dt,2))*
    exp((4*dt)/sqrt(5)) + 5*(-225 + 30*(-24 + sqrt(5))*dt + 
    24*(-5 + 7*sqrt(5))*pow(dt,2) + 32*(75 + 23*sqrt(5))*pow(dt,3) + 
    256*(5 + sqrt(5))*pow(dt,4))*exp((2 + 6/sqrt(5))*dt) + 
    5*(45*(5 + sqrt(5)) - 135*(-5 + sqrt(5))*dt - 
    8*(-15 + 77*sqrt(5))*pow(dt,2) - 64*(29 + 2*sqrt(5))*pow(dt,3) + 
    512*(-3 + sqrt(5))*pow(dt,4))*exp((2*(5 + sqrt(5))*dt)/5.) + 
    (225*(5 + sqrt(5)) + 330*(5 + 3*sqrt(5))*dt + 400*(-9 + sqrt(5))*pow(dt,2) - 
    64*(145 + 207*sqrt(5))*pow(dt,3) + 20480*pow(dt,4) + 
    4096*sqrt(5)*pow(dt,5))*exp((4*(5 + sqrt(5))*dt)/5.) + 
    5*dt*(15*(-7 + 3*sqrt(5)) + 8*(-45 + 19*sqrt(5))*dt + 
    64*(-5 + 2*sqrt(5))*pow(dt,2))*exp((6*(5 + sqrt(5))*dt)/5.) + 
    5*(-45*sqrt(5) - 30*(-5 + 6*sqrt(5))*dt + (840 - 344*sqrt(5))*pow(dt,2) + 
    224*(9 + 5*sqrt(5))*pow(dt,3) - 256*(-3 + 13*sqrt(5))*pow(dt,4) + 
    2048*(-1 + sqrt(5))*pow(dt,5))*exp(6*dt + (2*dt)/sqrt(5)) - 
    (-15 + 7*sqrt(5))*dt*(-15 - 80*dt - 192*pow(dt,2) + 256*pow(dt,3))*
    exp(8*dt + (4*dt)/sqrt(5)))/((-5*(3*sqrt(5) - (3 + 7*sqrt(5))*dt + 8*(-1 + sqrt(5))*pow(dt,2))*
    exp(2*dt) + (15*sqrt(5) - 6*(5 + 2*sqrt(5))*dt - 16*sqrt(5)*pow(dt,2))*
    exp((2*dt)/sqrt(5)) + (-15 + 7*sqrt(5))*dt*exp(4*dt + (2*dt)/sqrt(5))
    )*(5*(3*sqrt(5) - 4*(1 + 3*sqrt(5))*dt + 16*(1 + sqrt(5))*pow(dt,2))*
    exp(2*dt) - (15 - 80*dt + 192*pow(dt,2) + 256*pow(dt,3))*
    exp((2*dt)/sqrt(5)) - 5*
    (3*sqrt(5) + (4 - 12*sqrt(5))*dt + 16*(-1 + sqrt(5))*pow(dt,2))*
    exp((2 + 4/sqrt(5))*dt) + 
    (15 - 40*dt + 32*pow(dt,2))*exp(4*dt + (2*dt)/sqrt(5))))

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
        I_R = np.multiply(r*y[n]*dt,np.diff(quad_pts))
        
        # stupid predictor
        y_prev = [y[n]]*N_A
        y_curr = [y[n]]*N_A
        
        # less stupid predictor
        
        for m in range(0, N_A-1):
            dtp = dt*(quad_pts[m+1] - quad_pts[m])
            y_AD = (y_prev[m] + dtp*a*y[n])/(1 - d*dtp)
            y_prev[m+1] = (y_prev[m] + dtp*(a*y_prev[m] + d*y_AD))/(1 - r*dtp)
        
        I_R[0] = r*int_4(y_prev, 0, dt)
        I_R[1] = r*int_4(y_prev, 1, dt)
        I_R[2] = r*int_4(y_prev, 2, dt)
        
        
        # corrector iterations
        for k in range(max_iter):
            y_AD = y[n]
            for m in range(0,N_A-1):
                dtp = dt*(quad_pts[m+1] - quad_pts[m])
                
                #I_R[m] = r*int_4(y_prev, m, dt)
                
                y_AD = (y_AD + dtp*(
                          a*y_curr[m] - a*y_prev[m] 
                        - d*y_prev[m+1])
                        + (a+d)*int_4(y_prev, m, dt)
                        + I_R[m])/(1 - dtp*d)
                
                f_const = d*y_AD - d*y_prev[m+1] + a*y_curr[m] - a*y_prev[m]
                
                # compute the forcing term
                #D = y_prev[m+1]*d
                #A_n = a*y_prev[m]
                #D_n = d*y_prev[m]
                #A = a*y_prev[m+1]
                
                #f_lin_prev = A_n + D_n
                #f_lin_next = A + D
                #c1 = (f_lin_next - f_lin_prev)/dtp
                #c0 = f_const + f_lin_prev
                
                # solve the ODE here
                # solving y_next' = r*y_next + c + mt
                #y_curr[m+1] = ((exp(r*dtp)*(c1+c0*r) - c1 - c0*r - c1*r*dtp)/r**2
                #               + exp(r*dtp)*y_prev[m])
                
                
                #print 'miscd iteration: ', k, ' residual: ', f_const
                
                y_curr[m+1] = ode_solve(np.multiply((a+d),y_prev), 
                                        r, f_const, y_curr[m],
                                        dt, dtp)
                
                #y_curr[m+1] = (y_curr[m] + dtp*(#a*y_curr[m] - a*y_prev[m]
                #               + d*y_AD
                #               -  d*y_prev[m+1] - r*y_prev[m+1])
                #               + (a+d+r)*int_4(y_prev, m, dt))/(1 - r*dtp)
                               # + I_R[m])/(1 - r*dtp)
                
                I_R[m] = ((y_curr[m+1] - y_curr[m])
                             - f_const*dtp
                             - (a+d)*int_4(y_prev, m, dt))
               
            # move on to the next MISDC iteration...
            y_prev = y_curr
        
        print 'dt: ', dt, ' fixed point soln: ', fixed_pt(dt)
        #print 'difference: ', abs(fixed_pt(dt) - y_prev[N_A - 1])
        
        y[n+1] = y_prev[N_A - 1]
        #print 'y_prev final: ', y_prev
    
    exact = np.exp(C*t)
    l1exact = np.sum(exact)
    l1err = np.sum(np.abs(y-exact))
    
    if plot_it:
        print 'relative L^1 error: ', l1err/l1exact
        plt.plot(t,y,t,exact)
        plt.show()
    
    print 'L^1 error: ', l1err/l1exact
    return l1err/l1exact
