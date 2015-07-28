from __future__ import division
import numpy as np
import sys
import pdb
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
from scipy.interpolate import interp1d
from math import *

def solve_pde(a, eps, r, Nx=300, dt=None, FT=0.05,
              do_linear=False, do_modified=True, N_A=4, max_iter=2):
   endpt = 20.0
   Nx = Nx+1
   # spacing
   h = endpt/(Nx-1)
   
   # timestep
   if dt==None:
      # if not specified by a parameter, use the default (to ensure stability)
      dt = h/5.0
   # number of timesteps
   Nt = int(ceil(FT/dt))
   
   # this is unfortunately hard coded according to the number of substeps...
   if N_A == 4:
       # Gauss-Lobatto quadrature points
       quad_pts = [0, (1-sqrt(1/5.0))*0.5, (1+sqrt(1/5.0))*0.5, 1]
       # quadrature weights computed using Mathematica
       weights = [[(1.0/120.0)*(11 + sqrt(5)), (1.0/120.0)*(25 - sqrt(5)), (1.0/120.0)*(25 - 13*sqrt(5)), (1.0/120.0)*(-1+sqrt(5))],
                  [-1.0/(12*sqrt(5)), 7.0/(12*sqrt(5)), 7.0/(12*sqrt(5)), -1.0/(12*sqrt(5))],
                  [(1.0/120.0)*(-1 + sqrt(5)), (1.0/120.0)*(25 - 13*sqrt(5)), (1.0/120.0)*(25 - sqrt(5)), (1.0/120.0)*(11 + sqrt(5))]]
   elif N_A == 2:
       quad_pts = [0, 1]
   
   # perform our Gaussian sub-interval integration
   def integ4(f, m, l):
       return np.apply_along_axis(lambda f_j: np.dot(f_j, weights[m])*l, 0, f)
   
   # compute the integral using trapezoid rule
   def integ2(f, m, l):
       return 0.5*(f[0]+f[1])*l
   
   # return the integral (depending on how many substeps used)
   def integ(f, m, l):
       if N_A==4:
           return integ4(f, m, l)
       elif N_A==2:
           return integ2(f, m, l)

   # calculate the value of the interpolating cubic at an arbitrary time t
   def interp_4(f, t):
       return (2*pow(dt,3)*f[0] + 10*pow(t,3)*(-f[0] + sqrt(5)*f[1] - sqrt(5)*f[2] + f[3]) + 
       pow(dt,2)*t*(-12*f[0] + 5*(1 + sqrt(5))*f[1] + 5*f[2] - 5*sqrt(5)*f[2] + 2*f[3]) - 
       5*dt*pow(t,2)*(-4*f[0] + f[1] + 3*sqrt(5)*f[1] + f[2] - 3*sqrt(5)*f[2] + 2*f[3]))/(2.*pow(dt,3))

   # return the value of the linear interpolant
   def interp_linear(f, t):
       return (1-t/dt)*f[0] + (t/dt)*f[1]
   
   # return the value of the interpolating polynomial at time t
   def interp(f, t):
       if N_A == 4:
           if do_linear:
               return interp_4(f, t)
           else:
               return f
       elif N_A == 2:
           if do_linear:
               return interp_linear(f, t)
           else:
               return f

   # compute the right-hand side used in the Runge-Kutta integrator
   def rhs(f, z, t, t_n):
       if do_modified:
           return (interp(f,t) + FR(z))/(1-(t-t_n)*FRp(z))
       else:
           return interp(f,t) + FR(z)

   # use 4th order Runge-Kutta to compute the solution to the correction ODE
   def rk4(f, t_n, y_n, dtp):
       M = 1
       h_rk = dtp/float(M)
       
       z = y_n
       
       for t in np.linspace(t_n, t_n+dtp-h_rk, M):
           k1 = h_rk*rhs(f, z,        t,          t_n)
           k2 = h_rk*rhs(f, z+0.5*k1, t+0.5*h_rk, t_n)
           k3 = h_rk*rhs(f, z+0.5*k2, t+0.5*h_rk, t_n)
           k4 = h_rk*rhs(f, z+k3,     t+h_rk,     t_n)
           
           z  = z + (k1 + 2*k2 + 2*k3 + k4)/6.0;
       
       return z

   # return the reaction term
   def FR(z):
       return r*z*(z-1)*(z-0.5)
       
   # return the first derivative of the reaction term
   def FRp(z):
       return 0.5*r*(1-6*z+6*z**2)

   # return the advection term
   def FA(z):
       return np.dot(A, z) + bca

   # return the diffusion term
   def FD(z):
       return np.dot(D, z) + bcd
   
   def AD_RHS(z):
       return np.apply_along_axis(lambda z_m: FA(z_m) + FD(z_m), 1, z)
   
   def R_RHS(z):
       return np.apply_along_axis(lambda z_m: FR(z_m), 1, z)
   
   # compute the integral of the advection and diffusion terms for a subinterval
   def I_AD(z, m):    
       return integ(AD_RHS(z), m, dt)

   # compute the integral of the reaction term over a subinterval
   def I_R(z, m):
       return integ(R_RHS(z), m, dt)
   
   # advance the solution by one timestep
   def advance(n):
       # boundary values
       y[n+1][0] = y[n][0]
       y[n+1][-1] = y[n][-1]
           
       # stupid predictor
       y_prev = np.array([y[n][1:-1]]*N_A)
       y_curr = np.copy(y_prev)
       
       # corrector iterations
       for k in range(max_iter):
           y_AD = y[n][1:-1]
           for m in range(0,N_A-1):
               dtp = dt*(quad_pts[m+1] - quad_pts[m])
               
               # compute the intermediate solution using explicit advection
               # and implicit diffusion
               y_AD = spsolve(Dinv[m],
                          y_AD + dtp*(FA(y_curr[m]) - FA(y_prev[m])
                                    + bcd           - FD(y_prev[m+1]))
                               + I_AD(y_prev, m)
                               + I_R(y_prev, m))
               
               # compute the constant forcing term in the correction ODE
               f_const = (FA(y_curr[m]) - FA(y_prev[m])
                        + FD(y_AD)      - FD(y_prev[m+1]) 
                                        - FR(y_prev[m+1]))
               
               # add the interpolating polynomial given by the quadrature
               f = f_const + AD_RHS(y_prev) + R_RHS(y_prev)
               
               if not do_linear:
                  f = integ(f,m,dt)/dtp
               
               # compute the solution to the correction ODE
               soln = rk4(f, quad_pts[m]*dt, y_curr[m], dtp)
               y_curr[m+1] = soln
           
           # move on to the next MISDC iteration...
           y_prev = y_curr
       y[n+1][1:-1] = y_curr[N_A - 1]
   
   x = np.linspace(0, endpt, num=Nx)
   y = np.zeros((Nt+2, Nx))

   y[0] = 0.5*(np.tanh(10-2*x)+1)

   # identity
   I = np.eye(Nx-2)
   # create the differentiation matrices
   
   fourth_order = True
   
   bca = np.zeros(Nx-2)
   if fourth_order:
      Dx = (np.diag([1.0/12.0]*(Nx-4),-2)
          + np.diag([-2.0/3.0]*(Nx-3),-1)
          + np.diag([2.0/3.0]*(Nx-3),  1)
          + np.diag([-1.0/12.0]*(Nx-4)   ,2))/h
      Dx[0] = np.array([-13.0/12.0, 2.0, -1.0, 1.0/3.0, -1.0/20.0] + [0]*(Nx-7))/h
      Dx[-1] = np.array([0]*(Nx-7) + [1.0/20.0, -1.0/3.0, 1.0, -2.0, 13.0/12.0])/h
      bca[0] = (-1.0/5.0)*a/h
      bca[1] = (1.0/12.0)*a/h
   else:
      Dx = (np.diag([1]*(Nx-3),1) - np.diag([1]*(Nx-3),-1))/(2*h)
      bca[0] = -0.5*a/h
   
   A = a*Dx
   
   bcd = np.zeros(Nx-2)
   # laplacian operator
   if fourth_order:
      Dxx = (I*(-5.0/2.0) + np.diag([-1.0/12.0]*(Nx-4),-2)
                          + np.diag([4.0/3.0]*(Nx-3),  -1)
                          + np.diag([4.0/3.0]*(Nx-3),   1)
                          + np.diag([-1.0/12.0]*(Nx-4), 2))/h**2
      Dxx[0]  = np.array([-5.0/3.0, 1.0/2.0, 1.0/3.0, -1.0/12.0] + [0]*(Nx-6))/h**2
      Dxx[-1] = np.array([0]*(Nx-6) + [-1.0/12.0, 1.0/3.0, 1.0/2.0, -5.0/3.0])/h**2
      bcd[0] = (11.0/12.0)*eps/(h**2)
      bcd[1] =  (-1.0/12.0)*eps/(h**2)
   else:
      Dxx = (-I*2 + np.diag([1]*(Nx-3),1) + np.diag([1]*(Nx-3),-1))/h**2
      bcd[0] = eps/h**2
   
   D = eps*Dxx

   # we treat diffusion implicitly
   # for efficiency, we create the sparse matrices now
   # which we later use in the sparse solve
   # (I - dtp D) y = RHS
   Dinv = np.empty(N_A - 1, dtype=object)
   for m in range(N_A - 1):
       dtp = dt*(quad_pts[m+1] - quad_pts[m])
       Dinv[m] = csc_matrix(I - dtp*D)
   
   plot_it = True
   if plot_it:
      plt.ion()
   step = 10
   for n in range(Nt):
       advance(n)
#       if (n+1)%(int(1+Nt/10)) == 0:
#         print 't = ', dt*(n+1)
       if ((n+1)%step == 0 or n == Nt-1) and plot_it:
          print 't = ', dt*(n+1)
          plt.cla()
          plt.plot(x,y[n+1])
          plt.ylim((0,1.1))
          plt.draw()
    
#    np.savetxt('ft', np.array([x, y[Nt]]).T)
   
   #plt.plot(x,y[Nt])
   #plt.ylim((-0.1, 1.1))
   #plt.show(block=True)
   
   print 'solved. final time=',dt*Nt
   print 'Nx=',Nx,' dt=',dt
   return (x, y[Nt])
