from __future__ import division
import numpy as np
import sys
import pdb
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
from scipy.interpolate import interp1d
from math import *

def solve_pde(a, eps, r, Nx=300, dt=None, FT=0.05, max_iter=2):
   endpt = 20.0
   Nx = Nx+1
   # spacing
   h = endpt/(Nx-1)
   
   # timestep
   if dt==None:
      # if not specified by a parameter, use the default (to ensure stability)
      dt = h/5.0
   
   dt1 = dt/2.0
   dt2 = dt/2.0
   
   # number of timesteps
   Nt = int(ceil(FT/dt))

   # compute the right-hand side used in the Runge-Kutta integrator
   def rhs(f, z, t, t_n):
        return (f + FR(z))/(1-(t-t_n)*FRp(z))

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
   
   def int_1(f1, f2, f3):
      return (5*f1/24.0 + f2/3.0 - f3/24.0)*dt
   def int_2(f1, f2, f3):
      return (-f1/24.0 + f2/3.0 + 5*f3/24.0)*dt
   
   # advance the solution by one timestep
   def advance(n):
       # boundary values
       y[n+1][0] = y[n][0]
       y[n+1][-1] = y[n][-1]
           
       # stupid predictor
       yn      = np.copy(y[n][1:-1])         
       y_prev1 = np.copy(y[n][1:-1])
       y_prev2 = np.copy(y[n][1:-1])
       y_curr1 = np.copy(y[n][1:-1])
       y_curr2 = np.copy(y[n][1:-1])
       
       F_0 = FA(yn) + FD(yn) + FR(yn)
       # corrector iterations
       for k in range(max_iter):
           y_AD = y[n][1:-1]
           
           F_1 = FA(y_prev1) + FD(y_prev1) + FR(y_prev1)
           F_2 = FA(y_prev2) + FD(y_prev2) + FR(y_prev2)
           
           I_1 = int_1(F_0, F_1, F_2)
           I_2 = int_2(F_0, F_1, F_2)
           
           # compute the intermediate solution using explicit advection
           # and implicit diffusion
           y_AD1 = spsolve(Dinv, yn + dt1*(bcd - FD(y_prev1)) + I_1)
            
           # compute the constant forcing term in the correction ODE
           # and add the term given by the quadrature
           f = (FD(y_AD1) - FD(y_prev1) - FR(y_prev1)) + I_1/dt1
            
           # compute the solution to the correction ODE
           y_curr1 = rk4(f, 0, yn, dt1)
           
           y_AD2 = spsolve(Dinv, y_AD1 + dt2*(FA(y_curr1) - FA(y_prev1)
                                            + bcd         - FD(y_prev2))
                                       + I_2)
           f = (FA(y_curr1) - FA(y_prev1)
              + FD(y_AD2)   - FD(y_prev2) 
              - FR(y_prev2)) + I_2/dt2
           
           y_curr2 = rk4(f, dt1, y_curr1, dt2)
        
           # move on to the next MISDC iteration...
           y_prev1 = y_curr1
           y_prev2 = y_curr2
           
       y[n+1][1:-1] = y_curr2
   
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
   Dinv = csc_matrix(I - dt*0.5*D)
   
#   plt.ion()
   step = 20
   for n in range(Nt):
       advance(n)
       if (n+1)%int(Nt/10 + 1)==0:
          print 't = ', dt*(n+1)
#       if (n+1)%step == 0 or n == Nt-1:
#          print 't = ', dt*(n+1)
#          plt.cla()
#          plt.plot(x,y[n+1])
#          plt.ylim((0,1.1))
#          plt.draw()
    
#    np.savetxt('ft', np.array([x, y[Nt]]).T)
#   
#   plt.plot(x,y[Nt])
#   plt.ylim((-0.1, 1.1))
#   plt.show(block=True)
   
   print 'solved. final time=',dt*Nt
   print 'Nx=',Nx,' dt=',dt
   return (x, y[Nt])
