      subroutine update_vel(nx,vel_old,vel_new,gp,rhohalf,
     &                      macvel,sedge,alpha,beta,
     &                      Rhs,dx,dt,be_cn_theta)
      implicit none
      include 'spec.h'
      integer nx
      real*8 vel_old(-1:nx  )
      real*8 vel_new(-1:nx  )
      real*8      gp(0 :nx-1)
      real*8   alpha(0 :nx-1)
      real*8    beta(-1:nx  )
      real*8     Rhs(0 :nx-1)
      real*8  rhohalf(0 :nx-1)
      real*8   macvel(0 :nx)
      real*8    sedge(0 :nx)
      real*8 dx
      real*8 dt
      real*8 be_cn_theta
      
      real*8  dth,dxsqinv,aofs
      real*8  beta_lo,beta_hi
      real*8  visc_term, RWRK
      integer i,n,is, IWRK

      dth = be_cn_theta * dt * 2
      dxsqinv = 1.d0/(dx*dx)

c     rho.DU/Dt + G(pi) = D(tau), here D(tau) = 4/3 d/dx ( mu . du/dx )
      do i = 0,nx-1
         
         if (coef_avg_harm .eq. 1) then
            beta_lo = 2.d0 / (1.d0/beta(i)+1.d0/beta(i-1))
            beta_hi = 2.d0 / (1.d0/beta(i)+1.d0/beta(i+1))
         else
            beta_lo = 0.5*(beta(i) + beta(i-1))
            beta_hi = 0.5*(beta(i) + beta(i+1))
         endif

         visc_term =  
     $        ( beta_hi*(vel_old(i+1)-vel_old(i  ))
     $        - beta_lo*(vel_old(i  )-vel_old(i-1)) )*dxsqinv

         aofs = ( (macvel(i+1)*sedge(i+1) - macvel(i)*sedge(i))  -
     $        (macvel(i+1)-macvel(i)) * 
     $        0.5d0 * (sedge(i)+sedge(i+1)) ) /dx

         vel_new(i) = vel_old(i) + dt * ( -aofs + gp(i)/rhohalf(i) )

         Rhs(i) = vel_new(i) 
     &        + (4.d0/3.d0)*be_cn_theta*visc_term/rhohalf(i)

         alpha(i) = 1.d0
      enddo

      end
      
