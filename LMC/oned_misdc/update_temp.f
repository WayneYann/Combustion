      subroutine update_temp(nx,scal_old,scal_new,
     $                       aofs,alpha,beta,Rhs,dx,dt,be_cn_theta)
      implicit none
      include 'spec.h'
      integer nx
      real*8 scal_old(-1:nx  ,nscal)
      real*8 scal_new(-1:nx  ,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8    alpha(0 :nx-1)
      real*8     beta(-1:nx  ,nscal)
      real*8      Rhs(0 :nx-1)
      real*8 dx
      real*8 dt
      real*8 be_cn_theta
      
      real*8  dth,dxsqinv
      real*8  Ymid(Nspec), rho_old, rho_new, rho_mid, cpmix
      real*8  beta_lo,beta_hi
      real*8  visc_term, RWRK
      integer i,n,is, IWRK

      dth = be_cn_theta * dt
      dxsqinv = 1.d0/(dx*dx)
      
c*************************************************************************
c       Create RHS = time n diffusive term.
c*************************************************************************

      do i = 0,nx-1
         
         if (coef_avg_harm .eq. 1) then
            beta_lo = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i-1,Temp))
            beta_hi = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i+1,Temp))
         else
            beta_lo = 0.5*(beta(i,Temp) + beta(i-1,Temp))
            beta_hi = 0.5*(beta(i,Temp) + beta(i+1,Temp))
         endif

         visc_term =  
     $        beta_hi*(scal_old(i+1,Temp)-scal_old(i  ,Temp)) -
     $        beta_lo*(scal_old(i  ,Temp)-scal_old(i-1,Temp))
         Rhs(i) = visc_term * dxsqinv * dth

      enddo

c*************************************************************************
c     Construct alpha from half-time level data
c     (here, build rho to ensure correct extract of Y regardless of whether
c     rho was updated independently - ie, assume no reset_rho_in_rho_states))
c*************************************************************************

      do i = 0,nx-1
         rho_old = 0.d0
         rho_new = 0.d0
         do n = 1,Nspec
            is = FirstSpec + n - 1
            rho_old = rho_old + scal_old(i,is)
            rho_new = rho_new + scal_new(i,is)
         enddo
         do n = 1,Nspec
            is = FirstSpec + n - 1
            Ymid(n) =
     &           0.5d0*(scal_old(i,is)/rho_old + scal_new(i,is)/rho_new)
         enddo
         call CKCPBS(scal_old(i,Temp),Ymid,IWRK,RWRK,cpmix)

c     The correct rho_new to use in alpha is the one in the new state.  Depending 
c     on where we are in the algorithm, this may or may not agree with rho_new above
         alpha(i) = 0.5d0 * (rho_old + scal_new(i,Density)) * cpmix

      enddo


c*************************************************************************
c     Set Rhs.  Also, put advection terms in new state.
c*************************************************************************

      do i = 0,nx-1
         scal_new(i,Temp) = scal_old(i,Temp) + dt*aofs(i,Temp)
         Rhs(i) = alpha(i) * scal_new(i,Temp)
      enddo
      
      end
      
