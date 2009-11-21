      subroutine update_temp(nx,scal_old,scal_new,aofs,
     &                       alpha,beta_old,beta_new,I_R,Rhs,
     &                       dx,dt,be_cn_theta,time)
      implicit none
      include 'spec.h'
      integer nx
      real*8 scal_old(-1:nx  ,nscal)
      real*8 scal_new(-1:nx  ,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8    alpha(0 :nx-1)
      real*8 beta_old(-1:nx  ,nscal)
      real*8 beta_new(-1:nx  ,nscal)
      real*8      Rhs(0 :nx-1)
      real*8      I_R(0 :nx-1)
      real*8 dx, dt, time
      real*8 be_cn_theta
      
      real*8  dxsqinv
      real*8  Ymid(Nspec), rho_old, rho_new, rho_mid, cpmix
      real*8  beta_lo,beta_hi
      real*8  visc(0:nx-1)
      real*8  RWRK
      integer i,n,is, IWRK

c*************************************************************************
c     Add adv terms to old state prior to doing stuff below
c*************************************************************************
      do i = 0,nx-1
         scal_new(i,Temp) = scal_old(i,Temp) + dt*aofs(i,Temp)
      enddo

c*************************************************************************
c     Initialize RHS = (1-theta) * Div( lambda Grad(T) ) at old time
c*************************************************************************      
      call set_bc_grow_s(nx,scal_old,dx,time)
      dxsqinv = 1.d0/(dx*dx)
      do i = 0,nx-1
         Rhs(i) = I_R(i)
         if (coef_avg_harm .eq. 1) then
            beta_lo = 2.d0 /
     &           (1.d0/beta_old(i,Temp)+1.d0/beta_old(i-1,Temp))
            beta_hi = 2.d0 /
     &           (1.d0/beta_old(i,Temp)+1.d0/beta_old(i+1,Temp))
         else
            beta_lo = 0.5*(beta_old(i,Temp)+beta_old(i-1,Temp))
            beta_hi = 0.5*(beta_old(i,Temp)+beta_old(i+1,Temp))
         endif
         
         Rhs(i) =  Rhs(i) + (1.d0 - be_cn_theta)*dt*dxsqinv*(
     $        beta_hi*(scal_old(i+1,Temp)-scal_old(i  ,Temp)) -
     $        beta_lo*(scal_old(i  ,Temp)-scal_old(i-1,Temp)) )
      enddo

c*************************************************************************
c     Add rho.D.Grad(Y).Grad(H)  [with appropriate theta time centering]
c*************************************************************************      
      call rhoDgradHgradY(nx,scal_old,beta_old,visc,dx,time)
      do i = 0,nx-1
         Rhs(i) = Rhs(i) + (1.d0 - be_cn_theta)*visc(i) 
      enddo
      call rhoDgradHgradY(nx,scal_new,beta_new,visc,dx,time+dt)
      do i = 0,nx-1
         Rhs(i) = Rhs(i) + be_cn_theta*visc(i) 
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
         alpha(i) = 0.5d0 * (rho_old + scal_new(i,Density)) * cpmix
         Rhs(i) = Rhs(i) + scal_new(i,Temp)*alpha(i)
      enddo

      end
      
