      subroutine update_rhoh(nx,scal_old,scal_new,aofs,alpha,
     &                       beta,Rhs,dx,dt,be_cn_theta,time)
      implicit none
      include 'spec.h'
      integer nx
      real*8 scal_old(-1:nx  ,nscal)
      real*8 scal_new(-1:nx  ,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8    alpha(0 :nx-1)
      real*8     beta(-1:nx  ,nscal)
      real*8      Rhs(0 :nx-1)
      real*8 dx,dt,be_cn_theta,time
      
      real*8  dxsqinv
      real*8  beta_lo,beta_hi
      real*8  visc_term, RWRK
      integer i,n,is, IWRK

c     FIXME: Add NULN terms
      if (LeEQ1 .ne. 1) then
         print *,'update_rhoh does yet support non-unity Le'
         stop
      endif

      call set_bc_grow_s(nx,scal_old,dx,time)
      dxsqinv = 1.d0/(dx*dx)
      do i = 0,nx-1
         if (coef_avg_harm .eq. 1) then
            beta_lo = 2.d0 / (1.d0/beta(i,RhoH)+1.d0/beta(i-1,RhoH))
            beta_hi = 2.d0 / (1.d0/beta(i,RhoH)+1.d0/beta(i+1,RhoH))
         else
            beta_lo = 0.5*(beta(i,RhoH) + beta(i-1,RhoH))
            beta_hi = 0.5*(beta(i,RhoH) + beta(i+1,RhoH))
         endif

         visc_term = dxsqinv * dt * (1.d0 - be_cn_theta) * (
     $        beta_hi*(scal_old(i+1,RhoH)-scal_old(i  ,RhoH)) -
     $        beta_lo*(scal_old(i  ,RhoH)-scal_old(i-1,RhoH)) )

         scal_new(i,RhoH) = scal_old(i,RhoH) + dt*aofs(i,RhoH)
         Rhs(i) = scal_new(i,RhoH) + visc_term
         alpha(i) = scal_new(i,Density)
      enddo
      
      end
