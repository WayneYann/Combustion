      subroutine get_rhoh_visc_terms(nx,scal,beta,visc,dx)
      implicit none
      include 'spec.h'
      integer nx
      real*8 scal(-1:nx  ,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(0 :nx-1)
      real*8 dx
      
      integer i
      real*8 beta_lo,beta_hi
      real*8 flux_lo,flux_hi
      real*8 dxsqinv, h_lo, h_mid, h_hi
      
      dxsqinv = 1.d0/(dx*dx)

      do i = 0,nx-1
         
c     Compute Div( lambda/cp Grad(hmix) )
         if (coef_avg_harm.eq.1) then
            beta_lo = 2.d0 / (1.d0/beta(i,RhoH)+1.d0/beta(i-1,RhoH))
            beta_hi = 2.d0 / (1.d0/beta(i,RhoH)+1.d0/beta(i+1,RhoH))
         else
            beta_lo = 0.5*(beta(i,RhoH) + beta(i-1,RhoH))
            beta_hi = 0.5*(beta(i,RhoH) + beta(i+1,RhoH))
         endif

         h_hi  = scal(i+1,RhoH) / scal(i+1,Density)
         h_mid = scal(i  ,RhoH) / scal(i  ,Density)
         h_lo  = scal(i-1,RhoH) / scal(i-1,Density)

         flux_hi = beta_hi*(h_hi - h_mid)
         flux_lo = beta_lo*(h_mid - h_lo)
         visc(i) = (flux_hi - flux_lo) * dxsqinv
         
      end do
      end

