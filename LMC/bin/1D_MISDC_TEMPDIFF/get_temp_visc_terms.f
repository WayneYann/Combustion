      subroutine addDivLambdaGradT(scal,beta,visc,dx,lo,hi)
      implicit none
      include 'spec.h'
      real*8 scal(-2:nfine+1,nscal)
      real*8 beta(-1:nfine  ,nscal)
      real*8 visc(-1:nfine)
      real*8 dx
      integer lo,hi
      
      integer i
      real*8 beta_lo,beta_hi
      real*8 flux_lo,flux_hi
      real*8 dxsqinv

      dxsqinv = 1.d0/(dx*dx)
      do i=lo,hi
         if (coef_avg_harm.eq.1) then
            beta_lo = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i-1,Temp))
            beta_hi = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i+1,Temp))
         else
            beta_lo = 0.5*(beta(i,Temp) + beta(i-1,Temp))
            beta_hi = 0.5*(beta(i,Temp) + beta(i+1,Temp))
         endif
         
         flux_hi = beta_hi*(scal(i+1,Temp) - scal(i  ,Temp)) 
         flux_lo = beta_lo*(scal(i  ,Temp) - scal(i-1,Temp)) 
         visc(i) = visc(i) + (flux_hi - flux_lo) * dxsqinv

      enddo

      end
