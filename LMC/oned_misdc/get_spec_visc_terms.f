      subroutine get_spec_visc_terms(scal,beta,visc,dx,time)
      implicit none
      include 'spec.h'
      real*8 scal(-1:nx  ,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(0 :nx-1,*)
      real*8 dx,time
      
      integer i,n,is,IWRK
      real*8 beta_lo,beta_hi
      real*8 flux_lo(maxspec),flux_hi(maxspec)
      real*8 dxsqinv,RWRK
      real*8 Y(-1:nx,maxspec), sum_lo, sum_hi, sumRhoY_lo, sumRhoY_hi
      real*8 RhoYe_lo, RhoYe_hi

      real*8 beta_edge(-1:nx,Nspec)

      call set_bc_grow_s(scal,dx,time)
      do i = -1,nx
         do n=1,Nspec
            Y(i,n) = scal(i,FirstSpec+n-1)/scal(i,Density)
         enddo
      enddo

CCCCCCCCCCC debugging FIXME
C 1006 FORMAT((I5,1X),6(E22.15,1X))      
C      open(UNIT=11, FILE='corr.dat', STATUS = 'REPLACE')
C      write(11,*)'# 256 12'
CCCCCCCCCCCC      

      dxsqinv = 1.d0/(dx*dx)
      do i = 0,nx-1
         sum_lo = 0.d0
         sum_hi = 0.d0
         sumRhoY_lo = 0
         sumRhoY_hi = 0
         do n=1,Nspec
            is = FirstSpec + n - 1
            if (coef_avg_harm.eq.1) then
               beta_lo = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i-1,is))
               beta_hi = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i+1,is))
            else
               beta_lo = 0.5d0*(beta(i,is) + beta(i-1,is))
               beta_hi = 0.5d0*(beta(i,is) + beta(i+1,is))
            endif

CCCCCCCCCC
C            beta_edge(i,n) = beta_lo
C            beta_edge(i+1,n) = beta_hi
CCCCCCCCCC
            flux_hi(n) = beta_hi*(Y(i+1,n) - Y(i  ,n)) 
            flux_lo(n) = beta_lo*(Y(i  ,n) - Y(i-1,n)) 
            
            visc(i,n) = (flux_hi(n) - flux_lo(n))*dxsqinv
            
            if (LeEQ1 .eq. 0) then
               sum_lo = sum_lo + flux_lo(n)
               sum_hi = sum_hi + flux_hi(n)
               
               sumRhoY_lo = sumRhoY_lo+0.5d0*(scal(i-1,is)+scal(i,is))
               sumRhoY_hi = sumRhoY_hi+0.5d0*(scal(i,is)+scal(i+1,is))
            endif
         enddo

         if (LeEQ1 .eq. 0) then
            do n=1,Nspec
               is = FirstSpec + n - 1
               RhoYe_lo = .5d0*(scal(i-1,is)+scal(i,is))
               RhoYe_hi = .5d0*(scal(i,is)+scal(i+1,is))
               
               flux_lo(n) = flux_lo(n) - sum_lo*RhoYe_lo/sumRhoY_lo
               flux_hi(n) = flux_hi(n) - sum_hi*RhoYe_hi/sumRhoY_hi
               
               visc(i,n) =  (flux_hi(n) - flux_lo(n))*dxsqinv
            end do
         endif
         
      end do

C FIXME
C      do i=0,nx-1
C         write(11,1006) i,
C     &        (beta_edge(i,n)*1.d-1,n=1,Nspec)
C      enddo
C      close(11)
C      write(*,*)'spec diff'
C      stop
        
      end
