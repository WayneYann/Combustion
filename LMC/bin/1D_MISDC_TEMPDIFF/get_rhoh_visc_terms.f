      subroutine get_diffdiff_terms(scal_for_coeff,scal_for_grad,
     $                              gamma_lo,gamma_hi,
     $                              beta,diffdiff,dx,lo,hi)

      implicit none
      include 'spec.h'

      real*8 scal_for_coeff(-2:nfine+1,nscal)
      real*8 scal_for_grad (-2:nfine+1,nscal)
      real*8 gamma_lo  ( 0:nfine-1,Nspec)
      real*8 gamma_hi  ( 0:nfine-1,Nspec)
      real*8 beta          (-1:nfine  ,nscal)
      real*8 diffdiff      (-1:nfine)
      real*8 dx
      integer lo,hi

      real*8 dxsqinv,RWRK
      integer i,is,n,IWRK
      real*8 hm(Nspec,-1:nfine)
      real*8 flux_lo(Nspec),flux_hi(Nspec)
      real*8 Y(Nspec,-1:nfine)
      real*8 rho

      dxsqinv = 1.d0/(dx*dx)

      diffdiff = 0.d0

      do i=lo-1,hi+1
c        compute cell-centered h_m
         call CKHMS(scal_for_coeff(i,Temp),IWRK,RWRK,hm(1,i))
      end do

      do i=lo,hi
         do n=1,Nspec
            is = FirstSpec + n - 1

c     set face fluxes to h_m * gamma_m
            flux_lo(n) = gamma_lo(i,n)*(hm(n,i-1)+hm(n,i))/2.d0
            flux_hi(n) = gamma_hi(i,n)*(hm(n,i+1)+hm(n,i))/2.d0
 
c     differential diffusion is divergence of face fluxes
            diffdiff(i) = diffdiff(i) + 
     $           (flux_hi(n) - flux_lo(n))*dxsqinv

         end do
      end do

      end

