      subroutine get_rhoh_visc_terms(scal,beta,visc,dx)

      implicit none
      include 'spec.h'
      real*8 scal(-2:nx+1,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(-1:nx)
      real*8 dx
      
      integer i
      real*8 beta_lo,beta_hi
      real*8 flux_lo,flux_hi
      real*8 dxsqinv, h_lo, h_mid, h_hi

            dxsqinv = 1.d0/(dx*dx)
      do i = 0,nx-1         
         if (coef_avg_harm.eq.1) then
            beta_lo = 2.d0 / (1.d0/beta(i,RhoH)+1.d0/beta(i-1,RhoH))
            beta_hi = 2.d0 / (1.d0/beta(i,RhoH)+1.d0/beta(i+1,RhoH))
         else
            beta_lo = 0.5d0*(beta(i,RhoH) + beta(i-1,RhoH))
            beta_hi = 0.5d0*(beta(i,RhoH) + beta(i+1,RhoH))
         endif

         h_hi  = scal(i+1,RhoH) / scal(i+1,Density)
         h_mid = scal(i  ,RhoH) / scal(i  ,Density)
         h_lo  = scal(i-1,RhoH) / scal(i-1,Density)

         flux_hi = beta_hi*(h_hi - h_mid)
         flux_lo = beta_lo*(h_mid - h_lo)
         visc(i) = (flux_hi - flux_lo) * dxsqinv 
      enddo

      end


      subroutine get_diffdiff_terms(scal_for_coeff,scal_for_grad,
     $                              spec_flux_lo,spec_flux_hi,
     $                              beta,diffdiff,dx)

      implicit none
      include 'spec.h'

      real*8 scal_for_coeff(-2:nx+1,*)
      real*8 scal_for_grad (-2:nx+1,*)
      real*8 spec_flux_lo  ( 0:nx-1,*)
      real*8 spec_flux_hi  ( 0:nx-1,*)
      real*8 beta          (-1:nx  ,*)
      real*8 diffdiff      (-1:nx)
      real*8 dx

      real*8 dxsqinv,RWRK
      integer i,is,n,IWRK
      real*8 hm(Nspec,-1:nx)
      real*8 flux_lo(Nspec),flux_hi(Nspec)
      real*8 Y(Nspec,-1:nx)
      real*8 beta_lo, beta_hi, rho

      dxsqinv = 1.d0/(dx*dx)

      diffdiff = 0.d0

      do i=-1,nx
         rho = 0.d0
c        compute density
         do n=1,Nspec
            rho = rho + scal_for_grad(i,FirstSpec+n-1)
         enddo
c        compute Y = (rho*Y)/rho
         do n=1,Nspec
            Y(n,i) = scal_for_grad(i,FirstSpec+n-1)/rho
         enddo
c        compute cell-centered h_m
         call CKHMS(scal_for_coeff(i,Temp),IWRK,RWRK,hm(1,i))
      end do

      do i=0,nx-1
         do n=1,Nspec
            is = FirstSpec + n - 1

c     compute -lambda/cp on faces
            beta_lo = (-beta(i  ,RhoH)-beta(i-1,RhoH)) /2.d0
            beta_hi = (-beta(i+1,RhoH)-beta(i  ,RhoH)) /2.d0

c     set face fluxes to -lambda/cp * grad Y_m
            flux_lo(n) = beta_lo*(Y(n  ,i) - Y(n,i-1))
            flux_hi(n) = beta_hi*(Y(n,i+1) - Y(n  ,i))

c     set face fluxes to h_m * (rho D_m - lambda/cp) grad Y_m
            flux_lo(n) = (flux_lo(n) + spec_flux_lo(i,n))*
     $           (hm(n,i-1)+hm(n,i))/2.d0
            flux_hi(n) = (flux_hi(n) + spec_flux_hi(i,n))*
     $           (hm(n,i+1)+hm(n,i))/2.d0
 
c     differential diffusion is divergence of face fluxes
            diffdiff(i) = diffdiff(i) + 
     $           (flux_hi(n) - flux_lo(n))*dxsqinv

         end do
      end do

      end

