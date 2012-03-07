      subroutine get_spec_visc_terms(scal,beta,visc,spec_flux_lo,
     &                               spec_flux_hi,dx)

      implicit none
      include 'spec.h'
      real*8 scal(-2:nx+1,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(0 :nx-1,*)
      real*8 spec_flux_lo(0:nx-1,*)
      real*8 spec_flux_hi(0:nx-1,*)
      real*8 dx
      
      integer i,n,is
      real*8 beta_lo,beta_hi
      real*8 dxsqinv
      real*8 Y(-1:nx,maxspec), sum_lo, sum_hi, sumRhoY_lo, sumRhoY_hi
      real*8 RhoYe_lo, RhoYe_hi

      do i = -1,nx
         do n=1,Nspec
            Y(i,n) = scal(i,FirstSpec+n-1)/scal(i,Density)
         enddo
      enddo

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

            spec_flux_hi(i,n) = beta_hi*(Y(i+1,n) - Y(i  ,n)) 
            spec_flux_lo(i,n) = beta_lo*(Y(i  ,n) - Y(i-1,n)) 
 
            visc(i,n) = (spec_flux_hi(i,n)-spec_flux_lo(i,n))*dxsqinv

            if (LeEQ1 .eq. 0) then

c              need to correct fluxes so they add to zero on each face
c              build up the sum of species fluxes on lo and hi faces
c              this will be "rho * V_c"
               sum_lo = sum_lo + spec_flux_lo(i,n)
               sum_hi = sum_hi + spec_flux_hi(i,n)
               
c              build up the sum of rho*Y_m
c              this will be the density
               sumRhoY_lo = sumRhoY_lo+0.5d0*(scal(i-1,is)+scal(i,is))
               sumRhoY_hi = sumRhoY_hi+0.5d0*(scal(i,is)+scal(i+1,is))
               
            end if

         enddo

         if (LeEQ1 .eq. 0) then
c           correct the fluxes so they add up to zero before computing visc
            do n=1,Nspec
               is = FirstSpec + n - 1

c              compute rho*Y_m on each face
               RhoYe_lo = .5d0*(scal(i-1,is)+scal(i,is))
               RhoYe_hi = .5d0*(scal(i,is)+scal(i+1,is))

c              set flux = flux - (rho*V_c)*(rho*Y_m)/rho
               spec_flux_lo(i,n) = spec_flux_lo(i,n) 
     $              - sum_lo*RhoYe_lo/sumRhoY_lo
               spec_flux_hi(i,n) = spec_flux_hi(i,n) 
     $              - sum_hi*RhoYe_hi/sumRhoY_hi
               
               visc(i,n) = (spec_flux_hi(i,n)-spec_flux_lo(i,n))*dxsqinv
            end do
         end if

      end do
        
      end
