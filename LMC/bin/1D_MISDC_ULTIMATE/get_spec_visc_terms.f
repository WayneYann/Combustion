      subroutine get_spec_visc_terms(scal,beta,visc,gamma_lo,
     &                               gamma_hi,dx,lo,hi)

      implicit none
      include 'spec.h'
      real*8         scal(-2:nfine+1,nscal)
      real*8         beta(-1:nfine  ,nscal)
      real*8         visc(-1:nfine  ,Nspec)
      real*8 gamma_lo( 0:nfine-1,Nspec)
      real*8 gamma_hi( 0:nfine-1,Nspec)
      real*8 dx
      integer lo,hi
      
      integer i,n,is
      real*8 beta_lo,beta_hi
      real*8 dxsqinv
      real*8 Y(-1:nfine,Nspec), sum_lo, sum_hi, sumRhoY_lo, sumRhoY_hi
      real*8 RhoYe_lo, RhoYe_hi

      do i=lo-1,hi+1
         do n=1,Nspec
            Y(i,n) = scal(i,FirstSpec+n-1)/scal(i,Density)
         enddo
      enddo

      dxsqinv = 1.d0/(dx*dx)
      do i=lo,hi
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

            gamma_hi(i,n) = beta_hi*(Y(i+1,n) - Y(i  ,n)) 
            gamma_lo(i,n) = beta_lo*(Y(i  ,n) - Y(i-1,n)) 
 
            visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv

            if (LeEQ1 .eq. 0) then

c              need to correct fluxes so they add to zero on each face
c              build up the sum of species fluxes on lo and hi faces
c              this will be "rho * V_c"
               sum_lo = sum_lo + gamma_lo(i,n)
               sum_hi = sum_hi + gamma_hi(i,n)
               
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
               gamma_lo(i,n) = gamma_lo(i,n) 
     $              - sum_lo*RhoYe_lo/sumRhoY_lo
               gamma_hi(i,n) = gamma_hi(i,n) 
     $              - sum_hi*RhoYe_hi/sumRhoY_hi
               
               visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv
            end do
         end if

      end do
        
      end
