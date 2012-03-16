      subroutine get_temp_visc_terms(scal,beta,visc,dx,lo,hi)
      implicit none
      include 'spec.h'
      real*8 scal(-2:nfine+1,nscal)
      real*8 beta(-1:nfine  ,nscal)
      real*8 visc(-1:nfine)
      real*8 dx
      integer lo,hi

c     Compute Div(lambda.Grad(T)) + rho.D.Grad(Hi).Grad(Yi)
      call rhoDgradHgradY(scal,beta,visc,dx,lo,hi)

c     Add Div( lambda Grad(T) )
      call addDivLambdaGradT(scal,beta,visc,dx)

      end

      subroutine addDivLambdaGradT(scal,beta,visc,dx)
      implicit none
      include 'spec.h'
      real*8 scal(-2:nx+1,nscal)
      real*8 beta(-1:nx  ,nscal)
      real*8 visc(-1:nx)
      real*8 dx
      
      integer i
      real*8 beta_lo,beta_hi
      real*8 flux_lo,flux_hi
      real*8 dxsqinv

      dxsqinv = 1.d0/(dx*dx)
      do i = 0,nx-1
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

      subroutine rhoDgradHgradY(scal,beta,visc,dx,lo,hi)
      implicit none
      include 'spec.h'
      real*8 scal(-2:nfine+1,nscal)
      real*8 beta(-1:nfine  ,nscal)
      real*8 visc(-1:nfine)
      real*8 dx
      integer lo,hi
      
      integer i,n,is,IWRK
      real*8 beta_lo,beta_hi
      real*8 rdgydgh_lo,rdgydgh_hi
      real*8 dxsqinv,RWRK,rho,dv
      real*8 hm(Nspec,-1:nfine)
      real*8 Y(Nspec,-1:nfine)

      real*8 spec_flux_lo(0:nfine-1,Nspec)
      real*8 spec_flux_hi(0:nfine-1,Nspec)

      real*8 sum_lo,sum_hi
      real*8 sumRhoY_lo,sumRhoY_hi
      real*8 RhoYe_lo,RhoYe_hi

c     Compute rhoD Grad(Yi).Grad(hi) terms

      dxsqinv = 1.d0/(dx*dx)

c     Get Hi, Yi at cell centers
      do i=lo-1,hi+1
         rho = 0.d0
         do n=1,Nspec
            rho = rho + scal(i,FirstSpec+n-1)
         enddo
         call CKHMS(scal(i,Temp),IWRK,RWRK,hm(1,i))
         do n=1,Nspec
            Y(n,i) = scal(i,FirstSpec+n-1)/rho
         enddo
      enddo

c     Compute differences
      do i=lo,hi
         dv = 0.d0
         sum_lo = 0.d0
         sum_hi = 0.d0
         sumRhoY_lo = 0.d0
         sumRhoY_hi = 0.d0
         do n=1,Nspec
            is = FirstSpec + n - 1
            if (coef_avg_harm.eq.1) then
               beta_lo = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i-1,is))
               beta_hi = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i+1,is))
            else
               beta_lo = 0.5d0*(beta(i,is) + beta(i-1,is))
               beta_hi = 0.5d0*(beta(i,is) + beta(i+1,is))
            endif

            spec_flux_lo(i,n) = beta_lo*(Y(n,i)-Y(n,i-1))
            spec_flux_hi(i,n) = beta_hi*(Y(n,i+1)-Y(n,i))

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

         end do

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

            end do
         end if

         do n=1,Nspec
            rdgydgh_lo = spec_flux_lo(i,n)*(hm(n,i)-hm(n,i-1))
            rdgydgh_hi = spec_flux_hi(i,n)*(hm(n,i+1)-hm(n,i))

            dv = dv + (rdgydgh_hi + rdgydgh_lo)*0.5d0
         enddo
         
         visc(i) = dv*dxsqinv

        end do
      end
