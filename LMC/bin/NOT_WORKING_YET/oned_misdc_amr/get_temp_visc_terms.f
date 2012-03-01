      subroutine get_temp_visc_terms(scal,beta,visc,dx,time)
      implicit none
      include 'spec.h'
      real*8 scal(-1:nx  ,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(0 :nx-1)
      real*8 dx, time


c     Compute Div(lambda.Grad(T)) + rho.D.Grad(Hi).Grad(Yi)
      call rhoDgradHgradY(scal,beta,visc,dx,time)

c     Add Div( lambda Grad(T) )
      call addDivLambdaGradT(scal,beta,visc,dx,time)

      end

      subroutine addDivLambdaGradT(scal,beta,visc,dx,time)
      implicit none
      include 'spec.h'
      real*8 scal(-1:nx  ,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(0 :nx-1)
      real*8 dx, time
      
      integer i
      real*8 beta_lo,beta_hi
      real*8 flux_lo,flux_hi
      real*8 dxsqinv


c     this function will sometimes be called independently from 
c     get_temp_visc_terms, so need this here
      call set_bc_grow_s(scal,dx,time)

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

 1001 FORMAT(E22.15,1X)      
            

      end

      subroutine rhoDgradHgradY(scal,beta,visc,dx,time)
      implicit none
      include 'spec.h'
      real*8 scal(-1:nx  ,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(0 :nx-1)
      real*8 dx,time
      
      integer i,n,is,IWRK
      real*8 beta_lo,beta_hi
      real*8 rdgydgh_lo,rdgydgh_hi
      real*8 dxsqinv,RWRK,rho,dv
      real*8 hi(maxspec,-1:nx)
      real*8 Y(maxspec,-1:nx)

      real*8 spec_flux_lo(0:nx-1,maxspec)
      real*8 spec_flux_hi(0:nx-1,maxspec)

      real*8 sum_lo,sum_hi
      real*8 sumRhoY_lo,sumRhoY_hi
      real*8 RhoYe_lo,RhoYe_hi

c     Compute rhoD Grad(Yi).Grad(hi) terms

      dxsqinv = 1.d0/(dx*dx)

      call set_bc_grow_s(scal,dx,time)

c     Get Hi, Yi at cell centers
      do i = -1,nx
         rho = 0.d0
         do n=1,Nspec
            rho = rho + scal(i,FirstSpec+n-1)
         enddo
         call CKHMS(scal(i,Temp),IWRK,RWRK,hi(1,i))
         do n=1,Nspec
            Y(n,i) = scal(i,FirstSpec+n-1)/rho
         enddo
      enddo
C      do n = 1,Nspec
C         hi(n,-1) = hi(n,0)
C      enddo

c     Compute differences
      do i = 0,nx-1
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
            rdgydgh_lo = spec_flux_lo(i,n)*(hi(n,i)-hi(n,i-1))
            rdgydgh_hi = spec_flux_hi(i,n)*(hi(n,i+1)-hi(n,i))

            dv = dv + (rdgydgh_hi + rdgydgh_lo)*0.5d0
         enddo
         
         visc(i) = dv*dxsqinv

        end do
      end


      subroutine rhoDgradHgradY_nosetbc(scal,beta,visc,dx,time)
      implicit none
      include 'spec.h'
      real*8 scal(-1:nx  ,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(0 :nx-1)
      real*8 dx,time
      
      integer i,n,is,IWRK
      real*8 beta_lo,beta_hi
      real*8 rdgydgh_lo,rdgydgh_hi
      real*8 dxsqinv,RWRK,rho,dv
      real*8 hi(maxspec,-1:nx)
      real*8 Y(maxspec,-1:nx)

      real*8 spec_flux_lo(0:nx-1,maxspec)
      real*8 spec_flux_hi(0:nx-1,maxspec)

      real*8 sum_lo,sum_hi
      real*8 sumRhoY_lo,sumRhoY_hi
      real*8 RhoYe_lo,RhoYe_hi

c     Compute rhoD Grad(Yi).Grad(hi) terms

      dxsqinv = 1.d0/(dx*dx)

c     Get Hi, Yi at cell centers
      do i = -1,nx
         rho = 0.d0
         do n=1,Nspec
            rho = rho + scal(i,FirstSpec+n-1)
         enddo
         call CKHMS(scal(i,Temp),IWRK,RWRK,hi(1,i))
         do n=1,Nspec
            Y(n,i) = scal(i,FirstSpec+n-1)/rho
         enddo
      enddo
C      do n = 1,Nspec
C         hi(n,-1) = hi(n,0)
C      enddo

c     Compute differences
      do i = 0,nx-1
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
            rdgydgh_lo = spec_flux_lo(i,n)*(hi(n,i)-hi(n,i-1))
            rdgydgh_hi = spec_flux_hi(i,n)*(hi(n,i+1)-hi(n,i))

            dv = dv + (rdgydgh_hi + rdgydgh_lo)*0.5d0
         enddo
         
CCCCCCCCCCCCC
C         write(11,*)
C         write(11,*)
CCCCCCCCCCCCC

         visc(i) = dv*dxsqinv

        end do

CCCCCCCCCCCCCC
C         close(11)
C         write(*,*)'end of step'
C         stop
CCCCCCCCCCCCC      
      end
