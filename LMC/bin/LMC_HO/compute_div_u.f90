module div_u_module

   use diffusion_correction_module
   
   implicit none
   
   private
   public :: compute_div_u
contains

   subroutine get_temp_visc_terms(scal,beta,visc,dx)
      implicit none
      include 'spec.h'
      double precision, intent(in ) :: scal(-2:nx+1,nscal)
      double precision, intent(in ) :: beta(-1:nx  ,nscal)
      double precision, intent(out) :: visc(-1:nx)
      double precision, intent(in ) :: dx

      ! Compute Div(lambda.Grad(T)) + rho.D.Grad(Hi).Grad(Yi)
      call gamma_dot_gradh(scal,beta,visc,dx)
      ! Add Div( lambda Grad(T) )
      call addDivLambdaGradT(scal,beta,visc,dx)
      
   end subroutine get_temp_visc_terms

   subroutine addDivLambdaGradT(scal,beta,visc,dx)
      implicit none
      include 'spec.h'
      double precision, intent(in   ) :: scal(-2:nx+1,nscal)
      double precision, intent(in   ) :: beta(-1:nx  ,nscal)
      double precision, intent(inout) :: visc(-1:nx)
      double precision, intent(in   ) ::  dx
   
      integer i
      double precision :: beta_lo,beta_hi
      double precision :: flux_lo,flux_hi
      double precision :: dxsqinv

      dxsqinv = 1.d0/(dx*dx)
      do i=0,nx-1
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
      end do
   end subroutine addDivLambdaGradT

   subroutine gamma_dot_gradh(scal,beta,visc,dx)
      implicit none
      include 'spec.h'
      double precision, intent(in ) :: scal(-2:nx+1,nscal)
      double precision, intent(in ) :: beta(-1:nx  ,nscal)
      double precision, intent(out) :: visc(-1:nx)
      double precision, intent(in ) :: dx
   
      integer i,n,is,IWRK
      double precision :: beta_lo,beta_hi
      double precision :: gamma_dot_gradh_lo,gamma_dot_gradh_hi
      double precision :: dxsqinv,RWRK,rho,dv
      double precision :: hm(Nspec,-1:nx)
      double precision :: Y(Nspec,-1:nx)

      double precision :: gamma_lo(0:nx-1,Nspec)
      double precision :: gamma_hi(0:nx-1,Nspec)

      double precision :: sum_lo,sum_hi
      double precision :: sumRhoY_lo,sumRhoY_hi
      double precision :: RhoYe_lo,RhoYe_hi

      !    Compute rhoD Grad(Yi).Grad(hi) terms
      dxsqinv = 1.d0/(dx*dx)

      !  Get Hi, Yi at cell centers
      do i=-1,nx
         rho = 0.d0
         do n=1,Nspec
            rho = rho + scal(i,FirstSpec+n-1)
         enddo
         call CKHMS(scal(i,Temp),IWRK,RWRK,hm(1,i))
         do n=1,Nspec
            Y(n,i) = scal(i,FirstSpec+n-1)/rho
         enddo
      enddo

      !  Compute differences
      do i=0,nx-1
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

            gamma_lo(i,n) = beta_lo*(Y(n,i)-Y(n,i-1))
            gamma_hi(i,n) = beta_hi*(Y(n,i+1)-Y(n,i))

            ! need to correct fluxes so they add to zero on each face
            ! build up the sum of species fluxes on lo and hi faces
            ! this will be "rho * V_c"
            sum_lo = sum_lo + gamma_lo(i,n)
            sum_hi = sum_hi + gamma_hi(i,n)
         
            ! build up the sum of rho*Y_m
            ! this will be the density
            sumRhoY_lo = sumRhoY_lo+0.5d0*(scal(i-1,is)+scal(i,is))
            sumRhoY_hi = sumRhoY_hi+0.5d0*(scal(i,is)+scal(i+1,is))
         end do

         ! correct the fluxes so they add up to zero before computing visc
         do n=1,Nspec
            is = FirstSpec + n - 1

            ! compute rho*Y_m on each face
            RhoYe_lo = .5d0*(scal(i-1,is)+scal(i,is))
            RhoYe_hi = .5d0*(scal(i,is)+scal(i+1,is))

            ! set flux = flux - (rho*V_c)*(rho*Y_m)/rho
            gamma_lo(i,n) = gamma_lo(i,n) - sum_lo*RhoYe_lo/sumRhoY_lo
            gamma_hi(i,n) = gamma_hi(i,n) - sum_hi*RhoYe_hi/sumRhoY_hi
            
            gamma_dot_gradh_lo = gamma_lo(i,n)*(hm(n,i)-hm(n,i-1))
            gamma_dot_gradh_hi = gamma_hi(i,n)*(hm(n,i+1)-hm(n,i))
            dv = dv + (gamma_dot_gradh_hi + gamma_dot_gradh_lo)*0.5d0
          enddo
     
         visc(i) = dv*dxsqinv
      end do
   end

   subroutine compute_div_u(S, scal, beta, dx)
      implicit none
      include 'spec.h'
      ! parameters
      double precision, intent (out  ) :: S(0:nx-1)
      double precision, intent (in   ) :: scal(-2:nx+1,nscal)
      double precision, intent (in   ) :: beta(-1:nx,  nscal)
      double precision, intent (in   ) :: dx
   
      ! local variables
      double precision :: diff(-1:nx,nscal)
      double precision :: gamma_lo(0:nx-1,Nspec)
      double precision :: gamma_hi(0:nx-1,Nspec)
   
      double precision ::    Y(Nspec)
      double precision ::   HK(Nspec)
      double precision :: wdot(Nspec)
   
      double precision :: rho
      double precision :: T
      double precision :: mwmix
      double precision :: cpmix
   
      double precision :: rwrk
      integer :: iwrk
      integer :: i, n
   
      call get_temp_visc_terms(scal, beta, diff(:,Temp), dx)
      call get_spec_visc_terms(scal, beta, diff(:,FirstSpec:), &
                               gamma_lo,gamma_hi,dx)

      do i=0,nx-1
         rho = scal(i,Density)
         do n = 1,Nspec
            Y(n) = scal(i,FirstSpec + n - 1) / rho
         enddo
         T = scal(i,Temp)
      
         ! compute the production rates, wdot
         call CKWYR(rho, T, Y, IWRK, RWRK, wdot)
      
         call CKMMWY(Y,IWRK,RWRK,mwmix)
         call CKCPBS(T,Y,IWRK,RWRK,cpmix)
         call CKHMS(T,IWRK,RWRK,HK)

         ! construct S (according to the formula from DB99 paper)
         ! add the temperature contribution
         S(i) = diff(i,Temp)/(rho*cpmix*T)
         
         ! remove this!!!
         ! wdot = 0.d0
         
         ! add each of the species contributions
         do n=1,Nspec
            S(i) = S(i) + (diff(i,FirstSpec+n-1) + wdot(n)*mwt(n))*invmwt(n)*mwmix/rho &
                        - HK(n)*wdot(n)*mwt(n)/(rho*cpmix*T)
         enddo
       enddo
    end subroutine compute_div_u
    
end module div_u_module
