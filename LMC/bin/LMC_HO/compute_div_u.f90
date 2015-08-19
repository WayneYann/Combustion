module div_u_module

   use cell_conversions_module
   use ghost_cells_module
   use diffusion_correction_module
   
   implicit none
   
   private
   public :: compute_div_u
contains
   
   ! compute the temperature viscous terms
   ! div(lambda grad T) + sum_m (rho D_m grad Y_m . grad h_m)
   ! from equation (7), DB99
   subroutine get_temp_visc_terms(visc,scal,beta,gamma_face,dx)
      implicit none
      include 'spec.h'
      double precision, intent(out) ::       visc( 0:nx-1)
      double precision, intent(in ) ::       scal(-2:nx+1,nscal)
      double precision, intent(in ) ::       beta(-2:nx+1,nscal)
      double precision, intent(in ) :: gamma_face( 0:nx,  Nspec)
      double precision, intent(in ) ::         dx

      ! Compute Div(lambda.Grad(T)) + rho.D.Grad(Hi).Grad(Yi)
      call gamma_dot_gradh(visc,scal,beta,gamma_face,dx)
      !call gamma_dot_gradh_old(visc,scal,beta,dx)
      ! Add Div( lambda Grad(T) )
      call addDivLambdaGradT(visc,scal,beta,dx)
      !call addDivLambdaGradT_old(visc,scal,beta,dx)
   end subroutine get_temp_visc_terms

   subroutine addDivLambdaGradT(visc,scal,beta,dx)
      implicit none
      include 'spec.h'
      double precision, intent(inout) :: visc( 0:nx-1)
      double precision, intent(in   ) :: scal(-2:nx+1,nscal)
      double precision, intent(in   ) :: beta(-2:nx+1,nscal)
      double precision, intent(in   ) ::  dx
   
      integer i
      double precision :: beta_face(0:nx)
      double precision :: grad_T(0:nx)
      
      call cc_to_face(beta_face, beta(:,Temp))
      call cc_to_grad(grad_T, scal(:,Temp), dx)
      
      do i=0,nx-1
         visc(i) = (beta_face(i+1)*grad_T(i+1) - beta_face(i)*grad_T(i))/dx
      end do
   end subroutine addDivLambdaGradT
   
   subroutine addDivLambdaGradT_old(scal,beta,visc,dx)
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
   end
   
   subroutine gamma_dot_gradh(visc,scal,beta,gamma_face,dx)
      implicit none
      include 'spec.h'
      double precision, intent(out) ::       visc( 0:nx-1)
      double precision, intent(in ) ::       scal(-2:nx+1,nscal)
      double precision, intent(in ) ::       beta(-2:nx+1,nscal)
      double precision, intent(in ) :: gamma_face( 0:nx,  Nspec)
      double precision, intent(in ) ::         dx
      
      double precision ::        h_cc(-2:nx+1,Nspec)
      double precision ::       h_avg(-2:nx+1,Nspec)
      double precision ::      h_face( 0:nx,  Nspec)
      double precision ::   div_gamma(-2:nx+1,Nspec)
      double precision :: h_div_gamma( 0:nx-1,Nspec)
      
      double precision :: h_bdry(Nspec)
      
      double precision :: rwrk
      integer :: i, n
      integer :: iwrk
      
      ! get the boundary condition for the enthalpies
      ! using the Dirichlet condition for temperature
      call ckhms(T_bc(on_lo), iwrk, rwrk, h_bdry)
      
      ! compute the enthalpies at every cell center
      do i=-2,nx+1
         call ckhms(scal(i,Temp),iwrk,rwrk,h_cc(i,:))
      end do
      
      do n=1,Nspec
         ! convert the enthalpies to face values
         call cc_to_face(h_face(:,n), h_cc(:,n))
         ! convert the enthalpies to cell-averaged values
         call cc_to_avg(h_avg(:,n), h_cc(:,n), h_bdry(n))
         
         ! compute div(gamma) as a cell average
         do i=0,nx-1
            div_gamma(i,n) = (gamma_face(i+1,n) - gamma_face(i,n))/dx
         end do
         ! fill in the ghost cells for div(gamma) using extrapolation
         call extrapolate_avg_ghost_cells(div_gamma(:,n))
         call mult_avgs(h_div_gamma(:,n), h_avg(:,n), div_gamma(:,n))
      end do
      
      visc = 0
      do i=0,nx-1
         do n=1,Nspec
            visc(i) = visc(i) &
               + (h_face(i+1,n)*gamma_face(i+1,n) - h_face(i,n)*gamma_face(i,n))/dx & 
               - h_div_gamma(i,n)
         end do
      end do
   end subroutine gamma_dot_gradh

   subroutine compute_div_u(S, scal, beta, dx)
      implicit none
      include 'spec.h'
      ! parameters
      double precision, intent (out  ) :: S(0:nx-1)
      double precision, intent (in   ) :: scal(-2:nx+1,nscal)
      double precision, intent (in   ) :: beta(-2:nx+1,nscal)
      double precision, intent (in   ) :: dx
   
      ! local variables
      double precision :: diff(0:nx-1,nscal)
      double precision :: gamma_face(0:nx,Nspec)
   
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
      
      ! compute the viscous terms
      call get_spec_visc_terms(diff(:,FirstSpec:), scal, beta, &
                               gamma_face,dx)
      call get_temp_visc_terms(diff(:,Temp), scal, beta, gamma_face, dx)
      
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
         
         ! add each of the species contributions
         do n=1,Nspec
            S(i) = S(i) + (diff(i,FirstSpec+n-1) + wdot(n)*mwt(n))*invmwt(n)*mwmix/rho &
                        - HK(n)*wdot(n)*mwt(n)/(rho*cpmix*T)
         end do
       end do
    end subroutine compute_div_u
    
end module div_u_module
