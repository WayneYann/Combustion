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
   subroutine get_temp_visc_terms(visc_avg,scal_cc,beta_cc,gamma_face,dx)
      implicit none
      include 'spec.h'
      double precision, intent(out) ::       visc_avg( 0:nx-1)
      double precision, intent(in ) ::       scal_cc(-2:nx+1,nscal)
      double precision, intent(in ) ::       beta_cc(-2:nx+1,nscal)
      double precision, intent(in ) :: gamma_face( 0:nx,  Nspec)
      double precision, intent(in ) ::         dx

      ! Compute Div(lambda.Grad(T)) + rho.D.Grad(Hi).Grad(Yi)
      call gamma_dot_gradh(visc_avg,scal_cc,beta_cc,gamma_face,dx)

      ! Add Div( lambda Grad(T) )
      call addDivLambdaGradT(visc_avg,scal_cc,beta_cc,dx)

   end subroutine get_temp_visc_terms

   subroutine addDivLambdaGradT(visc_avg,scal_cc,beta_cc,dx)
      implicit none
      include 'spec.h'
      double precision, intent(inout) :: visc_avg( 0:nx-1)
      double precision, intent(in   ) :: scal_cc(-2:nx+1,nscal)
      double precision, intent(in   ) :: beta_cc(-2:nx+1,nscal)
      double precision, intent(in   ) ::  dx
   
      integer i
      double precision :: beta_face(0:nx)
      double precision :: grad_T(0:nx)
      
      call cc_to_face(beta_face, beta_cc(:,Temp))
      call cc_to_grad(grad_T, scal_cc(:,Temp), dx)
      
      do i=0,nx-1
         visc_avg(i) = visc_avg(i) + (beta_face(i+1)*grad_T(i+1) - beta_face(i)*grad_T(i))/dx
      end do
   end subroutine addDivLambdaGradT
   
   subroutine gamma_dot_gradh(visc_avg,scal_cc,beta_cc,gamma_face,dx)
      implicit none
      include 'spec.h'
      double precision, intent(out) ::       visc_avg( 0:nx-1)
      double precision, intent(in ) ::       scal_cc(-2:nx+1,nscal)
      double precision, intent(in ) ::       beta_cc(-2:nx+1,nscal)
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
         call ckhms(scal_cc(i,Temp),iwrk,rwrk,h_cc(i,:))
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
      
      do i=0,nx-1
         do n=1,Nspec
            visc_avg(i) = &
                 (h_face(i+1,n)*gamma_face(i+1,n) - h_face(i,n)*gamma_face(i,n))/dx & 
               - h_div_gamma(i,n)
         end do
      end do
   end subroutine gamma_dot_gradh

   subroutine compute_div_u(S_cc, scal_cc, beta_cc, dx)
      implicit none
      include 'spec.h'
      ! parameters
      double precision, intent (out  ) :: S_cc(0:nx-1)
      double precision, intent (in   ) :: scal_cc(-2:nx+1,nscal)
      double precision, intent (in   ) :: beta_cc(-2:nx+1,nscal)
      double precision, intent (in   ) :: dx
   
      ! local variables
      double precision :: diff_avg(0:nx-1,nscal)
      double precision :: diff_cc(0:nx-1,nscal)
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

      diff_avg = 0.d0
      
      ! compute the viscous terms
      call get_spec_visc_terms(diff_avg(:,FirstSpec:), scal_cc, beta_cc, &
                               gamma_face,dx)
      call get_temp_visc_terms(diff_avg(:,Temp), scal_cc, beta_cc, gamma_face, dx)
      
      do n=1,nscal
         call extrapolate_avg_to_cc(diff_cc(:,n),diff_avg(:,n))
      end do

      do i=0,nx-1
         rho = scal_cc(i,Density)
         do n = 1,Nspec
            Y(n) = scal_cc(i,FirstSpec + n - 1) / rho
         enddo
         T = scal_cc(i,Temp)
      
         ! compute the production rates, wdot
         call CKWYR(rho, T, Y, IWRK, RWRK, wdot)
      
         call CKMMWY(Y,IWRK,RWRK,mwmix)
         call CKCPBS(T,Y,IWRK,RWRK,cpmix)
         call CKHMS(T,IWRK,RWRK,HK)

         ! construct S (according to the formula from DB99 paper)
         ! add the temperature contribution
         S_cc(i) = diff_cc(i,Temp)/(rho*cpmix*T)
         
         ! add each of the species contributions
         do n=1,Nspec
            S_cc(i) = S_cc(i) + (diff_cc(i,FirstSpec+n-1) + wdot(n)*mwt(n))*invmwt(n)*mwmix/rho &
                 - HK(n)*wdot(n)*mwt(n)/(rho*cpmix*T)
         end do
       end do
       
    end subroutine compute_div_u
    
end module div_u_module
