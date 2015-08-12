module diffusion_correction_module
   
   implicit none
   private
   include 'spec.h'
   public :: species_AD_correction, enthalpy_AD_correction, get_spec_visc_terms
   
contains

   subroutine get_spec_visc_terms(scal,beta,visc,gamma_lo,gamma_hi,dx)
      implicit none
      include 'spec.h'
      double precision, intent(in ) ::     scal(-2:nx+1,nscal)
      double precision, intent(in ) ::     beta(-1:nx  ,nscal)
      double precision, intent(out) ::     visc(-1:nx  ,Nspec)
      double precision, intent(out) :: gamma_lo( 0:nx-1,Nspec)
      double precision, intent(out) :: gamma_hi( 0:nx-1,Nspec)
      double precision, intent(in ) :: dx
      
      integer i,n,is
      double precision :: beta_lo,beta_hi
      double precision :: dxsqinv
      double precision :: Y(-1:nx,Nspec), sum_lo, sum_hi, sumRhoY_lo, sumRhoY_hi
      double precision :: RhoYe_lo, RhoYe_hi
      
      ! this subroutine computes the term div(rho D_m grad Y_m)
      
      do i=-1,nx
         do n=1,Nspec
            Y(i,n) = scal(i,FirstSpec+n-1)/scal(i,Density)
         enddo
      enddo

      dxsqinv = 1.d0/(dx*dx)
      do i=0,nx-1
         sum_lo = 0.d0
         sum_hi = 0.d0
         sumRhoY_lo = 0
         sumRhoY_hi = 0
         do n=1,Nspec
            is = FirstSpec + n - 1
            
            beta_lo = 0.5d0*(beta(i,is) + beta(i-1,is))
            beta_hi = 0.5d0*(beta(i,is) + beta(i+1,is))
            
            gamma_hi(i,n) = beta_hi*(Y(i+1,n) - Y(i  ,n)) 
            gamma_lo(i,n) = beta_lo*(Y(i  ,n) - Y(i-1,n)) 
 
            !visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv
            
            ! need to correct fluxes so they add to zero on each face
            ! build up the sum of species fluxes on lo and hi faces
            ! this will be "rho * V_c"
            sum_lo = sum_lo + gamma_lo(i,n)
            sum_hi = sum_hi + gamma_hi(i,n)
            
            ! build up the sum of rho*Y_m
            ! this will be the density
            sumRhoY_lo = sumRhoY_lo+0.5d0*(scal(i-1,is)+scal(i,is))
            sumRhoY_hi = sumRhoY_hi+0.5d0*(scal(i,is)+scal(i+1,is))

         enddo

            ! correct the fluxes so they add up to zero before computing visc
         do n=1,Nspec
            is = FirstSpec + n - 1

            ! compute rho*Y_m on each face
            RhoYe_lo = .5d0*(scal(i-1,is)+scal(i,is))
            RhoYe_hi = .5d0*(scal(i,is)+scal(i+1,is))

            ! set flux = flux - (rho*V_c)*(rho*Y_m)/rho
            gamma_lo(i,n) = gamma_lo(i,n) - sum_lo*RhoYe_lo/sumRhoY_lo
            gamma_hi(i,n) = gamma_hi(i,n) - sum_hi*RhoYe_hi/sumRhoY_hi
            
            visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv
         end do
      end do
   end subroutine get_spec_visc_terms
   

   subroutine species_AD_correction(scal_AD_avg, scal_m_avg, advection_kp1, &
                                    advection_k, diffusion_k, I_k, dtm)
      implicit none
      double precision, intent(inout) ::  scal_AD_avg(-2:nx+1,nscal)
      double precision, intent(in   ) ::   scal_m_avg(-2:nx+1,nscal)
      double precision, intent(in   ) :: advection_kp1(0 :nx-1,nscal)
      double precision, intent(in   ) ::   advection_k(0 :nx-1,nscal)
      double precision, intent(in   ) ::   diffusion_k(-1:nx, nscal)
      double precision, intent(in   ) ::           I_k(0:nx-1, nscal)
      double precision, intent(in   ) ::           dtm
      
      integer :: n, ispec
      
      do n=1,Nspec
         ispec = FirstSpec+n-1
         
         call implicit_AD_solve(scal_AD_avg(:, ispec), scal_m_avg(:, ispec), &
            advection_kp1(:, ispec), advection_k(:, ispec), &
            diffusion_k(:, ispec), I_k(:, ispec), dtm)
         
      end do
      
   end subroutine species_AD_correction

   subroutine enthalpy_AD_correction(scal_AD_avg, scal_m_avg, advection_kp1, &
                                       advection_k, diffusion_k, I_k, dtm)
      implicit none
      double precision, intent(inout) ::   scal_AD_avg(-2:nx+1,nscal)
      double precision, intent(in   ) ::    scal_m_avg(-2:nx+1,nscal)
      double precision, intent(in   ) :: advection_kp1( 0:nx-1,nscal)
      double precision, intent(in   ) ::   advection_k( 0:nx-1,nscal)
      double precision, intent(in   ) ::   diffusion_k(-1:nx, nscal)
      double precision, intent(in   ) ::           I_k( 0:nx-1, nscal)
      double precision, intent(in   ) ::           dtm
      
      call implicit_AD_solve(scal_AD_avg(:, RhoH), scal_m_avg(:, RhoH), &
         advection_kp1(:, RhoH), advection_k(:, RhoH), &
         diffusion_k(:, RhoH), I_k(:, RhoH), dtm)
      
   end subroutine enthalpy_AD_correction
   
   subroutine implicit_AD_solve(phi_AD_avg, phi_m_avg, advection_kp1, &
                                advection_k, diffusion_k, I_k, dtm)
      implicit none
      double precision, intent(out  ) ::   phi_AD_avg(-2:nx+1)
      double precision, intent(in   ) ::    phi_m_avg(-2:nx+1)
      double precision, intent(in   ) :: advection_kp1(0:nx-1)
      double precision, intent(in   ) ::   advection_k(0:nx-1)
      double precision, intent(in   ) ::  diffusion_k(-1:nx)
      double precision, intent(in   ) ::           I_k(0:nx-1)
      double precision, intent(in   ) ::           dtm
      
      double precision :: rhs(nx)
      !double precision :: phi(nx)
      integer :: i
      
      
      ! construct the right hand side
      do i=0,nx-1
         rhs(i) = phi_m_avg(i) &
            + dtm*(advection_kp1(i) - advection_k(i) - diffusion_k(i)) + I_k(i)
      end do
         
      ! construct the matrix we need to invert
      
      ! perform the sparse linear solve for phi
      
      ! fill in scal_AD_avg with the values from phi
      
      ! fill in the ghost cells
      
   end subroutine implicit_AD_solve
   
end module diffusion_correction_module
