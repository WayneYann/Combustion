!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   advance.f90                          !!!
!!!   August 6, 2015                       !!!
!!!   Framework for MISDC timestep routine !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module advance_module

   use ghost_cells_module
   use div_u_module
   use cell_conversions_module
   use quadrature_module
   use diffusion_correction_module
   
   implicit none
   private
   public :: advance
   
contains
   
   subroutine advance(scal_n_cc, scal_np1_cc, vel, dt, dx)
      implicit none
      include 'spec.h'
      ! parameters
      ! two ghost cells
      double precision, intent(in   ) ::   scal_n_cc(-2:nx+1, nscal)
      double precision, intent(out  ) :: scal_np1_cc(-2:nx+1, nscal)
      ! face-value, no ghost cells
      double precision, intent(out  ) ::         vel( 0:nx)
      double precision, intent(in   ) :: dt, dx
      ! variable declarations

      ! we need dt for each timestep, denoted dtm
      double precision :: dtm(0:nnodes-2)

      ! we need to store our state variables
      ! we store the 'next' and 'previous' MISDC iterations
      ! indicated by k and kp1 (k plus 1)
      ! additionally, we need to store the state variables at 
      ! each temporal node
      ! (two ghost cells)
      double precision ::    scal_k_cc(0:nnodes-1, -2:nx+1, nscal)
      double precision ::  scal_kp1_cc(0:nnodes-1, -2:nx+1, nscal)
      double precision ::   scal_m_avg(-2:nx+1, nscal)

      ! we also need the provisional 'AD' solution
      ! for the next iterate, at the next temporal node
      double precision ::  scal_AD_avg(-2:nx+1, nscal)

      ! for each time substep, we have a delta chi prediction
      ! and correction term
      double precision :: delta_chi_pred(0:nnodes-1, 0:nx-1)
      double precision :: delta_chi_corr(0:nnodes-1, 0:nx-1)

      ! we also need to know the value of the contraint S
      ! (one ghost cell)
      double precision ::  S_cc(0:nx-1)
      double precision :: S_avg(0:nx-1)

      ! we need to store the advection terms for the next 
      ! and previous iterates at all temporal nodes
      double precision ::   advection_k(0:nnodes-1, 0:nx-1, nscal)
      double precision :: advection_kp1(0:nnodes-1, 0:nx-1, nscal)

      ! we need to store the diffusion terms for the next
      ! and previous iterates at all temporal nodes
      double precision ::   diffusion_k(0:nnodes-1, 0:nx-1, nscal)
      double precision :: diffusion_kp1(0:nnodes-1, 0:nx-1, nscal)

      ! we need to store the production rates for the next
      ! and previous iterates at all temporal nodes
      double precision ::   wdot_k(0:nnodes-1, 0:nx-1, Nspec)
      double precision :: wdot_kp1(0:nnodes-1, 0:nx-1, Nspec)
      
      ! we also store the integral approximation from the previous iterate
      double precision :: I_k(0:nnodes-2, 0:nx-1, nscal)
      
      ! beta stores the diffusion coefficients
      ! two ghost cells
      double precision :: beta(-2:nx+1, nscal)
      ! gammas are used to compute the differential diffusion
      double precision :: gamma_lo(0:nx-1, Nspec)
      double precision :: gamma_hi(0:nx-1, Nspec)
      
      ! k is the MISDC iteration
      integer :: k
      ! m is the temporal node
      integer :: m
      ! i is the gridpoint
      ! integer :: i
      
      ! time substeps (determined by Gauss-Lobatto rule)
      dtm(0) = 0.5*dt
      dtm(1) = 0.5*dt
      
      ! initialize delta chi
      delta_chi_pred = 0.d0
      delta_chi_corr = 0.d0
      
      ! initialize the 0th iterate to the initial condition
      scal_k_cc(0,:,:) = scal_n_cc
      
      call fill_scal_cc_ghost_cells(scal_k_cc(0,:,:))
      
      ! compute the advection term
      delta_chi_pred(0,:) = 0.d0
      call increment_delta_chi(delta_chi_pred(0,:), scal_k_cc(0,:,:), dtm(0), dx)
      call compute_diffusion_coefficients(beta, scal_k_cc(0,:,:))
      ! need to compute div u to fourth order
      call compute_div_u(S_cc, scal_k_cc(0,:,:), beta, dx)
      S_cc = S_cc + delta_chi_pred(0,:)
      call extrapolate_cc_to_avg(S_avg, S_cc)
      call compute_velocity(vel, S_avg, dx)
      ! compute the advection term
      ! todo: make sure this is right...
      call compute_advection(advection_k(0,:,:), scal_k_cc(0,:,:), vel, dx)
      ! compute the diffusion term
      call compute_diffusion(diffusion_k(0,:,:), scal_k_cc(0,:,:), beta, gamma_lo, gamma_hi, dx)
      call get_diffdiff_terms(advection_k(0,:,RhoH), scal_k_cc(0,:,:), beta, gamma_lo, gamma_hi, dx)
      ! compute the reaction term
      call compute_production_rate(wdot_k(0,:,:), scal_k_cc(0,:,:))
      
      do m = 1,nnodes-1
         scal_k_cc(m,:,:)   =   scal_k_cc(m-1,:,:)
         advection_k(m,:,:) = advection_k(m-1,:,:)
         diffusion_k(m,:,:) = diffusion_k(m-1,:,:)
         wdot_k(m,:,:)      =      wdot_k(m-1,:,:)
      end do
      
      call compute_integrals(I_k, advection_k, diffusion_k, wdot_k, dt)
      
      scal_kp1_cc(0,:,:) = scal_k_cc(0,:,:)
      advection_kp1(0,:,:) = advection_k(0,:,:)
      diffusion_kp1(0,:,:) = diffusion_k(0,:,:)
      wdot_kp1(0,:,:) = wdot_k(0,:,:)
      
      do k=1,misdc_iterMAX
         do m=0,nnodes-2
            ! convert the cell-centered scalars to cell-average
            call scal_cc_to_avg(scal_m_avg, scal_kp1_cc(m,:,:))
            
            ! compute delta chi prediction
            delta_chi_pred(m,:) = 0.d0
            call increment_delta_chi(delta_chi_pred(m,:), scal_kp1_cc(m,:,:), dtm(m), dx)
            
            ! compute the diffusion coefficients
            call compute_diffusion_coefficients(beta, scal_kp1_cc(m,:,:))
            ! compute the divergence constraint, S
            call compute_div_u(S_cc, scal_kp1_cc(m,:,:), beta, dx)
            
            ! add delta chi to the constraint
            S_cc = S_cc + delta_chi_pred(m,:) + delta_chi_corr(m,:)
            
            ! convert from cell-centered S to cell-average
            call extrapolate_cc_to_avg(S_avg, S_cc)
            ! compute velocity by integrating the constraint
            ! the velocity is computed at faces
            call compute_velocity(vel, S_avg, dx)
            
            ! compute the advection terms A^{m,(k+1)}
            call compute_advection(advection_kp1(m,:,:), scal_kp1_cc(m,:,:), vel, dx)
            call get_diffdiff_terms(advection_kp1(m,:,RhoH), scal_kp1_cc(m,:,:), beta, gamma_lo, gamma_hi, dx)
            
            ! update the density
            call update_density(scal_AD_avg, scal_m_avg, &
                                advection_kp1(m,:,:), advection_k(m,:,:), &
                                I_k(m,:,:), dtm(m))
            
            ! compute the iteratively-lagged diffusion coefficients
            call compute_diffusion_coefficients(beta, scal_k_cc(m+1,:,:))
            
            ! perform the diffusion correction for species
            call species_AD_correction(scal_AD_avg, scal_m_avg, beta, &
                                       advection_kp1(m,:,:), advection_k(m,:,:), &
                                       diffusion_k(m+1,:,:), I_k(m,:,:), dtm(m), dx)
            ! update the provisional enthalpy
            call enthalpy_AD_correction(scal_AD_avg, scal_m_avg, beta, &
                                        advection_kp1(m,:,:), advection_k(m,:,:), &
                                        diffusion_k(m+1,:,:), I_k(m,:,:), dtm(m), dx)
            
            ! convert provisional solution to cell-centered values
            call scal_avg_to_cc(scal_kp1_cc(m+1,:,:), scal_AD_avg(:,:))
            
            ! compute gamma, conservatively corrected diffusion term
            ! call compute_div_gamma
            call compute_diffusion(diffusion_kp1(m+1,:,:), scal_kp1_cc(m+1,:,:), &
                                   beta, gamma_lo, gamma_hi, dx)
            
            ! also need to convert the advection and diffusion to cell-centered 
            ! quantities
            
            !call write_plt(vel, scal_AD_avg, S_cc, dx, 2*m + 1, dt*m/2.0)
            
            ! call the chemistry solver to solve the correction equation
            call reaction_correction(scal_kp1_cc(m+1,:,:),   scal_kp1_cc(m,:,:), &
                                     advection_kp1(m,:,:),   advection_k(m,:,:), &
                                     diffusion_kp1(m+1,:,:), diffusion_k(m+1,:,:), &
                                     wdot_k(m+1,:,:), I_k(m,:,:), dtm(m))
            
            call write_plt(vel, scal_kp1_cc(m+1,:,:), S_cc, dx, 2*k + m, dt*m/2.0)
            
            ! need to also compute the new advection term
            ! new diffusion term is already computed
            ! compute the new production rates
            call compute_production_rate(wdot_kp1(m+1,:,:), scal_kp1_cc(m+1,:,:))
            
            ! increment the delta chi correction
            call increment_delta_chi(delta_chi_pred(m,:), scal_kp1_cc(m+1,:,:), dtm(m), dx)
         end do
         
         m = nnodes-1
         delta_chi_pred(nnodes-1,:) = 0
         call increment_delta_chi(delta_chi_pred(m,:), scal_kp1_cc(m,:,:), dtm(0), dx)
         
         do m=0,nnodes-1
            call compute_diffusion_coefficients(beta, scal_kp1_cc(m,:,:))
            call compute_div_u(S_cc, scal_kp1_cc(m,:,:), beta, dx)
            S_cc = S_cc + delta_chi_pred(m,:) + delta_chi_corr(m,:)
            call extrapolate_cc_to_avg(S_avg, S_cc)
            call compute_velocity(vel, S_avg, dx)
            call compute_diffusion(diffusion_kp1(m,:,:), scal_kp1_cc(m,:,:), &
                                   beta, gamma_lo, gamma_hi, dx)
            call compute_advection(advection_kp1(m,:,:), scal_kp1_cc(m,:,:), vel, dx)
            call get_diffdiff_terms(advection_kp1(m,:,RhoH), scal_kp1_cc(m,:,:), beta, gamma_lo, gamma_hi, dx)
         end do
         
         scal_k_cc = scal_kp1_cc
         advection_k = advection_kp1
         diffusion_k = diffusion_kp1
         wdot_k = wdot_kp1
         ! recompute the integral I_k
         call compute_integrals(I_k, advection_k, diffusion_k, wdot_k, dt)
      end do
      
      stop
      
      ! advance the solution
      scal_np1_cc = scal_kp1_cc(nnodes-1,:,:)
   end subroutine advance
   
end module advance_module
