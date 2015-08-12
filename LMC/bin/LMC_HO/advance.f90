      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!   advance.f                            !!!
      !!!   August 6, 2015                       !!!
      !!!   Framework for MISDC timestep routine !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine advance(scal_n_cc, scal_np1_cc, dt, dx, lo, hi)
      implicit none
      ! parameters
      double precision, intent(in   ) ::   scal_n_cc(-2:nx+1, nscal)
      double precision, intent(out  ) :: scal_np1_cc(-2:nx+1, nscal)
      double precision, intent(in   ) :: dt, dx
      integer,          intent(in   ) :: lo, hi
      ! variable declarations

      ! we need dt for each timestep, denoted dtm
      double precision :: dtm(0:nnodes-2)

      ! we need to store our state variables
      ! we store the 'next' and 'previous' MISDC iterations
      ! indicated by k and kp1 (k plus 1)
      ! additionally, we need to store the state variables at 
      ! each temporal node
      ! (two ghost cells)
      double precision ::   scal_k_avg(0:nnodes-1, -2:nx+1, nscal)
      double precision ::    scal_k_cc(0:nnodes-1, -2:nx+1, nscal)
      double precision :: scal_kp1_avg(0:nnodes-1, -2:nx+1, nscal)
      double precision ::  scal_kp1_cc(0:nnodes-1, -2:nx+1, nscal)

      ! we also need the provisional 'AD' solution
      ! for the next iterate, at the next temporal node
      double precision ::  scal_AD_avg(-2:nx+1, nscal)
      double precision ::   scal_AD_cc(-2:nx+1, nscal)

      ! for each time substep, we have a delta chi prediction
      ! and correction term
      double precision :: delta_chi_pred(0:nnodes-2, 0:nx-1)
      double precision :: delta_chi_corr(0:nnodes-2, 0:nx-1)
      
      ! we need the velocity of the next and previous iterates 
      ! (computed by integrating the constraint div U = S)
      ! (no ghost cells)
      double precision ::   vel_k(0:nnodes-1, 0:nx-1)
      double precision :: vel_kp1(0:nnodes-1, 0:nx-1)

      ! we also need to know the value of the contraint S
      ! (one ghost cell)
      double precision ::  S_cc(-1:nx)
      double precision :: S_avg(-1:nx)

      ! we need to store the advection terms for the next 
      ! and previous iterates at all temporal nodes
      double precision ::   advection_k(0:nnodes-1, 0:nx-1, nscal)
      double precision :: advection_kp1(0:nnodes-1, 0:nx-1, nscal)

      ! we need to store the diffusion terms for the next
      ! and previous iterates at all temporal nodes
      double precision ::   diffusion_k(0:nnodes-1, -1:nx, nscal)
      double precision :: diffusion_kp1(0:nnodes-1, -1:nx, nscal)
      double precision ::  diffusion_AD(-1:nx, nscal)

      ! we need to store the production rates for the next
      ! and previous iterates at all temporal nodes
      double precision ::   wdot_k(0:nnodes-1, 0:nx-1, Nspec)
      double precision :: wdot_kp1(0:nnodes-1, 0:nx-1, Nspec)
      
      ! we also store the integral approximation from the previous iterate
      double precision :: I_k(0:nnodes-2, 0:nx-1, nscal)
      
      ! beta is the diffusion coefficients
      double precision :: beta(-1:nx, nscal)
      double precision ::   mu(-1:nx)
      
      ! k is the MISDC iteration
      integer :: k
      ! m is the temporal node
      integer :: m
      ! i is the gridpoint
      integer :: i
      
      ! time substeps (determined by Gauss-Lobatto rule)
      dtm(0) = 0.5*dt
      dtm(1) = 0.5*dt
      
      ! initialize delta chi
      delta_chi_pred = 0.d0
      delta_chi_corr = 0.d0
      
      ! initialize the 0th iterate to the initial condition
      do m=0,nnodes-1
         scal_k_cc(m,:,:) = scal_n_cc
         
         ! compute the advection term
         delta_chi_pred(m,:) = 0.d0
         call increment_delta_chi(delta_chi_pred(m,:), scal_k_cc(m,:,:), dtm(m), dx, lo, hi)
         call compute_diffusion_coefficients(beta, mu, scal_k_cc(m,:,:), lo, hi)
         call compute_div_u(S_cc, scal_k_cc(m,:,:), beta, dx, lo, hi)
         S_cc = S_cc + delta_chi_pred(m,:)
         call cell_center_to_avg(S_avg, S_cc)
         call compute_velocity(vel_k(m,:), S_avg, dx, lo, hi)
         call compute_advection(advection_k(m,:,:), scal_k_cc(m,:,:), vel_k(m,:), dx, lo, hi)
         
         ! compute the diffusion term
         
         ! compute the reaction term
         call compute_production_rate(wdot_k(m,:,:), scal_k_cc(m,:,:), lo, hi)
      end do
      call compute_integrals(I_k, advection_k, diffusion_k, wdot_k, dt)
      
      scal_kp1_cc(0,:,:) = scal_k_cc(0,:,:)
      advection_kp1(0,:,:) = advection_k(0,:,:)
      diffusion_kp1(0,:,:) = diffusion_k(0,:,:)
      wdot_kp1(0,:,:) = wdot_k(0,:,:)
      
      do k=1,MISDC_iters
         do m=0,nnodes-2
            ! compute delta chi prediction
            delta_chi_pred(m,:) = 0.d0
            call increment_delta_chi(delta_chi_pred(m,:), scal_kp1_cc(m,:,:), dtm(m), dx, lo, hi)
            
            ! compute the diffusion coefficients
            call compute_diffusion_coefficients(beta, mu, scal_kp1_cc(m,:,:), lo, hi)
            ! compute the divergence constraint, S
            call compute_div_u(S_cc, scal_kp1_cc(m,:,:), beta, dx, lo, hi)
            
            ! add delta chi to the constraint
            S_cc = S_cc + delta_chi_pred(m,:) + delta_chi_corr(m,:)
            
            ! convert from cell-centered S to cell-average
            call cell_center_to_avg(S_avg, S_cc)
            
            ! compute velocity by integrating the constraint
            ! the velocity is computed at faces
            call compute_velocity(vel_kp1(m,:), S_avg, dx, lo, hi)
            
            ! compute the advection terms A^{m,(k+1)}
            call compute_advection(advection_kp1(m,:,:), scal_kp1_cc(m,:,:), vel_kp1(m,:), dx, lo, hi)
            
            ! update the density
            call update_density(scal_kp1_avg(m+1,:,:), scal_kp1_avg(m,:,:), &
                                advection_kp1(m,:,:), advection_k(m,:,:), &
                                I_k(m,:,:), dtm(m), lo, hi)
            
            ! perform the diffusion correction for species
            scal_AD_avg = scal_kp1_avg(m+1,:,:)
            call species_AD_correction(scal_AD_avg, scal_kp1_avg(m,:,:), &
                                       advection_kp1(m,:,:), advection_k(m,:,:), &
                                       diffusion_k(m+1,:,:), I_k(m,:,:), dtm(m), lo, hi)
            
            ! compute gamma, conservatively corrected diffusion term
            call compute_div_gamma
            
            ! update the provisional enthalpy
            call enthalpy_AD_correction(scal_AD_avg, scal_kp1_avg(m,:,:), &
                                        advection_kp1(m,:,:), advection_k(m,:,:), &
                                        diffusion_k(m+1,:,:), I_k(m,:,:), dtm(m), lo, hi)
            ! compute the AD diffusion term
            
            ! convert density to cell-centered value
            
            ! call the chemistry solver to solve the correction equation
            call reaction_correction(scal_kp1_cc(m+1,:,:), scal_kp1_cc(m,:,:), &
                                     advection_kp1(m,:,:), advection_k(m,:,:), &
                                     diffusion_kp1(m+1,:,:), diffusion_k(m+1,:,:), &
                                     wdot_k(m+1,:,:), I_k, dtm(m))
            
            ! need to also compute the new advection term
            ! new diffusion term is already computed
            ! compute the new production rates
            call compute_production_rate(wdot_kp1(m+1,:,:), scal_kp1_cc(m+1,:,:), lo, hi)
            
            ! increment the delta chi correction
            call increment_delta_chi(delta_chi_pred(m,:), scal_kp1_cc(m+1,:,:), dtm(m), dx, lo, hi)
         end do
         
         scal_k_avg = scal_kp1_avg
         scal_k_cc = scal_k_cc
         advection_k = advection_kp1
         diffusion_k = diffusion_kp1
         wdot_k = wdot_kp1
         ! recompute the integral I_k
         call compute_integral(I_k, advection_k, diffusion_k, wdot_k, dt)
      end do
      
      ! advance the solution
      scal_np1_cc = scal_kp1_cc(nnodes-1,:,:)

      end subroutine advance
