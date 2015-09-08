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
   use compute_advection_module
   
   implicit none
   private
   include 'spec.h'
   public :: advance
   
contains
   
   subroutine get_temp(scal_avg)
      double precision, intent(inout) :: scal_avg(-2:nx+1, nscal)
      
      double precision :: scal_cc(-2:nx+1, nscal)
      double precision :: Y(NSpec)
      double precision :: hmix
      integer :: n, i
      
      ! T from H parameters
      integer NiterMax, Niter
      parameter (NiterMax = 30)
      double precision res(NiterMax), errMax
      
      errMax = hmix_Typ*1.e-20
      
      call scal_avg_to_cc(scal_cc, scal_avg)
      
      do i=0,nx-1
         ! compute the temperature from h and Y
         do n = 1,Nspec
            Y(n) = scal_cc(i,FirstSpec+n-1)/scal_cc(i,Density)
         enddo
         hmix = scal_cc(i,RhoH) / scal_cc(i,Density)
         
         ! get the new value for the temperature
         call FORT_TfromHYpt(scal_cc(i,Temp), hmix, Y, &
                             Nspec, errMax, NiterMax, res, Niter)
         if (Niter.lt.0) then
            print *,'vodeF_T_RhoY: H to T solve failed'
            print *,'Niter=',Niter
            stop
         endif
      end do
      
      call fill_cc_ghost_cells(scal_cc(:,Temp), T_bc(on_lo))
      call cc_to_avg(scal_avg(:,Temp), scal_cc(:,Temp), T_bc(on_lo))
   end subroutine get_temp
   
   subroutine advance(scal_n_avg, scal_np1_avg, vel, S_avg, wdot, dt, dx)
      implicit none
      double precision, intent(in   ) ::   scal_n_avg(-2:nx+1, nscal)
      double precision, intent(out  ) :: scal_np1_avg(-2:nx+1, nscal)
      double precision, intent(out  ) ::          vel( 0:nx)
      double precision, intent(out  ) ::        S_avg( 0:nx-1)
      double precision, intent(inout) ::         wdot( 0:nx-1, Nspec)
      double precision, intent(in   ) :: dt, dx

      ! we need dt for each timestep, denoted dtm
      double precision :: dtm(0:nnodes-2)

      ! we need to store our state variables
      ! we store the 'next' and 'previous' MISDC iterations
      ! indicated by k and kp1 (k plus 1)
      ! additionally, we need to store the state variables at 
      ! each temporal node
      ! (two ghost cells)
      double precision ::    scal_k_avg(0:nnodes-1, -2:nx+1, nscal)
      double precision ::  scal_kp1_avg(0:nnodes-1, -2:nx+1, nscal)
      double precision ::   scal_kp1_cc(0:nnodes-1, -2:nx+1, nscal)
      double precision ::     scal_k_cc(0:nnodes-1, -2:nx+1, nscal)

      ! we also need the provisional 'AD' solution
      ! for the next iterate, at the next temporal node
      double precision ::  scal_AD_avg(-2:nx+1, nscal)

      ! for each time substep, we have a delta chi prediction
      ! and correction term
      double precision :: delta_chi(0:nnodes-1, 0:nx-1)

      ! we also need to know the value of the contraint S
      ! (no ghost cells)
      double precision ::  S_cc(0:nx-1)

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
      double precision :: I_k_avg(0:nnodes-2, 0:nx-1, nscal)
      double precision ::  I_k_cc(0:nnodes-2, 0:nx-1, nscal)
      
      ! beta stores the diffusion coefficients
      ! two ghost cells
      double precision :: beta(-2:nx+1, nscal)
      ! gamma stores rho D_j grad Y_j
      double precision :: gamma_face(0:nx, Nspec)
      
      ! k is the MISDC iteration
      integer :: k
      ! m is the temporal node
      integer :: m
      ! i is the gridpoint
      integer :: i

      ! time the routine
      double precision :: wallclock_time1, wallclock_time2
      
      call cpu_time(wallclock_time1)
      
      ! time substeps (determined by Gauss-Lobatto rule)
      if (nnodes .eq. 2) then
         dtm(0) = dt
      else if (nnodes .eq. 3) then
         dtm(0) = 0.5d0*dt
         dtm(1) = 0.5d0*dt
      end if
      
      ! initialize everything to zero
      delta_chi = 0
      advection_k = 0
      diffusion_kp1 = 0
      diffusion_k = 0
      diffusion_kp1 = 0
      wdot_k = 0
      wdot_kp1 = 0
      
      ! initialize the 0th iterate to the initial condition
      scal_k_avg(0,:,:) = scal_n_avg
      call fill_scal_avg_ghost_cells(scal_k_avg(0,:,:))
      call scal_avg_to_cc(scal_k_cc(0,:,:), scal_k_avg(0,:,:))
      
      ! compute the advection term
      call increment_delta_chi(delta_chi(0,:), scal_k_cc(0,:,:), dtm(0), dx)
      call compute_diffusion_coefficients(beta, scal_k_cc(0,:,:))
      
      ! compute the reaction term
      wdot_k(0,:,:) = wdot
      
      ! need to compute div u to fourth order
      call compute_div_u(S_cc, scal_k_avg(0,:,:), scal_k_cc(0,:,:), beta, wdot_k(0,:,:), dx)
      S_cc = S_cc + delta_chi(0,:)
      call extrapolate_cc_to_avg(S_avg, S_cc)
      
      call compute_velocity(vel, S_avg, dx)
      ! compute the advection term
      call compute_advection(advection_k(0,:,:), scal_k_avg(0,:,:), vel, dx)
      
      ! compute the diffusion term
      call compute_diffusion(diffusion_k(0,:,:), scal_k_avg(0,:,:), beta, gamma_face, dx)
      call add_diffdiff_terms(advection_k(0,:,RhoH), scal_k_avg(0,:,:), beta, gamma_face, dx)
      
      do m = 1,nnodes-1
         scal_k_avg(m,:,:)  =  scal_k_avg(m-1,:,:)
         advection_k(m,:,:) = advection_k(m-1,:,:)
         diffusion_k(m,:,:) = diffusion_k(m-1,:,:)
         wdot_k(m,:,:)      =      wdot_k(m-1,:,:)
      end do
      
      call compute_integrals_avg(I_k_avg, advection_k, diffusion_k, wdot_k, dt)
      call compute_integrals_cc(I_k_cc,  advection_k, diffusion_k, wdot_k, dt)
      
      scal_kp1_avg(0,:,:) = scal_k_avg(0,:,:)
      advection_kp1(0,:,:) = advection_k(0,:,:)
      diffusion_kp1(0,:,:) = diffusion_k(0,:,:)
      wdot_kp1(0,:,:) = wdot_k(0,:,:)
            
      do k=1,misdc_iterMax
         do m=0,nnodes-2
            write(*,'(AI2AI2)') 'k = ',k,'   m = ',m
            
            ! convert the cell-averaged scalars to cell centers
            call scal_avg_to_cc(scal_kp1_cc(m,:,:), scal_kp1_avg(m,:,:))
            
            ! compute delta chi prediction
            if (m .ne. 0) then
               call increment_delta_chi(delta_chi(m,:), scal_kp1_cc(m,:,:), dtm(m), dx)
            end if
            
            ! compute the diffusion coefficients
            call compute_diffusion_coefficients(beta, scal_kp1_cc(m,:,:))
            ! compute the divergence constraint, S
            call compute_div_u(S_cc, scal_kp1_avg(m,:,:), scal_kp1_cc(m,:,:), beta, wdot_kp1(m,:,:), dx)
            
            ! add delta chi to the constraint
            S_cc = S_cc + delta_chi(m,:)
            ! convert from cell-centered S to cell-average
            call extrapolate_cc_to_avg(S_avg, S_cc)
            
            ! compute velocity by integrating the constraint
            ! the velocity is computed at faces
            call compute_velocity(vel, S_avg, dx)
            
            ! compute the advection terms A^{m,(k+1)}
            call compute_advection(advection_kp1(m,:,:), scal_kp1_avg(m,:,:), vel, dx)
            call compute_diffusion(diffusion_kp1(m,:,:), scal_kp1_avg(m,:,:), &
                                   beta, gamma_face, dx)
            call add_diffdiff_terms(advection_kp1(m,:,RhoH), scal_kp1_avg(m,:,:), beta, gamma_face, dx)
            
            ! update the density
            call update_density(scal_AD_avg, scal_kp1_avg(m,:,:), &
                                advection_kp1(m,:,:), advection_k(m,:,:), &
                                I_k_avg(m,:,:), dtm(m))
            
            ! compute the iteratively-lagged diffusion coefficients
            call scal_avg_to_cc(scal_k_cc(m+1,:,:), scal_k_avg(m+1,:,:))
            call compute_diffusion_coefficients(beta, scal_k_cc(m+1,:,:))
            
            ! perform the diffusion correction for species
            call species_AD_correction(scal_AD_avg, scal_kp1_avg(m,:,:), beta, &
                                       advection_kp1(m,:,:), advection_k(m,:,:), &
                                       diffusion_k(m+1,:,:), I_k_avg(m,:,:), dtm(m), dx)
            
            ! update the provisional enthalpy
            call enthalpy_AD_correction(scal_AD_avg, scal_kp1_avg(m,:,:), beta, &
                                        advection_kp1(m,:,:), advection_k(m,:,:), &
                                        diffusion_k(m+1,:,:), I_k_avg(m,:,:), dtm(m), dx)
            
            scal_kp1_avg(m+1,:,:) = scal_AD_avg
            scal_kp1_avg(m+1,:,Temp) = scal_k_avg(m+1,:,Temp)
            
            ! compute gamma, conservatively corrected diffusion term
            !call compute_diffusion(diffusion_kp1(m+1,:,:), scal_kp1_avg(m+1,:,:), &
            !                       beta, gamma_face, dx)
            do i=0,nx-1
               diffusion_kp1(m+1,i,:) = (scal_kp1_avg(m+1,i,:) - scal_kp1_avg(m,i,:) - I_k_avg(m,i,:))/dtm(m) & 
                  - (advection_kp1(m,i,:) - advection_k(m,i,:) - diffusion_k(m+1,i,:))
            end do
            
            ! call the chemistry solver to solve the correction equation
            call reaction_correction(scal_kp1_avg(m+1,:,:),  scal_kp1_avg(m,:,:), &
                                     scal_k_cc(m+1,:,:),     k, &
                                     advection_kp1(m,:,:),   advection_k(m,:,:), &
                                     diffusion_kp1(m+1,:,:), diffusion_k(m+1,:,:), &
                                     wdot_k(m+1,:,:),        wdot_kp1(m+1,:,:), &
                                     I_k_avg(m,:,:), I_k_cc(m,:,:), dtm(m))
            call get_temp(scal_kp1_avg(m+1,:,:))
            call scal_avg_to_cc(scal_kp1_cc(m+1,:,:), scal_kp1_avg(m+1,:,:))
         end do
                  
         m = nnodes-1
         call increment_delta_chi(delta_chi(m,:), scal_kp1_cc(m,:,:), dtm(m-1), dx)
         
         ! recompute all the A, D, R terms
         do m=0,nnodes-1
            
            call compute_diffusion_coefficients(beta, scal_kp1_cc(m,:,:))
            call compute_div_u(S_cc, scal_kp1_avg(m,:,:), scal_kp1_cc(m,:,:), beta, wdot_kp1(m,:,:), dx)
            S_cc = S_cc + delta_chi(m,:)
            call extrapolate_cc_to_avg(S_avg, S_cc)
            call compute_velocity(vel, S_avg, dx)
            call compute_advection(advection_kp1(m,:,:), scal_kp1_avg(m,:,:), vel, dx)
            
            ! compute the iteratively-lagged diffusion coefficients
            call compute_diffusion_coefficients(beta, scal_k_cc(m,:,:))
            call compute_diffusion(diffusion_kp1(m,:,:), scal_kp1_avg(m,:,:), &
                                   beta, gamma_face, dx)
            
            call add_diffdiff_terms(advection_kp1(m,:,RhoH), scal_kp1_avg(m,:,:), beta, gamma_face, dx)

         end do
         
         scal_k_avg  = scal_kp1_avg
         advection_k = advection_kp1
         diffusion_k = diffusion_kp1
         wdot_k = wdot_kp1
         
         ! recompute the integral I_k
         call compute_integrals_avg(I_k_avg, advection_k, diffusion_k, wdot_k, dt)
         call compute_integrals_cc(I_k_cc, advection_k, diffusion_k, wdot_k, dt)
      end do
      
      wdot = wdot_kp1(nnodes-1,:,:)
      ! advance the solution
      scal_np1_avg = scal_kp1_avg(nnodes-1,:,:)

      call cpu_time(wallclock_time2)

      print*,'Time to advance time step',wallclock_time2-wallclock_time1

   end subroutine advance
   
end module advance_module
