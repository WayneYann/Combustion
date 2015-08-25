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
   
   subroutine temp_plot(v1,v2,v3,v4,dx,nsteps)
      implicit none
      include 'spec.h'

      integer nsteps
      double precision, dimension(0:nx-1) :: v1,v2,v3,v4
      real*8    dx
      
      character pltfile*(8)
      character char_of_int*(5)
      
      integer i
      pltfile(1:3) = 'plt'
      write(char_of_int,1005) nsteps
      pltfile(4:8) = char_of_int
 1005 format(i5.5)
 1006 FORMAT(200(E23.15E3,1X))

      open(10,file=pltfile,form='formatted')
      print *,'...writing data to ',pltfile

      do i=0,nx-1
         write(10,1006) dx*i, v1(i),v2(i),v3(i),v4(i)
      enddo

      close(10)

   end subroutine temp_plot
   
   subroutine get_temp(scal)
      double precision, intent(inout) :: scal(-2:nx+1, nscal)
      
      double precision :: Y(NSpec)
      double precision :: hmix
      integer :: n, i
      ! chemsolve stuff
      integer NiterMax, Niter
      parameter (NiterMax = 30)
      double precision res(NiterMax), errMax
      
      errMax = hmix_Typ*1.e-20
      
      do i=0,nx-1
         ! compute the temperature from h and Y
         do n = 1,Nspec
            Y(n) = scal(i,FirstSpec+n-1)/scal(i,Density)
         enddo
         hmix = scal(i,RhoH) / scal(i,Density)
         
         ! get the new value for the temperature
         call FORT_TfromHYpt(scal(i,Temp), hmix, Y, &
                             Nspec, errMax, NiterMax, res, Niter)
         if (Niter.lt.0) then
            print *,'vodeF_T_RhoY: H to T solve failed'
            print *,'Niter=',Niter
            stop
         endif
      end do
      call fill_scal_cc_ghost_cells(scal)
   end subroutine get_temp
   
   subroutine advance(scal_n_cc, scal_np1_cc, vel, S_avg, dt, dx)
      implicit none
      ! parameters
      ! two ghost cells
      double precision, intent(in   ) ::   scal_n_cc(-2:nx+1, nscal)
      double precision, intent(out  ) :: scal_np1_cc(-2:nx+1, nscal)
      ! face-value, no ghost cells
      double precision, intent(out  ) ::         vel( 0:nx)
      double precision, intent(out  ) ::       S_avg( 0:nx-1)
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
      double precision :: delta_chi(0:nnodes-1, 0:nx-1)

      ! we also need to know the value of the contraint S
      ! (one ghost cell)
      double precision ::  S_cc(0:nx-1)
      !double precision :: S_avg(0:nx-1)

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
      
      ! time substeps (determined by Gauss-Lobatto rule)
      dtm(0) = 0.5d0*dt
      dtm(1) = 0.5d0*dt
      
      ! initialize everything to zero
      delta_chi = 0
      advection_k = 0
      diffusion_kp1 = 0
      diffusion_k = 0
      diffusion_kp1 = 0
      
      ! initialize the 0th iterate to the initial condition
      scal_k_cc(0,:,:) = scal_n_cc
      call fill_scal_cc_ghost_cells(scal_k_cc(0,:,:))
      
      call fill_cc_ghost_cells(scal_k_cc(0,:,FirstSpec), Y_bc(1,0)*rho_bc(0))
      
      ! compute the advection term
      call increment_delta_chi(delta_chi(0,:), scal_k_cc(0,:,:), dtm(0), dx)
      call compute_diffusion_coefficients(beta, scal_k_cc(0,:,:))
      
      ! need to compute div u to fourth order
      call compute_div_u(S_cc, scal_k_cc(0,:,:), beta, dx)
      S_cc = S_cc + delta_chi(0,:)
      
      call extrapolate_cc_to_avg(S_avg, S_cc)
      call compute_velocity(vel, S_avg, dx)
      ! compute the advection term
      call compute_advection(advection_k(0,:,:), scal_k_cc(0,:,:), vel, dx)
      
      ! compute the diffusion term
      call compute_diffusion(diffusion_k(0,:,:), scal_k_cc(0,:,:), beta, gamma_face, dx)
      call add_diffdiff_terms(advection_k(0,:,RhoH), scal_k_cc(0,:,:), beta, gamma_face, dx)
      
      ! compute the reaction term
      call compute_production_rate(wdot_k(0,:,:), scal_k_cc(0,:,:))
      
      do m = 1,nnodes-1
         scal_k_cc(m,:,:)   =   scal_k_cc(m-1,:,:)
         advection_k(m,:,:) = advection_k(m-1,:,:)
         diffusion_k(m,:,:) = diffusion_k(m-1,:,:)
         wdot_k(m,:,:)      =      wdot_k(m-1,:,:)
      end do
      
      call compute_integrals_avg(I_k_avg, advection_k, diffusion_k, wdot_k, dt)
      call compute_integrals_cc(I_k_cc, advection_k, diffusion_k, wdot_k, dt)
      
      scal_kp1_cc(0,:,:) = scal_k_cc(0,:,:)
      advection_kp1(0,:,:) = advection_k(0,:,:)
      diffusion_kp1(0,:,:) = diffusion_k(0,:,:)
      wdot_kp1(0,:,:) = wdot_k(0,:,:)
      
      do k=1,misdc_iterMax
         do m=0,nnodes-2
            write(*,'(AI2AI2)') 'k = ',k,'   m = ',m
            ! convert the cell-centered scalars to cell-average
            call scal_cc_to_avg(scal_m_avg, scal_kp1_cc(m,:,:))
            
            ! compute delta chi prediction
            if (m .ne. 0) then
               call increment_delta_chi(delta_chi(m,:), scal_kp1_cc(m,:,:), dtm(m), dx)
            end if

            ! compute the diffusion coefficients
            call compute_diffusion_coefficients(beta, scal_kp1_cc(m,:,:))
            ! compute the divergence constraint, S
            call compute_div_u(S_cc, scal_kp1_cc(m,:,:), beta, dx)
            
            ! add delta chi to the constraint
            S_cc = S_cc + delta_chi(m,:)
            ! convert from cell-centered S to cell-average
            call extrapolate_cc_to_avg(S_avg, S_cc)
            ! compute velocity by integrating the constraint
            ! the velocity is computed at faces
            call compute_velocity(vel, S_avg, dx)
            
            ! compute the advection terms A^{m,(k+1)}
            call compute_advection(advection_kp1(m,:,:), scal_kp1_cc(m,:,:), vel, dx)
            call compute_diffusion(diffusion_kp1(m,:,:), scal_kp1_cc(m,:,:), &
                                   beta, gamma_face, dx)
               
            call add_diffdiff_terms(advection_kp1(m,:,RhoH), scal_kp1_cc(m,:,:), beta, gamma_face, dx)
            
            ! update the density
            call update_density(scal_AD_avg, scal_m_avg, &
                                advection_kp1(m,:,:), advection_k(m,:,:), &
                                I_k_avg(m,:,:), dtm(m))
            
            call fill_avg_ghost_cells(scal_m_avg(:,Density), rho_bc(on_lo))
            
            ! compute the iteratively-lagged diffusion coefficients
            call compute_diffusion_coefficients(beta, scal_k_cc(m+1,:,:))
            
            ! perform the diffusion correction for species
            call species_AD_correction(scal_AD_avg, scal_m_avg, beta, &
                                       advection_kp1(m,:,:), advection_k(m,:,:), &
                                       diffusion_k(m+1,:,:), I_k_avg(m,:,:), dtm(m), dx)
            
            ! update the provisional enthalpy
            call enthalpy_AD_correction(scal_AD_avg, scal_m_avg, beta, &
                                        advection_kp1(m,:,:), advection_k(m,:,:), &
                                        diffusion_k(m+1,:,:), I_k_avg(m,:,:), dtm(m), dx)
            
            ! convert provisional solution to cell-centered values
            call scal_avg_to_cc(scal_kp1_cc(m+1,:,:), scal_AD_avg(:,:))
            
            ! recompute the diffusion coefficients
            call compute_diffusion_coefficients(beta, scal_k_cc(m+1,:,:))
            ! compute gamma, conservatively corrected diffusion term
            ! call compute_div_gamma
            call compute_diffusion(diffusion_kp1(m+1,:,:), scal_kp1_cc(m+1,:,:), &
                                   beta, gamma_face, dx)
            
            !call write_plt(vel, scal_kp1_cc(m+1,:,:), S_cc, dx, 2*k + m - 1, dt*m/2.d0)
            
            ! call the chemistry solver to solve the correction equation
            call reaction_correction(scal_kp1_cc(m+1,:,:),   scal_kp1_cc(m,:,:), &
                                     advection_kp1(m,:,:),   advection_k(m,:,:), &
                                     diffusion_kp1(m+1,:,:), diffusion_k(m+1,:,:), &
                                     wdot_k(m+1,:,:), I_k_cc(m,:,:), dtm(m))
            
            !call temp_plot(diffusion_kp1(m+1,:,FirstSpec), scal_AD_avg(:,FirstSpec), &
            !   I_k_cc(m,:,FirstSpec), diffusion_k(m+1,:,FirstSpec), dx, 2*k + m - 1)
            
            !call write_plt(vel, scal_kp1_cc(m+1,:,:), S_cc, dx, 2*k + m - 1, dt*m/2.d0)
            
            ! todo: remove this
            !scal_kp1_cc(m+1,:,Temp) = scal_k_cc(m+1,:,Temp)
            !call get_temp(scal_kp1_cc(m+1,:,:))
            
         end do
         
         m = nnodes-1
         call increment_delta_chi(delta_chi(m,:), scal_kp1_cc(m,:,:), dtm(0), dx)
         
         ! recompute all the A, D, R terms 
         do m=0,nnodes-1
            call compute_diffusion_coefficients(beta, scal_kp1_cc(m,:,:))
            call compute_div_u(S_cc, scal_kp1_cc(m,:,:), beta, dx)
            
            S_cc = S_cc + delta_chi(m,:)
            call extrapolate_cc_to_avg(S_avg, S_cc)
            call compute_velocity(vel, S_avg, dx)
            call compute_diffusion(diffusion_kp1(m,:,:), scal_kp1_cc(m,:,:), &
                                   beta, gamma_face, dx)
            call compute_advection(advection_kp1(m,:,:), scal_kp1_cc(m,:,:), vel, dx)
            call add_diffdiff_terms(advection_kp1(m,:,RhoH), scal_kp1_cc(m,:,:), beta, gamma_face, dx)
            call compute_production_rate(wdot_kp1(m,:,:), scal_kp1_cc(m,:,:))
         end do
         
         scal_k_cc = scal_kp1_cc
         advection_k = advection_kp1
         diffusion_k = diffusion_kp1
         wdot_k = wdot_kp1
         
         ! recompute the integral I_k
         call compute_integrals_avg(I_k_avg, advection_k, diffusion_k, wdot_k, dt)
         call compute_integrals_cc(I_k_cc, advection_k, diffusion_k, wdot_k, dt)
         
      end do
      
      ! advance the solution
      scal_np1_cc = scal_kp1_cc(nnodes-1,:,:)
   end subroutine advance
   
end module advance_module
