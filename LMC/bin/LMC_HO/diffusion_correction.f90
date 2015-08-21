module diffusion_correction_module
   
   use ghost_cells_module
   use cell_conversions_module
   
   implicit none
   private
   include 'spec.h'
   public :: species_AD_correction, enthalpy_AD_correction, get_spec_visc_terms
   public :: compute_diffusion
contains
   
   ! compute the diffusion terms for rho Y_j and rho h
   subroutine compute_diffusion(diffusion, scal_cc, beta, gamma_face, dx)
      implicit none
      double precision, intent(out) ::  diffusion( 0:nx-1,nscal)
      double precision, intent(in ) ::    scal_cc(-2:nx+1,nscal)
      double precision, intent(in ) ::       beta(-2:nx+1,nscal)
      double precision, intent(out) :: gamma_face( 0:nx,  Nspec)

      double precision, intent(in ) :: dx
      
      ! compute the diffusion term for rho h
      call get_rhoh_visc_terms(diffusion(:,RhoH), scal_cc, beta, dx)
      ! compute the diffusion term for rho Y_j
      call get_spec_visc_terms(diffusion(:,FirstSpec:), scal_cc, beta, gamma_face, dx)
      
      ! temp: turn off diffusion:
      ! diffusion = 0
   end subroutine compute_diffusion
   
   subroutine get_rhoh_visc_terms(visc, scal, beta, dx)
      implicit none
      double precision, intent(out) :: visc( 0:nx-1)
      double precision, intent(in ) :: scal(-2:nx+1,nscal)
      double precision, intent(in ) :: beta(-2:nx+1,nscal)
      
      double precision, intent(in ) :: dx
      
      integer :: i
      
      double precision ::         h(-2:nx+1)
      double precision :: beta_face( 0:nx)
      double precision ::    grad_h( 0:nx)
      
      ! compute h at cell centers by dividing by rho
      do i=-2,nx+1
         h(i) = scal(i,RhoH)/scal(i,Density)
      end do
      
      ! obtain the gradient at faces
      call cc_to_grad(grad_h, h, dx)
      ! obtain the diffusion coefficients at faces
      call cc_to_face(beta_face, beta(:,RhoH))
      
      ! compute the viscous term
      do i=0,nx-1
         visc(i) = (beta_face(i+1)*grad_h(i+1) - beta_face(i)*grad_h(i))/dx
      end do
      
   end subroutine get_rhoh_visc_terms
   
   subroutine get_spec_visc_terms(visc_avg, scal_cc, beta_cc,gamma_face,dx)
      implicit none
      include 'spec.h'
      double precision, intent(out) ::       visc_avg( 0:nx-1,Nspec)
      double precision, intent(in ) ::       scal_cc(-2:nx+1,nscal)
      double precision, intent(in ) ::       beta_cc(-2:nx+1,nscal)
      double precision, intent(out) :: gamma_face( 0:nx,  Nspec)
      double precision, intent(in ) :: dx
      
      integer i,n,is
      double precision :: beta_face( 0:nx)
      double precision ::    grad_Y( 0:nx)
      double precision ::    Y_face( 0:nx,  Nspec)
      double precision ::         Y(-2:nx+1,Nspec)
      
      double precision :: sum_g
      double precision :: sum_y
      
      ! this subroutine computes the term div(rho D_m grad Y_m)
      
      ! get cell-centered values for Y_n by dividing by rho
      do n=1,Nspec
         do i=-2,nx+1
            Y(i,n) = scal_cc(i,FirstSpec+n-1)/scal_cc(i,Density)
         end do
         call fill_cc_ghost_cells(Y(:,n), Y_bc(n,on_lo))
      end do
      
      do n=1,Nspec
         is = FirstSpec + n - 1
         ! obtain the diffusion coefficient at faces
         call cc_to_face(beta_face, beta_cc(:, is))
         ! obtain rho Y_n at faces
         call cc_to_face(Y_face(:,n), Y(:,n))
         Y_face(0,n) = Y_bc(n,on_lo)
         ! compute the gradient of Y_n at faces
         call cc_to_grad(grad_Y, Y(:,n), dx)
         ! compute the flux, gamma
         do i=0,nx
            gamma_face(i,n) = beta_face(i)*grad_Y(i)
         end do
      end do
      
      ! correct the fluxes so they add up to zero before computing visc
      ! we correct according to a weighting of rho Y_n / rho
      ! where rho is taken to be the sum of rho Y_j on faces
      do i=0,nx
         sum_g = 0
         sum_y = 0
         do n=1,Nspec
            sum_g = sum_g + gamma_face(i,n)
            sum_y = sum_y + Y_face(i,n)
         end do
         do n=1,Nspec
            gamma_face(i,n) = gamma_face(i,n) - sum_g*Y_face(i,n)/sum_y
         end do
      end do
      
      ! compute the viscous term according to the divergence theorem
      do i=0,nx-1
         do n=1,Nspec
            visc_avg(i,n) = (gamma_face(i+1,n)-gamma_face(i,n))/dx
         end do
      end do
   end subroutine get_spec_visc_terms
   
   ! obtain the provision solution rho Y_AD
   ! this is just a wrapper for implicit_AD_solve
   subroutine species_AD_correction(scal_AD_avg, scal_m_avg, beta, advection_kp1, &
                                    advection_k, diffusion_k, I_k, dtm, dx)
      implicit none
      double precision, intent(inout) ::   scal_AD_avg(-2:nx+1,nscal)
      double precision, intent(in   ) ::    scal_m_avg(-2:nx+1,nscal)
      double precision, intent(in   ) ::          beta(-2:nx+1,nscal)
      double precision, intent(in   ) :: advection_kp1( 0:nx-1,nscal)
      double precision, intent(in   ) ::   advection_k( 0:nx-1,nscal)
      double precision, intent(in   ) ::   diffusion_k( 0:nx-1,  nscal)
      double precision, intent(in   ) ::           I_k( 0:nx-1,nscal)
      double precision, intent(in   ) ::       dtm, dx
      
      integer :: n, ispec
      
      do n=1,Nspec
         ispec = FirstSpec+n-1
         
         call implicit_AD_solve(scal_AD_avg(:, ispec), scal_m_avg(:, ispec), &
            scal_AD_avg(:,Density), beta(:, ispec), &
            advection_kp1(:, ispec), advection_k(:, ispec), &
            diffusion_k(:, ispec), I_k(:, ispec), dtm, dx, Y_bc(n,on_lo))
      end do
   end subroutine species_AD_correction

   ! obtain the provisional solution rho h_AD
   subroutine enthalpy_AD_correction(scal_AD_avg, scal_m_avg, beta, advection_kp1, &
                                       advection_k, diffusion_k, I_k, dtm, dx)
      implicit none
      double precision, intent(inout) ::   scal_AD_avg(-2:nx+1,nscal)
      double precision, intent(in   ) ::    scal_m_avg(-2:nx+1,nscal)
      double precision, intent(in   ) ::          beta(-2:nx+1,nscal)
      double precision, intent(in   ) :: advection_kp1( 0:nx-1,nscal)
      double precision, intent(in   ) ::   advection_k( 0:nx-1,nscal)
      double precision, intent(in   ) ::   diffusion_k( 0:nx-1,nscal)
      double precision, intent(in   ) ::           I_k( 0:nx-1,nscal)
      double precision, intent(in   ) ::       dtm, dx
      
      call implicit_AD_solve(scal_AD_avg(:, RhoH), scal_m_avg(:, RhoH), &
         scal_AD_avg(:,Density), beta(:, RhoH), &
         advection_kp1(:, RhoH), advection_k(:, RhoH), &
         diffusion_k(:, RhoH), I_k(:, RhoH), dtm, dx, h_bc(on_lo))
   end subroutine enthalpy_AD_correction
   
   ! this subroutine performs the implicit diffusion solve for a scalar quantity 'phi'
   ! in practice, phi is either a mass fraction, Y_j, or enthalpy, h.
   ! the result returned is the cell average of rho times phi
   subroutine implicit_AD_solve(rhophi_AD_avg, rhophi_m_avg, rho_mp1_avg, beta, &
                                advection_kp1, advection_k, diffusion_k, &
                                I_k, dtm, dx, phi_bdry)
      implicit none
      double precision, intent(out  ) :: rhophi_AD_avg(-2:nx+1)
      double precision, intent(in   ) ::  rhophi_m_avg(-2:nx+1)
      double precision, intent(in   ) ::   rho_mp1_avg(-2:nx+1)
      double precision, intent(in   ) ::          beta(-2:nx+1)
      double precision, intent(in   ) :: advection_kp1(0:nx-1)
      double precision, intent(in   ) ::   advection_k(0:nx-1)
      double precision, intent(in   ) ::   diffusion_k(0:nx-1)
      double precision, intent(in   ) ::           I_k(0:nx-1)
      double precision, intent(in   ) :: dtm, dx, phi_bdry
      
      ! the matrix system we solve is 'almost' pentadiagonal
      ! the first and last rows have a bandwidth of 3
      ! and therefore the matrix is technically 'heptadiagonal'
      ! 3 bands above, and 3 bands below the main diagonal
      integer, parameter :: ml = 3, mu = 3, lda = 2*ml + mu + 1, d = ml+mu+1
      
      double precision :: rhs(0:nx-1)
      double precision :: rho_G
      double precision :: abd(lda, 0:nx-1)
      double precision :: beta_face(0:nx)
      double precision :: phi(-2:nx+1)
      double precision :: dxsq
      
      integer :: i, info, ipvt(nx)
      
      ! construct the right hand side
      do i=0,nx-1
         rhs(i) = rhophi_m_avg(i) &
            + dtm*(advection_kp1(i) - advection_k(i) - diffusion_k(i)) + I_k(i)
      end do
      
      ! get the diffusion coefficients as face values rather than cell centers
      call cc_to_face(beta_face, beta)
      
      dxsq = dx*dx
      
      abd = 0
      ! construct the matrix we need to invert
      ! note that the matrix is stored in 'banded' form, which mean
      ! abd(i,j) is the entry in the jth column of the matrix, on the ith diagonal
      ! (where 1 is the uppermost diagonal, and 1+ml+mu us the bottommost)
      do i=0,nx-1
         rho_G = (5*rho_mp1_avg(i-2) - 34*rho_mp1_avg(i-1) &
                + 34*rho_mp1_avg(i+1) - 5*rho_mp1_avg(i+2))/27648.d0
         
         ! here we construct the banded matrix ("almost pentadiagonal..")
         ! these terms come from the product rule for cell averages
         ! we also add the term that comes from the diffusion operator
         ! computed using divergence theorem and the 4th order gradient stencil
         if (i .ge. 2) abd(d+2,i-2) = 5*rho_G + dtm*beta_face(i)/(12*dxsq)
         if (i .ge. 1) abd(d+1,i-1) = -34*rho_G - dtm*(beta_face(i+1) + 15*beta_face(i))/(12*dxsq)
         abd(d, i) = rho_mp1_avg(i) + dtm*(15*beta_face(i+1) + 15*beta_face(i))/(12*dxsq)
         if (i .le. nx-2) abd(d-1,i+1) = 34*rho_G - dtm*(15*beta_face(i+1) + beta_face(i))/(12*dxsq)
         if (i .le. nx-3) abd(d-2,i+2) = -5*rho_G + dtm*beta_face(i+1)/(12*dxsq)
      end do
      
      ! take care of the boundary condition/ghost cells
      ! the inflow ghost cells include the Dirichlet boundary condition, 
      ! whose term is added to the right-hand side
      ! all the remaining terms involve cell-averages in the 'valid region'
      ! and are added to the corresponding entries in the matrix
      
      ! the necessary terms to add are computed simply by substituting 
      ! the formulas for the cell-average ghost cells
      
      ! inflow:
      ! i = 0
      rho_G = (5*rho_mp1_avg(-2) - 34*rho_mp1_avg(-1) &
            + 34*rho_mp1_avg(1) - 5*rho_mp1_avg(2))/27648.d0
      
      rhs(0)     = rhs(0) + phi_bdry*(dtm*5*(10*beta_face(0) + beta_face(1))/(12.d0*dxsq) + 45*rho_G)
      abd(d,  0) = abd(d,  0) + dtm*(650*beta_face(0) + 77*beta_face(1))/(144.d0*dxsq) + 31*rho_G/4.d0
      abd(d-1,1) = abd(d-1,1) - dtm*(310*beta_face(0) + 43*beta_face(1))/(144.d0*dxsq) + 71*rho_G/4.d0
      abd(d-2,2) = abd(d-2,2) + dtm*(110*beta_face(0) + 17*beta_face(1))/(144.d0*dxsq) - 49*rho_G/4.d0
      abd(d-3,3) = abd(d-3,3) - dtm*(  6*beta_face(0) +    beta_face(1))/(48.d0*dxsq) + 11*rho_G/4.d0
      
      ! i = 1
      rho_G = (5*rho_mp1_avg(-1) - 34*rho_mp1_avg(0) &
          + 34*rho_mp1_avg(2) - 5*rho_mp1_avg(3))/27648.d0
      
      rhs(1)     = rhs(1) - phi_bdry*(dtm*5*beta_face(1)/(12.d0*dxsq) + 25*rho_G)
      abd(d+1,0) = abd(d+1,0) - dtm*beta_face(1)*77/(144.d0*dxsq) - 385*rho_G/12.d0
      abd(d,  1) = abd(d,  1) + dtm*beta_face(1)*43/(144.d0*dxsq) + 215*rho_G/12.d0
      abd(d-1,2) = abd(d-1,2) - dtm*beta_face(1)*17/(144.d0*dxsq) - 85*rho_G/12.d0 
      abd(d-2,3) = abd(d-2,3) + dtm*beta_face(1)/(48.d0*dxsq)  + 5*rho_G/4.d0
      
      ! for outflow, the ghost cells are filled to satisfy the Neumann condition, and 
      ! therefore we do not need to add anything to the right-hand side
      ! outflow:
      ! i = nx-2
      rho_G = (5*rho_mp1_avg(nx-4) - 34*rho_mp1_avg(nx-3) &
           + 34*rho_mp1_avg(nx-1) - 5*rho_mp1_avg(nx))/27648.d0
      
      abd(d+2,nx-4) = abd(d+2,nx-4) + dtm*beta_face(nx-1)/(120.d0*dxsq) - rho_G/2.d0
      abd(d+1,nx-3) = abd(d+1,nx-3) - dtm*beta_face(nx-1)/(24.d0*dxsq) + 5*rho_G/2.d0
      abd(d,  nx-2) = abd(d,  nx-2) + dtm*3*beta_face(nx-1)/(40.d0*dxsq) - 9*rho_G/2.d0
      abd(d-1,nx-1) = abd(d-1,nx-1) + dtm*beta_face(nx-1)/(24.d0*dxsq) - 5*rho_G/2.d0
      
      ! i = nx-1
      rho_G = (5*rho_mp1_avg(nx-3) - 34*rho_mp1_avg(nx-2) &
            + 34*rho_mp1_avg(nx) - 5*rho_mp1_avg(nx+1))/27648.d0
      
      abd(d+3,nx-4) = abd(d+3,nx-4) - dtm*beta_face(nx-1)/(120.d0*dxsq) - 41*rho_G/10.d0
      abd(d+2,nx-3) = abd(d+2,nx-3) + dtm*beta_face(nx-1)/(24.d0*dxsq) + 41*rho_G/2.d0
      abd(d+1,nx-2) = abd(d+1,nx-2) - dtm*(9*beta_face(nx-1) - 10*beta_face(nx))/(120.d0*dxsq) - 419*rho_G/10.d0
      abd(d,  nx-1) = abd(d,  nx-1) - dtm*(beta_face(nx-1) + 30*beta_face(nx))/(24.d0*dxsq) + 109*rho_G/2.d0
      
      ! perform the banded linear solve for phi
      ! call linpack to do the factorization
      call dgbfa(abd, lda, nx, ml, mu, ipvt, info)
      ! call linpack to do the solve -- the solution is returned in rhs
      call dgbsl(abd, lda, nx, ml, mu, ipvt, rhs, 0)
      
      do i=0,nx-1
         phi(i) = rhs(i)
      end do
      
      ! take the solution obtained from linpack
!      do i=0,nx-1
!         phi(i) = rhophi_m_avg(i) &
!            + dtm*(advection_kp1(i) - advection_k(i) - diffusion_k(i)) + I_k(i)
!      end do
!      
!      ! temporary!!!
!      rhophi_AD_avg = phi
!      call fill_avg_ghost_cells(rhophi_AD_avg, phi_bdry*rho_bc(on_lo))
      
      ! fill in the ghost cells
      call fill_avg_ghost_cells(phi, phi_bdry)
      
      ! obtain average of product of phi and rho (and fill in the ghost cells)
      call mult_avgs_bdry(rhophi_AD_avg, phi, rho_mp1_avg, rho_bc(on_lo)*phi_bdry)
   end subroutine implicit_AD_solve

end module diffusion_correction_module
