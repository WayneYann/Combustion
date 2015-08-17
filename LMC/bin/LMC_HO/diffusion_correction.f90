module diffusion_correction_module
   
   use ghost_cells_module
   use cell_conversions_module
   
   implicit none
   private
   include 'spec.h'
   public :: species_AD_correction, enthalpy_AD_correction, get_spec_visc_terms
   public :: compute_diffusion
contains
   
   subroutine compute_diffusion(diffusion, scal_cc, beta, gamma_lo, gamma_hi, dx)
      implicit none
      double precision, intent(out) :: diffusion( 0:nx-1,nscal)
      double precision, intent(in ) ::   scal_cc(-2:nx+1,nscal)
      double precision, intent(in ) ::      beta(-2:nx+1,nscal)
      double precision, intent(out) ::  gamma_lo( 0:nx-1,Nspec)
      double precision, intent(out) ::  gamma_hi( 0:nx-1,Nspec)

      double precision, intent(in ) :: dx
      
      call get_rhoh_visc_terms(diffusion(:,RhoH), scal_cc, beta, dx)
      call get_spec_visc_terms(diffusion(:,FirstSpec:), scal_cc, beta, gamma_lo, gamma_hi, dx)
   end subroutine compute_diffusion
   
   subroutine get_rhoh_visc_terms(visc, scal, beta, dx)
      implicit none
      double precision, intent(out) :: visc( 0:nx-1)
      double precision, intent(in ) :: scal(-2:nx+1,nscal)
      double precision, intent(in ) :: beta(-2:nx+1,nscal)
      
      double precision, intent(in ) :: dx
      
      integer :: i
      double precision :: beta_lo,beta_hi
      double precision :: flux_lo,flux_hi
      double precision :: dxsqinv, h_lo, h_mid, h_hi

      dxsqinv = 1.d0/(dx*dx)

      do i=0,nx-1
         ! todo: make this fourth order
         beta_lo = 0.5d0*(beta(i,RhoH) + beta(i-1,RhoH))
         beta_hi = 0.5d0*(beta(i,RhoH) + beta(i+1,RhoH))

         h_hi  = scal(i+1,RhoH) / scal(i+1,Density)
         h_mid = scal(i  ,RhoH) / scal(i  ,Density)
         h_lo  = scal(i-1,RhoH) / scal(i-1,Density)

         flux_hi = beta_hi*(h_hi - h_mid)
         flux_lo = beta_lo*(h_mid - h_lo)
         visc(i) = (flux_hi - flux_lo) * dxsqinv 
      enddo
   end subroutine get_rhoh_visc_terms
   
   subroutine get_spec_visc_terms(visc, scal, beta,gamma_lo,gamma_hi,dx)
      implicit none
      include 'spec.h'
      double precision, intent(out) ::     visc( 0:nx-1,Nspec)
      double precision, intent(in ) ::     scal(-2:nx+1,nscal)
      double precision, intent(in ) ::     beta(-2:nx+1,nscal)
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
            
            ! todo: fourth order avg to face of diffusion coefficient
            beta_lo = 0.5d0*(beta(i,is) + beta(i-1,is))
            beta_hi = 0.5d0*(beta(i,is) + beta(i+1,is))
            
            ! todo: fourth order gradient of species
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
         
         call fill_avg_ghost_cells(scal_AD_avg(:,ispec), Y_bc(n,on_lo)*rho_bc(on_lo))
      end do
      
   end subroutine species_AD_correction

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
      
      call fill_avg_ghost_cells(scal_AD_avg(:,RhoH), h_bc(on_lo)*rho_bc(on_lo))
      
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
                + 34*rho_mp1_avg(i+1) - 5*rho_mp1_avg(i+2))/27648.0
         
         ! here we construct the banded matrix ("almost pentadiagonal..")
         ! these terms come from the product rule for cell averages
         ! we also add the term that comes from the diffusion operator
         ! computed using divergence theorem and the 4th order gradient stencil
         if (i .ge. 2) abd(d+2,i-2) = 5*rho_G + dtm*beta_face(i)/(12.0*dxsq)
         if (i .ge. 1) abd(d+1,i-1) = -34*rho_G - dtm*(beta_face(i+1) + 15*beta_face(i))/(12.0*dxsq)
         abd(d, i) = rho_mp1_avg(i) + dtm*(15*beta_face(i+1) + 15*beta_face(i))/(12.0*dxsq)
         if (i .le. nx-2) abd(d-1,i+1) = 34*rho_G - dtm*(15*beta_face(i+1) + beta_face(i))/(12.0*dxsq)
         if (i .le. nx-3) abd(d-2,i+2) = -5*rho_G + dtm*beta_face(i+1)/(12.0*dxsq)
      end do
      
      ! take care of the boundary condition/ghost cells
      ! the inflow ghost cells include the Dirichlet boundary condition, 
      ! whose term is added to the right-hand side
      ! all the remaining terms involve cell-averages in the 'valid region'
      ! and are added to the corresponding entries in the matrix
      
      ! inflow:
      ! i = 0
      rho_G = (5*rho_mp1_avg(-2) - 34*rho_mp1_avg(-1) &
            + 34*rho_mp1_avg(1) - 5*rho_mp1_avg(2))/27648.0
      
      rhs(0)     = rhs(0) + phi_bdry*(dtm*5*(10*beta_face(0) + beta_face(1))/(12.0*dxsq) + 45*rho_G)
      abd(d,  0) = abd(d,  0) + dtm*(650*beta_face(0) + 77*beta_face(1))/(144.0*dxsq) + 31*rho_G/4.0
      abd(d-1,1) = abd(d-1,1) - dtm*(310*beta_face(0) + 43*beta_face(1))/(144.0*dxsq) + 71*rho_G/4.0
      abd(d-2,2) = abd(d-2,2) + dtm*(110*beta_face(0) + 17*beta_face(1))/(144.0*dxsq) - 49*rho_G/4.0
      abd(d-3,3) = abd(d-3,3) - dtm*(  6*beta_face(0) +    beta_face(1))/(48.0*dxsq) + 11*rho_G/4.0
      
      ! i = 1
      rho_G = (5*rho_mp1_avg(-1) - 34*rho_mp1_avg(0) &
          + 34*rho_mp1_avg(2) - 5*rho_mp1_avg(3))/27648.0
      
      rhs(1)     = rhs(1) - phi_bdry*(dtm*5*beta_face(1)/(12.0*dxsq) + 25*rho_G)
      abd(d+1,0) = abd(d+1,0) - dtm*beta_face(1)*77/(144.0*dxsq) - 385*rho_G/12.0
      abd(d,  1) = abd(d,  1) + dtm*beta_face(1)*43/(144.0*dxsq) + 215*rho_G/12.0
      abd(d-1,2) = abd(d-1,2) - dtm*beta_face(1)*17/(144.0*dxsq) - 85*rho_G/12.0 
      abd(d-2,3) = abd(d-2,3) + dtm*beta_face(1)/(48.0*dxsq)  + 5*rho_G/4.0
      
      ! for outflow, the ghost cells are filled to satisfy the Neumann condition, and 
      ! therefore we do not need to add anything to the right-hand side
      ! outflow:
      ! i = nx-2
      rho_G = (5*rho_mp1_avg(nx-4) - 34*rho_mp1_avg(nx-3) &
           + 34*rho_mp1_avg(nx-1) - 5*rho_mp1_avg(nx))/27648.0
      
      abd(d+2,nx-4) = abd(d+2,nx-4) + dtm*beta_face(nx-1)/(120.0*dxsq) - rho_G/2.0
      abd(d+1,nx-3) = abd(d+1,nx-3) - dtm*beta_face(nx-1)/(24.0*dxsq) + 5*rho_G/2.0
      abd(d,  nx-2) = abd(d,  nx-2) + dtm*3*beta_face(nx-1)/(40.0*dxsq) - 9*rho_G/2.0
      abd(d-1,nx-1) = abd(d-1,nx-1) + dtm*beta_face(nx-1)/(24.0*dxsq) - 5*rho_G/2.0
      
      ! i = nx-1
      rho_G = (5*rho_mp1_avg(nx-3) - 34*rho_mp1_avg(nx-2) &
            + 34*rho_mp1_avg(nx) - 5*rho_mp1_avg(nx+1))/27648.0
      
      abd(d+3,nx-4) = abd(d+3,nx-4) - dtm*beta_face(nx-1)/(120.0*dxsq) - 41*rho_G/10.0
      abd(d+2,nx-3) = abd(d+2,nx-3) + dtm*beta_face(nx-1)/(24.0*dxsq) + 41*rho_G/2.0
      abd(d+1,nx-2) = abd(d+1,nx-2) - dtm*(9*beta_face(nx-1) - 10*beta_face(nx))/(120.0*dxsq) - 419*rho_G/10.0
      abd(d,  nx-1) = abd(d,  nx-1) - dtm*(beta_face(nx-1) + 30*beta_face(nx))/(24.0*dxsq) + 109*rho_G/2.0
      
      ! perform the banded linear solve for phi
      ! call linpack to do the factorization
      call dgbfa(abd, lda, nx, ml, mu, ipvt, info)
      ! call linpack to do the solve -- the solution is returned in rhs
      call dgbsl(abd, lda, nx, ml, mu, ipvt, rhs, 0)
      
      do i=0,nx-1
         phi(i) = rhs(i)
      end do
      
      ! fill in the ghost cells
      call fill_avg_ghost_cells(phi, phi_bdry)
      
      ! obtain average of product of phi and rho
      call mult_avgs(rhophi_AD_avg, phi, rho_mp1_avg, rho_bc(on_lo)*phi_bdry)
      
   end subroutine implicit_AD_solve
   
end module diffusion_correction_module
