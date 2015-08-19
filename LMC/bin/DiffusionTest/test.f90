program test_it
   double precision :: error1, error2
   
   call get_error(32, error1)
   call get_error(64, error2)
   
   print *,log(error2/error1)
   print *,''
   
   call get_error(64, error1)
   call get_error(128, error2)
   
   print *,log(error2/error1)
   print *,''
   
   call get_error(128, error1)
   call get_error(256, error2)
   
   print *,log(error2/error1)
   print *,''
   
contains
   subroutine get_error(nx, error)
      integer, intent(in )          :: nx
      double precision, intent(out) :: error
      
      double precision :: rhs(0:nx-1)
      double precision :: rho(-2:nx+1)
      double precision :: phi(-2:nx+1)
      double precision :: dx
      double precision :: exact(-2:nx+1)
      
      double precision :: endpt, pi
      
      integer :: i
      
      pi = 4.d0*datan(1.d0)
      endpt = pi
      
      dx = endpt/nx
      
      do i=0,nx-1
         rhs(i) = (sin(dx*(i+1))-sin(dx*i))/dx
         exact(i) = (sin(dx*(i+1))-sin(dx*i))/dx
         rho(i) = 0.d0
      end do

      call fill_avg_ghost_cells(rho, 0.d0, nx)
      
      call implicit_AD_solve(phi, rho, rhs, dx, 1.d0, nx)

      error = 0.d0
      
      do i=0,nx-1
         !error = max(error, abs(exact(i)-phi(i)))
         error = error + abs(exact(i) - phi(i))*dx
         
         !print *,abs(exact(i) - phi(i))/dx
         !print *,dx*(i+0.5d0),exact(i),phi(i)
         
      end do
      print *,error
   end subroutine get_error
   
   ! fill cell-averaged ghost cells (Dirichlet inflow, Neumann outflow)
   subroutine fill_avg_ghost_cells(avg, bdry, nx)
      double precision, intent(inout) :: avg(-2:nx+1)
      double precision, intent(in   ) :: bdry
      integer,          intent(in   ) :: nx
      
      avg(-1) = (60*bdry - 77*avg(0) + 43*avg(1) - 17*avg(2) + 3*avg(3))/12.d0
      avg(-2) = (300*bdry - 505*avg(0) + 335*avg(1) - 145*avg(2) + 27*avg(3))/12.d0
      
      avg(nx) = (5*avg(nx-1) + 9*avg(nx-2) - 5*avg(nx-3) + avg(nx-4))/10.d0
      avg(nx+1) = (-15*avg(nx-1) + 29*avg(nx-2) - 15*avg(nx-3) + 3*avg(nx-4))/2.d0
   end subroutine fill_avg_ghost_cells
   
   subroutine implicit_AD_solve(rhophi_AD_avg, rho_mp1_avg, &
                                rhs, dx, phi_bdry, nx)
      implicit none
      integer,          intent(in   ) :: nx
      double precision, intent(out  ) :: rhophi_AD_avg(-2:nx+1)
      double precision, intent(in   ) ::   rho_mp1_avg(-2:nx+1)
      double precision, intent(inout) ::           rhs(0:nx-1)
      double precision, intent(in   ) :: dx, phi_bdry
      
      integer, parameter :: ml = 3, mu = 3, lda = 2*ml + mu + 1, d = ml+mu+1
      
      double precision :: rho_G
      double precision :: abd(lda, 0:nx-1)
      double precision :: beta_face(0:nx)
      double precision :: dxsq, dtm
      
      integer :: i, info, ipvt(nx)
      
      beta_face = 1
      
      dxsq = dx*dx
      dtm = 1.d0
      
      abd = 0
      ! construct the matrix we need to invert
      ! note that the matrix is stored in 'banded' form, which mean
      ! abd(i,j) is the entry in the jth column of the matrix, on the ith diagonal
      ! (where 1 is the uppermost diagonal, and 1+ml+mu us the bottommost)
      do i=0,nx-1
         ! here we construct the banded matrix ("almost pentadiagonal..")
         ! these terms come from the product rule for cell averages
         ! we also add the term that comes from the diffusion operator
         ! computed using divergence theorem and the 4th order gradient stencil
         if (i .ge. 2) abd(d+2,i-2) = 1.d0/(12.d0*dxsq)
         if (i .ge. 1) abd(d+1,i-1) = - 4.d0/(3.d0*dxsq)
         abd(d, i) = 5.d0/(2.d0*dxsq)
         if (i .le. nx-2) abd(d-1,i+1) = - 4.d0/(3.d0*dxsq)
         if (i .le. nx-3) abd(d-2,i+2) = 1.d0/(12.d0*dxsq)
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
!      abd(d+2,nx-4) = abd(d+2,nx-4) - 1.d0/(120.d0*dxsq)
!      abd(d+1,nx-3) = abd(d+1,nx-3) + 1.d0/(24.d0*dxsq)
!      abd(d,  nx-2) = abd(d,  nx-2) - 3.d0/(40.d0*dxsq)
!      abd(d-1,nx-1) = abd(d-1,nx-1) - 1.d0/(24.d0*dxsq)
!      
!      ! i = nx-1
!      abd(d+3,nx-4) = abd(d+3,nx-4) + 1.d0/(120.d0*dxsq)
!      abd(d+2,nx-3) = abd(d+2,nx-3) - 1.d0/(24.d0*dxsq)
!      abd(d+1,nx-2) = abd(d+1,nx-2) - 1.d0/(120.d0*dxsq)
!      abd(d,  nx-1) = abd(d,  nx-1) + 31.d0/(24.d0*dxsq)
      
            ! for outflow, the ghost cells are filled to satisfy the Neumann condition, and 
      ! therefore we do not need to add anything to the right-hand side
      ! outflow:
      ! i = nx-2
      abd(d+2,nx-4) = abd(d+2,nx-4) + beta_face(nx-1)/(120.d0*dxsq)
      abd(d+1,nx-3) = abd(d+1,nx-3) - beta_face(nx-1)/(24.d0*dxsq)
      abd(d,  nx-2) = abd(d,  nx-2) + 3.d0*beta_face(nx-1)/(40.d0*dxsq)
      abd(d-1,nx-1) = abd(d-1,nx-1) + beta_face(nx-1)/(24.d0*dxsq)
      
      ! i = nx-1
      abd(d+3,nx-4) = abd(d+3,nx-4) - beta_face(nx-1)/(120.d0*dxsq)
      abd(d+2,nx-3) = abd(d+2,nx-3) + beta_face(nx-1)/(24.d0*dxsq)
      abd(d+1,nx-2) = abd(d+1,nx-2) - (9*beta_face(nx-1) - 10*beta_face(nx))/(120.d0*dxsq)
      abd(d,  nx-1) = abd(d,  nx-1) - (beta_face(nx-1) + 30*beta_face(nx))/(24.d0*dxsq)
      
      
      ! perform the banded linear solve for phi
      ! call linpack to do the factorization
      call dgbfa(abd, lda, nx, ml, mu, ipvt, info)
      ! call linpack to do the solve -- the solution is returned in rhs
      call dgbsl(abd, lda, nx, ml, mu, ipvt, rhs, 0)
      
      ! take the solution obtained from linpack
      do i=0,nx-1
         rhophi_AD_avg(i) = rhs(i)
      end do
      
      !call fill_avg_ghost_cells(rhophi_AD_avg, phi_bdry, nx)
   end subroutine implicit_AD_solve   
end program test_it
