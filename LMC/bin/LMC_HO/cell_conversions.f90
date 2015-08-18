module cell_conversions_module
   
   use ghost_cells_module
   
   implicit none
   private
   
   include 'spec.h'
   
   public :: cc_to_avg, extrapolate_cc_to_avg, avg_to_cc, scal_cc_to_avg
   public :: cc_to_face, mult_avgs, scal_avg_to_cc, extrapolate_avg_to_cc
   public :: cc_to_grad

contains

   ! convert cell-centered quantity with 2 ghost cells
   ! to a cell-averaged quantity with 2 ghost cells
   ! filling the ghostcells according to Dirichlet inflow
   ! and Neumann outflow boundary conditions
   subroutine cc_to_avg(avg, cc, bdry)
      implicit none
      double precision, intent(out) ::  avg(-2:nx+1)
      double precision, intent(in ) ::   cc(-2:nx+1)
      double precision, intent(in ) :: bdry
      
      integer :: i
      
      do i=0,nx-1
         avg(i) = cc(i) + (cc(i-1) - 2.0*cc(i) + cc(i+1))/24.0
      end do
      
      call fill_avg_ghost_cells(avg, bdry)
   end subroutine cc_to_avg
   
   ! given a cell-centered quantity with no ghost cells, compute cell averages
   ! in the valid region, using a fourth-order extrapolation, and the 
   ! conversion stencil
   subroutine extrapolate_cc_to_avg(avg, cc)
      implicit none
      double precision, intent(out) :: avg(0:nx-1)
      double precision, intent(in ) ::  cc(0:nx-1)
      
      integer :: i
      
      ! one-sided stencils
      avg(0) = (27*cc(0) - 9*cc(1) + 10*cc(2) - 5*cc(3) + cc(4))/24.0
      avg(nx-1) = (501*cc(nx-1) + 31*cc(nx-2) - 5*cc(nx-3) + cc(nx-4))/528.0
      
      !avg(nx-1) = (cc(nx-5) - 5.0*cc(nx-4) + 10.0*cc(nx-3) &
      !           - 9.0*cc(nx-2) + 27.0*cc(nx-2))/24.0
      
      do i=1,nx-2
         avg(i) = cc(i) + (cc(i-1) - 2.0*cc(i) + cc(i+1))/24.0
      end do
   end subroutine extrapolate_cc_to_avg
   
   subroutine extrapolate_avg_to_cc(cc, avg)
      implicit none
      double precision, intent(out) ::  cc(0:nx-1)
      double precision, intent(in ) :: avg(0:nx-1)
      
      integer :: i
      
      cc(0) = (21*avg(0) + 9*avg(1) - 10*avg(2) + 5*avg(3) - avg(4))/24.0
      ! use the outflow condition
      cc(nx-1) = (255*avg(nx-1) - 19*avg(nx-2) + 5*avg(nx-3) - avg(nx-4))/240.0
      
      do i=1,nx-2
         cc(i) = avg(i) - (avg(i-1) - 2*avg(i) + avg(i+1))/24.0
      end do
   end subroutine
   
   ! convert from cell-averaged quantities to cell-centered quantities 
   ! and fill in the ghost cells according to the boundary condition
   subroutine avg_to_cc(cc, avg, bdry)
      implicit none
      double precision, intent(out) ::   cc(-2:nx+1)
      double precision, intent(in ) ::  avg(-2:nx+1)
      double precision, intent(in ) :: bdry
      
      integer :: i
      
      do i=0,nx-1
         cc(i) = avg(i) - (avg(i-1) - 2.0*avg(i) + avg(i+1))/24.0
      end do
      
      call fill_cc_ghost_cells(cc, bdry)
   end subroutine avg_to_cc
   
   ! convert the array of scalar quantities from cell-centered to cell-averaged
   ! fill in ghost cells for rho, rho Y, rho h, T according to the boundary
   ! conditions. fill in pthermo using the fourth order extrapolation
   subroutine scal_cc_to_avg(scal_avg, scal_cc)
      implicit none
      double precision, intent(out) :: scal_avg(-2:nx+1, nscal)
      double precision, intent(in ) ::  scal_cc(-2:nx+1, nscal)
      
      integer :: n, is
      
      do n=1,Nspec
         is = FirstSpec+n-1
         call cc_to_avg(scal_avg(:,is), scal_cc(:,is), rho_bc(on_lo)*Y_bc(n, on_lo))
      end do
      
      call cc_to_avg(scal_avg(:,Density), scal_cc(:,Density), rho_bc(on_lo))
      call cc_to_avg(scal_avg(:,Temp), scal_cc(:,Temp), T_bc(on_lo))
      call cc_to_avg(scal_avg(:,RhoH), scal_cc(:,RhoH), rho_bc(on_lo)*h_bc(on_lo))
      
      call extrapolate_cc_to_avg(scal_avg(:,RhoRT), scal_cc(:,RhoRT))
   end subroutine scal_cc_to_avg
   
   ! convert the array of scalar quantities from cell-averaged to cell-centered
   ! fill in ghost cells for rho, rho Y, rho h, T according to the boundary
   ! conditions. fill in pthermo using the fourth order extrapolation
   subroutine scal_avg_to_cc(scal_cc, scal_avg)
      implicit none
      double precision, intent(out) ::  scal_cc(-2:nx+1, nscal)
      double precision, intent(in ) :: scal_avg(-2:nx+1, nscal)
      
      integer :: n, is
      
      do n=1,Nspec
         is = FirstSpec+n-1
         call avg_to_cc(scal_cc(:,is), scal_avg(:,is), rho_bc(on_lo)*Y_bc(n, on_lo))
      end do
      
      call avg_to_cc(scal_cc(:,Density), scal_avg(:,Density), rho_bc(on_lo))
      call avg_to_cc(scal_cc(:,Temp), scal_avg(:,Temp), T_bc(on_lo))
      call avg_to_cc(scal_cc(:,RhoH), scal_avg(:,RhoH), rho_bc(on_lo)*h_bc(on_lo))
      call extrapolate_avg_to_cc(scal_cc(:,RhoRT), scal_avg(:,RhoRT))
      
   end subroutine scal_avg_to_cc
   
   ! convert a cell-centered quantity to a face-value quantity
   subroutine cc_to_face(face, cc)
      implicit none
      double precision, intent(out) :: face(0:nx)
      double precision, intent(in ) ::   cc(-2:nx+1)
      
      integer :: i
      
      do i=0,nx
         face(i) = 0.5d0*(cc(i) + cc(i-1))
         !face(i) = (-cc(i-2) + 9*cc(i-1) + 9*cc(i) - cc(i+1))/16.0
      end do
   end subroutine cc_to_face
   
   ! compute the cell-average of a product of two cell-averaged quantities
   ! according to the product rule, and fill the ghost cells according to 
   ! the supplied boundary condition
   subroutine mult_avgs(prod, avg_1, avg_2, bdry)
      implicit none
      double precision, intent(out) ::  prod(-2:nx+1)
      double precision, intent(in ) :: avg_1(-2:nx+1)
      double precision, intent(in ) :: avg_2(-2:nx+1)
      double precision, intent(in ) ::  bdry
      
      integer :: i
      
      do i=0,nx-1
         prod(i) = avg_1(i)*avg_2(i) &
            + (5*avg_1(i-2) - 34*avg_1(i-1) + 34*avg_1(i+1) - 5*avg_1(i+2)) &
             *(5*avg_2(i-2) - 34*avg_2(i-1) + 34*avg_2(i+1) - 5*avg_2(i+2))/27648.0
      end do
      
      call fill_avg_ghost_cells(prod, bdry)
   end subroutine mult_avgs
   
   ! compute the gradient given a cell-centered quantity
   subroutine cc_to_grad(grad, cc, dx)
      double precision, intent(out) :: grad( 0:nx)
      double precision, intent(in ) ::   cc(-2:nx+1)
      double precision, intent(in ) :: dx
      
      integer :: i
      
      do i=0,nx
         grad(i) = (cc(i-2) - 27*cc(i-1) + 27*cc(i) - cc(i+1))/(24*dx)
      end do
   end subroutine cc_to_grad
   
end module cell_conversions_module
