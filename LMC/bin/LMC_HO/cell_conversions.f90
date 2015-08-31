module cell_conversions_module
   
   use ghost_cells_module
   
   implicit none
   private
   
   include 'spec.h'
   
   public :: cc_to_avg, extrapolate_cc_to_avg, avg_to_cc, scal_cc_to_avg
   public :: cc_to_face, mult_avgs_bdry, scal_avg_to_cc, extrapolate_avg_to_cc
   public :: cc_to_grad, mult_avgs, scal_cc_to_face, avg_to_grad
   public :: avg_to_face, scal_avg_to_face, divide_avgs

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
         avg(i) = cc(i) + (cc(i-1) - 2.d0*cc(i) + cc(i+1))/24.d0
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
      avg(0) = (27*cc(0) - 9*cc(1) + 10*cc(2) - 5*cc(3) + cc(4))/24.d0
      avg(nx-1) = (501*cc(nx-1) + 31*cc(nx-2) - 5*cc(nx-3) + cc(nx-4))/528.d0
      
      !avg(nx-1) = (cc(nx-5) - 5*cc(nx-4) + 10*cc(nx-3) &
      !           - 9*cc(nx-2) + 27*cc(nx-2))/24.d0
      
      do i=1,nx-2
         avg(i) = cc(i) + (cc(i-1) - 2*cc(i) + cc(i+1))/24.d0
      end do
   end subroutine extrapolate_cc_to_avg
   
   subroutine extrapolate_avg_to_cc(cc, avg)
      implicit none
      double precision, intent(out) ::  cc(0:nx-1)
      double precision, intent(in ) :: avg(0:nx-1)
      
      integer :: i
      
      cc(0) = (21*avg(0) + 9*avg(1) - 10*avg(2) + 5*avg(3) - avg(4))/24.d0
      ! use the outflow condition
      cc(nx-1) = (255*avg(nx-1) - 19*avg(nx-2) + 5*avg(nx-3) - avg(nx-4))/240.d0
      
      do i=1,nx-2
         cc(i) = avg(i) - (avg(i-1) - 2*avg(i) + avg(i+1))/24.d0
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
         cc(i) = avg(i) - (avg(i-1) - 2*avg(i) + avg(i+1))/24.d0
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
      
      call extrapolate_cc_to_avg(scal_avg(0:,RhoRT), scal_cc(0:,RhoRT))
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
      
      call extrapolate_avg_to_cc(scal_cc(0:,RhoRT), scal_avg(0:,RhoRT))
   end subroutine scal_avg_to_cc
   
   ! convert a cell-centered quantity to a face-value quantity
   subroutine cc_to_face(face, cc)
      implicit none
      double precision, intent(out) :: face(0:nx)
      double precision, intent(in ) ::   cc(-2:nx+1)
      
      integer :: i
      
      do i=0,nx
         face(i) = (-cc(i-2) + 9*cc(i-1) + 9*cc(i) - cc(i+1))/16.d0
      end do
   end subroutine cc_to_face
   
   ! convert a cell average to a face-value quantity
   subroutine avg_to_face(face, avg)
      implicit none
      double precision, intent(out) :: face(0:nx)
      double precision, intent(in ) ::  avg(-2:nx+1)
      
      integer :: i
      
      do i=0,nx
         face(i) = (-avg(i-2) + 7*avg(i-1) + 7*avg(i) - avg(i+1))/12.d0
      end do
   end subroutine avg_to_face
   
   subroutine scal_cc_to_face(scal_face, scal_cc)
      implicit none
      double precision, intent(out) :: scal_face( 0:nx,   nscal)
      double precision, intent(in ) ::   scal_cc(-2:nx+1, nscal)
      
      integer :: n, is
      
      do n=1,Nspec
         is = FirstSpec+n-1
         call cc_to_face(scal_face(:,is), scal_cc(:,is))
         scal_face(0,is) = rho_bc(on_lo)*Y_bc(n,on_lo)
      end do
      
      call cc_to_face(scal_face(:,Density), scal_cc(:,Density))
      scal_face(0,Density) = rho_bc(on_lo)
      call cc_to_face(scal_face(:,Temp), scal_cc(:,Temp))
      scal_face(0,Temp) = T_bc(on_lo)
      call cc_to_face(scal_face(:,RhoH), scal_cc(:,RhoH))
      scal_face(0,RhoH) = rho_bc(on_lo)*h_bc(on_lo)
      call cc_to_face(scal_face(:,RhoRT), scal_cc(:,RhoRT))
   end subroutine scal_cc_to_face
   
   subroutine scal_avg_to_face(scal_face, scal_avg)
      implicit none
      double precision, intent(out) :: scal_face( 0:nx,   nscal)
      double precision, intent(in ) ::  scal_avg(-2:nx+1, nscal)
      
      integer :: n, is
      
      do n=1,Nspec
         is = FirstSpec+n-1
         call avg_to_face(scal_face(:,is), scal_avg(:,is))
         scal_face(0,is) = rho_bc(on_lo)*Y_bc(n,on_lo)
      end do
      
      call avg_to_face(scal_face(:,Density), scal_avg(:,Density))
      scal_face(0,Density) = rho_bc(on_lo)
      call avg_to_face(scal_face(:,Temp), scal_avg(:,Temp))
      scal_face(0,Temp) = T_bc(on_lo)
      call avg_to_face(scal_face(:,RhoH), scal_avg(:,RhoH))
      scal_face(0,RhoH) = rho_bc(on_lo)*h_bc(on_lo)
      call avg_to_face(scal_face(:,RhoRT), scal_avg(:,RhoRT))
   end subroutine scal_avg_to_face
   
   ! compute the cell-average of a product of two cell-averaged quantities
   ! according to the product rule, and fill the ghost cells according to 
   ! the supplied boundary condition
   subroutine mult_avgs_bdry(prod, avg_1, avg_2, bdry)
      implicit none
      double precision, intent(out) ::  prod(-2:nx+1)
      double precision, intent(in ) :: avg_1(-2:nx+1)
      double precision, intent(in ) :: avg_2(-2:nx+1)
      double precision, intent(in ) ::  bdry
      
      integer :: i
      
      do i=0,nx-1
         prod(i) = avg_1(i)*avg_2(i) &
            + (5*avg_1(i-2) - 34*avg_1(i-1) + 34*avg_1(i+1) - 5*avg_1(i+2)) &
             *(5*avg_2(i-2) - 34*avg_2(i-1) + 34*avg_2(i+1) - 5*avg_2(i+2))/27648.d0
      end do
      
      call fill_avg_ghost_cells(prod, bdry)
   end subroutine mult_avgs_bdry
   
   ! compute the cell-average of a product of two cell-averaged quantities
   ! according to the product rule, only in the valid region (no ghost cells!)
   ! the supplied boundary condition
   subroutine mult_avgs(prod, avg_1, avg_2)
      implicit none
      double precision, intent(out) ::  prod( 0:nx-1)
      double precision, intent(in ) :: avg_1(-2:nx+1)
      double precision, intent(in ) :: avg_2(-2:nx+1)
      
      integer :: i
      
      do i=0,nx-1
         prod(i) = avg_1(i)*avg_2(i) &
            + (5*avg_1(i-2) - 34*avg_1(i-1) + 34*avg_1(i+1) - 5*avg_1(i+2)) &
             *(5*avg_2(i-2) - 34*avg_2(i-1) + 34*avg_2(i+1) - 5*avg_2(i+2))/27648.d0
      end do
   end subroutine mult_avgs
   
   
   subroutine divide_avgs(q, avg_1, avg_2, bdry)
      implicit none
      double precision, intent(out) ::     q(-2:nx+1)
      double precision, intent(in ) :: avg_1(-2:nx+1)
      double precision, intent(in ) :: avg_2(-2:nx+1)
      double precision, intent(in ) ::  bdry
      
      double precision :: p_1, r_1
      
      integer :: i
      
      do i=0,nx-1
         p_1 = 5*avg_1(i-2) - 34*avg_1(i-1) + 34*avg_1(i+1) - 5*avg_1(i+2)
         r_1 = 5*avg_2(i-2) - 34*avg_2(i-1) + 34*avg_2(i+1) - 5*avg_2(i+2)
         
         q(i) = avg_1(i)/avg_2(i) &
            + (avg_1(i)*r_1*r_1/(avg_2(i)**3) - p_1*r_1/(avg_2(i)**2))/27648.d0
      end do
      
      call fill_avg_ghost_cells(q, bdry)
   end subroutine divide_avgs
   
   ! compute the gradient given a cell-centered quantity
   subroutine cc_to_grad(grad, cc, dx)
      double precision, intent(out) :: grad( 0:nx)
      double precision, intent(in ) ::   cc(-2:nx+1)
      double precision, intent(in ) :: dx
      
      integer :: i
      
      do i=0,nx
         grad(i) = (cc(i-2) - 27*cc(i-1) + 27*cc(i) - cc(i+1))/(24.d0*dx)
      end do
   end subroutine cc_to_grad
   
   ! compute the gradient given a cell-averaged quantity
   subroutine avg_to_grad(grad, avg, dx)
      double precision, intent(out) :: grad( 0:nx)
      double precision, intent(in ) ::  avg(-2:nx+1)
      double precision, intent(in ) ::   dx
      
      integer :: i
      
      do i=0,nx
         grad(i) = (avg(i-2) - 15*avg(i-1) + 15*avg(i) - avg(i+1))/(12.d0*dx)
      end do
   end subroutine avg_to_grad
   
end module cell_conversions_module
