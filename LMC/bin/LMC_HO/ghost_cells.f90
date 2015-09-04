module ghost_cells_module
   implicit none
   
   private
   include 'spec.h'
   
   public :: fill_avg_ghost_cells, fill_cc_ghost_cells, fill_scal_cc_ghost_cells
   public :: fill_scal_avg_ghost_cells, extrapolate_avg_ghost_cells
contains
   
   ! fill cell-averaged ghost cells (Dirichlet inflow, Neumann outflow)
   subroutine fill_avg_ghost_cells(avg, bdry)
      double precision, intent(inout) :: avg(-2:nx+1)
      double precision, intent(in   ) :: bdry
      
      avg(-1) = (60*bdry - 77*avg(0) + 43*avg(1) - 17*avg(2) + 3*avg(3))/12
      avg(-2) = (300*bdry - 505*avg(0) + 335*avg(1) - 145*avg(2) + 27*avg(3))/12
      
      avg(nx) = (5*avg(nx-1) + 9*avg(nx-2) - 5*avg(nx-3) + avg(nx-4))/10.d0
      avg(nx+1) = (-15*avg(nx-1) + 29*avg(nx-2) - 15*avg(nx-3) + 3*avg(nx-4))/2.d0
   end subroutine fill_avg_ghost_cells
   
   ! fill cell-averaged ghost cells (extrapolation at inflow, Neumann outflow)
   subroutine extrapolate_avg_ghost_cells(avg)
      double precision, intent(inout) :: avg(-2:nx+1)
      
      avg(-1) = 5*avg(0) - 10*avg(1) + 10*avg(2) - 5*avg(3) + avg(4)
      avg(-2) = 15*avg(0) - 40*avg(1) + 45*avg(2) - 24*avg(3) + 5*avg(4)
      
      avg(nx) = (5*avg(nx-1) + 9*avg(nx-2) - 5*avg(nx-3) + avg(nx-4))/10.d0
      avg(nx+1) = (-15*avg(nx-1) + 29*avg(nx-2) - 15*avg(nx-3) + 3*avg(nx-4))/2.d0
   end subroutine extrapolate_avg_ghost_cells
   
   ! fill cell-centered ghost cells (Dirichlet inflow, Neumann outflow)
   subroutine fill_cc_ghost_cells(cc, bdry)
      double precision, intent(inout) :: cc(-2:nx+1)
      double precision, intent(in   ) :: bdry
      
      cc(-1) = (128*bdry - 140*cc(0) + 70*cc(1) - 28*cc(2) + 5*cc(3))/35.d0
      cc(-2) = (128*bdry - 210*cc(0) + 140*cc(1) - 63*cc(2) + 12*cc(3))/7.d0
      
      cc(nx) = (17*cc(nx-1) + 9*cc(nx-2) - 5*cc(nx-3) + cc(nx-4))/22.d0
      cc(nx+1) = (-135*cc(nx-1) + 265*cc(nx-2) - 135*cc(nx-3) + 27*cc(nx-4))/22.d0
   end subroutine fill_cc_ghost_cells
   
   ! fill the ghost cells for all the scalar quantities (cell-centered)
   subroutine fill_scal_cc_ghost_cells(scal_cc)
      double precision, intent(inout) :: scal_cc(-2:nx+1, nscal)
      integer :: n
      
      call fill_cc_ghost_cells(scal_cc(:,Density), rho_bc(on_lo))
      call fill_cc_ghost_cells(scal_cc(:,Temp), T_bc(on_lo))
      call fill_cc_ghost_cells(scal_cc(:,RhoH), h_bc(on_lo)*rho_bc(on_lo))
      
      do n=1,Nspec
         call fill_cc_ghost_cells(scal_cc(:,FirstSpec+n-1), Y_bc(n,on_lo)*rho_bc(on_lo))
      end do
   end subroutine fill_scal_cc_ghost_cells
   
   ! fill the ghost cells for all the scalar quantities (cell-averaged)
   subroutine fill_scal_avg_ghost_cells(scal_avg)
      double precision, intent(inout) :: scal_avg(-2:nx+1, nscal)
      integer :: n
      
      call fill_avg_ghost_cells(scal_avg(:,Density), rho_bc(on_lo))
      call fill_avg_ghost_cells(scal_avg(:,Temp), T_bc(on_lo))
      call fill_avg_ghost_cells(scal_avg(:,RhoH), h_bc(on_lo)*rho_bc(on_lo))
      
      do n=1,Nspec
         call fill_avg_ghost_cells(scal_avg(:,FirstSpec+n-1), Y_bc(n,on_lo)*rho_bc(on_lo))
      end do
   end subroutine fill_scal_avg_ghost_cells
   
end module ghost_cells_module
