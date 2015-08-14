module ghost_cells_module
   implicit none
   
   private
   include 'spec.h'
   
   public :: fill_avg_ghost_cells, fill_cc_ghost_cells, fill_scal_cc_ghost_cells
   public :: fill_scal_avg_ghost_cells
contains
   
   subroutine fill_avg_ghost_cells(avg, bdry)
      double precision, intent(inout) :: avg(-2:nx+1)
      double precision, intent(in   ) :: bdry
      
      avg(-1) = (60.0*bdry - 77.0*avg(0) + 43.0*avg(1) - 17.0*avg(2) + 3.0*avg(3))/12.0
      avg(-2) = (300.0*bdry - 505.0*avg(0) + 335.0*avg(1) - 145.0*avg(2) + 27.0*avg(3))/12.0
      
      avg(nx) = (5.0*avg(nx-1) + 9.0*avg(nx-2) - 5.0*avg(nx-3) + avg(nx-4))/10.0
      avg(nx+1) = (-15.0*avg(nx-1) + 29.0*avg(nx-2) - 15.0*avg(nx-3) + 3.0*avg(nx-4))/2.0
   end subroutine fill_avg_ghost_cells
   
   subroutine fill_cc_ghost_cells(cc, bdry)
      double precision, intent(inout) :: cc(-2:nx+1)
      double precision, intent(in   ) :: bdry
      
      cc(-1) = (128.0*bdry - 140.0*cc(0) + 70.0*cc(1) - 28.0*cc(2) + 5.0*cc(3))/35.0
      cc(-2) = (128.0*bdry - 210.0*cc(0) + 140.0*cc(1) - 63.0*cc(2) + 12.0*cc(3))/7.0
      
      cc(nx) = (17.0*cc(nx-1) + 9.0*cc(nx-2) - 5.0*cc(nx-3) + cc(nx-4))/22.0
      cc(nx+1) = (-135.0*cc(nx-1) + 265.0*cc(nx-2) - 135.0*cc(nx-3) + 27.0*cc(nx-4))/22.0
   end subroutine fill_cc_ghost_cells
   
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
