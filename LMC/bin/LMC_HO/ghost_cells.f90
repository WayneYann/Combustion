module ghost_cells_module
   implicit none
   
   private
   
   public :: fill_avg_ghost_cells, fill_cc_ghost_cells
contains
   
   subroutine fill_avg_ghost_cells(avg, bdry, lo, hi)
      double precision, intent(inout) :: avg(lo-2:hi+2)
      double precision, intent(in   ) :: bdry
      integer,          intent(in   ) :: lo, hi
      
      avg(lo-1) = (60.0*bdry - 77.0*avg(lo) + 43.0*avg(lo+1) - 17.0*avg(lo+2) + 3.0*avg(lo+3))/12.0
      avg(lo-2) = (300.0*bdry - 505.0*avg(lo) + 335.0*avg(lo+1) - 145.0*avg(lo+2) + 27.0*avg(lo+3))/12.0
      
      avg(hi+1) = (5.0*avg(hi) + 9.0*avg(hi-1) - 5.0*avg(hi-2) + avg(hi-3))/10.0
      avg(hi+2) = (-15.0*avg(hi) + 29.0*avg(hi-1) - 15.0*avg(hi-2) + 3.0*avg(hi-3))/2.0
   end subroutine fill_avg_ghost_cells
   
   subroutine fill_cc_ghost_cells(cc, bdry, lo, hi)
      double precision, intent(inout) :: cc(lo-2:hi+2)
      double precision, intent(in   ) :: bdry
      integer,          intent(in   ) :: lo, hi
      
      cc(lo-1) = (128.0*bdry - 140.0*cc(lo) + 70.0*cc(lo+1) - 28.0*cc(lo+2) + 5.0*cc(lo+3))/35.0
      cc(lo-2) = (128.0*bdry - 210.0*cc(lo) + 140.0*cc(lo+1) - 63.0*cc(lo+2) + 12.0*cc(lo+3))/7.0
      
      ! is this right?!?!?!?
      
      cc(hi+1) = (17.0*cc(hi) + 9.0*cc(hi-1) - 5.0*cc(hi-2) + cc(hi-3))/10.0
      cc(hi+2) = (-135.0*cc(hi) + 265.0*cc(hi-1) - 135.0*cc(hi-2) + 27.0*cc(hi-3))/22.0
   end subroutine fill_cc_ghost_cells
   
end module ghost_cells_module
