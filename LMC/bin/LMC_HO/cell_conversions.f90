module cell_conversions_module
   
   use ghost_cells_module
   
   implicit none
   private
   
   include 'spec.h'
   
   public :: cc_to_avg, extrapolate_cc_to_avg, avg_to_cc

contains

   subroutine cc_to_avg(avg, cc)
      implicit none
      double precision, intent(in ) ::  cc(-2:nx+1)
      double precision, intent(out) :: avg(-2:nx+1)

      integer :: i

      ! note: we do not fill in the ghost cells here
      do i=0,nx-1
         avg(i) = cc(i) + (cc(i-1) - 2.0*cc(i) + cc(i+1))/24.0
      end do
      
   end subroutine cc_to_avg
         
   subroutine extrapolate_cc_to_avg(avg, cc)
      implicit none
      double precision, intent(in ) ::  cc(0:nx-1)
      double precision, intent(out) :: avg(0:nx-1)
      
      integer :: i
      
      ! one-sided stencils
      avg(0) = (27.0*cc(0) - 9.0*cc(1) + 10.0*cc(2) - 5.0*cc(3) + cc(4))/24.0
      avg(nx-1) = (cc(nx-5) - 5.0*cc(nx-4) + 10.0*cc(nx-3) &
                 - 9.0*cc(nx-2) + 27.0*cc(nx-2))/24.0
      
      do i=1,nx-2
         avg(i) = cc(i) + (cc(i-1) - 2.0*cc(i) + cc(i+1))/24.0
      end do
   end subroutine extrapolate_cc_to_avg
   
   subroutine avg_to_cc(cc, avg)
      implicit none
      double precision, intent(out) ::  cc(-2:nx+1)
      double precision, intent(in ) :: avg(-2:nx+1)
      
      integer :: i
      
      ! note: we do not fill in the ghost cells here
      do i=0,nx-1
         cc(i) = avg(i) - (avg(i-1) - 2.0*avg(i) + avg(i+1))/24.0
      end do
      
   end subroutine avg_to_cc

end module cell_conversions_module
