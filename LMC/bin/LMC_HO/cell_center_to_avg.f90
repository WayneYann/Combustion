      subroutine cell_center_to_avg(avg, cc, lo, hi)
      implicit none
      double precision, intent(in ) ::  cc(*)
      double precision, intent(out) :: avg(*)
      integer,          intent(in ) :: lo, hi
      
      integer :: i
      
      do i=lo,hi
         avg(i) = cc(i) + (cc(i-1) - 2.0*cc(i) + cc(i+1))/24.0
      end do
      
      end subroutine cell_center_to_avg
