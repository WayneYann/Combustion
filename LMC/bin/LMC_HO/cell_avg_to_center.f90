      subroutine cell_avg_to_center(cc, avg, lo, hi)
      implicit none
      double precision, intent(in ) ::  cc(*)
      double precision, intent(out) :: avg(*)
      integer,          intent(in ) :: lo, hi
      
      integer :: i
      
      do i=lo,hi
         cc(i) = avg(i) - (avg(i-1) - 2.0*avg(i) + avg(i+1))/24.0
      end do
      
      end subroutine cell_avg_to_center
