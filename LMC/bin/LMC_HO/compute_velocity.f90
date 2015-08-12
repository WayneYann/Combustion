      subroutine compute_velocity(vel, S, dx, lo, hi)
      implicit none
      include 'spec.h'
      double precision, intent(out  ) :: vel(0 :nx)
      double precision, intent(in   ) ::   S(-1:nx)
      double precision, intent(in   ) ::  dx
      integer,          intent(in   ) ::  lo, hi
      
      integer :: i
      
      ! need to actually implement the inflow velocity
      vel(lo) = vel_inflow
      do i=lo+1,hi
         vel(i) = dx*S(i-1) + vel(i-1)
      end do
      
      end
