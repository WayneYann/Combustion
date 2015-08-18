      ! integrate the constraint S to compute the velocity at faces
      subroutine compute_velocity(vel, S, dx)
         implicit none
         include 'spec.h'
         double precision, intent(out  ) :: vel(0:nx  )
         double precision, intent(in   ) ::   S(0:nx-1)
         double precision, intent(in   ) ::  dx
         
         integer :: i
         
         ! set the velocity at the inflow condition
         vel(0) = u_bc(0)
         ! integrate S
         do i=1,nx
            vel(i) = vel(i-1) + dx*S(i-1)
         end do
      end
