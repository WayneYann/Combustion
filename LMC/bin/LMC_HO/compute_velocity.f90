      ! integrate the constraint S to compute the velocity at faces
      subroutine compute_velocity(vel, S_avg, dx)
         implicit none
         include 'spec.h'
         double precision, intent(out  ) ::   vel(0:nx  )
         double precision, intent(in   ) :: S_avg(0:nx-1)
         double precision, intent(in   ) ::    dx
         
         integer :: i
         double precision :: x
         
         ! set the velocity at the inflow condition
         vel(0) = u_bc(0)
         ! integrate S
         do i=1,nx
            vel(i) = vel(i-1) + dx*S_avg(i-1)
         end do
         
!         do i=0,nx
!            x = dble(i)*dx
!            vel(i) = 15.d0 + 25.1d0*(1+tanh((x-0.615d0)/.039d0))
!         end do
      end
