      subroutine increment_delta_chi(delta_chi, scal_cc, dt)
      implicit none
      include 'spec.h'
      double precision, intent(in   ) ::   scal_cc(-2:nx+1,nscal)
      double precision, intent(inout) :: delta_chi( 0:nx-1)
      double precision, intent(in   ) :: dt
      
      double precision :: dpdt
      integer :: i
      
      call compute_pthermo(scal_cc)
      
      do i=0,nx-1
         dpdt = (scal_cc(i, RhoRT) - Pcgs) / dt
         delta_chi(i) = delta_chi(i) + dpdt*dpdt_factor/Pcgs
      end do
      
      end subroutine increment_delta_chi
