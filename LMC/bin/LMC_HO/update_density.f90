      subroutine update_density(scal_mp1_avg, scal_m_avg, advection_kp1, advection_k, I_k, dtm)
      implicit none
      include 'spec.h'
      double precision, intent(inout) ::  scal_mp1_avg(-2:nx+1,nscal)
      double precision, intent(in   ) ::    scal_m_avg(-2:nx+1,nscal)
      double precision, intent(in   ) :: advection_kp1( 0:nx-1,nscal)
      double precision, intent(in   ) ::   advection_k( 0:nx-1,nscal)
      double precision, intent(in   ) ::           I_k( 0:nx-1,nscal)
      double precision, intent(in   ) ::           dtm
      
      integer :: i
      
      do i=0,nx-1
         scal_mp1_avg(i, Density) = scal_m_avg(i, Density) &
            + dtm*(advection_kp1(i, Density) - advection_k(i, Density)) &
            + I_k(i, Density)
      end do
      
      end subroutine update_density
