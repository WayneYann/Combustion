      subroutine implicit_AD_correction(phi_AD_avg, phi_m_avg, advection_kp1, &
                                       advection_k  diffusion_k, I_k, dtm, lo, hi)
      implicit none
      double precision, intent(out  ) ::   phi_AD_avg(-2:nfine+1)
      double precision, intent(in   ) ::    phi_m_avg(-2:nfine+1)
      double precision, intent(in   ) :: advection_kp1(0 :nfine-1)
      double precision, intent(in   ) ::   advection_k(0 :nfine-1)
      double precision, intent(in   ) ::   diffusion_k(-1:nfine)
      double precision, intent(in   ) ::           I_k(0:nfine-1)
      double precision, intent(in   ) ::           dtm
      integer,          intnet(in   ) ::        lo, hi
      
      double precision :: rhs(nfine)
      double precision :: phi(nfine)
      integer :: i
      
      
      ! construct the right hand side
      do i=lo,hi
         rhs(i) = scal_m_avg(i) &
            + dtm*(advection_kp1(i) - advection_k(i) - diffusion_k(i)) + I_k(i)
      end do
         
      ! construct the matrix we need to invert
      
      ! perform the sparse linear solve for phi
      
      ! fill in scal_AD_avg with the values from phi
      
      ! fill in the ghost cells
      
      end subroutine implicit_AD_correction
