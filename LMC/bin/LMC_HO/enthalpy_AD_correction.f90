      subroutine enthalpy_AD_correction(scal_AD_avg, scal_m_avg, advection_kp1, &
                                       advection_k  diffusion_k, I_k, dtm, lo, hi)
      implicit none
      double precision, intent(inout) ::   scal_AD_avg(-2:nx+1,nscal)
      double precision, intent(in   ) ::    scal_m_avg(-2:nx+1,nscal)
      double precision, intent(in   ) :: advection_kp1( 0:nx-1,nscal)
      double precision, intent(in   ) ::   advection_k( 0:nx-1,nscal)
      double precision, intent(in   ) ::   diffusion_k(-1:nx, nscal)
      double precision, intent(in   ) ::           I_k( 0:nx-1, nscal)
      double precision, intent(in   ) ::           dtm
      integer,          intnet(in   ) ::        lo, hi
      
      call implicit_AD_solve(scal_AD_avg(:, RhoH), scal_m_avg(:, RhoH), &
         advection_kp1(:, RhoH), advection_k(:, isepc), &
         diffusion_k(:, RhoH), I_k(:, RhoH), dtm, lo, hi)
      
      end subroutine enthalpy_AD_correction
