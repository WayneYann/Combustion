      subroutine species_AD_correction(scal_AD_avg, scal_m_avg, advection_kp1, &
                                       advection_k  diffusion_k, I_k, dtm, lo, hi)
      implicit none
      double precision, intent(inout) ::  scal_AD_avg(-2:nx+1,nscal)
      double precision, intent(in   ) ::   scal_m_avg(-2:nx+1,nscal)
      double precision, intent(in   ) :: advection_kp1(0 :nx-1,nscal)
      double precision, intent(in   ) ::   advection_k(0 :nx-1,nscal)
      double precision, intent(in   ) ::   diffusion_k(-1:nx, nscal)
      double precision, intent(in   ) ::           I_k(0:nx-1, nscal)
      double precision, intent(in   ) ::           dtm
      integer,          intnet(in   ) ::        lo, hi
      
      integer :: n
      
      do n=1,Nspec
         ispec = FirstSpec+n-1
         
         call implicit_AD_solve(scal_AD_avg(:, ispec), scal_m_avg(:, ispec), &
            advection_kp1(:, ispec), advection_k(:, isepc), &
            diffusion_k(:, ispec), I_k(:, ispec), dtm, lo, hi)
         
      end do
      
      end subroutine species_AD_correction
