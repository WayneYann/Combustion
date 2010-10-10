      subroutine ca_fill_temp_cond(lo,hi,state,state_l1,state_l2,state_l3, &
                                               state_h1,state_h2,state_h3, &
                                         tau,  tau_l1, tau_l2, tau_l3,     &
                                               tau_h1, tau_h2, tau_h3,     &
                                         coefx,coefx_l1,coefx_l2,coefx_l3, &
                                               coefx_h1,coefx_h2,coefx_h3, &
                                         coefy,coefy_l1,coefy_l2,coefy_l3, &
                                               coefy_h1,coefy_h2,coefy_h3, &
                                         coefz,coefz_l1,coefz_l2,coefz_l3, &
                                               coefz_h1,coefz_h2,coefz_h3, &
                                         dx)
      use network
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS, UFX
      use probdata_module   
      use interpolate_module

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
      integer         , intent(in   ) :: tau_l1,tau_l2,tau_l3,tau_h1,tau_h2,tau_h3
      integer         , intent(in   ) :: coefx_l1,coefx_l2,coefx_l3,coefx_h1,coefx_h2,coefx_h3
      integer         , intent(in   ) :: coefy_l1,coefy_l2,coefy_l3,coefy_h1,coefy_h2,coefy_h3
      integer         , intent(in   ) :: coefz_l1,coefz_l2,coefz_l3,coefz_h1,coefz_h2,coefz_h3
      double precision, intent(in   ) :: state(state_l1:state_h1,state_l2:state_h2,&
                                               state_l3:state_h3,NVAR)
      double precision, intent(in   ) :: tau( tau_l1: tau_h1, tau_l2: tau_h2, &
                                               tau_l3:tau_h3)  
      double precision, intent(inout) :: coefx(coefx_l1:coefx_h1,coefx_l2:coefx_h2, &
                                               coefx_l3:coefx_h3)
      double precision, intent(inout) :: coefy(coefy_l1:coefy_h1,coefy_l2:coefy_h2, &
                                               coefy_l3:coefy_h3)
      double precision, intent(inout) :: coefz(coefz_l1:coefz_h1,coefz_l2:coefz_h2, &
                                               coefz_l3:coefz_h3)
      double precision, intent(in   ) :: dx(3)

      ! Local variables
      integer          :: i,j,k,n
      integer,save          :: ic12,iye
      double precision :: opac
      double precision :: coef_cc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
      double precision :: temp,den,u,v,w,c,eden,xn(nspec)
!      double precision :: dummy_gam,dummy_pres,dummy_dpdr,dummy_dpde
      double precision :: flame_thickness
      double precision :: x12,rho_fuel,factor_w,dummy

      double precision, parameter :: aconst = 7.5646d-15, lightspeed = 3.0d10

      logical, save :: firstCall = .true.

      !$OMP CRITICAL (fill_temp_lock)
      if (firstCall) then
         if (.NOT. network_initialized) then
            call bl_error("ERROR in EOS: must initialize network first")
         endif

         ic12= network_species_index("carbon-12")
         iye = network_species_index("Ye")

         firstCall = .false.
      endif
      !$OMP END CRITICAL (fill_temp_lock)

      !$OMP PARALLEL DO PRIVATE(i,j,k,temp,den,x12,u,v,w,eden,rho_fuel,dummy,factor_w,flame_thickness,opac)
      do k = lo(3)-1,hi(3)+1
      do j = lo(2)-1,hi(2)+1
      do i = lo(1)-1,hi(1)+1

         temp = state(i,j,k,UTEMP)
         den  = state(i,j,k,URHO)
         x12  = state(i,j,k,UFS+ic12-1) / den

         u    = state(i,j,k,UMX) / den
         v    = state(i,j,k,UMY) / den
         w    = state(i,j,k,UMZ) / den
         eden = state(i,j,k,UEDEN)/ den - 0.5d0 * (u**2+v**2+w**2)


!         do n=1,nspec
!            xn(n) = state(i,j,k,UFS+n-1) / den
!         enddo

!         call eos_given_ReX(dummy_gam, dummy_pres , c, temp, &
!                           dummy_dpdr, dummy_dpde, &
!                           den, eden, xn)

     ! define rho_fuel from burning table
!        call networkburn(log10(den),x12,rho_fuel,dummy,dummy,&
!                        dummy,dummy,dummy,dummy,&
!                        dummy,dummy,dummy,dummy) 

         coef_cc(i,j,k) = 0.0d0
      end do
      enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL PRIVATE(i,j,k)
      !$OMP DO
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
          coefx(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i-1,j,k))
        end do
      end do
      end do
      !$OMP END DO NOWAIT

      !$OMP DO
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)+1
        do i = lo(1),hi(1)
          coefy(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j-1,k))
        end do
      end do
      end do
      !$OMP END DO NOWAIT

      !$OMP DO
      do k = lo(3),hi(3)+1
      do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          coefz(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j,k-1))
        end do
      end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      end subroutine ca_fill_temp_cond
