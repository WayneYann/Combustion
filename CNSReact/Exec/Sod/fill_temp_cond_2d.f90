      subroutine ca_fill_temp_cond(lo,hi,state,state_l1,state_l2, &
                                               state_h1,state_h2, &
                                         tau,  tau_l1, tau_l2,     &
                                               tau_h1, tau_h2,     &
                                         coefx,coefx_l1,coefx_l2, &
                                               coefx_h1,coefx_h2, &
                                         coefy,coefy_l1,coefy_l2, &
                                               coefy_h1,coefy_h2, &
                                         dx)
      use network
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UTEMP, UFS, UFX
      use probdata_module   
      use interpolate_module

      implicit none
      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: state_l1,state_l2,state_h1,state_h2
      integer         , intent(in   ) :: tau_l1,tau_l2,tau_h1,tau_h2
      integer         , intent(in   ) :: coefx_l1,coefx_l2,coefx_h1,coefx_h2
      integer         , intent(in   ) :: coefy_l1,coefy_l2,coefy_h1,coefy_h2
      double precision, intent(in   ) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
      double precision, intent(in   ) :: tau( tau_l1: tau_h1, tau_l2: tau_h2)
      double precision, intent(inout) :: coefx(coefx_l1:coefx_h1,coefx_l2:coefx_h2)
      double precision, intent(inout) :: coefy(coefy_l1:coefy_h1,coefy_l2:coefy_h2)
      double precision, intent(in   ) :: dx(2)

      ! Local variables
      integer          :: i,j,n
      integer,save          :: ic12,iye
      double precision :: opac
      double precision :: coef_cc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)
      double precision :: temp,den,u,v,c,eden,xn(nspec)
!      double precision :: dummy_gam,dummy_pres,dummy_dpdr,dummy_dpde
      double precision :: flame_thickness
      double precision :: x12,rho_fuel,factor_w,dummy

      double precision, parameter :: aconst = 7.5646d-15, lightspeed = 3.0d10

      logical, save :: firstCall = .true.

      if (firstCall) then
         if (.NOT. network_initialized) then
            call bl_error("ERROR in EOS: must initialize network first")
         endif

         ic12= network_species_index("carbon-12")
         iye = network_species_index("Ye")

         firstCall = .false.
      endif

      do j = lo(2)-1,hi(2)+1
      do i = lo(1)-1,hi(1)+1

         temp = state(i,j,UTEMP)
         den  = state(i,j,URHO)
         x12  = state(i,j,UFS+ic12-1) / den

         u    = state(i,j,UMX) / den
         v    = state(i,j,UMY) / den
         eden = state(i,j,UEDEN)/ den - 0.5d0 * (u**2+v**2)

!         do n=1,nspec
!            xn(n) = state(i,j,UFS+n-1) / den
!         enddo

!         call eos_given_ReX(dummy_gam, dummy_pres , c, temp, &
!                           dummy_dpdr, dummy_dpde, &
!                           den, eden, xn)

     ! define rho_fuel from burning table
!        call networkburn(log10(den),x12,rho_fuel,dummy,dummy,&
!                        dummy,dummy,dummy,dummy,&
!                        dummy,dummy,dummy,dummy) 

         coef_cc(i,j) = 0.0d0
      enddo
      enddo

      do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
          coefx(i,j) = 0.5d0 * (coef_cc(i,j) + coef_cc(i-1,j))
        end do
      end do

      do j = lo(2),hi(2)+1
        do i = lo(1),hi(1)
          coefy(i,j) = 0.5d0 * (coef_cc(i,j) + coef_cc(i,j-1))
        end do
      end do

      end subroutine ca_fill_temp_cond
