module diff_flux_module

  implicit none

  private

  public :: diffFlux, fluxtosrc, diffup
 
contains

  subroutine diffFlux(lof,hif, &
       q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
       flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
       flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
       flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
       dx,dy,dz,df_timer)

    use meth_params_module, only : NVAR, QVAR
    use chemistry_module, only : nspec=>nspecies
    
    implicit none
    
    integer,intent(in):: lof(3),hif(3), df_timer
    integer,intent(in)::    q_l1,    q_l2,    q_l3,    q_h1,    q_h2,    q_h3
    integer,intent(in)::flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
    integer,intent(in)::flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
    integer,intent(in)::flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
    double precision,intent(in) :: dx,dy,dz
    double precision,intent(in   )::    q(    q_l1:    q_h1,    q_l2:    q_h2,    q_l3:    q_h3,QVAR)
    double precision,intent(inout)::flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision,intent(inout)::flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision,intent(inout)::flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)

    integer :: loD(3), hiD(3)
    double precision, allocatable :: rhoYD(:,:,:,:), lamcp(:,:,:), mu(:,:,:)

    if (df_timer .eq. 1) then
       flux1(:,:,:,:) = 0.d0
       flux2(:,:,:,:) = 0.d0
       flux3(:,:,:,:) = 0.d0
    end if

    loD = lof - 1
    hiD = hif + 1
    allocate(rhoYD(loD(1):hiD(1),loD(2):hiD(2),loD(3):hiD(3),nspec))
    allocate(lamcp(loD(1):hiD(1),loD(2):hiD(2),loD(3):hiD(3)))
    allocate(mu   (loD(1):hiD(1),loD(2):hiD(2),loD(3):hiD(3)))

    call get_transCoef(q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
         rhoYD, lamcp, mu, loD, hiD)

    call comp_dflux(lof,hif, &
         q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
         rhoYD, lamcp, mu, loD, hiD, &
         flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
         flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
         flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
         dx,dy,dz)

    deallocate(rhoYD, lamcp, mu)

  end subroutine diffFlux


  subroutine get_transCoef(q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
         rhoYD, lamcp, mu, loD, hiD)

    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QREINT, QPRES, QTEMP, QFS
    use chemistry_module, only : nspec=>nspecies, inv_mwt
    use eglib_module

    implicit none

    integer, intent(in) :: q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
    integer, intent(in) :: loD(3), hiD(3)
    double precision,intent(in )::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
    double precision,intent(out)::rhoYD(loD(1):hiD(1),loD(2):hiD(2),loD(3):hiD(3),nspec)
    double precision,intent(out)::lamcp(loD(1):hiD(1),loD(2):hiD(2),loD(3):hiD(3))
    double precision,intent(out)::   mu(loD(1):hiD(1),loD(2):hiD(2),loD(3):hiD(3))

    integer :: i, j, k, n, iwrk
    double precision :: rwrk, rhoD(nspec), Cpt(nspec), Xt(nspec), Yt(nspec), Tt, Wtm, &
         lam1, lam2, cpmix

    integer, parameter :: NP=1, ITLS=1, IFLAG=3
    double precision, parameter :: TMIN_TRANS=300.d0

    call eglib_init(nspec, NP, ITLS, IFLAG)

    do       k = loD(3), hiD(3)
       do    j = loD(2), hiD(2)
          do i = loD(1), hiD(1)

             Tt = max(TMIN_TRANS, q(i,j,k,QTEMP))
             Yt = q(i,j,k,QFS:QFS+nspec-1)
             call ckytx(Yt,iwrk,rwrk,Xt)
             call ckcpms(Tt, iwrk, rwrk, Cpt)

             call EGSPAR(Tt, Xt, Yt, Cpt, egwork, egiwork)

             call EGSE3(Tt, Yt, egwork, mu(i,j,k))

             call EGSVR1(Tt, Yt, egwork, rhoD)
             call ckmmwy(Yt, iwrk, rwrk, Wtm)
             rhoYD(i,j,k,:) = Yt * rhoD * Wtm * inv_mwt

             call EGSL1( 1.d0, Tt, Xt, egwork, lam1) 
             call EGSL1(-1.d0, Tt, Xt, egwork, lam2) 
             call ckcpbs(Tt, Yt, iwrk, rwrk, cpmix)

             lamcp(i,j,k) = 0.5d0*(lam1+lam2) / cpmix
          end do
       end do
    end do

  end subroutine get_transCoef


  subroutine comp_dflux( lof, hif, &
       q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
       rhoYD, lamcp, mu, loD, hiD, &
       flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
       flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
       flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
       dx,dy,dz)

    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QREINT, QPRES, QFS, &
         NVAR, UMX, UMY, UMZ, UEDEN, UEINT, UFS
    use chemistry_module, only : nspec=>nspecies

    implicit none

    integer, intent(in) :: lof(3), hif(3), loD(3), hiD(3)
    double precision, intent(in) :: dx, dy, dz
    integer,intent(in)::    q_l1,    q_l2,    q_l3,    q_h1,    q_h2,    q_h3
    integer,intent(in)::flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
    integer,intent(in)::flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
    integer,intent(in)::flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
    double precision,intent(in   )::rhoYD(loD(1):hiD(1),loD(2):hiD(2),loD(3):hiD(3),nspec)
    double precision,intent(in   )::lamcp(loD(1):hiD(1),loD(2):hiD(2),loD(3):hiD(3))
    double precision,intent(in   )::   mu(loD(1):hiD(1),loD(2):hiD(2),loD(3):hiD(3))
    double precision,intent(in   )::    q(    q_l1:    q_h1,    q_l2:    q_h2,    q_l3:    q_h3,QVAR)
    double precision,intent(inout)::flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision,intent(inout)::flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision,intent(inout)::flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)

    integer :: i, j, k, n
    double precision :: divu, dxinv, dyinv, dzinv
    double precision :: dudx, dvdy, dwdz, dudy, dudz, dvdx, dvdz, dwdx, dwdy
    double precision :: tauxx, tauxy, tauxz
    double precision :: tauyy, tauyx, tauyz
    double precision :: tauzz, tauzx, tauzy
    double precision :: dhdx, dhdy, dhdz
    double precision :: mum, rhoeflux
    double precision, allocatable :: h(:,:,:)

    double precision, parameter :: two3rd = 2.d0/3.d0

    dxinv = 1.d0/dx
    dyinv = 1.d0/dy
    dzinv = 1.d0/dz

    allocate(h(lof(1)-1:hif(1)+1, lof(2)-1:hif(2)+1, lof(3)-1:hif(3)+1))

    do k = lof(3)-1, hif(3)+1
    do j = lof(2)-1, hif(2)+1
    do i = lof(1)-1, hif(1)+1
       h(i,j,k) = (q(i,j,k,QREINT) + q(i,j,k,QPRES)) / q(i,j,k,QRHO)
    end do
    end do
    end do

    ! x-direction
    do k = lof(3), hif(3)
    do j = lof(2), hif(2)
    do i = lof(1), hif(1)+1

       dudx =        ( q(i, j, k,  QU) - q(i-1,j,k,QU)) * dxinv
       dudy = 0.25d0*( q(i-1,j+1,k,QU) + q(i,j+1,k,QU) &
            &        - q(i-1,j-1,k,QU) - q(i,j-1,k,QU)) * dyinv
       dudz = 0.25d0*( q(i-1,j,k+1,QU) + q(i,j,k+1,QU) &
            &        - q(i-1,j,k-1,QU) - q(i,j,k-1,QU)) * dzinv

       dvdx =        ( q(i, j, k,  QV) - q(i-1,j,k,QV)) * dxinv
       dvdy = 0.25d0*( q(i-1,j+1,k,QV) + q(i,j+1,k,QV) &
            &        - q(i-1,j-1,k,QV) - q(i,j-1,k,QV)) * dyinv

       dwdx =        ( q(i, j, k,  QW) - q(i-1,j,k,QW)) * dxinv
       dwdz = 0.25d0*( q(i-1,j,k+1,QW) + q(i,j,k+1,QW) &
            &        - q(i-1,j,k-1,QW) - q(i,j,k-1,QW)) * dzinv

       divu = dudx + dvdy + dwdz

       mum = 0.5d0 * (mu(i-1,j,k)+mu(i,j,k))

       tauxx = mum*(2.d0*dudx - two3rd*divu)       
       tauxy = mum*(dudy + dvdx)
       tauxz = mum*(dudz + dwdx)

       flux1(i,j,k,UMX) = flux1(i,j,k,UMX) - tauxx
       flux1(i,j,k,UMY) = flux1(i,j,k,UMY) - tauxy
       flux1(i,j,k,UMZ) = flux1(i,j,k,UMZ) - tauxz
       flux1(i,j,k,UEDEN) = flux1(i,j,k,UEDEN) - &
            ( tauxx * 0.5d0*(q(i-1,j,k,QU)+q(i,j,k,QU)) &
            + tauxy * 0.5d0*(q(i-1,j,k,QV)+q(i,j,k,QV)) &
            + tauxz * 0.5d0*(q(i-1,j,k,QW)+q(i,j,k,QW)) )

       dhdx = (h(i,j,k)-h(i-1,j,k))*dxinv
       rhoeflux = 0.5d0*(lamcp(i-1,j,k)+lamcp(i,j,k))*dhdx
       flux1(i,j,k,UEDEN) = flux1(i,j,k,UEDEN) + rhoeflux
       flux1(i,j,k,UEINT) = flux1(i,j,k,UEINT) + rhoeflux
    end do
    end do
    end do

    ! x-direction
    do n = 1, nspec
       do k = lof(3), hif(3)
       do j = lof(2), hif(2)
       do i = lof(1), hif(1)+1
          flux1(i,j,k,UFS+n-1) = flux1(i,j,k,UFS+n-1) &
               - 0.5d0*(rhoYD(i-1,j,k,n) + rhoYD(i,j,k,n)) &
               * (q(i,j,k,QFS+n-1) - q(i-1,j,k,QFS+n-1)) * dxinv
       end do
       end do
       end do
    end do

    ! y-direction
    do k = lof(3), hif(3)
    do j = lof(2), hif(2)+1
    do i = lof(1), hif(1)

       dudx = 0.25d0*( q(i+1,j-1,k,QU) + q(i+1,j,k,QU) &
            &        - q(i-1,j-1,k,QU) - q(i-1,j,k,QU)) * dxinv
       dudy =        ( q(i, j, k,  QU) - q(i,j-1,k,QU)) * dyinv

       dvdx = 0.25d0*( q(i+1,j-1,k,QV) + q(i+1,j,k,QV) &
            &        - q(i-1,j-1,k,QV) - q(i-1,j,k,QV)) * dxinv
       dvdy =        ( q(i, j, k,  QV) - q(i,j-1,k,QV)) * dyinv
       dvdz = 0.25d0*( q(i,j-1,k+1,QV) + q(i,j,k+1,QV) &
            &        - q(i,j-1,k-1,QV) - q(i,j,k-1,QV)) * dzinv

       dwdy =        ( q(i, j, k,  QW) - q(i,j-1,k,QW)) * dyinv
       dwdz = 0.25d0*( q(i,j-1,k+1,QW) + q(i,j,k+1,QW) &
            &        - q(i,j-1,k-1,QW) - q(i,j,k-1,QW)) * dzinv

       divu = dudx + dvdy + dwdz

       mum = 0.5d0 * (mu(i,j-1,k)+mu(i,j,k))

       tauyx = mum*(dvdx + dudy)
       tauyy = mum*(2.d0*dvdy - two3rd*divu)
       tauyz = mum*(dvdz + dwdy)

       flux2(i,j,k,UMX) = flux2(i,j,k,UMX) - tauyx
       flux2(i,j,k,UMY) = flux2(i,j,k,UMY) - tauyy
       flux2(i,j,k,UMZ) = flux2(i,j,k,UMZ) - tauyz
       flux2(i,j,k,UEDEN) = flux2(i,j,k,UEDEN) - &
            ( tauyx * 0.5d0*(q(i,j-1,k,QU)+q(i,j,k,QU)) &
            + tauyy * 0.5d0*(q(i,j-1,k,QV)+q(i,j,k,QV)) &
            + tauyz * 0.5d0*(q(i,j-1,k,QW)+q(i,j,k,QW)) )

       dhdy = (h(i,j,k)-h(i,j-1,k))*dyinv
       rhoeflux = 0.5d0*(lamcp(i,j-1,k)+lamcp(i,j,k))*dhdy
       flux2(i,j,k,UEDEN) = flux2(i,j,k,UEDEN) + rhoeflux
       flux2(i,j,k,UEINT) = flux2(i,j,k,UEINT) + rhoeflux
    end do
    end do
    end do

    ! y-direction
    do n = 1, nspec
       do k = lof(3), hif(3)
       do j = lof(2), hif(2)+1
       do i = lof(1), hif(1)
          flux2(i,j,k,UFS+n-1) = flux2(i,j,k,UFS+n-1) &
               - 0.5d0*(rhoYD(i,j-1,k,n) + rhoYD(i,j,k,n)) &
               * (q(i,j,k,QFS+n-1) - q(i,j-1,k,QFS+n-1)) * dyinv
       end do
       end do
       end do
    end do

    ! z-direction
    do k = lof(3), hif(3)+1
    do j = lof(2), hif(2)
    do i = lof(1), hif(1)

       dudx = 0.25d0*( q(i+1,j,k-1,QU) + q(i+1,j,k,QU) &
            &        - q(i-1,j,k-1,QU) - q(i-1,j,k,QU)) * dxinv
       dudz =        ( q(i, j, k,  QU) - q(i,j,k-1,QU)) * dzinv

       dvdy = 0.25d0*( q(i,j+1,k-1,QV) + q(i,j+1,k,QV) &
            &        - q(i,j-1,k-1,QV) - q(i,j-1,k,QV)) * dyinv
       dvdz =        ( q(i, j, k,  QV) - q(i,j,k-1,QV)) * dzinv

       dwdx = 0.25d0*( q(i+1,j,k-1,QW) + q(i+1,j,k,QW) &
            &        - q(i-1,j,k-1,QW) - q(i-1,j,k,QW)) * dxinv
       dwdy = 0.25d0*( q(i,j+1,k-1,QW) + q(i,j+1,k,QW) &
            &        - q(i,j-1,k-1,QW) - q(i,j-1,k,QW)) * dyinv
       dwdz =        ( q(i, j, k,  QW) - q(i,j,k-1,QW)) * dzinv

       divu = dudx + dvdy + dwdz

       mum = 0.5d0 * (mu(i,j,k-1)+mu(i,j,k))

       tauzx = mum*(dwdx + dudz)
       tauzy = mum*(dwdy + dvdz)
       tauzz = mum*(2.d0*dwdz - two3rd*divu)

       flux3(i,j,k,UMX) = flux3(i,j,k,UMX) - tauzx
       flux3(i,j,k,UMY) = flux3(i,j,k,UMY) - tauzy
       flux3(i,j,k,UMZ) = flux3(i,j,k,UMZ) - tauzz
       flux3(i,j,k,UEDEN) = flux3(i,j,k,UEDEN) - &
            ( tauzx * 0.5d0*(q(i,j,k-1,QU)+q(i,j,k,QU)) &
            + tauzy * 0.5d0*(q(i,j,k-1,QV)+q(i,j,k,QV)) &
            + tauzz * 0.5d0*(q(i,j,k-1,QW)+q(i,j,k,QW)) )

       dhdz = (h(i,j,k)-h(i,j,k-1))*dzinv
       rhoeflux = 0.5d0*(lamcp(i,j,k-1)+lamcp(i,j,k))*dhdz
       flux3(i,j,k,UEDEN) = flux3(i,j,k,UEDEN) + rhoeflux
       flux3(i,j,k,UEINT) = flux3(i,j,k,UEINT) + rhoeflux
    end do
    end do
    end do

    ! z-direction
    do n = 1, nspec
       do k = lof(3), hif(3)+1
       do j = lof(2), hif(2)
       do i = lof(1), hif(1)
          flux3(i,j,k,UFS+n-1) = flux3(i,j,k,UFS+n-1) &
               - 0.5d0*(rhoYD(i,j,k-1,n) + rhoYD(i,j,k,n)) &
               * (q(i,j,k,QFS+n-1) - q(i,j,k-1,QFS+n-1)) * dzinv
       end do
       end do
       end do
    end do

    deallocate(h)

  end subroutine comp_dflux

! ::
! :: ----------------------------------------------------------
! ::

  subroutine transCoef(lo,hi,q,TEMP,CP, &
       q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
       D, &
       D_l1,D_l2,D_l3,D_h1,D_h2,D_h3)

    use meth_params_module


      implicit none

!      include "cdwrk.H"
      ! xxxxxxxxxx
      integer, parameter :: Nspec = 0

      integer          :: lo(3), hi(3)
      integer          :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer          :: D_l1, D_l2, D_l3, D_h1, D_h2, D_h3
      double precision ::    q(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3,QVAR)
      double precision :: TEMP(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3)
      double precision ::   CP(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3) 
      double precision ::    D(D_l1:D_h1, D_l2:D_h2, D_l3:D_h3,Nspec+3)

      ! Local variables
      integer          :: i,j,k,doVelVisc,doTemp
!     double precision :: Dmax

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               TEMP(i,j,k) = MIN(MAX(q(i,j,k,QTHERM+1),300.d0),6500.d0)
            end do
         end do
      end do
      
      doTemp = 1
      doVelVisc = 1

      call mixavg_rhodiff_temp_pres(lo, hi, &
           D,  &
             D_l1,D_l2,D_l3,D_h1,D_h2,D_h3, &
           TEMP, &
             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
           q(q_l1,q_l2,q_l3,QTHERM+NADV+1), &
             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
           q(q_l1,q_l2,q_l3,QPRES), &
             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
           doTemp, doVelVisc)

!     call FORT_CPMIXfromTY(lo, hi, &
!      call dcpmty(lo, hi, &
!           CP, &
!             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
!           TEMP, &
!             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
!           q(q_l1,q_l2,q_l3,QTHERM+NADV+1), &
!             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3)

      ! Replace lambda with lambda/Cp
      D(:,:,:,Nspec+1) = D(:,:,:,Nspec+1)/CP(:,:,:)

!     Dmax = -1.d0
!     do k = lo(3),hi(3)
!        do j = lo(2),hi(2)
!           do i = lo(1),hi(1)
!              do n=1,Nspec
!                 Dmax = MAX(Dmax,D(i,j,k,n) / q(i,j,k,QRHO))
!              enddo
!           end do
!        end do
!     end do
!     print*, '********************      Dmax: ', Dmax

      end subroutine transCoef
! ::
! :: ----------------------------------------------------------
! ::
      subroutine myViscFlux(lo,hi,&
                            q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                            D,D_l1,D_l2,D_l3,D_h1,D_h2,D_h3, &
                            flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                            flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                            flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                            dx,dy,dz,dt)

        use meth_params_module

      implicit none

!      include "cdwrk.H"

! xxxxxx
      integer, parameter :: NSPEC = 0
      integer, parameter :: maxspec = 0

      integer          :: lo(3),hi(3)
      integer          :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer          :: D_l1, D_l2, D_l3, D_h1, D_h2, D_h3
      integer          :: flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3
      integer          :: flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3
      integer          :: flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3
      double precision :: q(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3, QVAR)
      double precision :: D(D_l1:D_h1, D_l2:D_h2, D_l3:D_h3, Nspec+3)
      double precision :: flux1(flux1_l1:flux1_h1, flux1_l2:flux1_h2, flux1_l3:flux1_h3, NVAR)
      double precision :: flux2(flux2_l1:flux2_h1, flux2_l2:flux2_h2, flux2_l3:flux2_h3, NVAR)
      double precision :: flux3(flux3_l1:flux3_h1, flux3_l2:flux3_h2, flux3_l3:flux3_h3, NVAR)

      double precision :: dx,dy,dz,dt

      ! Local variables
      integer          :: i,j,k,n,ifirstSp
      double precision :: tauxxm,tauxym,tauxzm
      double precision :: tauyym,tauyxm,tauyzm
      double precision :: tauzzm,tauzxm,tauzym
      double precision :: divxm,divym,divzm
      double precision :: muxm,muym,muzm
      double precision :: lamxm,lamym,lamzm
      double precision :: kxm,kym,kzm
      double precision :: phiflx

      double precision :: rhoDm(maxspec)
      double precision, parameter ::    two3rd = 2.d0/3.d0

      ifirstSp  = NTHERM+NADV+1

      ! compute x fluxes
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               
               muxm =  0.5d0*(D(i-1,j,k,Nspec+2) + D(i,j,k,Nspec+2))
               lamxm = -two3rd*muxm
               kxm =   0.5d0*(D(i-1,j,k,Nspec+1) + D(i,j,k,Nspec+1))
               do n=1,Nspec
                  rhoDm(n) = 0.5d0*(D(i-1,j,k,n) + D(i,j,k,n))
               end do
               
               tauxxm = 2.d0*muxm*(q(i,j,k,QU) - q(i-1,j,k,QU))/dx 
               
               divxm = lamxm*(q(i,j,k,QU) - q(i-1,j,k,QU))/dx  &
                    +.25d0*lamxm*(q(i-1,j+1,k,QV)+q(i,j+1,k,QV) &
                    -             q(i-1,j-1,k,QV)-q(i,j-1,k,QV))/dy &
                    +.25d0*lamxm*(q(i-1,j,k+1,QW)+q(i,j,k+1,QW) &
                    -             q(i-1,j,k-1,QW)-q(i,j,k-1,QW))/dz
               
               tauyxm =  muxm*(q(i,j,k,QV) - q(i-1,j,k,QV))/dx  &
                    +.25d0*muxm*(q(i,j+1,k,QU)+q(i-1,j+1,k,QU) &
                    -            q(i,j-1,k,QU)-q(i-1,j-1,k,QU))/dy
               
               tauzxm =  muxm*(q(i,j,k,QW) - q(i-1,j,k,QW))/dx &
                    +.25d0*muxm*(q(i,j,k+1,QU)+q(i-1,j,k+1,QU)&
                    -            q(i,j,k-1,QU)-q(i-1,j,k-1,QU))/dz

               phiflx =  tauxxm*(q(i,j,k,QU)+q(i-1,j,k,QU)) &
                    +     divxm*(q(i,j,k,QU)+q(i-1,j,k,QU)) &
                    +    tauyxm*(q(i,j,k,QV)+q(i-1,j,k,QV)) &
                    +    tauzxm*(q(i,j,k,QW)+q(i-1,j,k,QW))
               
               flux1(i,j,k,2) = flux1(i,j,k,2) - dt*dy*dz*(tauxxm+divxm)
               flux1(i,j,k,3) = flux1(i,j,k,3) - dt*dy*dz*tauyxm
               flux1(i,j,k,4) = flux1(i,j,k,4) - dt*dy*dz*tauzxm
               flux1(i,j,k,5) = flux1(i,j,k,5) - dt*dy*dz*(0.5d0*phiflx  &
                    + kxm*( (q(i  ,j,k,QREINT)+q(i  ,j,k,QPRES))/q(i  ,j,k,QRHO) &
                    -       (q(i-1,j,k,QREINT)+q(i-1,j,k,QPRES))/q(i-1,j,k,QRHO) )/dx)

               do n=1,Nspec
                  flux1(i,j,k,n+ifirstSp-1) = flux1(i,j,k,n+ifirstSp-1) &
                       - dt*dy*dz*rhoDm(n)* &
                       (q(i,j,k,QTHERM+NADV+n)-q(i-1,j,k,QTHERM+NADV+n))/dx
               end do
               
            enddo
         enddo
      enddo

      ! compute y fluxes
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)
               
               muym =  0.5d0*(D(i,j-1,k,Nspec+2) + D(i,j,k,Nspec+2))
               lamym = -two3rd*muym
               kym =   0.5d0*(D(i,j-1,k,Nspec+1) + D(i,j,k,Nspec+1))
               do n=1,Nspec
                  rhoDm(n) = 0.5d0*(D(i,j-1,k,n) + D(i,j,k,n))
               end do
               
               tauxym =  muym*(q(i,j,k,QU) - q(i,j-1,k,QU))/dy  &
                    +.25d0*muym*(q(i+1,j,k,QV)+q(i+1,j-1,k,QV) &
                    -            q(i-1,j,k,QV)-q(i-1,j-1,k,QV))/dx
               
               tauyym = 2.d0*muym*(q(i,j,k,QV) - q(i,j-1,k,QV))/dy 
               
               divym = lamym*(q(i,j,k,QV) - q(i,j-1,k,QV))/dy  &
                    +.25d0*lamym*(q(i+1,j-1,k,QU)+q(i+1,j,k,QU) &
                    -             q(i-1,j-1,k,QU)-q(i-1,j,k,QU))/dx &
                    +.25d0*lamym*(q(i,j-1,k+1,QW)+q(i,j,k+1,QW) &
                    -             q(i,j-1,k-1,QW)-q(i,j,k-1,QW))/dz
               
               tauzym =  muym*(q(i,j,k,QW) - q(i,j-1,k,QW))/dy  &
                    +.25d0*muym*(q(i,j,k+1,QV)+q(i,j-1,k+1,QV) &
                    -            q(i,j,k-1,QV)-q(i,j-1,k-1,QV))/dz

               phiflx = tauxym*(q(i,j,k,QU)+q(i,j-1,k,QU)) &
                    +   tauyym*(q(i,j,k,QV)+q(i,j-1,k,QV)) &
                    +    divym*(q(i,j,k,QV)+q(i,j-1,k,QV)) &
                    +   tauzym*(q(i,j,k,QW)+q(i,j-1,k,QW))
               
               flux2(i,j,k,2) = flux2(i,j,k,2) - dt*dx*dz*tauxym
               flux2(i,j,k,3) = flux2(i,j,k,3) - dt*dx*dz*(tauyym+divym)
               flux2(i,j,k,4) = flux2(i,j,k,4) - dt*dx*dz*tauzym
               flux2(i,j,k,5) = flux2(i,j,k,5) - dt*dy*dz*(0.5d0*phiflx  &
                    + 2.d0*kym*( (q(i,j  ,k,QREINT)+q(i,j  ,k,QPRES))/q(i,j  ,k,QRHO) &
                    -           (q(i,j-1,k,QREINT)+q(i,j-1,k,QPRES))/q(i,j-1,k,QRHO) )/dy)

               do n=1,Nspec
                  flux2(i,j,k,n+ifirstSp-1) = flux2(i,j,k,n+ifirstSp-1) &
                       - dt*dx*dz*rhoDm(n)* &
                       (q(i,j,k,QTHERM+NADV+n)-q(i,j-1,k,QTHERM+NADV+n))/dy
               end do
               
            enddo
         enddo
      enddo

      ! compute z fluxes
      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               
               muzm =  0.5d0*(D(i,j,k-1,Nspec+2) + D(i,j,k,Nspec+2))
               lamzm = -two3rd*muzm
               kzm =   0.5d0*(D(i,j,k-1,Nspec+1) + D(i,j,k,Nspec+1))
               do n=1,Nspec
                  rhoDm(n) = 0.5d0*(D(i,j,k-1,n) + D(i,j,k,n))
               end do
               
               tauxzm =  muzm*(q(i,j,k,QU) - q(i,j,k-1,QU))/dz  &
                    +.25d0*muzm*(q(i+1,j,k,QW)+q(i+1,j,k-1,QW) &
                    -            q(i-1,j,k,QW)-q(i-1,j,k-1,QW))/dx

               tauyzm =  muzm*(q(i,j,k,QV) - q(i,j,k-1,QV))/dz  &
                    +.25d0*muzm*(q(i,j+1,k,QW)+q(i,j+1,k-1,QW) &
                    -            q(i,j-1,k,QW)-q(i,j-1,k-1,QW))/dy
               
               tauzzm = 2.d0*muzm*(q(i,j,k,QW) - q(i,j,k-1,QW))/dz 

               divzm = lamzm*(q(i,j,k,QW) - q(i,j,k-1,QW))/dz  &
                    +.25d0*lamzm*(q(i+1,j,k-1,QU)+q(i+1,j,k,QU) &
                    -             q(i-1,j,k-1,QU)-q(i-1,j,k,QU))/dx &
                    +.25d0*lamzm*(q(i,j+1,k-1,QV)+q(i,j+1,k,QV) &
                    -             q(i,j-1,k-1,QV)-q(i,j-1,k,QV))/dy
               
               phiflx =tauxzm*(q(i,j,k,QU)+q(i,j,k-1,QU)) &
                    +  tauyzm*(q(i,j,k,QV)+q(i,j,k-1,QV)) &
                    +  tauzzm*(q(i,j,k,QW)+q(i,j,k-1,QW)) &
                    +   divzm*(q(i,j,k,QW)+q(i,j,k-1,QW))

               flux3(i,j,k,2) = flux3(i,j,k,2) - dt*dx*dy*tauxzm
               flux3(i,j,k,3) = flux3(i,j,k,3) - dt*dx*dy*tauyzm
               flux3(i,j,k,4) = flux3(i,j,k,4) - dt*dx*dy*(tauzzm+divzm)
               flux3(i,j,k,5) = flux3(i,j,k,5) - dt*dx*dy*(0.5d0*phiflx &
                    + 2.d0*kzm*( (q(i,j,k  ,QREINT)+q(i,j,k  ,QPRES))/q(i,j,k  ,QRHO) &
                    -            (q(i,j,k-1,QREINT)+q(i,j,k-1,QPRES))/q(i,j,k-1,QRHO) )/dz)

               do n=1,Nspec
                  flux3(i,j,k,n+ifirstSp -1) = flux3(i,j,k,n+ifirstSp -1) &
                       - dt*dx*dy*rhoDm(n)* &
                       (q(i,j,k,QTHERM+NADV+n)-q(i,j,k-1,QTHERM+NADV+n))/dz
               end do
               
            enddo
         enddo
      enddo

      end subroutine myViscFlux
! ::
! :: ----------------------------------------------------------
! ::
      subroutine mixavg_rhodiff_temp_pres(lo, hi, &
                          rd,rd_l1,rd_l2,rd_l3,rd_h1,rd_h2,rd_h3, &
                           T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, &
                           Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3, &
                           P, P_l1, P_l2, P_l3, P_h1, P_h2, P_h3, &
                           do_temp, do_VelVisc)

        use meth_params_module

      implicit none

!      include "cdwrk.H"
!xxxxx
      integer, parameter :: maxspec = 0
      integer, parameter :: Nspec = 0

      integer          :: lo(3),hi(3)
      integer          :: do_temp, do_VelVisc
      integer          :: rd_l1,rd_l2,rd_l3,rd_h1,rd_h2,rd_h3
      integer          ::  T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer          ::  Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      integer          ::  P_l1, P_l2, P_l3, P_h1, P_h2, P_h3
      double precision :: rd(rd_l1:rd_h1, rd_l2:rd_h2, rd_l3:rd_h3,*)
      double precision ::  T( T_l1: T_h1,  T_l2: T_h2,  T_l3: T_h3)
      double precision ::  Y( Y_l1: Y_h1,  Y_l2: Y_h2,  Y_l3: Y_h3,*)
      double precision ::  P( P_l1: P_h1,  P_l2: P_h2,  P_l3: P_h3)

      ! Local variables
      integer          :: i, j, k, n
      double precision :: Yt(maxspec), Dt(maxspec)
      double precision :: scal, tscal, Wavg, RHO, Tt, Pt, invmwt(maxspec)
      double precision :: alpha, l1, l2, X(maxspec), CPMS(maxspec)

!      call CKWT(IWRK(ckbi),RWRK(ckbr),invmwt)

      do n=1,Nspec
         invmwt(n) = 1.d0 / invmwt(n)
      end do

      scal  = 0.1d0
      tscal = 1.d-5
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
!xxx               Tt = MAX(T(i,j,k),TMIN_TRANS)
               Pt = P(i,j,k)*1.d1
!               CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
!               CALL CKCPMS(Tt,IWRK(ckbi),RWRK(ckbr),CPMS)
!               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
!wqz               CALL EGSPAR(Tt,X,Yt,CPMS,RWRK(egbr),IWRK(egbi))
!               CALL EGSPAR(Tt,X,Yt,CPMS,RWRK,IWRK)
!wqz               CALL EGSV1(Pt,Tt,Yt,Wavg,RWRK(egbr),Dt)
!               CALL EGSV1(Pt,Tt,Yt,Wavg,RWRK,Dt)
!               CALL CKRHOY(Pt,Tt,Yt,IWRK(ckbi),RWRK(ckbr),RHO)
               do n=1,Nspec
                  rd(i,j,k,n) = RHO * Wavg * invmwt(n) * Dt(n) * scal
               end do

               if (do_temp .ne. 0) then
                  alpha = 1
!wqz                  CALL EGSL1(alpha,Tt,X,RWRK(egbr),l1)
!                  CALL EGSL1(alpha,Tt,X,RWRK,l1)
                  alpha = -1
!wqz                  CALL EGSL1(alpha,Tt,X,RWRK(egbr),l2)
!                  CALL EGSL1(alpha,Tt,X,RWRK,l2)
                  rd(i,j,k,Nspec+1) = 0.5d0 * (l1 + l2) * tscal
               endif

               if (do_VelVisc .ne. 0) then
!                  CALL EGSE3(Tt,Yt,RWRK(egbr),rd(i,j,k,Nspec+2))
!                  CALL EGSE3(Tt,Yt,RWRK,rd(i,j,k,Nspec+2))
                  rd(i,j,k,Nspec+2) = rd(i,j,k,Nspec+2) * scal
               endif

            end do
         end do
      end do

      end subroutine mixavg_rhodiff_temp_pres


  subroutine fluxtosrc(lo_work, hi_work, &
       src  , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
       flux1,f1_l1,f1_l2,f1_l3,f1_h1,f1_h2,f1_h3, &
       flux2,f2_l1,f2_l2,f2_l3,f2_h1,f2_h2,f2_h3, &
       flux3,f3_l1,f3_l2,f3_l3,f3_h1,f3_h2,f3_h3, &
       area1,a1_l1,a1_l2,a1_l3,a1_h1,a1_h2,a1_h3, &
       area2,a2_l1,a2_l2,a2_l3,a2_h1,a2_h2,a2_h3, &
       area3,a3_l1,a3_l2,a3_l3,a3_h1,a3_h2,a3_h3, &
       vol  , v_l1, v_l2, v_l3, v_h1, v_h2, v_h3)

    use meth_params_module, only : NVAR, UTEMP

    integer, intent(in) :: lo_work(3), hi_work(3)
    integer, intent(in) :: f1_l1,f1_l2,f1_l3,f1_h1,f1_h2,f1_h3
    integer, intent(in) :: f2_l1,f2_l2,f2_l3,f2_h1,f2_h2,f2_h3
    integer, intent(in) :: f3_l1,f3_l2,f3_l3,f3_h1,f3_h2,f3_h3
    integer, intent(in) :: a1_l1,a1_l2,a1_l3,a1_h1,a1_h2,a1_h3
    integer, intent(in) :: a2_l1,a2_l2,a2_l3,a2_h1,a2_h2,a2_h3
    integer, intent(in) :: a3_l1,a3_l2,a3_l3,a3_h1,a3_h2,a3_h3
    integer, intent(in) ::  v_l1, v_l2, v_l3, v_h1, v_h2, v_h3
    integer, intent(in) ::  s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    double precision, intent(in) :: flux1(f1_l1:f1_h1,f1_l2:f1_h2,f1_l3:f1_h3,NVAR)
    double precision, intent(in) :: flux2(f2_l1:f2_h1,f2_l2:f2_h2,f2_l3:f2_h3,NVAR)
    double precision, intent(in) :: flux3(f3_l1:f3_h1,f3_l2:f3_h2,f3_l3:f3_h3,NVAR)
    double precision, intent(in) :: area1(a1_l1:a1_h1,a1_l2:a1_h2,a1_l3:a1_h3)
    double precision, intent(in) :: area2(a2_l1:a2_h1,a2_l2:a2_h2,a2_l3:a2_h3)
    double precision, intent(in) :: area3(a3_l1:a3_h1,a3_l2:a3_h2,a3_l3:a3_h3)
    double precision, intent(in) :: vol  ( v_l1: v_h1, v_l2: v_h2, v_l3: v_h3)
    double precision, intent(out):: src  ( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3,NVAR)

    integer :: i, j, k, n

    do n = 1, NVAR
       if (n .eq. UTEMP) then
          src(:,:,:,n) = 0.d0
       else
          !$omp parallel do private(i,j,k)
          do k       = lo_work(3), hi_work(3)
             do j    = lo_work(2), hi_work(2)
                do i = lo_work(1), hi_work(1)
                   
                   src(i,j,k,n) =  &
                        ( (flux1(i,j,k,n)*area1(i,j,k)-flux1(i+1,j,k,n)*area1(i+1,j,k)) &
                        + (flux2(i,j,k,n)*area2(i,j,k)-flux2(i,j+1,k,n)*area2(i,j+1,k)) &
                        + (flux3(i,j,k,n)*area3(i,j,k)-flux3(i,j,k+1,n)*area3(i,j,k+1))) &
                        / vol(i,j,k)
                end do
             end do
          end do
          !$omp end parallel do
       end if
    end do

  end subroutine fluxtosrc


  subroutine diffup(lo,hi, &
       uout,uout_l1,uout_l2,uout_l3, uout_h1,uout_h2,uout_h3, &
       src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
       flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
       flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
       flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
       area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
       area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
       area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
       vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
       dx,dy,dz,dt_flux,dt_src)

    use advection_module, only : normalize_species_fluxes
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UTEMP, normalize_species

    implicit none

    integer,intent(in):: lo(3), hi(3)
    integer,intent(in)::  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
    integer,intent(in)::   src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3 
    integer,intent(in):: flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
    integer,intent(in):: flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
    integer,intent(in):: flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
    integer,intent(in):: area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3
    integer,intent(in):: area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3
    integer,intent(in):: area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3
    integer,intent(in):: vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3

    double precision,intent(inout):: uout( uout_l1: uout_h1, uout_l2: uout_h2, uout_l3: uout_h3,NVAR)
    double precision,intent(in   )::  src(  src_l1:  src_h1,  src_l2:  src_h2,  src_l3:  src_h3,NVAR)
    double precision,intent(inout)::flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision,intent(inout)::flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision,intent(inout)::flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
    double precision,intent(in   )::area1(area1_l1:area1_h1,area1_l2:area1_h2,area1_l3:area1_h3)
    double precision,intent(in   )::area2(area2_l1:area2_h1,area2_l2:area2_h2,area2_l3:area2_h3)
    double precision,intent(in   )::area3(area3_l1:area3_h1,area3_l2:area3_h2,area3_l3:area3_h3)
    double precision,intent(in   )::  vol(vol_l1:vol_h1,vol_l2:vol_h2,vol_l3:vol_h3)
    double precision,intent(in) :: dx, dy, dz, dt_flux, dt_src

    integer :: i,j,k,n

    do n = 1, NVAR
         
       if ( n.eq.UTEMP ) then
          
          flux1(:,:,:,n) = 0.d0
          flux2(:,:,:,n) = 0.d0
          flux3(:,:,:,n) = 0.d0
          
       else
          
          !$OMP PARALLEL PRIVATE(i,j,k,div1)
          !$OMP DO
          do k       = lo(3),hi(3)
             do j    = lo(2),hi(2)
                do i = lo(1),hi(1)+1
                   flux1(i,j,k,n) = flux1(i,j,k,n) * area1(i,j,k) * dt_flux
                enddo
             enddo
          enddo
          !$OMP END DO NOWAIT
          !$OMP DO
          do k       = lo(3),hi(3)
             do j    = lo(2),hi(2)+1
                do i = lo(1),hi(1)
                   flux2(i,j,k,n) = flux2(i,j,k,n) * area2(i,j,k) * dt_flux
                enddo
             enddo
          enddo
          !$OMP END DO NOWAIT
          !$OMP DO
          do k       = lo(3),hi(3)+1
             do j    = lo(2),hi(2)
                do i = lo(1),hi(1)
                   flux3(i,j,k,n) = flux3(i,j,k,n) * area3(i,j,k) * dt_flux
                enddo
             enddo
          enddo
          !$OMP END DO
          !$OMP END PARALLEL
          
       endif
       
    enddo

    if (normalize_species .eq. 1) &
         call normalize_species_fluxes( &
         flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
         flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
         flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
         lo,hi)
    
    do n = 1, NVAR
       
       ! pass temperature through
       if (n .ne. UTEMP) then
          ! update everything else with fluxes and source terms
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k       = lo(3),hi(3)
             do j    = lo(2),hi(2)
                do i = lo(1),hi(1)
                   uout(i,j,k,n) = uout(i,j,k,n) &
                        + ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                        +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                        +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) / vol(i,j,k) &
                        +   dt_src * src(i,j,k,n)
                  enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         endif
         
      enddo

    end subroutine diffup

end module diff_flux_module
