module DME_jet_module

  use variables_module, only : irho, imx,imy,iene,iry1,ncons, iN2, iO2
  use chemistry_module, only : nspecies, get_species_index

  use probin_module, only : pamb, phi_in, T_in, vn_in, T_co, vn_co

  implicit none

  double precision, allocatable, save :: fuel_state(:), air_state(:)
  double precision, allocatable, save :: fuel_Y(:), air_Y(:)
  logical, save :: initialized = .false.

  private

  public :: fuel_state, air_state, fuel_Y, air_Y, initialized, init_DME_jet

contains

  subroutine init_DME_jet
    double precision, dimension(nspecies) :: Xt, Yt
    integer :: iwrk, iCH3OCH3, n
    double precision :: rwrk, Tt, rhot, Pt, et, vt, ek

    allocate(fuel_state(ncons))
    allocate( air_state(ncons))
    allocate(fuel_Y(nspecies))
    allocate( air_Y(nspecies))

    iCH3OCH3 = get_species_index("CH3OCH3")

    ! ----- Fuel -----
    Xt = 0.d0
    Xt(iCH3OCH3) = phi_in
    Xt(iN2) = 1.d0-Xt(iCH3OCH3)

    Pt = pamb
    Tt = T_in
    vt = vn_in

    call ckxty (Xt,iwrk,rwrk,Yt)
    call ckrhoy(Pt,Tt,Yt,iwrk,rwrk,rhot)
    call ckubms(Tt,Yt,iwrk,rwrk,et)

    ek = 0.5d0*vt**2

    fuel_state(irho) = rhot
    fuel_state(imx ) = rhot*vt
    fuel_state(imy ) = 0.d0
    fuel_state(iene) = rhot*(et+ek)
    do n=1,nspecies
       fuel_state(iry1+n-1) = rhot*Yt(n)
       fuel_Y(n) = Yt(n)
    end do

    ! ----- Air -----
    Xt = 0.d0
    Xt(iN2) = 0.79d0
    Xt(iO2) = 0.21d0

    Pt = pamb
    Tt = T_co
    vt = vn_co

    call ckxty (Xt,iwrk,rwrk,Yt)
    call ckrhoy(Pt,Tt,Yt,iwrk,rwrk,rhot)
    call ckubms(Tt,Yt,iwrk,rwrk,et)

    ek = 0.5d0*vt**2

    air_state(irho) = rhot
    air_state(imx ) = rhot*vt
    air_state(imy ) = 0.d0
    air_state(iene) = rhot*(et+ek)
    do n=1,nspecies
       air_state(iry1+n-1) = rhot*Yt(n)
       air_Y(n) = Yt(n)
    end do

    initialized = .true.
  end subroutine init_DME_jet

end module DME_jet_module

