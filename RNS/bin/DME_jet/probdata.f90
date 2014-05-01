module probdata_module

  implicit none

  ! parameters for DME_jet
  integer, save :: prob_type
  character (len=128), save :: turbfile = ""  
  double precision, save :: pamb, phi_in, T_in, vn_in, T_co, vn_co
  double precision, save :: splitx, xfrontw, Tfrontw
  double precision, save :: blobr, blobx, bloby, blobT
  double precision, save :: inflow_period, inflow_vnmag
  double precision, save :: splity, yfrontw
  double precision, allocatable, save :: fuel_Y(:), air_Y(:)
  double precision, allocatable, save :: fuel_state(:), air_state(:)

  logical, save :: dmejet_initialized = .false.

  double precision, save :: center(3)

  ! These determine the refinement criteria
  double precision, save :: denerr,   dengrad
  double precision, save :: velerr,   velgrad
  double precision, save :: presserr, pressgrad
  double precision, save :: temperr,  tempgrad
  double precision, save :: vorterr,  vortgrad
  double precision, save :: tracerr
  integer         , save :: max_denerr_lev    = -1
  integer         , save :: max_dengrad_lev   = -1
  integer         , save :: max_velerr_lev    = -1
  integer         , save :: max_velgrad_lev   = -1
  integer         , save :: max_presserr_lev  = -1 
  integer         , save :: max_pressgrad_lev = -1
  integer         , save :: max_temperr_lev   = -1
  integer         , save :: max_tempgrad_lev  = -1
  integer         , save :: max_vorterr_lev   = -1
  integer         , save :: max_vortgrad_lev  = -1
  integer         , save :: max_tracerr_lev   = -1

contains

  subroutine init_DME_jet

    use meth_params_module, only : ndim, NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS
    use chemistry_module, only : nspecies, get_species_index

    double precision, dimension(nspecies) :: Xt, Yt
    integer :: iwrk, iN2, iO2, iCH3OCH3, n
    double precision :: rwrk, Tt, rhot, Pt, et, vt, ek

    allocate(fuel_state(NVAR))
    allocate( air_state(NVAR))
    allocate(fuel_Y(nspecies))
    allocate( air_Y(nspecies))

    iN2 = get_species_index("N2")
    iO2 = get_species_index("O2")
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

    fuel_state(URHO ) = rhot
    fuel_state(UMX  ) = 0.d0
    if (ndim .eq. 3) then
       fuel_state(UMY  ) = 0.d0
       fuel_state(UMZ  ) = rhot*vt
    else
       fuel_state(UMY  ) = rhot*vt
    end if
    fuel_state(UEDEN) = rhot*(et+ek)
    do n=1,nspecies
       fuel_state(UFS+n-1) = rhot*Yt(n)
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

    air_state(URHO ) = rhot
    air_state(UMX  ) = 0.d0
    if (ndim .eq. 3) then
       air_state(UMY  ) = 0.d0
       air_state(UMZ  ) = rhot*vt
    else
       air_state(UMY  ) = rhot*vt
    end if
    air_state(UEDEN) = rhot*(et+ek)
    do n=1,nspecies
       air_state(UFS+n-1) = rhot*Yt(n)
       air_Y(n) = Yt(n)
    end do
    
    dmejet_initialized = .true.
  end subroutine init_DME_jet

end module probdata_module
