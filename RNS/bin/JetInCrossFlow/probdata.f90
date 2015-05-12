module probdata_module

  implicit none

  double precision, save :: pamb
  double precision, save :: v_cf, T_cf
  double precision, save :: v_jet, T_jet
  double precision, save :: X_H2_jet
  double precision, save :: r_jet

  ! boundary
  logical, save :: probdata_initialized = .false.
  double precision, allocatable, save :: state_cf(:), state_jet(:)

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

  subroutine init_probdata()
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS, NSPEC
    use chemistry_module, only : get_species_index

    integer :: n, iwrk
    integer :: iN2, iO2, iH2
    double precision :: rwrk, Xt(nspec), Yt(nspec), rhot, et

    iN2 = get_species_index("N2")
    iO2 = get_species_index("O2")
    iH2 = get_species_index("H2")

    allocate(state_cf(NVAR))
    allocate(state_jet(NVAR))

    ! cross flow
    Xt = 0.d0
    Xt(iN2) = 0.79d0
    Xt(iO2) = 0.21d0
    
    call CKXTY(Xt, iwrk, rwrk, Yt)
    call CKRHOX(pamb, T_cf, Xt, iwrk, rwrk, rhot)
    call CKUBMS(T_cf, Yt, iwrk, rwrk, et)
    
    state_cf(URHO)  = rhot
    state_cf(UMX)   = rhot*v_cf
    state_cf(UMY)   = 0.d0
    state_cf(UMY)   = 0.d0
    state_cf(UEDEN) = rhot*(et+0.5d0*v_cf**2)
    state_cf(UTEMP) = T_cf
    state_cf(UFS:UFS+nspec-1) = rhot*Yt

    ! jet
    Xt = 0.d0
    Xt(iH2) = 0.70d0
    Xt(iN2) = 0.30d0
    
    call CKXTY(Xt, iwrk, rwrk, Yt)
    call CKRHOX(pamb, T_jet, Xt, iwrk, rwrk, rhot)
    call CKUBMS(T_jet, Yt, iwrk, rwrk, et)
    
    state_jet(URHO)  = rhot
    state_jet(UMX)   = 0.d0
    state_jet(UMY)   = rhot*v_jet
    state_jet(UMY)   = 0.d0
    state_jet(UEDEN) = rhot*(et+0.5d0*v_jet**2)
    state_jet(UTEMP) = T_jet
    state_jet(UFS:UFS+nspec-1) = rhot*Yt

    probdata_initialized = .true.

  end subroutine init_probdata

end module probdata_module

