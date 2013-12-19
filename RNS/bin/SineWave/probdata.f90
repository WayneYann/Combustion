module probdata_module

  ! adiabatic index of the gas
  double precision, parameter :: gamma_gas = 1.4d0

  double precision, parameter :: p0     = 1.0d0
  double precision, parameter :: rho0   = 1.4d0
  double precision, parameter :: drho0  = 1.4d-1
  double precision, parameter :: u0     = 0.5d0

  integer, save :: idir = 0

  double precision, save :: center(3), length(3)

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

end module probdata_module
