module probdata_module

  ! adiabatic index of the gas
  double precision, parameter :: gamma_gas = 1.4d0

  ! unshocked gas
  double precision, parameter :: rho0 = 1.4d0
  double precision, parameter ::   p0 = 1.0d0
  double precision, parameter ::   v0 = 0.0d0

  ! shocked gas
  double precision, parameter :: rho1 = 8.d0
  double precision, parameter ::   p1 = 116.5d0
  double precision, parameter ::   v1 = 8.25d0

  ! shock velocity
  double precision, parameter :: vshock = 10.0d0

  ! Shock angle: 60 degree
  double precision, parameter :: thetashock = &
       3.141592653589793238462643383279502884197d0 / 3.d0

  ! Initial shock position at YLO
  double precision, parameter :: xshock0_lo = 1.d0/6.d0

  double precision, save :: xshock0_hi, vshock_x
  double precision, allocatable, save :: state1(:), state0(:)

end module probdata_module
