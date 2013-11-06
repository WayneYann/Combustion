module probin
  implicit none

  integer, save :: nx
  integer, save :: lo, hi
  integer, save :: nsteps
  integer, save :: bc(2)

  double precision, save :: dx
  double precision, save :: dt
  double precision, save :: prob_lo
  double precision, save :: prob_hi
  double precision, save :: v_in
  double precision, save :: flame_offset

contains

  subroutine read_probin(filename)
    character(len=*), intent(in   ) :: filename

    namelist /prbin/ nx, nsteps, dt, prob_hi, prob_lo, v_in, flame_offset

    ! defaults
    nx     = 256
    dt     = 1.d-10
    nsteps = 20
    bc     = [ 1, 2 ]

    prob_lo = 0.0d0
    prob_hi = 1.2d0

    v_in         = 1.2d0
    flame_offset = 0.d0

    ! read
    open(unit=66, file=filename, status='old', action='read')
    read(unit=66, nml=prbin)
    close(unit=66)

    lo = 0
    hi = nx-1
    dx = (prob_hi - prob_lo) / nx

  end subroutine read_probin

end module probin
