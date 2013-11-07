module probin
  implicit none

  integer, save :: nx
  integer, save :: lo, hi
  integer, save :: nsteps
  integer, save :: bc(2)
  integer, save :: probtype
  integer, save :: LeEQ1, coef_avg_harm

  double precision, save :: dx
  double precision, save :: dt
  double precision, save :: prob_lo
  double precision, save :: prob_hi
  double precision, save :: v_in
  double precision, save :: flame_offset
  double precision, save :: Pr, Sc, TMIN_TRANS


contains

  subroutine read_probin(filename)
    character(len=*), intent(in   ) :: filename

    namelist /prbin/ nx, nsteps, dt, prob_hi, prob_lo, v_in, flame_offset, probtype, LeEQ1, coef_avg_harm, Pr, Sc, TMIN_TRANS

    ! defaults
    nx     = 256
    dt     = 1.d-10
    nsteps = 20
    bc     = [ 1, 2 ]
    probtype = 1
    LeEQ1 = 0
    coef_avg_harm = 0

    prob_lo = 0.0d0
    prob_hi = 1.2d0

    Pr = 0.7d0
    Sc = 0.7d0
    TMIN_TRANS = 0.d0

    v_in         = 1.d20
    flame_offset = prob_hi/2

    ! read
    open(unit=66, file=filename, status='old', action='read')
    read(unit=66, nml=prbin)
    close(unit=66)

    lo = 0
    hi = nx-1
    dx = (prob_hi - prob_lo) / nx

  end subroutine read_probin

end module probin
