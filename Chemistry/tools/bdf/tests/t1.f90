module feval
  use bdf
  implicit none
  integer, parameter :: neq = 3
contains
  ! from example1.f90 in the f90 version of vode
  subroutine f(neq, y, t, ydot)
    integer,       intent(in)  :: neq
    real(dp),      intent(in)  :: y(neq), t
    real(dp),      intent(out) :: ydot(neq)
    YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
    YDOT(3) = 3.E7*Y(2)*Y(2)
    YDOT(2) = -YDOT(1) - YDOT(3)
  end subroutine f
  subroutine J(neq, y, t, pd)
    integer,       intent(in)  :: neq
    real(dp),      intent(in)  :: y(neq), t
    real(dp),      intent(out) :: pd(neq,neq)
    PD(1,1) = -.04D0
    PD(1,2) = 1.D4*Y(3)
    PD(1,3) = 1.D4*Y(2)
    PD(2,1) = .04D0
    PD(2,3) = -PD(1,3)
    PD(3,2) = 6.E7*Y(2)
    PD(2,2) = -PD(1,2) - PD(3,2)
  end subroutine J
end module feval

program test
  use bdf
  use feval
  implicit none

  type(bdf_ts)  :: ts
  double precision :: rtol(neq), atol(neq), dt
  double precision :: y0(neq), t0, y1(neq), t1

  integer :: i, ierr

  y0 = [ 1.d0, 0.d0, 0.d0 ]
  t0 = 0.d0
  t1 = 0.4d0

  rtol = 1.d-4
  atol = [ 1.d-8, 1.d-14, 1.d-6 ]
  dt   = 1.d-8

  call bdf_ts_build(ts, neq, rtol, atol, max_order=3)
  ! ts%verbose = 1

  do i = 1, 11
     ! call bdf_advance(ts, f, J, neq, y0, t0, y1, t1, dt, i/=1, ierr)
     !call bdf_advance(ts, f, J, neq, y0, t0, y1, t1, dt, .true., i/=1, ierr)
     call bdf_advance(ts, f, J, neq, y0, t0, y1, t1, dt, .true., .true., ierr)
     print *, t1, y1
     y0 = y1
     t0 = t1
     t1 = 10*t1
     dt = 2*ts%dt
  end do

  print *, ts%n, ts%nfe, ts%nje, ts%nit, ts%nse


  call bdf_ts_destroy(ts)

end program test
