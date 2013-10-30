module feval
  use bdf
  implicit none
  integer, parameter :: neq = 3
  integer, parameter :: npt = 2
contains
  ! from example1.f90 in the f90 version of vode
  subroutine f(neq, npt, y, t, ydot)
    integer,  intent(in   ) :: neq, npt
    real(dp), intent(in   ) :: y(neq,npt), t
    real(dp), intent(  out) :: ydot(neq,npt)
    integer :: p
    do p = 1, npt
       ydot(1,p) = -.04d0*y(1,p) + 1.d4*y(2,p)*y(3,p)
       ydot(3,p) = 3.e7*y(2,p)*y(2,p)
       ydot(2,p) = -ydot(1,p) - ydot(3,p)
    end do
  end subroutine f
  subroutine J(neq, y, t, pd)
    integer,  intent(in   ) :: neq
    real(dp), intent(in   ) :: y(neq), t
    real(dp), intent(  out) :: pd(neq,neq)
    pd(1,1) = -.04d0
    pd(1,2) = 1.d4*y(3)
    pd(1,3) = 1.d4*y(2)
    pd(2,1) = .04d0
    pd(2,3) = -pd(1,3)
    pd(3,2) = 6.e7*y(2)
    pd(2,2) = -pd(1,2) - pd(3,2)
  end subroutine J
end module feval

program test
  use bdf
  use feval
  implicit none

  type(bdf_ts)  :: ts
  double precision :: rtol(neq), atol(neq), dt
  double precision :: y0(neq,npt), t0, y1(neq,npt), t1

  integer :: i, ierr

  y0(:,1) = [ 1.d0, 0.d0, 0.d0 ]
!  y0(:,2) = [ 1.d0, 0.d0, 0.d0 ]
  y0(:,2) = [ 0.98516927747181138d0, 3.3863452485889568d-5, 1.4796859075703273d-2 ]

  t0 = 0.d0
  t1 = 0.4d0

  rtol = 1.d-4
  atol = [ 1.d-8, 1.d-14, 1.d-6 ]
  dt   = 1.d-8

  call bdf_ts_build(ts, neq, npt, rtol, atol, max_order=3)
  ! ts%verbose = 1

  do i = 1, 11
     !call bdf_advance(ts, f, J, neq, y0, t0, y1, t1, dt, i/=1, ierr)
     !call bdf_advance(ts, f, J, neq, y0, t0, y1, t1, dt, .true., i/=1, ierr)
     call bdf_advance(ts, f, J, neq, npt, y0, t0, y1, t1, dt, .true., .false., ierr)
     print *, t1, ierr, y1(:,1)
     print *, t1, ierr, y1(:,2)
     y0 = y1
     t0 = t1
     t1 = 10*t1
     dt = 2*ts%dt
  end do

  print *, ts%n, ts%nfe, ts%nje, ts%nit, ts%nse

  call bdf_ts_destroy(ts)

end program test
