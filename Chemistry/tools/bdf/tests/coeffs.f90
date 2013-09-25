
module bdf_params
  type :: bdf_ctx
  end type bdf_ctx
end module bdf_params

module feval
  use bdf
  implicit none
  integer, parameter :: neq = 1
contains
  subroutine f(neq, y, t, ydot, ctx)
    integer,          intent(in)  :: neq
    real(dp),         intent(in)  :: y(neq), t
    real(dp),         intent(out) :: ydot(neq)
    type(bdf_ctx),    intent(in)  :: ctx
  end subroutine f
  subroutine J(neq, y, t, pd, ctx)
    integer,          intent(in)  :: neq
    real(dp),         intent(in)  :: y(neq), t
    real(dp),         intent(out) :: pd(neq,neq)
    type(bdf_ctx),    intent(in)  :: ctx
  end subroutine J
end module feval
program test
  use bdf
  use feval

  double precision, parameter :: tol = 1.d-12

  type(bdf_ts)  :: ts
  type(bdf_ctx) :: ctx
  double precision :: rtol(neq), atol(neq)
  double precision :: a, v(0:6)

  rtol = 1
  atol = 1

  call bdf_ts_build(ts, neq, f, J, rtol, atol, max_order=6)

  print *, "====> order 1"
  ts%k = 1
  call random_number(ts%h)
  call bdf_ts_update(ts)
  print *, 'l', ts%l(0:1)
  call assert(all(ts%l(0:1) == [ 1.d0, 1.d0 ]), "error in l")

  ts%h = 1
  call bdf_ts_update(ts)
  print *, 'tq', ts%tq
  call assert(all(abs(ts%tq - [ 1.d0, 0.5d0, 0.2222222222222d0, 2.0d0 ]) < tol), "error in tq")

  print *, "====> order 2"
  ts%k = 2
  call random_number(ts%h)
  call bdf_ts_update(ts)
  print *, 'l', ts%l(0:2)
  call assert(all(abs(ts%l(0:2) - [ 1.d0, 1.5d0, 0.5d0 ]) < tol), "error in l")

  ts%h = 1
  call bdf_ts_update(ts)
  print *, 'tq', ts%tq
  call assert(all(abs(ts%tq - [ 1.0d0, 0.222222222222222d0, 0.13636363636363638d0, 6.0d0 ]) < tol), "error in tq")

  print *, "====> order 3"
  ts%k = 3
  ts%h = 1
  call bdf_ts_update(ts)
  print *, 'l', ts%l(0:3)
  v(0:3) = [ 1.d0, 1.8333333333333d0, 1.d0, 0.1666666666666d0 ]
  call assert(all(abs(ts%l(0:3) - v(0:3)) < tol), "error in l")

  ts%h(1) = 4
  call bdf_ts_update(ts)
  print *, 'l', ts%l(0:3)
  v(0:3) = [ 1.d0, 1.8333333333333d0, 0.96d0, 0.12666666666666d0 ]
  call assert(all(abs(ts%l(0:3) - v(0:3)) < tol), "error in l")

  ts%h = 1
  call bdf_ts_update(ts)
  print *, 'tq', ts%tq
  call assert(all(abs(ts%tq - [ 1.3333333333333d0, 0.13636363636363638d0, 9.6d-2, 24.0d0 ]) < tol), "error in tq")
  
  print *, "====> order 5"
  ts%k = 5
  ts%h = 1
  a = -2.2833333333333332d0
  print *, 'alpha0, alphahat0', a
  call assert(abs(alpha0(ts%k)         - a) < tol, "error in alpha0")
  call assert(abs(alphahat0(ts%k,ts%h) - a) < tol, "error in alphahat0")

  call bdf_ts_destroy(ts)

contains

  subroutine assert(cond, message)
    logical, intent(in) :: cond
    character(*), intent(in) :: message

    if (.not. cond) then
       print *, "ASSERTION FAILED: ", message
       call abort()
    end if
  end subroutine assert

end program test
