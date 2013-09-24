
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

  type(bdf_ts)  :: ts
  type(bdf_ctx) :: ctx
  double precision :: rtol(neq), atol(neq), dt
  ! double precision :: y0(neq), t0, y1(neq), t1

  rtol = 1
  atol = 1

  call bdf_ts_build(ts, neq, f, J, rtol, atol, max_order=6)

  ts%k = 5
  dt   = 0.021388d0
  ts%t = [ 5*dt, 4*dt, 3*dt, 2*dt, dt ]

  ! call nordsieck_update_coeffs(ts)

  call assert(abs(alpha0(ts%k) - (-2.2833333333333332d0)) < 1d-12, "error in alpha0")
  call assert(abs(alphahat0(ts%k,ts%t) - (-2.2833333333333332d0)) < 1d-12, "error in alphahat0")

contains

  subroutine assert(cond, message)
    logical, intent(in) :: cond
    character(*), intent(in) :: message

    if (.not. cond) then
       print *, "ASSERTION FAILED: ", message
       call abort
    end if
  end subroutine assert

end program test
