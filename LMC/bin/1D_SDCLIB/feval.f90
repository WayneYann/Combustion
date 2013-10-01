module feval
  use iso_c_binding
  ! use sdc_mod_farray
  implicit none
contains

  subroutine f1eval(fptr, qptr, t, state, ctx) bind(c)
    type(c_ptr),    intent(in), value :: fptr, qptr, state, ctx
    real(c_double), intent(in), value :: t

  end subroutine f1eval

  subroutine f2eval(fptr, qptr, t, state, ctx) bind(c)
    type(c_ptr),    intent(in), value :: fptr, qptr, state, ctx
    real(c_double), intent(in), value :: t

  end subroutine f2eval

  subroutine f2comp(fptr, qptr, t, dt, rhsptr, state, ctx) bind(c)
    type(c_ptr),    intent(in), value :: fptr, qptr, rhsptr, state, ctx
    real(c_double), intent(in), value :: t, dt

  end subroutine f2comp

end module feval
