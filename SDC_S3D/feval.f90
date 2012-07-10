! RHS routines for PFASST

module feval
  use iso_c_binding
  use layout_module
  implicit none
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(y0, t, nvars, yex)
    integer,      intent(in)  :: nvars
    real(kind=8), intent(in)  :: y0(nvars), t
    real(kind=8), intent(out) :: yex(nvars)

  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit function at y, t.
  subroutine eval_f1(y, t, nvars, level, ctxp, f1)
    use advance_module
    use multifab_module

    integer, intent(in)       :: nvars, level
    real(kind=8), intent(in)  :: y(nvars), t
    real(kind=8), intent(out) :: f1(nvars)
    type(c_ptr),  intent(in)  :: ctxp

    type(feval_ctx_t), pointer :: ctx

    type(multifab) :: U, F

    call c_f_pointer(ctxp, ctx)

    call build(U, ctx%la, ctx%nc, ctx%ng)
    call build(F, ctx%la, ctx%nc, 0)

    ! XXX: copy y into U

    call dUdt(U, F, ctx)

    ! XXX: extract f1 from F

    call destroy(U)
    call destroy(F)

  end subroutine eval_f1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine eval_f2(y, t, nvars, level, ctx, f2)
    integer, intent(in) :: nvars, level
    real(kind=8), intent(in ) :: y(nvars), t
    real(kind=8), intent(out) :: f2(nvars)
    type(c_ptr),  intent(in)  :: ctx

    !$OMP PARALLEL WORKSHARE
    f2 = 0.d0
    !$OMP END PARALLEL WORKSHARE
  end subroutine eval_f2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine comp_f2(y, t, dt, rhs, nvars, level, ctx, f2)
    integer, intent(in) :: nvars, level
    real(kind=8), intent(inout) :: y(nvars)
    real(kind=8), intent(in)    :: rhs(nvars), Dt, t
    real(kind=8), intent(out)   :: f2(nvars)
    type(c_ptr),  intent(in)    :: ctx

    !$OMP PARALLEL WORKSHARE
    f2 = 0.d0
    y  = rhs
    !$OMP END PARALLEL WORKSHARE
  end subroutine comp_f2

end module feval
