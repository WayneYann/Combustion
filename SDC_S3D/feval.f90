! RHS routines for PFASST

module feval
  use iso_c_binding
  use layout_module
  use multifab_module
  use encap
  implicit none
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(y0, t, nvars, yex)
    integer,      intent(in)  :: nvars
    real(kind=8), intent(in)  :: y0(nvars), t
    real(kind=8), intent(out) :: yex(nvars)

  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! subroutine copy_to_mfab(mf, y)
  !   type(multifab),   intent(inout) :: mf
  !   double precision, intent(in   ) :: y(:)

  !   double precision, pointer :: dp(:,:,:,:)
  !   integer :: lo(3), hi(3)

  !   dp => dataptr(mf,1)
  !   lo = lwb(get_box(mf,1))
  !   hi = upb(get_box(mf,1))

  !   dp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = reshape(y, &
  !        [ hi(1)-lo(1)+1, hi(2)-lo(2)+1, hi(3)-lo(3)+1, size(dp, 4) ] )

  ! end subroutine copy_to_mfab

  subroutine copy_from_mfab(y, mf)
    type(multifab),   intent(inout) :: mf
    double precision, intent(out  ) :: y(:)

    double precision, pointer :: dp(:,:,:,:)
    integer :: lo(3), hi(3)

    dp => dataptr(mf,1)
    lo = lwb(get_box(mf,1))
    hi = upb(get_box(mf,1))

    y = reshape(dp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:), [ size(y) ])
  end subroutine copy_from_mfab


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit function at y, t.
  subroutine eval_f1(y, t, level, ctxp, f1)
    use advance_module
    use multifab_module

    integer,          intent(in)    :: level
    real(kind=8),     intent(in)    :: t
    type(pf_encap_t), intent(inout) :: y
    type(pf_encap_t), intent(inout) :: f1
    type(c_ptr),      intent(in)    :: ctxp

    type(feval_ctx_t), pointer :: ctx

    call c_f_pointer(ctxp, ctx)

    call dUdt(y%q, f1%q, ctx)
  end subroutine eval_f1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine eval_f2(y, t, level, ctxp, f2)
    integer,          intent(in)    :: level
    real(kind=8),     intent(in)    :: t
    type(pf_encap_t), intent(inout) :: y
    type(pf_encap_t), intent(inout) :: f2
    type(c_ptr),      intent(in)    :: ctxp

    call setval(f2%q, 0.0d0)
  end subroutine eval_f2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine comp_f2(y, t, dt, rhs, level, ctxp, f2)
    integer,          intent(in)    :: level
    real(kind=8),     intent(in)    :: t, dt
    type(pf_encap_t), intent(inout) :: y, rhs
    type(pf_encap_t), intent(inout) :: f2
    type(c_ptr),      intent(in)    :: ctxp

    call setval(f2%q, 0.0d0)
    call copy(y%q, rhs%q)
  end subroutine comp_f2

end module feval
