! RHS routines for PFASST

module feval
  use iso_c_binding
  use layout_module
  use multifab_module
  implicit none
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(y0, t, nvars, yex)
    integer,      intent(in)  :: nvars
    real(kind=8), intent(in)  :: y0(nvars), t
    real(kind=8), intent(out) :: yex(nvars)

  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine copy_to_mfab(mf, y)
    type(multifab),   intent(inout) :: mf
    double precision, intent(in   ) :: y(:)

    integer :: i, j, k, m, n

    double precision, pointer :: dp(:,:,:,:)
    integer :: lo(3), hi(3)

    dp => dataptr(mf,1)
    lo = lwb(get_box(mf,1))
    hi = upb(get_box(mf,1))

    dp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = reshape(y, &
         [ hi(1)-lo(1)+1, hi(2)-lo(2)+1, hi(3)-lo(3)+1, size(dp, 4) ] )

  end subroutine copy_to_mfab

  subroutine copy_from_mfab(y, mf)
    type(multifab),   intent(inout) :: mf
    double precision, intent(out  ) :: y(:)

    integer :: i, j, k, m, n

    double precision, pointer :: dp(:,:,:,:)
    integer :: lo(3), hi(3)

    dp => dataptr(mf,1)
    lo = lwb(get_box(mf,1))
    hi = upb(get_box(mf,1))

    y = reshape(dp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:), [ size(y) ])
  end subroutine copy_from_mfab


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

    if (nboxes(U) > 1) then
       stop "Error: multibox PFASST not supported yet."
    end if

    call copy_to_mfab(U, y)
    call dUdt(U, F, ctx)
    call copy_from_mfab(f1, F)

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
