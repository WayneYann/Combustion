! Transfer (interpolate, restrict) routines for PFASST.

module transfer
  use iso_c_binding
  use encap
  implicit none
contains

  !
  ! Interpolate coarse qG to fine qF.
  !
  subroutine interpolate(qF, qG, levelF, ctxF, levelG, ctxG)
    use feval
    integer,          intent(in)    :: levelF, levelG
    type(pf_encap_t), intent(inout) :: qF, qG
    type(c_ptr),      intent(in)    :: ctxF, ctxG

    call copy(qF, qG)
  end subroutine interpolate


  !
  ! Restrict fine qF to coarse qG
  !
  subroutine restrict(qF, qG, levelF, ctxF, levelG, ctxG)
    integer,          intent(in)    :: levelF, levelG
    type(pf_encap_t), intent(inout) :: qF, qG
    type(c_ptr),      intent(in)    :: ctxF, ctxG

    call copy(qG, qF)
  end subroutine restrict

end module transfer
