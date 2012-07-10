! Transfer (interpolate, restrict) routines for PFASST.

module transfer
  use iso_c_binding
  implicit none
contains

  !
  ! Interpolate coarse qG to fine qF.
  !
  subroutine interpolate(qF, qG, nvarF, nvarG, levelF, ctxF, levelG, ctxG)
    use feval
    integer,      intent(in)    :: nvarF, nvarG, levelF, levelG
    real(kind=8), intent(inout) :: qF(nvarF)
    real(kind=8), intent(in)    :: qG(nvarG)
    type(c_ptr),  intent(in)    :: ctxF, ctxG

    !$OMP PARALLEL WORKSHARE
    qF = qG
    !$OMP END PARALLEL WORKSHARE
  end subroutine interpolate


  !
  ! Restrict fine qF to coarse qG
  !
  subroutine restrict(qF, qG, nvarF, nvarG, levelF, ctxF, levelG, ctxG)
    integer,      intent(in)    :: nvarF, nvarG, levelF, levelG
    real(kind=8), intent(in)    :: qF(nvarF)
    real(kind=8), intent(inout) :: qG(nvarG)
    type(c_ptr),  intent(in)    :: ctxF, ctxG

    !$OMP PARALLEL WORKSHARE
    qG = qF
    !$OMP END PARALLEL WORKSHARE
  end subroutine restrict

end module transfer
