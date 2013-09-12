module difterm_module
  use meth_params_module, only : NVAR
  implicit none
  private
  public :: difterm
contains

  subroutine difterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dxinv)
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    double precision, intent(in   ) :: dxinv(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),NVAR)
    return
  end subroutine difterm

end module difterm_module
