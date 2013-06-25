module difterm_module

  use meth_params_module, only : NVAR

  implicit none

  private

  public :: difterm

contains

  subroutine difterm(lo,hi,U,Ulo,Uhi,flx)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(in ) ::   U(Ulo(1):Uhi(1)  ,NVAR)
    double precision, intent(out) :: flx( lo(1): hi(1)+1,NVAR)

    flx = 0.d0

  end subroutine difterm

end module difterm_module
