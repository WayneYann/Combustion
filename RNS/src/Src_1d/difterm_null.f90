module difterm_module
  implicit none
  private
  public :: difterm
contains

  subroutine difterm(lo,hi,domlo,domhi,U,Ulo,Uhi,flx, dx)
    use meth_params_module, only : NVAR
    integer, intent(in) :: lo(1), hi(1), domlo(1), domhi(1), Ulo(1), Uhi(1)
    double precision, intent(in ) ::   U(Ulo(1):Uhi(1)  ,NVAR)
    double precision, intent(out) :: flx( lo(1): hi(1)+1,NVAR)
    double precision, intent(in) :: dx(1)
    flx = 0.d0
  end subroutine difterm
end module difterm_module
