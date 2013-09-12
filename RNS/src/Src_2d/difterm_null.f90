module difterm_module
  use meth_params_module, only : NVAR
  implicit none
  private
  public :: difterm
contains

  subroutine difterm(lo,hi,U,Ulo,Uhi,fx,fy,dxinv)
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(in   ) :: dxinv(2)
    double precision, intent(in   ) ::  U(Ulo(1):Uhi(1)  ,Ulo(2):Uhi(2)  ,NVAR)
    double precision, intent(inout) :: fx( lo(1): hi(1)+1, lo(2): hi(2)  ,NVAR)
    double precision, intent(inout) :: fy( lo(1): hi(1)  , lo(2): hi(2)+1,NVAR)
    return
  end subroutine difterm

end module difterm_module
