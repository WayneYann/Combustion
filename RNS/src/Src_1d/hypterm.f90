module hypterm_module

  use meth_params_module, only : NVAR
  use weno_module, only : reconstruct
  use riemann_module, only : riemann

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,flx)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(in ) ::   U(Ulo(1):Uhi(1)  ,NVAR)
    double precision, intent(out) :: flx( lo(1): hi(1)+1,NVAR)

    double precision, allocatable :: UL(:,:), UR(:,:)

    allocate(UL(lo(1):hi(1)+1,NVAR))
    allocate(UR(lo(1):hi(1)+1,NVAR))

    call reconstruct(lo, hi, U, Ulo, Uhi, UL, UR)
    
    call riemann(lo, hi, UL, UR, flx)
    
    deallocate(UL,UR)

  end subroutine hypterm

end module hypterm_module

