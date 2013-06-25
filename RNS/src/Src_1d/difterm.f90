module difterm_module

  use meth_params_module, only : NVAR

  implicit none

  private

  public :: difterm

contains

  subroutine difterm(lo,hi,U,Ulo,Uhi,flx)
    use weno_module, only : cellavg2face_1d
    use convert_module,only : cellavg2cc_1d
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(in ) ::   U(Ulo(1):Uhi(1)  ,NVAR)
    double precision, intent(out) :: flx( lo(1): hi(1)+1,NVAR)

    double precision, allocatable :: Uc(:,:), Uf(:,:)  ! cell-center and face
    integer :: n

    allocate(Uc(lo(1)-2:hi(1)+2, NVAR))
    allocate(Uf(lo(1)  :hi(1)+1, NVAR))

    do n=1,NVAR
       call cellavg2cc_1d(U, Ulo(1), Uhi(1), Uc, lo(1)-2, hi(1)+2)
    end do

    do n=1,NVAR
       call cellavg2face_1d(U, Ulo(1), Uhi(1), Uf, lo(1), hi(1)+1)
    end do

    flx = 0.d0

    deallocate(Uc,Uf)

  end subroutine difterm

end module difterm_module
