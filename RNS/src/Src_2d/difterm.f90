module difterm_module

  use meth_params_module, only : NVAR, NSPEC, QCVAR, QFVAR
  use weno_module, only : cellavg2gausspt_1d, cellavg2face_1d

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

    double precision, dimension(:,:,:), pointer :: Uag => Null()
    double precision, allocatable, target :: U1(:,:,:), U2(:,:,:)
    double precision, allocatable :: Qc(:,:,:), Qf(:,:,:)
    double precision, allocatable :: mu(:), xi(:), lam(:), Ddia(:,:)
    integer :: i, j, n, g

    ! ----- compute x-direction flux first -----

    allocate(U1(lo(1)-2:hi(1)+2,lo(2):hi(2),NVAR))
    allocate(U2(lo(1)-2:hi(1)+2,lo(2):hi(2),NVAR))

    allocate(Qf(lo(1)-1:hi(1)+1,lo(2):hi(2),QFVAR))


    ! cell-average => cell-avg-in-x and Gauss-point-in-y
    do n=1,NVAR
       do i=lo(1)-2,hi(1)+2
!          call cellavg2gausspt_1d(U(i,lo(2)-2:hi(2)+2,n), lo(2)-2, hi(2)+2, &
!               U1, U2, lo(2), hi(2))
       end do
    end do

    do g=1,2
       
       if (g .eq. 1) then
          Uag => U1
       else
          Uag => U2
       end if

       ! cell-avg-in-x and Gauss-point-in-y => xface and Gauss-point-in-y
       do n=1,NVAR
          do j=lo(2),hi(2)
!             call cellavg2face_1d(U1(:,j,n), lo(1)-2, hi(1)+2, &
!                  Qf1, lo(1), hi(1))
          end do
       end do

       Nullify(Uag)
    end do

    ! ----- compute y-direction flux -----

  end subroutine difterm

end module difterm_module
