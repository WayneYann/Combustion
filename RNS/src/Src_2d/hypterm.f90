module hypterm_module

  use meth_params_module, only : NVAR
  use reconstruct_module, only : reconstruct
  use riemann_module, only : riemann

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,fx,fy)
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(in ) ::  U(Ulo(1):Uhi(1)  ,Ulo(2):Uhi(2)  ,NVAR)
    double precision, intent(out) :: fx( lo(1): hi(1)+1, lo(2): hi(2)  ,NVAR)
    double precision, intent(out) :: fy( lo(1): hi(1)  , lo(2): hi(2)+1,NVAR)

    double precision, allocatable :: UL_y(:,:,:), UR_y(:,:,:), UG1_y(:,:,:), UG2_y(:,:,:)
    double precision, allocatable :: UL(:,:), UR(:,:), flux(:,:)
    integer :: tlo(2), thi(2), i, j, n

    do n=1,NVAR
       do j=lo(2), hi(2)
          do i=lo(1), hi(1)+1
             fx(i,j,n) = 0.d0
          end do
       end do

       do j=lo(2), hi(2)+1
          do i=lo(1), hi(1)
             fy(i,j,n) = 0.d0
          end do
       end do
    end do

    tlo(1) = lo(1)-3
    tlo(2) = lo(2)
    thi(1) = hi(1)+3
    thi(2) = hi(2)
    allocate( UL_y(tlo(1):thi(1),tlo(2):thi(2)+1,NVAR))
    allocate( UR_y(tlo(1):thi(1),tlo(2):thi(2)+1,NVAR))
    allocate(UG1_y(tlo(1):thi(1),tlo(2):thi(2)  ,NVAR))
    allocate(UG2_y(tlo(1):thi(1),tlo(2):thi(2)  ,NVAR))

    ! Given cell averages, reconstruct in y-direction
    ! Note that they are still averges in x-direction
    do i=tlo(1),thi(1)
       call reconstruct(tlo(2),thi(2), U(i,:,:), Ulo(2), Uhi(2), &
            UL=UL_y(i,:,:), UR=UR_y(i,:,:), UG1=UG1_y(i,:,:), UG2=UG2_y(i,:,:), dir=2)
    end do

    allocate(UL  (lo(1):hi(1)+1,NVAR))
    allocate(UR  (lo(1):hi(1)+1,NVAR))
    allocate(flux(lo(1):hi(1)+1,NVAR))

    ! flux in x-direction for two Gauss points in y-direction
    do j=lo(2),hi(2)
       call reconstruct(lo(1),hi(1), UG1_y(:,j,:), tlo(1), thi(1), &
            UL=UL, UR=UR, dir=1)
       call riemann(lo(1), hi(1), UL, UR, flux, dir=1)
       do n=1,NVAR
          do i=lo(1),hi(1)+1
             fx(i,j,n) = fx(i,j,n) + 0.5d0*flux(i,n)
          end do
       end do

       call reconstruct(lo(1),hi(1), UG2_y(:,j,:), tlo(1), thi(1), &
            UL=UL, UR=UR, dir=1)
       call riemann(lo(1), hi(1), UL, UR, flux, dir=1)
       do n=1,NVAR
          do i=lo(1),hi(1)+1
             fx(i,j,n) = fx(i,j,n) + 0.5d0*flux(i,n)
          end do
       end do
    end do

    deallocate(UL,UR, flux)
    
    ! flux in y-direction

    deallocate(UL_y,UR_y, UG1_y, UG2_y)

  end subroutine hypterm

end module hypterm_module

