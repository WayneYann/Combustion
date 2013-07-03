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

    double precision, allocatable :: ULy(:,:,:), URy(:,:,:), UG1y(:,:,:), UG2y(:,:,:)
    double precision, allocatable :: UGyL(:,:), UGyR(:,:), flux(:,:)
    double precision, allocatable :: ULyG1(:,:,:), ULyG2(:,:,:)
    double precision, allocatable :: URyG1(:,:,:), URyG2(:,:,:)
    double precision, allocatable :: U0(:,:)
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
    allocate( ULy(tlo(1):thi(1),tlo(2):thi(2)+1,NVAR))
    allocate( URy(tlo(1):thi(1),tlo(2):thi(2)+1,NVAR))
    allocate(UG1y(tlo(1):thi(1),tlo(2):thi(2)  ,NVAR))
    allocate(UG2y(tlo(1):thi(1),tlo(2):thi(2)  ,NVAR))

    ! Given cell averages, reconstruct in y-direction
    ! Note that they are still averges in x-direction
    do i=tlo(1),thi(1)
       call reconstruct(tlo(2),thi(2), U(i,:,:), Ulo(2), Uhi(2), &
            UL=ULy(i,:,:), UR=URy(i,:,:), &
            UG1=UG1y(i,:,:), UG2=UG2y(i,:,:), &
            dir=2)
    end do

    allocate(UGyL(lo(1):hi(1)+1,NVAR))
    allocate(UGyR(lo(1):hi(1)+1,NVAR))
    allocate(flux(lo(1):hi(1)+1,NVAR))
    allocate(U0(tlo(1):thi(1),NVAR))

    ! flux in x-direction for two Gauss points in y-direction
    do j=lo(2),hi(2)
       U0 = U(tlo(1):thi(1),j,:)

       call reconstruct(lo(1),hi(1), UG1y(:,j,:), tlo(1), thi(1), &
            UL=UGyL, UR=UGyR,  &
            U0=U0,  &
            dir=1)

       call riemann(lo(1), hi(1), UGyL, UGyR, flux, dir=1)
       do n=1,NVAR
          do i=lo(1),hi(1)+1
             fx(i,j,n) = fx(i,j,n) + 0.5d0*flux(i,n)
          end do
       end do

       call reconstruct(lo(1),hi(1), UG2y(:,j,:), tlo(1), thi(1), &
            UL=UGyL, UR=UGyR,  &
            U0=U0,  &
            dir=1)

       call riemann(lo(1), hi(1), UGyL, UGyR, flux, dir=1)
       do n=1,NVAR
          do i=lo(1),hi(1)+1
             fx(i,j,n) = fx(i,j,n) + 0.5d0*flux(i,n)
          end do
       end do
    end do

    deallocate(UGyL, UGyR, flux, U0)
    deallocate(UG1y, UG2y)

    ! flux in y-direction
    allocate(flux(lo(2):hi(2)+1,NVAR))
    allocate(ULyG1(lo(1):hi(1),lo(2):hi(2)+1,NVAR))
    allocate(ULyG2(lo(1):hi(1),lo(2):hi(2)+1,NVAR))
    allocate(URyG1(lo(1):hi(1),lo(2):hi(2)+1,NVAR))
    allocate(URyG2(lo(1):hi(1),lo(2):hi(2)+1,NVAR))
    allocate(U0(lo(1)-3:hi(1)+3,NVAR))

    ! reconstruct in x-direction for states on y-faces
    do j=lo(2), hi(2)+1

       U0 = U(lo(1)-3:hi(1)+3,j-1,:)
       call reconstruct(lo(1),hi(1), ULy(:,j,:), lo(1)-3, hi(1)+3, &
            UG1=ULyG1(:,j,:), UG2=ULyG2(:,j,:),  &
            U0=U0,  &
            dir=1)

       U0 = U(lo(1)-3:hi(1)+3,j  ,:)
       call reconstruct(lo(1),hi(1), URy(:,j,:), lo(1)-3, hi(1)+3, &
            UG1=URyG1(:,j,:), UG2=URyG2(:,j,:),  &
            U0=U0,  &
            dir=1)
    end do

    ! flux in y-direction for two Gauss points in x-direction
    do i=lo(1), hi(1)

       call riemann(lo(2), hi(2), ULyG1(i,:,:), URyG1(i,:,:), flux, dir=2)
       do n=1,NVAR
          do j=lo(2),hi(2)+1
             fy(i,j,n) = fy(i,j,n) + 0.5d0*flux(j,n)
          end do
       end do

       call riemann(lo(2), hi(2), ULyG2(i,:,:), URyG2(i,:,:), flux, dir=2)
       do n=1,NVAR
          do j=lo(2),hi(2)+1
             fy(i,j,n) = fy(i,j,n) + 0.5d0*flux(j,n)
          end do
       end do
    end do

    deallocate(ULy,URy, URyG1, ULyG2, ULyG1, URyG2, flux, U0)

  end subroutine hypterm

end module hypterm_module

