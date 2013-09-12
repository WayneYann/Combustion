module hypterm_module

  use meth_params_module, only : NVAR, URHO, UMX, UMY, UTEMP, difmag
  use reconstruct_module, only : reconstruct
  use riemann_module, only : riemann

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,fx,fy,dx)
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(in   ) :: dx(2)
    double precision, intent(in   ) ::  U(Ulo(1):Uhi(1)  ,Ulo(2):Uhi(2)  ,NVAR)
    double precision, intent(inout) :: fx( lo(1): hi(1)+1, lo(2): hi(2)  ,NVAR)
    double precision, intent(inout) :: fy( lo(1): hi(1)  , lo(2): hi(2)+1,NVAR)

    double precision, dimension(:,:,:), pointer :: UU
    double precision, dimension(:,:,:), allocatable, target :: UG1y, UG2y
    double precision, allocatable :: ULy(:,:,:), URy(:,:,:)
    double precision, allocatable :: UGyL(:,:), UGyR(:,:), flux(:,:)
    double precision, allocatable :: ULyG1(:,:,:), ULyG2(:,:,:)
    double precision, allocatable :: URyG1(:,:,:), URyG2(:,:,:)
    double precision, allocatable :: U0(:,:)
    integer :: tlo(2), thi(2), i, j, n, g

    double precision, allocatable :: divv(:,:), vx(:,:), vy(:,:)
    double precision :: rhoInv, dvdx, dvdy, div1, dxinv(2)

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
       call reconstruct(tlo(2),thi(2), Ulo(2),Uhi(2), tlo(2),thi(2)+1, tlo(2),thi(2), 0,0,&
            U(i,:,:), UL=ULy(i,:,:), UR=URy(i,:,:), UG1=UG1y(i,:,:), UG2=UG2y(i,:,:), dir=2)
    end do

    allocate(UGyL(lo(1):hi(1)+1,NVAR))
    allocate(UGyR(lo(1):hi(1)+1,NVAR))
    allocate(flux(lo(1):hi(1)+1,NVAR))
    allocate(U0(tlo(1):thi(1),NVAR))

    ! flux in x-direction for two Gauss points in y-direction
    do g=1,2

       if (g.eq.1) then
          UU => UG1y
       else
          UU => UG2y
       end if

       do j=lo(2),hi(2)
          U0 = U(tlo(1):thi(1),j,:)

          call reconstruct(lo(1),hi(1), tlo(1),thi(1), lo(1),hi(1)+1, 0,0, tlo(1),thi(1), &
               UU(:,j,:), UL=UGyL, UR=UGyR, U0=U0, dir=1)

          call riemann(lo(1), hi(1), UGyL, UGyR, flux, dir=1)
          do n=1,NVAR
             do i=lo(1),hi(1)+1
                fx(i,j,n) = fx(i,j,n) + 0.5d0*flux(i,n)
             end do
          end do
       end do

       Nullify(UU)

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
    U0 = U(lo(1)-3:hi(1)+3,lo(2)-1,:)
    !
    do j=lo(2), hi(2)+1

       call reconstruct(lo(1),hi(1), lo(1)-3,hi(1)+3, 0,0, lo(1),hi(1), lo(1)-3,hi(1)+3, &
            ULy(:,j,:), UG1=ULyG1(:,j,:), UG2=ULyG2(:,j,:), U0=U0, dir=1)

       U0 = U(lo(1)-3:hi(1)+3,j  ,:)
       call reconstruct(lo(1),hi(1), lo(1)-3,hi(1)+3, 0,0, lo(1),hi(1), lo(1)-3,hi(1)+3, &
            URy(:,j,:), UG1=URyG1(:,j,:), UG2=URyG2(:,j,:), U0=U0, dir=1)
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

    if (difmag .gt. 0.0d0) then

       dxinv = 1.d0/dx

       allocate(vx  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
       allocate(vy  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
       allocate(divv(lo(1)  :hi(1)+1,lo(2)  :hi(2)+1))

       do    j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             rhoInv = 1.d0/U(i,j,URHO)
             vx(i,j) = U(i,j,UMX)*rhoInv
             vy(i,j) = U(i,j,UMY)*rhoInv
          end do
       end do

       do    j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1
             
             dvdx = 0.5d0*(vx(i,j)-vx(i-1,j)+vx(i,j-1)-vx(i-1,j-1))*dxinv(1)
             dvdy = 0.5d0*(vy(i,j)-vy(i,j-1)+vy(i-1,j)-vy(i-1,j-1))*dxinv(2)
             divv(i,j) = dvdx + dvdy

          end do
       end do

       do n=1,NVAR
          if (n.ne.UTEMP) then

             do   j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  div1 = .5d0*(divv(i,j) + divv(i,j+1))
                  div1 = difmag*min(0.d0,div1)
                  fx(i,j,n) = fx(i,j,n) + dx(1)*div1*(U(i,j,n) - U(i-1,j,n))
               enddo
            enddo

            do    j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  div1 = .5d0*(divv(i,j) + divv(i+1,j))
                  div1 = difmag*min(0.d0,div1)
                  fy(i,j,n) = fy(i,j,n) + dx(2)*div1*(U(i,j,n) - U(i,j-1,n))
               enddo
            enddo

          end if
       end do

       deallocate(vx,vy,divv)
    end if

  end subroutine hypterm

end module hypterm_module

