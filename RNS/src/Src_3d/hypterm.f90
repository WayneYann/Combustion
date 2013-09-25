module hypterm_module

  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, difmag
  use reconstruct_module, only : reconstruct
  use riemann_module, only : riemann

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fxlo(3), fxhi(3), &
         fylo(3), fyhi(3), fzlo(3), fzhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   ):: U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    integer :: tlo(3), thi(3), i, j, k, n, g
    double precision, dimension(:,:,:,:), allocatable :: Ulz, URz, UG1z, UG2z

    tlo(1) = lo(1)-3
    tlo(2) = lo(2)-3
    tlo(3) = lo(3)
    thi(1) = hi(1)+3
    thi(2) = hi(2)+3
    thi(3) = hi(3)

    allocate( ULz(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)+1,NVAR))
    allocate( URz(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)+1,NVAR))
    allocate(UG1z(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)  ,NVAR))
    allocate(UG2z(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)  ,NVAR))

    ! Given cell averages, reconstruct in z-direction
    ! Note that they are still averges in x and y-direction
    do j=tlo(2),thi(2)
    do i=tlo(1),thi(1)
       call reconstruct(tlo(3),thi(3), Ulo(3),Uhi(3), tlo(3),thi(3)+1, tlo(3),thi(3), 0,0,&
            U(i,j,:,:), UL=ULz(i,j,:,:), UR=URz(i,j,:,:),  &
            UG1=UG1z(i,j,:,:), UG2=UG2z(i,j,:,:), dir=3)
    end do
    end do

    call hypterm_xy(lo,hi,UG1z,tlo,thi,fx,fxlo,fxhi,fy,fylo,fyhi,dx, &
         U,Ulo,Uhi,0.5d0)
    call hypterm_xy(lo,hi,UG2z,tlo,thi,fx,fxlo,fxhi,fy,fylo,fyhi,dx, &
         U,Ulo,Uhi,0.5d0)

    deallocate(UG1z,UG2z)

    ! z-direction flux


    ! difmag

    deallocate(ULz,URz)

  end subroutine hypterm


  subroutine hypterm_xy(lo,hi,Uz,Uzlo,Uzhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx,U,Ulo,Uhi,fac)
    integer, intent(in) :: lo(3), hi(3), Uzlo(3), Uzhi(3), fxlo(3), fxhi(3), &
         fylo(3), fyhi(3), Ulo(3), Uhi(3)
    double precision,intent(in   )::dx(3), fac
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,intent(in   )::Uz(Uzlo(1):Uzhi(1),Uzlo(2):Uzhi(2),Uzlo(3):Uzhi(3),NVAR)
    double precision,intent(in   )::U ( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)

    integer :: i, j, k, n, g, tlo(2), thi(2)
    double precision, dimension(:,:,:), allocatable, target :: UG1y, UG2y
    double precision, dimension(:,:,:), pointer :: UU
    double precision, dimension(:,:,:), allocatable :: ULy, URy, ULyG1, ULyG2, URyG1, URyG2
    double precision, dimension(:,:), allocatable :: U0x, U0y, UGyL, UGyR,flxx,flxy

    tlo(1) = lo(1)-3
    tlo(2) = lo(2)
    thi(1) = hi(1)+3
    thi(2) = hi(2)
    allocate(UG1y (tlo(1):thi(1),tlo(2):thi(2)  ,NVAR))
    allocate(UG2y (tlo(1):thi(1),tlo(2):thi(2)  ,NVAR))
    allocate( ULy (tlo(1):thi(1),tlo(2):thi(2)+1,NVAR))
    allocate( URy (tlo(1):thi(1),tlo(2):thi(2)+1,NVAR))
    allocate(ULyG1( lo(1): hi(1), lo(2): hi(2)+1,NVAR))
    allocate(ULyG2( lo(1): hi(1), lo(2): hi(2)+1,NVAR))
    allocate(URyG1( lo(1): hi(1), lo(2): hi(2)+1,NVAR))
    allocate(URyG2( lo(1): hi(1), lo(2): hi(2)+1,NVAR))

    allocate( U0x(tlo(1):thi(1)  ,NVAR))
    allocate( U0y(tlo(2):thi(2)  ,NVAR))
    allocate(flxy( lo(2): hi(2)+1,NVAR))

    allocate(UGyL(lo(1):hi(1)+1,NVAR))
    allocate(UGyR(lo(1):hi(1)+1,NVAR))
    allocate(flxx(lo(1):hi(1)+1,NVAR))

    do k=lo(3),hi(3)
       ! Given cell averages, reconstruct in y-direction
       ! Note that they are still averges in x-direction
       do i=tlo(1),thi(1)
          U0y(tlo(2):thi(2),:) = U(i,tlo(2):thi(2),k,:)
          call reconstruct(tlo(2),thi(2), Uzlo(2),Uzhi(2), &
               tlo(2),thi(2)+1, tlo(2),thi(2), Ulo(2),Uhi(2), &
               Uz(i,:,k,:), UL=ULy(i,:,:), UR=URy(i,:,:), UG1=UG1y(i,:,:), UG2=UG2y(i,:,:), &
               U0=U0y, dir=2)
       end do

       ! flux in x-direction for two Gauss points in y-direction
       do g=1,2
          if (g.eq.1) then
             UU => UG1y
          else
             UU => UG2y
          end if

          do j=lo(2),hi(2)
             U0x = U(tlo(1):thi(1),j,k,:)

             call reconstruct(lo(1),hi(1), tlo(1),thi(1), lo(1),hi(1)+1, 0,0, tlo(1),thi(1), &
                  UU(:,j,:), UL=UGyL, UR=UGyR, U0=U0x, dir=1)

             call riemann(lo(1),hi(1), UGyL, UGyR, lo(1),hi(1)+1, flxx, lo(1),hi(1)+1, dir=1)
             do n=1,NVAR
                do i=lo(1),hi(1)+1
                   fx(i,j,k,n) = fx(i,j,k,n) + fac*0.5d0*flxx(i,n)
                end do
             end do
          end do
          
          Nullify(UU)
       end do

       ! y-direction flux

       ! reconstruct in x-direction for states on y-faces
       U0x = U(lo(1)-3:hi(1)+3,lo(2)-1,k,:)
       !
       do j=lo(2), hi(2)+1
          call reconstruct(lo(1),hi(1), lo(1)-3,hi(1)+3, 0,0, lo(1),hi(1), lo(1)-3,hi(1)+3, &
               ULy(:,j,:), UG1=ULyG1(:,j,:), UG2=ULyG2(:,j,:), U0=U0x, dir=1)

          U0x = U(lo(1)-3:hi(1)+3,j,k,:)
          call reconstruct(lo(1),hi(1), lo(1)-3,hi(1)+3, 0,0, lo(1),hi(1), lo(1)-3,hi(1)+3, &
               URy(:,j,:), UG1=URyG1(:,j,:), UG2=URyG2(:,j,:), U0=U0x, dir=1)
       end do

       ! flux in y-direction for two Gauss points in x-direction
       do i=lo(1), hi(1)

          call riemann(lo(2),hi(2), ULyG1(i,:,:), URyG1(i,:,:), lo(2),hi(2)+1, &
               flxy, lo(2),hi(2)+1, dir=2)

          do n=1,NVAR
             do j=lo(2),hi(2)+1
                fy(i,j,k,n) = fy(i,j,k,n) + fac*0.5d0*flxy(j,n)
             end do
          end do

          call riemann(lo(2),hi(2), ULyG2(i,:,:), URyG2(i,:,:), lo(2),hi(2)+1, &
               flxy, lo(2),hi(2)+1, dir=2)

          do n=1,NVAR
             do j=lo(2),hi(2)+1
                fy(i,j,k,n) = fy(i,j,k,n) + fac*0.5d0*flxy(j,n)
             end do
          end do
       end do
    end do
    
    deallocate(UG1y,UG2y)
    deallocate(ULy,URy,ULyG1,ULyG2,URyG1,URyG2)
    deallocate(U0x,U0y,UGyL,UGyR,flxx,flxy)
    
  end subroutine hypterm_xy

end module hypterm_module
