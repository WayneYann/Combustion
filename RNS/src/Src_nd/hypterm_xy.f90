module hypterm_xy_module

  implicit none

contains

  subroutine hypterm_xy(fac,lo,hi,Uz,Uzlo,Uzhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx,U0lo,U0hi,U0)

    use meth_params_module, only : NVAR
    use reconstruct_module, only : reconstruct
    use riemann_module, only : riemann

    integer, intent(in) :: lo(3), hi(3), Uzlo(3), Uzhi(3), fxlo(3), fxhi(3), &
         fylo(3), fyhi(3), U0lo(3), U0hi(3)
    double precision,intent(in   )::dx(3), fac
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,target,intent(in):: &
         &                          Uz(Uzlo(1):Uzhi(1),Uzlo(2):Uzhi(2),Uzlo(3):Uzhi(3),NVAR)
    double precision,target,optional,intent(in)::  &
         &                          U0(U0lo(1):U0hi(1),U0lo(2):U0hi(2),U0lo(3):U0hi(3),NVAR)

    integer :: i, j, k, n, g
    double precision, dimension(:,:,:), pointer :: UU
    double precision, dimension(:,:,:,:), pointer :: Ubase
    double precision, dimension(:,:,:), allocatable, target :: UG1y, UG2y
    double precision, dimension(:,:,:), allocatable :: ULy, URy, ULyG1, ULyG2, URyG1, URyG2
    double precision, dimension(:,:), allocatable :: U0x, U0y, UGyL, UGyR,flxx,flxy

    if (present(U0)) then
       Ubase => U0
    else
       Ubase => Uz
    end if

    allocate(UG1y (lo(1)-3:hi(1)+3,lo(2):hi(2)  ,NVAR))
    allocate(UG2y (lo(1)-3:hi(1)+3,lo(2):hi(2)  ,NVAR))
    allocate(ULy  (lo(1)-3:hi(1)+3,lo(2):hi(2)+1,NVAR))
    allocate(URy  (lo(1)-3:hi(1)+3,lo(2):hi(2)+1,NVAR))
    allocate(ULyG1(lo(1)  :hi(1)  ,lo(2):hi(2)+1,NVAR))
    allocate(ULyG2(lo(1)  :hi(1)  ,lo(2):hi(2)+1,NVAR))
    allocate(URyG1(lo(1)  :hi(1)  ,lo(2):hi(2)+1,NVAR))
    allocate(URyG2(lo(1)  :hi(1)  ,lo(2):hi(2)+1,NVAR))

    allocate(U0x (lo(1)-1:hi(1)+1,NVAR))
    allocate(flxx(lo(1)  :hi(1)+1,NVAR))

    allocate(U0y (lo(2)-1:hi(2)+1,NVAR))
    allocate(flxy(lo(2)  :hi(2)+1,NVAR))

    allocate(UGyL(lo(1):hi(1)+1,NVAR))
    allocate(UGyR(lo(1):hi(1)+1,NVAR))

    do k=lo(3),hi(3)
       ! Given cell averages, reconstruct in y-direction
       ! Note that they are still averges in x-direction
       do i=lo(1)-3,hi(1)+3
          U0y = Ubase(i,lo(2)-1:hi(2)+1,k,:)
          call reconstruct(lo(2),hi(2), & 
               Uzlo(2),Uzhi(2),   &  ! for Uz
               lo(2)  ,  hi(2)+1, &  ! for UL & UR
               lo(2)  ,  hi(2),   &  ! for UG1 & UG2
               lo(2)-1,  hi(2)+1, &  ! for U0
               Uz(i,:,k,:), &
               UL =ULy (i,:,:), UR =URy (i,:,:), &
               UG1=UG1y(i,:,:), UG2=UG2y(i,:,:), &
               U0 =U0y, &
               dir=2)
       end do

       ! flux in x-direction for two Gauss points in y-direction
       do g=1,2
          if (g.eq.1) then
             UU => UG1y
          else
             UU => UG2y
          end if

          do j=lo(2),hi(2)
             U0x = Ubase(lo(1)-1:hi(1)+1,j,k,:)

             call reconstruct(lo(1),hi(1),  &
                  lo(1)-3,hi(1)+3,  &  ! for UU
                  lo(1)  ,hi(1)+1,  &  ! for UL & UR
                  0,0,              &  ! UG1 & UG2 are not present
                  lo(1)-1,hi(1)+1,  &  ! for U0
                  UU(:,j,:), &
                  UL=UGyL, UR=UGyR, &
                  U0=U0x, &
                  dir=1)

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
       U0x = Ubase(lo(1)-1:hi(1)+1,lo(2)-1,k,:)
       !
       do j=lo(2), hi(2)+1
          call reconstruct(lo(1),hi(1), &
               lo(1)-3,hi(1)+3, &  ! for ULy
               0,0,             &  ! UL & UR are not present
               lo(1)  ,hi(1),   &  ! for UG1 & UG2
               lo(1)-1,hi(1)+1, &  ! for U0
               ULy(:,j,:), &
               UG1=ULyG1(:,j,:), UG2=ULyG2(:,j,:), &
               U0=U0x, &
               dir=1)

          U0x = Ubase(lo(1)-1:hi(1)+1,j,k,:)
          call reconstruct(lo(1),hi(1), &
               lo(1)-3,hi(1)+3, &  ! for URy
               0,0,             &  ! UL & UR are not present
               lo(1)  ,hi(1),   &  ! for UG1 & UG2
               lo(1)-1,hi(1)+1, &  ! for U0
               URy(:,j,:), &
               UG1=URyG1(:,j,:), UG2=URyG2(:,j,:), &
               U0=U0x, &
               dir=1)
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

    Nullify(Ubase)
    
  end subroutine hypterm_xy

end module hypterm_xy_module
