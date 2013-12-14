module hypterm_module

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)

    use meth_params_module, only : NVAR, difmag
    use reconstruct_module, only : reconstruct_3d
    use hypterm_xy_module, only : hypterm_xy

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fxlo(3), fxhi(3), &
         fylo(3), fyhi(3), fzlo(3), fzhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   ):: U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    integer :: Glo(3), Ghi(3), Flo(3), Fhi(3), tlo(3), thi(3)
    double precision, dimension(:,:,:,:), allocatable :: Ulz, URz, UG1z, UG2z

    Glo(1:2) = lo(1:2)-3
    Glo(3)   = lo(3)
    Ghi(1:2) = hi(1:2)+3
    Ghi(3)   = hi(3)

    Flo(1:2) = lo(1:2)-3
    Flo(3)   = lo(3)
    Fhi(1:2) = hi(1:2)+3
    Fhi(3)   = hi(3)  +1    

    allocate( ULz(Flo(1):Fhi(1),Flo(2):Fhi(2),Flo(3):Fhi(3),NVAR))
    allocate( URz(Flo(1):Fhi(1),Flo(2):Fhi(2),Flo(3):Fhi(3),NVAR))
    allocate(UG1z(Glo(1):Ghi(1),Glo(2):Ghi(2),Glo(3):Ghi(3),NVAR))
    allocate(UG2z(Glo(1):Ghi(1),Glo(2):Ghi(2),Glo(3):Ghi(3),NVAR))

    ! Given cell averages, reconstruct in z-direction
    ! Note that they are still averges in x and y-direction
    call reconstruct_3d(3, Glo, Ghi, &
         Ulo, Uhi,   &  ! for U
         Flo, Fhi,   &  ! for UL & UR
         Glo, Ghi,   &  ! for UG1 & UG2
         Ulo, Ulo,   &  ! U0 is not present
         U,          &
         UL=ULz, UR=URz,  &
         UG1=UG1z, UG2=UG2z)

    call hypterm_xy(0.5d0,lo,hi,UG1z,Glo,Ghi,fx,fxlo,fxhi,fy,fylo,fyhi,dx, &
         Ulo,Uhi,U)
    call hypterm_xy(0.5d0,lo,hi,UG2z,Glo,Ghi,fx,fxlo,fxhi,fy,fylo,fyhi,dx, &
         Ulo,Uhi,U)

    deallocate(UG1z,UG2z)

    call hypterm_z(lo,hi,U,Ulo,Uhi,ULz,URz,Flo,Fhi,fz,fzlo,fzhi,dx)

    deallocate(ULz,URz)

    if (difmag .gt. 0.0d0) then
       call add_artifical_viscocity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)
    end if

  end subroutine hypterm


  subroutine hypterm_z(lo,hi,U,Ulo,Uhi,UL,UR,zlo,zhi,fz,fzlo,fzhi,dx)

    use meth_params_module, only : NVAR
    use reconstruct_module, only : reconstruct_3d
    use riemann_module, only : riemann

    integer, intent(in) :: lo(3),hi(3),Ulo(3),Uhi(3),zlo(3),zhi(3),fzlo(3),fzhi(3)
    double precision,intent(in) :: dx(3)
    double precision,intent(in   )::U ( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(in   )::UL( zlo(1): zhi(1), zlo(2): zhi(2), zlo(3): zhi(3),NVAR)
    double precision,intent(in   )::UR( zlo(1): zhi(1), zlo(2): zhi(2), zlo(3): zhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    integer :: i,j,k,n, Flo(3), Fhi(3), Glo(3), Ghi(3), U0lo(3), U0hi(3)
    double precision, dimension(:,:,:,:), allocatable :: UL1,UL2,UL3,UL4,UR1,UR2,UR3,UR4
    double precision, dimension(:,:,:,:), allocatable :: UG1, UG2
    double precision, dimension(:,:), allocatable :: flx

    Flo = lo
    Fhi(1:2) = hi(1:2)
    Fhi(3)   = hi(3) + 1

    Glo(1) = lo(1) - 2
    Glo(2) = lo(2)
    Glo(3) = lo(3)
    Ghi(1) = hi(1) + 2
    Ghi(2) = hi(2)
    Ghi(3) = hi(3) + 1

    allocate(UG1(Glo(1):Ghi(1),Glo(2):Ghi(2),Glo(3):Ghi(3),NVAR))
    allocate(UG2(Glo(1):Ghi(1),Glo(2):Ghi(2),Glo(3):Ghi(3),NVAR))

    allocate(UL1(Flo(1):Fhi(1),Flo(2):Fhi(2),Flo(3):Fhi(3),NVAR))
    allocate(UL2(Flo(1):Fhi(1),Flo(2):Fhi(2),Flo(3):Fhi(3),NVAR))
    allocate(UL3(Flo(1):Fhi(1),Flo(2):Fhi(2),Flo(3):Fhi(3),NVAR))
    allocate(UL4(Flo(1):Fhi(1),Flo(2):Fhi(2),Flo(3):Fhi(3),NVAR))

    allocate(UR1(Flo(1):Fhi(1),Flo(2):Fhi(2),Flo(3):Fhi(3),NVAR))
    allocate(UR2(Flo(1):Fhi(1),Flo(2):Fhi(2),Flo(3):Fhi(3),NVAR))
    allocate(UR3(Flo(1):Fhi(1),Flo(2):Fhi(2),Flo(3):Fhi(3),NVAR))
    allocate(UR4(Flo(1):Fhi(1),Flo(2):Fhi(2),Flo(3):Fhi(3),NVAR))

    ! work on UL

    U0lo(1:2) = Ulo(1:2)
    U0lo(3)   = Ulo(3) + 1 ! hack: shift U to be used as base (i.e., U0))
    U0hi(1:2) = Uhi(1:2)
    U0hi(3)   = Uhi(3) + 1 ! hack

    ! obtain Gauss points in y-direction
    call reconstruct_3d(2, Glo, Ghi, &
         zlo, zhi,   &  ! for input data array
         lo, hi,     &  ! L & R not present
         Glo, Ghi,   &  ! for G1 & G2
         U0lo, U0hi, &  ! for U0
         UL,  &
         UG1=UG1, UG2=UG2, &
         U0 = U)

    ! obtain Gauss points in x-direction
    call reconstruct_3d(1, Flo,Fhi, & 
         Glo, Ghi(1), &  ! for input data array
         lo, hi,      &  ! L & R not present
         Flo, Fhi,    &  ! for G1 & G2
         U0lo, U0hi,  &  ! for U0
         UG1,  &
         UG1=UL1, UG2=UL2, &
         U0 = U)
    call reconstruct_3d(1, Flo,Fhi, & 
         Glo, Ghi(1), &  ! for input data array
         lo, hi,      &  ! L & R not present
         Flo, Fhi,    &  ! for G1 & G2
         U0lo, U0hi,  &  ! for U0
         UG2,  &
         UG1=UL3, UG2=UL4, &
         U0 = U)

    ! work on UR

    U0lo = Ulo
    U0hi = Uhi

    ! obtain Gauss points in y-direction
    call reconstruct_3d(2, Glo, Ghi, &
         zlo, zhi,   &  ! for input data array
         lo, hi,     &  ! L & R not present
         Glo, Ghi,   &  ! for G1 & G2
         U0lo, U0hi, &  ! for U0
         UR,  &
         UG1=UG1, UG2=UG2, &
         U0 = U)

    ! obtain Gauss points in x-direction
    call reconstruct_3d(1, Flo,Fhi, & 
         Glo, Ghi(1), &  ! for input data array
         lo, hi,      &  ! L & R not present
         Flo, Fhi,    &  ! for G1 & G2
         U0lo, U0hi,  &  ! for U0
         UG1,  &
         UG1=UR1, UG2=UR2, &
         U0 = U)
    call reconstruct_3d(1, Flo,Fhi, & 
         Glo, Ghi(1), &  ! for input data array
         lo, hi,      &  ! L & R not present
         Flo, Fhi,    &  ! for G1 & G2
         U0lo, U0hi,  &  ! for U0
         UG2,  &
         UG1=UR3, UG2=UR4, &
         U0 = U)

    allocate(flx(lo(3):hi(3)+1,NVAR))

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          ! 1
          call riemann(lo(3),hi(3), UL1(i,j,:,:), UR1(i,j,:,:), lo(3),hi(3)+1, &
               flx, lo(3),hi(3)+1, dir=3)
          do n=1,NVAR
             do k=lo(3),hi(3)+1
                fz(i,j,k,n) = fz(i,j,k,n) + 0.25d0*flx(k,n)
             end do
          end do

          ! 2
          call riemann(lo(3),hi(3), UL2(i,j,:,:), UR2(i,j,:,:), lo(3),hi(3)+1, &
               flx, lo(3),hi(3)+1, dir=3)
          do n=1,NVAR
             do k=lo(3),hi(3)+1
                fz(i,j,k,n) = fz(i,j,k,n) + 0.25d0*flx(k,n)
             end do
          end do

          ! 3
          call riemann(lo(3),hi(3), UL3(i,j,:,:), UR3(i,j,:,:), lo(3),hi(3)+1, &
               flx, lo(3),hi(3)+1, dir=3)
          do n=1,NVAR
             do k=lo(3),hi(3)+1
                fz(i,j,k,n) = fz(i,j,k,n) + 0.25d0*flx(k,n)
             end do
          end do

          ! 4
          call riemann(lo(3),hi(3), UL4(i,j,:,:), UR4(i,j,:,:), lo(3),hi(3)+1, &
               flx, lo(3),hi(3)+1, dir=3)
          do n=1,NVAR
             do k=lo(3),hi(3)+1
                fz(i,j,k,n) = fz(i,j,k,n) + 0.25d0*flx(k,n)
             end do
          end do

       end do
    end do

    deallocate(UL1,UL2,UL3,UL4,UR1,UR2,UR3,UR4,UG1,UG2,flx)

  end subroutine hypterm_z


  subroutine add_artifical_viscocity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, difmag

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fxlo(3), fxhi(3), &
         fylo(3), fyhi(3), fzlo(3), fzhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   ):: U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    double precision, allocatable :: divv(:,:,:), vx(:,:,:), vy(:,:,:), vz(:,:,:)
    double precision :: rhoInv, dvdx, dvdy, dvdz, div1, dxinv(3)
    integer :: i, j, k, n

    dxinv = 1.d0/dx

    allocate(vx  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(vy  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(vz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(divv(lo(1)  :hi(1)+1,lo(2)  :hi(2)+1,lo(3)  :hi(3)+1))

    do       k=lo(3)-1,hi(3)+1
       do    j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             rhoInv = 1.d0/U(i,j,k,URHO)
             vx(i,j,k) = U(i,j,k,UMX)*rhoInv
             vy(i,j,k) = U(i,j,k,UMY)*rhoInv
             vz(i,j,k) = U(i,j,k,UMZ)*rhoInv
          end do
       end do
    end do

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1

             dvdx = .25d0*dxinv(1)*( &
                    + vx(i  ,j  ,k  ) - vx(i-1,j  ,k  ) &
                    + vx(i  ,j  ,k-1) - vx(i-1,j  ,k-1) &
                    + vx(i  ,j-1,k  ) - vx(i-1,j-1,k  ) &
                    + vx(i  ,j-1,k-1) - vx(i-1,j-1,k-1) )

             dvdy = .25d0*dxinv(2)*( &
                    + vy(i  ,j  ,k  ) - vy(i  ,j-1,k  ) &
                    + vy(i  ,j  ,k-1) - vy(i  ,j-1,k-1) &
                    + vy(i-1,j  ,k  ) - vy(i-1,j-1,k  ) &
                    + vy(i-1,j  ,k-1) - vy(i-1,j-1,k-1) )

             dvdz = .25d0*dxinv(3)*( &
                    + vz(i  ,j  ,k  ) - vz(i  ,j  ,k-1) &
                    + vz(i  ,j-1,k  ) - vz(i  ,j-1,k-1) &
                    + vz(i-1,j  ,k  ) - vz(i-1,j  ,k-1) &
                    + vz(i-1,j-1,k  ) - vz(i-1,j-1,k-1) )

             divv(i,j,k) = dvdx + dvdy + dvdz

          enddo
       enddo
    enddo

    do n = 1, NVAR
       if ( n.ne.UTEMP ) then
          
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)+1
                   div1 = .25d0*(divv(i,j,k) + divv(i,j+1,k) + divv(i,j,k+1) + divv(i,j+1,k+1))
                   div1 = difmag*min(0.d0,div1)
                   fx(i,j,k,n) = fx(i,j,k,n) + dx(1)*div1*(U(i,j,k,n)-U(i-1,j,k,n))
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)+1
                do i = lo(1),hi(1)
                   div1 = .25d0*(divv(i,j,k) + divv(i+1,j,k) + divv(i,j,k+1) + divv(i+1,j,k+1))
                   div1 = difmag*min(0.d0,div1)
                   fy(i,j,k,n) = fy(i,j,k,n) + dx(2)*div1*(U(i,j,k,n)-U(i,j-1,k,n))
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)+1
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   div1 = .25d0*(divv(i,j,k) + divv(i+1,j,k) + divv(i,j+1,k) + divv(i+1,j+1,k))
                   div1 = difmag*min(0.d0,div1)
                   fz(i,j,k,n) = fz(i,j,k,n) + dx(3)*div1*(U(i,j,k,n)-U(i,j,k-1,n))
                enddo
             enddo
          enddo
          
       endif
    enddo

  end subroutine add_artifical_viscocity

end module hypterm_module
