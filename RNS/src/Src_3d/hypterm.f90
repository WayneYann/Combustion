module hypterm_module

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)

    use meth_params_module, only : NVAR, difmag
    use reconstruct_module, only : reconstruct
    use hypterm_xy_module, only : hypterm_xy

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fxlo(3), fxhi(3), &
         fylo(3), fyhi(3), fzlo(3), fzhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   ):: U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    integer :: tlo4(4), thi4(4), tlo(3), thi(3), i, j, k, n, g
    double precision, dimension(:,:,:,:), allocatable :: Ulz, URz, UG1z, UG2z

    allocate( ULz(lo(1)-3:hi(1)+3,lo(2)-3:hi(2)+3,lo(3):hi(3)+1,NVAR))
    allocate( URz(lo(1)-3:hi(1)+3,lo(2)-3:hi(2)+3,lo(3):hi(3)+1,NVAR))
    allocate(UG1z(lo(1)-3:hi(1)+3,lo(2)-3:hi(2)+3,lo(3):hi(3)  ,NVAR))
    allocate(UG2z(lo(1)-3:hi(1)+3,lo(2)-3:hi(2)+3,lo(3):hi(3)  ,NVAR))

    ! Given cell averages, reconstruct in z-direction
    ! Note that they are still averges in x and y-direction
    do j=lo(2)-3,hi(2)+3
    do i=lo(1)-3,hi(1)+3
       call reconstruct(lo(3),hi(3), &
            Ulo(3),Uhi(3),   &  ! for U
            lo (3), hi(3)+1, &  ! for UL & UR
            lo (3), hi(3),   &  ! for UG1 & UG2
            0,0,             &  ! U0 is not present
            U(i,j,:,:), &
            UL=ULz(i,j,:,:), UR=URz(i,j,:,:),  &
            UG1=UG1z(i,j,:,:), UG2=UG2z(i,j,:,:), &
            dir=3)
    end do
    end do

    tlo4 = lbound(UG1z)
    thi4 = ubound(UG1z)
    tlo = tlo4(1:3)
    thi = thi4(1:3)

    call hypterm_xy(0.5d0,lo,hi,UG1z,tlo,thi,fx,fxlo,fxhi,fy,fylo,fyhi,dx, &
         Ulo,Uhi,U)
    call hypterm_xy(0.5d0,lo,hi,UG2z,tlo,thi,fx,fxlo,fxhi,fy,fylo,fyhi,dx, &
         Ulo,Uhi,U)

    deallocate(UG1z,UG2z)

    tlo4 = lbound(ULz)
    thi4 = ubound(ULz)
    tlo = tlo4(1:3)
    thi = thi4(1:3)

    call hypterm_z(lo,hi,U,Ulo,Uhi,ULz,URz,tlo,thi,fz,fzlo,fzhi,dx)

    deallocate(ULz,URz)

    if (difmag .gt. 0.0d0) then
       call add_artifical_viscocity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)
    end if

  end subroutine hypterm


  subroutine hypterm_z(lo,hi,U,Ulo,Uhi,UL,UR,zlo,zhi,fz,fzlo,fzhi,dx)

    use meth_params_module, only : NVAR
    use reconstruct_module, only : reconstruct
    use riemann_module, only : riemann

    integer, intent(in) :: lo(3),hi(3),Ulo(3),Uhi(3),zlo(3),zhi(3),fzlo(3),fzhi(3)
    double precision,intent(in) :: dx(3)
    double precision,intent(in   )::U ( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(in   )::UL( zlo(1): zhi(1), zlo(2): zhi(2), zlo(3): zhi(3),NVAR)
    double precision,intent(in   )::UR( zlo(1): zhi(1), zlo(2): zhi(2), zlo(3): zhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    integer :: i,j,k,n
    double precision, dimension(:,:,:,:), allocatable :: UL11,UL12,UL21,UL22,UR11,UR12,UR21,UR22
    double precision, dimension(:,:,:), allocatable :: UG1y, UG2y
    double precision, dimension(:,:), allocatable :: U0, flx

    allocate(UL11(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,NVAR))
    allocate(UL12(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,NVAR))
    allocate(UL21(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,NVAR))
    allocate(UL22(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,NVAR))

    allocate(UR11(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,NVAR))
    allocate(UR12(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,NVAR))
    allocate(UR21(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,NVAR))
    allocate(UR22(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,NVAR))

    allocate(UG1y(lo(1)-2:hi(1)+2,lo(2):hi(2),NVAR))
    allocate(UG2y(lo(1)-2:hi(1)+2,lo(2):hi(2),NVAR))

    allocate(U0(lo(1):hi(1),NVAR))

    allocate(flx(lo(3):hi(3)+1,NVAR))

    do k=lo(3),hi(3)+1
       
       ! ----- UL -----

       ! obtain Gauss points in y-direction
       do i=lo(1)-2,hi(1)+2
          call reconstruct(lo(2),hi(2), & 
               zlo(2),zhi(2),   &  ! for input data array
               0,0,             &  ! L & R not present
               lo (2), hi(2),   &  ! for G1 & G2
               lo (2), hi(2),   &  ! for U0
               UL(i,:,k,:), &
               UG1=UG1y(i,:,:), UG2=UG2y(i,:,:), &
               U0 = U(i,lo(2):hi(2),k-1,:), &
               dir=2)
       end do

       ! obtain Gauss points in x-direction
       do j=lo(2),hi(2)
          U0 = U(lo(1):hi(1),j,k-1,:)
          call reconstruct(lo(1),hi(1), & 
               lo(1)-2,hi(1)+2, &  ! for input data array
               0,0,             &  ! L & R not present
               lo (1), hi(1),   &  ! for G1 & G2
               lo (1), hi(1),   &  ! for U0
               UG1y(:,j,:), &
               UG1=UL11(:,j,k,:), UG2=UL12(:,j,k,:), &
               U0 = U0, &
               dir=1)          
          call reconstruct(lo(1),hi(1), & 
               lo(1)-2,hi(1)+2, &  ! for input data array
               0,0,             &  ! L & R not present
               lo (1), hi(1),   &  ! for G1 & G2
               lo (1), hi(1),   &  ! for U0
               UG2y(:,j,:), &
               UG1=UL21(:,j,k,:), UG2=UL22(:,j,k,:), &
               U0 = U0, &
               dir=1)          
       end do

       ! ----- UR -----

       ! obtain Gauss points in y-direction
       do i=lo(1)-2,hi(1)+2
          call reconstruct(lo(2),hi(2), & 
               zlo(2),zhi(2),   &  ! for input data array
               0,0,             &  ! L & R not present
               lo (2), hi(2),   &  ! for G1 & G2
               lo (2), hi(2),   &  ! for U0
               UR(i,:,k,:), &
               UG1=UG1y(i,:,:), UG2=UG2y(i,:,:), &
               U0 = U(i,lo(2):hi(2),k,:), &
               dir=2)
       end do

       ! obtain Gauss points in x-direction
       do j=lo(2),hi(2)
          U0 = U(lo(1):hi(1),j,k,:)
          call reconstruct(lo(1),hi(1), & 
               lo(1)-2,hi(1)+2, &  ! for input data array
               0,0,             &  ! L & R not present
               lo (1), hi(1),   &  ! for G1 & G2
               lo (1), hi(1),   &  ! for U0
               UG1y(:,j,:), &
               UG1=UR11(:,j,k,:), UG2=UR12(:,j,k,:), &
               U0 = U0, &
               dir=1)          
          call reconstruct(lo(1),hi(1), & 
               lo(1)-2,hi(1)+2, &  ! for input data array
               0,0,             &  ! L & R not present
               lo (1), hi(1),   &  ! for G1 & G2
               lo (1), hi(1),   &  ! for U0
               UG2y(:,j,:), &
               UG1=UR21(:,j,k,:), UG2=UR22(:,j,k,:), &
               U0 = U0, &
               dir=1)          
       end do

    end do

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          ! 11
          call riemann(lo(3),hi(3), UL11(i,j,:,:), UR11(i,j,:,:), lo(3),hi(3)+1, &
               flx, lo(3),hi(3)+1, dir=3)
          do n=1,NVAR
             do k=lo(3),hi(3)+1
                fz(i,j,k,n) = fz(i,j,k,n) + 0.25d0*flx(k,n)
             end do
          end do

          ! 12
          call riemann(lo(3),hi(3), UL12(i,j,:,:), UR12(i,j,:,:), lo(3),hi(3)+1, &
               flx, lo(3),hi(3)+1, dir=3)
          do n=1,NVAR
             do k=lo(3),hi(3)+1
                fz(i,j,k,n) = fz(i,j,k,n) + 0.25d0*flx(k,n)
             end do
          end do

          ! 21
          call riemann(lo(3),hi(3), UL21(i,j,:,:), UR21(i,j,:,:), lo(3),hi(3)+1, &
               flx, lo(3),hi(3)+1, dir=3)
          do n=1,NVAR
             do k=lo(3),hi(3)+1
                fz(i,j,k,n) = fz(i,j,k,n) + 0.25d0*flx(k,n)
             end do
          end do

          ! 22
          call riemann(lo(3),hi(3), UL22(i,j,:,:), UR22(i,j,:,:), lo(3),hi(3)+1, &
               flx, lo(3),hi(3)+1, dir=3)
          do n=1,NVAR
             do k=lo(3),hi(3)+1
                fz(i,j,k,n) = fz(i,j,k,n) + 0.25d0*flx(k,n)
             end do
          end do

       end do
    end do

    deallocate(UL11,UL12,UL21,UL22,UR11,UR12,UR21,UR22,UG1y,UG2y,U0,flx)

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
