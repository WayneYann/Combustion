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
       call reconstruct(tlo(3),thi(3), &
            Ulo(3),Uhi(3),   &  ! for U
            tlo(3),thi(3)+1, &  ! for UL & UR
            tlo(3),thi(3),   &  ! for UG1 & UG2
            0,0,             &  ! U0 is not present
            U(i,j,:,:), &
            UL=ULz(i,j,:,:), UR=URz(i,j,:,:),  &
            UG1=UG1z(i,j,:,:), UG2=UG2z(i,j,:,:), &
            dir=3)
    end do
    end do

    call hypterm_xy(0.5d0,lo,hi,UG1z,tlo,thi,fx,fxlo,fxhi,fy,fylo,fyhi,dx, &
         Ulo,Uhi,U)
    call hypterm_xy(0.5d0,lo,hi,UG2z,tlo,thi,fx,fxlo,fxhi,fy,fylo,fyhi,dx, &
         Ulo,Uhi,U)

    deallocate(UG1z,UG2z)

    ! z-direction flux


    if (difmag .gt. 0.0d0) then
       call add_numerical_viscosity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)
    end if

    deallocate(ULz,URz)

  end subroutine hypterm


  subroutine add_numerical_viscosity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)

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


  end subroutine add_numerical_viscosity

end module hypterm_module
