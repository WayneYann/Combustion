module hypterm_module

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)

    use meth_params_module, only : NVAR
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


    ! difmag

    deallocate(ULz,URz)

  end subroutine hypterm

end module hypterm_module
