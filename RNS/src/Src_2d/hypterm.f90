module hypterm_module

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx)

    use meth_params_module, only : NVAR, difmag
    use hypterm_xy_module, only : hypterm_xy

    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    double precision, intent(in   ) :: dx(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),NVAR)

    integer :: tlo(3), thi(3), tUlo(3), tUhi(3), tfxlo(3), tfxhi(3), tfylo(3), tfyhi(3)
    double precision :: tdx(3)

    tlo(1:2) = lo
    tlo(3) = 1
    thi(1:2) = hi
    thi(3) = 1
    
    tUlo(1:2) = Ulo
    tUlo(3) = 1
    tUhi(1:2) = Uhi
    tUhi(3) = 1
    
    tfxlo(1:2) = fxlo
    tfxlo(3) = 1
    tfxhi(1:2) = fxhi
    tfxhi(3) = 1
    
    tfylo(1:2) = fylo
    tfylo(3) = 1
    tfyhi(1:2) = fyhi
    tfyhi(3) = 1
    
    tdx(1:2) = dx
    tdx(3) = 0.d0
    
    call hypterm_xy(1.d0,tlo,thi,U,tUlo,tUhi,fx,tfxlo,tfxhi,fy,tfylo,tfyhi,tdx,tlo,thi)

    if (difmag .gt. 0.0d0) then
       call add_artifical_viscocity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx)
    end if

  end subroutine hypterm


  subroutine add_artifical_viscocity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx)

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UTEMP, UEDEN, UFS, NSPEC, difmag
    use eos_module, only : eos_get_T, eos_get_c

    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    double precision, intent(in   ) :: dx(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),NVAR)
    
    double precision, allocatable :: vx(:,:), vy(:,:), cs(:,:)
    double precision :: rhoInv, hdvdx, hdvdy, hdivv, dxinv(2), e, T, Y(nspec)
    double precision :: nu, fac, cmin
    integer :: i, j, n

    dxinv = 1.d0/dx
    
    allocate(vx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(vy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(cs(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    
    do    j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          rhoInv = 1.d0/U(i,j,URHO)
          vx(i,j) = U(i,j,UMX)*rhoInv
          vy(i,j) = U(i,j,UMY)*rhoInv
          e  = U(i,j,UEDEN)*rhoInv - 0.5d0*(vx(i,j)**2+vy(i,j)**2)
          Y = U(i,j,UFS:UFS+NSPEC-1)*rhoInv
          T = U(i,j,UTEMP)
          call eos_get_T(T, e, Y)
          call eos_get_c(cs(i,j), U(i,j,URHO), T, Y)
       end do
    end do

    ! x-direction

    fac = 0.25d0*dx(1)*dxinv(2)
    
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          hdvdx = vx(i,j)-vx(i-1,j)
          hdvdy = fac*(-vy(i-1,j-1)-vy(i,j-1)+vy(i-1,j+1)+vy(i,j+1))
          hdivv = hdvdx + hdvdy
          if (hdivv .lt. 0.d0) then
             cmin = min(cs(i-1,j),cs(i,j))
             nu = difmag * hdivv * min(hdivv*hdivv/(cmin*cmin*difmag), 1.d0)
             do n=1,NVAR
                if (n.ne.UTEMP) then
                   fx(i,j,n) = fx(i,j,n) - nu*(U(i-1,j,n)-U(i,j,n))
                end if
             end do
          end if
       end do
    end do

    ! y-direction

    fac = 0.25d0*dxinv(1)*dx(2)
    
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          hdvdx = fac*(-vx(i-1,j-1)-vx(i-1,j)+vx(i+1,j-1)+vx(i+1,j))
          hdvdy = vy(i,j)-vy(i,j-1)
          hdivv = hdvdx + hdvdy
          if (hdivv .lt. 0.d0) then
             cmin = min(cs(i,j-1),cs(i,j))
             nu = difmag * hdivv * min(hdivv*hdivv/(cmin*cmin*difmag), 1.d0)
             do n=1,NVAR
                if (n.ne.UTEMP) then
                   fy(i,j,n) = fy(i,j,n) - nu*(U(i,j-1,n)-U(i,j,n))
                end if
             end do
          end if
       end do
    end do

    deallocate(vx,vy,cs)

  end subroutine add_artifical_viscocity

end module hypterm_module

