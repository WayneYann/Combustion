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

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UTEMP, difmag

    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    double precision, intent(in   ) :: dx(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),NVAR)
    
    double precision, allocatable :: divv(:,:), vx(:,:), vy(:,:)
    double precision :: rhoInv, dvdx, dvdy, div1, dxinv(2)
    integer :: i, j, n

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

  end subroutine add_artifical_viscocity

end module hypterm_module

