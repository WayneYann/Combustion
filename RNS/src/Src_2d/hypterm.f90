module hypterm_module

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx)

    use meth_params_module, only : NVAR, difmag, do_quadrature_weno
    use hypterm_xy_module, only : hypterm_xy

    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    double precision, intent(in   ) :: dx(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),NVAR)

    integer :: tlo(3), thi(3), tUlo(3), tUhi(3), tfxlo(3), tfxhi(3), tfylo(3), tfyhi(3)
    double precision :: tdx(3)

    if (do_quadrature_weno) then

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

    else
       
       call hypterm_nq(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx)       

    end if

    if (difmag .gt. 0.0d0) then
       call add_artifical_viscocity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx)
    end if

  end subroutine hypterm


  subroutine hypterm_nq(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx)

    use meth_params_module, only : NVAR
    use reconstruct_module, only : reconstruct
    use convert_module, only : cellavg2cc_1d, cc2cellavg_1d
    use riemann_module, only : riemann

    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    double precision, intent(in   ) :: dx(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),NVAR)

    integer :: i, j, n, tlo(2), thi(2), dir
    double precision, allocatable, dimension(:,:,:) :: UL_a, UR_a, UL_c, UR_c, f_c
    double precision, allocatable, dimension(:,:) :: f_a

    allocate(UL_a(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,NVAR))
    allocate(UR_a(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,NVAR))
    allocate(UL_c(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,NVAR))
    allocate(UR_c(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,NVAR))
    allocate( f_c(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,NVAR))
    allocate( f_a(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))

    !----- x-direction first -----
    dir = 1
    
    do j=lo(2)-2, hi(2)+2
       call reconstruct(lo(1),hi(1), & 
            Ulo(1)  , Uhi(1),   &  ! for input data array
             lo(1)-2,  hi(1)+2, &  ! for UL & UR
            0, 0,               &  ! for UG1 & UG2
            0, 0,               &  ! for U0
            U(:,j,:), &
            UL = UL_a(:,j,:), UR = UR_a (:,j,:), &
            dir=dir)       
    end do

    ! y-average --> y-cell-center
    tlo(1) = lo(1)
    tlo(2) = lo(2)-1
    thi(1) = hi(1)+1
    thi(2) = hi(2)+1
    do n=1, NVAR
       call cellavg2cc_1d(tlo,thi,UL_a(:,:,n),lo-2,hi+2,UL_c(:,:,n),lo-1,hi+1,dir)
    end do
    do n=1, NVAR
       call cellavg2cc_1d(tlo,thi,UR_a(:,:,n),lo-2,hi+2,UR_c(:,:,n),lo-1,hi+1,dir)
    end do

    do j=lo(2)-1, hi(2)+1
       call riemann(lo(1),hi(1),UL_c(:,j,:),UR_c(:,j,:),lo(1)-1,hi(1)+1, &
            f_c(:,j,:), lo(1)-1,hi(1)+1, dir=dir)
    end do

    ! x-flux: y-cell-center --> y-average
    tlo(1) = lo(1)
    tlo(2) = lo(2)
    thi(1) = hi(1)+1
    thi(2) = hi(2)
    do n=1, NVAR
       call cc2cellavg_1d(tlo,thi,f_c(:,:,n),lo-1,hi+1,f_a,lo-1,hi+1,dir)
       fx(lo(1):hi(1)+1,lo(2):hi(2),n) = fx(lo(1):hi(1)+1,lo(2):hi(2),n) &
            +                           f_a(lo(1):hi(1)+1,lo(2):hi(2))
    end do

    !----- y-direction -----
    dir = 2

    do i=lo(1)-2, hi(1)+2
       call reconstruct(lo(2),hi(2), &
            Ulo(2)  , Uhi(2),   &  ! for input data array
             lo(2)-2,  hi(2)+2, &  ! for UL & UR
            0, 0,               &  ! for UG1 & UG2
            0, 0,               &  ! for U0
            U(i,:,:), &
            UL = UL_a(i,:,:), UR = UR_a (i,:,:), &
            dir=dir)       
    end do

    ! x-average --> x-cell-center
    tlo(1) = lo(1)-1
    tlo(2) = lo(2)
    thi(1) = hi(1)+1
    thi(2) = hi(2)+1
    do n=1, NVAR
       call cellavg2cc_1d(tlo,thi,UL_a(:,:,n),lo-2,hi+2,UL_c(:,:,n),lo-1,hi+1,dir)
    end do
    do n=1, NVAR
       call cellavg2cc_1d(tlo,thi,UR_a(:,:,n),lo-2,hi+2,UR_c(:,:,n),lo-1,hi+1,dir)
    end do
    
    do i=lo(1)-1, hi(1)+1
       call riemann(lo(2),hi(2),UL_c(i,:,:),UR_c(i,:,:),lo(2)-1,hi(2)+1, &
            f_c(i,:,:), lo(2)-1,hi(2)+1, dir=dir)
    end do

    ! y-flux: x-cell-center --> x-average
    tlo(1) = lo(1)
    tlo(2) = lo(2)
    thi(1) = hi(1)
    thi(2) = hi(2)+1
    do n=1, NVAR
       call cc2cellavg_1d(tlo,thi,f_c(:,:,n),lo-1,hi+1,f_a,lo-1,hi+1,dir)
       fy(lo(1):hi(1),lo(2):hi(2)+1,n) = fy(lo(1):hi(1),lo(2):hi(2)+1,n) &
            +                           f_a(lo(1):hi(1),lo(2):hi(2)+1)
    end do

    deallocate(UL_a, UR_a, UL_c, UR_c, f_c, f_a)
    
  end subroutine hypterm_nq


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

