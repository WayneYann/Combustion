module hypterm_module

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)

    use meth_params_module, only : NVAR, difmag, do_quadrature_weno

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fxlo(3), fxhi(3), &
         fylo(3), fyhi(3), fzlo(3), fzhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   ):: U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    if (do_quadrature_weno) then
       call hypterm_q(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)
    else
       call hypterm_nq(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)
    end if
    
    if (difmag .gt. 0.0d0) then
       call add_artifical_viscocity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)
    end if

  end subroutine hypterm

  subroutine hypterm_q(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)

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

    integer :: tlo4(4), thi4(4), tlo(3), thi(3), i, j
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

  end subroutine hypterm_q


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


  subroutine hypterm_nq(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)

    use meth_params_module, only : NVAR
    use reconstruct_module, only : reconstruct
    use riemann_module, only : riemann
    use convert_module, only : cc2cellavg_2d, cellavg2cc_2d

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fxlo(3), fxhi(3), &
         fylo(3), fyhi(3), fzlo(3), fzhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   ):: U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    integer :: dir, i, j, k, n, tlo(3), thi(3)
    double precision, allocatable, dimension(:,:,:,:) :: UL_a, UR_a, UL_c, UR_c, f_c
    double precision, allocatable, dimension(:,:,:) :: f_a
    
    allocate(UL_a(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2,NVAR))
    allocate(UR_a(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2,NVAR))
    allocate(UL_c(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,NVAR))
    allocate(UR_c(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,NVAR))
    allocate( f_c(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,NVAR))
    allocate( f_a(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    !----- x-direction first -----
    dir = 1

    do k=lo(3)-2, hi(3)+2
    do j=lo(2)-2, hi(2)+2
       if ( j.ge.lo(2).and.j.le.hi(2) .or. &
            k.ge.lo(3).and.k.le.hi(3) .or. &
            j.eq.lo(2)-1.and.k.eq.lo(3)-1 .or. &
            j.eq.lo(2)-1.and.k.eq.hi(3)+1 .or. &
            j.eq.hi(2)+1.and.k.eq.lo(3)-1 .or. &
            j.eq.hi(2)+1.and.k.eq.hi(3)+1 ) then

          call reconstruct(lo(1),hi(1), & 
               Ulo(1)  , Uhi(1),   &  ! for input data array
                lo(1)-2,  hi(1)+2, &  ! for UL & UR
               0, 0,               &  ! for UG1 & UG2
               0, 0,               &  ! for U0
               U(:,j,k,:), &
               UL = UL_a(:,j,k,:), UR = UR_a (:,j,k,:), &
               dir=dir)       
       else
          cycle
!          UL_a(:,j,k,:) = 0.d0
!          UR_a(:,j,k,:) = 0.d0
       end if
    end do
    end do

    ! x-face average -> x-face cell-center
    tlo(1) = lo(1)
    tlo(2) = lo(2)-1
    tlo(3) = lo(3)-1
    thi(1) = hi(1)+1
    thi(2) = hi(2)+1
    thi(3) = hi(3)+1
    do n=1, NVAR
       call cellavg2cc_2d(tlo,thi,UL_a(:,:,:,n),lo-2,hi+2,UL_c(:,:,:,n),lo-1,hi+1,dir)
    end do
    do n=1, NVAR
       call cellavg2cc_2d(tlo,thi,UR_a(:,:,:,n),lo-2,hi+2,UR_c(:,:,:,n),lo-1,hi+1,dir)
    end do

    do k=lo(3)-1, hi(3)+1
    do j=lo(2)-1, hi(2)+1
       if ( j.eq.lo(2)-1.and.k.eq.lo(3)-1 .or. &
            j.eq.lo(2)-1.and.k.eq.hi(3)+1 .or. &
            j.eq.hi(2)+1.and.k.eq.lo(3)-1 .or. &
            j.eq.hi(2)+1.and.k.eq.hi(3)+1 ) then
          cycle
          !          f_c(:,j,k,:) = 0.d0
       else
          call riemann(lo(1),hi(1),UL_c(:,j,k,:),UR_c(:,j,k,:),lo(1)-1,hi(1)+1, &
               f_c(:,j,k,:), lo(1)-1,hi(1)+1, dir=dir)
       end if
    end do
    end do

    ! x-flux: x-face cell-center --> x-face average
    tlo(1) = lo(1)
    tlo(2) = lo(2)
    tlo(3) = lo(3)
    thi(1) = hi(1)+1
    thi(2) = hi(2)
    thi(3) = hi(3)
    do n=1, NVAR
       call cc2cellavg_2d(tlo,thi,f_c(:,:,:,n),lo-1,hi+1,f_a,lo-1,hi+1,dir)
       fx      (lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = &
            fx (lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) + &
            f_a(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3))
    end do

    !----- y-direction -----
    dir = 2

    do k=lo(3)-2, hi(3)+2
    do i=lo(1)-2, hi(1)+2
       if ( i.ge.lo(1).and.i.le.hi(1) .or. &
            k.ge.lo(3).and.k.le.hi(3) .or. &
            i.eq.lo(1)-1.and.k.eq.lo(3)-1 .or. &
            i.eq.lo(1)-1.and.k.eq.hi(3)+1 .or. &
            i.eq.hi(1)+1.and.k.eq.lo(3)-1 .or. &
            i.eq.hi(1)+1.and.k.eq.hi(3)+1 ) then

          call reconstruct(lo(2),hi(2), & 
               Ulo(2)  , Uhi(2),   &  ! for input data array
                lo(2)-2,  hi(2)+2, &  ! for UL & UR
               0, 0,               &  ! for UG1 & UG2
               0, 0,               &  ! for U0
               U(i,:,k,:), &
               UL = UL_a(i,:,k,:), UR = UR_a (i,:,k,:), &
               dir=dir)       
       else
          cycle
!          UL_a(i,:,k,:) = 0.d0
!          UR_a(i,:,k,:) = 0.d0
       end if
    end do
    end do

    ! y-face average -> y-face cell-center
    tlo(1) = lo(1)-1
    tlo(2) = lo(2)
    tlo(3) = lo(3)-1
    thi(1) = hi(1)+1
    thi(2) = hi(2)+1
    thi(3) = hi(3)+1
    do n=1, NVAR
       call cellavg2cc_2d(tlo,thi,UL_a(:,:,:,n),lo-2,hi+2,UL_c(:,:,:,n),lo-1,hi+1,dir)
    end do
    do n=1, NVAR
       call cellavg2cc_2d(tlo,thi,UR_a(:,:,:,n),lo-2,hi+2,UR_c(:,:,:,n),lo-1,hi+1,dir)
    end do

    do k=lo(3)-1, hi(3)+1
    do i=lo(1)-1, hi(1)+1
       if ( i.eq.lo(1)-1.and.k.eq.lo(3)-1 .or. &
            i.eq.lo(1)-1.and.k.eq.hi(3)+1 .or. &
            i.eq.hi(1)+1.and.k.eq.lo(3)-1 .or. &
            i.eq.hi(1)+1.and.k.eq.hi(3)+1 ) then
          cycle
          !          f_c(i,:,k,:) = 0.d0
       else
          call riemann(lo(2),hi(2),UL_c(i,:,k,:),UR_c(i,:,k,:),lo(2)-1,hi(2)+1, &
               f_c(i,:,k,:), lo(2)-1,hi(2)+1, dir=dir)
       end if
    end do
    end do

    ! y-flux: y-face cell-center --> y-face average
    tlo(1) = lo(1)
    tlo(2) = lo(2)
    tlo(3) = lo(3)
    thi(1) = hi(1)
    thi(2) = hi(2)+1
    thi(3) = hi(3)
    do n=1, NVAR
       call cc2cellavg_2d(tlo,thi,f_c(:,:,:,n),lo-1,hi+1,f_a,lo-1,hi+1,dir)
       fy      (lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = &
            fy (lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) + &
            f_a(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3))
    end do

    !----- z-direction -----
    dir = 3

    do j=lo(2)-2, hi(2)+2
    do i=lo(1)-2, hi(1)+2
       if ( i.ge.lo(1).and.i.le.hi(1) .or. &
            j.ge.lo(2).and.j.le.hi(2) .or. &
            i.eq.lo(1)-1.and.j.eq.lo(2)-1 .or. &
            i.eq.lo(1)-1.and.j.eq.hi(2)+1 .or. &
            i.eq.hi(1)+1.and.j.eq.lo(2)-1 .or. &
            i.eq.hi(1)+1.and.j.eq.hi(2)+1 ) then

          call reconstruct(lo(3),hi(3), & 
               Ulo(3)  , Uhi(3),   &  ! for input data array
                lo(3)-2,  hi(3)+2, &  ! for UL & UR
               0, 0,               &  ! for UG1 & UG2
               0, 0,               &  ! for U0
               U(i,j,:,:), &
               UL = UL_a(i,j,:,:), UR = UR_a (i,j,:,:), &
               dir=dir)       
       else
          cycle
!          UL_a(i,j,:,:) = 0.d0
!          UR_a(i,j,:,:) = 0.d0
       end if
    end do
    end do

    ! z-face average -> z-face cell-center
    tlo(1) = lo(1)-1
    tlo(2) = lo(2)-1
    tlo(3) = lo(3)
    thi(1) = hi(1)+1
    thi(2) = hi(2)+1
    thi(3) = hi(3)+1
    do n=1, NVAR
       call cellavg2cc_2d(tlo,thi,UL_a(:,:,:,n),lo-2,hi+2,UL_c(:,:,:,n),lo-1,hi+1,dir)
    end do
    do n=1, NVAR
       call cellavg2cc_2d(tlo,thi,UR_a(:,:,:,n),lo-2,hi+2,UR_c(:,:,:,n),lo-1,hi+1,dir)
    end do

    do j=lo(2)-1, hi(2)+1
    do i=lo(1)-1, hi(1)+1
       if ( i.eq.lo(1)-1.and.j.eq.lo(2)-1 .or. &
            i.eq.lo(1)-1.and.j.eq.hi(2)+1 .or. &
            i.eq.hi(1)+1.and.j.eq.lo(2)-1 .or. &
            i.eq.hi(1)+1.and.j.eq.hi(2)+1 ) then
          cycle
          !          f_c(i,j,:,:) = 0.d0
       else
          call riemann(lo(3),hi(3),UL_c(i,j,:,:),UR_c(i,j,:,:),lo(3)-1,hi(3)+1, &
               f_c(i,j,:,:), lo(3)-1,hi(3)+1, dir=dir)
       end if
    end do
    end do

    ! z-flux: z-face cell-center --> z-face average
    tlo(1) = lo(1)
    tlo(2) = lo(2)
    tlo(3) = lo(3)
    thi(1) = hi(1)
    thi(2) = hi(2)
    thi(3) = hi(3)+1
    do n=1, NVAR
       call cc2cellavg_2d(tlo,thi,f_c(:,:,:,n),lo-1,hi+1,f_a,lo-1,hi+1,dir)
       fz      (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = &
            fz (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) + &
            f_a(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1)
    end do

    deallocate(UL_a, UR_a, UL_c, UR_c, f_c, f_a)

  end subroutine hypterm_nq


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
