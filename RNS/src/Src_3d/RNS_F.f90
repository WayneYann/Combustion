
subroutine rns_dudt_ad (lo, hi, &
     U   , U_l1, U_l2, U_l3, U_h1, U_h2, U_h3, &
     dUdt,Ut_l1,Ut_l2,Ut_l3,Ut_h1,Ut_h2,Ut_h3, &
     xflx,xf_l1,xf_l2,xf_l3,xf_h1,xf_h2,xf_h3, &
     yflx,yf_l1,yf_l2,yf_l3,yf_h1,yf_h2,yf_h3, &
     zflx,zf_l1,zf_l2,zf_l3,zf_h1,zf_h2,zf_h3, &
     dx)
  use prob_params_module, only : physbc_lo, physbc_hi, NoSlipWall
  use meth_params_module, only : NVAR, gravity_dir, gravity, URHO, UMX, UEDEN, do_weno, &
       xblksize, yblksize, zblksize, nthreads
  use hypterm_module, only : hypterm
  use difterm_module, only : difterm
  use threadbox_module, only : build_threadbox_3d, get_lo_hi
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) ::  U_l1,  U_h1,  U_l2,  U_h2,  U_l3,  U_h3
  integer, intent(in) :: Ut_l1, Ut_h1, Ut_l2, Ut_h2, Ut_l3, Ut_h3
  integer, intent(in) :: xf_l1, xf_h1, xf_l2, xf_h2, xf_l3, xf_h3
  integer, intent(in) :: yf_l1, yf_h1, yf_l2, yf_h2, yf_l3, yf_h3
  integer, intent(in) :: zf_l1, zf_h1, zf_l2, zf_h2, zf_l3, zf_h3
  double precision,intent(in) ::   U( U_l1: U_h1, U_l2: U_h2, U_l3: U_h3,NVAR)
  double precision,intent(out)::dUdt(Ut_l1:Ut_h1,Ut_l2:Ut_h2,Ut_l3:Ut_h3,NVAR)
  double precision,intent(out)::xflx(xf_l1:xf_h1,xf_l2:xf_h2,xf_l3:xf_h3,NVAR)
  double precision,intent(out)::yflx(yf_l1:yf_h1,yf_l2:yf_h2,yf_l3:yf_h3,NVAR)
  double precision,intent(out)::zflx(zf_l1:zf_h1,zf_l2:zf_h2,zf_l3:zf_h3,NVAR)
  double precision,intent(in) :: dx(3)

  integer :: igrav
  integer :: Ulo(3),Uhi(3),fxlo(3),fxhi(3),fylo(3),fyhi(3),fzlo(3),fzhi(3),tlo(3),thi(3)
  integer :: iblock, nblocks, nblocksxy, iblockxy, i, j, k, n, ib, jb, kb, nb(3), boxsize(3)
  double precision :: dxinv(3)
  integer, allocatable :: bxlo(:), bxhi(:), bylo(:), byhi(:), bzlo(:), bzhi(:)
  double precision, allocatable :: bxflx_h(:,:,:,:), byflx_h(:,:,:,:), bzflx_h(:,:,:,:)
  double precision, allocatable :: bxflx_d(:,:,:,:), byflx_d(:,:,:,:), bzflx_d(:,:,:,:)

  integer, parameter :: blocksize_min = 4

  dxinv = 1.d0/dx

  Ulo(1) = U_l1
  Ulo(2) = U_l2
  Ulo(3) = U_l3
  Uhi(1) = U_h1
  Uhi(2) = U_h2
  Uhi(3) = U_h3

  boxsize = hi-lo+1

  if (nthreads > 1) then
     call build_threadbox_3d(nthreads, boxsize, blocksize_min, nb)
     if (nb(1).eq.0) then
        nb = boxsize/blocksize_min
     end if
  else
     nb(1) = max(boxsize(1)/xblksize, 1)
     nb(2) = max(boxsize(2)/yblksize, 1)
     nb(3) = max(boxsize(3)/zblksize, 1)
  end if

  allocate(bxlo(0:nb(1)-1))
  allocate(bxhi(0:nb(1)-1))
  allocate(bylo(0:nb(2)-1))
  allocate(byhi(0:nb(2)-1))
  allocate(bzlo(0:nb(3)-1))
  allocate(bzhi(0:nb(3)-1))

  call get_lo_hi(boxsize(1), nb(1), bxlo, bxhi)
  call get_lo_hi(boxsize(2), nb(2), bylo, byhi)
  call get_lo_hi(boxsize(3), nb(3), bzlo, bzhi)

  nblocksxy = nb(1)*nb(2)
  nblocks   = nb(1)*nb(2)*nb(3)

  !$omp parallel private(fxlo,fxhi,fylo,fyhi,fzlo,fzhi,tlo,thi) &
  !$omp private(iblock,iblockxy,i,j,k,n,ib,jb,kb,bxflx_h,byflx_h,bzflx_h,bxflx_d,byflx_d,bzflx_d)

  !$omp do
  do iblock = 0, nblocks-1

     kb = iblock / nblocksxy
     iblockxy = iblock - kb*nblocksxy
     jb = iblockxy / nb(1)
     ib = iblockxy - jb*nb(1)

     tlo(1) = lo(1) + bxlo(ib)
     thi(1) = lo(1) + bxhi(ib)
     
     tlo(2) = lo(2) + bylo(jb)
     thi(2) = lo(2) + byhi(jb)
     
     tlo(3) = lo(3) + bzlo(kb)
     thi(3) = lo(3) + bzhi(kb)
     
     fxlo = tlo
     fxhi(1) = thi(1)+1
     fxhi(2) = thi(2)
     fxhi(3) = thi(3)
     
     fylo = tlo
     fyhi(1) = thi(1)
     fyhi(2) = thi(2)+1
     fyhi(3) = thi(3)
     
     fzlo = tlo
     fzhi(1) = thi(1)
     fzhi(2) = thi(2)
     fzhi(3) = thi(3)+1
     
     allocate(bxflx_h(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR))
     allocate(byflx_h(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR))
     allocate(bzflx_h(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR))

     allocate(bxflx_d(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR))
     allocate(byflx_d(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR))
     allocate(bzflx_d(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR))
     
     bxflx_h = 0.d0
     byflx_h = 0.d0
     bzflx_h = 0.d0
     
     if (do_weno) then
        call hypterm(tlo,thi,U,Ulo,Uhi,bxflx_h,fxlo,fxhi,byflx_h,fylo,fyhi,bzflx_h,fzlo,fzhi,dx)
     end if

     bxflx_d = 0.d0
     byflx_d = 0.d0
     bzflx_d = 0.d0

     call difterm(tlo,thi,U,Ulo,Uhi,bxflx_d,fxlo,fxhi,byflx_d,fylo,fyhi,bzflx_d,fzlo,fzhi,dxinv)
     
     if (tlo(1) .eq. 0 .and. physbc_lo(1) .eq. NoSlipWall) then
        bxflx_d(tlo(1),:,:,:) = 0.d0
     end if

     if (tlo(2) .eq. 0 .and. physbc_lo(2) .eq. NoSlipWall) then
        bxflx_d(:,tlo(2),:,:) = 0.d0
     end if

     if (tlo(3) .eq. 0 .and. physbc_lo(3) .eq. NoSlipWall) then
        bxflx_d(:,:,tlo(3),:) = 0.d0
     end if

     ! Note that fluxes are on faces.  So don't double count!
     if (thi(1) .ne. hi(1)) fxhi(1) = fxhi(1) - 1
     if (thi(2) .ne. hi(2)) fyhi(2) = fyhi(2) - 1
     if (thi(3) .ne. hi(3)) fzhi(3) = fzhi(3) - 1
     
     do n=1,NVAR
        do    k=fxlo(3),fxhi(3)
           do j=fxlo(2),fxhi(2)
              do i=fxlo(1),fxhi(1)
                 xflx(i,j,k,n) = bxflx_h(i,j,k,n) + bxflx_d(i,j,k,n)
              end do
           end do
        end do
        
        do    k=fylo(3),fyhi(3)
           do j=fylo(2),fyhi(2)
              do i=fylo(1),fyhi(1)
                 yflx(i,j,k,n) = byflx_h(i,j,k,n) + byflx_d(i,j,k,n)
              end do
           end do
        end do
        
        do    k=fzlo(3),fzhi(3)
           do j=fzlo(2),fzhi(2)
              do i=fzlo(1),fzhi(1)
                 zflx(i,j,k,n) = bzflx_h(i,j,k,n) + bzflx_d(i,j,k,n)
              end do
           end do
        end do
     end do
     
     deallocate(bxflx_h,byflx_h,bzflx_h,bxflx_d,byflx_d,bzflx_d)
     
  end do
  !$omp end do

  !$omp do
  do n=1, NVAR
     do    k=lo(3),hi(3)
        do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           dUdt(i,j,k,n) = dxinv(1)*(xflx(i,j,k,n)-xflx(i+1,j,k,n)) &
                +          dxinv(2)*(yflx(i,j,k,n)-yflx(i,j+1,k,n)) &
                +          dxinv(3)*(zflx(i,j,k,n)-zflx(i,j,k+1,n))
        end do
        end do
     end do
  end do
  !$omp end do

  if (gravity .ne. 0.d0) then
     igrav = UMX + gravity_dir-1
     !$omp do collapse(2)
     do    k=lo(3),hi(3)
        do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           dUdt(i,j,k,igrav) = dUdt(i,j,k,igrav) + U(i,j,k,URHO)*gravity
           dUdt(i,j,k,UEDEN) = dUdt(i,j,k,UEDEN) + U(i,j,k,igrav)*gravity
        end do
        end do
     end do
     !$omp end do
  end if
  
  !$omp end parallel

end subroutine rns_dudt_ad

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_advchem(lo,hi,U,U_l1,U_l2,U_l3,U_h1,U_h2,U_h3,  &
     st,st_l1,st_l2,st_l3,st_h1,st_h2,st_h3, dt)
  use meth_params_module, only : NVAR
  use chemterm_module, only : chemterm
  implicit none
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) ::  U_l1, U_l2, U_l3, U_h1, U_h2, U_h3, &
       st_l1,st_l2,st_l3,st_h1,st_h2,st_h3
  double precision, intent(inout) :: U(U_l1:U_h1,U_l2:U_h2,U_l3:U_h3,NVAR)
  double precision, intent(inout) :: st(st_l1:st_h1,st_l2:st_h2,st_l3:st_h3)
  double precision, intent(in) :: dt

  integer :: Ulo(3), Uhi(3), stlo(3), sthi(3)

  Ulo(1) = U_l1
  Ulo(2) = U_l2
  Ulo(3) = U_l3
  Uhi(1) = U_h1
  Uhi(2) = U_h2
  Uhi(3) = U_h3
  stlo(1) = st_l1
  stlo(2) = st_l2
  stlo(3) = st_l3
  sthi(1) = st_h1
  sthi(2) = st_h2
  sthi(3) = st_h3
  call chemterm(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt)
end subroutine rns_advchem

subroutine rns_advchem2(lo,hi,U,U_l1,U_l2,U_l3,U_h1,U_h2,U_h3, &
     st,st_l1,st_l2,st_l3,st_h1,st_h2,st_h3, &
     Up,Up_l1,Up_l2,Up_l3,Up_h1,Up_h2,Up_h3,dt)
  use meth_params_module, only : NVAR
  use chemterm_module, only : chemterm
  implicit none
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) ::  U_l1, U_l2, U_l3, U_h1, U_h2, U_h3, &
       st_l1,st_l2,st_l3,st_h1,st_h2,st_h3
  integer, intent(in) ::  Up_l1, Up_l2, Up_l3, Up_h1, Up_h2, Up_h3
  double precision, intent(inout) :: U(U_l1:U_h1,U_l2:U_h2,U_l3:U_h3,NVAR)
  double precision, intent(inout) :: st(st_l1:st_h1,st_l2:st_h2,st_l3:st_h3)
  double precision, intent(in   ) :: Up(Up_l1:Up_h1,Up_l2:Up_h2,Up_l3:Up_h3,NVAR)
  double precision, intent(in) :: dt

  integer :: Ulo(3), Uhi(3), stlo(3), sthi(3)

  Ulo(1) = U_l1
  Ulo(2) = U_l2
  Ulo(3) = U_l3
  Uhi(1) = U_h1
  Uhi(2) = U_h2
  Uhi(3) = U_h3
  stlo(1) = st_l1
  stlo(2) = st_l2
  stlo(3) = st_l3
  sthi(1) = st_h1
  sthi(2) = st_h2
  sthi(3) = st_h3
  call chemterm(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt, Up)
end subroutine rns_advchem2

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_dUdt_chem(lo,hi, &
     U , U_l1, U_l2, U_l3, U_h1, U_h2, U_h3, &
     Ut,Ut_l1,Ut_l2,Ut_l3,Ut_h1,Ut_h2,Ut_h3, &
     st,st_l1,st_l2,st_l3,st_h1,st_h2,st_h3 )
  use meth_params_module, only : NVAR
  use chemterm_module, only : dUdt_chem
  implicit none
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) ::  U_l1, U_h1, Ut_l1, Ut_h1, U_l2, U_h2, Ut_l2, Ut_h2, &
       U_l3, U_h3, Ut_l3, Ut_h3, st_l1,st_l2,st_l3,st_h1,st_h2,st_h3
  double precision, intent(in ) ::  U( U_l1: U_h1, U_l2: U_h2, U_l3: U_h3,NVAR)
  double precision, intent(out) :: Ut(Ut_l1:Ut_h1,Ut_l2:Ut_h2,Ut_l3:Ut_h3,NVAR)
  double precision, intent(inout) :: st(st_l1:st_h1,st_l2:st_h2,st_l3:st_h3)

  integer :: Ulo(3), Uhi(3), Utlo(3), Uthi(3), stlo(3), sthi(3)

  Ulo(1) = U_l1
  Ulo(2) = U_l2
  Ulo(3) = U_l3
  Uhi(1) = U_h1
  Uhi(2) = U_h2
  Uhi(3) = U_h3
  Utlo(1) = Ut_l1
  Utlo(2) = Ut_l2
  Utlo(3) = Ut_l3
  Uthi(1) = Ut_h1
  Uthi(2) = Ut_h2
  Uthi(3) = Ut_h3
  stlo(1) = st_l1
  stlo(2) = st_l2
  stlo(3) = st_l3
  sthi(1) = st_h1
  sthi(2) = st_h2
  sthi(3) = st_h3
  call dUdt_chem(lo, hi, U, Ulo, Uhi, Ut, Utlo, Uthi, st, stlo, sthi)
end subroutine rns_dUdt_chem

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_compute_temp(lo,hi,U,U_l1,U_l2,U_l3,U_h1,U_h2,U_h3)
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS, NSPEC
  use eos_module, only : eos_get_T
  implicit none
  
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: U_l1, U_l2, U_l3, U_h1, U_h2, U_h3
  double precision, intent(inout) :: U(U_l1:U_h1,U_l2:U_h2,U_l3:U_h3,NVAR)

  integer :: i, j, k, pt_index(3), ierr
  double precision :: rhoInv, e, vx, vy, vz, Y(NSPEC)

  !$omp parallel do private(i,j,k,rhoInv,e,vx,vy,vz,Y,pt_index,ierr) collapse(2)
  do k=lo(3),hi(3)
  do j=lo(2),hi(2)
  do i=lo(1),hi(1)
     pt_index(1) = i
     pt_index(2) = j
     pt_index(3) = k

     rhoInv = 1.0d0/U(i,j,k,URHO)

     vx = U(i,j,k,UMX)*rhoInv     
     vy = U(i,j,k,UMY)*rhoInv     
     vz = U(i,j,k,UMZ)*rhoInv     
     e  = U(i,j,k,UEDEN)*rhoInv - 0.5d0*(vx**2+vy**2+vz**2)

     Y = U(i,j,k,UFS:UFS+NSPEC-1)*rhoInv

     call eos_get_T(U(i,j,k,UTEMP), e, Y, pt_index,ierr)

     if (ierr .ne. 0) then
        print *, 'rns_compute_temp failed at ', i,j,k,U(i,j,k,:)
        call bl_error("rns_compute_temp failed")
     end if
  end do
  end do
  end do
  !$omp end parallel do
end subroutine rns_compute_temp

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_enforce_consistent_Y(lo,hi,U,U_l1,U_l2,U_l3,U_h1,U_h2,U_h3)
  use meth_params_module, only : NVAR, URHO, UFS, NSPEC
  implicit none
  
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: U_l1, U_l2, U_l3, U_h1, U_h2, U_h3
  double precision, intent(inout) :: U(U_l1:U_h1,U_l2:U_h2,U_l3:U_h3,NVAR)

  ! Local variables
  integer          :: i,j,k,n
  integer          :: int_dom_spec
  logical          :: any_negative
  double precision :: dom_spec,x,rhoInv, sumrY, fac

  double precision, parameter :: eps = -1.d-16

  !$omp parallel do private(i,j,k,n,int_dom_spec,any_negative,dom_spec) &
  !$omp private(x,rhoInv,sumrY,fac) collapse(2)
  do k = lo(3),hi(3)
  do j = lo(2),hi(2)
  do i = lo(1),hi(1)

     any_negative = .false.

     rhoInv = 1.d0/U(i,j,k,URHO)

     sumrY = 0.d0

     ! First deal with tiny undershoots by just setting them to zero
     do n = UFS, UFS+nspec-1
        if (U(i,j,k,n) .lt. 0.d0) then
           x = U(i,j,k,n) * rhoInv
           if (x .gt. eps) then
              U(i,j,k,n) = 0.d0
           else
              any_negative = .true.
           end if
        end if

        sumrY = sumrY + U(i,j,k,n)
     end do

     fac = U(i,j,k,URHO)/sumrY
     do n = UFS, UFS+nspec-1
        U(i,j,k,n) = U(i,j,k,n)*fac
     end do

     ! We know there are one or more undershoots needing correction 
     if (any_negative) then

        ! Find the dominant species
        dom_spec = 0.d0
        int_dom_spec = 0
        do n = UFS,UFS+nspec-1
           if (U(i,j,k,n) .gt. dom_spec) then
              dom_spec = U(i,j,k,n)
              int_dom_spec = n
           end if
        end do

        ! Now take care of undershoots greater in magnitude than 1e-16.
        do n = UFS, UFS+nspec-1
           
           if (U(i,j,k,n) .lt. 0.d0) then
              
              x = U(i,j,k,n)*rhoInv
              
              ! ! Here we only print the bigger negative values
              ! if (x .lt. -1.d-2) then
              !    print *,'Correcting negative species   ',n-UFS+1
              !    print *,'   at cell (i,j,k)            ',i,j
              !    print *,'Negative (rho*Y) is           ',U(i,j,k,n)
              !    print *,'Negative      Y  is           ',x
              !    print *,'Filling from dominant species ',int_dom_spec-UFS+1
              !    print *,'  which had Y =               ',&
              !         U(i,j,k,int_dom_spec) / U(i,j,k,URHO)
              ! end if

              ! Take enough from the dominant species to fill the negative one.
              U(i,j,k,int_dom_spec) = U(i,j,k,int_dom_spec) + U(i,j,k,n)
   
              ! Test that we didn't make the dominant species negative
              if (U(i,j,k,int_dom_spec) .lt. 0.d0) then 
                 print *,' Just made dominant species negative ',int_dom_spec-UFS+1,' at ',i,j,k
                 print *,'We were fixing species ',n-UFS+1,' which had value ',x
                 print *,'Dominant species became ',U(i,j,k,int_dom_spec)*rhoInv
                 call bl_error("rns_enforce_consistent_Y")
              end if

              ! Now set the negative species to zero
              U(i,j,k,n) = 0.d0

           end if

        end do

     end if
     
  end do
  end do
  end do
  !$omp end parallel do

end subroutine rns_enforce_consistent_Y


subroutine rns_sum_cons ( &
     U  ,U_l1,U_l2,U_l3,U_h1,U_h2,U_h3, &
     msk,m_l1,m_l2,m_l3,m_h1,m_h2,m_h3, &
     vol,v_l1,v_l2,v_l3,v_h1,v_h2,v_h3, &
     s)
  use meth_params_module, only : NVAR
  implicit none
  
  integer, intent(in) :: U_l1,U_l2,U_l3,U_h1,U_h2,U_h3
  integer, intent(in) :: m_l1,m_l2,m_l3,m_h1,m_h2,m_h3
  integer, intent(in) :: v_l1,v_l2,v_l3,v_h1,v_h2,v_h3
  double precision, intent(in) :: U  (U_l1:U_h1,U_l2:U_h2,U_l3:U_h3,NVAR)
  double precision, intent(in) :: msk(m_l1:m_h1,m_l2:m_h2,m_l3:m_h3)
  double precision, intent(in) :: vol(v_l1:v_h1,v_l2:v_h2,v_l3:v_h3)
  double precision, intent(inout) :: s(5)

  integer :: i, j, k, n

  !$omp parallel do private(i,j,k,n)
  do n=1,5
     do k=m_l3,m_h3
        do j=m_l2,m_h2
           do i=m_l1,m_h1
              s(n) = s(n) + msk(i,j,k)*vol(i,j,k)*U(i,j,k,n)
           end do
        end do
     end do
  end do
  !$omp end parallel do

end subroutine rns_sum_cons


! :: ----------------------------------------------------------
! :: Volume-weight average the fine grid data onto the coarse
! :: grid.  Overlap is given in coarse grid coordinates.
! ::
! :: INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  nvar	 => number of components in arrays
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
! :: NOTE:
! ::  Assumes all data cell centered
! :: ----------------------------------------------------------
! ::
      subroutine rns_avgdown(crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,nvar, &
                             cv,cv_l1,cv_l2,cv_l3,cv_h1,cv_h2,cv_h3, &
                             fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                             fv,fv_l1,fv_l2,fv_l3,fv_h1,fv_h2,fv_h3,lo,hi,lrat)

      implicit none
      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      integer cv_l1,cv_l2,cv_l3,cv_h1,cv_h2,cv_h3
      integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      integer fv_l1,fv_l2,fv_l3,fv_h1,fv_h2,fv_h3
      integer lo(3), hi(3)
      integer nvar, lrat(3)
      double precision crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,nvar)
      double precision cv(cv_l1:cv_h1,cv_l2:cv_h2,cv_l3:cv_h3)
      double precision fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3,nvar)
      double precision fv(fv_l1:fv_h1,fv_l2:fv_h2,fv_l3:fv_h3)

      integer i, j, k, n, ic, jc, kc, ioff, joff, koff
      integer lratx, lraty, lratz
      double precision   volfrac

      lratx   = lrat(1)
      lraty   = lrat(2)
      lratz   = lrat(3)
      volfrac = 1.d0/dble(lrat(1)*lrat(2)*lrat(3))

      !$omp parallel do private(i,j,k,n,ic,jc,kc,ioff,joff,koff)
      do n = 1, nvar
         !
         ! Set coarse grid to zero on overlap.
         !
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = 0.d0
               enddo
            enddo
         enddo
         !
         ! Sum fine data.
         !
         do koff = 0, lratz-1
            do kc = lo(3),hi(3)
               k = kc*lratz + koff
               do joff = 0, lraty-1
                  do jc = lo(2), hi(2)
                     j = jc*lraty + joff
                     do ioff = 0, lratx-1
                        do ic = lo(1), hi(1)
                           i = ic*lratx + ioff
                           crse(ic,jc,kc,n) = crse(ic,jc,kc,n) + fine(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !
         ! Divide out by volume weight.
         !
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = volfrac*crse(ic,jc,kc,n)
               enddo
            enddo
         enddo

      enddo
      !$omp end parallel do

      end subroutine rns_avgdown

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine rns_estdt(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi,dx,dt)
        use eos_module, only : eos_get_c
        use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS, NSPEC
        implicit none

        integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
        integer lo(3), hi(3)
        double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
        double precision dx(3), dt

        integer :: i, j, k
        double precision :: rhoInv, vx, vy, vz, T, e, c, Y(NSPEC)

        !$omp parallel do private(i,j,k,rhoInv,vx,vy,vz,T,e,c,Y) reduction(min:dt) &
        !$omp collapse(2)
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           rhoInv = 1.d0/u(i,j,k,URHO)

           vx = u(i,j,k,UMX)*rhoInv
           vy = u(i,j,k,UMY)*rhoInv
           vz = u(i,j,k,UMZ)*rhoInv
           T  = u(i,j,k,UTEMP)
           
           e = u(i,j,k,UEDEN)*rhoInv - 0.5d0*(vx**2+vy**2+vz**2)
           
           Y = u(i,j,k,UFS:UFS+NSPEC-1)*rhoInv
           
           call eos_get_c(c,u(i,j,k,URHO),T,Y)

           dt = min(dt, dx(1)/(abs(vx)+c+1.d-50), &
                &       dx(2)/(abs(vy)+c+1.d-50), &
                &       dx(3)/(abs(vz)+c+1.d-50) )
        end do
        end do
        end do
        !$omp end parallel do

      end subroutine rns_estdt
