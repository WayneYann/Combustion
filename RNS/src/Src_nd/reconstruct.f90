module reconstruct_module

  use meth_params_module, only : ndim, NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, &
       UFS, NSPEC, NCHARV, CFS, do_component_weno
  use weno_module, only : weno5, vweno5, weno_3d
  use eos_module, only : eos_given_ReY, eos_get_eref

  implicit none

  private

  public :: reconstruct, reconstruct_3d

contains

  ! L and R in UL and UR are relative to face
  ! UG1 and UG2 are at two Gauss points
  !
  ! Minimal ranges in the z-direction:
  !   U         : lo-2:hi+2 if UL & UR are not present; 
  !               lo-3:hi+3 if UL & UR are present
  !   UL  & UR  : lo  :hi+1
  !   UG1 & UG2 : lo  :hi
  !   U0        : lo  :hi   if UL & UR are not present; 
  !               lo-1:hi+1 if UL & UR are present

  subroutine reconstruct_3d(idir, lo, hi, &
       Ulo, Uhi, &
       ULRlo, ULRhi, &
       UGlo, UGhi, &
       U0lo, U0hi, &
       U, UL, UR, UG1, UG2, U0)

    integer, intent(in) :: idir, lo(3), hi(3), Ulo(3), Uhi(3), ULRlo(3), ULRhi(3), &
         UGlo(3), UGhi(3), U0lo(3), U0hi(3) 
    double precision, intent(in),target           :: &
         &            U  (  Ulo(1):  Uhi(1),  Ulo(2):  Uhi(2),  Ulo(3):  Uhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UL (ULRlo(1):ULRhi(1),ULRlo(2):ULRhi(2),ULRlo(3):ULRhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UR (ULRlo(1):ULRhi(1),ULRlo(2):ULRhi(2),ULRlo(3):ULRhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UG1( UGlo(1): UGhi(1), UGlo(2): UGhi(2), UGlo(3): UGhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UG2( UGlo(1): UGhi(1), UGlo(2): UGhi(2), UGlo(3): UGhi(3),NVAR)
    double precision, intent(in),target, optional :: &
         &            U0 ( U0lo(1): U0hi(1), U0lo(2): U0hi(2), U0lo(3): U0hi(3),NVAR)

    if (do_component_weno) then
       call reconstruct_comp_3d(idir, lo, hi, &
            Ulo, Uhi, &
            ULRlo, ULRhi, &
            UGlo, UGhi, &
            U, UL, UR, UG1, UG2)
    else
       call reconstruct_char_3d(idir, lo, hi, &
            Ulo, Uhi, &
            ULRlo, ULRhi, &
            UGlo, UGhi, &
            U0lo, U0hi, &
            U, UL, UR, UG1, UG2, U0)
    end if

  end subroutine reconstruct_3d


  subroutine reconstruct_comp_3d(idir, lo, hi, &
       Ulo, Uhi, &
       ULRlo, ULRhi, &
       UGlo, UGhi, &
       U, UL, UR, UG1, UG2)

    integer, intent(in) :: idir, lo(3), hi(3), Ulo(3), Uhi(3), ULRlo(3), ULRhi(3), &
         UGlo(3), UGhi(3)
    double precision, intent(in),target           :: &
         &            U  (  Ulo(1):  Uhi(1),  Ulo(2):  Uhi(2),  Ulo(3):  Uhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UL (ULRlo(1):ULRhi(1),ULRlo(2):ULRhi(2),ULRlo(3):ULRhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UR (ULRlo(1):ULRhi(1),ULRlo(2):ULRhi(2),ULRlo(3):ULRhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UG1( UGlo(1): UGhi(1), UGlo(2): UGhi(2), UGlo(3): UGhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UG2( UGlo(1): UGhi(1), UGlo(2): UGhi(2), UGlo(3): UGhi(3),NVAR)

    integer :: i, j, k, n, elo(3), ehi(3)
    double precision :: rho, rhoInv, Y(NSPEC)
    double precision, allocatable :: rhoE(:,:,:)
    logical :: do_gauss, do_face

    elo = lo
    elo(idir) = elo(idir) - 2
    ehi = hi
    ehi(idir) = ehi(idir) + 2

    do_face = .false.
    do_gauss = .false.
    
    if (present(UL) .and. present(UR)) then
       do_face = .true.
       elo(idir) = elo(idir)-1
       ehi(idir) = ehi(idir)+1
    end if
    
    if (present(UG1) .and. present(UG2)) then
       do_gauss = .true.
    end if

    allocate(rhoE(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3)))

    do k      =elo(3), ehi(3)
       do j   =elo(2), ehi(2)
          do i=elo(1), ehi(1)
             
             rho = 0.d0
             do n=1,nspec
                Y(n) = U(i,j,k,UFS+n-1)
                rho = rho + Y(n)
             end do
             
             rhoInv = 1.d0/rho
             Y = Y * rhoInv

             rhoE(i,j,k) = U(i,j,k,UEDEN) - rho*eos_get_eref(Y)

          end do
       end do
    end do

    if (do_face .and. do_gauss) then

       do n=1,NVAR
          if (n.eq.UEDEN) then
             call weno_3d(idir,lo,hi, rhoE,elo,ehi, ULRlo, ULRhi, UGlo, UGhi, &
                  UL, UR, UG1, UG2)
          else
             call weno_3d(idir,lo,hi, U(:,:,:,n),Ulo,Uhi, ULRlo, ULRhi, UGlo, UGhi, &
                  UL(:,:,:,n), UR(:,:,:,n), UG1(:,:,:,n), UG2(:,:,:,n))
          end if
       end do

    else if (do_face) then

       do n=1,NVAR
          if (n.eq.UEDEN) then
             call weno_3d(idir,lo,hi, rhoE,elo,ehi, ULRlo, ULRhi, UGlo, UGhi, &
                  vp=UL, vm=UR)
          else
             call weno_3d(idir,lo,hi, U(:,:,:,n),Ulo,Uhi, ULRlo, ULRhi, UGlo, UGhi, &
                  vp=UL(:,:,:,n), vm=UR(:,:,:,n))
          end if
       end do
       
    else ! do_gauss

       do n=1,NVAR
          if (n.eq.UEDEN) then
             call weno_3d(idir,lo,hi, rhoE,elo,ehi, ULRlo, ULRhi, UGlo, UGhi, &
                  vg1=UG1, vg2=UG2)
          else
             call weno_3d(idir,lo,hi, U(:,:,:,n),Ulo,Uhi, ULRlo, ULRhi, UGlo, UGhi, &
                  vg1=UG1(:,:,:,n), vg2=UG2(:,:,:,n))
          end if
       end do

    end if

    if (do_face) then
       elo(idir) = elo(idir) + 3
       ehi(idir) = ehi(idir) - 2
       call normalize_3d(elo,ehi, UL, ULRlo,ULRhi)
       call normalize_3d(elo,ehi, UR, ULRlo,ULRhi)
    end if

    if (do_gauss) then
       call normalize_3d(lo,hi, UG1, UGlo,UGhi)
       call normalize_3d(lo,hi, UG2, UGlo,UGhi)
    end if

    deallocate(rhoE)

  end subroutine reconstruct_comp_3d


  subroutine reconstruct_char_3d(idir, lo, hi, &
       Ulo, Uhi, &
       ULRlo, ULRhi, &
       UGlo, UGhi, &
       U0lo, U0hi, &
       U, UL, UR, UG1, UG2, U0)

    integer, intent(in) :: idir, lo(3), hi(3), Ulo(3), Uhi(3), ULRlo(3), ULRhi(3), &
         UGlo(3), UGhi(3), U0lo(3), U0hi(3) 
    double precision, intent(in),target           :: &
         &            U  (  Ulo(1):  Uhi(1),  Ulo(2):  Uhi(2),  Ulo(3):  Uhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UL (ULRlo(1):ULRhi(1),ULRlo(2):ULRhi(2),ULRlo(3):ULRhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UR (ULRlo(1):ULRhi(1),ULRlo(2):ULRhi(2),ULRlo(3):ULRhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UG1( UGlo(1): UGhi(1), UGlo(2): UGhi(2), UGlo(3): UGhi(3),NVAR)
    double precision, intent(out),       optional :: &
         &            UG2( UGlo(1): UGhi(1), UGlo(2): UGhi(2), UGlo(3): UGhi(3),NVAR)
    double precision, intent(in),target, optional :: &
         &            U0 ( U0lo(1): U0hi(1), U0lo(2): U0hi(2), U0lo(3): U0hi(3),NVAR)

    integer :: UM1, UM2, UM3  
    double precision :: vflag(3)
    integer :: dd, i, j, k, m, n, elo(3), ehi(3), iflag(3), ii, jj, kk
    double precision :: rho, rhoinv, v(3), ek, e, T, p, c, gamc, dpde, eref, &
         H, gt, cinv, b, gtinv
    double precision, dimension(NSPEC) :: Y, dpdr, d, Yref
    double precision, dimension(NCHARV) :: Udd, vp, vm, vg1, vg2
    double precision, dimension(NCHARV,NCHARV) :: egv, egvt
    double precision :: charv(-2:2,NCHARV) ! characteristic variables
    logical :: do_face, do_gauss
    double precision, pointer :: Ubase(:,:,:,:)

    if (present(U0)) then
       Ubase => U0
    else
       Ubase => U
    end if

    elo = lo
    ehi = hi

    do_face = .false.
    do_gauss = .false.
    
    if (present(UL) .and. present(UR)) then
       do_face = .true.
       elo(idir) = elo(idir)-1
       ehi(idir) = ehi(idir)+1
       UL = 0.d0
       UR = 0.d0
    end if
    
    if (present(UG1) .and. present(UG2)) then
       do_gauss = .true.
       UG1 = 0.d0
       UG2 = 0.d0
    end if

    v = 0.d0
    iflag = 0
    iflag(idir) = 1

    call set_vel_3d(idir, UM1, UM2, UM3, vflag)

    do k      =elo(3), ehi(3)
       do j   =elo(2), ehi(2)
          do i=elo(1), ehi(1)

             rho = 0.d0
             do n=1,nspec
                Y(n) = Ubase(i,j,k,UFS+n-1)
                rho = rho + Y(n)
             end do
       
             rhoInv = 1.d0/rho
       
             Y = Y*rhoInv

             v(1) = Ubase(i,j,k,UM1)*rhoInv*vflag(1)
             v(2) = Ubase(i,j,k,UM2)*rhoInv*vflag(2)
             v(3) = Ubase(i,j,k,UM3)*rhoInv*vflag(3)

             ek = 0.5d0*(v(1)**2+v(2)**2+v(3)**2)
             e  = Ubase(i,j,k,UEDEN)*rhoInV - ek
             T  = Ubase(i,j,k,UTEMP)

             call eos_given_ReY(p,c,gamc,T,dpdr,dpde,rho,e,Y)

             eref = eos_get_eref(Y)
             e = e - eref

             H = e + p*rhoInv + ek

             gt = dpde*rhoInv
             cinv = 1.d0/c
             b = gt*cinv*cinv
             gtinv = 1.d0/gt
             do n=1,nspec
                d(n) = b*(ek - e + dpdr(n)*gtinv)
             end do

             ! assemble left vectors
             egv(1,1) = -0.5d0*(cinv + b*v(1))
             egv(2,1) = -0.5d0*b*v(2)
             egv(3,1) = -0.5d0*b*v(3)
             egv(4,1) =  0.5d0*b
             do n=1,nspec
                egv(CFS+n-1,1) = 0.5d0*(v(1)*cinv + d(n))
             end do

             egv(1,2) =  0.5d0*(cinv - b*v(1))
             egv(2,2) = -0.5d0*b*v(2)
             egv(3,2) = -0.5d0*b*v(3)
             egv(4,2) =  0.5d0*b
             do n=1,nspec
                egv(CFS+n-1,2) = 0.5d0*(-v(1)*cinv + d(n))
             end do

             egv(1,3) = 0.d0
             egv(2,3) = 1.d0
             egv(3,3) = 0.d0
             egv(4,3) = 0.d0
             do n=1,nspec
                egv(CFS+n-1,3) = -v(2)
             end do

             egv(1,4) = 0.d0
             egv(2,4) = 0.d0
             egv(3,4) = 1.d0
             egv(4,4) = 0.d0
             do n=1,nspec
                egv(CFS+n-1,4) = -v(3)
             end do

             do m=1,nspec
                egv(1,CFS+m-1) =  Y(m)*b*v(1)
                egv(2,CFS+m-1) =  Y(m)*b*v(2)
                egv(3,CFS+m-1) =  Y(m)*b*v(3)
                egv(4,CFS+m-1) = -Y(m)*b
                do n=1,nspec
                   egv(CFS+n-1,CFS+m-1) = -Y(m)*d(n)
                end do
                egv(CFS+m-1,CFS+m-1) = egv(CFS+m-1,CFS+m-1) + 1.d0
             end do

             ! convert conserved variables to characteristic variables
             do dd=-2,2

                ii = i+dd*iflag(1)
                jj = j+dd*iflag(2)
                kk = k+dd*iflag(3)

                rhoinV = 1.d0 / U(ii,jj,kk,URHO)

                Udd(1) = U(ii,jj,kk,UM1)
                Udd(2) = U(ii,jj,kk,UM2)*vflag(2)
                Udd(3) = U(ii,jj,kk,UM3)*vflag(3)
                Udd(4) = U(ii,jj,kk,UEDEN)
          
                do n=1,nspec
                   Udd(CFS+n-1) = U(ii,jj,kk,UFS+n-1)
                   Yref(n) = U(ii,jj,kk,UFS+n-1) * rhoInv
                end do
                eref = eos_get_eref(Yref)
                Udd(4) = Udd(4) - U(ii,jj,kk,URHO) * eref

                do n=1,NCHARV
                   charv(dd,n) = dot_product(egv(:,n),Udd)
                end do

             end do

             if (do_face .and. do_gauss) then
                do n=1,NCHARV
                   call weno5(charv(:,n), vp(n), vm(n), vg1(n), vg2(n))
                end do
             else if (do_face) then
                do n=1,NCHARV
                   call weno5(charv(:,n), vp=vp(n), vm=vm(n))
                end do
             else  ! do_gauss
                do n=1,NCHARV
                   call weno5(charv(:,n), vg1=vg1(n), vg2=vg2(n))
                end do
             end if

             ! assemble right vectors
             egv(1,1) = v(1) - c
             egv(2,1) = v(2)
             egv(3,1) = v(3)
             egv(4,1) = H - v(1)*c
             do n=1,nspec
                egv(CFS+n-1,1) = Y(n)
             end do

             egv(1,2) = v(1) + c
             egv(2,2) = v(2)
             egv(3,2) = v(3)
             egv(4,2) = H + v(1)*c
             do n=1,nspec
                egv(CFS+n-1,2) = Y(n)
             end do

             egv(1,3) = 0.d0
             egv(2,3) = 1.d0
             egv(3,3) = 0.d0
             egv(4,3) = v(2)
             do n=1,nspec
                egv(CFS+n-1,3) = 0.d0
             end do
             
             egv(1,4) = 0.d0
             egv(2,4) = 0.d0
             egv(3,4) = 1.d0
             egv(4,4) = v(3)
             do n=1,nspec
                egv(CFS+n-1,4) = 0.d0
             end do

             do n=1,nspec
                egv(1,CFS+n-1) = v(1)
                egv(2,CFS+n-1) = v(2)
                egv(3,CFS+n-1) = v(3)
                egv(4,CFS+n-1) = e + ek - dpdr(n)*gtinv
                do m=1,nspec
                   egv(CFS+m-1,CFS+n-1) = 0.d0
                end do
                egv(CFS+n-1,CFS+n-1) = 1.d0
             end do
             
             egvt = transpose(egv)

             ii = i + iflag(1)
             jj = j + iflag(2)
             kk = k + iflag(3)
             if (do_face .and. ii.ne.hi(1)+2 .and. jj.ne.hi(2)+2 .and. kk.ne.hi(3)+2) then
                UL(ii,jj,kk,UM1)   = UL(ii,jj,kk,UM1)   + dot_product(vp,egvt(:,1))
                UL(ii,jj,kk,UM2)   = UL(ii,jj,kk,UM2)   + dot_product(vp,egvt(:,2))*vflag(2)
                UL(ii,jj,kk,UM3)   = UL(ii,jj,kk,UM3)   + dot_product(vp,egvt(:,3))*vflag(3)
                UL(ii,jj,kk,UEDEN) = UL(ii,jj,kk,UEDEN) + dot_product(vp,egvt(:,4))
                do m=1,nspec
                   UL(ii,jj,kk,UFS+m-1) = UL(ii,jj,kk,UFS+m-1) + dot_product(vp,egvt(:,CFS+m-1))
                   UL(ii,jj,kk,URHO) = UL(ii,jj,kk,URHO) + UL(ii,jj,kk,UFS+m-1)
                end do

                UL(ii,jj,kk,UTEMP) = U(i,j,k,UTEMP)
          
                rhoInv = 1.d0/UL(ii,jj,kk,URHO)
                do n=1,nspec
                   Yref(n) = UL(ii,jj,kk,UFS+n-1) * rhoInv
                end do
                eref = eos_get_eref(Yref)
                UL(ii,jj,kk,UEDEN) = UL(ii,jj,kk,UEDEN) + UL(ii,jj,kk,URHO) * eref
             end if

             if (do_face .and. i.ne.lo(1)-1 .and. j.ne.lo(2)-1 .and. k.ne.lo(3)-1) then
                UR(i,j,k,UM1)   = UR(i,j,k,UM1)   + dot_product(vm,egvt(:,1))
                UR(i,j,k,UM2)   = UR(i,j,k,UM2)   + dot_product(vm,egvt(:,2))*vflag(2)
                UR(i,j,k,UM3)   = UR(i,j,k,UM3)   + dot_product(vm,egvt(:,3))*vflag(3)
                UR(i,j,k,UEDEN) = UR(i,j,k,UEDEN) + dot_product(vm,egvt(:,4))
                do m=1,nspec
                   UR(i,j,k,UFS+m-1) = UR(i,j,k,UFS+m-1) + dot_product(vm, egvt(:,CFS+m-1))
                   UR(i,j,k,URHO) = UR(i,j,k,URHO) + UR(i,j,k,UFS+m-1)
                end do

                UR(i,j,k,UTEMP) = U(i,j,k,UTEMP)

                rhoInv = 1.d0/UR(i,j,k,URHO)
                do n=1,nspec
                   Yref(n) = UR(i,j,k,UFS+n-1) * rhoInv
                end do
                eref = eos_get_eref(Yref)
                UR(i,j,k,UEDEN) = UR(i,j,k,UEDEN) + UR(i,j,k,URHO) * eref
             end if

             if (do_gauss .and.  &
                  i.ne.lo(1)-1 .and. i.ne. hi(1)+1 .and. &
                  j.ne.lo(2)-1 .and. j.ne. hi(2)+1 .and. &
                  k.ne.lo(3)-1 .and. k.ne. hi(3)+1) then

                ! Gauss point 1
                UG1(i,j,k,UM1)   = UG1(i,j,k,UM1)   + dot_product(vg1,egvt(:,1))
                UG1(i,j,k,UM2)   = UG1(i,j,k,UM2)   + dot_product(vg1,egvt(:,2))*vflag(2)
                UG1(i,j,k,UM3)   = UG1(i,j,k,UM3)   + dot_product(vg1,egvt(:,3))*vflag(3)
                UG1(i,j,k,UEDEN) = UG1(i,j,k,UEDEN) + dot_product(vg1,egvt(:,4))
                do m=1,nspec
                   UG1(i,j,k,UFS+m-1) = UG1(i,j,k,UFS+m-1) + dot_product(vg1, egvt(:,CFS+m-1))
                   UG1(i,j,k,URHO) = UG1(i,j,k,URHO) + UG1(i,j,k,UFS+m-1)
                end do
          
                UG1(i,j,k,UTEMP) = U(i,j,k,UTEMP)
       
                rhoInv = 1.d0/UG1(i,j,k,URHO)
                do n=1,nspec
                   Yref(n) = UG1(i,j,k,UFS+n-1) * rhoInv
                end do
                eref = eos_get_eref(Yref)
                UG1(i,j,k,UEDEN) = UG1(i,j,k,UEDEN) + UG1(i,j,k,URHO) * eref
          
                ! Gauss point 2
                UG2(i,j,k,UM1)   = UG2(i,j,k,UM1)   + dot_product(vg2,egvt(:,1))
                UG2(i,j,k,UM2)   = UG2(i,j,k,UM2)   + dot_product(vg2,egvt(:,2))*vflag(2)
                UG2(i,j,k,UM3)   = UG2(i,j,k,UM3)   + dot_product(vg2,egvt(:,3))*vflag(3)
                UG2(i,j,k,UEDEN) = UG2(i,j,k,UEDEN) + dot_product(vg2,egvt(:,4))
                do m=1,nspec
                   UG2(i,j,k,UFS+m-1) = UG2(i,j,k,UFS+m-1) + dot_product(vg2, egvt(:,CFS+m-1))
                   UG2(i,j,k,URHO) = UG2(i,j,k,URHO) + UG2(i,j,k,UFS+m-1)
                end do
          
                UG2(i,j,k,UTEMP) = U(i,j,k,UTEMP)
          
                rhoInv = 1.d0/UG2(i,j,k,URHO)
                do n=1,nspec
                   Yref(n) = UG2(i,j,k,UFS+n-1) * rhoInv
                end do
                eref = eos_get_eref(Yref)
                UG2(i,j,k,UEDEN) = UG2(i,j,k,UEDEN) + UG2(i,j,k,URHO) * eref

             end if
             
          end do
       end do
    end do

    Nullify(Ubase)

  end subroutine reconstruct_char_3d


  subroutine reconstruct(lo, hi, &
       Ulo, Uhi, &
       ULRlo, ULRhi, &
       UGlo, UGhi, &
       U0lo, U0hi, &
       U, UL, UR, UG1, UG2, U0, dir)

    integer, intent(in) :: lo, hi, Ulo, Uhi, ULRlo, ULRhi, UGlo, UGhi, U0lo, U0hi 
    integer, intent(in), optional :: dir
    double precision, intent(in),target           :: U ( Ulo: Uhi,NVAR)
    double precision, intent(in),target, optional :: U0(U0lo:U0hi,NVAR)
    double precision, intent(out), dimension(ULRlo:ULRhi,NVAR), optional :: UL, UR
    double precision, intent(out), dimension( UGlo: UGhi,NVAR), optional :: UG1, UG2

    if (do_component_weno) then
       call reconstruct_comp(lo, hi, &
            Ulo, Uhi, &
            ULRlo, ULRhi, &
            UGlo, UGhi, &
            U, UL, UR, UG1, UG2)
    else
       call reconstruct_char(lo, hi, &
            Ulo, Uhi, &
            ULRlo, ULRhi, &
            UGlo, UGhi, &
            U0lo, U0hi, &
            U, UL, UR, UG1, UG2, U0, dir)
    end if

  end subroutine reconstruct

  subroutine reconstruct_comp(lo, hi, &
       Ulo, Uhi, &
       ULRlo, ULRhi, &
       UGlo, UGhi, &
       U, UL, UR, UG1, UG2)

    integer, intent(in) :: lo, hi, Ulo, Uhi, ULRlo, ULRhi, UGlo, UGhi
    double precision, intent(in) :: U(Ulo:Uhi,NVAR)
    double precision, intent(out), dimension(ULRlo:ULRhi,NVAR), optional :: UL, UR
    double precision, intent(out), dimension( UGlo: UGhi,NVAR), optional :: UG1, UG2

    integer :: i, n, iextra
    double precision, dimension(lo-1:hi+1) :: vp, vm
    double precision, dimension(lo  :hi  ) :: vg1, vg2
    double precision :: rhoE(Ulo:Uhi), rho, rhoInv, Y(nspec)
    logical :: do_gauss, do_face

    iextra = 0
    do_face = .false.
    do_gauss = .false.

    if (present(UL) .and. present(UR)) then
       do_face = .true.
       iextra  = 1
    end if

    if (present(UG1) .and. present(UG2)) then
       do_gauss = .true.
    end if

    do i = lo-iextra-2, hi+iextra+2

       rho = 0.d0
       do n=1,nspec
          Y(n) = U(i,UFS+n-1)
          rho = rho + Y(n)
       end do
       
       rhoInv = 1.d0/rho
       
       do n=1,nspec
          Y(n) = Y(n) * rhoInv
       end do

       rhoE(i) = U(i,UEDEN) - rho*eos_get_eref(Y)

    end do

    if (do_face .and. do_gauss) then

       do n=1,NVAR
          if (n.eq.UEDEN) then
             call vweno5(lo-1,hi+1, rhoE,Ulo,Uhi, lo,hi, vp,vm, vg1,vg2)
          else
             call vweno5(lo-1,hi+1, U(:,n),Ulo,Uhi, lo,hi, vp,vm,vg1,vg2)
          end if

          do i=lo, hi+1
             UL(i,n) = vp(i-1)
             UR(i,n) = vm(i)
          end do

          do i=lo, hi
             UG1(i,n) = vg1(i)
             UG2(i,n) = vg2(i)
          end do
       end do

    else if (do_face) then

       do n=1,NVAR
          if (n.eq.UEDEN) then
             call vweno5(lo-1,hi+1, rhoE,Ulo,Uhi, lo,hi, vp=vp,vm=vm)
          else
             call vweno5(lo-1,hi+1, U(:,n),Ulo,Uhi, lo,hi, vp=vp,vm=vm)
          end if

          do i=lo, hi+1
             UL(i,n) = vp(i-1)
             UR(i,n) = vm(i)
          end do
       end do

    else ! do_gauss

       do n=1,NVAR
          if (n.eq.UEDEN) then
             call vweno5(lo,hi, rhoE,Ulo,Uhi, lo,hi, vg1=vg1,vg2=vg2)
          else
             call vweno5(lo,hi, U(:,n),Ulo,Uhi, lo,hi, vg1=vg1,vg2=vg2)
          end if

          do i=lo, hi
             UG1(i,n) = vg1(i)
             UG2(i,n) = vg2(i)
          end do
       end do

    end if

    if (do_face) then
       call normalize(lo,hi+1, UL, ULRlo,ULRhi)
       call normalize(lo,hi+1, UR, ULRlo,ULRhi)
    end if

    if (do_gauss) then
       call normalize(lo,hi, UG1, UGlo,UGhi)
       call normalize(lo,hi, UG2, UGlo,UGhi)
    end if

  end subroutine reconstruct_comp
    

  subroutine reconstruct_char(lo, hi, &
       Ulo, Uhi, &
       ULRlo, ULRhi, &
       UGlo, UGhi, &
       U0lo, U0hi, &
       U, UL, UR, UG1, UG2, U0, dir)

    integer, intent(in) :: lo, hi, Ulo, Uhi, ULRlo, ULRhi, UGlo, UGhi, U0lo, U0hi 
    integer, intent(in), optional :: dir
    double precision, intent(in),target           :: U ( Ulo: Uhi,NVAR)
    double precision, intent(in),target, optional :: U0(U0lo:U0hi,NVAR)
    double precision, intent(out), dimension(ULRlo:ULRhi,NVAR), optional :: UL, UR
    double precision, intent(out), dimension( UGlo: UGhi,NVAR), optional :: UG1, UG2

    integer :: i, ii, ivar, m, n, ivel(3), idir, iextra
    double precision :: egv(NCHARV,NCHARV), egvt(NCHARV,NCHARV)
    double precision :: gt, b, d(NSPEC)
    double precision :: rho, rhoInv, p, c, gamc, T, dpdr(NSPEC), dpde, e, ek, H, Y(NSPEC)
    double precision :: gtinv, cinv
    double precision :: vel(3), vflag(3)
    double precision :: charv(-2:2,NCHARV) ! characteristic variables
    double precision, dimension(NCHARV) :: vp, vm, vg1, vg2
    double precision :: eref, Yref(NSPEC), Uii(NCHARV)
    logical :: do_gauss, do_face
    double precision, pointer :: Ubase(:,:)

    if (present(U0)) then
       Ubase => U0
    else
       Ubase => U
    end if

    if (present(dir)) then
       idir = dir
    else
       idir = 1
    end if

    iextra = 0
    do_face = .false.
    do_gauss = .false.

    if (present(UL) .and. present(UR)) then
       do_face = .true.
       iextra  = 1
       do n=1,NVAR
          do i=lo, hi+1
             UL(i,n) = 0.d0
             UR(i,n) = 0.d0
          end do
       end do
    end if

    if (present(UG1) .and. present(UG2)) then
       do_gauss = .true.
       do n=1,NVAR
          do i=lo, hi
             UG1(i,n) = 0.d0
             UG2(i,n) = 0.d0
          end do
       end do
    end if

    call set_vel(idir, ivel, vflag)


    do i = lo-iextra, hi+iextra
       
       rho = 0.d0
       do n=1,nspec
          Y(n) = Ubase(i,UFS+n-1)
          rho = rho + Y(n)
       end do
       
       rhoInv = 1.d0/rho
       
       do n=1,nspec
          Y(n) = Y(n) * rhoInv
       end do
       
       vel(1) = vflag(1) * Ubase(i,ivel(1))*rhoInv
       vel(2) = vflag(2) * Ubase(i,ivel(2))*rhoInv
       vel(3) = vflag(3) * Ubase(i,ivel(3))*rhoInv

       ek = 0.5d0*(vel(1)**2 + vel(2)**2 + vel(3)**2)
       e = Ubase(i,UEDEN)*rhoInV - ek
       T = Ubase(i,UTEMP)

       call eos_given_ReY(p,c,gamc,T,dpdr,dpde,rho,e,Y)

       eref = eos_get_eref(Y)
       e = e - eref

       H = e + p*rhoInv + ek

       gt = dpde*rhoInv
       cinv = 1.d0/c
       b = gt*cinv*cinv
       gtinv = 1.d0/gt
       do n=1,nspec
          d(n) = b*(ek - e + dpdr(n)*gtinv)
       end do

       ! assemble left vectors
       egv(1,1) = -0.5d0*(cinv + b*vel(1))
       egv(2,1) = -0.5d0*b*vel(2)
       egv(3,1) = -0.5d0*b*vel(3)
       egv(4,1) =  0.5d0*b
       do n=1,nspec
          egv(CFS+n-1,1) = 0.5d0*(vel(1)*cinv + d(n))
       end do

       egv(1,2) =  0.5d0*(cinv - b*vel(1))
       egv(2,2) = -0.5d0*b*vel(2)
       egv(3,2) = -0.5d0*b*vel(3)
       egv(4,2) =  0.5d0*b
       do n=1,nspec
          egv(CFS+n-1,2) = 0.5d0*(-vel(1)*cinv + d(n))
       end do

       egv(1,3) = 0.d0
       egv(2,3) = 1.d0
       egv(3,3) = 0.d0
       egv(4,3) = 0.d0
       do n=1,nspec
          egv(CFS+n-1,3) = -vel(2)
       end do

       egv(1,4) = 0.d0
       egv(2,4) = 0.d0
       egv(3,4) = 1.d0
       egv(4,4) = 0.d0
       do n=1,nspec
          egv(CFS+n-1,4) = -vel(3)
       end do

       do m=1,nspec
          egv(1,CFS+m-1) =  Y(m)*b*vel(1)
          egv(2,CFS+m-1) =  Y(m)*b*vel(2)
          egv(3,CFS+m-1) =  Y(m)*b*vel(3)
          egv(4,CFS+m-1) = -Y(m)*b
          do n=1,nspec
             egv(CFS+n-1,CFS+m-1) = -Y(m)*d(n)
          end do
          egv(CFS+m-1,CFS+m-1) = egv(CFS+m-1,CFS+m-1) + 1.d0
       end do

       ! convert conserved variables to characteristic variables
       do ii=-2,2

          rhoinV = 1.d0 / U(i+ii,URHO)

          Uii(1) = U(i+ii,ivel(1)) * vflag(1)
          Uii(2) = U(i+ii,ivel(2)) * vflag(2)
          Uii(3) = U(i+ii,ivel(3)) * vflag(3)
          Uii(4) = U(i+ii,UEDEN)
          
          do n=1,nspec
             Uii(CFS+n-1) = U(i+ii,UFS+n-1)
             Yref(n) = U(i+ii,UFS+n-1) * rhoInv
          end do
          eref = eos_get_eref(Yref)
          Uii(4) = Uii(4) - U(i+ii,URHO) * eref

          do n=1,NCHARV
             charv(ii,n) = dot_product(egv(:,n),Uii)
          end do

       end do

       if (do_face .and. do_gauss) then
          do ivar=1,NCHARV
             call weno5(charv(:,ivar), vp(ivar), vm(ivar), vg1(ivar), vg2(ivar))
          end do
       else if (do_face) then
          do ivar=1,NCHARV
             call weno5(charv(:,ivar), vp=vp(ivar), vm=vm(ivar))
          end do
       else  ! do_gauss
          do ivar=1,NCHARV
             call weno5(charv(:,ivar), vg1=vg1(ivar), vg2=vg2(ivar))
          end do
       end if

       ! assemble right vectors
       egv(1,1) = vel(1) - c
       egv(2,1) = vel(2)
       egv(3,1) = vel(3)
       egv(4,1) = H - vel(1)*c
       do n=1,nspec
          egv(CFS+n-1,1) = Y(n)
       end do

       egv(1,2) = vel(1) + c
       egv(2,2) = vel(2)
       egv(3,2) = vel(3)
       egv(4,2) = H + vel(1)*c
       do n=1,nspec
          egv(CFS+n-1,2) = Y(n)
       end do

       egv(1,3) = 0.d0
       egv(2,3) = 1.d0
       egv(3,3) = 0.d0
       egv(4,3) = vel(2)
       do n=1,nspec
          egv(CFS+n-1,3) = 0.d0
       end do

       egv(1,4) = 0.d0
       egv(2,4) = 0.d0
       egv(3,4) = 1.d0
       egv(4,4) = vel(3)
       do n=1,nspec
          egv(CFS+n-1,4) = 0.d0
       end do

       do n=1,nspec
          egv(1,CFS+n-1) = vel(1)
          egv(2,CFS+n-1) = vel(2)
          egv(3,CFS+n-1) = vel(3)
          egv(4,CFS+n-1) = e + ek - dpdr(n)*gtinv
          do m=1,nspec
             egv(CFS+m-1,CFS+n-1) = 0.d0
          end do
          egv(CFS+n-1,CFS+n-1) = 1.d0
       end do

       egvt = transpose(egv)

       if (do_face .and. i.ne.hi+1) then
          do n=1,NCHARV
             UL(i+1,ivel(1)) = UL(i+1,ivel(1)) + vp(n)*egv(1,n)*vflag(1)
             UL(i+1,ivel(2)) = UL(i+1,ivel(2)) + vp(n)*egv(2,n)*vflag(2)
             UL(i+1,ivel(3)) = UL(i+1,ivel(3)) + vp(n)*egv(3,n)*vflag(3)
             UL(i+1,UEDEN  ) = UL(i+1,UEDEN  ) + vp(n)*egv(4,n)
          end do

          do m=1,nspec
             UL(i+1,UFS+m-1) = UL(i+1,UFS+m-1) + dot_product(vp, egvt(:,CFS+m-1))
             UL(i+1,URHO) = UL(i+1,URHO) + UL(i+1,UFS+m-1)
          end do

          UL(i+1,UTEMP) = U(i,UTEMP)
          
          rhoInv = 1.d0/UL(i+1,URHO)
          do n=1,nspec
             Yref(n) = UL(i+1,UFS+n-1) * rhoInv
          end do
          eref = eos_get_eref(Yref)
          UL(i+1,UEDEN) = UL(i+1,UEDEN) + UL(i+1,URHO) * eref

       end if

       if (do_face .and. i.ne.lo-1) then
          do n=1,NCHARV
             UR(i,ivel(1)) = UR(i,ivel(1)) + vm(n)*egv(1,n)*vflag(1)
             UR(i,ivel(2)) = UR(i,ivel(2)) + vm(n)*egv(2,n)*vflag(2)
             UR(i,ivel(3)) = UR(i,ivel(3)) + vm(n)*egv(3,n)*vflag(3)
             UR(i,UEDEN  ) = UR(i,UEDEN  ) + vm(n)*egv(4,n)
          end do

          do m=1,nspec
             UR(i,UFS+m-1) = UR(i,UFS+m-1) + dot_product(vm, egvt(:,CFS+m-1))
             UR(i,URHO) = UR(i,URHO) + UR(i,UFS+m-1)
          end do

          UR(i,UTEMP) = U(i,UTEMP)

          rhoInv = 1.d0/UR(i,URHO)
          do n=1,nspec
             Yref(n) = UR(i,UFS+n-1) * rhoInv
          end do
          eref = eos_get_eref(Yref)
          UR(i,UEDEN) = UR(i,UEDEN) + UR(i,URHO) * eref
       end if

       if (do_gauss .and. i.ne.lo-1 .and. i.ne. hi+1) then

          ! Gauss point 1
          do n=1,NCHARV
             UG1(i,ivel(1)) = UG1(i,ivel(1)) + vg1(n)*egv(1,n)*vflag(1)
             UG1(i,ivel(2)) = UG1(i,ivel(2)) + vg1(n)*egv(2,n)*vflag(2)
             UG1(i,ivel(3)) = UG1(i,ivel(3)) + vg1(n)*egv(3,n)*vflag(3)
             UG1(i,UEDEN  ) = UG1(i,UEDEN  ) + vg1(n)*egv(4,n)
          end do

          do m=1,nspec
             UG1(i,UFS+m-1) = UG1(i,UFS+m-1) + dot_product(vg1, egvt(:,CFS+m-1))
             UG1(i,URHO) = UG1(i,URHO) + UG1(i,UFS+m-1)
          end do
          
          UG1(i,UTEMP) = U(i,UTEMP)
       
          rhoInv = 1.d0/UG1(i,URHO)
          do n=1,nspec
             Yref(n) = UG1(i,UFS+n-1) * rhoInv
          end do
          eref = eos_get_eref(Yref)
          UG1(i,UEDEN) = UG1(i,UEDEN) + UG1(i,URHO) * eref
          
          ! Gauss point 2
          do n=1,NCHARV
             UG2(i,ivel(1)) = UG2(i,ivel(1)) + vg2(n)*egv(1,n)*vflag(1)
             UG2(i,ivel(2)) = UG2(i,ivel(2)) + vg2(n)*egv(2,n)*vflag(2)
             UG2(i,ivel(3)) = UG2(i,ivel(3)) + vg2(n)*egv(3,n)*vflag(3)
             UG2(i,UEDEN  ) = UG2(i,UEDEN  ) + vg2(n)*egv(4,n)
          end do

          do m=1,nspec
             UG2(i,UFS+m-1) = UG2(i,UFS+m-1) + dot_product(vg2, egvt(:,CFS+m-1))
             UG2(i,URHO) = UG2(i,URHO) + UG2(i,UFS+m-1)
          end do
          
          UG2(i,UTEMP) = U(i,UTEMP)
          
          rhoInv = 1.d0/UG2(i,URHO)
          do n=1,nspec
             Yref(n) = UG2(i,UFS+n-1) * rhoInv
          end do
          eref = eos_get_eref(Yref)
          UG2(i,UEDEN) = UG2(i,UEDEN) + UG2(i,URHO) * eref

       end if

    end do

    Nullify(Ubase)

  end subroutine reconstruct_char


  subroutine set_vel(idir, ivel, vflag)
    integer, intent(in) :: idir
    integer, intent(out) :: ivel(3)
    double precision, intent(out) :: vflag(3)
    if (ndim .eq. 1) then
       vflag(1) = 1.d0
       vflag(2) = 0.d0
       vflag(3) = 0.d0
       ivel(1) = UMX
       ivel(2) = UMX
       ivel(3) = UMX
    else if (ndim .eq. 2) then
       vflag(1) = 1.d0
       vflag(2) = 1.d0
       vflag(3) = 0.d0   
       if (idir .eq. 1) then
          ivel(1) = UMX
          ivel(2) = UMY
          ivel(3) = UMX
       else
          ivel(1) = UMY
          ivel(2) = UMX
          ivel(3) = UMX
       end if
    else
       vflag(1) = 1.d0
       vflag(2) = 1.d0
       vflag(3) = 1.d0 
       if (idir .eq. 1) then
          ivel(1) = UMX
          ivel(2) = UMY
          ivel(3) = UMZ
       else if (idir .eq. 2) then
          ivel(1) = UMY
          ivel(2) = UMZ
          ivel(3) = UMX
       else
          ivel(1) = UMZ
          ivel(2) = UMX
          ivel(3) = UMY
       end if
    end if
  end subroutine set_vel


  subroutine set_vel_3d(idir, UM1, UM2, UM3, vflag)
    integer, intent(in) :: idir
    integer, intent(out) :: UM1, UM2, UM3
    double precision, intent(out) :: vflag(3)
    if (ndim .eq. 1) then
       vflag(1) = 1.d0
       vflag(2) = 0.d0
       vflag(3) = 0.d0
       UM1 = UMX
       UM2 = UMX
       UM3 = UMX
    else if (ndim .eq. 2) then
       vflag(1) = 1.d0
       vflag(2) = 1.d0
       vflag(3) = 0.d0   
       if (idir .eq. 1) then
          UM1 = UMX
          UM2 = UMY
          UM3 = UMX
       else
          UM1 = UMY
          UM2 = UMX
          UM3 = UMX
       end if
    else
       vflag(1) = 1.d0
       vflag(2) = 1.d0
       vflag(3) = 1.d0 
       if (idir .eq. 1) then
          UM1 = UMX
          UM2 = UMY
          UM3 = UMZ
       else if (idir .eq. 2) then
          UM1 = UMY
          UM2 = UMZ
          UM3 = UMX
       else
          UM1 = UMZ
          UM2 = UMX
          UM3 = UMY
       end if
    end if
  end subroutine set_vel_3d


  subroutine normalize(lo,hi,U,Ulo,Uhi)
    integer, intent(in) :: lo, hi, Ulo, Uhi
    double precision, intent(inout) :: U(Ulo:Uhi,NVAR)

    integer :: i, n
    double precision :: rhoInv, Y(nspec)
    
    do i=lo,hi
       U(i,URHO) = 0.d0
       do n=1, nspec
          Y(n) = U(i,UFS+n-1)
          U(i,URHO) = U(i,URHO) + Y(n)
       end do
       rhoInv = 1.d0/U(i,URHO)

       do n=1, nspec
          Y(n) = Y(n) * rhoInv
       end do
       
       U(i,UEDEN) = U(i,UEDEN) + eos_get_eref(Y)*U(i,URHO)
    end do
  end subroutine normalize


  subroutine normalize_3d(lo,hi,U,Ulo,Uhi)
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),Ulo(3):Uhi(3),NVAR)

    integer :: i, j, k, n
    double precision :: rhoInv, Y(nspec)
    
    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       U(i,j,k,URHO) = 0.d0
       do n=1, nspec
          Y(n) = U(i,j,k,UFS+n-1)
          U(i,j,k,URHO) = U(i,j,k,URHO) + Y(n)
       end do
       rhoInv = 1.d0/U(i,j,k,URHO)

       do n=1, nspec
          Y(n) = Y(n) * rhoInv
       end do
       
       U(i,j,k,UEDEN) = U(i,j,k,UEDEN) + eos_get_eref(Y)*U(i,j,k,URHO)
    end do
    end do
    end do
  end subroutine normalize_3d

end module reconstruct_module
