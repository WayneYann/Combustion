module reconstruct_module

  use meth_params_module, only : ndim, NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, &
       UFS, NSPEC, NCHARV, CFS, do_component_weno, do_mp5
  use weno_module, only : mp5, weno5, vweno5, weno5_center
  use eos_module, only : eos_given_RTY, eos_given_ReY, eos_get_eref
  use renorm_module, only : floor_species
  use passinfo_module, only : level

  implicit none

  private

  public :: reconstruct, reconstruct_comp, reconstruct_center, get_eigen_matrices_q

contains

  ! L and R in UL and UR are relative to face
  ! UG1 and UG2 are at two Gauss points
  !
  ! Minimal ranges:
  !   U         : lo-2:hi+2 if UL & UR are not present; 
  !               lo-3:hi+3 if UL & UR are present
  !   UL  & UR  : lo  :hi+1
  !   UG1 & UG2 : lo  :hi
  !   U0        : lo  :hi   if UL & UR are not present; 
  !               lo-1:hi+1 if UL & UR are present
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
    double precision :: rhoInv, vflag(3)
    double precision :: egv1(NCHARV,NCHARV), egv2(NCHARV,NCHARV)
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

       ! egv1: left matrix;  egv2: right matrix
       call get_eigen_matrices(Ubase(i,:), egv1, egv2, ivel, vflag)

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
             charv(ii,n) = dot_product(egv1(:,n),Uii)
          end do
       end do

       if (do_face .and. do_gauss) then
          do ivar=1,NCHARV
             call weno5(charv(:,ivar), vp(ivar), vm(ivar), vg1(ivar), vg2(ivar))
          end do
       else if (do_face) then
          if (do_mp5) then
             do ivar=1,NCHARV
                call mp5(charv(:,ivar), vp(ivar), vm(ivar))
             end do
          else
             do ivar=1,NCHARV
                call weno5(charv(:,ivar), vp=vp(ivar), vm=vm(ivar))
             end do
          end if
       else  ! do_gauss
          do ivar=1,NCHARV
             call weno5(charv(:,ivar), vg1=vg1(ivar), vg2=vg2(ivar))
          end do
       end if

       egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

       if (do_face .and. i.ne.hi+1) then
          do n=1,NCHARV
             UL(i+1,ivel(1)) = UL(i+1,ivel(1)) + vp(n)*egv2(1,n)*vflag(1)
             UL(i+1,ivel(2)) = UL(i+1,ivel(2)) + vp(n)*egv2(2,n)*vflag(2)
             UL(i+1,ivel(3)) = UL(i+1,ivel(3)) + vp(n)*egv2(3,n)*vflag(3)
             UL(i+1,UEDEN  ) = UL(i+1,UEDEN  ) + vp(n)*egv2(4,n)
          end do

          do m=1,nspec
             UL(i+1,UFS+m-1) = UL(i+1,UFS+m-1) + dot_product(vp, egv1(:,CFS+m-1))
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
             UR(i,ivel(1)) = UR(i,ivel(1)) + vm(n)*egv2(1,n)*vflag(1)
             UR(i,ivel(2)) = UR(i,ivel(2)) + vm(n)*egv2(2,n)*vflag(2)
             UR(i,ivel(3)) = UR(i,ivel(3)) + vm(n)*egv2(3,n)*vflag(3)
             UR(i,UEDEN  ) = UR(i,UEDEN  ) + vm(n)*egv2(4,n)
          end do

          do m=1,nspec
             UR(i,UFS+m-1) = UR(i,UFS+m-1) + dot_product(vm, egv1(:,CFS+m-1))
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
             UG1(i,ivel(1)) = UG1(i,ivel(1)) + vg1(n)*egv2(1,n)*vflag(1)
             UG1(i,ivel(2)) = UG1(i,ivel(2)) + vg1(n)*egv2(2,n)*vflag(2)
             UG1(i,ivel(3)) = UG1(i,ivel(3)) + vg1(n)*egv2(3,n)*vflag(3)
             UG1(i,UEDEN  ) = UG1(i,UEDEN  ) + vg1(n)*egv2(4,n)
          end do

          do m=1,nspec
             UG1(i,UFS+m-1) = UG1(i,UFS+m-1) + dot_product(vg1, egv1(:,CFS+m-1))
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
             UG2(i,ivel(1)) = UG2(i,ivel(1)) + vg2(n)*egv2(1,n)*vflag(1)
             UG2(i,ivel(2)) = UG2(i,ivel(2)) + vg2(n)*egv2(2,n)*vflag(2)
             UG2(i,ivel(3)) = UG2(i,ivel(3)) + vg2(n)*egv2(3,n)*vflag(3)
             UG2(i,UEDEN  ) = UG2(i,UEDEN  ) + vg2(n)*egv2(4,n)
          end do

          do m=1,nspec
             UG2(i,UFS+m-1) = UG2(i,UFS+m-1) + dot_product(vg2, egv1(:,CFS+m-1))
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


  subroutine reconstruct_center(lo,hi,U,Ulo,Uhi,Uc,clo,chi,idir)
    integer, intent(in) :: lo,hi,Ulo,Uhi,clo,chi,idir
    double precision, intent(in) :: U(Ulo:Uhi,NVAR)
    double precision, intent(out) :: Uc(clo:chi,NVAR)
    
    integer :: i, ii, ivar, m, n, ivel(3)
    double precision :: rhoInv, vflag(3)
    double precision :: egv1(NCHARV,NCHARV), egv2(NCHARV,NCHARV)
    double precision :: charv(-2:2,NCHARV) ! characteristic variables
    double precision, dimension(NCHARV) :: vc
    double precision :: eref, Yref(NSPEC), Uii(NCHARV)

    do n=1,NVAR
       do i=lo,hi
          Uc(i,n) = 0.d0
       end do
    end do

    call set_vel(idir, ivel, vflag)

    do i=lo,hi

       ! egv1: left matrix;  egv2: right matrix
       call get_eigen_matrices(U(i,:), egv1, egv2, ivel, vflag)

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
             charv(ii,n) = dot_product(egv1(:,n),Uii)
          end do
       end do

       do ivar=1,NCHARV
          call weno5_center(charv(:,ivar), vc(ivar))
       end do

       egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

       do n=1,NCHARV
          Uc(i,ivel(1)) = Uc(i,ivel(1)) + vc(n)*egv2(1,n)*vflag(1)
          Uc(i,ivel(2)) = Uc(i,ivel(2)) + vc(n)*egv2(2,n)*vflag(2)
          Uc(i,ivel(3)) = Uc(i,ivel(3)) + vc(n)*egv2(3,n)*vflag(3)
          Uc(i,UEDEN  ) = Uc(i,UEDEN  ) + vc(n)*egv2(4,n)
       end do
       
       do m=1,nspec
          Uc(i,UFS+m-1) = Uc(i,UFS+m-1) + dot_product(vc, egv1(:,CFS+m-1))
          Uc(i,URHO) = Uc(i,URHO) + Uc(i,UFS+m-1)
       end do
       
       Uc(i,UTEMP) = U(i,UTEMP)
       
       rhoInv = 1.d0/Uc(i,URHO)
       do n=1,nspec
          Yref(n) = Uc(i,UFS+n-1) * rhoInv
       end do
       eref = eos_get_eref(Yref)
       Uc(i,UEDEN) = Uc(i,UEDEN) + Uc(i,URHO) * eref
    end do

  end subroutine reconstruct_center


  subroutine reconstruct_center_comp(lo,hi,U,Ulo,Uhi,Uc,clo,chi,idir)
    integer, intent(in) :: lo,hi,Ulo,Uhi,clo,chi,idir
    double precision, intent(in) :: U(Ulo:Uhi,NVAR)
    double precision, intent(out) :: Uc(clo:chi,NVAR)
    
    integer :: i, n
    double precision :: rhoE(Ulo:Uhi), rho, rhoInv, Y(nspec)

    do i=lo,hi
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

    do i=lo,hi
       do n=1,NVAR
          if (n.eq.UEDEN) then
             call weno5_center(rhoE(i-2:i+2), Uc(i,n))
          else
             call weno5_center(U(i-2:i+2,n), Uc(i,n))
          end if
       end do
    end do

    call normalize(lo,hi,Uc,clo,chi)

  end subroutine reconstruct_center_comp


  subroutine get_eigen_matrices(U,lv,rv,ivel,vflag)
    double precision, intent(in) :: U(NVAR), vflag(3)
    integer, intent(in) :: ivel(3)
    double precision, intent(out) :: lv(NCHARV,NCHARV), rv(NCHARV,NCHARV)
    
    integer :: m, n, ierr
    double precision :: rho, rhoInv, p, c, gamc, T, dpdr(NSPEC), dpde, e, ek, H, Y(NSPEC)
    double precision :: gt, b, d(NSPEC), gtinv, cinv, vel(3)
    double precision :: eref

    rho = 0.d0
    do n=1,nspec
       Y(n) = U(UFS+n-1)
       rho = rho + Y(n)
    end do
    
    rhoInv = 1.d0/rho
    
    do n=1,nspec
       Y(n) = Y(n) * rhoInv
    end do
    
    vel(1) = vflag(1) * U(ivel(1))*rhoInv
    vel(2) = vflag(2) * U(ivel(2))*rhoInv
    vel(3) = vflag(3) * U(ivel(3))*rhoInv

    ek = 0.5d0*(vel(1)**2 + vel(2)**2 + vel(3)**2)
    e = U(UEDEN)*rhoInV - ek
    T = U(UTEMP)

    call floor_species(nspec, Y)
    
    call eos_given_ReY(p,c,gamc,T,dpdr,dpde,rho,e,Y, ierr=ierr)
    if (ierr .ne. 0) then
       print *, 'get_eigen_matrices: eos_given_ReY failed at ', &
            level, rho, e, Y, T, U(UTEMP)
       call bl_error('get_eigen_matrices: eos_given_ReY failed')
    end if
    
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
    lv(1,1) = -0.5d0*(cinv + b*vel(1))
    lv(2,1) = -0.5d0*b*vel(2)
    lv(3,1) = -0.5d0*b*vel(3)
    lv(4,1) =  0.5d0*b
    do n=1,nspec
       lv(CFS+n-1,1) = 0.5d0*(vel(1)*cinv + d(n))
    end do
    
    lv(1,2) =  0.5d0*(cinv - b*vel(1))
    lv(2,2) = -0.5d0*b*vel(2)
    lv(3,2) = -0.5d0*b*vel(3)
    lv(4,2) =  0.5d0*b
    do n=1,nspec
       lv(CFS+n-1,2) = 0.5d0*(-vel(1)*cinv + d(n))
    end do
    
    lv(1,3) = 0.d0
    lv(2,3) = 1.d0
    lv(3,3) = 0.d0
    lv(4,3) = 0.d0
    do n=1,nspec
       lv(CFS+n-1,3) = -vel(2)
    end do
    
    lv(1,4) = 0.d0
    lv(2,4) = 0.d0
    lv(3,4) = 1.d0
    lv(4,4) = 0.d0
    do n=1,nspec
       lv(CFS+n-1,4) = -vel(3)
    end do
    
    do m=1,nspec
       lv(1,CFS+m-1) =  Y(m)*b*vel(1)
       lv(2,CFS+m-1) =  Y(m)*b*vel(2)
       lv(3,CFS+m-1) =  Y(m)*b*vel(3)
       lv(4,CFS+m-1) = -Y(m)*b
       do n=1,nspec
          lv(CFS+n-1,CFS+m-1) = -Y(m)*d(n)
       end do
       lv(CFS+m-1,CFS+m-1) = lv(CFS+m-1,CFS+m-1) + 1.d0
    end do
    
    ! assemble right vectors
    rv(1,1) = vel(1) - c
    rv(2,1) = vel(2)
    rv(3,1) = vel(3)
    rv(4,1) = H - vel(1)*c
    do n=1,nspec
       rv(CFS+n-1,1) = Y(n)
    end do
    
    rv(1,2) = vel(1) + c
    rv(2,2) = vel(2)
    rv(3,2) = vel(3)
    rv(4,2) = H + vel(1)*c
    do n=1,nspec
       rv(CFS+n-1,2) = Y(n)
    end do
    
    rv(1,3) = 0.d0
    rv(2,3) = 1.d0
    rv(3,3) = 0.d0
    rv(4,3) = vel(2)
    do n=1,nspec
       rv(CFS+n-1,3) = 0.d0
    end do
    
    rv(1,4) = 0.d0
    rv(2,4) = 0.d0
    rv(3,4) = 1.d0
    rv(4,4) = vel(3)
    do n=1,nspec
       rv(CFS+n-1,4) = 0.d0
    end do
    
    do n=1,nspec
       rv(1,CFS+n-1) = vel(1)
       rv(2,CFS+n-1) = vel(2)
       rv(3,CFS+n-1) = vel(3)
       rv(4,CFS+n-1) = e + ek - dpdr(n)*gtinv
       do m=1,nspec
          rv(CFS+m-1,CFS+n-1) = 0.d0
       end do
       rv(CFS+n-1,CFS+n-1) = 1.d0
    end do

  end subroutine get_eigen_matrices


  subroutine get_eigen_matrices_q(rho, Y, T, vel, lv,rv)
    double precision, intent(in) :: rho, Y(NSPEC), T, vel(3)
    double precision, intent(out) :: lv(NCHARV,NCHARV), rv(NCHARV,NCHARV)
    
    integer :: m, n
    double precision :: rhoInv, p, c, dpdr(NSPEC), dpde, e, ek, H
    double precision :: gt, b, d(NSPEC), gtinv, cinv
    double precision :: eref

    rhoInv = 1.d0/rho
    
    ek = 0.5d0*(vel(1)**2 + vel(2)**2 + vel(3)**2)

    call eos_given_RTY(e, p, c, dpdr, dpde, rho, T, Y)
    
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
    lv(1,1) = -0.5d0*(cinv + b*vel(1))
    lv(2,1) = -0.5d0*b*vel(2)
    lv(3,1) = -0.5d0*b*vel(3)
    lv(4,1) =  0.5d0*b
    do n=1,nspec
       lv(CFS+n-1,1) = 0.5d0*(vel(1)*cinv + d(n))
    end do
    
    lv(1,2) =  0.5d0*(cinv - b*vel(1))
    lv(2,2) = -0.5d0*b*vel(2)
    lv(3,2) = -0.5d0*b*vel(3)
    lv(4,2) =  0.5d0*b
    do n=1,nspec
       lv(CFS+n-1,2) = 0.5d0*(-vel(1)*cinv + d(n))
    end do
    
    lv(1,3) = 0.d0
    lv(2,3) = 1.d0
    lv(3,3) = 0.d0
    lv(4,3) = 0.d0
    do n=1,nspec
       lv(CFS+n-1,3) = -vel(2)
    end do
    
    lv(1,4) = 0.d0
    lv(2,4) = 0.d0
    lv(3,4) = 1.d0
    lv(4,4) = 0.d0
    do n=1,nspec
       lv(CFS+n-1,4) = -vel(3)
    end do
    
    do m=1,nspec
       lv(1,CFS+m-1) =  Y(m)*b*vel(1)
       lv(2,CFS+m-1) =  Y(m)*b*vel(2)
       lv(3,CFS+m-1) =  Y(m)*b*vel(3)
       lv(4,CFS+m-1) = -Y(m)*b
       do n=1,nspec
          lv(CFS+n-1,CFS+m-1) = -Y(m)*d(n)
       end do
       lv(CFS+m-1,CFS+m-1) = lv(CFS+m-1,CFS+m-1) + 1.d0
    end do
    
    ! assemble right vectors
    rv(1,1) = vel(1) - c
    rv(2,1) = vel(2)
    rv(3,1) = vel(3)
    rv(4,1) = H - vel(1)*c
    do n=1,nspec
       rv(CFS+n-1,1) = Y(n)
    end do
    
    rv(1,2) = vel(1) + c
    rv(2,2) = vel(2)
    rv(3,2) = vel(3)
    rv(4,2) = H + vel(1)*c
    do n=1,nspec
       rv(CFS+n-1,2) = Y(n)
    end do
    
    rv(1,3) = 0.d0
    rv(2,3) = 1.d0
    rv(3,3) = 0.d0
    rv(4,3) = vel(2)
    do n=1,nspec
       rv(CFS+n-1,3) = 0.d0
    end do
    
    rv(1,4) = 0.d0
    rv(2,4) = 0.d0
    rv(3,4) = 1.d0
    rv(4,4) = vel(3)
    do n=1,nspec
       rv(CFS+n-1,4) = 0.d0
    end do
    
    do n=1,nspec
       rv(1,CFS+n-1) = vel(1)
       rv(2,CFS+n-1) = vel(2)
       rv(3,CFS+n-1) = vel(3)
       rv(4,CFS+n-1) = e + ek - dpdr(n)*gtinv
       do m=1,nspec
          rv(CFS+m-1,CFS+n-1) = 0.d0
       end do
       rv(CFS+n-1,CFS+n-1) = 1.d0
    end do

  end subroutine get_eigen_matrices_q


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

end module reconstruct_module
