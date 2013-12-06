module reconstruct_module

  use meth_params_module, only : ndim, NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, &
       UFS, NSPEC, NCHARV, CFS

  implicit none

  private

  public :: reconstruct

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

    use weno_module, only : weno5
    use eos_module, only : eos_given_ReY, eos_get_eref

    integer, intent(in) :: lo, hi, Ulo, Uhi, ULRlo, ULRhi, UGlo, UGhi, U0lo, U0hi 
    integer, intent(in), optional :: dir
    double precision, intent(in),target           :: U ( Ulo: Uhi,NVAR)
    double precision, intent(in),target, optional :: U0(U0lo:U0hi,NVAR)
    double precision, intent(out), dimension(ULRlo:ULRhi,NVAR), optional :: UL, UR
    double precision, intent(out), dimension( UGlo: UGhi,NVAR), optional :: UG1, UG2

    integer :: i, ii, ivar, m, n, ivel(3), idir, iextra
    double precision :: egv(NCHARV,NCHARV)
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

       if (do_face .and. i.ne.hi+1) then
          do n=1,NCHARV
             UL(i+1,ivel(1)) = UL(i+1,ivel(1)) + vp(n)*egv(1,n)*vflag(1)
             UL(i+1,ivel(2)) = UL(i+1,ivel(2)) + vp(n)*egv(2,n)*vflag(2)
             UL(i+1,ivel(3)) = UL(i+1,ivel(3)) + vp(n)*egv(3,n)*vflag(3)
             UL(i+1,UEDEN  ) = UL(i+1,UEDEN  ) + vp(n)*egv(4,n)
             do m=1,nspec
                UL(i+1,UFS+m-1) = UL(i+1,UFS+m-1) + vp(n)*egv(CFS+m-1,n)
             end do
          end do

          do m=1,nspec
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
             do m=1,nspec
                UR(i,UFS+m-1) = UR(i,UFS+m-1) + vm(n)*egv(CFS+m-1,n)
             end do
          end do

          do m=1,nspec
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
             do m=1,nspec
                UG1(i,UFS+m-1) = UG1(i,UFS+m-1) + vg1(n)*egv(CFS+m-1,n)
             end do
          end do
          
          do m=1,nspec
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
             do m=1,nspec
                UG2(i,UFS+m-1) = UG2(i,UFS+m-1) + vg2(n)*egv(CFS+m-1,n)
             end do
          end do
          
          do m=1,nspec
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

  end subroutine reconstruct
    

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

end module reconstruct_module
