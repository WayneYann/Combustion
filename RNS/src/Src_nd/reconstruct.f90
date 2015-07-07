module reconstruct_module

  implicit none

  private

  public :: reconstruct_comp

contains

  ! L and R in UL and UR are relative to face
  ! UG1 and UG2 are at two Gauss points
  !
  ! Minimal ranges:
  !   U         : lo-2:hi+2 if UL & UR are not present; 
  !               lo-3:hi+3 if UL & UR are present
  !   UL  & UR  : lo  :hi+1
  !   UG1 & UG2 : lo  :hi
  subroutine reconstruct_comp(lo, hi, &
       Ulo, Uhi, &
       ULRlo, ULRhi, &
       UGlo, UGhi, &
       U, UL, UR, UG1, UG2)

    use meth_params_module, only : NVAR, UEDEN, UFS, NSPEC
    use weno_module, only : vweno5
    use eos_module, only : eos_get_eref

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
    

  subroutine normalize(lo,hi,U,Ulo,Uhi)

    use meth_params_module, only : NVAR, URHO, UEDEN, UFS, NSPEC
    use eos_module, only : eos_get_eref

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
