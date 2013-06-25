module chemterm_module

  use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UTEMP, UFS, NSPEC
  use burner_module, only : burn
  use eos_module, only : eos_get_T
  use weno_module, only : cellavg2gausspt_1d

  implicit none

  private

  public :: chemterm

contains

  subroutine chemterm(lo, hi, U, Ulo, Uhi, dt)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(in) :: dt

    integer :: i, n
    double precision :: rho1, rho2, rhoinv, ei
    double precision :: YT1(nspec+1), YT2(nspec+1)
    double precision, allocatable :: U1(:,:), U2(:,:)

    allocate(U1(lo(1):hi(1),NVAR))
    allocate(U2(lo(1):hi(1),NVAR))

    do n=1,NVAR
       call cellavg2gausspt_1d(U(:,n), Ulo(1), Uhi(1), U1(:,n), U2(:,n), lo(1), hi(1))
    end do

    do i=lo(1),hi(1)

       ! Guass point 1
       rho1 = 0.d0
       do n=UFS,UFS+nspec-1
          YT1(n-UFS+1) = U1(i,n)
          rho1 = rho1 + U1(i,n)
       end do
       rhoinv = 1.d0/rho1

       YT1(1:nspec) = YT1(1:nspec) * rhoinv
       YT1(nspec+1) = U1(i,UTEMP)

       ei = rhoinv*( U1(i,UEDEN) - 0.5d0*rhoinv*U1(i,UMX)**2 )

       call eos_get_T(YT1(nspec+1), ei, YT1(1:nspec))
       
       call burn(rho1, YT1, dt)

       ! Guass point 2
       rho2 = 0.d0
       do n=UFS,UFS+nspec-1
          YT2(n-UFS+1) = U2(i,n)
          rho2 = rho2 + U2(i,n)
       end do
       rhoinv = 1.d0/rho2

       YT2(1:nspec) = YT2(1:nspec) * rhoinv
       YT2(nspec+1) = U2(i,UTEMP)

       ei = rhoinv*( U2(i,UEDEN) - 0.5d0*rhoinv*U2(i,UMX)**2 )

       call eos_get_T(YT2(nspec+1), ei, YT2(1:nspec))
       
       call burn(rho2, YT2, dt)

       do n=1,nspec
          U(i,UFS+n-1) = 0.5d0*(rho1*YT1(n) + rho2*YT2(n))
       end do

    end do

    deallocate(U1, U2)

  end subroutine chemterm

end module chemterm_module
