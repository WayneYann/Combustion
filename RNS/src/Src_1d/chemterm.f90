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

    integer :: i, n, g
    logical :: force_new_J
    double precision :: rho(2), rhoinv, ei
    double precision :: YT(nspec+1,2)
    double precision, allocatable :: UG(:,:,:)

    allocate(UG(lo(1):hi(1),NVAR,2))

    do n=1,NVAR
       call cellavg2gausspt_1d(lo(1),hi(1), U(:,n), Ulo(1),Uhi(1), UG(:,n,1), UG(:,n,2), lo(1),hi(1))
    end do
    
    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts
    
    do i=lo(1),hi(1)

       do g=1,2

          rho(g) = 0.d0
          do n=UFS,UFS+nspec-1
             YT(n-UFS+1,g) = UG(i,n,g)
             rho(g) = rho(g) + UG(i,n,g)
          end do
          rhoinv = 1.d0/rho(g)

          YT(1:nspec,g) = YT(1:nspec,g) * rhoinv
          YT(nspec+1,g) = UG(i,UTEMP,g)

          ei = rhoinv*( UG(i,UEDEN,g) - 0.5d0*rhoinv*UG(i,UMX,g)**2 )

          call eos_get_T(YT(nspec+1,g), ei, YT(1:nspec,g))
       
       end do

       call burn(2, rho, YT, dt, force_new_J)

       force_new_J = .false.

       do n=1,nspec
          U(i,UFS+n-1) = 0.5d0*(rho(1)*YT(n,1) + rho(2)*YT(n,2))
       end do

    end do

    deallocate(UG)

  end subroutine chemterm

end module chemterm_module
