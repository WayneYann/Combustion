module chemterm_module

  use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UTEMP, UFS, NSPEC
  use burner_module, only : burn
  use eos_module, only : eos_get_T

  implicit none

  private

  public :: chemterm

contains

  subroutine chemterm(lo, hi, U, Ulo, Uhi, dUdt, Utlo, Uthi, dt)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1), Utlo(1), Uthi(1)
    double precision, intent(in   ) ::    U( Ulo(1): Uhi(1),NVAR)
    double precision, intent(inout) :: dUdt(Utlo(1):Uthi(1),NVAR)
    double precision, intent(in) :: dt

    integer :: i, n, iwrk
    double precision :: rho, rhoinv, ei, rwrk, dtinv
    double precision :: YT(nspec+1), YTout(nspec+1)

    dtinv = 1.d0/dt

    do i=lo(1),hi(1)

       rho = u(i,URHO)
       rhoinv = 1.d0/rho

       ei = rhoinv*( u(i,UEDEN) - 0.5d0*rhoinv*u(i,UMX)**2 )
       YT(1:nspec)  = u(i,UFS:UFS+nspec-1) * rhoinv

       call eos_get_T(YT(nspec+1), ei, YT(1:nspec))
       
       call burn(rho, YT, YTout, dt)

       do n=1,nspec
          dUdt(i,UFS+n-1) = dUdt(i,UFS+n-1) + rho*(YTout(n) - YT(n)) * dtinv
       end do

       dUdt(i,UTEMP) = dUdt(i,UTEMP) + (YTout(nspec+1) - YT(nspec+1)) * dtinv

    end do

  end subroutine chemterm

end module chemterm_module
