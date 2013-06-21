module chemterm_module

  use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UTEMP, UFS, NSPEC
  use burner_module, only : burn
  use eos_module, only : eos_get_T

  implicit none

  private

  public :: chemterm

contains

  subroutine chemterm(lo, hi, U, Ulo, Uhi, dt)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(in) :: dt

    integer :: i
    double precision :: rho, rhoinv, ei
    double precision :: YT(nspec+1)

    do i=lo(1),hi(1)

       rho = u(i,URHO)
       rhoinv = 1.d0/rho

       ei = rhoinv*( u(i,UEDEN) - 0.5d0*rhoinv*u(i,UMX)**2 )
       YT(1:nspec)  = u(i,UFS:UFS+nspec-1) * rhoinv

       call eos_get_T(YT(nspec+1), ei, YT(1:nspec))
       
       call burn(rho, YT, dt)

       U(i,UFS:UFS+nspec-1) = rho * YT(1:nspec)

    end do

  end subroutine chemterm

end module chemterm_module
