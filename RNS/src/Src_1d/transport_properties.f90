module transport_properties

  use meth_params_module, only : NVAR, NSPEC, URHO, UMX, UEDEN, UTEMP, UFS
  use eglib_module

  implicit none

  ! eglib parameters
  integer, save :: ITLS=-1, IFLAG=-1

  private

  public :: get_transport_properties

contains

  subroutine get_transport_properties(lo, hi, U, mu, xi, lam, Ddiag)
    integer, intent(in) :: lo, hi
    double precision, intent(in) :: U(lo:hi,NVAR)
    double precision, intent(out) ::    mu(lo:hi)
    double precision, intent(out) ::    xi(lo:hi)
    double precision, intent(out) ::   lam(lo:hi)
    double precision, intent(out) :: Ddiag(lo:hi,NSPEC)

    integer :: i
    double precision :: rhoInv, vx, et, Tt, Xt(nspec), Yt(nspec), Cpt(nspec)
    
    logical, save :: first_call = .true.

    if (first_call) then
       first_call = .false.
!       if (use_bulk_viscosity) then
       if (.true.) then
          ITLS  = 1 
          IFLAG = 5
       else
          ITLS  = 1
          IFLAG = 3
       end if
    end if

    call eglib_init(nspec, 1, ITLS, IFLAG)

    do i = lo,hi

       rhoInv = 1.0d0/U(i,URHO)
       vx = U(i,UMX)*rhoInv 
       et = U(i,UEDEN)*rhoInv - 0.5d0*vx*vx
       
       Yt = U(i,UFS:UFS+NSPEC-1)*rhoInv
!       call eos_get_T(U(i,UTEMP), e, Y)
       

       call egspar(Tt, Xt, Yt, Cpt, egwork, egiwork)

    end do

  end subroutine get_transport_properties

end module transport_properties
