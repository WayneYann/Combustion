module burner_module

  implicit none

  private

  public :: burn

contains

  subroutine burn(np, rho, YT, dt, force_new_J)
    integer, intent(in) :: np
    double precision, intent(in   ) :: rho(*), dt
    double precision, intent(inout) :: YT(*)
    logical, intent(in) :: force_new_J

    return

  end subroutine burn


  subroutine compute_rhodYdt(np, rho, T, Y, rdYdt)
    integer, intent(in) :: np
    double precision, intent(in) :: rho(np), T(np), Y(np,nspecies)
    double precision, intent(out) :: rdYdt(np,nspecies)

    return

  end subroutine compute_rhodYdt

end module burner_module

