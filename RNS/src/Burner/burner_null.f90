module burner_module

  implicit none

  private

  public :: burn, compute_rhodYdt

contains

  subroutine burn(np, rho, YT, dt, force_new_J)
    integer, intent(in) :: np
    double precision :: rho(*), dt
    double precision :: YT(*)
    logical, intent(in) :: force_new_J

    return

  end subroutine burn


  subroutine compute_rhodYdt(np, rho, T, Y, rdYdt)
    integer, intent(in) :: np
    double precision :: rho(np), T(np), Y(*)
    double precision :: rdYdt(*)

    return

  end subroutine compute_rhodYdt

end module burner_module

