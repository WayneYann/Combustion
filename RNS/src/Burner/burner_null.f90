module burner_module

  implicit none

  private

  public :: burn, compute_rhodYdt, init_burn_linear, burn_linear

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


  subroutine init_burn_linear(rho0, Y0, dt)
    double precision, intent(in) :: rho0, Y0(*), dt
    return
  end subroutine init_burn_linear


  subroutine burn_linear(dY)
    double precision, intent(inout) :: dY(*)
    return;
  end subroutine burn_linear

end module burner_module

