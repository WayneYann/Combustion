module burner_module

  implicit none

  private

  public :: burn, compute_rhodYdt, splitburn, beburn

contains

  subroutine burn(np, rho, YT, dt, force_new_J, ierr)
    integer, intent(in) :: np
    double precision :: rho(*), dt
    double precision :: YT(*)
    logical, intent(in) :: force_new_J
    integer, intent(out), optional :: ierr
    if (present(ierr)) ierr=0
    return
  end subroutine burn


  subroutine compute_rhodYdt(np, rho, T, Y, rdYdt)
    integer, intent(in) :: np
    double precision :: rho(np), T(np), Y(*)
    double precision :: rdYdt(*)
    return
  end subroutine compute_rhodYdt


  subroutine splitburn(np, rho0, Y0, rho, YT, dt)
    integer, intent(in) :: np
    double precision, intent(in   ) :: rho0, rho(np), dt
    double precision, intent(in   ) :: Y0(*)
    double precision, intent(inout) :: YT(*)
    return
  end subroutine splitburn


  subroutine beburn(rho0, Y0, rho, YT, dt, g, ierr)
    integer, intent(in) :: g
    double precision, intent(in   ) :: rho0, rho, dt
    double precision, intent(in   ) :: Y0(*)
    double precision, intent(inout) :: YT(*)
    integer, intent(out), optional :: ierr
    if (present(ierr)) ierr=0
    return
  end subroutine beburn

end module burner_module

