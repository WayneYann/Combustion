module burner_module

  implicit none

  private

  public :: burn

contains

  subroutine burn(rho, YT, dt)
    double precision, intent(in   ) :: rho, dt
    double precision, intent(inout) :: YT(*)

    return

  end subroutine burn

end module burner_module


subroutine setfirst(frst)
  logical, intent(in) :: frst
end subroutine setfirst
