module burner_module

  use chemistry_module, only : nspecies
  use vode_module, only : itol, rtol, atol, MF_NUMERICAL_JAC, &
       voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar

  implicit none

  private

  public :: burn

contains

  subroutine burn(rho, YT, YTout, dt)
    double precision, intent(in ) :: rho, YT(nspecies+1), dt
    double precision, intent(out) :: YTout(nspecies+1)

    external jac, f_rhs, dvode
    
    ! vode stuff
    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate

    double precision :: time

    MF = MF_NUMERICAL_JAC

    voderpar(1) = rho
       
    YTout = YT

    istate = 1
    time = 0.d0

    call dvode(f_rhs, nspecies+1, YT, time, dt, itol, rtol, atol, itask, &
         istate, iopt, voderwork, lvoderwork, vodeiwork, lvodeiwork, &
         jac, MF, voderpar, vodeipar)

    if (istate < 0) then
       print *, 'chemsolv: VODE failed'
       print *, 'istate = ', istate, ' time =', time
       call bl_error("ERROR in burn: VODE failed")
    end if

  end subroutine burn

end module burner_module
