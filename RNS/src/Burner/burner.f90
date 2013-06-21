module burner_module

  use chemistry_module, only : nspecies, spec_names
  use vode_module, only : verbose, itol, rtol, atol, MF, &
       voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar

  implicit none

  private

  public :: burn

contains

  subroutine burn(rho, YT, dt)
    double precision, intent(in   ) :: rho, dt
    double precision, intent(inout) :: YT(nspecies+1)

    external jac, f_rhs, dvode
    
    ! vode stuff
    integer, parameter :: itask=1, iopt=1
    integer :: istate, ifail

    double precision :: time

    voderpar(1) = rho
       
    istate = 1
    time = 0.d0

    call dvode(f_rhs, nspecies+1, YT, time, dt, itol, rtol, atol, itask, &
         istate, iopt, voderwork, lvoderwork, vodeiwork, lvodeiwork, &
         jac, MF, voderpar, vodeipar)

    if (verbose .ge. 1) then
       write(6,*) '......dvode done:'
       write(6,*) ' last successful step size = ',voderwork(11)
       write(6,*) '          next step to try = ',voderwork(12)
       write(6,*) '   integrated time reached = ',voderwork(13)
       write(6,*) '      number of time steps = ',vodeiwork(11)
       write(6,*) '              number of fs = ',vodeiwork(12)
       write(6,*) '              number of Js = ',vodeiwork(13)
       write(6,*) '    method order last used = ',vodeiwork(14)
       write(6,*) '   method order to be used = ',vodeiwork(15)
       write(6,*) '            number of LUDs = ',vodeiwork(19)
       write(6,*) ' number of Newton iterations ',vodeiwork(20)
       write(6,*) ' number of Newton failures = ',vodeiwork(21)
       if (ISTATE.eq.-4 .or. ISTATE.eq.-5) then
          ifail = vodeiwork(16)
          if (ifail .eq. nspecies+1) then
             write(6,*) '   T has the largest error'
          else
             write(6,*) '   spec with largest error = ', trim(spec_names(ifail))
          end if
          call flush(6)
       end if
    end if

    if (istate < 0) then
       print *, 'chemsolv: VODE failed'
       print *, 'istate = ', istate, ' time =', time
       call bl_error("ERROR in burn: VODE failed")
    end if

  end subroutine burn

end module burner_module
