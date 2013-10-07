module burner_module

  use chemistry_module, only : nspecies, spec_names

  implicit none

  private

  public :: burn

contains

  subroutine burn(rho, YT, dt)
    use meth_params_module, only : use_vode
    double precision, intent(in   ) :: rho, dt
    double precision, intent(inout) :: YT(nspecies+1)

    if (use_vode) then
       call burn_vode(rho, YT, dt)
    else
       call burn_bdf(rho, YT, dt)
    end if

  end subroutine burn


  subroutine burn_vode(rho, YT, dt)
    use vode_module, only : verbose, itol, rtol, atol, vode_MF=>MF, always_new_j, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar
    double precision, intent(in   ) :: rho, dt
    double precision, intent(inout) :: YT(nspecies+1)

    external f_jac, f_rhs, dvode
    
    ! vode stuff
    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate, ifail

    double precision :: time

    voderpar(1) = rho
       
    istate = 1
    time = 0.d0

    if (always_new_j) call setfirst(.true.)

    MF = vode_MF  ! vode might change its sign!
    call dvode(f_rhs, nspecies+1, YT, time, dt, itol, rtol, atol, itask, &
         istate, iopt, voderwork, lvoderwork, vodeiwork, lvodeiwork, &
         f_jac, MF, voderpar, vodeipar)

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

  end subroutine burn_vode


  subroutine burn_bdf(rho_in, YT, dt)
    use bdf
    use bdf_data, only : ts, reuse_jac
    use feval, only : f_rhs, f_jac, rho
    double precision, intent(in   ) :: rho_in, dt
    double precision, intent(inout) :: YT(nspecies+1)

    double precision :: t0, t1, y1(nspecies+1)
    integer :: neq, ierr
    logical :: reset

    rho = rho_in

    neq = nspecies+1
    t0 = 0.d0
    t1 = dt

    reset = .true.

    call bdf_advance(ts, f_rhs, f_jac, neq, 1, YT, t0, y1, t1, dt, reset, reuse_jac, ierr)

    if (ierr .ne. 0) then
       print *, 'chemsolv: BDF failed'
       call bl_error("ERROR in burn: BDF failed")       
    end if

  end subroutine burn_bdf

end module burner_module
