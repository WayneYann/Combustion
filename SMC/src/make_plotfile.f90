module make_plotfile_module

  use bl_types
  use fabio_module
  use make_plot_variables_module
  use multifab_module
  use threadbox_module, only : tb_multifab_setval
  use variables_module

  use chemistry_module, only : nspecies, spec_names, get_species_index
  use probin_module, only : dm_in, plot_eint, plot_h, plot_rhoh, plot_divu, plot_magvort, &
       plot_Y, plot_X, plot_wbar, plot_hspec, plot_omegadot, plot_dYdt, plot_heatRelease, &
       plot_fuelConsumption, fuel_name, &
       nOutFiles, lUsingNFiles, single_prec_plotfiles, prob_lo, prob_hi

  use plotvar_index_module

  implicit none

  private
  public :: n_plot_comps, get_plot_names, make_plotfile, init_plot_variables

contains

  function get_next_plot_index(num) result (next)

    ! return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num
    integer :: num, next

    next = n_plot_comps + 1
    n_plot_comps = n_plot_comps + num

    return
  end function get_next_plot_index


  subroutine init_plot_variables()

    icomp_rho  = get_next_plot_index(1)
    icomp_vel  = get_next_plot_index(dm_in)
    icomp_pres = get_next_plot_index(1)
    icomp_temp = get_next_plot_index(1)

    if (plot_eint) then
       icomp_eint = get_next_plot_index(1)
    end if

    if (plot_h) then
       icomp_h = get_next_plot_index(1)
    end if

    if (plot_rhoh) then
       icomp_rhoh = get_next_plot_index(1)
    end if

    if (plot_divu) then
       icomp_divu = get_next_plot_index(1)
    end if

    if (plot_magvort) then
       icomp_magvort = get_next_plot_index(1)
    end if

    if (plot_Y) then
       icomp_Y = get_next_plot_index(nspecies)
    end if

    if (plot_X) then
       icomp_X = get_next_plot_index(nspecies)
    end if

    if (plot_wbar) then
       icomp_wbar = get_next_plot_index(1)
    end if

    if (plot_hspec) then
       icomp_hspec = get_next_plot_index(nspecies)
    end if

    nburn = 0
    if (plot_omegadot) then
       icomp_omegadot = get_next_plot_index(nspecies)
       icomp_burn = icomp_omegadot
       ib_omegadot = 1
       nburn = nburn + 1
    end if

    if (plot_dYdt) then
       icomp_dYdt = get_next_plot_index(nspecies)
       if (nburn > 0) then
          ib_dYdt = icomp_dYdt - icomp_burn + 1
       else
          icomp_burn = icomp_dYdt
          ib_dYdt = 1
       end if
       nburn = nburn + 1
    end if

    if (plot_heatRelease) then
       icomp_heatRelease = get_next_plot_index(1)
       if (nburn > 0) then
          ib_heatRelease = icomp_heatRelease - icomp_burn + 1
       else
          icomp_burn = icomp_heatRelease
          ib_heatRelease = 1
       end if
       nburn = nburn + 1
    end if

    if (plot_fuelConsumption) then
       ifuel = get_species_index(fuel_name)
       if (ifuel > 0) then
          icomp_fuelConsumption = get_next_plot_index(1)
          if (nburn > 0) then
             ib_fuelConsumption = icomp_fuelConsumption - icomp_burn + 1
          else
             icomp_burn = icomp_fuelConsumption
             ib_fuelConsumption = 1
          end if
          nburn = nburn + 1
       end if
    end if

  end subroutine init_plot_variables


  subroutine get_plot_names(plot_names)
    character(len=20), intent(inout) :: plot_names(:)

    integer :: i

    plot_names(icomp_rho  ) = "density"
    plot_names(icomp_vel  ) = "x_vel"
    plot_names(icomp_vel+1) = "y_vel"
    if (dm_in .eq. 3) then
       plot_names(icomp_vel+2) = "z_vel"
    end if
    plot_names(icomp_pres) = "pressure"
    plot_names(icomp_temp) = "temperature"

    if (plot_eint) then
       plot_names(icomp_eint) = "eint"
    end if

    if (plot_h) then
       plot_names(icomp_h) = "h"
    end if

    if (plot_rhoh) then
       plot_names(icomp_rhoh) = "rhoh"
    end if

    if (plot_divu) then
       plot_names(icomp_divu) = "divu"
    end if

    if (plot_magvort) then
       plot_names(icomp_magvort) = "magvort"
    end if

    if (plot_Y) then
       do i=1,nspecies
          plot_names(icomp_Y+i-1) = "Y("//trim(spec_names(i))//")"
       end do
    end if

    if (plot_X) then
       do i=1,nspecies
          plot_names(icomp_X+i-1) = "X("//trim(spec_names(i))//")"
       end do
    end if

    if (plot_wbar) then
       plot_names(icomp_wbar) = "wbar"
    end if

    if (plot_hspec) then
       do i=1,nspecies
          plot_names(icomp_hspec+i-1) = "h("//trim(spec_names(i))//")"
       end do
    end if

    if (plot_omegadot) then
       do i=1,nspecies
          plot_names(icomp_omegadot+i-1) = "rho*omgdot("//trim(spec_names(i))//")"
       end do
    end if

    if (plot_dYdt) then
       do i=1,nspecies
          plot_names(icomp_dYdt+i-1) = "Ydot("//trim(spec_names(i))//")"
       end do
    end if

    if (plot_heatRelease) then
       plot_names(icomp_heatRelease) = "HeatRelease"
    end if

    if (plot_fuelConsumption .and. ifuel>0) then
       plot_names(icomp_fuelConsumption) = "FuelConsumptionRate"
    end if

  end subroutine get_plot_names


  subroutine make_plotfile(dirname, la, U, plot_names, time, dx, write_pf_time)

    character(len=*) , intent(in   ) :: dirname
    type(layout)     , intent(in   ) :: la
    type(multifab)   , intent(inout) :: U
    character(len=20), intent(in   ) :: plot_names(:)
    real(dp_t)       , intent(in   ) :: time, dx(3)
    real(dp_t)       , intent(  out) :: write_pf_time

    ! dimensioned as an array of size 1 for fabio_ml_multifab_write_d
    type(multifab) :: plotdata(1), Q

    ! dimensioned as an array of size 0 for fabio_ml_multifab_write_d
    integer :: rr(0), prec, ngu, ngq, dm
    real(dp_t) :: writetime1, writetime2

    if (single_prec_plotfiles) then
       prec = FABIO_SINGLE
    else
       prec = FABIO_DOUBLE
    endif

    dm  = U%dim
    ngu = nghost(U)

    if (plot_divu .or. plot_magvort) then
       ngq = ngu
       call multifab_fill_boundary(U)
    else
       ngq = 0
    end if

    call multifab_build(Q,la,nprim, ngq)
    call tb_multifab_setval(Q, 0.d0, all=.true.)

    call ctoprim(U, Q, ngq)

    call multifab_build(plotdata(1),la,n_plot_comps,0)

    ! copy up to temperature
    call multifab_copy_c(plotdata(1),1,Q,1, qtemp)

    if (plot_eint) then
       call multifab_copy_c(plotdata(1),icomp_eint, Q,qe, 1)
    end if    

    if (plot_divu) then
       call make_plotvar(plotdata(1),icomp_divu, Q, dx)
    end if

    if (plot_magvort) then
       call make_plotvar(plotdata(1),icomp_magvort, Q, dx)
    end if

    if (plot_h) then
       call make_plotvar(plotdata(1),icomp_h, Q, dx)
    end if

    if (plot_rhoh) then
       call make_plotvar(plotdata(1),icomp_rhoh, Q, dx)
    end if

    if (plot_Y) then
       call multifab_copy_c(plotdata(1),icomp_Y, Q,qy1, nspecies)       
    end if

    if (plot_X) then
       call multifab_copy_c(plotdata(1),icomp_X, Q,qx1, nspecies)       
    end if

    if (plot_wbar) then
       call make_plotvar(plotdata(1),icomp_wbar, Q, dx)
    end if

    if (plot_hspec) then
       call multifab_copy_c(plotdata(1),icomp_hspec, Q,qh1, nspecies)       
    end if

    if (nburn > 0) then
       call make_plotvar(plotdata(1),icomp_burn, Q, dx)
    end if

    if (parallel_IOProcessor()) then
       write(6,*) ''
       write(6,*) 'Writing state to plotfile ',trim(dirname)
    end if

    writetime1 = parallel_wtime()

    call fabio_ml_multifab_write_d(plotdata, rr, dirname, plot_names, &
         la%lap%pd, prob_lo(1:dm), prob_hi(1:dm), time, dx(1:dm), &
         nOutFiles = nOutFiles, &
         lUsingNFiles = lUsingNFiles, prec = prec)

    writetime2 = parallel_wtime() - writetime1
    call parallel_reduce(writetime1, writetime2, MPI_MAX, proc=parallel_IOProcessorNode())
    if (parallel_IOProcessor()) then
       print*,'Time to write plotfile: ',writetime1,' seconds'
       print*,''
    end if

    write_pf_time = writetime1

    call multifab_destroy(plotdata(1))
    call multifab_destroy(Q)

  end subroutine make_plotfile

end module make_plotfile_module
