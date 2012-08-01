module make_plotfile_module

  use bl_types
  use variables
  use fabio_module
  use multifab_module
  use chemistry_module, only : nspecies, spec_names
  use probin_module, only : dm_in, plot_Y, plot_X, plot_h, nOutFiles, lUsingNFiles, &
       single_prec_plotfiles, prob_lo, prob_hi

  implicit none

  private
  public :: get_plot_names, make_plotfile

contains

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

  end subroutine get_plot_names


  subroutine make_plotfile(dirname, la, U, plot_names, time, dx, write_pf_time)
    character(len=*) , intent(in   ) :: dirname
    type(layout)     , intent(in   ) :: la
    type(multifab)   , intent(in   ) :: U
    character(len=20), intent(in   ) :: plot_names(:)
    real(dp_t)       , intent(in   ) :: time, dx(U%dim)
    real(dp_t)       , intent(  out) :: write_pf_time

    ! dimensioned as an array of size 1 for fabio_ml_multifab_write_d
    type(multifab) :: plotdata(1), Q

    ! dimensioned as an array of size 0 for fabio_ml_multifab_write_d
    integer :: rr(0), prec
    real(dp_t) :: writetime1, writetime2

    if (single_prec_plotfiles) then
       prec = FABIO_SINGLE
    else
       prec = FABIO_DOUBLE
    endif

    call multifab_build(Q,la,nprim, nghost(U))

    call ctoprim(U, Q, ng=0)

    call multifab_build(plotdata(1),la,n_plot_comps,0)

    ! copy density, velocities, pressure, and temperature
    ! that's 6 variables
    call multifab_copy_c(plotdata(1),1,Q,1, 6)

    if (plot_Y) then
       call multifab_copy_c(plotdata(1),icomp_Y, Q,qy1, nspecies)       
    end if

    if (plot_X) then
       call multifab_copy_c(plotdata(1),icomp_X, Q,qx1, nspecies)       
    end if

    if (plot_h) then
       call multifab_copy_c(plotdata(1),icomp_h, Q,qh1, nspecies)       
    end if

    if (parallel_IOProcessor()) then
       write(6,*) 'Writing state to plotfile ',trim(dirname)
    end if

    writetime1 = parallel_wtime()

    call fabio_ml_multifab_write_d(plotdata, rr, dirname, plot_names, &
         la%lap%pd, prob_lo, prob_hi, time, dx, &
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
