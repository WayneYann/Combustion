module make_plotfile_module

  use bl_types
  use variables
  use fabio_module
  use multifab_module
  use chemistry_module, only : nspecies, spec_names
  use probin_module, only : dm_in, plot_eint, plot_h, plot_divu, &
       plot_Y, plot_X, plot_hspec, plot_omegadot, &
       nOutFiles, lUsingNFiles, single_prec_plotfiles, prob_lo, prob_hi

  use make_plot_variables_module

  implicit none

  ! the total number of plot components
  integer, save :: n_plot_comps = 0
  integer, save :: icomp_rho=0, icomp_vel=0, icomp_pres=0, icomp_temp=0, icomp_eint=0, &
       icomp_h=0, icomp_divu=0, icomp_Y=0, icomp_X=0, icomp_hspec=0, icomp_omegadot=0

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

    if (plot_divu) then
       icomp_divu = get_next_plot_index(1)
    end if

    if (plot_Y) then
       icomp_Y = get_next_plot_index(nspecies)
    end if

    if (plot_X) then
       icomp_X = get_next_plot_index(nspecies)
    end if

    if (plot_hspec) then
       icomp_hspec = get_next_plot_index(nspecies)
    end if

    if (plot_omegadot) then
       icomp_omegadot = get_next_plot_index(nspecies)
    end if

    call make_plotvar_init(icomp_h, icomp_divu, icomp_omegadot)

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

    if (plot_divu) then
       plot_names(icomp_divu) = "divu"
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

  end subroutine get_plot_names


  subroutine make_plotfile(dirname, la, U, plot_names, time, dx, write_pf_time)

    character(len=*) , intent(in   ) :: dirname
    type(layout)     , intent(in   ) :: la
    type(multifab)   , intent(inout) :: U
    character(len=20), intent(in   ) :: plot_names(:)
    real(dp_t)       , intent(in   ) :: time, dx(U%dim)
    real(dp_t)       , intent(  out) :: write_pf_time

    ! dimensioned as an array of size 1 for fabio_ml_multifab_write_d
    type(multifab) :: plotdata(1), Q

    ! dimensioned as an array of size 0 for fabio_ml_multifab_write_d
    integer :: rr(0), prec, ngu, ngq
    real(dp_t) :: writetime1, writetime2

    if (single_prec_plotfiles) then
       prec = FABIO_SINGLE
    else
       prec = FABIO_DOUBLE
    endif

    ngu = nghost(U)

    if (plot_divu) then
       ngq = ngu
       call multifab_fill_boundary(U)
    else
       ngq = 0
    end if

    call multifab_build(Q,la,nprim, ngq)

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

    if (plot_h) then
       call make_plotvar(plotdata(1),icomp_h, Q, dx)
    end if

    if (plot_Y) then
       call multifab_copy_c(plotdata(1),icomp_Y, Q,qy1, nspecies)       
    end if

    if (plot_X) then
       call multifab_copy_c(plotdata(1),icomp_X, Q,qx1, nspecies)       
    end if

    if (plot_hspec) then
       call multifab_copy_c(plotdata(1),icomp_hspec, Q,qh1, nspecies)       
    end if

    if (plot_omegadot) then
       call make_plotvar(plotdata(1),icomp_omegadot, Q, dx)
    end if

    if (parallel_IOProcessor()) then
       write(6,*) ''
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
