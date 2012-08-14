subroutine smc()

  use advance_module
  use checkpoint_module
  use chemistry_module
  use derivative_stencil_module
  use initialize_module
  use layout_module
  use make_plotfile_module
  use multifab_module
  use omp_module
  use probin_module
  use runtime_init_module
  use sdcquad_module
  use time_module
  use variables_module

  use cputime_module, only: start_cputime_clock


  implicit none

  integer :: dm, numcell, i
  integer :: init_step, istep

  real(dp_t) :: dt
  real(dp_t)  , pointer     :: dx(:)

  integer :: last_plt_written,last_chk_written
  character(len=5)               :: plot_index, check_index
  character(len=6)               :: plot_index6, check_index6
  character(len=256)             :: plot_file_name, check_file_name
  character(len=20), allocatable :: plot_names(:)

!  logical :: dump_plotfile, dump_checkpoint
  real(dp_t) :: write_pf_time
  
  type(layout)   :: la
  type(multifab) :: U
  type(sdcquad)  :: sdc

  ! keep track of cputime
  call start_cputime_clock()

  last_plt_written = -1
  last_chk_written = -1

  call runtime_init()

  call stencil_init()

  call chemistry_init()
  if (verbose .ge. 1) then
     if (parallel_IOProcessor()) then
        print *, ''
        write(*,'(A,1X,I0,1X,A)', advance='no') "Chemistry model has", nspecies, "species:"
        do i=1,nspecies
           write(*, '(3X,A)', advance='no') trim(spec_names(i))
        end do
        print *, ''
     end if
  end if

  call init_variables()
  call init_plot_variables()

  allocate(plot_names(n_plot_comps))
  call get_plot_names(plot_names)

  if (restart >= 0) then

     if (restart <= 99999) then
        write(unit=check_index,fmt='(i5.5)') restart
        check_file_name = trim(check_base_name) // check_index
     else
        write(unit=check_index6,fmt='(i6.6)') restart
        check_file_name = trim(check_base_name) // check_index6
     endif
     
     if (parallel_IOProcessor() .and. verbose.ge.1) then
        print*,""
        print*,"Restarting from", check_file_name
     end if

     call initialize_from_restart(check_file_name, la,dt,dx,U)

  else 

     if (parallel_IOProcessor() .and. verbose.ge.1) then
        print*,""
        print*,"Starting from scratch"
     end if

     call initialize_from_scratch(la,dt,dx,U)

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! error checking
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dm = dm_in
  if (dm .ne. 3) then 
     call bl_error('SMC can only do 3D')
  end if

  ! check to make sure dimensionality is consistent in the inputs file
  if (dm .ne. get_dim(la)) then 
     call bl_error('dm_in not properly set in inputs file')
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! print processor and grid info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (parallel_IOProcessor()) then
     print *, ' '     
     print *, 'number of MPI processes = ', parallel_nprocs()
     print *, 'number of threads       = ', omp_get_max_threads()
     print *, ' '
     print *, 'number of dimensions    = ', dm
     print *, 'number of boxes         = ', nboxes(U)
     print *, ' '
  end if
  
  if (verbose .ge. 1) then
     numcell = multifab_volume(U,.false.)
     if (parallel_IOProcessor()) then
        print*,"Number of valid cells:              ",numcell
     end if
     numcell = multifab_volume(U,.true.)
     if (parallel_IOProcessor()) then
        print*,"Number of valid cells + ghost cells:",numcell
        print*,""
     end if
  end if


  if (restart < 0) then

     istep = 0

     if (chk_int > 0) then
        write(unit=check_index,fmt='(i5.5)') istep
        check_file_name = trim(check_base_name) // check_index

        call checkpoint_write(check_file_name, la, U, dt)
        
        last_chk_written = istep
     end if

     if (plot_int > 0 .or. plot_deltat > ZERO) then
        write(unit=plot_index,fmt='(i5.5)') istep
        plot_file_name = trim(plot_base_name) // plot_index
        
        call make_plotfile(plot_file_name,la,U,plot_names,time,dx,write_pf_time)

        call write_job_info(plot_file_name, la, write_pf_time)
        
        last_plt_written = istep
     end if
  end if

  if (restart < 0) then
     init_step = 1 
  else
     init_step = restart + 1
  end if

  if (advance_method == 2) then
     call create(sdc, 1, 3)  ! XXX: use probin here...
     sdc%tol_residual = 1.0d-8
     call build(sdc)
  end if

  if ( parallel_IOProcessor()) then
     print*,""
     print*,"BEGIN MAIN EVOLUTION LOOP"
     print*,""
  end if

  if ( (max_step >= init_step) .and. (time < stop_time .or. stop_time < 0.d0) ) then

     do istep = init_step, max_step

        if ( verbose .ge. 1 ) then
           if ( parallel_IOProcessor() ) then
              print*, ' '
              print *, 'MEMORY STATS AT START OF TIMESTEP ', istep
              print*, ' '
           end if
           call print(multifab_mem_stats(),    "    multifab")
           call print(fab_mem_stats(),         "         fab")
           call print(boxarray_mem_stats(),    "    boxarray")
           call print(layout_mem_stats(),      "      layout")
           call print(boxassoc_mem_stats(),    "    boxassoc")
!           call print(fgassoc_mem_stats(),     "     fgassoc")
!           call print(syncassoc_mem_stats(),   "   syncassoc")
           call print(copyassoc_mem_stats(),   "   copyassoc")
!           call print(fluxassoc_mem_stats(),   "   fluxassoc")
           if ( parallel_IOProcessor() ) print*, ''
        end if

        if (parallel_IOProcessor()) then
           print*,'Advancing time step',istep,'time = ',time
           print*, ""
        end if

        call advance(U,dt,dx,sdc,istep)

        time = time + dt

        if ( parallel_IOProcessor() ) then
           print *, 'End of step', istep,'time = ', time
        end if

        if (chk_int > 0) then
           if ( mod(istep,chk_int) .eq. 0 ) then
              
              if (istep <= 99999) then
                 write(unit=check_index,fmt='(i5.5)') istep
                 check_file_name = trim(check_base_name) // check_index
              else
                 write(unit=check_index6,fmt='(i6.6)') istep
                 check_file_name = trim(check_base_name) // check_index6
              endif
              
              call checkpoint_write(check_file_name, la, U, dt)
              
              last_chk_written = istep
              
           end if
        end if

        if (plot_int > 0 .or. plot_deltat > ZERO) then
           if ( (plot_int > 0 .and. mod(istep,plot_int) .eq. 0) .or. &
                (plot_deltat > ZERO .and. &
                mod(time - dt,plot_deltat) > mod(time,plot_deltat)) ) then

              if (istep <= 99999) then
                 write(unit=plot_index,fmt='(i5.5)') istep
                 plot_file_name = trim(plot_base_name) // plot_index
              else
                 write(unit=plot_index6,fmt='(i6.6)') istep
                 plot_file_name = trim(plot_base_name) // plot_index6
              endif

              call make_plotfile(plot_file_name,la,U,plot_names,time,dx,write_pf_time)

              call write_job_info(plot_file_name, la, write_pf_time)

              last_plt_written = istep

           end if
        end if

        if (parallel_IOProcessor() .and. verbose .ge. 2) then
           call flush()
        end if

        ! have we reached the stop time?
        if (stop_time >= 0.d0) then
           if (time >= stop_time) then
              goto 999
           end if
        end if

     end do

999  continue
     if (istep > max_step) istep = max_step

     if ( chk_int > 0 .and. last_chk_written .ne. istep ) then

        if (istep <= 99999) then
           write(unit=check_index,fmt='(i5.5)') istep
           check_file_name = trim(check_base_name) // check_index
        else
           write(unit=check_index6,fmt='(i6.6)') istep
           check_file_name = trim(check_base_name) // check_index6
        endif
        
        call checkpoint_write(check_file_name, la, U, dt)
              
     end if

     if ( plot_int > 0 .and. last_plt_written .ne. istep ) then

        if (istep <= 99999) then
           write(unit=plot_index,fmt='(i5.5)') istep
           plot_file_name = trim(plot_base_name) // plot_index
        else
           write(unit=plot_index6,fmt='(i6.6)') istep
           plot_file_name = trim(plot_base_name) // plot_index6
        endif

        call make_plotfile(plot_file_name,la,U,plot_names,time,dx,write_pf_time)

        call write_job_info(plot_file_name, la, write_pf_time)
     end if
  end if


  call destroy(U)

  call destroy(la)

  call chemistry_close()

  call runtime_close()

  deallocate(plot_names)
  deallocate(dx)

end subroutine smc
