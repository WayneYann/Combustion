subroutine smc()

  use advance_module
  use cputime_module, only: start_cputime_clock
  use initialize_module
  use layout_module
  use multifab_module
  use omp_module
  use probin_module
  use runtime_init_module
  use time_module
  use variables

  implicit none

  integer :: dm, numcell
  integer :: init_step, istep
  integer :: last_plt_written,last_chk_written
  real(dp_t) :: dt

  type(layout) :: la

  real(dp_t)  , pointer     :: dx(:)

  type(multifab) :: U

  ! keep track of cputime
  call start_cputime_clock()

  last_plt_written = -1
  last_chk_written = -1

  call runtime_init()

  call init_variables()
  call init_plot_variables()

!  allocate(plot_names(n_plot_comps))
!  call get_plot_names(plot_names)

  if (restart >= 0) then

     ! call initialize_from_restart
     ! ...

  else 

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

  ! xxxxx write plotfile


  init_step = 1  ! xxxx if restart set it to restart + 1

  if ( (max_step >= init_step) .and. (time < stop_time .or. stop_time < 0.d0) ) then

     do istep = init_step, max_step

        if ( verbose .ge. 1 ) then
           if ( parallel_IOProcessor() ) then
              print *, 'MEMORY STATS AT START OF TIMESTEP ', istep
              print*, ' '
           end if
           call print(multifab_mem_stats(),    "    multifab")
           call print(fab_mem_stats(),         "         fab")
           call print(boxarray_mem_stats(),    "    boxarray")
           call print(layout_mem_stats(),      "      layout")
           call print(boxassoc_mem_stats(),    "    boxassoc")
           call print(fgassoc_mem_stats(),     "     fgassoc")
           call print(syncassoc_mem_stats(),   "   syncassoc")
           call print(copyassoc_mem_stats(),   "   copyassoc")
           call print(fluxassoc_mem_stats(),   "   fluxassoc")
           if ( parallel_IOProcessor() ) print*, ''
        end if

        if (parallel_IOProcessor()) then
           print*,'Advancing time step',istep,'time = ',time
        end if

        call advance(U,dt,dx)

        ! write checkpoint? plotfile?

     end do

  end if


  call destroy(U)

  call destroy(la)

  call runtime_close()

!  deallocate(n_plot_comps)
  deallocate(dx)

end subroutine smc
