program main

  use BoxLib
  use layout_module
  use bl_prof_module
  use multifab_module
  use probin_module, only: bl_prof_name

  implicit none

  real(dp_t) :: r1, r2

!  call boxlib_initialize(thread_support_level=MPI_THREAD_FUNNELED)
  call boxlib_initialize()

  r1 = parallel_wtime()

  call bl_prof_initialize(on = .true.)

  call smc()

  call layout_flush_copyassoc_cache()

  call bl_prof_glean(bl_prof_name)
  call bl_prof_finalize()

  r2 = parallel_wtime() - r1

  call parallel_reduce(r1, r2, MPI_MAX, proc = parallel_IOProcessorNode())

  if ( parallel_IOProcessor() ) then
     print*, 'MEMORY STATS AT END OF PROGRAM'
     print*, ' '
  end if
  call print(multifab_mem_stats(),    "    multifab")
  call print(fab_mem_stats(),         "         fab")
  call print(boxarray_mem_stats(),    "    boxarray")
  call print(layout_mem_stats(),      "      layout")
  call print(boxassoc_mem_stats(),    "    boxassoc")
!  call print(fgassoc_mem_stats(),     "     fgassoc")
!  call print(syncassoc_mem_stats(),   "   syncassoc")
  call print(copyassoc_mem_stats(),   "   copyassoc")
!  call print(fluxassoc_mem_stats(),   "   fluxassoc")

  if (parallel_IOProcessor()) then
     print*, ' '
     print*, 'Total Run Time = ', r1
     call flush(6)
  end if

  call boxlib_finalize()

end program main
