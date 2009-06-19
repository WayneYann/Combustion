program main

  use BoxLib
  use parallel
  use layout_module
  use bl_prof_module
  use convergence

  implicit none

  real(dp_t) :: r1, r2
  integer    :: loc(1,4) = (/9032,9064,9128,9256/) 
  integer    :: files(1,2) = (/10,11/)        
  integer    :: i

  call boxlib_initialize()

  r1 = parallel_wtime()

  call bl_prof_initialize(on = .true.)

  !call layout_set_verbosity(1)

  open(10, FILE='l1_strang.dat', STATUS='NEW')
  open(11, FILE='l2_strang.dat', STATUS='NEW')

  do i = 1, 1
     call conv(loc(i,:),files(i,:))
  end do

  do i = 10,11
     close(i)
  end do

  call bl_prof_glean("bl_prof_res")

  call bl_prof_finalize()

  r2 = parallel_wtime() - r1

  call parallel_reduce(r1, r2, MPI_MAX, proc = parallel_IOProcessorNode())

  if (parallel_IOProcessor()) print*, 'Run Time = ', r1

  call boxlib_finalize()

end program main
