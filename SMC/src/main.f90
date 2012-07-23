program main

  use BoxLib
  use parallel
  use layout_module
  use bl_prof_module
  use multifab_module

  implicit none

  call boxlib_initialize()


  call boxlib_finalize()

end program main
