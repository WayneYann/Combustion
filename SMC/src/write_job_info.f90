subroutine write_job_info(dirname, la, write_pf_time)

  ! write out some basic information about the way the job was run
  ! to a file called job_info in the directory dir_name.  Usually
  ! dir_name will be the name of the checkpoint or plotfile toplevel
  ! directory

  use bl_types
  use parallel
  use probin_module, only: job_name, inputs_file_used, bcx_lo, bcy_lo, bcz_lo, &
       bcx_hi, bcy_hi, bcz_hi
  use runtime_init_module, only: probin
  use bl_system_module, only: BL_CWD_SIZE, get_cwd 
  use bc_module
  use layout_module
  use build_info_module, only: build_date, build_dir, build_machine, boxlib_dir, &
                               NUM_MODULES, modules, FCOMP, FCOMP_version, &
                               f90_compile_line, f_compile_line, &
                               C_compile_line, link_line, &
                               source_git_hash, boxlib_git_hash
  use omp_module
  use cputime_module, only: get_cputime

  implicit none

  character (len=*), intent(in) :: dirname
  type(layout), intent(in) :: la
  real(kind=dp_t)  , intent(in) :: write_pf_time

  character (len=256) :: out_name
  character (len=16) :: date_in, time_in
  integer, dimension(8) :: values
  character (len=BL_CWD_SIZE) :: cwd

  integer :: i, dm
  type(box) :: pbx

  call date_and_time(date_in, time_in, VALUES=values)
  call get_cwd(cwd)

  dm = get_dim(la)
  pbx = get_pd(la)

  out_name = trim(dirname) // "/job_info"

999  format(79('='))
1001 format(a,a)
1002 format(a,i6)
1003 format(a,i4.4,'-',i2.2,'-',i2.2)
1004 format(a,i2.2,':',i2.2,':',i2.2)
1005 format(a,g20.10)

  if (parallel_IOProcessor()) then
     open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")
     
     write (99,999)
     write (99,*) "Job Information"
     write (99,999)
     write (99,1001) "job name:    ", trim(job_name)
     write (99,1001) "inputs file: ", trim(inputs_file_used)
     write (99,*) " "     
     write (99,1002) "number of MPI processes ", parallel_nprocs()
     write (99,1002) "number of threads       ", omp_get_max_threads()
     write (99,*) " "
     write (99,1005) "CPU time used since start of simulation (CPU-hours) ", get_cputime()/3600.0_dp_t

     write (99,*) " "
     write (99,*) " "

     write (99,999)
     write (99,*) "Plotfile Information"
     write (99,999)
     write (99,1003) "output date:              ", values(1), values(2), values(3)
     write (99,1004) "output time:              ", values(5), values(6), values(7)
     write (99,1001) "output dir:               ", trim(cwd)
     write (99,1005) "time to write plotfile (s)", write_pf_time

     write (99,*) " "
     write (99,*) " "


     write (99,999)
     write (99,*) "Build Information"
     write (99,999)
     write (99,1001) "build date:    ", trim(build_date)
     write (99,1001) "build machine: ", trim(build_machine)
     write (99,1001) "build dir:     ", trim(build_dir)
     write (99,1001) "BoxLib dir:    ", trim(boxlib_dir)
     write (99,*) " "
     write (99,1001) "Combustion git hash: ", trim(source_git_hash)
     write (99,1001) "BoxLib     git hash: ", trim(boxlib_git_hash)
     write (99,*) " "
     write (99,1001) "modules used:  ", " "
     do i=1, NUM_MODULES
        write (99,1001) "  ", trim(modules(i))
     enddo
     write (99,*) " "
     write (99,1001) "FCOMP:            ", trim(FCOMP)
     write (99,1001) "FCOMP version:    ", trim(FCOMP_version)
     write (99,*) " "
     write (99,1001) "F90 compile line: ", trim(f90_compile_line)
     write (99,*) " "
     write (99,1001) "F77 compile line: ", trim(f_compile_line)
     write (99,*) " "     
     write (99,1001) "C compile line:   ", trim(C_compile_line)
     write (99,*) " "
     write (99,1001) "linker line:      ", trim(link_line)

     write (99,*) " "
     write (99,*) " "


     write (99,999)
     write (99,*) "Grid Information"
     write (99,999)
     write (99,*) "   number of boxes = ", nboxes(la)
     write (99,*) "   maximum zones   = ", (extent(pbx,i),i=1,dm)

     write (99,*) " "
     write (99,*) "Boundary Conditions"     
     write (99,*) "  -x: ", bc_integer_to_string(bcx_lo)
     write (99,*) "  +x: ", bc_integer_to_string(bcx_hi)

     if (dm >= 2) then
        write (99,*) " "
        write (99,*) "  -y: ", bc_integer_to_string(bcy_lo)
        write (99,*) "  +y: ", bc_integer_to_string(bcy_hi)
     endif

     if (dm == 3) then
        write (99,*) " "
        write (99,*) "  -z: ", bc_integer_to_string(bcz_lo)
        write (99,*) "  +z: ", bc_integer_to_string(bcz_hi)
     endif

     write (99,*) " "
     write (99,*) " "


!     write (99,999)
!     write (99,*) "Chemistry Model:", 
!
!     write (99,*) " "
!     write (99,*) " "


     write (99,999)
     write (99,*) "Runtime Parameter Information"
     write (99,999)
     write (99,nml=probin)
     close(99)
  endif

end subroutine write_job_info
  
