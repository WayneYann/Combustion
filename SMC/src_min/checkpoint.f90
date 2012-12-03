module checkpoint_module
  
  use bl_types, only: dp_t
  use multifab_module

  implicit none

  private

  public :: checkpoint_write, checkpoint_read

contains

  subroutine checkpoint_write(dirname, U, dt, courno)

    use parallel, only: parallel_IOProcessor, parallel_barrier
    use bl_IO_module, only: unit_new
    use fabio_module, only: fabio_mkdir, fabio_ml_multifab_write_d
    use bl_prof_module, only: bl_prof_timer, build, destroy
    use probin_module, only: verbose, nOutFiles, lUsingNFiles, &
         prob_lo_x, prob_lo_y, prob_lo_z, prob_hi_x, prob_hi_y, prob_hi_z, &
         n_cellx, n_celly, n_cellz
    use time_module, only: time
    use cputime_module, only: get_cputime

    character(len=*), intent(in) :: dirname
    type(multifab)  , intent(in) :: U
    real(kind=dp_t) , intent(in) :: dt, courno

    type(multifab) :: chkdata(1)
    character(len=256) :: sd_name
    integer :: un
    integer, dimension(0) :: rrs
    real(dp_t) :: writetime1, writetime2
    type(bl_prof_timer), save :: bpt
    
    call build(bpt, "checkpoint_write")

    chkdata(1) = U

    if ( parallel_IOProcessor() ) then
       call fabio_mkdir(dirname)
    end if

    call parallel_barrier()

    writetime1 = parallel_wtime()

    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_write_d(chkdata, rrs, sd_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) ''
       write(6,*) 'Writing state to checkpoint file ',trim(dirname)
    end if

    if (parallel_IOProcessor()) then
       un = unit_new()
       open(unit=un, &
            file = trim(dirname) // "/Header", &
            form = "formatted", access = "sequential", &
            status = "replace", action = "write")
       write(unit=un,fmt=1000) dt, courno
       write(unit=un,fmt=1000) time
       write(unit=un,fmt=1000) prob_lo_x, prob_lo_y, prob_lo_z
       write(unit=un,fmt=1000) prob_hi_x, prob_hi_y, prob_hi_z
       write(unit=un,fmt=1001) n_cellx, n_celly, n_cellz
       close(un)
    end if

    if (parallel_IOProcessor()) then
       un = unit_new()
       open(unit=un, file=trim(dirname) // "/CPUtime", &
            form="formatted", action="write", status="replace")
       write(unit=un,fmt=*) get_cputime()
       close(un)
    endif

    writetime2 = parallel_wtime() - writetime1
    call parallel_reduce(writetime1, writetime2, MPI_MAX, proc=parallel_IOProcessorNode())
    if (parallel_IOProcessor()) then
       print*,'Time to write checkpoint: ',writetime1,' seconds'
    end if

    call destroy(bpt)

1000 format(32(e30.20,1x))
1001 format(3(i0,1x))

  end subroutine checkpoint_write


  subroutine checkpoint_read(chkdata, dirname, dt, courno, prob_lo_chk, prob_hi_chk, ncell)

    use bl_IO_module, only: unit_new
    use fabio_module, only: fabio_ml_multifab_read_d
    use bl_prof_module, only: bl_prof_timer, build, destroy
    use time_module, only: time

    type(multifab  ), pointer :: chkdata(:)
    character(len=*), intent(in   ) :: dirname
    real(kind=dp_t) , intent(  out) :: dt, courno
    real(kind=dp_t), intent(out) :: prob_lo_chk(3), prob_hi_chk(3)
    integer, intent(out) :: ncell(3)

    ! local
    integer            :: un
    character(len=256) :: header, sd_name

    type(bl_prof_timer), save :: bpt

    call build(bpt, "checkpoint_read")

!   First read the header information
    header = "Header"
    un = unit_new()
    open(unit=un, &
         file = trim(dirname) // "/" // trim(header), &
         status = "old", &
         action = "read")
    read(unit=un,fmt=*) dt, courno
    read(unit=un,fmt=*) time
    read(unit=un,fmt=*) prob_lo_chk(1), prob_lo_chk(2), prob_lo_chk(3)
    read(unit=un,fmt=*) prob_hi_chk(1), prob_hi_chk(2), prob_hi_chk(3)
    read(unit=un,fmt=*) ncell(1), ncell(2), ncell(3)
    close(un)

    !   Read the state data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_read_d(chkdata, sd_name)

    call destroy(bpt)

  end subroutine checkpoint_read

end module checkpoint_module
