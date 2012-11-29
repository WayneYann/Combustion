module checkpoint_module
  
  use bl_types, only: dp_t
  use bc_module
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
         bcx_lo,bcx_hi,bcy_lo,bcy_hi,bcz_lo,bcz_hi, &
         n_cellx, n_celly, n_cellz
    use time_module, only: time
    use cputime_module, only: get_cputime
    use nscbc_module

    character(len=*), intent(in) :: dirname
    type(multifab)  , intent(in) :: U
    real(kind=dp_t) , intent(in) :: dt, courno

    type(multifab) :: chkdata(1), chkxlo(1), chkxhi(1), &
         chkylo(1), chkyhi(1), chkzlo(1), chkzhi(1)
    character(len=256) :: sd_name, xlo_name, ylo_name, zlo_name, &
         xhi_name, yhi_name, zhi_name
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

    if (bcx_lo .eq. INLET) then
       write(unit=xlo_name, fmt='(a,"/qinxlo")') trim(dirname)
       call build_chk_for_qin(chkxlo(1), qin_xlo)
       call fabio_ml_multifab_write_d(chkxlo, rrs, xlo_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)
       call destroy(chkxlo(1))
    end if

    if (bcy_lo .eq. INLET) then
       write(unit=ylo_name, fmt='(a,"/qinylo")') trim(dirname)
       call build_chk_for_qin(chkylo(1), qin_ylo)
       call fabio_ml_multifab_write_d(chkylo, rrs, ylo_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)
       call destroy(chkylo(1))
    end if

    if (bcz_lo .eq. INLET) then
       write(unit=zlo_name, fmt='(a,"/qinzlo")') trim(dirname)
       call build_chk_for_qin(chkzlo(1), qin_zlo)
       call fabio_ml_multifab_write_d(chkzlo, rrs, zlo_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)
       call destroy(chkzlo(1))
    end if

    if (bcx_hi .eq. INLET) then
       write(unit=xhi_name, fmt='(a,"/qinxhi")') trim(dirname)
       call build_chk_for_qin(chkxhi(1), qin_xhi)
       call fabio_ml_multifab_write_d(chkxhi, rrs, xhi_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)
       call destroy(chkxhi(1))
    end if

    if (bcy_hi .eq. INLET) then
       write(unit=yhi_name, fmt='(a,"/qinyhi")') trim(dirname)
       call build_chk_for_qin(chkyhi(1), qin_yhi)
       call fabio_ml_multifab_write_d(chkyhi, rrs, yhi_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)
       call destroy(chkyhi(1))
    end if

    if (bcz_hi .eq. INLET) then
       write(unit=zhi_name, fmt='(a,"/qinzhi")') trim(dirname)
       call build_chk_for_qin(chkzhi(1), qin_zhi)
       call fabio_ml_multifab_write_d(chkzhi, rrs, zhi_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)
       call destroy(chkzhi(1))
    end if

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
       write(unit=un,fmt=1001) bcx_lo, bcy_lo, bcz_lo
       write(unit=un,fmt=1001) bcx_hi, bcy_hi, bcz_hi
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


  subroutine checkpoint_read(chkdata, chkxlo, chkxhi, chkylo, chkyhi, chkzlo, chkzhi, &
       dirname, dt, courno, prob_lo_chk, prob_hi_chk, bc_lo_chk, bc_hi_chk, ncell)

    use bl_IO_module, only: unit_new
    use fabio_module, only: fabio_ml_multifab_read_d
    use bl_prof_module, only: bl_prof_timer, build, destroy
    use time_module, only: time

    type(multifab  ), pointer :: chkdata(:), chkxlo(:), chkxhi(:), &
         chkylo(:), chkyhi(:), chkzlo(:), chkzhi(:)
    character(len=*), intent(in   ) :: dirname
    real(kind=dp_t) , intent(  out) :: dt, courno
    real(kind=dp_t), intent(out) :: prob_lo_chk(3), prob_hi_chk(3)
    integer, intent(out) :: bc_lo_chk(3), bc_hi_chk(3)
    integer, intent(out) :: ncell(3)

    ! local
    integer            :: un
    character(len=256) :: header, sd_name, xlo_name, ylo_name, zlo_name, &
         xhi_name, yhi_name, zhi_name

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
    read(unit=un,fmt=*)   bc_lo_chk(1),   bc_lo_chk(2),   bc_lo_chk(3)
    read(unit=un,fmt=*)   bc_hi_chk(1),   bc_hi_chk(2),   bc_hi_chk(3)
    read(unit=un,fmt=*) ncell(1), ncell(2), ncell(3)
    close(un)

    !   Read the state data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_read_d(chkdata, sd_name)

    if (bc_lo_chk(1) .eq. INLET) then
       write(unit=xlo_name, fmt='(a,"/qinxlo")') trim(dirname)
       call fabio_ml_multifab_read_d(chkxlo, xlo_name)       
    end if

    if (bc_lo_chk(2) .eq. INLET) then
       write(unit=ylo_name, fmt='(a,"/qinylo")') trim(dirname)
       call fabio_ml_multifab_read_d(chkylo, ylo_name)       
    end if

    if (bc_lo_chk(3) .eq. INLET) then
       write(unit=zlo_name, fmt='(a,"/qinzlo")') trim(dirname)
       call fabio_ml_multifab_read_d(chkzlo, zlo_name)       
    end if

    if (bc_hi_chk(1) .eq. INLET) then
       write(unit=xhi_name, fmt='(a,"/qinxhi")') trim(dirname)
       call fabio_ml_multifab_read_d(chkxhi, xhi_name)       
    end if

    if (bc_hi_chk(2) .eq. INLET) then
       write(unit=yhi_name, fmt='(a,"/qinyhi")') trim(dirname)
       call fabio_ml_multifab_read_d(chkyhi, yhi_name)       
    end if

    if (bc_hi_chk(3) .eq. INLET) then
       write(unit=zhi_name, fmt='(a,"/qinzhi")') trim(dirname)
       call fabio_ml_multifab_read_d(chkzhi, zhi_name)       
    end if

    call destroy(bpt)

  end subroutine checkpoint_read


  subroutine build_chk_for_qin(chk, qin)
    use physbndry_reg_module
    type(multifab), intent(inout) :: chk
    type(physbndry_reg), intent(in) :: qin

    integer :: i,j,k,n
    integer :: lo(3), hi(3)
    double precision, pointer, dimension(:,:,:,:) :: qp, cp
    
    call build(chk, qin%la, qin%nc, ng=0, stencil=.false.)
    
    do n=1,nfabs(chk)

       cp => dataptr(chk,n)

       if (isValid(qin,n)) then
          lo = lwb(get_box(chk,n))
          hi = upb(get_box(chk,n))

          qp => dataptr(qin%data, n)

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   cp(i,j,k,:) = qp(:,i,j,k)
                end do
             end do
          end do

       else

          cp = 0.d0  ! set to 0 because fabio doesn't nans

       end if

    end do

  end subroutine build_chk_for_qin

end module checkpoint_module
