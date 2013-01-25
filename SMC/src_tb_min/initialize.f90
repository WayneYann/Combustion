module initialize_module

  use multifab_module
  use variables_module
  use threadbox_module

  implicit none

  private
  public :: initialize_from_restart, initialize_from_scratch

contains

  subroutine initialize_from_restart(dirname,la,dt,courno,dx,U)
    use checkpoint_module, only : checkpoint_read
    use probin_module, only: n_cellx, n_celly, n_cellz, prob_lo, prob_hi, dm_in, &
         max_grid_size, change_max_grid_size, pmask, t_trylayout, overlap_testing
    use derivative_stencil_module, only : stencil_ng

    character(len=*), intent(in) :: dirname
    type(layout),intent(inout) :: la
    real(dp_t), intent(out) :: dt,courno
    real(dp_t), pointer :: dx(:)
    type(multifab), intent(inout) :: U

    ! local
    type(multifab), pointer :: chkdata(:)
    type(layout) :: lachk
    type(boxarray) :: bachk, ba
    type(box) :: bx
    integer :: dm, ncell(3), idim, lo(3), hi(3), ng
    real(dp_t) :: prob_lo_chk(3), prob_hi_chk(3)

    call checkpoint_read(chkdata, dirname, dt, courno, prob_lo_chk, prob_hi_chk, ncell) 

    dm = chkdata(1)%dim
    if (dm .ne. dm_in) then
       call bl_error("Inconsistent dimensionality in checkpoint file and probin")
    end if

    if (ncomp(chkdata(1)) .ne. ncons) then
       call bl_error("Inconsistent number of components in checkpoint file and code")
    end if
    
    if (ncell(1) .ne. n_cellx) then
       call bl_error("Inconsistent n_cellx in checkpoint file and probin")
    end if
    if (ncell(2) .ne. n_celly) then
       call bl_error("Inconsistent n_celly in checkpoint file and probin")
    end if
    if (dm .eq. 3) then
       if (ncell(3) .ne. n_cellz) then
          call bl_error("Inconsistent n_cellz in checkpoint file and probin")
       end if
    end if

    do idim=1,dm
       if (prob_lo(idim) .ne. prob_lo_chk(idim) .or. prob_hi(idim) .ne. prob_hi_chk(idim)) then
          call bl_error("Inconsistent prob_lo or prob_hi in checkpoint file and probin")
       end if
    end do

    allocate(dx(dm))
    dx(1) = (prob_hi(1)-prob_lo(1)) / n_cellx
    dx(2) = (prob_hi(2)-prob_lo(2)) / n_celly
    if (dm > 2) then
       dx(3) = (prob_hi(3)-prob_lo(3)) / n_cellz
    end if

    ng = stencil_ng

    lachk = get_layout(chkdata(1))

    bachk = get_boxarray(chkdata(1))

    if (change_max_grid_size) then
       lo = 0
       hi(1) = n_cellx-1
       hi(2) = n_celly-1
       if (dm > 2) then
          hi(3) = n_cellz - 1
       end if
       bx = make_box(lo,hi)
       call boxarray_build_bx(ba,bx)
       call boxarray_maxsize(ba,max_grid_size)
    else
       call boxarray_build_copy(ba, bachk)
    end if

    call layout_build_ba(la,ba,boxarray_bbox(ba),pmask=pmask)

    call build_threadbox(la,ng)

    call multifab_build(U,la,ncons,ng)
    call multifab_copy_c(U,1,chkdata(1),1,ncons)

    call destroy(chkdata(1))
    call destroy(lachk)

    if (parallel_nprocs() > 1 .and. t_trylayout > 0.d0) then
       call build_better_layout(U,la,ba,boxarray_bbox(ba),pmask,ncons,ng,t_trylayout)
    end if

    call destroy(ba)

    if (overlap_testing) then
       call overlap_vs_nooverlap(U)
    end if

  end subroutine initialize_from_restart


  subroutine initialize_from_scratch(la,dt,courno,dx,U)

    use init_data_module, only : init_data
    use time_module, only : time

    use probin_module, only: n_cellx, n_celly, n_cellz, prob_lo, prob_hi, dm_in, &
         max_grid_size, pmask, t_trylayout, overlap_testing
    use derivative_stencil_module, only : stencil_ng

    type(layout),intent(inout) :: la
    real(dp_t), intent(inout) :: dt,courno
    real(dp_t), pointer :: dx(:)
    type(multifab), intent(inout) :: U

    ! local
    integer :: lo(dm_in), hi(dm_in), dm, ng
    type(box)          :: bx
    type(boxarray)     :: ba

    time = ZERO
    dt   = 1.d20
    courno = -1.d20

    dm = dm_in
    lo = 0
    hi(1) = n_cellx-1
    hi(2) = n_celly-1
    if (dm > 2) then
       hi(3) = n_cellz - 1
    end if

    allocate(dx(dm))
    dx(1) = (prob_hi(1)-prob_lo(1)) / n_cellx
    dx(2) = (prob_hi(2)-prob_lo(2)) / n_celly
    if (dm > 2) then
       dx(3) = (prob_hi(3)-prob_lo(3)) / n_cellz
    end if

    ng = stencil_ng

    bx = make_box(lo,hi)
    
    call boxarray_build_bx(ba,bx)
    call boxarray_maxsize(ba,max_grid_size)
    call layout_build_ba(la,ba,boxarray_bbox(ba),pmask=pmask)

    call build_threadbox(la,ng)

    call multifab_build(U,la,ncons,ng)
    call tb_multifab_setval(U, 0.d0, ALL=.true.)
  
    call init_data(U,dx,prob_lo,prob_hi)

    if (parallel_nprocs() > 1 .and. t_trylayout > 0.d0) then
       call build_better_layout(U,la,ba,boxarray_bbox(ba),pmask,ncons,ng,t_trylayout)
    end if

    call destroy(ba)

    if (overlap_testing) then
       call overlap_vs_nooverlap(U)
    end if

  end subroutine initialize_from_scratch

  
  subroutine build_better_layout(U, la, ba, pd, pmask, nc, ng, ttry)
    use probin_module, only : verbose, overlap_comm_comp, overlap_in_trial
    use advance_module, only : overlapped_part
    use bl_prof_module
    type(multifab),intent(inout) :: U
    type(layout),intent(inout) :: la
    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: pd
    logical, intent(in) :: pmask(:)
    integer, intent(in) :: nc, ng 
    real(dp_t), intent(in) :: ttry

    integer :: i, nb, ntry
    integer, dimension(:), allocatable :: dmap0, dmapt, dmapbest
    real(dp_t) :: t1, t2, timespent, tcom1, tcom2, tbest, tdefault
    type(layout) :: lat
    type(multifab) :: Ut
    type(mf_fb_data) :: Ut_fb_data

    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp_t), allocatable :: r(:)

    logical :: bp_state

    bp_state = bl_prof_get_state()
    call bl_prof_set_state(.false.) ! turn profiler off temporarily

    nb = nboxes(la)

    allocate(dmap0(nb))
    allocate(dmapt(nb))
    allocate(dmapbest(nb))
    allocate(r(nb))

    dmap0 = la%lap%prc(:)

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    do i=1,seed_size
       seed(i) = 134527 + 59*i
    end do
    call random_seed(put=seed)

    ntry = 0
    timespent = 0.0_dp_t
    tbest = 1.0d50
    do while (timespent < ttry)

       call parallel_barrier()
       t1 = parallel_wtime()

       ntry = ntry + 1

       if (ntry .eq. 1) then
          dmapt = dmap0
       else
          call random_number(r)
          call Fisher_Yates_shuffle(dmap0, dmapt, r, nb)
       end if

       call layout_build_ba(lat, ba, pd, pmask=pmask, &
            mapping = LA_EXPLICIT, explicit_mapping = dmapt)

       call build_threadbox(lat,ng)
       
       call multifab_build(Ut, lat, nc, ng)
       call multifab_copy_c(Ut,1,U,1,nc)

       call parallel_barrier()
       tcom1 = parallel_wtime()

       call multifab_fill_boundary_nowait(Ut, Ut_fb_data)
       if (overlap_comm_comp .and. overlap_in_trial) then
          call overlapped_part(Ut, Ut_fb_data)
       end if
       call multifab_fill_boundary_finish(Ut, Ut_fb_data)

       call parallel_barrier()
       tcom2 = parallel_wtime() - tcom1

       if (tcom2 < tbest) then
          tbest = tcom2
          dmapbest = dmapt
       end if

       if (ntry .eq. 1) then
          tdefault = tcom2
       end if

       call destroy(Ut)
       call destroy(lat)

       call parallel_barrier()
       t2 = parallel_wtime() - t1
       timespent = timespent + t2
       call parallel_bcast(timespent)
    end do

    call parallel_bcast(dmapbest)

    call layout_build_ba(lat, ba, pd, pmask=pmask, &
         mapping = LA_EXPLICIT, explicit_mapping = dmapbest)

    call build_threadbox(lat,ng)

    call multifab_build(Ut, lat, nc, ng)
    call multifab_copy_c(Ut,1,U,1,nc)

    call destroy(U)
    call destroy(la)

    la = lat
    U = Ut

    deallocate(dmap0,dmapt,dmapbest,seed,r)

    if (verbose > 0 .and. parallel_IOProcessor()) then
       print *, 'Tried', ntry, 'layouts for filling multifab boundaries.'
       if (overlap_comm_comp .and. overlap_in_trial) then
          print *, 'with overlapped computation'
       end if
       print *, '   The average time in second is', timespent/ntry
       print *, '   Using the default layout, it is', tdefault
       print *, '   The best time is', tbest
    end if

    call bl_prof_set_state(bp_state)

    contains

      subroutine Fisher_Yates_shuffle(s, a, r, n)
        integer, intent(in) :: n
        integer, intent(in) :: s(n)
        integer, intent(out) :: a(n)
        real(dp_t), intent(in) :: r(n)
        integer :: i, j
        a(1) = s(1)
        do i = 2, n
           j = int(r(i)*i)+1
           a(i) = a(j)
           a(j) = s(i)
        end do
      end subroutine Fisher_Yates_shuffle

  end subroutine build_better_layout


  subroutine overlap_vs_nooverlap(U)
    use probin_module, only : overlap_comm_comp
    use advance_module, only : overlapped_part
    use bl_prof_module
    type(multifab), intent(inout) :: U

    type(mf_fb_data) :: mfd1, mfd2
    double precision :: t1, t2, t3, t4
    logical :: bp_state

    bp_state = bl_prof_get_state()
    call bl_prof_set_state(.false.) ! turn profiler off temporarily

    call parallel_barrier()
    t1 = parallel_wtime()

    call multifab_fill_boundary_nowait(U, mfd1)
    call overlapped_part(U, mfd1)
    call multifab_fill_boundary_finish(U, mfd1)

    call parallel_barrier()
    t2 = parallel_wtime()

    call multifab_fill_boundary_nowait(U, mfd2)
    call multifab_fill_boundary_finish(U, mfd2)
    call overlapped_part(U, mfd2)

    call parallel_barrier()
    t3 = parallel_wtime()

    call overlapped_part(U, mfd2)

    call parallel_barrier()
    t4 = parallel_wtime()

    if (parallel_IOProcessor()) then
       print *, ''
       print *, 'Testing communication and computation overlapping:'
       print *, '   Computation time:', t4-t3
       print *, '   Communication time with overlapping:', t2-t1 - (t4-t3)
       print *, '   Communication time w/o  overlapping:', t3-t2 - (t4-t3)
       if (overlap_comm_comp .and. (t2-t1) > (t3-t2)) then
          print *, 'Turn off overlapping because the test shows it does not work.'
          overlap_comm_comp = .false.
       end if
       print *, ''
    end if

    call parallel_bcast(overlap_comm_comp)

    call bl_prof_set_state(bp_state)

  end subroutine overlap_vs_nooverlap

end module initialize_module
