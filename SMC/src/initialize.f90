module initialize_module

  use multifab_module
  use bc_module

  use smc_bc_module
  use nscbc_module
  use physbndry_reg_module
  use variables_module

  implicit none

  private
  public :: initialize_from_restart, initialize_from_scratch

contains

  subroutine initialize_from_restart(dirname,la,dt,dx,U)
    use checkpoint_module, only : checkpoint_read
    use probin_module, only: n_cellx, n_celly, n_cellz, prob_lo, prob_hi, &
         bcx_lo, bcy_lo, bcz_lo, bcx_hi, bcy_hi, bcz_hi, &
         dm_in, max_grid_size, change_max_grid_size, pmask, t_trylayout
    use derivative_stencil_module, only : stencil_ng

    character(len=*), intent(in) :: dirname
    type(layout),intent(inout) :: la
    real(dp_t), intent(out) :: dt
    real(dp_t), pointer :: dx(:)
    type(multifab), intent(inout) :: U

    ! local
    type(multifab), pointer :: chkdata(:), chkxlo(:), chkxhi(:), &
         chkylo(:), chkyhi(:), chkzlo(:), chkzhi(:)
    type(physbndry_reg) :: chqin_xlo, chqin_xhi, chqin_ylo, chqin_yhi, chqin_zlo, chqin_zhi
    type(layout) :: lachk
    type(boxarray) :: bachk, ba
    type(box) :: bx
    integer :: dm, ncell(3), idim, lo(3), hi(3), ng
    real(dp_t) :: prob_lo_chk(3), prob_hi_chk(3)
    integer :: bc_lo_chk(3), bc_hi_chk(3)

    call checkpoint_read(chkdata, chkxlo, chkxhi, chkylo, chkyhi, chkzlo, chkzhi, &
         dirname, dt, prob_lo_chk, prob_hi_chk, bc_lo_chk, bc_hi_chk, ncell) 

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

    if (bcx_lo .ne. bc_lo_chk(1)) then
       call bl_error("Inconsistent bcx_lo in checkpoint file and probin")
    end if
    if (bcy_lo .ne. bc_lo_chk(2)) then
       call bl_error("Inconsistent bcy_lo in checkpoint file and probin")
    end if
    if (bcz_lo .ne. bc_lo_chk(3)) then
       call bl_error("Inconsistent bcz_lo in checkpoint file and probin")
    end if
    if (bcx_hi .ne. bc_hi_chk(1)) then
       call bl_error("Inconsistent bcx_hi in checkpoint file and probin")
    end if
    if (bcy_hi .ne. bc_hi_chk(2)) then
       call bl_error("Inconsistent bcy_hi in checkpoint file and probin")
    end if
    if (bcz_hi .ne. bc_hi_chk(3)) then
       call bl_error("Inconsistent bcz_hi in checkpoint file and probin")
    end if

    allocate(dx(dm))
    dx(1) = (prob_hi(1)-prob_lo(1)) / n_cellx
    dx(2) = (prob_hi(2)-prob_lo(2)) / n_celly
    if (dm > 2) then
       dx(3) = (prob_hi(3)-prob_lo(3)) / n_cellz
    end if

    ng = stencil_ng

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

    if (parallel_nprocs() > 1 .and. t_trylayout > 0.d0) then
       call build_better_layout(la,ba,boxarray_bbox(ba),pmask,ncons,ng,t_trylayout)
    end if

    call destroy(ba)

    call multifab_build(U,la,ncons,ng)
    call setval(U, 0.d0, .true.)
    call multifab_copy_c(U,1,chkdata(1),1,ncons)

    lachk = get_layout(chkdata(1))
    call destroy(chkdata(1))
    call destroy(lachk)

    call smc_bc_init(la, U)

    call nscbc_init()
    call nscbc_build_registers(la)

    if (bcx_lo .eq. INLET) then
       call physbndry_reg_build(chqin_xlo,la,qin_xlo%nc,qin_xlo%idim,qin_xlo%iface,.false.)
       call multifab_copy_c(chqin_xlo%data,1,chkxlo(1),1,chqin_xlo%nc)
       call physbndry_reg_copy(chqin_xlo, qin_xlo)
       call physbndry_reg_destroy(chqin_xlo)
       lachk = get_layout(chkxlo(1))
       call destroy(chkxlo(1))
       call destroy(lachk)
    end if

    if (bcy_lo .eq. INLET) then
       call physbndry_reg_build(chqin_ylo,la,qin_ylo%nc,qin_ylo%idim,qin_ylo%iface,.false.)
       call multifab_copy_c(chqin_ylo%data,1,chkylo(1),1,chqin_ylo%nc)
       call physbndry_reg_copy(chqin_ylo, qin_ylo)
       call physbndry_reg_destroy(chqin_ylo)
       lachk = get_layout(chkylo(1))
       call destroy(chkylo(1))
       call destroy(lachk)
    end if

    if (bcz_lo .eq. INLET) then
       call physbndry_reg_build(chqin_zlo,la,qin_zlo%nc,qin_zlo%idim,qin_zlo%iface,.false.)
       call multifab_copy_c(chqin_zlo%data,1,chkzlo(1),1,chqin_zlo%nc)
       call physbndry_reg_copy(chqin_zlo, qin_zlo)
       call physbndry_reg_destroy(chqin_zlo)
       lachk = get_layout(chkzlo(1))
       call destroy(chkzlo(1))
       call destroy(lachk)
    end if

    if (bcx_hi .eq. INLET) then
       call physbndry_reg_build(chqin_xhi,la,qin_xhi%nc,qin_xhi%idim,qin_xhi%iface,.false.)
       call multifab_copy_c(chqin_xhi%data,1,chkxhi(1),1,chqin_xhi%nc)
       call physbndry_reg_copy(chqin_xhi, qin_xhi)
       call physbndry_reg_destroy(chqin_xhi)
       lachk = get_layout(chkxhi(1))
       call destroy(chkxhi(1))
       call destroy(lachk)
    end if

    if (bcy_hi .eq. INLET) then
       call physbndry_reg_build(chqin_yhi,la,qin_yhi%nc,qin_yhi%idim,qin_yhi%iface,.false.)
       call multifab_copy_c(chqin_yhi%data,1,chkyhi(1),1,chqin_yhi%nc)
       call physbndry_reg_copy(chqin_yhi, qin_yhi)
       call physbndry_reg_destroy(chqin_yhi)
       lachk = get_layout(chkyhi(1))
       call destroy(chkyhi(1))
       call destroy(lachk)
    end if

    if (bcz_hi .eq. INLET) then
       call physbndry_reg_build(chqin_zhi,la,qin_zhi%nc,qin_zhi%idim,qin_zhi%iface,.false.)
       call multifab_copy_c(chqin_zhi%data,1,chkzhi(1),1,chqin_zhi%nc)
       call physbndry_reg_copy(chqin_zhi, qin_zhi)
       call physbndry_reg_destroy(chqin_zhi)
       lachk = get_layout(chkzhi(1))
       call destroy(chkzhi(1))
       call destroy(lachk)
    end if

  end subroutine initialize_from_restart


  subroutine initialize_from_scratch(la,dt,dx,U)

    use init_data_module, only : init_data
    use time_module, only : time

    use probin_module, only: n_cellx, n_celly, n_cellz, prob_lo, prob_hi, dm_in, &
         max_grid_size, pmask, t_trylayout
    use derivative_stencil_module, only : stencil_ng

    type(layout),intent(inout) :: la
    real(dp_t), intent(inout) :: dt
    real(dp_t), pointer :: dx(:)
    type(multifab), intent(inout) :: U

    ! local
    integer :: lo(dm_in), hi(dm_in), dm, ng
    type(box)          :: bx
    type(boxarray)     :: ba

    time = ZERO
    dt   = 1.d20

    dm = dm_in
    lo = 0
    hi(1) = n_cellx-1
    hi(2) = n_celly-1
    if (dm > 2) then
       hi(3) = n_cellz - 1
    end if

    ng = stencil_ng

    bx = make_box(lo,hi)
    
    call boxarray_build_bx(ba,bx)
    call boxarray_maxsize(ba,max_grid_size)
    call layout_build_ba(la,ba,boxarray_bbox(ba),pmask=pmask)

    if (parallel_nprocs() > 1 .and. t_trylayout > 0.d0) then
       call build_better_layout(la,ba,boxarray_bbox(ba),pmask,ncons,ng,t_trylayout)
    end if

    call destroy(ba)

    allocate(dx(dm))
    dx(1) = (prob_hi(1)-prob_lo(1)) / n_cellx
    dx(2) = (prob_hi(2)-prob_lo(2)) / n_celly
    if (dm > 2) then
       dx(3) = (prob_hi(3)-prob_lo(3)) / n_cellz
    end if

    call multifab_build(U,la,ncons,ng)
    call setval(U, 0.d0, .true.)

    call init_data(U,dx,prob_lo,prob_hi)

    call smc_bc_init(la, U)

    call nscbc_init()
    call nscbc_build_registers(la)

    call nscbc_init_inlet_reg_from_scratch(U)

  end subroutine initialize_from_scratch


  subroutine build_better_layout(la, ba, pd, pmask, nc, ng, ttry)
    use probin_module, only : verbose
    type(layout),intent(inout) :: la
    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: pd
    logical, intent(in) :: pmask(:)
    integer, intent(in) :: nc, ng 
    real(dp_t), intent(in) :: ttry

    integer :: i, nb, ntry
    integer, dimension(:), allocatable :: dmap0, dmapt, dmapbest
    real(dp_t) :: t1, t2, tthis, timespent, tbest, tdefault
    type(layout) :: lat
    type(multifab) :: mf

    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp_t), allocatable :: r(:)

    nb = nboxes(la)

    allocate(dmap0(nb))
    allocate(dmapt(nb))
    allocate(dmapbest(nb))
    allocate(r(nb))

    dmap0 = la%lap%prc(:)
    call destroy(la)

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

       ntry = ntry + 1

       if (ntry .eq. 1) then
          dmapt = dmap0
       else
          call random_number(r)
          call Fisher_Yates_shuffle(dmap0, dmapt, r, nb)
       end if

       call layout_build_ba(lat, ba, pd, pmask=pmask, &
            mapping = LA_EXPLICIT, explicit_mapping = dmapt)
       
       call multifab_build(mf, lat, nc, ng)
       call setval(mf, 0.0_dp_t)

       call parallel_barrier()
       t1 = parallel_wtime()
       call multifab_fill_boundary(mf)
       call parallel_barrier()
       t2 = parallel_wtime()

       tthis = t2-t1
       timespent = timespent + tthis
       call parallel_bcast(timespent)

       if (tthis < tbest) then
          tbest = tthis
          dmapbest = dmapt
       end if

       if (ntry .eq. 1) then
          tdefault = tthis
       end if

       call destroy(mf)
       call destroy(lat)
    end do

    call parallel_bcast(dmapbest)

    call layout_build_ba(la, ba, pd, pmask=pmask, &
         mapping = LA_EXPLICIT, explicit_mapping = dmapbest)

    deallocate(dmap0,dmapt,dmapbest,seed,r)

    if (verbose > 0 .and. parallel_IOProcessor()) then
       print *, 'Tried', ntry, 'layouts for filling multifab boundaries.'
       print *, '   The average time in second is', timespent/ntry
       print *, '   Using the default layout, it is', tdefault
       print *, '   The best time is', tbest
    end if

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

end module initialize_module
