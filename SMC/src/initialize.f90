module initialize_module

  use multifab_module
  use variables_module

  implicit none

  private
  public :: initialize_from_restart, initialize_from_scratch

contains

  subroutine initialize_from_restart(dirname,la,dt,dx,U)
    use checkpoint_module, only : checkpoint_read
    use probin_module, only: n_cellx, n_celly, n_cellz, prob_lo, prob_hi, dm_in, &
         max_grid_size, change_max_grid_size, pmask
    use derivative_stencil_module, only : stencil_ng

    character(len=*), intent(in) :: dirname
    type(layout),intent(inout) :: la
    real(dp_t), intent(out) :: dt
    real(dp_t), pointer :: dx(:)
    type(multifab), intent(inout) :: U

    ! local
    type(multifab), pointer :: chkdata(:)
    type(layout) :: lachk
    type(boxarray) :: bachk, ba
    type(box) :: bx
    integer :: dm, ncell(3), idim, lo(3), hi(3), ng
    real(dp_t) :: prob_lo_chk(3), prob_hi_chk(3)

    call checkpoint_read(chkdata, dirname, dt, prob_lo_chk, prob_hi_chk, ncell) 

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
    call destroy(ba)

    ng = stencil_ng

    call multifab_build(U,la,ncons,ng)
    call multifab_copy_c(U,1,chkdata(1),1,ncons)

    call destroy(lachk)
    call destroy(chkdata(1))

  end subroutine initialize_from_restart


  subroutine initialize_from_scratch(la,dt,dx,U)

    use init_data_module, only : init_data
    use time_module, only : time

    use probin_module, only: n_cellx, n_celly, n_cellz, prob_lo, prob_hi, dm_in, &
         max_grid_size, pmask
    use derivative_stencil_module, only : stencil_ng

    type(layout),intent(inout) :: la
    real(dp_t), intent(inout) :: dt
    real(dp_t), pointer :: dx(:)
    type(multifab), intent(inout) :: U

    ! local
    integer :: lo(dm_in), hi(dm_in), dm, ng
    type(box)          :: bx
    type(boxarray)     :: ba

    if (n_cellx .le. 8) then
       call bl_error("n_cellx must be greater than 8")
    end if

    if (n_celly .le. 8) then
       call bl_error("n_celly must be greater than 8")
    end if

    if (n_cellz .le. 8) then
       call bl_error("n_cellz must be greater than 8")
    end if

    time = ZERO
    dt   = 1.d20

    dm = dm_in
    lo = 0
    hi(1) = n_cellx-1
    hi(2) = n_celly-1
    if (dm > 2) then
       hi(3) = n_cellz - 1
    end if

    bx = make_box(lo,hi)
    
    call boxarray_build_bx(ba,bx)
    call boxarray_maxsize(ba,max_grid_size)
    call layout_build_ba(la,ba,boxarray_bbox(ba),pmask=pmask)
    call destroy(ba)

    allocate(dx(dm))
    dx(1) = (prob_hi(1)-prob_lo(1)) / n_cellx
    dx(2) = (prob_hi(2)-prob_lo(2)) / n_celly
    if (dm > 2) then
       dx(3) = (prob_hi(3)-prob_lo(3)) / n_cellz
    end if

    ng = stencil_ng

    call multifab_build(U,la,ncons,ng)
  
    call init_data(U,dx,prob_lo,prob_hi)

  end subroutine initialize_from_scratch

end module initialize_module
