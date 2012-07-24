module initialize_module

  use bl_constants_module
  use multifab_module
  use variables

  implicit none

  private

  public :: initialize_from_restart, initialize_from_scratch

contains

  subroutine initialize_from_restart

  end subroutine initialize_from_restart


  subroutine initialize_from_scratch(la,dt,dx,U)

    use init_data_module, only : init_data
    use time_module, only : time

    use probin_module, only: n_cellx, n_celly, n_cellz, prob_lo, prob_hi, dm_in, &
         max_grid_size, pmask

    type(layout),intent(inout) :: la
    real(dp_t), intent(inout) :: dt
    real(dp_t), pointer :: dx(:)
    type(multifab), intent(inout) :: U

    ! local
    integer :: lo(dm_in), hi(dm_in), dm
    type(box)          :: bx
    type(boxarray)     :: ba
    integer, parameter :: NG=4

    time = ZERO
    dt = 1.d20

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

    !
    ! Compute initial condition.
    !
    call multifab_build(U,la,ncons,NG)
  
    call init_data(U,dx,prob_lo,prob_hi)

  end subroutine initialize_from_scratch

end module initialize_module
