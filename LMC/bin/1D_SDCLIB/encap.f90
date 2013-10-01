module encap
  use sdclib
  implicit none
  type :: lmc_encap
      real(8), pointer :: vel(:,:), scal(:,:,:), divu(:,:), press(:,:)
  end type lmc_encap
  type :: lmc_encap_ctx
     integer :: nfine, nscal
  end type lmc_encap_ctx
contains

  subroutine lmc_encap_build(enc, nfine, nscal)
    type(sdc_encap), intent(inout) :: enc
    integer,         intent(in)    :: nfine, nscal

    type(lmc_encap_ctx), pointer :: ctx

    enc%id = SDC_ID_ENCAP

    allocate(ctx)
    ctx%nfine = nfine
    ctx%nscal = nscal
    enc%ctx = c_loc(ctx)

    enc%create = c_funloc(lmc_encap_create)
    enc%destroy = c_funloc(lmc_encap_destroy)
    enc%getinfo = c_funloc(lmc_encap_info)
    enc%setval = c_funloc(lmc_encap_setval)
    enc%copy = c_funloc(lmc_encap_copy)
    enc%saxpy = c_funloc(lmc_encap_saxpy)
  end subroutine lmc_encap_build

  subroutine lmc_encap_create_simple(sol, nfine, nscal)
    type(lmc_encap), intent(inout) :: sol
    integer,         intent(in)    :: nfine, nscal
    allocate(sol%vel(0:0,-2:nfine+1))
    allocate(sol%scal(0:0,-2:nfine+1,nscal))
    allocate(sol%divu(0:0,-1:nfine))
    allocate(sol%press(0:0,-1:nfine+1))
  end subroutine lmc_encap_create_simple

  type(c_ptr) function lmc_encap_create(type, ctxptr) &
       result(r) bind(c, name='lmc_encap_create')
    integer(c_int), intent(in), value :: type
    type(c_ptr),    intent(in), value :: ctxptr
    type(lmc_encap_ctx), pointer :: ctx
    type(lmc_encap),     pointer :: sol
    call c_f_pointer(ctxptr, ctx)
    allocate(sol)
    call lmc_encap_create_simple(sol, ctx%nfine, ctx%nscal)
    r = c_loc(sol)
  end function lmc_encap_create

  ! Deallocate/destroy solution.
  subroutine lmc_encap_destroy(solptr) bind(c, name='lmc_encap_destroy')
    type(c_ptr), intent(in), value :: solptr
    type(lmc_encap), pointer :: sol
    call c_f_pointer(solptr, sol)
    deallocate(sol%vel,sol%scal,sol%divu,sol%press)
    deallocate(sol)
  end subroutine lmc_encap_destroy

  ! Get info.
  subroutine lmc_encap_info(solptr, dofs, buf) bind(c, name='lmc_encap_info')
    type(c_ptr),    intent(in), value :: solptr
    integer(c_int), intent(out)       :: dofs
    type(c_ptr),    intent(out)       :: buf
    ! type(lmc_encap), pointer :: sol
    ! call c_f_pointer(solptr, sol)
    ! dofs = size(sol%array)
    ! buf  = sol%aptr
  end subroutine lmc_encap_info

  ! Set solution value.
  subroutine lmc_encap_setval(solptr, val) bind(c, name='lmc_encap_setval')
    type(c_ptr), intent(in), value :: solptr
    real(c_double),    intent(in), value :: val
    type(lmc_encap), pointer :: sol
    call c_f_pointer(solptr, sol)
    sol%vel = val
    sol%scal = val
    sol%divu = val
    sol%press = val
  end subroutine lmc_encap_setval

  ! Copy solution value.
  subroutine lmc_encap_copy(dstptr, srcptr) bind(c, name='lmc_encap_copy')
    type(c_ptr), intent(in), value :: dstptr, srcptr
    type(lmc_encap), pointer :: dst, src
    call c_f_pointer(dstptr, dst)
    call c_f_pointer(srcptr, src)
    dst%vel = src%vel
    dst%scal = src%scal
    dst%divu = src%divu
    dst%press = src%press
  end subroutine lmc_encap_copy

  ! Compute y = a x + y where a is a scalar and x and y are solutions.
  subroutine lmc_encap_saxpy(yptr, a, xptr) bind(c, name='lmc_encap_saxpy')
    type(c_ptr), intent(in), value :: yptr, xptr
    real(c_double),    intent(in), value :: a
    type(lmc_encap), pointer :: y, x
    call c_f_pointer(yptr, y)
    call c_f_pointer(xptr, x)
    y%vel = a * x%vel + y%vel
    y%scal = a * x%scal + y%scal
    y%divu = a * x%divu + y%divu
    y%press = a * x%press + y%press
  end subroutine lmc_encap_saxpy

end module encap
