module encap
  use iso_c_binding
  use multifab_module
  use advance_module, only: feval_ctx_t
  implicit none

  ! Data precision
  integer, parameter :: pfdp = kind(0.0d0)

  ! Data encapsulation: multifab
  type :: pf_encap_t

     type(multifab) :: q

  end type pf_encap_t

  interface create
     module procedure encap_create
  end interface create

  interface destroy
     module procedure encap_destroy
  end interface destroy

  interface setval
     module procedure encap_setval
  end interface setval

  interface copy
     module procedure encap_copy
  end interface copy

  interface axpy
     module procedure encap_axpy
  end interface axpy
 
  interface pack
     module procedure encap_pack
  end interface pack

  interface unpack
     module procedure encap_unpack
  end interface unpack
 
contains

  ! Allocate/create solution (spatial data set) for the given level.
  !
  ! This is called for each SDC node.
  subroutine encap_create(sol, nvars, feval, level, ctxp)
    type(pf_encap_t),  intent(inout) :: sol
    integer,           intent(in)    :: nvars, level
    logical,           intent(in)    :: feval
    type(c_ptr),       intent(in)    :: ctxp

    type(feval_ctx_t), pointer :: ctx

    call c_f_pointer(ctxp, ctx)

    if (feval) then
       call build(sol%q, ctx%la, ctx%nc, 0)
    else
       call build(sol%q, ctx%la, ctx%nc, ctx%ng)
    end if
  end subroutine encap_create

  ! Deallocate/destroy solution.
  subroutine encap_destroy(sol)
    type(pf_encap_t), intent(inout) :: sol

    call destroy(sol%q)
  end subroutine encap_destroy

  ! Set solution value.
  subroutine encap_setval(sol, val)
    type(pf_encap_t), intent(inout) :: sol
    double precision, intent(in)    :: val

    call setval(sol%q, val)
  end subroutine encap_setval

  ! Copy solution value.
  subroutine encap_copy(dst, src)
    type(pf_encap_t), intent(inout) :: dst
    type(pf_encap_t), intent(in)    :: src

    integer :: lo(src%q%dim), hi(src%q%dim), i, j, k, m, n, nc
    double precision, pointer, dimension(:,:,:,:) :: sp, dp

    nc = ncomp(src%q)

    do n=1,nboxes(src%q)
       if ( remote(src%q,n) ) cycle

       sp => dataptr(src%q,n)
       dp => dataptr(dst%q,n)

       lo = lwb(get_box(src%q,n))
       hi = upb(get_box(src%q,n))

       do m = 1, nc
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   dp(i,j,k,m) = sp(i,j,k,m)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do


  end subroutine encap_copy

  ! Pack solution q into a flat array.
  subroutine encap_pack(arr, sol)
    type(pf_encap_t), intent(in)  :: sol
    double precision, intent(out) :: arr(:)

    integer :: n, ic, lo(sol%q%dim), hi(sol%q%dim)
    double precision, pointer, dimension(:,:,:,:) :: dp

    ic = 1

    do n=1,nboxes(sol%q)
       if ( remote(sol%q,n) ) cycle

       dp => dataptr(sol%q,n)

       lo = lwb(get_box(sol%q,n))
       hi = upb(get_box(sol%q,n))

       call reshape_d_4_1(arr,ic,dp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:))

       ic = ic + volume(get_box(sol%q, n))
    end do

    ! dp => dataptr(sol%q,1)
    ! lo = lwb(get_box(sol%q,1))
    ! hi = upb(get_box(sol%q,1))

    ! arr = reshape(dp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:), [ size(arr) ])

  end subroutine encap_pack

  ! Unpack solution from a flat array.
  subroutine encap_unpack(sol, arr)
    type(pf_encap_t), intent(inout) :: sol
    double precision, intent(in)    :: arr(:)

    integer :: n, ic, lo(sol%q%dim), hi(sol%q%dim), sh(4)
    double precision, pointer, dimension(:,:,:,:) :: dp


    ! dp => dataptr(sol%q,1)
    ! lo = lwb(get_box(sol%q,1))
    ! hi = upb(get_box(sol%q,1))

    ! dp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = reshape(arr, &
    !      [ hi(1)-lo(1)+1, hi(2)-lo(2)+1, hi(3)-lo(3)+1, size(dp, 4) ] )

    ic = 1

    do n=1,nboxes(sol%q)
       if ( remote(sol%q,n) ) cycle

       dp => dataptr(sol%q,n)

       lo = lwb(get_box(sol%q,n))
       hi = upb(get_box(sol%q,n))

       call reshape_d_1_4(dp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:),arr,ic,sh)

       ic = ic + volume(get_box(sol%q, n))
    end do

  end subroutine encap_unpack

  ! Compute y = a x + y where a is a scalar and x and y are solutions.
  subroutine encap_axpy(y, a, x)
    double precision, intent(in)    :: a
    type(pf_encap_t), intent(in)    :: x
    type(pf_encap_t), intent(inout) :: y

    call saxpy(y%q, a, x%q)
  end subroutine encap_axpy

end module encap
