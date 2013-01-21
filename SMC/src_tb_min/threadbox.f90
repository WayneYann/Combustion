module threadbox_module

  use bl_error_module
  use omp_module
  use layout_module
  use multifab_module

  implicit none

  integer, save :: numthreads, nthreads_d(3)
  integer, allocatable, save :: tb_lo(:,:,:), tb_hi(:,:,:) 
  integer, allocatable, save :: tb_glo(:,:,:), tb_ghi(:,:,:) 

  private

  public :: numthreads, build_threadbox, destroy_threadbox,  &
       tb_get_valid_lo, tb_get_valid_hi, tb_get_grown_lo, tb_get_grown_hi, &
       tb_multifab_setval

contains

  subroutine destroy_threadbox()
    deallocate(tb_lo, tb_hi, tb_glo, tb_ghi)
  end subroutine destroy_threadbox

  subroutine build_threadbox(la, ng)
    implicit none
    type(layout), intent(in) :: la
    integer, intent(in) :: ng

    type(box) :: bx
    integer :: ilocal, iglobal, idim, ndim, my_box_size(la%lap%dim), my_box_lo(la%lap%dim)
    integer, allocatable :: xsize(:), ysize(:), zsize(:)
    integer, allocatable :: xstart(:), ystart(:), zstart(:)
    integer :: i,j,k, ithread

    ! a limitation of current threadbox
    ! make sure all boxes have the same size
    
    ilocal = 1
    iglobal = global_index(la, ilocal)
    bx = get_box(la, iglobal)
    my_box_size = box_extent(bx)

    ndim = la%lap%dim

    do idim = 2, ndim
       if ( my_box_size(1) .ne. my_box_size(idim)) then
          call bl_error('Box must be a cube')
       end if
    end do

    do ilocal = 2, nlocal(la)
       iglobal = global_index(la, ilocal)
       bx = get_box(la, iglobal)
       do idim = 1, ndim
          if (my_box_size(idim) .ne. box_extent_d(bx,idim)) then
             call bl_error('All boxes in the same MPI task must have the same size')
          end if
       end do
    end do

    numthreads = omp_get_max_threads()

    if (ndim .eq. 2) then
       call bl_error('2D not supported')
    end if

    call init_thread_topology(numthreads, nthreads_d)

    if (numthreads .ne. nthreads_d(1)*nthreads_d(2)*nthreads_d(3)) then
       call bl_error("threadbox_module: How did this happen?")
    end if

    do idim=1, ndim
       if (my_box_size(idim) < nthreads_d(idim)) then
          print *, 'idim =', idim, ' box size = ', my_box_size(idim), '  threads in this direction', nthreads_d(idim)
          call bl_error("threadbox_module: Too many threads for such small box") 
       end if
    end do


    allocate(tb_lo (3,0:numthreads-1,nlocal(la)))
    allocate(tb_hi (3,0:numthreads-1,nlocal(la)))
    allocate(tb_glo(3,0:numthreads-1,nlocal(la)))
    allocate(tb_ghi(3,0:numthreads-1,nlocal(la)))

    allocate(xsize (nthreads_d(1)))
    allocate(ysize (nthreads_d(2)))
    allocate(zsize (nthreads_d(3)))
    allocate(xstart(nthreads_d(1)))
    allocate(ystart(nthreads_d(2)))
    allocate(zstart(nthreads_d(3)))

    ! valid box

    call split_domain(my_box_size(1), nthreads_d(1), xsize, xstart)
    call split_domain(my_box_size(2), nthreads_d(2), ysize, ystart)
    call split_domain(my_box_size(3), nthreads_d(3), zsize, zstart)

    ithread = 0
    do k = 1, nthreads_d(3)
       do j = 1, nthreads_d(2)
          do i = 1, nthreads_d(1)
             tb_lo(1,ithread,:) = xstart(i)
             tb_lo(2,ithread,:) = ystart(j)
             tb_lo(3,ithread,:) = zstart(k)
             tb_hi(1,ithread,:) = xstart(i) + xsize(i) - 1
             tb_hi(2,ithread,:) = ystart(j) + ysize(j) - 1
             tb_hi(3,ithread,:) = zstart(k) + zsize(k) - 1
             ithread = ithread + 1
          end do
       end do
    end do

    do ilocal = 1, nlocal(la)
       iglobal = global_index(la, ilocal)
       bx = get_box(la, iglobal)
       my_box_lo = box_lwb(bx)

       do ithread = 0, numthreads-1
          tb_lo(:,ithread,ilocal) = tb_lo(:,ithread,ilocal) + my_box_lo
          tb_hi(:,ithread,ilocal) = tb_hi(:,ithread,ilocal) + my_box_lo
       end do
    end do

    ! grown box

    call split_domain(my_box_size(1)+2*ng, nthreads_d(1), xsize, xstart)
    call split_domain(my_box_size(2)+2*ng, nthreads_d(2), ysize, ystart)
    call split_domain(my_box_size(3)+2*ng, nthreads_d(3), zsize, zstart)

    ithread = 0
    do k = 1, nthreads_d(3)
       do j = 1, nthreads_d(2)
          do i = 1, nthreads_d(1)
             tb_glo(1,ithread,:) = xstart(i)
             tb_glo(2,ithread,:) = ystart(j)
             tb_glo(3,ithread,:) = zstart(k)
             tb_ghi(1,ithread,:) = xstart(i) + xsize(i) - 1
             tb_ghi(2,ithread,:) = ystart(j) + ysize(j) - 1
             tb_ghi(3,ithread,:) = zstart(k) + zsize(k) - 1
             ithread = ithread + 1
          end do
       end do
    end do

    do ilocal = 1, nlocal(la)
       iglobal = global_index(la, ilocal)
       bx = get_box(la, iglobal)
       my_box_lo = box_lwb(bx) - ng

       do ithread = 0, numthreads-1
          tb_glo(:,ithread,ilocal) = tb_glo(:,ithread,ilocal) + my_box_lo
          tb_ghi(:,ithread,ilocal) = tb_ghi(:,ithread,ilocal) + my_box_lo
       end do
    end do

    deallocate(xsize, ysize, zsize, xstart, ystart, zstart)

  end subroutine build_threadbox


  subroutine init_thread_topology(n, nd)
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: nd(3)

    integer, allocatable :: facs(:)
    integer :: lfacs, f, nfac, rem, i, j, nmin, nmax

    lfacs = int(log(dble(n))/log(2.d0))
    allocate(facs(lfacs))
    
    nfac = 0
    f = 2
    rem = n
    do while (rem .ne. 1)
       if (mod(rem, f) .eq. 0) then
          rem = rem/f
          nfac = nfac+1
          facs(nfac) = f
       else
          f = f + 1
       end if
    end do

    nd = 1
    do i = nfac, 1, -1
       j = minloc(nd,1)
       nd(j) = nd(j) * facs(i)
    end do

    nmin = minval(nd)
    nmax = maxval(nd)

    nd(1) = nmin
    nd(2) = n / (nmin*nmax)
    nd(3) = nmax

    deallocate(facs)

  end subroutine init_thread_topology


  subroutine split_domain(totalsize, n, isize, start)
    implicit none
    integer, intent(in) :: totalsize, n
    integer, intent(out) :: isize(n), start(n)

    integer :: i, iavg, ileft

    iavg = int(totalsize/n)
    ileft = totalsize - iavg*n

    start(1) = 0
    isize(1) = iavg
    do i=2, n
       start(i) = start(i-1) + isize(i-1)
       if (ileft > 0) then
          isize(i) = iavg + 1
          ileft = ileft - 1
       else
          isize(i) = iavg
       end if
    end do

  end subroutine split_domain


  pure function tb_get_valid_lo(ithread, ilocal) result (lo)
    implicit none
    integer, intent(in) :: ithread, ilocal
    integer, dimension(3) :: lo
    lo = tb_lo(:,ithread,ilocal)
  end function tb_get_valid_lo

  pure function tb_get_valid_hi(ithread, ilocal) result (hi)
    implicit none
    integer, intent(in) :: ithread, ilocal
    integer, dimension(3) :: hi
    hi = tb_hi(:,ithread,ilocal)
  end function tb_get_valid_hi


  pure function tb_get_grown_lo(ithread, ilocal) result (lo)
    implicit none
    integer, intent(in) :: ithread, ilocal
    integer, dimension(3) :: lo
    lo = tb_glo(:,ithread,ilocal)
  end function tb_get_grown_lo

  pure function tb_get_grown_hi(ithread, ilocal) result (hi)
    implicit none
    integer, intent(in) :: ithread, ilocal
    integer, dimension(3) :: hi
    hi = tb_ghi(:,ithread,ilocal)
  end function tb_get_grown_hi


  subroutine tb_multifab_setval(mf, val, all)
    type(multifab), intent(inout) :: mf
    real(dp_t), intent(in) :: val
    logical, intent(in), optional :: all

    logical :: lall
    integer :: ib, i, j, k, tid, wlo(3), whi(3)
    double precision, pointer :: p(:,:,:,:)

    if (present(all)) then
       lall = all
    else
       lall = .false.
    end if

    if (lall) then
       !$omp parallel private(tid, ib, i, j, k, wlo, whi, p)
       tid = omp_get_thread_num()
       do ib = 1, nfabs(mf)
          p => dataptr(mf, ib)
          wlo = tb_get_grown_lo(tid, ib)
          whi = tb_get_grown_hi(tid, ib)
          p(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),:) = val
       end do
       !$omp end parallel
    else
       !$omp parallel private(tid, ib, i, j, k, wlo, whi)
       tid = omp_get_thread_num()
       do ib = 1, nlocal(mf%la)
          p => dataptr(mf, ib)
          wlo = tb_get_valid_lo(tid, ib)
          whi = tb_get_valid_hi(tid, ib)
          p(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),:) = val
       end do
       !$omp end parallel
    end if

  end subroutine tb_multifab_setval

end module threadbox_module
