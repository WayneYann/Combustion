module threadbox_module

  use bl_error_module
  use omp_module
  use layout_module
  use multifab_module

  implicit none

  integer, save :: numthreads, nthreads_d(3), ng
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

  subroutine build_threadbox(la, ng_in)
    use probin_module, only : tb_split_dim
    implicit none
    type(layout), intent(in) :: la
    integer, intent(in) :: ng_in

    type(box) :: bx
    integer :: ilocal, iglobal, idim, ndim, nb
    integer :: my_box_size(la%lap%dim), my_box_lo(la%lap%dim)
    integer, allocatable :: xsize(:), ysize(:), zsize(:)
    integer, allocatable :: xstart(:), ystart(:), zstart(:)
    integer, allocatable :: zero_lo(:,:), zero_hi(:,:)
    integer :: i,j,k, ithread, it2, itstride, itstart

    ng = ng_in
    nb = nlocal(la) ! number of local boxes

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

    do ilocal = 2, nb
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

    call init_thread_topology(numthreads, nthreads_d, tb_split_dim)

    if (numthreads .ne. nthreads_d(1)*nthreads_d(2)*nthreads_d(3)) then
       call bl_error("threadbox_module: How did this happen?")
    end if

    do idim=1, ndim
       if (my_box_size(idim) < nthreads_d(idim)) then
          print *, 'idim =', idim, ' box size = ', my_box_size(idim), '  threads in this direction', nthreads_d(idim)
          call bl_error("threadbox_module: Too many threads for such small box") 
       end if
    end do


    allocate(tb_lo (3,0:numthreads-1,nb))
    allocate(tb_hi (3,0:numthreads-1,nb))
    allocate(tb_glo(3,0:numthreads-1,nb))
    allocate(tb_ghi(3,0:numthreads-1,nb))

    allocate(zero_lo(3,0:numthreads-1))
    allocate(zero_hi(3,0:numthreads-1))

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
             zero_lo(1,ithread) = xstart(i)
             zero_lo(2,ithread) = ystart(j)
             zero_lo(3,ithread) = zstart(k)
             zero_hi(1,ithread) = xstart(i) + xsize(i) - 1
             zero_hi(2,ithread) = ystart(j) + ysize(j) - 1
             zero_hi(3,ithread) = zstart(k) + zsize(k) - 1
             ithread = ithread + 1
          end do
       end do
    end do

    itstride = numthreads / nb

    do ilocal = 1, nb
       iglobal = global_index(la, ilocal)
       bx = get_box(la, iglobal)
       my_box_lo = box_lwb(bx)

       itstart = (ilocal-1)*itstride
       do ithread = 0, numthreads-1
          it2 = mod(ithread+itstart, numthreads)
          tb_lo(:,ithread,ilocal) = zero_lo(:,it2) + my_box_lo
          tb_hi(:,ithread,ilocal) = zero_hi(:,it2) + my_box_lo
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
             zero_lo(1,ithread) = xstart(i)
             zero_lo(2,ithread) = ystart(j)
             zero_lo(3,ithread) = zstart(k)
             zero_hi(1,ithread) = xstart(i) + xsize(i) - 1
             zero_hi(2,ithread) = ystart(j) + ysize(j) - 1
             zero_hi(3,ithread) = zstart(k) + zsize(k) - 1
             ithread = ithread + 1
          end do
       end do
    end do

    do ilocal = 1, nb
       iglobal = global_index(la, ilocal)
       bx = get_box(la, iglobal)
       my_box_lo = box_lwb(bx) - ng

       itstart = (ilocal-1)*itstride
       do ithread = 0, numthreads-1
          it2 = mod(ithread+itstart, numthreads)
          tb_glo(:,ithread,ilocal) = zero_lo(:,it2) + my_box_lo
          tb_ghi(:,ithread,ilocal) = zero_hi(:,it2) + my_box_lo
       end do
    end do

    deallocate(xsize, ysize, zsize, xstart, ystart, zstart,zero_lo,zero_hi)

  end subroutine build_threadbox


  subroutine init_thread_topology(n, n3d, d)
    implicit none
    integer, intent(in) :: n, d
    integer, intent(out) :: n3d(3)

    integer, allocatable :: facs(:)
    integer :: lfacs, f, nfac, rem, i, j, nmin, nmax
    integer :: n2d(2)

    if (d.eq.1) then
       n3d(1:2) = 1
       n3d(3) = n
       return
    end if

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

    if (d .eq. 3) then
       n3d = 1
       do i = nfac, 1, -1
          j = minloc(n3d,1)
          n3d(j) = n3d(j) * facs(i)
       end do
       
       nmin = minval(n3d)
       nmax = maxval(n3d)
       
       n3d(1) = nmin
       n3d(2) = n / (nmin*nmax)
       n3d(3) = nmax
    else
       n2d = 1
       do i = nfac, 1, -1
          j = minloc(n2d,1)
          n2d(j) = n2d(j) * facs(i)
       end do
       
       n3d(1) = 1
       n3d(2:3) = n2d
    end if

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
    integer :: ib, tid, wlo(3), whi(3), ngmf
    double precision, pointer :: p(:,:,:,:)

    if (present(all)) then
       lall = all
    else
       lall = .false.
    end if

    ngmf = nghost(mf)

    if (lall) then
       if (ngmf .ne. ng .and. ngmf .ne. 0) then
          call bl_error("threadbox: do not know how to do setval in this case.")
       end if
    end if

    if (lall .and. ngmf.eq.ng) then
       !$omp parallel private(tid, ib, wlo, whi, p)
       tid = omp_get_thread_num()
       do ib = 1, nfabs(mf)
          p => dataptr(mf, ib)
          wlo = tb_get_grown_lo(tid, ib)
          whi = tb_get_grown_hi(tid, ib)
          p(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),:) = val
       end do
       !$omp end parallel
    else
       !$omp parallel private(tid, ib, wlo, whi, p)
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
