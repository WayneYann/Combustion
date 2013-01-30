module threadbox_module

  use bl_error_module
  use omp_module
  use layout_module
  use multifab_module

  implicit none

  integer, save :: numthreads, nthreads_d(3), ng, nb, ndim
  integer, save :: numboxgroups, boxgroupsize
  integer, save :: nthreadsperbox

  logical, allocatable, save :: worktodo(:,:)
  integer, allocatable, save :: tb_lo(:,:,:), tb_hi(:,:,:) 
  integer, allocatable, save :: tb_glo(:,:,:), tb_ghi(:,:,:) 

  private

  public ::build_threadbox, destroy_threadbox,  &
       tb_get_valid_lo, tb_get_valid_hi, tb_get_grown_lo, tb_get_grown_hi, &
       tb_multifab_setval, tb_worktodo

contains

  subroutine destroy_threadbox()
    if (allocated(tb_lo)) deallocate(tb_lo)
    if (allocated(tb_hi)) deallocate(tb_hi)
    if (allocated(tb_glo)) deallocate(tb_glo)
    if (allocated(tb_ghi)) deallocate(tb_ghi)
    if (allocated(worktodo)) deallocate(worktodo)
  end subroutine destroy_threadbox

  subroutine build_threadbox(la, ng_in)
    use probin_module, only : tb_split_dim, tb_collapse_boxes, tb_idim_more, tb_idim_less
    implicit none
    type(layout), intent(in) :: la
    integer, intent(in) :: ng_in

    integer :: ibox
    type(box) :: bx, bxg

    call destroy_threadbox()

    ndim = la%lap%dim
    if (ndim .eq. 2) then
       call bl_error('2D not supported')
    end if
    ng = ng_in
    nb = nlocal(la) ! number of local boxes
    numthreads = omp_get_max_threads()

    allocate(tb_lo (3,0:numthreads-1,nb))
    allocate(tb_hi (3,0:numthreads-1,nb))
    allocate(tb_glo(3,0:numthreads-1,nb))
    allocate(tb_ghi(3,0:numthreads-1,nb))

    allocate(worktodo(0:numthreads-1,nb))

    call setup_boxgroups(tb_collapse_boxes)

    call init_worktodo()

    if (tb_idim_more .eq. tb_idim_less) then
       call bl_error("threadbox_module: tb_idim_more .eq. tb_idim_less")
    end if
    if (tb_idim_more < 1 .or. tb_idim_more > 3) then
       call bl_error("threadbox_module: invalid tb_idim_more")
    end if
    if (tb_idim_less < 1 .or. tb_idim_less > 3) then
       call bl_error("threadbox_module: invalid tb_idim_less")
    end if
    call init_thread_topology(nthreadsperbox, nthreads_d, tb_split_dim, &
         tb_idim_more, tb_idim_less)

    call check_boxsize(la)

    do ibox=1, nb
       bx = get_box(la, global_index(la, ibox))
       call init_threadbox(ibox, bx, tb_lo(:,:,ibox), tb_hi(:,:,ibox))

       bxg = grow(bx, ng)
       call init_threadbox(ibox, bxg, tb_glo(:,:,ibox), tb_ghi(:,:,ibox))
    end do

  end subroutine build_threadbox


  ! Boxes are divided into groups.
  ! Then boxes in a group are divided by the number of threads.
  subroutine setup_boxgroups(collapse)
    implicit none
    logical, intent(in) :: collapse
    if (collapse) then
       boxgroupsize = greatest_common_factor(nb, numthreads)
    else
       boxgroupsize = 1
    end if
    numboxgroups = nb / boxgroupsize
    nthreadsperbox = numthreads / boxgroupsize
  contains
    recursive function greatest_common_factor(a, b) result(c)
      implicit none
      integer :: c
      integer, intent(in) :: a, b
      if (a.eq.b) then
         c = a
      else if (a.gt.b) then
         c = greatest_common_factor(a-b,b)
      else
         c = greatest_common_factor(a,b-a)
      end if
    end function greatest_common_factor
  end subroutine setup_boxgroups


  ! All threads will work on a group, but not necessarily on a box.
  subroutine init_worktodo()
    implicit none
    integer :: ithread, ibox, ib2, ib3
    do ibox=1,nb
       ib2 = mod(ibox-1,boxgroupsize) + 1 ! box index in a group
       do ithread=0,numthreads-1
          ib3 = (ithread/nthreadsperbox) + 1 ! box index for this thread
          if (ib2 .eq. ib3) then
             worktodo(ithread, ibox) = .true.
          else
             worktodo(ithread, ibox) = .false.
          end if
       end do
    end do
  end subroutine init_worktodo


  subroutine init_thread_topology(n, n3d, d, imore, iless)
    implicit none
    integer, intent(in) :: n, d, imore, iless
    integer, intent(out) :: n3d(3)

    integer, allocatable :: facs(:)
    integer :: lfacs, f, nfac, rem, i, j, nmin, nmax
    integer :: n2d(2)

    if (d.eq.1) then
       n3d = 1
       n3d(imore) = n
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
    else
       n2d = 1
       do i = nfac, 1, -1
          j = minloc(n2d,1)
          n2d(j) = n2d(j) * facs(i)
       end do
       n3d(1) = 1
       n3d(2:3) = n2d
    end if

    nmin = minval(n3d)
    nmax = maxval(n3d)

    n3d = n / (nmin*nmax)
    n3d(iless) = nmin
    n3d(imore) = nmax

    deallocate(facs)

  end subroutine init_thread_topology


  subroutine check_boxsize(la)
    implicit none
    type(layout),intent(in) :: la

    integer :: ilocal, iglobal, idim
    integer :: box_size_1(la%lap%dim), box_size_i(la%lap%dim)
    type(box) :: bx

    ilocal = 1
    iglobal = global_index(la, ilocal)
    bx = get_box(la, iglobal)
    box_size_1 = box_extent(bx)

    do ilocal = 2, nb

       iglobal = global_index(la, ilocal)
       bx = get_box(la, iglobal)
       box_size_i = box_extent(bx)

       do idim = 1, ndim
          if (box_size_1(idim) .ne. box_size_i(idim)) then
             call bl_error('All boxes in the same MPI task must have the same size')
             ! make sure all boxes have the same size
             ! otherwise, the load is unbalanced
          end if
       end do

    end do

    do ilocal = 1, nb

       iglobal = global_index(la, ilocal)
       bx = get_box(la, iglobal)
       box_size_i = box_extent(bx)

       do idim = 1, ndim
          if (box_size_i(idim) < nthreads_d(idim)) then
             print *, 'Box #', iglobal, 'idim =', idim, ' box size = ', box_size_i(idim), &
                  '  threads in this direction', nthreads_d(idim)
             call bl_error("threadbox_module: Too many threads for such small box") 
          end if
       end do

    end do

  end subroutine check_boxsize


  subroutine init_threadbox(ibox, bx, tbx_lo, tbx_hi)
    implicit none
    integer, intent(in) :: ibox
    type(box), intent(in) :: bx
    integer :: tbx_lo(3,0:numthreads-1), tbx_hi(3,0:numthreads-1)
    
    integer, allocatable :: xsize(:), ysize(:), zsize(:)
    integer, allocatable :: xstart(:), ystart(:), zstart(:)
    integer, allocatable :: zero_lo(:,:), zero_hi(:,:)
    integer :: my_box_size(3), my_box_lo(3)
    integer :: i,j,k, itbox, ithread

    allocate(zero_lo(3,0:nthreadsperbox-1))
    allocate(zero_hi(3,0:nthreadsperbox-1))

    allocate(xsize (nthreads_d(1)))
    allocate(ysize (nthreads_d(2)))
    allocate(zsize (nthreads_d(3)))
    allocate(xstart(nthreads_d(1)))
    allocate(ystart(nthreads_d(2)))
    allocate(zstart(nthreads_d(3)))

    my_box_size = box_extent(bx)
    my_box_lo = box_lwb(bx)

    ! valid box

    call split_domain(my_box_size(1), nthreads_d(1), xsize, xstart)
    call split_domain(my_box_size(2), nthreads_d(2), ysize, ystart)
    call split_domain(my_box_size(3), nthreads_d(3), zsize, zstart)

    itbox = 0
    do k = 1, nthreads_d(3)
       do j = 1, nthreads_d(2)
          do i = 1, nthreads_d(1)
             zero_lo(1,itbox) = xstart(i)
             zero_lo(2,itbox) = ystart(j)
             zero_lo(3,itbox) = zstart(k)
             zero_hi(1,itbox) = xstart(i) + xsize(i) - 1
             zero_hi(2,itbox) = ystart(j) + ysize(j) - 1
             zero_hi(3,itbox) = zstart(k) + zsize(k) - 1
             itbox = itbox + 1
          end do
       end do
    end do

    do ithread = 0, numthreads-1
       if (tb_worktodo(ithread, ibox)) then
          itbox = mod(ithread, nthreadsperbox)
          tbx_lo(:,ithread) = zero_lo(:,itbox) + my_box_lo
          tbx_hi(:,ithread) = zero_hi(:,itbox) + my_box_lo
       else
          tbx_lo(:,ithread) = HUGE(0)
          tbx_hi(:,ithread) = -(HUGE(0)-1)
       end if
    end do

    deallocate(xsize, ysize, zsize, xstart, ystart, zstart,zero_lo,zero_hi)

  contains

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

  end subroutine init_threadbox


  function tb_get_valid_lo(ithread, ilocal) result (lo)
    implicit none
    integer, intent(in) :: ithread, ilocal
    integer, dimension(3) :: lo
    lo = tb_lo(:,ithread,ilocal)
  end function tb_get_valid_lo

  function tb_get_valid_hi(ithread, ilocal) result (hi)
    implicit none
    integer, intent(in) :: ithread, ilocal
    integer, dimension(3) :: hi
    hi = tb_hi(:,ithread,ilocal)
  end function tb_get_valid_hi


  function tb_get_grown_lo(ithread, ilocal) result (lo)
    implicit none
    integer, intent(in) :: ithread, ilocal
    integer, dimension(3) :: lo
    lo = tb_glo(:,ithread,ilocal)
  end function tb_get_grown_lo

  function tb_get_grown_hi(ithread, ilocal) result (hi)
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
          if (.not.tb_worktodo(tid,ib)) cycle
          p => dataptr(mf, ib)
          wlo = tb_get_grown_lo(tid, ib)
          whi = tb_get_grown_hi(tid, ib)
          p(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),:) = val
       end do
       !$omp end parallel
    else
       !$omp parallel private(tid, ib, wlo, whi, p)
       tid = omp_get_thread_num()
       do ib = 1, nfabs(mf)
          if (.not.tb_worktodo(tid,ib)) cycle
          p => dataptr(mf, ib)
          wlo = tb_get_valid_lo(tid, ib)
          whi = tb_get_valid_hi(tid, ib)
          p(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),:) = val
       end do
       !$omp end parallel
    end if

  end subroutine tb_multifab_setval


  function tb_worktodo(ithread, ibox) result(r)
    implicit none
    logical :: r
    integer,intent(in) :: ithread, ibox
    r = worktodo(ithread,ibox)
  end function tb_worktodo


end module threadbox_module
