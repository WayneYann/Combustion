module advance_module

  use bl_constants_module
  use bl_error_module
  use multifab_module
  use probin_module
  use variables
  use time_module

  implicit none

  private

  public advance

contains

  subroutine advance(U, dt, dx)

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt
    double precision,  intent(in   ) :: dx(U%dim) 

    if (advance_method == 2) then
       call bl_error("call advance_sdc")
    else
       call advance_rk3(U, dt, dx)
    end if

  end subroutine advance


  subroutine advance_rk3 (U,dt,dx)

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt
    double precision,  intent(in   ) :: dx(U%dim)

    integer          :: nc, ng
    double precision :: courno, courno_proc
    type(layout)     :: la
    type(multifab)   :: Uprime, Unew

    nc = ncomp(U)
    ng = nghost(U)
    la = get_layout(U)

    call multifab_build(Uprime, la, nc, 0)
    call multifab_build(Unew,   la, nc, ng)

    ! RK Step 1
    courno_proc = 1.0d-50

    call dUdt(U,Uprime,courno=courno_proc)

    call parallel_reduce(courno, courno_proc, MPI_MAX)

  end subroutine advance_rk3


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute dU/dt given U.
  !
  ! The Courant number (courno) is also computed if passed.
  !
  subroutine dUdt (U,Uprime,courno)

    type(multifab),    intent(inout) :: U, Uprime
    double precision,  intent(inout), optional :: courno

    type(multifab) :: mu, xi ! viscosity
    type(multifab) :: lam ! thermal conductivity

    integer          :: lo(U%dim), hi(U%dim), i, j, k, m, n, nc, ng
    type(layout)     :: la
    type(multifab)   :: D, F, Q

    double precision, pointer, dimension(:,:,:,:) :: up, dp, fp, qp, upp, mup, xip, lamp

    nc = ncomp(U)
    ng = nghost(U)
    la = get_layout(U)

    !
    ! Sync U prior to calculating D & F
    !
    call multifab_fill_boundary(U)
    call multifab_build(D, la, nc,   0)
    call multifab_build(F, la, nc,   0)
    call multifab_build(Q, la, nc+1, ng)

    call multifab_build(mu , la, 1, ng)
    call multifab_build(xi , la, 1, ng)
    call multifab_build(lam, la, 1, ng)
    
!    call multifab_setval(mu ,  ctx%eta, all=.true.)
    call multifab_setval(xi ,     0.d0, all=.true.)
!    call multifab_setval(lam, ctx%alam, all=.true.)

    !
    ! Calculate primitive variables based on U
    !
    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

       up => dataptr(U,n)
       qp => dataptr(Q,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

!       call ctoprim(lo,hi,up,qp,ctx%dx,ng,courno=courno)
    end do

    !
    ! Calculate D
    !
    do n=1,nboxes(D)
       if ( remote(D,n) ) cycle

       qp => dataptr(Q,n)
       dp => dataptr(D,n)

       lo = lwb(get_box(D,n))
       hi = upb(get_box(D,n))

       mup  => dataptr(mu , n)
       xip  => dataptr(xi , n)
       lamp => dataptr(lam, n)
       !          call compact_diffterm(lo,hi,ng,ctx%dx,qp,dp,mup,xip,lamp)
    end do

    !
    ! Calculate F
    !
    do n=1,nboxes(F)
       if ( remote(F,n) ) cycle

       up => dataptr(U,n)
       qp => dataptr(Q,n)
       fp => dataptr(F,n)

       lo = lwb(get_box(F,n))
       hi = upb(get_box(F,n))

!       call hypterm(lo,hi,ng,ctx%dx,up,qp,fp)
    end do

    !
    ! Calculate U'
    !
    do n=1,nboxes(U)
       if ( remote(U,n) ) cycle

       dp  => dataptr(D,     n)
       fp  => dataptr(F,     n)
       upp => dataptr(Uprime,n)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       do m = 1, nc
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   upp(i,j,k,m) = dp(i,j,k,m) + fp(i,j,k,m)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do

    call destroy(D)
    call destroy(F)
    call destroy(Q)

    call destroy(mu)
    call destroy(xi)
    call destroy(lam)

  end subroutine dUdt


end module advance_module
