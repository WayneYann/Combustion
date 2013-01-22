module advance_module

  use bl_error_module
  use kernels_module
  use multifab_module
  use omp_module
  use threadbox_module
  use time_module
  use transport_properties
  use variables_module

  implicit none

  private
  public advance, overlapped_part

contains

  subroutine advance(U, dt, courno, dx, istep)

    use smcdata_module
    implicit none

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt, courno
    double precision,  intent(in   ) :: dx(U%dim)
    integer, intent(in) :: istep

    type(bl_prof_timer), save :: bpt_rkstep1, bpt_rkstep2, bpt_rkstep3

    call tb_multifab_setval(Unew, 0.d0, .true.)

    ! RK Step 1
    call build(bpt_rkstep1, "rkstep1")   !! vvvvvvvvvvvvvvvvvvvvvvv timer

    call dUdt(U, Uprime, Q, mu, xi, lam, Ddiag, dx, courno, istep)
    call set_dt(dt, courno, istep)
    call update_rk3(Zero,Unew, One,U, dt,Uprime)
    call reset_density(Unew)

    call destroy(bpt_rkstep1)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    ! RK Step 2
    call build(bpt_rkstep2, "rkstep2")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    call dUdt(Unew, Uprime, Q, mu, xi, lam, Ddiag, dx)
    call update_rk3(OneQuarter, Unew, ThreeQuarters, U, OneQuarter*dt, Uprime)
    call reset_density(Unew)
    call destroy(bpt_rkstep2)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    ! RK Step 3
    call build(bpt_rkstep3, "rkstep3")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    call dUdt(Unew, Uprime, Q, mu, xi, lam, Ddiag, dx)
    call update_rk3(OneThird, U, TwoThirds, Unew, TwoThirds*dt, Uprime)
    call reset_density(U)
    call destroy(bpt_rkstep3)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    if (contains_nan(U)) then
       call bl_error("U contains nan")
    end if

  end subroutine advance


  !
  ! Compute new time-step size
  !
  subroutine set_dt(dt, courno, istep)

    use probin_module, only : cflfac, fixed_dt, init_shrink, max_dt, max_dt_growth, small_dt, stop_time

    double precision, intent(inout) :: dt
    double precision, intent(in   ) :: courno
    integer,          intent(in   ) :: istep

    double precision :: dtold

    if (fixed_dt > 0.d0) then

       dt = fixed_dt

       if (parallel_IOProcessor()) then
          print*, ""
          print*, "Setting fixed dt =",dt
          print*, ""
       end if

    else

       dtold = dt
       dt    = cflfac / courno

       if (parallel_IOProcessor()) then
          print*, "CFL: dt =", dt
       end if

       if (istep .eq. 1) then
          dt = dt * init_shrink
          if (parallel_IOProcessor()) then
             print*,'Limited by init_shrink: dt =',dt
          end if
       else
          if (dt .gt. dtold * max_dt_growth) then
             dt = dtold * max_dt_growth
             if (parallel_IOProcessor()) then
                print*,'Limited by dt_growth: dt =',dt
             end if
          end if
       end if

       if(dt .gt. max_dt) then
          if (parallel_IOProcessor()) then
             print*,'Limited by max_dt: dt =',max_dt
          end if
          dt = max_dt
       end if

       if (dt < small_dt) then
          call bl_error("ERROR: timestep < small_dt")
       endif

       if (stop_time > 0.d0) then
          if (time + dt > stop_time) then
             dt = stop_time - time
             if (parallel_IOProcessor()) then
                print*, "Limited by stop_time: dt =",dt
             end if
          end if
       end if
       
       if (parallel_IOProcessor()) then
          print *, ""
       end if

    end if

  end subroutine set_dt



  !
  ! Compute U1 = a U1 + b U2 + c Uprime.
  !
  subroutine update_rk3 (a,U1,b,U2,c,Uprime)

    type(multifab),   intent(in   ) :: U2, Uprime
    type(multifab),   intent(inout) :: U1
    double precision, intent(in   ) :: a, b, c

    integer :: lo(U1%dim), hi(U1%dim), i, j, k, m, n, nc, tid
    double precision, pointer, dimension(:,:,:,:) :: u1p, u2p, upp

    nc = ncomp(U1)

    !$omp parallel private(tid,i,j,k,m,n,lo,hi,u1p,u2p,upp)
    tid = omp_get_thread_num()
    do n=1,nfabs(U1)
       u1p => dataptr(U1,    n)
       u2p => dataptr(U2,    n)
       upp => dataptr(Uprime,n)

       lo = tb_get_valid_lo(tid, n)
       hi = tb_get_valid_hi(tid, n)

       do m = 1, nc
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   u1p(i,j,k,m) = a * u1p(i,j,k,m) + b * u2p(i,j,k,m) + c * upp(i,j,k,m)
                end do
             end do
          end do
       end do
    end do
    !$omp end parallel

  end subroutine update_rk3


  !
  ! Compute dU/dt given U using the narrow stencil.
  !
  ! The Courant number (courno) is also computed if passed.
  !
  subroutine dUdt (U, Uprime, Q, mu, xi, lam, Ddiag, &
       dx, courno, istep)

    use probin_module, only : overlap_comm_comp, overlap_comm_gettrans, cfl_int, fixed_dt

    type(multifab),   intent(inout) :: U, Uprime, Q, mu, xi, lam, Ddiag
    double precision, intent(in   ) :: dx(U%dim)
    integer,          intent(in   ), optional :: istep
    double precision, intent(inout), optional :: courno

    integer :: lo(U%dim), hi(U%dim)
    integer :: n, ng, tid
    integer :: ng_ctoprim, ng_gettrans

    logical :: update_courno
    double precision :: courno_proc

    type(mf_fb_data) :: U_fb_data

    integer :: qlo(4), qhi(4), uplo(4), uphi(4), ulo(4), uhi(4)
    double precision, pointer, dimension(:,:,:,:) :: up, qp, mup, xip, lamp, Ddp, upp

    type(bl_prof_timer), save :: bpt_ctoprim, bpt_gettrans, bpt_hypdiffterm
    type(bl_prof_timer), save :: bpt_chemterm, bpt_courno

    ng = nghost(U)

    update_courno = .false.
    if (present(courno) .and. present(istep) .and. fixed_dt.le.0.d0) then
       if (mod(istep,cfl_int).eq.1 .or. cfl_int.le.1) then
          update_courno = .true.
       end if
    end if

    if (overlap_comm_comp) then
       call bl_error("overlap_comm_comp not supported yet!")
       call multifab_fill_boundary_nowait(U, U_fb_data)
       call multifab_fill_boundary_test(U, U_fb_data)
    else
       call multifab_fill_boundary(U)
       U_fb_data%rcvd = .true.
    end if

    call tb_multifab_setval(Uprime, ZERO)

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    if (U_fb_data%rcvd) then
       ng_ctoprim = ng
    else
       ng_ctoprim = 0
    end if

    !
    ! Calculate primitive variables based on U
    !
    call build(bpt_ctoprim, "ctoprim")    !! vvvvvvvvvvvvvvvvvvvvvvv timer
    call ctoprim(U, Q, ng_ctoprim)
    call destroy(bpt_ctoprim)             !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    if (update_courno) then
       courno_proc = -1.d50
       call compute_courno(Q, dx, courno_proc)
    end if

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    ! 
    ! chemistry
    !
    call build(bpt_chemterm, "chemterm")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    !$omp parallel private(tid,n,qp,upp,qlo,qhi,uplo,uphi,lo,hi)
    tid = omp_get_thread_num()
    do n=1,nfabs(Q)
       qp  => dataptr(Q,n)
       upp => dataptr(Uprime,n)

       qlo = lbound(qp)
       qhi = ubound(qp)
       uplo = lbound(upp)
       uphi = ubound(upp)

       lo = tb_get_valid_lo(tid,n)
       hi = tb_get_valid_hi(tid,n)

       call chemterm_3d(lo,hi,qp,qlo(1:3),qhi(1:3),upp,uplo(1:3),uphi(1:3))
    end do
    !$omp end parallel 
    call destroy(bpt_chemterm)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    if (overlap_comm_comp) then
       if (overlap_comm_gettrans) then
          call multifab_fill_boundary_test(U, U_fb_data)
       else
          call multifab_fill_boundary_waitrecv(U, U_fb_data)
       end if
    end if

    if (U_fb_data%rcvd) then
       ng_gettrans = ng
    else
       ng_gettrans = 0
    end if

    ! Fill ghost cells here for get_transport_properties
    if (ng_gettrans .eq. ng .and. ng_ctoprim .eq. 0) then
       call build(bpt_ctoprim, "ctoprim")    !! vvvvvvvvvvvvvvvvvvvvvvv timer
       call ctoprim(U, Q, ghostcells_only=.true.)
       call destroy(bpt_ctoprim)             !! ^^^^^^^^^^^^^^^^^^^^^^^ timer
       ng_ctoprim = ng 
    end if

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    !
    ! transport coefficients
    !
    call build(bpt_gettrans, "gettrans")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    call get_transport_properties(Q, mu, xi, lam, Ddiag, ng_gettrans)
    call destroy(bpt_gettrans)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    if (overlap_comm_comp) then
       call multifab_fill_boundary_waitrecv(U, U_fb_data)

       if (ng_ctoprim .eq. 0) then
          call build(bpt_ctoprim, "ctoprim")    !! vvvvvvvvvvvvvvvvvvvvvvv timer
          call ctoprim(U, Q, ghostcells_only=.true.)
          call destroy(bpt_ctoprim)             !! ^^^^^^^^^^^^^^^^^^^^^^^ timer
       end if

       if (ng_gettrans .eq. 0) then
          call build(bpt_gettrans, "gettrans")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
          call get_transport_properties(Q, mu, xi, lam, Ddiag, ghostcells_only=.true.)
          call destroy(bpt_gettrans)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer
       end if

       call multifab_fill_boundary_finish(U, U_fb_data)
    end if


    !
    ! Hyperbolic and Transport terms
    !
    call build(bpt_hypdiffterm, "hypdiffterm")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    !$omp parallel private(tid,n,lo,hi,up,ulo,uhi,upp,uplo,uphi,qp,qlo,qhi) &
    !$omp private(mup,xip,lamp,Ddp)
    tid = omp_get_thread_num()
    do n=1,nfabs(Q)
       up => dataptr(U,n)
       upp=> dataptr(Uprime,n)
       qp => dataptr(Q,n)
       mup  => dataptr(mu   , n)
       xip  => dataptr(xi   , n)
       lamp => dataptr(lam  , n)
       Ddp  => dataptr(Ddiag, n)

       ulo = lbound(up)
       uhi = ubound(up)
       qlo = lbound(qp)
       qhi = ubound(qp)
       uplo = lbound(upp)
       uphi = ubound(upp)

       lo = tb_get_valid_lo(tid,n)
       hi = tb_get_valid_hi(tid,n)

       call narrow_diffterm_3d(lo,hi,dx,qp,qlo(1:3),qhi(1:3),upp,uplo(1:3),uphi(1:3), &
            mup,xip,lamp,Ddp)

       call hypterm_3d(lo,hi,dx,up,ulo(1:3),uhi(1:3),qp,qlo(1:3),qhi(1:3),&
            upp,uplo(1:3),uphi(1:3))

    end do
    !$omp end parallel
    call destroy(bpt_hypdiffterm)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer


    if (update_courno) then
       call build(bpt_courno, "courno")
       call parallel_reduce(courno, courno_proc, MPI_MAX)
       call destroy(bpt_courno)
    end if

  end subroutine dUdt

  subroutine compute_courno(Q, dx, courno)
    type(multifab), intent(in) :: Q
    double precision, intent(in) :: dx(Q%dim)
    double precision, intent(inout) :: courno

    integer :: n, lo(Q%dim), hi(Q%dim), qlo(4), qhi(4), tid
    double precision :: courno_thread
    double precision, pointer :: qp(:,:,:,:)

    !$omp parallel private(tid, n, lo, hi, qlo, qhi, qp, courno_thread) &
    !$omp reduction(max:courno)
    tid = omp_get_thread_num()
    do n=1,nfabs(Q)
       qp => dataptr(Q,n)
       qlo = lbound(qp)
       qhi = ubound(qp)

       lo = tb_get_valid_lo(tid, n)
       hi = tb_get_valid_hi(tid, n)
       
       courno_thread = 0.d0

       call comp_courno_3d(lo,hi,dx,qp,qlo(1:3),qhi(1:3),courno_thread)

       courno = max(courno, courno_thread)
    end do
    !$omp end parallel
  end subroutine compute_courno


  ! only for testing communication and computation overlapping
  subroutine overlapped_part(U, U_fb_data)

    use chemistry_module, only : nspecies
    use probin_module, only : overlap_comm_gettrans

    type(multifab),   intent(inout) :: U
    type(mf_fb_data), intent(inout) :: U_fb_data

    integer :: tid, ng, ng_ctoprim, ng_gettrans, n, lo(U%dim), hi(U%dim)
    integer :: qlo(4), qhi(4), uplo(4), uphi(4)
    type(layout)     :: la
    type(multifab)   :: Q, Uprime, mu, xi, lam, Ddiag
    double precision, pointer, dimension(:,:,:,:) :: qp, upp

    call multifab_fill_boundary_test(U, U_fb_data)

    ng = nghost(U)
    la = get_layout(U)

    call multifab_build(Q, la, nprim, ng)

    call multifab_build(Uprime, la, ncons, 0)
    call multifab_setval(Uprime, 0.d0)

    call multifab_fill_boundary_test(U, U_fb_data)

    if (overlap_comm_gettrans) then
       call multifab_build(mu , la, 1, ng)
       call multifab_build(xi , la, 1, ng)
       call multifab_build(lam, la, 1, ng)
       call multifab_build(Ddiag, la, nspecies, ng)
    end if

    call multifab_fill_boundary_test(U, U_fb_data)
    
    ng_ctoprim = 0
    call ctoprim(U, Q, ng_ctoprim)

    call multifab_fill_boundary_test(U, U_fb_data)

    !$omp parallel private(tid,n,qp,upp,qlo,qhi,uplo,uphi,lo,hi)
    tid = omp_get_thread_num()
    do n=1,nfabs(Q)
       qp  => dataptr(Q,n)
       upp => dataptr(Uprime,n)

       qlo = lbound(qp)
       qhi = ubound(qp)
       uplo = lbound(upp)
       uphi = ubound(upp)

       lo = tb_get_valid_lo(tid,n)
       hi = tb_get_valid_hi(tid,n)

       call chemterm_3d(lo,hi,qp,qlo(1:3),qhi(1:3),upp,uplo(1:3),uphi(1:3))
    end do
    !$omp end parallel 
    
    call multifab_fill_boundary_test(U, U_fb_data)

    if (overlap_comm_gettrans) then
       ng_gettrans = 0
       call get_transport_properties(Q, mu, xi, lam, Ddiag, ng_gettrans)
    end if

    call destroy(Q)
    call destroy(Uprime)
    if (overlap_comm_gettrans) then
       call destroy(mu)
       call destroy(xi)
       call destroy(lam)
       call destroy(Ddiag)
    end if

  end subroutine overlapped_part

end module advance_module
