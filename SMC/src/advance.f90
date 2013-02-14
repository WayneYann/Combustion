module advance_module

  use bl_error_module
  use derivative_stencil_module
  use kernels_module
  use multifab_module
  use omp_module
  use nscbc_module
  use smc_bc_module
  use threadbox_module
  use time_module
  use transport_properties
  use variables_module
  use sdclib
  use sdclib_multifab
  use sdcquad_module

  use chemistry_module, only : nspecies

  implicit none

  private
  public advance, overlapped_part, srf1eval, srf1post, mrf1eval, mrf2eval

contains

  subroutine advance(U, dt, courno, dx, sdc, istep)

    use probin_module, only : advance_method

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt, courno
    double precision,  intent(in   ) :: dx(U%dim)
    integer,           intent(in   ) :: istep
    type(sdc_t),       intent(in   ) :: sdc

    select case(advance_method)
    case(1)
       call advance_rk3(U,dt,courno,dx,istep)
    case (2)
       call advance_sdc(U,dt,courno,dx,sdc,istep)
    case(3)
       call advance_multi_sdc(U,dt,courno,dx,sdc,istep)
    case default
       call bl_error("Invalid advance_method.")
    end select

    if (contains_nan(U)) then
       call bl_error("U contains nan")
    end if

  end subroutine advance

  !
  ! Advance U using SSP RK3
  !
  subroutine advance_rk3 (U,dt,courno,dx,istep)

    use smcdata_module, only : Unew, Uprime
    implicit none

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt, courno
    double precision,  intent(in   ) :: dx(U%dim)
    integer, intent(in) :: istep

    type(bl_prof_timer), save :: bpt_rkstep1, bpt_rkstep2, bpt_rkstep3

    call tb_multifab_setval(Unew, 0.d0, .true.)

    ! RK Step 1
    call build(bpt_rkstep1, "rkstep1")   !! vvvvvvvvvvvvvvvvvvvvvvv timer

    call dUdt(U, Uprime, dx, courno=courno, istep=istep)
    call set_dt(dt, courno, istep)
    call update_rk3(Zero,Unew, One,U, dt,Uprime)
    call reset_density(Unew)
    call impose_hard_bc(Unew)

    call destroy(bpt_rkstep1)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    ! RK Step 2
    call build(bpt_rkstep2, "rkstep2")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    call dUdt(Unew, Uprime, dx)
    call update_rk3(OneQuarter, Unew, ThreeQuarters, U, OneQuarter*dt, Uprime)
    call reset_density(Unew)
    call impose_hard_bc(Unew)
    call destroy(bpt_rkstep2)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    ! RK Step 3
    call build(bpt_rkstep3, "rkstep3")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    call dUdt(Unew, Uprime, dx)
    call update_rk3(OneThird, U, TwoThirds, Unew, TwoThirds*dt, Uprime)
    call reset_density(U)
    call impose_hard_bc(U)
    call destroy(bpt_rkstep3)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

  end subroutine advance_rk3


  !
  ! Advance U using single-rate SDC time-stepping
  !
  subroutine advance_sdc(U, dt, courno, dx, sdc, istep)

    use smcdata_module, only : Q
    use probin_module, only : cfl_int, fixed_dt

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt, courno
    double precision,  intent(in   ) :: dx(U%dim)
    type(sdc_t),       intent(in   ) :: sdc
    integer,           intent(in   ) :: istep

    logical :: update_courno
    double precision :: courno_proc

    integer :: k
    double precision :: res
    type(layout) :: la
    type(multifab), target :: R

    type(bl_prof_timer), save :: bpt_sdc_prep, bpt_sdc_iter

    !
    ! set dt
    !

    ! this really belongs in the first feval of an sdc sweep

    update_courno = .false.
    if (fixed_dt.le.0.d0) then
       if (mod(istep,cfl_int).eq.1 .or. cfl_int.le.1) then
          update_courno = .true.
       end if
    end if

    if (update_courno) then
       courno_proc = -1.d50
       call ctoprim(U, Q, 0)
       call compute_courno(Q, dx, courno_proc)
       call parallel_reduce(courno, courno_proc, MPI_MAX)
    end if

    call set_dt(dt, courno, istep)

    !
    ! advance (pass control to sdclib)
    !
    call build(bpt_sdc_prep, "sdc_prep")
    call sdc_srset_set_q0(sdc%srset, mfptr(U))
    call sdc_srset_spread(sdc%srset, 0.0d0)
    call destroy(bpt_sdc_prep)

    if (sdc%tol_residual > 0.d0) then
       la = get_layout(U)
    end if

    call build(bpt_sdc_iter, "sdc_iter")
    do k = 1, sdc%iters
       call sdc_srset_sweep(sdc%srset, 0.0d0, dt)

       ! check residual
       if (sdc%tol_residual > 0.d0) then
          call build(R, la, ncons, 0)
          call sdc_srset_residual(sdc%srset, dt, mfptr(R))
          call parallel_reduce(res, norm_inf(R), MPI_MAX)
          call destroy(R)

          if (parallel_IOProcessor()) then
             print *, "SDC: iter:", k, "residual:", res
          end if

          if (res < sdc%tol_residual) &
               exit
       end if
    end do
    call destroy(bpt_sdc_iter)

    call sdc_srset_get_qend(sdc%srset, mfptr(U))    

  end subroutine advance_sdc


  !
  ! SDCLib callbacks
  !
  subroutine srf1eval(Fptr, Uptr, t, ctxptr) bind(c)
    type(c_ptr),    intent(in), value :: Fptr, Uptr, ctxptr
    real(c_double), intent(in), value :: t

    type(multifab), pointer :: U, Uprime
    type(ctx_t), pointer    :: ctx

    call c_f_pointer(Uptr, U)
    call c_f_pointer(Fptr, Uprime)
    call c_f_pointer(ctxptr, ctx)

    call dUdt(U, Uprime, ctx%dx)
  end subroutine srf1eval

  subroutine srf1post(Uptr, Fptr, stateptr, ctxptr) bind(c)
    type(c_ptr), intent(in), value :: Uptr, Fptr, stateptr, ctxptr

    type(multifab), pointer :: U
    call c_f_pointer(Uptr, U)

    call reset_density(U)
    call impose_hard_bc(U)
  end subroutine srf1post

  subroutine mrf1eval(Fptr, Uptr, t, ctxptr) bind(c)
    type(c_ptr),    intent(in), value :: Fptr, Uptr, ctxptr
    real(c_double), intent(in), value :: t

    type(multifab), pointer :: U, Uprime
    type(ctx_t), pointer    :: ctx

    call c_f_pointer(Uptr, U)
    call c_f_pointer(Fptr, Uprime)
    call c_f_pointer(ctxptr, ctx)

    call dUdt_AD(U, Uprime, ctx%dx)
  end subroutine mrf1eval

  subroutine mrf2eval(Fptr, Uptr, t, ctxptr) bind(c)
    type(c_ptr),    intent(in), value :: Fptr, Uptr, ctxptr
    real(c_double), intent(in), value :: t

    type(multifab), pointer :: U, Uprime
    type(ctx_t), pointer    :: ctx

    call c_f_pointer(Uptr, U)
    call c_f_pointer(Fptr, Uprime)
    call c_f_pointer(ctxptr, ctx)

    call dUdt_R(U, Uprime, ctx%dx)
  end subroutine mrf2eval


  !
  ! Advance U using multi-rate SDC time-stepping
  !
  subroutine advance_multi_sdc(U, dt, courno, dx, sdc, istep)

    use smcdata_module, only : Q
    use probin_module, only : cfl_int, fixed_dt

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt, courno
    double precision,  intent(in   ) :: dx(U%dim)
    type(sdc_t),       intent(in   ) :: sdc
    integer,           intent(in   ) :: istep

    logical :: update_courno
    double precision :: courno_proc

    ! ideally we would have a preallocated workspace for the residual,
    ! and the computation of dt would be done in the first feval...

    !
    ! set dt
    !

    ! this really belongs in the first feval of an sdc sweep

    update_courno = .false.
    if (fixed_dt.le.0.d0) then
       if (mod(istep,cfl_int).eq.1 .or. cfl_int.le.1) then
          update_courno = .true.
       end if
    end if

    if (update_courno) then
       courno_proc = -1.d50
       call ctoprim(U, Q, 0)
       call compute_courno(Q, dx, courno_proc)
       call parallel_reduce(courno, courno_proc, MPI_MAX)
    end if

    call set_dt(dt, courno, istep)

    !
    ! advance
    !
    call sdc_mrset_set_q0(sdc%mrset, mfptr(U))
    call sdc_mrset_advance(sdc%mrset, 5, 0.0d0, dt)
    call sdc_mrset_get_qend(sdc%mrset, mfptr(U))

  end subroutine advance_multi_sdc


  !
  ! Compute new time-step size
  !
  subroutine set_dt(dt, courno, istep)

    use probin_module, only : fixed_dt, cflfac, init_shrink, max_dt_growth, &
         max_dt, small_dt, stop_time

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

       ! if (istep .eq. 1) then
       !    dt = dt * init_shrink
       !    if (parallel_IOProcessor()) then
       !       print*,'Limited by init_shrink: dt =',dt
       !    end if
       ! end if

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

       if (.not.tb_worktodo(tid,n)) cycle

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
  ! Compute advection/diffusion part of dU/dt given U.
  !
  subroutine dUdt_AD (U, Uprime, dx)
    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)

    call dUdt(U, Uprime, dx, include_r=.false.)
  end subroutine dUdt_AD


  !
  ! Compute reaction part of dU/dt given U.
  !
  subroutine dUdt_R (U, Uprime, dx)
    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)

    call dUdt(U, Uprime, dx, include_ad=.false.)
  end subroutine dUdt_R


  !
  ! Compute dU/dt given U using the narrow stencil.
  !
  ! The Courant number (courno) is also computed if passed.
  !
  subroutine dUdt (U, Uprime, dx, courno, istep, include_ad, include_r)

    use smcdata_module, only : Q, mu, xi, lam, Ddiag, Fdif
    use probin_module, only : overlap_comm_comp, overlap_comm_gettrans, cfl_int, fixed_dt

    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)
    double precision, intent(inout), optional :: courno
    integer,          intent(in   ), optional :: istep
    logical,          intent(in   ), optional :: include_ad, include_r

    integer ::    lo(U%dim),    hi(U%dim)
    integer ::   dlo(U%dim),   dhi(U%dim)
    integer ::   blo(U%dim),   bhi(U%dim)
    integer :: n, ng, tid
    integer :: ng_ctoprim, ng_gettrans

    logical :: update_courno
    double precision :: courno_proc

    type(mf_fb_data) :: U_fb_data

    logical :: inc_ad, inc_r

    integer :: qlo(4), qhi(4), uplo(4), uphi(4), ulo(4), uhi(4), flo(4), fhi(4)
    double precision, pointer, dimension(:,:,:,:) :: up, qp, mup, xip, lamp, Ddp, upp, fp

    type(bl_prof_timer), save :: bpt_ctoprim, bpt_gettrans, bpt_hypdiffterm
    type(bl_prof_timer), save :: bpt_chemterm, bpt_courno, bpt_nscbc

    ng = nghost(U)

    inc_ad = .true.; if (present(include_ad)) inc_ad = include_ad
    inc_r  = .true.; if (present(include_r))  inc_r  = include_r

    update_courno = .false.
    if (present(courno) .and. present(istep) .and. fixed_dt.le.0.d0) then
       if (mod(istep,cfl_int).eq.1 .or. cfl_int.le.1) then
          update_courno = .true.
       end if
    end if

    if (inc_ad) then
       call multifab_fill_boundary_nowait(U, U_fb_data)

       if (overlap_comm_comp) then
          call multifab_fill_boundary_test(U, U_fb_data)
       else
          call multifab_fill_boundary_finish(U, U_fb_data)
       end if
    end if

    call tb_multifab_setval(Uprime, ZERO)

    if (inc_ad .and. overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    if (.not. inc_ad) then
       ng_ctoprim = 0
    else if (U_fb_data%rcvd) then
       ng_ctoprim = ng
    else
       ng_ctoprim = 0
    end if

    !
    ! Calculate primitive variables based on U
    !
    call build(bpt_ctoprim, "ctoprim")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    call ctoprim(U, Q, ng_ctoprim)
    call destroy(bpt_ctoprim)            !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    if (inc_ad .and. overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    if (update_courno) then
       courno_proc = -1.d50
       call compute_courno(Q, dx, courno_proc)
    end if

    if (inc_ad .and. overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    !
    ! R
    !
    if (inc_r) then
       !
       ! chemistry
       !
       call build(bpt_chemterm, "chemterm")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
       !$omp parallel private(tid,n,qp,upp,qlo,qhi,uplo,uphi,lo,hi)
       tid = omp_get_thread_num()
       do n=1,nfabs(Q)

          if (.not.tb_worktodo(tid,n)) cycle

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
    end if

    !
    ! AD
    !
    if (inc_ad) then       

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
       call destroy(bpt_gettrans)               !! ^^^^^^^^^^^^^^^^^^^^^^^ timer       

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
       !$omp private(fp,flo,fhi,mup,xip,lamp,Ddp,dlo,dhi,blo,bhi)
       tid = omp_get_thread_num()
       do n=1,nfabs(Q)
          
          if (.not.tb_worktodo(tid,n)) cycle
          
          up => dataptr(U,n)
          upp=> dataptr(Uprime,n)
          fp => dataptr(Fdif,n)
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
          flo = lbound(fp)
          fhi = ubound(fp)
          
          lo = tb_get_valid_lo(tid,n)
          hi = tb_get_valid_hi(tid,n)

          call get_data_lo_hi(n,dlo,dhi)
          call get_boxbc(n,blo,bhi)

          call hypterm_3d(lo,hi,dx,up,ulo(1:3),uhi(1:3),qp,qlo(1:3),qhi(1:3),&
               upp,uplo(1:3),uphi(1:3),dlo,dhi,blo,bhi)
          
          call narrow_diffterm_3d(lo,hi,dx,qp,qlo(1:3),qhi(1:3),upp,uplo(1:3),uphi(1:3), &
               fp,flo(1:3),fhi(1:3),mup,xip,lamp,Ddp,dlo,dhi,blo,bhi)
          
       end do
       !$omp end parallel
       call destroy(bpt_hypdiffterm)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer


       !
       ! NSCBC boundary
       !
       call build(bpt_nscbc, "nscbc")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
       call nscbc(Q, U, Fdif, Uprime, dx)
       call destroy(bpt_nscbc)          !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    end if


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

       if (.not.tb_worktodo(tid,n)) cycle

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
    call tb_multifab_setval(Uprime, 0.d0)

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

       if (.not.tb_worktodo(tid,n)) cycle

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
