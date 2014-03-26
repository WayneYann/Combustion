module advance_module

  use bl_error_module
  use derivative_stencil_module
  use kernels_module
  use kernels_2d_module
  use kernels_s3d_module
  use multifab_module
  use nscbc_module
  use smc_bc_module
  use threadbox_module
  use transport_properties
  use variables_module
  use sdclib
  use sdclib_multifab
  use sdcquad_module

  use chemistry_module, only : nspecies

  implicit none

  private
  public advance, overlapped_part
  public single_sdc_feval, sdc_post_step_cb
  public multi_sdc_feval_slow, multi_sdc_feval_fast
  public sdc_get_q0

  logical, save, private :: trans_called, Mach_computed
  integer, save, private :: istep_first = -1
  integer, save, private :: istep_this
  double precision, save, private :: dt

  integer, public :: count_ad = 0, count_r = 0

contains

  subroutine advance(U, dtio, courno, dx, sdc, istep)
    use probin_module

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dtio, courno
    double precision,  intent(in   ) :: dx(3)
    integer,           intent(in   ) :: istep
    type(sdc_ctx),     intent(inout) :: sdc

    if (istep_first < 0) then
       istep_first = istep
    end if
    istep_this = istep
    trans_called = .false.
    Mach_computed = .false.

    dt = dtio

    select case (method)
    case (SMC_ADVANCE_SDC)
       call advance_sdc(U,courno,dx,sdc,istep==istep_first)
    case (SMC_ADVANCE_MRSDC)
       call advance_multi_sdc(U,courno,dx,sdc,istep==istep_first)
    case (SMC_ADVANCE_RK)
       if (rk_order == 3) then
          call advance_rk33(U,courno,dx)
       else if (rk_order == 4) then
          call advance_rk64(U,courno,dx)
       else
          call bl_error("Invalid rk_order.")
       end if
    case (SMC_ADVANCE_CUSTOM)
       call bl_error("Not implemented yet.")
    end select

    if (check_nans) then
       if (contains_nan(U)) then
          call bl_error("U contains nan")
       end if
    end if

    dtio = dt

  end subroutine advance

  !
  ! Advance U using SSP RK3
  !
  subroutine advance_rk33 (U,courno,dx)
    use time_module, only : time
    use smcdata_module, only : Unew, Uprime
    implicit none

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: courno
    double precision,  intent(in   ) :: dx(3)

    type(bl_prof_timer), save :: bpt_rkstep1, bpt_rkstep2, bpt_rkstep3

    call tb_multifab_setval(Unew, 0.d0, .true.)

    call build(bpt_rkstep1, "rkstep1")
    call dUdt(U, Uprime, time, dt, dx, courno=courno)
    call update_rk(Zero,Unew, One,U, dt, Uprime)
    call reset_density(Unew)
    call impose_hard_bc(Unew, time+OneThird*dt, dx)
    call destroy(bpt_rkstep1)

    call build(bpt_rkstep2, "rkstep2")
    call dUdt(Unew, Uprime, time+OneThird*dt, OneThird*dt, dx)
    call update_rk(OneQuarter, Unew, ThreeQuarters, U, OneQuarter*dt, Uprime)
    call reset_density(Unew)
    call impose_hard_bc(Unew, time+TwoThirds*dt, dx)
    call destroy(bpt_rkstep2)

    call build(bpt_rkstep3, "rkstep3")
    call dUdt(Unew, Uprime, time+TwoThirds*dt, TwoThirds*dt, dx)
    call update_rk(OneThird, U, TwoThirds, Unew, TwoThirds*dt, Uprime)
    call reset_density(U)
    call impose_hard_bc(U, time+dt, dx)
    call destroy(bpt_rkstep3)

  end subroutine advance_rk33

  !
  ! Advance U using 6-stage RK4 (same as S3D)
  !
  subroutine advance_rk64 (U,courno,dx)
    use time_module, only : time
    use smcdata_module, only : Unew, Uprime
    implicit none

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: courno
    double precision,  intent(in   ) :: dx(3)

    double precision, parameter :: rk64_alpha(6) = [ &
         +  3296351145737.d0/15110423921029.d0, &
         +  1879360555526.d0/ 7321162733569.d0, &
         + 10797097731880.d0/20472212111779.d0, &
         +   754636544611.d0/15563872110659.d0, &
         +  3260218886217.d0/ 2618290685819.d0, &
         +  5069185909380.d0/12292927838509.d0 ]

    double precision, parameter :: rk64_beta(6) = [ &
         -  1204558336989.d0/10607789004752.d0, &
         -  3028468927040.d0/14078136890693.d0, &
         -   455570672869.d0/ 8930094212428.d0, &
         - 17275898420483.d0/15997285579755.d0, &
         -  2453906524165.d0/ 9868353053862.d0, &
         +      0.d0 ]

    double precision, parameter :: rk64_time(6) = [ &
         +  3296351145737.d0/15110423921029.d0, &
         +  2703592154963.d0/ 7482974295227.d0, &
         +  7876405563010.d0/11693293791847.d0, &
         + 12920213460229.d0/19253602032679.d0, &
         +  7527523127717.d0/ 9001003553970.d0, &
         +      1.d0 ]

    character(len=7), save :: bpt_names(6) = [ &
         "rkstep1", "rkstep2", "rkstep3", "rkstep4", "rkstep5", "rkstep6" ]

    type(bl_prof_timer), save :: bpt_rkstep(6)
    double precision :: t
    integer :: j

    ! init carry-over
    call copy(Unew, U)

    do j = 1, 6
       call build(bpt_rkstep(j), bpt_names(j))
       if (j == 1) then
          call dUdt(U, Uprime, time, rk64_time(j)*dt, dx, courno=courno)
       else
          call dUdt(U, Uprime, t, rk64_time(j)*dt, dx)
       end if

       t = time + rk64_time(j) * dt

       ! update solution
       call update_rk(0.0d0, U, 1.0d0, Unew, rk64_alpha(j)*dt, Uprime)
       call reset_density(U)
       call impose_hard_bc(U, t, dx)

       ! update carry-over
       if (j < 6) then
          call update_rk(0.0d0, Unew, 1.0d0, U, rk64_beta(j)*dt, Uprime)
          call reset_density(Unew)
          call impose_hard_bc(Unew, t, dx)
       end if

       call destroy(bpt_rkstep(j))
    end do

    call reset_density(U)
    call impose_hard_bc(U, time+dt, dx)

  end subroutine advance_rk64

  !
  ! Advance U using single-rate SDC time-stepping
  !
  subroutine advance_sdc(U, courno, dx, sdc, first_step)

    use time_module, only : time
    use smcdata_module, only : Q
    use probin_module, only : cfl_int, fixed_dt

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: courno
    double precision,  intent(in   ) :: dx(3)
    type(sdc_ctx),     intent(inout) :: sdc
    logical,           intent(in   ) :: first_step

    logical :: update_courno
    double precision :: courno_proc

    integer :: k
    double precision :: res0, res1
    type(layout) :: la
    type(multifab), target :: R

    type(bl_prof_timer), save :: bpt_sdc_prep, bpt_sdc_iter

    !
    ! set dt
    !

    ! this really belongs in the first feval of an sdc sweep

    update_courno = .false.
    if (fixed_dt.le.0.d0) then
       if (mod(istep_this,cfl_int).eq.1 .or. cfl_int.le.1) then
          update_courno = .true.
       end if
    end if

    if (update_courno) then
       courno_proc = -1.d50
       call ctoprim(U, Q, 0)
       call compute_courno(Q, dx, courno_proc)
       call parallel_reduce(courno, courno_proc, MPI_MAX)
    end if

    call set_dt(courno, istep_this)

    !
    ! advance (pass control to sdclib)
    !
    call build(bpt_sdc_prep, "sdc_prep")
    if (first_step) then
       call sdc_imex_set_q0(sdc%imex, mfptr(U))
       call sdc_imex_spread(sdc%imex, time)
    else
       call sdc_imex_spread_qend(sdc%imex)
    end if
    call destroy(bpt_sdc_prep)

    if (sdc%tol_residual > 0.d0) then
       la = get_layout(U)
       call build(R, la, ncons, 0)
    end if

    res0 = -1.0d0

    call build(bpt_sdc_iter, "sdc_iter")
    do k = 1, sdc%iters
       call sdc_imex_sweep(sdc%imex, time, dt, 0)

       ! check residual
       if (sdc%tol_residual > 0.d0) then
          call sdc_imex_residual(sdc%imex, dt, mfptr(R))
          res1 = multifab_norm_l2(R)

          if (parallel_IOProcessor()) then
             if (res0 > 0.0d0) then
                print *, "SDC: iter:", k, "residual:", res1, res0/res1
             else
                print *, "SDC: iter:", k, "residual:", res1
             end if
          end if

          if (res0 > 0.0d0) then
             if (abs(res0 / res1 - 1.0d0) < sdc%tol_residual) &
               exit
          end if

          res0 = res1
       end if
    end do
    call destroy(bpt_sdc_iter)

    call sdc_imex_get_qend(sdc%imex, mfptr(U))

    call reset_density(U)
    call impose_hard_bc(U, time+dt, dx)

    if (sdc%tol_residual > 0.d0) then
       call destroy(R)
    end if

  end subroutine advance_sdc


  !
  ! SDCLib callbacks
  !
  subroutine single_sdc_feval(Fptr, Uptr, t, state, ctxptr) bind(c)
    type(c_ptr),     intent(in), value :: Fptr, Uptr, ctxptr
    type(sdc_state), intent(in)        :: state
    real(c_double),  intent(in), value :: t

    type(multifab), pointer :: U, Uprime
    type(sdc_ctx),  pointer :: ctx

    type(sdc_nset), pointer :: nset
    real(c_double), pointer :: nodes(:)
    real(c_double)          :: dt_m
    integer                 :: node

    call c_f_pointer(Uptr, U)
    call c_f_pointer(Fptr, Uprime)
    call c_f_pointer(ctxptr, ctx)

    ! hack: compute sub-step dt and wrap around appropriately
    call c_f_pointer(ctx%imex%nset, nset)
    call c_f_pointer(nset%nodes, nodes, [ nset%nnodes ])

    node = state%node + 1
    if (node >= nset%nnodes) then
       dt_m = dt * (nodes(2) - nodes(1))
    else
       dt_m = dt * (nodes(node+1) - nodes(node))
    end if

    call dUdt(U, Uprime, t, dt_m, ctx%dx)
  end subroutine single_sdc_feval

  subroutine sdc_post_step_cb(Uptr, state, ctxptr) bind(c)
    type(c_ptr),     intent(in), value :: Uptr, ctxptr
    type(sdc_state), intent(in)        :: state

    type(multifab), pointer :: U
    type(sdc_ctx),  pointer :: ctx

    call c_f_pointer(Uptr, U)
    call c_f_pointer(ctxptr, ctx)

    call reset_density(U)
    call impose_hard_bc(U,state%t,ctx%dx)
  end subroutine sdc_post_step_cb

  subroutine multi_sdc_feval_slow(Fptr, Uptr, t, state, ctxptr) bind(c)
    use probin_module, only : sdc_multirate_explicit

    type(c_ptr),     intent(in), value :: Fptr, Uptr, ctxptr
    type(sdc_state), intent(in)        :: state
    real(c_double),  intent(in), value :: t

    type(multifab), pointer :: U, Uprime, Uprime_chem
    type(sdc_ctx),  pointer :: ctx

    type(sdc_nodes), pointer :: nds
    real(c_double),  pointer :: nodes(:)
    real(c_double)           :: dt_m
    integer                  :: node

    call c_f_pointer(Uptr, U)
    call c_f_pointer(Fptr, Uprime)
    call c_f_pointer(ctxptr, ctx)

    ! hack: compute sub-step dt and wrap around appropriately
    nds => ctx%nodes1
    call c_f_pointer(nds%nodes, nodes, [ nds%nnodes ])

    node = state%node + 1
    if (node >= nds%nnodes) then
       dt_m = dt * (nodes(2) - nodes(1))
    else
       dt_m = dt * (nodes(node+1) - nodes(node))
    end if

    if (sdc_multirate_explicit) then
       Uprime_chem => sdc_get_chemterm(ctx, state%node)
       call dUdt(U, Uprime, t, dt_m, ctx%dx, include_r=.false., Uprime_c=Uprime_chem)
    else
       call dUdt(U, Uprime, t, dt_m, ctx%dx, include_ad=.false.)
    end if
  end subroutine multi_sdc_feval_slow

  subroutine multi_sdc_feval_fast(Fptr, Uptr, t, state, ctxptr) bind(c)
    use probin_module, only : sdc_multirate_explicit
    type(c_ptr),     intent(in), value :: Fptr, Uptr, ctxptr
    type(sdc_state), intent(in)        :: state
    real(c_double),  intent(in), value :: t

    type(multifab), pointer :: U, Uprime, Uprime_chem
    type(sdc_ctx),  pointer :: ctx

    type(sdc_nset), pointer :: nset
    real(c_double), pointer :: nodes(:)
    real(c_double)          :: dt_m
    integer                 :: node

    call c_f_pointer(Uptr, U)
    call c_f_pointer(Fptr, Uprime)
    call c_f_pointer(ctxptr, ctx)

    ! hack: compute sub-step dt and wrap around appropriately
    call c_f_pointer(ctx%mrex%nset, nset)
    call c_f_pointer(nset%nodes, nodes, [ nset%nnodes ])

    node = state%node + 1
    if (node >= nset%nnodes) then
       dt_m = dt * (nodes(2) - nodes(1))
    else
       dt_m = dt * (nodes(node+1) - nodes(node))
    end if

    if (sdc_multirate_explicit) then
       call dUdt(U, Uprime, t, dt_m, ctx%dx, include_ad=.false.)
    else
       Uprime_chem => sdc_get_chemterm(ctx, state%node)
       call dUdt(U, Uprime, t, dt_m, ctx%dx, include_r=.false., Uprime_c=Uprime_chem)
    end if
  end subroutine multi_sdc_feval_fast


  !
  ! Advance U using multi-rate SDC time-stepping
  !
  subroutine advance_multi_sdc(U, courno, dx, sdc, first_step)

    use time_module, only : time
    use smcdata_module, only : Q
    use probin_module, only : cfl_int, fixed_dt

    type(multifab),   intent(inout) :: U
    double precision, intent(inout) :: courno
    double precision, intent(in   ) :: dx(3)
    type(sdc_ctx),    intent(inout) :: sdc
    logical,          intent(in   ) :: first_step

    logical :: update_courno
    double precision :: courno_proc

    integer :: k
    double precision :: res0, res1
    type(layout) :: la
    type(multifab), target :: R

    type(bl_prof_timer), save :: bpt_sdc_prep, bpt_sdc_iter


    ! ideally we would have a preallocated workspace for the residual,
    ! and the computation of dt would be done in the first feval...

    !
    ! set dt
    !

    ! this really belongs in the first feval of an sdc sweep

    update_courno = .false.
    if (fixed_dt.le.0.d0) then
       if (mod(istep_this,cfl_int).eq.1 .or. cfl_int.le.1) then
          update_courno = .true.
       end if
    end if

    if (update_courno) then
       courno_proc = -1.d50
       call ctoprim(U, Q, 0)
       call compute_courno(Q, dx, courno_proc)
       call parallel_reduce(courno, courno_proc, MPI_MAX)
    end if

    call set_dt(courno, istep_this)

    !
    ! advance
    !
    call build(bpt_sdc_prep, "sdc_prep")
    if (first_step) then
       call sdc_mrex_set_q0(sdc%mrex, mfptr(U))
       call sdc_mrex_spread(sdc%mrex, time)
    else
       call sdc_mrex_spread_qend(sdc%mrex)
    end if
    call destroy(bpt_sdc_prep)

    if (sdc%tol_residual > 0.d0) then
       la = get_layout(U)
       call build(R, la, ncons, 0)
    end if

    res0 = -1.0d0

    call build(bpt_sdc_iter, "sdc_iter")
    do k = 1, sdc%iters
       call sdc_mrex_sweep(sdc%mrex, time, dt, 0);

       ! check residual
       if (sdc%tol_residual > 0.d0) then
          call sdc_mrex_residual(sdc%mrex, dt, mfptr(R))
          res1 = multifab_norm_l2(R)

          if (parallel_IOProcessor()) then
             if (res0 > 0.0d0) then
                print *, "SDC: iter:", k, "residual:", res1, res0/res1
             else
                print *, "SDC: iter:", k, "residual:", res1
             end if
          end if

          if (res0 > 0.0d0) then
             if (abs(res0 / res1 - 1.0d0) < sdc%tol_residual) &
               exit
          end if

          res0 = res1
       end if
    end do

    call sdc_mrex_get_qend(sdc%mrex, mfptr(U))

    call reset_density(U)
    call impose_hard_bc(U, time+dt, dx)

    if (sdc%tol_residual > 0.d0) then
       call destroy(R)
    end if
  end subroutine advance_multi_sdc


  !
  ! Compute new time-step size
  !
  subroutine set_dt(courno, istep)

    use time_module, only : time
    use probin_module, only : fixed_dt, cflfac, init_shrink, max_dt_growth, &
         max_dt, small_dt, stop_time

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
  subroutine update_rk(a,U1,b,U2,c,Uprime)

    type(multifab),   intent(in   ) :: U2, Uprime
    type(multifab),   intent(inout) :: U1
    double precision, intent(in   ) :: a, b, c

    integer :: lo(3), hi(3), i, j, k, m, n, nc, dm
    double precision, pointer, dimension(:,:,:,:) :: u1p, u2p, upp

    dm = U1%dim
    nc = ncomp(U1)

    !$omp parallel private(i,j,k,m,n,lo,hi,u1p,u2p,upp)
    lo(3) = 1
    hi(3) = 1
    do n=1,nfabs(U1)

       if (.not.tb_worktodo(n)) cycle

       u1p => dataptr(U1,    n)
       u2p => dataptr(U2,    n)
       upp => dataptr(Uprime,n)

       lo(1:dm) = tb_get_valid_lo(n)
       hi(1:dm) = tb_get_valid_hi(n)

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

  end subroutine update_rk


  !
  ! Compute dU/dt given U.
  !
  ! The Courant number (courno) is also computed if passed.
  !
  subroutine dUdt (U, Uprime, t, dt_m, dx, courno, include_ad, include_r, Uprime_c)

    use derivative_stencil_module, only : stencil, narrow, s3d
    use smcdata_module, only : Q, mu, xi, lam, Ddiag, Fdif, Upchem, qx, qy, qz
    use probin_module, only : overlap_comm_comp, overlap_comm_gettrans, cfl_int, fixed_dt, &
         trans_int, mach_int

    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: t, dt_m, dx(3)
    double precision, intent(inout), optional :: courno
    logical,          intent(in   ), optional :: include_ad, include_r
    type(multifab),   intent(in   ), optional :: Uprime_c

    integer ::    lo(U%dim),    hi(U%dim)
    integer ::   dlo(U%dim),   dhi(U%dim)
    integer ::   blo(U%dim),   bhi(U%dim)
    integer :: dm, n, ng, iblock
    integer :: ng_ctoprim, ng_gettrans

    logical :: update_trans, update_mach

    logical :: update_courno
    double precision :: courno_proc

    type(mf_fb_data) :: U_fb_data, qx_fb_data, qy_fb_data, qz_fb_data

    logical :: inc_ad, inc_r, rYt_only

    integer :: qlo(4), qhi(4), uplo(4), uphi(4), ulo(4), uhi(4), flo(4), fhi(4), &
         upclo(4), upchi(4), qxlo(4), qxhi(4), qylo(4), qyhi(4), qzlo(4), qzhi(4)
    double precision, pointer, dimension(:,:,:,:) :: up, qp, mup, xip, lamp, Ddp, upp, fp, &
         upcp, qxp, qyp, qzp

    type(bl_prof_timer), save :: bpt_ctoprim, bpt_gettrans, bpt_hypdiffterm
    type(bl_prof_timer), save :: bpt_chemterm, bpt_nscbc

    dm = U%dim
    ng = nghost(U)

    inc_ad = .true.; if (present(include_ad)) inc_ad = include_ad
    inc_r  = .true.; if (present(include_r))  inc_r  = include_r

    update_courno = .false.
    if (present(courno) .and. fixed_dt.le.0.d0) then
       if (mod(istep_this,cfl_int).eq.1 .or. cfl_int.le.1) then
          update_courno = .true.
       end if
    end if

    if (trans_int .le. 0) then
       update_trans = .true.
    else
       if (trans_called) then
          update_trans = .false.
       else if (mod((istep_this-istep_first), trans_int) .eq. 0) then
          update_trans = .true.
       else
          update_trans = .false.
       end if
    end if

    if (mach_int .le. 0) then
       update_mach = .true.
    else
       if (Mach_computed) then
          update_mach = .false.
       else if (mod((istep_this-istep_first), Mach_int) .eq. 0) then
          update_mach = .true.
       else
          update_mach = .false.
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

    if (inc_ad) count_ad = count_ad + 1
    if (inc_r)  count_r  = count_r  + 1

    !
    ! Calculate primitive variables based on U
    !
    call build(bpt_ctoprim, "ctoprim")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    rYt_only = .not. (inc_ad .or. update_courno)
    call ctoprim(U, Q, ng_ctoprim, rYT_only=rYt_only)
    call destroy(bpt_ctoprim)            !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    if (inc_ad .and. overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    if (update_courno) then
       courno_proc = -1.d50
       call compute_courno(Q, dx, courno_proc)
       call parallel_reduce(courno, courno_proc, MPI_MAX)
    end if

    if (present(courno)) then
       call set_dt(courno, istep_this)
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
       !$omp parallel private(n,qp,upp,upcp,qlo,qhi,uplo,uphi,lo,hi,upclo,upchi)
       do n=1,nfabs(Q)

          if (.not.tb_worktodo(n)) cycle

          qp   => dataptr(Q,n)
          upp  => dataptr(Uprime,n)
          upcp => dataptr(Upchem,n)

          qlo = lbound(qp)
          qhi = ubound(qp)
          uplo = lbound(upp)
          uphi = ubound(upp)
          upclo = lbound(upcp)
          upchi = ubound(upcp)

          lo = tb_get_valid_lo(n)
          hi = tb_get_valid_hi(n)

          if (dm .eq. 2) then
             call chemterm_2d(lo,hi,qp,qlo(1:2),qhi(1:2),upp,uplo(1:2),uphi(1:2), &
                  upcp,upclo(1:2),upchi(1:2), dt_m)
          else
             call chemterm_3d(lo,hi,qp,qlo(1:3),qhi(1:3),upp,uplo(1:3),uphi(1:3), &
                  upcp,upclo(1:3),upchi(1:3), dt_m)
          end if
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
       if (update_trans) then
          call build(bpt_gettrans, "gettrans")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
          call get_transport_properties(Q, mu, xi, lam, Ddiag, ng_gettrans)
          trans_called = .true.
          call destroy(bpt_gettrans)               !! ^^^^^^^^^^^^^^^^^^^^^^^ timer
       end if

       if (overlap_comm_comp) then
          call multifab_fill_boundary_waitrecv(U, U_fb_data)

          if (ng_ctoprim .eq. 0) then
             call build(bpt_ctoprim, "ctoprim")    !! vvvvvvvvvvvvvvvvvvvvvvv timer
             call ctoprim(U, Q, ghostcells_only=.true.)
             call destroy(bpt_ctoprim)             !! ^^^^^^^^^^^^^^^^^^^^^^^ timer
          end if

          if (ng_gettrans .eq. 0 .and. update_trans) then
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
       if (stencil .eq. narrow) then

          !$omp parallel private(n,iblock,lo,hi,up,ulo,uhi,upp,uplo,uphi,qp,qlo,qhi) &
          !$omp private(fp,flo,fhi,mup,xip,lamp,Ddp,dlo,dhi,blo,bhi)
          do n=1,nfabs(Q)

             if (.not.tb_worktodo(n)) cycle

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

             call get_data_lo_hi(n,dlo,dhi)
             call get_boxbc(n,blo,bhi)

             do iblock = 1, tb_get_nblocks(n)
                lo = tb_get_block_lo(iblock,n)
                hi = tb_get_block_hi(iblock,n)

                if (dm .eq. 2) then
                   call hypterm_2d(lo,hi,dx,up,ulo(1:2),uhi(1:2),qp,qlo(1:2),qhi(1:2),&
                        upp,uplo(1:2),uphi(1:2),dlo,dhi,blo,bhi)

                   call narrow_diffterm_2d(lo,hi,dx,qp,qlo(1:2),qhi(1:2),upp,uplo(1:2),uphi(1:2), &
                        fp,flo(1:2),fhi(1:2),mup,xip,lamp,Ddp,dlo,dhi,blo,bhi)
                else
                   call hypterm_3d(lo,hi,dx,up,ulo(1:3),uhi(1:3),qp,qlo(1:3),qhi(1:3),&
                        upp,uplo(1:3),uphi(1:3),dlo,dhi,blo,bhi)

                   call narrow_diffterm_3d(lo,hi,dx,qp,qlo(1:3),qhi(1:3),upp,uplo(1:3),uphi(1:3), &
                        fp,flo(1:3),fhi(1:3),mup,xip,lamp,Ddp,dlo,dhi,blo,bhi)
                end if
             end do

          end do
          !$omp end parallel

       else if (stencil .eq. s3d) then

          !$omp parallel private(n,iblock,lo,hi,qp,qlo,qhi) &
          !$omp private(qxp,qxlo,qxhi,qyp,qylo,qyhi,qzp,qzlo,qzhi) &
          !$omp private(fp,flo,fhi,mup,xip,dlo,dhi,blo,bhi)
          do n=1,nfabs(Q)
             if (.not.tb_worktodo(n)) cycle

             qp => dataptr(Q,n)
             qxp => dataptr(qx, n)
             qyp => dataptr(qy, n)
             if (dm .eq. 3) qzp => dataptr(qz, n)

             fp => dataptr(Fdif,n)

             mup  => dataptr(mu   , n)
             xip  => dataptr(xi   , n)

             qlo = lbound(qp)
             qhi = ubound(qp)
             qxlo = lbound(qxp)
             qxhi = ubound(qxp)
             qylo = lbound(qyp)
             qyhi = ubound(qyp)
             if (dm .eq. 3) then
                qzlo = lbound(qzp)
                qzhi = ubound(qzp)
             end if

             flo = lbound(fp)
             fhi = ubound(fp)

             call get_data_lo_hi(n,dlo,dhi)
             call get_boxbc(n,blo,bhi)

             do iblock = 1, tb_get_nblocks(n)
                lo = tb_get_block_lo(iblock,n)
                hi = tb_get_block_hi(iblock,n)

                if (dm .eq. 2) then
                   call bl_error("2D not supported for S3D mode")
                else
                   call s3d_diffterm_1_3d(lo,hi,dx,qp,qlo(1:3),qhi(1:3),  &
                        fp,flo(1:3),fhi(1:3),    &
                        qxp,qxlo(1:3),qxhi(1:3), &
                        qyp,qylo(1:3),qyhi(1:3), &
                        qzp,qzlo(1:3),qzhi(1:3), &
                        mup,xip,dlo,dhi,blo,bhi)
                end if
             end do
          end do
          !$omp end parallel

          call multifab_fill_boundary_nowait(qx, qx_fb_data, idim=1)
          call multifab_fill_boundary_nowait(qy, qy_fb_data, idim=2)
          if (dm .eq. 3) then
             call multifab_fill_boundary_nowait(qz, qz_fb_data, idim=3)
          end if
          if (overlap_comm_comp) then
             call multifab_fill_boundary_test(qx, qx_fb_data, idim=1)
             call multifab_fill_boundary_test(qy, qy_fb_data, idim=2)
             if (dm .eq. 3) then
                call multifab_fill_boundary_test(qz, qz_fb_data, idim=3)
             end if
          else
             call multifab_fill_boundary_finish(qx, qx_fb_data, idim=1)
             call multifab_fill_boundary_finish(qy, qy_fb_data, idim=2)
             if (dm .eq. 3) then
                call multifab_fill_boundary_finish(qz, qz_fb_data, idim=3)
             end if
          end if

          !$omp parallel private(n,iblock,lo,hi,up,ulo,uhi,upp,uplo,uphi,qp,qlo,qhi) &
          !$omp private(dlo,dhi,blo,bhi)
          do n=1,nfabs(Q)
             if (.not.tb_worktodo(n)) cycle

             up => dataptr(U,n)
             upp=> dataptr(Uprime,n)
             qp => dataptr(Q,n)

             ulo = lbound(up)
             uhi = ubound(up)
             qlo = lbound(qp)
             qhi = ubound(qp)
             uplo = lbound(upp)
             uphi = ubound(upp)

             call get_data_lo_hi(n,dlo,dhi)
             call get_boxbc(n,blo,bhi)

             do iblock = 1, tb_get_nblocks(n)
                lo = tb_get_block_lo(iblock,n)
                hi = tb_get_block_hi(iblock,n)

                if (dm .eq. 2) then
                   call bl_error("2D not supported for S3D mode")
                else
                   call hypterm_3d(lo,hi,dx,up,ulo(1:3),uhi(1:3),qp,qlo(1:3),qhi(1:3),&
                        upp,uplo(1:3),uphi(1:3),dlo,dhi,blo,bhi)
                end if
             end do
          end do
          !$omp end parallel

          if (overlap_comm_comp) then
             call multifab_fill_boundary_finish(qx, qx_fb_data, idim=1)
             call multifab_fill_boundary_finish(qy, qy_fb_data, idim=2)
             if (dm .eq. 3) then
                call multifab_fill_boundary_finish(qz, qz_fb_data, idim=3)
             end if
          end if

          !$omp parallel private(n,iblock,lo,hi,upp,uplo,uphi,qp,qlo,qhi) &
          !$omp private(qxp,qxlo,qxhi,qyp,qylo,qyhi,qzp,qzlo,qzhi) &
          !$omp private(fp,flo,fhi,mup,xip,lamp,Ddp,dlo,dhi,blo,bhi)
          do n=1,nfabs(Q)
             if (.not.tb_worktodo(n)) cycle

             qp  => dataptr(Q , n)
             qxp => dataptr(qx, n)
             qyp => dataptr(qy, n)
             if (dm .eq. 3) qzp => dataptr(qz, n)

             upp=> dataptr(Uprime,n)
             fp => dataptr(Fdif,n)

             mup  => dataptr(mu   , n)
             xip  => dataptr(xi   , n)
             lamp => dataptr(lam  , n)
             Ddp  => dataptr(Ddiag, n)

             qlo = lbound(qp)
             qhi = ubound(qp)
             qxlo = lbound(qxp)
             qxhi = ubound(qxp)
             qylo = lbound(qyp)
             qyhi = ubound(qyp)
             if (dm .eq. 3) then
                qzlo = lbound(qzp)
                qzhi = ubound(qzp)
             end if

             uplo = lbound(upp)
             uphi = ubound(upp)
             flo = lbound(fp)
             fhi = ubound(fp)

             call get_data_lo_hi(n,dlo,dhi)
             call get_boxbc(n,blo,bhi)

             do iblock = 1, tb_get_nblocks(n)
                lo = tb_get_block_lo(iblock,n)
                hi = tb_get_block_hi(iblock,n)

                if (dm .eq. 2) then
                   call bl_error("2D not supported for S3D mode")
                else
                   call s3d_diffterm_2_3d(lo,hi,dx,qp,qlo(1:3),qhi(1:3),  &
                        upp,uplo(1:3),uphi(1:3), fp,flo(1:3),fhi(1:3),    &
                        qxp,qxlo(1:3),qxhi(1:3), &
                        qyp,qylo(1:3),qyhi(1:3), &
                        qzp,qzlo(1:3),qzhi(1:3), &
                        mup,xip,lamp, Ddp, &
                        dlo,dhi,blo,bhi)
                end if
             end do
          end do
          !$omp end parallel

       else
          call bl_error("dUdt: unknown stencil type")
       end if
       call destroy(bpt_hypdiffterm)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

       !
       ! NSCBC boundary
       !
       call build(bpt_nscbc, "nscbc")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
       if (present(Uprime_c)) then
          call nscbc(Q, U, Fdif, Uprime_c, Uprime, t, dx, update_mach)
       else
          call nscbc(Q, U, Fdif, Upchem  , Uprime, t, dx, update_mach)
       end if
       if (update_mach) mach_computed = .true.
       call destroy(bpt_nscbc)          !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    end if

  end subroutine dUdt


  subroutine compute_courno(Q, dx, courno)
    type(multifab), intent(in) :: Q
    double precision, intent(in) :: dx(3)
    double precision, intent(inout) :: courno

    integer :: dm, n, lo(Q%dim), hi(Q%dim), qlo(4), qhi(4)
    double precision :: courno_thread
    double precision, pointer :: qp(:,:,:,:)

    dm = Q%dim

    !$omp parallel private(n, lo, hi, qlo, qhi, qp, courno_thread) &
    !$omp reduction(max:courno)
    do n=1,nfabs(Q)

       if (.not.tb_worktodo(n)) cycle

       qp => dataptr(Q,n)
       qlo = lbound(qp)
       qhi = ubound(qp)

       lo = tb_get_valid_lo(n)
       hi = tb_get_valid_hi(n)

       courno_thread = 0.d0

       if (dm .eq. 2) then
          call comp_courno_2d(lo,hi,dx,qp,qlo(1:2),qhi(1:2),courno_thread)
       else
          call comp_courno_3d(lo,hi,dx,qp,qlo(1:3),qhi(1:3),courno_thread)
       end if

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

    integer :: dm, ng, ng_ctoprim, ng_gettrans, n, lo(U%dim), hi(U%dim)
    integer :: qlo(4), qhi(4), uplo(4), uphi(4), upclo(4), upchi(4)
    type(layout)     :: la
    type(multifab)   :: Q, Uprime, Upchem, mu, xi, lam, Ddiag
    double precision, pointer, dimension(:,:,:,:) :: qp, upp, upcp

    call multifab_fill_boundary_test(U, U_fb_data)

    dm = U%dim
    ng = nghost(U)
    la = get_layout(U)

    call multifab_build(Q, la, nprim, ng)
    call tb_multifab_setval(Q, 0.d0, .true.)

    call multifab_build(Uprime, la, ncons, 0)
    call tb_multifab_setval(Uprime, 0.d0)

    call multifab_build(Upchem, la, nspecies, 0)

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

    !$omp parallel private(n,qp,upp,qlo,qhi,uplo,uphi,lo,hi)
    do n=1,nfabs(Q)

       if (.not.tb_worktodo(n)) cycle

       qp  => dataptr(Q,n)
       upp => dataptr(Uprime,n)
       upcp => dataptr(Upchem,n)

       qlo = lbound(qp)
       qhi = ubound(qp)
       uplo = lbound(upp)
       uphi = ubound(upp)
       upclo = lbound(upcp)
       upchi = ubound(upcp)

       lo = tb_get_valid_lo(n)
       hi = tb_get_valid_hi(n)

       dt = 1.d-10
       if (dm .eq. 2) then
          call chemterm_2d(lo,hi,qp,qlo(1:2),qhi(1:2),upp,uplo(1:2),uphi(1:2), &
               upcp,upclo(1:2),upchi(1:2), dt)
       else
          call chemterm_3d(lo,hi,qp,qlo(1:3),qhi(1:3),upp,uplo(1:3),uphi(1:3), &
               upcp,upclo(1:3),upchi(1:3), dt)
       end if
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

  subroutine sdc_get_q0(U0, sdc)
    use probin_module
    type(multifab), intent(inout) :: U0
    type(sdc_ctx),  intent(inout) :: sdc

    if (method .eq. SMC_ADVANCE_SDC) then
       call sdc_imex_get_q0(sdc%imex, mfptr(U0))
    else
       call sdc_mrex_get_q0(sdc%mrex, mfptr(U0))
    end if
  end subroutine sdc_get_q0

end module advance_module
