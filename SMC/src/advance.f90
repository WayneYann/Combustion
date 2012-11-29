module advance_module

  use bl_error_module
  use derivative_stencil_module
  use kernels_module
  use multifab_module
  use omp_module
  use nscbc_module
  use sdcquad_module
  use smc_bc_module
  use time_module
  use transport_properties
  use variables_module

  use chemistry_module, only : nspecies

  implicit none

  private
  public advance

contains

  subroutine advance(U, dt, courno, dx, sdc, istep)

    use probin_module, only : advance_method

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt, courno
    double precision,  intent(in   ) :: dx(U%dim)
    integer,           intent(in   ) :: istep
    type(sdcquad),     intent(in   ) :: sdc

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
    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt, courno
    double precision,  intent(in   ) :: dx(U%dim)
    integer, intent(in) :: istep

    integer          :: ng
    type(layout)     :: la
    type(multifab)   :: Uprime, Unew

    type(bl_prof_timer), save :: bpt_rkstep1, bpt_rkstep2, bpt_rkstep3, bpt_setdt

    ng = nghost(U)
    la = get_layout(U)

    call multifab_build(Uprime, la, ncons, 0)
    call multifab_build(Unew,   la, ncons, ng)
    call multifab_setval(Unew, 0.d0, .true.)

    ! RK Step 1
    call build(bpt_rkstep1, "rkstep1")   !! vvvvvvvvvvvvvvvvvvvvvvv timer

    call dUdt(U, Uprime, dx, courno=courno, istep=istep)
    call build(bpt_setdt, "setdt")
    call set_dt(dt, courno, istep)
    call destroy(bpt_setdt)
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

    call destroy(Unew)
    call destroy(Uprime)

  end subroutine advance_rk3

  !
  ! Advance U using SDC time-stepping
  !
  subroutine advance_sdc(U, dt, courno, dx, sdc, istep)

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt, courno
    double precision,  intent(in   ) :: dx(U%dim)
    type(sdcquad),     intent(in   ) :: sdc
    integer,           intent(in   ) :: istep

    integer          :: k, m, n, ng
    double precision :: res_proc, res
    type(layout)     :: la
    type(multifab)   :: uSDC(sdc%nnodes), fSDC(sdc%nnodes), S(sdc%nnodes-1)

    double precision :: dtsdc(sdc%nnodes-1)

    ng = nghost(U)
    la = get_layout(U)

    ! build u and u' multifabs for each node
    do m = 1, sdc%nnodes
       call build(uSDC(m), la, ncons, ng)
       call build(fSDC(m), la, ncons, 0)
    end do

    ! build S multifab (node to node integrals)
    do m = 1, sdc%nnodes-1
       call build(S(m), la, ncons, 0)
    end do

    ! set provisional solution, compute dt

    call copy(uSDC(1), U)
    call dUdt(uSDC(1), fSDC(1), dx, courno, istep)

    do m = 2, sdc%nnodes
       call copy(uSDC(m), uSDC(1))
       call copy(fSDC(m), fSDC(1))
    end do

    call set_dt(dt, courno, istep)

    ! perform sdc iterations
    res = 0.0d0
    dtsdc = dt * (sdc%nodes(2:sdc%nnodes) - sdc%nodes(1:sdc%nnodes-1))

    do k = 1, sdc%iters

       ! compute integrals (compact forward euler)
       do m = 1, sdc%nnodes-1
          call setval(S(m), 0.0d0)
          do n = 1, sdc%nnodes
             call saxpy(S(m), sdc%smats(m,n,1), fSDC(n))
          end do
       end do

       ! perform sub-step correction
       do m = 1, sdc%nnodes-1

          ! U(m+1) = U(m) + dt dUdt(m) + dt S(m)

          call copy(uSDC(m+1), uSDC(m))
          call saxpy(uSDC(m+1), dtsdc(m), fSDC(m))
          call saxpy(uSDC(m+1), dt, S(m))
          call reset_density(uSDC(m+1))
          call impose_hard_bc(uSDC(m+1))

          call dUdt(uSDC(m+1), fSDC(m+1), dx)

       end do

       ! check residual
       if (sdc%tol_residual > 0.d0) then
          res_proc = sdc_residual(uSDC, fSDC, S(1), dt, sdc)
          call parallel_reduce(res, res_proc, MPI_MAX)

          if (parallel_IOProcessor()) then
             print *, "SDC: iter:", k, "residual:", res
          end if

          if (res < sdc%tol_residual) exit
       end if
    end do

    call copy(U, uSDC(sdc%nnodes))

    ! destroy
    do m = 1, sdc%nnodes
       call destroy(uSDC(m))
       call destroy(fSDC(m))
    end do

    do m = 1, sdc%nnodes-1
       call destroy(S(m))
    end do

  end subroutine advance_sdc


  !
  ! Compute SDC residual
  !
  function sdc_residual (uSDC,fSDC,R,dt,sdc) result(res)
    real(dp_t)                      :: res
    type(sdcquad),    intent(in   ) :: sdc
    type(multifab),   intent(inout) :: uSDC(sdc%nnodes), fSDC(sdc%nnodes), R
    real(dp_t),       intent(in   ) :: dt

    integer :: m, n

    ! compute integral
    call copy(R, uSDC(1))

    do m = 1, sdc%nnodes-1
       do n = 1, sdc%nnodes
          call saxpy(R, dt*sdc%smat(m,n), fSDC(n))
       end do
    end do

    call saxpy(R, -1.0d0, uSDC(sdc%nnodes))

    res = norm_inf(R)

  end function sdc_residual

  !
  ! Advance U using multi-rate SDC time-stepping
  !
  subroutine advance_multi_sdc(U, dt, courno, dx, sdc, istep)

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt, courno
    double precision,  intent(in   ) :: dx(U%dim)
    type(sdcquad),     intent(in   ) :: sdc
    integer,           intent(in   ) :: istep

    integer          :: k, m, mm, ng
    double precision :: res_proc, res
    type(layout)     :: la
    type(multifab)   :: uAD(sdc%nnodes), fAD(sdc%nnodes), SAD(sdc%nnodes-1)
    type(multifab)   :: uR(sdc%nnodes), fR(sdc%nnodes), SR(sdc%nnodes-1)

    ng = nghost(U)
    la = get_layout(U)

    ! XXX: this is a work in progress
    print *, '*** MULTIRATE SDC IS A WORK IN PROGRESS ***'

    ! XXX: this just does normal SDC for now

    ! build u and u' multifabs for each node
    do m = 1, sdc%nnodes
       call build(uAD(m), la, ncons, ng)
       call build(fAD(m), la, ncons, 0)
       call build(uR(m), la, ncons, ng)
       call build(fR(m), la, ncons, 0)
    end do

    ! build S multifab (node to node integrals)
    do m = 1, sdc%nnodes-1
       call build(SAD(m), la, ncons, 0)
       call build(SR(m), la, ncons, 0)
    end do

    ! set provisional solution, compute dt

    call copy(uAD(1), U)
    call dUdt(uAD(1), fAD(1), dx)

    do m = 2, sdc%nnodes
       call copy(uAD(m), uAD(1))
       call copy(fAD(m), fAD(1))
    end do

    call set_dt(dt, courno, istep)

    ! perform sdc iterations
    res = 0.0d0

    do k = 1, sdc%iters
       ! call sdc_sweep(uAD, fAD, SAD, dx, dt, sdc)
    end do

    call copy(U, uAD(sdc%nnodes))

    ! destroy
    do m = 1, sdc%nnodes
       call destroy(uAD(m))
       call destroy(fAD(m))
       call destroy(uR(m))
       call destroy(fR(m))
    end do

    do m = 1, sdc%nnodes-1
       call destroy(SAD(m))
       call destroy(SR(m))
    end do

  end subroutine advance_multi_sdc


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

    integer :: lo(U1%dim), hi(U1%dim), i, j, k, m, n, nc

    double precision, pointer, dimension(:,:,:,:) :: u1p, u2p, upp

    nc = ncomp(U1)

    do n=1,nfabs(U1)
       u1p => dataptr(U1,    n)
       u2p => dataptr(U2,    n)
       upp => dataptr(Uprime,n)

       lo = lwb(get_box(U1,n))
       hi = upb(get_box(U1,n))

       do m = 1, nc
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   u1p(i,j,k,m) = a * u1p(i,j,k,m) + b * u2p(i,j,k,m) + c * upp(i,j,k,m)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do

  end subroutine update_rk3


  !
  ! Compute dU/dt given U.
  !
  ! The Courant number (courno) might also be computed if passed.
  !
  subroutine dUdt (U, Uprime, dx, courno, istep)
    use derivative_stencil_module, only : stencil, compact, s3d

    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)
    double precision, intent(inout), optional :: courno
    integer,          intent(in   ), optional :: istep

    if (stencil .eq. compact) then
       call dUdt_compact(U, Uprime, dx, courno, istep)
    else if (stencil .eq. s3d) then
       call dUdt_S3D(U, Uprime, dx, courno, istep)
    else
       call bl_error("advance: unknown stencil type")
    end if
  end subroutine dUdt


  !
  ! Compute advection/diffusion part of dU/dt given U.
  !
  subroutine dUdt_AD (U, Uprime, dx)
    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)

    call dUdt_compact(U, Uprime, dx, include_r=.false.)
  end subroutine dUdt_AD


  !
  ! Compute reaction part of dU/dt given U.
  !
  subroutine dUdt_R (U, Uprime, dx)
    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)

    call dUdt_compact(U, Uprime, dx, include_ad=.false.)
  end subroutine dUdt_R


  !
  ! Compute dU/dt given U using the compact stencil.
  !
  ! The Courant number (courno) is also computed if passed.
  !
  subroutine dUdt_compact (U, Uprime, dx, courno, istep, include_ad, include_r)

    use probin_module, only : overlap_comm_comp, overlap_comm_gettrans, cfl_int

    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)
    double precision, intent(inout), optional :: courno
    integer,          intent(in   ), optional :: istep
    logical,          intent(in   ), optional :: include_ad, include_r

    type(multifab) :: mu, xi ! viscosity
    type(multifab) :: lam ! partial thermal conductivity
    type(multifab) :: Ddiag ! diagonal components of rho * Y_k * D

    integer ::    lo(U%dim),    hi(U%dim)
    integer ::   dlo(U%dim),   dhi(U%dim)
    integer ::   blo(U%dim),   bhi(U%dim)
    integer :: i, j, k, m, n, ng, dm
    integer :: ng_ctoprim, ng_gettrans

    logical :: update_courno
    double precision :: courno_proc

    type(layout)     :: la
    type(multifab)   :: Q, Fhyp, Fdif
    type(mf_fb_data) :: U_fb_data

    logical :: inc_ad, inc_r

    double precision, pointer, dimension(:,:,:,:) :: up, fhp, fdp, qp, mup, xip, lamp, Ddp, upp

    type(bl_prof_timer), save :: bpt_ctoprim,  bpt_gettrans, bpt_hypterm, bpt_courno
    type(bl_prof_timer), save :: bpt_diffterm, bpt_calcU, bpt_chemterm, bpt_nscbc

    inc_ad = .true.; if (present(include_ad)) inc_ad = include_ad
    inc_r  = .true.; if (present(include_r))  inc_r  = include_r

    update_courno = .false.
    if (present(courno) .and. present(istep)) then
       if (mod(istep,cfl_int).eq.1 .or. cfl_int.le.1) then
          update_courno = .true.
       end if
    end if

    ! On hopper MPI_Test encourages the overlap of communication and compution.
    ! That's why we have so many calls to multifab_fill_boundary_test.

    if (inc_ad) then
       call multifab_fill_boundary_nowait(U, U_fb_data)
       if (overlap_comm_comp) then
          call multifab_fill_boundary_test(U, U_fb_data)
       else
          call multifab_fill_boundary_finish(U, U_fb_data)
       end if
    end if

    call setval(Uprime, ZERO)

    if (inc_ad .and. overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    dm = U%dim
    ng = nghost(U)
    la = get_layout(U)

    if (inc_ad .and. overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    call multifab_build(Q, la, nprim, ng)

    if (inc_ad) then
       if (overlap_comm_comp) then
          call multifab_fill_boundary_test(U, U_fb_data)
       end if

       call multifab_build(Fhyp, la, ncons, 0)
       call multifab_build(Fdif, la, ncons, 0)

       if (overlap_comm_comp) then
          call multifab_fill_boundary_test(U, U_fb_data)
       end if

       call multifab_build(mu , la, 1, ng)
       call multifab_build(xi , la, 1, ng)
       call multifab_build(lam, la, 1, ng)
       call multifab_build(Ddiag, la, nspecies, ng)

       if (overlap_comm_comp) then
          call multifab_fill_boundary_test(U, U_fb_data)
       end if
    end if

    if (U_fb_data%rcvd) then
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
       call build(bpt_courno, "courno")
       courno_proc = -1.d50
       call compute_courno(Q, dx, courno_proc)
       call destroy(bpt_courno)
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
       do n=1,nfabs(Q)
          qp  => dataptr(Q,n)
          upp => dataptr(Uprime,n)
          
          lo = lwb(get_box(Q,n))
          hi = upb(get_box(Q,n))
          
          if (dm .ne. 3) then
             call bl_error("Only 3D chemsitry_term is supported")
          else
             call chemterm_3d(lo,hi,ng,qp,upp)
          end if
       end do
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
       ! Transport terms
       !
       call build(bpt_diffterm, "diffterm")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
       do n=1,nfabs(Q)
          qp  => dataptr(Q,n)
          fdp => dataptr(Fdif,n)

          mup  => dataptr(mu   , n)
          xip  => dataptr(xi   , n)
          lamp => dataptr(lam  , n)
          Ddp  => dataptr(Ddiag, n)

          lo = lwb(get_box(Q,n))
          hi = upb(get_box(Q,n))

          call get_data_lo_hi(n,dlo,dhi)
          call get_boxbc(n,blo,bhi)

          if (dm .ne. 3) then
             call bl_error("Only 3D compact_diffterm is supported")
          else
             call compact_diffterm_3d(lo,hi,ng,dx,qp,fdp,mup,xip,lamp,Ddp,dlo,dhi,blo,bhi)
          end if
       end do
       call destroy(bpt_diffterm)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

       !
       ! Hyperbolic terms
       !
       call build(bpt_hypterm, "hypterm")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
       do n=1,nfabs(Fhyp)
          up => dataptr(U,n)
          qp => dataptr(Q,n)
          fhp=> dataptr(Fhyp,n)

          lo = lwb(get_box(Fhyp,n))
          hi = upb(get_box(Fhyp,n))

          call get_data_lo_hi(n,dlo,dhi)
          call get_boxbc(n,blo,bhi)

          if (dm .ne. 3) then
             call bl_error("Only 3D hypterm is supported")
          else
             call hypterm_3d(lo,hi,ng,dx,up,qp,fhp,dlo,dhi,blo,bhi)
          end if
       end do
       call destroy(bpt_hypterm)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

       !
       ! Calculate U'
       !
       call build(bpt_calcU, "calcU")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
       do n=1,nfabs(U)
          fhp => dataptr(Fhyp,  n)
          fdp => dataptr(Fdif,  n)
          upp => dataptr(Uprime,n)

          lo = lwb(get_box(U,n))
          hi = upb(get_box(U,n))

          do m = 1, ncons
             !$OMP PARALLEL DO PRIVATE(i,j,k)
             do k = lo(3),hi(3)
                do j = lo(2),hi(2)
                   do i = lo(1),hi(1)
                      upp(i,j,k,m) = upp(i,j,k,m) + fhp(i,j,k,m) + fdp(i,j,k,m)
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          end do
       end do
       call destroy(bpt_calcU)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    end if

    !
    ! NSCBC boundary
    !
    if (inc_ad) then
       ! XXX: MWE: not sure if this is reasonable...

       call build(bpt_nscbc, "nscbc")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
       call nscbc(Q, U, Fdif, Uprime, dx)
       call destroy(bpt_nscbc)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer
    end if


    !
    ! Destroy
    !
    call destroy(Q)

    if (inc_ad) then
       call destroy(Fhyp)
       call destroy(Fdif)

       call destroy(mu)
       call destroy(xi)
       call destroy(lam)
       call destroy(Ddiag)
    end if

    if (update_courno) then
       call build(bpt_courno, "courno")
       call parallel_reduce(courno, courno_proc, MPI_MAX)
       call destroy(bpt_courno)
    end if

  end subroutine dUdt_compact

  subroutine compute_courno(Q, dx, courno)
    type(multifab), intent(in) :: Q
    double precision, intent(in) :: dx(Q%dim)
    double precision, intent(inout) :: courno

    integer :: n, ng, dm, lo(Q%dim), hi(Q%dim)
    double precision, pointer :: qp(:,:,:,:)

    dm = Q%dim
    ng = nghost(Q)

    do n=1,nfabs(Q)
       qp => dataptr(Q,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       if (dm .ne. 3) then
          call bl_error("Only 3D compute_courno is supported")
       else
          call comp_courno_3d(lo,hi,ng,dx,qp,courno)
       end if
    end do
  end subroutine compute_courno


  !
  ! Compute dU/dt given U using the compact stencil.
  !
  ! The Courant number (courno) might also be computed if passed.
  !
  subroutine dUdt_S3D (U, Uprime, dx, courno, istep)

    use probin_module, only : overlap_comm_comp, overlap_comm_gettrans, cfl_int

    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)
    double precision, intent(inout), optional :: courno
    integer,          intent(in   ), optional :: istep

    integer, parameter :: ng = 4

    type(multifab) :: mu, xi    ! viscosity
    type(multifab) :: lam       ! partial thermal conductivity
    type(multifab) :: Ddiag     ! diagonal components of rho * Y_k * D

    integer ::    lo(U%dim),    hi(U%dim)
    integer ::   dlo(U%dim),   dhi(U%dim)
    integer ::   blo(U%dim),   bhi(U%dim)
    integer :: i,j,k,m,n, dm
    integer :: ndq, ng_ctoprim, ng_gettrans

    logical :: update_courno
    double precision :: courno_proc

    type(layout)     :: la
    type(multifab)   :: Q, Fhyp, Fdif
    type(multifab)   :: qx, qy, qz
    type(mf_fb_data) :: U_fb_data, qx_fb_data, qy_fb_data, qz_fb_data

    double precision, pointer, dimension(:,:,:,:) :: up, fhp, fdp, qp, mup, xip, lamp, &
         Ddp, upp, qxp, qyp, qzp

    type(bl_prof_timer), save :: bpt_ctoprim, bpt_gettrans, bpt_hypterm, bpt_courno
    type(bl_prof_timer), save :: bpt_diffterm_1, bpt_diffterm_2, bpt_calcU, bpt_chemterm

    update_courno = .false.
    if (present(courno) .and. present(istep)) then
       if (mod(istep,cfl_int).eq.1 .or. cfl_int.le.1) then
          update_courno = .true.
       end if
    end if

    call multifab_fill_boundary_nowait(U, U_fb_data)

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    else
       call multifab_fill_boundary_finish(U, U_fb_data)
    end if

    call setval(Uprime, ZERO)

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    ndq = idX1+nspecies-1

    dm = U%dim
    la = get_layout(U)

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    call multifab_build(Q, la, nprim, ng)

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    call multifab_build(Fhyp, la, ncons, 0)
    call multifab_build(Fdif, la, ncons, 0)

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    call multifab_build(mu , la, 1, ng)
    call multifab_build(xi , la, 1, ng)
    call multifab_build(lam, la, 1, ng)
    call multifab_build(Ddiag, la, nspecies, ng)

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    ! these multifabs are used to store first-derivatives
    call multifab_build(qx, la, ndq, ng)
    call multifab_build(qy, la, ndq, ng)
    call multifab_build(qz, la, ndq, ng)

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
    call build(bpt_ctoprim, "ctoprim")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    call ctoprim(U, Q, ng_ctoprim)
    call destroy(bpt_ctoprim)            !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if

    if (update_courno) then
       call build(bpt_courno, "courno")
       courno_proc = -1.d50
       call compute_courno(Q, dx, courno_proc)
       call destroy(bpt_courno)
    end if

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(U, U_fb_data)
    end if
    
    !
    ! chemistry
    !
    call build(bpt_chemterm, "chemterm")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    do n=1,nfabs(Q)
       qp  => dataptr(Q,n)
       upp => dataptr(Uprime,n)
       
       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))
       
       if (dm .ne. 3) then
          call bl_error("Only 3D chemsitry_term is supported")
       else
          call chemterm_3d(lo,hi,ng,qp,upp)
       end if
    end do
    call destroy(bpt_chemterm)              !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

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
    call destroy(bpt_gettrans)             !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

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
          call destroy(bpt_gettrans)             !! ^^^^^^^^^^^^^^^^^^^^^^^ timer
       end if

       call multifab_fill_boundary_finish(U, U_fb_data)
    end if

    !
    ! S3D_diffterm1: first derivative terms
    !
    call build(bpt_diffterm_1, "diffterm_1")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    do n=1,nfabs(Q)
       qp  => dataptr(Q,n)
       fdp => dataptr(Fdif,n)

       mup  => dataptr(mu, n)
       xip  => dataptr(xi, n)

       qxp => dataptr(qx, n)
       qyp => dataptr(qy, n)
       qzp => dataptr(qz, n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       if (dm .ne. 3) then
          call bl_error("Only 3D S3D_diffterm is supported")
       else
          call S3D_diffterm_1(lo,hi,ng,ndq,dx,qp,fdp,mup,xip,qxp,qyp,qzp)
       end if
    end do
    call destroy(bpt_diffterm_1)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    qx_fb_data%tag = 1001
    call multifab_fill_boundary_nowait(qx, qx_fb_data, idim=1)
    qy_fb_data%tag = 1002
    call multifab_fill_boundary_nowait(qy, qy_fb_data, idim=2)
    qz_fb_data%tag = 1003
    call multifab_fill_boundary_nowait(qz, qz_fb_data, idim=3)

    if (overlap_comm_comp) then
       call multifab_fill_boundary_test(qx, qx_fb_data, idim=1)
       call multifab_fill_boundary_test(qy, qy_fb_data, idim=2)
       call multifab_fill_boundary_test(qz, qz_fb_data, idim=3)
    else
       call multifab_fill_boundary_finish(qx, qx_fb_data, idim=1)
       call multifab_fill_boundary_finish(qy, qy_fb_data, idim=2)
       call multifab_fill_boundary_finish(qz, qz_fb_data, idim=3)
    end if

    !
    ! Hyperbolic terms
    !
    call build(bpt_hypterm, "hypterm")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    do n=1,nfabs(Fhyp)
       up => dataptr(U,n)
       qp => dataptr(Q,n)
       fhp=> dataptr(Fhyp,n)
       
       lo = lwb(get_box(Fhyp,n))
       hi = upb(get_box(Fhyp,n))
       
       call get_data_lo_hi(n,dlo,dhi)
       call get_boxbc(n,blo,bhi)
       
       if (dm .ne. 3) then
          call bl_error("Only 3D hypterm is supported")
       else
          call hypterm_3d(lo,hi,ng,dx,up,qp,fhp,dlo,dhi,blo,bhi)
       end if
    end do
    call destroy(bpt_hypterm)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    call multifab_fill_boundary_finish(qx, qx_fb_data, idim=1)
    call multifab_fill_boundary_finish(qy, qy_fb_data, idim=2)
    call multifab_fill_boundary_finish(qz, qz_fb_data, idim=3)
 
    !
    ! S3D_diffterm2: d(a du/dx)/dx terms
    !
    call build(bpt_diffterm_2, "diffterm_2")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    do n=1,nfabs(Q)
       qp  => dataptr(Q,n)
       fdp => dataptr(Fdif,n)

       mup  => dataptr(mu   , n)
       xip  => dataptr(xi   , n)
       lamp => dataptr(lam  , n)
       Ddp  => dataptr(Ddiag, n)

       qxp => dataptr(qx, n)
       qyp => dataptr(qy, n)
       qzp => dataptr(qz, n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       if (dm .ne. 3) then
          call bl_error("Only 3D S3D_diffterm is supported")
       else
          call S3D_diffterm_2(lo,hi,ng,ndq,dx,qp,fdp,mup,xip,lamp,Ddp,qxp,qyp,qzp)
       end if
    end do
    call destroy(bpt_diffterm_2)            !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    !
    ! Calculate U'
    !
    call build(bpt_calcU, "calcU")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    do n=1,nfabs(U)
       fhp => dataptr(Fhyp,  n)
       fdp => dataptr(Fdif,  n)
       upp => dataptr(Uprime,n)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       do m = 1, ncons
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   upp(i,j,k,m) = upp(i,j,k,m) + fhp(i,j,k,m) + fdp(i,j,k,m)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do
    call destroy(bpt_calcU)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    call destroy(Q)

    call destroy(Fhyp)
    call destroy(Fdif)

    call destroy(mu)
    call destroy(xi)
    call destroy(lam)
    call destroy(Ddiag)

    call destroy(qx)
    call destroy(qy)
    call destroy(qz)

    if (update_courno) then
       call build(bpt_courno, "courno")
       call parallel_reduce(courno, courno_proc, MPI_MAX)
       call destroy(bpt_courno)
    end if

  end subroutine dUdt_S3D

end module advance_module
