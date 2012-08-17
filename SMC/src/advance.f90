module advance_module

  use bl_error_module
  use derivative_stencil_module
  use kernels_module
  use multifab_module
  use sdcquad_module
  use time_module
  use transport_properties
  use variables_module

  use chemistry_module, only : nspecies

  implicit none

  private
  public advance

contains

  subroutine advance(U, dt, dx, sdc, istep)

    use probin_module, only : advance_method

    type(multifab),    intent(inout) :: U
    double precision,  intent(  out) :: dt
    double precision,  intent(in   ) :: dx(U%dim) 
    integer,           intent(in   ) :: istep
    type(sdcquad),     intent(in   ) :: sdc

    if (advance_method == 2) then
       call advance_sdc(U,dt,dx,sdc,istep)
    else
       call advance_rk3(U,dt,dx,istep)
    end if

  end subroutine advance


  !
  ! Advance U using SSP RK3
  !
  subroutine advance_rk3 (U,dt,dx,istep)
    type(multifab),    intent(inout) :: U
    double precision,  intent(  out) :: dt
    double precision,  intent(in   ) :: dx(U%dim)
    integer, intent(in) :: istep

    integer          :: ng
    double precision :: courno_proc
    type(layout)     :: la
    type(multifab)   :: Uprime, Unew

    ng = nghost(U)
    la = get_layout(U)

    call multifab_build(Uprime, la, ncons, 0)
    call multifab_build(Unew,   la, ncons, ng)

    ! RK Step 1
    courno_proc = 1.0d-50

    call dUdt(U, Uprime, dx, courno=courno_proc)
    call set_dt(dt, courno_proc, istep)
    call update_rk3(Zero,Unew, One,U, dt,Uprime)
    call reset_density(Unew)

    ! RK Step 2
    call dUdt(Unew, Uprime, dx)

    call update_rk3(OneQuarter, Unew, ThreeQuarters, U, OneQuarter*dt, Uprime)
    call reset_density(Unew)

    ! RK Step 3
    call dUdt(Unew, Uprime, dx)
    call update_rk3(OneThird, U, TwoThirds, Unew, TwoThirds*dt, Uprime)
    call reset_density(U)

    call destroy(Unew)
    call destroy(Uprime)

  end subroutine advance_rk3

  !
  ! Advance U using SDC time-stepping
  !
  subroutine advance_sdc(U, dt, dx, sdc, istep)

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt
    double precision,  intent(in   ) :: dx(U%dim)
    type(sdcquad),     intent(in   ) :: sdc
    integer,           intent(in   ) :: istep

    integer          :: k, m, ng
    double precision :: courno_proc, res_proc, res
    type(layout)     :: la
    type(multifab)   :: uSDC(sdc%nnodes), fSDC(sdc%nnodes), S(sdc%nnodes-1)

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
    courno_proc = 1.0d-50

    call copy(uSDC(1), U)
    call dUdt(uSDC(1), fSDC(1), dx, courno_proc)

    do m = 2, sdc%nnodes
       call copy(uSDC(m), uSDC(1))
       call copy(fSDC(m), fSDC(1))
    end do

    call set_dt(dt, courno_proc, istep)

    ! perform sdc iterations
    res = 0.0d0

    do k = 1, sdc%iters
       call sdc_sweep(uSDC, fSDC, S, dx, dt, sdc)

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
  ! Perform one SDC sweep.
  !
  subroutine sdc_sweep (uSDC,fSDC,S,dx,dt,sdc)
    type(sdcquad),     intent(in   ) :: sdc
    type(multifab),    intent(inout) :: uSDC(sdc%nnodes), fSDC(sdc%nnodes), S(sdc%nnodes-1)
    double precision,  intent(in   ) :: dt, dx(uSDC(1)%dim)

    integer :: m, n
    double precision :: dtsdc(sdc%nnodes-1)

    ! compute integrals (compact forward euler)
    do m = 1, sdc%nnodes-1
       call setval(S(m), 0.0d0)
       do n = 1, sdc%nnodes
          call saxpy(S(m), sdc%smats(m,n,1), fSDC(n))
       end do
    end do

    ! perform sub-step correction
    dtsdc = dt * (sdc%nodes(2:sdc%nnodes) - sdc%nodes(1:sdc%nnodes-1))
    do m = 1, sdc%nnodes-1

       ! U(m+1) = U(m) + dt dUdt(m) + dt S(m)

       call copy(uSDC(m+1), uSDC(m))
       call saxpy(uSDC(m+1), dtsdc(m), fSDC(m))
       call saxpy(uSDC(m+1), dt, S(m))

       call dUdt(uSDC(m+1), fSDC(m+1), dx)

    end do

  end subroutine sdc_sweep


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute SDC residual.
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
  ! Compute new time-step size
  !
  subroutine set_dt(dt, courno_proc, istep)

    use probin_module, only : cflfac, fixed_dt, init_shrink, max_dt, max_dt_growth, small_dt, stop_time

    double precision, intent(inout) :: dt
    double precision, intent(in   ) :: courno_proc
    integer,          intent(in   ) :: istep

    double precision :: dtold, courno

    if (fixed_dt > 0.d0) then

       dt = fixed_dt

       if (parallel_IOProcessor()) then
          print*, ""
          print*, "Setting fixed dt =",dt
          print*, ""
       end if

    else

       call parallel_reduce(courno, courno_proc, MPI_MAX)

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

    do n=1,nboxes(U1)
       if ( remote(U1,n) ) cycle

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
  ! The Courant number (courno) is also computed if passed.
  !
  subroutine dUdt (U, Uprime, dx, courno)
    use derivative_stencil_module, only : stencil, compact, s3d

    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)
    double precision, intent(inout), optional :: courno

    if (stencil .eq. compact) then
       call dUdt_compact(U, Uprime, dx, courno)
    else if (stencil .eq. s3d) then
       call dUdt_S3D(U, Uprime, dx, courno)
    else
       call bl_error("advance: unknown stencil type")
    end if       

  end subroutine dUdt


  !
  ! Compute dU/dt given U using the compact stencil.
  !
  ! The Courant number (courno) is also computed if passed.
  !
  subroutine dUdt_compact (U, Uprime, dx, courno)

    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)
    double precision, intent(inout), optional :: courno

    type(multifab) :: mu, xi ! viscosity
    type(multifab) :: lam ! partial thermal conductivity
    type(multifab) :: Ddiag ! diagonal components of rho * Y_k * D

    integer          :: lo(U%dim), hi(U%dim), i,j,k,m,n, ng, dm
    type(layout)     :: la
    type(multifab)   :: Q, Fhyp, Fdif

    double precision, pointer, dimension(:,:,:,:) :: up, fhp, fdp, qp, mup, xip, lamp, Ddp, upp

    dm = U%dim
    ng = nghost(U)
    la = get_layout(U)

    call multifab_fill_boundary(U)

    call multifab_build(Q, la, nprim, ng)

    call multifab_build(Fhyp, la, ncons, 0)
    call multifab_build(Fdif, la, ncons, 0)

    call multifab_build(mu , la, 1, ng)
    call multifab_build(xi , la, 1, ng)
    call multifab_build(lam, la, 1, ng)
    call multifab_build(Ddiag, la, nspecies, ng)

    !
    ! Calculate primitive variables based on U
    !
    call ctoprim(U, Q, ng)

    if (present(courno)) then
       call compute_courno(Q, dx, courno)
    end if

    call get_transport_properties(Q, mu, xi, lam, Ddiag)

    !
    ! Hyperbolic terms
    !
    do n=1,nboxes(Fhyp)
       if ( remote(Fhyp,n) ) cycle

       up => dataptr(U,n)
       qp => dataptr(Q,n)
       fhp=> dataptr(Fhyp,n)

       lo = lwb(get_box(Fhyp,n))
       hi = upb(get_box(Fhyp,n))

       if (dm .ne. 3) then
          call bl_error("Only 3D hypterm is supported")
       else
          call hypterm_3d(lo,hi,ng,dx,up,qp,fhp)
       end if
    end do

    
    !
    ! Transport terms
    !
    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

       qp  => dataptr(Q,n)
       fdp => dataptr(Fdif,n)

       mup  => dataptr(mu   , n)
       xip  => dataptr(xi   , n)
       lamp => dataptr(lam  , n)
       Ddp  => dataptr(Ddiag, n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       if (dm .ne. 3) then
          call bl_error("Only 3D compact_diffterm is supported")
       else
          call compact_diffterm_3d(lo,hi,ng,dx,qp,fdp,mup,xip,lamp,Ddp)
       end if
    end do

    !
    ! Calculate U'
    !
    do n=1,nboxes(U)
       if ( remote(U,n) ) cycle
       
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
                   upp(i,j,k,m) = fhp(i,j,k,m) + fdp(i,j,k,m)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do

    ! 
    ! Add chemistry
    !
    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

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

    call destroy(Q)

    call destroy(Fhyp)
    call destroy(Fdif)

    call destroy(mu)
    call destroy(xi)
    call destroy(lam)
    call destroy(Ddiag)

  end subroutine dUdt_compact

  subroutine compute_courno(Q, dx, courno)
    type(multifab), intent(in) :: Q
    double precision, intent(in) :: dx(Q%dim)
    double precision, intent(inout) :: courno

    integer :: n, ng, dm, lo(Q%dim), hi(Q%dim)
    double precision, pointer :: qp(:,:,:,:)

    dm = Q%dim
    ng = nghost(Q)

    do n=1,nboxes(Q)
       if (remote(Q,n)) cycle

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
  ! The Courant number (courno) is also computed if passed.
  !
  subroutine dUdt_S3D (U, Uprime, dx, courno)

    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)
    double precision, intent(inout), optional :: courno

    integer, parameter :: ng = 4

    type(multifab) :: mu, xi    ! viscosity
    type(multifab) :: lam       ! partial thermal conductivity
    type(multifab) :: Ddiag     ! diagonal components of rho * Y_k * D

    integer          :: lo(U%dim), hi(U%dim), i,j,k,m,n, dm
    type(layout)     :: la
    type(multifab)   :: Q, Fhyp, Fdif
    type(multifab)   :: qx, qy, qz

    double precision, pointer, dimension(:,:,:,:) :: up, fhp, fdp, qp, mup, xip, lamp, &
         Ddp, upp, qxp, qyp, qzp

    integer :: ndq

    ndq = idX1+nspecies-1

    dm = U%dim
    la = get_layout(U)

    call multifab_build(Q, la, nprim, ng)

    call multifab_build(Fhyp, la, ncons, 0)
    call multifab_build(Fdif, la, ncons, 0)

    call multifab_build(mu , la, 1, ng)
    call multifab_build(xi , la, 1, ng)
    call multifab_build(lam, la, 1, ng)
    call multifab_build(Ddiag, la, nspecies, ng)

    ! these multifabs are used to store first-derivatives
    call multifab_build(qx, la, ndq, ng)
    call multifab_build(qy, la, ndq, ng)
    call multifab_build(qz, la, ndq, ng)

    !
    ! Calculate primitive variables based on U
    !
    call ctoprim(U, Q, 0)

    call multifab_fill_boundary(Q)
    call multifab_fill_boundary(U)

    if (present(courno)) then
       call compute_courno(Q, dx, courno)
    end if

    call get_transport_properties(Q, mu, xi, lam, Ddiag)

    !
    ! Hyperbolic terms
    !
    do n=1,nboxes(Fhyp)
       if ( remote(Fhyp,n) ) cycle

       up => dataptr(U,n)
       qp => dataptr(Q,n)
       fhp=> dataptr(Fhyp,n)

       lo = lwb(get_box(Fhyp,n))
       hi = upb(get_box(Fhyp,n))

       if (dm .ne. 3) then
          call bl_error("Only 3D hypterm is supported")
       else
          call hypterm_3d(lo,hi,ng,dx,up,qp,fhp)
       end if
    end do

    !
    ! Transport terms
    ! S3D_diffterm1: first derivative terms
    ! S3D_diffterm2: d(a du/dx)/dx terms
    !
    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

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

    call multifab_fill_boundary(qx)
    call multifab_fill_boundary(qy)
    call multifab_fill_boundary(qz)

    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

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

    !
    ! Calculate U'
    !
    do n=1,nboxes(U)
       if ( remote(U,n) ) cycle
       
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
                   upp(i,j,k,m) = fhp(i,j,k,m) + fdp(i,j,k,m)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do

    ! 
    ! Add chemistry
    !
    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

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

  end subroutine dUdt_S3D

end module advance_module
