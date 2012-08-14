module advance_module

  use bl_error_module
  use derivative_stencil_module
  use kernels_module
  use multifab_module
  use time_module
  use transport_properties
  use variables_module

  use chemistry_module, only : nspecies
  use probin_module, only : advance_method, cflfac, fixed_dt, init_shrink, max_dt, &
       max_dt_growth, small_dt, stop_time, verbose

  implicit none

  private
  public advance

contains

  subroutine advance(U, dt, dx, istep)

    type(multifab),    intent(inout) :: U
    double precision,  intent(  out) :: dt
    double precision,  intent(in   ) :: dx(U%dim) 
    integer, intent(in) :: istep

    if (advance_method == 2) then
       call bl_error("call advance_sdc")
    else
       call advance_rk3(U, dt, dx, istep)
    end if

  end subroutine advance


  subroutine advance_rk3 (U,dt,dx,istep)
    use derivative_stencil_module, only : stencil, compact, S3D

    type(multifab),    intent(inout) :: U
    double precision,  intent(  out) :: dt
    double precision,  intent(in   ) :: dx(U%dim)
    integer, intent(in) :: istep

    integer          :: ng
    double precision :: courno, courno_proc, dtold
    type(layout)     :: la
    type(multifab)   :: Uprime, Unew

    ng = nghost(U)
    la = get_layout(U)

    call multifab_build(Uprime, la, ncons, 0)
    call multifab_build(Unew,   la, ncons, ng)

    ! RK Step 1
    courno_proc = 1.0d-50

    if (stencil .eq. compact) then
       call dUdt(U, Uprime, dx, courno=courno_proc)
    else if (stencil .eq. S3D) then
       call dUdt_S3D(U, Uprime, dx, courno=courno_proc)
    else
       call bl_error("advance_rk3: unknown stencil type")
    end if

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
       dt = cflfac / courno

       if (parallel_IOProcessor()) then
          print*, "CFL gives dt =", dt
       end if

       if (istep .eq. 1) then
          dt = dt * init_shrink
          if (parallel_IOProcessor()) then
             print*,'init_shrink factor limits the first dt =',dt
          end if
       else
          if (dt .gt. dtold * max_dt_growth) then
             dt = dtold * max_dt_growth
             if (parallel_IOProcessor()) then
                print*,'dt_growth factor limits the new dt =',dt
             end if
          end if
       end if

       if(dt .gt. max_dt) then
          if (parallel_IOProcessor()) then
             print*,'max_dt limits the new dt =',max_dt
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
                print*, "Stop time limits dt =",dt
             end if
          end if
       end if
       
       if (parallel_IOProcessor()) then
          print *, ""
       end if

    end if

    call update_rk3(Zero,Unew, One,U, dt,Uprime)
    call reset_density(Unew)

    ! RK Step 2
    if (stencil .eq. compact) then
       call dUdt(Unew,Uprime,dx)
    else if (stencil .eq. S3D) then
       call dUdt_S3D(Unew,Uprime,dx)
    else
       call bl_error("advance_rk3: unknown stencil type")
    end if

    call update_rk3(OneQuarter,Unew, ThreeQuarters,U, OneQuarter*dt,Uprime)
    call reset_density(Unew)

    ! RK Step 3
    if (stencil .eq. compact) then
       call dUdt(Unew,Uprime,dx)
    else if (stencil .eq. S3D) then
       call dUdt_S3D(Unew,Uprime,dx)
    else
       call bl_error("advance_rk3: unknown stencil type")
    end if       

    call update_rk3(OneThird,U, TwoThirds,Unew, TwoThirds*dt,Uprime)
    call reset_density(U)

    call destroy(Unew)
    call destroy(Uprime)

  end subroutine advance_rk3


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute dU/dt given U.
  !
  ! The Courant number (courno) is also computed if passed.
  !
  subroutine dUdt (U, Uprime, dx, courno)

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

  end subroutine dUdt



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



  ! S3D style
  subroutine dUdt_S3D (U, Uprime, dx, courno)

    type(multifab),   intent(inout) :: U, Uprime
    double precision, intent(in   ) :: dx(U%dim)
    double precision, intent(inout), optional :: courno

    integer, parameter :: ng = 4

    type(multifab) :: mu, xi ! viscosity
    type(multifab) :: lam ! partial thermal conductivity
    type(multifab) :: Ddiag ! diagonal components of rho * Y_k * D

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
