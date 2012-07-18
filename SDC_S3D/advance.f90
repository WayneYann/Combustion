 module advance_module

  use bl_error_module
  use multifab_module
  use bl_prof_module
  use sdcquad_module

  implicit none

  private

  !
  ! These index constants are shared with the initial data routine.
  !
  integer, parameter, public :: irho = 1
  integer, parameter, public :: imx  = 2
  integer, parameter, public :: imy  = 3
  integer, parameter, public :: imz  = 4
  integer, parameter, public :: iene = 5

  integer, parameter :: qu    = 2
  integer, parameter :: qv    = 3
  integer, parameter :: qw    = 4
  integer, parameter :: qpres = 5

  double precision, parameter :: ALP =  0.8d0
  double precision, parameter :: BET = -0.2d0
  double precision, parameter :: GAM =  4.d0/105.d0
  double precision, parameter :: DEL = -1.d0/280.d0

  type :: feval_ctx_t

     type(layout) :: la
     integer      :: nc, ng

     double precision :: dx(3), eta, alam

  end type feval_ctx_t

  public :: advance, dUdt, feval_ctx_t

contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Advance one time-step.
  !
  ! The time step is adjusted based on the CFL number and the computed
  ! Courant number.
  !
  subroutine advance (U,dt,dx,cfl,time,tfinal,method,ctx,sdc)

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt
    double precision,  intent(in   ) :: dx(U%dim), cfl, time, tfinal
    type(feval_ctx_t), intent(in   ) :: ctx 
    integer,           intent(in   ) :: method 
    type(sdcquad),     intent(in   ) :: sdc

    type(bl_prof_timer), save :: bpt_advance

    call build(bpt_advance, "bpt_advance")

    if (method == 2) then
       call advance_sdc(U,dt,dx,cfl,time,tfinal,ctx,sdc)
    else
       call advance_rk3(U,dt,dx,cfl,time,tfinal,ctx)
    end if

    call destroy(bpt_advance)

  end subroutine advance


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Advance one time-step using SDC.
  !
  ! The time step is adjusted based on the CFL number and the computed
  ! Courant number.
  !
  subroutine advance_sdc (U,dt,dx,cfl,time,tfinal,ctx,sdc)

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt
    double precision,  intent(in   ) :: dx(U%dim), cfl, time, tfinal
    type(sdcquad),     intent(in   ) :: sdc
    type(feval_ctx_t), intent(in   ) :: ctx

    integer          :: k, m, nc, ng
    double precision :: courno, courno_proc, res
    type(layout)     :: la
    type(multifab)   :: uSDC(sdc%nnodes), fSDC(sdc%nnodes)

    nc = ncomp(U)
    ng = nghost(U)
    la = get_layout(U)

    !
    ! Set provisional solution, compute Courant number, set dt
    !
    do m = 1, sdc%nnodes
       call build(uSDC(m), la, nc, ng)
       call build(fSDC(m), la, nc, 0)
    end do

    courno_proc = 1.0d-50

    call sdc_copy(uSDC(1), U)
    call dUdt(uSDC(1),fSDC(1),ctx,courno=courno_proc)

    call parallel_reduce(courno, courno_proc, MPI_MAX)

    if (cfl > 0.0d0) then
       dt = cfl / courno
    end if

    if (time + dt > tfinal) then
       dt = tfinal - time
    end if

    if ( parallel_IOProcessor() ) then
       print*, "dt,courno", dt, courno
    end if

    do m = 2, sdc%nnodes
       call sdc_copy(uSDC(m), uSDC(1))
       call sdc_copy(fSDC(m), fSDC(1))
    end do

    !
    ! Perform SDC iterations
    !
    do k = 1, sdc%iters
       call sdc_sweep(uSDC, fSDC, dx, dt, ctx, sdc)

       if (sdc%tol_residual > 0.d0) then
          res = sdc_residual(uSDC, fSDC, dt, sdc)

          print *, 'SDC iteration', k, 'residual = ', res

          if (res < sdc%tol_residual) then
             print *, 'SDC RESIDUAL CONDITION MET'
             exit
          end if
       end if
    end do

    call sdc_copy(U, uSDC(sdc%nnodes))
    
    do m = 1, sdc%nnodes
       call destroy(uSDC(m))
       call destroy(fSDC(m))
    end do

  end subroutine advance_sdc




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Perform one SDC sweep.
  !
  subroutine sdc_sweep (uSDC,fSDC,dx,dt,ctx,sdc)
    type(sdcquad),     intent(in   ) :: sdc
    type(multifab),    intent(inout) :: uSDC(sdc%nnodes), fSDC(sdc%nnodes)
    double precision,  intent(in   ) :: dt, dx(uSDC(1)%dim)
    type(feval_ctx_t), intent(in   ) :: ctx

    integer        :: m, n, nc
    type(multifab) :: S(sdc%nnodes-1)
    type(layout)   :: la

    double precision :: dtsdc(sdc%nnodes-1)

    la = get_layout(uSDC(1))
    nc = ncomp(uSDC(1))

    !
    ! Compute integrals (compact forward Euler)
    ! 
    do m = 1, sdc%nnodes-1
       call build(S(m), la, nc, 0)
       call setval(S(m), 0.0d0)
       do n = 1, sdc%nnodes
          call saxpy(S(m), sdc%smats(m,n,1), fSDC(n))
       end do
    end do

    !
    ! Perform sub-step correction
    !
    dtsdc = dt * (sdc%nodes(2:sdc%nnodes) - sdc%nodes(1:sdc%nnodes-1))
    do m = 1, sdc%nnodes-1

       ! U(m+1) = U(m) + dt dUdt(m) + dt S(m)

       call sdc_copy(uSDC(m+1), uSDC(m))
       call saxpy(uSDC(m+1), dtsdc(m), fSDC(m))
       call saxpy(uSDC(m+1), dt, S(m))

       call dUdt(uSDC(m+1), fSDC(m+1), ctx)

    end do


    do m = 1, sdc%nnodes-1
       call destroy(S(m))
    end do

  end subroutine sdc_sweep


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute SDC residual.
  !
  function sdc_residual (uSDC,fSDC,dt,sdc) result(res)
    real(dp_t)                      :: res
    type(sdcquad),    intent(in   ) :: sdc
    type(multifab),   intent(inout) :: uSDC(sdc%nnodes), fSDC(sdc%nnodes)
    real(dp_t),       intent(in   ) :: dt      

    integer        :: m, n, nc
    type(multifab) :: R
    type(layout)   :: la

    la = get_layout(uSDC(1))
    nc = ncomp(uSDC(1))

    !
    ! Compute integral
    ! 
    call build(R, la, nc, 0)
    call sdc_copy(R, uSDC(1))

    do m = 1, sdc%nnodes-1
       do n = 1, sdc%nnodes
          call saxpy(R, dt*sdc%smat(m,n), fSDC(n))
       end do
    end do

    call saxpy(R, -1.0d0, uSDC(sdc%nnodes))
    
    res = norm_inf(R)

    call destroy(R)

  end function sdc_residual


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Copy
  !
  subroutine sdc_copy (dst, src)

    type(multifab),   intent(in   ) :: src
    type(multifab),   intent(inout) :: dst

    integer :: lo(src%dim), hi(src%dim), i, j, k, m, n, nc

    double precision, pointer, dimension(:,:,:,:) :: sp, dp

    nc = ncomp(src)

    do n=1,nboxes(src)
       if ( remote(src,n) ) cycle

       sp => dataptr(src,n)
       dp => dataptr(dst,n)

       lo = lwb(get_box(src,n))
       hi = upb(get_box(src,n))

       do m = 1, nc
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   dp(i,j,k,m) = sp(i,j,k,m)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do

  end subroutine sdc_copy


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Advance one time-step using RK3.
  !
  ! The time step is adjusted based on the CFL number and the computed
  ! Courant number.
  !
  subroutine advance_rk3 (U,dt,dx,cfl,time,tfinal,ctx)

    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt
    double precision,  intent(in   ) :: dx(U%dim), cfl, time, tfinal
    type(feval_ctx_t), intent(in   ) :: ctx

    integer          :: nc, ng
    double precision :: courno, courno_proc
    type(layout)     :: la
    type(multifab)   :: Uprime, Unew

    !
    ! Some arithmetic constants.
    !
    double precision, parameter :: Zero          = 0.d0
    double precision, parameter :: One           = 1.d0
    double precision, parameter :: OneThird      = 1.d0/3.d0
    double precision, parameter :: TwoThirds     = 2.d0/3.d0
    double precision, parameter :: OneQuarter    = 1.d0/4.d0
    double precision, parameter :: ThreeQuarters = 3.d0/4.d0


    nc = ncomp(U)
    ng = nghost(U)
    la = get_layout(U)

    call multifab_build(Uprime, la, nc, 0)
    call multifab_build(Unew,   la, nc, ng)

    !
    ! Calculate U at time N+1/3
    !
    ! Also calculate Courant number and set dt
    !
    courno_proc = 1.0d-50

    call dUdt(U,Uprime,ctx,courno=courno_proc)

    call parallel_reduce(courno, courno_proc, MPI_MAX)

    if (cfl > 0.0d0) then
       dt = cfl / courno
    end if

    if (time + dt > tfinal) then
       dt = tfinal - time
    end if

    if ( parallel_IOProcessor() ) then
       print*, "dt,courno", dt, courno
    end if

    call update_rk3(One,U,Zero,Unew,dt,Uprime,Unew)

    !
    ! Calculate U at time N+2/3
    !
    call dUdt(Unew,Uprime,ctx)
    call update_rk3(ThreeQuarters,U,OneQuarter,Unew,OneQuarter*dt,Uprime,Unew)

    !
    ! Calculate U at time N+1
    !
    call dUdt(Unew,Uprime,ctx)
    call update_rk3(OneThird,U,TwoThirds,Unew,TwoThirds*dt,Uprime,U)

    call destroy(Unew)
    call destroy(Uprime)

  end subroutine advance_rk3


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute Unew = a U1 + b U2 + c Uprime.
  !
  subroutine update_rk3 (a,U1,b,U2,c,Uprime,Unew)

    type(multifab),   intent(in   ) :: U1, U2, Uprime
    type(multifab),   intent(inout) :: Unew
    double precision, intent(in   ) :: a, b, c

    integer :: lo(U1%dim), hi(U1%dim), i, j, k, m, n, nc

    double precision, pointer, dimension(:,:,:,:) :: u1p, u2p, upp, unp

    nc = ncomp(U1)

    do n=1,nboxes(U1)
       if ( remote(U1,n) ) cycle

       u1p => dataptr(U1,    n)
       u2p => dataptr(U2,    n)
       unp => dataptr(Unew,  n)
       upp => dataptr(Uprime,n)

       lo = lwb(get_box(Unew,n))
       hi = upb(get_box(Unew,n))

       do m = 1, nc
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   unp(i,j,k,m) = a * u1p(i,j,k,m) + b * u2p(i,j,k,m) + c * upp(i,j,k,m)
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
  subroutine dUdt (U,Uprime,ctx,courno)

    type(multifab),    intent(inout) :: U, Uprime
    type(feval_ctx_t), intent(in   ) :: ctx
    double precision,  intent(inout), optional :: courno

    integer          :: lo(U%dim), hi(U%dim), i, j, k, m, n, nc, ng
    type(layout)     :: la
    type(multifab)   :: D, F, Q

    double precision, pointer, dimension(:,:,:,:) :: up, dp, fp, qp, upp

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

    !
    ! Calculate primitive variables based on U
    !
    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

       up => dataptr(U,n)
       qp => dataptr(Q,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call ctoprim(lo,hi,up,qp,ctx%dx,ng,courno=courno)
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

       call diffterm(lo,hi,ng,ctx%dx,qp,dp,ctx%eta,ctx%alam)
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

       call hypterm(lo,hi,ng,ctx%dx,up,qp,fp)
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

  end subroutine dUdt


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute primitive variables based on conservative variables.
  !
  ! The Courant number (courno) is also computed if passed.
  !
  subroutine ctoprim (lo,hi,u,q,dx,ng,courno)

    use bl_prof_module

    integer,          intent(in ) :: lo(3), hi(3), ng
    double precision, intent(in ) :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,5)
    double precision, intent(in ) :: dx(3)
    double precision, intent(out) :: q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,6)

    double precision, intent(inout), optional :: courno

    integer          :: i, j, k
    double precision :: c, eint, courx, coury, courz, courmx, courmy, courmz, rhoinv
    double precision :: dx1inv, dx2inv, dx3inv, CVinv

    double precision, parameter :: GAMMA = 1.4d0
    double precision, parameter :: CV    = 8.3333333333d6

    type(bl_prof_timer), save :: bpt_ctoprim, bpt_ctoprim_loop1, bpt_ctoprim_loop2

    call build(bpt_ctoprim, "bpt_ctoprim")

    CVinv = 1.0d0 / CV

    call build(bpt_ctoprim_loop1, "bpt_ctoprim_loop1")
    !$OMP PARALLEL DO PRIVATE(i,j,k,eint,rhoinv)
    do k = lo(3)-ng,hi(3)+ng
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng

             rhoinv     = 1.0d0/u(i,j,k,1)
             q(i,j,k,1) = u(i,j,k,1)
             q(i,j,k,2) = u(i,j,k,2)*rhoinv
             q(i,j,k,3) = u(i,j,k,3)*rhoinv
             q(i,j,k,4) = u(i,j,k,4)*rhoinv

             eint = u(i,j,k,5)*rhoinv - 0.5d0*(q(i,j,k,2)**2 + q(i,j,k,3)**2 + q(i,j,k,4)**2)

             q(i,j,k,5) = (GAMMA-1.d0)*eint*u(i,j,k,1)
             q(i,j,k,6) = eint * CVinv

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_ctoprim_loop1)

    if ( present(courno) ) then

       courmx = -Huge(courmx)
       courmy = -Huge(courmy)
       courmz = -Huge(courmz)

       dx1inv = 1.0d0 / dx(1)
       dx2inv = 1.0d0 / dx(2)
       dx3inv = 1.0d0 / dx(3)

       call build(bpt_ctoprim_loop2, "bpt_ctoprim_loop2")
       !$OMP PARALLEL DO PRIVATE(i,j,k,c,courx,coury,courz) REDUCTION(max:courmx,courmy,courmz)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                c     = sqrt(GAMMA*q(i,j,k,5)/q(i,j,k,1))
                courx = ( c+abs(q(i,j,k,2)) ) * dx1inv
                coury = ( c+abs(q(i,j,k,3)) ) * dx2inv
                courz = ( c+abs(q(i,j,k,4)) ) * dx3inv

                courmx = max( courmx, courx )
                courmy = max( courmy, coury )
                courmz = max( courmz, courz )

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
       call destroy(bpt_ctoprim_loop2)
       !
       ! Compute running max of Courant number over grids.
       !
       courno = max( courmx, courmy, courmz , courno )

    end if

    call destroy(bpt_ctoprim)

  end subroutine ctoprim


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute hyperbolic part (due to boundary fluxes) of dU/dt.
  !
  subroutine hypterm (lo,hi,ng,dx,cons,q,flux)

    use bl_prof_module

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,5)
    double precision, intent(in ) ::    q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,6)
    double precision, intent(out) :: flux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)

    integer          :: i,j,k
    double precision :: unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4
    double precision :: dxinv(3)

    type(bl_prof_timer), save :: bpt_hypterm, bpt_hypterm_loop1, bpt_hypterm_loop2, bpt_hypterm_loop3

    call build(bpt_hypterm, "bpt_hypterm")


    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do


    call build(bpt_hypterm_loop1, "bpt_hypterm_loop1")
    !$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = q(i+1,j,k,qu)
             unp2 = q(i+2,j,k,qu)
             unp3 = q(i+3,j,k,qu)
             unp4 = q(i+4,j,k,qu)

             unm1 = q(i-1,j,k,qu)
             unm2 = q(i-2,j,k,qu)
             unm3 = q(i-3,j,k,qu)
             unm4 = q(i-4,j,k,qu)

             flux(i,j,k,irho)= - &
                   (ALP*(cons(i+1,j,k,imx)-cons(i-1,j,k,imx)) &
                  + BET*(cons(i+2,j,k,imx)-cons(i-2,j,k,imx)) &
                  + GAM*(cons(i+3,j,k,imx)-cons(i-3,j,k,imx)) &
                  + DEL*(cons(i+4,j,k,imx)-cons(i-4,j,k,imx)))*dxinv(1)

             flux(i,j,k,imx)= - &
                   (ALP*(cons(i+1,j,k,imx)*unp1-cons(i-1,j,k,imx)*unm1 &
                  + (q(i+1,j,k,qpres)-q(i-1,j,k,qpres)))               &
                  + BET*(cons(i+2,j,k,imx)*unp2-cons(i-2,j,k,imx)*unm2 &
                  + (q(i+2,j,k,qpres)-q(i-2,j,k,qpres)))               &
                  + GAM*(cons(i+3,j,k,imx)*unp3-cons(i-3,j,k,imx)*unm3 &
                  + (q(i+3,j,k,qpres)-q(i-3,j,k,qpres)))               &
                  + DEL*(cons(i+4,j,k,imx)*unp4-cons(i-4,j,k,imx)*unm4 &
                  + (q(i+4,j,k,qpres)-q(i-4,j,k,qpres))))*dxinv(1)

             flux(i,j,k,imy)= - &
                   (ALP*(cons(i+1,j,k,imy)*unp1-cons(i-1,j,k,imy)*unm1) &
                  + BET*(cons(i+2,j,k,imy)*unp2-cons(i-2,j,k,imy)*unm2) &
                  + GAM*(cons(i+3,j,k,imy)*unp3-cons(i-3,j,k,imy)*unm3) &
                  + DEL*(cons(i+4,j,k,imy)*unp4-cons(i-4,j,k,imy)*unm4))*dxinv(1)

             flux(i,j,k,imz)= - &
                   (ALP*(cons(i+1,j,k,imz)*unp1-cons(i-1,j,k,imz)*unm1) &
                  + BET*(cons(i+2,j,k,imz)*unp2-cons(i-2,j,k,imz)*unm2) &
                  + GAM*(cons(i+3,j,k,imz)*unp3-cons(i-3,j,k,imz)*unm3) &
                  + DEL*(cons(i+4,j,k,imz)*unp4-cons(i-4,j,k,imz)*unm4))*dxinv(1)

             flux(i,j,k,iene)= - &
                   (ALP*(cons(i+1,j,k,iene)*unp1-cons(i-1,j,k,iene)*unm1 &
                  + (q(i+1,j,k,qpres)*unp1-q(i-1,j,k,qpres)*unm1))       &
                  + BET*(cons(i+2,j,k,iene)*unp2-cons(i-2,j,k,iene)*unm2 &
                  + (q(i+2,j,k,qpres)*unp2-q(i-2,j,k,qpres)*unm2))       &
                  + GAM*(cons(i+3,j,k,iene)*unp3-cons(i-3,j,k,iene)*unm3 &
                  + (q(i+3,j,k,qpres)*unp3-q(i-3,j,k,qpres)*unm3))       &
                  + DEL*(cons(i+4,j,k,iene)*unp4-cons(i-4,j,k,iene)*unm4 &
                  + (q(i+4,j,k,qpres)*unp4-q(i-4,j,k,qpres)*unm4)))*dxinv(1) 
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_hypterm_loop1)

    call build(bpt_hypterm_loop2, "bpt_hypterm_loop2")
    !$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = q(i,j+1,k,qv)
             unp2 = q(i,j+2,k,qv)
             unp3 = q(i,j+3,k,qv)
             unp4 = q(i,j+4,k,qv)

             unm1 = q(i,j-1,k,qv)
             unm2 = q(i,j-2,k,qv)
             unm3 = q(i,j-3,k,qv)
             unm4 = q(i,j-4,k,qv)

             flux(i,j,k,irho)=flux(i,j,k,irho) - &
                   (ALP*(cons(i,j+1,k,imy)-cons(i,j-1,k,imy)) &
                  + BET*(cons(i,j+2,k,imy)-cons(i,j-2,k,imy)) &
                  + GAM*(cons(i,j+3,k,imy)-cons(i,j-3,k,imy)) &
                  + DEL*(cons(i,j+4,k,imy)-cons(i,j-4,k,imy)))*dxinv(2)

             flux(i,j,k,imx)=flux(i,j,k,imx) - &
                   (ALP*(cons(i,j+1,k,imx)*unp1-cons(i,j-1,k,imx)*unm1) &
                  + BET*(cons(i,j+2,k,imx)*unp2-cons(i,j-2,k,imx)*unm2) &
                  + GAM*(cons(i,j+3,k,imx)*unp3-cons(i,j-3,k,imx)*unm3) &
                  + DEL*(cons(i,j+4,k,imx)*unp4-cons(i,j-4,k,imx)*unm4))*dxinv(2)

             flux(i,j,k,imy)=flux(i,j,k,imy) - &
                   (ALP*(cons(i,j+1,k,imy)*unp1-cons(i,j-1,k,imy)*unm1 &
                  + (q(i,j+1,k,qpres)-q(i,j-1,k,qpres)))               &
                  + BET*(cons(i,j+2,k,imy)*unp2-cons(i,j-2,k,imy)*unm2 &
                  + (q(i,j+2,k,qpres)-q(i,j-2,k,qpres)))               &
                  + GAM*(cons(i,j+3,k,imy)*unp3-cons(i,j-3,k,imy)*unm3 &
                  + (q(i,j+3,k,qpres)-q(i,j-3,k,qpres)))               &
                  + DEL*(cons(i,j+4,k,imy)*unp4-cons(i,j-4,k,imy)*unm4 &
                  + (q(i,j+4,k,qpres)-q(i,j-4,k,qpres))))*dxinv(2)

             flux(i,j,k,imz)=flux(i,j,k,imz) - &
                   (ALP*(cons(i,j+1,k,imz)*unp1-cons(i,j-1,k,imz)*unm1) &
                  + BET*(cons(i,j+2,k,imz)*unp2-cons(i,j-2,k,imz)*unm2) &
                  + GAM*(cons(i,j+3,k,imz)*unp3-cons(i,j-3,k,imz)*unm3) &
                  + DEL*(cons(i,j+4,k,imz)*unp4-cons(i,j-4,k,imz)*unm4))*dxinv(2)

             flux(i,j,k,iene)=flux(i,j,k,iene) - &
                   (ALP*(cons(i,j+1,k,iene)*unp1-cons(i,j-1,k,iene)*unm1 &
                  + (q(i,j+1,k,qpres)*unp1-q(i,j-1,k,qpres)*unm1))       &
                  + BET*(cons(i,j+2,k,iene)*unp2-cons(i,j-2,k,iene)*unm2 &
                  + (q(i,j+2,k,qpres)*unp2-q(i,j-2,k,qpres)*unm2))       &
                  + GAM*(cons(i,j+3,k,iene)*unp3-cons(i,j-3,k,iene)*unm3 &
                  + (q(i,j+3,k,qpres)*unp3-q(i,j-3,k,qpres)*unm3))       &
                  + DEL*(cons(i,j+4,k,iene)*unp4-cons(i,j-4,k,iene)*unm4 &
                  + (q(i,j+4,k,qpres)*unp4-q(i,j-4,k,qpres)*unm4)))*dxinv(2)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_hypterm_loop2)

    call build(bpt_hypterm_loop3, "bpt_hypterm_loop3")
    !$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = q(i,j,k+1,qw)
             unp2 = q(i,j,k+2,qw)
             unp3 = q(i,j,k+3,qw)
             unp4 = q(i,j,k+4,qw)

             unm1 = q(i,j,k-1,qw)
             unm2 = q(i,j,k-2,qw)
             unm3 = q(i,j,k-3,qw)
             unm4 = q(i,j,k-4,qw)

             flux(i,j,k,irho)=flux(i,j,k,irho) - &
                   (ALP*(cons(i,j,k+1,imz)-cons(i,j,k-1,imz)) &
                  + BET*(cons(i,j,k+2,imz)-cons(i,j,k-2,imz)) &
                  + GAM*(cons(i,j,k+3,imz)-cons(i,j,k-3,imz)) &
                  + DEL*(cons(i,j,k+4,imz)-cons(i,j,k-4,imz)))*dxinv(3)

             flux(i,j,k,imx)=flux(i,j,k,imx) - &
                   (ALP*(cons(i,j,k+1,imx)*unp1-cons(i,j,k-1,imx)*unm1) &
                  + BET*(cons(i,j,k+2,imx)*unp2-cons(i,j,k-2,imx)*unm2) &
                  + GAM*(cons(i,j,k+3,imx)*unp3-cons(i,j,k-3,imx)*unm3) &
                  + DEL*(cons(i,j,k+4,imx)*unp4-cons(i,j,k-4,imx)*unm4))*dxinv(3)

             flux(i,j,k,imy)=flux(i,j,k,imy) - &
                   (ALP*(cons(i,j,k+1,imy)*unp1-cons(i,j,k-1,imy)*unm1) &
                  + BET*(cons(i,j,k+2,imy)*unp2-cons(i,j,k-2,imy)*unm2) &
                  + GAM*(cons(i,j,k+3,imy)*unp3-cons(i,j,k-3,imy)*unm3) &
                  + DEL*(cons(i,j,k+4,imy)*unp4-cons(i,j,k-4,imy)*unm4))*dxinv(3)

             flux(i,j,k,imz)=flux(i,j,k,imz) - &
                   (ALP*(cons(i,j,k+1,imz)*unp1-cons(i,j,k-1,imz)*unm1 &
                  + (q(i,j,k+1,qpres)-q(i,j,k-1,qpres)))               &
                  + BET*(cons(i,j,k+2,imz)*unp2-cons(i,j,k-2,imz)*unm2 &
                  + (q(i,j,k+2,qpres)-q(i,j,k-2,qpres)))               &
                  + GAM*(cons(i,j,k+3,imz)*unp3-cons(i,j,k-3,imz)*unm3 &
                  + (q(i,j,k+3,qpres)-q(i,j,k-3,qpres)))               &
                  + DEL*(cons(i,j,k+4,imz)*unp4-cons(i,j,k-4,imz)*unm4 &
                  + (q(i,j,k+4,qpres)-q(i,j,k-4,qpres))))*dxinv(3)

             flux(i,j,k,iene)=flux(i,j,k,iene) - &
                   (ALP*(cons(i,j,k+1,iene)*unp1-cons(i,j,k-1,iene)*unm1 &
                  + (q(i,j,k+1,qpres)*unp1-q(i,j,k-1,qpres)*unm1))       &
                  + BET*(cons(i,j,k+2,iene)*unp2-cons(i,j,k-2,iene)*unm2 &
                  + (q(i,j,k+2,qpres)*unp2-q(i,j,k-2,qpres)*unm2))       &
                  + GAM*(cons(i,j,k+3,iene)*unp3-cons(i,j,k-3,iene)*unm3 &
                  + (q(i,j,k+3,qpres)*unp3-q(i,j,k-3,qpres)*unm3))       &
                  + DEL*(cons(i,j,k+4,iene)*unp4-cons(i,j,k-4,iene)*unm4 &
                  + (q(i,j,k+4,qpres)*unp4-q(i,j,k-4,qpres)*unm4)))*dxinv(3)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_hypterm_loop3)

    call destroy(bpt_hypterm)

  end subroutine hypterm


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute diffusive part of dU/dt.
  !
  subroutine diffterm (lo,hi,ng,dx,q,difflux,eta,alam)

    use bl_prof_module

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,6)
    double precision, intent(out) :: difflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)
    double precision, intent(in ) :: eta, alam

    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz

    double precision :: dxinv(3)
    double precision :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
    double precision :: divu, uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz,txx,tyy,tzz
    double precision :: mechwork, uxy,uxz,vyz,wzx,wzy,vyx
    integer          :: i,j,k

    double precision, parameter :: OneThird   = 1.0d0/3.0d0
    double precision, parameter :: TwoThirds  = 2.0d0/3.0d0
    double precision, parameter :: FourThirds = 4.0d0/3.0d0

    double precision, parameter :: CENTER = -205.d0/72.d0
    double precision, parameter :: OFF1   =    8.d0/5.d0
    double precision, parameter :: OFF2   =   -0.2d0
    double precision, parameter :: OFF3   =    8.d0/315.d0
    double precision, parameter :: OFF4   =   -1.d0/560.d0

    type(bl_prof_timer), save :: bpt_diffterm, bpt_diffterm_loop123, bpt_diffterm_loop4
    type(bl_prof_timer), save :: bpt_diffterm_loop5, bpt_diffterm_loop6, bpt_diffterm_loop7


    call build(bpt_diffterm, "bpt_diffterm")

    allocate(ux(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(uy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(uz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    difflux(:,:,:,irho) = 0.0d0

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    call build(bpt_diffterm_loop123, "bpt_diffterm_loop123")
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1),hi(1)

             ux(i,j,k)= &
                   (ALP*(q(i+1,j,k,qu)-q(i-1,j,k,qu)) &
                  + BET*(q(i+2,j,k,qu)-q(i-2,j,k,qu)) &
                  + GAM*(q(i+3,j,k,qu)-q(i-3,j,k,qu)) &
                  + DEL*(q(i+4,j,k,qu)-q(i-4,j,k,qu)))*dxinv(1)

             vx(i,j,k)= &
                   (ALP*(q(i+1,j,k,qv)-q(i-1,j,k,qv)) &
                  + BET*(q(i+2,j,k,qv)-q(i-2,j,k,qv)) &
                  + GAM*(q(i+3,j,k,qv)-q(i-3,j,k,qv)) &
                  + DEL*(q(i+4,j,k,qv)-q(i-4,j,k,qv)))*dxinv(1)

             wx(i,j,k)= &
                   (ALP*(q(i+1,j,k,qw)-q(i-1,j,k,qw)) &
                  + BET*(q(i+2,j,k,qw)-q(i-2,j,k,qw)) &
                  + GAM*(q(i+3,j,k,qw)-q(i-3,j,k,qw)) &
                  + DEL*(q(i+4,j,k,qw)-q(i-4,j,k,qw)))*dxinv(1)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)   
          do i=lo(1)-ng,hi(1)+ng

             uy(i,j,k)= &
                   (ALP*(q(i,j+1,k,qu)-q(i,j-1,k,qu)) &
                  + BET*(q(i,j+2,k,qu)-q(i,j-2,k,qu)) &
                  + GAM*(q(i,j+3,k,qu)-q(i,j-3,k,qu)) &
                  + DEL*(q(i,j+4,k,qu)-q(i,j-4,k,qu)))*dxinv(2)

             vy(i,j,k)= &
                   (ALP*(q(i,j+1,k,qv)-q(i,j-1,k,qv)) &
                  + BET*(q(i,j+2,k,qv)-q(i,j-2,k,qv)) &
                  + GAM*(q(i,j+3,k,qv)-q(i,j-3,k,qv)) &
                  + DEL*(q(i,j+4,k,qv)-q(i,j-4,k,qv)))*dxinv(2)

             wy(i,j,k)= &
                   (ALP*(q(i,j+1,k,qw)-q(i,j-1,k,qw)) &
                  + BET*(q(i,j+2,k,qw)-q(i,j-2,k,qw)) &
                  + GAM*(q(i,j+3,k,qw)-q(i,j-3,k,qw)) &
                  + DEL*(q(i,j+4,k,qw)-q(i,j-4,k,qw)))*dxinv(2)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3),hi(3)
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng

             uz(i,j,k)= &
                   (ALP*(q(i,j,k+1,qu)-q(i,j,k-1,qu)) &
                  + BET*(q(i,j,k+2,qu)-q(i,j,k-2,qu)) &
                  + GAM*(q(i,j,k+3,qu)-q(i,j,k-3,qu)) &
                  + DEL*(q(i,j,k+4,qu)-q(i,j,k-4,qu)))*dxinv(3)

             vz(i,j,k)= &
                   (ALP*(q(i,j,k+1,qv)-q(i,j,k-1,qv)) &
                  + BET*(q(i,j,k+2,qv)-q(i,j,k-2,qv)) &
                  + GAM*(q(i,j,k+3,qv)-q(i,j,k-3,qv)) &
                  + DEL*(q(i,j,k+4,qv)-q(i,j,k-4,qv)))*dxinv(3)

             wz(i,j,k)= &
                   (ALP*(q(i,j,k+1,qw)-q(i,j,k-1,qw)) &
                  + BET*(q(i,j,k+2,qw)-q(i,j,k-2,qw)) &
                  + GAM*(q(i,j,k+3,qw)-q(i,j,k-3,qw)) &
                  + DEL*(q(i,j,k+4,qw)-q(i,j,k-4,qw)))*dxinv(3)
          enddo
       enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call destroy(bpt_diffterm_loop123)

    call build(bpt_diffterm_loop4, "bpt_diffterm_loop4")
    !$OMP PARALLEL DO PRIVATE(i,j,k,uxx,uyy,uzz,vyx,wzx)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             uxx = (CENTER*q(i,j,k,qu) &
                  + OFF1*(q(i+1,j,k,qu)+q(i-1,j,k,qu)) &
                  + OFF2*(q(i+2,j,k,qu)+q(i-2,j,k,qu)) &
                  + OFF3*(q(i+3,j,k,qu)+q(i-3,j,k,qu)) &
                  + OFF4*(q(i+4,j,k,qu)+q(i-4,j,k,qu)))*dxinv(1)**2

             uyy = (CENTER*q(i,j,k,qu) &
                  + OFF1*(q(i,j+1,k,qu)+q(i,j-1,k,qu)) &
                  + OFF2*(q(i,j+2,k,qu)+q(i,j-2,k,qu)) &
                  + OFF3*(q(i,j+3,k,qu)+q(i,j-3,k,qu)) &
                  + OFF4*(q(i,j+4,k,qu)+q(i,j-4,k,qu)))*dxinv(2)**2

             uzz = (CENTER*q(i,j,k,qu) &
                  + OFF1*(q(i,j,k+1,qu)+q(i,j,k-1,qu)) &
                  + OFF2*(q(i,j,k+2,qu)+q(i,j,k-2,qu)) &
                  + OFF3*(q(i,j,k+3,qu)+q(i,j,k-3,qu)) &
                  + OFF4*(q(i,j,k+4,qu)+q(i,j,k-4,qu)))*dxinv(3)**2

             vyx = (ALP*(vy(i+1,j,k)-vy(i-1,j,k)) &
                  + BET*(vy(i+2,j,k)-vy(i-2,j,k)) &
                  + GAM*(vy(i+3,j,k)-vy(i-3,j,k)) &
                  + DEL*(vy(i+4,j,k)-vy(i-4,j,k)))*dxinv(1)

             wzx = (ALP*(wz(i+1,j,k)-wz(i-1,j,k)) &
                  + BET*(wz(i+2,j,k)-wz(i-2,j,k)) &
                  + GAM*(wz(i+3,j,k)-wz(i-3,j,k)) &
                  + DEL*(wz(i+4,j,k)-wz(i-4,j,k)))*dxinv(1)

             difflux(i,j,k,imx) = eta*(FourThirds*uxx + uyy + uzz + OneThird*(vyx+wzx))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_diffterm_loop4)

    call build(bpt_diffterm_loop5, "bpt_diffterm_loop5")
    !$OMP PARALLEL DO PRIVATE(i,j,k,vxx,vyy,vzz,uxy,wzy)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             vxx = (CENTER*q(i,j,k,qv) &
                  + OFF1*(q(i+1,j,k,qv)+q(i-1,j,k,qv)) &
                  + OFF2*(q(i+2,j,k,qv)+q(i-2,j,k,qv)) &
                  + OFF3*(q(i+3,j,k,qv)+q(i-3,j,k,qv)) &
                  + OFF4*(q(i+4,j,k,qv)+q(i-4,j,k,qv)))*dxinv(1)**2

             vyy = (CENTER*q(i,j,k,qv) &
                  + OFF1*(q(i,j+1,k,qv)+q(i,j-1,k,qv)) &
                  + OFF2*(q(i,j+2,k,qv)+q(i,j-2,k,qv)) &
                  + OFF3*(q(i,j+3,k,qv)+q(i,j-3,k,qv)) &
                  + OFF4*(q(i,j+4,k,qv)+q(i,j-4,k,qv)))*dxinv(2)**2

             vzz = (CENTER*q(i,j,k,qv) &
                  + OFF1*(q(i,j,k+1,qv)+q(i,j,k-1,qv)) &
                  + OFF2*(q(i,j,k+2,qv)+q(i,j,k-2,qv)) &
                  + OFF3*(q(i,j,k+3,qv)+q(i,j,k-3,qv)) &
                  + OFF4*(q(i,j,k+4,qv)+q(i,j,k-4,qv)))*dxinv(3)**2

             uxy = (ALP*(ux(i,j+1,k)-ux(i,j-1,k)) &
                  + BET*(ux(i,j+2,k)-ux(i,j-2,k)) &
                  + GAM*(ux(i,j+3,k)-ux(i,j-3,k)) &
                  + DEL*(ux(i,j+4,k)-ux(i,j-4,k)))*dxinv(2)

             wzy = (ALP*(wz(i,j+1,k)-wz(i,j-1,k)) &
                  + BET*(wz(i,j+2,k)-wz(i,j-2,k)) &
                  + GAM*(wz(i,j+3,k)-wz(i,j-3,k)) &
                  + DEL*(wz(i,j+4,k)-wz(i,j-4,k)))*dxinv(2)

             difflux(i,j,k,imy) = eta*(vxx + FourThirds*vyy + vzz + OneThird*(uxy+wzy))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_diffterm_loop5)

    call build(bpt_diffterm_loop6, "bpt_diffterm_loop6")
    !$OMP PARALLEL DO PRIVATE(i,j,k,wxx,wyy,wzz,uxz,vyz)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             wxx = (CENTER*q(i,j,k,qw) &
                  + OFF1*(q(i+1,j,k,qw)+q(i-1,j,k,qw)) &
                  + OFF2*(q(i+2,j,k,qw)+q(i-2,j,k,qw)) &
                  + OFF3*(q(i+3,j,k,qw)+q(i-3,j,k,qw)) &
                  + OFF4*(q(i+4,j,k,qw)+q(i-4,j,k,qw)))*dxinv(1)**2

             wyy = (CENTER*q(i,j,k,qw) &
                  + OFF1*(q(i,j+1,k,qw)+q(i,j-1,k,qw)) &
                  + OFF2*(q(i,j+2,k,qw)+q(i,j-2,k,qw)) &
                  + OFF3*(q(i,j+3,k,qw)+q(i,j-3,k,qw)) &
                  + OFF4*(q(i,j+4,k,qw)+q(i,j-4,k,qw)))*dxinv(2)**2

             wzz = (CENTER*q(i,j,k,qw) &
                  + OFF1*(q(i,j,k+1,qw)+q(i,j,k-1,qw)) &
                  + OFF2*(q(i,j,k+2,qw)+q(i,j,k-2,qw)) &
                  + OFF3*(q(i,j,k+3,qw)+q(i,j,k-3,qw)) &
                  + OFF4*(q(i,j,k+4,qw)+q(i,j,k-4,qw)))*dxinv(3)**2

             uxz = (ALP*(ux(i,j,k+1)-ux(i,j,k-1)) &
                  + BET*(ux(i,j,k+2)-ux(i,j,k-2)) &
                  + GAM*(ux(i,j,k+3)-ux(i,j,k-3)) &
                  + DEL*(ux(i,j,k+4)-ux(i,j,k-4)))*dxinv(3)

             vyz = (ALP*(vy(i,j,k+1)-vy(i,j,k-1)) &
                  + BET*(vy(i,j,k+2)-vy(i,j,k-2)) &
                  + GAM*(vy(i,j,k+3)-vy(i,j,k-3)) &
                  + DEL*(vy(i,j,k+4)-vy(i,j,k-4)))*dxinv(3)

             difflux(i,j,k,imz) = eta*(wxx + wyy + FourThirds*wzz + OneThird*(uxz+vyz))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_diffterm_loop6)

    call build(bpt_diffterm_loop7, "bpt_diffterm_loop7")
    !$OMP PARALLEL DO PRIVATE(i,j,k,txx,tyy,tzz) &
    !$OMP PRIVATE(divu,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz,mechwork)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             txx = (CENTER*q(i,j,k,6) &
                  + OFF1*(q(i+1,j,k,6)+q(i-1,j,k,6)) &
                  + OFF2*(q(i+2,j,k,6)+q(i-2,j,k,6)) &
                  + OFF3*(q(i+3,j,k,6)+q(i-3,j,k,6)) &
                  + OFF4*(q(i+4,j,k,6)+q(i-4,j,k,6)))*dxinv(1)**2

             tyy = (CENTER*q(i,j,k,6) &
                  + OFF1*(q(i,j+1,k,6)+q(i,j-1,k,6)) &
                  + OFF2*(q(i,j+2,k,6)+q(i,j-2,k,6)) &
                  + OFF3*(q(i,j+3,k,6)+q(i,j-3,k,6)) &
                  + OFF4*(q(i,j+4,k,6)+q(i,j-4,k,6)))*dxinv(2)**2

             tzz = (CENTER*q(i,j,k,6) &
                  + OFF1*(q(i,j,k+1,6)+q(i,j,k-1,6)) &
                  + OFF2*(q(i,j,k+2,6)+q(i,j,k-2,6)) &
                  + OFF3*(q(i,j,k+3,6)+q(i,j,k-3,6)) &
                  + OFF4*(q(i,j,k+4,6)+q(i,j,k-4,6)))*dxinv(3)**2

             divu  = TwoThirds*(ux(i,j,k)+vy(i,j,k)+wz(i,j,k))
             tauxx = 2.d0*ux(i,j,k) - divu
             tauyy = 2.d0*vy(i,j,k) - divu
             tauzz = 2.d0*wz(i,j,k) - divu
             tauxy = uy(i,j,k)+vx(i,j,k)
             tauxz = uz(i,j,k)+wx(i,j,k)
             tauyz = vz(i,j,k)+wy(i,j,k)

             mechwork = tauxx*ux(i,j,k) + &
                        tauyy*vy(i,j,k) + &
                        tauzz*wz(i,j,k) + tauxy**2+tauxz**2+tauyz**2

             mechwork = eta*mechwork &
                  + difflux(i,j,k,imx)*q(i,j,k,qu) &
                  + difflux(i,j,k,imy)*q(i,j,k,qv) &
                  + difflux(i,j,k,imz)*q(i,j,k,qw)

             difflux(i,j,k,iene) = alam*(txx+tyy+tzz) + mechwork
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_diffterm_loop7)

    call destroy(bpt_diffterm)

  end subroutine diffterm

end module advance_module

