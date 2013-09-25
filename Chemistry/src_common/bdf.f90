!
! BDF (backward differentiation formula) time-stepping routines.
!
! See
!
!   1. VODE: A variable-coefficient ODE solver; Brown, Byrne, and
!      Hindmarsh; SIAM J. Sci. Stat. Comput., vol. 10, no. 5, pp.
!      1035-1051, 1989.
!
!   2. An alternative implementation of variable step-size multistep
!      formulas for stiff ODES; Jackson and Sacks-Davis; ACM
!      Trans. Math. Soft., vol. 6, no. 3, pp. 295-318, 1980.
!
!   3. A polyalgorithm for the numerical solution of ODEs; Byrne and
!      Hindmarsh; ACM Trans. Math. Soft., vol. 1, no. 1, pp. 71-96,
!      1975.
!

module bdf
  use bdf_params
  implicit none

  integer, parameter  :: dp   = kind(1.d0)
  real(dp), parameter :: one  = 1.0_dp
  real(dp), parameter :: two  = 2.0_dp
  real(dp), parameter :: half = 0.5_dp

  !
  ! bdf time-stepper
  !
  type :: bdf_ts

     integer :: neq                       ! number of equations (degrees of freedom)
     integer :: max_order                 ! maximum order (1 to 6)
     integer :: max_steps                 ! maximum allowable number of steps
     integer :: max_iters                 ! maximum allowable number of newton iterations
     real(dp) :: dt_min                   ! minimum allowable step-size
     real(dp) :: eta_min                  ! minimum allowable step-size shrink factor
     real(dp) :: eta_max                  ! maximum allowable step-size growth factor
     real(dp) :: eta_thresh               ! step-size growth threshold

     procedure(f_proc), pointer, nopass :: f
     procedure(J_proc), pointer, nopass :: J
     type(bdf_ctx) :: ctx

     real(dp) :: t                        ! current time
     real(dp) :: dt                       ! current time step
     real(dp) :: dt_nwt                   ! dt used when building newton iteration matrix
     integer  :: k                        ! current order
     integer  :: n                        ! current step
     integer  :: age                      ! age of jacobian
     integer  :: max_age                  ! maximum age of jacobian
     integer  :: k_age                    ! number of steps taken at current order

     real(dp) :: tq(-1:2)                 ! error coefficients (test quality)
     ! indexed as: tq(-1) for order k-1 error est. coeff.
     !             tq(0)  for order k error est. coeff.
     !             tq(1)  for order k+1 error est. coeff.
     !             tq(2)  for order k+1 error est. coeff. used for second derivative

     real(dp), pointer :: rtol(:)         ! realtive tolerances
     real(dp), pointer :: atol(:)         ! absolute tolerances

     ! jacobian and newton matrices, may be resused
     real(dp), pointer :: Jac(:,:)        ! jacobian matrix
     real(dp), pointer :: P(:,:)          ! newton iteration matrix

     ! work-spaces
     real(dp), pointer :: z(:,:)          ! nordsieck histroy array, indexed as (dof, n)
     real(dp), pointer :: z0(:,:)         ! nordsieck predictor array
     real(dp), pointer :: h(:)            ! time steps, h = [ h_n, h_{n-1}, ..., h_{n-k} ]
     real(dp), pointer :: l(:)            ! predictor/corrector update coefficients
     real(dp), pointer :: y(:)            ! current y
     real(dp), pointer :: yd(:)           ! current \dot{y}
     real(dp), pointer :: rhs(:)          ! solver rhs
     real(dp), pointer :: e(:)            ! accumulated corrections
     real(dp), pointer :: ewt(:)          ! cached error weights
     real(dp), pointer :: b(:)            ! solver work space

     integer,  pointer :: ipvt(:)         ! pivots
     integer,  pointer :: A(:,:)          ! pascal matrix

     ! counters
     integer :: nfe                       ! number of function evaluations
     integer :: nje                       ! number of jacobian evaluations
     integer :: nit                       ! number of non-linear solver iterations
     integer :: nse                       ! number of non-linear solver errors

  end type bdf_ts

  interface
     subroutine f_proc(neq, y, t, yd, ctx)
       import dp, bdf_ctx
       integer,          intent(in)  :: neq
       real(dp),         intent(in)  :: y(neq), t
       real(dp),         intent(out) :: yd(neq)
       type(bdf_ctx), intent(in)  :: ctx
     end subroutine f_proc

     subroutine J_proc(neq, y, t, J, ctx)
       import dp, bdf_ctx
       integer,          intent(in)  :: neq
       real(dp),         intent(in)  :: y(neq), t
       real(dp),         intent(out) :: J(neq, neq)
       type(bdf_ctx), intent(in)  :: ctx
     end subroutine J_proc
  end interface

  private :: &
!       alpha0, alphahat0, &
       xi_j, xi_star, &
       ewts, norm, eye, eye_r, eye_i, factorial

  interface eye
     module procedure eye_r
     module procedure eye_i
  end interface

  interface build
     module procedure bdf_ts_build
  end interface build

  interface destroy
     module procedure bdf_ts_destroy
  end interface destroy

contains

  !
  ! Advance system from t0 to t1.
  !
  subroutine bdf_advance(ts, neq, y0, t0, y1, t1, dt0, restart, reuse)
    type(bdf_ts),     intent(inout) :: ts
    integer,          intent(in)    :: neq
    real(dp),         intent(in)    :: y0(neq), t0, t1, dt0
    real(dp),         intent(out)   :: y1(neq)
    logical,          intent(in)    :: restart, reuse

    include 'LinAlg.inc'

    integer  :: i, j, k, iter, info, nse
    real(dp) :: c, dt_adj, error, eta(-1:1), etamax

    logical :: rebuild
    logical :: verbose = .false.


    nse = 0

    if (.not. restart) then
       ts%nfe = 0
       ts%nje = 0
       ts%nit = 0
       ts%nse = 0

       ts%y  = y0
       ts%dt = dt0
       ts%n  = 1
    end if

    if (.not. restart) ts%age = 666
       
    ts%t = t0
    do k = 1, 666666666

       if (ts%n > ts%max_steps) exit

       if (verbose) print *, 'BDF: stepping: ', ts%n, ts%k, ts%dt

       call ewts(ts, ts%y, ts%ewt)

       if (k == 1 .and. .not. restart) then
          ts%k = 1
          ts%h = ts%dt
          call bdf_ts_update(ts)

          call ts%f(neq, ts%y, ts%t, ts%yd, ts%ctx)
          ts%nfe = ts%nfe + 1

          ts%z(:,0) = ts%y
          ts%z(:,1) = ts%dt * ts%yd

          ts%k_age = 0
       end if

       !
       ! predict
       !

       do i = 0, ts%k
          ts%z0(:,i) = 0          
          do j = i, ts%k
             ts%z0(:,i) = ts%z0(:,i) + ts%A(i,j) * ts%z(:,j)
          end do
       end do

       !
       ! solve y_n - dt f(y_n,t) = y - dt yd for y_n
       !

       ! newton iteration general form is:
       !   solve:   P x = -c G(y(k)) for x
       !   update:  y(k+1) = y(k) + x
       ! where
       !   G(y) = y - dt * f(y,t) - rhs

       ts%e   = 0
       ts%rhs = ts%z0(:,0) - ts%z0(:,1) / ts%l(1)
       dt_adj = ts%dt / ts%l(1)
       ts%y   = ts%z0(:,0)

       if (ts%age > ts%max_age) rebuild = .true.

       do iter = 1, ts%max_iters

          if (rebuild) then
             ! build iteration matrix and factor
             call eye(ts%P)
             call ts%J(neq, ts%y, ts%t, ts%Jac, ts%ctx)
             ts%P = ts%P - dt_adj * ts%Jac
             call dgefa(ts%P, neq, neq, ts%ipvt, info)

             ts%nje    = ts%nje + 1
             ts%dt_nwt = dt_adj
             ts%age    = 0
             rebuild   = .false.
          end if

          ! solve using factorized iteration matrix
          call ts%f(neq, ts%y, ts%t, ts%yd, ts%ctx)
          ts%nfe = ts%nfe + 1
          ts%nit = ts%nit + 1

          c    = 2 * ts%dt_nwt / (dt_adj + ts%dt_nwt)
          ts%b = c * (ts%rhs - ts%y + dt_adj * ts%yd)
          call dgesl(ts%P, neq, neq, ts%ipvt, ts%b, 0)

          ts%e = ts%e + ts%b
          if (norm(ts%b, ts%ewt) < one) exit
          ts%y = ts%z0(:,0) + ts%e
       end do

       ts%age = ts%age + 1


       !
       ! retry if the solver didn't converge or the error estimate is too large
       !

       ! if solver failed many times, bail...
       if (iter >= ts%max_iters .and. nse > 7) then
          ! XXX: signal an error of some kind
          stop "BDF SOLVER FAILED LOTS OF TIMES IN A ROW"
       end if

       ! if solver failed to converge, shrink dt and try again
       if (iter >= ts%max_iters) then
          rebuild = .true.
          ts%nse = ts%nse + 1
          nse    = nse + 1
          eta(0) = 0.25d0
          call rescale_timestep(ts, eta(0))
          if (verbose) print *, 'BDF: solver failed'
          cycle
       else
          nse = 0
       end if

       error = ts%tq(0) * norm(ts%e, ts%ewt)

       ! if local error is fairly large, shrink dt and try again
       if (error > one) then
          eta(0) = one / ( (6.d0 * error) ** (one / ts%k) + 1.d-6 )
          if (verbose) print *, 'BDF: error tolerance exceeded'
          call rescale_timestep(ts, eta(0))
          cycle
       end if

       !
       ! correct
       !

       do i = 0, ts%k
          ts%z(:,i) = ts%z0(:,i) + ts%e * ts%l(i)
       end do

       if (ts%t >= t1) exit

       ts%h    = eoshift(ts%h, -1)
       ts%h(0) = ts%dt
       ts%t    = ts%t + ts%dt
       ts%n    = ts%n + 1

       !
       ! compute eta values for changing step-size and/or order
       !

       eta = 0
       if (ts%max_order == 1 .or. ts%k_age < ts%k) then
          eta(0) = one / ( (6.d0 * error) ** (one / ts%k) + 1.d-6 )
       else
          if (ts%k > 1) then
             error   = ts%tq(-1) * norm(ts%z(:,ts%k), ts%ewt)
             eta(-1) = one / ( (6.d0 * error) ** (one / ts%k) + 1.d-6 )
          end if
          if (ts%k < ts%max_order) then
             c = ts%tq(2) * (ts%h(0) / ts%h(2)) ** (ts%k+1)
             error  = ts%tq(1) * norm(ts%e + c * ts%z(:,2), ts%ewt)
             eta(1) = one / ( (10.d0 * error) ** (one / (ts%k+2)) + 1.d-6 )
          end if
       end if

       ts%k_age = ts%k_age + 1

       etamax = maxval(eta)
       if (etamax > ts%eta_thresh) then
          if (etamax == eta(-1)) then
             call decrease_order(ts)
          else if (etamax == eta(1)) then
             call increase_order(ts)
          end if
          call rescale_timestep(ts, etamax)
       end if

       !
       ! final adjustments to time step...
       !

       ! if (t + dt > t1) dt = t1 - t

    end do

    if (verbose) print *, 'BDF: done    : ', ts%n, k, ts%nfe, ts%nje, ts%nit, ts%nse
    y1 = ts%y

  end subroutine bdf_advance

  !
  ! Rescale time-step.
  !
  ! This consists of:
  !   1. bound eta to honor eta_min, eta_max, and dt_min
  !   2. scale dt and adjust time array t accordingly
  !   3. rescalel Nordsieck history array
  !   4. recompute Nordsieck update coefficients
  !   5. recompute local error coefficients
  !
  subroutine rescale_timestep(ts, eta)
    type(bdf_ts), intent(inout) :: ts
    real(dp),     intent(inout) :: eta
    integer :: i

    eta = max(eta, ts%dt_min / ts%dt, ts%eta_min)
    eta = min(eta, ts%eta_max)

    ts%dt   = eta * ts%dt
    ts%h(0) = ts%dt

    do i = 1, ts%k
       ts%z(:,i) = eta**i * ts%z(:,i)
    end do

    call bdf_ts_update(ts)
  end subroutine rescale_timestep

  !
  ! Decrease order.
  !
  subroutine decrease_order(ts)
    type(bdf_ts), intent(inout) :: ts
    integer  :: j
    real(dp) :: c(0:ts%k)

    if (ts%k > 2) then
       c = 0
       c(2) = 1
       do j = 1, ts%k-2
          c = eoshift(c, -1) + c * xi_j(ts%h, j)
       end do

       do j = 2, ts%k-1
          ts%z(:,j) = ts%z(:,j) - c(j) * ts%z(:,ts%k)
       end do
    end if

    ts%z(:,ts%k) = 0
    ts%k         = ts%k - 1
    ts%k_age     = 0
  end subroutine decrease_order

  !
  ! Increase order.
  !
  subroutine increase_order(ts)
    type(bdf_ts), intent(inout) :: ts
    integer  :: j
    real(dp) :: c(0:ts%k)

    c = 0
    c(2) = 1
    do j = 1, ts%k-2
       c = eoshift(c, -1) + c * xi_j(ts%h, j)
    end do

    ts%z(:,ts%k+1) = 0
    do j = 2, ts%k+1
       ts%z(:,j) = ts%z(:,j) + c(j) * ts%e
    end do

    ts%k = ts%k + 1
    ts%k_age = 0
  end subroutine increase_order

  !
  ! Compute Nordsieck update coefficients l and error coefficients tq.
  !
  ! Regarding the l coefficients, see section 5, and in particular
  ! eqn. 5.2, of Jackson and Sacks-Davis (1980).
  !
  ! Regarding the error coefficients tq, these have been adapted from cvode.
  !
  ! Note: 
  !
  !   1. The input vector t = [ t_n, t_{n-1}, ... t_{n-k} ] where we
  !      are advancing from step n-1 to step n.
  ! 
  !   2. The step size h_n = t_n - t_{n-1}.
  !
  subroutine bdf_ts_update(ts)
    type(bdf_ts), intent(inout) :: ts

    integer  :: j
    real(dp) :: a0, a0hat, a1, a2, a3, a4, a5, a6, xistar_inv, xi_inv, c

    ts%l  = 0
    ts%tq = 0

    ! compute l vector
    ts%l(0) = 1
    ts%l(1) = xi_j(ts%h, 1)
    if (ts%k > 1) then
       do j = 2, ts%k-1
          ts%l = ts%l + eoshift(ts%l, -1) / xi_j(ts%h, j)
       end do
       ts%l = ts%l + eoshift(ts%l, -1) / xi_star(ts%k, ts%h)
    end if

    ! compute error coefficients (adapted from cvode)
    a0hat = alphahat0(ts%k, ts%h)
    a0    = alpha0(ts%k)

    xi_inv     = one
    xistar_inv = one
    if (ts%k > 1) then
       xi_inv     = one / xi_j(ts%h, ts%k)
       xistar_inv = one / xi_star(ts%k, ts%h)
    end if

    a1 = one - a0hat + a0
    a2 = one + ts%k * a1
    ts%tq(0) = abs(a1 / (a0 * a2))
    ts%tq(2) = abs(a2 * xistar_inv / (ts%l(ts%k) * xi_inv))
    if (ts%k > 1) then
       c  = xistar_inv / ts%l(ts%k)
       a3 = a0 + one / ts%k
       a4 = a0hat + xi_inv
       ts%tq(-1) = abs(c * (one - a4 + a3) / a3)
    else
       ts%tq(-1) = one
    end if

    xi_inv = ts%h(0) / sum(ts%h(0:ts%k))
    a5 = a0 - one / (ts%k+1)
    a6 = a0hat - xi_inv
    ts%tq(1) = abs((one - a6 + a5) / a2 / (xi_inv * (ts%k+2) * a5))

    ! tq[4] = nlscoef / tq[2];

  end subroutine bdf_ts_update

  !
  ! Return $\alpha_0$.
  !
  function alpha0(k) result(a0)
    integer,  intent(in) :: k
    real(dp) :: a0
    integer  :: i
    a0 = -1
    do i = 2, k
       a0 = a0 - 1._dp/i
    end do
  end function alpha0

  !
  ! Return $\hat{\alpha}_{n,0}$.
  !
  function alphahat0(k, h) result(a0)
    integer,  intent(in) :: k
    real(dp), intent(in) :: h(0:k)
    real(dp) :: a0
    integer  :: i
    a0 = -1
    do i = 2, k
       a0 = a0 - h(0) / sum(h(0:i-1))
    end do
  end function alphahat0

  !
  ! Return $\xi^*_k$.
  !
  function xi_star(k, h) result(xi)
    integer,  intent(in) :: k
    real(dp), intent(in) :: h(0:)
    real(dp) :: xi, mu
    mu = -alpha0(k) / h(0) - qd(k-1, h) / q(k-1, h)
    xi = one / (mu * h(0))
  contains
    !
    ! Return $q_{k}(t_n)$.
    !
    function q(k, h) result(r)
      integer,  intent(in) :: k
      real(dp), intent(in) :: h(0:k)
      real(dp) :: r
      integer  :: j
      r = 1
      if (k > 0) then
         do j = 1, k
            r = r * sum(h(0:j-1))
         end do
      end if
    end function q

    !
    ! Return $\dot{q}_{k}(t_n)$.
    !
    function qd(k, h) result(r)
      integer,  intent(in) :: k
      real(dp), intent(in) :: h(0:k)
      real(dp) :: r, a
      integer  :: j, i
      r = 0
      do i = 1, k
         a = 1
         do j = 1, k
            if (j == i) cycle
            a = a * sum(h(0:j-1))
         end do
         r = r + a
      end do
    end function qd
  end function xi_star

  !
  ! Return $\xi_j$.
  !
  function xi_j(h, j) result(xi)
    integer,  intent(in) :: j
    real(dp), intent(in) :: h(0:)
    real(dp) :: xi
    xi = sum(h(0:j-1)) / h(0)
  end function xi_j


  !
  ! Estimate initial step-size.  See sec. 3.2 of Brown, Byrne, and
  ! Hindmarsh.
  !
  ! function bdf_estimate_dt(ts, f, y0, t0, t1, ctx) result(dt)
  !   type(bdf_ts),     intent(inout) :: ts
  !   real(dp),         intent(in)    :: y0(:), t0, t1
  !   type(bdf_ctx), intent(in) :: ctx
  !   interface 
  !      subroutine f(y, t, yd, ctx)
  !        import dp, bdf_ctx
  !        real(dp), intent(in)  :: y(:), t
  !        real(dp), intent(out) :: yd(:)
  !        type(bdf_ctx), intent(in) :: ctx
  !      end subroutine f
  !   end interface

  !   real(dp) :: dt
  !   real(dp) :: h, hnew, hl, hu, wrms, yd(size(y0)), ydd(size(y0)), ewt(size(y0))
  !   integer  :: i, k

  !   hl = 100.0_dp * epsilon(hl) * max(abs(t0), abs(t1))
  !   hu =   0.1_dp * abs(t0 - t1)

  !   ! reduce hu
  !   call f(y0, t0, yd, ctx)
  !   if (any(hu * abs(yd) > 0.1_dp * abs(y0) + ts%atol)) then
  !      hu = minval((0.1_dp * abs(y0) + ts%atol) / abs(yd))
  !   end if

  !   h = sqrt(hl*hu)
  !   if (hu < hl) then
  !      dt = h
  !      return
  !   end if

  !   call bdf_error_weight(ts, y0, ewt)
  !   do k = 1, 4
  !      call f(y0 + dt*yd, t0 + dt, ydd, ctx)
  !      wrms = bdf_norm(half * dt * (ydd - yd) / h, ewt)
  !      if (wrms*hu**2 > two) then
  !         hnew = sqrt(two / wrms)
  !      else
  !         hnew = sqrt(h * hu)
  !      end if
  !      if (hnew/h > half .and. hnew/h < two) exit
  !      h = hnew
  !   end do
  !   dt = half * hnew 
  ! end function bdf_estimate_dt
    
  ! 
  ! Pre-compute error weights.
  !
  subroutine ewts(ts, y, ewt)
    type(bdf_ts), intent(in)  :: ts
    real(dp),     intent(in)  :: y(:)
    real(dp),     intent(out) :: ewt(:)
    ewt = one / (ts%rtol * abs(y) + ts%atol)
  end subroutine ewts

  !
  ! Compute weighted norm of y.
  !
  function norm(y, ewt) result(r)
    real(dp), intent(in) :: y(:), ewt(:)
    real(dp) :: r
    r = sqrt(sum((y * ewt)**2)/size(y))
  end function norm

  !
  ! Build/destroy BDF time-stepper.
  !

  subroutine bdf_ts_build(ts, neq, f, J, rtol, atol, max_order)
    type(bdf_ts), intent(inout) :: ts
    integer,      intent(in   ) :: max_order, neq
    real(dp),     intent(in   ) :: rtol(neq), atol(neq)
    procedure(f_proc)           :: f
    procedure(J_proc)           :: J

    integer :: k, U(max_order+1, max_order+1), Uk(max_order+1, max_order+1)

    allocate(ts%rtol(neq))
    allocate(ts%atol(neq))
    allocate(ts%z(neq, 0:max_order))
    allocate(ts%z0(neq, 0:max_order))
    allocate(ts%l(0:max_order))
    allocate(ts%h(0:max_order))
    allocate(ts%A(0:max_order, 0:max_order))
    allocate(ts%P(neq, neq))
    allocate(ts%Jac(neq, neq))
    allocate(ts%y(neq))
    allocate(ts%yd(neq))
    allocate(ts%rhs(neq))
    allocate(ts%e(neq))
    allocate(ts%ewt(neq))
    allocate(ts%b(neq))
    allocate(ts%ipvt(neq))

    ts%max_order  = max_order
    ts%max_steps  = 1000000
    ts%max_iters  = 10
    ts%dt_min     = epsilon(ts%dt_min)
    ts%eta_min    = 0.2_dp
    ! ts%eta_max    = 1.5_dp
    ts%eta_max    = 1000000.0_dp
    ! ts%eta_thresh = 1.5_dp
    ts%eta_thresh = 1.1_dp
    ts%max_age    = 30

    ts%k = -1
    ts%f => f
    ts%J => J
    ts%rtol = rtol
    ts%atol = atol

    ! build pascal matrix A using A = exp(U)
    U = 0
    do k = 1, max_order
       U(k,k+1) = k
    end do
    Uk = U

    call eye(ts%A)
    do k = 1, max_order+1
       ts%A  = ts%A + Uk / factorial(k)
       Uk = matmul(U, Uk)
    end do

  end subroutine bdf_ts_build

  subroutine bdf_ts_destroy(ts)
    type(bdf_ts), intent(inout) :: ts
    deallocate(ts%z,ts%z0,ts%h,ts%l,ts%A,ts%rtol,ts%atol,ts%P,ts%Jac,ts%y,ts%yd,ts%rhs,ts%e,ts%ewt,ts%b,ts%ipvt)
  end subroutine bdf_ts_destroy

  !
  ! Various misc. helper functions
  !

  subroutine eye_r(A)
    real(dp), intent(inout) :: A(:,:)
    integer :: i
    A = 0
    do i = 1, size(A, 1)
       A(i,i) = 1
    end do
  end subroutine eye_r

  subroutine eye_i(A)
    integer, intent(inout) :: A(:,:)
    integer :: i
    A = 0
    do i = 1, size(A, 1)
       A(i,i) = 1
    end do
  end subroutine eye_i

  recursive function factorial(n) result(r)
    integer, intent(in) :: n
    integer :: r
    if (n == 1) then
       r = 1
    else
       r = n * factorial(n-1)
    end if
  end function factorial

end module bdf
