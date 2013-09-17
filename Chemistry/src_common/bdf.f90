!
! BDF (backward differentiation formula) time-stepping routines.
!
! 
!

module bdf_module
  use bdf_params_module
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
     integer :: order                     ! maximum order (1 to 6)
     integer :: max_steps                 ! maximum allowable number of steps
     integer :: max_iters                 ! maximum allowable number of newton iterations

     procedure(f_proc), pointer, nopass :: f
     procedure(J_proc), pointer, nopass :: J

     type(bdf_params) :: ctx

     real(dp) :: dt                       ! last step-size used
     real(dp), pointer :: rtol(:)         ! realtive tolerances
     real(dp), pointer :: atol(:)         ! absolute tolerances

     real(dp), pointer :: z(:,:)          ! nordsieck histroy array, indexed as (dof, order)
     real(dp), pointer :: z0(:,:)         ! nordsieck predictor array
     real(dp), pointer :: t(:)            ! times, t = [ t_n, t_{n-1}, ..., t_{n-k} ]
     real(dp), pointer :: l(:)            ! predictor/corrector update coefficients
     integer,  pointer :: A(:,:)          ! pascal matrix

     real(dp), pointer :: Jac(:,:)        ! jacobian matrix
     real(dp), pointer :: P(:,:)          ! newton iteration matrix
     real(dp), pointer :: b(:)            ! rhs workspace
  end type bdf_ts

  interface eye
     module procedure eye_r
     module procedure eye_i
  end interface

  interface
     subroutine f_proc(neq, y, t, yd, ctx)
       import dp, bdf_params
       integer,          intent(in)  :: neq
       real(dp),         intent(in)  :: y(neq), t
       real(dp),         intent(out) :: yd(neq)
       type(bdf_params), intent(in)  :: ctx
     end subroutine f_proc

     subroutine J_proc(neq, y, t, J, ctx)
       import dp, bdf_params
       integer,          intent(in)  :: neq
       real(dp),         intent(in)  :: y(neq), t
       real(dp),         intent(out) :: J(neq, neq)
       type(bdf_params), intent(in)  :: ctx
     end subroutine J_proc
  end interface

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

    integer  :: i, j, k, step
    real(dp) :: t, dt
    real(dp) :: y(neq), yd(neq), rhs(neq)

    !
    ! do time-stepping
    !
    y  = y0
    t  = t0
    dt = dt0

    do step = 1, ts%max_steps

       if (step == 1 .and. .not. restart) then
          k = 1
          ts%t(0:1) = [ t, t+dt ]
          call bdf_nordsieck_update_coeffs(k, ts%t(0:1), ts%l(0:1))
          call ts%f(neq, y, t, yd, ts%ctx)
          ts%z(:,0) = y
          ts%z(:,1) = dt * yd
       end if

       ! adjust dt
       ! if (t + dt > t1) dt = t1 - t

       ! update coeffs
       ! XXX

       ! predict next nordsieck array
       ! XXX: blas?
       do i = 0, k
          ts%z0(:,i) = 0          
          do j = i, k
             ts%z0(:,i) = ts%z0(:,i) + ts%A(i,j) * ts%z(:,j)
          end do
       end do

       ! solve y_n - dt f(y_n,t) = y - dt yd for y_n
       ! XXX: blas?
       rhs = ts%z0(:,0) - ts%z0(:,1) / ts%l(1)
       call bdf_solve(ts, neq, y, t+dt, dt/ts%l(1), rhs, restart, reuse)

       ! update nordsieck array
       do i = 0, k
          ts%z(:,i) = ts%z0(:,i) + (y - ts%z0(:,0)) * ts%l(i)
       end do

       t = t + dt

       if (t >= t1) exit

    end do

    y1 = y

  end subroutine bdf_advance

  subroutine bdf_solve(ts, neq, y, t, dt, a, restart, reuse)
    type(bdf_ts),     intent(inout) :: ts
    integer,          intent(in)    :: neq
    real(dp),         intent(inout) :: y(neq)
    real(dp),         intent(in)    :: dt, t, a(neq)
    logical,          intent(in)    :: restart, reuse

    include 'LinAlg.inc'

    integer  :: k, ipvt(neq), info
    ! XXX: pass all of these in... or just inline this???
    real(dp) :: dt_hat, c, s(neq), f(neq), delta, ewt(neq)

    real(dp), pointer :: P(:,:), J(:,:), b(:)
    P => ts%P; J => ts%Jac; b => ts%b

    dt_hat = dt

    call bdf_error_weight(ts, y, ewt)

    ! XXX: this doesn't save jacobians...

    !
    ! newton iteration
    !
    ! general form is:
    !   solve:   P s(k) = -c G(y(k)) for s(k)
    !   update:  y(k+1) = y(k) + s(k)
    ! where
    !   G(y) = y - dt * f(y,t) - a
    !
    do k = 1, ts%max_iters
      
       ! build iteration matrix and factor
       call eye(P)
       call ts%J(neq, y, t, J, ts%ctx)
       ! XXX: blas
       P = P - dt * J
       c = 2 * dt_hat / (dt + dt_hat)
       call dgefa(P, neq, neq, ipvt, info)

       ! solve using factorized iteration matrix
       call ts%f(neq, y, t, f, ts%ctx)
       ! XXX: blas
       b = -c * (y - dt * f - a)
       call dgesl(P, neq, neq, ipvt, b, 0)

       ! XXX: use an accumulated correction instead of updating y directly
       y = y + b

       delta = bdf_norm(b, ewt)

       if (delta < one) exit
    end do

  end subroutine bdf_solve

  !
  ! Compute Nordsieck update coefficients l based on times t.
  !
  ! See eqn. 5.2 of Jackson and Sacks-Davis (1980).
  !
  ! Note: 
  !
  !   1. The input vector t = [ t_n, t_{n-1}, ... t_{n-k} ] where we
  !      are advancing from step n-1 to step n.
  ! 
  !   2. The step size h_n = t_n - t_{n-1}.
  !
  subroutine bdf_nordsieck_update_coeffs(k, t, l)
    integer,  intent(in)  :: k
    real(dp), intent(in)  :: t(0:k)
    real(dp), intent(out) :: l(0:k)

    integer :: j
 
    l(0) = 1
    l(1) = xi_j(k, t, 1)
    do j = 2, k
       l = l + eoshift(l, 1) / xi_j(k, t, j)
    end do
  contains

    !
    ! Compute \xi_j, handle \xi_k specially.
    !
    function xi_j(k, t, j) result(xi)
      integer,  intent(in) :: k, j
      real(dp), intent(in) :: t(0:k)

      integer  :: i
      real(dp) :: xi, h, mu, a0

      h = t(0) - t(1)
      if (j == k) then
         a0 = 0
         do i = 1, k
            a0 = a0 - one/i
         end do
         mu = -a0 / h - qd(k, t) / q(k, t)
         xi = one / (mu * h)
      else
         xi = (t(0) - t(j)) / h
      end if
    end function xi_j

    !
    ! Compute q_{k-1}(t_n).
    !
    function q(k, t) result(r)
      integer,  intent(in) :: k
      real(dp), intent(in) :: t(0:k)
      real(dp) :: r
      integer  :: j
      r = 1
      if (k /= 1) then
         do j = 1, k-1
            r = r * (t(0) - t(j))
         end do
      end if
    end function q

    !
    ! Compute \dot{q}_{k-1}(t_n).
    !
    function qd(k, t) result(r)
      integer,  intent(in) :: k
      real(dp), intent(in) :: t(0:k)
      real(dp) :: r
      integer  :: j
      r = 0
      if (k /= 1) then
         do j = 1, k-1
            r = r + q(k, t) * (t(0) - t(j))
         end do
      end if
    end function qd
  end subroutine bdf_nordsieck_update_coeffs


  !
  ! Estimate initial step-size.  See sec. 3.2 of Brown, Byrne, and
  ! Hindmarsh.
  !
  ! function bdf_estimate_dt(ts, f, y0, t0, t1, ctx) result(dt)
  !   type(bdf_ts),     intent(inout) :: ts
  !   real(dp),         intent(in)    :: y0(:), t0, t1
  !   type(bdf_params), intent(in) :: ctx
  !   interface 
  !      subroutine f(y, t, yd, ctx)
  !        import dp, bdf_params
  !        real(dp), intent(in)  :: y(:), t
  !        real(dp), intent(out) :: yd(:)
  !        type(bdf_params), intent(in) :: ctx
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
    
  subroutine bdf_error_weight(ts, y, ewt)
    type(bdf_ts), intent(in)  :: ts
    real(dp),     intent(in)  :: y(:)
    real(dp),     intent(out) :: ewt(:)
    ewt = one / (ts%rtol * abs(y) + ts%atol)
  end subroutine bdf_error_weight

  function bdf_norm(y, ewt) result(r)
    real(dp), intent(in) :: y(:), ewt(:)
    real(dp) :: r
    r = sqrt(sum((y * ewt)**2)/size(y))
  end function bdf_norm

  !
  ! Build/destroy BDF time-stepper.
  !
  subroutine bdf_ts_build(ts, neq, f, J, rtol, atol, order)
    type(bdf_ts), intent(inout) :: ts
    integer,      intent(in   ) :: order, neq
    real(dp),     intent(in   ) :: rtol(neq), atol(neq)
    procedure(f_proc)           :: f
    procedure(J_proc)           :: J

    integer :: U(order, order)

    ! set options
    ts%neq  = neq
    ts%order = order

    ts%f => f
    ts%J => J

    ! set default options
    ts%max_steps = 10000
    ts%max_iters = 66

    ! pre-allocate work arrays etc
    allocate(ts%rtol(neq))
    allocate(ts%atol(neq))
    allocate(ts%z(neq, 0:order))
    allocate(ts%z0(neq, 0:order))
    allocate(ts%l(0:order))
    allocate(ts%t(0:order))
    allocate(ts%A(0:order, 0:order))
    allocate(ts%P(neq, neq))
    allocate(ts%Jac(neq, neq))
    allocate(ts%b(neq))

    ts%rtol = rtol
    ts%atol = atol

    ! build pascal matrix and update coefficients
    call bdf_upper_pascal(order+1, ts%A)

    ! XXX: l
  end subroutine bdf_ts_build

  subroutine bdf_ts_destroy(ts)
    type(bdf_ts), intent(inout) :: ts
    deallocate(ts%z,ts%z0,ts%t,ts%l,ts%A,ts%rtol,ts%atol,ts%P,ts%Jac,ts%b)
  end subroutine bdf_ts_destroy

  subroutine bdf_upper_pascal(q, A)
    integer, intent(in)  :: q
    integer, intent(out) :: A(q,q)

    integer :: k, U(q, q), Uk(q, q)

    ! build off diagonal matrix
    U = 0
    do k = 1, q-1
       U(k,k+1) = k
    end do

    ! compute A = exp(U)
    call eye(A)

    Uk = U
    do k = 1, q
       A  = A + Uk / factorial(k)
       Uk = matmul(U, Uk)
    end do

    ! do k = 1, q
    !    print *, A(k,:)
    ! end do

  contains

    recursive function factorial(n) result(r)
      integer, intent(in) :: n
      integer :: r
      if (n == 1) then
         r = 1
      else
         r = n * factorial(n-1)
      end if
    end function factorial

  end subroutine bdf_upper_pascal

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

end module bdf_module
