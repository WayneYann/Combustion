!
! Method-of-lines discretization for LMC.
!

module feval
  use encap
  use debug
  use lmc
  use probin
  use kernels
  implicit none
contains

  !
  ! Compute MOL dU/dt for the advection piece of LMC.
  !
  subroutine f1eval(fptr, qptr, t, state, ctx) bind(c)
    type(c_ptr),    intent(in   ), value :: fptr, qptr, state, ctx
    real(c_double), intent(in   ), value :: t

    type(lmc_encap), pointer :: q, f

    integer          :: i
    double precision :: gp(lo-1:hi+1)
    double precision :: divu(lo-1:hi+1)
    double precision :: omegadot(lo:hi,nspec)
    double precision :: edgevel(lo:hi+1)
    double precision :: macvel(lo:hi+1)
    double precision :: aofs(lo:hi,nscal)
    double precision :: aofv(lo:hi)
    double precision :: beta(lo-1:hi+1,nscal)
    double precision :: beta_for_Y(lo-1:hi+1,nscal)
    double precision :: beta_for_Wbar(lo-1:hi+1,nscal)
    double precision :: mu(lo-1:hi+1)

    call c_f_pointer(qptr, q)
    call c_f_pointer(fptr, f)

    print *, 'explicit eval ...'

    f%scal = 0; f%vel = 0; f%press = 0; aofs = 0; aofv = 0

    print *, '   ... computing divu'
    call calc_diffusivities(q%scal, beta, beta_for_Y, beta_for_Wbar, mu, lo, hi)
    call calc_omega_dot(q%scal, omegadot, lo, hi)
    call calc_divu(q%scal, beta, omegadot, divu, dx, lo, hi)

    print *, '   ... projecting'
    call calc_edge_vel(edgevel, q%vel, lo, hi, bc)
    macvel = edgevel
    call mac_project(macvel, q%scal(:,density), divu, dx, lo, hi, bc)

    print *, '   ... computing advective scalar flux'
    call scal_aofs(q%scal, macvel, aofs, divu, dx, lo, hi, bc)

    print *, '   ... computing advective velocity flux'
    call calc_grad_pi(gp, q%press, lo, hi, dx)
    do i = lo, hi
       aofv(i) = - ( (macvel(i+1)*edgevel(i+1) - macvel(i)*edgevel(i)) &
                     - 0.5d0 * (macvel(i+1)-macvel(i))*(edgevel(i)+edgevel(i+1)) &
                   ) / dx - gp(i) / q%scal(i,Density)
    end do

    ! XXX: set_bc_s, set_bc_v, need to add a hook to take care of this...

    f%scal(lo:hi,:) = aofs
    f%vel(lo:hi) = aofv

    call plot1(q%scal(:,density), 1, 'with lines title "density"', .false.)
    call plot1(q%press(:), 2, 'with lines title "pressure"', .false.)
    call plot2(q%vel(lo:hi), macvel(lo:hi), 4, &
         "u 1:2 w l title 'velocity', '' u 1:3 w l title 'macvel'", .true.)
  end subroutine f1eval

  !
  ! Compute MOL dU/dt for the diffusion and reaction pieces of LMC.
  !
  subroutine f2eval(fptr, qptr, t, state, ctx) bind(c)
    type(c_ptr),    intent(in   ), value :: fptr, qptr, state, ctx
    real(c_double), intent(in   ), value :: t

    type(lmc_encap), pointer :: q, f

    integer          :: i
    double precision :: diff(lo-1:hi+1,nscal)
    double precision :: diffdiff(lo-1:hi+1)
    double precision :: divu(lo-1:hi+1)
    double precision :: dofv(lo-1:hi+1)
    double precision :: omegadot(lo:hi,nspec)
    double precision :: beta(lo-1:hi+1,nscal)
    double precision :: beta_for_Y(lo-1:hi+1,nscal)
    double precision :: beta_for_Wbar(lo-1:hi+1,nscal)
    double precision :: mu(lo-1:hi+1)
    double precision :: gamma_lo(lo:hi,nspec), gamma_hi(lo:hi,nspec)

    print *, 'implicit eval ...'

    call c_f_pointer(qptr, q)
    call c_f_pointer(fptr, f)

    f%vel = 0; f%scal = 0; f%press = 0

    print *,'   ... computing divu'

    ! compute transport coefficients
    !    rho D_m     (for species)
    !    lambda / cp (for enthalpy)
    !    lambda      (for temperature)
    call calc_diffusivities(q%scal, beta, beta_for_Y, beta_for_Wbar, mu, lo, hi)
    call calc_omega_dot(q%scal, omegadot, lo, hi)
    call calc_divu(q%scal,beta,omegadot,divu,dx,lo,hi)

    print *,'   ... computing grad T and diffdiff'

    ! compute div lambda grad T
    diff(:,Temp) = 0.d0
    call addDivLambdaGradT(q%scal, beta, diff(:,Temp), dx, lo, hi)

    ! ! XXX: this eventually needs to end up in RhoH?

    ! compute a conservative div gamma_m
    ! save gamma_m for differential diffusion computation
    call get_spec_visc_terms(q%scal, beta, diff(:,FirstSpec:), gamma_lo, gamma_hi, dx, lo, hi)

    do i = lo, hi
       f%scal(i,FirstSpec:) = diff(i,FirstSpec:)
    end do

    ! "diffdiff" means "differential diffusion", which corresponds to
    ! sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ] in equation (3)

    ! compute div h_m Gamma_m
    ! we pass in conservative gamma_m via gamma
    ! we compute h_m using T from the scalar argument
    call get_diffdiff_terms(q%scal, gamma_lo, gamma_hi, diffdiff, dx, lo, hi)

    do i = lo, hi
       f%scal(i,RhoH) = f%scal(i,RhoH) + diffdiff(i)
    end do

    ! viscosity term
    call get_vel_visc_terms(q%vel, mu, dofv, dx, lo, hi)

    f%vel(lo-1:hi+1) = dofv

  end subroutine f2eval

  subroutine f2comp(fptr, qptr, t, dt, rhsptr, state, ctx) bind(c)
    type(c_ptr),    intent(in   ), value :: fptr, qptr, rhsptr, state, ctx
    real(c_double), intent(in   ), value :: t, dt

    type(lmc_encap), pointer :: q, f, rhs

    print *, 'implicit solve ...'

    call c_f_pointer(qptr, q)
    call c_f_pointer(fptr, f)
    call c_f_pointer(rhsptr, rhs)

    q%vel  = rhs%vel
    q%scal = rhs%scal

    f%vel  = 0
    f%scal = 0
  end subroutine f2comp

end module feval
