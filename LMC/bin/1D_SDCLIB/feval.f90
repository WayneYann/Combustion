!
! Method-of-lines discretization for LMC.
!
! 
!

module feval
  use iso_c_binding
  use encap
  use debug_module
  implicit none
contains


  !
  ! Compute MOL dU/dt for the advection piece of LMC.
  !
  subroutine advection(aofv, aofs, vel_old, scal_old, divu_old, press_old, dx, lo, hi, bc)
    use spec_module
    real*8  dx(0:nlevs-1)
    integer lo(0:nlevs-1)
    integer hi(0:nlevs-1)
    integer bc(0:nlevs-1,2)
    real*8  press_old(0:nlevs-1,-1:nfine+1)
    real*8    vel_old(0:nlevs-1,-2:nfine+1)
    real*8   scal_old(0:nlevs-1,-2:nfine+1,nscal)
    real*8   divu_old(0:nlevs-1,-1:nfine)

    real*8  dt(0:nlevs-1)
    real*8           gp(0:nlevs-1,-1:nfine)
    real*8       macvel(0:nlevs-1, 0:nfine  )
    real*8      veledge(0:nlevs-1, 0:nfine  )
    ! real*8     diff_old(0:nlevs-1,-1:nfine,  nscal)
    ! real*8     diff_new(0:nlevs-1,-1:nfine,  nscal)
    ! real*8     diff_hat(0:nlevs-1,-1:nfine,  nscal)
    ! real*8     diff_tmp(0:nlevs-1,-1:nfine,  nscal)
    real*8       tforce(0:nlevs-1,-1:nfine,  nscal)

    real*8, intent(inout) :: aofv(0:nlevs-1, 0:nfine-1)
    real*8, intent(inout) :: aofs(0:nlevs-1, 0:nfine-1,nscal)

    integer :: i

    aofv = 0
    aofs = 0

    print *,'... projecting'

    ! compute cell-centered grad pi from nodal pi
    do i=lo(0)-1,hi(0)+1
       gp(0,i) = (press_old(0,i+1) - press_old(0,i)) / dx(0)
    enddo

    ! compute edge velocities at current state
    dt = 0
    call pre_mac_predict(vel_old(0,:),scal_old(0,:,:),gp(0,:), &
         macvel(0,:),dx(0),dt(0),lo(0),hi(0),bc(0,:))

    ! project
    call macproj(macvel(0,:),scal_old(0,:,Density), &
         divu_old(0,:),dx,lo(0),hi(0),bc(0,:))

    call dsend(macvel(0,:), .false.)

    ! compute scalar advective flux
    print *,'... computing advective scalar flux'
    tforce = 0
    call scal_aofs(scal_old(0,:,:),macvel(0,:),aofs(0,:,:), &
                       divu_old(0,:),tforce(0,:,:),dx(0),dt(0), &
                       lo(0),hi(0),bc(0,:))

    ! compute velocity advective flux 
    print *,'... computing advective velocity flux'

    ! note that tforce is 0
    call vel_edge_states(vel_old(0,:),scal_old(0,:,Density),gp(0,:), &
         macvel(0,:),veledge(0,:),dx(0),dt(0), &
         tforce(0,0,:),lo(0),hi(0),bc(0,:))

    do i=lo(0),hi(0)
       aofv(0,i) = - ( (macvel(0,i+1)*veledge(0,i+1) - macvel(0,i)*veledge(0,i)) - &
            0.5d0*(macvel(0,i+1)-macvel(0,i))*(veledge(0,i)+veledge(0,i+1)) ) / dx(0) &
            - gp(0,i)/scal_old(0,i,Density)
    enddo

  end subroutine advection

  subroutine f1eval(fptr, qptr, t, state, ctx) bind(c)
    use spec_module
    type(c_ptr),    intent(in), value :: fptr, qptr, state, ctx
    real(c_double), intent(in), value :: t

    type(lmc_encap), pointer :: q, f
    integer :: lo(0:nlevs-1), hi(0:nlevs-1), bc(0:nlevs-1,2)
    real(8) :: dx(0:nlevs-1)
    
    call c_f_pointer(qptr, q)
    call c_f_pointer(fptr, f)

    print *, 'f1eval...'

    ! sick hacks...
    lo = 0
    hi = nx-1
    dx = (3.5)/dble(nx)
    bc(0,1) = 1
    bc(0,2) = 2

    f%scal = 0
    f%vel = 0

    call dsend(q%scal(0,:,density), .false.)
    call dsend(q%press(0,:), .false.)
    call dsend(q%vel(0,:), .false.)

    call advection(f%vel, f%scal, q%vel, q%scal, q%divu, q%press, dx, lo, hi, bc)

    ! XXX: set_bc_s, set_bc_v

    call dsend(f%scal(0,:,density), .false.)
    call dsend(f%vel(0,:), .true.)

  end subroutine f1eval

  subroutine f2eval(fptr, qptr, t, state, ctx) bind(c)
    use spec_module
    type(c_ptr),    intent(in), value :: fptr, qptr, state, ctx
    real(c_double), intent(in), value :: t

    type(lmc_encap), pointer :: q, f

    integer :: lo(0:nlevs-1), hi(0:nlevs-1), bc(0:nlevs-1,2)
    real(8) :: dx(0:nlevs-1)

    integer :: i

    real*8   beta(0:nlevs-1,-1:nfine  ,nscal)
    real*8   beta_for_Y(0:nlevs-1,-1:nfine  ,nscal)
    real*8   beta_for_Wbar(0:nlevs-1,-1:nfine  ,nscal)
    real*8       mu(0:nlevs-1,-1:nfine)

      real*8     diff(0:nlevs-1,-1:nfine,  nscal)
      real*8      gamma_lo(0:nlevs-1, 0:nfine-1,Nspec)
      real*8      gamma_hi(0:nlevs-1, 0:nfine-1,Nspec)
      real*8 diffdiff(0:nlevs-1,-1:nfine)
      real*8         visc(0:nlevs-1,-1:nfine)

    call c_f_pointer(qptr, q)
    call c_f_pointer(fptr, f)

    print *, 'f2eval...'

    f%vel = 0
    f%scal = 0

    ! sick hacks...
    lo = 0
    hi = nx-1
    dx = (3.5)/dble(nx)
    bc(0,1) = 1
    bc(0,2) = 2

    print *,'... compute lagged diff, D^{n+1,(k-1)}'

     ! compute transport coefficients
     !    rho D_m     (for species)
     !    lambda / cp (for enthalpy)
     !    lambda      (for temperature)
    call calc_diffusivities(q%scal(0,:,:),beta(0,:,:), &
         beta_for_Y(0,:,:), &
         beta_for_Wbar(0,:,:), &
         mu(0,:),lo(0),hi(0))

    ! compute div lambda grad T
    diff(0,:,Temp) = 0.d0
    call addDivLambdaGradT(q%scal(0,:,:),beta(0,:,:), &
         diff(0,:,Temp),dx(0),lo(0),hi(0))

    ! XXX: this eventually needs to end up in RhoH?

    ! compute a conservative div gamma_m
    ! save gamma_m for differential diffusion computation
    call get_spec_visc_terms(q%scal(0,:,:),beta(0,:,:), &
         diff(0,:,FirstSpec:), &
         gamma_lo(0,:,:),gamma_hi(0,:,:), &
         dx(0),lo(0),hi(0))

    do i = 0, nfine
       f%scal(0,i,FirstSpec:) = diff(0,i,FirstSpec:)
    end do

    ! compute div h_m Gamma_m
    ! we pass in conservative gamma_m via gamma
    ! we compute h_m using T from the scalar argument
    call get_diffdiff_terms(q%scal(0,:,:), &
         gamma_lo(0,:,:),gamma_hi(0,:,:), &
         diffdiff(0,:),dx(0),lo(0),hi(0))

    do i = 0, nfine
       f%scal(0,i,RhoH) = f%scal(0,i,RhoH) + diffdiff(0,i)
    end do


    !  ! instantaneous omegadot for divu calc
    ! do i=lo(0),hi(0)
    !    do n=1,Nspec
    !       C(n) = scal(0,i,FirstSpec+n-1)*invmwt(n)
    !    end do
    !    call CKWC(scal(0,i,Temp),C,IWRK,RWRK,WDOTK)
    !    do n=1,Nspec
    !       I_R_instant(0,i,n) = WDOTK(n)*mwt(n)
    !    end do
    ! end do

    ! ! divu
    ! call calc_divu(scal(0,:,:),beta(0,:,:),I_R_instant(0,:,:), &
    !      divu(0,:),dx(0),lo(0),hi(0))


    ! viscosity term
    call get_vel_visc_terms(q%vel(0,:),mu(0,:),f%vel(0,:), dx(0),lo(0),hi(0))

  end subroutine f2eval

  subroutine f2comp(fptr, qptr, t, dt, rhsptr, state, ctx) bind(c)
    type(c_ptr),    intent(in), value :: fptr, qptr, rhsptr, state, ctx
    real(c_double), intent(in), value :: t, dt

    type(lmc_encap), pointer :: q, f, rhs

    call c_f_pointer(qptr, q)
    call c_f_pointer(fptr, f)
    call c_f_pointer(rhsptr, rhs)

    print *, 'f2comp'

    q%vel = rhs%vel
    q%scal = rhs%scal

    f%vel = 0
    f%scal = 0

    
  end subroutine f2comp

end module feval
