!
! Method-of-lines discretization for LMC.
!

module feval
  use iso_c_binding
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
    double precision :: edgevel(lo:hi+1)
    double precision :: macvel(lo:hi+1)
    double precision :: aofs(lo:hi,0:nscal)
    double precision :: aofv(lo:hi)

    call c_f_pointer(qptr, q)
    call c_f_pointer(fptr, f)

    print *, 'explicit eval ...'

    aofs = 0
    aofv = 0

    call dsend(q%scal(:,density), .false.)
    call dsend(q%press(:), .false.)
    call dsend(q%vel(:), .false.)

    call calc_grad_pi(gp, q%press, lo, hi, dx)
    call calc_edge_vel(edgevel, q%vel, lo, hi, bc)

    print *, '   ... computing divu'
    ! call calc_divu()

    print *, '   ... projecting'
    macvel = edgevel
    call mac_project(macvel, q%scal(:,density), divu, dx, lo, hi, bc)

    print *, '   ... computing advective scalar flux'
    call scal_aofs(q%scal, macvel, aofs, divu, dx, lo, hi, bc)

    print *, '   ... computing advective velocity flux'
    do i = lo, hi
       aofv(i) = - ( (macvel(i+1)*edgevel(i+1) - macvel(i)*edgevel(i)) &
                     - 0.5d0 * (macvel(i+1)-macvel(i))*(edgevel(i)+edgevel(i+1)) &
                   ) / dx - gp(i) / q%scal(i,Density)
    end do

    ! XXX: set_bc_s, set_bc_v, need to add a hook to take care of this...

    f%scal = 0
    f%scal(lo:hi,0:nscal) = aofs

    f%vel = 0
    f%vel(lo:hi) = aofv

    call dsend(macvel, .false.)
    call dsend(f%scal(:,density), .false.)
    call dsend(f%vel, .true.)
  end subroutine f1eval

  subroutine f2eval(fptr, qptr, t, state, ctx) bind(c)
    type(c_ptr),    intent(in   ), value :: fptr, qptr, state, ctx
    real(c_double), intent(in   ), value :: t

    type(lmc_encap), pointer :: q, f


    ! integer                  :: lo(0:nlevs-1), hi(0:nlevs-1), bc(0:nlevs-1,2)
    ! real(8)                  :: dx(0:nlevs-1)

    ! integer :: i

    ! real*8   beta(0:nlevs-1,-1:nfine  ,nscal)
    ! real*8   beta_for_Y(0:nlevs-1,-1:nfine  ,nscal)
    ! real*8   beta_for_Wbar(0:nlevs-1,-1:nfine  ,nscal)
    ! real*8       mu(0:nlevs-1,-1:nfine)

    !   real*8     diff(0:nlevs-1,-1:nfine,  nscal)
    !   real*8      gamma_lo(0:nlevs-1, 0:nfine-1,Nspec)
    !   real*8      gamma_hi(0:nlevs-1, 0:nfine-1,Nspec)
    !   real*8 diffdiff(0:nlevs-1,-1:nfine)
    !   real*8         visc(0:nlevs-1,-1:nfine)

    print *, 'implicit eval ...'

    call c_f_pointer(qptr, q)
    call c_f_pointer(fptr, f)

    f%vel  = 0
    f%scal = 0

    ! ! sick hacks...
    ! lo = 0
    ! hi = nx-1
    ! dx = (3.5)/dble(nx)
    ! bc(0,1) = 1
    ! bc(0,2) = 2

    print *,'   ... computing diff'

    !  ! compute transport coefficients
    !  !    rho D_m     (for species)
    !  !    lambda / cp (for enthalpy)
    !  !    lambda      (for temperature)
    ! call calc_diffusivities(q%scal(0,:,:),beta(0,:,:), &
    !      beta_for_Y(0,:,:), &
    !      beta_for_Wbar(0,:,:), &
    !      mu(0,:),lo(0),hi(0))

    ! ! compute div lambda grad T
    ! diff(0,:,Temp) = 0.d0
    ! call addDivLambdaGradT(q%scal(0,:,:),beta(0,:,:), &
    !      diff(0,:,Temp),dx(0),lo(0),hi(0))

    ! ! XXX: this eventually needs to end up in RhoH?

    ! ! compute a conservative div gamma_m
    ! ! save gamma_m for differential diffusion computation
    ! call get_spec_visc_terms(q%scal(0,:,:),beta(0,:,:), &
    !      diff(0,:,FirstSpec:), &
    !      gamma_lo(0,:,:),gamma_hi(0,:,:), &
    !      dx(0),lo(0),hi(0))

    ! do i = 0, nfine
    !    f%scal(0,i,FirstSpec:) = diff(0,i,FirstSpec:)
    ! end do

    ! ! compute div h_m Gamma_m
    ! ! we pass in conservative gamma_m via gamma
    ! ! we compute h_m using T from the scalar argument
    ! call get_diffdiff_terms(q%scal(0,:,:), &
    !      gamma_lo(0,:,:),gamma_hi(0,:,:), &
    !      diffdiff(0,:),dx(0),lo(0),hi(0))

    ! do i = 0, nfine
    !    f%scal(0,i,RhoH) = f%scal(0,i,RhoH) + diffdiff(0,i)
    ! end do


    ! !  ! instantaneous omegadot for divu calc
    ! ! do i=lo(0),hi(0)
    ! !    do n=1,Nspec
    ! !       C(n) = scal(0,i,FirstSpec+n-1)*invmwt(n)
    ! !    end do
    ! !    call CKWC(scal(0,i,Temp),C,IWRK,RWRK,WDOTK)
    ! !    do n=1,Nspec
    ! !       I_R_instant(0,i,n) = WDOTK(n)*mwt(n)
    ! !    end do
    ! ! end do

    ! ! ! divu
    ! ! call calc_divu(scal(0,:,:),beta(0,:,:),I_R_instant(0,:,:), &
    ! !      divu(0,:),dx(0),lo(0),hi(0))


    ! ! viscosity term
    ! call get_vel_visc_terms(q%vel(0,:),mu(0,:),f%vel(0,:), dx(0),lo(0),hi(0))

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





  ! subroutine calc_divu(scal,beta,I_R,divu,dx,lo,hi)
  !   implicit none
  !   real*8 scal(-2:nfine+1,nscal)
  !   real*8 beta(-1:nfine  ,nscal)
  !   real*8  I_R(-1:nfine  ,0:Nspec)
  !   real*8 divu(-1:nfine)
  !   real*8 dx
  !   integer lo,hi

  !   real*8 Y(Nspec)
  !   real*8 HK(Nspec)
  !   real*8 cpmix,mwmix

  !   real*8 diff(-1:nfine,nscal)
  !   real*8 diffdiff(-1:nfine)

  !   real*8 RWRK,rho,T
  !   integer IWRK,i,n

  !   real*8 gamma_lo(0:nfine-1,Nspec)
  !   real*8 gamma_hi(0:nfine-1,Nspec)

  !   ! compute Gamma_m
  !   call get_spec_visc_terms(scal,beta,diff(:,FirstSpec:), gamma_lo,gamma_hi,dx,lo,hi)

  !   ! compute div lambda grad T
  !   diff(:,Temp) = 0.d0
  !   call addDivLambdaGradT(scal,beta,diff(:,Temp),dx,lo,hi)

  !   ! compute div h_m Gamma_m
  !   call get_diffdiff_terms(scal,gamma_lo,gamma_hi, diffdiff,dx,lo,hi)

  !   ! combine div lambda grad T + div h_m Gamma_m
  !   do i=lo,hi
  !      diff(i,Temp) = diff(i,Temp) + diffdiff(i)
  !   end do

  !   do i=lo,hi
  !      rho = scal(i,Density)
  !      do n = 1,Nspec
  !         Y(n) = scal(i,FirstSpec + n - 1) / rho
  !      enddo
  !      T = scal(i,Temp)
  !      call CKMMWY(Y,IWRK,RWRK,mwmix)
  !      call CKCPBS(T,Y,IWRK,RWRK,cpmix)
  !      call CKHMS(T,IWRK,RWRK,HK)

  !      divu(i) = diff(i,Temp)/(rho*cpmix*T)

  !      do n=1,Nspec
  !         divu(i) = divu(i) + (diff(i,FirstSpec+n-1) + I_R(i,n))*(invmwt(n)*mwmix/rho - HK(n)/(rho*cpmix*T))
  !      enddo
  !   enddo

  ! end subroutine calc_divu



end module feval
