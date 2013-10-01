!
! Method-of-lines discretization for LMC.
!
! 
!

module feval
  use iso_c_binding
  use encap
  implicit none
  
contains


  !
  ! Compute MOL dU/dt for the advection piece of LMC.
  !
  subroutine advection(vel_old, scal_old, divu_old, press_old, dx, lo, hi, bc)
    use spec_module

    real*8  dx(0:nlevs-1)
    integer lo(0:nlevs-1)
    integer hi(0:nlevs-1)
    real*8  press_old(0:nlevs-1,-1:nfine+1)
    real*8    vel_old(0:nlevs-1,-2:nfine+1)
    real*8   scal_old(0:nlevs-1,-2:nfine+1,nscal)
    real*8   divu_old(0:nlevs-1,-1:nfine)
    integer bc(0:nlevs-1,2)

    real*8  dt(0:nlevs-1)
    real*8           gp(0:nlevs-1,-1:nfine)
    real*8       macvel(0:nlevs-1, 0:nfine  )
    real*8     diff_old(0:nlevs-1,-1:nfine,  nscal)
    ! real*8     diff_new(0:nlevs-1,-1:nfine,  nscal)
    ! real*8     diff_hat(0:nlevs-1,-1:nfine,  nscal)
    ! real*8     diff_tmp(0:nlevs-1,-1:nfine,  nscal)
    real*8       tforce(0:nlevs-1,-1:nfine,  nscal)
    real*8 diffdiff_old(0:nlevs-1,-1:nfine)

    real*8          aofs(0:nlevs-1, 0:nfine-1,nscal)

    integer i, is, n

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


    print *,'... computing advective scalar update'

    ! compute advective flux divergence
    tforce = 0
    call scal_aofs(scal_old(0,:,:),macvel(0,:),aofs(0,:,:), &
                       divu_old(0,:),tforce(0,:,:),dx(0),dt(0), &
                       lo(0),hi(0),bc(0,:))

  end subroutine advection

  subroutine f1eval(fptr, qptr, t, state, ctx) bind(c)
    use spec_module
    type(c_ptr),    intent(in), value :: fptr, qptr, state, ctx
    real(c_double), intent(in), value :: t

    type(lmc_encap), pointer :: q, f
    integer :: lo(0:nlevs-1), hi(0:nlevs-1)
    
    call c_f_pointer(qptr, q)
    call c_f_pointer(fptr, f)

    lo = 0
    hi = nx-1

    ! call advection()

    ! XXX: set_bc_s, set_bc_v

  end subroutine f1eval

  subroutine f2eval(fptr, qptr, t, state, ctx) bind(c)
    type(c_ptr),    intent(in), value :: fptr, qptr, state, ctx
    real(c_double), intent(in), value :: t

    ! print *, 'F2EVAL'
  end subroutine f2eval

  subroutine f2comp(fptr, qptr, t, dt, rhsptr, state, ctx) bind(c)
    type(c_ptr),    intent(in), value :: fptr, qptr, rhsptr, state, ctx
    real(c_double), intent(in), value :: t, dt

    ! print *, 'F2COMP'
  end subroutine f2comp

end module feval
