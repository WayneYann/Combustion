module feval
  use bdf
  implicit none
  double precision :: rho
  !$omp threadprivate(rho)
contains

  subroutine f_rhs(neq, y, t, ydot)
    use chemistry_module, only : molecular_weight, nspecies
    integer,       intent(in)  :: neq
    real(dp),      intent(in)  :: y(neq), t
    real(dp),      intent(out) :: ydot(neq)

    integer :: iwrk, i
    double precision :: rwrk, Temp, cv, u(nspecies), Tdot, rYdot, rhoInv

    Temp = y(neq)

    call ckwyr(rho, Temp, y, iwrk, rwrk, ydot)

    call ckcvbs(Temp, y, iwrk, rwrk, cv)

    call ckums(Temp, iwrk, rwrk, u)

    rhoInv = 1.d0/rho
    Tdot = 0.d0
    do i=1,nspecies
       rYdot = ydot(i) * molecular_weight(i)
       Tdot = Tdot + u(i)* rYdot
       ydot(i) = rYdot * rhoInv
    end do

    ydot(neq) = -Tdot/(rho*cv)

  end subroutine f_rhs

  subroutine f_jac(neq, y, t, pd)
    use chemistry_module, only : molecular_weight, inv_mwt, nspecies
    integer,       intent(in)  :: neq
    real(dp),      intent(in)  :: y(neq), t
    real(dp),      intent(out) :: pd(neq,neq)

    ! local variables
    integer :: iwrk, i, j
    double precision :: rwrk, rhoinv, Temp, C(neq-1)
    integer, parameter :: consP = 0

    Temp = y(neq)

    rhoinv = 1.d0/rho
    
    call ckytcr(rho, Temp, y, iwrk, rwrk, C)
    call DWDOT(PD, C, Temp, consP)
    
    do j=1,neq-1
       do i=1,neq-1
          pd(i,j) = pd(i,j) * molecular_weight(i) * inv_mwt(j)
       end do
       i=neq
       pd(i,j) = pd(i,j) * inv_mwt(j) * rho
    end do
    
    j = neq
    do i=1,neq-1
       pd(i,j) = pd(i,j) * molecular_weight(i) * rhoinv
    enddo
    
  end subroutine f_jac

end module feval
