module feval
  use bdf
  implicit none
  double precision :: rho(8)
  !$omp threadprivate(rho)
contains

  subroutine f_rhs(neq, npt, y, t, ydot)
    use chemistry_module, only : molecular_weight, nspecies
    integer,       intent(in)  :: neq, npt
    real(dp),      intent(in)  :: y(neq,npt), t
    real(dp),      intent(out) :: ydot(neq,npt)

    integer :: iwrk, i, n
    double precision :: rwrk, cv, u(nspecies), Tdot, rYdot, rhoInv
    double precision :: Temp(npt), Ytmp(npt,nspecies), Ydottmp(npt,nspecies)

    do i=1,npt
       do n=1,nspecies
          Ytmp(i,n) = y(n,i)
       end do
       Temp(i) = y(neq,i)
    end do

    call vckwyr(npt, rho, Temp, Ytmp, iwrk, rwrk, ydottmp)

    do i=1,npt

       call ckcvbs(Temp(i), y(1,i), iwrk, rwrk, cv)

       call ckums(Temp(i), iwrk, rwrk, u)

       rhoInv = 1.d0/rho(i)
       Tdot = 0.d0
       do n=1,nspecies
          rYdot = ydottmp(i,n) * molecular_weight(n)
          Tdot = Tdot + u(n)* rYdot
          ydot(n,i) = rYdot * rhoInv
       end do

       ydot(neq,i) = -Tdot/(rho(i)*cv)

    end do

  end subroutine f_rhs

  subroutine f_jac(neq, npt, y, t, pd)
    use chemistry_module, only : molecular_weight, inv_mwt, nspecies
    integer,       intent(in)  :: neq, npt
    real(dp),      intent(in)  :: y(neq,npt), t
    real(dp),      intent(out) :: pd(neq,neq)

    ! local variables
    integer :: iwrk, i, j
    double precision :: rwrk, rho0, rhoinv, Temp, nptinv
    double precision :: YT0(neq), C(neq-1)
    integer, parameter :: consP = 0

    nptinv = 1.d0/dble(npt)

    rho0 = 0.d0
    do i=1, npt
       rho0 = rho0 + rho(i)*nptinv
    end do
    
    Yt0 = 0.d0
    do i=1, npt
       do j=1,neq
          YT0(j) = YT0(j) + y(j,i)*nptinv
       end do
    end do

    Temp = YT0(neq)

    call ckytcr(rho0, Temp, YT0, iwrk, rwrk, C)
    call DWDOT(PD, C, Temp, consP)
    
    do j=1,neq-1
       do i=1,neq-1
          pd(i,j) = pd(i,j) * molecular_weight(i) * inv_mwt(j)
       end do
       i=neq
       pd(i,j) = pd(i,j) * inv_mwt(j) * rho0
    end do
    
    j = neq
    rhoinv = 1.d0/rho0    
    do i=1,neq-1
       pd(i,j) = pd(i,j) * molecular_weight(i) * rhoinv
    enddo
    
  end subroutine f_jac

end module feval
