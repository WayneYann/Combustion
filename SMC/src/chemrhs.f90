subroutine f_rhs(n, t, y, ydot, rpar, ipar)
  use chemistry_module, only : molecular_weight, nspecies
  implicit none

  integer, intent(in) :: n, ipar
  double precision, intent(in) :: t, y(n), rpar(*)
  double precision, intent(out) :: ydot(n)

  integer :: iwrk, i
  double precision :: rwrk, Temp, rho, cv, u(nspecies), Tdot, rYdot, rhoInv

  rho = rpar(1)
  Temp = y(n)

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

  ydot(n) = -Tdot/(rho*cv)

end subroutine f_rhs


subroutine f_jac(neq, time, y, ml, mu, pd, nrpd, rpar, ipar)
  implicit none
  integer :: neq
  integer :: ml, mu, nrpd

  double precision :: y(neq)
  double precision :: pd(nrpd,neq)

  double precision :: rpar(*)
  integer :: ipar(*)

  double precision :: time

  ! local variables
  integer :: iwrk
  double precision :: rwrk, rho, T, C(neq-1)
  integer, parameter :: consP = 0

  rho = rpar(1)
  T = y(neq)

  call ckytcr(rho, T, y, iwrk, rwrk, C)
  call DWDOT(PD, C, T, consP)

  return
end subroutine f_jac
