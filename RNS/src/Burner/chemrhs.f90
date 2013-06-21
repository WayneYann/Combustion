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


subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  integer :: neq
  integer :: ml, mu, nrpd

  double precision :: y(neq)
  double precision :: pd(neq,neq)

  double precision :: rpar
  integer :: ipar

  double precision :: t

  return
end subroutine jac
