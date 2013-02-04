subroutine f_rhs(n, t, y, ydot, rpar, ipar)
  use chemistry_module, only : molecular_weight
  implicit none

  integer, intent(in) :: n, ipar
  double precision, intent(in) :: t, y(n), rpar(*)
  double precision, intent(out) :: ydot(n)

  integer :: iwrk
  double precision :: rwrk, Temp, rho, ei

  rho = rpar(1)
  ei  = rpar(2)

  call feeytt(ei, Y, iwrk, rwrk, Temp)

  call ckwyr(rho, Temp, y, iwrk, rwrk, ydot)

  ydot = ydot * molecular_weight / rho

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
