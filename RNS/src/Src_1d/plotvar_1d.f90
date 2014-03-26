
subroutine rns_ctoprim(lo, hi, &
     cons, c_l1, c_h1, &
     prim, p_l1, p_h1 )
  use meth_params_module, only : NVAR
  use convert_module, only : cellavg2cc_1d
  use variables_module, only : ctoprim
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: c_l1, c_h1
  integer, intent(in) :: p_l1, p_h1
  double precision, intent(in ) :: cons(c_l1:c_h1,NVAR) 
  double precision, intent(out) :: prim(p_l1:p_h1,NVAR) 

  integer :: n, plo(3), phi(3), tlo(3), thi(3)

  do n=1,NVAR
     call cellavg2cc_1d(lo(1),hi(1), cons(:,n), c_l1, c_h1, prim(:,n), p_l1, p_h1)
  end do

  tlo(1) = lo(1)
  tlo(2:3) = 1
  thi(1) = hi(1)
  thi(2:3) = 1

  plo(1) = p_l1
  plo(2:3) = 1
  phi(1) = p_h1
  phi(2:3) = 1

  call ctoprim(tlo, thi, prim, plo, phi, NVAR)

end subroutine rns_ctoprim

