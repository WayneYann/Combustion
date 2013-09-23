
subroutine rns_ctoprim(lo, hi, &
     cons, c_l1, c_l2, c_h1, c_h2, &
     prim, p_l1, p_l2, p_h1, p_h2 )
  use meth_params_module, only : NVAR
  use convert_2d_module, only : cellavg2cc_2d
  use variables_module, only : ctoprim
  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: c_l1, c_h1, c_l2, c_h2
  integer, intent(in) :: p_l1, p_h1, p_l2, p_h2
  double precision, intent(in ) :: cons(c_l1:c_h1,c_l2:c_h2,NVAR) 
  double precision, intent(out) :: prim(p_l1:p_h1,p_l2:p_h2,NVAR) 

  integer :: n, clo(3), chi(3), plo(3), phi(3), tlo(3), thi(3)

  clo(1) = c_l1
  clo(2) = c_l2
  clo(3) = 1
  chi(1) = c_h1
  chi(2) = c_h2
  chi(3) = 1

  plo(1) = p_l1
  plo(2) = p_l2
  plo(3) = 1
  phi(1) = p_h1
  phi(2) = p_h2
  phi(3) = 1

  do n=1,NVAR
     call cellavg2cc_2d(lo, hi, cons(:,:,n), clo(1:2), chi(1:2), prim(:,:,n), plo, phi)
  end do

  tlo(1:2) = lo
  tlo(3) = 1
  thi(1:2) = hi
  thi(3) = 1

  call ctoprim(tlo, thi, prim, plo, phi, NVAR)

end subroutine rns_ctoprim

