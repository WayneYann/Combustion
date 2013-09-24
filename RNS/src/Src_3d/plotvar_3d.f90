! xxxxxxxxxxxxxxxx todo

subroutine rns_ctoprim(lo, hi, &
     cons, c_l1, c_l2, c_l3, c_h1, c_h2, c_h3, &
     prim, p_l1, p_l2, p_l3, p_h1, p_h2, p_h3 )
  use meth_params_module, only : NVAR, xblksize, yblksize, zblksize
!xxx  use convert_2d_module, only : cellavg2cc_2d
  use variables_module, only : ctoprim
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: c_l1, c_h1, c_l2, c_h2, c_l3, c_h3
  integer, intent(in) :: p_l1, p_h1, p_l2, p_h2, p_l3, p_h3
  double precision,intent(in )::cons(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,NVAR) 
  double precision,intent(out)::prim(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3,NVAR) 

  integer :: n, clo(3), chi(3), plo(3), phi(3), tlo(3), thi(3)
  integer :: ib, jb, kb, nb(3), blocksize(3)


end subroutine rns_ctoprim

