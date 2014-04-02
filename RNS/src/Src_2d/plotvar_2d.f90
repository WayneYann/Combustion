
subroutine rns_ctoprim(lo, hi, &
     cons, c_l1, c_l2, c_h1, c_h2, &
     prim, p_l1, p_l2, p_h1, p_h2 )
  use meth_params_module, only : NVAR, xblksize, yblksize
  use convert_module, only : cellavg2cc_2d
  use variables_module, only : ctoprim
  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: c_l1, c_h1, c_l2, c_h2
  integer, intent(in) :: p_l1, p_h1, p_l2, p_h2
  double precision, intent(in ) :: cons(c_l1:c_h1,c_l2:c_h2,NVAR) 
  double precision, intent(out) :: prim(p_l1:p_h1,p_l2:p_h2,NVAR) 

  integer :: n, clo(3), chi(3), plo(3), phi(3), tlo(3), thi(3)
  integer :: ib, jb, nb(2), blocksize(2)

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

  !$omp parallel do private(n)
  do n=1,NVAR
     call cellavg2cc_2d(lo, hi, cons(:,:,n), clo(1:2), chi(1:2), prim(:,:,n), plo(1:2), phi(1:2))
  end do
  !$omp end parallel do

  blocksize(1) = xblksize
  blocksize(2) = yblksize

  nb = (hi-lo+blocksize)/blocksize

  !$omp parallel do private(ib,jb,tlo,thi) collapse(2)
  do jb=0,nb(2)-1
     do ib=0,nb(1)-1

        tlo(1) = lo(1) + ib*blocksize(1)
        tlo(2) = lo(2) + jb*blocksize(2)
        tlo(3) = 1

        thi(1) = min(tlo(1)+blocksize(1)-1, hi(1))
        thi(2) = min(tlo(2)+blocksize(2)-1, hi(2))
        thi(3) = 1

        call ctoprim(tlo, thi, prim, plo, phi, NVAR)
     end do
  end do
  !$omp end parallel do

end subroutine rns_ctoprim

