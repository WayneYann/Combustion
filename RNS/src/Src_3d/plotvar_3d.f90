subroutine rns_ctoprim(lo, hi, &
     cons, c_l1, c_l2, c_l3, c_h1, c_h2, c_h3, &
     prim, p_l1, p_l2, p_l3, p_h1, p_h2, p_h3 )
  use meth_params_module, only : NVAR, xblksize, yblksize, zblksize
  use convert_3d_module, only : cellavg2cc_3d
  use variables_module, only : ctoprim
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: c_l1, c_h1, c_l2, c_h2, c_l3, c_h3
  integer, intent(in) :: p_l1, p_h1, p_l2, p_h2, p_l3, p_h3
  double precision,intent(in )::cons(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,NVAR) 
  double precision,intent(out)::prim(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3,NVAR) 

  integer :: n, clo(3), chi(3), plo(3), phi(3), tlo(3), thi(3)
  integer :: ib, jb, kb, nb(3), blocksize(3)

  clo(1) = c_l1
  clo(2) = c_l2
  clo(3) = c_l3
  chi(1) = c_h1
  chi(2) = c_h2
  chi(3) = c_h3

  plo(1) = p_l1
  plo(2) = p_l2
  plo(3) = p_l3
  phi(1) = p_h1
  phi(2) = p_h2
  phi(3) = p_h3

  !$omp parallel do private(n)
  do n=1,NVAR
     call cellavg2cc_3d(lo, hi, cons(:,:,:,n), clo, chi, prim(:,:,:,n), plo, phi)
  end do
  !$omp end parallel do


  blocksize(1) = xblksize
  blocksize(2) = yblksize
  blocksize(3) = zblksize

  nb = (hi-lo+blocksize)/blocksize

  !$omp parallel do private(ib,jb,kb,tlo,thi) collapse(3)
  do kb=0,nb(3)-1
     do jb=0,nb(2)-1
        do ib=0,nb(1)-1

           tlo(1) = lo(1) + ib*blocksize(1)
           tlo(2) = lo(2) + jb*blocksize(2)
           tlo(3) = lo(3) + kb*blocksize(3)

           thi(1) = min(tlo(1)+blocksize(1)-1, hi(1))
           thi(2) = min(tlo(2)+blocksize(2)-1, hi(2))
           thi(3) = min(tlo(3)+blocksize(3)-1, hi(3))

           call ctoprim(tlo, thi, prim, plo, phi, NVAR)
        end do
     end do
  end do
  !$omp end parallel do

end subroutine rns_ctoprim

