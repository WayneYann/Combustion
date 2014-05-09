subroutine rns_ctoprim(lo, hi, &
     cons, c_l1, c_l2, c_l3, c_h1, c_h2, c_h3, &
     prim, p_l1, p_l2, p_l3, p_h1, p_h2, p_h3 )
  use meth_params_module, only : NVAR, xblksize, yblksize, zblksize, nthreads
  use convert_module, only : cellavg2cc_3d
  use variables_module, only : ctoprim
  use threadbox_module, only : build_threadbox_3d
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: c_l1, c_h1, c_l2, c_h2, c_l3, c_h3
  integer, intent(in) :: p_l1, p_h1, p_l2, p_h2, p_l3, p_h3
  double precision,intent(in )::cons(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,NVAR) 
  double precision,intent(out)::prim(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3,NVAR) 

  integer :: n, clo(3), chi(3), plo(3), phi(3), tlo(3), thi(3)
  integer :: ib, jb, kb, nb(3), blocksize(3), boxsize(3), nleft(3)
  integer, parameter :: blocksize_min = 4

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

  boxsize = hi-lo+1

  if (nthreads > 1) then
     call build_threadbox_3d(nthreads, boxsize, blocksize_min, nb)
     if (nb(1).eq.0) then
        nb = boxsize/blocksize_min
     end if
     blocksize = boxsize/nb
  else
     blocksize(1) = xblksize 
     blocksize(2) = yblksize
     blocksize(3) = zblksize
     nb = boxsize/blocksize
  end if

  nleft = boxsize - blocksize*nb

  !$omp parallel private(n,ib,jb,kb,tlo,thi)
  !$omp do
  do n=1,NVAR
     call cellavg2cc_3d(lo, hi, cons(:,:,:,n), clo, chi, prim(:,:,:,n), plo, phi)
  end do
  !$omp end do

  !$omp do collapse(3)
  do kb=0,nb(3)-1
     do jb=0,nb(2)-1
        do ib=0,nb(1)-1

           tlo(1) = lo(1) + ib*blocksize(1) + min(nleft(1),ib)
           tlo(2) = lo(2) + jb*blocksize(2) + min(nleft(2),jb)
           tlo(3) = lo(3) + kb*blocksize(3) + min(nleft(3),kb)
           
           thi(1) = tlo(1)+blocksize(1)-1
           thi(2) = tlo(2)+blocksize(2)-1
           thi(3) = tlo(3)+blocksize(3)-1
           
           if (ib < nleft(1)) thi(1) = thi(1) + 1
           if (jb < nleft(2)) thi(2) = thi(2) + 1
           if (kb < nleft(3)) thi(3) = thi(3) + 1

           thi(1) = min(hi(1), thi(1))
           thi(2) = min(hi(2), thi(2))
           thi(3) = min(hi(3), thi(3))

           call ctoprim(tlo, thi, prim, plo, phi, NVAR)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine rns_ctoprim

