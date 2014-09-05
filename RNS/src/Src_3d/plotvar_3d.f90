subroutine rns_ctoprim(lo, hi, &
     cons, c_l1, c_l2, c_l3, c_h1, c_h2, c_h3, &
     prim, p_l1, p_l2, p_l3, p_h1, p_h2, p_h3 )
  use meth_params_module, only : NVAR, xblksize, yblksize, zblksize, nthreads
  use convert_module, only : cellavg2cc_3d
  use variables_module, only : ctoprim
  use threadbox_module, only : build_threadbox_3d, get_lo_hi
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: c_l1, c_h1, c_l2, c_h2, c_l3, c_h3
  integer, intent(in) :: p_l1, p_h1, p_l2, p_h2, p_l3, p_h3
  double precision,intent(in )::cons(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,NVAR) 
  double precision,intent(out)::prim(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3,NVAR) 

  integer :: n, clo(3), chi(3), plo(3), phi(3), tlo(3), thi(3)
  integer :: ib, jb, kb, nb(3), boxsize(3)
  integer, allocatable :: bxlo(:), bxhi(:), bylo(:), byhi(:), bzlo(:), bzhi(:)
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
  else
     nb(1) = max(boxsize(1)/xblksize, 1)
     nb(2) = max(boxsize(2)/yblksize, 1)
     nb(3) = max(boxsize(3)/zblksize, 1)
  end if

  allocate(bxlo(0:nb(1)-1))
  allocate(bxhi(0:nb(1)-1))
  allocate(bylo(0:nb(2)-1))
  allocate(byhi(0:nb(2)-1))
  allocate(bzlo(0:nb(3)-1))
  allocate(bzhi(0:nb(3)-1))

  call get_lo_hi(boxsize(1), nb(1), bxlo, bxhi)
  call get_lo_hi(boxsize(2), nb(2), bylo, byhi)
  call get_lo_hi(boxsize(3), nb(3), bzlo, bzhi)

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

           tlo(1) = lo(1) + bxlo(ib)
           thi(1) = lo(1) + bxhi(ib)
           
           tlo(2) = lo(2) + bylo(jb)
           thi(2) = lo(2) + byhi(jb)

           tlo(3) = lo(3) + bzlo(kb)
           thi(3) = lo(3) + bzhi(kb)

           call ctoprim(tlo, thi, prim, plo, phi, NVAR)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine rns_ctoprim

