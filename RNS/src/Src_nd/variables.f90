module variables_module

  use meth_params_module
  use eos_module, only : eos_get_T

  implicit none

contains
  
  subroutine ctoprim(lo, hi, Q, Qlo, Qhi, QVAR)
    integer, intent(in) :: lo(3), hi(3), Qlo(3), Qhi(3), QVAR
    double precision, intent(inout) :: Q(Qlo(1):Qhi(1),Qlo(2):Qhi(2),Qlo(3):Qhi(3),QVAR)

    integer :: i,j,k,n,iwrk,pt_index(3),ierr
    double precision :: rwrk, rhoInv, ek, ei, Yt(NSPEC), Xt(NSPEC), ht(NSPEC)

    do       k = lo(3),hi(3)
       do    j = lo(2),hi(2)
          do i = lo(1),hi(1)

             pt_index(1) = i
             pt_index(2) = j
             pt_index(3) = k

             rhoInv = 1.d0/Q(i,j,k,QRHO)

             Q(i,j,k,QU) = Q(i,j,k,QU) * rhoInv
             ek = 0.5d0*Q(i,j,k,QU)**2
             if (ndim .ge. 2) then
                Q(i,j,k,QV) = Q(i,j,k,QV) * rhoInv
                ek = ek + 0.5d0*Q(i,j,k,QV)**2
             end if
             if (ndim .eq. 3) then
                Q(i,j,k,QW) = Q(i,j,k,QW) * rhoInv
                ek = ek + 0.5d0*Q(i,j,k,QW)**2
             end if

             ei = Q(i,j,k,UEDEN)*rhoInv - ek

             rhoInv = 1.d0/sum(Q(i,j,k,QFY:QFY+nspec-1))
             do n=1,nspec
                Q(i,j,k,QFY+n-1) = Q(i,j,k,QFY+n-1)*rhoInv
                Yt(n) = Q(i,j,k,QFY+n-1)
             end do

             call eos_get_T(Q(i,j,k,QTEMP), ei, Yt, pt_index, ierr)

             if (ierr .ne. 0) then
                print *, 'ctoprim failed at ', i,j,k,Q(i,j,k,1:QFY+nspec-1)
                call bl_error("ctoprim failed")
             end if

             call ckpy(Q(i,j,k,QRHO), Q(i,j,k,QTEMP), Yt, iwrk, rwrk, Q(i,j,k,qpres))

             if (QVAR .ge. QFX) then
                call ckytx(Yt, iwrk, rwrk, Xt)
                do n=1,nspec
                   Q(i,j,k,QFX+n-1) = Xt(n)
                end do
             end if

             if (QVAR .ge. QFH) then
                call ckhms(Q(i,j,k,QTEMP), iwrk, rwrk, ht)
                do n=1,nspec
                   Q(i,j,k,QFH+n-1) = ht(n)
                end do
             end if

          end do
       end do
    end do
  end subroutine ctoprim

end module variables_module
