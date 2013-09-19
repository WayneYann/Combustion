module variables_module

  use meth_params_module
  use eos_module, only : eos_get_T, eos_get_p

  implicit none

contains
  
  subroutine ctoprim(lo, hi, Q, Qlo, Qhi, QVAR)
    integer, intent(in) :: lo(3), hi(3), Qlo(3), Qhi(3), QVAR
    double precision, intent(inout) :: Q(Qlo(1):Qhi(1),Qlo(2):Qhi(2),Qlo(3):Qhi(3),QVAR)

    integer :: i,j,k,n
    double precision :: rhoInv, ek, ei, Yt(NSPEC)

    if (QVAR .ge. QFX) then
       print *, 'variables_module: QVAR > QFX'
       stop
    end if

    do       k = lo(3),hi(3)
       do    j = lo(2),hi(2)
          do i = lo(1),hi(1)

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

             do n=1,nspec
                Q(i,j,k,QFY+n-1) = Q(i,j,k,QFY+n-1)*rhoInv
                Yt(n) = Q(i,j,k,QFY+n-1)
             end do

             call eos_get_T(Q(i,j,k,QTEMP), ei, Yt)

             call eos_get_p(Q(i,j,k,qpres), Q(i,j,k,QRHO), Q(i,j,k,QTEMP), Yt)

          end do
       end do
    end do
  end subroutine ctoprim

end module variables_module
