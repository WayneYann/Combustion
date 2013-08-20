module chemterm_module

  use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UMY, UTEMP, UFS, NSPEC
  use burner_module, only : burn
  use eos_module, only : eos_get_T
  use weno_module, only : cellavg2gausspt_2d

  implicit none

  private

  public :: chemterm

contains

  subroutine chemterm(lo, hi, U, Ulo, Uhi, dt)
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),NVAR)
    double precision, intent(in) :: dt

    integer :: i, j, n, g
    double precision :: rhot, rhoinv, ei
    double precision :: Yt(nspec+1)
    double precision, allocatable :: UG(:,:,:,:)

    allocate(UG(lo(1):hi(1),lo(2):hi(2),4,NVAR))

    do n=1,NVAR
       call cellavg2gausspt_2d(U(:,:,n), Ulo, Uhi, UG(:,:,:,n), lo, hi)
    end do

    call setfirst(.true.)

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          do n=1,nspec
             U(i,j,UFS+n-1) = 0.d0
          end do

          do g=1,4

             rhot = 0.d0
             do n=1,NSPEC
                Yt(n) = UG(i,j,g,UFS+n-1)
                rhot = rhot + Yt(n)
             end do
             rhoinv = 1.d0/rhot

             Yt(1:nspec) = Yt(1:nspec) * rhoinv
             Yt(nspec+1) = UG(i,j,g,UTEMP)

             ei = rhoinv*( UG(i,j,g,UEDEN) - 0.5d0*rhoinv*(UG(i,j,g,UMX)**2 &
                  + UG(i,j,g,UMY)**2) )

             call eos_get_T(Yt(nspec+1), ei, Yt(1:nspec))
       
             call burn(rhot, Yt, dt)

             do n=1,nspec
                U(i,j,UFS+n-1) = U(i,j,UFS+n-1) + 0.25d0 * rhot*Yt(n)
             end do

          end do
       end do
    end do

    deallocate(UG)

  end subroutine chemterm

end module chemterm_module
