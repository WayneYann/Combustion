module chemterm_module

  use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UMY, UTEMP, UFS, NSPEC, &
       do_cc_burning
  use burner_module, only : burn
  use eos_module, only : eos_get_T

  implicit none

  private

  public :: chemterm

contains

  subroutine chemterm(lo, hi, U, Ulo, Uhi, dt)
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),NVAR)
    double precision, intent(in) :: dt

    if (do_cc_burning) then
       call chemterm_cellcenter(lo, hi, U, Ulo, Uhi, dt)
    else
       call chemterm_gauss(lo, hi, U, Ulo, Uhi, dt)
    end if

  end subroutine chemterm

  subroutine chemterm_gauss(lo, hi, U, Ulo, Uhi, dt)
    use weno_module, only : cellavg2gausspt_2d
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),NVAR)
    double precision, intent(in) :: dt

    integer :: i, j, n, g
    double precision :: rhot, rhoinv, ei
    double precision :: Yt(nspec+1)
    double precision, allocatable :: UG(:,:,:,:)

    allocate(UG(lo(1):hi(1),lo(2):hi(2),4,NVAR))

    !$omp parallel private(i,j,n,g,rhot,rhoinv,ei,Yt)

    !$omp do
    do n=1,NVAR
       call cellavg2gausspt_2d(lo,hi, U(:,:,n), Ulo,Uhi, UG(:,:,:,n), lo,hi)
    end do
    !$omp end do

    call setfirst(.true.)

    !$omp do collapse(2)
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
    !$omp end do
    !$omp end parallel

    deallocate(UG)

  end subroutine chemterm_gauss


  subroutine chemterm_cellcenter(lo, hi, U, Ulo, Uhi, dt)
    use convert_2d_module, only : cellavg2cc_2d, cellcenter2ca_2d
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),NVAR)
    double precision, intent(in) :: dt

    integer :: i, j, n
    double precision :: rhot, rhoinv, ei
    double precision :: Yt(nspec+1)
    double precision, allocatable :: Ucc(:,:,:)

    allocate(Ucc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,NVAR))

    !$omp parallel private(i,j,n,rhot,rhoinv,ei,Yt)

    !$omp do
    do n=1,NVAR
       call cellavg2cc_2d(lo-1,hi+1, U(:,:,n), Ulo,Uhi, Ucc(:,:,n), lo-1,hi+1)
    end do
    !$omp end do

    call setfirst(.true.)

    !$omp do collapse(2)
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1

          rhot = 0.d0
          do n=1,NSPEC
             Yt(n) = Ucc(i,j,UFS+n-1)
             rhot = rhot + Yt(n)
          end do
          rhoinv = 1.d0/rhot

          Yt(1:nspec) = Yt(1:nspec) * rhoinv
          Yt(nspec+1) = Ucc(i,j,UTEMP)

          ei = rhoinv*( Ucc(i,j,UEDEN) - 0.5d0*rhoinv*(Ucc(i,j,UMX)**2 &
               + Ucc(i,j,UMY)**2) )

          call eos_get_T(Yt(nspec+1), ei, Yt(1:nspec))
       
          call burn(rhot, Yt, dt)

          do n=1,nspec
             Ucc(i,j,UFS+n-1) = Yt(n)
          end do

       end do
    end do
    !$omp end do

    !$omp do
    do n=UFS,UFS+nspec-1
       call cellcenter2ca_2d(lo,hi, Ucc(:,:,n), lo-1,hi+1, U(:,:,n), Ulo,Uhi)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             U(i,j,n) = U(i,j,n) * U(i,j,URHO)
          end do
       end do
    end do
    !$omp end do

    !$omp end parallel

    deallocate(Ucc)

  end subroutine chemterm_cellcenter

end module chemterm_module
