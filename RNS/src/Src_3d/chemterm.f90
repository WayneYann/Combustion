module chemterm_module

  use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UMY, UMZ, UTEMP, UFS, NSPEC, &
       do_cc_burning
  use burner_module, only : burn
  use eos_module, only : eos_get_T

  implicit none

  private

  public :: chemterm

contains

  subroutine chemterm(lo, hi, U, Ulo, Uhi, dt)
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),Ulo(3):Uhi(3),NVAR)
    double precision, intent(in) :: dt

    if (do_cc_burning) then
       call chemterm_cellcenter(lo, hi, U, Ulo, Uhi, dt)
    else
       call chemterm_gauss(lo, hi, U, Ulo, Uhi, dt)
    end if

  end subroutine chemterm

  subroutine chemterm_gauss(lo, hi, U, Ulo, Uhi, dt)
    use weno_module, only : cellavg2gausspt_3d
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),Ulo(3):Uhi(3),NVAR)
    double precision, intent(in) :: dt

    integer :: i, j, k, n, g
    logical :: force_new_J
    double precision :: rhot(8), rhoinv, ei
    double precision :: Yt(nspec+1,8)
    double precision, allocatable :: UG(:,:,:,:,:)

    allocate(UG(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),8,NVAR))

    !$omp parallel private(i,j,k,n,g,rhot,rhoinv,ei,Yt,force_new_J)

    !$omp do
    do n=1,NVAR
       call cellavg2gausspt_3d(lo,hi, U(:,:,:,n), Ulo,Uhi, UG(:,:,:,:,n), lo,hi)
    end do
    !$omp end do

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    !$omp do collapse(2)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             do g=1,8

                rhot = 0.d0
                do n=1,NSPEC
                   Yt(n,g) = UG(i,j,k,g,UFS+n-1)
                   rhot(g) = rhot(g) + Yt(n,g)
                end do
                rhoinv = 1.d0/rhot(g)
                
                Yt(1:nspec,g) = Yt(1:nspec,g) * rhoinv
                Yt(nspec+1,g) = UG(i,j,k,g,UTEMP)

                ei = rhoinv*( UG(i,j,k,g,UEDEN) - 0.5d0*rhoinv*(UG(i,j,k,g,UMX)**2 &
                     + UG(i,j,k,g,UMY)**2 + UG(i,j,k,g,UMZ)**2) )

                call eos_get_T(Yt(nspec+1,g), ei, Yt(1:nspec,g))
       
             end do

             call burn(8, rhot, Yt, dt, force_new_J)

             force_new_J = .false.

             do n =1,nspec
                U(i,j,k,UFS+n-1) = 0.d0
                do g=1,8
                   U(i,j,k,UFS+n-1) = U(i,j,k,UFS+n-1) + rhot(g)*Yt(n,g)
                end do
                U(i,j,k,UFS+n-1) = U(i,j,k,UFS+n-1) * 0.125d0
             end do

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(UG)

  end subroutine chemterm_gauss


  subroutine chemterm_cellcenter(lo, hi, U, Ulo, Uhi, dt)
    use convert_3d_module, only : cellavg2cc_3d, cc2cellavg_3d
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),Ulo(3):Uhi(3),NVAR)
    double precision, intent(in) :: dt

    integer :: i, j, k, n
    logical :: force_new_J
    double precision :: rhot(1), rhoinv, ei, fac
    double precision :: Yt(nspec+1)
    double precision, allocatable :: Ucc(:,:,:,:)

    allocate(Ucc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,NVAR))

    !$omp parallel private(i,j,k,n,rhot,rhoinv,ei,fac,Yt,force_new_J)

    !$omp do
    do n=1,NVAR
       call cellavg2cc_3d(lo-1,hi+1, U(:,:,:,n), Ulo,Uhi, Ucc(:,:,:,n), lo-1,hi+1)
    end do
    !$omp end do

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    !$omp do collapse(2)
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1

             rhot = 0.d0
             do n=1,NSPEC
                Yt(n) = Ucc(i,j,k,UFS+n-1)
                rhot(1) = rhot(1) + Yt(n)
             end do
             rhoinv = 1.d0/rhot(1)
                
             Yt(1:nspec) = Yt(1:nspec) * rhoinv
             Yt(nspec+1) = Ucc(i,j,k,UTEMP)

             ei = rhoinv*( Ucc(i,j,k,UEDEN) - 0.5d0*rhoinv*(Ucc(i,j,k,UMX)**2 &
                  + Ucc(i,j,k,UMY)**2 + Ucc(i,j,k,UMZ)**2) )

             call eos_get_T(Yt(nspec+1), ei, Yt(1:nspec))
       
             call burn(1, rhot, Yt, dt, force_new_J)

             force_new_J = .false.

             do n=1,nspec
                Ucc(i,j,k,UFS+n-1) = rhot(1)*Yt(n)
             end do
             U(i,j,k,UTEMP) = Yt(nspec+1)

          end do
       end do
    end do
    !$omp end do

    !$omp do
    do n=UFS,UFS+nspec-1
       call cc2cellavg_3d(lo,hi, Ucc(:,:,:,n), lo-1,hi+1, U(:,:,:,n), Ulo,Uhi)
    end do
    !$omp end do

    !$omp do collapse(2)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhot = 0.d0
             do n=1,NSPEC
                rhot(1) = rhot(1) + U(i,j,k,UFS+n-1)
             end do
             fac = U(i,j,k,URHO)/rhot(1)
             do n=1,NSPEC
                U(i,j,k,UFS+n-1) = U(i,j,k,UFS+n-1) * fac
             end do
          end do
       end do
    end do
    !$omp end do

    !$omp end parallel

    deallocate(Ucc)

  end subroutine chemterm_cellcenter

end module chemterm_module
