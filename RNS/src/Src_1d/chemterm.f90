module chemterm_module

  use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UTEMP, UFS, NSPEC
  use burner_module, only : burn, compute_rhodYdt
  use eos_module, only : eos_get_T
  use weno_module, only : cellavg2gausspt_1d

  implicit none

  private

  public :: chemterm, dUdt_chem

contains

  subroutine chemterm(lo, hi, U, Ulo, Uhi, dt)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(in) :: dt

    integer :: i, n, g
    logical :: force_new_J
    double precision :: rho(2), rhoinv, ei
    double precision :: YT0(nspec+1,2), YT(nspec+1,2), dry(nspec)
    double precision, allocatable :: UG(:,:,:)

    allocate(UG(lo(1):hi(1),NVAR,2))

    do n=1,NVAR
       call cellavg2gausspt_1d(lo(1),hi(1), U(:,n), Ulo(1),Uhi(1), UG(:,n,1), UG(:,n,2), lo(1),hi(1))
    end do

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    do i=lo(1),hi(1)

       do g=1,2

          rho(g) = 0.d0
          do n=UFS,UFS+nspec-1
             YT(n-UFS+1,g) = UG(i,n,g)
             rho(g) = rho(g) + UG(i,n,g)
          end do
          rhoinv = 1.d0/rho(g)

          YT(1:nspec,g) = YT(1:nspec,g) * rhoinv
          YT(nspec+1,g) = UG(i,UTEMP,g)

          ei = rhoinv*( UG(i,UEDEN,g) - 0.5d0*rhoinv*UG(i,UMX,g)**2 )

          call eos_get_T(YT(nspec+1,g), ei, YT(1:nspec,g))

       end do

       YT0 = YT
       call burn(2, rho, YT, dt, force_new_J)

       force_new_J = .false.

       dry = 0.d0
       do g=1,2
          do n=1,nspec
             dry(n) = dry(n) + rho(g)*(Yt(n,g)-Yt0(n,g))
          end do
       end do
       
       ! note that the sum of dry is zero
       do n=1,nspec
          U(i,UFS+n-1) = U(i,UFS+n-1) + 0.5d0*dry(n)
       end do

    end do

    deallocate(UG)

  end subroutine chemterm


  subroutine dUdt_chem(lo, hi, U, Ulo, Uhi, Ut, Utlo, Uthi)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1), Utlo(1), Uthi(1)
    double precision, intent(in ) ::  U( Ulo(1): Uhi(1),NVAR)
    double precision, intent(out) :: Ut(Utlo(1):Uthi(1),NVAR)

    integer :: i, n, g, np
    double precision :: rhoinv, ei
    double precision :: rho(lo(1):hi(1)), T(lo(1):hi(1))
    double precision :: Ytmp(nspec)
    double precision :: Y(lo(1):hi(1),nspec), rdYdt(lo(1):hi(1),nspec)
    double precision, allocatable :: UG(:,:,:)

    np = hi(1)-lo(1)+1

    allocate(UG(lo(1):hi(1),NVAR,2))

    do n=1,NVAR
       do i=lo(1),hi(1)
          Ut(i,n) = 0.d0
       end do
    end do

    do n=1,NVAR
       call cellavg2gausspt_1d(lo(1),hi(1), U(:,n), Ulo(1),Uhi(1), UG(:,n,1), UG(:,n,2), lo(1),hi(1))
    end do

    do g=1,2

       do i=lo(1),hi(1)
          rho(i) = 0.d0
          do n=1,nspec
             Y(i,n) = UG(i,UFS+n-1,g)
             rho(i) = rho(i) + Y(i,n)
          end do
          rhoinv = 1.d0/rho(i)

          do n=1,nspec
             Y(i,n) = Y(i,n) * rhoinv
             Ytmp(n) = Y(i,n)
          end do

          ei = rhoinv*( UG(i,UEDEN,g) - 0.5d0*rhoinv*UG(i,UMX,g)**2 )

          T(i) = UG(i,UTEMP,g)
          call eos_get_T(T(i), ei, Ytmp)
       end do

       call compute_rhodYdt(np, rho,T,Y,rdYdt)

       do n=1,nspec
          do i=lo(1),hi(1)
             Ut(i,UFS+n-1) = Ut(i,UFS+n-1) + 0.5d0*rdYdt(i,n)
          end do
       end do

    end do

    deallocate(UG)

  end subroutine dUdt_chem

end module chemterm_module
