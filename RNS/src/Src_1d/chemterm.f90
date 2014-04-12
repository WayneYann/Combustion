module chemterm_module

  use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UTEMP, UFS, NSPEC, &
       do_cc_burning, split_burning, new_J_cell
  use burner_module, only : burn, compute_rhodYdt, splitburn
  use eos_module, only : eos_get_T
  use weno_module, only : cellavg2gausspt_1d
  use convert_module, only : cellavg2cc_1d, cc2cellavg_1d

  implicit none

  private

  public :: chemterm, dUdt_chem

contains

  subroutine chemterm(lo, hi, U, Ulo, Uhi, dt)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(in) :: dt

    if (split_burning) then
       call chemterm_split(lo, hi, U, Ulo, Uhi, dt)
    else if (do_cc_burning) then
       call chemterm_cellcenter(lo, hi, U, Ulo, Uhi, dt)
    else
       call chemterm_gauss(lo, hi, U, Ulo, Uhi, dt)
    end if

  end subroutine chemterm

  subroutine chemterm_gauss(lo, hi, U, Ulo, Uhi, dt)
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

    force_new_J = .true.  ! always recompute Jacobian when a new FAB starts

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

       force_new_J = new_J_cell

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

  end subroutine chemterm_gauss


  subroutine chemterm_cellcenter(lo, hi, U, Ulo, Uhi, dt)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(in) :: dt

    integer :: i, n
    logical :: force_new_J
    double precision :: rhot(1), rhoinv, ei
    double precision :: Yt0(nspec), Yt(nspec+1)
    double precision, allocatable :: Ucc(:,:)

    allocate(Ucc(lo(1)-1:hi(1)+1,NVAR))

    do n=1,NVAR
       call cellavg2cc_1d(lo(1)-1,hi(1)+1, U(:,n), Ulo(1),Uhi(1), Ucc(:,n), lo(1)-1,hi(1)+1)
    end do

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    do i=lo(1)-1,hi(1)+1

       rhot = 0.d0
       do n=1,NSPEC
          Yt(n) = Ucc(i,UFS+n-1)
          rhot(1) = rhot(1) + Yt(n)
       end do
       rhoinv = 1.d0/rhot(1)
       
       Yt(1:nspec) = Yt(1:nspec) * rhoinv
       Yt(nspec+1) = Ucc(i,UTEMP)

       ei = rhoinv*( Ucc(i,UEDEN) - 0.5d0*rhoinv*Ucc(i,UMX)**2 )

       call eos_get_T(Yt(nspec+1), ei, Yt(1:nspec))

       Yt0 = Yt(1:nspec)
       call burn(1, rhot, Yt, dt, force_new_J)

       force_new_J = new_J_cell

       do n=1,nspec
          Ucc(i,UFS+n-1) = rhot(1)*(Yt(n)-Yt0(n))
       end do

    end do

    do n=UFS,UFS+nspec-1
       call cc2cellavg_1d(lo(1),hi(1), Ucc(:,n), lo(1)-1,hi(1)+1, Ucc(:,UTEMP), lo(1)-1,hi(1)+1)
       do i=lo(1),hi(1)
          U(i,n) = U(i,n) + Ucc(i,UTEMP)
       end do
    end do

    deallocate(Ucc)

  end subroutine chemterm_cellcenter


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


  subroutine chemterm_split(lo, hi, U, Ulo, Uhi, dt)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(in) :: dt

    integer :: i, n, g
    logical :: force_new_J
    double precision :: rhot(2), rhoinv, ei, rho0(1)
    double precision :: Yt(nspec+1,2), Y0(nspec+1)
    double precision, allocatable :: UG(:,:,:)

    allocate(UG(lo(1):hi(1),NVAR,2))

    do n=1,NVAR
       call cellavg2gausspt_1d(lo(1),hi(1), U(:,n), Ulo(1),Uhi(1), UG(:,n,1), UG(:,n,2), lo(1),hi(1))
    end do

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    do i=lo(1),hi(1)
       
       Y0 = 0.d0
       rho0(1) = 0.d0
       
       do g=1,2
          
          rhot(g) = 0.d0
          do n=1,NSPEC
             Yt(n,g) = UG(i,g,UFS+n-1)
             rhot(g) = rhot(g) + Yt(n,g)
          end do
          rhoinv = 1.d0/rhot(g)

          Yt(1:nspec,g) = Yt(1:nspec,g) * rhoinv
          Yt(nspec+1,g) = UG(i,g,UTEMP)

          ei = rhoinv*( UG(i,g,UEDEN) - 0.5d0*rhoinv*UG(i,g,UMX)**2 )
          
          call eos_get_T(Yt(nspec+1,g), ei, Yt(1:nspec,g))
          
          Y0 = Y0 + 0.5d0*Yt(:,g)
          rho0(1) = rho0(1) + 0.5d0*rhot(g)
          
       end do
       
       call burn(1, rho0(1), Y0, dt, force_new_J)
       
       force_new_J = new_J_cell
       
       call splitburn(2, rho0(1), Y0, rhot, Yt, dt)
       
       do g=1,2
          do n=1,nspec
             U(i,UFS+n-1) = U(i,UFS+n-1) + 0.5d0*rhot(g)*Yt(n,g)
          end do
       end do
       
    end do

    deallocate(UG)

  end subroutine chemterm_split

end module chemterm_module
