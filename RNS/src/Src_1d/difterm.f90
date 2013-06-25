module difterm_module

  implicit none

  private

  public :: difterm

contains

  subroutine difterm(lo,hi,U,Ulo,Uhi,flx)

    use meth_params_module, only : NVAR, NSPEC, QCVAR, QFVAR
    use weno_module, only : cellavg2face_1d
    use convert_module, only : cellavg2cc_1d
    use variables_module, only : ctoprim
    use transport_properties, only : get_transport_properties

    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(in ) ::   U(Ulo(1):Uhi(1)  ,NVAR)
    double precision, intent(out) :: flx( lo(1): hi(1)+1,NVAR)

    double precision, allocatable :: Qc(:,:), Qf(:,:)  ! cell-center and face
    double precision, allocatable :: mu(:), xi(:), lam(:), Ddiag(:,:)
    integer :: Qclo(3), Qchi(3)
    integer :: Qflo(3), Qfhi(3)
    integer :: n

    Qclo = 1;  Qchi = 1;  Qflo = 1;  Qfhi = 1

    Qclo(1) = lo(1)-2   ! need 2 ghost cells for fourth-order derivatives on face
    Qchi(1) = hi(1)+2

    Qflo(1) = lo(1)
    Qfhi(1) = hi(1)+1

    allocate(   Qc(Qclo(1):Qchi(1), QCVAR))
    allocate(   Qf(Qflo(1):Qfhi(1), QFVAR))
    allocate(   mu(Qflo(1):Qfhi(1)))
    allocate(   xi(Qflo(1):Qfhi(1)))
    allocate(  lam(Qflo(1):Qfhi(1)))
    allocate(Ddiag(Qflo(1):Qfhi(1), NSPEC))
 
    ! cell-centered variables
    do n=1,NVAR
       call cellavg2cc_1d(U(:,n), Ulo(1), Uhi(1), Qc(:,n), Qclo(1), Qchi(1))
    end do
    !
    call ctoprim(Qc, Qclo, Qchi, QCVAR)

    ! face variables
    do n=1,NVAR
       call cellavg2face_1d(U(:,n), Ulo(1), Uhi(1), Qf(:,n), Qflo(1), Qfhi(1))
    end do
    !
    call ctoprim(Qf, Qflo, Qfhi, QFVAR)

    ! transport coefficients on face
    call get_transport_properties(Qf, Qflo, Qfhi, QFVAR, mu, xi, lam, Ddiag)

    flx = 0.d0

    deallocate(Qc,Qf)
    deallocate(mu,xi,lam,Ddiag)

  end subroutine difterm

end module difterm_module
