module difterm_module

  implicit none

  private

  public :: difterm

contains

  subroutine difterm(lo,hi,U,Ulo,Uhi,flx, dxinv)

    use meth_params_module, only : NVAR, NSPEC, QCVAR, QFVAR
    use weno_module, only : cellavg2face_1d
    use convert_module, only : cellavg2cc_1d
    use variables_module, only : ctoprim
    use transport_properties, only : get_transport_properties

    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(in ) ::   U(Ulo(1):Uhi(1)  ,NVAR)
    double precision, intent(out) :: flx( lo(1): hi(1)+1,NVAR)
    double precision, intent(in) :: dxinv(1)

    double precision, allocatable :: Qc(:,:), Qf(:,:)  ! cell-center and face
    double precision, allocatable :: mu(:), xi(:), lam(:), Ddia(:,:)
    integer :: Qclo(3), Qchi(3)
    integer :: Qflo(3), Qfhi(3)
    integer :: n

    Qclo = 1;  Qchi = 1;  Qflo = 1;  Qfhi = 1

    Qclo(1) = lo(1)-2   ! need 2 ghost cells for fourth-order derivatives on face
    Qchi(1) = hi(1)+2

    Qflo(1) = lo(1)
    Qfhi(1) = hi(1)+1

    allocate(  Qc(Qclo(1):Qchi(1), QCVAR))
    allocate(  Qf(Qflo(1):Qfhi(1), QFVAR))
    allocate(  mu(Qflo(1):Qfhi(1)))
    allocate(  xi(Qflo(1):Qfhi(1)))
    allocate( lam(Qflo(1):Qfhi(1)))
    allocate(Ddia(Qflo(1):Qfhi(1), NSPEC))
 
    ! cell-centered variables
    do n=1,NVAR
       call cellavg2cc_1d(U(:,n), Ulo(1), Uhi(1), Qc(:,n), Qclo(1), Qchi(1))
    end do
    !
    call ctoprim(Qclo,Qchi, Qc, Qclo,Qchi,QCVAR)

    ! face variables
    do n=1,NVAR
       call cellavg2face_1d(Qflo(1),Qfhi(1), U(:,n), Ulo(1),Uhi(1), Qf(:,n), Qflo(1), Qfhi(1))
    end do
    !
    call ctoprim(Qflo,Qfhi, Qf, Qflo,Qfhi,QFVAR)

    ! transport coefficients on face
    call get_transport_properties(Qflo,Qfhi, Qf, Qflo,Qfhi,QFVAR, mu, xi, lam, Ddia, Qflo,Qfhi)

    call comp_diff_flux(flx, Qf, mu, xi, lam, Ddia, Qflo, Qfhi, Qc, Qclo, Qchi, dxinv)

    deallocate(Qc,Qf)
    deallocate(mu,xi,lam,Ddia)

  end subroutine difterm


  subroutine comp_diff_flux(flx, Qf, mu, xi, lam, Ddia, Qflo, Qfhi, Qc, Qclo, Qchi, dxinv)

    use meth_params_module
    use derivative_stencil_module, only : FD4

    integer, intent(in) :: Qflo(1), Qfhi(1), Qclo(1), Qchi(1)
    double precision, intent(in) :: dxinv(1)
    double precision, intent(in)  ::   Qc(Qclo(1):Qchi(1),QCVAR)
    double precision, intent(in)  ::   Qf(Qflo(1):Qfhi(1),QFVAR)
    double precision, intent(in)  ::   mu(Qflo(1):Qfhi(1))
    double precision, intent(in)  ::   xi(Qflo(1):Qfhi(1))
    double precision, intent(in)  ::  lam(Qflo(1):Qfhi(1))
    double precision, intent(in)  :: Ddia(Qflo(1):Qfhi(1),NSPEC)
    double precision, intent(out) ::  flx(Qflo(1):Qfhi(1),NVAR)

    integer :: i, n, UYN, QYN, QXN, QHN
    double precision :: tauxx, dudx, dTdx, dXdx, Vd
    double precision, dimension(Qflo(1):Qfhi(1)) :: dlnpdx, rVc
    double precision, parameter :: fourThirds = 4.d0/3.d0

    flx(:,URHO) = 0.d0
    flx(:,UTEMP) = 0.d0

    do i = Qflo(1), Qfhi(1)

       ! viscous stress 
       dudx = dxinv(1) * (FD4(-2)*Qc(i-2,QU) + FD4(-1)*Qc(i-1,QU) &
            + FD4(0)*Qc(i,QU) + FD4(1)*Qc(i+1,QU))
       tauxx = (fourThirds*mu(i) + xi(i)) * dudx
       flx(i,UMX) = -tauxx
       flx(i,UEDEN) = -tauxx * Qf(i,QU)

       ! thermal conduction
       dTdx = dxinv(1) * (FD4(-2)*Qc(i-2,QTEMP) + FD4(-1)*Qc(i-1,QTEMP) &
            + FD4(0)*Qc(i,QTEMP) + FD4(1)*Qc(i+1,QTEMP))
       flx(i,UEDEN) = flx(i,UEDEN) - lam(i)*dTdx

       ! compute dpdx
       dlnpdx(i) = dxinv(1) * (FD4(-2)*Qc(i-2,QPRES) + FD4(-1)*Qc(i-1,QPRES) &
            + FD4(0)*Qc(i,QPRES) + FD4(1)*Qc(i+1,QPRES)) / Qf(i,QPRES)

    end do

    rVc = 0.d0

    do n=1,NSPEC

       UYN = UFS+n-1
       QYN = QFY+n-1
       QXN = QFX+n-1

       do i = Qflo(1), Qfhi(1)

          dXdx = dxinv(1) * (FD4(-2)*Qc(i-2,QXN) + FD4(-1)*Qc(i-1,QXN) &
            + FD4(0)*Qc(i,QXN) + FD4(1)*Qc(i+1,QXN))

          Vd = -Ddia(i,n)*(dXdx + (Qf(i,QXN)-Qf(i,QYN))*dlnpdx(i))
          
          flx(i,UYN) = Qf(i,QRHO)*Qf(i,QYN)*Vd
          rVc(i) = rVc(i) + flx(i,UYN)
       end do
    end do

    do n=1,NSPEC

       UYN = UFS+n-1
       QYN = QFY+n-1
       QHN = QFH+n-1

       do i = Qflo(1), Qfhi(1)
          flx(i,UYN) = flx(i,UYN) - Qf(i,QYN)*rVc(i)
          flx(i,UEDEN) = flx(i,UEDEN) + flx(i,UYN)*Qf(i,QHN)
       end do
    end do

  end subroutine comp_diff_flux

end module difterm_module
