module difterm_module

  implicit none

  private

  public :: difterm

contains

  subroutine difterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dxinv)

    use meth_params_module, only : NVAR, NSPEC, QCVAR, QFVAR
    use weno_module, only : cellavg2gausspt_1d, cellavg2face_1d, cellavg2dergausspt_1d
    use convert_2d_module, only : cellavg2cc_2d
    use variables_module, only : ctoprim
    use transport_properties, only : get_transport_properties

    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    double precision, intent(in   ) :: dxinv(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),NVAR)

    double precision, dimension(:,:,:), pointer ::  Uag
    double precision, dimension(:,:,:), pointer :: dUag
    double precision, allocatable, target ::  U1(:,:,:),  U2(:,:,:)
    double precision, allocatable, target :: dU1(:,:,:), dU2(:,:,:)
    double precision, allocatable :: Qc(:,:,:), Qf(:,:,:), dmom(:,:,:)
    double precision, allocatable :: mu(:,:), xi(:,:), lam(:,:), Ddia(:,:,:)
    integer :: i, j, n, g
    integer :: g2lo(2), g2hi(2), g3lo(2), g3hi(2)
    integer :: tlo(3), thi(3), Qclo(3), Qchi(3), Qflo(3), Qfhi(3)
    double precision, parameter :: fac = 0.5d0 ! due to Gauss quadrature

    g3lo = lo-3;  g3hi = hi+3
    allocate(U1(g3lo(1):g3hi(1),g3lo(2):g3hi(2),NVAR))
    allocate(U2(g3lo(1):g3hi(1),g3lo(2):g3hi(2),NVAR))

    g2lo = lo-2;  g2hi = hi+2
    allocate(dU1(g2lo(1):g2hi(1),g2lo(2):g2hi(2),3))
    allocate(dU2(g2lo(1):g2hi(1),g2lo(2):g2hi(2),3))    

    tlo = 1; thi = 1; Qclo = 1;  Qchi = 1;  Qflo = 1;  Qfhi = 1

    Qclo(1:2) = lo(1:2)-2
    Qchi(1:2) = hi(1:2)+2

    Qflo(1:2) = lo(1:2)
    Qfhi(1:2) = hi(1:2)+1

    allocate(Qc  (Qclo(1):Qchi(1),Qclo(2):Qchi(2),QCVAR))
    allocate(Qf  (Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR))
    allocate(mu  (Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(xi  (Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(lam (Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(Ddia(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC))
    allocate(dmom(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),3))

    ! ----- compute x-direction flux first -----

    ! cell-average of ? => cell-avg-in-x and Gauss-point-in-y of d?/dy
    do n=1,3
       do i=lo(1)-2,hi(1)+2
          call cellavg2dergausspt_1d(lo(2),hi(2), U(i,:,n), Ulo(2), Uhi(2), &
               dU1(i,:,n), dU2(i,:,n), g2lo(2), g2hi(2))
       end do
    end do

    ! cell-average => cell-avg-in-x and Gauss-point-in-y
    do n=1,NVAR
       do i=lo(1)-3,hi(1)+3
          call cellavg2gausspt_1d(lo(2),hi(2), U(i,:,n), Ulo(2), Uhi(2), &
               U1(i,:,n), U2(i,:,n), g3lo(2), g3hi(2))
       end do
    end do
    
    do g=1,2
       
       if (g .eq. 1) then
          Uag  =>  U1
          dUag => dU1
       else
          Uag  =>  U2
          dUag => dU2
       end if

       do n=1,3
          do j=lo(2),hi(2)
             ! cell-avg-in-x and Gauss-point-in-y => xface and Gauss-point-in-y
             call cellavg2face_1d(lo(1),hi(1)+1, dUag(:,j,n),g2lo(1),g2hi(1), &
                  dmom(:,j,n),Qflo(1),Qfhi(1))
          end do
       end do

       do n=1,NVAR
          do j=lo(2),hi(2)
             ! cell-avg-in-x and Gauss-point-in-y => xface and Gauss-point-in-y
             call cellavg2face_1d(lo(1),hi(1)+1, Uag(:,j,n),g3lo(1),g3hi(1), &
                  Qf(:,j,n),Qflo(1),Qfhi(1))
          end do
       end do

       tlo(1) = lo(1)-2
       tlo(2) = lo(2)
       thi(1) = hi(1)+2
       thi(2) = hi(2)
       ! cell-avg-in-x and Gauss-point-in-y => cell-center-in-x and Gauss-point-in-y
       do n=1,NVAR
          call cellavg2cc_2d(tlo(1:2),thi(1:2), Uag(:,:,n),g3lo,g3hi, &
               Qc(:,:,n),Qclo(1:2),Qchi(1:2),idir=1)
       end do

       tlo(1:2) = lo(1:2)
       thi(1) = hi(1)+1
       thi(2) = hi(2)
       call ctoprim(tlo,thi, Qf, Qflo,Qfhi,QFVAR)

       ! transport coefficients on face
       call get_transport_properties(tlo,thi, Qf,Qflo,Qfhi,QFVAR, mu,xi,lam,Ddia,Qflo,Qfhi)

       tlo(1) = lo(1)-2
       thi(1) = hi(1)+2
       call ctoprim(tlo,thi, Qc, Qclo,Qchi,QCVAR)

       tlo(1:2) = lo
       thi(1) = hi(1)+1
       thi(2) = hi(2)
       call comp_diff_flux_x(tlo(1:2), thi(1:2), fx, fxlo, fxhi, &
            Qf, mu, xi, lam, Ddia, dmom, Qflo, Qfhi, &
            Qc, Qclo, Qchi, dxinv, fac)

       Nullify(Uag,dUag)
    end do

    ! ----- compute y-direction flux -----

    ! cell-average of ? => cell-avg-in-y and Gauss-point-in-x of d?/dx
    do n=1,3
       do j=lo(2)-2,hi(2)+2
          call cellavg2dergausspt_1d(lo(1),hi(1), U(:,j,n), Ulo(1), Uhi(1), &
               dU1(:,j,n), dU2(:,j,n), g2lo(1), g2hi(1))
       end do
    end do

    ! cell-average => cell-avg-in-y and Gauss-point-in-x
    do n=1,NVAR
       do j=lo(2)-3,hi(2)+3
          call cellavg2gausspt_1d(lo(1),hi(1), U(:,j,n), Ulo(1), Uhi(1), &
               U1(:,j,n), U2(:,j,n), g3lo(1), g3hi(1))
       end do
    end do

    do g=1,2
       
       if (g .eq. 1) then
          Uag  =>  U1
          dUag => dU1
       else
          Uag  =>  U2
          dUag => dU2
       end if

       do n=1,3
          do i=lo(1),hi(1)
             ! cell-avg-in-y and Gauss-point-in-x => yface and Gauss-point-in-x
             call cellavg2face_1d(lo(2),hi(2)+1, dUag(i,:,n),g2lo(2),g2hi(2), &
                  dmom(i,:,n),Qflo(2),Qfhi(2))
          end do
       end do

       do n=1,NVAR
          do i=lo(1),hi(1)
             ! cell-avg-in-y and Gauss-point-in-x => yface and Gauss-point-in-x
             call cellavg2face_1d(lo(2),hi(2)+1, Uag(i,:,n),g3lo(2),g3hi(2), &
                  Qf(i,:,n),Qflo(2),Qfhi(2))
          end do
       end do

       tlo(1) = lo(1)
       tlo(2) = lo(2)-2
       thi(1) = hi(1)
       thi(2) = hi(2)+2
       do n=1,NVAR
          ! cell-avg-in-y and Gauss-point-in-x => cell-center-in-y and Gauss-point-in-x
          call cellavg2cc_2d(tlo(1:2),thi(1:2), Uag(:,:,n),g3lo,g3hi, &
               Qc(:,:,n),Qclo(1:2),Qchi(1:2),idir=2)
       end do

       tlo(1:2) = lo(1:2)
       thi(1) = hi(1)
       thi(2) = hi(2)+1
       call ctoprim(tlo,thi, Qf, Qflo,Qfhi,QFVAR)

       ! transport coefficients on face
       call get_transport_properties(tlo,thi, Qf,Qflo,Qfhi,QFVAR, mu,xi,lam,Ddia,Qflo,Qfhi)

       tlo(2) = lo(2)-2
       thi(2) = hi(2)+2
       call ctoprim(tlo,thi, Qc, Qclo,Qchi,QCVAR)

       tlo(1:2) = lo
       thi(1) = hi(1)
       thi(2) = hi(2)+1
       call comp_diff_flux_y(tlo(1:2), thi(1:2), fy, fylo, fyhi, &
            Qf, mu, xi, lam, Ddia, dmom, Qflo, Qfhi, &
            Qc, Qclo, Qchi, dxinv, fac)

       Nullify(Uag,dUag)
    end do

    deallocate(U1,U2,dU1,dU2,Qc,Qf,mu,xi,lam,Ddia,dmom)

  end subroutine difterm


  subroutine comp_diff_flux_x(lo, hi, flx, flo, fhi, &
       Qf, mu, xi, lam, Ddia, dmom, Qflo, Qfhi, &
       Qc, Qclo, Qchi, dxinv, fac)

    use meth_params_module
    use derivative_stencil_module, only : FD4

    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), Qflo(2), Qfhi(2), Qclo(2), Qchi(2)
    double precision, intent(in) :: dxinv(2), fac
    double precision, intent(inout) ::  flx( flo(1): fhi(1), flo(2): fhi(2),NVAR)
    double precision, intent(in   ) ::   Qf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR)
    double precision, intent(in   ) ::   mu(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::   xi(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::  lam(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) :: Ddia(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC)
    double precision, intent(in   ) :: dmom(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),3)
    double precision, intent(in   ) ::   Qc(Qclo(1):Qchi(1),Qclo(2):Qchi(2),QCVAR)

    integer :: i, j, n, UYN, QYN, QXN, QHN
    double precision :: tauxx, tauxy, dudx, dudy, dvdx, dvdy, divu, rhoinv
    double precision :: dTdx, dXdx, Vd
    double precision, dimension(lo(1):hi(1)) :: dlnpdx, Vc
    double precision, parameter :: twoThirds = 2.d0/3.d0

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          ! viscous stress
          dudx = dxinv(1)*(FD4(-2)*Qc(i-2,j,QU) + FD4(-1)*Qc(i-1,j,QU) &
               + FD4(0)*Qc(i,j,QU) + FD4(1)*Qc(i+1,j,QU))
          dvdx = dxinv(1)*(FD4(-2)*Qc(i-2,j,QV) + FD4(-1)*Qc(i-1,j,QV) &
               + FD4(0)*Qc(i,j,QV) + FD4(1)*Qc(i+1,j,QV))
          rhoinv = 1.d0/Qf(i,j,QRHO)
          dudy = dxinv(2)*rhoinv*(dmom(i,j,2)-Qf(i,j,QU)*dmom(i,j,1))
          dvdy = dxinv(2)*rhoinv*(dmom(i,j,3)-Qf(i,j,QV)*dmom(i,j,1))
          divu = dudx + dvdy
          tauxx = mu(i,j)*(2.d0*dudx-twoThirds*divu) + xi(i,j)*divu
          tauxy = mu(i,j)*(dudy+dvdx)
          flx(i,j,UMX)   = flx(i,j,UMX)   - fac*tauxx
          flx(i,j,UMY)   = flx(i,j,UMY)   - fac*tauxy
          flx(i,j,UEDEN) = flx(i,j,UEDEN) - fac*(tauxx*Qf(i,j,QU)+tauxy*Qf(i,j,QV))

          ! thermal conduction
          dTdx = dxinv(1) * (FD4(-2)*Qc(i-2,j,QTEMP) + FD4(-1)*Qc(i-1,j,QTEMP) &
               + FD4(0)*Qc(i,j,QTEMP) + FD4(1)*Qc(i+1,j,QTEMP))
          flx(i,j,UEDEN) = flx(i,j,UEDEN) - fac*lam(i,j)*dTdx

          ! compute dpdx
          dlnpdx(i) = dxinv(1) * (FD4(-2)*Qc(i-2,j,QPRES) + FD4(-1)*Qc(i-1,j,QPRES) &
               + FD4(0)*Qc(i,j,QPRES) + FD4(1)*Qc(i+1,j,QPRES)) / Qf(i,j,QPRES)
          Vc(i) = 0.d0
       end do

       do n=1,NSPEC
          UYN = UFS+n-1
          QYN = QFY+n-1
          QXN = QFX+n-1
          QHN = QFH+n-1
          do i = lo(1), hi(1)
             dXdx = dxinv(1) * (FD4(-2)*Qc(i-2,j,QXN) + FD4(-1)*Qc(i-1,j,QXN) &
                  + FD4(0)*Qc(i,j,QXN) + FD4(1)*Qc(i+1,j,QXN))
             Vd = -Ddia(i,j,n)*(dXdx + (Qf(i,j,QXN)-Qf(i,j,QYN))*dlnpdx(i))
             
             flx(i,j,UYN) = flx(i,j,UYN) + fac*Vd
             Vc(i) = Vc(i) + Vd
             flx(i,j,UEDEN) = flx(i,j,UEDEN) + fac*Vd*Qf(i,j,QHN)
          end do
       end do

       do n=1,NSPEC
          UYN = UFS+n-1
          QYN = QFY+n-1
          QHN = QFH+n-1
          do i = lo(1), hi(1)
             flx(i,j,UYN )  = flx(i,j,UYN  ) - (fac*Qf(i,j,QYN)*Vc(i))
             flx(i,j,UEDEN) = flx(i,j,UEDEN) - (fac*Qf(i,j,QYN)*Vc(i))*Qf(i,j,QHN)
          end do
       end do
    end do

  end subroutine comp_diff_flux_x


  subroutine comp_diff_flux_y(lo, hi, flx, flo, fhi, &
       Qf, mu, xi, lam, Ddia, dmom, Qflo, Qfhi, &
       Qc, Qclo, Qchi, dxinv, fac)

    use meth_params_module
    use derivative_stencil_module, only : FD4

    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), Qflo(2), Qfhi(2), Qclo(2), Qchi(2)
    double precision, intent(in) :: dxinv(2), fac
    double precision, intent(inout) ::  flx( flo(1): fhi(1), flo(2): fhi(2),NVAR)
    double precision, intent(in   ) ::   Qf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR)
    double precision, intent(in   ) ::   mu(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::   xi(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::  lam(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) :: Ddia(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC)
    double precision, intent(in   ) :: dmom(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),3)
    double precision, intent(in   ) ::   Qc(Qclo(1):Qchi(1),Qclo(2):Qchi(2),QCVAR)

    integer :: i, j, n, UYN, QYN, QXN, QHN
    double precision :: tauyy, tauxy, dudx, dudy, dvdx, dvdy, divu, rhoinv
    double precision :: dTdy, dXdy, Vd
    double precision, allocatable :: dlnpdy(:,:), Vc(:,:)
    double precision, parameter :: twoThirds = 2.d0/3.d0

    allocate(dlnpdy(lo(1):hi(1),lo(2):hi(2)))
    allocate(    Vc(lo(1):hi(1),lo(2):hi(2)))

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          ! viscous stress
          dudy = dxinv(2)*(FD4(-2)*Qc(i,j-2,QU) + FD4(-1)*Qc(i,j-1,QU) &
               + FD4(0)*Qc(i,j,QU) + FD4(1)*Qc(i,j+1,QU))
          dvdy = dxinv(2)*(FD4(-2)*Qc(i,j-2,QV) + FD4(-1)*Qc(i,j-1,QV) &
               + FD4(0)*Qc(i,j,QV) + FD4(1)*Qc(i,j+1,QV))
          rhoinv = 1.d0/Qf(i,j,QRHO)
          dudx = dxinv(1)*rhoinv*(dmom(i,j,2)-Qf(i,j,QU)*dmom(i,j,1))
          dvdx = dxinv(1)*rhoinv*(dmom(i,j,3)-Qf(i,j,QV)*dmom(i,j,1))
          divu = dudx + dvdy
          tauyy = mu(i,j)*(2.d0*dvdy-twoThirds*divu) + xi(i,j)*divu
          tauxy = mu(i,j)*(dudy+dvdx)
          flx(i,j,UMX)   = flx(i,j,UMX)   - fac*tauxy
          flx(i,j,UMY)   = flx(i,j,UMY)   - fac*tauyy
          flx(i,j,UEDEN) = flx(i,j,UEDEN) - fac*(tauxy*Qf(i,j,QU)+tauyy*Qf(i,j,QV))

          ! thermal conduction
          dTdy = dxinv(2) * (FD4(-2)*Qc(i,j-2,QTEMP) + FD4(-1)*Qc(i,j-1,QTEMP) &
               + FD4(0)*Qc(i,j,QTEMP) + FD4(1)*Qc(i,j+1,QTEMP))
          flx(i,j,UEDEN) = flx(i,j,UEDEN) - fac*lam(i,j)*dTdy

          ! compute dpdy
          dlnpdy(i,j) = dxinv(2) * (FD4(-2)*Qc(i,j-2,QPRES) + FD4(-1)*Qc(i,j-1,QPRES) &
               + FD4(0)*Qc(i,j,QPRES) + FD4(1)*Qc(i,j+1,QPRES)) / Qf(i,j,QPRES)
          Vc(i,j) = 0.d0
       end do
    end do

    do n=1,NSPEC
       UYN = UFS+n-1
       QYN = QFY+n-1
       QXN = QFX+n-1
       QHN = QFH+n-1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dXdy = dxinv(2) * (FD4(-2)*Qc(i,j-2,QXN) + FD4(-1)*Qc(i,j-1,QXN) &
                  + FD4(0)*Qc(i,j,QXN) + FD4(1)*Qc(i,j+1,QXN))
             Vd = -Ddia(i,j,n)*(dXdy + (Qf(i,j,QXN)-Qf(i,j,QYN))*dlnpdy(i,j))
             
             flx(i,j,UYN) = flx(i,j,UYN) + fac*Vd
             Vc(i,j) = Vc(i,j) + Vd
             flx(i,j,UEDEN) = flx(i,j,UEDEN) + fac*Vd*Qf(i,j,QHN)
          end do
       end do
    end do

    do n=1,NSPEC
       UYN = UFS+n-1
       QYN = QFY+n-1
       QHN = QFH+n-1
       do j= lo(2), hi(2)
          do i = lo(1), hi(1)
             flx(i,j,UYN )  = flx(i,j,UYN  ) - (fac*Qf(i,j,QYN)*Vc(i,j))
             flx(i,j,UEDEN) = flx(i,j,UEDEN) - (fac*Qf(i,j,QYN)*Vc(i,j))*Qf(i,j,QHN)
          end do
       end do
    end do

    deallocate(dlnpdy,Vc)

  end subroutine comp_diff_flux_y

end module difterm_module
