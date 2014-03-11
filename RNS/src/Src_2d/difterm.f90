module difterm_module

  implicit none

  private

  public :: difterm

contains

  subroutine difterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dxinv)

    use meth_params_module, only : NVAR, NSPEC, QCVAR, QFVAR, QU, QV
    use convert_module, only : cellavg2cc_2d
    use polyinterp_module, only : cc2xface_2d, cc2yface_2d, cc2DxYface_2d, cc2DyXface_2d
    use variables_module, only : ctoprim
    use transport_properties, only : get_transport_properties

    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    double precision, intent(in   ) :: dxinv(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),NVAR)

    double precision, allocatable :: Qcc(:,:,:),mucc(:,:),xicc(:,:),lamcc(:,:),Ddiacc(:,:,:)
    double precision, allocatable :: Qc1(:,:,:), Qc2(:,:,:), Qf1(:,:,:), Qf2(:,:,:)
    double precision, allocatable :: dvel1(:,:,:), dvel2(:,:,:)
    double precision, allocatable :: mu1(:,:), xi1(:,:), lam1(:,:), Ddia1(:,:,:)
    double precision, allocatable :: mu2(:,:), xi2(:,:), lam2(:,:), Ddia2(:,:,:)
    double precision, allocatable :: tmp1(:,:), tmp2(:,:)
    integer :: i, j, n, g
    integer :: g2lo(2), g2hi(2), Qflo(2), Qfhi(2), flo(2), fhi(2), tlo(3), thi(3)

    g2lo = lo-2
    g2hi = hi+2

    allocate(Qcc   (g2lo(1):g2hi(1),g2lo(2):g2hi(2),QFVAR))
    allocate(mucc  (g2lo(1):g2hi(1),g2lo(2):g2hi(2)))
    allocate(xicc  (g2lo(1):g2hi(1),g2lo(2):g2hi(2)))
    allocate(lamcc (g2lo(1):g2hi(1),g2lo(2):g2hi(2)))
    allocate(Ddiacc(g2lo(1):g2hi(1),g2lo(2):g2hi(2),NSPEC))

    allocate(Qc1(g2lo(1):g2hi(1),g2lo(2):g2hi(2),QCVAR))
    allocate(Qc2(g2lo(1):g2hi(1),g2lo(2):g2hi(2),QCVAR))    

    allocate(tmp1(g2lo(1):g2hi(1),g2lo(2):g2hi(2)))
    allocate(tmp2(g2lo(1):g2hi(1),g2lo(2):g2hi(2)))    

    do n=1,NVAR
       call cellavg2cc_2d(g2lo,g2hi, U(:,:,n), Ulo,Uhi, Qcc(:,:,n), g2lo,g2hi)
    end do

    tlo(1:2) = g2lo
    tlo(3) = 1
    thi(1:2) = g2hi
    thi(3) = 1
    call ctoprim(tlo, thi, Qcc, tlo, thi, QFVAR)

    ! transport coefficients at cell centers
    call get_transport_properties(tlo,thi, Qcc,tlo,thi,QFVAR, &
         mucc,xicc,lamcc,Ddiacc, tlo,thi)

    Qflo = lo
    Qfhi = hi+1

    allocate(dvel1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2))
    allocate(dvel2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2))

    allocate(Qf1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR))
    allocate(Qf2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR))

    allocate(  mu1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(  xi1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate( lam1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(Ddia1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC))

    allocate(  mu2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(  xi2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate( lam2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(Ddia2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC))

    ! ----- compute x-direction flux first -----

    ! cell center => Gauss points on x-face
    call cc2xface_2d(lo,hi,  mucc, g2lo, g2hi,  mu1,  mu2, Qflo, Qfhi,& 
         tmp1, tmp2, g2lo, g2hi)
    call cc2xface_2d(lo,hi,  xicc, g2lo, g2hi,  xi1,  xi2, Qflo, Qfhi,& 
         tmp1, tmp2, g2lo, g2hi)
    call cc2xface_2d(lo,hi, lamcc, g2lo, g2hi, lam1, lam2, Qflo, Qfhi,& 
         tmp1, tmp2, g2lo, g2hi)
    do n=1,NSPEC
       call cc2xface_2d(lo,hi, Ddiacc(:,:,n), g2lo, g2hi, &
            Ddia1(:,:,n), Ddia2(:,:,n), Qflo, Qfhi,& 
            tmp1, tmp2, g2lo, g2hi)
    end do

    ! cell-center => Qc: center-in-x and Gauss-point-in-y 
    !                Qf: xface and Gauss-point-in-y
    do n=1,QCVAR
       call cc2xface_2d(lo,hi, Qcc(:,:,n), g2lo, g2hi, &
            Qf1(:,:,n), Qf2(:,:,n), Qflo, Qfhi, &
            Qc1(:,:,n), Qc2(:,:,n), g2lo, g2hi)
    end do
    do n=QCVAR+1,QFVAR
       call cc2xface_2d(lo,hi, Qcc(:,:,n), g2lo, g2hi, &
            Qf1(:,:,n), Qf2(:,:,n), Qflo, Qfhi, &
            tmp1, tmp2, g2lo, g2hi)
    end do

    ! cell-average of ? => xface & Gauss-point-in-y of d?/dy
    call cc2DyXface_2d(lo,hi, Qcc(:,:,QU), g2lo, g2hi, &
         dvel1(:,:,1), dvel2(:,:,1), Qflo, Qfhi, &
         tmp1, tmp2, g2lo, g2hi)
    call cc2DyXface_2d(lo,hi, Qcc(:,:,QV), g2lo, g2hi, &
         dvel1(:,:,2), dvel2(:,:,2), Qflo, Qfhi, &
         tmp1, tmp2, g2lo, g2hi)

    flo = lo
    fhi(1) = hi(1)+1
    fhi(2) = hi(2)
    call comp_diff_flux_x(flo, fhi, fx, fxlo, fxhi, &
         Qf1, mu1, xi1, lam1, Ddia1, dvel1, Qflo, Qfhi, &
         Qc1, g2lo, g2hi, dxinv, 0.5d0)
    call comp_diff_flux_x(flo, fhi, fx, fxlo, fxhi, &
         Qf2, mu2, xi2, lam2, Ddia2, dvel2, Qflo, Qfhi, &
         Qc2, g2lo, g2hi, dxinv, 0.5d0)
    

    ! ----- compute y-direction flux -----

    ! cell center => Gauss points on x-face
    call cc2yface_2d(lo,hi,  mucc, g2lo, g2hi,  mu1,  mu2, Qflo, Qfhi,& 
         tmp1, tmp2, g2lo, g2hi)
    call cc2yface_2d(lo,hi,  xicc, g2lo, g2hi,  xi1,  xi2, Qflo, Qfhi,& 
         tmp1, tmp2, g2lo, g2hi)
    call cc2yface_2d(lo,hi, lamcc, g2lo, g2hi, lam1, lam2, Qflo, Qfhi,& 
         tmp1, tmp2, g2lo, g2hi)
    do n=1,NSPEC
       call cc2yface_2d(lo,hi, Ddiacc(:,:,n), g2lo, g2hi, &
            Ddia1(:,:,n), Ddia2(:,:,n), Qflo, Qfhi,& 
            tmp1, tmp2, g2lo, g2hi)
    end do

    ! cell-center => Qc: center-in-x and Gauss-point-in-y 
    !                Qf: xface and Gauss-point-in-y
    do n=1,QCVAR
       call cc2yface_2d(lo,hi, Qcc(:,:,n), g2lo, g2hi, &
            Qf1(:,:,n), Qf2(:,:,n), Qflo, Qfhi, &
            Qc1(:,:,n), Qc2(:,:,n), g2lo, g2hi)
    end do
    do n=QCVAR+1,QFVAR
       call cc2yface_2d(lo,hi, Qcc(:,:,n), g2lo, g2hi, &
            Qf1(:,:,n), Qf2(:,:,n), Qflo, Qfhi, &
            tmp1, tmp2, g2lo, g2hi)
    end do

    ! cell-average of ? => yface and Gauss-point-in-x of d?/dx
    call cc2DxYface_2d(lo,hi, Qcc(:,:,QU), g2lo, g2hi, &
         dvel1(:,:,1), dvel2(:,:,1), Qflo, Qfhi, &
         tmp1, tmp2, g2lo, g2hi)
    call cc2DxYface_2d(lo,hi, Qcc(:,:,QV), g2lo, g2hi, &
         dvel1(:,:,2), dvel2(:,:,2), Qflo, Qfhi, &
         tmp1, tmp2, g2lo, g2hi)

    flo = lo
    fhi(1) = hi(1)
    fhi(2) = hi(2)+1
    call comp_diff_flux_y(flo, fhi, fy, fylo, fyhi, &
         Qf1, mu1, xi1, lam1, Ddia1, dvel1, Qflo, Qfhi, &
         Qc1, g2lo, g2hi, dxinv, 0.5d0)
    call comp_diff_flux_y(flo, fhi, fy, fylo, fyhi, &
         Qf2, mu2, xi2, lam2, Ddia2, dvel2, Qflo, Qfhi, &
         Qc2, g2lo, g2hi, dxinv, 0.5d0)

    deallocate(Qcc,mucc,xicc,lamcc,Ddiacc)
    deallocate(Qc1,Qc2,tmp1,tmp2)
    deallocate(Qf1,Qf2,dvel1,dvel2)
    deallocate(mu1 ,xi1 ,lam1 ,Ddia1)
    deallocate(mu2 ,xi2 ,lam2 ,Ddia2)

  end subroutine difterm


  subroutine comp_diff_flux_x(lo, hi, flx, flo, fhi, &
       Qf, mu, xi, lam, Ddia, dvel, Qflo, Qfhi, &
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
    double precision, intent(in   ) :: dvel(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2)
    double precision, intent(in   ) ::   Qc(Qclo(1):Qchi(1),Qclo(2):Qchi(2),QCVAR)

    integer :: i, j, n, UYN, QYN, QXN, QHN
    double precision :: tauxx, tauxy, dudx, dudy, dvdx, dvdy, divu
    double precision :: dTdx, dXdx, Vd
    double precision :: ek, rhovn
    double precision, dimension(lo(1):hi(1)) :: dlnpdx, Vc
    double precision, parameter :: twoThirds = 2.d0/3.d0

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          ! viscous stress
          dudx = dxinv(1)*(FD4(-2)*Qc(i-2,j,QU) + FD4(-1)*Qc(i-1,j,QU) &
               + FD4(0)*Qc(i,j,QU) + FD4(1)*Qc(i+1,j,QU))
          dvdx = dxinv(1)*(FD4(-2)*Qc(i-2,j,QV) + FD4(-1)*Qc(i-1,j,QV) &
               + FD4(0)*Qc(i,j,QV) + FD4(1)*Qc(i+1,j,QV))
          dudy = dxinv(2)*dvel(i,j,1)
          dvdy = dxinv(2)*dvel(i,j,2)
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

    if (.not. do_weno) then
       ! compute hyperbolic flux
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhovn = fac*Qf(i,j,QRHO)*Qf(i,j,QU)
             flx(i,j,URHO) = flx(i,j,URHO) + rhovn
             flx(i,j,UMX ) = flx(i,j,UMX ) + rhovn*Qf(i,j,QU) + fac*Qf(i,j,QPRES)
             flx(i,j,UMY ) = flx(i,j,UMY ) + rhovn*Qf(i,j,QV)

             ek = 0.5d0*(Qf(i,j,QU)**2+Qf(i,j,QV)**2)
             flx(i,j,UEDEN) = flx(i,j,UEDEN) + rhovn*ek
             do n=1,NSPEC
                flx(i,j,UEDEN) = flx(i,j,UEDEN) + rhovn*Qf(i,j,QFH+n-1)
                flx(i,j,UFS+n-1) = flx(i,j,UFS+n-1) + rhovn*Qf(i,j,QFY+n-1)
             end do
          end do
       end do
    end if

  end subroutine comp_diff_flux_x


  subroutine comp_diff_flux_y(lo, hi, flx, flo, fhi, &
       Qf, mu, xi, lam, Ddia, dvel, Qflo, Qfhi, &
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
    double precision, intent(in   ) :: dvel(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2)
    double precision, intent(in   ) ::   Qc(Qclo(1):Qchi(1),Qclo(2):Qchi(2),QCVAR)

    integer :: i, j, n, UYN, QYN, QXN, QHN
    double precision :: tauyy, tauxy, dudx, dudy, dvdx, dvdy, divu, rhoinv
    double precision :: dTdy, dXdy, Vd
    double precision :: ek, rhovn
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
          dudx = dxinv(1)*dvel(i,j,1)
          dvdx = dxinv(1)*dvel(i,j,2)
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
    
    if (.not. do_weno) then
       ! compute hyperbolic flux
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhovn = fac*Qf(i,j,QRHO)*Qf(i,j,QV)
             flx(i,j,URHO) = flx(i,j,URHO) + rhovn
             flx(i,j,UMX ) = flx(i,j,UMX ) + rhovn*Qf(i,j,QU)
             flx(i,j,UMY ) = flx(i,j,UMY ) + rhovn*Qf(i,j,QV) + fac*Qf(i,j,QPRES)

             ek = 0.5d0*(Qf(i,j,QU)**2+Qf(i,j,QV)**2)
             flx(i,j,UEDEN) = flx(i,j,UEDEN) + rhovn*ek
             do n=1,NSPEC
                flx(i,j,UEDEN) = flx(i,j,UEDEN) + rhovn*Qf(i,j,QFH+n-1)
                flx(i,j,UFS+n-1) = flx(i,j,UFS+n-1) + rhovn*Qf(i,j,QFY+n-1)
             end do
          end do
       end do
    end if

  end subroutine comp_diff_flux_y

end module difterm_module
