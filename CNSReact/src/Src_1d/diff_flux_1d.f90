module diff_flux_module

  use cdwrk_module
  use meth_params_module
 
  implicit none
 
  public :: diffFlux
 
contains

      subroutine diffFlux(lo,hi, &
                          q,q_l1,q_h1,
                          flux1,flux1_l1,flux1_h1,&
                          dx,dt,initFlux)

      ! Get diffusion flux on faces surrounding Box(lo,hi), requires 
      ! valid data on region Box(lo,hi).grow(1)

      implicit none

      integer :: lo(1),hi(1)
      integer :: q_l1, q_h1
      integer :: flux1_l1, flux1_h1
      double precision ::     q(    q_l1:    q_h1, QVAR)
      double precision :: flux1(flux1_l1:flux1_h1,NVAR)
      double precision :: dx,dt
      integer          :: initFlux

      ! Local variables
      integer          :: loD(1),hiD(1)
      double precision, allocatable ::    D(:,:)
      double precision, allocatable :: TEMP(:)
      double precision, allocatable ::   CP(:)

      loD(1) = lo(1)-1
      hiD(1) = hi(1)+1

      allocate(   D(loD(1):hiD(1),Nspec+3))
      allocate(TEMP(q_l1:q_h1))
      allocate(  CP(q_l1:q_h1))

      if (initFlux.eq.1) then
         flux1(:,:) = 0.
      end if

      call transCoef(loD,hiD,q,TEMP,CP,q_l1,q_h1, &
                     D,loD(1),hiD(1))

      call myViscFlux(lo,hi, &
                      q,  q_l1,  q_h1,  &
                      D,loD(1),hiD(1), &
                      flux1,flux1_l1,flux1_h1,
                      dx,dt)

      deallocate(D,TEMP,CP)

      end subroutine diffFlux
! ::
! :: ----------------------------------------------------------
! ::

      subroutine transCoef(lo,hi,q,TEMP,CP,q_l1,q_h1,D,D_l1,D_h1)

      ! Note that lo,hi in this routine correspond to lo-1,hi+1 from the calling routine


      implicit none

      integer          :: lo(1), hi(1)
      integer          :: q_l1, q_h1
      integer          :: D_l1, D_h1
      double precision ::    q(q_l1:q_h1, QVAR)
      double precision :: TEMP(q_l1:q_h1)
      double precision ::   CP(q_l1:q_h1) 
      double precision ::    D(D_l1:D_h1, Nspec+3)

      ! Local variables
      integer          :: i,doVelVisc,doTemp
!     double precision :: Dmax

      do i = lo(1),hi(1)
         TEMP(i) = MIN(MAX(q(i,QTHERM+1),300.d0),6500.d0)
      end do
      
      doTemp = 1
      doVelVisc = 1

      call mixavg_rhodiff_temp_pres(lo, hi, &
           D,  &
             D_l1,D_h1,&
           TEMP, &
             q_l1,q_h1,&
           q(q_l1,QTHERM+NADV+1), &
             q_l1,q_h1, &
           q(q_l1,QPRES), &
             q_l1,q_h1, &
           doTemp, doVelVisc)

      call FORT_CPMIXfromTY(lo, hi, &
           CP, &
             q_l1,q_h1,&
           TEMP, &
             q_l1,q_h1, &
           q(q_l1,QTHERM+NADV+1), &
             q_l1,q_h1)

      ! Replace lambda with lambda/Cp
      D(:,Nspec+1) = D(:,Nspec+1)/CP(:)

!     Dmax = -1.d0
!     do i = lo(1),hi(1)
!        do n=1,Nspec
!           Dmax = MAX(Dmax,D(i,n) / q(i,QRHO))
!        enddo
!     end do
!     print*, '********************      Dmax: ', Dmax

      end subroutine transCoef
! ::
! :: ----------------------------------------------------------
! ::
      subroutine myViscFlux(lo,hi,&
                            q,q_l1,q_h1, &
                            D,D_l1,D_h1, &
                            flux1,flux1_l1,flux1_h1, &
                            dx,dt)


      implicit none

      integer          :: lo(1),hi(1)
      integer          :: q_l1, q_h1
      integer          :: D_l1, D_h1
      integer          :: flux1_l1, flux1_h1
      double precision :: q(q_l1:q_h1,,QVAR)
      double precision :: D(D_l1:D_h1,,Nspec+3)
      double precision :: flux1(flux1_l1:flux1_h1,NVAR)

      double precision :: dx,dt

      ! Local variables
      integer          :: i,n,ifirstSp
      double precision :: tauxxm
      double precision :: divxm
      double precision :: muxm
      double precision :: lamxm
      double precision :: kxm
      double precision :: phiflx

      double precision :: rhoDm(maxspec)
      double precision, parameter ::    two3rd = 2.d0/3.d0

      ifirstSp  = NTHERM+NADV+1

      ! compute x fluxes
      do i = lo(1),hi(1)+1
               
         muxm =  0.5d0*(D(i-1,Nspec+2) + D(i,Nspec+2))
         lamxm = -two3rd*muxm
         kxm =   0.5d0*(D(i-1,Nspec+1) + D(i,Nspec+1))
         do n=1,Nspec
            rhoDm(n) = 0.5d0*(D(i-1,n) + D(i,n))
         end do
         
         tauxxm = 2.d0*muxm*(q(i,QU) - q(i-1,QU))/dx 
         
         divxm = lamxm*(q(i,QU) - q(i-1,QU))/dx  

         phiflx =  tauxxm*(q(i,QU)+q(i-1,QU)) &
              +     divxm*(q(i,QU)+q(i-1,QU))
         
         flux1(i,2) = flux1(i,2) - dt*(tauxxm+divxm)
         flux1(i,3) = flux1(i,3) - dt*(0.5d0*phiflx  &
              + kxm*( (q(i  ,QREINT)+q(i  ,QPRES))/q(i  ,QRHO) &
              -       (q(i-1,QREINT)+q(i-1,QPRES))/q(i-1,QRHO) )/dx)

         do n=1,Nspec
            flux1(i,n+ifirstSp-1) = flux1(i,n+ifirstSp-1) &
                 - dt*dy*dz*rhoDm(n)* &
                 (q(i,QTHERM+NADV+n)-q(i-1,QTHERM+NADV+n))/dx
         end do
         
      enddo


      end subroutine myViscFlux
! ::
! :: ----------------------------------------------------------
! ::
      subroutine mixavg_rhodiff_temp_pres(lo, hi, &
                          rd,rd_l1,rd_h1, &
                           T, T_l1, T_h1, &
                           Y, Y_l1, Y_h1, &
                           P, P_l1, P_h1, &
                           do_temp, do_VelVisc)


      implicit none

      integer          :: lo(1),hi(1)
      integer          :: do_temp, do_VelVisc
      integer          :: rd_l1,rd_h1
      integer          ::  T_l1, T_h1
      integer          ::  Y_l1, Y_h1
      integer          ::  P_l1, P_h1
      double precision :: rd(rd_l1:rd_h1,:)
      double precision ::  T( T_l1: T_h1)
      double precision ::  Y( Y_l1: Y_h1,:)
      double precision ::  P( P_l1: P_h1)

      ! Local variables
      integer          :: i,  n
      double precision :: Yt(maxspec), Dt(maxspec)
      double precision :: scal, tscal, Wavg, RHO, Tt, Pt, invmwt(maxspec)
      double precision :: alpha, l1, l2, X(maxspec), CPMS(maxspec)

      call CKWT(IWRK(ckbi),RWRK(ckbr),invmwt)

      do n=1,Nspec
         invmwt(n) = 1.d0 / invmwt(n)
      end do

      scal  = 0.1d0
      tscal = 1.d-5
      do i=lo(1),hi(1)
         do n=1,Nspec
            Yt(n) = Y(i,n)
         end do
         Tt = MAX(T(i),TMIN_TRANS)
         Pt = P(i)*1.d1
         CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
         CALL CKCPMS(Tt,IWRK(ckbi),RWRK(ckbr),CPMS)
         CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
         CALL EGSPAR(Tt,X,Yt,CPMS,RWRK(egbr),IWRK(egbi))
         CALL EGSV1(Pt,Tt,Yt,Wavg,RWRK(egbr),Dt)
         CALL CKRHOY(Pt,Tt,Yt,IWRK(ckbi),RWRK(ckbr),RHO)
         do n=1,Nspec
            rd(i,n) = RHO * Wavg * invmwt(n) * Dt(n) * scal
         end do

         if (do_temp .ne. 0) then
            alpha = 1
            CALL EGSL1(alpha,Tt,X,RWRK(egbr),l1)
            alpha = -1
            CALL EGSL1(alpha,Tt,X,RWRK(egbr),l2)
            rd(i,Nspec+1) = 0.5d0 * (l1 + l2) * tscal
         endif

         if (do_VelVisc .ne. 0) then
            CALL EGSE3(Tt,Yt,RWRK(egbr),rd(i,Nspec+2))
            rd(i,Nspec+2) = rd(i,Nspec+2) * scal
         endif

      end do

      end subroutine mixavg_rhodiff_temp_pres

end module diff_flux_module
