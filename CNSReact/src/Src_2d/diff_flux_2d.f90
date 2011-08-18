module diff_flux_module

  use cdwrk_module
  use meth_params_module
 
  implicit none
 
  public :: diffFlux
 
contains

      subroutine diffFlux(lo,hi, &
                          q,q_l1,q_l2,q_h1,q_h2,&
                          flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                          flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                          dx,dy,dt,initFlux)

      ! Get diffusion flux on faces surrounding Box(lo,hi), requires 
      ! valid data on region Box(lo,hi).grow(1)

      implicit none

      integer lo(2),hi(2)
      integer :: q_l1, q_l2, q_h1, q_h2
      integer :: flux1_l1, flux1_l2, flux1_h1, flux1_h2
      integer :: flux2_l1, flux2_l2, flux2_h1, flux2_h2
      double precision :: q(q_l1:q_h1, q_l2:q_h2, QVAR)
      double precision :: flux1(flux1_l1:flux1_h1, flux1_l2:flux1_h2, NVAR)
      double precision :: flux2(flux2_l1:flux2_h1, flux2_l2:flux2_h2, NVAR)
      double precision :: dx,dy,dt
      integer          :: initFlux

      ! Local variables
      integer          :: loD(2),hiD(2)
      double precision, allocatable ::    D(:,:,:,:)
      double precision, allocatable :: TEMP(:,:,:)
      double precision, allocatable ::   CP(:,:,:)

      loD(:) = lo(:)-1
      hiD(:) = hi(:)+1

      allocate(   D(loD(1):hiD(1),loD(2):hiD(2),Nspec+3))
      allocate(TEMP(q_l1:q_h1,q_l2:q_h2))
      allocate(  CP(q_l1:q_h1,q_l2:q_h2))

      if (initFlux.eq.1) then
         flux1(:,:,:,:) = 0.
         flux2(:,:,:,:) = 0.
      end if

      call transCoef(loD,hiD,q,TEMP,CP,q_l1,q_l2,q_h1,q_h2, &
                     D,loD(1),loD(2),hiD(1),hiD(2))

      call myViscFlux(lo,hi, &
                      q,  q_l1,  q_l2, q_h1,  q_h2,  &
                      D,loD(1),loD(2),hiD(1),hiD(2), &
                      flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2,&
                      flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2,&
                      dx,dy,dt)

      deallocate(D,TEMP,CP)

      end subroutine diffFlux
! ::
! :: ----------------------------------------------------------
! ::

      subroutine transCoef(lo,hi,q,TEMP,CP, &
                           q_l1,q_l2,q_h1,q_h2, &
                           D, &
                           D_l1,D_l2,D_h1,D_h2)

      ! Note that lo,hi in this routine correspond to lo-1,hi+1 from the calling routine


      implicit none

      integer          :: lo(2), hi(2)
      integer          :: q_l1, q_l2, q_h1, q_h2
      integer          :: D_l1, D_l2, D_h1, D_h2
      double precision ::    q(q_l1:q_h1, q_l2:q_h2, QVAR)
      double precision :: TEMP(q_l1:q_h1, q_l2:q_h2)
      double precision ::   CP(q_l1:q_h1, q_l2:q_h2) 
      double precision ::    D(D_l1:D_h1, D_l2:D_h2, Nspec+3)

      ! Local variables
      integer          :: i,j,doVelVisc,doTemp
!     double precision :: Dmax

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            TEMP(i,j) = MIN(MAX(q(i,j,QTHERM+1),300.d0),6500.d0)
         end do
      end do
      
      doTemp = 1
      doVelVisc = 1

      call mixavg_rhodiff_temp_pres(lo, hi, &
           D,  &
             D_l1,D_l2,D_h1,D_h2,&
           TEMP, &
             q_l1,q_l2,q_h1,q_h2,&
           q(q_l1,q_l2,QTHERM+NADV+1), &
             q_l1,q_l2,q_h1,q_h2, &
           q(q_l1,q_l2,QPRES), &
             q_l1,q_l2,q_h1,q_h2, &
           doTemp, doVelVisc)

      call FORT_CPMIXfromTY(lo, hi, &
           CP, &
             q_l1,q_l2,q_h1,q_h2,&
           TEMP, &
             q_l1,q_l2,q_h1,q_h2, &
           q(q_l1,q_l2,QTHERM+NADV+1), &
             q_l1,q_l2,q_h1,q_h2)

      ! Replace lambda with lambda/Cp
      D(:,:,:,Nspec+1) = D(:,:,:,Nspec+1)/CP(:,:,:)

!     Dmax = -1.d0
!     do j = lo(2),hi(2)
!        do i = lo(1),hi(1)
!           do n=1,Nspec
!              Dmax = MAX(Dmax,D(i,j,n) / q(i,j,QRHO))
!           enddo
!        end do
!     end do
!     print*, '********************      Dmax: ', Dmax

      end subroutine transCoef
! ::
! :: ----------------------------------------------------------
! ::
      subroutine myViscFlux(lo,hi,&
                            q,q_l1,q_l2,q_h1,q_h2, &
                            D,D_l1,D_l2,D_h1,D_h2, &
                            flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                            flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                            dx,dy,dt)


      implicit none

      integer          :: lo(2),hi(2)
      integer          :: q_l1, q_l2, q_h1, q_h2
      integer          :: D_l1, D_l2, D_h1, D_h2
      integer          :: flux1_l1, flux1_l2, flux1_h1, flux1_h2
      integer          :: flux2_l1, flux2_l2, flux2_h1, flux2_h2
      double precision :: q(q_l1:q_h1, q_l2:q_h2, QVAR)
      double precision :: D(D_l1:D_h1, D_l2:D_h2, Nspec+3)
      double precision :: flux1(flux1_l1:flux1_h1, flux1_l2:flux1_h2, NVAR)
      double precision :: flux2(flux2_l1:flux2_h1, flux2_l2:flux2_h2, NVAR)

      double precision :: dx,dy,dt

      ! Local variables
      integer          :: i,j,n,ifirstSp
      double precision :: tauxym,tauyym
      double precision :: tauyxm
      double precision :: divxm,divym
      double precision :: muxm,muym
      double precision :: lamxm,lamym
      double precision :: kxm,kym
      double precision :: phiflx

      double precision :: rhoDm(maxspec)
      double precision, parameter ::    two3rd = 2.d0/3.d0

      ifirstSp  = NTHERM+NADV+1

      ! compute x fluxes
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)+1
               
               muxm =  0.5d0*(D(i-1,j,Nspec+2) + D(i,j,Nspec+2))
               lamxm = -two3rd*muxm
               kxm =   0.5d0*(D(i-1,j,Nspec+1) + D(i,j,Nspec+1))
               do n=1,Nspec
                  rhoDm(n) = 0.5d0*(D(i-1,j,n) + D(i,j,n))
               end do
               
               tauxxm = 2.d0*muxm*(q(i,j,QU) - q(i-1,j,QU))/dx 
               
               divxm = lamxm*(q(i,j,QU) - q(i-1,j,QU))/dx  &
                    +.25d0*lamxm*(q(i-1,j+1,QV)+q(i,j+1,QV) &
                    -             q(i-1,j-1,QV)-q(i,j-1,QV))/dy
               
               tauyxm =  muxm*(q(i,j,QV) - q(i-1,j,QV))/dx  &
                    +.25d0*muxm*(q(i,j+1,QU)+q(i-1,j+1,QU) &
                    -            q(i,j-1,QU)-q(i-1,j-1,QU))/dy

               phiflx =  tauxxm*(q(i,j,QU)+q(i-1,j,QU)) &
                    +     divxm*(q(i,j,QU)+q(i-1,j,QU)) &
                    +    tauyxm*(q(i,j,QV)+q(i-1,j,QV)) 
               
               flux1(i,j,2) = flux1(i,j,2) - dt*dy*(tauxxm+divxm)
               flux1(i,j,3) = flux1(i,j,3) - dt*dy*tauyxm
               flux1(i,j,4) = flux1(i,j,4) - dt*dy*(0.5d0*phiflx  &
                    + kxm*( (q(i  ,j,QREINT)+q(i  ,j,QPRES))/q(i  ,j,QRHO) &
                    -       (q(i-1,j,QREINT)+q(i-1,j,QPRES))/q(i-1,j,QRHO) )/dx)

               do n=1,Nspec
                  flux1(i,j,n+ifirstSp-1) = flux1(i,j,n+ifirstSp-1) &
                       - dt*dy*rhoDm(n)* &
                       (q(i,j,QTHERM+NADV+n)-q(i-1,j,QTHERM+NADV+n))/dx
               end do
               
         enddo
      enddo

      ! compute y fluxes
      do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)
               
               muym =  0.5d0*(D(i,j-1,Nspec+2) + D(i,j,Nspec+2))
               lamym = -two3rd*muym
               kym =   0.5d0*(D(i,j-1,Nspec+1) + D(i,j,Nspec+1))
               do n=1,Nspec
                  rhoDm(n) = 0.5d0*(D(i,j-1,n) + D(i,j,n))
               end do
               
               tauxym =  muym*(q(i,j,QU) - q(i,j-1,QU))/dy  &
                    +.25d0*muym*(q(i+1,j,QV)+q(i+1,j-1,QV) &
                    -            q(i-1,j,QV)-q(i-1,j-1,QV))/dx
               
               tauyym = 2.d0*muym*(q(i,j,QV) - q(i,j-1,QV))/dy 
               
               divym = lamym*(q(i,j,QV) - q(i,j-1,QV))/dy  &
                    +.25d0*lamym*(q(i+1,j-1,QU)+q(i+1,j,QU) &
                    -             q(i-1,j-1,QU)-q(i-1,j,QU))/dx
               

               phiflx = tauxym*(q(i,j,QU)+q(i,j-1,QU)) &
                    +   tauyym*(q(i,j,QV)+q(i,j-1,QV)) &
                    +    divym*(q(i,j,QV)+q(i,j-1,QV)) 
               
               flux2(i,j,2) = flux2(i,j,2) - dt*dx*tauxym
               flux2(i,j,3) = flux2(i,j,3) - dt*dx*(tauyym+divym)
               flux2(i,j,4) = flux2(i,j,4) - dt*dy*(0.5d0*phiflx  &
                    + 2.d0*kym*( (q(i,j  ,QREINT)+q(i,j  ,QPRES))/q(i,j  ,QRHO) &
                    -            (q(i,j-1,QREINT)+q(i,j-1,QPRES))/q(i,j-1,QRHO) )/dy)

               do n=1,Nspec
                  flux2(i,j,n+ifirstSp-1) = flux2(i,j,n+ifirstSp-1) &
                       - dt*dx*rhoDm(n)* &
                       (q(i,j,QTHERM+NADV+n)-q(i,j-1,QTHERM+NADV+n))/dy
               end do
               
         enddo
      enddo

      end subroutine myViscFlux
! ::
! :: ----------------------------------------------------------
! ::
      subroutine mixavg_rhodiff_temp_pres(lo, hi, &
                          rd,rd_l1,rd_l2,rd_l3,rd_h1,rd_h2,rd_h3, &
                           T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, &
                           Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3, &
                           P, P_l1, P_l2, P_l3, P_h1, P_h2, P_h3, &
                           do_temp, do_VelVisc)


      implicit none

      integer          :: lo(3),hi(3)
      integer          :: do_temp, do_VelVisc
      integer          :: rd_l1,rd_l2,rd_l3,rd_h1,rd_h2,rd_h3
      integer          ::  T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer          ::  Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      integer          ::  P_l1, P_l2, P_l3, P_h1, P_h2, P_h3
      double precision :: rd(rd_l1:rd_h1, rd_l2:rd_h2, rd_l3:rd_h3, :)
      double precision ::  T( T_l1: T_h1,  T_l2: T_h2,  T_l3: T_h3)
      double precision ::  Y( Y_l1: Y_h1,  Y_l2: Y_h2,  Y_l3: Y_h3, :)
      double precision ::  P( P_l1: P_h1,  P_l2: P_h2,  P_l3: P_h3)

      ! Local variables
      integer          :: i, j, k, n
      double precision :: Yt(maxspec), Dt(maxspec)
      double precision :: scal, tscal, Wavg, RHO, Tt, Pt, invmwt(maxspec)
      double precision :: alpha, l1, l2, X(maxspec), CPMS(maxspec)

      call CKWT(IWRK(ckbi),RWRK(ckbr),invmwt)

      do n=1,Nspec
         invmwt(n) = 1.d0 / invmwt(n)
      end do

      scal  = 0.1d0
      tscal = 1.d-5
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,n)
               end do
               Tt = MAX(T(i,j),TMIN_TRANS)
               Pt = P(i,j)*1.d1
               CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
               CALL CKCPMS(Tt,IWRK(ckbi),RWRK(ckbr),CPMS)
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
               CALL EGSPAR(Tt,X,Yt,CPMS,RWRK(egbr),IWRK(egbi))
               CALL EGSV1(Pt,Tt,Yt,Wavg,RWRK(egbr),Dt)
               CALL CKRHOY(Pt,Tt,Yt,IWRK(ckbi),RWRK(ckbr),RHO)
               do n=1,Nspec
                  rd(i,j,n) = RHO * Wavg * invmwt(n) * Dt(n) * scal
               end do

               if (do_temp .ne. 0) then
                  alpha = 1
                  CALL EGSL1(alpha,Tt,X,RWRK(egbr),l1)
                  alpha = -1
                  CALL EGSL1(alpha,Tt,X,RWRK(egbr),l2)
                  rd(i,j,Nspec+1) = 0.5d0 * (l1 + l2) * tscal
               endif

               if (do_VelVisc .ne. 0) then
                  CALL EGSE3(Tt,Yt,RWRK(egbr),rd(i,j,Nspec+2))
                  rd(i,j,Nspec+2) = rd(i,j,Nspec+2) * scal
               endif

            end do
         end do
      end do

      end subroutine mixavg_rhodiff_temp_pres

end module diff_flux_module
