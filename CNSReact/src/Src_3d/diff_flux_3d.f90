module diff_flux_module

  use cdwrk_module
  use meth_params_module
 
  implicit none
 
  public :: diffFlux
 
contains

      subroutine diffFlux(lo,hi, &
                          q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                          flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                          flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                          flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                          dx,dy,dz,dt,initFlux)

      ! Get diffusion flux on faces surrounding Box(lo,hi), requires 
      ! valid data on region Box(lo,hi).grow(1)

      implicit none

      integer lo(3),hi(3)
      integer :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer :: flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3
      integer :: flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3
      integer :: flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3
      double precision :: q(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3, QVAR)
      double precision :: flux1(flux1_l1:flux1_h1, flux1_l2:flux1_h2, flux1_l3:flux1_h3, NVAR)
      double precision :: flux2(flux2_l1:flux2_h1, flux2_l2:flux2_h2, flux2_l3:flux2_h3, NVAR)
      double precision :: flux3(flux3_l1:flux3_h1, flux3_l2:flux3_h2, flux3_l3:flux3_h3, NVAR)
      double precision :: dx,dy,dz,dt
      integer          :: initFlux

      ! Local variables
      integer          :: loD(3),hiD(3)
      double precision, allocatable ::    D(:,:,:,:)
      double precision, allocatable :: TEMP(:,:,:)
      double precision, allocatable ::   CP(:,:,:)

      allocate(   D(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,Nspec+3))
      allocate(TEMP(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3))
      allocate(  CP(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3))

      loD(:) = lo(:)-1
      hiD(:) = hi(:)+1

      if (initFlux.eq.1) then
         flux1(:,:,:,:) = 0.
         flux2(:,:,:,:) = 0.
         flux3(:,:,:,:) = 0.
      end if

      call transCoef(loD,hiD,q,TEMP,CP,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                     D,loD(1),loD(2),loD(3),hiD(1),hiD(2),hiD(3))

      call myViscFlux(lo,hi, &
                      q,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                      D,loD(1),loD(2),loD(3),hiD(1),hiD(2),hiD(3), &
                      flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                      flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                      flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                      dx,dy,dz,dt)

      deallocate(D,TEMP,CP)

      end subroutine diffFlux
! ::
! :: ----------------------------------------------------------
! ::

      subroutine transCoef(lo,hi,q,TEMP,CP, &
                           q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                           D, &
                           D_l1,D_l2,D_l3,D_h1,D_h2,D_h3)

      ! Note that lo,hi in this routine correspond to lo-1,hi+1 from the calling routine


      implicit none

      integer          :: lo(3), hi(3)
      integer          :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer          :: D_l1, D_l2, D_l3, D_h1, D_h2, D_h3
      double precision ::    q(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3,QVAR)
      double precision :: TEMP(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3)
      double precision ::   CP(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3) 
      double precision ::    D(D_l1:D_h1, D_l2:D_h2, D_l3:D_h3,Nspec+3)

      ! Local variables
      integer          :: i,j,k,doVelVisc,doTemp
!     double precision :: Dmax

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               TEMP(i,j,k) = MIN(MAX(q(i,j,k,QTHERM+1),300.d0),6500.d0)
            end do
         end do
      end do
      
      doTemp = 1
      doVelVisc = 1

      call mixavg_rhodiff_temp_pres(lo, hi, &
           D,  &
             D_l1,D_l2,D_l3,D_h1,D_h2,D_h3, &
           TEMP, &
             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
           q(q_l1,q_l2,q_l3,QTHERM+NADV+1), &
             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
           q(q_l1,q_l2,q_l3,QPRES), &
             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
           doTemp, doVelVisc)

      call FORT_CPMIXfromTY(lo, hi, &
           CP, &
             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
           TEMP, &
             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
           q(q_l1,q_l2,q_l3,QTHERM+NADV+1), &
             q_l1,q_l2,q_l3,q_h1,q_h2,q_h3)

      ! Replace lambda with lambda/Cp
      D(:,:,:,Nspec+1) = D(:,:,:,Nspec+1)/CP(:,:,:)

!     Dmax = -1.d0
!     do k = lo(3),hi(3)
!        do j = lo(2),hi(2)
!           do i = lo(1),hi(1)
!              do n=1,Nspec
!                 Dmax = MAX(Dmax,D(i,j,k,n) / q(i,j,k,QRHO))
!              enddo
!           end do
!        end do
!     end do
!     print*, '********************      Dmax: ', Dmax

      end subroutine transCoef
! ::
! :: ----------------------------------------------------------
! ::
      subroutine myViscFlux(lo,hi,&
                            q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                            D,D_l1,D_l2,D_l3,D_h1,D_h2,D_h3, &
                            flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                            flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                            flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                            dx,dy,dz,dt)


      implicit none

      integer          :: lo(3),hi(3)
      integer          :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer          :: D_l1, D_l2, D_l3, D_h1, D_h2, D_h3
      integer          :: flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3
      integer          :: flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3
      integer          :: flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3
      double precision :: q(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3, QVAR)
      double precision :: D(D_l1:D_h1, D_l2:D_h2, D_l3:D_h3, Nspec+3)
      double precision :: flux1(flux1_l1:flux1_h1, flux1_l2:flux1_h2, flux1_l3:flux1_h3, NVAR)
      double precision :: flux2(flux2_l1:flux2_h1, flux2_l2:flux2_h2, flux2_l3:flux2_h3, NVAR)
      double precision :: flux3(flux3_l1:flux3_h1, flux3_l2:flux3_h2, flux3_l3:flux3_h3, NVAR)

      double precision :: dx,dy,dz,dt

      ! Local variables
      integer          :: i,j,k,n,ifirstSp
      double precision :: tauxxm,tauxym,tauxzm
      double precision :: tauyym,tauyxm,tauyzm
      double precision :: tauzzm,tauzxm,tauzym
      double precision :: divxm,divym,divzm
      double precision :: muxm,muym,muzm
      double precision :: lamxm,lamym,lamzm
      double precision :: kxm,kym,kzm
      double precision :: phiflx

      double precision :: rhoDm(maxspec)
      double precision, parameter ::    two3rd = 2.d0/3.d0

      ifirstSp  = NTHERM+NADV+1

      ! compute x fluxes
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               
               muxm =  0.5d0*(D(i-1,j,k,Nspec+2) + D(i,j,k,Nspec+2))
               lamxm = -two3rd*muxm
               kxm =   0.5d0*(D(i-1,j,k,Nspec+1) + D(i,j,k,Nspec+1))
               do n=1,Nspec
                  rhoDm(n) = 0.5d0*(D(i-1,j,k,n) + D(i,j,k,n))
               end do
               
               tauxxm = 2.d0*muxm*(q(i,j,k,QU) - q(i-1,j,k,QU))/dx 
               
               divxm = lamxm*(q(i,j,k,QU) - q(i-1,j,k,QU))/dx  &
                    +.25d0*lamxm*(q(i-1,j+1,k,QV)+q(i,j+1,k,QV) &
                    -             q(i-1,j-1,k,QV)-q(i,j-1,k,QV))/dy &
                    +.25d0*lamxm*(q(i-1,j,k+1,QW)+q(i,j,k+1,QW) &
                    -             q(i-1,j,k-1,QW)-q(i,j,k-1,QW))/dz
               
               tauyxm =  muxm*(q(i,j,k,QV) - q(i-1,j,k,QV))/dx  &
                    +.25d0*muxm*(q(i,j+1,k,QU)+q(i-1,j+1,k,QU) &
                    -            q(i,j-1,k,QU)-q(i-1,j-1,k,QU))/dy
               
               tauzxm =  muxm*(q(i,j,k,QW) - q(i-1,j,k,QW))/dx &
                    +.25d0*muxm*(q(i,j,k+1,QU)+q(i-1,j,k+1,QU)&
                    -            q(i,j,k-1,QU)-q(i-1,j,k-1,QU))/dz

               phiflx =  tauxxm*(q(i,j,k,QU)+q(i-1,j,k,QU)) &
                    +     divxm*(q(i,j,k,QU)+q(i-1,j,k,QU)) &
                    +    tauyxm*(q(i,j,k,QV)+q(i-1,j,k,QV)) &
                    +    tauzxm*(q(i,j,k,QW)+q(i-1,j,k,QW))
               
               flux1(i,j,k,2) = flux1(i,j,k,2) - dt*dy*dz*(tauxxm+divxm)
               flux1(i,j,k,3) = flux1(i,j,k,3) - dt*dy*dz*tauyxm
               flux1(i,j,k,4) = flux1(i,j,k,4) - dt*dy*dz*tauzxm
               flux1(i,j,k,5) = flux1(i,j,k,5) - dt*dy*dz*(0.5d0*phiflx  &
                    + kxm*( (q(i  ,j,k,QREINT)+q(i  ,j,k,QPRES))/q(i  ,j,k,QRHO) &
                    -       (q(i-1,j,k,QREINT)+q(i-1,j,k,QPRES))/q(i-1,j,k,QRHO) )/dx)

               do n=1,Nspec
                  flux1(i,j,k,n+ifirstSp-1) = flux1(i,j,k,n+ifirstSp-1) &
                       - dt*dy*dz*rhoDm(n)* &
                       (q(i,j,k,QTHERM+NADV+n)-q(i-1,j,k,QTHERM+NADV+n))/dx
               end do
               
            enddo
         enddo
      enddo

      ! compute y fluxes
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)
               
               muym =  0.5d0*(D(i,j-1,k,Nspec+2) + D(i,j,k,Nspec+2))
               lamym = -two3rd*muym
               kym =   0.5d0*(D(i,j-1,k,Nspec+1) + D(i,j,k,Nspec+1))
               do n=1,Nspec
                  rhoDm(n) = 0.5d0*(D(i,j-1,k,n) + D(i,j,k,n))
               end do
               
               tauxym =  muym*(q(i,j,k,QU) - q(i,j-1,k,QU))/dy  &
                    +.25d0*muym*(q(i+1,j,k,QV)+q(i+1,j-1,k,QV) &
                    -            q(i-1,j,k,QV)-q(i-1,j-1,k,QV))/dx
               
               tauyym = 2.d0*muym*(q(i,j,k,QV) - q(i,j-1,k,QV))/dy 
               
               divym = lamym*(q(i,j,k,QV) - q(i,j-1,k,QV))/dy  &
                    +.25d0*lamym*(q(i+1,j-1,k,QU)+q(i+1,j,k,QU) &
                    -             q(i-1,j-1,k,QU)-q(i-1,j,k,QU))/dx &
                    +.25d0*lamym*(q(i,j-1,k+1,QW)+q(i,j,k+1,QW) &
                    -             q(i,j-1,k-1,QW)-q(i,j,k-1,QW))/dz
               
               tauzym =  muym*(q(i,j,k,QW) - q(i,j-1,k,QW))/dy  &
                    +.25d0*muym*(q(i,j,k+1,QV)+q(i,j-1,k+1,QV) &
                    -            q(i,j,k-1,QV)-q(i,j-1,k-1,QV))/dz

               phiflx = tauxym*(q(i,j,k,QU)+q(i,j-1,k,QU)) &
                    +   tauyym*(q(i,j,k,QV)+q(i,j-1,k,QV)) &
                    +    divym*(q(i,j,k,QV)+q(i,j-1,k,QV)) &
                    +   tauzym*(q(i,j,k,QW)+q(i,j-1,k,QW))
               
               flux2(i,j,k,2) = flux2(i,j,k,2) - dt*dx*dz*tauxym
               flux2(i,j,k,3) = flux2(i,j,k,3) - dt*dx*dz*(tauyym+divym)
               flux2(i,j,k,4) = flux2(i,j,k,4) - dt*dx*dz*tauzym
               flux2(i,j,k,5) = flux2(i,j,k,5) - dt*dy*dz*(0.5d0*phiflx  &
                    + 2.d0*kym*( (q(i,j  ,k,QREINT)+q(i,j  ,k,QPRES))/q(i,j  ,k,QRHO) &
                    -           (q(i,j-1,k,QREINT)+q(i,j-1,k,QPRES))/q(i,j-1,k,QRHO) )/dy)

               do n=1,Nspec
                  flux2(i,j,k,n+ifirstSp-1) = flux2(i,j,k,n+ifirstSp-1) &
                       - dt*dx*dz*rhoDm(n)* &
                       (q(i,j,k,QTHERM+NADV+n)-q(i,j-1,k,QTHERM+NADV+n))/dy
               end do
               
            enddo
         enddo
      enddo

      ! compute z fluxes
      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               
               muzm =  0.5d0*(D(i,j,k-1,Nspec+2) + D(i,j,k,Nspec+2))
               lamzm = -two3rd*muzm
               kzm =   0.5d0*(D(i,j,k-1,Nspec+1) + D(i,j,k,Nspec+1))
               do n=1,Nspec
                  rhoDm(n) = 0.5d0*(D(i,j,k-1,n) + D(i,j,k,n))
               end do
               
               tauxzm =  muzm*(q(i,j,k,QU) - q(i,j,k-1,QU))/dz  &
                    +.25d0*muzm*(q(i+1,j,k,QW)+q(i+1,j,k-1,QW) &
                    -            q(i-1,j,k,QW)-q(i-1,j,k-1,QW))/dx

               tauyzm =  muzm*(q(i,j,k,QV) - q(i,j,k-1,QV))/dz  &
                    +.25d0*muzm*(q(i,j+1,k,QW)+q(i,j+1,k-1,QW) &
                    -            q(i,j-1,k,QW)-q(i,j-1,k-1,QW))/dy
               
               tauzzm = 2.d0*muzm*(q(i,j,k,QW) - q(i,j,k-1,QW))/dz 

               divzm = lamzm*(q(i,j,k,QW) - q(i,j,k-1,QW))/dz  &
                    +.25d0*lamzm*(q(i+1,j,k-1,QU)+q(i+1,j,k,QU) &
                    -             q(i-1,j,k-1,QU)-q(i-1,j,k,QU))/dx &
                    +.25d0*lamzm*(q(i,j+1,k-1,QV)+q(i,j+1,k,QV) &
                    -             q(i,j-1,k-1,QV)-q(i,j-1,k,QV))/dy
               
               phiflx =tauxzm*(q(i,j,k,QU)+q(i,j,k-1,QU)) &
                    +  tauyzm*(q(i,j,k,QV)+q(i,j,k-1,QV)) &
                    +  tauzzm*(q(i,j,k,QW)+q(i,j,k-1,QW)) &
                    +   divzm*(q(i,j,k,QW)+q(i,j,k-1,QW))

               flux3(i,j,k,2) = flux3(i,j,k,2) - dt*dx*dy*tauxzm
               flux3(i,j,k,3) = flux3(i,j,k,3) - dt*dx*dy*tauyzm
               flux3(i,j,k,4) = flux3(i,j,k,4) - dt*dx*dy*(tauzzm+divzm)
               flux3(i,j,k,5) = flux3(i,j,k,5) - dt*dx*dy*(0.5d0*phiflx &
                    + 2.d0*kzm*( (q(i,j,k  ,QREINT)+q(i,j,k  ,QPRES))/q(i,j,k  ,QRHO) &
                    -            (q(i,j,k-1,QREINT)+q(i,j,k-1,QPRES))/q(i,j,k-1,QRHO) )/dz)

               do n=1,Nspec
                  flux3(i,j,k,n+ifirstSp -1) = flux3(i,j,k,n+ifirstSp -1) &
                       - dt*dx*dy*rhoDm(n)* &
                       (q(i,j,k,QTHERM+NADV+n)-q(i,j,k-1,QTHERM+NADV+n))/dz
               end do
               
            enddo
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
      double precision :: rd(rd_l1:rd_h1, rd_l2:rd_h2, rd_l3:rd_h3,*)
      double precision ::  T( T_l1: T_h1,  T_l2: T_h2,  T_l3: T_h3)
      double precision ::  Y( Y_l1: Y_h1,  Y_l2: Y_h2,  Y_l3: Y_h3,*)
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
                  Yt(n) = Y(i,j,k,n)
               end do
               Tt = MAX(T(i,j,k),TMIN_TRANS)
               Pt = P(i,j,k)*1.d1
               CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
               CALL CKCPMS(Tt,IWRK(ckbi),RWRK(ckbr),CPMS)
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
               CALL EGSPAR(Tt,X,Yt,CPMS,RWRK(egbr),IWRK(egbi))
               CALL EGSV1(Pt,Tt,Yt,Wavg,RWRK(egbr),Dt)
               CALL CKRHOY(Pt,Tt,Yt,IWRK(ckbi),RWRK(ckbr),RHO)
               do n=1,Nspec
                  rd(i,j,k,n) = RHO * Wavg * invmwt(n) * Dt(n) * scal
               end do

               if (do_temp .ne. 0) then
                  alpha = 1
                  CALL EGSL1(alpha,Tt,X,RWRK(egbr),l1)
                  alpha = -1
                  CALL EGSL1(alpha,Tt,X,RWRK(egbr),l2)
                  rd(i,j,k,Nspec+1) = 0.5d0 * (l1 + l2) * tscal
               endif

               if (do_VelVisc .ne. 0) then
                  CALL EGSE3(Tt,Yt,RWRK(egbr),rd(i,j,k,Nspec+2))
                  rd(i,j,k,Nspec+2) = rd(i,j,k,Nspec+2) * scal
               endif

            end do
         end do
      end do

      end subroutine mixavg_rhodiff_temp_pres

end module diff_flux_module
