#include "DDOp_F.H"
#include "ChemDriver_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"

#define SDIM 3

      subroutine FORT_DDCOEFS(lo, hi, coefs, DIMS(coefs), S, DIMS(S), Cp, DIMS(Cp), Full0_Mix1)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(coefs)
      integer DIMDEC(S)
      integer DIMDEC(Cp)
      REAL_T coefs(DIMV(coefs),Ncoefs)
      REAL_T S(DIMV(S), Nspec+1)
      REAL_T Cp(DIMV(Cp),Nspec)
      integer Full0_Mix1
      REAL_T Y(maxspec), Cpi(maxspec), mwtInv(maxspec), T
      integer i,j,k,n
      REAL_T rhoDi(maxspec),PTC,rhoTD(maxspec),rhoDij(maxspec*maxspec)

c     FIXME: These will never change...maybe should be set in common from the DDOp ctr
      call CKWT(IWRK(ckbi), RWRK(ckbr), mwtInv)
      do n=1,Nspec
         mwtInv(n) = 1.d0 / mwtInv(n)
      enddo

c     FIXME: perhaps calling the EGFlib or EGMlib routines would make this faster,
c     however, they do require the largest "pencil" size as part of the workspace setup
c     Personally, I hate dealing with that, so I'm going to punt and call point-by-point
c     Also, the actual code is in DDOp_F.F called by the _ND.F variants
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do n = 1,Nspec
                  Y(n) = S(i,j,k,n)
                  Cpi(n) = Cp(i,j,k,n)*1.d4
               enddo
               T = S(i,j,k,Nspec+1)
               call calc_coefs(T,Y,Cpi,mwtInv,rhoDi,Full0_Mix1,PTC,rhoTD,rhoDij)
               do n = 1,Nspec
                  coefs(i,j,k,DMIXoffset+n) = rhoDi(n)
               enddo
               coefs(i,j,k,PTCoffset+1) = PTC
               if (Full0_Mix1.eq.0) then
                  do n = 1,Nspec
                     coefs(i,j,k,TDoffset+n) = rhoTD(n)
                  enddo
                  do n = 1,Nspec*Nspec
                     coefs(i,j,k,Doffset+n) = rhoDij(n)
                  enddo
               endif            
            enddo
         enddo
      enddo
      end

      subroutine FORT_DDFLUX(lo, hi, delta, dir,
     &                       fluxYH, DIMS(fluxYH), FcpDTe, DIMS(FcpDTe),
     &                       YTc, DIMS(YTc), Xc, DIMS(Xc), coefc, DIMS(coefc),
     &                       CPic, DIMS(CPic), Ae, DIMS(Ae),
     &                       for_T0_H1, Hic, DIMS(Hic), fillAlpha, alphaYH, DIMS(alphaYH),
     &                       Full0_Mix1)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM),dir
      integer DIMDEC(fluxYH)
      integer DIMDEC(FcpDTe)
      integer DIMDEC(YTc)
      integer DIMDEC(Xc)
      integer DIMDEC(coefc)
      integer DIMDEC(CPic)
      integer DIMDEC(Ae)
      integer DIMDEC(Hic)
      integer for_T0_H1, fillAlpha
      integer DIMDEC(alphaYH)
      integer Full0_Mix1
      REAL_T  delta
      REAL_T  fluxYH(DIMV(fluxYH),*)
      REAL_T  FcpDTe(DIMV(FcpDTe))
      REAL_T     YTc(DIMV(YTc),*)
      REAL_T      Xc(DIMV(Xc)  ,*)
      REAL_T    CPic(DIMV(CPic) ,*)
      REAL_T      Ae(DIMV(Ae))
      REAL_T   coefc(DIMV(coefc) ,*)
      REAL_T     Hic(DIMV(Hic)   ,*)
      REAL_T  alphaYH(DIMV(alphaYH),*)

      integer i, j, k, n, m, idx
      REAL_T Te, Ye, Hie, CPie, dX, dY, dT, mwtInv(maxspec)
      REAL_T RU, RUC, P1atm, dxInv, MMWe
      REAL_T rhoTDe, rhoY_Dij, rDmixe, PTCe

      REAL_T, allocatable :: tMMwe(:,:,:)

      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RU = RU * 1.d-4
      call CKWT(IWRK(ckbi), RWRK(ckbr), mwtInv)
      do n=1,Nspec
         mwtInv(n) = 1.d0 / mwtInv(n)
      enddo

      dxInv = 1.d0/delta

      fluxYH(:,:,:,Nspec+1) = 0.0d0

      if ( dir.eq.0 ) then

         if ( Full0_Mix1 .eq.0 ) then
            allocate(tMMwe(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))

            tMMwe = 0.0d0

            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)+1
                        Ye = ( YTc(i,j,k,n) + YTc(i-1,j,k,n)  ) * 0.5d0
                        tMMwe(i,j,k) = tMMwe(i,j,k) + Ye*mwtInv(n)
                     end do
                  end do
               end do
            end do

            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)+1
                        Ye = ( YTc(i,j,k,n) + YTc(i-1,j,k,n)  ) * 0.5d0

                        if ( Ye .le. 1.d-30 ) then
                           rhoTDe = 0.5d0*(coefc(i-1,j,k,TDoffset+n) + coefc(i,j,k,TDoffset+n))
                        else
                           rhoTDe = 0.5d0*(coefc(i-1,j,k,TDoffset+n)*YTc(i-1,j,k,n)
     &                          +          coefc(i,j  ,k,TDoffset+n)*YTc(i  ,j,k,n))/Ye
                        endif

                        dX = (  Xc(i,j,k,n) -  Xc(i-1,j,k,n)  ) * dxInv
                        dT    = ( YTc(i,j,k,Nspec+1) - YTc(i-1,j,k,Nspec+1) ) * dxInv
                        Te    = 0.5d0 * ( YTc(i,j,k,Nspec+1) + YTc(i-1,j,k,Nspec+1) )
                        fluxYH(i,j,k,n) = - Ye * rhoTDe * dT / Te
                        fluxYH(i,j,k,Nspec+1) = fluxYH(i,j,k,Nspec+1) - (RU*Te*tMMwe(i,j,k))*rhoTDe*dX
                     end do
                  end do
               end do
            end do

            do n = 1,Nspec
               do m = 1,Nspec
                  idx = Doffset + (m-1)*Nspec + n
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        do i = lo(1), hi(1)+1
                           rhoY_Dij = 0.5d0*(coefc(i-1,j,k,idx)+coefc(i,j,k,idx))
                           dX = (  Xc(i,j,k,m) -  Xc(i-1,j,k,m)  ) * dxInv
                           fluxYH(i,j,k,n) = fluxYH(i,j,k,n) - rhoY_Dij * dX
                        end do
                     end do
                  end do
               end do
            end do

         end if

         if ( Full0_Mix1.eq.1 ) then 
            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)+1
                        dY = ( YTc(i,j,k,n) - YTc(i-1,j,k,n)  ) * dxInv
                        rDmixe = 0.5d0*(coefc(i-1,j,k,DMIXoffset+n)+coefc(i,j,k,DMIXoffset+n))
                        fluxYH(i,j,k,n) = - rDmixe * dY
                     end do
                  end do
               end do
            end do
         end if

         if ( fillAlpha.eq.1 ) then
            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)+1
                        rDmixe = 0.5d0*(coefc(i-1,j,k,DMIXoffset+n)+coefc(i,j,k,DMIXoffset+n))
                        alphaYH(i-1,j,k,n) = alphaYH(i-1,j,k,n)  +  rDmixe * dxInv * Ae(i,j,k)
                        alphaYH(i  ,j,k,n) = alphaYH(i  ,j,k,n)  +  rDmixe * dxInv * Ae(i,j,k)
                     end do
                  end do
               end do
            end do

            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                     PTCe = 0.5d0*(coefc(i-1,j,k,PTCoffset+1)+coefc(i,j,k,PTCoffset+1))
                     alphaYH(i-1,j,k,Nspec+1) = alphaYH(i-1,j,k,Nspec+1)  +  PTCe * dxInv * Ae(i,j,k)
                     alphaYH(i  ,j,k,Nspec+1) = alphaYH(i  ,j,k,Nspec+1)  +  PTCe * dxInv * Ae(i,j,k)
                  end do
               end do
            end do
         end if

         if ( for_T0_H1.eq.1 ) then
            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)+1
                        Hie = 0.5d0*(Hic(i-1,j,k,n)+Hic(i,j,k,n))
                        fluxYH(i,j,k,Nspec+1) = fluxYH(i,j,k,Nspec+1) + Hie * fluxYH(i,j,k,n)
                     end do
                  end do
               end do
            enddo
         else
            FcpDTe(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = 0.0d0

            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)+1
                        dT = ( YTc(i,j,k,Nspec+1) - YTc(i-1,j,k,Nspec+1) ) * dxInv
                        CPie = 0.5d0*(CPic(i-1,j,k,n)+CPic(i,j,k,n))
                        FcpDTe(i,j,k) = FcpDTe(i,j,k)  +  CPie * fluxYH(i,j,k,n) * dT
                     end do
                  end do
               end do
            enddo
         endif

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1
                  dT = ( YTc(i,j,k,Nspec+1) - YTc(i-1,j,k,Nspec+1) ) * dxInv
                  PTCe = 0.5d0*(coefc(i-1,j,k,PTCoffset+1)+coefc(i,j,k,PTCoffset+1))
                  fluxYH(i,j,k,Nspec+1) = fluxYH(i,j,k,Nspec+1)  -  PTCe * dT
               end do
            end do
         end do

         do n = 1,Nspec+1
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                     fluxYH(i,j,k,n) = fluxYH(i,j,k,n) * Ae(i,j,k)
                  end do
               end do
            end do
         enddo

      else if (dir.eq.1) then

         if ( Full0_Mix1 .eq.0 ) then
            allocate(tMMwe(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))

            tMMwe = 0.0d0

            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)+1
                     do i = lo(1), hi(1)
                        Ye = ( YTc(i,j,k,n) + YTc(i,j-1,k,n)  ) * 0.5d0
                        tMMwe(i,j,k) = tMMwe(i,j,k) + Ye*mwtInv(n)
                     end do
                  end do
               end do
            end do

            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)+1
                     do i = lo(1), hi(1)
                        Ye = ( YTc(i,j,k,n) + YTc(i,j-1,k,n)  ) * 0.5d0

                        if ( Ye .le. 1.d-30 ) then
                           rhoTDe = 0.5d0*(coefc(i,j-1,k,TDoffset+n) + coefc(i,j,k,TDoffset+n))
                        else
                           rhoTDe = 0.5d0*(coefc(i,j-1,k,TDoffset+n)*YTc(i,j-1,k,n)
     &                          +          coefc(i,j  ,k,TDoffset+n)*YTc(i,j  ,k,n))/Ye
                        endif

                        dX = (  Xc(i,j,k,n) -  Xc(i,j-1,k,n)  ) * dxInv
                        dT = ( YTc(i,j,k,Nspec+1) - YTc(i,j-1,k,Nspec+1) ) * dxInv
                        Te = 0.5d0 * ( YTc(i,j,k,Nspec+1) + YTc(i,j-1,k,Nspec+1) )
                        fluxYH(i,j,k,n) = - Ye * rhoTDe * dT / Te
                        fluxYH(i,j,k,Nspec+1) = fluxYH(i,j,k,Nspec+1) - (RU*Te*tMMwe(i,j,k))*rhoTDe*dX
                     end do
                  end do
               end do
            end do

            do n = 1,Nspec
               do m = 1,Nspec
                  idx = Doffset + (m-1)*Nspec + n
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)+1
                        do i = lo(1), hi(1)
                           rhoY_Dij = 0.5d0*(coefc(i,j-1,k,idx)+coefc(i,j,k,idx))
                           dX = (  Xc(i,j,k,m) -  Xc(i,j-1,k,m)  ) * dxInv
                           fluxYH(i,j,k,n) = fluxYH(i,j,k,n)  -  rhoY_Dij * dX
                        end do
                     end do
                  end do
               end do
            end do
         end if

         if ( Full0_Mix1.eq.1 ) then 
            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)+1
                     do i = lo(1), hi(1)
                        rDmixe = 0.5d0*(coefc(i,j-1,k,DMIXoffset+n)+coefc(i,j,k,DMIXoffset+n))
                        dY = ( YTc(i,j,k,n) - YTc(i,j-1,k,n)  ) * dxInv
                        fluxYH(i,j,k,n) = - rDmixe * dY
                     end do
                  end do
               end do
            end do
         end if

         if ( fillAlpha.eq.1 ) then
            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)+1
                     do i = lo(1), hi(1)
                        rDmixe = 0.5d0*(coefc(i,j-1,k,DMIXoffset+n)+coefc(i,j,k,DMIXoffset+n))
                        alphaYH(i,j-1,k,n) = alphaYH(i,j-1,k,n)  +  rDmixe * dxInv * Ae(i,j,k)
                        alphaYH(i  ,j,k,n) = alphaYH(i  ,j,k,n)  +  rDmixe * dxInv * Ae(i,j,k)
                     end do
                  end do
               end do
            end do

            do k = lo(3), hi(3)
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     PTCe = 0.5d0*(coefc(i,j-1,k,PTCoffset+1)+coefc(i,j,k,PTCoffset+1))
                     alphaYH(i,j-1,k,Nspec+1) = alphaYH(i,j-1,k,Nspec+1)  +  PTCe * dxInv * Ae(i,j,k)
                     alphaYH(i  ,j,k,Nspec+1) = alphaYH(i  ,j,k,Nspec+1)  +  PTCe * dxInv * Ae(i,j,k)
                  end do
               end do
            end do
         end if

         if ( for_T0_H1.eq.1 ) then
            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)+1
                     do i = lo(1), hi(1)
                        Hie = 0.5d0*(Hic(i,j-1,k,n)+Hic(i,j,k,n))
                        fluxYH(i,j,k,Nspec+1) = fluxYH(i,j,k,Nspec+1)  +  Hie * fluxYH(i,j,k,n)
                     end do
                  end do
               end do
            enddo
         else
            FcpDTe(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)) = 0.0d0

            do n = 1,Nspec
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)+1
                     do i = lo(1), hi(1)
                        dT = ( YTc(i,j,k,Nspec+1) - YTc(i,j-1,k,Nspec+1) ) * dxInv
                        CPie = 0.5d0*(CPic(i,j-1,k,n)+CPic(i,j,k,n))
                        FcpDTe(i,j,k) = FcpDTe(i,j,k)  +  CPie * fluxYH(i,j,k,n) * dT
                     end do
                  end do
               end do
            enddo
         endif

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)
                  dT = ( YTc(i,j,k,Nspec+1) - YTc(i,j-1,k,Nspec+1) ) * dxInv
                  PTCe = 0.5d0*(coefc(i,j-1,k,PTCoffset+1)+coefc(i,j,k,PTCoffset+1))
                  fluxYH(i,j,k,Nspec+1) = fluxYH(i,j,k,Nspec+1)  -  PTCe * dT
               end do
            end do
         end do

         do n = 1,Nspec+1
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     fluxYH(i,j,k,n) = fluxYH(i,j,k,n) * Ae(i,j,k)
                  end do
               end do
            end do
         enddo

      else if (dir.eq.2) then

         if ( Full0_Mix1 .eq.0 ) then
            allocate(tMMwe(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))

            tMMwe = 0.0d0

            do n = 1,Nspec
               do k = lo(3), hi(3)+1
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        Ye = ( YTc(i,j,k,n) + YTc(i,j,k-1,n)  ) * 0.5d0
                        tMMwe(i,j,k) = tMMwe(i,j,k) + Ye*mwtInv(n)
                     end do
                  end do
               end do
            end do

            do n = 1,Nspec
               do k = lo(3), hi(3)+1
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        Ye = ( YTc(i,j,k,n) + YTc(i,j,k-1,n)  ) * 0.5d0

                        if ( Ye.le.1.d-30 ) then
                           rhoTDe = 0.5d0*(coefc(i,j,k-1,TDoffset+n) + coefc(i,j,k,TDoffset+n))
                        else
                           rhoTDe = 0.5d0*(coefc(i,j,k-1,TDoffset+n)*YTc(i,j,k-1,n)
     &                          +          coefc(i,j  ,k,TDoffset+n)*YTc(i,j  ,k,n))/Ye
                        endif

                        dX = (  Xc(i,j,k,n) -  Xc(i,j,k-1,n)  ) * dxInv
                        dT = ( YTc(i,j,k,Nspec+1) - YTc(i,j,k-1,Nspec+1) ) * dxInv
                        Te = 0.5d0 * ( YTc(i,j,k,Nspec+1) + YTc(i,j,k-1,Nspec+1) )

                        fluxYH(i,j,k,n) = - Ye * rhoTDe * dT / Te
                        fluxYH(i,j,k,Nspec+1) = fluxYH(i,j,k,Nspec+1) - (RU*Te*tMMwe(i,j,k))*rhoTDe*dX
                     end do
                  end do
               end do
            end do

            do n = 1,Nspec
               do m = 1,Nspec
                  idx = Doffset + (m-1)*Nspec + n
                  do k = lo(3), hi(3)+1
                     do j = lo(2), hi(2)
                        do i = lo(1), hi(1)
                           rhoY_Dij = 0.5d0*(coefc(i,j,k-1,idx)+coefc(i,j,k,idx))
                           dX = (  Xc(i,j,k,m) -  Xc(i,j,k-1,m)  ) * dxInv
                           fluxYH(i,j,k,n) = fluxYH(i,j,k,n)  -  rhoY_Dij * dX
                        end do
                     end do
                  end do
               end do
            end do
         end if

         if ( Full0_Mix1.eq.1 ) then 
            do n = 1,Nspec
               do k = lo(3), hi(3)+1
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        rDmixe = 0.5d0*(coefc(i,j,k-1,DMIXoffset+n)+coefc(i,j,k,DMIXoffset+n))
                        dY = ( YTc(i,j,k,n) - YTc(i,j,k-1,n)  ) * dxInv
                        fluxYH(i,j,k,n) = - rDmixe * dY
                     end do
                  end do
               end do
            end do
         end if

         if ( fillAlpha.eq.1 ) then
            do n = 1,Nspec
               do k = lo(3), hi(3)+1
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        rDmixe = 0.5d0*(coefc(i,j,k-1,DMIXoffset+n)+coefc(i,j,k,DMIXoffset+n))
                        alphaYH(i,j,k-1,n) = alphaYH(i,j,k-1,n)  +  rDmixe * dxInv * Ae(i,j,k)
                        alphaYH(i,j,k  ,n) = alphaYH(i,j,k  ,n)  +  rDmixe * dxInv * Ae(i,j,k)
                     end do
                  end do
               end do
            end do

            do k = lo(3), hi(3)+1
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     PTCe = 0.5d0*(coefc(i,j,k-1,PTCoffset+1)+coefc(i,j,k,PTCoffset+1))
                     alphaYH(i,j,k-1,Nspec+1) = alphaYH(i,j,k-1,Nspec+1)  +  PTCe * dxInv * Ae(i,j,k)
                     alphaYH(i,j,k  ,Nspec+1) = alphaYH(i,j,k  ,Nspec+1)  +  PTCe * dxInv * Ae(i,j,k)
                  end do
               end do
            end do
         end if

         if ( for_T0_H1.eq.1 ) then
            do n = 1,Nspec
               do k = lo(3), hi(3)+1
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        Hie = 0.5d0*(Hic(i,j,k-1,n)+Hic(i,j,k,n))
                        fluxYH(i,j,k,Nspec+1) = fluxYH(i,j,k,Nspec+1)  +  Hie * fluxYH(i,j,k,n)
                     end do
                  end do
               end do
            enddo
         else
            FcpDTe(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1) = 0.0d0

            do n = 1,Nspec
               do k = lo(3), hi(3)+1
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        dT = ( YTc(i,j,k,Nspec+1) - YTc(i,j,k-1,Nspec+1) ) * dxInv
                        CPie = 0.5d0*(CPic(i,j,k-1,n)+CPic(i,j,k,n))
                        FcpDTe(i,j,k) = FcpDTe(i,j,k)  +  CPie * fluxYH(i,j,k,n) * dT
                     end do
                  end do
               end do
            enddo
         endif

         do k = lo(3), hi(3)+1
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  dT = ( YTc(i,j,k,Nspec+1) - YTc(i,j,k-1,Nspec+1) ) * dxInv
                  PTCe = 0.5d0*(coefc(i,j,k-1,PTCoffset+1)+coefc(i,j,k,PTCoffset+1))
                  fluxYH(i,j,k,Nspec+1) = fluxYH(i,j,k,Nspec+1)  -  PTCe * dT
               end do
            end do
         end do

         do n = 1,Nspec+1
            do k = lo(3), hi(3)+1
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     fluxYH(i,j,k,n) = fluxYH(i,j,k,n) * Ae(i,j,k)
                  end do
               end do
            end do
         enddo

      endif
      end


      subroutine FORT_DDC2E(lo, hi, 
     &                      DIMS(cfab), cfab,
     &                      DIMS(efab), efab, nc, dir)
      implicit none
      integer lo(SDIM), hi(SDIM), nc, dir
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)

      integer i,j,k,n

c     Here, the cc data is assumed to be computable over the box (lo,hi) 
c     surroundingNodes of grow((lo,hi),-1) in the dir direction.
      if (dir .EQ. 0) then
         do n = 1,nc
            do k = lo(3)+1, hi(3)-1
               do j = lo(2)+1, hi(2)-1
                  do i = lo(1)+1, hi(1)
                     efab(i,j,k,n) = half*(cfab(i-1,j,k,n) + cfab(i,j,k,n))
                  enddo
               enddo
            enddo
         enddo
      else if (dir .EQ. 1) then
         do n = 1,nc
            do k = lo(3)+1, hi(3)-1
               do j = lo(2)+1, hi(2)
                  do i = lo(1)+1, hi(1)-1
                     efab(i,j,k,n) = half*(cfab(i,j-1,k,n) + cfab(i,j,k,n))
                  enddo
               enddo
            enddo
         enddo
      else
         do n = 1,nc
            do k = lo(3)+1, hi(3)
               do j = lo(2)+1, hi(2)-1
                  do i = lo(1)+1, hi(1)-1
                     efab(i,j,k,n) = half*(cfab(i,j,k-1,n) + cfab(i,j,k,n))
                  enddo
               enddo
            enddo
         enddo
      end if
      end

      subroutine FORT_DDETC(lo, hi, div, DIMS(div), V, DIMS(V),
     &                      a, dir, nc, d0a1);
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM), dir, nc, d0a1
      integer DIMDEC(div)
      integer DIMDEC(V)
      REAL_T div(DIMV(div) ,*)
      REAL_T   V(DIMV(V)   ,*)
      REAL_T a

      integer i,j,k,n

      if (d0a1.eq.0) then
         if (dir .EQ. 0) then
            do n= 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        div(i,j,k,n) = div(i,j,k,n)
     &                       + a*(V(i+1,j,k,n) - V(i,j,k,n))
                     enddo
                  enddo
               enddo               
            enddo
         elseif (dir .EQ. 1) then
            do n= 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        div(i,j,k,n) = div(i,j,k,n)
     &                       + a*(V(i,j+1,k,n) - V(i,j,k,n))
                     enddo
                  enddo
               enddo               
            enddo
         else 
            do n= 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        div(i,j,k,n) = div(i,j,k,n)
     &                       + a*(V(i,j,k+1,n) - V(i,j,k,n))
                     enddo
                  enddo
               enddo               
            enddo
         end if
      else
         if (dir .EQ. 0) then
            do n= 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        div(i,j,k,n) = div(i,j,k,n)
     &                       + a*half*(V(i+1,j,k,n) + V(i,j,k,n))
                     enddo
                  enddo
               enddo               
            enddo
         elseif (dir .EQ. 1) then
            do n= 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        div(i,j,k,n) = div(i,j,k,n)
     &                       + a*half*(V(i,j+1,k,n) + V(i,j,k,n))
                     enddo
                  enddo
               enddo               
            enddo
         else 
            do n= 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        div(i,j,k,n) = div(i,j,k,n)
     &                       + a*half*(V(i,j,k+1,n) + V(i,j,k,n))
                     enddo
                  enddo
               enddo               
            enddo
         endif
      endif
      end

      subroutine FORT_CRSNCCBND(lo, hi, fine, DIMS(fine),
     &                          crse, DIMS(crse), nc, face, ratio)
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(fine)
      integer DIMDEC(crse)
      REAL_T fine(DIMV(fine),nc)
      REAL_T crse(DIMV(crse),nc)
      integer face, ratio(SDIM)
      
      integer i,j,k,n,ii,jj,kk, clo1, chi1, clo2, chi2, clo3, chi3, offset, dir
      
      if ((face.EQ.0).OR.(face.EQ.3)) then
         dir = 0
      elseif ((face.EQ.1).OR.(face.EQ.4)) then
         dir = 1
      else
         dir = 2
      end if

      if (face .LT. SDIM) then
         offset = 1
      else
         offset = 0
      end if
      
      if (dir .eq. 0) then
         clo1 = (lo(1)-offset)/2
         clo2 = lo(2)/2
         clo3 = lo(3)/2
         chi1 = (hi(1)-offset)/2
         chi2 = (hi(2)+1)/2 - 1
         chi3 = (hi(3)+1)/2 - 1
         do k=clo3,chi3
            do j=clo2,chi2
               do i=clo1,chi1
                  do n=1,nc
                     crse(i,j,k,n) = 0.d0
                     do jj=1,ratio(2)
                        do kk=1,ratio(3)
                           crse(i,j,k,n) = crse(i,j,k,n) + fine(2*i+offset,2*j+jj-1,2*k+kk-1,n)
                        enddo
                     enddo
                     crse(i,j,k,n) = crse(i,j,k,n) / (ratio(2)*ratio(3))
                  enddo
               enddo
            enddo
         enddo
      elseif (dir .eq. 1) then
         clo1 = lo(1)/2
         clo2 = (lo(2)-offset)/2
         clo3 = lo(3)/2
         chi1 = (hi(1)+1)/2 - 1
         chi2 = (hi(2)-offset)/2
         chi3 = (hi(3)+1)/2 - 1
         do k=clo3,chi3
            do j=clo2,chi2
               do i=clo1,chi1
                  do n=1,nc
                     crse(i,j,k,n) = 0.d0
                     do ii=1,ratio(1)
                        do kk=1,ratio(3)
                           crse(i,j,k,n) = crse(i,j,k,n) + fine(2*i+ii-1,2*j+offset,2*k+kk-1,n)
                        enddo
                     enddo
                     crse(i,j,k,n) = crse(i,j,k,n) / (ratio(1)*ratio(3))
                  enddo
               enddo
            enddo
         enddo
      else
         clo1 = lo(1)/2
         clo2 = lo(2)/2
         clo3 = (lo(3)-offset)/2
         chi1 = (hi(1)+1)/2 - 1
         chi2 = (hi(2)+1)/2 - 1
         chi3 = (hi(3)-offset)/2
         do k=clo3,chi3
            do j=clo2,chi2
               do i=clo1,chi1
                  do n=1,nc
                     crse(i,j,k,n) = 0.d0
                     do ii=1,ratio(1)
                        do jj=1,ratio(2)
                           crse(i,j,k,n) = crse(i,j,k,n) + fine(2*i+ii-1,2*j+jj-1,2*k+offset,n)
                        enddo
                     enddo
                     crse(i,j,k,n) = crse(i,j,k,n) / (ratio(1)*ratio(2))
                  enddo
               enddo
            enddo
         enddo
      end if
      end

      subroutine FORT_DDCCAVG(c, DIMS(c), f, DIMS(f), lo, hi, nc, ratio)
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
      integer ratio(SDIM)
c
      integer i,j,k,n,ii,jj,kk
      REAL_T fac

      fac = 1.d0/(ratio(1)*ratio(2)*ratio(3))
c     NOTE: generalize to rz
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  c(i,j,k,n) = 0.d0
                  do kk = 1,ratio(3)
                     do jj = 1,ratio(2)
                         do ii = 1,ratio(1)
                            c(i,j,k,n) = c(i,j,k,n) + f(2*i+ii-1,2*j+jj-1,2*k+kk-1,n)
                         enddo
                     enddo
                  enddo
                  c(i,j,k,n) = c(i,j,k,n)*fac
               enddo
            enddo
         enddo
      enddo
      end

      subroutine FORT_DDCCINT(f, DIMS(f), c, DIMS(c), lo, hi, nc, ratio)
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
      integer ratio(SDIM)
c
      integer i,j,k,n,ii,jj,kk

c     NOTE: generalize to rz
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do kk = 1,ratio(3)
                     do jj = 1,ratio(2)
                        do ii = 1,ratio(1)
                           f(2*i+ii-1,2*j+jj-1,2*k+kk-1,n) = c(i,j,k,n) + f(2*i+ii-1,2*j+jj-1,2*k+kk-1,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      end
      
