#include "DDOp_F.H"
#include "ChemDriver_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"

#define SDIM 2

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
      integer i,j,n
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
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 1,Nspec
               Y(n) = S(i,j,n)
               Cpi(n) = Cp(i,j,n)*1.d4
            enddo
            T = S(i,j,Nspec+1)
            call calc_coefs(T,Y,Cpi,mwtInv,rhoDi,Full0_Mix1,PTC,rhoTD,rhoDij)
            do n = 1,Nspec
               coefs(i,j,DMIXoffset+n) = rhoDi(n)
            enddo
            coefs(i,j,PTCoffset+1) = PTC
            if (Full0_Mix1.eq.0) then
               do n = 1,Nspec
                  coefs(i,j,TDoffset+n) = rhoTD(n)
               enddo
               do n = 1,Nspec*Nspec
                  coefs(i,j,Doffset+n) = rhoDij(n)
               enddo
            endif            
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

      integer i, j, n, m, idx
      REAL_T Te, Ye, Hie, CPie, dX, dY, dT, mwtInv(maxspec)
      REAL_T RU, RUC, P1atm, dxInv
      REAL_T rhoTDe, rhoY_Dij, rDmixe, PTCe

      REAL_T, allocatable :: tMMwe(:,:)

      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RU = RU * 1.d-4
      call CKWT(IWRK(ckbi), RWRK(ckbr), mwtInv)
      do n=1,Nspec
         mwtInv(n) = 1.d0 / mwtInv(n)
      enddo

      dxInv = 1.d0/delta

      fluxYH(:,:,Nspec+1) = 0.0d0

      if (dir.eq.0) then         

         if ( Full0_Mix1 .eq.0 ) then
            allocate(tMMwe(lo(1):hi(1)+1,lo(2):hi(2)))

            tMMwe = 0.0d0

            do n = 1,Nspec
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                     Ye = ( YTc(i,j,n) + YTc(i-1,j,n)  ) * 0.5d0
                     tMMwe(i,j) = tMMwe(i,j) + Ye*mwtInv(n)
                  end do
               end do
            end do

            do n = 1,Nspec
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                     Ye = ( YTc(i,j,n) + YTc(i-1,j,n)  ) * 0.5d0

                     if ( Ye .le. 1.d-30 ) then
                        rhoTDe = 0.5d0*(coefc(i-1,j,TDoffset+n) + coefc(i,j  ,TDoffset+n))
                     else
                        rhoTDe = 0.5d0*(coefc(i-1,j,TDoffset+n)*YTc(i-1,j,n)
     &                       +          coefc(i,j  ,TDoffset+n)*YTc(i,j  ,n))/Ye
                     endif

                     dX = (  Xc(i,j,n) -  Xc(i-1,j,n)  ) * dxInv
                     dT = ( YTc(i,j,Nspec+1) - YTc(i-1,j,Nspec+1) ) * dxInv
                     Te = 0.5d0 * ( YTc(i,j,Nspec+1) + YTc(i-1,j,Nspec+1) )

                     fluxYH(i,j,n) = - Ye * rhoTDe * dT / Te
                     fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1) - (RU*Te*tMMWe(i,j)*rhoTDe*dX)
                  end do
               end do
            end do

            do n = 1,Nspec
               do m = 1,Nspec
                  idx = Doffset + (m-1)*Nspec + n
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)+1
                        rhoY_Dij = 0.5d0*(coefc(i-1,j,idx)+coefc(i,j,idx))
                        dX = (  Xc(i,j,m) -  Xc(i-1,j,m)  ) * dxInv
                        fluxYH(i,j,n) = fluxYH(i,j,n)  -  rhoY_Dij * dX
                     end do
                  end do
               end do
            end do

         end if

         if ( Full0_Mix1.eq.1 ) then 
            do n = 1,Nspec
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                     dY = ( YTc(i,j,n) - YTc(i-1,j,n)  ) * dxInv
                     rDmixe = 0.5d0*(coefc(i-1,j,DMIXoffset+n)+coefc(i,j,DMIXoffset+n))
                     fluxYH(i,j,n) = - rDmixe * dY
                  end do
               end do
            end do
         end if

         if ( fillAlpha.eq.1 ) then
            do n = 1,Nspec
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                     rDmixe = 0.5d0*(coefc(i-1,j,DMIXoffset+n)+coefc(i,j,DMIXoffset+n))
                     alphaYH(i-1,j,n) = alphaYH(i-1,j,n)  +  rDmixe * dxInv * Ae(i,j)
                     alphaYH(i  ,j,n) = alphaYH(i  ,j,n)  +  rDmixe * dxInv * Ae(i,j)
                  end do
               end do
            end do

            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1
                  PTCe = 0.5d0*(coefc(i-1,j,PTCoffset+1)+coefc(i,j,PTCoffset+1))
                  alphaYH(i-1,j,Nspec+1) = alphaYH(i-1,j,Nspec+1)  +  PTCe * dxInv * Ae(i,j)
                  alphaYH(i  ,j,Nspec+1) = alphaYH(i  ,j,Nspec+1)  +  PTCe * dxInv * Ae(i,j)
               end do
            end do
         end if

         if ( for_T0_H1.eq.1 ) then
            do n = 1,Nspec
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                     Hie = 0.5d0*(Hic(i-1,j,n)+Hic(i,j,n))
                     fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1)  +  Hie * fluxYH(i,j,n)
                  end do
               end do
            enddo
         else
            FcpDTe(lo(1):hi(1)+1,lo(2):hi(2)) = 0.0d0
            do n = 1,Nspec
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                     dT = ( YTc(i,j,Nspec+1) - YTc(i-1,j,Nspec+1) ) * dxInv
                     CPie = 0.5d0*(CPic(i-1,j,n)+CPic(i,j,n))
                     FcpDTe(i,j) = FcpDTe(i,j)  +  CPie * fluxYH(i,j,n) * dT
                  end do
               end do
            enddo
         endif

         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               dT = ( YTc(i,j,Nspec+1) - YTc(i-1,j,Nspec+1) ) * dxInv
               PTCe = 0.5d0*(coefc(i-1,j,PTCoffset+1)+coefc(i,j,PTCoffset+1))
               fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1)  -  PTCe * dT
            end do
         end do

         do n = 1,Nspec+1
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1
                  fluxYH(i,j,n) = fluxYH(i,j,n) * Ae(i,j)
               end do
            end do
         enddo

      else

         if ( Full0_Mix1 .eq.0 ) then
            allocate(tMMwe(lo(1):hi(1),lo(2):hi(2)+1))

            tMMwe = 0.0d0

            do n = 1,Nspec
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     Ye = ( YTc(i,j,n) + YTc(i,j-1,n)  ) * 0.5d0
                     tMMwe(i,j) = tMMwe(i,j) + Ye*mwtInv(n)
                  end do
               end do
            end do

            do n = 1,Nspec
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     Ye = ( YTc(i,j,n) + YTc(i,j-1,n)  ) * 0.5d0

                     if ( Ye.le.1.d-30 ) then
                        rhoTDe = 0.5d0*(coefc(i,j-1,TDoffset+n) + coefc(i,j  ,TDoffset+n))
                     else
                        rhoTDe = 0.5d0*(coefc(i,j-1,TDoffset+n)*YTc(i,j-1,n)
     &                       +          coefc(i,j  ,TDoffset+n)*YTc(i,j  ,n))/Ye
                     endif

                     dX = (  Xc(i,j,n) -  Xc(i,j-1,n)  ) * dxInv
                     dT = ( YTc(i,j,Nspec+1) - YTc(i,j-1,Nspec+1) ) * dxInv
                     Te = 0.5d0 * ( YTc(i,j,Nspec+1) + YTc(i,j-1,Nspec+1) )

                     fluxYH(i,j,n) = - Ye * rhoTDe * dT / Te
                     fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1) - (RU*Te*tMMWe(i,j)*rhoTDe*dX)
                  end do
               end do
            end do

            do n = 1,Nspec
               do m = 1,Nspec
                  idx = Doffset + (m-1)*Nspec + n
                  do j = lo(2), hi(2)+1
                     do i = lo(1), hi(1)
                        rhoY_Dij = 0.5d0*(coefc(i,j-1,idx)+coefc(i,j,idx))
                        dX = (  Xc(i,j,m) -  Xc(i,j-1,m)  ) * dxInv
                        fluxYH(i,j,n) = fluxYH(i,j,n)  -  rhoY_Dij * dX
                     end do
                  end do
               end do
            end do

         end if

         if ( Full0_Mix1.eq.1 ) then 
            do n = 1,Nspec
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     dY = ( YTc(i,j,n) - YTc(i,j-1,n)  ) * dxInv
                     rDmixe = 0.5d0*(coefc(i,j-1,DMIXoffset+n)+coefc(i,j,DMIXoffset+n))
                     fluxYH(i,j,n) = - rDmixe * dY
                  end do
               end do
            end do
         end if

         if ( fillAlpha.eq.1 ) then
            do n = 1,Nspec
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     rDmixe = 0.5d0*(coefc(i,j-1,DMIXoffset+n)+coefc(i,j,DMIXoffset+n))
                     alphaYH(i,j-1,n) = alphaYH(i,j-1,n)  +  rDmixe * dxInv * Ae(i,j)
                     alphaYH(i  ,j,n) = alphaYH(i  ,j,n)  +  rDmixe * dxInv * Ae(i,j)
                  end do
               end do
            end do

            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)
                  PTCe = 0.5d0*(coefc(i,j-1,PTCoffset+1)+coefc(i,j,PTCoffset+1))
                  alphaYH(i,j-1,Nspec+1) = alphaYH(i,j-1,Nspec+1)  +  PTCe * dxInv * Ae(i,j)
                  alphaYH(i  ,j,Nspec+1) = alphaYH(i  ,j,Nspec+1)  +  PTCe * dxInv * Ae(i,j)
               end do
            end do
         end if

         if ( for_T0_H1.eq.1 ) then
            do n = 1,Nspec
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     Hie = 0.5d0*(Hic(i,j-1,n)+Hic(i,j,n))
                     fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1)  +  Hie * fluxYH(i,j,n)
                  end do
               end do
            enddo
         else
            FcpDTe(lo(1):hi(1),lo(2):hi(2)+1) = 0.0d0
            do n = 1,Nspec
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     dT = ( YTc(i,j,Nspec+1) - YTc(i,j-1,Nspec+1) ) * dxInv
                     CPie = 0.5d0*(CPic(i,j-1,n)+CPic(i,j,n))
                     FcpDTe(i,j) = FcpDTe(i,j)  +  CPie * fluxYH(i,j,n) * dT
                  end do
               end do
            enddo
         endif

         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
               dT = ( YTc(i,j,Nspec+1) - YTc(i,j-1,Nspec+1) ) * dxInv
               PTCe = 0.5d0*(coefc(i,j-1,PTCoffset+1)+coefc(i,j,PTCoffset+1))
               fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1)  -  PTCe * dT
            end do
         end do

         do n = 1,Nspec+1
            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)
                  fluxYH(i,j,n) = fluxYH(i,j,n) * Ae(i,j)
               end do
            end do
         enddo

      endif

      end

      subroutine FORT_DDC2E(lo, hi, 
     &                      cfab, DIMS(cfab),
     &                      efab, DIMS(efab), nc, dir)
      implicit none
      integer lo(SDIM), hi(SDIM), nc, dir
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)

      integer i,j,n

      if (dir .EQ. 0) then
         do n = 1,nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1
                  efab(i,j,n) = half*(cfab(i-1,j,n) + cfab(i,j,n))
               enddo
            enddo
         enddo
      else
         do n = 1,nc
            do i = lo(1), hi(1)
               do j = lo(2), hi(2)+1
                  efab(i,j,n) = half*(cfab(i,j-1,n) + cfab(i,j,n))
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

      integer i,j,n

      if (d0a1.eq.0) then
         if (dir .EQ. 0) then
            do n= 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     div(i,j,n) = div(i,j,n)
     &                    +a*(V(i+1,j,n) - V(i,j,n))
                  enddo
               enddo               
            enddo
         else 
            do n= 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     div(i,j,n) = div(i,j,n)
     &                    + a*(V(i,j+1,n) - V(i,j,n))
                  enddo
               enddo               
            enddo
         endif
      else
         if (dir .EQ. 0) then
            do n= 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     div(i,j,n) = div(i,j,n)
     &                    + a*half*(V(i+1,j,n) + V(i,j,n))
                  enddo
               enddo               
            enddo
         else 
            do n= 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     div(i,j,n) = div(i,j,n)
     &                    + a*half*(V(i,j+1,n) + V(i,j,n))
                  enddo
               enddo               
            enddo
         endif
      endif
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

      integer i,j,n,ii,jj
      REAL_T fac

      fac = 1.d0 / (ratio(1)*ratio(2))

c     NOTE: generalize to rz
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j,n) = 0.d0
               do jj = 1,ratio(2)
                  do ii = 1,ratio(1)
                     c(i,j,n) = c(i,j,n) + f(2*i+ii-1,2*j+jj-1,n)
                  enddo
               enddo
               c(i,j,n) = c(i,j,n)*fac
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
      integer i,j,n,ii,jj

c     NOTE: generalize to rz
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do jj = 1,ratio(2)
                  do ii = 1,ratio(1)
                     f(2*i+ii-1,2*j+jj-1,n) = c(i,j,n) + f(2*i+ii-1,2*j+jj-1,n)
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      end
      
