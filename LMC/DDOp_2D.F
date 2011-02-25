#include "DDOp_F.H"
#include "ChemDriver_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"

#define SDIM 2

      subroutine FORT_DDCOEFS(lo, hi, coefs, DIMS(coefs), T, DIMS(T), Y, DIMS(Y), Cp, DIMS(Cp), Full0_Mix1)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(coefs)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer DIMDEC(Cp)
      REAL_T coefs(DIMV(coefs),Ncoefs)
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y), Nspec)
      REAL_T Cp(DIMV(Cp),Nspec)
      integer Full0_Mix1
      REAL_T Ypt(maxspec), Cpi(maxspec), mwtInv(maxspec), Tpt
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
               Ypt(n) = Y(i,j,n)
               Cpi(n) = Cp(i,j,n)*1.d4
            enddo
            Tpt = T(i,j)
            call calc_coefs(Tpt,Ypt,Cpi,mwtInv,rhoDi,Full0_Mix1,PTC,rhoTD,rhoDij)
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
     &                       fluxYH, DIMS(fluxYH), F_dTe, DIMS(F_dTe),
     &                       YTc, DIMS(YTc), Xc, DIMS(Xc),
     &                       coefe, DIMS(coefe), Ae, DIMS(Ae),
     &                       for_T0_H1, Hic, DIMS(Hic), fillAlpha, alphaYH, DIMS(alphaYH),
     &                       Full0_Mix1)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM),dir
      integer DIMDEC(fluxYH)
      integer DIMDEC(F_dTe)
      integer DIMDEC(YTc)
      integer DIMDEC(Xc)
      integer DIMDEC(coefe)
      integer DIMDEC(Ae)
      integer DIMDEC(Hic)
      integer for_T0_H1, fillAlpha
      integer DIMDEC(alphaYH)
      integer Full0_Mix1
      REAL_T  delta(0:SDIM-1)
      REAL_T  fluxYH(DIMV(fluxYH),*)
      REAL_T   F_dTe(DIMV(F_dTe),*)
      REAL_T     YTc(DIMV(YTc),*)
      REAL_T      Xc(DIMV(Xc)  ,*)
      REAL_T      Ae(DIMV(Ae))
      REAL_T   coefe(DIMV(coefe) ,*)
      REAL_T     Hic(DIMV(Hic)   ,*)
      REAL_T  alphaYH(DIMV(alphaYH),*)

      integer i, j, n, m, idx
      REAL_T Te, Ye, Hie, dX, dY, mwtInv(maxspec)
      REAL_T RU, RUC, P1atm, dxInv
      REAL_T rhoTDe, rhoY_Dij, rDmixe, PTCe

      REAL_T, allocatable :: tMMwe(:,:)

      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RU = RU * 1.d-4
      call CKWT(IWRK(ckbi), RWRK(ckbr), mwtInv)
      do n=1,Nspec
         mwtInv(n) = 1.d0 / mwtInv(n)
      enddo

      dxInv = 1.d0/delta(dir)

      fluxYH(:,:,Nspec+1) = 0.0d0

      if (dir.eq.0) then         

c     Put Grad(T) in F_dTe slot 1 until after we have computed F
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               F_dTe(i,j,1) = ( YTc(i,j,Nspec+1) - YTc(i-1,j,Nspec+1) ) * dxInv
               PTCe = coefe(i,j,PTCoffset+1)
               fluxYH(i,j,Nspec+1) = -  PTCe * F_dTe(i,j,1)
            end do
         end do

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
                     rhoTDe = coefe(i,j,TDoffset+n)
                     dX = (  Xc(i,j,n) -  Xc(i-1,j,n)  ) * dxInv
                     Te = 0.5d0 * ( YTc(i,j,Nspec+1) + YTc(i-1,j,Nspec+1) )
                     fluxYH(i,j,n) = - Ye * rhoTDe * F_dTe(i,j,1) / Te
                     fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1) - (RU*Te*tMMWe(i,j)*rhoTDe*dX)
                  end do
               end do
            end do

            do n = 1,Nspec
               do m = 1,Nspec
                  idx = Doffset + (m-1)*Nspec + n
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)+1
                        rhoY_Dij = coefe(i,j,idx)
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
                     rDmixe = coefe(i,j,DMIXoffset+n)
                     fluxYH(i,j,n) = - rDmixe * dY
                  end do
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
c     Go backwards on n so that we have slot 1 until the end
            do n = Nspec,1,-1
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                     F_dTe(i,j,n) = fluxYH(i,j,n) * F_dTe(i,j,1)
                  end do
               end do
            enddo
         endif

         do n = 1,Nspec+1
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1
                  fluxYH(i,j,n) = fluxYH(i,j,n) * Ae(i,j)
               end do
            end do
         enddo

         if ( fillAlpha.eq.1 ) then
            do n = 1,Nspec
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                     rDmixe = coefe(i,j,DMIXoffset+n)
                     alphaYH(i-1,j,n) = alphaYH(i-1,j,n)  +  rDmixe * dxInv * Ae(i,j)
                     alphaYH(i  ,j,n) = alphaYH(i  ,j,n)  +  rDmixe * dxInv * Ae(i,j)
                  end do
               end do
            end do

            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1
                  PTCe = coefe(i,j,PTCoffset+1)
                  alphaYH(i-1,j,Nspec+1) = alphaYH(i-1,j,Nspec+1)  +  PTCe * dxInv * Ae(i,j)
                  alphaYH(i  ,j,Nspec+1) = alphaYH(i  ,j,Nspec+1)  +  PTCe * dxInv * Ae(i,j)
               end do
            end do
         end if

      else

c     Put Grad(T) in F_dTe slot 1 until after we have computed F
         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
               F_dTe(i,j,1) = ( YTc(i,j,Nspec+1) - YTc(i,j-1,Nspec+1) ) * dxInv
               PTCe = coefe(i,j,PTCoffset+1)
               fluxYH(i,j,Nspec+1) = PTCe * F_dTe(i,j,1)
            end do
         end do

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
                     rhoTDe = coefe(i,j,TDoffset+n)
                     dX = (  Xc(i,j,n) -  Xc(i,j-1,n)  ) * dxInv
                     Te = 0.5d0 * ( YTc(i,j,Nspec+1) + YTc(i,j-1,Nspec+1) )
                     fluxYH(i,j,n) = - Ye * rhoTDe * F_dTe(i,j,1) / Te
                     fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1) - (RU*Te*tMMWe(i,j)*rhoTDe*dX)
                  end do
               end do
            end do

            do n = 1,Nspec
               do m = 1,Nspec
                  idx = Doffset + (m-1)*Nspec + n
                  do j = lo(2), hi(2)+1
                     do i = lo(1), hi(1)
                        rhoY_Dij = coefe(i,j,idx)
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
                     rDmixe = coefe(i,j,DMIXoffset+n)
                     fluxYH(i,j,n) = - rDmixe * dY
                  end do
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
c     Go backwards on n so that we have slot 1 until the end
            do n = Nspec,1,-1
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     F_dTe(i,j,n) = fluxYH(i,j,n) * F_dTe(i,j,1)
                  end do
               end do
            enddo
         endif

         do n = 1,Nspec+1
            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)
                  fluxYH(i,j,n) = fluxYH(i,j,n) * Ae(i,j)
               end do
            end do
         enddo

         if ( fillAlpha.eq.1 ) then
            do n = 1,Nspec
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     rDmixe = coefe(i,j,DMIXoffset+n)
                     alphaYH(i,j-1,n) = alphaYH(i,j-1,n)  +  rDmixe * dxInv * Ae(i,j)
                     alphaYH(i  ,j,n) = alphaYH(i  ,j,n)  +  rDmixe * dxInv * Ae(i,j)
                  end do
               end do
            end do

            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)
                  PTCe = coefe(i,j,PTCoffset+1)
                  alphaYH(i,j-1,Nspec+1) = alphaYH(i,j-1,Nspec+1)  +  PTCe * dxInv * Ae(i,j)
                  alphaYH(i  ,j,Nspec+1) = alphaYH(i  ,j,Nspec+1)  +  PTCe * dxInv * Ae(i,j)
               end do
            end do
         end if

      endif
      end

      subroutine FORT_DDC2E(lo, hi, 
     &                      cfab, DIMS(cfab),
     &                      efab, DIMS(efab), nc, dir, avg_H0_A1)
      implicit none
      integer lo(SDIM), hi(SDIM), nc, dir, avg_H0_A1
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)

      integer i,j,n
      REAL_T den

      if (avg_H0_A1 .eq. 0) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                     den = cfab(i-1,j,n) + cfab(i,j,n)
                     if (den .ne. 0.d0) then
                        efab(i,j,n) = 2*(cfab(i-1,j,n) * cfab(i,j,n))/den
                     else
                        efab(i,j,n) = 0.d0
                     endif
                  enddo
               enddo
            enddo
         else
            do n = 1,nc
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     den = cfab(i,j-1,n) + cfab(i,j,n)
                     if (den .ne. 0.d0) then
                        efab(i,j,n) = 2*(cfab(i,j-1,n) * cfab(i,j,n))/den
                     else
                        efab(i,j,n) = 0.d0
                     endif
                  enddo
               enddo
            enddo
         endif
      else
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
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     efab(i,j,n) = half*(cfab(i,j-1,n) + cfab(i,j,n))
                  enddo
               enddo
            enddo
         endif
      endif
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

      subroutine FORT_DDECAVG(lo, hi, c, DIMS(c), f, DIMS(f), nc, ratio, dir, avg_H0_A1)
      implicit none
      integer lo(SDIM), hi(SDIM), nc, ratio(SDIM), dir, avg_H0_A1
      integer DIMDEC(c)
      integer DIMDEC(f)
      REAL_T c(DIMV(c),nc)
      REAL_T f(DIMV(f),nc)

      integer i,j,n,ii,jj,fi,fj
      REAL_T den, num, fac

      if (dir .eq. 0 ) then
         if (avg_H0_A1 .eq. 0) then 
            do n = 1, nc
               do j = lo(2), hi(2)
                  fj = ratio(2)*j
                  do i = lo(1), hi(1)
                     fi = ratio(1)*i
                     den = f(fi,fj,n)
                     do jj=1,ratio(2)-1
                        den = den + f(fi,fj+jj,n)
                     enddo
                     if (den .ne. 0.d0) then
                        num = ratio(2) * f(fi,fj,n)
                        do jj=1,ratio(2)-1
                           num = num * f(fi,fj+jj,n)
                        enddo
                        c(i,j,n) = num / den
                     else
                        c(i,j,n) = 0.d0
                     endif
                  enddo
               enddo
            enddo
         else
            fac = 1.d0/DBLE(ratio(2))
            do n = 1, nc
               do j = lo(2), hi(2)
                  fj = ratio(2)*j
                  do i = lo(1), hi(1)
                     fi = ratio(1)*i
                     c(i,j,n) = f(fi,fj,n)
                     do jj=1,ratio(2)-1
                        c(i,j,n) = c(i,j,n) + f(fi,fj+jj,n)
                     enddo
                     c(i,j,n) = c(i,j,n) * fac
                  enddo
               enddo
            enddo
         endif
      else
         if (avg_H0_A1 .eq. 0) then 
            do n = 1, nc
               do j = lo(2), hi(2)
                  fj = ratio(2)*j
                  do i = lo(1), hi(1)
                     fi = ratio(1)*i
                     den = ratio(1) * f(fi,fj,n)
                     do ii=1,ratio(1)-1
                        den = den + f(fi+ii,fj,n)
                     enddo
                     if (den .ne. 0.d0) then
                        num = f(fi,fj,n)
                        do ii=1,ratio(1)-1
                           num = num * f(fi+ii,fj,n)
                        enddo
                        c(i,j,n) = num / den
                     else
                        c(i,j,n) = 0.d0
                     endif
                  enddo
               enddo
            enddo
         else
            fac = 1.d0/DBLE(ratio(1))
            do n = 1, nc
               do j = lo(2), hi(2)
                  fj = ratio(2)*j
                  do i = lo(1), hi(1)
                     fi = ratio(1)*i
                     c(i,j,n) = f(fi,fj,n)
                     do ii=1,ratio(1)-1
                        c(i,j,n) = c(i,j,n) + f(fi+ii,fj,n)
                     enddo
                     c(i,j,n) = c(i,j,n) * fac
                  enddo
               enddo
            enddo
         endif
      endif
      end

      subroutine FORT_DDCCAVG(lo, hi, c, DIMS(c), f, DIMS(f), nc, ratio)
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
               do jj = 0,ratio(2)-1
                  do ii = 0,ratio(1)-1
                     c(i,j,n) = c(i,j,n) + f(ratio(1)*i+ii,ratio(2)*j+jj,n)
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

      integer i,j,n,ii,jj,fi,fj

c     NOTE: generalize to rz
      do n = 1, nc
         do j = lo(2), hi(2)
            fj = ratio(2)*j
            do i = lo(1), hi(1)
               fi = ratio(1)*i
               do jj = 0,ratio(2)-1
                  do ii = 0,ratio(1)-1
                     f(fi+ii,fj+jj,n) = c(i,j,n) + f(fi+ii,fj+jj,n)
                  enddo
               enddo
            enddo
         enddo
      enddo
      end
      
