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
     &                       YTc, DIMS(YTc), Xc, DIMS(Xc), coefc, DIMS(coefc), CPic, DIMS(CPic),
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
      REAL_T   coefc(DIMV(coefc) ,*)
      REAL_T     Hic(DIMV(Hic)   ,*)
      REAL_T  alphaYH(DIMV(alphaYH),*)

      integer i, j, n, m, idx
      REAL_T Te, Ye(maxspec), Hie, CPie, dX(maxspec), dY(maxspec), dT, mwtInv(maxspec)
      REAL_T RU, RUC, P1atm, dxInv, dx2Inv, MMWe
      REAL_T rhoTDe, rhoY_Dij, rDmixe, PTCe

      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RU = RU * 1.d-4
      call CKWT(IWRK(ckbi), RWRK(ckbr), mwtInv)
      do n=1,Nspec
         mwtInv(n) = 1.d0 / mwtInv(n)
      enddo

      dxInv = 1.d0/delta
      dx2Inv = dxInv*dxInv

      if (dir.eq.0) then         
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               do n = 1,Nspec
                  Ye(n) = ( YTc(i,j,n) + YTc(i-1,j,n)  ) * 0.5d0
                  dX(n) = (  Xc(i,j,n) -  Xc(i-1,j,n)  ) * dxInv
                  dY(n) = ( YTc(i,j,n) - YTc(i-1,j,n)  ) * dxInv
               enddo
               dT = ( YTc(i,j,Nspec+1) - YTc(i-1,j,Nspec+1) ) * dxInv
               Te = 0.5d0 * ( YTc(i,j,Nspec+1) + YTc(i-1,j,Nspec+1) )
               call CKMMWY(Ye, IWRK(ckbi), RWRK(ckbr), MMWe)

               fluxYH(i,j,Nspec+1) = 0.d0
               if (Full0_Mix1 .eq.0) then
                  do n = 1,Nspec
                     if (Ye(n).le.1.d-30) then
                        rhoTDe = 0.5d0*(coefc(i-1,j,TDoffset+n) + coefc(i,j  ,TDoffset+n))
                     else
                        rhoTDe = 0.5d0*(coefc(i-1,j,TDoffset+n)*YTc(i-1,j,n)
     &                       +          coefc(i,j  ,TDoffset+n)*YTc(i,j  ,n))/Ye(n)
                     endif
                     fluxYH(i,j,n) = - Ye(n) * rhoTDe * dT / Te
c     Note: Forming D_(n,m) . Grad(X_m), and assuming D was loaded with n varying fastest
                     do m = 1,Nspec
                        idx = Doffset + (m-1)*Nspec + n
                        rhoY_Dij = 0.5d0*(coefc(i-1,j,idx)+coefc(i,j,idx))
                        fluxYH(i,j,n) = fluxYH(i,j,n)  -  rhoY_Dij * dX(m)
                     enddo
                     fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1) - (RU*Te/MMWe)*rhoTDe*dX(n)
                  enddo
               endif
               do n = 1,Nspec
                  rDmixe = 0.5d0*(coefc(i-1,j,DMIXoffset+n)+coefc(i,j,DMIXoffset+n))
                  if (Full0_Mix1.eq.1) then 
                     fluxYH(i,j,n) = - rDmixe * dY(n)
                  endif
                  if (fillAlpha.eq.1) then
                     if (i-1.ge.lo(1)) then
                        alphaYH(i-1,j,n) = alphaYH(i-1,j,n)  +  rDmixe * dx2Inv
                     endif
                     if (i.le.hi(1)) then
                        alphaYH(i  ,j,n) = alphaYH(i  ,j,n)  +  rDmixe * dx2Inv
                     endif
                  endif
               enddo

               if (for_T0_H1.eq.1) then
                  do n = 1,Nspec
                     Hie = 0.5d0*(Hic(i-1,j,n)+Hic(i,j,n))
                     fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1)  +  Hie * fluxYH(i,j,n)
                  enddo
               else
                  FcpDTe(i,j) = 0.d0               
                  do n = 1,Nspec
                     CPie = 0.5d0*(CPic(i-1,j,n)+CPic(i,j,n))
                     FcpDTe(i,j) = FcpDTe(i,j)  +  CPie * fluxYH(i,j,n) * dT
                  enddo
               endif

               PTCe = 0.5d0*(coefc(i-1,j,PTCoffset+1)+coefc(i,j,PTCoffset+1))
               fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1)  -  PTCe * dT
               if (fillAlpha.eq.1) then
                  if (i-1.ge.lo(1)) then
                     alphaYH(i-1,j,Nspec+1) = alphaYH(i-1,j,Nspec+1)  +  PTCe * dx2Inv
                  endif
                  if (i.le.hi(1)) then
                     alphaYH(i  ,j,Nspec+1) = alphaYH(i  ,j,Nspec+1)  +  PTCe * dx2Inv
                  endif
               endif
            enddo
         enddo
      else
         do i = lo(1), hi(1)
            do j = lo(2), hi(2)+1
               do n = 1,Nspec
                  Ye(n) = ( YTc(i,j,n) + YTc(i,j-1,n)  ) * 0.5d0
                  dX(n) = (  Xc(i,j,n) -  Xc(i,j-1,n)  ) * dxInv
                  dY(n) = ( YTc(i,j,n) - YTc(i,j-1,n)  ) * dxInv
               enddo
               dT = ( YTc(i,j,Nspec+1) - YTc(i,j-1,Nspec+1) ) * dxInv
               Te = 0.5d0 * ( YTc(i,j,Nspec+1) + YTc(i,j-1,Nspec+1) )
               call CKMMWY(Ye, IWRK(ckbi), RWRK(ckbr), MMWe)

               fluxYH(i,j,Nspec+1) = 0.d0
               if (Full0_Mix1 .eq.0) then
                  do n = 1,Nspec
                     if (Ye(n).le.1.d-30) then
                        rhoTDe = 0.5d0*(coefc(i,j-1,TDoffset+n) + coefc(i,j  ,TDoffset+n))
                     else
                        rhoTDe = 0.5d0*(coefc(i,j-1,TDoffset+n)*YTc(i,j-1,n)
     &                       +          coefc(i,j  ,TDoffset+n)*YTc(i,j  ,n))/Ye(n)
                     endif
                     fluxYH(i,j,n) = - Ye(n) * rhoTDe * dT / Te
c     Note: Forming D_(n,m) . Grad(X_m), and assuming D was loaded with n varying fastest
                     do m = 1,Nspec
                        idx = Doffset + (m-1)*Nspec + n
                        rhoY_Dij = 0.5d0*(coefc(i,j-1,idx)+coefc(i,j,idx))
                        fluxYH(i,j,n) = fluxYH(i,j,n)  -  rhoY_Dij * dX(m)
                     enddo
                     fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1) - (RU*Te/MMWe)*rhoTDe*dX(n)
                  enddo
               endif
               do n = 1,Nspec
                  rDmixe = 0.5d0*(coefc(i,j-1,DMIXoffset+n)+coefc(i,j,DMIXoffset+n))
                  if (Full0_Mix1.eq.1) then 
                     fluxYH(i,j,n) = - rDmixe * dY(n)
                  endif
                  if (fillAlpha.eq.1) then
                     if (j-1.ge.lo(2)) then
                        alphaYH(i,j-1,n) = alphaYH(i,j-1,n)  +  rDmixe * dx2Inv
                     endif
                     if (j.le.hi(2)) then
                        alphaYH(i  ,j,n) = alphaYH(i  ,j,n)  +  rDmixe * dx2Inv
                     endif
                  endif
               enddo
               if (for_T0_H1.eq.1) then
                  do n = 1,Nspec
                     Hie = 0.5d0*(Hic(i,j-1,n)+Hic(i,j,n))
                     fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1)  +  Hie * fluxYH(i,j,n)
                  enddo
               else
                  FcpDTe(i,j) = 0.d0               
                  do n = 1,Nspec
                     CPie = 0.5d0*(CPic(i,j-1,n)+CPic(i,j,n))
                     FcpDTe(i,j) = FcpDTe(i,j)  +  CPie * fluxYH(i,j,n) * dT
                  enddo
               endif

               PTCe = 0.5d0*(coefc(i,j-1,PTCoffset+1)+coefc(i,j,PTCoffset+1))
               fluxYH(i,j,Nspec+1) = fluxYH(i,j,Nspec+1)  -  PTCe * dT
               if (fillAlpha.eq.1) then
                  if (j-1.ge.lo(2)) then
                     alphaYH(i,j-1,Nspec+1) = alphaYH(i,j-1,Nspec+1)  +  PTCe * dx2Inv
                  endif
                  if (j  .le.hi(2)) then
                     alphaYH(i  ,j,Nspec+1) = alphaYH(i  ,j,Nspec+1)  +  PTCe * dx2Inv
                  endif
               endif
            enddo
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
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               do n = 1,nc
                  efab(i,j,n) = half*(cfab(i-1,j,n) + cfab(i,j,n))
               enddo
            enddo
         enddo
      else
         do i = lo(1), hi(1)
            do j = lo(2), hi(2)+1
               do n = 1,nc
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
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n= 1,nc
                     div(i,j,n) = div(i,j,n)
     &                    +a*(V(i+1,j,n) - V(i,j,n))
                  enddo
               enddo               
            enddo
         else 
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n= 1,nc
                     div(i,j,n) = div(i,j,n)
     &                    + a*(V(i,j+1,n) - V(i,j,n))
                  enddo
               enddo               
            enddo
         endif
      else
         if (dir .EQ. 0) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n= 1,nc
                     div(i,j,n) = div(i,j,n)
     &                    + a*half*(V(i+1,j,n) + V(i,j,n))
                  enddo
               enddo               
            enddo
         else 
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n= 1,nc
                     div(i,j,n) = div(i,j,n)
     &                    + a*half*(V(i,j+1,n) + V(i,j,n))
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
      
      integer icrse,jcrse,n,ii,jj,ifine,jfine, clo1, chi1, clo2, chi2, offset, dir
      
      if ((face.EQ.0).OR.(face.EQ.2)) then
         dir = 0
      else
         dir = 1
      end if

      if (face .LT. SDIM) then
         offset = 1
      else
         offset = 0
      end if
      
c     NOTE: need to generalize to rz
      if (dir .eq. 0) then
         clo1 = (lo(1)+offset)/ratio(1)-offset
         clo2 = lo(2)/ratio(2)
         chi1 = (hi(1)+offset)/ratio(1)-offset
         chi2 = (hi(2)+1)/ratio(2) - 1
         do jcrse=clo2,chi2
            do icrse=clo1,chi1
               do n=1,nc
                  crse(icrse,jcrse,n) = 0.d0
                  ifine = ratio(1)*(icrse + offset) - offset
                  do jj = 1,ratio(2)
                     jfine = ratio(2)*jcrse + jj - 1
                     crse(icrse,jcrse,n) = crse(icrse,jcrse,n) + fine(ifine,jfine,n)
                  enddo
                  crse(icrse,jcrse,n) = crse(icrse,jcrse,n)/ratio(2)
               enddo
            enddo
         enddo
      else
         clo1 = lo(1)/ratio(1)
         clo2 = (lo(2)+offset)/ratio(2)-offset
         chi1 = (hi(1)+1)/ratio(1) - 1
         chi2 = (hi(2)+offset)/ratio(2)-offset
         do jcrse=clo2,chi2
            do icrse=clo1,chi1
               do n=1,nc
                  crse(icrse,jcrse,n) = 0.d0
                  jfine = ratio(2)*(jcrse + offset) - offset
                  do ii = 1,ratio(1)
                     ifine = ratio(1)*icrse + ii - 1
                     crse(icrse,jcrse,n) = crse(icrse,jcrse,n) + fine(ifine,jfine,n) 
                  enddo
                  crse(icrse,jcrse,n) = crse(icrse,jcrse,n)/ratio(1)
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
      
