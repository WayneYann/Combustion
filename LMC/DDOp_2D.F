#include "DDOp_F.H"
#include "ChemDriver_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"

#define SDIM 2

      subroutine FORT_DDFLUX(lo, hi, delta, dir,
     &                       fluxYH, DIMS(fluxYH),
     &                       YTc, DIMS(YTc), cpe, DIMS(cpe),
     &                       for_T0_H1, fillAlpha, alphaYH, DIMS(alphaYH) )
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(fluxYH)
      integer DIMDEC(YTc)
      integer DIMDEC(cpe)
      integer for_T0_H1, fillAlpha
      integer DIMDEC(alphaYH)
      integer dir
      REAL_T  delta
      REAL_T  fluxYH(DIMV(fluxYH),*)
      REAL_T     YTc(DIMV(YTc)   ,*)
      REAL_T     cpe(DIMV(cpe)   ,*)
      REAL_T  alphaYH(DIMV(alphaYH),*)

      integer i, j, n
      REAL_T XM(maxspec), XP(maxspec), YM(maxspec), YP(maxspec), TM, TP, CP(maxspec)
      REAL_T f(maxspec+1), alpha(maxspec+1)

      if (dir.eq.0) then

         do j = lo(2), hi(2)
            do n = 1,Nspec
               YP(n)  = YTc(lo(1)-1,j,n)
            enddo
            CALL CKYTX(YP,IWRK(ckbi),RWRK(ckbr),XP)
            TP = YTc(lo(1)-1,j,Nspec+1)            
            do i = lo(1), hi(1)+1
               do n = 1,Nspec
                  YM(n) = YP(n)
                  YP(n) = YTc(i,j,n)
                  XM(n) = XP(n)
                  CP(n) = cpe(i,j,n)
               enddo
               CALL CKYTX(YP,IWRK(ckbi),RWRK(ckbr),XP)
               TM = TP
               TP = YTc(i,j,Nspec+1)

               call DDFLUX_helper(XM,XP,YM,YP,TM,TP,CP,delta,for_T0_H1,fillAlpha,f,alpha)

               do n = 1,Nspec+1
                  fluxYH(i,j,n) = f(n)
               enddo

               if (for_T0_H1.eq.0) then
                  do n=1,Nspec
                     cpe(i,j,n) = CP(n)
                  enddo
               endif

               if (fillAlpha.eq.1) then
                  if (i.le.hi(1)) then
                     do n=1,Nspec+1
                        alphaYH(i,j,n) = alphaYH(i,j,n) + alpha(n)
                     enddo
                  endif
                  if (i.gt.lo(1)) then
                     do n=1,Nspec+1
                        alphaYH(i-1,j,n) = alphaYH(i-1,j,n) + alpha(n)
                     enddo
                  endif
               endif

            enddo
         enddo
      else if (dir.eq.1) then
         do i = lo(1), hi(1)
            do n = 1,Nspec
               YP(n)  = YTc(i,lo(2)-1,n)
            enddo
            CALL CKYTX(YP,IWRK(ckbi),RWRK(ckbr),XP)
            TP = YTc(i,lo(2)-1,Nspec+1)            
            do j = lo(2), hi(2)+1
               do n = 1,Nspec
                  YM(n) = YP(n)
                  YP(n) = YTc(i,j,n)
                  XM(n) = XP(n)
                  CP(n) = cpe(i,j,n)
               enddo
               CALL CKYTX(YP,IWRK(ckbi),RWRK(ckbr),XP)
               TM = TP
               TP = YTc(i,j,Nspec+1)
               
               call DDFLUX_helper(XM,XP,YM,YP,TM,TP,CP,delta,for_T0_H1,fillAlpha,f,alpha)

               do n = 1,Nspec+1
                  fluxYH(i,j,n) = f(n)
               enddo

               if (for_T0_H1.eq.0) then
                  do n=1,Nspec
                     cpe(i,j,n) = CP(n)
                  enddo
               endif

               if (fillAlpha.eq.1) then
                  if (j.le.hi(2)) then
                     do n=1,Nspec+1
                        alphaYH(i,j,n) = alphaYH(i,j,n) + alpha(n)
                     enddo
                  endif
                  if (j.gt.lo(2)) then
                     do n=1,Nspec+1
                        alphaYH(i,j-1,n) = alphaYH(i,j-1,n) + alpha(n)
                     enddo
                  endif
               endif

            enddo
         enddo

c         print *,'for_T0_H1:',for_T0_H1
c         open(unit=21,file="junk")
c         do j=115,125
c            write(21,'(i4,20g13.5)') j,(j+0.5d0)*delta*100,(fluxYH(15,j,n),n=1,Nspec+1)
c         enddo
c         stop

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
      
