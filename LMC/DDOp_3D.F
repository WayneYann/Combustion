#include "DDOp_F.H"
#include "ChemDriver_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"

#define SDIM 3

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

      integer i, j, k, n
      REAL_T XM(maxspec), XP(maxspec), YM(maxspec), YP(maxspec), TM, TP, CP(maxspec)
      REAL_T f(maxspec+1), alpha(maxspec+1)

      if (dir.eq.0) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do n = 1,Nspec
                  YP(n)  = YTc(lo(1)-1,j,k,n)
               enddo
               CALL CKYTX(YP,IWRK(ckbi),RWRK(ckbr),XP)
               TP = YTc(lo(1)-1,j,k,Nspec+1)            
               do i = lo(1), hi(1)+1
                  do n = 1,Nspec
                     YM(n) = YP(n)
                     YP(n) = YTc(i,j,k,n)
                     XM(n) = XP(n)
                     CP(n) = cpe(i,j,k,n)
                  enddo
                  CALL CKYTX(YP,IWRK(ckbi),RWRK(ckbr),XP)
                  TM = TP
                  TP = YTc(i,j,k,Nspec+1)
                  
                  call DDFLUX_helper(XM,XP,YM,YP,TM,TP,CP,delta,for_T0_H1,fillAlpha,f,alpha)
                  
                  do n = 1,Nspec+1
                     fluxYH(i,j,k,n) = f(n)
                  enddo
                  
                  if (for_T0_H1.eq.0) then
                     do n=1,Nspec
                        cpe(i,j,k,n) = CP(n)
                     enddo
                  endif
                  
                  if (fillAlpha.eq.1) then
                     if (i.le.hi(1)) then
                        do n=1,Nspec+1
                           alphaYH(i,j,k,n) = alphaYH(i,j,k,n) + alpha(n)
                        enddo
                     endif
                     if (i.gt.lo(1)) then
                        do n=1,Nspec+1
                           alphaYH(i-1,j,k,n) = alphaYH(i-1,j,k,n) + alpha(n)
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)
               do n = 1,Nspec
                  YP(n)  = YTc(i,lo(2)-1,k,n)
               enddo
               CALL CKYTX(YP,IWRK(ckbi),RWRK(ckbr),XP)
               TP = YTc(i,lo(2)-1,k,Nspec+1)            
               do j = lo(2), hi(2)+1
                  do n = 1,Nspec
                     YM(n) = YP(n)
                     YP(n) = YTc(i,j,k,n)
                     XM(n) = XP(n)
                     CP(n) = cpe(i,j,k,n)
                  enddo
                  CALL CKYTX(YP,IWRK(ckbi),RWRK(ckbr),XP)
                  TM = TP
                  TP = YTc(i,j,k,Nspec+1)
                  
                  call DDFLUX_helper(XM,XP,YM,YP,TM,TP,CP,delta,for_T0_H1,fillAlpha,f,alpha)
                  
                  do n = 1,Nspec+1
                     fluxYH(i,j,k,n) = f(n)
                  enddo
                  
                  if (for_T0_H1.eq.0) then
                     do n=1,Nspec
                        cpe(i,j,k,n) = CP(n)
                     enddo
                  endif
                  
                  if (fillAlpha.eq.1) then
                     if (j.le.hi(2)) then
                        do n=1,Nspec+1
                           alphaYH(i,j,k,n) = alphaYH(i,j,k,n) + alpha(n)
                        enddo
                     endif
                     if (j.gt.lo(2)) then
                        do n=1,Nspec+1
                           alphaYH(i,j-1,k,n) = alphaYH(i,j-1,k,n) + alpha(n)
                        enddo
                     endif
                  endif
            enddo
         enddo
      else if (dir.eq.2) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do n = 1,Nspec
                  YP(n)  = YTc(i,j,lo(3)-1,n)
               enddo
               CALL CKYTX(YP,IWRK(ckbi),RWRK(ckbr),XP)
               TP = YTc(i,k,lo(3)-1,Nspec+1)            
               do j = lo(3), hi(3)+1
                  do n = 1,Nspec
                     YM(n) = YP(n)
                     YP(n) = YTc(i,j,k,n)
                     XM(n) = XP(n)
                     CP(n) = cpe(i,j,k,n)
                  enddo
                  CALL CKYTX(YP,IWRK(ckbi),RWRK(ckbr),XP)
                  TM = TP
                  TP = YTc(i,j,k,Nspec+1)
                  
                  call DDFLUX_helper(XM,XP,YM,YP,TM,TP,CP,delta,for_T0_H1,fillAlpha,f,alpha)
                  
                  do n = 1,Nspec+1
                     fluxYH(i,j,k,n) = f(n)
                  enddo
                  
                  if (for_T0_H1.eq.0) then
                     do n=1,Nspec
                        cpe(i,j,k,n) = CP(n)
                     enddo
                  endif
                  
                  if (fillAlpha.eq.1) then
                     if (k.le.hi(3)) then
                        do n=1,Nspec+1
                           alphaYH(i,j,k,n) = alphaYH(i,j,k,n) + alpha(n)
                        enddo
                     endif
                     if (j.gt.lo(3)) then
                        do n=1,Nspec+1
                           alphaYH(i,j-1,k,n) = alphaYH(i,j-1,k,n) + alpha(n)
                        enddo
                     endif
                  endif
               enddo
            enddo
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
      REAL_T div(DIMV(div) ,1)
      REAL_T   V(DIMV(V)   ,1)
      REAL_T a

      integer i,j,k,n

      if (d0a1.eq.0) then
         if (dir .EQ. 0) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     do n= 1,nc
                        div(i,j,k,n) = div(i,j,k,n)
     &                       + a*(V(i+1,j,k,n) - V(i,j,k,n))
                     enddo
                  enddo
               enddo               
            enddo
         elseif (dir .EQ. 1) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     do n= 1,nc
                        div(i,j,k,n) = div(i,j,k,n)
     &                       + a*(V(i,j+1,k,n) - V(i,j,k,n))
                     enddo
                  enddo
               enddo               
            enddo
         else 
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     do n= 1,nc
                        div(i,j,k,n) = div(i,j,k,n)
     &                       + a*(V(i,j,k+1,n) - V(i,j,k,n))
                     enddo
                  enddo
               enddo               
            enddo
         end if
      else
         if (dir .EQ. 0) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     do n= 1,nc
                        div(i,j,k,n) = div(i,j,k,n)
     &                       + a*half*(V(i+1,j,k,n) + V(i,j,k,n))
                     enddo
                  enddo
               enddo               
            enddo
         elseif (dir .EQ. 1) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     do n= 1,nc
                        div(i,j,k,n) = div(i,j,k,n)
     &                       + a*half*(V(i,j+1,k,n) + V(i,j,k,n))
                     enddo
                  enddo
               enddo               
            enddo
         else 
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     do n= 1,nc
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
      
