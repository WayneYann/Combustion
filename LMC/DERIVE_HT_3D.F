#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROB_F.H"
#include "ArrayLim.H"
#include "ChemKinDriver_F.H"

#   if   BL_SPACEDIM==1
#       define  ARLIM(x)  x(1)
#   elif BL_SPACEDIM==2
#       define  ARLIM(x)  x(1),x(2)
#   elif BL_SPACEDIM==3
#       define  ARLIM(x)  x(1),x(2),x(3)
#   endif

#define SDIM 3

      subroutine FORT_DERRHOMINUSSUMRHOY (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc,
     $                         level,grid_no)
      implicit none
c
c ::: This routine will computes rho - sum (rho*Y)
c

#include "htdata.H"

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k,n
      integer    nxlo,nxhi,nylo,nyhi,nzlo,nzhi
      
c     ::::: lets punt if not in domain interior
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
	 call bl_abort("FORT_DERRHOMINUSSUMRHOY: outside domain")
      endif

      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,1)
            enddo
         enddo
      enddo
      
      do n=2,ncomp
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  e(i,j,k,1) = e(i,j,k,1)-dat(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
c     
c     end of routine
c
      return
      end

c=========================================================

      subroutine FORT_DERSUMYDOT (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc,
     $                         level,grid_no)
      implicit none
c
c ::: This routine will computes sum (Ydot)
c

#include "htdata.H"

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k,n
      integer    nxlo, nxhi, nylo, nyhi, nzlo,nzhi
      
c     ::::: lets punt if not in domain interior
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
	 call bl_abort("FORT_DERSUMYDOT: outside domain")
      endif

      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = 0.0
            enddo
         enddo
      enddo
      
      do n=1,ncomp
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  e(i,j,k,1) = e(i,j,k,1)+dat(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
c     
c     end of routine
c
      return
      end

c=========================================================

      subroutine FORT_DERRHORT (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc,
     $                         level,grid_no)
      implicit none
c     
c ::: This routine will derive rho*R*T
c

#include "htdata.H"
#include "cdwrk.H"
      
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i, j, k, n, rho, T, fS
      integer    nxlo,nxhi,nylo,nyhi,nzlo,nzhi
      REAL_T     Yt(maxspec)
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/
      
c     ::::: lets punt if not in domain interior
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
	 call bl_abort("FORT_DERRHORT: outside domain")
      endif
      
c     Set pointers into state (these must agree with setup for this derived quant)
      rho = 1
      T = 2
      fS = 3

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = dat(i,j,k,fS+n-1) / dat(i,j,k,rho)
               end do
               call FORT_PfromRTY(lo_chem, hi_chem, 
     &              e(i,j,k,1),     ARLIM(lo_chem),ARLIM(hi_chem), 
     &              dat(i,j,k,rho), ARLIM(lo_chem),ARLIM(hi_chem),
     &              dat(i,j,k,T),   ARLIM(lo_chem),ARLIM(hi_chem), 
     &              Yt,             ARLIM(lo_chem),ARLIM(hi_chem))
            end do
         end do
      end do
      end

      subroutine FORT_DERMOLEFRAC (x,DIMS(x),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc,
     $                         level,grid_no)
      implicit none
#include "cdwrk.H"
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(x)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     x(DIMV(x),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer i,j,k,n
      REAL_T Yt(maxspec),Xt(maxspec)
      integer fS,rho
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      rho = 1 
      fS = 2

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = dat(i,j,k,fS+n-1)/dat(i,j,k,rho) 
               enddo
               call FORT_MASSTOMOLE(lo_chem, hi_chem,
     &              Yt, ARLIM(lo_chem),ARLIM(hi_chem),
     &              Xt, ARLIM(lo_chem),ARLIM(hi_chem))
               do n = 1,Nspec
                  x(i,j,k,n) = Xt(n)
               enddo
            enddo
         enddo
      enddo

      end


      subroutine FORT_DERCONCENTRATION(C,DIMS(C),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc,
     $                         level,grid_no)
      implicit none
#include "cdwrk.H"
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(C)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     C(DIMV(C),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer i,j,k,n
      REAL_T Yt(maxspec),Ct(maxspec)
      integer fS,rho,T
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      rho = 1 
      T   = 2
      fS  = 3

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = dat(i,j,k,fS+n-1)/dat(i,j,k,rho) 
               enddo
            call FORT_MASSR_TO_CONC(lo_chem,hi_chem,
     &           Yt,             ARLIM(lo_chem),ARLIM(hi_chem),
     &           dat(i,j,k,T),   ARLIM(lo_chem),ARLIM(hi_chem),
     &           dat(i,j,k,rho), ARLIM(lo_chem),ARLIM(hi_chem),
     &           Ct,             ARLIM(lo_chem),ARLIM(hi_chem))
               do n = 1,Nspec
                  C(i,j,k,n) = Ct(n)
               enddo
            enddo
         enddo
      enddo

      end

