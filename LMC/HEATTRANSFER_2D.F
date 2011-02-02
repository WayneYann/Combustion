c
c $Id: HEATTRANSFER_2D.F,v 1.34 2011-02-02 01:36:44 marc Exp $
c
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "HEATTRANSFER_F.H"
#include "ChemDriver_F.H"
#include "ArrayLim.H"

#   if   BL_SPACEDIM==1
#       define  ARLIM(x)  x(1)
#   elif BL_SPACEDIM==2
#       define  ARLIM(x)  x(1),x(2)
#   elif BL_SPACEDIM==3
#       define  ARLIM(x)  x(1),x(2),x(3)
#   endif

#define SDIM 2

      subroutine FORT_COMPUTE_RHODGRADHDOTGRADY(dx,
     &         lo, hi, DIMS(species), species,
     &         DIMS(h), h, DIMS(betax), betax,
     &         DIMS(betay), betay, DIMS(rdghdgy), rdghdgy)

      implicit none

c ... inputs

      integer lo(SDIM), hi(SDIM)
      REAL_T  dx(SDIM)
      integer DIMDEC(species)
      integer DIMDEC(h)
      REAL_T  species(DIMV(species))
      REAL_T  h(DIMV(h))
      integer DIMDEC(betax)
      integer DIMDEC(betay)
      REAL_T betax(DIMV(betax))
      REAL_T betay(DIMV(betay))
      integer DIMDEC(rdghdgy)

c ... outputs

      REAL_T rdghdgy(DIMV(rdghdgy))

c ... local

      integer i,j
      REAL_T  betaxcen, betaycen, dxsqr, dysqr
      REAL_T  betadotleft, betadotright
      REAL_T  betadottop,  betadotbot

c -------------------- begin here ----------------------------

      dxsqr = dx(1)**2
      dysqr = dx(2)**2
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
#if 0
          betaxcen = half*(betax(i,j)+betax(i+1,j))
          betaycen = half*(betay(i,j)+betay(i,j+1))
          rdghdgy(i,j) = betaxcen*(h(i+1,j)-h(i-1,j))*
     &                            (species(i+1,j)-species(i-1,j))/dxsqr +
     &                   betaycen*(h(i,j+1)-h(i,j-1))*
     &                             (species(i,j+1)-species(i,j-1))/dysqr
#else
          betadotleft  = betax(i,j)*(h(i,j)-h(i-1,j))*(species(i,j)-species(i-1,j))
          betadotright = betax(i+1,j)*(h(i+1,j)-h(i,j))*(species(i+1,j)-species(i,j))
          betadotbot   = betay(i,j)*(h(i,j)-h(i,j-1))*(species(i,j)-species(i,j-1))
          betadottop   = betay(i,j+1)*(h(i,j+1)-h(i,j))*(species(i,j+1)-species(i,j))
          rdghdgy(i,j) =  half*((betadotleft + betadotright)/dxsqr +
     &                          (betadottop  + betadotbot)/dysqr)
#endif
        enddo
      enddo
      end

      subroutine FORT_VELVISC(lo,hi,
     &                        DIMS(T), T,
     &                        DIMS(Y), Y,
     &                        DIMS(mu), mu)
      implicit none
#include "visc.H"
#include "cdwrk.H"      
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T  T(DIMV(T))      
      REAL_T  Y(DIMV(Y),*)      
      integer DIMDEC(mu)
      REAL_T  mu(DIMV(mu))

      integer i, j

      if (.not.use_constant_mu) then
         if (.not. LeEQ1) then
            call FORT_MIX_SHEAR_VISC(lo, hi, mu, DIMS(mu),
     &                               T, DIMS(T), Y, DIMS(Y))
         else 
            do j=lo(2), hi(2)
               do i=lo(1), hi(1)
                  mu(i,j) = 1.85e-5*(T(i,j)/298.0)**.7
               end do
            end do
         end if
      else
         do j=lo(2), hi(2)
            do i=lo(1), hi(1)
               mu(i,j) = constant_mu_val
            end do
         end do
      end if
      end

      subroutine FORT_SPECTEMPVISC(lo,hi,
     &     DIMS(T), T,
     &     DIMS(Y), Y,
     &     DIMS(rhoD), rhoD,
     &     ncompd, P1ATM_MKS, do_temp, do_VelVisc)
      implicit none
#include "cdwrk.H"
#include "visc.H"
#include "htdata.H"
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer DIMDEC(rhoD)
      integer ncompd, do_temp, do_VelVisc
      REAL_T  T(DIMV(T))
      REAL_T  Y(DIMV(Y),*)
      REAL_T  rhoD(DIMV(rhoD),ncompd)
      
      integer i, j
      integer n, nspecies
      REAL_T P1ATM_MKS, Patm
      REAL_T Yt(maxspec)
      REAL_T cpmix, Tfac, Yfac
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1/
      data hi_chem /1,1/

      nspecies = ncompd - 2

c     ncompd-1 is Temp stuff
c     ncompd   is Velocity stuff 

c     Warning, FORT_VELVISC is called separately from this routine, so if there's
c     any hacking to be done on viscosity, be sure to catch it there as well.
      Tfac = thickFacTR / Pr
      Yfac = thickFacTR / Sc

      if (.not.use_constant_rhoD) then
         if (.not. LeEQ1) then
            Patm = pamb / P1ATM_MKS
            call FORT_MIXAVG_RHODIFF_TEMP(lo, hi, rhoD, DIMS(rhoD),
     &           T, DIMS(T), Y, DIMS(Y), Patm, do_temp, do_VelVisc)
            if (thickFacTR.ne.1.d0) then
               do n=1,ncompd-1
                  do j=lo(2), hi(2)
                     do i=lo(1), hi(1)
                        rhoD(i,j,n) = rhoD(i,j,n)*thickFacTR
                     end do
                  end do
               end do
            endif
         else
            call FORT_VELVISC(lo,hi,DIMS(T),T,DIMS(Y),Y,
     &           DIMS(rhoD),rhoD)
            do j=lo(2), hi(2)
               do i=lo(1), hi(1)
                  rhoD(i,j,ncompd-1) = rhoD(i,j,1)
                  rhoD(i,j,ncompd  ) = rhoD(i,j,1)
                  rhoD(i,j,1)        = rhoD(i,j,1) * Yfac
               end do
            end do
            do n=2,nspecies
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     rhoD(i,j,n) = rhoD(i,j,1)
                  end do
               end do
            end do

            if (do_temp .ne. 0) then
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     do n=1,nspecies
                        Yt(n) = Y(i,j,n)
                     end do
                     CALL FORT_CPMIXfromTY(lo_chem, hi_chem, 
     &                    cpmix,  ARLIM(lo_chem), ARLIM(hi_chem),
     &                    T(i,j), ARLIM(lo_chem), ARLIM(hi_chem),
     &                    Yt,     ARLIM(lo_chem), ARLIM(hi_chem))
                     rhoD(i,j,ncompd-1) = rhoD(i,j,ncompd-1)*cpmix*Tfac
                  end do
               end do
            end if
         end if
      else
         do n=1,nspecies
            do j=lo(2), hi(2)
               do i=lo(1), hi(1)
                  rhoD(i,j,n) = constant_rhoD_val
               end do
            end do
         end do
      end if
      end

      subroutine FORT_EST_DIVU_DT(flag, dtfactor, delta, divu, DIMS(divu), 
     &                           dsdt, rho, DIMS(rho), 
     &                           u, DIMS(u), 
     &                           volume, DIMS(volume), 
     &                           areax,  DIMS(areax),
     &                           areay,  DIMS(areay),
     &                           lo, hi, dt, rhomin)

      implicit none

      integer flag
      integer lo(SDIM), hi(SDIM)
      REAL_T  delta(SDIM)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
      integer DIMDEC(u)
      REAL_T  rho(DIMV(rho))      
      REAL_T  u(DIMV(u),BL_SPACEDIM)      
      REAL_T  rhomin, dtfactor
      REAL_T  divu(DIMV(divu))
      REAL_T  dsdt(DIMV(divu))

      integer DIMDEC(volume)
      integer DIMDEC(areax)
      integer DIMDEC(areay)
      REAL_T  volume(DIMV(volume))
      REAL_T  areax(DIMV(areax))
      REAL_T  areay(DIMV(areay))

      REAL_T dt

      integer i,j
      REAL_T  dtcell, dtcell2, denom, rhominij, rhoij
      REAL_T  fluxtop, fluxbot, fluxleft, fluxright
      REAL_T  a,b,c

      dt = 1.0e20
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          dtcell = dt
          if (flag.eq.1) then
            if(divu(i,j).gt.zero) then
              if(rho(i,j).gt.rhomin) then
                dtcell = dtfactor*(one-rhomin/rho(i,j))/divu(i,j)
              else
                dtcell = dtfactor*.5/divu(i,j)
              endif
              if (dsdt(i,j).gt.1.0e-20) then
                if (abs(rho(i,j)).gt.rhomin) then
                  rhominij = rhomin
                else
                  rhominij = .9*abs(rho(i,j)) 
                endif
                rhoij = abs(rho(i,j))
#if 0
                dtcell2 = (-rhoij*divu(i,j) + 
     &                     sqrt((rhoij*divu(i,j))**2+
     &                          two*dsdt(i,j)*rhoij*abs(rhoij-rhominij)))/
     &                    (rhoij*dsdt(i,j))
#else
c ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
c           We use the latter because it is more robust
                a = rhoij*dsdt(i,j)*half
                b = rhoij*divu(i,j)
                c = rhominij - rhoij
                dtcell2 = two*c/(-b-sqrt(b**2-four*a*c))
#endif
                dtcell2 = dtfactor*dtcell2
                dtcell = min(dtcell,dtcell2)
              endif
            endif
            if(dtcell.le.0.0)then
              write(6,*)'aha'
            endif
          else if (flag.eq.2) then
            denom = rho(i,j)*divu(i,j)+
     &              u(i,j,1)*(rho(i+1,j)-rho(i-1,j))/delta(1) +
     &              u(i,j,2)*(rho(i,j+1)-rho(i,j-1))/delta(2)
            if(denom.gt.zero)then
              if(rho(i,j).gt.rhomin) then
                dtcell = dtfactor*(rho(i,j)-rhomin)/denom
              else
                dtcell = dtfactor*abs(rho(i,j))/denom
              endif
            endif
          else if (flag.eq.3) then
            fluxtop   = fourth*(rho(i,j)+rho(i,j+1))*(u(i,j,2)+u(i,j+1,2))
            fluxbot   = fourth*(rho(i,j)+rho(i,j-1))*(u(i,j,2)+u(i,j-1,2))
            fluxleft  = fourth*(rho(i,j)+rho(i-1,j))*(u(i,j,1)+u(i-1,j,1))
            fluxright = fourth*(rho(i,j)+rho(i+1,j))*(u(i,j,1)+u(i+1,j,1))
            denom = ((areax(i+1,j)*fluxright-areax(i,j)*fluxleft)+
     &               (areay(i,j+1)*fluxtop-areay(i,j)*fluxbot))/volume(i,j)

            if(denom.gt.zero)then
              if(rho(i,j).gt.rhomin) then
                dtcell = dtfactor*(rho(i,j)-rhomin)/denom
              else
                dtcell = dtfactor*abs(rho(i,j))/denom
              endif
            endif
          endif
#if 0
          write(6,*)'i,j,dtcell=',i,j,dtcell
#endif
          dt = min(dtcell,dt)
        enddo
      enddo
      end

      subroutine FORT_CHECK_DIVU_DT(flag, dtfactor, delta, divu, DIMS(divu), 
     &                           dsdt, rho, DIMS(rho), 
     &                           u, DIMS(u), 
     &                           volume, DIMS(volume), 
     &                           areax,  DIMS(areax),
     &                           areay,  DIMS(areay),
     &                           lo, hi, dt, rhomin)

      implicit none

      integer flag
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
      integer DIMDEC(u)
      REAL_T  delta(SDIM)
      REAL_T  rho(DIMV(rho))      
      REAL_T  u(DIMV(u),BL_SPACEDIM)      
      REAL_T  rhomin
      REAL_T  divu(DIMV(divu))
      REAL_T  dsdt(DIMV(divu))
      REAL_T  dt, dtfactor

      integer DIMDEC(volume)
      integer DIMDEC(areax)
      integer DIMDEC(areay)
      REAL_T  volume(DIMV(volume))
      REAL_T  areax(DIMV(areax))
      REAL_T  areay(DIMV(areay))

      integer i,j
      REAL_T  dtcell, denom
      REAL_T  fluxtop, fluxbot, fluxleft, fluxright
      REAL_T a,b,c,dtcell2,rhominij,rhoij

      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          dtcell = bigreal
          if (flag.eq.1) then
            if(divu(i,j).gt.zero) then
              if(rho(i,j).gt.rhomin) then
                dtcell = (one-rhomin/rho(i,j))/divu(i,j)
              else
                dtcell = one/divu(i,j)
              endif
              if (dsdt(i,j).gt.1.0e-20) then
                if (abs(rho(i,j)).gt.rhomin) then
                  rhominij = rhomin
                else
                  rhominij = .9*abs(rho(i,j)) 
                endif
                rhoij = abs(rho(i,j))
#if 0
                dtcell2 = (-rhoij*divu(i,j) + 
     &                     sqrt((rhoij*divu(i,j))**2+
     &                          two*dsdt(i,j)*rhoij*abs(rhoij-rhominij)))/
     &                    (rhoij*dsdt(i,j))
#else
c ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
c           We use the latter because it is more robust
                a = rhoij*dsdt(i,j)*half
                b = rhoij*divu(i,j)
                c = rhominij - rhoij
                dtcell2 = two*c/(-b-sqrt(b**2-four*a*c))
#endif
                dtcell = min(dtcell,dtcell2)
              endif

            endif
          else if (flag.eq.2) then
            denom = rho(i,j)*divu(i,j)+
     &              u(i,j,1)*(rho(i+1,j)-rho(i-1,j))/delta(1) +
     &              u(i,j,2)*(rho(i,j+1)-rho(i,j-1))/delta(2)
            if(denom.gt.zero)then
              if(rho(i,j).gt.rhomin) then
                dtcell = (rho(i,j)-rhomin)/denom
              else
                dtcell = abs(rho(i,j))/denom
              endif
            endif
          else if (flag.eq.3) then
            fluxtop = fourth*(rho(i,j)+rho(i,j+1))*(u(i,j,2)+u(i,j+1,2))
            fluxbot = fourth*(rho(i,j)+rho(i,j-1))*(u(i,j,2)+u(i,j-1,2))
            fluxleft = fourth*(rho(i,j)+rho(i-1,j))*(u(i,j,1)+u(i-1,j,1))
            fluxright = fourth*(rho(i,j)+rho(i+1,j))*(u(i,j,1)+u(i+1,j,1))
            denom = ((areax(i+1,j)*fluxright-areax(i,j)*fluxleft)+
     &               (areay(i,j+1)*fluxtop-areay(i,j)*fluxbot))/volume(i,j)
            if(denom.gt.zero)then
              if(rho(i,j).gt.rhomin) then
                dtcell = (rho(i,j)-rhomin)/denom
              else
                dtcell = abs(rho(i,j))/denom
              endif
            endif
          endif
          if (dt.gt.dtcell) then
            write(6,*)'ERROR: FORT_CHECK_DIVU_DT : i,j,dt>dtcell = ',
     &            i,j,dt,dtcell
          else if (dt.gt.dtcell*dtfactor) then
            write(6,*)'WARNING: ',
     &            'FORT_CHECK_DIVU_DT : i,j,dt>dtcell*dtfactor = ',
     &            i,j,dt,dtcell*dtfactor
          endif
        enddo
      enddo
      end

c===================================================================

c ... The routines YDOTFILL, DIVUFILL, DQRADFILL, and DSDTFILL
c     are here instead of in the problem dependent code because
c     we always fill the quantitities ydot, divu, dqrad, and dsdt
c     the same way

c===================================================================

      subroutine FORT_DQRADFILL (dqrad,DIMS(dqrad),domlo,domhi,delta,
     &                         xlo,time,bc )

      integer    DIMDEC(dqrad), bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     dqrad(DIMV(dqrad))

      integer    i, j
      integer    ilo, ihi, jlo, jhi

      lo(1) = ARG_L1(dqrad)
      hi(1) = ARG_H1(dqrad)
      lo(2) = ARG_L2(dqrad)
      hi(2) = ARG_H2(dqrad)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))

      call filcc (dqrad,DIMS(dqrad),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
           if(jlo.le.jhi)then
             do j = jlo, jhi
               do i = lo(1), domlo(1)-1
                 dqrad(i,j) = dqrad(domlo(1),j)
               enddo
             enddo
           endif
           if (lo(2).lt.domlo(2)) then
             do j = lo(2), domlo(2)-1
               do i = lo(1), domlo(1)-1
                 dqrad(i,j) = dqrad(domlo(1),domlo(2))
               enddo
             enddo
           endif
           if(hi(2).gt.domhi(2))then
             do j = domhi(2)+1, hi(2)
               do i = lo(1), domlo(1)-1
                 dqrad(i,j) = dqrad(domlo(1),domhi(2))
               enddo
             enddo
           endif

      endif            

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
           if(jlo.le.jhi)then
             do j = jlo,jhi
               do i = domhi(1)+1,hi(1)
                 dqrad(i,j) = dqrad(domhi(1),j)
               enddo
             enddo
           endif
           if (lo(2).lt.domlo(2)) then
             do j = lo(2), domlo(2)-1
               do i = domhi(1)+1,hi(1)
                 dqrad(i,j) = dqrad(domhi(1),domlo(2))
               enddo
             enddo
           endif
           if(hi(2).gt.domhi(2))then
             do j = domhi(2)+1, hi(2)
               do i = domhi(1)+1,hi(1)
                 dqrad(i,j) = dqrad(domhi(1),domhi(2))
               enddo
             enddo
           endif
      endif            

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
           if(ilo.le.ihi)then
             do j = lo(2), domlo(2)-1
               do i = ilo,ihi
                 dqrad(i,j) = dqrad(i,domlo(2))
               enddo
             enddo
           endif
           if (lo(1).lt.domlo(1)) then
             do j = lo(2), domlo(2)-1
               do i = lo(1), domlo(1)-1
                 dqrad(i,j) = dqrad(domlo(1),domlo(2))
               enddo
             enddo
           endif
           if(hi(1).gt.domhi(1))then
             do j = lo(2), domlo(2)-1
               do i = domhi(1)+1, hi(1)
                 dqrad(i,j) = dqrad(domhi(1),domlo(2))
               enddo
             enddo
           endif

      endif            

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
           if(ilo.le.ihi)then
             do j = domhi(2)+1, hi(2)
               do i = ilo,ihi
                 dqrad(i,j) = dqrad(i,domhi(2))
               enddo
             enddo
           endif
           if (lo(1).lt.domlo(1)) then
             do j = domhi(2)+1, hi(2)
               do i = lo(1), domlo(1)-1
                 dqrad(i,j) = dqrad(domlo(1),domhi(2))
               enddo
             enddo
           endif
           if(hi(1).gt.domhi(1))then
             do j = domhi(2)+1, hi(2)
               do i = domhi(1)+1, hi(1)
                 dqrad(i,j) = dqrad(domhi(1),domhi(2))
               enddo
             enddo
           endif
      endif            
      end

      subroutine FORT_DIVUFILL (divu,DIMS(divu),domlo,domhi,delta,
     &                        xlo,time,bc )

      integer    DIMDEC(divu)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     divu(DIMV(divu))

      integer    i, j
      integer    ilo, ihi, jlo, jhi

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(divu)
      hi(1) = ARG_H1(divu)
      lo(2) = ARG_L2(divu)
      hi(2) = ARG_H2(divu)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))

      call filcc (divu,DIMS(divu),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
           if(jlo.le.jhi)then
             do j = jlo, jhi
               do i = lo(1), domlo(1)-1
                 divu(i,j) = divu(domlo(1),j)
               enddo
             enddo
           endif
           if (lo(2).lt.domlo(2)) then
             do j = lo(2), domlo(2)-1
               do i = lo(1), domlo(1)-1
                 divu(i,j) = divu(domlo(1),domlo(2))
               enddo
             enddo
           endif
           if(hi(2).gt.domhi(2))then
             do j = domhi(2)+1, hi(2)
               do i = lo(1), domlo(1)-1
                 divu(i,j) = divu(domlo(1),domhi(2))
               enddo
             enddo
           endif

      endif            

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
           if(jlo.le.jhi)then
             do j = jlo,jhi
               do i = domhi(1)+1,hi(1)
                 divu(i,j) = divu(domhi(1),j)
               enddo
             enddo
           endif
           if (lo(2).lt.domlo(2)) then
             do j = lo(2), domlo(2)-1
               do i = domhi(1)+1,hi(1)
                 divu(i,j) = divu(domhi(1),domlo(2))
               enddo
             enddo
           endif
           if(hi(2).gt.domhi(2))then
             do j = domhi(2)+1, hi(2)
               do i = domhi(1)+1,hi(1)
                 divu(i,j) = divu(domhi(1),domhi(2))
               enddo
             enddo
           endif
      endif            

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
           if(ilo.le.ihi)then
             do j = lo(2), domlo(2)-1
               do i = ilo,ihi
                 divu(i,j) = divu(i,domlo(2))
               enddo
             enddo
           endif
           if (lo(1).lt.domlo(1)) then
             do j = lo(2), domlo(2)-1
               do i = lo(1), domlo(1)-1
                 divu(i,j) = divu(domlo(1),domlo(2))
               enddo
             enddo
           endif
           if(hi(1).gt.domhi(1))then
             do j = lo(2), domlo(2)-1
               do i = domhi(1)+1, hi(1)
                 divu(i,j) = divu(domhi(1),domlo(2))
               enddo
             enddo
           endif

      endif            

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
           if(ilo.le.ihi)then
             do j = domhi(2)+1, hi(2)
               do i = ilo,ihi
                 divu(i,j) = divu(i,domhi(2))
               enddo
             enddo
           endif
           if (lo(1).lt.domlo(1)) then
             do j = domhi(2)+1, hi(2)
               do i = lo(1), domlo(1)-1
                 divu(i,j) = divu(domlo(1),domhi(2))
               enddo
             enddo
           endif
           if(hi(1).gt.domhi(1))then
             do j = domhi(2)+1, hi(2)
               do i = domhi(1)+1, hi(1)
                 divu(i,j) = divu(domhi(1),domhi(2))
               enddo
             enddo
           endif
      endif            
      end

      subroutine FORT_DSDTFILL (dsdt,DIMS(dsdt),domlo,domhi,delta,
     &                         xlo,time,bc )

      integer    DIMDEC(dsdt)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     dsdt(DIMV(dsdt))

      integer    i, j
      integer    ilo, ihi, jlo, jhi

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(dsdt)
      hi(1) = ARG_H1(dsdt)
      lo(2) = ARG_L2(dsdt)
      hi(2) = ARG_H2(dsdt)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))

      call filcc (dsdt,DIMS(dsdt),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
           do i = lo(1), domlo(1)-1
             do j = lo(2), hi(2)
               dsdt(i,j) = zero
             enddo
           enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
           do i = domhi(1)+1, hi(1)
             do j = lo(2), hi(2)
               dsdt(i,j) = zero
             enddo
           enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
c                               inflow for burner in a can (bic, biac)

           do j = lo(2), domlo(2)-1
              do i = lo(1), hi(1)
                 dsdt(i,j) = zero
              enddo
           enddo
      endif            

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
           do j = domhi(2)+1, hi(2)
              do i = lo(1), hi(1)
                 dsdt(i,j) = zero
              enddo
           enddo
      endif            
      end

      subroutine FORT_YDOTFILL (ydot,DIMS(ydot),domlo,domhi,delta,
     &                          xlo,time,bc)

      integer    DIMDEC(ydot), bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     ydot(DIMV(ydot))

      integer    i, j
      integer    jlo, jhi, ilo, ihi
      REAL_T     x

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(ydot)
      hi(1) = ARG_H1(ydot)
      lo(2) = ARG_L2(ydot)
      hi(2) = ARG_H2(ydot)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))

      call filcc (ydot,DIMS(ydot),domlo,domhi,delta,xlo,bc)

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
c                               inflow for pipe
           do j = lo(2), domlo(2)-1
              do i = lo(1), hi(1)
                 ydot(i,j) = 0.0
              enddo
           enddo
      endif            

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
           do j = domhi(2)+1, hi(2)
              do i = lo(1), hi(1)
                 ydot(i,j) = 0.0
              enddo
           enddo
      endif            

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
           do i = lo(1), domlo(1)-1
             do j = lo(2), hi(2)
               ydot(i,j) = zero
             enddo
           enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
           do j = lo(2),hi(2)
                x = (float(i-lo(1))+.5)*delta(1)+xlo(1)
	       ydot(i,j) = zero
	    enddo
	 enddo
      endif            
      end

      subroutine FORT_FABMINMAX(lo, hi, 
     &                          fab, DIMS(fab),
     &                          fmin, fmax, nc)
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(fab)
      REAL_T  fab(DIMV(fab),nc)
      REAL_T  fmin, fmax

      integer i,j,n

      do n = 1,nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               fab(i,j,n) = MAX( fmin, MIN( fmax, fab(i,j,n) ) )
            end do
         end do
      end do
      end

      subroutine FORT_REPAIR_FLUX(lo, hi,
     &                            flux, DIMS(flux),
     &                            state, DIMS(state), nc, dir)
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM), nc, dir
      integer DIMDEC(flux)
      integer DIMDEC(state)
      REAL_T flux(DIMV(flux),nc)
      REAL_T state(DIMV(state),nc)
      
      integer i, j, n
      REAL_T sum, yedge, sumRhoY

      if (dir.eq.0) then
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               sum = zero
               sumRhoY = zero
               do n=1,nc
                  sum = sum + flux(i,j,n)
                  sumRhoY = sumRhoY + half*(state(i-1,j,n)+state(i,j,n))
               end do               
               do n=1,nc
                  yedge = half*(state(i-1,j,n)+state(i,j,n))
                  flux(i,j,n) = flux(i,j,n) - sum*yedge/sumRhoY
               end do
            end do
         end do
      else
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)
               sum = zero
               sumRhoY = zero
               do n=1,nc
                  sum = sum + flux(i,j,n)
                  sumRhoY = sumRhoY + half*(state(i,j-1,n)+state(i,j,n))
               end do
               do n=1,nc
                  yedge = half*(state(i,j-1,n)+state(i,j,n))
                  flux(i,j,n) = flux(i,j,n) - sum*yedge/sumRhoY
               end do
            end do
         end do
      end if
      end

      subroutine  FORT_RECOMP_UPDATE(lo, hi,
     &                               update, DIMS(update),
     &                               xflux,  DIMS(xflux),
     &                               yflux,  DIMS(yflux),
     &                               vol,    DIMS(vol),
     &                               nc)
      implicit none
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(update)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(vol)
      REAL_T update(DIMV(update),nc)
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T vol(DIMV(vol))

      integer i, j, n
      
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            do n=1,nc
               update(i,j,n)=-((xflux(i+1,j,n)-xflux(i,j,n))
     &              +          (yflux(i,j+1,n)-yflux(i,j,n)))/vol(i,j)
            end do
         end do
      end do
      end
      
      subroutine  FORT_INCRWEXTFLXDIV(lo, hi,
     &                                xflux,  DIMS(xflux),
     &                                yflux,  DIMS(yflux),
     &                                stateo, DIMS(stateo),
     &                                staten, DIMS(staten),
     &                                vol,    DIMS(vol),
     &                                nc, dt)
      implicit none
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(stateo)
      integer DIMDEC(staten)
      integer DIMDEC(vol)
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T stateo(DIMV(stateo))
      REAL_T staten(DIMV(staten))
      REAL_T vol(DIMV(vol))
      REAL_T dt

      integer i, j, n
      REAL_T dF
      
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            dF = zero
            do n=1,nc
               dF = dF + ( xflux(i+1,j,n) - xflux(i,j,n) )
     &              +    ( yflux(i,j+1,n) - yflux(i,j,n) ) 
            end do
            staten(i,j) = stateo(i,j) + dF*dt/vol(i,j)
         end do
      end do
      end
      
      subroutine  FORT_COMPUTE_UGRADP(p, DIMS(p), ugradp, DIMS(ugp),
     &                                umac,  DIMS(umac),
     &                                vmac,  DIMS(vmac),
     &                                lo, hi, dx)

      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(p)
      integer DIMDEC(ugp)
      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T      p(DIMV(p))
      REAL_T ugradp(DIMV(ugp))
      REAL_T dx(SDIM)

      integer i, j
      REAL_T uadv, vadv
      REAL_T p_x_lo, p_x_hi
      REAL_T p_y_lo, p_y_hi
      
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            uadv = half*(umac(i,j) + umac(i+1,j))
            vadv = half*(vmac(i,j) + vmac(i,j+1))
            p_x_hi = cvmgp(p(i  ,j),p(i+1,j),umac(i+1,j))
            p_x_lo = cvmgp(p(i-1,j),p(i  ,j),umac(i  ,j))
            p_y_hi = cvmgp(p(i,j  ),p(i,j+1),vmac(i,j+1))
            p_y_lo = cvmgp(p(i,j-1),p(i,j  ),vmac(i,j  ))
            ugradp(i,j) = uadv * (p_x_hi - p_x_lo) / dx(1) +
     &                    vadv * (p_y_hi - p_y_lo) / dx(2)      
         end do
      end do
      end

      integer function FORT_CONSERVATIVE_T_FLOOR(
     &     loF,hiF,state,DIMS(state),
     &     min_T, Tcomp, Rcomp, first_spec, last_spec, RhoH,
     &     ratio, tmp, nt)
      implicit none

      integer loF(SDIM), hiF(SDIM)
      integer DIMDEC(state)
      REAL_T  state(DIMV(state),0:*)      
      integer Tcomp, Rcomp, first_spec, last_spec, RhoH, ratio(SDIM), nt
      REAL_T  min_T, tmp(0:nt-1)
      integer n,i,j, loC(SDIM),hiC(SDIM),ii,jj,iii,jjj,ncells
      Real ncellsInv
      logical bad_T

c     Returns the number of fine cells fixed up
      FORT_CONSERVATIVE_T_FLOOR = 0

      ncells = 1
      do n=1,SDIM
         loC(n) = loF(n)/ratio(n)
         hiC(n) = (hiF(n)+1)/ratio(n) - 1
         ncells = ncells*ratio(n)
      enddo
      ncellsInv = 1.d0 / ncells

      do j=loC(2),hiC(2)
         do i=loC(1),hiC(1)
            
            bad_T = .false.
            do jj=0,ratio(2)-1
               jjj = ratio(2)*j + jj
               do ii=0,ratio(1)-1
                  iii = ratio(1)*i + ii
                  if (state(iii,jjj,Tcomp).lt.min_T) then
                     bad_T = .true.
                  endif
               enddo
            enddo
            
            if (bad_T .eqv. .true.) then

               tmp(Rcomp) = 0.d0
               do n=first_spec,last_spec
                  tmp(n) = 0.d0
               enddo
               tmp(RhoH) = 0.d0

               do jj=0,ratio(2)-1
                  jjj = ratio(2)*j + jj
                  do ii=0,ratio(1)-1
                     iii = ratio(1)*i + ii

                     tmp(Rcomp) = tmp(Rcomp) + state(iii,jjj,Rcomp)
                     do n=first_spec,last_spec
                        tmp(n) = tmp(n) + state(iii,jjj,n)
                     enddo
                     tmp(RhoH) = tmp(RhoH) + state(iii,jjj,RhoH)
                     
                  enddo
               enddo

               FORT_CONSERVATIVE_T_FLOOR = FORT_CONSERVATIVE_T_FLOOR + ncells

               tmp(Rcomp) = tmp(Rcomp) * ncellsInv
               do n=first_spec,last_spec
                  tmp(n) = tmp(n) * ncellsInv
               enddo
               tmp(RhoH) = tmp(RhoH)* ncellsInv
               
               do jj=0,ratio(2)-1
                  jjj = ratio(2)*j + jj
                  do ii=0,ratio(1)-1
                     iii = ratio(1)*i + ii
                     
                     state(iii,jjj,Rcomp) = tmp(Rcomp)
                     do n=first_spec,last_spec
                        state(iii,jjj,n) = tmp(n)
                     enddo
                     state(iii,jjj,RhoH) = tmp(RhoH)
                     
                  enddo
               enddo
               

            endif

         enddo
      enddo
      end

      subroutine FORT_HTDDRELAX(lo, hi, S, DIMS(S), yc, Tc, hc, rc,
     &     L, DIMS(L), a, DIMS(a), R, DIMS(R), dt, fac, maxRes, maxCor,
     &     for_T0_H1, res_only, mult)
      implicit none
#include "cdwrk.H"      
      integer lo(SDIM), hi(SDIM), yc, Tc, hc, rc
      integer DIMDEC(S)
      integer DIMDEC(L)
      integer DIMDEC(a)
      integer DIMDEC(R)
      REAL_T S(DIMV(s),0:*)
      REAL_T L(DIMV(L),0:*)
      REAL_T a(DIMV(a),0:*)
      REAL_T R(DIMV(R),0:*)
      REAL_T dt, fac(0:*), maxRes(0:*), maxCor(0:*)
      integer for_T0_H1, res_only
      REAL_T mult, cor
      integer i, j, n

#ifdef HT_SKIP_NITROGEN
      REAL_T, allocatable :: tfab(:,:)
#endif

c
c     Helper function to compute:
c     (a) dR = Rhs + (rho.phi - dt.Lphi)
c     (b) Res = Rhs - A(S) = Rhs - (rho.phi - dt.Lphi)
c               and then relax: phi = phi + Res/alpha
c     
c     NOTE: Assumes maxRes and maxCor have been initialized properly
c
      do n=0,Nspec-1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               L(i,j,n) = R(i,j,n) + mult*(S(i,j,rc)*S(i,j,yc+n) - dt*L(i,j,n))
               maxRes(n) = MAX(maxRes(n),ABS(L(i,j,n)))
            enddo
         enddo
      enddo
#ifdef HT_SKIP_NITROGEN
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            L(i,j,Nspec-1) = 0.d0
         enddo
      enddo
      maxRes(Nspec-1) = 0.d0
#endif

      if (for_T0_H1.eq.0) then
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               L(i,j,Nspec) = R(i,j,Nspec) + mult*(S(i,j,Tc) - dt*L(i,j,Nspec)/S(i,j,rc))
               maxRes(Nspec) = MAX(maxRes(Nspec),ABS(L(i,j,n)))
            enddo
         enddo
      else
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               L(i,j,Nspec) = R(i,j,Nspec) + mult*(S(i,j,rc)*S(i,j,hc) - dt*L(i,j,Nspec))
               maxRes(Nspec) = MAX(maxRes(Nspec),ABS(L(i,j,n)))
            enddo
         enddo
      endif

      if (res_only.ne.1) then
         do n=0,Nspec-1
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  cor = fac(n) * L(i,j,n) / (1.d0 + dt*a(i,j,n))
                  maxCor(n) = MAX(maxCor(n),ABS(cor))
                  S(i,j,yc+n) = S(i,j,yc+n)  +  cor
               enddo
            enddo
         enddo
#ifdef HT_SKIP_NITROGEN
         allocate(tfab(lo(1):hi(1),lo(2):hi(2)))
         tfab = 0.0d0
         do n=0,Nspec-2
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  tfab(i,j) = tfab(i,j) + S(i,j,yc+n)
               enddo
            enddo
         end do
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               S(i,j,Nspec-1) = 1.d0 - tfab(i,j)
            enddo
         enddo
         maxCor(Nspec-1) = 0.d0
         deallocate(tfab)
#endif

         if (for_T0_H1.eq.0) then
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  cor = fac(Nspec) * L(i,j,Nspec) / (1.d0 + dt*a(i,j,Nspec)/S(i,j,rc))
                  maxCor(Nspec) = MAX(maxCor(Nspec),ABS(cor))
                  S(i,j,Tc) = S(i,j,Tc)  +  cor
               enddo
            enddo
         else
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  cor = fac(Nspec) * L(i,j,Nspec) / (1.d0 + dt*a(i,j,Nspec))
                  maxCor(Nspec) = MAX(maxCor(Nspec),ABS(cor))
                  S(i,j,Tc) = S(i,j,Tc)  +  cor
               enddo
            enddo
         endif
      endif
      end
