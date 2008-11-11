c
c $Id: HEATTRANSFER_2D.F,v 1.14 2008-11-11 02:02:22 marc Exp $
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

      subroutine FORT_FLOOR_AND_SCALE_SPEC(nspecies,lo,hi,
     &  DIMS(spec),sum,spec)

      implicit none

      integer nspecies
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(spec)
      REAL_T  sum(DIMV(spec))      
      REAL_T  spec(DIMV(spec),nspecies)      

      integer i,j,k

      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          sum(i,j) = zero
        enddo
      enddo
      do k=1,nspecies
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            spec(i,j,k) = max(zero,spec(i,j,k))
            sum(i,j)    = sum(i,j)+spec(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nspecies
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            spec(i,j,k) = spec(i,j,k)/sum(i,j)
          enddo
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

      integer i, j, n

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

      subroutine FORT_TEMPVISC(lo,hi,
     &                         DIMS(T), T,
     &                         DIMS(Y), Y,
     &                         DIMS(lambda), lambda)
      implicit none
#include "visc.H"
#include "cdwrk.H"
      integer nspecies
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer DIMDEC(lambda)
      REAL_T  T(DIMV(T))
      REAL_T  Y(DIMV(Y),*)
      REAL_T  lambda(DIMV(lambda))
      
      integer i, j, n
      REAL_T Yt(maxspec)
      REAL_T cpmix
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1/
      data hi_chem /1,1/

      print *,'**** in FORT_TEMPVISC'
      call bl_abort()

      if (.not. LeEQ1) then
         call FORT_MIX_COND(lo, hi, lambda, DIMS(lambda),
     &        T, DIMS(T), Y, DIMS(Y))
      else
         call FORT_VELVISC(lo,hi,DIMS(T),T,DIMS(Y),Y,
     &        DIMS(lambda),lambda)
         do j=lo(2), hi(2)
            do i=lo(1), hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,n)
               end do
               CALL FORT_CPMIXfromTY(lo_chem, hi_chem, 
     &              CPMIX,  ARLIM(lo_chem), ARLIM(hi_chem),
     &              T(i,j), ARLIM(lo_chem), ARLIM(hi_chem),
     &              Yt,     ARLIM(lo_chem), ARLIM(hi_chem))
               lambda(i,j) = lambda(i,j) * cpmix / Pr
            end do
         end do
      end if
      end

      subroutine FORT_SPECVISC(lo,hi,
     &                         DIMS(T), T,
     &                         DIMS(Y), Y,
     &                         DIMS(rhoD), rhoD,
     &                         ncompd, P1ATM_MKS)
      implicit none
#include "visc.H"
#include "htdata.H"
#include "cdwrk.H"
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer DIMDEC(rhoD)
      integer ncompd
      REAL_T  T(DIMV(T))
      REAL_T  Y(DIMV(Y),*)
      REAL_T  rhoD(DIMV(rhoD),ncompd)
      
      integer i, j
      integer n
      REAL_T P1ATM_MKS, Ptmp, val, Patm

      print *,'**** in FORT_SPECVISC'
      call bl_abort()

      if (.not.use_constant_rhoD) then
         if (.not. LeEQ1) then
            Patm = pamb / P1ATM_MKS
            call FORT_MIXAVG_RHODIFF(lo, hi, rhoD, DIMS(rhoD), T, DIMS(T),
     &                               Y, DIMS(Y), Patm)
         else
            call FORT_VELVISC(lo,hi,DIMS(T),T,DIMS(Y),Y,
     &                        DIMS(rhoD),rhoD)
            do j=lo(2), hi(2)
               do i=lo(1), hi(1)
                  val = rhoD(i,j,1) / Sc
                  do n=1,ncompd
                     rhoD(i,j,n) = val
                  end do
               end do
            end do
         end if
      else
         do j=lo(2), hi(2)
            do i=lo(1), hi(1)
               do n=1,ncompd
                  rhoD(i,j,n) = constant_rhoD_val
               end do
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
#include "visc.H"
#include "htdata.H"
#include "cdwrk.H"
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
      REAL_T P1ATM_MKS, Ptmp, val, Patm
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
      REAL_T     y

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
      REAL_T     y

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
      REAL_T     y

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
      REAL_T     x, joff

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

      
#define VCOR    0
#define N2DOM   1
#define FINDDOM 2
#define JBBMOD  3
#define METH VCOR
      
#if !(METH == JBBMOD)
      
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
      
      character*(maxspnml) name
      integer i, j, n, pm, maxs, iN2
      REAL_T sum, val, maxv
      
#if (METH == VCOR)
      REAL_T yedge, sumRhoY
#endif

      iN2 = -1
      do n = 1,nc
         call get_spec_name(name,n)
         if (name .eq. 'N2' ) iN2 = n
      end do
         
      if (dir.eq.0) then
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               sum = zero
               sumRhoY = zero
               do n=1,nc
                  sum = sum + flux(i,j,n)
                  sumRhoY = sumRhoY + half*(state(i-1,j,n)+state(i,j,n))
               end do               
#if (METH == VCOR)
               do n=1,nc
                  yedge = half*(state(i-1,j,n)+state(i,j,n))
                  flux(i,j,n) = flux(i,j,n) - sum*yedge/sumRhoY
               end do
#else
#  if (METH == N2DOM)
               maxs = iN2
#  elif (METH == FINDDOM)
               pm = cvmgt(0,-1,sum.gt.zero)
               maxs=1
               maxv=state(i+pm,j,maxs)
               do n=2,nc
                  val = state(i+pm,j,n)
                  if (val.gt.maxv) then
                     maxv=val
                     maxs=n
                  end if
               end do
#  endif
               flux(i,j,maxs) = flux(i,j,maxs) - sum
#endif
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
#if (METH == VCOR)
               do n=1,nc
                  yedge = half*(state(i,j-1,n)+state(i,j,n))
                  flux(i,j,n) = flux(i,j,n) - sum*yedge/sumRhoY
               end do
#else
#  if (METH == N2DOM)
               maxs = iN2
#  elif (METH == FINDDOM)
               pm = cvmgt(0,-1,sum.gt.zero)
               maxs=1
               maxv=state(i,j+pm,maxs)
               do n=2,nc
                  val = state(i,j+pm,n)
                  if (val.gt.maxv) then
                     maxv=val
                     maxs=n
                  end if
               end do
#  endif
               flux(i,j,maxs) = flux(i,j,maxs) - sum
#endif
            end do
         end do
      end if
      end
#else
      subroutine FORT_REPAIR_FLUX(lo, hi,
     &                            flux, DIMS(flux),
     &                            state, DIMS(state), nc, dir)
      implicit none
      integer lo(SDIM), hi(SDIM), nc, dir
      integer DIMDEC(flux)
      integer DIMDEC(state)
      REAL_T flux(DIMV(flux),nc)
      REAL_T state(DIMV(state),nc)
#include "cdwrk.H"
      integer i,j,n,hi_e(SDIM),ii,jj
      REAL_T sump,sumn,delta,alphabar,sumy,beta,y(maxspec)
      character*(maxspnml) name
      
      if (nc.lt.Nspec) then
         call bl_abort('FORT_REPAIR_FLUX: too few flux comps for repair')
      end if

c      
c     The scheme:
c     We want the fluxes to sum to zero, take f+ and f- as +ve and -ve fluxes
c     at a face.  If the net flux is positive, we want alphabar*f+ + f- = 0.
c     We will reduce each fl+ (ie give flux back to originator) by an amount
c     alphal proportional to the mass frac of the downstream cell, thus
c                        (1 - alphal)*fl+ = beta*Yl
c     Summing over l, we get (1-alphabar)*sum(fl+) = beta*sum(Yl) as an equation
c     for beta given alphabar (alphabar computed from initial data so that 
c     flux sums to zero).  We exclude any fluxes from this mod if the downstream 
c     mass frac is <= 0.
c      
      if (dir .EQ. 0) then
         hi_e(1) = hi(1) + 1
         hi_e(2) = hi(2)
      else if (dir .EQ. 1) then
         hi_e(1) = hi(1)
         hi_e(2) = hi(2) + 1
      else
         write(6,*) 'FORT_REPAIR_FLUX: Bad dir'
      end if
      
      do j = lo(2),hi_e(2)
         do i = lo(1),hi_e(1)
            sump=zero
            sumn=zero
            do n=1,Nspec
               if (flux(i,j,n) .GT. zero) then
                  sump = sump + flux(i,j,n)
               else if (flux(i,j,n) .LT. zero) then
                  sumn = sumn + flux(i,j,n)
               end if
            end do
            delta = sump + sumn
            if (delta .NE. zero) then
               sumy = zero
               if (delta .LT. zero) then
                  do n=1,Nspec
                     if (dir .EQ. 0) then
                        y(n) = MAX(zero,state(i-1,j,n))
                     else
                        y(n) = MAX(zero,state(i,j-1,n))
                     endif
                     if (flux(i,j,n).GE.zero) y(n) = zero
                     sumy = sumy + y(n)
                  end do
                  alphabar = - sump/sumn
                  beta=cvmgt(sumn*(one-alphabar)/sumy,zero,sumy.NE.zero)
               else
                  do n=1,Nspec
                     y(n) = MAX(zero,state(i,j,n))
                     if (flux(i,j,n).LE.zero) y(n) = zero
                     sumy = sumy + y(n)
                  end do
                  alphabar = - sumn/sump
                  beta=cvmgt(sump*(one-alphabar)/sumy,zero,sumy.NE.zero)
               end if
               if (alphabar .GT. one) goto 100
               do n=1,Nspec
                  flux(i,j,n) = flux(i,j,n) - beta*y(n)
               end do
            end if
         end do
      end do
      return
      
 100  continue 
      write(6,*) 'FORT_REPAIR_FLUX: alphabar > 1'
      write(6,*) 'fluxes: (dir=',dir,')'
      do n=1,Nspec
         call get_spec_name(name,n)
         write(6,996) name,flux(i,j,n)
      end do
      write(6,*) 'left state:'
      do n=1,Nspec
         call get_spec_name(name,n)
         write(6,996) name,state(i,j,n)
      end do
      write(6,*) 'right state:'
      ii = i
      jj = j
      if (dir .EQ. 0) then
         ii = i + 1
      else
         jj = j + 1
      end if
      do n=1,Nspec
         call get_spec_name(name,n)
         write(6,996) name,state(ii,jj,n)
      end do
 996  format(a16,1x,e30.22)
      call bl_abort(" ")
      end
#endif
      
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

c///
c///  start of FORT_RK_MIXTURE_AVERAGED
c///

c     evaluate the "operator" (extensive fluxes and their divergences) for the 
c     Runge-Kutta implementation of the diffusion update using the mixture
c     averaged formulation.

c     the divergences are scaled by 1/vol so they do approximate the analytic
c     divergence.  when used in the Runge-Kutta formula their additional scaling,
c     given by scale_factor, should be set to - dt.  the minus occurs because
c     the standard ODE is written y_prime = f(y), that is, the divergences are
c     moved to the opposite side of the equation from the time derivative.  when
c     used in getViscTerms the scale_factor should be set to -1.

c     to do:
c     1) must adjust diffusion coefficients by thickness factors
c     2) unsure of LMC support for r-z coordinates

c     arguments are alphabetical, mostly:
c
c     domain_lo, domain_hi,             ! INPUT limits of valid region of the domain
c     lo, hi,                           ! INPUT limits of valid region of the box
c   * areax, DIMS(areax),               ! INPUT areas of the faces perendicular to x axis
c   * areay, DIMS(areay),               ! INPUT areas of the faces perpendicular to y axis
c     bc,                               ! INPUT boundary condition array for all comps
c     dx,                               ! INPUT physical dimensions of grid cells
c     index_of_firstY,                  ! INPUT index of rho Y for the first species in the state
c     index_of_lastY,                   ! INPUT index of rho Y for the last species in the state
c     index_of_rho,                     ! INPUT index of rho in the state
c     index_of_rhoH,                    ! INPUT index of rho H in the state
c     index_of_T,                       ! INPUT index of T in the state
c     maximum_error,                    ! INPUT maximum error in calculation of T
c     maximum_iterations,               ! INPUT maximum iterations in calculation of T
c     ncomps,                           ! INPUT total number of components in the state
c     nspecies,                         ! INPUT total number of species in the state
c     scale_factor,                     ! INPUT scale factor for the divergences of the fluxes
c     state, DIMS(state),               ! INPUT all variables in the state
c   * volume, DIMS(volume),             ! INPUT volumes of the cells
c     div_of_flux_for_H, DIMS(div_of_flux_for_H), ! OUTPUT divergence of the flux for rho H
c     div_of_flux_for_Y, DIMS(div_of_flux_for_Y), ! OUTPUT divergences of the fluxes for rho Y
c     xflux_for_H, DIMS(xflux_for_H),   ! OUTPUT extensive x fluxes for rho H
c     xflux_for_Y, DIMS(xflux_for_Y),   ! OUTPUT extensive x fluxes for rho Y
c     yflux_for_H, DIMS(yflux_for_H),   ! OUTPUT extensive y fluxes for rho H
c     yflux_for_Y, DIMS(yflux_for_Y),   ! OUTPUT extensive y fluxes for rho Y
c
c     * these arguments are not used

      subroutine FORT_RK_MIXTURE_AVERAGED
     &     (
     &     domain_lo, domain_hi,
     &     lo, hi,
     &     areax, DIMS(areax),
     &     areay, DIMS(areay),
     &     bc,
     &     dx,
     &     index_of_firstY,
     &     index_of_lastY,
     &     index_of_rho,
     &     index_of_rhoH,
     &     index_of_T,
     &     maximum_error,
     &     maximum_iterations,
     &     ncomps,
     &     nspecies,
     &     scale_factor,
     &     state, DIMS(state),
     &     volume, DIMS(volume),
     &     div_of_flux_for_H, DIMS(div_of_flux_for_H),
     &     div_of_flux_for_Y, DIMS(div_of_flux_for_Y),
     &     xflux_for_H, DIMS(xflux_for_H),
     &     xflux_for_Y, DIMS(xflux_for_Y),
     &     yflux_for_H, DIMS(yflux_for_H),
     &     yflux_for_Y, DIMS(yflux_for_Y)
     &     )

      implicit none

#include "BC_TYPES.H"
#include "cdwrk.H"

c     arguments
      integer ncomps
      integer nspecies
      integer domain_lo(BL_SPACEDIM), domain_hi(BL_SPACEDIM)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer DIMDEC(areax)
      integer DIMDEC(areay)
      integer bc(BL_SPACEDIM,2,ncomps)
      integer DIMDEC(state)
      integer DIMDEC(div_of_flux_for_H)
      integer DIMDEC(div_of_flux_for_Y)
      integer DIMDEC(xflux_for_H)
      integer DIMDEC(yflux_for_H)
      integer DIMDEC(xflux_for_Y)
      integer DIMDEC(yflux_for_Y)
      integer DIMDEC(volume)
      integer index_of_firstY
      integer index_of_lastY
      integer index_of_rho
      integer index_of_rhoH
      integer index_of_T
      integer maximum_iterations
      REAL_T areax(DIMV(areax))
      REAL_T areay(DIMV(areay))
      REAL_T dx(BL_SPACEDIM)
      REAL_T maximum_error
      REAL_T div_of_flux_for_H(DIMV(div_of_flux_for_H))
      REAL_T div_of_flux_for_Y(DIMV(div_of_flux_for_Y),nspecies)
      REAL_T scale_factor
      REAL_T state(DIMV(state),ncomps)
      REAL_T volume(DIMV(volume))
      REAL_T xflux_for_H(DIMV(xflux_for_H))
      REAL_T yflux_for_H(DIMV(yflux_for_H))
      REAL_T xflux_for_Y(DIMV(xflux_for_Y),ncomps)
      REAL_T yflux_for_Y(DIMV(yflux_for_Y),ncomps)

c     local variables
      external FORT_TfromHYpt
      integer FORT_TfromHYpt
      integer ell
      integer i
      integer iterations
      integer j
      integer n
      integer ni
      integer nj
      integer offset
      logical at_hi
      logical at_lo
      logical error
      logical hi_dirichlet_check
      logical hi_dirichlet_for_T
      logical hi_dirichlet_for_Y
      logical hi_internal
      logical hi_internal_check
      logical lo_dirichlet_check
      logical lo_dirichlet_for_T
      logical lo_dirichlet_for_Y
      logical lo_internal
      logical lo_internal_check
      REAL_T alpha
      REAL_T area_for_x
      REAL_T area_for_y
      REAL_T ambient_pressure_cgs
      REAL_T ambient_pressure_si
      REAL_T coefficient
      REAL_T derivative_of_T
      REAL_T dpdt_factor
      REAL_T FORT_P1ATMMKS
      REAL_T h
      REAL_T hi_coefficient_0
      REAL_T hi_coefficient_1
      REAL_T hi_coefficient_2
      REAL_T lambda_on_edge_cgs
      REAL_T lambda_on_edge_si
      REAL_T lo_coefficient_0
      REAL_T lo_coefficient_1
      REAL_T lo_coefficient_2
      REAL_T minus_value_cgs
      REAL_T mmw_on_edge
      REAL_T one_atmosphere_cgs
      REAL_T one_atmosphere_si
      REAL_T one_over_dx
      REAL_T one_over_dy
      REAL_T one_over_rho
      REAL_T plus_value_cgs
      REAL_T RUC
      REAL_T scale_factor_over_volume
      REAL_T sum_h_xflux_for_Y
      REAL_T sum_h_yflux_for_Y
      REAL_T sum_of_fluxes
      REAL_T sum_of_rho_Y
      REAL_T T_on_edge
      REAL_T universal_gas_constant_cgs

#ifdef BL_CYGWIN_NT
      print *,'Need to rebuild with f90'
      stop
#else
c     local arrays
      REAL_T, allocatable :: contiguous_Y (:)
      REAL_T, allocatable :: cpms_on_edge (:)
      REAL_T, allocatable :: derivative_of_Y (:)
      REAL_T, allocatable :: hms_on_edge_cgs (:)
      REAL_T, allocatable :: hms_on_edge_si (:)
      REAL_T, allocatable :: one_over_wt (:)
      REAL_T, allocatable :: residual (:)
      REAL_T, allocatable :: rho_Dmix_on_edge_si (:)
      REAL_T, allocatable :: rho_Upsilon_on_edge_cgs (:)
      REAL_T, allocatable :: rho_Y (:,:,:)
      REAL_T, allocatable :: rho_Y_on_edge (:)
      REAL_T, allocatable :: T (:,:)
      REAL_T, allocatable :: wt (:)
      REAL_T, allocatable :: X (:,:,:)
      REAL_T, allocatable :: X_on_edge (:)
      REAL_T, allocatable :: Y (:,:,:)
      REAL_T, allocatable :: Y_on_edge (:)

      allocate (contiguous_Y (nspecies))
      allocate (cpms_on_edge (nspecies))
      allocate (derivative_of_Y (nspecies))
      allocate (hms_on_edge_cgs (nspecies))
      allocate (hms_on_edge_si (nspecies))
      allocate (one_over_wt (nspecies))
      allocate (residual (maximum_iterations))
      allocate (rho_Dmix_on_edge_si (nspecies))
      allocate (rho_Upsilon_on_edge_cgs (nspecies))
      allocate (rho_Y (DIMV(state), nspecies))
      allocate (rho_Y_on_edge (nspecies))
      allocate (T (DIMV(state)))
      allocate (wt (nspecies))
      allocate (X (DIMV(state), nspecies))
      allocate (X_on_edge (nspecies))
      allocate (Y (DIMV(state), nspecies))
      allocate (Y_on_edge (nspecies))

c///
c///  initialization
c///

c     offset of species in the state
      offset = index_of_firstY - 1

c     some constants, at least in Cartesian coordinates
      area_for_x = dx(2)
      area_for_y = dx(1)
      one_over_dx = 1.0d0 / dx(1)
      one_over_dy = 1.0d0 / dx(2)
      scale_factor_over_volume = scale_factor / (dx(1) * dx(2))

c     coefficients for difference formulas at the hi and low edges
c     note still must divide by dx or dy to get the approximate derivative
      hi_coefficient_0 = + 1.0d0 / 3.0d0
      hi_coefficient_1 = - 3.0d0
      hi_coefficient_2 = + 8.0d0 / 3.0d0
      lo_coefficient_0 = - 8.0d0 / 3.0d0
      lo_coefficient_1 = + 3.0d0
      lo_coefficient_2 = - 1.0d0 / 3.0d0

c     get ambient pressure
c     one_atmosphere_cgs = 1013250.0 (dyne / cm**2) = 10.0 (Pascal)
c     one_atmosphere_si = 101325.0 (Pascal)
      call CKRP (IWRK(ckbi), RWRK(ckbr), universal_gas_constant_cgs, RUC, one_atmosphere_cgs)
      call FORT_GETPAMB (ambient_pressure_si, dpdt_factor)
      one_atmosphere_si = FORT_P1ATMMKS ()
      ambient_pressure_cgs = one_atmosphere_cgs * ambient_pressure_si / one_atmosphere_si

c     get reciprocals of molecular weights
      call CKWT (IWRK(ckbi), RWRK(ckbr), wt)
      do n=1,nspecies
         one_over_wt(n) = 1.0d0 / wt(n)
      end do

c     form rho_Y, T, X, and Y from rho, rho H and rho Y in the state
      do j=ARG_L2(state),ARG_H2(state)
         do i=ARG_L1(state),ARG_H1(state)
            one_over_rho = 1.0d0 / state(i, j, index_of_rho)
            h = one_over_rho * state(i, j, index_of_rhoH)
            do n=1,nspecies
               rho_Y(i, j, n) = state(i, j, offset + n)
               Y(i, j, n) = one_over_rho * state(i, j, offset + n)
            end do
            do n=1,nspecies
               contiguous_Y(n) = Y(i, j, n)
            end do
c           use the state value for T as a guess for the Newton iteration
            T(i, j) = state(i, j, index_of_T)
            iterations = FORT_TfromHYpt
     &           (T(i, j),
     &           h, 
     &           contiguous_Y,
     &           maximum_error,
     &           maximum_iterations,
     &           residual)
         end do
      end do

c///
c///  x fluxes
c///

c     set flags for special treatment of the lo boundary

      lo_dirichlet_for_T = .false.
      lo_dirichlet_for_Y = .false.
      lo_internal        = .true.
      if (domain_lo(1) .eq. lo(1)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         lo_dirichlet_for_T = EXT_DIR .eq. bc(1, 1, index_of_T)
         lo_dirichlet_for_Y = EXT_DIR .eq. bc(1, 1, index_of_firstY)
         lo_dirichlet_check = EXT_DIR .eq. bc(1, 1, index_of_firstY)
         lo_internal        = INT_DIR .eq. bc(1, 1, index_of_T)
         lo_internal_check  = INT_DIR .eq. bc(1, 1, index_of_T)
         do n=1,nspecies
            lo_dirichlet_for_Y = lo_dirichlet_for_Y .and. (EXT_DIR .eq. bc(1, 1, offset + n))
            lo_dirichlet_check = lo_dirichlet_check  .or. (EXT_DIR .eq. bc(1, 1, offset + n))
            lo_internal        = lo_internal        .and. (INT_DIR .eq. bc(1, 1, offset + n))
            lo_internal_check  = lo_internal_check   .or. (INT_DIR .eq. bc(1, 1, offset + n))
         end do
         error = .not. (
     &        (lo_dirichlet_for_Y .eqv. lo_dirichlet_check) .and.
     &        (lo_internal .eqv. lo_internal_check) .and.
     &        (.not. lo_dirichlet_for_Y .or. lo_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_MIXTURE_AVERAGED: inconsistent bc array at low x")
      end if

c     set flags for special treatment of the hi boundary

      hi_dirichlet_for_T = .false.
      hi_dirichlet_for_Y = .false.
      hi_internal        = .true.
      if (domain_hi(1) .eq. hi(1)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         hi_dirichlet_for_T = EXT_DIR .eq. bc(1, 2, index_of_T)
         hi_dirichlet_for_Y = EXT_DIR .eq. bc(1, 2, index_of_firstY)
         hi_dirichlet_check = EXT_DIR .eq. bc(1, 2, index_of_firstY)
         hi_internal        = INT_DIR .eq. bc(1, 2, index_of_T)
         hi_internal_check  = INT_DIR .eq. bc(1, 2, index_of_T)
         do n=1,nspecies
            hi_dirichlet_for_Y = hi_dirichlet_for_Y .and. (EXT_DIR .eq. bc(1, 2, offset + n))
            hi_dirichlet_check = hi_dirichlet_check  .or. (EXT_DIR .eq. bc(1, 2, offset + n))
            hi_internal        = hi_internal        .and. (INT_DIR .eq. bc(1, 2, offset + n))
            hi_internal_check  = hi_internal_check   .or. (INT_DIR .eq. bc(1, 2, offset + n))
         end do
         error = .not. (
     &        (hi_dirichlet_for_Y .eqv. hi_dirichlet_check) .and.
     &        (hi_internal .eqv. hi_internal_check) .and.
     &        (.not. hi_dirichlet_for_Y .or. hi_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_MIXTURE_AVERAGED: inconsistent bc array at high x")
      end if

c     calculate fluxes in the x direction

      do j=lo(2),hi(2)
c        fluxes in x direction are at edges indexed from lo(1) to hi(1)+1
         do i=lo(1),hi(1)+1

c           set flags indicating whether the edge is at the left or right boundary
            at_lo = i .eq. domain_lo(1)
            at_hi = i .eq. domain_hi(1) + 1

c           evaluate the value for T on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            if (at_lo .and. lo_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i-1, j, index_of_T)
            else if (at_hi .and. hi_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i, j, index_of_T)
            else
               T_on_edge = 0.5d0 * (T(i-1,j) + T(i,j))
            end if

c           use a floored T just in case
            T_on_edge = MAX (TMIN_TRANS, T_on_edge)

c           evaluate the values for Y and rho_Y on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            do n = 1,nspecies
               if (at_lo .and. lo_dirichlet_for_Y) then
                  rho_Y_on_edge(n) = rho_Y(i-1, j, n)
                  Y_on_edge(n) = Y(i-1, j, n)
               else if (at_hi .and. hi_dirichlet_for_Y) then
                  rho_Y_on_edge(n) = rho_Y(i, j, n)
                  Y_on_edge(n) = Y(i, j, n)
               else
                  rho_Y_on_edge(n) = 0.5d0 * (rho_Y(i-1,j,n) + rho_Y(i,j,n))
                  Y_on_edge(n) = 0.5d0 * (Y(i-1,j,n) + Y(i,j,n))
               end if
            end do

c           evaluate the values for X on the edge
            call CKYTX 
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), X_on_edge)

c           evaluate the derivative for T on the edge
            if (at_lo .and. .not. lo_internal) then
               if (lo_dirichlet_for_T) then
                  derivative_of_T 
     &                 = one_over_dx
     &                 * (lo_coefficient_0 * T_on_edge
     &                 + lo_coefficient_1 * T(i,j)
     &                 + lo_coefficient_2 * T(i+1,j))
               else
                  derivative_of_T = 0.0d0
               end if
            else if (at_hi .and. .not. hi_internal) then
               if (hi_dirichlet_for_T) then
                  derivative_of_T 
     &                 = one_over_dx
     &                 * (hi_coefficient_0 * T(i-2,j)
     &                 + hi_coefficient_1 * T(i-1,j)
     &                 + hi_coefficient_2 * T_on_edge)
               else
                  derivative_of_T = 0.0d0
               end if
            else
               derivative_of_T = (T(i,j) - T(i-1,j)) * one_over_dx
            end if

c           evaluate the derivatives for Y on the edge
            do n = 1,nspecies
               if (at_lo .and. .not. lo_internal) then
                  if (lo_dirichlet_for_Y) then
                     derivative_of_Y(n)
     &                    = one_over_dx
     &                    * (lo_coefficient_0 * Y_on_edge(n)
     &                    + lo_coefficient_1 * Y(i,j,n)
     &                    + lo_coefficient_2 * Y(i+1,j,n))
                  else
                     derivative_of_Y(n) = 0.0d0
                  end if
               else if (at_hi .and. .not. hi_internal) then
                  if (hi_dirichlet_for_Y) then
                     derivative_of_Y(n)
     &                    = one_over_dx
     &                    * (hi_coefficient_0 * Y(i-2,j,n)
     &                    + hi_coefficient_1 * Y(i-1,j,n)
     &                    + hi_coefficient_2 * Y_on_edge(n))
                  else
                     derivative_of_Y(n) = 0.0d0
                  end if
               else
                  derivative_of_Y(n) = (Y(i,j,n) - Y(i-1,j,n)) * one_over_dx
               end if
            end do

c           evaluate diffusion coefficients and other things on the edge
            call CKCPMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), cpms_on_edge)
            call EGSPAR 
     &           (T_on_edge, X_on_edge, Y_on_edge, cpms_on_edge, 
     &           RWRK(egbr), IWRK(egbi))
            call CKMMWY 
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), mmw_on_edge)
            CALL EGSVR1 
     &           (T_on_edge, Y_on_edge, RWRK(egbr), rho_Upsilon_on_edge_cgs)
            call CKHMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms_on_edge_cgs)
c           the following averaging is equation 50 in the chemkin transport package 
c           manual.  the EGLIB manual does not mention the averaging, but it does have
c           provisions to compute it through the alpha parameter in its equation 17.
            alpha = 1.0d0
            CALL EGSL1 (alpha, T_on_edge, X_on_edge, RWRK(egbr), plus_value_cgs)
            alpha = - 1.0d0
            CALL EGSL1 (alpha, T_on_edge, X_on_edge, RWRK(egbr), minus_value_cgs)
            lambda_on_edge_cgs = 0.5d0 * (plus_value_cgs + minus_value_cgs)

c           convert the library values from cgs units to si units

c           per page 98 of the chemkin manual, the enthalpies for the species,
c           hms_on_edge_cgs, have units
c
c                            erg / g
c                 =          (J / 10**7) / (Kg / 10**3)
c                 = 1.0d-4 * J / Kg
c
c           the si unit of thermal conductivity is J / Kg, so the EGLIB values must
c           be multiplied by 1.0d-4 to convert them to si units.
            do n=1,nspecies
               hms_on_edge_si(n) = 1.0d-4 * hms_on_edge_cgs(n)
            end do

c           per page 23 in the EGLIB manual, the rescaled flux diffusion vector,
c           rho_Upsilon_cgs, has units
c
c                            g / (cm s)
c                 =          (Kg / 10**3) / ((m / 10**2) s)
c                 = 1.0d-1 * Kg / (m s)
c
c           so the EGLIB values must be multiplied by 1.0d-1 to convert them to si 
c           units. per equation 12 in the EGLIB manual, the scaled mixture averaged 
c           diffusion coefficient rho_Dmix is rho_Upsilon times the dimensionless 
c           mmw / wt(n).  
            do n=1,nspecies
                rho_Dmix_on_edge_si(n) 
     &              = 1.0d-1 *  rho_Upsilon_on_edge_cgs(n) * mmw_on_edge * one_over_wt(n)
            end do

c           per page 15 of the EGLIB manual, these thermal conductivities are returned
c           with units (erg / (cm s K)) = ((J / 10**7) / ((m / 10**2) s K)) =
c           = ((W s) / (10**5 m s K)) = (W / (m K 10**5)).  the si unit of thermal
c           conductivity is (W / (M K)), so these values must be divided by 10**5
c           to convert them to si units.
            lambda_on_edge_si = 1.0d-5 * lambda_on_edge_cgs

c           assemble the fluxes
            sum_h_xflux_for_Y = 0
            do n = 1,nspecies
               xflux_for_Y(i,j,n)
     &              = - area_for_x
     &              * rho_Dmix_on_edge_si(n)
     &              * derivative_of_Y(n)
               sum_h_xflux_for_Y
     &              = sum_h_xflux_for_Y
     &              + hms_on_edge_si(n)
     &              * xflux_for_Y(i,j,n)
            end do
            xflux_for_H(i,j)
     &           = sum_h_xflux_for_Y
     &           -  area_for_x
     &           * lambda_on_edge_si 
     &           * derivative_of_T

c           flux correction
            sum_of_fluxes = 0
            sum_of_rho_Y = 0
            do n=1,nspecies
               sum_of_fluxes = sum_of_fluxes + xflux_for_Y(i,j,n)
               sum_of_rho_Y = sum_of_rho_Y + rho_Y_on_edge(n)
            end do
            coefficient = sum_of_fluxes / sum_of_rho_Y
            do n=1,nspecies
               xflux_for_Y(i,j,n) 
     &              = xflux_for_Y(i,j,n) 
     &              - rho_Y_on_edge(n)
     &              * coefficient
            end do

         end do
      end do

c///
c///  y fluxes
c///

c     set flags for special treatment of the lo boundary

      lo_dirichlet_for_T = .false.
      lo_dirichlet_for_Y = .false.
      lo_internal        = .true.
      if (domain_lo(2) .eq. lo(2)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         lo_dirichlet_for_T = EXT_DIR .eq. bc(2, 1, index_of_T)
         lo_dirichlet_for_Y = EXT_DIR .eq. bc(2, 1, index_of_firstY)
         lo_dirichlet_check = EXT_DIR .eq. bc(2, 1, index_of_firstY)
         lo_internal        = INT_DIR .eq. bc(2, 1, index_of_T)
         lo_internal_check  = INT_DIR .eq. bc(2, 1, index_of_T)
         do n=1,nspecies
            lo_dirichlet_for_Y = lo_dirichlet_for_Y .and. (EXT_DIR .eq. bc(2, 1, offset + n))
            lo_dirichlet_check = lo_dirichlet_check  .or. (EXT_DIR .eq. bc(2, 1, offset + n))
            lo_internal        = lo_internal        .and. (INT_DIR .eq. bc(2, 1, offset + n))
            lo_internal_check  = lo_internal_check   .or. (INT_DIR .eq. bc(2, 1, offset + n))
         end do
         error = .not. (
     &        (lo_dirichlet_for_Y .eqv. lo_dirichlet_check) .and.
     &        (lo_internal .eqv. lo_internal_check) .and.
     &        (.not. lo_dirichlet_for_Y .or. lo_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_MIXTURE_AVERAGED: inconsistent bc array at low y")
      end if

c     set flags for special treatment of the hi boundary

      hi_dirichlet_for_T = .false.
      hi_dirichlet_for_Y = .false.
      hi_internal        = .true.
      if (domain_hi(2) .eq. hi(2)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         hi_dirichlet_for_T = EXT_DIR .eq. bc(2, 2, index_of_T)
         hi_dirichlet_for_Y = EXT_DIR .eq. bc(2, 2, index_of_firstY)
         hi_dirichlet_check = EXT_DIR .eq. bc(2, 2, index_of_firstY)
         hi_internal        = INT_DIR .eq. bc(2, 2, index_of_T)
         hi_internal_check  = INT_DIR .eq. bc(2, 2, index_of_T)
         do n=1,nspecies
            hi_dirichlet_for_Y = hi_dirichlet_for_Y .and. (EXT_DIR .eq. bc(2, 2, offset + n))
            hi_dirichlet_check = hi_dirichlet_check  .or. (EXT_DIR .eq. bc(2, 2, offset + n))
            hi_internal        = hi_internal        .and. (INT_DIR .eq. bc(2, 2, offset + n))
            hi_internal_check  = hi_internal_check   .or. (INT_DIR .eq. bc(2, 2, offset + n))
         end do
         error = .not. (
     &        (hi_dirichlet_for_Y .eqv. hi_dirichlet_check) .and.
     &        (hi_internal .eqv. hi_internal_check) .and.
     &        (.not. hi_dirichlet_for_Y .or. hi_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_MIXTURE_AVERAGED: inconsistent bc array at high y")
      end if

c     calculate fluxes in the y direction

c     fluxes in y direction are at edges indexed from lo(2) to hi(2)+1
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)

c           set flags indicating whether the edge is at the left or right boundary
            at_lo = j .eq. domain_lo(2)
            at_hi = j .eq. domain_hi(2) + 1

c           evaluate the value for T on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            if (at_lo .and. lo_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i, j-1, index_of_T)
            else if (at_hi .and. hi_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i, j, index_of_T)
            else
               T_on_edge = 0.5d0 * (T(i,j-1) + T(i,j))
            end if

c           use a floored T just in case
            T_on_edge = MAX (TMIN_TRANS, T_on_edge)

c           evaluate the values for Y on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            do n = 1,nspecies
               if (at_lo .and. lo_dirichlet_for_Y) then
                  Y_on_edge(n) = Y(i, j-1, n)
               else if (at_hi .and. hi_dirichlet_for_Y) then
                  Y_on_edge(n) = Y(i, j, n)
               else
                  Y_on_edge(n) = 0.5d0 * (Y(i,j-1,n) + Y(i,j,n))
               end if
            end do

c           evaluate the values for Y and rho_Y on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            do n = 1,nspecies
               if (at_lo .and. lo_dirichlet_for_Y) then
                  rho_Y_on_edge(n) = rho_Y(i, j-1, n)
                  Y_on_edge(n) = Y(i, j-1, n)
               else if (at_hi .and. hi_dirichlet_for_Y) then
                  rho_Y_on_edge(n) = rho_Y(i, j, n)
                  Y_on_edge(n) = Y(i, j, n)
               else
                  rho_Y_on_edge(n) = 0.5d0 * (rho_Y(i,j-1,n) + rho_Y(i,j,n))
                  Y_on_edge(n) = 0.5d0 * (Y(i,j-1,n) + Y(i,j,n))
               end if
            end do

c           evaluate the values for X on the edge
            call CKYTX
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), X_on_edge)

c           evaluate the derivative for T on the edge
            if (at_lo .and. .not. lo_internal) then
               if (lo_dirichlet_for_T) then
                  derivative_of_T
     &                 = one_over_dy
     &                 * (lo_coefficient_0 * T_on_edge
     &                 + lo_coefficient_1 * T(i,j)
     &                 + lo_coefficient_2 * T(i,j+1))
               else
                  derivative_of_T = 0
               end if
            else if (at_hi .and. .not. hi_internal) then
               if (hi_dirichlet_for_T) then
                  derivative_of_T
     &                 = one_over_dy
     &                 * (hi_coefficient_0 * T(i,j-2)
     &                 + hi_coefficient_1 * T(i,j-1)
     &                 + hi_coefficient_2 * T_on_edge)
               else
                  derivative_of_T = 0
               end if
            else
               derivative_of_T = (T(i,j) - T(i,j-1)) * one_over_dy
            end if

c           evaluate the derivatives for X on the edge
            do n = 1,nspecies
               if (at_lo .and. .not. lo_internal) then
                  if (lo_dirichlet_for_Y) then
                     derivative_of_Y(n)
     &                    = one_over_dy
     &                    * (lo_coefficient_0 * Y_on_edge(n)
     &                    + lo_coefficient_1 * Y(i,j,n)
     &                    + lo_coefficient_2 * Y(i,j+1,n))
                  else
                     derivative_of_Y(n) = 0.0d0
                  end if                     
               else if (at_hi .and. .not. hi_internal) then
                  if (hi_dirichlet_for_Y) then
                     derivative_of_Y(n)
     &                    = one_over_dy
     &                    * (hi_coefficient_0 * Y(i,j-2,n)
     &                    + hi_coefficient_1 * Y(i,j-1,n)
     &                    + hi_coefficient_2 * Y_on_edge(n))
                  else
                     derivative_of_Y(n) = 0.0d0
                  end if
               else
                  derivative_of_Y(n) = (Y(i,j,n) - Y(i,j-1,n)) * one_over_dy
               end if
            end do

c           evaluate diffusion coefficients and other things on the edge
            call CKCPMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), cpms_on_edge)
            call EGSPAR 
     &           (T_on_edge, X_on_edge, Y_on_edge, cpms_on_edge, 
     &           RWRK(egbr), IWRK(egbi))
            call CKMMWY 
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), mmw_on_edge)
            CALL EGSVR1 
     &           (T_on_edge, Y_on_edge, RWRK(egbr), rho_Upsilon_on_edge_cgs)
            call CKHMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms_on_edge_cgs)
c           the following averaging is equation 50 in the chemkin transport package 
c           manual.  the EGLIB manual does not mention the averaging, but it does have
c           provisions to compute it through the alpha parameter in its equation 17.
            alpha = 1.0d0
            CALL EGSL1 (alpha, T_on_edge, X_on_edge, RWRK(egbr), plus_value_cgs)
            alpha = - 1.0d0
            CALL EGSL1 (alpha, T_on_edge, X_on_edge, RWRK(egbr), minus_value_cgs)
            lambda_on_edge_cgs = 0.5d0 * (plus_value_cgs + minus_value_cgs)

c           convert the library values from cgs units to si units

c           per page 98 of the chemkin manual, the enthalpies for the species,
c           hms_on_edge_cgs, have units
c
c                            erg / g
c                 =          (J / 10**7) / (Kg / 10**3)
c                 = 1.0d-4 * J / Kg
c
c           the si unit of thermal conductivity is J / Kg, so the EGLIB values must
c           be multiplied by 1.0d-4 to convert them to si units.
            do n=1,nspecies
               hms_on_edge_si(n) = 1.0d-4 * hms_on_edge_cgs(n)
            end do

c           per page 23 in the EGLIB manual, the rescaled flux diffusion vector,
c           rho_Upsilon_cgs, has units
c
c                            g / (cm s)
c                 =          (Kg / 10**3) / ((m / 10**2) s)
c                 = 1.0d-1 * Kg / (m s)
c
c           so the EGLIB values must be multiplied by 1.0d-1 to convert them to si 
c           units. per equation 12 in the EGLIB manual, the scaled mixture averaged 
c           diffusion coefficient rho_Dmix is rho_Upsilon times the dimensionless 
c           mmw / wt(n).  
            do n=1,nspecies
                rho_Dmix_on_edge_si(n) 
     &              = 1.0d-1 *  rho_Upsilon_on_edge_cgs(n) * mmw_on_edge * one_over_wt(n)
            end do

c           per page 15 of the EGLIB manual, these thermal conductivities are returned
c           with units (erg / (cm s K)) = ((J / 10**7) / ((m / 10**2) s K)) =
c           = ((W s) / (10**5 m s K)) = (W / (m K 10**5)).  the si unit of thermal
c           conductivity is (W / (M K)), so these values must be divided by 10**5
c           to convert them to si units.
            lambda_on_edge_si = 1.0d-5 * lambda_on_edge_cgs

c           assemble the fluxes
            sum_h_yflux_for_Y = 0
            do n = 1,nspecies
               yflux_for_Y(i,j,n)
     &              = - area_for_x
     &              * rho_Dmix_on_edge_si(n)
     &              * derivative_of_Y(n)
               sum_h_yflux_for_Y
     &              = sum_h_yflux_for_Y
     &              + hms_on_edge_si(n)
     &              * yflux_for_Y(i,j,n)
            end do
            yflux_for_H(i,j)
     &           = sum_h_yflux_for_Y
     &           - area_for_y
     &           * lambda_on_edge_si 
     &           * derivative_of_T

c           flux correction
            sum_of_fluxes = 0
            sum_of_rho_Y = 0
            do n=1,nspecies
               sum_of_fluxes = sum_of_fluxes + yflux_for_Y(i,j,n)
               sum_of_rho_Y = sum_of_rho_Y + rho_Y_on_edge(n)
            end do
            coefficient = sum_of_fluxes / sum_of_rho_Y
            do n=1,nspecies
               yflux_for_Y(i,j,n) 
     &              = yflux_for_Y(i,j,n) 
     &              - rho_Y_on_edge(n)
     &              * coefficient
            end do
         end do
      end do

c///
c///  form the updates weighted by scale_factor / vol
c///

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            div_of_flux_for_H(i, j) 
     &           = scale_factor_over_volume
     &           * ((xflux_for_H(i+1,j) - xflux_for_H(i,j)) 
     &           + (yflux_for_H(i,j+1) - yflux_for_H(i,j)))
            do n=1,nspecies
               div_of_flux_for_Y(i, j, n) 
     &              = scale_factor_over_volume
     &              * ((xflux_for_Y(i+1,j,n) - xflux_for_Y(i,j,n))
     &              + (yflux_for_Y(i,j+1,n) - yflux_for_Y(i,j,n)))
            end do
         end do
      end do

c///
c///  release temporary arrays (not strictly necessary)
c///

      deallocate (contiguous_Y)
      deallocate (cpms_on_edge)
      deallocate (derivative_of_Y)
      deallocate (hms_on_edge_cgs)
      deallocate (hms_on_edge_si)
      deallocate (one_over_wt)
      deallocate (residual)
      deallocate (rho_Dmix_on_edge_si)
      deallocate (rho_Upsilon_on_edge_cgs)
      deallocate (rho_Y)
      deallocate (rho_Y_on_edge)
      deallocate (T)
      deallocate (wt)
      deallocate (X)
      deallocate (X_on_edge)
      deallocate (Y)
      deallocate (Y_on_edge)

c///
c///  end of FORT_RK_MIXTURE_AVERAGED
c///
#endif 
      return
      end

c///
c///  start of FORT_RK_MULTICOMPONENT
c///

c     evaluate the "operator" (extensive fluxes and their divergences) for the 
c     Runge-Kutta implementation of the diffusion update using the multicomponent
c     formulation with Soret and Dufour effects

c     the divergences are scaled by 1/vol so they do approximate the analytic
c     divergence.  when used in the Runge-Kutta formula their additional scaling,
c     given by scale_factor, should be set to - dt.  the minus occurs because
c     the standard ODE is written y_prime = f(y), that is, the divergences are
c     moved to the opposite side of the equation from the time derivative.  when
c     used in getViscTerms the scale_factor should be set to -1.

c     to do:
c     1) must adjust diffusion coefficients by thickness factors
c     2) unsure of LMC support for r-z coordinates

c     arguments are alphabetical, mostly:
c
c     domain_lo, domain_hi,             ! INPUT limits of valid region of the domain
c     lo, hi,                           ! INPUT limits of valid region of the box
c   * areax, DIMS(areax),               ! INPUT areas of the faces perendicular to x axis
c   * areay, DIMS(areay),               ! INPUT areas of the faces perpendicular to y axis
c     bc,                               ! INPUT boundary condition array for all comps
c     dx,                               ! INPUT physical dimensions of grid cells
c     index_of_firstY,                  ! INPUT index of rho Y for the first species in the state
c     index_of_lastY,                   ! INPUT index of rho Y for the last species in the state
c     index_of_rho,                     ! INPUT index of rho in the state
c     index_of_rhoH,                    ! INPUT index of rho H in the state
c     index_of_T,                       ! INPUT index of T in the state
c     maximum_error,                    ! INPUT maximum error in calculation of T
c     maximum_iterations,               ! INPUT maximum iterations in calculation of T
c     ncomps,                           ! INPUT total number of components in the state
c     nspecies,                         ! INPUT total number of species in the state
c     scale_factor,                     ! INPUT scale_factor
c     state, DIMS(state),               ! INPUT all variables in the state
c   * volume, DIMS(volume),             ! INPUT volumes of the cells
c     div_of_flux_for_H, DIMS(div_of_flux_for_H), ! OUTPUT divergence of the flux for rho H
c     div_of_flux_for_Y, DIMS(div_of_flux_for_Y), ! OUTPUT divergences of the fluxes for rho Y
c     xflux_for_H, DIMS(xflux_for_H),   ! OUTPUT extensive x fluxes for rho H
c     xflux_for_Y, DIMS(xflux_for_Y),   ! OUTPUT extensive x fluxes for rho Y
c     yflux_for_H, DIMS(yflux_for_H),   ! OUTPUT extensive y fluxes for rho H
c     yflux_for_Y, DIMS(yflux_for_Y),   ! OUTPUT extensive y fluxes for rho Y
c
c     * these arguments are not used

      subroutine FORT_RK_MULTICOMPONENT
     &     (
     &     domain_lo, domain_hi,
     &     lo, hi,
     &     areax, DIMS(areax),
     &     areay, DIMS(areay),
     &     bc,
     &     dx,
     &     index_of_firstY,
     &     index_of_lastY,
     &     index_of_rho,
     &     index_of_rhoH,
     &     index_of_T,
     &     maximum_error,
     &     maximum_iterations,
     &     ncomps,
     &     nspecies,
     &     scale_factor,
     &     state, DIMS(state),
     &     volume, DIMS(volume),
     &     div_of_flux_for_H, DIMS(div_of_flux_for_H),
     &     div_of_flux_for_Y, DIMS(div_of_flux_for_Y),
     &     xflux_for_H, DIMS(xflux_for_H),
     &     xflux_for_Y, DIMS(xflux_for_Y),
     &     yflux_for_H, DIMS(yflux_for_H),
     &     yflux_for_Y, DIMS(yflux_for_Y)
     &     )

      implicit none

#include "BC_TYPES.H"
#include "cdwrk.H"

c     arguments
      integer ncomps
      integer nspecies
      integer domain_lo(BL_SPACEDIM), domain_hi(BL_SPACEDIM)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer DIMDEC(areax)
      integer DIMDEC(areay)
      integer bc(BL_SPACEDIM,2,ncomps)
      integer DIMDEC(state)
      integer DIMDEC(div_of_flux_for_H)
      integer DIMDEC(div_of_flux_for_Y)
      integer DIMDEC(xflux_for_H)
      integer DIMDEC(yflux_for_H)
      integer DIMDEC(xflux_for_Y)
      integer DIMDEC(yflux_for_Y)
      integer DIMDEC(volume)
      integer index_of_firstY
      integer index_of_lastY
      integer index_of_rho
      integer index_of_rhoH
      integer index_of_T
      integer maximum_iterations
      REAL_T areax(DIMV(areax))
      REAL_T areay(DIMV(areay))
      REAL_T div_of_flux_for_H(DIMV(div_of_flux_for_H))
      REAL_T div_of_flux_for_Y(DIMV(div_of_flux_for_Y),nspecies)
      REAL_T dx(BL_SPACEDIM)
      REAL_T maximum_error
      REAL_T scale_factor
      REAL_T state(DIMV(state),ncomps)
      REAL_T volume(DIMV(volume))
      REAL_T xflux_for_H(DIMV(xflux_for_H))
      REAL_T yflux_for_H(DIMV(yflux_for_H))
      REAL_T xflux_for_Y(DIMV(xflux_for_Y),ncomps)
      REAL_T yflux_for_Y(DIMV(yflux_for_Y),ncomps)

c     local variables
      external FORT_TfromHYpt 
      integer FORT_TfromHYpt
      integer ell
      integer i
      integer iterations
      integer j
      integer n
      integer ni
      integer nj
      integer offset
      logical at_hi
      logical at_lo
      logical error
      logical hi_dirichlet_check
      logical hi_dirichlet_for_T
      logical hi_dirichlet_for_Y
      logical hi_internal
      logical hi_internal_check
      logical lo_dirichlet_check
      logical lo_dirichlet_for_T
      logical lo_dirichlet_for_Y
      logical lo_internal
      logical lo_internal_check
      REAL_T area_for_x
      REAL_T area_for_y
      REAL_T ambient_pressure_cgs
      REAL_T ambient_pressure_si
      REAL_T derivative_of_T
      REAL_T dpdt_factor
      REAL_T floored_T
      REAL_T FORT_P1ATMMKS
      REAL_T H
      REAL_T hi_coefficient_0
      REAL_T hi_coefficient_1
      REAL_T hi_coefficient_2
      REAL_T lambda_prime_on_edge_cgs
      REAL_T lambda_prime_on_edge_si
      REAL_T lo_coefficient_0
      REAL_T lo_coefficient_1
      REAL_T lo_coefficient_2
      REAL_T minus_value
      REAL_T mmw_on_edge
      REAL_T one_atmosphere_cgs
      REAL_T one_atmosphere_si
      REAL_T one_over_dx
      REAL_T one_over_dy
      REAL_T one_over_rho
      REAL_T plus_value
      REAL_T pressure_over_rho_on_edge_cgs
      REAL_T pressure_over_rho_on_edge_si
      REAL_T RUC
      REAL_T scale_factor_over_volume
      REAL_T sum_h_xflux_for_Y
      REAL_T sum_h_yflux_for_Y
      REAL_T sum_rho_theta_derivative_of_X
      REAL_T sum_rho_Y_D_derivative_of_X
      REAL_T T_on_edge
      REAL_T universal_gas_constant_cgs

#ifdef BL_CYGWIN_NT
      print *,'Need to rebuild with f90'
      stop
#else
c     local arrays
      REAL_T, allocatable :: contiguous_X (:)
      REAL_T, allocatable :: contiguous_Y (:)
      REAL_T, allocatable :: cpms_on_edge (:)
      REAL_T, allocatable :: derivative_of_X (:)
      REAL_T, allocatable :: hms_on_edge_cgs (:)
      REAL_T, allocatable :: hms_on_edge_si (:)
      REAL_T, allocatable :: residual (:)
      REAL_T, allocatable :: rho_theta_on_edge_cgs (:)
      REAL_T, allocatable :: rho_theta_on_edge_si (:)
      REAL_T, allocatable :: rho_Y_D_on_edge_cgs (:,:)
      REAL_T, allocatable :: rho_Y_D_on_edge_si (:,:)
      REAL_T, allocatable :: T (:,:)
      REAL_T, allocatable :: X (:,:,:)
      REAL_T, allocatable :: X_on_edge (:)
      REAL_T, allocatable :: Y (:,:,:)
      REAL_T, allocatable :: Y_on_edge (:)

      allocate (contiguous_X (nspecies))
      allocate (contiguous_Y (nspecies))
      allocate (cpms_on_edge (nspecies))
      allocate (derivative_of_X (nspecies))
      allocate (hms_on_edge_cgs (nspecies))
      allocate (hms_on_edge_si (nspecies))
      allocate (residual (maximum_iterations))
      allocate (rho_theta_on_edge_cgs (nspecies))
      allocate (rho_theta_on_edge_si (nspecies))
      allocate (rho_Y_D_on_edge_cgs (nspecies, nspecies))
      allocate (rho_Y_D_on_edge_si (nspecies, nspecies))
      allocate (T (DIMV(state)))
      allocate (X (DIMV(state), nspecies))
      allocate (X_on_edge (nspecies))
      allocate (Y (DIMV(state), nspecies))
      allocate (Y_on_edge (nspecies))

c///
c///  initialization
c///

c     offset of species in the state
      offset = index_of_firstY - 1

c     some constants, at least in Cartesian coordinates
      area_for_x = dx(2)
      area_for_y = dx(1)
      one_over_dx = 1.0d0 / dx(1)
      one_over_dy = 1.0d0 / dx(2)
      scale_factor_over_volume = scale_factor / (dx(1) * dx(2))

c     coefficients for difference formulas at the hi and low edges
c     note still must divide by dx or dy to get the approximate derivative
      hi_coefficient_0 = + 1.0d0 / 3.0d0
      hi_coefficient_1 = - 3.0d0
      hi_coefficient_2 = + 8.0d0 / 3.0d0
      lo_coefficient_0 = - 8.0d0 / 3.0d0
      lo_coefficient_1 = + 3.0d0
      lo_coefficient_2 = - 1.0d0 / 3.0d0

c     get ambient pressure
c     one_atmosphere_cgs = 1013250.0 (dyne / cm**2) = 10.0 (Pascal)
c     one_atmosphere_si = 101325.0 (Pascal)
      call CKRP (IWRK(ckbi), RWRK(ckbr), universal_gas_constant_cgs, RUC, one_atmosphere_cgs)
      call FORT_GETPAMB (ambient_pressure_si, dpdt_factor)
      one_atmosphere_si = FORT_P1ATMMKS ()
      ambient_pressure_cgs = one_atmosphere_cgs * ambient_pressure_si / one_atmosphere_si

c     calculate T, X, and Y from rho, rho H and rho Y in the state
      do j=ARG_L2(state),ARG_H2(state)
         do i=ARG_L1(state),ARG_H1(state)
            one_over_rho = 1.0d0 / state(i, j, index_of_rho)
            H = one_over_rho * state(i, j, index_of_rhoH)
            do n=1,nspecies
               Y(i, j, n) = one_over_rho * state(i, j, offset + n)
            end do
            do n=1,nspecies
               contiguous_Y(n) = Y(i, j, n)
            end do
c           use the state value for T as a guess for the Newton iteration
            T(i, j) = state(i, j, index_of_T)
            iterations = FORT_TfromHYpt 
     &           (T(i, j),
     &           H, 
     &           contiguous_Y,
     &           maximum_error,
     &           maximum_iterations,
     &           residual)
            CALL CKYTX (contiguous_Y, IWRK(ckbi), RWRK(ckbr), contiguous_X)
            do n=1,nspecies
               X(i, j, n) = contiguous_X(n)
            end do
         end do
      end do

c///
c///  x fluxes
c///

c     set flags for special treatment of the lo boundary

      lo_dirichlet_for_T = .false.
      lo_dirichlet_for_Y = .false.
      lo_internal        = .true.
      if (domain_lo(1) .eq. lo(1)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         lo_dirichlet_for_T = EXT_DIR .eq. bc(1, 1, index_of_T)
         lo_dirichlet_for_Y = EXT_DIR .eq. bc(1, 1, index_of_firstY)
         lo_dirichlet_check = EXT_DIR .eq. bc(1, 1, index_of_firstY)
         lo_internal        = INT_DIR .eq. bc(1, 1, index_of_T)
         lo_internal_check  = INT_DIR .eq. bc(1, 1, index_of_T)
         do n=1,nspecies
            lo_dirichlet_for_Y = lo_dirichlet_for_Y .and. (EXT_DIR .eq. bc(1, 1, offset + n))
            lo_dirichlet_check = lo_dirichlet_check  .or. (EXT_DIR .eq. bc(1, 1, offset + n))
            lo_internal        = lo_internal        .and. (INT_DIR .eq. bc(1, 1, offset + n))
            lo_internal_check  = lo_internal_check   .or. (INT_DIR .eq. bc(1, 1, offset + n))
         end do
         error = .not. (
     &        (lo_dirichlet_for_Y .eqv. lo_dirichlet_check) .and.
     &        (lo_internal .eqv. lo_internal_check) .and.
     &        (.not. lo_dirichlet_for_Y .or. lo_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_MULTICOMPONENT: inconsistent bc array at low x")
      end if

c     set flags for special treatment of the hi boundary

      hi_dirichlet_for_T = .false.
      hi_dirichlet_for_Y = .false.
      hi_internal        = .true.
      if (domain_hi(1) .eq. hi(1)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         hi_dirichlet_for_T = EXT_DIR .eq. bc(1, 2, index_of_T)
         hi_dirichlet_for_Y = EXT_DIR .eq. bc(1, 2, index_of_firstY)
         hi_dirichlet_check = EXT_DIR .eq. bc(1, 2, index_of_firstY)
         hi_internal        = INT_DIR .eq. bc(1, 2, index_of_T)
         hi_internal_check  = INT_DIR .eq. bc(1, 2, index_of_T)
         do n=1,nspecies
            hi_dirichlet_for_Y = hi_dirichlet_for_Y .and. (EXT_DIR .eq. bc(1, 2, offset + n))
            hi_dirichlet_check = hi_dirichlet_check  .or. (EXT_DIR .eq. bc(1, 2, offset + n))
            hi_internal        = hi_internal        .and. (INT_DIR .eq. bc(1, 2, offset + n))
            hi_internal_check  = hi_internal_check   .or. (INT_DIR .eq. bc(1, 2, offset + n))
         end do
         error = .not. (
     &        (hi_dirichlet_for_Y .eqv. hi_dirichlet_check) .and.
     &        (hi_internal .eqv. hi_internal_check) .and.
     &        (.not. hi_dirichlet_for_Y .or. hi_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_MULTICOMPONENT: inconsistent bc array at high x")
      end if

c     calculate fluxes in the x direction

      do j=lo(2),hi(2)
c        fluxes in x direction are at edges indexed from lo(1) to hi(1)+1
         do i=lo(1),hi(1)+1

c           set flags indicating whether the edge is at the left or right boundary
            at_lo = i .eq. domain_lo(1)
            at_hi = i .eq. domain_hi(1) + 1

c           evaluate the value for T on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            if (at_lo .and. lo_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i-1, j, index_of_T)
            else if (at_hi .and. hi_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i, j, index_of_T)
            else
               T_on_edge = 0.5d0 * (T(i-1,j) + T(i,j))
            end if

c           use a floored T just in case
            T_on_edge = MAX (TMIN_TRANS, T_on_edge)

c           evaluate the values for Y on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            do n = 1,nspecies
               if (at_lo .and. lo_dirichlet_for_Y) then
                  Y_on_edge(n) = Y(i-1, j, n)
               else if (at_hi .and. hi_dirichlet_for_Y) then
                  Y_on_edge(n) = Y(i, j, n)
               else
                  Y_on_edge(n) = 0.5d0 * (Y(i-1,j,n) + Y(i,j,n))
               end if
            end do

c           evaluate the values for X on the edge
            call CKYTX 
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), X_on_edge)

c           evaluate the derivative for T on the edge
            if (at_lo .and. .not. lo_internal) then
               if (lo_dirichlet_for_T) then
                  derivative_of_T 
     &                 = one_over_dx
     &                 * (lo_coefficient_0 * T_on_edge
     &                 + lo_coefficient_1 * T(i,j)
     &                 + lo_coefficient_2 * T(i+1,j))
               else
                  derivative_of_T = 0.0d0
               end if
            else if (at_hi .and. .not. hi_internal) then
               if (hi_dirichlet_for_T) then
                  derivative_of_T 
     &                 = one_over_dx
     &                 * (hi_coefficient_0 * T(i-2,j)
     &                 + hi_coefficient_1 * T(i-1,j)
     &                 + hi_coefficient_2 * T_on_edge)
               else
                  derivative_of_T = 0.0d0
               end if
            else
               derivative_of_T = (T(i,j) - T(i-1,j)) * one_over_dx
            end if

c           evaluate the derivatives for X on the edge
            do n = 1,nspecies
               if (at_lo .and. .not. lo_internal) then
                  if (lo_dirichlet_for_Y) then
                     derivative_of_X(n)
     &                    = one_over_dx
     &                    * (lo_coefficient_0 * X_on_edge(n)
     &                    + lo_coefficient_1 * X(i,j,n)
     &                    + lo_coefficient_2 * X(i+1,j,n))
                  else
                     derivative_of_X(n) = 0.0d0
                  end if
               else if (at_hi .and. .not. hi_internal) then
                  if (hi_dirichlet_for_Y) then
                     derivative_of_X(n)
     &                    = one_over_dx
     &                    * (hi_coefficient_0 * X(i-2,j,n)
     &                    + hi_coefficient_1 * X(i-1,j,n)
     &                    + hi_coefficient_2 * X_on_edge(n))
                  else
                     derivative_of_X(n) = 0.0d0
                  end if
               else
                  derivative_of_X(n) = (X(i,j,n) - X(i-1,j,n)) * one_over_dx
               end if
            end do

c           evaluate diffusion coefficients and other things on the edge
            call CKCPMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), cpms_on_edge)
            call EGSPAR 
     &           (T_on_edge, X_on_edge, Y_on_edge, 
     &           cpms_on_edge, RWRK(egbr), IWRK(egbi))
            call CKMMWY 
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), mmw_on_edge)
c           SUBROUTINE EGSLTDR5 (T, Y, WW, WEG, IWEG, PTC, THETA, D)
            call EGSLTDR5
     &           (T_on_edge, Y_on_edge, mmw_on_edge, RWRK(egbr), IWRK(egbi),
     &           lambda_prime_on_edge_cgs, rho_theta_on_edge_cgs, rho_Y_D_on_edge_cgs)
            call CKHMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms_on_edge_cgs)
c           there does not seem to be a chemkin call to get this ratio from T and Y
            pressure_over_rho_on_edge_cgs 
     &           = universal_gas_constant_cgs * T_on_edge / mmw_on_edge

c           convert the library values from cgs units to si units

c           per page 98 of the chemkin manual, the enthalpies for the species,
c           hms_on_edge_cgs, have units
c
c                            erg / g
c                 =          (J / 10**7) / (Kg / 10**3)
c                 = 1.0d-4 * J / Kg
c
c           the si unit of thermal conductivity is J / Kg, so the EGLIB values must
c           be multiplied by 1.0d-4 to convert them to si units.
            do n=1,nspecies
               hms_on_edge_si(n) = 1.0d-4 * hms_on_edge_cgs(n)
            end do

c           per page 23 of the EGLIB manual, the partial thermal conductivity, 
c           lambda_prime_on_edge_cgs, has units
c
c                            erg / (cm s K)
c                 =          (J / 10**7) / ((m / 10**2) s K)
c                 = 1.0d-5 * ((W s) / (m s K))
c                 = 1.0d-5 * W / (m K).
c
c           the si unit of thermal conductivity is W / (m K), so the EGLIB value must
c           be multiplied by 1.0d-5 to convert it to si units.
            lambda_prime_on_edge_si = 1.0d-5 * lambda_prime_on_edge_cgs

c           per page 23 in the EGLIB manual, the rescaled thermal diffusion vector,
c           rho_theta_on_edge_cgs, has units
c
c                            g / (cm s)
c                 =          (Kg / 10**3) / ((m / 10**2) s)
c                 = 1.0d-1 * Kg / (m s)
c
c           so the EGLIB values must be multiplied by 1.0d-1 to convert them to si units.
            do n=1,nspecies
               rho_theta_on_edge_si(n) = 1.0d-1 * rho_theta_on_edge_cgs(n)
            end do

c           per page 23 in the EGLIB manual, the rescaled flux diffusion matrix,
c           rho_Y_D_on_edge_cgs, has units
c
c                            g / (cm s)
c                 =          (Kg / 10**3) / ((m / 10**2) s)
c                 = 1.0d-1 * Kg / (m s)
c
c           so the EGLIB values must be multiplied by 1.0d-1 to convert them to si units.
            do ni=1,nspecies
               do nj=1,nspecies
                  rho_Y_D_on_edge_si(ni,nj) = 1.0d-1 * rho_Y_D_on_edge_cgs(ni,nj)
               end do
            end do

c           the cgs unit for pressure is the barye which equals 0.1 Pascal.
c           the si unit for pressure is the Pascal which equals 1 Kg m s**2.

c           the pressure and density ratio, pressure_over_rho_on_edge_cgs, has units
c
c                            barye / (g / cm**3)
c                 =          (1.0d-1 Pa) / ((1.0d-3 Kg) / (1.0d-2 m)**3)
c                 = 1.0d-4 * Pa / (Kg / m**2)
c
c           so the cgs value must be multiplied by 1.0d-4 to convert it to si units.
            pressure_over_rho_on_edge_si = 1.0d-4 * pressure_over_rho_on_edge_cgs 

c           assemble the fluxes
            sum_h_xflux_for_Y = 0
            sum_rho_theta_derivative_of_X = 0
            do n = 1,nspecies
               sum_rho_Y_D_derivative_of_X = 0
               do ell = 1,nspecies
                  sum_rho_Y_D_derivative_of_X 
     &                 = sum_rho_Y_D_derivative_of_X 
     &                 + rho_Y_D_on_edge_si(n,ell) 
     &                 * derivative_of_X(ell)
               end do
               xflux_for_Y(i,j,n) 
     &              = - area_for_x
     &              * (sum_rho_Y_D_derivative_of_X
     &              + Y_on_edge(n) * rho_theta_on_edge_si(n) * derivative_of_T / T_on_edge)
               sum_h_xflux_for_Y 
     &              = sum_h_xflux_for_Y 
     &              + hms_on_edge_si(n) 
     &              * xflux_for_Y(i,j,n)
               sum_rho_theta_derivative_of_X 
     &              = sum_rho_theta_derivative_of_X 
     &              + rho_theta_on_edge_si(n)
     &              * derivative_of_X(n)
            end do
            xflux_for_H(i,j)
     &           = sum_h_xflux_for_Y
     &           - area_for_x
     &           * (lambda_prime_on_edge_si * derivative_of_T
     &           + pressure_over_rho_on_edge_si * sum_rho_theta_derivative_of_X)
            end do
         end do

c///
c///  y fluxes
c///

c     set flags for special treatment of the lo boundary

      lo_dirichlet_for_T = .false.
      lo_dirichlet_for_Y = .false.
      lo_internal        = .true.
      if (domain_lo(2) .eq. lo(2)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         lo_dirichlet_for_T = EXT_DIR .eq. bc(2, 1, index_of_T)
         lo_dirichlet_for_Y = EXT_DIR .eq. bc(2, 1, index_of_firstY)
         lo_dirichlet_check = EXT_DIR .eq. bc(2, 1, index_of_firstY)
         lo_internal        = INT_DIR .eq. bc(2, 1, index_of_T)
         lo_internal_check  = INT_DIR .eq. bc(2, 1, index_of_T)
         do n=1,nspecies
            lo_dirichlet_for_Y = lo_dirichlet_for_Y .and. (EXT_DIR .eq. bc(2, 1, offset + n))
            lo_dirichlet_check = lo_dirichlet_check  .or. (EXT_DIR .eq. bc(2, 1, offset + n))
            lo_internal        = lo_internal        .and. (INT_DIR .eq. bc(2, 1, offset + n))
            lo_internal_check  = lo_internal_check   .or. (INT_DIR .eq. bc(2, 1, offset + n))
         end do
         error = .not. (
     &        (lo_dirichlet_for_Y .eqv. lo_dirichlet_check) .and.
     &        (lo_internal .eqv. lo_internal_check) .and.
     &        (.not. lo_dirichlet_for_Y .or. lo_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_MULTICOMPONENT: inconsistent bc array at low y")
      end if

c     set flags for special treatment of the hi boundary

      hi_dirichlet_for_T = .false.
      hi_dirichlet_for_Y = .false.
      hi_internal        = .true.
      if (domain_hi(2) .eq. hi(2)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         hi_dirichlet_for_T = EXT_DIR .eq. bc(2, 2, index_of_T)
         hi_dirichlet_for_Y = EXT_DIR .eq. bc(2, 2, index_of_firstY)
         hi_dirichlet_check = EXT_DIR .eq. bc(2, 2, index_of_firstY)
         hi_internal        = INT_DIR .eq. bc(2, 2, index_of_T)
         hi_internal_check  = INT_DIR .eq. bc(2, 2, index_of_T)
         do n=1,nspecies
            hi_dirichlet_for_Y = hi_dirichlet_for_Y .and. (EXT_DIR .eq. bc(2, 2, offset + n))
            hi_dirichlet_check = hi_dirichlet_check  .or. (EXT_DIR .eq. bc(2, 2, offset + n))
            hi_internal        = hi_internal        .and. (INT_DIR .eq. bc(2, 2, offset + n))
            hi_internal_check  = hi_internal_check   .or. (INT_DIR .eq. bc(2, 2, offset + n))
         end do
         error = .not. (
     &        (hi_dirichlet_for_Y .eqv. hi_dirichlet_check) .and.
     &        (hi_internal .eqv. hi_internal_check) .and.
     &        (.not. hi_dirichlet_for_Y .or. hi_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_MULTICOMPONENT: inconsistent bc array at high y")
      end if

c     calculate fluxes in the y direction

c     fluxes in y direction are at edges indexed from lo(2) to hi(2)+1
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)

c           set flags indicating whether the edge is at the left or right boundary
            at_lo = j .eq. domain_lo(2)
            at_hi = j .eq. domain_hi(2) + 1

c           evaluate the value for T on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            if (at_lo .and. lo_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i, j-1, index_of_T)
            else if (at_hi .and. hi_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i, j, index_of_T)
            else
               T_on_edge = 0.5d0 * (T(i,j-1) + T(i,j))
            end if

c           use a floored T just in case
            T_on_edge = MAX (TMIN_TRANS, T_on_edge)

c           evaluate the values for Y on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            do n = 1,nspecies
               if (at_lo .and. lo_dirichlet_for_Y) then
                  Y_on_edge(n) = Y(i, j-1, n)
               else if (at_hi .and. hi_dirichlet_for_Y) then
                  Y_on_edge(n) = Y(i, j, n)
               else
                  Y_on_edge(n) = 0.5d0 * (Y(i,j-1,n) + Y(i,j,n))
               end if
            end do

c           evaluate the values for X on the edge
            call CKYTX
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), X_on_edge)

c           evaluate the derivative for T on the edge
            if (at_lo .and. .not. lo_internal) then
               if (lo_dirichlet_for_T) then
                  derivative_of_T
     &                 = one_over_dy
     &                 * (lo_coefficient_0 * T_on_edge
     &                 + lo_coefficient_1 * T(i,j)
     &                 + lo_coefficient_2 * T(i,j+1))
               else
                  derivative_of_T = 0
               end if
            else if (at_hi .and. .not. hi_internal) then
               if (hi_dirichlet_for_T) then
                  derivative_of_T
     &                 = one_over_dy
     &                 * (hi_coefficient_0 * T(i,j-2)
     &                 + hi_coefficient_1 * T(i,j-1)
     &                 + hi_coefficient_2 * T_on_edge)
               else
                  derivative_of_T = 0
               end if
            else
               derivative_of_T = (T(i,j) - T(i,j-1)) * one_over_dy
            end if

c           evaluate the derivatives for X on the edge
            do n = 1,nspecies
               if (at_lo .and. .not. lo_internal) then
                  if (lo_dirichlet_for_Y) then
                     derivative_of_X(n)
     &                    = one_over_dy
     &                    * (lo_coefficient_0 * X_on_edge(n)
     &                    + lo_coefficient_1 * X(i,j,n)
     &                    + lo_coefficient_2 * X(i,j+1,n))
                  else
                     derivative_of_X(n) = 0.0d0
                  end if                     
               else if (at_hi .and. .not. hi_internal) then
                  if (hi_dirichlet_for_Y) then
                     derivative_of_X(n)
     &                    = one_over_dy
     &                    * (hi_coefficient_0 * X(i,j-2,n)
     &                    + hi_coefficient_1 * X(i,j-1,n)
     &                    + hi_coefficient_2 * X_on_edge(n))
                  else
                     derivative_of_X(n) = 0.0d0
                  end if
               else
                  derivative_of_X(n) = (X(i,j,n) - X(i,j-1,n)) * one_over_dy
               end if
            end do

c           evaluate diffusion coefficients and other things on the edge
            call CKCPMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), cpms_on_edge)
            call EGSPAR 
     &           (T_on_edge, X_on_edge, Y_on_edge, 
     &           cpms_on_edge, RWRK(egbr), IWRK(egbi))
            call CKMMWY 
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), mmw_on_edge)
c           SUBROUTINE EGSLTDR5 (T, Y, WW, WEG, IWEG, PTC, THETA, D)
            call EGSLTDR5 
     &           (T_on_edge, Y_on_edge, mmw_on_edge, RWRK(egbr), IWRK(egbi),
     &           lambda_prime_on_edge_cgs, rho_theta_on_edge_cgs, rho_Y_D_on_edge_cgs)
            call CKHMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms_on_edge_cgs)
c           there does not seem to be a chemkin call to get this ratio from T and Y
            pressure_over_rho_on_edge_cgs 
     &           = universal_gas_constant_cgs * T_on_edge / mmw_on_edge

c           convert the library values from cgs units to si units

c           per page 98 of the chemkin manual, the enthalpies for the species,
c           hms_on_edge_cgs, have units
c
c                            erg / g
c                 =          (J / 10**7) / (Kg / 10**3)
c                 = 1.0d-4 * J / Kg
c
c           the si unit of thermal conductivity is J / Kg, so the EGLIB values must
c           be multiplied by 1.0d-4 to convert them to si units.
            do n=1,nspecies
               hms_on_edge_si(n) = 1.0d-4 * hms_on_edge_cgs(n)
            end do

c           per page 23 of the EGLIB manual, the partial thermal conductivity, 
c           lambda_prime_on_edge_cgs, has units
c
c                            erg / (cm s K)
c                 =          (J / 10**7) / ((m / 10**2) s K)
c                 = 1.0d-5 * ((W s) / (m s K))
c                 = 1.0d-5 * W / (m K).
c
c           the si unit of thermal conductivity is W / (m K), so the EGLIB value must
c           be multiplied by 1.0d-5 to convert it to si units.
            lambda_prime_on_edge_si = 1.0d-5 * lambda_prime_on_edge_cgs

c           per page 23 in the EGLIB manual, the rescaled thermal diffusion vector,
c           rho_theta_on_edge_cgs, has units
c
c                            g / (cm s)
c                 =          (Kg / 10**3) / ((m / 10**2) s)
c                 = 1.0d-1 * Kg / (m s)
c
c           so the EGLIB values must be multiplied by 1.0d-1 to convert them to si units.
            do n=1,nspecies
               rho_theta_on_edge_si(n) = 1.0d-1 * rho_theta_on_edge_cgs(n)
            end do

c           per page 23 in the EGLIB manual, the rescaled flux diffusion matrix,
c           rho_Y_D_on_edge_cgs, has units
c
c                            g / (cm s)
c                 =          (Kg / 10**3) / ((m / 10**2) s)
c                 = 1.0d-1 * Kg / (m s)
c
c           so the EGLIB values must be multiplied by 1.0d-1 to convert them to si units.
            do ni=1,nspecies
               do nj=1,nspecies
                  rho_Y_D_on_edge_si(ni,nj) = 1.0d-1 * rho_Y_D_on_edge_cgs(ni,nj)
               end do
            end do

c           the cgs unit for pressure is the barye which equals 0.1 Pascal.
c           the si unit for pressure is the Pascal which equals 1 Kg m s**2.

c           the pressure and density ratio, pressure_over_rho_on_edge_cgs, has units
c
c                            barye / (g / cm**3)
c                 =          (1.0d-1 Pa) / ((1.0d-3 Kg) / (1.0d-2 m)**3)
c                 = 1.0d-4 * Pa / (Kg / m**2)
c
c           so the cgs value must be multiplied by 1.0d-4 to convert it to si units.
            pressure_over_rho_on_edge_si = 1.0d-4 * pressure_over_rho_on_edge_cgs 

c           assemble the fluxes
            sum_h_yflux_for_Y = 0
            sum_rho_theta_derivative_of_X = 0
            do n = 1,nspecies
               sum_rho_Y_D_derivative_of_X = 0
               do ell = 1,nspecies
                  sum_rho_Y_D_derivative_of_X 
     &                 = sum_rho_Y_D_derivative_of_X 
     &                 + rho_Y_D_on_edge_si(n,ell) 
     &                 * derivative_of_X(ell)
               end do
               yflux_for_Y(i,j,n) 
     &              = - area_for_y
     &              * (sum_rho_Y_D_derivative_of_X
     &              + Y_on_edge(n) * rho_theta_on_edge_si(n) * derivative_of_T / T_on_edge)
               sum_h_yflux_for_Y 
     &              = sum_h_yflux_for_Y 
     &              + hms_on_edge_si(n) 
     &              * yflux_for_Y(i,j,n)
               sum_rho_theta_derivative_of_X 
     &              = sum_rho_theta_derivative_of_X 
     &              + rho_theta_on_edge_si(n) 
     &              * derivative_of_X(n)
            end do
            yflux_for_H(i,j)
     &           = sum_h_yflux_for_Y
     &           - area_for_y
     &           * (lambda_prime_on_edge_si * derivative_of_T
     &           + pressure_over_rho_on_edge_si * sum_rho_theta_derivative_of_X)
         end do
      end do

c///
c///  form the updates weighted by scale_factor / vol
c///

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            div_of_flux_for_H(i, j) 
     &           = scale_factor_over_volume
     &           * ((xflux_for_H(i+1,j) - xflux_for_H(i,j)) 
     &           + (yflux_for_H(i,j+1) - yflux_for_H(i,j)))
            do n=1,nspecies
               div_of_flux_for_Y(i, j, n) 
     &              = scale_factor_over_volume
     &              * ((xflux_for_Y(i+1,j,n) - xflux_for_Y(i,j,n))
     &              + (yflux_for_Y(i,j+1,n) - yflux_for_Y(i,j,n)))
            end do
         end do
      end do

c///
c///  release temporary arrays (not strictly necessary)
c///

      deallocate (contiguous_X)
      deallocate (contiguous_Y)
      deallocate (cpms_on_edge)
      deallocate (derivative_of_X)
      deallocate (hms_on_edge_cgs)
      deallocate (hms_on_edge_si)
      deallocate (residual)
      deallocate (rho_theta_on_edge_cgs)
      deallocate (rho_theta_on_edge_si)
      deallocate (rho_Y_D_on_edge_cgs)
      deallocate (rho_Y_D_on_edge_si)
      deallocate (T)
      deallocate (X)
      deallocate (X_on_edge)
      deallocate (Y)
      deallocate (Y_on_edge)

c///
c///  end of FORT_RK_MULTICOMPONENT
c///
#endif
      return
      end

c///
c///  start of FORT_RK_STEP_SELECTION
c///

c     estimate the time step for the Runge-Kutta implementation of the diffusion
c     update.  this estimation is based on the diffusion coefficients of the 
c     species into the mixture, independent of whether or not the "operator" is 
c     based on the mixture averaged or the multicomponent formulation with Soret 
c     and Dufour effects.
c
c                               space_step ** 2
c                 time_step  =  ---------------
c                               D_mix * 2 * DIM
c
c     the minimum of this formula is taken over all cells and species

c     to do:
c     1) must adjust diffusion coefficients by thickness factors
c     2) unsure of LMC support for r-z coordinates
c     3) could include an estimate for thermal diffusion

c     arguments are alphabetical, mostly:
c
c     lo, hi             ! INPUT limits of valid region of the box
c     dx                 ! INPUT physical dimensions of grid cells
c     index_of_firstY    ! INPUT index of rho Y for the first species in the state
c     index_of_lastY     ! INPUT index of rho Y for the last species in the state
c     index_of_rho       ! INPUT index of rho in the state
c     index_of_rhoH      ! INPUT index of rho H in the state
c     index_of_T         ! INPUT index of T in the state
c     maximum_error      ! INPUT maximum error in calculation of T
c     maximum_iterations ! INPUT maximum iterations in calculation of T
c     ncomps             ! INPUT total number of components in the state
c     nspecies           ! INPUT total number of species in the state
c     smallest_time_step ! OUTPUT smallest time step over all cells and species
c     state, DIMS(state) ! INPUT all variables in the state

      subroutine FORT_RK_STEP_SELECTION
     &     (
     &     lo, hi,
     &     dx,
     &     index_of_firstY,
     &     index_of_lastY,
     &     index_of_rho,
     &     index_of_rhoH,
     &     index_of_T,
     &     maximum_error,
     &     maximum_iterations,
     &     ncomps,
     &     nspecies,
     &     smallest_time_step,
     &     state, DIMS(state)
     &     )

      implicit none

#include "REAL.H"
#include "HEATTRANSFER_F.H"
#include "ChemDriver_F.H"
#include "cdwrk.H"

c     arguments
      integer ncomps
      integer nspecies
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer DIMDEC(state)
      integer index_of_firstY
      integer index_of_lastY
      integer index_of_rho
      integer index_of_rhoH
      integer index_of_T
      integer maximum_iterations
      REAL_T maximum_error
      REAL_T dx(BL_SPACEDIM)
      REAL_T smallest_time_step
      REAL_T state(DIMV(state),ncomps)

c     local variables
      external FORT_TfromHYpt 
      integer FORT_TfromHYpt 
      integer i
      integer iterations
      integer j
      integer n
      integer offset
      logical first
      REAL_T ambient_pressure_cgs
      REAL_T ambient_pressure_si
      REAL_T D
      REAL_T dpdt_factor
      REAL_T floored_T
      REAL_T FORT_P1ATMMKS
      REAL_T H
      REAL_T mmw
      REAL_T one_atmosphere_cgs
      REAL_T one_atmosphere_si
      REAL_T one_over_rho
      REAL_T rho_cgs
      REAL_T RU
      REAL_T RUC
      REAL_T space_step
      REAL_T T
      REAL_T time_step

#ifdef BL_CYGWIN_NT
      print *,'Need to rebuild with f90'
      stop
#else
c     local arrays
      REAL_T, allocatable :: contiguous_X (:)
      REAL_T, allocatable :: contiguous_Y (:)
      REAL_T, allocatable :: cpms (:)
      REAL_T, allocatable :: one_over_wt (:)
      REAL_T, allocatable :: residual (:)
      REAL_T, allocatable :: Upsilon (:)
      REAL_T, allocatable :: wt (:)

      allocate (contiguous_X (nspecies))
      allocate (contiguous_Y (nspecies))
      allocate (cpms (nspecies))
      allocate (one_over_wt (nspecies))
      allocate (residual (maximum_iterations))
      allocate (Upsilon (nspecies))
      allocate (wt (nspecies))

c///
c///  initialization
c///

c     offset of species in the state
      offset = index_of_firstY - 1

c     get ambient pressure in cgs units
c     one_atmosphere_cgs = 1013250.0 (dyne / cm**2)
c     one_atmosphere_si = 101325.0 (pascal)
      call CKRP (IWRK(ckbi), RWRK(ckbr), RU, RUC, one_atmosphere_cgs)
      call FORT_GETPAMB (ambient_pressure_si, dpdt_factor)
      one_atmosphere_si = FORT_P1ATMMKS ()
      ambient_pressure_cgs = one_atmosphere_cgs * ambient_pressure_si / one_atmosphere_si

c     reciprocal of molecular weights
      call CKWT (IWRK(ckbi),RWRK(ckbr),wt)
      do n=1,nspecies
         one_over_wt(n) = 1.0d0 / wt(n)
      end do

c     smallest space step
      space_step = min (dx(1), dx(2))

c///
c///  find the smallest time step for the given part of the domain
c///

      first = .true.
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

c           calculate T and Y from rho, rho H, and rho Y in the state
            one_over_rho = 1.0d0 / state(i, j, index_of_rho)
            H = one_over_rho * state(i, j, index_of_rhoH)
            do n=1,nspecies
               contiguous_Y(n) = one_over_rho * state(i, j, offset + n)
            end do
c           use the state value for T as a guess for the Newton iteration
            T = state(i, j, index_of_T)
            iterations = FORT_TfromHYpt 
     &           (T, 
     &           H, 
     &           contiguous_Y,
     &           maximum_error,
     &           maximum_iterations,
     &           residual)
            floored_T = MAX (T, TMIN_TRANS) 

c           evaluate the diffusion coefficients
            CALL CKCPMS 
     &           (floored_T, IWRK(ckbi), RWRK(ckbr), cpms)
            CALL CKYTX 
     &           (contiguous_Y, IWRK(ckbi), RWRK(ckbr), contiguous_X)
            CALL EGSPAR 
     &           (floored_T, contiguous_X, contiguous_Y, cpms, RWRK(egbr), 
     &           IWRK(egbi))
            CALL CKMMWY 
     &           (contiguous_Y, IWRK(ckbi), RWRK(ckbr), mmw)
            CALL EGSV1 
     &           (ambient_pressure_cgs, floored_T, contiguous_Y, mmw, 
     &           RWRK(egbr), Upsilon)
c           get a rho_cgs consistent with floored_T and the rest
            CALL CKRHOY 
     &           (ambient_pressure_cgs, floored_T, contiguous_Y, IWRK(ckbi),
     &           RWRK(ckbr), rho_cgs)

c           select the smallest time step
            do n=1,nspecies
c              per page 17 in the EGLIB manual, Upsilon has units (cm**2 / s).  
c              per equation 12 in the EGLIB manual, the mixture averaged diffusion 
c              coefficient is Upsilon times the dimensionless mmw / wt(n).  thus 
c              the resulting number with units (cm**2 / s) = 1.0d-4 (m**2 / s) must
c              be multiplied by 1.0d-4 to obtain si units.
               D = mmw * one_over_wt(n) * Upsilon(n) * 1.0d-4
               time_step = space_step ** 2 / (D * 2.0d0 * BL_SPACEDIM)
               if (first) then
                  smallest_time_step = time_step
                  first = .false.
               else
                  smallest_time_step = min (time_step, smallest_time_step)
               end if
            end do

         end do
      end do

c///
c///  release temporary arrays (not strictly necessary)
c///

      deallocate (contiguous_X)
      deallocate (contiguous_Y)
      deallocate (cpms)
      deallocate (one_over_wt)
      deallocate (residual)
      deallocate (Upsilon)
      deallocate (wt)

c///
c///  end of FORT_RK_STEP_SELECTION
c///
#endif
      return
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
