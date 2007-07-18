c
c $Id: HEATTRANSFER_2D.F,v 1.4 2007-07-18 20:19:38 sepp Exp $
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
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     do n=1,ncompd-1
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
                  val = rhoD(i,j,1) * Yfac
                  do n=1,nspecies
                     rhoD(i,j,n) = val
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
         do j=lo(2), hi(2)
            do i=lo(1), hi(1)
               do n=1,nspecies
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

#include "FORT_RK_MIXTURE_AVERAGED.F"

#include "FORT_RK_MULTICOMPONENT.F"

#include "FORT_RK_STEP_SELECTION.F"
