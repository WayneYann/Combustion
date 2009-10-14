c
c $Id: HEATTRANSFER_3D.F,v 1.14 2009-10-14 23:24:08 lijewski Exp $
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

#define SDIM 3

      subroutine FORT_COMPUTE_RHODGRADHDOTGRADY(dx,
     &         lo, hi, DIMS(species), species,
     &         DIMS(h), h, DIMS(betax), betax,
     &         DIMS(betay), betay, DIMS(betaz), betaz, 
     &         DIMS(rdghdgy), rdghdgy)

      implicit none

      integer lo(SDIM), hi(SDIM)
      REAL_T  dx(SDIM)
      integer DIMDEC(species)
      integer DIMDEC(h)
      REAL_T  species(DIMV(species))
      REAL_T  h(DIMV(h))
      integer DIMDEC(betax)
      integer DIMDEC(betay)
      integer DIMDEC(betaz)
      REAL_T betax(DIMV(betax))
      REAL_T betay(DIMV(betay))
      REAL_T betaz(DIMV(betaz))
      integer DIMDEC(rdghdgy)

      REAL_T rdghdgy(DIMV(rdghdgy))

      integer i,j,k
      REAL_T  betaxcen, betaycen, betazcen, dxsqr, dysqr, dzsqr
      REAL_T  betadotxlo, betadotxhi, betadotylo, betadotyhi
      REAL_T  betadotzlo, betadotzhi

      dxsqr = dx(1)**2
      dysqr = dx(2)**2
      dzsqr = dx(3)**2
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
#if 0
            betaxcen = half*(betax(i,j,k)+betax(i+1,j,k))
            betaycen = half*(betay(i,j,k)+betay(i,j+1,k))
            betazcen = half*(betaz(i,j,k)+betaz(i,j,k+1))
            rdghdgy(i,j,k) = betaxcen*(h(i+1,j,k)-h(i-1,j,k))*
     &                         (species(i+1,j,k)-species(i-1,j,k))/dxsqr +
     &                   betaycen*(h(i,j+1,k)-h(i,j-1,k))*
     &                         (species(i,j+1,k)-species(i,j-1,k))/dysqr +
     &                   betazcen*(h(i,j,k+1)-h(i,j,k-1))*
     &                          (species(i,j,k+1)-species(i,j,k-1))/dzsqr
#else 
            betadotxlo = betax(i,j,k  )*(h(i,j,k)-h(i-1,j,k))
     &                        *(species(i,j,k)-species(i-1,j,k))
            betadotxhi = betax(i+1,j,k)*(h(i+1,j,k)-h(i,j,k))
     &                        *(species(i+1,j,k)-species(i,j,k))
            betadotylo = betay(i,j,k  )*(h(i,j,k)-h(i,j-1,k))
     &                        *(species(i,j,k)-species(i,j-1,k))
            betadotyhi = betay(i,j+1,k)*(h(i,j+1,k)-h(i,j,k))
     &                        *(species(i,j+1,k)-species(i,j,k))
            betadotzlo = betaz(i,j,k  )*(h(i,j,k)-h(i,j,k-1))
     &                        *(species(i,j,k)-species(i,j,k-1))
            betadotzhi = betaz(i,j,k+1)*(h(i,j,k+1)-h(i,j,k))
     &                        *(species(i,j,k+1)-species(i,j,k))
            rdghdgy(i,j,k) =  half*((betadotxlo + betadotxhi)/dxsqr +
     &                            (betadotylo + betadotyhi)/dysqr +
     &                            (betadotzlo + betadotzhi)/dysqr)

#endif

          enddo
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

      integer i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               sum(i,j,k) = zero
            enddo
         enddo
      enddo
      do n=1,nspecies
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  spec(i,j,k,n) = max(zero,spec(i,j,k,n))
                  sum(i,j,k)    = sum(i,j,k)+spec(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
      do n=1,nspecies
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  spec(i,j,k,n) = spec(i,j,k,n)/sum(i,j,k)
               enddo
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
      
      integer i, j, k, n

      if (.not.use_constant_mu) then
         if (.not. LeEQ1) then
            call FORT_MIX_SHEAR_VISC(lo, hi, mu, DIMS(mu),
     &                               T, DIMS(T), Y, DIMS(Y))       
         else 
            do k=lo(3), hi(3)
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     mu(i,j,k) = 1.85e-5*(T(i,j,k)/298.0)**.7
                  end do
               end do
            end do
         end if
      else
         do k=lo(3), hi(3)
            do j=lo(2), hi(2)
               do i=lo(1), hi(1)
                  mu(i,j,k) = constant_mu_val
               end do
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
      
      integer i, j, k, n
      REAL_T Yt(maxspec)
      REAL_T cpmix
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/
      
      print *,'**** in FORT_TEMPVISC'
      call bl_abort()

      if (.not.use_constant_lambda) then
         if (.not. LeEQ1) then
            call FORT_MIX_COND(lo, hi, lambda, DIMS(lambda),
     &           T, DIMS(T), Y, DIMS(Y))
         else
            call FORT_VELVISC(lo,hi,DIMS(T),T,DIMS(Y),Y,
     &           DIMS(lambda),lambda)
            do k=lo(3), hi(3)
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     do n=1,Nspec
                        Yt(n) = Y(i,j,k,n)
                     end do
                     CALL FORT_CPMIXfromTY(lo_chem, hi_chem, 
     &                    CPMIX,    ARLIM(lo_chem), ARLIM(hi_chem),
     &                    T(i,j,k), ARLIM(lo_chem), ARLIM(hi_chem),
     &                    Yt,       ARLIM(lo_chem), ARLIM(hi_chem))
                     lambda(i,j,k) = lambda(i,j,k) * cpmix / Pr
                  end do
               end do
            end do
         end if
      else
         do k=lo(3), hi(3)
            do j=lo(2), hi(2)
               do i=lo(1), hi(1)
                  lambda(i,j,k) = constant_lambda_val
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
 
      integer i, j, k, n, nspecies
      REAL_T P1ATM_MKS, Ptmp, val, Patm
      REAL_T Yt(maxspec)
      REAL_T cpmix, Tfac, Yfac
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      nspecies = ncompd - 2

c     ncompd-1 is Temp stuff
c     ncompd   is velocity stuff 

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
                  do k=lo(3), hi(3)
                     do j=lo(2), hi(2)
                        do i=lo(1), hi(1)
                           rhoD(i,j,k,n) = rhoD(i,j,k,n)*thickFacTR
                        end do
                     end do
                  end do
               end do
            endif
         else
            call FORT_VELVISC(lo,hi,DIMS(T),T,DIMS(Y),Y,
     &           DIMS(rhoD),rhoD)
            do k=lo(3), hi(3)
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     rhoD(i,j,k,ncompd-1) = rhoD(i,j,k,1)
                     rhoD(i,j,k,ncompd  ) = rhoD(i,j,k,1)
                     rhoD(i,j,k,1)        = rhoD(i,j,k,1) * Yfac
                  end do
               end do
            end do
            do n=2,nspecies
               do k=lo(3), hi(3)
                  do j=lo(2), hi(2)
                     do i=lo(1), hi(1)
                        rhoD(i,j,k,n) = rhoD(i,j,k,1)
                     end do
                  end do
               end do
            end do

            if (do_temp .ne. 0) then
               do k=lo(3), hi(3)
                  do j=lo(2), hi(2)
                     do i=lo(1), hi(1)
                        do n=1,Nspec
                           Yt(n) = Y(i,j,k,n)
                        end do
                        CALL FORT_CPMIXfromTY(lo_chem, hi_chem, 
     &                       CPMIX,    ARLIM(lo_chem), ARLIM(hi_chem),
     &                       T(i,j,k), ARLIM(lo_chem), ARLIM(hi_chem),
     &                       Yt,       ARLIM(lo_chem), ARLIM(hi_chem))
                        rhoD(i,j,k,ncompd-1) = rhoD(i,j,k,ncompd-1)*cpmix*Tfac
                     end do
                  end do
               end do
            endif
         end if
      else
         do n=1,nspecies
            do k=lo(3), hi(3)
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     rhoD(i,j,k,n) = constant_rhoD_val
                  end do
               end do
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
      
      integer i, j, k
      integer n
      REAL_T P1ATM_MKS, Ptmp, val, Patm

      print *,'**** in FORT_SPECVISC'
      call bl_abort()

      end

      subroutine FORT_EST_DIVU_DT(flag, dtfactor, delta, divu, DIMS(divu), 
     &                           dsdt, rho, DIMS(rho), 
     &                           u, DIMS(u), 
     &                           volume, DIMS(volume), 
     &                           areax,  DIMS(areax),
     &                           areay,  DIMS(areay),
     &                           areaz,  DIMS(areaz),
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
      integer DIMDEC(areaz)
      REAL_T  volume(DIMV(volume))
      REAL_T  areax(DIMV(areax))
      REAL_T  areay(DIMV(areay))
      REAL_T  areaz(DIMV(areaz))

      REAL_T dt

      integer i,j,k
      REAL_T  dtcell, dtcell2, denom, rhominij, rhoij
      REAL_T  fluxxlo, fluxxhi, fluxylo, fluxyhi, fluxzlo, fluxzhi
      REAL_T  a,b,c

      dt = 1.0e20
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               dtcell = dt
               if (flag.eq.1) then
                  if(divu(i,j,k).gt.zero) then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = dtfactor*(one-rhomin/rho(i,j,k))/divu(i,j,k)
                     else
                        dtcell = dtfactor*.5/divu(i,j,k)
                     endif
                     if (dsdt(i,j,k).gt.1.0e-20) then
                        if (abs(rho(i,j,k)).gt.rhomin) then
                           rhominij = rhomin
                        else
                           rhominij = .9*abs(rho(i,j,k)) 
                        endif
                        rhoij = abs(rho(i,j,k))
#if 0 
                        dtcell2 = (-rhoij*divu(i,j,k) + 
     &                       sqrt((rhoij*divu(i,j,k))**2+
     &                       two*dsdt(i,j,k)*rhoij*abs(rhoij-rhominij)))/
     &                       (rhoij*dsdt(i,j,k))
#else 
c     ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
c     We use the latter because it is more robust
                        a = rhoij*dsdt(i,j,k)*half
                        b = rhoij*divu(i,j,k)
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
                  denom = rho(i,j,k)*divu(i,j,k)+
     &                 u(i,j,k,1)*(rho(i+1,j,k)-rho(i-1,j,k))/delta(1) +
     &                 u(i,j,k,2)*(rho(i,j+1,k)-rho(i,j-1,k))/delta(2) +
     &                 u(i,j,k,3)*(rho(i,j,k+1)-rho(i,j,k-1))/delta(3)
                  if(denom.gt.zero)then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = dtfactor*(rho(i,j,k)-rhomin)/denom
                     else
                        dtcell = dtfactor*abs(rho(i,j,k))/denom
                     endif
                  endif
               else if (flag.eq.3) then
                  fluxxlo = fourth*(rho(i,j,k)+rho(i-1,j,k))
     &                              *(u(i,j,k,1)+u(i-1,j,k,1))
                  fluxxhi = fourth*(rho(i,j,k)+rho(i+1,j,k))
     &                              *(u(i,j,k,1)+u(i+1,j,k,1))
                  fluxylo = fourth*(rho(i,j,k)+rho(i,j-1,k))
     &                              *(u(i,j,k,2)+u(i,j-1,k,2))
                  fluxyhi = fourth*(rho(i,j,k)+rho(i,j+1,k))
     &                              *(u(i,j,k,2)+u(i,j+1,k,2))
                  fluxzhi = fourth*(rho(i,j,k)+rho(i,j,k-1))
     &                              *(u(i,j,k,3)+u(i,j,k-1,3))
                  fluxzlo = fourth*(rho(i,j,k)+rho(i,j,k+1))
     &                              *(u(i,j,k,3)+u(i,j,k+1,3))
                  denom = ((areax(i+1,j,k)*fluxxhi-areax(i,j,k)*fluxxlo)+
     &                     (areay(i,j+1,k)*fluxyhi-areay(i,j,k)*fluxylo)+
     &                     (areaz(i,j,k+1)*fluxzhi-areaz(i,j,k)*fluxzlo))
     &                      /volume(i,j,k)
                  
                  if(denom.gt.zero)then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = dtfactor*(rho(i,j,k)-rhomin)/denom
                     else
                        dtcell = dtfactor*abs(rho(i,j,k))/denom
                     endif
                  endif
               endif
#if 0 
               write(6,*)'i,j,k,dtcell=',i,j,k,dtcell
#endif
               dt = min(dtcell,dt)
            enddo
         enddo
      enddo
      end

      subroutine FORT_CHECK_DIVU_DT(flag, dtfactor, delta, divu, DIMS(divu), 
     &                           dsdt, rho, DIMS(rho), 
     &                           u, DIMS(u), 
     &                           volume, DIMS(volume), 
     &                           areax,  DIMS(areax),
     &                           areay,  DIMS(areay),
     &                           areaz,  DIMS(areaz),
     &                           lo, hi,
     &                           dt, rhomin)

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
      integer DIMDEC(areaz)
      REAL_T  volume(DIMV(volume))
      REAL_T  areax(DIMV(areax))
      REAL_T  areay(DIMV(areay))
      REAL_T  areaz(DIMV(areaz))

      integer i,j,k
      REAL_T  dtcell, denom
      REAL_T  fluxxlo,fluxxhi,fluxylo,fluxyhi,fluxzlo,fluxzhi
      REAL_T  a,b,c,dtcell2,rhominij,rhoij

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               dtcell = bigreal
               if (flag.eq.1) then
                  if(divu(i,j,k).gt.zero) then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = (one-rhomin/rho(i,j,k))/divu(i,j,k)
                     else
                        dtcell = one/divu(i,j,k)
                     endif
                     if (dsdt(i,j,k).gt.1.0e-20) then
                        if (abs(rho(i,j,k)).gt.rhomin) then
                           rhominij = rhomin
                        else
                           rhominij = .9*abs(rho(i,j,k)) 
                        endif
                        rhoij = abs(rho(i,j,k))
#if 0 
                        dtcell2 = (-rhoij*divu(i,j,k) + 
     &                       sqrt((rhoij*divu(i,j,k))**2+
     &                       two*dsdt(i,j,k)*rhoij*abs(rhoij-rhominij)))/
     &                       (rhoij*dsdt(i,j,k))
#else 
c     ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
c     We use the latter because it is more robust
                        a = rhoij*dsdt(i,j,k)*half
                        b = rhoij*divu(i,j,k)
                        c = rhominij - rhoij
                        dtcell2 = two*c/(-b-sqrt(b**2-four*a*c))
#endif
                        dtcell = min(dtcell,dtcell2)
                     endif
                  endif
               else if (flag.eq.2) then
                  denom = rho(i,j,k)*divu(i,j,k)+
     &                 u(i,j,k,1)*(rho(i+1,j,k)-rho(i-1,j,k))/delta(1) +
     &                 u(i,j,k,2)*(rho(i,j+1,k)-rho(i,j-1,k))/delta(2) +
     &                 u(i,j,k,3)*(rho(i,j,k+1)-rho(i,j,k-1))/delta(3)
                  if(denom.gt.zero)then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = (rho(i,j,k)-rhomin)/denom
                     else
                        dtcell = abs(rho(i,j,k))/denom
                     endif
                  endif
               else if (flag.eq.3) then
                  fluxxlo = fourth*(rho(i,j,k)+rho(i-1,j,k))
     &                              *(u(i,j,k,1)+u(i-1,j,k,1))
                  fluxxhi = fourth*(rho(i,j,k)+rho(i+1,j,k))
     &                              *(u(i,j,k,1)+u(i+1,j,k,1))
                  fluxylo = fourth*(rho(i,j,k)+rho(i,j-1,k))
     &                              *(u(i,j,k,2)+u(i,j-1,k,2))
                  fluxyhi = fourth*(rho(i,j,k)+rho(i,j+1,k))
     &                              *(u(i,j,k,2)+u(i,j+1,k,2))
                  fluxzlo = fourth*(rho(i,j,k)+rho(i,j,k-1))
     &                              *(u(i,j,k,3)+u(i,j,k-1,3))
                  fluxzhi = fourth*(rho(i,j,k)+rho(i,j,k+1))
     &                              *(u(i,j,k,3)+u(i,j,k+1,3))
                  denom = ((areax(i+1,j,k)*fluxxhi-areax(i,j,k)*fluxxlo)+
     &                     (areay(i,j+1,k)*fluxyhi-areay(i,j,k)*fluxylo)+
     &                     (areaz(i,j,k+1)*fluxzhi-areaz(i,j,k)*fluxzlo))
     &                      /volume(i,j,k)
                  if(denom.gt.zero)then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = (rho(i,j,k)-rhomin)/denom
                     else
                        dtcell = abs(rho(i,j,k))/denom
                     endif
                  endif
               endif
               if (dt.gt.dtcell) then
                  write(6,*)'ERROR: FORT_CHECK_DIVU_DT : i,j,k,dt>dtcell = ',
     &                 i,j,k,dt,dtcell
               else if (dt.gt.dtcell*dtfactor) then
                  write(6,*)'WARNING: FORT_CHECK_DIVU_DT : i,j,k,dt>dtcell*dtfactor = ',
     &                 i,j,k,dt,dtcell*dtfactor
               endif
            enddo
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

      integer    DIMDEC(dqrad)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     dqrad(DIMV(dqrad))

      call filcc (dqrad,DIMS(dqrad),domlo,domhi,delta,xlo,bc)
      call fillEdges(dqrad,DIMS(dqrad),domlo,domhi,bc)
      end

      subroutine FORT_DIVUFILL (divu,DIMS(divu),domlo,domhi,delta,
     &                        xlo,time,bc )

      integer    DIMDEC(divu)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     divu(DIMV(divu))

      call filcc (divu,DIMS(divu),domlo,domhi,delta,xlo,bc)

      call fillEdges(divu,DIMS(divu),domlo,domhi,bc)

      end

      subroutine FORT_DSDTFILL (dsdt,DIMS(dsdt),domlo,domhi,delta,
     &                         xlo,time,bc )

      integer    DIMDEC(dsdt)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     dsdt(DIMV(dsdt))

      call filcc (dsdt,DIMS(dsdt),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(dsdt,DIMS(dsdt),domlo,domhi,bc)

      end

      subroutine FORT_YDOTFILL (ydot,DIMS(ydot),domlo,domhi,delta,
     &                          xlo,time,bc)

      integer    DIMDEC(ydot), bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     ydot(DIMV(ydot))

      call filcc (ydot,DIMS(ydot),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(ydot,DIMS(ydot),domlo,domhi,bc)

      end

      subroutine FORT_FABMINMAX(lo, hi, 
     &                          fab, DIMS(fab),
     &                          fmin, fmax, nc)
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(fab)
      REAL_T  fab(DIMV(fab),nc)
      REAL_T  fmin, fmax

      integer i,j,k,n

      do n = 1,nc
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  fab(i,j,k,n) = MAX( fmin, MIN( fmax, fab(i,j,k,n) ) )
               end do
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
      integer i, j, k, n, pm, maxs, iN2
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
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  sum = zero
                  sumRhoY = zero
                  do n=1,nc
                     sum = sum + flux(i,j,k,n)
                     sumRhoY = sumRhoY + half*(state(i-1,j,k,n)+state(i,j,k,n))
                  end do               
#if (METH == VCOR)
                  do n=1,nc
                     yedge = half*(state(i-1,j,k,n)+state(i,j,k,n))
                     flux(i,j,k,n) = flux(i,j,k,n) - sum*yedge/sumRhoY
                  end do
#else 
#  if (METH == N2DOM)
                  maxs = iN2
#  elif (METH == FINDDOM)
                  pm = cvmgt(0,-1,sum.gt.zero)
                  maxs=1
                  maxv=state(i+pm,j,k,maxs)
                  do n=2,nc
                     val = state(i+pm,j,k,n)
                     if (val.gt.maxv) then
                        maxv=val
                        maxs=n
                     end if
                  end do
#  endif
                  flux(i,j,k,maxs) = flux(i,j,k,maxs) - sum
#endif
               end do
            end do
         end do
      else if (dir.eq.1) then
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  sum = zero
                  sumRhoY = zero
                  do n=1,nc
                     sum = sum + flux(i,j,k,n)
                     sumRhoY = sumRhoY + half*(state(i,j-1,k,n)+state(i,j,k,n))
                  end do
#if (METH == VCOR)
                  do n=1,nc
                     yedge = half*(state(i,j-1,k,n)+state(i,j,k,n))
                     flux(i,j,k,n) = flux(i,j,k,n) - sum*yedge/sumRhoY
                  end do
#else 
#  if (METH == N2DOM)
                  maxs = iN2
#  elif (METH == FINDDOM)
                  pm = cvmgt(0,-1,sum.gt.zero)
                  maxs=1
                  maxv=state(i,j+pm,k,maxs)
                  do n=2,nc
                     val = state(i,j+pm,k,n)
                     if (val.gt.maxv) then
                        maxv=val
                        maxs=n
                     end if
                  end do
#  endif
                  flux(i,j,k,maxs) = flux(i,j,k,maxs) - sum
#endif
               end do
            end do
         end do
      else
         do k = lo(3),hi(3)+1
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  sum = zero
                  sumRhoY = zero
                  do n=1,nc
                     sum = sum + flux(i,j,k,n)
                     sumRhoY = sumRhoY + half*(state(i,j,k-1,n)+state(i,j,k,n))
                  end do
#if (METH == VCOR)
                  do n=1,nc
                     yedge = half*(state(i,j,k-1,n)+state(i,j,k,n))
                     flux(i,j,k,n) = flux(i,j,k,n) - sum*yedge/sumRhoY
                  end do
#else 
#  if (METH == N2DOM)
                  maxs = iN2
#  elif (METH == FINDDOM)
                  pm = cvmgt(0,-1,sum.gt.zero)
                  maxs=1
                  maxv=state(i,j,k+pm,maxs)
                  do n=2,nc
                     val = state(i,j,k+pm,n)
                     if (val.gt.maxv) then
                        maxv=val
                        maxs=n
                     end if
                  end do
#  endif
                  flux(i,j,k,maxs) = flux(i,j,k,maxs) - sum
#endif
               end do
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
      integer i,j,k,n,hi_e(SDIM),ii,jj,kk
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
         hi_e(3) = hi(3)
      else if (dir .EQ. 1) then
         hi_e(1) = hi(1)
         hi_e(2) = hi(2) + 1
         hi_e(3) = hi(3)
      else if (dir .EQ. 2) then
         hi_e(1) = hi(1)
         hi_e(2) = hi(2)
         hi_e(3) = hi(3) + 1
      else
         write(6,*) 'FORT_REPAIR_FLUX: Bad dir'
      end if
      
      do k = lo(3),hi_e(3)
         do j = lo(2),hi_e(2)
            do i = lo(1),hi_e(1)
               sump=zero
               sumn=zero
               do n=1,Nspec
                  if (flux(i,j,k,n) .GT. zero) then
                     sump = sump + flux(i,j,k,n)
                  else if (flux(i,j,k,n) .LT. zero) then
                     sumn = sumn + flux(i,j,k,n)
                  end if
               end do
               delta = sump + sumn
               if (delta .NE. zero) then
                  sumy = zero
                  if (delta .LT. zero) then
                     do n=1,Nspec
                        if (dir .EQ. 0) then
                           y(n) = MAX(zero,state(i-1,j,k,n))
                        else if (dir .EQ. 1) then
                           y(n) = MAX(zero,state(i,j-1,k,n))
                        else
                           y(n) = MAX(zero,state(i,j,k-1,n))
                        endif
                        if (flux(i,j,k,n).GE.zero) y(n) = zero
                        sumy = sumy + y(n)
                     end do
                     alphabar = - sump/sumn
                     beta=cvmgt(sumn*(one-alphabar)/sumy,zero,sumy.NE.zero)
                  else
                     do n=1,Nspec
                        y(n) = MAX(zero,state(i,j,k,n))
                        if (flux(i,j,k,n).LE.zero) y(n) = zero
                        sumy = sumy + y(n)
                     end do
                     alphabar = - sumn/sump
                     beta=cvmgt(sump*(one-alphabar)/sumy,zero,sumy.NE.zero)
                  end if
                  if (alphabar .GT. one) goto 100
                  do n=1,Nspec
                     flux(i,j,k,n) = flux(i,j,k,n) - beta*y(n)
                  end do
               end if
            end do
         end do
      end do
      return
      
 100  continue 
      write(6,*) 'FORT_REPAIR_FLUX: alphabar > 1'
      write(6,*) 'fluxes: (dir=',dir,')'
      do n=1,Nspec
         call get_spec_name(name,n)
         write(6,996) name,flux(i,j,k,n)
      end do
      write(6,*) 'left state:'
      do n=1,Nspec
         call get_spec_name(name,n)
         write(6,996) name,state(i,j,k,n)
      end do
      write(6,*) 'right state:'
      ii = i
      jj = j
      kk = k
      if (dir .EQ. 0) then
         ii = i + 1
      else if (dir .EQ. 1) then
         jj = j + 1
      else
         kk = k + 1
      end if
      do n=1,Nspec
         call get_spec_name(name,n)
         write(6,996) name,state(ii,jj,kk,n)
      end do
 996  format(a16,1x,e30.22)
      call bl_abort(" ")
      end
#endif
      
      subroutine  FORT_RECOMP_UPDATE(lo, hi,
     &                               update, DIMS(update),
     &                               xflux,  DIMS(xflux),
     &                               yflux,  DIMS(yflux),
     &                               zflux,  DIMS(zflux),
     &                               vol,    DIMS(vol),
     &                               nc)
      implicit none
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(update)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(vol)
      REAL_T update(DIMV(update),nc)
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T zflux(DIMV(zflux),nc)
      REAL_T vol(DIMV(vol))

      integer i, j, k, n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,nc
                  update(i,j,k,n)=-((xflux(i+1,j,k,n)-xflux(i,j,k,n))
     &                 +            (yflux(i,j+1,k,n)-yflux(i,j,k,n))
     &                 +            (zflux(i,j,k+1,n)-zflux(i,j,k,n)))
     &                     /vol(i,j,k)
               end do
            end do
         end do
      end do
      end


      subroutine  FORT_INCRWEXTFLXDIV(lo, hi,
     &                                xflux,  DIMS(xflux),
     &                                yflux,  DIMS(yflux),
     &                                zflux,  DIMS(zflux),
     &                                stateo, DIMS(stateo),
     &                                staten, DIMS(staten),
     &                                vol,    DIMS(vol),
     &                                nc, dt)
      implicit none
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(stateo)
      integer DIMDEC(staten)
      integer DIMDEC(vol)
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T zflux(DIMV(zflux),nc)
      REAL_T stateo(DIMV(stateo))
      REAL_T staten(DIMV(staten))
      REAL_T vol(DIMV(vol))
      REAL_T dt

      integer i, j, k, n
      REAL_T dF
      
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               dF = zero
               do n=1,nc
                  dF = dF + ( xflux(i+1,j,k,n) - xflux(i,j,k,n) )
     &                 +    ( yflux(i,j+1,k,n) - yflux(i,j,k,n) )
     &                 +    ( zflux(i,j,k+1,n) - zflux(i,j,k,n) )
               end do
               staten(i,j,k) = stateo(i,j,k) + dF*dt/vol(i,j,k)
            end do
         end do
      end do
      end

      subroutine fillEdges(dat,DIMS(dat),domlo,domhi,bc )
c  This routine fills ghost cells with the value from the nearest
c  interior cell.  
      integer    DIMDEC(dat)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dat(DIMV(dat))

      integer    lo(SDIM), hi(SDIM)
      integer    i, j, k
      integer    ilo, ihi, jlo, jhi, klo, khi
      logical    do_xlo,do_xhi,do_ylo,do_yhi,do_zlo,do_zhi

      lo(1) = ARG_L1(dat)
      hi(1) = ARG_H1(dat)
      lo(2) = ARG_L2(dat)
      hi(2) = ARG_H2(dat)
      lo(3) = ARG_L3(dat)
      hi(3) = ARG_H3(dat)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      klo = max(lo(3),domlo(3))
      khi = min(hi(3),domhi(3))

      do_xlo = lo(1).lt.domlo(1)
      do_xhi = hi(1).gt.domhi(1)
      do_ylo = lo(2).lt.domlo(2)
      do_yhi = hi(2).gt.domhi(2)
      do_zlo = lo(3).lt.domlo(3)
      do_zhi = hi(3).gt.domhi(3)

      if (bc(1,1).eq.EXT_DIR.and. do_xlo) then
         do k = klo,khi
            do j = jlo, jhi
               do i = lo(1), domlo(1)-1
                  dat(i,j,k) = dat(domlo(1),j,k)
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = klo, khi
               do j = lo(2), domlo(2)-1
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domlo(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = lo(2), domlo(2)-1
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = lo(2), domlo(2)-1
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_yhi) then
            do k = klo, khi
               do j = domhi(2)+1, hi(2)
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domhi(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = domhi(2)+1,hi(2)
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = domhi(2)+1,hi(2)
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_zlo) then
            do k = lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, hi(3)
               do j = jlo, jhi
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
      endif            


      if (bc(1,2).eq.EXT_DIR.and. do_xhi) then
         do k = klo, khi
            do j = jlo, jhi
               do i = domhi(1)+1,hi(1)
                  dat(i,j,k) = dat(domhi(1),j,k)
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = klo, khi
               do j = lo(2), domlo(2)-1
                  do i = domhi(1)+1,hi(1)
                     dat(i,j,k) = dat(domhi(1),domlo(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = lo(2), domlo(2)-1
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = lo(2), domlo(2)-1
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_yhi) then
            do k = klo, khi
               do j = domhi(2)+1, hi(2)
                  do i = domhi(1)+1,hi(1)
                     dat(i,j,k) = dat(domhi(1),domhi(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = domhi(2)+1,hi(2)
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = domhi(2)+1,hi(2)
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_zlo) then
            do k = lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = domhi(1)+1,hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, hi(3)
               do j = jlo, jhi
                  do i = domhi(1)+1,hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
      endif            

      if (bc(2,1).eq.EXT_DIR.and.do_ylo) then
         do k = klo, khi
            do j = lo(2), domlo(2)-1
               do i = ilo,ihi
                  dat(i,j,k) = dat(i,domlo(2),k)
               enddo
            enddo
         enddo
         if (do_xlo) then
            do k = klo, khi
               do j = lo(2), domlo(2)-1
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domlo(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = lo(2), domlo(2)-1
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = lo(2), domlo(2)-1
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if(do_xhi)then
            do k = klo, khi
               do j = lo(2), domlo(2)-1
                  do i = domhi(1)+1, hi(1)
                     dat(i,j,k) = dat(domhi(1),domlo(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = lo(2),domlo(2)-1
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = lo(2),domlo(2)-1
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_zlo) then
            do k = lo(3), domlo(3)-1
               do j = lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, hi(3)
               do j = lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domhi(3))
                  enddo
               enddo
            enddo
         endif
      endif            
      
      if (bc(2,2).eq.EXT_DIR.and.do_yhi) then
         do k = klo, khi
            do j = domhi(2)+1, hi(2)
               do i = ilo, ihi
                  dat(i,j,k) = dat(i,domhi(2),k)
               enddo
            enddo
         enddo
         if (do_xlo) then
            do k = klo, khi
               do j = domhi(2)+1, hi(2)
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domhi(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = domhi(2)+1, hi(2)
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = domhi(2)+1, hi(2)
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if(do_xhi)then
            do k = klo, khi
               do j = domhi(2)+1, hi(2)
                  do i = domhi(1)+1, hi(1)
                     dat(i,j,k) = dat(domhi(1),domhi(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = domhi(2)+1, hi(2)
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = domhi(2)+1, hi(2)
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_zlo) then
            do k = lo(3), domlo(3)-1
               do j = domhi(2)+1, hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, hi(3)
               do j = domhi(2)+1, hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domhi(3))
                  enddo
               enddo
            enddo
         endif
      endif            

      if (bc(3,1).eq.EXT_DIR.and. do_zlo) then
         do k = lo(3), domlo(3)-1
            do j = jlo, jhi
               do i = ilo, ihi
                  dat(i,j,k) = dat(i,j,domlo(3))
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = lo(3), domlo(3)-1
               do j = lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domlo(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = lo(3), domlo(3)-1
                  do j = lo(2), domlo(2)-1
                     do i = lo(1),domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = lo(3), domlo(3)-1
                  do j = lo(2), domlo(2)-1
                     do i = domhi(1)+1, hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_yhi) then
            do k = lo(3), domlo(3)-1
               do j = domhi(2)+1, hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domlo(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = lo(3), domlo(3)-1
                  do j = domhi(2)+1,hi(2)
                     do i = lo(1),domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = lo(3), domlo(3)-1
                  do j = domhi(2)+1,hi(2)
                     do i = domhi(1)+1, hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_xlo) then
            do k = lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_xhi)then
            do k = lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = domhi(1)+1, hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
      endif            

      if (bc(3,2).eq.EXT_DIR.and. do_zhi) then
         do k = domhi(3)+1,hi(3)
            do j = jlo, jhi
               do i = ilo,ihi
                  dat(i,j,k) = dat(i,j,domhi(3))
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = domhi(3)+1,hi(3)
               do j = lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domhi(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = domhi(3)+1,hi(3)
                  do j = lo(2), domlo(2)-1
                     do i = lo(1),domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = domhi(3)+1,hi(3)
                  do j = lo(2), domlo(2)-1
                     do i = domhi(1)+1, hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_yhi) then
            do k = domhi(3)+1,hi(3)
               do j = domhi(2)+1, hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domhi(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = domhi(3)+1,hi(3)
                  do j = domhi(2)+1,hi(2)
                     do i = lo(1),domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = domhi(3)+1,hi(3)
                  do j = domhi(2)+1,hi(2)
                     do i = domhi(1)+1, hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_xlo) then
            do k = domhi(3)+1,hi(3)
               do j = jlo, jhi
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
         if(do_xhi)then
            do k = domhi(3)+1,hi(3)
               do j = jlo, jhi
                  do i = domhi(1)+1, hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
      endif            
      end

      subroutine fillWithZeros(dat,DIMS(dat),domlo,domhi,bc)

      integer    DIMDEC(dat)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dat(DIMV(dat))

      integer    lo(SDIM), hi(SDIM)
      integer    i, j, k

      lo(1) = ARG_L1(dat)
      hi(1) = ARG_H1(dat)
      lo(2) = ARG_L2(dat)
      hi(2) = ARG_H2(dat)
      lo(3) = ARG_L3(dat)
      hi(3) = ARG_H3(dat)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), domlo(1)-1
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif            
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = domhi(1)+1,hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif            
      
      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do k=lo(3),hi(3)
            do j = lo(2), domlo(2)-1
               do i = lo(1), hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif            
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do k=lo(3),hi(3)
            do j = domhi(2)+1, hi(2)
               do i = lo(1), hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif            
      
      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k=lo(3),domlo(3)-1
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif            
      
      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k=domhi(3)+1,hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif            
      end

      subroutine  FORT_COMPUTE_UGRADP(p, DIMS(p), ugradp, DIMS(ugp),
     &                                umac,  DIMS(umac),
     &                                vmac,  DIMS(vmac),
     &                                wmac,  DIMS(wmac),
     &                                lo, hi, dx)

      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(p)
      integer DIMDEC(ugp)
      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(wmac)
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T  wmac(DIMV(wmac))
      REAL_T      p(DIMV(p))
      REAL_T ugradp(DIMV(ugp))
      REAL_T dx(SDIM)

      integer i, j, k
      REAL_T uadv, vadv, wadv
      REAL_T p_x_lo, p_x_hi
      REAL_T p_y_lo, p_y_hi
      REAL_T p_z_lo, p_z_hi
      
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            uadv = half*(umac(i,j,k) + umac(i+1,j,k))
            vadv = half*(vmac(i,j,k) + vmac(i,j+1,k))
            wadv = half*(wmac(i,j,k) + wmac(i,j,k+1))
            p_x_hi = cvmgp(p(i  ,j,k),p(i+1,j,k),umac(i+1,j,k))
            p_x_lo = cvmgp(p(i-1,j,k),p(i  ,j,k),umac(i  ,j,k))
            p_y_hi = cvmgp(p(i,j  ,k),p(i,j+1,k),vmac(i,j+1,k))
            p_y_lo = cvmgp(p(i,j-1,k),p(i,j  ,k),vmac(i,j  ,k))
            p_z_hi = cvmgp(p(i,j,k  ),p(i,j,k+1),wmac(i,j,k+1))
            p_z_lo = cvmgp(p(i,j,k-1),p(i,j,k  ),wmac(i,j,k  ))
            ugradp(i,j,k) = uadv * (p_x_hi - p_x_lo) / dx(1) +
     &                      vadv * (p_y_hi - p_y_lo) / dx(2) +    
     &                      wadv * (p_z_hi - p_z_lo) / dx(3)
         end do
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
      integer n,i,j,k, loC(SDIM),hiC(SDIM),ii,jj,kk,iii,jjj,kkk,ncells
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

      do k=loC(3),hiC(3)
         do j=loC(2),hiC(2)
            do i=loC(1),hiC(1)

               bad_T = .false.
               do kk=0,ratio(3)-1
                  kkk = ratio(3)*k + kk
                  do jj=0,ratio(2)-1
                     jjj = ratio(2)*j + jj
                     do ii=0,ratio(1)-1
                        iii = ratio(1)*i + ii
                        if (state(iii,jjj,kkk,Tcomp).lt.min_T) then
                           bad_T = .true.
                        endif
                     enddo
                  enddo
               enddo

               if (bad_T .eqv. .true.) then

                  tmp(Rcomp) = 0.d0
                  do n=first_spec,last_spec
                     tmp(n) = 0.d0
                  enddo
                  tmp(RhoH) = 0.d0


                  do kk=0,ratio(3)-1
                     kkk = ratio(3)*k + kk
                     do jj=0,ratio(2)-1
                        jjj = ratio(2)*j + jj
                        do ii=0,ratio(1)-1
                           iii = ratio(1)*i + ii
                           
                           tmp(Rcomp) = tmp(Rcomp) + state(iii,jjj,kkk,Rcomp)
                           do n=first_spec,last_spec
                              tmp(n) = tmp(n) + state(iii,jjj,kkk,n)
                           enddo
                           tmp(RhoH) = tmp(RhoH) + state(iii,jjj,kkk,RhoH)
                           
                        enddo
                     enddo
                  enddo

                  FORT_CONSERVATIVE_T_FLOOR = FORT_CONSERVATIVE_T_FLOOR + ncells
                  tmp(Rcomp) = tmp(Rcomp) * ncellsInv
                  do n=first_spec,last_spec
                     tmp(n) = tmp(n) * ncellsInv
                  enddo
                  tmp(RhoH) = tmp(RhoH)* ncellsInv
                  
                  do kk=0,ratio(3)-1
                     kkk = ratio(3)*k + kk
                     do jj=0,ratio(2)-1
                        jjj = ratio(2)*j + jj
                        do ii=0,ratio(1)-1
                           iii = ratio(1)*i + ii
                           
                           state(iii,jjj,kkk,Rcomp) = tmp(Rcomp)
                           do n=first_spec,last_spec
                              state(iii,jjj,kkk,n) = tmp(n)
                           enddo
                           state(iii,jjj,kkk,RhoH) = tmp(RhoH)
                           
                        enddo
                     enddo
                  enddo


               endif

            enddo
         enddo
      enddo
      end
