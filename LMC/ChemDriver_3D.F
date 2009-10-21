#include "ChemDriver_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"

c     Series of runtime hacks
#define DO_JBB_HACK
#undef DO_JBB_HACK_TEMP
#define DO_JBB_HACK_TEMP
#define HACK_TEMP_MIN 273.d0
#define HACK_TEMP_MAX 1900.d0
#define TRIGGER_NEW_J
#undef ALWAYS_NEW_J
#undef SOLN_IS_1D

#define CONPF_FILE conpFY
#define CONPJ_FILE conpJY

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define twothousand 2000
#define one100th    0.01
#define ten2minus19 1.e-19
#define million     1.e6
#define one2minus3  1.e-3
#else
#define twothousand 2000d0
#define one100th    0.01d0
#define ten2minus19 1.d-19
#define million     1.d6
#define one2minus3  1.d-3
#endif

#define SDIM 3

      subroutine FORT_NORMMASS(lo, hi, xsID,
     &                         Y, DIMS(Y), Ynorm, DIMS(YNORM))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(Ynorm)
      integer xsID
      REAL_T Y(DIMV(Y),*)
      REAL_T Ynorm(DIMV(Ynorm),*)

      integer i, j, k, n
      REAL_T sum

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               sum = zero
               do n=1,Nspec
                  Ynorm(i,j,k,n) =  MAX( Y(i,j,k,n),zero)
                  sum = sum + Ynorm(i,j,k,n)
               end do
               Ynorm(i,j,k,xsID) = Y(i,j,k,xsID)+ one - sum
            end do
         end do
      end do
      end

      subroutine FORT_FRrateXTP(lo,hi,X,DIMS(X),T,DIMS(T),
     &                          FwdK,DIMS(FwdK),RevK,DIMS(RevK),
     &                          Patm,rxns,Nrxns)
      implicit none
#include "cdwrk.H"
#include "conp.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(X)
      integer DIMDEC(T)
      integer DIMDEC(FwdK)
      integer DIMDEC(RevK)
      integer Nrxns
      integer rxns(Nrxns)
      REAL_T X(DIMV(X),*)
      REAL_T T(DIMV(T))
      REAL_T FwdK(DIMV(FwdK),*)
      REAL_T RevK(DIMV(RevK),*)
      REAL_T Patm, scale

      REAL_T Xt(maxspec),FwdKt(maxreac),RevKt(maxreac)
      REAL_T sum, Yt(maxspec)
      integer i,j,k,n
      REAL_T P1atm,RU,RUC,Pdyne

      parameter(scale = million)

      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      Pdyne = Patm * P1atm

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Xt(n) = X(i,j,k,n)
               end do
#ifdef DO_JBB_HACK
               CALL CKXTY(Xt,IWRK(ckbi),RWRK(ckbr),Yt)
               sum = zero
               do n=1,Nspec
                  Yt(n) =MAX( Yt(n),zero)
                  sum = sum+Yt(n)
               end do
               if (iN2 .gt. 0) then
                  Yt(iN2) = Yt(iN2)+one-sum
               endif
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),Xt)
#endif
#ifdef MIKE1
               CALL CKKFKR(Pdyne,T(i,j,k),Xt,IWRK(ckbi),RWRK(ckbr),FwdKt,RevKt)
#else
               call bl_abort('FORT_FRrateXTP: CKKFKR not available')
#endif
               do n=1,Nrxns
                  FwdK(i,j,k,n) = FwdKt(rxns(n)+1)*scale
                  RevK(i,j,k,n) = RevKt(rxns(n)+1)*scale
               end do
            end do
         end do
      end do
      end

      subroutine FORT_HTRLS(lo,hi,Y,DIMS(Y),T,DIMS(T),
     &                      Q,DIMS(Q),Patm)
      implicit none
#include "cdwrk.H"
#include "conp.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      integer DIMDEC(Q)
      REAL_T Y(DIMV(Y),*)
      REAL_T T(DIMV(T))
      REAL_T Q(DIMV(Q))
      REAL_T Patm

      REAL_T Zt(maxspec+1),Zdott(maxspec+1)
      integer i,j,k,n
      integer ndummy
      REAL_T tdummy,P1atm,RU,RUC
      REAL_T RHO, CPB, scal
c
c     NOTE: scal converts result from assumed cgs to MKS (1 erg/s.cm^3 = .1 J/s.m^3)
c
      parameter(scal = tenth)

      ndummy = Nspec
      tdummy = 0.0D0
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RWRK(NP) = Patm * P1atm

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Zt(1) = T(i,j,k)
               do n=1,Nspec
                  Zt(n+1) = Y(i,j,k,n)
               end do
               call conpFY(ndummy,tdummy,Zt,Zdott,RWRK,IWRK)
               CALL CKRHOY(RWRK(NP),Zt(1),Zt(2),IWRK(ckbi),RWRK(ckbr),RHO)
               CALL CKCPBS(Zt(1),Zt(2),IWRK(ckbi),RWRK(ckbr),CPB)
               Q(i,j,k) = Zdott(1) * RHO * CPB * scal
            end do
         end do
      end do
      end

      subroutine FORT_RRATEY(lo,hi,Y,DIMS(Y),T,DIMS(T),
     &                       Ydot,DIMS(Ydot),Patm)
      implicit none
#include "cdwrk.H"
#include "conp.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      integer DIMDEC(Ydot)
      REAL_T Y(DIMV(Y),*)
      REAL_T T(DIMV(T))
      REAL_T Ydot(DIMV(Ydot),*)
      REAL_T Patm

      REAL_T Zt(maxspec+1),Zdott(maxspec+1)
      integer i,j,k,n
      integer ndummy
      REAL_T tdummy,P1atm,RU,RUC

      ndummy = Nspec
      tdummy = 0.
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RWRK(NP) = Patm * P1atm

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Zt(1) = T(i,j,k)
               do n=1,Nspec+1
                  Zt(n+1) = Y(i,j,k,n)
               end do
               call conpFY(ndummy,tdummy,Zt,Zdott,RWRK,IWRK)
               do n=1,Nspec
                  Ydot(i,j,k,n) = Zdott(n+1)
               end do
            end do
         end do
      end do
      end

      subroutine FORT_RRATEC(lo,hi,C,DIMS(C),T,DIMS(T),
     &                       Cdot,DIMS(Cdot),Patm)
      implicit none
#include "cdwrk.H"
#include "conp.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(C)
      integer DIMDEC(T)
      integer DIMDEC(Cdot)
      REAL_T C(DIMV(C),*)
      REAL_T T(DIMV(T))
      REAL_T Cdot(DIMV(Cdot),*)
      REAL_T Patm

      REAL_T Zt(maxspec+1),Zdott(maxspec+1)
      integer i,j,k,n
      integer ndummy
      REAL_T tdummy,P1atm,RU,RUC

      ndummy = Nspec
      tdummy = 0.
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RWRK(NP) = Patm * P1atm

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Zt(1) = T(i,j,k)
               do n=1,Nspec+1
                  Zt(n+1) = C(i,j,k,n)
               end do
               call conpFC(ndummy,tdummy,Zt,Zdott,RWRK,IWRK)
               do n=1,Nspec
                  Cdot(i,j,k,n) = Zdott(n+1)
               end do
            end do
         end do
      end do
      end

      subroutine FORT_MASSTOMOLE(lo, hi, Y, DIMS(Y), X, DIMS(X))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(X)
      REAL_T Y(DIMV(Y),*)
      REAL_T X(DIMV(X),*)

      REAL_T Xt(maxspec), Yt(maxspec)
      integer i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),Xt)
               do n = 1,Nspec
                  X(i,j,k,n) = Xt(n)
               end do
            end do
         end do
      end do
      end
      
      subroutine FORT_MOLETOMASS(lo, hi, X, DIMS(X), Y, DIMS(Y))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(X)
      integer DIMDEC(Y)
      REAL_T X(DIMV(X),*)
      REAL_T Y(DIMV(Y),*)
      
      REAL_T Xt(maxspec), Yt(maxspec)
      integer i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Xt(n) = X(i,j,k,n)
               end do
               CALL CKXTY(Xt,IWRK(ckbi),RWRK(ckbr),Yt)
               do n = 1,Nspec
                  Y(i,j,k,n) = Yt(n)
               end do
            end do
         end do
      end do
      end

      subroutine FORT_MASSTP_TO_CONC(lo, hi, Patm,
     &                           Y, DIMS(Y), T, DIMS(T), C, DIMS(C))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      integer DIMDEC(C)
      REAL_T Patm
      REAL_T Y(DIMV(Y),*)
      REAL_T T(DIMV(T))
      REAL_T C(DIMV(C),*)
      
      REAL_T Yt(maxspec), Ct(maxspec), RU, RUC, P1ATM, Ptmp
      integer i,j,k,n

      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKYTCP(Ptmp,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),Ct)
               do n = 1,Nspec
                  C(i,j,k,n) = Ct(n)*million
               end do
            end do
         end do
      end do
      end

      subroutine FORT_MASSR_TO_CONC(lo, hi, Y, DIMS(Y), 
     &                              T, DIMS(T), RHO, DIMS(RHO), C, DIMS(C))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      integer DIMDEC(C)
      integer DIMDEC(RHO)
      REAL_T Y(DIMV(Y),*)
      REAL_T T(DIMV(T))
      REAL_T C(DIMV(C),*)
      REAL_T RHO(DIMV(RHO))

      REAL_T Yt(maxspec), Ct(maxspec), rhoScl
      integer i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               rhoScl = RHO(i,j,k)*one2minus3
               CALL CKYTCR(rhoScl,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),Ct)
               do n = 1,Nspec
                  C(i,j,k,n) = Ct(n)*million
               end do
            end do
         end do
      end do
      end

      subroutine FORT_CONC_TO_MOLE(lo, hi,
     &                             C, DIMS(C), X, DIMS(X))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(C)
      integer DIMDEC(X)
      REAL_T C(DIMV(C),*)
      REAL_T X(DIMV(X),*)

      REAL_T Ct(maxspec), Xt(maxspec), scale
      integer i,j,k,n

      parameter(scale = one/million)
      
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Ct(n) = C(i,j,k,n)*scale
               end do
               CALL CKCTX(Ct,IWRK(ckbi),RWRK(ckbr),Xt)
               do n = 1,Nspec
                  X(i,j,k,n) = Xt(n)
               end do
            end do
         end do
      end do
      end

      subroutine FORT_MOLPROD(lo, hi, id, 
     &                        Q, DIMS(Q), C, DIMS(C), T, DIMS(T) )
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM), id
      integer DIMDEC(Q)
      integer DIMDEC(C)
      integer DIMDEC(T)
      REAL_T Q(DIMV(Q),*)
      REAL_T C(DIMV(C),*)
      REAL_T T(DIMV(T))

      REAL_T Ct(maxspec), Qt(maxreac), Qkt(maxreac), millionth
      integer i,j,k,n

      millionth = one/million
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Ct(n) = C(i,j,k,n)*millionth
               end do
               CALL CKQC(T(i,j,k),Ct,IWRK(ckbi),RWRK(ckbr),Qt)
#ifdef MIKE
               CALL CKCONT(id,Qt,IWRK(ckbi),RWRK(ckbr),Qkt)
#else
               call bl_abort('FORT_MOLPROD: CKCONT not available')
#endif
               do n = 1,Nreac
                  Q(i,j,k,n) = Qkt(n)*million
               end do
            end do
         end do
      end do
      end

      subroutine FORT_GETELTMOLES(namenc, namlen, lo, hi,
     &                            Celt, DIMS(Celt), C, DIMS(C))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer namlen, maxlen
      integer namenc(namlen)
      integer DIMDEC(Celt)
      integer DIMDEC(C)
      REAL_T Celt(DIMV(Celt))
      REAL_T C(DIMV(C),*)
      integer thenames(maxelts*2)
      logical match
      integer i, j, k, theidx, n, lout
      integer NCF(Nelt,Nspec)
c     Find index of desired element
      CALL CKSYME(thenames,2)
      theidx = -1
      do i=1,Nelt
         match = .true.
         do j=1,namlen               
            if (namenc(j) .NE. thenames((i-1)*2+j)) match = .false.
         enddo
         if (match .eqv. .true.) theidx = i
      end do
      if (theidx.lt.0) then
         call bl_pd_abort()
      endif
c     Get the matrix of elements versus species
      call CKNCF(Nelt,IWRK,RWRK,NCF)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               Celt(i,j,k) = zero
               do n = 1,Nspec
                  Celt(i,j,k) = Celt(i,j,k) +C(i,j,k,n)*NCF(theidx,n)
               end do
            end do
         end do
      end do
      end

      subroutine FORT_CONPSOLV(lo, hi,
     &                         Ynew, DIMS(Ynew), Tnew, DIMS(Tnew),
     &                         Yold, DIMS(Yold), Told, DIMS(Told),
     &                         FuncCount, DIMS(FuncCount),
     &                         Patm,
     &                         dt,
     &                         diag,do_diag)
      implicit none

#include "cdwrk.H"
#include "conp.H"
#include "vode.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(Yold)
      integer DIMDEC(Told)
      integer DIMDEC(Ynew)
      integer DIMDEC(Tnew)
      integer DIMDEC(FuncCount)
      integer do_diag
      REAL_T Yold(DIMV(Yold),*)
      REAL_T Told(DIMV(Told))
      REAL_T Ynew(DIMV(Ynew),*)
      REAL_T Tnew(DIMV(Tnew))
      REAL_T FuncCount(DIMV(FuncCount))
      REAL_T Patm, dt
      REAL_T diag(DIMV(FuncCount),*)
      
      integer ITOL, IOPT, ITASK, open_vode_failure_file
      parameter (ITOL=1, IOPT=1, ITASK=1)
      REAL_T RTOL, ATOL(maxspec+1),ATOLEPS
      REAL_T spec_scalT, NEWJ_TOL
      parameter (RTOL=1.0D-8, ATOLEPS=1.0D-8)
      parameter (spec_scalT=twothousand, NEWJ_TOL=one100th)
      external CONPF_FILE, CONPJ_FILE, open_vode_failure_file
      REAL_T TT1, TT2, RU, RUC, P1atm, sum, atoln
      integer i, j, k, m, MF, ISTATE, LOUTCK, lout

      integer nsubchem, nsub, node, tid, nthrds
      REAL_T dtloc, weight, TT1save, scale, tspecies(0:maxspec)
      REAL_T Ct(maxspec),Qt(maxreac),Ytemp(maxspec),Yres(maxspec)

      character*(maxspnml) name
      parameter (LOUTCK=6)

      logical newJ_triggered, bad_soln

#ifdef BL_USE_OMP
      include "omp_lib.h"
      REAL_T,  pointer, save :: dvrwrk(:,:)
      integer, pointer, save :: dviwrk(:,:)
      integer,          save :: dvrlen
      integer,          save :: dvilen
#endif

      REAL_T YJ_SAVE(80)
      LOGICAL FIRST
      COMMON /VHACK/ YJ_SAVE, FIRST
      SAVE   /VHACK/

!$omp threadprivate(/VHACK/,/DVOD01/,/DVOD02/)

c      
c     Set IOPT=1 parameter settings for VODE
c
c     TODO -- these really only need to be set once.
c
      RWRK(dvbr+4) = 0
      RWRK(dvbr+5) = 0
      RWRK(dvbr+6) = ten2minus19
      IWRK(dvbi+4) = 0
      IWRK(dvbi+5) = max_vode_subcycles
      IWRK(dvbi+6) = 0
      
      if (do_diag.eq.1) nsubchem = nchemdiag
c
c     Set molecular weights and pressure in area accessible by conpF
c
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)

      RWRK(NP) = Patm * P1atm

      if (nstiff .eq. 1) then
         MF = 22  ! finite difference jacobian
      else
         MF = 10
      endif
c      
c     Set up ATOL
c
      if (ITOL.eq.2) then
         ATOL(1) = spec_scalT*ATOLEPS
         do m=1,Nspec
            ATOL(m+1) = ATOLEPS*spec_scalY(m)
         end do
      else
         ATOL(1) = ATOLEPS
      endif

      if (do_diag.eq.1) then
         nsub  = nsubchem
         dtloc = dt/nsubchem
      else
         nsub  = 1
         dtloc = dt
      endif

!$omp parallel private(nthrds) copyin(/VHACK/)

#ifdef BL_USE_OMP

!$omp single
      if (.not.associated(dvrwrk)) then
         dvrlen = 22 + 9*(Nspec+1) + 2*(Nspec+1)**2
         dvilen = 30 + Nspec + 1
         nthrds = omp_get_num_threads()

         allocate(dvrwrk(dvrlen,nthrds),dviwrk(dvilen,nthrds))

         do i = 1,nthrds
            dvrwrk(5,i) = 0
            dvrwrk(6,i) = 0
            dvrwrk(7,i) = ten2minus19
            dviwrk(5,i) = 0
            dviwrk(6,i) = max_vode_subcycles
            dviwrk(7,i) = 0
         end do
      endif
!$omp end single

#define DVRWK   dvrwrk(1,tid)
#define DVIWK   dviwrk(1,tid)
#define DVIFCNT dviwrk(12,tid)
#define DVRLEN  dvrlen
#define DVILEN  dvilen

#else

#define DVRWK   RWRK(dvbr)
#define DVIWK   IWRK(dvbi)
#define DVIFCNT IWRK(dvbi+11)
#define DVRLEN  dvr
#define DVILEN  dvi

#endif /*BL_USE_OMP*/

!$omp do private(i,j,k,m,TT1,ISTATE,sum,atoln,Ytemp,Yres,newJ_triggered)
!$omp&private(scale,bad_soln,Ct,Qt,node,weight,TT1save,TT2,tid,tspecies)
!$omp&private(lout,name)
      do k=lo(3),hi(3)

#ifdef BL_USE_OMP
         tid = omp_get_thread_num() + 1
#endif
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               TT1    = zero
               ISTATE = 1

#ifdef SOLN_IS_1D
            if (i.ne.lo(1) .and. j.ne.lo(2)) then
               do m=1,Nspec
                  Ynew(i,j,k,m) = Ynew(lo(1),lo(2),k,m)
               end do
               Tnew(i,j,k) = Tnew(lo(1),lo(2),k)
            else
#endif
               tspecies(0) = Told(i,j,k)

#ifdef DO_JBB_HACK
#ifdef DO_JBB_HACK_TEMP

               tspecies(0) = MAX(HACK_TEMP_MIN, MIN(HACK_TEMP_MAX, tspecies(0)))
#endif
            sum = zero
            do m=1,Nspec
               Ytemp(m) = Yold(i,j,k,m)
               atoln = ATOLEPS
               if (ITOL.eq.2) atoln = atol(m+1)
               Ytemp(m) = MAX(Yold(i,j,k,m),zero)
               sum = sum+Ytemp(m)
            end do
            if (iN2 .gt. 0) Ytemp(iN2) = Ytemp(iN2)+one-sum
#else
            do m=1,Nspec
               Ytemp(m) = Yold(i,j,k,m)
            end do
#endif
            do m = 1,Nspec
               tspecies(m) = Ytemp(m)
            end do

#ifdef TRIGGER_NEW_J
            newJ_triggered = .FALSE.
            sum = zero
            do m=1,NEQ
               scale = spec_scalT
               if (m.ne.1) scale = spec_scalY(m-1)
               sum = sum + ABS(tspecies(m-1)-YJ_SAVE(m))/scale
            end do
            if (sum .gt. NEWJ_TOL) then
               FIRST = .TRUE.
               newJ_triggered = .TRUE.
            end if
#endif

#ifdef ALWAYS_NEW_J
            FIRST = .TRUE.
#endif
            if (do_diag.eq.1) then
               FuncCount(i,j,k) = 0
               CALL CKYTCP(RWRK(NP),tspecies(0),tspecies(1),IWRK(ckbi),RWRK(ckbr),Ct)
               CALL CKQC(tspecies(0),Ct,IWRK(ckbi),RWRK(ckbr),Qt)
               do m=1,Nreac
                  diag(i,j,k,m) = diag(i,j,k,m)+half*dtloc*Qt(m)*million
               enddo
            endif

            do node = 1,nsub
               if (node.lt.nsub) then
                  weight = one
               else
                  weight = half
               endif

               TT1save = TT1
               TT2     = TT1 + dtloc

#if !defined(BL_USE_DOUBLE) || defined(BL_T3E)
               CALL SVODE
#else 
               CALL DVODE
#endif
     &              (CONPF_FILE, NEQ, tspecies, TT1, TT2, ITOL, RTOL, ATOL,
     &              ITASK, ISTATE, IOPT, DVRWK, DVRLEN, DVIWK, DVILEN,
     &              CONPJ_FILE, MF, RWRK, IWRK)
c
c     If the step was bad, and we reused an old Jacobian, try again from scratch
c               
#if defined(TRIGGER_NEW_J) && defined(DO_JBB_HACK)
            if (newJ_triggered .EQV. .FALSE.) then 
               bad_soln = .FALSE.
               do m=1,Nspec
                  if (tspecies(m) .lt. -1.D-6*spec_scalY(m))
     &                 bad_soln = .TRUE.
               end do
               if (bad_soln .EQV. .TRUE.) then
                  FIRST = .TRUE. 
                  tspecies(0) = Told(i,j,k)
                  do m=1,Nspec
                     tspecies(m) = Ytemp(m)
                  end do
                  TT1    = zero
                  ISTATE = 1
#if !defined(BL_USE_DOUBLE) || defined(BL_T3E)
                  CALL SVODE
#else 
                  CALL DVODE
#endif
     &                 (CONPF_FILE, NEQ, tspecies, TT1, TT2, ITOL, RTOL, ATOL,
     &                 ITASK, ISTATE, IOPT, DVRWK, DVRLEN, DVIWK, DVILEN,
     &                 CONPJ_FILE, MF, RWRK, IWRK)
               end if
            end if
#endif
               TT1 = TT2

               if (do_diag.eq.1) then
                  CALL CKYTCP(RWRK(NP),tspecies(0),tspecies(1),IWRK(ckbi),RWRK(ckbr),Ct)
                  CALL CKQC(tspecies(0),Ct,IWRK(ckbi),RWRK(ckbr),Qt)
                  do m=1,Nreac
                     diag(i,j,k,m) = diag(i,j,k,m)+weight*dtloc*Qt(m)*million
                  enddo
                  FuncCount(i,j,k) = FuncCount(i,j,k) + DVIFCNT
               else
                  FuncCount(i,j,k) = DVIFCNT
               endif

               if (verbose_vode .eq. 1) then
#ifdef BL_USE_OMP
                  write(6,*) '......dvode done:'
                  write(6,*) ' last successful step size = ',dvrwrk(11,tid)
                  write(6,*) '          next step to try = ',dvrwrk(12,tid)
                  write(6,*) '   integrated time reached = ',dvrwrk(13,tid)
                  write(6,*) '      number of time steps = ',dviwrk(11,tid)
                  write(6,*) '              number of fs = ',dviwrk(12,tid)
                  write(6,*) '              number of Js = ',dviwrk(13,tid)
                  write(6,*) '    method order last used = ',dviwrk(14,tid)
                  write(6,*) '   method order to be used = ',dviwrk(15,tid)
                  write(6,*) '            number of LUDs = ',dviwrk(19,tid)
                  write(6,*) ' number of Newton iterations ',dviwrk(20,tid)
                  write(6,*) ' number of Newton failures = ',dviwrk(21,tid)
                  if (ISTATE.eq.-4 .or. ISTATE.eq.-5) then
                     call get_spec_name(name,dviwrk(16,tid))
                     write(6,*) '   spec with largest error = ',name
                  end if
#else
                  write(6,*) '......dvode done:'
                  write(6,*) ' last successful step size = ',RWRK(dvbr+10)
                  write(6,*) '          next step to try = ',RWRK(dvbr+11)
                  write(6,*) '   integrated time reached = ',RWRK(dvbr+12)
                  write(6,*) '      number of time steps = ',IWRK(dvbi+10)
                  write(6,*) '              number of fs = ',IWRK(dvbi+11)
                  write(6,*) '              number of Js = ',IWRK(dvbi+12)
                  write(6,*) '    method order last used = ',IWRK(dvbi+13)
                  write(6,*) '   method order to be used = ',IWRK(dvbi+14)
                  write(6,*) '            number of LUDs = ',IWRK(dvbi+18)
                  write(6,*) ' number of Newton iterations ',IWRK(dvbi+19)
                  write(6,*) ' number of Newton failures = ',IWRK(dvbi+20)
                  if (ISTATE.eq.-4 .or. ISTATE.eq.-5) then
                     call get_spec_name(name,IWRK(dvbi+15))
                     write(6,*) '   spec with largest error = ',name
                  end if
#endif
               end if

               if (ISTATE .LE. -1) then
#ifdef BL_USE_OMP
                  call CONPF_FILE(NEQ, TT1, tspecies, Qt, RWRK, IWRK)
                  lout = open_vode_failure_file()                  
                  write(lout,*)
                  write(lout,995) 'VODE Failed at (i,j,k) = (',i,',',j,
     &                 ',',k,'),   Return code = ',ISTATE
                  write(lout,996) 'time(T2,Tl,dt)  ',dt, TT1, dt-TT1
                write(lout,995) 'State ID, old, last, dY/dt, dY/dt*(dt)'
                  write(lout,996) 'T               ',
     &                 Told(i,j,k),tspecies(0),Qt(1),Qt(1)*(dt-TT1)
                  do m=1,Nspec
                     call get_spec_name(name,m)
                     write(lout,996) name,Yold(i,j,k,m),
     &                    tspecies(m),Qt(1+m),Qt(1+m)*(dt-TT1)
                  end do
#else
                  call CONPF_FILE(NEQ, TT1, tspecies, RWRK(NWDOT), RWRK, IWRK)
                  lout = open_vode_failure_file()                  
                  write(lout,*)
                  write(lout,995) 'VODE Failed at (i,j,k) = (',i,',',j,
     &                 ',',k,'),   Return code = ',ISTATE
                  write(lout,996) 'time(T2,Tl,dt)  ',dt, TT1, dt-TT1
                write(lout,995) 'State ID, old, last, dY/dt, dY/dt*(dt)'
                  write(lout,996) 'T               ',
     &                 Told(i,j,k),tspecies(0),RWRK(NWDOT),RWRK(NWDOT)*(dt-TT1)
                  do m=1,Nspec
                     call get_spec_name(name,m)
                     write(lout,996) name,Yold(i,j,k,m),
     &                    tspecies(m),RWRK(NWDOT+m),RWRK(NWDOT+m)*(dt-TT1)
                  end do
#endif
 995              format(a,4(i4,a))
 996              format(a16,1x,4e30.22)
                  close(lout)
                  call bl_abort('VODE failed, see drop file, exiting')
               end if
             enddo

#ifdef DO_JBB_HACK_POST

               tspecies(0) = MAX(HACK_TEMP_MIN, MIN(HACK_TEMP_MAX, tspecies(0)))
#endif
               Tnew(i,j,k) = tspecies(0)

               do m= 1,Nspec
                  Yres(m) = tspecies(m)
               end do

#if defined(DO_JBB_HACK)
               do m=1,Nspec
                  Ynew(i,j,k,m) = Yold(i,j,k,m)+Yres(m)-Ytemp(m)
               end do
#elif  defined(DO_JBB_HACK_POST)
            sum = zero
            do m=1,Nspec
               Ytemp(m) = Yres(i,j,k,m)
               atoln = ATOLEPS
               if (ITOL.eq.2) atoln = atol(m+1)
               Ytemp(m) = MAX(Yres(i,j,k,m),zero)
               sum = sum+Ytemp(m)
            end do
            if (iN2 .gt. 0) then
               Ytemp(iN2) = Ytemp(iN2)+one-sum
            endif
#else
               do m=1,Nspec
                  Ynew(i,j,k,m) = Yres(m)
               end do
#endif

#ifdef SOLN_IS_1D
               endif
#endif
            end do
         end do
      end do
!$omp end do

!$omp end parallel


#undef DVSRWK
#undef DVSIWK
#undef DVRLEN
#undef DVILEN

      end

      subroutine FORT_CHEMEQ(lo, hi,
     &     Ynew, DIMS(Ynew), 
     &     Tnew, DIMS(Tnew),
     &     Yold, DIMS(Yold), 
     &     Told, DIMS(Told),
     &     FuncCount, DIMS(FuncCount),
     &     Patm,
     &     dt)
      implicit none
#include "cdwrk.H"
#include "conp.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(Yold)
      integer DIMDEC(Told)
      integer DIMDEC(Ynew)
      integer DIMDEC(Tnew)
      integer DIMDEC(FuncCount)
      REAL_T Yold(DIMV(Yold),*)
      REAL_T Told(DIMV(Told))
      REAL_T Ynew(DIMV(Ynew),*)
      REAL_T Tnew(DIMV(Tnew))
      REAL_T FuncCount(DIMV(FuncCount))
      REAL_T Patm, dt

      integer i,j,k,n
      REAL_T Yin(maxspec), Temp, Yout(maxspec), sum, mwt(maxspec)
      REAL_T RU, RUC, P1atm, Perg, Tin, rho
      integer gc,rc,tgc,trc

c     Set molecular weights and pressure in area accessible by conpF
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      Perg = Patm * P1atm

      call CKWT(IWRK(ckbi),RWRK(ckbr),mwt)
      
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

#if 0
               rho = zero
               do n=1,Nspec
                  rho = rho + Yold(i,j,k,n)
               enddo

               sum = 0.d0
               do n=1,Nspec
                  Yin(n) = MAX(Yold(i,j,k,n)/rho,zero)
                  Yin(n) = Yold(i,j,k,n)/rho
                  sum = sum + Yin(n)/mwt(n)
               end do
               Tin = Perg / (RU * rho*1.d-3 * sum)

               if (ABS(Tin-Told(i,j,k)).gt.20.d0) then
                  write(6,100)'********',i,j,Tin,Told(i,j,k)
 100              format(a,2i4,2e20.10)
               endif

               call chemeqpt(Perg,Tin,Tnew(i,j,k),Yin,Yout,dt)
               call chemct(1.0,gc,rc,tgc,trc)
               FuncCount(i,j,k) = gc
              
               do n=1,Nspec
                  Ynew(i,j,k,n) = rho*(Yold(i,j,k,n)+Yout(n)-Yin(n))
               end do
#else
               sum = zero
               do n=1,Nspec
                  Yin(n) = MAX(Yold(i,j,k,n),zero)
                  sum = sum + Yin(n)
               end do
               if (iN2 .gt. 0) then
                  Yin(iN2) = Yin(iN2)+one-sum
               endif

               call chemeqpt(Perg,Told(i,j,k),Tnew(i,j,k),Yin,Yout,dt)
               call chemct(1.0,gc,rc,tgc,trc)
               FuncCount(i,j,k) = gc

               do n=1,Nspec
                  Ynew(i,j,k,n) = Yold(i,j,k,n)+Yout(n)-Yin(n)
               end do
#endif
            end do
         end do
      end do
      
      end
 
      subroutine FORT_MIXAVG_RHODIFF_TEMP(lo, hi, RD, DIMS(RD), T,
     &     DIMS(T), Y, DIMS(Y), Patm, do_temp, do_VelVisc)
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM), do_temp, do_VelVisc
      integer DIMDEC(RD)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T RD(DIMV(RD),*)
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T Patm

      integer i, j, k, n
      REAL_T RU, RUC, P1ATM, Ptmp, Yt(maxspec), Dt(maxspec)
      REAL_T SCAL, TSCAL, Wavg, RHO, Tt, invmwt(maxspec)
      REAL_T alpha, l1, l2, X(maxspec), CPMS(maxspec)

#ifdef BL_USE_OMP
      include "omp_lib.h"
      REAL_T,  pointer :: egrwrk(:,:)
      integer, pointer :: egiwrk(:,:)
      common /egwork/ egrwrk,egiwrk
      save /egwork/
#endif
      integer egrlen, egilen, tid, nthrds

      parameter(SCAL = tenth, TSCAL = one / 100000.0D0)

      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM
      call CKWT(IWRK(ckbi),RWRK(ckbr),invmwt)

      do n=1,Nspec
         invmwt(n) = one / invmwt(n)
      end do

!$omp parallel private(egrlen,egilen,nthrds)

#ifdef BL_USE_OMP

!$omp single
      if (.not.associated(egrwrk)) then
         egrlen = 23 + 14*Nspec + 32*Nspec**2 + 13*eg_nodes
     &        + 30*eg_nodes*Nspec + 5*eg_nodes*Nspec**2
         egilen = Nspec
         nthrds = omp_get_num_threads()

         allocate(egrwrk(egrlen,nthrds),egiwrk(egilen,nthrds))

         do i = 1,nthrds
            do n=1,egrlen
               egrwrk(n,i) = RWRK(egbr+n-1)
            enddo
         enddo
         do i = 1,nthrds
            do n=1,egilen
               egiwrk(n,i) = IWRK(egbi+n-1)
            enddo
         end do
      endif
!$omp end single

#define EGSRWK egrwrk(1,tid)
#define EGSIWK egiwrk(1,tid)

#else

#define EGSRWK RWRK(egbr)
#define EGSIWK IWRK(egbi)

#endif /*BL_USE_OMP*/

!$omp do private(i,j,k,n,Yt,Tt,alpha,Dt,Wavg)
!$omp&private(CPMS,X,RHO,l1,l2,tid)

      do k=lo(3),hi(3)

#ifdef BL_USE_OMP
         tid = omp_get_thread_num() + 1
#endif
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               Tt = MAX(T(i,j,k),TMIN_TRANS)
               CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
               CALL CKCPMS(Tt,IWRK(ckbi),RWRK(ckbr),CPMS)
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
               CALL EGSPAR(Tt,X,Yt,CPMS,EGSRWK,EGSIWK)
               CALL EGSV1(Ptmp,Tt,Yt,Wavg,EGSRWK,Dt)
               CALL CKRHOY(Ptmp,Tt,Yt,IWRK(ckbi),RWRK(ckbr),RHO)
               do n=1,Nspec
                  RD(i,j,k,n) = RHO * Wavg * invmwt(n) * Dt(n) * SCAL
               end do

               if (do_temp .ne. 0) then
                  alpha = 1
                  CALL EGSL1(alpha,Tt,X,EGSRWK,l1)
                  alpha = -1
                  CALL EGSL1(alpha,Tt,X,EGSRWK,l2)
                  RD(i,j,k,Nspec+1) = half * (l1 + l2) * TSCAL
               endif

               if (do_VelVisc .ne. 0) then
                  CALL EGSE3(Tt,Yt,EGSRWK,RD(i,j,k,Nspec+2))
                  RD(i,j,k,Nspec+2) = RD(i,j,k,Nspec+2) * SCAL
               endif

            end do
         end do
      end do
!$omp end do

!$omp end parallel

#undef EGSRWK
#undef EGSIWK
      end

      subroutine FORT_MIX_SHEAR_VISC(lo, hi, eta, DIMS(eta),
     &                               T, DIMS(T), Y, DIMS(Y))
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(eta)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T eta(DIMV(eta))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T SCAL
      
      integer i, j, k, n
      REAL_T X(maxspec), Yt(maxspec), CPMS(maxspec), Tt

#ifdef BL_USE_OMP
      include "omp_lib.h"
      REAL_T,  pointer :: egrwrk(:,:)
      integer, pointer :: egiwrk(:,:)
      common /egwork/ egrwrk,egiwrk
      save /egwork/
#endif
      integer egrlen, egilen, tid, nthrds

      parameter(SCAL = tenth)

!$omp parallel private(egrlen,egilen,nthrds)

#ifdef BL_USE_OMP

!$omp single

      if (.not.associated(egrwrk)) then
         egrlen = 23 + 14*Nspec + 32*Nspec**2 + 13*eg_nodes
     &        + 30*eg_nodes*Nspec + 5*eg_nodes*Nspec**2
         egilen = Nspec
         nthrds = omp_get_num_threads()

         allocate(egrwrk(egrlen,nthrds),egiwrk(egilen,nthrds))

         do i = 1,nthrds
            do n=1,egrlen
               egrwrk(n,i) = RWRK(egbr+n-1)
            enddo
         enddo
         do i = 1,nthrds
            do n=1,egilen
               egiwrk(n,i) = IWRK(egbi+n-1)
            enddo
         end do
      endif

!$omp end single

#define EGSRWK egrwrk(1,tid)
#define EGSIWK egiwrk(1,tid)

#else

#define EGSRWK RWRK(egbr)
#define EGSIWK IWRK(egbi)

#endif /*BL_USE_OMP*/
!
! The following computes the mixture averaged shear viscosity using EGLib
! Note that SCAL converts assumed cgs units to MKS (1 g/cm.s = .1 kg/m.s)
!
!$omp do private(i,j,k,n,Yt,Tt,CPMS,X,tid)
      do k=lo(3),hi(3)

#ifdef BL_USE_OMP
         tid = omp_get_thread_num() + 1
#endif
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               Tt = MAX(T(i,j,k),TMIN_TRANS) 
               CALL CKCPMS(Tt,IWRK(ckbi),RWRK(ckbr),CPMS)
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
               CALL EGSPAR(Tt,X,Yt,CPMS,EGSRWK,EGSIWK)
               CALL EGSE3(Tt,Yt,EGSRWK,eta(i,j,k))
               eta(i,j,k) = eta(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end do

!$omp end parallel

#undef EGSRWK
#undef EGSIWK
      end

      subroutine FORT_RHOfromPTY(lo, hi, RHO, DIMS(RHO), T, DIMS(T),
     &                           Y, DIMS(Y), Patm)
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T Patm
      
      integer i, j, k, n
      REAL_T RU, RUC, P1ATM, Ptmp, Yt(maxspec), SCAL
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
c
      parameter(SCAL = one * 1000.0D0)
      
      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM

!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKRHOY(Ptmp,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),RHO(i,j,k))
               RHO(i,j,k) = RHO(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do

      end
      
      subroutine FORT_RHOfromPvTY(lo, hi, RHO, DIMS(RHO), T, DIMS(T),
     &                           Y, DIMS(Y), P, DIMS(P))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer DIMDEC(P)
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T P(DIMV(P))
      
      integer i, j, k, n
      REAL_T RU, RUC, P1ATM, Ptmp, Yt(maxspec), SCAL
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
c
      parameter(SCAL = one * 1000)
      
      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               Ptmp = P(i,j,k) * P1ATM
               CALL CKRHOY(Ptmp,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),RHO(i,j,k))
               RHO(i,j,k) = RHO(i,j,k) * SCAL
            end do
         end do
      end do
      end
      
      subroutine FORT_PfromRTY(lo, hi, P, DIMS(P), RHO, DIMS(RHO),
     &                         T, DIMS(T), Y, DIMS(Y))
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(P)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T P(DIMV(P))
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(maxspec), RHOt, SCAL, SCAL1
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 dyne/cm^2 = .1 Pa)
c           SCAL1 converts density (1 kg/m^3 = 1.e-3 g/cm^3)
c
      parameter(SCAL = tenth, SCAL1 = tenth**3)

!$omp parallel do private(i,j,k,n,Yt,RHOt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               RHOt = RHO(i,j,k) * SCAL1
               CALL CKPY(RHOt,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),P(i,j,k))
               P(i,j,k) = P(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do

      end
      
      subroutine FORT_TfromPRY(lo, hi, T, DIMS(T), RHO, DIMS(RHO),
     &                         Y, DIMS(Y), Patm)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T Patm
      
      integer i, j, k, n
      REAL_T RU, RUC, P1ATM, Ptmp, Yt(maxspec), SCAL, Wavg, RHOt
c
c     NOTE: SCAL converts density (1 kg/m^3 = 1.e-3 g/cm^3)
c
      parameter(SCAL = tenth**3)
      
      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
               RHOt = RHO(i,j,k) * SCAL
               T(i,j,k) = Ptmp / (RHOt * RU / Wavg)
            end do
         end do
      end do
      end
      
      subroutine FORT_CPMIXfromTY(lo, hi, CPMIX, DIMS(CPMIX), T, DIMS(T),
     &                            Y, DIMS(Y))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(CPMIX)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T CPMIX(DIMV(CPMIX))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(maxspec), SCAL
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
c
      parameter(SCAL = tenth**4)

!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKCPBS(T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),CPMIX(i,j,k))
               CPMIX(i,j,k) = CPMIX(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do
      end
      
      subroutine FORT_CVMIXfromTY(lo, hi, CVMIX, DIMS(CVMIX), T, DIMS(T),
     &                            Y, DIMS(Y))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(CVMIX)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T CVMIX(DIMV(CVMIX))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(maxspec), SCAL
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
c
      parameter(SCAL = tenth**4)

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKCVBS(T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),CVMIX(i,j,k))
               CVMIX(i,j,k) = CVMIX(i,j,k) * SCAL
            end do
         end do
      end do
      end
      
      subroutine FORT_HMIXfromTY(lo, hi, HMIX, DIMS(HMIX), T, DIMS(T),
     &                           Y, DIMS(Y))
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(HMIX)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T HMIX(DIMV(HMIX))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(maxspec), SCAL
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
c
      parameter(SCAL = tenth**4)

!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKHBMS(T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),HMIX(i,j,k))
               HMIX(i,j,k) = HMIX(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do

      end
      
      subroutine FORT_MWMIXfromY(lo, hi, MWMIX, DIMS(MWMIX), Y, DIMS(Y))
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(MWMIX)
      integer DIMDEC(Y)
      REAL_T MWMIX(DIMV(MWMIX))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k,n
      REAL_T Yt(maxspec)
!
!     Returns mean molecular weight in kg/kmole
!
!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),MWMIX(i,j,k))
            end do
         end do
      end do
!$omp end parallel do

      end
      
      subroutine FORT_CPfromT(lo, hi, CP, DIMS(CP), T, DIMS(T))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(CP)
      integer DIMDEC(T)
      REAL_T CP(DIMV(CP),*)
      REAL_T T(DIMV(T))
      
      integer i, j, k, n
      REAL_T SCAL, CPt(maxspec)
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
c
      parameter(SCAL = tenth**4)
      
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               CALL CKCPMS(T(i,j,k),IWRK(ckbi),RWRK(ckbr),CPt)
               do n=1,Nspec
                  CP(i,j,k,n) = CPt(n) * SCAL
               end do
            end do
         end do
      end do
      end
      
      subroutine FORT_HfromT(lo, hi, H, DIMS(H), T, DIMS(T))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(H)
      integer DIMDEC(T)
      REAL_T H(DIMV(H),*)
      REAL_T T(DIMV(T))
      
      integer i, j, k, n
      REAL_T SCAL, Ht(maxspec)
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
c
      parameter(SCAL = tenth**4)

!$omp parallel do private(i,j,k,n,Ht)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               CALL CKHMS(T(i,j,k),IWRK(ckbi),RWRK(ckbr),Ht)
               do n=1,Nspec
                  H(i,j,k,n) = Ht(n) * SCAL
               end do
            end do
         end do
      end do
!$omp end parallel do

      end

      integer function FORT_TfromHY(lo, hi, T, DIMS(T),
     &                              HMIX, DIMS(HMIX), Y, DIMS(Y),
     &                              errMax, NiterMAX, res)
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer NiterMAX
      integer DIMDEC(T)
      integer DIMDEC(HMIX)
      integer DIMDEC(Y)
      REAL_T T(DIMV(T))
      REAL_T HMIX(DIMV(HMIX))
      REAL_T Y(DIMV(Y),*)
      REAL_T errMAX
      REAL_T res(0:NiterMAX-1)
      REAL_T Yt(maxspec)
      integer i,j,k,n,FORT_TfromHYpt,Niter,MAXiters

      MAXiters = 0

!$omp parallel do private(i,j,k,n,Yt,res,Niter) reduction(max:MAXiters)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do

               Niter = FORT_TfromHYpt(T(i,j,k),HMIX(i,j,k),Yt,errMax,NiterMAX,res)

               if (Niter .lt. 0) then
                  write(6,*) 'T from h,y solve in FORT_TfromHY failed',Niter
                  call bl_abort(" ")
               end if
               
               if (Niter .gt. MAXiters) MAXiters = Niter

            end do
         end do
      end do
!$omp end parallel do

c     Set max iters taken during this solve, and exit
      FORT_TfromHY = MAXiters
      return
      end
c
c     Optically thin radiation model, specified at
c            http://www.ca.sandia.gov/tdf/Workshop/Submodels.html
c     
c     Q(T,species) = 4*sigma*SUM{pi*aP,i} *(T4-Tb4) 
c     
c     sigma=5.669e-08 W/m2K4 is the Steffan-Boltzmann constant, 
c     SUM{ } represents a summation over the species in the radiation calculation, 
c     pi is partial pressure of species i in atm (Xi times local pressure)
c     aP,i is the Planck mean absorption coefficient of species i, 1/[m.atm]
c     T is the local flame temperature (K)
c     Tb is the background temperature (300K or as spec. in expt)
c
c     For H2O and CO2,
c         aP = exp{c0 + c1*ln(T) + c2*{ln(T)}2 + c3*{ln(T)}3 + c4*{ln(T)}4} 
c     
c                            H2O                  CO2
c              c0       0.278713E+03         0.96986E+03
c              c1      -0.153240E+03        -0.58838E+03
c              c2       0.321971E+02         0.13289E+03
c              c3      -0.300870E+01        -0.13182E+02
c              c4       0.104055E+00         0.48396E+00
c     For CH4:
c     
c     aP,ch4 = 6.6334 - 0.0035686*T + 1.6682e-08*T2 + 2.5611e-10*T3 - 2.6558e-14*T4
c
c     For CO:   aP,co = c0+T*(c1 + T*(c2 + T*(c3 + T*c4)))
c
c           T <= 750                 else
c      
c         c0   4.7869              10.09       
c         c1  -0.06953             -0.01183    
c         c2   2.95775e-4          4.7753e-6   
c         c3  -4.25732e-7          -5.87209e-10
c         c4   2.02894e-10         -2.5334e-14 
c      
      
      subroutine FORT_OTrad_TDF(lo, hi, Qloss, DIMS(Qloss),
     &                          T, DIMS(T), X, DIMS(X), Patm, T_bg) 
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(Qloss)
      integer DIMDEC(T)
      integer DIMDEC(X)
      REAL_T Qloss(DIMV(Qloss))
      REAL_T T(DIMV(T))
      REAL_T X(DIMV(X),*)
      REAL_T Patm, T_bg
      
      character*(maxspnml) name      
      integer n, i, j, k, iH2O, iCO2, iCH4, iCO
      REAL_T aP, c0, c1, c2, c3, c4
      REAL_T T1,T2,T3,T4,lnT1,lnT2,lnT3,lnT4,Tb4,sigma

      parameter (sigma = 5.669D-08)
      
      iH2O = 0
      iCO2 = 0
      iCH4 = 0
      iCO  = 0
      
      do n = 1,Nspec
         call get_spec_name(name, n)
         if (name .EQ. 'H20') iH2O = n
         if (name .EQ. 'CO2') iCO2 = n
         if (name .EQ. 'CH4') iCH4 = n
         if (name .EQ. 'CO')  iCO  = n
      end do
      
      Tb4 = T_bg**4

!$omp parallel do private(i,j,k,T1,T2,T3,T4,lnT1,lnT2,lnT3,lnT4)
!$omp&private(aP,c0,c1,c2,c3,c4)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               
               T1 = T(i,j,k)
               T2 = T1*T1
               T3 = T2*T1
               T4 = T3*T1
               
               if ( (iH2O.gt.0) .or. (iCO2.gt.0) ) then
                  lnT1 = LOG(T1)
                  lnT2 = lnT1*lnT1
                  lnT3 = lnT2*lnT1
                  lnT4 = lnT3*lnT1
               end if
               
               aP = zero
               
               if ((iH2O.gt.0).and.(X(i,j,k,iH2O).gt.zero)) then
                  aP = aP + X(i,j,k,iH2O)*EXP(
     &                 + 0.278713D+03
     &                 - 0.153240D+03*lnT1
     &                 + 0.321971D+02*lnT2
     &                 - 0.300870D+01*lnT3
     &                 + 0.104055D+00*lnT4 )
               end if
               
               if ((iCO2.gt.0).and.(X(i,j,k,iCO2).gt.zero)) then            
                  aP = aP + X(i,j,k,iCO2)*EXP(
     &                 + 0.96986D+03
     &                 - 0.58838D+03*lnT1
     &                 + 0.13289D+03*lnT2
     &                 - 0.13182D+02*lnT3
     &                 + 0.48396D+00*lnT4 )
               end if

               if ((iCH4.gt.0).and.(X(i,j,k,iCH4).gt.zero)) then
                  aP = aP + X(i,j,k,iCH4)*
     &                 ( 6.6334D0
     &                 - 0.0035686D0 *T1
     &                 + 1.6682D-08*T2
     &                 + 2.5611D-10*T3
     &                 - 2.6558D-14*T4 )         
               end if
               
               if ((iCO.gt.0).and.(X(i,j,k,iCO).gt.zero)) then
                  if ( T1 .le. 750.0D0 ) then
                     c0 =  4.7869D0
                     c1 = -0.06953D0
                     c2 =  2.95775D-4
                     c3 = -4.25732D-7
                     c4 =  2.02894D-10
                  else
                     c0 =  10.09D0
                     c1 = -0.01183D0
                     c2 =  4.7753D-6
                     c3 = -5.87209D-10
                     c4 = -2.5334D-14
                  endif
                  aP = aP + X(i,j,k,iCO)*(c0 + c1*T1 + c2*T2 + c3*T3 + c4*T4)
               end if

               Qloss(i,j,k) = four*sigma*Patm*(T4-Tb4)*aP
               
            end do
         end do
      end do
!$omp end parallel do

      end
