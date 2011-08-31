C     CVS $Revision: 1.1.1.1 $ reposited $Date: 2006/05/26 19:09:33 $

C///////////////////////////////////////////////////////////////////////
C
C     CHEMKIN-III file dasac.f V.3.4 October 1997
C
C     the public-domain (?) subroutines known as DASAC from the
C     University of Wisconsin.  DASAC is a modification of DASSL,
C     for use with SENKIN
C
C
C     v3.4:  97/10/29  E. Meeks
C        1. Removed function D1MACH, to use D1MACH available in
C           mach.f
C        2. Added change blocks around 2 definitions of UROUND
C           to use D1MACH when double, R1MACH when single.
C     v3.3:  97/10/29  Joseph Grcar
C        1. Replaced SIGN by SIGN77.
C     v3.2:
C        1. Corrected error in calls to XERRWV in DDASAC.  The
C           argument list was missing a 0 after IDID variable
C           (two occurences). (F. Rupley)
C        2. Inital CVS version 3.2.
C
C
C///////////////////////////////////////////////////////////////////////
C
C*****precision > single
C      SUBROUTINE SDASAC (RES,NSYS,T,Y,YPRIME,TOUT,
C     *  INFO,ISEN,RTOL,ATOL,IDID,SWORK,LSW,
C     *  RWORK,LRW,IWORK,LIW,RPAR,IPAR,
C     *  JAC,DRES,DFDYP)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      SUBROUTINE DDASAC (RES,NSYS,T,Y,YPRIME,TOUT,
     *  INFO,ISEN,RTOL,ATOL,IDID,SWORK,LSW,
     *  RWORK,LRW,IWORK,LIW,RPAR,IPAR,
     *  JAC,DRES,DFDYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C-------------------------------------------------------------------
C     This code solves a system of differential/algebraic
C     equations of the form f(t,y,yprime;p)=0. where p is a
C     parameter vector of time independent constants.
C     If requested the code will also solve the following system
C     of equations for the parametric sensitivity functions W(t) :
C
C        (df/dyprime)*W'(t)+(df/dy)*W(t)+df/dp=0.
C
C     where   W(t):=dy/dp
C--------------------------------------------------------------------
      SAVE
      LOGICAL DONE
      EXTERNAL RES,JAC,DRES,DFDYP,SIGN77
      DIMENSION Y(*), YPRIME(*), INFO(15),ISEN(5)
      DIMENSION RWORK(*), IWORK(*), SWORK(*), RTOL(*), ATOL(*)
      DIMENSION RPAR(*), IPAR(*)
      COMMON/DDA001/NPD,NTEMP,
     *   LML,LMU,LMXORD,LMTYPE,
     *   LNST,LNSTL,LNRE,LNJE,LETF,LCTF,LIPVT
      COMMON/DDA002/INDEX,NALG,IDFDP,ICALC,NPAR
      DATA LTSTOP,LHMAX,LH,LTN,
     *   LCJ,LCJOLD,LHOLD,LS,LROUND,
     *   LALPHA,LBETA,LGAMMA,
     *   LPSI,LSIGMA,LDELTA
     *   /1,2,3,4,
     *   5,6,7,8,9,
     *   11,17,23,
     *   29,35,41/
      DATA ZERO/0.0/
C
C     First see if parametric sensitivity calculations
C     have been requested by the user.
C
      ICALC=ISEN(1)
      IF(ICALC.EQ.0) THEN
         LDTEM = 1
         LEMAT = 1
         GO TO 5
      ENDIF
C
C     Rename the elements of the isen array. Store them
C     in common block DDA002 for later use. Also set
C     pointers into swork array.
C
      INDEX=ISEN(2)
      NALG=ISEN(3)
      IDFDP=ISEN(4)
      NPAR=ISEN(5)
      LDTEM=1
      LEMAT=NSYS+1
C
C     Compute the total number of equations.
C
      NEQ=NSYS*(NPAR+1)
5     CONTINUE
C
C     If sensitivity analysis has not been requested, set
C     total number of equations equal to the number of
C     state equations.
C
      IF(ICALC.EQ.0)NEQ=NSYS
      IF(INFO(1).NE.0)GO TO 100
C-----------------------------------------------------------------------
C     The following block is executed for the initial call only.
C     It checks the inputs and initializations.
C-----------------------------------------------------------------------
C
C     Compute length of swork and check for illegal input.
C
      IF(ICALC.EQ.0)GO TO 8
      LENSW=NSYS
      IF(INDEX.EQ.2)LENSW=NSYS*(NSYS+1)
      IF(LSW.LT.LENSW)GO TO 7045
8     CONTINUE
C
C     Check info array to make sure all elements
C     are either zero or one.
C
      DO 10 I=2,11
         IF(INFO(I).NE.0.AND.INFO(I).NE.1)GO TO 701
10       CONTINUE
      IF(NEQ.LE.0)GO TO 702
C
C     Set pointers into iwork.
C
      LML=1
      LMU=2
      LMXORD=3
      LMTYPE=4
      LJCALC=5
      LPHASE=6
      LK=7
      LKOLD=8
      LNS=9
      LNSTL=10
      LNST=11
      LNRE=12
      LNJE=13
      LETF=14
      LCTF=15
      LIPVT=21
      LIWM=1
C
C     Check and store maximum integration order.
C
      MXORD=5
      IF(INFO(9).EQ.0)GO TO 20
         MXORD=IWORK(LMXORD)
         IF(MXORD.LT.1.OR.MXORD.GT.5)GO TO 703
20       IWORK(LMXORD)=MXORD
C
C     Compute mtype,lenpd,lenrw.  Check ml and mu
C
      IF(INFO(6).NE.0)GO TO 40
         LENPD=NSYS**2
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
         IF(INFO(5).NE.0)GO TO 30
            IWORK(LMTYPE)=2
            GO TO 60
30          IWORK(LMTYPE)=1
            GO TO 60
40    IF(IWORK(LML).LT.0.OR.IWORK(LML).GE.NSYS)GO TO 717
      IF(IWORK(LMU).LT.0.OR.IWORK(LMU).GE.NSYS)GO TO 718
      LENPD=(2*IWORK(LML)+IWORK(LMU)+1)*NSYS
      IF(INFO(5).NE.0)GO TO 50
         IWORK(LMTYPE)=5
         MBAND=IWORK(LML)+IWORK(LMU)+1
         MSAVE=(NSYS/MBAND)+1
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD+2*MSAVE
         GO TO 60
50       IWORK(LMTYPE)=4
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
C
C     Check lengths of rwork and iwork.
C
60    LENIW=20+NEQ
      IF(LRW.LT.LENRW)GO TO 704
      IF(LIW.LT.LENIW)GO TO 705
C
C     Check to see that tout is different from t.
C
      IF(TOUT .EQ. T)GO TO 719
C
C     Check hmax.
C
      IF(INFO(7).EQ.0)GO TO 70
         HMAX=RWORK(LHMAX)
         IF(HMAX.LE.0.0E0)GO TO 710
70    CONTINUE
C
C     Initialize counters.
C
      IWORK(LNST)=0
      IWORK(LNRE)=0
      IWORK(LNJE)=0
      IWORK(LNSTL)=0
      IDID=1
      GO TO 200
C-----------------------------------------------------------------------
C     This block is for continuation calls
C     only. Here we check info(1),and if the
C     last step was interrupted we check whether
C     appropriate action was taken.
C-----------------------------------------------------------------------
100   CONTINUE
      IF(INFO(1).EQ.1)GO TO 110
      IF(INFO(1).NE.-1)GO TO 701
C
C     If control reaches here, the last step was interrupted
C     by an error condition from ddastp,and
C     appropriate action was not taken. This
C     is a fatal error.
C
      CALL XERRWV('DASSL--  THE LAST STEP TERMINATED WITH A NEGATIVE',
     1             49,201,0,0,0,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  VALUE (=I1) OF IDID AND NO APPROPRIATE',
     1             47,202,0,1,IDID,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  ACTION WAS TAKEN. RUN TERMINATED',
     1             41,203,1,0,0,0,0,ZERO,ZERO)
      RETURN
110   CONTINUE
      IWORK(LNSTL)=IWORK(LNST)
C-----------------------------------------------------------------------
C     This block is executed on all calls.
C     The error tolerance parameters are
C     checked, and the work array pointers
C     are set.
C-----------------------------------------------------------------------
200   CONTINUE
      NZFLG=0
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 210 I=1,NEQ
         IF(INFO(2).EQ.1)RTOLI=RTOL(I)
         IF(INFO(2).EQ.1)ATOLI=ATOL(I)
         IF(RTOLI.GT.0.0E0.OR.ATOLI.GT.0.0E0)NZFLG=1
         IF(RTOLI.LT.0.0E0)GO TO 706
         IF(ATOLI.LT.0.0E0)GO TO 707
210      CONTINUE
      IF(NZFLG.EQ.0)GO TO 708
C
C     Set up rwork storage. The iwork storage is fixed
C     in a data statement.
C
      LE=LDELTA+NEQ
      LWT=LE+NEQ
      LPHI=LWT+NEQ
      LPD=LPHI+(IWORK(LMXORD)+1)*NEQ
      LWM=LPD
      NPD=1
      NTEMP=NPD+LENPD
      IF(INFO(1).EQ.1)GO TO 400
C-----------------------------------------------------------------------
C     This block is executed on the initial call
C     only. Set the initial step size,
C     the error weight vector, and phi.
C     Compute initial yprime, if necessary.
C-----------------------------------------------------------------------
300   CONTINUE
      TN=T
      IDID=1
C
C     Set error weight vector wt.
C
      CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      DO 305 I = 1,NEQ
         IF(RWORK(LWT+I-1).LE.0.0E0) GO TO 713
305      CONTINUE
C
C     Compute unit roundoff and hmin.
C
C*****precision > single
C      UROUND = R1MACH(4)
C*****END precision > single
C*****precision > double
      UROUND = D1MACH(4)
C*****END precision > double
      RWORK(LROUND) = UROUND
      HMIN = 4.0E0*UROUND*MAX(ABS(T),ABS(TOUT))
C
C     Check initial interval to see that it is long enough.
C
      TDIST = ABS(TOUT - T)
      IF(TDIST .LT. HMIN) GO TO 714
C
C     Check ho, if this was input.
C
      IF (INFO(8) .EQ. 0) GO TO 310
         HO = RWORK(LH)
         IF ((TOUT - T)*HO .LT. 0.0E0) GO TO 711
         IF (HO .EQ. 0.0E0) GO TO 712
         GO TO 320
310    CONTINUE
C
C     Compute initial stepsize, to be used by either
C     ddastp or ddaini, depending on info(11).
C
      HO = 0.001E0*TDIST
      YPNORM = DDANRM(NEQ,YPRIME,RWORK(LWT),RPAR,IPAR)
      IF (YPNORM .GT. 0.5E0/HO) HO = 0.5E0/YPNORM
      HO = SIGN77(HO,TOUT-T)
C
C     Adjust ho if necessary to meet hmax bound.
C
320   IF (INFO(7) .EQ. 0) GO TO 330
         RH = ABS(HO)/HMAX
         IF (RH .GT. 1.0E0) HO = HO/RH
C
C     Compute tstop, if applicable.
C
330   IF (INFO(4) .EQ. 0) GO TO 340
         TSTOP = RWORK(LTSTOP)
         IF ((TSTOP - T)*HO .LT. 0.0E0) GO TO 715
         IF ((T + HO - TSTOP)*HO .GT. 0.0E0) HO = TSTOP - T
         IF ((TSTOP - TOUT)*HO .LT. 0.0E0) GO TO 709
C
C     Compute initial derivative, if applicable.
C
340   IF (INFO(11) .EQ. 0) GO TO 350
      CALL DDAINI(T,Y,YPRIME,NSYS,
     *  RES,JAC,HO,RWORK(LWT),IDID,RPAR,IPAR,
     *  RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
     *  RWORK(LWM),IWORK(LIWM),HMIN,RWORK(LROUND),
     *  SWORK(LDTEM), INFO(10))
      IF (IDID .LT. 0) GO TO 390
C
C     Load h with ho.  Store h in rwork(lh).
C
350   H = HO
      RWORK(LH) = H
C
C     Load y and h*yprime into phi(*,1) and phi(*,2).
C
360   ITEMP = LPHI + NEQ
      DO 370 I = 1,NEQ
         RWORK(LPHI + I - 1) = Y(I)
370      RWORK(ITEMP + I - 1) = H*YPRIME(I)
390   GO TO 500
C-------------------------------------------------------
C     This block is for continuation calls only. Its
C     purpose is to check stop conditions before
C     taking a step.
C     Adjust h if necessary to meet hmax bound.
C-------------------------------------------------------
400   CONTINUE
      DONE = .FALSE.
      TN=RWORK(LTN)
      H=RWORK(LH)
      IF(INFO(7) .EQ. 0) GO TO 410
         RH = ABS(H)/HMAX
         IF(RH .GT. 1.0E0) H = H/RH
410   CONTINUE
      IF(T .EQ. TOUT) GO TO 719
      IF((T - TOUT)*H .GT. 0.0E0) GO TO 711
      IF(INFO(4) .EQ. 1) GO TO 430
      IF(INFO(3) .EQ. 1) GO TO 420
      IF((TN-TOUT)*H.LT.0.0E0)GO TO 490
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
420   IF((TN-T)*H .LE. 0.0E0) GO TO 490
      IF((TN - TOUT)*H .GT. 0.0E0) GO TO 425
      CALL DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
425   CONTINUE
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
430   IF(INFO(3) .EQ. 1) GO TO 440
      TSTOP=RWORK(LTSTOP)
      IF((TN-TSTOP)*H.GT.0.0E0) GO TO 715
      IF((TSTOP-TOUT)*H.LT.0.0E0)GO TO 709
      IF((TN-TOUT)*H.LT.0.0E0)GO TO 450
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *   RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
440   TSTOP = RWORK(LTSTOP)
      IF((TN-TSTOP)*H .GT. 0.0E0) GO TO 715
      IF((TSTOP-TOUT)*H .LT. 0.0E0) GO TO 709
      IF((TN-T)*H .LE. 0.0E0) GO TO 450
      IF((TN - TOUT)*H .GT. 0.0E0) GO TO 445
      CALL DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
445   CONTINUE
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
450   CONTINUE
C
C     Check whether we are with in roundoff of tstop.
C
      IF(ABS(TN-TSTOP).GT.100.0E0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 460
      IDID=2
      T=TSTOP
      DONE = .TRUE.
      GO TO 490
460   TNEXT=TN+H*(1.0E0+4.0E0*UROUND)
      IF((TNEXT-TSTOP)*H.LE.0.0E0)GO TO 490
      H=(TSTOP-TN)*(1.0E0-4.0E0*UROUND)
      RWORK(LH)=H
490   IF (DONE) GO TO 590
C-------------------------------------------------------
C     The next block contains the call to the
C     one-step integrator ddastp.
C     This is a looping point for the integration
C     steps.
C     Check for too many steps.
C     Update wt.
C     Check for too much accuracy requested.
C     Compute minimum stepsize.
C-------------------------------------------------------
500   CONTINUE
C
C     Check for failure to compute initial yprime.
C
      IF (IDID .EQ. -12) GO TO 527
C
C     Check for too many steps.
C
      IF((IWORK(LNST)-IWORK(LNSTL)).LT.500)
     *   GO TO 510
           IDID=-1
           GO TO 527
C
C     Update wt.
C
510   CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,RWORK(LPHI),
     *  RWORK(LWT),RPAR,IPAR)
      DO 520 I=1,NEQ
         IF(RWORK(I+LWT-1).GT.0.0E0)GO TO 520
           IDID=-3
           GO TO 527
520   CONTINUE
C
C     Test for too much accuracy requested.
C
      R=DDANRM(NEQ,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)*
     *   100.0E0*UROUND
      IF(R.LE.1.0E0)GO TO 525
C
C     Multiply rtol and atol by r and return.
C
      IF(INFO(2).EQ.1)GO TO 523
           RTOL(1)=R*RTOL(1)
           ATOL(1)=R*ATOL(1)
           IDID=-2
           GO TO 527
523   DO 524 I=1,NEQ
           RTOL(I)=R*RTOL(I)
524        ATOL(I)=R*ATOL(I)
      IDID=-2
      GO TO 527
525   CONTINUE
C
C     Compute minimum stepsize.
C
      HMIN=4.0E0*UROUND*MAX(ABS(TN),ABS(TOUT))
      CALL DDASTP(TN,Y,YPRIME,NEQ,NSYS,
     *   RES,JAC,H,RWORK(LWT),INFO(1),IDID,RPAR,IPAR,
     *   RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
     *   RWORK(LWM),IWORK(LIWM),
     *   RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),
     *   RWORK(LPSI),RWORK(LSIGMA),
     *   RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),
     *   RWORK(LS),HMIN,RWORK(LROUND),
     *   IWORK(LPHASE),IWORK(LJCALC),IWORK(LK),
     *   IWORK(LKOLD),IWORK(LNS),INFO(10),
     *   SWORK(LDTEM),SWORK(LEMAT),DRES,DFDYP)
527   IF(IDID.LT.0)GO TO 600
C------------------------------------------------------
C     This block handles the case of a successful
C     return from ddastp (idid=1) test for
C     stop conditions.
C--------------------------------------------------------
C
      IF(INFO(4).NE.0)GO TO 540
           IF(INFO(3).NE.0)GO TO 530
             IF((TN-TOUT)*H.LT.0.0E0)GO TO 500
             CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
530          IF((TN-TOUT)*H.GE.0.0E0)GO TO 535
             T=TN
             IDID=1
             GO TO 580
535          CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
540   IF(INFO(3).NE.0)GO TO 550
      IF((TN-TOUT)*H.LT.0.0E0)GO TO 542
         CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
         T=TOUT
         IDID=3
         GO TO 580
542   IF(ABS(TN-TSTOP).LE.100.0E0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 545
      TNEXT=TN+H*(1.0E0+4.0E0*UROUND)
      IF((TNEXT-TSTOP)*H.LE.0.0E0)GO TO 500
      H=(TSTOP-TN)*(1.0E0-4.0E0*UROUND)
      GO TO 500
545   IDID=2
      T=TSTOP
      GO TO 580
550   IF((TN-TOUT)*H.GE.0.0E0)GO TO 555
      IF(ABS(TN-TSTOP).LE.100.0E0*UROUND*(ABS(TN)+ABS(H)))GO TO 552
      T=TN
      IDID=1
      GO TO 580
552   IDID=2
      T=TSTOP
      GO TO 580
555   CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *   IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID=3
580   CONTINUE
C--------------------------------------------------------
C     All successful returns from ddassl are made from
C     this block.
C--------------------------------------------------------
590   CONTINUE
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C-----------------------------------------------------------------------
C     This block handles all unsuccessful
C     returns other than for illegal input.
C-----------------------------------------------------------------------
600   CONTINUE
      ITEMP=-IDID
      GO TO (610,620,630,690,690,640,650,660,670,675,
     *  680,685), ITEMP
C
C     The maximum number of steps was taken before
C     reaching tout.
C
610   CALL XERRWV('DASSL--  AT CURRENT T (=R1)  500 STEPS',
     1             38,610,0,0,0,0,1,TN,ZERO)
      CALL XERRWV('DASSL--  TAKEN ON THIS CALL BEFORE REACHING TOUT',
     1             48,611,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     Tolerances too small for machine precision.
C
620   CALL XERRWV('DASSL--  AT T (=R1) TOO MUCH ACCURACY REQUESTED',
     1             47,620,0,0,0,0,1,TN,ZERO)
      CALL XERRWV('DASSL--  FOR PRECISION OF MACHINE. RTOL AND ATOL',
     1             48,621,0,0,0,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  WERE INCREASED TO APPROPRIATE VALUES',
     1             45,622,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     wt(i) .le. 0.0E0 for some i (not at start of problem).
C
630   CALL XERRWV('DASSL--  AT T (=R1) SOME ELEMENT OF WT',
     1             38,630,0,0,0,0,1,TN,ZERO)
      CALL XERRWV('DASSL--  HAS BECOME .LE. 0.0',
     1             28,631,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     Error test failed repeatedly or with h=hmin.
C
640   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     1             44,640,0,0,0,0,2,TN,H)
      CALL XERRWV(
     1'DASSL--  ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN',
     2 57,641,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     Corrector convergence failed repeatedly or with h=hmin.
C
650   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     1             44,650,0,0,0,0,2,TN,H)
      CALL XERRWV('DASSL--  CORRECTOR FAILED TO CONVERGE REPEATEDLY',
     2             48,651,0,0,0,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  OR WITH ABS(H)=HMIN',
     1             28,652,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     The iteration matrix is singular.
C
660   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     1             44,660,0,0,0,0,2,TN,H)
      CALL XERRWV('DASSL--  ITERATION MATRIX IS SINGULAR',
     1             37,661,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     Corrector failure preceeded by error test failures.
C
670   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     1             44,670,0,0,0,0,2,TN,H)
      CALL XERRWV('DASSL--  CORRECTOR COULD NOT CONVERGE.  ALSO, THE',
     1             49,671,0,0,0,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  ERROR TEST FAILED REPEATEDLY.',
     1             38,672,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     Corrector failure because ires = -1.
C
675   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     1             44,675,0,0,0,0,2,TN,H)
      CALL XERRWV('DASSL--  CORRECTOR COULD NOT CONVERGE BECAUSE',
     1             45,676,0,0,0,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  IRES WAS EQUAL TO MINUS ONE',
     1             36,677,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     Failure because ires = -2.
C
680   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2)',
     1             40,680,0,0,0,0,2,TN,H)
      CALL XERRWV('DASSL--  IRES WAS EQUAL TO MINUS TWO',
     1             36,681,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     DDAINI subroutine failed to compute initial yprime.
C
685   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     1             44,685,0,0,0,0,2,TN,HO)
      CALL XERRWV('DASSL--  INITIAL YPRIME COULD NOT BE COMPUTED',
     1             45,686,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
690   CONTINUE
      INFO(1)=-1
      T=TN
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C-----------------------------------------------------------------------
C     This block handles all error returns due
C     to illegal input, as detected before calling
C     ddastp. First the error message routine is
C     called. If this happens twice in
C     succession, execution is terminated.
C-----------------------------------------------------------------------
701   CALL XERRWV(
     1'DASSL--  SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR ONE',
     2 55,1,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
C
702   CALL XERRWV('DASSL--  NEQ (=I1) .LE. 0',
     1             25,2,0,1,NEQ,0,0,ZERO,ZERO)
      GO TO 750
C
703   CALL XERRWV('DASSL--  MAXORD (=I1) NOT IN RANGE',
     1             34,3,0,1,MXORD,0,0,ZERO,ZERO)
      GO TO 750
C
704   CALL XERRWV(
     1'DASSL--  RWORK LENGTH NEEDED, LENRW (=I1), EXCEEDS LRW (=I2)',
     2 60,4,0,2,LENRW,LRW,0,ZERO,ZERO)
      GO TO 750
C
7045  CALL XERRWV(
     1'DASSL--  SWORK LENGTH NEEDED, LENSW (=I1), EXCEEDS LSW (=I2)',
     2 60,4,0,2,LENSW,LSW,0,ZERO,ZERO)
      GO TO 750
C
705   CALL XERRWV(
     1'DASSL--  IWORK LENGTH NEEDED, LENIW (=I1), EXCEEDS LIW (=I2)',
     2 60,5,0,2,LENIW,LIW,0,ZERO,ZERO)
      GO TO 750
C
706   CALL XERRWV('DASSL--  SOME ELEMENT OF RTOL IS .LT. 0',
     1             39,6,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
C
707   CALL XERRWV('DASSL--  SOME ELEMENT OF ATOL IS .LT. 0',
     1             39,7,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
C
708   CALL XERRWV('DASSL--  ALL ELEMENTS OF RTOL AND ATOL ARE ZERO',
     1             47,8,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
C
709   CALL XERRWV(
     1'DASSL--  INFO(4) = 1 AND TSTOP (=R1) BEHIND TOUT (=R2)',
     2 54,9,0,0,0,0,2,TSTOP,TOUT)
      GO TO 750
C
710   CALL XERRWV('DASSL--  HMAX (=R1) .LT. 0.0',
     1             28,10,0,0,0,0,1,HMAX,ZERO)
      GO TO 750
C
711   CALL XERRWV('DASSL--  TOUT (=R1) BEHIND T (=R2)',
     1             34,11,0,0,0,0,2,TOUT,T)
      GO TO 750
C
712   CALL XERRWV('DASSL--  INFO(8)=1 AND H0=0.0',
     1             29,12,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
C
713   CALL XERRWV('DASSL--  SOME ELEMENT OF WT IS .LE. 0.0',
     1             39,13,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
C
714   CALL XERRWV(
     1'DASSL--  TOUT (=R1) TOO CLOSE TO T (=R2) TO START INTEGRATION',
     2 61,14,0,0,0,0,2,TOUT,T)
      GO TO 750
C
715   CALL XERRWV('DASSL--  INFO(4)=1 AND TSTOP (=R1) BEHIND T (=R2)',
     1             49,15,0,0,0,0,2,TSTOP,T)
      GO TO 750
C
717   CALL XERRWV(
     1'DASSL--  ML (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ',
     2 52,17,0,1,IWORK(LML),0,0,ZERO,ZERO)
      GO TO 750
C
718   CALL XERRWV(
     1'DASSL--  MU (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ',
     2 52,18,0,1,IWORK(LMU),0,0,ZERO,ZERO)
      GO TO 750
C
719   CALL XERRWV('DASSL--  TOUT (=R1) IS EQUAL TO T (=R2)',
     1             39,19,0,0,0,0,2,TOUT,T)
      GO TO 750
C
750   IF(INFO(1).EQ.-1) GO TO 760
      INFO(1)=-1
      IDID=-33
      RETURN
C
760   CALL XERRWV('DASSL--  REPEATED OCCURRENCES OF ILLEGAL INPUT',
     1             46,801,0,0,0,0,0,ZERO,ZERO)
C
770   CALL XERRWV('DASSL--  RUN TERMINATED. APPARENT INFINITE LOOP',
     1             47,802,1,0,0,0,0,ZERO,ZERO)
      RETURN
C
      END
      SUBROUTINE DDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,DTEM,
     1  IER,WT,E,WM,IWM,RES,IRES,UROUND,JAC,RPAR,IPAR,JFACT)
C---------------------------------------------------------------------
C     This routine computes the iteration matrix
C     pd=dg/dy+cj*dg/dyprime (where g(x,y,yprime)=0).
C     Here pd is computed by the user-supplied
C     routine jac if iwm(mtype) is 1 or 4, and
C     it is computed by numerical finite differencing
C     if iwm(mtype)is 2 or 5.
C-----------------------------------------------------------------------
C*****precision > double
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
C*****END precision > double
      EXTERNAL RES,JAC,SIGN77
      DIMENSION Y(*), YPRIME(*), DELTA(*), WT(*), E(*)
      DIMENSION WM(*), IWM(*), RPAR(*), IPAR(*), DTEM(*)
      COMMON/DDA001/NPD,NTEMP,
     *  LML,LMU,LMXORD,LMTYPE,
     *  LNST,LNSTL,LNRE,LNJE,LETF,LCTF,LIPVT
      IER = 0
      NPDM1=NPD-1
      MTYPE=IWM(LMTYPE)
      GO TO (100,200,300,400,500),MTYPE
C
C     Dense user-supplied matrix.
C
100   LENPD=NEQ*NEQ
      DO 110 I=1,LENPD
110      WM(NPDM1+I)=0.0E0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      GO TO 230
C
C     Dense finite-difference-generated matrix.
C
200   IRES=0
      NROW=NPDM1
      SQUR = SQRT(UROUND)
      DO 210 I=1,NEQ
         DEL=SQUR*MAX(ABS(Y(I)),ABS(H*YPRIME(I)),
     *     ABS(WT(I)))
         DEL=SIGN77(DEL,H*YPRIME(I))
         DEL=(Y(I)+DEL)-Y(I)
         YSAVE=Y(I)
         YPSAVE=YPRIME(I)
         Y(I)=Y(I)+DEL
         YPRIME(I)=YPRIME(I)+CJ*DEL
         IF(JFACT.EQ.0)GO TO 600
         CALL RES(X,Y,YPRIME,DTEM,IRES,RPAR,IPAR)
         IF (IRES .LT. 0) RETURN
         DELINV=1.0E0/DEL
         DO 220 L=1,NEQ
220         WM(NROW+L)=(DTEM(L)-DELTA(L))*DELINV
         GO TO 610
600      CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
         IF(IRES.LT.0)RETURN
         DELINV=1.0E0/DEL
         DO 620 L=1,NEQ
620      WM(NROW+L)=(E(L)-DELTA(L))*DELINV
610        CONTINUE
      NROW=NROW+NEQ
      Y(I)=YSAVE
      YPRIME(I)=YPSAVE
210   CONTINUE
C
C     Do dense-matrix lu decomposition on pd.
C
C*****precision > double
230      CALL DGEFA(WM(NPD),NEQ,NEQ,IWM(LIPVT),IER)
C*****END precision > double
C*****precision > single
C230      CALL SGEFA(WM(NPD),NEQ,NEQ,IWM(LIPVT),IER)
C*****END precision > single
      RETURN
C
C     Dummy section for iwm(mtype)=3.
C
300   RETURN
C
C     Banded user-supplied matrix.
C
400   LENPD=(2*IWM(LML)+IWM(LMU)+1)*NEQ
      DO 410 I=1,LENPD
410      WM(NPDM1+I)=0.0E0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 550
C
C     Banded finite-difference-generated matrix.
C
500   MBAND=IWM(LML)+IWM(LMU)+1
      MBA=MIN0(MBAND,NEQ)
      MEBAND=MBAND+IWM(LML)
      MEB1=MEBAND-1
      MSAVE=(NEQ/MBAND)+1
      ISAVE=NTEMP-1
      IPSAVE=ISAVE+MSAVE
      IRES=0
      SQUR=SQRT(UROUND)
      DO 540 J=1,MBA
         DO 510 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          WM(ISAVE+K)=Y(N)
          WM(IPSAVE+K)=YPRIME(N)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),
     *      ABS(WT(N)))
          DEL=SIGN77(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          Y(N)=Y(N)+DEL
510       YPRIME(N)=YPRIME(N)+CJ*DEL
      IF(JFACT.EQ.0)GO TO 700
      CALL RES(X,Y,YPRIME,DTEM,IRES,RPAR,IPAR)
      IF(IRES.LT.0)RETURN
      DO 730 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          Y(N)=WM(ISAVE+K)
          YPRIME(N)=WM(IPSAVE+K)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),
     *      ABS(WT(N)))
          DEL=SIGN77(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          DELINV=1.0E0/DEL
          I1=MAX0(1,(N-IWM(LMU)))
          I2=MIN0(NEQ,(N+IWM(LML)))
          II=N*MEB1-IWM(LML)+NPDM1
          DO 720 I=I1,I2
720         WM(II+I)=(DTEM(I)-DELTA(I))*DELINV
730      CONTINUE
      GO TO 540
700   CONTINUE
      CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) RETURN
      DO 530 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          Y(N)=WM(ISAVE+K)
          YPRIME(N)=WM(IPSAVE+K)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),
     *      ABS(WT(N)))
          DEL=SIGN77(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          DELINV=1.0E0/DEL
          I1=MAX0(1,(N-IWM(LMU)))
          I2=MIN0(NEQ,(N+IWM(LML)))
          II=N*MEB1-IWM(LML)+NPDM1
          DO 520 I=I1,I2
520         WM(II+I)=(E(I)-DELTA(I))*DELINV
530      CONTINUE
540   CONTINUE
C
C     Do lu decomposition of banded pd.
C
C*****precision > double
550   CALL DGBFA(WM(NPD),MEBAND,NEQ,
     *    IWM(LML),IWM(LMU),IWM(LIPVT),IER)
C*****END precision > double
C*****precision > single
C550   CALL SGBFA(WM(NPD),MEBAND,NEQ,
C     *    IWM(LML),IWM(LMU),IWM(LIPVT),IER)
C*****END precision > single
      RETURN
      END
       SUBROUTINE DDASTP(X,Y,YPRIME,NEQ,NSYS,
     *  RES,JAC,H,WT,JSTART,IDID,RPAR,IPAR,
     *  PHI,DELTA,E,WM,IWM,
     *  ALPHA,BETA,GAMMA,PSI,SIGMA,
     *  CJ,CJOLD,HOLD,S,HMIN,UROUND,
     *  IPHASE,JCALC,K,KOLD,NS,NONNEG,
     *  DTEM,EMAT,DRES,DFDYP)
C-----------------------------------------------------------------
C     DDASTP solves a system of differential and algebraic
C     equations for one step, normally from t to t+h.
C-----------------------------------------------------------------
C*****precision > double
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER (I-N)
C*****END precision > double
      LOGICAL CONVGD
      DIMENSION Y(*), YPRIME(*), WT(*), PHI(NEQ,*), DELTA(*), E(*)
      DIMENSION DTEM(*), EMAT(*), WM(*), IWM(*)
      DIMENSION PSI(*), ALPHA(*), BETA(*), GAMMA(*), SIGMA(*)
      DIMENSION RPAR(*), IPAR(*)
      EXTERNAL RES,JAC,DRES,DFDYP
      COMMON/DDA001/NPD,NTEMP,
     *   LML,LMU,LMXORD,LMTYPE,
     *   LNST,LNSTL,LNRE,LNJE,LETF,LCTF,LIPVT
      COMMON/DDA002/INDEX,NALG,IDFDP,ICALC,NPAR
      DATA MAXIT/4/
      DATA XRATE/0.25E0/
      DATA ZERO/0.0E0/, PT25/0.25E0/, PT5/0.5E0/, PT9/0.9E0/
C-----------------------------------------------------------------------
C     Block 1.
C     Initialize. On the first call,set
C     the order to 1 and initialize
C     other variables.
C-----------------------------------------------------------------------
C
C     Initializations for all calls.
C
      IDID=1
      XOLD=X
      NCF=0
      NSF=0
      NEF=0
      IF(JSTART .NE. 0) GO TO 120
C
C     If this is the first step,perform
C     other initializations.
C
      IWM(LETF) = 0
      IWM(LCTF) = 0
      K=1
      KOLD=0
      HOLD=0.0E0
      JSTART=1
      PSI(1)=H
      CJOLD = 1.0E0/H
      CJ = CJOLD
      S = 100.E0
      JCALC = -1
      DELNRM=1.0E0
      IPHASE = 0
      NS=0
120   CONTINUE
C-----------------------------------------------------------------------
C     Block 2.
C     Compute coefficients of formulas for
C     this step.
C-----------------------------------------------------------------------
200   CONTINUE
      KP1=K+1
      KP2=K+2
      KM1=K-1
      XOLD=X
      IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0
      NS=MIN0(NS+1,KOLD+2)
      NSP1=NS+1
      IF(KP1 .LT. NS)GO TO 230
      BETA(1)=1.0E0
      ALPHA(1)=1.0E0
      TEMP1=H
      GAMMA(1)=0.0E0
      SIGMA(1)=1.0E0
      DO 210 I=2,KP1
         TEMP2=PSI(I-1)
         PSI(I-1)=TEMP1
         BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
         TEMP1=TEMP2+H
         ALPHA(I)=H/TEMP1
         SIGMA(I)=FLOAT(I-1)*SIGMA(I-1)*ALPHA(I)
         GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
210      CONTINUE
      PSI(KP1)=TEMP1
230   CONTINUE
      ALPHAS = 0.0E0
      ALPHA0 = 0.0E0
      DO 240 I = 1,K
        ALPHAS = ALPHAS - 1.0E0/FLOAT(I)
        ALPHA0 = ALPHA0 - ALPHA(I)
240     CONTINUE
C
C     Compute leading coefficient cj.
C
      CJLAST = CJ
      CJ = -ALPHAS/H
C
C     Compute variable stepsize error coefficient ck.
C
      CK = ABS(ALPHA(KP1) + ALPHAS - ALPHA0)
      CK = MAX(CK,ALPHA(KP1))
C
C     Decide whether new jacobian is needed.
C
      TEMP1 = (1.0E0 - XRATE)/(1.0E0 + XRATE)
      TEMP2 = 1.0E0/TEMP1
      IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
      IF (CJ .NE. CJLAST) S = 100.E0
C
C     Change phi to phi star.
C
      IF(KP1 .LT. NSP1) GO TO 280
      DO 270 J=NSP1,KP1
         DO 260 I=1,NEQ
260         PHI(I,J)=BETA(J)*PHI(I,J)
270      CONTINUE
280   CONTINUE
C
C     Update time.
C
      X=X+H
C-----------------------------------------------------------------------
C     Block 3.
C     Predict the solution and derivative,
C     and solve the corrector equation.
C-----------------------------------------------------------------------
C
C     First,predict the solution and derivative.
C
300   CONTINUE
      DO 310 I=1,NEQ
         Y(I)=PHI(I,1)
310      YPRIME(I)=0.0E0
      DO 330 J=2,KP1
         DO 320 I=1,NEQ
            Y(I)=Y(I)+PHI(I,J)
320         YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
330   CONTINUE
      PNORM = DDANRM (NSYS,Y,WT,RPAR,IPAR)
C
C     Solve the corrector equation using a
C     modified Newton scheme.
C
      CONVGD= .TRUE.
      M=0
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
C
C     If indicated,reevaluate the
C     iteration matrix pd = dg/dy + cj*dg/dyprime
C     (where g(x,y,yprime)=0). Set
C     jcalc to 0 as an indicator that
C     this has been done.
C
      IF(JCALC .NE. -1)GO TO 340
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      JFACT=0
      CALL DDAJAC(NSYS,X,Y,YPRIME,DELTA,CJ,H,DTEM,
     * IER,WT,E,WM,IWM,RES,IRES,UROUND,JAC,RPAR,IPAR,JFACT)
      CJOLD=CJ
      S = 100.E0
      IF (IRES .LT. 0) GO TO 380
      IF(IER .NE. 0)GO TO 380
      NSF=0
C
C     Initialize the error accumulation vector e.
C
340   CONTINUE
      DO 345 I=1,NEQ
345      E(I)=0.0E0
      S = 100.E0
C
C     Corrector loop.
C
350   CONTINUE
C
C     Multiply residual by temp1 to accelerate convergence.
C
      TEMP1 = 2.0E0/(1.0E0 + CJ/CJOLD)
      DO 355 I = 1,NSYS
355     DELTA(I) = DELTA(I) * TEMP1
C
C     Compute a new iterate (back-substitution).
C     Store the correction in delta.
C
      CALL DDASLV(NSYS,DELTA,WM,IWM)
C
C     Update y,e,and yprime.
C
      DO 360 I=1,NSYS
         Y(I)=Y(I)-DELTA(I)
         E(I)=E(I)-DELTA(I)
360      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
C
C     Test for convergence of the iteration.
C
      DELNRM=DDANRM(NSYS,DELTA,WT,RPAR,IPAR)
      IF (DELNRM .LE. 100.E0*UROUND*PNORM) GO TO 375
      IF (M .GT. 0) GO TO 365
         OLDNRM = DELNRM
         GO TO 367
365   RATE = (DELNRM/OLDNRM)**(1.0E0/FLOAT(M))
      IF (RATE .GT. 0.90E0) GO TO 370
      S = RATE/(1.0E0 - RATE)
367   IF (S*DELNRM .LE. 0.33E0) GO TO 375
C
C     The corrector has not yet converged.
C     Update m and test whether the
C     maximum number of iterations have
C     been tried.
C
      M=M+1
      IF(M.GE.MAXIT)GO TO 370
C
C     Evaluate the residual
C     and go back to do another iteration.
C
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL RES(X,Y,YPRIME,DELTA,IRES,
     *  RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
      GO TO 350
C
C     The corrector failed to converge in maxit
C     iterations. If the iteration matrix
C     is not current,re-do the step with
C     a new iteration matrix.
C
370   CONTINUE
      IF(JCALC.EQ.0)GO TO 380
      JCALC=-1
      GO TO 300
375   CONTINUE
C
C      If nonnegativity of solution is required,set the
C      solution nonnegative,if the perturbation to do it
C      is small enough. If the change is too large,then
C      consider the corrector iteration to have failed.
C
      IF(NONNEG.EQ.0)GO TO 376
      DO 377 I = 1,NSYS
377      DELTA(I) = MIN(Y(I),ZERO)
      DELNRM = DDANRM(NSYS,DELTA,WT,RPAR,IPAR)
      IF(DELNRM .GT. 0.33E0) GO TO 380
      DO 378 I = 1,NSYS
378      E(I) = E(I) - DELTA(I)
C--------------------------------------------------------
C     Block for parametric sensitivities.
C     The modified Newton method has converged ; proceed
C     by calculating the parametric sensitivities
C     if requested by the user.
C-------------------------------------------------------------
376   IF(NSYS.EQ.NEQ)GO TO 390
C
C     Update the right hand side of the state equations
C     and the Jacobian matrix only if necessary.
C     Store the update of the right hand side in delta(i,1).
C     Set mcalc equal to zero to indicate that the update
C     has been done.
C
      MCALC=1
      MTYPE=IWM(LMTYPE)
      IF(IDFDP.EQ.1)GO TO 391
      MCALC=0
      IWM(LNRE)=IWM(LNRE)+1
      IRES=0
      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF(IRES.LT.0)GO TO 380
391   CONTINUE
      IF(JCALC.EQ.0)GO TO 393
      IF(MTYPE.EQ.1.OR.MTYPE.EQ.4)GO TO 392
      IF(MCALC.EQ.0)GO TO 392
      MCALC=0
      IWM(LNRE)=IWM(LNRE)+1
      IRES=0
      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF(IRES.LT.0)GO TO 380
392   CONTINUE
      IWM(LNJE)=IWM(LNJE)+1
      JFACT=1
      CALL DDAJAC(NSYS,X,Y,YPRIME,DELTA,CJ,H,DTEM,
     * IER,WT,E,WM,IWM,RES,IRES,UROUND,JAC,RPAR,IPAR,JFACT)
      IF(IRES.LT.0)GO TO 380
      IF(IER.NE.0)GO TO 380
393   CONTINUE
      CALL DDSENS(X,Y,YPRIME,NSYS,NEQ,E,DELTA,DTEM,
     *     WM,EMAT,RPAR,IPAR,CJ,IWM,UROUND,IRES,
     *     RES,DRES,DFDYP)
      IF(IRES.LT.0)GO TO 380
      GO TO 390
C
C     Singular iteration matrix or no convergence
C     with current iteration matrix.
C
380   CONVGD= .FALSE.
390   JCALC = 1
      IF(.NOT.CONVGD)GO TO 600
C-----------------------------------------------------------------------
C     Block 4.
C     Estimate the errors at orders k,k-1,k-2
C     as if constant stepsize was used. Estimate
C     the local error at order k and test
C     whether the current step is successful.
C-----------------------------------------------------------------------
      ENORM = DDANRM(NEQ,E,WT,RPAR,IPAR)
      ERK = SIGMA(K+1)*ENORM
      TERK = FLOAT(K+1)*ERK
      EST = ERK
      KNEW=K
      IF(K .EQ. 1)GO TO 430
      DO 405 I = 1,NEQ
405     DELTA(I) = PHI(I,KP1) + E(I)
      ERKM1=SIGMA(K)*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM1 = FLOAT(K)*ERKM1
      IF(K .GT. 2)GO TO 410
      IF(TERKM1 .LE. 0.5*TERK)GO TO 420
      GO TO 430
410   CONTINUE
      DO 415 I = 1,NEQ
415     DELTA(I) = PHI(I,K) + DELTA(I)
      ERKM2=SIGMA(K-1)*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM2 = FLOAT(K-1)*ERKM2
      IF(MAX(TERKM1,TERKM2).GT.TERK)GO TO 430
C
C     Lower the order.
C
420   CONTINUE
      KNEW=K-1
      EST = ERKM1
C
C     Calculate the local error for the current step
C     to see if the step was successful.
C
430   CONTINUE
      ERR = CK * ENORM
      IF(ERR .GT. 1.0E0)GO TO 600
C-----------------------------------------------------------------------
C     Block 5.
C     The step is successful. Determine
C     the best order and stepsize for
C     the next step. Update the differences
C     for the next step.
C-----------------------------------------------------------------------
      IDID=1
      IWM(LNST)=IWM(LNST)+1
      KDIFF=K-KOLD
      KOLD=K
      HOLD=H
C
C     Estimate the error at order k+1 unless
C     already decided to lower order, or
C     already using maximum order, or
C     stepsize not constant, or
C     order raised in previous step.
C
      IF(KNEW.EQ.KM1.OR.K.EQ.IWM(LMXORD))IPHASE=1
      IF(IPHASE .EQ. 0)GO TO 545
      IF(KNEW.EQ.KM1)GO TO 540
      IF(K.EQ.IWM(LMXORD)) GO TO 550
      IF(KP1.GE.NS.OR.KDIFF.EQ.1)GO TO 550
      DO 510 I=1,NEQ
510      DELTA(I)=E(I)-PHI(I,KP2)
      ERKP1 = (1.0E0/FLOAT(K+2))*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKP1 = FLOAT(K+2)*ERKP1
      IF(K.GT.1)GO TO 520
      IF(TERKP1.GE.0.5E0*TERK)GO TO 550
      GO TO 530
520   IF(TERKM1.LE.MIN(TERK,TERKP1))GO TO 540
      IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD))GO TO 550
C
C     Raise order.
C
530   K=KP1
      EST = ERKP1
      GO TO 550
C
C     Lower order.
C
540   K=KM1
      EST = ERKM1
      GO TO 550
C
C     If iphase = 0, increase order by one and multiply stepsize by
C     factor two.
C
545   K = KP1
      HNEW = H*2.0E0
      H = HNEW
      GO TO 575
C
C     Determine the appropriate stepsize for
C     the next step.
C
550   HNEW=H
      TEMP2=K+1
      R=(2.0E0*EST+0.0001E0)**(-1.0E0/TEMP2)
      IF(R .LT. 2.0E0) GO TO 555
      HNEW = 2.0E0*H
      GO TO 560
555   IF(R .GT. 1.0E0) GO TO 560
      R = MAX(PT5,MIN(PT9,R))
      HNEW = H*R
560   H=HNEW
C
C     Update differences for next step.
C
575   CONTINUE
      IF(KOLD.EQ.IWM(LMXORD))GO TO 585
      DO 580 I=1,NEQ
580      PHI(I,KP2)=E(I)
585   CONTINUE
      DO 590 I=1,NEQ
590      PHI(I,KP1)=PHI(I,KP1)+E(I)
      DO 595 J1=2,KP1
         J=KP1-J1+1
         DO 595 I=1,NEQ
595      PHI(I,J)=PHI(I,J)+PHI(I,J+1)
      RETURN
C-----------------------------------------------------------------------
C     Block 6.
C     The step is unsuccessful. Restore x,psi,phi.
C     Determine appropriate stepsize for
C     continuing the integration, or exit with
C     an error flag if there have been many
C     failures.
C-----------------------------------------------------------------------
600   IPHASE = 1
C
C     Restore x,phi,psi.
C
      X=XOLD
      IF(KP1.LT.NSP1)GO TO 630
      DO 620 J=NSP1,KP1
         TEMP1=1.0E0/BETA(J)
         DO 610 I=1,NEQ
610         PHI(I,J)=TEMP1*PHI(I,J)
620      CONTINUE
630   CONTINUE
      DO 640 I=2,KP1
640      PSI(I-1)=PSI(I)-H
C
C     Test whether failure is due to corrector iteration
C     or error test.
C
      IF(CONVGD)GO TO 660
      IWM(LCTF)=IWM(LCTF)+1
C
C     The newton iteration failed to converge with
C     a current iteration matrix.  Determine the cause
C     of the failure and take appropriate action.
C
      IF(IER.EQ.0)GO TO 650
C
C     The iteration matrix is singular. Reduce
C     the stepsize by a factor of 4. If
C     this happens three times in a row on
C     the same step, return with an error flag.
C
      NSF=NSF+1
      R = 0.25E0
      H=H*R
      IF (NSF .LT. 3 .AND. ABS(H) .GE. HMIN) GO TO 690
      IDID=-8
      GO TO 675
C
C     The newton iteration failed to converge for a reason
C     other than a singular iteration matrix.  If ires = -2, then
C     return.  Otherwise, reduce the stepsize and try again, unless
C     too many failures have occured.
C
650   CONTINUE
      IF (IRES .GT. -2) GO TO 655
      IDID = -11
      GO TO 675
655   NCF = NCF + 1
      R = 0.25E0
      H = H*R
      IF (NCF .LT. 10 .AND. ABS(H) .GE. HMIN) GO TO 690
      IDID = -7
      IF (IRES .LT. 0) IDID = -10
      IF (NEF .GE. 3) IDID = -9
      GO TO 675
C
C     The newton scheme converged,and the cause
C     of the failure was the error estimate
C     exceeding the tolerance.
C
660   NEF=NEF+1
      IWM(LETF)=IWM(LETF)+1
      IF (NEF .GT. 1) GO TO 665
C
C     On first error test failure, keep current order or lower
C     order by one.  Compute new stepsize based on differences
C     of the solution.
C
      K = KNEW
      TEMP2 = K + 1
      R = 0.90E0*(2.0E0*EST+0.0001E0)**(-1.0E0/TEMP2)
      R = MAX(PT25,MIN(PT9,R))
      H = H*R
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     On second error test failure, use the current order or
C     decrease order by one.  Reduce the stepsize by a factor of
C     one quarter.
C
665   IF (NEF .GT. 2) GO TO 670
      K = KNEW
      H = 0.25E0*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     On third and subsequent error test failures, set the order to
C     one and reduce the stepsize by a factor of one quarter.
C
670   K = 1
      H = 0.25E0*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     For all crashes, restore y to its last value,
C     interpolate to find yprime at last x, and return.
C
675   CONTINUE
      CALL DDATRP(X,X,Y,YPRIME,NEQ,K,PHI,PSI)
      RETURN
C
C     Go back and try this step again.
C
690   GO TO 200
      END
      SUBROUTINE DDASLV(NEQ,DELTA,WM,IWM)
C-----------------------------------------------------------------------
C     This routine manages the solution of the linear
C     system arising in the newton iteration.
C     Matrices and real temporary storage and
C     real information are stored in the array wm.
C     Integer matrix information is stored in
C     the array iwm.
C     For a dense matrix, the linpack routine
C     dgesl is called.
C     For a banded matrix,the linpack routine
C     dgbsl is called.
C-----------------------------------------------------------------------
C*****precision > double
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
C*****END precision > double
      DIMENSION DELTA(*), WM(*), IWM(*)
      COMMON/DDA001/NPD,NTEMP,LML,LMU,
     *  LMXORD,LMTYPE,
     *  LNST,LNSTL,LNRE,LNJE,LETF,LCTF,LIPVT
      MTYPE=IWM(LMTYPE)
      GO TO(100,100,300,400,400),MTYPE
C
C     Dense matrix.
C
C*****precision > double
100   CALL DGESL(WM(NPD),NEQ,NEQ,IWM(LIPVT),DELTA,0)
C*****END precision > double
C*****precision > single
C100   CALL SGESL(WM(NPD),NEQ,NEQ,IWM(LIPVT),DELTA,0)
C*****END precision > single
      RETURN
C
C     Dummy section for mtype=3.
C
300   CONTINUE
      RETURN
C
C     Banded matrix.
C
400   MEBAND=2*IWM(LML)+IWM(LMU)+1
C*****precision > double
      CALL DGBSL(WM(NPD),MEBAND,NEQ,IWM(LML),
     *  IWM(LMU),IWM(LIPVT),DELTA,0)
C*****END precision > double
C*****precision > single
C      CALL DGBSL(WM(NPD),MEBAND,NEQ,IWM(LML),
C     *  IWM(LMU),IWM(LIPVT),DELTA,0)
C*****END precision > single
      RETURN
      END
C
      SUBROUTINE DDAINI(X,Y,YPRIME,NEQ,
     *   RES,JAC,H,WT,IDID,RPAR,IPAR,
     *   PHI,DELTA,E,WM,IWM,
     *   HMIN,UROUND,DTEM,NONNEG)
C--------------------------------------------------------------------
C     DDAINI takes one step of size h or smaller
C     with the backward Euler method, to
C     find yprime at the initial time x. A modified
C     damped Newton iteration is used to
C     solve the corrector iteration.
C     The initial guess yprime is used in the
C     prediction, and in forming the iteration
C     matrix, but is not involved in the
C     error test. This may have trouble
C     converging if the initial guess is poor,
C     or if g(x,y,yprime) depends
C     nonlinearly on yprime.
C-------------------------------------------------------------------
C*****precision > double
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER (I-N)
C*****END precision > double
      LOGICAL CONVGD
      DIMENSION Y(*), YPRIME(*), WT(*), PHI(NEQ,*), DELTA(*), E(*)
      DIMENSION WM(*), IWM(*), RPAR(*), IPAR(*), DTEM(*)
      EXTERNAL RES,JAC
      COMMON/DDA001/NPD,NTEMP,
     *  LML,LMU,LMXORD,LMTYPE,
     *  LNST,LNSTL,LNRE,LNJE,LETF,LCTF,LIPVT
      DATA MAXIT/10/,MJAC/5/
      DATA DAMP/0.75E0/
      DATA ZERO/0.0E0/, PT5/0.5E0/, PT1/0.1E0/
C---------------------------------------------------
C     Block 1.
C     Initializations.
C---------------------------------------------------
      IDID=1
      NEF=0
      NCF=0
      NSF=0
      YNORM=DDANRM(NEQ,Y,WT,RPAR,IPAR)
C
C     Save y and yprime in phi.
C
      DO 100 I=1,NEQ
         PHI(I,1)=Y(I)
100      PHI(I,2)=YPRIME(I)
C----------------------------------------------------
C     Block 2.
C     Do one backward Euler step.
C----------------------------------------------------
C
C     Set up for start of corrector iteration.
C
200   CJ=1.0E0/H
      XNEW=X+H
C
C     Predict solution and derivative.
C
      DO 250 I=1,NEQ
250     Y(I)=Y(I)+H*YPRIME(I)
      JCALC=-1
      M=0
      CONVGD=.TRUE.
C
C     Corrector loop.
C
300   IWM(LNRE)=IWM(LNRE)+1
      IRES=0
      CALL RES(XNEW,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES.LT.0) GO TO 430
C
C     Evaluate the iteration matrix.
C
      IF (JCALC.NE.-1) GO TO 310
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      JFACT=0
      CALL DDAJAC(NEQ,XNEW,Y,YPRIME,DELTA,CJ,H,DTEM,
     *   IER,WT,E,WM,IWM,RES,IRES,
     *   UROUND,JAC,RPAR,IPAR,JFACT)
      S=1000000.E0
      IF (IRES.LT.0) GO TO 430
      IF (IER.NE.0) GO TO 430
      NSF=0
C
C     Multiply residual by damping factor.
C
310   CONTINUE
      DO 320 I=1,NEQ
320      DELTA(I)=DELTA(I)*DAMP
C
C     Compute a new iterate (back substitution)
C     and store the correction in delta.
C
      CALL DDASLV(NEQ,DELTA,WM,IWM)
C
C     Update y and yprime.
C
      DO 330 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
330      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
C
C     Test for convergence of the iteration.
C
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.LE.100.E0*UROUND*YNORM)
     *   GO TO 400
      IF (M.GT.0) GO TO 340
         OLDNRM=DELNRM
         GO TO 350
340   RATE=(DELNRM/OLDNRM)**(1.0E0/FLOAT(M))
      IF (RATE.GT.0.90E0) GO TO 430
      S=RATE/(1.0E0-RATE)
350   IF (S*DELNRM .LE. 0.33E0) GO TO 400
C
C     The corrector has not yet converged. Update
C     m and and test whether the maximum
C     number of iterations have been tried.
C     Every mjac iterations, get a new
C     iteration matrix.
C
      M=M+1
      IF (M.GE.MAXIT) GO TO 430
      IF ((M/MJAC)*MJAC.EQ.M) JCALC=-1
      GO TO 300
C
C     The iteration has converged.
C     Check nonnegativity constraints.
C
400   IF (NONNEG.EQ.0) GO TO 450
      DO 410 I=1,NEQ
410      DELTA(I)=MIN(Y(I),ZERO)
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.GT.0.33E0) GO TO 430
      DO 420 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
420      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
      GO TO 450
C
C     Exits from corrector loop.
C
430   CONVGD=.FALSE.
450   IF (.NOT.CONVGD) GO TO 600
C-----------------------------------------------------
C     Block 3.
C     The corrector iteration converged.
C     Do error test.
C-----------------------------------------------------
      DO 510 I=1,NEQ
510      E(I)=Y(I)-PHI(I,1)
      ERR=DDANRM(NEQ,E,WT,RPAR,IPAR)
      IF (ERR.LE.1.0E0) RETURN
C--------------------------------------------------------
C     Block 4.
C     The backward Euler step failed. Restore y
C     and yprime to their original values.
C     Reduce stepsize and try again, if
C     possible.
C---------------------------------------------------------
600   CONTINUE
      DO 610 I=1,NEQ
         Y(I)=PHI(I,1)
610      YPRIME(I)=PHI(I,2)
      IF (CONVGD) GO TO 640
      IF (IER.EQ.0) GO TO 620
         NSF=NSF+1
         H=H*0.25E0
         IF (NSF.LT.3.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
620   IF (IRES.GT.-2) GO TO 630
         IDID=-12
         RETURN
630   NCF=NCF+1
      H=H*0.25E0
      IF (NCF.LT.10.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
640   NEF=NEF+1
      R=0.90E0/(2.0E0*ERR+0.0001E0)
      R = MAX(PT1,MIN(PT5,R))
      H=H*R
      IF (ABS(H).GE.HMIN.AND.NEF.LT.10) GO TO 690
         IDID=-12
         RETURN
690      GO TO 200
      END
C*****precision > single
C      REAL FUNCTION DDANRM(NEQ,V,WT,RPAR,IPAR)
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      DOUBLE PRECISION FUNCTION DDANRM(NEQ,V,WT,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C-----------------------------------------------------------------------
C     This function routine computes the weighted
C     Euclidean norm of the vector of length
C     neq contained in the array v,with weights
C     contained in the array wt of length neq.
C        ddanrm=sqrt((1/neq)*sum(v(i)/wt(i))**2)
C-----------------------------------------------------------------------
      DIMENSION V(NEQ), WT(NEQ), RPAR(*), IPAR(*)
      DDANRM = 0.0E0
      VMAX = 0.0E0
      DO 10 I = 1,NEQ
10      IF(ABS(V(I)/WT(I)) .GT. VMAX) VMAX = ABS(V(I)/WT(I))
      IF(VMAX .LE. 0.0E0) GO TO 30
      SUM = 0.0E0
      DO 20 I = 1,NEQ
20      SUM = SUM + ((V(I)/WT(I))/VMAX)**2
      DDANRM = VMAX*SQRT(SUM/FLOAT(NEQ))
30    CONTINUE
      RETURN
      END
      SUBROUTINE DDATRP(X,XOUT,YOUT,YPOUT,NEQ,KOLD,PHI,PSI)
C-----------------------------------------------------------------------
C     The methods in subroutine dastep use polynomials
C     to approximate the solution. DDATRP approximates the
C     solution and its derivative at time xout by evaluating
C     one of these polynomials,and its derivative,there.
C     Information defining this polynomial is passed from
C     dastep, so ddatrp cannot be used alone.
C-----------------------------------------------------------------------
C*****precision > double
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER (I-N)
C*****END precision > double
      DIMENSION YOUT(*), YPOUT(*), PHI(NEQ,*), PSI(*)
      KOLDP1=KOLD+1
      TEMP1=XOUT-X
      DO 10 I=1,NEQ
         YOUT(I)=PHI(I,1)
10       YPOUT(I)=0.0E0
      C=1.0E0
      D=0.0E0
      GAMMA=TEMP1/PSI(1)
      DO 30 J=2,KOLDP1
         D=D*GAMMA+C/PSI(J-1)
         C=C*GAMMA
         GAMMA=(TEMP1+PSI(J-1))/PSI(J)
         DO 20 I=1,NEQ
            YOUT(I)=YOUT(I)+C*PHI(I,J)
20          YPOUT(I)=YPOUT(I)+D*PHI(I,J)
30       CONTINUE
      RETURN
      END
      SUBROUTINE DDAWTS(NEQ,IWT,RTOL,ATOL,Y,WT,RPAR,IPAR)
C-----------------------------------------------------------------------
C     This subroutine sets the error weight vector
C     wt according to wt(i)=rtol(i)*abs(y(i))+atol(i),
C     i=1,-,n.
C     RTOL and ATOL are scalars if iwt = 0,
C     and vectors if iwt = 1.
C-----------------------------------------------------------------------
C*****precision > double
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
C*****END precision > double
      DIMENSION RTOL(*), ATOL(*), Y(*), WT(*), RPAR(*), IPAR(*)
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 20 I=1,NEQ
         IF (IWT .EQ.0) GO TO 10
           RTOLI=RTOL(I)
           ATOLI=ATOL(I)
10         WT(I)=RTOLI*ABS(Y(I))+ATOLI
20         CONTINUE
      RETURN
      END
C*****precision > single
C      SUBROUTINE SDINIT(NSYS,NPAR,T,RPAR,IPAR,YPRIME,
C     *                    Y,IRES,ICALC,IDFDP,RES,DRES)
C*****END precision > single
C*****precision > double
      SUBROUTINE DDINIT(NSYS,NPAR,T,RPAR,IPAR,YPRIME,
     *                    Y,IRES,ICALC,IDFDP,RES,DRES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
C*****END precision > double
C---------------------------------------------------------------------
C      This subroutine may be used to compute the initial
C      values of yprime of the state equations as well as
C      the initial values of the sensitivity functions and
C      the initial derivatives of the sensitivities. The
C      subroutine will give accurate results only when the
C      initial conditions of the state equations are
C      inependent of the parameter vector p and the matrix
C      df/dyprime is a unit matrix.
C---------------------------------------------------------------------
C      Modified by A. Lutz to declare RES, DRES as external routines
C      so names may be passed thru argument list.  (4/16/87)
C---------------------------------------------------------------------
        DIMENSION Y(NSYS,*), YPRIME(NSYS,*), RPAR(*), IPAR(*)
        EXTERNAL RES, DRES
C
C       First calculate the initial value of the derivatives
C       of the state equations .
C
        DO 90 I=1,NSYS
        YPRIME(I,1)=0.0E0
90        CONTINUE
        CALL RES(T,Y(1,1),YPRIME(1,1),YPRIME(1,1),IRES,RPAR,IPAR)
        DO 91 I=1,NSYS
91        YPRIME(I,1)=-YPRIME(I,1)
C
C         If sensitivity calculations have been requested,
C            compute the initial values of the sensitivities
C       and their time derivatives.
C
        IF(ICALC.EQ.0)GO TO 31
        DO 92 I=1,NSYS
        DO 92 J=2,NPAR+1
92        Y(I,J)=0.0E0
        IF(IDFDP.EQ.1)GO TO 100
C*****precision > single
C        UROUND = R1MACH(4)
C*****END precision > single
C*****precision > double
        UROUND = D1MACH(4)
C*****END precision > double
        SQUR=SQRT(UROUND)
        DO 210 I=1,NPAR
        DEL=SQUR*RPAR(I)
        IF(RPAR(I).EQ.0.0E0)DEL=SQUR
        RSAVE=RPAR(I)
        RPAR(I)=RPAR(I)+DEL
        CALL RES(T,Y(1,1),YPRIME(1,1),YPRIME(1,I+1),IRES,RPAR,IPAR)
        DELINV=1.0E0/DEL
        DO 220 L=1,NSYS
220        YPRIME(L,I+1)=YPRIME(L,I+1)*DELINV
        RPAR(I)=RSAVE
210        CONTINUE
        GO TO 110
100        DO 93 I=1,NPAR
93        CALL DRES(T,Y(1,1),YPRIME(1,1),YPRIME(1,1+I),I,RPAR,IPAR)
110        DO 94 I=1,NSYS
           DO 94 J=2,NPAR+1
94        YPRIME(I,J)=-YPRIME(I,J)
31        CONTINUE
        RETURN
        END
C
        SUBROUTINE DDSENS(T,Y,YPRIME,NSYS,NEQ,E,DELTA,
     *        DTEM,WM,EMAT,RPAR,IPAR,CJ,IWM,UROUND,IRES,
     *  RES,DRES,DFDYP)
C---------------------------------------------------------------------
C       This routine controls the strategy for the solution
C        of the sensitivity equations.
C--------------------------------------------------------------------
C*****precision > double
        IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
C*****END precision > double
        DIMENSION Y(NSYS,*), YPRIME(NSYS,*), E(NSYS,*), DELTA(NSYS,*)
        DIMENSION EMAT(*), WM(*), RPAR(*), IWM(*), DTEM(*), IPAR(*)
        EXTERNAL RES,DRES,DFDYP
        COMMON/DDA001/NPD,NTEMP,
     *    LML,LMU,LMXORD,LMTYPE,
     *    LNST,LNSTL,LNRE,LNJE,LEFT,LCTF,LIPVT
        COMMON/DDA002/INDEX,NALG,IDFDP,ICALC,NPAR
C
C         Before the update initialize the arrays where the sensitivity
C       corrections are to be stored.
C
        DO 10 I=1,NSYS
        DO 10 J=2,NPAR+1
        DELTA(I,J)=0.0E0
10        CONTINUE
C
C        Compute the partial derivatives of f(t,y,yprime,p)
C        withe respect to parameter vector p. Store the
C        computation in delta(i,j) temporarily.
C
        IF(IDFDP.EQ.1)GO TO 30
        DO 20 IPARM=1,NPAR
        CALL DFDPJ(T,Y(1,1),YPRIME(1,1),NSYS,DTEM,IPARM,
     *  DELTA(1,1),DELTA(1,IPARM+1),RPAR,IPAR,UROUND,IRES,RES)
        IF(IRES.LT.0)RETURN
20        CONTINUE
        GO TO 50
30        DO 40 IPARM=1,NPAR
        CALL DRES(T,Y(1,1),YPRIME(1,1),DELTA(1,IPARM+1),
     *  IPARM,RPAR,IPAR)
40        CONTINUE
50        CONTINUE
C
C         Now compute the right hand side of the sensitivity
C        equations. Store the result in delta(i,j).
C
        IF(INDEX.EQ.2)
     *        CALL DFDYP(T,Y(1,1),YPRIME(1,1),EMAT,RPAR,IPAR)
        DO 60 IPARM=1,NPAR
        CALL DDSRHS(T,Y(1,IPARM+1),YPRIME(1,IPARM+1),NSYS,EMAT,
     *  DELTA(1,IPARM+1),CJ,RPAR,IPAR)
60        CONTINUE
C
C        Now use the LU decomposition of the iteration
C        matrix to compute the solution of the sensitivity
C        equations. Store the solution elements in delta(i,j)
C
        DO 70 IPARM=1,NPAR
        CALL DDASLV(NSYS,DELTA(1,IPARM+1),WM,IWM)
70        CONTINUE
C
C        Update the sensitivity coefficients and the error
C        vector and return.
C
        DO 80 I=1,NSYS
        DO 80 J=2,NPAR+1
        TEMP=DELTA(I,J)
        DELTA(I,J)=Y(I,J)+TEMP
        Y(I,J)=-TEMP
        E(I,J)=E(I,J)-DELTA(I,J)
        YPRIME(I,J)=YPRIME(I,J)-CJ*DELTA(I,J)
80        CONTINUE
        RETURN
        END
        SUBROUTINE DFDPJ(T,Y,YPRIME,NSYS,DTEM,IPARM,DNOM,
     *  DELTA,RPAR,IPAR,UROUND,IRES,RES)
C---------------------------------------------------------------------
C        This routine computes the partial derivatives
C        of the function f(t,y,yprime,p) with respect to
C        the parameter vector p.
C----------------------------------------------------------------------
C*****precision > double
        IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
C*****END precision > double
        DIMENSION Y(*), YPRIME(*), DTEM(*), DNOM(*), DELTA(*)
        DIMENSION RPAR(*), IPAR(*)
        EXTERNAL RES
C
C        Compute the perturbation step and perturb the parameter
C        rpar(iparm) to compute the derivatives.
C
        SQUR=SQRT(UROUND)
        DEL=SQUR*RPAR(IPARM)
        IF(RPAR(IPARM).EQ.0.0E0)DEL=SQUR
        SAVE=RPAR(IPARM)
        RPAR(IPARM)=RPAR(IPARM)+DEL
        IRES=0
        CALL RES(T,Y,YPRIME,DTEM,IRES,RPAR,IPAR)
        IF(IRES.LT.0)RETURN
        DELINV=1.0E0/DEL
        DO 10 I=1,NSYS
           DELTA(I)=(DTEM(I)-DNOM(I))*DELINV
10        CONTINUE
        RPAR(IPARM)=SAVE
        RETURN
        END
        SUBROUTINE DDSRHS(T,Y,YPRIME,NSYS,EMAT,DELTA,
     *        CJ,RPAR,IPAR)
C-----------------------------------------------------------------------
C        This routine computes the right hand side of the
C        sensitivity equations and returns the results
C        through the array delta(i,j).
C-----------------------------------------------------------------------
C*****precision > double
        IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
C*****END precision > double
        DIMENSION Y(*), YPRIME(*), DELTA(*), EMAT(NSYS,*)
        COMMON/DDA002/INDEX,NALG,IDFDP,ICALC,NPAR
        INP1=INDEX+1
        GO TO (10,10,20)INP1
10        NDIF=NSYS-NALG
        DO 40 I=1,NSYS
        FAC=0.0E0
        IF(I.LE.NDIF)FAC=1.0E0
        DELTA(I)=DELTA(I)+FAC*(YPRIME(I)-CJ*Y(I))
40        CONTINUE
        RETURN
20        DO 60 I=1,NSYS
        SUM=0.0E0
        DO 70 J=1,NSYS
        SUM=SUM+EMAT(I,J)*(YPRIME(J)-CJ*Y(J))
70        CONTINUE
        DELTA(I)=DELTA(I)+SUM
60        CONTINUE
        RETURN
        END
