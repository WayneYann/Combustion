C     CVS $Revision: 1.1.1.1 $ repositied $Date: 2006/05/26 19:09:33 $

C///////////////////////////////////////////////////////////////////////
C
C     CHEMKIN-III file dassl.f
C
C     the public-domain (?) subroutines known as DASSL
C
C
C     97/01/xx
C
C        1. Original chemkin version 3.0.
C
C     97/10/29  Joseph Grcar
C
C        1. Replaced SIGN by SIGN77.
C
C///////////////////////////////////////////////////////////////////////

C*****precision > double
      SUBROUTINE DDAINI (X,Y, YPRIME, NEQ, RES, JAC, H, WT, IDID, RPAR,
     +   IPAR, PHI, DELTA, E, WM, IWM, HMIN, UROUND, NONNEG, NTEMP)
C***BEGIN PROLOGUE  DDAINI
C***SUBSIDIARY
C***PURPOSE  Initialization routine for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDAINI-S, DDAINI-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------
C     DDAINI TAKES ONE STEP OF SIZE H OR SMALLER
C     WITH THE BACKWARD EULER METHOD, TO
C     FIND YPRIME.  X AND Y ARE UPDATED TO BE CONSISTENT WITH THE
C     NEW STEP.  A MODIFIED DAMPED NEWTON ITERATION IS USED TO
C     SOLVE THE CORRECTOR ITERATION.
C
C     THE INITIAL GUESS FOR YPRIME IS USED IN THE
C     PREDICTION, AND IN FORMING THE ITERATION
C     MATRIX, BUT IS NOT INVOLVED IN THE
C     ERROR TEST. THIS MAY HAVE TROUBLE
C     CONVERGING IF THE INITIAL GUESS IS NO
C     GOOD, OR IF G(X,Y,YPRIME) DEPENDS
C     NONLINEARLY ON YPRIME.
C
C     THE PARAMETERS REPRESENT:
C     X --         INDEPENDENT VARIABLE
C     Y --         SOLUTION VECTOR AT X
C     YPRIME --    DERIVATIVE OF SOLUTION VECTOR
C     NEQ --       NUMBER OF EQUATIONS
C     H --         STEPSIZE. IMDER MAY USE A STEPSIZE
C                  SMALLER THAN H.
C     WT --        VECTOR OF WEIGHTS FOR ERROR
C                  CRITERION
C     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS
C                  IDID= 1 -- YPRIME WAS FOUND SUCCESSFULLY
C                  IDID=-12 -- DDAINI FAILED TO FIND YPRIME
C     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS
C                  THAT ARE NOT ALTERED BY DDAINI
C     PHI --       WORK SPACE FOR DDAINI
C     DELTA,E --   WORK SPACE FOR DDAINI
C     WM,IWM --    REAL AND INTEGER ARRAYS STORING
C                  MATRIX INFORMATION
C
C-----------------------------------------------------------------
C***ROUTINES CALLED  DDAJAC, DDANRM, DDASLV
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C   901030  Minor corrections to declarations.  (FNF)
C***END PROLOGUE  DDAINI
C
      INTEGER  NEQ, IDID, IPAR(*), IWM(*), NONNEG, NTEMP
      DOUBLE PRECISION
     *   X, Y(*), YPRIME(*), H, WT(*), RPAR(*), PHI(NEQ,*), DELTA(*),
     *   E(*), WM(*), HMIN, UROUND
      EXTERNAL  RES, JAC
C
      EXTERNAL  DDAJAC, DDANRM, DDASLV
      DOUBLE PRECISION  DDANRM
C
      INTEGER  I, IER, IRES, JCALC, LNJE, LNRE, M, MAXIT, MJAC, NCF,
     *   NEF, NSF
      DOUBLE PRECISION
     *   CJ, DAMP, DELNRM, ERR, OLDNRM, R, RATE, S, XOLD, YNORM
      LOGICAL  CONVGD
C
      PARAMETER (LNRE=12)
      PARAMETER (LNJE=13)
C
      DATA MAXIT/10/,MJAC/5/
      DATA DAMP/0.75D0/
C
C
C---------------------------------------------------
C     BLOCK 1.
C     INITIALIZATIONS.
C---------------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT  DDAINI
      IDID=1
      NEF=0
      NCF=0
      NSF=0
      XOLD=X
      YNORM=DDANRM(NEQ,Y,WT,RPAR,IPAR)
C
C     SAVE Y AND YPRIME IN PHI
      DO 100 I=1,NEQ
         PHI(I,1)=Y(I)
100      PHI(I,2)=YPRIME(I)
C
C
C----------------------------------------------------
C     BLOCK 2.
C     DO ONE BACKWARD EULER STEP.
C----------------------------------------------------
C
C     SET UP FOR START OF CORRECTOR ITERATION
200   CJ=1.0D0/H
      X=X+H
C
C     PREDICT SOLUTION AND DERIVATIVE
      DO 250 I=1,NEQ
250     Y(I)=Y(I)+H*YPRIME(I)
C
      JCALC=-1
      M=0
      CONVGD=.TRUE.
C
C
C     CORRECTOR LOOP.
300   IWM(LNRE)=IWM(LNRE)+1
      IRES=0
C
      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES.LT.0) GO TO 430
C
C
C     EVALUATE THE ITERATION MATRIX
      IF (JCALC.NE.-1) GO TO 310
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      CALL DDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,
     *   IER,WT,E,WM,IWM,RES,IRES,
     *   UROUND,JAC,RPAR,IPAR,NTEMP)
C
      S=1000000.D0
      IF (IRES.LT.0) GO TO 430
      IF (IER.NE.0) GO TO 430
      NSF=0
C
C
C
C     MULTIPLY RESIDUAL BY DAMPING FACTOR
310   CONTINUE
      DO 320 I=1,NEQ
320      DELTA(I)=DELTA(I)*DAMP
C
C     COMPUTE A NEW ITERATE (BACK SUBSTITUTION)
C     STORE THE CORRECTION IN DELTA
C
      CALL DDASLV(NEQ,DELTA,WM,IWM)
C
C     UPDATE Y AND YPRIME
      DO 330 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
330      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
C
C     TEST FOR CONVERGENCE OF THE ITERATION.
C
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.LE.100.D0*UROUND*YNORM)
     *   GO TO 400
C
      IF (M.GT.0) GO TO 340
         OLDNRM=DELNRM
         GO TO 350
C
340   RATE=(DELNRM/OLDNRM)**(1.0D0/M)
      IF (RATE.GT.0.90D0) GO TO 430
      S=RATE/(1.0D0-RATE)
C
350   IF (S*DELNRM .LE. 0.33D0) GO TO 400
C
C
C     THE CORRECTOR HAS NOT YET CONVERGED. UPDATE
C     M AND AND TEST WHETHER THE MAXIMUM
C     NUMBER OF ITERATIONS HAVE BEEN TRIED.
C     EVERY MJAC ITERATIONS, GET A NEW
C     ITERATION MATRIX.
C
      M=M+1
      IF (M.GE.MAXIT) GO TO 430
C
      IF ((M/MJAC)*MJAC.EQ.M) JCALC=-1
      GO TO 300
C
C
C     THE ITERATION HAS CONVERGED.
C     CHECK NONNEGATIVITY CONSTRAINTS
400   IF (NONNEG.EQ.0) GO TO 450
      DO 410 I=1,NEQ
410      DELTA(I)=MIN(Y(I),0.0D0)
C
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.GT.0.33D0) GO TO 430
C
      DO 420 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
420      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
      GO TO 450
C
C
C     EXITS FROM CORRECTOR LOOP.
430   CONVGD=.FALSE.
450   IF (.NOT.CONVGD) GO TO 600
C
C
C
C-----------------------------------------------------
C     BLOCK 3.
C     THE CORRECTOR ITERATION CONVERGED.
C     DO ERROR TEST.
C-----------------------------------------------------
C
      DO 510 I=1,NEQ
510      E(I)=Y(I)-PHI(I,1)
      ERR=DDANRM(NEQ,E,WT,RPAR,IPAR)
C
      IF (ERR.LE.1.0D0) RETURN
C
C
C
C--------------------------------------------------------
C     BLOCK 4.
C     THE BACKWARD EULER STEP FAILED. RESTORE X, Y
C     AND YPRIME TO THEIR ORIGINAL VALUES.
C     REDUCE STEPSIZE AND TRY AGAIN, IF
C     POSSIBLE.
C---------------------------------------------------------
C
600   CONTINUE
      X = XOLD
      DO 610 I=1,NEQ
         Y(I)=PHI(I,1)
610      YPRIME(I)=PHI(I,2)
C
      IF (CONVGD) GO TO 640
      IF (IER.EQ.0) GO TO 620
         NSF=NSF+1
         H=H*0.25D0
         IF (NSF.LT.3.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
620   IF (IRES.GT.-2) GO TO 630
         IDID=-12
         RETURN
630   NCF=NCF+1
      H=H*0.25D0
      IF (NCF.LT.10.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
C
640   NEF=NEF+1
      R=0.90D0/(2.0D0*ERR+0.0001D0)
      R=MAX(0.1D0,MIN(0.5D0,R))
      H=H*R
      IF (ABS(H).GE.HMIN.AND.NEF.LT.10) GO TO 690
         IDID=-12
         RETURN
690      GO TO 200
C
C-------------END OF SUBROUTINE DDAINI----------------------
      END
      SUBROUTINE DDAJAC (NEQ, X, Y, YPRIME, DELTA, CJ, H,
     +   IER, WT, E, WM, IWM, RES, IRES, UROUND, JAC, RPAR,
     +   IPAR, NTEMP)
C***BEGIN PROLOGUE  DDAJAC
C***SUBSIDIARY
C***PURPOSE  Compute the iteration matrix for DDASSL and form the
C            LU-decomposition.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDAJAC-S, DDAJAC-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***DESCRIPTION
C----------------------------------------------------------------------
C     THIS ROUTINE COMPUTES THE ITERATION MATRIX
C     PD=DG/DY+CJ*DG/DYPRIME (WHERE G(X,Y,YPRIME)=0).
C     HERE PD IS COMPUTED BY THE USER-SUPPLIED
C     ROUTINE JAC IF IWM(MTYPE) IS 1 OR 4, AND
C     IT IS COMPUTED BY NUMERICAL FINITE DIFFERENCING
C     IF IWM(MTYPE)IS 2 OR 5
C     THE PARAMETERS HAVE THE FOLLOWING MEANINGS.
C     Y        = ARRAY CONTAINING PREDICTED VALUES
C     YPRIME   = ARRAY CONTAINING PREDICTED DERIVATIVES
C     DELTA    = RESIDUAL EVALUATED AT (X,Y,YPRIME)
C                (USED ONLY IF IWM(MTYPE)=2 OR 5)
C     CJ       = SCALAR PARAMETER DEFINING ITERATION MATRIX
C     H        = CURRENT STEPSIZE IN INTEGRATION
C     IER      = VARIABLE WHICH IS .NE. 0
C                IF ITERATION MATRIX IS SINGULAR,
C                AND 0 OTHERWISE.
C     WT       = VECTOR OF WEIGHTS FOR COMPUTING NORMS
C     E        = WORK SPACE (TEMPORARY) OF LENGTH NEQ
C     WM       = REAL WORK SPACE FOR MATRICES. ON
C                OUTPUT IT CONTAINS THE LU DECOMPOSITION
C                OF THE ITERATION MATRIX.
C     IWM      = INTEGER WORK SPACE CONTAINING
C                MATRIX INFORMATION
C     RES      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
C                TO EVALUATE THE RESIDUAL FUNCTION G(X,Y,YPRIME)
C     IRES     = FLAG WHICH IS EQUAL TO ZERO IF NO ILLEGAL VALUES
C                IN RES, AND LESS THAN ZERO OTHERWISE.  (IF IRES
C                IS LESS THAN ZERO, THE MATRIX WAS NOT COMPLETED)
C                IN THIS CASE (IF IRES .LT. 0), THEN IER = 0.
C     UROUND   = THE UNIT ROUNDOFF ERROR OF THE MACHINE BEING USED.
C     JAC      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
C                TO EVALUATE THE ITERATION MATRIX (THIS ROUTINE
C                IS ONLY USED IF IWM(MTYPE) IS 1 OR 4)
C----------------------------------------------------------------------
C***ROUTINES CALLED  DGBFA, DGEFA
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901010  Modified three MAX calls to be all on one line.  (FNF)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C   901101  Corrected PURPOSE.  (FNF)
C***END PROLOGUE  DDAJAC
C
      INTEGER  NEQ, IER, IWM(*), IRES, IPAR(*), NTEMP
      DOUBLE PRECISION
     *   X, Y(*), YPRIME(*), DELTA(*), CJ, H, WT(*), E(*), WM(*),
     *   UROUND, RPAR(*), SIGN77
      EXTERNAL  RES, JAC, SIGN77
C
      EXTERNAL  DGBFA, DGEFA
C
      INTEGER  I, I1, I2, II, IPSAVE, ISAVE, J, K, L, LENPD, LIPVT,
     *   LML, LMTYPE, LMU, MBA, MBAND, MEB1, MEBAND, MSAVE, MTYPE, N,
     *   NPD, NPDM1, NROW
      DOUBLE PRECISION  DEL, DELINV, SQUR, YPSAVE, YSAVE
C
      PARAMETER (NPD=1)
      PARAMETER (LML=1)
      PARAMETER (LMU=2)
      PARAMETER (LMTYPE=4)
      PARAMETER (LIPVT=21)
C
C***FIRST EXECUTABLE STATEMENT  DDAJAC
      IER = 0
      NPDM1=NPD-1
      MTYPE=IWM(LMTYPE)
      GO TO (100,200,300,400,500),MTYPE
C
C
C     DENSE USER-SUPPLIED MATRIX
100   LENPD=NEQ*NEQ
      DO 110 I=1,LENPD
110      WM(NPDM1+I)=0.0D0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      GO TO 230
C
C
C     DENSE FINITE-DIFFERENCE-GENERATED MATRIX
200   IRES=0
      NROW=NPDM1
      SQUR = SQRT(UROUND)
      DO 210 I=1,NEQ
         DEL=SQUR*MAX(ABS(Y(I)),ABS(H*YPRIME(I)),ABS(WT(I)))
         DEL=SIGN77(DEL,H*YPRIME(I))
         DEL=(Y(I)+DEL)-Y(I)
         YSAVE=Y(I)
         YPSAVE=YPRIME(I)
         Y(I)=Y(I)+DEL
         YPRIME(I)=YPRIME(I)+CJ*DEL
         CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
         IF (IRES .LT. 0) RETURN
         DELINV=1.0D0/DEL
         DO 220 L=1,NEQ
220      WM(NROW+L)=(E(L)-DELTA(L))*DELINV
      NROW=NROW+NEQ
      Y(I)=YSAVE
      YPRIME(I)=YPSAVE
210   CONTINUE
C
C
C     DO DENSE-MATRIX LU DECOMPOSITION ON PD
230      CALL DGEFA(WM(NPD),NEQ,NEQ,IWM(LIPVT),IER)
      RETURN
C
C
C     DUMMY SECTION FOR IWM(MTYPE)=3
300   RETURN
C
C
C     BANDED USER-SUPPLIED MATRIX
400   LENPD=(2*IWM(LML)+IWM(LMU)+1)*NEQ
      DO 410 I=1,LENPD
410      WM(NPDM1+I)=0.0D0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 550
C
C
C     BANDED FINITE-DIFFERENCE-GENERATED MATRIX
500   MBAND=IWM(LML)+IWM(LMU)+1
      MBA=MIN(MBAND,NEQ)
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
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
          DEL=SIGN77(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          Y(N)=Y(N)+DEL
510       YPRIME(N)=YPRIME(N)+CJ*DEL
      CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) RETURN
      DO 530 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          Y(N)=WM(ISAVE+K)
          YPRIME(N)=WM(IPSAVE+K)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
          DEL=SIGN77(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          DELINV=1.0D0/DEL
          I1=MAX(1,(N-IWM(LMU)))
          I2=MIN(NEQ,(N+IWM(LML)))
          II=N*MEB1-IWM(LML)+NPDM1
          DO 520 I=I1,I2
520         WM(II+I)=(E(I)-DELTA(I))*DELINV
530      CONTINUE
540   CONTINUE
C
C
C     DO LU DECOMPOSITION OF BANDED PD
550   CALL DGBFA(WM(NPD),MEBAND,NEQ,
     *    IWM(LML),IWM(LMU),IWM(LIPVT),IER)
      RETURN
C------END OF SUBROUTINE DDAJAC------
      END
      DOUBLE PRECISION FUNCTION DDANRM (NEQ, V, WT, RPAR, IPAR)
C***BEGIN PROLOGUE  DDANRM
C***SUBSIDIARY
C***PURPOSE  Compute vector norm for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDANRM-S, DDANRM-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***DESCRIPTION
C----------------------------------------------------------------------
C     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED
C     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH
C     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS
C     CONTAINED IN THE ARRAY WT OF LENGTH NEQ.
C        DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
C----------------------------------------------------------------------
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  DDANRM
C
      INTEGER  NEQ, IPAR(*)
      DOUBLE PRECISION  V(NEQ), WT(NEQ), RPAR(*)
C
      INTEGER  I
      DOUBLE PRECISION  SUM, VMAX
C
C***FIRST EXECUTABLE STATEMENT  DDANRM
      DDANRM = 0.0D0
      VMAX = 0.0D0
      DO 10 I = 1,NEQ
        IF(ABS(V(I)/WT(I)) .GT. VMAX) VMAX = ABS(V(I)/WT(I))
10      CONTINUE
      IF(VMAX .LE. 0.0D0) GO TO 30
      SUM = 0.0D0
      DO 20 I = 1,NEQ
20      SUM = SUM + ((V(I)/WT(I))/VMAX)**2
      DDANRM = VMAX*SQRT(SUM/NEQ)
30    CONTINUE
      RETURN
C------END OF FUNCTION DDANRM------
      END
      SUBROUTINE DDASLV (NEQ, DELTA, WM, IWM)
C***BEGIN PROLOGUE  DDASLV
C***SUBSIDIARY
C***PURPOSE  Linear system solver for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDASLV-S, DDASLV-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***DESCRIPTION
C----------------------------------------------------------------------
C     THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR
C     SYSTEM ARISING IN THE NEWTON ITERATION.
C     MATRICES AND REAL TEMPORARY STORAGE AND
C     REAL INFORMATION ARE STORED IN THE ARRAY WM.
C     INTEGER MATRIX INFORMATION IS STORED IN
C     THE ARRAY IWM.
C     FOR A DENSE MATRIX, THE LINPACK ROUTINE
C     DGESL IS CALLED.
C     FOR A BANDED MATRIX,THE LINPACK ROUTINE
C     DGBSL IS CALLED.
C----------------------------------------------------------------------
C***ROUTINES CALLED  DGBSL, DGESL
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  DDASLV
C
      INTEGER  NEQ, IWM(*)
      DOUBLE PRECISION  DELTA(*), WM(*)
C
      EXTERNAL  DGBSL, DGESL
C
      INTEGER  LIPVT, LML, LMU, LMTYPE, MEBAND, MTYPE, NPD
      PARAMETER (NPD=1)
      PARAMETER (LML=1)
      PARAMETER (LMU=2)
      PARAMETER (LMTYPE=4)
      PARAMETER (LIPVT=21)
C
C***FIRST EXECUTABLE STATEMENT  DDASLV
      MTYPE=IWM(LMTYPE)
      GO TO(100,100,300,400,400),MTYPE
C
C     DENSE MATRIX
100   CALL DGESL(WM(NPD),NEQ,NEQ,IWM(LIPVT),DELTA,0)
      RETURN
C
C     DUMMY SECTION FOR MTYPE=3
300   CONTINUE
      RETURN
C
C     BANDED MATRIX
400   MEBAND=2*IWM(LML)+IWM(LMU)+1
      CALL DGBSL(WM(NPD),MEBAND,NEQ,IWM(LML),
     *  IWM(LMU),IWM(LIPVT),DELTA,0)
      RETURN
C------END OF SUBROUTINE DDASLV------
      END
      SUBROUTINE DDASSL (RES,NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
     +   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
C***BEGIN PROLOGUE  DDASSL
C***PURPOSE  This code solves a system of differential/algebraic
C            equations of the form G(T,Y,YPRIME) = 0.
C***LIBRARY   SLATEC (DASSL)
C***CATEGORY  I1A2
C***TYPE      DOUBLE PRECISION (SDASSL-S, DDASSL-D)
C***KEYWORDS  DIFFERENTIAL/ALGEBRAIC, BACKWARD DIFFERENTIATION
C             FORMULAS, IMPLICIT DIFFERENTIAL SYSTEMS
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C             COMPUTING AND MATHEMATICS RESEARCH DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             L - 316, P.O. BOX 808,
C             LIVERMORE, CA.    94550
C***DESCRIPTION
C
C *Usage:
C
C      EXTERNAL RES, JAC
C      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR
C      DOUBLE PRECISION T, Y(NEQ), YPRIME(NEQ), TOUT, RTOL, ATOL,
C     *   RWORK(LRW), RPAR
C
C      CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
C     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
C
C
C *Arguments:
C  (In the following, all real arrays should be type DOUBLE PRECISION.)
C
C  RES:EXT     This is a subroutine which you provide to define the
C              differential/algebraic system.
C
C  NEQ:IN      This is the number of equations to be solved.
C
C  T:INOUT     This is the current value of the independent variable.
C
C  Y(*):INOUT  This array contains the solution components at T.
C
C  YPRIME(*):INOUT  This array contains the derivatives of the solution
C              components at T.
C
C  TOUT:IN     This is a point at which a solution is desired.
C
C  INFO(N):IN  The basic task of the code is to solve the system from T
C              to TOUT and return answer at TOUT.  INFO is an integer
C              array which is used to communicate exactly how you want
C              this task to be carried out.  (See below for details.)
C              N must be greater than or equal to 15.
C
C  RTOL,ATOL:INOUT  These quantities represent relative and absolute
C              error tolerances which you provide to indicate how
C              accurately you wish the solution to be computed.  You
C              may choose them to be both scalars or else both vectors.
C              Caution:  In Fortran 77, a scalar is not the same as an
C                        array of length 1.  Some compilers may object
C                        to using scalars for RTOL,ATOL.
C
C  IDID:OUT    This scalar quantity is an indicator reporting what the
C              code did.  You must monitor this integer variable to
C              decide  what action to take next.
C
C  RWORK:WORK  A real work array of length LRW which provides the
C              code with needed storage space.
C
C  LRW:IN      The length of RWORK.  (See below for required length.)
C
C  IWORK:WORK  An integer work array of length LIW which probides the
C              code with needed storage space.
C
C  LIW:IN      The length of IWORK.  (See below for required length.)
C
C  RPAR,IPAR:IN  These are real and integer parameter arrays which
C              you can use for communication between your calling
C              program and the RES subroutine (and the JAC subroutine)
C
C  JAC:EXT     This is the name of a subroutine which you may choose
C              to provide for defining a matrix of partial derivatives
C              described below.
C
C  Quantities which may be altered by DDASSL are:
C     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL,
C     IDID, RWORK(*) AND IWORK(*)
C
C *Description
C
C  Subroutine DDASSL uses the backward differentiation formulas of
C  orders one thru five to solve a system of the above form for Y and
C  YPRIME.  Values Y and YPRIME at the initial time must be given as
C  input.  These values must be consistent, (that is, if T,Y,YPRIME are
C  the initial values, they must satisfy G(T,Y,YPRIME) = 0.).  The
C  subroutine solves the system from T to TOUT.  It is easy to continue
C  solution to get results at additional TOUT.  This is the interval
C  mode of operation.  Intermediate results can also be obtained easily
C  by using the intermediate-output capability.
C
C  The following detailed description is divided into subsections:
C    1. Input required for the first call to DDASSL.
C    2. Output after any return from DDASSL.
C    3. What to do to continue the integration.
C    4. Error messages.
C
C
C  ----- INPUT -- WHAT TO DO ON THE FIRST CALL TO DDASSL ------------
C
C  The first call of the code is defined to be the start of each new
C  problem. Read through the descriptions of all the following items,
C  provide sufficient storage space for designated arrays, set
C  appropriate variables for the initialization of the problem, and
C  give information about how you want the problem to be solved.
C
C
C  RES -- Provide a subroutine of the form
C             SUBROUTINE RES(T,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
C         to define the system of differential/algebraic
C         equations which is to be solved. For the given values
C         of T,Y and YPRIME, the subroutine should
C         return the residual of the defferential/algebraic
C         system
C             DELTA = G(T,Y,YPRIME)
C         (DELTA(*) is a vector of length NEQ which is
C         output for RES.)
C
C         Subroutine RES must not alter T,Y or YPRIME.
C         You must declare the name RES in an external
C         statement in your program that calls DDASSL.
C         You must dimension Y,YPRIME and DELTA in RES.
C
C         IRES is an integer flag which is always equal to
C         zero on input. Subroutine RES should alter IRES
C         only if it encounters an illegal value of Y or
C         a stop condition. Set IRES = -1 if an input value
C         is illegal, and DDASSL will try to solve the problem
C         without getting IRES = -1. If IRES = -2, DDASSL
C         will return control to the calling program
C         with IDID = -11.
C
C         RPAR and IPAR are real and integer parameter arrays which
C         you can use for communication between your calling program
C         and subroutine RES. They are not altered by DDASSL. If you
C         do not need RPAR or IPAR, ignore these parameters by treat-
C         ing them as dummy arguments. If you do choose to use them,
C         dimension them in your calling program and in RES as arrays
C         of appropriate length.
C
C  NEQ -- Set it to the number of differential equations.
C         (NEQ .GE. 1)
C
C  T -- Set it to the initial point of the integration.
C         T must be defined as a variable.
C
C  Y(*) -- Set this vector to the initial values of the NEQ solution
C         components at the initial point. You must dimension Y of
C         length at least NEQ in your calling program.
C
C  YPRIME(*) -- Set this vector to the initial values of the NEQ
C         first derivatives of the solution components at the initial
C         point.  You must dimension YPRIME at least NEQ in your
C         calling program. If you do not know initial values of some
C         of the solution components, see the explanation of INFO(11).
C
C  TOUT -- Set it to the first point at which a solution
C         is desired. You can not take TOUT = T.
C         integration either forward in T (TOUT .GT. T) or
C         backward in T (TOUT .LT. T) is permitted.
C
C         The code advances the solution from T to TOUT using
C         step sizes which are automatically selected so as to
C         achieve the desired accuracy. If you wish, the code will
C         return with the solution and its derivative at
C         intermediate steps (intermediate-output mode) so that
C         you can monitor them, but you still must provide TOUT in
C         accord with the basic aim of the code.
C
C         The first step taken by the code is a critical one
C         because it must reflect how fast the solution changes near
C         the initial point. The code automatically selects an
C         initial step size which is practically always suitable for
C         the problem. By using the fact that the code will not step
C         past TOUT in the first step, you could, if necessary,
C         restrict the length of the initial step size.
C
C         For some problems it may not be permissible to integrate
C         past a point TSTOP because a discontinuity occurs there
C         or the solution or its derivative is not defined beyond
C         TSTOP. When you have declared a TSTOP point (SEE INFO(4)
C         and RWORK(1)), you have told the code not to integrate
C         past TSTOP. In this case any TOUT beyond TSTOP is invalid
C         input.
C
C  INFO(*) -- Use the INFO array to give the code more details about
C         how you want your problem solved.  This array should be
C         dimensioned of length 15, though DDASSL uses only the first
C         eleven entries.  You must respond to all of the following
C         items, which are arranged as questions.  The simplest use
C         of the code corresponds to answering all questions as yes,
C         i.e. setting all entries of INFO to 0.
C
C       INFO(1) - This parameter enables the code to initialize
C              itself. You must set it to indicate the start of every
C              new problem.
C
C          **** Is this the first call for this problem ...
C                Yes - Set INFO(1) = 0
C                 No - Not applicable here.
C                      See below for continuation calls.  ****
C
C       INFO(2) - How much accuracy you want of your solution
C              is specified by the error tolerances RTOL and ATOL.
C              The simplest use is to take them both to be scalars.
C              To obtain more flexibility, they can both be vectors.
C              The code must be told your choice.
C
C          **** Are both error tolerances RTOL, ATOL scalars ...
C                Yes - Set INFO(2) = 0
C                      and input scalars for both RTOL and ATOL
C                 No - Set INFO(2) = 1
C                      and input arrays for both RTOL and ATOL ****
C
C       INFO(3) - The code integrates from T in the direction
C              of TOUT by steps. If you wish, it will return the
C              computed solution and derivative at the next
C              intermediate step (the intermediate-output mode) or
C              TOUT, whichever comes first. This is a good way to
C              proceed if you want to see the behavior of the solution.
C              If you must have solutions at a great many specific
C              TOUT points, this code will compute them efficiently.
C
C          **** Do you want the solution only at
C                TOUT (and not at the next intermediate step) ...
C                 Yes - Set INFO(3) = 0
C                  No - Set INFO(3) = 1 ****
C
C       INFO(4) - To handle solutions at a great many specific
C              values TOUT efficiently, this code may integrate past
C              TOUT and interpolate to obtain the result at TOUT.
C              Sometimes it is not possible to integrate beyond some
C              point TSTOP because the equation changes there or it is
C              not defined past TSTOP. Then you must tell the code
C              not to go past.
C
C           **** Can the integration be carried out without any
C                restrictions on the independent variable T ...
C                 Yes - Set INFO(4)=0
C                  No - Set INFO(4)=1
C                       and define the stopping point TSTOP by
C                       setting RWORK(1)=TSTOP ****
C
C       INFO(5) - To solve differential/algebraic problems it is
C              necessary to use a matrix of partial derivatives of the
C              system of differential equations. If you do not
C              provide a subroutine to evaluate it analytically (see
C              description of the item JAC in the call list), it will
C              be approximated by numerical differencing in this code.
C              although it is less trouble for you to have the code
C              compute partial derivatives by numerical differencing,
C              the solution will be more reliable if you provide the
C              derivatives via JAC. Sometimes numerical differencing
C              is cheaper than evaluating derivatives in JAC and
C              sometimes it is not - this depends on your problem.
C
C           **** Do you want the code to evaluate the partial
C                derivatives automatically by numerical differences ...
C                   Yes - Set INFO(5)=0
C                    No - Set INFO(5)=1
C                  and provide subroutine JAC for evaluating the
C                  matrix of partial derivatives ****
C
C       INFO(6) - DDASSL will perform much better if the matrix of
C              partial derivatives, DG/DY + CJ*DG/DYPRIME,
C              (here CJ is a scalar determined by DDASSL)
C              is banded and the code is told this. In this
C              case, the storage needed will be greatly reduced,
C              numerical differencing will be performed much cheaper,
C              and a number of important algorithms will execute much
C              faster. The differential equation is said to have
C              half-bandwidths ML (lower) and MU (upper) if equation i
C              involves only unknowns Y(J) with
C                             I-ML .LE. J .LE. I+MU
C              for all I=1,2,...,NEQ. Thus, ML and MU are the widths
C              of the lower and upper parts of the band, respectively,
C              with the main diagonal being excluded. If you do not
C              indicate the equation has a banded matrix of partial
C              derivatives, the code works with a full matrix of NEQ**2
C              elements (stored in the conventional way). Computations
C              with banded matrices cost less time & storage than with
C              full matrices if 2*ML+MU .LT. NEQ. If you tell the
C              code that the matrix of partial derivatives has a banded
C              structure and you want to provide subroutine JAC to
C              compute partial derivatives, then you must be careful
C              to store the elements of the matrix in the special form
C              indicated in the description of JAC.
C
C          **** Do you want to solve the problem using a full
C               (dense) matrix (and not a special banded
C               structure) ...
C                Yes - Set INFO(6)=0
C                 No - Set INFO(6)=1
C                       and provide the lower (ML) and upper (MU)
C                       bandwidths by setting
C                       IWORK(1)=ML
C                       IWORK(2)=MU ****
C
C
C        INFO(7) -- You can specify a maximum (absolute value of)
C              stepsize, so that the code
C              will avoid passing over very
C              large regions.
C
C          ****  Do you want the code to decide
C                on its own maximum stepsize?
C                Yes - Set INFO(7)=0
C                 No - Set INFO(7)=1
C                      and define HMAX by setting
C                      RWORK(2)=HMAX ****
C
C        INFO(8) -- Differential/algebraic problems
C              may occaisionally suffer from
C              severe scaling difficulties on the
C              first step. If you know a great deal
C              about the scaling of your problem, you can
C              help to alleviate this problem by
C              specifying an initial stepsize HO.
C
C          ****  Do you want the code to define
C                its own initial stepsize?
C                Yes - Set INFO(8)=0
C                 No - Set INFO(8)=1
C                      and define HO by setting
C                      RWORK(3)=HO ****
C
C        INFO(9) -- If storage is a severe problem,
C              you can save some locations by
C              restricting the maximum order MAXORD.
C              the default value is 5. for each
C              order decrease below 5, the code
C              requires NEQ fewer locations, however
C              it is likely to be slower. In any
C              case, you must have 1 .LE. MAXORD .LE. 5
C          ****  Do you want the maximum order to
C                default to 5?
C                Yes - Set INFO(9)=0
C                 No - Set INFO(9)=1
C                      and define MAXORD by setting
C                      IWORK(3)=MAXORD ****
C
C        INFO(10) --If you know that the solutions to your equations
C               will always be nonnegative, it may help to set this
C               parameter. However, it is probably best to
C               try the code without using this option first,
C               and only to use this option if that doesn't
C               work very well.
C           ****  Do you want the code to solve the problem without
C                 invoking any special nonnegativity constraints?
C                  Yes - Set INFO(10)=0
C                   No - Set INFO(10)=1
C
C        INFO(11) --DDASSL normally requires the initial T,
C               Y, and YPRIME to be consistent. That is,
C               you must have G(T,Y,YPRIME) = 0 at the initial
C               time. If you do not know the initial
C               derivative precisely, you can let DDASSL try
C               to compute it.
C          ****   Are the initialHE INITIAL T, Y, YPRIME consistent?
C                 Yes - Set INFO(11) = 0
C                  No - Set INFO(11) = 1,
C                       and set YPRIME to an initial approximation
C                       to YPRIME.  (If you have no idea what
C                       YPRIME should be, set it to zero. Note
C                       that the initial Y should be such
C                       that there must exist a YPRIME so that
C                       G(T,Y,YPRIME) = 0.)
C
C  RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL
C         error tolerances to tell the code how accurately you
C         want the solution to be computed.  They must be defined
C         as variables because the code may change them.  You
C         have two choices --
C               Both RTOL and ATOL are scalars. (INFO(2)=0)
C               Both RTOL and ATOL are vectors. (INFO(2)=1)
C         in either case all components must be non-negative.
C
C         The tolerances are used by the code in a local error
C         test at each step which requires roughly that
C               ABS(LOCAL ERROR) .LE. RTOL*ABS(Y)+ATOL
C         for each vector component.
C         (More specifically, a root-mean-square norm is used to
C         measure the size of vectors, and the error test uses the
C         magnitude of the solution at the beginning of the step.)
C
C         The true (global) error is the difference between the
C         true solution of the initial value problem and the
C         computed approximation.  Practically all present day
C         codes, including this one, control the local error at
C         each step and do not even attempt to control the global
C         error directly.
C         Usually, but not always, the true accuracy of the
C         computed Y is comparable to the error tolerances. This
C         code will usually, but not always, deliver a more
C         accurate solution if you reduce the tolerances and
C         integrate again.  By comparing two such solutions you
C         can get a fairly reliable idea of the true error in the
C         solution at the bigger tolerances.
C
C         Setting ATOL=0. results in a pure relative error test on
C         that component.  Setting RTOL=0. results in a pure
C         absolute error test on that component.  A mixed test
C         with non-zero RTOL and ATOL corresponds roughly to a
C         relative error test when the solution component is much
C         bigger than ATOL and to an absolute error test when the
C         solution component is smaller than the threshhold ATOL.
C
C         The code will not attempt to compute a solution at an
C         accuracy unreasonable for the machine being used.  It will
C         advise you if you ask for too much accuracy and inform
C         you as to the maximum accuracy it believes possible.
C
C  RWORK(*) --  Dimension this real work array of length LRW in your
C         calling program.
C
C  LRW -- Set it to the declared length of the RWORK array.
C               You must have
C                    LRW .GE. 40+(MAXORD+4)*NEQ+NEQ**2
C               for the full (dense) JACOBIAN case (when INFO(6)=0), or
C                    LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
C               for the banded user-defined JACOBIAN case
C               (when INFO(5)=1 and INFO(6)=1), or
C                     LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
C                           +2*(NEQ/(ML+MU+1)+1)
C               for banded finite-difference-generated JACOBIAN case
C               (when INFO(5)=0 and INFO(6)=1)
C
C  IWORK(*) --  Dimension this integer work array of length LIW in
C         your calling program.
C
C  LIW -- Set it to the declared length of the IWORK array.
C               You must have LIW .GE. 20+NEQ
C
C  RPAR, IPAR -- These are parameter arrays, of real and integer
C         type, respectively.  You can use them for communication
C         between your program that calls DDASSL and the
C         RES subroutine (and the JAC subroutine).  They are not
C         altered by DDASSL.  If you do not need RPAR or IPAR,
C         ignore these parameters by treating them as dummy
C         arguments.  If you do choose to use them, dimension
C         them in your calling program and in RES (and in JAC)
C         as arrays of appropriate length.
C
C  JAC -- If you have set INFO(5)=0, you can ignore this parameter
C         by treating it as a dummy argument.  Otherwise, you must
C         provide a subroutine of the form
C             SUBROUTINE JAC(T,Y,YPRIME,PD,CJ,RPAR,IPAR)
C         to define the matrix of partial derivatives
C             PD=DG/DY+CJ*DG/DYPRIME
C         CJ is a scalar which is input to JAC.
C         For the given values of T,Y,YPRIME, the
C         subroutine must evaluate the non-zero partial
C         derivatives for each equation and each solution
C         component, and store these values in the
C         matrix PD.  The elements of PD are set to zero
C         before each call to JAC so only non-zero elements
C         need to be defined.
C
C         Subroutine JAC must not alter T,Y,(*),YPRIME(*), or CJ.
C         You must declare the name JAC in an EXTERNAL statement in
C         your program that calls DDASSL.  You must dimension Y,
C         YPRIME and PD in JAC.
C
C         The way you must store the elements into the PD matrix
C         depends on the structure of the matrix which you
C         indicated by INFO(6).
C               *** INFO(6)=0 -- Full (dense) matrix ***
C                   Give PD a first dimension of NEQ.
C                   When you evaluate the (non-zero) partial derivative
C                   of equation I with respect to variable J, you must
C                   store it in PD according to
C                   PD(I,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
C               *** INFO(6)=1 -- Banded JACOBIAN with ML lower and MU
C                   upper diagonal bands (refer to INFO(6) description
C                   of ML and MU) ***
C                   Give PD a first dimension of 2*ML+MU+1.
C                   when you evaluate the (non-zero) partial derivative
C                   of equation I with respect to variable J, you must
C                   store it in PD according to
C                   IROW = I - J + ML + MU + 1
C                   PD(IROW,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
C
C         RPAR and IPAR are real and integer parameter arrays
C         which you can use for communication between your calling
C         program and your JACOBIAN subroutine JAC. They are not
C         altered by DDASSL. If you do not need RPAR or IPAR,
C         ignore these parameters by treating them as dummy
C         arguments. If you do choose to use them, dimension
C         them in your calling program and in JAC as arrays of
C         appropriate length.
C
C
C  OPTIONALLY REPLACEABLE NORM ROUTINE:
C
C     DDASSL uses a weighted norm DDANRM to measure the size
C     of vectors such as the estimated error in each step.
C     A FUNCTION subprogram
C       DOUBLE PRECISION FUNCTION DDANRM(NEQ,V,WT,RPAR,IPAR)
C       DIMENSION V(NEQ),WT(NEQ)
C     is used to define this norm. Here, V is the vector
C     whose norm is to be computed, and WT is a vector of
C     weights.  A DDANRM routine has been included with DDASSL
C     which computes the weighted root-mean-square norm
C     given by
C       DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
C     this norm is suitable for most problems. In some
C     special cases, it may be more convenient and/or
C     efficient to define your own norm by writing a function
C     subprogram to be called instead of DDANRM. This should,
C     however, be attempted only after careful thought and
C     consideration.
C
C
C  ----- OUTPUT -- AFTER ANY RETURN FROM DDASSL ---------------------
C
C  The principal aim of the code is to return a computed solution at
C  TOUT, although it is also possible to obtain intermediate results
C  along the way. To find out whether the code achieved its goal
C  or if the integration process was interrupted before the task was
C  completed, you must check the IDID parameter.
C
C
C  T -- The solution was successfully advanced to the
C               output value of T.
C
C  Y(*) -- Contains the computed solution approximation at T.
C
C  YPRIME(*) -- Contains the computed derivative
C               approximation at T.
C
C  IDID -- Reports what the code did.
C
C                     *** Task completed ***
C                Reported by positive values of IDID
C
C           IDID = 1 -- A step was successfully taken in the
C                   intermediate-output mode. The code has not
C                   yet reached TOUT.
C
C           IDID = 2 -- The integration to TSTOP was successfully
C                   completed (T=TSTOP) by stepping exactly to TSTOP.
C
C           IDID = 3 -- The integration to TOUT was successfully
C                   completed (T=TOUT) by stepping past TOUT.
C                   Y(*) is obtained by interpolation.
C                   YPRIME(*) is obtained by interpolation.
C
C                    *** Task interrupted ***
C                Reported by negative values of IDID
C
C           IDID = -1 -- A large amount of work has been expended.
C                   (About 500 steps)
C
C           IDID = -2 -- The error tolerances are too stringent.
C
C           IDID = -3 -- The local error test cannot be satisfied
C                   because you specified a zero component in ATOL
C                   and the corresponding computed solution
C                   component is zero. Thus, a pure relative error
C                   test is impossible for this component.
C
C           IDID = -6 -- DDASSL had repeated error test
C                   failures on the last attempted step.
C
C           IDID = -7 -- The corrector could not converge.
C
C           IDID = -8 -- The matrix of partial derivatives
C                   is singular.
C
C           IDID = -9 -- The corrector could not converge.
C                   there were repeated error test failures
C                   in this step.
C
C           IDID =-10 -- The corrector could not converge
C                   because IRES was equal to minus one.
C
C           IDID =-11 -- IRES equal to -2 was encountered
C                   and control is being returned to the
C                   calling program.
C
C           IDID =-12 -- DDASSL failed to compute the initial
C                   YPRIME.
C
C
C
C           IDID = -13,..,-32 -- Not applicable for this code
C
C                    *** Task terminated ***
C                Reported by the value of IDID=-33
C
C           IDID = -33 -- The code has encountered trouble from which
C                   it cannot recover. A message is printed
C                   explaining the trouble and control is returned
C                   to the calling program. For example, this occurs
C                   when invalid input is detected.
C
C  RTOL, ATOL -- These quantities remain unchanged except when
C               IDID = -2. In this case, the error tolerances have been
C               increased by the code to values which are estimated to
C               be appropriate for continuing the integration. However,
C               the reported solution at T was obtained using the input
C               values of RTOL and ATOL.
C
C  RWORK, IWORK -- Contain information which is usually of no interest
C               to the user but necessary for subsequent calls.
C               However, you may find use for
C
C               RWORK(3)--Which contains the step size H to be
C                       attempted on the next step.
C
C               RWORK(4)--Which contains the current value of the
C                       independent variable, i.e., the farthest point
C                       integration has reached. This will be different
C                       from T only when interpolation has been
C                       performed (IDID=3).
C
C               RWORK(7)--Which contains the stepsize used
C                       on the last successful step.
C
C               IWORK(7)--Which contains the order of the method to
C                       be attempted on the next step.
C
C               IWORK(8)--Which contains the order of the method used
C                       on the last step.
C
C               IWORK(11)--Which contains the number of steps taken so
C                        far.
C
C               IWORK(12)--Which contains the number of calls to RES
C                        so far.
C
C               IWORK(13)--Which contains the number of evaluations of
C                        the matrix of partial derivatives needed so
C                        far.
C
C               IWORK(14)--Which contains the total number
C                        of error test failures so far.
C
C               IWORK(15)--Which contains the total number
C                        of convergence test failures so far.
C                        (includes singular iteration matrix
C                        failures.)
C
C
C  ----- INPUT -- WHAT TO DO TO CONTINUE THE INTEGRATION ------------
C                    (CALLS AFTER THE FIRST)
C
C  This code is organized so that subsequent calls to continue the
C  integration involve little (if any) additional effort on your
C  part. You must monitor the IDID parameter in order to determine
C  what to do next.
C
C  Recalling that the principal task of the code is to integrate
C  from T to TOUT (the interval mode), usually all you will need
C  to do is specify a new TOUT upon reaching the current TOUT.
C
C  Do not alter any quantity not specifically permitted below,
C  in particular do not alter NEQ,T,Y(*),YPRIME(*),RWORK(*),IWORK(*)
C  or the differential equation in subroutine RES. Any such
C  alteration constitutes a new problem and must be treated as such,
C  i.e., you must start afresh.
C
C  You cannot change from vector to scalar error control or vice
C  versa (INFO(2)), but you can change the size of the entries of
C  RTOL, ATOL. Increasing a tolerance makes the equation easier
C  to integrate. Decreasing a tolerance will make the equation
C  harder to integrate and should generally be avoided.
C
C  You can switch from the intermediate-output mode to the
C  interval mode (INFO(3)) or vice versa at any time.
C
C  If it has been necessary to prevent the integration from going
C  past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
C  code will not integrate to any TOUT beyond the currently
C  specified TSTOP. Once TSTOP has been reached you must change
C  the value of TSTOP or set INFO(4)=0. You may change INFO(4)
C  or TSTOP at any time but you must supply the value of TSTOP in
C  RWORK(1) whenever you set INFO(4)=1.
C
C  Do not change INFO(5), INFO(6), IWORK(1), or IWORK(2)
C  unless you are going to restart the code.
C
C                 *** Following a completed task ***
C  If
C     IDID = 1, call the code again to continue the integration
C                  another step in the direction of TOUT.
C
C     IDID = 2 or 3, define a new TOUT and call the code again.
C                  TOUT must be different from T. You cannot change
C                  the direction of integration without restarting.
C
C                 *** Following an interrupted task ***
C               To show the code that you realize the task was
C               interrupted and that you want to continue, you
C               must take appropriate action and set INFO(1) = 1
C  If
C    IDID = -1, The code has taken about 500 steps.
C                  If you want to continue, set INFO(1) = 1 and
C                  call the code again. An additional 500 steps
C                  will be allowed.
C
C    IDID = -2, The error tolerances RTOL, ATOL have been
C                  increased to values the code estimates appropriate
C                  for continuing. You may want to change them
C                  yourself. If you are sure you want to continue
C                  with relaxed error tolerances, set INFO(1)=1 and
C                  call the code again.
C
C    IDID = -3, A solution component is zero and you set the
C                  corresponding component of ATOL to zero. If you
C                  are sure you want to continue, you must first
C                  alter the error criterion to use positive values
C                  for those components of ATOL corresponding to zero
C                  solution components, then set INFO(1)=1 and call
C                  the code again.
C
C    IDID = -4,-5  --- Cannot occur with this code.
C
C    IDID = -6, Repeated error test failures occurred on the
C                  last attempted step in DDASSL. A singularity in the
C                  solution may be present. If you are absolutely
C                  certain you want to continue, you should restart
C                  the integration. (Provide initial values of Y and
C                  YPRIME which are consistent)
C
C    IDID = -7, Repeated convergence test failures occurred
C                  on the last attempted step in DDASSL. An inaccurate
C                  or ill-conditioned JACOBIAN may be the problem. If
C                  you are absolutely certain you want to continue, you
C                  should restart the integration.
C
C    IDID = -8, The matrix of partial derivatives is singular.
C                  Some of your equations may be redundant.
C                  DDASSL cannot solve the problem as stated.
C                  It is possible that the redundant equations
C                  could be removed, and then DDASSL could
C                  solve the problem. It is also possible
C                  that a solution to your problem either
C                  does not exist or is not unique.
C
C    IDID = -9, DDASSL had multiple convergence test
C                  failures, preceeded by multiple error
C                  test failures, on the last attempted step.
C                  It is possible that your problem
C                  is ill-posed, and cannot be solved
C                  using this code. Or, there may be a
C                  discontinuity or a singularity in the
C                  solution. If you are absolutely certain
C                  you want to continue, you should restart
C                  the integration.
C
C    IDID =-10, DDASSL had multiple convergence test failures
C                  because IRES was equal to minus one.
C                  If you are absolutely certain you want
C                  to continue, you should restart the
C                  integration.
C
C    IDID =-11, IRES=-2 was encountered, and control is being
C                  returned to the calling program.
C
C    IDID =-12, DDASSL failed to compute the initial YPRIME.
C                  This could happen because the initial
C                  approximation to YPRIME was not very good, or
C                  if a YPRIME consistent with the initial Y
C                  does not exist. The problem could also be caused
C                  by an inaccurate or singular iteration matrix.
C
C    IDID = -13,..,-32  --- Cannot occur with this code.
C
C
C                 *** Following a terminated task ***
C
C  If IDID= -33, you cannot continue the solution of this problem.
C                  An attempt to do so will result in your
C                  run being terminated.
C
C
C  ----- ERROR MESSAGES ---------------------------------------------
C
C      The SLATEC error print routine XERMSG is called in the event of
C   unsuccessful completion of a task.  Most of these are treated as
C   "recoverable errors", meaning that (unless the user has directed
C   otherwise) control will be returned to the calling program for
C   possible action after the message has been printed.
C
C   In the event of a negative value of IDID other than -33, an appro-
C   priate message is printed and the "error number" printed by XERMSG
C   is the value of IDID.  There are quite a number of illegal input
C   errors that can lead to a returned value IDID=-33.  The conditions
C   and their printed "error numbers" are as follows:
C
C   Error number       Condition
C
C        1       Some element of INFO vector is not zero or one.
C        2       NEQ .le. 0
C        3       MAXORD not in range.
C        4       LRW is less than the required length for RWORK.
C        5       LIW is less than the required length for IWORK.
C        6       Some element of RTOL is .lt. 0
C        7       Some element of ATOL is .lt. 0
C        8       All elements of RTOL and ATOL are zero.
C        9       INFO(4)=1 and TSTOP is behind TOUT.
C       10       HMAX .lt. 0.0
C       11       TOUT is behind T.
C       12       INFO(8)=1 and H0=0.0
C       13       Some element of WT is .le. 0.0
C       14       TOUT is too close to T to start integration.
C       15       INFO(4)=1 and TSTOP is behind T.
C       16       --( Not used in this version )--
C       17       ML illegal.  Either .lt. 0 or .gt. NEQ
C       18       MU illegal.  Either .lt. 0 or .gt. NEQ
C       19       TOUT = T.
C
C   If DDASSL is called again without any action taken to remove the
C   cause of an unsuccessful return, XERMSG will be called with a fatal
C   error flag, which will cause unconditional termination of the
C   program.  There are two such fatal errors:
C
C   Error number -998:  The last step was terminated with a negative
C       value of IDID other than -33, and no appropriate action was
C       taken.
C
C   Error number -999:  The previous call was terminated because of
C       illegal input (IDID=-33) and there is illegal input in the
C       present call, as well.  (Suspect infinite loop.)
C
C  -----------------------------------------------------------------
C
C***REFERENCES  A DESCRIPTION OF DASSL: A DIFFERENTIAL/ALGEBRAIC
C                 SYSTEM SOLVER, L. R. PETZOLD, SAND82-8637,
C                 SANDIA NATIONAL LABORATORIES, SEPTEMBER 1982.
C***ROUTINES CALLED  D1MACH, DDAINI, DDANRM, DDASTP, DDATRP, DDAWTS,
C                    XERMSG
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   880387  Code changes made.  All common statements have been
C           replaced by a DATA statement, which defines pointers into
C           RWORK, and PARAMETER statements which define pointers
C           into IWORK.  As well the documentation has gone through
C           grammatical changes.
C   881005  The prologue has been changed to mixed case.
C           The subordinate routines had revision dates changed to
C           this date, although the documentation for these routines
C           is all upper case.  No code changes.
C   890511  Code changes made.  The DATA statement in the declaration
C           section of DDASSL was replaced with a PARAMETER
C           statement.  Also the statement S = 100.D0 was removed
C           from the top of the Newton iteration in DDASTP.
C           The subordinate routines had revision dates changed to
C           this date.
C   890517  The revision date syntax was replaced with the revision
C           history syntax.  Also the "DECK" comment was added to
C           the top of all subroutines.  These changes are consistent
C           with new SLATEC guidelines.
C           The subordinate routines had revision dates changed to
C           this date.  No code changes.
C   891013  Code changes made.
C           Removed all occurrances of FLOAT or DBLE.  All operations
C           are now performed with "mixed-mode" arithmetic.
C           Also, specific function names were replaced with generic
C           function names to be consistent with new SLATEC guidelines.
C           In particular:
C              Replaced DSQRT with SQRT everywhere.
C              Replaced DABS with ABS everywhere.
C              Replaced DMIN1 with MIN everywhere.
C              Replaced MIN0 with MIN everywhere.
C              Replaced DMAX1 with MAX everywhere.
C              Replaced MAX0 with MAX everywhere.
C              Replaced DSIGN with SIGN everywhere.
C           Also replaced REVISION DATE with REVISION HISTORY in all
C           subordinate routines.
C  901004  Miscellaneous changes to prologue to complete conversion
C          to SLATEC 4.0 format.  No code changes.  (F.N.Fritsch)
C  901009  Corrected GAMS classification code and converted subsidiary
C          routines to 4.0 format.  No code changes.  (F.N.Fritsch)
C  901010  Converted XERRWV calls to XERMSG calls.  (R.Clemens,AFWL)
C  901019  Code changes made.
C          Merged SLATEC 4.0 changes with previous changes made
C          by C. Ulrich.  Below is a history of the changes made by
C          C. Ulrich. (Changes in subsidiary routines are implied
C          by this history)
C          891228  Bug was found and repaired inside the DDASSL
C                  and DDAINI routines.  DDAINI was incorrectly
C                  returning the initial T with Y and YPRIME
C                  computed at T+H.  The routine now returns T+H
C                  rather than the initial T.
C                  Cosmetic changes made to DDASTP.
C          900904  Three modifications were made to fix a bug (inside
C                  DDASSL) re interpolation for continuation calls and
C                  cases where TN is very close to TSTOP:
C
C                  1) In testing for whether H is too large, just
C                     compare H to (TSTOP - TN), rather than
C                     (TSTOP - TN) * (1-4*UROUND), and set H to
C                     TSTOP - TN.  This will force DDASTP to step
C                     exactly to TSTOP under certain situations
C                     (i.e. when H returned from DDASTP would otherwise
C                     take TN beyond TSTOP).
C
C                  2) Inside the DDASTP loop, interpolate exactly to
C                     TSTOP if TN is very close to TSTOP (rather than
C                     interpolating to within roundoff of TSTOP).
C
C                  3) Modified IDID description for IDID = 2 to say
C                     the solution is returned by stepping exactly to
C                     TSTOP, rather than TOUT.  (In some cases the
C                     solution is actually obtained by extrapolating
C                     over a distance near unit roundoff to TSTOP,
C                     but this small distance is deemed acceptable in
C                     these circumstances.)
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue, removed unreferenced labels,
C           and improved XERMSG calls.  (FNF)
C   901030  Added ERROR MESSAGES section and reworked other sections to
C           be of more uniform format.  (FNF)
C   910624  Fixed minor bug related to HMAX (five lines ending in
C           statement 526 in DDASSL).   (LRP)
C
C***END PROLOGUE  DDASSL
C
C**End
C
C     Declare arguments.
C
      INTEGER  NEQ, INFO(15), IDID, LRW, IWORK(*), LIW, IPAR(*)
      DOUBLE PRECISION
     *   T, Y(*), YPRIME(*), TOUT, RTOL(*), ATOL(*), RWORK(*),
     *   RPAR(*), SIGN77
      EXTERNAL  RES, JAC, SIGN77
C
C     Declare externals.
C
      EXTERNAL  D1MACH, DDAINI, DDANRM, DDASTP, DDATRP, DDAWTS, XERMSG
      DOUBLE PRECISION  D1MACH, DDANRM
C
C     Declare local variables.
C
      INTEGER  I, ITEMP, LALPHA, LBETA, LCJ, LCJOLD, LCTF, LDELTA,
     *   LENIW, LENPD,LENRW,LE, LETF, LGAMMA, LH, LHMAX, LHOLD, LIPVT,
     *   LJCALC, LK, LKOLD, LIWM, LML, LMTYPE, LMU, LMXORD, LNJE, LNPD,
     *   LNRE, LNS, LNST, LNSTL, LPD, LPHASE, LPHI, LPSI, LROUND, LS,
     *   LSIGMA,LTN,LTSTOP, LWM, LWT, MBAND, MSAVE, MXORD, NPD, NTEMP,
     *   NZFLG
      DOUBLE PRECISION
     *   ATOLI, H, HMAX, HMIN, HO, R, RH, RTOLI, TDIST, TN, TNEXT,
     *   TSTOP, UROUND, YPNORM
      LOGICAL  DONE
C       Auxiliary variables for conversion of values to be included in
C       error messages.
      CHARACTER*8  XERN1, XERN2
      CHARACTER*16 XERN3, XERN4
C
C     SET POINTERS INTO IWORK
      PARAMETER (LML=1, LMU=2, LMXORD=3, LMTYPE=4, LNST=11,
     *  LNRE=12, LNJE=13, LETF=14, LCTF=15, LNPD=16,
     *  LIPVT=21, LJCALC=5, LPHASE=6, LK=7, LKOLD=8,
     *  LNS=9, LNSTL=10, LIWM=1)
C
C     SET RELATIVE OFFSET INTO RWORK
      PARAMETER (NPD=1)
C
C     SET POINTERS INTO RWORK
      PARAMETER (LTSTOP=1, LHMAX=2, LH=3, LTN=4,
     *  LCJ=5, LCJOLD=6, LHOLD=7, LS=8, LROUND=9,
     *  LALPHA=11, LBETA=17, LGAMMA=23,
     *  LPSI=29, LSIGMA=35, LDELTA=41)
C
C***FIRST EXECUTABLE STATEMENT  DDASSL
      IF(INFO(1).NE.0)GO TO 100
C
C----------------------------------------------------------------------
C     THIS BLOCK IS EXECUTED FOR THE INITIAL CALL ONLY.
C     IT CONTAINS CHECKING OF INPUTS AND INITIALIZATIONS.
C----------------------------------------------------------------------
C
C     FIRST CHECK INFO ARRAY TO MAKE SURE ALL ELEMENTS OF INFO
C     ARE EITHER ZERO OR ONE.
      DO 10 I=2,11
         IF(INFO(I).NE.0.AND.INFO(I).NE.1)GO TO 701
10       CONTINUE
C
      IF(NEQ.LE.0)GO TO 702
C
C     CHECK AND COMPUTE MAXIMUM ORDER
      MXORD=5
      IF(INFO(9).EQ.0)GO TO 20
         MXORD=IWORK(LMXORD)
         IF(MXORD.LT.1.OR.MXORD.GT.5)GO TO 703
20       IWORK(LMXORD)=MXORD
C
C     COMPUTE MTYPE,LENPD,LENRW.CHECK ML AND MU.
      IF(INFO(6).NE.0)GO TO 40
         LENPD=NEQ**2
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
         IF(INFO(5).NE.0)GO TO 30
            IWORK(LMTYPE)=2
            GO TO 60
30          IWORK(LMTYPE)=1
            GO TO 60
40    IF(IWORK(LML).LT.0.OR.IWORK(LML).GE.NEQ)GO TO 717
      IF(IWORK(LMU).LT.0.OR.IWORK(LMU).GE.NEQ)GO TO 718
      LENPD=(2*IWORK(LML)+IWORK(LMU)+1)*NEQ
      IF(INFO(5).NE.0)GO TO 50
         IWORK(LMTYPE)=5
         MBAND=IWORK(LML)+IWORK(LMU)+1
         MSAVE=(NEQ/MBAND)+1
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD+2*MSAVE
         GO TO 60
50       IWORK(LMTYPE)=4
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
C
C     CHECK LENGTHS OF RWORK AND IWORK
60    LENIW=20+NEQ
      IWORK(LNPD)=LENPD
      IF(LRW.LT.LENRW)GO TO 704
      IF(LIW.LT.LENIW)GO TO 705
C
C     CHECK TO SEE THAT TOUT IS DIFFERENT FROM T
      IF(TOUT .EQ. T)GO TO 719
C
C     CHECK HMAX
      IF(INFO(7).EQ.0)GO TO 70
         HMAX=RWORK(LHMAX)
         IF(HMAX.LE.0.0D0)GO TO 710
70    CONTINUE
C
C     INITIALIZE COUNTERS
      IWORK(LNST)=0
      IWORK(LNRE)=0
      IWORK(LNJE)=0
C
      IWORK(LNSTL)=0
      IDID=1
      GO TO 200
C
C----------------------------------------------------------------------
C     THIS BLOCK IS FOR CONTINUATION CALLS
C     ONLY. HERE WE CHECK INFO(1),AND IF THE
C     LAST STEP WAS INTERRUPTED WE CHECK WHETHER
C     APPROPRIATE ACTION WAS TAKEN.
C----------------------------------------------------------------------
C
100   CONTINUE
      IF(INFO(1).EQ.1)GO TO 110
      IF(INFO(1).NE.-1)GO TO 701
C
C     IF WE ARE HERE, THE LAST STEP WAS INTERRUPTED
C     BY AN ERROR CONDITION FROM DDASTP,AND
C     APPROPRIATE ACTION WAS NOT TAKEN. THIS
C     IS A FATAL ERROR.
      WRITE (XERN1, '(I8)') IDID
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'THE LAST STEP TERMINATED WITH A NEGATIVE VALUE OF IDID = ' //
     *   XERN1 // ' AND NO APPROPRIATE ACTION WAS TAKEN.  ' //
     *   'RUN TERMINATED', -998, 2)
      RETURN
110   CONTINUE
      IWORK(LNSTL)=IWORK(LNST)
C
C----------------------------------------------------------------------
C     THIS BLOCK IS EXECUTED ON ALL CALLS.
C     THE ERROR TOLERANCE PARAMETERS ARE
C     CHECKED, AND THE WORK ARRAY POINTERS
C     ARE SET.
C----------------------------------------------------------------------
C
200   CONTINUE
C     CHECK RTOL,ATOL
      NZFLG=0
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 210 I=1,NEQ
         IF(INFO(2).EQ.1)RTOLI=RTOL(I)
         IF(INFO(2).EQ.1)ATOLI=ATOL(I)
         IF(RTOLI.GT.0.0D0.OR.ATOLI.GT.0.0D0)NZFLG=1
         IF(RTOLI.LT.0.0D0)GO TO 706
         IF(ATOLI.LT.0.0D0)GO TO 707
210      CONTINUE
      IF(NZFLG.EQ.0)GO TO 708
C
C     SET UP RWORK STORAGE.IWORK STORAGE IS FIXED
C     IN DATA STATEMENT.
      LE=LDELTA+NEQ
      LWT=LE+NEQ
      LPHI=LWT+NEQ
      LPD=LPHI+(IWORK(LMXORD)+1)*NEQ
      LWM=LPD
      NTEMP=NPD+IWORK(LNPD)
      IF(INFO(1).EQ.1)GO TO 400
C
C----------------------------------------------------------------------
C     THIS BLOCK IS EXECUTED ON THE INITIAL CALL
C     ONLY. SET THE INITIAL STEP SIZE, AND
C     THE ERROR WEIGHT VECTOR, AND PHI.
C     COMPUTE INITIAL YPRIME, IF NECESSARY.
C----------------------------------------------------------------------
C
      TN=T
      IDID=1
C
C     SET ERROR WEIGHT VECTOR WT
      CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      DO 305 I = 1,NEQ
         IF(RWORK(LWT+I-1).LE.0.0D0) GO TO 713
305      CONTINUE
C
C     COMPUTE UNIT ROUNDOFF AND HMIN
      UROUND = D1MACH(4)
      RWORK(LROUND) = UROUND
      HMIN = 4.0D0*UROUND*MAX(ABS(T),ABS(TOUT))
C
C     CHECK INITIAL INTERVAL TO SEE THAT IT IS LONG ENOUGH
      TDIST = ABS(TOUT - T)
      IF(TDIST .LT. HMIN) GO TO 714
C
C     CHECK HO, IF THIS WAS INPUT
      IF (INFO(8) .EQ. 0) GO TO 310
         HO = RWORK(LH)
         IF ((TOUT - T)*HO .LT. 0.0D0) GO TO 711
         IF (HO .EQ. 0.0D0) GO TO 712
         GO TO 320
310    CONTINUE
C
C     COMPUTE INITIAL STEPSIZE, TO BE USED BY EITHER
C     DDASTP OR DDAINI, DEPENDING ON INFO(11)
      HO = 0.001D0*TDIST
      YPNORM = DDANRM(NEQ,YPRIME,RWORK(LWT),RPAR,IPAR)
      IF (YPNORM .GT. 0.5D0/HO) HO = 0.5D0/YPNORM
      HO = SIGN77(HO,TOUT-T)
C     ADJUST HO IF NECESSARY TO MEET HMAX BOUND
320   IF (INFO(7) .EQ. 0) GO TO 330
         RH = ABS(HO)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) HO = HO/RH
C     COMPUTE TSTOP, IF APPLICABLE
330   IF (INFO(4) .EQ. 0) GO TO 340
         TSTOP = RWORK(LTSTOP)
         IF ((TSTOP - T)*HO .LT. 0.0D0) GO TO 715
         IF ((T + HO - TSTOP)*HO .GT. 0.0D0) HO = TSTOP - T
         IF ((TSTOP - TOUT)*HO .LT. 0.0D0) GO TO 709
C
C     COMPUTE INITIAL DERIVATIVE, UPDATING TN AND Y, IF APPLICABLE
340   IF (INFO(11) .EQ. 0) GO TO 350
      CALL DDAINI(TN,Y,YPRIME,NEQ,
     *  RES,JAC,HO,RWORK(LWT),IDID,RPAR,IPAR,
     *  RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
     *  RWORK(LWM),IWORK(LIWM),HMIN,RWORK(LROUND),
     *  INFO(10),NTEMP)
      IF (IDID .LT. 0) GO TO 390
C
C     LOAD H WITH HO.  STORE H IN RWORK(LH)
350   H = HO
      RWORK(LH) = H
C
C     LOAD Y AND H*YPRIME INTO PHI(*,1) AND PHI(*,2)
      ITEMP = LPHI + NEQ
      DO 370 I = 1,NEQ
         RWORK(LPHI + I - 1) = Y(I)
370      RWORK(ITEMP + I - 1) = H*YPRIME(I)
C
390   GO TO 500
C
C------------------------------------------------------
C     THIS BLOCK IS FOR CONTINUATION CALLS ONLY. ITS
C     PURPOSE IS TO CHECK STOP CONDITIONS BEFORE
C     TAKING A STEP.
C     ADJUST H IF NECESSARY TO MEET HMAX BOUND
C------------------------------------------------------
C
400   CONTINUE
      UROUND=RWORK(LROUND)
      DONE = .FALSE.
      TN=RWORK(LTN)
      H=RWORK(LH)
      IF(INFO(7) .EQ. 0) GO TO 410
         RH = ABS(H)/RWORK(LHMAX)
         IF(RH .GT. 1.0D0) H = H/RH
410   CONTINUE
      IF(T .EQ. TOUT) GO TO 719
      IF((T - TOUT)*H .GT. 0.0D0) GO TO 711
      IF(INFO(4) .EQ. 1) GO TO 430
      IF(INFO(3) .EQ. 1) GO TO 420
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 490
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
420   IF((TN-T)*H .LE. 0.0D0) GO TO 490
      IF((TN - TOUT)*H .GT. 0.0D0) GO TO 425
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
      IF((TN-TSTOP)*H.GT.0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H.LT.0.0D0)GO TO 709
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 450
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *   RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
440   TSTOP = RWORK(LTSTOP)
      IF((TN-TSTOP)*H .GT. 0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H .LT. 0.0D0) GO TO 709
      IF((TN-T)*H .LE. 0.0D0) GO TO 450
      IF((TN - TOUT)*H .GT. 0.0D0) GO TO 445
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
C     CHECK WHETHER WE ARE WITHIN ROUNDOFF OF TSTOP
      IF(ABS(TN-TSTOP).GT.100.0D0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 460
      CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      DONE = .TRUE.
      GO TO 490
460   TNEXT=TN+H
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 490
      H=TSTOP-TN
      RWORK(LH)=H
C
490   IF (DONE) GO TO 580
C
C------------------------------------------------------
C     THE NEXT BLOCK CONTAINS THE CALL TO THE
C     ONE-STEP INTEGRATOR DDASTP.
C     THIS IS A LOOPING POINT FOR THE INTEGRATION STEPS.
C     CHECK FOR TOO MANY STEPS.
C     UPDATE WT.
C     CHECK FOR TOO MUCH ACCURACY REQUESTED.
C     COMPUTE MINIMUM STEPSIZE.
C------------------------------------------------------
C
500   CONTINUE
C     CHECK FOR FAILURE TO COMPUTE INITIAL YPRIME
      IF (IDID .EQ. -12) GO TO 527
C
C     CHECK FOR TOO MANY STEPS
      IF((IWORK(LNST)-IWORK(LNSTL)).LT.500)
     *   GO TO 510
           IDID=-1
           GO TO 527
C
C     UPDATE WT
510   CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,RWORK(LPHI),
     *  RWORK(LWT),RPAR,IPAR)
      DO 520 I=1,NEQ
         IF(RWORK(I+LWT-1).GT.0.0D0)GO TO 520
           IDID=-3
           GO TO 527
520   CONTINUE
C
C     TEST FOR TOO MUCH ACCURACY REQUESTED.
      R=DDANRM(NEQ,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)*
     *   100.0D0*UROUND
      IF(R.LE.1.0D0)GO TO 525
C     MULTIPLY RTOL AND ATOL BY R AND RETURN
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
C     COMPUTE MINIMUM STEPSIZE
      HMIN=4.0D0*UROUND*MAX(ABS(TN),ABS(TOUT))
C
C     TEST H VS. HMAX
      IF (INFO(7) .EQ. 0) GO TO 526
         RH = ABS(H)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) H = H/RH
526   CONTINUE
C
      CALL DDASTP(TN,Y,YPRIME,NEQ,
     *   RES,JAC,H,RWORK(LWT),INFO(1),IDID,RPAR,IPAR,
     *   RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
     *   RWORK(LWM),IWORK(LIWM),
     *   RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),
     *   RWORK(LPSI),RWORK(LSIGMA),
     *   RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),
     *   RWORK(LS),HMIN,RWORK(LROUND),
     *   IWORK(LPHASE),IWORK(LJCALC),IWORK(LK),
     *   IWORK(LKOLD),IWORK(LNS),INFO(10),NTEMP)
527   IF(IDID.LT.0)GO TO 600
C
C-------------------------------------------------------
C     THIS BLOCK HANDLES THE CASE OF A SUCCESSFUL RETURN
C     FROM DDASTP (IDID=1).  TEST FOR STOP CONDITIONS.
C-------------------------------------------------------
C
      IF(INFO(4).NE.0)GO TO 540
           IF(INFO(3).NE.0)GO TO 530
             IF((TN-TOUT)*H.LT.0.0D0)GO TO 500
             CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
530          IF((TN-TOUT)*H.GE.0.0D0)GO TO 535
             T=TN
             IDID=1
             GO TO 580
535          CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
540   IF(INFO(3).NE.0)GO TO 550
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 542
         CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
         T=TOUT
         IDID=3
         GO TO 580
542   IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 545
      TNEXT=TN+H
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 500
      H=TSTOP-TN
      GO TO 500
545   CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,
     *  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      GO TO 580
550   IF((TN-TOUT)*H.GE.0.0D0)GO TO 555
      IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*(ABS(TN)+ABS(H)))GO TO 552
      T=TN
      IDID=1
      GO TO 580
552   CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,
     *  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      GO TO 580
555   CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *   IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID=3
      GO TO 580
C
C-------------------------------------------------------
C     ALL SUCCESSFUL RETURNS FROM DDASSL ARE MADE FROM
C     THIS BLOCK.
C-------------------------------------------------------
C
580   CONTINUE
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C
C----------------------------------------------------------------------
C     THIS BLOCK HANDLES ALL UNSUCCESSFUL
C     RETURNS OTHER THAN FOR ILLEGAL INPUT.
C----------------------------------------------------------------------
C
600   CONTINUE
      ITEMP=-IDID
      GO TO (610,620,630,690,690,640,650,660,670,675,
     *  680,685), ITEMP
C
C     THE MAXIMUM NUMBER OF STEPS WAS TAKEN BEFORE
C     REACHING TOUT
610   WRITE (XERN3, '(1P,D15.6)') TN
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT CURRENT T = ' // XERN3 // ' 500 STEPS TAKEN ON THIS ' //
     *   'CALL BEFORE REACHING TOUT', IDID, 1)
      GO TO 690
C
C     TOO MUCH ACCURACY FOR MACHINE PRECISION
620   WRITE (XERN3, '(1P,D15.6)') TN
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' TOO MUCH ACCURACY REQUESTED FOR ' //
     *   'PRECISION OF MACHINE. RTOL AND ATOL WERE INCREASED TO ' //
     *   'APPROPRIATE VALUES', IDID, 1)
      GO TO 690
C
C     WT(I) .LE. 0.0 FOR SOME I (NOT AT START OF PROBLEM)
630   WRITE (XERN3, '(1P,D15.6)') TN
      CALL XERMSG ('SLATEC', 'DDASSL',
     * 'AT T = ' // XERN3 // ' SOME ELEMENT OF WT HAS BECOME .LE. ' //
     *   '0.0', IDID, 1)
      GO TO 690
C
C     ERROR TEST FAILED REPEATEDLY OR WITH H=HMIN
640   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN',
     *   IDID, 1)
      GO TO 690
C
C     CORRECTOR CONVERGENCE FAILED REPEATEDLY OR WITH H=HMIN
650   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE CORRECTOR FAILED TO CONVERGE REPEATEDLY OR WITH ' //
     *   'ABS(H)=HMIN', IDID, 1)
      GO TO 690
C
C     THE ITERATION MATRIX IS SINGULAR
660   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE ITERATION MATRIX IS SINGULAR', IDID, 1)
      GO TO 690
C
C     CORRECTOR FAILURE PRECEEDED BY ERROR TEST FAILURES.
670   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE CORRECTOR COULD NOT CONVERGE.  ALSO, THE ERROR TEST ' //
     *   'FAILED REPEATEDLY.', IDID, 1)
      GO TO 690
C
C     CORRECTOR FAILURE BECAUSE IRES = -1
675   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE CORRECTOR COULD NOT CONVERGE BECAUSE IRES WAS EQUAL ' //
     *   'TO MINUS ONE', IDID, 1)
      GO TO 690
C
C     FAILURE BECAUSE IRES = -2
680   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' IRES WAS EQUAL TO MINUS TWO', IDID, 1)
      GO TO 690
C
C     FAILED TO COMPUTE INITIAL YPRIME
685   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') HO
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE INITIAL YPRIME COULD NOT BE COMPUTED', IDID, 1)
      GO TO 690
C
690   CONTINUE
      INFO(1)=-1
      T=TN
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C
C----------------------------------------------------------------------
C     THIS BLOCK HANDLES ALL ERROR RETURNS DUE
C     TO ILLEGAL INPUT, AS DETECTED BEFORE CALLING
C     DDASTP. FIRST THE ERROR MESSAGE ROUTINE IS
C     CALLED. IF THIS HAPPENS TWICE IN
C     SUCCESSION, EXECUTION IS TERMINATED
C
C----------------------------------------------------------------------
701   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR ONE', 1, 1)
      GO TO 750
C
702   WRITE (XERN1, '(I8)') NEQ
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'NEQ = ' // XERN1 // ' .LE. 0', 2, 1)
      GO TO 750
C
703   WRITE (XERN1, '(I8)') MXORD
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'MAXORD = ' // XERN1 // ' NOT IN RANGE', 3, 1)
      GO TO 750
C
704   WRITE (XERN1, '(I8)') LENRW
      WRITE (XERN2, '(I8)') LRW
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'RWORK LENGTH NEEDED, LENRW = ' // XERN1 //
     *   ', EXCEEDS LRW = ' // XERN2, 4, 1)
      GO TO 750
C
705   WRITE (XERN1, '(I8)') LENIW
      WRITE (XERN2, '(I8)') LIW
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'IWORK LENGTH NEEDED, LENIW = ' // XERN1 //
     *   ', EXCEEDS LIW = ' // XERN2, 5, 1)
      GO TO 750
C
706   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'SOME ELEMENT OF RTOL IS .LT. 0', 6, 1)
      GO TO 750
C
707   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'SOME ELEMENT OF ATOL IS .LT. 0', 7, 1)
      GO TO 750
C
708   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'ALL ELEMENTS OF RTOL AND ATOL ARE ZERO', 8, 1)
      GO TO 750
C
709   WRITE (XERN3, '(1P,D15.6)') TSTOP
      WRITE (XERN4, '(1P,D15.6)') TOUT
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'INFO(4) = 1 AND TSTOP = ' // XERN3 // ' BEHIND TOUT = ' //
     *   XERN4, 9, 1)
      GO TO 750
C
710   WRITE (XERN3, '(1P,D15.6)') HMAX
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'HMAX = ' // XERN3 // ' .LT. 0.0', 10, 1)
      GO TO 750
C
711   WRITE (XERN3, '(1P,D15.6)') TOUT
      WRITE (XERN4, '(1P,D15.6)') T
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'TOUT = ' // XERN3 // ' BEHIND T = ' // XERN4, 11, 1)
      GO TO 750
C
712   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'INFO(8)=1 AND H0=0.0', 12, 1)
      GO TO 750
C
713   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'SOME ELEMENT OF WT IS .LE. 0.0', 13, 1)
      GO TO 750
C
714   WRITE (XERN3, '(1P,D15.6)') TOUT
      WRITE (XERN4, '(1P,D15.6)') T
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'TOUT = ' // XERN3 // ' TOO CLOSE TO T = ' // XERN4 //
     *   ' TO START INTEGRATION', 14, 1)
      GO TO 750
C
715   WRITE (XERN3, '(1P,D15.6)') TSTOP
      WRITE (XERN4, '(1P,D15.6)') T
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'INFO(4)=1 AND TSTOP = ' // XERN3 // ' BEHIND T = ' // XERN4,
     *   15, 1)
      GO TO 750
C
717   WRITE (XERN1, '(I8)') IWORK(LML)
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'ML = ' // XERN1 // ' ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ',
     *   17, 1)
      GO TO 750
C
718   WRITE (XERN1, '(I8)') IWORK(LMU)
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'MU = ' // XERN1 // ' ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ',
     *   18, 1)
      GO TO 750
C
719   WRITE (XERN3, '(1P,D15.6)') TOUT
      CALL XERMSG ('SLATEC', 'DDASSL',
     *  'TOUT = T = ' // XERN3, 19, 1)
      GO TO 750
C
750   IDID=-33
      IF(INFO(1).EQ.-1) THEN
         CALL XERMSG ('SLATEC', 'DDASSL',
     *      'REPEATED OCCURRENCES OF ILLEGAL INPUT$$' //
     *      'RUN TERMINATED. APPARENT INFINITE LOOP', -999, 2)
      ENDIF
C
      INFO(1)=-1
      RETURN
C----------END OF SUBROUTINE DDASSL------------------------------------
      END
      SUBROUTINE DDASTP (X, Y, YPRIME, NEQ, RES, JAC, H, WT, JSTART,
     +   IDID, RPAR, IPAR, PHI, DELTA, E, WM, IWM, ALPHA, BETA, GAMMA,
     +   PSI, SIGMA, CJ, CJOLD, HOLD, S, HMIN, UROUND, IPHASE, JCALC,
     +   K, KOLD, NS, NONNEG, NTEMP)
C***BEGIN PROLOGUE  DDASTP
C***SUBSIDIARY
C***PURPOSE  Perform one step of the DDASSL integration.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDASTP-S, DDASTP-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***DESCRIPTION
C----------------------------------------------------------------------
C     DDASTP SOLVES A SYSTEM OF DIFFERENTIAL/
C     ALGEBRAIC EQUATIONS OF THE FORM
C     G(X,Y,YPRIME) = 0,  FOR ONE STEP (NORMALLY
C     FROM X TO X+H).
C
C     THE METHODS USED ARE MODIFIED DIVIDED
C     DIFFERENCE,FIXED LEADING COEFFICIENT
C     FORMS OF BACKWARD DIFFERENTIATION
C     FORMULAS. THE CODE ADJUSTS THE STEPSIZE
C     AND ORDER TO CONTROL THE LOCAL ERROR PER
C     STEP.
C
C
C     THE PARAMETERS REPRESENT
C     X  --        INDEPENDENT VARIABLE
C     Y  --        SOLUTION VECTOR AT X
C     YPRIME --    DERIVATIVE OF SOLUTION VECTOR
C                  AFTER SUCCESSFUL STEP
C     NEQ --       NUMBER OF EQUATIONS TO BE INTEGRATED
C     RES --       EXTERNAL USER-SUPPLIED SUBROUTINE
C                  TO EVALUATE THE RESIDUAL.  THE CALL IS
C                  CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
C                  X,Y,YPRIME ARE INPUT.  DELTA IS OUTPUT.
C                  ON INPUT, IRES=0.  RES SHOULD ALTER IRES ONLY
C                  IF IT ENCOUNTERS AN ILLEGAL VALUE OF Y OR A
C                  STOP CONDITION.  SET IRES=-1 IF AN INPUT VALUE
C                  OF Y IS ILLEGAL, AND DDASTP WILL TRY TO SOLVE
C                  THE PROBLEM WITHOUT GETTING IRES = -1.  IF
C                  IRES=-2, DDASTP RETURNS CONTROL TO THE CALLING
C                  PROGRAM WITH IDID = -11.
C     JAC --       EXTERNAL USER-SUPPLIED ROUTINE TO EVALUATE
C                  THE ITERATION MATRIX (THIS IS OPTIONAL)
C                  THE CALL IS OF THE FORM
C                  CALL JAC(X,Y,YPRIME,PD,CJ,RPAR,IPAR)
C                  PD IS THE MATRIX OF PARTIAL DERIVATIVES,
C                  PD=DG/DY+CJ*DG/DYPRIME
C     H --         APPROPRIATE STEP SIZE FOR NEXT STEP.
C                  NORMALLY DETERMINED BY THE CODE
C     WT --        VECTOR OF WEIGHTS FOR ERROR CRITERION.
C     JSTART --    INTEGER VARIABLE SET 0 FOR
C                  FIRST STEP, 1 OTHERWISE.
C     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS:
C                  IDID= 1 -- THE STEP WAS COMPLETED SUCCESSFULLY
C                  IDID=-6 -- THE ERROR TEST FAILED REPEATEDLY
C                  IDID=-7 -- THE CORRECTOR COULD NOT CONVERGE
C                  IDID=-8 -- THE ITERATION MATRIX IS SINGULAR
C                  IDID=-9 -- THE CORRECTOR COULD NOT CONVERGE.
C                             THERE WERE REPEATED ERROR TEST
C                             FAILURES ON THIS STEP.
C                  IDID=-10-- THE CORRECTOR COULD NOT CONVERGE
C                             BECAUSE IRES WAS EQUAL TO MINUS ONE
C                  IDID=-11-- IRES EQUAL TO -2 WAS ENCOUNTERED,
C                             AND CONTROL IS BEING RETURNED TO
C                             THE CALLING PROGRAM
C     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS THAT
C                  ARE USED FOR COMMUNICATION BETWEEN THE
C                  CALLING PROGRAM AND EXTERNAL USER ROUTINES
C                  THEY ARE NOT ALTERED BY DDASTP
C     PHI --       ARRAY OF DIVIDED DIFFERENCES USED BY
C                  DDASTP. THE LENGTH IS NEQ*(K+1),WHERE
C                  K IS THE MAXIMUM ORDER
C     DELTA,E --   WORK VECTORS FOR DDASTP OF LENGTH NEQ
C     WM,IWM --    REAL AND INTEGER ARRAYS STORING
C                  MATRIX INFORMATION SUCH AS THE MATRIX
C                  OF PARTIAL DERIVATIVES,PERMUTATION
C                  VECTOR,AND VARIOUS OTHER INFORMATION.
C
C     THE OTHER PARAMETERS ARE INFORMATION
C     WHICH IS NEEDED INTERNALLY BY DDASTP TO
C     CONTINUE FROM STEP TO STEP.
C
C---------------------------------------------------------------------
C***ROUTINES CALLED  DDAJAC, DDANRM, DDASLV, DDATRP
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  DDASTP
C
      INTEGER  NEQ, JSTART, IDID, IPAR(*), IWM(*), IPHASE, JCALC, K,
     *   KOLD, NS, NONNEG, NTEMP
      DOUBLE PRECISION
     *   X, Y(*), YPRIME(*), H, WT(*), RPAR(*), PHI(NEQ,*), DELTA(*),
     * E(*), WM(*), ALPHA(*), BETA(*), GAMMA(*), PSI(*), SIGMA(*), CJ,
     *   CJOLD, HOLD, S, HMIN, UROUND
      EXTERNAL  RES, JAC
C
      EXTERNAL  DDAJAC, DDANRM, DDASLV, DDATRP
      DOUBLE PRECISION  DDANRM
C
      INTEGER  I, IER, IRES, J, J1, KDIFF, KM1, KNEW, KP1, KP2, LCTF,
     *   LETF, LMXORD, LNJE, LNRE, LNST, M, MAXIT, NCF, NEF, NSF, NSP1
      DOUBLE PRECISION
     *   ALPHA0, ALPHAS, CJLAST, CK, DELNRM, ENORM, ERK, ERKM1,
     *   ERKM2, ERKP1, ERR, EST, HNEW, OLDNRM, PNORM, R, RATE, TEMP1,
     *   TEMP2, TERK, TERKM1, TERKM2, TERKP1, XOLD, XRATE
      LOGICAL  CONVGD
C
      PARAMETER (LMXORD=3)
      PARAMETER (LNST=11)
      PARAMETER (LNRE=12)
      PARAMETER (LNJE=13)
      PARAMETER (LETF=14)
      PARAMETER (LCTF=15)
C
      DATA MAXIT/4/
      DATA XRATE/0.25D0/
C
C
C
C
C
C----------------------------------------------------------------------
C     BLOCK 1.
C     INITIALIZE. ON THE FIRST CALL,SET
C     THE ORDER TO 1 AND INITIALIZE
C     OTHER VARIABLES.
C----------------------------------------------------------------------
C
C     INITIALIZATIONS FOR ALL CALLS
C***FIRST EXECUTABLE STATEMENT  DDASTP
      IDID=1
      XOLD=X
      NCF=0
      NSF=0
      NEF=0
      IF(JSTART .NE. 0) GO TO 120
C
C     IF THIS IS THE FIRST STEP,PERFORM
C     OTHER INITIALIZATIONS
      IWM(LETF) = 0
      IWM(LCTF) = 0
      K=1
      KOLD=0
      HOLD=0.0D0
      JSTART=1
      PSI(1)=H
      CJOLD = 1.0D0/H
      CJ = CJOLD
      S = 100.D0
      JCALC = -1
      DELNRM=1.0D0
      IPHASE = 0
      NS=0
120   CONTINUE
C
C
C
C
C
C--------------------------------------------------------------------
C     BLOCK 2
C     COMPUTE COEFFICIENTS OF FORMULAS FOR
C     THIS STEP.
C--------------------------------------------------------------------
200   CONTINUE
      KP1=K+1
      KP2=K+2
      KM1=K-1
      XOLD=X
      IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0
      NS=MIN(NS+1,KOLD+2)
      NSP1=NS+1
      IF(KP1 .LT. NS)GO TO 230
C
      BETA(1)=1.0D0
      ALPHA(1)=1.0D0
      TEMP1=H
      GAMMA(1)=0.0D0
      SIGMA(1)=1.0D0
      DO 210 I=2,KP1
         TEMP2=PSI(I-1)
         PSI(I-1)=TEMP1
         BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
         TEMP1=TEMP2+H
         ALPHA(I)=H/TEMP1
         SIGMA(I)=(I-1)*SIGMA(I-1)*ALPHA(I)
         GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
210      CONTINUE
      PSI(KP1)=TEMP1
230   CONTINUE
C
C     COMPUTE ALPHAS, ALPHA0
      ALPHAS = 0.0D0
      ALPHA0 = 0.0D0
      DO 240 I = 1,K
        ALPHAS = ALPHAS - 1.0D0/I
        ALPHA0 = ALPHA0 - ALPHA(I)
240     CONTINUE
C
C     COMPUTE LEADING COEFFICIENT CJ
      CJLAST = CJ
      CJ = -ALPHAS/H
C
C     COMPUTE VARIABLE STEPSIZE ERROR COEFFICIENT CK
      CK = ABS(ALPHA(KP1) + ALPHAS - ALPHA0)
      CK = MAX(CK,ALPHA(KP1))
C
C     DECIDE WHETHER NEW JACOBIAN IS NEEDED
      TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)
      TEMP2 = 1.0D0/TEMP1
      IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
      IF (CJ .NE. CJLAST) S = 100.D0
C
C     CHANGE PHI TO PHI STAR
      IF(KP1 .LT. NSP1) GO TO 280
      DO 270 J=NSP1,KP1
         DO 260 I=1,NEQ
260         PHI(I,J)=BETA(J)*PHI(I,J)
270      CONTINUE
280   CONTINUE
C
C     UPDATE TIME
      X=X+H
C
C
C
C
C
C--------------------------------------------------------------------
C     BLOCK 3
C     PREDICT THE SOLUTION AND DERIVATIVE,
C     AND SOLVE THE CORRECTOR EQUATION
C-------------------------------------------------------------------
C
C     FIRST,PREDICT THE SOLUTION AND DERIVATIVE
300   CONTINUE
      DO 310 I=1,NEQ
         Y(I)=PHI(I,1)
310      YPRIME(I)=0.0D0
      DO 330 J=2,KP1
         DO 320 I=1,NEQ
            Y(I)=Y(I)+PHI(I,J)
320         YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
330   CONTINUE
      PNORM = DDANRM (NEQ,Y,WT,RPAR,IPAR)
C
C
C
C     SOLVE THE CORRECTOR EQUATION USING A
C     MODIFIED NEWTON SCHEME.
      CONVGD= .TRUE.
      M=0
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
C
C
C     IF INDICATED,REEVALUATE THE
C     ITERATION MATRIX PD = DG/DY + CJ*DG/DYPRIME
C     (WHERE G(X,Y,YPRIME)=0). SET
C     JCALC TO 0 AS AN INDICATOR THAT
C     THIS HAS BEEN DONE.
      IF(JCALC .NE. -1)GO TO 340
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      CALL DDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,
     * IER,WT,E,WM,IWM,RES,IRES,UROUND,JAC,RPAR,
     * IPAR,NTEMP)
      CJOLD=CJ
      S = 100.D0
      IF (IRES .LT. 0) GO TO 380
      IF(IER .NE. 0)GO TO 380
      NSF=0
C
C
C     INITIALIZE THE ERROR ACCUMULATION VECTOR E.
340   CONTINUE
      DO 345 I=1,NEQ
345      E(I)=0.0D0
C
C
C     CORRECTOR LOOP.
350   CONTINUE
C
C     MULTIPLY RESIDUAL BY TEMP1 TO ACCELERATE CONVERGENCE
      TEMP1 = 2.0D0/(1.0D0 + CJ/CJOLD)
      DO 355 I = 1,NEQ
355     DELTA(I) = DELTA(I) * TEMP1
C
C     COMPUTE A NEW ITERATE (BACK-SUBSTITUTION).
C     STORE THE CORRECTION IN DELTA.
      CALL DDASLV(NEQ,DELTA,WM,IWM)
C
C     UPDATE Y,E,AND YPRIME
      DO 360 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
         E(I)=E(I)-DELTA(I)
360      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
C
C     TEST FOR CONVERGENCE OF THE ITERATION
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM .LE. 100.D0*UROUND*PNORM) GO TO 375
      IF (M .GT. 0) GO TO 365
         OLDNRM = DELNRM
         GO TO 367
365   RATE = (DELNRM/OLDNRM)**(1.0D0/M)
      IF (RATE .GT. 0.90D0) GO TO 370
      S = RATE/(1.0D0 - RATE)
367   IF (S*DELNRM .LE. 0.33D0) GO TO 375
C
C     THE CORRECTOR HAS NOT YET CONVERGED.
C     UPDATE M AND TEST WHETHER THE
C     MAXIMUM NUMBER OF ITERATIONS HAVE
C     BEEN TRIED.
      M=M+1
      IF(M.GE.MAXIT)GO TO 370
C
C     EVALUATE THE RESIDUAL
C     AND GO BACK TO DO ANOTHER ITERATION
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL RES(X,Y,YPRIME,DELTA,IRES,
     *  RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
      GO TO 350
C
C
C     THE CORRECTOR FAILED TO CONVERGE IN MAXIT
C     ITERATIONS. IF THE ITERATION MATRIX
C     IS NOT CURRENT,RE-DO THE STEP WITH
C     A NEW ITERATION MATRIX.
370   CONTINUE
      IF(JCALC.EQ.0)GO TO 380
      JCALC=-1
      GO TO 300
C
C
C     THE ITERATION HAS CONVERGED.  IF NONNEGATIVITY OF SOLUTION IS
C     REQUIRED, SET THE SOLUTION NONNEGATIVE, IF THE PERTURBATION
C     TO DO IT IS SMALL ENOUGH.  IF THE CHANGE IS TOO LARGE, THEN
C     CONSIDER THE CORRECTOR ITERATION TO HAVE FAILED.
375   IF(NONNEG .EQ. 0) GO TO 390
      DO 377 I = 1,NEQ
377      DELTA(I) = MIN(Y(I),0.0D0)
      DELNRM = DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF(DELNRM .GT. 0.33D0) GO TO 380
      DO 378 I = 1,NEQ
378      E(I) = E(I) - DELTA(I)
      GO TO 390
C
C
C     EXITS FROM BLOCK 3
C     NO CONVERGENCE WITH CURRENT ITERATION
C     MATRIX,OR SINGULAR ITERATION MATRIX
380   CONVGD= .FALSE.
390   JCALC = 1
      IF(.NOT.CONVGD)GO TO 600
C
C
C
C
C
C------------------------------------------------------------------
C     BLOCK 4
C     ESTIMATE THE ERRORS AT ORDERS K,K-1,K-2
C     AS IF CONSTANT STEPSIZE WAS USED. ESTIMATE
C     THE LOCAL ERROR AT ORDER K AND TEST
C     WHETHER THE CURRENT STEP IS SUCCESSFUL.
C------------------------------------------------------------------
C
C     ESTIMATE ERRORS AT ORDERS K,K-1,K-2
      ENORM = DDANRM(NEQ,E,WT,RPAR,IPAR)
      ERK = SIGMA(K+1)*ENORM
      TERK = (K+1)*ERK
      EST = ERK
      KNEW=K
      IF(K .EQ. 1)GO TO 430
      DO 405 I = 1,NEQ
405     DELTA(I) = PHI(I,KP1) + E(I)
      ERKM1=SIGMA(K)*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM1 = K*ERKM1
      IF(K .GT. 2)GO TO 410
      IF(TERKM1 .LE. 0.5D0*TERK)GO TO 420
      GO TO 430
410   CONTINUE
      DO 415 I = 1,NEQ
415     DELTA(I) = PHI(I,K) + DELTA(I)
      ERKM2=SIGMA(K-1)*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM2 = (K-1)*ERKM2
      IF(MAX(TERKM1,TERKM2).GT.TERK)GO TO 430
C     LOWER THE ORDER
420   CONTINUE
      KNEW=K-1
      EST = ERKM1
C
C
C     CALCULATE THE LOCAL ERROR FOR THE CURRENT STEP
C     TO SEE IF THE STEP WAS SUCCESSFUL
430   CONTINUE
      ERR = CK * ENORM
      IF(ERR .GT. 1.0D0)GO TO 600
C
C
C
C
C
C------------------------------------------------------------------
C     BLOCK 5
C     THE STEP IS SUCCESSFUL. DETERMINE
C     THE BEST ORDER AND STEPSIZE FOR
C     THE NEXT STEP. UPDATE THE DIFFERENCES
C     FOR THE NEXT STEP.
C------------------------------------------------------------------
      IDID=1
      IWM(LNST)=IWM(LNST)+1
      KDIFF=K-KOLD
      KOLD=K
      HOLD=H
C
C
C     ESTIMATE THE ERROR AT ORDER K+1 UNLESS:
C        ALREADY DECIDED TO LOWER ORDER, OR
C        ALREADY USING MAXIMUM ORDER, OR
C        STEPSIZE NOT CONSTANT, OR
C        ORDER RAISED IN PREVIOUS STEP
      IF(KNEW.EQ.KM1.OR.K.EQ.IWM(LMXORD))IPHASE=1
      IF(IPHASE .EQ. 0)GO TO 545
      IF(KNEW.EQ.KM1)GO TO 540
      IF(K.EQ.IWM(LMXORD)) GO TO 550
      IF(KP1.GE.NS.OR.KDIFF.EQ.1)GO TO 550
      DO 510 I=1,NEQ
510      DELTA(I)=E(I)-PHI(I,KP2)
      ERKP1 = (1.0D0/(K+2))*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKP1 = (K+2)*ERKP1
      IF(K.GT.1)GO TO 520
      IF(TERKP1.GE.0.5D0*TERK)GO TO 550
      GO TO 530
520   IF(TERKM1.LE.MIN(TERK,TERKP1))GO TO 540
      IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD))GO TO 550
C
C     RAISE ORDER
530   K=KP1
      EST = ERKP1
      GO TO 550
C
C     LOWER ORDER
540   K=KM1
      EST = ERKM1
      GO TO 550
C
C     IF IPHASE = 0, INCREASE ORDER BY ONE AND MULTIPLY STEPSIZE BY
C     FACTOR TWO
545   K = KP1
      HNEW = H*2.0D0
      H = HNEW
      GO TO 575
C
C
C     DETERMINE THE APPROPRIATE STEPSIZE FOR
C     THE NEXT STEP.
550   HNEW=H
      TEMP2=K+1
      R=(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      IF(R .LT. 2.0D0) GO TO 555
      HNEW = 2.0D0*H
      GO TO 560
555   IF(R .GT. 1.0D0) GO TO 560
      R = MAX(0.5D0,MIN(0.9D0,R))
      HNEW = H*R
560   H=HNEW
C
C
C     UPDATE DIFFERENCES FOR NEXT STEP
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
C
C
C
C
C
C------------------------------------------------------------------
C     BLOCK 6
C     THE STEP IS UNSUCCESSFUL. RESTORE X,PSI,PHI
C     DETERMINE APPROPRIATE STEPSIZE FOR
C     CONTINUING THE INTEGRATION, OR EXIT WITH
C     AN ERROR FLAG IF THERE HAVE BEEN MANY
C     FAILURES.
C------------------------------------------------------------------
600   IPHASE = 1
C
C     RESTORE X,PHI,PSI
      X=XOLD
      IF(KP1.LT.NSP1)GO TO 630
      DO 620 J=NSP1,KP1
         TEMP1=1.0D0/BETA(J)
         DO 610 I=1,NEQ
610         PHI(I,J)=TEMP1*PHI(I,J)
620      CONTINUE
630   CONTINUE
      DO 640 I=2,KP1
640      PSI(I-1)=PSI(I)-H
C
C
C     TEST WHETHER FAILURE IS DUE TO CORRECTOR ITERATION
C     OR ERROR TEST
      IF(CONVGD)GO TO 660
      IWM(LCTF)=IWM(LCTF)+1
C
C
C     THE NEWTON ITERATION FAILED TO CONVERGE WITH
C     A CURRENT ITERATION MATRIX.  DETERMINE THE CAUSE
C     OF THE FAILURE AND TAKE APPROPRIATE ACTION.
      IF(IER.EQ.0)GO TO 650
C
C     THE ITERATION MATRIX IS SINGULAR. REDUCE
C     THE STEPSIZE BY A FACTOR OF 4. IF
C     THIS HAPPENS THREE TIMES IN A ROW ON
C     THE SAME STEP, RETURN WITH AN ERROR FLAG
      NSF=NSF+1
      R = 0.25D0
      H=H*R
      IF (NSF .LT. 3 .AND. ABS(H) .GE. HMIN) GO TO 690
      IDID=-8
      GO TO 675
C
C
C     THE NEWTON ITERATION FAILED TO CONVERGE FOR A REASON
C     OTHER THAN A SINGULAR ITERATION MATRIX.  IF IRES = -2, THEN
C     RETURN.  OTHERWISE, REDUCE THE STEPSIZE AND TRY AGAIN, UNLESS
C     TOO MANY FAILURES HAVE OCCURED.
650   CONTINUE
      IF (IRES .GT. -2) GO TO 655
      IDID = -11
      GO TO 675
655   NCF = NCF + 1
      R = 0.25D0
      H = H*R
      IF (NCF .LT. 10 .AND. ABS(H) .GE. HMIN) GO TO 690
      IDID = -7
      IF (IRES .LT. 0) IDID = -10
      IF (NEF .GE. 3) IDID = -9
      GO TO 675
C
C
C     THE NEWTON SCHEME CONVERGED,AND THE CAUSE
C     OF THE FAILURE WAS THE ERROR ESTIMATE
C     EXCEEDING THE TOLERANCE.
660   NEF=NEF+1
      IWM(LETF)=IWM(LETF)+1
      IF (NEF .GT. 1) GO TO 665
C
C     ON FIRST ERROR TEST FAILURE, KEEP CURRENT ORDER OR LOWER
C     ORDER BY ONE.  COMPUTE NEW STEPSIZE BASED ON DIFFERENCES
C     OF THE SOLUTION.
      K = KNEW
      TEMP2 = K + 1
      R = 0.90D0*(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      R = MAX(0.25D0,MIN(0.9D0,R))
      H = H*R
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     ON SECOND ERROR TEST FAILURE, USE THE CURRENT ORDER OR
C     DECREASE ORDER BY ONE.  REDUCE THE STEPSIZE BY A FACTOR OF
C     FOUR.
665   IF (NEF .GT. 2) GO TO 670
      K = KNEW
      H = 0.25D0*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     ON THIRD AND SUBSEQUENT ERROR TEST FAILURES, SET THE ORDER TO
C     ONE AND REDUCE THE STEPSIZE BY A FACTOR OF FOUR.
670   K = 1
      H = 0.25D0*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C
C
C
C     FOR ALL CRASHES, RESTORE Y TO ITS LAST VALUE,
C     INTERPOLATE TO FIND YPRIME AT LAST X, AND RETURN
675   CONTINUE
      CALL DDATRP(X,X,Y,YPRIME,NEQ,K,PHI,PSI)
      RETURN
C
C
C     GO BACK AND TRY THIS STEP AGAIN
690   GO TO 200
C
C------END OF SUBROUTINE DDASTP------
      END
      SUBROUTINE DDATRP (X, XOUT, YOUT, YPOUT, NEQ, KOLD, PHI, PSI)
C***BEGIN PROLOGUE  DDATRP
C***SUBSIDIARY
C***PURPOSE  Interpolation routine for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDATRP-S, DDATRP-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***DESCRIPTION
C------------------------------------------------------------------
C     THE METHODS IN SUBROUTINE DDASTP USE POLYNOMIALS
C     TO APPROXIMATE THE SOLUTION. DDATRP APPROXIMATES THE
C     SOLUTION AND ITS DERIVATIVE AT TIME XOUT BY EVALUATING
C     ONE OF THESE POLYNOMIALS,AND ITS DERIVATIVE,THERE.
C     INFORMATION DEFINING THIS POLYNOMIAL IS PASSED FROM
C     DDASTP, SO DDATRP CANNOT BE USED ALONE.
C
C     THE PARAMETERS ARE:
C     X     THE CURRENT TIME IN THE INTEGRATION.
C     XOUT  THE TIME AT WHICH THE SOLUTION IS DESIRED
C     YOUT  THE INTERPOLATED APPROXIMATION TO Y AT XOUT
C           (THIS IS OUTPUT)
C     YPOUT THE INTERPOLATED APPROXIMATION TO YPRIME AT XOUT
C           (THIS IS OUTPUT)
C     NEQ   NUMBER OF EQUATIONS
C     KOLD  ORDER USED ON LAST SUCCESSFUL STEP
C     PHI   ARRAY OF SCALED DIVIDED DIFFERENCES OF Y
C     PSI   ARRAY OF PAST STEPSIZE HISTORY
C------------------------------------------------------------------
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  DDATRP
C
      INTEGER  NEQ, KOLD
      DOUBLE PRECISION  X, XOUT, YOUT(*), YPOUT(*), PHI(NEQ,*), PSI(*)
C
      INTEGER  I, J, KOLDP1
      DOUBLE PRECISION  C, D, GAMMA, TEMP1
C
C***FIRST EXECUTABLE STATEMENT  DDATRP
      KOLDP1=KOLD+1
      TEMP1=XOUT-X
      DO 10 I=1,NEQ
         YOUT(I)=PHI(I,1)
10       YPOUT(I)=0.0D0
      C=1.0D0
      D=0.0D0
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
C
C------END OF SUBROUTINE DDATRP------
      END
      SUBROUTINE DDAWTS (NEQ, IWT, RTOL, ATOL, Y, WT, RPAR, IPAR)
C***BEGIN PROLOGUE  DDAWTS
C***SUBSIDIARY
C***PURPOSE  Set error weight vector for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDAWTS-S, DDAWTS-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***DESCRIPTION
C------------------------------------------------------------------
C     THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR
C     WT ACCORDING TO WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
C     I=1,-,N.
C     RTOL AND ATOL ARE SCALARS IF IWT = 0,
C     AND VECTORS IF IWT = 1.
C------------------------------------------------------------------
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  DDAWTS
C
      INTEGER  NEQ, IWT, IPAR(*)
      DOUBLE PRECISION  RTOL(*), ATOL(*), Y(*), WT(*), RPAR(*)
C
      INTEGER  I
      DOUBLE PRECISION  ATOLI, RTOLI
C
C***FIRST EXECUTABLE STATEMENT  DDAWTS
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 20 I=1,NEQ
         IF (IWT .EQ.0) GO TO 10
           RTOLI=RTOL(I)
           ATOLI=ATOL(I)
10         WT(I)=RTOLI*ABS(Y(I))+ATOLI
20         CONTINUE
      RETURN
C------END OF SUBROUTINE DDAWTS------------------------------------
      END
C*****END precision > double
C
C*****precision > single
C      SUBROUTINE SDAINI (X,Y,YPRIME, NEQ, RES, JAC, H, WT, IDID, RPAR,
C     +   IPAR, PHI, DELTA, E, WM, IWM, HMIN, UROUND, NONNEG, NTEMP)
CC***BEGIN PROLOGUE  SDAINI
CC***SUBSIDIARY
CC***PURPOSE  Initialization routine for SDASSL.
CC***LIBRARY   SLATEC (DASSL)
CC***TYPE      SINGLE PRECISION (SDAINI-S, DDAINI-D)
CC***AUTHOR  PETZOLD, LINDA R., (LLNL)
CC***DESCRIPTION
CC------------------------------------------------------------
CC     SDAINI TAKES ONE STEP OF SIZE H OR SMALLER
CC     WITH THE BACKWARD EULER METHOD, TO
CC     FIND YPRIME.  X AND Y ARE UPDATED TO BE CONSISTENT WITH THE
CC     NEW STEP.  A MODIFIED DAMPED NEWTON ITERATION IS USED TO
CC     SOLVE THE CORRECTOR ITERATION.
CC
CC     THE INITIAL GUESS FOR YPRIME IS USED IN THE
CC     PREDICTION, AND IN FORMING THE ITERATION
CC     MATRIX, BUT IS NOT INVOLVED IN THE
CC     ERROR TEST. THIS MAY HAVE TROUBLE
CC     CONVERGING IF THE INITIAL GUESS IS NO
CC     GOOD, OR IF G(X,Y,YPRIME) DEPENDS
CC     NONLINEARLY ON YPRIME.
CC
CC     THE PARAMETERS REPRESENT:
CC     X --         INDEPENDENT VARIABLE
CC     Y --         SOLUTION VECTOR AT X
CC     YPRIME --    DERIVATIVE OF SOLUTION VECTOR
CC     NEQ --       NUMBER OF EQUATIONS
CC     H --         STEPSIZE. IMDER MAY USE A STEPSIZE
CC                  SMALLER THAN H.
CC     WT --        VECTOR OF WEIGHTS FOR ERROR
CC                  CRITERION
CC     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS
CC                  IDID= 1 -- YPRIME WAS FOUND SUCCESSFULLY
CC                  IDID=-12 -- SDAINI FAILED TO FIND YPRIME
CC     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS
CC                  THAT ARE NOT ALTERED BY SDAINI
CC     PHI --       WORK SPACE FOR SDAINI
CC     DELTA,E --   WORK SPACE FOR SDAINI
CC     WM,IWM --    REAL AND INTEGER ARRAYS STORING
CC                  MATRIX INFORMATION
CC
CC------------------------------------------------------------
CC***ROUTINES CALLED  SDAJAC, SDANRM, SDASLV
CC***REVISION HISTORY  (YYMMDD)
CC   830315  DATE WRITTEN
CC   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
CC   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
CC   901026  Added explicit declarations for all variables and minor
CC           cosmetic changes to prologue.  (FNF)
CC   901030  Minor corrections to declarations.  (FNF)
CC***END PROLOGUE  SDAINI
CC
C      INTEGER  NEQ, IDID, IPAR(*), IWM(*), NONNEG, NTEMP
C      REAL  X,Y(*),YPRIME(*), H, WT(*), RPAR(*), PHI(NEQ,*), DELTA(*),
C     *   E(*), WM(*), HMIN, UROUND
C      EXTERNAL  RES, JAC
CC
C      EXTERNAL  SDAJAC, SDANRM, SDASLV
C      REAL  SDANRM
CC
C      INTEGER  I, IER, IRES, JCALC, LNJE, LNRE, M, MAXIT, MJAC, NCF,
C     *   NEF, NSF
C      REAL  CJ, DAMP, DELNRM, ERR, OLDNRM, R, RATE, S, XOLD, YNORM
C      LOGICAL  CONVGD
CC
C      PARAMETER (LNRE=12)
C      PARAMETER (LNJE=13)
CC
C      DATA MAXIT/10/,MJAC/5/
C      DATA DAMP/0.75E0/
CC
CC
CC---------------------------------------------------
CC     BLOCK 1.
CC     INITIALIZATIONS.
CC---------------------------------------------------
CC
CC***FIRST EXECUTABLE STATEMENT  SDAINI
C      IDID=1
C      NEF=0
C      NCF=0
C      NSF=0
C      XOLD=X
C      YNORM=SDANRM(NEQ,Y,WT,RPAR,IPAR)
CC
CC     SAVE Y AND YPRIME IN PHI
C      DO 100 I=1,NEQ
C         PHI(I,1)=Y(I)
C100      PHI(I,2)=YPRIME(I)
CC
CC
CC----------------------------------------------------
CC     BLOCK 2.
CC     DO ONE BACKWARD EULER STEP.
CC----------------------------------------------------
CC
CC     SET UP FOR START OF CORRECTOR ITERATION
C200   CJ=1.0E0/H
C      X=X+H
CC
CC     PREDICT SOLUTION AND DERIVATIVE
C      DO 250 I=1,NEQ
C250     Y(I)=Y(I)+H*YPRIME(I)
CC
C      JCALC=-1
C      M=0
C      CONVGD=.TRUE.
CC
CC
CC     CORRECTOR LOOP.
C300   IWM(LNRE)=IWM(LNRE)+1
C      IRES=0
CC
C      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
C      IF (IRES.LT.0) GO TO 430
CC
CC
CC     EVALUATE THE ITERATION MATRIX
C      IF (JCALC.NE.-1) GO TO 310
C      IWM(LNJE)=IWM(LNJE)+1
C      JCALC=0
C      CALL SDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,
C     *   IER,WT,E,WM,IWM,RES,IRES,
C     *   UROUND,JAC,RPAR,IPAR,NTEMP)
CC
C      S=1000000.E0
C      IF (IRES.LT.0) GO TO 430
C      IF (IER.NE.0) GO TO 430
C      NSF=0
CC
CC
CC
CC     MULTIPLY RESIDUAL BY DAMPING FACTOR
C310   CONTINUE
C      DO 320 I=1,NEQ
C320      DELTA(I)=DELTA(I)*DAMP
CC
CC     COMPUTE A NEW ITERATE (BACK SUBSTITUTION)
CC     STORE THE CORRECTION IN DELTA
CC
C      CALL SDASLV(NEQ,DELTA,WM,IWM)
CC
CC     UPDATE Y AND YPRIME
C      DO 330 I=1,NEQ
C         Y(I)=Y(I)-DELTA(I)
C330      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
CC
CC     TEST FOR CONVERGENCE OF THE ITERATION.
CC
C      DELNRM=SDANRM(NEQ,DELTA,WT,RPAR,IPAR)
C      IF (DELNRM.LE.100.E0*UROUND*YNORM)
C     *   GO TO 400
CC
C      IF (M.GT.0) GO TO 340
C         OLDNRM=DELNRM
C         GO TO 350
CC
C340   RATE=(DELNRM/OLDNRM)**(1.0E0/M)
C      IF (RATE.GT.0.90E0) GO TO 430
C      S=RATE/(1.0E0-RATE)
CC
C350   IF (S*DELNRM .LE. 0.33E0) GO TO 400
CC
CC
CC     THE CORRECTOR HAS NOT YET CONVERGED. UPDATE
CC     M AND AND TEST WHETHER THE MAXIMUM
CC     NUMBER OF ITERATIONS HAVE BEEN TRIED.
CC     EVERY MJAC ITERATIONS, GET A NEW
CC     ITERATION MATRIX.
CC
C      M=M+1
C      IF (M.GE.MAXIT) GO TO 430
CC
C      IF ((M/MJAC)*MJAC.EQ.M) JCALC=-1
C      GO TO 300
CC
CC
CC     THE ITERATION HAS CONVERGED.
CC     CHECK NONNEGATIVITY CONSTRAINTS
C400   IF (NONNEG.EQ.0) GO TO 450
C      DO 410 I=1,NEQ
C410      DELTA(I)=MIN(Y(I),0.0E0)
CC
C      DELNRM=SDANRM(NEQ,DELTA,WT,RPAR,IPAR)
C      IF (DELNRM.GT.0.33E0) GO TO 430
CC
C      DO 420 I=1,NEQ
C         Y(I)=Y(I)-DELTA(I)
C420      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
C      GO TO 450
CC
CC
CC     EXITS FROM CORRECTOR LOOP.
C430   CONVGD=.FALSE.
C450   IF (.NOT.CONVGD) GO TO 600
CC
CC
CC
CC-----------------------------------------------------
CC     BLOCK 3.
CC     THE CORRECTOR ITERATION CONVERGED.
CC     DO ERROR TEST.
CC-----------------------------------------------------
CC
C      DO 510 I=1,NEQ
C510      E(I)=Y(I)-PHI(I,1)
C      ERR=SDANRM(NEQ,E,WT,RPAR,IPAR)
CC
C      IF (ERR.LE.1.0E0) RETURN
CC
CC
CC
CC--------------------------------------------------------
CC     BLOCK 4.
CC     THE BACKWARD EULER STEP FAILED. RESTORE X, Y
CC     AND YPRIME TO THEIR ORIGINAL VALUES.
CC     REDUCE STEPSIZE AND TRY AGAIN, IF
CC     POSSIBLE.
CC---------------------------------------------------------
CC
C600   CONTINUE
C      X = XOLD
C      DO 610 I=1,NEQ
C         Y(I)=PHI(I,1)
C610      YPRIME(I)=PHI(I,2)
CC
C      IF (CONVGD) GO TO 640
C      IF (IER.EQ.0) GO TO 620
C         NSF=NSF+1
C         H=H*0.25E0
C         IF (NSF.LT.3.AND.ABS(H).GE.HMIN) GO TO 690
C         IDID=-12
C         RETURN
C620   IF (IRES.GT.-2) GO TO 630
C         IDID=-12
C         RETURN
C630   NCF=NCF+1
C      H=H*0.25E0
C      IF (NCF.LT.10.AND.ABS(H).GE.HMIN) GO TO 690
C         IDID=-12
C         RETURN
CC
C640   NEF=NEF+1
C      R=0.90E0/(2.0E0*ERR+0.0001E0)
C      R=MAX(0.1E0,MIN(0.5E0,R))
C      H=H*R
C      IF (ABS(H).GE.HMIN.AND.NEF.LT.10) GO TO 690
C         IDID=-12
C         RETURN
C690      GO TO 200
CC
CC-------------END OF SUBROUTINE SDAINI----------------------
C      END
C      SUBROUTINE SDAJAC (NEQ, X, Y, YPRIME, DELTA, CJ, H,
C     +   IER, WT, E, WM, IWM, RES, IRES, UROUND, JAC, RPAR,
C     +   IPAR, NTEMP)
CC***BEGIN PROLOGUE  SDAJAC
CC***SUBSIDIARY
CC***PURPOSE  Compute the iteration matrix for SDASSL and form the
CC            LU-decomposition.
CC***LIBRARY   SLATEC (DASSL)
CC***TYPE      SINGLE PRECISION (SDAJAC-S, DDAJAC-D)
CC***AUTHOR  PETZOLD, LINDA R., (LLNL)
CC***DESCRIPTION
CC------------------------------------------------------------------
CC     THIS ROUTINE COMPUTES THE ITERATION MATRIX
CC     PD=DG/DY+CJ*DG/DYPRIME (WHERE G(X,Y,YPRIME)=0).
CC     HERE PD IS COMPUTED BY THE USER-SUPPLIED
CC     ROUTINE JAC IF IWM(MTYPE) IS 1 OR 4, AND
CC     IT IS COMPUTED BY NUMERICAL FINITE DIFFERENCING
CC     IF IWM(MTYPE)IS 2 OR 5
CC     THE PARAMETERS HAVE THE FOLLOWING MEANINGS.
CC     Y        = ARRAY CONTAINING PREDICTED VALUES
CC     YPRIME   = ARRAY CONTAINING PREDICTED DERIVATIVES
CC     DELTA    = RESIDUAL EVALUATED AT (X,Y,YPRIME)
CC                (USED ONLY IF IWM(MTYPE)=2 OR 5)
CC     CJ       = SCALAR PARAMETER DEFINING ITERATION MATRIX
CC     H        = CURRENT STEPSIZE IN INTEGRATION
CC     IER      = VARIABLE WHICH IS .NE. 0
CC                IF ITERATION MATRIX IS SINGULAR,
CC                AND 0 OTHERWISE.
CC     WT       = VECTOR OF WEIGHTS FOR COMPUTING NORMS
CC     E        = WORK SPACE (TEMPORARY) OF LENGTH NEQ
CC     WM       = REAL WORK SPACE FOR MATRICES. ON
CC                OUTPUT IT CONTAINS THE LU DECOMPOSITION
CC                OF THE ITERATION MATRIX.
CC     IWM      = INTEGER WORK SPACE CONTAINING
CC                MATRIX INFORMATION
CC     RES      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
CC                TO EVALUATE THE RESIDUAL FUNCTION G(X,Y,YPRIME)
CC     IRES     = FLAG WHICH IS EQUAL TO ZERO IF NO ILLEGAL VALUES
CC                IN RES, AND LESS THAN ZERO OTHERWISE.  (IF IRES
CC                IS LESS THAN ZERO, THE MATRIX WAS NOT COMPLETED)
CC                IN THIS CASE (IF IRES .LT. 0), THEN IER = 0.
CC     UROUND   = THE UNIT ROUNDOFF ERROR OF THE MACHINE BEING USED.
CC     JAC      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
CC                TO EVALUATE THE ITERATION MATRIX (THIS ROUTINE
CC                IS ONLY USED IF IWM(MTYPE) IS 1 OR 4)
CC------------------------------------------------------------------
CC***ROUTINES CALLED  SGBFA, SGEFA
CC***REVISION HISTORY  (YYMMDD)
CC   830315  DATE WRITTEN
CC   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
CC   901010  Modified three MAX calls to be all on one line.  (FNF)
CC   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
CC   901026  Added explicit declarations for all variables and minor
CC           cosmetic changes to prologue.  (FNF)
CC   901101  Corrected PURPOSE.  (FNF)
CC***END PROLOGUE  SDAJAC
CC
C      INTEGER  NEQ, IER, IWM(*), IRES, IPAR(*), NTEMP
C      REAL  X, Y(*), YPRIME(*), DELTA(*), CJ, H, WT(*), E(*), WM(*),
C     *   UROUND, RPAR(*), SIGN77
C      EXTERNAL  RES, JAC, SIGN77
CC
C      EXTERNAL  SGBFA, SGEFA
CC
C      INTEGER  I, I1, I2, II, IPSAVE, ISAVE, J, K, L, LENPD, LIPVT,
C     *   LML, LMTYPE, LMU, MBA, MBAND, MEB1, MEBAND, MSAVE, MTYPE, N,
C     *   NPD, NPDM1, NROW
C      REAL  DEL, DELINV, SQUR, YPSAVE, YSAVE
CC
C      PARAMETER (NPD=1)
C      PARAMETER (LML=1)
C      PARAMETER (LMU=2)
C      PARAMETER (LMTYPE=4)
C      PARAMETER (LIPVT=21)
CC
CC***FIRST EXECUTABLE STATEMENT  SDAJAC
C      IER = 0
C      NPDM1=NPD-1
C      MTYPE=IWM(LMTYPE)
C      GO TO (100,200,300,400,500),MTYPE
CC
CC
CC     DENSE USER-SUPPLIED MATRIX
C100   LENPD=NEQ*NEQ
C      DO 110 I=1,LENPD
C110      WM(NPDM1+I)=0.0E0
C      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
C      GO TO 230
CC
CC
CC     DENSE FINITE-DIFFERENCE-GENERATED MATRIX
C200   IRES=0
C      NROW=NPDM1
C      SQUR = SQRT(UROUND)
C      DO 210 I=1,NEQ
C         DEL=SQUR*MAX(ABS(Y(I)),ABS(H*YPRIME(I)),ABS(WT(I)))
C         DEL=SIGN77(DEL,H*YPRIME(I))
C         DEL=(Y(I)+DEL)-Y(I)
C         YSAVE=Y(I)
C         YPSAVE=YPRIME(I)
C         Y(I)=Y(I)+DEL
C         YPRIME(I)=YPRIME(I)+CJ*DEL
C         CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
C         IF (IRES .LT. 0) RETURN
C         DELINV=1.0E0/DEL
C         DO 220 L=1,NEQ
C220      WM(NROW+L)=(E(L)-DELTA(L))*DELINV
C      NROW=NROW+NEQ
C      Y(I)=YSAVE
C      YPRIME(I)=YPSAVE
C210   CONTINUE
CC
CC
CC     DO DENSE-MATRIX LU DECOMPOSITION ON PD
C230      CALL SGEFA(WM(NPD),NEQ,NEQ,IWM(LIPVT),IER)
C      RETURN
CC
CC
CC     DUMMY SECTION FOR IWM(MTYPE)=3
C300   RETURN
CC
CC
CC     BANDED USER-SUPPLIED MATRIX
C400   LENPD=(2*IWM(LML)+IWM(LMU)+1)*NEQ
C      DO 410 I=1,LENPD
C410      WM(NPDM1+I)=0.0E0
C      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
C      MEBAND=2*IWM(LML)+IWM(LMU)+1
C      GO TO 550
CC
CC
CC     BANDED FINITE-DIFFERENCE-GENERATED MATRIX
C500   MBAND=IWM(LML)+IWM(LMU)+1
C      MBA=MIN(MBAND,NEQ)
C      MEBAND=MBAND+IWM(LML)
C      MEB1=MEBAND-1
C      MSAVE=(NEQ/MBAND)+1
C      ISAVE=NTEMP-1
C      IPSAVE=ISAVE+MSAVE
C      IRES=0
C      SQUR=SQRT(UROUND)
C      DO 540 J=1,MBA
C         DO 510 N=J,NEQ,MBAND
C          K= (N-J)/MBAND + 1
C          WM(ISAVE+K)=Y(N)
C          WM(IPSAVE+K)=YPRIME(N)
C          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
C          DEL=SIGN77(DEL,H*YPRIME(N))
C          DEL=(Y(N)+DEL)-Y(N)
C          Y(N)=Y(N)+DEL
C510       YPRIME(N)=YPRIME(N)+CJ*DEL
C      CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
C      IF (IRES .LT. 0) RETURN
C      DO 530 N=J,NEQ,MBAND
C          K= (N-J)/MBAND + 1
C          Y(N)=WM(ISAVE+K)
C          YPRIME(N)=WM(IPSAVE+K)
C          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
C          DEL=SIGN77(DEL,H*YPRIME(N))
C          DEL=(Y(N)+DEL)-Y(N)
C          DELINV=1.0E0/DEL
C          I1=MAX(1,(N-IWM(LMU)))
C          I2=MIN(NEQ,(N+IWM(LML)))
C          II=N*MEB1-IWM(LML)+NPDM1
C          DO 520 I=I1,I2
C520         WM(II+I)=(E(I)-DELTA(I))*DELINV
C530      CONTINUE
C540   CONTINUE
CC
CC
CC     DO LU DECOMPOSITION OF BANDED PD
C550   CALL SGBFA(WM(NPD),MEBAND,NEQ,
C     *    IWM(LML),IWM(LMU),IWM(LIPVT),IER)
C      RETURN
CC------END OF SUBROUTINE SDAJAC------
C      END
C      REAL FUNCTION SDANRM (NEQ, V, WT, RPAR, IPAR)
CC***BEGIN PROLOGUE  SDANRM
CC***SUBSIDIARY
CC***PURPOSE  Compute vector norm for SDASSL.
CC***LIBRARY   SLATEC (DASSL)
CC***TYPE      SINGLE PRECISION (SDANRM-S, DDANRM-D)
CC***AUTHOR  PETZOLD, LINDA R., (LLNL)
CC***DESCRIPTION
CC------------------------------------------------------------------
CC     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED
CC     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH
CC     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS
CC     CONTAINED IN THE ARRAY WT OF LENGTH NEQ.
CC        SDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
CC------------------------------------------------------------------
CC***ROUTINES CALLED  (NONE)
CC***REVISION HISTORY  (YYMMDD)
CC   830315  DATE WRITTEN
CC   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
CC   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
CC   901026  Added explicit declarations for all variables and minor
CC           cosmetic changes to prologue.  (FNF)
CC***END PROLOGUE  SDANRM
CC
C      INTEGER  NEQ, IPAR(*)
C      REAL  V(NEQ), WT(NEQ), RPAR(*)
CC
C      INTEGER  I
C      REAL  SUM, VMAX
CC
CC***FIRST EXECUTABLE STATEMENT  SDANRM
C      SDANRM = 0.0E0
C      VMAX = 0.0E0
C      DO 10 I = 1,NEQ
C        IF(ABS(V(I)/WT(I)) .GT. VMAX) VMAX = ABS(V(I)/WT(I))
C10      CONTINUE
C      IF(VMAX .LE. 0.0E0) GO TO 30
C      SUM = 0.0E0
C      DO 20 I = 1,NEQ
C20      SUM = SUM + ((V(I)/WT(I))/VMAX)**2
C      SDANRM = VMAX*SQRT(SUM/NEQ)
C30    CONTINUE
C      RETURN
CC------END OF FUNCTION SDANRM------
C      END
C      SUBROUTINE SDASLV (NEQ, DELTA, WM, IWM)
CC***BEGIN PROLOGUE  SDASLV
CC***SUBSIDIARY
CC***PURPOSE  Linear system solver for SDASSL.
CC***LIBRARY   SLATEC (DASSL)
CC***TYPE      SINGLE PRECISION (SDASLV-S, DDASLV-D)
CC***AUTHOR  PETZOLD, LINDA R., (LLNL)
CC***DESCRIPTION
CC------------------------------------------------------------------
CC     THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR
CC     SYSTEM ARISING IN THE NEWTON ITERATION.
CC     MATRICES AND REAL TEMPORARY STORAGE AND
CC     REAL INFORMATION ARE STORED IN THE ARRAY WM.
CC     INTEGER MATRIX INFORMATION IS STORED IN
CC     THE ARRAY IWM.
CC     FOR A DENSE MATRIX, THE LINPACK ROUTINE
CC     SGESL IS CALLED.
CC     FOR A BANDED MATRIX,THE LINPACK ROUTINE
CC     SGBSL IS CALLED.
CC------------------------------------------------------------------
CC***ROUTINES CALLED  SGBSL, SGESL
CC***REVISION HISTORY  (YYMMDD)
CC   830315  DATE WRITTEN
CC   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
CC   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
CC   901026  Added explicit declarations for all variables and minor
CC           cosmetic changes to prologue.  (FNF)
CC***END PROLOGUE  SDASLV
CC
C      INTEGER  NEQ, IWM(*)
C      REAL  DELTA(*), WM(*)
CC
C      EXTERNAL  SGBSL, SGESL
CC
C      INTEGER  LIPVT, LML, LMU, LMTYPE, MEBAND, MTYPE, NPD
C      PARAMETER (NPD=1)
C      PARAMETER (LML=1)
C      PARAMETER (LMU=2)
C      PARAMETER (LMTYPE=4)
C      PARAMETER (LIPVT=21)
CC
CC***FIRST EXECUTABLE STATEMENT  SDASLV
C      MTYPE=IWM(LMTYPE)
C      GO TO(100,100,300,400,400),MTYPE
CC
CC     DENSE MATRIX
C100   CALL SGESL(WM(NPD),NEQ,NEQ,IWM(LIPVT),DELTA,0)
C      RETURN
CC
CC     DUMMY SECTION FOR MTYPE=3
C300   CONTINUE
C      RETURN
CC
CC     BANDED MATRIX
C400   MEBAND=2*IWM(LML)+IWM(LMU)+1
C      CALL SGBSL(WM(NPD),MEBAND,NEQ,IWM(LML),
C     *  IWM(LMU),IWM(LIPVT),DELTA,0)
C      RETURN
CC------END OF SUBROUTINE SDASLV------
C      END
C      SUBROUTINE SDASSL (RES,NEQ,T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
C     +   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
CC***BEGIN PROLOGUE  SDASSL
CC***PURPOSE  This code solves a system of differential/algebraic
CC            equations of the form G(T,Y,YPRIME) = 0.
CC***LIBRARY   SLATEC (DASSL)
CC***CATEGORY  I1A2
CC***TYPE      SINGLE PRECISION (SDASSL-S, DDASSL-D)
CC***KEYWORDS  DIFFERENTIAL/ALGEBRAIC, BACKWARD DIFFERENTIATION
CC             FORMULAS, IMPLICIT DIFFERENTIAL SYSTEMS
CC***AUTHOR  PETZOLD, LINDA R., (LLNL)
CC             COMPUTING AND MATHEMATICS RESEARCH DIVISION
CC             LAWRENCE LIVERMORE NATIONAL LABORATORY
CC             L - 316, P.O. BOX 808,
CC             LIVERMORE, CA.    94550
CC***DESCRIPTION
CC
CC *Usage:
CC
CC      EXTERNAL RES, JAC
CC      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR
CC      REAL T, Y(NEQ), YPRIME(NEQ), TOUT, RTOL, ATOL,
CC     *   RWORK(LRW), RPAR
CC
CC      CALL SDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
CC     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
CC
CC
CC *Arguments:
CC
CC  RES:EXT     This is a subroutine which you provide to define the
CC              differential/algebraic system.
CC
CC  NEQ:IN      This is the number of equations to be solved.
CC
CC  T:INOUT     This is the current value of the independent variable.
CC
CC  Y(*):INOUT  This array contains the solution components at T.
CC
CC  YPRIME(*):INOUT  This array contains the derivatives of the solution
CC              components at T.
CC
CC  TOUT:IN     This is a point at which a solution is desired.
CC
CC  INFO(N):IN  The basic task of the code is to solve the system from T
CC              to TOUT and return answer at TOUT.  INFO is an integer
CC              array which is used to communicate exactly how you want
CC              this task to be carried out.  (See below for details.)
CC              N must be greater than or equal to 15.
CC
CC  RTOL,ATOL:INOUT  These quantities represent relative and absolute
CC              error tolerances which you provide to indicate how
CC              accurately you wish the solution to be computed.  You
CC              may choose them to be both scalars or else both vectors.
CC              Caution:  In Fortran 77, a scalar is not the same as an
CC                        array of length 1.  Some compilers may object
CC                        to using scalars for RTOL,ATOL.
CC
CC  IDID:OUT    This scalar quantity is an indicator reporting what the
CC              code did.  You must monitor this integer variable to
CC              decide  what action to take next.
CC
CC  RWORK:WORK  A real work array of length LRW which provides the
CC              code with needed storage space.
CC
CC  LRW:IN      The length of RWORK.  (See below for required length.)
CC
CC  IWORK:WORK  An integer work array of length LIW which probides the
CC              code with needed storage space.
CC
CC  LIW:IN      The length of IWORK.  (See below for required length.)
CC
CC  RPAR,IPAR:IN  These are real and integer parameter arrays which
CC              you can use for communication between your calling
CC              program and the RES subroutine (and the JAC subroutine)
CC
CC  JAC:EXT     This is the name of a subroutine which you may choose
CC              to provide for defining a matrix of partial derivatives
CC              described below.
CC
CC  Quantities which may be altered by SDASSL are:
CC     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL,
CC     IDID, RWORK(*) AND IWORK(*)
CC
CC *Description
CC
CC  Subroutine SDASSL uses the backward differentiation formulas of
CC  orders one thru five to solve a system of the above form for Y and
CC  YPRIME.  Values for Y and YPRIME at initial time must be given as
CC  input.  These values must be consistent, (that is, if T,Y,YPRIME are
CC  the initial values, they must satisfy G(T,Y,YPRIME) = 0.).  The
CC  subroutine solves the system from T to TOUT.  It is easy to continue
CC  solution to get results at additional TOUT.  This is theinterval
CC  mode of operation.  Intermediate results can also be obtained easily
CC  by using the intermediate-output capability.
CC
CC  The following detailed description is divided into subsections:
CC    1. Input required for the first call to SDASSL.
CC    2. Output after any return from SDASSL.
CC    3. What to do to continue the integration.
CC    4. Error messages.
CC
CC
CC  ----- INPUT -- WHAT TO DO ON THE FIRST CALL TO SDASSL ------------
CC
CC  The first call of the code is defined to be the start of each new
CC  problem. Read through the descriptions of all the following items,
CC  provide sufficient storage space for designated arrays, set
CC  appropriate variables for the initialization of the problem, and
CC  give information about how you want the problem to be solved.
CC
CC
CC  RES -- Provide a subroutine of the form
CC             SUBROUTINE RES(T,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
CC         to define the system of differential/algebraic
CC         equations which is to be solved. For the given values
CC         of T,Y and YPRIME, the subroutine should
CC         return the residual of the defferential/algebraic
CC         system
CC             DELTA = G(T,Y,YPRIME)
CC         (DELTA(*) is a vector of length NEQ which is
CC         output for RES.)
CC
CC         Subroutine RES must not alter T,Y or YPRIME.
CC         You must declare the name RES in an external
CC         statement in your program that calls SDASSL.
CC         You must dimension Y,YPRIME and DELTA in RES.
CC
CC         IRES is an integer flag which is always equal to
CC         zero on input. Subroutine RES should alter IRES
CC         only if it encounters an illegal value of Y or
CC         a stop condition. Set IRES = -1 if an input value
CC         is illegal, and SDASSL will try to solve the problem
CC         without getting IRES = -1. If IRES = -2, SDASSL
CC         will return control to the calling program
CC         with IDID = -11.
CC
CC         RPAR and IPAR are real and integer parameter arrays which
CC         you can use for communication between your calling program
CC         and subroutine RES. They are not altered by SDASSL. If you
CC         do not need RPAR or IPAR, ignore these parameters by treat-
CC         ing them as dummy arguments. If you do choose to use them,
CC         dimension them in your calling program and in RES as arrays
CC         of appropriate length.
CC
CC  NEQ -- Set it to the number of differential equations.
CC         (NEQ .GE. 1)
CC
CC  T -- Set it to the initial point of the integration.
CC         T must be defined as a variable.
CC
CC  Y(*) -- Set this vector to the initial values of the NEQ solution
CC         components at the initial point. You must dimension Y of
CC         length at least NEQ in your calling program.
CC
CC  YPRIME(*) -- Set this vector to the initial values of the NEQ
CC         first derivatives of the solution components at the initial
CC         point.  You must dimension YPRIME at least NEQ in your
CC         calling program. If you do not know initial values of some
CC         of the solution components, see the explanation of INFO(11).
CC
CC  TOUT -- Set it to the first point at which a solution
CC         is desired. You can not take TOUT = T.
CC         integration either forward in T (TOUT .GT. T) or
CC         backward in T (TOUT .LT. T) is permitted.
CC
CC         The code advances the solution from T to TOUT using
CC         step sizes which are automatically selected so as to
CC         achieve the desired accuracy. If you wish, the code will
CC         return with the solution and its derivative at
CC         intermediate steps (intermediate-output mode) so that
CC         you can monitor them, but you still must provide TOUT in
CC         accord with the basic aim of the code.
CC
CC         The first step taken by the code is a critical one
CC         because it must reflect how fast the solution changes near
CC         the initial point. The code automatically selects an
CC         initial step size which is practically always suitable for
CC         the problem. By using the fact that the code will not step
CC         past TOUT in the first step, you could, if necessary,
CC         restrict the length of the initial step size.
CC
CC         For some problems it may not be permissible to integrate
CC         past a point TSTOP because a discontinuity occurs there
CC         or the solution or its derivative is not defined beyond
CC         TSTOP. When you have declared a TSTOP point (SEE INFO(4)
CC         and RWORK(1)), you have told the code not to integrate
CC         past TSTOP. In this case any TOUT beyond TSTOP is invalid
CC         input.
CC
CC  INFO(*) -- Use the INFO array to give the code more details about
CC         how you want your problem solved.  This array should be
CC         dimensioned of length 15, though SDASSL uses only the first
CC         eleven entries.  You must respond to all of the following
CC         items, which are arranged as questions.  The simplest use
CC         of the code corresponds to answering all questions as yes,
CC         i.e. setting all entries of INFO to 0.
CC
CC       INFO(1) - This parameter enables the code to initialize
CC              itself. You must set it to indicate the start of every
CC              new problem.
CC
CC          **** Is this the first call for this problem ...
CC                Yes - Set INFO(1) = 0
CC                 No - Not applicable here.
CC                      See below for continuation calls.  ****
CC
CC       INFO(2) - How much accuracy you want of your solution
CC              is specified by the error tolerances RTOL and ATOL.
CC              The simplest use is to take them both to be scalars.
CC              To obtain more flexibility, they can both be vectors.
CC              The code must be told your choice.
CC
CC          **** Are both error tolerances RTOL, ATOL scalars ...
CC                Yes - Set INFO(2) = 0
CC                      and input scalars for both RTOL and ATOL
CC                 No - Set INFO(2) = 1
CC                      and input arrays for both RTOL and ATOL ****
CC
CC       INFO(3) - The code integrates from T in the direction
CC              of TOUT by steps. If you wish, it will return the
CC              computed solution and derivative at the next
CC              intermediate step (the intermediate-output mode) or
CC              TOUT, whichever comes first. This is a good way to
CC              proceed if you want to see the behavior of the solution.
CC              If you must have solutions at a great many specific
CC              TOUT points, this code will compute them efficiently.
CC
CC          **** Do you want the solution only at
CC                TOUT (and not at the next intermediate step) ...
CC                 Yes - Set INFO(3) = 0
CC                  No - Set INFO(3) = 1 ****
CC
CC       INFO(4) - To handle solutions at a great many specific
CC              values TOUT efficiently, this code may integrate past
CC              TOUT and interpolate to obtain the result at TOUT.
CC              Sometimes it is not possible to integrate beyond some
CC              point TSTOP because the equation changes there or it is
CC              not defined past TSTOP. Then you must tell the code
CC              not to go past.
CC
CC           **** Can the integration be carried out without any
CC                restrictions on the independent variable T ...
CC                 Yes - Set INFO(4)=0
CC                  No - Set INFO(4)=1
CC                       and define the stopping point TSTOP by
CC                       setting RWORK(1)=TSTOP ****
CC
CC       INFO(5) - To solve differential/algebraic problems it is
CC              necessary to use a matrix of partial derivatives of the
CC              system of differential equations. If you do not
CC              provide a subroutine to evaluate it analytically (see
CC              description of the item JAC in the call list), it will
CC              be approximated by numerical differencing in this code.
CC              although it is less trouble for you to have the code
CC              compute partial derivatives by numerical differencing,
CC              the solution will be more reliable if you provide the
CC              derivatives via JAC. Sometimes numerical differencing
CC              is cheaper than evaluating derivatives in JAC and
CC              sometimes it is not - this depends on your problem.
CC
CC           **** Do you want the code to evaluate the partial
CC                derivatives automatically by numerical differences ...
CC                   Yes - Set INFO(5)=0
CC                    No - Set INFO(5)=1
CC                  and provide subroutine JAC for evaluating the
CC                  matrix of partial derivatives ****
CC
CC       INFO(6) - SDASSL will perform much better if the matrix of
CC              partial derivatives, DG/DY + CJ*DG/DYPRIME,
CC              (here CJ is a scalar determined by SDASSL)
CC              is banded and the code is told this. In this
CC              case, the storage needed will be greatly reduced,
CC              numerical differencing will be performed much cheaper,
CC              and a number of important algorithms will execute much
CC              faster. The differential equation is said to have
CC              half-bandwidths ML (lower) and MU (upper) if equation i
CC              involves only unknowns Y(J) with
CC                             I-ML .LE. J .LE. I+MU
CC              for all I=1,2,...,NEQ. Thus, ML and MU are the widths
CC              of the lower and upper parts of the band, respectively,
CC              with the main diagonal being excluded. If you do not
CC              indicate the equation has a banded matrix of partial
CC              derivatives, the code works with a full matrix of NEQ**2
CC              elements (stored in the conventional way). Computations
CC              with banded matrices cost less time & storage than with
CC              full matrices if 2*ML+MU .LT. NEQ. If you tell the
CC              code that the matrix of partial derivatives has a banded
CC              structure and you want to provide subroutine JAC to
CC              compute partial derivatives, then you must be careful
CC              to store the elements of the matrix in the special form
CC              indicated in the description of JAC.
CC
CC          **** Do you want to solve the problem using a full
CC               (dense) matrix (and not a special banded
CC               structure) ...
CC                Yes - Set INFO(6)=0
CC                 No - Set INFO(6)=1
CC                       and provide the lower (ML) and upper (MU)
CC                       bandwidths by setting
CC                       IWORK(1)=ML
CC                       IWORK(2)=MU ****
CC
CC
CC        INFO(7) -- You can specify a maximum (absolute value of)
CC              stepsize, so that the code
CC              will avoid passing over very
CC              large regions.
CC
CC          ****  Do you want the code to decide
CC                on its own maximum stepsize?
CC                Yes - Set INFO(7)=0
CC                 No - Set INFO(7)=1
CC                      and define HMAX by setting
CC                      RWORK(2)=HMAX ****
CC
CC        INFO(8) -- Differential/algebraic problems
CC              may occaisionally suffer from
CC              severe scaling difficulties on the
CC              first step. If you know a great deal
CC              about the scaling of your problem, you can
CC              help to alleviate this problem by
CC              specifying an initial stepsize HO.
CC
CC          ****  Do you want the code to define
CC                its own initial stepsize?
CC                Yes - Set INFO(8)=0
CC                 No - Set INFO(8)=1
CC                      and define HO by setting
CC                      RWORK(3)=HO ****
CC
CC        INFO(9) -- If storage is a severe problem,
CC              you can save some locations by
CC              restricting the maximum order MAXORD.
CC              the default value is 5. for each
CC              order decrease below 5, the code
CC              requires NEQ fewer locations, however
CC              it is likely to be slower. In any
CC              case, you must have 1 .LE. MAXORD .LE. 5
CC          ****  Do you want the maximum order to
CC                default to 5?
CC                Yes - Set INFO(9)=0
CC                 No - Set INFO(9)=1
CC                      and define MAXORD by setting
CC                      IWORK(3)=MAXORD ****
CC
CC        INFO(10) --If you know that the solutions to your equations
CC               will always be nonnegative, it may help to set this
CC               parameter. However, it is probably best to
CC               try the code without using this option first,
CC               and only to use this option if that doesn't
CC               work very well.
CC           ****  Do you want the code to solve the problem without
CC                 invoking any special nonnegativity constraints?
CC                  Yes - Set INFO(10)=0
CC                   No - Set INFO(10)=1
CC
CC        INFO(11) --SDASSL normally requires the initial T,
CC               Y, and YPRIME to be consistent. That is,
CC               you must have G(T,Y,YPRIME) = 0 at the initial
CC               time. If you do not know the initial
CC               derivative precisely, you can let SDASSL try
CC               to compute it.
CC          ****   Are the initialHE INITIAL T, Y, YPRIME consistent?
CC                 Yes - Set INFO(11) = 0
CC                  No - Set INFO(11) = 1,
CC                       and set YPRIME to an initial approximation
CC                       to YPRIME.  (If you have no idea what
CC                       YPRIME should be, set it to zero. Note
CC                       that the initial Y should be such
CC                       that there must exist a YPRIME so that
CC                       G(T,Y,YPRIME) = 0.)
CC
CC  RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL
CC         error tolerances to tell the code how accurately you
CC         want the solution to be computed.  They must be defined
CC         as variables because the code may change them.  You
CC         have two choices --
CC               Both RTOL and ATOL are scalars. (INFO(2)=0)
CC               Both RTOL and ATOL are vectors. (INFO(2)=1)
CC         in either case all components must be non-negative.
CC
CC         The tolerances are used by the code in a local error
CC         test at each step which requires roughly that
CC               ABS(LOCAL ERROR) .LE. RTOL*ABS(Y)+ATOL
CC         for each vector component.
CC         (More specifically, a root-mean-square norm is used to
CC         measure the size of vectors, and the error test uses the
CC         magnitude of the solution at the beginning of the step.)
CC
CC         The true (global) error is the difference between the
CC         true solution of the initial value problem and the
CC         computed approximation.  Practically all present day
CC         codes, including this one, control the local error at
CC         each step and do not even attempt to control the global
CC         error directly.
CC         Usually, but not always, the true accuracy of the
CC         computed Y is comparable to the error tolerances. This
CC         code will usually, but not always, deliver a more
CC         accurate solution if you reduce the tolerances and
CC         integrate again.  By comparing two such solutions you
CC         can get a fairly reliable idea of the true error in the
CC         solution at the bigger tolerances.
CC
CC         Setting ATOL=0. results in a pure relative error test on
CC         that component.  Setting RTOL=0. results in a pure
CC         absolute error test on that component.  A mixed test
CC         with non-zero RTOL and ATOL corresponds roughly to a
CC         relative error test when the solution component is much
CC         bigger than ATOL and to an absolute error test when the
CC         solution component is smaller than the threshhold ATOL.
CC
CC         The code will not attempt to compute a solution at an
CC         accuracy unreasonable for the machine being used.  It will
CC         advise you if you ask for too much accuracy and inform
CC         you as to the maximum accuracy it believes possible.
CC
CC  RWORK(*) --  Dimension this real work array of length LRW in your
CC         calling program.
CC
CC  LRW -- Set it to the declared length of the RWORK array.
CC               You must have
CC                    LRW .GE. 40+(MAXORD+4)*NEQ+NEQ**2
CC               for the full (dense) JACOBIAN case (when INFO(6)=0), or
CC                    LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
CC               for the banded user-defined JACOBIAN case
CC               (when INFO(5)=1 and INFO(6)=1), or
CC                     LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
CC                           +2*(NEQ/(ML+MU+1)+1)
CC               for banded finite-difference-generated JACOBIAN case
CC               (when INFO(5)=0 and INFO(6)=1)
CC
CC  IWORK(*) --  Dimension this integer work array of length LIW in
CC         your calling program.
CC
CC  LIW -- Set it to the declared length of the IWORK array.
CC               You must have LIW .GE. 20+NEQ
CC
CC  RPAR, IPAR -- These are parameter arrays, of real and integer
CC         type, respectively.  You can use them for communication
CC         between your program that calls SDASSL and the
CC         RES subroutine (and the JAC subroutine).  They are not
CC         altered by SDASSL.  If you do not need RPAR or IPAR,
CC         ignore these parameters by treating them as dummy
CC         arguments.  If you do choose to use them, dimension
CC         them in your calling program and in RES (and in JAC)
CC         as arrays of appropriate length.
CC
CC  JAC -- If you have set INFO(5)=0, you can ignore this parameter
CC         by treating it as a dummy argument.  Otherwise, you must
CC         provide a subroutine of the form
CC             SUBROUTINE JAC(T,Y,YPRIME,PD,CJ,RPAR,IPAR)
CC         to define the matrix of partial derivatives
CC             PD=DG/DY+CJ*DG/DYPRIME
CC         CJ is a scalar which is input to JAC.
CC         For the given values of T,Y,YPRIME, the
CC         subroutine must evaluate the non-zero partial
CC         derivatives for each equation and each solution
CC         component, and store these values in the
CC         matrix PD.  The elements of PD are set to zero
CC         before each call to JAC so only non-zero elements
CC         need to be defined.
CC
CC         Subroutine JAC must not alter T,Y,(*),YPRIME(*), or CJ.
CC         You must declare the name JAC in an EXTERNAL statement in
CC         your program that calls SDASSL.  You must dimension Y,
CC         YPRIME and PD in JAC.
CC
CC         The way you must store the elements into the PD matrix
CC         depends on the structure of the matrix which you
CC         indicated by INFO(6).
CC               *** INFO(6)=0 -- Full (dense) matrix ***
CC                   Give PD a first dimension of NEQ.
CC                   When you evaluate the (non-zero) partial derivative
CC                   of equation I with respect to variable J, you must
CC                   store it in PD according to
CC                   PD(I,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
CC               *** INFO(6)=1 -- Banded JACOBIAN with ML lower and MU
CC                   upper diagonal bands (refer to INFO(6) description
CC                   of ML and MU) ***
CC                   Give PD a first dimension of 2*ML+MU+1.
CC                   when you evaluate the (non-zero) partial derivative
CC                   of equation I with respect to variable J, you must
CC                   store it in PD according to
CC                   IROW = I - J + ML + MU + 1
CC                   PD(IROW,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
CC
CC         RPAR and IPAR are real and integer parameter arrays
CC         which you can use for communication between your calling
CC         program and your JACOBIAN subroutine JAC. They are not
CC         altered by SDASSL. If you do not need RPAR or IPAR,
CC         ignore these parameters by treating them as dummy
CC         arguments. If you do choose to use them, dimension
CC         them in your calling program and in JAC as arrays of
CC         appropriate length.
CC
CC
CC  OPTIONALLY REPLACEABLE NORM ROUTINE:
CC
CC     SDASSL uses a weighted norm SDANRM to measure the size
CC     of vectors such as the estimated error in each step.
CC     A FUNCTION subprogram
CC       REAL FUNCTION SDANRM(NEQ,V,WT,RPAR,IPAR)
CC       DIMENSION V(NEQ),WT(NEQ)
CC     is used to define this norm. Here, V is the vector
CC     whose norm is to be computed, and WT is a vector of
CC     weights.  A SDANRM routine has been included with SDASSL
CC     which computes the weighted root-mean-square norm
CC     given by
CC       SDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
CC     this norm is suitable for most problems. In some
CC     special cases, it may be more convenient and/or
CC     efficient to define your own norm by writing a function
CC     subprogram to be called instead of SDANRM. This should,
CC     however, be attempted only after careful thought and
CC     consideration.
CC
CC
CC  ----- OUTPUT -- AFTER ANY RETURN FROM SDASSL ---------------------
CC
CC  The principal aim of the code is to return a computed solution at
CC  TOUT, although it is also possible to obtain intermediate results
CC  along the way. To find out whether the code achieved its goal
CC  or if the integration process was interrupted before the task was
CC  completed, you must check the IDID parameter.
CC
CC
CC  T -- The solution was successfully advanced to the
CC               output value of T.
CC
CC  Y(*) -- Contains the computed solution approximation at T.
CC
CC  YPRIME(*) -- Contains the computed derivative
CC               approximation at T.
CC
CC  IDID -- Reports what the code did.
CC
CC                     *** Task completed ***
CC                Reported by positive values of IDID
CC
CC           IDID = 1 -- A step was successfully taken in the
CC                   intermediate-output mode. The code has not
CC                   yet reached TOUT.
CC
CC           IDID = 2 -- The integration to TSTOP was successfully
CC                   completed (T=TSTOP) by stepping exactly to TSTOP.
CC
CC           IDID = 3 -- The integration to TOUT was successfully
CC                   completed (T=TOUT) by stepping past TOUT.
CC                   Y(*) is obtained by interpolation.
CC                   YPRIME(*) is obtained by interpolation.
CC
CC                    *** Task interrupted ***
CC                Reported by negative values of IDID
CC
CC           IDID = -1 -- A large amount of work has been expended.
CC                   (About 500 steps)
CC
CC           IDID = -2 -- The error tolerances are too stringent.
CC
CC           IDID = -3 -- The local error test cannot be satisfied
CC                   because you specified a zero component in ATOL
CC                   and the corresponding computed solution
CC                   component is zero. Thus, a pure relative error
CC                   test is impossible for this component.
CC
CC           IDID = -6 -- SDASSL had repeated error test
CC                   failures on the last attempted step.
CC
CC           IDID = -7 -- The corrector could not converge.
CC
CC           IDID = -8 -- The matrix of partial derivatives
CC                   is singular.
CC
CC           IDID = -9 -- The corrector could not converge.
CC                   there were repeated error test failures
CC                   in this step.
CC
CC           IDID =-10 -- The corrector could not converge
CC                   because IRES was equal to minus one.
CC
CC           IDID =-11 -- IRES equal to -2 was encountered
CC                   and control is being returned to the
CC                   calling program.
CC
CC           IDID =-12 -- SDASSL failed to compute the initial
CC                   YPRIME.
CC
CC
CC
CC           IDID = -13,..,-32 -- Not applicable for this code
CC
CC                    *** Task terminated ***
CC                Reported by the value of IDID=-33
CC
CC           IDID = -33 -- The code has encountered trouble from which
CC                   it cannot recover. A message is printed
CC                   explaining the trouble and control is returned
CC                   to the calling program. For example, this occurs
CC                   when invalid input is detected.
CC
CC  RTOL, ATOL -- These quantities remain unchanged except when
CC               IDID = -2. In this case, the error tolerances have been
CC               increased by the code to values which are estimated to
CC               be appropriate for continuing the integration. However,
CC               the reported solution at T was obtained using the input
CC               values of RTOL and ATOL.
CC
CC  RWORK, IWORK -- Contain information which is usually of no
CC               interest to user but necessary for subsequent calls.
CC               However, you may find use for
CC
CC               RWORK(3)--Which contains the step size H to be
CC                       attempted on the next step.
CC
CC               RWORK(4)--Which contains the current value of the
CC                       independent variable, i.e., the farthest point
CC                       integration has reached. This will be different
CC                       from T only when interpolation has been
CC                       performed (IDID=3).
CC
CC               RWORK(7)--Which contains the stepsize used
CC                       on the last successful step.
CC
CC               IWORK(7)--Which contains the order of the method to
CC                       be attempted on the next step.
CC
CC               IWORK(8)--Which contains the order of the method used
CC                       on the last step.
CC
CC               IWORK(11)--Which contains the number of steps taken so
CC                        far.
CC
CC               IWORK(12)--Which contains the number of calls to RES
CC                        so far.
CC
CC               IWORK(13)--Which contains the number of evaluations of
CC                        the matrix of partial derivatives needed so
CC                        far.
CC
CC               IWORK(14)--Which contains the total number
CC                        of error test failures so far.
CC
CC               IWORK(15)--Which contains the total number
CC                        of convergence test failures so far.
CC                        (includes singular iteration matrix
CC                        failures.)
CC
CC
CC  ----- INPUT -- WHAT TO DO TO CONTINUE THE INTEGRATION ------------
CC                    (CALLS AFTER THE FIRST)
CC
CC  This code is organized so that subsequent calls to continue the
CC  integration involve little (if any) additional effort on your
CC  part. You must monitor the IDID parameter in order to determine
CC  what to do next.
CC
CC  Recalling that the principal task of the code is to integrate
CC  from T to TOUT (the interval mode), usually all you will need
CC  to do is specify a new TOUT upon reaching the current TOUT.
CC
CC  Do not alter any quantity not specifically permitted below,
CC  in particular do not alter NEQ,T,Y(*),YPRIME(*),RWORK(*),IWORK(*)
CC  or the differential equation in subroutine RES. Any such
CC  alteration constitutes a new problem and must be treated as such,
CC  i.e., you must start afresh.
CC
CC  You cannot change from vector to scalar error control or vice
CC  versa (INFO(2)), but you can change the size of the entries of
CC  RTOL, ATOL. Increasing a tolerance makes the equation easier
CC  to integrate. Decreasing a tolerance will make the equation
CC  harder to integrate and should generally be avoided.
CC
CC  You can switch from the intermediate-output mode to the
CC  interval mode (INFO(3)) or vice versa at any time.
CC
CC  If it has been necessary to prevent the integration from going
CC  past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
CC  code will not integrate to any TOUT beyond the currently
CC  specified TSTOP. Once TSTOP has been reached you must change
CC  the value of TSTOP or set INFO(4)=0. You may change INFO(4)
CC  or TSTOP at any time but you must supply the value of TSTOP in
CC  RWORK(1) whenever you set INFO(4)=1.
CC
CC  Do not change INFO(5), INFO(6), IWORK(1), or IWORK(2)
CC  unless you are going to restart the code.
CC
CC                 *** Following a completed task ***
CC  If
CC     IDID = 1, call the code again to continue the integration
CC                  another step in the direction of TOUT.
CC
CC     IDID = 2 or 3, define a new TOUT and call the code again.
CC                  TOUT must be different from T. You cannot change
CC                  the direction of integration without restarting.
CC
CC                 *** Following an interrupted task ***
CC               To show the code that you realize the task was
CC               interrupted and that you want to continue, you
CC               must take appropriate action and set INFO(1) = 1
CC  If
CC    IDID = -1, The code has taken about 500 steps.
CC                  If you want to continue, set INFO(1) = 1 and
CC                  call the code again. An additional 500 steps
CC                  will be allowed.
CC
CC    IDID = -2, The error tolerances RTOL, ATOL have been
CC                  increased to values the code estimates appropriate
CC                  for continuing. You may want to change them
CC                  yourself. If you are sure you want to continue
CC                  with relaxed error tolerances, set INFO(1)=1 and
CC                  call the code again.
CC
CC    IDID = -3, A solution component is zero and you set the
CC                  corresponding component of ATOL to zero. If you
CC                  are sure you want to continue, you must first
CC                  alter the error criterion to use positive values
CC                  for those components of ATOL corresponding to zero
CC                  solution components, then set INFO(1)=1 and call
CC                  the code again.
CC
CC    IDID = -4,-5  --- Cannot occur with this code.
CC
CC    IDID = -6, Repeated error test failures occurred on the
CC                  last attempted step in SDASSL. A singularity in the
CC                  solution may be present. If you are absolutely
CC                  certain you want to continue, you should restart
CC                  the integration. (Provide initial values of Y and
CC                  YPRIME which are consistent)
CC
CC    IDID = -7, Repeated convergence test failures occurred
CC                  on the last attempted step in SDASSL. An inaccurate
CC                  or ill-conditioned JACOBIAN may be the problem. If
CC                  you are absolutely certain you want to continue, you
CC                  should restart the integration.
CC
CC    IDID = -8, The matrix of partial derivatives is singular.
CC                  Some of your equations may be redundant.
CC                  SDASSL cannot solve the problem as stated.
CC                  It is possible that the redundant equations
CC                  could be removed, and then SDASSL could
CC                  solve the problem. It is also possible
CC                  that a solution to your problem either
CC                  does not exist or is not unique.
CC
CC    IDID = -9, SDASSL had multiple convergence test
CC                  failures, preceeded by multiple error
CC                  test failures, on the last attempted step.
CC                  It is possible that your problem
CC                  is ill-posed, and cannot be solved
CC                  using this code. Or, there may be a
CC                  discontinuity or a singularity in the
CC                  solution. If you are absolutely certain
CC                  you want to continue, you should restart
CC                  the integration.
CC
CC    IDID =-10, SDASSL had multiple convergence test failures
CC                  because IRES was equal to minus one.
CC                  If you are absolutely certain you want
CC                  to continue, you should restart the
CC                  integration.
CC
CC    IDID =-11, IRES=-2 was encountered, and control is being
CC                  returned to the calling program.
CC
CC    IDID =-12, SDASSL failed to compute the initial YPRIME.
CC                  This could happen because the initial
CC                  approximation to YPRIME was not very good, or
CC                  if a YPRIME consistent with the initial Y
CC                  does not exist. The problem could also be caused
CC                  by an inaccurate or singular iteration matrix.
CC
CC    IDID = -13,..,-32  --- Cannot occur with this code.
CC
CC
CC                 *** Following a terminated task ***
CC
CC  If IDID= -33, you cannot continue the solution of this problem.
CC                  An attempt to do so will result in your
CC                  run being terminated.
CC
CC
CC  ----- ERROR MESSAGES ---------------------------------------------
CC
CC      The SLATEC error print routine XERMSG is called in the event of
CC   unsuccessful completion of a task.  Most of these are treated as
CC   "recoverable errors", meaning that (unless the user has directed
CC   otherwise) control will be returned to the calling program for
CC   possible action after the message has been printed.
CC
CC   In the event of a negative value of IDID other than -33, an appro-
CC   priate message is printed and the "error number" printed by XERMSG
CC   is the value of IDID.  There are quite a number of illegal input
CC   errors that can lead to a returned value IDID=-33.  The conditions
CC   and their printed "error numbers" are as follows:
CC
CC   Error number       Condition
CC
CC        1       Some element of INFO vector is not zero or one.
CC        2       NEQ .le. 0
CC        3       MAXORD not in range.
CC        4       LRW is less than the required length for RWORK.
CC        5       LIW is less than the required length for IWORK.
CC        6       Some element of RTOL is .lt. 0
CC        7       Some element of ATOL is .lt. 0
CC        8       All elements of RTOL and ATOL are zero.
CC        9       INFO(4)=1 and TSTOP is behind TOUT.
CC       10       HMAX .lt. 0.0
CC       11       TOUT is behind T.
CC       12       INFO(8)=1 and H0=0.0
CC       13       Some element of WT is .le. 0.0
CC       14       TOUT is too close to T to start integration.
CC       15       INFO(4)=1 and TSTOP is behind T.
CC       16       --( Not used in this version )--
CC       17       ML illegal.  Either .lt. 0 or .gt. NEQ
CC       18       MU illegal.  Either .lt. 0 or .gt. NEQ
CC       19       TOUT = T.
CC
CC   If SDASSL is called again without any action taken to remove the
CC   cause of an unsuccessful return, XERMSG will be called with a fatal
CC   error flag, which will cause unconditional termination of the
CC   program.  There are two such fatal errors:
CC
CC   Error number -998:  The last step was terminated with a negative
CC       value of IDID other than -33, and no appropriate action was
CC       taken.
CC
CC   Error number -999:  The previous call was terminated because of
CC       illegal input (IDID=-33) and there is illegal input in the
CC       present call, as well.  (Suspect infinite loop.)
CC
CC  ----------------------------------------------------------------
CC
CC***REFERENCES  A DESCRIPTION OF DASSL: A DIFFERENTIAL/ALGEBRAIC
CC                 SYSTEM SOLVER, L. R. PETZOLD, SAND82-8637,
CC                 SANDIA NATIONAL LABORATORIES, SEPTEMBER 1982.
CC***ROUTINES CALLED  R1MACH, SDAINI, SDANRM, SDASTP, SDATRP, SDAWTS,
CC                    XERMSG
CC***REVISION HISTORY  (YYMMDD)
CC   830315  DATE WRITTEN
CC   880387  Code changes made.  All common statements have been
CC           replaced by a DATA statement, which defines pointers into
CC           RWORK, and PARAMETER statements which define pointers
CC           into IWORK.  As well the documentation has gone through
CC           grammatical changes.
CC   881005  The prologue has been changed to mixed case.
CC           The subordinate routines had revision dates changed to
CC           this date, although the documentation for these routines
CC           is all upper case.  No code changes.
CC   890511  Code changes made.  The DATA statement in the declaration
CC           section of SDASSL was replaced with a PARAMETER
CC           statement.  Also the statement S = 100.E0 was removed
CC           from the top of the Newton iteration in SDASTP.
CC           The subordinate routines had revision dates changed to
CC           this date.
CC   890517  The revision date syntax was replaced with the revision
CC           history syntax.  Also the "DECK" comment was added to
CC           the top of all subroutines.  These changes are consistent
CC           with new SLATEC guidelines.
CC           The subordinate routines had revision dates changed to
CC           this date.  No code changes.
CC   891013  Code changes made.
CC           Removed all occurrances of FLOAT.  All operations
CC           are now performed with "mixed-mode" arithmetic.
CC           Also, specific function names were replaced with generic
CC           function names to be consistent with new SLATEC guidelines.
CC           In particular:
CC              Replaced AMIN1 with MIN everywhere.
CC              Replaced MIN0 with MIN everywhere.
CC              Replaced AMAX1 with MAX everywhere.
CC              Replaced MAX0 with MAX everywhere.
CC           Also replaced REVISION DATE with REVISION HISTORY in all
CC           subordinate routines.
CC  901004  Miscellaneous changes to prologue to complete conversion
CC          to SLATEC 4.0 format.  No code changes.  (F.N.Fritsch)
CC  901009  Corrected GAMS classification code and converted subsidiary
CC          routines to 4.0 format.  No code changes.  (F.N.Fritsch)
CC  901010  Converted XERRWV calls to XERMSG calls.  (R.Clemens,AFWL)
CC  901019  Code changes made.
CC          Merged SLATEC 4.0 changes with previous changes made
CC          by C. Ulrich.  Below is a history of the changes made by
CC          C. Ulrich. (Changes in subsidiary routines are implied
CC          by this history)
CC          891228  Bug was found and repaired inside the SDASSL
CC                  and SDAINI routines.  SDAINI was incorrectly
CC                  returning the initial T with Y and YPRIME
CC                  computed at T+H.  The routine now returns T+H
CC                  rather than the initial T.
CC                  Cosmetic changes made to SDASTP.
CC          900904  Three modifications were made to fix a bug (inside
CC                  SDASSL) re interpolation for continuation calls and
CC                  cases where TN is very close to TSTOP:
CC
CC                  1) In testing for whether H is too large, just
CC                     compare H to (TSTOP - TN), rather than
CC                     (TSTOP - TN) * (1-4*UROUND), and set H to
CC                     TSTOP - TN.  This will force SDASTP to step
CC                     exactly to TSTOP under certain situations
CC                     (i.e. when H returned from SDASTP would otherwise
CC                     take TN beyond TSTOP).
CC
CC                  2) Inside the SDASTP loop, interpolate exactly to
CC                     TSTOP if TN is very close to TSTOP (rather than
CC                     interpolating to within roundoff of TSTOP).
CC
CC                  3) Modified IDID description for IDID = 2 to say
CC                     the solution is returned by stepping exactly to
CC                     TSTOP, rather than TOUT.  (In some cases the
CC                     solution is actually obtained by extrapolating
CC                     over a distance near unit roundoff to TSTOP,
CC                     but this small distance is deemed acceptable in
CC                     these circumstances.)
CC   901026  Added explicit declarations for all variables and minor
CC           cosmetic changes to prologue, removed unreferenced labels,
CC           and improved XERMSG calls.  (FNF)
CC   901030  Added ERROR MESSAGES section and reworked other sections to
CC           be of more uniform format.  (FNF)
CC   910624  Fixed minor bug related to HMAX (five lines ending
CC           in label 526).  (LRP)
CC
CC***END PROLOGUE  SDASSL
CC
CC**End
CC
CC     Declare arguments.
CC
C      INTEGER  NEQ, INFO(15), IDID, LRW, IWORK(*), LIW, IPAR(*)
C      REAL  T, Y(*), YPRIME(*), TOUT, RTOL(*), ATOL(*), RWORK(*),
C     *   RPAR(*), SIGN77
C      EXTERNAL  RES, JAC, SIGN77
CC
CC     Declare externals.
CC
C      EXTERNAL  R1MACH, SDAINI, SDANRM, SDASTP, SDATRP, SDAWTS, XERMSG
C      REAL  R1MACH, SDANRM
CC
CC     Declare local variables.
CC
C      INTEGER  I, ITEMP, LALPHA, LBETA, LCJ, LCJOLD, LCTF, LDELTA,
C     *   LENIW,LENPD,LENRW, LE, LETF, LGAMMA, LH, LHMAX, LHOLD, LIPVT,
C     *   LJCALC, LK, LKOLD, LIWM, LML, LMTYPE, LMU, LMXORD, LNJE, LNPD,
C     *   LNRE, LNS, LNST, LNSTL, LPD, LPHASE, LPHI, LPSI, LROUND, LS,
C     *   LSIGMA,LTN,LTSTOP, LWM, LWT, MBAND, MSAVE, MXORD, NPD, NTEMP,
C     *   NZFLG
C      REAL  ATOLI, H, HMAX, HMIN, HO, R, RH, RTOLI, TDIST, TN, TNEXT,
C     *   TSTOP, UROUND, YPNORM
C      LOGICAL  DONE
CC       Auxiliary variables for conversion of values to be included in
CC       error messages.
C      CHARACTER*8  XERN1, XERN2
C      CHARACTER*16 XERN3, XERN4
CC
CC     SET POINTERS INTO IWORK
C      PARAMETER (LML=1, LMU=2, LMXORD=3, LMTYPE=4, LNST=11,
C     *  LNRE=12, LNJE=13, LETF=14, LCTF=15, LNPD=16,
C     *  LIPVT=21, LJCALC=5, LPHASE=6, LK=7, LKOLD=8,
C     *  LNS=9, LNSTL=10, LIWM=1)
CC
CC     SET RELATIVE OFFSET INTO RWORK
C      PARAMETER (NPD=1)
CC
CC     SET POINTERS INTO RWORK
C      PARAMETER (LTSTOP=1, LHMAX=2, LH=3, LTN=4,
C     *  LCJ=5, LCJOLD=6, LHOLD=7, LS=8, LROUND=9,
C     *  LALPHA=11, LBETA=17, LGAMMA=23,
C     *  LPSI=29, LSIGMA=35, LDELTA=41)
CC
CC***FIRST EXECUTABLE STATEMENT  SDASSL
C      IF(INFO(1).NE.0)GO TO 100
CC
CC------------------------------------------------------------------
CC     THIS BLOCK IS EXECUTED FOR THE INITIAL CALL ONLY.
CC     IT CONTAINS CHECKING OF INPUTS AND INITIALIZATIONS.
CC------------------------------------------------------------------
CC
CC     FIRST CHECK INFO ARRAY TO MAKE SURE ALL ELEMENTS OF INFO
CC     ARE EITHER ZERO OR ONE.
C      DO 10 I=2,11
C         IF(INFO(I).NE.0.AND.INFO(I).NE.1)GO TO 701
C10       CONTINUE
CC
C      IF(NEQ.LE.0)GO TO 702
CC
CC     CHECK AND COMPUTE MAXIMUM ORDER
C      MXORD=5
C      IF(INFO(9).EQ.0)GO TO 20
C         MXORD=IWORK(LMXORD)
C         IF(MXORD.LT.1.OR.MXORD.GT.5)GO TO 703
C20       IWORK(LMXORD)=MXORD
CC
CC     COMPUTE MTYPE,LENPD,LENRW.CHECK ML AND MU.
C      IF(INFO(6).NE.0)GO TO 40
C         LENPD=NEQ**2
C         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
C         IF(INFO(5).NE.0)GO TO 30
C            IWORK(LMTYPE)=2
C            GO TO 60
C30          IWORK(LMTYPE)=1
C            GO TO 60
C40    IF(IWORK(LML).LT.0.OR.IWORK(LML).GE.NEQ)GO TO 717
C      IF(IWORK(LMU).LT.0.OR.IWORK(LMU).GE.NEQ)GO TO 718
C      LENPD=(2*IWORK(LML)+IWORK(LMU)+1)*NEQ
C      IF(INFO(5).NE.0)GO TO 50
C         IWORK(LMTYPE)=5
C         MBAND=IWORK(LML)+IWORK(LMU)+1
C         MSAVE=(NEQ/MBAND)+1
C         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD+2*MSAVE
C         GO TO 60
C50       IWORK(LMTYPE)=4
C         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
CC
CC     CHECK LENGTHS OF RWORK AND IWORK
C60    LENIW=20+NEQ
C      IWORK(LNPD)=LENPD
C      IF(LRW.LT.LENRW)GO TO 704
C      IF(LIW.LT.LENIW)GO TO 705
CC
CC     CHECK TO SEE THAT TOUT IS DIFFERENT FROM T
C      IF(TOUT .EQ. T)GO TO 719
CC
CC     CHECK HMAX
C      IF(INFO(7).EQ.0)GO TO 70
C         HMAX=RWORK(LHMAX)
C         IF(HMAX.LE.0.0E0)GO TO 710
C70    CONTINUE
CC
CC     INITIALIZE COUNTERS
C      IWORK(LNST)=0
C      IWORK(LNRE)=0
C      IWORK(LNJE)=0
CC
C      IWORK(LNSTL)=0
C      IDID=1
C      GO TO 200
CC
CC------------------------------------------------------------------
CC     THIS BLOCK IS FOR CONTINUATION CALLS
CC     ONLY. HERE WE CHECK INFO(1),AND IF THE
CC     LAST STEP WAS INTERRUPTED WE CHECK WHETHER
CC     APPROPRIATE ACTION WAS TAKEN.
CC------------------------------------------------------------------
CC
C100   CONTINUE
C      IF(INFO(1).EQ.1)GO TO 110
C      IF(INFO(1).NE.-1)GO TO 701
CC
CC     IF WE ARE HERE, THE LAST STEP WAS INTERRUPTED
CC     BY AN ERROR CONDITION FROM SDASTP,AND
CC     APPROPRIATE ACTION WAS NOT TAKEN. THIS
CC     IS A FATAL ERROR.
C      WRITE (XERN1, '(I8)') IDID
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'THE LAST STEP TERMINATED WITH A NEGATIVE VALUE OF IDID = ' //
C     *   XERN1 // ' AND NO APPROPRIATE ACTION WAS TAKEN.  ' //
C     *   'RUN TERMINATED', -998, 2)
C      RETURN
C110   CONTINUE
C      IWORK(LNSTL)=IWORK(LNST)
CC
CC------------------------------------------------------------------
CC     THIS BLOCK IS EXECUTED ON ALL CALLS.
CC     THE ERROR TOLERANCE PARAMETERS ARE
CC     CHECKED, AND THE WORK ARRAY POINTERS
CC     ARE SET.
CC------------------------------------------------------------------
CC
C200   CONTINUE
CC     CHECK RTOL,ATOL
C      NZFLG=0
C      RTOLI=RTOL(1)
C      ATOLI=ATOL(1)
C      DO 210 I=1,NEQ
C         IF(INFO(2).EQ.1)RTOLI=RTOL(I)
C         IF(INFO(2).EQ.1)ATOLI=ATOL(I)
C         IF(RTOLI.GT.0.0E0.OR.ATOLI.GT.0.0E0)NZFLG=1
C         IF(RTOLI.LT.0.0E0)GO TO 706
C         IF(ATOLI.LT.0.0E0)GO TO 707
C210      CONTINUE
C      IF(NZFLG.EQ.0)GO TO 708
CC
CC     SET UP RWORK STORAGE.IWORK STORAGE IS FIXED
CC     IN DATA STATEMENT.
C      LE=LDELTA+NEQ
C      LWT=LE+NEQ
C      LPHI=LWT+NEQ
C      LPD=LPHI+(IWORK(LMXORD)+1)*NEQ
C      LWM=LPD
C      NTEMP=NPD+IWORK(LNPD)
C      IF(INFO(1).EQ.1)GO TO 400
CC
CC------------------------------------------------------------------
CC     THIS BLOCK IS EXECUTED ON THE INITIAL CALL
CC     ONLY. SET THE INITIAL STEP SIZE, AND
CC     THE ERROR WEIGHT VECTOR, AND PHI.
CC     COMPUTE INITIAL YPRIME, IF NECESSARY.
CC------------------------------------------------------------------
CC
C      TN=T
C      IDID=1
CC
CC     SET ERROR WEIGHT VECTOR WT
C      CALL SDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
C      DO 305 I = 1,NEQ
C         IF(RWORK(LWT+I-1).LE.0.0E0) GO TO 713
C305      CONTINUE
CC
CC     COMPUTE UNIT ROUNDOFF AND HMIN
C      UROUND = R1MACH(4)
C      RWORK(LROUND) = UROUND
C      HMIN = 4.0E0*UROUND*MAX(ABS(T),ABS(TOUT))
CC
CC     CHECK INITIAL INTERVAL TO SEE THAT IT IS LONG ENOUGH
C      TDIST = ABS(TOUT - T)
C      IF(TDIST .LT. HMIN) GO TO 714
CC
CC     CHECK HO, IF THIS WAS INPUT
C      IF (INFO(8) .EQ. 0) GO TO 310
C         HO = RWORK(LH)
C         IF ((TOUT - T)*HO .LT. 0.0E0) GO TO 711
C         IF (HO .EQ. 0.0E0) GO TO 712
C         GO TO 320
C310    CONTINUE
CC
CC     COMPUTE INITIAL STEPSIZE, TO BE USED BY EITHER
CC     SDASTP OR SDAINI, DEPENDING ON INFO(11)
C      HO = 0.001E0*TDIST
C      YPNORM = SDANRM(NEQ,YPRIME,RWORK(LWT),RPAR,IPAR)
C      IF (YPNORM .GT. 0.5E0/HO) HO = 0.5E0/YPNORM
C      HO = SIGN77(HO,TOUT-T)
CC     ADJUST HO IF NECESSARY TO MEET HMAX BOUND
C320   IF (INFO(7) .EQ. 0) GO TO 330
C         RH = ABS(HO)/RWORK(LHMAX)
C         IF (RH .GT. 1.0E0) HO = HO/RH
CC     COMPUTE TSTOP, IF APPLICABLE
C330   IF (INFO(4) .EQ. 0) GO TO 340
C         TSTOP = RWORK(LTSTOP)
C         IF ((TSTOP - T)*HO .LT. 0.0E0) GO TO 715
C         IF ((T + HO - TSTOP)*HO .GT. 0.0E0) HO = TSTOP - T
C         IF ((TSTOP - TOUT)*HO .LT. 0.0E0) GO TO 709
CC
CC     COMPUTE INITIAL DERIVATIVE, UPDATING TN AND Y, IF APPLICABLE
C340   IF (INFO(11) .EQ. 0) GO TO 350
C      CALL SDAINI(TN,Y,YPRIME,NEQ,
C     *  RES,JAC,HO,RWORK(LWT),IDID,RPAR,IPAR,
C     *  RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
C     *  RWORK(LWM),IWORK(LIWM),HMIN,RWORK(LROUND),
C     *  INFO(10),NTEMP)
C      IF (IDID .LT. 0) GO TO 390
CC
CC     LOAD H WITH HO.  STORE H IN RWORK(LH)
C350   H = HO
C      RWORK(LH) = H
CC
CC     LOAD Y AND H*YPRIME INTO PHI(*,1) AND PHI(*,2)
C      ITEMP = LPHI + NEQ
C      DO 370 I = 1,NEQ
C         RWORK(LPHI + I - 1) = Y(I)
C370      RWORK(ITEMP + I - 1) = H*YPRIME(I)
CC
C390   GO TO 500
CC
CC-------------------------------------------------------
CC     THIS BLOCK IS FOR CONTINUATION CALLS ONLY. ITS
CC     PURPOSE IS TO CHECK STOP CONDITIONS BEFORE
CC     TAKING A STEP.
CC     ADJUST H IF NECESSARY TO MEET HMAX BOUND
CC-------------------------------------------------------
CC
C400   CONTINUE
C      UROUND=RWORK(LROUND)
C      DONE = .FALSE.
C      TN=RWORK(LTN)
C      H=RWORK(LH)
C      IF(INFO(7) .EQ. 0) GO TO 410
C         RH = ABS(H)/RWORK(LHMAX)
C         IF(RH .GT. 1.0E0) H = H/RH
C410   CONTINUE
C      IF(T .EQ. TOUT) GO TO 719
C      IF((T - TOUT)*H .GT. 0.0E0) GO TO 711
C      IF(INFO(4) .EQ. 1) GO TO 430
C      IF(INFO(3) .EQ. 1) GO TO 420
C      IF((TN-TOUT)*H.LT.0.0E0)GO TO 490
C      CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
C     *  RWORK(LPHI),RWORK(LPSI))
C      T=TOUT
C      IDID = 3
C      DONE = .TRUE.
C      GO TO 490
C420   IF((TN-T)*H .LE. 0.0E0) GO TO 490
C      IF((TN - TOUT)*H .GT. 0.0E0) GO TO 425
C      CALL SDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
C     *  RWORK(LPHI),RWORK(LPSI))
C      T = TN
C      IDID = 1
C      DONE = .TRUE.
C      GO TO 490
C425   CONTINUE
C      CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
C     *  RWORK(LPHI),RWORK(LPSI))
C      T = TOUT
C      IDID = 3
C      DONE = .TRUE.
C      GO TO 490
C430   IF(INFO(3) .EQ. 1) GO TO 440
C      TSTOP=RWORK(LTSTOP)
C      IF((TN-TSTOP)*H.GT.0.0E0) GO TO 715
C      IF((TSTOP-TOUT)*H.LT.0.0E0)GO TO 709
C      IF((TN-TOUT)*H.LT.0.0E0)GO TO 450
C      CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
C     *   RWORK(LPHI),RWORK(LPSI))
C      T=TOUT
C      IDID = 3
C      DONE = .TRUE.
C      GO TO 490
C440   TSTOP = RWORK(LTSTOP)
C      IF((TN-TSTOP)*H .GT. 0.0E0) GO TO 715
C      IF((TSTOP-TOUT)*H .LT. 0.0E0) GO TO 709
C      IF((TN-T)*H .LE. 0.0E0) GO TO 450
C      IF((TN - TOUT)*H .GT. 0.0E0) GO TO 445
C      CALL SDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
C     *  RWORK(LPHI),RWORK(LPSI))
C      T = TN
C      IDID = 1
C      DONE = .TRUE.
C      GO TO 490
C445   CONTINUE
C      CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
C     *  RWORK(LPHI),RWORK(LPSI))
C      T = TOUT
C      IDID = 3
C      DONE = .TRUE.
C      GO TO 490
C450   CONTINUE
CC     CHECK WHETHER WE ARE WITHIN ROUNDOFF OF TSTOP
C      IF(ABS(TN-TSTOP).GT.100.0E0*UROUND*
C     *   (ABS(TN)+ABS(H)))GO TO 460
C      CALL SDATRP(TN,TSTOP,Y,YPRIME,NEQ,IWORK(LKOLD),
C     *  RWORK(LPHI),RWORK(LPSI))
C      IDID=2
C      T=TSTOP
C      DONE = .TRUE.
C      GO TO 490
C460   TNEXT=TN+H
C      IF((TNEXT-TSTOP)*H.LE.0.0E0)GO TO 490
C      H=TSTOP-TN
C      RWORK(LH)=H
CC
C490   IF (DONE) GO TO 580
CC
CC-------------------------------------------------------
CC     THE NEXT BLOCK CONTAINS THE CALL TO THE
CC     ONE-STEP INTEGRATOR SDASTP.
CC     THIS IS A LOOPING POINT FOR THE INTEGRATION STEPS.
CC     CHECK FOR TOO MANY STEPS.
CC     UPDATE WT.
CC     CHECK FOR TOO MUCH ACCURACY REQUESTED.
CC     COMPUTE MINIMUM STEPSIZE.
CC-------------------------------------------------------
CC
C500   CONTINUE
CC     CHECK FOR FAILURE TO COMPUTE INITIAL YPRIME
C      IF (IDID .EQ. -12) GO TO 527
CC
CC     CHECK FOR TOO MANY STEPS
C      IF((IWORK(LNST)-IWORK(LNSTL)).LT.500)
C     *   GO TO 510
C           IDID=-1
C           GO TO 527
CC
CC     UPDATE WT
C510   CALL SDAWTS(NEQ,INFO(2),RTOL,ATOL,RWORK(LPHI),
C     *  RWORK(LWT),RPAR,IPAR)
C      DO 520 I=1,NEQ
C         IF(RWORK(I+LWT-1).GT.0.0E0)GO TO 520
C           IDID=-3
C           GO TO 527
C520   CONTINUE
CC
CC     TEST FOR TOO MUCH ACCURACY REQUESTED.
C      R=SDANRM(NEQ,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)*
C     *   100.0E0*UROUND
C      IF(R.LE.1.0E0)GO TO 525
CC     MULTIPLY RTOL AND ATOL BY R AND RETURN
C      IF(INFO(2).EQ.1)GO TO 523
C           RTOL(1)=R*RTOL(1)
C           ATOL(1)=R*ATOL(1)
C           IDID=-2
C           GO TO 527
C523   DO 524 I=1,NEQ
C           RTOL(I)=R*RTOL(I)
C524        ATOL(I)=R*ATOL(I)
C      IDID=-2
C      GO TO 527
C525   CONTINUE
CC
CC     COMPUTE MINIMUM STEPSIZE
C      HMIN=4.0E0*UROUND*MAX(ABS(TN),ABS(TOUT))
CC
CC     TEST H VS. HMAX
C      IF (INFO(7) .EQ. 0) GO TO 526
C         RH = ABS(H)/RWORK(LHMAX)
C         IF (RH .GT. 1.0E0) H = H/RH
C526   CONTINUE
CC
C      CALL SDASTP(TN,Y,YPRIME,NEQ,
C     *   RES,JAC,H,RWORK(LWT),INFO(1),IDID,RPAR,IPAR,
C     *   RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
C     *   RWORK(LWM),IWORK(LIWM),
C     *   RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),
C     *   RWORK(LPSI),RWORK(LSIGMA),
C     *   RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),
C     *   RWORK(LS),HMIN,RWORK(LROUND),
C     *   IWORK(LPHASE),IWORK(LJCALC),IWORK(LK),
C     *   IWORK(LKOLD),IWORK(LNS),INFO(10),NTEMP)
C527   IF(IDID.LT.0)GO TO 600
CC
CC--------------------------------------------------------
CC     THIS BLOCK HANDLES THE CASE OF A SUCCESSFUL RETURN
CC     FROM SDASTP (IDID=1).  TEST FOR STOP CONDITIONS.
CC--------------------------------------------------------
CC
C      IF(INFO(4).NE.0)GO TO 540
C           IF(INFO(3).NE.0)GO TO 530
C             IF((TN-TOUT)*H.LT.0.0E0)GO TO 500
C             CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,
C     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
C             IDID=3
C             T=TOUT
C             GO TO 580
C530          IF((TN-TOUT)*H.GE.0.0E0)GO TO 535
C             T=TN
C             IDID=1
C             GO TO 580
C535          CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,
C     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
C             IDID=3
C             T=TOUT
C             GO TO 580
C540   IF(INFO(3).NE.0)GO TO 550
C      IF((TN-TOUT)*H.LT.0.0E0)GO TO 542
C         CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,
C     *     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
C         T=TOUT
C         IDID=3
C         GO TO 580
C542   IF(ABS(TN-TSTOP).LE.100.0E0*UROUND*
C     *   (ABS(TN)+ABS(H)))GO TO 545
C      TNEXT=TN+H
C      IF((TNEXT-TSTOP)*H.LE.0.0E0)GO TO 500
C      H=TSTOP-TN
C      GO TO 500
C545   CALL SDATRP(TN,TSTOP,Y,YPRIME,NEQ,
C     *  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
C      IDID=2
C      T=TSTOP
C      GO TO 580
C550   IF((TN-TOUT)*H.GE.0.0E0)GO TO 555
C      IF(ABS(TN-TSTOP).LE.100.0E0*UROUND*(ABS(TN)+ABS(H)))GO TO 552
C      T=TN
C      IDID=1
C      GO TO 580
C552   CALL SDATRP(TN,TSTOP,Y,YPRIME,NEQ,
C     *  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
C      IDID=2
C      T=TSTOP
C      GO TO 580
C555   CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,
C     *   IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
C      T=TOUT
C      IDID=3
C      GO TO 580
CC
CC--------------------------------------------------------
CC     ALL SUCCESSFUL RETURNS FROM SDASSL ARE MADE FROM
CC     THIS BLOCK.
CC--------------------------------------------------------
CC
C580   CONTINUE
C      RWORK(LTN)=TN
C      RWORK(LH)=H
C      RETURN
CC
CC------------------------------------------------------------------
CC     THIS BLOCK HANDLES ALL UNSUCCESSFUL
CC     RETURNS OTHER THAN FOR ILLEGAL INPUT.
CC------------------------------------------------------------------
CC
C600   CONTINUE
C      ITEMP=-IDID
C      GO TO (610,620,630,690,690,640,650,660,670,675,
C     *  680,685), ITEMP
CC
CC     THE MAXIMUM NUMBER OF STEPS WAS TAKEN BEFORE
CC     REACHING TOUT
C610   WRITE (XERN3, '(1P,E15.6)') TN
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'AT CURRENT T = ' // XERN3 // ' 500 STEPS TAKEN ON THIS ' //
C     *   'CALL BEFORE REACHING TOUT', IDID, 1)
C      GO TO 690
CC
CC     TOO MUCH ACCURACY FOR MACHINE PRECISION
C620   WRITE (XERN3, '(1P,E15.6)') TN
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'AT T = ' // XERN3 // ' TOO MUCH ACCURACY REQUESTED FOR ' //
C     *   'PRECISION OF MACHINE. RTOL AND ATOL WERE INCREASED TO ' //
C     *   'APPROPRIATE VALUES', IDID, 1)
C      GO TO 690
CC
CC     WT(I) .LE. 0.0 FOR SOME I (NOT AT START OF PROBLEM)
C630   WRITE (XERN3, '(1P,E15.6)') TN
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *'AT T = ' // XERN3 // ' SOME ELEMENT OF WT HAS BECOME .LE. ' //
C     *   '0.0', IDID, 1)
C      GO TO 690
CC
CC     ERROR TEST FAILED REPEATEDLY OR WITH H=HMIN
C640   WRITE (XERN3, '(1P,E15.6)') TN
C      WRITE (XERN4, '(1P,E15.6)') H
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
C     *   ' THE ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN',
C     *   IDID, 1)
C      GO TO 690
CC
CC     CORRECTOR CONVERGENCE FAILED REPEATEDLY OR WITH H=HMIN
C650   WRITE (XERN3, '(1P,E15.6)') TN
C      WRITE (XERN4, '(1P,E15.6)') H
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
C     *   ' THE CORRECTOR FAILED TO CONVERGE REPEATEDLY OR WITH ' //
C     *   'ABS(H)=HMIN', IDID, 1)
C      GO TO 690
CC
CC     THE ITERATION MATRIX IS SINGULAR
C660   WRITE (XERN3, '(1P,E15.6)') TN
C      WRITE (XERN4, '(1P,E15.6)') H
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
C     *   ' THE ITERATION MATRIX IS SINGULAR', IDID, 1)
C      GO TO 690
CC
CC     CORRECTOR FAILURE PRECEEDED BY ERROR TEST FAILURES.
C670   WRITE (XERN3, '(1P,E15.6)') TN
C      WRITE (XERN4, '(1P,E15.6)') H
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
C     *   ' THE CORRECTOR COULD NOT CONVERGE.  ALSO, THE ERROR TEST ' //
C     *   'FAILED REPEATEDLY.', IDID, 1)
C      GO TO 690
CC
CC     CORRECTOR FAILURE BECAUSE IRES = -1
C675   WRITE (XERN3, '(1P,E15.6)') TN
C      WRITE (XERN4, '(1P,E15.6)') H
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
C     *   ' THE CORRECTOR COULD NOT CONVERGE BECAUSE IRES WAS EQUAL ' //
C     *   'TO MINUS ONE', IDID, 1)
C      GO TO 690
CC
CC     FAILURE BECAUSE IRES = -2
C680   WRITE (XERN3, '(1P,E15.6)') TN
C      WRITE (XERN4, '(1P,E15.6)') H
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
C     *   ' IRES WAS EQUAL TO MINUS TWO', IDID, 1)
C      GO TO 690
CC
CC     FAILED TO COMPUTE INITIAL YPRIME
C685   WRITE (XERN3, '(1P,E15.6)') TN
C      WRITE (XERN4, '(1P,E15.6)') HO
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
C     *   ' THE INITIAL YPRIME COULD NOT BE COMPUTED', IDID, 1)
C      GO TO 690
CC
C690   CONTINUE
C      INFO(1)=-1
C      T=TN
C      RWORK(LTN)=TN
C      RWORK(LH)=H
C      RETURN
CC
CC------------------------------------------------------------------
CC     THIS BLOCK HANDLES ALL ERROR RETURNS DUE
CC     TO ILLEGAL INPUT, AS DETECTED BEFORE CALLING
CC     SDASTP. FIRST THE ERROR MESSAGE ROUTINE IS
CC     CALLED. IF THIS HAPPENS TWICE IN
CC     SUCCESSION, EXECUTION IS TERMINATED
CC
CC------------------------------------------------------------------
C701   CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR ONE', 1, 1)
C      GO TO 750
CC
C702   WRITE (XERN1, '(I8)') NEQ
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'NEQ = ' // XERN1 // ' .LE. 0', 2, 1)
C      GO TO 750
CC
C703   WRITE (XERN1, '(I8)') MXORD
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'MAXORD = ' // XERN1 // ' NOT IN RANGE', 3, 1)
C      GO TO 750
CC
C704   WRITE (XERN1, '(I8)') LENRW
C      WRITE (XERN2, '(I8)') LRW
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'RWORK LENGTH NEEDED, LENRW = ' // XERN1 //
C     *   ', EXCEEDS LRW = ' // XERN2, 4, 1)
C      GO TO 750
CC
C705   WRITE (XERN1, '(I8)') LENIW
C      WRITE (XERN2, '(I8)') LIW
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'IWORK LENGTH NEEDED, LENIW = ' // XERN1 //
C     *   ', EXCEEDS LIW = ' // XERN2, 5, 1)
C      GO TO 750
CC
C706   CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'SOME ELEMENT OF RTOL IS .LT. 0', 6, 1)
C      GO TO 750
CC
C707   CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'SOME ELEMENT OF ATOL IS .LT. 0', 7, 1)
C      GO TO 750
CC
C708   CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'ALL ELEMENTS OF RTOL AND ATOL ARE ZERO', 8, 1)
C      GO TO 750
CC
C709   WRITE (XERN3, '(1P,E15.6)') TSTOP
C      WRITE (XERN4, '(1P,E15.6)') TOUT
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'INFO(4) = 1 AND TSTOP = ' // XERN3 // ' BEHIND TOUT = ' //
C     *   XERN4, 9, 1)
C      GO TO 750
CC
C710   WRITE (XERN3, '(1P,E15.6)') HMAX
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'HMAX = ' // XERN3 // ' .LT. 0.0', 10, 1)
C      GO TO 750
CC
C711   WRITE (XERN3, '(1P,E15.6)') TOUT
C      WRITE (XERN4, '(1P,E15.6)') T
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'TOUT = ' // XERN3 // ' BEHIND T = ' // XERN4, 11, 1)
C      GO TO 750
CC
C712   CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'INFO(8)=1 AND H0=0.0', 12, 1)
C      GO TO 750
CC
C713   CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'SOME ELEMENT OF WT IS .LE. 0.0', 13, 1)
C      GO TO 750
CC
C714   WRITE (XERN3, '(1P,E15.6)') TOUT
C      WRITE (XERN4, '(1P,E15.6)') T
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'TOUT = ' // XERN3 // ' TOO CLOSE TO T = ' // XERN4 //
C     *   ' TO START INTEGRATION', 14, 1)
C      GO TO 750
CC
C715   WRITE (XERN3, '(1P,E15.6)') TSTOP
C      WRITE (XERN4, '(1P,E15.6)') T
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'INFO(4)=1 AND TSTOP = ' // XERN3 // ' BEHIND T = ' // XERN4,
C     *   15, 1)
C      GO TO 750
CC
C717   WRITE (XERN1, '(I8)') IWORK(LML)
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'ML = ' // XERN1 // ' ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ',
C     *   17, 1)
C      GO TO 750
CC
C718   WRITE (XERN1, '(I8)') IWORK(LMU)
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *   'MU = ' // XERN1 // ' ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ',
C     *   18, 1)
C      GO TO 750
CC
C719   WRITE (XERN3, '(1P,E15.6)') TOUT
C      CALL XERMSG ('SLATEC', 'SDASSL',
C     *  'TOUT = T = ' // XERN3, 19, 1)
C      GO TO 750
CC
C750   IDID=-33
C      IF(INFO(1).EQ.-1) THEN
C         CALL XERMSG ('SLATEC', 'SDASSL',
C     *      'REPEATED OCCURRENCES OF ILLEGAL INPUT$$' //
C     *      'RUN TERMINATED. APPARENT INFINITE LOOP', -999, 2)
C      ENDIF
CC
C      INFO(1)=-1
C      RETURN
CC------END OF SUBROUTINE SDASSL------------------------------------
C      END
C      SUBROUTINE SDASTP (X, Y, YPRIME, NEQ, RES, JAC, H, WT, JSTART,
C     +   IDID, RPAR, IPAR, PHI, DELTA, E, WM, IWM, ALPHA, BETA, GAMMA,
C     +   PSI, SIGMA, CJ, CJOLD, HOLD, S, HMIN, UROUND, IPHASE, JCALC,
C     +   K, KOLD, NS, NONNEG, NTEMP)
CC***BEGIN PROLOGUE  SDASTP
CC***SUBSIDIARY
CC***PURPOSE  Perform one step of the SDASSL integration.
CC***LIBRARY   SLATEC (DASSL)
CC***TYPE      SINGLE PRECISION (SDASTP-S, DDASTP-D)
CC***AUTHOR  PETZOLD, LINDA R., (LLNL)
CC***DESCRIPTION
CC------------------------------------------------------------------
CC     SDASTP SOLVES A SYSTEM OF DIFFERENTIAL/
CC     ALGEBRAIC EQUATIONS OF THE FORM
CC     G(X,Y,YPRIME) = 0,  FOR ONE STEP (NORMALLY
CC     FROM X TO X+H).
CC
CC     THE METHODS USED ARE MODIFIED DIVIDED
CC     DIFFERENCE,FIXED LEADING COEFFICIENT
CC     FORMS OF BACKWARD DIFFERENTIATION
CC     FORMULAS. THE CODE ADJUSTS THE STEPSIZE
CC     AND ORDER TO CONTROL THE LOCAL ERROR PER
CC     STEP.
CC
CC
CC     THE PARAMETERS REPRESENT
CC     X  --        INDEPENDENT VARIABLE
CC     Y  --        SOLUTION VECTOR AT X
CC     YPRIME --    DERIVATIVE OF SOLUTION VECTOR
CC                  AFTER SUCCESSFUL STEP
CC     NEQ --       NUMBER OF EQUATIONS TO BE INTEGRATED
CC     RES --       EXTERNAL USER-SUPPLIED SUBROUTINE
CC                  TO EVALUATE THE RESIDUAL.  THE CALL IS
CC                  CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
CC                  X,Y,YPRIME ARE INPUT.  DELTA IS OUTPUT.
CC                  ON INPUT, IRES=0.  RES SHOULD ALTER IRES ONLY
CC                  IF IT ENCOUNTERS AN ILLEGAL VALUE OF Y OR A
CC                  STOP CONDITION.  SET IRES=-1 IF AN INPUT VALUE
CC                  OF Y IS ILLEGAL, AND SDASTP WILL TRY TO SOLVE
CC                  THE PROBLEM WITHOUT GETTING IRES = -1.  IF
CC                  IRES=-2, SDASTP RETURNS CONTROL TO THE CALLING
CC                  PROGRAM WITH IDID = -11.
CC     JAC --       EXTERNAL USER-SUPPLIED ROUTINE TO EVALUATE
CC                  THE ITERATION MATRIX (THIS IS OPTIONAL)
CC                  THE CALL IS OF THE FORM
CC                  CALL JAC(X,Y,YPRIME,PD,CJ,RPAR,IPAR)
CC                  PD IS THE MATRIX OF PARTIAL DERIVATIVES,
CC                  PD=DG/DY+CJ*DG/DYPRIME
CC     H --         APPROPRIATE STEP SIZE FOR NEXT STEP.
CC                  NORMALLY DETERMINED BY THE CODE
CC     WT --        VECTOR OF WEIGHTS FOR ERROR CRITERION.
CC     JSTART --    INTEGER VARIABLE SET 0 FOR
CC                  FIRST STEP, 1 OTHERWISE.
CC     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS:
CC                  IDID= 1 -- THE STEP WAS COMPLETED SUCCESSFULLY
CC                  IDID=-6 -- THE ERROR TEST FAILED REPEATEDLY
CC                  IDID=-7 -- THE CORRECTOR COULD NOT CONVERGE
CC                  IDID=-8 -- THE ITERATION MATRIX IS SINGULAR
CC                  IDID=-9 -- THE CORRECTOR COULD NOT CONVERGE.
CC                             THERE WERE REPEATED ERROR TEST
CC                             FAILURES ON THIS STEP.
CC                  IDID=-10-- THE CORRECTOR COULD NOT CONVERGE
CC                             BECAUSE IRES WAS EQUAL TO MINUS ONE
CC                  IDID=-11-- IRES EQUAL TO -2 WAS ENCOUNTERED,
CC                             AND CONTROL IS BEING RETURNED TO
CC                             THE CALLING PROGRAM
CC     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS THAT
CC                  ARE USED FOR COMMUNICATION BETWEEN THE
CC                  CALLING PROGRAM AND EXTERNAL USER ROUTINES
CC                  THEY ARE NOT ALTERED BY SDASTP
CC     PHI --       ARRAY OF DIVIDED DIFFERENCES USED BY
CC                  SDASTP. THE LENGTH IS NEQ*(K+1),WHERE
CC                  K IS THE MAXIMUM ORDER
CC     DELTA,E --   WORK VECTORS FOR SDASTP OF LENGTH NEQ
CC     WM,IWM --    REAL AND INTEGER ARRAYS STORING
CC                  MATRIX INFORMATION SUCH AS THE MATRIX
CC                  OF PARTIAL DERIVATIVES,PERMUTATION
CC                  VECTOR,AND VARIOUS OTHER INFORMATION.
CC
CC     THE OTHER PARAMETERS ARE INFORMATION
CC     WHICH IS NEEDED INTERNALLY BY SDASTP TO
CC     CONTINUE FROM STEP TO STEP.
CC
CC------------------------------------------------------------------
CC***ROUTINES CALLED  SDAJAC, SDANRM, SDASLV, SDATRP
CC***REVISION HISTORY  (YYMMDD)
CC   830315  DATE WRITTEN
CC   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
CC   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
CC   901026  Added explicit declarations for all variables and minor
CC           cosmetic changes to prologue.  (FNF)
CC***END PROLOGUE  SDASTP
CC
C      INTEGER  NEQ, JSTART, IDID, IPAR(*), IWM(*), IPHASE, JCALC, K,
C     *   KOLD, NS, NONNEG, NTEMP
C      REAL  X,Y(*),YPRIME(*), H, WT(*), RPAR(*), PHI(NEQ,*), DELTA(*),
C     * E(*), WM(*), ALPHA(*), BETA(*), GAMMA(*), PSI(*), SIGMA(*), CJ,
C     *   CJOLD, HOLD, S, HMIN, UROUND
C      EXTERNAL  RES, JAC
CC
C      EXTERNAL  SDAJAC, SDANRM, SDASLV, SDATRP
C      REAL  SDANRM
CC
C      INTEGER  I, IER, IRES, J, J1, KDIFF, KM1, KNEW, KP1, KP2, LCTF,
C     *   LETF, LMXORD, LNJE, LNRE, LNST, M, MAXIT, NCF, NEF, NSF, NSP1
C      REAL  ALPHA0, ALPHAS, CJLAST, CK, DELNRM, ENORM, ERK, ERKM1,
C     *   ERKM2, ERKP1, ERR, EST, HNEW, OLDNRM, PNORM, R, RATE, TEMP1,
C     *   TEMP2, TERK, TERKM1, TERKM2, TERKP1, XOLD, XRATE
C      LOGICAL  CONVGD
CC
C      PARAMETER (LMXORD=3)
C      PARAMETER (LNST=11)
C      PARAMETER (LNRE=12)
C      PARAMETER (LNJE=13)
C      PARAMETER (LETF=14)
C      PARAMETER (LCTF=15)
CC
C      DATA MAXIT/4/
C      DATA XRATE/0.25E0/
CC
CC
CC
CC
CC
CC------------------------------------------------------------------
CC     BLOCK 1.
CC     INITIALIZE. ON THE FIRST CALL,SET
CC     THE ORDER TO 1 AND INITIALIZE
CC     OTHER VARIABLES.
CC------------------------------------------------------------------
CC
CC     INITIALIZATIONS FOR ALL CALLS
CC***FIRST EXECUTABLE STATEMENT  SDASTP
C      IDID=1
C      XOLD=X
C      NCF=0
C      NSF=0
C      NEF=0
C      IF(JSTART .NE. 0) GO TO 120
CC
CC     IF THIS IS THE FIRST STEP,PERFORM
CC     OTHER INITIALIZATIONS
C      IWM(LETF) = 0
C      IWM(LCTF) = 0
C      K=1
C      KOLD=0
C      HOLD=0.0E0
C      JSTART=1
C      PSI(1)=H
C      CJOLD = 1.0E0/H
C      CJ = CJOLD
C      S = 100.E0
C      JCALC = -1
C      DELNRM=1.0E0
C      IPHASE = 0
C      NS=0
C120   CONTINUE
CC
CC
CC
CC
CC
CC------------------------------------------------------------------
CC     BLOCK 2
CC     COMPUTE COEFFICIENTS OF FORMULAS FOR
CC     THIS STEP.
CC------------------------------------------------------------------
C200   CONTINUE
C      KP1=K+1
C      KP2=K+2
C      KM1=K-1
C      XOLD=X
C      IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0
C      NS=MIN(NS+1,KOLD+2)
C      NSP1=NS+1
C      IF(KP1 .LT. NS)GO TO 230
CC
C      BETA(1)=1.0E0
C      ALPHA(1)=1.0E0
C      TEMP1=H
C      GAMMA(1)=0.0E0
C      SIGMA(1)=1.0E0
C      DO 210 I=2,KP1
C         TEMP2=PSI(I-1)
C         PSI(I-1)=TEMP1
C         BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
C         TEMP1=TEMP2+H
C         ALPHA(I)=H/TEMP1
C         SIGMA(I)=(I-1)*SIGMA(I-1)*ALPHA(I)
C         GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
C210      CONTINUE
C      PSI(KP1)=TEMP1
C230   CONTINUE
CC
CC     COMPUTE ALPHAS, ALPHA0
C      ALPHAS = 0.0E0
C      ALPHA0 = 0.0E0
C      DO 240 I = 1,K
C        ALPHAS = ALPHAS - 1.0E0/I
C        ALPHA0 = ALPHA0 - ALPHA(I)
C240     CONTINUE
CC
CC     COMPUTE LEADING COEFFICIENT CJ
C      CJLAST = CJ
C      CJ = -ALPHAS/H
CC
CC     COMPUTE VARIABLE STEPSIZE ERROR COEFFICIENT CK
C      CK = ABS(ALPHA(KP1) + ALPHAS - ALPHA0)
C      CK = MAX(CK,ALPHA(KP1))
CC
CC     DECIDE WHETHER NEW JACOBIAN IS NEEDED
C      TEMP1 = (1.0E0 - XRATE)/(1.0E0 + XRATE)
C      TEMP2 = 1.0E0/TEMP1
C      IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
C      IF (CJ .NE. CJLAST) S = 100.E0
CC
CC     CHANGE PHI TO PHI STAR
C      IF(KP1 .LT. NSP1) GO TO 280
C      DO 270 J=NSP1,KP1
C         DO 260 I=1,NEQ
C260         PHI(I,J)=BETA(J)*PHI(I,J)
C270      CONTINUE
C280   CONTINUE
CC
CC     UPDATE TIME
C      X=X+H
CC
CC
CC
CC
CC
CC------------------------------------------------------------------
CC     BLOCK 3
CC     PREDICT THE SOLUTION AND DERIVATIVE,
CC     AND SOLVE THE CORRECTOR EQUATION
CC------------------------------------------------------------------
CC
CC     FIRST,PREDICT THE SOLUTION AND DERIVATIVE
C300   CONTINUE
C      DO 310 I=1,NEQ
C         Y(I)=PHI(I,1)
C310      YPRIME(I)=0.0E0
C      DO 330 J=2,KP1
C         DO 320 I=1,NEQ
C            Y(I)=Y(I)+PHI(I,J)
C320         YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
C330   CONTINUE
C      PNORM = SDANRM (NEQ,Y,WT,RPAR,IPAR)
CC
CC
CC
CC     SOLVE THE CORRECTOR EQUATION USING A
CC     MODIFIED NEWTON SCHEME.
C      CONVGD= .TRUE.
C      M=0
C      IWM(LNRE)=IWM(LNRE)+1
C      IRES = 0
C      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
C      IF (IRES .LT. 0) GO TO 380
CC
CC
CC     IF INDICATED,REEVALUATE THE
CC     ITERATION MATRIX PD = DG/DY + CJ*DG/DYPRIME
CC     (WHERE G(X,Y,YPRIME)=0). SET
CC     JCALC TO 0 AS AN INDICATOR THAT
CC     THIS HAS BEEN DONE.
C      IF(JCALC .NE. -1)GO TO 340
C      IWM(LNJE)=IWM(LNJE)+1
C      JCALC=0
C      CALL SDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,
C     * IER,WT,E,WM,IWM,RES,IRES,UROUND,JAC,RPAR,
C     * IPAR,NTEMP)
C      CJOLD=CJ
C      S = 100.E0
C      IF (IRES .LT. 0) GO TO 380
C      IF(IER .NE. 0)GO TO 380
C      NSF=0
CC
CC
CC     INITIALIZE THE ERROR ACCUMULATION VECTOR E.
C340   CONTINUE
C      DO 345 I=1,NEQ
C345      E(I)=0.0E0
CC
CC
CC     CORRECTOR LOOP.
C350   CONTINUE
CC
CC     MULTIPLY RESIDUAL BY TEMP1 TO ACCELERATE CONVERGENCE
C      TEMP1 = 2.0E0/(1.0E0 + CJ/CJOLD)
C      DO 355 I = 1,NEQ
C355     DELTA(I) = DELTA(I) * TEMP1
CC
CC     COMPUTE A NEW ITERATE (BACK-SUBSTITUTION).
CC     STORE THE CORRECTION IN DELTA.
C      CALL SDASLV(NEQ,DELTA,WM,IWM)
CC
CC     UPDATE Y,E,AND YPRIME
C      DO 360 I=1,NEQ
C         Y(I)=Y(I)-DELTA(I)
C         E(I)=E(I)-DELTA(I)
C360      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
CC
CC     TEST FOR CONVERGENCE OF THE ITERATION
C      DELNRM=SDANRM(NEQ,DELTA,WT,RPAR,IPAR)
C      IF (DELNRM .LE. 100.E0*UROUND*PNORM) GO TO 375
C      IF (M .GT. 0) GO TO 365
C         OLDNRM = DELNRM
C         GO TO 367
C365   RATE = (DELNRM/OLDNRM)**(1.0E0/M)
C      IF (RATE .GT. 0.90E0) GO TO 370
C      S = RATE/(1.0E0 - RATE)
C367   IF (S*DELNRM .LE. 0.33E0) GO TO 375
CC
CC     THE CORRECTOR HAS NOT YET CONVERGED.
CC     UPDATE M AND TEST WHETHER THE
CC     MAXIMUM NUMBER OF ITERATIONS HAVE
CC     BEEN TRIED.
C      M=M+1
C      IF(M.GE.MAXIT)GO TO 370
CC
CC     EVALUATE THE RESIDUAL
CC     AND GO BACK TO DO ANOTHER ITERATION
C      IWM(LNRE)=IWM(LNRE)+1
C      IRES = 0
C      CALL RES(X,Y,YPRIME,DELTA,IRES,
C     *  RPAR,IPAR)
C      IF (IRES .LT. 0) GO TO 380
C      GO TO 350
CC
CC
CC     THE CORRECTOR FAILED TO CONVERGE IN MAXIT
CC     ITERATIONS. IF THE ITERATION MATRIX
CC     IS NOT CURRENT,RE-DO THE STEP WITH
CC     A NEW ITERATION MATRIX.
C370   CONTINUE
C      IF(JCALC.EQ.0)GO TO 380
C      JCALC=-1
C      GO TO 300
CC
CC
CC     THE ITERATION HAS CONVERGED.  IF NONNEGATIVITY OF SOLUTION IS
CC     REQUIRED, SET THE SOLUTION NONNEGATIVE, IF THE PERTURBATION
CC     TO DO IT IS SMALL ENOUGH.  IF THE CHANGE IS TOO LARGE, THEN
CC     CONSIDER THE CORRECTOR ITERATION TO HAVE FAILED.
C375   IF(NONNEG .EQ. 0) GO TO 390
C      DO 377 I = 1,NEQ
C377      DELTA(I) = MIN(Y(I),0.0E0)
C      DELNRM = SDANRM(NEQ,DELTA,WT,RPAR,IPAR)
C      IF(DELNRM .GT. 0.33E0) GO TO 380
C      DO 378 I = 1,NEQ
C378      E(I) = E(I) - DELTA(I)
C      GO TO 390
CC
CC
CC     EXITS FROM BLOCK 3
CC     NO CONVERGENCE WITH CURRENT ITERATION
CC     MATRIX,OR SINGULAR ITERATION MATRIX
C380   CONVGD= .FALSE.
C390   JCALC = 1
C      IF(.NOT.CONVGD)GO TO 600
CC
CC
CC
CC
CC
CC------------------------------------------------------------------
CC     BLOCK 4
CC     ESTIMATE THE ERRORS AT ORDERS K,K-1,K-2
CC     AS IF CONSTANT STEPSIZE WAS USED. ESTIMATE
CC     THE LOCAL ERROR AT ORDER K AND TEST
CC     WHETHER THE CURRENT STEP IS SUCCESSFUL.
CC------------------------------------------------------------------
CC
CC     ESTIMATE ERRORS AT ORDERS K,K-1,K-2
C      ENORM = SDANRM(NEQ,E,WT,RPAR,IPAR)
C      ERK = SIGMA(K+1)*ENORM
C      TERK = (K+1)*ERK
C      EST = ERK
C      KNEW=K
C      IF(K .EQ. 1)GO TO 430
C      DO 405 I = 1,NEQ
C405     DELTA(I) = PHI(I,KP1) + E(I)
C      ERKM1=SIGMA(K)*SDANRM(NEQ,DELTA,WT,RPAR,IPAR)
C      TERKM1 = K*ERKM1
C      IF(K .GT. 2)GO TO 410
C      IF(TERKM1 .LE. 0.5E0*TERK)GO TO 420
C      GO TO 430
C410   CONTINUE
C      DO 415 I = 1,NEQ
C415     DELTA(I) = PHI(I,K) + DELTA(I)
C      ERKM2=SIGMA(K-1)*SDANRM(NEQ,DELTA,WT,RPAR,IPAR)
C      TERKM2 = (K-1)*ERKM2
C      IF(MAX(TERKM1,TERKM2).GT.TERK)GO TO 430
CC     LOWER THE ORDER
C420   CONTINUE
C      KNEW=K-1
C      EST = ERKM1
CC
CC
CC     CALCULATE THE LOCAL ERROR FOR THE CURRENT STEP
CC     TO SEE IF THE STEP WAS SUCCESSFUL
C430   CONTINUE
C      ERR = CK * ENORM
C      IF(ERR .GT. 1.0E0)GO TO 600
CC
CC
CC
CC
CC
CC------------------------------------------------------------------
CC     BLOCK 5
CC     THE STEP IS SUCCESSFUL. DETERMINE
CC     THE BEST ORDER AND STEPSIZE FOR
CC     THE NEXT STEP. UPDATE THE DIFFERENCES
CC     FOR THE NEXT STEP.
CC------------------------------------------------------------------
C      IDID=1
C      IWM(LNST)=IWM(LNST)+1
C      KDIFF=K-KOLD
C      KOLD=K
C      HOLD=H
CC
CC
CC     ESTIMATE THE ERROR AT ORDER K+1 UNLESS:
CC        ALREADY DECIDED TO LOWER ORDER, OR
CC        ALREADY USING MAXIMUM ORDER, OR
CC        STEPSIZE NOT CONSTANT, OR
CC        ORDER RAISED IN PREVIOUS STEP
C      IF(KNEW.EQ.KM1.OR.K.EQ.IWM(LMXORD))IPHASE=1
C      IF(IPHASE .EQ. 0)GO TO 545
C      IF(KNEW.EQ.KM1)GO TO 540
C      IF(K.EQ.IWM(LMXORD)) GO TO 550
C      IF(KP1.GE.NS.OR.KDIFF.EQ.1)GO TO 550
C      DO 510 I=1,NEQ
C510      DELTA(I)=E(I)-PHI(I,KP2)
C      ERKP1 = (1.0E0/(K+2))*SDANRM(NEQ,DELTA,WT,RPAR,IPAR)
C      TERKP1 = (K+2)*ERKP1
C      IF(K.GT.1)GO TO 520
C      IF(TERKP1.GE.0.5E0*TERK)GO TO 550
C      GO TO 530
C520   IF(TERKM1.LE.MIN(TERK,TERKP1))GO TO 540
C      IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD))GO TO 550
CC
CC     RAISE ORDER
C530   K=KP1
C      EST = ERKP1
C      GO TO 550
CC
CC     LOWER ORDER
C540   K=KM1
C      EST = ERKM1
C      GO TO 550
CC
CC     IF IPHASE = 0, INCREASE ORDER BY ONE AND MULTIPLY STEPSIZE BY
CC     FACTOR TWO
C545   K = KP1
C      HNEW = H*2.0E0
C      H = HNEW
C      GO TO 575
CC
CC
CC     DETERMINE THE APPROPRIATE STEPSIZE FOR
CC     THE NEXT STEP.
C550   HNEW=H
C      TEMP2=K+1
C      R=(2.0E0*EST+0.0001E0)**(-1.0E0/TEMP2)
C      IF(R .LT. 2.0E0) GO TO 555
C      HNEW = 2.0E0*H
C      GO TO 560
C555   IF(R .GT. 1.0E0) GO TO 560
C      R = MAX(0.5E0,MIN(0.9E0,R))
C      HNEW = H*R
C560   H=HNEW
CC
CC
CC     UPDATE DIFFERENCES FOR NEXT STEP
C575   CONTINUE
C      IF(KOLD.EQ.IWM(LMXORD))GO TO 585
C      DO 580 I=1,NEQ
C580      PHI(I,KP2)=E(I)
C585   CONTINUE
C      DO 590 I=1,NEQ
C590      PHI(I,KP1)=PHI(I,KP1)+E(I)
C      DO 595 J1=2,KP1
C         J=KP1-J1+1
C         DO 595 I=1,NEQ
C595      PHI(I,J)=PHI(I,J)+PHI(I,J+1)
C      RETURN
CC
CC
CC
CC
CC
CC------------------------------------------------------------------
CC     BLOCK 6
CC     THE STEP IS UNSUCCESSFUL. RESTORE X,PSI,PHI
CC     DETERMINE APPROPRIATE STEPSIZE FOR
CC     CONTINUING THE INTEGRATION, OR EXIT WITH
CC     AN ERROR FLAG IF THERE HAVE BEEN MANY
CC     FAILURES.
CC------------------------------------------------------------------
C600   IPHASE = 1
CC
CC     RESTORE X,PHI,PSI
C      X=XOLD
C      IF(KP1.LT.NSP1)GO TO 630
C      DO 620 J=NSP1,KP1
C         TEMP1=1.0E0/BETA(J)
C         DO 610 I=1,NEQ
C610         PHI(I,J)=TEMP1*PHI(I,J)
C620      CONTINUE
C630   CONTINUE
C      DO 640 I=2,KP1
C640      PSI(I-1)=PSI(I)-H
CC
CC
CC     TEST WHETHER FAILURE IS DUE TO CORRECTOR ITERATION
CC     OR ERROR TEST
C      IF(CONVGD)GO TO 660
C      IWM(LCTF)=IWM(LCTF)+1
CC
CC
CC     THE NEWTON ITERATION FAILED TO CONVERGE WITH
CC     A CURRENT ITERATION MATRIX.  DETERMINE THE CAUSE
CC     OF THE FAILURE AND TAKE APPROPRIATE ACTION.
C      IF(IER.EQ.0)GO TO 650
CC
CC     THE ITERATION MATRIX IS SINGULAR. REDUCE
CC     THE STEPSIZE BY A FACTOR OF 4. IF
CC     THIS HAPPENS THREE TIMES IN A ROW ON
CC     THE SAME STEP, RETURN WITH AN ERROR FLAG
C      NSF=NSF+1
C      R = 0.25E0
C      H=H*R
C      IF (NSF .LT. 3 .AND. ABS(H) .GE. HMIN) GO TO 690
C      IDID=-8
C      GO TO 675
CC
CC
CC     THE NEWTON ITERATION FAILED TO CONVERGE FOR A REASON
CC     OTHER THAN A SINGULAR ITERATION MATRIX.  IF IRES = -2, THEN
CC     RETURN.  OTHERWISE, REDUCE THE STEPSIZE AND TRY AGAIN, UNLESS
CC     TOO MANY FAILURES HAVE OCCURED.
C650   CONTINUE
C      IF (IRES .GT. -2) GO TO 655
C      IDID = -11
C      GO TO 675
C655   NCF = NCF + 1
C      R = 0.25E0
C      H = H*R
C      IF (NCF .LT. 10 .AND. ABS(H) .GE. HMIN) GO TO 690
C      IDID = -7
C      IF (IRES .LT. 0) IDID = -10
C      IF (NEF .GE. 3) IDID = -9
C      GO TO 675
CC
CC
CC     THE NEWTON SCHEME CONVERGED,AND THE CAUSE
CC     OF THE FAILURE WAS THE ERROR ESTIMATE
CC     EXCEEDING THE TOLERANCE.
C660   NEF=NEF+1
C      IWM(LETF)=IWM(LETF)+1
C      IF (NEF .GT. 1) GO TO 665
CC
CC     ON FIRST ERROR TEST FAILURE, KEEP CURRENT ORDER OR LOWER
CC     ORDER BY ONE.  COMPUTE NEW STEPSIZE BASED ON DIFFERENCES
CC     OF THE SOLUTION.
C      K = KNEW
C      TEMP2 = K + 1
C      R = 0.90E0*(2.0E0*EST+0.0001E0)**(-1.0E0/TEMP2)
C      R = MAX(0.25E0,MIN(0.9E0,R))
C      H = H*R
C      IF (ABS(H) .GE. HMIN) GO TO 690
C      IDID = -6
C      GO TO 675
CC
CC     ON SECOND ERROR TEST FAILURE, USE THE CURRENT ORDER OR
CC     DECREASE ORDER BY ONE.  REDUCE THE STEPSIZE BY A FACTOR OF
CC     FOUR.
C665   IF (NEF .GT. 2) GO TO 670
C      K = KNEW
C      H = 0.25E0*H
C      IF (ABS(H) .GE. HMIN) GO TO 690
C      IDID = -6
C      GO TO 675
CC
CC     ON THIRD AND SUBSEQUENT ERROR TEST FAILURES, SET THE ORDER TO
CC     ONE AND REDUCE THE STEPSIZE BY A FACTOR OF FOUR.
C670   K = 1
C      H = 0.25E0*H
C      IF (ABS(H) .GE. HMIN) GO TO 690
C      IDID = -6
C      GO TO 675
CC
CC
CC
CC
CC     FOR ALL CRASHES, RESTORE Y TO ITS LAST VALUE,
CC     INTERPOLATE TO FIND YPRIME AT LAST X, AND RETURN
C675   CONTINUE
C      CALL SDATRP(X,X,Y,YPRIME,NEQ,K,PHI,PSI)
C      RETURN
CC
CC
CC     GO BACK AND TRY THIS STEP AGAIN
C690   GO TO 200
CC
CC------END OF SUBROUTINE SDASTP------
C      END
C      SUBROUTINE SDATRP (X, XOUT, YOUT, YPOUT, NEQ, KOLD, PHI, PSI)
CC***BEGIN PROLOGUE  SDATRP
CC***SUBSIDIARY
CC***PURPOSE  Interpolation routine for SDASSL.
CC***LIBRARY   SLATEC (DASSL)
CC***TYPE      SINGLE PRECISION (SDATRP-S, DDATRP-D)
CC***AUTHOR  PETZOLD, LINDA R., (LLNL)
CC***DESCRIPTION
CC------------------------------------------------------------------
CC     THE METHODS IN SUBROUTINE SDASTP USE POLYNOMIALS
CC     TO APPROXIMATE THE SOLUTION. SDATRP APPROXIMATES THE
CC     SOLUTION AND ITS DERIVATIVE AT TIME XOUT BY EVALUATING
CC     ONE OF THESE POLYNOMIALS,AND ITS DERIVATIVE,THERE.
CC     INFORMATION DEFINING THIS POLYNOMIAL IS PASSED FROM
CC     SDASTP, SO SDATRP CANNOT BE USED ALONE.
CC
CC     THE PARAMETERS ARE:
CC     X     THE CURRENT TIME IN THE INTEGRATION.
CC     XOUT  THE TIME AT WHICH THE SOLUTION IS DESIRED
CC     YOUT  THE INTERPOLATED APPROXIMATION TO Y AT XOUT
CC           (THIS IS OUTPUT)
CC     YPOUT THE INTERPOLATED APPROXIMATION TO YPRIME AT XOUT
CC           (THIS IS OUTPUT)
CC     NEQ   NUMBER OF EQUATIONS
CC     KOLD  ORDER USED ON LAST SUCCESSFUL STEP
CC     PHI   ARRAY OF SCALED DIVIDED DIFFERENCES OF Y
CC     PSI   ARRAY OF PAST STEPSIZE HISTORY
CC------------------------------------------------------------------
CC***ROUTINES CALLED  (NONE)
CC***REVISION HISTORY  (YYMMDD)
CC   830315  DATE WRITTEN
CC   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
CC   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
CC   901026  Added explicit declarations for all variables and minor
CC           cosmetic changes to prologue.  (FNF)
CC***END PROLOGUE  SDATRP
CC
C      INTEGER  NEQ, KOLD
C      REAL  X, XOUT, YOUT(*), YPOUT(*), PHI(NEQ,*), PSI(*)
CC
C      INTEGER  I, J, KOLDP1
C      REAL  C, D, GAMMA, TEMP1
CC
CC***FIRST EXECUTABLE STATEMENT  SDATRP
C      KOLDP1=KOLD+1
C      TEMP1=XOUT-X
C      DO 10 I=1,NEQ
C         YOUT(I)=PHI(I,1)
C10       YPOUT(I)=0.0E0
C      C=1.0E0
C      D=0.0E0
C      GAMMA=TEMP1/PSI(1)
C      DO 30 J=2,KOLDP1
C         D=D*GAMMA+C/PSI(J-1)
C         C=C*GAMMA
C         GAMMA=(TEMP1+PSI(J-1))/PSI(J)
C         DO 20 I=1,NEQ
C            YOUT(I)=YOUT(I)+C*PHI(I,J)
C20          YPOUT(I)=YPOUT(I)+D*PHI(I,J)
C30       CONTINUE
C      RETURN
CC
CC------END OF SUBROUTINE SDATRP------
C      END
C      SUBROUTINE SDAWTS (NEQ, IWT, RTOL, ATOL, Y, WT, RPAR, IPAR)
CC***BEGIN PROLOGUE  SDAWTS
CC***SUBSIDIARY
CC***PURPOSE  Set error weight vector for SDASSL.
CC***LIBRARY   SLATEC (DASSL)
CC***TYPE      SINGLE PRECISION (SDAWTS-S, DDAWTS-D)
CC***AUTHOR  PETZOLD, LINDA R., (LLNL)
CC***DESCRIPTION
CC------------------------------------------------------------------
CC     THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR
CC     WT ACCORDING TO WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
CC     I=1,-,N.
CC     RTOL AND ATOL ARE SCALARS IF IWT = 0,
CC     AND VECTORS IF IWT = 1.
CC------------------------------------------------------------------
CC***ROUTINES CALLED  (NONE)
CC***REVISION HISTORY  (YYMMDD)
CC   830315  DATE WRITTEN
CC   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
CC   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
CC   901026  Added explicit declarations for all variables and minor
CC           cosmetic changes to prologue.  (FNF)
CC***END PROLOGUE  SDAWTS
CC
C      INTEGER  NEQ, IWT, IPAR(*)
C      REAL  RTOL(*), ATOL(*), Y(*), WT(*), RPAR(*)
CC
C      INTEGER  I
C      REAL  ATOLI, RTOLI
CC
CC***FIRST EXECUTABLE STATEMENT  SDAWTS
C      RTOLI=RTOL(1)
C      ATOLI=ATOL(1)
C      DO 20 I=1,NEQ
C         IF (IWT .EQ.0) GO TO 10
C           RTOLI=RTOL(I)
C           ATOLI=ATOL(I)
C10         WT(I)=RTOLI*ABS(Y(I))+ATOLI
C20         CONTINUE
C      RETURN
CC------END OF SUBROUTINE SDAWTS------------------------------------
C      END
C*****END precision > single
