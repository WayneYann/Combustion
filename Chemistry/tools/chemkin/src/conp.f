C  CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:32 $
      SUBROUTINE CONP (LINKCK, LIN, LOUT, LENIWK, LENRWK, LENCWK,  
     1                  KMAX, IWORK, RWORK, CWORK, KSYM, X, Z, 
     2                  RTOL, ATOL)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (ZERO=0.0, NK=5, NLMAX=55, ITOL=1, IOPT=0, ITASK=1)
      DIMENSION IWORK(LENIWK), RWORK(LENRWK), X(KMAX), Z(KMAX)
      CHARACTER*16 CWORK(LENCWK), KSYM(KMAX), PRVERS, PRDATE, PREC
      CHARACTER*80 LINE
      LOGICAL KERR, IERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH, FUN
      COMMON /ICONS/ KK, NP, NWT, NH, NWDOT
C
C     V.1.6 98/02/23 (E. Meeks)
C     1. Fix Bug#0151: Remove lines that cannot be reached in 
C        IF (ISTATE .LE. -2) block after call to VODE. 
C     V.1.5 97/03/01
C     1. implement main "driver" program
C     VERSION 1.3, 96/05/24
C     1. initial sccs version
C     VERSION 1.2:
C     1.  Implement new VODE solver
C     VERSION 1.3:
C     1.  Add IFLAG to CKLEN, CKINIT argument list.
C     VERSION 1.4, 97/01/29:
C     1.  Change character index "(:" to "(1:"
C
      DATA PRVERS/'1.6'/, PRDATE/'98/02/23'/,
C*****precision > double
     1PREC/'DOUBLE'/
C*****END precision > double
C*****precision > single
C     1PREC/'SINGLE'/
C*****END precision > single
C
      WRITE (LOUT, '(/A, /1X,A, A, A, A, /A, /A, //)')
     1' CONP: CHEMKIN-III Constant Pressure Kinetics Code,',
     2  PREC(1:CKLSCH(PREC)), ' PRECISION VERS. ',
     3  PRVERS(1:CKLSCH(PRVERS)+1), PRDATE,
     4  ' Copyright 1995, Sandia Corporation.',
     5' The U.S. Government retains a limited license in this software.'
C
      KERR = .FALSE.
      DO 05 K = 1, KMAX
         X(K) = 0.0
         KSYM(K) = ' '
   05 CONTINUE
C
      CALL CKLEN  (LINKCK, LOUT, LENI, LENR, LENC, IFLAG)
      IF (IFLAG .GT. 0) GO TO 7000
      CALL CKINIT (LENIWK, LENRWK, LENCWK, LINKCK, LOUT, IWORK,
     1             RWORK, CWORK, IFLAG)
      IF (IFLAG .GT. 0) GO TO 7000
      CALL CKINDX (IWORK, RWORK, MM, KK, II, NFIT)
C
      NEQ   = KK + 1
      LRW   = 22 + 9*NEQ + 2*NEQ**2
      NVODE = LENR + 1
      NP    = NVODE + LRW
      NWT   = NP + 1
      NH    = NWT  + KK
      NWDOT = NH   + KK
      NTOT  = NWDOT+ KK - 1
C
      LIW   = 30 + NEQ
      IVODE = LENI + 1
      ITOT  = IVODE + LIW - 1
C
      IF (KK.GT.KMAX .OR. LENRWK.LT.NTOT .OR. LENIWK.LT.ITOT) THEN
         IF (KK .GT. KMAX)  WRITE (LOUT, *)
     1   ' Error...KMAX too small...must be at least ',KK
         IF (LENRWK .LT. NTOT) WRITE (LOUT, *)
     1   ' Error...LENRWK too small...must be at least', NTOT
         IF (LENIWK .LT. ITOT) WRITE (LOUT, *)
     1   ' Error...LENIWK too small...must be at least', ITOT
         GO TO 7000
      ENDIF
C
      CALL CKSYMS (CWORK, LOUT, KSYM, IERR)
      IF (IERR) KERR = .TRUE.
      CALL CKWT   (IWORK, RWORK, RWORK(NWT))
      CALL CKRP   (IWORK, RWORK, RU, RUC, PATM)
C
C     Pressure and temperature
C
      WRITE (LOUT, '(/A,//A)')
     1' ADIABATIC FIXED PRESSURE PROBLEM,',
     2' INPUT PRESSURE(ATM) AND TEMPERATURE(K):'
      READ  (LIN,    *) PA, T
      WRITE (LOUT,7105) PA, T
      RWORK(NP) = PA*PATM
C
   40 CONTINUE
C     Initial non-zero moles
C
      LINE = ' '
      WRITE (LOUT, '(/A)') ' INPUT MOLES OF NEXT SPECIES'
      READ  (LIN,  '(A)', END=45)   LINE
      WRITE (LOUT, '(1X,A)') LINE
      CALL CKDTAB (LINE)
      ILEN = INDEX (LINE, '!')
      IF (ILEN .EQ. 1) GO TO 40
C
      ILEN = ILEN - 1
      IF (ILEN .LE. 0) ILEN = LEN(LINE)
      IF (INDEX(LINE(1:ILEN), 'END') .EQ. 0) THEN
         IF (LINE(1:ILEN) .NE. ' ') THEN
            CALL CKSNUM (LINE(1:ILEN), 1, LOUT, KSYM, KK, KNUM,
     1                   NVAL, VAL, IERR)
            IF (IERR) THEN
               WRITE (LOUT,*) ' Error reading moles...'
               KERR = .TRUE.
            ELSE
               X(KNUM) = VAL
            ENDIF
         ENDIF
         GO TO 40
      ENDIF
C
   45 CONTINUE
C     Final time and print interval
C
      WRITE (LOUT, '(/A)') ' INPUT FINAL TIME AND DT'
      READ  (LIN,    *) T2, DT
      WRITE (LOUT,7105) T2, DT
C
      IF (KERR) GO TO 7000
C
C     Normalize the mole fractions
      CALL CKNORM (X, KK)
C
C     Initial conditions and mass fractions
      TT1  = 0.0
      Z(1) = T
      CALL CKXTY (X, IWORK, RWORK, Z(2))
C
C     Integration control parameters for VODE
C
      TT2   = TT1
      MF = 22
      ISTATE= 1
      NLINES=NLMAX + 1
C
  250 CONTINUE
C     Integration loop
C
      IF (NLINES .GE. NLMAX) THEN
C        Print page heading
C
         WRITE (LOUT, 7003)
         WRITE (LOUT, 7100) (KSYM(K)(1:10), K=1,MIN(NK,KK))
         NLINES = 1
C
         DO 200 K1 = NK+1, KK, NK
            WRITE (LOUT, 7110) (KSYM(K)(1:10), K=K1, MIN(K1+NK-1, KK))
            NLINES = NLINES + 1
  200    CONTINUE
      ENDIF
C
      T = Z(1)
      CALL CKYTX (Z(2), IWORK, RWORK, X)
C
C     Print the solution
      WRITE (LOUT, 7105) TT1, T, (X(K), K=1,MIN(NK,KK))
      NLINES = NLINES + 1
C
      DO 300 K1 = NK+1, KK, NK
         WRITE (LOUT, 7115) (X(K), K=K1, MIN(K1+NK-1,KK))
         NLINES = NLINES + 1
  300 CONTINUE
C
      IF (TT2 .GE. T2) GO TO 7000
      TT2 = MIN(TT2 + DT, T2)
C
  350 CONTINUE
C
C     Call the differential equation solver
C
C*****precision > single
C      CALL SVODE
C*****END precision > single
C*****precision > double
      CALL DVODE
C*****END precision > double
     *           (FUN, NEQ, Z, TT1, TT2, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK(NVODE), LRW, IWORK(IVODE),
     2            LIW, JAC, MF, RWORK, IWORK)
C
      IF (ISTATE .LE. -2) THEN
         WRITE (LOUT,*) ' ISTATE=',ISTATE
         GO TO 7000
      ENDIF
      GO TO 250
C
C     FORMATS
 7003 FORMAT (1H1)
 7100 FORMAT (2X, 'T(SEC)', 6X, 'TMP(K)', 6X, 5(1X,A10))
 7105 FORMAT (12E11.3)
 7110 FORMAT (26X, 5(1X,A10))
 7115 FORMAT (22X, 10E11.3)
C
 7000 CONTINUE
C     end of SUBROUTINE CPRUN
      RETURN
      END
C
      SUBROUTINE FUN (N, TIME, Z, ZP, RPAR, IPAR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER(I-N)
C*****END precision > single
C
      COMMON /ICONS/ KK, NP, NWT, NH, NWDOT
      DIMENSION Z(*), ZP(*), RPAR(*), IPAR(*)
C
C     Variables in Z are:  Z(1)   = T
C                          Z(K+1) = Y(K)
C
C     Call CHEMKIN subroutines
C
      CALL CKRHOY (RPAR(NP), Z(1), Z(2), IPAR, RPAR, RHO)
      CALL CKCPBS (Z(1), Z(2), IPAR, RPAR, CPB)
      CALL CKWYP  (RPAR(NP), Z(1), Z(2), IPAR, RPAR, RPAR(NWDOT))
      CALL CKHMS  (Z(1), IPAR, RPAR, RPAR(NH))
C
C     Form governing equation
C
      SUM = 0.0
      DO 100 K = 1, KK
         H    = RPAR(NH    + K - 1)
         WDOT = RPAR(NWDOT + K - 1)
         WT   = RPAR(NWT   + K - 1)
         ZP(K+1) = WDOT * WT / RHO
         SUM = SUM + H * WDOT * WT
 100  CONTINUE
      ZP(1) = -SUM / (RHO*CPB)
C
C     end of SUBROUTINE FUN
      RETURN
      END
