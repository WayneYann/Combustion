C   CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      SUBROUTINE SKSAMP (LIN, LOUT, LINCK, LINSK, ATOL, RTOL, LIWORK,
     1                   IWORK, LRWORK, RWORK, LCWORK, CWORK)
C
C     V. 3.7 98/03/03 (E. Meeks)
C     1. Action #0114: Remove unused parameter ZERO in SKSAMP.
C     V. 3.6 97/03/20
C     1. checks for NNBULK, NNSURF and add these variables to IPAR
C        (per M. Coltrin)
C     V. 3.5 97/03/01
C     1. make new main "driver" program to set up arrays,
C        to satisfy Reaction Design requirement of providing object
C        files instead of source code.
C     V. 3.4, 96/05/24
C     1. initial sccs version
C     VERSION 3.3 (2/27/95 F. Rupley)
C     1.  Change character index "(:" to "(1:"
C       Gas-phase and surface reaction in a constant volume isothermal
C       container with fixed surface area.
C     VERSION 3.2 (1/19/95 F. Rupley)
C     1.  Add integer error flag to CKLEN,CKINIT,SKLEN,SKINIT
C         call lists.
C     VERSION 3.1 (7/14/94 F. Rupley)
C     1.  Add ISKWRK to some CALL SK lists.
C     VERSION 3.0 (3/15/94 F. Rupley)
C     1.  DOS/PC compatibility effort includes adding file names to
C         OPEN statements, removing unused variables in CALL lists,
C         unusued but possibly initialized variables.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER(I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (NLMAX=50, ITOL=1, IOPT=0, ITASK=1)
      DIMENSION IWORK(LIWORK), RWORK(LRWORK)
      COMMON /RPAR/T, AVRAT, RU
      COMMON /IPAR/KKGAS, KKSURF, KKBULK, KKTOT, NFSURF, NLSURF, NFBULK,
     1             NLBULK, NISK, NIPKK, NIPKF, NIPKL, NRSK, NSDEN,
     2             NRCOV, NX, NWDOT, NWT, NACT, NSDOT, NSITDT, NNSURF,
     3             NNBULK
      CHARACTER CWORK(LCWORK)*16, LINE*80
      LOGICAL KERR, IERR
      EXTERNAL FUN
      DATA KERR/.FALSE./
C
C     Find lengths necessary for arrays
      CALL CKLEN (LINCK, LOUT, LENI, LENR, LENC, IFLAG1)
      CALL SKLEN (LINSK, LOUT, LENIS, LENRS, LENCS, IFLAG2)
      IF (IFLAG1.GT.0 .OR. IFLAG2.GT.0) RETURN
C
      LITOT = LENI + LENIS
      LRTOT = LENR + LENRS
      LCTOT = MAX(LENC, LENCS)
      IF (LIWORK.GE.LITOT .AND. LRWORK.GE.LRTOT .AND. LCWORK.GE.LCTOT)
     1THEN
         CALL CKINIT (LENI, LENR, LENC, LINCK, LOUT, IWORK, RWORK,
     1                CWORK, IFLAG1)
         NISK = LENI + 1
         NRSK = LENR + 1
         CALL SKINIT (LENIS, LENRS, LENCS, LINSK, LOUT, IWORK(NISK),
     1                RWORK(NRSK), CWORK, IFLAG2)
         IF (IFLAG1.GT.0 .OR. IFLAG2.GT.0) RETURN
         CALL CKINDX (IWORK, RWORK, MM, KKGAS, II, NFIT)
         CALL SKINDX (IWORK(NISK), NELEM, KKGAS, KKSUR, KKBULK, 
     1                KKTOT, NPHASE, NNSURF, NFSURF, NLSURF, 
     2                NNBULK, NFBULK, NLBULK, IISUR)
         IKSYM = LCTOT + 1
         IPSYM = IKSYM + KKTOT
         LCTOT = IPSYM + NPHASE - 1
      ENDIF
C
      NIPKK = NISK  + LENIS
      NIPKF = NIPKK + NPHASE
      NIPKL = NIPKF + NPHASE
      NICOV = NIPKL + NPHASE
      NEQ   = KKTOT + 1 + NNSURF
      NIODE = NICOV + KKTOT
      LIW   = 30 + NEQ
      ITOT  = NIODE + LIW - 1
C
      NSDEN = NRSK + LENRS
      NRCOV = NSDEN + NPHASE
      NX    = NRCOV  + KKTOT
      NZ    = NX     + KKTOT
      NWDOT = NZ     + KKTOT + 1 + NNSURF
      NWT   = NWDOT  + KKGAS
      NACT  = NWT    + KKTOT
      NSDOT = NACT   + KKTOT
      NSITDT= NSDOT  + KKTOT
      NRODE = NSITDT + NPHASE
      LRW   = 22 + 9*NEQ + 2*NEQ**2
      NTOT  = NRODE + LRW - 1
C
      IF (LIWORK.LT.ITOT .OR. LRWORK.LT.NTOT .OR. LCWORK.LT.LCTOT)
     1   THEN
         IF (LIWORK .LT. ITOT) WRITE (LOUT, *)
     1   ' ERROR: IWORK must be at least ', ITOT
         IF (LRWORK .LT. NTOT) WRITE (LOUT, *)
     1   ' ERROR: RWORK must be at least ', NTOT
         IF (LCWORK .LT. LCTOT) WRITE (LOUT, *)
     1   ' ERROR: CWORK must be at least ', LCTOT
         RETURN
      ENDIF
C
      CALL SKPKK  (IWORK(NISK), IWORK(NIPKK), IWORK(NIPKF),
     1             IWORK(NIPKL))
      CALL SKSDEN (IWORK(NISK), RWORK(NRSK), RWORK(NSDEN))
      CALL SKCOV  (IWORK(NISK), IWORK(NICOV))
C
      DO 30 K = 1, KKTOT
         RWORK(NRCOV + K - 1) = IWORK(NICOV + K - 1)
         RWORK(NX + K - 1) = 0.0
   30 CONTINUE
C
      CALL SKSYMS (IWORK(NISK), CWORK, LOUT, CWORK(IKSYM), IERR)
      KERR = KERR.OR.IERR
      CALL SKSYMP (IWORK(NISK), CWORK, LOUT, CWORK(IPSYM), IERR)
      KERR = KERR.OR.IERR
      CALL SKWT   (IWORK(NISK), RWORK(NRSK), RWORK(NWT))
      CALL SKRP   (IWORK(NISK), RWORK(NRSK), RU, RUC, PATM)
      IF (KERR) THEN
         WRITE (LOUT, *)
     1   'STOP...ERROR INITIALIZING CONSTANTS...'
         RETURN
      ENDIF
C
C     Pressure and temperature
      WRITE (LOUT, '(/A)')
     1                ' INPUT INITIAL PRESSURE(ATM) AND TEMPERATURE(K)'
      READ  (LIN,    *) PA, T
      WRITE (LOUT,7105) PA, T
      P = PA*PATM
C
C     Initial non-zero moles
   40 CONTINUE
      LINE = ' '
      WRITE (LOUT, '(/A)') ' INPUT INITIAL ACTIVITY OF NEXT SPECIES'
      READ  (LIN,  '(A)', END=45)   LINE
      WRITE (LOUT, '(1X,A)') LINE
      ILEN = INDEX (LINE, '!')
      IF (ILEN .EQ. 1) GO TO 40
C
      ILEN = ILEN - 1
      IF (ILEN .LE. 0) ILEN = LEN(LINE)
      IF (INDEX(LINE(1:ILEN), 'END') .EQ. 0) THEN
         IF (LINE(1:ILEN) .NE. ' ') THEN
            CALL SKSNUM (LINE(1:ILEN), 1, LOUT, CWORK(IKSYM), 
     1                   KKTOT, CWORK(IPSYM), NPHASE, IWORK(NIPKK),
     1                   KNUM, NKF, NVAL, VAL, IERR)
            IF (IERR) THEN
               WRITE (LOUT,*) ' Error reading moles...'
               KERR = .TRUE.
            ELSE
               RWORK(NX + KNUM - 1) = VAL
            ENDIF
         ENDIF
         GO TO 40
      ENDIF
C
   45 CONTINUE
      IF (KERR) THEN
         WRITE (LOUT, *) 'STOP...ERROR INITIALIZING USER INPUT...'
         RETURN
      ENDIF
C
C     Surface area to volume ratio
C
      WRITE (LOUT, '(/A)') ' INPUT SURFACE AREA TO VOLUME RATIO'
      READ  (LIN,    *) AVRAT
      WRITE (LOUT,7105) AVRAT
C
C     Final time and print interval
C
      WRITE (LOUT, '(/A)') ' INPUT FINAL TIME AND DT'
      READ  (LIN,    *) T2, DT
      WRITE (LOUT,7105) T2, DT
C
C     Normalize the mole fractions for each phase
C
      DO 60 N = 1, NPHASE
         XTOT = 0.0
         KFIRST = IWORK(NIPKF + N - 1)
         KLAST  = IWORK(NIPKL + N - 1)
         DO 50 K = KFIRST, KLAST
            XTOT = XTOT + RWORK(NX + K - 1)
   50    CONTINUE
         IF (XTOT .NE. 0.0) THEN
            DO 55 K = KFIRST, KLAST
               RWORK(NX + K - 1) = RWORK(NX + K - 1) / XTOT
   55       CONTINUE
         ELSE
            WRITE (LOUT, *)
     1      ' ERROR...NO SPECIES WERE INPUT FOR PHASE ',
     2      CWORK(IPSYM+N-1)
            KERR = .TRUE.
         ENDIF
   60 CONTINUE
      IF (KERR) THEN
         WRITE (LOUT, *) 'STOP...ERROR INITIALIZING SOLUTION...'
         RETURN
      ENDIF
C
C     Initial conditions
C
      TT1  = 0.0
C     Initial gas-phase mass fractions
      CALL CKXTY (RWORK(NX), IWORK, RWORK, RWORK(NZ))
      IF (NNSURF .GT. 0) THEN
C        Initial surface site fractions
         KFIRST = IWORK(NIPKF + NFSURF - 1)
         KLAST  = IWORK(NIPKL + NLSURF - 1)
         DO 110 K = KFIRST, KLAST
            RWORK(NZ+K-1) = RWORK(NX + K - 1)
  110    CONTINUE
      ENDIF
      IF (NNBULK .GT. 0) THEN
C        Initial bulk deposit amounts
         KFIRST = IWORK(NIPKF + NFBULK - 1)
         KLAST  = IWORK(NIPKL + NLBULK - 1)
         DO 120 K = KFIRST, KLAST
            RWORK(NZ+K-1) = 0.0
  120    CONTINUE
      ENDIF
C     Initial gas-phase mass density
      CALL CKRHOY (P, T, RWORK(NZ), IWORK, RWORK, RWORK(NZ+KKTOT))
      IF (NNSURF .GT. 0) THEN
C        Initial surface site densities
         DO 130 N = NFSURF, NLSURF
            IZ = NZ + KKTOT + 1 + N - NFSURF
            RWORK(IZ) = RWORK(NSDEN + N - 1)
C           Z(KKTOT+1+N-NFSURF+1) = RWORK(NSDEN + N - 1)
  130    CONTINUE
      ENDIF
C
C     Integration control parameters for LSODE
C
      TT2   = TT1
      MF    = 22
      ISTATE= 1
C
C     Integration loop
C
  250 CONTINUE
C
C           Print the solution
C
      CALL CKPY  (RWORK(NZ+KKTOT), T, RWORK(NZ), IWORK, RWORK, P)
      WRITE (LOUT,*) ' '
      WRITE (LOUT,*) ' TIME = ', TT2
      WRITE (LOUT, 7100) P, T, RWORK(NZ+KKTOT)
      WRITE (LOUT, *) ' GAS-PHASE MOLE FRACTIONS'
      CALL CKYTX (RWORK(NZ), IWORK, RWORK, RWORK(NX))
      CALL PRT1  (KKGAS, CWORK(IKSYM), LOUT, RWORK(NX))
C
      IF (NNSURF .GT. 0) THEN
         DO 190 N = NFSURF, NLSURF
            WRITE (LOUT, *) 
     1      ' SURFACE SITE FRACTIONS ON PHASE (SITE) ', N
            KKPHAS = IWORK(NIPKK + N - 1)
            KFIRST = IWORK(NIPKF + N - 1)
            CALL PRT1 (KKPHAS, CWORK(IKSYM+KFIRST-1), LOUT, 
     1                 RWORK(NZ+KFIRST-1))
C
            SUM = 0.0
            KFIRST = IWORK(NIPKF + N - 1)
            KLAST  = IWORK(NIPKL + N - 1)
            DO 185 K = KFIRST, KLAST
               SUM = SUM + RWORK(NZ+K-1)
  185       CONTINUE
            WRITE (LOUT,*)'  SUM OF SURFACE SITE FRACTIONS', SUM
            IZ = NZ + KKTOT + 1 + N - NFSURF
            WRITE (LOUT,*)'  SURFACE SITE DENSITY  ', RWORK(IZ)
  190    CONTINUE
      ENDIF
C
      IF (NNBULK .GT. 0) THEN
         DO 195 N = NFBULK, NLBULK
            WRITE (LOUT, *) ' BULK DEPOSITION (GM/CM**2) IN PHASE ', N
            KKPHAS = IWORK(NIPKK + N - 1)
            KFIRST = IWORK(NIPKF + N - 1)
            CALL PRT1 (KKPHAS, CWORK(IKSYM+KFIRST-1), LOUT, 
     1                 RWORK(NZ+KFIRST-1))
  195    CONTINUE
      ENDIF
C
      IF (TT2 .GE. T2) THEN
         WRITE (LOUT, *) 'STOP...TIME LIMIT REACHED...'
         RETURN
      ENDIF
      TT2 = MIN(TT2 + DT, T2)
C
C     Call the differential equation solver
C
  350 CONTINUE
C*****precision > single
C      CALL SVODE
C*****END precision > single
C*****precision > double
      CALL DVODE
C*****END precision > double
     *           (FUN, NEQ, RWORK(NZ), TT1, TT2, ITOL, RTOL, ATOL, 
     1            ITASK, ISTATE, IOPT, RWORK(NRODE), LRW,
     2            IWORK(NIODE), LIW, JAC, MF, RWORK, IWORK)
C
      IF (ISTATE .LE. -1) THEN
         IF (ISTATE .EQ. -1) THEN
            ISTATE = 2
            GO TO 350
         ELSE
            WRITE (LOUT,*) 'ERROR, ISTATE=',ISTATE
            RETURN
         ENDIF
      ENDIF
      GO TO 250
C
 7003 FORMAT (1H1)
 7100 FORMAT (1H , ' GAS-PHASE STATE', /,
     1   '  P = ', 1PE12.4, ' T = ', 1PE12.4, ' DENSITY = ', 1PE12.4)
 7105 FORMAT (12E11.3)
 7110 FORMAT (26X, 5(1X,A10))
 7115 FORMAT (22X, 10E11.3)
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE FUN (NEQ, TIME, Z, ZP, RWORK, IWORK)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER(I-N)
C*****END precision > single
C
      DIMENSION Z(NEQ), ZP(NEQ), RWORK(*), IWORK(*)
C
      COMMON /RPAR/T, AVRAT, RU
      COMMON /IPAR/KKGAS, KKSURF, KKBULK, KKTOT, NFSURF, NLSURF, NFBULK,
     1             NLBULK, NISK, NIPKK, NIPKF, NIPKL, NRSK, NSDEN, 
     2             NRCOV, NX, NWDOT, NWT, NACT, NSDOT, NSITDT, NNSURF, 
     3             NNBULK
C
C
C     Variables in Z are:  Z(K) = Y(K), K=1,KKGAS
C                          Z(K) = SURFACE SITE FRACTIONS,
C                                 K=KFIRST(NFSURF), KLAST(NLSURF)
C                          Z(K) = BULK SPECIES MASS,
C                                 K=KFIRST(NFBULK), KLAST(NLBULK)
C                          Z(K) = GAS-PHASE MASS DENSITY, K=KKTOT+1
C                          Z(K) = SURFACE SITE MOLAR DENSITIES,
C                                 K=KKTOT+2, KKTOT+1+NNSURF
C
C        Call CHEMKIN and SURFACE CHEMKIN subroutines
C
      CALL CKPY  (Z(KKTOT+1), T, Z(1), IWORK, RWORK, P)
      CALL CKWYP (P, T, Z(1), IWORK, RWORK, RWORK(NWDOT))
      CALL CKYTX (Z, IWORK, RWORK, RWORK(NACT))
C
      IF (NNSURF .GT. 0) THEN
         KFIRST = IWORK(NIPKF + NFSURF - 1)
         KLAST  = IWORK(NIPKL + NLSURF - 1)
         DO 100 K = KFIRST, KLAST
            RWORK(NACT + K - 1) = Z(K)
  100    CONTINUE
      ENDIF
C
      IF (NNBULK .GT. 0) THEN
         KFIRST = IWORK(NIPKF + NFBULK - 1)
         KLAST  = IWORK(NIPKL + NLBULK - 1)
         DO 150 K = KFIRST, KLAST
            RWORK(NACT + K - 1) = RWORK(NX + K - 1)
150      CONTINUE
      ENDIF
C
      IF (NNSURF .GT. 0) THEN
         DO 175 N = NFSURF, NLSURF
            RWORK(NSDEN + N - 1) = Z(KKTOT+1+N-NFSURF+1)
  175    CONTINUE
      ENDIF
C
      CALL SKRAT (P, T, RWORK(NACT), RWORK(NSDEN), IWORK(NISK),
     1            RWORK(NRSK), RWORK(NSDOT), RWORK(NSITDT))
C
C        Form mass density equation
C
      SUM = 0.0
      DO 200 K = 1, KKGAS
         SUM = SUM + AVRAT * RWORK(NSDOT+K-1) * RWORK(NWT+K-1)
 200  CONTINUE
      ZP(KKTOT+1) = SUM
C
C        Form the gas-phase mass conservation equation
C
      DO 300 K = 1, KKGAS
         WDOT = RWORK(NWDOT + K - 1)
         WT   = RWORK(NWT   + K - 1)
         SDOT = RWORK(NSDOT + K - 1)
         ZP(K) = ( - Z(K) * ZP(KKTOT+1) + WDOT * WT
     1             + AVRAT * SDOT * WT ) / Z(KKTOT+1)
 300  CONTINUE
C
      IF (NNSURF .GT. 0) THEN
C        Form the surface mass equations
         DO 410 N = NFSURF, NLSURF
            SITDOT = RWORK(NSITDT + N - 1)
            SDEN0  = RWORK(NSDEN + N - 1)
            KFIRST = IWORK(NIPKF + N - 1)
            KLAST  = IWORK(NIPKL + N - 1)
            DO 400 K = KFIRST, KLAST
             SDOT = RWORK(NSDOT + K - 1)
             RCOV = RWORK(NRCOV + K - 1)
             ZP(K) = (SDOT*RCOV - Z(K) * SITDOT) / SDEN0
  400       CONTINUE
  410    CONTINUE
      ENDIF
C
      IF (NNBULK .GT. 0) THEN
C        Form the bulk mass equations
         KFIRST = IWORK(NIPKF + NFBULK - 1)
         KLAST  = IWORK(NIPKL + NLBULK - 1)
         DO 500 K = KFIRST, KLAST
            ZP(K) = RWORK(NSDOT + K - 1)*RWORK(NWT + K - 1) * AVRAT
  500    CONTINUE
      ENDIF
C
      IF (NNSURF .GT. 0) THEN
C        Form the surface site number-density equations
         DO 575 N = NFSURF, NLSURF
            ZP(KKTOT+1+N-NFSURF+1) = RWORK(NSITDT + N - 1)
  575    CONTINUE
      ENDIF
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE PRT1 (KK, KSYM, LOUT, X)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION X(KK)
      CHARACTER*(*) KSYM(KK)
C
      DO 10 K = 1, KK, 3
         WRITE (LOUT, 6010) (KSYM(L), X(L), L=K, MIN(K+2, KK))
   10 CONTINUE
 6010 FORMAT (3X, 3(A12,'=', 1PE10.3, 4X))
C
      RETURN
      END
