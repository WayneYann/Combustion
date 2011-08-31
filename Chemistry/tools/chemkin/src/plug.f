C    CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
C
      SUBROUTINE PLUG (LIN, LOUT, LINCK, LINSK, LIPAR, IPAR, LRPAR,
     1                 RPAR, LCPAR, CPAR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /PTR/ NICK, NISK, NKFIRS, NKLAST, NKKPHA, NKDOM,
     1             NKOCC, NID1, NID0, NITOT, NRCK, NRSK, NSDEN,
     2             NSIDOT, NF, NFP, NWDOT, NX, NXMOL, NY, NZ, NZOLD,
     3             NZP, NWT, NHML, NACT, NSDOT, NRD1, NRD0, NRTOT
      COMMON /PARAM/ RU, PATM, KG, KG4, KS, KB, NEQ, IENERG, KGS,
     1               KT, NNPHAS, NFSURF, NLSURF, NFBULK, NLBULK,
     2               LENICK, LENRCK, LENISK, LENRSK, LIDWK0, LRDWK0,
     3               LIDWK1, LRDWK1
      COMMON /PROPS/ P, T
C
      DIMENSION  IPAR(LIPAR), RPAR(LRPAR), INFO(15)
      CHARACTER*16  CPAR(LCPAR)
      LOGICAL KERR
      EXTERNAL RES1
C
C         PRINT HEADER
C
      WRITE (LOUT,373)
  373 FORMAT ('   PLUG:  Plug Flow Reactor Code')
      WRITE (LOUT,3370)
 3370 FORMAT ('          Copyright 1996, Sandia Corporation.')
      WRITE (LOUT,3371)
 3371 FORMAT ('          Under the terms of Contract DE-AC04-94AL85000,
     1there is a')
      WRITE (LOUT,3372)
 3372 FORMAT ('          non-exclusive license for use of this work by o
     1r on behalf of')
      WRITE (LOUT,3373)
 3373 FORMAT ('          the U. S. Government.')
      WRITE (LOUT,374)
  374 FORMAT ('          Version 2.7, 98/04/20 ')
      WRITE (LOUT,375)
C*****precision > double
  375 FORMAT ('          DOUBLE PRECISION')
C*****END precision > double
C*****precision > single
C  375 FORMAT ('          SINGLE PRECISION')
C*****END precision > single
C
C         SET CONTROL PARAMETERS FOR DASSL
C
      DO 236 K=1,15
         INFO(K) = 0
  236 CONTINUE
C
C         SET VARIOUS CONSTANTS AND CHECK FOR SUFFICIENT SPACE IN ARRAYS
C
      KERR = .FALSE.
C     get gas-phase storage requirments
      CALL CKLEN (LINCK, LOUT, LENICK, LENRCK, LENCCK, IFLAG)
      IF (IFLAG .GT. 0) THEN
         WRITE (LOUT, *) 'Error initializing gas-phase linkfile...'
         KERR = .TRUE.
      ENDIF
C
C     get surface mechanism storage requirements
      CALL SKLEN (LINSK, LOUT, LENISK, LENRSK, LENCSK, IFLAG)
      IF (IFLAG .GT. 0) THEN
         WRITE (LOUT, *) 'Error initializing surface linkfile...'
         KERR = .TRUE.
      ENDIF
      IF (LIPAR .LT. LENICK+LENISK) THEN
         WRITE (LOUT, *) 'Error...LIPAR must be at least ',LENICK+LENISK
         KERR = .TRUE.
      ENDIF
      IF (LRPAR .LT. LENRCK+LENRSK) THEN
         WRITE (LOUT, *) 'Error...LRPAR must be at least ',LENRCK+LENRSK
         KERR = .TRUE.
      ENDIF
      IF (KERR) RETURN
C
C     store ICKWRK here
      NICK = 1
C     store RCKWRK here
      NRCK = 1
C
      CALL CKINIT (LENICK, LENRCK, LENCCK, LINCK, LOUT, IPAR(NICK),
     1             RPAR(NRCK), CPAR, IFLAG)
      CALL CKRP   (IPAR(NICK), RPAR(NRCK), RU, RUC, PATM)
C
C     store ISKWRK here
      NISK = LENICK + 1
C     store RSKWRK here
      NRSK = LENRCK + 1
C
      CALL SKINIT (LENISK, LENRSK, LENCSK, LINSK, LOUT, IPAR(NISK),
     1             RPAR(NRSK), CPAR, IFLAG)
      CALL SKINDX (IPAR(NISK), NELEM, KG, KS, KB, KT, NNPHAS, NNSURF, 
     1             NFSURF, NLSURF, NNBULK, NFBULK, NLBULK, IISUR)
C
      NEQ  = KG + KS + 5
      LRDWK1 = 40 + 9*NEQ + NEQ**2
      LIDWK1 = 20 + NEQ
      LRDWK0 = 40 + 9*KS + KS**2
      LIDWK0 = 20 + KS
C
C     store KFIRST of phases here
      NKFIRS = NISK + LENISK 
C     store KLAST of phases here
      NKLAST = NKFIRS + NNPHAS
C     store KKPHAS of phase here
      NKKPHA = NKLAST + NNPHAS
C     space for KDOM for phases
      NKDOM  = NKKPHA + NNPHAS
C     space for KCOV for species
      NKOCC  = NKDOM + NNPHAS
C     RES1 integer workspace
      NID1 = NKOCC + KT
C     RES0 integer workspace
      NID0 = NID1 + LIDWK1
      NITOT  = NID0 + LIDWK0 - 1  
C
C     store SDEN of phases here
      NSDEN  = NRSK + LENRSK
C     store SITDOT of phases here
      NSIDOT = NSDEN + NNPHAS
C     store F of NEQ here
      NF     = NSIDOT + NNPHAS
C     store Fprime of NEQ here
      NFP    = NF     + NEQ
C     space for WDOT for gases
      NWDOT  = NFP    + NEQ
C     space for X for gases
      NX     = NWDOT  + KG
C     space for XMOL for gases
      NXMOL  = NX     + KG
C     space for Y for gases
      NY     = NXMOL  + KG
C     space for Z for site species
      NZ     = NY     + KG
C     space for ZOLD for site species
      NZOLD  = NZ     + KS
C     space for Zprime for site species
      NZP    = NZOLD  + KS
C     space for WT for KTot
      NWT    = NZP    + KS
C     space for HML for KTot
      NHML   = NWT    + KT
C     space for ACT for KTot
      NACT   = NHML   + KT
C     space for SDOT for KTot
      NSDOT  = NACT   + KT
C     RES1 real workspace
      NRD1 = NSDOT  + KT
C     RES0 real workspace
      NRD0 = NRD1 + LRDWK1 
      NRTOT  = NRD0 + LRDWK0 - 1
C
      IF (NITOT .GT. LIPAR) THEN
         WRITE (LOUT, *) 'Error...LIPAR must be at least ',NITOT
         KERR = .TRUE.
      ENDIF
      IF (NRTOT .GT. LRPAR) THEN
         WRITE (LOUT, *) 'Error...LRPAR must be at least ',NRTOT
         KERR = .TRUE.
      ENDIF 
      NCSK = 1
      NKNAME = NCSK + LENCSK
      NCTOT  = NKNAME + KT - 1
      IF (NCTOT .GT. LCPAR) THEN
         WRITE (LOUT, *) 'Error...LCPAR must be at least ',NCTOT
         KERR = .TRUE.
      ENDIF
      IF (KERR) RETURN
C
      CALL SKPKK  (IPAR(NISK), IPAR(NKKPHA), IPAR(NKFIRS),
     1             IPAR(NKLAST))
      CALL SKSYMS (IPAR(NISK), CPAR, LOUT, CPAR(NKNAME), KERR)
      CALL SKWT   (IPAR(NISK), RPAR(NRSK), RPAR(NWT))
      CALL SKSDEN (IPAR(NISK), RPAR(NRSK), RPAR(NSDEN))
      CALL SKCOV  (IPAR(NISK), IPAR(NKOCC))
      KG4  = KG + 4
      KGS  = KG + KS
C
C         READ INPUT FILE
C
      CALL READIN (LIN, LOUT, CPAR(NKNAME), RPAR(NX), RPAR(NY), 
     1             RPAR(NZ), XSTR, XEND, DX, ATOL, RTOL, T, P, U, 
     2             IPAR(NICK), RPAR(NRCK), NNEG, DT, RCHG, 
     3             RPAR(NACT + KGS), IPAR(NKFIRS), IPAR(NKLAST),
     4             IPAR(NKDOM), KERR)
      IF (KERR) RETURN
C
C         COMPUTE GAS DENSITY AND SITE FRACTIONS AT REACTOR INLET
C
      CALL CKRHOY (P, T, RPAR(NY), IPAR(NICK), RPAR(NRCK), RHO)
      INFO(10) = NNEG
      IF (KS .GT. 0) THEN
         CALL PLRES0 (LOUT, TIME, DT, TOUT, P, T, RPAR(NACT),
     1                RPAR(NSDEN), IPAR(NISK), RPAR(NRSK), RPAR(NSDOT),
     2                RPAR(NSIDOT), IPAR(NKFIRS), IPAR(NKLAST), 
     2                RPAR(NX), RPAR(NZ), RPAR(NZP), RPAR(NZOLD),
     3                IPAR(NKDOM), INFO, RTOL, ATOL, RPAR(NRD0),
     3                IPAR(NID0), RPAR, IPAR, RCHG, KERR)
         IF (KERR) RETURN
      ENDIF
  770 CONTINUE
C
C         SET INITIAL CONDITIONS FOR INTEGRATION
C
      CALL PLINIT (INFO, XSTR, TAU, T, RHO, P, U, RPAR(NY), RPAR(NZ),
     1             RPAR(NF), RPAR(NFP))
      X1 = XSTR
      X2 = XSTR
C
  250 CONTINUE
C
C         COMPUTE DEPOSITION RATE
C
      CALL PLDEP (RPAR(NX), RPAR(NZ), RPAR(NACT), P, T, RPAR(NWT),
     1            RPAR(NSDEN), IPAR(NISK), RPAR(NRSK), RPAR(NSDOT),
     2            RPAR(NSIDOT), DEP)
C
C         PRINT CURRENT RESULTS
C
      CALL PLPRNT (LOUT, X1, TAU, T, RHO, P, U, DEP, RPAR(NX),
     1             CPAR(NKNAME), RPAR(NZ), RPAR(NSDOT))
C
C         INTEGRATE ONE STEP USING DASSL
C
      X2 = X2 + DX
      IF (X2 .GT. XEND + 1.E-5) RETURN
C
C*****precision > double
      CALL DDASSL 
C*****END precision > double
C*****precision > single
C      CALL SDASSL 
C*****END precision > single
     1            (RES1, NEQ, X1, RPAR(NF), RPAR(NFP), X2, INFO, RTOL,
     2             ATOL, IDID, RPAR(NRD1), LRDWK1, IPAR(NID1), LIDWK1, 
     3             RPAR, IPAR, JAC)
C
      IF (IDID .GT. 0) THEN
         CALL PLREST (RPAR(NF), IPAR(NICK), RPAR(NRCK), T, RHO, P,
     1                U, TAU, RPAR(NX), RPAR(NZ))
         GO TO 250
      ELSE
         WRITE (LOUT,7000) IDID
         WRITE (LOUT,5000)
         RETURN
      END IF
C
 5000 FORMAT (' EXECUTION TERMINATED DUE TO FATAL ERROR(S)')
 7000 FORMAT (' FAILURE IN DASSL WITH IDID =',I3)
      END
C-----------------------------------------------------------------------
C
      SUBROUTINE READIN (LIN, LOUT, KNAME, X, Y, Z, XSTR, XEND, DX,
     1                   ATOL, RTOL, T, P, U, ICKWK, RCKWK, NNEG, DT,
     2                   RCHG, BACT, KFIRST, KLAST, KDIM, KERR)
C
C         THIS SUBROUTINE READS THE INPUT FILE FOR A PARTICULAR
C         APPLICATION AND CHECKS FOR ERRORS.
C
C         SUMMARY OF AVAILABLE KEYWORDS:
C
C          'ATOL'   ABSOLUTE ERROR TOLERANCE FOR DASSL --
C                   DEFAULT = 1.E-8
C
C          'RTOL'   RELATIVE ERROR TOLERANCE FOR DASSL --
C                   DEFAULT = 1.E-6
C
C          'XSTR'   INLET AXIAL POSITION (CM) -- DEFAULT = 0.
C
C          'XEND'   OUTLET AXIAL POSITION (CM) -- MUST BE INPUT (NO
C                   DEFAULT VALUE)
C
C          'DX'     OUTPUT DISTANCE INTERVAL (CM) -- DEFAULT = 0.01
C
C          'DIAM'   TUBE DIAMETER (CM) -- TO BE INPUT ONLY IF THE
C                   REACTOR IS A ROUND, EMPTY TUBE OF CONSTANT CROSS
C                   SECTION; OTHERWISE SUBROUTINE GEOM MUST BE USED
C
C          'TEMP'   INLET GAS TEMPERATURE (KELVIN) -- WILL BE IGNORED IF
C                   THE AXIAL TEMPERATURE PROFILE IS SPECIFIED (KEYWORD
C                   'TFIX'), BUT MUST BE INPUT OTHERWISE
C
C          'PRES'   INLET PRESSURE (ATM) -- MUST BE INPUT (NO DEFAULT
C                   VALUE)
C
C          'VEL'    INLET VELOCITY (CM/S) -- MUST BE INPUT IF AND ONLY
C                   IF INLET VOLUMETRIC FLOW RATE IS NOT SPECIFIED
C
C          'VDOT'   INLET VOLUMETRIC FLOW RATE (CM**3/S) -- MUST BE
C                   INPUT IF AND ONLY IF INLET VELOCITY IS NOT SPECIFIED
C
C          UNITS FOR INPUT GAS COMPOSITION -- DEFAULT = MOLE FRACTIONS
C
C          'MASS'   MASS FRACTIONS
C          'MOLE'   MOLE FRACTIONS
C
C          'GAS'    INLET GAS-PHASE MOLE OR MASS FRACTION FOR THE GIVEN
C                   SPECIES -- DEFAULT = 0.; INPUT VALUES WILL BE
C                   NORMALIZED, BUT THEIR SUM MUST STILL BE WITHIN 0.01
C                   OF UNITY
C
C          'SURF'   ESTIMATE (INITIAL GUESS) FOR INLET SURFACE SITE
C                   FRACTION OF THE GIVEN SPECIES -- DEFAULT = 0.; INPUT
C          'BULK'   BULK PHASE ACTIVITY FOR THE GIVEN SPECIES (ASSUMED
C                   CONSTANT) -- DEFAULT = 1.; INPUT VALUES FOR EACH
C                   PHASE WILL BE NORMALIZED, BUT THEIR SUM MUST STILL
C                   BE WITHIN 0.01 OF UNITY
C
C          TYPE OF REACTOR -- DEFAULT = ADIABATIC
C
C          'ISO'    ISOTHERMAL REACTOR
C          'ADIA'   ADIABATIC REACTOR
C          'TFIX'   REACTOR WITH SPECIFIED AXIAL TEMPERATURE PROFILE --
C                   FUNCTION TSPEC(X) MUST BE PROVIDED
C          'QFIX'   REACTOR WITH SPECIFIED EXTERNAL HEAT FLUX
C                   PROFILE -- FUNCTION QE(X) MUST BE PROVIDED
C          'HEAT'   REACTOR WITH SPECIFIED AMBIENT TEMPERATURE (TINF)
C                   AND OVERALL HEAT-TRANSFER COEFFICIENT (BIGU)
C
C          'TINF'   AMBIENT TEMPERATURE (K) -- WILL BE USED ONLY IF
C                   KEYWORD 'HEAT' IS SPECIFIED; DEFAULT = 298
C
C          'BIGU'   OVERALL HEAT-TRANSFER COEFFICIENT BASED ON INTERNAL
C                   SURFACE AREA (CGS UNITS) -- MUST BE INPUT IF KEYWORD
C                   'HEAT' IS SPECIFIED
C
C          'VIS'    VISCOSITY OF INLET MIXTURE (POISE) -- DEFAULT = 0.
C
C          'NNEG'   FLAG INSTRUCTING DASSL TO TRY TO CONSTRAIN ALL
C                   COMPONENTS OF THE SOLUTION VECTOR TO BE NON-
C                   NEGATIVE
C
C          'TSTP'   INITIAL TIME STEP FOR INTEGRATION OF FICTITIOUS
C                   TRANSIENT EQUATIONS (SEC) -- DEFAULT = 1.
C
C          'RCHG'   MAXIMUM RELATIVE CHANGE IN Z(K) INDICATING
C                   CONVERGENCE OF FICTITIOUS TRANSIENT EQUATIONS TO
C                   STEADY STATE -- DEFAULT = 1.E-6
C
C          'END'    END OF INPUT FILE -- REQUIRED
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      COMMON /PARAM/ RU, PATM, KG, KG4, KS, KB, NEQ, IENERG, KGS,
     1               KT, NNPHAS, NFSURF, NLSURF, NFBULK, NLBULK,
     2               LENICK, LENRCK, LENISK, LENRSK, LIDWK0, LRDWK0,
     3               LIDWK1, LRDWK1
      COMMON /SIZE/  DIAM, A, DADX, AI, AE, ISIZE
      COMMON /SWAP/  TIN, TINF, BIGU, VISC, IVISC
      DIMENSION Y(KG), X(KG), Z(KS), ICKWK(LENICK), RCKWK(LENRCK), 
     1          VALUE(2), NEC(9),
     1          BACT(*), KFIRST(NNPHAS), KLAST(NNPHAS)
      CHARACTER KNAME(KT)*(*), KEY*4, LINE*80, CKCHUP*4
      EXTERNAL CKCHUP
      LOGICAL KERR, IERR
C
      DATA NEC/9*0/, PI/3.141592654/
C
C         SET DEFAULT VALUES
C
      KERR = .FALSE.
      XSTR   = 0.
      DX     = 0.01
      ISIZE  = 0
      ATOL   = 1.E-8
      RTOL   = 1.E-6
      IENERG = 1
      IMOLF  = 1
      TINF   = 298.
      VISC   = 0.
      IVISC  = 0
      NNEG   = 0
      DT     = 1.
      RCHG   = 1.E-6
      DO 10 K=1,KG
         X(K) = 0.
   10 CONTINUE
      IF (KS .GT. 0) THEN
         DO 11 K=1,KS
            Z(K) = 0.
   11    CONTINUE
      END IF
      IF (KB .GT. 0) THEN
         DO 12 K=1,KB
            BACT(K) = 1.
   12    CONTINUE
      END IF
C
      WRITE (LOUT,9903)
   90 CONTINUE
C
C         READ AND ECHO INPUT LINE
C
      KEY = ' '
      LINE   = ' '
      READ  (LIN,'(A)')  LINE
      WRITE (LOUT,'(10X,A)') LINE(1:70)
      KEY = CKCHUP(LINE,4)
      LINE(1:4) = ' '
C
C         ASSIGN INPUT VALUE TO APPROPRIATE VARIABLE AFTER CHECKING FOR
C         ERRORS
C
      IF (KEY .EQ. 'XSTR') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XSTR, IERR)
C
      ELSE IF (KEY .EQ. 'XEND') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XEND, IERR)
         IF (IERR .EQV. .FALSE.) NEC(1) = 1
C
      ELSE IF (KEY .EQ. 'DX') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DX, IERR)
C
      ELSE IF (KEY .EQ. 'DIAM') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DIAM, IERR)
         IF (IERR .EQV. .FALSE.) THEN
            A     = PI*DIAM**2/4.
            DADX  = 0.
            AI    = PI*DIAM
            AE    = AI
            ISIZE = 1
         END IF
C
      ELSE IF (KEY .EQ. 'ATOL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, ATOL, IERR)
C
      ELSE IF (KEY .EQ. 'RTOL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RTOL, IERR)
C
      ELSE IF (KEY .EQ. 'TEMP') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, T, IERR)
         TIN = T
         IF (IERR .EQV. .FALSE.) NEC(2) = 1
C
      ELSE IF (KEY .EQ. 'PRES') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, P, IERR)
         IF (IERR .EQV. .FALSE.) THEN
            NEC(3) = 1
            P = P * PATM
         END IF
C
      ELSE IF (KEY .EQ. 'VEL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, U, IERR)
         IF (IERR .EQV. .FALSE.) NEC(4) = 1
C
      ELSE IF (KEY .EQ. 'VDOT') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VOL, IERR)
         IF (IERR .EQV. .FALSE.) NEC(7) = 1
C
      ELSE IF (KEY .EQ. 'ISO') THEN
         IENERG = 0
C
      ELSE IF (KEY .EQ. 'ADIA') THEN
         IENERG = 1
C
      ELSE IF (KEY .EQ. 'TFIX') THEN
         IENERG = 2
C
      ELSE IF (KEY .EQ. 'QFIX') THEN
         IENERG = 3
C
      ELSE IF (KEY .EQ. 'HEAT') THEN
         IENERG = 4
C
      ELSE IF (KEY .EQ. 'TINF') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TINF, IERR)
C
      ELSE IF (KEY .EQ. 'BIGU') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, BIGU, IERR)
         IF (IERR .EQV. .FALSE.) NEC(8) = 1
C
      ELSE IF (KEY .EQ. 'VIS') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VISC, IERR)
         IF (IERR .EQV. .FALSE. .AND. VISC.NE.0) IVISC = 1
C
      ELSE IF (KEY .EQ. 'MASS') THEN
         IMOLF = 0
C
      ELSE IF (KEY .EQ. 'MOLE') THEN
         IMOLF = 1
C
      ELSE IF (KEY .EQ. 'GAS') THEN
         CALL CKSNUM (LINE, 1, LOUT, KNAME, KT, KSPEC, NVAL, VALUE,
     1                IERR)
         IF (KSPEC .LE. 0) THEN
            WRITE (LOUT,8202)
            IERR = .TRUE.
         ELSE IF (KSPEC .GT. KG) THEN
            WRITE (LOUT,8303)
            IERR = .TRUE.
         ELSEIF (IERR .EQV. .FALSE.) THEN
            X(KSPEC) = VALUE(1)
         END IF
C
      ELSE IF (KEY .EQ. 'SURF') THEN
         CALL CKSNUM (LINE, 1, LOUT, KNAME, KT, KSPEC, NVAL, VALUE,
     1                IERR)
         IF (KSPEC .LE. 0) THEN
            WRITE (LOUT,8202)
            IERR = .TRUE.
         ELSE IF (KSPEC .LE. KG .OR. KSPEC .GT. KGS) THEN
            WRITE (LOUT,8303)
            IERR = .TRUE.
         ELSEIF (IERR .EQV. .FALSE.) THEN
            Z(KSPEC-KG) = VALUE(1)
         END IF
C
      ELSE IF (KEY .EQ. 'BULK') THEN
         CALL CKSNUM (LINE, 1, LOUT, KNAME, KT, KSPEC, NVAL, VALUE,
     1                IERR)
         IF (KSPEC .LE. 0) THEN
            WRITE (LOUT,8202)
            IERR = .TRUE.
         ELSE IF (KSPEC .LE. KGS .OR. KSPEC .GT. KT) THEN
            WRITE (LOUT,8303)
            IERR = .TRUE.
         ELSEIF (IERR .EQV. .FALSE.) THEN
            BACT(KSPEC-KGS) = VALUE(1)
         END IF
C
      ELSE IF (KEY .EQ. 'NNEG') THEN
         NNEG = 1
C
      ELSE IF (KEY .EQ. 'TSTP') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DT, IERR)
C
      ELSE IF (KEY .EQ. 'RCHG') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RCHG, IERR)
C
      ELSE IF (KEY .EQ. 'END') THEN
         GO TO 6000
C
      ELSE
         WRITE (LOUT,8203)
         IERR = .TRUE.
      END IF
      KERR = KERR.OR.IERR
      GO TO 90
C
 6000 CONTINUE
C
C         CHECK TO SEE IF INPUT IS COMPLETE AND CONSISTENT
C
      SUM = 0.
      DO 5100 K=1,KG
         SUM = SUM + X(K)
 5100 CONTINUE
      IF (ABS(SUM-1.0) .GE. 0.01) THEN
         WRITE (LOUT,8300)
         KERR = .TRUE.
      ELSE
         DO 5103 K=1,KG
            X(K) = X(K)/SUM
 5103    CONTINUE
      END IF
C
      IF (KS .LE. 0) GO TO 771
      DO 5101 N=NFSURF,NLSURF
         SUM = 0.
         DO 5102 K=KFIRST(N),KLAST(N)
            SUM = SUM + Z(K-KG)
 5102    CONTINUE
         IF (ABS(SUM-1.0) .GE. 0.01) THEN
            WRITE (LOUT,8301) N
            KERR = .TRUE.
         ELSE
            DO 5104 K=KFIRST(N),KLAST(N)
               Z(K-KG) = Z(K-KG)/SUM
 5104       CONTINUE
         END IF
 5101 CONTINUE
  771 CONTINUE
C
      IF (KB .LE. 0) GO TO 772
      DO 6101 N=NFBULK,NLBULK
         SUM = 0.
         DO 6102 K=KFIRST(N),KLAST(N)
            SUM = SUM + BACT(K-KGS)
 6102    CONTINUE
         IF (ABS(SUM-1.0) .GE. 0.01) THEN
            WRITE (LOUT,8302) N
            KERR = .TRUE.
         ELSE
            DO 6104 K=KFIRST(N),KLAST(N)
               BACT(K-KGS) = BACT(K-KGS)/SUM
 6104       CONTINUE
         END IF
 6101 CONTINUE
  772 CONTINUE
C
      IF (NEC(1) .NE. 1) THEN
         WRITE (LOUT,8400)
         KERR = .TRUE.
      END IF
      IF (NEC(2) .NE. 1 .AND. IENERG .NE. 2) THEN
         WRITE (LOUT,8450)
         KERR = .TRUE.
      END IF
      IF (NEC(2) .EQ. 1 .AND. IENERG .EQ. 2) THEN
         WRITE (LOUT,8490)
         WRITE (LOUT,8491)
      END IF
      IF (NEC(8) .NE. 1 .AND. IENERG .EQ. 4) THEN
         WRITE (LOUT,9901)
         KERR = .TRUE.
      END IF
      IF (NEC(8) .EQ. 1 .AND. IENERG .NE. 4) WRITE (LOUT,9902)
      IF (NEC(3) .NE. 1) THEN
         WRITE (LOUT,8460)
         KERR = .TRUE.
      END IF
      IF (NEC(4) .NE. 1 .AND. NEC(7) .NE. 1) THEN
         WRITE (LOUT,8470)
         KERR = .TRUE.
      END IF
      IF (NEC(4) .EQ. 1 .AND. NEC(7) .EQ. 1) THEN
         WRITE (LOUT,9123)
         KERR = .TRUE.
      END IF
      IF (ISIZE .NE. 1) THEN
         WRITE (LOUT,8480)
         WRITE (LOUT,8481)
      END IF
      IF (KERR) THEN
         WRITE (LOUT,500)
         RETURN
      END IF
C
C         PUT MASS FRACTIONS INTO Y AND, IF NECESSARY, MOLE FRACTIONS
C         INTO X
C
      IF (IMOLF .NE. 1) THEN
         DO 5350 K=1,KG
            Y(K)=X(K)
 5350    CONTINUE
         CALL CKYTX (Y, ICKWK, RCKWK, X)
      ELSE
         CALL CKXTY (X, ICKWK, RCKWK, Y)
      END IF
C
C         COMPUTE INITIAL TEMPERATURE IF APPROPRIATE
C
      IF (IENERG .EQ. 2) T = TSPEC(XSTR)
C
C         COMPUTE INLET VELOCITY IF NECESSARY
C
      IF (NEC(7) .EQ. 1) THEN
         IF (ISIZE .EQ. 0) THEN
            ZERO = 0.
            CALL GEOM (ZERO, A, DADX, AI, AE, DIAM)
         END IF
         U = VOL/A
      END IF
C
C         FORMAT STATEMENTS
C
  500 FORMAT (' EXECUTION TERMINATED DUE TO FATAL ERROR(S)')
 8200 FORMAT (A4,A76)
 8201 FORMAT (10X,A4,A66)
 8202 FORMAT (' (UNDECLARED SPECIES)')
 8203 FORMAT (' PLUG ERROR: ILLEGAL KEYWORD')
 8300 FORMAT (' PLUG ERROR: GAS-PHASE MOLE OR MASS FRACTIONS ',
     1        ' DO NOT SUM TO UNITY')
 8301 FORMAT (' PLUG ERROR: SITE FRACTIONS DO NOT SUM TO UNITY',
     1        ' FOR PHASE #',I2)
 8302 FORMAT (' PLUG ERROR: BULK ACTIVITIES DO NOT SUM TO UNITY',
     1        ' FOR PHASE #',I2)
 8303 FORMAT (' PLUG ERROR: KEYWORD AND SPECIES NAME ARE INCONSISTENT')
 8400 FORMAT (' PLUG ERROR: FINAL DISTANCE NOT SPECIFIED')
 8450 FORMAT (' PLUG ERROR: INITIAL TEMPERATURE NOT SPECIFIED')
 8460 FORMAT (' PLUG ERROR: INITIAL PRESSURE NOT SPECIFIED')
 8470 FORMAT (' PLUG ERROR: NEITHER INITIAL VELOCITY NOR VOLUMETRIC',
     1        ' FLOW RATE IS SPECIFIED')
 8480 FORMAT (' NOTE: DIAMETER NOT SPECIFIED IN INPUT FILE;')
 8481 FORMAT ('       SUBROUTINE GEOM WILL BE USED')
 8490 FORMAT (' WARNING: INPUT TEMPERATURE WILL BE IGNORED;')
 8491 FORMAT ('          FUNCTION TSPEC WILL BE USED INSTEAD')
 9123 FORMAT (' PLUG ERROR: INLET VELOCITY AND VOLUMETRIC FLOW RATE',
     1        ' CANNOT BOTH BE SPECIFIED')
 9901 FORMAT (' PLUG ERROR: HEAT-TRANSFER COEFFICIENT NOT SPECIFIED')
 9902 FORMAT (' WARNING: HEAT-TRANSFER COEFFICIENT WILL BE IGNORED')
 9903 FORMAT (' ')
C
C     end of SUBROUTINE READIN
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RES0 (TIME, Z, ZP, ZDELTA, IRES, RPAR, IPAR)
C
C         THIS SUBROUTINE PROVIDES THE (FICTITIOUS) TRANSIENT EQUATIONS
C         THAT ARE INTEGRATED BY DASSL TO FIND THE STEADY STATE SURFACE
C         SITE FRACTIONS AT THE REACTOR INLET.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z)
C*****END precision > single
      COMMON /PTR/ NICK, NISK, NKFIRS, NKLAST, NKKPHA, NKDOM,
     1             NKOCC, NID1, NID0, NITOT, NRCK, NRSK, NSDEN,
     2             NSIDOT, NF, NFP, NWDOT, NX, NXMOL, NY, NZ, NZOLD,
     3             NZP, NWT, NHML, NACT, NSDOT, NRD1, NRD0, NRTOT
      COMMON /PARAM/ RU, PATM, KG, KG4, KS, KB, NEQ, IENERG, KGS,
     1               KT, NNPHAS, NFSURF, NLSURF, NFBULK, NLBULK,
     2               LENICK, LENRCK, LENISK, LENRSK, LIDWK0, LRDWK0,
     3               LIDWK1, LRDWK1
      COMMON /PROPS/ P, T
      DIMENSION Z(KS), ZP(KS), ZDELTA(KS), IPAR(*), RPAR(*) 
C
      DO 253 K=1,KG
         RPAR(NACT + K - 1) = RPAR(NX + K - 1)
  253 CONTINUE
      DO 254 K=1,KS
         RPAR(NACT + KG + K - 1) = Z(K)
  254 CONTINUE
      CALL SKRAT (P, T, RPAR(NACT), RPAR(NSDEN), IPAR(NISK),  
     1                  RPAR(NRSK), RPAR(NSDOT), RPAR(NSIDOT))
      DO 600 N=NFSURF,NLSURF
         SUM = 0.
         SDEN = RPAR(NSDEN + N - 1)
         KFIRST = IPAR(NKFIRS + N - 1)
         KLAST  = IPAR(NKLAST + N - 1)
         DO 601 K = KFIRST, KLAST
            SDOT = RPAR(NSDOT + K - 1)
            KOCC = IPAR(NKOCC + K - 1)
            ZDELTA(K - KG) = ZP(K - KG) - SDOT*KOCC/SDEN
            SUM = SUM + Z(K - KG)
  601    CONTINUE
         ZDELTA(KFIRST - KG) = SUM - 1.0
  600 CONTINUE
C
C     end of SUBROUTINE RES0
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RES1 (X, F, FP, DELTA, IRES, RPAR, IPAR)
C
C         THIS SUBROUTINE FORMULATES THE GOVERNING EQUATIONS FOR THE
C         REACTOR IN A FORM SUITABLE FOR USE BY DASSL.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      COMMON /PTR/ NICK, NISK, NKFIRS, NKLAST, NKKPHA, NKDOM,
     1             NKOCC, NID1, NID0, NITOT, NRCK, NRSK, NSDEN,
     2             NSIDOT, NF, NFP, NWDOT, NX, NXMOL, NY, NZ, NZOLD,
     3             NZP, NWT, NHML, NACT, NSDOT, NRD1, NRD0, NRTOT
      COMMON /PARAM/ RU, PATM, KG, KG4, KS, KB, NEQ, IENERG, KGS,
     1               KT, NNPHAS, NFSURF, NLSURF, NFBULK, NLBULK,
     2               LENICK, LENRCK, LENISK, LENRSK, LIDWK0, LRDWK0,
     3               LIDWK1, LRDWK1
      COMMON /SIZE/  DIAM, A, DADX, AI, AE, ISIZE
      COMMON /SWAP/  TIN, TINF, BIGU, VISC, IVISC
      DIMENSION F(NEQ), FP(NEQ), DELTA(NEQ), IPAR(*), RPAR(*)
C
C         TABULATE LOCAL VARIABLES
C
      T   = F(1)
      RHO = F(2)
      P   = F(3)
      U   = F(4)
      US  = U*U
      IF (ISIZE .EQ. 0) CALL GEOM (X, A, DADX, AI, AE, DIAM)
      RHOUA = RHO*U*A
C
C         COMPUTE THERMODYNAMIC PROPERTIES AND SPECIES GENERATION RATES
C
      CALL CKMMWY (F(5), IPAR(NICK), RPAR(NRCK), WTM)
      CALL CKCPBS (T, F(5), IPAR(NICK), RPAR(NRCK), CPB)
      CALL SKHML  (T, IPAR(NISK), RPAR(NRSK), RPAR(NHML))
      CALL CKWYP  (P, T, F(5), IPAR(NICK), RPAR(NRCK), RPAR(NWDOT))
      CALL CKYTX  (F(5), IPAR(NICK), RPAR(NRCK), RPAR(NXMOL))
      DO 255 K=1,KG
         RPAR(NACT + K - 1) = RPAR(NXMOL + K - 1)
  255 CONTINUE
      IF (KS .GT. 0) THEN
         DO 256 K=1,KS
            RPAR(NACT + KG + K - 1) = F(KG4 + K)
  256    CONTINUE
      END IF
      CALL SKRAT  (P, T, RPAR(NACT), RPAR(NSDEN), IPAR(NISK), 
     1                   RPAR(NRSK), RPAR(NSDOT), RPAR(NSIDOT))
C
C         COMPUTE LOCAL DEPOSITION RATE
C
      SUM3 = 0.
      DO 10 K=1,KG
         SDOT = RPAR(NSDOT + K - 1)
         WT   = RPAR(NWT   + K - 1)
         SUM3 = SUM3 + SDOT * WT
   10 CONTINUE
      DEP = AI*SUM3
C
C         COMPONENTS OF SOLUTION VECTOR:
C
C       F(1)     = TEMPERATURE (T)
C       F(2)     = GAS MASS DENSITY (RHO)
C       F(3)     = PRESSURE (P)
C       F(4)     = VELOCITY (U)
C       F(K+4)   = GAS MASS FRACTION Y(K)
C       F(K+KG4) = SURFACE SITE FRACTION Z(K)
C       F(KGS+5) = LOCAL RESIDENCE TIME (TAU)
C
C       FP(K)    = D(F(K))/DX
C
C         ENERGY EQUATION (ADIABATIC REACTOR OR REACTOR WITH SPECIFIED
C                          HEAT FLUX OR HEAT-TRANSFER COEFFICIENT)
C
      IF (IENERG .EQ. 1 .OR. IENERG .EQ. 3 .OR. IENERG .EQ. 4) THEN
         SUM1 = 0.
         SUM2 = 0.
         DO 50 K=1,KG
            HML = RPAR(NHML + K - 1)
            WT = RPAR(NWT + K - 1)
            SUM1 = SUM1 + (HML/WT)*FP(K+4) 
            SUM2 = SUM2 + (HML/WT)*F(K+4) 
   50    CONTINUE
         HEAT = 0.
         IF (IENERG .EQ. 3) THEN
            HEAT = QE(X)*AE
         ELSE IF (IENERG .EQ. 4) THEN
            HEAT = BIGU*(TINF - T)*AI
         END IF
         SUM4 = 0.
         IF (IENERG .EQ. 1 .OR. IENERG .EQ. 3) THEN
            IF (KB .GT. 0) THEN
               DO 20 K=1,KB
                  SDOT = RPAR(NSDOT + KGS + K - 1)
                  WT   = RPAR(NWT   + KGS + K - 1)
                  HML  = RPAR(NHML  + KGS + K - 1)
                  SUM4 = SUM4 + SDOT * HML
   20          CONTINUE
            END IF
         ELSE
            DO 120 K=1,KG
               SDOT = RPAR(NSDOT + K - 1)
               WT   = RPAR(NWT   + K - 1)
               HML  = RPAR(NHML  + K - 1)
               SUM4 = SUM4 - SDOT * HML
  120       CONTINUE
         END IF
         DELTA(1) = RHOUA*(SUM1 + CPB*FP(1) + U*FP(4))
     1              + (SUM2 + US/2.)*DEP + AI*SUM4 - HEAT
C
C         ENERGY EQUATION (ISOTHERMAL REACTOR)
C
      ELSE IF (IENERG .EQ. 0) THEN
         DELTA(1) = FP(1)
C
C         ENERGY EQUATION (REACTOR WITH SPECIFIED TEMPERATURE PROFILE)
C
      ELSE
         DELTA(1) = F(1) - TSPEC(X)
      END IF
C
C         CONTINUITY EQUATION
C
      DELTA(2) = RHO*U*DADX + RHO*A*FP(4) + U*A*FP(2) - DEP
C
C         EQUATION OF STATE
C
      DELTA(3) = P*WTM - RHO*RU*T
C
C         MOMENTUM EQUATION
C
      IF (IVISC .EQ. 0) THEN
         DRAG = 0.
      ELSE
         RENO = DIAM*U*RHO/(VISC*SQRT(T/TIN))
         DRAG = 0.5*AI*RHO*US*FRIC(RENO)
      END IF
      DELTA(4) = A*FP(3) + RHOUA*FP(4) + U*DEP + DRAG
C
C         GAS-PHASE SPECIES EQUATIONS
C
      DO 30 K=1,KG
         WDOT = RPAR(NWDOT + K - 1)
         SDOT = RPAR(NSDOT + K - 1)
         WT   = RPAR(NWT   + K - 1)
         DELTA(K+4) = RHOUA*FP(K+4) + F(K+4)*DEP -
     1                WT * (SDOT*AI + WDOT*A)
   30 CONTINUE
C
C         SURFACE SPECIES EQUATIONS
C
      IF (KS .GT. 0) THEN
         DO 258 N=NFSURF,NLSURF
            SUM = 0.
            DO 40 K=IPAR(NKFIRS+N-1), IPAR(NKLAST+N-1)
               DELTA(K + 4) = RPAR(NSDOT + K - 1)
               SUM = SUM + F(K + 4)
   40       CONTINUE
            DELTA(IPAR(NKDOM + N - 1) + 4) = SUM - 1.0
  258    CONTINUE
      END IF
C
C         RESIDENCE TIME EQUATION
C
      DELTA(NEQ) = FP(NEQ) - 1./F(4)
C
C     end of SUBROUTINE RES1
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GEOM (X, A, DADX, AI, AE, DIAM)
C
C         THIS SUBROUTINE IS USED TO SPECIFY THE REACTOR GEOMETRY;
C         HOWEVER, ITS USE IS NOT REQUIRED IF THE REACTOR IS A ROUND,
C         EMPTY TUBE OF CONSTANT CROSS SECTION (AND NEGLIGIBLE WALL
C         THICKNESS IF KEYWORD 'HEAT' OR 'QFIX' IS SPECIFIED).
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C         CROSS-SECTIONAL (FLOW) AREA (CM**2)
C
      A = 139.
C
C         DERIVATIVE OF CROSS-SECTIONAL AREA WITH RESPECT TO LENGTH (CM)
C
      DADX = 0.
C
C         DEPOSITION AREA PER UNIT LENGTH (CM)
C
      AI = 43.9
      IF (X .LE. 0.1) AI = 2073.9
C
C         EXTERNAL SURFACE AREA PER UNIT LENGTH (CM)
C
      AE = 47.9
C
C         HYDRAULIC DIAMETER (CM)
C
      DIAM = 4.*A/AI
C
C     end of SUBROUTINE GEOM
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION QE(X)
C
C         THIS FUNCTION IS USED TO SPECIFY THE POSITION-DEPENDENT
C         HEAT FLUX (ERG/CM**2/S) TO THE OUTSIDE REACTOR WALL IF
C         THE KEYWORD 'QFIX' APPEARS IN THE INPUT FILE.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      QE = 0.
C
C     end of FUNCTION QE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION TSPEC(X)
C
C         THIS FUNCTION IS USED TO SPECIFY THE AXIAL TEMPERATURE
C         PROFILE IN THE REACTOR IF THE KEYWORD 'TFIX' APPEARS IN
C         THE INPUT FILE.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      TSPEC = 50630./(32.49 - LOG(X + 12.07))
C
C     end of FUNCTION TSPEC
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION FRIC(RENO)
C
C         FANNING FRICTION FACTOR AS A FUNCTION OF REYNOLDS NUMBER
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C         LAMINAR FLOW
C
      IF (RENO .LT. 2100.) THEN
         FRIC = 16./RENO
C
C         TURBULENT FLOW
C
      ELSE
         FRIC = 0.0791/RENO**0.25
C
      END IF
C
C     end of FUNCTION FRIC
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PLABS
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C///////////////////////////////////////////////////////////////////////
C
C                          PLUG FLOW REACTOR CODE
C         WRITTEN BY:
C                          RICHARD S. LARSON
C                          THERMAL AND PLASMA PROCESSES DEPARTMENT 8345
C                          MAIL STOP 9042
C                          SANDIA NATIONAL LABORATORIES
C                          P. O. BOX 969
C                          LIVERMORE, CALIFORNIA  94551-0969
C                          510-294-3008
C
C///////////////////////////////////////////////////////////////////////
C
C
C         THIS PROGRAM SOLVES THE COUPLED STEADY STATE CONTINUITY,
C         MOMENTUM, ENERGY, AND SPECIES BALANCE EQUATIONS FOR A PLUG
C         FLOW REACTOR. BOTH HOMOGENEOUS (GAS-PHASE) AND HETEROGENEOUS
C         (SURFACE) REACTIONS CAN BE ACCOMMODATED. THE REACTOR MAY BE
C         EITHER ISOTHERMAL OR ADIABATIC OR MAY HAVE A SPECIFIED AXIAL
C         TEMPERATURE OR HEAT FLUX PROFILE; ALTERNATIVELY, AN AMBIENT
C         TEMPERATURE AND AN OVERALL HEAT-TRANSFER COEFFICIENT CAN BE
C         SPECIFIED. THE CROSS-SECTIONAL AREA AND SURFACE AREA MAY VARY
C         WITH AXIAL POSITION, AND VISCOUS DRAG IS INCLUDED. IDEAL GAS
C         BEHAVIOR AND SURFACE SITE CONSERVATION ARE ASSUMED.
C
C         VERSION 2.7, 98/04/20 E. Meeks
C         1. Fixed bug#0162: Replaced call to SKHMS with call to 
C            SKHML in subroutine RES1.  Renamed variable HMS to HML
C            to avoid confusion about units. Renamed pointer NHMS
C            to NHML globally.  Added divide by WT(K) in loop 50 in
C            subroutine RES1 (only loops over gas-phase species).
C            Removed multiplication by WT in definition of SUM4 in
C            loops 20 and 120 in RES1.
C         VERSION 2.6, 97/11/18 F. Rupley
C         1. Fixed bug#132a:  Replaced DIMENSION (1) with (*) for
C            BACT, IPAR, and RPAR. 
C         VERSION 2.5, 97/11/11 E. Meeks
C         1. Add dummy subroutine DASERR to allow optional link with 
C            cdassl, as per Action Request #092.
C         2. Remove unused EXTERNAL RES0 in routine PLUG; bug #109.
C         VERSION 2.4, 97/04/21 F. Rupley
C         1. move X1=XSTR time start and X2=XSTR current time outside
C            of PLINIT 
C         VERSION 2.3, 97/03/01 F. Rupley
C         1. make new main "driver" program to set up arrays,
C            to satisfy Reaction Design requirement of providing object
C            files instead of source code - move OPENs to driver;
C            modularize by creating subroutines
C         VERSION 2.2, 96/08/06 F. Rupley
C         1. change LIN=5, LOUT=6 in order to use redirect
C         2. change several WRITE (6 to WRITE (LOUT
C         CHANGES FOR VERSION 2.1 (3-28-96):
C               1. LINKING FILES MAY NOW BE EITHER BINARY OR ASCII
C         CHANGES FOR VERSION 2.0 (1-9-96):
C               1. KEYWORD 'INIT' IS REPLACED BY 'GAS'
C               2. BULK SPECIES ACTIVITIES ARE NOW INPUT RATHER THAN
C                  BEING SET EQUAL TO UNITY
C               3. SPECIES NAMES USED WITH 'GAS', 'SURF', AND 'BULK' ARE
C                  CHECKED TO BE SURE THAT THE TYPE OF PHASE IS CORRECT
C         CHANGES FOR VERSION 1.9 (1-8-96):
C               1. ALLOWANCE IS MADE FOR PROBLEMS WITH NO SURFACE AND/OR
C                  BULK PHASES
C         CHANGES FOR VERSION 1.8 (1-4-96):
C               1. FICTITIOUS TRANSIENT EQUATIONS NOW INVOLVE SPECIES
C                  COVERAGE PARAMETERS (SITE OCCUPANCY NUMBERS)
C         CHANGES FOR VERSION 1.7 (12-21-95):
C               1. MIXTURE VISCOSITY IS NOW INPUT AT INLET TEMPERATURE
C                  RATHER THAN ROOM TEMPERATURE
C               2. OVERALL HEAT-TRANSFER COEFFICIENT IS NOW BASED ON
C                  INTERNAL RATHER THAN EXTERNAL SURFACE AREA
C               3. KEYWORD 'RCHG' ADDED
C         CHANGES FOR VERSION 1.6:
C               1. FICTITIOUS TRANSIENT EQUATIONS NOW INVOLVE SDEN(N)
C               2. EXTENSIVE CHECKING FOR SUFFICIENT SPACE IN ARRAYS
C               3. KEYWORD 'TSTP' ADDED
C         CHANGES FOR VERSION 1.5:
C               1. FANNING FRICTION FACTOR IS NOW COMPUTED FROM
C                  CORRELATIONS RATHER THAN BEING INPUT AS A CONSTANT
C               2. MIXTURE VISCOSITY AT ROOM TEMPERATURE IS NOW INPUT
C         CHANGES FOR VERSION 1.4:
C               1. KEYWORD 'NNEG' ADDED
C
C         THIS PROGRAM IS CURRENT WITH THE FOLLOWING CODE VERSIONS:
C               CHEMKIN-III GAS-PHASE MECHANISM INTERPRETER 5.0
C               CHEMKIN-III GAS-PHASE CHEMICAL KINETICS LIBRARY 5.0
C               CHEMKIN-III SURFACE MECHANISM INTERPRETER 5.0
C               CHEMKIN-III SURFACE KINETICS LIBRARY 6.0
      RETURN
      END
C
      SUBROUTINE PLRES0 (LOUT, TIME, DT, TOUT, P, T, ACT, SDEN, ISKWK,
     1                   RSKWK, SDOT, SITDOT, KFIRST, KLAST, X, Z, ZP,
     2                   ZOLD, KDOM, INFO, RTOL, ATOL, RDWK0, IDWK0,
     3                   RPAR, IPAR, RCHG, KERR)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      COMMON /PARAM/ RU, PATM, KG, KG4, KS, KB, NEQ, IENERG, KGS,
     1               KT, NNPHAS, NFSURF, NLSURF, NFBULK, NLBULK,
     2               LENICK, LENRCK, LENISK, LENRSK, LIDWK0, LRDWK0,
     3               LIDWK1, LRDWK1
      DIMENSION INFO(15), X(KG), Z(KS), ZP(KS), ZOLD(KS), KDOM(NNPHAS),
     1          IDWK0(LIDWK0), RDWK0(LRDWK0), IPAR(*), RPAR(*), 
     2          ACT(KT), ISKWK(LENISK), RSKWK(LENRSK), SDEN(NNPHAS),
     3          SDOT(KT), SITDOT(NNPHAS), KFIRST(NNPHAS), KLAST(NNPHAS)
      LOGICAL KERR
      EXTERNAL RES0
C
      KOUNT = -1
      KKOUNT = -1
      TIME = 0.
      DO 253 K=1,KG
         ACT(K) = X(K)
  253 CONTINUE
      DO 254 K=1,KS
         ACT(K+KG) = Z(K)
  254 CONTINUE
      CALL SKRAT (P, T, ACT, SDEN, ISKWK, RSKWK, SDOT, SITDOT)
      DO 600 K=1,KS
         ZP(K) = SDOT(K+KG)
  600 CONTINUE
      DO 900 N=NFSURF,NLSURF
         ZMAX = -1.
         DO 901 K=KFIRST(N),KLAST(N)
            ZK = Z(K-KG)
            IF (ZK .GT. ZMAX) THEN
               ZMAX = ZK
               KDOM(N) = K
            END IF
  901    CONTINUE
  900 CONTINUE
  206 CONTINUE
      KOUNT = KOUNT + 1
      KKOUNT = KKOUNT + 1
      IF (KOUNT .GT. 200) THEN
         WRITE (LOUT,7001)
         WRITE (LOUT,7004)
         WRITE (LOUT,5000)
         KERR = .TRUE.
         RETURN
      ELSE IF (KKOUNT .GT. 9) THEN
         KKOUNT = 0
         DT = DT*10.
      END IF
      DO 207 K=1,KS
         ZOLD(K) = Z(K)
  207 CONTINUE
      TOUT = TIME + DT
C
C*****precision > double
      CALL DDASSL (RES0, KS, TIME, Z, ZP, TOUT, INFO, RTOL, ATOL,
C*****END precision > double
C*****precision > single
C      CALL SDASSL (RES0, KS, TIME, Z, ZP, TOUT, INFO, RTOL, ATOL,
C*****END precision > single
     1             IDID, RDWK0, LRDWK0, IDWK0, LIDWK0, RPAR, IPAR, JAC)
C
      IF (IDID .LT. 0) THEN
         WRITE (LOUT,7001)
         WRITE (LOUT,7000) IDID
         WRITE (LOUT,5000)
         KERR = .TRUE.
         RETURN
      END IF
      DO 205 K=1,KS
         DEV = ABS(1. - ZOLD(K)/Z(K))
         IF (DEV .GT. RCHG) GO TO 206
  205 CONTINUE
      DO 400 N=NFSURF,NLSURF
         ZMAX = -1.
         DO 401 K=KFIRST(N),KLAST(N)
            ZK = Z(K-KG)
            IF (ZK .GT. ZMAX) THEN
               ZMAX = ZK
               KDOM(N) = K
            END IF
  401    CONTINUE
  400 CONTINUE
      WRITE (LOUT,7002)
C
 5000 FORMAT (' EXECUTION TERMINATED DUE TO FATAL ERROR(S)')
 7000 FORMAT (' FAILURE IN DASSL WITH IDID =',I3)
 7001 FORMAT (/,' SITE FRACTIONS AT REACTOR INLET CANNOT BE FOUND')
 7002 FORMAT (/,' SITE FRACTIONS AT REACTOR INLET FOUND SUCCESSFULLY')
 7004 FORMAT ('   (TOO MANY TIME STEPS TAKEN)')
      RETURN
      END
C
      SUBROUTINE PLDEP (X, Z, ACT, P, T, WT, SDEN, ISKWK, RSKWK,
     1                  SDOT, SITDOT, DEP)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      COMMON /PARAM/ RU, PATM, KG, KG4, KS, KB, NEQ, IENERG, KGS,
     1               KT, NNPHAS, NFSURF, NLSURF, NFBULK, NLBULK,
     2               LENICK, LENRCK, LENISK, LENRSK, LIDWK0, LRDWK0,
     3               LIDWK1, LRDWK1
      DIMENSION X(KG), Z(KS), ACT(KT), WT(KT), SDEN(NNPHAS), 
     1          ISKWK(LENISK), RSKWK(LENRSK), SDOT(KT), 
     2          SITDOT(NNPHAS)
C
      DO 251 K=1,KG
         ACT(K) = X(K)
  251 CONTINUE
      IF (KS .GT. 0) THEN
         DO 252 K=1,KS
            ACT(K+KG) = Z(K)
  252    CONTINUE
      END IF
      CALL SKRAT (P, T, ACT, SDEN, ISKWK, RSKWK, SDOT, SITDOT)
      DEP = 0.
      IF (KB .GT. 0) THEN
         DO 237 K=1,KB
            DEP = DEP + SDOT(K+KGS)*WT(K+KGS)
  237    CONTINUE
      END IF
      RETURN
      END
C
      SUBROUTINE PLPRNT (LOUT, X1, TAU, T, RHO, P, U, DEP, X,
     1                   KNAME, Z, SDOT)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      COMMON /PARAM/ RU, PATM, KG, KG4, KS, KB, NEQ, IENERG, KGS,
     1               KT, NNPHAS, NFSURF, NLSURF, NFBULK, NLBULK,
     2               LENICK, LENRCK, LENISK, LENRSK, LIDWK0, LRDWK0,
     3               LIDWK1, LRDWK1
      DIMENSION X(KG), Z(KS), SDOT(KT)
      CHARACTER*16 KNAME(KT)
C
      WRITE (LOUT,7200) X1
      WRITE (LOUT,7209) TAU
      WRITE (LOUT,7201) T
      WRITE (LOUT,7202) RHO
      PTORR = P*7.5006E-4
      WRITE (LOUT,7203) PTORR
      WRITE (LOUT,7204) U
      IF (KB .GT. 0) WRITE (LOUT,7205) DEP
      WRITE (LOUT,7206)
      DO 300 K1=1,KG,3
         K2 = MIN(K1+2, KG)
         WRITE (LOUT,7115) (KNAME(K),X(K),K=K1,K2)
  300 CONTINUE
      IF (KS .GT. 0) THEN
         WRITE (LOUT,7207)
         DO 301 K1=1,KS,3
            K2 = MIN(K1+2, KS)
            WRITE (LOUT,7115) (KNAME(K+KG),Z(K),K=K1,K2)
  301    CONTINUE
      END IF
      IF (KB .GT. 0) THEN
         WRITE (LOUT,7208)
         DO 302 K1=1,KB,3
            K2 = MIN(K1+2, KB)
            WRITE (LOUT,7115) (KNAME(K+KGS),SDOT(K+KGS),K=K1,K2)
  302    CONTINUE
      END IF
C
 7200 FORMAT (/,' X = ',1PE9.3,' CM')
 7201 FORMAT (' TEMPERATURE = ',F6.1,' K')
 7202 FORMAT (' GAS DENSITY = ',1PE9.3,' G/CM**3')
 7203 FORMAT (' PRESSURE = ',1PE9.3,' TORR')
 7204 FORMAT (' VELOCITY = ',1PE9.3,' CM/S')
 7205 FORMAT (' TOTAL DEPOSITION RATE = ',1PE10.3,' G/CM**2/S')
 7206 FORMAT (' GAS-PHASE MOLE FRACTIONS')
 7207 FORMAT (' SURFACE SITE FRACTIONS')
 7208 FORMAT (' BULK SPECIES DEPOSITION RATES (MOL/CM**2/S)')
 7209 FORMAT (' RESIDENCE TIME = ', 1PE9.3, ' S')
 7115 FORMAT (3(3X,A12,'=',1PE10.3))
      RETURN
      END
C
      SUBROUTINE PLINIT (INFO, XSTR, TAU, T, RHO, P, U, Y, Z, F, FP)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      COMMON /PARAM/ RU, PATM, KG, KG4, KS, KB, NEQ, IENERG, KGS,
     1               KT, NNPHAS, NFSURF, NLSURF, NFBULK, NLBULK,
     2               LENICK, LENRCK, LENISK, LENRSK, LIDWK0, LRDWK0,
     3               LIDWK1, LRDWK1
      DIMENSION INFO(15), F(NEQ), FP(NEQ), Y(KG), Z(KS)
C
      INFO(1) = 0
      INFO(11) = 1
C      X1 = XSTR
C      X2 = XSTR
      TAU = 0.
      F(1) = T
      F(2) = RHO
      F(3) = P
      F(4) = U
      DO 233 K=1,KG
         F(K+4) = Y(K)
  233 CONTINUE
      IF (KS .GT. 0) THEN
         DO 234 K=1,KS
            F(K+KG4) = Z(K)
  234    CONTINUE
      END IF
      F(NEQ) = TAU
      DO 235 K=1,NEQ
         FP(K) = 0.
  235 CONTINUE
      RETURN
      END
C
      SUBROUTINE PLREST (F, ICKWK, RCKWK, T, RHO, P, U, TAU, X, Z)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      COMMON /PARAM/ RU, PATM, KG, KG4, KS, KB, NEQ, IENERG, KGS,
     1               KT, NNPHAS, NFSURF, NLSURF, NFBULK, NLBULK,
     2               LENICK, LENRCK, LENISK, LENRSK, LIDWK0, LRDWK0,
     3               LIDWK1, LRDWK1
      DIMENSION F(NEQ), ICKWK(LENICK), RCKWK(LENRCK), X(KG), Z(KS)
C
      T   = F(1)
      RHO = F(2)
      P   = F(3)
      U   = F(4)
      CALL CKYTX (F(5), ICKWK, RCKWK, X)
      IF (KS .GT. 0) THEN
         DO 238 K=1,KS
            Z(K) = F(K+KG4)
  238    CONTINUE
      END IF
      TAU = F(NEQ)
      RETURN
      END
C 
C the following dummy subroutine allows use of cdassl in place
C of dassl without having an unresolved external.
C
      SUBROUTINE DASERR (T,NEQ,H,E,NFUNCT,Y,YPRIME,WT,DELTA,IPAR,RPAR)
      RETURN
      END
