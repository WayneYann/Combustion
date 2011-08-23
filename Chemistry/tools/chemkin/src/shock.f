C  CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      SUBROUTINE SHRUN (LINKCK, LIN, LOUT, LSAVE, LIPAR, LRPAR,
     1                  LCPAR, IPAR, RPAR, CPAR)
C
C     V.2.6 98/07/20
C     1. Action#0182: Corrected equation of state for reflected shock
C        conditions, defining RHO3, in routine SHKIC.
C     2. Action#0183: Corrected Equation for T3OT1 in routine SHKIC;
C        Had multiplication by (xi-eta) instead of divide by (xi-eta).
C     3. Action#0184: Corrected (again) equation for BETA in routine 
C        SHKIC.  Action request #013 was implemented incorrectly.
C        Should be:  BETA=(2.0*GAMMA1*XM12 - (GAMMA1-1)+1.0)/GAMP1
C     V.2.5 98/03/03
C     1. Action#113: Remove unused parameter, ZERO, in SHRUN. (E.Meeks)
C     V.2.4 97/08/27
C     1. Bug#013:  In SHKIC correct first equation for BETA to be:
C        BETA = (2.0*GAMMA1*XM12 - GAMP1 + 1.0)/GAMP1 (was missing +1.0)
C     V.2.3 97/03/28
C     1. initialize I1,I2,I5 to zero in SUBROUTINE READIN
C     V. 2.2 97/03/01
C     1. make new main "driver" programto set up arrays,
C        to satisfy Reaction Design requirement of providing object
C        files instead of source code.
C     V. 2.1, 96/05/24
C     1. initial sccs version
C     VERSION 2.0 (4/17/96 F. Rupley)
C     1. initial CHEMKIN_III version
C     2. set up for ascii linkfile
C     CHANGES FOR VERSION 1.9 (2/27/95 F. Rupley)
C     1.  Change character index "(:" to "(1:"
C     CHANGES FOR VERSION 1.8 (1/19/95 F. Rupley)
C     1.  Add integer error flag to CKLEN,CKINIT call lists.
C     CHANGES FOR VERSION 1.7 (3/15/94 F. Rupley)
C     1.  DOS/PC compatibility effort includes adding file names to
C         OPEN statements, removing unused variables in CALL lists,
C         unusued but possibly initialized variables.
C     CHANGES FOR VERSION 1.6 (6/17/93 F. Rupley)
C     1.  Add binary file data to solution file to facilitate
C         post-processing
C     CHANGES FOR VERSION 1.5
C     1.  Implement new VODE solver
C     CHANGES FOR VERSION 1.4
C     1.  FORMAT 290 in RUNIT changed to reduce length of line below
C         132 characters.
C     CHANGES FOR VERSION 1.3
C     1.  READ (LINKCK) VERS, PREC, KERR
C     CHANGES SINCE VERSION 1.1:
C     1.  modified OPEN statements for unix
C     CHANGES SINCE VERSION 1.0:
C     1.  binary file now contains array lengths
C     2.  changed REAL*8 to DOUBLE PRECISION
C     3.  allow upper/lower case input
C     4.  add cray ctss option
C     5.  error handling in character string routines
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /RCONS/ HO, XLM, DIA, VISCOF, BETA, RU, PATM, VIS,
     1               P1, T1, RHO1, H1, V1, XM1, P2, T2, RHO2, H2, V2,
     2               XM2, P3, T3, RHO3, H3, V3, XM3, VRS, ASHOCK,
     3               ATOL, RTOL, TINIT, DT, TEND, TSTART
      COMMON /ICONS/ IPRB, IGOT, IMOLF, KK, NZZ, NT, NRHO, NKK, NTL,
     1               NV, NA, NCWK, NXX, NYY, NWT, NHMS, NGW, LRW,
     2               ICWK, IKS, ICC, IVODE, NVODE, LIW
C
C     all storage space needed is in RPAR, IPAR and CPAR
      DIMENSION RPAR(LRPAR), IPAR(LIPAR)
      CHARACTER CPAR(LCPAR)*16, PRVERS*16, PRDATE*16, PREC*16, ITITL*76
      LOGICAL LGSCAL, KERR, LCNTUE
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C*****licensing > FLEXlm
CC      CHARACTER FEATUR*16
C      INTEGER FEATUR
C      CHARACTER ERRSTR*512
C*****END licensing > FLEXlm
C
      DATA PRVERS/'2.5'/, PRDATE/'98/03/03'/
C
C     initlize all SHOCK variables to zero
      DATA P1A, P2A, P3A/ 3*0.0/
C
      KERR = .FALSE.
C*****licensing > FLEXlm
CC      FEATUR = 'SHOCK'
C      FEATUR = 7
C      ERRSTR = ''
C      IHFLEX = 0
C      IFLXFL = 0
C      CALL CHKOUT ( FEATUR, IHFLEX, IFLXFL, ERRSTR )
C      IF ( IFLXFL .NE. 0 ) THEN
CC        Report the licensing error and exit
C         WRITE (LOUT,'(/A,I2,/A,//A)')
C     &         'ERROR: could not get a license for Feature #', FEATUR,
C     &         'Flexlm error report: ', ERRSTR
C         RETURN
C      ENDIF
C*****END licensing > FLEXlm
      T1  = 0.0
      T2  = 0.0
      T3  = 0.0
      RHO1 = 0.0
      RHO2 = 0.0
      RHO3 = 0.0
      V1  = 0.0
      V2  = 0.0
      V3  = 0.0
      XM1 = 0.0
      XM2 = 0.0
      XM3 = 0.0
      VRS = 0.0
C
C*****precision > double
      PREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C      PREC = 'SINGLE'
C*****END precision > single
C
      WRITE (LOUT, '(/A, /1X,A, A, A, A, /A, /A, //)')
     1   ' SHOCK: CHEMKIN-III Shock Tube Code,',
     2   PREC(1:CKLSCH(PREC)), ' PRECISION Vers. ',
     3   PRVERS(1:CKLSCH(PRVERS)+1), PRDATE,
     4   ' Copyright 1995, Sandia Corporation.',
     5' The U.S. Government retains a limited license in this software.'
C
C     Read the CHEMKIN linkfile to determine mechanism size
C
C     Set up pointers for the work arrays
      CALL SETWRK (LINKCK, LSAVE, LOUT, LIPAR, LRPAR, LCPAR,
     1             IPAR, RPAR, CPAR, NEQ, KERR)
C*****licensing > FLEXlm
C      IF (KERR) CALL CHKIN ( IHFLEX )
C*****END licensing > FLEXlm
      IF (KERR) RETURN
  110 CONTINUE
C
C     Read keyword input
      CALL READIN (LIN, LOUT, LGSCAL, CPAR(IKS), IPAR(ICWK),
     1             RPAR(NCWK), RPAR(NXX), ITITL, KERR, LCNTUE)
C*****licensing > FLEXlm
C      IF (KERR) CALL CHKIN ( IHFLEX )
C*****END licensing > FLEXlm
      IF (KERR) RETURN
C
C     Write problem information to file LSAVE for later processing,
C     for example, plotting
      WRITE (LSAVE) 'SOLUTION        '
      WRITE (LSAVE) ITITL, IPRB, IMOLF
C
C     Print heading indicating problem type;
C     IPRB is a flag indicating the type of problem to be solved,
C     IPRB = 1, incident shock problem with boundary layer correction
C     IPRB = 2, incident shock problem without boundary layer correction
C     IPRB = 3, reflected shock problem
C
      IF (IPRB .EQ. 1) THEN
         WRITE (LOUT, '(/1X,A)')
     1   'INCIDENT SHOCK WITH BOUNDARY LAYER PROBLEM CORRECTION'
      ELSEIF (IPRB .EQ. 2) THEN
         WRITE (LOUT, '(/1X,A)')
     1   'INCIDENT SHOCK WITHOUT BOUNDARY LAYER PROBLEM CORRECTION'
      ELSEIF (IPRB .EQ. 3) THEN
         WRITE (LOUT, '(/1X,A)')
     1   'REFLECTED SHOCK PROBLEM'
      ENDIF
C
C     Print integration control and initial mole fractions
      WRITE (LOUT, '(/1X,A,4(/5X,A,1PE12.4))')
     1'INTEGRATION PARAMETERS:',
     2'RTOL            ', RTOL,
     3'ATOL            ', ATOL,
     4'T1              ', TINIT,
     5'T2              ', TEND
C
      IF (LGSCAL) THEN
         WRITE (LOUT,'(5X,A,1PE12.4)')
     6   'LGDT            ',DT
      ELSE
         WRITE (LOUT,'(5X,A,1PE12.4)')
     6   'DT              ',  DT
      ENDIF
C
      WRITE (LOUT, '(/1X,A,/)') 'MOLE FRACTIONS'
      WRITE (LOUT, '(5X,A,1PE12.4)')
     1      (CPAR(IKS+K-1), RPAR(NXX+K-1), K = 1, KK)
C
C     Set up initial conditions; convert to mass fractions
      CALL SHKIC (LOUT, IPAR, RPAR, KERR)
C*****licensing > FLEXlm
C      IF (KERR) CALL CHKIN ( IHFLEX )
C*****END licensing > FLEXlm
      IF (KERR) RETURN
      P1A = P1/PATM
      P2A = P2/PATM
C
      IF (IPRB .EQ. 3) THEN
         P3A = P3/PATM
C
         WRITE (LOUT, '(/2(/T24,A,T44,A,T64,A),/)')
     1   'CONDITIONS BEFORE','CONDITIONS BEHIND','CONDITIONS BEHIND',
     2   'INCIDENT SHOCK'   ,'INCIDENT SHOCK',   'REFLECTED SHOCK'
C
         WRITE (LOUT,'(5(/1X, A, T24,1PE10.4,T44,1PE10.4,T64,1PE10.4))')
     1   'PRESSURE    (atm)  ', P1A,  P2A,  P3A,
     2   'TEMPERATURE (K)    ', T1,   T2,   T3,
     3   'DENSITY     (g/cc) ', RHO1, RHO2, RHO3,
     4   'VELOCITY    (cm/s) ', V1,   V2,   V3,
     5   'MACH No.           ', XM1,  XM2,  XM3
C
         WRITE (LOUT, '(/1X,A,1PE10.4)')
     1  'REFLECTED SHOCK VELOCITY (cm/s): ',VRS

      ELSE
         WRITE (LOUT, '(/2(/T24,A,T44,A),/)')
     1   'CONDITONS',        'CONDITONS',
     2   'BEFORE THE SHOCK', 'AFTER THE SHOCK'
C
         WRITE (LOUT, '(1X, A, T24, 1PE10.4, T44,1PE10.4)')
     1   'PRESSURE (atm)  ', P1A,  P2A,
     2   'TEMPERATURE (K) ', T1,   T2,
     3   'DENSITY (g/cc)  ', RHO1, RHO2,
     4   'VELOCITY (cm/s) ', V1,   V2,
     5   'MACH NO.        ', XM1,  XM2
C
         IF (IPRB .EQ. 1) WRITE (LOUT, '(3(/1X,A,1PE10.4))')
     1   'GAS VISCOSITY BEFORE SHOCK (g/cm/s) ', VIS,
     2   'BOUNDARY LAYER PARAMETER BETA       ', BETA,
     3   'LIMITING SEPARATION BETWEEN SHOCK AND CONTACT SURFACE (cm) ',
     4   XLM
      ENDIF
C
C     Call subroutine that solves the problem
      CALL RUNIT (NEQ, LGSCAL, LOUT, LSAVE, IPAR, RPAR, CPAR(IKS),
     1            KERR)
C*****licensing > FLEXlm
C      IF (KERR) CALL CHKIN ( IHFLEX )
C*****END licensing > FLEXlm
      IF (KERR) RETURN
C
C     Read another case
      IF (LCNTUE) GO TO 110
C*****licensing > FLEXlm
C      CALL CHKIN ( IHFLEX )
C*****END licensing > FLEXlm
      RETURN
C
C     end of PROGRAM SHOCK
      END
C
      SUBROUTINE SETWRK (LINKCK, LSAVE, LOUT, LIPAR, LRPAR, LCPAR,
     1                   IPAR, RPAR, CPAR, NEQ, KERR)
C
C  START PROLOGUE
C     SETWRK sets pointers for storage RPAR, IPAR, CPAR
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /ICONS/ IPRB, IGOT, IMOLF, KK, NZZ, NT, NRHO, NKK, NTL,
     1               NV, NA, NCWK, NXX, NYY, NWT, NHMS, NGW, LRW,
     2               ICWK, IKS, ICC, IVODE, NVODE, LIW
      COMMON /RCONS/ HO, XLM, DIA, VISCOF, BETA, RU, PATM, VIS,
     1               P1, T1, RHO1, H1, V1, XM1, P2, T2, RHO2, H2, V2,
     2               XM2, P3, T3, RHO3, H3, V3, XM3, VRS, ASHOCK,
     3               ATOL, RTOL, TINIT, DT, TEND, TSTART
C
      DIMENSION IPAR(*), RPAR(*)
      CHARACTER CPAR(*)*(*)
      LOGICAL KERR
C
      KERR = .FALSE.
      CALL CKLEN (LINKCK, LOUT, LENICK, LENRCK, LENCCK, IFLAG)
      IF (IFLAG .GT. 0) THEN
         KERR = .TRUE.
         RETURN
      ENDIF
C
      IF (LIPAR.GE.LENICK .AND. LRPAR.GE.LENRCK .AND. LCPAR.GE.LENCCK)
     1   THEN
         CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT,
     1                IPAR, RPAR, CPAR, IFLAG)
         IF (IFLAG .GT. 0) THEN
            KERR = .TRUE.
            RETURN
         ENDIF
         CALL CKINDX (IPAR, RPAR, MM, KK, II, NFIT)
      ELSE
         WRITE (LOUT, '(/1X,A,I6)')
     1   'ERROR...RPAR MUST BE DIMENSIONED AT LEAST ', LENICK,
     2   '        IPAR MUST BE DIMENSIONED AT LEAST ', LENRCK,
     3   '        CPAR MUST BE DIMENSIONED AT LEAST ', LENCCK
         KERR = .TRUE.
         RETURN
      ENDIF
C
C     REAL WORK SPACE, RPAR
C
      NZZ  = 1
C            points to the dependent variable vector
      NT   = NZZ
C            points to temperature, the first element
      NRHO = NT + 1
C            points to density, the second element
      NKK  = NRHO + 1
C            points to the species elements
      NTL  = NKK + KK
C            TL pointer
      NV   = NTL + 1
C            volume
      NA   = NV + 1
C            area
      NEQ  = KK + 5
      NCWK = NZZ + NEQ
C            points to the CHEMKIN workspace
      NXX  = NCWK + LENRCK
C            points to the mole fraction vector
      NYY  = NXX + KK
C            points to the mass fraction vector
      NWT  = NYY + KK
C            points to the molecular weight vector
      NHMS = NWT + KK
C            points to the enthalpies in mass units
      NVODE = NHMS + KK
C             points to real VODE workspace
      LRW = 22 + 9*NEQ + 2*NEQ**2
      NTOT = NVODE + LRW - 1
C
C     INTEGER SPACE IPAR
C
      ICWK = 1
C            points to integer CHEMKIN workspace
C
      IVODE = ICWK + LENICK
C             points to integer VODE workspace
      LIW = 30 + NEQ
      ITOT = IVODE + LIW - 1
C
      IKS = 1
C           points to character species names
      ICC = IKS + KK
C           points to character CHEMKIN workspace
      ICTOT = ICC + LENCCK
C
      IF (LRPAR.LE.NTOT .OR. LIPAR.LE.ITOT .OR. LCPAR.LE.ICTOT) THEN
         WRITE (LOUT, '(/1X,A,I6)')
     1   'ERROR...RPAR MUST BE DIMENSIONED AT LEAST ', NTOT,
     2   '        IPAR MUST BE DIMENSIONED AT LEAST ', ITOT,
     3   '        CPAR MUST BE DIMENSIONED AT LEAST ', ICTOT
         KERR = .TRUE.
         RETURN
      ENDIF
C
C     Store CHEMKIN workspace
C
      REWIND (LINKCK)
      CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT, IPAR(ICWK),
     1             RPAR(NCWK), CPAR(ICC), IFLAG)
      IF (IFLAG .GT. 0) THEN
         KERR = .TRUE.
         RETURN
      ENDIF
C
C     Store CHEMKIN variables
      CALL CKSYMS (CPAR(ICC),  LOUT, CPAR(IKS), KERR)
      CALL CKRP   (IPAR(ICWK), RPAR(NCWK), RU, RUC, PATM)
      CALL CKWT   (IPAR(ICWK), RPAR(NCWK), RPAR(NWT))
C
      WRITE (LSAVE) 'CKLINK          '
      CALL CKSAVE (LOUT, LSAVE, IPAR(ICWK), RPAR(NCWK), CPAR(ICC))
C
C     end of SUBROUTINE SETWRK
      RETURN
      END
C
      SUBROUTINE RUNIT (NEQ, LGSCAL, LOUT, LSAVE, IPAR, RPAR, KSYM,
     1                  KERR)
C
C  START PROLOGUE
C     RUNIT loops over the requested time steps, calls the stiff ODE
C     integrator VODE, and prints the results.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /ICONS/ IPRB, IGOT, IMOLF, KK, NZZ, NT, NRHO, NKK, NTL,
     1               NV, NA, NCWK, NXX, NYY, NWT, NHMS, NGW, LRW,
     2               ICWK, IKS, ICC, IVODE, NVODE, LIW
      COMMON /RCONS/ HO, XLM, DIA, VISCOF, BETA, RU, PATM, VIS,
     1               P1, T1, RHO1, H1, V1, XM1, P2, T2, RHO2, H2, V2,
     2               XM2, P3, T3, RHO3, H3, V3, XM3, VRS, ASHOCK,
     3               ATOL, RTOL, TINIT, DT, TEND, TSTART
C
      DIMENSION RPAR(*), IPAR(*)
      CHARACTER KSYM(*)*(*)
      LOGICAL LGSCAL, KERR
      EXTERNAL FUNSHK
C
      KERR = .FALSE.
      ITOL = 1
      MF = 22
      HO = 1.0E-25
C
C     Initialize optional inputs for VODE
      DO 110 N = 5, 10
         RPAR(NVODE + N - 1) = 0.0
         IPAR(IVODE + N - 1) = 0
110   CONTINUE
C
      ISTATE   = 1
      ITASK    = 1
      IOPT     = 1
      RPAR(NVODE+4) = HO
      NLINES   = 100
C
      TT1 = TINIT
      TT2 = TINIT
      IF (LGSCAL) THEN
         TT2   = TINIT
         TINIT = 0.0
         TT2L  = LOG10(TT2)
         TT1L  = TT2L - DT
         TT1   = 10**TT1L
      ENDIF
C
C     Begin timestep loop
130   CONTINUE
C
      IF (NLINES .GE. 55) THEN
         NLINES = 0
C
         WRITE (LOUT, '(1H1,1X,8(A))')
     1   'T(sec)     ',
     2   'TEMP(K)    ',
     3   'PRES(atm)  ',
     4   'RHO(g/cc)  ',
     5   'MEAN WT    ',
     6   'AREA(cm**2)',
     7   'VEL(cm/s)  ',
     8   'LAB TIME(sec)'
C
         NLINES = NLINES + 1
         DO 140 K1 = 1, KK, 7
            WRITE (LOUT,'(13X,7(A))') (KSYM(K)(1:11),K=K1,MIN(K1+6,KK))
            NLINES = NLINES + 1
140      CONTINUE
         WRITE (LOUT, '(1H )')
         NLINES = NLINES + 1
      ENDIF
C
      TT2 = TT1 + DT
      IF (LGSCAL) THEN
         TT1L = LOG10(TT1)
         TT2L = TT1L+DT
         TT2  = 10.0**TT2L
         IF (ISTATE.EQ.1) TT1 = 0.0
      ENDIF
160   CONTINUE
C
C     Retrieve local variables from solution
      T   = RPAR(NT)
      RHO = RPAR(NRHO)
      TL  = RPAR(NTL)
      CALL CKPY (RHO, T, RPAR(NKK), IPAR(ICWK), RPAR(NCWK), P)
C
      IF (IMOLF .EQ. 1) THEN
C        compute mole fractions
         CALL CKYTX (RPAR(NKK), IPAR(ICWK), RPAR(NCWK), RPAR(NXX))
      ELSE
C        compute molar concentrations
         CALL CKYTCP (P, T, RPAR(NKK), IPAR(ICWK), RPAR(NCWK),
     1                RPAR(NXX))
      ENDIF
C
      V  = RPAR(NV)
      CALL CKMMWY (RPAR(NKK), IPAR(ICWK), RPAR(NCWK), WTM)
      AREA = ASHOCK
      IF (IPRB .EQ. 1) THEN
         XX = RPAR(NA)
         AREA = ASHOCK/(1.0-SQRT(XX/XLM))
         IF (XX/XLM .GE. 0.90) THEN
            WRITE (LOUT, '(//10X,A,A)')
     1      '  MAXIMUM FLOW DURATION IS APPROACHED.',
     2      '  PROBLEM IS TERMINATED.'
            RETURN
         ENDIF
      ENDIF
C
C     Write solution
      WRITE (LOUT,  '(8(1PE11.3))')
     1       TT1, T, P/PATM, RHO, WTM, AREA, V, TL
      NLINES = NLINES + 1
C
C     Write post-processing file
      WRITE (LSAVE) TT1, T, P/PATM, RHO, WTM, AREA, V, TL,
     1              (RPAR(NXX + K - 1),K=1,KK)
C
      DO 200 K1 = 1, KK, 7
         WRITE (LOUT, '(11X,7(1PE11.3))')
     1   (RPAR(NXX + K - 1), K=K1,MIN(K1+6,KK))
         NLINES = NLINES + 1
200   CONTINUE
      WRITE (LOUT, '(1H )')
      NLINES = NLINES + 1
         IF (TT1 .GE. TEND) RETURN
210   CONTINUE
C
C     Call integrator
C
C*****precision > single
C      CALL SVODE
C*****END precision > single
C*****precision > double
      CALL DVODE
C*****END precision > double
     *          (FUNSHK, NEQ, RPAR(NZZ), TT1, TT2, ITOL, RTOL, ATOL,
     1           ITASK, ISTATE, IOPT, RPAR(NVODE), LRW, IPAR(IVODE),
     2           LIW, JAC, MF, RPAR, IPAR)
C
      IF (ISTATE .EQ. 2) GO TO 130
      IF (ISTATE .NE. -1) THEN
         WRITE (LOUT, '(/5X,A,I3)') 'ISTATE = ',ISTATE
         KERR = .TRUE.
         RETURN
      ENDIF
      ISTATE = 2
      GO TO 210
C
      END
C
      SUBROUTINE FUNSHK (NEQ, TIME, Z, ZP, RPAR, IPAR)
C
C  START PROLOGUE
C     FUNSHK is where the governing equations are defined;
C     the dependent variable vector Z is passed in from VODE, and
C     the right hand sides of the governing equations are evaluated
C     and passed back through ZP.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /ICONS/ IPRB, IGOT, IMOLF, KK, NZZ, NT, NRHO, NKK, NTL,
     1               NV, NA, NCWK, NXX, NYY, NWT, NHMS, NGW, LRW,
     2               ICWK, IKS, ICC, IVODE, NVODE, LIW
      COMMON /RCONS/ HO, XLM, DIA, VISCOF, BETA, RU, PATM, VIS,
     1               P1, T1, RHO1, H1, V1, XM1, P2, T2, RHO2, H2, V2,
     2               XM2, P3, T3, RHO3, H3, V3, XM3, VRS, ASHOCK,
     3               ATOL, RTOL, TINIT, DT, TEND, TSTART
C
      DIMENSION Z(*), ZP(*), RPAR(*), IPAR(*)
C
      AREA(AUG1) = ASHOCK / (1.0-SQRT(AUG1/XLM))
      DAREA(AUG1,AUG2) = AUG2**2 /2.0 /XLM /ASHOCK /SQRT(AUG1/XLM)
C
C        Z(1)=TEMPERATURE    ZP(1)=D(T)/DT
C        Z(2)=RHO            ZP(2)=D(RHO)/DT
C        Z(K+2)=Y(K)         ZP(K+2)=D(Y(K))/DT
C        Z(KK+3)=LAB TIME    ZP(KK+3)=D(LT)/DT
C        Z(KK+4)=VELOCITY    ZP(KK+4)=D(U)/DT
C        Z(KK+5)=DISTANCE    ZP(KK+5) = D(X)/DT
C
C     Compute local variables:
      T    = Z(NT)
      RHO  = Z(NRHO)
      U    = Z(NV)
      X    = Z(NA)
      A    = ASHOCK
      DADX = 0.0
      U2   = U**2
      CALL CKPY   (RHO, T, Z(NKK), IPAR(ICWK), RPAR(NCWK), P)
      CALL CKMMWY (Z(NKK), IPAR(ICWK), RPAR(NCWK), WTM)
      CALL CKWYP  (P, T, Z(NKK), IPAR(ICWK), RPAR(NCWK), ZP(NKK))
      CALL CKCPBS (T, Z(NKK), IPAR(ICWK), RPAR(NCWK), CPB)
      CALL CKHMS  (T, IPAR(ICWK), RPAR(NCWK), RPAR(NHMS))
      CPBT = CPB*T
C
      IF (IPRB .EQ. 1) THEN
C        Area profile
C
         A = AREA(X)
         IF (X .NE. 0.0) THEN
            DADX = DAREA(X,A)
         ELSE
C
C           at time=0 the derivative function is singular;
C           evaluate DADX at a time of H) seconds (the smallest
C           time step employed).
C
            X0 = HO * V2
            A = AREA(X0)
            DADX = DAREA(X0,A)
         ENDIF
      ENDIF
C
      SUM = 0.0
      DO 130 K = 1, KK
         WT  = RPAR(NWT + K - 1)
         HMS = RPAR(NHMS + K - 1)
         SUM = SUM + WT * ZP(K+2) * (HMS - (WTM/WT)*CPBT)
130   CONTINUE
      T1    = RU * RHO * SUM / (WTM * CPB)
      T2    = (1.0 - RU/(WTM*CPB)) * U2 * RHO**2 * U * DADX/A
      TTHREE    = P * (1.0+U2/CPBT) - RHO*U2
      ZP(NRHO) = (T1+T2) / TTHREE
C
      SUM = 0.0
      DO 140 K = 1, KK
         WT  = RPAR(NWT + K - 1)
         HMS = RPAR(NHMS + K - 1)
         SUM = SUM + WT * ZP(NKK + K - 1) * HMS
         ZP(NKK + K - 1) = WT * ZP(NKK + K - 1) / RHO
140   CONTINUE
      ZP(NT) = (U2*ZP(NRHO) - SUM)/(CPB*RHO) + U2*U*DADX/(A*CPB)
C
      IF (IPRB .EQ. 3) THEN
         ZP(NTL) = 1.0
      ELSE
         ZP(NTL) = RHO1 * ASHOCK/RHO/A
      ENDIF
C
      ZP(NV) = -U * ZP(NRHO)/RHO - U2*DADX/A
      ZP(NA) = U
C
C     end of SUBROUTINE FUNSHK
      RETURN
      END
C
      SUBROUTINE READIN (LIN, LOUT, LGSCAL, KSYM, ICKWRK, RCKWRK,
     1                   X, ITITL, KERR, LCNTUE)
C
C  START PROLOGUE
C     READIN reads keyword input and passes the problem parameters
C     back to the main code.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /ICONS/ IPRB, IGOT, IMOLF, KK, NZZ, NT, NRHO, NKK, NTL,
     1               NV, NA, NCWK, NXX, NYY, NWT, NHMS, NGW, LRW,
     2               ICWK, IKS, ICC, IVODE, NVODE, LIW
      COMMON /RCONS/ HO, XLM, DIA, VISCOF, BETA, RU, PATM, VIS,
     1               P1, T1, RHO1, H1, V1, XM1, P2, T2, RHO2, H2, V2,
     2               XM2, P3, T3, RHO3, H3, V3, XM3, VRS, ASHOCK,
     3               ATOL, RTOL, TINIT, DT, TEND, TSTART
C
      DIMENSION ICKWRK(*), RCKWRK(*), X(*)
      CHARACTER LINE*80, KEY*4, CKCHUP*4, KSYM(*)*(*), ITITL*(*)
      LOGICAL IERR, KERR, LGSCAL, NEC(5), N1(4), N2(4), N3(4), NB(2),
     1        LCNTUE
      EXTERNAL CKCHUP
C
      DATA NEC,N1,N2,N3,NB/19*.FALSE./
C
C     Initialize variables
      KERR = .FALSE.
      LCNTUE = .FALSE.
      DO 110 K = 1, KK
         X(K) = 0.0
110   CONTINUE
C
      IMOLF  = 1
      TINIT  = 0.0
      DT     = 1.0E-6
      LGSCAL = .FALSE.
      RTOL   = 1.E-6
      ATOL   = 1.E-10
      ASHOCK = 1.0
      DIA    = SQRT(4.0*ASHOCK/3.1415926535)
C
      WRITE(LOUT, '(/1X,A,/)') ' KEYWORD INPUT:'
C
      I1 = 0
      I2 = 0
      I5 = 0
C
C     Read input
130   CONTINUE
      KEY = ' '
      LINE = ' '
      IERR = .FALSE.
      READ  (LIN,  '(A)', END=135) LINE
      WRITE (LOUT, '(5X,A)') LINE
      CALL CKDTAB (LINE)
      KEY = CKCHUP(LINE,4)
      LINE(1:4) = ' '
C
C     Determine if a keyword was read
      IF (KEY .EQ. 'TSTR') THEN
C        Initial time
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TINIT, IERR)
C
      ELSEIF (KEY .EQ. 'TEND') THEN
C        Final time
         NEC(1) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TEND, IERR)
C
      ELSEIF (KEY(1:2) .EQ. 'DT') THEN
C        Output time interval
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DT, IERR)
C
      ELSEIF (KEY .EQ. 'LGDT') THEN
C        Log output time interval
         LGSCAL = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DT, IERR)
C
      ELSEIF (KEY .EQ. 'ATOL') THEN
C        Floor value for integration
         CALL CKXNUM (LINE, 1, LOUT, NVAL, ATOL, IERR)
C
      ELSEIF (KEY .EQ. 'RTOL') THEN
C        Error tolerance for integration
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RTOL, IERR)
C
      ELSEIF (KEY(1:2) .EQ. 'T1') THEN
C        Temperature before incident shock
         N1(1) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, T1, IERR)
C
      ELSEIF (KEY(1:3) .EQ. 'P1A') THEN
C        Pressure before incident shock (atm)
         N1(2) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, P1, IERR)
         P1 = P1*PATM
C
      ELSEIF (KEY .EQ. 'RHO1') THEN
C        Density before incident shock
         N1(3) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RHO1, IERR)
C
      ELSEIF (KEY .EQ. 'VSHK') THEN
C        Incident shock velocity
         N1(4) = .TRUE.
         N2(4) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, V1, IERR)
C
      ELSEIF (KEY(1:2) .EQ. 'T2') THEN
C        Temperature behind incident shock
         N2(1) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, T2, IERR)
C
      ELSEIF (KEY(1:3) .EQ. 'P2A') THEN
C        Pressure behind incident shock (atm)
         N2(2) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, P2, IERR)
         P2 = P2*PATM
C
      ELSEIF (KEY .EQ. 'RHO2') THEN
C        Density behind incident shock
         N2(3) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RHO2, IERR)
C
      ELSEIF (KEY(1:2) .EQ. 'T3') THEN
C        Temperature behind reflected shock
         N3(1) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, T3, IERR)
C
      ELSEIF (KEY(1:3) .EQ. 'P3A') THEN
C        Pressure behind reflected shock (atm)
         N3(2) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, P3, IERR)
         P3 = P3*PATM
C
      ELSEIF (KEY .EQ. 'RHO3') THEN
C        Density behind reflected shock
         N3(3) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RHO3, IERR)
C
      ELSEIF (KEY(1:3) .EQ. 'VRS') THEN
C        Reflected shock velocity
         N3(4) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VRS, IERR)
C
      ELSEIF (KEY .EQ. 'CONC') THEN
C        Output in concentration units (mole/cc)
         IMOLF = 0
C
      ELSEIF (KEY .EQ. 'ISHK') THEN
C        Incident shock problem
         NEC(2) = .TRUE.
         IPRB = 2
C
      ELSEIF (KEY .EQ. 'ISKB') THEN
C        Incident shock with boundary layer correction
         NEC(2) = .TRUE.
         IPRB = 1
C
      ELSEIF (KEY .EQ. 'RSHK') THEN
C        Reflected shock problem
         NEC(2) = .TRUE.
         IPRB = 3
C
      ELSEIF (KEY .EQ. 'INIT') THEN
C        Initial mole fractions
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL, VALUE,
     1                IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)') ' ERROR READING KEYWORD '//KEY
         ELSE
            X(KSPEC) = VALUE
         ENDIF
C
      ELSEIF (KEY(1:3) .EQ. 'DIA') THEN
C        Diameter of tube
         NB(1) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DIA, IERR)
         ASHOCK = 3.1415926535*DIA**2/4.0
C
      ELSEIF (KEY .EQ. 'VISC') THEN
C        Viscosity coefficient
         NB(2) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VISCOF, IERR)
C
      ELSEIF (KEY .EQ. 'TITL') THEN
C        Title for problem
         ITITL = ' '
         ITITL = LINE
C
      ELSEIF (KEY(1:3) .EQ. 'END') THEN
C        End of input
         GO TO 480
C
      ELSEIF (KEY .EQ. 'CNTN') THEN
         LCNTUE = .TRUE.
C
      ELSE
         WRITE (LOUT, 750) 'ERROR...ILLEGAL KEYWORD'
         IERR = .TRUE.
      ENDIF
C
C     Read next input
      KERR = KERR.OR.IERR
      GO TO 130
C
C     Check to see if input is complete
C
480   CONTINUE
      KERR = KERR.OR.IERR
C
C     Normalize the mole fractions
      CALL CKNORM (X, KK)
      IF (.NOT. NEC(1)) THEN
         WRITE (LOUT, 750) 'ERROR...FINAL TIME NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
      IF (.NOT. NEC(2)) THEN
         WRITE (LOUT, 750) 'ERROR...PROBLEM TYPE NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF (IPRB .NE. 3) THEN
         IF (IPRB .EQ. 1) THEN
C           Check incident shock with boundary layer
C
            IF (.NOT. NB(1)) THEN
               WRITE (LOUT, 750) 'ERROR...DIAMETER NOT SPECIFIED'
               KERR = .TRUE.
            ENDIF
            IF (.NOT. NB(2)) THEN
               WRITE (LOUT, 750)
     1         'ERROR...VISCOSITY COEFFICIENT NOT SPECIFIED'
               KERR = .TRUE.
            ENDIF
         ENDIF
C
         DO 570 N = 1, 3
            IF (N1(N)) I1 = 1
            IF (N2(N)) I2 = 1
570      CONTINUE
         IF (I1.EQ.1 .AND. I2.EQ.1) THEN
            WRITE (LOUT, 750)
     1'ERROR...CONDITIONS BEFORE AND BEHIND INCIDENT SHOCK SPECIFIED'
            KERR = .TRUE.
         ENDIF
         IF (.NOT. N1(4)) THEN
            WRITE (LOUT, 850)
            KERR = .TRUE.
         ENDIF
C
         IF (KERR) RETURN
C
C        Determine P, T, and RHO for conditions 1 or 2
         IGOT = 1
         CALL SETCON (N1, X, RU, ICKWRK, RCKWRK, T1, P1, RHO1, IER)
         IF (IER .NE. 0) THEN
            CALL SETCON (N2, X, RU, ICKWRK, RCKWRK, T2, P2, RHO2, IER)
            IF (IER .NE. 0) THEN
               WRITE (LOUT, 860)
               KERR = .TRUE.
               RETURN
            ENDIF
            IGOT = 2
         ENDIF
         RETURN
      ENDIF
C
C     IPRB=3, Reflected shock
C
      DO 650 N = 1, 3
         IF (N1(N)) I1 = 1
         IF (N3(N)) I5 = 1
650   CONTINUE
      IF (I1.EQ.1 .AND. I5.EQ.1) THEN
         WRITE (LOUT, 870)
         KERR = .TRUE.
      ENDIF
C
      IF (I5.NE.1 .AND. (.NOT.N1(4))) THEN
         WRITE (LOUT, 850)
         KERR = .TRUE.
      ENDIF
C
      IF (KERR) RETURN
C
      IGOT = 1
      CALL SETCON (N1, X, RU, ICKWRK, RCKWRK, T1, P1, RHO1, IER)
      IF (IER .NE. 0) THEN
         CALL SETCON (N3, X, RU, ICKWRK, RCKWRK, T3, P3, RHO3, IER)
         IF (IER .NE. 0) THEN
            WRITE (LOUT, 880)
            KERR = .TRUE.
            RETURN
         ENDIF
         IGOT = 5
      ENDIF
C
750   FORMAT (10X, A)
850   FORMAT (10X,'ERROR...SHOCK SPEED NOT SPECIFIED')
860   FORMAT (10X,'ERROR...ANY TWO OF (T1,P1,RHO1) ',
     1                 'OR ANY TWO OF (T2,P2,RHO2) MUST BE SPECIFIED')
870   FORMAT (10X,'ERROR...CONDITIONS BEFORE INCIDENT SHOCK AND BEHIND',
     1 ' REFLECTED SHOCK GIVEN')
880   FORMAT (10X,'ERROR...ANY TWO OF (T1,P1,RHO1) ',
     1                 'OR ANY TWO OF (T3,P3,RHO3) MUST BE SPECIFIED')
C
C     end of SUBROUTINE READIN
      RETURN
  135 CONTINUE
      WRITE (LOUT, *) ' END OF KEYWORD INPUT '
      RETURN
      END
C
      SUBROUTINE SETCON (NF, X, RU, ICKWRK, RCKWRK, T, P, RHO, IER)
C
C  START PROLOGUE
C     SETCON determines which state variables were defined by the
C     keyword input and computes the remaining ones; for example,
C     if pressure and temperature were given density is calculated.
C     The routine also does error checks to make sure sufficient data
C     was provided.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION X(*), ICKWRK(*), RCKWRK(*)
      LOGICAL NF(*)
C
      IER = 0
C
      IF (.NOT. NF(1)) THEN
C        Temperature not given
         IF (.NOT.NF(2) .AND. .NOT.NF(3)) THEN
            IER = 1
         ELSE
C           Pressure and density given
            CALL CKMMWX (X, ICKWRK, RCKWRK, WTM)
            T = P * WTM/RHO/RU
         ENDIF
         RETURN
      ENDIF
C
C     Temperature given
      IF (NF(2)) THEN
C        Temperature and pressure given
         CALL CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
      ELSE
         IF (.NOT. NF(3)) THEN
C           Density not given
            IER = 1
         ELSE
C           Temperature and density given
            CALL CKPX (RHO, T, X, ICKWRK, RCKWRK, P)
          ENDIF
      ENDIF
C
C     end of SUBROUTINE SETCON
      RETURN
      END
C
C
      SUBROUTINE SHKIC (LOUT, IPAR, RPAR, KERR)
C
C  START PROLOGUE
C     SHKIC computes conditions immediately after either an
C     incident or a reflected shock from the Rankine-Hugoniot
C     relations;  if needed, it also sets up parameters required
C     to make the boundary layer corrections for incident shock
C     problems.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /ICONS/ IPRB, IGOT, IMOLF, KK, NZZ, NT, NRHO, NKK, NTL,
     1               NV, NA, NCWK, NXX, NYY, NWT, NHMS, NGW, LRW,
     2               ICWK, IKS, ICC, IVODE, NVODE, LIW
      COMMON /RCONS/ HO, XLM, DIA, VISCOF, BETA, RU, PATM, VIS,
     1               P1, T1, RHO1, H1, V1, XM1, P2, T2, RHO2, H2, V2,
     2               XM2, P3, T3, RHO3, H3, V3, XM3, VRS, ASHOCK,
     3               ATOL, RTOL, TINIT, DT, TEND, TSTART
      LOGICAL KERR
      DIMENSION IPAR(*), RPAR(*)
      EXTERNAL FUN1, FUN2, FUN3, FUN3A
C
      KERR = .FALSE.
      RE = 1.0E-06
      AE = 1.0E-06
      CALL CKXTY  (RPAR(NXX), IPAR(ICWK), RPAR(NCWK), RPAR(NYY))
      CALL CKXTY  (RPAR(NXX), IPAR(ICWK), RPAR(NCWK), RPAR(NKK))
      CALL CKMMWX (RPAR(NXX), IPAR(ICWK), RPAR(NCWK), WTM)
C
      IF (IGOT .EQ. 1) THEN
C        Conditions before incident shock given
C
         CALL CKHBMS (T1, RPAR(NKK), IPAR(ICWK), RPAR(NCWK), H1)
         CALL CKCPBL (T1, RPAR(NXX), IPAR(ICWK), RPAR(NCWK), CPB)
         GAMMA1 = CPB / (CPB-RU)
         XM1 = V1 / SQRT(GAMMA1 * RU * T1/WTM)
         IF (XM1 .LE. 1.0) THEN
            WRITE (LOUT, '(//20X,A)')
     1      'SHOCK WAVE SPEED IS LESS THAN MACH 1'
            KERR = .TRUE.
            RETURN
         ENDIF
C
         ALPHA = 1.0
C        ALPHA is the temperature ratio across the shock
C
         XM12 = XM1**2
         GAMP1 = GAMMA1 + 1.0
C EM BETA is the presure ratio across the shock.
C EM        BETA  = (2.0*GAMMA1*XM12 - GAMP1 + 1.0)/GAMP1
C EM 7/15/98:  Correct the following equation according to Eq.5.18a
C EM           in Shapiro's Gas Dynamics book.
C
         BETA = (2.0*GAMMA1*XM12 - (GAMMA1-1.0))/GAMP1
C EM The following equation is based on Eq. 37 in the SHOCK manual
C EM and substituting the above definition of Beta:
         ALPHI = BETA*(2.0/XM12 + GAMMA1-1.0)/GAMP1
         CALL ZROIN (FUN1, ALPHA, ALPHI, RE, AE, IFLAG, IPAR, RPAR)
C
         IF (IFLAG .NE. 1) THEN
            WRITE (LOUT, 230) 'FUN1',IFLAG
            IF (IFLAG .NE. 4) THEN
               KERR = .TRUE.
               RETURN
            ENDIF
         ENDIF
C
         V12 = V1**2
         T2     = ALPHA * T1
C EM The following three lines determine the solution of the quadratic
C EM equation for BETA, based on Eq. 36 of the SHOCK manual.
         B      = -(1.0 + RHO1*V12/P1)
         C      = RHO1 * V12/P1 * ALPHA
         BETA   = 0.5 * (-B + SQRT(B**2 - 4.0*C))
         P2     = BETA * P1
         RHO2   = P2 * WTM/RU/T2
         V2     = RHO1 * V1/RHO2
         CALL CKCPBL (T2, RPAR(NXX), IPAR(ICWK), RPAR(NCWK), CPB)
         GAMMA2 = CPB / (CPB-RU)
         XM2    = V2 / SQRT(GAMMA2 * RU * T2/WTM)
         CALL CKHBMS (T2, RPAR(NKK), IPAR(ICWK), RPAR(NCWK), H2)
C
         IF (IPRB .NE. 3) GO TO 140
C
C        Conditions behind reflected shock
C
         ALPHA = 1.0
         ETA   = RHO2/RHO1
         ETAM1 = ETA - 1.0
         XM12  = XM1**2
C EM The following is Eq. 53 in the SHOCK manual:
         XI    = ETA * (XM12*ETAM1*GAMMA1 + ETA) /
     1           (XM12*ETAM1*(GAMMA1-1.0) + ETA)
C EM Corrected this equation according to Eq. 55 in the SHOCK manual:
         T3OT1 = 1.0 + XM12*(GAMMA1-1.0) * (XI-1.0) * ETAM1
     1               / (XI-ETA)/ETA
         ALPHI = T3OT1*T1/T2
C
         IF (VRS .EQ. 0.0) THEN
            CALL ZROIN (FUN3, ALPHA, ALPHI, RE, AE, IFLAG, IPAR, RPAR)
            IF (IFLAG .NE. 1) THEN
               WRITE (LOUT, 230) 'FUN3',IFLAG
               IF (IFLAG .NE. 4) THEN
                  KERR = .TRUE.
                  RETURN
               ENDIF
            ENDIF
            B    = -(1.0 + RHO2 * (V1-V2)**2 / P2 + ALPHA)
            BETA = 0.5 * (-B + SQRT(B**2 - 4.0*ALPHA))
            VRS  = ALPHA/BETA * (V1-V2)/(1.0 - ALPHA/BETA)
         ELSE
C
            CALL ZROIN (FUN3A, ALPHA, ALPHI, RE, AE, IFLAG, IPAR, RPAR)
            IF (IFLAG .NE. 1) THEN
               WRITE (LOUT, 230) 'FUN3',IFLAG
               IF (IFLAG .NE. 4) THEN
                  KERR = .TRUE.
                  RETURN
               ENDIF
            ENDIF
C
            V2P  = (VRS + V1 - V2)**2
            B    = -(1.0 + RHO2/P2 * V2P)
            BETA = 0.5*(-B +SQRT(B**2 - 4.0*ALPHA*RHO2/P2*V2P))
         ENDIF
C
         T3 = ALPHA * T2
         P3 = BETA  * P2
C         RHO3   = P3*WTM*T3/RU
         RHO3 = P3*WTM/(RU*T3)
         V3     = RHO2/RHO3 * (VRS+V1-V2)
         CALL CKCPBL (T3, RPAR(NXX), IPAR(ICWK), RPAR(NCWK), CPB)
         CALL CKHBMS (T3, RPAR(NKK), IPAR(ICWK), RPAR(NCWK), H3)
         GAMMA5 = CPB / (CPB-RU)
         XM3    = V3 / SQRT(GAMMA5 * RU * T3/WTM)
C
      ELSEIF (IGOT .EQ. 2) THEN
C        Incident shock velocity and conditions behind incident shock
C        given
C
         CALL CKRHOX (P2, T2, RPAR(NXX), IPAR(ICWK), RPAR(NCWK), RHO2)
         CALL CKHBMS (T2, RPAR(NKK), IPAR(ICWK), RPAR(NCWK), H2)
         CALL CKCPBL (T2, RPAR(NXX), IPAR(ICWK), RPAR(NCWK), CPB)
         GAMMA2 = CPB/(CPB-RU)
         ALPHA  = 1.0
         ALPHI  = T2/298.0
C
         CALL ZROIN (FUN2, ALPHA, ALPHI, RE, AE, IFLAG, IPAR, RPAR)
         IF (IFLAG .NE. 1) THEN
            WRITE (LOUT, 230) 'FUN2',IFLAG
            IF (IFLAG .NE. 4) THEN
               KERR = .TRUE.
               RETURN
            ENDIF
         ENDIF
C
         V12 = V1**2
         T1   = T2 / ALPHA
         B    = -(1.0 + RHO2/P2 * V12 * ALPHA)
         C    = RHO2/P2 * V12 * ALPHA**2
         BETA = 0.5 * (-B + SQRT(B**2 - 4.0*C))
         P1   = P2 / BETA
         RHO1 = P1*WTM/RU/T1
         V2   = RHO1 * V1/RHO2
         XM2  = V2 / SQRT(GAMMA2 * RU * T2/WTM)
         CALL CKCPBL (T1, RPAR(NXX), IPAR(ICWK), RPAR(NCWK), CPB)
         GAMMA1 = CPB / (CPB-RU)
         XM1  = V1 / SQRT(GAMMA1 * RU * T1/WTM)
         CALL CKHBMS (T1, RPAR(NKK), IPAR(ICWK), RPAR(NCWK), H1)
         GO TO 140
C
      ELSE
C        Conditions behind reflected shock given
C
         CALL CKHBMS (T3, RPAR(NKK), IPAR(ICWK), RPAR(NCWK), H3)
C
C        Conditions before and behind incident shock not needed;
C        values defaulted to zero
         XM2  = 0.0
         P1   = 0.0
         T1   = 0.0
         RHO1 = 0.0
         H1   = 0.0
         XM1  = 0.0
         V1   = 0.0
         P2   = 0.0
         T2   = 0.0
         RHO2 = 0.0
         H2   = 0.0
         V2   = 0.0
         XM2  = 0.0
      ENDIF
C
C     Initial conditions for integration behind reflected shock
      RPAR(NT)   = T3
      RPAR(NRHO) = RHO3
      RPAR(NTL)  = 0.0
      RPAR(NV)   = 0.0
      RPAR(NA)   = 0.0
C
      RETURN
C
  140 CONTINUE
C     Initial conditions for integration behind incident shock
      RPAR(NT)  = T2
      RPAR(NRHO)= RHO2
      RPAR(NTL) = 0.0
      RPAR(NV)  = V2
      RPAR(NA)  = 0.0
      IF (IPRB .EQ. 2) RETURN
C
C     Boundary layer parameters
C
      ASHOCK = 3.1415926535*DIA**2 / 4.0
      VIS    = VISCOF*(T1/300.0)**0.6756
      CALL CKRHOX (P2, T1, RPAR(NXX), IPAR(ICWK), RPAR(NCWK), RHOW)
      WW     = RHO2/RHO1
      VIS2   = VISCOF*(T2/300.0)**0.6756
      CEO    = (RHO2*VIS2/RHOW/VIS)**0.37
      ZETA   = (GAMMA1+1.0) / (GAMMA1-1.0)
      IF (WW .GE. ZETA) ZETA = WW
      BETA   = 1.59*(1.0+(1.796+0.802*WW)/(ZETA*WW-1.0))*CEO
      XLM    = (DIA*RHO2/4.0/BETA/RHOW)**2/(WW-1.0)*V2*RHOW/VIS
C
      RETURN
C
230   FORMAT (//20X,'IFLAG FOR ',A,' IS',I4)
C
C     end of SUBROUTINE SHKIC
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION FUN1 (ALPHA, IPAR, RPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION FUN1 (ALPHA, IPAR, RPAR)
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C     FUN1 is used by ZROIN as part of the iteration process
C     to determine conditions behind an incident shock given
C     conditions befor the shock.
C
      COMMON /ICONS/ IPRB, IGOT, IMOLF, KK, NZZ, NT, NRHO, NKK, NTL,
     1               NV, NA, NCWK, NXX, NYY, NWT, NHMS, NGW, LRW,
     2               ICWK, IKS, ICC, IVODE, NVODE, LIW
      COMMON /RCONS/ HO, XLM, DIA, VISCOF, BETA, RU, PATM, VIS,
     1               P1, T1, RHO1, H1, V1, XM1, P2, T2, RHO2, H2, V2,
     2               XM2, P3, T3, RHO3, H3, V3, XM3, VRS, ASHOCK,
     3               ATOL, RTOL, TINIT, DT, TEND, TSTART
C
      DIMENSION IPAR(*), RPAR(*)
C
      V12 = V1**2
      T    = ALPHA*T1
      CALL CKHBMS (T, RPAR(NYY), IPAR(ICWK), RPAR(NCWK), H)
      B    = -(1.0 + RHO1 * V12/P1)
      C    = RHO1 * V12/P1 * ALPHA
      BETA = 0.5 * (-B + SQRT(B**2 - 4.0*C))
      FUN1 = H1 + 0.5*V12 * (1.0 - ALPHA**2/BETA**2) - H
C
C     end of FUNCTION FUN1
      RETURN
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION FUN2 (ALPHA, IPAR, RPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION FUN2 (ALPHA, IPAR, RPAR)
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C     FUN2 is used by ZROIN as part of the iteration process to
C     determine conditions about an incident shock given the
C     shock velocity and state variables behind the shock.
C
      COMMON /ICONS/ IPRB, IGOT, IMOLF, KK, NZZ, NT, NRHO, NKK, NTL,
     1               NV, NA, NCWK, NXX, NYY, NWT, NHMS, NGW, LRW,
     2               ICWK, IKS, ICC, IVODE, NVODE, LIW
      COMMON /RCONS/ HO, XLM, DIA, VISCOF, BETA, RU, PATM, VIS,
     1               P1, T1, RHO1, H1, V1, XM1, P2, T2, RHO2, H2, V2,
     2               XM2, P3, T3, RHO3, H3, V3, XM3, VRS, ASHOCK,
     3               ATOL, RTOL, TINIT, DT, TEND, TSTART
C
      DIMENSION IPAR(*), RPAR(*)
C
      V12 = V1**2
      T    = T2/ALPHA
      CALL CKHBMS (T, RPAR(NYY), IPAR(ICWK), RPAR(NCWK), H)
      B    = -(1.0 + RHO2 * V12 * ALPHA/P2)
      C    = RHO2 * V12 * ALPHA**2 / P2
      BETA = 0.5 * (-B + SQRT(B**2 - 4.0*C))
      FUN2 = H + 0.5 * V12 * (1.0 - ALPHA**2/BETA**2) - H2
C
C     end of FUNCTION FUN2
      RETURN
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION FUN3 (ALPHA, IPAR, RPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION FUN3 (ALPHA, IPAR, RPAR)
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C     FUN3 is used by ZROIN as part of the iteration process
C     to determine conditions behind a reflected shock given
C     conditions before the indicent shock.
C
      COMMON /ICONS/ IPRB, IGOT, IMOLF, KK, NZZ, NT, NRHO, NKK, NTL,
     1               NV, NA, NCWK, NXX, NYY, NWT, NHMS, NGW, LRW,
     2               ICWK, IKS, ICC, IVODE, NVODE, LIW
      COMMON /RCONS/ HO, XLM, DIA, VISCOF, BETA, RU, PATM, VIS,
     1               P1, T1, RHO1, H1, V1, XM1, P2, T2, RHO2, H2, V2,
     2               XM2, P3, T3, RHO3, H3, V3, XM3, VRS, ASHOCK,
     3               ATOL, RTOL, TINIT, DT, TEND, TSTART
C
      DIMENSION IPAR(*), RPAR(*)
C
      V12 = (V1 - V2)**2
      T    = ALPHA*T2
      CALL CKHBMS (T, RPAR(NYY), IPAR(ICWK), RPAR(NCWK), H)
      B    = -(1.0 + RHO2/P2 * V12 + ALPHA)
      BETA = 0.5 * (-B + SQRT(B**2 - 4.0*ALPHA))
      AB = ALPHA/BETA
      FUN3 = H2 + 0.5 * V12 / (1.0-AB)*(1.0+AB) - H
C
C     end of FUNCTION FUN3
      RETURN
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION FUN3A (ALPHA, IPAR, RPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION FUN3A (ALPHA, IPAR, RPAR)
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C     FUN3A is used by ZROIN as part of the iteration process
C     to determine conditions behind a reflected shock given
C     conditions befor ethe incident shock, and the reflected
C     shock velocity.
C
      COMMON /ICONS/ IPRB, IGOT, IMOLF, KK, NZZ, NT, NRHO, NKK, NTL,
     1               NV, NA, NCWK, NXX, NYY, NWT, NHMS, NGW, LRW,
     2               ICWK, IKS, ICC, IVODE, NVODE, LIW
      COMMON /RCONS/ HO, XLM, DIA, VISCOF, BETA, RU, PATM, VIS,
     1               P1, T1, RHO1, H1, V1, XM1, P2, T2, RHO2, H2, V2,
     2               XM2, P3, T3, RHO3, H3, V3, XM3, VRS, ASHOCK,
     3               ATOL, RTOL, TINIT, DT, TEND, TSTART
C
      DIMENSION IPAR(*), RPAR(*)
C
      T = ALPHA * T2
      CALL CKHBMS (T, RPAR(NYY), IPAR(ICWK), RPAR(NCWK), H)
      V2P = (VRS + V1 - V2)**2
      B = -(1.0 + RHO2/P2 * V2P)
      BETA = 0.5 * (-B + SQRT(B**2 - 4.0*ALPHA*RHO2/P2*V2P))
      FUN3A = H2 + 0.5 * V2P * (1.0 - ALPHA**2/BETA**2) - H
C
C     end of FUNCTION FUN3A
      RETURN
      END
C
      SUBROUTINE ZROIN (F, B, C, RE, AE, IFLAG, IPAR, RPAR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION IPAR(*),RPAR(*)
      EXTERNAL F
      DATA ONE,ZERO/1.0,0.0/
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2646
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 8.1  AUGUST 1980
C                   *************************
C                   *       ISSUED BY       *
C                   *  SANDIA LABORATORIES, *
C                   *   A PRIME CONTRACTOR  *
C                   ********     TO THE     *
C                          *  UNITED STATES *
C                          *   DEPARTMENT   *
C                          *       OF       *
C                          *     ENERGY     *
C      *********************  ---NOTICE---  *********************
C      *THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED*
C      *  BY THE UNITED STATES GOVERNMENT.  NEITHER THE UNITED  *
C      *   STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,   *
C      *               NOR ANY OF THEIR EMPLOYEES,              *
C      * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR *
C      * EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR  *
C      * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE  *
C      *          **********    ACCURACY,   **********          *
C      *          *        *  COMPLETENESS  *        *          *
C      *          *        *  OR USEFULNESS *        *          *
C      *          *        *     OF ANY     *        *          *
C      *          *        *  INFORMATION,  *        *          *
C      *          *        *   APPARATUS,   *        *          *
C      *       ****        *     PRODUCT    *        ****       *
C      *       *           *   OR PROCESS   *           *       *
C      *       *           *   DISCLOSED,   *           *       *
C      *       *           *  OR REPRESENTS *           *       *
C      *       *          **    THAT ITS    **          *       *
C      *       *          **  USE WOULD NOT **          *       *
C      *********          **    INFRINGE    **          *********
C                         **    PRIVATELY   **
C                         **      OWNED     **
C                         **     RIGHTS.    **
C                         **                **
C                         **                **
C                         **                **
C                         ********************
C
C     BASED ON A METHOD BY T J DEKKER
C     WRITTEN BY L F SHAMPINE AND H A WATTS
C     MODIFIED FOR THE MATH LIBRARY BY C B BAILEY
C
C     ABSTRACT
C        ZROIN SEARCHES FOR A ZERO OF A FUNCTION F(X) BETWEEN
C        THE GIVEN VALUES B AND C UNTIL THE WIDTH OF THE INTERVAL
C        (B,C) HAS COLLAPSED TO WITHIN A TOLERANCE SPECIFIED BY
C        THE STOPPING CRITERION, ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).
C        THE METHOD USED IS AN EFFICIENT COMBINATION OF BISECTION AND
C        THE SECANT RULE.  IN ORDER TO INSURE THAT ZROIN WILL CONVERGE
C        TO A ZERO, THE USER SHOULD PICK VALUES FOR B AND C AT WHICH
C        THE FUNCTION DIFFERS IN SIGN.
C
C     DESCRIPTION OF ARGUMENTS
C     F,B,C,RE AND AE ARE INPUT PARAMETERS
C     B,C AND IFLAG ARE OUTPUT PARAMETERS
C        F     - NAME OF THE REAL VALUED EXTERNAL FUNCTION.  THIS NAME
C                MUST BE IN AN EXTERNAL STATEMENT IN THE CALLING
C                PROGRAM.  F MUST BE A FUNCTION OF ONE REAL ARGUMENT.
C        B     - ONE END OF THE INTERVAL (B,C).  THE VALUE RETURNED FOR
C                B USUALLY IS THE BETTER APPROXIMATION TO A ZERO OF F.
C        C     - THE OTHER END OF THE INTERVAL (B,C)
C        RE    - RELATIVE ERROR USED FOR RW IN THE STOPPING CRITERION.
C                IF THE REQUESTED RE IS LESS THAN MACHINE PRECISION,
C                THEN RW IS SET TO APPROXIMATELY MACHINE PRECISION.
C        AE    - ABSOLUTE ERROR USED IN THE STOPPING CRITERION.  IF THE
C                GIVEN INTERVAL (B,C) CONTAINS THE ORIGIN, THEN A
C                NONZERO VALUE SHOULD BE CHOSEN FOR AE.
C        IFLAG - A STATUS CODE.  USER MUST CHECK IFLAG AFTER EACH CALL.
C                CONTROL RETURNS TO THE USER FROM ZROIN IN ALL CASES.
C                XERROR DOES NOT PROCESS DIAGNOSTICS IN THESE CASES.
C                 1 B IS WITHIN THE REQUESTED TOLERANCE OF A ZERO.
C                   THE INTERVAL (B,C) COLLAPSED TO THE REQUESTED
C                   TOLERANCE, THE FUNCTION CHANGES SIGN IN (B,C), AND
C                   F(X) DECREASED IN MAGNITUDE AS (B,C) COLLAPSED.
C                 2 F(B) = 0.  HOWEVER, THE INTERVAL (B,C) MAY NOT HAVE
C                   COLLAPSED TO THE REQUESTED TOLERANCE.
C                 3 B MAY BE NEAR A SINGULAR POINT OF F(X).
C                   THE INTERVAL (B,C) COLLAPSED TO THE REQUESTED
C                   TOLERANCE AND THE FUNCTION CHANGES SIGN IN (B,C) BUT
C                   F(X) INCREASED IN MAGNITUDE AS (B,C) COLLAPSED,I.E.
C                     ABS(F(B OUT)) .GT. MAX(ABS(F(B IN)),ABS(F(C IN)))
C                 4 NO CHANGE IN SIGN OF F(X) WAS FOUND ALTHOUGH THE
C                   INTERVAL (B,C) COLLAPSED TO THE REQUESTED TOLERANCE.
C                   THE USER MUST EXAMINE THIS CASE AND DECIDE WHETHER
C                   B IS NEAR A LOCAL MINIMUM OF F(X), OR B IS NEAR A
C                   ZERO OF EVEN MULTIPLICITY, OR NEITHER OF THESE.
C                 5 TOO MANY (.GT. 500) FUNCTION EVALUATIONS USED.
C
C     REFERENCES
C       1.  L F SHAMPINE AND H A WATTS, ZEROIN, A ROOT-SOLVING CODE,
C           SC-TM-70-631, SEPT 1970.
C       2.  T J DEKKER, FINDING A ZERO BY MEANS OF SUCCESSIVE LINEAR
C           INTERPOLATION, *CONSTRUCTIVE ASPECTS OF THE FUNDAMENTAL
C           THEOREM OF ALGEBRA*, EDITED BY B DEJON AND P HENRICI, 1969.
C
C     ER IS TWO TIMES THE COMPUTER UNIT ROUNDOFF VALUE WHICH IS
C     DEFINED HERE BY THE FUNCTION D1MACH.
C
C*****precision > double
      ER = 2.0 * D1MACH(4)
C*****END precision > double
C*****precision > single
C      ER = 2.0 * R1MACH(4)
C*****END precision > single
C
C     INITIALIZE
C
      RW = MAX(RE,ER)
      AW = MAX(AE,ZERO)
      IC = 0
      ACBS = ABS(B-C)
      A = C
      T = A
      FA = F(T,IPAR,RPAR)
      T = B
      FB = F(T,IPAR,RPAR)
      FC = FA
      KOUNT = 2
      FX = MAX(ABS(FB),ABS(FC))
C
    1 CONTINUE
      IF (ABS(FC) .LT. ABS(FB)) THEN
C        PERFORM INTERCHANGE
         A  = B
         FA = FB
         B  = C
         FB = FC
         C  = A
         FC = FA
      ENDIF
C
      IF (FB .EQ. ZERO) THEN
         IFLAG = 2
         RETURN
      ENDIF
      CMB  = 0.5*(C-B)
      ACMB = ABS(CMB)
      TOL  = RW*ABS(B)+AW
C
C     TEST STOPPING CRITERION
      IF (ACMB .LE. TOL) GO TO 10
C
C     CALCULATE NEW ITERATE IMPLICITLY AS B+P/Q
C     WHERE WE ARRANGE P .GE. 0.
C     THE IMPLICIT FORM IS USED TO PREVENT OVERFLOW.
      P = (B-A)*FB
      Q = FA - FB
      IF (P .LE. ZERO) THEN
         P = -P
         Q = -Q
      ENDIF
C
C     UPDATE A AND CHECK FOR SATISFACTORY REDUCTION
C     IN THE SIZE OF OUR BOUNDING INTERVAL.
      A  = B
      FA = FB
      IC = IC + 1
      IF (IC .GE. 4) THEN
         IF (8.0*ACMB .GE. ACBS) GO TO 6
         IC   = 0
         ACBS = ACMB
      ENDIF
C
C     TEST FOR TOO SMALL A CHANGE
C
      IF (P .LE. ABS(Q)*TOL) THEN
C
C     INCREMENT BY TOLERANCE
         B = B + SIGN(TOL,CMB)
         GO TO 7
      ENDIF
C
C     ROOT OUGHT TO BE BETWEEN B AND (C+B)/2.0
C
      IF (P .LT. CMB*Q) THEN
C
C     INTERPOLATE
         B = B + P/Q
         GO TO 7
      ENDIF
C
    6 CONTINUE
      B = 0.5*(C+B)
C     BISECT
C
C     HAVE COMPLETED COMPUTATION FOR NEW ITERATE B
    7 CONTINUE
      T = B
      FB = F(T,IPAR,RPAR)
      IF (FB .EQ. ZERO) THEN
         IFLAG = 2
         RETURN
      ENDIF
C
C     DECIDE WHETHER NEXT STEP IS INTERPOLATION OR EXTRAPOLATION
      IF (SIGN(ONE,FB) .EQ. SIGN(ONE,FC)) THEN
         C = A
         FC = FA
      ENDIF
C
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. 500) THEN
         IFLAG = 5
         RETURN
      ENDIF
      GO TO 1
C
C     FINISHED. PROCESS RESULTS FOR PROPER SETTING OF IFLAG
C
   10 CONTINUE
      IF (SIGN(ONE,FB) .EQ. SIGN(ONE,FC)) THEN
         IFLAG = 4
         RETURN
      ENDIF
      IF (ABS(FB) .GT. FX) THEN
         IFLAG = 3
         RETURN
      ENDIF
      IFLAG = 1
C
C     end of SUBROUTINE ZEROIN
      RETURN
      END
