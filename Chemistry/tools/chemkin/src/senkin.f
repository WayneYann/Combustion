C  CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      SUBROUTINE SENKIN (LIN, LOUT, LINKCK, LSAVE, LREST,
     1                   LR, R, LI, I, LC, C)
C
C  START PROLOGUE
C     This subroutine computes the spontaneous ignition delay
C     for an adiabatic gas mixture in a closed system with either:
C     (1) pressure constant,
C     (2) volume a function of time,
C     (3) pressure constant and temperature a function of time.
C     The CHEMKIN package is used to handle chemical kinetics.
C     The time integration and kinetic sensitivity analysis
C     employ program DASAC.

C     Author:  Andy Lutz,
C              MS-9042
C              Computational Mechanics Division
C              Sandia National Laboratories, Livermore, CA, 94551-0969
C              (510) 294-2761
C              aelutz@ca.sandia.gov
C
C     V. 3.9 98/01/24 (E. Meeks)
C     1. Fix bug#148: Remove definition of P after call to TEMPT in
C        subroutine RTEMP, since PA is undefined here.  Add definition 
C        of PA = P/PATM just before call to TEMPT.
C     V. 3.8 97/10/27
C     1. Fix bug#006: move DATA statement after COMMON
C     2. Fix bug#112: eliminate unused variables ZERO and IREAC
C     V. 3.7 97/03/01
C     1. move OPEN statements to main "driver" program,
C        replace STOP with RETURN, etc.,
C        to satisfy Reaction Design requirement of providing object
C        files instead of source code.
C     V. 3.6, 96/11/05
C     1. merged with A. Lutz's latest version 3.4a, which uses
C        keyword TRST to select a restart solution time.
C     2. added SUBROUTINE SETPAR to contain tolerance initialization
C     3. initialize TLIM=0.0 in SUBROUTINE REDKEY
C     V. 3.5, 96/05/24
C     1. initial sccs version
C     VERSION 3.4 (4/10/96 F. Rupley)
C     1. CHEMKIN-III initial version; set up for ascii linkfile
C     2. use explicit dimensioning if possible
C     3. initialize values in work arrays
C     4. use SUBROUTINE RDSAVE for restarts
C     5. simplify error processing in REDKEY
C     CHANGES FOR VERSION 3.3 (1/16/96 A. Lutz)
C     1.  Add ICASE=6 for constant T,V.
C     CHANGES FOR VERSION 3.2 (3/6/95 A. Lutz)
C     1.  Merge version 3.1 into 2.0a
C     2.  Generalize RTEMP to allow P(t).
C     CHANGES FOR VERSION 3.1 (1/19/95 F. Rupley)
C     1.  Add integer error flag to CKLEN, CKINIT call lists.
C     CHANGES FOR VERSION 3.0 (3/15/94 F. Rupley)
C     1.  DOS/PC compatibility effort includes adding file names to
C         OPEN statements, removing unused variables in CALL lists,
C         unusued but possibly initialized variables.
C     CHANGES FOR VERSION 2.0a
C     1.  Allow restart using input T,P (keyword USET). (2/16/95)
C     CHANGES FOR VERSION 2.0
C     1.  Place Chemkin linking file in binary save file.
C     2.  Removed extra printing to unit LIGN.  Note elimination from
C         SENKIN call list requires change to driver.
C     3.  Removed UNICOS blocks for time remaining.
C     4.  Add keyword DTSV for a time increment on printing to the
C         binary file.                                    (9/29/94)
C     CHANGES FOR VERSION 1.8
C     1.  Place reaction constant (D) at front of RPAR (remove A,B,E).
C     2.  Call CKRDEX instead of CKA to perturb reaction constant.
C     CHANGES FOR VERSION 1.7
C     1.  Call CKLEN and CKINDX instead of reading the linking file,
C         and changed position of some pointers in I, R, and C
C         arrays
C     CHANGES FOR VERSION 1.6
C     1.  Call CKRHOY to compute variable density in RCONT
C     CHANGES FOR VERSION 1.5
C     1.  CALL CKLEN will check linking file version and precision
C     2.  Read statements for new linking file
C     CHANGES FROM VERSION 1.3
C     1.  modified OPEN, 'tremain' statements for unicos
C     CHANGES FROM VERSION 1.2
C     1.  read KERR and MAXTB from linking file
C     CHANGES FROM VERSION 1.1
C     1.  changed REAL*8 to DOUBLE PRECISION
C     2.  read linking file for array lengths
C     3.  remove need for COMMON/CKSTRT/
C     4.  allow upper/lower case input
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION R(LR), I(LI)
      CHARACTER C(LC)*(*)
      LOGICAL LSENS, KERR
      COMMON /POINT/ IPICK, IPRCK, IPWT, IPWDOT, IPU, IPRD
C
      WRITE (LOUT, 15)
   15 FORMAT(
     1/' SENKIN:  Sensitity Analysis',
     2/'          Author: Andy Lutz',
     3/'          CHEMKIN-III Version 3.9, 98/01/24',
C*****precision > double
     4/'          DOUBLE PRECISION')
C*****END precision > double
C*****precision > single
C     4/'          SINGLE PRECISION')
C*****END precision > single
C
C     determine mechanism size
      IPICK = 20
      CALL CKLEN (LINKCK, LOUT, LENI, LENR, LENC, IFLAG)
      IF (LENI .LE. LI .AND. LENR.LE.LR .AND. LENC.LE.LC) THEN
         CALL CKINIT (LI, LR, LC, LINKCK, LOUT, I(IPICK), R, C, IFLAG)
         IF (IFLAG .GT. 0) RETURN
         CALL CKINDX (I, R, MM, KK, II, NFIT)
      ELSE
         LITOT = LENI
         LRTOT = LENR
         LCTOT = LENC
         GO TO 555
      ENDIF
C
C     read sensitivity and problem choice keywords
      KERR = .FALSE.
      CALL REDSEN (ICASE, LIN, LOUT, LSENS, KERR)
      IF (KERR) RETURN
C
C     compute DASAC workspace
      IF (ICASE.LT.1 .OR. ICASE.GT.6) THEN
         WRITE (LOUT, '(/1X,A,/)') ' Stop, ICASE not found in SENKIN.'
         RETURN
      ELSEIF (ICASE .GT. 3) THEN
         NSYS = KK
      ELSE
         NSYS  = KK + 1
      ENDIF
C
      IF (LSENS) THEN
         NEQ    = NSYS * (II + 1)
         LSDAS  = NSYS
      ELSE
         NEQ = NSYS
         LSDAS  = 1
      ENDIF
      LIDAS  = 20 + NEQ
      LRDAS  = 40 + NEQ*9 + NSYS*NSYS
C
C     set RPAR pointers (note:  must have IPRD 1st in RPAR!)
      IPRD   = 1
      IPRCK  = IPRD  + II
      IPWT   = IPRCK + LENR
      IPWDOT = IPWT  + KK
      IPU    = IPWDOT+ KK
      IPTOT  = IPU   + KK
C
C     apportion real workspace
      NRPAR  = 1
      NRDAS  = NRPAR + IPTOT
      NSDAS  = NRDAS + LRDAS
      NATOL  = NSDAS + LSDAS
      NRTOL  = NATOL + NEQ
      NXMOL  = NRTOL + NEQ
      NZ     = NXMOL + KK
      NZP    = NZ    + NEQ
      LRTOT  = NZP   + NEQ
C
C     apportion integer workspace
      NIPAR  = 1
      NIDAS  = NIPAR + LENI + IPICK
      LITOT  = NIDAS + LIDAS
C
C     apportion character workspace
      IPCCK = 1
      NKSYM = IPCCK + LENC
      LCTOT = NKSYM + KK
C
  555 CONTINUE
C     check for sufficient space
      WRITE (LOUT, 7020) LI, LITOT, LR, LRTOT, LC, LCTOT
 7020 FORMAT (/, '                Working Space Requirements',
     1        /, '                 Provided        Required ',
     2        /, ' Integer  ', 2I15,
     3        /, ' Real     ', 2I15,
     4        /, ' Character', 2I15, /)
C
      IF (LRTOT.GT.LR .OR. LITOT.GT.LI .OR. LCTOT.GT.LC) THEN
         WRITE (LOUT, *) '  Stop, not enough work space provided.'
         RETURN
      ENDIF
C
C      go to main level
      CALL BEGIN (NSYS, ICASE, II, KK, LENI, LENR, LENC, LINKCK,
     1            LIN, LOUT, LSAVE, LREST, LSENS, LIDAS, LRDAS,
     2            LSDAS, I(NIDAS), R(NRDAS), R(NSDAS), R(NRPAR),
     3            I(NIPAR), R(NZ), R(NZP), R(NRTOL), R(NATOL),
     4            R(NXMOL), C(NKSYM), C(IPCCK))
C
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE BEGIN (NSYS, ICASE, II, KK, LENICK, LENRCK, LENCCK,
     1                  LINKCK, LIN, LOUT, LSAVE, LREST, LSENS, LIDAS,
     2                  LRDAS, LSDAS, IDWORK, DWORK, SDWORK, RPAR, IPAR,
     3                  Z, ZP, RTOL, ATOL, XMOL, KSYM, CCKWRK)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), XMOL(KK), DWORK(LRDAS), IDWORK(LIDAS), 
     1          SDWORK(LSDAS), RTOL(*), ATOL(*), RPAR(*), IPAR(*),
     2          TOLS(4), INFO(15), ISEN(5)
      CHARACTER CCKWRK(LENCCK)*(*), KSYM(KK)*(*), ILINK*16, CSENK*16
      PARAMETER (ILINK='CKLINK', CSENK='SENKIN SOLUTION')
      EXTERNAL RCONP, RCONV, RCONT, RVOLT, RTEMP, RCONTV
      LOGICAL LSENS, RESTRT, IERR, USEINP, KERR
C
C     physical variables in COMMON
      COMMON /RES1/ P
      COMMON /RES2/ TOTMAS
      COMMON /RES3/ T
      COMMON /RES4/ RHO
      COMMON /ATM/ PATM
C
C     pointers in COMMON
      COMMON /POINT/ IPICK, IPRCK, IPWT, IPWDOT, IPU, IPRD
      IPAR(1) = KK
      IPAR(2) = IPRCK
      IPAR(3) = IPRD
      IPAR(6) = IPWT
      IPAR(7) = IPWDOT
      IPAR(8) = IPU
      IPAR(9) = IPICK
      IPAR(10) = LOUT
      IPAR(11) = II
C
      DATA ISEN/5*0/, INFO/15*0/
C
C     initialize CHEMKIN (reset CHEMKIN workspace)
      CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT,
     1             IPAR(IPICK), RPAR(IPRCK), CCKWRK, IFLAG )
      IF (IFLAG .GT. 0) RETURN
C
      CALL CKSYMS (CCKWRK, LOUT, KSYM, IERR)
      CALL CKWT   (IPAR(IPICK), RPAR(IPRCK), RPAR(IPWT))
      CALL CKRP   (IPAR(IPICK), RPAR(IPRCK), RU, RUC, PATM)
C
C     write the CHEMKIN linking info to binary file
      WRITE (LSAVE) ILINK
      CALL CKSAVE (LOUT, LSAVE, IPAR(IPICK), RPAR(IPRCK), CCKWRK)
C
C     load the reaction parameters into RPAR
      DO 15 I = 1, II
        CALL CKRDEX (I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C     read keyword input
      KERR = .FALSE.
      CALL REDKEY (DTPRNT, DTSAVE, KK, KSYM, LIN, LOUT, PA, RESTRT,
     1             T, TLIM, TOLS, TRES, TSTOP, XMOL, USEINP, TRST,
     2             KERR)
      IF (KERR) RETURN
C
C       PHYSICAL INITIAL CONDITIONS
C
C     CASES:
C     ICASE = 1 == CONSTANT PRESSURE
C     ICASE = 2 == CONSTANT VOLUME
C     ICASE = 3 == VOLUME GIVEN IN TIME
C     ICASE = 4 == CONSTANT TEMPERATURE, PRESSURE
C     ICASE = 5 == TEMPERATURE GIVEN IN TIME
C     ICASE = 6 == CONSTANT TEMPERATURE, VOLUME
C
      IF (ICASE .LT. 1 .OR. ICASE .GT. 6) THEN
         WRITE (LOUT, '(/1X,A,/)') ' Stop, ICASE not found in BEGIN.'
         RETURN
      ENDIF
C
      IF (RESTRT) THEN
C        read restart data
C
         CALL RDSAVE (LREST, TRST, LOUT, KK, ICASE, ILINK, CSENK,  
     1                IPAR(IPICK), RPAR(IPRCK), XMOL, Z, TFIL, TMPFIL, 
     2                PFIL, IERR)
         IF (IERR) RETURN
         IF (TRES .LT. 0.) THEN
            TIM = TFIL
         ELSE
            TIM = TRES
         ENDIF
         IF (USEINP) THEN
C           option to use input (T,P) for restart
            WRITE (LOUT,'(/5X,A/)')
     1      'Restart solution will use input (T,P).'
         ELSE
            T = TMPFIL
            P = PFIL
         ENDIF
C
      ELSE
C        original job
C
         TIM = 0.E+0
         IF (ICASE .LE. 3) THEN
            Z(1) = T
            CALL CKXTY  (XMOL, IPAR(IPICK), RPAR(IPRCK), Z(2))
         ELSE
            CALL CKXTY  (XMOL, IPAR(IPICK), RPAR(IPRCK), Z(1))
         ENDIF
      ENDIF
C
      IF (ICASE .EQ. 5) CALL TEMPT (TIM, T, PA)
C
      P = PA * PATM
      IF (ICASE .LE. 3) THEN
         CALL CKRHOY (P, T, Z(2), IPAR(IPICK), RPAR(IPRCK), RHO)
      ELSE
         CALL CKRHOY (P, T, Z(1), IPAR(IPICK), RPAR(IPRCK), RHO)
      ENDIF
C
      IF (RHO .LE. 0.0) THEN
         WRITE (LOUT, '(/1X,A,/)') ' Stop, initial density < 0.'
         RETURN
      ENDIF
C
C     print initial conditions
      IF (RESTRT) WRITE (LOUT, 7000)
C
      IF (ICASE .EQ. 1) THEN
         WRITE (LOUT, 7111)
      ELSEIF (ICASE .EQ. 2) THEN
         WRITE (LOUT, 7112)
      ELSEIF (ICASE .EQ. 3) THEN
         WRITE (LOUT, 7113)
      ELSEIF (ICASE .EQ. 4) THEN
         WRITE (LOUT, 7114)
      ELSEIF (ICASE .EQ. 5) THEN
         WRITE (LOUT, 7115)
      ELSEIF (ICASE .EQ. 6) THEN
         WRITE (LOUT, 7116)
      ENDIF
C
      WRITE (LOUT, 7103)
      WRITE (LOUT, 7100) PA, T, RHO
      IF (ICASE .EQ. 3) THEN
         CALL VOLT (TIM, VZERO, DVDT)
         TOTMAS = RHO * VZERO
         WRITE (LOUT, 7104) TOTMAS
      ENDIF
C
      WRITE (LOUT, 7101)
      DO 130 K = 1, KK
         WRITE (LOUT, 7102) KSYM(K), XMOL(K)
130   CONTINUE
C
      CALL SETPAR (INFO, NSYS, TOLS, RTOL, ATOL, LSENS, ISEN, II)
C     integration routine to run problem
C
      IF (ICASE .EQ. 1) THEN
         CALL SENS13 (RCONP, NSYS, KK, II, DTPRNT, DTSAVE, TSTOP,
     1                TIM, PATM, TLIM, LOUT, LSAVE, LSENS, LIDAS,
     2                LRDAS, LSDAS, Z, ZP, DWORK, IDWORK, SDWORK, RPAR,
     3                IPAR, ATOL, RTOL, XMOL, KSYM, CSENK, INFO, ISEN)
      ELSEIF (ICASE .EQ. 2) THEN
         CALL SENS13 (RCONV, NSYS, KK, II, DTPRNT, DTSAVE, TSTOP,
     1                TIM, PATM, TLIM, LOUT, LSAVE, LSENS, LIDAS,
     2                LRDAS, LSDAS, Z, ZP, DWORK, IDWORK, SDWORK, RPAR,
     3                IPAR, ATOL, RTOL, XMOL, KSYM, CSENK, INFO, ISEN)
      ELSEIF (ICASE .EQ. 3) THEN
         CALL SENS13 (RVOLT, NSYS, KK, II, DTPRNT, DTSAVE, TSTOP,
     1                TIM, PATM, TLIM, LOUT, LSAVE, LSENS, LIDAS,
     2                LRDAS, LSDAS, Z, ZP, DWORK, IDWORK, SDWORK, RPAR,
     3                IPAR, ATOL, RTOL, XMOL, KSYM, CSENK, INFO, ISEN)
      ELSEIF (ICASE .EQ. 4) THEN
         CALL SENS46 (RCONT, NSYS, KK, II, DTPRNT, DTSAVE, TSTOP, TIM,
     1                PATM, TLIM, LOUT, LSAVE, LSENS, LIDAS, LRDAS,
     2                LSDAS, Z, ZP, DWORK, IDWORK, SDWORK, RPAR, IPAR,
     3                ATOL, RTOL, XMOL, KSYM, CSENK, INFO, ISEN)
      ELSEIF (ICASE .EQ. 5) THEN
         CALL SENS46 (RTEMP, NSYS, KK, II, DTPRNT, DTSAVE, TSTOP, TIM,
     1                PATM, TLIM, LOUT, LSAVE, LSENS, LIDAS, LRDAS,
     2                LSDAS, Z, ZP, DWORK, IDWORK, SDWORK, RPAR,
     3                IPAR, ATOL, RTOL, XMOL, KSYM, CSENK, INFO, ISEN)
      ELSEIF (ICASE .EQ. 6) THEN
         CALL SENS46 (RCONTV, NSYS, KK, II, DTPRNT, DTSAVE, TSTOP, TIM,
     1                PATM, TLIM, LOUT, LSAVE, LSENS, LIDAS, LRDAS,
     2                LSDAS, Z, ZP, DWORK, IDWORK, SDWORK, RPAR,
     3                IPAR, ATOL, RTOL, XMOL, KSYM, CSENK, INFO, ISEN)
      ENDIF
C
      RETURN
C
7000  FORMAT (/5X,'Restart calculation from previous solution.'/)
7100  FORMAT ('  Pressure (atm)  =', 1PE12.4,/,
     1        '  Temperature (K) =', 1PE12.4,/,
     2        '  Density (gm/cc) =', 1PE12.4)
7101  FORMAT (/,'  Mole Fractions:')
7102  FORMAT (1X,A10,2H =,1PE11.4)
7103  FORMAT (/,'  Initial Conditions:',/)
7104  FORMAT ('  Mass (gm)       =', 1PE12.4)
7111  FORMAT (/5X,'Pressure is constant.'/)
7112  FORMAT (/5X,'Volume is constant.'/)
7113  FORMAT (/5X,'Volume is a user-specified function of time.'/)
7114  FORMAT (/5X,'Temperature & pressure are constant.'/)
7115  FORMAT (/5X,'Temperature is a user-specified function of time.'/)
7116  FORMAT (/5X,'Temperature & volume are constant.'/)
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE RDSAVE (LREST, TRST, LOUT, KK, ICASE, ILINK, CSENK, 
     1                   IPAR, RPAR, XMOL, Z, TIM, T, P, IERR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION XMOL(KK), Z(*), IPAR(*), RPAR(*)
      CHARACTER*16 ILINK, CSENK, ICHR
      LOGICAL LSENSR, IERR
C
      IERR = .FALSE.
      TFIL = -1.0
      TIM  = TFIL
   18 CONTINUE
      ICHR = ' '
      READ (LREST, END=50, ERR=50) ICHR
      IF (ICHR .EQ. ILINK) THEN
C
C        skip CHEMKIN linking info
         DO 20 L = 1, 4
            READ (LREST, END=50, ERR=50)
20       CONTINUE
         GO TO 18
C
      ELSEIF (ICHR .EQ. CSENK ) THEN
C
        READ (LREST, END=50, ERR=50) LSENSR
        READ (LREST, END=50, ERR=50) NSYSR, KKR, IIR
C
        IF (KKR .NE. KK) THEN
            WRITE (LOUT, '(/3X, A, /9X, A, /9X, A, I6, /9X, A, I6)')
     2      'Stop! Number of species in restart file does not match',
     3      'the number from the CHEMKIN link file.',
     4      'No. species in restart file = ', KKR,
     5      'No. species in linking file = ', KK
            IERR = .TRUE.
            RETURN
        ENDIF
C
C       loop through restart file; store mass fractions in XMOL
        DO 30 J = 1, 10000
           READ (LREST, END=50, ERR=50) 
     1           TFIL, PFIL, TMPFIL, (XMOL(K),K=1,KK)
C          save latest solution
           TIM = TFIL
           P = PFIL
           IF (ICASE .LE. 3) THEN
              T    = TMPFIL
              Z(1) = T
              DO 25 K = 2, KK+1
                 Z(K) = XMOL(K-1)
   25         CONTINUE
           ELSE
              DO 26 K = 1, KK
                 Z(K) = XMOL(K)
   26         CONTINUE
           ENDIF
           
C          check for requested time
           IF (TRST .GE. 0.0) THEN
              IF (TIM .GE. TRST) GO TO 50
           ENDIF
C          skip sensitivities if present
           IF (LSENSR) READ (LREST, END=50, ERR=50)
   30    CONTINUE
      ENDIF
C
   50 CONTINUE
C
      IF (TIM .LT. 0.0) THEN
C        Error reading the binary file for a solution
         WRITE (LOUT,'(//5X,A//)') 
     1   'Stop! Error reading restart file.'
         IERR = .TRUE.
         RETURN
      ELSE
         IF (TIM .LT. TRST) THEN
C           Error or end-of-file, but can proceed with last time
            WRITE (LOUT, '(2(/5X,A),/5X,A,1PE13.6)')
     1      'DID NOT FIND RESTART TIME.',
     2      'WILL USE LAST TIME FOUND.',
     3      'Time (s) = ', TIM
         ENDIF
         IF (ICASE .LE. 3) THEN
            CALL CKYTX (Z(2), IPAR, RPAR, XMOL)
         ELSE
            CALL CKYTX (Z(1), IPAR, RPAR, XMOL)
         ENDIF
      ENDIF
C     end of SUBROUTINE SENKIN
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SETPAR (INFO, NSYS, TOLS, RTOL, ATOL, LSENS, ISEN, II)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION INFO(15), ISEN(5), RTOL(NSYS,*), ATOL(NSYS,*), TOLS(*)
      LOGICAL LSENS
C
C     set parameters for DASAC
      INFO(3) = 1
      INFO(2) = 1
C
C     tolerances on dependent variables
      DO 10 I = 1, NSYS
         RTOL(I,1) = TOLS(1)
         ATOL(I,1) = TOLS(2)
10    CONTINUE
C
C     tolerances on sensitivity coefficients
      IF (LSENS) THEN
         ISEN(1) = 1
         ISEN(5) = II
c         DO 20 I = 1, NSYS
c            DO 15 J = 2, II+1
          do 20 j = 2, II+1
             do 15 i = 1, nsys
               RTOL(I,J) = TOLS(3)
               ATOL(I,J) = TOLS(4)
15          CONTINUE
20       CONTINUE
      ENDIF
C
C     end of SUBROUTINE SETPAR
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SENS13 (RES, NSYS, KK, II, DTPRNT, DTSAVE, TSTOP,
     1                   TIM, PATM, TLIM, LOUT, LSAVE, LSENS, LIW,
     2                   LRW, LSW, Z, ZP, ELWRK, IELWRK, SWORK, RPAR,
     3                   IPAR, ATOL, RTOL, XMOL, KSYM, CSENK, INFO,
     4                   ISEN)
C
C  START PROLOGUE
C  This module directs the integration for cases 1-3, where temperature
C  is not known and the energy equation is included.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(NSYS,*), ZP(NSYS,*), ELWRK(*), IELWRK(*), SWORK(*),
     1          ISEN(*), RPAR(*), IPAR(*), INFO(*), RTOL(NSYS,*),
     2          ATOL(NSYS,*), XMOL(KK)
C
      EXTERNAL RES, DRES
      LOGICAL LSENS
      CHARACTER KSYM(KK)*(*), CSENK*16
C
      COMMON /RES1/ P
      DATA NOSAV /1/
C 
C     compute initial derivatives
C*****precision > single
C      CALL SDINIT 
C*****END precision > single
C*****precision > double
      CALL DDINIT 
C*****END precision > double
     1   (NSYS, ISEN(5), TIM, RPAR, IPAR, ZP, Z,
     2    IRES, ISEN(1), ISEN(4), RES, DRES)
C
C     initial outputs
      WRITE (LSAVE) CSENK
      WRITE (LSAVE) LSENS
      WRITE (LSAVE) NSYS, KK, II
      WRITE (LSAVE) TIM, P, (Z(I,1), I = 1, NSYS)
      IF (LSENS) WRITE (LSAVE) (( Z(I,J), I=1,NSYS), J = 2, II+1)
      WRITE (LOUT, '(//A/)') '  Time Integration:'
      CALL TEXT (IPAR, KK, KSYM, LOUT,
     1           P, PATM, RPAR, TIM, XMOL, Z(1,1), Z(2,1) )
C
C     integration loop
      TPRINT = DTPRNT + TIM
      TSAVE  = DTSAVE + TIM
      IFLG = 0
      TIGN = 0.0
250   CONTINUE
C
C     call the ODE solver
      TIMOLD = TIM
      TOLD   = Z(1,1)
C
C*****precision > single
C      CALL SDASAC 
C*****END precision > single
C*****precision > double
      CALL DDASAC 
C*****END precision > double
     1  (RES, NSYS, TIM, Z, ZP, TSTOP, INFO, ISEN,
     2   RTOL, ATOL, IDID, SWORK, LSW, ELWRK, LRW,
     3   IELWRK, LIW, RPAR, IPAR, JAC, DRES, DFDYP)
C
      IF (IDID .LT. 0) THEN
         WRITE (LOUT, '(1X,A,I3)') 'IDID =', IDID
         RETURN
      ENDIF
C
C     check for thermal runaway
      IF (IFLG .EQ. 0) THEN
         T = Z(1,1)
         IF (T .GE. TLIM) THEN
            DT = ELWRK(7)
            TIGN = TIMOLD + DT *(TLIM-TOLD)/(T-TOLD)
            IFLG = 1
         ENDIF
      ENDIF
C
C     save solution to binary file
      IF (DTSAVE .GT. 0.) THEN
         IF (TIM .LT. TSAVE) THEN
C            wait for a later solution to save
             GO TO 125
         ELSE
C            set next solution time
             TSAVE = TSAVE + DTSAVE
         ENDIF
      ENDIF
C
      WRITE (LSAVE) TIM, P, (Z(I,1), I = 1, NSYS)
      IF (LSENS) WRITE (LSAVE) ((Z(I,J), I = 1, NSYS), J = 2, II+1)
      NOSAV = NOSAV + 1
      TLASTS = TIM
C
  125 CONTINUE
C
C     print solution to text file
      IF (TIM .GE. TPRINT) THEN
         CALL TEXT (IPAR, KK, KSYM, LOUT,
     1              P, PATM, RPAR, TIM, XMOL, Z(1,1), Z(2,1) )
         TLASTP = TIM
         TPRINT = TPRINT + DTPRNT
      ENDIF
C
      IF (TIM .LT. TSTOP) GO TO 250
C
C     last print
      IF (TLASTP .NE. TIM) CALL TEXT (IPAR, KK, KSYM,
     1  LOUT, P, PATM, RPAR, TIM, XMOL, Z(1,1), Z(2,1) )

      IF (TLASTS .NE. TIM) THEN
        WRITE (LSAVE) TIM, P, (Z(I,1), I = 1, NSYS)
        IF (LSENS) WRITE (LSAVE) ((Z(I,J),I = 1,NSYS), J = 2,II+1)
        NOSAV = NOSAV + 1
      ENDIF
C
      WRITE (LOUT, 7040) TIGN
      WRITE (LOUT, 7045) TLIM
      WRITE (LOUT, 7050) NOSAV
C
C     FORMATS
 7040 FORMAT(//1X,'  Integration completed:',/2X,
     1' Ignition Time (sec) = ',1PE12.4)
 7045 FORMAT(2X,' Temp criteria (K) = ',1PE12.4)
 7050 FORMAT(/1X,'  Binary file has ',I6,' time datasets.')
C
C     end of SUBROUTINE SENS13
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SENS46 (RES, NSYS, KK, II, DTPRNT, DTSAVE, TSTOP,
     1                   TIM, PATM, TLIM, LOUT, LSAVE, LSENS, LIW,
     2                   LRW, LSW, Z, ZP, ELWRK, IELWRK, SWORK, RPAR,
     3                   IPAR, ATOL, RTOL, XMOL, KSYM, CSENK, INFO,
     4                   ISEN)
C
C  START PROLOGUE
C  This module directs the integration for cases 4-6, where temperature
C  is specified and so there is no energy equation. Here  NSYS = KK.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(NSYS,*), ZP(NSYS,*), ELWRK(*), IELWRK(*), SWORK(*),
     1          ISEN(*), RPAR(*), IPAR(*), INFO(*), RTOL(NSYS,*),
     2          ATOL(NSYS,*), XMOL(KK)
C
      EXTERNAL RES, DRES
      LOGICAL LSENS
      CHARACTER KSYM(KK)*(*), CSENK*16
      COMMON /RES1/ P
      COMMON /RES3/ T
      DATA NOSAV /1/
C
C     compute initial derivatives
C*****precision > single
C      CALL SDINIT
C*****END precision > single
C*****precision > double
      CALL DDINIT 
C*****END precision > double
     1  (NSYS, ISEN(5), TIM, RPAR, IPAR, ZP, Z,
     2   IRES, ISEN(1), ISEN(4), RES, DRES)
C
C     print
      WRITE (LSAVE) CSENK
      WRITE (LSAVE) LSENS
      WRITE (LSAVE) NSYS, KK, II
      WRITE (LSAVE) TIM, P, T, (Z(I,1), I = 1, NSYS)
      IF (LSENS) WRITE (LSAVE) (( Z(I,J), I=1,NSYS), J = 2, II+1)
      WRITE (LOUT, '(//A/)') '  Time Integration:'
      CALL TEXT (IPAR, KK, KSYM, LOUT,
     1           P, PATM, RPAR, TIM, XMOL, T, Z)
C
C     integration loop
      TPRINT = DTPRNT + TIM
      TSAVE  = DTSAVE + TIM
      IFLG = 0
250   CONTINUE
C
C     call the ODE solver
C*****precision > single
C      CALL SDASAC 
C*****END precision > single
C*****precision > double
      CALL DDASAC 
C*****END precision > double
     1  (RES, NSYS, TIM, Z, ZP, TSTOP, INFO, ISEN,
     2   RTOL, ATOL, IDID, SWORK, LSW, ELWRK, LRW,
     3   IELWRK, LIW, RPAR, IPAR, JAC, DRES, DFDYP)
C
      IF (IDID .LT. 0) THEN
         WRITE (LOUT, '(1X,A,I3)') 'IDID =', IDID
         RETURN
      ENDIF
C
C     save solution to binary file
      IF (DTSAVE .GT. 0.) THEN
        IF (TIM .LT. TSAVE) THEN
C          wait for a later solution to save
           GO TO 125
        ELSE
C          set next save time
           TSAVE = TSAVE + DTSAVE
        ENDIF
      ENDIF
C
      WRITE (LSAVE) TIM, P, T, (Z(I,1), I = 1, NSYS)
      IF (LSENS) WRITE (LSAVE) ((Z(I,J),I = 1,NSYS), J = 2,II+1)
      NOSAV = NOSAV + 1
      TLASTS = TIM
C
  125 CONTINUE
C     print solution to text file
      IF (TIM .GE. TPRINT) THEN
         CALL TEXT (IPAR, KK, KSYM, LOUT,
     1              P, PATM, RPAR, TIM, XMOL, T, Z)
         TLASTP = TIM
         TPRINT = TPRINT + DTPRNT
      ENDIF
C
      IF (TIM .LT. TSTOP) GO TO 250
C
C     last text print
      IF (TLASTP .NE. TIM) CALL TEXT (IPAR, KK, KSYM,
     1                LOUT, P, PATM, RPAR, TIM, XMOL, T, Z)
C
      IF (TLASTS .NE. TIM) THEN
        WRITE (LSAVE) TIM, P, T, (Z(I,1), I = 1, NSYS)
        IF (LSENS) WRITE (LSAVE) ((Z(I,J),I = 1,NSYS), J = 2,II+1)
        NOSAV = NOSAV + 1
      ENDIF
C
      WRITE (LOUT, 7050) NOSAV
C
C     FORMATS
 7040 FORMAT(//1X,'  Integration completed:',/2x,
     1' Ignition Time (sec) = ',1PE12.4)
 7050 FORMAT(/1X,'  Binary file has ',I6,' time datasets.')
C
C     end of SENS46
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RTEMP (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
C
C  START PROLOGUE
C  Residual of differential equations for case where pressure
C  and temperature are user specified functions of time.
C
C  Variables:
C    Z(K) = species mass fractions
C    RPAR = array of reaction pre-exponential constants
C    T    = temperature (Kelvin)
C    P    = pressure (dyne/cm2)
C    PA   = pressure (atm)
C    PATM = atmospheric pressure (dyne/cm2)
C    TIME = time (sec)

C  User supplys a subroutine for the temperature & pressure
C    SUBROUTINE TEMPT (TIME, T, PA)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*)
      COMMON /RES1/ P
      COMMON /RES3/ T
      COMMON /ATM/ PATM
C
      KK     = IPAR(1)
      IPRCK  = IPAR(2)
      IPRD   = IPAR(3)
      IPWT   = IPAR(6)
      IPWDOT = IPAR(7)
C      IPH    = IPAR(8)
      IPICK  = IPAR(9)
      LOUT   = IPAR(10)
      II     = IPAR(11)
C
C     modify CHEMKIN work array for pre-exponential
      DO 15 I = 1, II
        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C     user supplies temperature
      PA = P / PATM
      CALL TEMPT (TIME, T, PA)
C
C     call CHEMKIN subroutines
      CALL CKRHOY (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RHO)
      CALL CKWYP  (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RPAR(IPWDOT))
      IF (RHO .EQ. 0.0) THEN
         WRITE (LOUT, '(/1X,A)') 'Stop, zero density in RTEMP.'
         STOP
      ENDIF
      VOLSP = 1. / RHO
C
C     species equations
      DO 200 K = 1, KK
         K1 = K-1
         DELTA(K) = ZP(K) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
 200  CONTINUE
C
C     end of SUBROUTINE RTEMP
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RCONT (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
C
C  START PROLOGUE
C  Residual of differential equations for case where temperature,
C  and pressure are constant.  Density varies with the mean molecular
C  weight of the mixture.
C
C  Variables:
C    Z(K) = species mass fractions
C    T    = temperature (Kelvin) - constant in time
C    P    = pressure  (dyne/cm2) - constant in time
C    TIME = time (sec)
C    RHO  = density of mixture (gm/cm3)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*)
      COMMON /RES1/ P
      COMMON /RES3/ T
C
      KK     = IPAR(1)
      IPRCK  = IPAR(2)
      IPRD   = IPAR(3)
      IPWT   = IPAR(6)
      IPWDOT = IPAR(7)
C      IPH    = IPAR(8)
      IPICK  = IPAR(9)
C      LOUT   = IPAR(10)
      II     = IPAR(11)
C
C     modify CHEMKIN work array for pre-exponential
      DO 15 I = 1, II
        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C     call CHEMKINK subroutines
      CALL CKRHOY (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RHO)
      VOLSP = 1./ RHO
      CALL CKWYP (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RPAR(IPWDOT))
C
C     species equations
      DO 200 K = 1, KK
         K1 = K-1
         DELTA(K) = ZP(K) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
 200  CONTINUE
C
C     end of SUBROUTINE RCONT
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RCONTV (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
C
C  START PROLOGUE
C  Residual of differential equations for case where temperature,
C  and volume are constant.  Pressure varies with the mean molecular
C  weight of the mixture.
C
C  Variables:
C    Z(K) = species mass fractions
C    T    = temperature (Kelvin) - constant in time
C    P    = pressure  (dyne/cm2) 
C    TIME = time (sec)
C    RHO  = density of mixture (gm/cm3) - constant in time
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*)
      COMMON /RES3/ T
      COMMON /RES4/ RHO
C
      KK     = IPAR(1)
      IPRCK  = IPAR(2)
      IPRD   = IPAR(3)
      IPWT   = IPAR(6)
      IPWDOT = IPAR(7)
      IPICK  = IPAR(9)
      II     = IPAR(11)
C
C     modify CHEMKIN work array for pre-exponential
      DO 15 I = 1, II
        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C     call CHEMKIN subroutines
      CALL CKPY (RHO, T, Z, IPAR(IPICK), RPAR(IPRCK), P)
      VOLSP = 1./ RHO
      CALL CKWYP (P, T, Z, IPAR(IPICK), RPAR(IPRCK), RPAR(IPWDOT))
C
C     species equations
      DO 200 K = 1, KK
         K1 = K-1
         DELTA(K) = ZP(K) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
 200  CONTINUE
C
C     end of SUBROUTINE RCONTV
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RCONP (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
C
C  START PROLOGUE
C  Residual of differential equations for constant pressure case
C
C  Variables:
C    Z(1)   = temperature (Kelvin)
C    Z(K+1) = species mass fractions
C    P      = pressure (dyne/cm2) - constant in time
C    RHO    = density (gm/cm3)
C    RPAR   = array of reaction pre-exponential constants
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*)
      COMMON /RES1/ P
C
      KK     = IPAR(1)
      IPRCK  = IPAR(2)
      IPRD   = IPAR(3)
      IPWT   = IPAR(6)
      IPWDOT = IPAR(7)
      IPH    = IPAR(8)
      IPICK  = IPAR(9)
      LOUT   = IPAR(10)
      II     = IPAR(11)
C
C     modify CHEMKIN work array for pre-exponential
      DO 15 I = 1, II
        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C     call CHEMKIN subroutines
      CALL CKRHOY (P, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), RHO)
      CALL CKCPBS (Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), CPB)
      CALL CKWYP  (P, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK),
     1             RPAR(IPWDOT))
      CALL CKHMS  (Z(1), IPAR(IPICK), RPAR(IPRCK), RPAR(IPH))
      IF (RHO .EQ. 0.0) THEN
         WRITE (LOUT, '(/1X,A)') 'Stop, zero density in RCONP.'
         STOP
      ENDIF
      VOLSP = 1. / RHO
C
C     energy equation
      SUM = 0.
      DO 100 K = 0, KK-1
         SUM = SUM + RPAR(IPH+K) * RPAR(IPWDOT+K) * RPAR(IPWT+K)
 100  CONTINUE
      DELTA(1) = ZP(1) + VOLSP *SUM /CPB
C
C     species equations
      DO 200 K = 1, KK
         K1 = K-1
         DELTA(K+1) = ZP(K+1) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
 200  CONTINUE
C
C     end of SUBROUTINE RCONP
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RVOLT (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
C
C  START PROLOGUE
C  Residual of differential equations for case where volume is
C  a user-specified function of time.
C
C  Variables:
C    Z(1)   = temperature (Kelvin)
C    Z(K+1) = species mass fractions
C    P      = pressure (dyne/cm2)
C    RHO    = density (gm/cm3)
C    TOTMAS = mass (gm) - constant for closed system
C    TIME   = time (sec)
C
C  User supplys a subroutine for the volume:
C    SUBROUTINE VOLT (TIME, VOL, DVDT)
C     VOL    = volume of system
C     DVDT   = time derivative of system volume
C     VOLSP  = specific volume
C     VDOT   = time derivative of specific volume
C  Note: Consistent units for volume are (cm3), but the volume can
C        be considered to be normalized and therefore unitless.
C        The problem is solved using intensive variables, so the
C        solution is independent of extensive variables such as
C        the volume and total mass.  Subroutine VOLT is called at
C        time zero and a total mass is computed, but the solution
C        only depends on the density, so the total mass and total
C        volume are not important.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*)
      EXTERNAL VOLT
      COMMON /RES1/ P
      COMMON /RES2/ TOTMAS
C
      KK     = IPAR(1)
      IPRCK  = IPAR(2)
      IPRD   = IPAR(3)
      IPWT   = IPAR(6)
      IPWDOT = IPAR(7)
      IPU    = IPAR(8)
      IPICK  = IPAR(9)
      LOUT   = IPAR(10)
      II     = IPAR(11)
C
C     modify CHEMKIN work array for pre-exponential
      DO 15 I = 1, II
        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C     get volume as function of time
      CALL VOLT (TIME, VOL, DVDT)
      VOLSP = VOL  / TOTMAS
      VDOT  = DVDT / TOTMAS
      IF (VOLSP .EQ. 0.0) THEN
         WRITE (LOUT, '(/1X,A)') 'Stop, zero volume in RVOLT.'
         STOP
      ENDIF
      RHO = 1./ VOLSP
C
C     call CHEMKIN subroutines
      CALL CKCVBS (Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), CVB)
      CALL CKWYR  (RHO, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK),
     1             RPAR(IPWDOT))
      CALL CKUMS  (Z(1), IPAR(IPICK), RPAR(IPRCK), RPAR(IPU))
      CALL CKPY   (RHO, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), P)
C
C     energy equation
      SUM = 0.
      DO 100 K = 0, KK-1
         SUM = SUM + RPAR(IPU+K) * RPAR(IPWDOT+K) * RPAR(IPWT+K)
 100  CONTINUE
      DELTA(1) = ZP(1) + (P *VDOT + VOLSP *SUM) /CVB
C
C
C     species equations
      DO 200 K = 1, KK
         K1 = K-1
         DELTA(K+1) = ZP(K+1) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
 200  CONTINUE
C
C     end of SUBROUTINE RVOLT
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE RCONV (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
C
C  START PROLOGUE
C  Residual of differential equations when volume is constant
C
C  Variables:
C    Z(1)   = temperature (Kelvin)
C    Z(K+1) = species mass fractions
C    P      = pressure (dyne/cm2)
C    RHO    = density (gm/cm3) - constant in time
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Z(*), ZP(*), DELTA(*), RPAR(*), IPAR(*)
      COMMON /RES1/ P
      COMMON /RES4/ RHO
C
      KK     = IPAR(1)
      IPRCK  = IPAR(2)
      IPRD   = IPAR(3)
      IPWT   = IPAR(6)
      IPWDOT = IPAR(7)
      IPU    = IPAR(8)
      IPICK  = IPAR(9)
C      LOUT   = IPAR(10)
      II     = IPAR(11)
C
C     modify CHEMKIKN work array for pre-exponential
      DO 15 I = 1, II
        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
15    CONTINUE
C
C     constant density
      VOLSP = 1. / RHO
C
C     call CHEMKIN subroutines
      CALL CKCVBS (Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), CVB)
      CALL CKWYR  (RHO, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK),
     1             RPAR(IPWDOT))
      CALL CKUMS  (Z(1), IPAR(IPICK), RPAR(IPRCK), RPAR(IPU))
      CALL CKPY   (RHO, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), P)
C     
C     energy equation
      SUM = 0.
      DO 100 K = 0, KK-1
         SUM = SUM + RPAR(IPU+K) * RPAR(IPWDOT+K) * RPAR(IPWT+K)
 100  CONTINUE
      DELTA(1) = ZP(1) + VOLSP *SUM /CVB
C
C     species equations
      DO 200 K = 1, KK
         K1 = K-1
         DELTA(K+1) = ZP(K+1) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
 200  CONTINUE
C
C     end of SUBROUTINE RCONV
      RETURN
      END

C----------------------------------------------------------------------

      SUBROUTINE REDKEY (DTPRNT, DTSAVE, KK, KSYM, LIN, LOUT, P, RESTRT,
     1                   T, TLIM, TOLS, TRES, TSTOP, REAC, USEINP, TRST,
     1                   KERR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION REAC(KK), TOLS(*), VALUE(5)
C
      CHARACTER KEY*4, KSYM(*)*(*), LINE*80, CKCHUP*4
      LOGICAL KERR, IERR, RESTRT, USEINP
      EXTERNAL CKCHUP
      ICASE = 0
C
C     defaults
      KERR = .FALSE.
      DTSAVE = 0.
      TLIM = 0.0
      TRES = -1.
      RESTRT = .FALSE.
      TRST = -1.
      USEINP = .FALSE.
      DO 10 K = 1, KK
         REAC(K) = 0.
10    CONTINUE
      TOLS(1) = 1.E-8
      TOLS(2) = 1.E-20
      TOLS(3) = 1.E-5
      TOLS(4) = 1.E-5
C
      WRITE (LOUT, '(/A,/)') '     Keyword input:'
C
C     read next input line
90    CONTINUE
      LINE = ' '
      IERR = .FALSE.
      READ  (LIN,  '(A)',END=100) LINE
      WRITE (LOUT, '(A)') LINE
      CALL CKDTAB (LINE)
      KEY = ' '
      KEY = CKCHUP(LINE,4)
      LINE(1:4) = ' '
C
C     is this a keyword comment?
      IF (KEY(1:1) .EQ. '.' .OR. KEY(1:1) .EQ. '/' .OR.
     1    KEY(1:1) .EQ. '!') GO TO 90
C
C     tolerances (optional)
      IF (KEY .EQ. 'RTOL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TOLS(1), IERR)

      ELSEIF (KEY .EQ. 'ATOL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TOLS(2), IERR)

      ELSEIF (KEY .EQ. 'RTLS') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TOLS(3), IERR)

      ELSEIF (KEY .EQ. 'ATLS') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TOLS(4), IERR)
C
      ELSEIF (KEY .EQ. 'PRES') THEN
C        pressure
         CALL CKXNUM (LINE, 1, LOUT, NVAL, P, IERR)
C
      ELSEIF (KEY .EQ. 'TEMP') THEN
C        temperature
         CALL CKXNUM (LINE, 1, LOUT, NVAL, T, IERR)
C
      ELSEIF (KEY .EQ. 'TLIM') THEN
C        ignition temperature
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TLIM, IERR)
C
      ELSEIF (KEY .EQ. 'REAC') THEN
C        species
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (.NOT. IERR) REAC(KSPEC) = VALUE(1)
C
      ELSEIF (KEY .EQ. 'TIME') THEN
C        integration time
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TSTOP, IERR)
C
      ELSEIF (KEY .EQ. 'DELT') THEN
C        text output time interval
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DTPRNT, IERR)
C
      ELSEIF (KEY .EQ. 'DTSV') THEN
C        binary output time interval
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DTSAVE, IERR)
C
      ELSEIF (KEY .EQ. 'REST') THEN
C        restart job
         RESTRT = .TRUE.
C
      ELSEIF (KEY .EQ. 'USET') THEN
C        use input T,P for restart
         USEINP = .TRUE.
C
      ELSEIF (KEY .EQ. 'TRES') THEN
C        restart time for new job
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TRES, IERR)
C
      ELSEIF (KEY .EQ. 'TRST') THEN
C        restart time from file
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TRST, IERR)
C
      ELSEIF (KEY .EQ. 'END ') THEN
C        last line
         GO TO 100

      ELSE
C       keyword is bogus
C
        WRITE (LOUT, '(/2X,A,A,A)')
     1  ' Warning, keyword ', KEY, ' was not found.'
      ENDIF
C
      IF (IERR) THEN
         KERR =.TRUE.
         WRITE (LOUT, '(2X,A,A,A,/)')
     1    ' Error reading data for keyword ', KEY,'.'
      ENDIF
      GO TO 90
C
  100 CONTINUE
C
C     END OF INPUT
C
C     set limit temperature if not input
      IF (TLIM.EQ.0.) TLIM = T + 400.
C
      IF (RESTRT) RETURN
C
C     check for necessary input
      IF (P .LE. 0) THEN
         WRITE (LOUT, *) ' Error, "PRES" not properly specified.'
         KERR = .TRUE.
      ENDIF
      IF (T .LE. 0) THEN
         WRITE (LOUT, *) ' Error, "TEMP" not properly specified.'
         KERR = .TRUE.
      ENDIF
      IF (TSTOP .LE. 0) THEN
         WRITE (LOUT, *) ' Error, "TIME" not properly specified.'
         KERR = .TRUE.
      ENDIF
      IF (DTPRNT .LE. 0) THEN
         WRITE (LOUT, *) ' Error, "DELT" not properly specified.'
         KERR = .TRUE.
      ENDIF
C
C     check the reactant sum
      SUMR = 0.
      DO 6100 K = 1, KK
         SUMR = SUMR + REAC(K)
6100  CONTINUE
      IF (SUMR .LE. 0) THEN
         WRITE (LOUT, *) ' Error, "REAC" not properly specified.'
         KERR = .TRUE.
      ENDIF
C
      IF (KERR) RETURN
C
C     normalize reactant mole fractions
      DO 6200 K = 1, KK
         REAC(K) = REAC(K) / SUMR
6200  CONTINUE
C
C     end of SUBROUTINE REDKEY
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE REDSEN (ICASE, LIN, LOUT, LSENS, KERR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*4 ISTR, KEY, CKCHUP
      EXTERNAL CKCHUP
      LOGICAL LSENS, KERR
      LSENS = .FALSE.
C
      KERR = .FALSE.
C     issue a prompt
20    WRITE (LOUT, '(//5X,A)') ' Enter keyword: '
C
C     read next input line
50    READ  (LIN,  '(A)')     ISTR
      WRITE (LOUT, '(10X,A4)')ISTR
      CALL CKDTAB (ISTR)
      KEY = CKCHUP(ISTR,4)
C
C     is this a keyword comment?
      IF (KEY(1:1) .EQ. '.' .OR. KEY(1:1) .EQ. '/' .OR.
     1    KEY(1:1) .EQ. '!') GO TO 50
C
C     sensitivity keyword
      IF (KEY .EQ. 'SENS') THEN
         WRITE (LOUT, 9000)
         LSENS = .TRUE.
         GO TO 20
C
C      problem choice keyword
       ELSEIF (KEY .EQ. 'CONP') THEN
         ICASE = 1
       ELSEIF (KEY .EQ. 'CONV') THEN
         ICASE = 2
       ELSEIF (KEY .EQ. 'VTIM') THEN
         ICASE = 3
       ELSEIF (KEY .EQ. 'CONT') THEN
         ICASE = 4
       ELSEIF (KEY .EQ. 'TTIM') THEN
         ICASE = 5
       ELSEIF (KEY .EQ. 'CNTV') THEN
         ICASE = 6
       ELSE
          WRITE (LOUT, 9002) 'Warning, keyword ',KEY, ' not found.'
       ENDIF
C
C     was sensitivity found?
      IF (.NOT. LSENS) WRITE (LOUT, 9001)
C
C     check for problem selection
      IF (ICASE.GE.1 .AND. ICASE.LE.6) RETURN
C
      WRITE (LOUT, '(/10X, A, /5X, A, /5X, A)')
     1'Fatal error.',
     2'The problem choice keyword must appear on one of the first',
     3'two input lines (CONP, CONV, CONT, VTIM, TTIM, or CNTV).'
      KERR = .TRUE.
C
9000  FORMAT (/5X, 'Sensitivity analysis will be performed.')
9001  FORMAT (/5X, 'Sensitivity analysis will not be performed.')
9002  FORMAT (5X,A)
9003  FORMAT (5X,A,A,A)
C
C     end of SUBROUTINE REDSEN
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE TEXT (IPAR, KK, KSYM, LOUT,
     1                   P, PATM, RPAR, TIM, XMOL, T, Y)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION Y(KK), RPAR(*), IPAR(*), XMOL(KK)
      CHARACTER*(*) KSYM(KK)
      DATA NPL /3/
C
C     compute mole fractions
      IPRCK  = IPAR(2)
      IPICK  = IPAR(9)
      CALL CKYTX (Y, IPAR(IPICK), RPAR(IPRCK), XMOL)
C
      WRITE (LOUT, 7700) TIM, P/PATM, T
      DO 30 K = 1, KK, NPL
        WRITE (LOUT, 7710) (KSYM(I),XMOL(I),I=K,MIN(K+NPL-1,KK))
30    CONTINUE
C
 7700 FORMAT(/,' t(sec)=',1PE11.4,'   P(atm)=',1PE11.4,
     1'   T(K)=',1PE11.4)
 7710 FORMAT(3(2X,A10,'=',1PE9.2))
C
C     end of SUBROUTINE TEXT
      RETURN
      END
C
C----------------------------------------------------------
C
      SUBROUTINE DRES (T, Y, YP, D, IM, RP, IP)
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
      DIMENSION Y(*), YP(*), D(*), RP(*), IP(*)
C     dummy routine to avoid missing external
      RETURN
      END
