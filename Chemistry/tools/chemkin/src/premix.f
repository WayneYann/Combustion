C     CVS $Revision: 1.28 $ created $Date: 2008-07-28 01:08:00 $
      SUBROUTINE PMABS
C
C///////////////////////////////////////////////////////////////////////
C
C     ONE DIMENSIONAL LAMINAR PREMIXED FLAME CODE
C
C     WRITTEN BY:
C         ROBERT J. KEE
C         COMPUTATIONAL MECHANICS DIVISION
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94550
C         (510) 294-3272
C
C     MODIFIED BY:
C         JOSEPH F. GRCAR
C         COMPUTATIONAL MECHANICS DIVISION
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94551-0969
C         (510) 294-2662
C         na.grcar@na-net.ornl.gov
C         sepp@snll-arparw.ornl.gov
C
C     AND BY:
C         FRAN RUPLEY
C         COMPUTATIONAL MECHANICS DIVISION
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94551-0969
C         (510) 294-3657
C
C/////////////////////////////////////////////////////////////////////
C
C     Revision: 3.15, Date: 1998/3/3, Author: E. Meeks.
C     1)  Fixed action#110: removed unused parameter ZERO from FLDRIV.
C     Revision: 3.14, Date: 1998/2/23, Author: E. Meeks
C     1)  Fixed bug#152:  Moved FUNCTION AREA to prdriv.f.
C     Revision: 3.13, Date: 1998/2/8, Author: E. Meeks
C     1)  Fixed bug#150:  Changed logic in line 2946 in
C         subroutine RDKEY to not require TEMP keywords on a restart.
C     Revision: 3.12, Date: 1997/10/29, Author: E. Meeks
C     1)  Fixed bug #098: Set ABSOL=RELAT=SQRT(D1MACH(4)).
C     2)  Removed UROUND function.
C     3)  Put double-precision change blocks around above, and
C         add single-precision change block using R1MACH(4).
C     Revision: 3.11, Date: 1997/07/25, Author: J. Grcar
C     1)  Fixed bug #064:  Correct IERR = IERR .OR. NP+1 .LE. JMAX
C         to be IERR = IERR .OR. .NOT. (NP+1 .LE. JMAX) on line 2576.
C     2)  On line 2942, TWSETL argument should be ERROR, not KERR;
C         add line following call:  KERR = KERR .OR. ERROR.
C     Revision: 3.10, Date: 1997/07/23, Author: J. Grcar
C     1)  Fix bug 43: QTDEC should be TDEC near line 2537.
C     Revision: 3.9, Date: 1997/06/10, Author: J. Grcar
C     1)  Fixed bug #022.  Write "continuing to new problem" on many
C         separate lines rather than many times on one line.
C     2)  Edited "CVS $Revision" line to fit in 72 columns.
C     3)  Updated output version number and creation date.
C   V3.8 97/03/01
C   1. new main "driver" program to set up arrays and open files,
C      to satisfy Reaction Design requirement of providing object
C      files instead of source code.
C   V3.7, 96/09/04
C   1. dimension COND(JMAX) added to FLDRIV
C   V. 3.6, 96/05/24
C   1. initial sccs version
C   VERSION 3.5 (4/17/95 F. Rupley)
C   1.  initial CHEMKIN-III version
C   2.  implement twopnt v.3.22
C   3.  clean up some loops by setting JM1, JP1, etc.
C   4.  set up START/END PROLOGUE sections
C   5.  converted some subroutines to function
C
C   CHANGES FOR V.3.4 (11/95 F. Rupley)
C     1. COMMON blocks added
C     2. Debug REGRID and TWOPNT upgrades
C   CHANGES FOR V.3.3 (11/24/92 F. Rupley)
C     1. Regridding must be done BEFORE the CALL RESTART, for purposes
C        of setting JFIXT, X(JFIXT), etc.
C     2. In START, some variables were not being initialized when a
C        point is added.
C   CHANGES FOR V.3.2 (11/12/92 F. Rupley)
C     1. Bugfix add SPOS to PRETWO arguments
C     2. Update for CKLIB.40: do not need linking file unit numbers
C        for CKSAVE, MCSAVE
C     3. Add regridding scheme; pointer NREG for an array of length
C        PMAX, real variables JJREGD, PCTACP and RATGTC, logical LREGRD
C     4. Some character variables (KSYM, CNTRL) cannot be (*)*(*) for
C        SGI compiler - need actual lengths in declarations
C   CHANGES FOR V.3.0 (6/4/92 F. Rupley)
C     1. SUBROUTINE PRETWO for TWOPNT input/output
C     2. SUBROUTINE PRTINIT for TWOPNT array initialization
C     3. Remove options SCAL, SCLT, ESEN, DSEN, etc.
C   CHANGES FOR V.2.9 (4/28/92 F. Rupley)
C     1. SUBROUTINE PRREAD reads solution for restart
C     2. Remove variable and keyword UFAC
C   CHANGES FOR VERSION 2.8 (4/03/92 by J. F. Grcar)
C     1. Change FUN to avoid rate evaluation at points where nothing
C        is perturbed during Jacobian preparation. This reduces Cray YMP
C        computing time by about 25 percent for the hydrogen flame
C        example (10.04 seconds to 7.59).
C     2. Change NMAX to PMAX in all subroutines.
C     3. Change subroutine POINT to POINTR to avoid naming conflicts.
C     4. Changed names of some change blocks.
C   CHANGES FOR VERSION 2.7 (3/23/92 by J. F. Grcar)
C     1. Change FUN to use CKLIB version 3.7 routines CKKFRT and CKWYPK.
C        This reduces Cray YMP computing time by about 20 percent for
C        the hydrogen flame example (12.43 YMP seconds to 10.04).
C   CHANGES FOR VERSION 2.6 (2/17/92 by J. F. Grcar)
C     1. Install Twopnt version 3.08 with changes to subroutines
C        FLDRIV, POINT, PREMIX, RDKEY, REASEN and REHSEN.
C     2. In some subroutines, variables KK, COMPS and NMAX have been
C        renamed POINTS, COMPS and PMAX to conform to TWOPNT.
C     3. In FLDRIV, modify a printing message to include the time
C        of computing sensitivity coefficients.
C     4. In RDKEY, omit trailing blanks when echoing keyword input.
C     5. In PRINT, reposition column headers and add a change block
C        for 80 column output.
C     6. Remove subroutine JACOB.
C   CHANGES FOR VERSION 2.5 (7/16/91 F. Rupley per R. Kee)
C     1. Correction for velocity calculation in Subroutine PRINT
C   CHANGES FOR VERSION 2.4 (5/9/91 F. Rupley per R. Kee)
C     1. Subroutine REHSEN and keyword "HSEN" for sensitivity to
C        heats of formation.
C   CHANGES FOR VERSION 2.3 (4/29/91 F. Rupley per Bob Kee)
C     1. Correction in MTNRPR for TDIF option.
C   CHANGES FOR VERSION 2.1 (1/15/91 F. Rupley per Bob Kee)
C     1. Set default value of LVCOR to .TRUE. in Subroutine RDKEY
C   CHANGES FOR VERSION 2.0
C     1. Do not use linking file information stored in a restart
C        solution, as mechanism may have changed.
C     2. Call list for TWOPNT requires additional input; optional
C        use of new keywords reset default values:
C        'ISTP' n - sets NINIT initial time steps before newton
C                   (default is 0)
C        'IRET' n - set retirement age IRETIR of old time step
C                   (default 50)
C        'NJAC' n - set retirement age NJAC of Jacobian during
C                   steady state newton (default 20)
C        'TJAC' n - set retirement age ITJAC of Jacobian during
C                   time stepping (default 20)
C        'DTMX' x - set maximum time step DTMAX (default 1.0E-4)
C     3. Keyword SPOS sets a minimum value of species solutions.
C     4. Keyword NOFT allows skipping of the fixed temperature
C        problem.
C     5. Keywork GFAC allows scaling of production rates.
C   CHANGES FOR VERSION 1.9
C     1. Restore file has additioal work array records, as in solution
C        file
C   CHANGES FOR VERSION 1.8
C     1. Refine some DO loops
C   CHANGES FOR VERSION 1.7
C     1. Multicomponent version
C     2. Solution has additional work array records
C   CHANGES FOR VERSION 1.6
C     1. Modify POINT to use new subroutines CKLEN and TPLEN instead
C        of reading linking file.
C   CHANGES FOR VERSION 1.4
C     1.  CONVERSION TO MULTICOMPONENT GAS TRANSPORT
C     2.  INCLUDE SPECIES V CORRECTION AS OPTION
C   CHANGES FOR VERSION 1.3
C     1.  LINKING FILE HAS ADDITIONAL VARIABLES KERR, MAXTB
C   CHANGES FOR VERSION 1.2
C     1.  ALLOW MIXED CASE SPECIES NAMES TO REMAIN MIXED
C   CHANGES FOR VERSION 1.1
C     1.  READ LENICK, LENCK, LENCCK FROM LINKING FILE
C         (SUBROUTINE POINT) INSTEAD OF CALCULATING
C     2.  ALLOW UPPER/LOWER CASE INPUT
C   CHANGES FOR VERSION 1.0
C     1.  CHANGED REAL*8 TO DOUBLE PRECISION
C
C///////////////////////////////////////////////////////////////////////
C
C     end of SUBROUTINE PMABS
      END
C
      SUBROUTINE FLDRIV (LIN, LOUT, LREST, LSAVE, LRCRVR, LINKCK,
     1                   LINKMC, JMAX, RCKWRK, RMCWRK, EPS, WT, REAC,
     2                   SCRTCH, X, COND, REG, TGIVEN, XGIVEN, D, DKJ,
     3                   TDR, YV, ABOVE, BELOW, BUFFER, S, SN, F, FN,
     4                   DS, A6, A, ICKWRK, IMCWRK, KSYM, CCKWRK, KR,
     5                   KI, KP, IPIVOT, ACTIVE, MARK, NAME, ITWWRK,
     6                   RTWWRK, SSAVE, RKFT, RKRT, RSAVE)
C
C  START PROLOGUE
C
C  LIN        - integer scalar, formatted input file unit number
C  LOUT       - integer scalar, formatted output file unit number
C  LREST      - integer scalar, binary input restart file unit number
C  LSAVE      - integer scalar, binary output solution file unit number
C  LRCRVR     - integer scalar, binary output scratch file unit number
C  LINKCK     - integer scalar, CHEMKIN input linkfile unit number
C  LINKMC     - integer scalar, TRANSPORT input linkfile unit number
C  JMAX       - integer scalar, maximum number of gridpoints
C  RCKWRK(*)  - real array, CHEMKIN workspace
C  RMCWRK(*)  - real array, TRANSPORT workspace
C  EPS(*)     - real array, mass flux fractions at the burner
C  WT(*)      - real array, species molecular weights (gm/mole)
C  REAC(*)    - real array, reactant input mole/mass fractions
C  SCRTCH(*,*)- real matrix, workspace
C  X(*)       - real array, mesh point locations (cm)
C  COND(*)    - real array, thermal conductivities at mesh midpoints
C               (erg/cm-K-sec);
C               if LVARMC=.TRUE., these are computed each time function
C               is called, otherwise the stored values are used
C  REG(*)     - real array, mesh locations of a regrid operation (cm)
C  TGIVEN(*)  - real array, temperatures at given distances (K)
C  XGIVEN(*)  - real array, distances from the burner for TGIVEN
C  D(*,*)     - real matrix, species diff coeff'nts at mesh midpoints
C               (cm^2/sec);
C               if LVARMC=.TRUE., these are computed each time function
C               is called, otherwise the stored values are used
C  DKJ(*,*,*) - real 3-dimensional array, multicomponent diff coeff'nts
C               at mesh midpoints (cm^2/sec)
C  TDR(*,*)   - real matrix, species thermal diff ratios at mesh
C               midpoints;
C               if LVARMC=.TRUE., these are computed each time fucntion
C               is called, otherwise the stored values are used.
C  YV(*,*)    - real matrix, mass fractions time diffusion velocities at
C               mesh midpoints;
C               YV(K,J) is flux of species K between J and J+1.
C  ABOVE(*)   - real array, upper bounds for solving dependent variables
C  BELOW(*)   - real array, lower bounds for solving dependent variables
C  BUFFER(*,*)-
C  S(*,*)     - real matrix, dependent variables;
C                            S(NT,J) = T(J), temperature
C                            S(NYS+K,J) = Y(K,J), mass fractions
C                            S(NM,J) = FLRT(J), flow rate
C  SN(*,*)    - real matrix, dependent variables at previous timestep
C  F(*,*)     - real matrix, governing equations residuals of S(*,*);
C               F(NT,J) is the energy equation residual,
C               F(NM,J) is the mass equation residual,
C               F(NYS+K,J) is the species equation for species K
C  FN(*,*)    -
C  DS(*)      -
C  A6(*)      -
C  A(*)       -
C  ICKWRK(*)  - integer array, CHEMKIN workspace
C  IMCWRK(*)  - real array, TRANSPORT workspace
C  KSYM(*)    - character-array, species names
C  CCKWRK(*)  - character-array, CHEMKIN workspace
C  KR(*)      - integer array, reactant species indices
C  KI(*)      - integer array, intermediate species indices
C  KP(*)      - integer array, product species indices
C  IPIVOT(*)
C  ACTIVE(*)
C  MARK(*)
C  NAME(*)    - character-string array
C  ITWWRK(*)  - integer array, TWOPNT workspace
C  RTWWRK(*)  - real array, TWOPNT workspace
C  SSAVE(*)   - real array, temperature values at which stored reaction
C               rates were evaluated.
C  RKFT(*)    - real array, forward reaction rates at SSAVE.
C  RKRT(*)    - real array, reverse reaction rates at SSAVE.
C  RSAVE(*)   - real matrix, for ICASE=2, save species production rates,
C               for ICASE=3, use RSAVE for species production rates
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
C
      CHARACTER CCKWRK(LENCCK)*(*), KSYM(KK)*(*), NAME(NATJ)*(*),
     1          ISOLUT*16, ISENSI*16, IHSENS*16, ID*9
C
      LOGICAL ACTIVE(NATJ), IERR, LASEN, LBURNR, LCNTUE, KERR,
     1        LENRGY, LHSEN, LMOLE, LMULTI, LRSTRT, LTDIF, LTIME,
     2        LUMESH, LUSTGV, LVARMC, LVCOR, MARK(JMAX), RSTCNT, LREGRD,
     3        linflow
C
C     Integer arrays
      DIMENSION ICKWRK(LENICK), IMCWRK(LENIMC), ITWWRK(LENITW), KI(KK),
     1          KP(KK), KR(KK), IPIVOT(NATJ,JMAX)
C     Real arrays
      DIMENSION A(IASIZE), A6(NTR), ABOVE(NATJ), BELOW(NATJ),
     1          BUFFER(NATJ,JMAX), COND(JMAX), REG(JMAX), D(KK,JMAX),
     2          DKJ(KK,KK,JMAX), DS(JMAX), F(NATJ,JMAX), FN(NATJ,JMAX),
     3          RCKWRK(LENRCK), REAC(KK), RKFT(II,JMAX), RKRT(II,JMAX),
     4          RMCWRK(LENRMC), RTWWRK(LENRTW), S(NATJ,JMAX),
     4          SCRTCH(KK, 6), EPS(KK), SN(NATJ,JMAX), TDR(KK,JMAX),
     5          TGIVEN(JMAX), SSAVE(NATJ,JMAX), WT(KK), X(JMAX),
     6          XGIVEN(JMAX), YV(KK,JMAX), RSAVE(KK,JMAX)
C
      PARAMETER (ID='FLDRIV:  ', ISOLUT='SOLUTION',
     1           ISENSI='SENSITIVITY', IHSENS='HSENSITIITY')
      DATA LCNTUE /.FALSE./
C
C///  INITIALIZATION.
C
      KERR = .FALSE.
      ONE = 1.0
C*****precision > double
      ABSOL = SQRT(D1MACH(4))
C*****END precision > double
C*****precision > single
C      ABSOL = SQRT(R1MACH(4))
C*****END precision > single
      RELAT = ABSOL
      RSTCNT = .FALSE.
      NT = 1
      NYS = 1
      NY = 2
      NM = KK + 2
C
      CALL CKSYMS (CCKWRK, LOUT, KSYM, IERR)
      CALL CKWT (ICKWRK, RCKWRK, WT)
      CALL CKRP (ICKWRK, RCKWRK, RU, RUC, PATM)
C
C///////////////////////////////////////////////////////////////////////
C
C     (1) TOP OF THE LOOP OVER PROBLEM CONTINUATIONS.
C
C///////////////////////////////////////////////////////////////////////
C
C///  TOP OF THE LOOP.
C
1010  CONTINUE
C      CALL TWTIME(XSTART)
C
C///  READ THE KEYWORDS.
C
      CALL RDKEY (JMAX, LIN, LOUT, KSYM, LBURNR, LMOLE, LUSTGV, LENRGY,
     1            LMULTI, LVCOR, LTDIF, LUMESH, LRSTRT, LCNTUE, MFILE,
     2            LASEN, LHSEN, NTOT, X, REAC, SCRTCH(1,2),
     3            SCRTCH(1, 3), KR, KI, KP, XGIVEN, TGIVEN, N1CALL,
     4            LREGRD, PCTADP, RATGTC, KERR, linflow, tinflow)
      IF (KERR) RETURN
C     RDKEY sets JJ=6
C
      IF (LRSTRT) THEN
C        Read the restart data file for starting profiles.
C
         IF (.NOT. RSTCNT) THEN
C           read a restart file solution (sets JJ to old solution)
            CALL PRREST (JMAX, LREST, LOUT, MFILE, X, S, KERR)
            IF (KERR) RETURN
         ENDIF
C
c         IF (LREGRD .AND. (RSTCNT.OR.LRSTRT) .AND. (JJ.GT.JJREGD))
         IF (LREGRD .AND. (RSTCNT.OR.LRSTRT) )
C        (sets JJ=JJREGD)
     1   CALL PRREGD (SN, S, REG, X, COND, PCTADP, RATGTC, XGIVEN,
     3                TGIVEN, LOUT, LBURNR, SCRTCH(1,4))
C        (adds to JJ if required)
c         CALL RESTRT (JMAX, LOUT, LMOLE, LBURNR, LUSTGV, REAC,
C        JFG: this appears to be a bug. shouldn't the restart file
C        always be written in the same units?
         CALL RESTRT (JMAX, LOUT, .true., LBURNR, LUSTGV, REAC,
     1                XGIVEN, TGIVEN, ICKWRK, RCKWRK, EPS, X, S,
     2                KERR, linflow, tinflow)
         IF (KERR) RETURN
C
      ELSE
C        Form the starting profiles;
C        in the following call, COND is length KK scratch space.
C
         CALL START (LOUT, LMOLE, LUMESH, LBURNR, REAC, SCRTCH(1,2),
     1               SCRTCH(1, 3), KR, KI, KP, XGIVEN, TGIVEN,
     2               ICKWRK, RCKWRK, SCRTCH(1, 4), COND, EPS, X, S,
     3               JMAX, KERR)
         IF (KERR) RETURN
      ENDIF
C
C///  CALL TWOPNT
C
      CALL PRETWO (LBURNR, LENRGY, LMULTI, LVCOR, LTDIF, LVARMC, LTIME,
     1             LMOLE, WT, EPS, XGIVEN, TGIVEN, X, S, SN, BUFFER,
     2             SCRTCH, YV, COND, D, DKJ, TDR, ICKWRK, RCKWRK,
     3             CCKWRK, IMCWRK, RMCWRK, ITWWRK, RTWWRK, F, FN,
     4             SSAVE, RKFT, RKRT, RSAVE, LOUT, LRCRVR, A, CONDIT,
     5             IPIVOT, NAME, KSYM, ABOVE, BELOW, JMAX, MARK,
     6             ACTIVE, N1CALL, KERR)
      IF (KERR) RETURN
C
C     WRITE TO LSAVE WHEN SOLUTION IS COMPLETE
C
      WRITE (LSAVE) ISOLUT
      WRITE (LSAVE) NATJ, JJ, P, S(NM, 1)
      WRITE (LSAVE) (X(J), J = 1, JJ)
      WRITE (LSAVE) ((S(N, J), N = 1, NATJ), J = 1, JJ)
C
C///////////////////////////////////////////////////////////////////////
C
C     (3) PERFORM SENSITIVITY ANALYSES.
C
C///////////////////////////////////////////////////////////////////////
C
      IF (LASEN .OR. LHSEN) THEN
C
         IF (LASEN)
C        Sensitivity to reaction rates
     1   CALL PRSENS (ISENSI, LBURNR, LENRGY, LMULTI, LTDIF, LMOLE,
     2                LSAVE, LRCRVR, LOUT, LVARMC, LTIME, WT, EPS,
     3                XGIVEN, TGIVEN, X, SN, S, SCRTCH, YV, COND, D,
     4                DKJ, TDR, ICKWRK, RCKWRK, IMCWRK, RMCWRK, F, FN,
     5                DS, A, A6, IPIVOT, BUFFER, SSAVE, RKFT, RKRT,
     6                RSAVE, KERR)

         IF (LHSEN)
C        Sensitivity of enthalpies
     1   CALL PRSENS (IHSENS, LBURNR, LENRGY, LMULTI, LTDIF, LMOLE,
     2                LSAVE, LRCRVR, LOUT, LVARMC, LTIME, WT, EPS,
     3                XGIVEN, TGIVEN, X, SN, S, SCRTCH, YV, COND, D,
     4                DKJ, TDR, ICKWRK, RCKWRK, IMCWRK, RMCWRK, F, FN,
     5                DS, A, A6, IPIVOT, BUFFER, SSAVE, RKFT, RKRT,
     6                RSAVE, KERR)
         IF (KERR) RETURN
C
      ENDIF
C
C      CALL TWTIME (xSTOP)
C      telaps = xstop - xstart
C      if (telaps .gt. 60) then
C      WRITE (STRING, '(F10.3, A)') telaps/60., ' MINUTES FOR SOLUTION.'
C      else
C      write (string, '(F10.3, a)') telaps, ' SECONDS FOR SOLUTION.'
C      endif
C      CALL TWSQEZ (LENGTH, STRING)
C      WRITE (LOUT, 10001) ID, STRING (1 : LENGTH)
10001 FORMAT (/1X, A9, A)
C
      IF (.NOT. LCNTUE) RETURN
C
      WRITE (LOUT, '( /////)')
      DO 1020 L = 1, 5
         WRITE (LOUT, *)
     1      ' //////////// CONTINUING TO NEW PROBLEM ////////////'
1020  CONTINUE
      WRITE (LOUT, '( /////)')
C
      RSTCNT = .TRUE.
      LRSTRT = .TRUE.
      GO TO 1010
C
      END

      SUBROUTINE FUN (LBURNR, LENRGY, LMULTI, LVCOR, LTDIF, LVARMC,
     1                LTIME, WT, EPS, XGIVEN, TGIVEN, X, SN, S, YAV,
     2                YV, WDOT, CP, H, COND, D, DKJ, TDR, ICKWRK,
     3                RCKWRK, IMCWRK, RMCWRK, F, XAV, SSAVE, RKFT,
     4                RKRT, ICASE, RSAVE, MARCDBG)
C
C  START PROLOGUE
C
C  INPUT
C  LBURNR     - logical, T  burner stabilized flame,
C                        F, freely propagating adiabatic flame
C  LENRGY     - logical, T  solve energy equation,
C                        F, use fixed temperature profile
C  LMULTI     - logical, T, use multicomponent formulas,
C                        F, use mixture-averaged formulas
C  LVCOR      - logical, T, use correction velocity formulism
C                        F, use 'trace' approximation, lumping
C                           all transport errors into final species
C  LTDIF      - logical, T, evaluate thermal diffusion ratios as
C                           well as diffusion coefficients
C  LVARMC     - logical, T, compute new transport properties
C                        F, use previously stored values
C  LTIME      - logical, T, time step of DT will be added to residual
C  WT(*)      - real array, species molecular weights (gm/mole)
C  EPS(*)     - real array, mass flux fractions at the burner
C  XGIVEN(*)  - real array, distances from the burner at which
C                           temperatures are specified (cm)
C  TGIVEN(*)  - real array, temperatures at given distances (K)
C  X(*)       - real array, mesh point locations (cm)
C  SN(*,*)    - real matrix, dependent variables at previous timestep
C  S(*,*)     - real matrix, dependent variables;
C                            S(NT,J) = T(J), temperature
C                            S(NYS+K,J) = Y(K,J), mass fractions
C                            S(NM,J) = FLRT(J), flow rate
C  YAV(*)    - real array, mass fractions at mesh midpoints;
C              YAV(K) is between J and J+1
C  YV(*,*)   - real matrix, mass fractions times diffusion velocities at
C              mesh midpoints;
C              YV(K,J) is flux of species K between J and J+1.
C  WDOT(*)   - real array, species chemical production rates
C              (moles/cm^3-sec)
C  CP(*)     - real array, species specific heats (ergs/gm-K)
C  H(*)      - real array, species enthalpies (ergs/mole)
C  COND(*)   - real array, thermal conductivities at mesh midpoints
C              (erg/cm-K-sec);
C              if LVARMC=.TRUE., these are computed each time function
C              is called, otherwise the stored values are used
C  D(*,*)    - real matrix, species diff coeff'nts at mesh midpoints
C              (cm^2/sec);
C              if LVARMC=.TRUE., these are computed each time function
C              is called, otherwise the stored values are used.
C  DKJ(*,*,*)- real 3-dimensional array, multicomponent diff coeff'nts
C              at mesh midpoints (cm^2/sec)
C  TDR(*,*)  - real matrix, species thermal diff ratios at mesh
C              midpoints;
C              if LVARMC=.TRUE., these are computed each time function
C              is called, otherwise the stored values are used.
C  ICKWRK(*) - integer array, CHEMKIN workspace
C  RCKWRK(*) - real array, CHEMKIN workspace
C  IMCWRK(*) - integer array, TRANSPORT workspace
C  RMCWRK(*) - integer array, TRANSPORT workspace
C  F(*,*)    - real matrix, governing equations residuals of S(*,*);
C              F(NT,J) is the energy equation residual,
C              F(NM,J) is the mass equation residual,
C              F(NYS+K,J) is the species equation for species K
C  XAV(*)    - real array, mole fractions at mesh midpoints
C  SSAVE(*)  - real array, temperature values at which stored reaction
C              rates were evaluated.
C  RKFT(*)   - real array, forward reaction rates at SSAVE.
C  RKRT(*)   - real array, reverse reaction rates at SSAVE.
C  ICASE     - integer scalar, flag for production rate evaluation;
C              =1, CALL CKWYP for production rates,
C              =2, save temperatures in SSAVE, CALL CKKFRT for RKFT
C                  and RKRT, and CALL CKWYPK for production rates,
C              =3, compare new temperature values to SSAVE, and if
C                  same, CALL CKKFRT as in (2), else do (1)
C  RSAVE(*,*)- real matrix, for ICASE=2, save species production rates,
C              for ICASE=3, use RSAVE for species production rates
C
C  END PROLOGUE
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
C     Integer arrays
      DIMENSION ICKWRK(LENICK), IMCWRK(LENIMC)
C     Real arrays
      DIMENSION WT(KK), EPS(KK), XGIVEN(NTEMP), TGIVEN(NTEMP), X(JJ),
     1          S(NATJ,JJ), SN(NATJ,JJ), YAV(KK), YV(KK,JJ), WDOT(KK),
     2          CP(KK), H(KK), COND(JJ), D(KK,JJ), TDR(KK,JJ),
     3          RCKWRK(LENRCK), RMCWRK(LENRMC), F(NATJ,JJ),
     4          DKJ(KK,KK,JJ), XAV(KK), SSAVE(NATJ,JJ), RKFT(II,JJ),
     5          RKRT(II,JJ), RSAVE(KK,JJ)
C
      LOGICAL LBURNR, LENRGY, LTDIF, LVARMC, LTIME, LMULTI, LVCOR, MATCH
      EXTERNAL CKBSEC, AREA, VTRACE
      LOGICAL MARCDBG

      logical LULNUM
      common / MARC / DFTH, CHTH, LULNUM

      double precision CONV,DIFF,REAC,DIFFE
      integer iBudget
      iBudget = 1

c      if (MARCDBG.eqv..TRUE.) print *,'In FUN MARCDBG = ',MARCDBG
C
C     Evaluate and store the transport coefficients
      IF (LVARMC) THEN
         IF (LMULTI) THEN
            CALL MCMULT (LENRGY, LTDIF, X, S, YAV, ICKWRK, RCKWRK,
     1                   IMCWRK, RMCWRK, WT, H, CP, XAV, COND, D,
     1                   TDR, DKJ)
         ELSE
            CALL MCMIX (LENRGY, LTDIF, X, S, YAV, ICKWRK, RCKWRK,
     1                  IMCWRK, RMCWRK, WT, H, CP, XAV, COND, D,
     2                  TDR, DKJ)
         ENDIF
      ENDIF

c     Unity Lewis number:   K = (1 / PRANDTL) * mu * CPMIX
      if ( LULNUM ) then
         do J=1,JJ-1
            CALL CKCPBS (S(NT,J),   S(NY,J),   ICKWRK, RCKWRK, CPB)
            CALL CKCPBS (S(NT,J+1), S(NY,J+1), ICKWRK, RCKWRK, CPBJ1)
            TAV = 0.5d0 * ( S(NT,J+1) + S(NT,J) )
            CPB = 0.5d0 * ( CPB + CPBJ1 )
            VISC = .185d0*(TAV/298.d0)**0.7
c            COND(J) = CPB * ( 1. / 0.7d0 ) * VISC
         enddo
      endif

C
C     Evaluate and store the diffusion velocities
C     (in the following call h(*) and cp(*) are used temporarily
C     for scratch space)
      CALL MDIFV (LMULTI, LVCOR, LTDIF, X, S, WT, YAV, H,
     1            CP, D, TDR, ICKWRK, RCKWRK, YV, DKJ)

c      open(unit=4411,file="budget.dat",status='unknown')
c      rewind(unit=4411)
c      write(4411,'(a)'),'VARIABLES=X T uGradY DivFD FcpDT Reac FDp'

c      print *,'ICASE, GFAC:',ICASE,GFAC
C
C     Left boundary
      XMID = 0.5 * (X(1)+X(2))
      TAV = 0.5 * (S(NT,1) + S(NT,2))
      CALL CKAVG (KK, S(NYS+1,1), S(NYS+1,2), YAV)
      AREAP = AREA (XMID)
      CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOP)
C
      RHOPA = RHOP * AREAP
      RHOPF = RHOPA / S(NM,1)
      DO 0020 K = 1, KK-1
         N = NYS + K
C         F(N,1) = EPS(K) - S(N,1) - RHOPA * YV(K,1) / S(NM,1)
         F(N,1) = EPS(K) - S(N,1) - YV(K,1) * RHOPF
0020  CONTINUE
      k = 1
C
      N = NYS + KK
      IF (LVCOR) THEN
C         F(N,1) = EPS(KK) - S(N,1) - RHOPA * YV(KK,1) / S(NM,1)
         F(N,1) = EPS(KK) - S(N,1) - YV(KK,1) * RHOPF
      ELSE
         F(N,1) = VTRACE (S(NYS+1,1), KK)
      ENDIF
C
C     Inflow boundary condition for temperature
      F(NT,1) = S(NT,1) - TGIVEN(1)
C
      IF (LBURNR) THEN
         F(NM,1) = S(NM,1) - FLRT
      ELSE
         F(NM,1) = S(NM,2) - S(NM,1)
      ENDIF
C
      DO 0120 J = 2, JJ-1
C        Interior mesh points
C
        JP1 = J + 1
        JM1 = J - 1
        TAV = 0.5 * (S(NT,J) + S(NT,JP1))
        CALL CKAVG (KK, S(NYS+1,J), S(NYS+1,JP1), YAV)
C
        XMDOT = S(NM,J)
        RHOM = RHOP
        CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOP)
C
        XMID = 0.5 * (X(J)+X(JP1))
        AREAM = AREAP
        AREAP = AREA(XMID)
        AREAJ = 0.5 * (AREAP+AREAM)
C
C       Form the chemical rate terms
        IF (ICASE .EQ. 1) THEN
           CALL CKWYP (P, S(NT,J), S(NY,J), ICKWRK, RCKWRK, WDOT)
C
        ELSEIF (ICASE .EQ. 2) THEN
           CALL CKCOPY (NATJ, S(1,J), SSAVE(1,J))
           CALL CKKFRT (P, S(NT,J), ICKWRK, RCKWRK, RKFT(1,J),
     1                  RKRT(1,J))
           CALL CKWYPK (P, S(NT,J), S(NY,J), RKFT(1,J), RKRT(1,J),
     1                  ICKWRK, RCKWRK, WDOT)
           CALL CKCOPY (KK, WDOT, RSAVE(1,J))
C
        ELSEIF (ICASE .EQ. 3) THEN
           MATCH = S(NT, J) .EQ. SSAVE(NT, J)
           DO 0060 K = 1, KK
              N = NYS + K
              MATCH = MATCH .AND. (S(N, J) .EQ. SSAVE(N, J))
0060       CONTINUE
C
           IF (MATCH) THEN
              CALL CKCOPY (KK, RSAVE(1,J), WDOT)
           ELSEIF (S(NT, J) .EQ. SSAVE(NT, J)) THEN
              CALL CKWYPK (P, S(NT,J), S(NY,J), RKFT(1,J), RKRT(1,J),
     1                     ICKWRK, RCKWRK, WDOT)
           ELSE
              CALL CKWYP (P, S(NT,J), S(NY,J), ICKWRK, RCKWRK, WDOT)
           ENDIF
        ENDIF
C
        IF (GFAC .NE. 1.0) THEN
           DO 0080 K = 1, KK
C             Damp the rates
              WDOT(K) = WDOT(K)*GFAC
0080       CONTINUE
        ENDIF

        IF (CHTH .NE. 1.0) THEN
           do K = 1, KK
C             Damp the rates
              WDOT(K) = WDOT(K) / CHTH
           enddo
        ENDIF
C
C       Evaluate enthalpies and heats
        CALL CKHML (S(NT,J), ICKWRK, RCKWRK, H)
        CALL CKCPBS (S(NT,J), S(NY,J), ICKWRK, RCKWRK, CPB)
        CALL CKCPMS (S(NT,J), ICKWRK, RCKWRK, CP)
C
C       Form the mesh differences
        DXP =        (X(JP1) - X(J)  )
        DXM =        (X(J)   - X(JM1))
        DXAV = 0.5 * (X(JP1) - X(JM1))
        DXPM =       (X(JP1) - X(JM1))
C
C       Form the coefficients for central differences
        CENDFM = - DXP / (DXM*DXPM)
        CENDFC =   (DXP-DXM) / (DXP*DXM)
        CENDFP =   DXM / (DXP*DXPM)
C
C       Species conservation equation
        RHOPA = RHOP * AREAP
        RHOMA = RHOM * AREAM
        IF (WNDFAC .EQ. 1.0) THEN
          XMDXM = XMDOT / DXM
          DO 0090 K = 1, KK-1
             N = NYS + K
             F(N,J) = XMDXM * ( S(N,J) - S(N,JM1) ) +
     1             (RHOPA*YV(K,J) - RHOMA*YV(K,JM1)) / DXAV -
     2               AREAJ * WDOT(K) * WT(K)
0090      CONTINUE
          N = NYS + KK
          IF (LVCOR) THEN
             F(N,J) = XMDXM * (S(N,J) - S(N,JM1)) +
     1            (RHOPA*YV(KK,J) - RHOMA*YV(KK,JM1))/ DXAV -
     2               AREAJ * WDOT(KK) * WT(KK)
          ELSE
             F(N,J) = VTRACE (S(NYS+1,J), KK)
          ENDIF

          CONV = XMDXM * ( S(NYS+iBudget,J) - S(NYS+iBudget,JM1) )

        ELSE
          DO 0100 K = 1, KK-1
             N = NYS + K
             F(N,J) = XMDOT *
     1                    (CENDFP*S(N,JP1) + CENDFC*S(N,J) +
     2                     CENDFM*S(N,JM1) )  +
     3             (RHOPA*YV(K,J) - RHOMA*YV(K,JM1)) / DXAV -
     4               AREAJ * WDOT(K) * WT(K)
0100      CONTINUE
          N = NYS + KK
          IF (LVCOR) THEN
             F(N,J) = XMDOT *
     1                    (CENDFP*S(N,JP1) + CENDFC*S(N,J) +
     2                     CENDFM*S(N,JM1) )  +
     3           (RHOPA*YV(KK,J) - RHOMA*YV(KK,JM1)) / DXAV -
     4               AREAJ * WDOT(KK) * WT(KK)
          ELSE
             F(N,J) = VTRACE (S(NYS+1,J), KK)
          ENDIF

          CONV = XMDOT *  (CENDFP*S(NYS+iBudget,JP1)
     &                   + CENDFC*S(NYS+iBudget,J) +
     2                     CENDFM*S(NYS+iBudget,JM1) )

        ENDIF


        DIFF = (RHOPA*YV(iBudget,J) - RHOMA*YV(iBudget,JM1)) / DXAV
        DIFFE = 0.d0
        REAC = AREAJ * WDOT(iBudget) * WT(iBudget)
c        write(4411,9988) X(J),S(NT,J),CONV,DIFF,DIFFE,REAC,
c     &       RHOP*AREAP*YV(iBudget,J)
            
C
C       Mass flow rate equation
        IF (LBURNR) THEN
           F(NM,J) = S(NM,J) - S(NM,JM1)
        ELSE
           IF (J .GT. JFIXT) THEN
              F(NM,J) = S(NM,J) - S(NM,JM1)
           ELSE
              F(NM,J) = S(NM,JP1) - S(NM,J)
           ENDIF
           IF (J .EQ. JFIXT) F(NM,J) = S(NT,J) - TFIXT
        ENDIF
C
C       Energy equation
        IF (LENRGY) THEN
C
           SUM = 0.0
           TDOT = 0.0
           CENT = CENDFP * S(NT,JP1) +
     1            CENDFC * S(NT,J)   +
     2            CENDFM * S(NT,JM1)

           DO 0110 K = 1, KK
              TDOT = TDOT + WDOT(K)*H(K)
              SUM = SUM + 0.5 * (RHOP*YV(K,J) + RHOM*YV(K,JM1)) *
     1                    CP(K) * CENT
0110       CONTINUE
C
           CONJ =   (COND(J)  * AREAP*(S(NT,JP1) - S(NT,J)) / DXP -
     1               COND(JM1)* AREAM*(S(NT,J) - S(NT,JM1)) / DXM) /
     2              (CPB*DXAV) - AREAJ * (SUM + TDOT) / CPB
C
           IF (WNDFAC .EQ. 1.0) THEN
             F(NT,J) = XMDOT * ( S(NT,J) - S(NT,JM1) ) / DXM - CONJ
           ELSE
             F(NT,J) = XMDOT * CENT - CONJ
           ENDIF
C
        ELSE
           F(NT, J) = S(NT,J) - CKBSEC(NTEMP, X(J), XGIVEN, TGIVEN)
        ENDIF
C
c        CONV = XMDOT * CENT
c        DIFF = ( COND(J)  * AREAP*(S(NT,JP1) - S(NT,J)) / DXP
c     1       -   COND(JM1)* AREAM*(S(NT,J) - S(NT,JM1)) / DXM) /
c     1       (DXAV * CPB)
c        DIFFE =  - AREADJ * SUM / CPB
c        REAC = AREAJ * TDOT / CPB
c        write(4411,9988) X(J),S(NT,J),CONV,DIFF,DIFFE,REAC,
c     &       COND(J)  * AREAP*(S(NT,JP1) - S(NT,J)) / DXP

0120  CONTINUE
C
C     Right boundary
      JJM1 = JJ - 1
      DO 0130 K = 1, KK-1
         N = NYS + K
         F(N,JJ) = S(N,JJ) - S(N,JJM1)
0130  CONTINUE
C
      N = NYS + KK
      IF (LVCOR) THEN
         F(N,JJ) = S(N,JJ) - S(N,JJM1)
      ELSE
         F(N,JJ) = VTRACE (S(NYS+1,JJ), KK)
      ENDIF
C
      IF (LENRGY) THEN
         F(NT, JJ) = S(NT,JJ) - S(NT,JJM1)
      ELSE
         F(NT, JJ) = S(NT,JJ) - CKBSEC(NTEMP, X(JJ), XGIVEN, TGIVEN)
      ENDIF
C
      F(NM, JJ) = S(NM,JJ) - S(NM,JJM1)
C
 9988 format(10e14.6)
c      close(4411)
c      print *,'stopped in FUN'
c      stop

      IF (.NOT. LTIME) RETURN
C
C     Add the time step, if needed
      IF (LVCOR) THEN
         KK1 = KK
      ELSE
         KK1 = KK - 1
      ENDIF
C
      DO 0150 J = 2, JJ-1
         CALL CKRHOY (P, S(NT,J), S(NY,J), ICKWRK, RCKWRK, RHO)
         RHODT = RHO * AREA(X(J)) / DT
         IF (LENRGY) F(NT,J) = F(NT,J) + RHODT * (S(NT,J) - SN(NT,J))
         DO 0150 K = 1, KK1
            N = NYS + K
            F(N,J) = F(N,J) + RHODT * (S(N,J) - SN(N,J))
0150  CONTINUE

C
C     end of SUBROUTINE FLDRIV
      RETURN
      END
C
      SUBROUTINE MDIFV (LMULTI, LVCOR, LTDIF, X, S, WT, YAV, XMF, XMFP,
     1                  D, TDR, ICKWRK, RCKWRK, YV, DKJ)
C
C  START PROLOGUE
C
C  LMULTI     - logical, T, use multicomponent formulas,
C                        F, use mixture-averaged formulas
C  LVCOR      - logical, T, use correction velocity formulism
C                        F, use 'trace' aproximation, lumping
C                           all transport errors into final species
C  LTDIF      - logical, T, evaluate thermal diffusion ratios
C                           as well as diffusion coefficients
C  X(*)       - real array, mesh point locations (cm)
C  S(*,*)     - real matrix, dependent variables;
C                            S(NT,J) = T(J), temperature
C                            S(NYS+K,J) = Y(K,J), mass fractions
C                            S(NM,J) = FLRT(J), flow rate
C  WT(*)      - real array, species molecular weights (gm/mole)
C  YAV(*)     - real array, mass fractions at mesh midpoints;
C               YAV(K) is between J and J+1
C  XMF(*)     - real array, mole fractions at mesh points
C  XMFP(*)    - real array,
C  D(*,*)     - real matrix, species diffusion coeff'nts at mesh
C               midpoints (cm^2/sec);
C               if LVARMC=.TRUE., these are computed each time function
C               is called, otherwise the stored values are used.
C  TDR(*,*)   - real matrix, species thermal diffusion ratios at mesh
C               midpoints;
C               if LVARMC=.TRUE., these are comptued each time function
C               is called, otherwise the stored values are used.
C  ICKWRK(*)  - integer array, CHEMKIN workspace
C  RCKWRK(*)  - real array, CHEMKIN workspace
C  YV(*,*)    - real matrix, mass fractions times diffusion velocities
C               at mesh midpoints;
C               YV(K,J) is flux of species K between J and J+1.
C  DKJ(*,*)   - real 3-dimensional array, multicomponent diffusion
C               coefficients at mesh midpoints (cm^2/sec)
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
C
      logical LULNUM
      common / MARC / DFTH, CHTH, LULNUM

      DIMENSION X(JJ), S(NATJ,JJ), WT(KK), YAV(KK), XMF(KK), XMFP(KK),
     1          D(KK,JJ), TDR(KK,JJ), YV(KK,JJ), DKJ(KK,KK,JJ), XAV(KK)
C
      LOGICAL LTDIF, LMULTI, LVCOR
C
      CALL CKYTX (S(NY,1), ICKWRK, RCKWRK, XMFP)
C
C     Loop over all mesh points, computing the diffusion velocity
C     at the midpoints; the indexing is such that YV(K,J) is the
C     diffusion velocity of the Kth species midway between nodes
C     J and J+1.
C
      DO 1000 J = 1, JJ-1
         JP1 = J + 1
         JM1 = J - 1
C
         TAV = 0.5 * (S(NT,J) + S(NT,JP1))
         CALL CKAVG (KK, S(NYS+1,J), S(NYS+1,JP1), YAV)
         CALL CKCOPY (KK, XMFP, XMF)
         CALL CKMMWY (YAV, ICKWRK, RCKWRK, WTM)
         CALL CKYTX (S(NY,JP1), ICKWRK, RCKWRK, XMFP)
         XDIF = X(JP1) - X(J)
C

         if ( LULNUM ) then
c     Unity Lewis number.
c     First, compute D = (DFTH / SCHMIDT) * mu   (same for all species)
            VISC = .185d0*(TAV/298.d0)**0.7
            D(1,J) = ( DFTH / 0.7d0 ) * VISC
            
c     Then, form diffusion velocities
            do K = 1, KK
               YV(K,J) = - D(1,J) * ( S(NYS+K,JP1) - S(NYS+K,J) ) / XDIF
            enddo
            
         else

         IF ( LMULTI ) THEN
C           Evaluate the multicomponent diffusion velocity directly,
C           rather than use the mixture-averaged form for D(K,J)
C
            WTM2 = WTM**2
            DO 475 K = 1, KK
               SUM = 0.0
               DO 450 L = 1, KK
                  SUM = SUM + WT(L) * DKJ(K,L,J) *
     1                     (XMFP(L)-XMF(L)) / XDIF
  450          CONTINUE
               YV(K,J) = (WT(K)/WTM2) * SUM
  475       CONTINUE
         ELSE
C           Use mixture-averaged form for Fickian diffusion,
C           whether we are using the multicomponent formalism
C           or mixture-averaged
C
            DO 500 K = 1, KK
               YV(K,J) = - D(K,J) * (WT(K)/WTM) *
     1                     (XMFP(K)-XMF(K)) / XDIF
500         CONTINUE
         ENDIF
C
         IF (LTDIF) THEN
C           add thermal diffusion, if requested
C
            CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOAV)
            TRHO = TAV * RHOAV
            TDIF = (S(NT, JP1) - S(NT,J)) / (XDIF * TRHO)
            DO 600 K = 1, KK
               YV(K,J) = YV(K,J) - TDIF * TDR(K,J)
600         CONTINUE
C
         ENDIF
         endif
C
         IF (LVCOR) THEN
C           compute and add the correction velocity
C
            SUM = 0.0
            DO 700 K = 1, KK
               SUM = SUM + YV(K,J)
700         CONTINUE
C
            VC = - SUM
C
            DO 800 K = 1, KK
               YV(K,J) = YV(K,J) + YAV(K)*VC
800         CONTINUE
         ENDIF
C
1000  CONTINUE
C
C     end of SUBROUTINE MDIFV
      RETURN
      END
C
      SUBROUTINE MCMIX (LENRGY, LTDIF, X, S, YAV, ICKWRK,
     1                   RCKWRK, IMCWRK, RMCWRK, WT, XMF, XMFP,
     2                   XAV, COND, D, TDR, DKJ)
C
C  START PROLOGUE
C
C  LENRGY     - logical, T  solve energy equation,
C                        F, use fixed temperature profile
C  LTDIF      - logical, T, evaluate thermal diffusion ratios as
C                           well as diffusion coefficients
C  X(*)       - real array, mesh point locations (cm)
C  S(*,*)     - real matrix, dependent variables;
C                            S(NT,J) = T(J), temperature
C                            S(NYS+K,J) = Y(K,J), mass fractions
C                            S(NM,J) = FLRT(J), flow rate
C  YAV(*)    - real array, mass fractions at mesh midpoints;
C              YAV(K) is between J and J+1
C  ICKWRK(*) - integer array, CHEMKIN workspace
C  RCKWRK(*) - real array, CHEMKIN workspace
C  IMCWRK(*) - integer array, TRANSPORT workspace
C  RMCWRK(*) - integer array, TRANSPORT workspace
C  WT(*)      - real array, species molecular weights (gm/mole)
C  XMF(*)     - real array, mole fractions at mesh points
C  XMFP(*)    -
C  XAV(*)    - real array, mole fractions at mesh midpoints
C  COND(*)   - real array, thermal conductivities at mesh midpoints
C              (erg/cm-K-sec);
C              if LVARMC=.TRUE., these are computed each time function
C              is called, otherwise the stored values are used
C  D(*,*)    - real matrix, species diffusion coefficients at mesh
C              midpoints (cm^2/sec);
C              if LVARMC=.TRUE., these are computed each time function
C              is called, otherwise the stored values are used.
C  TDR(*,*)  - real matrix, species thermal diffusion ratios at mesh
C              midpoints;
C              if LVARMC=.TRUE., these are computed each time function
C              is called, otherwise the stored values are used.
C  DKJ(*,*,*)- real 3-dimensional array, multicomponent diffusion
C              coefficients at mesh midpoints (cm^2/sec)
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
C
      DIMENSION X(JJ), S(NATJ,JJ), YAV(KK), ICKWRK(LENICK),
     1          RCKWRK(LENRCK),
     1          IMCWRK(LENIMC), RMCWRK(LENRMC),
     2          COND(JJ), D(KK,JJ), TDR(KK,JJ),
     2          DKJ(KK,KK,JJ), XAV(KK), WT(KK), XMF(KK), XMFP(KK)
C
      LOGICAL LENRGY, LTDIF
C
      DO 400 J = 1, JJ-1
         JP1 = J + 1
         TAV = 0.5 * (S(NT,J) + S(NT,JP1))
C
C        Dimensional temperature at the grid points
         CALL CKAVG (KK, S(NYS+1,J), S(NYS+1,JP1), YAV)
         CALL CKYTX (YAV, ICKWRK, RCKWRK, XAV)
         CALL MCADIF(P, TAV, XAV, RMCWRK, D(1,J) )
C
C        Determine the mixture conductivity at J
         IF (LENRGY) CALL MCACON( TAV, XAV, RMCWRK, COND(J) )
C
         IF (LTDIF) THEN
            CALL MCATDR( TAV, XAV, IMCWRK, RMCWRK, TDR(1,J) )
            CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOAV)
            CALL CKMMWY (YAV, ICKWRK, RCKWRK, WTM)
            DO 350 K = 1, KK
               TDR(K,J) = D(K,J) * TDR(K,J) * RHOAV * WT(K)/WTM
  350       CONTINUE
         ENDIF
C
400   CONTINUE
C
C     end of SUBROUTINE MCMIX
      RETURN
      END
      SUBROUTINE MCMULT (LENRGY, LTDIF, X, S, YAV, ICKWRK, RCKWRK,
     1                   IMCWRK, RMCWRK, WT, XMF, XMFP, XAV, COND,
     2                   D, TDR, DKJ)
C
C  START PROLOGUE
C
C  LENRGY     - logical, T  solve energy equation,
C                        F, use fixed temperature profile
C  LTDIF      - logical, T, evaluate thermal diffusion ratios as
C                           well as diffusion coefficients
C  X(*)       - real array, mesh point locations (cm)
C  S(*,*)     - real matrix, dependent variables;
C                            S(NT,J) = T(J), temperature
C                            S(NYS+K,J) = Y(K,J), mass fractions
C                            S(NM,J) = FLRT(J), flow rate
C  YAV(*)    - real array, mass fractions at mesh midpoints;
C              YAV(K) is between J and J+1
C  ICKWRK(*) - integer array, CHEMKIN workspace
C  RCKWRK(*) - real array, CHEMKIN workspace
C  IMCWRK(*) - integer array, TRANSPORT workspace
C  RMCWRK(*) - integer array, TRANSPORT workspace
C  WT(*)      - real array, species molecular weights (gm/mole)
C  XMF(*)     - real array, mole fractions at mesh points
C  XMFP(*)    -
C  XAV(*)    - real array, mole fractions at mesh midpoints
C  COND(*)   - real array, thermal conductivities at mesh midpoints
C              (erg/cm-K-sec);
C              if LVARMC=.TRUE., these are computed each time function
C              is called, otherwise the stored values are used
C  D(*,*)    - real matrix, species diffusion coefficients at mesh
C              midpoints (cm^2/sec);
C              if LVARMC=.TRUE., these are computed each time function
C              is called, otherwise the stored values are used.
C  TDR(*,*)  - real matrix, species thermal diffusion ratios at mesh
C              midpoints;
C              if LVARMC=.TRUE., these are computed each time function
C              is called, otherwise the stored values are used.
C  DKJ(*,*,*)- real 3-dimensional array, multicomponent diffusion
C              coefficients at mesh midpoints (cm^2/sec)
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
C
      DIMENSION X(JJ), S(NATJ,JJ), YAV(KK), ICKWRK(LENICK),
     1          RCKWRK(LENRCK),
     1          IMCWRK(LENIMC),
     2          RMCWRK(LENRMC), COND(JJ), D(KK,JJ), TDR(KK,JJ),
     2          DKJ(KK,KK,JJ), XAV(KK), WT(KK), XMF(KK), XMFP(KK)
C
      LOGICAL LENRGY, LTDIF
C
      DATA EPS/1.0E-30/
C
      CALL CKYTX (S(NY,1), ICKWRK, RCKWRK, XMFP)
      DO 200 J = 1, JJ-1
         JP1 = J + 1
         TAV = 0.5 * (S(NT,J) + S(NT,JP1))
         XDIF = X(JP1) - X(J)
C
C        Dimensional temperature at the gridpoints
         CALL CKAVG (KK, S(NYS+1,J), S(NYS+1,JP1), YAV)
         CALL CKCOPY (KK, XMF, XMFP)
         CALL CKYTX (YAV, ICKWRK, RCKWRK, XAV)
         CALL CKYTX (S(NY,JP1), ICKWRK, RCKWRK, XMFP)
         CALL CKMMWX ( XAV, ICKWRK, RCKWRK, WTMAV)
         CALL MCMDIF(P, TAV, XAV, KK, IMCWRK, RMCWRK, DKJ(1,1,J) )
C
         DO 75 K = 1, KK
            SUMN = 0.0
            DO 50 L = 1, KK
               N = NYS + L
               SUMN = SUMN + DKJ(K,L,J) *
     1              (S(N,JP1) - S(N,J)) / XDIF
   50       CONTINUE
            N = NYS + K
            DENOM = - (S(N,JP1) - S(N,J)) / XDIF
            D(K,J) = (SUMN + EPS) / ( WTMAV * (DENOM + EPS))
   75    CONTINUE

C
C        Determine the mixture conductivity and thermal diffusion
C        coefficient at J
C
         IF (LENRGY .OR. LTDIF) CALL MCMCDT
     1      ( P, TAV, XAV, IMCWRK, RMCWRK, ICKWRK, RCKWRK,
     2        TDR(1,J), COND(J) )
C
200   CONTINUE
C
C     end of SUBROUTINE MCMULT
      RETURN
      END
      SUBROUTINE LINWMX (WMIX, XCEN, XNODE, XRE, XPD)
C
C  START PROLOGUE
C
C  WMIX     - real scalar, width of the mixing region over which
C             starting estimates are fit (cm)
C  XCEN     - real scalar, X position where the initial starting
C             estimates are centered (cm)
C  XNODE    - real scalar
C  XRE      - real scalar
C  XPD      - real scalar
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      ZERO = 0.0
      WMID = 0.5*WMIX
C
      IF (XNODE .LE. XCEN-WMID) THEN
         XPD = ZERO
      ELSE
         IF (XNODE .LT. XCEN+WMID) THEN
            XPD = (XNODE-XCEN)/WMIX + 0.5
         ELSE
            XPD = 1.0
         ENDIF
      ENDIF
C
      XRE = 1.0 - XPD
C
C     end of SUBROUTINE LINWMX
      RETURN
      END
C
      SUBROUTINE POINTR (LINKCK, LINKMC, LENIWK, LENRWK, LENCWK, JMAX,
     1                   LOUT, LSAVE, LTOT, ITOT, NTOT, ICTOT, I, R, C,
     2                   NIWK, NRWK)
C
C  START PROLOGUE
C
C  LINKCK     - integer scalar, CHEMKIN linkfile input unit number
C  LINKMC     - integer scalar, TRANSPORT linkfile input unit number
C  LENIWK     - integer scalar, length of integer workspace provided
C  LENRWK     - integer scalar, length of real workspace provided
C  LENCWK     - integer scalar, length of character-string workspace
C  JMAX       - integer scalar, maximum gridpoints
C  LOUT       - integer scalar, formatted output file unit number
C  LSAVE      - integer scalar, binary solution file unit number
C  LTOT       - integer scalar, amount of logical workspace required
C  ITOT       - integer scalar, amount of integer workspace required
C  NTOT       - integer scalar, amount of real workspace required
C  ICTOT      - integer scalar, amount of character workspace required
C  I(*)       - integer array, problem workspace
C  R(*)       - real array, problem workspace
C  C(*)       - character-string array, problem workspace
C  NIWK       - integer scalar, size of twopnt integer workspace
C  NRWK       - integer scalar, size of twopnt real workspace
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
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
      COMMON /FLFLFL/
     +   ICC, ICKW, IIP, IKI, IKP, IKR, IKS, IMCW, INAME, LAC, LMK, NA,
     +   NABV, NBLW, NBUF, NCKW, NCON, NREG, ND, NDKJ, NDS, NEPS, NF,
     +   NFN,  NKA6, NMCW, NRE,  NRKFT,NRKRT,NRSAVE, NS, NSCH, NSN,
     +   NTDR, NTGV, NSSAVE, NWT, NX, NXGV, NYV
C
      DIMENSION I(LENIWK), R(LENRWK)
      CHARACTER C(LENCWK)*(*)
C
C      SET THE POINTERS INTO THE SOLUTION VECTOR
C
C      NT  = 1
C      NYS = 1
C      NY  = 2
C      NM  = KK+2
C
      CALL CKLEN (LINKCK, LOUT, LENICK, LENRCK, LENCCK, IFLAG)
      CALL MCLEN (LINKMC, LOUT, LENIMC, LENRMC, IFLAG)
C
C     real chemkin work space
      NCKW = 1
C     real transport work space
      NMCW = NCKW + LENRCK
      NTOT = NMCW + LENRMC
C
C     integer chemkin work space
      ICKW = 1
C     integer transport work space
      IMCW = ICKW + LENICK
      ITOT = IMCW + LENIMC
C
C     character chemkin work space
      ICC = 1
      ICTOT = ICC + LENCCK
C
C
      IF (ITOT.LT.LENIWK .AND. NTOT.LT.LENRWK .AND. ICTOT.LT.LENCWK)
     1   THEN
         CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT, I, R, C,
     1        IFLAG)
         CALL CKINDX (I, R, MM, KK, II, NFIT)
C
         CALL CKMXTP (I, MAXTP)
         NTR = MAXTP - 1
C
         CALL MCINIT (LINKMC, LOUT, LENIMC, LENRMC, I(IMCW), R(NMCW),
     1        IFLAG)
         REWIND LSAVE
         CALL PRSAVE (I, R, C, I(IMCW), R(NMCW), LOUT, LSAVE)
      ENDIF
C
      NATJ = KK+2
C
C          APPORTION THE BALANCE OF THE FLOATING POINT SPACE
C
      NEPS = NTOT
      NWT  = NEPS + KK
      NRE  = NWT  + KK
      NSCH = NRE  + KK
      NX   = NSCH + KK*6
      NCON = NX   + JMAX
      NREG = NCON + JMAX
      NTGV = NREG + JMAX
      NXGV = NTGV + JMAX
      ND   = NXGV + JMAX
      NDKJ = ND   + KK*JMAX
      NTDR = NDKJ + KK*KK*JMAX
      NYV  = NTDR + KK*JMAX
      NABV = NYV  + KK*JMAX
      NBLW = NABV + NATJ
      NBUF = NBLW + NATJ
      NS   = NBUF + NATJ * JMAX
      NSN  = NS   + NATJ * JMAX
      NF   = NSN  + NATJ * JMAX
      NFN  = NF   + NATJ * JMAX
      NDS  = NFN  + NATJ * JMAX
C     thermodynamic coefficients a6
      NKA6 = NDS  + JMAX
C     jacobian space
      NA   = NKA6 + NTR
      IASIZE = (6 * NATJ - 1) * (NATJ * JMAX)
C     twopnt real work space
      NRWK = NA + IASIZE
      LENRTW = 3 * JMAX + 9 * NATJ * JMAX
C     saved solution in FUN
      NSSAVE = NRWK + LENRTW - 1
C     saved forward rates
      NRKFT = NSSAVE + NATJ * JMAX
C     saved reverse rates
      NRKRT = NRKFT + II * JMAX
C     saved production rates
      NRSAVE = NRKRT + II * JMAX
      NTOT = NRSAVE + KK * JMAX
C
C           APPORTION THE BALANCE OF THE INTEGER SPACE
C
      IKR  = ITOT
      IKI  = IKR  + KK
      IKP  = IKI  + KK
      IIP  = IKP  + KK
C     twopnt real work space
      NIWK = IIP  + NATJ * JMAX
      LENITW = 3 * JMAX
      ITOT = NIWK + LENITW - 1
C
C           APPORTION THE LOGICAL SPACE
C
      LAC  = 1
      LMK  = LAC  + NATJ
      LTOT = LMK  + JMAX - 1
C
C           APPORTION THE BALANCE OF THE CHARACTER SPACE
C
      IKS = ICTOT
      INAME = IKS + KK + 2
      ICTOT = INAME + KK - 1
C
C     end of SUBROUTINE POINTR
      RETURN
      END
C
      SUBROUTINE PREMIX (JMAX, LIN, LOUT, LINKCK, LINKMC, LREST, LSAVE,
     1                   LRCRVR, LENLWK, L, LENIWK, I, LENRWK, R,
     2                   LENCWK, C)
C
C  START PROLOGUE
C
C  JMAX     - integer scalar, maximum number of gridpoints allowed
C  LIN      - integer scalar, formatted input file unit number
C  LOUT     - integer scalar, formatted output file unit number
C  LINKCK   - integer scalar, CHEMKIN linkfile input unit number
C  LINKMC   - integer scalar, TRANSPORT linkfile input unit number
C  LREST    - integer scalar, binary restart input file unit number
C  LSAVE    - integer scalar, binary solution output file unit number
C  LRCRVR   - integer scalar, binary scratch file unit number
C  LENLWK   - integer scalar, size of logical problem workspace
C  L(*)     - logical array, problem workspace
C  LENIWK   - integer scalar, size of integer problem workspace
C  I(*)     - integer array, problem workspace
C  LENRWK   - integer scalar, size of real problem workspace
C  R(*)     - real array, problem workspace
C  LENCWK   - integer scalar, size of character-string problem workspace
C  C(*)     - character-string array, problem workspace
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
      COMMON /FLFLFL/
     +   ICC, ICKW, IIP, IKI, IKP, IKR, IKS, IMCW, INAME, LAC, LMK, NA,
     +   NABV, NBLW, NBUF, NCKW, NCON, NREG, ND, NDKJ, NDS, NEPS, NF,
     +   NFN,  NKA6, NMCW, NRE,  NRKFT,NRKRT,NRSAVE, NS, NSCH, NSN,
     +   NTDR, NTGV, NSSAVE, NWT, NX, NXGV, NYV
C
      DIMENSION I(LENIWK), R(LENRWK)
      CHARACTER C(LENCWK)*(*), PRVERS*16, PRDATE*16, PREC*16
      LOGICAL L(LENLWK)
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      DATA PRVERS/'3.15'/, PRDATE/'98/03/03'/
C
C*****precision > double
      PREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C      PREC = 'SINGLE'
C*****END precision > single
C
      WRITE (LOUT, '(/A, A,/1X,A, A, A, A, /A, /A, //)')
     1' PREMIX: ',
     2'CHEMKIN-III One-dimensional steady premixed laminar flame code,',
     3PREC(1:CKLSCH(PREC)), ' PRECISION Vers. ',
     4PRVERS(1:CKLSCH(PRVERS)+1), PRDATE,
     5' Copyright 1995, Sandia Corporation.',
     6' The U.S. Government retains a limited license in this software.'
C
C     Set up internal work pointers
      CALL POINTR (LINKCK, LINKMC, LENIWK, LENRWK, LENCWK, JMAX, LOUT,
     1             LSAVE, LTOT, ITOT, NTOT, ICTOT, I, R, C, NIWK, NRWK)
C
C     Check for enough space
      WRITE (LOUT, 7000) LENLWK, LTOT, LENIWK, ITOT, LENRWK, NTOT,
     1                   LENCWK, ICTOT
7000  FORMAT (/,'                WORKING SPACE REQUIREMENTS',
     1        /,'                 PROVIDED        REQUIRED ',
     2        /,' LOGICAL  ' , 2I15,
     3        /,' INTEGER  ' , 2I15,
     4        /,' REAL     ' , 2I15,
     5        /,' CHARACTER' , 2I15,/)
C
      IF (LTOT.GT.LENLWK .OR. ITOT.GT.LENIWK .OR. NTOT.GT.LENRWK
     1                   .OR. ICTOT.GT.LENCWK) THEN
         WRITE (LOUT, *) '  FATAL ERROR, NOT ENOUGH WORK SPACE PROVIDED'
         RETURN
      ENDIF
C
      CALL FLDRIV (LIN, LOUT, LREST, LSAVE, LRCRVR, LINKCK, LINKMC,
     1             JMAX, R(NCKW), R(NMCW), R(NEPS), R(NWT), R(NRE),
     2             R(NSCH), R(NX), R(NCON), R(NREG), R(NTGV), R(NXGV),
     3             R(ND), R(NDKJ), R(NTDR), R(NYV), R(NABV), R(NBLW),
     4             R(NBUF), R(NS), R(NSN), R(NF), R(NFN), R(NDS),
     5             R(NKA6), R(NA), I(ICKW), I(IMCW), C(IKS), C(ICC),
     6             I(IKR), I(IKI), I(IKP), I(IIP), L(LAC), L(LMK),
     7             C(INAME), I(NIWK), R(NRWK), R(NSSAVE), R(NRKFT),
     8             R(NRKRT), R(NRSAVE))
C
C     end of SUBROUTINE PREMIX
      RETURN
      END
C
      SUBROUTINE PRETWO (LBURNR, LENRGY, LMULTI, LVCOR, LTDIF, LVARMC,
     1                   LTIME, LMOLE, WT, EPS, XGIVEN, TGIVEN, X, S,
     2                   SN, BUFFER, SCRTCH, YV, COND, D, DKJ, TDR,
     3                   ICKWRK, RCKWRK, CCKWRK, IMCWRK, RMCWRK, ITWWRK,
     4                   RTWWRK, F, FN, SSAVE, RKFT, RKRT, RSAVE, LOUT,
     5                   LRCRVR, A, CONDIT, IPIVOT, NAME, KSYM, ABOVE,
     6                   BELOW, JMAX, MARK, ACTIVE, N1CALL, KERR)
C
C  START PROLOGUE
C
C  LBURNR     - logical, T  burner stabilized flame,
C                        F, freely propagating adiabatic flame
C  LENRGY     - logical, T  solve energy equation,
C                        F, use fixed temperature profile
C  LMULTI     - logical, T, use multicomponent formulas,
C                        F, use mixture-averaged formulas
C  LVCOR      - logical, T, use correction velocity formulism
C                        F, use 'trace' approximation, lumping
C                           all transport errors into final species
C  LTDIF      - logical, T, evaluate thermal diffusion ratios as
C                           well as diffusion coefficients
C  LVARMC     - logical, T, compute new transport properties
C                        F, use previously stored values
C  LTIME      - logical, T, time step of DT will be added to residual
C  LMOLE     - logical, TRUE for input/output in mole fraction,
C                      FALSE for input/output in mass fraction.
C  WT(*)      - real array, species molecular weights (gm/mole)
C  EPS(*)     - real array, mass flux fractions at the burner
C  XGIVEN(*)  - real array, distances from the burner for TGIVEN
C  TGIVEN(*)  - real array, temperatures at given distances (K)
C  X(*)       - real array, mesh point locations (cm)
C  S(*,*)     - real matrix, dependent variables;
C                            S(NT,J) = T(J), temperature
C                            S(NYS+K,J) = Y(K,J), mass fractions
C                            S(NM,J) = FLRT(J), flow rate
C  SN(*,*)    -
C  BUFFER(*,*)-
C  SCRTCH(*,*)-
C  YV(*,*)    - real matrix, mass fractions time diffusion velocities at
C               mesh midpoints;
C               YV(K,J) is flux of species K between J and J+1.
C  COND(*)    - real array, thermal conductivities at mesh midpoints
C               (erg/cm-K-sec);
C               if LVARMC=.TRUE., these are computed each time function
C               is called, otherwise the stored values are used
C  D(*,*)     - real matrix, species diffusion coefficients at mesh
C               midpoints (cm^2/sec);
C               if LVARMC=.TRUE., these are computed each time function
C               is called, otherwise the stored values are used
C  DKJ(*,*,*) - real 3-dimensional array, multicomponent diffusion
C               coefficients at mesh midpoints (cm^2/sec)
C  TDR(*,*)   - real matrix, species thermal diffusion ratios at mesh
C               midpoints;
C               if LVARMC=.TRUE., these are computed each time fucntion
C               is called, otherwise the stored values are used.
C  ICKWRK(*)  - integer array, CHEMKIN workspace
C  RCKWRK(*)  - real array, CHEMKIN workspace
C  CCKWRK(*)  - character-string array, CHEMKIN workspace
C  IMCWRK(*)  - integer array, TRANSPORT workspace
C  RMCWRK(*)  - real array, TRANSPORT workspace
C  ITWWRK(*)  - integer array, TWOPNT workspace
C  RTWWRK(*)  - real array, TWPONT workspace
C  F(*,*)     - real matrix, governing equations residuals of S(*,*);
C               F(NT,J) is the energy equation residual,
C               F(NM,J) is the mass equation residual,
C               F(NYS+K,J) is the species equation for species K
C  FN(*,*)    -
C  SSAVE(*)   - real array, temperature values at which stored reaction
C               rates were evaluated.
C  RKFT(*)    - real array, forward reaction rates at SSAVE.
C  RKRT(*)    - real array, reverse reaction rates at SSAVE.
C  RSAVE(*)   - real matrix, for ICASE=2, save species production rates,
C               for ICASE=3, use RSAVE for species production rates
C  LOUT       - integer scalar, formatted output file unit number
C  LRCRVR     - integer scalar, binary scratch file unit number
C  A(*)       -
C  CONDIT     -
C  IPIVOT(*)  -
C  NAME(*)    - character-string array, dependent variable names
C  KSYM(*)    - character-string array, species names
C  ABOVE(*)   - real array, dependent variable upper bounds
C  BELOW(*)   - real array, dependent variable lower bounds
C  JMAX       - integer scalar, maximum gridpoints
C  MARK(*)    -
C  ACTIVE(*)  -
C  N1CALL     -
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
      INTEGER CALL, CALLS
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
C
      CHARACTER VERSIO*80, SIGNAL*16, REPORT*16, NAME(NATJ)*(*),
     1          KSYM(KK)*(*), CCKWRK(LENCCK)*(*), ISOLUT*16
      LOGICAL LENRGY, ENERGY, LBURNR, BURNER, LMULTI, LVCOR, LVARMC,
     1        LTDIF, LTIME, LMOLE, ACTIVE(NATJ), ERROR, MARK(JMAX),
     2        RETURN, KERR, latom, lfuel, lpath, lrrop, lsolve, ldflx,
     +        lewis, lheat, ldiff
C     Integer arrays
      DIMENSION ICKWRK(LENICK), IMCWRK(LENIMC), ITWWRK(LENITW),
     1          IPIVOT(NATJ,JMAX)
C     Real arrays
      DIMENSION X(JMAX), S(NATJ,JMAX), ABOVE(NATJ), BELOW(NATJ),
     1          BUFFER(NATJ,JMAX), RCKWRK(LENRCK), RMCWRK(LENIMC),
     3          RTWWRK(LENRTW), A(IASIZE), RKFT(II,JMAX),
     4          RKRT(II,JMAX), D(KK,JMAX), F(NATJ,JMAX),
     5          FN(NATJ,JMAX), SN(NATJ,JMAX), TDR(KK,JMAX),
     6          TGIVEN(JMAX), SSAVE(NATJ,JMAX), WT(KK), XGIVEN(JMAX),
     7          YV(KK,JMAX), RSAVE(KK,JMAX), SCRTCH(KK,6), EPS(KK),
     8          DKJ(KK,KK,JMAX), COND(JMAX)
C
      EXTERNAL CKBSEC
      PARAMETER (ISOLUT='SOLUTION')

      LOGICAL MARCDBG

      common / prlcon / 
     +   latom, ldflx, ldiff, lewis, lfuel, lheat, lpath, lrrop, lsolve
C
C///////////////////////////////////////////////////////////////////////
C
C     (2) TOP OF THE BLOCKS TO CALL TWOPNT.
C
C///////////////////////////////////////////////////////////////////////
C
      KERR = .FALSE.
C///  CHOOSE THE VERSION.
C
C*****precision > double
      VERSIO = 'DOUBLE PRECISION VERSION 3.18'
C*****END precision > double
C*****precision > single
C      VERSIO = 'SINGLE PRECISION VERSION 3.18'
C*****END precision > single
C
      NAME(1) = 'TEMPERATURE'
      DO 0020 K = 1, KK
         NAME(1 + K) = KSYM(K)
0020  CONTINUE
      NAME(KK + 2) = 'MASS FLOW RATE'
C
C     Set the solution bounds.
      ZERO = 0.0
      BELOW(NT) = 50.0
      ABOVE(NT) = 6000.0
C     JFG: CHANGED LOWER BOUND OF MASS FLOW RATE TO POSTIVE VALUE
C     TO AVOID DIVISION BY ZERO WHEN LOWER BOUND IS USED AS VALUE
      BELOW(NM) = ZERO
      BELOW(NM) = 1.0E-6
      ABOVE(NM) = 10000.0
C
      DO 1030 K = 1, KK
         N = NYS + K
         BELOW (N) = SFLR
         ABOVE (N) = 10.0
1030  CONTINUE
C
C     Set the variables on which to adapt.
      ACTIVE(NM) = .FALSE.
      ACTIVE(NT) = LENRGY
      DO 1040 K = 1, KK
         N = NYS + K
         ACTIVE(N) = .TRUE.
1040  CONTINUE
C
C     Decide how many times to call TWOPNT.
      ENERGY = LENRGY
      BURNER = LBURNR
      IF (ENERGY) THEN
         CALLS = 2
      ELSE
         CALLS = 1
      ENDIF
      if (.not. lsolve) calls = 0
C
C     Set XFIXT at TFIXT.
      IF (.NOT. LBURNR) XFIXT = X(JFIXT)
C
C     Top of the loop over calls to TWOPNT.
      IGRPA = 0
      IGRPB = 0
      ERROR = .FALSE.
C
      DO 2070 CALL = N1CALL, CALLS
C
         IF (CALL .EQ. 2) WRITE (LOUT, '( / A / )')
     1   ' FLDRIV: FINISHED FIXED TEMPERATURE, ADDING ENERGY EQUATION'
C
         IF (ENERGY) THEN
            IF (CALL .EQ. 1) THEN
               LENRGY = .FALSE.
               LBURNR = .TRUE.
            ELSE
               LENRGY = .TRUE.
               LBURNR = BURNER
            ENDIF
         ENDIF
C
         ACTIVE(NT) = LENRGY
C
         SIGNAL = ' '
2010     CONTINUE
C
         IF (CALL .GT. 1) THEN
            CALL TWSETL (ERROR, LOUT, 'ADAPT', .TRUE.)
            KERR = KERR.OR.ERROR
            CALL TWSETI (ERROR, LOUT, 'STEPS1', NUMDT2)
            KERR = KERR.OR.ERROR
            CALL TWSETR (ERROR, LOUT, 'STRID0', DT2)
            KERR = KERR.OR.ERROR
            IF (KERR) RETURN
         ENDIF
C
         CALL TWOPNT (ERROR, LOUT, VERSIO, ABOVE, ACTIVE, BELOW, BUFFER,
     1                NATJ, CONDIT, IGRPA, IGRPB, LENITW, ITWWRK, MARK,
     2                NAME, NATJ, JMAX, JJ, REPORT, LENRTW, RTWWRK,
     3                SIGNAL, DT, LTIME, S, X)
         KERR = KERR.OR.ERROR
         IF (KERR) RETURN
C
C        Top of the block to service requests from TWOPNT.

         IF (SIGNAL .NE. ' ') THEN
C
            IF (SIGNAL .EQ. 'RESIDUAL') THEN
C              Evaluate the residual.
C
               LVARMC = .TRUE.
               ICASE = 1
               MARCDBG = .FALSE.
               CALL FUN (LBURNR, LENRGY, LMULTI, LVCOR, LTDIF, LVARMC,
     1                   LTIME, WT, EPS, XGIVEN, TGIVEN, X, SN, BUFFER,
     3                   SCRTCH(1, 1), YV, SCRTCH(1, 2), SCRTCH(1, 3),
     4                   SCRTCH(1, 4), COND, D, DKJ, TDR, ICKWRK,
     5                   RCKWRK, IMCWRK, RMCWRK, F, SCRTCH(1, 5),
     6                   SSAVE, RKFT, RKRT, ICASE, RSAVE, MARCDBG)
               CALL CKCOPY (NATJ * JJ, F, BUFFER)
C
            ELSEIF (SIGNAL .EQ. 'PREPARE') THEN
C              Prepare the Jacobian matrix.
C
               ICASE = 2
               LVARMC = .TRUE.
               RETURN = .FALSE.
C
2020           CONTINUE
               ERROR = .FALSE.
               CALL TWPREP (ERROR, LOUT, A, IASIZE, BUFFER, NATJ, 
     1                      CONDIT, IGRPA, IGRPB, IPIVOT, JJ, RETURN)
               KERR = KERR.OR.ERROR
               IF (KERR) RETURN
C
               MARCDBG = .FALSE.
               IF (RETURN) THEN
                  CALL FUN (LBURNR, LENRGY, LMULTI, LVCOR, LTDIF,
     1                      LVARMC, LTIME, WT, EPS, XGIVEN, TGIVEN,
     3                      X, SN, BUFFER, SCRTCH(1, 1), YV,
     4                      SCRTCH(1, 2), SCRTCH(1, 3), SCRTCH(1, 4),
     5                      COND, D, DKJ, TDR, ICKWRK, RCKWRK, IMCWRK,
     6                      RMCWRK, F, SCRTCH(1, 5), SSAVE, RKFT,
     7                      RKRT, ICASE, RSAVE, MARCDBG)
C
                  ICASE = 3
                  LVARMC = .FALSE.
                  CALL CKCOPY (NATJ * JJ, F, BUFFER)
                  GO TO 2020
               ENDIF
C
            ELSEIF (SIGNAL .EQ. 'SOLVE') THEN
C              Solve the linear equations.
C
               CALL TWSOLV (ERROR, LOUT, A,IASIZE, BUFFER, NATJ,
     1                      IGRPA, IGRPB, IPIVOT, JJ)
               KERR = KERR.OR.ERROR
               IF (KERR) RETURN
C
            ELSEIF (SIGNAL .EQ. 'RETAIN') THEN
C              Retain the solution for time integration
C
               CALL CKCOPY (NATJ * JJ, BUFFER, SN)
               IF (SPOS .GE. 0.0) THEN
                  DO 2030 J = 1, JJ
                     CALL PRPOS (KK, SN(NYS+1,J), SPOS)
2030              CONTINUE
               ENDIF
C
            ELSEIF (SIGNAL .EQ. 'SHOW') THEN
C              Show the solution

c     HACK
               print *,'showing...'
               ICASE = 1
               LVARMC = .TRUE.
               MARCDBG = .TRUE.
               LENRGY = .TRUE.
               CALL FUN (LBURNR, LENRGY, LMULTI, LVCOR, LTDIF,
     1                      LVARMC, LTIME, WT, EPS, XGIVEN, TGIVEN,
     3                      X, SN, BUFFER, SCRTCH(1, 1), YV,
     4                      SCRTCH(1, 2), SCRTCH(1, 3), SCRTCH(1, 4),
     5                      COND, D, DKJ, TDR, ICKWRK, RCKWRK, IMCWRK,
     6                      RMCWRK, F, SCRTCH(1, 5), SSAVE, RKFT,
     7                      RKRT, ICASE, RSAVE, MARCDBG)
c     HACK
               CALL PRINT (LOUT, LENRGY, LMOLE, LMULTI, LTDIF, LVCOR, 
     +              X, BUFFER, SCRTCH(1, 1), SCRTCH(1, 2), KSYM, WT,
     +              cckwrk, ICKWRK, RCKWRK, IMCWRK, RMCWRK)
C
            ELSEIF (SIGNAL .EQ. 'SAVE') THEN

C              Save the solution
               REWIND LRCRVR
               CALL PRSAVE (ICKWRK, RCKWRK, CCKWRK,
     1                      IMCWRK, RMCWRK, LOUT, LRCRVR)
               WRITE (LRCRVR) ISOLUT
               WRITE (LRCRVR) NATJ, JJ, P, BUFFER(NM, 1)
               WRITE (LRCRVR) (X(J), J = 1, JJ)
               WRITE (LRCRVR) ((BUFFER(N, J), N=1,NATJ), J=1,JJ)
C
            ELSEIF (SIGNAL .EQ. 'UPDATE') THEN
C              Update the solution to a new grid.
               IF (.NOT. LENRGY) THEN
                  DO 2040 J = 1, JJ
                     BUFFER(NT, J) = CKBSEC (NTEMP,X(J),XGIVEN,TGIVEN)
2040              CONTINUE
               ENDIF
               IF (.NOT. LBURNR) THEN
                  DO 2050 J = 1, JJ
                     IF (X(J) .EQ. XFIXT) JFIXT = J
2050              CONTINUE
               ENDIF
               DO 2060 J = 2, JJ
                  BUFFER(NM, J) = BUFFER(NM, 1)
2060           CONTINUE
C
            ENDIF
C           Bottom of the block to service requests from TWOPNT.
C
            GO TO 2010
         ENDIF
C
C///  TWOPNT FINISHES.
C
2070  CONTINUE

c     just print the solution

      if (.not. lsolve) then
               print *,'showing...'
               MARCDBG = .TRUE.
               CALL FUN (LBURNR, LENRGY, LMULTI, LVCOR, LTDIF,
     1                      LVARMC, LTIME, WT, EPS, XGIVEN, TGIVEN,
     3                      X, SN, S, SCRTCH(1, 1), YV,
     4                      SCRTCH(1, 2), SCRTCH(1, 3), SCRTCH(1, 4),
     5                      COND, D, DKJ, TDR, ICKWRK, RCKWRK, IMCWRK,
     6                      RMCWRK, F, SCRTCH(1, 5), SSAVE, RKFT,
     7                      RKRT, ICASE, RSAVE, MARCDBG)
c     HACK
               CALL PRINT (LOUT, LENRGY, LMOLE, LMULTI, LTDIF, LVCOR, 
     +              X, s, SCRTCH(1, 1), SCRTCH(1, 2), KSYM, WT, 
     +              cckwrk, ICKWRK, RCKWRK, IMCWRK, RMCWRK)
               write (lout, "( )")
      end if
C
C     end of SUBROUTINE PRETWO
      RETURN
      END
C
      SUBROUTINE PRINT 
     +   (LOUT, LENRGY, LMOLE, LMULTI, LTDIF, LVCOR, 
     +   X, S, Y, WDOT, KSYM, WT,
     +   cckwrk, ICKWRK, RCKWRK, IMCWRK, RMCWRK)
C
C  START PROLOGUE
C  Print the latest solution.
C
C  LOUT      - integer scalar, formatted file output unit number.
C  LMOLE     - logical, TRUE for input/output in mole fraction,
C                      FALSE for input/output in mass fraction.
C  X(*)      - real array, mesh point locations (cm)
C  S(*,*)    - real matrix, depdent variables, where at X(J),
C                 S(NT,J) is temperature,
C                 S(NM,J) is flow rate, and
C                 S(NYS+K,J) is species K mass fraction.
C  Y(*)      - real array, species mass fractions.
C  KSYM(*)   - character string array, species names.
c  cckwrk(*) - character array, chemkin workspace.
C  ICKWRK(*) - integer array, CHEMKIN workspace.
C  RCKWRK(*) - real array, CHEMKIN workspace.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      double precision 
     +   AREA, cond, d, diffusionCoefficients, dkj, elementFlux, hml, 
     +   molarDensity, moleFractions,
     +   negativeIntegral, p, positiveIntegral,
     +   rho, sum, tav, tdr, thermalConductivity, thisarea, 
     +   wdot, xav, xmf, mixtureViscosity,
     +   xmfp, yav, yv, z, z0
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
c      real 
C     +   AREA, cond, d, diffusionCoefficients, dkj, elementFlux, hml, 
C     +   molarDensity, moleFractions,
C     +   negativeIntegral, p, positiveIntegral,
C     +   rho, sum, tav, tdr, thermalConductivity, thisarea, xav, xmf, 
C     +   xmfp, yav, yv, z, z0, mixtureViscosity
C*****END precision > single
C
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
      LOGICAL
     +   kerr, lewis, lfuel, LMOLE, lpath, lrrop, lsolve, latom, ldflx, 
     +   LMULTI, LENRGY, LTDIF, LVCOR, lheat, ldiff

      common / prlcon / 
     +   latom, ldflx, ldiff, lewis, lfuel, lheat, lpath, lrrop, lsolve
C
      allocatable cond
      allocatable d
      allocatable diffusionCoefficients
      allocatable dkj
      allocatable elementConcentration
      allocatable elementName
      allocatable hml
      allocatable ki
      allocatable labels
      allocatable moleFractions
      allocatable negativeIntegral
      allocatable positiveIntegral
      allocatable quantity
      allocatable reactionRates
      allocatable speciesConcentration
      allocatable tdr
      allocatable values
      allocatable xav
      allocatable xmf
      allocatable xmfp
      allocatable yav
      allocatable yv
      allocatable z
      allocatable z0
      DIMENSION 
     +   CCKWRK(LENCCK), cond(:), d(:,:), diffusionCoefficients(:),
     +   dkj(:,:,:), 
     +   elementConcentration(:), elementName(:), hml(:), ki(:),
     +   ICKWRK(LENICK), labels(:), moleFractions(:),
     +   negativeIntegral(:),
     +   positiveIntegral(:), quantity(:, :), RCKWRK(LENRCK),
     +   reactionRates(:), S(NATJ, JJ), speciesConcentration(:), 
     +   tdr(:,:), values(:), wdot(kk), X(JJ), Y(KK), z(:), z0(:), 
     +   WT(KK), xmf(:), xmfp(:), xav(:), yav(:), yv(:,:)
      CHARACTER 
     +   CCKWRK*(*), elementName*20, KSYM(KK)*(*), labels*20,
     +   string*100, values*20
      integer 
     +   element, elements, field, fields, itemp1, itemp2, itemp3,
     +   ndim, nspec, quantity
      intrinsic 
     +   len_trim

c     get the molecular weights (in case they have not been gotten)

      CALL CKWT (ICKWRK, RCKWRK, WT)

c     count the output fields and prepare to write them

      call ckindx (ickwrk, rckwrk, elements, itemp1, itemp2, itemp3)

      field = 0
      fields = 0
c     always include x, t, v, rho
      fields = fields + 4
c     heat release
      if (lheat) fields = fields + 1
c     fuel consumption
      if (lfuel) fields = fields + 1
c     molar fractions of atoms and percent change thereof
      if (latom) fields = fields + elements + elements
c     always include species
      fields = fields + kk
c     diffusion flux
      if (ldflx) fields = fields + 1 + kk 
      if (ldflx .and. latom) fields = fields + elements
c     reaction rates of progress
      if (lrrop) fields = fields + ii
c     lewis numbers
c     includes heat capacity, mixture viscosity, and Prandtl number
      if (lewis) fields = fields + 3 + kk
c     diffusivities
      if (ldiff) fields = fields + 1 + kk

      allocate (labels(fields))

c     labels for x, t, v, rho

      field = field + 1
      labels(field) = 'X                   '
      field = field + 1
      labels(field) = 'T                   '
      field = field + 1
      labels(field) = 'V                   '
      field = field + 1
      labels(field) = 'RHO                 '

c     label for heat release

      if (lheat) then
         field = field + 1
         labels(field) = 'HEAT_RELEASE        '
      end if

c     label for fuel consumption

      if (lfuel) then
         field = field + 1
         labels(field) = 'FUEL_CONS           '
      end if

c     labels for molar fractions of atoms and percent change thereof

      if (latom) then
         allocate (quantity(elements, kk))
         allocate (speciesConcentration(kk))
         allocate (elementConcentration(elements))
         allocate (elementName(elements))
         allocate (z(elements))
         allocate (z0(elements))

         call cksyme (cckwrk, lout, elementName, kerr)
         if (kerr) stop
         call ckncf (elements, ickwrk, rckwrd, quantity)

         do element = 1, elements
            field = field + 1
            labels(field) = 
     +         'z(' // elementName(element)(1 : 
     +         len_trim (elementName(element))) // ')'
         end do
         do element = 1, elements
            field = field + 1
            labels(field) = 
     +         'change_z(' // elementName(element)(1 : 
     +         len_trim (elementName(element))) // ')'
         end do
      end if

c     labels for species

      do k = 1, kk
         field = field + 1
         labels(field) = ksym(k)
      end do

c     labels for reaction rates of progress

      if (lrrop) then
         do i = 1, ii
            field = field + 1
            write (labels(field), "('Rxn', i3.3)") i
         end do
      end if

c     labels for diffusion flux

      if (ldflx) then
         field = field + 1
         labels(field) = 'X-half              '
         do k = 1, kk
            field = field + 1
            labels(field) = 
     +         'dflx(' // ksym(k) (1 : len_trim (ksym(k))) // ')'
         end do

         if (latom) then
            do element = 1, elements
               field = field + 1
               labels(field) = 
     +              'dflx(elem_' // elementName(element)(1 : 
     +              len_trim (elementName(element))) // ')'
            end do
         end if
      end if

c     labels for Lewis number
c     (includes mixture viscosity and Prandtl number)

      if (lewis) then
         field = field + 1
         labels(field) = '     cp'
         field = field + 1
         labels(field) = 'mixvisc'
         field = field + 1
         labels(field) = 'Prandtl'
         do k = 1, kk
            field = field + 1
            labels(field) = 
     +         'lewis(' // ksym(k) (1 : len_trim (ksym(k))) // ')'
         end do
      end if

c     labels for diffusivities

      if (ldiff) then
         field = field + 1
         labels(field) = 'diff(T)'
         do k = 1, kk
            field = field + 1
            labels(field) = 
     +         'diff(' // ksym(k) (1 : len_trim (ksym(k))) // ')'
         end do
      end if

c     write the labels

      write (string, "('(/8x,', i10, '(a20))')") fields
      write (LOUT, string) labels

c     prepare to evalute reaction rates

      if (lpath .or. lrrop) then
         allocate (reactionRates(ii))
      end if

      if (lpath) then
         allocate (negativeIntegral(ii))
         allocate (positiveIntegral(ii))

         do i = 1, ii
            negativeIntegral(i) = 0
            positiveIntegral(i) = 0
         end do
      end if

c     values for diffusion flux

      if (ldflx) then
         allocate (cond(jj))
         allocate (d(kk, jj))
         allocate (dkj(kk, kk, jj))
         allocate (tdr(kk, jj))
         allocate (xav(kk))
         allocate (xmf(kk))
         allocate (xmfp(kk))
         allocate (yav(kk))
         allocate (yv(kk, jj))

         IF (LMULTI) THEN
            CALL MCMULT (LENRGY, LTDIF, X, S, yav, ICKWRK, RCKWRK,
     1                   IMCWRK, RMCWRK, WT, xmf, xmfp, xav, cond, d,
     1                   tdr, dkj)
         ELSE
            CALL MCMIX (LENRGY, LTDIF, X, S, yav, ICKWRK, RCKWRK,
     1                  IMCWRK, RMCWRK, WT, xmf, xmfp, xav, cond, d,
     2                  tdr, dkj)
         ENDIF

         CALL MDIFV (LMULTI, LVCOR, LTDIF, X, S, WT, yav, xmf, xmfp,
     1        d, tdr, ICKWRK, RCKWRK, yv, dkj)

         do j = 1, jj - 1
            tav = 0.5 * (S(NT, J) + S(NT, J + 11))
            CALL CKAVG (KK, S(NYS + 1, J), S(NYS + 1, J + 1), yav)
            CALL CKRHOY (p, tav, yav, ICKWRK, RCKWRK, rho)
            thisarea = AREA((X(J) + X(J + 1)) / 2)
            do k = 1, kk
C              comments in the code suggest yv is already scaled by rho
C              but the comment block for MDIFV suggests yv is just the
C              diffusion velocity.  the following computes a molar flux.
               yv(k, j) = yv(k, j) * thisarea * rho / WT(k)
            end do
         end do
         do k = 1, kk
            yv(k, jj) = yv(k, jj - 1)
         end do

         deallocate (cond)
         deallocate (d)
         deallocate (dkj)
         deallocate (tdr)
         deallocate (xav)
         deallocate (xmf)
         deallocate (xmfp)
         deallocate (yav)
      end if

c     top of the loop to write the data

      allocate (values(fields))
      write (string, "('(8x,', i10, '(1x, a20))')") fields

      DO J = 1, JJ
         field = 0

c     values for x, t, v, and rho

      CALL CKRHOY (P, S(NT, J), S(NY,J), ICKWRK, RCKWRK, RHO)
      FLMSPD = S(NM, J) / (RHO * AREA(X(J)))

      field = field + 1
      write (values(field), "(1pe20.12)") X(J)
      field = field + 1
      write (values(field), "(1pe20.12)") S(NT, J)
      field = field + 1
      write (values(field), "(1pe20.12)") FLMSPD
      field = field + 1
      write (values(field), "(1pe20.12)") RHO

c     values for heat release

      if (lheat) then
         ndim = 100
         allocate (hml(kk))

         call CKHML (S(NT, J), ICKWRK, RCKWRK, hml)
         call CKWYP (P, S(NT, J), S(NY, J), ICKWRK, RCKWRK, wdot)

         sum = 0
         do k = 1, kk
c           from CKWYP wdot is mole / (cm**3 s)
c           from CKHML hml  is erg / mole
c           sum = sum - wdot(k) * hml(k) * wt(k)
            sum = sum - wdot(k) * hml(k)
         end do

         field = field + 1
         write (values(field), "(1pe20.12)") sum

         deallocate (hml)
      end if

c     values for fuel consumption (gr / cm**3)

      if (lfuel) then
         CALL CKWYP (P, S(NT,J), S(NY,J), ICKWRK, RCKWRK, WDOT)
         field = field + 1
         write (values(field), "(1pe20.12)") - wdot(ifuel) * WT(ifuel)
      end if

c     values for molar concentrations of atoms

      if (latom) then
         call ckytcp 
     +      (p, s(nt, j), s(ny, j), ickwrk, rckwrk, 
     +      speciesConcentration)
         do element = 1, elements
            elementConcentration(element) = 0
            do k = 1, kk
               elementConcentration(element) =
     +            elementConcentration(element) +
     +            speciesConcentration(k) * quantity(element, k)
            end do
         end do
         molarDensity = 0
         do element = 1, elements
            molarDensity = molarDensity + elementConcentration(element)
         end do
         do element = 1, elements
            z(element) = elementConcentration(element) / molarDensity
            if (j .eq. 1) then
               z0(element) = z(element)
            end if
         end do
         do element = 1, elements
            field = field + 1
            write (values(field), "(1pe20.12)") z(element)
         end do
         do element = 1, elements
            field = field + 1
            if (z0(element) .ne. 0) then
               write (values(field), "(1pe20.12)") 
     +              (z(element) / z0(element) - 1) * 100
            else
               write (values(field), "(1pe20.12)") 0
            end if
         end do
      end if

c     values for species

      IF (LMOLE) THEN
         CALL CKYTX (S(NY, J), ICKWRK, RCKWRK, Y)
      ELSE
         CALL CKCOPY (KK, S(NYS+1,J), Y)
      ENDIF
      do k = 1, kk
         field = field + 1
         write (values(field), "(1pe20.12)") Y(k)
      end do

c     values for diffusion flux

      if (ldflx) then
         if (j .lt. jj) then
            field = field + 1
            write (values(field), "(1pe20.12)") (X(J) + X(J + 1)) / 2
            do k = 1, kk
               field = field + 1
               write (values(field), "(1pe20.12)") yv(k, j)
            end do
         else
            field = field + 1
            write (values(field), "(1pe20.12)") X(J)
            do k = 1, kk
               field = field + 1
               write (values(field), "(1pe20.12)") yv(k, j - 1)
            end do
         end if

         if (latom) then
            do element = 1, elements
               elementFlux = 0
               do k = 1, kk
                  elementFlux =
     +                 elementFlux + yv(k, j) * quantity(element, k)
               end do
               field = field + 1
               write (values(field), "(1pe20.12)") elementFlux
            end do
         end if
      end if

c     values for reaction rates of progress

      if (lpath .or. lrrop) then
         call ckqyp 
     +        (p, S(NT, J), S(NY, J), ICKWRK, RCKWRK, reactionRates)
      end if

      if (lrrop) then
         do i = 1, ii
            field = field + 1
            write (values(field), "(1pe20.12)") reactionRates(i)
         end do
      end if

c     values for Lewis numbers

      if (lewis) then
         allocate (diffusionCoefficients(kk))
         allocate (moleFractions(kk))

         call CKYTX (S(NY, J), ICKWRK, RCKWRK, moleFractions)
         call MCAVIS 
     +        (S(NT, J), moleFractions, RMCWRK, mixtureViscosity) 
         call MCACON 
     +        (S(NT, J), moleFractions, RMCWRK, thermalConductivity)
         call MCADIF 
     +        (P, S(NT, J), moleFractions, RMCWRK, 
     +        diffusionCoefficients)
         CALL CKRHOY (P, S(NT, J), S(NY, J), ICKWRK, RCKWRK, RHO)
         call CKCPBS (S(NT, J), S(NY, J), ICKWRK, RCKWRK, cp)

c        from MCACON    themalconductivity is erg / cm K s
c        from CKRHOY                   rho is g / cm**3
c        from CKCPBS                    cp is erg / g K
c        from MCADIF diffusionCoefficients is cm**2 /s
c        from MCAVIS      mixtureViscosity is g / cm s
         field = field + 1
         write (values(field), "(1pe20.12)") 
     +      cp
         field = field + 1
         write (values(field), "(1pe20.12)") 
     +      mixtureViscosity
         field = field + 1
         write (values(field), "(1pe20.12)") 
     +      mixtureViscosity * cp / thermalconductivity
         do k = 1, kk
            field = field + 1
            write (values(field), "(1pe20.12)") 
     +           thermalConductivity 
     +           / (RHO * cp * diffusionCoefficients(k))
         end do

         deallocate (diffusionCoefficients)
         deallocate (moleFractions)
      end if

c     values for diffusivities

      if (ldiff) then
         allocate (diffusionCoefficients(kk))
         allocate (moleFractions(kk))

         call CKYTX (S(NY, J), ICKWRK, RCKWRK, moleFractions)
         call MCACON 
     +        (S(NT, J), moleFractions, RMCWRK, thermalConductivity)
         call MCADIF 
     +        (P, S(NT, J), moleFractions, RMCWRK, 
     +        diffusionCoefficients)
         CALL CKRHOY (P, S(NT, J), S(NY, J), ICKWRK, RCKWRK, RHO)
         call CKCPBS (S(NT, J), S(NY, J), ICKWRK, RCKWRK, cp)

c        from MCACON    themalconductivity is erg / cm K s
c        from CKRHOY                   rho is g / cm**3
c        from CKCPBS                    cp is erg / g K
c        from MCADIF diffusionCoefficients is cm**2 /s
         field = field + 1
         write (values(field), "(1pe20.12)") 
     +      thermalConductivity / (RHO * cp)
         do k = 1, kk
            field = field + 1
            write (values(field), "(1pe20.12)") 
     +         diffusionCoefficients(k)
         end do

         deallocate (diffusionCoefficients)
         deallocate (moleFractions)
      end if

c     write the values

      write (string, "('(i8,', i10, '(a20))')") fields
      write (LOUT, string) j, values

c     gather reaction data for drawing path diagrams.

      if (lpath) then
         weight = 0
         if (j .gt. 1) then
            weight = weight + 0.5 * (x(j) - x(j - 1))
         else if (j .lt. jj) then
            weight = weight + 0.5 * (x(j + 1) - x(j))
         end if
         do i = 1, ii
            if (reactionRates(i) .le. 0) then
               negativeIntegral(i) = 
     +              negativeIntegral(i) + weight * reactionRates(i)
            else
               positiveIntegral(i) = 
     +              positiveIntegral(i) + weight * reactionRates(i)
            end if
         end do
      end if

c     bottom of the loop to write the data

      end do

c     Write reaction data to path.dat for drawing path diagrams.

      if (lpath) then
         open (unit = 73, file = "path.dat")
         write (73, '(1x, 1pe20.12)') 
     +      (negativeIntegral(i), i = 1, ii),
     +      (positiveIntegral(i), i = 1, ii)
         close (unit = 73)
      end if

c     dellocate before exiting.

      if (latom) deallocate (elementConcentration)
      if (latom) deallocate (elementName)
      deallocate (labels)
      if (lpath) deallocate (negativeIntegral)
      if (lpath) deallocate (positiveIntegral)
      if (latom) deallocate (quantity)
      if (lpath .or. lrrop) deallocate (reactionRates)
      if (latom) deallocate (speciesConcentration)
      if (ldflx) deallocate (yv)
      if (latom) deallocate (z)
      if (latom) deallocate (z0)
      deallocate (values)

C     end of SUBROUTINE PRINT
      RETURN
      END
C
      SUBROUTINE PRREST (JMAX, LREST, LOUT, MFILE, X, S, KERR)
C
C  START PROLOGUE
C  Fill mesh with a previous solution from a restart file.
C
C  JMAX      - integer scalar, maximum gridpoints
C  LREST     - integer scalar, binary solution input file unit number.
C  LOUT      - integer scalar, formatted output file unit number.
C  MFILE     - integer scalar, index of solution on restart file.
C  X(*)      - real array, mesh point locations (cm)
C  S(*,*)    - real matrix, dependent variables;
C                           S(NT,J) = T(J), temperature
C                           S(NYS+K,J) = Y(K,J), mass fractions
C                           S(NM,J) = FLRT(J), flow rate
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
C
      DIMENSION X(JMAX), S(NATJ,JMAX)
      CHARACTER*16 ICHR, ICKLNK, IMCLNK, ISOLUT, ISENSI, IHSENS
      PARAMETER (ICKLNK='CKLINK', IMCLNK='MCLINK', ISOLUT='SOLUTION',
     1           ISENSI='SENSITIVITY', IHSENS='HSENSITIVITY')
      LOGICAL KERR
      KERR = .FALSE.
C
C     Read the restart data file for starting profiles.
      NREST = 0
1050  CONTINUE
      READ (LREST, END = 1000, ERR = 1000) ICHR
      IF (ICHR .EQ. ICKLNK) THEN
         DO 1060 L = 1, 4
            READ (LREST, END = 1000, ERR = 1000)
1060     CONTINUE
      ELSEIF (ICHR .EQ. IMCLNK) THEN
         DO 1070 L = 1, 3
            READ (LREST, END = 1000, ERR = 1000)
1070     CONTINUE
      ELSEIF (ICHR .EQ. ISOLUT) THEN
         READ (LREST, END = 1000, ERR = 1000) NNNN, JJ, P, DMFLRT
         READ (LREST, END = 1000, ERR = 1000) (X(J), J = 1, JJ)
         READ (LREST, END = 1000, ERR = 1000)
     +      ((S(N, J), N = 1, NNNN), J = 1, JJ)
         IF (NNNN .NE. NATJ) THEN
            WRITE (LOUT, *) ' FATAL ERROR, INCOMPATIBLE RESTART FILE'
            KERR = .TRUE.
            RETURN
         ENDIF
         NREST = NREST + 1
         IF (NREST .EQ. MFILE) GO TO 1000
      ELSEIF (ICHR .EQ. ISENSI) THEN
         DO 1080 I = 1, II
            READ (LREST, END = 1000, ERR = 1000)
1080     CONTINUE
      ELSEIF (ICHR .EQ. IHSENS) THEN
         DO 1090 K = 1, KK
            READ (LREST, END = 1000, ERR = 1000)
1090     CONTINUE
      ELSE
         WRITE (LOUT, *)
     1   'FATAL ERROR, NOT A SOLUTION ON RESTART FILE'
         KERR = .TRUE.
         RETURN
      ENDIF
      GO TO 1050
C
1000  CONTINUE
      IF (NREST .NE. MFILE) THEN
         WRITE (LOUT, *) ' Error reading solution file...'
         KERR = .TRUE.
      ENDIF
C
C     end of SUBROUTINE PRREST
      RETURN
      END
C
      SUBROUTINE PRSAVE (ICKWRK, RCKWRK, CCKWRK,
     1                   IMCWRK, RMCWRK, LOUT, LSAVE)
C
C  START PROLOGUE
C  Write CHEMKIN and TRANSPORT work arrays to binary file.
C
C  ICKWRK(*) -  integer array, CHEMKIN workspace
C  RCKWRK(*) -  real array, CHEMKIN workspace
C  CCKWRK(*) -  character-string array, CHEMKIN workspace
C  IMCWRK(*) -  integer array, TRANSPORT workspace
C  RMCWRK(*) -  real array, TRANSPORT workspace
C  LOUT      -  integer scalar, formatted output file unit number
C  LSAVE     -  integer scalar, binary output file unit number
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
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
C
      DIMENSION ICKWRK(LENICK), RCKWRK(LENRCK),
     1          IMCWRK(LENIMC), RMCWRK(LENRMC)
      CHARACTER*(*) CCKWRK(LENCCK)
      CHARACTER*16 ILINK
C
      ILINK = 'CKLINK          '
      WRITE (LSAVE) ILINK
      CALL CKSAVE (LOUT, LSAVE, ICKWRK, RCKWRK, CCKWRK)
C
      ILINK = 'MCLINK          '
      WRITE (LSAVE) ILINK
      CALL MCSAVE (LOUT, LSAVE, IMCWRK, RMCWRK)
C
C     end of SUBROUTINE PRSAVE
      RETURN
      END
C
      SUBROUTINE PRSENS (ISENS, LBURNR, LENRGY, LMULTI, LTDIF, LMOLE,
     1                   LSAVE, LRCRVR, LOUT, LVARMC, LTIME, WT, EPS,
     3                   XGIVEN, TGIVEN, X, SN, S, SCRTCH, YV, COND,
     4                   D, DKJ, TDR, ICKWRK, RCKWRK, IMCWRK, RMCWRK,
     5                   F, FN, DS, XA, A6, IPIVOT, BUFFER, SSAVE,
     6                   RKFT, RKRT, RSAVE, KERR)
C
C  START PROLOGUE
C
C  ISENS      - character-string, type of sensitivity analysis
C  LBURNR     - logical, T  burner stabilized flame,
C                        F, freely propagating adiabatic flame
C  LENRGY     - logical, T  solve energy equation,
C                        F, use fixed temperature profile
C  LMULTI     - logical, T, use multicomponent formulas,
C                        F, use mixture-averaged formulas
C  LTDIF      - logical, T, evaluate thermal diffusion ratios as
C                           well as diffusion coefficients
C  LMOLE     - logical, TRUE for input/output in mole fraction,
C                      FALSE for input/output in mass fraction.
C  LSAVE     -  integer scalar, binary output file unit number
C  LRCRVR     - integer scalar, binary output scratch file unit number
C  LOUT       - integer scalar, formatted output file unit number
C  LVARMC     - logical, T, compute new transport properties
C                        F, use previously stored values
C  LTIME      - logical, T, time step of DT will be added to residual
C  WT(*)      - real array, species molecular weights (gm/mole)
C  EPS(*)     - real array, mass flux fractions at the burner
C  XGIVEN(*)  - real array, distances from the burner for TGIVEN
C  TGIVEN(*)  - real array, temperatures at given distances (K)
C  X(*)
C  SN(*,*)    - real matrix, dependent variables at previous timestep
C  S(*,*)     - real matrix, dependent variables;
C                            S(NT,J) = T(J), temperature
C                            S(NYS+K,J) = Y(K,J), mass fractions
C                            S(NM,J) = FLRT(J), flow rate
C  SCRTCH(*,*)
C  YV(*,*)    - real matrix, mass fractions time diffusion velocities at
C               mesh midpoints;
C               YV(K,J) is flux of species K between J and J+1.
C  COND(*)    - real array, thermal conductivities at mesh midpoints
C               (erg/cm-K-sec);
C               if LVARMC=.TRUE., these are computed each time function
C               is called, otherwise the stored values are used
C  D(*,*)    - real matrix, species diffusion coefficients at mesh
C              midpoints (cm^2/sec);
C              if LVARMC=.TRUE., these are computed each time function
C              is called, otherwise the stored values are used.
C  DKJ(*,*,*)- real 3-dimensional array, multicomponent diffusion
C              coefficients at mesh midpoints (cm^2/sec)
C  TDR(*,*)  - real matrix, species thermal diffusion ratios at mesh
C              midpoints;
C              if LVARMC=.TRUE., these are computed each time function
C              is called, otherwise the stored values are used.
C  ICKWRK(*) - integer array, CHEMKIN workspace
C  RCKWRK(*) - real array, CHEMKIN workspace
C  IMCWRK(*) - integer array, TRANSPORT workspace
C  RMCWRK(*) - integer array, TRANSPORT workspace
C  F(*,*)    - real matrix, governing equations residuals of S(*,*);
C              F(NT,J) is the energy equation residual,
C              F(NM,J) is the mass equation residual,
C              F(NYS+K,J) is the species equation for species K
C  VN(*,*)
C  DS(*)
C  XA(*)
C  A6(*)
C  IPIVOT(*)
C  BUFFER(*,*)
C  SSAVE(*)  - real array, temperature values at which stored reaction
C              rates were evaluated.
C  RKFT(*)   - real array, forward reaction rates at SSAVE.
C  RKRT(*)   - real array, reverse reaction rates at SSAVE.
C  ICASE     - integer scalar, flag for production rate evaluation;
C              =1, CALL CKWYP for production rates,
C              =2, save temperatures in SSAVE, CALL CKKFRT for RKFT
C                  and RKRT, and CALL CKWYPK for production rates,
C              =3, compare new temperature values to SSAVE, and if
C                  same, CALL CKKFRT as in (2), else do (1)
C  RSAVE(*,*)- real matrix, for ICASE=2, save species production rates,
C              for ICASE=3, use RSAVE for species production rates
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
      CHARACTER ISENS*(*)
C
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
C
      LOGICAL
     +   ERROR, LBURNR, LENRGY, LMOLE, LMULTI, LTDIF, LTIME, LVARMC,
     +   LVCOR, RETURN, KERR
C     Integer arrays
      DIMENSION ICKWRK(LENICK), IMCWRK(LENIMC), IPIVOT(2000)
C     Real arrays
      DIMENSION BUFFER(NATJ,JJ), COND(JJ), D(KK,JJ), DKJ(KK,KK,JJ),
     1          DS(JJ), EPS(KK), F(NATJ,JJ), FN(NATJ,JJ),
     3          RCKWRK(LENRCK), RKFT(II,JJ), RKRT(II,JJ),
     4          RMCWRK(LENRMC), S(NATJ,JJ), SCRTCH(KK,JJ),
     5          SN(NATJ,JJ), TDR(KK,JJ), TGIVEN(NTEMP), SSAVE(JJ),
     6          A6(NTR), WT(KK), X(JJ), XA(IASIZE), XGIVEN(NTEMP),
     7          YV(KK,JJ), RSAVE(KK,JJ)

      LOGICAL MARCDBG
C
C///////////////////////////////////////////////////////////////////////
C
C     EVALUATE THE JACOBIAN MATRIX.
C
C///////////////////////////////////////////////////////////////////////
C
      KERR = .FALSE.
      CALL CKCOPY (NATJ * JJ, S, BUFFER)
C
      ICASE = 2
      IGRPA = 0
      IGRPB = 0
      LTIME = .FALSE.
      LVARMC = .TRUE.
      RETURN = .FALSE.
0100  CONTINUE
C
      ERROR = .FALSE.
      CALL TWPREP (ERROR, LOUT, XA,IASIZE, BUFFER, NATJ, CONDIT,
     1             IGRPA, IGRPB, IPIVOT, JJ, RETURN)
      KERR = KERR.OR.ERROR
      IF (KERR) RETURN
C
      IF (RETURN) THEN
C
         MARCDBG = .FALSE.
         CALL FUN (LBURNR, LENRGY, LMULTI, LVCOR, LTDIF, LVARMC, LTIME,
     1             WT, EPS, XGIVEN, TGIVEN, X, SN, BUFFER, SCRTCH(1, 1),
     3             YV, SCRTCH(1, 2), SCRTCH(1, 3), SCRTCH(1, 4), COND,
     4             D, DKJ, TDR, ICKWRK, RCKWRK, IMCWRK, RMCWRK, F,
     5             SCRTCH(1, 5), SSAVE, RKFT, RKRT, ICASE, RSAVE,
     6             MARCDBG)
C
         ICASE = 3
         LVARMC = .FALSE.
         CALL CKCOPY (NATJ * JJ, F, BUFFER)
         GO TO 0100
      ENDIF
C
C///////////////////////////////////////////////////////////////////////
C
C     EVALUATE THE RESIDUAL AT THE SOLUTION.
C
C///////////////////////////////////////////////////////////////////////
C
      CALL CKCOPY (NATJ * JJ, S, BUFFER)
C
      ICASE = 1
      LVARMC = .TRUE.
C
      MARCDBG = .FALSE.
      CALL FUN (LBURNR, LENRGY, LMULTI, LVCOR, LTDIF, LVARMC, LTIME,
     1          WT, EPS, XGIVEN, TGIVEN, X, SN, BUFFER, SCRTCH(1,1),
     3          YV, SCRTCH(1, 2), SCRTCH(1, 3), SCRTCH(1, 4), COND,
     4          D, DKJ, TDR, ICKWRK, RCKWRK, IMCWRK, RMCWRK, FN,
     5          SCRTCH(1, 5), SSAVE, RKFT, RKRT, ICASE, RSAVE, MARCDBG)
C
C///////////////////////////////////////////////////////////////////////
C
C     EVALUATE SENSITIVITY COEFFICIENTS D(MASS FRACTION) / D(IND)
C
C///////////////////////////////////////////////////////////////////////
C
C     Write the binary solution file header.
      WRITE (LSAVE) ISENS
C
C     Top of the loop over reactions.
      IF (ISENS .EQ. 'SENSITIVITY') THEN
         ITOT = II
      ELSEIF (ISENS .EQ. 'HSENSITIVITY') THEN
         ITOT = KK
      ENDIF
C
      DO 1000 IND = 1, ITOT
C
         IF (ISENS .EQ. 'SENSITIVITY') THEN
C           Perturb the reaction rate
C
            CALL CKRDEX (IND, RCKWRK, SAVEP)
            DP = RELAT * SAVEP + ABSOL
            CALL CKRDEX (- IND, RCKWRK, SAVEP + DP)
         ELSEIF (ISENS .EQ. 'HSENSITIVITY') THEN
C           Perturb the species' specific heat
C
            CALL CKRHEX (IND, RCKWRK, A6)
            DP = RELAT * A6(1) + ABSOL
            DO 0200 L = 1, NTR
               A6(L) = A6(L) + DP
0200        CONTINUE
            CALL CKRHEX (-IND, RCKWRK, A6)
         ENDIF
C
C        Evaluate the perturbed residual.
         ICASE = 1
         MARCDBG = .FALSE.
         CALL FUN (LBURNR, LENRGY, LMULTI, LVCOR, LTDIF, LVARMC, LTIME,
     1             WT, EPS, XGIVEN, TGIVEN, X, SN, S, SCRTCH(1, 1), YV,
     2             SCRTCH(1,2), SCRTCH(1,3), SCRTCH(1, 4), COND, D,
     3             DKJ, TDR, ICKWRK, RCKWRK, IMCWRK, RMCWRK, F,
     4             SCRTCH(1,5), SSAVE, RKFT, RKRT, ICASE, RSAVE,MARCDBG)
C
         IF (ISENS .EQ. 'SENSITIVITY') THEN
C           Restore the reaction rate.
C
            CALL CKRDEX (- IND, RCKWRK, SAVEP)
         ELSE
C           Restore the species' specific heat
C
            DO 0300 L = 1, NTR
               A6(L) = A6(L) - DP
0300        CONTINUE
            CALL CKRHEX (-IND, RCKWRK, A6)
            SAVEP = A6(1)
         ENDIF
C
C///  FORM THE DERIVATIVE OF THE RESIDUAL WITH RESPECT TO THE RATE.
C
         DO 0400 J = 1, JJ
            DO 0400 N = 1, NATJ
            SN(N, J) =  - (F(N, J) - FN(N, J)) / DP
0400     CONTINUE
C
C///  APPLY THE INVERSE OF THE JACOBIAN MATRIX TO OBTAIN THE THE RAW
C///  SENSIVITIES, DS(N, J) / DA(IND).
C
         CALL TWSOLV (ERROR, LOUT, XA,IASIZE, SN, NATJ, IGRPA, IGRPB,
     1                IPIVOT, JJ)
         KERR = KERR.OR. ERROR
         IF (KERR) RETURN
C
C        normalize the sensitivity coefficients
         IF (LMOLE) THEN
            DO 0700 J = 1, JJ
               SUM = 0.0
               DO 0600 L = 1, KK
                  SUM = SUM + SN(NYS+L, J) / WT(L)
0600           CONTINUE
               CALL CKMMWY (S(NY, J), ICKWRK, RCKWRK, WTM)
               WTMSUM = WTM * SUM
               DO 0700 L = NYS+1, NYS+KK
                  SN(L, J) = SAVEP * (SN(L,J)/S(L,J)-WTMSUM)
0700        CONTINUE
         ELSE
            DO 0800 J = 1, JJ
               DO 0800 L = NYS+1, NYS+KK
                  SN(L, J) = SAVEP * SN(L,J) / S(L,J)
0800        CONTINUE
         ENDIF
C
         DO 0900 J = 1, JJ
            SN(NT, J) = SAVEP * SN(NT, J) / S(NT, J)
            SN(NM, J) = SAVEP * SN(NM, J) / S(NM, J)
0900     CONTINUE
C
C        write out the sensitivities
         WRITE (LSAVE) IND, ((SN(N, J), N=1,NATJ), J=1,JJ)
         WRITE (LRCRVR) IND,((SN(N, J), N=1,NATJ), J=1,JJ)
C
C        bottom of the loop over reactions
1000  CONTINUE
C
C     end of SUBROUTINE PRSENS
      RETURN
      END
C
      SUBROUTINE RDKEY (JMAX, LIN, LOUT, KSYM, LBURNR, LMOLE, LUSTGV,
     1                  LENRGY, LMULTI, LVCOR, LTDIF, LUMESH, LRSTRT,
     2                  LCNTUE, MFILE, LASEN, LHSEN, NTOT, X, REAC,
     3                  XINTM, PROD, KR, KI, KP, XX, TT, N1CALL,
     4                  LREGRD, PCTADP, RATCTC, KERR, linflow, tinflow)
C
C  START PROLOGUE
C
C  JMAX     - integer scalar, maximum number of grid points; JMAX
C             is used for dimensioning and for dynamic storage
C             allocation; it can only be changed in the main code.
C  LIN      - integer scalar, formatted input file unit number.
C  KSYM(*)  - character-string array, species names.
C  PATM     - real scalar, pressure of one atmosphere (dynes/cm^2)
C  LBURNR   - logical, .TRUE. for burner stabilized flame problem,
C                      .FALSE. for freely propagating adiabatic flame
C  LTIME    - logical, .TRUE. when time-stepping is used,
C                      .FALSE. when only Newton's method is used
C  LTIME2   - logical, .TRUE. when time-stepping used on each new mesh,
C                      .FALSE. when only Newton's method is used
C  LMOLE    - logical, .TRUE. for mole fraction input and output,
C                      .FALSE. for mass fraction input and output
C  LUSTGV   - logical, .TRUE., use given temperature profile on restart,
C                      .FALSE., use the restart file profile
C  LENRGY   - logical, .TRUE., solve the energy equation,
C                      .FALSE., use XX(*), TT(*) temperature profile
C  LTDIF    - logical, .TRUE., include thermal diffusion,
C                      .FALSE., neglect thermal diffusion
C  LMULTI   - logical, .TRUE., use multicomponent formulas,
C                      .FALSE., use mixture-averaged formulas
C  LVCOR    - logical, .TRUE., use correction velocity formulism,
C                      .FALSE., use 'trace' approximation, lumping
C                               all transport error into final species
C  LUMESH   - logical, .TRUE.,  start on a uniform mesh,
C                      .FALSE., start on specified non-uniform X(*)
C                               specified by GRID keyword
C  LRST     - logical, .TRUE.,  start from a previous profile,
C                      .FALSE., start fresh
C  LCNTUE   - logical, .TRUE.,  return to RDKEY for continuation,
C                      .FALSE., this is the final problem input
C  LASEN    - logical, .TRUE.,  compute sensitivies for reactions,
C                      .FALSE., do not compute
C  LHSEN    - logical, .TRUE.,  compute heat-of-formation sensitivities,
C                      .FALSE., do not compute
C  NTOT     - integer scalar, maximum number of grid points for problem,
C  X(*)     - real array, mesh point locations (cm)
C  REAC(*)  - real array, reactant input mole/mass fractions.
C  XINTM(*) - real array, intermediate input mole/mass fractions.
C  KR(*)    - integer array, reactant input species indices.
C  KI(*)    - integer array, intermediate input species indices.
C  KP(*)    - integer array, product input species indices.
C  XX(*)    - real array, locations for initial temp. profiles (cm)
C  TT(*)    - real array, initial temperatures (K)
C  N1CALL   -
C  LREGRD   -
C  PCTADP   -
C  RATCTC   -
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
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
      common / prlcon / 
     +   latom, ldflx, ldiff, lewis, lfuel, lheat, lpath, lrrop, lsolve

      logical LULNUM
      common / MARC / DFTH, CHTH, LULNUM
C
      LOGICAL
     +   CNTNUD, IERR, KERR, ERROR, LASEN, LBURNR, LCNTUE,
     +   LENRGY, LFIRST, LHSEN, LMOLE, LMULTI,
     +   LRSTRT, LTDIF, LTIME, LTIME2, LUMESH, LUSTGV,
     +   LVCOR, NEC, NOPT, LREGRD, lpath, lrrop, lfuel, latom,
     +   lsolve, linflow, ldflx, lheat, lewis, ldiff
      INTEGER 
     +   CKLSCH, FIRST, IFUEL, LAST
      EXTERNAL CKLSCH
C
      CHARACTER ID*9, KEY*4, CKCHUP*4, KSYM(*)*16, LINE*80, fuel*16
      EXTERNAL CKCHUP
      PARAMETER (ID = ' RDKEY:  ')
C
      DIMENSION KI(KK), KP(KK), KR(KK), NEC(8), NOPT(5), PROD(KK),
     1          REAC(KK), TT(JMAX), VALUE(5), X(JMAX), XINTM(KK),
     2          XX(JMAX)
C
C///////////////////////////////////////////////////////////////////////
C
C     INITIALIZE VARIABLES
C
C///////////////////////////////////////////////////////////////////////
C
      DATA NEC,NOPT,CNTNUD/ 14*.FALSE./
C
      KERR = .FALSE.
      IF (LCNTUE) THEN
         CNTNUD = .TRUE.
         NP = 0
      ELSE
         DO 0100 K = 1, KK
            REAC(K) = 0.
            XINTM(K) = 0.
            PROD(K)  = 0.
0100     CONTINUE
C
         NREAC  = 0
         NINTM  = 0
         NPROD  = 0
         JJ = 6
         NTOT   = JMAX
         SFLR   = - 1.E-3
         DT2    = 1.0E-6
         N1CALL = 1
         GFAC   = 1.0
         DFTH   = 1.0
         CHTH   = 1.0
         LULNUM = .FALSE.
         SPOS   = - 1.0
         MFILE  = 1
         NTEMP  = 0
         NP     = 0
         WNDFAC = 1.0
         XSTR   = 0.0
         PCTADP = 0.75
         RATGTC = 1.0
         TRFAC  = 0.d0
         LUMESH = .TRUE.
         LUSTGV = .FALSE.
         LRSTRT = .FALSE.
         LTDIF  = .FALSE.
         LMULTI = .FALSE.
         LVCOR  = .TRUE.
         LTIME  = .FALSE.
         LTIME2 = .FALSE.
         LENRGY = .FALSE.
         LASEN  = .FALSE.
         LHSEN  = .FALSE.
         LREGRD = .FALSE.
         ifuel  = 0
         lpath  = .false.
         lrrop  = .false.
         latom  = .false.
         ldflx  = .false.
         lewis  = .false.
         ldiff  = .false.
         lfuel  = .false.
         lheat  = .false.
         lsolve = .true.
         linflow = .false.
C
C///  SET VALUES FOR TWOPNT'S CONTROLS FOR DEFAULT
C///  FIXED TEMPERATURE PROBLEM
C
         ERROR = .FALSE.
         CALL TWSETL (ERROR, LOUT, 'ADAPT', .FALSE.)
         CALL TWSETI (ERROR, LOUT, 'LEVELD',   1)
         CALL TWSETI (ERROR, LOUT, 'LEVELM',   1)
         CALL TWSETI (ERROR, LOUT, 'PADD',  JMAX)
         ATOL = 1.0E-9
         CALL TWSETR (ERROR, LOUT, 'SSABS', ATOL)
         CALL TWSETI (ERROR, LOUT, 'SSAGE',   20)
         RTOL = 1.0E-4
         CALL TWSETR (ERROR, LOUT, 'SSREL', RTOL)
         CALL TWSETL (ERROR, LOUT, 'STEADY', .TRUE.)
         CALL TWSETI (ERROR, LOUT, 'STEPS0',   0)
         CALL TWSETI (ERROR, LOUT, 'STEPS1', 100)
         DT1 = 1.0E-6
         CALL TWSETR (ERROR, LOUT, 'STRID0', DT1)
         CALL TWSETR (ERROR, LOUT, 'TDABS', ATOL)
         CALL TWSETI (ERROR, LOUT, 'TDAGE',   10)
         TDEC = 2.2
         CALL TWSETR (ERROR, LOUT, 'TDEC',  TDEC)
         CALL TWSETR (ERROR, LOUT, 'TDREL', RTOL)
         TINC = 2.0
         CALL TWSETR (ERROR, LOUT, 'TINC',  TINC)
         TMAX = 1.0E-4
         CALL TWSETR (ERROR, LOUT, 'TMAX',  TMAX)
         TMIN = 1.0E-10
         CALL TWSETR (ERROR, LOUT, 'TMIN',  TMIN)
         TOLER0 = 1.0E-8
         CALL TWSETR (ERROR, LOUT, 'TOLER0', TOLER0)
         TOLER1 = 0.1
         CALL TWSETR (ERROR, LOUT, 'TOLER1', TOLER1)
         TOLER2 = 0.5
         CALL TWSETR (ERROR, LOUT, 'TOLER2', TOLER2)
C
         KERR = KERR.OR.ERROR
      ENDIF
      LFIRST = .TRUE.
      LCNTUE = .FALSE.
C
C     Top of the loop over input lines
      WRITE (LOUT, '(/10X, A /)') 'KEYWORD INPUT'
C
C     Read next input line
0200  CONTINUE
      KEY = ' '
      LINE = ' '
      IERR  = .FALSE.
      ERROR = .FALSE.
      READ (LIN, '(A)') LINE
      WRITE (LOUT, '(10X, A)') LINE (1 : CKLSCH(LINE))
      CALL CKDTAB (LINE)
      KEY = CKCHUP(LINE(1:4),4)
      LINE(1:4) = ' '
C
C     Comment line if '.' or '/' is first character
      IF (KEY(1 : 1) .EQ. '.' .OR. KEY(1 : 1) .EQ. '/') GO TO 0200
      IND = INDEX(LINE, '(')
      IF (IND .GT. 0) LINE(IND : ) = ' '
C
C///////////////////////////////////////////////////////////////////////
C
C     PROBLEM TYPE KEYWORDS
C
C///////////////////////////////////////////////////////////////////////
C
      IF (KEY .EQ. 'BURN') THEN
C        Burner-stabilized flame problem
         LBURNR = .TRUE.
         NEC(8) = .TRUE.
C
      ELSEIF (KEY .EQ. 'FREE') THEN
C        Adiabatic, freely-propagating flame problem
         LBURNR = .FALSE.
         LENRGY = .TRUE.
         NEC(8) = .TRUE.
C
      ELSEIF (KEY .EQ. 'MOLE') THEN
C        Mole fraction input/output
         NEC(1) = .TRUE.
         LMOLE = .TRUE.
C
      ELSEIF (KEY .EQ. 'MASS') THEN
C        Mass fraction input/output
         NEC(1) = .TRUE.
         LMOLE = .FALSE.
C
      ELSEIF (KEY .EQ. 'TGIV') THEN
C        Energy equation is not included
         NEC(2) = .TRUE.
         LENRGY = .FALSE.
C
      ELSEIF (KEY .EQ. 'ENRG') THEN
C        Energy equation is included
         NEC(2) = .TRUE.
         LENRGY = .TRUE.

      ELSEIF (KEY .EQ. 'NSLV') THEN
C        Do not solve the problem, just print the solution
         lsolve = .false.
C
C///////////////////////////////////////////////////////////////////////
C
C     METHOD OPTIONS KEYWORDS
C
C///////////////////////////////////////////////////////////////////////
C
      ELSEIF (KEY .EQ. 'ATOL') THEN
C        Absolute Newton iteration convergence criteria
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'SSABS', VALUE(1))
C
      ELSEIF (KEY .EQ. 'NJAC') THEN
C        Retirement age of Jacobian during steady-state Newton
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETI (ERROR, LOUT, 'SSAGE', INT(VALUE(1)))
C
      ELSEIF (KEY .EQ. 'RTOL') THEN
C        Relative Newton iteration convergence criteria
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'SSREL', VALUE(1))
C
      ELSEIF (KEY .EQ. 'ATIM') THEN
C        Absolute Newton convergence criteria for timesteps
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TDABS', VALUE(1))
C
      ELSEIF (KEY .EQ. 'RTIM') THEN
C        Relative Newton convergence criteria for timesteps
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TDREL', VALUE(1))
C
      ELSEIF (KEY .EQ. 'TIME') THEN
C        Timestep starting procedure
         LTIME = .TRUE.
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         CALL TWSETI (ERROR, LOUT, 'STEPS1', INT(VALUE(1)))
         CALL TWSETR (ERROR, LOUT, 'STRID0', VALUE(2))
C
      ELSEIF (KEY .EQ. 'ISTP') THEN
C        Number of initial timesteps before Newton
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETI (ERROR, LOUT, 'STEPS0', INT(VALUE(1)))
C
      ELSEIF (KEY .EQ. 'IRET') THEN
C        Retirement age of old time step (default=50)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETI (ERROR, LOUT, 'STEPS2', INT(VALUE(1)))
C
      ELSEIF (KEY .EQ. 'TIM2') THEN
C        Timestepping, after adding the energy equation
         LTIME2 = .TRUE.
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
C        save dt2 for later
         DT2 = VALUE(2)
         NUMDT2 = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'DFAC') THEN
C        Timestep decrease when Newton fails convergence on timestep
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TDEC', VALUE(1))
C
      ELSEIF (KEY .EQ. 'DTMN') THEN
C        Minimum timestep
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TMIN', VALUE(1))
C
      ELSEIF (KEY .EQ. 'DTMX') THEN
C        Maximum timestep
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TMAX', VALUE(1))
C
      ELSEIF (KEY .EQ. 'TJAC') THEN
C        Retirement age of Jacobian during timestepping
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETI (ERROR, LOUT, 'TDAGE', INT(VALUE(1)))
C
C///////////////////////////////////////////////////////////////////////
C
C     GRID PARAMETER KEYWORDS
C
C///////////////////////////////////////////////////////////////////////
C
      ELSEIF (KEY .EQ. 'NPTS') THEN
C        Number of initial meshpoints (this is overwritten "GRID" input0
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         JJ = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'GRID') THEN
C        Initial mesh
         LUMESH = .FALSE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IERR = IERR .OR. .NOT. (NP + 1 .LE. JMAX)
         IF (.NOT. IERR) THEN
            NP = NP + 1
            X(NP) = VALUE(1)
         ENDIF
C
      ELSEIF (KEY .EQ. 'GRAD') THEN
C        Gradient mesh adaption parameter
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TOLER1', VALUE(1))
C
      ELSEIF (KEY .EQ. 'CURV') THEN
C        Curvature mesh adaption parameter
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TOLER2', VALUE(1))
C
      ELSEIF (KEY .EQ. 'XSTR') THEN
C        Point for left boundary condition
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XSTR, IERR)
C
      ELSEIF (KEY .EQ. 'XCEN') THEN
C        Center of mixing region
         NOPT(2) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XCEN, IERR)
C
      ELSEIF (KEY .EQ. 'XEND') THEN
C        Distance at which end boundary condition is applied
         NEC(5) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XEND, IERR)
C
      ELSEIF (KEY .EQ. 'WMIX') THEN
C        Width of mixing zone
         NOPT(1) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, WMIX, IERR)
C
      ELSEIF (KEY .EQ. 'WDIF') THEN
C        Windward differencing
         WNDFAC = 1.0
C
      ELSEIF (KEY .EQ. 'TRFA') THEN
C        if MULT, flux=TRFAC*MC + (1-TRFAC)*MA
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TRFAC, IERR)         
         TRFAC = MIN( 1.d0, MAX(0.d0, TRFAC) )
C
      ELSEIF (KEY .EQ. 'CDIF') THEN
C        Central differencing
         WNDFAC = 0.0
C
      ELSEIF (KEY .EQ. 'SFLR') THEN
C        Floor value for the species bounds
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SFLR, IERR)
C
C///////////////////////////////////////////////////////////////////////
C
C     FLAME DEFINITION KEYWORDS
C
C///////////////////////////////////////////////////////////////////////
C
      ELSEIF (KEY .EQ. 'PRES') THEN
C        Pressure
         NEC(3) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, P, IERR)
         P = P * PATM
C
      ELSEIF (KEY .EQ. 'FLRT') THEN
C        Mass flow rate (gm/sec)
         NEC(4) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, FLRT, IERR)
C
      ELSEIF (KEY .EQ. 'REAC') THEN
C        Reactant
         IF (LFIRST) THEN
            LFIRST = .FALSE.
            NREAC = 0
            DO 0400 K = 1, KK
               REAC(K) = 0.
0400        CONTINUE
         ENDIF
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     +                VALUE, IERR)
         IF (IERR .OR. NREAC + 1 .GT. KK) THEN
            WRITE (LOUT, *) ' ERROR READING DATA FOR KEYWORD ' // KEY
         ELSE
            NREAC = NREAC + 1
            KR(NREAC) = KSPEC
            REAC(KSPEC) = VALUE(1)
         ENDIF
C
      ELSEIF (KEY .EQ. 'INTM') THEN
C        Intermediate
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     +                VALUE, IERR)
         IF (IERR .OR. NINTM + 1 .GT. KK) THEN
            WRITE (LOUT, *) ' ERROR READING DATA FOR KEYWORD ' // KEY
         ELSE
            NINTM = NINTM + 1
            KI(NINTM) = KSPEC
            XINTM(KSPEC) = VALUE(1)
         ENDIF
C
      ELSEIF (KEY .EQ. 'PROD') THEN
C        Product
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     +                VALUE, IERR)
         IF (IERR .OR. NPROD + 1 .GT. KK) THEN
            WRITE (LOUT, *) ' ERROR READING DATA FOR KEYWORD ' // KEY
         ELSE
            NPROD = NPROD + 1
            KP(NPROD) = KSPEC
            PROD(KSPEC) = VALUE(1)
         ENDIF
C
      ELSEIF (KEY .EQ. 'TEMP') THEN
C        Read specified temperature profile (X, T) pairs
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         IF (NTEMP + 1 .GT. JMAX) THEN
            WRITE (LOUT, *)
     +         ' ERROR... THE PROBLEM IS ONLY DIMENSIONED FOR ', JMAX,
     +         ' (X, T) PAIRS'
            IERR = .TRUE.
         ELSE
            NTEMP = NTEMP + 1
            XX(NTEMP) = VALUE(1)
            TT(NTEMP) = VALUE(2)
         ENDIF
C
      ELSEIF (KEY .EQ. 'USTG') THEN
C        If a restart, use given temperature profile, not the
C        one one the restart file
         LUSTGV = .TRUE.
C
      ELSEIF (KEY .EQ. 'TINF') THEN
C        JFG: On restart without USTG, set a new inflow temperature.
         linflow = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, tinflow, IERR)
C
      ELSEIF (KEY .EQ. 'TFIX') THEN
C        Temperature which is to be held fixed for a free flame
         NEC(6) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TFIXT, IERR)
C
C///////////////////////////////////////////////////////////////////////
C
C     TRANSPORT OPTIONS KEYWORDS
C
C///////////////////////////////////////////////////////////////////////
C
      ELSEIF (KEY .EQ. 'MULT') THEN
C        Multicomponent formulas used
         LMULTI = .TRUE.
C
      ELSEIF (KEY .EQ. 'MIX') THEN
C        Mixture-averaged formulas used
         LMULTI = .FALSE.
C
      ELSEIF (KEY .EQ. 'TDIF') THEN
C        Thermal diffusion included
         LTDIF = .TRUE.
C
      ELSEIF (KEY .EQ. 'VCOR') THEN
C        Use correction velocity formalism
         LVCOR = .TRUE.
C
      ELSEIF (KEY .EQ. 'TRCE') THEN
C        Use "trace" approximation, lump all transport errors into
C        the "last" species
         LVCOR = .FALSE.
C
C///////////////////////////////////////////////////////////////////////
C
C     SENSITIVITY KEYWORDS
C
C///////////////////////////////////////////////////////////////////////
C
      ELSEIF (KEY .EQ. 'ASEN') THEN
C        All reaction sensitivity
         LASEN = .TRUE.
C
      ELSEIF (KEY .EQ. 'HSEN') THEN
C        Sensitivity to heats of formation
         LHSEN = .TRUE.
C
C///////////////////////////////////////////////////////////////////////
C
C     PRINTING AND RESTARTING KEYWORDS
C
C///////////////////////////////////////////////////////////////////////
C
      elseif (key .eq. "ATOM") then
c        Write molar concentration of atoms.
         latom = .true.

      elseif (key .eq. "DFLX") then
c        Write diffusion fluxes.
         ldflx = .true.

      elseif (key .eq. "LWSN") then
c        Write Lewis numbers.
         lewis = .true.

      elseif (key .eq. "DIFF") then
c        Write diffusivities.
         ldiff = .true.

      elseif (key .eq. "HEAT") then
c        Write heat release.
         lheat = .true.

      elseif (key .eq. "PATH") then
c        Write reaction data to path.dat for drawing path diagrams.
         lpath = .true.

      elseif (key .eq. "RROP") then
c        Write reaction rates of progress to rrop.dat.
         lrrop = .true.

      ELSEIF (KEY .EQ. "FUEL") THEN
C        Name of fuel species for fuel consumption output
         lfuel = .true.
         first = 0
         last = 0
         do i = 1, len (line)
            if (line (i : i) .eq. " ") then
               if (first .ne. 0 .and. last .eq. 0) last = i
            else
               if (first .eq. 0) first = i
            end if
         end do
         fuel = line (first : last)
         do k = 1, kk
            if (fuel .eq. ksym(k)) ifuel = k
         end do
C
      ELSEIF (KEY .EQ. 'PRNT') THEN
C        Print level control
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IPRNT = INT (VALUE(1))
         IF (IPRNT .LT. 10) THEN
            IPRNT = MIN (2, IPRNT)
            CALL TWSETI (ERROR, LOUT, 'LEVELM', IPRNT)
            CALL TWSETI (ERROR, LOUT, 'LEVELD', IPRNT)
         ELSE
            IPRNT1 = IPRNT / 10
            IPRNT2 = IPRNT - 10 * IPRNT1
            CALL TWSETI (ERROR, LOUT, 'LEVELD', IPRNT2)
            IPRNT3 = MAX (IPRNT1, IPRNT2)
            CALL TWSETI (ERROR, LOUT, 'LEVELM', IPRNT3)
         ENDIF
C
      ELSEIF (KEY .EQ. 'SKIP') THEN
C        Restart skips
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         MFILE = INT(VALUE(1)) + 1
C
      ELSEIF (KEY .EQ. 'RSTR') THEN
C        Restart check
         LRSTRT = .TRUE.
C
      ELSEIF (KEY .EQ. 'CNTN') THEN
C        Continuiation flag
         LCNTUE = .TRUE.
C
      ELSEIF (KEY .EQ. 'JJRG') THEN
C        Number of mesh points in a regrid
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         JJREGD = INT(VALUE(1))
         LREGRD = .TRUE.
C
      ELSEIF (KEY .EQ. 'PCAD') THEN
C        Percentage of regrid points dedicated to adaption
         CALL CKXNUM (LINE, 1, LOUT, NVAL, PCTADP, IERR)
C
      ELSEIF (KEY .EQ. 'RGTC') THEN
C        Ratio of gradient regrid points to curvate points
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RATGTC, IERR)
C
      ELSEIF (KEY .EQ. 'END ') THEN
C        End of input for this job
         KERR = KERR.OR.IERR.OR.ERROR
         GO TO 0500
C
      ELSEIF (KEY .EQ. 'NTOT') THEN
C        Maximum points allowed for this job
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         NTOT = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'NADP') THEN
C        Maximum points that can be added per adaption
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETI (ERROR, LOUT, 'PADD', INT(VALUE(1)))
C
      ELSEIF (KEY .EQ. 'NOFT') THEN
C        Do not do the fixed temperature problem
         N1CALL = 2
C
      ELSEIF (KEY .EQ. 'SPOS') THEN
C        Convert negative species solutions
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SPOS, IERR)
C
      ELSEIF (KEY .EQ. 'GFAC') THEN
C        Scalar for gas-phase rate constants
         CALL CKXNUM (LINE, 1, LOUT, NVAL, GFAC, IERR)
c
      ELSEIF (KEY .EQ. 'DFTH') THEN
C        Scalar for unity Lewis number diffusion velocities for 
C        flame thickening
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DFTH, IERR)
c
      ELSEIF (KEY .EQ. 'CHTH') THEN
C        Scalar for reaction chemistry for flame thickening
         CALL CKXNUM (LINE, 1, LOUT, NVAL, CHTH, IERR)
c
      ELSEIF (KEY .EQ. 'LULN') THEN
C        Scalar for unity Lewis number
         LULNUM = .TRUE.
C
C///////////////////////////////////////////////////////////////////////
C
C     BOTTOM OF THE LOOP OVER INPUT LINES
C
C///////////////////////////////////////////////////////////////////////
C
      ELSE
C        To get here, an invalid keyword was read
         WRITE (LOUT, *) ' ERROR...ILLEGAL KEYWORD:',KEY
         KERR = .TRUE.
      ENDIF
C
C     BOTTOM OF THE LOOP
C
      KERR = KERR.OR.IERR.OR.ERROR
      GO TO 0200
0500  CONTINUE
C
C///////////////////////////////////////////////////////////////////////
C
C     CHECK AND PREPARE THE DATA SUPPLIED BY THE KEYWORKDS
C
C///////////////////////////////////////////////////////////////////////
C
C     Check the reactant and product sums
      SUMR = 0.0
      SUMP = 0.0
      DO 0600 K = 1, KK
         SUMR = SUMR + REAC(K)
         SUMP = SUMP + PROD(K)
0600  CONTINUE
C
C     Normalize reactant and product fractions
      DO 0700 K = 1, KK
         REAC(K) = REAC(K) / SUMR
         PROD(K) = PROD(K) / SUMP
0700  CONTINUE
C
      IF (ABS (SUMR - 1.0) .GT. 1.E-6)
     +   WRITE (LOUT, *) ' CAUTION...REACTANT FRACTIONS SUM TO ', SUMR
      IF ((.NOT.CNTNUD) .AND. ABS (SUMP - 1.0) .GT. 1.E-6)
     +   WRITE (LOUT, *) ' CAUTION...PRODUCT FRACTIONS SUM TO ',  SUMP
C
C     Check for necessary input
      IF (.NOT. NEC(8) )THEN
         WRITE (LOUT, *) ' ERROR..."BURN" OR "FREE" NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(1) ) THEN
         WRITE (LOUT, *)
     +      ' ERROR...MUST SPECIFY EITHER "MOLE" OR "MASS" '
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. LBURNR) LENRGY = .TRUE.
C
      IF ((.NOT. NEC(2)) .AND. (LBURNR)) THEN
         WRITE (LOUT, *)
     +   ' ERROR..."ENRG" OR "TGIV" MUST BE PROVIDED FOR A BURNER FLAME'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(3) ) THEN
         WRITE (LOUT, *) ' ERROR...PRESSURE NOT GIVEN'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(4) .AND. .NOT. LRSTRT) THEN
         WRITE (LOUT, *) ' ERROR...MASS FLOW RATE NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(5)) THEN
         WRITE (LOUT, *) ' ERROR..."XEND" NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF ((.NOT. NEC(6)) .AND. (.NOT. LBURNR)) THEN
         WRITE (LOUT, *) ' ERROR..."TFIX" NOT GIVEN FOR A FREE FLAME'
         KERR = .TRUE.
      ENDIF
C
C     MAKE SURE THE (X, T) PAIRS ARE IN ORDER
C
      DO 0800 N = 2, NTEMP
         IF (XX(N - 1) .GE. XX(N)) THEN
            WRITE (LOUT, *)
     +         ' ERROR...SPECIFIED TEMPERATURES ARE OUT OF ORDER'
            KERR = .TRUE.
         ENDIF
0800  CONTINUE
C
C     Make sure the initial gridpoints are in order.
      IF ((.NOT.CNTNUD) .AND. (.NOT.LUMESH)) THEN
         JJ = NP
         DO 0900 N = 2, JJ
            IF (X(N - 1) .GE. X(N)) THEN
               WRITE (LOUT, *) ' ERROR...INITIAL GRID IS OUT OF ORDER'
               KERR = .TRUE.
            ENDIF
0900     CONTINUE
      ENDIF
C
C     Make sure the given temperatures span the XEND-XSTR domain
      IF ((.NOT.LRSTRT .AND. .NOT.CNTNUD) .OR. LUSTGV) THEN
         IF (XX(1) .GT. XSTR .OR. XX(NTEMP) .LT. XEND) THEN
            WRITE (LOUT, *)
     +      ' ERROR...GIVEN TEMPERATURE PROFILE DOES NOT SPAN XEND-XSTR'
            KERR = .TRUE.
         ENDIF
      ENDIF
C
C     Set optional input if needed
      IF (.NOT. NOPT(1)) WMIX = (XEND - XSTR) * 0.50
      IF (.NOT. NOPT(2)) XCEN = (XEND - XSTR) * 0.35
C
C     Set adaption.
      CALL TWSETL (ERROR, LOUT, 'ADAPT', .NOT. LENRGY)
      KERR = KERR .OR. ERROR
C
C///////////////////////////////////////////////////////////////////////
C
C     EXIT.
C
C///////////////////////////////////////////////////////////////////////
C
99999 CONTINUE
C
C     end of SUBROUTINE RDKEY
      RETURN
      END
C
      SUBROUTINE RESTRT (JMAX, LOUT, LMOLE, LBURNR, LUSTGV, REAC,
     1                   XGIVEN, TGIVEN, ICKWRK, RCKWRK, EPS, X, S,
     2                   KERR, linflow, tinflow)
C
C  START PROLOGUE
C
C  JMAX     - integer scalar, maximum number of mesh points
C  LOUT     - integer scalar, formatted output file unit number
C  LMOLE    - logical, T, input and output is mole fractions,
C                      F, input/poutput is mass fractions
C  LBURNR   - logical, T, burner stabilized flame problem,
C                      F, freely propagating adiabatic flame
C  LUSTGV   - logical, T, use given temperature profile,
C                      F, use the restart file profile
C  REAC(*)  - real array, species input mole/mass fractions.
C  XGIVEN(*)- real array, distances from the burner for specified
C             temperatures (cm)
C  TGIVEN(*)- real array, specified temperature profile (K)
C  ICKWRK(*)- integer array, CHEMKIN workspace
C  RCKWRK(*)- real array, CHEMKIN workspace
C  EPS(*)   - real array, inlet mass flux fractions
C  X(*)     - real array, starting mesh locations (cm)
C  S(*,*)   - real matrix, starting solution.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
C     Integer arrays
      DIMENSION ICKWRK(LENICK)
C     Real arrays
      DIMENSION REAC(KK), RCKWRK(LENRCK), EPS(KK), XGIVEN(JMAX),
     1          TGIVEN(JMAX), X(JMAX), S(NATJ,JMAX)
      LOGICAL LMOLE, LBURNR, LUSTGV, KERR, linflow
C
      KERR = .FALSE.
C     Initialize mass flux fractions
      DO 100 K = 1, KK
         EPS(K) = 0.0
100   CONTINUE
C
      IF (XSTR .LT. X(1)) THEN
C        Shift the solution to the right if a new XSTR is < X(1)
         JJ = JJ + 1
         IF (JJ .GT. JMAX) THEN
            WRITE (LOUT, *) ' ERROR...NEW XSTR NEEDS TOO MANY POINTS'
            KERR = .TRUE.
            RETURN
         ENDIF
         DO 200 I = 2, JJ
            J = JJ + 2 - I
            JM1 = J - 1
            X(J) = X(JM1)
            CALL CKCOPY (NATJ, S(1,JM1), S(1,J))
200      CONTINUE
         X(1) = XSTR
         CALL CKCOPY (NATJ, S(1,2), S(1,1))
      ENDIF
C
      IF (XEND .GT. (X(JJ)+1.E-4) ) THEN
C        If a new XEND > X(JJ) then add a point at JJ+1
C        if a new XEND <=X(JJ) then reduce points and X(JJ)=XEND
C
         JJ = JJ + 1
         IF (JJ .GT. JMAX) THEN
            WRITE (LOUT, *) ' ERROR...NEW XEND NEEDS TOO MANY POINTS'
            KERR = .TRUE.
            RETURN
         ENDIF
         X(JJ) = XEND
         JJM1 = JJ - 1
         CALL CKCOPY (NATJ, S(1,JJM1), S(1,JJ))
      ENDIF
C
C     Set the max flux fraction boundary conditions
      IF (LMOLE) THEN
         CALL CKXTY (REAC, ICKWRK, RCKWRK, EPS)
      ELSE
         CALL CKCOPY (NREAC, EPS, REAC)
      ENDIF
C
C     If not a burner, overwrite FLRT from the restart file,
C     otherwise, set FLRT from the new keyword.
      IF (.NOT. LBURNR) FLRT = S(NM,1)
C
C     Set the new flow rate profile
      DO 1200 J = 1, JJ
         S(NM,J) = FLRT
1200  CONTINUE
C
C     Set XGIVEN and TGIVEN to the old solution
      IF (.NOT. LUSTGV) THEN
         NTEMP = JJ
         DO 1250 N = 1, NTEMP
            XGIVEN(N) = X(N)
            if (linflow) then
               if (n .eq. 1) then
                  TGIVEN(N) = tinflow
               else
                  TGIVEN(N) = max (tinflow, S(NT,N))
               end if
            else
               TGIVEN(N) = S(NT,N)
            end if
1250     CONTINUE
      ENDIF
C
      IF (LBURNR) RETURN
C
C     For a free flame set the mesh point for the fixed temperature;
C     since on the restart file, TFIXT may not exactly equal S(NT,J),
C     we check for a point with 2K of TFIXT and change TFIXT if
C     necessary.
C
      DO 1300 J = 1, JJ
         IF ( ABS(S(NT,J)-TFIXT) .LT. 2.0) THEN
            JFIXT = J
            TFIXT = S(NT,J)
         ENDIF
1300  CONTINUE
C
C     end of SUBROUTINE RESTRT
      RETURN
      END
C
      SUBROUTINE START (LOUT, LMOLE, LUMESH, LBURNR, REAC, XINTM,
     1                  PROD, KR, KI, KP, XGIVEN, TGIVEN, ICKWRK,
     3                  RCKWRK, Y, SI, EPS, X, S, JMAX, KERR)
C
C  START PROLOGUE
C
C  LOUT     - integer scalar, formatted output file unit number
C  LMOLE    - logical, T, input/output mole fractions,
C                      F, input/output mass fractions
C  LUMESH   - logical, T, use uniform starting mesh
C  LBURNR   - logical, T, burner stabilized flame probem,
C                      F, freely propagating adiabatic flame
C  REAC(*)  - real array, reactant species fractions
C  XINTM(*) - real array, intermediate species fractions
C  PROD(*)  - real array, product species fractions
C  KR(*)    - integer array, reactant species indices
C  KI(*)    - integer array, intermediate species indices
C  KP(*)    - integer array, product species indices
C  XGIVEN(*)- real array, distances from the burner of given
C             temperatures
C  TGIVEN(*)- real array, given temperatures
C  ICKWRK(*)- integer array, CHEMKIN workspace
C  RCKWRK(*)- real array, CHEMKIN workspace
C  Y(*)     - real array, mass fractions
C  SI(*)    - real array, sum of intermediates at gridpoints
C  EPS(*)   - real array, inlet mass flux fractions.
C  X(*)     - real array, grid point distances from burner
C  S(*,*)   - real matrix, starting solution estimates
C  JMAX     - integer scalar, maximum number of gridpoints
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
C
      LOGICAL LBURNR, LMOLE, LUMESH, KERR
C     Integer arrays
      DIMENSION ICKWRK(LENICK), KI(NINTM), KP(NPROD), KR(NREAC)
C     Real arrays
      DIMENSION EPS(KK), PROD(KK), RCKWRK(LENRCK), REAC(KK),
     2          S(NATJ,JMAX), SI(JMAX), TGIVEN(JMAX), X(JMAX),
     3          XGIVEN(JMAX), XINTM(KK), Y(KK)
C
C     Initialize mass flux fractions and mass fractions.
      ZERO = 0.0
      DO 100 K = 1, KK
         EPS(K) = ZERO
         DO 95 J = 1, JJ
            S(NYS+K, J) = ZERO
 95      CONTINUE
100   CONTINUE
C
      IF (LUMESH) THEN
C        Set uniform X mesh coordinates
C
         DX = (XEND-XSTR) / (JJ-1)
         DO 200 J = 1, JJ
            X(J) = XSTR + DX*(J-1)
200      CONTINUE
      ENDIF
C
      IF (.NOT. LBURNR) THEN
C        For free flames, add the fixed temperature point to the mesh
C
         DO 300 N = 1, NTEMP
            NFIXT = N
            IF (TGIVEN(N) .GE. TFIXT) GO TO 350
300      CONTINUE
         WRITE (LOUT,*) ' ERROR...NO USABLE MESH LOCATION FOR TFIXT'
         KERR = .TRUE.
         RETURN
350      CONTINUE
C
         IF (TGIVEN(NFIXT) .EQ. TFIXT) THEN
            XFIXT = XGIVEN(NFIXT)
         ELSE
            XFIXT = XGIVEN(NFIXT-1) + (TFIXT-TGIVEN(NFIXT-1)) *
     1              (XGIVEN(NFIXT) - XGIVEN(NFIXT-1)) /
     2              (TGIVEN(NFIXT) - TGIVEN(NFIXT-1))
         ENDIF
C
         DO 400 J = 1, JJ
            IF (XFIXT .EQ. X(J)) THEN
               JFIXT = J
               GO TO 700
            ENDIF
400      CONTINUE
C
         DO 500 J = 1, JJ-1
            JP1 = J + 1
            IF (XFIXT.GT.X(J) .AND. XFIXT.LT.X(JP1)) JFIXT = JP1
500      CONTINUE
C
         IF (JMAX .LE. JJ) THEN
            WRITE (LOUT,*) ' ERROR...THE GRID CANOT BE ENLARGED'
            KERR = .TRUE.
            RETURN
         ENDIF
         JJ = JJ + 1
         DO 600 J = JJ, JFIXT+1, -1
            JM1 = J - 1
            X(J) = X(JM1)
            CALL CKCOPY (NATJ, S(1,JM1), S(1,J))
  600    CONTINUE
C
         X(JFIXT) = XFIXT
         S(NT, JFIXT) = TFIXT
700      CONTINUE
C
      ENDIF
C
C     Set intermediate gaussians
      W = -LOG(0.15) / (WMIX/2.)**2
      DO 800 N = 1, NINTM
         K = KI(N)
         DO 800 J = 1, JJ
            S(NYS+K, J) = XINTM(K) * EXP(-W * (X(J)-XCEN)**2)
800   CONTINUE
C
      DO 1000 J = 1, JJ
C
C        Sum the intermediates
         SI(J) = ZERO
         DO 900 N = 1, NINTM
            SI(J) = SI(J) + S(NYS+KI(N), J)
900      CONTINUE
C
C        Set starting species profiles
         CALL LINWMX (WMIX, XCEN, X(J), XRE, XPD)
         FAC = 1.0 - SI(J)
         DO 1100 N = 1, NREAC
            S(NYS+KR(N), J) = (XPD*PROD(KR(N)) + XRE*REAC(KR(N)))*FAC
1100     CONTINUE
C
         DO 1050 N = 1, NPROD
            K = KP(N)
            DO 1025 L = 1, NREAC
               IF (K .EQ. KR(L)) GO TO 1050
1025        CONTINUE
            S(NYS+K, J) = (XPD*PROD(K) + XRE*REAC(K))*FAC
1050     CONTINUE
1000  CONTINUE
C
C     Set the mass flux fraction boundary conditions
      IF (LMOLE) THEN
         CALL CKXTY (REAC, ICKWRK, RCKWRK, EPS)
C
C        Convert starting extimates to mass fraction
         DO 1600 J = 1, JJ
            CALL CKXTY (S(NY, J), ICKWRK, RCKWRK, Y)
            CALL CKCOPY (KK, Y, S(NYS+1,J))
1600     CONTINUE
      ELSE
         CALL CKCOPY (NREAC, REAC, EPS)
      ENDIF
C
C     Set the temperature and flow rate profiles
      DO 1700 J = 1, JJ
         S(NT, J) = CKBSEC(NTEMP, X(J), XGIVEN, TGIVEN)
         S(NM,J) = FLRT
1700  CONTINUE
      IF (.NOT. LBURNR) S(NT,JFIXT)=TFIXT
C
C     end of SUBROUTINE START
      RETURN
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION VTRACE (S, KK)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION VTRACE (S, KK)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C  START PROLOGUE
C  In the case of LVCOR=.FALSE., use 'trace' approximation,
C  lumping all transport errors into final species
C
C  S(*) - real array, species at solution gridpoint
C  KK   - integer scalar, species count
C
      DIMENSION S(KK)
      SUMYK = 0.0
      DO 10 K = 1, KK
         SUMYK = SUMYK + S(K)
  10  CONTINUE
      VTRACE = 1.0 - SUMYK
C
C     end of FUNCTION VTRACE
      RETURN
      END
C
      SUBROUTINE REGRID (SNEW, SOLD, NVAR, JJNEW, JJOLD, XNEW, XOLD,
     1                   WORK, P, R, N )
C
C  START PROLOGUE
C
C  SNEW(*,*) - real matrix
C  SOLD(*,*) - real matrix
C  NVAR      - integer scalar
C  JJNEW     - integer scalar
C  JJOLD     - integer scalar
C  XNEW(*)   - real array
C  XOLD(*)   - real array
C  WORK(*)   - real array
C  P         - real scalar
C  R         - real scalar
C  N         - integer scalar
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION SNEW(NVAR,JJNEW), SOLD(NVAR,JJOLD), XNEW(JJNEW),
     1          XOLD(JJOLD), WORK(*)
C
C     interpolate the solution onto an equidistributing mesh
C
C       compute weight coefficients
C
      R0 = 1. - P
      R1 = P *R /(R+1.)
      R2 = P - R1
      TV1 = 0.
      DO 10 I = 2, JJOLD
        TV1 = TV1 + ABS( SOLD(N,I) - SOLD(N,I-1) )
10    CONTINUE
      TV2 = 0.
      DO 20 I = 2, JJOLD-1
        TV2 = TV2 + ABS( (SOLD(N,I+1)-SOLD(N,I))/(XOLD(I+1)-XOLD(I))
     1                 - (SOLD(N,I)-SOLD(N,I-1))/(XOLD(I)-XOLD(I-1)) )
20    CONTINUE
      XLEN = XOLD(JJOLD) - XOLD(1)
      B1 = R1 * XLEN /((1.-P) * TV1)
      B2 = R2 * XLEN /((1.-P) * TV2)
C
C     compute partial sums of weight function
C
      WORK(1) = 0.
      DO 50 I = 2, JJOLD
         DX = XOLD(I) -XOLD(I-1)
         WORK(I) = DX + B1*ABS( SOLD(N,I) - SOLD(N,I-1) ) + WORK(I-1)
     1           + B2*ABS( (SOLD(N,I+1)-SOLD(N,I))/(XOLD(I+1)-XOLD(I))
     2                   - (SOLD(N,I)-SOLD(N,I-1))/(XOLD(I)-XOLD(I-1)) )
50    CONTINUE
      DO 65 I = 2, JJOLD
         WORK(I) = WORK(I)/WORK(JJOLD)
65    CONTINUE
C
C      interpolate onto uniform eta grid to find new x
C
      XNEW(1) = XOLD(1)
      XNEW(JJNEW) = XOLD(JJOLD)
      ISTART  = 2
      DETA = 1./FLOAT(JJNEW-1)
      DO 80 J = 2, JJNEW-1
         ETAJ = (J-1)*DETA
         DO 70 I = ISTART, JJOLD
            IF (ETAJ .LE. WORK(I)) THEN
               IM1 = I - 1
               DEL = (ETAJ-WORK(IM1))/(WORK(I)-WORK(IM1))
               XNEW(J) = XOLD(IM1)+(XOLD(I)-XOLD(IM1))*DEL
               GO TO 80
            ELSE
               ISTART = I
            ENDIF
70       CONTINUE
         WRITE (6, *) ' *** VALUE OF ETA NOT FOUND ***'
80    CONTINUE
C
C        interpolate solution...
C
      do j=2,jjnew-1
         xnew(j) = xnew(1) + (j-1)*(xnew(jjnew)-xnew(1))/(jjnew-1)
      end do
      CALL INTPL8 (SOLD, SNEW, XOLD, XNEW, NVAR, JJOLD, JJNEW)
C
C     end of SUBROUTINE REGRID
      RETURN
      END
C
C------------------------------------------------------------------
      SUBROUTINE INTPL8 (F1, F2, X1, X2, MVAR, N1, N2)
C
C  START PROLOGUE
C
C  F1   -
C  F2   -
C  X1   -
C  X2   -
C  MVAR -
C  N1   -
C  N2   -
C
C  END PROLOGUE
C
C------------------------------------------------------------------
C
C       interpolate to get f2(x2) from f1(x1)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION F1(MVAR,*), F2(MVAR,*), X1(*), X2(*)
      DO 100 J = 2, N2-1
        XVAL = X2(J)
        DO 50 I = 2, N1
          IF (XVAL .LE. X1(I)) THEN
            DEL = (XVAL-X1(I-1)) / (X1(I)-X1(I-1))
            DO 40 L = 1, MVAR
              F2(L,J) = F1(L,I-1) + (F1(L,I) - F1(L,I-1))*DEL
40          CONTINUE
            GO TO 80
          ENDIF
50      CONTINUE
        WRITE(7,*) ' *** STOP...INTERPOLATION ERROR.'
80      CONTINUE
100   CONTINUE
C
C     endpoints..
C
      DO 110 L = 1, MVAR
        F2(L,1)  = F1(L,1)
        F2(L,N2) = F1(L,N1)
110   CONTINUE
C
C     end of SUBROUTINE INTPL8
      RETURN
      END
C
      SUBROUTINE PRREGD (SN, S, REG, X, COND, PCTADP, RATGTC,
     1                   XGIVEN, TGIVEN, LOUT, LBURNR, SCRK)
C
C  START PROLOGUE
C
C  SN(*,*)    - real matrix,
C  S(*,*)     - real matrix,
C  REG(*)     - real array,
C  X(*)       - real array, mesh point locations (cm)
C  COND(*)    - real array,
C  PCTADP     - real scalar,
C  RATGTC     - real scalar,
C  XGIVEN(*)  - real array, distances from the burner for TGIVEN
C  TGIVEN(*)  - real array, temperatures at given distances (K)
C  LOUT
C  LBURNR
C  SCRK(*,*)
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
      COMMON / PRICON / KK, II, NATJ, JJ, IASIZE, LENICK, LENRCK,
     1                  LENCCK, LENIMC, LENRMC,
     2                  LENITW, LENRTW, NT, NM, NYS, NY, NTR, NTEMP,
     3                  NREAC, NPROD, NINTM, JJREGD, JFIXT, NUMDT2,
     +                  ifuel
      COMMON / PRRCON / RELAT, ABSOL, SPOS, GFAC, PATM, FLRT, SFLR,
     1                  XSTR, XCEN, XEND, WMIX, WNDFAC, P, DT, DT2,
     2                  TFIXT, TRFAC
C
      LOGICAL LBURNR, LFIXT
      DIMENSION SN(NATJ,JJ), S(NATJ,JJ), REG(JJ), X(JJ), COND(JJ),
     1          XGIVEN(NTEMP), TGIVEN(NTEMP), SCRK(KK)
C
      WRITE (LOUT, *)
     1'REGRIDDING FROM ',JJ,' TO ', JJREGD, ' POINTS'
C
      LFIXT = .TRUE.
C
      IF (.NOT. LBURNR) THEN
	 JFIXT = 0
	 DO 10 J = 1, JJ
	    IF (S(1,J) .EQ. TFIXT) JFIXT = J
   10    CONTINUE
	 IF (JFIXT .NE. 0) THEN
	    XFIXT = X(JFIXT)
            CALL CKCOPY (KK, S(NYS+1,JFIXT), SCRK)
	    FLRT = S(KK+2, JFIXT)
         ELSE
	    WRITE (LOUT, *)
     1      ' CANNOT FIND FIXED TEMPERATURE IN REGRID...REGRID FAILS '
	    RETURN
         ENDIF
      ENDIF
C
  100 CONTINUE
C
      CALL REGRID (SN, S, NATJ, JJREGD, JJ, REG, X, COND, PCTADP,
     1             RATGTC, NT)
C
      IF (.NOT. LBURNR) THEN
         LFIXT = .FALSE.
         DO 150 J = 1, JJREGD
            IF (SN(1,J) .EQ. TFIXT) LFIXT = .TRUE.
  150    CONTINUE
      ENDIF
C
      IF (LFIXT) THEN
C
         CALL CKCOPY (JJREGD, REG, X)
         DO 380 J = 1, JJREGD
            CALL CKCOPY (NATJ, SN(1,J), S(1,J))
  380    CONTINUE
	 JJ = JJREGD
C
      ELSE
C
	 JFIXT = 0
	 DO 390 J = 1, JJREGD-1
	    IF (REG(J).LE.XFIXT .AND. REG(J+1).GE.XFIXT) JFIXT=J
  390    CONTINUE
C
	 IF (JFIXT .EQ. 0) THEN
	    WRITE (LOUT, *)
     1      ' CANNOT FIND FIXED TEMPERATURE ON NEW GRID...REGRID FAILS '
	    RETURN
C
         ELSE
            CALL CKCOPY (JFIXT, REG, X)
	    DO 400 J = 1, JFIXT
               CALL CKCOPY (NATJ, SN(1,J), S(1,J))
  400       CONTINUE
C
	    JFIXT = JFIXT + 1
	    X(JFIXT) = XFIXT
            S(1, JFIXT) = TFIXT
	    S(KK+2,JFIXT) = FLRT
            CALL CKCOPY (KK, SCRK, S(NYS+1,JFIXT))
C
	    DO 420 J = JFIXT, JJREGD
               JP1 = J + 1
	       X(JP1) = REG(J)
               CALL CKCOPY (NATJ, SN(1,J), S(1,JP1))
  420       CONTINUE
C
	    JJ = JJREGD + 1
         ENDIF
C
      ENDIF
C
      IF (SPOS .GT. 0.0) THEN
         DO 376 J = 1, JJ
            CALL PRPOS (KK, S(NYS+1,J), SPOS)
  376    CONTINUE
      ENDIF
C
      NTEMP = JJ
      DO 450 N = 1, NTEMP
	 XGIVEN(N) = X(N)
	 TGIVEN(N) = S(NT,N)
  450 CONTINUE
C
C     end of SUBROUTINE PRREGD
      RETURN
      END
C
      SUBROUTINE PRPOS (IDIM, S, SPOS)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER(I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER(I-N)
C*****END precision > single
C
      DIMENSION S(IDIM)
      DO 10 I = 1, IDIM
         S(I) = MAX (S(I), SPOS)
   10 CONTINUE
C
C     end of SUBROUTINE PRPOS
      RETURN
      END
