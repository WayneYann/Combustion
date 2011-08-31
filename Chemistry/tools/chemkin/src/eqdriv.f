C     CVS: $Revision: 1.1.1.1 $ deposited $Date: 2006/05/26 19:09:33 $

      PROGRAM EQDRIV

C///////////////////////////////////////////////////////////////////////
C
C     This is the driver routine for EQUIL, the equilibrium code.
C     The file is used to make 'equil.exe'
C
C     The parameters and unit numbers that may be changed
C     by the user are described below:
C
C     unit numbers:
C
C     LIN     user Keyword input
C     LOUT    formatted solution and diagnostic printing (output)
C     LINCK   gas-phase Chemkin Linking File
C     LINSK   Surface Chemkin Linking File
C     LSAVE   binary solution file (output)
C
C     dimensions:
C
C     LIWORK  maximum integer workspace available for Equil
C     LRWORK  maximum real  workspace available for Equil
C     LCWORK  maximum character workspace available for Equil
C
C     Other parameters controlling equil calculations:
C
C     LSURF   an integer flag indicating whether or not the equilibrium
C             calculations will involve surface (bulk) species.  This
C             driver is currently set up to look for a Surface Chemkin
C             Linking File, 'surf.asc'.  If the file is not found, a
C             gas-phase only equilibrium calculation is assumed.
C             The LSURF parameter is set accordingly:
C
C             = 0   gas-phase only
C             = 1   surface (bulk) species may be included
C
C///////////////////////////////////////////////////////////////////////

C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single

      CHARACTER*16 C, FORM, SUFFIX
      EXTERNAL CKTIME
      LOGICAL FOUND

      PARAMETER
     +  (LIN = 5,
     +   LINCK = 25,
     +   LINSK = 26,
     +   LOUT = 6,
     +   LSAVE = 7,
     +   LIWORK = 5000,
     +   LRWORK = 6000,
     +   LCWORK = 500)

      DIMENSION C(LCWORK), I(LIWORK), R(LRWORK)

      ZERO = 0
      TSTART = CKTIME (ZERO)
C
C///////////////////////////////////////////////////////////////////////
C   Note:  The ascii-formatted linking file option is the recommended 
C          standard, as shown below in the "linking files > ascii"
C          section.  Any change to these open statements requires
C          systematic changes to all of the Chemkin driver routines.
C///////////////////////////////////////////////////////////////////////
C
C*****linking files > ascii
      FORM = 'FORMATTED'
      SUFFIX = '.asc'
C*****end linking files > ascii
C*****linking files > binary
C      FORM = 'UNFORMATTED'
C      SUFFIX = '.bin'
C*****end linking files > binary

      OPEN
     +  (FILE = 'chem' // SUFFIX,
     +   FORM = FORM,
     +   STATUS = 'UNKNOWN',
     +   UNIT = LINCK)

      INQUIRE (FILE = 'surf' // SUFFIX, EXIST = FOUND)

      IF (FOUND) THEN
         OPEN
     +     (FILE = 'surf' // SUFFIX,
     +      FORM = FORM,
     +      STATUS = 'OLD',
     +      UNIT = LINSK)
         LSURF = 1
      ELSE
         LSURF = 0
      ENDIF

      OPEN
     +  (FILE = 'equil.bin',
     +   FORM = 'UNFORMATTED',
     +   STATUS = 'UNKNOWN',
     +   UNIT = LSAVE)

C     SUBROUTINE EQINTP (LINCK, LINSK, LIN, LOUT, LSAVE, LSURF,
C    +                   LIWORK, I, LRWORK, R, LCWORK, C)

      CALL EQINTP (LINCK, LINSK, LIN, LOUT, LSAVE, LSURF,
     +             LIWORK, I, LRWORK, R, LCWORK, C)

      TEND = CKTIME (TSTART)
      IF (TEND .GT. 60) THEN
         WRITE (LOUT, '(/A, 1PE15.2)')
     +      ' Total CPU time (min): ', TEND / 60
      ELSE
         WRITE (LOUT, '(/A, 1PE15.2)')
     +      ' Total CPU time (sec): ', TEND
      ENDIF

      END
