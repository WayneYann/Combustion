C     CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      PROGRAM SHOCK
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for SHOCK
C  The file is used to make 'shock.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     LIN     user Keyword input
C     LOUT    solution and diagnostic printing (output)
C     LINCK   gas-phase Chemkin Linking File (input)
C     LSAVE   binary solution Save File (output)
C
C     dimensions:
C
C     LIPAR   maximum integer workspace available for Shock
C     LRPAR   maximum real workspace available for Shock
C     LCPAR   maximum character  workspace available for Shock
C
C///////////////////////////////////////////////////////////////////////
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (LRPAR=4000, LIPAR=4000, LCPAR=100, LIN=5, LOUT=6,
     1           LINCK=25, LSAVE=55, ZERO=0.0)
C
      DIMENSION RPAR(LRPAR), IPAR(LIPAR)
      CHARACTER*16 CPAR(LCPAR)
C
C     use appropriate machine-dependent cktime.f
C
      EXTERNAL CKTIME
      TSTART = CKTIME (ZERO)
C
C///////////////////////////////////////////////////////////////////////
C   Note:  The ascii-formatted linking file option is the recommended 
C          standard, as shown below in the "linking files > ascii"
C          section.  Any change to these open statements requires
C          systematic changes to all of the Chemkin driver routines.
C///////////////////////////////////////////////////////////////////////
C
C*****linking files > binary
C      OPEN (LINCK,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE='./chem.bin')
C*****end linking files > binary
C*****linking files > ascii
      OPEN (LINCK,STATUS='UNKNOWN',FORM='FORMATTED',FILE='./chem.asc') 
C*****end linking files > ascii
C
      OPEN (LSAVE,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./save.bin')
C
      CALL SHRUN (LINCK,LIN,LOUT,LSAVE,LIPAR,LRPAR,LCPAR,IPAR,RPAR,CPAR)
C
      TEND = CKTIME (TSTART)
      IF (TEND .GT. 60.) THEN
         WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (min): ',TEND/60.
      ELSE
          WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (sec): ',TEND
      ENDIF
C
      CLOSE (LINCK)
      CLOSE (LSAVE)
      STOP
      END
