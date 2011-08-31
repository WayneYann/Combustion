C     CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      PROGRAM PSRMN
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for AURORA, the perfectly stirred reactor
C  code for use in modeling thermal and plasma reactors, with and 
C  without surface chemistry.
C  The file is used to make 'aurora.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     LIN     user Keyword input
C     LOUT    formatted solution and diagnostic printing (output)
C     LINCK   gas-phase Chemkin Linking File (input)
C     LINSK   Surface Chemkin Linking File (input)
C     LSAVE   binary solution Save File (output)
C     LREAD   binary Restart File (input)
C     LRECOV  binary solution Recover File (output)
C     LOUTSH  formatted printing of plasma sheath variables
C             (only used in plasma cases when the RF sheath option
C              is selected) (output)
C
C     dimensions:
C
C     LNIWRK  maximum integer workspace available for Aurora
C     LNRWRK  maximum real  workspace available for Aurora
C     LNCWRK  maximum character workspace available for Aurora
C     LNLWRK  maximum logical workspace available for Aurora
C     MAXPSR  maximum number of perfectly stirred reactors in series
C     MAXMAT  maximum number of different materials where surface
C             chemistry mechanisms are specified.
C
C///////////////////////////////////////////////////////////////////////
C
C*****precision > quad
C      IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (LNIWRK=20000, LNRWRK=45000, LNCWRK=2000, LNLWRK=1000,
     1           LIN=5, LOUT=6, LINCK=25, LINSK=26, LREAD=14, 
     2           LSAVE=15, LRECOV=16, LOUTSH=31, MAXPSR=3, MAXMAT=3, 
     3           ZERO=0.0)
C
      COMMON IWORK(LNIWRK), RWORK(LNRWRK)
      CHARACTER*16 CWORK(LNCWRK)
      LOGICAL LWORK(LNLWRK)
      EXTERNAL CKTIME
      TSTART = CKTIME(ZERO)
C
C///////////////////////////////////////////////////////////////////////
C   Note:  The ascii-formatted linking file option is the recommended 
C          standard, as shown below in the "linking files > ascii"
C          section.  Any change to these open statements requires
C          systematic changes to all of the Chemkin driver routines.
C///////////////////////////////////////////////////////////////////////
C
C*****linking files > binary
C      OPEN(LINCK,STATUS='OLD',FORM='UNFORMATTED',FILE='./chem.bin')
C      OPEN(LINSK,STATUS='OLD',FORM='UNFORMATTED',FILE='./surf.bin')
C*****end linking files > binary
C*****linking files > ascii
      OPEN(LINCK,STATUS='OLD',FORM='FORMATTED',FILE='./chem.asc')
      OPEN(LINSK,STATUS='OLD',FORM='FORMATTED',FILE='./surf.asc')
C*****end linking files > ascii
C
      OPEN(LREAD,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE='./rest.bin')
      OPEN(LSAVE,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE='./save.bin')
      OPEN(LRECOV,STATUS='UNKNOWN',FORM='UNFORMATTED',
     1     FILE='./recov.bin')
      OPEN(LOUTSH,STATUS='UNKNOWN',FORM='FORMATTED',FILE='./last.s')
C
      CALL PSRUN (LIN, LOUT, LINCK, LINSK, LREAD, LSAVE, LRECOV, LOUTSH,
     1            LNLWRK, LWORK, LNIWRK, IWORK, LNRWRK, RWORK, LNCWRK, 
     2            CWORK, MAXPSR, MAXMAT)
C
      TEND = CKTIME (TSTART)
      IF (TEND .GT. 60.) THEN
         WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (min): ',TEND/60.
      ELSE
          WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (sec): ',TEND
      ENDIF
C
      CLOSE (LINCK)
      CLOSE (LINSK)
      CLOSE (LREAD)
      CLOSE (LSAVE)
      CLOSE (LRECOV)
      STOP
      END
