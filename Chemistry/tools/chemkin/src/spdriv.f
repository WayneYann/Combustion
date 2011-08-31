C     CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      PROGRAM DRIVER
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for SPIN, the rotating-disk, 
C  stagnation-flow code.
C  The file is used to make 'spin.exe'
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
C     LINMC   Transport Linking File (input)
C     LROUT   binary solution Save File (output)
C     LRIN    binary Restart File (input)
C     LRCVR   binary solution Recover File (output)
C     LFLUX   formatted printing of resulting species fluxes at the
C             disk surface (output)
C
C     dimensions:
C
C     LENIWK  maximum integer workspace available for Spin
C     LENRWK  maximum real  workspace available for Spin
C     LENCWK  maximum character workspace available for Spin
C     LENLWK  maximum logical workspace available for Spin
C     NMAX    total number of grid points allowed for grid refinement.
C
C///////////////////////////////////////////////////////////////////////
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (LENRWK=900000, LENIWK=8000, LENLWK=300, LENCWK=200,
     1           LIN=5, LOUT=6, LFLUX=10, LINCK=25, LINMC=35,
     2           LINSK=26, LRIN=14, LROUT=15, LRCVR=16, NMAX=200,
     3           ZERO=0.0)
C
      DIMENSION RWORK(LENRWK), IWORK(LENIWK)
      LOGICAL LWORK(LENLWK)
      CHARACTER*16 CWORK(LENCWK)
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
C      OPEN(LINCK,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./chem.bin')
C      OPEN(LINMC,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./tran.bin')
C      OPEN(LINSK,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./surf.bin')
C*****end linking files > binary
C*****linking files > ascii
      OPEN(LINCK,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./chem.asc')
      OPEN(LINMC,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./tran.asc')
      OPEN(LINSK,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./surf.asc')
C*****end linking files > ascii
C
      OPEN(LRIN,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./rest.bin')
      OPEN(LROUT,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./save.bin')
      OPEN(LRCVR,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./recov.bin')
      OPEN(LFLUX,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./flux.out')
C
      CALL SPIN (NMAX, LIN, LOUT, LINCK, LINSK, LINMC, LRIN, LROUT,
     1           LRCVR, LFLUX, LENLWK, LWORK, LENIWK, IWORK, LENRWK,
     2           RWORK, LENCWK, CWORK)
C
      TEND = CKTIME (TSTART)
      IF (TEND .GT. 60.) THEN
         WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (min): ',TEND/60.
      ELSE
          WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (sec): ',TEND
      ENDIF
C
      CLOSE(LINCK)
      CLOSE(LINMC)
      CLOSE(LINSK)
      CLOSE(LFLUX)
      STOP
      END
