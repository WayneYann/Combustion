C     CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      PROGRAM SDRIV
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for SKSAMPLE, the example code for
C  Surface Chemkin.
C  The file is used to make 'sksample.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     LIN     user Keyword input
C     LOUT    formatted solution and diagnostic printing (output)
C     LINCK   gas-phase Chemkin Linking File 
C     LINSK   Surface Chemkin Linking File
C
C     dimensions:
C
C     LIWORK  maximum integer workspace available for SKSAMPLE
C     LRWORK  maximum real  workspace available for SKSAMPLE
C     LCWORK  maximum character workspace available for SKSAMPLE
C
C     numerical tolerances for solution convergence:
C    
C     ATOL    the absolute tolerance for solution values
C     RTOL    the relative tolerance for solution values
C
C///////////////////////////////////////////////////////////////////////
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER(I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (LIN=5, LOUT=6, LINCK=25, LINSK=26, RTOL=1.0E-6,
     1           ATOL=1.0E-15, LIWORK=4000, LRWORK=4500, LCWORK=200,
     2           ZERO=0.0)
C
      DIMENSION IWORK(LIWORK), RWORK(LRWORK)
      CHARACTER*16 CWORK(LCWORK)
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
C      OPEN(LINSK,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./surf.bin')
C*****end linking files > binary
C*****linking files > ascii
      OPEN(LINCK,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./chem.asc')
      OPEN(LINSK,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./surf.asc')
C*****end linking files > ascii
C
      CALL SKSAMP (LIN, LOUT, LINCK, LINSK, ATOL, RTOL, LIWORK, IWORK,
     1             LRWORK, RWORK, LCWORK, CWORK)
C
      TEND = CKTIME (TSTART)
      IF (TEND .GT. 60.) THEN
         WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (min): ',TEND/60.
      ELSE
          WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (sec): ',TEND
      ENDIF
C
      CLOSE(LINCK)
      CLOSE(LINSK)
      STOP
      END
