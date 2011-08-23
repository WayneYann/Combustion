C     CVS $Revision: 1.1.1.1 $ created $Date: 2006/05/26 19:09:33 $
      PROGRAM OPPDR
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for OPPDIF, the opposed-flow flame code.
C  The file is used to make 'oppdif.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     LIN     user Keyword input
C     LOUT    formatted solution and diagnostic printing (output)
C     LINCK   gas-phase Chemkin Linking File (input)
C     LINMC   Transport Linking File (input)
C     LRIN    binary restart file (input)
C     LROUT   binary solution Save File (output)
C     LRCVR   binary solution Recover File (output)
C
C     dimensions:
C
C     LENLWK  maximum logical workspace available for Oppdif
C     LENIWK  maximum integer workspace available for Oppdif
C     LENRWK  maximum real  workspace available for Oppdif
C     LENCWK  maximum character workspace available for Oppdif
C     NMAX    total number of grid points allowed for grid refinement
C
C  The following parameter should *NOT* be changed from its
C  value of 16:
C
C     LENSYM  the number of characters in a Chemkin character string
C              = 16
C
C///////////////////////////////////////////////////////////////////////
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (LENLWK=200, LENIWK=20500, LENRWK=2359000, LENCWK=200,
     1           LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=90)
C
      DIMENSION IWORK(LENIWK), RWORK(LENRWK)
      CHARACTER CWORK(LENCWK)*(LENSYM)
      LOGICAL LWORK(LENLWK)
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
C*****end linking files > binary
C*****linking files > ascii
      OPEN(LINCK,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./chem.asc')
      OPEN(LINMC,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./tran.asc')
C*****end linking files > ascii
C
      OPEN(LRIN,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./rest.bin')
      OPEN(LROUT,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./save.bin')
      OPEN(LRCVR,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./recov.bin')
C
      CALL OPPDIF
     1 (LIN, LOUT, LINCK, LINMC, LRIN, LROUT, LRCVR, NMAX,
     2  LENCWK, LENIWK, LENLWK, LENRWK, CWORK, IWORK, LWORK, RWORK)
C
      CLOSE(LRIN)
      CLOSE(LROUT)
      CLOSE(LRCVR)
      CLOSE(LINCK)
      CLOSE(LINMC)
C
      STOP
      END
