C     CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:32 $
      PROGRAM CPDRIV
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for CONP, the gas-phase Chemkin Example.
C  The file is used to make 'conp.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     LIN     user Keyword input
C     LOUT    solution and diagnostic printing (output)
C     LINCK   gas-phase Chemkin Linking File (input)
C
C     dimensions:
C
C     KMAX    maximum number of species
C     LENIWK  maximum integer workspace available for CONP
C     LENRWK  maximum real workspace available for CONP
C     LENCWK  maximum character workspace available for CONP
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
      PARAMETER (LENIWK=4000, LENRWK=4000, LENCWK=500, LIN=5, LOUT=6,
     1           LINCK=25, KMAX=50, RTOL=1.0E-6, ATOL=1.0E-15, ZERO=0.0)
C
      DIMENSION IWORK(LENIWK), RWORK(LENRWK), X(KMAX), Z(KMAX)
      CHARACTER*16 CWORK(LENCWK), KSYM(KMAX)
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
      CALL CONP (LINCK, LIN, LOUT, LENIWK, LENRWK, LENCWK, KMAX,
     1            IWORK, RWORK, CWORK, KSYM, X, Z, RTOL, ATOL)
C
      TEND = CKTIME (TSTART)
      IF (TEND .GT. 60.) THEN
         WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (min): ',TEND/60.
      ELSE
          WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (sec): ',TEND
      ENDIF
C
      CLOSE (LINCK)
      STOP
      END
