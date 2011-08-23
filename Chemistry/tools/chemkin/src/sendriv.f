C     CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      PROGRAM DRIVER
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for SENKIN, the sensitivity analysis code.
C  The file is used to make 'senkin.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     LIN     user Keyword input
C     LOUT    formatted solution and diagnostic printing (output)
C     LINCK   gas-phase Chemkin Linking File (input)
C     LREST   binary restart file (input)
C     LSAVE   binary solution Save File (output)
C
C     dimensions:
C
C     LENIWK  maximum integer workspace available for Senkin
C     LENRWK  maximum real  workspace available for Senkin
C     LENCWK  maximum character workspace available for Senkin
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
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      PARAMETER (LENIWK=7500, LENRWK=60000, LENCWK=500, LENSYM=16,
     1           LIN=5, LOUT=6, LSAVE=7, LREST=10, LINCK=25, ZERO=0.0)
      DIMENSION IWORK (LENIWK), RWORK (LENRWK)
      CHARACTER CWORK(LENCWK)*(LENSYM)
C
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
C*****end linking files > binary
C*****linking files > ascii
      OPEN(LINCK,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./chem.asc')
C*****end linking files > ascii
C
      OPEN(LSAVE,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./save.bin')
      OPEN(LREST,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./rest.bin')
C
      CALL SENKIN (LIN, LOUT, LINCK, LSAVE, LREST, LENRWK, RWORK,
     1             LENIWK, IWORK, LENCWK, CWORK)
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
      CLOSE (LREST)
      END
C
      SUBROUTINE TEMPT (A, B, C)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      RETURN
      END
C
      SUBROUTINE VOLT (A, B, C)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      RETURN
      END
