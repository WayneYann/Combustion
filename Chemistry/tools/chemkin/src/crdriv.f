C     CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:32 $
      PROGRAM CRDRIV
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for CRESLAF, the chemically reacting,
C  shear-layer flow code.
C  The file is used to make 'creslaf.exe'
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
C     LSAVE   binary solution Save File (output)
C     LREST   binary Restart File (input)
C     LDAS    binary dump of DASSL solver workspace used internally
C     LPRO    binary file w/ user-supplied initial solution profile 
C             (input)
C
C     dimensions:
C
C     LNIPAR  maximum integer workspace available for Creslaf
C     LNRPAR  maximum real  workspace available for Creslaf
C     LNCWRK  maximum character workspace available for Creslaf
C     LNLWRK  maximum logical workspace available for Creslaf
C
C///////////////////////////////////////////////////////////////////////
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (LNIPAR=5000, LNRPAR=500000, LNCWRK=100, LNLWRK=100,
     1           LIN=5, LOUT=6, LINCK=25, LINSK=26, LINMC=35,
     2           LSAVE=37, LDAS=38, LPRO=33, LREST=14, ZERO=0.0)
C
      DIMENSION IPAR(LNIPAR), RPAR(LNRPAR)
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
      OPEN(LSAVE,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./save.bin')
      OPEN(LREST,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./rest.bin')
      OPEN(LDAS,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./ldas.bin')
      OPEN(LPRO,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE='./cres.pro')
C
      CALL CRES (LIN, LOUT, LINCK, LINSK, LINMC, LREST, LSAVE, LDAS, 
     1           LPRO, LNIPAR, IPAR, LNRPAR, RPAR,
     1           LNCWRK, CWORK, LNLWRK, LWORK)
C
      TEND = CKTIME (TSTART)
      IF (TEND .GT. 60.) THEN
         WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (min): ',TEND/60.
      ELSE
          WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (sec): ',TEND
      ENDIF
C
      CLOSE (LINCK)
      CLOSE(LINMC)
      CLOSE(LINSK)
      CLOSE(LDAS)
      CLOSE(LSAVE)
      CLOSE(LPRO)
      STOP
      END
