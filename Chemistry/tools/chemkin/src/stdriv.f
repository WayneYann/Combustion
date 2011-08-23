C     CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
       PROGRAM STDRIV
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for SURFTHERM, the code for analyzing
C  gas-phase and surface Chemkin mechanisms and thermodynamic data.
C  The file is used to make 'surftherm.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     LOUT    formatted solution and diagnostic printing (output)
C     LCOM    user Keyword input file (input)
C     LINKCK  gas-phase Chemkin Linking File (input)
C     LINKSK  Surface Chemkin Linking File (input)
C     LINKMC  Transport Linking File (input)
C     LIN     This unit is not currently used by Surftherm
C
C     dimensions:
C
C     LIWORK  maximum integer workspace available for Surftherm
C     LRWORK  maximum real  workspace available for Surftherm
C     LCWORK  maximum character workspace available for Surftherm
C     LWRK64
C     LWRK48
C
C  Other parameters controlling Surftherm calculations:
C
C     LTRAN  a logical flag indicating whether or not transport
C            properties will be included as part of the Surftherm
C            calculations.  If LTRAN = .TRUE. the driver will try
C            to open a Transport Linking File, 'tran.asc'.
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
      PARAMETER (LIWORK=8000, LRWORK=8000, LCWORK=500,
     1           LWRK64=100, LWRK48=10, ZERO=0.0)
C
      DIMENSION IWORK(LIWORK), RWORK(LRWORK)
      CHARACTER CWORK(LCWORK)*16, CWRK64(LWRK64)*64, CWRK48(LWRK48)*48
      LOGICAL LCOM_EXIST, LTRAN
      COMMON /SF_UIO/ LIN, LOUT, LCOM, LINKCK, LINKSK, LINKMC,
     1                LCOM_EXIST
      COMMON /SF_KEY/ LTRAN
      EXTERNAL CKTIME
      LTRAN = .TRUE.
C     (set LTRAN .FALSE. if no transport in input)
C
      TSTART = CKTIME(ZERO)
      LIN = 5
      LOUT = 6
      LCOM = 15
      LINKCK=25
      LINKSK=26
      LINKMC=35
C
C///////////////////////////////////////////////////////////////////////
C   Note:  The ascii-formatted linking file option is the recommended 
C          standard, as shown below in the "linking files > ascii"
C          section.  Any change to these open statements requires
C          systematic changes to all of the Chemkin driver routines.
C///////////////////////////////////////////////////////////////////////
C
C*****linking files > binary
C      OPEN(LINKCK,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./chem.bin')
C      OPEN(LINKMC,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./tran.bin')
C      OPEN(LINKSK,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./surf.bin')
C*****end linking files > binary
C*****linking files > ascii
      OPEN(LINKCK,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./chem.asc')
      OPEN(LINKMC,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./tran.asc')
      OPEN(LINKSK,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./surf.asc')
C*****end linking files > ascii
C
      OPEN(LCOM,STATUS='UNKNOWN',FORM='FORMATTED',
     1     FILE='./surftherm.inp')
C
      CALL STHERM (LIWORK, IWORK, LRWORK, RWORK, LCWORK, CWORK,
     1             LWRK64, CWRK64, LWRK48, CWRK48)
C
      TEND = CKTIME (TSTART)
      IF (TEND .GT. 60.) THEN
         WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (min): ',TEND/60.
      ELSE
          WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (sec): ',TEND
      ENDIF
C
      CLOSE(LINKCK)
      CLOSE(LINKMC)
      CLOSE(LCOM)
      CLOSE(LINKSK)
      STOP
      END
