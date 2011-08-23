C     CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      PROGRAM PRDRIV
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for PREMIX, the premixed-flame code.
C  The file is used to make 'premix.exe'
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
C     LENLWK  maximum logical workspace available for Premix
C     LENIWK  maximum integer workspace available for Premix
C     LENRWK  maximum real  workspace available for Premix
C     LENCWK  maximum character workspace available for Premix
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
c     Used these for fine uniform drm-noAR run
c      PARAMETER (LENLWK=5034, LENIWK=190531, LENRWK=44658792, LENCWK=103
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=5000, ZERO=0.0)
c      PARAMETER (LENLWK=2067, LENIWK=153184, LENRWK=66609550, LENCWK=202
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=2000, ZERO=0.0)
c      These for Hunter mech
c      PARAMETER (LENLWK=549, LENIWK=34452, LENRWK=9077330, LENCWK=148
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=500, ZERO=0.0)
c      Used these for all the GRI phi runs...
c      PARAMETER (LENLWK=1034, LENIWK=42531, LENRWK=8950792, LENCWK=103
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=500, ZERO=0.0)
c      Used these for all the GRI phi runs...
c      PARAMETER (LENLWK=1034, LENIWK=42531, LENRWK=11360799, LENCWK=166
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=500, ZERO=0.0)
cc      PARAMETER (LENLWK=2055, LENIWK=125993, LENRWK=45260799, LENCWK=166
cc     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
cc     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=2000, ZERO=0.0)
c-commented out 6/2/06 to test large stiff propane mech with parameters in following set
c      PARAMETER (LENLWK=2067, LENIWK=153184, LENRWK=66609550, LENCWK=202
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=2000, ZERO=0.0)
c     3-lean-flames settings
c      PARAMETER (LENLWK=4038, LENIWK=169189, LENRWK=43953437, LENCWK=202
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=4000, ZERO=0.0)
      PARAMETER (LENLWK=4055, LENIWK=241933, LENRWK=90460799, LENCWK=202
     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=4000, ZERO=0.0)
c      Used these for the DME runs...
c      PARAMETER (LENLWK=581, LENIWK=52704, LENRWK=23997821, LENCWK=244
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=500, ZERO=0.0)
c      Used these for the DME runs...
c      PARAMETER (LENLWK=281, LENIWK=27504, LENRWK=9676421, LENCWK=244
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=200, ZERO=0.0)
c     Used these for the refined H runs
c      PARAMETER (LENLWK=2511, LENIWK=35883, LENRWK=2652459, LENCWK=103
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=2500, ZERO=0.0)
c     Used these for the small 2-step runs
c      PARAMETER (LENLWK=259, LENIWK=3241, LENRWK=177461, LENCWK=28
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=250, ZERO=0.0)
c     Used these for the small 2-step runs
c      PARAMETER (LENLWK=284, LENIWK=14781, LENRWK=2255542, LENCWK=103
c     1           ,LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
c     2           LRCVR=16, LINCK=25, LINMC=35, NMAX=250, ZERO=0.0)

      DIMENSION IWORK(LENIWK), RWORK(LENRWK)
      CHARACTER CWORK(LENCWK)*(LENSYM)
      LOGICAL LWORK(LENLWK)
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
C*****end linking files > binary
C*****linking files > ascii
c      OPEN(LIN,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./premix.inp')
      OPEN(LINCK,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./chem.asc')
      OPEN(LINMC,FORM='FORMATTED',STATUS='UNKNOWN',FILE='./tran.asc')
C*****end linking files > ascii
C
      OPEN(LRIN,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./rest.bin')
      OPEN(LROUT,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./save.bin')
      OPEN(LRCVR,FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='./recov.bin')
C
      CALL PREMIX (NMAX, LIN, LOUT, LINCK, LINMC, LRIN, LROUT, LRCVR,
     1             LENLWK, LWORK, LENIWK, IWORK, LENRWK, RWORK, LENCWK,
     2             CWORK)
C
      TEND = CKTIME (TSTART)
      IF (TEND .GT. 60.) THEN
         WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (min): ',TEND/60.
      ELSE
          WRITE (LOUT, '(A,1PE15.2)') ' Total CPUtime (sec): ',TEND
      ENDIF
      CLOSE(LRIN)
      CLOSE(LROUT)
      CLOSE(LRCVR)
      CLOSE(LINCK)
      CLOSE(LINMC)
      STOP
      END
C
C*****precision > double
C
        DOUBLE PRECISION FUNCTION AREA(X)
C
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C
C        REAL FUNCTION AREA(X)
C
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C  START PROLOGUE
C
C  User Function to define flame tube area as a function of position. 
C
C  X - real scalar, distance from the burner in (cm)
C  AREA is the area of the flame tube relative to the burner area
C       as a function of X;  AREA = A(X)/Aburner.  (UNITS = none)
C
C  END PROLOGUE
C
      AREA = 1.0
C
C     end of FUNCTION AREA
      RETURN
      END
