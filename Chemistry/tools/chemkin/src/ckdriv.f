C     CVS: $Revision: 1.1.1.1 $ created $Date: 2006/05/26 19:09:32 $
      PROGRAM CKDRIV
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for the Chemkin Interpreter;
C  The file is used to make 'chem.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     LIN     gas-phase species and reaction info (input)
C     LINCK   gas-phase Chemkin Linking File (output)
C     LOUT    diagnostic printing (output)
C     LTHRM   Chemkin Thermodynamics Database file (input)
C
C     dimensions:
C
C     IDIM    maximum number of reactions
C     KDIM    maximum number of species
C     MDIM    maximum number of elements
C
C     MAXORD  maximum number of times in a reaction that 
C             the order of a species in the reaction rate is
C             changed to an arbitrary value (non-mass-action kinetics)
C             (This value should not exceed MAXSP=12)
C     MAXTB   maximum number of auxiliary 3rd-body enhancement entries
C             allowd for a single reaction
C             (Usually a value of 10 is sufficient)
C     MAXTP   maximum number of temperatures used to bound temperature
C             ranges for thermodynamic property fits of each species
C             (Usually this is 3, allowing two temperature ranges)
C
C  The following parameter should *NOT* be changed from its
C  value of 12:
C
C     MAXSP   maximum number of reactants + products in a 
C             single reaction = 12
C
C  The driver calculates the following maximum sizes for the work arrays
C  needed by the chemkin interpreter, based on the dimensions above:
C
C     LIWORK  the length of the integer work array
C     LRWORK  the length of the real work array
C     LCWORK  the length of the character work array
C     LLWORK  the length of the logical work array
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
      CHARACTER*16 C
      EXTERNAL CKINTP, CKTIME
C      INTEGER LCWORK, LIWORK, LRWORK, LLWORK
C
      PARAMETER (LIN = 5, LINCK = 25, LOUT = 6, LTHRM = 17)
C
      PARAMETER
     +  (IDIM = 1000,
     +   KDIM = 200,
     +   MDIM = 20,
     +   MAXORD = 10,
     +   MAXSP = 12,
     +   MAXTB = 10,
     +   MAXTP = 3)
C
      PARAMETER
     +  (LCWORK = KDIM + MDIM,
     +   LIWORK = IDIM * (27 + MAXORD + 2 * MAXSP + MAXTB)
     +      + KDIM * (4 + MDIM),
     +   LLWORK = KDIM,
     +   LRWORK = IDIM * (33 + MAXORD + MAXSP + MAXTB)
     +      + 2 * KDIM * (4 * MAXTP - 3) + MDIM)
C
      DIMENSION C(LCWORK), I(LIWORK), L(LLWORK), R(LRWORK)
C
      ZERO = 0
      TSTART = CKTIME (ZERO)
C
C///////////////////////////////////////////////////////////////////////
C   Note:  The ascii-formatted linking file option is the recommended 
C          standard, as shown below in the "linking files > ascii"
C          section.  Any change to these open statements requires
C          systematic changes to all of the Chemkin driver routines.
C///////////////////////////////////////////////////////////////////////
C
C*****linking files > ascii
      OPEN (LINCK, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = './chem.asc')
C*****end linking files > ascii
C*****linking files > binary
C      OPEN (LINCK, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED',
C     +   FILE = './chem.bin')
C*****end linking files > binary
C
      OPEN (LTHRM, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = './therm.dat')
C
      CALL CKINTP
     +  (MDIM, KDIM, IDIM, MAXTP, MAXSP, MAXTB, LIN, LOUT, LTHRM, LINCK,
     +   MAXORD, LIWORK, I, LRWORK, R, LCWORK, C, L)

      TEND = CKTIME (TSTART)
      IF (60 .LT. TEND) THEN
         TEND = TEND / 60
         WRITE (LOUT, '(A, 1PE15.2)') ' Total CPUtime (min): ', TEND
      ELSE
         WRITE (LOUT, '(A, 1PE15.2)') ' Total CPUtime (sec): ', TEND
      ENDIF

      CLOSE (LINCK)
      CLOSE (LTHRM)

      STOP
      END
