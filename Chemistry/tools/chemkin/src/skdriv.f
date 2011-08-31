C     CVS: $Revision: 1.1.1.1 $ created $Date: 2006/05/26 19:09:33 $

      PROGRAM SKDRIV
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for the Surface Chemkin Interpreter;
C  The file is used to make 'surf.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     unit numbers:
C
C     LIN     surface species and reaction info (input)
C     LINCK   gas-phase Chemkin Linking File (input)
C     LINSK   Surface Chemkin Linking File (output)
C     LOUT    diagnostic printing (output)
C     LTHRM   Chemkin Thermodynamics Database file (input)
C
C     dimensions:
C     dimensions:
C
C     IDIM    maximum number of reactions
C     KDIM    maximum number of species
C     MDIM    maximum number of elements
C     MXPHSE  maximum number of surface phases for each material
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
C             single reaction
C
C  The driver calculates the following maximum sizes for the work arrays
C  needed by the surface chemkin interpreter, based on the dimensions 
C  above:
C
C     LIWORK  max length of the gas-phase chemkin integer workspace
C     LRWORK  max length of the gas-phase chemkin real workspace
C     LCWORK  max length of the gas-phase chemkin character workspace
C     LLWORK  max length of the gas-phase chemkin logical workspace
C     LISWRK  the length of the surface integer work array
C     LRSWRK  the length of the surface real work array
C     LCSWRK  the length of the surface character work array
C     LLSWRK  the length of the surface logical work array
C
C///////////////////////////////////////////////////////////////////////

C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*16 C
      EXTERNAL CKTIME, SKINTP
C      INTEGER LCWORK, LCSWRK, LIWORK, LISWRK, LLWORK, LLSWRK, LRWORK,
C     +   LRSWRK
      LOGICAL L
C
      PARAMETER (LIN = 5, LINCK = 25, LINSK = 26, LOUT = 6, LTHRM = 17)
C
      PARAMETER
     +  (IDIM = 1000,
     +   KDIM = 200,
     +   MDIM = 20,
     +   MAXORD = 10,
     +   MAXSP = 12,
     +   MAXTB = 10,
     +   MAXTP = 3,
     +   MXPHSE = 20)
C
      PARAMETER
     +  (LCWORK = KDIM + MDIM,
     +   LIWORK = IDIM * (27 + MAXORD + 2 * MAXSP + MAXTB)
     +      + KDIM * (4 + MDIM),
     +   LLWORK = KDIM,
     +   LRWORK = IDIM * (33 + MAXORD + MAXSP + MAXTB)
     +      + 2 * KDIM * (4 * MAXTP - 3) + MDIM)
C
      PARAMETER
     +  (LCSWRK = LCWORK + MDIM + KDIM + MXPHSE,
     +   LISWRK = LIWORK + IDIM * (17 + MAXORD + 3 * MAXSP)
     +      + KDIM * (5 + MDIM) + 3 * MXPHSE,
     +   LLSWRK = LLWORK + KDIM + MXPHSE,
     +   LRSWRK = LRWORK + KDIM * (8 * MAXTP - 5) + MDIM + MXPHSE
     +      + IDIM * (17 + MAXORD + MAXSP + 2 * MXPHSE))
C
      DIMENSION C(LCSWRK), I(LISWRK), L(LLSWRK), R(LRSWRK)
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
      OPEN (LINSK, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = './surf.asc')
C*****end linking files > ascii
C*****linking files > binary
C      OPEN (LINCK, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED',
C     +   FILE = './chem.bin')
C      OPEN (LINSK, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED',
C     +   FILE = './surf.bin')
C*****end linking files > binary
C
      OPEN (LTHRM, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = './therm.dat')
C
      CALL SKINTP
     +  (MDIM, KDIM, IDIM, LIN, LOUT, LTHRM, LINCK, LINSK, MAXTP, MAXSP,
     +   MXPHSE, MAXORD, LISWRK, I, LRSWRK, R, LCSWRK, C, LLSWRK, L)
C
      TEND = CKTIME (TSTART)
      IF (60 .LT. TEND) THEN
         TEND = TEND / 60
         WRITE (LOUT, '(A, 1PE15.2)') ' Total CPUtime (min): ', TEND
      ELSE
         WRITE (LOUT, '(A, 1PE15.2)') ' Total CPUtime (sec): ', TEND
      ENDIF
C
      CLOSE (LINCK)
      CLOSE (LINSK)
      CLOSE (LTHRM)
C
      STOP
      END

