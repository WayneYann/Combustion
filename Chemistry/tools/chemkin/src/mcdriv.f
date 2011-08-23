C     CVS: $Revision: 1.1.1.1 $ created $Date: 2006/05/26 19:09:33 $

      PROGRAM MCDRIV
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for the Chemkin Interpreter;
C  The file is used to make 'chem.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     LINCK   gas-phase Chemkin Linking File (input)
C     LINMC   gas-phase Transport Linking File (output)
C     LTRAN   gas-phase Transport Database file (input)
C     LOUT    diagnostic printing (output)
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
C     LIWORK  max length of the gas-phase chemkin integer workspace
C     LRWORK  max length of the gas-phase chemkin real workspace
C     LCWORK  max length of the gas-phase chemkin character workspace
C     LLWORK  max length of the gas-phase chemkin logical workspace
C     LIMCWK  the length of the transport integer work array
C     LRMCWK  the length of the transport real work array
C     LCMCWK  the length of the transport character work array
C
C///////////////////////////////////////////////////////////////////////

C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single

      CHARACTER*16 C
      EXTERNAL CKTIME, TRANFT
C      INTEGER LCWORK, LCMCWK, LIWORK, LIMCWK, LLWORK, LRWORK, LRMCWK


      PARAMETER (LINCK = 25, LINMC = 35, LOUT = 6, LTRAN = 31)


      PARAMETER
     +  (IDIM = 1000,
     +   KDIM = 200,
     +   MDIM = 20,
     +   MAXORD = 10,
     +   MAXSP = 12,
     +   MAXTB = 10,
     +   MAXTP = 3)

      PARAMETER
     +  (LCWORK = KDIM + MDIM,
     +   LIWORK = IDIM * (27 + MAXORD + 2 * MAXSP + MAXTB)
     +      + KDIM * (4 + MDIM),
     +   LLWORK = KDIM,
     +   LRWORK = IDIM * (33 + MAXORD + MAXSP + MAXTB)
     +      + 2 * KDIM * (4 * MAXTP - 3) + MDIM)

      PARAMETER
     +  (LCMCWK = LCWORK + KDIM,
     +   LIMCWK = LIWORK + 3 * KDIM,
     +   LRMCWK = LRWORK + 300 + KDIM * (48 + 4 * KDIM + 8 * MAXTP))

      DIMENSION C(LCMCWK), I(LIMCWK), R(LRMCWK)

      ZERO = 0
      TSTART = CKTIME (ZERO)
C
C///////////////////////////////////////////////////////////////////////
C   Note:  The ascii-formatted linking file option is the recommended 
C          standard, as shown below in the "linking files > ascii"
C          section.  Any change to these open statements requires
C          systematic changes to all of the Chemkin driver routines.
C///////////////////////////////////////////////////////////////////////

C*****linking files > ascii
      OPEN (LINCK, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = './chem.asc')
      OPEN (LINMC, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = './tran.asc')
C*****end linking files > ascii
C*****linking files > binary
C      OPEN (LINCK, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED',
C     +   FILE = './chem.bin')
C      OPEN (LINMC, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED',
C     +   FILE = './tran.bin')
C*****end linking files > binary

C     OPEN (LDATA, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED',
C    +   FILE = './tran.inp')

      OPEN (LTRAN, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = './tran.dat')

C     SUBROUTINE TRANFT
C    +  (LINCK, LINMC, LTRAN, LOUT, MAXTP, LIWORK, LRWORK, LCWORK, I, R,
C    +   C)

      CALL TRANFT
     +  (LINCK, LINMC, LTRAN, LOUT, MAXTP, LIMCWK, LRMCWK, LCMCWK, I, R,
     +   C)

      TEND = CKTIME (TSTART)
      IF (60 .LT. TEND) THEN
         TEND = TEND / 60
         WRITE (LOUT, '(A, 1PE15.2)') ' Total CPUtime (min): ', TEND
      ELSE
         WRITE (LOUT, '(A, 1PE15.2)') ' Total CPUtime (sec): ', TEND
      ENDIF

      CLOSE (LINCK)
      CLOSE (LINMC)
      CLOSE (LTRAN)

      STOP
      END
