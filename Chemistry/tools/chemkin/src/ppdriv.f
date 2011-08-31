C     CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      PROGRAM PSROUT
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for EPOST, the postprocessor for AURORA,
C  the perfectly stirred reactor code..
C  The file is used to make 'epost.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     LIN     user Keyword input
C     LOUT    diagnostic printing (output)
C     LSAVE   the Aurora binary solution Save File (input)
C     LPRNT   postprocessing file output containing colums
C             of solution data as requested by user (output)
C
C     dimensions:
C
C     MXSPEC  maximum total number of species in the Aurora solution
C     MXPTS   maximum number of solutions that can be read in from
C             the Aurora Save File
C     MXREAC  maximum number of gas or surface reactions considered
C     MXPHSE  maximum number of surface phases for each material
C     MAXION  maximum number of ionic species
C     MAXMAT  maximum number of different materials where surface
C             chemistry mechanisms are specified.
C     MAXPSR  maximum number of perfectly stirred reactors in series
C     MAXBUL  maximum number of bulk species
C     MAXPAR  maximum number of parameters the user can request in
C             the solution output
C     LREAC   maximum total of gas+surface reactions on all materials
C     NCMAX   maximum number of columns of data output 
C     LENICK  maximum length of the gas-phase Chemkin integer array
C     LENRCK  maximum length of the gas-phase Chemkin real array
C     LENCCK  maximum length of the gas-phase Chemkin character array
C     LENISK  maximum length of the Surface Chemkin integer array
C     LENRSK  maximum length of the Surface Chemkin real array
C     LENCSK  maximum length of the Surface Chemkin character array
C     LIPAR   maximum integer workspace available for Epost
C     LRPAR   maximum real workspace available for Epost
C     LCPAR   maximum character workspace available for Epost
C     LLPAR   maximum logical workspace available for Epost
C
C///////////////////////////////////////////////////////////////////////
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
      PARAMETER (LIN=5, LOUT=6, LSAVE=45, LPRNT=8, MXSPEC=100, MXPTS=25,
     1           MXREAC=300, MXPHSE=5, MAXION=20, MAXMAT=4, MAXPSR=3,
     2           MAXBUL=5, MAXPAR=14, LREAC=2500, NCMAX=500,
     3           LENICK=8000, LENRCK=5000, LENCCK=100,
     4           LENISK=8000, LENRSK=8000, LENCSK=200,
     5           LIPAR=5000, LRPAR=8000000,  LCPAR=2500, LLPAR=300)
C
      DIMENSION ICKWRK(LENICK), RCKWRK(LENRCK),
     1          ISKWRK(LENISK), RSKWRK(LENRSK),
     2          IPAR(LIPAR),    RPAR(LRPAR)
      CHARACTER*16 CCKWRK(LENCCK), CSKWRK(LENCSK), CPAR(LCPAR)
      LOGICAL LPAR(LLPAR)
      CHARACTER*60 REACTS(LREAC)
C
      OPEN (LSAVE,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE='./save.bin')
C
      CALL PPSR (LIN, LOUT, LSAVE, LPRNT, MXSPEC, MXREAC, MAXMAT,
     1           MAXPSR, MAXBUL, MXPTS, MXPHSE, MAXION, MAXPAR,
     2           LREAC, NCMAX, REACTS, 
     2           LENICK, ICKWRK, LENRCK, RCKWRK, LENCCK, CCKWRK,
     3           LENISK, ISKWRK, LENRSK, RSKWRK, LENCSK, CSKWRK,
     4           LIPAR,  IPAR,   LRPAR,  RPAR,   LCPAR,  CPAR,
     5           LLPAR,  LPAR)
      STOP
      END
