C     CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      PROGRAM SENKOUT
C///////////////////////////////////////////////////////////////////////
C  This is the driver routine for PSENK, the postprocessor for SENKIN,
C  the sensitivity analysis code.
C  The file is used to make 'psenk.exe'
C
C  The parameters and unit numbers that may be changed
C  by the user are described below:
C
C     unit numbers:
C
C     LIN     user Keyword input
C     LOUT    diagnostic printing (output)
C     LSAVE   the Senkin binary solution Save File (input)
C     LTXT    postprocessing file output containing colums
C             of solution data as requested by user (output)
C
C     dimensions:
C
C     KMAX    maximum total number of species in the Senkin solution
C     JMAX    maximum number of points that can be read in from
C             the Senkin Save File
C     IMAX    maximum number of gas or surface reactions considered
C     LENICK  maximum length of the gas-phase Chemkin integer array
C     LENRCK  maximum length of the gas-phase Chemkin real array
C     LENCCK  maximum length of the gas-phase Chemkin character array
C     LNIPAR  maximum integer workspace available for PSENK
C     LNRPAR  maximum real workspace available for PSENK
C     LNCPAR  maximum character workspace available for PSENK
C
C///////////////////////////////////////////////////////////////////////
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (LIN=5, LOUT=6, LSAVE=55, LTXT=7, KMAX=50, JMAX=300,
     1           IMAX=100, LENICK=5000, LENRCK=5000, LENCCK=250,
     2           LNIPAR=5000, LNRPAR=150000, LNCPAR=250)
C
      DIMENSION ICKWRK(LENICK), RCKWRK(LENRCK),
     1          IPAR(LNIPAR), RPAR(LNRPAR)
      CHARACTER CCKWRK(LENCCK)*16, CPAR(LNCPAR)*16, REACTS(IMAX)*32
C
      CALL PSENK (LIN, LOUT, LSAVE, LTXT, KMAX, JMAX, IMAX,
     1            LENICK, ICKWRK, LENRCK, RCKWRK, LENCCK, CCKWRK,
     2            LNIPAR, IPAR, LNRPAR, RPAR, LNCPAR, CPAR, REACTS)
C
      CLOSE (LSAVE)
      CLOSE (LTXT)
      STOP
      END
