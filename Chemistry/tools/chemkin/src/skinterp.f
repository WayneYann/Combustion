C    CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
C
      SUBROUTINE SKINTP (MDIM, KDIM, IDIM, LIN, LOUT, LTHRM, LINCK,
     1                   LINSK, MAXTP, MAXSP, MXPHSE, MAXORD, LIWORK,
     3                   I, LRWORK, R, LCWORK, C, LLWORK, L)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (NSPAR=3, NCP=5, NFIT=NCP+2, NSCOV=3, NEDPAR=3,
     1           NYPAR=4)
      DIMENSION I(LIWORK), R(LRWORK)
      CHARACTER*16 C(LCWORK), PRVERS, FILVER, PRDATE, PREC
      LOGICAL L(LLWORK), KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
      COMMON /SKVERS/ PRVERS, PREC, FILVER
C
      NTR = MAXTP - 1
C----------------------------------------------------------------------C
C
C     Write to standard output information about the skinterp program
C     and linkfile.
C
      PRVERS = '7.10'
      FILVER= '1.0'
      PRDATE = '98/03/25'
C*****precision > double
      PREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C      PREC = 'SINGLE'
C*****END precision > single
C
      WRITE (LOUT, '(/A, /1X,A, A, A, A, /A, /A, //)')
     1   ' CHEMKIN-III SURFACE MECHANISM INTERPRETER:',
     2    PREC(1:CKLSCH(PREC)), ' PRECISION Vers. ',
     3    PRVERS(1:CKLSCH(PRVERS)+1), PRDATE,
     4   ' Copyright 1995, Sandia Corporation.',
     5' The U.S. Government retains a limited license in this software.'
C
      KERR = .FALSE.
C     Integer arrays
      IICK = 1
      CALL CKLEN (LINCK, LOUT, LIWK, LRWK, LCWK, IFLAG)
      IF (IFLAG .GT. 0) THEN
         KERR = .TRUE.
         RETURN
      ENDIF
      INT    = IICK   + LIWK
      IKCHRG = INT    + KDIM
      IKPHSE = IKCHRG + KDIM
      IKNCF  = IKPHSE + KDIM
      IKCOV  = IKNCF  + MDIM * KDIM
      IKFIRS = IKCOV  + KDIM
      IKLAST = IKFIRS + MXPHSE
      IKKPHA = IKLAST + MXPHSE
      INSPEC = IKKPHA + MXPHSE
      INU    = INSPEC + IDIM
      INUNK  = INU    + MAXSP*IDIM
      INREAC = INUNK  + MAXSP*IDIM
      INUSUM = INREAC + IDIM
      IIDUP  = INUSUM + IDIM
      IICOV  = IIDUP  + IDIM
      IIKCOV = IICOV  + IDIM
      IIREV  = IIKCOV + IDIM
      IISTK  = IIREV  + IDIM
      IIMOTZ = IISTK  + IDIM
      IIRNU  = IIMOTZ + IDIM
      IIORD  = IIRNU  + IDIM
      IKORD  = IIORD  + IDIM
      IIBHM  = IKORD  + MAXORD * IDIM
      IKBHM  = IIBHM  + IDIM
      IKION  = IKBHM  + IDIM
      IIEDP  = IKION  + KDIM
      IKEDP  = IIEDP  + IDIM
      IIYLD  = IKEDP  + IDIM
      IKYION = IIYLD  + IDIM
      IKYLD  = IKYION + IDIM
      IITOT  = IKYLD  + MAXSP * IDIM - 1
      IF (LIWORK .LT. IITOT) THEN
         WRITE (LOUT, *)
     1   'SKINTP ERROR: IWORK needs to be at least ', IITOT
         KERR = .TRUE.
      ENDIF
C     Real arrays
      IRCK   = 1
      ITMP   = IRCK   + LRWK
      IA     = ITMP   + MAXTP * KDIM
      IAWT   = IA     + NFIT * NTR * KDIM
      IWT    = IAWT   + MDIM
      IDEN   = IWT    + KDIM
      IPDEN  = IDEN   + KDIM
      IRNCF  = IPDEN  + MXPHSE
      ISPAR  = IRNCF  + MXPHSE * IDIM
      ICPAR  = ISPAR  + NSPAR  * IDIM
      IRSPAR = ICPAR  + NSCOV  * IDIM
      IRNU   = IRSPAR + NSPAR  * IDIM
      IRORD  = IRNU   + MAXSP  * IDIM
      IPEDP  = IRORD  + MAXORD * IDIM
      IPYLD  = IPEDP  + NEDPAR * IDIM
      IYNCF  = IPYLD  + NYPAR  * IDIM
      IEQFAC = IYNCF  + MXPHSE * IDIM
      IRTOT  = IEQFAC + IDIM - 1
      IF (LRWORK .LT. IRTOT) THEN
         WRITE (LOUT, *)
     1   'SKINTP ERROR: RWORK needs to be at least ', IRTOT
         KERR = .TRUE.
      ENDIF
      ICCK = 1
      ICKNAM = ICCK + LCWK
      ICENAM = ICKNAM + KDIM
      ICPNAM = ICENAM + MDIM
      ICTOT  = ICPNAM + MXPHSE - 1
      IF (LCWORK .LT. ICTOT) THEN
         WRITE (LOUT, *)
     1   'SKINTP ERROR: CWORK needs to be at least ', ICTOT
         KERR = .TRUE.
      ENDIF
      LITHRM = 1
      LIPDEN = LITHRM + KDIM
      LIKDEN = LIPDEN + MXPHSE
      LTOT   = LIKDEN + KDIM - 1
      IF (LLWORK .LT. LTOT) THEN
         WRITE (LOUT, *)
     1   'SKINTP ERROR: LWORK needs to be at least ', LTOT
         KERR = .TRUE.
      ENDIF
      IF (KERR) RETURN
C
      CALL SKMECH
     1(MDIM, KDIM, IDIM, LIWK, LRWK, LCWK, LIN, LOUT, LTHRM, LINCK,
     2 LINSK, NSPAR, NCP, NFIT, MAXTP, NTR, MAXSP, NSCOV, NEDPAR,
     3 NYPAR, MXPHSE, MAXORD, I(IICK), I(INT), I(IKCHRG), I(IKPHSE),
     4 I(IKNCF), I(IKCOV), I(IKFIRS), I(IKLAST), I(IKKPHA),
     5 I(INSPEC), I(INU), I(INUNK), I(INREAC), I(INUSUM), I(IIDUP),
     6 I(IICOV), I(IIKCOV), I(IIREV), I(IISTK), I(IIMOTZ), I(IIRNU),
     7 I(IIORD), I(IKORD), I(IIBHM), I(IKBHM), I(IKION), I(IIEDP),
     8 I(IKEDP), I(IIYLD), I(IKYION), I(IKYLD), R(IRCK), R(ITMP),
     9 R(IA), R(IAWT), R(IWT), R(IDEN), R(IPDEN), R(IRNCF),
     * R(ISPAR), R(ICPAR), R(IRSPAR), R(IRNU), R(IRORD), R(IPEDP),
     1 R(IPYLD), R(IYNCF), R(IEQFAC), C(ICKNAM), C(ICENAM), C(ICCK),
     2 C(ICPNAM), L(LITHRM), L(LIPDEN), L(LIKDEN))
      RETURN
      END

      SUBROUTINE SKMECH
     1(MDIM, KDIM, IDIM, LIWK, LRWK, LCWK, LIN, LOUT, LTHRM, LINCK,
     2 LINSK, NSPAR, NCP, NFIT, MAXTP, NTR, MAXSPR, NSCOV, NEDPAR,
     3 NYPAR, MXPHSE, MAXORD, IWORK, NT, KCHRG, KPHSE, KNCF, KCOV,
     4 KFIRST, KLAST, KKPHAS, NSPEC, NU, NUNK, NREAC, NUSUMK, IDUP,
     5 ICOV, IKCOV, IREV, ISTK, MOTZ, IRNU, IORD, KORD, IBHM, KBHM,
     6 KION, IEDP, KEDP, IYLD, KYION, KYLD, WORK, TMP, A, AWT, WT,
     7 DEN, PDEN, RNCF, SPAR, CPAR, RSPAR, RNU, RORD, PEDP, PYLD,
     8 YNCF, EQFAC, KNAME, ENAME, CWORK, PNAME, ITHRM, LPDEN, LKDEN)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (SKMIN=0.001)
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
      COMMON /SKVERS/ PRVERS, PREC, FILVER
C
C     Integer arrays
      DIMENSION IWORK(LIWK), NT(KDIM), KCHRG(KDIM), KPHSE(KDIM),
     1          KNCF(MDIM,KDIM), KCOV(KDIM), KFIRST(MXPHSE),
     2          KLAST(MXPHSE), KKPHAS(MXPHSE), NSPEC(IDIM),
     3          NU(MAXSPR,IDIM), NUNK(MAXSPR,IDIM),
     4          NREAC(IDIM), NUSUMK(IDIM), IDUP(IDIM),
     5          ICOV(IDIM), IKCOV(IDIM), IREV(IDIM),
     6          ISTK(IDIM), MOTZ(IDIM), IRNU(IDIM),
     7          IORD(IDIM), KORD(MAXORD,IDIM), IBHM(IDIM),
     8          KBHM(IDIM), KION(KDIM), IEDP(IDIM), KEDP(IDIM),
     9          IYLD(IDIM), KYION(IDIM), KYLD(MAXSPR,IDIM)
C     Real arrays
      DIMENSION WORK(LRWK), TMP(MAXTP,KDIM), A(NFIT, NTR, KDIM),
     1          AWT(MDIM), WT(KDIM), DEN(KDIM), PDEN(MXPHSE),
     2          RNCF(MXPHSE,IDIM), SPAR(NSPAR,IDIM),
     3          CPAR(NSCOV,IDIM), RSPAR(NSPAR,IDIM),
     4          RNU(MAXSPR,IDIM), RORD(MAXORD,IDIM),
     5          PEDP(NEDPAR,IDIM), PYLD(NYPAR,IDIM),
     6          YNCF(MXPHSE,IDIM), EQFAC(IDIM)
C
      CHARACTER*16 MATNAM, KNAME(KDIM), ENAME(MDIM), CWORK(LCWK),
     1             PNAME(MXPHSE), PRVERS, PREC, FILVER
C
      LOGICAL KERR, ITHRM(KDIM), NONCON, LPDEN(MXPHSE), LKDEN(KDIM),
     1        LMOTZ
C
C----------------------------------------------------------------------C
C
C     Keep a record of which "MATErial" (separate surface chemical
C     mechanism) we are currently processing.
C
      NMAT = 0
  100 CONTINUE
C
C     Increment the number of the current material being processed
C
      NMAT = NMAT + 1
C
C     Set relevant variables to appropriate defaults
C
      CALL SKSET
     1   (LRWK, WORK, LIWK, IWORK, LCWK, CWORK, MDIM, KDIM, IDIM,
     2    MXPHSE, MAXSPR, MAXTP, NFIT, NTR, NSPAR, NSCOV, NEDPAR,
     3    NYPAR, MAXORD, KION, NT, TMP, A, KCHRG, KPHSE, AWT, WT,
     4    DEN, KNCF, KCOV, KFIRST, KLAST, KKPHAS,PDEN, NSPEC, NU,
     5    NUNK, RNCF, NREAC, NUSUMK, SPAR, IDUP, ICOV, IKCOV, CPAR,
     6    IREV, RSPAR, ISTK, MOTZ, IRNU, RNU, IORD, KORD, RORD, IBHM,
     7    KBHM, IEDP, KEDP, PEDP, KERR, ITHRM, NONCON, LPDEN, LKDEN,
     8    LMOTZ, KNAME, ENAME, PNAME, IYLD, KYION, KYLD, PYLD, YNCF,
     9    EQFAC)
C
      IF (NMAT .EQ. 1) THEN
C
C        we are working on the first material, so we need to get
C        gas-phase CHEMKIN information
C
C        initialize gas-phase chemistry
C
C*****linkfile (gas) > binary
C         READ (LINCK, END=111)
C*****END linkfile (gas) > binary
C*****linkfile (gas) > ascii
         READ (LINCK, *, END=111)
C*****END linkfile (gas) > ascii
C
C        rewind the unit so that we can read linkfile from
C        the beginning
         REWIND (LINCK)
C
C        read CHEMKIN linkfile and fill up many arrays
C        containing gas-phase mechanism information
         CALL SKBGIN (LIWK, LRWK, LCWK, LINCK, LOUT, MDIM, KDIM,
     1                MAXTP, NFIT, NTR, IWORK, WORK, CWORK, AWT,
     2                ENAME, KNAME, WT, ITHRM, KCHRG, KPHSE, NT,
     3                TMP, KNCF, KION, A, KERR)
         IF (KERR) RETURN
C
      ENDIF
C
C     Interpret surface-phase mechanism
C
      CALL SKKEY
     1   (LIN, LTHRM, LOUT, MDIM, KDIM, IDIM, MXPHSE, MAXSPR, MAXTP,
     2    NFIT, NSPAR, NSCOV, NEDPAR, NYPAR, MAXORD, ENAME, AWT,
     3    KNAME, ITHRM, WT, KNCF, NONCON, KCOV, KPHSE, KCHRG, NT,
     4    NTR, TMP, A, PNAME, PDEN, LPDEN, KFIRST, KLAST, KKPHAS,
     5    DEN, LKDEN, SPAR, NSPEC, NREAC, NU, NUNK, NUSUMK, RNCF,
     6    IDUP, IREV, RSPAR, ICOV, IKCOV, CPAR, ISTK, KERR, IRNU,
     7    RNU, SKMIN, IORD, KORD, RORD, IBHM, KBHM, LMOTZ, MOTZ,
     8    IEDP, KEDP, PEDP, IYLD, KYION, KYLD, PYLD, YNCF, MATNAM,
     9    EQFAC)
C
      MXTP = 0
      DO 120 K = 1, KKTOT
         MXTP = MAX (MXTP, NT(K))
  120 CONTINUE
      CALL SKSIZE (NSPAR, MAXSPR, MXTP, MXTP-1, NFIT, NSCOV, NEDPAR,
     1             NYPAR, MAXORD)
C
C     Write all of the material information to linkfile
C*****linkfile (surface) > binary
C      CALL SKBIN
C*****END linkfile (surface) > binary
C*****linkfile (surface) > ascii
      CALL SKFORM
C*****END linkfile (surface) > ascii
     1   (LINSK, LOUT, SKMIN, MDIM, KDIM, IDIM, MXPHSE, KERR, MAXSPR,
     2    MAXTP, NCP, NSPAR, NSCOV, NEDPAR, NYPAR, MAXORD, MATNAM,
     3    ENAME, AWT, KNAME, WT, KNCF, KCHRG, NT, KPHSE, KCOV, DEN, TMP,
     4    NFIT, NTR, A, KION, PNAME, KFIRST, KLAST, KKPHAS, PDEN, NSPEC,
     5    NREAC, NU, NUNK, NUSUMK, SPAR, RNCF, ICOV, IKCOV, CPAR, IREV,
     6    RSPAR, ISTK, MOTZ, IBHM, KBHM, IRNU, RNU, IORD, KORD, RORD,
     7    IEDP, KEDP, PEDP, IYLD, KYION, KYLD, YNCF, PYLD, EQFAC)
C
C     tell user how much integer, real and character workspace
C     will be needed by the SURFACE CHEMKIN library
C
      WRITE (LOUT, '(/A,3(/A,I6))')
     1   ' WORKING SPACE REQUIREMENTS ARE',
     2   '    INTEGER:   ',LENISK,
     3   '    REAL:      ',LENRSK,
     4   '    CHARACTER: ',LENCSK
C
C     if another material is to be processed, go back to the top
C
      IF (MORE .GT. 0) GO TO 100
C
C     have finished processing all materials; close all files and quit
C
      RETURN
C
C     we branch to these statements upon errors in reading the
C     CHEMKIN linkfile and/or the interpreter input file
C
  111 CONTINUE
      WRITE (LOUT,*) ' Error...cannot read chem.bin...'
      RETURN
      END
C                                                                      C
      SUBROUTINE SKBIN
     1   (LINSK, LOUT, SKMIN, MDIM, KDIM, IDIM, MXPHSE, KERR, MAXSPR,
     2    MAXTP, NCP, NSPAR, NSCOV, NEDPAR, NYPAR, MAXORD, MATNAM,
     3    ENAME, AWT, KNAME, WT, KNCF, KCHRG, NT, KPHSE, KCOV, DEN, TMP,
     4    NFIT, NTR, A, KION, PNAME, KFIRST, KLAST, KKPHAS, PDEN, NSPEC,
     5    NREAC, NU, NUNK, NUSUMK, SPAR, RNCF, ICOV, IKCOV, CPAR, IREV,
     6    RSPAR, ISTK, MOTZ, IBHM, KBHM, IRNU, RNU, IORD, KORD, RORD,
     7    IEDP, KEDP, PEDP, IYLD, KYION, KYLD, YNCF, PYLD, EQFAC)
C
C  START PROLOGUE
C
C  SUBROUTINE SKBIN
C  writes mechanism data to unformatted surface linkfile.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
      COMMON /SKVERS/ PRVERS, PREC, FILVER
C
C     Integer arrays
      DIMENSION KNCF(MDIM,KDIM), KCHRG(KDIM), NT(KDIM), KPHSE(KDIM),
     1          KCOV(KDIM), KFIRST(MXPHSE), KLAST(MXPHSE),
     2          KKPHAS(MXPHSE), KION(KDIM), NSPEC(IDIM),
     3          NREAC(IDIM), NUSUMK(IDIM), NU(MAXSPR, IDIM),
     4          NUNK(MAXSPR, IDIM), IREV(IDIM), ICOV(IDIM),
     5          IKCOV(IDIM), ISTK(IDIM), MOTZ(IDIM), IBHM(IDIM),
     6          KBHM(IDIM), IRNU(IDIM), IORD(IDIM), KORD(MAXORD,IDIM),
     7          IEDP(IDIM), KEDP(IDIM), IYLD(IDIM), KYION(IDIM),
     8          KYLD(MAXSPR,IDIM)
C     Real arrays
      DIMENSION AWT(MDIM), WT(KDIM), DEN(KDIM), TMP(MAXTP,KDIM),
     1          A(NFIT,NTR,KDIM), PDEN(MXPHSE), SPAR(NSPAR,IDIM),
     2          RNCF(MXPHSE,IDIM), CPAR(NSCOV,IDIM),
     3          RSPAR(NSPAR,IDIM), RNU(MAXSPR,IDIM), RORD(MAXORD,IDIM),
     3          PEDP(NEDPAR,IDIM), YNCF(MXPHSE,IDIM), PYLD(NYPAR,IDIM),
     5          EQFAC(IDIM)
C
      CHARACTER*16 MATNAM, ENAME(MDIM), KNAME(KDIM), PNAME(MXPHSE),
     1             PRVERS, PREC, FILVER
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C     Linkfile version (should be 1.0)
      WRITE (LINSK, ERR=5003) FILVER
C
C     Skinterp program version
      WRITE (LINSK, ERR=5003) PRVERS
C
C     Precision
      WRITE (LINSK, ERR=5003) PREC
C
C     Whether skinterp successfully completed or not
      WRITE (LINSK, ERR=5001) KERR
C
      IF (KERR) THEN
         WRITE (LOUT, '(//A,/A)')
     1   'WARNING...THERE IS AN ERROR IN THE SURFACE LINKFILE',
     2   '          DUE TO SURFACE MECHANISM INPUT ERRORS.'
         RETURN
      ENDIF
C
C     Work space sizes
      WRITE (LINSK, ERR=5001) LENISK, LENRSK, LENCSK
C
C     Constants
      MXTK = 0
      DO 120 K = 1, KKTOT
         MXTK = MAX (MXTK, NT(K))
  120 CONTINUE
      WRITE (LINSK, ERR=5001) MAXSPR, MXTK, NCP, NSPAR, NSCOV,
     1                         NEDPAR, NYPAR, MAXORD
C
C     Problem specifications
      WRITE (LINSK, ERR=5001) NELEM, KKGAS, KKSURF, KKBULK, KKTOT,
     1                         NPHASE, NFSURF, NLSURF, NNSURF, NFBULK,
     2                         NLBULK, NNBULK, NIISUR, NCOV, NREV,
     3                         NSTK, NCON, NBHM, NRNU, NORD, NEDP,
     4                         NYLD, NKION, KELECT, MORE
C
C     Elemental, charge, and site balance tolerance parameter
      WRITE (LINSK, ERR=5002) SKMIN
C
C     material name for this surface
      WRITE (LINSK, ERR=5003) MATNAM
C
C     element names
      WRITE (LINSK, ERR=5003) (ENAME(N), N = 1, NELEM)
C
C     elemental atomic weights
      WRITE (LINSK, ERR=5002) (AWT(N), N = 1, NELEM)
C
C     species names
      WRITE (LINSK, ERR=5003) (KNAME(K), K = 1, KKTOT)
C
C     species molecular weights
      WRITE (LINSK, ERR=5002) (WT(K), K = 1, KKTOT)
C
C     species elemental compositions
      WRITE (LINSK, ERR=5001) ((KNCF(N,K), N=1,NELEM), K=1,KKTOT)
C
C     species electronic charges
      WRITE (LINSK, ERR=5001) (KCHRG(K), K = 1, KKTOT)
C
C     the number of temperatures regions used to fit species
C     thermodynamics
      WRITE (LINSK, ERR=5001) (NT(K), K = 1, KKTOT)
C
C     an integer representing physical state for each species
      WRITE (LINSK, ERR=5001) (KPHSE(K), K = 1, KKTOT)
C
C     site coverage of each species
      WRITE (LINSK, ERR=5001) (KCOV(K), K = 1, KKTOT)
C
C     Densities for each species
      WRITE (LINSK, ERR=5002) (DEN(K), K = 1, KKTOT)
C
C     array of temperatures used to specify boundaries of
C     temperature regions used in fits to thermodynamic functions
      WRITE (LINSK, ERR=5002) ((TMP(L,K), L=1,MXTK), K=1,KKTOT)
C
C     fits to thermodynamic functions for each species;
C     NTR=MAXTP-1, NFIT=NCP+2
      WRITE (LINSK, ERR=5002)
     1      (((A(N,L,K), N=1,NFIT), L=1,MXTK-1), K=1,KKTOT)
C
C     Integer species indices for ionic species
      IF (NKION .GT. 0) THEN
         WRITE (LINSK, ERR=5001) NKION
         WRITE (LINSK, ERR=5001) (KION(N), N = 1, NKION)
      ENDIF
C
C     Phase names
      WRITE (LINSK, ERR=5003) (PNAME(N), N=1,NPHASE)
C
C     Phase species indices
      WRITE (LINSK, ERR=5001) (KFIRST(N), N=1,NPHASE)
      WRITE (LINSK, ERR=5001) (KLAST(N),  N=1,NPHASE)
      WRITE (LINSK, ERR=5001) (KKPHAS(N), N=1,NPHASE)
C
C     Phase densities
      WRITE (LINSK, ERR=5002) (PDEN(N), N = 1, NPHASE)
C
      IF (NIISUR .LE. 0) THEN
         WRITE (LOUT, '(/A,/A)')
     1   ' WARNING...NO SURFACE REACTIONS FOUND, ',
     2   ' LINKFILE HAS NO REACTION INFORMATION ON IT.'
         GO TO 500
      ENDIF
C
C     SURFACE REACTION OUTPUT
C     count of species defined as reactants and products
      WRITE (LINSK, ERR=5001) (NSPEC(I), I=1,NIISUR)
C     count of reactants defined for each reaction, also,
C     a negative number flags an irreversible reaction
      WRITE (LINSK, ERR=5001) (NREAC(I), I=1,NIISUR)
C     stoichiometric information for each reaction
      WRITE (LINSK, ERR=5001)
     1      ((NU(N,I), NUNK(N,I), N=1,MAXSPR), I=1,NIISUR)
C     stoichiometric coefficients totals
      WRITE (LINSK, ERR=5001) (NUSUMK(I), I=1,NIISUR)
C     forward Arrhenius coefficients for each reaction
      WRITE (LINSK, ERR=5002) ((SPAR(N,I),N=1,NSPAR),I=1,NIISUR)
C     multiplicative equilibrium factor
      WRITE (LINSK, ERR=5002) (EQFAC(I), I=1,NIISUR)
C     phase balance for each reaction
      WRITE (LINSK, ERR=5002) ((RNCF(N,I),N=1,NPHASE),I=1,NIISUR)
C
C     SPECIAL REACTIONS:
C
      IF (NCOV .GT. 0) THEN
C        coverage reactions
         WRITE (LINSK, ERR=5001) NCOV, NSCOV
C        coverage reaction indices
         WRITE (LINSK, ERR=5001) (ICOV(N), N = 1, NCOV)
C        coverage species indices
         WRITE (LINSK, ERR=5001) (IKCOV(N), N = 1, NCOV)
C        coverage expression parameters
         WRITE (LINSK, ERR=5002) ((CPAR(L,N),L=1,NSCOV),N=1,NCOV)
      ENDIF
C
      IF (NREV .GT. 0) THEN
C        reactions with explicit reverse rate constants
         WRITE (LINSK, ERR=5001) NREV, NSPAR
C        reaction indices
         WRITE (LINSK, ERR=5001) (IREV(N), N=1,NREV)
C        Arrhenius reverse rate constant coefficients
         WRITE (LINSK, ERR=5002) ((RSPAR(L,N),L=1,NSPAR), N=1,NREV)
      ENDIF
C
      IF (NSTK .GT. 0) THEN
C        Reactions for which parameters act as sticking coefficients
         WRITE (LINSK, ERR=5001) NSTK
C        reaction indices
         WRITE (LINSK, ERR=5001) (ISTK(N), N=1,NSTK)
C        Motz-Wise correction indicator flag
         WRITE (LINSK, ERR=5001) (MOTZ(N), N=1,NSTK)
      ENDIF
C
      IF (NBHM .GT. 0) THEN
C        Bohm ion-collision correction reactions
         WRITE (LINSK, ERR=5001) NBHM
C        reaction indices
         WRITE (LINSK, ERR=5001) (IBHM(N), N=1,NBHM)
C        ionic species indices
         WRITE (LINSK, ERR=5001) (KBHM(N), N=1,NBHM)
      ENDIF
C
      IF (NRNU .GT. 0) THEN
C        real stoichiometry reactions
         WRITE (LINSK, ERR=5001) NRNU
C        reaction indices
         WRITE (LINSK, ERR=5001) (IRNU(N), N=1,NRNU)
C        stoichiometric coefficient array for that reaction
         WRITE (LINSK, ERR=5002) ((RNU(L,N), L=1,MAXSPR), N=1,NRNU)
      ENDIF
C
      IF (NORD .GT. 0) THEN
C        Modified species order reactions
         WRITE (LINSK, ERR=5001) NORD, MAXORD
C        reaction indicies
         WRITE (LINSK, ERR=5001) (IORD(N), N = 1, NORD)
C        order dependency species indices
         WRITE (LINSK, ERR=5001) ((KORD(L,N),L=1,MAXORD),N=1,NORD)
C        order dependency values for the species assigned by KORD(L,N)
         WRITE (LINSK, ERR=5002) ((RORD(L,N),L=1,MAXORD),N=1,NORD)
      ENDIF
C
      IF (NEDP .GT. 0) THEN
C        reactions for which ENRGDEP keyword was given
         WRITE (LINSK, ERR=5001) NEDP, NEDPAR
C        reaction indices
         WRITE (LINSK, ERR=5001) (IEDP(N), N = 1, NEDP)
C        species indices
         WRITE (LINSK, ERR=5001) (KEDP(N), N = 1, NEDP)
C        parameters
         WRITE (LINSK, ERR=5002) ((PEDP(L,N),L=1,NEDPAR),N=1,NEDP)
      ENDIF
C
      IF (NYLD .GT. 0) THEN
C        modified yield reactions
         WRITE (LINSK, ERR=5001) NYLD, NYPAR
C        reaction indices
         WRITE (LINSK, ERR=5001) (IYLD(N), N=1,NYLD)
C        ion indices
         WRITE (LINSK, ERR=5001) (KYION(N), N=1,NYLD)
C        yield flags for reaction species
         WRITE (LINSK, ERR=5001) ((KYLD(L,N),L=1,MAXSPR),N=1,NYLD)
C        yield parameters
         WRITE (LINSK, ERR=5002) ((PYLD(L,N),L=1,NYPAR), N=1,NYLD)
C        net phase balance affect after yield-modify
         WRITE (LINSK, ERR=5002) ((YNCF(L,N),L=1,NPHASE),N=1,NYLD)
      ENDIF
C
  500 CONTINUE
C
C     CLEANUP
      WRITE (LOUT, '(///A,/A,A,A)')
     1' NO ERRORS FOUND ON INPUT:',
     2' BINARY Version ',FILVER(1:CKLSCH(FILVER)),
     3' surface linkfile surf.bin written.'
C
      RETURN
C
C     ERROR HANDLING
 5001 CONTINUE
      WRITE (LOUT, '(////A)')
     1   'ERROR: Failure to write binary integer data to surf.bin'
      RETURN
 5002 CONTINUE
      WRITE (LOUT, '(///A)')
     1   'ERROR: Failure to write binary real data to surf.bin'
      RETURN
 5003 CONTINUE
      WRITE (LOUT, '(///A)')
     1   'ERROR: Failure to write binary character data to surf.bin'
      RETURN
      END
C
      SUBROUTINE SKFORM
     1   (LINSK, LOUT, SKMIN, MDIM, KDIM, IDIM, MXPHSE, KERR, MAXSPR,
     2    MAXTP, NCP, NSPAR, NSCOV, NEDPAR, NYPAR, MAXORD, MATNAM,
     3    ENAME, AWT, KNAME, WT, KNCF, KCHRG, NT, KPHSE, KCOV, DEN, TMP,
     4    NFIT, NTR, A, KION, PNAME, KFIRST, KLAST, KKPHAS, PDEN, NSPEC,
     5    NREAC, NU, NUNK, NUSUMK, SPAR, RNCF, ICOV, IKCOV, CPAR, IREV,
     6    RSPAR, ISTK, MOTZ, IBHM, KBHM, IRNU, RNU, IORD, KORD, RORD,
     7    IEDP, KEDP, PEDP, IYLD, KYION, KYLD, YNCF, PYLD, EQFAC)
C
C  START PROLOGUE
C
C  SUBROUTINE SKFORM
C  writes surface mechanism data to formatted surface linkfile.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
      COMMON /SKVERS/ PRVERS, PREC, FILVER
C
C     Integer arrays
      DIMENSION KNCF(MDIM,KDIM), KCHRG(KDIM), NT(KDIM), KPHSE(KDIM),
     1          KCOV(KDIM), KFIRST(MXPHSE), KLAST(MXPHSE),
     2          KKPHAS(MXPHSE), KION(NKION), NSPEC(IDIM),
     3          NREAC(IDIM), NUSUMK(IDIM), NU(MAXSPR, IDIM),
     4          NUNK(MAXSPR, IDIM), IREV(IDIM), ICOV(IDIM),
     5          IKCOV(IDIM), ISTK(IDIM), MOTZ(IDIM), IBHM(IDIM),
     6          KBHM(IDIM), IRNU(IDIM), IORD(IDIM), KORD(MAXORD,IDIM),
     7          IEDP(IDIM), KEDP(IDIM), IYLD(IDIM), KYION(IDIM),
     8          KYLD(MAXSPR,IDIM)
C     Real arrays
      DIMENSION AWT(MDIM), WT(KDIM), DEN(KDIM), TMP(MAXTP,KDIM),
     1          A(NFIT,NTR,KDIM), PDEN(MXPHSE), SPAR(NSPAR,IDIM),
     2          RNCF(MXPHSE,IDIM), CPAR(NSCOV,IDIM),
     3          RSPAR(NSPAR,IDIM), RNU(MAXSPR,IDIM), RORD(MAXORD,IDIM),
     3          PEDP(NEDPAR,IDIM), YNCF(MXPHSE,IDIM), PYLD(NYPAR,IDIM),
     5          EQFAC(IDIM)
C
      CHARACTER*16 MATNAM, ENAME(MDIM), KNAME(KDIM), PNAME(MXPHSE),
     1             PRVERS, PREC, FILVER
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      CHARACTER*16 CFMT, IFMT, LFMT, RFMT
      PARAMETER
     1(CFMT='(8A16)', IFMT='(10I12)', LFMT='(L8)', RFMT='(1P,5E24.16)')
C
C     Linkfile version (should be 1.0)
      WRITE (LINSK, '(A16)', ERR=5003) FILVER
C
C     Skinterp program version
      WRITE (LINSK, '(A16)', ERR=5003) PRVERS
C
C     Precision
      WRITE (LINSK, '(A16)', ERR=5003) PREC
C
C     Whether skinterp successfully completed or not
      WRITE (LINSK, '(L8)', ERR=5001) KERR
C
      IF (KERR) THEN
         WRITE (LOUT, '(//A,/A)')
     1   'WARNING...THERE IS AN ERROR IN THE SURFACE LINKFILE',
     2   '          DUE TO SURFACE MECHANISM INPUT ERRORS.'
         RETURN
      ENDIF
C
C     Work space sizes
      WRITE (LINSK, '(4I12)', ERR=5001) LENISK,LENRSK,LENCSK,MORE
C
C     Constants
      MXTK = 0
      DO 120 K = 1, KKTOT
         MXTK = MAX (MXTK, NT(K))
  120 CONTINUE
      WRITE (LINSK, IFMT, ERR=5001) MAXSPR, MXTK, NCP, NSPAR, NSCOV,
     1                         NEDPAR, NYPAR, MAXORD
C
C     Problem specifications
      WRITE (LINSK, IFMT, ERR=5001)
     1     NELEM, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE, NFSURF,
     2     NLSURF, NNSURF, NFBULK, NLBULK, NNBULK, NIISUR, NCOV,
     3     NREV, NSTK, NCON, NBHM, NRNU, NORD, NEDP, NYLD, NKION,
     4     KELECT, MORE
C
C     Elemental, charge, and site balance tolerance parameter
      WRITE (LINSK, RFMT, ERR=5002) SKMIN
C
C     material name for this surface
      WRITE (LINSK, '(A16)', ERR=5003) MATNAM
C
C     element names
      WRITE (LINSK, CFMT, ERR=5003) (ENAME(N), N = 1, NELEM)
C
C     elemental atomic weights
      WRITE (LINSK, RFMT, ERR=5002) (AWT(N), N = 1, NELEM)
C
C     species names
      WRITE (LINSK, CFMT, ERR=5003) (KNAME(K), K = 1, KKTOT)
C
C     species molecular weights
      WRITE (LINSK, RFMT, ERR=5002) (WT(K), K = 1, KKTOT)
C
C     species elemental compositions
      WRITE (LINSK, IFMT, ERR=5001) ((KNCF(N,K), N=1,NELEM), K=1,KKTOT)
C
C     species electronic charges
      WRITE (LINSK, IFMT, ERR=5001) (KCHRG(K), K = 1, KKTOT)
C
C     the number of temperatures regions used to fit thermodynamics
C     for each species
      WRITE (LINSK, IFMT, ERR=5001) (NT(K), K = 1, KKTOT)
C
C     an integer representing physical state for each species
      WRITE (LINSK, IFMT, ERR=5001) (KPHSE(K), K = 1, KKTOT)
C
C     site coverage of each species
      WRITE (LINSK, IFMT, ERR=5001) (KCOV(K), K = 1, KKTOT)
C
C     densities for each species
      WRITE (LINSK, RFMT, ERR=5002) (DEN(K), K = 1, KKTOT)
C
C     array of temperatures used to specify boundaries of
C     temperature regions used in the fits to thermodynamic functions
      WRITE (LINSK, RFMT, ERR=5002) ((TMP(L,K),L=1,MXTK),K=1,KKTOT)
C
C     fits to thermodynamic functions for each species;
C     NTR=MAXTP-1, NFIT=NCP+2
      WRITE (LINSK, RFMT, ERR=5002)
     1      (((A(N,L,K), N=1,NFIT), L=1,MXTK-1), K=1,KKTOT)
C
C     integer species indices for ionic species
      IF (NKION .GT. 0) THEN
         WRITE (LINSK, IFMT, ERR=5001) NKION
         WRITE (LINSK, IFMT, ERR=5001) (KION(N), N = 1, NKION)
      ENDIF
C
C     phase names
      WRITE (LINSK, CFMT, ERR=5003) (PNAME(N), N=1,NPHASE)
C
C     phase first species indices
      WRITE (LINSK, IFMT, ERR=5001) (KFIRST(N), N = 1, NPHASE)
C     phase last species indices
      WRITE (LINSK, IFMT, ERR=5001) (KLAST(N), N = 1, NPHASE)
C     phase species counts
      WRITE (LINSK, IFMT, ERR=5001) (KKPHAS(N), N = 1, NPHASE)
C
C     phase densities
      WRITE (LINSK, RFMT, ERR=5002) (PDEN(N), N = 1, NPHASE)
C
      IF (NIISUR .LE. 0) THEN
         WRITE (LOUT, '(/A,/A)')
     1   ' WARNING...NO SURFACE REACTIONS FOUND, ',
     2   ' LINKFILE HAS NO REACTION INFORMATION ON IT.'
         GO TO 500
      ENDIF
C
C     SURFACE REACTION OUTPUT
C     count of species defined as reactants and products
      WRITE (LINSK, IFMT, ERR=5001) (NSPEC(I), I=1,NIISUR)
C     count of reactants defined for each reaction, also,
C     a negative number flags an irreversible reaction
      WRITE (LINSK, IFMT, ERR=5001) (NREAC(I), I=1,NIISUR)
C     stoichiometric information for each reaction
      WRITE (LINSK, IFMT, ERR=5001)
     1      ((NU(N,I), NUNK(N,I), N=1,MAXSPR), I=1,NIISUR)
C     stoichiometric coefficients totals
      WRITE (LINSK, IFMT, ERR=5001) (NUSUMK(I), I=1,NIISUR)
C     forward Arrhenius coefficients for each reaction
      WRITE (LINSK, RFMT, ERR=5002) ((SPAR(N,I),N=1,NSPAR),I=1,NIISUR)
C     Multiplicative equilibrium factor
      WRITE (LINSK, RFMT, ERR=5002) (EQFAC(I), I=1,NIISUR)
C     Phase balance for each reaction
      WRITE (LINSK, RFMT, ERR=5002) ((RNCF(N,I),N=1,NPHASE),I=1,NIISUR)
C
C     SPECIAL REACTIONS:
C
      IF (NCOV .GT. 0) THEN
C        coverage reactions
         WRITE (LINSK, IFMT, ERR=5001) NCOV, NSCOV
C        reaction indices
         WRITE (LINSK, IFMT, ERR=5001) (ICOV(N), N = 1, NCOV)
C        coverage species indices
         WRITE (LINSK, IFMT, ERR=5001) (IKCOV(N), N = 1, NCOV)
C        coverage expression parameters
         WRITE (LINSK, RFMT, ERR=5002) ((CPAR(L,N),L=1,NSCOV),N=1,NCOV)
      ENDIF
C
      IF (NREV .GT. 0) THEN
C        reactions with explicit reverse rate constants
         WRITE (LINSK, IFMT, ERR=5001) NREV, NSPAR
C        reaction indices
         WRITE (LINSK, IFMT, ERR=5001) (IREV(N), N=1,NREV)
C        Arrhenius reverse rate constant coefficients
         WRITE (LINSK, RFMT, ERR=5002)
     1   ((RSPAR(L,N),L=1,NSPAR),N=1,NREV)
      ENDIF
C
      IF (NSTK .GT. 0) THEN
C        reactions for which parameters act as sticking coefficients
         WRITE (LINSK, IFMT, ERR=5001) NSTK
C        reaction indices
         WRITE (LINSK, IFMT, ERR=5001) (ISTK(N), N=1,NSTK)
C        Motz-Wise correction indicator flag
         WRITE (LINSK, IFMT, ERR=5001) (MOTZ(N), N=1,NSTK)
      ENDIF
C
      IF (NBHM .GT. 0) THEN
C        Bohm ion-collision correction reactions
         WRITE (LINSK, IFMT, ERR=5001) NBHM
C        reaction indices
         WRITE (LINSK, IFMT, ERR=5001) (IBHM(N), N=1,NBHM)
C        ionic species indices
         WRITE (LINSK, IFMT, ERR=5001) (KBHM(N), N=1,NBHM)
      ENDIF
C
      IF (NRNU .GT. 0) THEN
C        real stoichiometry reactions
         WRITE (LINSK, IFMT, ERR=5001) NRNU
C        reaction indices
         WRITE (LINSK, IFMT, ERR=5001) (IRNU(N), N=1,NRNU)
C        stoichiometric coefficient array for that reaction
         WRITE (LINSK, RFMT, ERR=5002)
     1   ((RNU(L,N), L=1,MAXSPR), N=1,NRNU)
      ENDIF
C
      IF (NORD .GT. 0) THEN
C        Modified species order reactions
         WRITE (LINSK, IFMT, ERR=5001) NORD, MAXORD
C        reaction indices
         WRITE (LINSK, IFMT, ERR=5001) (IORD(N), N = 1, NORD)
C        order dependency species indices
         WRITE (LINSK, IFMT, ERR=5001)
     1   ((KORD(L,N),L=1,MAXORD),N=1,NORD)
C        order dependency values for the species assigned by KORD(L,N)
         WRITE (LINSK, RFMT, ERR=5002)
     1   ((RORD(L,N),L=1,MAXORD),N=1,NORD)
      ENDIF
C
      IF (NEDP .GT. 0) THEN
C        reactions for which ENRGDEP keyword was given
         WRITE (LINSK, IFMT, ERR=5001) NEDP, NEDPAR
C        reaction indices
         WRITE (LINSK, IFMT, ERR=5001) (IEDP(N), N = 1, NEDP)
C        species indices
         WRITE (LINSK, IFMT, ERR=5001) (KEDP(N), N = 1, NEDP)
C        parameters
         WRITE (LINSK, RFMT, ERR=5002)
     1   ((PEDP(L,N),L=1,NEDPAR),N=1,NEDP)
      ENDIF
C
      IF (NYLD .GT. 0) THEN
C        modified yield reactions
         WRITE (LINSK, IFMT, ERR=5001) NYLD, NYPAR
C        reaction indices
         WRITE (LINSK, IFMT, ERR=5001) (IYLD(N), N=1,NYLD)
C        ion indices
         WRITE (LINSK, IFMT, ERR=5001) (KYION(N), N=1,NYLD)
C        yield flags for reaction species
         WRITE (LINSK, IFMT, ERR=5001)
     1   ((KYLD(L,N),L=1,MAXSPR),N=1,NYLD)
C        yield parameters
         WRITE (LINSK, RFMT, ERR=5002)
     1   ((PYLD(L,N),L=1,NYPAR), N=1,NYLD)
C        net phase balance affect after yield-modify
         WRITE (LINSK, RFMT, ERR=5002)
     1   ((YNCF(L,N),L=1,NPHASE),N=1,NYLD)
      ENDIF
C
  500 CONTINUE
C
C     CLEANUP
      WRITE (LOUT, '(///A,/A,A,A)')
     1' NO ERRORS FOUND ON INPUT:',
     2' ASCII Version ',FILVER(1:CKLSCH(FILVER)),
     3' surface linkfile surf.asc written.'
      RETURN
C
C     ERROR HANDLING
 5001 CONTINUE
      WRITE (LOUT, '(////A)')
     1   'ERROR: Failure to write ascii integer data to surf.asc'
      RETURN
 5002 CONTINUE
      WRITE (LOUT, '(///A)')
     1   'ERROR: Failure to write ascii real data to surf.asc'
      RETURN
 5003 CONTINUE
      WRITE (LOUT, '(///A)')
     1   'ERROR: Failure to write ascii character data to surf.asc'
      RETURN
C
C     end of SUBROUTINE CKFORM
      END
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKAUXL
     1         (LOUT, KDIM, MDIM, IDIM, MXPHSE, NSPAR, MAXSPR, NSCOV,
     2          NEDPAR, NYPAR, MAXORD, SUB, NSUB, NSPEC, NREAC, NU,
     3          NUNK, KNAME, PNAME, KFIRST, KLAST, KKPHAS, KCHRG, IDUP,
     4          IREV, RSPAR, ICOV, IKCOV, CPAR, ISTK, LIMOTZ, MOTZ,
     5          IRNU, RNU, IORD, KORD, RORD, IBHM, IEDP, PEDP, IYLD,
     6          PYLD, AXUNIT, EXUNIT, SKMIN, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKAUXL
C  processes character strings in lines which follow a reaction to
C  to specify additional options;
C  strings are in the form of a "KEYword" (3 characters are usually
C  sufficient, but more may be used to better describe the option),
C  and may be followed by slash(/)-delimited arguments:
C
C  DUPlicate               - reaction is allowed to be duplicated.
C  STIck                   - rate parameters are        CALL SKNSTK.
C                            sticking coefficients,
C  BOHm                    - this is a Bohm reaction,   CALL SKNBHM.
C  MWO[FF],MWO[N]          - set Motz-Wise correction   CALL SKMOTZ.
C  REV/val(n), n=1,3/      - explicit reverse rate,     CALL SKNREV.
C  UNIt/string/            - rate parameter units,      CALL SKUNIT.
C  COV/KNAME val(n),n=1,3) - coverage reaction,         CALL SKNCOV.
C  ENR/val(n), n=1,3/      - energy-dependence,         CALL SKNEDP.
C  FORd/KNAME val1/        - species forward order,     CALL SKNORD.
C  RORd/KNAME val1/        - species reverse order,     CALL SKNORD.
C  YIEld/val(n),n=1,4/     - species yield-modify,      CALL SKNYLD.
C                            *required "+#" in reaction
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  KDIM     - Integer scalar, maximum number of species.
C  MDIM     - Integer scalar, maximum number of elements.
C  IDIM     - Integer scalar, maximum number of reactions.
C  MXPHSE   - Integer scalar, maximum number of phases.
C  NSPAR    - Integer scalar, required number of Arrhenius parameters.
C  MAXSPR   - Integer scalar, maximum number of reaction species.
C  NSCOV    - Integer scalar, required number of coverage parameters.
C  NEDPAR   - Integer scalar, required number of T-dep. parameters.
C  NYPAR    - Integer scalar, required number of yield parameters.
C  SUB(*)   - Character string array of option keyword and parameters.
C  NSUB     - Integer scalar, total count of character strings.
C  NSPEC(*) - Integer array, reaction's count of reactants+products.
C  NREAC(*) - Integer array, reaction's count of reactants only.
C  NU(*)    - Integer array, reaction's stoichiometric coefficients.
C  NUNK(*)  - Integer array, reaction's species indices.
C  KNAME(*) - Character string array, species names.
C  PNAME(*) - Character string array, phase names.
C  KFIRST(*)- Integer array, starting species indices for phases.
C  KLAST(*) - Integer array, ending species indices for phases.
C  KKPHAS(*)- Integer array, total count of species for the phases.
C  KCHRG(*)  - Integer array, species electronic charges.
C  NEDPAR   - Integer scalar, required number of energy-dep parameters.
C  IDUP(*)  - Integer array, a flag to allow duplicate reactions.
C             -1, reaction is allowed to be duplicated,
C              0, not.
C  IREV(*)  - Integer array, explicit reverse parameter rxn indices.
C  RSPAR(*,*)-Real matrix, reverse Arrhenius rate parameters.
C  ICOV(*)  - Integer array, coverage reaction indices.
C  IKCOV(*) - Integer array, coverage reaction coverage species indices.
C  CPAR(*,*)- Real matrix, coverage reaction parameters.
C  ISTK(*)  - Integer array, sticking reaction indices.
C  LIMOTZ   - Logical, Motz-Wise correction flag for sticking reactions.
C  MOTZ(*)  - Integer array, Motz-Wise 0/1 correction flags.
C  IRNU(*)  - Integer array, real stoichiometry reaction indices.
C  RNU(*,*) - Real matrix, real stoichiometric coefficients.
C  IORD(*)  - Integer array, changed-order reaction indices.
C  KORD(*,*)- Integer matrix, changed-order species indices.
C             < 0, change forward order for species KORD
C             > 0, change reverse order for species KORD
C  RORD(*,*)- Real matrix, order values for change-order species.
C  KERR     - Logical error flag.
C  IBHM(*) - Integer array, Bohm-correction reaction indices.
C  IEDP(*)  - Integer array, ion-energy dependent reaction indices.
C  PEDP(*,*)-Real matrix, coefficients for energy-dependent reactions.
C  IYLD(*)  - Integer array, yield-modified reaction indices.
C  PYLD(*,*)- Real matrix, yield-modify parameters.
C  AXUNIT   - Character string, description of pre-exponential A units.
C  EXUNIT   - Character string, description of activation energy E units
C  SKMIN    - Real scalar, error tolerance.
C
C  END PROLOGUE
C
C----------------------------------------------------------------------C
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
C
C     Integer arrays
      DIMENSION NU(MAXSPR), NUNK(MAXSPR), IDUP(IDIM), IREV(IDIM),
     1          ICOV(IDIM), IKCOV(IDIM), ISTK(IDIM), MOTZ(IDIM),
     2          KFIRST(MXPHSE), KLAST(MXPHSE), KKPHAS(MXPHSE),
     4          IORD(IDIM), KORD(MAXORD,IDIM),
     5          IRNU(IDIM), IBHM(IDIM), KCHRG(KDIM), IEDP(IDIM),
     6          IYLD(IDIM)
C     Real arrays
      DIMENSION RSPAR(NSPAR,IDIM), CPAR(NSCOV,IDIM),
     1          RORD(MAXORD,IDIM), RNU(MAXSPR,IDIM),
     2          PEDP(NEDPAR,IDIM), PYLD(NYPAR,IDIM)
      CHARACTER*(*) SUB(NSUB), AXUNIT, EXUNIT
      CHARACTER*16 KNAME(KDIM), PNAME(MXPHSE)
      CHARACTER RSTR*80, IUNITS*80, KEY*3, CKCHUP*3
      LOGICAL KERR, LDUM1, LDUM2, LSTK, LREV, LCOV, LBHM, LEDP,
     1        LORD, LRNU, LIMOTZ, LYLD
      EXTERNAL CKCHUP
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C     NSUB sets of auxiliary information
      DO 100 N = 1, NSUB
         ILEN = CKLSCH(SUB(N))
         IF (ILEN .LE. 0) GO TO 100
C        check status of this reaction so far
         LSTK = NIISUR .EQ. ISTK(MAX(NSTK,1))
         LREV = NIISUR .EQ. IREV(MAX(NREV,1))
         LCOV = NIISUR .EQ. ICOV(MAX(NCOV,1))
         LBHM = NIISUR .EQ. IBHM(MAX(NBHM,1))
         LEDP = NIISUR .EQ. IEDP(MAX(NEDP,1))
         LORD = NIISUR .EQ. IORD(MAX(NORD,1))
         LYLD = NIISUR .EQ. IYLD(MAX(NYLD,1))
         LRNU = NIISUR .EQ. IRNU(MAX(NRNU,1))
C
         KEY = ' '
         KEY = CKCHUP(SUB(N),3)
C
C        DUP, STI, BOH, MWO[N],[F] do not need additional parameters
C
         IF (KEY .EQ. 'DUP') THEN
            IDUP(NIISUR) = -1
            WRITE (LOUT, '(6X,A)') 'Declared duplicate reaction...'
            GO TO 100
C
         ELSEIF (KEY .EQ. 'STI') THEN
C
C           Arrhenius parameters are sticking coefficients;
            CALL SKNSTK (LOUT, IDIM, LSTK, NSTK, ISTK, LIMOTZ, MOTZ,
     1                   NIISUR, KERR)
            GO TO 100
C
         ELSEIF (KEY .EQ. 'BOH') THEN
C           use Bohm velocity correction for collisions of
C           positive ions with a surface
            CALL SKNBHM (LOUT, IDIM, LBHM, NBHM, IBHM, KKGAS, KELECT,
     1                   NSPEC, NIISUR, KERR)
            GO TO 100
C
         ELSEIF (KEY .EQ. 'MWO') THEN
C           Motz-Wise correction to sticking coefficients
            RSTR = ' '
            RSTR = CKCHUP(SUB(N)(4:4), 1)
            CALL SKMOTZ (LOUT, IDIM, LSTK, NSTK, SUB(N), RSTR, MOTZ)
            GO TO 100
         ENDIF
C
C        other options require /-delimited parameters
         CALL CKDLIM (SUB(N), '/', I1, I2)
         IF (I1 .GE. I2) THEN
            KERR = .TRUE.
            WRITE (LOUT, *)
     1      ' Error searching for /-delimited parameters...',
     2        SUB(N)(1:ILEN)
            GO TO 100
         ENDIF
C
C        contents of /-delimited data
         RSTR = ' '
         RSTR = SUB(N)(I1+1:I2-1)
C
         IF (KEY .EQ. 'COV') THEN
C           coverage-dependence of the rate constant was specified:
C
C           COV / species parameters / is the simplest construct,
C           KEY='COV' and RSTR=SUB(N)(I1+1:I2-1)='species parameters',
C           BUT...it might look like two unrelated constructs,
C           COV / species / phase / parameters /, in which case
C           KEY='COV' and RSTR='species',
C           with next SUB(N+1) = phase / parameters /,
C           so need to find out if only one substring in RSTR
C
            ILAST = CKLSCH(RSTR)
            IF (INDEX(RSTR(1:ILAST), ' ') .GT. 0) THEN
C
C             there is more than one substring in RSTR, so assume
C             it does contain the species and the parameters
            ELSE
C             it is likely that SUB(N+1) contains phase / parameters /,
C             so make RSTR look like species/phase/ parameters
              IF (N .EQ. NSUB) THEN
C                can't do any more
              ELSE
                 CALL CKDLIM (SUB(N+1), '/', I1, I2)
                 IF (I1.GT.0 .AND. I2.GT.I1) THEN
                    RSTR(ILAST+1:ILAST+1) = '/'
                    RSTR(ILAST+2:) = SUB(N+1)
C                   RSTR now looks like "species/phase/ parameters";
C                   blank out next substring
                    SUB(N+1) = ' '
                 ENDIF
               ENDIF
            ENDIF
C
            CALL SKNCOV (LOUT, KDIM, IDIM, MXPHSE, RSTR, NIISUR,  LCOV,
     1                   NCOV, ICOV, KNAME, KKTOT, PNAME, KKPHAS,
     2                   KFIRST(NFSURF), KLAST(NLSURF), NSCOV, CPAR,
     3                   IKCOV, KERR)
            GO TO 100
         ENDIF
C
         IF (KEY .EQ. 'REV') THEN
C
C           reverse arrhenius parameters given
            CALL SKNREV (LOUT, IDIM, NIISUR, RSTR, LREV, NSPEC, NSPAR,
     1                   NREV, IREV, RSPAR, KERR)
C
         ELSEIF (KEY .EQ. 'ENR') THEN
C
C           ion energy-dependence of the rate coefficient
            CALL SKNEDP (LOUT, IDIM, RSTR, NIISUR, LEDP, NEDP, IEDP,
     1                   NEDPAR, PEDP, KERR)
C
         ELSEIF (KEY .EQ. 'UNI') THEN
C
C           units conversion for reaction parameters is required:
            IF (CKCHUP(RSTR,3) .EQ. 'NON') THEN
C              site NONConservation keyword only useful on REACTION line
               WRITE (LOUT, '(6X,A,A)')
     1         'Warning...non-conservation of sites ignored here...',
     2         SUB(N)(1:ILEN)
            ELSE
               IUNITS = ' '
               CALL SKUNIT (RSTR, 1, LOUT, AXUNIT, EXUNIT, IUNITS,
     1                      LDUM1, LDUM2)
               ILAST = CKLSCH(IUNITS)
               WRITE (LOUT, '(6X,A,A)')
     1         'Units for this reaction are...',IUNITS(1:ILAST)
            ENDIF
C
         ELSEIF (KEY .EQ. 'FOR' .OR. KEY .EQ. 'ROR') THEN
C
C           non-standard reaction order declared
            CALL SKNORD (LOUT, KDIM, IDIM, KEY, RSTR, NIISUR, KKTOT,
     1                   KNAME, PNAME, NPHASE, KKPHAS,
     2                   NORD, MAXORD, IORD, KORD, RORD, MAXSPR,
     2                   NUNK, NU, NREAC, NSPEC, LORD, LRNU, NRNU, RNU,
     3                   KERR)
C
         ELSEIF (KEY .EQ. 'YIE') THEN
C
            CALL SKNYLD (LOUT, IDIM, RSTR, NIISUR, NSPEC, KKGAS,
     1                   LYLD, NYLD, IYLD, NYPAR, PYLD, KERR)
C
         ELSE
            WRITE (LOUT, '(6X,A,A)')
     1      'Unrecognized option ignored...', SUB(N)(1:ILEN)
         ENDIF
  100 CONTINUE
C
C     end of SUBROUTINE SKAUXL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKBAL (LOUT, MDIM, KDIM, IDIM, MXPHSE, MAXSPR, NU,
     1                  NUNK, KCHRG, KCOV, NONCON, IRNU, RNU, IYLD,
     2                  KYLD, KFIRST, KLAST, KNCF, RNCF, YNCF, PNAME,
     3                  SKMIN, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKBAL
C  checks a reaction to ensure that its reactants' mass and electronic
C  charge equals its product's mass and electronic charge.
C
C  Arguments:
C  LOUT      - Integer scalar, formatted output file unit number.
C  MDIM      - Integer scalar, maximum number of elements.
C  KDIM      - Integer scalar, maximum number of species.
C  IDIM      - Integer scalar, maximum number of reactions.
C  MAXPHSE   - Integer scalar, maximum number of phases.
C  MAXSP     - Integer scalar, reactions maximum number of species.
C  NU(*)     - Integer array, reaction species stoichiometric coeff'nts.
C  NUNK(*)   - Integer array, reaction species indices.
C  KCHRG(*)  - Integer array, species electronic charges;
C  KCOV(*)   - Integer array, species site coverage.
C  NONCON    - Logical, site conservation flag.
C  IRNU(*)   - Integer array, real stoichiometry reaction indices.
C  RNU(*,*)  - Real matrix, real stoichiometric coefficients.
C  IYLD(*)   - Integer array, yield-modify reaction indices.
C  KYLD(*,*) - Integer array, reaction species yield-modify flags.
C  KFIRST(*) - Integer array, phases first species indices.
C  KLAST(*)  - Integer array, phases last species indices.
C  KNCF(*,*) - Integer matrix, elemental composition of species.
C  RNCF(*,*) - Real matrix, reactions phase balances.
C  YNCF(*,*) - Real matrix, yield-modify reactions phase balances.
C  PNAME(*)  - Character string array, phase names.
C  SKMIN     - Real scalar, an error tolerance for balancing.
C  KERR      - Logical error flag.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
C     Integer arrays
      DIMENSION NU(MAXSPR), NUNK(MAXSPR), IRNU(IDIM), IYLD(IDIM),
     1          KYLD(MAXSPR,IDIM), KNCF(MDIM,KDIM), KCHRG(KDIM),
     2          KCOV(KDIM), KFIRST(MXPHSE), KLAST(MXPHSE)
C     Real arrays
      DIMENSION RNU(MAXSPR,IDIM), YNCF(MXPHSE,IDIM), RNCF(MXPHSE)
      CHARACTER*16 PNAME(MXPHSE)
      LOGICAL KERR, NONCON, LRNU, LYLD, IERR, YERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      LRNU = NIISUR .EQ. IRNU(MAX(NRNU,1))
      LYLD = NIISUR .EQ. IYLD(MAX(NYLD,1))
C
      IERR = .FALSE.
      YERR = .FALSE.
      DO 60 M = 1, NELEM
C        non-# coefficient element/mass balance of reaction
         SUMM = 0.0
C        #-modified coefficient element/mass balance of reaction
         YSUMM = 0.0
         DO 50 N = 1, MAXSPR
C           index number for this species
            NK = NUNK(N)
            IF (NK .NE. 0) THEN
C              stoichiometric coefficient for this species
               IF (LRNU) THEN
                  COEF = RNU(N, NRNU)
               ELSE
                  COEF = NU(N)
               ENDIF
C
               YCOEF = 0.0
               IF (LYLD) THEN
C                 #-modified stoichiometric coefficients exist
                  IF (KYLD(N, NYLD) .GT. 0) THEN
C                    #-modified stoichiometric coefficient 
                     YCOEF = COEF
                     COEF = 0.0
                  ENDIF
               ENDIF
C              sum of this element over non-# coefficients
               SUMM = SUMM + COEF * KNCF(M, NK)
C              sum of this element over #-modified coefficients
               YSUMM = YSUMM + YCOEF * KNCF(M, NK)
            ENDIF
   50    CONTINUE
         IF (ABS(SUMM) .GT. SKMIN) IERR = .TRUE.
         IF (ABS(YSUMM).GT. SKMIN) YERR = .TRUE.
   60 CONTINUE
      IF (IERR) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...reaction does not balance in elements/mass...'
      ENDIF
      IF (YERR) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...reaction does not balance in #-modified elements/',
     2   'mass...'
      ENDIF
C
C     non-# coefficient charge balance of reaction
      SUMCH = 0.0
C     #-modified coefficient charge balance of reaction
      YSUMCH = 0.0
      DO 70 N = 1, MAXSPR
C        index number for this species
         NK = NUNK(N)
         IF (NK .GT. 0) THEN
C           non-# stoichiometric coefficient for this species
            IF (LRNU) THEN
               COEF = RNU(N,NRNU)
            ELSE
               COEF = NU(N)
            ENDIF
C
            YCOEF = 0.0
            IF (LYLD) THEN
C              #-modified stoichiometric coefficients exist
               IF (KYLD(N,NYLD) .GT. 0) THEN
C                 #-modified stoichiometric coefficient for this species
                  YCOEF = COEF
                  COEF = 0.0
               ENDIF
            ENDIF
C           sum of charge over non-# coefficients
            SUMCH = SUMCH + COEF * KCHRG(NK)
C           sum of charge over #-modified coefficients
            YSUMCH = YSUMCH + YCOEF * KCHRG(NK)
         ENDIF
   70 CONTINUE   
      IF (ABS(SUMCH) .GT. SKMIN) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1   'Error...reaction does not balance in electronic charge...'
      ENDIF
      IF (ABS(YSUMCH) .GT. SKMIN) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...reaction does not balance in #-modified electronic ',
     2   'charge...'
      ENDIF
C
C     site balance
      DO 80 J = 1, NPHASE
C        balance sites for a phase
         DO 75 N = 1, MAXSPR
C           index number for this species
            NK = NUNK(N)
            IF (NK.GE.KFIRST(J) .AND. NK.LE.KLAST(J)) THEN
C              non-# stoichiometric coefficient for this species
               IF (LRNU) THEN
                  COEF = RNU(N,NRNU)
               ELSE
                  COEF = NU(N)
               ENDIF
C
               YCOEF = 0.0
               IF (LYLD) THEN
C                 #-modified stoichiometric coefficients exist
                  IF (KYLD(N,NYLD) .GT. 0) THEN
C                    #-modified stoich. coefficient for this species
                     YCOEF = COEF
C                    site sum for this phase, #-modified coefficients
                     YNCF(J,NYLD) = YNCF(J,NYLD) + YCOEF * KCOV(NK)
                     COEF = 0.0
                  ENDIF
               ENDIF
C              sum of sites for this phase over non-# coefficients
               RNCF(J) = RNCF(J) + COEF * KCOV(NK)
            ENDIF
   75    CONTINUE
   80 CONTINUE
C
      DO 90 J = 2, NNSURF + 1
C        site phase conservation test
         IF (ABS(RNCF(J)) .GT. SKMIN .OR.
     1       (LYLD .AND. ABS(YNCF(J,MAX(NYLD,1))) .GT. SKMIN)) THEN
C           a site does not balance
            IF (NONCON) THEN
               NCON = NCON + 1
               RETURN
            ELSE
               KERR = .TRUE.
               ILEN = CKLSCH(PNAME(J))
               WRITE (LOUT, '(6X,A,A)')
     1         'Error...reaction does not conserve site ',
     2         PNAME(J)(1:ILEN)
            ENDIF
         ENDIF
   90 CONTINUE
C
C     end of SUBROUTINE SKBAL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKBULK (SUB, NSUB, NDIM, STRAY, RAY, LRAY, NN, KSTART,
     1                   KERR, LOUT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKBULK
C  processes an array of character-strings, storing character
C  information in array STRAY, real-value conversions into real array
C  RAY, and setting logical flags in logical array LRAY.
C
C  Arguments:
C  SUB(*)     - Character-string array.
C  NSUB       - Integer scalar, count of character strings.
C  NDIM       - Integer scalar, size of STRAY, RAY, and LRAY.
C  STRAY(*)   - Character-string array.
C  RAY(*)     - Real array.
C  LRAY(*)    - Logical array.
C  NN         - Integer scalar, latest count of STRAY, RAY, and LRAY.
C  KSTART     - Integer scalar.
C  KERR       - Logical error flag.
C  LOUT       - Integer scalar, formatted output file unit number.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION RAY(NDIM)
      CHARACTER SUB(NSUB)*(*), STRAY(NDIM)*(*), ISTR*80, CKCHUP*4
      LOGICAL KERR, LRAY(NDIM)
      INTEGER CKLSCH
      EXTERNAL CKCHUP, CKLSCH
C
      ILEN = LEN(STRAY(1))
C
      DO 200 N = 1, NSUB
         IF (CKCHUP(SUB(N), 3) .EQ. 'END') RETURN
         ILAST = CKLSCH(SUB(N))
C
         ISTR = ' '
         I1 = INDEX(SUB(N),'/')
         IF (I1 .EQ. 1) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...misplaced value...',SUB(N)(1:ILAST)
         ELSE
            IF (I1 .LE. 0) THEN
               ISTR = SUB(N)
            ELSE
               ISTR = SUB(N)(1:I1-1)
            ENDIF
            CALL CKNCMP (ISTR, STRAY, NN, KNUM1, NF)
C
            IF (NN .GT. KSTART) THEN
               CALL CKNCMP (ISTR, STRAY(KSTART), NN-KSTART, KNUM2, NF)
            ELSE
               KNUM2 = 0
            ENDIF
C
            IF (KNUM1 .GT. 0) THEN
               WRITE (LOUT, '(6X,A,A,A)')
     1         'Error...bulk species name duplicates ',
     2         'gas or site species...',SUB(N)(1:ILAST)
               KERR = .TRUE.
            ELSEIF (KNUM2 .GT. 0) THEN
               WRITE (LOUT, '(6X,A,A)')
     1         'Warning...duplicate bulk species name ignored...',
     2         SUB(N)(1:ILAST)
            ELSE
               IF (NN .LT. NDIM) THEN
                  IF (ISTR(ILEN+1:) .NE. ' ') THEN
                     WRITE (LOUT, '(6X,A,A)')
     1               'Error...bulk species name too long...',
     2               SUB(N)(1:ILAST)
                     KERR = .TRUE.
                  ELSE
                     NN = NN + 1
                     STRAY(NN) = ' '
                     STRAY(NN) = ISTR(1:ILEN)
                     IF (I1 .GT. 0) THEN
                        I2 = I1 + INDEX(SUB(N)(I1+1:),'/')
                        ISTR = ' '
                        ISTR = SUB(N)(I1+1:I2-1)
                        CALL CKPARR (ISTR,1,1,RAY(NN), NVAL, IER, LOUT)
                        KERR = KERR .OR. (IER.NE.0)
                        LRAY(NN) = .TRUE.
                     ENDIF
                  ENDIF
               ELSE
                  WRITE (LOUT, '(6X,A,A)')
     1            'Error...species array size too small for...',
     2            SUB(N)(1:ILAST)
                  KERR = .TRUE.
               ENDIF
            ENDIF
         ENDIF
  200 CONTINUE
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKDUP (LOUT, IDIM, NIISUR, MAXSPR, NSPEC, NREAC, NU,
     1                  NUNK, IDUP, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKDUP
C  checks all components of reaction NIISUR against the (NIISUR-1)
C  reactions for illegal duplication.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  IDIM     - Integer scalar, maximum number of reactions.
C  NIISUR   - Integer scalar, total reaction count, and index of
C             current reaction.
C  MAXSPR   - Integer scalar, maximum reaction species count.
C  NSPEC(*) - Integer array, reactions counts of products+reactants.
C  NREAC(*) - Integer array, reactions counts of reactants onlyu.
C  NU(*,*)  - Integer matrix, reactions stoichiometric coefficients.
C  NUNK(*,*)- Integer matrix, reactions participant species indices.
C  IDUP(*)  - Logical array, flag for reactions duplication.
C  KERR     - Logical error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NSPEC(IDIM), NREAC(IDIM), NU(MAXSPR,IDIM),
     1          NUNK(MAXSPR,IDIM), IDUP(IDIM)
      LOGICAL KERR
C
      ISAME = 0
      KHALF = MAXSPR/2
      II = NIISUR
C
C     II is the current reaction
C     NRI is II's reactant species count
C     NPI is II's product species count
      NRI = NREAC(II)
      NPI = ABS(NSPEC(II)) - NREAC(II)
C
C     check II against the previous reactions
      DO 500 J = 1, II-1
C
C        NRJ is J's reactant species count
C        NPJ is J's product species count
         NRJ = NREAC(J)
         NPJ = ABS(NSPEC(J)) - NREAC(J)
C
C        II is not like J if a species count is different
         IF (NRI.NE.NRJ .OR. NPI.NE.NPJ) GO TO 100
C
C        does J have same reactants as II
         DO 20 N = 1, NRI
            KI = NUNK(N,II)
            IF (KI .NE. 0) THEN
               CALL SKINUM (KI, NUNK(1,J), NRI, KJ)
C              II is not like J if KI reactant not in J,
C              or if it's coefficient is different
               IF (KJ .LE. 0) GO TO 100
               IF (NU(KJ,J) .NE. NU(N,II)) GO TO 100
            ENDIF
   20    CONTINUE
C
C        does J have same products as II
         DO 25 N = KHALF+1, KHALF + NPI
            KI = NUNK(N,II)
            IF (KI .NE. 0) THEN
               CALL SKINUM (KI, NUNK(KHALF+1,J), NPI, KJ)
               KJ = KHALF + KJ
C              II is not like J if KI product not in J,
C              or if it's coefficient is different
               IF (KJ .LE. KHALF) GO TO 100
               IF (NU(KJ,J) .NE. NU(N,II)) GO TO 100
            ENDIF
   25    CONTINUE
C
C        same products, reactants, coefficients, pres dependency, and
C        third-body relationship
C
         ISAME = J
         GO TO 600
C
  100    CONTINUE
C
C        II not same as J in forward direction; check reverse,
C        if II has same number of reactants as J has products,
C        and same number of products as J has reactants
C
         IF (NPI.NE.NRJ .OR. NPJ.NE.NRI) GO TO 500
C
C        check I reactants against J products
         DO 30 N = 1, NRI
            KI = NUNK(N,II)
            IF (KI .NE. 0) THEN
               CALL SKINUM (KI, NUNK(KHALF+1,J), NPJ, KJ)
               KJ = KHALF + KJ
C              II is not like J if KI reactant not in J products,
C              or if it's coefficient is different
               IF (KJ .LE. KHALF) GO TO 500
               IF (NU(N,II) .NE. -NU(KJ,J)) GO TO 500
            ENDIF
   30    CONTINUE
C
C        check I products against J reactants
         DO 35 N = KHALF+1, KHALF + NPI
            KI = NUNK(N,II)
            IF (KI .NE. 0) THEN
               CALL SKINUM (KI, NUNK(1,J), NRJ, KJ)
C              II is not like J if KI product not in J reactants,
C              or if it's coefficient is different
               IF (KJ .LE. 0) GO TO 500
               IF (-NU(N,II) .NE. NU(KJ,J)) GO TO 500
            ENDIF
   35    CONTINUE
C
C        J products same as II reactants, and
C        J reactants same as II products
         IF (NSPEC(J).LT.0 .AND. NSPEC(II).LT.0) THEN
C           only OK if J and II are both irreversible
         ELSE
            ISAME = J
            GO TO 600
         ENDIF
C
  500 CONTINUE
C
  600 CONTINUE
C
      IF (ISAME .EQ. 0) RETURN
C
C     Reaction #ISAME is a duplicate
      IF (IDUP(ISAME).NE.0 .AND. IDUP(II).NE.0) THEN
C        Reaction #ISAME is a legal duplicate
         IDUP(ISAME) = ABS(IDUP(ISAME))
         IDUP(II) = ABS(IDUP(II))
         RETURN
      ENDIF
C
      KERR = .TRUE.
      WRITE (LOUT, 1050) ISAME
 1050 FORMAT (6X,'Error...undeclared duplicate to reaction number ',I3)
C
C     end of SUBROUTINE SKDUP
      RETURN
      END
C
      SUBROUTINE SKIEQ (KDIM, IDIM, MXPHSE, MAXSPR, KFIRST, KLAST,
     1                  KCOV, NU, NUNK, PDEN, IRNU, RNU, SKMIN, EQFAC)
C
C  START PROLOGUE
C
C  SUBROUTINE SKIEQ
C  sets an equilibrium constant scalar for reactions.
C
C  Arguments:
C  KDIM     - Integer scalar, maximum number of species.
C  IDIM     - Integer scalar, maximum number of reactions.
C  MXPHSE   - Integer scalar, maximum number of phases.
C  MAXSPR   - Integer scalar, maximum number of reaction species.
C  KFIRST(*)- Integer array, first species indices of phases.
C  KLAST(*) - Integer array, last species indices of phases.
C  KCOV(*)  - Integer array, species site coverages.
C  NU(*)    - Integer array, reactions stoichiometric coefficients.
C  NUNK(*)  - Integer array, reactions participant species.
C  PDEN(*)  - Real array, phase densities.
C  IRNU(*)  - Integer array, real stoichiometry reaction indices.
C  RNU(*,*) - Real matrix, real stoichiometric coefficients.
C  SKMIN    - Real scalar, small difference for scaling, balancing.
C  EQFAC    - Real scalar, equilibrium constant scalar.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
C
C     Integer arrays
      DIMENSION KFIRST(MXPHSE), KLAST(MXPHSE), KCOV(KDIM), NU(MAXSPR),
     1          NUNK(MAXSPR), IRNU(IDIM)
C     Real arrays
      DIMENSION PDEN(MXPHSE), RNU(MAXSPR,IDIM)
      LOGICAL LRNU
C
      LRNU = NIISUR .EQ. IRNU(MAX(NRNU,1))
      IF (LRNU) THEN
         DO 20 N = NFSURF, NLSURF
            RNUSUM = 0.0
            DO 10 K = KFIRST(N), KLAST(N)
               RKCOV = KCOV(K)
               DO 5 M = 1, MAXSPR
                  IF (NUNK(M).EQ.K .AND. ABS(RNU(M,NRNU)).GT.SKMIN)
     1            THEN
                     EQFAC = EQFAC * RKCOV**(-RNU(M,NRNU))
                     RNUSUM = RNUSUM + RNU(M,NRNU)
                  ENDIF
    5          CONTINUE
   10       CONTINUE
            EQFAC = EQFAC * PDEN(N)**RNUSUM
   20    CONTINUE
      ELSE
         DO 40 N = NFSURF, NLSURF
            NUSUM = 0
            DO 30 K = KFIRST(N), KLAST(N)
               RKCOV = KCOV(K)
               DO 25 M = 1, MAXSPR
                  IF (NUNK(M).EQ.K .AND. NU(M).NE.0) THEN
                     EQFAC = EQFAC * RKCOV**(-NU(M))
                     NUSUM = NUSUM + NU(M)
                  ENDIF
   25          CONTINUE
   30       CONTINUE
            EQFAC = EQFAC * PDEN(N)**NUSUM
   40    CONTINUE
      ENDIF
C
C     end of SUBROUTINE SKIEQ
      RETURN
      END
C
      SUBROUTINE SKINUM (I, IARRAY, NI, IND)
C
C  START PROLOGUE
C
C  SUBROUTINE SKINUM
C  searches integer array IARRAY of length NI for the
C  first occurrence of integer I, and returns its index IND.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION IARRAY(NI)
C
      IND = 0
      DO 50 N = 1, NI
         IF (I .EQ. IARRAY(N)) THEN
            IND = N
            RETURN
         ENDIF
   50 CONTINUE
C
C     end of SUBROUTINE SKINUM
      RETURN
      END
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKISUB (LINE, SUB, NSUB)
C
C  START PROLOGUE
C
C  SUBROUTINE SKISUB
C  splits a longer character string into substrings, using blanks as
C  the delimited; any portion of input preceeded by a '!' is ignored.
C
C  Arguments:
C  LINE     - Character string.
C  SUB(*)   - Character-string array.
C  NSUB     - Integer scalar, count of substrings.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) SUB(*), LINE
      INTEGER CKFRCH, CKLSCH, CKSLEN
      EXTERNAL CKFRCH, CKLSCH, CKSLEN
C
C     initialize substring count
      NSUB = 0
C
C     if input line is zero length, just return
      IF (CKSLEN(LINE) .LE. 0) RETURN
C
C     index of last non-blank character in string
      ILEN = CKLSCH(LINE)
C
C     index of first non-blank character in string
      NSTART = CKFRCH(LINE)
C
   10 CONTINUE
C
      ISTART = NSTART
      NSUB = NSUB + 1
      SUB(NSUB) = ' '
C
      ILAST = INDEX(LINE(ISTART:),' ') - 1
      IF (ILAST .GT. 0) THEN
         ILAST = ISTART + ILAST - 1
      ELSE
         ILAST = ILEN
      ENDIF
      SUB(NSUB) = LINE(ISTART:ILAST)
C
C     end of line
      IF (ILAST .EQ. ILEN) GO TO 50
C
      NSTART = ILAST + CKFRCH(LINE(ILAST+1:))
C
C     Does SUB have any slashes?
C
      I1 = INDEX(SUB(NSUB),'/')
      IF (I1 .LE. 0) THEN
C        no slash found
         IF (LINE(NSTART:NSTART) .NE. '/') GO TO 10
         NEND = NSTART + INDEX(LINE(NSTART+1:),'/')
         IND = INDEX(SUB(NSUB),' ')
         SUB(NSUB)(IND:) = LINE(NSTART:NEND)
C        end of line, so quit
         IF (NEND .EQ. ILEN) GO TO 50
         NSTART = NEND + CKFRCH(LINE(NEND+1:))
         GO TO 10
      ENDIF
C
C     one slash found; does SUB have 2 slashes?
C
      I2 = INDEX(SUB(NSUB)(I1+1:),'/')
      IF (I2 .GT. 0) GO TO 10
C
C     two slashes
      NEND = NSTART + INDEX(LINE(NSTART+1:),'/')
      IND = INDEX(SUB(NSUB),' ') + 1
      SUB(NSUB)(IND:) = LINE(NSTART:NEND)
C
C     end of line
      IF (NEND .LT. ILEN) THEN
         NSTART = NEND + CKFRCH(LINE(NEND+1:))
         GO TO 10
      ENDIF
C
  50  CONTINUE
C
C     end of SUBROUTINE SKISUB
      RETURN
      END
C                                                                      C
      SUBROUTINE SKKEL (MDIM, MELECT, NELEM, KNCF, KEL)
C
C  START PROLOGUE
C
C  SUBROUTINE SKKEL
C  flags an electron species by examining its elemental composition.
C
C  Arguments:
C  MELECT   - Integer scalar, element index of the electron.
C  NELEM    - Integer scalar, total element count.
C  KNCF(*)  - Integer array, elemental composition of this species.
C  KEL      - Integer scalar, 1 if species contains only an electron.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KNCF(MDIM)
C
      KEL = 0
C     return if does not contain one and only one electron
      IF (KNCF(MELECT) .NE. 1) RETURN
      DO 10 N = 1, NELEM
C        return if has element other than the electron
         IF (N.NE.MELECT .AND. KNCF(N).GT.0) RETURN
   10 CONTINUE
C     has only one element, and it is one electron
      KEL = 1
C
C     end of SUBROUTINE SKKEL
      RETURN
      END
C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKKEY
     1   (LIN, LTHRM, LOUT, MDIM, KDIM, IDIM, MXPHSE, MAXSPR, MAXTP,
     2    NFIT, NSPAR, NSCOV, NEDPAR, NYPAR, MAXORD, ENAME, AWT, KNAME,
     3    ITHRM, WT, KNCF, NONCON, KCOV, KPHSE, KCHRG, NT, NTR, TMP, A,
     4    PNAME, PDEN, LPDEN, KFIRST, KLAST, KKPHAS, DEN, LKDEN, SPAR,
     5    NSPEC, NREAC, NU, NUNK, NUSUMK, RNCF, IDUP, IREV, RSPAR,
     6    ICOV, IKCOV, CPAR, ISTK, KERR, IRNU, RNU, SKMIN, IORD, KORD,
     7    RORD, IBHM, KBHM, LMOTZ, MOTZ, IEDP, KEDP, PEDP, IYLD, KYION,
     8    KYLD, PYLD, YNCF, MATNAM, EQFAC)
C
C  START PROLOGUE
C
C  SUBROUTINE SKKEY
C  is the main parsing routine for character-string input of a surface
C  mechanism.  Depending on input flags "SITE", "BULK", "THERMO",
C  "REACTION", and specialized auxiliary reaction line options,
C  character strings are passed to other subroutines for processing
C  of mechanism information.
C
C  Arguments:
C  LIN      - Integer scalar, formatted input file unit number
C             for the surface mechanism.
C  LTHRM    - Integer scalar, formatted input file unit number
C             for thermodynamic data.
C  LOUT     - Integer scalar, formatted output file unit number.
C  MDIM     - Integer scalar, maximum number of elements.
C  KDIM     - Integer scalar, maximum number of species.
C  IDIM     - Integer scalar, maximum number of reactions.
C  MXPHSE   - Integer scalar, maximum number of phases.
C  MAXSPR   - Integer scalar, maximum number of reaction species.
C  MAXTP    - Integer scalar, maximum number of fit temperatures.
C  NFIT     - Integer scalar, maximum number of fit coefficients.
C  NSPAR    - Integer scalar, required number of Arrhenius parameters.
C  NSCOV    - Integer scalar, required number of coverage parameters.
C  NEDPAR   - Integer scalar, required number of ion-energy-dependence
C             parameters.
C  NYPAR    - Integer scalar, required number of yield-modify params.
C  MAXORD   - Integer scalar, maximum number of reaction change-orders.
C  ENAME(*) - Character-string array, element names.
C  AWT(*)   - Real array, element atomic weights.
C  KNAME(*) - Character-string array, species names.
C  ITHRM(*) - Logical array, species thermodynamic flags.
C  WT(*)    - Real array, species molecular weights.
C  KNCF(*,*)- Integer matrix, species elemental composition.
C  NONCON   - Logical, site conservation flag.
C  KCOV(*)  - Integer array, species site coverage.
C  KPHSE(*) - Integer array, species physical state.
C  KCHRG(*) - Integer array, species electronic charge.
C  NT(*)    - Integer array, species number of fit temperatures.
C  NTR      - Integer scalar, number of fit temperature ranges.
C  TMP(*,*) - Real matrix, species fit temperatures.
C  A(*,*,*) - Real matrix, species thermodynamic polynomial coeff'nts.
C  PNAME(*) - Character-string array, phase names.
C  PDEN(*)  - Real array, phase densities.
C  LPDEN(*) - Logical array, phase density flag.
C  KFIRST(*)- Integer array, phases first species index.
C  KLAST(*) - Integer array, phases last species index.
C  KKPHAS(*)- Integer array, phases total species count.
C  DEN(*)   - Real array, species densities.
C  LKDEN(*) - Logical array, species density flag.
C  SPAR(*,*)- Real matrix, reactions Arrhenius parameters.
C  NSPEC(*) - Integer array, count of reactions reactants+products.
C  NREAC(*) - Integer array, count of reactions reactants only.
C  NU(*,*)  - Integer matrix, reactions stoichiometric coefficients.
C  NUNK(*,*)- Integer matrix, reactions species indices.
C  NUSUMK(*)- Integer array, reactions sums of stoichiometric coeff'nts.
C  RNCF(*,*)- Real matrix, reactions net phase changes.
C  IDUP(*)  - Integer array, reaction duplication flag.
C  IREV(*)  - Integer array, reverse parameter reaction indices.
C  RSPAR(*,*)-Real matrix, reverse reaction parameters.
C  ICOV(*)  - Integer array, coverage reaction indices.
C  IKCOV(*) - Integer array, coverage reaction species indices.
C  CPAR(*,*)- Real matrix, coverage reaction parameters.
C  ISTK(*)  - Integer array, stiking reaction indices.
C  KERR     - Logical, error flag.
C  IRNU(*)  - Integer array, real stoichiometry reaction indices.
C  RNU(*,*) - Real matrix, real stoichiometric coefficients.
C  SKMIN    - Real scalar, minimum value for content and balancing.
C  IORD(*)  - Integer array, change-order reaction indices.
C  KORD(*,*)- Integer matrix, change-order reaction species indices.
C  RORD(*,*)- Real matrix, change-order values.
C  IBHM(*)  - Integer array, Bohm-correction reaction indices.
C  KBHM(*)  - Integer array, Bohm-correction reaction species indices.
C  LMOTZ    - Logical, sticking reaction Motz-Wise correction flag.
C  MOTZ(*)  - Integer array, sticking reaction Motz-Wise option.
C  IEDP(*)  - Integer array, ion-energy reaction indices.
C  KEDP(*)  - Integer array, ion-energy reaction species indices.
C  PEDP(*,*)- Real matrix, ion-energy reaction parameters.
C  IYLD(*)  - Integer array, yield-modify reaction indices.
C  KYION(*) - Integer array, yield-modify reaction ion indices.
C  KYLD(*,*)- Integer matrix, reaction species yield-modify flags.
C  PYLD(*,*)- Real matrix, yield-modify reaction parameters.
C  YNCF(*,*)- Real matrix, yield-modify reaction phase balances.
C  MATNAM   - Character string, material/surface name.
C  EQFAC(*) - Real array, reaction equilibrium rate scalar.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C     MAXIMUM NUMBER OF CONTINUATION LINES ALLOWED
      PARAMETER (MXCONT=10)
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
C
C     Integer arrays
      DIMENSION AWT(KDIM), WT(KDIM), KNCF(MDIM,KDIM), KCOV(KDIM),
     1          KPHSE(KDIM), KCHRG(KDIM), NT(KDIM), KKPHAS(MXPHSE),
     3          KFIRST(MXPHSE), KLAST(MXPHSE),
     4          NSPEC(IDIM), NU(MAXSPR,IDIM), NUNK(MAXSPR,IDIM),
     6          NREAC(IDIM), NUSUMK(IDIM), IDUP(IDIM), IREV(IDIM),
     7          ICOV(IDIM), IKCOV(IDIM), ISTK(IDIM), MOTZ(IDIM),
     8          IRNU(IDIM), IORD(IDIM), KORD(MAXORD,IDIM), IBHM(IDIM),
     *          KBHM(IDIM), IEDP(IDIM), KEDP(IDIM), IYLD(IDIM),
     1          KYLD(MAXSPR,IDIM), KYION(IDIM)
C     Real arrays
      DIMENSION A(NFIT,NTR,KDIM), PDEN(MXPHSE), DEN(KDIM),
     1          TMP(MAXTP,KDIM), SPAR(NSPAR,IDIM), RNCF(MXPHSE,IDIM),
     1          RSPAR(NSPAR,IDIM), CPAR(NSCOV,IDIM), RNU(MAXSPR,IDIM),
     2          RORD(MAXORD,IDIM), PEDP(NEDPAR,IDIM),
     3          PYLD(NYPAR,IDIM), YNCF(MXPHSE,IDIM), EQFAC(IDIM)
      CHARACTER*16 KNAME(KDIM), ENAME(MDIM), PNAME(MXPHSE), MATNAM
      CHARACTER ISTR*80, SUB(80)*80, IUNITS*80, LINE(MXCONT)*80,
     1          ILINE*80, CKCHUP*4, SKEY*4, IST(20)*2, EUNITS*4,
     2          AUNITS*4, AXUNIT*4, EXUNIT*4, IAUXL*3
      LOGICAL KERR, IERR, THERMO, ITHRM(KDIM), NONCON, 
     1        LPDEN(MXPHSE), LKDEN(MXPHSE), LTHERM, LMOTZ, LIMOTZ,
     2        LMATE, SKFILE
      INTEGER CKLSCH, CKSLEN
      EXTERNAL CKLSCH, CKSLEN, CKCHUP, SKFILE
      DATA IST/'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ','10',
     1         '11','12','13','14','15','16','17','18','19','20'/
C
C     if THERMO is true, we need to read information from a separate
C     thermo database file
      THERMO = .TRUE.
C
C     LTHERM will be false unless we encounter the 'THERM' keyword
C     in the input file
      LTHERM = .FALSE.
C
      IF (NMAT .EQ. 1) THEN
C
C        if we are processing the first material, then we will need
C        to print out the gas-phase information;
C        make a column header with element names:
         WRITE (LOUT, 1500)
         WRITE (LOUT, 1600) (ENAME(M)(1:2),M=1,NELEM)
         WRITE (LOUT, 1500)
C
C        list the gas species names and elemental compositions
         WRITE (LOUT, '(/1X,A)') 'Gas phase species:'
         DO 75 K = 1, KKGAS
            WRITE (LOUT, 1650)
     1      K, KNAME(K), WT(K), (KNCF(M,K),M=1,NELEM)
   75    CONTINUE
      ENDIF
C
C     the first phase is always the gas
      NPHASE = 1
C
C     the first species index in the gas is always 1
      KFIRST(NPHASE) = 1
C
C     the index of the last species in the phase
      KLAST (NPHASE) = KKGAS
C
C     the total number of species in this phase
      KKPHAS(NPHASE) = KKGAS
C
C     the name of the phase
      PNAME (NPHASE) = 'GAS'
C
C     total number of species processed in the mechanism so far
      KKTOT = KKGAS
C
C     no more materials follow this one
      MORE = 0
C
C     logical flag says we already have read a line containing the
C     keyword "MATE"
      LMATE = .FALSE.
C
C     the name of the current material
      MATNAM = 'MATERIAL'//IST(NMAT)
C
C     flag indicating what task we are currently doing;
C     ITASK of 0 indicates that there are no tasks currently pending
      ITASK = 0
C
C     initialize reaction lines
      NLINES = 0
      DO 90 N = 1, MXCONT
         LINE(N) = ' '
   90 CONTINUE
C
  100 CONTINUE
C
C     read a line of input from the input file
      ILINE = ' '
      READ (LIN, '(A)', END=9999) ILINE
      CALL CKDTAB (ILINE)
C
  105 CONTINUE
C     get the length of the information read (non-blank)
      ILEN = CKSLEN(ILINE)
C
C     if there was nothing on the line, go read the next line
      IF (ILEN .LE. 0) GO TO 100
C
      CALL SKISUB (ILINE(1:ILEN), SUB, NSUB)
C
C-----look for 'keyword' input
C
C     convert the first four characters of the first substring to
C     upper case and put the result into the variable SKEY
      SKEY = CKCHUP(SUB(1), 4)
C
      IF (SKEY.EQ.'MATE') THEN
C        we have encountered a new material line
C
         IF (LMATE) THEN
C           the LMATE flag has been set already; we are free to finish
C           processing the current material and return from the
C           subroutine
C
C           we will backspace the file, so that the line will be read
C           again the first time through when this subroutine is called
C           again (called again because another material follows)
            BACKSPACE (LIN)
C
C           set the flag to indicate that more materials follow
            MORE = 1
C
C           finished processing the mechanism and then exit the
C           subroutine
            GO TO 9999
         ENDIF
C
C        the logical flag has not been set, so this is the first time
C        that we have read the line with the material name; set the
C        flag, so that the next time a material line is read, we know
C        that we have finished with the current material
         LMATE = .TRUE.
C
C        get the name of the material, if one was supplied
         IF (NSUB .GT. 1) MATNAM = SUB(2)
C
C        print out header information with the material name,
C        and columns of the element names
         WRITE (LOUT, *)
         WRITE (LOUT, 1500)
         WRITE (LOUT, '(1X,A,A)') 'MATERIAL:',MATNAM
         WRITE (LOUT, 1600) (ENAME(M)(1:2),M=1,NELEM)
         WRITE (LOUT,1500)
C
C        get another line of input, and continue processing keywords
C        for the current material
         GO TO 100
      ENDIF
C
      IF (SKEY.EQ.'SITE' .OR. SKEY.EQ.'BULK') THEN
C
C--------new phase
C
         IF (NPHASE+1 .GT. MXPHSE) THEN
C           we have not procided enough space for another phase;
C           issue a warning and set the error flag
            WRITE (LOUT, '(A)')' Error...too many phases...'
            KERR = .TRUE.
C
         ELSE
C
C           increment the number of phases
            NPHASE = NPHASE + 1
C           find a '/'-delimited string
            CALL CKDLIM (SUB(1), '/', I1, I2)
            IF (I2 .GT. I1) THEN
C              the phase name appears between the two '/' delimiters
               PNAME(NPHASE) = SUB(1)(I1+1:I2-1)
            ENDIF
C
C           no errors encountered so far
            IERR = .FALSE.
C
C           the number of the first species in the phase will be zero
C           if we are just starting to process the phase; if so,
C           set it to the index of the next species to be read
            IF (KFIRST(NPHASE) .EQ. 0) KFIRST(NPHASE) = KKTOT+1
C
C           save the current total number of species
            KOLD = KKTOT
C
            IF (SKEY .EQ. 'SITE') THEN
C
C              a surface site phase is to be processed;
C              set an appropriate task flag
               ITASK = 1
C
C              increment the number of surface phases found
               NNSURF = NNSURF + 1
C
C              if a name for the phase was not supplied, then
C              provide a default
               IF (PNAME(NPHASE) .EQ. ' ')
     1             PNAME(NPHASE)='SITE'//IST(NNSURF)
C
C              make sure that this name is not duplicated by a phase
C              name already read in
               CALL CKNCMP (PNAME(NPHASE), PNAME, NPHASE-1, I, NF)
               IF (I .GT. 0) THEN
C                 a duplicate name was found; tell the user and set
C                 an error flag
                  ILS = CKLSCH(PNAME(NPHASE))
                  ISTR = ' '
                  ISTR = PNAME(NPHASE)(1:ILS)
                  WRITE (LOUT, '(A)')
     1            ' Warning...duplicate phase name'//ISTR(1:ILS)
                  KERR = .TRUE.
               ENDIF
C
C              if we have just found our first surface phase, set the
C              pointer to the first surface phase
               IF (NNSURF .EQ. 1) NFSURF = NPHASE
C
C              the current number of the surface phase is the last
C              surface phase (at least until another one is encountered)
               NLSURF = NPHASE
C
               IF (NNBULK .GT. 0) THEN
C                 we have already encountered a bulk phase, so we
C                 should not be reading another surface phase; issue
C                 a warning and set the error flag
                  WRITE (LOUT, '(A)')
     1            ' Error...sites must precede bulk phases...',
     2            PNAME(NPHASE)
                  KERR = .TRUE.
               ENDIF
C
               IF (NIISUR .GT. 0) THEN
C                 we have already encountered a surface reaction, so we
C                 should not be reading another surface phase; issue
C                 a warning and set the error flag
                  WRITE (LOUT, '(A)')
     1            ' Error...site must precede reactions...',
     2            PNAME(NPHASE)
                  KERR = .TRUE.
               ENDIF
C
               IF (LTHERM) THEN
C                 we have already encountered thermo information, so we
C                 should not be reading another surface phase; issue
C                 a warning and set the error flag
                  WRITE (LOUT, '(A)')
     1            ' Error...site must precede THERMO data...',
     2            PNAME(NPHASE)
                  KERR = .TRUE.
               ENDIF
C
C              read information about the surface phase; this may be the
C              site density or the names of surface species; this
C              routine adds species names to KNAME and increment KKTOT
C              for any species found
               CALL SKSURF (SUB(2), NSUB-1, KDIM, KNAME, KCOV, KKTOT,
     1                      KKGAS, KFIRST(NPHASE), LPDEN(NPHASE),
     2                      PDEN(NPHASE), IERR, LOUT)
C
            ELSE
C
C              processing a bulk phase; set an appropriate task flag
               ITASK = 2
C
C              increment the number of bulk phases found
               NNBULK = NNBULK + 1
C
C              if no name for this phase was provided, supply a
C              default name
               IF (PNAME(NPHASE) .EQ. ' ')
     1             PNAME(NPHASE)='BULK'//IST(NNBULK)
C
C              make sure that the name of this phase doesn't
C              duplicate one already encountered
               CALL CKNCMP (PNAME(NPHASE), PNAME, NPHASE-1, I, NF)
C
               IF (I .GT. 0) THEN
C                 the name for this phase was already used;
C                 issue a warning and set an error flag
                  ILS = CKLSCH(PNAME(NPHASE))
                  ISTR = ' '
                  ISTR = PNAME(NPHASE)(1:ILS)
                  WRITE (LOUT, '(A)')
     1            ' Warning...duplicate phase name',ISTR(1:ILS)
                  KERR = .TRUE.
               ENDIF
C
C              if this is the first bulk phase that we have encountered,
C              set the pointer to the first bulk phase
               IF (NNBULK .EQ. 1) NFBULK = NPHASE
C
C              the current bulk phase is also the last one (so far)
               NLBULK = NPHASE
C
               IF (NIISUR .GT. 0) THEN
C                 we have already read some surface reactions, so we
C                 should not be encountering any more bulk phases;
C                 issue a warning and set the error flag
                  WRITE (LOUT, '(A)')
     1            ' Error...bulk phase must precede reactions...',
     2            PNAME(NPHASE)
                  KERR = .TRUE.
               ENDIF
C
               IF (LTHERM) THEN
C                 we have already read some thermo data, so we
C                 should not be encountering any more bulk phses;
C                 issue a warning and set the error flag
                  WRITE (LOUT, '(A)')
     1            ' Error...bulk phase must precede THERMO data...',
     2            PNAME(NPHASE)
                  KERR = .TRUE.
               ENDIF
C
C              read information about the bulk phase; this may be the
C              bulk densities or the names of bulk species; this
C              routine will add species names to KNAME and increment
C              KKTOT for any species found
               CALL SKBULK (SUB(2), NSUB-1, KDIM, KNAME, DEN, LKDEN,
     1                      KKTOT, KFIRST(NPHASE), IERR, LOUT)
            ENDIF
C
C           find the number of new species read by either SKSURF or
C           SKBULK
            KNEW   = KKTOT - KOLD
C
C           increment the number of species in this phase by the
C           number of species just found
            KKPHAS(NPHASE) = KKPHAS(NPHASE) + KNEW
C
C           increment the total number of surface or bulk species
C           by the number just found
            IF (ITASK .EQ. 1) THEN
C              we are working on surface species
               KKSURF = KKSURF + KNEW
            ELSE
C              we are working on bulk species
               KKBULK = KKBULK + KNEW
            ENDIF
C
C           update the pointer to the last species in the current phase
            KLAST(NPHASE) = KKTOT
C           update the error flag
            KERR = KERR.OR.IERR
         ENDIF
         GO TO 100
C
      ELSEIF (SKEY .EQ. 'REAC') THEN
C
C--------reactions
C
C        set the task flag to an appropriate value for continuing
C        to process reactions
         ITASK = 3
C
C        User may specify rate parameter units on same line as REACTION,
C        or 'NONCON' options to allow non-conservation of sites,
C        or 'MWOFF', 'MWON' options (Motz-Wise correction off/on) for
C        sticking coefficients;
C        a character-string units legend is created, and
C        printed after all reactions.
C        User may later use the auxiliary reaction keyword UNITS/.../
C        to set other rate parameter units for the particular reaction,
C        or to toggle MWOFF/MWON for a particular STICKing reaction.
C
         IUNITS = ' '
         EUNITS = ' '
         AUNITS = ' '
         CALL SKUNIT (SUB(2), NSUB-1, LOUT, AUNITS, EUNITS, IUNITS,
     1                NONCON, LMOTZ)
         IF (.NOT. LMOTZ) THEN
C           Motz-wise correction is turned off
            DO 111 I = NIISUR+1, IDIM
               MOTZ(I) = 0
  111       CONTINUE
         ELSE
            DO 112 I = NIISUR+1, IDIM
               MOTZ(I) = 1
  112       CONTINUE
         ENDIF
         GO TO 100
C
      ELSEIF (SKEY .EQ. 'THER') THEN
C
C--------thermodynamic data
C
C        this flag says that we have encountered thermo data;
C        if we should encounter a SITE or BULK keyword after this,
C        it would be an error
         LTHERM = .TRUE.
C
C        set an appropriate task flag meaning that nothing is pending
         ITASK = 0
C
         IF (KKTOT .EQ. KKGAS) THEN
C           do not need thermo data if no additional species, so
C           just read and skip through user's thermo data
            THERMO = .FALSE.
            GO TO 100
         ENDIF
C
C        if not THERMO ALL, need LTHRM at least for temperature ranges,
C        else use LIN for temperature ranges; in either case, use
C        LIN for thermo data
C
C        process user's thermodynamic data
         IF (NSUB.GT.1 .AND. CKCHUP(SUB(2),3).EQ.'ALL') THERMO=.FALSE.
C
C        if not THERMO, don't use LTHRM at all, not even for
C        getting temperature ranges;
C        if THERMO, use LTHRM to get temperature ranges, then use
C        LIN to get thermodynamic data
C
         IF (THERMO) THEN
            IF (.NOT. SKFILE(LTHRM,'FORM')) THEN
               WRITE (LOUT, 334)
               KERR = .TRUE.
               RETURN
            ENDIF
         ENDIF
         CALL SKTHRM (LIN, LTHRM, THERMO, MDIM, KDIM, ENAME, AWT,
     1                KNAME, KNCF, KPHSE, KCHRG, WT, MAXTP, NT, NTR,
     2                TMP, NFIT, A, ITHRM, KERR, LOUT, ILINE)
C
C        use last ILINE read for finding next ITASK
         GO TO 105
C
      ELSEIF (INDEX(SKEY, 'END') .GT. 0) THEN
C        this task flag indicates that no tasks are currently pending
         ITASK = 0
         GO TO 100
      ENDIF
C
C     to get here, no new keyword; process other input phase based on
C     an existing ITASK
C
  555 CONTINUE
C
      IF (ITASK .EQ. 3) THEN
C        we are processing reaction data; is it an auxiliary line?
C
         IND = 0
C        does the line contain any reaction auxiliary keywords?
         DO 400 N = 1, NSUB
            IAUXL = CKCHUP(SUB(N), 3)
            IF (IAUXL.EQ.'COV' .OR. IAUXL.EQ.'DUP' .OR.
     1          IAUXL.EQ.'REV' .OR. IAUXL.EQ.'STI' .OR.
     2          IAUXL.EQ.'FOR' .OR. IAUXL.EQ.'ROR' .OR.
     3          IAUXL.EQ.'BOH' .OR. IAUXL.EQ.'ENR' .OR.
     4          IAUXL.EQ.'UNI' .OR. IAUXL.EQ.'YIE') IND=MAX(IND,N)
  400    CONTINUE
C
         IF (IND .GT. 0) THEN
C           this is auxiliary reaction data;
            IF (NIISUR .LE. 0) THEN
               KERR = .TRUE.
               WRITE (LOUT, '(6X,A)')
     1         'Error...no reaction to modify...'
            ELSE
C              extract the auxiliary information
               CALL SKAUXL
     1         (LOUT, KDIM, MDIM, IDIM, MXPHSE, NSPAR, MAXSPR, NSCOV,
     2          NEDPAR, NYPAR, MAXORD, SUB, NSUB, NSPEC(NIISUR),
     2          NREAC(NIISUR), NU(1,NIISUR), NUNK(1,NIISUR), KNAME,
     2          PNAME, KFIRST, KLAST, KKPHAS, KCHRG, IDUP, IREV,
     3          RSPAR, ICOV, IKCOV, CPAR, ISTK, LIMOTZ, MOTZ, IRNU,
     4          RNU, IORD, KORD, RORD, IBHM, IEDP, PEDP, IYLD, PYLD,
     5          AXUNIT, EXUNIT, SKMIN, KERR)
               GO TO 100
            ENDIF
         ENDIF
C
C        does this reaction have a continuation sign '&'
         IS = INDEX(ILINE(1:ILEN), '&')
         IF (IS .GT. 0) THEN
C           is there any information on the line?
            ILEN = IS - 1
            IF (ILEN.LE.0 .OR. ILINE(1:ILEN).EQ.' ') GO TO 100
         ENDIF
C
C        array of reaction lines
         IF (NLINES+1 .GT. MXCONT) THEN
C           reached maximum number of reaction lines; issue a warning
C           message and set the error flag
            KERR = .TRUE.
            WRITE (LOUT, *)
     1      ' Error...more than 10 continuation lines...'
C
         ELSE
C           increment number of lines and add the input line
            NLINES = NLINES + 1
            LINE(NLINES) = ILINE(1:ILEN)
         ENDIF
C
C        if continuations expected, get another line
         IF (IS .GT. 0) GO TO 100
C
C        have all lines for this reaction
         IF (NIISUR .LE. 0) THEN
C
C           this is the first surface reaction
C           assume end of site and bulk data; need thermo data if not
C           THERMO ALL and more than KKGAS species
C
            IF (KKTOT .EQ. KKGAS) THERMO = .FALSE.
            IF (THERMO) THEN
C              not THERMO ALL, so use LTHRM both for temperature
C              ranges and for thermo data
               IF (.NOT. SKFILE(LTHRM,'FORM')) THEN
                  WRITE (LOUT, 334)
                  KERR = .TRUE.
                  RETURN
               ENDIF
               CALL SKTHRM (LTHRM, LTHRM, THERMO, MDIM, KDIM, ENAME,
     1                      AWT, KNAME, KNCF, KPHSE, KCHRG, WT, MAXTP,
     2                      NT, NTR, TMP, NFIT, A, ITHRM, KERR, LOUT,
     3                      ILINE)
C              done with thermodynamic data
               THERMO = .FALSE.
            ENDIF
C
C           in any THERMO case, check and print species
            CALL SKPRNT (LOUT, KDIM, MXPHSE, MAXTP, KNAME, ITHRM, TMP,
     1                   NT, TMID, WT, PDEN, LPDEN, DEN, LKDEN, MDIM,
     2                   KNCF, KCOV, KFIRST, KLAST, KKPHAS, PNAME, KERR)
            WRITE (LOUT, *)
            WRITE (LOUT, 1500)
            WRITE (LOUT, 1800)
C
         ELSE
C
C           check previous reaction for balance and duplication, and
C           convert units of rate constant coefficients if necessary;
C           we are doing this checking for the previous reaction and not
C           the one just read because we can't be be sure that we have
C           finished with the one just read as there may be more
C           auxiliary information
C
            CALL SPREAC
     1      (LOUT, MDIM, KDIM, IDIM, MXPHSE, MAXSPR, NSPAR, MAXORD,
     2       NEDPAR, NYPAR, NSCOV, PNAME, NSPEC, NREAC, SPAR, RSPAR,
     3       AXUNIT, EXUNIT, NUNK, NU, RNCF, KCHRG, ISTK, ICOV, CPAR,
     4       KCOV, IREV, KNCF, NONCON, IDUP, KFIRST, KLAST, IRNU, RNU,
     6       IEDP, KEDP, PEDP, IYLD, KYLD, PYLD, YNCF, SKMIN, IBHM,
     7       KBHM, IORD, KORD, RORD, PDEN, EQFAC, KNAME, KERR)
         ENDIF
C
C        this is a reaction with one to 10 lines
         IF (NIISUR .GE. IDIM) THEN
C           we have read more reactions than the interpreter can
C           hold; issue a warning message and set an error flag
            WRITE (LOUT, *) ' Reaction dimension too small,',
     1                      ' IDIM must be increased...'
            KERR = .TRUE.
            GO TO 100
         ENDIF
C
C        increment counter for the total number of surface reactions
         NIISUR = NIISUR + 1
C
C        find the reactants, products, stoichiometric coefficients and
C        rate coefficients for the current reaction;
C        if real stoichiometric coefficients, set NRNU, IRNU, and RNU;
C        if #-modified stoichiometry, set NYLD, IYLD, and KYLD;
C        balancing and units conversion will
C        be done after any auxiliary reaction lines are processed.
C
         CALL SKREAC
     1        (LOUT, KDIM, IDIM, MXPHSE, MAXSPR, NSPAR, MXCONT,
     2         LINE, NLINES, KNAME, PNAME, KKPHAS, NSPEC(NIISUR),
     3         NREAC(NIISUR), NUNK(1,NIISUR), NU(1,NIISUR),
     4         NUSUMK(NIISUR), SPAR(1,NIISUR), IRNU, RNU, IYLD,
     5         KYION, KYLD, KCHRG, KERR)
C
C        set auxiliary units declaration string for this reaction
         AXUNIT = AUNITS
         EXUNIT = EUNITS
C        set default Motz-Wise correction for this reaction
         LIMOTZ = LMOTZ
C
C        clear reaction lines
         DO 55 N = 1, NLINES
            LINE(N) = ' '
   55    CONTINUE
         NLINES = 0
         GO TO 100
C
      ELSEIF (ITASK .GT. 0) THEN
C
C        we are reading species data (for a SITE or BULK phase)
C        set the error flag for this task to FALSE
         IERR = .FALSE.
C
C        KOLD holds the number of species found in this phase so far
         KOLD = KKTOT
C
C        if no species have been registered for this phase so far,
C        set the first species in the phase to be the index one
C        greater than have been counted so far
C        (IN POINT OF FACT, COLTRIN CAN'T SEE HOW WE COULD HAVE
C        GOTTEN HERE WITHOUT HAVING SET KFIRST(NPHASE)...)
C
         IF (KFIRST(NPHASE) .EQ. 0) KFIRST(NPHASE) = KOLD + 1
         IF (ITASK .EQ. 1) THEN
C           we are processing surface species; read species names
            CALL SKSURF (SUB, NSUB, KDIM, KNAME, KCOV, KKTOT, KKGAS,
     1                   KFIRST(NPHASE), LPDEN(NPHASE), PDEN(NPHASE),
     2                   IERR, LOUT)
         ELSE
C              we are processing BULK species; read species names
            CALL SKBULK (SUB, NSUB, KDIM, KNAME, DEN, LKDEN, KKTOT,
     1                   KFIRST(NPHASE), IERR, LOUT)
         ENDIF
C
         KERR = KERR.OR.IERR
C        see how many new species we just encountered
         KNEW   = KKTOT - KOLD
C        add that number to the runing total number of species
C        in the phase
         KKPHAS(NPHASE) = KKPHAS(NPHASE) + KNEW
C        increment the total number of surface or bulk species
C        by the number just found
         IF (ITASK .EQ. 1) THEN
            KKSURF = KKSURF + KNEW
         ELSE
            KKBULK = KKBULK + KNEW
         ENDIF
C        update the species number of the last species in the phase
         KLAST(NPHASE) = KKTOT
      ENDIF
C
C     read another line of input
      GO TO 100
C
 9999 CONTINUE
C
C     end of input
C
      IF (NIISUR .GT. 0) THEN
C
C        check final reaction for balance and for duplication
C
         CALL SPREAC
     1      (LOUT, MDIM, KDIM, IDIM, MXPHSE, MAXSPR, NSPAR, MAXORD,
     2       NEDPAR, NYPAR, NSCOV, PNAME, NSPEC, NREAC, SPAR, RSPAR,
     3       AXUNIT, EXUNIT, NUNK, NU, RNCF, KCHRG, ISTK, ICOV, CPAR,
     5       KCOV, IREV, KNCF, NONCON, IDUP, KFIRST, KLAST, IRNU, RNU,
     6       IEDP, KEDP, PEDP, IYLD, KYLD, PYLD, YNCF, SKMIN, IBHM,
     7       KBHM, IORD, KORD, RORD, PDEN, EQFAC, KNAME, KERR)
C
         DO 600 I = 1, NIISUR
C
C           check reactions declared as duplicates
            IF (IDUP(I) .LT. 0) THEN
C              reaction was declared to be a duplicate, but its partner
C              was never found (which would have switched its sign
C              to positive)
               KERR = .TRUE.
               WRITE (LOUT, '(6X,A,I3)')
     1         'Error...no duplicated declared for reaction no.',I
            ENDIF
  600    CONTINUE
C
         WRITE (LOUT, '(/1X,A)') ' NOTE:  '//IUNITS(1:CKLSCH(IUNITS))
         IF (NONCON) WRITE (LOUT, '(10X,A)')
     1               'Site non-conservation specified'
         IF (NSTK .GT. 0) THEN
            IF (.NOT.LMOTZ) THEN
            WRITE (LOUT, '(10X,A,A)')
     1      'Default Motz-Wise correction to sticking coefficients ',
     2      'is turned OFF.'
            ELSE
            WRITE (LOUT, '(10X,A,A)')
     1      'Default Motz-Wise correction to sticking coefficients ',
     2      'is turned ON.'
            ENDIF
         ENDIF
C
C        Done with surface mechanism...
C
         RETURN
      ENDIF
C
      IF (.NOT.THERMO .OR. KKTOT.EQ.KKGAS) RETURN
C
C     the THERMO flag is still true, meaning that we still may need
C     to get information from the thermo data file on unit LTHRM;
C     there were no reactions (thermo data must be complete to do
C     reactions), and user did not use "THERMO ALL"; only need this
C     if user added surface or bulk species.
C
C     read thermodynamic polynomial coefficients from the data
C     file LTHRM
      IF (.NOT. SKFILE(LTHRM,'FORM')) THEN
         WRITE (LOUT, 334)
         KERR = .TRUE.
         RETURN
      ENDIF
      CALL SKTHRM (LTHRM, LTHRM, THERMO, MDIM, KDIM, ENAME, AWT,
     1             KNAME, KNCF, KPHSE, KCHRG, WT, MAXTP, NT, NTR,
     2             TMP, NFIT, A, ITHRM, KERR, LOUT, ILINE)
C
C     print the names of all species in this material and
C     their elemental compositions
      CALL SKPRNT (LOUT, KDIM, MXPHSE, MAXTP, KNAME, ITHRM, TMP,
     1             NT, TMID, WT, PDEN, LPDEN, DEN, LKDEN, MDIM,
     2             KNCF, KCOV, KFIRST, KLAST, KKPHAS, PNAME, KERR)
      WRITE (LOUT, *)
      WRITE (LOUT, 1500)
C
C     Formats
  334 FORMAT (/6X, 'Error accessing thermodynamic data file...')
 1500 FORMAT (1X, 79('-'))
 1600 FORMAT (1X, 'SPECIES',16X,'MOLECULAR',24X,'ELEMENT COUNT',/
     1        1X, 'CONSIDERED',13X,'WEIGHT',7X,'Density',4X,'Nsites',
     2        3X, 15A3)
 1650 FORMAT (I4,'. ', A16, F11.5, 23X, 10I3)
 1800 FORMAT (54X, '(k = A T**b exp(-E/RT))',/,
     1         6X, 'SURFACE REACTIONS CONSIDERED',
     2        22X, 'A',8X,'b',8X,'E',/)
C
C     end of SUBROUTINE SKKEY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNBHM (LOUT, IDIM, LBHM, NBHM, IBHM, KKGAS, KELECT,
     1                   NSPEC, NIISUR, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNBHM
C  initializes Bohm-correction reactions, but additional checking must
C  be done by SKPBHM after all information on the reaction is complete.
C
C  Arguments:
C  LOUT   - Integer scalar, formatted output file unit number.
C  LBHM   - Logical, reaction Bohm-correction status.
C  NBHM   - Integer scalar, count of Bohm-correction reactions
C  IBHM(*)- Integer array, Bohm-correction reaction indices.
C  KKGAS  - Integer scalar, gas-phase species count.
C  KELECT - Integer scalar, species index of the electron species.
C  NSPEC  - Integer scalar, reaction count of reactants+products.
C  NIISUR - Integer scalar, total reaction count and the index of
C           the current reaction.
C  KERR   - Logical, error flag.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C     Bohm velocity correction for collisions of positive ions
C     with a surface
C
      DIMENSION IBHM(IDIM)
      LOGICAL LBHM, KERR
C
      IF (LBHM) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple BOHM declarations...'
         RETURN
      ENDIF
C     increment the counter for the Bohm correction reactions
      NBHM = NBHM + 1
C     add this reaction's index to the NBHM reactions
      IBHM(NBHM) = NIISUR
C
      IF (NSPEC .GT. 0) THEN
          KERR = .TRUE.
          WRITE (LOUT,'(6X,A)')
     1    'Error...BOHM and reversible reaction mutually exclusive...'
      ENDIF
C
      IF (KELECT .LE. 0 .OR. KELECT .GT. KKGAS) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)') 'Error...',
     1   'BOHM requires a gas-phase electron species in mechanism...'
      ENDIF
C
C     end of SUBROUTINE SKNBHM
      RETURN
      END
C
      SUBROUTINE SKPBHM (LOUT, KDIM, IDIM, MAXSPR, MAXORD, NIISUR,
     1                   KKGAS, KNAME, KELECT, KCHRG, SKMIN, NREAC,
     2                   NUNK, NU, NRNU, IRNU, RNU, NORD, IORD,
     2                   KORD, RORD, KBHM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKPBHM
C  checks Bohm-correction reactions after all information about the
C  reaction is complete.
C
C  Arguments:
C  LOUT    - Integer scalar, formatted output file unit number.
C  NIISUR  - Integer scalar, total reaction count and the index
C            of the current reaction.
C  KKGAS   - Integer scalar, gas-phase species count.
C  KNAME(*)- Character-string array, species names.
C  KELECT  - Integer scalar, species index of the electron species.
C  KCHRG(*)- Integer array, species electronic charges.
C  SKMIN   - Real scalar, minimum value for content and balance.
C  NREAC   - Integer scalar, reaction of reactants+products.
C  NUNK(*) - Integer array, reaction species indices.
C  NU(*)   - Integer array, reaction stoichiometric coefficients.
C  MAXSPR  - Integer scalar, maximum reaction species.
C  NRNU    - Integer scalar, real stoichiometry reaction count.
C  IRNU(*) - Integer array, real stoichiometry reaction indices.
C  RNU(*,*)- Real matrix, real stoichiometric coefficients.
C  NORD    - Integer scalar, change-order reaction count.
C  IORD(*) - Integer array, change-order reaction indices.
C  MAXORD  - Integer scalar, maximum reaction change-orders.
C  KORD(*,*)-Integer scalar, reaction change-order species indices.
C  KBHM    - Integer scalar, Bohm-correction species index.
C  KERR    - Logical, error flag.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NUNK(MAXSPR), NU(MAXSPR), IRNU(IDIM),
     1          RNU(MAXSPR,IDIM), IORD(IDIM), KORD(MAXORD,IDIM),
     2          RORD(MAXORD,IDIM), KCHRG(KDIM)
      CHARACTER KNAME(KDIM)*16
      LOGICAL LRNU, LORD, KERR, IERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C     need exactly 1 non-electron positive ion (KBHM) of order approx 1,
      KBHM = 0
      COEF = 0.0
      IERR = .FALSE.
C
      LORD = NIISUR .EQ. IORD(MAX(NORD,1))
      LRNU = NIISUR .EQ. IRNU(MAX(NRNU,1))
C
      IF (.NOT. LORD) THEN
         DO 300 N = 1, NREAC
            IF (NUNK(N).GT.KKGAS .OR. NUNK(N).EQ.KELECT) GO TO 300
C           this is a gas-phase reactant and not the electron species
            IF (KCHRG(NUNK(N)) .GT. 0) THEN
C              this is the gas-phase ion stoichiometric coefficient
               IF (LRNU) THEN
                 COEF = ABS(RNU(N, NRNU))
               ELSE
                  COEF = ABS(NU(N))
               ENDIF
               IF (COEF.LT.1.0+SKMIN .AND. COEF.GT. 1.0-SKMIN) THEN
                  IF (KBHM .NE. 0) THEN
                     IERR = .TRUE.
                  ELSE
                     KBHM = NUNK(N)
                  ENDIF
               ENDIF
            ENDIF
  300    CONTINUE
C
      ELSE
C
C        check change-of-order reactants
         DO 350 N = 1, MAXORD
C           only need to check FORD/reactant change-of-order species
            IF (KORD(N, NORD) .GE. 0) GO TO 350
C
            NK = ABS(KORD(N, NORD))
            IF (NK .LE. KKGAS .AND. NK .NE. KELECT) THEN
C              this is a gas-phase reactant but not the electron
               IF (KCHRG(NK) .GT. 0) THEN
C                 this is a gas-phase ion reactant
                  COEF = RORD(N, NORD)
                  IF (COEF .LT. 1.0+SKMIN .AND. COEF.GT. 1.0-SKMIN) THEN
                     IF (KBHM .NE. 0) THEN
                        IERR = .TRUE.
                     ELSE
                        KBHM = NUNK(N)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
  350    CONTINUE
      ENDIF
C
      IF (KBHM.GT.0 .AND. .NOT.IERR) THEN
         ILEN = CKLSCH(KNAME(KBHM))
         WRITE (LOUT, '(6X,A,A)')
     1   'Bohm correction for ionic species...',KNAME(KBHM)(1:ILEN)
      ELSE
         WRITE (LOUT, 500)
         KERR = .TRUE.
      ENDIF
C
  500 FORMAT (6X,
     2    'Error...BOHM requires exactly one gas-phase positive ',
     3    'ion participant, and its order must be 1.0...')
C
C     end of SUBROUTINE SKPBHM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNCOV (LOUT, KDIM, IDIM, MXPHSE, RSTR, NIISUR, LCOV,
     1                   NCOV, ICOV, KNAME, KKTOT, PNAME, KKPHAS,
     2                   KFSURF, KLSURF, NSCOV, CPAR, IKCOV, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNCOV
C  initializes information about a coverage reaction.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  RSTR     - Character string.
C  NIISUR   - Integer scalar, total reaction count and the index of
C             the current reaction.
C  LCOV     - Logical, coverage reaction status.
C  ICOV(*)  - Integer array, coverage reaction indices.
C  KNAME(*) - Character-string array, species names.
C  KKTOT    - Integer scalar, total species count.
C  PNAME(*) - Character-string array, phase names.
C  KKPHAS(*)- Integer array, phase total species count.
C  KFSURF   - Integer scalar, species index of first site species.
C  KLSURF   - Integer scalar, species index of last site species.
C  NSCOV    - Integer scalar, required number of coverage parameters.
C  CPAR(*,*)- Real matrix, coverage reaction parameters.
C  IKCOV(*) - Integer array, coverage species indices.
C  KERR     - Logical, error flag.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER RSTR*(*), PNAME(MXPHSE)*16, KNAME(KDIM)*16
      DIMENSION ICOV(IDIM), KKPHAS(MXPHSE), CPAR(NSCOV,IDIM),
     1          IKCOV(IDIM), VALUE(5)
      LOGICAL LCOV, KERR, IERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IERR = .FALSE.
C
C      IF (LCOV) THEN
C         IERR = .TRUE.
C         WRITE (LOUT,'(A)') 'Error...multiple COV declarations...'
C      ELSE
C        increment the number of reactions with coverage dependence
         NCOV = NCOV+1
C        add the index for this reaction to the NCOV reactions
         ICOV(NCOV) = NIISUR
C      ENDIF
C
C     find species name and the coverage parameters; we must see if
C     the species' phase name has been appended, e.g. SIH3/PLANE/
C
C     "species/phase/parameters" or "species/parameters/other_key"?
      I2 = INDEX(RSTR,' ')
      IF (INDEX(RSTR,'/') .GT. 0) CALL CKDLIM (RSTR, '/', I1, I2)
      CALL SKPCMP (RSTR(1:I2), KNAME, KKTOT, PNAME, NPHASE,
     1             KKPHAS, KNUM, NF)
C
      ISTART = I1
      IEND = I2
      IF (KNUM .NE. 0) THEN
C        "species/phase/" found, "parameters" in RSTR(I2+1:)
         ISTART = I2+1
         IEND = CKLSCH(RSTR)
C
      ELSE
         CALL SKPCMP(RSTR(1:I1-1), KNAME, KKTOT, PNAME, NPHASE,
     1               KKPHAS, KNUM, NF)
         IF (KNUM .NE. 0) THEN
C           "species" found, "/parameters/" in RSTR(I1:I2)
            ISTART = I1 + 1
            IEND = I2 - 1
         ENDIF
      ENDIF
C
      IF (KNUM .LE. 0) THEN
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...unrecognized COV species...',RSTR(1:I2)
C
      ELSE
         IF (KNUM .LT. KFSURF .OR. KNUM .GT. KLSURF) THEN
            IERR = .TRUE.
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...COV species must be be site-phase...',
     2      RSTR(1:ISTART-1)
         ENDIF
      ENDIF
C
      CALL CKPARR(RSTR(ISTART:IEND), 1, NSCOV, VALUE, NVAL, IER, LOUT)
      IF (IER.NE.0 .OR. NVAL.NE.NSCOV) THEN
         IERR = .TRUE.
         WRITE (LOUT, '(6X, I3, A, A)')
     1   'Error...COV requires ',NSCOV,' parameters...',
     2   RSTR(ISTART:IEND)
      ENDIF
C
      IF (IERR) THEN
         KERR = .TRUE.
      ELSE
C        we have found the species name and its coverage
C        parameters successfully; save the species index
C        and write the information to the output file
C
         IKCOV(NCOV) = KNUM
         ILEN = CKLSCH(KNAME(KNUM))
         DO 10 N = 1, NSCOV
            CPAR(N,NCOV) = VALUE(N)
   10    CONTINUE
         WRITE (LOUT,200) KNAME(KNUM)(1:ILEN),
     1                    (CPAR(J,NCOV),J=1,NSCOV)
      ENDIF
C
  200 FORMAT (6X,'Coverage parameters for species ',A,': ',/,10X,
     1            3(1PE10.3))
C
C     end of SUBROUTINE SKNCOV
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNEDP (LOUT, IDIM, RSTR, NIISUR, LEDP, NEDP, IEDP,
     1                   NEDPAR, PEDP, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNEDP
C  initializes information about an ion-energy dependent reaction;
C  further checking is done in SKPEDP after all information about the
C  reaction has been processed.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  RSTR     - Character string.
C  NIISUR   - Integer scalar, total reaction count and index of the
C             current reaction.
C  LEDP     - Logical, ion-energy reaction status.
C  NEDP     - Integer scalar, ion-energy reactions count.
C  IEDP(*)  - Integer array, ion-energy reactions indices.
C  NEDPAR   - Integer scalar, required number of ion-energy params.
C  PEDP(*,*)- Real matrix, ion-energy reaction parameters.
C  KERR     - Logical, error flag.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) RSTR
      DIMENSION IEDP(IDIM), PEDP(NEDPAR,IDIM)
      LOGICAL KERR, LEDP
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IF (LEDP) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple ENR declarations...'
         RETURN
      ENDIF
C
C     increment the total of reactions having ion-energy dependence
      NEDP = NEDP + 1
C     add the index for this reaction to the NEDP reactons
      IEDP(NEDP) = NIISUR
C
C     get the parameters from the user's input string
      CALL CKPARR (RSTR, 1, NEDPAR, PEDP(1,NEDP), NVAL, IER, LOUT)
      IF (IER.NE.0 .OR. NVAL.NE.NEDPAR) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(RSTR)
         WRITE (LOUT, '(6X,A,I3,A)')
     1   'Error...ENR requires ',NEDPAR,' parameters...',RSTR(1:ILEN)
      ELSE
         WRITE (LOUT, 600)
         WRITE (LOUT, 900) 'Energy-dependence coefficients:',
     1                     (PEDP(N,NEDP), N=1,NEDPAR)
      ENDIF
C
  600 FORMAT (6X,'Energy-dependence of reaction rate for ionic species',
     1           ' will be applied...')
  900 FORMAT (6X,A, T53, 1PE8.2, 2X, 0PF5.1, 2X, F9.1)
C
C     end of SUBROUTINE SKNEDP
      RETURN
      END
C
      SUBROUTINE SKPEDP (LOUT, KDIM, IDIM, NIISUR, KKGAS, KNAME, KCHRG,
     1                   SKMIN, NREAC, NUNK, NU, MAXSPR, NRNU, IRNU,
     2                   RNU, NORD, IORD, MAXORD, KORD, RORD, KEDP,
     3                   KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKPEDP
C  checks ion-energy-dependence reactions for completeness after all
C  information about the reaction is processed.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  NIISUR   - Integer scalar, total reaction count and the index
C             of current reaction.
C  KKGAS    - Integer scalar, gas-phase species count.
C  KNAME(*) - Character-string array, species names.
C  KCHRG(*) - Integer array, species electronic charges.
C  SKMIN    - Real scalar, minimum value for content and balancing.
C  NREAC    - Integer scalar, reaction count of reactants+products.
C  NUNK(*)  - Integer array, reaction species indices.
C  NU(*)    - Integer array, reaction stoichiometric coefficients.
C  MAXSPR   - Integer scalar, maximum number of reaction species.
C  NRNU     - Integer scalar, real stoichiometry reaction count.
C  IRNU(*)  - Integer array, real stoichiometry reaction indices.
C  RNU(*,*) - Real matrix, real stoichiometric coefficients.
C  NORD     - Integer scalar, change-order reaction count.
C  IORD(*)  - Integer array, change-order reaction indices.
C  MAXORD   - Integer scalar, maximum reaction change-orders.
C  KORD(*,*)- Integer matrix, reaction change-order species indices.
C  RORD(*,*)- Real matrix, reaction change-order species values.
C  KEDP     - Integer scalar, temperature-dependent species.
C  KERR     - Logical, error flag.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KCHRG(KDIM), NUNK(MAXSPR), NU(MAXSPR), IRNU(IDIM),
     1          RNU(MAXSPR,IDIM), IORD(IDIM), KORD(MAXORD,IDIM),
     2          RORD(MAXORD,IDIM)
      CHARACTER KNAME(KDIM)*16
      LOGICAL KERR, LRNU, LORD
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C     there must be exactly one gas-phase reactant and it must
C     be an ionic species (KEDP)
C
      NION = 0
      KEDP = 0
      COEF = 0.0
C
      LRNU = NIISUR .EQ. IRNU(MAX(NRNU,1))
      LORD = NIISUR .EQ. IORD(MAX(NORD,1))
C
      DO 300 N = 1, NREAC
         IF (NUNK(N).GT.KKGAS .OR. KCHRG(N).LE.0) GO TO 300
C        this is a gas-phase ion species index
         KEDP = NUNK(N)
         NION = NION + 1
C        this is a gas-phase ion stoichiometric coefficient
         IF (LRNU) THEN
            COEF = ABS(RNU(N, NRNU))
         ELSE
            COEF = ABS(NU(N))
         ENDIF
  300 CONTINUE
C
      IF (LORD) THEN
C
C        check change-of-order reactants
         DO 350 N = 1, MAXORD
C           only check FORD/reactant change-of-order species
            IF (KORD(N, NORD) .GE. 0) GO TO 350
            NK = ABS(KORD(N, NORD))
            IF (NK.LE.KKGAS .AND. KCHRG(NK).GT.0) THEN
C              this is a gas-phase ion reactant
               NION = NION + 1
               IF (NK .EQ. KEDP) THEN
C                 new coefficient for the gas-phase ion reactant
                  COEF = RORD(N, NORD)
               ELSEIF (KEDP .EQ. 0) THEN
C                this is the first gas-phase ion reactant found
                  KEDP = NK
C                 this is its coefficient
                  COEF = RORD(N, NORD)
               ENDIF
            ENDIF
  350    CONTINUE
      ENDIF
C
      IF (NION .NE. 1) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1   'Error...ENR requires exactly one gas-phase ion reactant...'
      ELSE
         IF (COEF.GT.1.0+SKMIN .OR. COEF.LT.1.0-SKMIN) THEN
            KERR = .TRUE.
            ILEN = CKLSCH(KNAME(KEDP))
            WRITE (LOUT, '(6X,A,A,A)')
     1      'Error...ENR gas-phase ion ',KNAME(KEDP)(1:ILEN),
     2      ' must have stoichiometric coefficient of 1...'
         ENDIF
      ENDIF
C
C     end of SUBROUTINE SKPEDP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNORD (LOUT, KDIM, IDIM, KEY, RSTR, NIISUR, KKTOT,
     1                   KNAME, PNAME, NPHASE, KKPHAS, NORD, MAXORD,
     2                   IORD, KORD, RORD, MAXSPR, NUNK, NU, NREAC,
     2                   NSPEC, LORD, LRNU, NRNU, RNU, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNORD
C  processes information about changed-order species in reactions.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  KEY      - Character string, FORD or NORD change-order option.
C  RSTR     - Character string, species and change-order value.
C  NIISUR   - Integer scalar, total reaction count and the index
C             of this reaction.
C  KKTOT    - Integer scalar, total species count.
C  KNAME(*) - Character-string array, species names.
C  PNAME(*) - Character-string array, phase names.
C  NPHASE   - Integer scalar, total phase count.
C  KKPHAS(*)- Integer array, phase species cïcounts.
C  NORD     - Integer scalar, total count of change-order reactions.
C  MAXORD   - Integer scalar, maximum number of reaction change-orders.
C  IORD(*)  - Integer array, change-order reaction indices.
C  KORD(*,*)- Integer matrix, reactions change-order species.
C  RORD(*,*)- Real matrix, reactions change-order values.
C  MAXSPR   - Integer scalar, maximum number of reaction species.
C  NUNK(*)  - Integer array, the reactions species indices.
C  NU(*)    - Integer array, the reactions stoichiometry coeff'nts.
C  NREAC    - Integer scalar, reaction count of reactants+products.
C  NSPEC    - Integer scalar, reaction count of reactants only.
C  LORD     - Logical, flag of reaction's change-order status.
C  LRNU     - Logical, flag of reaction's real stoichiometry status.
C  NRNU     - Integer scalar, total count of real stoichiometry reaxs.
C  RNU(*,*) - Real matrix, real stoichiometric coefficients.
C  KERR     - Logical error flag.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*16 KNAME(KDIM), PNAME(NPHASE), KEY*(*), RSTR*(*)
      DIMENSION IORD(IDIM), NUNK(MAXSPR), NU(MAXSPR), KORD(MAXORD,IDIM),
     1          RORD(MAXORD,IDIM), RNU(MAXSPR,IDIM), KKPHAS(NPHASE)
      LOGICAL LORD, LFORD, LRORD, LRNU, KERR, IERR
      INTEGER CKLSCH, CKSLEN
      EXTERNAL CKLSCH, CKSLEN
C
      IERR = .FALSE.
C
C     LRORD is reverse change-of-order,
C     LFORD is forward reaction change-of-order
C
      LFORD = KEY(1:1) .EQ. 'F'
      LRORD = KEY(1:1) .EQ. 'R'

      IF (.NOT. LORD) THEN
C
C        need to initialize this change-of-order reaction;
C        increment the number of reactions with non-standard orders
         NORD = NORD + 1
C        add this reaction index to the list of NORD reactions
         IORD(NORD) = NIISUR
C
C        copy original stoichiometry
         NKORD = 0
         NKORD = ABS(NREAC)
         DO 100 N = 1, NKORD
C           reactants
            KORD(N, NORD) = -NUNK(N)
            IF (LRNU) THEN
               RORD(N, NORD) = ABS(RNU(N,NRNU))
            ELSE
               RORD(N, NORD) = IABS(NU(N))
            ENDIF
  100    CONTINUE
         DO 110 N = MAXSPR/2 + 1, MAXSPR
C           products
            IF (NUNK(N) .NE. 0) THEN
               NKORD = NKORD + 1
               KORD(NKORD, NORD) = NUNK(N)
               IF (LRNU) THEN
                  RORD(NKORD, NORD) = RNU(N,NRNU)
               ELSE
                  RORD(NKORD, NORD) = NU(N)
               ENDIF
            ENDIF
  110    CONTINUE
      ENDIF
C
C     need to split species name from change-order value
      KEND = 0
      ILEN = CKSLEN(RSTR)
      DO 55 I = 2, ILEN
         IF (KEND .EQ. 0) THEN
            IF (RSTR(I:I).EQ.' ' .AND. RSTR(I-1:I-1).NE.' ') KEND=I-1
            IF (RSTR(I:I).NE.' ' .AND. I.EQ.ILEN) KEND=I
         ENDIF
   55 CONTINUE
C
      KNUM = 0
      IF (KEND .GT. 0)
     1    CALL SKPCMP (RSTR(1:KEND), KNAME, KKTOT, PNAME, NPHASE,
     2             KKPHAS, KNUM, NF)
C
      IF (KNUM .EQ. 0) THEN
C        could not get species; issue error message and set flag
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A,A,A,A)')
     1      'Error...unrecognized ',KEY,' species...', RSTR(1:ILEN)
      ENDIF
C
      NVAL = 0
      IF (KEND .GT. 0)
     1   CALL CKPARR (RSTR(KEND+1:), 1, 1, VAL, NVAL, IER, LOUT)
C
      IF (IER.NE.0 .OR. NVAL.NE.1) THEN
C        could not get good order value for species
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A,A,A,A)')
     1   'Error...problem finding ',KEY,' parameter...', RSTR(1:ILEN)
      ENDIF
C
      IF (LRORD .AND. NSPEC.LT.0) THEN
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...RORD incompatible with irreversible reaction...',
     2   RSTR(1:ILEN)
      ENDIF
C
      IF (KNUM .NE. 0) THEN
C        store species index; negative for reactants
         IF (LFORD) KNUM = -KNUM
C
         NPOS =0
         DO 200 N = 1, MAXORD
            NK = KORD(N, NORD)
            IF (NPOS.EQ.0 .AND. (NK.EQ.KNUM .OR. NK.EQ.0))
C              this is the place to put the change-order
     1         NPOS = N
  200    CONTINUE
C
         IF (NPOS .EQ. 0) THEN
C           did not find the species or an empty space
            IERR = .TRUE.
            WRITE (LOUT, '(6X,A,I5,A)')
     1      'Error...more than ',MAXORD,' change-orders...',
     2      RSTR(1:ILEN)
         ENDIF
      ENDIF
C
      IF (IERR) THEN
         KERR = .TRUE.
      ELSE
         KORD(NPOS, NORD) = KNUM
         RORD(NPOS, NORD) = VAL
         KNUM = ABS(KNUM)
         ILEN = CKLSCH(KNAME(KNUM))
         IF (LFORD) THEN
            WRITE (LOUT, '(6X, A, A, 1PE12.3)')
     1      'Forward order ', KNAME(KNUM)(1:ILEN), VAL
         ELSE
            WRITE (LOUT, '(6X, A, A, 1PE12.3)')
     1      'Reverse order ', KNAME(KNUM)(1:ILEN), VAL
         ENDIF
      ENDIF
C
C     end of SUBROUTINE SKNORD
      RETURN
      END
C                                                                      C
      SUBROUTINE SKMOTZ (LOUT, IDIM, LSTK, NSTK, SUB, RSTR, MOTZ)
C
C  START PROLOGUE
C
C  SUBROUTINE SKMOTZ
C  checks Motz-Wise correction requested for reaction.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  LSTK     - Logical flag, sticking reaction status.
C  NSTK     - Integer scalar, total count of sticking reactions.
C  SUB      - Character string, the Motz-Wise option input string.
C  RSTR     - Character string, the one-character N/F (on/off) option.
C  MOTZ(*)  - Integer array, sticking reactions Motz-Wise flags.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION MOTZ(IDIM)
      CHARACTER SUB*(*), RSTR*(*)
      LOGICAL LSTK
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      ILEN = CKLSCH(SUB)
      IF (RSTR  .EQ. 'N') THEN
         IF (.NOT. LSTK) THEN
            WRITE (LOUT, '(6X,A,A)')
     1      'Warning...Motz-Wise correction will be ON only if ',
     2      'STICK found for this reaction...',SUB(1:ILEN)
         ELSE
            WRITE (LOUT, '(6X,A,A)')
     1      'Motz-Wise sticking rate correction turned ON...',
     2      SUB(1:ILEN)
            MOTZ(NSTK) = 1
         ENDIF
      ELSEIF (RSTR .EQ. 'F') THEN
C        Motz-Wise correction to sticking coefficients turned-off
C        for this reaction
         IF (.NOT. LSTK) THEN
            WRITE (LOUT, '(6X,A,A)')
     1      'Warning...Motz-Wise correction will be OFF only if ',
     2      'STICK found for this reaction...',SUB(1:ILEN)
         ELSE
            WRITE (LOUT, '(6X,A,A)')
     1      'Motz-Wise sticking rate correction turned OFF...',
     2      SUB(1:ILEN)
            MOTZ(NSTK) = 0
         ENDIF
      ENDIF
C
C     end of SUBROUTINE SKMOTZ
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNREV (LOUT, IDIM, NIISUR, RSTR, LREV, NSPEC, NSPAR,
     1                   NREV, IREV, RSPAR, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNREV
C  initializes information about a reaction when explicit reverse
C  Arrhenius parameters are given.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  NIISUR   - Integer scalar, total reaction count and the index
C             of this reaction.
C  RSTR     - Character string, representing reverse parameters.
C  LREV     - Logical, reaction reverse parameter status.
C  NSPEC    - Integer scalar, reaction count of reactants+products.
C  NSPAR    - Integer scalar, required number of Arrhenius parameters.
C  NREV     - Integer scalar, total reaction count of those with
C             explicit reverse parameters.
C  IREV(*)  - Integer array, reaction indices for those with
C             explicit reverse parameters.
C  RSPAR(*,*)-Real matrix, reactions explicit reverse parameters.
C  KERR     - Logical error flag.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) RSTR
      DIMENSION IREV(IDIM), RSPAR(NSPAR,IDIM)
      LOGICAL KERR, LREV
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IF (LREV) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple REV declarations...'
         RETURN
      ENDIF
C     increment number of reactions for which reverse
C     rate coefficients have been specified
      NREV = NREV + 1
C     add this reaction's index to the list of NREV reactions
      IREV(NREV) = NIISUR
C
      IF (NSPEC .LT. 0) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1   'Error...REV and irreversible reaction mutually exclusive...'
      ENDIF
C
C     extract the value of the coefficients
C
      CALL CKPARR (RSTR, 1, NSPAR, RSPAR(1,NREV), NVAL, IER, LOUT)
      IF (IER.NE.0 .OR. NVAL.NE.NSPAR) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(RSTR)
         WRITE (LOUT, '(6X,A,I3,A)')
     1   'Error...REV requires ',NSPAR,' parameters...',RSTR(1:ILEN)
      ELSE
         WRITE (LOUT, 100) (RSPAR(N,NREV), N=1,NSPAR)
      ENDIF
      IF (NSPAR.GE.1 .AND. RSPAR(1,NREV).LT.0.0) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(RSTR)
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...REV pre-exponential must be positive...',RSTR(1:ILEN)
      ENDIF
C
  100 FORMAT (6X,'Reverse Arrhenius coefficients:',
     1           T53, 1PE8.2, 2X, 0PF5.1, 2X, F9.1)
C
C     end of SUBROUTINE SKNREV
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNSTK (LOUT, IDIM, LSTK, NSTK, ISTK, LIMOTZ, MOTZ,
     1                   NIISUR, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNSTK
C  intializes information about a sticking reaction;
C  further checking will be done in SKPSTK after all the reactions
C  auxiliary information is processed.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  LSTK     - Logical, sticking reaction status.
C  NSTK     - Integer scalar, total sticking reaction count.
C  ISTK(*)  - Integer array, sticking reaction indices.
C  LIMOTZ   - Logical, mechanism Motz-Wise correction flag.
C  MOTZ(*)  - Integer array, sticking reactions Motz-Wise flags.
C  NIISUR   - Integer scalar, total reaction count and the index
C             of the current reaction.
C  KERR     - Logical error flag.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISTK(IDIM), MOTZ(IDIM)
      LOGICAL LSTK, LIMOTZ, KERR
C
      IF (LSTK) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple STICK declarations...'
C        have previously checked reactants
      ELSE
         NSTK = NSTK + 1
         ISTK(NSTK) = NIISUR
         IF (LIMOTZ) THEN
            MOTZ(NSTK) = 1
         ELSE
            MOTZ(NSTK) = 0
         ENDIF
      ENDIF
C
C     end of SUBROUTINE SKNSTK
      RETURN
      END
C
      SUBROUTINE SKPSTK (LOUT, KDIM, IDIM, NIISUR, KKGAS, KNAME, SKMIN,
     1                  NREAC, NUNK, NU, MAXSPR, NRNU, IRNU, RNU, NORD,
     2                  IORD, MAXORD, KORD, RORD, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKPSTK
C  completes checking of a sticking reaction after all of its
C  auxiliary data is processed.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  NIISUR   - Integer scalar, total reaction count and the index of
C             the current reaction.
C  KKGAS    - Integer scalar, total gas-phase species count.
C  KNAME(*) - Character-string array, species names.
C  SKMIN    - Real scalar, minimum value for balance, content.
C  NREAC    - Integer scalar, reaction's count of reactants+products.
C  NUNK(*)  - Integer array, reaction's stoichiometric coefficients.
C  NU(*)    - Integer array, reaction's species indices.
C  MAXSPR   - Integer scalar, maximum number of reaction species.
C  NRNU     - Integer scalar, count of real stoichiometry reactions.
C  IRNU(*)  - Integer array, real stoichiometry reaction indices.
C  RNU(*,*) - Real matrix, real stoichiometry coefficients.
C  NORD     - Integer scalar, count of change-order reactions.
C  IORD(*)  - Integer array, change-order reaction indices.
C  MAXORD   - Integer scalar, maximum reaction species change-orders.
C  KORD(*,*)- Integer matrix, reaction change-order species indices.
C  RORD(*,*)- Real matrix, reaction change-order values.
C  KERR     - Logical error flag.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NUNK(MAXSPR), NU(MAXSPR), IRNU(IDIM), RNU(MAXSPR,IDIM),
     1          IORD(IDIM), KORD(MAXORD,IDIM), RORD(MAXORD,IDIM)
      CHARACTER KNAME(KDIM)*16
      LOGICAL LRNU, LORD, KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      LRNU = NIISUR .EQ. IRNU(MAX(NRNU,1))
      LORD = NIISUR .EQ. IORD(MAX(NORD,1))
C
C     need to cross-check every gas-phase reactant and every FORD
C     gas-phase change-order species to find exactly one gas species
C     (KGAS) with a coefficient approx 1

      KGAS = 0
      COEF = 0.0
C
      DO 300 N = 1, NREAC
C
         KTEST = NUNK(N)
         IF (KTEST.EQ.0 .OR. KTEST.GT.KKGAS) GO TO 300
C
C        this is a gas-phase reactant
         IF (LRNU) THEN
C           this is it's coefficient
            CTEST = ABS(RNU(N,NRNU))
         ELSE
            CTEST = ABS(NU(N))
         ENDIF
C
C        does KTEST change order; if so, get the last CTEST?
         IF (LORD) THEN
            DO 310 M = 1, MAXORD
C              at this point, only need to look for KTEST in the
C              FORDs for changes to CTEST
               IF (KORD(M,NORD) .LT. 0 .AND.
     1             ABS(KORD(M,NORD)) .EQ. KTEST) CTEST = RORD(M,NORD)
  310       CONTINUE
         ENDIF
C
C        is the final order of KTEST approx. 1.0?
         IF (CTEST.LE.1.0+SKMIN .AND. CTEST.GE.1.0-SKMIN) THEN
C
C           KTEST is a gas-phase species with coef approx. 1.0;
            IF (KGAS.GT.0 .AND. COEF.NE.0.0) THEN
C              previously-found KTEST and CTEST; could be same
               IF (KTEST .EQ. KGAS) GO TO 300
C
               KERR = .TRUE.
               WRITE (LOUT, '(6X,A,A)')
     1         'Error...STICK requires exactly one active gas-phase ',
     2         'participant...'
            ELSE
               KGAS = KTEST
               COEF = CTEST
               KTEST = 0
               CTEST = 0.0
            ENDIF
         ENDIF
  300 CONTINUE
C
C     may or may not have a KGAS>0 and COEF approx. 1.0;
C     in any case, need to look at other gas-phase change-orders
C
      IF (LORD) THEN
C
         DO 350 N = MAXORD, 1, -1
            KTEST = KORD(N,NORD)
C           don't need KTEST if 0, if KTEST>0 (product), or if site/bulk
            IF (KTEST .GE. 0 .OR. ABS(KTEST) .GT. KKGAS) GO TO 350
C           don't need KTEST if matches KGAS
            IF (ABS(KTEST) .EQ. KGAS) GO TO 350
C
            KTEST = ABS(KTEST)
            CTEST = RORD(N,NORD)
            IF (CTEST.LE.1.0+SKMIN .AND. CTEST.GE.1.0-SKMIN) THEN
C
C              this is a gas-phase species with coef approx. 1.0;
C              is it the same as previously-found?
C
               IF (KGAS .GT. 0 .AND. COEF.NE.0.0) THEN
C                 this is a not-allowed additional gas-phase reactant
                  KERR = .TRUE.
                  WRITE (LOUT, '(6X,A,A)')
     1            'Error...STICK requires exactly one active ',
     2            'gas-phase participant...'
               ELSE
                  KGAS = KTEST
                  COEF = CTEST
                  KTEST = 0
                  CTEST = 0.0
               ENDIF
            ENDIF
  350    CONTINUE
      ENDIF
C
      IF (KGAS .EQ. 0) THEN
C        no gas-phase found
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...STICK requires a gas-phase reactant...'
C
      ELSE
          IF (COEF.GT.1.0+SKMIN .OR. COEF.LT.1.0-SKMIN) THEN
C              gas-phase reactant coefficient must be 1
               KERR = .TRUE.
               ILEN = CKLSCH(KNAME(KGAS))
               WRITE (LOUT, '(6X,A,A)')
     1         'Error...STICK gas-phase reactant ',
     2         KNAME(KGAS)(1:ILEN), ' must have a coefficient of 1...'
         ELSE
            WRITE (LOUT, '(6X,A)')
     1      'Coefficients are sticking parameters...'
         ENDIF
      ENDIF
C
C     end of SUBROUTINE SKPSTK
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNYLD (LOUT, IDIM, RSTR, NIISUR, NSPEC, KKGAS,
     1                   LYLD, NYLD, IYLD, NYPAR, PYLD, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNYLD
C  initializes information about a reaction with species yield
C  modify; final check will be done by SKPYLD after all its
C  auxiliary data has been processed.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  RSTR     - Character string, representing yield parameters.
C  NIISUR   - Integer scalar, total reaction count and the index
C             of the current reaction.
C  NSPEC    - Integer scalar, reaction's count of reactants+products.
C  LYLD     - Logical, reaction's yield-modify status.
C  NYLD     - Integer scalar, yield-modify reaction count.
C  IYLD(*)  - Integer array, yield-modify reaction indices.
C  NYPAR    - Integer scalar, required number of yield parameters.
C  PYLD(*,*)- Real matrix, yield-modify reaction yield parameters.
C  KERR     - Logical error flag.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) RSTR
      DIMENSION IYLD(IDIM), PYLD(NYPAR,IDIM)
      LOGICAL LYLD, KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IF (.NOT. LYLD) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(RSTR)
         WRITE (LOUT,'(6X,A,A)')
     1   'Error...YIELD requires # before effected reaction species...',
     2   RSTR(1:ILEN)
         NYLD = NYLD + 1
         IYLD(NYLD) = NIISUR
C
C        the following are checked in SKREAC if #-sign had been found
C        need irreversible reaction
         IF (NSPEC .GT. 0) THEN
            KERR = .TRUE.
            ILEN = CKLSCH(RSTR)
            WRITE (LOUT, '(6X,A,A)')
     1   'Error...YIELD and reversible reaction mutually exclusive...',
     2      RSTR(1:ILEN)
         ENDIF
      ENDIF
C
      CALL CKPARR (RSTR, 1, NYPAR, PYLD(1,NYLD), NVAL, IER, LOUT)
      IF (IER .NE. 0 .OR. NVAL.NE.NYPAR) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(RSTR)
         WRITE (LOUT, '(6X,A,I3,A,A)')
     1   'Error...YIELD requires ',NYPAR, ' parameters...',
     2   RSTR(1:ILEN)
      ELSE
C
C        PYLD(1) = A, PYLD(2) = Eth[eV], PYLD(3) = a, PYLD(4) = b
C        # = A(Ei^a - Eth^a)^b
C        then nu(K) = # * nu(K)
C
         WRITE (LOUT, 100) (PYLD(J,NYLD),J=1,NYPAR)
      ENDIF
C
  100 FORMAT (6X,'#Yield parameters: ',
     1    'A=',1PE10.3,', Eth=',1PE10.3,', a=',0PF4.1,', b=',F4.1)
C
C     end of SUBROUTINE SKNYLD
      RETURN
      END
C                                                                      C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKPYLD (LOUT, KDIM, IDIM, NIISUR, KKGAS, KNAME, KCHRG,
     1                   SKMIN, NREAC, NUNK, NU, MAXSPR, NRNU, IRNU,
     2                   RNU, NORD, IORD, MAXORD, KORD, RORD, KELECT,
     3                   KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKPYLD
C  finalizes and checks the information about a yield-modify reaction,
C  after all of its auxiliary data has been processed.
C
C  Arguments:
C  LOUT     - Integer scalar, formated output file unit number.
C  NIISUR   - Integer scalar, total reaction count and the index
C             of the current reaction.
C  KKGAS    - Integer scalar, total gas-phase species count.
C  KNAME(*) - Character-string array, species names.
C  KCHRG(*) - Integer array, species physical state.
C  SKMIN    - Real scalar, minimum value for content and balancing.
C  NREAC    - Integer scalar, reaction's count of reactants only.
C  NUNK(*)  - Integer array, reaction's species indices.
C  NU(*)    - Integer array, reaction's stoichiometric coefficients.
C  MAXSPR   - Integer scalar, maximum number of reaction species.
C  NRNU     - Integer scalar, count of real stoichiometry reactions.
C  IRNU(*)  - Integer array, real stoichiometry reactions indices.
C  RNU(*,*) - Real matrix, real stoichiometry coefficients.
C  NORD     - Integer scalar, total count of change-order reactions.
C  IORD(*)  - Integer array, change-order reaction indices.
C  MAXORD   - Integer scalar, maximum reaction change-order species.
C  KORD(*,*)- Integer matrix, reactions change-order species indices.
C  RORD(*,*)- Real matrix, reactions change-order values.
C  KELECT   - Integer scalar, species index of an electron species.
C  KERR     - Logical error flag.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KCHRG(KDIM), NUNK(MAXSPR), NU(MAXSPR), IRNU(IDIM),
     1          RNU(MAXSPR,IDIM), IORD(IDIM), KORD(MAXORD,IDIM),
     2          RORD(MAXORD,IDIM)
      CHARACTER KNAME(KDIM)*16
      LOGICAL LRNU, LORD, KERR
C
C     need gas-phase non-electron ion
      NION = 0
      LRNU = NIISUR .EQ. IRNU(MAX(NRNU,1))
      LORD = NIISUR .EQ. IORD(MAX(NORD,1))
C
      DO 50 N = 1, NREAC
         NK = NUNK(N)
C        if there is no electron species in the mechanism, KELECT is 0
C        if there is an electron species in the mechanism, the
C        ion we want is NOT KELECT
         IF (NK.NE.0 .AND. NK.NE.KELECT) THEN
            IF (KCHRG(NK) .NE. 0) NION = NION + 1
         ENDIF
   50 CONTINUE
C
      IF (NION .LT. 1) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...YIELD requires a gas-phase ion reactant...'
      ENDIF
C
C     end of SUBROUTINE SKPYLD
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKPRNT (LOUT, KDIM, MXPHSE, MAXTP, KNAME, ITHRM, TMP,
     1                   NT, TMID, WT, PDEN, LPDEN, DEN, LKDEN, MDIM,
     3                   KNCF, KCOV, KFIRST, KLAST, KKPHAS, PNAME, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKPRNT
C  outputs formatted surface and solid species data.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  KDIM     - Integer scalar, maximum number of species.
C  MXPHSE   - Integer scalar, maximum number of phases.
C  MAXTP    - Integer scalar, maximum number of fit temperatures.
C  KNAME(*) - Character-string array, species names.
C  ITHRM(*) - Logical array, species thermodynamics flags.
C  TMP(*,*) - Real matrix, species fit temperatures.
C  NT(*)    - Integer array, species number of fit temperatures.
C  TMID     - Real scalar, default middle fit temperature.
C  WT(*)    - Real array, species molecular weights.
C  PDEN(*)  - Real array, phase densities.
C  LPDEN(*) - Logical array, phase density flags.
C  DEN(*)   - Real array, species densities.
C  LKDEN(*) - Logical array, species density flags.
C  MDIM     - Integer scalar, maximum number of elements.
C  KNCF(*,*)- Integer matrix, species elemental composition.
C  KCOV(*)  - Integer array, species site coverage.
C  KFIRST(*)- Integer array, phases first species indices.
C  KLAST(*) - Integer array, phases last species indices.
C  KKPHAS(*)- Integer array, phases total species count.
C  PNAME(*) - Character-string array, phase names.
C  KERR     - Logical error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
C
C     Integer arrays
      DIMENSION KNCF(MDIM,KDIM), KCOV(KDIM), NT(KDIM), KFIRST(MXPHSE),
     1          KLAST(MXPHSE), KKPHAS(MXPHSE)
C     Real arrays
      DIMENSION WT(KDIM), PDEN(MXPHSE), DEN(KDIM), TMP(MAXTP,KDIM)
C     Character arrays
      CHARACTER*16 PNAME(MXPHSE), KNAME(KDIM)
C     Logical arrays
      LOGICAL ITHRM(KDIM), LPDEN(MXPHSE), LKDEN(KDIM), KERR
C
      IF (NNSURF .GT. 0) THEN
C        print information about each surface phase
         DO 10 N = NFSURF,NLSURF
C
            IF (LPDEN(N)) THEN
C              print the phase name and site density
               WRITE (LOUT,'(/2(1X,A),10X,E10.3,A)')
     1         'SITE:', PNAME(N), PDEN(N), ' moles/cm**2'
               IF (PDEN(N) .LE. 0.0) THEN
C                 site density found was not a positive number;
                  WRITE (LOUT, '(6X,A)')
     1            'Error...site phase density must be > 0'
                  KERR = .TRUE.
               ENDIF
            ELSE
C              the site density was not supplied;
               WRITE (LOUT,'(/2(1X,A))') 'SITE:', PNAME(N)
               WRITE (LOUT, '(6X,A)')
     1         'Error...site phase density is required input'
               KERR = .TRUE.
            ENDIF
C
            IF (KKPHAS(N) .LE. 0) THEN
C              this phase has no species, which is considered an error;
C              issue an error message and set the error flag
               WRITE (LOUT,'(6X,A)') 'Error...site phase has no species'
               KERR = .TRUE.
            ELSE
C              print out the species names, molecular weight, number
C              of sites occupied, and the elemental composition
               DO 5 K = KFIRST(N), KLAST(N)
                  WRITE (LOUT, 1655) K, KNAME(K), WT(K), KCOV(K),
     1                               (KNCF(M,K),M=1,NELEM)
                  IF (KCOV(K) .LE. 0) THEN
C                    species site coverage must be positive
                     WRITE (LOUT, '(6X,A)')
     1               'Error...number of sites occupied must be > 0'
                     KERR = .TRUE.
                  ENDIF
C                 do several error checks on the validity of the
C                 information accumulated on this species
                  CALL SKSPEC (KDIM, K, KNAME, KKTOT, ITHRM, MAXTP,
     1                         TMP, NT, TMID, LOUT, KERR)
    5          CONTINUE
            ENDIF
   10    CONTINUE
      ENDIF
C
      IF (NNBULK .LE. 0) RETURN
C
C     print information about each bulk phase
      DO 20 N = NFBULK, NLBULK
         WRITE (LOUT, '(/2(1X,A))') 'BULK:', PNAME(N)
         IF (KKPHAS(N) .LE. 0) THEN
C           this phase has no species, which is considered an error;
C           issue an error message and set the error flag
            WRITE (LOUT, '(6X,A)') 'Error...bulk phase has no species'
            KERR = .TRUE.
         ELSE
C
C           print out the species names, molecular weight, density,
C           and the elemental composition
            DO 15 K = KFIRST(N), KLAST(N)
               WRITE (LOUT, 1660) K, KNAME(K), WT(K),
     1                            DEN(K), (KNCF(M,K), M=1,NELEM)
               IF (LKDEN(K) .AND. DEN(K).LE.0.0) THEN
C                 density must be positive for bulk species
                  WRITE (LOUT, '(6X,A)')
     1            'Error...bulk species density must be > 0'
                  KERR = .TRUE.
               ENDIF
C              do several error checks on the validity of the
C              information accumulated on this species
               CALL SKSPEC (KDIM, K, KNAME, KKTOT, ITHRM, MAXTP, TMP,
     1                      NT, TMID, LOUT, KERR)
   15       CONTINUE
         ENDIF
   20 CONTINUE
C
 1655 FORMAT (I4,'. ', A16, F11.5, 17X, I3, 3X, 10I3)
 1660 FORMAT (I4,'. ', A16, F11.5, E10.3,' g/cm**3', 5X, 10I3)
C
C     end of SUBROUTINE SKPRNT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKREAC
     1        (LOUT, KDIM, IDIM, MXPHSE, MAXSPR, NSPAR, MXCONT, LINE,
     2         NLINES, KNAME, PNAME, KKPHAS, NSPEC, NREAC, NUNK, NU,
     3         NUSUMK, SPAR, IRNU, RNU, IYLD, KYION, KYLD, KCHRG, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKREAC
C  processes character strings (lines of text) which make up a surface
C  reaction, to find reactants, products, their coefficients,
C  Arrhenius parameters, and symbols which may indicate a type of
C  rate formulation or other options.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  KDIM     - Integer scalar, maximum number of species.
C  IDIM     - Integer scalar, maximum number of reactions.
C  MXPHSE   - Integer scalar, maximum number of phases.
C  MAXSPR   - Integer scalar, maximum number of reaction species.
C  NSPAR    - Integer scalar, required number of Arrhenius parameters.
C  MXCONT   - Integer scalar, maximum number of character-strings.
C  LINE(*)  - Character string array, the text of a reaction.
C  NLINES   - Integer scalar, character-string line count.
C  KNAME(*) - Character string array, species names.
C  PNAME(*) - Character string array, phase names.
C  KKPHAS(*)- Integer array, phase species counts.
C  NSPEC    - Integer scalar, reaction's count of reactants+products.
C  NREAC    - Integer scalar, reaction's count of reactants only.
C  NUNK(*)  - Integer array,  reaction's species indices.
C  NU(*)    - Integer array, reaction's stoichiometric coefficients.
C  NUSUMK   - Integer scalar, sum of the stoichiometric coefficients.
C  SPAR(*)  - Real array, reaction's Arrhenius parameters.
C  IRNU(*)  - Integer array, reaction indices for those with real
C             stoichiometry.
C  RNU(*,*) - Real matrix, stoichiometric coefficients for reactions
C             with real stoichiometry.
C  IYLD(*)  - Integer array, reaction indices for yield-modify
C             reacitons.
C  KYION(*) - Integer array, ion species indices for yield-modify
C             reactions.
C  KYLD(*,*)- Integer matrix, yield-modify flags for the species in
C             yield-modify reactions.
C  KCHRG(*) - Integer array, species electronic charges.
C
C  END PROLOGUE
C----------------------------------------------------------------------C
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C----------------------------------------------------------------------C
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
C
C     Integer arrays
      DIMENSION NUNK(MAXSPR), NU(MAXSPR), KKPHAS(MXPHSE), IRNU(IDIM),
     1          IYLD(IDIM), KYION(IDIM), KYLD(MAXSPR,IDIM), KCHRG(KDIM),
     2          IPLUS(20)
C     Real arrays
      DIMENSION SPAR(NSPAR), RNU(MAXSPR,IDIM)
C     Character arrays
      CHARACTER*16 KNAME(KDIM), PNAME(MXPHSE)
      CHARACTER*80 RSTR, PSTR, ISTR, ITEMP, ISPEC, LINE(MXCONT)
      LOGICAL KERR, IERR, LRNU, LREV, LYLD
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C     find last 3 character strings and try to convert to the
C     3 Arrhenius parameters
C
      IERR = .FALSE.
      NWANT = NSPAR
C     loop backwards over number of lines for reaction
      DO 20 L = NLINES, 1, -1
C
C        keep track of this line number
         LASTLINE = L
C
C        loop backwards over number of characters in the line
         DO 20 N = CKLSCH(LINE(L)), 1, -1
            RSTR = ' '
            IF (N .EQ. 1 .AND. (LINE(L)(N:) .NE. ' ') .OR.
     1          (LINE(L)(N:).NE.' ' .AND. LINE(L)(N-1:N-1).EQ.' '))
     2      THEN
               RSTR = LINE(L)(N:)
               CALL CKPARR (RSTR, -1, 1, SPAR(NWANT), NVAL, IER, LOUT)
               IF (IER .NE. 0) IERR = .TRUE.
C              done with these characters
               LINE(L)(N:) = ' '
               NWANT = NWANT - 1
               IF (NWANT .EQ. 0) THEN
C                 done if found 3
                  IF (LINE(L)(1:) .EQ. ' ') LASTLINE = LASTLINE-1
C                 done with this line if blank
                  GO TO 25
               ENDIF
            ENDIF
   20 CONTINUE
C
   25 CONTINUE
C
C     assume reaction reversible
      LREV = .TRUE.
C     assume no real coefficients
      LRNU = .FALSE.
C     assume no #-modified species coefficient
      LYLD = .FALSE.
C
C     find delimeter and set reactant string and product string
      NDLIM = 0
      RSTR = ' '
      PSTR = ' '
C
C     loop over remaining lines for reaction
      DO 50 L = 1, LASTLINE
C
         IREAC = 0
         IPROD = 1
         IF (NDLIM .LE. 0) THEN
C           does this line have the delimiter
            NDLIM = INDEX(LINE(L), '<=>')
            IF (NDLIM .GT. 0) THEN
C              reversible reaction
               LREV = .TRUE.
               IREAC = NDLIM - 1
               IPROD = NDLIM + 3
C
            ELSE
               NDLIM = INDEX(LINE(L), '=>')
               IF (NDLIM .GT. 0) THEN
C                 irreversible reaction
                  LREV = .FALSE.
                  IREAC = NDLIM - 1
                  IPROD = NDLIM + 2
C
               ELSE
                  NDLIM = INDEX(LINE(L), '=')
                  IF (NDLIM .GT. 0) THEN
C                    reversible reaction
                     LREV = .TRUE.
                     IREAC = NDLIM - 1
                     IPROD = NDLIM + 1
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
         IF (NDLIM .EQ. 0) THEN
C           this whole line belongs to RSTR reactant string
            IREAC = CKLSCH(LINE(L))
            IPROD = IREAC + 1
         ENDIF
C
C        add non-blank characters to RSTR reactant string
C        IREAC counts whole line if no delimeter found yet
C        IREAC counts part of this line if delimeter in this line
C        IREAC is 0 otherwise
C
         ITEMP = ' '
         ITEMP = LINE(L)(1:IREAC)
         IR = 0
         DO 45 I = 1, IREAC
            IF (ITEMP(I:I) .EQ. ' ') GO TO 45
            IR = IR + 1
            RSTR(IR:IR) = ITEMP(I:I)
   45    CONTINUE
C
C        add non-blank characters to PSTR product string
C        IPROD starts after delimeter if in this line
C        IPROD is greater than line length if no delimeter found yet
C        IPROD is 1 so need whole line if delimeter found previously
C
         ITEMP = ' '
         ITEMP = LINE(L)(IPROD:CKLSCH(LINE(L)))
         IP = 0
         DO 46 I = 1, CKLSCH(ITEMP)
            IF (ITEMP(I:I) .EQ. ' ') GO TO 46
            IP = IP + 1
            PSTR(IP:IP) = ITEMP(I:I)
   46    CONTINUE
C        go to previous reaction line
   50 CONTINUE
C
      IF (NDLIM .LE. 0) THEN
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1   'Error...no delimeter =, <=>, or => found....'
      ELSE
C        write reaction string
         IF ((IR + IP + 3) .GE. 38) THEN
            WRITE (LOUT, 900) NIISUR, RSTR(1:IR),
     1         (SPAR(N), N=1,NSPAR)
            IF (LREV) THEN
               WRITE (LOUT, '(6X,A,A)') '<=>',PSTR(1:IP)
            ELSE
               WRITE (LOUT, '(6X,A,A)') '=>',PSTR(1:IP)
            ENDIF
         ELSE
            IF (LREV) THEN
               WRITE (LOUT, 910) NIISUR, RSTR(1:IR), '<=>', PSTR(1:IP),
     1                           (SPAR(N), N=1,NSPAR)
            ELSE
               WRITE (LOUT, 910) NIISUR, RSTR(1:IR), '=>', PSTR(1:IP),
     1                           (SPAR(N), N=1,NSPAR)
            ENDIF
         ENDIF
      ENDIF
C
      IF (NWANT .NE. 0) THEN
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1   'Error...problem finding Arrhenius parameters...'
      ELSEIF (SPAR(1) .LE. 0.0) THEN
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1   'Error...pre-exponential must be positive...'
      ENDIF
      KERR = KERR.OR.IERR
C
C     find reactants (J=1), products (J=2), and their coefficients
C
      NREAC = 0
      NPROD = 0
      DO 600 J = 1, 2
C
         IF (J .EQ. 1) THEN
            ISTR = RSTR
         ELSE
            ISTR = PSTR
         ENDIF
C
C        '+'-sign is delimeter between species
         NPLUS = 0
         ILAST = CKLSCH(ISTR)
         DO 100 N = 1, ILAST
            IF (N .EQ. ILAST) THEN
               NPLUS = NPLUS + 1
               IPLUS(NPLUS) = ILAST+1
            ELSEIF
     1      (ISTR(N:N).EQ.'+' .AND. ISTR(N+1:N+1).NE.'+') THEN
C              index of the character after the end of a species
               NPLUS = NPLUS + 1
               IPLUS(NPLUS) = N
            ENDIF
  100    CONTINUE
C
C        loop over the '+'delimited substrings
         DO 200 N = 1, NPLUS
C
            IERR = .FALSE.
C           a reactant or product substring
            IF (N .EQ. 1) THEN
               ISTART = 1
            ELSE
               ISTART = IPLUS(N-1)+1
            ENDIF
C
            ISPEC = ' '
            ISPEC = ISTR(ISTART:IPLUS(N)-1)
C
            ISTART = 0
            DO 120 I = 1, CKLSCH(ISPEC)
C              does this species have #-sign
               IF (ISPEC(I:I) .EQ. '#') THEN
                  CALL SKPCMP(ISPEC(I+1:), KNAME, KKTOT, PNAME,
     1                        NPHASE, KKPHAS, KNUM, NF)
               ELSE
                  CALL SKPCMP (ISPEC(I:), KNAME, KKTOT, PNAME, NPHASE,
     1                         KKPHAS, KNUM, NF)
               ENDIF
C              if species found, go to find a coefficient
               IF (KNUM .GT. 0) THEN
                  ISTART = I
                  GO TO 125
               ENDIF
  120       CONTINUE
C
  125       CONTINUE
C
C           default integer coefficient
            ICOEF = 1
C           default real coefficient
            RCOEF = 1.0
C
            IER = 0
            NVAL = 1
C
C           species name starts at ISTART; if ISTART>1, coefficient
            IF (ISTART .GT. 1) THEN
               IND = INDEX(ISPEC(1:ISTART-1), '.')
               IF (IND .GT. 0) THEN
C                 real coefficient
                  CALL CKPARR (ISPEC(1:ISTART-1), -1, 1, RCOEF, NVAL,
     1                         IER, LOUT)
                  IF (.NOT. LRNU) THEN
C                    initialize this real-coefficient reaction
                     NRNU = NRNU + 1
                     IRNU(NRNU) = NIISUR
                     LRNU = .TRUE.
C                    convert any previous coefficients
                     DO 111 L = 1, MAXSPR
                        RNU(L, NRNU) = NU(L)
                        NU(L) = 0
  111                CONTINUE
                  ENDIF
C
               ELSE
C                 integer coefficient
                  CALL CKPARI (ISPEC(1:ISTART-1), -1, 1, ICOEF, NVAL,
     1                         IER, LOUT)
                  RCOEF = ICOEF
              ENDIF
           ENDIF
C
           IF (IER.NE.0 .OR. NVAL.NE.1) THEN
              IERR = .TRUE.
              WRITE (LOUT, '(6X,A,A)')
     1        'Error finding stoichiometric coefficient...',
     2        ISPEC(1:CKLSCH(ISPEC))
           ENDIF
C
           IF (NF .GT. 1) THEN
              IERR = .TRUE.
              WRITE (LOUT, '(6X,A,A,/a)')
     1        'Error...non-unique species name...',
     2        ISPEC(1:CKLSCH(ISPEC)),
     2        'possible correction is "species/phase/"...'
           ENDIF
C
           IF (KNUM .EQ. 0) THEN
C             this was not a good species name
              IERR = .TRUE.
              WRITE (LOUT, '(6X,A,A)')
     1        'Error...undeclared species...', ISPEC(1:CKLSCH(ISPEC))
           ENDIF
C
           KERR = KERR.OR.IERR
           IF (IERR) GO TO 200
C
C          is this a #-modified species?
           IF (ISPEC(ISTART:ISTART) .EQ. '#' .AND. .NOT. LYLD) THEN
C             this is first species with '#' for this reaction
              LYLD = .TRUE.
              NYLD = NYLD + 1
              IYLD(NYLD) = NIISUR
           ENDIF
C
C          to get here, species=KNUM, coef=ICOEF and/or RCOEF
C          find place to put this species
           NPOS = 0
           IF (J .EQ. 1) THEN
C             check previous reactant species
              DO 221 L = 1, NREAC
                 IF (KNUM .EQ. NUNK(L)) THEN
C                   species match, but # must be treated separately
                    IF (ISPEC(ISTART:ISTART).EQ.'#') THEN
C                      the new species is #-modified (NYLD > 0)
                       IF (KYLD(L,NYLD) .EQ. 1) THEN
C                         the previous species is also #-modified
C                         OK to add to this species coefficients
                          NPOS = L
                       ENDIF
                    ELSE
C                      the new species is not #-modified
                       IF (.NOT. LYLD) THEN
                          NPOS = L
                       ELSE 
                          IF (KYLD(L,NYLD) .NE. 1) THEN
C                            neither is the previous occurence
C                            OK to add to this species coefficients
                             NPOS = L
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
  221         CONTINUE
C
              IF (NPOS .EQ. 0) THEN
C                this is a new reactant
                 IF (NREAC .EQ. MAXSPR/2) THEN
                    IERR = .TRUE.
                    WRITE (LOUT, '(6X,A,A)')
     1              'Error...reactant array full, cannot add...',
     2              ISPEC(1:CKLSCH(ISPEC))
                 ELSE
                    NREAC = NREAC + 1
                    NPOS = NREAC
                 ENDIF
              ENDIF
C
           ELSE
C             check previous product species
              DO 222 L = MAXSPR/2 + 1, MAXSPR/2 + NPROD
                 IF (KNUM .EQ. NUNK(L)) THEN
C                   species match, but # must be treated separately
                    IF (ISPEC(ISTART:ISTART).EQ.'#') THEN
C                      the new species is #-modified (NYLD > 0)
                       IF (KYLD(L,NYLD) .EQ. 1) THEN
C                         the previous species is also #-modified
C                         OK to add to this species coefficients
                          NPOS = L
                       ENDIF
                    ELSE
C                      the new species is not #-modified
                       IF (.NOT. LYLD) THEN
                          NPOS = L
                       ELSE
                          IF (KYLD(L,NYLD) .NE. 1) THEN
C                            neither is the previous occurence
C                            OK to add to this species coefficients
                             NPOS = L
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
  222         CONTINUE
C
              IF (NPOS .EQ. 0) THEN
C                this is a new product
                 IF (NPROD .EQ. MAXSPR/2) THEN
                    IERR = .TRUE.
                    WRITE (LOUT, '(6X,A,A)')
     1              'Error...product array full, cannot add...',
     1              ISPEC(1:CKLSCH(ISPEC))
                 ELSE
                    NPROD = NPROD + 1
                    NPOS = MAXSPR/2 + NPROD
                 ENDIF
              ENDIF
           ENDIF
C
           KERR = KERR.OR.IERR
           IF (.NOT. IERR) THEN
C
              NUNK(NPOS) = KNUM
              IF (LRNU) THEN
                 RNU(NPOS, NRNU) = RNU(NPOS, NRNU) + RCOEF
                 NU(NPOS) = 0
              ELSE
                 NU(NPOS) = NU(NPOS) + ICOEF
              ENDIF
C
              IF (LYLD) THEN
                 IF (ISPEC(ISTART:ISTART) .EQ. '#') THEN
C                   this species uses #-modification
                    KYLD(NPOS,NYLD) = 1
                 ELSE
                    KYLD(NPOS,NYLD) = 0
                 ENDIF
              ENDIF
           ENDIF
C
C          go to next '+'-delimited substring parsing
 200     CONTINUE
C
  600 CONTINUE
C
C     reactants have negative stoichiometric coefficients
      IF (LRNU) THEN
         DO 650 N = 1, NREAC
            RNU(N,NRNU) = -RNU(N,NRNU)
  650    CONTINUE
      ELSE
         DO 700 N = 1, NREAC
            NU(N) = -NU(N)
  700    CONTINUE
      ENDIF
C
C     sum of gas-phase reactant and product coefficients
      DO 750 N = 1, MAXSPR
         IF (NUNK(N) .LE. KKGAS) NUSUMK = NUSUMK + NU(N)
  750 CONTINUE
C
      NSPEC = NREAC + NPROD
C
C     irreversible reaction has negative NSPEC
      IF (.NOT. LREV) NSPEC = -NSPEC
C
      IF (LYLD) THEN
C        #-modified reaction must have exactly one gas-phase ion
C        reactant, we check here in SKREAC rather than in SKAUXL
C        in case user neglects to use YIELD/parameters/ option.
C
         IF (KELECT.LE.0 .OR. KELECT.GT.KKGAS) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...#-modify requires ',
     2      'gas-phase electron in mechanism SPECIES...'
         ENDIF
C
         IF (NSPEC .GT. 0) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A,A)') 'Error...#-modify requires ',
     1      'irreversible reaction...'
         ENDIF
C
         NGI = 0
         DO 925 M = 1, MAXSPR
            NK = NUNK(M)
            IF (NK.NE.0 .AND. NK.NE.KELECT .AND. NK.LE.KKGAS) THEN
               IF (KCHRG(NK).NE.0) THEN
C                  this is a gas-phase non-electron ion
                  NGI = NGI + 1
                  IF (KYLD(M, NYLD) .GT. 0) THEN
C                  this species has a #
                     KERR = .TRUE.
                     ILEN = CKLSCH(KNAME(NK))
                     WRITE (LOUT, '(6X,A,A)')
     1                  'Error...cannot #-modify the gas-phase ion...',
     2                  KNAME(NK)(1:ILEN)
                  ELSE
                     KYION(NYLD) = NK
                  ENDIF
               ENDIF
            ENDIF
  925    CONTINUE
         IF (NGI .LT. 1) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A,A)') 'Error...#-modify requires ',
     1      'one gas-phase non-electron ion in reaction...'
         ENDIF
         IF (NGI .GT. 1) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A,A)') 'Error...#-modify requires ',
     1      'exactly one gas-phase non-electron ion in reaction...'
         ENDIF
C
      ENDIF
C
  900 FORMAT (I4,'. ', A, T53, 1PE8.2, 2X, 0PF5.1, 2X, F9.1)
  910 FORMAT (I4,'. ', A, A, A, T53, 1PE8.2, 2X, 0PF5.1,2X,F9.1)
C
C     end of SUBROUTINE SKREAC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSET
     1   (LRWK, WORK, LIWK, IWORK, LCWK, CWORK, MDIM, KDIM, IDIM,
     2    MXPHSE, MAXSPR, MAXTP, NFIT, NTR, NSPAR, NSCOV, NEDPAR,
     3    NYPAR, MAXORD, KION, NT, TMP, A, KCHRG, KPHSE, AWT, WT, DEN,
     4    KNCF, KCOV, KFIRST, KLAST, KKPHAS, PDEN, NSPEC, NU, NUNK,
     5    RNCF, NREAC, NUSUMK, SPAR, IDUP, ICOV, IKCOV, CPAR, IREV,
     6    RSPAR, ISTK, MOTZ, IRNU, RNU, IORD, KORD, RORD, IBHM, KBHM,
     7    IEDP, KEDP, PEDP, KERR, ITHRM, NONCON, LPDEN, LKDEN, LMOTZ,
     8    KNAME, ENAME, PNAME, IYLD, KYION, KYLD, PYLD, YNCF, EQFAC)
C
C     SKSET sets initial values of arrays and matrices
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
C
C     Integer arrays
      DIMENSION IWORK(LIWK), NT(KDIM), KCHRG(KDIM), KPHSE(KDIM),
     1          KNCF(MDIM,KDIM), KCOV(KDIM), KION(KDIM),
     2          KFIRST(MXPHSE), KLAST(MXPHSE), KKPHAS(MXPHSE),
     2          NSPEC(IDIM), NU(MAXSPR,IDIM), NUNK(MAXSPR,IDIM),
     3          NREAC(IDIM), NUSUMK(IDIM), IDUP(IDIM), ICOV(IDIM),
     5          IKCOV(IDIM), IREV(IDIM), ISTK(IDIM), MOTZ(IDIM),
     6          IRNU(IDIM), IORD(IDIM), KORD(MAXORD,IDIM), IBHM(IDIM),
     7          KBHM(IDIM), IEDP(IDIM), KEDP(IDIM), IYLD(IDIM),
     8          KYION(IDIM), KYLD(MAXSPR,IDIM)
C
C     Real arrays
      DIMENSION WORK(LRWK), TMP(MAXTP,KDIM), A(NFIT,NTR,KDIM),
     1          AWT(KDIM), WT(KDIM), DEN(KDIM), PDEN(MXPHSE),
     2          RNCF(MXPHSE,IDIM), SPAR(NSPAR,IDIM), CPAR(NSCOV,IDIM),
     3          RSPAR(NSPAR,IDIM), RNU(MAXSPR,IDIM), RORD(MAXORD,IDIM),
     4          PEDP(NEDPAR,IDIM), PYLD(NYPAR,IDIM), YNCF(MXPHSE,IDIM),
     5          EQFAC(IDIM)
C
      LOGICAL KERR, ITHRM(KDIM), NONCON, LPDEN(MXPHSE), LKDEN(KDIM),
     1        LMOTZ
      CHARACTER*16 KNAME(KDIM), ENAME(MDIM), CWORK(LCWK), PNAME(MXPHSE)
C
      IF (NMAT .EQ. 1) THEN
C
C        we are processing the first material, so initialize chemkin
C        real, integer, and character work space
C
         KKGAS = 0
         DO 1 L = 1, LRWK
            WORK(L) = 0.0
    1    CONTINUE
         DO 2 L = 1, LIWK
            IWORK(L) = 0
    2    CONTINUE
         DO 3 L = 1, LCWK
            CWORK(L) = ' '
    3    CONTINUE
C
C        initialize atomic weights and element names
C
         DO 5 M = 1, MDIM
            AWT(M) = 0.0
            ENAME(M) = ' '
    5    CONTINUE
      ENDIF
C
C     initialize space associated with surface and bulk species
C     (all of the species past kkgas)
C
      DO 100 K = KKGAS+1, KDIM
C
C        species names
C
         KNAME(K) = ' '
C
C        number of temperatures used to define polynomial fits,
C        i.e., TLOW, TMID, THIGH for the case of MAXTP=3
C
         NT(K) = MAXTP
C
C        initialize temperature bounds to a negative
C        number as a flag to indicate that they haven't been set
C
         DO 10 M = 1, MAXTP
            TMP(M,K) = -1.0
   10    CONTINUE
C
C        "A(m,n,k)" contains thermo polynomial coefficients;
C        first dimension runs over the NFIT number of coefficients,
C        second dimension runs over the temperature ranges,
C        third dimension is for the species index;
C        initialize values to zero
C
         DO 15 M = 1, NFIT
            DO 15 N = 1, NTR
               A(M, N, K) = 0.0
   15    CONTINUE
C
C        species electronic charge
C
         KCHRG(K) = 0
C
C        phase in which the species resides
C
         KPHSE(K) = 0
C
C        species molecular weight
C
         WT(K) = 0.0
C
         KION(K) = 0
C
C        species density
C
         DEN(K) = -1.0
C
C        elemental composition; the first dimension runs over
C        elements, the second is the species index.
C
         DO 20 M = 1, MDIM
            KNCF(M,K) = 0
   20    CONTINUE
C
C        the number of surface sites that this species covers
C
         KCOV(K) = 1
C
C        thermo fits have not yet been read for this species
C
         ITHRM(K) = .FALSE.
C
C        bulk species will turn this flag to true when they are read
C
         LKDEN(K) = .FALSE.
  100 CONTINUE
C
C
C     initialize information having to do with phases
C
      DO 200 M = 1, MXPHSE
C        PHASE NAME
         PNAME(M) = ' '
C
C        arrays containing species indices of the first and last
C        species in a phase
C
         KFIRST(M) = 0
         KLAST(M) = 0
C
C        total count of species in each phase
C
         KKPHAS(M) = 0
C
C        phase density; primarily relevant to surface phases,
C        where it holds site density
C
         PDEN(M) = -1.0
C
C        flag will be set to true when a site density is read for
C        a given phase
C
         LPDEN(M) = .FALSE.
  200 CONTINUE
C
C     initialize information concerning reactions
C
      DO 300 I = 1, IDIM
C
C        total number of reactants plus products in a reaction
C        (is set to a negative number if reaction is irreversible)
C
         NREAC(I) = 0
C
C        MAXSPR is maximum number of species allowed per reaction
C
         DO 205 M = 1, MAXSPR
C
C           stoichiometric coefficient for the Mth reactant in the
C           Ith reaction
C
            NU(M,I) = 0
C
C           species number corresponding to NU(M,I)
C
            NUNK(M,I) = 0
  205    CONTINUE
C
         DO 210 M = 1, MXPHSE
C
C           net change in phase M due to reaction I
C
            RNCF(M,I) = 0.0
  210    CONTINUE
C
C        total number of reactants in reaction I
C
         NREAC(I) = 0
C
C        sum of stoichiometric coefficients of reactants and products
C
         NUSUMK(I) = 0
C
C        NSPAR rate coefficients for reaction I
C
         DO 215 N = 1, NSPAR
            SPAR(N,I) = 0.0
  215    CONTINUE
C
C        reaction I is duplicate reaction
C
         IDUP(I) = 0
C
C        coverage parameters were specified for reaction I
C
         ICOV(I) = 0
C
C        the species number for which reaction I has a coverage
C        dependence
C
         IKCOV(I) = 0
C
C        NSCOV coverage coefficients for reaction I
C
         DO 220 N = 1, NSCOV
            CPAR(N,I) = 0.0
  220    CONTINUE
C
C        reaction I is specified to be irreversible
         IREV(I) = 0
C
C        NSPAR reverse rate coefficients for reaction I
         DO 225 N = 1, NSPAR
            RSPAR(N,I) = 0.0
  225    CONTINUE
C
C        sticking coefficients were specified for reaction I
         ISTK(I) = 0
C
C        Motz-Wise correction ON is the default
         MOTZ(I) = 1
C
C        real stoichiometric coefficients are used for reaction I
         IRNU(I) = 0
C
         DO 230 M = 1, MAXSPR
C           real stoichiometric coefficient for Mth species in
C           reaction I
            RNU(M,I) = 0.0
  230    CONTINUE
C
C        a non-standard reaction order dependence was specified for
C        reaction I
         IORD(I) = 0
         DO 235 M = 1, MAXORD
C
C           non-standard reaction order species index
            KORD(M,I) = 0
C
C           reaction order for that species
            RORD(M,I) = 0.0
  235    CONTINUE
C
C        a Bohm keyword was read for reaction I
         IBHM(I) = 0
C
C        Bohm reaction ionic species index
         KBHM(I) = 0
C
         DO 240 N = 1, NEDPAR
C           coefficients to specify the energy dependence
C           in reactions for which ENRGDEP keyword was given
            PEDP(N,I) = 0.0
  240    CONTINUE
C
C        ENRGDEP keyword was read for reaction I
         IEDP(I) = 0
C
C        species index of the ionic species in a reaction for which
C        the ENRGDEP keyword was read
         KEDP(I) = 0
C
         IYLD(I) = 0
         KYION(I) = 0
         DO 290 N = 1, MAXSPR
            KYLD(N,I) = 0
  290    CONTINUE
         DO 295 N = 1, NYPAR
            PYLD(N,I) = 0.0
  295    CONTINUE
         DO 296 N = 1, MXPHSE
            YNCF(N,I) = 0.0
  296    CONTINUE
C
         EQFAC(I) = 1.0
  300 CONTINUE
C
C     no errors encountered so far
      KERR = .FALSE.
C
C     use default definition of a sticking coefficient, which
C     applies Motz-Wise correction
      LMOTZ = .TRUE.
C
C     surface sites are normally assumed to be conserved, so the
C     NONCON flag is initialized to FALSE
      NONCON = .FALSE.
C
C     number of phases, surface phases, bulk phases, surface species,
C     bulk species, total number of species, surface reactions:
C
      NPHASE = 0
      NNSURF = 0
      NNBULK = 0
      KKSURF = 0
      KKBULK = 0
      KKTOT = 0
      NIISUR = 0
C
C     number of reactions which include a coverage dependence
      NCOV = 0
C
C     number of reactions that are irreversible
      NREV = 0
C
C     number of reactions that have sticking coefficients
      NSTK = 0
C
C     number of reactions which do not conserve sites
      NCON = 0
C
C     number of reactions that have real stoichiometric coefficients
      NRNU = 0
C
C     number of reactions that have non-standard reaction orders
      NORD = 0
C
C     number of reactions that have Bohm correction
      NBHM = 0
C
C     number of reactions for which an ion-energy dependence was given
      NEDP = 0
C
      NYLD = 0
C
C     end of SUBROUTINE SKSET
      RETURN
      END
C                                                                      C
      SUBROUTINE SKSIZE (NSPAR, MAXSPR, MAXTP, NTR, NFIT, NSCOV,
     1                   NEDPAR, NYPAR, MAXORD)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE

C     Size of integer workspace for surface chemkin:
      LENISK = 88
C              (ISKWRK pointers)
     1       + NPHASE*3
C              (start, end, count of phase species)
     1       + KKTOT*(4 + NELEM)
C              (phase, charge, coverage, composition of species)
     2       + NIISUR*(3 + 2*MAXSPR)
C              (nreac, nspec, nsum, nunk, nu for reactions)
     3       + NCOV*2
C              (reaction indices, coverage species)
     4       + NREV
C              (reaction indices)
     5       + NSTK*2
C              (reaction indices, Motz flags)
     6       + NBHM*2
C              (reaction indices, Bohm species indices)
     7       + NRNU
C              (reaction indices)
     8       + NORD*(1 + MAXORD)
C              (reaction indices, species indices)
     9       + NEDP*2
C              (reaction indices, temperature flags)
     *       + NYLD*(2 + MAXSPR)
C              (reaction indices, ion species indices, yield indices)
     1       + NKION
C              (ion species)
     2       + KKGAS
C              (temperature flags)
     3       + 1
C              (MORE)
C
C
      LENRSK = 4
C             (constants)
     1       + NPHASE
C             (phase densities)
     2       + KKTOT*(MAXTP + NFIT*(MAXTP-1) + 1)
C             (fit temps, thermodynamic coefficients, species dens.)
     3       + NELEM
C             (atomic weights)
     4       + KKTOT
C             (molecular weights)
     5       + NIISUR*(NSPAR+1)
C             (Arrhenius coefficients)
     6       + NCOV*NSCOV
C             (coverage parameters)
     7       + NREV*(NSPAR+1)
C             (reverse Arrhenius coefficients)
     8       + NIISUR
C             (equilibrium scalar)
     9       + NRNU*MAXSPR
C             (real stoichiometric coefficients)
     *       + NIISUR*NPHASE
C             (reaction phase/site balances)
     1       + NYLD*(NYPAR + NPHASE)
C             (yield reaction phase/site balances)
     2       + NORD*MAXORD
C             (changed-order values)
     3       + NEDP*NEDPAR
C             (ion-energy reaction parameters
     4       + KKGAS
C             (species ion energy array)
     5       + NIISUR*2
C             (forward and reverse temp-dependent reaction rates)
     6       + KKTOT*2
C             (species scratch space)
     7       + NPHASE
C             (phase scratch space)
     8       + NIISUR*3
C             (reaction scratch space)
      LENCSK = NELEM
C              (element names)
     1       + KKTOT
C              (species names)
     2       + NPHASE
C              (phase names)
     3       + 1
C              (material name)
C
C     end of SUBROUTINE SKSIZE
      RETURN
      END
C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSPEC (KDIM, K, KNAME, KKTOT, ITHRM, MAXTP, TMP, NT,
     1                   TMID, LOUT, KERR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NT(KDIM), TMP(MAXTP,KDIM), IPLUS(10)
      CHARACTER KNAME(KDIM)*16, INUM(10)*1
      LOGICAL ITHRM(*), KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      DATA INUM/'0','1','2','3','4','5','6','7','8','9'/
C
C     each species must have thermodynamic data
C
      IF (.NOT.ITHRM(K)) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...no thermodynamic properties for species ',KNAME(K)
      ENDIF
C
C     thermodynamic temperatures must be in order
      DO 111 N = 2, NT(K)
         IF (TMP(N,K) .LT. TMP(N-1,K)) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A)')
     1      'Error...thermodynamics temperatures must be in order...'
         ENDIF
  111 CONTINUE
C
C     a species cannot start with a number
      CALL CKNCMP (KNAME(K)(1:1), INUM, 10, I, NF)
      IF (I .GT. 0) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1   'Error...species name starts with a number'
      ENDIF
C
C     is there another species name after a +
      NPLUS = 0
      DO 30 N = 1, CKLSCH(KNAME(K))
         IF (KNAME(K)(N:N) .EQ. '+') THEN
            NPLUS = NPLUS + 1
            IPLUS(NPLUS) = N
         ENDIF
   30 CONTINUE
      DO 40 N = 1, NPLUS
         I1 = IPLUS(N)
         IF (I1 .EQ. 1) THEN
            WRITE (LOUT, '(6X,A)')
     1      'Error...species name starts with a plus'
            KERR = .TRUE.
         ELSE
C
C           is there another species name after a +
            I1 = I1 + 1
            IF (N .LT. NPLUS) THEN
               DO 35 L = N+1, NPLUS
                  I2 = IPLUS(L)
                  IF (I2 .GT. I1) THEN
                     CALL CKNCMP (KNAME(K)(I1:I2-1),KNAME,
     1                            KKTOT, KNUM, NF)
                     IF (KNUM .GT. 0) THEN
                        WRITE (LOUT, 230)
                        KERR = .TRUE.
                     ENDIF
                  ENDIF
   35          CONTINUE
            ENDIF
            I2 = CKLSCH(KNAME(K))
            IF (I2 .GE. I1) THEN
               CALL CKNCMP (KNAME(K)(I1:I2), KNAME, KKTOT, KNUM, NF)
               IF (KNUM .GT. 0) THEN
                  WRITE (LOUT, 230)
                  KERR = .TRUE.
               ENDIF
            ENDIF
         ENDIF
   40 CONTINUE
C
  230 FORMAT (6X,'Error...illegal + in species name')
C
C     end of SUBROUTINE SKSPEC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKBGIN (LIWK, LRWK, LCWK, LINCK, LOUT, MDIM, KDIM,
     1                   MAXTP, NFIT, NTR, IWORK, WORK, CWORK, AWT,
     2                   ENAME, KNAME, WT, ITHRM, KCHRG, KPHSE, NT,
     3                   TMP, KNCF, KION, A, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKBGIN
C  intializes a CHEMKIN linkfile to get gas-phase chemistry.
C
C  Arguments:
C  LIWK     - Integer scalar, size of integer work array.
C  LRWK     - Integer scalar, size of real work array.
C  LCWK     - Integer scalar, size of character-string work array.
C  LINCK   - Integer scalar, CHEMKIN input linkfile unit number.
C  LOUT     - Integer scalar, formatted output file unit number.
C  MDIM     - Integer scalar, maximum number of elements.
C  KDIM     - Integer scalar, maximum number of species.
C  MAXTP    - Integer scalar, maximum number of fit temperatures.
C  NFIT     - Integer scalar, maximum number of fit coefficients.
C  NTR      - Integer scalar, maximum number of temperature ranges.
C  IWORK(*) - Integer array, CHEMKIN workspace.
C  WORK(*)  - Real array, CHEMKIN workspace.
C  CWORK(*) - Character-string array, CHEMKIN workspace.
C  AWT(*)   - Real array, atomic weights.
C  ENAME(*) - Character-string array, element names.
C  KNAME(*) - Character-string array, species names.
C  WT(*)    - Real array, molecular weights.
C  ITHRM(*) - Logical array, species thermodynamic data flag.
C  KCHRG(*) - Integer array, species electronic charges.
C  KPHSE(*) - Integer array, species physical state.
C  NT(*)    - Integer array, species fit temperature counts.
C  TMP(*,*) - Real matrix, species fit temperatures.
C  KNCF(*,*)- Integer matrix, species elemental compositions.
C  KION(*)  - Integer array, ion species indices.
C  A(*,*,*) - Real three-dimensional array, species thermodymic data
C             polynomial coefficients.
C  KERR     - Logical error flag
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
C
C     Integer arrays
      DIMENSION IWORK(LIWK), KNCF(MDIM,KDIM), KCHRG(KDIM), NT(KDIM),
     1          KPHSE(KDIM), KION(KDIM)
C     Real arrays
      DIMENSION WORK(LRWK), AWT(MDIM), WT(KDIM), TMP(MAXTP,KDIM),
     1          A(NFIT, NTR, KDIM)
C     Character arrays
      CHARACTER*16 KNAME(KDIM), ENAME(MDIM), CWORK(LCWK)
      LOGICAL ITHRM(KDIM), IERR, KERR
C
C     initialize the gas-phase CHEMKIN library
      CALL CKINIT (LIWK, LRWK, LCWK, LINCK, LOUT, IWORK, WORK, CWORK,
     1             IFLAG)
C
C     quit if an error occurred in CKINIT
      IF (IFLAG .GT. 0) THEN
         KERR = .TRUE.
         RETURN
      ENDIF
C     get gas-phase indices
      CALL CKINDX (IWORK, WORK, NELEM, KKGAS, IGAS, NCKFIT)
      CALL CKMXTP  (IWORK, MXKTMP)
      IF (MXKTMP .GT. MAXTP) THEN
         WRITE (LOUT, *) ' MAXTP too small for gas-phase ',
     1                   ' fit temperature arrays...must be at least ',
     2                   MXKTMP
         KERR = .TRUE.
      ENDIF
C
C     is enough room provided for the number of elements?
      IF (MDIM .LT. NELEM) THEN
C        no; issue an error message
         WRITE (LOUT, *) ' Element dimension too small,',
     1                   ' MDIM must be at least...',NELEM
         KERR = .TRUE.
      ENDIF
C
C     is enough space provided for the species?
      IF (KDIM .LT. KKGAS) THEN
C        no; issue an error message
         WRITE (LOUT, *) ' Species dimension too small,',
     1                   ' KDIM must be at least...', KK
         KERR = .TRUE.
      ENDIF
C
C     cannot process element or species data
      IF (KERR) RETURN
C
C     get names of elements and species
      CALL CKSYME (CWORK, LOUT, ENAME, IERR)
      KERR = KERR.OR.IERR
      CALL CKSYMS (CWORK, LOUT, KNAME, IERR)
      KERR = KERR.OR.IERR
C
      IF (KERR) RETURN
C
C     get element atomic weights and species molecular weights
      CALL CKAWT  (IWORK, WORK, AWT)
      CALL CKWT   (IWORK, WORK, WT)
C
C     get electronic charge of species
      CALL CKCHRG (IWORK, WORK, KCHRG)
C
C     find out what phase the species is declared to be in:
C        KPHASE(K)=-1, species K is solid
C        KPHASE(K)= 0, species K is gaseous
C        KPHASE(K)=+1, species K is liquid
      CALL CKPHAZ (IWORK, WORK, KPHSE)
C
C     get elemental composition of the species
      CALL CKNCF  (MDIM, IWORK, WORK, KNCF)
C
C     get thermodynamic coefficients for the species
      CALL CKATHM (NFIT, NTR, IWORK, WORK, MAXTP, NT, TMP, A)
      DO 101 K = 1, KKGAS
C        thermo data has been read for the gas-phase species
         ITHRM(K) = .TRUE.
  101 CONTINUE
C
      MELECT = 0
C     MELECT is the electron element, if present
      KELECT = 0
C     KELECT is the electron species, if present
      NKION = 0
C     NKION is the count of ionic species, if present

      DO 200 M = 1, NELEM
         IF (ENAME(M).EQ.'e'.OR. ENAME(M).EQ.'E') MELECT = M
  200 CONTINUE
C
      IF (MELECT .EQ. 0) RETURN
C     KELECT is the electron species, if present
      DO 215 K = 1, KKGAS
         CALL SKKEL (MDIM, MELECT, NELEM, KNCF(1,K), KEL)
         IF (KEL .LE. 0) GO TO 215
         IF (KELECT .GT. 0) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...duplicate electron-only species...',
     2      KNAME(K)
         ELSE
            KELECT = K
         ENDIF
  215    CONTINUE
C
      IF (KELECT .GT. 0 .AND. WT(KELECT) .GT. 1.0E-3) THEN
C        there must be some problem, because the electron's mass
C        should be very tiny; issue an error message
         WRITE (LOUT, *)' THERE IS A PROBLEM IN SKBGIN:'
         WRITE (LOUT, *) ' SPECIES NUMBER ', KELECT
         WRITE (LOUT, *) ' SPECIES NAME: ', KNAME(KELECT)
         WRITE (LOUT, *) ' WAS IDENTIFIED AS THE ELECTRON.'
         WRITE (LOUT, *) ' ITS MASS IS: ', WT(KELECT),
     1                   ', WHICH SEEMS TOO LARGE.'
         KERR = .TRUE.
      ENDIF
C
C     look for ionic species
      DO 120 K = 1, KKGAS
         IF (KCHRG(K) .NE. 0.0) THEN
C           have found an ion; increment counter and
C           save species index in array KION
            NKION = NKION + 1
            KION(NKION) = K
         ENDIF
  120 CONTINUE
C
C     end of SUBROUTINE SKBGIN
      IF (KERR) RETURN
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSURF (SUB, NSUB, NDIM, STRAY, IRAY, NN, KCHECK,
     1                   KSTART, LPDEN, PDEN, KERR, LOUT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSURF
C  processes an array of character strings and adds them to
C  either a character string array STRAY or an integer array IRAY.
C
C  Arguments:
C  SUB(*)     - Character string array.
C  NSUB       - Integer scalar, the count of strings to be processed.
C  NDIM       - Integer scalar, size of STRAY and IRAY arrays.
C  NN         - Integer scalar, count of items in STRAY, IRAY.
C  KCHECK     -
C  KSTART     -
C  LPDEN      - Logical
C  PDEN       - Real scalar,
C  KERR       - Logical error flag.
C  LOUT       - Integer scalar, formatted output file unit number.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION IRAY(*), IPAR(1)
      CHARACTER*(*) SUB(*), STRAY(*)
      CHARACTER ISTR*80, CKCHUP*4
      LOGICAL KERR, IERR, LPDEN
      INTEGER CKLSCH
      EXTERNAL CKLSCH, CKCHUP
C
      ILEN = LEN(STRAY(1))
C
      DO 200 N = 1, NSUB
         IF (CKCHUP(SUB(N), 3) .EQ. 'END') GO TO 1000
         ILAST = CKLSCH(SUB(N))
         IF (CKCHUP(SUB(N), 4) .EQ. 'SDEN') THEN
            I1 = INDEX(SUB(N),'/')
            IF (I1 .GT. 0) THEN
               I2 = I1 + INDEX(SUB(N)(I1+1:),'/')
               IF (I2 .GT. I1) THEN
                  ISTR = ' '
                  ISTR = SUB(N)(I1+1:I2-1)
                  CALL CKXNUM (ISTR, 1, LOUT, NVAL, PDEN, IERR)
               ENDIF
            ENDIF
            IF (IERR) WRITE (LOUT, '(6X,A,A)')
     1      'Error in density declaration...',SUB(N)(1:ILAST)
            KERR = KERR.OR.IERR
            IF (LPDEN) THEN
               WRITE (LOUT, '(6X,A,A)')
     1         'Error...multiple density declaration...',
     2         SUB(N)(1:ILAST)
            ELSE
               LPDEN = .TRUE.
            ENDIF
         ELSE
            ISTR = ' '
            I1 = INDEX(SUB(N),'/')
            IF (I1 .EQ. 1) THEN
               KERR = .TRUE.
               WRITE (LOUT, '(6X,A,A)')
     1         'Error...misplaced value...',SUB(N)(1:ILAST)
            ELSE
               IF (I1 .LE. 0) THEN
                  ISTR = SUB(N)
               ELSE
                  ISTR = SUB(N)(1:I1-1)
               ENDIF
               CALL CKNCMP (ISTR, STRAY, KCHECK, KNUM1, NF)
               IF (NN .GT. KSTART) THEN
                  CALL CKNCMP (ISTR, STRAY(KSTART), NN-KSTART, KNUM2,
     1                         NF)
               ELSE
                  KNUM2 = 0
               ENDIF
C
               IF (KNUM1 .GT. 0) THEN
                  WRITE (LOUT, '(6X,A,A,A)')
     1            'Error...site species name duplicates ',
     2            'gas species name...', SUB(N)(1:ILAST)
                  KERR = .TRUE.
               ELSEIF (KNUM2 .GT. 0) THEN
                  WRITE (LOUT, '(6X,A,A)')
     1            'Warning...duplicate site species name ignored...',
     2            SUB(N)(1:ILAST)
               ELSE
                  IF (NN .LT. NDIM) THEN
                     IF (ISTR(ILEN+1:) .NE. ' ') THEN
                        WRITE (LOUT, '(6X,A,A)')
     1                  'Error...site species name too long...',
     2                  SUB(N)(1:ILAST)
                        KERR = .TRUE.
                     ELSE
                        NN = NN + 1
                        STRAY(NN) = ' '
                        STRAY(NN) = ISTR(1:ILEN)
                        IF (I1 .GT. 0) THEN
                           I2 = I1 + INDEX(SUB(N)(I1+1:),'/')
                           ISTR = ' '
                           ISTR = SUB(N)(I1+1:I2-1)
                           CALL CKPARI (ISTR,1,1,IPAR,NVAL,IER,LOUT)
                           KERR = KERR .OR. (IER.NE.0)
                           IRAY(NN) = IPAR(1)
                        ENDIF
                     ENDIF
                  ELSE
                     WRITE (LOUT, '(6X,A,A)')
     1               'Error...species array size too small for...',
     2               SUB(N)(1:ILAST)
                     KERR = .TRUE.
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  200 CONTINUE
C
C     end of SUBROUTINE SKSURF
 1000 RETURN
      END
C                                                                      C
      SUBROUTINE SKTHRM (LUNIT, LTHRM, THERMO, MDIM, KDIM, ENAME, AWT,
     1                   KNAME, KNCF, KPHSE, KCHRG, WT, MAXTP, NT, NTR,
     2                   T, NFIT, A, ITHRM, KERR, LOUT, ISTR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKTHRM
C  processes formatted thermodynamic properties inputs of species;
C  format -
C  Line 1: species name, optional comments, elemental composition,
C          phase, T(low), T(high), T(mid), additional elemental
C          composition, line number (col. 80);
C          format(A10,A14,4(A2,I3),A1,E10.0,E10.0,E8.0,(A2,I3),I1)
C  Line 2: coefficients a(1--5) for upper temperature range,
C          line number (col. 80); format (5(e15.0), I1)
C  Line 3: coefficients a(6--7) for upper temperature range,
C          coefficients a(1--3) for lower temperature range,
C          line number (col. 80); format (5(e15.0), I1)
C  Line 4: coefficients a(4--7) for lower temperature range,
C          line number (col. 80); format (4(e15.0), I1)
C
C  Arguments:
C  LUNIT     - Integer scalar, formatted input file unit number.
C  LTHRM     - Integer scalar, formatted input file unit number.
C  THERMO    - Logical, thermodynamic processing status flag.
C  MDIM      - Integer scalar, maximum number of elements.
C  ENAME(*)  - Character-string array, element names.
C  AWT(*)    - Real array, element atomic weights.
C  KNAME(*)  - Character-string array, species names.
C  KNCF(*,*) - Integer matrix, elemental composition of species.
C  KPHSE(*)  - Integer array, species physical states.
C  WT(*)     - Real array, species molecular weights.
C  MAXTP     - Integer scalar, maximum number of species fit temps.
C  NT(*)     - Integer array, count of species fit temperatures.
C  MAXTR     - Integer scalar, number of fit temperature ranges.
C  T(*,*)    - Real matrix, species fit temperatures.
C  NFIT      - Integer scalar, number of species fit coefficients.
C  A(*,*,*)  - Real three-dimensional array, species thermodynamic
C              polynomial coefficients.
C  ITHRM(*)  - Logical array, species thermodynamics status flag.
C  KERR      - Logical, error flag.
C  LOUT      - Integer scalar, formatted output file unit number.
C  ISTR      - Character string, the final string input to SKTHRM.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
C     Integer arrays
      DIMENSION NT(KDIM), KPHSE(KDIM), KNCF(MDIM,KDIM), KCHRG(KDIM)
C     Real arrays
      DIMENSION WT(KDIM), T(MAXTP,KDIM), A(NFIT, NTR, KDIM),
     1          AWT(MDIM), VALUE(10)
      CHARACTER ENAME(MDIM)*16, KNAME(MDIM)*16, ELEM*2,
     1          LINE(15)*80, ISTR*(*), SUB(80)*80, CKCHUP*4, KEY*4
      LOGICAL THERMO, KERR, ITHRM(KDIM), SETTMP
      INTEGER CKSLEN
      EXTERNAL CKCHUP, CKSLEN

C     local temperature array
      PARAMETER (MAXTMP=10)
      DIMENSION TFIT(MAXTMP)
C
      IF (NELEM .LE. 0) THEN
         WRITE (LOUT, '(6X,A)')
     1   'Error...cannot use THERM until ELEMents have been input...'
         KERR = .TRUE.
      ENDIF
      IF (KKTOT-KKGAS .LE. 0) THEN
         WRITE (LOUT, '(6X,A)')
     1   'Error...cannot use THERM until SPECies have been input...'
         KERR = .TRUE.
      ENDIF
C      IF (KERR) RETURN
C
      IF (THERMO) THEN
         REWIND LTHRM
         LTEMP = LTHRM
      ELSE
C        THERMO ALL option eliminates need for CHEMKIN thermo data
         LTEMP = LUNIT
      ENDIF
C
  300 CONTINUE
C     stepping through THERMO data
      ISTR = ' '
      READ (LTEMP,'(A)',ERR=22222, END=400) ISTR
      IF (CKSLEN(ISTR) .LE. 0) GO TO 300
C
C     CHEMKIN data has the THERM line to contend with
      IND = MAX(INDEX(ISTR,'THERM'), INDEX(ISTR,'therm'))
      IF (IND .GT. 0) GO TO 300
C
C     first non-blank line after THERM has NTFIT fit temperatures
C     (LTEMP here could be either user's input or CHEMKIN data file);
C     NTFIT is the number of temperatures (nominally 3) on the line,
C     must be <= MAXTP set by calling program, and in ascending order.
C     The species data which follows must then have 7 x NTFIT
C     coefficients, up to 5 per line, from HIGHEST to LOWEST
C     temperature ranges.
C     (To change the requirements temporarily for a specific species,
C      add a line TEMPS value(1) ... value(n), after Line 1.)
C
C     initial array of fit temperatures
      DO 5 N = 1, MAXTMP
         TFIT(N) = -1.0
    5 CONTINUE
      CALL CKDTAB(ISTR)
C
C     fit temperatures to be used as defaults
      CALL CKPARR (ISTR, -1, -MAXTMP, TFIT, NTFIT, IER, LOUT)
      IF (NTFIT .GT. MAXTP) THEN
         WRITE (LOUT, '(6X,A,I3)')
     1   'Error...must increase MAXTP parameter to at least ',NTFIT
         KERR = .TRUE.
         NTFIT = MAXTP
      ENDIF
      IF (IER .NE. 0) THEN
         WRITE (LOUT, '(6X,A,A)')
     1   'Error reading default temperatures...',ISTR
         KERR = .TRUE.
      ENDIF
C
C     temperatures given must be in ascending order
      DO 10 N = 2, NTFIT
         IF (TFIT(N) .LE. TFIT(N-1)) THEN
            WRITE (LOUT, '(6X,A)')
     1      'Error...fit temperatures are not in ascending order'
            KERR = .TRUE.
         ENDIF
   10 CONTINUE
C
C     Species thermodynamic data from here on,
C     LUNIT could be either user's input or CHEMKIN data file.
C
   25 CONTINUE
C     initialize species data storage
      NLINES = 0
      KNUM = 0
      SETTMP = .FALSE.
C
   50 CONTINUE
      ISTR = ' '
      READ (LUNIT,'(A)', ERR=22222, END=400) ISTR
C     skip blank lines
      CALL CKDTAB (ISTR)
      ILEN = CKSLEN(ISTR)
      IF (ILEN .LE. 0) GO TO 50
C
      CALL SKISUB (ISTR(1:ILEN), SUB, NSUB)
      KEY = ' '
      KEY = CKCHUP(SUB(1), 4)
C
      IND = MAX (INDEX(KEY,'END'), INDEX(KEY,'REAC'),
     1           INDEX(KEY,'SITE'),INDEX(KEY,'BULK'))
      IF (IND .GT. 0) THEN
C        end of thermo data
         IF (KNUM .EQ. 0) RETURN
C        still need to process previous species, KNUM
         BACKSPACE (LUNIT)
C        500 is the main species processing section
         GO TO 500
      ENDIF
C
      IF (ILEN.GE.80 .AND. ISTR(80:80) .EQ. '1') THEN
C        '1' in col. 80 starts a set of species data
         IF (KNUM .GT. 0) THEN
C           have previously-stored species data to process
C           before proceeding to check this species
            BACKSPACE (LUNIT)
            GO TO 500
         ENDIF
C        new species data; compare to those in mechanism
         CALL CKNCMP (SUB(1), KNAME, KKTOT, KNUM, KTIMES)
C        this species is not used
         IF (KNUM .LE. 0) GO TO 50
C        this species already has thermo data
         IF (ITHRM(KNUM)) GO TO 25
C        initialize data for the new species
         NLINES = 1
         LINE(NLINES) = ISTR(1:ILEN)
C        look for more input for this species
         GO TO 50
      ENDIF
C
      IF (KNUM .GT. 0) THEN
C        continue storing data for this species
         NLINES = NLINES + 1
         LINE(NLINES) = ISTR(1:ILEN)
         IF (KEY .EQ. 'TEMP') SETTMP = .TRUE.
C        look for more input for this species
      ENDIF
C     next input
      GO TO 50
C
  500 CONTINUE
C
C     Line 1 Input:
C     a) get/store elemental composition;
C     b) use that to calculate/store molecule mass and charge,
C     c) get/store integer representation of molecule phase
C     d) get/store T(1), T(3), T(2) for this particular species,
C        (used only if present and if NTFIT=3)
C
      ICOL = 25
      DO 110 N = 1, 5
         ELEM = ' '
         ELEM = LINE(1)(ICOL:ICOL+1)
         CALL CKPARR (LINE(1)(ICOL+2:ICOL+4), -1, 1, VALUE, NVAL,
     1                IER, LOUT)
         IF (ELEM.NE.' ' .AND. NVAL.GT.0) THEN
            IELEM = INT (VALUE(1))
            IF (IELEM .NE. 0) THEN
               CALL CKCOMP (ELEM, ENAME, NELEM, M)
               IF (M .GT. 0) THEN
C                 composition
                  KNCF(M,KNUM) = KNCF(M,KNUM) + IELEM
C                 molecular weight
                  WT(KNUM) = WT(KNUM) + AWT(M) * VALUE(1)
C                 electronic charge
                  IF (M .EQ. MELECT) KCHRG(KNUM)=KCHRG(KNUM)-VALUE(1)
               ELSE
                  WRITE (LOUT, '(6X,A,A,A)')
     1            'Error...element ', ELEM, 'not declared for species ',
     2            KNAME(KNUM)(1:10)
                  KERR = .TRUE.
               ENDIF
            ENDIF
         ENDIF
         IF (N .LT. 4) THEN
            ICOL = ICOL + 5
         ELSE
            ICOL = 74
         ENDIF
  110 CONTINUE
      IF (CKCHUP(LINE(1)(45:),1) .EQ. 'L') KPHSE(KNUM)=1
      IF (CKCHUP(LINE(1)(45:),1) .EQ. 'S') KPHSE(KNUM)=-1
C
      IF (SETTMP) THEN
C        line 2 has TEMPS followed by substitute fit temperatures
         CALL SKISUB (LINE(2), SUB, NSUB)
         IF (NSUB-1 .GT. MAXTP) THEN
            WRITE (LOUT, '(6X,A,I3,A,A)')
     1      'Error...MAXTP must be increased at least to ',NSUB-1,
     2      'for species ', KNAME(KNUM)
            KERR = .TRUE.
C           can only use MAXTP data for this species
            SETTMP = .FALSE.
         ELSE
            NT(KNUM) = NSUB - 1
            DO 116 N = 2, NSUB
               CALL CKPARR (SUB(N), -1, 1, T(N-1,KNUM), NVAL, IER, LOUT)
  116       CONTINUE
         ENDIF
         LINE(2) = ' '
      ENDIF
C
      IF (.NOT. SETTMP) THEN
         NT(KNUM) = NTFIT
         DO 115 N = 1, NTFIT
            T(N,KNUM) = TFIT(N)
  115    CONTINUE
      ENDIF
C
      IF (NT(KNUM).LE.3 .AND. LINE(1)(46:73).NE.' ') THEN
C        Line 1 fields OK for standard 2-range case, or 1 range
         IF (LINE(1)(46:55) .NE. ' ') CALL CKPARR
     1      (LINE(1)(46:55), 0, 1, T(1,KNUM), NVAL, IER, LOUT)
         IF (LINE(1)(66:73) .NE. ' ') CALL CKPARR
     1      (LINE(1)(66:73), 0, 1, T(2,KNUM), NVAL, IER, LOUT)
         IF (LINE(1)(56:65) .NE. ' ') CALL CKPARR
     1      (LINE(1)(56:65), 0, 1, T(3,KNUM), NVAL, IER, LOUT)
      ENDIF
C
C     Lines 2..NLINES:
C     15-column polynomial coefficients, HIGHEST to LOWEST temperatures
C
      NA = 0
      NRANGE = NT(KNUM) - 1
      DO 80 N = 2, NLINES
         ICOL = 1
         DO 75 L = 1, 5
            IF (LINE(N)(ICOL:ICOL+14) .NE. ' ') THEN
               NA = NA + 1
               IF (NA .EQ. 8) THEN
C                 back up to previous temperature range 7-polynomial set
                  NA = 1
                  NRANGE = NRANGE - 1
                  IF (NRANGE .EQ. 0) GO TO 90
               ENDIF
               READ (LINE(N)(ICOL:ICOL+14), '(E15.8)')
     1               A(NA, NRANGE, KNUM)
            ENDIF
            ICOL = ICOL + 15
   75    CONTINUE
   80 CONTINUE
      IF (NRANGE .GT. 1) THEN
            WRITE (LOUT, '(6X,A,I3,A)')
     1      'Error...data represents less than ',NT(KNUM)-1,
     2      ' temperature ranges...'
         KERR = .TRUE.
      ENDIF
C
   90 CONTINUE
C     set thermo flag for this species
      ITHRM(KNUM) = .TRUE.
C     this species data may be needed in another phase:
      IF (KTIMES .GT. 1) THEN
         DO 95 K = KNUM+1, KKTOT
            IF (KNAME(K) .EQ. KNAME(KNUM)) THEN
               ITHRM(K) = ITHRM(KNUM)
               KCHRG(K) = KCHRG(KNUM)
               KPHSE(K) = KPHSE(KNUM)
               WT(K)   = WT(KNUM)
               NT(K)    = NT(KNUM)
               DO 91 M = 1, NELEM
                  KNCF(M,K) = KNCF(M,KNUM)
   91          CONTINUE
               DO 93 N = 1, NT(K)
                  T(N,K) = T(N,KNUM)
                  DO 92 L = 1, NFIT
                    A(L, N, K) = A(L, N, KNUM)
   92             CONTINUE
   93          CONTINUE
            ENDIF
   95    CONTINUE
      ENDIF
C     reset defaults
      GO TO 25
C
C     end of SUBROUTINE CKTHRM thermo input
C
  400 CONTINUE
C     still have last species to process?
      IF (NLINES .GT. 0) GO TO 500
      RETURN
C
22222 CONTINUE
      WRITE (LOUT,*) ' Error reading thermodynamic data...'
      KERR = .TRUE.
C
      RETURN
      END
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKUNIT (SUB, NSUB, LOUT, AUNITS, EUNITS, IUNITS,
     1                   NONCON, LMOTZ)
C
C  START PROLOGUE
C
C  SUBROUTINE SKUNIT
C  processes a user's character string input to
C  (1)  toggle type of units for the input Arrhenius parameters
C  (2)  toggle conservation/non-conservation of sites
C  (3)  toggle Motz-Wise correction of stocking coefficients
C
C  Arguments:
C  SUB(*)     - Character string array
C  NSUB       - Integer scalar, count of character strings
C  LOUT       - Integer scalar, formatted output file unit number.
C  AUNITS     - Character string, representation of A units,
C               the pre-exponential factor.
C  EUNITS     - Character string, representation of E units,
C               the activation energy.
C  IUNITS     - Character string, reaction units legend for printing.
C  NONCON     - Logical, non-conservation flag.
C  LMOTZ      - Logical, Motz-Wise correction flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) SUB(*), IUNITS, EUNITS, AUNITS
      CHARACTER*4 CKCHUP, KEY
      LOGICAL NONCON, LMOTZ
      INTEGER CKLSCH
C
      DO 85 N = 1, NSUB
C
         KEY = ' '
         KEY = CKCHUP(SUB(N), 4)
C
C        units legend
         IND = CKLSCH(IUNITS) + 1
C
         IF (KEY .EQ. 'NONC') THEN
C           reactions which do not conserve surface sites
C           declared to be legitimate
            NONCON = .TRUE.
         ELSEIF (KEY .EQ. 'MWOF') THEN
C           user does not want to include Motz-Wise correction
C           in converting a sticking coefficient to a rate constant
            LMOTZ = .FALSE.
C
         ELSEIF (KEY .EQ. 'MWON') THEN
C           user wants Motz-Wise correction
            LMOTZ = .TRUE.
C
         ELSEIF (KEY .EQ. 'CAL/') THEN
C           input line contains the keyword for CAL/MOLE
            IF (IUNITS .EQ. ' ') THEN
C              the character variable to be printed is currently
C              empty; write energy units into IUNITS
               IUNITS = 'E units cal/mole'
            ELSE
C              there is already some information in the variable
C              IUNITS, so add energy information to the end
               IUNITS(IND:) = ', E units cal/mole'
            ENDIF
C           units for energy also get put into EUNITS
            EUNITS = 'CAL/'
C
C        check for other unit specifications for energy:
         ELSEIF (KEY .EQ. 'KCAL') THEN
            IF (IUNITS .EQ. ' ') THEN
               IUNITS = 'E units Kcal/mole'
            ELSE
               IUNITS(IND:) = ', E units Kcal/mole'
            ENDIF
            EUNITS = 'KCAL'
         ELSEIF (KEY .EQ. 'JOUL') THEN
            IF (IUNITS .EQ. ' ') THEN
               IUNITS = 'E units Joules/mole'
            ELSE
               IUNITS(IND:) = ', E units Joules/mole'
            ENDIF
            EUNITS = 'JOUL'
         ELSEIF (KEY .EQ. 'KJOU') THEN
            IF (IUNITS .EQ. ' ') THEN
               IUNITS = 'E units Kjoule/mole'
            ELSE
               IUNITS(IND:) = ', E units Kjoule/mole'
            ENDIF
            EUNITS = 'KJOU'
         ELSEIF (KEY .EQ. 'EVOL') THEN
            IF (IUNITS .EQ. ' ') THEN
               IUNITS = 'E units eV'
            ELSE
               IUNITS(IND:) = ', E units eV'
            ENDIF
            EUNITS = 'EVOL'
         ELSEIF (KEY .EQ. 'KELV') THEN
            IF (IUNITS .EQ. ' ') THEN
               IUNITS = 'E units Kelvins'
            ELSE
               IUNITS(IND:) = ', E units Kelvins'
            ENDIF
            EUNITS = 'KELV'
C
         ELSEIF (KEY .EQ. 'MOLE' .AND.
     1      (SUB(N)(5:5).EQ.'S' .OR. SUB(N)(5:5).EQ.'s')) THEN
C           user wants units of moles
            IF (IUNITS .EQ. ' ') THEN
C              the character variable to be printed is currently
C              empty, so put the units specification into it
               IUNITS = 'A units moles'
            ELSE
C              there is already some information in the string,
C              so append the unit information
               IUNITS(IND:) = ', A units moles'
            ENDIF
C           also put pre-exponential unit information into string
            AUNITS = 'MOLE'
C
         ELSEIF (KEY .EQ. 'MOLE' .AND.
     1      (SUB(N)(5:5).EQ.'C' .OR.SUB(N)(5:5).EQ.'c')) THEN
C           user wants "MOLECULES" pre-exponential units
            IF (IUNITS .EQ. ' ') THEN
C              the character variable to be printed is currently
C              empty, so put the units specification into it
               IUNITS = 'A units molecules'
            ELSE
C              there is already some information in the string,
C              so append the unit information
               IUNITS(IND:) = ', A units molecules'
            ENDIF
C           also put pre-exponential unit information into string
            AUNITS = 'MOLC'
         ELSE
C
            WRITE (LOUT, '(6X,A,A)')
     1      'Warning...unrecognized units expression ignored...',
     2      SUB(N)
         ENDIF
   85 CONTINUE
C
C     supply defaults if units were not specified
C
      IF (AUNITS .EQ. ' ') THEN
         AUNITS = 'MOLE'
         IND = MAX(CKLSCH(IUNITS), 1)
         IF (IND .GT. 1) THEN
            IUNITS(IND+1:) = ', A units moles'
         ELSE
            IUNITS(IND:) = ' A units moles'
         ENDIF
      ENDIF
C
      IF (EUNITS .EQ. ' ') THEN
         EUNITS = 'CAL/'
         IND = MAX(CKLSCH(IUNITS), 1)
         IF (IND .GT. 1) THEN
            IUNITS(IND+1:) = ', E units cal/mole'
         ELSE
            IUNITS(IND:) = ' E units cal/mole'
         ENDIF
      ENDIF
C
C     end of SUBROUTINE SKUNIT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SPREAC
     1      (LOUT, MDIM, KDIM, IDIM, MXPHSE, MAXSPR, NSPAR, MAXORD,
     2       NEDPAR, NYPAR, NSCOV, PNAME, NSPEC, NREAC, SPAR, RSPAR,
     3       AXUNIT, EXUNIT, NUNK, NU, RNCF, KCHRG, ISTK, ICOV, CPAR,
     5       KCOV, IREV, KNCF, NONCON, IDUP, KFIRST, KLAST, IRNU, RNU,
     6       IEDP, KEDP, PEDP, IYLD, KYLD, PYLD, YNCF, SKMIN, IBHM,
     7       KBHM, IORD, KORD, RORD, PDEN, EQFAC, KNAME, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SPREAC
C  coordinates checking of a reaction after all auxiliary lines have
C  been processed;
C  checking is done for mass and charge balance, for duplication,
C  for required reactants, options, parameters, etc., and
C  conversion of units is done here.
C
C  Arguments:
C  LOUT     - Integer scalar,formatted output file unit number.
C  MDIM     - Integer scalar, maximum number of elements.
C  KDIM     - Integer scalar, maximum number of species.
C  IDIM     - Integer scalar, maximum number of reactions
C  MXPHSE   - Integer scalar, maximum number of phases.
C  MAXSPR   - Integer scalar, maximum number of reaction species.
C  NSPAR    - Integer scalar, required number of Arrhenius parameters.
C  MAXORD   - Integer scalar, maximum number of reaction change-orders.
C  NEDPAR   - Integer scalar, required number of temperature-
C             dependence parameters.
C  NYPAR    - Integer scalar, required number of yield parameters.
C  NSCOV    - Integer scalar, required number of coverage parameters.
C  PNAME(*) - Character-string array, phase names.
C  NSPEC(*) - Integer array, reaction's count of reactants+products.
C  NREAC(*) - Integer array, reaction's count of reactants only.
C  SPAR(*,*)- Real matrix, reaction Arrhenius parameters.
C  RSPAR(*,*)-Real matrix, reverse parameters, if given.
C  AXUNIT   - Character string, represents A factor units type.
C  EXUNIT   - Character string, represents E factor units type.
C  NUNK(*,*)- Integer matrix, reaction's species indices.
C  NU(*,*)  - Integer matrix, reaction's stoichiometric coefficients.
C  RNCF(*,*)- Real matrix, reaction's net change for phases.
C  KCHRG(*) - Integer array, species electronic charges.
C  ISTK(*)  - Integer array, reaction indices of sticking reactions.
C  ICOV(*)  - Integer array, reaction indices of coverage reactions.
C  CPAR(*,*)- Real matrix, reaction coverage parameters, if given.
C  KCOV(*)  - Integer array, reaction coverage species, if given.
C  IREV(*)  - Integer array, reaction indices for those with explicit
C             reverse parameters.
C  KNCF(*,*)- Integer matrix, elemental composition of species.
C  NONCON   - Logical, flag for site non-conservation.
C  IDUP(*)  - Logical array, reactions duplication flag.
C  KFIRST(*)- Integer array, starting species indices for phases.
C  KLAST(*) - Integer array, ending species indices for phases.
C  IRNU(*)  - Integer array, real stoichiometry reaction indices.
C  RNU(*,*) - Real matrix, real stoichiometric coefficients.
C  IEDP(*)  - Integer array, energy-dependent reaction indices.
C  KEDP(*)  - Integer array, energy-dependence species indices.
C  PEDP(*)  - Integer array, energy-dependence parameters.
C  IYLD(*)  - Integer array, yield-modify reaction indices.
C  KYLD(*,*)- Integer matrix, yield-modify flags for reaction species.
C  PYLD(*,*)- Real matrix, yield-modify parameters.
C  YNCF(*)  - Real matrix, phase balance due to yield-modify.
C  SKMIN    - Real scalar, a minimum used in balancing.
C  IBHM(*)  - Integer array, Bohm-correction reaction indices.
C  KBHM(*)  - Integer array, ion species index for Bohm-correction.
C  IORD(*)  - Integer array, change-order reaction indices.
C  KORD(*,*)- Integer matrix, change-order species indices.
C  RORD(*,*)- Real matrix, change-order species values.
C  PDEN(*)  - Real array, phase densitities.
C  EQFAC(*) - Real array, reaction equilibrium constant scalar.
C  KNAME(*) - Character-string array, species names
C  KERR     - Logical error flag
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER          (RU_JOUL = 8.314510D0, AVOG = 6.0221367D23)
C*****END precision > double
C*****precision > single
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C       PARAMETER         (RU_JOUL = 8.314510E0, AVOG = 6.0221367E23)
C*****END precision > single
C
C     (Value of Avogadro's Constant and Gas Constant
C      from 1986 CODATA recommended values (1993 CRC)
C      J. Research National Bureau of Standards, 92, 95, 1987
C      6.0221367(39)E23 mol-1
C      8.314510 Joules / (Mole K)       )
C
      COMMON /SKINT/ LENISK, LENRSK, LENCSK, NMAT, NELEM, MELECT,
     1               KELECT, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NFSURF, NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     3               NIISUR, NKION, NCOV, NREV, NSTK, NCON, NRNU, NORD,
     4               NBHM, NEDP, NYLD, MORE
C     Integer arrays
      DIMENSION NSPEC(IDIM), NREAC(IDIM), NUNK(MAXSPR,IDIM),
     1          NU(MAXSPR,IDIM), KCHRG(KDIM), ISTK(IDIM), ICOV(IDIM),
     2          KCOV(KDIM), IREV(IDIM), KNCF(MDIM,KDIM), IDUP(IDIM),
     3          KFIRST(MXPHSE), KLAST(MXPHSE), IRNU(IDIM), IEDP(IDIM),
     4          IYLD(IDIM), KYLD(MAXSPR,IDIM), IBHM(IDIM), IORD(IDIM),
     5          KORD(MAXORD,IDIM), KEDP(IDIM), KBHM(IDIM)
C     Real arrays
      DIMENSION SPAR(NSPAR,IDIM), RSPAR(NSPAR,IDIM),
     1          RNCF(MXPHSE,IDIM), CPAR(NSCOV,IDIM), RNU(MAXSPR,IDIM),
     2          PEDP(NEDPAR,IDIM), PYLD(NYPAR,IDIM),
     3          YNCF(MXPHSE,IDIM), RORD(MAXORD,IDIM),
     4          PDEN(MXPHSE), EQFAC(IDIM)
C
      CHARACTER*(*) AXUNIT, EXUNIT
      CHARACTER*16 PNAME(MXPHSE), KNAME(KDIM)
      LOGICAL KERR, NONCON
C
      IF (NNSURF .GT. 0) THEN
C        calculate the multiplicative factor for surface equilibrium
C        constant
         CALL SKIEQ (KDIM, IDIM, MXPHSE, MAXSPR, KFIRST, KLAST, KCOV,
     1               NU(1,NIISUR), NUNK(1,NIISUR), PDEN, IRNU, RNU,
     2               SKMIN, EQFAC(NIISUR))
      ENDIF
C
C     check for reaction balance
      CALL SKBAL (LOUT, MDIM, KDIM, IDIM, MXPHSE, MAXSPR,
     1            NU(1,NIISUR), NUNK(1,NIISUR), KCHRG, KCOV,
     2            NONCON, IRNU, RNU, IYLD, KYLD, KFIRST, KLAST,
     3            KNCF, RNCF(1,NIISUR), YNCF, PNAME, SKMIN, KERR)
C
C     check for error relating to duplicate reactions
      CALL SKDUP (LOUT, IDIM, NIISUR, MAXSPR, NSPEC, NREAC, NU, NUNK,
     1            IDUP, KERR)
C
      IF (NIISUR .EQ. IBHM(MAX(NBHM,1))) THEN
         CALL SKPBHM (LOUT, KDIM, IDIM, MAXSPR, MAXORD, NIISUR,
     1                KKGAS, KNAME, KELECT, KCHRG, SKMIN,
     2                NREAC(NIISUR), NUNK(1,NIISUR), NU(1,NIISUR),
     3                NRNU, IRNU, RNU,  NORD, IORD,
     4                KORD, RORD, KBHM(NBHM), KERR)
      ENDIF
C
      IF (NIISUR .EQ. IEDP(MAX(NEDP,1))) THEN
         CALL SKPEDP (LOUT, KDIM, IDIM, NIISUR, KKGAS, KNAME, KCHRG,
     1               SKMIN, NREAC(NIISUR), NUNK(1,NIISUR), NU(1,NIISUR),
     2               MAXSPR, NRNU, IRNU, RNU, NORD, IORD, MAXORD, KORD,
     3               RORD, KEDP(NEDP), KERR)
      ENDIF
C
      IF (NIISUR .EQ. ISTK(MAX(NSTK,1))) THEN
         CALL SKPSTK (LOUT, KDIM, IDIM, NIISUR, KKGAS, KNAME, SKMIN,
     1                NREAC(NIISUR), NUNK(1,NIISUR), NU(1,NIISUR),
     2                MAXSPR, NRNU, IRNU, RNU, NORD, IORD, MAXORD, KORD,
     3                RORD, KERR)
      ENDIF
C
      IF (NIISUR .EQ. IYLD(MAX(NYLD,1))) THEN
         CALL SKPYLD (LOUT, KDIM, IDIM, NIISUR, KKGAS, KNAME, KCHRG,
     1                SKMIN, NREAC(NIISUR), NUNK(1,NIISUR),
     2                NU(1,NIISUR), MAXSPR, NRNU, IRNU, RNU, NORD, IORD,
     2                MAXORD, KORD, RORD, KELECT, KERR)
      ENDIF
C
      IF (AXUNIT .EQ. 'MOLC') THEN
C
C        convert units of pre-exponential to standard units
         IF (NIISUR .EQ. ISTK(MAX(NSTK,1))) THEN
C           no unit conversion if a sticking coefficient (no units)
         ELSEIF (NIISUR .EQ. IBHM(MAX(NBHM,1))) THEN
C           no unit conversion if Bohm velocity reaction
         ELSE
C           sum of stoichiometric coefficients of reactants
            NSTOR = 0
            DO 50 N = 1, NREAC(NIISUR)
               IF (NUNK(N,NIISUR) .LE. KKGAS+KKSURF)
     1         NSTOR = NSTOR + ABS(NU(N,NIISUR))
   50       CONTINUE
C           convert from molecules to moles
            NSTOR = NSTOR - 1
            IF (NSTOR .GT. 0)
     1         SPAR(1,NIISUR) = SPAR(1,NIISUR) * AVOG**NSTOR
         ENDIF
C
C        convert units for reverse rate coefficient, if given
         IF (NIISUR .EQ. IREV(MAX(NREV,1))) THEN
            NSTOP = 0
            DO 60 N = MAXSPR/2+1, MAXSPR
               IF (NUNK(N,NIISUR).GT.0 .AND.
     1             NUNK(N,NIISUR).LE.KKGAS+KKSURF)
     1         NSTOP = NSTOP + NU(N,NIISUR)
   60       CONTINUE
            NSTOP = NSTOP - 1
            IF (NSTOP .GT. 0) RSPAR(1,NREV) = RSPAR(1,NREV)
     1                        * AVOG**NSTOP
         ENDIF
      ENDIF

C     find the appropriate conversion factor between the energy unit
C     specified by the user and Kelvin units
C
      IF (EXUNIT .EQ. 'KELV') RETURN
      EFAC = 1.0
      IF (EXUNIT .EQ. 'CAL/') THEN
C        convert E from cal/mole to Kelvin
         EFAC = 4.184  / RU_JOUL
      ELSEIF (EXUNIT .EQ. 'KCAL') THEN
C        convert E from kcal/mole to Kelvin
         EFAC = 4184.0 / RU_JOUL
      ELSEIF (EXUNIT .EQ. 'JOUL') THEN
C        convert E from Joules/mole to Kelvin
         EFAC = 1.0    / RU_JOUL
      ELSEIF (EXUNIT .EQ. 'KJOU') THEN
C        convert E from Kjoule/mole to Kelvin
         EFAC = 1000.0 / RU_JOUL
      ELSEIF (EXUNIT .EQ. 'EVOL') THEN
C        convert E from eV to Kelvin
         EFAC = 1.60219E-19 * AVOG / RU_JOUL
      ENDIF
C
C     convert activation energies to Kelvin
C
      SPAR(3,NIISUR) = SPAR(3,NIISUR) * EFAC
C
      IF (NIISUR .EQ. IREV(MAX(NREV,1))) THEN
C       reverse coefficients were given; convert reverse
C       activation energy to Kelvin
        RSPAR(3,NREV) = RSPAR(3,NREV) * EFAC
      ENDIF
C
      IF (NIISUR .EQ. ICOV(MAX(NCOV,1))) THEN
C        conversion of coverage parameter units
         CPAR(3,NCOV) = CPAR(3,NCOV) * EFAC
      ENDIF
C
      IF (NIISUR .EQ. IEDP(MAX(NEDP,1))) THEN
C        ion-energy dependence of reaction rate declared;
C        convert energy units and reaction pre-exponential
         PEDP(1,NEDP) = PEDP(1,NEDP) * EFAC
         SPAR(1,NIISUR) = SPAR(1,NIISUR) *
     1       EFAC**(-PEDP(2,NEDP)*PEDP(3,NEDP))
      ENDIF
C
      IF (NIISUR .EQ. IYLD(MAX(NYLD,1))) THEN
C        ion-yield dependence of reaction rate declared;
C        convert yield threshold energy and scale factor
C        PYLD(1)=A, PYLD(2)=Eth[eV], PYLD(3)=a, PYLD(4_=b
C
         PYLD(2,NYLD) = PYLD(2,NYLD) * EFAC
         PYLD(1,NYLD) = PYLD(1,NYLD) *
     1                  EFAC**(-PYLD(3,NYLD)*PYLD(4,NYLD))
      ENDIF
C
C     end of SUBROUTINE SPREAC
      RETURN
      END
      SUBROUTINE SKIABS
C
C     The surface chemistry interpreter reads a CHEMKIN binary file
C     LINCK, containing information about gas-phase chemistry for
C     a problem, then reads and interprets a surface mechanism data
C     file and creates a binary file LINSK for surface chemistry.
C
C///////////////////////////////////////////////////////////////////
C
C     SKINTERP: SURFACE KINETICS INTERPRETER
C                      VERSION 7.10
C
C     WRITTEN BY:
C         FRAN M. RUPLEY
C         COMPUTATIONAL MECHANICS DEPARTMENT, 8745
C         MAIL STOP 9042
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94550
C         (415) 294-3657
C
C       Copyright 1990, Sandia Corporation.
C       The U.S. Goverment retains a limited license in this
C       software.
C
C       The U.S. Government retains, in this software, a paid-up,
C       nonexclusive, irrevocable worldwide license to reproduce,
C       prepare derivative works, perform publicly and display
C       publicly by or for the Government, including the right to
C       distribute to other Government contractors.
C
C       Neither the United States, the U.S. Dept. of Energy, nor
C       any of their employees, makes any warranty, express or
C       implied, or assumes any legal liability or responsibility
C       for the accuracy, completeness, or usefulness of any
C       information, apparatus, product, or process disclosed, or
C       represents that its use would not infringe privately owned
C       rights.
C
C/////////////////////////////////////////////////////////////////////
C       v.7.10 98/03/25
C       1. Action#161:  Changed the name of SKSTRT to SKBGIN, to avoid
C          conflict with the common block /SKSTRT/ in sklib.f
C       v.7.9 98/03/03
C       1. Action#119:  Removed duplicate SKPCMP routine from 
C          skinterp.f and added sklib.f to link for surf.exe in 
C          the Makefile. (E. Meeks)
C       v.7.8 97/08/18
C       1. Fix bug#074a: in SKTHRM check NVAL>0 to avoid possible
C          IELEM=VALUE(1) error
C       v.7.7 97/08/04
C       1. Fix bug#069: in SKBAL, in loops 221 and 222, put check
C          for whether or not new and old species are #-modified 
C          (i.e. part of YIELD auxiliary option) before determining
C          whether to add stoichiometric coefficients to previously 
C          found values.  Remove logic following these loops that
C          attempted to do the same thing.
C          (E.Meeks)
C       v.7.6 97/07/31
C       1. Fix bug#068: SKDUP, split IF tests to avoid out-of-bounds
C          errors where index KJ may be zero (F. Rupley)
C       v.7.5 97/07/23
C       1. Fix bug#012: remove spacing from logical operators
C          to conform to f77 standard. Spaces found in SKBULK,
C          SKSTRT (2), and SKSURF.
C       V.7.4 97/07/22
C       1. SKTHRM, comment out IF (KERR) RETURN in order to
C          continue through thermo data in spite of previous error
C       V.7.3 97/07/16
C       1. (2) comment lines shortened in order to use CHANGE program
C       V.7.1 97/04/16
C       1. need KERR set if CKLEN fails
C       V.7.0 97/04/15
C       1. logical function SKFILE is used in SKKEY to check for
C          existence of opened thermodynamics file;
C          CHEMKIN Committee bugfix #001
C       2. also need KERR set if CKINIT fails
C       V.6.9 97/03/27
C       1. SKBAL needs YERR initialized .FALSE.
C       DEVELOPMENT V.6.8 97/03/26
C       1. SKBAL needs INTEGER CKLSCH and EXTERNAL CKLSCH
C       2. SKBULK needs INTEGER CKLSCH and EXTERNAL CKLSCH
C       DEVELOPMENT V.6.7 97/03/25
C       1. have corrected bugs in SKBAL reaction balancing
C       2. add MORE to ascii linkfile WRITE in CKFORM
C       V.6.6 97/03/02 per M. Coltrin
C       1. in SKNORD fix DO 110 to start at MAXSPR/2+1
C       V.6.5 97/03/01 F. Rupley
C       1. make new main "driver" programto set up arrays,
C          to satisfy Reaction Design requirement of providing object
C          files instead of source code - move OPENs to driver
C       V.6.4, 97/01/25 F. Rupley
C       1. allow more than 3 fit temperatures
C       V.6.3, 97/01/21 F. Rupley
C       1. allow multiple COV declarations in a reaction
C       V.6.2, 96/12/11 F. Rupley
C       1. several bugfixes
C       V.6.1, 96/11/11 F. Rupley
C       1. renamed and incorporated some of the utility parsing
C          routines into cklib.f, so deleted from skinterp.f
C       V.6.0, 96/09/04 F. Rupley
C       1. delete VALUE(5) from SKNREV
C       VERSION 5.9, 96/08/14 F. Rupley
C       1. change SKDUP to search for differences, rather than sameness
C       2. change SKPBHM to fix logic looking for the ion species
C       3. change SKNORD to allow any species to change order, not only
C          a gas-phase species; this requires additional arguments in
C          the call list
C       VERSION 5.8, 96/08/06 F. Rupley
C       1. OPEN therm.dat as "OLD" and initialize ISTR=' '
C       VERSION 5.7, 96/08/05 per E. Meeks
C       1. in SKPRNT no error if MELECT.GT.0
C       VERSION 5.6, 96/05/24
C       1. initial sccs version
C       VERSION 5.5 (F. Rupley 5/14/96)
C       1. don't print Motz legend if no sticking reactions
C       VERSION 5.4 (F. Rupley 5/2/96)
C       1. declare LYLD logical in SKAUXL, LIMOTZ in SKKEY
C       VERSION 5.3 (F. Rupley 5/1/96)
C       1. Switch to MDIM, KDIM, IDIM for dimensioning in SKBIN/SKFORM,
C          instead of the previously utilized MM, KK, II (or Nxxx),
C          which might be zero.
C       VERSION 5.2 (F. Rupley 4/29/96, per E. Meeks)
C       1. Verified that Bohm option doe not have requirement for an
C          electron species reactant (fixed comments that would
C          indicate otherwise); maintain that Bohm option requires
C          exactly one non-electron ion with order approx. 1, and no
C          other gas-phase participant species.  Still need to
C          verify that an electron species is present in the mechanism
C          to use this option, since the rate calculation needs to
C          know its species index.
C       2. Yield-modify reaction does not require electron species
C          in the mechanism or in the reaction; it does require exactly
C          one ion reactant, but its coefficient is not restricted to 1,
C          and there may be other gas-phase species in the reaction.
C
C       VERSION 5.1 (F. Rupley April 24, 1996)
C       CAUTION...THIS WORK IS STILL IN PROGRESS:
C       1.  declare LORD logical in SKPEDP
C       2.  REWIND LTHRM in SKTHRM
C       3.  ascii linkfile name "surf.asc"
C       VERSION 5.0 (Initial CHEMKIN-III version, F. Rupley Jan. 1996)
C       1.  binary/ascii linkfile options.
C       2.  '#' yield-modify for reactions.
C       3.  linkfile version numbers separate from interpreter version.
C       4.  PROLOGUE documentation.
C       5.  separate SKN... subroutines for auxiliary reaction options.
C       6.  MOTZ flag is now an array over sticking reactions
C       7.  separate SKP... subroutines called by SPREAC to verify
C           existence of required species types after all auxiliary
C           reaction input
C       8.  fully dimensioned arrays in (most) subroutines
C       9.  COMMON /SKINT/ to carry counter variables
C
C        CHANGES FOR VERSION 4.9 (8/18/95 E. Meeks)
C        1.  Add auxiliary keyword UNITS/string/ to specify units
C            differently for a specific reaction.  /string/ can
C            contain any units descriptors from the REACTION
C            line.
C        2.  Moved initialization of AUNITS out of SKUNIT (passed in).
C        3.  Call SKUNIT from SKAUXL to set AXUNIT(NIISUR),
C            EXUNIT(NIISUR) unit expressions for individual reactions.
C        4.  Initialize AXUNIT and EXUNIT arrays in SKSET.
C        5.  Apply reaction-specific unit conversion info in SPREAC.
C        CHANGES FOR VERSION 4.8 (8/15/95 E. Meeks)
C        1.  Corrected calculation of NKION to include negative
C            ions as well as positive ions.
C        2.  Increased LENRWK (+KKGAS) and LENIWK (+1) to allow new
C            pointers in sklib.f (formed in SKINIT).
C        CHANGES FOR VERSION 4.7 (6/2/95 M. Coltrin)
C        1.  Add auxiliary keyword ENRGDEP and add associated variables
C            to linkfile (which changes with this version).
C        2.  Allow users to specify energies in eV (EVOL keyword on
C            REAC line).
C        3.  Fix mistake in parsing "KCAL/MOL" unit specifier;
C            program was comparing to string 'KCAL/' instead of
C            'KCAL'.
C        4.  Change all occurances of AVAG to AVOG (mis-spelling of
C            Mssr. Avogadro's name).
C        5.  Add KION(*), NKION, KELECT to linkfile.
C        6.  Delete all occurances of array IBT, which is no
C            longer needed.
C        7.  Change BOHM auxiliary keyword parsing logic; ionic
C            species name and temperture pointer are no longer
C            input parameters.
C        8.  Add ENRGDEP auxiliary keyword. This adds the following
C            variables to linkfile: NEDPAR, NEDP, IEDP(*),
C            KEDP(*), PEDP(NEDPAR,*).
C        9.  In SKSET, changed initialization of KORD(I) from 0.0 to 0
C        10. Loop 5 in SKPRNT used a format number 1660, which did not
C            apply to the error found. Changed it to a new format
C            1750.
C        11. IF-test for RSPAR(1,NREV) in SKAUXL compared to "0" instead
C            of "0.0".
C        12. Add KELECT to arguments of SKAUXL and SKKEY in order
C            to do an error check for BOHM reactions.
C            (See format 3097.)
C        CHANGES FOR VERSION 4.6 (2/27/95 F. Rupley)
C        1.  Change character index "(:" to "(1:"
C        2.  In SKAUXL add check to see if order has changed before
C            applying error check in STK option (per E. Meeks)
C        CHANGES FOR VERSION 4.5 (10/3/94 F. Rupley)
C        1.  Add NORD to argument list for SKSET and initialize to 0.
C        CHANGES FOR VERSION 4.31 (8/10/94 H. Moffat)
C        1.  Changed physical constants to conform to 1986 CODATA
C            recommendations
C        CHANGES FOR VERSION 4.3
C        1.  Lengthen ISKWRK by 3
C        CHANGES FOR VERSION 4.2
C        1.  Correct LENISK (NBHM > NBHM)
C        CHANGES FOR VERSION 4.10 (7/14/94 F. Rupley)
C        1.  Implement MATERIAL keyword to enable multiple materials
C        2.  Subroutine SKSET to initialize values
C        CHANGES FOR VERSION 4.09 (6/28/94 E. Meeks)
C        1.  Add REACTION line keyword MWOFF to turn off Motz-Wise
C            correction to sticking coefficient formulation.
C            Addition integer flag MOTZ passed through linkfile.
C        CHANGES FOR VERSION 4.08c (6/3/94 F. Rupley per H. Moffat)
C        1.  Add error checks/messages upon opening therm.dat,
C            chem.bin, surf.inp
C        2.  Allow comments (!) in thermo data
C        3.  Correct phase name logic
C        CHANGES FOR VERSION 4.08b (5/20/94 F. Rupley per E. Meeks)
C        1.  Incorporate plasma options.
C        CHANGES FOR VERSION 4.08 (4/29/94 F. Rupley)
C        1.  Cannot change RORD for an irreversible reaction.
C        CHANGES FOR VERSION 4.07 (4/19/94 F. Rupley)
C        1.  correct indexing in SKBAL, SKRBAL
C        CHANGES FOR VERSION 4.06 (3/15/94 F. Rupley)
C        1.  DOS/PC compatibility effort includes adding filenames to
C            OPEN statements, removing unused variables in CALL lists,
C            unusued but possibly initialized variables.
C         CHANGES FOR VERSION 4.05 (1/26/94 F. Rupley per R. Kee)
C         1. Allow real stoichometric coefficients; NRNU, IRNU(*),
C            RNU(MAXSPR,*)
C         CHANGES FOR VERSION 4.04 (4/13/92 F. Rupley per M. Coltrin)
C         1. Correct logic for SKDUP, requires additional argument.
C         CHANGES TO MAKE VERSION 4.03
C         1. TMID must be passed to SKSPEC
C         CHANGES TO MAKE VERSION 4.02
C         (7/17/91 F. Rupley per M. Coltrin)
C         1. SKPRNT was checking sites, bulks, and species and
C            printing error messages for an empty surface mechanism;
C            corrected to skip these checks if there are no sites
C            or bulks.
C         CHANGES TO MAKE VERSION 4.01
C         (7/11/91 F. Rupley per M. Coltrin)
C         1. Additional error checking:
C            no species on a site or bulk
C            site occupancy numbers > 0
C            site density required and > 0
C            bulk density > 0 if given
C            lack of Arhennius coefficients
C            duplicate phase names
C         CHANGES TO MAKE VERSION 4.0
C         1. Change units on surface-coverage modification of the
C            rate of progress (concentration units --> site fractions)
C         2. Add an extra NIISUR to length of real work array required
C            to account for new multiplicative factor (EQFAC) in SKLIB
C         CHANGES TO VERSION 3.78
C         1. Increase length LENRSK by NIISUR+NREV for perturbation
C            factor.
C         CHANGES TO VERSION 3.77
C         1. Change "LINE" to "ILINE" for reading thermodynamic
C            database.
C         CHANGES TO VERSION 3.76
C         1. Need to get TLO,THI,TMID from thermodynamic database
C            BEFORE reading user's THERMO data
C         CHANGES TO VERSION 3.75
C         1. Additional thermodynamic checks:
C            a) if TLO,TMID,THI not given on THERMO CARDS, use values
C               from THERMO.DAT
C            b) Check TLO < THI, TLO <= TMID, TMID <= THI
C         CHANGES TO VERSION 3.74
C         1. Correct error similar to V3.72 for RSPAR(1, NREV)*AFAC,etc.
C         CHANGES TO VERSION 3.73
C         1. Initialize logical variable LPDEN
C         CHANGES TO VERSION 3.72
C         1. Indexing in scaling of 3rd Reverse Arrhenius parameter
C            storage was wrong, should be RSPAR(3, NREV)*EFAC,
C            not RSPAR(3, NIISUR).
C         CHANGES TO VERSION 3.71
C         1. Error in UPCASE was causing auxiliary keywords to be
C            ignored.
C         CHANGES TO VERSION 3.7
C         1. INCF(*,*) previously calculated for sites only, is now
C            calculated for all phases (SUBROUTINE SKBAL).
C         CHANGES TO VERSION 3.64
C         1. Modify reaction interpretation to allow continuation
C            lines as an array (thus avoiding CHARACTER*160).
C         CHANGES TO VERSION 3.63
C         1. Set KERR=.TRUE. if coefficient greater than 1 for
C            gas phase species in a reaction with sticking coeff.
C         CHANGES TO VERSION 3.62
C         1. Need NU in argument list for SKAUXL
C         CHANGES TO VERSION 3.61
C         1. Error if site phase or bulk phas has no species declared
C         2. Change SKUNIT to parse LINE instead of SUB(*) to
C            correct misinterpretation of unit strings containing
C            a slash
C         3. Increase length of real work space by NPHASE in order
C            to have scratch space for phases
C         4. Modify SKTHRM such that if a species occurs in more than
C            one site or bulk phase, thermodynamic data will be
C            stored for each occurrence.
C         5. With sticking coefficients, the one gas-phase reactant
C            must have a stoichiometric coefficient of one.
C         6. A site must have a density delcaration.
C         CHANGES TO VERSION 3.6
C         1. Correct unit conversions
C         2. Bring up to date with latest version of manual
C         CHANGES FOR VERSION 3.5
C         1. Add IKCOV array of species for reactions with coverage
C            parameters
C         2. Add integer NCON, the number of reactions which do
C            not conserve sites, to first record of binary file
C         3. Add first record to binary file, which contains a
C            character*16 string VERS to name the version number of
C            the binary file, and PRES to name its precision
C            (DOUBLE or SINGLE)
C         4. Correction to SKBULK and SKSURF to initialize KNUM2=0.
C         5. Check that pre-exponential factor is positive.
C       CHANGES TO VERSION 3.4
C         1. 6 reactants, 6 products
C         2. 'NONCON'servation of sites option
C         3. default phase names
C         4. continuation character '&' for reaction lines
C         5. phase names must be unique
C         6. site-phase species names may not duplicate gas-phase
C            species names, but MAY duplicate other site-phase species
C            names not in same site phase.
C         7. bulk-phase species names may not duplicate either gas-
C            or site-phase species names, but MAY duplicate other
C            bulk-phase species names not in same bulk phase.
C         8. site densities are now moles/cm**2 instead of number
C            densities
C       CHANGES TO VERSION 3.3
C         1. Add auxiliary reaction keyword "STICK" for sticking
C            coefficients, in which case there must be one gas-phase
C            reactant.
C       CHANGES TO VERSION 3.2
C         1. Eliminate MIXSOL and PURE in favor of BULK.
C         2. Implement slash-delimited options.
C       CHANGES SINCE VERSION 1.5
C         1. Replace DEPOSIT species with MIXSOL and PURE species
C       CHANGES SINCE VERSION 1.4
C         1. Add "unix" change blocks
C       CHANGES SINCE VERSION 1.3
C         1. Change name to SINTRP to conform to ANSI standard
C         2. Change MAXSITE to MXSITE to conform to ANSI standard
C       CHANGES SINCE VERSION 1.1
C         1. Include thermodynamic properties for all species
C         2. Allow reversible/irreversible surface reactions
C         3. Allow reaction species to end with '=' or '-'
C         4. Allow real values of elemental composition in THERMO cards
C         5. Allow upper/upper case input
C       CHANGES SINCE VERSION 1.0
C         1. Replace "REAL*8" with "DOUBLE PRECISION"
C
C/////////////////////////////////////////////////////////////////////
      RETURN
      END
C
      LOGICAL FUNCTION SKFILE (LUNIT, LFORM)
      INTEGER LUNIT
      CHARACTER LFORM*(*), LTYPE*16, STR1*1
      LOGICAL LOPEN
C
C     INQUIRE returns logical OPENED, character LTYPE
      INQUIRE (LUNIT, OPENED=LOPEN, FORM=LTYPE)
      SKFILE = LOPEN .AND. (LFORM(1:1).EQ.LTYPE(1:1))
      IF (SKFILE) THEN
         IF (LTYPE(1:1) .EQ. 'U') THEN
            READ (LUNIT, ERR=500)
         ELSE
            READ (LUNIT, '(A)', ERR=500) STR1
         ENDIF
         BACKSPACE (LUNIT)
         RETURN
  500    CONTINUE
C        error reading file
         SKFILE = .FALSE.
      ENDIF
      RETURN
      END
