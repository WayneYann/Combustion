C CVS:$Revision: 1.1.1.1 $ created $Date: 2006/05/26 19:09:33 $
C----------------------------------------------------------------------C
      SUBROUTINE STHERM (LEN_IWORK, IWORK, LEN_RWORK, RWORK,
     1                   LEN_CWORK, CWORK,
     2                   LEN_CWRK64, CWRK64, LEN_CWRK48, CWRK48)
C
*=======================================================================
* $Id: surftherm.f,v 1.1.1.1 2006/05/26 19:09:33 marc Exp $
* $Source: /usr/local/cvsroot/chemkin/src/surftherm.f,v $
*
* $Log: surftherm.f,v $
* Revision 1.1.1.1  2006/05/26 19:09:33  marc
* This is part of the chemkin code suite.  Various routines in
* the src folder have been modified for CCSE-specific purposes,
* particularly in the use of premix.  If you need a complete
* set of chemkin functionality, you need to get the original
* chemkin distribution.
*
*
* Revision 1.16  1998/04/24 20:50:30  emeeks
* 1) Clean up change blocks and any long comments,
*    so that change will run.
*
* Revision 1.15  1998/04/24 18:41:48  emeeks
* 1) Shortened comment in line 17, to allow change.exe to not give 
*    error.
*   (Last check-in did not contain message due to remote editing 
*    errors.)
*
* Revision 1.14  1998/04/24 18:04:18  emeeks
* *** empty log message ***
*
* Revision 1.13  1997/07/25 21:05:15  emeeks
* 1.  Fix bug #012: remove extra spaces in logical operators.  Found 4
*     instances of . GT. and 6 instances of . NE.
* 2.  Moved DATA statement in ID_INT after declarations for f77 stndard.
* 3.  Consolidated change-history comments to driver routine and
*     eliminated some old RCS commands (keep $LOG$ for log history).
*
* Revision 1.12  1997/05/01 21:59:02  emeeks
* Set revision to 1.12; Created 97/03/01.
*
* Revision 1.1  1997/05/01 21:23:55  emeeks
* Initial CVS revision.
*
C======================================================================
C Change history previous to CVS installation:
C
C Version 1.12 (97/03/01 F. Rupley)
C 1) make new main "driver" program and add pointers to set up arrays,
C    to satisfy Reaction Design requirement of providing object
C    files instead of source code - move OPENs to driver.
C Version 1.11 (6/19/96: H. K. Moffat)
C  1) Initial Chemkin_III Version. Doesn't encompass any of the new
C     features introduced in Chemkin_III. However, it no longer core
C     dumps.
C  2) Fixed logic error with reversible surface rxns whose reverse
C     rate constant is given by REV keyword. More has to be done
C     for this subcase.
C Version 1.10B (96/11/11 F. Rupley)
C 1) delete SUBROUTINE CKxxx as later versions have been incorporated
C    into cklib.f
C Version 1.10 (H. K. Moffat)
C  1) Bug fix of third body section
C  2) Implemented global variables through the surftherm.h
C     include file. Expect this feature to grow
C Version 1.9 (H. K. Moffat)
C  1) Added static and dynamic memory allocation options
C  2) Changed default carrier gas to be the last species,
C     and the default major reactant to be the first species
C Version 1.8 (H.K. Moffat)
C  1) Compilation error fixed
C VERSION 1.6 (M. Coltrin)
C  1) Add IFLAG to the arguments of CKINIT, CKLEN, SKINIT, SKLEN,
C     MCINIT, MCLEN
C VERSION 1.5 (F. Rupley)
C  1) Added ISKWRK to some CALL SK subroutine calls.
C VERSION 1.45
C  1) Changed .h file to chemkin's ckstrt.h file
C VERSION 1.36:
C  1) Fixed a bug in calculating the Arrhenius parameters for fall-off
C     reaction rates as a function of pressure.
C VERSION 1.37:
C  1) Fixed a bug. The ACT vector was not being initialized before
C     calling NDCALC_SUR.
C VERSION 1.38
C  1) Added pointers and work space storage allocation to the program.
C VERSION 1.39
C  1) Shoved long arrays into blank common, and increased size of
C     fixed dimensions within the program.
C VERSION 1.40
C  1) Changed getarg to igetarg in HP version
C VERSION 1.41
C  1) Delete "EXTERNAL IGETARG" statement
C  2) Made minor changes in the following format statement numbers:
C     7002, 7005, 7008, 8356, 8600, 8601, 9599, 9600, 9700
C  3) Added format statement 8253
C  4) Changed print of bath bulk activities in loop 303 to use format
C     number 8253
C VERSION 1.42
C  1) Changed names of linking files in open statements to conform
C     with new versions of the interpreters
C VERSION 1.43
C  1) Added SKWRK to the arguments of SKNU.
C  2) Fixed unit number in print-out of variable RP using format
C     number 8200 in subroutine P_REAC
C VERSION 1.44
C  1) Fixed duplicate declaration statements.
C  2) Explicit typing option on the HP now works.  Had to move
C     the location of a common block
C  3) NaN's are fixed for HP version, ck_com.h needed to be updated.
C=======================================================================
C
C
C      PROGRAM SURFDRVR
C
C----------------------------------------------------------------------C
C          SURFTHERM
C                    Version 1.13
C
C          WRITTEN BY M. E. COLTRIN and H. K. MOFFAT, 1126
C
C----------------------------------------------------------------------C
C
C Surftherm:
C     WHAT IT DOES
C    ---------------
C
C    Surftherm analyses a gas/surface CHEMKIN mechanism.
C
C     HOW TO USE IT: unix
C   ----------------------
C
C    Surftherm expects to find a CHEMKIN linking file named "chem.bin",
C a SURFACE CHEMKIN linking file named "surf.bin", and a transport
C linking file named "tran.bin" in the current directory.  It will then
C produce an analysis of the mechanism represented in the linking
C files on standard output.
C    Optionally, it can take a command file of keywords to alter the
C default output.  The command file must be named "surftherm.inp" and
C reside in the current directory.  If no such file exists, a default
C set of output is produced.  The list of recognized keywords is
C contained in comments at the top of the SURFKEY subroutine.
C    On unix machines, the command file name can also be input on the
C command line (i.e. surftherm surftherm.inp).  The first command file
C argument is assumed to be the name of the command file.
C
C  F77 Extensions Needed
C ----------------------------
C  0) INCLUDE directive
C  1) Use of the underscore in variable names
C  2) 16 character length long global and local variable names
C  3) $ formatting feature
C  4) Automatic arrays for the optional dynamic memory allocation
C     change block
C----------------------------------------------------------------------C
C
C
C      CHARACTER VERSN*64
C
C Declarations for work arrays
C
      INTEGER    LEN_IWORK,       LEN_CWORK,      LEN_RWORK
C      PARAMETER (LEN_IWORK=10000, LEN_CWORK=3000, LEN_RWORK=40000)

      INTEGER      IWORK(LEN_IWORK)
      CHARACTER*16 CWORK(LEN_CWORK)
      INTEGER LEN_CWRK48, LEN_CWRK64
      CHARACTER*48 CWRK48(LEN_CWRK48)
      CHARACTER*64 CWRK64(LEN_CWRK64)
C*****B-2) precision > double
      DOUBLE PRECISION RWORK(LEN_RWORK)
C*****END B-2) precision > double
C*****B-1) precision > single
C      REAL            RWORK(LEN_RWORK)
C*****END B-1) precision > single
C
C      COMMON        IWORK, CWORK, RWORK
C
      INCLUDE 'surftherm.h'
C
C Declarations for communication with main subroutine
C
      LOGICAL         RESTRT
C
C
C Declarations for pointers into work arrays
C
      INTEGER         N_CKW, N_SKW, N_MCK, N_FPR, N_AKI, N_H__, N_S__,
     1                N_CP_, N_G__, N_DSG, N_DHG, N_DGG, N_DSS, N_DHS,
     2                N_DGS, N_DJK, N_XBT,
     3                I_CKW, I_SKW, I_MCK, I_NEL, I_KSG, I_KSF, I_KSS,
     4                I_KSU, I_NST,
     $                L_THM, L_GRX, L_SRX, L_STK, L_SCV, L_PFL, L_TFL,
     $                L_GTB,
     5                C_CKW, C_SKW,
     6                P_IW,  P_RW,  P_CW
      COMMON /SURPNT/ N_CKW, N_SKW, N_MCK, N_FPR, N_AKI, N_H__, N_S__,
     1                N_CP_, N_G__, N_DSG, N_DHG, N_DGG, N_DSS, N_DHS,
     2                N_DGS, N_DJK, N_XBT,
     3                I_CKW, I_SKW, I_MCK, I_NEL, I_KSG, I_KSF, I_KSS,
     4                I_KSU, I_NST,
     $                L_THM, L_GRX, L_SRX, L_STK, L_SCV, L_PFL, L_TFL,
     $                L_GTB,
     5                C_CKW, C_SKW,
     6                P_IW, P_RW, P_CW
C
C     Pointers extended for ReactionDesign by Fran Rupley
      INTEGER NSCOV
      INTEGER C_KNAM, C_ENAM, C_PNAM, C_GNAM, C_SNAM, P_C64, P_C48,
     1        I_KPHAS, I_KFIRS, I_KLAST, I_KCHRG, I_KOCC, I_NELPH,
     2        I_GREAC, I_GPROD, I_ITHB, I_GREV, I_IFOP, I_IFLO, I_KFAL,
     3        I_GDELT, I_ISTFL, I_MDELS, I_SDELS, I_BDELS, I_NRPPS,
     4        I_ICOV, I_SREV, I_GREACS, I_SREACS, I_BREACS,
     5        I_BPRODS, I_KGASR, I_KGASP, I_NCON, I_KCOVI, I_IMOTZ,
     6        I_IEDP, I_IBHM, I_IORD, I_IYLD, I_NLIN, L_LTIME,
     7        N_WT, N_DEN, N_XEST, N_EPS, N_SIG, N_DIP, N_POL,
     8        N_ZROT, N_VIS, N_CON, N_HBTH, N_SBTH, N_CPBTH,
     9        N_GBTH, N_ACT, N_SDEN, N_H298, N_GEQ, N_GEQD, N_RAG,
     *        N_RBG, N_REG, N_RPG, N_DRCF, N_DRCR, N_SEQ, N_SEQD,
     1        N_RAS, N_RBS, N_RES, N_RPS, N_SDRAT, N_CPARI
C
C Common block for Unit numbers for I/O
C
      INTEGER         LIN, LOUT, LCOM, LINKCK, LINKSK, LINKMC
      LOGICAL         LCOM_EXIST
      COMMON /SF_UIO/ LIN, LOUT, LCOM, LINKCK, LINKSK, LINKMC,
     1                LCOM_EXIST
C
C
C Common block for variables associated with the keyword interface
C
      LOGICAL         LTRAN
      COMMON /SF_KEY/ LTRAN
C
C Other variables associated with table lengths
C
      INTEGER MAX_TTABLE_NUM
C
C
C Variables Associated with command line arguments
C
      CHARACTER ARGV*256
C*****A-1): C-1) Command line arguments: sun
C      INTEGER NUM_ARGS
C      INTEGER IARGC
C      EXTERNAL IARGC, GETARG
C*****END A-1): C-1) Command line arguments: sun
C*****A-1): C-2) Command line arguments: hp
C      INTEGER NUM_ARGS, TMPVAR
C      INTEGER IARGC, IGETARG
C*****END A-1): C-2) Command line arguments: hp
C
C Declarations for IEEE exception handling
C
C*****A-1): B-1) sun debug
C      INTEGER  ieeer, ieee_handler
C      EXTERNAL myhandler,my_undfl_handler, ieee_handler
C*****END A-1): B-1) sun debug
C
C Misc Local Variables
C
      INTEGER LT, IFLAG
      LOGICAL TRAN_INIT
C
C Externals:
C
      INTEGER  CKLSCH
      EXTERNAL CKLSCH
C
C Data Statements
C
C----------------------------------------------------------------------C
C      DATA VERSN
C     1 /'@(#) surftherm.f: SCCS vers 1.1 05/06/96'/
C      DATA LIN/5/, LOUT/6/, LCOM/15/, LINKCK/25/, LINKSK/26/,
C     1     LINKMC/35/
C       DATA ARGV /'surftherm.inp'/, LCOM_EXIST /.FALSE./
        DATA ARGV /'surftherm.inp'/
C       DATA LCOM_EXIST /.TRUE./
C*****A-1): C-1) Command line arguments: sun
C      DATA NUM_ARGS /0/
C*****END A-1): C-1) Command line arguments: sun
C*****A-1): C-2) Command line arguments: hp
C      DATA NUM_ARGS /0/
C*****END A-1): C-2) Command line arguments: hp
      DATA  RESTRT /.FALSE./, TRAN_INIT/.FALSE./, ZERO /0.0/
C----------------------------------------------------------------------C
C    END OF DECLARATION STATEMENTS FOR INITIALIZATION ROUTINE
C----------------------------------------------------------------------C
*
       LCOM_EXIST = .TRUE.
C      LTRAN = .FALSE.
*
C-----------------------------------------------------------------------
C                             BLOCK 1
C                       SET UP ERROR HANDLING
C-----------------------------------------------------------------------
C
C*****A-1): B-1) sun debug
C      ieeer = ieee_handler('set', 'common', myhandler)
CC     ieeer = ieee_handler('set', 'underflow', my_undfl_handler)
C*****END A-1): B-1) sun debug
*
C-----------------------------------------------------------------------
C                              BLOCK 2
C          WRITE OUT HEADER INFORMATION ABOUT THE PROGRAM
C-----------------------------------------------------------------------
*
      WRITE ( LOUT, *)
     $ ' SURFTHERM: Program to analyze gas/surface reaction mechanisms'
      WRITE (LOUT,*)'            surftherm.f, ',
     $              'Version 1.13, 97/07/25'
      WRITE (LOUT,*)'                         SCCS version ',
     $    '1.1, 05/06/96'
C*****B-2) precision > double
      WRITE (LOUT,*)'            DOUBLE PRECISION'
C*****END B-2) precision > double
C*****B-1) precision > single
C      WRITE (LOUT,*)'            SINGLE PRECISION'
C*****END B-1) precision > single
*
C----------------------------------------------------------------------
C                             BLOCK 3
C       OPEN UP THE COMMAND FILE FOR READING
C----------------------------------------------------------------------
C
C*****A-1): C-1) Command line arguments: sun
CC     assign the first argument to the command file
C      NUM_ARGS = IARGC()
C      IF (NUM_ARGS .GT. 0) THEN
C        ARGV = ' '
C        CALL GETARG(1,ARGV)
C        IF (NUM_ARGS .GE. 2) THEN
C          WRITE(LOUT,*)'SURFTHERM WARNING: Only the first argument',
C     1                 ' will be used.'
C        END IF
C        OPEN(UNIT=LCOM, FILE=ARGV, STATUS='OLD',
C     1       FORM='FORMATTED',   ERR=3)
C        LCOM_EXIST = .TRUE.
C        LT = MAX(1, CKLSCH(ARGV))
C        WRITE(LOUT,9001) ARGV(1:LT)
C        GO TO 5
C3       CONTINUE
C        LT = MAX(1,CKLSCH(ARGV))
C        WRITE(LOUT,*)'SURFTHERM ERROR: Could not open command file',
C     $               ' input on the command line, ', ARGV(1:LT)
C        RETURN
C      END IF
C*****END A-1): C-1) Command line arguments: sun
C*****A-1): C-2) Command line arguments: hp
CC     assign the first argument to the command file
C      NUM_ARGS = IARGC()
C      IF (NUM_ARGS .GT. 0) THEN
C        ARGV = ' '
C        TMPVAR = IGETARG(1, ARGV, 256)
C        IF (NUM_ARGS .GE. 2) THEN
C          WRITE(LOUT,*)'SURFTHERM WARNING: Only the first argument',
C     $                 ' will be used.'
C        END IF
C        OPEN(UNIT=LCOM, FILE=ARGV, STATUS='OLD',
C     $       FORM='FORMATTED',   ERR=3)
C        LCOM_EXIST = .TRUE.
C        LT = MAX(1, CKLSCH(ARGV))
C        WRITE(LOUT,9001) ARGV(1:LT)
C        GO TO 5
C3       CONTINUE
C        LT = MAX(1,CKLSCH(ARGV))
C        WRITE(LOUT,*)'SURFTHERM ERROR: Could not open command file',
C     $               ' input on the command line, ', ARGV(1:LT)
C        RETURN
C      END IF
C*****END A-1): C-2) Command line arguments: hp
C*****A-1) unix
C      OPEN(UNIT=LCOM,   FILE=ARGV,
C     $     STATUS='OLD', FORM='FORMATTED',   ERR=5)
C*****END A-1) unix
      LCOM_EXIST = .TRUE.
      LT = MAX(1, CKLSCH(ARGV))
      WRITE(LOUT,9001) ARGV(1:LT)
5     CONTINUE
C
C-----------------------------------------------------------------------
C                              BLOCK 4
C                     OPEN CHEMKIN DATA FILES
C-----------------------------------------------------------------------
C
C               OPEN STATEMENTS
C      OPEN(UNIT=LINKCK, FILE='chem.bin',
C     $     STATUS='OLD', FORM='UNFORMATTED', ERR=1001)
C      OPEN(UNIT=LINKSK, FILE='surf.bin',
C     $     STATUS='OLD', FORM='UNFORMATTED', ERR=1002)
C
C-----------------------------------------------------------------------
C                           BLOCK 5
C         OPEN TRANSPORT PACKAGE DATA FILE IF REQUIRED
C-----------------------------------------------------------------------
C
500   CONTINUE
C      IF (LTRAN .AND. .NOT. TRAN_INIT) THEN
C        OPEN( UNIT=LINKMC, FILE='tran.bin',
C     $        STATUS='OLD', FORM='UNFORMATTED', ERR=1003)
C      END IF
C
C-----------------------------------------------------------------------
C                           BLOCK 5.5
C       INITIALIZE CHEMKIN AND QUERY FOR SIZE AND DIMENSIONS
C                 OF THE CHEMISTRY PROBLEM
C-----------------------------------------------------------------------
*
* Set Position Pointers into Work Spaces
*
      P_IW = 1
      P_RW = 1
      P_CW = 1
C
C Gas-phase chemkin
C
C        Find the length of the internal work spaces
C
      IF (.NOT. RESTRT) THEN
        CALL CKLEN  (LINKCK, LOUT, LENICK, LENRCK, LENCCK, IFLAG)
        IF ( IFLAG .NE. 0) THEN
           WRITE(LOUT,*)' QUITTING BECAUSE CKLEN IFLAG = ', IFLAG
C           STOP
           RETURN
        ENDIF
      ENDIF
C
C        Assign values to the pointers to chemkin workspaces
C
      N_CKW = P_RW
      P_RW  = N_CKW  + LENRCK
      IF (P_RW .GT. LEN_RWORK) THEN
        PRINT 2001, LEN_RWORK, P_RW
C        STOP 4
        RETURN
      END IF
      I_CKW = P_IW
      P_IW  = I_CKW  + LENICK
      IF (P_IW .GT. LEN_IWORK) THEN
        PRINT 2002, LEN_IWORK, P_IW
C        STOP 4
         RETURN
      END IF
      C_CKW = P_CW
      P_CW  = C_CKW  + LENCCK
      IF (P_CW .GT. LEN_CWORK) THEN
        PRINT 2003, LEN_CWORK, P_CW
C        STOP 4
        RETURN
      END IF
C
C         Initialize the Chemkin work spaces
C
      IF (.NOT. RESTRT) THEN
        CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT,
     $               IWORK(I_CKW), RWORK(N_CKW), CWORK(C_CKW), IFLAG)
        IF (IFLAG .NE. 0) THEN
           WRITE(LOUT,*)' QUITTING BECAUSE CKINIT IFLAG = ', IFLAG
C           STOP
           RETURN
        ENDIF
C
C         Find out the values of key length variables in mechanism
C
        CALL CKINDX (IWORK(I_CKW), RWORK(N_CKW), NELEM, KKGAS, IIGAS,
     $               NFIT)
C
C Surface Chemkin
C
C         Find out the lengths of the surface chemkin work spaces
C
        CALL SKLEN (LINKSK, LOUT, LENISK, LENRSK, LENCSK, IFLAG)
        IF ( IFLAG .NE. 0) THEN
           WRITE(LOUT,*)' QUITTING BECAUSE SKLEN IFLAG = ', IFLAG
C           STOP
           RETURN
        ENDIF
      ENDIF
C
C         Assign values to the pointers to the schemkin work space
C
      N_SKW = P_RW
      P_RW  = P_RW + LENRSK
      IF (P_RW .GT. LEN_RWORK) THEN
        PRINT 2001, LEN_RWORK, P_RW
C        STOP 4
        RETURN
      END IF
      I_SKW = P_IW
      P_IW  = P_IW + LENISK
      IF (P_IW .GT. LEN_IWORK) THEN
        PRINT 2002, LEN_IWORK, P_IW
C        STOP 4
        RETURN
      END IF
      C_SKW = P_CW
      P_CW  = P_CW + LENCSK
      IF (P_CW .GT. LEN_CWORK) THEN
        PRINT 2003, LEN_CWORK, P_CW
C        STOP 4
        RETURN
      END IF
C
C         Initialize the surface chemkin work space
C
      IF (.NOT. RESTRT) THEN
        CALL SKINIT (LENISK, LENRSK, LENCSK, LINKSK, LOUT,
     $               IWORK(I_SKW), RWORK(N_SKW), CWORK(C_SKW), IFLAG)
        IF ( IFLAG .NE. 0) THEN
           WRITE(LOUT,*)' QUITTING BECAUSE SKINIT IFLAG = ', IFLAG
C           STOP
           RETURN
        ENDIF
C
C         Find out the values of key length variables in mechanism
C
        CALL SKINDX (IWORK(I_SKW), NELEM, KKGAS, KKSURF, KKBULK,
     $               KKTOT, NNPHAS, NNSURF, NFSURF, NLSURF, NNBULK,
     $               NFBULK, NLBULK, IISUR)
      ENDIF
C
C-----------------------------------------------------------------------
C                           BLOCK 6
C         SETUP FIXED POINTERS FOR INTEGER AND LOGICAL WORK SPACE
C-----------------------------------------------------------------------
C
      I_NEL = P_IW
      I_KSG = I_NEL  +  NELEM * KKTOT
      I_KSF = I_KSG  +  KKGAS * IIGAS
      I_KSS = I_KSF  +  KKGAS * IIGAS
      I_KSU = I_KSS  +  KKTOT * IISUR
      I_NST = I_KSU  +  KKTOT * IISUR
      P_IW  = I_NST  +  IISUR * NNPHAS
*
*         Logical work space is assigned in the integer space
*
      L_THM = P_IW
      L_GRX = L_THM  +  KKTOT
      L_SRX = L_GRX  +  IIGAS
      L_STK = L_SRX  +  IISUR
      L_SCV = L_STK  +  IISUR
      L_PFL = L_SCV  +  IISUR
      L_TFL = L_PFL  +  IIGAS
      L_GTB = L_TFL  +  IIGAS
      P_IW  = L_GTB  +  IIGAS
*
C-----------------------------------------------------------------------
C                           BLOCK 7
C         SETUP FIXED POINTERS FOR REAL WORK SPACE
C-----------------------------------------------------------------------
*
      NFAR = 8
      N_FPR = P_RW
      N_AKI = N_FPR  +  NFAR  * IIGAS
      N_DJK = N_AKI  +  KKGAS * IIGAS
      N_XBT = N_DJK  +  KKGAS * KKGAS
      P_RW  = N_XBT  +  KKTOT
*
C     BLOCK 7.B Pointers extended by Fran Rupley for ReactionDesign
      C_KNAM = P_CW
      C_ENAM = C_KNAM + KKTOT
      C_PNAM = C_ENAM + NELEM
      P_CW   = C_PNAM + NNPHAS 
      IF (P_CW .GT. LEN_CWORK) THEN
        PRINT 2003, LEN_CWORK, P_CW
C        STOP 4
        RETURN
      END IF
C
      C_GNAM = 1
      C_SNAM = C_GNAM + IIGAS
      P_C64  = C_SNAM + IISUR
      IF (P_C64 .GT. LEN_CWRK64) THEN
        PRINT 3003, LEN_CWRK64, P_C64
C        STOP 4
        RETURN
      END IF
C
      P_C48  = NNPHAS+1
      IF (P_C48 .GT. LEN_CWRK48) THEN
         PRINT 3004, LEN_CWRK48, P_C48
         RETURN
      ENDIF
C
      I_KPHAS = P_IW
      I_KFIRS = I_KPHAS + NNPHAS
      I_KLAST = I_KFIRS + NNPHAS
      I_KCHRG = I_KLAST + NNPHAS
      I_KOCC  = I_KCHRG + KKTOT
      I_NELPH = I_KOCC  + KKTOT
      I_GREAC = I_NELPH + NNPHAS + 1
      I_GPROD = I_GREAC + IIGAS
      I_ITHB  = I_GPROD + IIGAS
      I_GREV  = I_ITHB  + IIGAS
      I_IFOP  = I_GREV  + IIGAS
      I_IFLO  = I_IFOP  + IIGAS
      I_KFAL  = I_IFLO  + IIGAS
      I_GDELT = I_KFAL  + IIGAS
      I_ISTFL = I_GDELT + IIGAS
      I_MDELS = I_ISTFL + IISUR
      I_SDELS = I_MDELS + IISUR
      I_BDELS = I_SDELS + IISUR
      I_NRPPS = I_BDELS + IISUR
      I_ICOV  = I_NRPPS + IISUR
      I_SREV  = I_ICOV  + IISUR
      I_GREACS= I_SREV  + IISUR
      I_SREACS= I_GREACS+ IISUR
      I_BREACS= I_SREACS+ IISUR
      I_GPRODS= I_BREACS+ IISUR
      I_SPRODS= I_GPRODS+ IISUR
      I_BPRODS= I_SPRODS+ IISUR
      I_KGASR = I_BPRODS+ IISUR
      I_KGASP = I_KGASR + IISUR
      I_NCON  = I_KGASP + IISUR
      I_KCOVI = I_NCON  + IISUR
      I_IMOTZ = I_KCOVI + KKTOT
      I_IEDP  = I_IMOTZ + IISUR
      I_IBHM  = I_IEDP  + IISUR
      I_IORD  = I_IBHM  + IISUR
      I_IYLD  = I_IORD  + IISUR
      I_NLIN  = I_IYLD  + IISUR
      L_LTIME = I_NLIN  + KKGAS
      P_IW    = L_LTIME  + NNPHAS
C
      N_WT = P_RW
      N_DEN = N_WT + KKTOT
      N_XEST = N_DEN + KKTOT
      N_EPS  = N_XEST + KKTOT
      N_SIG  = N_EPS  + KKGAS
      N_DIP  = N_SIG + KKGAS
      N_POL  = N_DIP + KKGAS
      N_ZROT = N_POL + KKGAS
      N_VIS  = N_ZROT+ KKGAS
      N_CON  = N_VIS + KKGAS
      N_HBTH = N_CON + KKGAS
      N_SBTH = N_HBTH+ KKTOT
      N_CPBTH= N_SBTH+ KKTOT
      N_GBTH = N_CPBTH+KKTOT
      N_ACT = N_GBTH + KKTOT
      N_SDEN = N_ACT + KKTOT
      N_H298 = N_SDEN + NNPHAS
      N_GEQ  = N_H298 + KKTOT
      N_GEQD = N_GEQ  + IIGAS
      N_RAG  = N_GEQD + IIGAS
      N_RBG  = N_RAG  + IIGAS
      N_REG  = N_RBG  + IIGAS
      N_RPG  = N_REG  + IIGAS
      N_DRCF = N_RPG  + IIGAS
      N_DRCR = N_DRCF + IIGAS
      N_SEQ  = N_DRCR + IIGAS
      N_SEQD = N_SEQ  + IISUR
      N_RAS  = N_SEQD + IISUR
      N_RBS  = N_RAS  + IISUR
      N_RES  = N_RBS  + IISUR
      N_RPS  = N_RES  + IISUR
      N_SDRAT= N_RPS  + IISUR
      N_CPARI= N_SDRAT+ IISUR 
      NSCOV = 3
      P_RW   = N_CPARI+ NSCOV * KKTOT
C
      
C-----------------------------------------------------------------------
C                           BLOCK 8
C         INITIALIZE TRANSPORT PACKAGE IF REQUESTED
C-----------------------------------------------------------------------
*
      IF (LTRAN) THEN
*
*           Find out the lengths of transport work spaces
*
        IF (.NOT. TRAN_INIT) THEN
          CALL MCLEN (LINKMC, LOUT, LENIMC, LENRMC, IFLAG)
          IF ( IFLAG .NE. 0) THEN
              WRITE(LOUT,*)' QUITTING BECAUSE MCLEN IFLAG = ', IFLAG
C              STOP
              RETURN
          ENDIF
        ENDIF
*
*           Assign values to pointers to transport work space
*
        N_MCK = P_RW
        P_RW  = P_RW + LENRMC
        IF (P_RW .GT. LEN_RWORK) THEN
          PRINT 2001, LEN_RWORK, P_RW
C          STOP 4
           RETURN
        END IF
        I_MCK = P_IW
        P_IW  = P_IW + LENIMC
        IF (P_IW .GT. LEN_IWORK) THEN
          PRINT 2002, LEN_IWORK, P_IW
C          STOP 4
          RETURN
        END IF
*
*           Initialize Transport work space
*
        IF (.NOT. TRAN_INIT) THEN
          CALL MCINIT (LINKMC, LOUT, LENIMC, LENRMC, IWORK(I_MCK),
     $                 RWORK(N_MCK), IFLAG)
          IF (IFLAG .NE. 0) THEN
             WRITE(LOUT,*)' QUITTING BECAUSE MCINIT IFLAG = ', IFLAG
C             STOP 4
             RETURN
          ENDIF
          TRAN_INIT = .TRUE.
        ENDIF
      ELSE
        N_MCK = 1
        I_MCK = 1
      END IF
*
C-----------------------------------------------------------------------
C                           BLOCK 9
C         SETUP OTHER POINTERS FOR INTEGER AND LOGICAL WORK SPACE
C-----------------------------------------------------------------------
*
      IF (P_IW .GT. LEN_IWORK) THEN
        PRINT 2012, LEN_IWORK, P_IW
C        STOP 4
        RETURN
      END IF
*
C-----------------------------------------------------------------------
C                           BLOCK 10
C         SETUP OTHER POINTERS FOR REAL WORK SPACE
C-----------------------------------------------------------------------
*
*         Calculate the maximum permissible size of MAX_TTABLE_NUM
*
      MAX_TTABLE_NUM = (LEN_RWORK - P_RW) /
     $                 (4 * KKTOT + 3 * IIGAS + 3 * IISUR)

      N_H__ = P_RW
      N_S__ = N_H__  +  KKTOT * MAX_TTABLE_NUM
      N_CP_ = N_S__  +  KKTOT * MAX_TTABLE_NUM
      N_G__ = N_CP_  +  KKTOT * MAX_TTABLE_NUM
      N_DSG = N_G__  +  KKTOT * MAX_TTABLE_NUM
      N_DHG = N_DSG  +  IIGAS * MAX_TTABLE_NUM
      N_DGG = N_DHG  +  IIGAS * MAX_TTABLE_NUM
      N_DSS = N_DGG  +  IIGAS * MAX_TTABLE_NUM
      N_DHS = N_DSS  +  IISUR * MAX_TTABLE_NUM
      N_DGS = N_DHS  +  IISUR * MAX_TTABLE_NUM
      P_RW  = N_DGS  +  IISUR * MAX_TTABLE_NUM
*
      IF (P_RW .GT. LEN_RWORK) THEN
        PRINT 2011, LEN_CRWORK, P_RW
C        STOP 4
        RETURN
      END IF
*
C-----------------------------------------------------------------------
C                           BLOCK 11
C                       CALL SURFTHERM
C-----------------------------------------------------------------------
*
      CALL SURFTHERM (
     1  RESTRT, LEN_IWORK, LEN_CWORK, LEN_RWORK, MAX_TTABLE_NUM, ARGV,
     3  RWORK(N_CKW), RWORK(N_SKW), RWORK(N_MCK), RWORK(N_FPR),
     4  RWORK(N_AKI), RWORK(N_H__), RWORK(N_S__), RWORK(N_CP_),
     5  RWORK(N_G__), RWORK(N_DSG), RWORK(N_DHG), RWORK(N_DGG),
     6  RWORK(N_DSS), RWORK(N_DHS), RWORK(N_DGS), RWORK(N_DJK),
     7  RWORK(N_XBT),
     8  IWORK(I_CKW), IWORK(I_SKW), IWORK(I_MCK), IWORK(I_NEL),
     9  IWORK(I_KSG), IWORK(I_KSF), IWORK(I_KSS), IWORK(I_KSU),
     *  IWORK(I_NST),
     1  IWORK(L_THM), IWORK(L_GRX), IWORK(L_SRX), IWORK(L_STK),
     2  IWORK(L_SCV), IWORK(L_PFL), IWORK(L_TFL), IWORK(L_GTB),
C     2  CWORK(C_CKW), CWORK(C_SKW))
     3  CWORK(C_CKW), CWORK(C_SKW), 
     4  CWORK(C_KNAM), CWORK(C_ENAM), CWORK(C_PNAM),
C       KNAME,         ENAME,         PNAME,
     5  CWRK64(C_GNAM), CWRK64(C_SNAM), CWRK48,    IWORK(I_KPHAS),
C       IGASNAME,       ISURNAME,       NELSTRING, KKPHAS,
     6  IWORK(I_KFIRS), IWORK(I_KLAST), IWORK(I_KCHRG),
C       KFIRST,         KLAST,          KCHARG,
     7  IWORK(I_KOCC), IWORK(I_NELPH), RWORK(N_WT), RWORK(N_DEN),
C       KOCC,          NELPHA,         WT,          DEN,
     8  IWORK(I_NLIN), RWORK(N_XEST), RWORK(N_EPS), RWORK(N_SIG), 
C       NLIN,          X_EST,         EPS,          SIG,
     9  RWORK(N_DIP),
C       NDIP,
     9  RWORK(N_POL), RWORK(N_ZROT), RWORK(N_VIS), RWORK(N_CON), 
C       POL,          ZROT,          VIS,          CON,
     *  RWORK(N_HBTH), RWORK(N_SBTH), RWORK(N_CPBTH), RWORK(N_GBTH), 
C       H_BATH,        S_BATH,        CP_BATH,        G_BATH,
     1  RWORK(N_ACT), RWORK(N_SDEN), RWORK(N_H298), RWORK(N_GEQ), 
C       ACT,          SDEN,          H298,          GEQKC,
     2  RWORK(N_GEQD), RWORK(N_RAG), RWORK(N_RBG), RWORK(N_REG), 
C       GEQKC_DELTA,    RA_GAS,       RB_GAS,       RE_GAS,
     3  RWORK(N_RPG), RWORK(N_DRCF), RWORK(N_DRCR), IWORK(I_GREAC), 
C       RP_GAS,       ND_RCF,        ND_RCR,        GAS_REAC_GAS,
     4  IWORK(I_GPROD), IWORK(I_ITHB), IWORK(I_GREV), IWORK(I_IFOP), 
C       GAS_PROD_GAS,   ITHB,          IREV_GAS,      IFOP,
     5  IWORK(I_IFLO), IWORK(I_KFAL), IWORK(I_GDELT),   
C       IFLO,          KFAL,          GAS_MOLE_CHANGE_GAS,
     6  RWORK(N_SEQ), RWORK(N_SEQD), RWORK(N_RAS), RWORK(N_RBS),
C       SEQKC,        SEQKC_DELTA,   RA_SUR,       RB_SUR,
     7  RWORK(N_RES), RWORK(N_RPS), RWORK(N_SDRAT), NSCOV, 
C       RE_SUR,       RP_SUR,       SDEN_RATIO,     NSCOV,
     8  RWORK(N_CPARI), IWORK(I_ISTFL), IWORK(I_MDELS),
C       CPARI,          ISTFL,          GAS_MOLE_CHANGE_SUR,
     9  IWORK(I_SDELS),      IWORK(I_BDELS),       IWORK(I_NRPPS), 
C       SUR_MOLE_CHANGE_SUR, BULK_MOLE_CHANGE_SUR, NRPP_SUR,
     *  IWORK(I_ICOV), IWORK(I_SREV), IWORK(I_GREACS),  
C       ICOV,          IREV_SUR,      GAS_REAC,
     1  IWORK(I_SREACS),  IWORK(I_BREACS), IWORK(I_GPRODS),  
C       SUR_REAC,         BULK_REAC,       GAS_PROD,
     2  IWORK(I_SPRODS), IWORK(I_BPRODS), IWORK(I_KGASR), 
C       SUR_PROD,        BULK_PROD,       K_GAS_REAC,
     3  IWORK(I_KGASP), IWORK(I_NCON), IWORK(I_KCOVI), 
C       K_GAS_PROD,     NCON,          KCOVI,
     4  IWORK(I_IMOTZ), IWORK(I_IEDP), IWORK(I_IBHM),
C       IMOTZ,          IEDP,          IBHM,
     5  IWORK(I_IORD), IWORK(I_IYLD), IWORK(L_LTIME))
C       IORD,          IYLD,          LESTIM
C
C-----------------------------------------------------------------------
C                           BLOCK 12
C                CHECK FOR RESTART CONDITIONS
C-----------------------------------------------------------------------
*
      IF (RESTRT) GO TO 500
*
C      STOP
      RETURN
*
C-----------------------------------------------------------------------
C                           BLOCK 13
C                        ERROR HANDLING
C-----------------------------------------------------------------------
*
C1001  WRITE(LOUT,9110)
C      STOP
C1002  WRITE(LOUT,9111)
C      STOP
C1003  WRITE(LOUT,9112)
C      STOP
*
C----------------------------------------------------------------------C
*
2001  FORMAT(' SURFTHERM ERROR: The length of the real work space,'
     1   ,I6,', is too small to even load the CHEMKIN data arrays.'/
     2    15X,'Increase size to at least ',I6,' and probably much more')
2002  FORMAT(' SURFTHERM ERROR: The length of the integer work space,'
     1   ,I6,', is too small to even load the CHEMKIN data arrays.'/
     2    15X,'Increase size to at least ',I6,' and probably much more')
2003  FORMAT(' SURFTHERM ERROR: The length of the character work space,'
     1   ,I6,', is too small to even load the CHEMKIN data arrays.'/
     2    15X,'Increase size to at least ',I6,' and probably much more')
3003  FORMAT(' SURFTHERM ERROR: The length of the CHARACTER*64 work ',
     1    'space,',I6,', is too small to load the reaction arrays.'/
     2    15X,'Increase size to at least ',I6)
3004  FORMAT(' SURFTHERM ERROR: The length of the CHARACTER*48 work ',
     1    'space,',I6,', is too small to load the phase array.'/
     2    15X,'Increase size to at least ',I6)
2011  FORMAT(' SURFTHERM ERROR: The length of the real work space,'
     1   ,I6,', is too small '/ 15X,'Increase its size to ',I6)
2012  FORMAT(' SURFTHERM ERROR: The length of the integer work space,'
     1   ,I6,', is too small '/ 15X,'Increase its size to ',I6)
9001  FORMAT (12X,'Command File = ', A)
9110  FORMAT (' SURFTHERM ERROR: Error openning chemkin linking file'/
     1        '                  Probably couldn''t find chem.bin in ',
     2        'current direction')
9111  FORMAT (' SURFTHERM ERROR: Error openning surface chemkin ',
     1        'linking file'/
     1        '                  Probably couldn''t find surf.bin in ',
     2        'current direction')
9112  FORMAT (' SURFTHERM ERROR: Error openning transport linking file'/
     1        '                  Probably couldn''t find tran.bin in ',
     2        'current direction')
      END
C----------------------------------------------------------------------C
      SUBROUTINE SURFTHERM (
     1   RESTRT,      LEN_IWORK,   LEN_CWORK,   LEN_RWORK,
     2   MAX_TTABLE_NUM, ARGV,
     3   CKWRK,       SKWRK,       RMCWRK,      FPAR,
     4   AKI,         H,           S,           CP,
     5   G,           DELTA_S_GAS, DELTA_H_GAS, DELTA_G_GAS,
     6   DELTA_S_SUR, DELTA_H_SUR, DELTA_G_SUR, DJK,
     $   X_BATH,
     7   ICKWRK,      ISKWRK,      IMCWRK,      NEL,
     8   KSTOIC_GAS,  KSTOICF_GAS, KSTOICF_SUR, KSTOIC_SUR,
     9   NSTOIC,
     *   LTHRM,       LGRXN,       LSRXN,       LSTCK,
     1   LSCOV,       LPFAL,       LTFAL,       LGTHB,
     2   CCKWRK,      CSKWRK,
C     3  CWORK(C_KNAM), CWORK(C_ENAM), CWORK(C_PNAM),
     3  KNAME,         ENAME,         PNAME,
C     4  CWRK64(C_GNAM), CWRK64(C_SNAM), CWRK48,    IWORK(I_KPHAS),
     4  IGASNAME,       ISURNAME,       NELSTRING, KKPHAS,
C     5  IWORK(I_KFIRS), IWORK(I_KLAST), IWORK(I_KCHRG),
     5  KFIRST,         KLAST,          KCHARG,
C     6  IWORK(I_KOCC), IWORK(I_NELPH), RWORK(N_WT), RWORK(N_DEN),
     6  KOCC,          NELPHA,         WT,           DEN,
C     7  IWORK(I_NLIN), RWORK(N_XEST), RWORK(N_EPS), RWORK(N_SIG),
     7  NLIN,          X_EST,         EPS,          SIG,
C     8  RWORK(N_DIP),
     8  DIP,
C     8  RWORK(N_POL), RWORK(N_ZROT), RWORK(N_VIS), RWORK(N_CON),
     8  POL,          ZROT,          VIS,          CON,
C     9  RWORK(N_HBTH), RWORK(N_SBTH), RWORK(N_CPBTH), RWORK(N_GBTH),
     9  H_BATH,        S_BATH,        CP_BATH,        G_BATH,
C     *  RWORK(N_ACT), RWORK(N_SDEN), RWORK(N_H298), RWORK(N_GEQ),
     *  ACT,          SDEN,          H298,          GEQKC,
C     1  RWORK(N_GEQD), RWORK(N_RAG), RWORK(N_RBG), RWORK(N_REG),
     1  GEQKC_DELTA,    RA_GAS,       RB_GAS,       RE_GAS,
C     2  RWORK(N_RPG), RWORK(N_DRCF), RWORK(N_DRCR), IWORK(I_GREAC),
     2  RP_GAS,       ND_RCF,        ND_RCR,        GAS_REAC_GAS,
C     3  IWORK(I_GPROD), IWORK(I_ITHB), IWORK(I_GREV), IWORK(I_IFOP),
     3  GAS_PROD_GAS,   ITHB,          IREV_GAS,      IFOP,
C     4  IWORK(I_IFLO), IWORK(I_KFAL), IWORK(I_GDELT),
     4  IFLO,          KFAL,          GAS_MOLE_CHANGE_GAS,
C     5  RWORK(N_SEQ), RWORK(N_SEQD), RWORK(N_RAS), RWORK(N_RBS),
     5  SEQKC,        SEQKC_DELTA,   RA_SUR,       RB_SUR,
C     6  RWORK(N_RES), RWORK(N_RPS), RWORK(N_SDRAT),
     6  RE_SUR,       RP_SUR,       SDEN_RATIO,          NSCOV,
C     7  RWORK(N_CPARI), IWORK(I_ISTFL), IWORK(I_MDELS), NSCOV,
     7  CPARI,          ISTFL_SUR,          GAS_MOLE_CHANGE_SUR,
C     8  IWORK(I_SDELS),      IWORK(I_BDELS),       IWORK(I_NRPPS),
     8  SUR_MOLE_CHANGE_SUR, BULK_MOLE_CHANGE_SUR, NRPP_SUR,
C     9  IWORK(I_ICOV), IWORK(I_SREV), IWORK(I_GREACS),
     9  ICOV_SUR,          IREV_SUR,      GAS_REAC,
C     *  IWORK(I_SREACS),  IWORK(I_BREACS), IWORK(I_GPRODS),
     *  SUR_REAC,         BULK_REAC,       GAS_PROD,
C     1  IWORK(I_SPRODS), IWORK(I_BPRODS), IWORK(I_KGASR),
     1  SUR_PROD,        BULK_PROD,       K_GAS_REAC,
C     2  IWORK(I_KGASP), IWORK(I_NCON), IWORK(I_KCOVI),
     2  K_GAS_PROD,     NCON,          KCOVI,
C     3  IWORK(I_IMOTZ), IWORK(I_IEDP), IWORK(I_IBHM),
     3  IMOTZ_SUR,          IEDP_SUR,          IBHM_SUR,
C     4  IWORK(I_IORD), IWORK(I_IYLD), IWORK(L_TIME))
     4  IORD_SUR,          IYLD_SUR,          LESTIM )
C
C----------------------------------------------------------------------C
      INCLUDE 'surftherm.h'
C----------------------------------------------------------------------C
C
C Dummy Variables
C
      LOGICAL   RESTRT
      INTEGER   LEN_IWORK, LEN_CWORK, LEN_RWORK, MAX_TTABLE_NUM,
     $          ICKWRK(LENICK), ISKWRK(LENISK), IMCWRK(LENIMC),
     $          NEL(NELEM,KKTOT),  KSTOIC_GAS(KKGAS,IIGAS),
     $          KSTOICF_GAS(KKGAS,IIGAS), KSTOICF_SUR(IISUR,KKTOT),
     $          KSTOIC_SUR(IISUR,KKTOT), NSTOIC(IISUR,NNPHAS)
      LOGICAL   LTHRM(KKTOT), LGRXN(IIGAS), LSRXN(IISUR), LSTCK(IISUR),
     $          LSCOV(IISUR), LPFAL(IIGAS), LTFAL(IIGAS), LGTHB(IIGAS)
      CHARACTER CCKWRK(LENCCK)*16, CSKWRK(LENCSK)*16, ARGV*(*)
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1      CKWRK(LENRCK), SKWRK(LENRSK), RMCWRK(LENRMC),
     2      FPAR(NFAR,IIGAS), AKI(KKGAS,IIGAS),
     3      H( KKTOT,MAX_TTABLE_NUM),  S(KKTOT,MAX_TTABLE_NUM),
     4      CP(KKTOT,MAX_TTABLE_NUM),  G(KKTOT,MAX_TTABLE_NUM),
     5      DELTA_S_GAS(IIGAS,MAX_TTABLE_NUM),
     6      DELTA_H_GAS(IIGAS,MAX_TTABLE_NUM),
     7      DELTA_G_GAS(IIGAS,MAX_TTABLE_NUM),
     8      DELTA_S_SUR(IISUR,MAX_TTABLE_NUM),
     9      DELTA_H_SUR(IISUR,MAX_TTABLE_NUM),
     *      DELTA_G_SUR(IISUR,MAX_TTABLE_NUM), DJK(KKGAS,KKGAS),
     1      X_BATH(KKTOT)
C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
C Variables for the Chemistry and Species set-up
C
C      INTEGER     MAXBULK, MAXSUR, MAXGAS, MAXTOT, MAXPHA, MAXELEM
C*****C-1) static memory allocation
C      PARAMETER (MAXBULK=4,  MAXSUR=60, MAXGAS=50,
C     $           MAXTOT=MAXGAS+MAXSUR+MAXBULK, MAXPHA=8, MAXELEM=8)
C*****END C-1) static memory allocation
C*****C-2) dynamic memory allocation
C      EQUIVALENCE (KKBULK, MAXBULK), (KKSURF, MAXSUR),
C     $            (KKGAS,  MAXGAS),  (KKTOT,  MAXTOT),
C     $            (NNPHAS, MAXPHA),  (NELEM,  MAXELEM)
C*****END C-2) dynamic memory allocation
*
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1   WT(KKTOT), ACT(KKTOT), SDEN(NNPHAS),  SDTOT, H298(KKTOT)
      INTEGER
     1   KKPHAS(NNPHAS), KFIRST(NNPHAS), KLAST(NNPHAS),  KCHARG(KKTOT),
     2   KOCC(KKTOT), NELPHA(0:NNPHAS)
      CHARACTER*16 KNAME(KKTOT), PNAME(NNPHAS),  ENAME(NELEM)
      CHARACTER    NELSTRING(0:NNPHAS)*48, RATE_UNITS(-1:3)*22
      SAVE         RATE_UNITS
C
C Variables for the Gas-phase Reactions
C
C      INTEGER    MAXIIGAS
C*****C-1) static memory allocation
C      PARAMETER (MAXIIGAS=200)
C*****END C-1) static memory allocation
C*****C-2) dynamic memory allocation
C      EQUIVALENCE (IIGAS, MAXIIGAS)
C*****END C-2) dynamic memory allocation
*
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1   GEQKC(IIGAS), GEQKC_DELTA(IIGAS),
     2   RA_GAS(IIGAS), RB_GAS(IIGAS), RE_GAS(IIGAS),
     3   RP_GAS(IIGAS), ND_RCF(IIGAS), ND_RCR(IIGAS),
     $   CTB, RKLOW, PR, FC, RCFINF, PCOR
      INTEGER
     2   GAS_REAC_GAS(IIGAS), GAS_PROD_GAS(IIGAS),
     3   ITHB(IIGAS), IREV_GAS(IIGAS),IFOP(IIGAS),
     $   IFLO(IIGAS), KFAL(IIGAS),
     $   GAS_MOLE_CHANGE_GAS(IIGAS)
      CHARACTER IGASNAME(IIGAS)*64
C
C Variables for the Surface Reactions
C
      INTEGER    NSCOV
C      PARAMETER (NSCOV=3)
C      INTEGER   MAXIISUR
C*****C-1) static memory allocation
C      PARAMETER (IISUR=125)
C*****END C-1) static memory allocation
C*****C-2) dynamic memory allocation
C      EQUIVALENCE (IISUR, IISUR)
C*****END C-2) dynamic memory allocation
*
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1   SEQKC(IISUR),  SEQKC_DELTA(IISUR),
     2   RA_SUR(IISUR), RB_SUR(IISUR), RE_SUR(IISUR),
     3   RP_SUR(IISUR), RA_STK, RB_STK, RE_STK, STK,
     4   VELOC_CRT, EFF_FLUX, SDEN_RATIO(IISUR), CPARI(NSCOV,KKTOT)
      INTEGER
     1   GAS_MOLE_CHANGE_SUR(IISUR),  SUR_MOLE_CHANGE_SUR(IISUR),
     2   BULK_MOLE_CHANGE_SUR(IISUR), SUR_COV_CHANGE,
     3   NRPP_SUR(IISUR),   ICOV_SUR(IISUR),  IREV_SUR(IISUR),
     4   ISTFL_SUR(IISUR),  IMOTZ_SUR(IISUR), IEDP_SUR(IISUR),
     5   IBHM_SUR(IISUR),   IORD_SUR(IISUR),  IYLD_SUR(IISUR),
     6   GAS_REAC(IISUR),   SUR_REAC(IISUR),  BULK_REAC(IISUR),
     7   GAS_PROD(IISUR),   SUR_PROD(IISUR),  BULK_PROD(IISUR),
     8   K_GAS_REAC(IISUR), K_GAS_PROD(IISUR),
     9   NCON(NNPHAS),  NCOVI, KCOVI(KKTOT)
      CHARACTER
     1   ISURNAME(IISUR)*64, SUR_FUNITS*24, SUR_RUNITS*24
C
C Variables for the "Solution"
C
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1   TTABLE(MAXTEMP), DEN(KKTOT)
C
C Variables associated with the "Bath Gas" and Transport Calculations
C
      INTEGER  NLIN(KKGAS), ID_T_BATH
C
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1   X_EST(KKTOT),
     1   EPS(KKGAS), SIG(KKGAS),  DIP(KKGAS),
     1   POL(KKGAS), ZROT(KKGAS), VIS(KKGAS), CON(KKGAS),
     1   H_BATH(KKTOT), S_BATH(KKTOT), CP_BATH(KKTOT),
     1   G_BATH(KKTOT), PTORR
C
C Variables for the keyword interface
C
      LOGICAL LCNTUE, CNTNUD, LGEN, LNDIM_GAS, LNDIM_SUR,
     $        LTSUM_SPECIES, LTSUM_GAS, LTSUM_SUR
      SAVE    LCNTUE, CNTNUD, LGEN, LNDIM_GAS, LNDIM_SUR,
     $        LTSUM_SPECIES, LTSUM_GAS, LTSUM_SUR
      LOGICAL LESTIM(NNPHAS)
C
C Common block for variables associated with the keyword interface
C
      LOGICAL         LTRAN
      COMMON /SF_KEY/ LTRAN
C
C Common block for Unit numbers for I/O
C
      INTEGER         LIN, LOUT, LCOM, LINKCK, LINKSK, LINKMC
      LOGICAL         LCOM_EXIST
      COMMON /SF_UIO/ LIN, LOUT, LCOM, LINKCK, LINKSK, LINKMC,
     1                LCOM_EXIST
C
C Misc. Variables
C
      INTEGER   ITEMP, LT, LT2, I, K, N, M, IPHASE, LINES, BRKLEN(4),
     1          F_RATE_UNITS, R_RATE_UNITS
      LOGICAL   IERR, FLAG, FLAG2
      CHARACTER BROKLINE(2)*40, BRKSTRNG(4)*2, STRING*80
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    WORK(MAXTEMP), T, T298, DELTAT, A_FACTOR, RF,
     2    RF_DELTA, EACT, T_OLD, RF_NEG1, RF_NEG1_DELTA, P,
     3    EACT_NEG1, A_FACTOR_NEG1, RAR, RBR, RER, RA_REV, RB_REV,
     4    RE_REV, RA_REV_ST, RB_REV_ST, RE_REV_ST, RA_GEN, RB_GEN,
     5    RE_GEN, F_RATE_SUR_UNITS, r_rate_sur, FAC_DELTA,
     6    A_FACTOR2, EACT2, FAC, RF_PRIME_NEG1, RF_PRIME, K_STAR,
     7    K_STAR_NEG1, GAS_DA_FACTOR, SURF_DA_FACTOR, F_DBLE,
     8    R_DBLE
C
C Local Variables
C
      LOGICAL TRAN_OPENED
C
C
C Externals
C
      INTEGER  PHASEID
      EXTERNAL BRKLINE, CALCRSTA, GENTOSTK, HIGH_PRESS, P_REAC,
     1         PHASEID, REVGEN, REVSTK, SRX_UNIT, STKACT,
     1         STKCALCA, STKTOGEN, NDCALC_SUR, COV_DEPA
      INTEGER  CKLSCH
      EXTERNAL CKLSCH
      EXTERNAL CKABE,  CKEQXP, CKFAL, CKINDX, CKITR,
     1         CKNU, CKNUF, CKRP, CKRDEX, CKSYMR, CKTHB
      EXTERNAL SKABE, SKCOV,  SKCPOR, SKDEN, SKEQ, SKFLGS, SKHORT,
     1         SKICOV, SKINDX, SKINIT, SKIREV, SKNCF, SKNCON, SKNU,
     1         SKNUF, SKPKK, SKRDEX, SKSDEN, SKSOR, SKSYME, SKSYMP,
     1         SKSYMR, SKSYMS, SKWT
      EXTERNAL MCLEN, MCINIT, MCPRAM, MCSCON, MCSDIF
C
C-----------------------------------------------------------------------
C
      SAVE BRKSTRNG, BRKLEN, DELTAT
      DATA BRKSTRNG /'<=', '=>', '=', '+'/
      DATA BRKLEN   / 2, 2, 1, 1/
      DATA DELTAT /1.0E-3/
      DATA CNTNUD/.FALSE./
      DATA RATE_UNITS/' (mole/(cm**3 sec))  ', '       (sec-1)       ',
     $                ' (cm**3/(mole sec))  ', '(cm**6/(mole**2 sec))',
     $                '(cm**9/(mole**3 sec))'/
C
C----------------------------------------------------------------------C
C
C-----------------------------------------------------------------------
C                               BLOCK 1
C                     INITIALIZATION OF CHEMKIN
C-----------------------------------------------------------------------
*
* Get gas constant and check sizes of chemkin static dimensions
*
      CALL CKRP (ICKWRK, CKWRK, RU, RUC, PATM)
*
* Calculate the gas constant in terms of kcal/moleK
*
      RKCAL = RUC * 1.0E-3
*
C*****C-1) static memory allocation
C      IF (KKGAS .GT. MAXGAS) THEN
C        PRINT *,'KKGAS is greater than MAXGAS,', KKGAS, MAXGAS
C        STOP
C      END IF
C      IF (IIGAS .GT. MAXIIGAS) THEN
C        PRINT *,'Number of gas reactions is greater than '
C     $         ,'storage dimensions'
C        STOP
C      END IF
C*****END C-1) static memory allocation
*
C
C Get the names of the gas-phase reactions
C
      DO 80 I = 1, IIGAS
        CALL CKSYMR(I, LOUT, ICKWRK, CKWRK, CCKWRK, LT, IGASNAME(I),
     1              IERR)
        IF (IERR) THEN
          PRINT *,'Error flag set in call to CKSYMR'
C          STOP
        RETURN
        END IF
80    CONTINUE
C
C Get the third body and reversiblity flags for gas-phase reactions
      CALL CKITR (ICKWRK, CKWRK, ITHB, IREV_GAS)
C Get the enhanced third body efficiencies for the gas-phase reactions
      CALL CKTHB (KKGAS, ICKWRK, CKWRK, AKI)
C Get the net Stoichiometric coefficients arrays for gas reactions
      CALL CKNU  (KKGAS, ICKWRK, CKWRK, KSTOIC_GAS)
C Get the forward Stoichiometric coefficients arrays for gas reactions
      CALL CKNUF (KKGAS, ICKWRK, CKWRK, KSTOICF_GAS)
C Get the Arrhenius parameters for the gas reactions - change to
      CALL CKABE (ICKWRK, CKWRK, RA_GAS, RB_GAS, RE_GAS)
      DO 81 I = 1, IIGAS
        CALL CKRDEX (I, CKWRK, RP_GAS(I))
81    CONTINUE
*
C Get the fall-off parameters for the gas reactions
C  IFOP(*) - flags indicating type of pressure-dependent behavior;
C  IFLO(*) - flags indication pressure-depencency;
C              0 - unimolecular fall-off,
C              1 - chemically activated bi-molecular.
C  KFAL(*) - flags indicating type of bath-gas concentration to be used
C            in expressions
C               0 - Use total concentration of gas mixture
C                    (with the added capability of using enhanced
C                     third body coefficients) (default)
C               K - Use the concentration of species K
C  FPAR(*,*) - Pressure dependency parameters;
*
      CALL CKFAL (NFAR, ICKWRK, CKWRK, IFOP, IFLO, KFAL, FPAR)
*
C-----------------------------------------------------------------------
C                                BLOCK 1.5
C                 INITIALIZATION OF SURFACE CHEMKIN
C-----------------------------------------------------------------------
*
C*****C-1) static memory allocation
C      IF (NELEM .GT. MAXELEM) THEN
C        PRINT *,'NELEM is greater than MAXELEM,', NELEM, MAXELEM
C        STOP
C      END IF
C      IF (KKSURF .GT. MAXSUR) THEN
C        PRINT *,'KKSURF is greater than MAXSUR,', KKSURF, MAXSUR
C        STOP
C      END IF
C      IF (KKBULK .GT. MAXBULK) THEN
C        PRINT *,'KKBULK is greater than MAXBULK,', KKBULK, MAXBULK
C        STOP
C      END IF
C      IF (KKTOT .GT. MAXTOT) THEN
C        PRINT *,'KKTOT is greater than MAXTOT,', KKTOT, MAXTOT
C        STOP
C      END IF
C      IF (NNPHAS .GT. NNPHAS) THEN
C        PRINT *,'NNPHAS is greater than MAXPHA,', NNPHAS, MAXPHA
C        STOP
C      END IF
C      IF (IISUR .GT. MAXIISUR) THEN
C        PRINT *,'IISUR is greater than MAXIISUR,', IISUR, MAXIISUR
C        STOP
C      END IF
C*****END C-1) static memory allocation
*
C Get Phase pointers
      CALL SKPKK  (ISKWRK, KKPHAS, KFIRST, KLAST)
C Get molecular weights for all species
      CALL SKWT   (ISKWRK, SKWRK, WT)
C Get the elemental compositions of all species
      CALL SKNCF(NELEM, ISKWRK, NEL)
C Get the site occupancy numbers for all species
      CALL SKCOV (ISKWRK, KOCC)
C Get Character strings for the element names
      CALL SKSYME (ISKWRK, CSKWRK, LOUT, ENAME, IERR)
      IF (IERR) THEN
        PRINT *,'Error flag set in call to SKSYME'
C        STOP
        RETURN
      END IF
C Get Character strings for the surface species
      CALL SKSYMS (ISKWRK, CSKWRK, LOUT, KNAME, IERR)
      IF (IERR) THEN
        PRINT *,'Error flag set in call to SKSYMS'
C        STOP
        RETURN
      END IF
C Get Character strings for the phase names
      CALL SKSYMP (ISKWRK, CSKWRK, LOUT, PNAME, IERR)
      IF (IERR) THEN
        PRINT *,'Error flag set in call to SKSYMP'
C        STOP
        RETURN
      END IF
C Get Character strings for the surface phase reactions
      DO 82 I = 1, IISUR
        CALL SKSYMR(I, LOUT, ISKWRK, SKWRK, CSKWRK, LT, ISURNAME(I),
     $              IERR)
        IF (IERR) THEN
          PRINT *,'Error flag set in call to CKSYMR'
C          STOP
          RETURN
        END IF
82    CONTINUE
C Get Standard state surface site densities
      CALL SKSDEN (ISKWRK, SKWRK, SDEN)
C Get Charges for the species
      CALL SKCHRG(ISKWRK, SKWRK, KCHARG)
C Get Stoichiometric coefficients arrays for surface reactions
      CALL SKNU  (IISUR, ISKWRK, SKWRK, KSTOIC_SUR, NSTOIC)
C Get the forward Stoichiometric coefficients arrays for surf reactions
      CALL SKNUF (IISUR, ISKWRK, KSTOICF_SUR)
C Get the total number of surface reactions which don't conserve sites
      CALL SKNCON (ISKWRK, SKWRK, NCON)
C Get Arrhenius parameters for surface reactions
      CALL SKABE (ISKWRK, SKWRK, RA_SUR, RB_SUR, RE_SUR, ISTFL_SUR)
      DO 83 I = 1, IISUR
        CALL SKRDEX(I, ISKWRK, SKWRK, RP_SUR(I))
83    CONTINUE
*
*        Get flags describing types of surface reactions:
*
C  NRPP_SUR  - Integer scalar, number of species (reactants+products)
C              for surface reaction IR, combined with reversibility
C              flag.
C              NRPP > 0, NRPP species, reversible surface reaction,
C                  < 0, ABS(NRPP) species, irreversible reaction.
C  IREV_SUR  - Integer scalar, flag for explicit reverse Arrhenius
C              parameters.
C              =1, reaction has explicit reverse Arrhenius parameters
C              =0, no (may or may not be reversible, see NRPP).
C  ISTFL_SUR - Integer scalar, flag for sticking coefficients;
C              =1, reaction does not use sticking coefficients
C              =0, no
C  IMOTZ_SUR - Integer scalar, flag for Motz-Wise correction of
C              sticking coefficients;
C              =1, sticking reaction with Motz-Wise correction
C              =0, no (may or may not be sticking reaction, see ISTFL)
C  ICOV_SUR  - Integer scalar, flag to indidicate that reaction has
C              coverage dependence;
C              =1, reaction has coverage dependence
C              =0, no.
C  IEDP_SUR  - Integer scalar, flag for energy-dependence;
C              =1, reaction is energy-dependent,
C              =0, no.
C  IBHM_SUR  - Integer scalar, flag for Bohm correction;
C              =1, Bohm reaction,
C              =0, no
C  IORD_SUR  - Integer scalar, flag for species order change;
C              =1, reaction has species order change,
C              =0, no
C  IYLD_SUR  - Integer scalar, flag for yield-modification;
C              =1, yield-modification in reaction;
C              =0, no
*
      DO 84 I = 1, IISUR
        CALL SKFLGS (I, ISKWRK, NRPP_SUR(I), IREV_SUR(I), ISTFL_SUR(I),
     $               ICOV_SUR(I), IMOTZ_SUR(I), IEDP_SUR(I),
     $               IBHM_SUR(I), IORD_SUR(I), IYLD_SUR(I))
84    CONTINUE
C-----------------------------------------------------------------------
C                         BLOCK 1.7
C              QUERY USER CONCERNING PREFERENCES
C      1. Calculate Heats of Formation at 298.15
C      2. Get keywords, if available
C-----------------------------------------------------------------------
*
* Calculate Heats of Formation (kcal/mole) at 298.15 K
*
      T298 = 298.15
      CALL SKHORT (T298, ISKWRK, SKWRK, H298)
      DO 90 K = 1, KKTOT
        H298(K) = H298(K) * RUC * T298 / 1000.
90    CONTINUE
C
C Get the keywords if this is not a restart
C
100   CONTINUE
      IF (RESTRT) THEN
        RESTRT = .FALSE.
      ELSE
        CALL SURFKEY
     1    (LCOM_EXIST, LCOM, LIN, LOUT, LCNTUE, CNTNUD,
     1     LTHRM, LTRAN, LGRXN, LGEN, LSRXN, LNDIM_GAS, LNDIM_SUR,
     3     LSTCK, LSCOV, LPFAL, LTFAL, LGTHB, LTSUM_SPECIES,
     4     LTSUM_GAS, LTSUM_SUR, X_BATH, KNAME, KFIRST, KLAST, KKPHAS,
     9     PNAME, IGASNAME, ISURNAME, LESTIM, X_EST)
      END IF
C
C Check to see if its necessary to go back and open transport
C data file - Set RESTRT if this is true
C
      IF (LTRAN) THEN
        INQUIRE(LINKMC, OPENED=TRAN_OPENED)
        IF (.NOT. TRAN_OPENED) THEN
          RESTRT =.TRUE.
        END IF
      END IF
C
C Calculate the number of entries in tables that span temperature
C
      TTABLE_NUM = (TTABLE_MAX - TTABLE_MIN)/TTABLE_DELTAT + 3
C
      IF (TTABLE_NUM .GT. MAXTEMP) THEN
        WRITE (LOUT,*) 'WARNING: TTABLE_NUM is greater than MAXTEMP.',
     1     'It will be reduced to ', MAXTEMP
        TTABLE_NUM = MAXTEMP
      END IF
      IF (TTABLE_NUM .GT. MAX_TTABLE_NUM )  THEN
        WRITE (LOUT,*)
     1           'WARNING: TTABLE_NUM is greater than MAX_TTABLE_NUM.'
        WRITE (LOUT,*)
     $           '         Tables will be truncated!'
        WRITE (LOUT,*)
     $           '         Desired value = ', TTABLE_NUM,
     $                   ' Actual  value = ', MAX_TTABLE_NUM
        WRITE (LOUT,*)
     $           '         Next time, increase size of RWORK'
      END IF
C
C Check to see if a restart of the subroutine is necessary
C
      IF (RESTRT) RETURN
C
C-----------------------------------------------------------------------
C                                BLOCK 1.8
C                 INITIALIZATION OF TRANSPORT PACKAGE
C-----------------------------------------------------------------------
*
      IF (LTRAN .AND. .NOT. CNTNUD) THEN
        CALL MCPRAM(IMCWRK, RMCWRK, EPS, SIG, DIP, POL, ZROT, NLIN)
      END IF
*
C-----------------------------------------------------------------------
C                         BLOCK 1.9
C     CALCULATE BATH GAS QUANTITIES AND VECTOR OF ACTIVITIES
C-----------------------------------------------------------------------
*
* Calculate the bath gas pressure in cgs units
*
      P_CGS = P_BATH * PATM / 760.
*
* Calculate the bath gas concentration in mole/cm**3
*
      CTOT_BATH = P_CGS / (RU * T_BATH)
C
C Get the nondimensional thermo quantities of all molecular species
C at the bath gas conditions
C
      CALL SKHORT (T_BATH, ISKWRK, SKWRK, H_BATH)
      CALL SKSOR  (T_BATH, ISKWRK, SKWRK, S_BATH)
      CALL SKCPOR (T_BATH, ISKWRK, SKWRK, CP_BATH)
C Convert standard enthalpies to kcal/mole
      DO 110 K = 1 , KKTOT
        H_BATH(K) = H_BATH(K) * RUC * T_BATH / 1000.
        S_BATH(K) = S_BATH(K) * RUC
        CP_BATH(K)= CP_BATH(K)* RUC
        G_BATH(K) = H_BATH(K) - T_BATH * S_BATH(K) / 1000.
110   CONTINUE
C
C Right now, the activity vector is assumed to be equal to bath
C gas mole fraction vector
C
      DO 120 K = 1, KKTOT
        ACT(K) = X_BATH(K)
120   CONTINUE
C-----------------------------------------------------------------------
C                         BLOCK 2
C              PROBLEM DIMENSIONS AND PHASE INFORMATION
C-----------------------------------------------------------------------
C Calculate the NELPHA vector
      DO 188 IPHASE = 0, NNPHAS
        NELPHA(IPHASE) = 0
        NELSTRING(IPHASE) = '('
        DO 186 M = 1, NELEM
          IF (IPHASE .EQ. 0) THEN
            LT = CKLSCH(NELSTRING(IPHASE))
            LT2 =CKLSCH(ENAME(M))
            NELSTRING(IPHASE) =
     1        NELSTRING(IPHASE)(1:LT)//' '//ENAME(M)(1:LT2)//' '
            GO TO 185
          END IF
          DO 184 K = KFIRST(IPHASE), KLAST(IPHASE)
            IF (NEL(M,K) .NE. 0 ) THEN
              NELPHA(IPHASE) = NELPHA(IPHASE) + 1
              LT = CKLSCH(NELSTRING(IPHASE))
              LT2 =CKLSCH(ENAME(M))
              NELSTRING(IPHASE) =
     1          NELSTRING(IPHASE)(1:LT)//' '//ENAME(M)(1:LT2)//' '
              GO TO 185
            END IF
184       CONTINUE
185       CONTINUE
186     CONTINUE
        IF ( NELSTRING(IPHASE) .EQ. '(' ) THEN
          NELSTRING(IPHASE) = ' '
        ELSE
          LT = CKLSCH(NELSTRING(IPHASE))
          NELSTRING(IPHASE) =  NELSTRING(IPHASE)(1:LT)//' )'
        END IF
188   CONTINUE
C
      IF (LGEN) THEN
        LT = CKLSCH(NELSTRING(0))
        WRITE(LOUT,6999) NELEM, NELSTRING(0)(1:LT),
     1                   KKTOT, NNPHAS, IIGAS, IISUR, RU, RUC, PATM
        LT = CKLSCH(NELSTRING(1))
        WRITE(LOUT,7000) KKGAS, NCON(1), NELPHA(1), NELSTRING(1)(1:LT)
      END IF
C
      SDTOT = 0.0
      IF (NNSURF .GT. 0) THEN
              IF (LGEN) WRITE(LOUT,7002)
        DO 190 IPHASE = NFSURF, NLSURF
          SDTOT = SDTOT + SDEN(IPHASE)
          LT = CKLSCH(NELSTRING(IPHASE))
          IF (LGEN)
     1    WRITE(LOUT,7003) IPHASE, PNAME(IPHASE), SDEN(IPHASE),
     1                     KKPHAS(IPHASE), NELPHA(IPHASE), NCON(IPHASE),
     2                     NELSTRING(IPHASE)(1:LT)
190     CONTINUE
        IF (LGEN) WRITE(LOUT,7008) SDTOT
      ELSE
        IF (LGEN) WRITE(LOUT,7004)
      END IF
C
      IF (LGEN) THEN
        IF (NNBULK .GT. 0) THEN
          WRITE(LOUT,7005)
          DO 191 IPHASE = NFBULK, NLBULK
            LT = CKLSCH(NELSTRING(IPHASE))
            WRITE(LOUT,7006) IPHASE, PNAME(IPHASE), KKPHAS(IPHASE),
     1                       NELPHA(IPHASE), NCON(IPHASE),
     2                       NELSTRING(IPHASE)(1:LT)
191       CONTINUE
        ELSE
          WRITE(LOUT,7007)
        END IF
      END IF
C-----------------------------------------------------------------------
C                            BLOCK 3
C            SHORT DESCRIPTION OF REACTION MECHANISM
C      1. Print out information on species in the mechanism with the
C         bath gas composition.
C      2. Print out a summary of the Species thermo at
C         bath gas conditions.
C      3. Calculate mole changes for each reaction
C      4. Surface reactions:
C         a) Calculate number and type of reactants and products
C         b) Identify the gas-phase species number for later use in
C            determining effusive velocities
C         c) Calculate the surface density ratio for later use in
C            analysing sticking coefficient reactions.
C-----------------------------------------------------------------------
*
      IF (LGEN .AND. KKTOT .GT. 0) THEN
        LT  = CKLSCH(KNAME(ID_CARR))
        LT2 = CKLSCH(KNAME(ID_MAJOR))
        WRITE(LOUT,8250) P_BATH, T_BATH, KNAME(ID_CARR)(1:LT),
     1                   KNAME(ID_MAJOR)(1:LT2)
        DO 301 K = 1, KKGAS
          WRITE(LOUT,8251) K, KNAME(K), X_BATH(K),
     1                     X_BATH(K)*(P_BATH/760.)/(82.05*T_BATH)
301     CONTINUE
        IF (NNSURF .GT. 0) THEN
          DO 302 IPHASE = NFSURF, NLSURF
            DO 302 K = KFIRST(IPHASE), KLAST(IPHASE)
              WRITE(LOUT,8252) K, KNAME(K), X_BATH(K),
     1                         X_BATH(K)*SDEN(IPHASE)/KOCC(K)
302       CONTINUE
        END IF
        IF (NNBULK .GT. 0) THEN
          DO 303 IPHASE = NFBULK, NLBULK
            DO 303 K = KFIRST(IPHASE), KLAST(IPHASE)
              WRITE(LOUT,8253) K, KNAME(K), X_BATH(K)
303       CONTINUE
        END IF
      END IF
C
C  Print out summary of thermodynamic functions for the species
C     at bath gas conditions
C
      IF (LTSUM_SPECIES .AND. KKTOT .GT. 0) THEN
        WRITE(LOUT,8260) T_BATH
        DO 311 K = 1, KKTOT
          WRITE(LOUT,8261) K, KNAME(K), H298(K), H_BATH(K), CP_BATH(K),
     1                     S_BATH(K), G_BATH(K)
311     CONTINUE
        WRITE(LOUT,8262)
      END IF
C
      IF (IIGAS .GT. 0) THEN
        IF (LGEN) WRITE(LOUT,8200)
        DO 350 I = 1, IIGAS
          GAS_MOLE_CHANGE_GAS(I) = 0
          GAS_REAC_GAS(I)        = 0
          DO 340 K = 1, KKGAS
            GAS_MOLE_CHANGE_GAS(I) = GAS_MOLE_CHANGE_GAS(I)
     1                             + KSTOIC_GAS(K,I)
            GAS_REAC_GAS(I) = GAS_REAC_GAS(I) - KSTOICF_GAS(K,I)
340       CONTINUE
          GAS_PROD_GAS(I) = GAS_REAC_GAS(I) + GAS_MOLE_CHANGE_GAS(I)
          IF (LGEN) THEN
            CALL BRKLINE(IGASNAME(I), BROKLINE, 2, 40, BRKSTRNG,
     1                   BRKLEN, 4, LINES)
            WRITE(LOUT,8201) I, BROKLINE(1), GAS_MOLE_CHANGE_GAS(I),
     1                       GAS_REAC_GAS(I)
            IF (LINES .GT. 1) WRITE(LOUT,8205) BROKLINE(2)
          END IF
350     CONTINUE
        IF (LGEN) WRITE(LOUT,8204)
      END IF
      IF (IISUR .GT. 0) THEN
      IF (LGEN) WRITE(LOUT,8202)
      DO 390 I = 1, IISUR
        GAS_MOLE_CHANGE_SUR(I) = 0
        GAS_REAC(I) = 0
        K_GAS_REAC(I) = 1
        K_GAS_PROD(I) = 1
        DO 360 K = 1, KKGAS
          GAS_MOLE_CHANGE_SUR(I) = GAS_MOLE_CHANGE_SUR(I)
     1                           + KSTOIC_SUR(I,K)
          IF (KSTOICF_SUR(I,K) .LT. 0) THEN
            GAS_REAC(I) = GAS_REAC(I) - KSTOICF_SUR(I,K)
            K_GAS_REAC(I) = K
          END IF
          IF (KSTOIC_SUR(I,K) .GT. 0) K_GAS_PROD(I) = K
360     CONTINUE
        SUR_MOLE_CHANGE_SUR(I) = 0
        SUR_REAC(I) = 0
        SDEN_RATIO(I) = 1.0
        SUR_COV_CHANGE = 0
        IF (KKSURF .GT. 0) THEN
        DO 370 IPHASE = NFSURF, NLSURF
        DO 370 K = KFIRST(IPHASE), KLAST(IPHASE)
          SUR_MOLE_CHANGE_SUR(I) = SUR_MOLE_CHANGE_SUR(I)
     1                           + KSTOIC_SUR(I,K)
          SUR_COV_CHANGE = SUR_COV_CHANGE
     1                    + KSTOIC_SUR(I,K) * KOCC(K)
          IF (KSTOICF_SUR(I,K) .LT. 0) THEN
            SUR_REAC(I) = SUR_REAC(I) - KSTOICF_SUR(I,K)
            SDEN_RATIO(I) = SDEN_RATIO(I) * (SDEN(IPHASE)/SDTOT)
          END IF
370     CONTINUE
        END IF
        BULK_MOLE_CHANGE_SUR(I) = 0
        BULK_REAC(I) = 0
        IF (KKBULK .GT. 0) THEN
        DO 380 K = KFIRST(NFBULK), KLAST(NLBULK)
          BULK_MOLE_CHANGE_SUR(I) = BULK_MOLE_CHANGE_SUR(I)
     1                           + KSTOIC_SUR(I,K)
          IF (KSTOICF_SUR(I,K) .LT. 0) THEN
            BULK_REAC(I) = BULK_REAC(I) - KSTOICF_SUR(I,K)
          END IF
380     CONTINUE
        END IF
        GAS_PROD(I)  = GAS_REAC(I)  + GAS_MOLE_CHANGE_SUR(I)
        SUR_PROD(I)  = SUR_REAC(I)  + SUR_MOLE_CHANGE_SUR(I)
        BULK_PROD(I) = BULK_REAC(I) + BULK_MOLE_CHANGE_SUR(I)
        IF (LGEN) THEN
          CALL BRKLINE(ISURNAME(I), BROKLINE, 2, 40, BRKSTRNG,
     1                    BRKLEN, 4, LINES)
          WRITE(LOUT,8203) I, BROKLINE(1), GAS_MOLE_CHANGE_SUR(I),
     1              SUR_MOLE_CHANGE_SUR(I), BULK_MOLE_CHANGE_SUR(I),
     1              SUR_COV_CHANGE
          IF (LINES .GT. 1) WRITE(LOUT,8205) BROKLINE(2)
        END IF
390   CONTINUE
      IF (LGEN) WRITE(LOUT,8206)
      END IF
C-----------------------------------------------------------------------
C                         BLOCK 3.2
C     TABLE OF NON-DIMENSIONAL GAS REACTION RATE CONSTANTS
C                 AT THE BATH GAS CONDITIONS
C-----------------------------------------------------------------------
      IF (LNDIM_GAS .AND. IIGAS .GT. 0) THEN
C
C Go get the matrix of binary diffusion coefficients
C (If the transport package is not used, calculate a fake diffusion
C  coefficient, using a reasonable value.  It is needed later for
C  the calculation of the Damkoeler numbers
C
        IF (LTRAN) THEN
          CALL MCSDIF((P_BATH/760.*PATM), T_BATH, KKGAS, RMCWRK, DJK)
        ELSE
          DJK(ID_CARR, ID_MAJOR) = 1.241E-5 * T_BATH**1.707
     1                                      * (760./P_BATH)
        END IF
C
C Calculate the nondimensionalizing factor used in the calculation
C of volumetric Damkoeler numbers.  This has units of mole/(cm**3*sec).
C
        GAS_DA_FACTOR  =  CTOT_BATH*DJK(ID_CARR, ID_MAJOR)
     $                  / LENGTH_SCALE**2
        WRITE(LOUT,8349) P_BATH, T_BATH
        CALL NDCALC_GAS (ND_RCF, ND_RCR, GAS_REAC_GAS,
     1                   P_CGS, T_BATH, X_BATH, ICKWRK, CKWRK,
     2                   GAS_PROD_GAS)
        DO 395 I = 1, IIGAS
          CALL BRKLINE(IGASNAME(I), BROKLINE, 2, 40, BRKSTRNG,
     1                 BRKLEN, 4, LINES)
          IF (IREV_GAS(I) .EQ. 0) THEN
            WRITE(LOUT,8351) I, BROKLINE(1), ND_RCF(I), ND_RCR(I),
     1                       ND_RCF(I)/GAS_DA_FACTOR,
     2                       ND_RCR(I)/GAS_DA_FACTOR
            FLAG = .TRUE.
          ELSE
            WRITE(LOUT,8352) I, BROKLINE(1), ND_RCF(I), ND_RCR(I),
     1                       ND_RCF(I)/GAS_DA_FACTOR,
     2                       ND_RCR(I)/GAS_DA_FACTOR
          END IF
          IF (LINES .GT. 1) WRITE(LOUT,8353) BROKLINE(2)
395     CONTINUE
        WRITE(LOUT,8354)
        IF (FLAG) WRITE(LOUT,8355)
        IF (LTRAN) THEN
          LT  = CKLSCH(KNAME(ID_CARR))
          LT2 = CKLSCH(KNAME(ID_MAJOR))
          WRITE(LOUT,8356) KNAME(ID_CARR)(1:LT), KNAME(ID_MAJOR)(1:LT2),
     1                     CTOT_BATH, DJK(ID_CARR, ID_MAJOR),
     2                     LENGTH_SCALE, GAS_DA_FACTOR
        ELSE
          WRITE(LOUT,8356) 'O2', 'N2', CTOT_BATH,
     1                     DJK(ID_CARR, ID_MAJOR), LENGTH_SCALE,
     2                     GAS_DA_FACTOR
        END IF
      END IF
*
C-----------------------------------------------------------------------
C                         BLOCK 3.5
C     TABLE OF NON-DIMENSIONAL SURFACE REACTION RATE CONSTANTS
C                 AT THE BATH GAS CONDITIONS
C-----------------------------------------------------------------------
*
      IF (LNDIM_SUR .AND. IISUR .GT. 0) THEN
      IF (.NOT. (LNDIM_GAS .AND. IIGAS .GT. 0)) THEN
        IF (LTRAN) THEN
          CALL MCSDIF((P_BATH/760.*PATM), T_BATH, KKGAS, RMCWRK, DJK)
        ELSE
          DJK(ID_CARR, ID_MAJOR)=1.241E-5 * T_BATH**1.707 *(760./P_BATH)
        END IF
      END IF
      SURF_DA_FACTOR  =  CTOT_BATH *DJK(ID_CARR, ID_MAJOR)
     $                 / LENGTH_SCALE
      WRITE(LOUT,8359) P_BATH, T_BATH
      FLAG = .FALSE.
      DO 400 I = 1, IISUR
      CALL NDCALC_SUR(K_STAR, K_STAR_NEG1, I, ISTFL_SUR(I), RA_SUR(I),
     1            RB_SUR(I), RE_SUR(I), RP_SUR(I), WT(K_GAS_REAC(I)),
     2            SDTOT, GAS_REAC(I), SUR_REAC(I), BULK_REAC(I), SDEN,
     3            KFIRST, KLAST, KSTOICF_SUR,
     4            T_BATH, X_BATH, ICOV_SUR(I),
     5            NCOVI, KCOVI, CPARI, NSCOV, ACT, ISKWRK, SKWRK,
     6            SEQKC, GAS_PROD(I), SUR_PROD(I), BULK_PROD(I),
     7            KSTOIC_SUR, KOCC)
        CALL BRKLINE(ISURNAME(I), BROKLINE, 2, 40, BRKSTRNG,
     1                    BRKLEN, 4, LINES)
        IF (NRPP_SUR(I) .LT. 0) THEN
          WRITE(LOUT,8361) I, BROKLINE(1), K_STAR, K_STAR_NEG1,
     1                     K_STAR/SURF_DA_FACTOR,
     2                     K_STAR_NEG1/SURF_DA_FACTOR
          FLAG = .TRUE.
        ELSE
          WRITE(LOUT,8362) I, BROKLINE(1), K_STAR, K_STAR_NEG1,
     1                     K_STAR/SURF_DA_FACTOR,
     2                     K_STAR_NEG1/SURF_DA_FACTOR
        END IF
        IF (LINES .GT. 1) WRITE(LOUT,8363) BROKLINE(2)
400   CONTINUE
      WRITE(LOUT,8364)
      IF (FLAG) WRITE(LOUT,8365)
      IF (LTRAN) THEN
        LT  = CKLSCH(KNAME(ID_CARR))
        LT2 = CKLSCH(KNAME(ID_MAJOR))
        WRITE(LOUT,8366) KNAME(ID_CARR)(1:LT), KNAME(ID_MAJOR)(1:LT2),
     1                   CTOT_BATH, DJK(ID_CARR, ID_MAJOR),
     2                   LENGTH_SCALE, SURF_DA_FACTOR
      ELSE
        WRITE(LOUT,8366) 'O2', 'N2', CTOT_BATH, DJK(ID_CARR, ID_MAJOR),
     1                   LENGTH_SCALE, SURF_DA_FACTOR
      END IF
      END IF
C-----------------------------------------------------------------------
C                         BLOCK 4
C   MOLECULAR SPECIES THERMODYNAMICS and TRANSPORT PROPERTIES
C                         PRINTOUT
C-----------------------------------------------------------------------
C
C Get thermo quantities at a range of temperatures
C
      T_OLD = TTABLE_MIN - TTABLE_DELTAT
      FLAG  = .FALSE.
      FLAG2 = .FALSE.
      ID_T_BATH = 0
      DO 450 ITEMP = 1, TTABLE_NUM
C
      IF( ( (TTABLE_MIN .GT. T298) .OR.
     1    (T_OLD .LT. T298 .AND. (T_OLD+TTABLE_DELTAT) .GT. T298))
     2             .AND. (.NOT. FLAG) ) THEN
        IF (T_BATH .LT. T298 .AND. .NOT. FLAG2) THEN
          T = T_BATH
          FLAG2 = .TRUE.
          ID_T_BATH = ITEMP
        ELSE
          T = T298
          FLAG  = .TRUE.
          IF (T_BATH .EQ. T298) THEN
            FLAG2 = .TRUE.
            ID_T_BATH = ITEMP
          END IF
        END IF
      ELSE IF( ( (TTABLE_MIN .GT. T_BATH) .OR.
     1      (T_OLD .LT. T_BATH .AND. (T_OLD+TTABLE_DELTAT) .GT. T_BATH))
     2              .AND. (.NOT. FLAG2)    ) THEN
        T = T_BATH
        FLAG2 = .TRUE.
        ID_T_BATH = ITEMP
      ELSE
        T     = T_OLD + TTABLE_DELTAT
        T_OLD = T_OLD + TTABLE_DELTAT
        IF (T_BATH .EQ. T) THEN
          FLAG2 = .TRUE.
          ID_T_BATH = ITEMP
        END IF
        IF (T298 .EQ. T) FLAG = .TRUE.
      END IF
C
      IF (T .GT. TTABLE_MAX) THEN
        IF (.NOT. FLAG) THEN
          FLAG = .TRUE.
          T = T298
          IF (T_BATH .EQ. T298) THEN
            FLAG2 = .TRUE.
            ID_T_BATH = ITEMP
          END IF
        ELSE IF (.NOT. FLAG2) THEN
          T = T_BATH
          FLAG2 = .TRUE.
          ID_T_BATH = ITEMP
        ELSE IF (ITEMP .GE. 4) THEN
          TTABLE_NUM = ITEMP - 1
          GO TO 449
        END IF
      END IF
C
      IF (T .EQ. 1000.) THEN
        TTABLE(ITEMP) = 1000.01
      ELSE
        TTABLE(ITEMP) = T
      END IF
C
C Get the nondimensional ss specific heat at constant pressure
      CALL SKCPOR (TTABLE(ITEMP), ISKWRK, SKWRK, CP(1,ITEMP))
C Get the nondimensional ss enthalpy at temperature T
      CALL SKHORT (TTABLE(ITEMP), ISKWRK, SKWRK, H(1,ITEMP))
C Get the nondimensional ss entropy at temperature T
      CALL SKSOR  (TTABLE(ITEMP), ISKWRK, SKWRK, S(1,ITEMP))
C
c Convert to dimensional units and calc G from H and T
c     S and CP is in cals/mole*k, while H and G are in kcals/mole
      DO 440 K = 1, KKTOT
         CP(K,ITEMP) = CP(K,ITEMP)  * RUC
         H(K,ITEMP)  = H(K,ITEMP)   * RUC * TTABLE(ITEMP) / 1000.
         S(K,ITEMP)  = S(K,ITEMP)   * RUC
         G(K,ITEMP)  = H(K,ITEMP) - TTABLE(ITEMP) * S(K,ITEMP) / 1000.
 440  CONTINUE
 449  CONTINUE
 450  CONTINUE
C
C Print out JANAF-like table for each species
C     (loop through species - only print out a table if LTHRM is true)
C
      DO 500 K = 1 , KKTOT
      IF (LTHRM(K)) THEN
C
      WRITE(LOUT,8000)
      IPHASE = PHASEID (K, NNPHAS, KFIRST, KLAST)
      LT  = CKLSCH(KNAME(K))
      LT2 = CKLSCH(PNAME(IPHASE))
      WRITE(LOUT,8005) KNAME(K)(1:LT), PNAME(IPHASE)(1:LT2)
      WRITE(LOUT,8015) K, K-KFIRST(IPHASE)+1, PNAME(IPHASE)(1:LT2)
      IF (KCHARG(K) .NE. 0) WRITE(LOUT,8006) KCHARG(K)
      WRITE(LOUT,8007)
      FLAG = .FALSE.
      DO 480 M = 1, NELEM
        IF (NEL(M,K) .NE. 0) THEN
          LT = MAX(2,CKLSCH(ENAME(M)))
          WRITE(LOUT,8008) ENAME(M)(1:LT), NEL(M,K)
          FLAG = .TRUE.
        END IF
 480  CONTINUE
      IF (.NOT. FLAG) WRITE(LOUT,8011)
      IF (NNSURF .GT. 0) THEN
        IF (K .GE. KFIRST(NFSURF) .AND. K .LE. KLAST(NLSURF))
     1    WRITE(LOUT, 8012) KOCC(K)
      END IF
      IF (NNBULK .GT. 0) THEN
        IF (K .GE. KFIRST(NFBULK)) THEN
       CALL SKDEN(PATM, T298, ACT, SDEN, ISKWRK, SKWRK, DEN)
       IF (DEN(K) .EQ. -1.0) THEN
         WRITE(LOUT,8013)
       ELSE
         WRITE(LOUT,8014) DEN(K)
         WRITE(LOUT,8023) ACT(K)
       END IF
      END IF
      END IF
C
C     Write out the gas-phase table with transport parameters
C
      IF (K .LE. KKGAS .AND. LTRAN) THEN
        WRITE(LOUT,8016) EPS(K), SIG(K), DIP(K), POL(K), ZROT(K)
        IF (NLIN(K) .EQ. 0) THEN
          WRITE(LOUT,8017)
        ELSE IF (NLIN(K) .EQ. 1) THEN
          WRITE(LOUT,8018)
        ELSE IF (NLIN(K) .EQ. 2) THEN
          WRITE(LOUT,8019)
        END IF
        LT = CKLSCH(KNAME(ID_CARR))
        STRING = ' '
        STRING = 'Dif_Co_with_'//KNAME(ID_CARR)(1:LT)
        CALL CENT_JST(STRING, 26)
        WRITE(LOUT,8020) H298(K), WT(K), STRING(1:26)
C
C       Loop through temperatures in the table
C
        DO 490 ITEMP = 1 , TTABLE_NUM
          CALL MCSVIS (TTABLE(ITEMP), RMCWRK, VIS)
          CALL MCSCON (TTABLE(ITEMP), RMCWRK, CON)
          CALL MCSDIF ((P_BATH/760.*PATM), TTABLE(ITEMP), KKGAS,
     1                 RMCWRK, DJK)
          WRITE(LOUT,8021) TTABLE(ITEMP), H(K,ITEMP)-H298(K),
     1                     G(K,ITEMP) - H298(K), CP(K,ITEMP),
     2                     S(K,ITEMP), VIS(K), CON(K), DJK(ID_CARR,K)
 490    CONTINUE
        WRITE(LOUT,8022) P_BATH
      ELSE
C
C     Write out the surface phase table without transport parameters
C
        WRITE(LOUT,8009) H298(K), WT(K)
C       Loop through temperatures in the table
        DO 495 ITEMP = 1 , TTABLE_NUM
          WRITE(LOUT,8010) TTABLE(ITEMP), H(K,ITEMP)-H298(K),
     1                     G(K,ITEMP) - H298(K), CP(K,ITEMP), S(K,ITEMP)
 495    CONTINUE
        WRITE(LOUT,8000)
      END IF
      WRITE(LOUT,*)
C
      END IF
 500  CONTINUE
*
C-----------------------------------------------------------------------
C                         BLOCK 5.0
C       CALCULATION OF THERMO QUANTITIES FOR GAS-PHASE REACTIONS
C (This section is done irrespective of whether a printout is requested)
C-----------------------------------------------------------------------
*
* Only do this section if there are gas-phase reactions
*
      IF (IIGAS .GT. 0) THEN
        DO 600 ITEMP = 1 , TTABLE_NUM
        DO 600 I = 1, IIGAS
          DELTA_S_GAS(I,ITEMP) = 0.0
          DELTA_H_GAS(I,ITEMP) = 0.0
          DO 550 K = 1, KKGAS
            DELTA_S_GAS(I,ITEMP) = DELTA_S_GAS(I,ITEMP) +
     $                             KSTOIC_GAS(K,I)*S(K,ITEMP)
            DELTA_H_GAS(I,ITEMP) = DELTA_H_GAS(I,ITEMP) +
     $                             KSTOIC_GAS(K,I)*H(K,ITEMP)
 550      CONTINUE
          DELTA_G_GAS(I,ITEMP) = DELTA_H_GAS(I,ITEMP) -
     $                     TTABLE(ITEMP) * DELTA_S_GAS(I,ITEMP) / 1000.
 600    CONTINUE
      END IF
*
C-----------------------------------------------------------------------
C                         BLOCK 5.1
C            SUMMARY TABLE OF THERMO FOR GAS-PHASE REACTIONS
C               Print out a table at the bath gas conditions
C               describing the overall thermo functions
C-----------------------------------------------------------------------
*
      IF (LTSUM_GAS .AND. (ID_T_BATH .GT. 0) .AND. (IIGAS .GT. 0)) THEN
        WRITE(LOUT,8210) T_BATH
        DO 610 I = 1, IIGAS
            CALL BRKLINE(IGASNAME(I), BROKLINE, 2, 40, BRKSTRNG,
     1                   BRKLEN, 4, LINES)
            WRITE(LOUT,8211) I, BROKLINE(1), DELTA_G_GAS(I,ID_T_BATH),
     1                DELTA_H_GAS(I,ID_T_BATH), DELTA_S_GAS(I,ID_T_BATH)
            IF (LINES .GT. 1) WRITE(LOUT,8212) BROKLINE(2)
610     CONTINUE
        WRITE(LOUT,8213)
      END IF
*
C-----------------------------------------------------------------------
C                         BLOCK 5.2
C            LONG DESCRIPTION OF GAS-PHASE REACTIONS
C               Print out a table for each individual reaction
C-----------------------------------------------------------------------
*
* Loop through the reactions, and print out a detailed description
* if LGRXN(I) is turned on from the input file
*
      DO 1000 I = 1, IIGAS
      IF (LGRXN(I)) THEN
*
* Print Out a detailed header, containing the symbol description
*
      WRITE(LOUT,5000)
      WRITE(LOUT,8381) I, IGASNAME(I)
      WRITE(LOUT,8400) GAS_MOLE_CHANGE_GAS(I)
      IF (RP_GAS(I) .NE. 1.0) WRITE(LOUT,8410) RP_GAS(I)
*
* Output Third body information, if any
*
      IF (ITHB(I) .EQ. -1) THEN
        WRITE(LOUT,8385)
      ELSE IF (ITHB(I) .EQ. 0) THEN
        WRITE(LOUT,8386)
      ELSE IF (ITHB(I) .GT. 0) THEN
        WRITE(LOUT,8387) ITHB(I)
        DO 620 K = 1 , KKGAS
        IF (AKI(K,I) .NE. 1.0) THEN
           LT = CKLSCH(KNAME(K))
           WRITE(LOUT,8388) KNAME(K)(1:LT), AKI(K,I)
        END IF
620     CONTINUE
      END IF
*
* Calculate the units for the forward rate constant
*
      F_RATE_UNITS = -1
      DO 640 K = 1, KKGAS
        IF (KSTOICF_GAS(K,I) .LT. 0) THEN
          F_RATE_UNITS = F_RATE_UNITS - KSTOICF_GAS(K,I)
        END IF
640   CONTINUE
      IF (ITHB(I) .GE. 0 .AND. IFOP(I) .EQ. 0) THEN
         F_RATE_UNITS = F_RATE_UNITS + 1
         STRING = ' (including the third body)'
         LT = 27
      ELSE
         STRING = ' '
         LT = 1
      END IF
*
* Calculate the units for the reverse rate constant
*
      R_RATE_UNITS = F_RATE_UNITS + GAS_MOLE_CHANGE_GAS(I)
*
* Handle the output of the reversibility vector
*
      IF (IREV_GAS(I) .GT. 0) THEN
        WRITE(6,8510) F_RATE_UNITS+1, R_RATE_UNITS+1, STRING(1:LT)
      ELSE IF(IREV_GAS(I) .LT. 0) THEN
        WRITE(6,8520) F_RATE_UNITS+1, R_RATE_UNITS+1, STRING(1:LT)
      ELSE
         PRINT *,'ERROR IREV_GAS(I) = 0'
C         STOP
         RETURN
      END IF
*
* Print out the high pressure arrhenius parameters:
*
      CALL P_REAC(LOUT, RATE_UNITS(F_RATE_UNITS), RA_GAS(I),
     1            RB_GAS(I), RE_GAS(I), RP_GAS(I), 1)
*
* Print out any fall-off information
*
      IF (IFOP(I) .NE. 0) THEN
        IF (IFOP(I) .EQ. 1) THEN
          WRITE(LOUT,8561) RATE_UNITS(F_RATE_UNITS+1),
     1      FPAR(1,I), FPAR(2,I), FPAR(3,I)*RKCAL
        ELSE IF (IFOP(I) .EQ. 2) THEN
          WRITE(LOUT,8562) RATE_UNITS(F_RATE_UNITS+1),
     1      FPAR(1,I), FPAR(2,I), FPAR(3,I)*RKCAL,
     2      FPAR(4,I), FPAR(5,I), FPAR(6,I), FPAR(7,I), FPAR(8,I)
        ELSE IF (IFOP(I) .EQ. 3) THEN
          WRITE(LOUT,8563) RATE_UNITS(F_RATE_UNITS+1),
     1      FPAR(1,I), FPAR(2,I), FPAR(3,I)*RKCAL,
     2      FPAR(4,I), FPAR(5,I), FPAR(6,I)
        ELSE IF (IFOP(I) .EQ. 4) THEN
          WRITE(LOUT,8564) RATE_UNITS(F_RATE_UNITS+1),
     1      FPAR(1,I), FPAR(2,I), FPAR(3,I)*RKCAL,
     2      FPAR(4,I), FPAR(5,I), FPAR(6,I), FPAR(7,I)
        END IF
        IF (KFAL(I) .NE. 0) THEN
          LT = CKLSCH(KNAME(KFAL(I)))
          WRITE(LOUT,8570) KNAME(KFAL(I))(1:LT)
        END IF
      END IF
*
* Print out headings for the high pressure rate constant table
*
      IF (IREV_GAS(I) .GT. 0) THEN
        WRITE(LOUT,8601)
      ELSE
        WRITE(LOUT,8600)
      END IF
      WRITE(LOUT,8610) RATE_UNITS(F_RATE_UNITS),
     $                 RATE_UNITS(R_RATE_UNITS)
*
* Loop through the body of the high pressure rate constant table
* calculating reaction information for each temperature in the table
*
      DO 900 ITEMP = 1 , TTABLE_NUM
        T = TTABLE(ITEMP)
*
*        Calculate the high pressure rate constant at two different
*        temperatures to get the arrhenius dependence
*
        CALL HIGH_PRESS(RF_DELTA, A_FACTOR, EACT, RA_GAS(I),
     $                  RB_GAS(I), RE_GAS(I), RP_GAS(I),
     $                  T + DELTAT)
        CALL HIGH_PRESS(RF, A_FACTOR, EACT, RA_GAS(I), RB_GAS(I),
     $                  RE_GAS(I), RP_GAS(I), T)
*
*        Calculate the equilibrium constant- doesn't depend on p or
*        X_bath
*
        CALL CKEQXP (PATM, T, X_BATH, ICKWRK, CKWRK, GEQKC)
        CALL CKEQXP (PATM, T+DELTAT, X_BATH, ICKWRK, CKWRK,
     $               GEQKC_DELTA)
*
*        Calculate the high pressure reverse rate constant:
*
        RF_NEG1       = RF       / MAX(GEQKC(I), SMALL)
        RF_NEG1_DELTA = RF_DELTA / MAX(GEQKC_DELTA(I), SMALL)
        CALL CALC_ARR(RF_NEG1, T, RF_NEG1_DELTA, (T+DELTAT),
     $                A_FACTOR_NEG1, EACT_NEG1)
*
*        Calculate Non-dimensional gas rate constants
*
        CALL NDCALC_GAS (ND_RCF, ND_RCR, GAS_REAC_GAS,
     $                   P_CGS, T, X_BATH, ICKWRK, CKWRK,
     $                   GAS_PROD_GAS)
*
*        Write all of the above out on a single line
*
        WRITE(6,8700) T, RF, A_FACTOR, EACT, DELTA_G_GAS(I,ITEMP),
     $                DELTA_H_GAS(I,ITEMP), DELTA_S_GAS(I,ITEMP),
     $                RF_NEG1, A_FACTOR_NEG1, EACT_NEG1, ND_RCF(I),
     $                ND_RCR(I)
900   CONTINUE
      WRITE(LOUT,5000)
*
C-----------------------------------------------------------------------
C                         BLOCK 5.3
C            FALL-OFF CURVES AS A FUNCTION OF TEMPERATURE
C            Print out a table of fall-off information
C-----------------------------------------------------------------------
*
*        Only print out table if rxn has fall-off effects, and
*        the fall-off logical variable, LTFAL, is turned on from
*        the input file
*
      IF (LTFAL(I) .AND. IFOP(I) .NE. 0) THEN
        WRITE(LOUT,8705) I, P_BATH
        WRITE(LOUT,8710) RATE_UNITS(F_RATE_UNITS),
     $                   RATE_UNITS(F_RATE_UNITS+1),
     $                   RATE_UNITS(R_RATE_UNITS)
        DO 920 ITEMP = 1, TTABLE_NUM
          T = TTABLE(ITEMP)
*
*            Calculate the reaction rate constant for rxn I at the
*            current conditions -> two temperatures to get
*            local arrhenius parameters
*
          CALL CKRCXP (P_CGS, T + DELTAT, X_BATH, ICKWRK, CKWRK,
     $                 ND_RCF, ND_RCR)
          RF_DELTA      = ND_RCF(I)
          RF_NEG1_DELTA = ND_RCR(I)
          CALL CKRCXP (P_CGS, T, X_BATH, ICKWRK, CKWRK, ND_RCF, ND_RCR)
          RF            = ND_RCF(I)
          RF_NEG1       = ND_RCR(I)
*
*  Calculate the high pressure reaction rate constant
*
          CALL HIGH_PRESS(RCFINF, A_FACTOR, EACT, RA_GAS(I), RB_GAS(I),
     $                    RE_GAS(I), RP_GAS(I), T)
*
*  Calculate the arrhenius parameters for the forward and reverse
*  directions of the rate constants
*
          CALL CALC_ARR (RF, T, RF_DELTA, (T+DELTAT),
     $                   A_FACTOR, EACT)
          CALL CALC_ARR (RF_NEG1, T, RF_NEG1_DELTA, (T+DELTAT),
     $                   A_FACTOR_NEG1, EACT_NEG1)
*
*  Obtain a detailed breakdown of the reaction rate constant
*  calculation for the current fall-off reaction
*
          CALL CKFALP (P_CGS, T, X_BATH, ICKWRK, CKWRK, I,
     $                 RKLOW, CTB, PR, FC, PCOR)
*
*  Print everything out on a single line
*
          WRITE(6,8720) T, RF, A_FACTOR, EACT, PCOR, RKLOW,
     $                  PR, FC, CTB, RF_NEG1, A_FACTOR_NEG1, EACT_NEG1
920     CONTINUE
        WRITE(LOUT,5000)
      END IF
*
C-----------------------------------------------------------------------
C                         BLOCK 5.4
C            FALL-OFF CURVES AS A FUNCTION OF PRESSURE
C            Print out a table of fall-off information
C-----------------------------------------------------------------------
*
      IF (LPFAL(I) .AND. IFOP(I) .NE. 0) THEN
*
* Loop over a range of pressures keeping the temperature fixed at the
* bath temperature
*
        T = T_BATH
*
*          Find the low pressure rate constant for printout
*
        CALL CKFALP (P_CGS, T, X_BATH, ICKWRK, CKWRK, I,
     $               RKLOW, CTB, PR, FC, PCOR)
        WRITE(LOUT,8725) I, T_BATH, RKLOW, RATE_UNITS(F_RATE_UNITS+1)
        WRITE(LOUT,8730) RATE_UNITS(F_RATE_UNITS),
     $                   RATE_UNITS(F_RATE_UNITS),
     $                   RATE_UNITS(R_RATE_UNITS)
        DO 940 ITEMP = 1 , PTABLE_NUM
          PTORR = 10.0**(
     $      LOG10(PTABLE_MAX)*(PTABLE_NUM-ITEMP)/(PTABLE_NUM-1.0) +
     $      LOG10(PTABLE_MIN)*(ITEMP-1.0)       /(PTABLE_NUM-1.0) )
          P = PTORR *PATM / 760.
*
*            Calculate the reaction rate constant for rxn I at the
*            current conditions -> two temperatures to get
*            local arrhenius parameters
*
          CALL CKRCXP (P, T + DELTAT, X_BATH, ICKWRK, CKWRK,
     $                 ND_RCF, ND_RCR)
          RF_DELTA      = ND_RCF(I)
          RF_NEG1_DELTA = ND_RCR(I)
          CALL CKRCXP (P, T, X_BATH, ICKWRK, CKWRK, ND_RCF, ND_RCR)
          RF            = ND_RCF(I)
          RF_NEG1       = ND_RCR(I)
*
*  Calculate the high pressure reaction rate constant
*
          CALL HIGH_PRESS(RCFINF, A_FACTOR, EACT, RA_GAS(I), RB_GAS(I),
     $                    RE_GAS(I), RP_GAS(I), T)
*
*  Calculate the arrhenius parameters for the forward and reverse
*  directions of the rate constants
*
          CALL CALC_ARR (RF, T, RF_DELTA, (T+DELTAT),
     $                   A_FACTOR, EACT)
          CALL CALC_ARR (RF_NEG1, T, RF_NEG1_DELTA, (T+DELTAT),
     $                   A_FACTOR_NEG1, EACT_NEG1)
*
*  Obtain a detailed breakdown of the reaction rate constant
*  calculation for the current fall-off reaction
*
          CALL CKFALP (P, T, X_BATH, ICKWRK, CKWRK, I,
     $                 RKLOW, CTB, PR, FC, PCOR)
*
*
        WRITE(6,8740) PTORR, RF, A_FACTOR, EACT, PCOR, RKLOW*CTB,
     1                PR, FC, CTB, RF_NEG1, A_FACTOR_NEG1, EACT_NEG1
940     CONTINUE
        WRITE(LOUT,5000)
      END IF
*
C-----------------------------------------------------------------------
C                         BLOCK 5.5
C         THIRD BODY REACTIONS - LUMPING "M" WITH REACTION RATE
C            Print out a table of information
C-----------------------------------------------------------------------
*
* Print an extra table of information for a rxn IF
*    Third body option is turned on
*    The rxn involves a third body
*    The rxn is not already handled by fall-off section
* The extra information will involve lumping in the third body term
* into the rate constant. This effective rate constant will then be
* reported. The bath gas conditions will be used to calculate the
* "effective concentration" of the bath gas. Effective, here, means
* the concentration of the bath gas after the effects of any
* enhanced third body efficiency factors have been added in.
*
      IF (LGTHB(I) .AND. ITHB(I) .GE. 0 .AND. IFOP(I) .EQ. 0) THEN
*
*         Write out the header
*         - Note that the units of the reaction rate constant have
*           changed by cm**3/mole factor since the bath gas
*           concentration is now lumped in.
*
        WRITE(LOUT,8745) I, P_BATH
        WRITE(LOUT,8750) RATE_UNITS(F_RATE_UNITS-1),
     $                   RATE_UNITS(R_RATE_UNITS-1)
*
*         Create a table wrt temperature
*
        DO 960 ITEMP = 1 , TTABLE_NUM
          T = TTABLE(ITEMP)
*
*           Calculate the reaction rate constant for rxn I at the
*           current conditions -> two temperatures to get
*           local arrhenius parameters
*
          CALL CKRCXP (P_CGS, T + DELTAT, X_BATH, ICKWRK, CKWRK,
     $                 ND_RCF, ND_RCR)
          RF_DELTA      = ND_RCF(I)
          RF_NEG1_DELTA = ND_RCR(I)
*
          CALL CKRCXP (P_CGS, T, X_BATH, ICKWRK, CKWRK, ND_RCF, ND_RCR)
          RF            = ND_RCF(I)
          RF_NEG1       = ND_RCR(I)
*
          CALL CALC_ARR(RF, T, RF_DELTA, (T+DELTAT),
     $                A_FACTOR, EACT)
          CALL CALC_ARR(RF_NEG1, T, RF_NEG1_DELTA, (T+DELTAT),
     $                A_FACTOR_NEG1, EACT_NEG1)
*
*           We call CKFALP in order to get the effective
*           concentration, CTB
*
          CALL CKFALP (P_CGS, T, X_BATH, ICKWRK, CKWRK, I,
     $                 RKLOW, CTB, PR, FC, PCOR)
*
*           Print out the results on a single line
*
          WRITE(6,8760) T, RF, A_FACTOR, EACT, P_BATH/(760.*82.05*T),
     $                  CTB, RF_NEG1, A_FACTOR_NEG1, EACT_NEG1
960     CONTINUE
        WRITE(LOUT,5000)
      END IF
*
*        End of logic block LGRXN(I)
*
      ENDIF
*
*        End of loop over gas phase reactions
*
1000  CONTINUE
*
C-----------------------------------------------------------------------
C                         BLOCK 6.0
C            CALCULATE SURFACE REACTION THERMO PARAMETERS
C-----------------------------------------------------------------------
*
      IF (IISUR .GT. 0) THEN
        DO 1600 ITEMP = 1 , TTABLE_NUM
        DO 1600 I = 1, IISUR
          T = TTABLE(ITEMP)
          DELTA_S_SUR(I,ITEMP) = 0.0
          DELTA_H_SUR(I,ITEMP) = 0.0
          DO 1550 K = 1, KKTOT
            DELTA_S_SUR(I,ITEMP) = DELTA_S_SUR(I,ITEMP) +
     1                             KSTOIC_SUR(I,K)*S(K,ITEMP)
            DELTA_H_SUR(I,ITEMP) = DELTA_H_SUR(I,ITEMP) +
     1                             KSTOIC_SUR(I,K)*H(K,ITEMP)
1550      CONTINUE
          DELTA_G_SUR(I,ITEMP) = DELTA_H_SUR(I,ITEMP) -
     1                           T * DELTA_S_SUR(I,ITEMP) / 1000.
1600    CONTINUE
      END IF
C
C-----------------------------------------------------------------------
C                         BLOCK 6.1
C          WRITE A SUMMARY TABLE OF SURFACE-PHASE REACTIONS
C-----------------------------------------------------------------------
C
      IF (LTSUM_SUR .AND. (ID_T_BATH .GT. 0) .AND. (IISUR .GT. 0)) THEN
        WRITE(LOUT,8209) T_BATH
        DO 1610 I = 1, IISUR
            CALL BRKLINE(ISURNAME(I), BROKLINE, 2, 40, BRKSTRNG,
     1                   BRKLEN, 4, LINES)
            WRITE(LOUT,8211) I, BROKLINE(1), DELTA_G_SUR(I,ID_T_BATH),
     1                DELTA_H_SUR(I,ID_T_BATH), DELTA_S_SUR(I,ID_T_BATH)
            IF (LINES .GT. 1) WRITE(LOUT,8212) BROKLINE(2)
1610    CONTINUE
        WRITE(LOUT,8213)
      END IF
C
C-----------------------------------------------------------------------
C                         BLOCK 6.2
C          WRITE A TABLE FOR EACH SURFACE-PHASE REACTION
C-----------------------------------------------------------------------
C
      IF (IISUR .GT. 0) THEN
      DO 2000 I = 1, IISUR
      IF (LSRXN(I)) THEN
C
      WRITE(LOUT,5000)
      WRITE(LOUT,9301) I, ISURNAME(I)
      WRITE(LOUT,9400) GAS_MOLE_CHANGE_SUR(I), SUR_MOLE_CHANGE_SUR(I),
     1                 BULK_MOLE_CHANGE_SUR(I)
      IF (RP_SUR(I) .NE. 1.0) WRITE(LOUT,8410) RP_SUR(I)
C
C Calculate the units for the forward rate constant.  Also determine
C   the identity of the gas-phase species for later use in calculation
C   of the effusive velocity.
C
      F_RATE_UNITS = GAS_REAC(I)
      F_RATE_SUR_UNITS = SUR_REAC(I) - 1.0
C
C Print out coverage dependent parameters, and calculate units
C     for the reverse direction
C
      F_DBLE = F_RATE_UNITS
      CALL SRX_UNIT(SUR_FUNITS, F_DBLE, F_RATE_SUR_UNITS,
     1              IERR)
      IF (IERR) THEN
        PRINT *,'Overwrite forward rate constant units string'
      ENDIF
      R_RATE_UNITS = F_RATE_UNITS + GAS_MOLE_CHANGE_SUR(I)
      r_rate_sur = F_RATE_SUR_UNITS + SUR_MOLE_CHANGE_SUR(I)
      R_DBLE = R_RATE_UNITS
      CALL SRX_UNIT(SUR_RUNITS, R_DBLE, r_rate_sur,
     1              IERR)
      IF (IERR) THEN
        PRINT *,'Overwrite reverse units string'
      ENDIF
      IF (ICOV_SUR(I) .NE. 0) THEN
        CALL SKICOV(I, NSCOV, ISKWRK, SKWRK, NCOVI, KCOVI, CPARI)
        WRITE(LOUT,9420) NCOVI
        DO 1760 N = 1, NCOVI
          K = KCOVI(N)
          LT = CKLSCH(KNAME(K))
          WRITE(LOUT,9421) KNAME(K)(1:LT),   CPARI(1,N),
     1                     CPARI(2,N),       CPARI(3,N) * RUC
1760    CONTINUE
      END IF
C
C Handle the output of the reversibility vector
C
      IF (NRPP_SUR(I) .GT. 0) THEN
        WRITE(LOUT,9510) GAS_REAC(I), SUR_REAC(I), BULK_REAC(I),
     1                   GAS_PROD(I), SUR_PROD(I), BULK_PROD(I)
      ELSE
        WRITE(LOUT,9520) GAS_REAC(I), SUR_REAC(I), BULK_REAC(I),
     1                   GAS_PROD(I), SUR_PROD(I), BULK_PROD(I)
      END IF
C
C Handle the output of whether this is a sticking coefficient reaction
C   or whether it could have been considered a sticking coefficient
C   reaction. - FORWARD DIRECTION
C
      IF (ISTFL_SUR(I) .GT. 0) THEN
        WRITE(LOUT,9528)
        IF (IMOTZ_SUR(I) .EQ. 1) THEN
          WRITE(LOUT,9529)
        ELSE
          WRITE(LOUT,9530)
        ENDIF
        CALL P_REAC(LOUT, SUR_FUNITS, RA_SUR(I),
     1              RB_SUR(I), RE_SUR(I), RP_SUR(I), 5)
C       Save sticking coefficient parameters -
        RA_STK = RA_SUR(I)
        RB_STK = RB_SUR(I)
        RE_STK = RE_SUR(I)
C
C       Convert sticking coefficient expression into general rate
C             constant expression
C
        CALL STKTOGEN(RA_STK, RB_STK, RE_STK, WT(K_GAS_REAC(I)),
     1                SUR_REAC(I), SDTOT, TTABLE, RA_GEN,
     2                RB_GEN, RE_GEN, WORK, IERR)
        IF (IERR) WRITE(LOUT, 9539)
        WRITE(LOUT,9536)
        CALL P_REAC(LOUT, SUR_FUNITS, RA_GEN, RB_GEN, RE_GEN,
     1              RP_SUR(I), 1)
      ELSE IF (GAS_REAC(I) .EQ. 1 .AND. GAS_PROD(I) .EQ. 1) THEN
        WRITE(LOUT,9531)
        CALL P_REAC(LOUT, SUR_FUNITS, RA_SUR(I),
     1              RB_SUR(I), RE_SUR(I), RP_SUR(I), 1)
        CALL GENTOSTK(RA_STK, RB_STK, RE_STK, WT(K_GAS_REAC(I)),
     1                SUR_REAC(I), SDTOT, TTABLE,
     1                RA_SUR(I), RB_SUR(I), RE_SUR(I), WORK, IERR)
        IF (IERR) THEN
          WRITE(LOUT,9534)
        ELSE
          WRITE(LOUT,9535)
        END IF
        CALL P_REAC(LOUT, SUR_FUNITS, RA_STK, RB_STK, RE_STK,
     1              RP_SUR(I), 3)
        IF (IERR) WRITE(LOUT, 9539)
      ELSE IF (GAS_REAC(I) .EQ. 1) THEN
        WRITE(LOUT,9532)
        CALL P_REAC(LOUT, SUR_FUNITS, RA_SUR(I),
     1              RB_SUR(I), RE_SUR(I), RP_SUR(I), 1)
        CALL GENTOSTK(RA_STK, RB_STK, RE_STK, WT(K_GAS_REAC(I)),
     1                SUR_REAC(I), SDTOT, TTABLE,
     2                RA_SUR(I), RB_SUR(I), RE_SUR(I), WORK, IERR)
        IF (IERR) THEN
          WRITE(LOUT,9534)
        ELSE
          WRITE(LOUT,9535)
        END IF
        CALL P_REAC(LOUT, SUR_FUNITS, RA_STK, RB_STK, RE_STK,
     1              RP_SUR(I), 3)
        IF (IERR) WRITE(LOUT, 9539)
      ELSE
        CALL P_REAC(LOUT, SUR_FUNITS, RA_SUR(I),
     1              RB_SUR(I), RE_SUR(I), RP_SUR(I), 1)
      END IF
C
C Calculate the reverse rate constant arrhenius parameters
C   Handle the output of whether this is a sticking coefficient reaction
C   or whether it could have been considered a sticking coefficient
C   reaction. - REVERSE DIRECTION
C
      CALL REVGEN(SKWRK, ISKWRK, I, RA_REV, RB_REV, RE_REV,
     1            WT(K_GAS_REAC(I)), SUR_REAC(I), SDTOT, TTABLE,
     2            ISTFL_SUR(I), RA_SUR(I), RB_SUR(I), RE_SUR(I),
     3            ACT, SDEN, SEQKC, WORK)
      IF (NRPP_SUR(I) .LE. 0) THEN
        WRITE(LOUT,9541)
      ELSE
        WRITE(LOUT,9537)
      END IF
      CALL P_REAC(LOUT, SUR_RUNITS, RA_REV, RB_REV, RE_REV,
     1            RP_SUR(I), 2)
      IF (GAS_PROD(I) .EQ. 1) THEN
        CALL REVSTK(SKWRK, ISKWRK, I, RA_REV_ST, RB_REV_ST, RE_REV_ST,
     1              WT(K_GAS_REAC(I)), SUR_REAC(I), WT(K_GAS_PROD(I)),
     2              SUR_PROD(I), SDTOT, TTABLE, ISTFL_SUR(I),
     3              RA_SUR(I), RB_SUR(I), RE_SUR(I),
     4              ACT, SDEN, SEQKC, WORK, IERR)
        WRITE(LOUT,9538)
        CALL P_REAC(LOUT, SUR_RUNITS, RA_REV_ST, RB_REV_ST,
     1              RE_REV_ST, RP_SUR(I), 4)
        IF (IERR) WRITE(LOUT, 9540)
      END IF
C
C Handle the header output for the case where reverse
C     arrhenius parameters are explicitly supplied:
C
      IF (IREV_SUR(I) .GT. 0) THEN
        CALL SKIREV(I, ISKWRK, SKWRK, IREV_SUR(I), RAR, RBR, RER)
        CALL P_REAC(LOUT, SUR_RUNITS, RAR, RBR, RER, RP_SUR(I), 2)
      END IF
C
C Print out headings for the high pressure rate constant table
C
      IF (NRPP_SUR(I) .LE. 0) THEN
        WRITE(LOUT,9599)
      ELSE
        WRITE(LOUT,9600)
      END IF
      WRITE(LOUT,9601) SUR_FUNITS, SUR_RUNITS
C
C Loop through the body of the high pressure rate constant table
C
      DO 1900 ITEMP = 1 , TTABLE_NUM
      T = TTABLE(ITEMP)
C     Calculate the high pressure rate constant:
      IF (ISTFL_SUR(I) .EQ. 0) THEN
        CALL HIGH_PRESS(RF_DELTA , A_FACTOR , EACT, RA_SUR(I),
     1                  RB_SUR(I) , RE_SUR(I), RP_SUR(I), T + DELTAT)
        CALL HIGH_PRESS(RF, A_FACTOR, EACT, RA_SUR(I), RB_SUR(I),
     1                  RE_SUR(I), RP_SUR(I), T)
      ELSE
        CALL STKCALCA(RA_STK, RB_STK, RE_STK, RP_SUR(I),
     1                WT(K_GAS_REAC(I)), SUR_REAC(I), SDTOT, T + DELTAT,
     2                RF_DELTA, A_FACTOR, EACT)
        CALL STKCALCA(RA_STK, RB_STK, RE_STK, RP_SUR(I),
     1                WT(K_GAS_REAC(I)), SUR_REAC(I), SDTOT, T, RF,
     2                A_FACTOR, EACT)
      END IF
C
C Calculate the equilibrium constant- doesn't depend on p or y
C  (note that this subroutine handles the irev_sur=1 case correctly,
C   ie. - it will return the equilibrium constant consistent with
C         the reaction rate input)
C
      CALL SKEQ (PATM, T, ACT, SDEN, ISKWRK, SKWRK, SEQKC)
      CALL SKEQ (PATM, T+DELTAT, ACT, SDEN, ISKWRK, SKWRK, SEQKC_DELTA)
C
C     Calculate the high pressure reverse rate constant:
C
      RF_NEG1       = RF       / MAX(SEQKC(I)      ,SMALL)
      RF_NEG1_DELTA = RF_DELTA / MAX(SEQKC_DELTA(I),SMALL)
      CALL CALC_ARR(RF_NEG1, T, RF_NEG1_DELTA, (T+DELTAT),
     1              A_FACTOR_NEG1, EACT_NEG1)
C
C Calculate k_star, the non-dimensional rate_constant
C
      CALL NDCALC_SUR(K_STAR, K_STAR_NEG1, I, ISTFL_SUR(I), RA_SUR(I),
     1            RB_SUR(I), RE_SUR(I), RP_SUR(I), WT(K_GAS_REAC(I)),
     2            SDTOT, GAS_REAC(I), SUR_REAC(I), BULK_REAC(I), SDEN,
     3            KFIRST, KLAST, KSTOICF_SUR,
     4            T, X_BATH, ICOV_SUR(I), NCOVI,
     5            KCOVI, CPARI, NSCOV, ACT, ISKWRK, SKWRK, SEQKC,
     6            GAS_PROD(I), SUR_PROD(I), BULK_PROD(I), KSTOIC_SUR,
     7            KOCC)
      WRITE(LOUT,9602) T, RF, A_FACTOR, EACT, DELTA_G_SUR(I,ITEMP),
     1         DELTA_H_SUR(I,ITEMP), DELTA_S_SUR(I,ITEMP),
     2         RF_NEG1, A_FACTOR_NEG1, EACT_NEG1, K_STAR, K_STAR_NEG1
1900  CONTINUE
C
C Special Section to analyse surface coverage dependent reaction rates
C
      IF (ICOV_SUR(I) .NE. 0 .AND. LSCOV(I)) THEN
C
C     Print out headings for the high pressure rate constant table
C
      WRITE(LOUT,9610)
      WRITE(LOUT,9611) SUR_FUNITS, SUR_FUNITS, SUR_RUNITS, SUR_RUNITS
C
C     Loop through the body of the high pressure rate constant table
C
      DO 1920 ITEMP = 1 , TTABLE_NUM
      T = TTABLE(ITEMP)
      IF (ISTFL_SUR(I) .EQ. 0) THEN
        CALL HIGH_PRESS(RF_DELTA , A_FACTOR , EACT, RA_SUR(I),
     1                  RB_SUR(I) , RE_SUR(I), RP_SUR(I), T + DELTAT)
        CALL HIGH_PRESS(RF, A_FACTOR, EACT, RA_SUR(I), RB_SUR(I),
     1                  RE_SUR(I), RP_SUR(I), T)
      ELSE
        CALL STKCALCA(RA_STK, RB_STK, RE_STK, RP_SUR(I),
     1                WT(K_GAS_REAC(I)), SUR_REAC(I), SDTOT, T + DELTAT,
     2                RF_DELTA, A_FACTOR, EACT)
        CALL STKCALCA(RA_STK, RB_STK, RE_STK, RP_SUR(I),
     1                WT(K_GAS_REAC(I)), SUR_REAC(I), SDTOT, T, RF,
     2                A_FACTOR, EACT)
      END IF
      CALL COV_DEPA(X_BATH, T + DELTAT, NCOVI, KCOVI, CPARI, NSCOV,
     1              FAC_DELTA, A_FACTOR2, EACT2)
      CALL COV_DEPA(X_BATH, T, NCOVI, KCOVI, CPARI, NSCOV,
     1              FAC, A_FACTOR2, EACT2)
      A_FACTOR = A_FACTOR * A_FACTOR2
      EACT = EACT + EACT2
      RF_PRIME = RF * FAC
      RF_DELTA = RF_DELTA * FAC_DELTA
C
C Calculate the equilibrium constant- doesn't depend on p or y
C  (note that this subroutine handles the irev_sur=1 case correctly,
C   ie. - it will return the equilibrium constant consistent with
C         the reaction rate input)
C
      CALL SKEQ (PATM, T, ACT, SDEN, ISKWRK, SKWRK, SEQKC)
      CALL SKEQ (PATM, T+DELTAT, ACT, SDEN, ISKWRK, SKWRK, SEQKC_DELTA)
C
C Calculate the high pressure reverse rate constant:
C
      RF_NEG1       = RF       / MAX(SEQKC(I)      ,SMALL)
      RF_PRIME_NEG1 = RF_PRIME / MAX(SEQKC(I)      ,SMALL)
      RF_NEG1_DELTA = RF_DELTA / MAX(SEQKC_DELTA(I),SMALL)
      CALL CALC_ARR(RF_PRIME_NEG1, T, RF_NEG1_DELTA, (T+DELTAT),
     1              A_FACTOR_NEG1, EACT_NEG1)
C
      WRITE(LOUT,9612) T, RF_PRIME, A_FACTOR, EACT, RF, FAC, RF_NEG1,
     1                 RF_PRIME_NEG1, A_FACTOR_NEG1, EACT_NEG1
C
1920  CONTINUE
      END IF
C
C Describe surface reaction in terms of sticking coefficient,
C     if possible.
C
      IF (LSTCK(I)) THEN
      IF (ISTFL_SUR(I) .NE. 0 .OR. GAS_REAC(I) .EQ. 1) THEN
        WRITE(LOUT,9700) SDTOT**(SUR_REAC(I)), SUR_REAC(I),
     1                   2*SUR_REAC(I)
        IF (RP_SUR(I) .NE. 1.0) WRITE(LOUT,9706)RP_SUR(I)
        WRITE(LOUT,9701) SUR_FUNITS
        DO 1950 ITEMP = 1, TTABLE_NUM
        T = TTABLE(ITEMP)
        IF (ISTFL_SUR(I) .GT. 0) THEN
          CALL STKACT(STK, A_FACTOR , EACT, RA_STK,
     1                  RB_STK , RE_STK, T, FLAG)
        ELSE
          CALL CALCSTKA(STK, WT(K_GAS_REAC(I)), SUR_REAC(I), SDTOT, T,
     1                  RA_SUR(I), RB_SUR(I), RE_SUR(I), IERR,
     2                  A_FACTOR, EACT)
        END IF
        VELOC_CRT = 1.0/(1.0-0.5*STK)
        EFF_FLUX = 3637.6011 * SQRT(T/WT(K_GAS_REAC(I)))
        K_STAR = STK * EFF_FLUX * VELOC_CRT * SDEN_RATIO(I)
        RF = STK * EFF_FLUX * VELOC_CRT / SDTOT**SUR_REAC(I) * RP_SUR(I)
        IF (.NOT. FLAG) THEN
          WRITE(LOUT,9702) T, STK, A_FACTOR, EACT, EFF_FLUX, VELOC_CRT,
     1                     SDEN_RATIO(I), K_STAR, RF
        ELSE
          WRITE(LOUT,9705) T, STK, A_FACTOR, EACT, EFF_FLUX, VELOC_CRT,
     1                     SDEN_RATIO(I), K_STAR, RF
        END IF
1950    CONTINUE
      END IF
C
C Describe reverse surface reaction in terms of sticking coefficient,
C     if possible.
C
      IF (GAS_PROD(I) .EQ. 1) THEN
        IF (NRPP_SUR(I) .LE. 0) THEN
          WRITE(LOUT,9703) SDTOT**(SUR_PROD(I)), SUR_PROD(I),
     1                     2*SUR_PROD(I)
        ELSE
          WRITE(LOUT,9704) SDTOT**(SUR_PROD(I)), SUR_PROD(I),
     1                     2*SUR_PROD(I)
        END IF
        IF (RP_SUR(I) .NE. 1.0) WRITE(LOUT,9706) RP_SUR(I)
        WRITE(LOUT,9701) SUR_RUNITS
        DO 1970 ITEMP = 1, TTABLE_NUM
        T = TTABLE(ITEMP)
        CALL CALCRSTA(SKWRK, ISKWRK, I, STK, WT(K_GAS_REAC(I)),
     1                SUR_REAC(I), WT(K_GAS_PROD(I)),SUR_PROD(I),
     2                SDTOT, T, ISTFL_SUR(I), RA_SUR(I), RB_SUR(I),
     3                RE_SUR(I), ACT, SDEN, SEQKC, IERR,
     4                A_FACTOR, EACT)
        IF (STK .LT. 1.0) THEN
          FLAG = .FALSE.
        ELSE
          FLAG = .TRUE.
        END IF
        VELOC_CRT = 1.0/(1.0-0.5*STK)
        EFF_FLUX = 3637.6011 * SQRT(T/WT(K_GAS_PROD(I)))
        K_STAR = STK * EFF_FLUX * VELOC_CRT * SDEN_RATIO(I)
        RF = STK * EFF_FLUX * VELOC_CRT / SDTOT**SUR_PROD(I) * RP_SUR(I)
        IF (.NOT. FLAG) THEN
          WRITE(LOUT,9702) T, STK, A_FACTOR, EACT, EFF_FLUX, VELOC_CRT,
     1                     SDEN_RATIO(I), K_STAR, RF
        ELSE
          WRITE(LOUT,9705) T, STK, A_FACTOR, EACT, EFF_FLUX, VELOC_CRT,
     1                     SDEN_RATIO(I), K_STAR, RF
        END IF
1970    CONTINUE
      END IF
      WRITE(LOUT,5000)
C
C End of loop over 'IF (LSTCK(I)) THEN;
C
      END IF
C
C End of Conditional statements for each surface reaction
C
      END IF
2000  CONTINUE
      END IF
C
C-----------------------------------------------------------------------
C                               BLOCK 8
C                  HANDLE THE CONTINUATION RUN CASE - Normal exits too
C-----------------------------------------------------------------------
*
*         If we want to go to a new problem, set the flag, CNTNUD,
*         which indicates that we are now on a new set of conditions
*         and go to the top of the file
*
      IF (LCNTUE) THEN
        CNTNUD = .TRUE.
        GO TO 100
      END IF
*
      RETURN
*
C-----------------------------------------------------------------------
C                               BLOCK 9
C                           ERROR HANDLING
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                               BLOCK 10
C                           FORMAT STATEMENTS
C-----------------------------------------------------------------------
C
C  Generic Format Statements:
C
5000  FORMAT (132('='))
C
6999  FORMAT (132('=')/
     1 ' GENERAL INFORMATION CONCERNING THE SURFACE CHEMKIN PROBLEM:'/
     2 10X,'Total number of elements declared =', I3,' ',A/
     3 10X,'Total number of species  =',I4/
     4 10X,'Total number of phases   =',I3/
     5 10X,'Total number of gas-phase reactions =',I4/
     6 10X,'Total number of surface-phase reactions =',I4/
     7 10X,'Universal gas constant = ', 1PG13.7,' ergs/(mole*K)'/
     8 10X,'Universal gas constant used for activation energies = ',
     9              1PG13.7,' cals/(mole*K)'/
     1 10x,'Pressure of one standard atmosphere = ',1PG13.7,
     2             ' dynes/cm**2'/)
7000  FORMAT (10X,'GAS phase (Always phase # 1, with the name, "GAS"):'/
     1        15X,'Number of species =',I4/
     2        15X,'Number of surface reactions where the # ',
     1 'of gas products is different than the # of gas reactants =',I3/
     3        15X,'Number of elements in the phase =',I3,' ',A)
7002  FORMAT (/10X,'SURFACE PHASES:'/
     1  15X,'Phase    Name of        SS Site Density   Number of  '
     1     ,'Number of  Site_changing    Elements:'/
     1  15X,'Number    Phase           (moles/cm**2)    Species   '
     1     ,'Elements   Surf_rxns'
     1  /15X,90('-'))
7003  FORMAT (15X,I3,6X,A16,' ',1PE12.4,'      ',I4,'      ',I4,
     1          '       ',I4,'        ',A)
7008  FORMAT (15X,90('-')/
     1    15X,'Tot SS Site Dens= ',1PE12.4,
     1        ' (mole/cm**2) (used in sticking coefficient express)')
7004  FORMAT (/10X,'SURFACE PHASES:'/
     2    15X,'There are no surface phases')
7005  FORMAT (/10X,'BULK PHASES:'/
     1  15X,'Phase    Name of              Number of   '
     1     ,'Number of   Mole_changing    Elements:'/
     1  15X,'Number    Phase               Species    '
     1     ,' Elements    Surf_rxns'
     1  /15X,90('-'))
7006  FORMAT (15X,I3,6X,A16,'      ',I4,'       ',I4,'         ',
     1        I4,'         ',A)
7007  FORMAT (/10X,'BULK PHASES:'/
     2    15X,'There are no bulk phases')
C
8000  FORMAT (1X,80('-'))
8005  FORMAT (1X,'THERMO TABLE FOR MOLECULE "',A,'" IN PHASE "',A,'"')
8006  FORMAT (20X,'Electronic charge =', I3)
8007  FORMAT (20X,'Elemental Composition:')
8008  FORMAT (25X,A,':',I3)
8009  FORMAT (20X,'Heat of Formation at 298 = ',F10.3
     1          , ' kcal/mole'/
     1  ,20X,'Molecular Weight = ',G13.5,' gm/mole'/1X,80('-')/
     1,' |     Temp  |   (H-H298)          (G-H298)     '
     1,'     Cp              S          |'/
     1,' |      (K)  |  (kcal/mole)      (kcal/mole)    '
     1,'(cal/mole*K)    (cal/mole*K)    |'
     1,/1X,'|', 11('-'),'|', 66('-'),'|')
8010  FORMAT (1X,'|',F10.2,1X,'|',1X,1PG13.5,3X,G13.5,3X,G13.5,
     1        3X,G13.5,'    |')
8011  FORMAT (25X,'There are no elements in this species!')
8012  FORMAT (20X,'Number of surface sites occupied by the species ='
     1        ,I2)
8013  FORMAT (20X,'No bulk density was input for this species')
8014  FORMAT (20X,'Bulk Density = ',1PG12.5,' gm/cm**3')
8015  FORMAT (20X,'Overall, this is the ',I3,
     1          'th species in the mechanism'/,
     1        20X,'It is the ',I3,'th species in phase ',A)
8016  FORMAT (
     1   20X,'L-J Potential well depth = ', 1PG11.3,' K'/
     1   20X,'L-J collision diameter   = ', 1PG11.3,' Angstroms'/
     1   20X,'Dipole Moment = ', 1PG11.3,' Debye'/
     1   20X,'Polarizability = ', 1PG11.3,' Angstroms**3'/
     1   20X,'Rotational Collision number at 298K = ',1PG11.3)
8017  FORMAT (20X,'This molecule is consists of a single atom')
8018  FORMAT (20X,'This molecule is linear')
8019  FORMAT (20X,'This molecule is non-linear')
8020  FORMAT (20X,'Heat of Formation at 298 = ',F10.3
     1          , ' kcal/mole'/
     1  ,20X,'Molecular Weight = ',G13.5,' gm/mole'/1X,128('-')/
     1 ' |     Temp  |   (H-H298)          (G-H298)     ',
     1 '     Cp              S       |  Viscosity  Therm_Cond ',
     1   A26,'|'/
     1 ' |      (K)  |  (kcal/mole)    (kcal/mole)    ',
     1 '(cal/mole*K)    (cal/mole*K)   | (gm/cm*sec) (erg/cm*sec*K) ',
     1 '     (cm**2/sec)      |'/
     1 1X,'|', 11('-'),'|', 63('-'),'|', 51('-') )
8021  FORMAT (1X,'|',F10.2, 1X, '|', 1X, 1PG13.5, 3X ,G13.5, 3X, G13.5,
     1        3X, G13.5, ' |', 1PG11.3, 2X, 1PG11.3, 6X, 1PG11.3,9X,'|')
8022  FORMAT (1X, 128('-')/
     1        30X,'[Pressure for binary diffusion coeff. calc. = ',
     1        4PG11.3,' torr]')
8023  FORMAT (20X,'Activity (bath gas dependent) = ',1PG12.5)
C
C Table for the short description of gas-phase and surface-phase
C reactions
C
8200  FORMAT (/132('=')/
     1    10X,'SHORT DESCRIPTION OF GAS-PHASE REACTIONS'//
     1    20X,'Number                      Description        '
     1       ,'    Gas_Mole         Gas_Mole'/
     1    20X,'                                               '
     1       ,'     Change          Reactants'/
     1    20X,80('-'))
8201  FORMAT (20X,I3,'.   ',A40,'      ', I3,13X,I3)
8202  FORMAT (/132('=')/
     1    10X,'SHORT DESCRIPTION OF SURFACE-PHASE REACTIONS'//
     1    20X,'Number              Description                '
     1       ,'Gas Mole  Surf Mole   Bulk Mole    Surf_Site'/
     1    20X,'                                               '
     1       ,' Change    Change      Change       Change'/
     1    20X,90('-'))
8203  FORMAT (20X,I3,'.   ',A40,2X,I3,8X,I3,8X,I3,10X,I3)
8204  FORMAT (20X,80('-'))
8205  FORMAT (27X, A40)
8206  FORMAT (20X,90('-'))
C
8209  FORMAT (/132('=')/
     1   10X,'SUMMARY OF STANDARD STATE THERMODYNAMICS FUNCTIONS FOR ',
     1       'SURFACE-PHASE REACTIONS:'//
     1    20X,'Bath Gas Temperture = ', 1PG11.3//
     1    10X,'Number                      Description        '
     1       ,'    Delta_G        Delta_H        Delta_S'/
     1    10X,'                                               '
     1       ,'  (kcal/mole)    (kcal/mole)    (cal/moleK)'/
     1    10X,90('-'))
8210  FORMAT (/132('=')/
     1   10X,'SUMMARY OF STANDARD STATE THERMODYNAMICS FUNCTIONS FOR ',
     1       'GAS-PHASE REACTIONS:'//
     1    20X,'Bath Gas Temperture = ', 1PG11.3//
     1    10X,'Number                      Description        '
     1       ,'    Delta_G        Delta_H        Delta_S'/
     1    10X,'                                               '
     1       ,'  (kcal/mole)    (kcal/mole)    (cal/moleK)'/
     1    10X,90('-'))
8211  FORMAT (10X,I3,'.   ',A40,3X,1PG12.5,3X,1PG12.5,3X,1PG12.5,3X)
8212  FORMAT (10X,7X,A40)
8213  FORMAT (10X,90('-'))
C
8250  FORMAT (/132('=')//10X,'SUMMARY OF SPECIES IN THE MECHANISM '
     1        ,' with a DESCRIPTION OF BATH GAS COMPOSITION:'//
     1        20X,'Total pressure = ', 4PG11.3, ' torr'/
     1        20X,'Temperature (where needed) = ',4PG11.3,' Kelvin'/
     1        20X,'Carrier Gas (used in diff. calcs) = ',A/
     1        20X,'Major Gas Species (used in nondim calcs) = ',A//
     1     20X,'Number   Name         Mole_fraction      Concentration'/
     1        20X,80('-'))
8251  FORMAT (20X, I3,'.    ',A16,'  ',G11.4,'  ',G12.4,' mole/cm**3')
8252  FORMAT (20X, I3,'.    ',A16,'  ',G11.4,'  ',G12.4,' mole/cm**2')
8253  FORMAT (20X, I3,'.    ',A16,'  ',G11.4,'  ',' activity(unit1ess)')
C
8260  FORMAT (132('=')//
     1 10X,'SUMMARY OF STANDARD STATE THERMODYNAMIC FUNCTIONS FOR ',
     1        'SPECIES AT BATH GAS CONDITIONS:'//
     1        20X,'Bath Gas Temperature = ',4PG11.4,' Kelvin'//
     1 10X,'Number   Name',13X,'H(298 K)    H(T_bath)     Cp(T_bath) ',
     1 '    S(T_bath)     G(T_bath) '/
     1 10X,'           ',13X,'(kcal/mole)  (kcal/mole)   (cal/moleK) ',
     1 '   (cal/moleK)   (kcal/mole)'/
     1        10X,100('-'))
8261  FORMAT (10X,I3,'.    ',A16,' ', 1PG12.5,2X, G12.5,2X,
     1        G12.5,2X, G12.5,2X, G12.5)
8262  FORMAT (132('='))
C
C
C  TABLE OF NON-DIMENSIONAL GAS REACTION RATE CONSTANTS
C    AT BATH GAS CONDITIONS
C
8349  FORMAT (/132('=')/
     1    10X,'NON-DIMENSIONAL GAS REACTION RATE CONSTANTS'/
     1    10X,'         AT THE BATH GAS CONDITIONS'/
     1    20X,'Total Pressure =  ',4PG11.3,' torr'/
     1    20X,'Temperature    = ',4PG11.3,' Kelvin'//
     1    10X,'Number              Description          ',
     1    12X,'k_star      k_star_rev    |  Gas_Da_For ',
     1        '   Gas_Da_Rev'/
     1    56X,'(mole/cm**3 sec) (mole/cm**3 sec)|'/
     1    10X,79('-'),'|',31('-'))
8351  FORMAT (10X,I3,'.   ',A40,2X,1PG11.3,2X,'[',1PG11.3,']',
     1        4X,'|',2X,1PG11.3,'  [',1PG11.3,']')
8352  FORMAT (10X,I3,'.   ',A40,2X,1PG11.3,3X,1PG11.3,
     1        5X,'|',2X,1PG11.3, 3X,  1PG11.3)
8353  FORMAT (17X, A40, 32X,'|')
8354  FORMAT (10X,111('-'))
8355  FORMAT (30X,
     1 '[  ] indicates that this reaction is not in mechanism')
8356  FORMAT (10X,'NOTE ON THE ABSOLUTE NUMBERS IN THIS TABLE:'/
     1  20X,'The rate rate constants (mole/cm**3*sec) should be ',
     1  'compared to rate of mass transport in order'/
     1  20X,'to characterize their values as being fast or slow.',
     1  '  The nondimensionalization of the mass transport '/20X,
     1  'involves the following multiplicative factor, which also',
     1  ' has the units of mole/cm**3*sec:'/
     1  30X,'Total_Concentration * Diffusivity / Length_scale**2'/
     1  20X,'Using the binary diffusion coefficient between ',A/
     1  20X,'and ',A,', the following factors are calculated',
     1      ' at bath gas conditions:'/
     1  30X,'Total Concentration          = ',1PG10.3,' mole/cm**3'/
     1  30X,'Binary Diffusion Coefficient = ',1PG10.3,' cm**2/sec'/
     1  30X,'Length scale                 = ',1PG9.3, ' cm'/
     1  20X,'Therefore, the non-dimensionalization factor for gas ',
     1      'reactions becomes:'/
     1  30X,'Conc * Diff / Length**2      = ',1PG10.3,' mole/cm**3*sec'/
     1  20X,'Note that this number is independent of pressure'/)
C
C  TABLE OF NON-DIMENSIONAL SURFACE REACTION RATE CONSTANTS
C    AT BATH GAS CONDITIONS
C
8359  FORMAT (/132('=')/
     1    10X,'NON-DIMENSIONAL SURFACE REACTION RATE CONSTANTS'/
     1    10X,'         AT THE BATH GAS CONDITIONS'/
     1    20X,'Total Pressure =  ',4PG11.3,' torr'/
     1    20X,'Temperature    = ',4PG11.3,' Kelvin'//
     1    10X,'Number              Description          ',
     1    12X,'k_star      k_star_rev    | Surf_Da_For ',
     1        '  Surf_Da_Rev'/
     1    56X,'(mole/cm**2 sec) (mole/cm**2 sec)|'/
     1    10X,79('-'),'|',31('-'))
8361  FORMAT (10X,I3,'.   ',A40,2X,1PG11.3,2X,'[',1PG11.3,']',
     1        4X,'|',2X,1PG11.3,'  [',1PG11.3,']')
8362  FORMAT (10X,I3,'.   ',A40,2X,1PG11.3,3X,1PG11.3,
     1        5X,'|',2X,1PG11.3, 3X,  1PG11.3)
8363  FORMAT (17X, A40, 32X,'|')
8364  FORMAT (10X,111('-'))
8365  FORMAT (30X,
     1 '[  ] indicates that this reaction is not in mechanism')
8366  FORMAT (10X,'NOTE ON THE ABSOLUTE NUMBERS IN THIS TABLE:'/
     1  20X,'The rate rate constants (mole/cm**2*sec) should be ',
     1  'compared to rate of mass transport to the surface in order'/
     1  20X,'to characterize their values as being fast or slow.',
     1  '  The nondimensionalization of the mass transport '/20X,
     1  'involves the following multiplicative factor, which also',
     1  ' has the units of mole/cm**2*sec:'/
     1  30X,'Total_Concentration * Diffusivity / Length_scale'/
     1  20X,'Using the binary diffusion coefficient between ',A/
     1  20X,'and ',A,', the following factors are calculated',
     1      ' at bath gas conditions:'/
     1  30X,'Total Concentration          = ',1PG10.3,' mole/cm**3'/
     1  30X,'Binary Diffusion Coefficient = ',1PG10.3,' cm**2/sec'/
     1  30X,'Length scale                 = ',1PG9.3, ' cm'/
     1  20X,'Therefore, the non-dimensionalization factor for surface ',
     1      'reactions becomes:'/
     1  30X,'Conc * Diff / Length         = ',1PG10.3,' mole/cm**2*sec'/
     1  20X,'Note that this number is independent of pressure'/)
C
C  TABLE FOR LONG DESCRIPTION OF REACTION RATE AS A FUNCTION OF T
C
8381  FORMAT(/1X,'Gas Reaction # ',I3,10X,A64/)
8385  FORMAT(20X,'This reaction does not have any third body effects')
8386  FORMAT(20X,'This reaction does have third body effects, but'
     1  ,'  no enhanced third body efficiencies were input')
8387  FORMAT(20X,'This reaction does have third body effects. '
     1     ,I2,' modified enhanced third body efficiencies were input')
8388  FORMAT(30X,'Species "',A,
     1   '", modified enhanced third body efficiency'
     1  ,' for the reaction = ',G11.4)
C
8400  FORMAT(20X,'Change in moles in the reaction =',I2)
8410  FORMAT(20X,'Perturbation factor for rate contant =',E12.5)
8510  FORMAT(20X,'This is a reversible reaction, having ',I1
     1          ,' reactant species'/
     1       52X,'and ',I1,' product species',A)
8520  FORMAT(20X,'This is a irreversible reaction, having ',I1
     1          ,' reactant species'/
     1       54X,'and ',I1,' product species',A)
C
8561  FORMAT(20X,'Reaction has a fall-off behavior with a ',
     1  'Lindeman function form:'/
     2 25X,'klow ',A,'= ',1PG11.3,' T**(',1PG11.4,') exp( -',1PG10.3,
     1           ' kcal/mole / RT)'/)
8562  FORMAT(20X,'Reaction has a fall-off behavior with a ',
     1  'Lindeman function form:'/
     2 25X,'klow ',A,'= ',1PG11.3,' T**(',G11.4,') exp( -',1PG10.3,
     1           ' kcal/mole / RT)'/
     1 30X,'a = ', 1PG11.3/
     1 30X,'b = ', 1PG11.3,' Kelvin'/
     1 30X,'c = ', 1PG11.3,' Kelvin'/
     1 30X,'d = ', 1PG11.3/
     1 30X,'e = ', 1PG11.3)
8563  FORMAT(20X,'Reaction has a fall-off behavior with a ',
     1  '6 parameter Troe function form:'/
     2 25X,'klow ',A,'= ',1PG11.3,' T**(',1PG11.4,') exp( -',3PG10.3,
     1           ' kcal/mole / RT)'/
     1 30X,'a = ', 1PG12.4/
     1 30X,'T*** = ',4PG11.3,' Kelvin'/
     1 30X,'T*   = ',4PG11.3,' Kelvin' )
8564  FORMAT(20X,'Reaction has a fall-off behavior with a ',
     1  '7 parameter Troe function form:'/
     2 25X,'klow ',A,'= ',1PG11.3,' T**(',1PG11.4,') exp( -',3PG10.3,
     1           ' kcal/mole / RT)'/
     1 30X,'a = ', 1PG12.4/
     1 30X,'T*** = ',4PG11.3,' Kelvin'/
     1 30X,'T*   = ',4PG11.3,' Kelvin'/
     1 30X,'T**  = ',4PG11.3,' Kelvin')
8570  FORMAT(25X,'Third body in fall-off calculation is a specific',
     1   ' species, ',A)
C
8600  FORMAT(132('-')/
     1,30X,'HIGH PRESSURE GAS REACTION RATE CONSTANTS AS A FUNCTION',
     1  ' OF TEMPERATURE'//
     1 50X,'(note: reverse rate constant is not in mechanism)',
     1    8X, '|UnifDimensnl Rate_Const|'/
     1 '   T     |    k      A_factor    Ea     |',
     1 '  DeltaG     DeltaH     DeltaS   ',
     1 '| k_rev   A_factor_rev  Ea_rev   | k_star    k_star_rev  |')
8601  FORMAT(132('-')/
     1,30X,'HIGH PRESSURE GAS REACTION RATE CONSTANTS AS A FUNCTION',
     1  ' OF TEMPERATURE'//
     1 107X,'|UnifDimensnl Rate_Const|'/
     1 '   T     |    k      A_factor    Ea     |',
     1 '  DeltaG     DeltaH     DeltaS   ',
     1 '| k_rev   A_factor_rev  Ea_rev   | k_star    k_star_rev  |')
8610  FORMAT('  (K)',A22,'  (kcal/mol)'
     1 ,' |(kcal/mol)(kcal/mole) (cal/moleK)|',A22,'(kcal/mol)|'
     1 ,'   (mole/cm**3*sec)    |'
     1 /1X,8('-'),'|',30('-'),'|',33('-'),'|',32('-'),'|',23('-'),'|'/
     1     9X    ,'|',30X    ,'|',33X    ,'|',32X    ,'|',23X    ,'|')
8700  FORMAT(1X,F7.2,' |',1PE9.2,1PG9.2, 1X, 1PG11.4,'|',
     1       1PG11.4,1PG11.4,1PG11.4,'|', 1PE10.2,1PG11.2,1PG11.4,'|',
     1        1PG11.2,1PG11.2,' |')
C
C         TABLE OF FALL-OFF BEHAVIOR AS A FUNCTION OF TEMPERATURE.
C
8705  FORMAT(132('-')/,25X,
     1'FALL-OFF BEHAVIOR AS A FUNCTION OF TEMPERATURE: Reaction # ',
     1    I4/
     1 35X,'Bath Gas Pressure = ', 4PG11.3,' torr'//
     1 ' Temperature   k         A_factor    Ea     |  '
     1,' k/kinf     klow     Reduc_Press   FC   Eff_Conc'
     1,'   |   k_rev   A_factor_rev  Ea_rev'
     1)
8710  FORMAT('   (K)   ',A22,' (kcal/mole) |'
     1 ,'       ',A22,'           (mole/cm**3) |',A22,'(kcal/mole)'
     1   /1X,43('-'),'|',53('-'),'|',33('-'))
8720  FORMAT(1X,4PG10.3, 1PG11.4, 1PG11.3, 1PG11.4,
     1        '|', 1PG11.3, 1PG11.3, 1PG11.3, 1PG10.3, 1PG9.2,
     1        ' |', 1PG11.4, 1PG11.3, 1PG11.4)
C
C         TABLE OF FALL-OFF BEHAVIOR AS A FUNCTION OF PRESSURE.
C
8725  FORMAT(132('-')/
     1,25X,'FALL-OFF BEHAVIOR AS A FUNCTION OF PRESSURE: Reaction # ',
     1    I4/
     1 35X,'Bath Gas Temperature = ', 4PG11.3,' K'/
     1 35X,'Low Pressure Limiting Reaction Rate = klow = ', 1PG11.4,
     1     ' ',A22//
     1 '   Pressure    k         A_factor    Ea     |  '
     1,' k/kinf klow*Eff_Conc Reduc_Press  FC   Eff_Conc'
     1,'   |   k_rev   A_factor_rev  Ea_rev'
     1)
8730  FORMAT('  (torr) ',A22,' (kcal/mole) |'
     1 ,'       ',A22,'           (mole/cm**3) |',A22,'(kcal/mole)'
     1   /1X,43('-'),'|',53('-'),'|',33('-'))
8740  FORMAT(1X,4PG10.3, 1PG11.4, 1PG11.3, 1PG11.4,
     1        '|', 1PG11.3, 1PG11.3, 1PG11.3, 1PG10.3, 1PG9.2,
     1        ' |', 1PG11.4, 1PG11.3, 1PG11.4)
C
C   TABLE OF THIRD BODY REACTION RATES - [M] LUMPED WITH RATE CONSTANT
C
8745  FORMAT(132('-')/,25X,
     1'ANALYSIS OF THIRD BODY REACTIONS: LUMPING [M] WITH RATE CNST'/
     1    30X,'Reaction # ',I3/
     1 35X,'Bath Gas Pressure = ', 4PG11.3,' torr'//
     1 ' Temperature   k         A_factor    Ea     |  '
     1,' Concentration         C_eff'
     1,'   |   k_rev   A_factor_rev  Ea_rev'
     1)
8750  FORMAT('   (K)   ',A22,' (kcal/mole) |'
     1 ,'   (mole/cm**3)     (mole/cm**3) |',A22,'(kcal/mole)  |'
     1   /1X,43('-'),'|',33('-'),'|',35('-'),'|')
8760  FORMAT(1X, F10.2, 1PG11.4, 1PG11.3, 1PG11.4,
     1        '|   ',  1PG11.3, 5X, 1PG11.3, 3X,
     1        '|  ', 1PG11.4, 1PG11.3, 1PG11.4,'|')
C
C
C      SECTION ON THE SURFACE REACTION RATE AS A FUNCTION OF TEMPERATURE
C
9301  FORMAT(/1X,'Surface Reaction # ',I3,10X,A64/)
9400  FORMAT(20X,'Change in gas     moles in the reaction =',I2/
     1       20X,'Change in surface moles in the reaction =',I2/
     1       20X,'Change in bulk    moles in the reaction =',I2)
9420  FORMAT(/20X,'This reaction has ',I2,' species whose surface ',
     1         'coverage modify the rate constant '/
     1       25X,'Each of these species has three parameters that',
     1       ' multiplicatively modify the rate constant as follows:'/
     1       30X,' k_prime = k * 10**(Z_k*nu_ki) * Z_k**mu_ki * ',
     1          ' exp[ - eps_ki*Z_k / Rc*T ]'/
     1       35X,'where Z_k = Site Fraction of species k')
9421  FORMAT(25X,'Species = ', A/
     1       30X,'nu_ki  = ', 1PG11.4,' (cm**2/mole)'/
     1       30X,'mu_ki  = ', 1PG11.4/
     1       30X,'eps_ki = ', 1PG11.4,' (cal*cm**2/mole**2)'/)
9510  FORMAT(20X,'This is a reversible surface reaction, having '
     1          ,'the following types of reactant species:'/
     1       35X,I1,' gas-phase species'/
     1       35X,I1,' surface-phase species'/
     1       35X,I1,' bulk-phase species'/
     1       25X,'and the following types  of product species:'/
     1       35X,I1,' gas-phase species'/
     1       35X,I1,' surface-phase species'/
     1       35X,I1,' bulk-phase species')
9520  FORMAT(20X,'This is an irreversible surface reaction, having '
     1          ,'the following types of reactant species:'/
     1       35X,I1,' gas-phase species'/
     1       35X,I1,' surface-phase species'/
     1       35X,I1,' bulk-phase species'/
     1       25X,'and the following types  of product species:'/
     1       35X,I1,' gas-phase species'/
     1       35X,I1,' surface-phase species'/
     1       35X,I1,' bulk-phase species')
9528  FORMAT(20X,'The reaction rate constant was input via ',
     1  'a sticking coefficient in the interpretor input file')
9529  FORMAT(20X,'Motz-Wise Correction factor is used')
9530  FORMAT(20X,'Motz-Wise Correction factor is not used')
9531  FORMAT(20X,'A sticking coefficient was not used though ',
     1   'the forward reaction could have used one')
9532  FORMAT(20X,'A sticking coefficient was not used though ',
     1   'the forward reaction could have used one')
9534  FORMAT(25X,'The forward rate constant can be fit to the ',
     1      'following sticking coefficient expression:'/
     1       25X,'(WARNING: STICKING COEFF. IS GREATER THAN ONE)')
9535  FORMAT(25X,'The forward rate constant can be fit to the ',
     1      'following sticking coefficient expression:')
9536  FORMAT(
     1 25X,'It can be fit to the following general rate constant form:')
9537  FORMAT(20X,'The reverse rate constant can be fit to the ',
     1      'following form:')
9541  FORMAT(20X,'Even though this reaction is IRREVERSIBLE, the '
     1        ,'reverse rate constant will also be analysed:'/
     1       20X,' The reverse rate constant can be fit to the ',
     1      'following form:')
9538  FORMAT(25X,'The reverse rate constant can also be expressed in a',
     1      ' sticking coefficient form:')
9539  FORMAT(30X,'WARNING: Sticking coefficient is greater than one ',
     1       'for some temperatures')
9540  FORMAT(30X,'WARNING: Sticking coefficient for reverse rate ',
     1           'constant is greater than one for some temperatures')
C
C       FORMATTING FOR SURFACE REACTION RATE CONSTANT TABLE
C
9599  FORMAT(1X,131('-')/
     1 35X,'FORWARD AND REVERSE SURFACE REACTION RATE CONSTANTS',
     1 21X,'|  Bath Gas Dependent   |'/,
     1 50X,'(note: reverse rate constant is not in mechanism)',
     1    8X, '|UnifDimensnl Rate_Const|'/,
     1 '   T     |    k      A_factor    Ea      |',
     1 '  DeltaG     DeltaH    DeltaS   ',
     1 '| k_rev   A_factor_rev  Ea_rev   | k_star    k_star_rev  |')
9600  FORMAT(1X,131('-')/
     1 35X,'FORWARD AND REVERSE SURFACE REACTION RATE CONSTANTS',
     1 21X,'|  Bath Gas Dependent   |',/
     1 107X,'|UnifDimensnl Rate_Const|',/
     1 '   T     |    k      A_factor    Ea      |',
     1 '  DeltaG     DeltaH    DeltaS   ',
     1 '| k_rev   A_factor_rev  Ea_rev   | k_star    k_star_rev  |')
9601  FORMAT('  (K)',A22,'  (kcal/mol)'
     1 ,'  |(kcal/mol)(kcal/mole)(cal/moleK)|',A22,'(kcal/mol)|'
     1 ,'   (mole/cm**2*sec)    |'
     1 /1X,8('-'),'|',31('-'),'|',32('-'),'|',32('-'),'|',23('-'),'|'/
     1     9X    ,'|',31X    ,'|',32X    ,'|',32X    ,'|',23X    ,'|')
9602  FORMAT(1X,F7.2,' |',1PG10.2,1PG10.2,1X,G10.4,'|',
     1       G11.4,G11.4,G10.4,'|', 1PG10.2,1PG11.2,G11.4,'|',
     1        1PG11.2,1PG11.2,' |')
C
9610  FORMAT(1X,131('-')//
     1 35X,'ANALYSIS OF FORWARD AND REVERSE COVERAGE DEPENDENT'/
     1 35X,'  SURFACE RATE CONSTANTS AT BATH GAS CONDITIONS:'//
     1 '   T       k_prime       A_factor    Ea     |   '
     1,'      k          Cov_fac           k_rev      | '
     1,'krev_prime   A_factor_rev   Ea_rev' )
9611  FORMAT('  (K) ',A22,'    (kcal/mole) |'
     1 ,A22,'(cgs)',A22,'|',A22,'    (kcal/mole)'
     1 /1X,43('-'),'|',49('-'),'|',37('-')/
     1  1X,43(' '),'|',49(' '),'|' )
9612  FORMAT(1X,F7.2,1X,1PG12.3, G12.3, G11.4,'|',
     1      2X, G12.4, 4X, G12.4, 4X, G12.4,'   | ',3G12.4)
C
C TABLE OF THE BREAKDOWN OF THE FORWARD AND REVERSE
C REACTION'S STICKING COEFFICIENT
C
9700  FORMAT(1X,131('-')//
     1,35X,'BREAKDOWN OF FORWARD REACTION''S STICKING COEFFICIENT'//
     1,40X,'Surface site density devisor = ',1PE12.4,' mole**',I1,
     1      '/cm**',I1/)
9706  FORMAT(40X,'Reaction Rate Perturbation factor = ',1PG11.3)
9701  FORMAT(
     1  '   T        Stck_Coeff    A_factor     Ea      '
     1 ,' Eff_Veloc   Veloc_Corr   Sden_Ratio   '
     1 ,'   k*          k'/
     1  '   (K)      (unitless)                        '
     1 ,' (cm/sec)                              '
     1 ,'(cm/sec)', A,   /1X,131('-')/)
9702  FORMAT(1X,F8.2,1X,7(1PG12.4),3X,G12.4)
9703  FORMAT(1X,131('-')//
     1,35X,'BREAKDOWN OF REVERSE REACTION''S STICKING COEFFICIENT'/
     1,40X,'(note: reverse reaction is not in the mechanism)'/
     1,40X,'Surface site density devisor = ',1PE12.4,' mole**',I1,
     1      '/cm**',I1/)
9704  FORMAT(1X,131('-')//
     1,35X,'BREAKDOWN OF REVERSE REACTION''S STICKING COEFFICIENT'//
     1,40X,'Surface site density devisor = ',1PE12.4,' mole**',I1,
     1      '/cm**',I1/)
9705  FORMAT(1X,F8.2,'*',1PG11.4,'*',6G12.4,3X,G12.4)
C
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE BACKSYM ( RNAM, BRNAM)
C
C Reverse the character string representing the reaction
C
      CHARACTER*(*) RNAM, BRNAM
      CHARACTER*8  ARROW
      INTEGER LEQ, LLEFT, LRIGHT, LTOT, FLEN, RLEN, IARROW
C
      LEQ    = INDEX(RNAM,'=')
      LLEFT  = INDEX(RNAM,'<')
      LRIGHT = INDEX(RNAM,'>')
      LTOT   = INDEX(RNAM,' ')-1
      IF ( LLEFT .GT. 0) THEN
*          Reaction contains a < symbol, find the length of the
*          forward reaction string.
         FLEN = LLEFT - 1
      ELSE
         FLEN = LEQ - 1
      ENDIF
      IF ( LRIGHT .GT. LEQ) THEN
C        REACTION CONTAINS A > SYMBOL
C        FIND LENGTH OF REVERSE REACTION
         RLEN = LTOT - LRIGHT
      ELSE
         RLEN = LTOT - LEQ
      ENDIF
C     MAKE REVERSE ARROW SYMBOL
      IARROW = 0
      IF ( LRIGHT .GT. 0) THEN
         IARROW = IARROW + 1
         ARROW(IARROW:IARROW) = '<'
      ENDIF
      IARROW = IARROW + 1
      ARROW(IARROW:IARROW) = '='
      IF ( LLEFT .GT. 0) THEN
         IARROW = IARROW + 1
         ARROW(IARROW:IARROW) = '>'
      ENDIF
      BRNAM = ' '
      BRNAM(1:RLEN) = RNAM(LTOT-RLEN+1:LTOT)
      BRNAM(RLEN+1:RLEN+IARROW) = ARROW(1:IARROW)
      BRNAM(RLEN+IARROW+1:RLEN+IARROW+FLEN) = RNAM(1:FLEN)
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE HIGH_PRESS(RF, A_FACTOR, EACT, RA, RB, RE, RP, T)
C
C
C    Calculate the reaction rate, A_factor, and activation energy
C given a standard 3 parameter Arrhenius expression with a
C 4th parameter as a perturbation factor.  EACT is returned
C with units in kcal/mole.
C
C note it is assumed that RE has units of Kelvin
*
      INCLUDE 'surftherm.h'
*
C Dummy Variables
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     $                 A_FACTOR, EACT, RA, RB, RE, T, RF, RP
*
      IF (RE .EQ. 0.0 .AND. RB .EQ. 0.0) THEN
        EACT = 0.0
        A_FACTOR = RP * RA
        RF = A_FACTOR
      ELSE
        EACT = (RE + RB * T) * RKCAL
        A_FACTOR = RP * RA * T **(RB) * EXP(RB)
        RF = A_FACTOR * EXP(-EACT / (RKCAL * T))
      END IF
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CALC_ARR(K1, T1, K2, T2, A_FACTOR, EACT)
C
C    This subroutine will calculate Arrhenius parameters, given
C a reaction rate at two different tmperatures.  EACT is returned with
C units of kcal/mole.
C    It pays special
C attention to avoiding an overflow or underflow condition.  It uses
C the MACH chemkin common blocks. It expects these to be initialized.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Variables
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     $       K1, T1, K2, T2, A_FACTOR, EACT
C
C Local Variables
C
C*****B-2) double precision
      DOUBLE PRECISION         TMP
C*****END B-2) double precision
C*****B-1) single precision
C      REAL                    TMP
C*****END B-1) single precision
*
* Check the Arguments
*
      IF (T2 .EQ. T1) THEN
        A_FACTOR = K1
        EACT = 0.0
        RETURN
      END IF
      IF (K1 .LE. 0.0 .AND. K2 .GE. 0.0) THEN
        A_FACTOR = K2
        EACT = 0.0
        RETURN
      END IF
      IF (K2 .LE. 0.0 .AND. K1 .GE. 0.0) THEN
        A_FACTOR = K1
        EACT = 0.0
        RETURN
      END IF
C
      IF (T2 .GT. T1) THEN
        TMP  =    T2 /(MAX(T2-T1,SMALL)) * LOG (K2/K1)
      ELSE
        TMP  = -  T2 /(MAX(T1-T2,SMALL)) * LOG (K2/K1)
      END IF
C
C Check EACT against EXPARG
C
      IF (TMP .GE. 0.0) THEN
        IF (TMP .GT. EXPARG) THEN
          TMP =   EXPARG
        END IF
        EACT = TMP * (RKCAL*T1)
        TMP = EXP(TMP)
        IF ((BIG/TMP) .GT. ABS(K1)) THEN
          A_FACTOR = K1 * TMP
        ELSE
          IF (K1 .GT. 0) THEN
            A_FACTOR = BIG
            EACT = RKCAL * T1 * (LOG(BIG) - LOG(K1))
          ELSE
            A_FACTOR = - BIG
            EACT = RKCAL * T1 * (LOG(BIG) - LOG(-K1))
          END IF
        END IF
      ELSE
        IF ( - TMP .GT. EXPARG) THEN
          TMP = - EXPARG
        END IF
        EACT = TMP * (RKCAL*T1)
        TMP = EXP(TMP)
        IF ((SMALL/MAX(TMP,SMALL)) .LT. ABS(K1)) THEN
          A_FACTOR = K1 * TMP
        ELSE
          IF (K1 .GT. 0.0) THEN
            A_FACTOR = SMALL
            EACT = RKCAL * T1 * (LOG(SMALL) - LOG(K1))
          ELSE
            A_FACTOR = - SMALL
            EACT = RKCAL * T1 * (LOG(SMALL) - LOG(-K1))
          END IF
        END IF
      END IF
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE STKACT(RF, A_FACTOR, EACT, RA, RB, RE, T, IERR)
C
C Calculate the sticking coefficient, A_factor, and activation energy
C given a standard 3 parameter Arrhenius expression .  EACT is returned
C with units in kcal/mole. Mins the sticking coefficient against 1.0
C Returns a positive error flag if raw sticking coefficient is greater
C than one.
C
C note it is assumed that RE has units of Kelvin. EACT is returned
C with units in kcal/mole.
*
      INCLUDE 'surftherm.h'
*
      LOGICAL IERR
C*****B-2) double precision
      DOUBLE PRECISION A_FACTOR, EACT, RA, RB, RE, T, RF
C*****END B-2) double precision
C*****B-1) single precision
C      REAL             A_FACTOR, EACT, RA, RB, RE, T, RF
C*****END B-1) single precision
*
      RF = RA * T**RB * EXP(- RE/T )
      IF (RF .GT. 1.0) THEN
        A_FACTOR = 1.0
        EACT = 0.0
        RF = 1.0
        IERR = .TRUE.
      ELSE
        IERR = .FALSE.
        IF (RE .EQ. 0.0 .AND. RB .EQ. 0.0) THEN
          EACT = 0.0
          A_FACTOR = RA
        ELSE
          EACT = (RE + RB * T) * RKCAL
          A_FACTOR = RA * T **(RB) * EXP(RB)
        END IF
      END IF
      RETURN
      END
*
C----------------------------------------------------------------------C
      INTEGER FUNCTION PHASEID(K, NNPHAS, KFIRST, KLAST)
C
C   Given a species number this integer function returns the phase
C  ID of the species. For
      INTEGER K, NNPHAS, KFIRST(NNPHAS), KLAST(NNPHAS)
      INTEGER IPHASE
      IF (K .LE. 0 ) THEN
        PHASEID = 0
        RETURN
      END IF
      DO 10 IPHASE = 1, NNPHAS
        IF (KFIRST(IPHASE) .LE. K) THEN
          IF (KLAST(IPHASE) .GE. K) THEN
            PHASEID = IPHASE
            RETURN
          END IF
        END IF
10    CONTINUE
      PHASEID = 0
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE BRKLINE(STRING, BROKLINE, MAXLINES, CHATOL, BRKSTRNG,
     1                    BRKLEN, NUMCOND, LINES)
C
C
C   This subroutine breaks up the string, STRING, into a vector of
C strings packed into the character vector, BROKLINE.  The number
C of lines in BROKLINE is returned in the integer variable, LINES.
C Whether STRING is broken up or not depends on the length of the
C string, as determined by a string length tolerance, CHATOL.
C Breakpoints for the string are input via the character vector,
C BRKSTRNG.  There are NUMCOND entries in BRKSTRNG.  Each one can
C have a distinct character length, input from the integer vector
C BRKLEN.  Note, this subroutine does work with ' ' as a
C delimiter.  Line breaking is strictly prioritized according to the
C position of entries in BRKSTRNG.  This program looks for the first
C two occurrences of BRKSTRNG in each line. It then picks the one
C that occurs closest to the middle of the line for the location
C of the next line break.
C
C Dummy Variables
      INTEGER LINES, NUMCOND, BRKLEN(NUMCOND), MAXLINES
      CHARACTER STRING*(*), BROKLINE(MAXLINES)*(*),
     1          BRKSTRNG(NUMCOND)*(*)
      INTEGER CHATOL
C Local Variables
      INTEGER I, J, K, OK, IND, IND1, IND2, LEN_BROK, BREAKPOINTS(0:10),
     1        CHALEN
C Externals
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C-----------------------------------------------------------------------
      OK = 0
      LINES = 1
      LEN_BROK = LEN(BROKLINE(1))
      BREAKPOINTS(0) = 0
      BREAKPOINTS(1) = CKLSCH(STRING)
C
   10 CONTINUE
      IF (LINES .LT. MAXLINES .AND. LINES .LT. 10) THEN
        DO 20 I = OK + 1, LINES
          CHALEN =  BREAKPOINTS(I)-BREAKPOINTS(I-1)
          IF ( (CHALEN .GT. CHATOL  .OR. CHALEN .GT. LEN_BROK) .AND.
     1         CHALEN .GT. 2 ) THEN
            DO 15 J = 1, NUMCOND
              IND1 = INDEX(STRING((BREAKPOINTS(I-1)+2):BREAKPOINTS(I)),
     1                     BRKSTRNG(J)(1:BRKLEN(J))  )
              IF (IND1 .NE. 0) THEN
                IND2 = IND1 + INDEX(
     1              STRING((BREAKPOINTS(I-1)+2+IND1):BREAKPOINTS(I)),
     1                     BRKSTRNG(J)(1:BRKLEN(J))  )
                IF (IND2 .GT. IND1) THEN
                  IF (ABS(IND2-CHALEN/2) .LT. ABS(IND1-CHALEN/2))
     1              THEN
                    IND = IND2
                  ELSE
                    IND = IND1
                  END IF
                ELSE
                  IND = IND1
                END IF
                LINES = LINES + 1
                DO 12 K = LINES, I + 1, -1
                  BREAKPOINTS(K) = BREAKPOINTS(K-1)
   12           CONTINUE
                BREAKPOINTS(I) = IND + BREAKPOINTS(I-1)
                GO TO 10
              END IF
   15       CONTINUE
          ELSE
            OK = OK + 1
          END IF
   20   CONTINUE
      END IF
C
      DO 100 I = 1, LINES
        BROKLINE(I) = STRING((BREAKPOINTS(I-1)+1):BREAKPOINTS(I))
  100 CONTINUE
      IF (LINES .LT. MAXLINES) THEN
        DO 110 I = LINES+1, MAXLINES
          BROKLINE(I) = ' '
  110   CONTINUE
      END IF
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE SRX_UNIT(TSTRING, F_RATE_UNITS, F_RATE_SUR_UNITS, IERR)
C
C    Packs the correct string for display of the units for
C surface reaction rate constants. The string is
C packed center-justified.  Note F_RATE_SUR_UNITS may be a non-integer
C due to coverage dependent surface reaction rates.  In those cases
C it will be truncated to units of 1/100.
C  F_RATE_UNIT     = power of (cm**3/mole) in the rate constant
C  F_RATE_SUR_UNIT = power of (cm**2/mole) in the rate constant
C  IERR is true if there isn't enough room in the string.
C
C The final string will have units of:
C
C (cm**3/mole)**F_RATE_UNIT * (cm**2/mole)**F_RATE_SUR_UNIT / sec)
C
      CHARACTER  TSTRING*(*), FSTRING*128, STRING4*4
      INTEGER    LT, LT2
      LOGICAL    IERR
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     $   F_RATE_UNITS, F_RATE_SUR_UNITS, CM_POW, MOLE_POW
C
C External Variables
C
      INTEGER   CHAR_EXP
      EXTERNAL  CHAR_EXP
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - C
C
      IERR = .FALSE.
*
* Work with a temporary string that we know will have enough room
*
      FSTRING = '('
      LT = 1
*
*       Determine the exponent on the cm and mole units
*
      CM_POW   =  F_RATE_UNITS*3 + F_RATE_SUR_UNITS*2
      MOLE_POW = -F_RATE_UNITS   - F_RATE_SUR_UNITS
*
* Write the numerator
*
      IF (CM_POW .LE. 0.0 .AND. MOLE_POW .LE. 0.0) THEN
        FSTRING = FSTRING(1:LT)//'1 '
        LT = LT + 2
      ELSE IF (MOLE_POW .EQ. 1.0) THEN
        FSTRING = FSTRING(1:LT) // 'mol '
        LT = LT + 4
      ELSE IF (MOLE_POW .GT. 0.0) THEN
        LT2 = CHAR_EXP(STRING4, MOLE_POW)
        IF (LT2 .GT. 2) THEN
          FSTRING = FSTRING(1:LT) // 'm**' // STRING4(1:LT2)
          LT = LT + 3 + LT2
        ELSE
          FSTRING = FSTRING(1:LT) // 'mol**' // STRING4(1:LT2)
          LT = LT + 5 + LT2
        ENDIF
      END IF
*
* Add a space, if necessary
*
      IF (FSTRING(LT:LT) .NE. ' ' .AND. FSTRING(LT:LT) .NE. '(') THEN
        LT = LT + 1
        FSTRING(LT:LT) =  ' '
      ENDIF
      IF (CM_POW .EQ. 1.0) THEN
        FSTRING = FSTRING(1:LT)//'cm '
        LT = LT + 3
      ELSE IF (CM_POW .GT. 0.0) THEN
        LT2 = CHAR_EXP(STRING4, CM_POW)
        FSTRING = FSTRING(1:LT) // 'cm**' // STRING4(1:LT2)
        LT = LT + 4 + LT2
      END IF
*
* Add the division symbol
*
      FSTRING = FSTRING(1:LT) // '/'
      LT = LT + 1
*
* Write the Denominator
*
      IF (MOLE_POW .EQ. -1.0) THEN
        FSTRING = FSTRING(1:LT)//'mol '
        LT = LT + 4
      ELSE IF (MOLE_POW .LT. 0.0) THEN
        LT2 = CHAR_EXP(STRING4, -MOLE_POW)
        FSTRING = FSTRING(1:LT) // 'mol**' // STRING4
        LT = LT + 5 + LT2
      ENDIF
*
* Add a space, if necessary
*
      IF (FSTRING(LT:LT) .NE. ' ' .AND. FSTRING(LT:LT) .NE. '(') THEN
        LT = LT + 1
        FSTRING(LT:LT) =  ' '
      ENDIF
      IF (CM_POW .EQ. -1.0) THEN
        FSTRING = FSTRING(1:LT) // 'cm '
        LT = LT + 3
      ELSE IF (CM_POW .LT. 0.0) THEN
        LT2 = CHAR_EXP(STRING4, -CM_POW)
        FSTRING = FSTRING(1:LT) // 'cm**' // STRING4
        LT = LT + 4 + LT2
      END IF
*
* Add a space, if necessary
*
      IF (FSTRING(LT:LT) .NE. ' ' .AND. FSTRING(LT:LT) .NE. '(') THEN
        LT = LT + 1
        FSTRING(LT:LT) =  ' '
      ENDIF
*
* Add sec units -> always assumed to be sec-1
*
      FSTRING = FSTRING(1:LT) // 'sec)'
      LT = LT + 4
*
* Check to see if TSTRING was long enough to hold the string just
* created. If it isn't, then we'll just cram what we can into it
*
      IF (LT .GT. LEN(TSTRING)) THEN
        DO 900 LT2 = 1, LEN(TSTRING)
          TSTRING(LT2:LT2) = FSTRING(LT2+1:LT2+1)
 900    CONTINUE
        IERR = .TRUE.
*
* Normal case -> TSTRING is long enough.
*
      ELSE
        TSTRING = FSTRING
      ENDIF
*
*           Center Justify the string in the space, LENT
*
      CALL CENT_JST(TSTRING, LEN(TSTRING))
*
      RETURN
      END
C----------------------------------------------------------------------C
*
      SUBROUTINE CENT_JST(STRING, SPACE)
*
      CHARACTER           STRING*(*)
      INTEGER             SPACE
C
C This subroutine will center-justify the string, STRING, in the space,
C SPACE.  If the length of STRING is such that this is not possible,
C then the subroutine will right justify STRING in the available space.
C The first and last significant characters in STRING, which are used
C in determining positions for center justification, are obtained from
C the CHEMKIN functions, CKFRCH and CKLSCH.
C
      INTEGER LENG, LT, FT, EXTRA_SPACE, SIG_LEN, I, SHIFT_C
      INTEGER  CKLSCH, CKFRCH
      EXTERNAL CKLSCH, CKFRCH
C
      LENG = LEN(STRING)
      IF (LENG .EQ. 0) RETURN
      FT = MAX(CKFRCH(STRING), 1)
      LT = MAX(CKLSCH(STRING), 1)
      SIG_LEN = LT - FT + 1
      EXTRA_SPACE = MAX( ((SPACE -  SIG_LEN) / 2) , 0 )
      EXTRA_SPACE = MIN( EXTRA_SPACE, (LENG - SIG_LEN) )
      SHIFT_C = EXTRA_SPACE - FT + 1
      IF (SHIFT_C .GT. 0) THEN
        DO 10 I = LT, FT, -1
          STRING((I+SHIFT_C):(I+SHIFT_C)) = STRING(I:I)
10      CONTINUE
        DO 11 I = FT, FT + SHIFT_C - 1
          STRING(I:I) = ' '
11      CONTINUE
      ELSE IF (SHIFT_C .LT. 0) THEN
        DO 20 I = FT, LT
          STRING((I+SHIFT_C):(I+SHIFT_C)) = STRING(I:I)
20      CONTINUE
        DO 21 I = LT + SHIFT_C + 1, LT
          STRING(I:I) = ' '
21      CONTINUE
      END IF
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE STKTOGEN(RA_ST, RB_ST, RE_ST, WT, NSURF, SDTOT,
     1                    TTABLE, RA, RB, RE, WORK, IERR)
C
C   Converts Arrhenius coefficients for a surface reaction rate
C constant expressed in terms of sticking coefficient parameters,
C            STK = RA_ST * T**(RB_ST) * exp ( -RE_ST / T)
C into general rate constant Arrhenius parameters (cgs units):
C            k = RA * T**(RB) * exp ( -RE / T)
C Requires a table of temperatures, over which the parameters
C are fit with a least squares algorithm (due to the velocity
C distribution correction factor).
C
C The formula for conversion of the sticking coefficient to the
C reaction rate uses the velocity distribution correction factor.
C Because this factor becomes singular when the sticking coefficient
C is 2 and because mechanisms may often inadvertently dictate
C sticking coefficients greater than two, I have implemented
C a change in the equations.  In this implementation, the
C sticking coefficient can never be greater than 1.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Parameters
      INTEGER  NSURF
      LOGICAL IERR
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    RA_ST, RB_ST, RE_ST, WT, SDTOT, TTABLE(TTABLE_NUM),
     1    RA, RB, RE, WORK(*)
C Local Variables
      INTEGER ITEMP
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    ALOGT, STK
C Externals
      EXTERNAL FITKT
      IERR = .FALSE.
      DO 10 ITEMP = 1, TTABLE_NUM
        ALOGT = LOG(TTABLE(ITEMP))
        STK   = RA_ST * EXP( RB_ST*ALOGT - RE_ST/TTABLE(ITEMP) )
        IF (STK .GT. 1.0) THEN
          STK = 1.0
          IERR = .TRUE.
        END IF
        IF (NSURF .EQ. 0) THEN
          WORK(ITEMP) = STK * 3637.6011 * SQRT(TTABLE(ITEMP)/WT)
     1              / (1.0 - 0.5*STK)
        ELSE
          WORK(ITEMP) = STK * 3637.6011 * SQRT(TTABLE(ITEMP)/WT)
     1              / ((1.0 - 0.5*STK) * SDTOT**NSURF)
        END IF
10    CONTINUE
      CALL FITKT (TTABLE_NUM, TTABLE, WORK, RA, RB, RE)
      RE = RE / RUC
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE STKCALCA(RA_ST, RB_ST, RE_ST, RP_SUR, WT, NSURF, SDTOT,
     $                    T, RATE_CNS, A_FACTOR, EACT)
C
C   Calculates a surface reaction rate constant, pre-exponential factor,
C and activation energy given  a sticking coefficient expressed
C in terms of a three parameter Arrhenius-like expression
C            STK = RA_ST * T**(RB_ST) * exp ( -RE_ST / T)
C The surface reaction rate constant has the appropriate cgs units.
C The final rate constant is also multiplied by a perturbation factor,
C           RF = RF * RP_SUR(I)
C
C The formula for conversion of the sticking coefficient to the
C reaction rate uses the velocity distribution correction factor.
C Because this factor becomes singular when the sticking coefficient
C is 2 and because mechanisms may often inadvertently dictate
C sticking coefficients greater than two, I have implemented
C a change in the equations.  In this implementation, the
C sticking coefficient can never be greater than 1.
C EACT is returned with units of kcal/mole.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Parameters
      INTEGER  NSURF
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    RA_ST, RB_ST, RE_ST, RP_SUR, WT, SDTOT, T, RATE_CNS,
     2    A_FACTOR, EACT
C Local Variables
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1             DELTA_T, RATE_CNS_DELTA, T2
      SAVE DELTA_T
C Externals
      EXTERNAL STKCALC
      DATA DELTA_T/0.01/
      CALL STKCALC(RA_ST, RB_ST, RE_ST, RP_SUR, WT, NSURF, SDTOT,
     1                    T, RATE_CNS)
C
      T2 = T + DELTA_T
      CALL STKCALC(RA_ST, RB_ST, RE_ST, RP_SUR, WT, NSURF, SDTOT,
     1                    T2, RATE_CNS_DELTA)
C
      CALL CALC_ARR(RATE_CNS, T, RATE_CNS_DELTA, T2, A_FACTOR, EACT)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE STKCALC(RA_ST, RB_ST, RE_ST, RP_SUR, WT, NSURF, SDTOT,
     1                    T, RATE_CNS)
C
C   Calculates a general surface reaction rate constant,
C given a sticking coefficient expressed
C in terms of a three parameter Arrhenius-like expression
C            STK = RA_ST * T**(RB_ST) * exp ( -RE_ST / T)
C The surface reaction rate constant has the appropriate cgs units.
C The final rate constant is also multiplied by a perturbation factor,
C           RF = RF * RP_SUR(I)
C
C The formula for conversion of the sticking coefficient to the
C reaction rate uses the velocity distribution correction factor.
C Because this factor becomes singular when the sticking coefficient
C is 2 and because mechanisms may often inadvertently dictate
C sticking coefficients greater than two, I have implemented
C a change in the equations.  In this implementation, the
C sticking coefficient can never be greater than 1.
C EACT is returned with units of kcal/mole.
C
C Dummy Parameters
      INTEGER  NSURF
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    RA_ST, RB_ST, RE_ST, RP_SUR, WT, SDTOT, T, RATE_CNS
C Local Variables
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1             STK
      STK   = RA_ST * EXP( RB_ST*LOG(T) - RE_ST/T )
C*****B-2) double precision
      STK = MIN (STK, 1.0D0)
      IF (NSURF .EQ. 0) THEN
        RATE_CNS = RP_SUR * STK * 3637.6011D0 * SQRT(T/WT)
     1               / (1.0D0 - 0.5D0*STK)
      ELSE
        RATE_CNS = RP_SUR * STK * 3637.6011D0 * SQRT(T/WT)
     1            / ((1.0D0 - 0.5D0*STK) * SDTOT**NSURF)
      END IF
C*****END B-2) double precision
C*****B-1) single precision
C      STK = MIN (STK, 1.0)
C      IF (NSURF .EQ. 0) THEN
C        RATE_CNS = RP_SUR * STK * 3637.6011 * SQRT(T/WT)
C     1               / (1.0 - 0.5*STK)
C      ELSE
C        RATE_CNS = RP_SUR * STK * 3637.6011 * SQRT(T/WT)
C     1            / ((1.0 - 0.5*STK) * SDTOT**NSURF)
C      END IF
C*****END B-1) single precision
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE GENTOSTK(RA_ST, RB_ST, RE_ST, WT, NSURF, SDTOT,
     $                    TTABLE, RA, RB, RE, WORK, IERR)
C
C   Converts a general rate constant expression (appropriate cgs units)
C            k = RA * T**(RB) * exp ( -RE / T)
C into sticking coefficient expression (unitless).
C            STK = RA_ST * T**(RB_ST) * exp ( -RE_ST / T)
C Requires a table of temperatures, over which the parameters
C are fit with a least squares algorithm.
C
C The formula for conversion of the rate constant to
C a sticking coefficient
C uses the velocity distribution correction factor.
C This factor becomes singular when the sticking coefficient
C is 2. Moreover, mechanisms may often inadvertently dictate
C sticking coefficients greater than two.  From this subroutine,
C the sticking coefficient will never be greater than 2
C due to this singularity.
C However, if the sticking coefficient determined by this
C subroutine is greater than 1.0, an error flag, IERR, is set.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Parameters
      INTEGER  NSURF
      LOGICAL IERR
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     $    RA_ST, RB_ST, RE_ST, WT, SDTOT,
     $    TTABLE(TTABLE_NUM), RA, RB, RE, WORK(*)
C Local Variables
      INTEGER ITEMP
      LOGICAL KERR
C Externals
      EXTERNAL FITKT, CALCSTK
      IERR = .FALSE.
      DO 10 ITEMP = 1, TTABLE_NUM
        CALL CALCSTK(WORK(ITEMP), WT, NSURF, SDTOT,
     1               TTABLE(ITEMP), RA, RB, RE, KERR)
        IERR = IERR .OR. KERR
10    CONTINUE
      CALL FITKT (TTABLE_NUM, TTABLE, WORK, RA_ST, RB_ST, RE_ST)
      RE_ST = RE_ST / RUC
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE CALCSTKA(STK, WT, NSURF, SDTOT, T, RA, RB, RE, IERR,
     1                    A_FACTOR, EACT)
C
C   Calculates a sticking coefficient given a rate constant given
C in terms of a general rate constant Arrhenius expression
C (appropriate cgs units):
C            k = RA * T**(RB) * exp ( -RE / T)
C Also calculates the A_factor and activation energy of the sticking
C coefficient.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Parameters
      INTEGER  NSURF
      LOGICAL IERR
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    STK, T, RA, RB, RE, WT, SDTOT, A_FACTOR, EACT
C Local Variables
      LOGICAL KERR
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    T2, DELTA_T, STK_DELTA
C Externals
      EXTERNAL CALCSTK
      DATA DELTA_T/0.01/
C
      CALL CALCSTK(STK, WT, NSURF, SDTOT, T, RA, RB, RE, IERR)
      T2 = T + DELTA_T
      CALL CALCSTK(STK_DELTA, WT, NSURF, SDTOT, T2, RA, RB, RE, KERR)
      IERR = IERR .OR. KERR
      CALL CALC_ARR(STK, T, STK_DELTA, T2, A_FACTOR, EACT)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE CALCSTK(STK, WT, NSURF, SDTOT, T, RA, RB, RE, IERR)
C
C   Calculates a sticking coefficient given a rate constant given
C in terms of a general rate constant Arrhenius expression
C (appropriate cgs units):
C            k = RA * T**(RB) * exp ( -RE / T)
C Also calculates the A_factor and activation energy of the sticking
C coefficient.
C
C The formula for conversion of the rate constant to
C a sticking coefficient
C uses the velocity distribution correction factor.
C This factor becomes singular when the sticking coefficient
C is 2. Moreover, mechanisms may often inadvertently dictate
C sticking coefficients greater than two.  From this subroutine,
C the sticking coefficient will never be greater than 2
C due to this singularity.
C However, if the sticking coefficient determined by this
C subroutine is greater than 1.0, an error flag, IERR, is set.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Parameters
      INTEGER  NSURF
      LOGICAL IERR
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    STK, T, RA, RB, RE, WT, SDTOT
C Local Variables
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    ALOGT, AK, FAC
C Externals
C*****B-2) double precision
      DOUBLE PRECISION D1MACH
      EXTERNAL         D1MACH
C*****END B-2) double precision
C*****B-1) single precision
C      REAL            R1MACH
C      EXTERNAL        R1MACH
C*****END B-1) single precision
      IERR = .FALSE.
      ALOGT = LOG(T)
      AK   = RA * EXP( RB*ALOGT - RE/T )
      IF (NSURF .EQ. 0) THEN
         FAC = 3637.6011 * SQRT(T/WT)
      ELSE
         FAC = 3637.6011 * SQRT(T/WT) / (SDTOT**NSURF)
      END IF
      STK =  MIN( AK/(0.5*AK + FAC) ,
C*****B-2) double precision
     1       (2.0D0 - 2.0D0 * D1MACH(3)) )
C*****END B-2) double precision
C*****B-1) single precision
C     1       (2.0 - 2.0 * R1MACH(3)) )
C*****END B-1) single precision
      IF (STK .GT. 1.001) IERR = .TRUE.
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE REVSTK(SKWRK, ISKWRK, I, RA_REV_ST, RB_REV_ST,
     1                  RE_REV_ST, WT_FOR, NSURF_FOR, WT_REV,
     1                  NSURF_REV, SDTOT, TTABLE, ISTFL_SUR, RA,
     1                  RB, RE, ACT, SDEN, SEQKC, WORK, IERR)
C
C   Converts a forward rate constant (appropriate cgs units)
C expressed either in terms of a general rate expression
C            k   = RA * T**(RB) * exp ( -RE / T)
C or a sticking coefficient expression
C            STK = RA * T**(RB) * exp ( -RE / T)
C into sticking coefficient expression for the reverse reaction.
C            STK = RA_ST * T**(RB_ST) * exp ( -RE_ST / T)
C Requires a table of temperatures, over which the parameters
C are fit with a least squares algorithm.
C
C The formula for conversion of the rate constant to
C a sticking coefficient
C uses the velocity distribution correction factor.
C This factor becomes singular when the sticking coefficient
C is 2. Moreover, mechanisms may often inadvertently dictate
C sticking coefficients greater than two.  From this subroutine,
C the sticking coefficient will never be greater than 2
C due to this singularity.
C However, if the sticking coefficient determined by this
C subroutine is greater than 1.0, an error flag, IERR, is set.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Parameters
      INTEGER  ISKWRK(*), I, ISTFL_SUR, NSURF_FOR, NSURF_REV
      LOGICAL IERR
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    SKWRK(*), RA_REV_ST, RB_REV_ST, RE_REV_ST, WT_FOR, WT_REV,
     2    SDTOT, TTABLE(TTABLE_NUM), RA, RB, RE, ACT(*),
     3    SDEN(*), SEQKC(*), WORK(*)
C
C Local Variables
      INTEGER ITEMP
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    ALOGT, AK, FAC, AK_REV, RP
C Externals
      EXTERNAL FITKT
      IERR = .FALSE.
      DO 10 ITEMP = 1, TTABLE_NUM
C     Calc the forward rate constant
        IF (ISTFL_SUR .GT. 0) THEN
          RP = 1.0
          CALL STKCALC(RA, RB, RE, RP, WT_FOR, NSURF_FOR, SDTOT,
     1                 TTABLE(ITEMP), AK)
        ELSE
          ALOGT = LOG(TTABLE(ITEMP))
          AK    = RA * EXP( RB*ALOGT - RE/TTABLE(ITEMP) )
        END IF
C     Calc the reverse rate constant from a call to SURFACE CHEMKIN
C       to find the equilibrium constant.
        CALL SKEQ (PATM, TTABLE(ITEMP), ACT, SDEN, ISKWRK, SKWRK, SEQKC)
        AK_REV = AK / MAX(SEQKC(I),SMALL)
C     Calc the reverse sticking coefficient
        IF (NSURF_REV .EQ. 0) THEN
          FAC = 3637.6011 * SQRT(TTABLE(ITEMP)/WT_REV)
        ELSE
          FAC = 3637.6011 * SQRT(TTABLE(ITEMP)/WT_REV)
     1                    / (SDTOT**NSURF_REV)
        END IF
        WORK(ITEMP) =  AK_REV/(0.5*AK_REV + FAC)
        IF (WORK(ITEMP) .GT. 1.001) IERR = .TRUE.
10    CONTINUE
      CALL FITKT (TTABLE_NUM, TTABLE, WORK, RA_REV_ST, RB_REV_ST,
     1            RE_REV_ST)
      RE_REV_ST = RE_REV_ST / RUC
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE CALCRSTA(SKWRK, ISKWRK, I, STK_REV, WT_FOR,
     1                    NSURF_FOR, WT_REV, NSURF_REV, SDTOT,
     1                    T, ISTFL_SUR, RA, RB, RE,
     1                    ACT, SDEN, SEQKC, IERR, A_FACTOR, EACT)
C
C    Calculates a reverse sticking coefficient given
C a forward rate constant (appropriate cgs units)
C expressed either in terms of a general rate expression
C            k   = RA * T**(RB) * exp ( -RE / T)
C or a sticking coefficient expression
C            STK = RA * T**(RB) * exp ( -RE / T)
C This subroutine also calculates the A factor and activation energy
C of the reverse sticking coefficient.
C
C The formula for conversion of the rate constant to
C a sticking coefficient
C uses the velocity distribution correction factor.
C This factor becomes singular when the sticking coefficient
C is 2. Moreover, mechanisms may often inadvertently dictate
C sticking coefficients greater than two.  From this subroutine,
C the sticking coefficient will never be greater than 2
C due to this singularity.
C However, if the sticking coefficient determined by this
C subroutine is greater than 1.0, an error flag, IERR, is set.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Parameters
      INTEGER  ISKWRK(*), I, ISTFL_SUR, NSURF_FOR, NSURF_REV
      LOGICAL IERR
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    SKWRK(*), STK_REV, WT_FOR, WT_REV, SDTOT, T, RA, RB, RE,
     1    ACT(*), SDEN(*), SEQKC(*), A_FACTOR, EACT
C Local Variables
      LOGICAL KERR
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1             DELTA_T, STK_REV_DELTA, T2
      SAVE DELTA_T
C Externals
      EXTERNAL CALCRST
      DATA DELTA_T/0.01/
      CALL CALCRST(SKWRK, ISKWRK, I, STK_REV, WT_FOR,
     1                    NSURF_FOR, WT_REV, NSURF_REV, SDTOT,
     1                    T, ISTFL_SUR, RA, RB, RE,
     1                    ACT, SDEN, SEQKC, IERR)
C
      T2 = T + DELTA_T
      CALL CALCRST(SKWRK, ISKWRK, I, STK_REV_DELTA, WT_FOR,
     1                    NSURF_FOR, WT_REV, NSURF_REV, SDTOT,
     1                    T2, ISTFL_SUR, RA, RB, RE,
     1                    ACT, SDEN, SEQKC, KERR)
      IERR = IERR .OR. KERR
C
      CALL CALC_ARR(STK_REV, T, STK_REV_DELTA, T2, A_FACTOR, EACT)
      RETURN
      END
*
C----------------------------------------------------------------------C
*
      SUBROUTINE CALCRST(SKWRK, ISKWRK, I, STK_REV, WT_FOR, NSURF_FOR,
     1                   WT_REV, NSURF_REV, SDTOT, T, ISTFL_SUR, RA, RB,
     1                   RE, ACT, SDEN, SEQKC, IERR)
C
C    Calculates a reverse sticking coefficient given
C a forward rate constant (appropriate cgs units)
C expressed either in terms of a general rate expression
C            k   = RA * T**(RB) * exp ( -RE / T)
C or a sticking coefficient expression
C            STK = RA * T**(RB) * exp ( -RE / T)
C
C The formula for conversion of the rate constant to
C a sticking coefficient
C uses the velocity distribution correction factor.
C This factor becomes singular when the sticking coefficient
C is 2. Moreover, mechanisms may often inadvertently dictate
C sticking coefficients greater than two.  From this subroutine,
C the sticking coefficient will never be greater than 2
C due to this singularity.
C However, if the sticking coefficient determined by this
C subroutine is greater than 1.0, an error flag, IERR, is set.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Parameters
      INTEGER  ISKWRK(*), I, ISTFL_SUR, NSURF_FOR, NSURF_REV
      LOGICAL  IERR
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    SKWRK(*), STK_REV, WT_FOR, WT_REV, SDTOT, T, RA, RB, RE,
     1    ACT(*), SDEN(*), SEQKC(*)
C
C Local Variables
C*****B-2) double precision
      DOUBLE PRECISION            ALOGT, AK, FAC, AK_REV, RP
C*****END B-2) double precision
C*****B-1) single precision
C      REAL                       ALOGT, AK, FAC, AK_REV, RP
C*****END B-1) single precision
C Externals
C*****B-2) double precision
      DOUBLE PRECISION D1MACH
      EXTERNAL         D1MACH
C*****END B-2) double precision
C*****B-1) single precision
C      REAL           R1MACH
C      EXTERNAL       R1MACH
C*****END B-1) single precision
      IERR = .FALSE.
C     Calc the forward rate constant
      IF (ISTFL_SUR .GT. 0) THEN
        RP = 1.0
        CALL STKCALC(RA, RB, RE, RP, WT_FOR, NSURF_FOR, SDTOT,
     1               T, AK)
      ELSE
        ALOGT = LOG(T)
        AK   = RA * EXP( RB*ALOGT - RE/T )
       END IF
C   Calc the reverse rate constant from a call to SURFACE CHEMKIN
C     to find the equilibrium constant.
      CALL SKEQ (PATM, T, ACT, SDEN, ISKWRK, SKWRK, SEQKC)
      AK_REV = AK / MAX(SEQKC(I), SMALL)
C     Calc the reverse sticking coefficient
      IF (NSURF_REV .EQ. 0) THEN
        FAC = 3637.6011 * SQRT(T/WT_REV)
      ELSE
        FAC = 3637.6011 * SQRT(T/WT_REV) / (SDTOT**NSURF_REV)
      END IF
C*****B-1) single precision
C      STK_REV =  MIN(AK_REV/(0.5*AK_REV + FAC)
C     1             , (2.0 - 2.0*R1MACH(3)))
C*****END B-1) single precision
C*****B-2) double precision
      STK_REV =  MIN(AK_REV/(0.5D0*AK_REV + FAC)
     1              , (2.0D0 - 2.0D0*D1MACH(3)))
C*****END B-2) double precision
      IF (STK_REV .GT. 1.001) IERR = .TRUE.
      RETURN
      END
*
C----------------------------------------------------------------------C
*
      SUBROUTINE REVGEN(SKWRK, ISKWRK, I, RA_REV, RB_REV, RE_REV,
     1                  WT_FOR, NSURF_FOR, SDTOT, TTABLE,
     1                  ISTFL_SUR, RA, RB, RE, ACT, SDEN, SEQKC,
     1                  WORK)
C
C   Converts a forward rate constant (appropriate cgs units)
C expressed either in terms of a general rate expression
C            k     = RA * T**(RB) * exp ( -RE / T)
C or in terms of a sticking coefficient
C            gamma =  RA * T**(RB) * exp ( -RE / T)
C into a reverse rate constant with the general rate expression
C functional form.  Requires a table of temperatures, over which the
C parameters are fit with a least squares algorithm.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Parameters
      INTEGER  ISKWRK(*), I, ISTFL_SUR, NSURF_FOR
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    SKWRK(*), RA_REV, RB_REV, RE_REV, WT_FOR, SDTOT,
     1    TTABLE(TTABLE_NUM), RA, RB, RE, ACT(*), SDEN(*),
     1    SEQKC(*), WORK(*)
C
C Local Variables
      INTEGER ITEMP
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    ALOGT, AK, RP
C Externals
      EXTERNAL FITKT, STKCALC
      DO 10 ITEMP = 1, TTABLE_NUM
C     Calc the forward rate constant
        IF (ISTFL_SUR .GT. 0) THEN
          RP = 1.0
          CALL STKCALC(RA, RB, RE, RP, WT_FOR, NSURF_FOR, SDTOT,
     1                 TTABLE(ITEMP), AK)
        ELSE
          ALOGT = LOG(TTABLE(ITEMP))
          AK    = RA * EXP( RB*ALOGT - RE/TTABLE(ITEMP) )
        END IF
C     Calc the reverse rate constant from a call to SURFACE CHEMKIN
C       to find the equilibrium constant.
        CALL SKEQ (PATM, TTABLE(ITEMP), ACT, SDEN, ISKWRK, SKWRK, SEQKC)
        WORK(ITEMP) = AK / MAX(SEQKC(I), SMALL)
10    CONTINUE
      CALL FITKT (TTABLE_NUM, TTABLE, WORK, RA_REV, RB_REV,
     1            RE_REV)
      RE_REV = RE_REV / RUC
      RETURN
      END
*
C----------------------------------------------------------------------C
*
      SUBROUTINE P_REAC(LOUT, UNITS, RA, RB, RE, RP, FORM)
C
C   Prints out a rate constant expression. The following FORMS
C are used.
C   1   = forward rate constant - 4 parameter arrhenius form
C   2   = reverse rate constnat - 4 parameter arrhenius form
C   3   = forward sticking coefficient - 4 parameter arrhenius form
C   4   = reverse sticking coefficient - 4 parameter arrhenius form
C   5   = forward sticking coefficient - 4 parameter arrhenius form
C                  min with respect to 1 is printed out
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Variables
      INTEGER LOUT, FORM
      CHARACTER UNITS*(*)
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1    RA, RB, RE, RP
C Local Variables
      INTEGER LT, FT
      CHARACTER RATE*96, SUFF*6, EXP_SGN*1
C*****B-2) double precision
      DOUBLE PRECISION              RP_LOCAL
C*****END B-2) double precision
C*****B-1) single precision
C      REAL                         RP_LOCAL
C*****END B-1) single precision
C Externals
      INTEGER  CKLSCH, CKFRCH
      EXTERNAL CKLSCH, CKFRCH
      RP_LOCAL = RP
      LT = MAX(CKLSCH(UNITS),1)
      FT = MAX(CKFRCH(UNITS),1)
      IF (FORM .EQ. 1) THEN
        RATE = 'k ' // UNITS(FT:LT) // ' = '
        SUFF = ' '
      ELSE IF (FORM .EQ. 2) THEN
        RATE = 'k(rev) ' // UNITS(FT:LT) // ' = '
        SUFF = ' '
      ELSE IF (FORM .EQ. 3) THEN
        RATE = 'Sticking Coeff = '
        SUFF = ' '
      ELSE IF (FORM .EQ. 4) THEN
        RATE = 'Sticking Coeff(rev) = '
        SUFF = ' '
      ELSE IF (FORM .EQ. 5) THEN
        IF (RP .NE. 1.0) THEN
          WRITE(LOUT,8200) RP
          RP_LOCAL = 1.0
        ELSE
          RATE = 'Sticking Coeff = MIN( '
        END IF
        SUFF = ' , 1 )'
      ELSE
        print *,'P_REAC: FORM IS UNRECOGNIZED', FORM
        STOP
      END IF
      LT = CKLSCH(RATE)
C Determine the Sign of the exponent
      IF (RE .LT. 0.0) THEN
        EXP_SGN = '+'
      ELSE
        EXP_SGN = '-'
      END IF
C
C Print out the high pressure arrhenius parameters:
      IF (RB .EQ. 0.0) THEN
        IF (RP_LOCAL .EQ. 1.0) THEN
          WRITE(LOUT,8540)  RATE(1:LT), RA, EXP_SGN, ABS(RE*RKCAL), SUFF
        ELSE
          WRITE(LOUT,8541)  RATE(1:LT), RP_LOCAL, RA, EXP_SGN,
     1                      ABS(RE*RKCAL), SUFF
        END IF
      ELSE
        IF (RP_LOCAL .EQ. 1.0) THEN
          WRITE(LOUT,8542)  RATE(1:LT), RA, RB, EXP_SGN, ABS(RE*RKCAL),
     1                      SUFF
        ELSE
          WRITE(LOUT,8543)  RATE(1:LT), RP_LOCAL, RA, RB, EXP_SGN,
     1                      ABS(RE*RKCAL), SUFF
        END IF
      END IF
      RETURN
8200  FORMAT( 'Sticking Coeff = ', 1PG10.2, ' * MIN( ' )
8540  FORMAT(/30X,A,1PG12.4,' exp( ',A1,G10.4, ' kcal/mole / RT) ',A/)
8541  FORMAT(/30X, A,1PG11.3,' * ',G12.4,' exp( ',A1,G10.4,
     1           ' kcal/mole / RT) ',A/)
8542  FORMAT(/30X, A,1PG12.4,' T**(',1PG11.4,') exp( ',A1,3PG10.4,
     1           ' kcal/mole / RT) ',A/)
8543  FORMAT(/30X, A,1PG11.3,' * ',1PG12.4,
     1   ' T**(',1PG11.4,') exp( ',A1,3PG10.4,' kcal/mole / RT) ',A/)
      END
C----------------------------------------------------------------------C
      SUBROUTINE COV_DEPA(X_BATH, T, NCOVI, KCOVI, CPARI, NSCOV,
     1                    FAC, A_FACTOR, EACT)
C
C    Calculates the coverage dependence factor for a surface reaction
C rate constant.  It also calculates the A_factor and activation
C energy due to these terms.  A bath gas is assumed.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Variables
      INTEGER NCOVI, KCOVI(NCOVI), NSCOV
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1      X_BATH(*), T, CPARI(NSCOV,NCOVI), FAC, A_FACTOR, EACT
C Local Variables
      INTEGER N, K
      A_FACTOR = 1.0
      EACT = 0.0
      DO 100 N = 1, NCOVI
        K = KCOVI(N)
        A_FACTOR = A_FACTOR * 10.0**(X_BATH(K)*CPARI(1,N))
     1                      * X_BATH(K)**CPARI(2,N)
        EACT = EACT + CPARI(3,N)*X_BATH(K)
100   CONTINUE
      FAC = A_FACTOR * EXP( - EACT / T)
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE NDCALC_GAS( ND_RCF, ND_RCR, GAS_REAC_GAS, P, T, X,
     $                       ICKWRK, CKWRK, GAS_PROD_GAS)
C
C    Calculates the "non-dimensional" forward and reverse reaction rate
C constants for gas reactions at a single temperature.
C All reaction rates will have units of mole cm-3 sec-1.
C
C Input
C ------
C   P    gas pressure in cgs units
C   T    temperature (kelvin)
C   X    Vector of mole fractions
C   Patm Pressure of one atmosphere in cgs units
C   GAS_REAC_GAS(IIGAS) Number of gas phase reactants
C   GAS_PROD_GAS(IIGAS) Number of gas phase products
C------------------------------------------------------------------
*
      INCLUDE 'surftherm.h'
*
* Dummy Variables
*
      INTEGER  GAS_REAC_GAS(IIGAS), ICKWRK(*), GAS_PROD_GAS(IIGAS)
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
     $      ND_RCF(IIGAS), ND_RCR(IIGAS), P, T, X(KKGAS), CKWRK(*)
*
* Local Variables
*
      INTEGER I
C*****B-2) double precision
      DOUBLE PRECISION               CTOT
C*****END B-2) double precision
C*****B-1) single precision
C      REAL                          CTOT
C*****END B-1) single precision
*
* Externals
*
      EXTERNAL CKRCXP
*
* Call an auxiliary routine that calculates the reaction rate
* constants
*
      CALL CKRCXP(P, T, X, ICKWRK, CKWRK, ND_RCF, ND_RCR)
      CTOT =  P/ (RU * T)
      DO 10 I = 1, IIGAS
        ND_RCF(I) = ND_RCF(I) * CTOT**GAS_REAC_GAS(I)
        ND_RCR(I) = ND_RCR(I) * CTOT**GAS_PROD_GAS(I)
10    CONTINUE
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE NDCALC_SUR(
     1              K_STAR, K_STAR_NEG1, I, ISTFL_SUR, RA_SUR, RB_SUR,
     1              RE_SUR, RP_SUR, WT, SDTOT, GAS_REAC, SUR_REAC,
     2              BULK_REAC, SDEN, KFIRST, KLAST,
     3              KSTOICF_SUR, T, X, ICOV_SUR,
     4              NCOVI, KCOVI, CPARI, NSCOV, ACT, ISKWRK,
     5              SKWRK, SEQKC, GAS_PROD, SUR_PROD, BULK_PROD,
     6              KSTOIC_SUR, KOCC)
C
C    Calculates the non-dimensional forward and reverse reaction rate
C for surface reactions.  All reaction rates will have units of
C mole/cm**2*sec.
C
*
      INCLUDE 'surftherm.h'
*
C Dummy Variables
      INTEGER  I, ISTFL_SUR, GAS_REAC, SUR_REAC, BULK_REAC,
     1         KFIRST(NNPHAS), KLAST(NNPHAS),
     2         KSTOICF_SUR(IISUR,KKTOT), ICOV_SUR, NCOVI,
     3         KCOVI(*), NSCOV, ISKWRK(*), GAS_PROD, SUR_PROD,
     4         BULK_PROD, KSTOIC_SUR(IISUR,KKTOT), KOCC(KKTOT)
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1      K_STAR, K_STAR_NEG1, RA_SUR, RB_SUR, RE_SUR, RP_SUR, WT,
     2      SDTOT, SDEN(NNPHAS), T, X(KKTOT), CPARI(NSCOV,*),
     3      ACT(KKTOT), SKWRK(*), SEQKC(*)
C
C Local Variables
      INTEGER N, K, IPHASE
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     $      A_FACTOR, EACT, FAC, CTOT
C
C Externals
      INTEGER  PHASEID
      EXTERNAL HIGH_PRESS, STKCALC, COV_DEPA, PHASEID
      EXTERNAL SKEQ, SKICOV
      IF (ISTFL_SUR .EQ. 0) THEN
        CALL HIGH_PRESS(K_STAR, A_FACTOR, EACT, RA_SUR, RB_SUR,
     1                  RE_SUR, RP_SUR, T)
      ELSE
        CALL STKCALC (RA_SUR, RB_SUR, RE_SUR, RP_SUR, WT,
     1                SUR_REAC, SDTOT, T, K_STAR)
      END IF
      IF (ICOV_SUR .NE. 0) THEN
        CALL SKICOV(I, NSCOV, ISKWRK, SKWRK, NCOVI, KCOVI, CPARI)
        CALL COV_DEPA(X, T, NCOVI, KCOVI, CPARI, NSCOV,
     1                FAC, A_FACTOR, EACT)
        K_STAR = K_STAR * FAC
      END IF
      CALL SKEQ (PATM, T, ACT, SDEN, ISKWRK, SKWRK, SEQKC)
      K_STAR_NEG1 = K_STAR / MAX(SEQKC(I), SMALL)
      CTOT =  (P_BATH/760.) / (82.05 * T)
      IF (GAS_REAC .GT. 0) THEN
        DO 10 N = 1, GAS_REAC
          K_STAR = K_STAR * CTOT
10      CONTINUE
      END IF
      IF (SUR_REAC .GT. 0) THEN
        DO 20 K = KFIRST(NFSURF), KLAST(NLSURF)
          IF (KSTOICF_SUR(I,K) .NE. 0) THEN
            IPHASE = PHASEID(K, NNPHAS, KFIRST, KLAST)
            K_STAR = K_STAR
     1             * (SDEN(IPHASE)/KOCC(K))**(ABS(KSTOICF_SUR(I,K)))
          END IF
20      CONTINUE
      END IF
      IF (GAS_PROD .GT. 0) THEN
        DO 30 N = 1, GAS_PROD
          K_STAR_NEG1 = K_STAR_NEG1 * CTOT
30      CONTINUE
      END IF
      IF (SUR_PROD .GT. 0) THEN
        DO 40 K = KFIRST(NFSURF), KLAST(NLSURF)
          IF ((KSTOIC_SUR(I,K)-KSTOICF_SUR(I,K)) .NE. 0) THEN
            IPHASE = PHASEID(K, NNPHAS, KFIRST, KLAST)
            K_STAR_NEG1 = K_STAR_NEG1 *
     1        (SDEN(IPHASE)/KOCC(K))**(KSTOIC_SUR(I,K)-KSTOICF_SUR(I,K))
          END IF
40      CONTINUE
      END IF
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE SURFKEY(LCOM_EXIST, LCOM, LIN, LOUT, LCNTUE, CNTNUD,
     1  LTHRM, LTRAN, LGRXN, LGEN, LSRXN, LNDIM_GAS, LNDIM_SUR, LSTCK,
     1  LSCOV, LPFAL, LTFAL, LGTHB, LTSUM_SPECIES, LTSUM_GAS, LTSUM_SUR,
     1  X_BATH, KNAME, KFIRST, KLAST, KKPHAS,
     1  PNAME, IGASNAME, ISURNAME, LESTIM, X_EST)
C-----------------------------------------------------------------------
C    SURFKEY collects keywords concerning one output data set
C It returns to the main program after it encounters the 'END '
C keyword or an end of file condition.  If it encounters the 'CNTR'
C keyword, it instructs the main program through the LCONT flag
C to return to this program after listing the requested data.
C
C
C                           LIST OF KEYWORDS
C-----------------------------------------------------------------------
C
C
C General keywords
C----------------------
C
C END
C       Signals end of the current keyword input.  An end of file
C       condition serves the same purpose.
C
C CNTR
C       Signals that another set of outputs will be requested.
C       After the current set of outpus is processed by the program,
C       the program will return to this subroutine to get more
C       keywords
C       The NONE keyword is automatically assumed on later SURFKEY
C       calls.
C
C NONE
C       Turns all default I/O into no output.  Use this keyword if
C       you want to turn on only one feature.  Also, turns off all
C       previously specified output.
C
C ALL
C       Turns on all I/O.  Use this keyword if you want to turn on all
C       features.  Note, this is the default, usually.
C
C
C  Keywords that Control Functionality
C--------------------------------------
C
C GEN   [ ALL ]  [ NONE ]
C       Controls the printing of general information about the
C       mechanism.
C
C TSUM  [ ALL ]  [ NONE ]  [ SPECIES ]  [ GAS ]  [ SUR ]
C       Controls the printing of summary tables for the thermodynamic
C       functions at bath gas conditions.  There are three tables:
C       one for the species, one for the gas reactions, and one
C       for the surface reactions.  The last three options turn
C       on each table individually.
C
C THRM  [ ALL ]  [ GAS ]  [ SUR ] [ BULK ] [ NONE ]  [ Species_Name... ]
C       [ Species_number... ]
C        Prints out individual thermodynamics tables for species in the
C        mechanism. The default is ALL.
C
C TRAN  [ ALL ]  [ NONE ]
C       Prints out the transport data base properties. Also
C       expands the thermo table to create a table of transport
C       properties as a function of temperature. If this parameter
C       is totally turned off, the transport data base, trandat,
C       does not have to exist in the current directory.
C
C NDIM  [ ALL ]  [ GAS ]  [ SUR ]  [ NONE ]
C       Prints out the nondimensional reaction rate constant for
C       gas and surface reactions. Note this operation may use bath
C       gas quantities.
C
C GRXN  [ ALL ]  [ NONE ]  [ Gas_Reaction_Number ... ]
C       [ Gas_Reaction_Expression ... ]
C       Prints out a table of reaction rates and other pertinant info
C       for a gas phase reaction.
C
C PFAL  [ ALL ]  [ NONE ]  [ Gas_Reaction_Number ... ]
C       [ Gas_Reaction_Expression ... ]
C       Analyse the fall-off of a gas phase reaction, i.e. create a
C       table of reaction rates versus total gas pressure at a
C       constant temperature.
C
C TFAL  [ ALL ]  [ NONE ]  [ Gas_Reaction_Number ... ]
C       [ Gas_Reaction_Expression ... ]
C       Analyse the fall-off of a gas-phase reaction, i.e. create a
C       table of reaction rates versus temperature at a constant
C       presssure.
C
C GTHB  [ ALL ]  [ NONE ]  [ Gas_Reaction_Number ... ]
C       [ Gas_Reaction_Expression ... ]
C       Create an extra table of the reaction rates for reactions which
C       have third bodies.  Use the bath gas to get effective
C       reaction rates.
C
C SRXN  [ ALL ]  [ NONE ]  [ Surface_Reaction_Number ... ]
C       [ Surface_Reaction_Expression ... ]
C       Prints out a table of reaction rates and other pertinant info
C       for a surface reaction.
C
C STCK  [ ALL ]  [ NONE ]  [ Surface_Reaction_Number ... ]
C       [ Surface_Reaction_Expression ... ]
C       Analyses the forward and reverse surface reaction's sticking
C       coefficient, if applicable.
C
C SCOV  [ ALL ]  [ NONE ]  [ Surface_Reaction_Number ... ]
C       [ Surface_Reaction_Expression ... ]
C       Analyse the coverage dependence of a surface reaction, i.e.
C       create a table of effective reaction rates.
C
C
C  Keywords that Control Bath Gas Quantities
C--------------------------------------------
C
C XBTH  Species_Name|Species_Number   Value
C       Sets the bath gas composition.  The Species_Name and mole
C       fraction are required parameters.  If one species in a
C       phase has been set with the XBATH keyword, then all specified
C       mole fractions for that phase are summed and normalized so
C       that they add up to one.  If no XBATH parameters have been
C       specified for any species in a phase, then mole fractions
C       for all species in that phase are set equal to one another.
C
C PBTH  Value
C       Set the total bath gas pressure in torr.
C
C TBTH  Value
C       Set the bath gas temperature.  This temperature is used
C       whereever a single temperature is needed.
C
C
C
C  Keywords that the Composition of Tables
C--------------------------------------------
C
C TLOW  Value
C       Set the lower limit of the temperature (K) in tables where the
C       temperature is varied.
C
C THIG  Value
C       Set the upper limit of the temperature (K) in tables where the
C       temperature is varied.
C
C TDEL  Value
C       Set the total number of temperatures in tables where the
C       temperature is varied.
C
C PLOW  Value
C       Set the lower limit of the pressure (torr) in tables where the
C       gas pressure is varied.
C
C PHIG  Value
C       Set the upper limit of the pressure (torr) in tables where the
C       gas pressure is varied.
C
C PNUM  Value
C       Set the total number of pressures in tables where the
C       gas pressure is varied.
C
C CARR  Gas_Species_Name|Gas_Species_Number
C       Set the carrier gas to a species.  This quantity is used to
C       identify the speciesin claculating binary diffusion
C       coefficients for tables and for non-dimensionalizations, which
C       require a binary diffusion coefficient.
C       Default = Use the gas species with the largest mole fraction.
C
C MAJ   Gas_Species_Name|Gas_Species_Number
C       Sets the "Major Species".  This is only used to calculate an
C       effective diffusion coefficient, when non-dimensionalizing
C       the reaction rate constants.
C
C LSCL  Value
C       Sets the length scale (cm) for the calculation of gas
C       and surface Damkoeler numbers. The default is 1. cm.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*
      INCLUDE 'surftherm.h'
*
C Dummy Variables
      INTEGER LCOM, LIN, LOUT,
     1        KFIRST(NNPHAS), KLAST(NNPHAS), KKPHAS(NNPHAS)
      LOGICAL LCOM_EXIST, CNTNUD, 
     1        LCNTUE, LTHRM(KKTOT), LTRAN, LGRXN(IIGAS), LGEN,
     1        LSRXN(IISUR), LNDIM_GAS, LNDIM_SUR, LSTCK(IISUR),
     1        LSCOV(IISUR), LPFAL(IIGAS), LTFAL(IIGAS), LGTHB(IIGAS),
     1        LTSUM_SPECIES, LTSUM_GAS, LTSUM_SUR, LESTIM(NNPHAS)
      CHARACTER KNAME(KKTOT)*16, PNAME(NNPHAS)*16, IGASNAME(IIGAS)*64,
     1          ISURNAME(IISUR)*64
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1   X_BATH(KKTOT), X_EST(KKTOT)
C Local Variables
      INTEGER    NDIM
      PARAMETER (NDIM=10)
*
      INTEGER I, K, IPHASE, KNUM, NFOUND, IWORD, LT, NT, NVAL
      LOGICAL KERR, IERR, FLAG, LCARR, LMAJ
      CHARACTER LINE*128, KEYWRD*4, SUB(NDIM)*128
C*****B-2) double precision
      DOUBLE PRECISION
C*****END B-2) double precision
C*****B-1) single precision
C      REAL
C*****END B-1) single precision
     1      VALUE(3), SUM
C Externals
      LOGICAL  ID_INT
      INTEGER  PHASEID
C      EXTERNAL SET_ALL, ASG_BATH, ID_INT, PHASEID
      EXTERNAL SET_ALL, ID_INT, PHASEID
      INTEGER  CKLSCH
      CHARACTER*4 CKCHUP
      EXTERNAL CKSUBS, CKLSCH, CKCHUP
C-----------------------------------------------------------------------
      SAVE LCARR, LMAJ
      DATA LCARR/.FALSE./, LMAJ/.FALSE./
C-----------------------------------------------------------------------
*
*        Initilialize variables
*
      KERR = .FALSE.
      DO 1 IPHASE = 1, NNPHAS
        LESTIM(IPHASE) = .FALSE.
1     CONTINUE
      DO 2 K = 1, KKTOT
        X_EST(K) = 0.0
2     CONTINUE
*
*        Do not assume that there will be a continuation run
*        after this run, unless a keyword dictates it
*
      LCNTUE = .FALSE.
*
*        The defaults for options change, when we are on a
*        continuation run or not. CNTNUD indicates whether we are.
*        For continuation runs, we want to turn off all printout
*        as the default. For the first time in the routine, all
*        printout is enabled as the default
*
      IF (CNTNUD) THEN
        WRITE(LOUT,9002)
        CALL SET_ALL(LTHRM, LTRAN, LGRXN, LGEN, LSRXN, LNDIM_GAS,
     $               LNDIM_SUR, LSTCK, LSCOV, LPFAL, LTFAL, LGTHB,
     $               LTSUM_SPECIES, LTSUM_GAS, LTSUM_SUR, .FALSE.)
      ELSE
        WRITE(LOUT,9001)
        CALL SET_ALL(LTHRM, LTRAN, LGRXN, LGEN, LSRXN, LNDIM_GAS,
     $               LNDIM_SUR, LSTCK, LSCOV, LPFAL, LTFAL, LGTHB,
     $               LTSUM_SPECIES, LTSUM_GAS, LTSUM_SUR, .TRUE.)
      END IF
*
*        If the keyword file doesn't exist, print a message and go
*        directly to the post processing of the keyword information
*        using the defaults assigned above.
*
      IF(.NOT. LCOM_EXIST) THEN
       WRITE(LOUT,1001)
       GO TO 5000
      END IF
*
C-----------------------------------------------------------------------
C
C         READ NEXT INPUT LINE
C
      WRITE (LOUT, '(/A)') '           KEYWORD INPUT '
C
   10 CONTINUE
      KEYWRD = ' '
      LINE = ' '
      READ (LCOM, '(A)', END=5000,ERR=6000) LINE
      KEYWRD = CKCHUP (LINE(1:4), 4)
      LINE(1:4) = ' '
C
C Is this a keyword comment ?
C
      IF (KEYWRD(1:1) .EQ. ' ' .OR.
     1    KEYWRD(1:1) .EQ. '.' .OR.
     3    KEYWRD(1:1) .EQ. '!'       )      GO TO 10
      WRITE(LOUT,'(1X,A,A)') KEYWRD, LINE(5:)
C
C Is this a keyword comment that is to appear in the output ?
C
      IF (KEYWRD(1:1) .EQ. '/' )            GO TO 10
C
C Strip Line of inline "!" Comments
C
C------------------------MAIN KEYWORDS----------------------------------
C
C
C####---------------------------------------------------------------NONE
C
      IF (KEYWRD .EQ. 'NONE') THEN
C
        CALL SET_ALL(LTHRM, LTRAN, LGRXN, LGEN, LSRXN, LNDIM_GAS,
     1               LNDIM_SUR, LSTCK, LSCOV, LPFAL, LTFAL, LGTHB,
     2               LTSUM_SPECIES, LTSUM_GAS, LTSUM_SUR, .FALSE.)
C
C####---------------------------------------------------------------ALL
C
      ELSEIF (KEYWRD .EQ. 'ALL') THEN
C
        CALL SET_ALL(LTHRM, LTRAN, LGRXN, LGEN, LSRXN, LNDIM_GAS,
     1               LNDIM_SUR, LSTCK, LSCOV, LPFAL, LTFAL, LGTHB,
     2               LTSUM_SPECIES, LTSUM_GAS, LTSUM_SUR, .TRUE.)
C
C####---------------------------------------------------------------END
C
      ELSEIF (KEYWRD .EQ. 'END') THEN
C
         GO TO 5000
C
C####---------------------------------------------------------------CNTR
C
      ELSEIF (KEYWRD .EQ. 'CNTR') THEN
C
         LCNTUE = .TRUE.
C
C####---------------------------------------------------------------GEN
C
      ELSEIF (KEYWRD .EQ. 'GEN') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 50 IWORD = 1, NFOUND
           IF (SUB(IWORD) .EQ. 'ALL') THEN
             LGEN = .TRUE.
           ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
             LGEN = .FALSE.
           ELSE
             LT = CKLSCH(SUB(IWORD))
             WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1             ' is not ALL or NONE'
             KERR = .TRUE.
           END IF
50       CONTINUE
C
C####---------------------------------------------------------------TSUM
C
      ELSEIF (KEYWRD .EQ. 'TSUM') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 80 IWORD = 1, NFOUND
           IF (SUB(IWORD) .EQ. 'ALL') THEN
             LTSUM_SPECIES = .TRUE.
             LTSUM_GAS     = .TRUE.
             LTSUM_SUR     = .TRUE.
           ELSE IF (SUB(IWORD) .EQ. 'SPECIES') THEN
             LTSUM_SPECIES = .TRUE.
           ELSE IF (SUB(IWORD) .EQ. 'GAS') THEN
             LTSUM_GAS     = .TRUE.
           ELSE IF (SUB(IWORD) .EQ. 'SUR' .OR.
     1              SUB(IWORD) .EQ. 'SURF') THEN
             LTSUM_SUR     = .TRUE.
           ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
             LTSUM_SPECIES = .FALSE.
             LTSUM_GAS     = .FALSE.
             LTSUM_SUR     = .FALSE.
           ELSE
             LT = CKLSCH(SUB(IWORD))
             WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1             ' is not ALL, NONE, SPECIES, GAS, or SUR'
             KERR = .TRUE.
           END IF
80       CONTINUE
C
C####---------------------------------------------------------------TRAN
C
      ELSEIF (KEYWRD .EQ. 'TRAN') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 85 IWORD = 1, NFOUND
           IF (SUB(IWORD) .EQ. 'ALL') THEN
             LTRAN = .TRUE.
           ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
             LTRAN = .FALSE.
           ELSE
             LT = CKLSCH(SUB(IWORD))
             WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1             ' is not ALL or NONE'
             KERR = .TRUE.
           END IF
85       CONTINUE
C
C####---------------------------------------------------------------THRM
C
      ELSEIF (KEYWRD .EQ. 'THRM') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 100 IWORD = 1, NFOUND
           FLAG = ID_INT(SUB(IWORD))
           IF (FLAG) THEN
             CALL CKXNUM(SUB(IWORD), 1, LOUT, NVAL, VALUE, IERR)
             IF (IERR) THEN
               KERR = .TRUE.
             ELSE
               IF (NVAL .EQ. 1) THEN
                 K = INT(VALUE(1))
                 IF (K .GT. 0 .AND. K .LE. KKTOT) THEN
                   LTHRM(K) = .TRUE.
                 ELSE
                   LT = CKLSCH(SUB(IWORD))
                   WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1               ' is not a correct species number'
                   KERR = .TRUE.
                 END IF
               END IF
             END IF
           ELSE
             IF (SUB(IWORD) .EQ. 'ALL') THEN
               DO 90 K = 1, KKTOT
                 LTHRM(K) = .TRUE.
90             CONTINUE
             ELSE IF (SUB(IWORD) .EQ. 'GAS') THEN
               DO 91 K = 1, KKGAS
                 LTHRM(K) = .TRUE.
91             CONTINUE
             ELSE IF (SUB(IWORD) .EQ. 'SUR'
     1                .OR. SUB(IWORD) .EQ. 'SURF') THEN
               IF (NNSURF .GT. 0) THEN
                 DO 92 K = KFIRST(NFSURF), KLAST(NLSURF)
                   LTHRM(K) = .TRUE.
92               CONTINUE
               END IF
             ELSE IF (SUB(IWORD) .EQ. 'BULK') THEN
               IF (NNBULK .GT. 0) THEN
                 DO 93 K = KFIRST(NFBULK), KLAST(NLBULK)
                   LTHRM(K) = .TRUE.
93               CONTINUE
               END IF
             ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
               DO 94 K = 1, KKTOT
                 LTHRM(K) = .FALSE.
94             CONTINUE
             ELSE
               CALL SKCOMP(SUB(IWORD), PNAME, NNPHAS, IPHASE, NT)
               IF (NT .GT. 0) THEN
                 DO 95 K = KFIRST(IPHASE), KLAST(IPHASE)
                   LTHRM(K) = .TRUE.
95               CONTINUE
               ELSE
                 CALL SKPCMP (SUB(IWORD), KNAME, KKTOT, PNAME, NNPHAS,
     1                        KKPHAS, KNUM, NT)
                 IF (KNUM .GT. 0 .AND. KNUM .LE. KKTOT) THEN
                   LTHRM(KNUM) = .TRUE.
                 ELSE
                   LT = CKLSCH(SUB(IWORD))
                   WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                           ' is not a species or phase name'
                   KERR = .TRUE.
                 END IF
               END IF
             END IF
           END IF
100      CONTINUE
C
C####---------------------------------------------------------------GRXN
C
      ELSEIF (KEYWRD .EQ. 'GRXN') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 200 IWORD = 1, NFOUND
           IF (SUB(IWORD) .EQ. 'ALL') THEN
             DO 190 I = 1, IIGAS
               LGRXN(I) = .TRUE.
190           CONTINUE
           ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
             DO 194 I = 1, IIGAS
               LGRXN(I) = .FALSE.
194          CONTINUE
           ELSE
             FLAG = ID_INT(SUB(IWORD))
             IF (FLAG) THEN
               CALL CKXNUM(SUB(IWORD), 1, LOUT, NVAL, VALUE, IERR)
               IF (IERR) THEN
                 KERR = .TRUE.
               ELSE
                 IF (NVAL .EQ. 1) THEN
                   I = INT(VALUE(1))
                   IF (I .GT. 0 .AND. I .LE. IIGAS) THEN
                     LGRXN(I) = .TRUE.
                   ELSE
                     LT = CKLSCH(SUB(IWORD))
                     WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                 ' is not a correct gas-phase reaction number'
                     KERR = .TRUE.
                   END IF
                 END IF
               END IF
             ELSE
               CALL SKCOMP(SUB(IWORD), IGASNAME, IIGAS, I, NT)
               IF (NT .GT. 0) THEN
                 LGRXN(I) = .TRUE.
               ELSE
                 LT = CKLSCH(SUB(IWORD))
                 WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                  ' is not the name of a gas-phase reaction '
                 KERR = .TRUE.
               END IF
             END IF
           END IF
200      CONTINUE
C
C####---------------------------------------------------------------PFAL
C
      ELSEIF (KEYWRD .EQ. 'PFAL') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 250 IWORD = 1, NFOUND
           IF (SUB(IWORD) .EQ. 'ALL') THEN
             DO 240 I = 1, IIGAS
               LPFAL(I) = .TRUE.
240           CONTINUE
           ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
             DO 245 I = 1, IIGAS
               LPFAL(I) = .FALSE.
245          CONTINUE
           ELSE
             FLAG = ID_INT(SUB(IWORD))
             IF (FLAG) THEN
               CALL CKXNUM(SUB(IWORD), 1, LOUT, NVAL, VALUE, IERR)
               IF (IERR) THEN
                 KERR = .TRUE.
               ELSE
                 IF (NVAL .EQ. 1) THEN
                   I = INT(VALUE(1))
                   IF (I .GT. 0 .AND. I .LE. IIGAS) THEN
                     LPFAL(I) = .TRUE.
                   ELSE
                     LT = CKLSCH(SUB(IWORD))
                     WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                 ' is not a correct gas-phase reaction number'
                     KERR = .TRUE.
                   END IF
                 END IF
               END IF
             ELSE
               CALL SKCOMP(SUB(IWORD), IGASNAME, IIGAS, I, NT)
               IF (NT .GT. 0) THEN
                 LPFAL(I) = .TRUE.
               ELSE
                 LT = CKLSCH(SUB(IWORD))
                 WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                  ' is not the name of a gas-phase reaction '
                 KERR = .TRUE.
               END IF
             END IF
           END IF
250      CONTINUE
C
C####---------------------------------------------------------------TFAL
C
      ELSEIF (KEYWRD .EQ. 'TFAL') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 260 IWORD = 1, NFOUND
           IF (SUB(IWORD) .EQ. 'ALL') THEN
             DO 252 I = 1, IIGAS
               LTFAL(I) = .TRUE.
252           CONTINUE
           ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
             DO 254 I = 1, IIGAS
               LTFAL(I) = .FALSE.
254          CONTINUE
           ELSE
             FLAG = ID_INT(SUB(IWORD))
             IF (FLAG) THEN
               CALL CKXNUM(SUB(IWORD), 1, LOUT, NVAL, VALUE, IERR)
               IF (IERR) THEN
                 KERR = .TRUE.
               ELSE
                 IF (NVAL .EQ. 1) THEN
                   I = INT(VALUE(1))
                   IF (I .GT. 0 .AND. I .LE. IIGAS) THEN
                     LTFAL(I) = .TRUE.
                   ELSE
                     LT = CKLSCH(SUB(IWORD))
                     WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                 ' is not a correct gas-phase reaction number'
                     KERR = .TRUE.
                   END IF
                 END IF
               END IF
             ELSE
               CALL SKCOMP(SUB(IWORD), IGASNAME, IIGAS, I, NT)
               IF (NT .GT. 0) THEN
                 LTFAL(I) = .TRUE.
               ELSE
                 LT = CKLSCH(SUB(IWORD))
                 WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                  ' is not the name of a gas-phase reaction '
                 KERR = .TRUE.
               END IF
             END IF
           END IF
260      CONTINUE
C
C####---------------------------------------------------------------GTHB
C
      ELSEIF (KEYWRD .EQ. 'GTHB') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 270 IWORD = 1, NFOUND
           IF (SUB(IWORD) .EQ. 'ALL') THEN
             DO 272 I = 1, IIGAS
               LGTHB(I) = .TRUE.
272           CONTINUE
           ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
             DO 274 I = 1, IIGAS
               LGTHB(I) = .FALSE.
274          CONTINUE
           ELSE
             FLAG = ID_INT(SUB(IWORD))
             IF (FLAG) THEN
               CALL CKXNUM(SUB(IWORD), 1, LOUT, NVAL, VALUE, IERR)
               IF (IERR) THEN
                 KERR = .TRUE.
               ELSE
                 IF (NVAL .EQ. 1) THEN
                   I = INT(VALUE(1))
                   IF (I .GT. 0 .AND. I .LE. IIGAS) THEN
                     LGTHB(I) = .TRUE.
                   ELSE
                     LT = CKLSCH(SUB(IWORD))
                     WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                 ' is not a correct gas-phase reaction number'
                     KERR = .TRUE.
                   END IF
                 END IF
               END IF
             ELSE
               CALL SKCOMP(SUB(IWORD), IGASNAME, IIGAS, I, NT)
               IF (NT .GT. 0) THEN
                 LGTHB(I) = .TRUE.
               ELSE
                 LT = CKLSCH(SUB(IWORD))
                 WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                  ' is not the name of a gas-phase reaction '
                 KERR = .TRUE.
               END IF
             END IF
           END IF
270      CONTINUE
C
C####---------------------------------------------------------------SRXN
C
      ELSEIF (KEYWRD .EQ. 'SRXN') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 300 IWORD = 1, NFOUND
           IF (SUB(IWORD) .EQ. 'ALL') THEN
             DO 290 I = 1, IISUR
               LSRXN(I) = .TRUE.
290           CONTINUE
           ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
             DO 294 I = 1, IISUR
               LSRXN(I) = .FALSE.
294          CONTINUE
           ELSE
             FLAG = ID_INT(SUB(IWORD))
             IF (FLAG) THEN
               CALL CKXNUM(SUB(IWORD), 1, LOUT, NVAL, VALUE, IERR)
               IF (IERR) THEN
                 KERR = .TRUE.
               ELSE
                 IF (NVAL .EQ. 1) THEN
                   I = INT(VALUE(1))
                   IF (I .GT. 0 .AND. I .LE. IISUR) THEN
                     LSRXN(I) = .TRUE.
                   ELSE
                     LT = CKLSCH(SUB(IWORD))
                     WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                 ' is not a correct surface-phase reaction number'
                     KERR = .TRUE.
                   END IF
                 END IF
               END IF
             ELSE
               CALL SKCOMP(SUB(IWORD), ISURNAME, IISUR, I, NT)
               IF (NT .GT. 0) THEN
                 LSRXN(I) = .TRUE.
               ELSE
                 LT = CKLSCH(SUB(IWORD))
                 WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                  ' is not the name of a surface-phase reaction '
                 KERR = .TRUE.
               END IF
             END IF
           END IF
300      CONTINUE
C
C####---------------------------------------------------------------NDIM
C
      ELSEIF (KEYWRD .EQ. 'NDIM') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 400 IWORD = 1, NFOUND
           IF (SUB(IWORD) .EQ. 'ALL') THEN
             LNDIM_GAS = .TRUE.
             LNDIM_SUR = .TRUE.
           ELSE IF (SUB(IWORD) .EQ. 'GAS') THEN
             LNDIM_GAS = .TRUE.
           ELSE IF (SUB(IWORD) .EQ. 'SUR' .OR.
     1              SUB(IWORD) .EQ. 'SURF') THEN
             LNDIM_SUR = .TRUE.
           ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
             LNDIM_GAS = .FALSE.
             LNDIM_SUR = .FALSE.
           ELSE
             LT = CKLSCH(SUB(IWORD))
             WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1             ' is not ALL GAS SUR or NONE'
             KERR = .TRUE.
           END IF
400      CONTINUE
C
C####---------------------------------------------------------------STCK
C
      ELSEIF (KEYWRD .EQ. 'STCK') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 500 IWORD = 1, NFOUND
           IF (SUB(IWORD) .EQ. 'ALL') THEN
             DO 490 I = 1, IISUR
               LSTCK(I) = .TRUE.
490           CONTINUE
           ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
             DO 494 I = 1, IISUR
               LSTCK(I) = .FALSE.
494          CONTINUE
           ELSE
             FLAG = ID_INT(SUB(IWORD))
             IF (FLAG) THEN
               CALL CKXNUM(SUB(IWORD), 1, LOUT, NVAL, VALUE, IERR)
               IF (IERR) THEN
                 KERR = .TRUE.
               ELSE
                 IF (NVAL .EQ. 1) THEN
                   I = INT(VALUE(1))
                   IF (I .GT. 0 .AND. I .LE. IISUR) THEN
                     LSTCK(I) = .TRUE.
                   ELSE
                     LT = CKLSCH(SUB(IWORD))
                     WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                 ' is not a correct surface-phase reaction number'
                     KERR = .TRUE.
                   END IF
                 END IF
               END IF
             ELSE
               CALL SKCOMP(SUB(IWORD), ISURNAME, IISUR, I, NT)
               IF (NT .GT. 0) THEN
                 LSTCK(I) = .TRUE.
               ELSE
                 LT = CKLSCH(SUB(IWORD))
                 WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                  ' is not the name of a surface-phase reaction '
                 KERR = .TRUE.
               END IF
             END IF
           END IF
500      CONTINUE
C
C####---------------------------------------------------------------SCOV
C
      ELSEIF (KEYWRD .EQ. 'SCOV') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .LE. 0) THEN
           SUB(1) = 'ALL'
           NFOUND = 1
         END IF
         DO 600 IWORD = 1, NFOUND
           IF (SUB(IWORD) .EQ. 'ALL') THEN
             DO 590 I = 1, IISUR
               LSCOV(I) = .TRUE.
590           CONTINUE
           ELSE IF (SUB(IWORD) .EQ. 'NONE') THEN
             DO 594 I = 1, IISUR
               LSCOV(I) = .FALSE.
594          CONTINUE
           ELSE
             FLAG = ID_INT(SUB(IWORD))
             IF (FLAG) THEN
               CALL CKXNUM(SUB(IWORD), 1, LOUT, NVAL, VALUE, IERR)
               IF (IERR) THEN
                 KERR = .TRUE.
               ELSE
                 IF (NVAL .EQ. 1) THEN
                   I = INT(VALUE(1))
                   IF (I .GT. 0 .AND. I .LE. IISUR) THEN
                     LSCOV(I) = .TRUE.
                   ELSE
                     LT = CKLSCH(SUB(IWORD))
                     WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                 ' is not a correct surface-phase reaction number'
                     KERR = .TRUE.
                   END IF
                 END IF
               END IF
             ELSE
               CALL SKCOMP(SUB(IWORD), ISURNAME, IISUR, I, NT)
               IF (NT .GT. 0) THEN
                 LSCOV(I) = .TRUE.
               ELSE
                 LT = CKLSCH(SUB(IWORD))
                 WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(IWORD)(1:LT),
     1                  ' is not the name of a surface-phase reaction '
                 KERR = .TRUE.
               END IF
             END IF
           END IF
600      CONTINUE
C
C####---------------------------------------------------------------XBTH
C
      ELSEIF (KEYWRD .EQ. 'XBTH') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 2) THEN
           WRITE(LOUT,*)'RDKEY ERROR: XBTH needs two arguments'
           KERR = .TRUE.
         END IF
         FLAG = ID_INT(SUB(1))
         IF (FLAG) THEN
           CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
           IF (IERR) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(2)(1:LT),
     1                    ' can not be interpreted as a species number'
             KERR = .TRUE.
           ELSE
             K = INT(VALUE(1))
             IF (K .GT. 0 .AND. K .LE. KKTOT) THEN
               CALL CKXNUM(SUB(2), 1, LOUT, NVAL, VALUE, IERR)
               IF (IERR) THEN
                 WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(2)(1:LT),
     1                   ' can not be interpreted as a mole fraction'
                 KERR = .TRUE.
               ELSE
                 IF (VALUE(1) .LT. 0.0 .OR. VALUE(1) .GT. 1.0) THEN
                   WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(2)(1:LT),
     1                     ' is greater than 1 or less than 0'
                   KERR = .TRUE.
                 ELSE
                   IPHASE = PHASEID (K, NNPHAS, KFIRST, KLAST)
                   LESTIM(IPHASE) = .TRUE.
                   X_EST(K) = VALUE(1)
                 END IF
               END IF
             ELSE
               LT = CKLSCH(SUB(1))
               WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1           ' is not a correct species number'
               KERR = .TRUE.
             END IF
           END IF
         ELSE
           CALL SKSNUM (LINE, 1, LOUT, KNAME, KKTOT, PNAME, NNPHAS,
     1                  KKPHAS, K, NT, NVAL, VALUE, IERR)
           IF (IERR) THEN
             WRITE(LOUT,'(A)') 'RDKEY ERROR: ',
     1          ' ERROR READING DATA FOR KEYWORD '//KEYWRD//':'//LINE
             KERR = .TRUE.
           ELSE
             IF (K .LE. 0 .OR. K .GT. KKTOT) THEN
               LT = CKLSCH(SUB(1))
               WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1           ' is not a correct species name'
               KERR = .TRUE.
             ELSE
               IF (VALUE(1) .LT. 0.0 .OR. VALUE(1) .GT. 1.0) THEN
                 WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(2)(1:LT),
     1                   ' is greater than 1 or less than 0'
                 KERR = .TRUE.
               ELSE
                 IPHASE = PHASEID (K, NNPHAS, KFIRST, KLAST)
                 LESTIM(IPHASE) = .TRUE.
                 X_EST(K) = VALUE(1)
               END IF
             END IF
           END IF
         END IF
C
C####---------------------------------------------------------------PBTH
C
      ELSEIF (KEYWRD .EQ. 'PBTH') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 1) THEN
           WRITE(LOUT,*)
     1       'RDKEY ERROR: PBTH needs one and only one argument'
           KERR = .TRUE.
         END IF
         CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
         IF (IERR .OR. NVAL .NE. 1) THEN
           LT = CKLSCH(SUB(1))
           WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                  ' can not be interpreted as a single number'
           KERR = .TRUE.
         ELSE
           P_BATH = VALUE(1)
           IF (P_BATH .LE. 0.0) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: P_BATH must be ',
     1                    'greater than zero.'
             KERR = .TRUE.
           END IF
         END IF
C
C####---------------------------------------------------------------TBTH
C
      ELSEIF (KEYWRD .EQ. 'TBTH') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 1) THEN
           WRITE(LOUT,*)
     1       'RDKEY ERROR: TBTH needs one and only one argument'
           KERR = .TRUE.
         END IF
         CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
         IF (IERR .OR. NVAL .NE. 1) THEN
           LT = CKLSCH(SUB(1))
           WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                  ' can not be interpreted as a single number'
           KERR = .TRUE.
         ELSE
           T_BATH = VALUE(1)
           IF (T_BATH .LE. 0.0) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: T_BATH must be ',
     1                    'greater than zero.'
             KERR = .TRUE.
           END IF
         END IF
C
C####---------------------------------------------------------------PLOW
C
      ELSEIF (KEYWRD .EQ. 'PLOW') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 1) THEN
           WRITE(LOUT,*)
     1       'RDKEY ERROR: PLOW needs one and only one argument'
           KERR = .TRUE.
         END IF
         CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
         IF (IERR .OR. NVAL .NE. 1) THEN
           LT = CKLSCH(SUB(1))
           WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                  ' can not be interpreted as a single number'
           KERR = .TRUE.
         ELSE
           PTABLE_MIN = VALUE(1)
           IF (PTABLE_MIN .LE. 0.0) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: PTABLE_MIN must be ',
     1                    'greater than zero.'
             KERR = .TRUE.
           END IF
         END IF
C
C####---------------------------------------------------------------PHIG
C
      ELSEIF (KEYWRD .EQ. 'PHIG') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 1) THEN
           WRITE(LOUT,*)
     1       'RDKEY ERROR: PHIG needs one and only one argument'
           KERR = .TRUE.
         END IF
         CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
         IF (IERR .OR. NVAL .NE. 1) THEN
           LT = CKLSCH(SUB(1))
           WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                  ' can not be interpreted as a single number'
           KERR = .TRUE.
         ELSE
           PTABLE_MAX = VALUE(1)
           IF (PTABLE_MAX .LE. 0.0) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: PTABLE_MAX must be ',
     1                    'greater than zero.'
             KERR = .TRUE.
           END IF
         END IF
C
C####---------------------------------------------------------------PNUM
C
      ELSEIF (KEYWRD .EQ. 'PNUM') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 1) THEN
           WRITE(LOUT,*)
     1       'RDKEY ERROR: PNUM needs one and only one argument'
           KERR = .TRUE.
         END IF
         CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
         IF (IERR .OR. NVAL .NE. 1) THEN
           LT = CKLSCH(SUB(1))
           WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                  ' can not be interpreted as a single number'
           KERR = .TRUE.
         ELSE
           PTABLE_NUM = INT(VALUE(1))
           IF (PTABLE_NUM .LE. 0) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: PTABLE_NUM must be ',
     1                    'greater than zero.'
             KERR = .TRUE.
           END IF
         END IF
C
C####---------------------------------------------------------------TLOW
C
      ELSEIF (KEYWRD .EQ. 'TLOW') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 1) THEN
           WRITE(LOUT,*)
     1       'RDKEY ERROR: TLOW needs one and only one argument'
           KERR = .TRUE.
         END IF
         CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
         IF (IERR .OR. NVAL .NE. 1) THEN
           LT = CKLSCH(SUB(1))
           WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                  ' can not be interpreted as a single number'
           KERR = .TRUE.
         ELSE
           TTABLE_MIN = VALUE(1)
           IF (TTABLE_MIN .LE. 0.0) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: TTABLE_MIN must be ',
     1                    'greater than zero.'
             KERR = .TRUE.
           END IF
         END IF
C
C####---------------------------------------------------------------THIG
C
      ELSEIF (KEYWRD .EQ. 'THIG') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 1) THEN
           WRITE(LOUT,*)
     1       'RDKEY ERROR: THIG needs one and only one argument'
           KERR = .TRUE.
         END IF
         CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
         IF (IERR .OR. NVAL .NE. 1) THEN
           LT = CKLSCH(SUB(1))
           WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                  ' can not be interpreted as a single number'
           KERR = .TRUE.
         ELSE
           TTABLE_MAX = VALUE(1)
           IF (TTABLE_MAX .LE. 0.0) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: TTABLE_MAX must be ',
     1                    'greater than zero.'
             KERR = .TRUE.
           END IF
         END IF
C
C####---------------------------------------------------------------TDEL
C
      ELSEIF (KEYWRD .EQ. 'TDEL') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 1) THEN
           WRITE(LOUT,*)
     1       'RDKEY ERROR: TDEL needs one and only one argument'
           KERR = .TRUE.
         END IF
         CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
         IF (IERR .OR. NVAL .NE. 1) THEN
           LT = CKLSCH(SUB(1))
           WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                  ' can not be interpreted as a single number'
           KERR = .TRUE.
         ELSE
           TTABLE_DELTAT = VALUE(1)
           IF (TTABLE_DELTAT .LE. 0.0) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: TTABLE_DELTAT must be ',
     1                    'greater than zero.'
             KERR = .TRUE.
           END IF
         END IF
C
C####---------------------------------------------------------------CARR
C
      ELSEIF (KEYWRD .EQ. 'CARR') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 1) THEN
           WRITE(LOUT,*)
     1      'RDKEY ERROR: CARR keyword takes one and only one argument'
         END IF
         FLAG = ID_INT(SUB(1))
         IF (FLAG) THEN
           CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
           IF (IERR .OR. NVAL .NE. 1) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1         ' is not a correct single gas species number'
             KERR = .TRUE.
           ELSE
             K = INT(VALUE(1))
             IF (K .GT. 0 .AND. K .LE. KKGAS) THEN
               ID_CARR = K
               LCARR = .TRUE.
             ELSE
               LT = CKLSCH(SUB(1))
               WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                       ' is not a correct gas species number'
               KERR = .TRUE.
             END IF
           END IF
         ELSE
           CALL SKPCMP (SUB(1), KNAME, KKGAS, PNAME, 1,
     1                  KKPHAS, KNUM, NT)
           IF (KNUM .GT. 0 .AND. KNUM .LE. KKGAS) THEN
             ID_CARR = KNUM
             LCARR = .TRUE.
           ELSE
             LT = CKLSCH(SUB(1))
             WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                     ' is not a gas species name'
             KERR = .TRUE.
           END IF
         END IF
C
C####---------------------------------------------------------------MAJ
C
      ELSEIF (KEYWRD .EQ. 'MAJ') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 1) THEN
           WRITE(LOUT,*)
     1      'RDKEY ERROR: CARR keyword takes one and only one argument'
         END IF
         FLAG = ID_INT(SUB(1))
         IF (FLAG) THEN
           CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
           IF (IERR .OR. NVAL .NE. 1) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1         ' is not a correct single gas species number'
             KERR = .TRUE.
           ELSE
             K = INT(VALUE(1))
             IF (K .GT. 0 .AND. K .LE. KKGAS) THEN
               ID_MAJOR = K
               LMAJ = .TRUE.
             ELSE
               LT = CKLSCH(SUB(1))
               WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                       ' is not a correct gas species number'
               KERR = .TRUE.
             END IF
           END IF
         ELSE
           CALL SKPCMP (SUB(1), KNAME, KKGAS, PNAME, 1,
     1                  KKPHAS, KNUM, NT)
           IF (KNUM .GT. 0 .AND. KNUM .LE. KKGAS) THEN
             ID_MAJOR = KNUM
             LMAJ = .TRUE.
           ELSE
             LT = CKLSCH(SUB(1))
             WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                     ' is not a gas species name'
             KERR = .TRUE.
           END IF
         END IF
C
C####---------------------------------------------------------------LSCL
C
      ELSEIF (KEYWRD .EQ. 'LSCL') THEN
C
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         KERR = KERR .OR. IERR
         IF (NFOUND .NE. 1) THEN
           WRITE(LOUT,*)
     1       'RDKEY ERROR: LSCL needs one and only one argument'
           KERR = .TRUE.
         END IF
         CALL CKXNUM(SUB(1), 1, LOUT, NVAL, VALUE, IERR)
         IF (IERR .OR. NVAL .NE. 1) THEN
           LT = CKLSCH(SUB(1))
           WRITE(LOUT,*) 'RDKEY ERROR: ',SUB(1)(1:LT),
     1                  ' can not be interpreted as a single number'
           KERR = .TRUE.
         ELSE
           LENGTH_SCALE = VALUE(1)
           IF (LENGTH_SCALE .LE. 0.0) THEN
             WRITE(LOUT,*) 'RDKEY ERROR: LENGTH_SCALE must be greater',
     1                  ' than zero.'
             KERR = .TRUE.
           END IF
         END IF
C
C####----------------------------------------------------UNKNOWN KEYWORD
C
      ELSE
C
        WRITE(LOUT,*) 'RDKEY ERROR: ',KEYWRD,' is an unknown keyword'
        KERR = .TRUE.
C
      END IF
C
      IF (.NOT. KERR) GO TO 10
C
C------------------------------------END OF KEYWORDS--------------------
C
 5000 CONTINUE
C
C
C        STOP IF ERRORS ENCOUNTERED
C
      IF (KERR) STOP
C
C-----------------------------POST PROCESSING OF KEYWORD INFORMATION----
C
      DO 5550 IPHASE = 1 , NNPHAS
*
*         If information about the bath gas was supplied for a
*         phase, Process it here. It was storred in X_EST(K).
*
        IF (LESTIM(IPHASE)) THEN
          SUM = 0.0
          DO 5400 K = KFIRST(IPHASE), KLAST(IPHASE)
             SUM = SUM + X_EST(K)
5400      CONTINUE
          IF (SUM .GT. 0.0) THEN
             DO 5500 K = KFIRST(IPHASE), KLAST(IPHASE)
                X_BATH(K) = X_EST(K) / SUM
5500         CONTINUE
             IF (ABS(SUM-1.0) .GT. 1.E-3) THEN
               WRITE (LOUT, *)
     $           'RDKEY CAUTION...XEST MOLE FRACTIONS SUM TO ', SUM,
     $           ' FOR PHASE', IPHASE,' = ',PNAME(IPHASE)
               WRITE(LOUT,*)'           They have been renormalized'
             END IF
          ELSE
            LT = CKLSCH(PNAME(IPHASE))
            WRITE(LOUT,*)
     $        'RDKEY ERROR: XBTH keywords were supplied for phase, ',
     $         PNAME(IPHASE)(1:LT),
     $        '.  However all mole fractions were zero.'
            STOP
          ENDIF
        ELSE
*
*           If this is the first time and no information was
*           provided, then assume mole fractions are equal
*
          IF ( .NOT. CNTNUD) THEN
            IF (KFIRST(IPHASE) .GT. 0) THEN
              DO 5520 K = KFIRST(IPHASE), KLAST(IPHASE)
C*****B-2) double precision
                X_BATH(K) =1.0D0/DBLE(KLAST(IPHASE)-KFIRST(IPHASE)+1)
C*****END B-2) double precision
C*****B-1) single precision
C               X_BATH(K) =1.0/FLOAT(KLAST(IPHASE)-KFIRST(IPHASE)+1)
C*****END B-1) single precision
5520          CONTINUE
            END IF
          END IF
        END IF
5550  CONTINUE
*
*        If the carrier gas has not been specified, set its default
*        to be the gas with the largest bath gas concentration
*        (or the last gas phase species, as a last resort)
*
      IF (.NOT. LCARR) THEN
        ID_CARR = KKGAS
        DO 5600 K = 1, KKGAS
          IF (X_BATH(K) .GT. X_BATH(ID_CARR)) ID_CARR = K
5600    CONTINUE
      END IF
*
*        Set the id of the major reactant in the gas to be
*        the species with the second to highest bath gas concentration
*        (or the first species, as a last resort)
*
      IF (.NOT. LMAJ) THEN
        ID_MAJOR = 1
        DO 5700 K = 1, KKGAS
          IF (X_BATH(K) .GT. X_BATH(ID_MAJOR)) THEN
            IF (ID_CARR .NE. K) ID_MAJOR = K
          END IF
5700    CONTINUE
      END IF
C
      IF (TTABLE_MAX .LT. TTABLE_MIN) THEN
        WRITE(LOUT,*)  'RDKEY ERROR: TTABLE_MAX is less than ',
     1                 'TTABLE_MIN.'
        STOP
      END IF
C
      IF (PTABLE_MAX .LT. PTABLE_MIN) THEN
        WRITE(LOUT,*)  'RDKEY ERROR: PTABLE_MAX is less than ',
     1                 'PTABLE_MIN.'
        STOP
      END IF
C
      WRITE(LOUT,9003)
      RETURN
C------------------------------------------------------------Error Block
C
C  Read Errors
C
 6000 CONTINUE
      WRITE (LOUT,*) 'RDKEY: ERROR ... error reading keyword file'
      STOP
C
C         FORMATS
C
9001  FORMAT(132('-')/10X,'Command lines read:'//)
9002  FORMAT(132('-')/10X,'Continuation command lines read:'//)
9003  FORMAT(132('-')/)
 1001 FORMAT(/' SURFKEY: Command input file, surftherm.inp, does not ',
     $       'exist.'/'          Default set of I/O will be done'/)
      END
C----------------------------------------------------------------------C
      BLOCK DATA DEFAULTS
*
      INCLUDE 'surftherm.h'
*
      DATA P_BATH /760./, T_BATH /298.15/,
     $     TTABLE_MIN/300./, TTABLE_MAX/1500./, TTABLE_DELTAT/100./,
     $     PTABLE_MIN/1.0/,  PTABLE_MAX/1000./, PTABLE_NUM/14/,
     $     LENGTH_SCALE/1.0/
*
      END
C----------------------------------------------------------------------C
      SUBROUTINE SET_ALL(LTHRM, LTRAN, LGRXN, LGEN, LSRXN, LNDIM_GAS,
     $                   LNDIM_SUR, LSTCK, LSCOV, LPFAL, LTFAL, LGTHB,
     $                   LTSUM_SPECIES, LTSUM_GAS, LTSUM_SUR, FLAG)
*
      INCLUDE 'surftherm.h'
*
* Set all I/O options either on or off
*
C Dummy Variables
      LOGICAL LTHRM(KKTOT), LTRAN, LGRXN(IIGAS), LGEN,
     $        LSRXN(IISUR), LNDIM_GAS, LNDIM_SUR, LSTCK(IISUR),
     $        LSCOV(IISUR), LPFAL(IIGAS), LTFAL(IIGAS), LGTHB(IIGAS),
     $        LTSUM_SPECIES, LTSUM_GAS, LTSUM_SUR, FLAG
*
C Local Variables
      INTEGER I
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LGEN          = FLAG
      LNDIM_GAS     = FLAG
      LNDIM_SUR     = FLAG
      LTSUM_SPECIES = FLAG
      LTSUM_GAS     = FLAG
      LTSUM_SUR     = FLAG
      LTRAN         = FLAG
      DO 10 I = 1, KKTOT
        LTHRM(I) = FLAG
10    CONTINUE
      DO 20 I = 1, IISUR
        LSRXN(I) = FLAG
        LSTCK(I) = FLAG
        LSCOV(I) = FLAG
20    CONTINUE
      DO 30 I = 1, IIGAS
        LGRXN(I) = FLAG
        LPFAL(I) = FLAG
        LTFAL(I) = FLAG
        LGTHB(I) = FLAG
30    CONTINUE
      RETURN
      END
C----------------------------------------------------------------------C
C
      LOGICAL FUNCTION NONINT(VALUE)
*
*       This small logical function will return true if its real*8
*       argument is an integer. False otherwise.
*
C*****B-2) double precision
      DOUBLE PRECISION        VALUE, TRUNC
      TRUNC      = DBLE ( INT ( VALUE ) )
C*****END B-2) double precision
C*****B-1) single precision
C      REAL                   VALUE, TRUNC
C      TRUNC      = REAL ( INT ( VALUE ) )
C*****END B-1) single precision
      IF (TRUNC .EQ. VALUE) THEN
        NONINT = .FALSE.
      ELSE
        NONINT = .TRUE.
      END IF
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      INTEGER FUNCTION CHAR_EXP (STRING_OUT, VALUE)
*
*  This integer function is meant to translate a real number into
*  a string variable. It is used for output of power law expressions
*  and is meant to use as few characters as possible. A hard limit
*  of 4 characters plus perhaps a minus sign is imposed.
*  If the length of STRING_OUT is even less, then the routine will
*  attempt to accommodate a smaller length.
*    If the routine can not do a reasonable job of translating
*  within the given amount of room, then a character string filled
*  with '#' is returned.
*
*      INPUT
*        VALUE           = real** value, usually a positive integer
*        LEN(STRING_OUT) = 1 to 5
*                            (4 or 5 is optimal)
*      OUTPUT
*        STRING_OUT   = string, no longer than 5 character or
*                       the length of string, whichever is
*                       smaller.
*
*   The subroutine returns the length of the string used.
*
*
      INTEGER   LENS, LT, LTMIN, ST
      LOGICAL   NEGYES
      CHARACTER STRING_OUT*(*), STRING1*4
C*****B-2) double precision
      DOUBLE PRECISION         EXPONENT, VALUE
C*****END B-2) double precision
C*****B-1) single precision
C      REAL                    EXPONENT, VALUE
C*****END B-1) single precision
      LOGICAL   NONINT
      EXTERNAL  NONINT
*
      LENS = LEN(STRING_OUT)
      IF (LENS .LT. 1) THEN
        CHAR_EXP = 0
        RETURN
      ENDIF
      LT = 1
*
* Handle zero and negative exponent
*
      IF (VALUE .EQ. 0.0) THEN
        STRING_OUT = '0'
        CHAR_EXP = 1
        RETURN
      ELSE IF (VALUE .LT. 0.0) THEN
        NEGYES = .TRUE.
        EXPONENT = - VALUE
      ELSE
        EXPONENT = VALUE
        NEGYES = .FALSE.
      ENDIF
*
*  First do a "best string*4" translation for the variable
*     LT = length of the string
*     LTMIN = Min length of string that would be meaningful
*             if the end of the string were whacked off
*
*  Handle the real exponent case
*
      IF (NONINT(EXPONENT)) THEN
        IF (EXPONENT .LT.  0.001) THEN
          STRING1 = '####'
          LT = 4
          LTMIN = 1
        ELSE IF (EXPONENT .LT.  1.0) THEN
          IF (NONINT(10.0 * EXPONENT)) THEN
            IF (NONINT(100.0 * EXPONENT)) THEN
              WRITE (STRING1, 3043) EXPONENT
              LT = 4
              IF (EXPONENT .GE. 0.1) THEN
                LTMIN = 2
              ELSE IF (EXPONENT .GE. 0.001) THEN
                LTMIN = 3
              ELSE
                LTMIN = 4
              ENDIF
            ELSE
              WRITE (STRING1, 3032) EXPONENT
              LT = 3
              IF (EXPONENT .GE. 0.1) THEN
                LTMIN = 2
              ELSE
                LTMIN = 3
              ENDIF
            ENDIF
          ELSE
            WRITE (STRING1,3021) EXPONENT
            LT = 2
            LTMIN = LT
          ENDIF
        ELSE IF (EXPONENT .LT. 10.0) THEN
          IF (NONINT(10.0 * EXPONENT)) THEN
            WRITE (STRING1,3042) EXPONENT
            LT = 4
          ELSE
            WRITE (STRING1,3031) EXPONENT
            LT = 3
          ENDIF
          LTMIN = 2
        ELSE IF (EXPONENT .LT. 100.0) THEN
          WRITE (STRING1,3041) EXPONENT
          LT =  4
          LTMIN = 2
        ELSE IF (EXPONENT .LT. 1000.0) THEN
          WRITE (STRING1,3040) EXPONENT
          LT =  4
          LTMIN = 3
        ELSE IF (EXPONENT .LT. 10000.0) THEN
          WRITE (STRING1,2004) INT(EXPONENT)
          LT = 4
          LTMIN = 4
        ELSE
          STRING1 = '####'
          LT = 4
          LTMIN = 1
        ENDIF
*
* Handle Integer Exponent case
*
      ELSE
        IF (EXPONENT .LT. 10.0) THEN
          WRITE (STRING1,2001) INT(EXPONENT)
          LT = 1
          LTMIN = LT
        ELSE IF (EXPONENT .LT. 100.0) THEN
          WRITE (STRING1,2002) INT(EXPONENT)
          LT = 2
          LTMIN = LT
        ELSE IF (EXPONENT .LT. 1000.0) THEN
          WRITE (STRING1,2003) INT(EXPONENT)
          LT = 3
          LTMIN = LT
        ELSE IF (EXPONENT .LT. 10000.0) THEN
          WRITE (STRING1,2004) INT(EXPONENT)
          LT =  4
          LTMIN = LT
        ELSE
          STRING1 = '####'
          LT = 4
          LTMIN = 1
        ENDIF
      ENDIF
*
*  Now try to fit the string STRING1 into the output string
*
      IF (NEGYES) THEN
        ST = 1
        IF (STRING1(1:1) .EQ. '#') THEN
          STRING_OUT = '#'
        ELSE
          STRING_OUT = '-'
        ENDIF
        LENS = LENS - 1
      ELSE
        ST = 0
        STRING_OUT = ' '
      ENDIF
      IF (LT .LE. LENS) THEN
        STRING_OUT(ST+1:ST+LT) =  STRING1(1:LT)
        CHAR_EXP = ST + LT
      ELSE IF (LTMIN .LE. LENS) THEN
        LT = LENS
        STRING_OUT(ST+1:ST+LT) =  STRING1(1:LT)
        CHAR_EXP = ST + LT
      ELSE
        LT = LENS
        STRING1 = '####'
        STRING_OUT(1:1) = '#'
        STRING_OUT(ST+1:ST+LT) =  STRING1(1:LT)
        CHAR_EXP = ST + LT
      ENDIF
      RETURN
*
2001  FORMAT(I1)
2002  FORMAT(I2)
2003  FORMAT(I3)
2004  FORMAT(I4)
*
3040  FORMAT(SS,F4.0)
3021  FORMAT(SS,F2.1)
3031  FORMAT(SS,F3.1)
3041  FORMAT(SS,F4.1)
3032  FORMAT(SS,F3.2)
3042  FORMAT(SS,F4.2)
3043  FORMAT(SS,F4.3)
*
      END
C
C----------------------------------------------------------------------C
C
      LOGICAL FUNCTION ID_INT(STRING)
C
C   ID_INT will try to identify a character string as either an integer
C or something else.  It will return TRUE if it determines that
C the string is an integer.  It will return FALSE otherwise.
C
C Dummy Variables
      CHARACTER STRING*(*)
*
C Local Variables
      INTEGER LENCH, I, K, IS, IF
      CHARACTER OK(10)*1
C
C Externals
      INTEGER  CKFRCH, CKLSCH
      EXTERNAL CKFRCH, CKLSCH
C
      DATA OK/'1','2','3','4','5','6','7','8','9','0'/
C
      ID_INT = .TRUE.
      LENCH = CKLSCH(STRING)
      IF    = CKFRCH(STRING)
      IF (LENCH .EQ. 0) THEN
        ID_INT = .FALSE.
        RETURN
      END IF
      IF (STRING(IF:IF) .EQ. '-') THEN
        IS = IF+1
        IF (LENCH .EQ. IF) THEN
          ID_INT = .FALSE.
          RETURN
        END IF
      ELSE
        IS = IF
      END IF
      DO 10 I = IS, LENCH
        DO 5 K = 1, 10
          IF (STRING(I:I) .EQ. OK(K)) GO TO 9
5       CONTINUE
        IF (I .EQ. LENCH) THEN
          IF(STRING(I:I) .EQ. '.') RETURN
        END IF
        ID_INT = .FALSE.
        RETURN
9       CONTINUE
10    CONTINUE
      RETURN
      END
C
      subroutine fitkt ( npts , t , ak , a , b , Ea )
C@(#)===================================================================
C@(#)   fitkt: Version 1.5
C@(#)===================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
c
c Given a reaction rate, ak, at npts temperatures, t, this subroutine
c will fit an expression of the form
c
c                         a T**(b) exp ( - Ea / R T )
c
c to the values.
c
c     where   Ea  is returned with units of cal/mole
c             b   is dimensionless
c             a   can have any dimension
c             T   is in Kelvin
c
C=======================================================================
C
C Change History:
C
C
C=======================================================================
C
c                      DUMMY VARIABLES
c number of points
      integer npts
C*****precision > double
c vector containing the temperatures
      double precision T (npts)
c vector containing the reaction rates
      double precision ak(npts)
c Preexponential factor returned ( has the units of ak)
      double precision a
c Beta parameter in the fit
      double precision b
c activation energy returned ( cals/mole)
      double precision Ea
C*****END precision > double
C*****precision > single
C      real T(npts), ak(npts), a, b, Ea
C*****END precision > single

c                          LOCAL VARIABLES
c work space
C*****precision > single
C      real wk(20), aa(3,3),bb(3), zero, one, maxT, minT
C*****END precision > single
C*****precision > double
      double precision wk(20), aa(3,3),bb(3), zero, one, maxT, minT
C*****END precision > double
      save zero, one
C*****slatec: precision > double
      integer iwk(3)
C*****END slatec: precision > double
C*****slatec: precision > single
C      integer iwk(3)
C*****END slatec: precision > single

c gas constant in cals/(mole*K)
C*****precision > single
C      real R
C*****END precision > single
C*****precision > double
      double precision R
C*****END precision > double
      save R

C*****precision > single
C      real slnt , sit , slnk , slnt2 , slntit , slnkt
C     1              ,  sit2 , slnkit
C      real aloga, explarge
C*****END precision > single
C*****precision > double
      double precision slnt , sit , slnk , slnt2 , slntit , slnkt
     1              ,  sit2 , slnkit
      double precision aloga, explarge
C*****END precision > double
      integer i
      logical ierr2, ierr3

c                    EXTERNALS CALLED
C*****imsl: precision > single
C      external leqt2f
C      integer ier, isig
C*****END imsl: precision > single
C*****slatec: precision > double
      external dgefs
      integer ind
C*****END slatec: precision > double
C*****slatec: precision > single
C      external sgefs
C      integer ind
C*****END slatec: precision > single

C String for the "what" command:
c
      character*40 versn
      data         versn
     1   /'@(#) fitkt.f V.1.5: sccs v1.2 02/07/94'/
C
c-----------------------------------------------------------------------
C
C*****precision > single
C      data R /1.987/, zero/0.0/, one/1.0/
C      data explarge /85./
C*****END precision > single
C*****precision > double
      data R /1.987d0/, zero/0.0d0/, one/1.0d0/
      data explarge /700.d0/
C*****END precision > double
C
c-----------------------------------------------------------------------
C
C*****printing > on
C      print *,'FITKT: number of points = ',npts
C      print *
C*****END printing > on
C
      ierr2 = .false.
      ierr3 = .false.
c collect sums for the matrix
      slnt  =zero
      sit   =zero
      slnk  =zero
      slnt2 =zero
      slntit=zero
      slnkt =zero
      sit2  =zero
      slnkit=zero
      maxT  = t(1)
      minT  = t(1)

      do 100 i=1,npts

        if (t(i) .le. 0.0) then
          print *,'fitkt: Error point ,',i, ' has a temperature ',
     1            'less than zero, ', t(i)
          stop
        end if

        if (ak(i) .le. 0.0) then
          print *,'fitkt: Error point ,',i, ' has a rate constant ',
     1            'less than zero, ', ak(i)
          stop
        end if

        maxT   = max (maxT, t(i))
        minT   = min (minT, t(i))
        slnt   = slnt   + log(t(i))
        sit    = sit    - one/(r*t(i))
        slnk   = slnk   + log(ak(i))
        slnt2  = slnt2  + (log(t(i)))**2
        slntit = slntit - log(t(i))/(r*t(i))
        slnkt  = slnkt  + log(ak(i)) * log(t(i))
        sit2   = sit2   + (one/(r*t(i)))**2
        slnkit = slnkit - log(ak(i))/(r*t(i))
100     continue

      aa(1,1) = npts
      aa(1,2) = slnt
      aa(1,3) = sit
      bb(1)   = slnk
      aa(2,1) = slnt
      aa(2,2) = slnt2
      aa(2,3) = slntit
      bb(2)   = slnkt
      aa(3,1) = sit
      aa(3,2) = slntit
      aa(3,3) = sit2
      bb(3)   = slnkit

C
C             Call a linear solver to solve the 3 by 3 system
C
C*****imsl: precision > single
C      isig=10
C      call leqt2f( aa, 1, 3, 3, bb, isig, wk, ier)
C*****END imsl: precision > single
C
C*****slatec: precision > double
      call dgefs( aa, 3, 3, bb, 1, ind, wk, iwk)
      if (ind .le. 2 .and. ind .ge. 0) then
        print *,'fitkt warning: solution to fitting problem only ',
     1          'has ', ind,' significant digits.'
      else if (ind .lt. 0) then
        print *,'fitkt error: dgefs reports error, ind = ', ind
        ierr3 = .true.
      end if
C*****END slatec: precision > double
C
C*****slatec: precision > single
C      call sgefs( aa, 3, 3, bb, 1, ind, wk, iwk)
C      if (ind .le. 2 .and. ind .ge. 0) then
C        print *,'fitkt warning: solution to fitting problem only ',
C     1          'has ', ind,' significant digits.'
C      else if (ind .lt. 0) then
C        print *,'fitkt error: sgefs reports error, ind = ', ind
C        ierr3 = .true.
C      end if
C*****END slatec: precision > single
C
      aloga = bb(1)
      if (aloga .gt. explarge) then
        if (.not. ierr3) then
          print *,'fitkt warning: A factor is too large for machine.'
          print *,'               Will try a two parameter fit.'
          ierr3 = .true.
        end if
        a  = exp (explarge)
      else
        a     = exp(aloga)
      end if
C
      if (bb(2) .gt. 0) then
        if (bb(2) .gt. (explarge/log(maxT))) then
          if (.not. ierr3) then
            print *,'fitkt warning: b is too large for machine'
            print *,'               Will try a two parameter fit.'
            ierr3 = .true.
          end if
          b = explarge/log(maxT)
        else
          b     = bb(2)
        end if
      else
        if (abs(bb(2)) .gt. (explarge/log(maxT))) then
          if (.not. ierr3) then
            print *,'fitkt warning: b is too small for machine'
            print *,'               Will try a two parameter fit.'
            ierr3 = .true.
          end if
          b = - explarge/log(maxT)
        else
          b     = bb(2)
        end if
      end if
C
      if (bb(3) .gt. 0) then
        if (bb(3) .gt. (explarge*R*minT)) then
          if (.not. ierr3) then
            print *,'fitkt warning: Ea is too large for machine'
            print *,'               Will try a two parameter fit.'
            ierr3 = .true.
          end if
          Ea    = explarge*R*minT
        else
          Ea    = bb(3)
        end if
      else
        if (abs(bb(3)) .gt. (explarge*R*minT)) then
          if (.not. ierr3) then
            print *,'fitkt warning: Ea is too small for machine'
            print *,'               Will try a two parameter fit.'
            ierr3 = .true.
          end if
          Ea    = - explarge*R*minT
        else
          Ea    = bb(3)
        end if
      end if
      if (.not. ierr3) go to 180
C
C-----------------------------------------------------------------------
C
C         Do two parameter fit if three parameter fit fails
C
C-----------------------------------------------------------------------
C
      aa(1,1)  =npts
      aa(1,2)  =sit
      bb(1)    =slnk
      aa(2,1)  =sit
      aa(2,2)  =sit2
      bb(2)    =slnkit
C
C             Call a linear solver to solve the 2 by 2 system
C
C*****imsl: precision > single
C      isig=10
C      call leqt2f( aa, 1, 3, 2, bb, isig, wk, ier)
C*****END imsl: precision > single
C
C*****slatec: precision > double
      call dgefs( aa, 3, 2, bb, 1, ind, wk, iwk)
      if (ind .le. 2 .and. ind .ge. 0) then
        print *,'fitkt warning: solution to fitting problem only ',
     1          'has ', ind,' significant digits.'
      else if (ind .lt. 0) then
        print *,'fitkt error: dgefs reports error, ind = ', ind
        ierr2 = .true.
      end if
C*****END slatec: precision > double
C
C*****slatec: precision > single
C      call sgefs( aa, 3, 2, bb, 1, ind, wk, iwk)
C      if (ind .le. 2 .and. ind .ge. 0) then
C        print *,'fitkt warning: solution to fitting problem only ',
C     1          'has ', ind,' significant digits.'
C      else if (ind .lt. 0) then
C        print *,'fitkt error: sgefs reports error, ind = ', ind
C        ierr2 = .true.
C      end if
C*****END slatec: precision > single
C

      if (ierr2) then
        print *,'fitkt: Will return the three parameter fit'
        go to 180
      end if
C
      aloga = bb(1)
      if (aloga .gt. explarge) then
        if (.not. ierr2) then
          print *,'fitkt warning: A factor is still ',
     1            'too large for machine.'
          print *,'               It will be adjusted downwards'
        end if
        a  = exp (explarge)
      else
        a  = exp (aloga)
      end if
C
      b = 0.0
C
      if (bb(2) .ge. 0.0) then
        if (bb(2) .gt. (explarge * R * minT)) then
          print *,'fitkt warning: Ea is still too large for machine'
          print *,'               It will be adjusted downwards.'
          Ea    = explarge*R*minT
        else
          Ea    = bb(2)
        end if
      else
        if (abs(bb(2)) .gt. (explarge * R * minT)) then
          print *,'fitkt warning: Ea is too small for machine'
          print *,'               Will try a two parameter fit.'
          Ea    = - explarge*R*minT
        else
          Ea    = bb(2)
        end if
      end if
C
C-----------------------------------------------------------------------
C
C                        Return from here
C
C-----------------------------------------------------------------------
C
180   continue
C
C*****printing > on
C      print *,'K(T) = ',a,' T**(',b,') exp ( (- ',
C     1        Ea,' cals/mole) / RT)'
C      print *,'FITKT: Comparison of fit to values:'
C      print *,'      Temperature       Values            Fit'
C      do 190 i=1,npts
C       print *,'FITKT:     ',T(i),'  ',ak(i), '     ',
C     1        a * T(i)**b * exp(-Ea/(R*T(i)))
C190   continue
C*****END printing > on
C
      return
      end
C*****precision > single
C*DECK SGEFS
C      SUBROUTINE SGEFS (A, LDA, N, V, ITASK, IND, WORK, IWORK)
CC***BEGIN PROLOGUE  SGEFS
CC***PURPOSE  Solve a general system of linear equations.
CC***LIBRARY   SLATEC
CC***CATEGORY  D2A1
CC***TYPE      SINGLE PRECISION (SGEFS-S, DGEFS-D, CGEFS-C)
CC***KEYWORDS  COMPLEX LINEAR EQUATIONS, GENERAL MATRIX,
CC             GENERAL SYSTEM OF LINEAR EQUATIONS
CC***AUTHOR  Voorhees, E. A., (LANL)
CC***DESCRIPTION
CC
CC    Subroutine SGEFS solves a general NxN system of single
CC    precision linear equations using LINPACK subroutines SGECO
CC    and SGESL.  That is, if A is an NxN real matrix and if X
CC    and B are real N-vectors, then SGEFS solves the equation
CC
CC                          A*X=B.
CC
CC    The matrix A is first factored into upper and lower tri-
CC    angular matrices U and L using partial pivoting.  These
CC    factors and the pivoting information are used to find the
CC    solution vector X.  An approximate condition number is
CC    calculated to provide a rough estimate of the number of
CC    digits of accuracy in the computed solution.
CC
CC    If the equation A*X=B is to be solved for more than one vector
CC    B, the factoring of A does not need to be performed again and
CC    the option to only solve (ITASK .GT. 1) will be faster for
CC    the succeeding solutions.  In this case, the contents of A,
CC    LDA, N and IWORK must not have been altered by the user follow-
CC    ing factorization (ITASK=1).  IND will not be changed by SGEFS
CC    in this case.
CC
CC  Argument Description ***
CC
CC    A      REAL(LDA,N)
CC             on entry, the doubly subscripted array with dimension
CC               (LDA,N) which contains the coefficient matrix.
CC             on return, an upper triangular matrix U and the
CC               multipliers necessary to construct a matrix L
CC               so that A=L*U.
CC    LDA    INTEGER
CC             the leading dimension of the array A.  LDA must be great-
CC             er than or equal to N.  (terminal error message IND=-1)
CC    N      INTEGER
CC             the order of the matrix A.  The first N elements of
CC             the array A are the elements of the first column of
CC             the  matrix A.  N must be greater than or equal to 1.
CC             (terminal error message IND=-2)
CC    V      REAL(N)
CC             on entry, the singly subscripted array(vector) of di-
CC               mension N which contains the right hand side B of a
CC               system of simultaneous linear equations A*X=B.
CC             on return, V contains the solution vector, X .
CC    ITASK  INTEGER
CC             If ITASK=1, the matrix A is factored and then the
CC               linear equation is solved.
CC             If ITASK .GT. 1, equation is solved using the existing
CC               factored matrix A and IWORK.
CC             If ITASK .LT. 1, then terminal error message IND=-3 is
CC               printed.
CC    IND    INTEGER
CC             GT. 0  IND is a rough estimate of the number of digits
CC                     of accuracy in the solution, X.
CC             LT. 0  see error message corresponding to IND below.
CC    WORK   REAL(N)
CC             a singly subscripted array of dimension at least N.
CC    IWORK  INTEGER(N)
CC             a singly subscripted array of dimension at least N.
CC
CC  Error Messages Printed ***
CC
CC    IND=-1  terminal   N is greater than LDA.
CC    IND=-2  terminal   N is less than 1.
CC    IND=-3  terminal   ITASK is less than 1.
CC    IND=-4  terminal   The matrix A is computationally singular.
CC                         A solution has not been computed.
CC    IND=-10 warning    The solution has no apparent significance.
CC                         The solution may be inaccurate or the matrix
CC                         A may be poorly scaled.
CC
CC               Note-  The above terminal(*fatal*) error messages are
CC                      designed to be handled by XERMSG in which
CC                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
CC                      for warning error messages from XERMSG.  Unless
CC                      the user provides otherwise, an error message
CC                      will be printed followed by an abort.
CC
CC***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
CC                 Stewart, LINPACK Users' Guide, SIAM, 1979.
CC***ROUTINES CALLED  R1MACH, SGECO, SGESL, XERMSG
CC***REVISION HISTORY  (YYMMDD)
CC   800317  DATE WRITTEN
CC   890531  Changed all specific intrinsics to generic.  (WRB)
CC   890831  Modified array declarations.  (WRB)
CC   890831  REVISION DATE from Version 3.2
CC   891214  Prologue converted to Version 4.0 format.  (BAB)
CC   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
CC   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
CC   920501  Reformatted the REFERENCES section.  (WRB)
CC***END PROLOGUE  SGEFS
CC
C      INTEGER LDA,N,ITASK,IND,IWORK(*)
C      REAL A(LDA,*),V(*),WORK(*),R1MACH
C      REAL RCOND
C      CHARACTER*8 XERN1, XERN2
CC***FIRST EXECUTABLE STATEMENT  SGEFS
C      IF (LDA.LT.N) THEN
C         IND = -1
C         WRITE (XERN1, '(I8)') LDA
C         WRITE (XERN2, '(I8)') N
C         CALL XERMSG ('SLATEC', 'SGEFS', 'LDA = ' // XERN1 //
C     *      ' IS LESS THAN N = ' // XERN2, -1, 1)
C         RETURN
C      ENDIF
CC
C      IF (N.LE.0) THEN
C         IND = -2
C         WRITE (XERN1, '(I8)') N
C         CALL XERMSG ('SLATEC', 'SGEFS', 'N = ' // XERN1 //
C     *      ' IS LESS THAN 1', -2, 1)
C         RETURN
C      ENDIF
CC
C      IF (ITASK.LT.1) THEN
C         IND = -3
C         WRITE (XERN1, '(I8)') ITASK
C         CALL XERMSG ('SLATEC', 'SGEFS', 'ITASK = ' // XERN1 //
C     *      ' IS LESS THAN 1', -3, 1)
C         RETURN
C      ENDIF
CC
C      IF (ITASK.EQ.1) THEN
CC
CC        FACTOR MATRIX A INTO LU
CC
C         CALL SGECO(A,LDA,N,IWORK,RCOND,WORK)
CC
CC        CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
CC
C         IF (RCOND.EQ.0.0) THEN
C            IND = -4
C            CALL XERMSG ('SLATEC', 'SGEFS',
C     *         'SINGULAR MATRIX A - NO SOLUTION', -4, 1)
C            RETURN
C         ENDIF
CC
CC        COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
CC        AND CHECK FOR IND GREATER THAN ZERO
CC
C         IND = -LOG10(R1MACH(4)/RCOND)
C         IF (IND.LE.0) THEN
C            IND=-10
C            CALL XERMSG ('SLATEC', 'SGEFS',
C     *         'SOLUTION MAY HAVE NO SIGNIFICANCE', -10, 0)
C         ENDIF
C      ENDIF
CC
CC     SOLVE AFTER FACTORING
CC
C      CALL SGESL(A,LDA,N,IWORK,V,0)
C      RETURN
C      END
C*****END precision > single
C*****precision > double
*DECK DGEFS
      SUBROUTINE DGEFS (A, LDA, N, V, ITASK, IND, WORK, IWORK)
C***BEGIN PROLOGUE  DGEFS
C***PURPOSE  Solve a general system of linear equations.
C***LIBRARY   SLATEC
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGEFS-S, DGEFS-D, CGEFS-C)
C***KEYWORDS  COMPLEX LINEAR EQUATIONS, GENERAL MATRIX,
C             GENERAL SYSTEM OF LINEAR EQUATIONS
C***AUTHOR  Voorhees, E. A., (LANL)
C***DESCRIPTION
C
C    Subroutine DGEFS solves a general NxN system of double
C    precision linear equations using LINPACK subroutines DGECO
C    and DGESL.  That is, if A is an NxN double precision matrix
C    and if X and B are double precision N-vectors, then DGEFS
C    solves the equation
C
C                          A*X=B.
C
C    The matrix A is first factored into upper and lower tri-
C    angular matrices U and L using partial pivoting.  These
C    factors and the pivoting information are used to find the
C    solution vector X.  An approximate condition number is
C    calculated to provide a rough estimate of the number of
C    digits of accuracy in the computed solution.
C
C    If the equation A*X=B is to be solved for more than one vector
C    B, the factoring of A does not need to be performed again and
C    the option to only solve (ITASK.GT.1) will be faster for
C    the succeeding solutions.  In this case, the contents of A,
C    LDA, N and IWORK must not have been altered by the user follow-
C    ing factorization (ITASK=1).  IND will not be changed by DGEFS
C    in this case.
C
C  Argument Description ***
C
C    A      DOUBLE PRECISION(LDA,N)
C             on entry, the doubly subscripted array with dimension
C               (LDA,N) which contains the coefficient matrix.
C             on return, an upper triangular matrix U and the
C               multipliers necessary to construct a matrix L
C               so that A=L*U.
C    LDA    INTEGER
C             the leading dimension of the array A.  LDA must be great-
C             er than or equal to N.  (terminal error message IND=-1)
C    N      INTEGER
C             the order of the matrix A.  The first N elements of
C             the array A are the elements of the first column of
C             the matrix A.  N must be greater than or equal to 1.
C             (terminal error message IND=-2)
C    V      DOUBLE PRECISION(N)
C             on entry, the singly subscripted array(vector) of di-
C               mension N which contains the right hand side B of a
C               system of simultaneous linear equations A*X=B.
C             on return, V contains the solution vector, X .
C    ITASK  INTEGER
C             If ITASK=1, the matrix A is factored and then the
C               linear equation is solved.
C             If ITASK .GT. 1, equation is solved using the existing
C               factored matrix A and IWORK.
C             If ITASK .LT. 1, then terminal error message IND=-3 is
C               printed.
C    IND    INTEGER
C             GT. 0  IND is a rough estimate of the number of digits
C                     of accuracy in the solution, X.
C             LT. 0  see error message corresponding to IND below.
C    WORK   DOUBLE PRECISION(N)
C             a singly subscripted array of dimension at least N.
C    IWORK  INTEGER(N)
C             a singly subscripted array of dimension at least N.
C
C  Error Messages Printed ***
C
C    IND=-1  terminal   N is greater than LDA.
C    IND=-2  terminal   N is less than 1.
C    IND=-3  terminal   ITASK is less than 1.
C    IND=-4  terminal   The matrix A is computationally singular.
C                         A solution has not been computed.
C    IND=-10 warning    The solution has no apparent significance.
C                         The solution may be inaccurate or the matrix
C                         A may be poorly scaled.
C
C               Note-  The above terminal(*fatal*) error messages are
C                      designed to be handled by XERMSG in which
C                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
C                      for warning error messages from XERMSG.  Unless
C                      the user provides otherwise, an error message
C                      will be printed followed by an abort.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  D1MACH, DGECO, DGESL, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800326  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGEFS
C
      INTEGER LDA,N,ITASK,IND,IWORK(*)
      DOUBLE PRECISION A(LDA,*),V(*),WORK(*),D1MACH
      DOUBLE PRECISION RCOND
      CHARACTER*8 XERN1, XERN2
C***FIRST EXECUTABLE STATEMENT  DGEFS
      IF (LDA.LT.N) THEN
         IND = -1
         WRITE (XERN1, '(I8)') LDA
         WRITE (XERN2, '(I8)') N
         CALL XERMSG ('SLATEC', 'DGEFS', 'LDA = ' // XERN1 //
     *      ' IS LESS THAN N = ' // XERN2, -1, 1)
         RETURN
      ENDIF
C
      IF (N.LE.0) THEN
         IND = -2
         WRITE (XERN1, '(I8)') N
         CALL XERMSG ('SLATEC', 'DGEFS', 'N = ' // XERN1 //
     *      ' IS LESS THAN 1', -2, 1)
         RETURN
      ENDIF
C
      IF (ITASK.LT.1) THEN
         IND = -3
         WRITE (XERN1, '(I8)') ITASK
         CALL XERMSG ('SLATEC', 'DGEFS', 'ITASK = ' // XERN1 //
     *      ' IS LESS THAN 1', -3, 1)
         RETURN
      ENDIF
C
      IF (ITASK.EQ.1) THEN
C
C        FACTOR MATRIX A INTO LU
C
         CALL DGECO(A,LDA,N,IWORK,RCOND,WORK)
C
C        CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
C
         IF (RCOND.EQ.0.0D0) THEN
            IND = -4
            CALL XERMSG ('SLATEC', 'DGEFS',
     *         'SINGULAR MATRIX A - NO SOLUTION', -4, 1)
            RETURN
         ENDIF
C
C        COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
C        AND CHECK FOR IND GREATER THAN ZERO
C
         IND = -LOG10(D1MACH(4)/RCOND)
         IF (IND.LE.0) THEN
            IND=-10
            CALL XERMSG ('SLATEC', 'DGEFS',
     *         'SOLUTION MAY HAVE NO SIGNIFICANCE', -10, 0)
         ENDIF
      ENDIF
C
C     SOLVE AFTER FACTORING
C
      CALL DGESL(A,LDA,N,IWORK,V,0)
      RETURN
      END
C*****END precision > double
