C     CVS $Revision: 1.1.1.1 $ created $Date: 2006/05/26 19:09:33 $
      SUBROUTINE OPPDIF 
     1 (LIN, LOUT, LINKCK, LINKMC, LRIN, LROUT, LRCRVR, NMAX,
     2  LENCWK, LENIWK, LENLWK, LENRWK, C, I, L, R)
C
C///////////////////////////////////////////////////////////////////////
C
C     MAIN PROGRAM
C
C///////////////////////////////////////////////////////////////////////
C
      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   AUTHOR*80, BANNER*80, C*16, ID*9, STRING*80, VDATE*16,
     +   VNUMBR*16
C*****precision > double
      DOUBLE PRECISION
C*****END precision > double
C*****precision > single
C      REAL
C*****END precision > single
     +   R
      EXTERNAL
     +   EXTENT, FLDRIV, POINT
      INTEGER
     +   I, ICC, ICKW, ICTOT, II, IIP, IKS,
     +   IMCW, ITOT, ITPW, J, K, KK, LAC, LENCWK, LENGTH,
     +   LENICK, LENITW, LENIWK, LENLWK, LENRCK, LENRWK,
     +   LENTWP, LIN, LINKCK, LINKMC, LMK, LOUT, LRCRVR, LRIN, LROUT,
     +   LTOT, MM, NA, NABV, NATJ, NBLW, NBUF, NCKW, NCON, ND, NDKJ,
     +   NDS, NFF, NFL, NFN, NMAX, NMCW, NOX, NPD, NS,
     +   NSCH, NSN, NSSV, NTDR, NTGV, NTOT, NTWP, NVIS, NWT, NX,
     +   NXGV, NYFL, NYOX, NYV, ASIZE, VMAX, GROUPA,
     +   GROUPB, NRKFT, NRKRT, NUSV, NRSV, CKSLEN
      LOGICAL L
      EXTERNAL CKSLEN
C
      PARAMETER (ID = 'OPPDIF:  ')
C
      PARAMETER (GROUPA = 0)
      PARAMETER (GROUPB = 0)
C
      PARAMETER (VDATE = 'DECEMBER 1997')
      PARAMETER (VNUMBR = '8.4')
C     VERSION 4.05 = 8.1 IN CVS
C     V. 8.2: CHANGE SIGN ON BC FOR AFUEL, AOXID (8/97)
C     V. 8.3: A.E. Lutz, 10-28-97
C     1. Fix bug #099: delete unused UROUND function.
C     2. Fix bug #107: remove unused variables.
C     3. Fix bug #108: remove extra parameter TI from call list
C        of CKBSEC in routine FUN.
C     4. Correct DO 2200 loop to go to KK instead of KK-1 in FUN.
C     V. 8.4: A.E. Lutz, J. Grcar 11-09-97
C     1. Fix bug #137: SN to S in call to OPPOS at line 383.
C
      DIMENSION
     +   AUTHOR(11, 2), BANNER(10), C(LENCWK), I(LENIWK), L(LENLWK),
     +   R(LENRWK)
C
      COMMON /FLFLFL/ NCKW, NMCW, NYOX, NYFL, NWT,  NFL,  NPD, NOX,
     1                NSCH, NX,   NVIS, NCON, NTGV, NXGV, ND,  NDKJ,
     2                NTDR, NYV,  NABV, NBLW, NBUF, NTWP, NS,  NSN,
     3                NFF,  NFN,  NDS,  NSSV, NA,
     4                NRKFT, NRKRT, NUSV, NRSV,
     5                ICKW, IMCW, ITPW, IKS, ICC, IIP, LAC, LMK
C      DATA I/LENIWK*0/, L/LENLWK*.FALSE./, C/LENCWK*' '/,
C     1     R/LENRWK*0.0/
C
C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////
C
C
C///  PRINT THE HEADER.
C
C*****precision > double
      STRING = 'DOUBLE'
C*****END precision > double
C*****precision > single
C      STRING = 'SINGLE'
C*****END precision > single
     +   // ' PRECISION VERSION ' // VNUMBR // ' OF ' // VDATE
C
      CALL TWSQEZ (LENGTH, STRING)
      WRITE (LOUT, 10001) ID, STRING (1 : LENGTH)
C
      WRITE (LOUT, 10002)
      WRITE (LOUT, 10002)
C
      DO 1010 J = 1, 9
         STRING = BANNER(J)
         CALL EXTENT (LENGTH, STRING)
         LENGTH = MAX (CKSLEN(STRING), 1)
         WRITE (LOUT, 10003) STRING (1 : LENGTH)
1010  CONTINUE
C
      DO 1020 J = 1, 2
         WRITE (LOUT, 10002)
         WRITE (LOUT, 10002)
         DO 1020 K = 1, 11
            STRING = AUTHOR(K, J)
            CALL EXTENT (LENGTH, STRING)
            LENGTH = MAX (CKSLEN(STRING), 1)
            WRITE (LOUT, 10003) STRING (1 : LENGTH)
1020  CONTINUE
C
C///////////////////////////////////////////////////////////////////////
C
C     (2) SET UP INTERNAL WORK POINTERS.
C
C///////////////////////////////////////////////////////////////////////
C
      WRITE (LOUT, 10004)
C
      REWIND LROUT
C
      CALL POINT (LINKCK, LINKMC, NMAX, LOUT, LROUT, KK, II, MM, NATJ,
     1            LENIWK, LENRWK, LENCWK, LENITW, LENTWP, LENICK,
     2            LENRCK, LTOT, ITOT, NTOT, ICTOT, I, R, C,
     3            ASIZE, GROUPA, GROUPB, VMAX)
C
      WRITE (LOUT, 10005) ID,
     +   ICTOT, ITOT, LTOT, NTOT,
     +   LENCWK - ICTOT, LENIWK - ITOT, LENLWK - LTOT, LENRWK - NTOT,
     +   LENCWK, LENIWK, LENLWK, LENRWK
C
      IF (LTOT .GT. LENLWK .OR. ITOT .GT. LENIWK .OR.
     +    NTOT .GT. LENRWK .OR. ICTOT .GT. LENCWK) GO TO 9201
C
      WRITE (LOUT, 10004)
C
      CALL FLDRIV (LENITW, LENTWP, LIN, LOUT, LINKCK, LINKMC, LRIN,
     1           LROUT, LRCRVR, KK, II, MM, NATJ, NMAX, LENICK, LENRCK,
     2           R(NCKW), R(NMCW), R(NYOX), R(NYFL), R(NWT),  R(NFL),
     3           R(NPD),  R(NOX),  R(NSCH), R(NX),   R(NVIS), R(NCON),
     4           R(NTGV), R(NXGV), R(ND),   R(NDKJ), R(NTDR), R(NYV),
     5           R(NABV), R(NBLW), R(NBUF), R(NTWP), R(NS),   R(NSN),
     6           R(NFF),  R(NFN),  R(NDS),  R(NSSV), R(NA),
     *           R(NUSV), R(NRSV), R(NRKFT), R(NRKRT),
     7           I(ICKW), I(IMCW), I(ITPW), C(IKS), C(ICC),
     8           I(IIP),  L(LAC), L(LMK),
     9           ASIZE, GROUPA, GROUPB, VMAX)
C
C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////
C
      DATA BANNER /
     + ' OOOOO  PPPPPP  PPPPPP  DDDDDD    III   FFFFFFF',
     + 'O     O P     P P     P D     D    I    F',
     + 'O     O P     P P     P D     D    I    F',
     + 'O     O PPPPPP  PPPPPP  D     D    I    FFFFF',
     + 'O     O P       P       D     D    I    F',
     + 'O     O P       P       D     D    I    F',
     + ' OOOOO  P       P       DDDDDD    III   F',
     +   ' ',
     +   ' ',
     +   'ONE DIMENSIONAL OPPOSED-FLOW DIFFUSION FLAME CODE' /
C
      DATA (AUTHOR(J, 1), J = 1, 11) /
     +   'WRITTEN BY:',
     +   ' ',
     +   '   DR. ANDREW E. LUTZ',
     +   '   SANDIA NATIONAL LABORATORIES',
     +   '   THERMAL & PLASMA PROCESSES DEPARTMENT 8345',
     +   '   MS-9042',
     +   '   LIVERMORE, CA 94551-0969 USA',
     +   '   aelutz@california.sandia.gov' ,' ', ' ', ' ' /
C
      DATA (AUTHOR(J, 2), J = 1, 11) /
     +   'WITH THE ASSISTANCE OF:',
     +   '   DR. JOSEPH F. GRCAR',
     +   '   DR. ROBERT J. KEE',
     +   '   MS. FRAN M. RUPLEY', ' ', ' ', ' ', ' ', ' ', ' ',' '/
C
10001 FORMAT
     +  (/9X, 35(' /'),
     +  //1X, A9, A)
C
10002 FORMAT
     +   ()
C
10003 FORMAT
     +  (10X, A)
C
10004 FORMAT
     +  (/9X, 35(' /'))
C
10005 FORMAT
     + (/1X, A9, 'THE WORK SPACE REQUIREMENTS ARE AS FOLLOWS.'
     + //10X, '           CHARACTER    INTEGER    LOGICAL       REAL'
     + //10X, '     USED', 4(2X, I9)
     +  /10X, '   EXCESS', 4(2X, I9)
     +  /10X, '    TOTAL', 4(2X, I9))
C
C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////
C
      GO TO 99999
C
9201  WRITE (LOUT, 99201) ID
      GO TO 99999
C
99201 FORMAT
     +   (/1X, A9, 'ERROR.  THERE IS NOT ENOUGH WORK SPACE.')
C
C///  EXIT.
C
99999 CONTINUE
C
C     end PROGRAM OPPDIF
      END
C
      SUBROUTINE EXTENT (LENGTH, STRING)
C
C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     EXTENT
C
C     RETURN THE LENGTH OF A STRING.  THAT IS, RETURN THE POSITION OF
C     THE LAST NONBLANK CHARACTER, OR 1 IF ALL ARE BLANK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER STRING*(*)
      INTEGER J, LENGTH

      LENGTH = 1
      DO 0100 J = 1, LEN (STRING)
         IF (STRING(J : J) .NE. ' ') LENGTH = J
0100  CONTINUE
C
C     end of SUBROUTINE EXTENT
      RETURN
      END
C
      SUBROUTINE FLDRIV
     +  (LENITW, LENTWP, LIN, LOUT, LINKCK, LINKMC, LRIN, LROUT, LRCRVR,
     +   KK, II, MM, NATJ, NMAX, LENICK, LENRCK, RCKWRK, RMCWRK, YOXID,
     +   YFUEL, WT, FUEL, PROD, OXID, SCRCHK, X, VISC, COND, TGIVEN,
     +   XGIVEN, D, DKJ, TDR, YV, ABOVE, BELOW, BUFFER, TWPWK, S, SN, F,
     +   FN, DS, SSAVE, A, USAVE, RSAVE, RKFT, RKRT,
     +   ICKWRK, IMCWRK, ITWPWK, KSYM, CCKWRK, IP,
     +   ACTIVE, MARK, ASIZE, GROUPA, GROUPB, VMAX)
C
C///////////////////////////////////////////////////////////////////////
C
C     D I F F U S : FLDRIV
C
C     MAIN SUBROUTINE
C
C///////////////////////////////////////////////////////////////////////
C
      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   CCKWRK*(*),
     +   ISOLUT*16, KSYM*(*), NAME*16, REPORT*16,
     +   SIGNAL*16, VERSIO*80, ID*9, STRING*80
C*****precision > double
      DOUBLE PRECISION
C*****END precision > double
C*****precision > single
C      REAL
C*****END precision > single
     +   A, ABOVE, AFUEL, AOXID, BELOW, BUFFER, COND,
     +   CONDIT, D, DKJ, DS, DT, F,
     +   FN, FUEL, OXID, P, PATM, PCTADP, PROD, RATEF,
     +   RATGTC, RCKWRK, RFAC, RMCWRK,
     +   RU, RUC, S, SCRCHK, SFLR, SN, SPOS, SSAVE, TDR, TFUEL,
     +   TGIVEN, TMAX, TOXID, TWPWK, VFUEL, VISC, VOXID,
     +   WMIX, WNDFAC, WT, X, XCEN, XEND, XGIVEN, YFUEL,
     +   YOXID, YV, RKFT, RKRT, USAVE, RSAVE,
     +   DT2, CKBSEC
      EXTERNAL
     +   CKRP, CKSYMS, CKWT, FUN, PRINT, RDKEY, REGRID, RESTRT,
     +   START, TWCOPY, TWOPNT, TWPREP, TWSOLV, CKBSEC
      INTEGER
     +   ASIZE, CALL, CALLS, GROUPA, GROUPB, ICKWRK,
     +   II, IMCWRK, INDEX, IP, IPROFL, IRETIR,
     +   ITWPWK, J, JJ, JJREGD, K, KK, KW,
     +   LENICK, LENITW, LENRCK, LENTWP,
     +   LIN, LINKCK, LINKMC, LOUT, LRCRVR, LRIN, LROUT, MM, N,
     +   N1CALL, NATJ, NMAX, NT,
     +   NTEMP, NW, NY, NYS, VMAX, VRBL, VRBLS, LENGTH,
     +   NDPR, NF, NG, NH, IRATE, JJNEW, JJOLD,
     +   ISTPS2
      INTRINSIC
     +   MAX, MIN, SQRT
      LOGICAL
     +   ACTIVE, ENERGY, ERROR, IERR,
     +   LENRGY, LMULTI, LREGRD, LRSTRT, LSEN, LTDIF,
     +   LTIME, LUMESH, LUSTGV, LVARMC, MARK,
     +   RETURN, TRYGRD
C
      PARAMETER (ID = 'FLDRIV:  ')
C
      DIMENSION
     +   A(ASIZE), ABOVE(GROUPA + NATJ + GROUPB), ACTIVE(*),
     +   BELOW(GROUPA + NATJ + GROUPB), BUFFER(VMAX),
     +   COND(*), D(KK,*), DKJ(KK,KK,*), DS(*), F(VMAX), FN(VMAX),
     +   FUEL(*), ICKWRK(*), IMCWRK(*), IP(VMAX),
     +   ITWPWK(*), KW(6),
     +   OXID(*), PROD(*), RCKWRK(*),
     +   RMCWRK(*), S(VMAX), SCRCHK(KK,5), SN(VMAX),
     +   SSAVE(*), TDR(KK,*), TGIVEN(*), TWPWK(*), VISC(*), WT(*), X(*),
     +   XGIVEN(*), YFUEL(KK), YOXID(KK), YV(KK,*),
     +   CCKWRK(*), KSYM(*), RKFT(*), RKRT(*), RSAVE(*), USAVE(*)
C
      COMMON /LOCS/ NT, NG, NF, NH, NYS, NY
C
      DATA ISOLUT /'SOLUTION        '/
C
      DATA RFAC/1.0/, ERROR/.FALSE./
C
C     get some things from CHEMKIN
      CALL CKSYMS (CCKWRK, LOUT, KSYM, IERR)
      CALL CKWT   (ICKWRK, RCKWRK, WT)
      CALL CKRP   (ICKWRK, RCKWRK, RU, RUC, PATM)
C
C     read the keyword input
      CALL RDKEY (KK, NMAX, NATJ, LIN, LOUT, KSYM, PATM, LUSTGV,
     1            LENRGY, LMULTI, LTDIF, LUMESH, LRSTRT,
     2            IPROFL, RATEF, N1CALL, LSEN, VOXID, VFUEL,
     3            AOXID, AFUEL, TFUEL, TOXID, TMAX, P, JJ, X, FUEL,
     4            PROD, OXID, XCEN, XEND, WMIX, IRETIR, SFLR, NTEMP,
     6            XGIVEN, TGIVEN, WNDFAC, LREGRD, JJREGD, PCTADP,
     7            RATGTC, KW, NW, SPOS, NDPR, DT2, ISTPS2)
C
C///  SET SOME TWOPNT ARRAYS.
C
C     solution bounds (these should take group A and B into account)
      BELOW(NT) = 200.
      ABOVE(NT) = 6000.
      BELOW(NG) = -10000.
      ABOVE(NG) =  10000.
      BELOW(NF) = -10000.
      ABOVE(NF) =  10000.
      BELOW(NH) = -1.E10
      ABOVE(NH) =  1.E10
      DO 1020 K = 1, KK
         N = NYS + K
         BELOW (N) = SFLR
         ABOVE (N) = 10.
1020  CONTINUE
C
C     variables on which to adapt
      DO 1030 N = 1, NATJ
         ACTIVE(N) = .TRUE.
1030  CONTINUE
C
C     logical for regrid attempts
      TRYGRD = .FALSE.
C
      IF (.NOT. LRSTRT) THEN
C        set up an original problem
         CALL START (KK, NMAX, NATJ, JJ, LOUT, LUMESH, LENRGY,
     1               IPROFL, TMAX, TFUEL, TOXID, VFUEL, VOXID,
     2               P, FUEL, PROD, OXID, NTEMP, XGIVEN,
     3               TGIVEN, XCEN, XEND, WMIX, ICKWRK, RCKWRK,
     4               SCRCHK(1,4), YOXID, YFUEL, X, S)
C
      ELSE

C///     READ THE RESTART FILE.

            CALL OPREAD (LRIN, LOUT, GROUPA, GROUPB, II,
     1                   JJ, NATJ, X, S, ERROR)
            IF (ERROR) GO TO 99999
C
C        limit species minimum on restart
         DO 3050 J = 1, JJ
            VRBL = GROUPA + NATJ * (J - 1) + NYS + 1
            CALL OPPOS (S(VRBL), KK, SPOS)
3050     CONTINUE
C
         CALL RESTRT (KK, NMAX, NATJ, JJ, LOUT, LUMESH, LUSTGV,
     1                FUEL, SCRCHK(1,2), OXID, NTEMP, XGIVEN,
     2                TGIVEN, XCEN, XEND, WMIX, ICKWRK, RCKWRK,
     3                SCRCHK(1,4), COND, YOXID, YFUEL, X, S)
C
         IF (LREGRD) THEN
C           regrid
            WRITE (STRING, '(A, I10, A)')
     +         'REGRIDDING TO ', JJREGD, ' POINTS'
            CALL TWSQEZ (LENGTH, STRING)
            WRITE (LOUT, 10001) ID, STRING (1 : LENGTH)
C
C           do not allow regrid to change 2 points at each boundary
            VRBL = GROUPA + NATJ + 1
            JJNEW = JJREGD - 2
            JJOLD = JJ - 2
            CALL REGRID (SN(VRBL), S(VRBL), NATJ, JJNEW, JJOLD,
     1                   VISC(2), X(2), COND, PCTADP, RATGTC, NT)
            JJ    = JJNEW + 2
            JJOLD = JJOLD + 2
C
C           COPY NEW GRID AND SOLUTION
C
C           copy last point
            X(JJ) = X(JJOLD)
            VRBL = GROUPA + NATJ * (JJ - 1)
            N    = GROUPA + NATJ * (JJOLD - 1)
            CALL TWCOPY (NYS+KK, S(N+1), S(VRBL+1))
C
C           copy interior points
            DO 3080 J = 2, JJ-1
               X(J) = VISC(J)
               VRBL = GROUPA + NATJ * (J - 1)
               CALL TWCOPY (NYS, SN(VRBL+1), S(VRBL+1))
               CALL TWCOPY (KK,  SN(VRBL+NYS+1), S(VRBL+NYS+1))
               CALL OPPOS  (S(VRBL+NYS+1), KK, SPOS)
3080        CONTINUE
C
C           correct negative species at endpoints
            VRBL = GROUPA + 1
            CALL OPPOS (S(VRBL), KK, SPOS)
            VRBL = GROUPA + NATJ * (JJ - 1) + 1
            CALL OPPOS (S(VRBL), KK, SPOS)
C
            IF (.NOT. LUSTGV) THEN
C              set given temperatures to solution
               NTEMP = JJ
               DO 3100 J = 1, NTEMP
                  XGIVEN(J) = X(J)
                  VRBL = GROUPA + NATJ * (J - 1) + NT
                 TGIVEN(J) = S(VRBL)
3100           CONTINUE
            ENDIF
         ENDIF
      ENDIF
C
C///  TOP OF THE LOOP OVER SUBPROBLEMS.
C
C     save the energy flag
      ENERGY = LENRGY
C
C     choose the number of subproblems
      IF (ENERGY) THEN
         CALLS = 2
      ELSE
         CALLS = 1
      ENDIF
C
      DO 5000 CALL = N1CALL, CALLS
C
         IF (CALL .EQ. 1) THEN
            LENRGY = .FALSE.
         ELSEIF (CALL .EQ. 2) THEN
            LENRGY = ENERGY
            WRITE (LOUT, 10002) ID
         ENDIF
C
         IF (LENRGY) RFAC = RATEF
C
C///     TOP OF THE LOOP TO RE-ENTER TWOPNT.
C
C        specify the expected TWOPNT version
C*****precision > double
         VERSIO = 'DOUBLE PRECISION VERSION 3.22'
C*****END precision > double
C*****precision > single
C         VERSIO = 'SINGLE PRECISION VERSION 3.22'
C*****END precision > single
C
         IF (CALL .LT. CALLS) THEN
            INDEX = 1
         ELSE
            INDEX = 2
            CALL TWSETI (ERROR, LOUT, 'STEPS1', ISTPS2)
            CALL TWSETR (ERROR, LOUT, 'STRID0', DT2)
            CALL TWSETL (ERROR, LOUT, 'ADAPT', .TRUE.)
         ENDIF
         SIGNAL = ' '
4020     CONTINUE
C
         CALL TWOPNT (ERROR, LOUT, VERSIO, ABOVE, ACTIVE, BELOW,
     1             BUFFER, NATJ, CONDIT, GROUPA, GROUPB, LENITW,
     3             ITWPWK, MARK, NAME, 1, NMAX, JJ, REPORT,
     5             LENTWP, TWPWK, SIGNAL, DT, LTIME, S, X)
         IF (ERROR) STOP
C
C        update the number of variables
         VRBLS = GROUPA + NATJ * JJ + GROUPB
C
         IF (SIGNAL .EQ. 'RESIDUAL') THEN
C           evaluate the residual
            IRATE = 1
            LVARMC = .TRUE.
            CALL FUN (II, KK, JJ, NATJ, LENRGY, LMULTI, LTDIF, LVARMC,
     1                LTIME, P, WT, YOXID, YFUEL, DT, NTEMP, XGIVEN,
     2                TGIVEN, X, SN, BUFFER, WNDFAC, SCRCHK(1,1),
     3               SCRCHK(1,5), YV, SCRCHK(1,2), SCRCHK(1,3),
     4               SCRCHK(1,4), TFUEL, TOXID, VFUEL, VOXID,
     5               AFUEL, AOXID, VISC, COND, D, DKJ, TDR,
     6               ICKWRK, RCKWRK, IMCWRK, RMCWRK, F, RFAC,
     7               RKFT, RKRT, USAVE, RSAVE, IRATE)
            CALL TWCOPY (VRBLS, F, BUFFER)
C
         ELSEIF (SIGNAL .EQ. 'PREPARE') THEN
C           prepare the Jacobian matrix
C
            IRATE = 2
            LVARMC = .TRUE.
            RETURN = .FALSE.
4030        CONTINUE

            CALL TWPREP (ERROR, LOUT, A, ASIZE, BUFFER, NATJ,
     1                   CONDIT, GROUPA, GROUPB, IP, JJ, RETURN)
            IF (ERROR) STOP
C
            IF (RETURN) THEN
               CALL FUN (II, KK, JJ, NATJ, LENRGY, LMULTI, LTDIF,
     1                   LVARMC, LTIME, P, WT, YOXID, YFUEL, DT,
     2                   NTEMP, XGIVEN, TGIVEN, X, SN, BUFFER,
     3                   WNDFAC, SCRCHK(1,1), SCRCHK(1,5), YV,
     4                   SCRCHK(1,2), SCRCHK(1,3), SCRCHK(1,4),
     5                   TFUEL, TOXID, VFUEL, VOXID, AFUEL,
     6                   AOXID, VISC, COND, D, DKJ, TDR, ICKWRK,
     7                   RCKWRK, IMCWRK, RMCWRK, F, RFAC, RKFT,
     8                   RKRT, USAVE, RSAVE, IRATE)
C
               IRATE = 3
               LVARMC = .FALSE.
               CALL TWCOPY (VRBLS, F, BUFFER)
               GO TO 4030
            ENDIF
C
         ELSEIF (SIGNAL .EQ. 'SOLVE') THEN
C           solve the linear equations
C
            CALL TWSOLV (ERROR, LOUT, A, ASIZE, BUFFER, NATJ,
     1                   GROUPA, GROUPB, IP, JJ)
            IF (ERROR) STOP
C
         ELSEIF (SIGNAL .EQ. 'RETAIN') THEN
C           retain the solution for time integration
            CALL TWCOPY (VRBLS, BUFFER, SN)
            IF (SPOS .GE. 0.0) THEN
               DO 4040 J = 1, JJ
                  VRBL = GROUPA + NATJ * (J - 1) + NYS
                  CALL OPPOS (SN(VRBL+1), KK, SPOS)
4040           CONTINUE
            ENDIF
C
         ELSEIF (SIGNAL .EQ. 'SHOW') THEN
C           show the solution
            CALL PRINT (LOUT, KK, JJ, NATJ, P, X, BUFFER, SCRCHK(1,1),
     1              KSYM, ICKWRK, RCKWRK, KW, NW)
C
         ELSEIF (SIGNAL .EQ. 'SAVE') THEN
C           save the solution
            REWIND LRCRVR
            WRITE (LRCRVR) ISOLUT
            WRITE (LRCRVR) NATJ, JJ, P, AFUEL, AOXID, VFUEL, VOXID
            WRITE (LRCRVR) (X(J), J=1,JJ)
            WRITE (LRCRVR) (BUFFER(J), J=1, VRBLS)
C
         ELSEIF (SIGNAL .EQ. 'UPDATE') THEN
C           update the solution to a new grid
            IF (.NOT. LENRGY) THEN
               DO 4050 J = 1, JJ
                  VRBL = GROUPA + NATJ * (J - 1) + NT
                  BUFFER(VRBL) = CKBSEC(NTEMP, X(J), XGIVEN, TGIVEN)
4050           CONTINUE
            ENDIF
C
C        bottom of the blocks to service requests from TWOPNT
      ENDIF

C///  BOTTOM OF THE LOOP TO RE-ENTRY TWOPNT.

      IF (SIGNAL .NE. ' ') GO TO 4020
C
C///  BOTTOM OF THE LOOP OVER SUBPROBLEMS.
C
C     check for solution to fixed T problem
      IF ((INDEX .EQ. 1) .AND. (REPORT .NE. ' ')) THEN
         WRITE (LOUT,'(5X,A)') 'DID NOT FINISH FIXED TEMPERATURE CASE!'
         STOP
      ENDIF
C
C     check for failure; regrid, and try again
      IF ( (INDEX .EQ. 2) .AND. ((REPORT .EQ. 'SOME SOLVED') .OR.
     1                           (REPORT .EQ. 'NOT SOLVED'))
     2                    .AND.  .NOT. TRYGRD ) THEN
C
            WRITE (STRING, '(A, I10, A)')
     +         'REGRIDDING TO ', JJ, ' POINTS'
            CALL TWSQEZ (LENGTH, STRING)
            WRITE (LOUT, 10001) ID, STRING (1 : LENGTH)
C
C           do not allow regrid to change 2 points at each boundary
            VRBL = GROUPA + NATJ + 1
            JJOLD = JJ - 2
            JJNEW = JJOLD
            CALL REGRID (SN(VRBL), S(VRBL), NATJ, JJNEW, JJOLD,
     1                   VISC(2), X(2), COND, PCTADP, RATGTC, NT)
C
C           copy new grid and solution
            JJ = JJNEW + 2
            DO 4080 J = 2, JJ-1
               X(J) = VISC(J)
               VRBL = GROUPA + NATJ * (J - 1)
               CALL TWCOPY (NYS, SN(VRBL+1), S(VRBL+1))
               DO 4080 K = NYS + 1, NYS + KK
                  S(VRBL+K) = MAX (SPOS, SN(VRBL+K))
4080        CONTINUE

            WRITE (LOUT, '(//5X,A//)')
     1      ' ATTEMPTING SOLUTION AFTER REGRID'
            TRYGRD = .TRUE.
            SIGNAL = ' '
            REPORT = ' '
            GO TO 4020
C
      ENDIF
C
5000  CONTINUE
C
C     write the completed solution
      WRITE (LROUT) ISOLUT
      WRITE (LROUT) NATJ, JJ, P, AFUEL, AOXID, VFUEL, VOXID
      WRITE (LROUT) (X(J), J=1,JJ)
      WRITE (LROUT) (S(J), J=1,VRBLS)
C
C
C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +  (/1X, A, A)

10002 FORMAT
     +  (/1X, A, 'FINISHED FIXED TEMPERATURE CASE, ',
     +           'ADDING THE ENERGY EQUATION.')

C
C///  EXIT.
C
99999 CONTINUE
C
C     end of SUBROUTINE FLDRIV
      RETURN
      END
C
      SUBROUTINE FUN (II, KK, JJ, NATJ, LENRGY, LMULTI, LTDIF, LVARMC,
     1                LTIME, P, WT, YOXID, YFUEL, DT, NTEMP, XGIVEN,
     2                TGIVEN, X, SN, S, WNDFAC, YAV, XAV, YV, WDOT, CP,
     3                H, TFUEL, TOXID, VFUEL, VOXID, AFUEL, AOXID,
     4                VISC, COND, D, DKJ, TDR, ICKWRK, RCKWRK, IMCWRK,
     5                RMCWRK, F, RFAC, RKFT, RKRT, USAVE, RSAVE, IRATE)
C
C  START PROLOGUE
C     This subroutine forms the function for the differential equations
C     Indicies:
C      NT : Temperature
C      NG : Velocity gradient
C      NF : Velocity
C      NH : Eigenvalue, pressure gradient
C      KK : no. of species
C      JJ : no. of gridpoints
C      NATJ : no. variables at each point
C      LENERGY = .TRUE. if solve energy equation
C      LMULTI
C      LTDIFF
C      LVARMC
C      LTIME
C      VFUEL, VOXID, TFUEL, TOXID, YFUEL, YOXID...
C  END PROLOGUE
C---------------------------------------------------------------------
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NG, NF, NH, NYS, NY
C
      DIMENSION WT(KK), YOXID(KK), YFUEL(KK), XGIVEN(*), TGIVEN(*),
     1          X(JJ),
     1          S(NATJ,JJ), SN(NATJ,JJ), YAV(KK), XAV(KK), YV(KK,*),
     2          WDOT(KK), CP(KK), H(KK), VISC(JJ), COND(*), D(KK,JJ),
     3          TDR(KK,*), DKJ(KK,KK,JJ), ICKWRK(*), RCKWRK(*),
     4          IMCWRK(*), RMCWRK(*), F(NATJ,JJ),
     5          RKFT(II,JJ), RKRT(II,JJ), USAVE(NATJ,JJ), RSAVE(KK,JJ)
      LOGICAL LENRGY, LMULTI, LTDIF, LVARMC, LTIME, MATCH
C
      IF (LVARMC .AND. .NOT.LTIME) THEN
C        evaluate and store the transport coefficients
         CALL MTRNPR (KK, JJ, NATJ, LENRGY, LMULTI, LTDIF, P, X, S,
     1                    YAV, ICKWRK, RCKWRK, IMCWRK, RMCWRK, WT,
     2                    H, CP, XAV, VISC, COND, D, TDR, DKJ)
      ENDIF
C
C     evaluate and store the diffusion velocities
C     (in the following call h(*) and cp(*) are used temporarily
C     for acratch space
      CALL MDIFV  (KK, JJ, NATJ, LMULTI, LTDIF, X, S, WT, YAV,
     1             H, CP, P, D, TDR, ICKWRK, RCKWRK, YV, DKJ)
C
C     fuel boundary (X=0)
      XMID = 0.5 * (X(1)+X(2))
      TAV = 0.5 * (S(NT,1) + S(NT,2))
      CALL CKAVG (KK, S(NYS+1,1), S(NYS+1,2), YAV)
      CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOP)
C
      RFUEL = RHOP * VFUEL
      DO 30 K = 1, KK
         N = NYS + K
         F(N, 1) = RFUEL * (YFUEL(K) - S(N,1)) - RHOP*YV(K,1)
30    CONTINUE
C
      IF (LENRGY) THEN
         F(NT,1) = TFUEL - S(NT,1)
      ELSE
         F(NT, 1) = S(NT,1) - CKBSEC(NTEMP, X(1), XGIVEN, TGIVEN)
      ENDIF
      F(NF,1) = S(NF,1) - VFUEL*RHOP/2.
      F(NG,1) = AFUEL*RHOP + S(NG,1)
      F(NH,1) = S(NH,2) - S(NH,1)
C
C     interior mesh points
      DO 1000 J = 2, JJ-1
        JP1 = J + 1
        JM1 = J - 1
        TAV = 0.5 * (S(NT,J) + S(NT,JP1))
        CALL CKAVG (KK, S(NYS+1,J), S(NYS+1,JP1), YAV)
C
        RHOM = RHOP
        CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOP)
        CALL CKRHOY (P, S(NT,J), S(NY,J), ICKWRK, RCKWRK, RHOJ)
C
        XMID = 0.5 * (X(J)+X(JP1))
        XP = XMID
        XM = 0.5 * (X(J) + X(JM1))
C
C       form the chemical rate terms
        IF (IRATE .EQ. 1) THEN
C          normal call to rate routine
           CALL CKWYP (P, S(NT,J), S(NY,J), ICKWRK, RCKWRK, WDOT)
C
        ELSEIF (IRATE .EQ. 2) THEN
C          save all of the solution components at this node
C
           CALL TWCOPY (NATJ, S(1,J), USAVE(1,J))
C
C          get the forward and reverse rate constants
           CALL CKKFRT (P, S(NT,J), ICKWRK, RCKWRK, RKFT(1,J),
     1                      RKRT(1,J) )
C          calculate molar production rates
           CALL CKWYPK (P, S(NT,J), S(NY,J), RKFT(1,J), RKRT(1,J),
     1                      ICKWRK, RCKWRK, WDOT)
C          save the production rates for all species at this node
           CALL TWCOPY (KK, WDOT, RSAVE(1,J))
C
        ELSE
C
           MATCH = S(NT,J) .EQ. USAVE(NT,J)
           DO 110 K = 1, KK
              N = NYS + K
              MATCH = MATCH .AND. (S(N,J) .EQ. USAVE(N,J))
110        CONTINUE
C
C          if MATCH is .TRUE. after these two tests, then neither
C          temperature nor the species concentrations have changed,
C          so don't need to recalculate the production rates
           IF (MATCH) THEN
              CALL TWCOPY (KK, RSAVE(1,J), WDOT)

           ELSEIF ( S(NT,J) .EQ. USAVE(NT,J) ) THEN
C             the temperature has not changed, so don't need to
C             recalculate the rate constants to get new production
C             rates
              CALL CKWYPK (P, S(NT,J), S(NY,J), RKFT(1,J),
     1                         RKRT(1,J), ICKWRK, RCKWRK, WDOT)
C
           ELSE
C             temperature has changed, so calculate new rate
C             constants and production rates
              CALL CKWYP (P, S(NT,J), S(NY,J), ICKWRK, RCKWRK, WDOT)
C
           ENDIF
        ENDIF
C
C       accelerate rates to facilitate solution
        IF (RFAC .NE. 1.0) THEN
           DO 120 K = 1, KK
              WDOT(K) = WDOT(K) *RFAC
120        CONTINUE
        ENDIF
C
C       form the mesh differences
        DXP =        (X(JP1) - X(J)  )
        DXM =        (X(J)   - X(JM1))
        DXAV = 0.5 * (X(JP1) - X(JM1))
        DXPM =       (X(JP1) - X(JM1))
C
C       form the coefficients for central differences
         CENDFM = - DXP / (DXM*DXPM)
         CENDFC =   (DXP-DXM) / (DXP*DXM)
         CENDFP =   DXM / (DXP*DXPM)
C
C       windward differencing factors
        IF (WNDFAC .EQ. 1.0) THEN
           IF (S(NF,J) .LE. 0.) THEN
              WINDA =  1.
              WINDB = -1.
              WINDC =  0.
           ELSE
              WINDA =  0.
              WINDB =  1.
              WINDC = -1.
           ENDIF
           DXWIND = WINDA*X(JP1) + WINDB*X(J) + WINDC*X(JM1)
C        ENDIF
C
C       species conservation equation
C        IF (WNDFAC .EQ. 1.0) THEN
          DO 200 K = 1, KK
             N = NYS + K
             DYK = WINDA*S(N,JP1) + WINDB*S(N,J) + WINDC*S(N,JM1)
             F(N,J) =  2.0*S(NF,J) * DYK / DXWIND
     1                  -  WDOT(K) * WT(K)
     2                  +  (RHOP*YV(K,J) - RHOM*YV(K,JM1)) / DXAV
200       CONTINUE
        ELSE
          DO 210 K = 1, KK
             N = NYS + K
             F(N,J) = 2.0*S(NF,J) *
     1                    (CENDFP*S(N,JP1) + CENDFC*S(N,J) +
     1                     CENDFM*S(N,JM1) )  -
     2                    WDOT(K) * WT(K) +
     3                      (RHOP*YV(K,J) - RHOM*YV(K,JM1)) / DXAV
210       CONTINUE
        ENDIF
C
C       streamfunction equation
        F(NF,J) = (S(NG,J)+S(NG,JM1))/2. - (S(NF,J) - S(NF,JM1))/DXM
C
C       momentum equation
        CALL CKRHOY (P, S(NT,JM1), S(NY,JM1), ICKWRK, RCKWRK, RHOJM1)
        CALL CKRHOY (P, S(NT,JP1), S(NY,JP1), ICKWRK, RCKWRK, RHOJP1)
        BIGCP = VISC(J)
        BIGCM = VISC(JM1)
        F(NG,J) = - 3.0*S(NG,J)**2/RHOJ
     3             - ( BIGCP * (S(NG,JP1)/RHOJP1 - S(NG,J)/RHOJ) / DXP -
     3                 BIGCM * (S(NG,J)/RHOJ -S(NG,JM1)/RHOJM1)/DXM ) /
     3                 DXAV - S(NH,J)
        IF (WNDFAC .EQ. 1.0) THEN
           DUGR = WINDA*( S(NF,JP1)*S(NG,JP1)/RHOJP1 )
     1          + WINDB*( S(NF,J  )*S(NG,J  )/RHOJ )
     2          + WINDC*( S(NF,JM1)*S(NG,JM1)/RHOJM1 )
           F(NG,J) = F(NG, J) + 2.0 * DUGR /DXWIND
        ELSE
           F(NG,J) = F(NG, J) +
     2               2.0 * (CENDFP*S(NG,JP1)*S(NF,JP1)/RHOJP1 +
     2                     CENDFC*S(NG,J)*S(NF,J)/RHOJ +
     2                     CENDFM*S(NG,JM1)*S(NF,JM1)/RHOJM1 )
        ENDIF
C
        IF (LENRGY) THEN
C          energy equation
C
C          thermodynamic properties
           CALL CKHML (S(NT,J), ICKWRK, RCKWRK, H)
           CALL CKCPBS (S(NT,J), S(NY,J), ICKWRK, RCKWRK, CPB)
           CALL CKCPMS (S(NT,J), ICKWRK, RCKWRK, CP)
C
           SUM = 0.0
           TDOT = 0.0
           CEND = CENDFP*S(NT,JP1) + CENDFC*S(NT,J) +
     1            CENDFM*S(NT,JM1)
           DO 400 K = 1, KK
              TDOT = TDOT + WDOT(K) * H(K)
              SUM = SUM + 0.5 * (YV(K,J) + YV(K,JM1)) * CP(K)*CEND
400        CONTINUE
C
           CFAC = -( COND(J)*(S(NT,JP1)-S(NT,J))/DXP -
     2              COND(JM1)*(S(NT,J)-S(NT,JM1))/DXM ) /
     3                                            (CPB*DXAV) +
     4               RHOJ * SUM / CPB + TDOT/CPB
C
           IF (WNDFAC .EQ. 1.0) THEN
             DTEMP = WINDA*S(NT,JP1) + WINDB*S(NT,J) + WINDC*S(NT,JM1)
             F(NT,J) = 2.0*S(NF,J) * DTEMP / DXWIND + CFAC
           ELSE
             F(NT,J) = 2.0*S(NF,J) * CEND + CFAC
           ENDIF
C
        ELSE
           F(NT, J) = S(NT,J) - CKBSEC(NTEMP, X(J), XGIVEN, TGIVEN)
        ENDIF
C
C       eigenvalue equation, H
        F(NH,J) = S(NH,JP1) - S(NH,J)
C
1000  CONTINUE
C
C     oxidizer boundary (X=L)
      CALL CKRHOY (P, S(NT,JJ), S(NY,JJ), ICKWRK, RCKWRK, RHOJ)
C
      JJM1 = JJ - 1
      RHOV = RHOJ * VOXID
      DO 1200 K = 1, KK
         N = NYS + K
         F(N, JJ) = RHOV * (YOXID(K) - S(N,JJ)) - RHOJ*YV(K,JJM1)
1200  CONTINUE
C
      IF (LENRGY) THEN
         F(NT, JJ) = S(NT,JJ) - TOXID
      ELSE
         F(NT, JJ) = S(NT,JJ) - CKBSEC(NTEMP, X(JJ), XGIVEN, TGIVEN)
      ENDIF
C
      F(NH,JJ) = 0.5*RHOV - S(NF,JJ)
      F(NG,JJ) = S(NG,JJ) + AOXID*RHOJ
      F(NF,JJ) = 0.5*(S(NG,JJ)+S(NG,JJM1)) - (S(NF,JJ)-S(NF,JJM1))/DXP
C
      IF (.NOT. LTIME) RETURN
C
C     add the time step, if needed
      DO 2500 J = 2, JJ-1
         CALL CKRHOY (P, S(NT,J), S(NY,J), ICKWRK, RCKWRK, RHO)
C
         DO 2200 K = 1, KK
            N = NYS + K
            DYDT = (S(N,J) - SN(N,J)) / DT
            F(N,J) = F(N,J) + RHO*DYDT
2200     CONTINUE

         DGDT = (S(NG,J) - SN(NG,J)) / DT
         F(NG,J) = F(NG,J) + DGDT
C
         IF (LENRGY) THEN
            DTDT = (S(NT,J) - SN(NT,J)) / DT
            F(NT,J) = F(NT,J) + RHO*DTDT
         ENDIF
2500  CONTINUE
C
C     end of FUNCTION FUN
      RETURN
      END
C
      SUBROUTINE INTPL8 (F1, F2, X1, X2, MVAR, N1, N2)
C
C  START PROLOGUE
C       interpolate to get f2(x2) from f1(x1)
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION F1(MVAR,*), F2(MVAR,*), X1(*), X2(*)
      DO 100 J = 2, N2-1
        XVAL = X2(J)
        DO 50 I = 2, N1
          IF(XVAL.LE.X1(I)) THEN
            DEL = (XVAL-X1(I-1))/(X1(I)-X1(I-1))
            DO 40 L = 1, MVAR
              F2(L,J)=F1(L,I-1)+(F1(L,I)-F1(L,I-1))*DEL
40          CONTINUE
            GOTO 80
          ENDIF
50      CONTINUE
        WRITE(7,*) ' *** STOP...INTERPOLATION ERROR.'
80      CONTINUE
100   CONTINUE
C
C     endpoints..
      DO 110 L = 1, MVAR
        F2(L,1)  = F1(L,1)
        F2(L,N2) = F1(L,N1)
110   CONTINUE
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION SLINE (WMIX, XCEN, FLOORL, FLOORR, X)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION SLINE (WMIX, XCEN, FLOORL, FLOORR, X)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C  START PROLOGUE
C  END PROLOGUE
C
      IF (X .LE. (XCEN-WMIX/2.) ) THEN
         SLINE = FLOORL
         RETURN
      ELSEIF ( X .LT. (XCEN+WMIX/2.) ) THEN
         SLINE = (FLOORR-FLOORL)/WMIX * (X-XCEN+WMIX*0.5) + FLOORL
         RETURN
      ELSE
         SLINE = FLOORR
         RETURN
      ENDIF
C
C     end of FUNCTION SLINE
      END
C
      SUBROUTINE MDIFV  (KK, JJ, NATJ, LMULTI, LTDIF, X, S, WT, YAV,
     1                   XMF, XMFP, P, D, TDR, ICKWRK, RCKWRK, YV, DKJ)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NG, NF, NH, NYS, NY
      DIMENSION X(*), S(NATJ,*), WT(*), YAV(*), XMF(*), XMFP(*),
     1          D(KK,*), TDR(KK,*), DKJ(KK,KK,*), YV(KK,*),
     2          ICKWRK(*), RCKWRK(*)
      LOGICAL LTDIF, LDIRCT, LMULTI
C
      LDIRCT = .TRUE.
      CALL CKYTX (S(NY,1), ICKWRK, RCKWRK, XMFP)
C
C     loop over all mesh points, computing the diffusion velocity
C     at the midpoints; the indexing is such that YV(K,J) is the
C     diffusion velocity of species K midway betwee nodes J and J+1
      DO 1000 J = 1, JJ-1
C
         JP1 = J + 1
         TMF   = S(NT, J)
         TMFP  = S(NT, JP1)
         TAV = 0.5 * (TMF + TMFP)
         CALL CKAVG (KK, S(NYS+1,J), S(NYS+1,JP1), YAV)
         CALL TWCOPY (KK, XMFP, XMF)
         CALL CKMMWY (YAV, ICKWRK, RCKWRK, WTM)
         CALL CKYTX  (S(NY,JP1), ICKWRK, RCKWRK, XMFP)
C
         XDIF = X(JP1) - X(J)
         IF (LMULTI. AND. LDIRCT) THEN
C           evaluate the multicomponent diffusion velocity directly,
C           rather than use the mixture-averaged form for D(K,J)
C
            WTM2 = WTM**2
            DO 475 K = 1, KK
               SUM = 0.0
               DO 450 L = 1, KK
                  SUM = SUM + WT(L) * DKJ(K,L,J) *
     1                     (XMFP(L)-XMF(L)) / XDIF
  450          CONTINUE
               YV(K,J) = (WT(K)/WTM2) * SUM
  475       CONTINUE
         ELSE
C           use mixture-averaged form for Fickian diffusion,
C           whether we are using the multicomponent formalism
C           or mixture-averaged
C
            CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOAV)
            DO 500 K = 1, KK
               YV(K,J) = - D(K,J) * (WT(K)/WTM) *
     1                     (XMFP(K)-XMF(K)) / XDIF
500         CONTINUE
         ENDIF
C
         IF (LTDIF) THEN
C           add the thermal diffusion, if requested
C
            TFAC = (TMFP-TMF) / (XDIF * TAV * RHOAV)
            TRHO = TAV * RHOAV
            DO 600 K = 1, KK
               YV(K,J) = YV(K,J) - TFAC * TDR(K,J)
600         CONTINUE
C
         ENDIF
C
C        add correction velocity to diffusion
         SUM = 0.
         DO 700 K = 1, KK
           SUM = SUM + YV(K,J)
700      CONTINUE
         VC = - SUM
         DO 800 K = 1, KK
           YV(K,J) = YV(K,J) + YAV(K)*VC
800      CONTINUE
C
1000  CONTINUE
C
C     end of SUBROUTINE MDIFV
      RETURN
      END
C
      SUBROUTINE MTRNPR (KK, JJ, NATJ, LENRGY, LMULTI, LTDIF, P, X, S,
     1                    YAV, ICKWRK, RCKWRK, IMCWRK, RMCWRK, WT,
     2                    XMF, XMFP, XAV, VISC, COND, D, TDR, DKJ)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NG, NF, NH, NYS, NY
      DIMENSION X(JJ), S(NATJ,JJ), YAV(KK), XAV(KK), WT(KK), XMF(KK),
     1          XMFP(KK), ICKWRK(*), RCKWRK(*), IMCWRK(*), RMCWRK(*),
     2          VISC(JJ), COND(JJ), DKJ(KK,KK,JJ), D(KK,JJ), TDR(KK,JJ)
      LOGICAL LENRGY, LTDIF, LMULTI
C
      DATA EPS/1.0E-30/
C
      IF (LMULTI) THEN
C        multicomponent formulas
C
         CALL CKYTX (S(NY,1), ICKWRK, RCKWRK, XMFP)
         DO 200 J = 1, JJ-1
            JP1 = J + 1
            TAV = 0.5 * (S(NT,J) + S(NT,JP1))
C
C           dimensional temperature at the grid points
            CALL CKAVG (KK, S(NYS+1,J), S(NYS+1,JP1), YAV)
            CALL TWCOPY (KK, XMFP, XMF)
            CALL CKYTX  (YAV, ICKWRK, RCKWRK, XAV)
            CALL CKYTX  (S(NY,JP1), ICKWRK, RCKWRK, XMFP)
            CALL CKMMWX (XAV, ICKWRK, RCKWRK, WTMAV)
            CALL MCMDIF (P, TAV, XAV, KK, IMCWRK, RMCWRK, DKJ(1,1,J))
            CALL MCAVIS (TAV, XAV, RMCWRK, VISC(J))
C
            XDIF = X(JP1) - X(J)
            DO 75 K = 1, KK
               SUMN = 0.0
               DO 50 L = 1, KK
                  N = NYS + L
                  SUMN = SUMN + DKJ(K,L,J) * (S(N,JP1)-S(N,J))/XDIF
   50          CONTINUE
               N = NYS + K
               DENOM = - (S(N,JP1) - S(N,J)) / XDIF
               D(K,J) = (SUMN + EPS) / ( WTMAV * (DENOM + EPS))
   75       CONTINUE
C
C           determine the mixture conductivity and thermal
C           diffusion coefficient at J
            IF (LENRGY .OR. LTDIF)
     1         CALL MCMCDT (P, TAV, XAV, IMCWRK, RMCWRK,
     2                      ICKWRK, RCKWRK, TDR(1,J), COND(J))
C
200      CONTINUE
         RETURN
      ENDIF
C
C     mixture-averaged formulas
      DO 400 J = 1, JJ-1
         JP1 = J + 1
         TAV = 0.5 * (S(NT,J) + S(NT,JP1))
C
C        dimensional temperature at the grid points
         CALL CKAVG (KK, S(NYS+1,J), S(NYS+1,JP1), YAV)
         CALL CKYTX  (YAV, ICKWRK, RCKWRK, XAV)
         CALL MCADIF (P, TAV, XAV, RMCWRK, D(1,J))
         CALL MCAVIS (TAV, XAV, RMCWRK, VISC(J))
C
C        determine the mixture conductivity at J
         IF (LENRGY) CALL MCACON (TAV, XAV, RMCWRK, COND(J))
         IF (LTDIF) CALL MCATDR (TAV, XAV, IMCWRK, RMCWRK, TDR(1,J))
C
400   CONTINUE
C
C     end of SUBROUTINE MTRNPR
      RETURN
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION PLATOW (W, XC, FL, FP, FR, X0, X, XL)
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION PLATOW (W, XC, FL, FP, FR, X0, X, XL)
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C  START PROLOGUE
C  END PROLOGUE
C
      WMID = 0.5*W
      IF (X .LE. XC-WMID) THEN
         PLATOW = FL + (FP-FL)*(X-X0)/(XC-WMID-X0)
C         RETURN
      ELSEIF ( X .LT. XC+WMID) THEN
         PLATOW = FP
C         RETURN
      ELSEIF ( X .LE. XL) THEN
         PLATOW = FP + (FR-FP)*(X-XC+WMID)/(XL-XC+WMID)
C         RETURN
      ELSE
         WRITE (6,*) 'STOP, X NOT FOUND IN -PLATOW-'
         STOP
      ENDIF
C
C     end of FUNCTION PLATOW
      RETURN
      END
C
      SUBROUTINE POINT (LINKCK, LINKMC, NMAX, LOUT, LROUT, KK, II, MM,
     1                  NATJ, LENIWK, LENRWK, LENCWK, LENITW, LENTWP,
     2                  LENICK, LENRCK, LTOT, ITOT, NTOT, ICTOT, I, R,
     4                  C, ASIZE, GROUPA, GROUPB, VMAX)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION I(LENIWK), R(LENRWK)
      CHARACTER C(LENCWK)*(*), ILINK*16
      INTEGER ASIZE, GROUPA, GROUPB, VMAX
C
      COMMON /FLFLFL/ NCKW, NMCW, NYOX, NYFL, NWT,  NFL,  NPD, NOX,
     1                NSCH, NX,   NVIS, NCON, NTGV, NXGV, ND,  NDKJ,
     2                NTDR, NYV,  NABV, NBLW, NBUF, NTWP, NS,  NSN,
     3                NFF,  NFN,  NDS,  NSSV, NA,
     4                NRKFT, NRKRT, NUSV, NRSV,
     5                ICKW, IMCW, ITPW, IKS, ICC, IIP, LAC, LMK
      COMMON /LOCS/ NT, NG, NF, NH, NYS, NY
C
C     set the pointers into the solution vector
      NT  = 1
      NG  = 2
      NF  = 3
      NH  = 4
      NYS = 4
      NY  = 5
C
C     call for CHEMKINK work array lengths
      CALL CKLEN (LINKCK, LOUT, LENICK, LENRCK, LENCCK, IFLAG)
      IF (IFLAG .NE. 0) THEN
         WRITE (LOUT,*) ' STOP.  BAD RETURN FROM CKLEN!'
         STOP
      ENDIF
      CALL MCLEN (LINKMC, LOUT, LENIMC, LENRMC, IFLAG)
      IF (IFLAG .NE. 0) THEN
         WRITE (LOUT,*) ' STOP. BAD RETURN FROM MCLEN!'
         STOP
      ENDIF

C     real chemkin work space
      NCKW = 1
C     real transport work space
      NMCW = NCKW + LENRCK
      NTOT = NMCW + LENRMC

C     integer chemkin space
      ICKW = 1
C     integer transport space
      IMCW = ICKW + LENICK
      ITOT = IMCW + LENIMC

C     characcter chemkin space
      ICC = 1
      ICTOT = ICC + LENCCK

      IF (ITOT.LT.LENIWK .AND. NTOT.LT.LENRWK .AND. ICTOT.LT.LENCWK)
     1   THEN
         CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT, I, R, C,
     1                IFLAG)
         IF (IFLAG .NE. 0) THEN
           WRITE (LOUT,*) ' STOP.  BAD RETURN FROM CKINIT!'
           STOP
         ENDIF
         ILINK = 'CKLINK           '
         WRITE (LROUT) ILINK
         CALL CKSAVE (LOUT, LROUT, I, R, C)
         CALL CKINDX (I, R, MM, KK, II, NFIT)
C
         CALL MCINIT (LINKMC, LOUT, LENIMC, LENRMC, I(IMCW), R(NMCW),
     1                IFLAG)
         IF (IFLAG .NE. 0) THEN
           WRITE (LOUT,*) ' STOP.  BAD RETURN FROM MCINIT!'
           STOP
         ENDIF
         ILINK = 'MCLINK          '
         WRITE (LROUT) ILINK
         CALL MCSAVE (LOUT, LROUT, I(IMCW), R(NMCW))
      ENDIF
C
C     number of TWOPNT variables
      NATJ = KK + NYS
      VMAX = GROUPA + NATJ * NMAX + GROUPB
C
C     TWOPNT work space lengths
      LENITW = 3 * NMAX
      LENTWP = 3 * NMAX + 9 * VMAX
C
C     apportion the floating point space
      NYOX = NTOT
      NYFL = NYOX + KK
      NWT  = NYFL + KK
      NFL  = NWT  + KK
      NPD  = NFL  + KK
      NOX  = NPD  + KK
      NSCH = NOX  + KK
      NX   = NSCH + KK*5
      NVIS = NX   + NMAX
      NCON = NVIS + NMAX
      NTGV = NCON + NMAX
      NXGV = NTGV + NMAX
      ND   = NXGV + NMAX
      NDKJ = ND   + KK*NMAX
      NTDR = NDKJ + KK*KK*NMAX
      NYV  = NTDR + KK*NMAX
      NABV = NYV  + KK*NMAX
      NBLW = NABV + GROUPA + NATJ + GROUPB
      NBUF = NBLW + GROUPA + NATJ + GROUPB
      NTWP = NBUF + VMAX
      NS   = NTWP + LENTWP
      NSN  = NS   + VMAX
      NFF  = NSN  + VMAX
      NFN  = NFF  + VMAX
      NDS  = NFN  + VMAX
      NSSV = NDS  + NMAX
      NRKFT = NSSV + NMAX
      NRKRT = NRKFT + II*NMAX
      NUSV = NRKRT + II*NMAX
      NRSV = NUSV + VMAX
      NA   = NRSV + KK*NMAX
      ASIZE = (6*NATJ-1) * VMAX
      NTOT = NA   + ASIZE - 1
C
C     apportion the integer space
      ITPW = ITOT
      IIP  = ITPW  + LENITW
      ITOT  = IIP  + VMAX - 1
C
C     apportion the logical space
      LAC  = 1
      LMK = LAC + NATJ
      LTOT = LMK  + NMAX - 1
C
C     apportion the character space
      IKS = ICTOT
      ICTOT = IKS + KK
C
C     end of SUBROUTINE POINT
      RETURN
      END
C
      SUBROUTINE PRINT (LOUT, KK, JJ, NATJ, P, X, S, XMF, KSYM,
     1                  ICKWRK, RCKWRK, KW, NW)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NG, NF, NH, NYS, NY
      DIMENSION X(*), S(NATJ, *), XMF(*), ICKWRK(*), RCKWRK(*),
     1          KW(*)
      CHARACTER KSYM(KK)*(*)
C
      DATA KPERLN /6/
C
      LS = 1
      K2 = MIN (KPERLN, KK)
C
C     print the flow variables
      WRITE (LOUT, 7070)
      DO 100 J = 1, JJ
         CALL CKRHOY (P, S(NT,J), S(NY,J), ICKWRK, RCKWRK, RHO)
         VEL = 2.0*S(NF,J)/RHO
         WRITE (LOUT, 7020) X(J), S(NF,J), VEL, S(NG,J), S(NH,J),
     1                   S(NT,J), RHO
100   CONTINUE
C
C     print the mole fractions
      IF (NW .EQ. 0) THEN
C
         DO 200 L = LS, KK, KPERLN
            K1   = L
            K2   = L+KPERLN-1
            K2 = MIN (K2, KK)
C
            WRITE (LOUT, 7060) (KSYM(K), K=K1,K2)
            DO 200 J = 1, JJ
               CALL CKYTX (S(NY,J), ICKWRK, RCKWRK, XMF)
               WRITE (LOUT, 7020) X(J), (XMF(K), K=K1,K2)
200      CONTINUE
C
      ELSE
C
         WRITE (LOUT, 7060) (KSYM(KW(K)), K=1,NW)
         DO 210 J = 1, JJ
            CALL CKYTX (S(NY,J), ICKWRK, RCKWRK, XMF)
            WRITE (LOUT, 7020) X(J), (XMF(KW(K)), K=1,NW)
210      CONTINUE
      ENDIF
C
7020  FORMAT(1X, F5.2, 6(1PE11.3))
7060  FORMAT(/2X, 'X(cm)' , 3X, 10(A10, 1X))
7070  FORMAT(/2X, 'X(cm)', 5X, 'F', 7X, 'V(cm/s)', 7X, 'G',
     1           8X, 'H', 10X, 'T(K)', 7X, 'RHO')
C
C     end of SUBROUTINE PRINT
      RETURN
      END
C
      SUBROUTINE RDKEY (KK, NMAX, NATJ, LIN, LOUT, KSYM, PATM, LUSTGV,
     1                  LENRGY, LMULTI, LTDIF, LUMESH, LRSTRT,
     2                  IPROFL, RATEF, N1CALL, LSEN, VOXID,
     3                  VFUEL, AFUEL, AOXID, TFUEL, TOXID, TMAX, P,
     4                  NPTS, X, FUEL, PROD, OXID, XCEN, XEND, WMIX,
     5                  IRETIR, SFLR, NTEMP, XX, TT, WNDFAC, LREGRD,
     6                  JJREGD, PCTADP, RATGTC, KW, NW, SPOS, NDPR,
     7                  DT2, ISTPS2)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NG, NF, NH, NYS, NY
      DIMENSION FUEL(*), OXID(*), PROD(*),
     1          XX(*), TT(*), X(*), VALUE(5), KW(*)
      CHARACTER KEY*4, KSYM(*)*(*), LINE*80, ID*9, CKCHUP*4
      LOGICAL LENRGY, LMULTI, LTDIF, LUMESH, LRSTRT,
     1        LUSTGV, LSEN, NEC(8), NOPT(5), CNTNUD, LFIRST,
     2        IERR, KERR, ERROR, LREGRD
      INTEGER LENGTH, LEVELM, LEVELD, CKSLEN
      EXTERNAL CKCHUP, CKSLEN
C
      PARAMETER (ID = ' RDKEY:  ')
      DATA NEC/8*.FALSE./, NOPT/5*.FALSE./, CNTNUD/.FALSE./
C
C     initialize variables
      KERR = .FALSE.
C
         DO 10 K = 1, KK
            FUEL(K) = 0.
            OXID(K) = 0.
            PROD(K)  = 0.
10       CONTINUE
         DO 11 K = 1, 6
           KW(K) = 0
 11      CONTINUE
C
         NPTS  = 6
         SFLR  = -1.E-4
         SPOS  = 1.E-10
         NTEMP = 0
         NP    = 0
         NW    = 0
         WNDFAC= 1.0
         VOXID = 0.
         N1CALL= 1
         P     = PATM
         TFUEL = 300.
         TOXID = 300.
         TMAX  = 2200.
         AFUEL = 0.
         AOXID = 0.
         JJREGD= 40
         PCTADP= .75
         RATGTC= 1.0
         IPRNT = 1
         RATEF = 1.0
         LFIRST= .TRUE.
         LUMESH= .TRUE.
         LUSTGV= .FALSE.
         LRSTRT= .FALSE.
         LTDIF = .FALSE.
         LMULTI= .FALSE.
         LENRGY= .FALSE.
         LSEN  = .FALSE.
         LREGRD= .FALSE.
         IPROFL= 3
         NDPR = 100

C        set default values for TWOPNT's controls, for two
C        different problems, (1) fixed temperature and
C        (2) energy equation
C
         ERROR = .FALSE.
         CALL TWSETL (ERROR, LOUT, 'ADAPT', .FALSE.)
         CALL TWSETI (ERROR, LOUT, 'LEVELD', 1)
         CALL TWSETI (ERROR, LOUT, 'LEVELM', 1)
         CALL TWSETI (ERROR, LOUT, 'PADD',  10)
         ATOL = 1.0E-9
         CALL TWSETR (ERROR, LOUT, 'SSABS', ATOL)
         CALL TWSETI (ERROR, LOUT, 'SSAGE', 20)
         RTOL = 1.0E-4
         CALL TWSETR (ERROR, LOUT, 'SSREL', RTOL)
         CALL TWSETL (ERROR, LOUT, 'STEADY', .TRUE.)
         CALL TWSETI (ERROR, LOUT, 'STEPS0', 0)
         CALL TWSETI (ERROR, LOUT, 'STEPS1', 100)
         CALL TWSETI (ERROR, LOUT, 'STEPS2', 50)
         STRID0 = 1.0E-6
         CALL TWSETR (ERROR, LOUT, 'STRID0', STRID0)
         ATIM = 1.0E-9
         CALL TWSETR (ERROR, LOUT, 'TDABS', ATIM)
         CALL TWSETI (ERROR, LOUT, 'TDAGE', 20)
         TDEC = 2.2
         CALL TWSETR (ERROR, LOUT, 'TDEC', TDEC)
         RTIM = 1.0E-4
         CALL TWSETR (ERROR, LOUT, 'TDREL', RTIM)
         TINC = 2.0
         CALL TWSETR (ERROR, LOUT, 'TINC', TINC)
         TMAX = 1.0E-4
         CALL TWSETR (ERROR, LOUT, 'TMAX', TMAX)
         TMIN = 1.0E-10
         CALL TWSETR (ERROR, LOUT, 'TMIN', TMIN)
         TOLER0 = 1.0E-8
         CALL TWSETR (ERROR, LOUT, 'TOLER0', TOLER0)
         TOLER1 = 0.1
         CALL TWSETR (ERROR, LOUT, 'TOLER1', TOLER1)
         TOLER2 = 0.5
         CALL TWSETR (ERROR, LOUT, 'TOLER2', TOLER2)
C
         KERR = KERR.OR.ERROR

C
C--------------------------------------------------------------
C
      WRITE (LOUT,'(/1X, A, A/)') ID, 'READING THE KEYWORD INPUT.'
C
90    CONTINUE
C     read the next input line
C
      KEY = ' '
      LINE = ' '
      READ  (LIN,  '(A)') LINE
      CALL EXTENT (LENGTH, LINE)
      LENGTH = MAX(CKSLEN(LINE), 1)
      WRITE (LOUT, 7000) LINE(1 : LENGTH)
      KEY = CKCHUP(LINE(1:4), 4)
      LINE(1:4) = ' '
C
C               IS THIS A KEYWORD COMMENT?
C
      IF (KEY(1:1) .EQ. '.' .OR. KEY(1:1) .EQ. '/'
     1                         .OR. KEY(1:1) .EQ. '!') GO TO 90
C
C-----PROBLEM TYPE KEYWORDS--------------------
C
      IERR = .FALSE.
      ERROR = .FALSE.
      IF (KEY .EQ. 'TGIV') THEN
C        energy equation is not included
         LENRGY   = .FALSE.
C
      ELSEIF (KEY .EQ. 'ENRG') THEN
C        energy equation is included
         LENRGY = .TRUE.
C
      ELSEIF (KEY .EQ. 'GFAC') THEN
C        factor on reaction rates
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RATEF, IERR)
C
      ELSEIF (KEY .EQ. 'LINE') THEN
C        linear profile on initial solution
         IPROFL = 2
C
      ELSEIF (KEY .EQ. 'PLAT') THEN
C        profile on initial solution is plateau
         IPROFL = 3
C
C-----ORIGINAL METHOD OPTIONS KEYWORDS--------------------
C
      ELSEIF (KEY .EQ. 'NOFT') THEN
C        do not do the fixed temperature problem
         N1CALL = 2
C
      ELSEIF (KEY .EQ. 'WDIF') THEN
C        windward differencing
         WNDFAC = 1.0
C
      ELSEIF (KEY .EQ. 'CDIF') THEN
C        central differencing
         WNDFAC = 0.0
C
      ELSEIF (KEY .EQ. 'SFLR') THEN
C        floor value for the species bounds
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SFLR, IERR)
C
      ELSEIF (KEY .EQ. 'SPOS') THEN
C        positive value for negative species
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SPOS, IERR)
C
      ELSEIF (KEY .EQ. 'KOUT') THEN
C        oxidizer
         CALL CKCRAY (LINE, KK, KSYM, LOUT, 6, KW, NFD, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING KEYWRD '//KEY
            IF (NFD .GT. 6) THEN
             WRITE (LOUT,'(A,I4,A)') 'FOUND ', NFD,' SPECIES'
             NW = 6
            ENDIF
         ELSE
           NW = NFD
         ENDIF
C
C///////////////////////////////////////////////////////////////////////
C
C     METHOD OPTIONS KEYWORDS
C
C///////////////////////////////////////////////////////////////////////
C
      ELSEIF (KEY .EQ. 'ATOL') THEN
C        absolute Newton iteration convergence criteria
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'SSABS', VALUE(1))
C
      ELSEIF (KEY .EQ. 'RTOL') THEN
C        relative Newton iteration convergence criteria
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'SSREL', VALUE(1))
C
      ELSEIF (KEY .EQ. 'ATIM') THEN
C        absolute Newton convergence criteria for timesteps
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TDABS', VALUE(1))
C
      ELSEIF (KEY .EQ. 'RTIM') THEN
C        relative Newton convergence criterial for timesteps
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TDREL', VALUE(1))
C
      ELSEIF (KEY .EQ. 'TIME') THEN
C        time step starting procedure
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         ISTPS1 = INT(VALUE(1))
         CALL TWSETI (ERROR, LOUT, 'STEPS1', ISTPS1)
         CALL TWSETR (ERROR, LOUT, 'STRID0', VALUE(2))
C
      ELSEIF (KEY .EQ. 'ISTP') THEN
C        number of initial time steps before Newton
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         ISTPS0 = INT(VALUE(1))
         CALL TWSETI (ERROR, LOUT, 'STEPS0', ISTPS0)
C
      ELSEIF (KEY .EQ. 'NDPR') THEN
C        interval for time solution output
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         NDPR = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'TRAN') THEN
C        do not try Newton for steady state
         CALL TWSETL (ERROR, LOUT, 'STEADY', .FALSE.)
C
      ELSEIF (KEY .EQ. 'IRET') THEN
C        retirement age of old time step (default=50)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         ISTPS2 = INT(VALUE(1))
         CALL TWSETI (ERROR, LOUT, 'STEPS2', ISTPS2)
C
      ELSEIF (KEY .EQ. 'TIM2') THEN
C        time stepping, after adding the energy equation
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         DT2 = VALUE(2)
         ISTPS2 = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'UFAC') THEN
C        timestep increase when timestep does not change solution
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TINC', VALUE(1))
C
      ELSEIF (KEY .EQ. 'DFAC') THEN
C        timestep decrease when Newton fails convergence on timestep
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TDEC', VALUE(1))
C
      ELSEIF (KEY .EQ. 'DTMN') THEN
C        minimum timestep
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TMIN', VALUE(1))
C
      ELSEIF (KEY .EQ. 'DTMX') THEN
C        maximum timestep
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TMAX', VALUE(1))
C
      ELSEIF (KEY .EQ. 'TJAC') THEN
C        retirement age of Jacobian during time stepping
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         ITAGE = INT(VALUE(1))
         CALL TWSETI (ERROR, LOUT, 'TDAGE', ITAGE)
C
      ELSEIF (KEY .EQ. 'NADP') THEN
C        number of mesh points to add at once
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IADD = INT(VALUE(1))
         CALL TWSETI (ERROR, LOUT, 'PADD', IADD)

C///////////////////////////////////////////////////////////////////////
C
C     GRID PARAMETER KEYWORDS
C
C///////////////////////////////////////////////////////////////////////
C
      ELSEIF (KEY .EQ. 'NPTS') THEN
C        number of initial mesh points (overwritten 'GRID' input)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         NPTS = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'GRID') THEN
C        initial mesh
         LUMESH = .FALSE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IERR = IERR .OR. (NP+1 .GT. NMAX)
         IF (.NOT. IERR) THEN
            NP = NP + 1
            X(NP) = VALUE(1)
         ENDIF
C
      ELSEIF (KEY .EQ. 'GRAD') THEN
C        gradient mesh adaption parameter
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TOLER1', VALUE(1))
C
      ELSEIF (KEY .EQ. 'CURV') THEN
C        curvature mesh adaption parameter
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         CALL TWSETR (ERROR, LOUT, 'TOLER2', VALUE(1))
C
      ELSEIF (KEY .EQ. 'XCEN') THEN
C        center of mixing region
         NOPT(2) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XCEN, IERR)
C
      ELSEIF (KEY .EQ. 'XEND') THEN
C        distance at which end boundary condition is applied
         NEC(1) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XEND, IERR)
C
      ELSEIF (KEY .EQ. 'WMIX') THEN
C        width of mixing zone
         NOPT(1) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, WMIX, IERR)
C
      ELSEIF (KEY .EQ. 'WDIF') THEN
C        windward differencing
         WNDFAC = 1.0
C
      ELSEIF (KEY .EQ. 'CDIF') THEN
C        central differencing
         WNDFAC = 0.0
C
C-----FLAME DEFINITION KEYWORDS--------------------
C
      ELSEIF (KEY .EQ. 'VOXI') THEN
C        reactant inlet velocities
         NEC(4)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VOXID, IERR)
         VOXID = - VOXID
C
      ELSEIF (KEY .EQ. 'VFUE') THEN
         NEC(5)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VFUEL, IERR)
C
      ELSEIF (KEY .EQ. 'AOXI') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, AOXID, IERR)
C
      ELSEIF (KEY .EQ. 'AFUE') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, AFUEL, IERR)
C
      ELSEIF (KEY .EQ. 'TOXI') THEN
C        inlet temperature
         NEC(3)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TOXID, IERR)
C
      ELSEIF (KEY .EQ. 'TFUE') THEN
         NEC(2)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TFUEL, IERR)
C
      ELSEIF (KEY .EQ. 'TMAX') THEN
C        maximum temperature
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TMAX, IERR)
C
      ELSEIF (KEY .EQ. 'PRES') THEN
C        pressure
         CALL CKXNUM (LINE, 1, LOUT, NVAL, P, IERR)
         P = P*PATM
C
      ELSEIF (KEY .EQ. 'FUEL') THEN
C        fuel
         IF (LFIRST) THEN
            LFIRST = .FALSE.
            DO 1100 K = 1, KK
               FUEL(K) = 0.
1100        CONTINUE
         ENDIF
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING KEYWRD '//KEY
         ELSE
            FUEL(KSPEC) = VALUE(1)
         ENDIF
C
      ELSEIF (KEY .EQ. 'PROD') THEN
C        product
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING KEYWRD '//KEY
         ELSE
            PROD(KSPEC) = VALUE(1)
         ENDIF
C
      ELSEIF (KEY .EQ. 'OXID') THEN
C        oxidizer
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING KEYWRD '//KEY
         ELSE
            OXID(KSPEC) = VALUE(1)
         ENDIF
C
      ELSEIF (KEY .EQ. 'TEMP') THEN
C        read specified temperature profile (X,T) pairs
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         IF (NTEMP+1 .GT. NMAX) THEN
            WRITE (LOUT, *)
     1       ' ERROR... THE PROBLEM IS ONLY DIMENSIONED FOR ', I4,
     2       ' (X,T) PAIRS'
            IERR = .TRUE.
         ENDIF
         IF (.NOT. IERR) THEN
            NTEMP = NTEMP+1
            XX(NTEMP) = VALUE(1)
            TT(NTEMP) = VALUE(2)
         ENDIF
C
      ELSEIF (KEY .EQ. 'USTG') THEN
C        on a restart use given temperature profile,
C        not the one on the restart file
         LUSTGV = .TRUE.
C
C-----TRANSPORT OPTIONS KEYWORDS--------------------
C
      ELSEIF (KEY .EQ. 'MULT') THEN
C        multicomponent formulas used
         LMULTI = .TRUE.
C
      ELSEIF (KEY .EQ. 'MIX') THEN
C        mixture-averaged formulas used
         LMULTI = .FALSE.
C
      ELSEIF (KEY .EQ. 'TDIF') THEN
C        thermal diffusion included
         LTDIF    = .TRUE.
C
C-----SENSITIVITY KEYWORDS--------------------
C
      ELSEIF (KEY .EQ. 'ASEN') THEN
C        all reaction sensitivity
         LSEN = .TRUE.
C
C-----PRINTING AND RESTARTING KEYWORDS--------------------
C
      ELSEIF (KEY .EQ. 'PRNT') THEN
C        print control
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IPRNT    = INT(VALUE(1))
         IF (IPRNT .LT. 10) THEN
            LEVELM = MIN (2, IPRNT)
            LEVELD = LEVELM
         ELSE
            LEVELM = IPRNT / 10
            LEVELD = IPRNT - 10 * LEVELM
            LEVELM = MAX( LEVELM, LEVELD )
         ENDIF
         CALL TWSETI (ERROR, LOUT, 'LEVELD', LEVELD)
         CALL TWSETI (ERROR, LOUT, 'LEVELM', LEVELM)
C
      ELSEIF (KEY .EQ. 'RSTR') THEN
C        restart check
         LRSTRT = .TRUE.
C
      ELSEIF (KEY .EQ. 'JJRG') THEN
C        number of mesh points in the regrid
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         JJREGD = INT(VALUE(1))
         LREGRD = .TRUE.
C
      ELSEIF (KEY .EQ. 'PCAD') THEN
C        percentage of regrid points dedicated to adaption
         CALL CKXNUM (LINE, 1, LOUT, NVAL, PCTADP, IERR)
C
      ELSEIF (KEY .EQ. 'RGTC') THEN
C        ratio of gradient regrid points to curvature points
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RATGTC, IERR)
C
      ELSEIF (KEY .EQ. 'END ') THEN
C        last line
         GO TO 6000
C
      ELSE
C
C-----THE KEYWORDS AFTER HERE ARE NOT OPERATIONAL------------
C
C-----END OF KEYWORDS--------------------
C
C        to get here, an invalid keyword was read
         WRITE (LOUT, *) ' ERROR...ILLEGAL KEYWORD'
         IERR = .TRUE.
      ENDIF
C
      KERR = KERR.OR.IERR.OR.ERROR
      GO TO 90
C
C     check the reactant and product sums
6000  CONTINUE
      IF (KERR) STOP
C
      SUMF       = 0.
      SUMO       = 0.
      SUMP       = 0.
      DO 6100 K = 1, KK
         SUMF = SUMF+FUEL(K)
         SUMO = SUMO+OXID(K)
         SUMP = SUMP+PROD(K)
6100  CONTINUE
C
C     normalize reactant and product fractions
      DO 6200 K = 1, KK
         FUEL(K) = FUEL(K)/SUMF
         OXID(K) = OXID(K)/SUMO
         PROD(K) = PROD(K)/SUMP
6200  CONTINUE
C
      IF (ABS(SUMF-1.0) .GT. 1.E-6) WRITE (LOUT, *)
     1                ' CAUTION...FUEL FRACTIONS SUM TO ', SUMF
      IF (ABS(SUMO-1.0) .GT. 1.E-6) WRITE (LOUT, *)
     1                ' CAUTION...OXIDIZER FRACTIONS SUM TO ', SUMO
      IF (ABS(SUMP-1.0) .GT. 1.E-6) WRITE (LOUT, *)
     1                ' CAUTION...PRODUCT FRACTIONS SUM TO ', SUMP
C
C     check for necessary input
      IF (.NOT. NEC(1)) THEN
         WRITE (LOUT, *) ' ERROR..."XEND" NOT SPECIFIED'
         STOP
      ENDIF
C
      IF (.NOT. NEC(2) ) THEN
         WRITE (LOUT, *) ' ERROR..."TFUEL" NOT GIVEN'
         STOP
      ENDIF
      IF (.NOT. NEC(3) ) THEN
         WRITE (LOUT, *) ' ERROR..."TOXID" NOT GIVEN'
         STOP
      ENDIF
C
      IF (.NOT. NEC(4) ) THEN
         WRITE (LOUT, *) ' ERROR..."VOXID" NOT GIVEN'
         STOP
      ENDIF
C
      IF (.NOT. NEC(5) ) THEN
         WRITE (LOUT, *) ' ERROR..."VFUEL" NOT GIVEN'
         STOP
      ENDIF
C
C     make sure the initial grid points are in order
      IF ((.NOT. CNTNUD) .AND. (.NOT. LUMESH)) THEN
         NPTS = NP
         DO 5550 N = 2, NPTS
            IF (X(N-1) .GE. X(N)) THEN
              WRITE (LOUT, *)
     1              ' ERROR...INITIAL GRID IS OUT OF ORDER'
              STOP
            ENDIF
5550     CONTINUE
      ENDIF
C
      IF (NTEMP.GT.0) THEN
C         make sure the (X,T) pairs are in order
C
          DO 5500 N = 2, NTEMP
            IF (XX(N-1) .GE. XX(N)) THEN
            WRITE (LOUT, *)
     1              ' ERROR...SPECIFIED TEMPERATURES ARE OUT OF ORDER'
            STOP
            ENDIF
5500      CONTINUE
C
C         make sure the given temperatures span the XEND-XSTR domain
          IF (.NOT.LRSTRT .OR. .NOT.CNTNUD .OR. LUSTGV) THEN
            IF (XX(1).GT.0.0 .OR. XX(NTEMP).LT.XEND) THEN
               WRITE (LOUT, *)
     1      ' ERROR...GIVEN TEMPERATURE PROFILE DOES NOT SPAN XEND-XSTR'
               STOP
            ENDIF
          ENDIF
      ENDIF
C
C     set optional input if needed
      IF (.NOT. NOPT(1)) WMIX = (XEND-0.0)*0.50
      IF (.NOT. NOPT(2)) XCEN = (XEND-0.0)*0.35
C
99999 CONTINUE
C
C     formats
7000  FORMAT(10X, A)
C
C     end of SUBROUTINE RDKEY
      RETURN
      END


      SUBROUTINE REGRID (SNEW, SOLD, NVAR, JJNEW, JJOLD, XNEW, XOLD,
     1                   WORK, P, R, N )
C
C  START PROLOGUE
C  INTERPOLATE THE SOLUTION ONTO AN EQUIDISTRIBUTING MESH
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION SNEW(NVAR,*), SOLD(NVAR,*), XNEW(*), XOLD(*), WORK(*)
C
C     compute weight coefficients
      R0 = 1. - P
      R1 = P *R /(R+1.)
      R2 = P - R1
      TV1 = 0.
      DO 10 I = 2, JJOLD
        TV1 = TV1 + ABS( SOLD(N,I) - SOLD(N,I-1) )
10    CONTINUE
      TV2 = 0.
      DO 20 I = 2, JJOLD-1
        TV2 = TV2 + ABS( (SOLD(N,I+1)-SOLD(N,I))/(XOLD(I+1)-XOLD(I))
     1                 - (SOLD(N,I)-SOLD(N,I-1))/(XOLD(I)-XOLD(I-1)) )
20    CONTINUE
      XLEN = XOLD(JJOLD) - XOLD(1)
      B1 = R1 * XLEN /((1.-P) * TV1)
      B2 = R2 * XLEN /((1.-P) * TV2)
C
C     compute partial sums of weight function
C
      WORK(1) = 0.
      DO 50 I = 2, JJOLD
         DX = XOLD(I) -XOLD(I-1)
         WORK(I) = DX + B1*ABS( SOLD(N,I) - SOLD(N,I-1) ) + WORK(I-1)
     1           + B2*ABS( (SOLD(N,I+1)-SOLD(N,I))/(XOLD(I+1)-XOLD(I))
     2                   - (SOLD(N,I)-SOLD(N,I-1))/(XOLD(I)-XOLD(I-1)) )
50    CONTINUE
      DO 65 I = 2, JJOLD
         WORK(I) = WORK(I)/WORK(JJOLD)
65    CONTINUE
C
C      interpolate onto uniform eta grid to find new x
C
      XNEW(1) = XOLD(1)
      XNEW(JJNEW) = XOLD(JJOLD)
      ISTART  = 2
      DETA = 1./FLOAT(JJNEW-1)
      DO 80 J = 2, JJNEW-1
         ETAJ = (J-1)*DETA
         DO 70 I = ISTART, JJOLD
            IF (ETAJ .LE. WORK(I)) THEN
               DEL = (ETAJ-WORK(I-1))/(WORK(I)-WORK(I-1))
               XNEW(J) = XOLD(I-1)+(XOLD(I)-XOLD(I-1))*DEL
               GO TO 80
            ELSE
               ISTART = I
            ENDIF
70       CONTINUE
         WRITE (6, *) ' *** VALUE OF ETA NOT FOUND ***'
80    CONTINUE
C
C        interpolate solution...
C
      CALL INTPL8 (SOLD, SNEW, XOLD, XNEW, NVAR, JJOLD, JJNEW)
C
C     end of SUBROUTINE REGRID
      RETURN
      END


      SUBROUTINE RESTRT (KK, NMAX, NATJ, JJ, LOUT, LUMESH, LUSTGV,
     1                   FUEL, PROD, OXID,
     2                   NTEMP, XGIVEN, TGIVEN, XCEN, XEND,
     3                   WMIX, ICKWRK, RCKWRK, Y, SI, YOXID, YFUEL,
     4                   X, S)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NG, NF, NH, NYS, NY

      DIMENSION FUEL(*), OXID(*), PROD(*),
     1          ICKWRK(*), RCKWRK(*), Y(*), YOXID(KK), YFUEL(KK),
     2          XGIVEN(*), TGIVEN(*), X(*), S(NATJ, *)
C
      LOGICAL LUMESH, LUSTGV
C
C     initialize mass flux fractions to zero
      DO 100 K = 1, KK
         YOXID(K) = 0.0
         YFUEL(K) = 0.0
100   CONTINUE
C
C     if a new XEND > X(JJ), then add a point at JJ + 1;
C     if a new XEND <=X(JJ), then reduce JJ and set X(JJ)=XEND
      IF (XEND .GT. (X(JJ)+1.E-4) ) THEN
         JJ = JJ + 1
         IF (JJ .GT. NMAX) THEN
            WRITE (LOUT, *) ' ERROR...NEW XEND NEEDS TOO MANY POINTS'
            STOP
         ENDIF
         X(JJ) = XEND
         CALL TWCOPY (NATJ, S(1,JJ-1), S(1,JJ))
      ELSE
         DO 500 J = 1, JJ
            IF (XEND .LE. X(J)) THEN
               X(J) = XEND
               GO TO 550
            ENDIF
500      CONTINUE
550      CONTINUE
         JJ = J
      ENDIF
C
C     set the mass flux fraction boundary conditions
      CALL CKXTY (FUEL, ICKWRK, RCKWRK, YFUEL)
      CALL CKXTY (OXID, ICKWRK, RCKWRK, YOXID)
C
C     set XGIVEN and TGIVEN to the old solution
      IF (.NOT. LUSTGV) THEN
         NTEMP = JJ
         DO 1250 N = 1, NTEMP
            XGIVEN(N) = X(N)
            TGIVEN(N) = S(NT,N)
1250     CONTINUE
      ENDIF
C
C     end of SUBROUTINE RESTRT
      RETURN
      END
C
      SUBROUTINE START (KK, NMAX, NATJ, JJ, LOUT, LUMESH, LENRGY,
     1                  IPROFL, TMAX, TFUEL, TOXID, VFUEL, VOXID,
     2                  P, FUEL, PROD, OXID, NTEMP, XGIVEN,
     3                  TGIVEN, XCEN, XEND, WMIX, ICKWRK, RCKWRK,
     4                  Y, YOXID, YFUEL, X, S)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NG, NF, NH, NYS, NY
      DIMENSION FUEL(*), OXID(*), PROD(*),
     1          ICKWRK(*), RCKWRK(*), Y(*), YOXID(KK), YFUEL(KK),
     2          XGIVEN(*), TGIVEN(*), X(*), S(NATJ, *)
C
      LOGICAL LUMESH, LENRGY
      EXTERNAL CKBSEC
C
C     initialize mass flux fractions and mass fractions to zero
      DO 100 K = 1, KK
         YOXID(K) = 0.0
         YFUEL(K) = 0.0
         DO 100 J = 1, NMAX
            S(NYS+K, J) = 0.0
100   CONTINUE
C
C     set the mass flux fraction boundary conditions
      CALL CKXTY (OXID, ICKWRK, RCKWRK, YOXID)
      CALL CKXTY (FUEL, ICKWRK, RCKWRK, YFUEL)
      CALL CKXTY (PROD, ICKWRK, RCKWRK, Y )
C
C     compute stagnation plane location
      CALL CKRHOY (P, TFUEL, YFUEL, ICKWRK, RCKWRK, RHOF)
      CALL CKRHOY (P, TOXID, YOXID, ICKWRK, RCKWRK, RHOO)
      IF (XCEN .EQ. 0.) THEN
         XCEN = XEND / 2.
         IF (LUMESH) THEN
           WMIX = 2.*MIN( XEND-XCEN, XCEN )
         ELSE
           WMIX = 2.*MIN( X(JJ)-XCEN, XCEN-X(1) )
         ENDIF
      ENDIF
C
      IF (LUMESH) THEN
C        set uniform X mesh coordinates
C
         DX = (XEND-0.0) / FLOAT(JJ-1)
         DO 200 J = 1, JJ
            X(J) = 0.0 + DX*FLOAT(J-1)
200      CONTINUE
      ENDIF
C
      IF (IPROFL .EQ. 2) THEN
C        set linear profiles
C
         DO 900 K = 1,KK
            N = NYS + K
            DO 900 J = 1,JJ
               S(N,J) = SLINE (WMIX, XCEN, YFUEL(K), YOXID(K), X(J))
900      CONTINUE
         IF (LENRGY) THEN
           DO 950 J = 1, JJ
             S(NT,J) = TFUEL
950        CONTINUE
         ENDIF
C
      ELSEIF (IPROFL .EQ. 3) THEN
C        set plateau profiles
C
         DO 960 K = 1,KK
            N = NYS + K
            DO 960 J = 1,JJ
               S(N,J) = PLATOW (WMIX, XCEN, YFUEL(K), Y(K), YOXID(K),
     1                      X(1), X(J), X(JJ))
960      CONTINUE
         IF (LENRGY) THEN
           DO 970 J = 1,JJ
             S(NT,J) = PLATOW (WMIX, XCEN, TFUEL, TMAX, TOXID,
     1                      X(1), X(J), X(JJ))
970        CONTINUE
         ENDIF
C
      ENDIF
C
C      CHECK IF TGIVEN (NOT ENERGY), SET TGIVEN OR S(NT,J)
C
      IF (LENRGY) THEN
C       set the 'given' temperatures to the profile
        NTEMP = JJ
        DO 1770 N = 1, NTEMP
          XGIVEN(N) = X(N)
          TGIVEN(N) = S(NT,N)
1770    CONTINUE
C
      ELSEIF ( (.NOT. LENRGY) .AND. (.NOT. LUMESH) ) THEN
C        set profile to given temperatures
         DO 1780 J = 1, JJ
           S(NT,J) = CKBSEC(NTEMP, X(J), XGIVEN, TGIVEN)
1780     CONTINUE
      ENDIF
C
      FFUEL = 0.5 * VFUEL * RHOF
      FOXID = 0.5 * VOXID * RHOO
      DFDX = (FOXID - FFUEL)/(X(JJ) - X(1))

      DO 1800 J = 1, JJ
C        set a linear F profile
         S(NF,J) = DFDX* (X(J)-X(1)) + FFUEL
C        set G as constant (to math linear profile in F)
         S(NG,J) = DFDX
C        set the pressure gradient eigenvalue
         S(NH,J) = -100.
1800  CONTINUE
C
C     end of SUBROUTINE START
      RETURN
      END
      SUBROUTINE OPPOS (S, IDIM, SPOS)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER(I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER(I-N)
C*****END precision > single
C
      DIMENSION S(IDIM)
      DO 10 I = 1, IDIM
         S(I) = MAX (S(I), SPOS)
   10 CONTINUE
C
C     end of SUBROUTINE OPPOS
      RETURN
      END
C
      SUBROUTINE OPREAD (LRIN, LOUT, GROUPA, GROUPB, II,
     1                   JJ, NATJ, X, S, ERROR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INTEGER GROUPA, GROUPB, VRBLS
      DIMENSION X(*), S(*)
      LOGICAL ERROR
      CHARACTER*16 ICHR, ICKLNK, IMCLNK, ISOLUT, ISENSI
      PARAMETER (ICKLNK='CKLINK', IMCLNK='MCLINK', ISOLUT='SOLUTION',
     1           ISENSI='SENSITIVITY')
C
C///  READ THE RESTART FILE.
C
3010  CONTINUE
      READ (LRIN) ICHR
C
C     don't read link arrays with solution
      IF (ICHR .EQ. ICKLNK) THEN
         DO 3020 L=1,4
            READ (LRIN)
3020     CONTINUE
C
      ELSEIF (ICHR .EQ. IMCLNK) THEN
         DO 3030 L=1,3
            READ (LRIN)
3030     CONTINUE
C
      ELSEIF (ICHR .EQ. ISOLUT) THEN
         READ (LRIN) NNNN, JJ, PDUM, AFDUM, AODUM, VWDUM, VWLDUM
         READ (LRIN) (X(J), J=1,JJ)
         VRBLS = GROUPA + NNNN*JJ + GROUPB
         READ (LRIN) (S(J), J=1, VRBLS)
         ERROR = NNNN .NE. NATJ
         IF (ERROR) GO TO 9301
         GO TO 3040
C
      ELSEIF (ICHR .EQ. ISENSI) THEN
         DO 3045 I = 1, II
            READ (LRIN)
3045     CONTINUE
C
      ELSE
         ERROR = .TRUE.
         IF (ERROR) GO TO 9302
      ENDIF
      GO TO 3010
3040  CONTINUE
      GO TO 99999
C
9301  WRITE (LOUT, 99301) ID
      GO TO 99999
C
9302  WRITE (LOUT, 99302) ID
      GO TO 99999
C
9303  WRITE (LOUT, 99303) ID
      GO TO 99999
C
99301 FORMAT
     +  (/1X, A, 'ERROR. INCOMPATIBLE RESTART FILE.')
C
99302 FORMAT
     +  (/1X, A, 'ERROR. NOT A SOLUTION ON RESTART FILE.')
C
99303 FORMAT
     +  (/1X, A, 'ERROR. ERROR READING SOLUTION FILE.')
C
C///  EXIT.
C
99999 CONTINUE
C
C     end of SUBROUTINE OPREAD
      RETURN
      END
