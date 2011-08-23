C  CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      SUBROUTINE PSENK (LIN, LOUT, LSAVE, LTXT, KMAX, JMAX, IMAX, 
     1                  LENIWK, IWORK, LENRWK, RWORK, LENCWK, CWORK,
     2                  LI, I, LR, R, LC, C, ISYM)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      DIMENSION IWORK(LENIWK), RWORK(LENRWK), I(LI), R(LR)
      CHARACTER*16 CWORK(LENCWK), C(LC)
      CHARACTER*32 ISYM(IMAX)
      LOGICAL KERR
      KERR = .FALSE.
C
      NS    = 1
      NTIM  = NS    + KMAX * JMAX
      NCONT = NTIM  + JMAX
      NSUMRT= NCONT + JMAX * IMAX
      NYMF  = NSUMRT+ JMAX
      NXMF  = NYMF  + JMAX * KMAX
      NXT   = NXMF  + JMAX * KMAX
      NTEMP = NXT   + KMAX
      NPRES = NTEMP + JMAX
      NWT   = NPRES + JMAX
      NQ    = NWT   + KMAX
      NWDOT = NQ    + JMAX
      NCIK  = NWDOT + KMAX
      NSENSJ= NCIK  + IMAX
      NSEN  = NSENSJ+ IMAX * KMAX
      NSSUM = NSEN  + IMAX * JMAX
      NWTM  = NSSUM + IMAX * JMAX
      NTOT  = NWTM  + JMAX - 1
      IF (NTOT .GT. LR) THEN
         WRITE (LOUT, *)
     1   'SENKOUT ERROR: RPAR needs to be at least ', NTOT
         KERR = .TRUE.
      ENDIF
C
      ISPEC = 1
      ICONT = ISPEC + KMAX
      ITOT  = ICONT + IMAX - 1 
      IF (ITOT .GT. LI) THEN
         WRITE (LOUT, *)
     1   'SENKOUT ERROR: IPAR needs to be at least ', ITOT
         KERR = .TRUE.
      ENDIF
      IF (KMAX .GT. LC) THEN
         WRITE (LOUT, *)
     1   'SENKOUT ERROR: CPAR needs to be at least ', KMAX
         KERR = .TRUE.
      ENDIF
      IF (KERR) RETURN
C
      CALL POSTSEN (LIN, LOUT, LSAVE, LTXT, KMAX, JMAX, IMAX,
     1              LENIWK, IWORK, LENRWK, RWORK, LENCWK, CWORK,
     2              R(NS), R(NTIM), R(NCONT), R(NSUMRT), R(NYMF),
     3              R(NXMF), R(NXT), R(NTEMP), R(NPRES), R(NWT),
     4              R(NQ), R(NWDOT), R(NCIK), R(NSENSJ), R(NSEN),
     5              R(NSSUM), R(NWTM),
     6              I(ISPEC), I(ICONT), C, ISYM)
      RETURN
      END
C
      SUBROUTINE POSTSEN (LIN, LOUT, LSAVE, LTXT, KMAX, JMAX, IMAX,
     1           LENIWK, IWORK, LENRWK, RWORK, LENCWK, CWORK,
     2           S, TIM, CONT, SUMRT, YMF, XMF, XT, T, P, WT,
     3           Q, WDOT, CIK, SENSJ, SEN, SSUM, WTM, ISPEC,
     4           ICONT, KSYM, ISYM)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single

C      TEXT POST-PROCESSING FOR SENKIN CALCULATIONS.
C      KEYWORD INPUT TO SELECT TEMPERATURE, PRESSURE, MOLE FRACTIONS,
C      NORMALIZED SENSITIVITY COEFFICIENTS, AND RATE-OF-PRODUCTION.
C      VERSION DATE: 04/15/97
C      AUTHOR:  ANDREW LUTZ
C               SANDIA NATIONAL LABORATORIES
C               LIVERMORE, CA 94551-0969
C               (510) 294-2761
C               aelutz@ca.sandia.gov
C
C      V.2.1: fix out-of-range subscript for -n32 compiler
C      V.2.2: fix bug#111: remove unused variables CKLSCH, REDKEY,
C             and ILINK.
C      V.2.3: fix bug#149: Add initialization of JNUM and JJ to 0
C             before Loop 80 in RDSAVE for case when sensitivities
C             were not stored. (E. Meeks)

C     DIMENSIONS:
C      KMAX  = NO. SPECIES + 1   = KK + 1
C      IMAX  = NO. REACTIONS     = II (IF NO SENSITIVITY COEFFICIENTS)
C            = NO. REACTION + 1  = II + 1 (IF SENSITIVITY )
C      JMAX  = MAX. NO. TIME DATASETS TO BE USED

C     ARRAY STRUCTURE:
C      T(J)     = TEMPERATURE AT J-TH TIME
C      S(K,J) = MASS FRACTION OF K-TH SPECIES AT J-TH TIME
C      SEN(L,J) = SENSITIVITY COEFFICIENT AT J-TH TIME
C                 W.R.T. THE L-TH REACTION
C      OTHER ARRAYS ARE WORK ARRAYS

C      PARAMETER ( JMAX=300, IMAX=450, KMAX=100, LENIWK=9500,
C     1            LENRWK=9500, LENCWK=100, LENSYM=16, LENIHL=32,
C     2            LIN=5, LOUT=6, LSAVE=55, LTXT=7)
      PARAMETER (LENSYM=16, LENIHL=32)

      DIMENSION  IWORK(LENIWK), RWORK(LENRWK), S(KMAX,JMAX), TIM(JMAX),
     1           CONT(JMAX,IMAX), SUMRT(JMAX), YMF(JMAX,KMAX),
     2           XMF(JMAX,KMAX), XT(KMAX), T(JMAX), P(JMAX), 
     3           WT(KMAX), Q(JMAX), WDOT(KMAX), CIK(IMAX),
     4           SENSJ(KMAX,IMAX), SEN(JMAX,IMAX), SSUM(JMAX,IMAX),
     5           ISPEC(KMAX), WTM(JMAX), ICONT(IMAX)

      CHARACTER*24 ISTR, FILNAM, CKCHLO

      LOGICAL LWANT(4), LMOLE, LINIT, KERR, IERR

      CHARACTER KSYM(KMAX)*(LENSYM), ISYM(IMAX)*(LENIHL),
C     1          CWORK(LENCWK)*(LENSYM), LIST(9)*60, SLIST(3)*(76),
     1          CWORK(LENCWK)*(LENSYM), SLIST(3)*(76),
     2          NFORM*16, RFORM*16, SFORM*16
      EXTERNAL CKCHLO

      DATA JJ/0/, LWANT, LMOLE, LINIT/6*.FALSE./,
     1     TMAX, TMIN /2* 0./, SLIST /3*' '/

C     NUMBER OF COLUMNS PER LINE AND FORMATS FOR TABLES
      DATA NPERLS/10/, NFORM/'(10(1X,1PE11.4))'/, SFORM/'(10(1X,A10))'/,
     1     NPERLR/8/, RFORM/'(8(1X,A15))'/

      NDIM = KMAX

      WRITE (LOUT,'(//,A,/,A,/)') ' SENKIN Post-Processor',
     1 '         (CHEMKIN-II) Version 2.1, 97/04/15'
      WRITE (LOUT,'(A)') '  Enter name for binary data file.'
      ISTR = ' '
      READ (LIN,'(A)') ISTR
      IND = INDEX(ISTR,'!')
      IF (IND .GT. 0) ISTR(IND:) = ' '
      FILNAM = CKCHLO(ISTR, 24)
      OPEN (LSAVE, FILE=FILNAM, STATUS='OLD',
     1      FORM='UNFORMATTED', ERR=10)
      WRITE (LOUT,'(A)') '  Reading file '//FILNAM
      GO TO 15
10    CONTINUE
      WRITE (LOUT,'(A)') '  Stop, cannot open file '//FILNAM
      RETURN
15    CONTINUE

      WRITE (LOUT,'(A)') '  Enter name for text output file.'
      ISTR = ' '
      READ (LIN,'(A)') ISTR
      IND = INDEX(ISTR,'!')
      IF (IND .GT. 0) ISTR(IND:) = ' '
      FILNAM = CKCHLO(ISTR, 24)
      OPEN (LTXT, STATUS='UNKNOWN', FORM='FORMATTED', 
     1      FILE=FILNAM, ERR=20)
      WRITE (LOUT,'(A)') ' Writing to file '//FILNAM
      GO TO 25
20    CONTINUE
        WRITE (LOUT,'(A)') '  Stop, cannot open file '//FILNAM
        RETURN
25    CONTINUE

C     READ KEYWORD INPUT

      CALL REDKEY ( LIN, LOUT, LWANT, SLIST, 
     1              TMIN, TMAX, LMOLE, KERR )
      IF (KERR) RETURN
      DELTIM = (TMAX - TMIN) / FLOAT(JMAX)

C     READ THE BINARY FILE (FIRST TIME)

      KNUM = -1
      CALL RDSAVE (LSAVE, LINIT, LOUT, IWORK, RWORK, KNUM, KMAX,
     1             JMAX, IMAX, WT, TIM, P, T, KK, II, JJ, TMIN, TMAX,
     2             DELTIM, S, SEN, SENSJ, SSUM, XT, YMF, XMF, WTM,
     3             CWORK, ISYM, PATM, LENIWK, LENRWK, LENCWK, KSYM,
     4             KERR)
      IF (KERR) RETURN

C     WRITE TEMPERATURE VS TIME

      IF ( LWANT(1) ) THEN
          IF (SLIST(1) .NE. ' ') THEN
            CALL CKCRAY (SLIST(1), KK, KSYM, LOUT, NDIM, 
     1                   ISPEC, NWANT, IERR)
            IF (IERR) THEN
              WRITE (LOUT,'(A)') 'ERROR FINDING SPECIES IN LIST!'
              NWANT = 0
            ENDIF
          ELSE
            NWANT = KK
            DO 100 K = 1, KK
100           ISPEC(K) = K
          ENDIF
          CALL TABLE (LOUT, LTXT, TIM, T, P, XMF, YMF, JJ, JMAX, KK, 
     1                WTM, ISPEC, NWANT, KSYM, LMOLE, NPERLS, NFORM,
     2                SFORM )
      ENDIF

C     WRITE RATES VS TIME

      IF (LWANT(4)) THEN
        CALL CKCRAY (SLIST(2), KK, KSYM, LOUT, NDIM, 
     1                   ISPEC, NWANT, IERR)
        IF (IERR) THEN
           WRITE (LOUT,'(A)') 'ERROR FINDING SPECIES IN LIST!'
        ELSE   
          CALL RATES (LOUT, LTXT, NWANT, ISPEC, KMAX, S, XT, T, P,
     1                IWORK, RWORK, CWORK, Q, CIK, SUMRT, JMAX,
     2                CONT, II, KSYM, ISYM, TIM, JJ, ICONT, 
     3                NPERLR, NFORM, RFORM )
        ENDIF
      ENDIF
C
C      SENSITIVITY OF TEMPERATURE
C
      IF ( LWANT(2) ) THEN

         KNUM = 0
         CALL RDSAVE (LSAVE, LINIT, LOUT, IWORK, RWORK, KNUM, KMAX,
     1                JMAX, IMAX, WT, TIM, P, T, KK, II, JJ, TMIN, TMAX,
     2                DELTIM, S, SEN, SENSJ, SSUM, XT, YMF, XMF, WTM,
     3                CWORK, ISYM, PATM, LENIWK, LENRWK, LENCWK, KSYM,
     4                KERR)
         IF (KERR) RETURN

         CALL TSENS (LOUT, LTXT, T, TIM, SEN, CONT, JMAX, JJ,
     1               ISYM, II, ICONT, NPERLR, NFORM, RFORM)

      ENDIF

C     SENSITIVITY OF SPECIES

      IF ( LWANT(3) ) THEN

         CALL CKCRAY (SLIST(3), KK, KSYM, LOUT, NDIM, 
     1                   ISPEC, NWANT, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)') 'ERROR FINDING SPECIES IN LIST!'
            NWANT = 0
         ENDIF

         DO 900 N = 1, NWANT

           KNUM = ISPEC(N)

           CALL RDSAVE (LSAVE, LINIT, LOUT, IWORK, RWORK, KNUM, KMAX,
     1                  JMAX, IMAX, WT, TIM, P, T, KK, II, JJ, 
     2                  TMIN, TMAX,
     2                  DELTIM, S, SEN, SENSJ, SSUM, XT, YMF, XMF, WTM,
     3                  CWORK, ISYM, PATM, LENIWK, LENRWK, LENCWK, KSYM,
     4                  KERR)
          IF (KERR) RETURN

           CALL KSENS (LMOLE, LOUT, LTXT, KNUM, YMF, XMF, SEN, SSUM,
     1                 WT, CONT, JMAX, JJ, KK, II, XT, IWORK, RWORK,
     2                 CWORK, ISYM, KSYM, TIM, ICONT, NPERLR,
     3                  NFORM, RFORM)

900       CONTINUE

      ENDIF

      RETURN
      END

C--------------------------------------------------------------------

      SUBROUTINE TABLE (LOUT, LTXT, TIM, T, P, XMF, YMF, JJ, JMAX, KK, 
     1                  WTM, KWANT, NW, KSYM, MOLE, NPERL, NFORM,
     2                  SFORM)

C     PRINT TABLE: T, P, X_K IN TIME INTERVAL FOR SOLUTION

C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double

      DIMENSION TIM(JMAX), T(JMAX), P(JMAX), XMF(JMAX,KK), 
     1          KWANT(NW), WTM(JMAX), YMF(JMAX, KK)

      LOGICAL MOLE

      CHARACTER KSYM(KK)*16, NFORM*(*), SFORM*(*)

C     WRITE FIRST GROUP, INCLUDING T,P

      K1 = 1
1000  CONTINUE
        K2 = MIN( K1 + NPERL - 4, NW )
        IF (MOLE) THEN
          WRITE (LTXT,'(1X,A)') 'MOL FRACTIONS:'
          WRITE (LTXT, SFORM) 't(s)', 'T(K)', 'P(atm)',
     1      (KSYM(KWANT(K)), K = K1, K2)
          DO 1100 J = 1, JJ          
            WRITE (LTXT, NFORM)
     1      TIM(J), T(J), P(J), (XMF(J,KWANT(K)), K = K1, K2)
1100      CONTINUE
        ELSE
          WRITE (LTXT,'(1X,A)') 'MASS FRACTIONS:'
          WRITE (LTXT, SFORM) 't(s)', 'T(K)', 'P(atm)',
     1      (KSYM(KWANT(K)), K = K1, K2)
          DO 1200 J = 1, JJ          
            WRITE (LTXT, NFORM)
     1      TIM(J), T(J), P(J), (YMF(J,KWANT(K)), K = K1, K2)
1200      CONTINUE
        ENDIF

      IF (K2 .EQ. NW) RETURN

C     LOOP FOR OTHER GROUPS OF SPECIES


1300  CONTINUE
        K1 = K2 + 1
        K2 = MIN( K1 + NPERL - 2, NW )
        IF (MOLE) THEN
          WRITE (LTXT,'(1X,A)') 'MOL FRACTIONS:'
          WRITE (LTXT, SFORM) 't(s)',
     1      (KSYM(KWANT(K)), K = K1, K2)
          DO 1400 J = 1, JJ          
            WRITE (LTXT, NFORM)
     1      TIM(J), (XMF(J,KWANT(K)), K = K1, K2)
1400      CONTINUE
        ELSE
          WRITE (LTXT,'(1X,A)') 'MASS FRACTIONS:'
          WRITE (LTXT, SFORM) 't(s)',
     1      (KSYM(KWANT(K)), K = K1, K2)
          DO 1500 J = 1, JJ          
            WRITE (LTXT, NFORM)
     1      TIM(J), (YMF(J,KWANT(K)), K = K1, K2)
1500      CONTINUE
        ENDIF

        IF (K2 .LT. NW) GO TO 1300

      RETURN         
      END

C--------------------------------------------------------------------

      SUBROUTINE MINMAX (V, N, VMIN, VMAX, IPK)

C      Returns min and max values of vector, and index
C      of the peak absolute value.

C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double

      DIMENSION V(N)

      IPK = 0
      VMAX = V(1)
      VMIN = V(1)
      VPK  = ABS(VMAX)
      DO 100 J = 2, N
         VMIN = MIN (V(J), VMIN)
         VMAX = MAX (V(J), VMAX)
         AMIN = ABS(VMIN)
         AMAX = ABS(VMAX)
         IF (AMIN .GT. VPK) THEN
            VPK = AMIN
            IPK = J
         ELSEIF (AMAX .GT. VPK) THEN
            VPK = AMAX
            IPK = J
         ENDIF
  100 CONTINUE
      RETURN
      END

C---------------------------------------------------------------------

      SUBROUTINE REDKEY (LIN, LOUT, LWANT,
     1                   SLIST, TMIN, TMAX, MOLE, KERR )

C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double

      CHARACTER KEYWRD*4, SLIST(*)*(*)
      CHARACTER ISTR*76, LINE*76, CKCHUP*76
      EXTERNAL CKCHUP
      LOGICAL LWANT(4), KERR, IERR, MOLE

      KERR = .FALSE.
      KRWANT = 0
      KXWANT = 0
      KSWANT = 0
C
C           ISSUE A PROMPT
C
      WRITE (LOUT,'(//2X,A,/,2X,A)') 'Available keywords:',
     1' TMIN, TMAX, TEMP, MOLS, MASS, SPEC, TSEN, XSEN, RATE, END'

C......
C
C         READ NEXT INPUT LINE
C
90    CONTINUE
      WRITE (LOUT,'(/2X,A)') 'Enter keyword:'
      READ  (LIN,  '(A76)') ISTR
      WRITE (LOUT, '(10X,A)') ISTR
C
C               IS THIS A KEYWORD COMMENT?
C
      IF (ISTR(1:1) .EQ. '.' .OR. ISTR(1:1) .EQ. '/' .OR.
     1    ISTR(1:1) .EQ. '!' ) GO TO 90
C
      IND = INDEX(ISTR,'!')
      IF (IND .GT. 0) ISTR(IND:) = ' '
      LINE = CKCHUP(ISTR, 76)
      KEYWRD = LINE(1:4)
      LINE(1:4) = ' '
C
      IF (KEYWRD .EQ. 'TEMP') THEN
         LWANT(1) = .TRUE.
      ELSEIF (KEYWRD .EQ. 'MOLS') THEN
         MOLE = .TRUE.
      ELSEIF (KEYWRD .EQ. 'MASS') THEN
         MOLE = .FALSE.

      ELSEIF (KEYWRD .EQ. 'SPEC') THEN
         LWANT(1) = .TRUE.
         SLIST(1) = LINE

      ELSE IF (KEYWRD .EQ. 'TMIN') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TMIN, IERR)
         KERR = KERR.OR.IERR

      ELSE IF (KEYWRD .EQ. 'TMAX') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TMAX, IERR)
         KERR = KERR.OR.IERR

      ELSE IF (KEYWRD .EQ. 'TSEN') THEN
         LWANT(2) = .TRUE.

      ELSE IF (KEYWRD .EQ. 'XSEN') THEN
         LWANT(3) = .TRUE.
         SLIST(3) = LINE

      ELSE IF (KEYWRD .EQ. 'RATE') THEN
         LWANT(4) = .TRUE.
         SLIST(2) = LINE

      ELSE IF (KEYWRD .EQ. 'END') THEN
         GO TO 6000
      ELSE
C
C..........END OF KEYWORDS
C
C        TO GET HERE, AN INVALID KEYWORD WAS READ
C
         WRITE (LOUT,*) ' Error, illegal keyword.'
      ENDIF
C
C        GO BACK UP AND READ THE NEXT LINE
C
      GO TO 90
C
C         DONE READING CARDS
C
6000  CONTINUE
C
C        STOP IF ERRORS ENCOUNTERED
C
      IF (KERR) RETURN
  555 FORMAT (' Sorry...sensitivity coefficients are not available.')
C
      RETURN
      END

C-------------------------------------------------------------------

      SUBROUTINE RATES (LOUT, LTXT, NWANT, KWANT, IS, S, XT, T, P,
     1                  IWORK, RWORK, CWORK, Q, CIK, SUMRT, IC, CONT,
     2                  II, KSYM, ISYM, TIME, JJ, ICONT, NPERL, NFORM,
     3                  RFORM)

C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double

      DIMENSION S(IS,JJ), IWORK(1), RWORK(1), XT(IS), Q(II), CIK(II),
     1          CONT(IC,II), SUMRT(II), T(JJ), TIME(JJ),  P(JJ),
     2          KWANT(NWANT), ICONT(II)

      CHARACTER KSYM(*)*16, CWORK(*)*16, ISYM(II)*32,
     1          NFORM*16, RFORM*16

C     LOOP OVER SPECIES WANTED

      DO 5000 N = 1, NWANT

         KNUM = KWANT(N)
         WRITE (LOUT,'(A,A)') ' Printing RATE for species ',KSYM(KNUM)

C        CONTRIBUTIONS OF REACTIONS
         DO 310 J = 1, JJ
            CALL CKYTX (S(1,J), IWORK, RWORK, XT)
            CALL CKQXP  (P(J), T(J), XT, IWORK, RWORK, Q)
            CALL CKCONT (KNUM, Q, IWORK, RWORK, CIK)
            CALL CKRHOX (P(J), T(J), XT, IWORK, RWORK, RHO)
            SUMRT(J) = 0.0
            DO 310 I=1,II
               CONT(J,I) = CIK(I)/RHO
               SUMRT(J)  = SUMRT(J) + CONT(J,I)
  310    CONTINUE

C        MAX AND MIN RATES
         YMAX = -1.E30
         YMIN = +1.E30
         DO 325 I = 1, II
            CALL MINMAX (CONT(1,I), JJ, YMN, YMX, IPK)
            YMIN = MIN (YMIN,YMN)
            YMAX = MAX (YMAX,YMX)
  325    CONTINUE
         CALL MINMAX (SUMRT, JJ, YMN, YMX, IPK)
         YMIN = MIN(YMIN, YMN)
         YMAX = MAX(YMAX, YMX)

         WRITE (LTXT,'(2X,A,A)') KSYM(KNUM), 'RATE (moles/gm-sec)'

C        COUNT AND INDEX THE CONTRIBUTING REACTIONS
         NCONT = 0
         DO 330 I = 1, II
            CALL MINMAX (CONT(1,I), JJ, CMIN, CMAX, IPK)
            IF (ABS(CMAX) .GT. .1*ABS(YMAX) .OR.
     1             ABS(CMIN) .GT. .1*ABS(YMIN)) THEN
                NCONT = NCONT + 1
                ICONT(NCONT) = I
            ENDIF
330      CONTINUE

C        PRINT TOTAL PRODUCTION AND LARGEST CONTRIBUTORS

        I2 = 0
1000    CONTINUE
          I1 = I2 + 1
          I2 = MIN( I1 + NPERL - 3, NCONT )
          WRITE (LTXT, RFORM) 't(s) ', 'Total ',
     1           (ISYM(ICONT(I)), I = I1, I2)
          DO 1100 J = 1, JJ          
            WRITE (LTXT, NFORM)
     1      TIME(J), SUMRT(J), (CONT(J,ICONT(I)), I = I1, I2)
1100      CONTINUE

        IF (I2 .EQ. NCONT) GO TO 5000

1200    CONTINUE
          I1 = I2 + 1
          I2 = MIN( I1 + NPERL - 2, NCONT )
          WRITE (LTXT, RFORM) 't(s) ',
     1           (ISYM(ICONT(I)), I = I1, I2)
          DO 1300 J = 1, JJ          
            WRITE (LTXT, NFORM)
     1      TIME(J), (CONT(J,ICONT(I)), I = I1, I2)
1300      CONTINUE
        IF (I2 .LT. NCONT) GO TO 1200

5000  CONTINUE
      RETURN
      END

C----------------------------------------------------------------------
      SUBROUTINE RDSAVE (LSAVE, LINIT, LOUT, ICKWRK, RCKWRK, KNUM, KMAX,
     1                   JMAX, IMAX, WT, TIM, P, T, KK, II, JJ, TMIN,
     3                   TMAX, DELTIM, S, SEN, SENSJ, SSUM, XT,
     2                   YMF, XMF, WTM, CCKWRK, IREAC, PATM,
     3                   LENICK, LENRCK, LENCCK, KSYM, KERR)

C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double

      DIMENSION ICKWRK(LENICK), RCKWRK(LENRCK), S(KMAX,JMAX), 
     1          SENSJ(KMAX,IMAX), SEN(JMAX,IMAX),
     1          SSUM(JMAX,IMAX), XT(KMAX), YMF(JMAX,KMAX),
     2          XMF(JMAX,KMAX), WT(KMAX), P(JMAX), T(JMAX), TIM(JMAX), 
     3          WTM(JMAX)

      LOGICAL LSENS, LINIT, KERR, IERR

      CHARACTER*16 ICHR, CSENK, VERS, PREC
      CHARACTER CCKWRK(LENCCK)*16, IREAC(IMAX)*32, KSYM(KMAX)*16

      DATA CSENK/'SENKIN SOLUTION'/

      LSENS = .FALSE.

      WRITE(LOUT,'(/A/)') '  Reading data.'

C     READ SOLUTION FROM UNIT LSAVE

      REWIND (LSAVE)

10    CONTINUE
      READ (LSAVE) ICHR

      IF (ICHR .EQ. 'CKLINK') THEN
C
C        Initialize Chemkin Library common block
C
         CALL CKPNT (LSAVE, LOUT, NPOINT, VERS, PREC, LENI, LENR, LENC,
     1               IERR)
         KERR = KERR.OR.IERR
C
C*****double precision
         IF (INDEX(PREC,'DOUB') .LE. 0) THEN
C*****END double precision
C*****single precision
C         IF (INDEX(PREC,'SING') .LE. 0) THEN
C*****END single precision
C
            WRITE (LOUT, '(//2X,A,A,A//)') ' Error...', ICHR,
     2      ' precision incompatible with post-processor'
            KERR = .TRUE.
            RETURN
        ENDIF
C
C       Initialize Chemkin work arrays
C
        IF (LENI.LE.LENICK .AND. LENR.LE.LENRCK .AND. LENC.LE.LENCCK)
     1     THEN
           READ (LSAVE, ERR=200, END=200) (ICKWRK(L), L = 1, LENI)
           READ (LSAVE, ERR=200, END=200) (RCKWRK(L), L = 1, LENR)
           READ (LSAVE, ERR=200, END=200) (CCKWRK(L), L = 1, LENC)
           CALL CKINDX (ICKWRK, RCKWRK, MM, KK, II, NFIT)
           CALL CKRP   (ICKWRK, RCKWRK, RU, RUC, PATM)
           CALL CKWT   (ICKWRK, RCKWRK, WT)

           IF (KK .LE. KMAX) THEN
              CALL CKSYMS (CCKWRK, LOUT, KSYM, IERR)
           ELSE
              WRITE (LOUT, *) ' Error...species dimension too small,',
     1                        ' KMAX should be at least ', KK
              KERR = .TRUE.
           ENDIF

           IF (II .LE. IMAX) THEN
              DO 20 I = 1, II
                 CALL CKSYMR (I, LOUT, ICKWRK, RCKWRK, CCKWRK, LT,
     1                        IREAC(I), IERR)               
                 IF (IERR) WRITE (LOUT,*) ' ERROR IN REACTION SYMBOL',I
   20         CONTINUE
           ELSE
              WRITE (LOUT, *) ' Error...reaction dimension too small,',
     1                        ' IMAX should be at least ',II
              KERR = .TRUE.
           ENDIF

           GO TO 10

        ELSE
           READ (LSAVE, ERR=200, END=200) (IDUM, L = 1, LENI)
           READ (LSAVE, ERR=200, END=200) (RDUM, L = 1, LENR)
           READ (LSAVE, ERR=200, END=200) (CDUM, L = 1, LENC)
           WRITE (LOUT, 700)
     1        ICHR, LENI, LENICK, LENR, LENRCK, LENC, LENCCK
 700  FORMAT(' Error...not enough ',A,
     1                      ' work space provided:',
     2        /,'                PROVIDED        REQUIRED',
     3        /,' INTEGER  ', 2I15,
     4        /,' REAL     ', 2I15,
     5        /,' CHARACTER', 2I15,/)
           KERR = .TRUE.
        ENDIF       

      ELSEIF (ICHR .EQ. CSENK) THEN

        READ (LSAVE) LSENS
        READ (LSAVE) NATJ, KKS, IIS

C       CHECK IF SENSITIVITY SOLUTION REQUESTED

        IF ( (KNUM .GE. 0) .AND. (.NOT. LSENS) ) THEN
           WRITE (LOUT,'(//5X,A//)') 'SENSITIVITY NOT AVAILABLE!'
           RETURN
        ENDIF
        IF ( (KNUM .EQ. 0) .AND. (NATJ .EQ. KKS) ) THEN
           WRITE (LOUT,'(//5X,A//)') 
     1     'TEMPERATURE SENSITIVITY NOT AVAILABLE!'
           RETURN
        ENDIF


C       LOOP OVER DATASETS

        IF (LSENS) THEN

C        READ SENSITIVITY COEFFICIENTS ALSO

         JJ = 0
         JNUM = 0
         NT = 0
         IF (NATJ .NE. KKS) NT = 1
         DO 50 J=1,JMAX

45          CONTINUE
            READ (LSAVE,END=90) TIM(J),P(J),T(J),(S(K,J),K=1,KK)
            READ (LSAVE)        ((SENSJ(N,I), N=1,NATJ), I=1,II)

            CALL CKYTX (S(1,J), ICKWRK, RCKWRK, XT)

C           GET MEAN WT...
            CALL CKMMWX (XT, ICKWRK, RCKWRK, WTM)

            DO 992 K=1,KK
C              STORE MASS FRACTIONS
               YMF(J,K) = S(K,J)
C              STORE MOLE FRACTIONS
               XMF(J,K) = XT(K)
  992       CONTINUE

            DO 991 I=1,II
               SSUM(J,I) = 0.0
               DO 991 K=1,KK
                  SSUM(J,I) = SSUM(J,I) + SENSJ(K+NT,I)/WT(K)
  991       CONTINUE

            IF (KNUM .GE. 0) THEN
C              STORE SENS COEF
               DO 8049 I=1,II
                  SEN(J,I) = SENSJ(KNUM+NT,I)
 8049          CONTINUE
            ENDIF

C            NDUM = NATJ*II
c            CALL MINMAX (SENSJ, NDUM, SJMIN, SJMAX, IPKS)
c            WRITE (LOUT,*) ' Max sens = ',SJMAX


C          CHECK FOR SAVING THE DATA POINT

            JNUM = JNUM + 1
            IF (DELTIM .GT. 0.) THEN
              IF (TIM(J) .LT. TMIN) GO TO 45
              IF (J .GT. 1) THEN
                  IF (TIM(J)-TIM(J-1) .LT. DELTIM) GOTO 45
              ENDIF
              JJ = JJ + 1
              IF (TIM(J) .GE. TMAX) GO TO 90
            ELSE
              JJ = JJ + 1
            ENDIF             
50       CONTINUE

        ELSE

C        NO SENSITIVITY COEFFICIENTS TO READ

         JJ = 0 
         JNUM = 0
         DO 80 J=1,JMAX

65          CONTINUE
            READ (LSAVE,END = 90) TIM(J), P(J), T(J), (S(K,J),K=1,KK)

            CALL CKYTX (S(1,J), ICKWRK, RCKWRK, XT)
            DO 993 K=1,KK
C              STORE MASS FRACTIONS
               YMF(J,K) = S(K,J)
C              STORE MOLE FRACTIONS
               XMF(J,K) = XT(K)
  993       CONTINUE

C          CHECK FOR SAVING THE DATA POINT

            JNUM = JNUM + 1
            IF (DELTIM .GT. 0.) THEN
              IF (TIM(J) .LT. TMIN) GO TO 65
              IF (J .GT. 1) THEN
                  IF (TIM(J)-TIM(J-1) .LT. DELTIM) GOTO 65
              ENDIF
              JJ = JJ + 1
              IF (TIM(J) .GE. TMAX) GO TO 90
            ELSE
              JJ = JJ + 1
            ENDIF
80       CONTINUE

        ENDIF
      ENDIF

C     IF HERE, DID NOT REACH END OF FILE
      WRITE (LOUT,'(//5X,A,I15/5X,A//)')
     1 'WARNING: DID NOT REACH END OF FILE. JMAX =', JMAX,
     2 '         INCREASE JMAX, OR ENTER A MAXIMUM TIME.'

90    CONTINUE

C     DONE READING DATA

      IF ( (TIM(JJ) .GT. TMAX) .AND. (DELTIM .GT. 0.) ) JJ = JJ - 1
      IF (.NOT. LINIT) THEN
      WRITE (LOUT,'(/,1X,A)') 'Finished reading data.'
      WRITE (LOUT,'(1X,A,I8)') 'Number of time datasets read = ', JNUM
      WRITE (LOUT,'(1X,A,I8)') 'Number kept = ', JJ
      LINIT = .TRUE.
      ENDIF

      RETURN

200   CONTINUE
      WRITE (LOUT, '(//5X,A)') 
     1 'STOP. PROBLEM READING CHEMKIN DATA FROM FILE.'
      KERR = .TRUE.
      RETURN
      END

C-------------------------------------------------------------------

      SUBROUTINE TSENS (LOUT, LTXT, T, TIM, SEN, CONT, JMAX, JJ, 
     1                  ISYM, II, ICONT, NPERL, NFORM, RFORM)

C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double

      DIMENSION SEN(JMAX,II), CONT(JMAX,II), T(JMAX), ICONT(II), 
     1          TIM(JMAX)

      CHARACTER ISYM(II)*32, NFORM*16, RFORM*16


      WRITE (LOUT,*) ' Temperature sensitivities ...'

      DO 8100 J = 1, JJ
        DO 8100 L = 1, II
          CONT(J,L) = SEN(J,L) /( T(J) + 1.E-8)
8100  CONTINUE
      YMAX = -1.E20
      YMIN = -YMAX
      DO 8200 I = 1, II
         CALL MINMAX (CONT(1,I), JJ, YMN, YMX, IPK)
         YMIN = MIN(YMIN, YMN)
         YMAX = MAX(YMAX, YMX)
8200  CONTINUE

      WRITE (LTXT,'(2X,A)') 'NORMALIZED TEMPERATURE SENSITIVITY'

C     COUNT AND INDEX THE CONTRIBUTING REACTIONS
      NCONT = 0
      DO 8210 I = 1, II
         CALL MINMAX (CONT(1,I), JJ, RMIN, RMAX, IPK)
         IF(ABS(RMAX) .GT. 0.1*ABS(YMAX) .OR.
     1      ABS(RMIN) .GT. 0.1*ABS(YMIN)) THEN
                NCONT = NCONT + 1
                ICONT(NCONT) = I
          ENDIF
8210  CONTINUE

C        PRINT MOST IMPORTANT REACTIONS

        I2 = 0
9000    CONTINUE
          I1 = I2 + 1
          I2 = MIN( I1 + NPERL - 1, NCONT )
          WRITE (LTXT, RFORM) 
     1  't(s)', (ISYM(ICONT(I)), I = I1, I2)
          DO 9100 J = 1, JJ          
            WRITE (LTXT, NFORM)
     1      TIM(J), (CONT(J,ICONT(I)), I = I1, I2)
9100      CONTINUE
        IF (I2 .LT. NCONT) GO TO 9000

      RETURN
      END

C----------------------------------------------------------------------

       SUBROUTINE KSENS (LMOLE, LOUT, LTXT, KNUM, YMF, XMF, SEN, SSUM,
     1                   WT, CONT, JMAX, JJ, KK, II, XT, IWORK, RWORK,
     2                   CWORK, ISYM, KSYM, TIM, ICONT, NPERL,
     3                   NFORM, RFORM )

C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double

      DIMENSION YMF(JMAX,KK), XMF(JMAX,KK), CONT(JMAX,II),
     1          SEN(JMAX,II), WT(KK), SSUM(JMAX,II), XT(KK),
     2          IWORK(*), RWORK(*), TIM(JMAX), ICONT(II)

      CHARACTER CWORK(*)*16, KSYM(KK)*16, ISYM(II)*16,
     1          NFORM*16, RFORM*16

      LOGICAL LMOLE

      WRITE (LOUT,'(A,A)') ' Sensitivities for ', KSYM(KNUM)

      IF (.NOT. LMOLE) THEN

         CALL MINMAX (YMF(1,KNUM), JJ, DUM, YMAXK, IPK)
         YMAXK = YMAXK + 1.E-14
         DO 9110 J = 1, JJ
            DO 9110 L = 1, II
               CONT(J,L) = SEN(J,L) /YMAXK
9110     CONTINUE

      ELSE

         CALL MINMAX (XMF(1,KNUM), JJ, DUM, XMAXK, IPK)
         XMAXK = XMAXK + 1.E-14
         DO 9111 J=1,JJ
            DO 992 K=1,KK
               XT(K) = XMF(J,K)
  992       CONTINUE
            CALL CKMMWX (XT, IWORK, RWORK, WTM)

            DO 9111 I=1,II
               CONT(J,I) = (WTM/(XMAXK*WT(KNUM))) *
     1                ( SEN(J,I) - WT(KNUM)*YMF(J,KNUM)*SSUM(J,I) )
 9111    CONTINUE

      ENDIF

      YMAX = -1.E20
      YMIN = -YMAX
      DO 9200 I = 1, II
         CALL MINMAX (CONT(1,I), JJ, YMN, YMX, IPK)
         YMIN = MIN(YMIN, YMN)
         YMAX = MAX(YMAX, YMX)
9200  CONTINUE

      IF (LMOLE) THEN
        WRITE (LTXT,'(2X,A,A)') 
     1  KSYM(KNUM), 'NORMALIZED MOLE FRACTION SENSITIVITY'
      ELSE
        WRITE (LTXT,'(2X,A,A)') 
     1  KSYM(KNUM), 'NORMALIZED MASS FRACTION SENSITIVITY'
      ENDIF

C     COUNT AND INDEX THE CONTRIBUTING REACTIONS
      NCONT = 0
      DO 9230 I = 1, II
         CALL MINMAX (CONT(1,I), JJ, RMIN, RMAX, IPK)
         IF (ABS(RMAX) .GT. 0.1*ABS(YMAX) .OR.
     1       ABS(RMIN) .GT. 0.1*ABS(YMIN)) THEN
                NCONT = NCONT + 1
                ICONT(NCONT) = I
         ENDIF
9230  CONTINUE
      IF (NCONT .EQ. 0) THEN
         WRITE (LOUT,'(//2X,A//)') 'NO CONTRIBUTING REACTIONS?'
         RETURN
      ENDIF

C        PRINT MOST IMPORTANT REACTIONS

        I2 = 0
9500    CONTINUE
          I1 = I2 + 1
          I2 = MIN( I1 + NPERL - 1, NCONT )
          WRITE (LTXT, RFORM) 't(s)', (ISYM(ICONT(I)), I = I1, I2)
          DO 9600 J = 1, JJ          
            WRITE (LTXT, NFORM)
     1      TIM(J), (CONT(J,ICONT(I)), I = I1, I2)
9600      CONTINUE
        IF (I2 .LT. NCONT) GO TO 9500

      RETURN
      END














