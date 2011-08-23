C     CVS $Revision: 1.1.1.1 $ reposited $Date: 2006/05/26 19:09:32 $

* * * * * * * * * * * * * * * * *
* Copyright (C) 1997            *
* Joseph F. Grcar               *
* Sandia National Laboratories  *
* sepp@california.sandia.gov    *
* * * * * * * * * * * * * * * * *

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C     VERSION 1.13 OF DECEMBER 1997
C
C     BLOCK SPARSE MATRIX PACKAGE
C
C     WRITTEN BY:  DR. JOSEPH F. GRCAR
C                  DEPARTMENT 8345
C                  MAIL STOP 9051
C                  SANDIA NATIONAL LABORATORIES
C                  LIVERMORE, CA 94551-0969 USA
C
C                  (510) 294-2662
C
C                  sepp@california.sandia.gov
C                  na.grcar@na-net.ornl.gov
C
C///////////////////////////////////////////////////////////////////////
C
C     CHANGES FROM THE PREVIOUS VERSION:
C
C     1) ADD COPYRIGHT NOTICE.
C
C     2) CORRECT ADDRESS.
C
C///////////////////////////////////////////////////////////////////////

      SUBROUTINE BLCH2
     +  (ERROR, TEXT, PRINT,
     +   ACOL1, ACOLQ1, BLKS1, BMISS, BNULL, CMISS, CNULL, COL1, COLS1,
     +   EXPECT, N, OCCUPY, RETURN, ROW1, ROWS1, X, Y, Y0)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLCH2
C
C>    CALLED BY BLCHE
C
C>    CHECK THE SPARSITY PATTERN.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9, LINE*80, STRING*80, WORD*80
      INTEGER
     +   ACOL1, ACOLQ1, BLOCK, BLKS1, BMISS, BNULL, C1, CHARS, CMISS,
     +   CNULL, COL, COL1, COLS1, FIRST, J, LENGTH, MISS, N, NULL,
     +   PRINT, R1, ROUTE, ROW, ROW1, ROWS1, TEXT, YR1, YR2
      LOGICAL ERROR, EXPECT, HEAD, OCCUPY, RETURN
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   ABSOL, EPS, FACTOR, RELAT, SAVE, X, Y, Y0

      PARAMETER (ID = ' BLCH2:  ')

      PARAMETER (QFRST = 1, QLAST = 2, QFRSTD = 3, QLASTD = 4)
      PARAMETER (QIND = 1, QLOC = 2)

      DIMENSION
     +   ACOL1(2, BLKS1), ACOLQ1(2, COLS1), COL1(2, COLS1),
     +   EXPECT(ROWS1), OCCUPY(ROWS1), ROW1(2, ROWS1), X(N), Y(N), Y0(N)

      SAVE

C///////////////////////////////////////////////////////////////////////
C
C     RETURN FROM REVERSE COMMUNICATION.
C
C///////////////////////////////////////////////////////////////////////

      IF (RETURN) THEN
         RETURN = .FALSE.
         GO TO (2010, 3030) ROUTE
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. BLKS1 .AND. 0 .LT. COLS1 .AND. 0 .LT. N
     +   .AND. 0 .LT. ROWS1)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) INITIALIZE.
C
C///////////////////////////////////////////////////////////////////////

C///  FORM MACHINE EPSILON AND THE ABSOLUTE AND RELATIVE PERTURBATIONS.

      CALL BLEPS (EPS)

      ABSOL = SQRT (EPS)
      RELAT = SQRT (EPS)

C///  EVALUATE THE FUNCTION AT THE UNPERTURBED X.

C     GO TO 2010 WHEN ROUTE = 1
      ROUTE = 1
      RETURN = .TRUE.
      GO TO 99999
2010  CONTINUE

      DO 2020 ROW = 1, N
         Y0(ROW) = Y(ROW)
2020  CONTINUE

C///  INITIALIZE THE COUNTERS.

      BMISS = 0
      BNULL = 0
      CMISS = 0
      CNULL = 0

C///////////////////////////////////////////////////////////////////////
C
C     TOP OF THE LOOP OVER COLUMN BLOCKS.
C
C///////////////////////////////////////////////////////////////////////

      C1 = 1
      HEAD = .FALSE.
2030  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (3) BUILD THE LIST OF OCCUPIED ROW BLOCKS.
C
C///////////////////////////////////////////////////////////////////////

C///  INITIALIZE THE LIST.

      DO 3010 R1 = 1, ROWS1
         OCCUPY(R1) = .FALSE.
3010  CONTINUE

C///  TOP OF THE LOOP OVER COLUMNS IN THE BLOCK.

      ERROR = .NOT. (COL1(QFRST, C1) .LE. COL1(QLAST, C1))
      IF (ERROR) GO TO 9301

      COL = COL1(QFRST, C1)
3020  CONTINUE

C///  EVALUATE THE COLUMN.

      SAVE = X(COL)
      X(COL) = X(COL) + (ABS (X(COL)) * RELAT + ABSOL)

C     GO TO 3030 WHEN ROUTE = 2
      ROUTE = 2
      RETURN = .TRUE.
      GO TO 99999
3030  CONTINUE

      X(COL) = SAVE

      FACTOR = 1.0 / (ABS (X(COL)) * RELAT + ABSOL)
      DO 3040 ROW = 1, N
         Y(ROW) = (Y(ROW) - Y0(ROW)) * FACTOR
3040  CONTINUE

C///  UPDATE THE LIST.

      DO 3060 R1 = 1, ROWS1
         YR1 = ROW1(QFRST, R1)
         YR2 = ROW1(QLAST, R1)
         DO 3050 ROW = YR1, YR2
            OCCUPY(R1) = OCCUPY(R1) .OR. Y(ROW) .NE. 0.0
3050     CONTINUE
3060  CONTINUE

C///  BOTTOM OF THE LOOP OVER COLUMNS IN THE BLOCK.

      COL = COL + 1
      IF (COL .LE. COL1(QLAST, C1)) GO TO 3020

C///////////////////////////////////////////////////////////////////////
C
C     (4) COMPARE LISTS.
C
C///////////////////////////////////////////////////////////////////////

C///  BUILD THE EXPECTED LIST.

      DO 4010 R1 = 1, ROWS1
         EXPECT(R1) = .FALSE.
4010  CONTINUE

      DO 4020 BLOCK = ACOLQ1(QFRST, C1), ACOLQ1(QLAST, C1)
         R1 = ACOL1(QIND, BLOCK)
         EXPECT(R1) = .TRUE.
4020  CONTINUE

C///  COMPARE THE LISTS.

      MISS = 0
      NULL = 0
      DO 4030 R1 = 1, ROWS1
         IF (.NOT. EXPECT(R1) .AND. OCCUPY(R1)) THEN
            MISS = MISS + 1
         ELSE IF (EXPECT(R1) .AND. .NOT. OCCUPY(R1)) THEN
            NULL = NULL + 1
         END IF
4030  CONTINUE

C///  UPDATE THE COUNTERS.

      BMISS = BMISS + MISS
      BNULL = BNULL + NULL

      IF (0 .LT. MISS) CMISS = CMISS + 1
      IF (0 .LT. NULL) CNULL = CNULL + 1

C///////////////////////////////////////////////////////////////////////
C
C     (5) PRINT THE ROW INDICES.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE PRINT BLOCK.

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT .AND.
     +   (0 .LT. MISS .OR. 0 .LT. NULL)) THEN

C///  PRINT THE HEADER.

      IF (.NOT. HEAD) THEN
         HEAD = .TRUE.
         WRITE (TEXT, 10001)
      END IF

C///  PRINT THE ROW INDICES.

      LENGTH = 12
      WRITE (LINE, '(I10, 2X)') C1
      DO 5020 R1 = 1, ROWS1
         IF (EXPECT(R1) .OR. OCCUPY(R1)) THEN
            IF (EXPECT(R1) .AND. .NOT. OCCUPY(R1)) THEN
               WRITE (WORD, '(I15, A1)') R1, '?'
            ELSE IF (.NOT. EXPECT(R1) .AND. OCCUPY(R1)) THEN
               WRITE (WORD, '(I15, A1)') R1, '!'
            ELSE
               WRITE (WORD, '(I16)') R1
            END IF

            DO 5010 J = 16, 1, - 1
               IF (WORD (J : J) .NE. ' ') FIRST = J
5010        CONTINUE
            CHARS = 16 - FIRST + 1

            IF (60 .LT. LENGTH + 1 + CHARS) THEN
               WRITE (TEXT, 10002) LINE (1 : LENGTH)
               LINE = ' '
               LENGTH = 12
            END IF

            STRING = LINE (1 : LENGTH) // ' ' // WORD (FIRST : 16)
            LINE = STRING
            LENGTH = LENGTH + 1 + CHARS
         END IF
5020  CONTINUE
      WRITE (TEXT, 10002) LINE (1 : LENGTH)

C///  BOTTOM OF THE PRINT BLOCK.

      END IF

C///////////////////////////////////////////////////////////////////////
C
C     BOTTOM OF THE LOOP OVER COLUMN BLOCKS.
C
C///////////////////////////////////////////////////////////////////////

      C1 = C1 + 1
      IF (C1 .LE. COLS1) GO TO 2030

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +   (/10X, 'THE FOLLOWING COLUMNS HAVE INCORRECT ROW INDICES.'
     +    /10X, '"!" INDICATES A MISSING ROW, "?" AN UNNECESSARY ROW.'
     +   //10X, '    COLUMN   ROW INDICES'
     +    /)

10002 FORMAT
     +   (10X, A)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, BLKS1, COLS1, N, ROWS1
      GO TO 99999

9301  IF (0 .LT. TEXT) WRITE (TEXT, 99301) ID,
     +   C1, COL1(QFRST, C1), COL1(QLAST, C1)
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  BLKS1'
     +    /10X, I10, '  COLS1'
     +    /10X, I10, '  N',
     +    /10X, I10, '  ROWS1')

99301 FORMAT
     +   (/1X, A9, 'ERROR.  BLJA5 FAILS.')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLCHE
     +  (ERROR, TEXT, PRINT,
     +   ILAST, IMAX, ISIZE, IWORK,
     +   LLAST, LMAX, LSIZE, LWORK,
     +   RLAST, RMAX, RSIZE, RWORK,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   N, BMISS, BNULL, CMISS, CNULL, RETURN, X, Y)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLCHE
C
C>    CALLED BY THE  U S E R
C
C>    CHECK THE MATRIX SPARSITY PATTERN.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   ASIZE2, BLKS1, BLKS2, BMISS, BNULL, BORDER, CMISS, CNULL,
     +   COLS1, COLS2, GLUES, IDATA, IDSIZ, ILAST, IMAX, ISAVE, ISIZE,
     +   IWORK, LLAST, LMAX, LSAVE, LSIZE, N, N1, N2, PRINT, RDSIZ,
     +   RLAST, RMAX, ROWS1, ROWS2, RSAVE, RSIZE, TEXT
      LOGICAL ERROR, LWORK, RETURN
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   RDATA, RWORK, X, Y

      PARAMETER (ID = ' BLCHE:  ')

      DIMENSION
     +   IDATA(IDSIZ), IWORK(ISIZE), LWORK(LSIZE), RDATA(RDSIZ),
     +   RWORK(RSIZE), X(N), Y(N)

      SAVE

C///////////////////////////////////////////////////////////////////////
C
C     RETURN FROM REVERSE COMMUNICATION.
C
C///////////////////////////////////////////////////////////////////////

      IF (RETURN) GO TO 2010

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  SAVE THE LEVELS OF THE WORK SPACES.

      ISAVE = ILAST
      LSAVE = LLAST
      RSAVE = RLAST

C///  PRINT.

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) WRITE (TEXT, 10001) ID

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. IDSIZ .AND. 0 .LT. ISIZE .AND. 0 .LT. LSIZE
     +   .AND. 0 .LT. N .AND. 0 .LT. RDSIZ .AND. 0 .LT. RSIZE)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) CALL BLCH2.
C
C///////////////////////////////////////////////////////////////////////

C///  UNPACK THE DATA SPACES.

      CALL BLDAT
     +  (ERROR, TEXT,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   ASIZE2, BLKS1, BLKS2, BORDER, COLS1, COLS2, GLUES, QA2, QAC1,
     +   QAC2, QACQ1, QACQ2, QAR2, QARQ2, QBMAP, QCMAP, QCOL1, QCOL2,
     +   QGMAP, QRMAP, QROW1, QROW2, QRPIV, QRSCAL, N1, N2, ROWS1,
     +   ROWS2)
      IF (ERROR) GO TO 9201

      ERROR = .NOT. (N .EQ. N1)
      IF (ERROR) GO TO 9202

C///  RESERVE SOME WORK SPACE.

C     SUBROUTINE RESERV
C    +  (ERROR, TEXT,
C    +   FATAL, QLAST, QMAX, QSIZE,
C    +   NAME, NUMBER, Q)

C     LOGICAL EXPECT(ROWS1)

      CALL RESERV (ERROR, TEXT, .FALSE., LLAST, LMAX, LSIZE,
     +   'QEXPEC', ROWS1, QEXPEC)
      IF (ERROR) GO TO 9203

C     LOGICAL OCCUPY(ROWS1)

      CALL RESERV (ERROR, TEXT, .FALSE., LLAST, LMAX, LSIZE,
     +   'QOCCUP', ROWS1, QOCCUP)
      IF (ERROR) GO TO 9203

C     REAL Y0(N)

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QY0', N, QY0)
      IF (ERROR) GO TO 9203

      ERROR = .NOT. (ILAST .LE. ISIZE .AND. LLAST .LE. LSIZE .AND.
     +   RLAST .LE. RSIZE)
      IF (ERROR) GO TO 9204

C///  CALL.

C     SUBROUTINE BLCH2
C    +  (ERROR, TEXT, PRINT,
C    +   ACOL1, ACOLQ1, BLKS1, BMISS, BNULL, CMISS, CNULL, COL1, COLS1,
C    +   EXPECT, N, OCCUPY, RETURN, ROW1, ROWS1, X, Y, Y0)

2010  CONTINUE

      CALL BLCH2
     +  (ERROR, TEXT, PRINT - 1,
     +   IDATA(QAC1), IDATA(QACQ1), BLKS1, BMISS, BNULL, CMISS, CNULL,
     +   IDATA(QCOL1), COLS1, LWORK(QEXPEC), N, LWORK(QOCCUP),
     +   RETURN, IDATA(QROW1), ROWS1, X, Y, RWORK(QY0))
      IF (ERROR) GO TO 9205

      IF (RETURN) GO TO 99999

C///////////////////////////////////////////////////////////////////////
C
C     (3) EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) THEN
         IF (0 .LT. CMISS .AND. 0 .LT. CNULL) THEN
            WRITE (TEXT, 10002)
         ELSE IF (0 .LT. CMISS) THEN
            WRITE (TEXT, 10003)
         ELSE IF (0 .LT. CNULL) THEN
            WRITE (TEXT, 10004)
         ELSE
            WRITE (TEXT, 10005)
         END IF

         IF (0 .LT. CMISS .OR. 0 .LT. CNULL)
     +      WRITE (TEXT, 10006) CMISS, BMISS, CNULL, BNULL
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +   (/1X, A, 'CHECK THE MATRIX SPARSITY PATTERN.')

10002 FORMAT
     +  (/10X, 'THE SPARSITY PATTERN IS DEFICIENT.  IT ALSO CONTAINS'
     +   /10X, 'SOME UNNECESSARY BLOCKS, BUT THESE MAY BE NEEDED WHEN'
     +   /10X, 'THE JACOBIAN MATRIX IS EVALUATED ELSEWHERE.')

10003 FORMAT
     +  (/10X, 'THE SPARSITY PATTERN IS DEFICIENT.')

10004 FORMAT
     +  (/10X, 'THE SPARSITY PATTERN CONTAINS SOME UNNECESSARY BLOCKS,'
     +   /10X, 'BUT THESE MAY BE NEEDED WHEN THE JACOBIAN MATRIX IS'
     +   /10X, 'EVALUATED ELSEWHERE.')

10005 FORMAT
     +  (/10X, 'THE SPARSITY PATTERN IS CORRECT.')

10006 FORMAT
     +  (/10X, I10, '  DEFICIENT COLUMNS'
     +   /10X, I10, '  TOTAL MISSING BLOCKS'
     +  //10X, I10, '  COLUMNS WITH UNNECESSARY BLOCKS'
     +   /10X, I10, '  TOTAL UNNECESSARY BLOCKS')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99998

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   IDSIZ, ISIZE, LSIZE, N, RDSIZ, RSIZE
      GO TO 99998

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID
      GO TO 99998

9202  IF (0 .LT. TEXT) WRITE (TEXT, 99202) ID, N, N1
      GO TO 99998

9203  IF (0 .LT. TEXT) WRITE (TEXT, 99203) ID
      GO TO 99998

9204  IF (0 .LT. TEXT) WRITE (TEXT, 99204) ID,
     +   0, ISIZE, LSIZE, RSIZE, 0, ILAST, LLAST, RLAST
      GO TO 99998

9205  IF (0 .LT. TEXT) WRITE (TEXT, 99205) ID
      GO TO 99998

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  IDSIZ'
     +    /10X, I10, '  ISIZE'
     +    /10X, I10, '  LSIZE'
     +    /10X, I10, '  N',
     +    /10X, I10, '  RDSIZ'
     +    /10X, I10, '  RSIZE')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  BLDAT FAILS.')

99202 FORMAT
     +   (/1X, A9, 'ERROR.  THE MATRIX ORDER IS INCONSISTENT.'
     +   //10X, I10,  '  ORDER IN ARGUMENT LIST'
     +    /10X, I10,  '  IN DATA FROM BLSYM')

99203 FORMAT
     +   (/1X, A9, 'ERROR.  RESERV FAILS.')

99204 FORMAT
     +   (/1X, A9, 'ERROR.  SOME WORKSPACES ARE TOO SMALL.'
C               123456789_  123456789_  123456789_  123456789_
     +  //25X, ' CHARACTER     INTEGER     LOGICAL        REAL'
C               123456789_12345
     +  //10X, ' PRESENT SIZE', 4I12
     +   /10X, 'REQUIRED SIZE', 4I12
     +  //10X, 'THIS SUBROUTINE WILL NEED NO MORE SPACE.')

99205 FORMAT
     +   (/1X, A9, 'ERROR.  BLCH2 FAILS.')

C///  EXIT.

99998 CONTINUE
      ILAST = ISAVE
      LLAST = LSAVE
      RLAST = RSAVE
99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLDAT
     +  (ERROR, TEXT,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   ASIZE2, BLKS1, BLKS2, BORDER, COLS1, COLS2, GLUES, QA2, QAC1,
     +   QAC2, QACQ1, QACQ2, QAR2, QARQ2, QBMAP, QCMAP, QCOL1, QCOL2,
     +   QGMAP, QRMAP, QROW1, QROW2, QRPIV, QRSCAL, N1, N2, ROWS1,
     +   ROWS2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLDAT
C
C>    CALLED BY BLFAC, BLJAC, BLSOL
C
C>    UNPACK THE DATA SPACE.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   ASIZE2, BLKS1, BLKS2, BORDER, COLS1, COLS2, GLUES, IDATA,
     +   IDSIZ, ILAST, IMAX, ISIZE, N1, N2, RDSIZ, RLAST, RMAX, ROWS1,
     +   ROWS2, RSIZE, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   RDATA

      PARAMETER (ID = ' BLDAT:  ')

      DIMENSION IDATA(IDSIZ), RDATA(RDSIZ)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. IDSIZ .AND. 0 .LT. RDSIZ)
      IF (ERROR) GO TO 9101

C///  INITIALIZE THE STACK VARIABLES.

      ILAST = 0
      IMAX = 0
      ISIZE = IDSIZ

      RLAST = 0
      RMAX = 0
      RSIZE = RDSIZ

C///////////////////////////////////////////////////////////////////////
C
C     (2) UNPACK THE INTEGER SCALARS.
C
C///////////////////////////////////////////////////////////////////////

C     SUBROUTINE RESERV
C    +  (ERROR, TEXT,
C    +   FATAL, QLAST, QMAX, QSIZE,
C    +   NAME, NUMBER, Q)

C     INTEGER IDATA'S 11 SCALARS

      CALL RESERV (ERROR, TEXT, .TRUE., ILAST, IMAX, ISIZE,
     +   'QISCLR', 11, QISCLR)
      IF (ERROR) GO TO 9201

      ASIZE2 = IDATA(QISCLR + 00)
      BLKS1  = IDATA(QISCLR + 01)
      BLKS2  = IDATA(QISCLR + 02)
      BORDER = IDATA(QISCLR + 03)
      COLS1  = IDATA(QISCLR + 04)
      COLS2  = IDATA(QISCLR + 05)
      GLUES  = IDATA(QISCLR + 06)
      N1     = IDATA(QISCLR + 07)
      N2     = IDATA(QISCLR + 08)
      ROWS1  = IDATA(QISCLR + 09)
      ROWS2  = IDATA(QISCLR + 10)

      ERROR = .NOT. (0 .LT. ASIZE2 .AND. 0 .LT. BLKS1 .AND. 0 .LT. BLKS2
     +   .AND. 0 .LE. BORDER .AND. 0 .LT. COLS1 .AND. 0 .LT. COLS2 .AND.
     +   0 .LE. GLUES .AND. 0 .LT. N1 .AND. 0 .LT. N2 .AND. 0 .LT. ROWS1
     +   .AND. 0 .LT. ROWS2)
      IF (ERROR) GO TO 9202

      IF (0 .EQ. BORDER) THEN
         ERROR = .NOT. (BLKS1 .LE. BLKS2 .AND. COLS1 .EQ. COLS2 .AND.
     +      COLS1 .EQ. ROWS1 .AND. COLS1 .EQ. ROWS2 .AND. GLUES .EQ. 0
     +      .AND. N1 .EQ. N2)
         IF (ERROR) GO TO 9203
      ELSE
         ERROR = .NOT. (BLKS1 .LT. BLKS2 .AND. COLS1 .EQ. ROWS1 .AND.
     +      COLS1 .EQ. ROWS2 .AND. 2 * COLS1 .EQ. COLS2 .AND.
     +      COLS1 - 1 .EQ. GLUES .AND. N1 + BORDER * ROWS1 .EQ. N2)
         IF (ERROR) GO TO 9204
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (3) UNPACK THE ARRAYS.
C
C///////////////////////////////////////////////////////////////////////

C     SUBROUTINE RESERV
C    +  (ERROR, TEXT,
C    +   FATAL, QLAST, QMAX, QSIZE,
C    +   NAME, NUMBER, Q)

C     INTEGER ACOL1(2, BLKS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QAC1', 2 * BLKS1, QAC1)
      IF (ERROR) GO TO 9201

C     INTEGER ACOLQ1(2, COLS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QACQ1', 2 * COLS1, QACQ1)
      IF (ERROR) GO TO 9201

C     INTEGER ACOLQ2(4, COLS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QACQ2', 4 * COLS2, QACQ2)
      IF (ERROR) GO TO 9201

C     INTEGER BMAP(COLS1 + 1)

      IF (0 .LT. BORDER) THEN
         CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +      'QBMAP', COLS1 + 1, QBMAP)
         IF (ERROR) GO TO 9201
      END IF

C     INTEGER CMAP(COLS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QCMAP', COLS1, QCMAP)
      IF (ERROR) GO TO 9201

C     INTEGER COL1(2, COLS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QCOL1', 2 * COLS1, QCOL1)
      IF (ERROR) GO TO 9201

C     INTEGER COL2(2, COLS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QCOL2', 2 * COLS2, QCOL2)
      IF (ERROR) GO TO 9201

C     INTEGER GMAP(3, 0 : GLUES)

      IF (0 .LT. BORDER) THEN
         CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +      'QGMAP', 3 * (1 + GLUES), QGMAP)
         IF (ERROR) GO TO 9201
      END IF

C     INTEGER RMAP(ROWS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QRMAP', ROWS1, QRMAP)
      IF (ERROR) GO TO 9201

C     INTEGER ROW1(2, ROWS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QROW1', 2 * ROWS1, QROW1)
      IF (ERROR) GO TO 9201

C     INTEGER ROW2(2, ROWS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QROW2', 2 * ROWS2, QROW2)
      IF (ERROR) GO TO 9201

C     INTEGER ACOL2(2, BLKS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QAC2', 2 * BLKS2, QAC2)
      IF (ERROR) GO TO 9201

C     REAL A2(ASIZE2)

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QA2', ASIZE2, QA2)
      IF (ERROR) GO TO 9201

C     INTEGER AROW2(2, BLKS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QAR2', 2 * BLKS2, QAR2)
      IF (ERROR) GO TO 9201

C     INTEGER AROWQ2(4, ROWS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QARQ2', 4 * ROWS2, QARQ2)
      IF (ERROR) GO TO 9201

C     INTEGER RPIVOT(N2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QRPIV', N2, QRPIV)
      IF (ERROR) GO TO 9201

C     REAL RSCALE(N2)

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QRSCAL', N2, QRSCAL)
      IF (ERROR) GO TO 9201

      ERROR = .NOT. (ILAST .EQ. IDSIZ .AND. RLAST .EQ. RDSIZ)
      IF (ERROR) GO TO 9301

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, IDSIZ, RDSIZ
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID
      GO TO 99999

9202  IF (0 .LT. TEXT) WRITE (TEXT, 99202) ID,
     +   ASIZE2, BLKS1, BLKS2, BORDER, COLS1, COLS2, GLUES, N1, N2,
     +   ROWS1, ROWS2
      GO TO 99999

9203  IF (0 .LT. TEXT) WRITE (TEXT, 99203) ID,
     +   BLKS1, BLKS2, COLS1, COLS2, ROWS1, ROWS2, GLUES, N1, N2
      GO TO 99999

9204  IF (0 .LT. TEXT) WRITE (TEXT, 99204) ID,
     +   BLKS1, BLKS2, COLS1, ROWS1, ROWS2, 2 * COLS1, COLS2, COLS1 - 1,
     +   GLUES, N1, N1 + BORDER * ROWS1, N2
      GO TO 99999

9301  IF (0 .LT. TEXT) WRITE (TEXT, 99301) ID,
     +   ILAST, RLAST, IDSIZ, RDSIZ
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  IDSIZ'
     +    /10X, I10, '  RDSIZ')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  RESERV FAILS.')

99202 FORMAT
     +  (/1X, A9, 'ERROR.  THE INTEGER DATA SPACE IS CORRUPT.'
     +  //10X, I10, '  ASIZE2'
     +   /10X, I10, '  BLKS1'
     +   /10X, I10, '  BLKS2'
     +   /10X, I10, '  BORDER'
     +   /10X, I10, '  COLS1'
     +   /10X, I10, '  COLS2'
     +   /10X, I10, '  GLUES (MAY BE ZERO)'
     +   /10X, I10, '  N1'
     +   /10X, I10, '  N2'
     +   /10X, I10, '  ROWS1'
     +   /10X, I10, '  ROWS2')

99203 FORMAT
     +  (/1X, A9, 'ERROR.  THE INTEGER DATA SPACE IS CORRUPT.'
     +  //10X, I10, '  BLKS1'
     +   /10X, I10, '  BLKS2'
     +  //10X, I10, '  COLS1'
     +   /10X, I10, '  COLS2'
     +   /10X, I10, '  ROWS1'
     +   /10X, I10, '  ROWS2'
     +  //10X, I10, '  GLUES'
     +  //10X, I10, '  N1'
     +   /10X, I10, '  N2')

99204 FORMAT
     +  (/1X, A9, 'ERROR.  THE INTEGER DATA SPACE IS CORRUPT.'
     +  //10X, I10, '  BLKS1'
     +   /10X, I10, '  BLKS2'
     +  //10X, I10, '  COLS1'
     +   /10X, I10, '  ROWS1'
     +   /10X, I10, '  ROWS2'
     +  //10X, I10, '  2 * COLS1'
     +   /10X, I10, '  COLS2'
     +  //10X, I10, '  COLS1 - 1'
     +   /10X, I10, '  GLUES'
     +  //10X, I10, '  N1'
     +   /10X, I10, '  N1 + BORDER * ROWS1'
     +   /10X, I10, '  N2')

99301 FORMAT
     +   (/1X, A9, 'ERROR.  THE DATA SPACES ARE CORRUPT.'
C               123456789_  123456789_
     +  //27X, '   INTEGER        REAL'
C               123456789_1234567
     +   /10X, '    DIMENSIONED  ', 2I12
     +   /10X, '  EXPECTED SIZE  ', 2I12)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLEP2
     +   (SAME, VALUE1, VALUE2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLEP2
C
C>    CALLED BY BLEPS
C
C>    TEST EQUALITY OF FLOATING POINT NUMBERS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   VALUE1, VALUE2
      LOGICAL SAME

      SAME = VALUE1 .EQ. VALUE2

      RETURN
      END
      SUBROUTINE BLEPS
     +   (EPS)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLEPS
C
C>    CALLED BY BLCH2 AND BLJA2.
C
C>    FIND MACHINE EPSILON.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   EPS, VALUE
      INTEGER
     +   SCRTCH
      LOGICAL
     +   SAME

      PARAMETER (SCRTCH = 98)

C///  IEEE STANDARD

C*****PRECISION > DOUBLE
      VALUE = 1.1102230246251565D-16
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      VALUE = 5.9604645E-08
C*****END PRECISION > SINGLE

C*****MACHINE EPSILON > IEEE STANDARD
      EPS = VALUE
C*****END MACHINE EPSILON > IEEE STANDARD

C///  COMPUTED

C*****MACHINE EPSILON > COMPUTED
C      OPEN (ACCESS = 'SEQUENTIAL', FORM = 'UNFORMATTED',
C     +   STATUS = 'SCRATCH', UNIT = SCRTCH)
C
C      EPS = 1
C1010  CONTINUE
C      EPS = 0.5 * EPS
C
C      VALUE = 1 + EPS
C
C      REWIND (SCRTCH)
C      WRITE (SCRTCH) VALUE
C
C      REWIND (SCRTCH)
C      READ (SCRTCH) VALUE
C
C      SAME = 1 .EQ. VALUE
C
C      IF (.NOT. SAME) GO TO 1010
C
C      CLOSE (UNIT = SCRTCH)
C*****END MACHINE EPSILON > COMPUTED

C///  EXIT.

      RETURN
      END
      SUBROUTINE BLFA2
     +  (ERROR, TEXT, PRINT,
     +   A2, ACOL1, ACOL2, ACOLQ1, ACOLQ2, AROW2, AROWQ2, ASIZE2, BLKS1,
     +   BLKS2, BORDER, CMAP, COL2, COLS1, COLS2, CSUM, N2, RMAP, ROW2,
     +   ROWS1, ROWS2, RPIVOT, RSCALE)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLFA2
C
C>    CALLED BY BLFAC
C
C>    FACTOR THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   AC, AC1, AC2, ACOL1, ACOL2, ACOLQ1, ACOLQ2, AR, AR1, AR2,
     +   AROW2, AROWQ2, ASIZE2, BLKS1, BLKS2, BLOCK, BORDER, C1, C2,
     +   CMAP, COL2, COLS1, COLS2, COUNT, COUNT1, COUNT2, FIRST, INDEX,
     +   J, LAST, LC, LC1, LC2, LR1, LR2, N2, PRINT, R1, R2, RMAP, ROW2,
     +   ROWS1, ROWS2, RPIVOT, TEXT, UC1, UC2, UR, UR1, UR2
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A2, CSUM, RSCALE

      PARAMETER (ID = ' BLFA2:  ')

      PARAMETER (QFRST = 1, QLAST = 2, QFRSTD = 3, QLASTD = 4)
      PARAMETER (QIND = 1, QLOC = 2)

      DIMENSION
     +   A2(ASIZE2), ACOL1(2, BLKS1), ACOL2(2, BLKS2), ACOLQ1(2, COLS1),
     +   ACOLQ2(4, COLS2), AROW2(2, BLKS2), AROWQ2(4, ROWS2),
     +   CMAP(COLS1), COL2(2, COLS2), COUNT(3), CSUM(N2), INDEX(5),
     +   RMAP(ROWS1), ROW2(2, ROWS2), RPIVOT(N2), RSCALE(N2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. ASIZE2 .AND. 0 .LT. BLKS1 .AND. 0 .LT. BLKS2
     +   .AND. 0 .LT. COLS1 .AND. 0 .LT. COLS2 .AND. 0 .LT. N2 .AND.
     +   0 .LT. ROWS1 .AND. 0 .LT. ROWS2)
      IF (ERROR) GO TO 9101

C///  GATHER SOME DATA.

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) THEN

      COUNT(1) = 0
      DO 1010 J = 1, ASIZE2
         IF (A2(J) .NE. 0.0) COUNT(1) = COUNT(1) + 1
1010  CONTINUE

      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (2) CHECK THE COLUMN SUMS.
C
C///////////////////////////////////////////////////////////////////////

C///  PREPARE THE SUMS.

      DO 2010 J = 1, N2
         CSUM(J) = 0.0
2010  CONTINUE

C///  TOP OF THE LOOP OVER THE A2 BLOCKS WITH A1 ENTRIES.

      DO 2030 C1 = 1, COLS1
         C2 = CMAP(C1)
         AC1 = COL2(QFRST, C2)
         AC2 = COL2(QLAST, C2)
         DO 2020 BLOCK = ACOLQ1(QFRST, C1), ACOLQ1(QLAST, C1)
            Q = ACOL1(QLOC, BLOCK)
            R1 = ACOL1(QIND, BLOCK)
            R2 = RMAP(R1)
            AR1 = ROW2(QFRST, R2)
            AR2 = ROW2(QLAST, R2)

C///  SUM THE COLUMNS OF THE A2 BLOCK.

C     SUBROUTINE BLFA3
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, AR1, AR2, CSUM)

      CALL BLFA3
     +  (ERROR, TEXT,
     +   A2(Q), AC1, AC2, AR1, AR2, CSUM(AC1))
      IF (ERROR) GO TO 9201

C///  BOTTOM OF THE LOOP OVER THE A2 BLOCKS WITH A1 ENTRIES.

2020     CONTINUE
2030  CONTINUE

C///  TOP OF THE LOOP OVER THE BORDER COLUMNS, IF ANY.

      IF (1 .LE. BORDER) THEN
         AC1 = COL2(QFRST, COLS2)
         AC2 = COL2(QLAST, COLS2)
         DO 2040 BLOCK = ACOLQ2(QFRST, COLS2), ACOLQ2(QLAST, COLS2)
            Q = ACOL2(QLOC, BLOCK)
            R2 = ACOL2(QIND, BLOCK)
            AR1 = ROW2(QFRST, R2)
            AR2 = ROW2(QLAST, R2)

C///  SUM THE BORDER COLUMNS.

C     SUBROUTINE BLFA3
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, AR1, AR2, CSUM)

      CALL BLFA3
     +  (ERROR, TEXT,
     +   A2(Q), AC1, AC2, AR1, AR2, CSUM(AC1))
      IF (ERROR) GO TO 9201

C///  BOTTOM OF THE LOOP OVER THE BORDER COLUMNS.

2040     CONTINUE
      END IF

C///  CHECK FOR NULL COLUMNS.

      COUNT1 = 0
      DO 2050 J = 1, N2
         IF (CSUM(J) .EQ. 0.0) COUNT1 = COUNT1 + 1
2050  CONTINUE

      ERROR = .NOT. (COUNT1 .EQ. 0)
      IF (ERROR) GO TO 9202

C///////////////////////////////////////////////////////////////////////
C
C     (3) SCALE THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

C///  PREPARE THE SCALE FACTORS.

      COUNT1 = 0
      COUNT2 = 0
      DO 3010 J = 1, N2
         IF (RSCALE(J) .LT. 0.0) THEN
            COUNT1 = COUNT1 + 1
         ELSE IF (RSCALE(J) .EQ. 0.0) THEN
            COUNT2 = COUNT2 + 1
         END IF
3010  CONTINUE
      ERROR = .NOT. (COUNT1 .EQ. 0 .AND. COUNT2 .EQ. 0)
      IF (ERROR) GO TO 9301

      DO 3020 J = 1, N2
         RSCALE(J) = 1.0 / RSCALE(J)
3020  CONTINUE

C///  TOP OF THE LOOP OVER THE A2 BLOCKS WITH A1 ENTRIES.

      DO 3040 C1 = 1, COLS1
         C2 = CMAP(C1)
         AC1 = COL2(QFRST, C2)
         AC2 = COL2(QLAST, C2)
         DO 3030 BLOCK = ACOLQ1(QFRST, C1), ACOLQ1(QLAST, C1)
            Q = ACOL1(QLOC, BLOCK)
            R1 = ACOL1(QIND, BLOCK)
            R2 = RMAP(R1)
            AR1 = ROW2(QFRST, R2)
            AR2 = ROW2(QLAST, R2)

C///  SCALE THE A2 BLOCK INCLUDING STRETCHED ROWS, IF ANY.

C     SUBROUTINE BLFA4
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, AR1, AR2, RSCALE)

      CALL BLFA4
     +  (ERROR, TEXT,
     +   A2(Q), AC1, AC2, AR1, AR2, RSCALE(AR1))
      IF (ERROR) GO TO 9302

C///  BOTTOM OF THE LOOP OVER THE A2 BLOCKS WITH A1 ENTRIES.

3030     CONTINUE
3040  CONTINUE

C///  TOP OF THE LOOP OVER THE BORDER COLUMNS, IF ANY.

      IF (1 .LE. BORDER) THEN
         AC1 = COL2(QFRST, COLS2)
         AC2 = COL2(QLAST, COLS2)
         DO 3050 BLOCK = ACOLQ2(QFRST, COLS2), ACOLQ2(QLAST, COLS2)
            Q = ACOL2(QLOC, BLOCK)
            R2 = ACOL2(QIND, BLOCK)
            AR1 = ROW2(QFRST, R2)
            AR2 = ROW2(QLAST, R2)

C///  SCALE THE BORDER COLUMNS.

C     SUBROUTINE BLFA4
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, AR1, AR2, RSCALE)

      CALL BLFA4
     +  (ERROR, TEXT,
     +   A2(Q), AC1, AC2, AR1, AR2, RSCALE(AR1))
      IF (ERROR) GO TO 9302

C///  BOTTOM OF THE LOOP OVER THE BORDER COLUMNS.

3050     CONTINUE
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (4) FACTOR THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP OVER THE BLOCKS TO BE CHANGED.

      DO 4030 AC = 1, COLS2
         UC1 = COL2(QFRST, AC)
         UC2 = COL2(QLAST, AC)
         DO 4020 QAR = ACOLQ2(QFRST, AC), ACOLQ2(QLAST, AC)
            AR = ACOL2(QIND, QAR)
            QA = ACOL2(QLOC, QAR)
            LR1 = ROW2(QFRST, AR)
            LR2 = ROW2(QLAST, AR)

C///  TOP OF THE LOOP OVER THE ELIMINATION BLOCKS.

      QC = AROWQ2(QFRST, AR)
      LC = AROW2(QIND, QC)

      QR = ACOLQ2(QFRST, AC)
      UR = ACOL2(QIND, QR)

4010  CONTINUE

      QL = AROW2(QLOC, QC)
      LC1 = COL2(QFRST, LC)
      LC2 = COL2(QLAST, LC)

      QU = ACOL2(QLOC, QR)
      UR1 = ROW2(QFRST, UR)
      UR2 = ROW2(QLAST, UR)

      FIRST = MAX (LC1, UR1)
      LAST = MIN (LC2, UR2)

      IF (FIRST .LE. LAST) THEN

C///  A, L, AND U BLOCKS.

      IF (LC .LT. AC .AND. UR .LT. AR) THEN

C     CALL BLFA5
C    +  (ERROR, TEXT,
C    +   A, L, LC1, LC2, LR1, LR2, U, UC1, UC2, UR1, UR2)

      CALL BLFA5
     +  (ERROR, TEXT,
     +   A2(QA), A2(QL), LC1, LC2, LR1, LR2, A2(QU), UC1, UC2, UR1, UR2)
      IF (ERROR) GO TO 9401

C///  A AND L BLOCKS.

      ELSE IF (LC .LT. AC .AND. UR .EQ. AR) THEN

C     SUBROUTINE BLFA6
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, L, LC1, LC2, LR1, LR2, RPIVOT)

      CALL BLFA6
     +  (ERROR, TEXT,
     +   A2(QA), UC1, UC2, A2(QL), LC1, LC2, LR1, LR2, RPIVOT(LR1))
      IF (ERROR) GO TO 9402

C///  A AND U BLOCKS.

      ELSE IF (LC .EQ. AC .AND. UR .LT. AR) THEN

C     SUBROUTINE BLFA7
C    +  (ERROR, TEXT,
C    +   A, AR1, AR2, U, UC1, UC2, UR1, UR2)

      CALL BLFA7
     +  (ERROR, TEXT,
     +   A2(QA), LR1, LR2, A2(QU), UC1, UC2, UR1, UR2)
      IF (ERROR) GO TO 9403

C///  A BLOCK ONLY.

      ELSE

C     SUBROUTINE BLFA8
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, AR1, AR2, RPIVOT)

      CALL BLFA8
     +  (ERROR, TEXT,
     +   A2(QA), UC1, UC2, LR1, LR2, RPIVOT(LR1))
      IF (ERROR) GO TO 9404

      END IF

C///  BOTTOM OF THE LOOP OVER THE ELIMINATION BLOCKS.

      END IF

      IF (LAST .EQ. LC2) QC = QC + 1
      IF (LAST .EQ. UR2) QR = QR + 1

      IF (QC .LE. AROWQ2(QLAST, AR) .AND.
     +    QR .LE. ACOLQ2(QLAST, AC)) THEN
         LC = AROW2(QIND, QC)
         UR = ACOL2(QIND, QR)
         IF (LC .LE. AC .AND. UR .LE. AR) GO TO 4010
      END IF

C///  BOTTOM OF THE LOOP OVER THE BLOCKS TO BE CHANGED.

4020     CONTINUE
4030  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (5) EPILOGUE
C
C///////////////////////////////////////////////////////////////////////

C///  GATHER SOME DATA.

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) THEN

      COUNT(2) = 0
      DO 5010 J = 1, ASIZE2
         IF (A2(J) .NE. 0.0) COUNT(2) = COUNT(2) + 1
5010  CONTINUE

      COUNT(3) = 0
      DO 5020 J = 1, N2
         IF (RPIVOT(J) .NE. J) COUNT(3) = COUNT(3) + 1
5020  CONTINUE

C///  PRINT.

      WRITE(TEXT, 10001)
     +   COUNT(1), COUNT(2), ASIZE2,
     +   NINT (100.0 * REAL (COUNT(2)) / REAL (ASIZE2)),
     +   COUNT(3), N2, NINT (100.0 * REAL (COUNT(3)) / REAL (N2))

      END IF

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +   (/10X, I10, '  NONZEROS IN MATRIX'
     +    /10X, I10, '  NONZEROS IN FACTORS'
     +    /10X, I10, '  MATRIX STORAGE SPACE'
     +    /10X, I10, '  PERCENT UTILIZATION'
     +   //10X, I10, '  ROW EXCHANGES'
     +    /10X, I10, '  MATRIX ORDER'
     +    /10X, I10, '  PERCENT OF POSSIBLE EXCHANGES')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE(TEXT, 99101) ID,
     +   ASIZE2, BLKS1, BLKS2, COLS1, COLS2, N2, ROWS1, ROWS2
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE(TEXT, 99201) ID
      GO TO 99999

9202  IF (0 .LT. TEXT) THEN
         IF (50 .LT. COUNT1) THEN
            WRITE(TEXT, 99202) ID, COUNT1, N2, 'A PARTIAL'
         ELSE
            WRITE(TEXT, 99202) ID, COUNT1, N2, 'THE'
         END IF
         COUNT1 = 0
         DO 9292 C2 = 1, COLS2
            AC1 = COL2(QFRST, C2)
            AC2 = COL2(QLAST, C2)
            DO 9291 J = AC1, AC2
               IF (CSUM(J) .EQ. 0.0) THEN
                  COUNT1 = COUNT1 + 1
                  IF (50 .LT. COUNT1) GO TO 9293
                  WRITE (TEXT, 99291) C2, AC2 - AC1 + 1, J - AC1 + 1, J
               END IF
9291        CONTINUE
9292     CONTINUE
9293     CONTINUE
      END IF
      GO TO 99999

9301  IF (0 .LT. TEXT) THEN
         WRITE(TEXT, 99301) ID, COUNT1, COUNT2, N2
         IF (0 .LT. COUNT2) THEN
            WRITE(TEXT, 99391)
            COUNT2 = 0
            DO 9391 J = 1, N2
               IF (RSCALE(J) .EQ. 0.0) THEN
                  IF (COUNT2 .EQ. 5) THEN
                     WRITE(TEXT, 99392) INDEX
                     COUNT2 = 0
                  END IF
                  COUNT2 = COUNT2 + 1
                  INDEX(COUNT2) = J
               END IF
9391        CONTINUE
            IF (0 .LT. COUNT2)
     +         WRITE(TEXT, 99392) (INDEX(J), J = 1, COUNT2)
         END IF
      END IF
      GO TO 99999

9302  IF (0 .LT. TEXT) WRITE(TEXT, 99302) ID
      GO TO 99999

9401  IF (0 .LT. TEXT) WRITE(TEXT, 99401) ID
      GO TO 99999

9402  IF (0 .LT. TEXT) WRITE(TEXT, 99402) ID
      GO TO 99999

9403  IF (0 .LT. TEXT) WRITE(TEXT, 99403) ID
      GO TO 99999

9404  IF (0 .LT. TEXT) WRITE(TEXT, 99404) ID
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  ASIZE2'
     +    /10X, I10, '  BLKS1'
     +    /10X, I10, '  BLKS2'
     +    /10X, I10, '  COLS1'
     +    /10X, I10, '  COLS2'
     +    /10X, I10, '  N2'
     +    /10X, I10, '  ROWS1'
     +    /10X, I10, '  ROWS2')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  BLFA3 FAILS.')

99202 FORMAT
     +   (/1X, A9, 'ERROR.  SOME COLUMNS ARE ZERO.'
     +   //10X, I10, '  ZERO COLUMNS'
     +    /10X, I10, '  N2'
     +   //10X, A, ' LIST OF ZERO COLUMNS FOLLOWS.'
C                1234567890  1234567890  1234567890  1234567890
     +   //10X, '                        COLUMN INDEX'
     +    /10X, '                        ----------------------'
     +    /10X, '     BLOCK  BLOCK SIZE    IN BLOCK   IN MATRIX'
     +    /)

99291 FORMAT
     +   (10X, I10, 4(2X, I10))

99301 FORMAT
     +   (/1X, A9, 'ERROR.  SOME ABSOLUTE ROW SUMS ARE NOT POSITIVE.'
     +   //10X, I10, '  NEGATIVE SUMS'
     +    /10X, I10, '  ZERO SUMS'
     +    /10X, I10, '  N2')

99391 FORMAT
     +   (/10X, 'ROWS WITH ZERO SUMS:'
     +    /)

99392 FORMAT
     +   (10X, I10, 4(2X, I10))

99302 FORMAT
     +   (/1X, A9, 'ERROR.  BLFA4 FAILS.')

99401 FORMAT
     +   (/1X, A9, 'ERROR.  BLFA5 FAILS.')

99402 FORMAT
     +   (/1X, A9, 'ERROR.  BLFA6 FAILS.')

99403 FORMAT
     +   (/1X, A9, 'ERROR.  BLFA7 FAILS.')

99404 FORMAT
     +   (/1X, A9, 'ERROR.  BLFA8 FAILS.')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLFA3
     +  (ERROR, TEXT,
     +   A, AC1, AC2, AR1, AR2, CSUM)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLFA3
C
C>    CALLED BY BLFA2
C
C>    SUM A BLOCK'S COLUMNS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER AC, AC1, AC2, AR, AR1, AR2, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, CSUM

      PARAMETER (ID = ' BLFA3:  ')

      DIMENSION A(AR1 : AR2, AC1 : AC2), CSUM(AC1 : AC2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (AC1 .LE. AC2 .AND. AR1 .LE. AR2)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) SUM THE BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      DO 2020 AC = AC1, AC2
         DO 2010 AR = AR1, AR2
            CSUM(AC) = CSUM(AC) + ABS (A(AR, AC))
2010     CONTINUE
2020  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, AC1, AC2, AR1, AR2
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLFA4
     +  (ERROR, TEXT,
     +   A, AC1, AC2, AR1, AR2, RSCALE)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLFA4
C
C>    CALLED BY BLFA2
C
C>    SCALE A BLOCK'S ROWS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER AC, AC1, AC2, AR, AR1, AR2, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, RSCALE

      PARAMETER (ID = ' BLFA4:  ')

      DIMENSION A(AR1 : AR2, AC1 : AC2), RSCALE(AR1 : AR2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (AC1 .LE. AC2 .AND. AR1 .LE. AR2)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) SCALE THE BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      DO 2020 AC = AC1, AC2
         DO 2010 AR = AR1, AR2
            A(AR, AC) = RSCALE(AR) * A(AR, AC)
2010     CONTINUE
2020  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, AC1, AC2, AR1, AR2
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLFA5
     +  (ERROR, TEXT,
     +   A, L, LC1, LC2, LR1, LR2, U, UC1, UC2, UR1, UR2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLFA5
C
C>    CALLED BY BLFA2
C
C>    PERFORM ELIMINATION WITH LOWER AND UPPER BLOCKS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   FIRST, J, LAST, LC1, LC2, LR1, LR2, TEXT, UC, UC1, UC2, UR,
     +   UR1, UR2
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, L, TEMP, U

      PARAMETER (ID = ' BLFA5:  ')

      DIMENSION
     +   A(LR1 : LR2, UC1 : UC2), L(LR1 : LR2, LC1 : LC2),
     +   U(UR1 : UR2, UC1 : UC2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (LC1 .LE. LC2 .AND. LR1 .LE. LR2 .AND. UC1 .LE. UC2
     +   .AND. UR1 .LE. UR2)
      IF (ERROR) GO TO 9101

C///  CHECK THE RANGE.

      FIRST = MAX (LC1, UR1)
      LAST = MIN (LC2, UR2)
      ERROR = .NOT. (FIRST .LE. LAST)
      IF (ERROR) GO TO 9102

C///////////////////////////////////////////////////////////////////////
C
C     (2) PERFORM THE ELIMINATION.
C
C///////////////////////////////////////////////////////////////////////

      DO 2030 UC = UC1, UC2
         DO 2020 UR = FIRST, LAST
            TEMP = U(UR, UC)
            DO 2010 J = LR1, LR2
               A(J, UC) = A(J, UC) - TEMP * L(J, UR)
2010        CONTINUE
2020     CONTINUE
2030  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   LC1, LC2, LR1, LR2, UC1, UC2, UR1, UR2
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, LC1, LC2, UR1, UR2
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  LC1'
     +    /10X, I10, '  LC2'
     +   //10X, I10, '  LR1'
     +    /10X, I10, '  LR2'
     +   //10X, I10, '  UC1'
     +    /10X, I10, '  UC2'
     +   //10X, I10, '  UR1'
     +    /10X, I10, '  UR2')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  THE L COLUMNS AND U ROWS ARE DISJOINT.'
     +   //10X, I10, '  LC1'
     +    /10X, I10, '  LC2'
     +   //10X, I10, '  UR1'
     +    /10X, I10, '  UR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLFA6
     +  (ERROR, TEXT,
     +   A, AC1, AC2, L, LC1, LC2, LR1, LR2, RPIVOT)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLFA6
C
C>    CALLED BY BLFA2
C
C>    PERFORM ELIMINATION WITH A LOWER BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   AC, AC1, AC2, AR, FIRST, J, LAST, LC1, LC2, LR1, LR2, RINDEX,
     +   RPIVOT, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, L, TEMP

      PARAMETER (ID = ' BLFA6:  ')

      DIMENSION
     +   A(LR1 : LR2, AC1 : AC2), L(LR1 : LR2, LC1 : LC2),
     +   RPIVOT(LR1 : LR2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (AC1 .LE. AC2 .AND. LC1 .LE. LC2 .AND. LR1 .LE. LR2)
      IF (ERROR) GO TO 9101

C///  CHECK THE RANGE.

      FIRST = MAX (LC1, LR1)
      LAST = MIN (LC2, LR2)
      ERROR = .NOT. (FIRST .LE. LAST)
      IF (ERROR) GO TO 9102

C///////////////////////////////////////////////////////////////////////
C
C     (2) PERFORM THE ELIMINATION.
C
C///////////////////////////////////////////////////////////////////////

      DO 2030 AC = AC1, AC2
         DO 2020 AR = FIRST, MIN (AC, LAST)
            RINDEX = RPIVOT(AR)
            TEMP = A(RINDEX, AC)
            A(RINDEX, AC) = A(AR, AC)
            A(AR, AC) = TEMP

            DO 2010 J = AR + 1, LR2
               A(J, AC) = A(J, AC) - TEMP * L(J, AR)
2010        CONTINUE
2020     CONTINUE
2030  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   AC1, AC2, LC1, LC2, LR1, LR2
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, LC1, LC2, LR1, LR2
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  LC1'
     +    /10X, I10, '  LC2'
     +   //10X, I10, '  LR1'
     +    /10X, I10, '  LR2')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  THE L COLUMNS AND U ROWS ARE DISJOINT.'
     +   //10X, I10, '  LC1'
     +    /10X, I10, '  LC2'
     +   //10X, I10, '  LR1'
     +    /10X, I10, '  LR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLFA7
     +  (ERROR, TEXT,
     +   A, AR1, AR2, U, UC1, UC2, UR1, UR2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLFA7
C
C>    CALLED BY BLFA2
C
C>    PERFORM ELIMINATION WITH AN UPPER BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   AR, AR1, AR2, FIRST, LAST, TEXT, UC, UC1, UC2, UR, UR1, UR2
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, TEMP, U

      PARAMETER (ID = ' BLFA7:  ')

      DIMENSION A(AR1 : AR2, UC1 : UC2), U(UR1 : UR2, UC1 : UC2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (AR1 .LE. AR2 .AND. UC1 .LE. UC2 .AND. UR1 .LE. UR2)
      IF (ERROR) GO TO 9101

C///  CHECK THE RANGE.

      FIRST = MAX (UC1, UR1)
      LAST = MIN (UC2, UR2)
      ERROR = .NOT. (FIRST .LE. LAST)
      IF (ERROR) GO TO 9102

C///////////////////////////////////////////////////////////////////////
C
C     (2) PERFORM THE ELIMINATION.
C
C///////////////////////////////////////////////////////////////////////

      DO 2040 UC = FIRST, LAST
         DO 2020 UR = FIRST, UC - 1
            TEMP = U(UR, UC)
            DO 2010 AR = AR1, AR2
               A(AR, UC) = A(AR, UC) - A(AR, UR) * TEMP
2010        CONTINUE
2020     CONTINUE

         TEMP = U(UC, UC)
         DO 2030 AR = AR1, AR2
            A(AR, UC) = A(AR, UC) * TEMP
2030     CONTINUE
2040  CONTINUE

      DO 2070 UC = LAST + 1, UC2
         DO 2060 UR = FIRST, LAST
            TEMP = U(UR, UC)
            DO 2050 AR = AR1, AR2
               A(AR, UC) = A(AR, UC) - A(AR, UR) * TEMP
2050        CONTINUE
2060     CONTINUE
2070  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   AR1, AR2, UC1, UC2, UR1, UR2
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, UC1, UC2, UR1, UR2
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2'
     +   //10X, I10, '  UC1'
     +    /10X, I10, '  UC2'
     +   //10X, I10, '  UR1'
     +    /10X, I10, '  UR2')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  THE U COLUMNS AND U ROWS ARE DISJOINT.'
     +   //10X, I10, '  UC1'
     +    /10X, I10, '  UC2'
     +   //10X, I10, '  UR1'
     +    /10X, I10, '  UR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLFA8
     +  (ERROR, TEXT,
     +   A, AC1, AC2, AR1, AR2, RPIVOT)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLFA8
C
C>    CALLED BY BLFA2
C
C>    FACTOR A BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   AC, AC1, AC2, AR, AR1, AR2, FIRST, J, LAST, RINDEX, RPIVOT,
     +   TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, SWAP, TEMP

      PARAMETER (ID = ' BLFA8:  ')

      DIMENSION
     +   A(AR1 : AR2, AC1 : AC2), RPIVOT(AR1 : AR2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (AC1 .LE. AC2 .AND. AR1 .LE. AR2)
      IF (ERROR) GO TO 9101

C///  CHECK THE RANGE.

      FIRST = MAX (AC1, AR1)
      LAST = MIN (AC2, AR2)
      ERROR = .NOT. (FIRST .LE. LAST)
      IF (ERROR) GO TO 9102

C///////////////////////////////////////////////////////////////////////
C
C     (2) PERFORM THE FACTORIZATION.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP OVER THE COLUMNS.

      DO 2050 AC = FIRST, AC2

C///  ELIMINATION.

      DO 2020 AR = FIRST, MIN (AC - 1, LAST)
         RINDEX = RPIVOT(AR)

         TEMP = A(RINDEX, AC)
         A(RINDEX, AC) = A(AR, AC)
         A(AR, AC) = TEMP

         DO 2010 J = AR + 1, AR2
            A(J, AC) = A(J, AC) - A(J, AR) * TEMP
2010     CONTINUE
2020  CONTINUE

C///  ROW AND COLUMN PIVOTING.

      IF (AC .LE. LAST) THEN
         RINDEX = AC
         TEMP = ABS (A(AC, AC))
         DO 2030 AR = AC, AR2
            IF (ABS (A(AR, AC)) .GT. TEMP) THEN
               RINDEX = AR
               TEMP = ABS (A(AR, AC))
            END IF
2030     CONTINUE

         RPIVOT(AC) = RINDEX
         IF (AC .NE. RINDEX) THEN
            SWAP = A(AC, AC)
            A(AC, AC) = A(RINDEX, AC)
            A(RINDEX, AC) = SWAP
         END IF

         TEMP = A(AC, AC)
         ERROR = TEMP .EQ. 0.0
         IF (ERROR) GO TO 9201
         TEMP = 1.0 / TEMP
         A(AC, AC) = TEMP

         DO 2040 J = AC + 1, AR2
            A(J, AC) = A(J, AC) * TEMP
2040     CONTINUE
      END IF

C///  BOTTOM OF THE LOOP OVER THE COLUMNS.

2050  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, AC1, AC2, AR1, AR2
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, AC1, AC2, AR1, AR2
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID, AC, AC1, AC2, AR1, AR2
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  THE COLUMNS AND ROWS ARE DISJOINT.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  A PIVOT VANISHES.'
     +   //10X, I10, '  PIVOT'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLFAC
     +  (ERROR, TEXT, PRINT,
     +   RLAST, RMAX, RSIZE, RWORK,
     +   IDATA, IDSIZ, RDATA, RDSIZ)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLFAC
C
C>    CALLED BY THE  U S E R
C
C>    FACTOR THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   ASIZE2, BLKS1, BLKS2, BORDER, COLS1, COLS2, GLUES, IDATA,
     +   IDSIZ, N1, N2, PRINT, RDSIZ, RLAST, RMAX, ROWS1, ROWS2, RSAVE,
     +   RSIZE, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   RDATA, RWORK

      PARAMETER (ID = ' BLFAC:  ')

      DIMENSION IDATA(IDSIZ), RDATA(RDSIZ), RWORK(RSIZE)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  SAVE THE LEVELS OF THE WORK SPACES.

      RSAVE = RLAST

C///  PRINT.

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) WRITE (TEXT, 10001) ID

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. IDSIZ .AND. 0 .LT. RDSIZ .AND. 0 .LT. RSIZE)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) CALL BLFA2.
C
C///////////////////////////////////////////////////////////////////////

C///  UNPACK THE DATA SPACES.

      CALL BLDAT
     +  (ERROR, TEXT,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   ASIZE2, BLKS1, BLKS2, BORDER, COLS1, COLS2, GLUES, QA2, QAC1,
     +   QAC2, QACQ1, QACQ2, QAR2, QARQ2, QBMAP, QCMAP, QCOL1, QCOL2,
     +   QGMAP, QRMAP, QROW1, QROW2, QRPIV, QRSCAL, N1, N2, ROWS1,
     +   ROWS2)
      IF (ERROR) GO TO 9201

C///  RESERVE SOME WORK SPACE.

C     SUBROUTINE RESERV
C    +  (ERROR, TEXT,
C    +   FATAL, QLAST, QMAX, QSIZE,
C    +   NAME, NUMBER, Q)

C     REAL CSUM(N2)

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QCSUM', N2, QCSUM)
      IF (ERROR) GO TO 9202

      ERROR = .NOT. (RLAST .LE. RSIZE)
      IF (ERROR) GO TO 9203

C///  CALL.

C     SUBROUTINE BLFA2
C    +  (ERROR, TEXT, PRINT,
C    +   A2, ACOL1, ACOL2, ACOLQ1, ACOLQ2, AROW2, AROWQ2, ASIZE2, BLKS1,
C    +   BLKS2, BORDER, CMAP, COL2, COLS1, COLS2, CSUM, N2, RMAP, ROW2,
C    +   ROWS1, ROWS2, RPIVOT, RSCALE)

      CALL BLFA2
     +  (ERROR, TEXT, PRINT - 1,
     +   RDATA(QA2), IDATA(QAC1), IDATA(QAC2), IDATA(QACQ1),
     +   IDATA(QACQ2), IDATA(QAR2), IDATA(QARQ2), ASIZE2, BLKS1, BLKS2,
     +   BORDER, IDATA(QCMAP), IDATA(QCOL2), COLS1, COLS2, RWORK(QCSUM),
     +   N2, IDATA(QRMAP), IDATA(QROW2), ROWS1, ROWS2, IDATA(QRPIV),
     +   RDATA(QRSCAL))
      IF (ERROR) GO TO 9204

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +   (/1X, A, 'FACTOR A SPARSE BLOCKED MATRIX.')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, IDSIZ, RDSIZ, RSIZE
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID
      GO TO 99999

9202  IF (0 .LT. TEXT) WRITE (TEXT, 99202) ID
      GO TO 99999

9203  IF (0 .LT. TEXT) WRITE (TEXT, 99203) ID, RSIZE, RLAST
      GO TO 99999

9204  IF (0 .LT. TEXT) WRITE (TEXT, 99204) ID
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  IDSIZ'
     +    /10X, I10, '  RDSIZ'
     +    /10X, I10, '  RSIZE')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  BLDAT FAILS.')

99202 FORMAT
     +   (/1X, A9, 'ERROR.  RESERV FAILS.')

99203 FORMAT
     +   (/1X, A9, 'ERROR.  SOME WORKSPACES ARE TOO SMALL.'
C               123456789_  123456789_  123456789_  123456789_
     +  //25X, ' CHARACTER     INTEGER     LOGICAL        REAL'
C               123456789_12345
     +  //10X, ' PRESENT SIZE', 36X, I12
     +   /10X, 'REQUIRED SIZE', 36X, I12
     +  //10X, 'THIS SUBROUTINE WILL NEED NO MORE SPACE.')

99204 FORMAT
     +   (/1X, A9, 'ERROR.  BLFA2 FAILS.')

C///  EXIT.

99999 CONTINUE
      RLAST = RSAVE
      RETURN
      END
      SUBROUTINE BLJA2
     +  (ERROR, TEXT, PRINT,
     +   A2, ACOL1, ACOLQ1, ASIZE2, BLKS1, CMAP, COL1, COLS1, N, N2,
     +   NEXT, OCCUPY, POINT, RETURN, RMAP, ROW1, ROW2, ROWS1, ROWS2,
     +   RSCALE, SAVE, X, Y)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLJA2
C
C>    CALLED BY BLJAC
C
C>    EVALUATE A JACOBIAN MATRIX.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   AC, AC1, AC2, ACOL1, ACOLQ1, AR1, AR2, ASIZE2, BLKS1, C1, C2,
     +   CMAP, COL1, COLS1, COUNT, EVALS, INDEX, J, N, N2, NEXT, POINT,
     +   PRINT, R1, R2, RMAP, ROUTE, ROW1, ROW2, ROWS1, ROWS2, TEXT,
     +   YR1, YR2
      LOGICAL ERROR, OCCUPY, RETURN, TOUCH
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A2, ABSOL, EPS, FACTOR, RELAT, RSCALE, SAVE, X, Y

      PARAMETER (ID = ' BLJA2:  ')

      PARAMETER (QFRST = 1, QLAST = 2, QFRSTD = 3, QLASTD = 4)
      PARAMETER (QIND = 1, QLOC = 2)

      DIMENSION
     +   A2(ASIZE2), ACOL1(2, BLKS1), ACOLQ1(2, COLS1), CMAP(COLS1),
     +   COL1(2, COLS1), NEXT(COLS1), OCCUPY(ROWS1), POINT(COLS1),
     +   RMAP(ROWS1), ROW1(2, ROWS1), ROW2(2, ROWS2), RSCALE(N2),
     +   SAVE(COLS1), X(N), Y(N)

      SAVE

C///////////////////////////////////////////////////////////////////////
C
C     RETURN FROM REVERSE COMMUNICATION.
C
C///////////////////////////////////////////////////////////////////////

      IF (RETURN) THEN
         RETURN = .FALSE.
         GO TO (2030, 3070) ROUTE
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. ASIZE2 .AND. 0 .LT. BLKS1 .AND. 0 .LT. COLS1
     +   .AND. 0 .LT. N .AND. 0 .LT. N2 .AND. 0 .LT. ROWS1 .AND.
     +   0 .LT. ROWS2)
      IF (ERROR) GO TO 9101

C///  INITIALIZE THE COUNT OF FUNCTION EVALUATIONS.

      EVALS = 0

C///////////////////////////////////////////////////////////////////////
C
C     (2) INITIALIZE.
C
C///////////////////////////////////////////////////////////////////////

C///  FORM MACHINE EPSILON AND THE ABSOLUTE AND RELATIVE PERTURBATIONS.

      CALL BLEPS (EPS)

      ABSOL = SQRT (EPS)
      RELAT = SQRT (EPS)

C///  CLEAR THE ROW SUMS.

      DO 2010 J = 1, N2
         RSCALE(J) = 0.0
2010  CONTINUE

C///  CLEAR THE MATRIX.

      DO 2020 J = 1, ASIZE2
         A2(J) = 0.0
2020  CONTINUE

C///  EVALUATE THE FUNCTION AT THE UNPERTURBED X.

      EVALS = EVALS + 1
C     GO TO 2030 WHEN ROUTE = 1
      ROUTE = 1
      RETURN = .TRUE.
      GO TO 99999
2030  CONTINUE

C///  PLACE THE FUNCTION VALUES IN THE MATRIX.

C     SUBROUTINE BLJA3
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, AR1, AR2, N, Y, YR1, YR2)

      DO 2050 C1 = 1, COLS1
         AC1 = COL1(QFRST, C1)
         AC2 = COL1(QLAST, C1)
         DO 2040 Q = ACOLQ1(QFRST, C1), ACOLQ1(QLAST, C1)
            QA2 = ACOL1(QLOC, Q)
            R1 = ACOL1(QIND, Q)
            YR1 = ROW1(QFRST, R1)
            YR2 = ROW1(QLAST, R1)
            R2 = RMAP(R1)
            AR1 = ROW2(QFRST, R2)
            AR2 = ROW2(QLAST, R2)

      CALL BLJA3
     +  (ERROR, TEXT,
     +   A2(QA2), AC1, AC2, AR1, AR2, N, Y, YR1, YR2)
      IF (ERROR) GO TO 9201

2040     CONTINUE
2050  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (3) BUILD THE JACOBIAN.
C
C///////////////////////////////////////////////////////////////////////

C///  INITIALIZE POINTERS TO THE UNFINISHED COLUMNS WITHIN EACH BLOCK.

      DO 3010 C1 = 1, COLS1
         NEXT(C1) = COL1(QFRST, C1)
3010  CONTINUE

C///  TOP OF THE LOOP OVER GROUPS OF DISJOINT COLUMNS.

3020  CONTINUE

C///  FIND THE GROUP AND PERTURB THE ENTRIES OF X.

      DO 3030 R1 = 1, ROWS1
         OCCUPY(R1) = .FALSE.
3030  CONTINUE

      COUNT = 0
      DO 3060 C1 = 1, COLS1
         IF (NEXT(C1) .LE. COL1(QLAST, C1)) THEN
            TOUCH = .FALSE.
            DO 3040 Q = ACOLQ1(QFRST, C1), ACOLQ1(QLAST, C1)
               R1 = ACOL1(QIND, Q)
               TOUCH = TOUCH .OR. OCCUPY(R1)
3040        CONTINUE
            IF (.NOT. TOUCH) THEN
               COUNT = COUNT + 1
               POINT(COUNT) = C1

               DO 3050 Q = ACOLQ1(QFRST, C1), ACOLQ1(QLAST, C1)
                  R1 = ACOL1(QIND, Q)
                  OCCUPY(R1) = .TRUE.
3050           CONTINUE

               QX = NEXT(C1)
               SAVE(COUNT) = X(QX)
               X(QX) = X(QX) + (ABS (X(QX)) * RELAT + ABSOL)
            END IF
         END IF
3060  CONTINUE

C///  QUIT IF NO COLUMNS WERE PERTURBED.

      IF (COUNT .EQ. 0) GO TO 3100

C///  EVALUATE THE FUNCTION AT THE PERTURBED X.

      EVALS = EVALS + 1
C     GO TO 3070 WHEN ROUTE = 2
      ROUTE = 2
      RETURN = .TRUE.
      GO TO 99999
3070  CONTINUE

C///  TOP OF THE LOOP OVER THE PERTURBED COLUMNS.

      DO 3090 INDEX = 1, COUNT
         C1 = POINT(INDEX)
         AC = NEXT(C1)
         FACTOR = 1.0 / (ABS (SAVE(INDEX)) * RELAT + ABSOL)
         C2 = CMAP(C1)
         AC1 = COL1(QFRST, C1)
         AC2 = COL1(QLAST, C1)

C        INCREMENT THE POINTER
         NEXT(C1) = AC + 1

C        RESTORE X
         X(AC) = SAVE(INDEX)

C///  TOP OF THE LOOP OVER THE ROWS.

      DO 3080 Q = ACOLQ1(QFRST, C1), ACOLQ1(QLAST, C1)
         QA2 = ACOL1(QLOC, Q)
         R1 = ACOL1(QIND, Q)
         YR1 = ROW1(QFRST, R1)
         YR2 = ROW1(QLAST, R1)
         R2 = RMAP(R1)
         AR1 = ROW2(QFRST, R2)
         AR2 = ROW2(QLAST, R2)

C///  EVALUATE THE DIVIDED DIFFERENCES.

C     SUBROUTINE BLJA4
C    +  (ERROR, TEXT,
C    +   A, AC, AC1, AC2, AR1, AR2, FACTOR, N, RSCALE, Y, YR1, YR2)

      CALL BLJA4
     +  (ERROR, TEXT,
     +   A2(QA2), AC, AC1, AC2, AR1, AR2, FACTOR, N, RSCALE(AR1), Y,
     +   YR1, YR2)
      IF (ERROR) GO TO 9301

C///  BOTTOM OF THE LOOPS OVER ROWS AND COLUMS.

3080     CONTINUE
3090  CONTINUE

C///  BOTTOM OF THE LOOP GROUPS.

      GO TO 3020
3100  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (4) EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) WRITE (TEXT, 10001)
     +   EVALS, N, NINT (100.0 * REAL (EVALS) / REAL (N))

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +   (/10X, I10, '  FUNCTION EVALUATIONS'
     +    /10X, I10, '  MATRIX ORDER'
     +    /10X, I10, '  EVALUATIONS AS PERCENT OF ORDER')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   ASIZE2, BLKS1, COLS1, N, N2, ROWS1, ROWS2
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID
      GO TO 99999

9301  IF (0 .LT. TEXT) WRITE (TEXT, 99301) ID
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  ASIZE2'
     +    /10X, I10, '  BLKS1'
     +    /10X, I10, '  COLS1'
     +    /10X, I10, '  N',
     +    /10X, I10, '  N2',
     +    /10X, I10, '  ROWS1'
     +    /10X, I10, '  ROWS2')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  BLJA3 FAILS.')

99301 FORMAT
     +   (/1X, A9, 'ERROR.  BLJA4 FAILS.')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLJA3
     +  (ERROR, TEXT,
     +   A, AC1, AC2, AR1, AR2, N, Y, YR1, YR2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLJA3
C
C>    CALLED BY BLJA2
C
C>    REPLICATE A COLUMN IN A BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER AC, AC1, AC2, AR1, AR2, N, TEXT, YR, YR1, YR2
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, Y

      PARAMETER (ID = ' BLJA3:  ')

      DIMENSION A(AR1 : AR2, AC1 : AC2), Y(N)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (AC1 .LE. AC2 .AND. AR1 .LE. AR2 .AND. 0 .LT. N)
      IF (ERROR) GO TO 9101

C///  CHECK THE PLACEMENT.

      ERROR = .NOT. (1 .LE. YR1 .AND. YR1 .LE. YR2 .AND. YR2 .LE. N)
      IF (ERROR) GO TO 9102

      ERROR = .NOT. (YR2 - YR1 .LE. AR2 - AR1)
      IF (ERROR) GO TO 9103

C///////////////////////////////////////////////////////////////////////
C
C     (2) COPY THE COLUMN.
C
C///////////////////////////////////////////////////////////////////////

      DO 2020 AC = AC1, AC2
         DO 2010 YR = YR1, YR2
            A(AR1 - YR1 + YR, AC) = Y(YR)
2010     CONTINUE
2020  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, AC1, AC2, AR1, AR2, N
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, YR1, YR2, N
      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID, AR1, AR2, YR1, YR2

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2'
     +   //10X, I10, '  N')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  THE COLUMN DOES NOT CONTAIN THE BLOCK TO'
     +    /10X, 'COPY.'
     +   //10X, I10, '  YR1'
     +    /10X, I10, '  YR2'
     +    /10X, I10, '  N')

99103 FORMAT
     +   (/1X, A9, 'ERROR.  THE COPIED BLOCK IS TOO LARGE.'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2'
     +   //10X, I10, '  YR1'
     +    /10X, I10, '  YR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLJA4
     +  (ERROR, TEXT,
     +   A, AC, AC1, AC2, AR1, AR2, FACTOR, N, RSCALE, Y, YR1, YR2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLJA4
C
C>    CALLED BY BLJA2
C
C>    EVALUATE A COLUMN IN A BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER AC, AC1, AC2, AR, AR1, AR2, N, TEXT, YR, YR1, YR2
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, FACTOR, RSCALE, Y

      PARAMETER (ID = ' BLJA4:  ')

      DIMENSION A(AR1 : AR2, AC1 : AC2), RSCALE(AR1 : AR2), Y(N)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (AC1 .LE. AC2 .AND. AR1 .LE. AR2 .AND. 0 .LT. N)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) FORM THE COLUMN AND UPDATE THE SCALE FACTORS.
C
C///////////////////////////////////////////////////////////////////////

      DO 2010 YR = YR1, YR2
         AR = AR1 - YR1 + YR
         A(AR, AC) = (Y(YR) - A(AR, AC)) * FACTOR
         RSCALE(AR) = RSCALE(AR) + ABS (A(AR, AC))
2010  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, AC1, AC2, AR1, AR2, N
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2'
     +   //10X, I10, '  N')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLJAC
     +  (ERROR, TEXT, PRINT,
     +   ILAST, IMAX, ISIZE, IWORK,
     +   LLAST, LMAX, LSIZE, LWORK,
     +   RLAST, RMAX, RSIZE, RWORK,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   N, RETURN, X, Y)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLJAC
C
C>    CALLED BY THE  U S E R
C
C>    EVALUATE A JACOBIAN MATRIX.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   ASIZE2, BLKS1, BLKS2, BORDER, COLS1, COLS2, GLUES, IDATA,
     +   IDSIZ, ILAST, IMAX, ISAVE, ISIZE, IWORK, LLAST, LMAX, LSAVE,
     +   LSIZE, N, N1, N2, PRINT, RDSIZ, RLAST, RMAX, ROWS1, ROWS2,
     +   RSAVE, RSIZE, TEXT
      LOGICAL ERROR, LWORK, RETURN
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   RDATA, RWORK, X, Y

      PARAMETER (ID = ' BLJAC:  ')

      DIMENSION
     +   IDATA(IDSIZ), IWORK(ISIZE), LWORK(LSIZE), RDATA(RDSIZ),
     +   RWORK(RSIZE), X(N), Y(N)

      SAVE

C///////////////////////////////////////////////////////////////////////
C
C     RETURN FROM REVERSE COMMUNICATION.
C
C///////////////////////////////////////////////////////////////////////

      IF (RETURN) GO TO 2010

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  SAVE THE LEVELS OF THE WORK SPACES.

      ISAVE = ILAST
      LSAVE = LLAST
      RSAVE = RLAST

C///  PRINT.

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) WRITE (TEXT, 10001) ID

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. IDSIZ .AND. 0 .LT. ISIZE .AND. 0 .LT. LSIZE
     +   .AND. 0 .LT. N .AND. 0 .LT. RDSIZ .AND. 0 .LT. RSIZE)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) CALL BLJA2.
C
C///////////////////////////////////////////////////////////////////////

C///  UNPACK THE DATA SPACES.

      CALL BLDAT
     +  (ERROR, TEXT,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   ASIZE2, BLKS1, BLKS2, BORDER, COLS1, COLS2, GLUES, QA2, QAC1,
     +   QAC2, QACQ1, QACQ2, QAR2, QARQ2, QBMAP, QCMAP, QCOL1, QCOL2,
     +   QGMAP, QRMAP, QROW1, QROW2, QRPIV, QRSCAL, N1, N2, ROWS1,
     +   ROWS2)
      IF (ERROR) GO TO 9201

      ERROR = .NOT. (N .EQ. N1)
      IF (ERROR) GO TO 9202

C///  RESERVE SOME WORK SPACE.

C     SUBROUTINE RESERV
C    +  (ERROR, TEXT,
C    +   FATAL, QLAST, QMAX, QSIZE,
C    +   NAME, NUMBER, Q)

C     INTEGER NEXT(COLS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QNEXT', COLS1, QNEXT)
      IF (ERROR) GO TO 9203

C     LOGICAL OCCUPY(ROWS1)

      CALL RESERV (ERROR, TEXT, .FALSE., LLAST, LMAX, LSIZE,
     +   'QOCCUP', ROWS1, QOCCUP)
      IF (ERROR) GO TO 9203

C     INTEGER POINT(COLS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QPOINT', COLS1, QPOINT)
      IF (ERROR) GO TO 9203

C     REAL SAVE(COLS1)

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QSAVE', COLS1, QSAVE)
      IF (ERROR) GO TO 9203

      ERROR = .NOT. (ILAST .LE. ISIZE .AND. LLAST .LE. LSIZE .AND.
     +   RLAST .LE. RSIZE)
      IF (ERROR) GO TO 9204

C///  CALL.

C     SUBROUTINE BLJA2
C    +  (ERROR, TEXT, PRINT,
C    +   A2, ACOL1, ACOLQ1, ASIZE2, BLKS1, CMAP, COL1, COLS1, N, N2,
C    +   NEXT, OCCUPY, POINT, RETURN, RMAP, ROW1, ROW2, ROWS1, ROWS2,
C    +   RSCALE, SAVE, X, Y)

2010  CONTINUE

      CALL BLJA2
     +  (ERROR, TEXT, PRINT - 1,
     +   RDATA(QA2), IDATA(QAC1), IDATA(QACQ1), ASIZE2, BLKS1,
     +   IDATA(QCMAP), IDATA(QCOL1), COLS1, N, N2, IWORK(QNEXT),
     +   LWORK(QOCCUP), IWORK(QPOINT), RETURN, IDATA(QRMAP),
     +   IDATA(QROW1), IDATA(QROW2), ROWS1, ROWS2, RDATA(QRSCAL),
     +   RWORK(QSAVE), X, Y)
      IF (ERROR) GO TO 9205

      IF (RETURN) GO TO 99999

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +   (/1X, A, 'EVALUATE A SPARSE BLOCKED JACOBIAN MATRIX.')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99998

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   IDSIZ, ISIZE, LSIZE, N, RDSIZ, RSIZE
      GO TO 99998

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID
      GO TO 99998

9202  IF (0 .LT. TEXT) WRITE (TEXT, 99202) ID, N, N1
      GO TO 99998

9203  IF (0 .LT. TEXT) WRITE (TEXT, 99203) ID
      GO TO 99998

9204  IF (0 .LT. TEXT) WRITE (TEXT, 99204) ID,
     +   0, ISIZE, LSIZE, RSIZE, 0, ILAST, LLAST, RLAST
      GO TO 99998

9205  IF (0 .LT. TEXT) WRITE (TEXT, 99205) ID
      GO TO 99998

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  IDSIZ'
     +    /10X, I10, '  ISIZE'
     +    /10X, I10, '  LSIZE'
     +    /10X, I10, '  N',
     +    /10X, I10, '  RDSIZ'
     +    /10X, I10, '  RSIZE')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  BLDAT FAILS.')

99202 FORMAT
     +   (/1X, A9, 'ERROR.  THE MATRIX ORDER IS INCONSISTENT.'
     +   //10X, I10,  '  ORDER IN ARGUMENT LIST'
     +    /10X, I10,  '  IN DATA FROM BLSYM')

99203 FORMAT
     +   (/1X, A9, 'ERROR.  RESERV FAILS.')

99204 FORMAT
     +   (/1X, A9, 'ERROR.  SOME WORKSPACES ARE TOO SMALL.'
C               123456789_  123456789_  123456789_  123456789_
     +  //25X, ' CHARACTER     INTEGER     LOGICAL        REAL'
C               123456789_12345
     +  //10X, ' PRESENT SIZE', 4I12
     +   /10X, 'REQUIRED SIZE', 4I12
     +  //10X, 'THIS SUBROUTINE WILL NEED NO MORE SPACE.')

99205 FORMAT
     +   (/1X, A9, 'ERROR.  BLJA2 FAILS.')

C///  EXIT.

99998 CONTINUE
      ILAST = ISAVE
      LLAST = LSAVE
      RLAST = RSAVE
99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLMA2
     +  (ERROR, TEXT, PRINT,
     +   ACOL1, ACOLQ1, ACOL2, ACOLQ2, ACTIVE, ASIZE2, BLKS1, BLKS2,
     +   BLOCKS, BMAX, CMAP, COL1, COL2, COLS1, COLS2, CPOINT, DIAG,
     +   FILL, ITEMP, OCCUPY, N1, N2, RINDEX, RMAP, ROW1, ROW2, ROWS1,
     +   ROWS2, SIZE, TOTAL)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLMA2
C
C>    CALLED BY BLMAK
C
C>    MAKE A DATA SPACE FOR A BORDERLESS MATRIX.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   ACOL1, ACOLQ1, ACOL2, ACOLQ2, ASIZE1, ASIZE2, BLKS1, BLKS2,
     +   BLOCK1, BLOCK2, BLOCKS, BMAX, C1, C2, CMAP, COL1, COL2, COLS1,
     +   COLS2, COUNT, COUNT1, COUNT2, CPOINT, DIAG, FILL, FIRST, FMAX,
     +   ITEMP, J, LAST, LOC, NEXT, N1, N2, PRINT, R1, R2, RINDEX, RMAP,
     +   ROW1, ROW2, ROWS1, ROWS2, SIZE, TEXT, TOTAL
      LOGICAL ACTIVE, ERROR, FOUND, OCCUPY

      PARAMETER (ID = ' BLMA2:  ')

      PARAMETER (QFRST = 1, QLAST = 2, QFRSTD = 3, QLASTD = 4)
      PARAMETER (QIND = 1, QLOC = 2)

      DIMENSION
     +   ACOL1(2, BLKS1), ACOLQ1(2, COLS1), ACOL2(2, BMAX),
     +   ACOLQ2(4, COLS2), ACTIVE(ROWS2), CMAP(COLS1), COL1(2, COLS1),
     +   COL2(2, COLS2), CPOINT(BLOCKS), DIAG(2, ROWS2), ITEMP(ROWS2),
     +   OCCUPY(ROWS1), RINDEX(TOTAL), RMAP(ROWS1), ROW1(2, ROWS1),
     +   ROW2(2, ROWS2), SIZE(BLOCKS)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. BLKS1 .AND. 0 .LT. BLOCKS .AND. 0 .LT. BMAX
     +   .AND. 0 .LT. COLS1 .AND. 0 .LT. COLS2 .AND. 0 .LT. ROWS1 .AND.
     +   0 .LT. ROWS2 .AND. 0 .LT. TOTAL)
      IF (ERROR) GO TO 9101

      ERROR = .NOT. (BLKS1 .EQ. TOTAL .AND. BLOCKS .EQ. COLS1 .AND.
     +   BLOCKS .EQ. COLS2 .AND. BLOCKS .EQ. ROWS1 .AND.
     +   BLOCKS .EQ. ROWS2)
      IF (ERROR) GO TO 9102

C///  CHECK THE FACTORIZATION FILL LEVEL.

      ERROR = .NOT. (0 .LE. FILL)
      IF (ERROR) GO TO 9103

C///////////////////////////////////////////////////////////////////////
C
C     (2) MAKE THE SIMPLE, DIMENSIONAL POINTERS.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE BLOCK SIZES.

      COUNT = 0
      DO 2010 C1 = 1, COLS1
         IF (0 .LT. SIZE(C1)) COUNT = COUNT + 1
2010  CONTINUE
      ERROR = .NOT. (COUNT .EQ. COLS1)
      IF (ERROR) GO TO 9201

C///  CMAP(COLS1) : A2 COL BLOCK CONTAINING AN A1 BLOCK

      DO 2020 C1 = 1, COLS1
         CMAP(C1) = C1
2020  CONTINUE

C///  COL1(2, COLS1) : POINTERS INTO A1 COLS

      LOC = 0
      DO 2030 C1 = 1, COLS1
         COL1(QFRST, C1) = LOC + 1
         LOC = LOC + SIZE(C1)
         COL1(QLAST, C1) = LOC
2030  CONTINUE

C///  COL2(2, COLS2) : POINTERS INTO A2 COLS

      DO 2040 C2 = 1, COLS2
         COL2(QFRST, C2) = COL1(QFRST, C2)
         COL2(QLAST, C2) = COL1(QLAST, C2)
2040  CONTINUE

C///  RMAP(ROWS1) : A2 ROW BLOCK CONTAINING AN A1 BLOCK

      DO 2050 R1 = 1, ROWS1
         RMAP(R1) = R1
2050  CONTINUE

C///  ROW1(2, ROWS1) : POINTERS INTO A1 ROWS

      LOC = 0
      DO 2060 R1 = 1, ROWS1
         ROW1(QFRST, R1) = LOC + 1
         LOC = LOC + SIZE(R1)
         ROW1(QLAST, R1) = LOC
2060  CONTINUE

C///  MAKE THE ORDERS.

      N1 = LOC
      N2 = N1

C///  ROW2(2, ROWS2) : POINTERS INTO A2 ROWS

      DO 2070 R2 = 1, ROWS2
         ROW2(QFRST, R2) = ROW1(QFRST, R2)
         ROW2(QLAST, R2) = ROW1(QLAST, R2)
2070  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (3) MAKE THE A1 COLUMN STRUCTURE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE POINTERS.

      ERROR = .NOT. (1 .EQ. CPOINT(1))
      IF (ERROR) GO TO 9301

      COUNT = 0
      DO 3010 C1 = 1, COLS1
         IF (1 .LE. CPOINT(C1) .AND. CPOINT(C1) .LE. TOTAL)
     +      COUNT = COUNT + 1
3010  CONTINUE
      ERROR = .NOT. (COUNT .EQ. COLS1)
      IF (ERROR) GO TO 9302

      COUNT1 = 0
      COUNT2 = 0
      DO 3020 C1 = 1, COLS1
         IF (C1 .LT. COLS1) THEN
            NEXT = CPOINT(C1 + 1)
         ELSE
            NEXT = TOTAL + 1
         END IF
         IF (CPOINT(C1) .GT. NEXT) COUNT1 = COUNT1 + 1
         IF (CPOINT(C1) .EQ. NEXT) COUNT2 = COUNT2 + 1
3020  CONTINUE

      ERROR = .NOT. (COUNT1 .EQ. 0)
      IF (ERROR) GO TO 9303

      ERROR = .NOT. (COUNT2 .EQ. 0)
      IF (ERROR) GO TO 9304

C///  TOP OF THE LOOP OVER THE COLUMNS.

      DO 3030 R1 = 1, ROWS1
         OCCUPY(R1) = .FALSE.
3030  CONTINUE

      COUNT = 0
      DO 3060 C1 = 1, COLS1

C///  FIND THE BEGINNING AND END OF THE LIST.

      FIRST = CPOINT(C1)
      IF (C1 .LT. COLS1) THEN
         LAST = CPOINT(C1 + 1) - 1
      ELSE
         LAST = TOTAL
      END IF

C///  MARK THE OCCUPIED ROWS.

      DO 3040 J = FIRST, LAST
         R1 = RINDEX(J)
         OCCUPY(R1) = .TRUE.
3040  CONTINUE

C///  SAVE THE OCCUPIED ROWS.

      ACOLQ1(QFRST, C1) = COUNT + 1
      DO 3050 R1 = 1, ROWS1
         IF (OCCUPY(R1)) THEN
            COUNT = COUNT + 1
            ACOL1(QIND, COUNT) = R1
            OCCUPY(R1) = .FALSE.
         END IF
3050  CONTINUE
      ACOLQ1(QLAST, C1) = COUNT

C///  CHECK FOR DUPLICATE ROWS.

      ERROR = .NOT. (LAST - FIRST .EQ.
     +   ACOLQ1(QLAST, C1) - ACOLQ1(QFRST, C1))
      IF (ERROR) GO TO 9305

C///  BOTTOM OF THE LOOP OVER THE COLUMNS.

3060  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (4) PERFORM SYMBOLIC FACTORIZATION AND MAKE MOST OF THE A2 COLUMN
C     STRUCTURE.
C
C///////////////////////////////////////////////////////////////////////

C///  INITIALIZE THE SCRATCH ARRAYS.

      DO 4010 R2 = 1, ROWS2
         ACTIVE(R2) = .FALSE.
         ITEMP(R2) = FILL + 1
4010  CONTINUE

C///  TOP OF THE LOOP OVER THE COLUMNS.

      BLKS2 = 0
      FMAX = 0
      DO 4030 C2 = 1, COLS2

C///  COPY THE A1 SPARSITY PATTERN.

      C1 = C2
      DO 4020 LOC = ACOLQ1(QFRST, C1), ACOLQ1(QLAST, C1)
         R1 = ACOL1(QIND, LOC)
         R2 = RMAP(R1)
         ITEMP(R2) = 0
4020  CONTINUE

C///  ELIMINATE AND BUILD.

C     SUBROUTINE BLMA4
C    +  (ERROR, TEXT,
C    +   ACOL2, ACOLQ2, ACTIVE, BLKS2, BMAX, COL, COL2, COLS2, DIAG,
C    +   FILL, FMAX, LEVEL, ROW2, ROWS2)

      CALL BLMA4
     +  (ERROR, TEXT,
     +   ACOL2, ACOLQ2, ACTIVE, BLKS2, BMAX, C2, COL2, COLS2, DIAG,
     +   FILL, FMAX, ITEMP, ROW2, ROWS2)
      IF (ERROR) GO TO 9401

C///  BOTTOM OF THE LOOP OVER THE COLUMNS.

4030  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (5) COMPLETE THE A1 AND A2 COLUMN POINTER STRUCTURES.
C
C///////////////////////////////////////////////////////////////////////

C///  MAKE THE BLOCK ASSIGNMENTS FOR A2.

      ASIZE2 = 0
      DO 5020 C2 = 1, COLS2
         DO 5010 BLOCK2 = ACOLQ2(QFRST, C2), ACOLQ2(QLAST, C2)
            ACOL2(QLOC, BLOCK2) = ASIZE2 + 1
            R2 = ACOL2(QIND, BLOCK2)
            ASIZE2 = ASIZE2
     +         + (COL2(QLAST, C2) - COL2(QFRST, C2) + 1)
     +         * (ROW2(QLAST, R2) - ROW2(QFRST, R2) + 1)
5010     CONTINUE
5020  CONTINUE

C///  COPY THE BLOCK LOCATIONS FOR A1.

      ASIZE1 = 0
      DO 5060 C1 = 1, COLS1
         C2 = CMAP(C1)

         DO 5050 BLOCK1 = ACOLQ1(QFRST, C1), ACOLQ1(QLAST, C1)
            R1 = ACOL1(QIND, BLOCK1)
            R2 = RMAP(R1)

            ASIZE1 = ASIZE1
     +         + (COL1(QLAST, C1) - COL1(QFRST, C1) + 1)
     +         * (ROW1(QLAST, R1) - ROW1(QFRST, R1) + 1)

            FOUND = .FALSE.
            DO 5030 BLOCK2 = ACOLQ2(QFRST, C2), ACOLQ2(QLAST, C2)
               IF (R2 .EQ. ACOL2(QIND, BLOCK2)) THEN
                  ACOL1(QLOC, BLOCK1) = ACOL2(QLOC, BLOCK2)
                  FOUND = .TRUE.
                  GO TO 5040
               END IF
5030        CONTINUE
5040        CONTINUE
            ERROR = .NOT. FOUND
            IF (ERROR) GO TO 9501
5050     CONTINUE
5060  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     PRINT.
C
C///////////////////////////////////////////////////////////////////////

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) THEN
         WRITE (TEXT, 10001)
     +      N1, 0, BLKS1, ASIZE1, N2, FMAX, BLKS2, ASIZE2
         IF (FMAX .LE. FILL) WRITE (TEXT, 10002)
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
C               123456789_  123456789_  123456789_  123456789_
     +  (/10X, '                  FILL       TOTAL      MATRIX'
     +   /10X, '     ORDER       LEVEL      BLOCKS       SPACE'
     +   /10X, 4(I10, 2X), 'ORIGINAL MATRIX'
     +   /10X, 4(I10, 2X), 'FACTORED MATRIX')

10002 FORMAT
     +  (/10X, 'THE FACTORIZATION COMPLETES AT THIS LEVEL.')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   BLKS1, BLOCKS, BMAX, COLS1, COLS2, ROWS1, ROWS2, TOTAL
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID,
     +   BLKS1, TOTAL, BLOCKS, COLS1, COLS2, ROWS1, ROWS2
      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID, FILL
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID, COLS1, COLS1 - COUNT
      GO TO 99999

9301  IF (0 .LT. TEXT) WRITE (TEXT, 99301) ID, CPOINT(1)
      GO TO 99999

9302  IF (0 .LT. TEXT) WRITE (TEXT, 99302) ID,
     +   TOTAL, COLS1, COLS1 - COUNT
      GO TO 99999

9303  IF (0 .LT. TEXT) WRITE (TEXT, 99303) ID, COLS1, COUNT1
      GO TO 99999

9304  IF (0 .LT. TEXT) WRITE (TEXT, 99304) ID, COLS1, COUNT2
      GO TO 99999

9305  IF (0 .LT. TEXT) WRITE (TEXT, 99305) ID,
     +   COLS1, C1, LAST - FIRST + 1,
     +   ACOLQ1(QLAST, C1) - ACOLQ1(QFRST, C1) + 1,
     +   (J, RINDEX(J), J = FIRST, LAST)
      GO TO 99999

9401  IF (0 .LT. TEXT) WRITE (TEXT, 99401) ID
      GO TO 99999

9501  IF (0 .LT. TEXT) WRITE (TEXT, 99501) ID, C1, C2, R1, R2
      GO TO 99999

99101 FORMAT
     +  (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +  //10X, I10, '  BLKS1'
     +   /10X, I10, '  BLOCKS'
     +   /10X, I10, '  BMAX'
     +   /10X, I10, '  COLS1'
     +   /10X, I10, '  COLS2'
     +   /10X, I10, '  ROWS1'
     +   /10X, I10, '  ROWS2'
     +   /10X, I10, '  TOTAL')

99102 FORMAT
     +  (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +  //10X, I10, '  BLKS1'
     +   /10X, I10, '  TOTAL'
     +  //10X, I10, '  BLOCKS'
     +   /10X, I10, '  COLS1'
     +   /10X, I10, '  COLS2'
     +   /10X, I10, '  ROWS1'
     +   /10X, I10, '  ROWS2')

99103 FORMAT
     +  (/1X, A9, 'ERROR.  THE FILL LEVEL MUST BE ZERO OR POSITIVE.'
     +  //10X, I10, '  FILL')

99201 FORMAT
     +  (/1X, A9, 'ERROR.  ALL VECTOR BLOCK SIZES MUST BE POSITIVE.'
     +  //10X, I10, '  VECTOR BLOCKS'
     +   /10X, I10, '  WITHOUT POSITIVE SIZE')

99301 FORMAT
     +  (/1X, A9, 'ERROR.  THE FIRST COLUMN POINTER MUST BE 1.'
     +  //10X, I10, '  CPOINT(1)')

99302 FORMAT
     +  (/1X, A9, 'ERROR.  COLUMN POINTERS TO ROW INDICES MUST RANGE'
     +   /10X, 'FROM 1 TO TOTAL MATRIX BLOCKS.'
     +  //10X, I10, '  TOTAL MATRIX BLOCKS'
     +   /10X, I10, '  COLUMNS'
     +   /10X, I10, '  POINTERS OUT OF RANGE')

99303 FORMAT
     +  (/1X, A9, 'ERROR.  COLUMN POINTERS TO ROW INDICES MUST'
     +   /10X, 'INCREASE WITH THE COLUMN NUMBER.'
     +  //10X, I10, '  COLUMNS'
     +   /10X, I10, '  POINTERS EXCEEDING THEIR SUCCESSORS')

99304 FORMAT
     +  (/1X, A9, 'ERROR.  SOME COLUMNS HAVE NO BLOCKS.'
     +  //10X, I10, '  COLUMNS'
     +   /10X, I10, '  COLUMN POINTERS MATCHING THEIR SUCCESSORS')

99305 FORMAT
     +  (/1X, A9, 'ERROR.  A COLUMN HAS DUPLICATE ROW INDICES.'
     +  //10X, I10, '  TOTAL COLUMNS'
     +   /10X, I10, '  INDEX OF THIS COLUMN'
     +  //10X, I10, '  ROW INDICES FOR THIS COLUMN'
     +   /10X, I10, '  DISTINCT INDICES'
C               123456789_  123456789_
     +  //10X, '     COUNT       INDEX'
     +  /(10X, I10, 2X, I10))

99401 FORMAT
     +  (/1X, A9, 'ERROR.  BLMA4 FAILS.')

99501 FORMAT
     +  (/1X, A9, 'ERROR.  A BLOCK FROM A1 IS NOT IN A2.'
     +  //10X, I10, '  C1'
     +   /10X, I10, '  C2'
     +  //10X, I10, '  R1'
     +   /10X, I10, '  R2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLMA3
     +  (ERROR, TEXT, PRINT,
     +   ACOL1, ACOLQ1, ACOL2, ACOLQ2, ACTIVE, ASIZE2, BLKS1, BLKS2,
     +   BLOCKS, BMAP, BMAX, BORDER, CMAP, COL1, COL2, COLS1, COLS2,
     +   CPOINT, DIAG, FILL, GLUES, GMAP, ITEMP, N1, N2, OCCUPY, RINDEX,
     +   RMAP, ROW1, ROW2, ROWS1, ROWS2, SIZE, TOTAL)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLMA3
C
C>    CALLED BY BLMAK
C
C>    MAKE A DATA SPACE FOR A BORDERED MATRIX.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   ACOL1, ACOLQ1, ACOL2, ACOLQ2, ASIZE1, ASIZE2, BLKS1, BLKS2,
     +   BLOCK1, BLOCK2, BLOCKS, BMAP, BMAX, BORDER, C1, C2, CMAP, COL1,
     +   COL2, COLS1, COLS2, COUNT, COUNT1, COUNT2, CPOINT, DIAG, FILL,
     +   FIRST, FMAX, GLUES, GMAP, ITEMP, J, LAST, LOC, N1, N2, NEXT,
     +   PRINT, R1, R2, RINDEX, RMAP, ROW1, ROW2, ROWS1, ROWS2, SIZE,
     +   TEXT, TOTAL
      LOGICAL ACTIVE, ERROR, FOUND, OCCUPY

      PARAMETER (ID = ' BLMA3:  ')

      PARAMETER (QFRST = 1, QLAST = 2, QFRSTD = 3, QLASTD = 4)
      PARAMETER (QIND = 1, QLOC = 2)
      PARAMETER (QCOL = 1, QNEG = 2, QPOS = 3)

      DIMENSION
     +   ACOL1(2, BLKS1), ACOLQ1(2, COLS1), ACOL2(2, BMAX),
     +   ACOLQ2(4, COLS2), ACTIVE(ROWS2), BMAP(COLS1 + 1), CMAP(COLS1),
     +   COL1(2, COLS1), COL2(2, COLS2), CPOINT(BLOCKS), DIAG(2, ROWS2),
     +   GMAP(3, 0 : GLUES), ITEMP(ROWS2), OCCUPY(ROWS1), RINDEX(TOTAL),
     +   RMAP(ROWS1), ROW1(2, ROWS1), ROW2(2, ROWS2), SIZE(BLOCKS)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. BLKS1 .AND. 0 .LT. BLOCKS .AND. 0 .LT. BMAX
     +   .AND. 0 .LT. BORDER .AND. 0 .LT. COLS1 .AND. 0 .LT. COLS2 .AND.
     +   0 .LE. GLUES .AND. 0 .LT. ROWS1 .AND. 0 .LT. ROWS2 .AND.
     +   0 .LT. TOTAL)
      IF (ERROR) GO TO 9101

      ERROR = .NOT. (BLKS1 .EQ. TOTAL .AND. BLOCKS .EQ. COLS1 .AND.
     +   BLOCKS .EQ. ROWS1 .AND. BLOCKS .EQ. ROWS2 .AND.
     +   2 * BLOCKS .EQ. COLS2 .AND. BLOCKS - 1 .EQ. GLUES)
      IF (ERROR) GO TO 9102

C///  CHECK THE FACTORIZATION FILL LEVEL.

      ERROR = .NOT. (0 .LE. FILL)
      IF (ERROR) GO TO 9103

C///////////////////////////////////////////////////////////////////////
C
C     (2) MAKE THE SIMPLE, DIMENSIONAL POINTERS.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE BLOCK SIZES.

      COUNT = 0
      DO 2010 C1 = 1, COLS1
         IF (0 .LT. SIZE(C1)) COUNT = COUNT + 1
2010  CONTINUE
      ERROR = .NOT. (COUNT .EQ. COLS1)
      IF (ERROR) GO TO 9201

C///  CMAP(COLS1) : A2 COL BLOCK CONTAINING AN A1 BLOCK

      DO 2020 C1 = 1, COLS1
         CMAP(C1) = 2 * C1 - 1
2020  CONTINUE

C///  COL1(2, COLS1) : POINTERS INTO A1 COLS

      LOC = 0
      DO 2030 C1 = 1, COLS1
         COL1(QFRST, C1) = LOC + 1
         LOC = LOC + SIZE(C1)
         COL1(QLAST, C1) = LOC
2030  CONTINUE

C///  COL2(2, COLS2) : POINTERS INTO A2 COLS

      C2 = 0
      LOC = 0
      DO 2040 C1 = 1, COLS1
         C2 = C2 + 1
         COL2(QFRST, C2) = LOC + 1
         LOC = LOC + SIZE(C1)
         COL2(QLAST, C2) = LOC

         C2 = C2 + 1
         COL2(QFRST, C2) = LOC + 1
         LOC = LOC + BORDER
         COL2(QLAST, C2) = LOC
2040  CONTINUE

C///  RMAP(ROWS1) : A2 ROW BLOCK CONTAINING AN A1 BLOCK

      DO 2050 R1 = 1, ROWS1
         RMAP(R1) = R1
2050  CONTINUE

C///  ROW1(2, ROWS1) : POINTERS INTO A1 ROWS

      LOC = 0
      DO 2060 R1 = 1, ROWS1
         ROW1(QFRST, R1) = LOC + 1
         LOC = LOC + SIZE(R1)
         ROW1(QLAST, R1) = LOC
2060  CONTINUE
      N1 = LOC

C///  ROW2(2, ROWS2) : POINTERS INTO A2 ROWS

      LOC = 0
      DO 2070 R1 = 1, ROWS1
         ROW2(QFRST, R1) = LOC + 1
         LOC = LOC + SIZE(R1) + BORDER
         ROW2(QLAST, R1) = LOC
2070  CONTINUE
      N2 = LOC

C///////////////////////////////////////////////////////////////////////
C
C     (3) MAKE THE A1 COLUMN STRUCTURE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE POINTERS.

      ERROR = .NOT. (1 .EQ. CPOINT(1))
      IF (ERROR) GO TO 9301

      COUNT = 0
      DO 3010 C1 = 1, COLS1
         IF (1 .LE. CPOINT(C1) .AND. CPOINT(C1) .LE. TOTAL)
     +      COUNT = COUNT + 1
3010  CONTINUE
      ERROR = .NOT. (COUNT .EQ. COLS1)
      IF (ERROR) GO TO 9302

      COUNT1 = 0
      COUNT2 = 0
      DO 3020 C1 = 1, COLS1
         IF (C1 .LT. COLS1) THEN
            NEXT = CPOINT(C1 + 1)
         ELSE
            NEXT = TOTAL + 1
         END IF
         IF (CPOINT(C1) .GT. NEXT) COUNT1 = COUNT1 + 1
         IF (CPOINT(C1) .EQ. NEXT) COUNT2 = COUNT2 + 1
3020  CONTINUE

      ERROR = .NOT. (COUNT1 .EQ. 0)
      IF (ERROR) GO TO 9303

      ERROR = .NOT. (COUNT2 .EQ. 0)
      IF (ERROR) GO TO 9304

C///  TOP OF THE LOOP OVER THE COLUMNS.

      DO 3030 R1 = 1, ROWS1
         OCCUPY(R1) = .FALSE.
3030  CONTINUE

      COUNT = 0
      DO 3060 C1 = 1, COLS1

C///  FIND THE BEGINNING AND END OF THE LIST.

      FIRST = CPOINT(C1)
      IF (C1 .LT. COLS1) THEN
         LAST = CPOINT(C1 + 1) - 1
      ELSE
         LAST = TOTAL
      END IF

C///  MARK THE OCCUPIED ROWS.

      DO 3040 J = FIRST, LAST
         R1 = RINDEX(J)
         OCCUPY(R1) = .TRUE.
3040  CONTINUE

C///  SAVE THE OCCUPIED ROWS.

      ACOLQ1(QFRST, C1) = COUNT + 1
      DO 3050 R1 = 1, ROWS1
         IF (OCCUPY(R1)) THEN
            COUNT = COUNT + 1
            ACOL1(QIND, COUNT) = R1
            OCCUPY(R1) = .FALSE.
         END IF
3050  CONTINUE
      ACOLQ1(QLAST, C1) = COUNT

C///  CHECK FOR DUPLICATE ROWS.

      ERROR = .NOT. (LAST - FIRST .EQ.
     +   ACOLQ1(QLAST, C1) - ACOLQ1(QFRST, C1))
      IF (ERROR) GO TO 9305

C///  BOTTOM OF THE LOOP OVER THE COLUMNS.

3060  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (4) PERFORM SYMBOLIC FACTORIZATION AND MAKE MOST OF THE A2 COLUMN
C     STRUCTURE.
C
C///////////////////////////////////////////////////////////////////////

C///  INITIALIZE THE SCRATCH ARRAYS.

      DO 4010 R2 = 1, ROWS2
         ACTIVE(R2) = .FALSE.
         ITEMP(R2) = FILL + 1
4010  CONTINUE

C///  TOP OF THE LOOP OVER THE COLUMNS.

      BLKS2 = 0
      C2 = 0
      FMAX = 0
      DO 4040 C1 = 1, COLS1

C///  COPY THE A1 SPARSITY PATTERN.

      DO 4020 LOC = ACOLQ1(QFRST, C1), ACOLQ1(QLAST, C1)
         R1 = ACOL1(QIND, LOC)
         R2 = RMAP(R1)
         ITEMP(R2) = 0
4020  CONTINUE

C///  ELIMINATE AND BUILD.

      C2 = C2 + 1

C     SUBROUTINE BLMA4
C    +  (ERROR, TEXT,
C    +   ACOL2, ACOLQ2, ACTIVE, BLKS2, BMAX, COL, COL2, COLS2, DIAG,
C    +   FILL, FMAX, LEVEL, ROW2, ROWS2)

      CALL BLMA4
     +  (ERROR, TEXT,
     +   ACOL2, ACOLQ2, ACTIVE, BLKS2, BMAX, C2, COL2, COLS2, DIAG,
     +   FILL, FMAX, ITEMP, ROW2, ROWS2)
      IF (ERROR) GO TO 9401

C///  INSERT THE PATTERN.

      IF (C1 .LT. COLS1) THEN
C        GLUE
         QD = ACOLQ2(QFRSTD, C2)
         R2 = ACOL2(QIND, QD)
         ITEMP(R2) = 0
         ITEMP(R2 + 1) = 1
      ELSE
C        BORDER COLUMN
         DO 4030 R1 = 1, ROWS1
            R2 = RMAP(R1)
            ITEMP(R2) = 0
4030     CONTINUE
      END IF

C///  ELIMINATE AND BUILD.

      C2 = C2 + 1

C     SUBROUTINE BLMA4
C    +  (ERROR, TEXT,
C    +   ACOL2, ACOLQ2, ACTIVE, BLKS2, BMAX, COL, COL2, COLS2, DIAG,
C    +   FILL, FMAX, LEVEL, ROW2, ROWS2)

      CALL BLMA4
     +  (ERROR, TEXT,
     +   ACOL2, ACOLQ2, ACTIVE, BLKS2, BMAX, C2, COL2, COLS2, DIAG,
     +   FMAX, FILL, ITEMP, ROW2, ROWS2)
      IF (ERROR) GO TO 9401

C///  BOTTOM OF THE LOOP OVER THE COLUMNS.

4040  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (5) COMPLETE THE A1 AND A2 COLUMN POINTER STRUCTURES.
C
C///////////////////////////////////////////////////////////////////////

C///  MAKE THE BLOCK ASSIGNMENTS FOR A2.

      ASIZE2 = 0
      DO 5020 C2 = 1, COLS2
         DO 5010 BLOCK2 = ACOLQ2(QFRST, C2), ACOLQ2(QLAST, C2)
            ACOL2(QLOC, BLOCK2) = ASIZE2 + 1
            R2 = ACOL2(QIND, BLOCK2)
            ASIZE2 = ASIZE2
     +         + (COL2(QLAST, C2) - COL2(QFRST, C2) + 1)
     +         * (ROW2(QLAST, R2) - ROW2(QFRST, R2) + 1)
5010     CONTINUE
5020  CONTINUE

C///  COPY THE BLOCK LOCATIONS FOR A1.

      ASIZE1 = 0
      DO 5060 C1 = 1, COLS1
         C2 = CMAP(C1)

         DO 5050 BLOCK1 = ACOLQ1(QFRST, C1), ACOLQ1(QLAST, C1)
            R1 = ACOL1(QIND, BLOCK1)
            R2 = RMAP(R1)

            ASIZE1 = ASIZE1
     +         + (COL1(QLAST, C1) - COL1(QFRST, C1) + 1)
     +         * (ROW1(QLAST, R1) - ROW1(QFRST, R1) + 1)

            FOUND = .FALSE.
            DO 5030 BLOCK2 = ACOLQ2(QFRST, C2), ACOLQ2(QLAST, C2)
               IF (R2 .EQ. ACOL2(QIND, BLOCK2)) THEN
                  ACOL1(QLOC, BLOCK1) = ACOL2(QLOC, BLOCK2)
                  FOUND = .TRUE.
                  GO TO 5040
               END IF
5030        CONTINUE
5040        CONTINUE
            ERROR = .NOT. FOUND
            IF (ERROR) GO TO 9501
5050     CONTINUE
5060  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (6) MAKE THE GLUE POINTERS.
C
C///////////////////////////////////////////////////////////////////////

C///  BMAP(COLS1 + 1) : LOC IN ACOL2 OF BLOCKS CONTIANING
C///                    BLOCKS FROM THE BORDER ROW

      DO 6010 C1 = 1, COLS1 + 1
         IF (C1 .LE. COLS1) THEN
            C2 = CMAP(C1)
         ELSE
            C2 = COLS2
         END IF
         BMAP(C1) = ACOLQ2(QFRSTD, C2)
6010  CONTINUE

C///  GMAP(3, 0 : GLUES) : POINTERS TO GLUE BLOCKS.

      DO 6020 C1 = 1, COLS1 - 1
         C2 = CMAP(C1) + 1
         GMAP(QCOL, C1) = C2
         GMAP(QNEG, C1) = ACOLQ2(QFRSTD, C2)
         GMAP(QPOS, C1) = ACOLQ2(QFRSTD, C2) + 1
6020  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     PRINT.
C
C///////////////////////////////////////////////////////////////////////

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) THEN
         WRITE (TEXT, 10001)
     +      N1, 0, BLKS1, ASIZE1,
     +      N1 + BORDER, 0, BLKS1 + 2 * BLOCKS + 1,
     +         ASIZE1 + 2 * BORDER * N1 + BORDER ** 2,
     +      N2, 0, BLKS1 + 3 * BLOCKS - 2,
     +      N2, FMAX, BLKS2, ASIZE2
         IF (FMAX .LE. FILL) WRITE (TEXT, 10002)
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
C               123456789_  123456789_  123456789_  123456789_
     +  (/10X, '                  FILL       TOTAL      MATRIX'
     +   /10X, '     ORDER       LEVEL      BLOCKS       SPACE'
     +   /10X, 4(I10, 2X), 'ORIGINAL MATRIX'
     +   /10X, 4(I10, 2X), 'BORDERED MATRIX'
     +   /10X, 3(I10, 2X), 12X, 'STRETCHED MATRIX'
     +   /10X, 4(I10, 2X), 'FACTORED MATRIX')

10002 FORMAT
     +  (/10X, 'THE FACTORIZATION COMPLETES AT THIS LEVEL.')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   BLKS1, BLOCKS, BMAX, BORDER, COLS1, COLS2, GLUES, ROWS1, ROWS2,
     +   TOTAL
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID,
     +   BLKS1, TOTAL, BLOCKS, COLS1, ROWS1, ROWS2, 2 * BLOCKS, COLS2,
     +   BLOCKS - 1, GLUES
      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID, FILL
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID, COLS1, COLS1 - COUNT
      GO TO 99999

9301  IF (0 .LT. TEXT) WRITE (TEXT, 99301) ID, CPOINT(1)
      GO TO 99999

9302  IF (0 .LT. TEXT) WRITE (TEXT, 99302) ID,
     +   TOTAL, COLS1, COLS1 - COUNT
      GO TO 99999

9303  IF (0 .LT. TEXT) WRITE (TEXT, 99303) ID, COLS1, COUNT1
      GO TO 99999

9304  IF (0 .LT. TEXT) WRITE (TEXT, 99304) ID, COLS1, COUNT2
      GO TO 99999

9305  IF (0 .LT. TEXT) WRITE (TEXT, 99305) ID,
     +   COLS1, C1, LAST - FIRST + 1,
     +   ACOLQ1(QLAST, C1) - ACOLQ1(QFRST, C1) + 1,
     +   (J, RINDEX(J), J = FIRST, LAST)
      GO TO 99999

9401  IF (0 .LT. TEXT) WRITE (TEXT, 99401) ID
      GO TO 99999

9501  IF (0 .LT. TEXT) WRITE (TEXT, 99501) ID, C1, C2, R1, R2
      GO TO 99999

99101 FORMAT
     +  (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +  //10X, I10, '  BLKS1'
     +   /10X, I10, '  BLOCKS'
     +   /10X, I10, '  BMAX'
     +   /10X, I10, '  BORDER'
     +   /10X, I10, '  COLS1'
     +   /10X, I10, '  COLS2'
     +   /10X, I10, '  GLUES (MAY BE ZERO)'
     +   /10X, I10, '  ROWS1'
     +   /10X, I10, '  ROWS2'
     +   /10X, I10, '  TOTAL')

99102 FORMAT
     +  (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +  //10X, I10, '  BLKS1'
     +   /10X, I10, '  TOTAL'
     +  //10X, I10, '  BLOCKS'
     +   /10X, I10, '  COLS1'
     +   /10X, I10, '  ROWS1'
     +   /10X, I10, '  ROWS2'
     +  //10X, I10, '  2 * BLOCKS'
     +   /10X, I10, '  COLS2'
     +  //10X, I10, '  BLOCKS - 1'
     +   /10X, I10, '  GLUES')

99103 FORMAT
     +  (/1X, A9, 'ERROR.  THE FILL LEVEL MUST BE ZERO OR POSITIVE.'
     +  //10X, I10, '  FILL')

99201 FORMAT
     +  (/1X, A9, 'ERROR.  ALL VECTOR BLOCK SIZES MUST BE POSITIVE.'
     +  //10X, I10, '  VECTOR BLOCKS'
     +   /10X, I10, '  WITHOUT POSITIVE SIZE')

99301 FORMAT
     +  (/1X, A9, 'ERROR.  THE FIRST COLUMN POINTER MUST BE 1.'
     +  //10X, I10, '  CPOINT(1)')

99302 FORMAT
     +  (/1X, A9, 'ERROR.  COLUMN POINTERS TO ROW INDICES MUST RANGE'
     +   /10X, 'FROM 1 TO TOTAL MATRIX BLOCKS.'
     +  //10X, I10, '  TOTAL MATRIX BLOCKS'
     +   /10X, I10, '  COLUMNS'
     +   /10X, I10, '  POINTERS OUT OF RANGE')

99303 FORMAT
     +  (/1X, A9, 'ERROR.  COLUMN POINTERS TO ROW INDICES MUST'
     +   /10X, 'INCREASE WITH THE COLUMN NUMBER.'
     +  //10X, I10, '  COLUMNS'
     +   /10X, I10, '  POINTERS EXCEEDING THEIR SUCCESSORS')

99304 FORMAT
     +  (/1X, A9, 'ERROR.  SOME COLUMNS HAVE NO BLOCKS.'
     +  //10X, I10, '  COLUMNS'
     +   /10X, I10, '  COLUMN POINTERS MATCHING THEIR SUCCESSORS')

99305 FORMAT
     +  (/1X, A9, 'ERROR.  A COLUMN HAS DUPLICATE ROW INDICES.'
     +  //10X, I10, '  TOTAL COLUMNS'
     +   /10X, I10, '  INDEX OF THIS COLUMN'
     +  //10X, I10, '  ROW INDICES FOR THIS COLUMN'
     +   /10X, I10, '  DISTINCT INDICES'
C               123456789_  123456789_
     +  //10X, '     COUNT       INDEX'
     +  /(10X, I10, 2X, I10))

99401 FORMAT
     +  (/1X, A9, 'ERROR.  BLMA4 FAILS.')

99501 FORMAT
     +  (/1X, A9, 'ERROR.  A BLOCK FROM A1 IS NOT IN A2.'
     +  //10X, I10, '  C1'
     +   /10X, I10, '  C2'
     +  //10X, I10, '  R1'
     +   /10X, I10, '  R2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLMA4
     +  (ERROR, TEXT,
     +   ACOL2, ACOLQ2, ACTIVE, BLKS2, BMAX, COL, COL2, COLS2, DIAG,
     +   FILL, FMAX, LEVEL, ROW2, ROWS2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLMA4
C
C>    CALLED BY BLMA2, BLMA3
C
C>    PERFORM SYMBOLIC ELIMINATION ON A COLUMN.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   ACOL2, ACOLQ2, BLKS2, BMAX, C2, COL, COL2, COLS2, DIAG, FILL,
     +   FMAX, LEVEL, R2, ROW, ROW2, ROWS2, TEXT
      LOGICAL ACTIVE, ERROR, FOUND

      PARAMETER (ID = ' BLMA4:  ')

      PARAMETER (QFRST = 1, QLAST = 2, QFRSTD = 3, QLASTD = 4)
C     PARAMETER (QIND = 1, QLOC = 2)
      PARAMETER (QIND = 1, QLEV = 2)

      DIMENSION
     +   ACOL2(2, BMAX), ACOLQ2(4, COLS2), ACTIVE(ROWS2),
     +   COL2(2, COLS2), DIAG(2, ROWS2), LEVEL(ROWS2), ROW2(2, ROWS2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. BMAX .AND. 0 .LT. COLS2 .AND. 0 .LT. ROWS2)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) MAKE THE ACOL POINTER STRUCTURE AND PART OF THE ACOLQ POINTER
C     STRUCTURE FOR THE NEW COLUMN.
C
C///////////////////////////////////////////////////////////////////////

C///  INITIALIZE THE POINTERS.

      ACOLQ2(QFRST, COL) = BLKS2 + 1
      ACOLQ2(QLAST, COL) = BLKS2

C///  TOP OF THE LOOP OVER THE ROWS.

      DO 2040 R2 = 1, ROWS2
         IF (LEVEL(R2) .LE. FILL) THEN

C///  PACK ACOL2.

      BLKS2 = BLKS2 + 1
      ERROR = .NOT. (BLKS2 .LE. BMAX)
      IF (ERROR) GO TO 9201

      ACOL2(QIND, BLKS2) = R2
      ACOL2(QLEV, BLKS2) = LEVEL(R2)
      ACOLQ2(QLAST, COL) = BLKS2

C///  SYMBOLIC ELIMINATION.

      IF (ACTIVE(R2) .AND. LEVEL(R2) .LT. FILL) THEN
         DO 2030 C2 = DIAG(QFRST, R2), DIAG(QLAST, R2)
            DO 2020 QD = ACOLQ2(QFRSTD, C2), ACOLQ2(QLASTD, C2)
               IF (ACOL2(QIND, QD) .EQ. R2) THEN
                  DO 2010 QR = QD + 1, ACOLQ2(QLAST, C2)
                     ROW = ACOL2(QIND, QR)
                     LEVEL(ROW) = MIN (LEVEL(ROW), 1 + MAX (
     +                  ACOL2(QLEV, QD), ACOL2(QLEV, QR),
     +                  LEVEL(R2)))
                     FMAX = MAX (FMAX, LEVEL(ROW))
2010              CONTINUE
               END IF
2020        CONTINUE
2030     CONTINUE
      END IF

C///  RESTORE THE LEVEL ARRAY FOR THE COLUMN.

      LEVEL(R2) = FILL + 1

C///  BOTTOM OF THE LOOP OVER THE ROWS.

         END IF
2040  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (3) COMPLETE THE ACOLQ POINTER STRUCTURE FOR THE NEW COLUMN.
C
C///////////////////////////////////////////////////////////////////////

C///  BUILD POINTERS TO THE DIAGONAL BLOCKS.

      ERROR = .NOT. (ACOLQ2(QFRST, COL) .LE. ACOLQ2(QLAST, COL))
      IF (ERROR) GO TO 9301

      FOUND = .FALSE.
      DO 3010 QR = ACOLQ2(QFRST, COL), ACOLQ2(QLAST, COL)
         R2 = ACOL2(QIND, QR)
         IF (MAX (COL2(QFRST, COL), ROW2(QFRST, R2)) .LE.
     +       MIN (COL2(QLAST, COL), ROW2(QLAST, R2))) THEN
            ACOLQ2(QLASTD, COL) = QR
            IF (.NOT. FOUND) THEN
               ACOLQ2(QFRSTD, COL) = QR
               FOUND = .TRUE.
            END IF
         END IF
3010  CONTINUE

      ERROR = .NOT. FOUND
      IF (ERROR) GO TO 9302

C///////////////////////////////////////////////////////////////////////
C
C     (4) UPDATE POINTERS TO THE FIRST AND LAST COLUMNS HAVING DIAGONAL
C     BLOCKS IN A GIVEN ROW.
C
C///////////////////////////////////////////////////////////////////////

      DO 4010 QD = ACOLQ2(QFRSTD, COL), ACOLQ2(QLASTD, COL)
         R2 = ACOL2(QIND, QD)
         IF (ACTIVE(R2)) THEN
            DIAG(QFRST, R2) = MIN (DIAG(QFRST, R2), COL)
            DIAG(QLAST, R2) = MAX (DIAG(QLAST, R2), COL)
         ELSE
            ACTIVE(R2) = .TRUE.
            DIAG(QFRST, R2) = COL
            DIAG(QLAST, R2) = COL
         END IF
4010  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, BMAX, COLS2, ROWS2
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID,
     +   2 * ((COLS2 - COL + 1) * ROWS2 - R2 + 1)
      GO TO 99999

9301  IF (0 .LT. TEXT) WRITE (TEXT, 99301) ID, COLS2
      GO TO 99999

9302  IF (0 .LT. TEXT) WRITE (TEXT, 99302) ID, COLS2
      GO TO 99999

99101 FORMAT
     +  (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +  //10X, I10, '  BMAX'
     +   /10X, I10, '  COLS2'
     +   /10X, I10, '  ROWS2')

99201 FORMAT
     +  (/1X, A9, 'ERROR.  THE INTEGER WORKSPACE IS TOO SMALL'
     +   /10X, 'TO COMPLETE THE SYMBOLIC FACTORIZATION.'
     +  //10X, I10, '  SUFFICIENT ADDITIONAL SPACE'
     +  //10X, 'THIS MAY BE A SEVERE OVERESTIMATE.')

99301 FORMAT
     +  (/1X, A9, 'ERROR.  A COLUMN HAS NO BLOCKS.'
     +  //10X, I10, '  COLUMN')

99302 FORMAT
     +  (/1X, A9, 'ERROR.  A COLUMN HAS NO DIAGONAL BLOCKS.'
     +  //10X, I10, '  COLUMN')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLMA5
     +  (ERROR, TEXT,
     +   ACOL2, ACOLQ2, AROW2, AROWQ2, BLKS2, COL2, COLS2, ROW2, ROWS2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLMA5
C
C>    CALLED BY BLMAK
C
C>    MAKE ROW INDICES FROM COLUMN INDICES.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER HEADER*80, ID*9
      INTEGER
     +   ACOL2, ACOLQ2, AROW2, AROWQ2, BLKS2, BLOCK, BLOCK2, C2, COL2,
     +   COLS2, COUNT, J, K, LAST, R2, ROW2, ROWS2, TEXT
      LOGICAL ERROR, FOUND

      PARAMETER (ID = ' BLMA5:  ')

      PARAMETER (QFRST = 1, QLAST = 2, QFRSTD = 3, QLASTD = 4)
      PARAMETER (QIND = 1, QLOC = 2)

      DIMENSION
     +   ACOL2(2, BLKS2), ACOLQ2(4, COLS2), AROW2(2, BLKS2),
     +   AROWQ2(4, ROWS2), COL2(2, COLS2), HEADER(2, 2), ROW2(2, ROWS2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. BLKS2 .AND. 0 .LT. COLS2 .AND. 0 .LT. ROWS2)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) MAKE THE AROW POINTER STRUCTURE AND PART OF THE AROWQ POINTER
C     STRUCTURE FOR THE NEW COLUMN.
C
C///////////////////////////////////////////////////////////////////////

C///  COUNT THE BLOCKS IN EACH ROW.

      DO 2010 R2 = 1, ROWS2
         AROWQ2(QFRST, R2) = 0
2010  CONTINUE

      DO 2030 C2 = 1, COLS2
         DO 2020 BLOCK2 = ACOLQ2(QFRST, C2), ACOLQ2(QLAST, C2)
            R2 = ACOL2(QIND, BLOCK2)
            AROWQ2(QFRST, R2) = AROWQ2(QFRST, R2) + 1
2020     CONTINUE
2030  CONTINUE

C///  BUILD THE FIRST POINTERS.

      LAST = 0
      DO 2040 R2 = 1, ROWS2
         COUNT = AROWQ2(QFRST, R2)
         AROWQ2(QFRST, R2) = LAST + 1
         AROWQ2(QLAST, R2) = LAST
         LAST = LAST + COUNT
2040  CONTINUE

      ERROR = .NOT. (BLKS2 .EQ. LAST)
      IF (ERROR) GO TO 9201

C///  BUILD THE BLOCK POINTERS AND THE LAST POINTERS.

      DO 2060 C2 = 1, COLS2
         DO 2050 BLOCK2 = ACOLQ2(QFRST, C2), ACOLQ2(QLAST, C2)
            R2 = ACOL2(QIND, BLOCK2)
            AROWQ2(QLAST, R2) = AROWQ2(QLAST, R2) + 1
            BLOCK = AROWQ2(QLAST, R2)
            AROW2(QIND, BLOCK) = C2
            AROW2(QLOC, BLOCK) = ACOL2(QLOC, BLOCK2)
2050     CONTINUE
2060  CONTINUE

C///  BUILD THE POINTERS TO THE DIAGONAL BLOCKS.

      DO 2080 R2 = 1, ROWS2
         ERROR = .NOT. (AROWQ2(QFRST, R2) .LE. AROWQ2(QLAST, R2))
         IF (ERROR) GO TO 9202

         FOUND = .FALSE.
         DO 2070 BLOCK2 = AROWQ2(QFRST, R2), AROWQ2(QLAST, R2)
            C2 = AROW2(QIND, BLOCK2)
            IF (MAX (COL2(QFRST, C2), ROW2(QFRST, R2)) .LE.
     +          MIN (COL2(QLAST, C2), ROW2(QLAST, R2))) THEN
               AROWQ2(QLASTD, R2) = BLOCK2
               IF (.NOT. FOUND) THEN
                  AROWQ2(QFRSTD, R2) = BLOCK2
                  FOUND = .TRUE.
               END IF
            END IF
2070     CONTINUE

         ERROR = .NOT. FOUND
         IF (ERROR) GO TO 9203
2080  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, BLKS2, COLS2, ROWS2
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID, BLKS2, LAST
      GO TO 99999

9202  IF (0 .LT. TEXT) WRITE (TEXT, 99202) ID, R2
      GO TO 99999

9203  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99203) ID,
     +      ROWS2, R2, AROWQ2(QLAST, R2) - AROWQ2(QFRST, R2) + 1, HEADER
         DO 9291 BLOCK2 = AROWQ2(QFRST, R2), AROWQ2(QLAST, R2)
            C2 = AROW2(QIND, BLOCK2)
            WRITE (TEXT, 99291)
     +         BLOCK2 - AROWQ2(QFRST, R2) + 1, C2, COL2(QFRST, C2),
     +         COL2(QLAST, C2), ROW2(QFRST, R2), ROW2(QLAST, R2)
9291     CONTINUE
      END IF
      GO TO 99999

99101 FORMAT
     +  (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +  //10X, I10, '  BLKS2'
     +   /10X, I10, '  COLS2'
     +   /10X, I10, '  ROWS2')

99201 FORMAT
     +  (/1X, A9, 'ERROR.  THE NUMBER OF BLOCKS FOUND IN THE COLUMN'
     +   /10X, 'POINTERS IS NOT THE NUMBER EXPECTED.'
     +  //10X, I10, '  FOUND'
     +   /10X, I10, '  EXPECTED')

99202 FORMAT
     +  (/1X, A9, 'ERROR.  A ROW HAS NO BLOCKS.'
     +  //10X, I10, '  ROW')

      DATA ((HEADER(J, K), K = 1, 2), J = 1, 2)
C         123456789_123456789_123456789_1
C         12345  123456789_  123456789_
     + / '           COLUMN       FIRST  ',
     +   'BLOCK       INDEX      COLUMN  ',
C         123456789_123456789_123456789_1234
C         123456789_  123456789_  123456789_
     +   '      LAST       FIRST        LAST',
     +   '    COLUMN         ROW         ROW' /

99203 FORMAT
     +  (/1X, A9, 'ERROR.  A ROW HAS NO DIAGONAL BLOCKS.'
     +  //10X, I10, '  TOTAL ROWS'
     +   /10X, I10, '  INDEX OF THIS ROW'
     +   /10X, I10, '  BLOCKS IN THIS ROW'
     +  //10X, A31, A34
     +   /10X, A31, A34)

99291 FORMAT
     +  (10X, I5, 5(2X, I10))

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLMAK
     +  (ERROR, TEXT, PRINT,
     +   ILAST, IMAX, ISIZE, IWORK,
     +   LLAST, LMAX, LSIZE, LWORK,
     +   RLAST, RMAX, RSIZE, RWORK,
     +   QIDAT, QRDAT,
     +   IDSIZ, RDSIZ,
     +   BLOCKS, BORDER, CPOINT, FILL, RINDEX, SIZE, TOTAL)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLMAK
C
C>    CALLED BY THE  U S E R
C
C>    MAKE THE DATA SPACE.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   ASIZE2, BLKS1, BLKS2, BLOCKS, BMAX, BORDER, CPOINT, COLS1,
     +   COLS2, FILL, GLUES, IDSIZ, ILAST, IMARK, IMAX, ISAVE, ISIZE,
     +   IWORK, J, LLAST, LMARK, LMAX, LSAVE, LSIZE, N1, N2, PRINT,
     +   RDSIZ, RINDEX, RLAST, RMARK, RMAX, ROWS1, ROWS2, RSAVE, RSIZE,
     +   SIZE, TEXT, TOTAL
      LOGICAL ERROR, LWORK
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   RWORK

      PARAMETER (ID = ' BLMAK:  ')

      DIMENSION
     +   CPOINT(BLOCKS), IWORK(ISIZE), LWORK(LSIZE), RINDEX(TOTAL),
     +   RWORK(RSIZE), SIZE(BLOCKS)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  SAVE THE LEVELS OF THE WORK SPACES.

      ISAVE = ILAST
      LSAVE = LLAST
      RSAVE = RLAST

C///  PRINT.

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) WRITE (TEXT, 10001) ID

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. BLOCKS .AND. 0 .LT. ISIZE .AND. 0 .LT. LSIZE
     +   .AND. 0 .LT. RSIZE .AND. 0 .LT. TOTAL)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) ALLOCATE SOME DATA AND WORK SPACE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHOOSE THE DIMENSIONS.

      BLKS1 = TOTAL
      COLS1 = BLOCKS
      COLS2 = BLOCKS
      ROWS1 = BLOCKS
      ROWS2 = BLOCKS

      IF (0 .EQ. BORDER) THEN
         COLS2 = BLOCKS
         GLUES = 0
      ELSE IF (0 .LT. BORDER) THEN
         COLS2 = 2 * BLOCKS
         GLUES = BLOCKS - 1
      ELSE
         ERROR = .TRUE.
         GO TO 9201
      END IF

C///  RESERVE SOME DATA SPACE.

C     SUBROUTINE RESERV
C    +  (ERROR, TEXT,
C    +   FATAL, QLAST, QMAX, QSIZE,
C    +   NAME, NUMBER, Q)

C     INTEGER IDATA'S 11 SCALARS

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QISCLR', 11, QISCLR)
      IF (ERROR) GO TO 9202

C     INTEGER ACOL1(2, BLKS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QAC1', 2 * BLKS1, QAC1)
      IF (ERROR) GO TO 9202

C     INTEGER ACOLQ1(2, COLS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QACQ1', 2 * COLS1, QACQ1)
      IF (ERROR) GO TO 9202

C     INTEGER ACOLQ2(4, COLS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QACQ2', 4 * COLS2, QACQ2)
      IF (ERROR) GO TO 9202

C     INTEGER BMAP(COLS1 + 1)

      IF (0 .LT. BORDER) THEN
         CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +      'QBMAP', COLS1 + 1, QBMAP)
         IF (ERROR) GO TO 9202
      END IF

C     INTEGER CMAP(COLS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QCMAP', COLS1, QCMAP)
      IF (ERROR) GO TO 9202

C     INTEGER COL1(2, COLS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QCOL1', 2 * COLS1, QCOL1)
      IF (ERROR) GO TO 9202

C     INTEGER COL2(2, COLS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QCOL2', 2 * COLS2, QCOL2)
      IF (ERROR) GO TO 9202

C     INTEGER GMAP(3, 0 : GLUES)

      IF (0 .LT. BORDER) THEN
         CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +      'QGMAP', 3 * (1 + GLUES), QGMAP)
         IF (ERROR) GO TO 9202
      END IF

C     INTEGER RMAP(ROWS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QRMAP', ROWS1, QRMAP)
      IF (ERROR) GO TO 9202

C     INTEGER ROW1(2, ROWS1)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QROW1', 2 * ROWS1, QROW1)
      IF (ERROR) GO TO 9202

C     INTEGER ROW2(2, ROWS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QROW2', 2 * ROWS2, QROW2)
      IF (ERROR) GO TO 9202

C///  MARK THE LEVELS OF THE WORK SPACES.

      IMARK = ILAST
      LMARK = LLAST
      RMARK = RLAST

C///  RESERVE SOME WORK SPACE.

C     SUBROUTINE RESERV
C    +  (ERROR, TEXT,
C    +   FATAL, QLAST, QMAX, QSIZE,
C    +   NAME, NUMBER, Q)

C     LOGICAL ACTIVE(ROWS2)

      CALL RESERV (ERROR, TEXT, .FALSE., LLAST, LMAX, LSIZE,
     +   'QACTIV', ROWS2, QACTIV)
      IF (ERROR) GO TO 9202

C     INTEGER DIAG(2, ROWS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QDIAG', 2 * ROWS2, QDIAG)
      IF (ERROR) GO TO 9202

C     INTEGER ITEMP(ROWS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QITEMP', ROWS2, QITEMP)
      IF (ERROR) GO TO 9202

C     LOGICAL OCCUPY(ROWS1)

      CALL RESERV (ERROR, TEXT, .FALSE., LLAST, LMAX, LSIZE,
     +   'QOCCUP', ROWS2, QOCCUP)
      IF (ERROR) GO TO 9202

C///  RESERVE SOME TEMPORARY DATA SPACE.

C     SUBROUTINE RESERV
C    +  (ERROR, TEXT,
C    +   FATAL, QLAST, QMAX, QSIZE,
C    +   NAME, NUMBER, Q)

C     INTEGER ACOL2(2, BMAX)

      BMAX = BLKS1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QAC2', 2 * BMAX, QAC2)
      IF (ERROR) GO TO 9202

C///  CHECK FOR ADEQUATE WORK SPACE.

      ERROR = .NOT. (ILAST .LE. ISIZE .AND. LLAST .LE. LSIZE .AND.
     +   RLAST .LE. RSIZE)
      IF (ERROR) GO TO 9203

C///  EXTEND THE TEMPORARY DATA SPACE TO THE ENTIRE INTEGER WORK SPACE.

      BMAX = BMAX + (ISIZE - ILAST) / 2

C///////////////////////////////////////////////////////////////////////
C
C     (3) CALL BLMA2 OR BLMA3.
C
C///////////////////////////////////////////////////////////////////////

      IF (0 .EQ. BORDER) THEN

C     SUBROUTINE BLMA2
C    +  (ERROR, TEXT, PRINT,
C    +   ACOL1, ACOLQ1, ACOL2, ACOLQ2, ACTIVE, ASIZE2, BLKS1, BLKS2,
C    +   BLOCKS, BMAX, CMAP, COL1, COL2, COLS1, COLS2, CPOINT, DIAG,
C    +   FILL, ITEMP, OCCUPY, N1, N2, RINDEX, RMAP, ROW1, ROW2, ROWS1,
C    +   ROWS2, SIZE, TOTAL)

      CALL BLMA2
     +  (ERROR, TEXT, PRINT - 1,
     +   IWORK(QAC1), IWORK(QACQ1), IWORK(QAC2), IWORK(QACQ2),
     +   LWORK(QACTIV), ASIZE2, BLKS1, BLKS2, BLOCKS, BMAX,
     +   IWORK(QCMAP), IWORK(QCOL1), IWORK(QCOL2), COLS1, COLS2, CPOINT,
     +   IWORK(QDIAG), FILL, IWORK(QITEMP), LWORK(QOCCUP), N1, N2,
     +   RINDEX, IWORK(QRMAP), IWORK(QROW1), IWORK(QROW2), ROWS1, ROWS2,
     +   SIZE, TOTAL)
      IF (ERROR) GO TO 9301

      ELSE

C     SUBROUTINE BLMA3
C    +  (ERROR, TEXT, PRINT,
C    +   ACOL1, ACOLQ1, ACOL2, ACOLQ2, ACTIVE, ASIZE2, BLKS1, BLKS2,
C    +   BLOCKS, BMAP, BMAX, BORDER, CMAP, COL1, COL2, COLS1, COLS2,
C    +   CPOINT, DIAG, FILL, GLUES, GMAP, ITEMP, N1, N2, OCCUPY, RINDEX,
C    +   RMAP, ROW1, ROW2, ROWS1, ROWS2, SIZE, TOTAL)

      CALL BLMA3
     +  (ERROR, TEXT, PRINT - 1,
     +   IWORK(QAC1), IWORK(QACQ1), IWORK(QAC2), IWORK(QACQ2),
     +   LWORK(QACTIV), ASIZE2, BLKS1, BLKS2, BLOCKS, IWORK(QBMAP),
     +   BMAX, BORDER, IWORK(QCMAP), IWORK(QCOL1), IWORK(QCOL2), COLS1,
     +   COLS2, CPOINT, IWORK(QDIAG), FILL, GLUES, IWORK(QGMAP),
     +   IWORK(QITEMP), N1, N2, LWORK(QOCCUP), RINDEX, IWORK(QRMAP),
     +   IWORK(QROW1), IWORK(QROW2), ROWS1, ROWS2, SIZE, TOTAL)
      IF (ERROR) GO TO 9302

      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (4) ALLOCATE MORE DATA SPACE.
C
C///////////////////////////////////////////////////////////////////////

C///  RELEASE THE TEMPORARY SPACE.

      ILAST = IMARK
      LLAST = LMARK
      RLAST = RMARK

C///  MOVE THE ACOL2 SPACE.

      Q = QAC2

C     SUBROUTINE RESERV
C    +  (ERROR, TEXT,
C    +   FATAL, QLAST, QMAX, QSIZE,
C    +   NAME, NUMBER, Q)

C     INTEGER ACOL2(2, BLKS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QAC2', 2 * BLKS2, QAC2)
      IF (ERROR) GO TO 9202

      DO 4010 J = 0, 2 * BLKS2 - 1
         IWORK(QAC2 + J) = IWORK(Q + J)
4010  CONTINUE

C///  RESERVE SOME DATA SPACE.

C     SUBROUTINE RESERV
C    +  (ERROR, TEXT,
C    +   FATAL, QLAST, QMAX, QSIZE,
C    +   NAME, NUMBER, Q)

C     REAL A2(ASIZE2)

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QA2', ASIZE2, QA2)
      IF (ERROR) GO TO 9202

C     INTEGER AROW2(2, BLKS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QAR2', 2 * BLKS2, QAR2)
      IF (ERROR) GO TO 9202

C     INTEGER AROWQ2(4, ROWS2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QARQ2', 4 * ROWS2, QARQ2)
      IF (ERROR) GO TO 9202

C     INTEGER RPIVOT(N2)

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QRPIV', N2, QRPIV)
      IF (ERROR) GO TO 9202

C     REAL RSCALE(N2)

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QRSCAL', N2, QRSCAL)
      IF (ERROR) GO TO 9202

C///  CHECK FOR ADEQUATE WORK SPACE.

      ERROR = .NOT. (ILAST .LE. ISIZE .AND. LLAST .LE. LSIZE .AND.
     +   RLAST .LE. RSIZE)
      IF (ERROR) GO TO 9401

C///////////////////////////////////////////////////////////////////////
C
C     (5) MAKE THE ROW POINTER STRUCTURE.
C
C///////////////////////////////////////////////////////////////////////

C     SUBROUTINE BLMA5
C    +  (ERROR, TEXT,
C    +   ACOL2, ACOLQ2, AROW2, AROWQ2, BLKS2, COL2, COLS2, ROW2, ROWS2)

      CALL BLMA5
     +  (ERROR, TEXT,
     +   IWORK(QAC2), IWORK(QACQ2), IWORK(QAR2), IWORK(QARQ2), BLKS2,
     +   IWORK(QCOL2), COLS2, IWORK(QROW2), ROWS2)
      IF (ERROR) GO TO 9501

C///////////////////////////////////////////////////////////////////////
C
C     (6) COMPLETE THE DATA STRUCTURE.
C
C///////////////////////////////////////////////////////////////////////

C///  PACK THE SCALARS.

      IWORK(QISCLR + 00) = ASIZE2
      IWORK(QISCLR + 01) = BLKS1
      IWORK(QISCLR + 02) = BLKS2
      IWORK(QISCLR + 03) = BORDER
      IWORK(QISCLR + 04) = COLS1
      IWORK(QISCLR + 05) = COLS2
      IWORK(QISCLR + 06) = GLUES
      IWORK(QISCLR + 07) = N1
      IWORK(QISCLR + 08) = N2
      IWORK(QISCLR + 09) = ROWS1
      IWORK(QISCLR + 10) = ROWS2

C///  CHOOSE THE POINTERS AND SIZES.

      QIDAT = ISAVE + 1
      IDSIZ = ILAST - ISAVE

      QRDAT = RSAVE + 1
      RDSIZ = RLAST - RSAVE

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +   (/1X, A, 'MAKE DATA STRUCTURES FOR A STRETCHED BORDERED SPARSE'
     +   /10X, 'BLOCKED MATRIX.')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   BLOCKS, ISIZE, LSIZE, RSIZE, TOTAL
      GO TO 99998

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID
      GO TO 99998

9202  IF (0 .LT. TEXT) WRITE (TEXT, 99202) ID
      GO TO 99998

9203  IF (0 .LT. TEXT) WRITE (TEXT, 99203) ID,
     +   0, ISIZE, LSIZE, RSIZE, 0, ILAST, LLAST, RLAST
      GO TO 99999

9301  IF (0 .LT. TEXT) WRITE (TEXT, 99301) ID
      GO TO 99998

9302  IF (0 .LT. TEXT) WRITE (TEXT, 99302) ID
      GO TO 99998

9401  IF (0 .LT. TEXT) WRITE (TEXT, 99401) ID,
     +   0, ISIZE, LSIZE, RSIZE, 0, ILAST, LLAST, RLAST
      GO TO 99998

9501  IF (0 .LT. TEXT) WRITE (TEXT, 99501) ID
      GO TO 99998

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  BLOCKS'
     +    /10X, I10, '  ISIZE'
     +    /10X, I10, '  LSIZE'
     +    /10X, I10, '  RSIZE'
     +    /10X, I10, '  TOTAL')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  THERE MUST BE ZERO OR MORE BORDERS.'
     +   //10X, I10, '  BORDER')

99202 FORMAT
     +   (/1X, A9, 'ERROR.  RESERV FAILS.')

99203 FORMAT
     +   (/1X, A9, 'ERROR.  SOME WORKSPACES ARE TOO SMALL.'
C               123456789_  123456789_  123456789_  123456789_
     +  //25X, ' CHARACTER     INTEGER     LOGICAL        REAL'
C               123456789_12345
     +  //10X, ' PRESENT SIZE', 4I12
     +   /10X, 'REQUIRED SIZE', 4I12
     +  //10X, 'THIS SUBROUTINE WILL NEED MORE SPACE LATER.')

99301 FORMAT
     +   (/1X, A9, 'ERROR.  BLMA2 FAILS.')

99302 FORMAT
     +   (/1X, A9, 'ERROR.  BLMA3 FAILS.')

99401 FORMAT
     +   (/1X, A9, 'ERROR.  SOME WORKSPACES ARE TOO SMALL.'
C               123456789_  123456789_  123456789_  123456789_
     +  //25X, ' CHARACTER     INTEGER     LOGICAL        REAL'
C               123456789_12345
     +  //10X, ' PRESENT SIZE', 4I12
     +   /10X, 'REQUIRED SIZE', 4I12
     +  //10X, 'THIS SUBROUTINE WILL NEED NO MORE SPACE.')

99501 FORMAT
     +   (/1X, A9, 'ERROR.  BLMA5 FAILS.')

C///  EXIT.

99998 CONTINUE
      ILAST = ISAVE
      LLAST = LSAVE
      RLAST = RSAVE
99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLPU2
     +  (ERROR, TEXT,
     +   A2, ACOL2, ACOLQ2, ASIZE2, BCMAX, BCOL, BLKS2, BMAP, BORDER,
     +   BRMAX, BROW, CMAP, COL1, COL2, COLS1, COLS2, GLUES, GMAP, N1,
     +   N2, RMAP, ROW1, ROW2, ROWS1, ROWS2, RSCALE)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLPU2
C
C>    CALLED BY BLPUT
C
C>    PUT A BORDER ON THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   AC1, AC2, ACOL2, ACOLQ2, QA, AR1, AR2, ASIZE2, BCMAX, BLKS2,
     +   BMAP, BORDER, BR1, BR2, BRMAX, C1, C2, CMAP, COL, COL1, COL2,
     +   COLS1, COLS2, GLUE, GLUES, GMAP, INDEX, N1, N2, R1, R2, RMAP,
     +   ROW, ROW1, ROW2, ROWS1, ROWS2, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A2, BCOL, BROW, RSCALE, VALUE

      PARAMETER (ID = ' BLPU2:  ')

      PARAMETER (QFRST = 1, QLAST = 2, QFRSTD = 3, QLASTD = 4)
      PARAMETER (QIND = 1, QLOC = 2)
      PARAMETER (QCOL = 1, QNEG = 2, QPOS = 3)

      DIMENSION
     +   A2(ASIZE2), ACOL2(2, BLKS2), ACOLQ2(4, COLS2),
     +   BCOL(BCMAX, BORDER), BMAP(COLS1 + 1), BROW(BRMAX, BORDER),
     +   CMAP(COLS1), COL1(2, COLS1), COL2(2, COLS2),
     +   GMAP(3, 0 : GLUES), RMAP(ROWS1), ROW1(2, ROWS1),
     +   ROW2(2, ROWS2), RSCALE(N2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. ASIZE2 .AND. 0 .LT. BCMAX .AND.
     +   0 .LT. BLKS2 .AND. 0 .LT. BORDER .AND. 0 .LT. BRMAX .AND.
     +   0 .LT. COLS1 .AND. 0 .LT. COLS2 .AND. 0 .LT. GLUES .AND.
     +   0 .LT. N2 .AND. 0 .LT. ROWS1)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) PLACE THE BORDER ROW.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP OVER THE COLUMNS IN THE BORDER.

      DO 2010 C1 = 1, COLS1 + 1
         Q = BMAP(C1)
         QA = ACOL2(QLOC, Q)
         R2 = ACOL2(QIND, Q)
         AR1 = ROW2(QFRST, R2)
         AR2 = ROW2(QLAST, R2)

C///  BORDER COLUMN.

      IF (C1 .EQ. COLS1 + 1) THEN
         BR1 = N1 + 1
         BR2 = N1 + BORDER
         C2 = COLS2

C///  COLUMNS OF A1.

      ELSE
         BR1 = COL1(QFRST, C1)
         BR2 = COL1(QLAST, C1)
         C2 = CMAP(C1)
      END IF

C///  COPY THE BLOCK.

      AC1 = COL2(QFRST, C2)
      AC2 = COL2(QLAST, C2)

C     SUBROUTINE BLPU3
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, AR1, AR2, BRMAX, BROW, BORDER, BR1, BR2)

      CALL BLPU3
     +  (ERROR, TEXT,
     +   A2(QA), AC1, AC2, AR1, AR2, BRMAX, BROW, BORDER, BR1, BR2)
      IF (ERROR) GO TO 9201

C///  BOTTOM OF THE LOOP OVER THE COLUMNS IN THE BORDER.

2010  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (3) FORM AND PLACE THE BORDER ROW SCALE FACTORS.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP OVER THE COLUMNS IN THE BORDER.

      DO 3040 C1 = 1, COLS1 + 1
         Q = BMAP(C1)
         R2 = ACOL2(QIND, Q)
         AR1 = ROW2(QFRST, R2)
         AR2 = ROW2(QLAST, R2)

C///  BORDER COLUMN.

      IF (C1 .EQ. COLS1 + 1) THEN
         BR1 = N1 + 1
         BR2 = N1 + BORDER
         C2 = COLS2

C///  COLUMNS OF A1.

      ELSE
         BR1 = COL1(QFRST, C1)
         BR2 = COL1(QLAST, C1)
         C2 = CMAP(C1)
      END IF

C///  SUM THE ROWS.

      IF (C1 .EQ. 1) THEN
         DO 3020 ROW = 1, BORDER
            INDEX = AR2 - BORDER + ROW
            RSCALE(INDEX) = 0
            DO 3010 COL = 1, N1 + BORDER
               RSCALE(INDEX) = RSCALE(INDEX) + ABS (BROW(COL, ROW))
3010        CONTINUE
3020     CONTINUE

         INDEX = AR2 - BORDER

C///  COPY THE SUMS.

      ELSE
         DO 3030 ROW = 1, BORDER
            RSCALE(AR2 - BORDER + ROW) = RSCALE(INDEX + BORDER)
3030     CONTINUE
      END IF

C///  BOTTOM OF THE LOOP OVER THE COLUMNS IN THE BORDER.

3040  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (4) PLACE THE COLUMN.
C
C///////////////////////////////////////////////////////////////////////

      ERROR = .NOT. (ROWS1 .EQ. ROWS2 .AND.
     +   ACOLQ2(QFRST, COLS2) + ROWS2 .EQ. ACOLQ2(QLAST, COLS2) + 1)
      IF (ERROR) GO TO 9401

C///  TOP OF THE LOOP OVER THE ROWS OF A1.

      AC1 = COL2(QFRST, COLS2)
      AC2 = COL2(QLAST, COLS2)
      DO 4010 R1 = 1, ROWS1
         R2 = RMAP(R1)

C///  GET THE POINTER TO THE BLOCK OF A2.

      Q = ACOLQ2(QFRST, COLS2) + R2 - 1
      ERROR = .NOT. (R2 .EQ. ACOL2(QIND, Q))
      IF (ERROR) GO TO 9402

C///  COPY THE BLOCK.

      QA = ACOL2(QLOC, Q)

      AR1 = ROW2(QFRST, R2)
      AR2 = ROW2(QLAST, R2)

      BR1 = ROW1(QFRST, R1)
      BR2 = ROW1(QLAST, R1)

C     SUBROUTINE BLPU4
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, AR1, AR2, BCOL, BORDER, BR1, BR2, NMAX, RSCALE)

      CALL BLPU4
     +  (ERROR, TEXT,
     +   A2(QA), AC1, AC2, AR1, AR2, BCOL, BORDER, BR1, BR2, BCMAX,
     +   RSCALE(AR1))
      IF (ERROR) GO TO 9403

C///  BOTTOM OF THE LOOP OVER THE ROWS OF A1.

4010  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (5) PLACE THE GLUE.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP OVER THE GLUE COLUMNS.

      DO 5010 GLUE = 1, GLUES
         C2 = GMAP(QCOL, GLUE)
         AC1 = COL2(QFRST, C2)
         AC2 = COL2(QLAST, C2)

C///  PLACE THE POSITIVE GLUE.

      Q = GMAP(QPOS, GLUE)
      QA = ACOL2(QLOC, Q)
      R2 = ACOL2(QIND, Q)
      AR1 = ROW2(QFRST, R2)
      AR2 = ROW2(QLAST, R2)

      VALUE = 1.0

C     SUBROUTINE BLPU5
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, AR1, AR2, BORDER, VALUE)

      CALL BLPU5
     +  (ERROR, TEXT,
     +   A2(QA), AC1, AC2, AR1, AR2, BORDER, VALUE)
      IF (ERROR) GO TO 9501

C///  PLACE THE NEGATIVE GLUE.

      Q = GMAP(QNEG, GLUE)
      QA = ACOL2(QLOC, Q)
      R2 = ACOL2(QIND, Q)
      AR1 = ROW2(QFRST, R2)
      AR2 = ROW2(QLAST, R2)

      VALUE = - 1.0

C     SUBROUTINE BLPU5
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, AR1, AR2, BORDER, VALUE)

      CALL BLPU5
     +  (ERROR, TEXT,
     +   A2(QA), AC1, AC2, AR1, AR2, BORDER, VALUE)
      IF (ERROR) GO TO 9501

C///  BOTTOM OF THE LOOP OVER THE GLUE COLUMNS.

5010  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   ASIZE2, BCMAX, BLKS2, BORDER, BRMAX, COLS1, COLS2, GLUES, N2,
     +   ROWS1
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID
      GO TO 99999

9401  IF (0 .LT. TEXT) WRITE (TEXT, 99401) ID,
     +   ROWS1, ROWS2, ACOLQ2(QFRST, COLS2) - ACOLQ2(QLAST, COLS2) + 1
      GO TO 99999

9402  IF (0 .LT. TEXT) WRITE (TEXT, 99402) ID, R2, ACOL2(QIND, Q)
      GO TO 99999

9403  IF (0 .LT. TEXT) WRITE (TEXT, 99403) ID
      GO TO 99999

9501  IF (0 .LT. TEXT) WRITE (TEXT, 99501) ID
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  ASIZE2',
     +    /10X, I10, '  BCMAX',
     +    /10X, I10, '  BLKS2',
     +    /10X, I10, '  BORDER',
     +    /10X, I10, '  BRMAX',
     +    /10X, I10, '  COLS1',
     +    /10X, I10, '  COLS2',
     +    /10X, I10, '  GLUES',
     +    /10X, I10, '  N2',
     +    /10X, I10, '  ROWS1')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  BLPU3 FAILS.')

99401 FORMAT
     +   (/1X, A9, 'ERROR.  A1 AND A2 MUST HAVE THE SAME ROW DIMENSION'
     +    /10X, 'AND THE FINAL COLUMN OF A2 MUST BE DENSE.'
     +   //10X, I10, '  ROWS1'
     +    /10X, I10, '  ROWS2'
     +    /10X, I10, '  ROWS IN FINAL COLUMN')

99402 FORMAT
     +   (/1X, A9, 'ERROR.  THE POINTERS FOR THE FINAL COLUMN OF A2 ARE'
     +    /10X, 'OUT OF ORDER.'
     +   //10X, I10, '  EXPECTED ROW'
     +    /10X, I10, '  ACTUAL ROW')

99403 FORMAT
     +   (/1X, A9, 'ERROR.  BLPU4 FAILS.')

99501 FORMAT
     +   (/1X, A9, 'ERROR.  BLPU5 FAILS.')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLPU3
     +  (ERROR, TEXT,
     +   A, AC1, AC2, AR1, AR2, BRMAX, BROW, BORDER, BR1, BR2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLPU3
C
C>    CALLED BY BLPU2
C
C>    COPY A ROW BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   AC1, AC2, AR1, AR2, BC, BORDER, BR, BR1, BR2, BRMAX, SC1, SC2,
     +   SR1, SR2, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, BROW

      PARAMETER (ID = ' BLPU3:  ')

      DIMENSION A(AR1 : AR2, AC1 : AC2), BROW(BRMAX, BORDER)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (AC1 .LE. AC2 .AND. AR1 .LE. AR2 .AND. 0 .LT. BORDER
     +   .AND. 0 .LT. BRMAX)
      IF (ERROR) GO TO 9101

C///  CHECK THE PLACEMENT.

      ERROR = .NOT. (1 .LE. BR1 .AND. BR1 .LE. BR2 .AND. BR2 .LE. BRMAX)
      IF (ERROR) GO TO 9102

      SC1 = AC1
      SC2 = AC1 + BR2 - BR1
      SR1 = AR2 - BORDER + 1
      SR2 = AR2
      ERROR = .NOT. (AC1 .EQ. SC1 .AND. SC2 .EQ. AC2 .AND. AR1 .LE. SR1
     +   .AND. SR2 .EQ. AR2)
      IF (ERROR) GO TO 9103

C///////////////////////////////////////////////////////////////////////
C
C     (2) PLACE THE BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      DO 2020 BC = 1, BORDER
         DO 2010 BR = BR1, BR2
            A(AR2 - BORDER + BC, AC1 - BR1 + BR) = BROW(BR, BC)
2010     CONTINUE
2020  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   AC1, AC2, AR1, AR2, BORDER, BRMAX
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, BR1, BR2, BRMAX
      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID,
     +   AC1, SC1, SC2, AC2, AR1, SR1, SR2, AR2

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2'
     +   //10X, I10, '  BORDER'
     +    /10X, I10, '  BRMAX')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  THE ARRAY OF BORDER ROWS DOES NOT CONTAIN'
     +    /10X, 'THE BLOCK TO COPY.'
     +   //10X, I10, '  BR1'
     +    /10X, I10, '  BR2'
     +    /10X, I10, '  BRMAX')

99103 FORMAT
     +   (/1X, A9, 'ERROR.  THE COPIED BLOCK IS OFF TARGET.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  SC1'
     +    /10X, I10, '  SC2'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  SR1'
     +    /10X, I10, '  SR2'
     +    /10X, I10, '  AR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLPU4
     +  (ERROR, TEXT,
     +   A, AC1, AC2, AR1, AR2, BCOL, BORDER, BR1, BR2, NMAX, RSCALE)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLPU4
C
C>    CALLED BY BLPU2
C
C>    COPY A COLUMN BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   AC, AC1, AC2, AR, AR1, AR2, BC, BORDER, BR, BR1, BR2, NMAX,
     +   SC1, SC2, SR1, SR2, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, BCOL, RSCALE

      PARAMETER (ID = ' BLPU4:  ')

      DIMENSION
     +   A(AR1 : AR2, AC1 : AC2), BCOL(NMAX, BORDER), RSCALE(AR1 : AR2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (AC1 .LE. AC2 .AND. AR1 .LE. AR2 .AND. 0 .LT. BORDER
     +   .AND. 0 .LT. NMAX)
      IF (ERROR) GO TO 9101

C///  CHECK THE PLACEMENT.

      ERROR = .NOT. (1 .LE. BR1 .AND. BR1 .LE. BR2 .AND. BR2 .LE. NMAX)
      IF (ERROR) GO TO 9102

      SC1 = AC1
      SC2 = AC1 - 1 + BORDER
      SR1 = AR1
      SR2 = AR1 - 1 + BR2 - BR1
      ERROR = .NOT. (AC1 .EQ. SC1 .AND. SC2 .EQ. AC2 .AND. AR1 .EQ. SR1
     +   .AND. SR2 .LE. AR2)
      IF (ERROR) GO TO 9103

C///////////////////////////////////////////////////////////////////////
C
C     (2) PLACE THE BLOCK AND UPDATE THE SCALE FACTORS.
C
C///////////////////////////////////////////////////////////////////////

      DO 2020 BC = 1, BORDER
         DO 2010 BR = BR1, BR2
            AC = AC1 - 1 + BC
            AR = AR1 - BR1 + BR
            A(AR, AC) = BCOL(BR, BC)
            RSCALE(AR) = RSCALE(AR) + ABS (A(AR, AC))
2010     CONTINUE
2020  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   AC1, AC2, AR1, AR2, BORDER, NMAX
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, BR1, BR2, NMAX
      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID,
     +   AC1, SC1, SC2, AC2, AR1, SR1, SR2, AR2

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2'
     +   //10X, I10, '  BORDER'
     +    /10X, I10, '  NMAX')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  THE ARRAY OF BORDER ROWS DOES NOT CONTAIN'
     +    /10X, 'THE BLOCK TO COPY.'
     +   //10X, I10, '  BR1'
     +    /10X, I10, '  BR2'
     +    /10X, I10, '  NMAX')

99103 FORMAT
     +   (/1X, A9, 'ERROR.  THE COPIED BLOCK IS OFF TARGET.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  SC1'
     +    /10X, I10, '  SC2'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  SR1'
     +    /10X, I10, '  SR2'
     +    /10X, I10, '  AR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLPU5
     +  (ERROR, TEXT,
     +   A, AC1, AC2, AR1, AR2, BORDER, VALUE)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLPU5
C
C>    CALLED BY BLPU2
C
C>    MAKE A GLUE BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER AC1, AC2, AR1, AR2, BC, BORDER, SC1, SC2, SR1, SR2, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, VALUE

      PARAMETER (ID = ' BLPU5:  ')

      DIMENSION A(AR1 : AR2, AC1 : AC2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (AC1 .LE. AC2 .AND. AR1 .LE. AR2 .AND.
     +   0 .LT. BORDER)
      IF (ERROR) GO TO 9101

C///  CHECK THE VALUE.

      ERROR = VALUE .EQ. 0.0
      IF (ERROR) GO TO 9102

C///  CHECK THE PLACEMENT.

      SC1 = AC1
      SC2 = AC1 + BORDER - 1
      SR1 = AR2 - BORDER + 1
      SR2 = AR2
      ERROR = .NOT. (AC1 .LE. SC1 .AND. SC2 .LE. AC2 .AND. AR1 .LE. SR1
     +   .AND. SR2 .EQ. AR2)
      IF (ERROR) GO TO 9103

C///////////////////////////////////////////////////////////////////////
C
C     (2) PLACE THE BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      DO 2010 BC = 1, BORDER
         A(AR2 - BORDER + BC, AC1 - 1 + BC) = VALUE
2010  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   AC1, AC2, AR1, AR2, BORDER
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID,
     +   AC1, SC1, SC2, AC2, AR1, SR1, SR2, AR2

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2'
     +   //10X, I10, '  BORDER')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  THE GLUE VALUE IS ZERO.')

99103 FORMAT
     +   (/1X, A9, 'ERROR.  THE GLUE BLOCK IS OFF TARGET.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  SC1'
     +    /10X, I10, '  SC2'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  SR1'
     +    /10X, I10, '  SR2'
     +    /10X, I10, '  AR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLPUT
     +  (ERROR, TEXT, PRINT,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   BCMAX, BCOL, BORDER, BRMAX, BROW)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLPUT
C
C>    CALLED BY THE  U S E R
C
C>    PUT A BORDER ON THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   ASIZE2, BCMAX, BLKS1, BLKS2, BORDER, BRMAX, COLS1, COLS2,
     +   GLUES, IDATA, IDSIZ, N1, N2, PRINT, RDSIZ, ROWS1, ROWS2, TEXT,
     +   XBRDR
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   BCOL, BROW, RDATA

      PARAMETER (ID = ' BLPUT:  ')

      DIMENSION
     +   BCOL(BCMAX, BORDER), BROW(BRMAX, BORDER), IDATA(IDSIZ),
     +   RDATA(RDSIZ)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  PRINT.

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) WRITE (TEXT, 10001) ID

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. BCMAX .AND. 0 .LT. BORDER .AND.
     +   0 .LT. BRMAX .AND. 0 .LT. IDSIZ .AND. 0 .LT. RDSIZ)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) CALL BLPU2.
C
C///////////////////////////////////////////////////////////////////////

C///  UNPACK THE DATA SPACES.

      CALL BLDAT
     +  (ERROR, TEXT,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   ASIZE2, BLKS1, BLKS2, XBRDR, COLS1, COLS2, GLUES, QA2, QAC1,
     +   QAC2, QACQ1, QACQ2, QAR2, QARQ2, QBMAP, QCMAP, QCOL1, QCOL2,
     +   QGMAP, QRMAP, QROW1, QROW2, QRPIV, QRSCAL, N1, N2, ROWS1,
     +   ROWS2)
      IF (ERROR) GO TO 9201

      ERROR = .NOT. (BORDER .EQ. XBRDR)
      IF (ERROR) GO TO 9202

      ERROR = .NOT. (N1 .LE. BCMAX .AND. N1 + BORDER .LE. BRMAX)
      IF (ERROR) GO TO 9203

C///  CALL.

C     SUBROUTINE BLPU2
C    +  (ERROR, TEXT,
C    +   A2, ACOL2, ACOLQ2, ASIZE2, BCMAX, BCOL, BLKS2, BMAP, BORDER,
C    +   BRMAX, BROW, CMAP, COL1, COL2, COLS1, COLS2, GLUES, GMAP, N1,
C    +   N2, RMAP, ROW1, ROW2, ROWS1, ROWS2, RSCALE)

2010  CONTINUE

      CALL BLPU2
     +  (ERROR, TEXT,
     +   RDATA(QA2), IDATA(QAC2), IDATA(QACQ2), ASIZE2, BCMAX, BCOL,
     +   BLKS2, IDATA(QBMAP), BORDER, BRMAX, BROW, IDATA(QCMAP),
     +   IDATA(QCOL1), IDATA(QCOL2), COLS1, COLS2, GLUES, IDATA(QGMAP),
     +   N1, N2, IDATA(QRMAP), IDATA(QROW1), IDATA(QROW2), ROWS1, ROWS2,
     +   RDATA(QRSCAL))
      IF (ERROR) GO TO 9204

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +   (/1X, A, 'PLACE THE BORDER AND GLUE IN A STRETCHED SPARSE'
     +   /10X, 'BLOCKED MATRIX.')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   BCMAX, BORDER, BRMAX, IDSIZ, RDSIZ
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID
      GO TO 99999

9202  IF (0 .LT. TEXT) WRITE (TEXT, 99204) ID, BORDER, XBRDR
      GO TO 99999

9203  IF (0 .LT. TEXT) WRITE (TEXT, 99203) ID,
     +   N1, BCMAX, N1 + BORDER, BRMAX
      GO TO 99999

9204  IF (0 .LT. TEXT) WRITE (TEXT, 99204) ID
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  BCMAX'
     +    /10X, I10, '  BORDER'
     +    /10X, I10, '  BRMAX'
     +    /10X, I10, '  IDSIZ'
     +    /10X, I10, '  RDSIZ')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  BLDAT FAILS.')

99202 FORMAT
     +   (/1X, A9, 'ERROR.  THE BORDER SIZE IS WRONG.'
     +   //10X, I10, '  BORDER SIZE IN ARGUMENT LIST'
     +    /10X, I10, '  EXPECTED SIZE')

99203 FORMAT
     +   (/1X, A9, 'ERROR.  THE LEADING DIMENSIONS OF BCOL OR BROW ARE'
     +    /10X, 'TOO SMALL.'
     +   //10X, I10, '  N1'
     +    /10X, I10, '  BCMAX'
     +   //10X, I10, '  N1 + BORDER'
     +    /10X, I10, '  BRMAX')

99204 FORMAT
     +   (/1X, A9, 'ERROR.  BLPU2 FAILS.')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLSH2
     +  (ERROR, TEXT,
     +   A2, ACOL1, ACOLQ1, ASIZE2, BLKS1, CMAP, COL1, COL2, COLS1,
     +   COLS2, COLUMN, RMAP, ROW, ROW1, ROW2, ROWS1, ROWS2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLSH2
C
C>    CALLED BY BLSHO
C
C>    PRINT A BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   AC1, AC2, ACOL1, ACOLQ1, AR1, AR2, ASIZE2, BC1, BC2, BLKS1,
     +   BLOCK, BR1, BR2, C2, CMAP, COFF, COL1, COL2, COLS1, COLS2,
     +   COLUMN, R2, RMAP, ROFF, ROW, ROW1, ROW2, ROWS1, ROWS2, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A2

      PARAMETER (ID = ' BLSH2:  ')

      PARAMETER (QFRST = 1, QLAST = 2, QFRSTD = 3, QLASTD = 4)
      PARAMETER (QIND = 1, QLOC = 2)

      DIMENSION
     +   A2(ASIZE2), ACOL1(2, BLKS1), ACOLQ1(2, COLS1), CMAP(COLS1),
     +   COL1(2, COLS1), COL2(2, COLS2), RMAP(ROWS1), ROW1(2, ROWS1),
     +   ROW2(2, ROWS2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. ASIZE2 .AND. 0 .LT. BLKS1 .AND. 0 .LT. COLS1
     +   .AND. 0 .LT. COLS2 .AND. 0 .LT. ROWS1 .AND. 0 .LT. ROWS2)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) PRINT THE BLOCK.
C
C///////////////////////////////////////////////////////////////////////

C///  SEARCH FOR THE BLOCK.

      ERROR = .NOT. (1 .LE. COLUMN .AND. COLUMN .LE. COLS1 .AND.
     +   1 .LE. ROW .AND. ROW .LE. ROWS1)
      IF (ERROR) GO TO 9201

      DO 2010 BLOCK = ACOLQ1(QFRST, COLUMN), ACOLQ1(QLAST, COLUMN)
         IF (ACOL1(QIND, BLOCK) .EQ. ROW) THEN
            QA2 = ACOL1(QLOC, BLOCK)
            GO TO 2020
         END IF
2010  CONTINUE
      ERROR = .TRUE.
      GO TO 9202

2020  CONTINUE

C///  PRINT THE BLOCK.

      C2 = CMAP(COLUMN)
      AC1 = COL2(QFRST, C2)
      AC2 = COL2(QLAST, C2)

      R2 = RMAP(ROW)
      AR1 = ROW2(QFRST, R2)
      AR2 = ROW2(QLAST, R2)

      BC1 = COL1(QFRST, COLUMN)
      BC2 = COL1(QLAST, COLUMN)

      BR1 = ROW1(QFRST, ROW)
      BR2 = ROW1(QLAST, ROW)

C     THE A1 BLOCKS LIE IN THE UPPER LEFT CORNERS OF THE A2 BLOCKS
      COFF = AC1 - BC1
      ROFF = AR1 - BR1

C     SUBROUTINE BLSH3
C    +  (ERROR, TEXT,
C    +   A, AC1, AC2, AR1, AR2, BC1, BC2, BR1, BR2, COFF, ROFF)

      CALL BLSH3
     +  (ERROR, TEXT,
     +   A2(QA2), AC1, AC2, AR1, AR2, BC1, BC2, BR1, BR2, COFF, ROFF)
      IF (ERROR) GO TO 9203

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   ASIZE2, BLKS1, COLS1, COLS2, ROWS1, ROWS2
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID,
     +   COLUMN, COLS1, ROW, ROWS1
      GO TO 99999

9202  IF (0 .LT. TEXT) WRITE (TEXT, 99202) ID,
     +   COLUMN, ROW, (ACOL1(QIND, BLOCK),
     +   BLOCK = ACOLQ1(QFRST, COLUMN), ACOLQ1(QLAST, COLUMN))
      GO TO 99999

9203  IF (0 .LT. TEXT) WRITE (TEXT, 99203) ID
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  ASIZE2'
     +    /10X, I10, '  BLKS1'
     +    /10X, I10, '  COLS1'
     +    /10X, I10, '  COLS2'
     +    /10X, I10, '  ROWS1'
     +    /10X, I10, '  ROWS2')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  THE BLOCK IS NOT IN THE MATRIX.'
     +   //10X, I10, '  BLOCK COLUMN'
     +    /10X, I10, '  COLUMNS'
     +   //10X, I10, '  BLOCK ROW'
     +    /10X, I10, '  ROWS')

99202 FORMAT
     +   (/1X, A9, 'ERROR.  THE BLOCK IS NOT IN THE MATRIX.'
     +   //10X, I10, '  BLOCK COLUMN'
     +    /10X, I10, '  BLOCK ROW'
     +   //10X, 'ROWS FOR THIS COLUMN'
     +  //(10X, 6I10))

99203 FORMAT
     +   (/1X, A9, 'ERROR.  BLSH3 FAILS.')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLSH3
     +  (ERROR, TEXT,
     +   A, AC1, AC2, AR1, AR2, BC1, BC2, BR1, BR2, COFF, ROFF)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLSH3
C
C>    CALLED BY BLSH2
C
C>    PRINT A BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9, WORD*80
      INTEGER
     +   AC1, AC2, AR1, AR2, BC1, BC2, BR1, BR2, COFF, COL, COUNT,
     +   FIRST, J, LAST, ROFF, ROW, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, TEMP

      PARAMETER (ID = ' BLSH3:  ')

      DIMENSION A(AR1 : AR2, AC1 : AC2), WORD(5)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (AC1 .LE. AC2 .AND. AR1 .LE. AR2)
      IF (ERROR) GO TO 9101

C///  CHECK THE PLACEMENT.

      ERROR = .NOT. (BC1 .LE. BC2 .AND. BR1 .LE. BR2)
      IF (ERROR) GO TO 9102

      ERROR = .NOT. (AC1 .LE. BC1 + COFF .AND. BC2 + COFF .LE. AC2 .AND.
     +   AR1 .LE. BR1 + ROFF .AND. BR2 + ROFF .LE. AR2)
      IF (ERROR) GO TO 9103

C///////////////////////////////////////////////////////////////////////
C
C     (2) PRINT THE BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IF (0 .LT. TEXT) THEN
         DO 2030 FIRST = BC1, BC2, 5
            LAST = MIN (FIRST + 4, BC2)

            WRITE (TEXT, 10001) (COL, COL = FIRST, LAST)
            WRITE (TEXT, 10002) ('____________', COL = FIRST, LAST)
            WRITE (TEXT, 10003)

            COUNT = LAST - FIRST + 1
            DO 2020 ROW = BR1, BR2
               DO 2010 J = 1, COUNT
                  TEMP = A(ROW + ROFF, J + FIRST - 1 + COFF)
                  IF (TEMP .EQ. 0.0) THEN
                     WORD(J) = '           0'
                  ELSE
                     WRITE (WORD(J), '(1P, E12.2)') TEMP
                  END IF
2010           CONTINUE
               WORD(1) (1 : 1) = '|'
               WRITE (TEXT, 10004) ROW, (WORD(J), J = 1, COUNT)
2020        CONTINUE
2030     CONTINUE
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT (/19X, 5I12)

10002 FORMAT (20X, A11, 4A12)

10003 FORMAT (19X, '|')

10004 FORMAT (6X, I12, 1X, 5A12)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, AC1, AC2, AR1, AR2
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, BC1, BC2, BR1, BR2
      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID,
     +   BC1, BC2, COFF, BC1 + COFF, BC2 + COFF, AC1, AC2, BR1, BR2,
     +   ROFF, BR1 + ROFF, BR2 + ROFF, AR1, AR2

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  BC1'
     +    /10X, I10, '  BC2'
     +   //10X, I10, '  BR1'
     +    /10X, I10, '  BR2')

99103 FORMAT
     +   (/1X, A9, 'ERROR.  THE DATA BLOCK DOES NOT CONTAIN THE PRINT'
     +    /10X, 'BLOCK.'
     +   //10X, I10, '  BC1'
     +    /10X, I10, '  BC2'
     +   //10X, I10, '  COFF'
     +   //10X, I10, '  BC1 + COFF'
     +    /10X, I10, '  BC2 + COFF'
     +   //10X, I10, '  AC1'
     +    /10X, I10, '  AC2'
     +   //10X, I10, '  BR1'
     +    /10X, I10, '  BR2'
     +   //10X, I10, '  ROFF'
     +   //10X, I10, '  BR1 + ROFF'
     +    /10X, I10, '  BR2 + ROFF'
     +   //10X, I10, '  AR1'
     +    /10X, I10, '  AR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLSHO
     +  (ERROR, TEXT,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   COLUMN, ROW)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLSHO
C
C>    CALLED BY THE  U S E R
C
C>    PRINT A BLOCK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   ASIZE2, BLKS1, BLKS2, BORDER, COLS1, COLS2, COLUMN, GLUES,
     +   IDATA, IDSIZ, N1, N2, RDSIZ, ROW, ROWS1, ROWS2, TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   RDATA

      PARAMETER (ID = ' BLSHO:  ')

      DIMENSION IDATA(IDSIZ), RDATA(RDSIZ)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  PRINT.

      IF (0 .LT. TEXT) WRITE (TEXT, 10001) ID, COLUMN, ROW

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. IDSIZ .AND. 0 .LT. RDSIZ)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) CALL BLSH2.
C
C///////////////////////////////////////////////////////////////////////

C///  UNPACK THE DATA SPACES.

      CALL BLDAT
     +  (ERROR, TEXT,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   ASIZE2, BLKS1, BLKS2, BORDER, COLS1, COLS2, GLUES, QA2, QAC1,
     +   QAC2, QACQ1, QACQ2, QAR2, QARQ2, QBMAP, QCMAP, QCOL1, QCOL2,
     +   QGMAP, QRMAP, QROW1, QROW2, QRPIV, QRSCAL, N1, N2, ROWS1,
     +   ROWS2)
      IF (ERROR) GO TO 9201

C///  CALL.

C     SUBROUTINE BLSH2
C    +  (ERROR, TEXT,
C    +   A2, ACOL1, ACOLQ1, ASIZE2, BLKS1, CMAP, COL1, COL2, COLS1,
C    +   COLS2, COLUMN, RMAP, ROW, ROW1, ROW2, ROWS1, ROWS2)

2010  CONTINUE

      CALL BLSH2
     +  (ERROR, TEXT,
     +   RDATA(QA2), IDATA(QAC1), IDATA(QACQ1), ASIZE2, BLKS1,
     +   IDATA(QCMAP), IDATA(QCOL1), IDATA(QCOL2), COLS1, COLS2, COLUMN,
     +   IDATA(QRMAP), ROW, IDATA(QROW1), IDATA(QROW2), ROWS1, ROWS2)
      IF (ERROR) GO TO 9204

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +   (/1X, A, 'PRINT A BLOCK.'
     +   //10X, I10, '  COLUMN BLOCK'
     +    /10X, I10, '  ROW BLOCK')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, IDSIZ, RDSIZ
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID
      GO TO 99999

9204  IF (0 .LT. TEXT) WRITE (TEXT, 99204) ID
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  IDSIZ'
     +    /10X, I10, '  RDSIZ')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  BLDAT FAILS.')

99204 FORMAT
     +   (/1X, A9, 'ERROR.  BLSH2 FAILS.')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLSO2
     +  (ERROR, TEXT,
     +   A2, AROW2, AROWQ2, ASIZE2, BLKS2, BORDER, CMAP, COL1, COL2,
     +   COLS1, COLS2, N1, N2, RMAP, ROW1, ROW2, ROWS1, ROWS2, RPIVOT,
     +   RSCALE, X2, YB)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLSO2
C
C>    CALLED BY BLSOL
C
C>    SOLVE EQUATIONS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   AROW2, AROWQ2, ASIZE2, BLKS2, BORDER, C1, C2, CMAP, COL1, COL2,
     +   COLS1, COLS2, J, LC1, LC2, LR1, LR2, N1, N2, OFFSET, R1, R2,
     +   RMAP, ROW1, ROW2, ROWS1, ROWS2, RPIVOT, TEXT, UC1, UC2, UR1,
     +   UR2
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A2, RSCALE, X2, YB

      PARAMETER (ID = ' BLSO2:  ')

      PARAMETER (QFRST = 1, QLAST = 2, QFRSTD = 3, QLASTD = 4)
      PARAMETER (QIND = 1, QLOC = 2)

      DIMENSION
     +   A2(ASIZE2), AROW2(2, BLKS2), AROWQ2(4, ROWS2), CMAP(COLS1),
     +   COL1(2, COLS1), COL2(2, COLS2), RMAP(ROWS1), ROW1(2, ROWS1),
     +   ROW2(2, ROWS2), RPIVOT(N2), X2(N2), RSCALE(N2), YB(N1 + BORDER)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. ASIZE2 .AND. 0 .LT. BLKS2 .AND.
     +   0 .LE. BORDER .AND. 0 .LT. COLS1 .AND. 0 .LT. COLS2 .AND.
     +   0 .LT. N1 .AND. 0 .LT. N2 .AND. 0 .LT. ROWS2)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) STRETCH AND SCALE THE RIGHT SIDE.
C
C///////////////////////////////////////////////////////////////////////

C///  CLEAR THE RIGHT SIDE.

      DO 2010 J = 1, N2
         X2(J) = 0.0
2010  CONTINUE

C///  COPY THE A1 BLOCKS.

      DO 2030 R1 = 1, ROWS1
         R2 = RMAP(R1)
         OFFSET = ROW2(QFRST, R2) - ROW1(QFRST, R1)
         DO 2020 J = ROW1(QFRST, R1), ROW1(QLAST, R1)
            X2(OFFSET + J) = YB(J)
2020     CONTINUE
2030  CONTINUE

C///  COPY THE BORDER BLOCK.

      DO 2040 J = 1, BORDER
         X2(N2 - BORDER + J) = YB(N1 + J)
2040  CONTINUE

C///  SCALE THE RIGHT SIDE.

      DO 2050 J = 1, N2
         X2(J) = RSCALE(J) * X2(J)
2050  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (3) BACKWARD SUBSTITUTION.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP OVER THE BLOCKS OF ROWS AND COLUMNS.

      DO 3020 R2 = 1, ROWS2
         LR1 = ROW2(QFRST, R2)
         LR2 = ROW2(QLAST, R2)
         DO 3010 QC2 = AROWQ2(QFRST, R2), AROWQ2(QLASTD, R2)
            C2 = AROW2(QIND, QC2)
            LC1 = COL2(QFRST, C2)
            LC2 = COL2(QLAST, C2)
            QA2 = AROW2(QLOC, QC2)

C///  PERFORM THE BACKWARD SUBSTITION.

C     SUBROUTINE BLSO3
C    +  (ERROR, TEXT,
C    +   L, LC1, LC2, LR1, LR2, N2, RPIVOT, X)

      CALL BLSO3
     +  (ERROR, TEXT,
     +   A2(QA2), LC1, LC2, LR1, LR2, N2, RPIVOT, X2)
      IF (ERROR) GO TO 9301

C///  BOTTOM OF THE LOOP OVER THE BLOCKS OF ROWS AND COLUMNS.

3010     CONTINUE
3020  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (4) FORWARD SUBSTITUTION.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP OVER THE BLOCKS OF ROWS AND COLUMNS.

      DO 4020 R2 = ROWS2, 1, - 1
         UR1 = ROW2(QFRST, R2)
         UR2 = ROW2(QLAST, R2)
         DO 4010 QC2 = AROWQ2(QLAST, R2), AROWQ2(QFRSTD, R2), - 1
            C2 = AROW2(QIND, QC2)
            UC1 = COL2(QFRST, C2)
            UC2 = COL2(QLAST, C2)
            QA2 = AROW2(QLOC, QC2)

C///  PERFORM THE FORWARD SUBSTITION.

C     SUBROUTINE BLSO4
C    +  (ERROR, TEXT,
C    +   N2, U, UC1, UC2, UR1, UR2, X)

      CALL BLSO4
     +  (ERROR, TEXT,
     +   N2, A2(QA2), UC1, UC2, UR1, UR2, X2)
      IF (ERROR) GO TO 9401

C///  BOTTOM OF THE LOOP OVER THE BLOCKS OF ROWS AND COLUMNS.

4010     CONTINUE
4020  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (5) SQUEEZE THE SOLUTION.
C
C///////////////////////////////////////////////////////////////////////

C///  COPY THE A1 BLOCKS.

      DO 5020 C1 = 1, COLS1
         C2 = CMAP(C1)
         OFFSET = COL2(QFRST, C2) - COL1(QFRST, C1)
         DO 5010 J = COL1(QFRST, C1), COL1(QLAST, C1)
            YB(J) = X2(OFFSET + J)
5010     CONTINUE
5020  CONTINUE

C///  COPY THE BORDER BLOCK.

      DO 5030 J = 1, BORDER
         YB(N1 + J) = X2(N2 - BORDER + J)
5030  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   ASIZE2, BLKS2, BORDER, COLS1, COLS2, N1, N2, ROWS2
      GO TO 99999

9301  IF (0 .LT. TEXT) WRITE (TEXT, 99301) ID
      GO TO 99999

9401  IF (0 .LT. TEXT) WRITE (TEXT, 99401) ID
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  ASIZE2'
     +    /10X, I10, '  BLKS2'
     +    /10X, I10, '  BORDER (MAY BE ZERO)'
     +    /10X, I10, '  COLS1'
     +    /10X, I10, '  COLS2'
     +    /10X, I10, '  N1'
     +    /10X, I10, '  N2'
     +    /10X, I10, '  ROWS2')

99301 FORMAT
     +   (/1X, A9, 'ERROR.  BLSO3 FAILS.')

99401 FORMAT
     +   (/1X, A9, 'ERROR.  BLSO4 FAILS.')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLSO3
     +  (ERROR, TEXT,
     +   L, LC1, LC2, LR1, LR2, N2, RPIVOT, X2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLSO3
C
C>    CALLED BY BLSO2
C
C>    PERFORM BACKWARD SUBSTITUTION.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   FIRST, LAST, LC, LC1, LC2, LR, LR1, LR2, N2, RINDEX, RPIVOT,
     +   TEXT
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   L, TEMP, X2

      PARAMETER (ID = ' BLSO3:  ')

      DIMENSION L(LR1 : LR2, LC1 : LC2), RPIVOT(N2), X2(N2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (LC1 .LE. LC2 .AND. LR1 .LE. LR2 .AND. 0 .LT. N2)
      IF (ERROR) GO TO 9101

C///  CHECK THE RANGE.

      ERROR = .NOT. (LC1 .LE. LR2)
      IF (ERROR) GO TO 9102

C///////////////////////////////////////////////////////////////////////
C
C     (2) PERFORM THE ELIMINATION.
C
C///////////////////////////////////////////////////////////////////////

      FIRST = MAX (LC1, LR1)
      LAST = MIN (LC2, LR2)

      DO 2020 LC = LC1, MIN (FIRST - 1, LC2)
         TEMP = X2(LC)
         DO 2010 LR = LR1, LR2
            X2(LR) = X2(LR) - TEMP * L(LR, LC)
2010     CONTINUE
2020  CONTINUE

      DO 2040 LC = FIRST, LAST
         RINDEX = RPIVOT(LC)
         TEMP = X2(RINDEX)
         X2(RINDEX) = X2(LC)
         X2(LC) = TEMP
         DO 2030 LR = LC + 1, LR2
            X2(LR) = X2(LR) - TEMP * L(LR, LC)
2030     CONTINUE
2040  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, LC1, LC2, LR1, LR2, N2
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, LC1, LR2
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  LC1'
     +    /10X, I10, '  LC2'
     +   //10X, I10, '  LR1'
     +    /10X, I10, '  LR2'
     +   //10X, I10, '  N2')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  THE BLOCK LIES ABOVE THE MAIN DIAGONAL.'
     +   //10X, I10, '  LC1'
     +    /10X, I10, '  LR2')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLSO4
     +  (ERROR, TEXT,
     +   N2, U, UC1, UC2, UR1, UR2, X2)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLSO4
C
C>    CALLED BY BLSO2
C
C>    PERFORM FORWARD SUBSTITUTION.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER N2, TEXT, UC, UC1, UC2, UR, UR1, UR2
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   TEMP, U, X2

      PARAMETER (ID = ' BLSO4:  ')

      DIMENSION U(UR1 : UR2, UC1 : UC2), X2(N2)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LT. N2 .AND. UC1 .LE. UC2 .AND. UR1 .LE. UR2)
      IF (ERROR) GO TO 9101

C///  CHECK THE RANGE.

      ERROR = .NOT. (UR1 .LE. UC2)
      IF (ERROR) GO TO 9102

C///////////////////////////////////////////////////////////////////////
C
C     (2) PERFORM THE ELIMINATION.
C
C///////////////////////////////////////////////////////////////////////

      DO 2020 UC = UC2, MAX (UR2 + 1, UC1), - 1
         TEMP = X2(UC)
         DO 2010 UR = UR1, UR2
            X2(UR) = X2(UR) - U(UR, UC) * TEMP
2010     CONTINUE
2020  CONTINUE

      DO 2040 UC = MIN (UC2, UR2), MAX (UC1, UR1), - 1
         TEMP = U(UC, UC) * X2(UC)
         X2(UC) = TEMP
         DO 2030 UR = UR1, UC - 1
            X2(UR) = X2(UR) - U(UR, UC) * TEMP
2030     CONTINUE
2040  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, N2, UC1, UC2, UR1, UR2
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, UC2, UR1
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE INCONSISTENT.'
     +   //10X, I10, '  N2'
     +   //10X, I10, '  UC1'
     +    /10X, I10, '  UC2'
     +   //10X, I10, '  UR1'
     +    /10X, I10, '  UR2')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  THE BLOCK LIES BELOW THE MAIN DIAGONAL.'
     +   //10X, I10, '  UC2'
     +    /10X, I10, '  UR1')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE BLSOL
     +  (ERROR, TEXT, PRINT,
     +   RLAST, RMAX, RSIZE, RWORK,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   BORDER, N, YB)

C///////////////////////////////////////////////////////////////////////
C
C     B L O C K
C
C>    BLSOL
C
C>    CALLED BY THE  U S E R
C
C>    SOLVE EQUATIONS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   ASIZE2, BLKS1, BLKS2, BORDER, COLS1, COLS2, GLUES, IDATA,
     +   IDSIZ, N, N1, N2, PRINT, ROWS1, ROWS2, RDSIZ, RLAST, RMAX,
     +   RSAVE, RSIZE, TEXT, XBRDR
      LOGICAL ERROR
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   RDATA, RWORK, YB

      PARAMETER (ID = ' BLSOL:  ')

      DIMENSION IDATA(IDSIZ), RDATA(RDSIZ), RWORK(RSIZE), YB(N + BORDER)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  SAVE THE LEVELS OF THE WORK SPACES.

      RSAVE = RLAST

C///  PRINT.

      IF (0 .LT. PRINT .AND. 0 .LT. TEXT) WRITE (TEXT, 10001) ID

C///  CHECK THE DIMENSIONAL ARGUMENTS.

      ERROR = .NOT. (0 .LE. BORDER .AND. 0 .LT. IDSIZ .AND.
     +   0 .LT. N .AND. 0 .LT. RDSIZ .AND. 0 .LT. RSIZE)
      IF (ERROR) GO TO 9101

C///////////////////////////////////////////////////////////////////////
C
C     (2) CALL BLSO2.
C
C///////////////////////////////////////////////////////////////////////

C///  UNPACK THE DATA SPACES.

      CALL BLDAT
     +  (ERROR, TEXT,
     +   IDATA, IDSIZ, RDATA, RDSIZ,
     +   ASIZE2, BLKS1, BLKS2, XBRDR, COLS1, COLS2, GLUES, QA2, QAC1,
     +   QAC2, QACQ1, QACQ2, QAR2, QARQ2, QBMAP, QCMAP, QCOL1, QCOL2,
     +   QGMAP, QRMAP, QROW1, QROW2, QRPIV, QRSCAL, N1, N2, ROWS1,
     +   ROWS2)
      IF (ERROR) GO TO 9201

      ERROR = .NOT. (N .EQ. N1 .AND. BORDER .EQ. XBRDR)
      IF (ERROR) GO TO 9202

C///  RESERVE SOME WORK SPACE.

C     SUBROUTINE RESERV
C    +  (ERROR, TEXT,
C    +   FATAL, QLAST, QMAX, QSIZE,
C    +   NAME, NUMBER, Q)

C     REAL Y2(N2)

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QX2', N2, QX2)
      IF (ERROR) GO TO 9203

      ERROR = .NOT. (RLAST .LE. RSIZE)
      IF (ERROR) GO TO 9204

C///  CALL.

C     SUBROUTINE BLSO2
C    +  (ERROR, TEXT,
C    +   A2, AROW2, AROWQ2, ASIZE2, BLKS2, BORDER, CMAP, COL1, COL2,
C    +   COLS1, COLS2, N1, N2, RMAP, ROW1, ROW2, ROWS1, ROWS2, RPIVOT,
C    +   RSCALE, X2, YB)

2010  CONTINUE

      CALL BLSO2
     +  (ERROR, TEXT,
     +   RDATA(QA2), IDATA(QAR2), IDATA(QARQ2), ASIZE2, BLKS2, BORDER,
     +   IDATA(QCMAP), IDATA(QCOL1), IDATA(QCOL2), COLS1, COLS2,
     +   N1, N2, IDATA(QRMAP), IDATA(QROW1), IDATA(QROW2), ROWS1, ROWS2,
     +   IDATA(QRPIV), RDATA(QRSCAL), RWORK(QX2), YB)
      IF (ERROR) GO TO 9205

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +   (/1X, A, 'SOLVE A SYSTEM OF EQUATIONS WITH A FACTORED SPARSE'
     +   /10X, 'BLOCKED MATRIX.')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID,
     +   BORDER, IDSIZ, N, RDSIZ, RSIZE
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID
      GO TO 99999

9202  IF (0 .LT. TEXT) WRITE (TEXT, 99202) ID, N, N1, BORDER, XBRDR
      GO TO 99999

9203  IF (0 .LT. TEXT) WRITE (TEXT, 99203) ID
      GO TO 99999

9204  IF (0 .LT. TEXT) WRITE (TEXT, 99204) ID,
     +   0, 0, 0, RSIZE, 0, 0, 0, RLAST
      GO TO 99999

9205  IF (0 .LT. TEXT) WRITE (TEXT, 99205) ID
      GO TO 99999

99101 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +   //10X, I10, '  BORDER'
     +    /10X, I10, '  IDSIZ'
     +    /10X, I10, '  N',
     +    /10X, I10, '  RDSIZ'
     +    /10X, I10, '  RSIZE')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  BLDAT FAILS.')

99202 FORMAT
     +   (/1X, A9, 'ERROR.  THE MATRIX ORDER IS INCONSISTENT.'
     +   //10X, I10,  '  ORDER IN ARGUMENT LIST'
     +    /10X, I10,  '  IN DATA FROM BLSYM'
     +   //10X, I10,  '  BORDER IN ARGUMENT LIST'
     +    /10X, I10,  '  IN DATA FROM BLSYM')

99203 FORMAT
     +   (/1X, A9, 'ERROR.  RESERV FAILS.')

99204 FORMAT
     +   (/1X, A9, 'ERROR.  SOME WORKSPACES ARE TOO SMALL.'
C               123456789_  123456789_  123456789_  123456789_
     +  //25X, ' CHARACTER     INTEGER     LOGICAL        REAL'
C               123456789_12345
     +  //10X, ' PRESENT SIZE', 4I12
     +   /10X, 'REQUIRED SIZE', 4I12
     +  //10X, 'THIS SUBROUTINE WILL NEED NO MORE SPACE.')

99205 FORMAT
     +   (/1X, A9, 'ERROR.  BLSO2 FAILS.')

C///  EXIT.

99999 CONTINUE
      RLAST = RSAVE
      RETURN
      END
