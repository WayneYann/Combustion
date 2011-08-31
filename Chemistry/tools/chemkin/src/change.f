C  CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:32 $
C///////////////////////////////////////////////////////////////////////
C
C     C H A N G E
C
C     VERSION 2.05 OF MARCH 1995
C
C     COMMENT-OUT OR UN-COMMENT-IN LINES OF A FORTRAN PROGRAM OR A UNIX
C     SHELL SCRIPT
C
C     WRITTEN BY:  DR. JOSEPH F. GRCAR
C                  SANDIA NATIONAL LABORATORIES
C                  DEPARTMENT 8745
C                  LIVERMORE, CA 94551-0969 USA
C
C                  (510) 294-2662
C
C                  na.grcar@na-net.ornl.gov
C                  sepp@sandia.california.gov
C
C///////////////////////////////////////////////////////////////////////
C
C     DOCUMENTATION:
C
C     J. F. Grcar, "The Change Tool for Changing Programs and Scripts,"
C     Sandia National Laboratories Report SAND92-8225C, Livermore,
C     California, December 1993.
C
C///////////////////////////////////////////////////////////////////////
C
C     CHANGES FROM THE PREVIOUS VERSION:
C
C     1) CORRECT ERROR 9106 ABOUT UNAVAILABLE COPY FILES.
C
C     2) CORRECT ERROR 9208 ABOUT INCONSISTENT CHANGE BLOCKS.
C
C     3) REWRITE THE VARIABLE DECLARATIONS.
C
C///////////////////////////////////////////////////////////////////////

      PROGRAM CHANGE

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)

C     THE FOLLOWING CHARACTER STRINGS SHOULD HAVE LENGTH EQUAL TO THE
C     MAXS PARAMETER:  LINE, NAME, SAVLIN, TEMP, TEMP2, TEMP3, TOP

      CHARACTER
     +   FILE1*80, FILE2*80, ITEM*8, LINE*200, LIST*80, MARK*8,
     +   NAME*200, REMARK*6, SAVLIN*200, SMARK*8, TEMP*200, TEMP2*200,
     +   TEMP3*200, TOP*200
      EXTERNAL
     +   EXTENT, SQUEEZ
      INTEGER
     +   CMND, CMNTS, COUNT, INDEX, ITEMP, J, K, LEN1, LEN2, LEN3, LEN4,
     +   LEN5, LENGTH, LINES, LONG, MAXL, MAXN, MAXS, NAMES, NUMBER,
     +   QADD, QCMND, QCOPY, QFORT, QSHELL, QSTRIP, REPORT, SAVLEN,
     +   SAVNUM, STATUS, TOPNUM, TYPE, UNIT1, UNIT2
      INTRINSIC
     +   LEN, MOD
      LOGICAL
     +   ACTIVE, ERROR, EXIST, FLAG, FOUND, INSIDE, MANUAL, MAY, NEW,
     +   OPEN, SAVED, SPECL, VARY

      PARAMETER (REPORT = 1000, MAXN = 100)
      PARAMETER (QADD = 1, QCOPY = 2, QSTRIP = 3)
      PARAMETER (QFORT = 1, QSHELL = 2)
      PARAMETER (MAXS = 200)

      DIMENSION
     +   ACTIVE(MAXN, 2), CMND(0 : MAXN), LONG(MAXN), NAME(MAXN),
     +   REMARK(MAXN, 2), SMARK(MAXN), TYPE(0 : MAXN)

C     FORCE STATIC ALLOCATION OF VARIABLES.  THIS IS NOT NEEDED BY THE
C     CHANGE PROGRAM, BUT IT SEEMS TO HELP WITH TROUBLESOME COMPILERS.

      SAVE

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHOOSE THE UNIT NUMBERS.

      UNIT1 = 20
      UNIT2 = 21

C///  DECIDE WHETHER THE COPY MAY OVERWRITE THE ORIGINAL.

      MAY = .FALSE.
C*****COPY MAY OVERWRITE ORIGINAL (RECOMMENDED)
      MAY = .TRUE.
C*****END COPY MAY OVERWRITE ORIGINAL (RECOMMENDED)

C///  DECIDE WHETHER TO LIMIT LINES TO 72 CHARACTERS.

C*****LINE LENGTH > AS LARGE AS POSSIBLE
C      MAXL = MAXS
C      VARY = .FALSE.
C*****END LINE LENGTH > AS LARGE AS POSSIBLE

C*****LINE LENGTH > 72 CHARACTERS (FORTRAN 77 STANDARD)
C      MAXL = 72
C      VARY = .FALSE.
C*****END LINE LENGTH > 72 CHARACTERS (FORTRAN 77 STANDARD)

C*****LINE LENGTH > 72 OR 73
      MAXL = 72
      VARY = .TRUE.
C*****END LINE LENGTH > 72 OR 73

C///  WRITE THE HEADER.

      WRITE (*, 10001)

C///  PRINT ERROR MESSAGES FOR THE USER MANUAL.

      MANUAL = .FALSE.
C      MANUAL = .TRUE.

      IF (MANUAL) THEN
         LINE = 'C*****ABRAHAM LINCOLN'
         NAME(1) = 'GEORGE WASHINGTON'
         NUMBER = 8701
         SAVLIN = '     +  (/13X, ''ERROR.  IN THE ACTIVE CHANGE BLOCK '
     +      // 'WHICH BEGINS ON LINE '', A, '','''
         SAVNUM = 4528
         STATUS = 43
         TOP = 'C*****GEORGE WASHINGTON'
         TOPNUM = 4521
         INDEX = 1
         NAMES = 1
         GO TO 9101
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (1) OPEN THE ORIGINAL FILE.
C
C///////////////////////////////////////////////////////////////////////

1010  CONTINUE

C///  PROMPT FOR THE FILE NAME.

      WRITE (*, 10002)
      READ
     +  (END = 99999,
     +   ERR = 9101,
     +   FMT = '(A)',
     +   IOSTAT = STATUS,
     +   UNIT = *)
     +   FILE1

      CALL SQUEEZ (LENGTH, FILE1)
      IF (FILE1 .EQ. ' ') GO TO 99999

C///  CHECK THE FILE.

      ERROR = .TRUE.
      INQUIRE
     +  (ERR = 9102,
     +   EXIST = EXIST,
     +   FILE = FILE1,
     +   IOSTAT = STATUS,
     +   OPENED = OPEN)
      ERROR = .FALSE.

      ERROR = .NOT. EXIST
      IF (ERROR) GO TO 9103

      ERROR = EXIST .AND. OPEN
      IF (ERROR) GO TO 9104

C///  OPEN THE FILE.

      ERROR = .TRUE.
      OPEN
     +  (ACCESS = 'SEQUENTIAL',
     +   ERR = 9105,
     +   FILE = FILE1,
     +   FORM = 'FORMATTED',
     +   IOSTAT = STATUS,
     +   STATUS = 'OLD',
     +   UNIT = UNIT1)
      ERROR = .FALSE.

C///////////////////////////////////////////////////////////////////////
C
C     (2) LOOP THROUGH THE ORIGINAL FILE.
C
C///////////////////////////////////////////////////////////////////////

      INSIDE = .FALSE.
      MARK = ' '
      NUMBER = 0
      NAMES = 0

C///  TOP OF THE LOOP.  READ THE NEXT LINE.

2010  CONTINUE
      READ
     +  (END = 2050,
     +   ERR = 9201,
     +   FMT = '(A)',
     +   IOSTAT = STATUS,
     +   UNIT = UNIT1)
     +   LINE
      NUMBER = NUMBER + 1

C///  GATHER INFORMATION ABOUT THE LINE.

      LENGTH = 0
      DO 2020 J = MAXS, 1, - 1
         IF (LINE (J : J) .NE. ' ') THEN
            LENGTH = J
            GO TO 2030
         END IF
2020  CONTINUE
2030  CONTINUE

      SPECL = LINE (1 : 6) .EQ. 'C*****'
     +   .OR. LINE (1 : 6) .EQ. 'c*****'
     +   .OR. LINE (1 : 6) .EQ. '******'
     +   .OR. LINE (1 : 6) .EQ. '#*****'

C///  TREAT SPECIAL LINES.

      IF (SPECL) THEN
         ERROR = MAXL .LT. LENGTH
         IF (ERROR) GO TO 9202

C///  TREAT SPECIAL LINES FOUND FROM OUTSIDE.

         IF (LINE (1 : 1) .EQ. '#') THEN
            TYPE(0) = QSHELL
         ELSE
            TYPE(0) = QFORT
         END IF

         IF (.NOT. INSIDE) THEN
            FOUND = .FALSE.
            DO 2040 J = 1, NAMES
               IF (TYPE(0) .EQ. TYPE(J) .AND.
     +            LINE (7 : ) .EQ. NAME(J)) THEN
                  FOUND = .TRUE.
                  INDEX = J
               END IF
2040        CONTINUE

            IF (.NOT. FOUND) THEN
               NAMES = NAMES + 1
               ERROR = NAMES .GT. MAXN
               IF (ERROR) GO TO 9203

               NEW = .TRUE.
               INDEX = NAMES
               NAME(INDEX) = LINE (7 : )
               TYPE(INDEX) = TYPE(0)
            END IF

            CMNTS = 0
            INSIDE = .TRUE.
            LINES = 0
            MARK = LINE (1 : 1)
            SAVED = .FALSE.
            TOP = LINE
            TOPNUM = NUMBER

C///  TREAT SPECIAL LINES FOUND FROM INSIDE.

         ELSE
            ERROR = .NOT. (
     +         (LINE (7 : 7) .EQ. 'E' .OR. LINE (7 : 7) .EQ. 'e') .AND.
     +         (LINE (8 : 8) .EQ. 'N' .OR. LINE (8 : 8) .EQ. 'n') .AND.
     +         (LINE (9 : 9) .EQ. 'D' .OR. LINE (9 : 9) .EQ. 'd') .AND.
     +          LINE (10 : 10) .EQ. ' ')
            IF (ERROR) GO TO 9204

            ERROR = .NOT. (NAME(INDEX) .EQ. LINE (11 : ))
            IF (ERROR) GO TO 9205

            ERROR = .NOT. (TOP (1 : 1) .EQ. LINE (1 : 1))
            IF (ERROR) GO TO 9206

            ERROR = LINES .EQ. 0
            IF (ERROR) GO TO 9207

            FLAG = CMNTS .LT. LINES
            IF (NEW) THEN
               ACTIVE(INDEX, 1) = FLAG
               NEW = .FALSE.
            ELSE
               ERROR = .NOT. (ACTIVE(INDEX, 1) .EQV. FLAG)
               IF (ERROR) GO TO 9208
            END IF

            IF (SAVED) THEN
               IF (VARY .EQV. ACTIVE(INDEX, 1)) THEN
                  ERROR = MAXL .LT. SAVLEN
                  IF (ERROR) GO TO 9209
               ELSE IF (.NOT. ACTIVE(INDEX, 1)) THEN
                  ERROR = MAXL + 1 .LT. SAVLEN
                  IF (ERROR) GO TO 9210
               ELSE
                  ERROR = MAXL - 1 .LT. SAVLEN
                  IF (ERROR) GO TO 9211
               END IF
            END IF

            INDEX = 0
            INSIDE = .FALSE.
         END IF

C///  TREAT SIMPLE LINES FROM OUTSIDE.

      ELSE
         IF (.NOT. INSIDE) THEN
            ERROR = MAXL .LT. LENGTH
            IF (ERROR) GO TO 9202

C///  TREAT SIMPLE LINES FROM INSIDE.

         ELSE
            LINES = LINES + 1

            IF (TYPE (INDEX) .EQ. QFORT) THEN
               IF (LINE (1 : 1) .EQ. 'C' .OR. LINE (1 : 1) .EQ. 'c'
     +            .OR. LINE (1 : 1) .EQ. '*') CMNTS = CMNTS + 1
            ELSE IF (TYPE (INDEX) .EQ. QSHELL) THEN
               IF (LINE (1 : 1) .EQ. '#') CMNTS = CMNTS + 1
            END IF

            IF (SAVED) THEN
               FLAG = SAVLEN .LT. LENGTH
            ELSE
               FLAG = .TRUE.
            END IF
            IF (FLAG) THEN
               SAVED = .TRUE.
               SAVLEN = LENGTH
               SAVLIN = LINE
               SAVNUM = NUMBER
            END IF
         END IF
      END IF

C///  BOTTOM OF THE LOOP THROUGH THE FILE.

      IF (MOD (NUMBER, REPORT) .EQ. 0) THEN
         IF (NUMBER .EQ. REPORT) WRITE (*, 10004)
         WRITE (TEMP, 10005) NUMBER
         CALL SQUEEZ (LENGTH, TEMP)
         WRITE (*, 10006) TEMP (1 : LENGTH)
      END IF

      GO TO 2010
2050  CONTINUE

      IF (MOD (NUMBER, REPORT) .NE. 0) THEN
         IF (NUMBER .LT. REPORT) WRITE (*, 10004)
         WRITE (TEMP, 10005) NUMBER
         CALL SQUEEZ (LENGTH, TEMP)
         WRITE (*, 10006) TEMP (1 : LENGTH)
      END IF

      ERROR = INSIDE
      IF (ERROR) GO TO 9212

C///////////////////////////////////////////////////////////////////////
C
C     (3) DETERMINE WHICH BLOCKS SHOULD BE ACTIVE.
C
C///////////////////////////////////////////////////////////////////////

C///  DECIDE WHETHER TO WRITE A COPY.

      IF (NAMES .EQ. 0) THEN
         WRITE (*, 10010)
         GO TO 99999
      END IF

C///  SORT THE CHANGE BLOCKS.

      DO 3010 INDEX = 1, NAMES
         IF (TYPE(INDEX) .EQ. QFORT) THEN
            SMARK(INDEX) = ')'
         ELSE
            SMARK(INDEX) = '#'
         END IF
3010  CONTINUE

      DO 3020 K = 1, NAMES
         DO 3020 J = K + 1, NAMES
            IF (NAME(J) .LT. NAME(K)) THEN
               ITEMP = TYPE(J)
               TYPE(J) = TYPE(K)
               TYPE(K) = ITEMP

               TEMP = SMARK(J)
               SMARK(J) = SMARK(K)
               SMARK(K) = TEMP

               TEMP = NAME(J)
               NAME(J) = NAME(K)
               NAME(K) = TEMP

               FLAG = ACTIVE(J, 1)
               ACTIVE(J, 1) = ACTIVE(K, 1)
               ACTIVE(K, 1) = FLAG
            END IF
3020  CONTINUE

C///  INITIALIZE THE CHOICES AND LENGTHS.

      DO 3030 J = 1, NAMES
         ACTIVE(J, 2) = ACTIVE(J, 1)
         CALL EXTENT (LONG(J), NAME(J))
3030  CONTINUE

C///  TOP OF THE QUERY LOOP.

      LIST = ' '
      WRITE (*, 10007)
3040  CONTINUE

C///  DISPLAY THE NAMES OF THE CHANGE BLOCKS.  ASK WHICH TO CHANGE.

      IF (LIST .EQ. ' ') THEN
         DO 3050 K = 1, 2
            DO 3050 J = 1, NAMES
                 IF (ACTIVE(J, K)) THEN
                    REMARK(J, K) = 'ACTIVE'
                 ELSE
                    REMARK(J, K) = '  --  '
                 END IF
3050     CONTINUE

         WRITE (*, 10008) (J, SMARK(J), REMARK(J, 1), REMARK(J, 2),
     +      NAME(J) (1 : LONG(J)), J = 1, NAMES)

         WRITE (*, 10009)
         READ
     +     (END = 3090,
     +      ERR = 9101,
     +      FMT = '(A)',
     +      IOSTAT = STATUS,
     +      UNIT = *)
     +      LIST
         IF (LIST .EQ. ' ') GO TO 3090
      END IF

C///  GET THE NUMBER OF THE NEXT BLOCK TO CHANGE.

      CALL SQUEEZ (LENGTH, LIST)
      COUNT = 1
      DO 3060 J = 1, LEN (LIST)
         IF (LIST (J : J) .EQ. ' ') GO TO 3070
         COUNT = J
3060  CONTINUE
3070  CONTINUE

      IF (COUNT .LE. 8) THEN
         ITEM = ' '
         ITEM (9 - COUNT : 8) = LIST (1 : COUNT)
         LIST (1 : COUNT) = ' '
         READ
     +     (ERR = 3040,
     +      FMT = '(I8)',
     +      IOSTAT = STATUS,
     +      UNIT = ITEM)
     +      J
      ELSE
         LIST (1 : COUNT) = ' '
         GO TO 3040
      END IF

C///  CHANGE THE BLOCK.

      IF (1 .LE. J .AND. J .LE. NAMES) ACTIVE(J, 2) = .NOT. ACTIVE(J, 2)

C///  BOTTOM OF THE QUERY LOOP.

         GO TO 3040
3090  CONTINUE

C///  DETERMINE THE ACTION TO BE TAKEN ON EACH BLOCK.

      DO 3100 J = 1, NAMES
         IF (ACTIVE(J, 1) .EQV. ACTIVE(J, 2)) THEN
            CMND(J) = QCOPY
         ELSE IF (.NOT. ACTIVE(J, 1) .AND. ACTIVE(J, 2)) THEN
            CMND(J) = QSTRIP
         ELSE IF (ACTIVE(J, 1) .AND. .NOT. ACTIVE(J, 2)) THEN
            CMND(J) = QADD
         END IF
3100  CONTINUE

C///  DECIDE WHETHER TO WRITE A COPY.

      FLAG = .FALSE.
      DO 3110 J = 1, NAMES
         FLAG = FLAG .OR. (CMND(J) .NE. QCOPY)
3110  CONTINUE

      IF (.NOT. FLAG) THEN
         WRITE (*, 10011)
         GO TO 99999
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (4) OPEN THE COPY FILE.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE BLOCK.

4010  CONTINUE

C///  PROMPT FOR THE FILE NAME.

      CALL SQUEEZ (LENGTH, FILE1)
      WRITE (*, 10003) FILE1 (1 : LENGTH)
      READ
     +  (END = 99999,
     +   ERR = 9101,
     +   FMT = '(A)',
     +   IOSTAT = STATUS,
     +   UNIT = *)
     +   FILE2

      CALL SQUEEZ (LENGTH, FILE2)
      IF (FILE2 .EQ. ' ') GO TO 99999

C///  USE THE ORIGINAL FILE.

      IF (FILE1 .EQ. FILE2) THEN
         IF (.NOT. MAY) GO TO 9401

C///  OPEN THE SCRATCH FILE.

         ERROR = .TRUE.
         OPEN
     +     (ACCESS = 'SEQUENTIAL',
     +      ERR = 9402,
     +      FORM = 'FORMATTED',
     +      IOSTAT = STATUS,
     +      STATUS = 'SCRATCH',
     +      UNIT = UNIT2)
         ERROR = .FALSE.

C///  COPY FROM THE ORIGINAL FILE TO THE SCRATCH FILE.

         NUMBER = 0
         REWIND UNIT1
4020     CONTINUE
            READ
     +        (END = 4050,
     +         ERR = 9201,
     +         FMT = '(A)',
     +         IOSTAT = STATUS,
     +         UNIT = UNIT1)
     +         LINE
            NUMBER = NUMBER + 1

            LENGTH = 1
            DO 4030 J = MAXS, 1, - 1
               IF (LINE (J : J) .NE. ' ') THEN
                  LENGTH = J
                  GO TO 4040
               END IF
4030        CONTINUE
4040        CONTINUE

            WRITE (UNIT2, '(A)') LINE (1 : LENGTH)
            GO TO 4020
4050     CONTINUE

C///  SWITCH THE UNIT NUMBERS.

         ITEMP = UNIT1
         UNIT1 = UNIT2
         UNIT2 = ITEMP

C///  REWIND THE FILES.

         REWIND UNIT1
         REWIND UNIT2

C///  USE A NEW FILE.

      ELSE

C///  CHECK THE FILE.

         ERROR = .TRUE.
         INQUIRE
     +     (ERR = 9102,
     +      EXIST = EXIST,
     +      FILE = FILE2,
     +      IOSTAT = STATUS,
     +      OPENED = OPEN)
         ERROR = .FALSE.

         ERROR = EXIST .AND. OPEN
         IF (ERROR) GO TO 9106

C///  OPEN THE FILE.

         IF (EXIST) THEN
            ERROR = .TRUE.
            OPEN
     +        (ACCESS = 'SEQUENTIAL',
     +         ERR = 9105,
     +         FILE = FILE2,
     +         FORM = 'FORMATTED',
     +         IOSTAT = STATUS,
     +         STATUS = 'OLD',
     +         UNIT = UNIT2)
            ERROR = .FALSE.

C///  CREATE THE FILE.

         ELSE
            ERROR = .TRUE.
            OPEN
     +        (ACCESS = 'SEQUENTIAL',
     +         ERR = 9107,
     +         FILE = FILE2,
     +         FORM = 'FORMATTED',
     +         IOSTAT = STATUS,
     +         STATUS = 'NEW',
     +         UNIT = UNIT2)
            ERROR = .FALSE.
         END IF

C///  REWIND THE FILES.

         REWIND UNIT1
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (5) MAKE THE COPY.
C
C///////////////////////////////////////////////////////////////////////

      NUMBER = 0

C///  TOP OF THE LOOP.  READ THE NEXT LINE.

5010  CONTINUE
      READ
     +  (END = 5060,
     +   ERR = 9201,
     +   FMT = '(A)',
     +   IOSTAT = STATUS,
     +   UNIT = UNIT1)
     +   LINE
      NUMBER = NUMBER + 1

C///  GATHER INFORMATION ABOUT THE LINE.

      LENGTH = 0
      DO 5020 J = MAXS, 1, - 1
         IF (LINE (J : J) .NE. ' ') THEN
            LENGTH = J
            GO TO 5030
         END IF
5020  CONTINUE
5030  CONTINUE

      SPECL = LINE (1 : 6) .EQ. 'C*****'
     +   .OR. LINE (1 : 6) .EQ. 'c*****'
     +   .OR. LINE (1 : 6) .EQ. '******'
     +   .OR. LINE (1 : 6) .EQ. '#*****'

C///  TREAT SPECIAL LINES FOUND FROM OUTSIDE.

      IF (SPECL) THEN
         IF (LINE (1 : 1) .EQ. '#') THEN
            TYPE(0) = QSHELL
         ELSE
            TYPE(0) = QFORT
         END IF

         IF (.NOT. INSIDE) THEN
            FOUND = .FALSE.
            DO 5040 J = 1, NAMES
               IF (TYPE(0) .EQ. TYPE(J) .AND.
     +            LINE (7 : ) .EQ. NAME(J)) THEN
                  FOUND = .TRUE.
                  INDEX = J
               END IF
5040        CONTINUE
            ERROR = .NOT. FOUND
            IF (ERROR) GO TO 9501

            INSIDE = .TRUE.
            MARK = LINE (1 : 1)
            QCMND = QCOPY

C///  TREAT SPECIAL LINES FOUND FROM INSIDE.

         ELSE
            INDEX = 0
            INSIDE = .FALSE.
            QCMND = QCOPY
         END IF

C///  TREAT SIMPLE LINES FROM OUTSIDE.

      ELSE
         IF (.NOT. INSIDE) THEN
            QCMND = QCOPY

C///  TREAT SIMPLE LINES FROM INSIDE.

         ELSE
            LINES = LINES + 1
            QCMND = CMND(INDEX)
         END IF
      END IF

C///  WRITE THE LINE.

5050  CONTINUE
      IF (QCMND .EQ. QADD) THEN
         IF (0 .LT. LENGTH) THEN
            TEMP = MARK (1 : 1) // LINE (1 : LENGTH)
         ELSE
            TEMP = MARK (1 : 1)
         END IF
         LINE = TEMP
         LENGTH = LENGTH + 1
      ELSE IF (QCMND .EQ. QSTRIP) THEN
         LENGTH = LENGTH - 1
         TEMP = LINE (2 : )
         LINE = TEMP
      END IF

      IF (1 .LE. LENGTH) THEN
         WRITE (UNIT2, '(A)') LINE (1 : LENGTH)
      ELSE
         WRITE (UNIT2, '()')
      END IF

C///  BOTTOM OF THE LOOP THROUGH THE FILE.

      IF (MOD (NUMBER, REPORT) .EQ. 0) THEN
         IF (NUMBER .EQ. REPORT) WRITE (*, 10004)
         WRITE (TEMP, 10005) NUMBER
         CALL SQUEEZ (LENGTH, TEMP)
         WRITE (*, 10006) TEMP (1 : LENGTH)
      END IF

      GO TO 5010
5060  CONTINUE

      IF (MOD (NUMBER, REPORT) .NE. 0) THEN
         IF (NUMBER .LT. REPORT) WRITE (*, 10004)
         WRITE (TEMP, 10005) NUMBER
         CALL SQUEEZ (LENGTH, TEMP)
         WRITE (*, 10006) TEMP (1 : LENGTH)
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

10001 FORMAT
     +   (/ '    CHANGE:  VERSION 2.05 OF MARCH 1995',
     +   ' BY JOSEPH GRCAR.')

C*****$ EDIT DESCRIPTOR > YES (RECOMMENDED)
10002 FORMAT
     +  (/13X, 'ENTER THE NAME OF THE ORIGINAL FILE.'
     +   //1X, '     FILE?  ', $)

10003 FORMAT
     +  (/13X, 'ENTER THE NAME FOR THE COPY,'
     +   /13X, 'OR ENTER BLANK TO QUIT.'
     +  //13X, A, ' (NAME OF THE ORIGINAL FILE)'
     +   //1X, '     FILE?  ', $)
C*****END $ EDIT DESCRIPTOR > YES (RECOMMENDED)

C*****$ EDIT DESCRIPTOR > NO  (FORTRAN 77 STANDARD)
C10002 FORMAT
C     +  (/13X, 'ENTER THE NAME OF THE ORIGINAL FILE.'
C     +   /)
C
C10003 FORMAT
C     +  (/13X, 'ENTER THE NAME FOR THE COPY,'
C     +   /13X, 'OR ENTER BLANK TO QUIT.'
C     +   //1X, A, ' (NAME OF THE ORIGINAL FILE)'
C     +   /)
C*****END $ EDIT DESCRIPTOR > NO  (FORTRAN 77 STANDARD)

10004 FORMAT
     +   ()

10005 FORMAT
     +   (I10, ' LINES')

10006 FORMAT
     +   (13X, A)

10007 FORMAT
     +  (/13X, 'STATUS OF THE CHANGE BLOCKS.')

10008 FORMAT
     +  (/12X, 'ORIGINAL  COPY   NAME OF THE CHANGE BLOCK'
     +   /12X, ' ------  ------  --------------------------'
     +   /(1X, I9, A1, '  ', A6, 2X, A6, 2X, A))

C*****$ EDIT DESCRIPTOR > YES (RECOMMENDED)
10009 FORMAT
     +  (/13X, 'ENTER THE NUMBERS OF BLOCKS TO CHANGE,'
     +   /13X, 'OR ENTER BLANK TO ACCEPT THINGS AS THEY ARE.'
     +   //1X, '   NUMBER?  ', $)
C*****END $ EDIT DESCRIPTOR > YES (RECOMMENDED)

C*****$ EDIT DESCRIPTOR > NO  (FORTRAN 77 STANDARD)
C10009 FORMAT
C     +  (/13X, 'ENTER THE NUMBERS OF BLOCKS TO CHANGE,'
C     +   /13X, 'OR ENTER BLANK TO ACCEPT THINGS AS THEY ARE.'
C     +   /)
C*****END $ EDIT DESCRIPTOR > NO  (FORTRAN 77 STANDARD)

10010 FORMAT
     +  (/13X, 'NO CHANGE BLOCKS.')

10011 FORMAT
     +  (/13X, 'NO CHANGES.')

C*****PAUSE AT FINISH > NO  (RECOMMENDED)
10012 FORMAT
     +  (/13X, 'FINISH.'
     +   /)
C*****END PAUSE AT FINISH > NO  (RECOMMENDED)

C*****PAUSE AT FINISH > YES (SOME WINDOW ENVIRONMENTS)
C10012 FORMAT
C     +  (/13X, 'PRESS RETURN TO FINISH.')
C*****END PAUSE AT FINISH > YES (SOME WINDOW ENVIRONMENTS)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

9101  WRITE (TEMP, '(I20)') STATUS
      CALL SQUEEZ (LENGTH, TEMP)
      WRITE (*, 99101) TEMP (1 : LENGTH)
      IF (.NOT. MANUAL) GO TO 99999

9102  WRITE (TEMP, '(I20)') STATUS
      CALL SQUEEZ (LENGTH, TEMP)
      WRITE (*, 99102) TEMP (1 : LENGTH)
      IF (.NOT. MANUAL) GO TO 99999

9103  WRITE (*, 99103)
      IF (.NOT. MANUAL) GO TO 1010

9104  WRITE (*, 99104)
      IF (.NOT. MANUAL) GO TO 1010

9105  WRITE (TEMP, '(I20)') STATUS
      CALL SQUEEZ (LENGTH, TEMP)
      WRITE (*, 99105) TEMP (1 : LENGTH)
      IF (.NOT. MANUAL) GO TO 99999

9106  WRITE (*, 99106)
      IF (.NOT. MANUAL) GO TO 4010

9107  WRITE (TEMP, '(I20)') STATUS
      CALL SQUEEZ (LENGTH, TEMP)
      WRITE (*, 99107) TEMP (1 : LENGTH)
      IF (.NOT. MANUAL) GO TO 99999

9201  IF (NUMBER .EQ. 1) THEN
         WRITE (TEMP, '(I20, A)') NUMBER, ' LINE'
      ELSE
         WRITE (TEMP, '(I20, A)') NUMBER, ' LINES'
      END IF
      CALL SQUEEZ (LEN1, TEMP)
      WRITE (TEMP2, '(I20)') STATUS
      CALL SQUEEZ (LEN2, TEMP2)
      WRITE (*, 99201) TEMP (1 : LEN1), TEMP2 (1 : LEN2)
      IF (.NOT. MANUAL) GO TO 99999

9202  WRITE (TEMP, '(I20)') NUMBER
      WRITE (TEMP2, '(I20)') MAXL
      CALL SQUEEZ (LEN1, TEMP)
      CALL SQUEEZ (LEN2, TEMP2)
      WRITE (*, 99202) TEMP (1 : LEN1), TEMP2 (1 : LEN2)
      IF (.NOT. MANUAL) GO TO 99999

9203  WRITE (TEMP, '(I20)') MAXN
      CALL SQUEEZ (LENGTH, TEMP)
      WRITE (*, 99203) TEMP (1 : LENGTH)
      IF (.NOT. MANUAL) GO TO 99999

9204  CALL EXTENT (LEN1, TOP)
      CALL EXTENT (LEN2, LINE)
      WRITE (*, 99204)
     +   TOPNUM, TOP (1 : LEN1), NUMBER, LINE (1 : LEN2)
      IF (.NOT. MANUAL) GO TO 99999

9205  CALL EXTENT (LEN1, TOP)
      CALL EXTENT (LEN2, LINE)
      WRITE (*, 99205)
     +   TOPNUM, TOP (1 : LEN1), NUMBER, LINE (1 : LEN2)
      IF (.NOT. MANUAL) GO TO 99999

9206  CALL EXTENT (LEN1, TOP)
      CALL EXTENT (LEN2, LINE)
      WRITE (*, 99206)
     +   TOPNUM, TOP (1 : LEN1), NUMBER, LINE (1 : LEN2)
      IF (.NOT. MANUAL) GO TO 99999

9207  CALL EXTENT (LEN1, TOP)
      CALL EXTENT (LEN2, LINE)
      WRITE (*, 99207)
     +   TOPNUM, TOP (1 : LEN1), NUMBER, LINE (1 : LEN2)
      IF (.NOT. MANUAL) GO TO 99999

9208  CALL SQUEEZ (LEN1, TOP)
      IF (ACTIVE(INDEX, 1)) THEN
         WRITE (*, 99208) 'ACTIVE', 'INACTIVE', TOPNUM, TOP (1 : LEN1)
      ELSE
         WRITE (*, 99208) 'INACTIVE', 'ACTIVE', TOPNUM, TOP (1 : LEN1)
      END IF
      IF (.NOT. MANUAL) GO TO 99999

9209  WRITE (TEMP, '(I20)') SAVNUM
      WRITE (TEMP2, '(I20)') MAXL
      CALL SQUEEZ (LEN1, TEMP)
      CALL SQUEEZ (LEN2, TEMP2)
      WRITE (*, 99209) TEMP (1 : LEN1), TEMP2 (1 : LEN2)
      IF (.NOT. MANUAL) GO TO 99999

9210  WRITE (TEMP, '(I20)') TOPNUM
      CALL SQUEEZ (LEN1, TEMP)
      WRITE (TEMP2, '(I20)') SAVNUM
      CALL SQUEEZ (LEN2, TEMP2)
      WRITE (TEMP3, '(I20)') MAXL + 1
      CALL SQUEEZ (LEN3, TEMP3)
      CALL EXTENT (LEN4, TOP)
      SAVLIN (58 : MAXS) = ' '
      CALL EXTENT (LEN5, SAVLIN)
      WRITE (*, 99210)
     +   TEMP (1 : LEN1), TEMP2 (1 : LEN2), TEMP3 (1 : LEN3),
     +   TOPNUM, TOP (1 : LEN4), SAVNUM, SAVLIN (1 : LEN5) // '...'
      IF (.NOT. MANUAL) GO TO 99999

9211  WRITE (TEMP, '(I20)') TOPNUM
      CALL SQUEEZ (LEN1, TEMP)
      WRITE (TEMP2, '(I20)') SAVNUM
      CALL SQUEEZ (LEN2, TEMP2)
      WRITE (TEMP3, '(I20)') MAXL - 1
      CALL SQUEEZ (LEN3, TEMP3)
      CALL EXTENT (LEN4, TOP)
      SAVLIN (58 : MAXS) = ' '
      CALL EXTENT (LEN5, SAVLIN)
      WRITE (*, 99211)
     +   TEMP (1 : LEN1), TEMP2 (1 : LEN2), TEMP3 (1 : LEN3),
     +   TOPNUM, TOP (1 : LEN4), SAVNUM, SAVLIN (1 : LEN5) // '...'
      IF (.NOT. MANUAL) GO TO 99999

9212  WRITE (*, 99212)
      IF (.NOT. MANUAL) GO TO 99999

9401  WRITE (*, 99401)
      IF (.NOT. MANUAL) GO TO 4010

9402  WRITE (TEMP, '(I20)') STATUS
      CALL SQUEEZ (LENGTH, TEMP)
      WRITE (*, 99402) TEMP (1 : LENGTH)
      IF (.NOT. MANUAL) GO TO 99999

9501  CALL EXTENT (LEN1, LINE)
      WRITE (*, 99501) NUMBER, LINE (1 : LEN1)
      QCMND = QCOPY
      IF (.NOT. MANUAL) GO TO 5050

99101 FORMAT
     +  (/13X, 'ERROR.  READ FAILS WITH I/O STATUS ', A, '.')

99102 FORMAT
     +  (/13X, 'ERROR.  INQUIRE FAILS WITH I/O STATUS ', A, '.')

99103 FORMAT
     +  (/13X, 'ERROR.  THE FILE DOES NOT EXIST.')

99104 FORMAT
     +  (/13X, 'ERROR.  THE FILE IS OPEN TO SOMEONE ELSE.')

99105 FORMAT
     +  (/13X, 'ERROR.  OPEN FAILS WITH I/O STATUS ', A, '.')

99106 FORMAT
     +  (/13X, 'ERROR.  THE FILE IS OPEN TO SOMEONE ELSE.')

99107 FORMAT
     +  (/13X, 'ERROR.  FILE CREATE FAILS WITH I/O STATUS ', A, '.')

99201 FORMAT
     +  (/13X, 'ERROR.  AFTER ', A, ' SUCCESSFULLY READ,'
     +   /21X, 'A FILE READ FAILS WITH I/O STATUS ', A, '.')

99202 FORMAT
     +  (/13X, 'ERROR.  LINE NUMBER ', A, ' HAS OVER ', A,
     +   ' CHARACTERS.')

99203 FORMAT
     +  (/13X, 'ERROR.  THERE ARE OVER ', A, ' CHANGE BLOCK NAMES.')

99204 FORMAT
     +  (/13X, 'ERROR.  THE LINE AT THE BOTTOM OF A CHANGE BLOCK'
     +   /21X, 'DOES NOT SAY "END".'
     +   //1X, I9, ':  ', A
     +    /1X, I9, ':  ', A)

99205 FORMAT
     +  (/13X, 'ERROR.  THE LINES AT THE TOP AND BOTTOM OF A CHANGE'
     +   /21X, 'BLOCK HAVE DIFFERENT NAMES.'
     +   //1X, I9, ':  ', A
     +    /1X, I9, ':  ', A)

99206 FORMAT
     +  (/13X, 'ERROR.  THE LINES AT THE TOP AND BOTTOM OF A CHANGE'
     +   /21X, 'BLOCK HAVE DIFFERENT COMMENT MARKS.'
     +   //1X, I9, ':  ', A
     +    /1X, I9, ':  ', A)

99207 FORMAT
     +  (/13X, 'ERROR.  A CHANGE BLOCK HAS NO LINES BETWEEN THE TOP AND'
     +   /21X, 'BOTTOM LINE.'
     +   //1X, I9, ':  ', A
     +    /1X, I9, ':  ', A)

99208 FORMAT
     +  (/13X, 'ERROR.  THE CHANGE BLOCKS ARE INCONSISTENT.  AMONG'
     +   /21X, 'BLOCKS WITH THE NAME BELOW, THE FIRST IS ', A, ','
     +   /21X, 'BUT THE BLOCK BEGINNING ON THE LINE BELOW IS ', A, '.'
     +   //1X, I9, ':  ', A)

99209 FORMAT
     +  (/13X, 'ERROR.  LINE NUMBER ', A, ' HAS OVER ', A,
     +   ' CHARACTERS.')

99210 FORMAT
     +  (/13X, 'ERROR.  IN THE INACTIVE CHANGE BLOCK WHICH BEGINS ON ',
     +         'LINE ', A, ','
     +   /21X, 'LINE NUMBER ', A, ' IS LONGER THAN ', A, ' CHARACTERS.'
     +   //1X, I9, ':  ', A
     +    /1X, I9, ':  ', A)

99211 FORMAT
     +  (/13X, 'ERROR.  IN THE ACTIVE CHANGE BLOCK WHICH BEGINS ON ',
     +         'LINE ', A, ','
     +   /21X, 'LINE NUMBER ', A, ' IS LONGER THAN ', A, ' CHARACTERS.'
     +   //1X, I9, ':  ', A
     +    /1X, I9, ':  ', A)

99212 FORMAT
     +  (/13X, 'ERROR.  THE FILE ENDS IN THE MIDDLE OF A CHANGE',
     +   ' BLOCK.')

99401 FORMAT
     +  (/13X, 'ERROR.  THE COPY MAY NOT OVERWRITE THE ORIGINAL FILE.')

99402 FORMAT
     +  (/13X, 'ERROR.  SCRATCH OPEN FAILS WITH I/O STATUS ', A, '.')

99501 FORMAT
     +  (/13X, 'ERROR.  A CHANGE BLOCK NAME IS NOT RECOGNIZED.'
     +   //1X, I9, ':  ', A)

C///  EXIT.

99999 CONTINUE
      WRITE (*, 10012)
C*****PAUSE AT FINISH > YES (SOME WINDOW ENVIRONMENTS)
C      PAUSE ' '
C*****END PAUSE AT FINISH > YES (SOME WINDOW ENVIRONMENTS)
      END

      SUBROUTINE EXTENT (LENGTH, STRING)

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

      CHARACTER
     +   STRING*(*)
      INTEGER
     +   J, LENGTH
      INTRINSIC
     +   LEN

C     FORCE STATIC ALLOCATION OF VARIABLES.  THIS IS NOT NEEDED BY THE
C     CHANGE PROGRAM, BUT IT SEEMS TO HELP WITH TROUBLESOME COMPILERS.

      SAVE

      LENGTH = 1
      DO 0100 J = 1, LEN (STRING)
         IF (STRING(J : J) .NE. ' ') LENGTH = J
0100  CONTINUE

      RETURN
      END

      SUBROUTINE SQUEEZ (LENGTH, STRING)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     SQUEEZ
C
C     SQUEEZE LEADING BLANKS AND MULTIPLE BLANKS FROM A CHARACTER
C     STRING.  RETURN THE LENGTH OF THE STRING (OR 1 IF ALL BLANK).
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)

      CHARACTER
     +   CHAR*1, STRING*(*)
      INTEGER
     +   END, J, LENGTH
      INTRINSIC
     +   LEN
      LOGICAL
     +   BLANK

C     FORCE STATIC ALLOCATION OF VARIABLES.  THIS IS NOT NEEDED BY THE
C     CHANGE PROGRAM, BUT IT SEEMS TO HELP WITH TROUBLESOME COMPILERS.

      SAVE

C///  SQUEEZE THE STRING.

      LENGTH = 0
      BLANK = .TRUE.
      DO 0100 J = 1, LEN (STRING)
         CHAR = STRING (J : J)
         IF (.NOT. BLANK .OR. CHAR .NE. ' ') THEN
            BLANK = CHAR .EQ. ' '
            LENGTH = LENGTH + 1
            STRING (LENGTH : LENGTH) = CHAR
         END IF
0100  CONTINUE

C///  ADJUST THE LENGTH AND PAD THE STRING.

      IF (0 .LT. LENGTH) THEN
         IF (STRING (LENGTH : LENGTH) .EQ. ' ') LENGTH = LENGTH - 1
C        THE ABSOFT FORTRAN 77 COMPILER VERSION 3.1.2 FOR THE MACINTOSH
C        MIS-COMPILES THE FOLLOWING STATEMENT.
C        IF (LENGTH .LT. LEN (STRING)) STRING (LENGTH + 1 : ) = ' '
         END = LEN (STRING)
         IF (LENGTH .LT. LEN (STRING)) STRING (LENGTH + 1 : END) = ' '
      ELSE
         LENGTH = 1
      END IF

C///  EXIT.

      RETURN
      END

