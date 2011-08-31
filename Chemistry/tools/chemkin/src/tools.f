C     CVS $Revision: 1.1.1.1 $ reposited $Date: 2006/05/26 19:09:33 $

* * * * * * * * * * * * * * * * *
* Copyright (C) 1997            *
* Joseph F. Grcar               *
* Sandia National Laboratories  *
* sepp@california.sandia.gov    *
* * * * * * * * * * * * * * * * *

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     VERSION 1.04 OF DECEMBER 1997
C
C     FORTRAN TOOLS FOR CHARACTER MANIPULATION, KEYWORD INPUT, OUTPUT OF
C     PLOTTING DATA AND WORKSPACE ACCOUNTING.
C
C     WRITTEN BY:  DR. JOSEPH F. GRCAR
C                  DIVISION 8345
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
C     SUBROUTINES FOR CHARACTER MANIPULATION:
C
C     EXTENT - EXTENT (LENGTH) OF A STRING
C     RIGHTJ - RIGHT JUSTIFY A STRING
C     SQUEEZ - SQUEEZ LEADING AND MULTIPLE BLANKS FROM A STRING
C
C     SUBROUTINES FOR KEYWORD INPUT:
C
C     READI  - READ AN INTEGER NUMBER
C     READL  - READ A LOGICAL VALUE
C     READR1 - READ A REAL NUMBER (SINGLE PRECISION)
C     READR2 - READ A REAL NUMBER (DOUBLE PRECISION)
C     READV1 - READ VALUES FOR A MENU OF KEYWORDS (SINGLE PRECISION)
C     READV2 - READ VALUES FOR A MENU OF KEYWORDS (DOUBLE PRECISION)
C     READW  - READ A WORD (A BLANK OR DOUBLE QUOTED DELIMITED STRING)
C     VERIFY - VERIFY THE NEXT WORDS IN THE INPUT FILE
C
C     SUBROUTINES FOR OUTPUT OF PLOTTING DATA:
C
C     WRT1D1 - WRITE A PLOT1D RECORD FROM SINGLE PRECISION DATA
C     WRT1D2 - WRITE A PLOT1D RECORD FROM DOUBLE PRECISION DATA
C     WRT2D1 - WRITE A PLOT2D RECORD FROM SINGLE PRECISION DATA
C     WRT2D2 - WRITE A PLOT2D RECORD FROM DOUBLE PRECISION DATA
C
C     SUBROUTINES FOR WORKSPACE ACCOUTING:
C
C     RESERV - RESERVE SPACE IN AN ARRAY
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

      IMPLICIT COMPLEX (A - Z)
      CHARACTER STRING*(*)
      INTEGER J, LENGTH

      LENGTH = 1
      DO 0100 J = 1, LEN (STRING)
         IF (STRING(J : J) .NE. ' ') LENGTH = J
0100  CONTINUE

      RETURN
      END
      SUBROUTINE READI (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, VALUE)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     READI
C
C     READ THE NEXT SPACE-DELIMITED WORD OR DOUBLE-QUOTE-DELIMITED 
C     STRING FROM A SCRIPT FILE AND INTERPRET IT AS AN INTEGER.
C
C     USES SUBROUTINES EXTENT AND READW.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER ID*9, LINE*(*), WORD*80
      INTEGER LENGTH, NUMBER, PRINT, SCRIPT, STATUS, TEXT, VALUE
      LOGICAL ERROR

      PARAMETER (ID = ' READI:  ')

C///  READ A STRING.

C     SUBROUTINE READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)

      CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)
      IF (ERROR) GO TO 9001

C///  CONVERT THE STRING.

      READ (WORD, '(BN, I80)', ERR = 9002, IOSTAT = STATUS) VALUE

C///  ERROR MESSAGES.

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID
      GO TO 99999

9002  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, WORD)
         IF (58 .LT. LENGTH) THEN
            WORD (54 : ) = ' ...'
            LENGTH = 58
         END IF
         WRITE (TEXT, 99002) ID, SCRIPT, NUMBER, WORD (1 : LENGTH)
      END IF
      ERROR = .TRUE.
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  READW FAILS.')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  THE NEXT ITEM IN THE SCRIPT FILE SHOULD BE'
     +   /10X, 'AN INTEGER.'
     +  //10X, I10, '  UNIT NUMBER'
     +   /10X, I10, '  LINE NUMBER'
     +  //10X, '   STRING:  ', A)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE READL (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, VALUE)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     READL
C
C     READ THE NEXT SPACE-DELIMITED WORD OR DOUBLE-QUOTE-DELIMITED 
C     STRING FROM A SCRIPT FILE AND INTERPRET IT AS A LOGICAL.
C
C     USES SUBROUTINES EXTENT AND READW.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER ID*9, LINE*(*), WORD*80
      INTEGER LENGTH, NUMBER, PRINT, SCRIPT, TEXT
      LOGICAL ERROR, VALUE

      PARAMETER (ID = ' READL:  ')

C///  READ A STRING.

C     SUBROUTINE READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)

      CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)
      IF (ERROR) GO TO 9001

C///  CONVERT THE STRING.

      IF (WORD .EQ. '.T.' .OR.
     +   WORD .EQ. '.TRUE.' .OR. WORD .EQ. 'TRUE' .OR. 
     +   WORD .EQ. '.YES.' .OR. WORD .EQ. 'YES') THEN
         VALUE = .TRUE.
      ELSE IF (WORD .EQ. '.F.' .OR.
     +   WORD .EQ. '.FALSE.' .OR. WORD .EQ. 'FALSE' .OR. 
     +   WORD .EQ. '.NO.' .OR. WORD .EQ. 'NO') THEN
         VALUE = .FALSE.
      ELSE
         ERROR = .TRUE.
         GO TO 9002
      END IF

C///  ERROR MESSAGES.

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID
      GO TO 99999

9002  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, WORD)
         IF (58 .LT. LENGTH) THEN
            WORD (54 : ) = ' ...'
            LENGTH = 58
         END IF
         WRITE (TEXT, 99002) ID, SCRIPT, NUMBER, WORD (1 : LENGTH)
      END IF
      ERROR = .TRUE.
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  READW FAILS.')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  THE NEXT ITEM IN THE SCRIPT FILE SHOULD BE'
     +   /10X, 'A LOGICAL.'
     +  //10X, I10, '  UNIT NUMBER'
     +   /10X, I10, '  LINE NUMBER'
     +  //10X, '   STRING:  ', A)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE READR1 
     +   (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, VALUE)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     READR1 
C
C     READ THE NEXT SPACE-DELIMITED WORD OR DOUBLE-QUOTE-DELIMITED
C     STRING FROM A SCRIPT FILE AND INTERPRET IT AS A SINGLE PRECISION
C     REAL NUMBER.
C
C     USES SUBROUTINES EXTENT AND READW.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER ID*9, LINE*(*), WORD*80
      INTEGER LENGTH, NUMBER, PRINT, SCRIPT, STATUS, TEXT
      LOGICAL ERROR
      REAL VALUE

      PARAMETER (ID = 'READR1:  ')

C///  READ A STRING.

C     SUBROUTINE READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)

      CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)
      IF (ERROR) GO TO 9001

C///  CONVERT THE STRING.

      READ (WORD, '(BN, F80.0)', ERR = 9002, IOSTAT = STATUS) VALUE

C///  ERROR MESSAGES.

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID
      GO TO 99999

9002  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, WORD)
         IF (58 .LT. LENGTH) THEN
            WORD (54 : ) = ' ...'
            LENGTH = 58
         END IF
         WRITE (TEXT, 99002) ID, SCRIPT, NUMBER, WORD (1 : LENGTH)
      END IF
      ERROR = .TRUE.
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  READW FAILS.')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  THE NEXT ITEM IN THE SCRIPT FILE SHOULD BE'
     +   /10X, 'A REAL NUMBER.'
     +  //10X, I10, '  UNIT NUMBER'
     +   /10X, I10, '  LINE NUMBER'
     +  //10X, '   STRING:  ', A)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE READR2
     +   (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, VALUE)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     READR2
C
C     READ THE NEXT SPACE-DELIMITED WORD OR DOUBLE-QUOTE-DELIMITED
C     STRING FROM A SCRIPT FILE AND INTERPRET IT AS A DOUBLE PRECISION
C     REAL NUMBER.
C
C     USES SUBROUTINES EXTENT AND READW.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER ID*9, LINE*(*), WORD*80
      INTEGER LENGTH, NUMBER, PRINT, SCRIPT, STATUS, TEXT
      LOGICAL ERROR
      DOUBLE PRECISION VALUE

      PARAMETER (ID = 'READR2:  ')

C///  READ A STRING.

C     SUBROUTINE READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)

      CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)
      IF (ERROR) GO TO 9001

C///  CONVERT THE STRING.

      READ (WORD, '(BN, F80.0)', ERR = 9002, IOSTAT = STATUS) VALUE

C///  ERROR MESSAGES.

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID
      GO TO 99999

9002  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, WORD)
         IF (58 .LT. LENGTH) THEN
            WORD (54 : ) = ' ...'
            LENGTH = 58
         END IF
         WRITE (TEXT, 99002) ID, SCRIPT, NUMBER, WORD (1 : LENGTH)
      END IF
      ERROR = .TRUE.
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  READW FAILS.')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  THE NEXT ITEM IN THE SCRIPT FILE SHOULD BE'
     +   /10X, 'A REAL NUMBER.'
     +  //10X, I10, '  UNIT NUMBER'
     +   /10X, I10, '  LINE NUMBER'
     +  //10X, '   STRING:  ', A)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE READV1
     +  (ERROR, TEXT,
     +   CVALUE, ESIGN, FOUND, KEY, KEYS, IVALUE, LINE, LVALUE, NEED,
     +   NUMBER, PRINT, RVALUE, SCRIPT, TYPE)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     READV1 
C
C     READ VALUES FOR KEYWORD PHRASES FROM A SCRIPT FILE, RETURNING
C     SINGLE PRECISION VALUES FOR REAL NUMBERS.
C
C     USES SUBROUTINES EXTENT, READI, READR1 AND READW.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   CVALUE*(*), ID*9, KEY*(*), LINE*(*), PHRASE*132, STRING*132,
     +   TYPE*(*)
      INTEGER
     +   COUNT1, COUNT2, INDEX, IVALUE, J, KEYS, LENGTH, NUMBER, PINDEX,
     +   PRINT, SCRIPT, TEXT
      LOGICAL ERROR, ESIGN, FOUND, LVALUE, NEED
      REAL RVALUE

      PARAMETER (ID = 'READV1:  ')

      DIMENSION
     +   CVALUE(KEYS), FOUND(KEYS), IVALUE(KEYS), KEY(KEYS),
     +   LVALUE(KEYS), NEED(KEYS), RVALUE(KEYS), TYPE(KEYS)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (0 .LT. KEYS)
      IF (ERROR) GO TO 9001

      DO 1100 J = 1, KEYS
         CALL EXTENT (LENGTH, KEY(J))
         ERROR = .NOT. (LENGTH .LE. LEN (PHRASE))
         IF (ERROR) GO TO 9002

         ERROR = KEY(J) .EQ. 'END'
         IF (ERROR) GO TO 9003
1100  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     READ THE PHRASES.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP OVER THE PHRASES.

      DO 1200 J = 1, KEYS
         FOUND(J) = .FALSE.
1200  CONTINUE

1300  CONTINUE

      PHRASE = ' '
      PINDEX = 0

C///  TOP OF THE LOOP OVER THE STRINGS.

1400  CONTINUE

C///  APPEND A STRING TO THE PHRASE.

C     SUBROUTINE READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)

      CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, STRING)
      IF (ERROR) GO TO 9004

      CALL EXTENT (LENGTH, STRING)

      IF (PINDEX .EQ. 0) THEN
         PHRASE = STRING
         PINDEX = LENGTH
      ELSE IF (PINDEX + 2 .LE. LEN (PHRASE)) THEN
         PHRASE (PINDEX + 2 : ) = STRING
         PINDEX = MIN (LEN (PHRASE), PINDEX + 1 + LENGTH)
      END IF

C///  COMPARE THE PHRASE WITH THE KEYWORDS.

      IF (PHRASE .EQ. 'END') GO TO 1600

      COUNT1 = 0
      COUNT2 = 0
C     ADDING 1 INCLUDES THE NEXT CHARACTER IN THE KEY AND PREVENTS
C     MATCHES TO PARTIAL WORDS.
      LENGTH = MIN (PINDEX + 1, LEN (KEY(1)))
      DO 1500 J = 1, KEYS
         IF (PHRASE .EQ. KEY(J) (1 : LENGTH)) THEN
            COUNT1 = COUNT1 + 1
            IF (PHRASE .EQ. KEY(J)) THEN
               COUNT2 = COUNT2 + 1
               INDEX = J
            END IF
         END IF
1500  CONTINUE

      ERROR = 0 .EQ. COUNT1
      IF (ERROR) GO TO 9005

      ERROR = (1 .LT. COUNT1 .AND. 0 .LT. COUNT2) .OR. (1 .LT. COUNT2)
      IF (ERROR) GO TO 9006

C///  READ THE VALUE.

      IF (COUNT2 .EQ. 1) THEN
         ERROR = FOUND(INDEX)
         IF (ERROR) GO TO 9007

         IF (ESIGN .AND. TYPE(INDEX) .NE. 'N') THEN
            CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, '=')
            IF (ERROR) GO TO 9008
         END IF

         FOUND(INDEX) = .TRUE.

         IF (TYPE(INDEX) .EQ. 'C') THEN
            CALL READW
     +         (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, CVALUE(INDEX))
            IF (ERROR) GO TO 9009
         ELSE IF (TYPE(INDEX) .EQ. 'I') THEN
            CALL READI
     +         (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, IVALUE(INDEX))
            IF (ERROR) GO TO 9010
         ELSE IF (TYPE(INDEX) .EQ. 'L') THEN
            CALL READW
     +         (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, STRING)
            IF (ERROR) GO TO 9004

            IF (STRING .EQ. 'F' .OR. STRING .EQ. 'FALSE' .OR.
     +         STRING .EQ. 'NO') THEN
               LVALUE(INDEX) = .FALSE.
            ELSE IF (STRING .EQ. 'T' .OR. STRING .EQ. 'TRUE' .OR.
     +         STRING .EQ. 'YES') THEN
               LVALUE(INDEX) = .TRUE.
            ELSE
               ERROR = .TRUE.
               GO TO 9011
            END IF
         ELSE IF (TYPE(INDEX) .EQ. 'N') THEN
            CONTINUE
         ELSE IF (TYPE(INDEX) .EQ. 'R') THEN
            CALL READR1
     +         (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, RVALUE(INDEX))
            IF (ERROR) GO TO 9012
         ELSE
            ERROR = .TRUE.
            GO TO 9013
         END IF

C///  BOTTOM OF THE LOOP OVER THE STRINGS.

      ELSE
         GO TO 1400
      END IF

C///  BOTTOM OF THE LOOP OVER THE PHRASES.

      GO TO 1300
1600  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE NEEDED KEYWORDS.

      DO 2100 J = 1, KEYS
         ERROR = (.NOT. FOUND(J)) .AND. NEED(J)
         IF (ERROR) GO TO 9014
2100  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, KEYS
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, LENGTH, LEN (PHRASE)
      GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID
      GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID
      GO TO 99999

9005  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, PHRASE)
         IF (58 .LT. LENGTH) PHRASE (54 : ) = ' ...'
         CALL EXTENT (LENGTH, PHRASE)
         WRITE (TEXT, 99005) ID, NUMBER, PHRASE (1 : LENGTH)
         DO 3100 J = 1, KEYS
            IF (.NOT. FOUND(J)) THEN
               STRING = KEY(J)
               CALL EXTENT (LENGTH, STRING)
               IF (58 .LT. LENGTH) STRING (54 : ) = ' ...'
               CALL EXTENT (LENGTH, STRING)
               WRITE (TEXT, '(10X, A, A)')
     +            ' EXPECTED:  ', STRING (1 : LENGTH)
            END IF
3100     CONTINUE
      END IF
      GO TO 99999

9006  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99006) ID
         DO 3200 J = 1, KEYS
            STRING = KEY(J)
            CALL EXTENT (LENGTH, STRING)
            IF (58 .LT. LENGTH) STRING (54 : ) = ' ...'
            CALL EXTENT (LENGTH, STRING)
            WRITE (TEXT, '(10X, A, A)')
     +         '  KEYWORD:  ', STRING (1 : LENGTH)
3200     CONTINUE
      END IF
      GO TO 99999

9007  IF (0 .LT. TEXT) THEN
         STRING = KEY(INDEX)
         CALL EXTENT (LENGTH, STRING)
         IF (58 .LT. LENGTH) STRING (54 : ) = ' ...'
         CALL EXTENT (LENGTH, STRING)
         WRITE (TEXT, 99007) ID, NUMBER, STRING (1 : LENGTH)
      END IF
      GO TO 99999

9008  IF (0 .LT. TEXT) WRITE (TEXT, 99008) ID
      GO TO 99999

9009  IF (0 .LT. TEXT) WRITE (TEXT, 99009) ID
      GO TO 99999

9010  IF (0 .LT. TEXT) WRITE (TEXT, 99010) ID
      GO TO 99999

9011  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, STRING)
         IF (58 .LT. LENGTH) STRING (54 : ) = ' ...'
         CALL EXTENT (LENGTH, STRING)
         WRITE (TEXT, 99011) ID, NUMBER, STRING (1 : LENGTH)
      END IF
      GO TO 99999

9012  IF (0 .LT. TEXT) WRITE (TEXT, 99012) ID
      GO TO 99999

9013  IF (0 .LT. TEXT) WRITE (TEXT, 99013) ID
      GO TO 99999

9014  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99014) ID
         DO 3300 J = 1, KEYS
            IF ((.NOT. FOUND(J)) .AND. NEED(J)) THEN
               STRING = KEY(J)
               CALL EXTENT (LENGTH, STRING)
               IF (58 .LT. LENGTH) STRING (54 : ) = ' ...'
               CALL EXTENT (LENGTH, STRING)
               WRITE (TEXT, '(10X, A, A)')
     +            '  KEYWORD:  ', STRING (1 : LENGTH)
            END IF
3300     CONTINUE
      END IF
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  THE NUMBER OF KEYWORDS MUST BE POSITIVE.'
     +  //10X, I10, '  KEYS')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  A KEYWORD PHRASE IS TOO LONG.'
     +  //10X, I10, '  LENGTH'
     +   /10X, I10, '  LIMIT')

99003 FORMAT
     +   (/1X, A9, 'ERROR.  "END" CANNOT BE A KEYWORD.')

99004 FORMAT
     +   (/1X, A9, 'ERROR.  READW FAILS.')

99005 FORMAT
     +   (/1X, A9, 'ERROR.  THE NEXT PHRASE IN THE SCRIPT FILE IS NOT'
     +   /10X, 'THE ONE EXPECTED.'
     +  //10X, I10, '  LINE NUMBER'
     +  //10X, '    FOUND:  ', A
     +  //10X, ' EXPECTED:  END')

99006 FORMAT
     +   (/1X, A9, 'ERROR.  SOME KEYWORD PHRASES DUPLICATE OR ENCOMPASS'
     +   /10X, 'OTHERS.'
     +   /)

99007 FORMAT
     +   (/1X, A9, 'ERROR.  THE INPUT SCRIPT REPEATS A KEYWORD.'
     +  //10X, I10, '  LINE NUMBER'
     +  //10X, '  KEYWORD:  ', A)

99008 FORMAT
     +   (/1X, A9, 'ERROR.  VERIFY FAILS.')

99009 FORMAT
     +   (/1X, A9, 'ERROR.  READW FAILS.')

99010 FORMAT
     +   (/1X, A9, 'ERROR.  READI FAILS.')

99011 FORMAT
     +   (/1X, A9, 'ERROR.  A LOGICAL VALUE IS EXPECTED.'
     +  //10X, I10, '  LINE NUMBER'
     +  //10X, '    FOUND:  ', A
     +  //10X, ' EXPECTED:  F, FALSE, NO, T, TRUE, YES')

99012 FORMAT
     +   (/1X, A9, 'ERROR.  READR1 FAILS.')

99013 FORMAT
     +   (/1X, A9, 'ERROR.  THE TYPE OF THE VALUE TO BE READ IS NOT'
     +   /10X, 'RECOGNIZED.')

99014 FORMAT
     +   (/1X, A9, 'ERROR.  ONE OR MORE NEEDED KEYWORDS ARE MISSING.'
     +   /)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE READV2
     +  (ERROR, TEXT,
     +   CVALUE, ESIGN, FOUND, KEY, KEYS, IVALUE, LINE, LVALUE, NEED,
     +   NUMBER, PRINT, RVALUE, SCRIPT, TYPE)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     READV2 
C
C     READ VALUES FOR KEYWORD PHRASES FROM A SCRIPT FILE, RETURNING
C     DOUBLE PRECISION VALUES FOR REAL NUMBERS.
C
C     USES SUBROUTINES EXTENT, READI, READR1 AND READW.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   CVALUE*(*), ID*9, KEY*(*), LINE*(*), PHRASE*132, STRING*132,
     +   TYPE*(*)
      INTEGER
     +   COUNT1, COUNT2, INDEX, IVALUE, J, KEYS, LENGTH, NUMBER, PINDEX,
     +   PRINT, SCRIPT, TEXT
      LOGICAL ERROR, ESIGN, FOUND, LVALUE, NEED
      DOUBLE PRECISION RVALUE

      PARAMETER (ID = 'READV2:  ')

      DIMENSION
     +   CVALUE(KEYS), FOUND(KEYS), IVALUE(KEYS), KEY(KEYS),
     +   LVALUE(KEYS), NEED(KEYS), RVALUE(KEYS), TYPE(KEYS)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (0 .LT. KEYS)
      IF (ERROR) GO TO 9001

      DO 1100 J = 1, KEYS
         CALL EXTENT (LENGTH, KEY(J))
         ERROR = .NOT. (LENGTH .LE. LEN (PHRASE))
         IF (ERROR) GO TO 9002

         ERROR = KEY(J) .EQ. 'END'
         IF (ERROR) GO TO 9003
1100  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     READ THE PHRASES.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP OVER THE PHRASES.

      DO 1200 J = 1, KEYS
         FOUND(J) = .FALSE.
1200  CONTINUE

1300  CONTINUE

      PHRASE = ' '
      PINDEX = 0

C///  TOP OF THE LOOP OVER THE STRINGS.

1400  CONTINUE

C///  APPEND A STRING TO THE PHRASE.

C     SUBROUTINE READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)

      CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, STRING)
      IF (ERROR) GO TO 9004

      CALL EXTENT (LENGTH, STRING)

      IF (PINDEX .EQ. 0) THEN
         PHRASE = STRING
         PINDEX = LENGTH
      ELSE IF (PINDEX + 2 .LE. LEN (PHRASE)) THEN
         PHRASE (PINDEX + 2 : ) = STRING
         PINDEX = MIN (LEN (PHRASE), PINDEX + 1 + LENGTH)
      END IF

C///  COMPARE THE PHRASE WITH THE KEYWORDS.

      IF (PHRASE .EQ. 'END') GO TO 1600

      COUNT1 = 0
      COUNT2 = 0
C     ADDING 1 INCLUDES THE NEXT CHARACTER IN THE KEY AND PREVENTS
C     MATCHES TO PARTIAL WORDS.
      LENGTH = MIN (PINDEX + 1, LEN (KEY(1)))
      DO 1500 J = 1, KEYS
         IF (PHRASE .EQ. KEY(J) (1 : LENGTH)) THEN
            COUNT1 = COUNT1 + 1
            IF (PHRASE .EQ. KEY(J)) THEN
               COUNT2 = COUNT2 + 1
               INDEX = J
            END IF
         END IF
1500  CONTINUE

      ERROR = 0 .EQ. COUNT1
      IF (ERROR) GO TO 9005

      ERROR = (1 .LT. COUNT1 .AND. 0 .LT. COUNT2) .OR. (1 .LT. COUNT2)
      IF (ERROR) GO TO 9006

C///  READ THE VALUE.

      IF (COUNT2 .EQ. 1) THEN
         ERROR = FOUND(INDEX)
         IF (ERROR) GO TO 9007

         IF (ESIGN .AND. TYPE(INDEX) .NE. 'N') THEN
            CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, '=')
            IF (ERROR) GO TO 9008
         END IF

         FOUND(INDEX) = .TRUE.

         IF (TYPE(INDEX) .EQ. 'C') THEN
            CALL READW
     +         (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, CVALUE(INDEX))
            IF (ERROR) GO TO 9009
         ELSE IF (TYPE(INDEX) .EQ. 'I') THEN
            CALL READI
     +         (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, IVALUE(INDEX))
            IF (ERROR) GO TO 9010
         ELSE IF (TYPE(INDEX) .EQ. 'L') THEN
            CALL READW
     +         (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, STRING)
            IF (ERROR) GO TO 9004

            IF (STRING .EQ. 'F' .OR. STRING .EQ. 'FALSE' .OR.
     +         STRING .EQ. 'NO') THEN
               LVALUE(INDEX) = .FALSE.
            ELSE IF (STRING .EQ. 'T' .OR. STRING .EQ. 'TRUE' .OR.
     +         STRING .EQ. 'YES') THEN
               LVALUE(INDEX) = .TRUE.
            ELSE
               ERROR = .TRUE.
               GO TO 9011
            END IF
         ELSE IF (TYPE(INDEX) .EQ. 'N') THEN
            CONTINUE
         ELSE IF (TYPE(INDEX) .EQ. 'R') THEN
            CALL READR2
     +         (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, RVALUE(INDEX))
            IF (ERROR) GO TO 9012
         ELSE
            ERROR = .TRUE.
            GO TO 9013
         END IF

C///  BOTTOM OF THE LOOP OVER THE STRINGS.

      ELSE
         GO TO 1400
      END IF

C///  BOTTOM OF THE LOOP OVER THE PHRASES.

      GO TO 1300
1600  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE NEEDED KEYWORDS.

      DO 2100 J = 1, KEYS
         ERROR = (.NOT. FOUND(J)) .AND. NEED(J)
         IF (ERROR) GO TO 9014
2100  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, KEYS
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, LENGTH, LEN (PHRASE)
      GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID
      GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID
      GO TO 99999

9005  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, PHRASE)
         IF (58 .LT. LENGTH) PHRASE (54 : ) = ' ...'
         CALL EXTENT (LENGTH, PHRASE)
         WRITE (TEXT, 99005) ID, NUMBER, PHRASE (1 : LENGTH)
         DO 3100 J = 1, KEYS
            IF (.NOT. FOUND(J)) THEN
               STRING = KEY(J)
               CALL EXTENT (LENGTH, STRING)
               IF (58 .LT. LENGTH) STRING (54 : ) = ' ...'
               CALL EXTENT (LENGTH, STRING)
               WRITE (TEXT, '(10X, A, A)')
     +            ' EXPECTED:  ', STRING (1 : LENGTH)
            END IF
3100     CONTINUE
      END IF
      GO TO 99999

9006  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99006) ID
         DO 3200 J = 1, KEYS
            STRING = KEY(J)
            CALL EXTENT (LENGTH, STRING)
            IF (58 .LT. LENGTH) STRING (54 : ) = ' ...'
            CALL EXTENT (LENGTH, STRING)
            WRITE (TEXT, '(10X, A, A)')
     +         '  KEYWORD:  ', STRING (1 : LENGTH)
3200     CONTINUE
      END IF
      GO TO 99999

9007  IF (0 .LT. TEXT) THEN
         STRING = KEY(INDEX)
         CALL EXTENT (LENGTH, STRING)
         IF (58 .LT. LENGTH) STRING (54 : ) = ' ...'
         CALL EXTENT (LENGTH, STRING)
         WRITE (TEXT, 99007) ID, NUMBER, STRING (1 : LENGTH)
      END IF
      GO TO 99999

9008  IF (0 .LT. TEXT) WRITE (TEXT, 99008) ID
      GO TO 99999

9009  IF (0 .LT. TEXT) WRITE (TEXT, 99009) ID
      GO TO 99999

9010  IF (0 .LT. TEXT) WRITE (TEXT, 99010) ID
      GO TO 99999

9011  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, STRING)
         IF (58 .LT. LENGTH) STRING (54 : ) = ' ...'
         CALL EXTENT (LENGTH, STRING)
         WRITE (TEXT, 99011) ID, NUMBER, STRING (1 : LENGTH)
      END IF
      GO TO 99999

9012  IF (0 .LT. TEXT) WRITE (TEXT, 99012) ID
      GO TO 99999

9013  IF (0 .LT. TEXT) WRITE (TEXT, 99013) ID
      GO TO 99999

9014  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99014) ID
         DO 3300 J = 1, KEYS
            IF ((.NOT. FOUND(J)) .AND. NEED(J)) THEN
               STRING = KEY(J)
               CALL EXTENT (LENGTH, STRING)
               IF (58 .LT. LENGTH) STRING (54 : ) = ' ...'
               CALL EXTENT (LENGTH, STRING)
               WRITE (TEXT, '(10X, A, A)')
     +            '  KEYWORD:  ', STRING (1 : LENGTH)
            END IF
3300     CONTINUE
      END IF
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  THE NUMBER OF KEYWORDS MUST BE POSITIVE.'
     +  //10X, I10, '  KEYS')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  A KEYWORD PHRASE IS TOO LONG.'
     +  //10X, I10, '  LENGTH'
     +   /10X, I10, '  LIMIT')

99003 FORMAT
     +   (/1X, A9, 'ERROR.  "END" CANNOT BE A KEYWORD.')

99004 FORMAT
     +   (/1X, A9, 'ERROR.  READW FAILS.')

99005 FORMAT
     +   (/1X, A9, 'ERROR.  THE NEXT PHRASE IN THE SCRIPT FILE IS NOT'
     +   /10X, 'THE ONE EXPECTED.'
     +  //10X, I10, '  LINE NUMBER'
     +  //10X, '    FOUND:  ', A
     +  //10X, ' EXPECTED:  END')

99006 FORMAT
     +   (/1X, A9, 'ERROR.  SOME KEYWORD PHRASES DUPLICATE OR ENCOMPASS'
     +   /10X, 'OTHERS.'
     +   /)

99007 FORMAT
     +   (/1X, A9, 'ERROR.  THE INPUT SCRIPT REPEATS A KEYWORD.'
     +  //10X, I10, '  LINE NUMBER'
     +  //10X, '  KEYWORD:  ', A)

99008 FORMAT
     +   (/1X, A9, 'ERROR.  VERIFY FAILS.')

99009 FORMAT
     +   (/1X, A9, 'ERROR.  READW FAILS.')

99010 FORMAT
     +   (/1X, A9, 'ERROR.  READI FAILS.')

99011 FORMAT
     +   (/1X, A9, 'ERROR.  A LOGICAL VALUE IS EXPECTED.'
     +  //10X, I10, '  LINE NUMBER'
     +  //10X, '    FOUND:  ', A
     +  //10X, ' EXPECTED:  F, FALSE, NO, T, TRUE, YES')

99012 FORMAT
     +   (/1X, A9, 'ERROR.  READR2 FAILS.')

99013 FORMAT
     +   (/1X, A9, 'ERROR.  THE TYPE OF THE VALUE TO BE READ IS NOT'
     +   /10X, 'RECOGNIZED.')

99014 FORMAT
     +   (/1X, A9, 'ERROR.  ONE OR MORE NEEDED KEYWORDS ARE MISSING.'
     +   /)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     READW
C
C     READ THE NEXT SPACE-DELIMITED WORD OR DOUBLE-QUOTE-DELIMITED 
C     STRING FROM A SCRIPT FILE.
C
C     USES SUBROUTINE EXTENT.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER ID*9, LINE*(*), TEMP*132, WORD*(*)
      INTEGER 
     +   ACT1, ACT2, ACT3, ACT4, LENGTH, LIMIT, LINDEX, NUMBER, PRINT, 
     +   SCRIPT, STATE, STATUS, SYMBOL, TEXT, WINDEX
      LOGICAL DONE, ERROR, HOLD

      PARAMETER (ID = ' READW:  ')

      DIMENSION ACT1(4, 3), ACT2(4, 3), ACT3(4, 3), ACT4(4, 3)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (0 .LE. NUMBER)
      IF (ERROR) GO TO 9001

      ERROR = .NOT. (0 .LT. SCRIPT)
      IF (ERROR) GO TO 9002

C///  EVALUATE THE LINE LENGTH LIMIT.

      LIMIT = LEN (LINE)

C///////////////////////////////////////////////////////////////////////
C
C     READ A WORD OR STRING.
C
C///////////////////////////////////////////////////////////////////////

C///  INITIALIZE.

C     STATE: 1 SEARCHING, 2 BUILDING WORD, 3 BUILDING STRING

      DONE = .FALSE.
      HOLD = .FALSE.
      LINDEX = 1
      STATE = 1
      WINDEX = 0
      WORD = ' '

C///  TOP OF THE LOOP OVER THE CHARACTERS.

0100  CONTINUE

C///  READ A LINE.

      IF (0 .EQ. NUMBER .OR. LIMIT .LT. LINDEX) THEN
         READ (SCRIPT, '(A)', END = 9003, ERR = 9004, IOSTAT = STATUS) 
     +      TEMP
         NUMBER = NUMBER + 1

         CALL EXTENT (LENGTH, TEMP)
 
         IF (0 .LT. PRINT .AND. 0 .LT. TEXT) THEN
            IF (LENGTH .LE. 70) THEN
               WRITE (TEXT, '(10X, A)') TEMP (1 : LENGTH)
            ELSE
               WRITE (TEXT, '(10X, A, A)') TEMP (1 : 66), ' ...'
            END IF
         END IF

         ERROR = .NOT. (LENGTH + 2 .LE. LEN (LINE))
         IF (ERROR) GO TO 9005

         LINDEX = 1
         LINE = ' ' // TEMP (1 : LENGTH) // '!'
      END IF

C///  INTERPRET THE CHARACTER.

C     SYMBOL: 1 BLANK, 2 END OF LINE, 3 QUOTES, 4 OTHER

      IF (LINE (LINDEX : LINDEX) .EQ. ' ') THEN
         SYMBOL = 1
      ELSE IF (LINE (LINDEX : LINDEX) .EQ. '!') THEN
         SYMBOL = 2
      ELSE IF (LINE (LINDEX : LINDEX) .EQ. '"') THEN
         SYMBOL = 3
      ELSE
         SYMBOL = 4
      END IF

C///  PERFORM ACTION 1.

C     SYMBOL: 1 BLANK, 2 END OF LINE, 3 QUOTES, 4 OTHER
C     STATE: 1 SEARCHING, 2 BUILDING WORD, 3 BUILDING STRING
C     ACTION 1: 0 IGNORE CHARACTER, 1 COPY CHARACTER

C        SYMBOL   1  2  3  4    1  2  3  4    1  2  3  4
C         STATE   1  1  1  1    2  2  2  2    3  3  3  3
      DATA ACT1 / 0, 0, 0, 1,   0, 0, 0, 1,   1, 0, 0, 1 /

      IF (ACT1(SYMBOL, STATE) .EQ. 1) THEN
         IF (LINE (LINDEX : LINDEX) .NE. ' ') THEN
            IF (HOLD) THEN
               ERROR = .NOT. (WINDEX .LT. LEN (WORD))
               IF (ERROR) GO TO 9006
               WINDEX = WINDEX + 1
               IF (WINDEX .LE. 80) WORD (WINDEX : WINDEX) = ' '
            END IF
            ERROR = .NOT. (WINDEX .LT. LEN (WORD))
            IF (ERROR) GO TO 9006
            WINDEX = WINDEX + 1
            WORD (WINDEX : WINDEX) = LINE (LINDEX : LINDEX)
         END IF
         HOLD = LINE (LINDEX : LINDEX) .EQ. ' ' .AND. 0 .LT. WINDEX
      END IF

C///  PERFORM ACTION 2.

C     SYMBOL: 1 BLANK, 2 END OF LINE, 3 QUOTES, 4 OTHER
C     STATE: 1 SEARCHING, 2 BUILDING WORD, 3 BUILDING STRING
C     ACTION 2: 0 OLD CHARACTER, 1 NEW CHARACTER

C        SYMBOL   1  2  3  4    1  2  3  4    1  2  3  4
C         STATE   1  1  1  1    2  2  2  2    3  3  3  3
      DATA ACT2 / 1, 1, 1, 1,   0, 0, 0, 1,   1, 1, 1, 1 /

      IF (ACT2(SYMBOL, STATE) .EQ. 1) THEN
         IF (SYMBOL .EQ. 2) THEN
            LINDEX = LIMIT + 1
         ELSE
            LINDEX = LINDEX + 1
         END IF
      END IF

C///  PERFORM ACTION 3.

C     SYMBOL: 1 BLANK, 2 END OF LINE, 3 QUOTES, 4 OTHER
C     STATE: 1 SEARCHING, 2 BUILDING WORD, 3 BUILDING STRING
C     ACTION 3: 0 CONTINUE, 1 DONE

C        SYMBOL   1  2  3  4    1  2  3  4    1  2  3  4
C         STATE   1  1  1  1    2  2  2  2    3  3  3  3
      DATA ACT3 / 0, 0, 0, 0,   1, 1, 1, 0,   0, 0, 1, 0 /

      DONE = ACT3(SYMBOL, STATE) .EQ. 1

C///  PERFORM ACTION 4.

C     SYMBOL: 1 BLANK, 2 END OF LINE, 3 QUOTES, 4 OTHER
C     STATE: 1 SEARCHING, 2 BUILDING WORD, 3 BUILDING STRING

C        SYMBOL   1  2  3  4    1  2  3  4    1  2  3  4
C         STATE   1  1  1  1    2  2  2  2    3  3  3  3
      DATA ACT4 / 1, 1, 3, 2,   2, 2, 2, 2,   3, 3, 3, 3 /

      STATE = ACT4(SYMBOL, STATE)

C///  BOTTOM OF THE LOOP OVER THE CHARACTERS.

      IF (.NOT. DONE) GO TO 0100

C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  POSITION THE LINE FOR THE NEXT READ.

      IF (LINDEX .LE. LEN (LINE)) THEN
         TEMP = LINE (LINDEX : )
         LINE = TEMP 
      ELSE
         LINE = ' '
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, NUMBER
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, SCRIPT
      GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, SCRIPT, NUMBER, STATUS
      ERROR = .TRUE.
      GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, SCRIPT, NUMBER, STATUS
      ERROR = .TRUE.
      GO TO 99999

9005  IF (0 .LT. TEXT) WRITE (TEXT, 99005) 
     +   ID, SCRIPT, NUMBER, LENGTH, LIMIT - 2
      GO TO 99999

9006  IF (0 .LT. TEXT) WRITE (TEXT, 99006) 
     +   ID, SCRIPT, NUMBER, LEN (WORD), WORD (1 : MIN (54, WINDEX))
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  THE NUMBER OF LINES READ IS NEGATIVE.'
     +  //10X, I10, '  LINES READ')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  THE SCRIPT UNIT NUMBER MUST BE POSITIVE.'
     +  //10X, I10, '  UNIT NUMBER')

99003 FORMAT
     +   (/1X, A9, 'ERROR.  THE SCRIPT FILE ENDS UNEXPECTEDLY.'
     +  //10X, I10, '  UNIT NUMBER'
     +   /10X, I10, '  LINES SUCCESSFULLY READ'
     +   /10X, I10, '  I/O STATUS')

99004 FORMAT
     +   (/1X, A9, 'ERROR.  A READ FROM THE SCRIPT FILE FAILS.'
     +  //10X, I10, '  UNIT NUMBER'
     +   /10X, I10, '  LINES SUCCESSFULLY READ'
     +   /10X, I10, '  I/O STATUS')

99005 FORMAT
     +   (/1X, A9, 'ERROR.  A LINE IN THE SCRIPT FILE IS TOO LONG.'
     +  //10X, I10, '  UNIT NUMBER'
     +   /10X, I10, '  LINE NUMBER'
     +   /10X, I10, '  LINE LENGTH'
     +   /10X, I10, '  LENGTH LIMIT')

99006 FORMAT
     +   (/1X, A9, 'ERROR.  A WORD OR STRING IN THE SCRIPT FILE IS TOO'
     +   /10X, 'LONG.'
     +  //10X, I10, '  UNIT NUMBER'
     +   /10X, I10, '  LINE NUMBER'
     +   /10X, I10, '  LENGTH LIMIT'
     +  //10X, '   STRING:  ', A, ' ...')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE RESERV 
     +  (ERROR, TEXT, 
     +   FATAL, QLAST, QMAX, QSIZE, 
     +   NAME, NUMBER, Q)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     RESERV
C
C     RESERVE SPACE IN A WORKSPACE.
C
C     USES SUBROUTINE EXTENT.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9, NAME*(*)
      INTEGER LENGTH, NUMBER, TEXT
      LOGICAL ERROR, FATAL

      PARAMETER (ID = 'RESERV:  ')

C///  CHECK THE ARGUMENTS.

      ERROR =
     +   .NOT. (0 .LE. QLAST .AND. QLAST .LE. QMAX .AND. 0 .LE. QSIZE)
     +   .OR. (FATAL .AND. .NOT. (QMAX .LE. QSIZE))
      IF (ERROR) GO TO 9001

      ERROR = .NOT. (1 .LE. NUMBER)
      IF (ERROR) GO TO 9002

      ERROR = FATAL .AND. .NOT. (QLAST + NUMBER .LE. QSIZE)
      IF (ERROR) GO TO 9003

C///  RESERV THE SPACE.

      Q = QLAST + 1
      QLAST = QLAST + NUMBER
      QMAX = MAX (QLAST, QMAX)

C///  ERROR MESSAGES.

      GO TO 99999

9001  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, NAME)
         WRITE (TEXT, 99001) ID, QLAST, QMAX, QSIZE, NAME (1 : LENGTH)
      END IF
      GO TO 99999

9002  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, NAME)
         WRITE (TEXT, 99002) ID, NUMBER, NAME (1 : LENGTH)
      END IF
      GO TO 99999

9003  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, NAME)
         WRITE (TEXT, 99003) 
     +      ID, QSIZE, QLAST + NUMBER, NAME (1 : LENGTH)
      END IF
      GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  THE WORKSPACE POINTERS ARE INCONSISTENT.'
     +  //10X, I10, '  LAST ENTRY IN USE'
     +   /10X, I10, '  MAXIMUM ENTRY EVER IN USE'
     +   /10X, I10, '  SIZE OF THE WORKSPACE'
     +  //10X, 'NAME OF THE REQUEST:  ', A)

99002 FORMAT
     +  (/1X, A9, 'ERROR.  THE NUMBER TO RESERV MUST BE POSITIVE.'
     +  //10X, I10, '  NUMBER'
     +  //10X, 'NAME OF THE REQUEST:  ', A)

99003 FORMAT
     +  (/1X, A9, 'ERROR.  THE WORKSPACE IS TOO SMALL.'
     +  //10X, I10, '  SIZE OF THE WORKSPACE'
     +   /10X, I10, '  NEEDED SIZE (MORE MAY BE NEEDED LATER)'
     +  //10X, 'NAME OF THE REQUEST:  ', A)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE RIGHTJ (STRNG1, STRNG2, WIDTH)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     RIGHTJ
C
C     RIGHT JUSTIFY ONE STRING IN ANOTHER WITH BLANKS SQUEEZED OUT.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER CHAR*1, STRNG1*(*), STRNG2*(*)
      INTEGER J, LENGTH, LIMIT, OFFSET, TOTAL, WIDTH

C///  COPY AND SQUEEZE THE STRING.

      LENGTH = 0
      LIMIT = MIN (LEN (STRNG2), WIDTH)
      OFFSET = 1
      STRNG2 = ' '
      TOTAL = 0
      DO 0100 J = 1, LEN (STRNG1)
         CHAR = STRNG1 (J : J)
         IF (CHAR .NE. ' ') THEN
            TOTAL = TOTAL + OFFSET
            IF (LENGTH + OFFSET .LE. LIMIT) THEN
               LENGTH = LENGTH + OFFSET
               OFFSET = 1
               STRNG2 (LENGTH : LENGTH) = CHAR 
            END IF
         ELSE
            OFFSET = 2
         END IF
0100  CONTINUE

C///  RIGHT JUSTIFY THE STRING.

      IF (0 .LT. TOTAL .AND. TOTAL .LT. LIMIT) THEN
         DO 0200 J = TOTAL, 1, - 1
            CHAR = STRNG2(J : J)
            STRNG2(J + LIMIT - TOTAL : J + LIMIT - TOTAL) = CHAR
0200     CONTINUE
         STRNG2(1 : LIMIT - TOTAL) = ' '
      ELSE IF (LIMIT .LT. TOTAL) THEN
         DO 0300 J = MAX (1, LIMIT - 2), LIMIT
            STRNG2(J : J) = '.'
0300     CONTINUE
      END IF
         
C///  EXIT.

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
      CHARACTER CHAR*1, STRING*(*)
      INTEGER J, LENGTH
      LOGICAL BLANK

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
         IF (LENGTH .LT. LEN (STRING)) STRING (LENGTH + 1 : ) = ' '
      ELSE
         LENGTH = 1
      END IF

C///  EXIT.

      RETURN
      END
      SUBROUTINE VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     VERIFY
C
C     VERIFY THAT THE NEXT PHRASE IN A SCRIPT FILE IS THE ONE EXPECTED.
C
C     USES SUBROUTINES EXTENT AND READW.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER ID*9, LINE*(*), PHRASE*132, TEMP*132, WORD*(*)
      INTEGER LENGTH, NUMBER, PINDEX, PRINT, SCRIPT, TEXT, WINDEX
      LOGICAL ERROR

      PARAMETER (ID = 'VERIFY:  ')

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE ARGUMENTS.

      CALL EXTENT (WINDEX, WORD)
      ERROR = .NOT. (WINDEX .LE. LEN (PHRASE))
      IF (ERROR) GO TO 9001

C///////////////////////////////////////////////////////////////////////
C
C     READ THE PHRASE.
C
C///////////////////////////////////////////////////////////////////////

C///  INITIALIZE.

      PHRASE = ' '
      PINDEX = 0

C///  TOP OF THE LOOP OVER THE STRINGS.

0100  CONTINUE

C///  APPEND A STRING TO THE PHRASE.

C     SUBROUTINE READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)

      CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, TEMP)
      IF (ERROR) GO TO 9002

      CALL EXTENT (LENGTH, TEMP)

      IF (PINDEX .EQ. 0) THEN
         PHRASE = TEMP
         PINDEX = LENGTH
      ELSE IF (PINDEX + 2 .LE. LEN (PHRASE)) THEN
         PHRASE (PINDEX + 2 : ) = TEMP
         PINDEX = MIN (LEN (PHRASE), PINDEX + 1 + LENGTH)
      END IF

C///  BOTTOM OF THE LOOP OVER THE STRINGS.

      IF (PINDEX .LT. WINDEX) GO TO 0100

C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

      ERROR = .NOT. (PHRASE .EQ. WORD)
      IF (ERROR) GO TO 9003

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, WINDEX, LEN (PHRASE)
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID
      GO TO 99999

9003  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99003) ID, SCRIPT, NUMBER
         IF (58 .LT. PINDEX) THEN
            PHRASE (54 : ) = ' ...'
            PINDEX = 58
         END IF
         WRITE (TEXT, '(10X, A, A)') '    FOUND:  ', PHRASE (1 : PINDEX)
         PHRASE = WORD
         IF (58 .LT. WINDEX) THEN
            PHRASE (54 : ) = ' ...'
            WINDEX = 58
         END IF
         WRITE (TEXT, '(10X, A, A)') ' EXPECTED:  ', PHRASE (1 : WINDEX)
      END IF
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  THE PHRASE IS TOO LONG.'
     +  //10X, I10, '  LENGTH'
     +   /10X, I10, '  LIMIT')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  READW FAILS.')

99003 FORMAT
     +   (/1X, A9, 'ERROR.  THE NEXT PHRASE IN THE SCRIPT FILE IS NOT'
     +   /10X, 'THE ONE EXPECTED.'
     +  //10X, I10, '  UNIT NUMBER'
     +   /10X, I10, '  LINE NUMBER'
     +  /)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE WRT1D1
     +  (ERROR, TEXT,
     +   DATA, LEGE, NAME, SIZE, SUBT, TITL, XCOORD, XINC, XLAB, YCOORD,
     +   YINC, YLAB)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     WRT1D1
C
C     WRITE SINGLE PRECISION DATA TO A RECORD IN A PLOT1D FILE.
C
C     USES SUBROUTINE EXTENT.
C
C///////////////////////////////////////////////////////////////////////
C
C     CHANGE HISTORY:
C
C     12 APR 90  CORRECTED ORDER OF ARGUMENT LIST.
C     13 APR 90  CHANGE TO ACCEPT CURVES WITH NO POINTS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   ID*9, LEGE*(*), NAME*(*), SUBT*(*), TITL*(*), XLAB*(*), 
     +   YLAB*(*)
      INTEGER DATA, J, LENGTH, SIZE, STATUS, TEXT, XINC, YINC
      LOGICAL ERROR
      REAL XCOORD, YCOORD

      PARAMETER (ID = 'WRT1D1:  ')

      DIMENSION XCOORD(1 + SIZE * XINC), YCOORD(1 + SIZE * YINC)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (0 .LE. SIZE)
      IF (ERROR) GO TO 9001

      ERROR = .NOT. (1 .LE. XINC)
      IF (ERROR) GO TO 9002

      ERROR = .NOT. (1 .LE. YINC)
      IF (ERROR) GO TO 9003

C///////////////////////////////////////////////////////////////////////
C
C     WRITE THE RECORD.
C
C///////////////////////////////////////////////////////////////////////

C///  SPACE.

      WRITE (DATA, 10001, ERR = 9004, IOSTAT = STATUS)

C///  NAME.

      CALL EXTENT (LENGTH, NAME)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   '    NAME', NAME (1 : LENGTH)

C///  TITLE.

      CALL EXTENT (LENGTH, TITL)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   '   TITLE', TITL (1 : LENGTH)

C///  SUBTITLE.

      CALL EXTENT (LENGTH, SUBT)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   'SUBTITLE', SUBT (1 : LENGTH)

C///  LEGEND.

      CALL EXTENT (LENGTH, LEGE)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   '  LEGEND', LEGE (1 : LENGTH)

C///  X LABEL.

      CALL EXTENT (LENGTH, XLAB)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   ' X LABEL', XLAB (1 : LENGTH)

C///  Y LABEL.

      CALL EXTENT (LENGTH, YLAB)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   ' Y LABEL', YLAB (1 : LENGTH)

C///  SIZE.

      WRITE (DATA, 10003, ERR = 9004, IOSTAT = STATUS)
     +   '    SIZE', SIZE

C///  COORDINATES.

      IF (0 .LT. SIZE) WRITE (DATA, 10004, ERR = 9004, IOSTAT = STATUS)
     +   (1 + J, REAL (XCOORD(1 + J * XINC)),
     +   REAL (YCOORD(1 + J * YINC)), J = 0, SIZE - 1)

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT ()

10002 FORMAT (1X, A8, ' : ', A)

10003 FORMAT (1X, A8, ' : ', I10)

10004 FORMAT (1X, I8, ' : ', 1P, E15.7, 2X, E15.7)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, SIZE
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, XINC
      GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, YINC
      GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, DATA, STATUS
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  THE NUMBER OF DATA POINTS IN THE CURVE MUST'
     +   /10X, 'BE POSITIVE OR ZERO.'
     +  //10X, I10, '  SIZE')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  THE X INCREMENT MUST BE POSITIVE.'
     +  //10X, I10, '  XINC')

99003 FORMAT
     +   (/1X, A9, 'ERROR.  THE Y INCREMENT MUST BE POSITIVE.'
     +  //10X, I10, '  YINC')

99004 FORMAT
     +   (/1X, A9, 'ERROR.  FORMATTED WRITE TO THE DATA FILE FAILS.'
     +   //10X, I10, '  UNIT NUMBER'
     +    /10X, I10, '  I/O STATUS')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE WRT1D2
     +  (ERROR, TEXT,
     +   DATA, LEGE, NAME, SIZE, SUBT, TITL, XCOORD, XINC, XLAB, YCOORD,
     +   YINC, YLAB)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     WRT1D2
C
C     WRITE DOUBLE PRECISION DATA TO A RECORD IN A PLOT1D FILE.
C
C     USES SUBROUTINE EXTENT.
C
C///////////////////////////////////////////////////////////////////////
C
C     CHANGE HISTORY:
C
C     12 APR 90  CORRECTED ORDER OF ARGUMENT LIST.
C     13 APR 90  CHANGE TO ACCEPT CURVES WITH NO POINTS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   ID*9, LEGE*(*), NAME*(*), SUBT*(*), TITL*(*), XLAB*(*), 
     +   YLAB*(*)
      INTEGER DATA, J, LENGTH, SIZE, STATUS, TEXT, XINC, YINC
      LOGICAL ERROR
      DOUBLE PRECISION XCOORD, YCOORD

      PARAMETER (ID = 'WRT1D2:  ')

      DIMENSION XCOORD(1 + SIZE * XINC), YCOORD(1 + SIZE * YINC)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (0 .LE. SIZE)
      IF (ERROR) GO TO 9001

      ERROR = .NOT. (1 .LE. XINC)
      IF (ERROR) GO TO 9002

      ERROR = .NOT. (1 .LE. YINC)
      IF (ERROR) GO TO 9003

C///////////////////////////////////////////////////////////////////////
C
C     WRITE THE RECORD.
C
C///////////////////////////////////////////////////////////////////////

C///  SPACE.

      WRITE (DATA, 10001, ERR = 9004, IOSTAT = STATUS)

C///  NAME.

      CALL EXTENT (LENGTH, NAME)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   '    NAME', NAME (1 : LENGTH)

C///  TITLE.

      CALL EXTENT (LENGTH, TITL)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   '   TITLE', TITL (1 : LENGTH)

C///  SUBTITLE.

      CALL EXTENT (LENGTH, SUBT)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   'SUBTITLE', SUBT (1 : LENGTH)

C///  LEGEND.

      CALL EXTENT (LENGTH, LEGE)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   '  LEGEND', LEGE (1 : LENGTH)

C///  X LABEL.

      CALL EXTENT (LENGTH, XLAB)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   ' X LABEL', XLAB (1 : LENGTH)

C///  Y LABEL.

      CALL EXTENT (LENGTH, YLAB)
      WRITE (DATA, 10002, ERR = 9004, IOSTAT = STATUS)
     +   ' Y LABEL', YLAB (1 : LENGTH)

C///  SIZE.

      WRITE (DATA, 10003, ERR = 9004, IOSTAT = STATUS)
     +   '    SIZE', SIZE

C///  COORDINATES.

      IF (0 .LT. SIZE) WRITE (DATA, 10004, ERR = 9004, IOSTAT = STATUS)
     +   (1 + J, REAL (XCOORD(1 + J * XINC)),
     +   REAL (YCOORD(1 + J * YINC)), J = 0, SIZE - 1)

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT ()

10002 FORMAT (1X, A8, ' : ', A)

10003 FORMAT (1X, A8, ' : ', I10)

10004 FORMAT (1X, I8, ' : ', 1P, E15.7, 2X, E15.7)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, SIZE
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, XINC
      GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, YINC
      GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, DATA, STATUS
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  THE NUMBER OF DATA POINTS IN THE CURVE MUST'
     +   /10X, 'BE POSITIVE OR ZERO.'
     +  //10X, I10, '  SIZE')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  THE X INCREMENT MUST BE POSITIVE.'
     +  //10X, I10, '  XINC')

99003 FORMAT
     +   (/1X, A9, 'ERROR.  THE Y INCREMENT MUST BE POSITIVE.'
     +  //10X, I10, '  YINC')

99004 FORMAT
     +   (/1X, A9, 'ERROR.  FORMATTED WRITE TO THE DATA FILE FAILS.'
     +   //10X, I10, '  UNIT NUMBER'
     +    /10X, I10, '  I/O STATUS')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE WRT2D1
     +  (ERROR, TEXT,
     +   DATA, LEGE, NAME, SUBT, TITL, XCOORD, XINC, XLAB, XSIZE,
     +   YCOORD, YINC, YLAB, YSIZE, ZCOORD, ZLAB)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     WRT2D1
C
C     WRITE SINGLE PRECISION DATA TO A RECORD IN A PLOT2D FILE.
C
C     USES SUBROUTINE EXTENT.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   ID*9, LEGE*(*), NAME*(*), SUBT*(*), TITL*(*), XLAB*(*),
     +   YLAB*(*), ZLAB*(*)
      INTEGER
     +   DATA, J, K, LAST, LENGTH, STATUS, TEXT, XINC, XSIZE, YINC,
     +   YSIZE
      LOGICAL ERROR
      REAL XCOORD, YCOORD, ZCOORD

      PARAMETER (ID = 'WRT2D1:  ')

      DIMENSION
     +   XCOORD(XSIZE), YCOORD(YSIZE),
     +   ZCOORD(1 + (XSIZE - 1) * XINC + (YSIZE - 1) * YINC)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (0 .LE. XSIZE)
      IF (ERROR) GO TO 9001

      ERROR = .NOT. (1 .LE. XINC)
      IF (ERROR) GO TO 9002

      ERROR = .NOT. (0 .LE. YSIZE)
      IF (ERROR) GO TO 9003

      ERROR = .NOT. (1 .LE. YINC)
      IF (ERROR) GO TO 9004

C///////////////////////////////////////////////////////////////////////
C
C     WRITE THE RECORD.
C
C///////////////////////////////////////////////////////////////////////

C///  SPACE.

      WRITE (DATA, 10001, ERR = 9005, IOSTAT = STATUS)

C///  NAME.

      CALL EXTENT (LENGTH, NAME)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   '    NAME', NAME (1 : LENGTH)

C///  TITLE.

      CALL EXTENT (LENGTH, TITL)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   '   TITLE', TITL (1 : LENGTH)

C///  SUBTITLE.

      CALL EXTENT (LENGTH, SUBT)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   'SUBTITLE', SUBT (1 : LENGTH)

C///  LEGEND.

      CALL EXTENT (LENGTH, LEGE)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   '  LEGEND', LEGE (1 : LENGTH)

C///  X LABEL.

      CALL EXTENT (LENGTH, XLAB)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   ' X LABEL', XLAB (1 : LENGTH)

C///  Y LABEL.

      CALL EXTENT (LENGTH, YLAB)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   ' Y LABEL', YLAB (1 : LENGTH)

C///  Z LABEL.

      CALL EXTENT (LENGTH, ZLAB)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   ' Z LABEL', ZLAB (1 : LENGTH)

C///  X SIZE.

      WRITE (DATA, 10003, ERR = 9005, IOSTAT = STATUS)
     +   '  X SIZE', XSIZE

C///  Y SIZE.

      WRITE (DATA, 10003, ERR = 9005, IOSTAT = STATUS)
     +   '  Y SIZE', YSIZE

C///  X, Y AND Z COORDINATES.

      LAST = MIN (XSIZE, YSIZE)

      IF (0 .LT. LAST) THEN
         WRITE (DATA, 10004, ERR = 9005, IOSTAT = STATUS)
     +   (J, REAL (XCOORD(J)), REAL (YCOORD(J)), J = 1, LAST)

         IF (YSIZE .LT. XSIZE)
     +      WRITE (DATA, 10005, ERR = 9005, IOSTAT = STATUS)
     +      (J, REAL (XCOORD(J)), J = YSIZE + 1, XSIZE)

         IF (XSIZE .LT. YSIZE)
     +      WRITE (DATA, 10006, ERR = 9005, IOSTAT = STATUS)
     +      (J, REAL (YCOORD(J)), J = XSIZE + 1, YSIZE)

         WRITE (DATA, 10007, ERR = 9005, IOSTAT = STATUS)
     +      ((J + 1, K + 1, REAL (ZCOORD(1 + J * XINC + K * YINC)),
     +      J = 0, XSIZE - 1), K = 0, YSIZE - 1)
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT ()

10002 FORMAT (1X, A8, ' : ', A)

10003 FORMAT (1X, A8, ' : ', I10)

10004 FORMAT (1X, I8, ' : ', 1P, E15.7, 2X, E15.7)

10005 FORMAT (1X, I8, ' : ', 1P, E15.7)

10006 FORMAT (1X, I8, ' : ', 1P, 17X, E15.7)

10007 FORMAT (1X, I8, 2X, I8, ' : ', 1P, E15.7)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, XSIZE
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, XINC
      GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, YSIZE
      GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, YINC
      GO TO 99999

9005  IF (0 .LT. TEXT) WRITE (TEXT, 99005) ID, DATA, STATUS
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  THE NUMBER OF DATA POINTS IN THE CURVE MUST'
     +   /10X, 'BE POSITIVE OR ZERO.'
     +  //10X, I10, '  X SIZE')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  THE X INCREMENT MUST BE POSITIVE.'
     +  //10X, I10, '  XINC')

99003 FORMAT
     +   (/1X, A9, 'ERROR.  THE NUMBER OF DATA POINTS IN THE CURVE MUST'
     +   /10X, 'BE POSITIVE OR ZERO.'
     +  //10X, I10, '  Y SIZE')

99004 FORMAT
     +   (/1X, A9, 'ERROR.  THE Y INCREMENT MUST BE POSITIVE.'
     +  //10X, I10, '  YINC')

99005 FORMAT
     +   (/1X, A9, 'ERROR.  FORMATTED WRITE TO THE DATA FILE FAILS.'
     +   //10X, I10, '  UNIT NUMBER'
     +    /10X, I10, '  I/O STATUS')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE WRT2D2
     +  (ERROR, TEXT,
     +   DATA, LEGE, NAME, SUBT, TITL, XCOORD, XINC, XLAB, XSIZE,
     +   YCOORD, YINC, YLAB, YSIZE, ZCOORD, ZLAB)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     WRT2D2
C
C     WRITE DOUBLE PRECISION DATA TO A RECORD IN A PLOT2D FILE.
C
C     USES SUBROUTINE EXTENT.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   ID*9, LEGE*(*), NAME*(*), SUBT*(*), TITL*(*), XLAB*(*),
     +   YLAB*(*), ZLAB*(*)
      INTEGER
     +   DATA, J, K, LAST, LENGTH, STATUS, TEXT, XINC, XSIZE, YINC,
     +   YSIZE
      LOGICAL ERROR
      DOUBLE PRECISION XCOORD, YCOORD, ZCOORD

      PARAMETER (ID = 'WRT2D2:  ')

      DIMENSION
     +   XCOORD(XSIZE), YCOORD(YSIZE),
     +   ZCOORD(1 + (XSIZE - 1) * XINC + (YSIZE - 1) * YINC)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (0 .LE. XSIZE)
      IF (ERROR) GO TO 9001

      ERROR = .NOT. (1 .LE. XINC)
      IF (ERROR) GO TO 9002

      ERROR = .NOT. (0 .LE. YSIZE)
      IF (ERROR) GO TO 9003

      ERROR = .NOT. (1 .LE. YINC)
      IF (ERROR) GO TO 9004

C///////////////////////////////////////////////////////////////////////
C
C     WRITE THE RECORD.
C
C///////////////////////////////////////////////////////////////////////

C///  SPACE.

      WRITE (DATA, 10001, ERR = 9005, IOSTAT = STATUS)

C///  NAME.

      CALL EXTENT (LENGTH, NAME)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   '    NAME', NAME (1 : LENGTH)

C///  TITLE.

      CALL EXTENT (LENGTH, TITL)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   '   TITLE', TITL (1 : LENGTH)

C///  SUBTITLE.

      CALL EXTENT (LENGTH, SUBT)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   'SUBTITLE', SUBT (1 : LENGTH)

C///  LEGEND.

      CALL EXTENT (LENGTH, LEGE)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   '  LEGEND', LEGE (1 : LENGTH)

C///  X LABEL.

      CALL EXTENT (LENGTH, XLAB)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   ' X LABEL', XLAB (1 : LENGTH)

C///  Y LABEL.

      CALL EXTENT (LENGTH, YLAB)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   ' Y LABEL', YLAB (1 : LENGTH)

C///  Z LABEL.

      CALL EXTENT (LENGTH, ZLAB)
      WRITE (DATA, 10002, ERR = 9005, IOSTAT = STATUS)
     +   ' Z LABEL', ZLAB (1 : LENGTH)

C///  X SIZE.

      WRITE (DATA, 10003, ERR = 9005, IOSTAT = STATUS)
     +   '  X SIZE', XSIZE

C///  Y SIZE.

      WRITE (DATA, 10003, ERR = 9005, IOSTAT = STATUS)
     +   '  Y SIZE', YSIZE

C///  X, Y AND Z COORDINATES.

      LAST = MIN (XSIZE, YSIZE)

      IF (0 .LT. LAST) THEN
         WRITE (DATA, 10004, ERR = 9005, IOSTAT = STATUS)
     +   (J, REAL (XCOORD(J)), REAL (YCOORD(J)), J = 1, LAST)

         IF (YSIZE .LT. XSIZE)
     +      WRITE (DATA, 10005, ERR = 9005, IOSTAT = STATUS)
     +      (J, REAL (XCOORD(J)), J = YSIZE + 1, XSIZE)

         IF (XSIZE .LT. YSIZE)
     +      WRITE (DATA, 10006, ERR = 9005, IOSTAT = STATUS)
     +      (J, REAL (YCOORD(J)), J = XSIZE + 1, YSIZE)

         WRITE (DATA, 10007, ERR = 9005, IOSTAT = STATUS)
     +      ((J + 1, K + 1, REAL (ZCOORD(1 + J * XINC + K * YINC)),
     +      J = 0, XSIZE - 1), K = 0, YSIZE - 1)
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT ()

10002 FORMAT (1X, A8, ' : ', A)

10003 FORMAT (1X, A8, ' : ', I10)

10004 FORMAT (1X, I8, ' : ', 1P, E15.7, 2X, E15.7)

10005 FORMAT (1X, I8, ' : ', 1P, E15.7)

10006 FORMAT (1X, I8, ' : ', 1P, 17X, E15.7)

10007 FORMAT (1X, I8, 2X, I8, ' : ', 1P, E15.7)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, XSIZE
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, XINC
      GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, YSIZE
      GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, YINC
      GO TO 99999

9005  IF (0 .LT. TEXT) WRITE (TEXT, 99005) ID, DATA, STATUS
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  THE NUMBER OF DATA POINTS IN THE CURVE MUST'
     +   /10X, 'BE POSITIVE OR ZERO.'
     +  //10X, I10, '  X SIZE')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  THE X INCREMENT MUST BE POSITIVE.'
     +  //10X, I10, '  XINC')

99003 FORMAT
     +   (/1X, A9, 'ERROR.  THE NUMBER OF DATA POINTS IN THE CURVE MUST'
     +   /10X, 'BE POSITIVE OR ZERO.'
     +  //10X, I10, '  Y SIZE')

99004 FORMAT
     +   (/1X, A9, 'ERROR.  THE Y INCREMENT MUST BE POSITIVE.'
     +  //10X, I10, '  YINC')

99005 FORMAT
     +   (/1X, A9, 'ERROR.  FORMATTED WRITE TO THE DATA FILE FAILS.'
     +   //10X, I10, '  UNIT NUMBER'
     +    /10X, I10, '  I/O STATUS')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
