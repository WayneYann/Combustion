C     CVS $Revision: 1.1.1.1 $ reposited $Date: 2006/05/26 19:09:34 $

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     VERSION 3.29 OF APRIL 1998
C
C     THE TWOPNT PROGRAM FOR BOUNDARY VALUE PROBLEMS
C
C     WRITTEN BY: DR. JOSEPH F. GRCAR
C                 SANDIA NATIONAL LABORATORIES
C                 MAIL STOP 9051
C                 LIVERMORE, CALIFORNIA  94551-0969  USA
C
C                 (925) 294-2662
C                 (FTS) 234-2662
C
C                 na.grcar@na-net.ornl.gov
C                 sepp@california.sandia.gov
C
C///////////////////////////////////////////////////////////////////////
C
C     DOCUMENTATION:
C
C     J. F. Grcar, "The Twopnt Program for Boundary Value Problems,"
C     Sandia National Laboratories Report SAND91-8230, Livermore,
C     California, April 1992.  Reprinted February 1996.
C
C///////////////////////////////////////////////////////////////////////
C
C     CHANGES FROM THE PREVIOUS VERSION:
C
C     1) PUT CHANGE BLOCK AROUND DECLARATION OF SAME IN TWEPS.
C
C     2) ALTER AREA CODE.
C
C///////////////////////////////////////////////////////////////////////

      SUBROUTINE EVOLVE
     +  (ERROR, TEXT,
     +   ABOVE, BELOW, BUFFER, COMPS, CONDIT, DESIRE, GROUPA, GROUPB,
     +   LEVELD, LEVELM, NAME, NAMES, POINTS, REPORT, S0, S1, SIGNAL,
     +   STEP, STEPS2, STRID0, STRIDE, SUCCES, TDABS, TDAGE, TDEC,
     +   TDREL, TIME, TINC, TMAX, TMIN, V0, V1, VSAVE, Y0, Y1, YNORM)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     EVOLVE
C
C     PERFORM TIME EVOLUTION.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER
     +   CWORD*80, HEADER*80, ID*9, JWORD*80, NAME*(*), REMARK*80,
     +   SIGNAL*(*), YWORD*80
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   ABOVE, BELOW, BUFFER, CHANGE, CONDIT, CSAVE, DUMMY, HIGH, LOW,
     +   S0, S1, STRID0, STRIDE, TDABS, TDEC, TDREL, TINC, TMAX, TMIN,
     +   V0, V1, VSAVE, Y0, Y1, YNORM
      EXTERNAL
     +   SEARCH, TWCOPY, TWLOGR, TWNORM
      INTEGER
     +   AGE, AGEJ, COMPS, COUNT, DESIRE, FIRST, GROUPA, GROUPB, J,
     +   LAST, LENGTH, LEVELD, LEVELM, NAMES, NUMBER, POINTS, QBNDS,
     +   QDVRG, QNULL, REPORT, ROUTE, STEP, STEPS2, TDAGE, TEXT,
     +   XREPOR
      INTRINSIC
     +   LOG10, MAX, MIN
      LOGICAL
     +   ERROR, EXIST, JACOB, MESS, SUCCES, TIME, XSUCCE

      PARAMETER (ID = 'EVOLVE:  ')

C     REPORT CODES
      PARAMETER (QNULL = 0, QBNDS = 1, QDVRG = 2)

      DIMENSION
     +   ABOVE(GROUPA + COMPS * POINTS + GROUPB),
     +   BELOW(GROUPA + COMPS * POINTS + GROUPB),
     +   BUFFER(GROUPA + COMPS * POINTS + GROUPB), HEADER(2, 3),
     +   NAME(NAMES), S0(GROUPA + COMPS * POINTS + GROUPB),
     +   V0(GROUPA + COMPS * POINTS + GROUPB),
     +   V1(GROUPA + COMPS * POINTS + GROUPB),
     +   VSAVE(GROUPA + COMPS * POINTS + GROUPB),
     +   Y0(GROUPA + COMPS * POINTS + GROUPB)

C///  SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.

      SAVE

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  INITIALIZE.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      MESS = .FALSE.

C     TURN OF REVERSE COMMUNICATION FLAGS.
      TIME = .FALSE.

C     TURN OFF ALL COMPLETION STATUS FLAGS.
      ERROR = .FALSE.
      REPORT = QNULL
      SUCCES = .FALSE.

C///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      IF (SIGNAL .NE. ' ') THEN
         GO TO (1010, 1020, 1060, 1080, 1090, 2020) ROUTE
         ERROR = .TRUE.
         GO TO 9001
      END IF

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (((0 .LT. COMPS) .EQV. (0 .LT. POINTS)) .AND.
     +   0 .LE. COMPS .AND. 0 .LE. POINTS .AND. 0 .LE. GROUPA .AND.
     +   0 .LE. GROUPB .AND. 0 .LT. GROUPA + COMPS * POINTS + GROUPB)
      IF (ERROR) GO TO 9002

      ERROR = .NOT. (0 .LT. DESIRE)
      IF (ERROR) GO TO 9003

      ERROR = .NOT. (1.0 .LE. TDEC .AND. 1.0 .LE. TINC)
      IF (ERROR) GO TO 9004

      ERROR = .NOT. (0.0 .LT. TMIN .AND. TMIN .LE. TMAX)
      IF (ERROR) GO TO 9005

      ERROR = .NOT. (TMIN .LE. STRID0 .AND. STRID0 .LE. TMAX)
      IF (ERROR) GO TO 9006

      ERROR = .NOT. (0 .LE. STEP)
      IF (ERROR) GO TO 9007

      ERROR = 1.0 .LT. TINC .AND. .NOT. 0 .LT. STEPS2
      IF (ERROR) GO TO 9008

C///  WRITE ALL MESSAGES.

C                     123456789_123456789_123456789_123456789_1234
C                     123456   123456   123456   123456   12345
      HEADER(1, 1) = '  TIME   LOG10                      NEWTON S'
      HEADER(1, 2) = ' POINT   ------------------------   --------'
      HEADER(1, 3) = 'NUMBER   NORM F   CHANGE   STRIDE   STEPS   '

C                     123456789_123456789_1
C                     123   123456   123456
      HEADER(2, 1) = 'EARCH                '
      HEADER(2, 2) = '---------------------'
      HEADER(2, 3) = 'J''S   COND J   REMARK'

      IF (MESS .AND. 0 .LT. TEXT) THEN
         ROUTE = 0
         YNORM = 1.0E-4
         CALL TWLOGR (YWORD, YNORM)

         WRITE (TEXT, 10001) ID, HEADER, STEP, YWORD
         WRITE (TEXT, 20001) ID, STEP, YWORD, LOG10 (STRID0)
         WRITE (TEXT, 10002) ID, HEADER, STEP, YWORD
         WRITE (TEXT, 20002) ID, STEP, YWORD, LOG10 (STRID0)
         WRITE (TEXT, 20003) ID, STEP, YWORD, LOG10 (STRID0)
         WRITE (TEXT, 10006) ID
         WRITE (TEXT, 10008) ID
         WRITE (TEXT, 20007) ID, STEP, YWORD
         WRITE (TEXT, 20004) ID, STEP, YWORD, LOG10 (STRID0)
         WRITE (TEXT, 10007) ID
         WRITE (TEXT, 20006) ID, STEP, YWORD
         WRITE (TEXT, 20008) ID
         WRITE (TEXT, 20005) ID, STEP, YWORD, LOG10 (STRID0)

         GO TO 9001
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     TIME EVOLUTION.
C
C///////////////////////////////////////////////////////////////////////

C///  0 < M?

      IF (0 .LT. STEP) GO TO 1010
         STRIDE = STRID0
         AGE = 0

C        RETAIN THE LATEST SOLUTION FOR USE BY THE FUNCTION
         CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, V0, BUFFER)
         SIGNAL = 'RETAIN'
C        GO TO 1010 WHEN ROUTE = 1
         ROUTE = 1
         GO TO 99999
1010  CONTINUE
      SIGNAL = ' '

C///  FIRST := STEP, LAST := STEP + DESIRE.

      EXIST = .FALSE.
      FIRST = STEP
      LAST = STEP + DESIRE

C///  PRINT.

      IF (.NOT. (0 .LT. LEVELM .AND. 0 .LT. TEXT)) GO TO 1030
         CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, V0, BUFFER)
         SIGNAL = 'RESIDUAL'
         TIME = .FALSE.
C        GO TO 1020 WHEN ROUTE = 2
         ROUTE = 2
         GO TO 99999
1020     CONTINUE
         SIGNAL = ' '
         CALL TWNORM (GROUPA + COMPS * POINTS + GROUPB, YNORM, BUFFER)
         CALL TWLOGR (YWORD, YNORM)

         IF (1 .EQ. LEVELM) THEN
            IF (STEP .EQ. 0) THEN
               WRITE (TEXT, 10001) ID, HEADER, STEP, YWORD
            ELSE
               WRITE (TEXT, 10002) ID, HEADER, STEP, YWORD
            END IF
         ELSE IF (1 .LT. LEVELM .AND. 0 .EQ. STEP) THEN
            IF (STEP .EQ. 0) THEN
               WRITE (TEXT, 20001) ID, STEP, YWORD, LOG10 (STRIDE)
            ELSE
               WRITE (TEXT, 20002) ID, STEP, YWORD, LOG10 (STRIDE)
            END IF
         END IF
1030  CONTINUE

C///  LOW := TMIN, HIGH := TMAX.

1040  CONTINUE

      LOW = TMIN
      HIGH = TMAX

C///  IF AGE = STEPS2 AND STRIDE < HIGH AND 1 < TINC, THEN INCREASE
C///  STRIDE.

      IF (AGE .EQ. STEPS2 .AND. STRIDE .LT. HIGH .AND. 1.0 .LT. TINC)
     +   THEN
         AGE = 0
         EXIST = .FALSE.
         LOW = STRIDE * TDEC
         STRIDE = MIN (HIGH, STRIDE * TINC)
         IF (1 .LT. LEVELM .AND. 0 .LT. TEXT)
     +      WRITE (TEXT, 20003) ID, STEP, YWORD, LOG10 (STRIDE)
      ELSE
         IF (1 .LT. LEVELM .AND. 0 .LT. TEXT .AND. 0 .LT. STEP)
     +      WRITE (TEXT, 20002) ID, STEP, YWORD, LOG10 (STRIDE)
      END IF

C///  NEWTON SEARCH.

1050  CONTINUE

C     STORE THE LATEST SOLUTION SHOULD THE SEARCH FAIL
      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, V0, VSAVE)

      COUNT = 0
      CSAVE = 0.0
      JWORD = ' '
      JACOB = .FALSE.

1060  CONTINUE

      IF (JACOB) THEN
         EXIST = .TRUE.
         COUNT = COUNT + 1
         CSAVE = MAX (CONDIT, CSAVE)
         IF (CSAVE .EQ. 0.0) THEN
            WRITE (JWORD, '(I3, 3X, A6)') COUNT, '    NA'
         ELSE
            WRITE (JWORD, '(I3, 3X, F6.2)') COUNT, LOG10 (CSAVE)
         END IF
      END IF

C     SUBROUTINE SEARCH
C    +  (ERROR, TEXT,
C    +   ABOVE, AGE, BELOW, BUFFER, COMPS, CONDIT, EXIST, GROUPA,
C    +   GROUPB, LEVELD, LEVELM, NAME, NAMES, POINTS, REPORT, S0, S1,
C    +   SIGNAL, STEPS, SUCCES, V0, V1, XXABS, XXAGE, XXREL, Y0, Y0NORM,
C    +   Y1)

      CALL SEARCH
     +  (ERROR, TEXT,
     +   ABOVE, AGEJ, BELOW, BUFFER, COMPS, CONDIT, EXIST, GROUPA,
     +   GROUPB, LEVELD - 1, LEVELM - 1, NAME, NAMES, POINTS, XREPOR,
     +   S0, S1, SIGNAL, NUMBER, XSUCCE, V0, V1, TDABS, TDAGE, TDREL,
     +   Y0, DUMMY, Y1)
      IF (ERROR) GO TO 9009

      IF (SIGNAL .NE. ' ') THEN
         JACOB = SIGNAL .EQ. 'PREPARE'
         TIME = .TRUE.
C        GO TO 1060 WHEN ROUTE = 3
         ROUTE = 3
         GO TO 99999
      END IF

C///  UNSUCCESSFUL?

      IF (.NOT. XSUCCE) THEN
         IF (1 .EQ. LEVELM .AND. 0 .LT. TEXT) THEN
            IF (XREPOR .EQ. QBNDS) THEN
               LENGTH = 6
               REMARK = 'BOUNDS'
            ELSE IF (XREPOR .EQ. QDVRG) THEN
               LENGTH = 7
               REMARK = 'DIVERGE'
            ELSE
               LENGTH = 1
               REMARK = ' '
            END IF
            WRITE (TEXT, 10003) STEP + 1, LOG10 (STRIDE), NUMBER, JWORD,
     +         REMARK (1 : LENGTH)
         END IF

C///  IF ALSO LOW < STRIDE AND 1 < TDEC, THEN DECREASE STRIDE.

         IF (LOW .LT. STRIDE .AND. 1.0 .LT. TDEC) THEN
            AGE = 0
            CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, VSAVE, V0)
            EXIST = .FALSE.
            HIGH = STRIDE / TINC
            STRIDE = MAX (LOW, STRIDE / TDEC)
            IF (1 .LT. LEVELM .AND. 0 .LT. TEXT) WRITE (TEXT, 20004)
     +         ID, STEP, YWORD, LOG10 (STRIDE)
            GO TO 1050
         END IF

C///  OTHERWISE END, FAILURE.

         GO TO 2010
      END IF

C///  IF NO CHANGE AND STRIDE .LT. HIGH AND 1.0 .LT. TINC, THEN
C///  INCREASE STRIDE.  OTHERWISE END, FAILURE.

      DO 1070 J = 1, GROUPA + COMPS * POINTS + GROUPB
         BUFFER(J) = V0(J) - VSAVE(J)
1070  CONTINUE
      CALL TWNORM (GROUPA + COMPS * POINTS + GROUPB, CHANGE, BUFFER)
      CALL TWLOGR (CWORD, CHANGE)

      IF (CHANGE .EQ. 0.0) THEN
         IF (1 .EQ. LEVELM .AND. 0 .LT. TEXT) THEN
            WRITE (TEXT, 10004)
     +         STEP + 1, '  ZERO', LOG10 (STRIDE), NUMBER, JWORD
         END IF

         IF (1.0 .LT. TINC .AND. STRIDE .LT. HIGH) THEN
            AGE = 0
            EXIST = .FALSE.
            LOW = STRIDE * TDEC
            STRIDE = MIN (HIGH, STRIDE * TINC)
            IF (1 .LT. LEVELM .AND. 0 .LT. TEXT)
     +         WRITE (TEXT, 20005) ID, STEP, YWORD, LOG10 (STRIDE)
            GO TO 1050
         END IF
         GO TO 2010
      END IF

C///  AGE := AGE + 1, M := M + 1.

      AGE = AGE + 1
      STEP = STEP + 1

C     RETAIN THE LATEST SOLUTION FOR USE BY THE FUNCTION.
C     GO TO 1080 WHEN ROUTE = 4
      ROUTE = 4
      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, V0, BUFFER)
      SIGNAL = 'RETAIN'
      GO TO 99999
1080  CONTINUE
      SIGNAL = ' '

C///  PRINT.

      IF (.NOT. (0 .LT. LEVELM .AND. 0 .LT. TEXT)) GO TO 1100
         CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, V0, BUFFER)
         SIGNAL = 'RESIDUAL'
         TIME = .FALSE.
C        GO TO 1090 WHEN ROUTE = 5
         ROUTE = 5
         GO TO 99999
1090     CONTINUE
         SIGNAL = ' '
         CALL TWNORM (GROUPA + COMPS * POINTS + GROUPB, YNORM, BUFFER)
         CALL TWLOGR (YWORD, YNORM)

         IF (1 .EQ. LEVELM) WRITE (TEXT, 10005)
     +      STEP, YWORD, CWORD, LOG10 (STRIDE), NUMBER, JWORD
1100  CONTINUE

C///  M < LAST?

      IF (STEP .LT. LAST) GO TO 1040

C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

2010  CONTINUE

C///  PRINT.

      IF (0 .LT. LEVELM .AND. 0 .LT. TEXT) THEN
         IF (1 .EQ. LEVELM) THEN
            IF (STEP .EQ. FIRST) THEN
               WRITE (TEXT, 10006) ID
            ELSE IF (STEP .EQ. LAST) THEN
               WRITE (TEXT, 10007) ID
            ELSE
               WRITE (TEXT, 10008) ID
            END IF
         ELSE IF (1 .LT. LEVELM) THEN
            IF (STEP .EQ. FIRST) THEN
               WRITE (TEXT, 10006) ID
            ELSE IF (STEP .EQ. LAST) THEN
               WRITE (TEXT, 20006) ID, STEP, YWORD
            ELSE
               WRITE (TEXT, 20007) ID, STEP, YWORD
            END IF
         END IF

         IF (FIRST .LT. LAST .AND. 1 .EQ. LEVELD) THEN
            WRITE (TEXT, 20008) ID
            CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, V0, BUFFER)
            SIGNAL = 'SHOW'
C           GO TO 2020 WHEN ROUTE = 6
            ROUTE = 6
            GO TO 99999
         END IF
      END IF

2020  CONTINUE
      SIGNAL = ' '

C///  SET THE COMPLETION STATUS FLAGS.

      SUCCES = FIRST .LT. STEP
      IF (STEP .LT. LAST) REPORT = XREPOR

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +  (/1X, A9, 'BEGIN TIME EVOLUTION.'
     +  /3(/10X, A44, A21)
     +  /10X, I6, 3X, A6)

10002 FORMAT
     +  (/1X, A9, 'CONTINUE TIME EVOLUTION.'
     +  /3(/10X, A44, A21)
     +  /10X, I6, 3X, A6)

10003 FORMAT
     +  (10X, I6, 21X, F6.2, 3X, I5, 3X, A12, 3X, A)

10004 FORMAT
     +  (10X, I6, 12X, A6, 3X, F6.2, 3X, I5, 3X, A12)

10005 FORMAT
     +  (10X, I6, 2(3X, A6), 3X, F6.2, 3X, I5, 3X, A12)

10006 FORMAT
     +  (/1X, A9, 'FAILURE.  NO TIME EVOLUTION.')

10007 FORMAT
     +  (/1X, A9, 'SUCCESS.  TIME EVOLUTION COMPLETED.')

10008 FORMAT
     +  (/1X, A9, 'PARTIAL SUCCESS.  TIME EVOLUTION INCOMPLETE.')

20001 FORMAT
     +  (/1X, A9, 'BEGIN TIME EVOLUTION.'
     + //10X, I10, '  LATEST TIME POINT'
     +  /14X, A6, '  LOG10 STEADY STATE RESIDUAL HERE'
     +  /10X, F10.2, '  LOG10 STRIDE TO NEXT TIME POINT'
     + //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')

20002 FORMAT
     +  (/1X, A9, 'CONTINUE TIME EVOLUTION.'
     + //10X, I10, '  LATEST TIME POINT'
     +  /14X, A6, '  LOG10 STEADY STATE RESIDUAL HERE'
     +  /10X, F10.2, '  LOG10 STRIDE TO NEXT TIME POINT'
     + //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')

20003 FORMAT
     +  (/1X, A9, 'CONTINUE TIME EVOLUTION WITH INCREASED STRIDE.'
     + //10X, I10, '  LATEST TIME POINT'
     +  /14X, A6, '  LOG10 STEADY STATE RESIDUAL HERE'
     +  /10X, F10.2, '  LOG10 INCREASED STRIDE TO NEXT TIME POINT'
     + //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')

20004 FORMAT
     +  (/1X, A9, 'RETRY THE STEP WITH A DECREASED TIME STRIDE.'
     + //10X, I10, '  LATEST TIME POINT'
     +  /14X, A6, '  LOG10 STEADY STATE RESIDUAL HERE'
     +  /10X, F10.2, '  LOG10 DECREASED STRIDE TO NEXT TIME POINT'
     + //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE, AGAIN.')

20005 FORMAT
     +  (/1X, A9, 'THE SOLUTION DID NOT CHANGE.  RETRYING THE STEP'
     +  /10X, 'WITH AN INCREASED TIME STRIDE.'
     + //10X, I10, '  LATEST TIME POINT'
     +  /14X, A6, '  LOG10 STEADY STATE RESIDUAL HERE'
     +  /10X, F10.2, '  LOG10 INCREASED STRIDE TO NEXT TIME POINT'
     + //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE, AGAIN.')

20006 FORMAT
     +  (/1X, A9, 'SUCCESS.  TIME EVOLUTION COMPLETED.'
     + //10X, I10, '  LAST TIME POINT'
     +  /14X, A6, '  LOG10 STEADY STATE RESIDUAL HERE')

20007 FORMAT
     +  (/1X, A9, 'PARTIAL SUCCESS.  TIME EVOLUTION INCOMPLETE.'
     + //10X, I10, '  LAST TIME POINT'
     +  /14X, A6, '  LOG10 STEADY STATE RESIDUAL HERE')

20008 FORMAT
     + (/1X, A9, 'THE LATEST SOLUTION:')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, ROUTE
      IF (.NOT. MESS) GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID,
     +   COMPS, POINTS, GROUPA, GROUPB, GROUPA + COMPS * POINTS + GROUPB
      IF (.NOT. MESS) GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, DESIRE
      IF (.NOT. MESS) GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, TDEC, TINC
      IF (.NOT. MESS) GO TO 99999

9005  IF (0 .LT. TEXT) WRITE (TEXT, 99005) ID, TMIN, TMAX
      IF (.NOT. MESS) GO TO 99999

9006  IF (0 .LT. TEXT) WRITE (TEXT, 99006) ID, TMIN, STRID0, TMAX
      IF (.NOT. MESS) GO TO 99999

9007  IF (0 .LT. TEXT) WRITE (TEXT, 99007) ID, STEP
      IF (.NOT. MESS) GO TO 99999

9008  IF (0 .LT. TEXT) WRITE (TEXT, 99008) ID, STEPS2
      IF (.NOT. MESS) GO TO 99999

9009  IF (0 .LT. TEXT) WRITE (TEXT, 99009) ID
      IF (.NOT. MESS) GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, I10, '  ROUTE')

99002 FORMAT
     +  (/1X, A9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES'
     +  /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  TOTAL UNKNOWNS')

99003 FORMAT
     +  (/1X, A9, 'ERROR.  THE NUMBER OF TIME STEPS MUST BE POSITIVE.'
     + //10X, I10, '  STEPS0 OR STEPS1, DESIRED NUMBER OF STEPS')

99004 FORMAT
     +  (/1X, A9, 'ERROR.  THE FACTORS FOR CHANGING THE TIME STRIDE'
     +  /10X, 'MUST BE NO SMALLER THAN 1.'
     + //10X, 1P, E10.2, '  TDEC, DECREASE FACTOR',
     +  /10X, 1P, E10.2, '  TINC, INCREASE FACTOR')

99005 FORMAT
     +  (/1X, A9, 'ERROR.  THE BOUNDS ON THE TIME STRIDE ARE OUT OF'
     +  /10X, 'ORDER.'
     + //10X, 1P, E10.2, '  TMIN, SHORTEST STRIDE'
     +  /10X, 1P, E10.2, '  TMAX, LONGEST STRIDE')

99006 FORMAT
     +  (/1X, A9, 'ERROR.  THE INITIAL TIME STRIDE MUST LIE BETWEEN'
     +  /10X, 'THE LOWER AND UPPER BOUNDS.'
     + //10X, 1P, E10.2, '  TMIN, SHORTEST STRIDE'
     +  /10X, 1P, E10.2, '  STRID0, INITIAL STRIDE'
     +  /10X, 1P, E10.2, '  TMAX, LONGEST STRIDE')

99007 FORMAT
     +  (/1X, A9, 'ERROR.  THE COUNT OF TIME STEPS MUST BE ZERO OR'
     +  /10X, 'POSITIVE.'
     + //10X, I10, '  STEP')

99008 FORMAT
     +  (/1X, A9, 'ERROR.  THE TIME STEPS BEFORE STRIDE INCREASES'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, I10, '  STEPS2, TIME STEPS BEFORE STRIDE INCREASES')

99009 FORMAT
     +  (/1X, A9, 'ERROR.  SEARCH FAILS.')

C///  EXIT.

      STOP
99999 CONTINUE
      RETURN
      END
      SUBROUTINE SEARCH
     +  (ERROR, TEXT,
     +   ABOVE, AGE, BELOW, BUFFER, COMPS, CONDIT, EXIST, GROUPA,
     +   GROUPB, LEVELD, LEVELM, NAME, NAMES, POINTS, REPORT, S0, S1,
     +   SIGNAL, STEPS, SUCCES, V0, V1, XXABS, XXAGE, XXREL, Y0, Y0NORM,
     +   Y1)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     SEARCH
C
C     PERFORM THE DAMPED, MODIFIED NEWTON'S SEARCH.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)

      CHARACTER
     +   COLUMN*16, CTEMP1*80, CTEMP2*80, HEADER*80, ID*9, NAME*(*),
     +   SIGNAL*(*), STRING*80
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   ABOVE, ABS0, ABS1, BELOW, BUFFER, CONDIT, DELTAB, DELTAD, REL0,
     +   REL1, S0, S0NORM, S1, S1NORM, SJ, TEMP, V0, V1, VALUE, VJ,
     +   XXABS, XXREL, Y0, Y0NORM, Y1, Y1NORM, ZERO
      EXTERNAL
     +   TWCOPY, TWLOGR, TWNORM, TWSQEZ
      INTEGER
     +   AGE, COMPS, COUNT, ENTRY, EXPONE, GROUPA, GROUPB, I, J, K,
     +   LEN1, LEN2, LENGTH, LEVELD, LEVELM, LINES, NAMES, NUMBER,
     +   POINTS, QBNDS, QDVRG, QNULL, REPORT, ROUTE, STEPS, TEXT, XXAGE
      INTRINSIC
     +   ABS, INT, LOG10, MAX, MIN, MOD
      LOGICAL
     +   ERROR, EXIST, FORCE, MESS, SUCCES

      PARAMETER (ID = 'SEARCH:  ')
      PARAMETER (LINES = 20)
      PARAMETER (ZERO = 0.0)

C     REPORT CODES
      PARAMETER (QNULL = 0, QBNDS = 1, QDVRG = 2)

      DIMENSION
     +   ABOVE(GROUPA + COMPS * POINTS + GROUPB),
     +   BELOW(GROUPA + COMPS * POINTS + GROUPB),
     +   BUFFER(GROUPA + COMPS * POINTS + GROUPB), COLUMN(7),
     +   HEADER(3, 2), NAME(NAMES),
     +   S0(GROUPA + COMPS * POINTS + GROUPB),
     +   S1(GROUPA + COMPS * POINTS + GROUPB),
     +   V0(GROUPA + COMPS * POINTS + GROUPB),
     +   V1(GROUPA + COMPS * POINTS + GROUPB),
     +   Y0(GROUPA + COMPS * POINTS + GROUPB),
     +   Y1(GROUPA + COMPS * POINTS + GROUPB)

C///  SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.

      SAVE

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  EVERY-TIME INITIALIZATION.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      MESS = .FALSE.

C     TURN OFF ALL COMPLETION STATUS FLAGS.
      ERROR = .FALSE.
      REPORT = QNULL
      SUCCES = .FALSE.

C///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      IF (SIGNAL .NE. ' ') THEN
         GO TO (2020, 2040, 2050, 2140, 2150, 2180) ROUTE
         ERROR = .TRUE.
         GO TO 9001
      END IF

C///  ONE-TIME INITIALIZATION.

      NUMBER = 0

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (((0 .LT. COMPS) .EQV. (0 .LT. POINTS)) .AND.
     +   0 .LE. COMPS .AND. 0 .LE. POINTS .AND. 0 .LE. GROUPA .AND.
     +   0 .LE. GROUPB .AND. 0 .LT. GROUPA + COMPS * POINTS + GROUPB)
      IF (ERROR) GO TO 9002

      ERROR = .NOT. (NAMES .EQ. 1 .OR.
     +   NAMES .EQ. GROUPA + COMPS + GROUPB)
      IF (ERROR) GO TO 9003

      COUNT = 0
      DO 1010 J = 1, GROUPA + COMPS * POINTS + GROUPB
         IF (.NOT. (BELOW(J) .LT. ABOVE(J))) COUNT = COUNT + 1
1010  CONTINUE
      ERROR = COUNT .NE. 0
      IF (ERROR) GO TO 9004

      COUNT = 0
      DO 1020 J = 1, GROUPA + COMPS * POINTS + GROUPB
         IF (.NOT. (BELOW(J) .LE. V0(J) .AND. V0(J) .LE. ABOVE(J)))
     +      COUNT = COUNT + 1
1020  CONTINUE
      ERROR = COUNT .NE. 0
      IF (ERROR) GO TO 9005

      ERROR = .NOT. (0.0 .LE. XXABS .AND. 0.0 .LE. XXREL)
      IF (ERROR) GO TO 9006

      ERROR = .NOT. (0 .LT. XXAGE)
      IF (ERROR) GO TO 9007

C///  WRITE ALL MESSAGES.

      IF (MESS .AND. 0 .LT. TEXT) THEN
         ROUTE = 0

         WRITE (TEXT, 10003) ID
         WRITE (TEXT, 10002) ID
         COUNT = 0
         DO 1030 J = 1, GROUPA + COMPS * POINTS + GROUPB
            COUNT = COUNT + 1
            IF (COUNT .LE. LINES) THEN
               IF (J .LE. GROUPA) THEN
                  I = J
               ELSE IF (J .LE. GROUPA + COMPS * POINTS) THEN
                  I = GROUPA + MOD (J - GROUPA - 1, COMPS) + 1
               ELSE
                  I = J - GROUPA - COMPS * POINTS
               END IF

               IF (NAMES .EQ. COMPS + GROUPA + GROUPB) THEN
                  CTEMP1 = NAME(I)
               ELSE
                  CTEMP1 = ' '
               END IF
               CALL TWSQEZ (LEN1, CTEMP1)

               IF (J .LE. GROUPA) THEN
                  WRITE (CTEMP2, 80001) 'A', I
               ELSE IF (J .LE. GROUPA + COMPS * POINTS) THEN
                  WRITE (CTEMP2, 80002) 'C', I,
     +               'P', INT ((J - GROUPA - 1) / COMPS) + 1
               ELSE
                  WRITE (CTEMP2, 80001) 'B', I
               END IF
               CALL TWSQEZ (LEN2, CTEMP2)

               IF (CTEMP1 .EQ. ' ') THEN
                  STRING = CTEMP2
                  LENGTH = LEN2
               ELSE IF (LEN1 + 2 + LEN2 .LE. 30) THEN
                  STRING = CTEMP1 (1 : LEN1) // '  ' // CTEMP2
                  LENGTH = LEN1 + 2 + LEN2
               ELSE IF (LEN1 + 1 + LEN2 .LE. 30) THEN
                  STRING = CTEMP1 (1 : LEN1) // ' ' // CTEMP2
                  LENGTH = LEN1 + 1 + LEN2
               ELSE
                  LEN1 = 30 - LEN2 - 4
                  STRING = CTEMP1 (1 : LEN1) // '... ' // CTEMP2
                  LENGTH = 30
               END IF

               WRITE (TEXT, 80003) 'LOWER', V0(J), STRING (1 : LENGTH)
            END IF
1030     CONTINUE
         IF (LINES .LT. COUNT) WRITE (TEXT, 80004)
         WRITE (TEXT, 10001) ID
         WRITE (TEXT, 10006) ID
         WRITE (TEXT, 10005) ID

         GO TO 9001
      END IF

C///  PRINT THE HEADER.

C                     123456789_123456789_123456789_123456789_1234
C                     123456   123456   123456   123456   123456
      HEADER(1, 1) = '         LOG10                              '
      HEADER(2, 1) = '  SLTN   -----------------------------------'
      HEADER(3, 1) = 'NUMBER   NORM F   COND J   NORM S      ABS A'

C                     123456789_123456789_123
C                     123456   123456  123456
      HEADER(1, 2) = '                       '
      HEADER(2, 2) = '-----------------------'
      HEADER(3, 2) = 'ND REL    DELTA B AND D'

      IF (LEVELM .GE. 1 .OR. MESS) THEN
         IF (0 .LT. TEXT) WRITE (TEXT, 10001)
     +      ID, ((HEADER(J, K), K = 1, 2), J = 1, 3)
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     SIR ISSAC NEWTON'S ALGORITHM.
C
C///////////////////////////////////////////////////////////////////////

C///  J EXIST?

      IF (.NOT. EXIST) GO TO 2010

C///  AGE < XXAGE?

      IF (AGE .LT. XXAGE) GO TO 2030

C///  EVALUATE J AT V0.  RE-EVALUATE Y0 := F(V0) IN CASE F CHANGES WHEN
C///  J DOES.  SOLVE J S0 = Y0.  EVAUATE ABS0 AND REL0.

2010  CONTINUE

      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, V0, BUFFER)
      SIGNAL = 'PREPARE'
C     GO TO 2020 WHEN ROUTE = 1
      ROUTE = 1
      GO TO 99999
2020  CONTINUE
      SIGNAL = ' '
      AGE = 0

C     JACOBIAN EVALUATION SHOULD RETURN A NEW RESIDUAL TOO.

      IF (0 .LT. LEVELM .AND. 0 .LT. TEXT) THEN
         IF (0.0 .LT. CONDIT) THEN
            WRITE (COLUMN(2), '(F6.2)') LOG10 (CONDIT)
         ELSE
            COLUMN(2) = '    NA'
         END IF
      END IF

C///  EVALUATE Y0 := F(V0).  SOLVE J S0 = Y0.  EVAUATE ABS0 AND REL0.

2030  CONTINUE

      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, V0, BUFFER)
      SIGNAL = 'RESIDUAL'
C     GO TO 2040 WHEN ROUTE = 2
      ROUTE = 2
      GO TO 99999
2040  CONTINUE
      SIGNAL = ' '
      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, BUFFER, Y0)
      CALL TWNORM (GROUPA + COMPS * POINTS + GROUPB, Y0NORM, Y0)

      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, Y0, BUFFER)
      SIGNAL = 'SOLVE'
C     GO TO 2050 WHEN ROUTE = 3
      ROUTE = 3
      GO TO 99999
2050  CONTINUE
      SIGNAL = ' '
      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, BUFFER, S0)
      CALL TWNORM (GROUPA + COMPS * POINTS + GROUPB, S0NORM, S0)

      ABS0 = 0.0
      REL0 = 0.0
      DO 2060 J = 1, GROUPA + COMPS * POINTS + GROUPB
         SJ = ABS (V0(J) - (V0(J) - S0(J)))
         VJ = ABS (V0(J))
         IF (XXREL * VJ .LT. SJ) ABS0 = MAX (ABS0, SJ)
         IF (XXABS .LT. SJ .AND. 0.0 .LT. VJ)
     +      REL0 = MAX (REL0, SJ / VJ)
2060  CONTINUE

C///  CHECK FOR SUCCESS.

      IF (ABS0 .LE. XXABS .AND. REL0 .LE. XXREL) GO TO 2170

C///  CHOOSE DELTAB.

2070  CONTINUE

C     DELTAB IS THE LARGEST DAMPING COEFFICIENT BETWEEN 0 AND 1 THAT
C     KEEPS V1 WITHIN BOUNDS.  IF V1 BELONGS ON THE BOUNDARY, THEN
C     PROVISIONS ARE MADE TO FORCE IT THERE DESPITE ROUNDING ERROR.

      DELTAB = 1.0
      FORCE = .FALSE.
      DO 2080 J = 1, GROUPA + COMPS * POINTS + GROUPB
         IF (S0(J) .GT. MAX (ZERO, V0(J) - BELOW(J))) THEN
            TEMP = (V0(J) - BELOW(J)) / S0(J)
            IF (TEMP .LT. DELTAB) THEN
               DELTAB = TEMP
               ENTRY = J
               FORCE = .TRUE.
               VALUE = BELOW(J)
            END IF
         ELSE IF (S0(J) .LT. MIN (ZERO, V0(J) - ABOVE(J))) THEN
            TEMP = (V0(J) - ABOVE(J)) / S0(J)
            IF (TEMP .LT. DELTAB) THEN
               DELTAB = TEMP
               ENTRY = J
               FORCE = .TRUE.
               VALUE = ABOVE(J)
            END IF
         END IF
2080  CONTINUE

      ERROR = DELTAB .LT. 0.0
      IF (ERROR) GO TO 9008

C///  0 < DELTAB?

      IF (.NOT. (0.0 .LT. DELTAB)) THEN
         IF (0 .LT. AGE) GO TO 2010

         IF (0 .LT. LEVELM .AND. 0 .LT. TEXT) THEN
            CALL TWLOGR (COLUMN(1), Y0NORM)
            CALL TWLOGR (COLUMN(3), S0NORM)
            CALL TWLOGR (COLUMN(4), ABS0)
            CALL TWLOGR (COLUMN(5), REL0)
            COLUMN(6) = ' '
            IF (DELTAB .NE. 1.0) CALL TWLOGR (COLUMN(6), DELTAB)
            COLUMN(7) = ' '
            IF (DELTAD .NE. 1.0) CALL TWLOGR (COLUMN(7), DELTAD)
            WRITE (TEXT, 10004) NUMBER, COLUMN
            WRITE (TEXT, 10002) ID

            COUNT = 0
            DO 2090 J = 1, GROUPA + COMPS * POINTS + GROUPB
               IF ((BELOW(J) .EQ. V0(J) .AND. 0.0 .LT. S0(J)) .OR.
     +            (V0(J) .EQ. ABOVE(J) .AND. S0(J) .LT. 0.0)) THEN
                  COUNT = COUNT + 1
                  IF (COUNT .LE. LINES) THEN
                     IF (J .LE. GROUPA) THEN
                        I = J
                     ELSE IF (J .LE. GROUPA + COMPS * POINTS) THEN
                        I = GROUPA + MOD (J - GROUPA - 1, COMPS) + 1
                     ELSE
                        I = J - GROUPA - COMPS * POINTS
                     END IF

                     IF (NAMES .EQ. COMPS + GROUPA + GROUPB) THEN
                        CTEMP1 = NAME(I)
                     ELSE
                        CTEMP1 = ' '
                     END IF
                     CALL TWSQEZ (LEN1, CTEMP1)

                     IF (J .LE. GROUPA) THEN
                        WRITE (CTEMP2, 80001) 'A', I
                     ELSE IF (J .LE. GROUPA + COMPS * POINTS) THEN
                        WRITE (CTEMP2, 80002) 'C', I,
     +                     'P', INT ((J - GROUPA - 1) / COMPS) + 1
                     ELSE
                        WRITE (CTEMP2, 80001) 'B', I
                     END IF
                     CALL TWSQEZ (LEN2, CTEMP2)

                     IF (CTEMP1 .EQ. ' ') THEN
                        STRING = CTEMP2
                        LENGTH = LEN2
                     ELSE IF (LEN1 + 2 + LEN2 .LE. 30) THEN
                        STRING = CTEMP1 (1 : LEN1) // '  ' // CTEMP2
                        LENGTH = LEN1 + 2 + LEN2
                     ELSE IF (LEN1 + 1 + LEN2 .LE. 30) THEN
                        STRING = CTEMP1 (1 : LEN1) // ' ' // CTEMP2
                        LENGTH = LEN1 + 1 + LEN2
                     ELSE
                        LEN1 = 30 - LEN2 - 4
                        STRING = CTEMP1 (1 : LEN1) // '... ' // CTEMP2
                        LENGTH = 30
                     END IF

                     IF (BELOW(J) .EQ. V0(J)) THEN
                        WRITE (TEXT, 80003)
     +                     'LOWER', V0(J), STRING (1 : LENGTH)
                     ELSE
                        WRITE (TEXT, 80003)
     +                     'UPPER', V0(J), STRING (1 : LENGTH)
                     END IF
                  END IF
               END IF
2090        CONTINUE
            IF (LINES .LT. COUNT) WRITE (TEXT, 80005)
         END IF

         REPORT = QBNDS
         SUCCES = .FALSE.
         GO TO 99999
      END IF

C///  DELTAD := 1.

      DELTAD = 1.0
      EXPONE = 0

C///  V1 := V0 - DELTAB DELTAD S0.  EVALUATE Y1 := F(V1).  SOLVE
C///  J S1 = Y1.  EVALUATE ABS1 AND REL1.

2100  CONTINUE

      TEMP = DELTAB * DELTAD
      DO 2110 J = 1, GROUPA + COMPS * POINTS + GROUPB
         V1(J) = V0(J) - TEMP * S0(J)
2110  CONTINUE

C     KEEP V1 IN BOUNDS DESPITE ROUNDING ERROR.

      DO 2120 J = 1, GROUPA + COMPS * POINTS + GROUPB
         V1(J) = MIN (V1(J), ABOVE(J))
2120  CONTINUE
      DO 2130 J = 1, GROUPA + COMPS * POINTS + GROUPB
         V1(J) = MAX (V1(J), BELOW(J))
2130  CONTINUE
      IF (EXPONE .EQ. 0 .AND. FORCE) V1(ENTRY) = VALUE

      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, V1, BUFFER)
      SIGNAL = 'RESIDUAL'
C     GO TO 2140 WHEN ROUTE = 4
      ROUTE = 4
      GO TO 99999
2140  CONTINUE
      SIGNAL = ' '
      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, BUFFER, Y1)
      CALL TWNORM (GROUPA + COMPS * POINTS + GROUPB, Y1NORM, Y1)

      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, Y1, BUFFER)
      SIGNAL = 'SOLVE'
C     GO TO 2150 WHEN ROUTE = 5
      ROUTE = 5
      GO TO 99999
2150  CONTINUE
      SIGNAL = ' '
      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, BUFFER, S1)
      CALL TWNORM (GROUPA + COMPS * POINTS + GROUPB, S1NORM, S1)

      ABS1 = 0.0
      REL1 = 0.0
      DO 2160 J = 1, GROUPA + COMPS * POINTS + GROUPB
         SJ = ABS (V1(J) - (V1(J) - S1(J)))
         VJ = ABS (V1(J))
         IF (XXREL * VJ .LT. SJ) ABS1 = MAX (ABS1, SJ)
         IF (XXABS .LT. SJ .AND. 0.0 .LT. VJ)
     +      REL1 = MAX (REL1, SJ / VJ)
2160  CONTINUE

C///  NORM S1 < OR = NORM S0?

      IF (S1NORM .LE. S0NORM) THEN
      ELSE
         DELTAD = 0.5 * DELTAD
         EXPONE = EXPONE + 1
         IF (EXPONE .LE. 5) GO TO 2100
            IF (0 .LT. AGE) GO TO 2010
               IF (0 .LT. LEVELM .AND. 0 .LT. TEXT) THEN
                  CALL TWLOGR (COLUMN(1), Y0NORM)
                  CALL TWLOGR (COLUMN(3), S0NORM)
                  CALL TWLOGR (COLUMN(4), ABS0)
                  CALL TWLOGR (COLUMN(5), REL0)
                  COLUMN(6) = ' '
                  IF (DELTAB .NE. 1.0) CALL TWLOGR (COLUMN(6), DELTAB)
                  COLUMN(7) = ' '
                  IF (DELTAD .NE. 1.0) CALL TWLOGR (COLUMN(7), DELTAD)
                  WRITE (TEXT, 10004) NUMBER, COLUMN
                  WRITE (TEXT, 10003) ID
               END IF
               REPORT = QDVRG
               SUCCES = .FALSE.
               GO TO 99999
      END IF

C///  PRINT.

      IF (0 .LT. LEVELM .AND. 0 .LT. TEXT) THEN
         CALL TWLOGR (COLUMN(1), Y0NORM)
         CALL TWLOGR (COLUMN(3), S0NORM)
         CALL TWLOGR (COLUMN(4), ABS0)
         CALL TWLOGR (COLUMN(5), REL0)
         COLUMN(6) = ' '
         IF (DELTAB .NE. 1.0) CALL TWLOGR (COLUMN(6), DELTAB)
         COLUMN(7) = ' '
         IF (DELTAD .NE. 1.0) CALL TWLOGR (COLUMN(7), DELTAD)
         WRITE (TEXT, 10004) NUMBER, COLUMN
         COLUMN(2) = ' '
      END IF

C///  S0 := S1, U := V1, Y0 := Y1, AGE := AGE + 1.

      AGE = AGE + 1
      NUMBER = NUMBER + 1
      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, S1, S0)
      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, V1, V0)
      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, Y1, Y0)
      S0NORM = S1NORM
      Y0NORM = Y1NORM
      ABS0 = ABS1
      REL0 = REL1

C///  S0 SMALL VS V0?

      IF (.NOT. (ABS0 .LE. XXABS .AND. REL0 .LE. XXREL)) THEN
         IF (AGE .LT. XXAGE) GO TO 2070
         GO TO 2010
      END IF

C///  SUCCESS.

2170  CONTINUE

C///  PRINT.

      IF (0 .LT. LEVELM .AND. 0 .LT. TEXT) THEN
         CALL TWLOGR (COLUMN(1), Y0NORM)
         CALL TWLOGR (COLUMN(3), S0NORM)
         CALL TWLOGR (COLUMN(4), ABS0)
         CALL TWLOGR (COLUMN(5), REL0)
         COLUMN(6) = ' '
         COLUMN(7) = ' '
         IF (0 .LT. LEVELD) THEN
            WRITE (TEXT, 10004) NUMBER, COLUMN
            WRITE (TEXT, 10005) ID
            SIGNAL = 'SHOW'
            CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, V0, BUFFER)
C           GO TO 2180 WHEN ROUTE = 6
            ROUTE = 6
            GO TO 99999
         ELSE
            WRITE (TEXT, 10004) NUMBER, COLUMN
            WRITE (TEXT, 10006) ID
         END IF
      END IF

2180  CONTINUE
      SIGNAL = ' '

      SUCCES = .TRUE.

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +  (/1X, A9, 'SOLVE NONLINEAR, NONDIFFERENTIAL EQUATIONS.'
     +  /4(/10X, A44, A23)/)

10002 FORMAT
     + (/1X, A9, 'FAILURE.  THE SEARCH FOR THE FOLLOWING UNKNOWNS GOES'
     +  /10X, 'OUT OF BOUNDS.'
C              12345  123456789_
     + //10X, 'BOUND       VALUE   UNKNOWN'
     +  /)

10003 FORMAT
     +  (/1X, A9, 'FAILURE.  THE SEARCH DIVERGES.')

10004 FORMAT
     +  (10X, I6, 3(3X, A6), 2(3X, A6, 2X, A6))

10005 FORMAT
     +  (/1X, A9, 'SUCCESS.  THE SOLUTION:')

10006 FORMAT
     +  (/1X, A9, 'SUCCESS.')

80001 FORMAT
     +   ('(', A, ' ', I10, ')')

80002 FORMAT
     +  ('(', A, ' ', I10, ' ', A, ' ', I10, ')')

80003 FORMAT
     +  (10X, A5, 2X, 1P, E10.2, 3X, A)

80004 FORMAT
     +  (30X, '... MORE')

80005 FORMAT
     +  (10X, '  ... MORE')

80006 FORMAT
     +  (10X, 1P, E10.2, 2X, E10.2, 3X, A)

80007 FORMAT
     +  (10X, 1P, E10.2, 2X, E10.2, 2X, E10.2, 3X, A)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, ROUTE
      IF (.NOT. MESS) GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID,
     +   COMPS, POINTS, GROUPA, GROUPB, GROUPA + COMPS * POINTS + GROUPB
      IF (.NOT. MESS) GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID,
     +   NAMES, COMPS, GROUPA, GROUPB, GROUPA + COMPS + GROUPB
      IF (.NOT. MESS) GO TO 99999

9004  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99004) ID,
     +      GROUPA, GROUPB, COMPS, GROUPA + COMPS + GROUPB, COUNT
         COUNT = 0
         DO 8010 J = 1, GROUPA + COMPS + GROUPB
            IF (.NOT. (BELOW(J) .LT. ABOVE(J)) .OR. MESS) THEN
               COUNT = COUNT + 1
               IF (COUNT .LE. LINES) THEN
                  IF (NAMES .EQ. COMPS + GROUPA + GROUPB) THEN
                     CTEMP1 = NAME(J)
                  ELSE
                     CTEMP1 = ' '
                  END IF
                  CALL TWSQEZ (LEN1, CTEMP1)

                  IF (J .LE. GROUPA) THEN
                     WRITE (CTEMP2, 80001) 'A', J
                  ELSE IF (J .LE. GROUPA + COMPS) THEN
                     WRITE (CTEMP2, 80001) 'C', J - GROUPA
                  ELSE
                     WRITE (CTEMP2, 80001) 'B', J - GROUPA - COMPS
                  END IF
                  CALL TWSQEZ (LEN2, CTEMP2)

                  IF (CTEMP1 .EQ. ' ') THEN
                     STRING = CTEMP2
                     LENGTH = LEN2
                  ELSE IF (LEN1 + 2 + LEN2 .LE. 40) THEN
                     STRING = CTEMP1 (1 : LEN1) // '  ' // CTEMP2
                     LENGTH = LEN1 + 2 + LEN2
                  ELSE IF (LEN1 + 1 + LEN2 .LE. 40) THEN
                     STRING = CTEMP1 (1 : LEN1) // ' ' // CTEMP2
                     LENGTH = LEN1 + 1 + LEN2
                  ELSE
                     LEN1 = 40 - LEN2 - 4
                     STRING = CTEMP1 (1 : LEN1) // '... ' // CTEMP2
                     LENGTH = 40
                  END IF

                  WRITE (TEXT, 80006)
     +               BELOW(J), ABOVE(J), STRING (1 : LENGTH)
               END IF
            END IF
8010     CONTINUE
         IF (LINES .LT. COUNT) WRITE (TEXT, 80005)
      END IF
      IF (.NOT. MESS) GO TO 99999

9005  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99005) ID, GROUPA, GROUPB, COMPS, POINTS,
     +      GROUPA + COMPS * POINTS + GROUPB, COUNT
         COUNT = 0
         DO 8020 J = 1, GROUPA + COMPS * POINTS + GROUPB
            IF (.NOT. (BELOW(J) .LE. V0(J) .AND. V0(J) .LE. ABOVE(J))
     +         .OR. MESS) THEN
               COUNT = COUNT + 1
               IF (COUNT .LE. LINES) THEN
                  IF (J .LE. GROUPA) THEN
                     I = J
                  ELSE IF (J .LE. GROUPA + COMPS * POINTS) THEN
                     I = GROUPA + MOD (J - GROUPA - 1, COMPS) + 1
                  ELSE
                     I = J - GROUPA - COMPS * POINTS
                  END IF

                  IF (NAMES .EQ. COMPS + GROUPA + GROUPB) THEN
                     CTEMP1 = NAME(I)
                  ELSE
                     CTEMP1 = ' '
                  END IF
                  CALL TWSQEZ (LEN1, CTEMP1)

                  IF (J .LE. GROUPA) THEN
                     WRITE (CTEMP2, 80001) 'A', I
                  ELSE IF (J .LE. GROUPA + COMPS * POINTS) THEN
                     WRITE (CTEMP2, 80002) 'C', I,
     +                  'P', INT ((J - GROUPA - 1) / COMPS) + 1
                  ELSE
                     WRITE (CTEMP2, 80001) 'B', I
                  END IF
                  CALL TWSQEZ (LEN2, CTEMP2)

                  IF (CTEMP1 .EQ. ' ') THEN
                     STRING = CTEMP2
                     LENGTH = LEN2
                  ELSE IF (LEN1 + 2 + LEN2 .LE. 30) THEN
                     STRING = CTEMP1 (1 : LEN1) // '  ' // CTEMP2
                     LENGTH = LEN1 + 2 + LEN2
                  ELSE IF (LEN1 + 1 + LEN2 .LE. 30) THEN
                     STRING = CTEMP1 (1 : LEN1) // ' ' // CTEMP2
                     LENGTH = LEN1 + 1 + LEN2
                  ELSE
                     LEN1 = 30 - LEN2 - 4
                     STRING = CTEMP1 (1 : LEN1) // '... ' // CTEMP2
                     LENGTH = 30
                  END IF

                  WRITE (TEXT, 80007)
     +               BELOW(J), V0(J), ABOVE(J), STRING (1 : LENGTH)
               END IF
            END IF
8020     CONTINUE
         IF (LINES .LT. COUNT) WRITE (TEXT, 80005)
      END IF
      IF (.NOT. MESS) GO TO 99999

9006  IF (0 .LT. TEXT) WRITE (TEXT, 99006) ID, XXABS, XXREL
      IF (.NOT. MESS) GO TO 99999

9007  IF (0 .LT. TEXT) WRITE (TEXT, 99007) ID, XXAGE
      IF (.NOT. MESS) GO TO 99999

9008  IF (0 .LT. TEXT) WRITE (TEXT, 99008) ID, DELTAB
      IF (.NOT. MESS) GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, I10, '  ROUTE')

99002 FORMAT
     +  (/1X, A9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES'
     +  /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  TOTAL UNKNOWNS')

99003 FORMAT
     +  (/1X, A9, 'ERROR.  THE NUMBER OF NAMES IS WRONG.'
     + //10X, I10, '  NAMES'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  TOTAL NUMBER')

99004 FORMAT
     +  (/1X, A9, 'ERROR.  THE LOWER AND UPPER BOUNDS ON SOME UNKNOWNS'
     +  /10X, 'ARE OUT OF ORDER.'
     + //10X, I10, '  GROUP A UNKNOWNS (A)'
     +  /10X, I10, '  GROUP B UNKNOWNS (B)'
     +  /10X, I10, '  COMPONENTS AT POINTS (C)'
     +  /10X, I10, '  TOTAL TYPES OF UNKNOWNS'
     +  /10X, I10, '  NUMBER OF BOUNDS OUT OF ORDER'
C              123456789_  123456789_
     + //10X, '     LOWER       UPPER'
     +  /10X, '     BOUND       BOUND   UNKNOWN'
     +  /)

99005 FORMAT
     +  (/1X, A9, 'ERROR.  THE GUESSES FOR SOME UNKNOWNS ARE OUT OF'
     +  /10X, 'BOUNDS.'
     + //10X, I10, '  GROUP A UNKNOWNS (A)'
     +  /10X, I10, '  GROUP B UNKNOWNS (B)'
     +  /10X, I10, '  COMPONENTS AT POINTS (C)'
     +  /10X, I10, '  POINTS (P)'
     +  /10X, I10, '  TOTAL UNKNOWNS'
     +  /10X, I10, '  NUMBER OUT OF BOUNDS'
C              123456789_  123456789_  123456789_
     + //10X, '     LOWER                   UPPER'
     +  /10X, '     BOUND       VALUE       BOUND   UNKNOWN'
     +  /)

99006 FORMAT
     +  (/1X, A9, 'ERROR.  THE BOUNDS FOR THE ABSOLUTE AND RELATIVE'
     +  /10X, 'CONVERGENCE TESTS MUST BE ZERO OR POSITIVE.'
     + //10X, 1P, E10.2, '  SSABS OR TDABS, ABSOLUTE ERROR'
     +  /10X, 1P, E10.2, '  SSREL OR TDREL, RELATIVE ERROR')

99007 FORMAT
     +  (/1X, A9, 'ERROR.  THE RETIREMENT AGE OF THE JACOBIAN MATRIX'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, I10, '  SSAGE OR TDAGE, MATRIX RETIREMENT AGE')

99008 FORMAT
     +  (/1X, A9, 'ERROR.  THE DAMPING COEFFICIENT FOR STAYING'
     +  /10X, 'IN BOUNDS IS NEGATIVE.'
     + //10X, 1P, E10.2, '  DELTA B')

C///  EXIT.

      STOP
99999 CONTINUE

C     COPY THE PROTECTED LOCAL VARIABLE
      STEPS = NUMBER

      RETURN
      END
      SUBROUTINE TWCOPY (N, X, Y)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWCOPY
C
C     COPY ONE VECTOR TO ANOTHER.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      INTEGER J, N
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   X, Y

      DIMENSION X(N), Y(N)

      DO 0100 J = 1, N
         Y(J) = X(J)
0100  CONTINUE

      RETURN
      END
      SUBROUTINE TWEPS (EPS)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWEPS
C
C     FIND MACHINE EPSILON.
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
C*****MACHINE EPSILON > COMPUTED
C      LOGICAL
C     +   SAME
C*****end MACHINE EPSILON > COMPUTED

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
      SUBROUTINE TWGBCO (A, LDA, N, LOWER, UPPER, PIVOT, RCOND, Z)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWGBCO
C
C     FACTOR A BANDED MATRIX AND ESTIMATE THE RECIPROCAL OF ITS
C     CONDITION NUMBER.  BASED ON _GBCO FROM THE LINPACK LIBRARY.
C
C///////////////////////////////////////////////////////////////////////

      DOUBLE PRECISION
     +   DSUM
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, ANORM, EK, RCOND, S, SM, SUM, T, WK, WKM, YNORM, Z
      EXTERNAL
     +   TWGBFA
      INTEGER
     +   FIRST, INFO, J, JDIAG, JU, K, LAST, LDA, LOWER, MM, N, PIVOT,
     +   UPPER
      INTRINSIC
     +   ABS, DBLE, MAX, MIN, SIGN

      DIMENSION
     +   A(LDA,N), PIVOT(N), Z(N)

      JDIAG = LOWER + UPPER + 1

C///  COMPUTE THE 1-NORM OF A

      ANORM = 0.0
      DO 1020 K = 1, N
         FIRST = MAX (LOWER + 1, JDIAG + 1 - K)
         LAST = MIN (JDIAG + LOWER, JDIAG + N - K)
         SUM = 0.0
         DO 1010 J = FIRST, LAST
            SUM = SUM + ABS (A(J, K))
1010     CONTINUE
         ANORM = MAX (ANORM, SUM)
1020  CONTINUE

C///  FACTOR A

      CALL TWGBFA (A, LDA, N, LOWER, UPPER, PIVOT, INFO)

C///  SOLVE TRANSPOSE(U) * W = E

      EK = 1.0
      DO 2010 J = 1, N
         Z(J) = 0.0
2010  CONTINUE

      JU = 0
      DO 2050 K = 1, N
         IF (Z(K) .NE. 0.0) EK = SIGN (EK, - Z(K))

         IF (ABS (EK - Z(K)) .GT. ABS (A(JDIAG, K))) THEN
            S = ABS (A(JDIAG, K)) / ABS (EK - Z(K))
            EK = S * EK
            DO 2020 J = 1, N
               Z(J) = S * Z(J)
2020        CONTINUE
         END IF

         WK = EK - Z(K)
         WKM = - EK - Z(K)
         S = ABS (WK)
         SM = ABS (WKM)
         IF (A(JDIAG, K) .NE. 0.0) THEN
            WK = WK / A(JDIAG, K)
            WKM = WKM / A(JDIAG, K)
         ELSE
            WK = 1.0
            WKM = 1.0
         END IF

         JU = MIN (MAX (JU, UPPER + PIVOT(K)), N)
         MM = JDIAG
         IF (K + 1 .LE. JU) THEN
            DO 2030 J = K + 1, JU
               MM = MM - 1
               SM = SM + ABS (Z(J) + WKM * A(MM, J))
               Z(J) = Z(J) + WK * A(MM, J)
               S = S + ABS (Z(J))
2030        CONTINUE

            IF (S .LT. SM) THEN
               T = WKM - WK
               WK = WKM
               MM = JDIAG
               DO 2040 J = K + 1, JU
                  MM = MM - 1
                  Z(J) = Z(J) + T * A(MM, J)
2040           CONTINUE
            END IF
         END IF

         Z(K) = WK
2050  CONTINUE

      SUM = 0.0
      DO 2060 J = 1, N
         SUM = SUM + ABS (Z(J))
2060  CONTINUE
      S = 1.0 / SUM

      DO 2070 J = 1, N
         Z(J) = S * Z(J)
2070  CONTINUE

C///  SOLVE TRANSPOSE(L) * Y = W

      DO 3030 K = N, 1, - 1
         DSUM = 0.0
         DO 3010 J = 1, MIN (LOWER, N - K)
            DSUM = DSUM + DBLE (A(JDIAG + J, K)) * DBLE (Z(K + J))
3010     CONTINUE
         Z(K) = Z(K) + DSUM

         IF (1.0 .LT. ABS (Z(K))) THEN
            S = 1.0 / ABS (Z(K))
            DO 3020 J = 1, N
               Z(J) = S * Z(J)
3020        CONTINUE
         END IF

         J = PIVOT(K)
         T = Z(J)
         Z(J) = Z(K)
         Z(K) = T
3030  CONTINUE

      SUM = 0.0
      DO 3040 J = 1, N
         SUM = SUM + ABS (Z(J))
3040  CONTINUE
      S = 1.0 / SUM

      DO 3050 J = 1, N
         Z(J) = S * Z(J)
3050  CONTINUE

      YNORM = 1.0

C///  SOLVE L * V = Y

      DO 4030 K = 1, N
         J = PIVOT(K)
         T = Z(J)
         Z(J) = Z(K)
         Z(K) = T

         DO 4010 J = 1, MIN (LOWER, N - K)
            Z(K + J) = T * A(JDIAG + J, K) + Z(K + J)
4010     CONTINUE

         IF (1.0 .LT. ABS (Z(K))) THEN
            S = 1.0 / ABS (Z(K))
            DO 4020 J = 1, N
               Z(J) = S * Z(J)
4020        CONTINUE

            YNORM = S * YNORM
         END IF
4030  CONTINUE

      SUM = 0.0
      DO 4040 J = 1, N
         SUM = SUM + ABS (Z(J))
4040  CONTINUE
      S = 1.0 / SUM

      DO 4050 J = 1, N
         Z(J) = S * Z(J)
4050  CONTINUE

      YNORM = S * YNORM

C///  SOLVE U * Z = W

      DO 5030 K = N, 1, - 1
         IF (ABS (Z(K)) .GT. ABS (A(JDIAG, K))) THEN
            S = ABS (A(JDIAG, K)) / ABS (Z(K))
            DO 5010 J = 1, N
               Z(J) = S * Z(J)
5010        CONTINUE
            YNORM = S*YNORM
         END IF

         IF (A(JDIAG, K) .NE. 0.0) THEN
            Z(K) = Z(K) / A(JDIAG, K)
         ELSE
            Z(K) = 1.0
         END IF

         T = - Z(K)
         DO 5020 J = 1, MIN (K, JDIAG) - 1
            Z(K - J) = T * A(JDIAG - J, K) + Z(K - J)
5020     CONTINUE
5030  CONTINUE

      SUM = 0.0
      DO 5040 J = 1, N
         SUM = SUM + ABS (Z(J))
5040  CONTINUE
      S = 1.0 / SUM

      DO 5050 J = 1, N
         Z(J) = S * Z(J)
5050  CONTINUE

      YNORM = S * YNORM

C///  FORM RCOND

      IF (ANORM .NE. 0.0) THEN
         RCOND = YNORM / ANORM
      ELSE
         RCOND = 0.0
      END IF

C///  EXIT

      RETURN
      END
      SUBROUTINE TWGBFA (A, LDA, N, LOWER, UPPER, PIVOT, INFO)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWGBFA
C
C     FACTOR A BANDED MATRIX FOR TWGBCO. BASED ON _GBFA FROM THE LINPACK
C     LIBRARY.
C
C///////////////////////////////////////////////////////////////////////

C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, T, VALUE
      INTEGER
     +   I, INFO, J, JK, K, LDA, LOWER, N, PACK, PDIAG, PIVOT, PJK,
     +   SAVE, TRUE, UPPER
      INTRINSIC
     +   ABS, MAX, MIN

      DIMENSION
     +   A(LDA, N), PIVOT(N)

C///  STATEMENT FUNTIONS

C     PACKED ROW POSITION OF ENTRY J IN COLUMN K
      PACK (J, K) = J - K + PDIAG

C     TRUE ROW POSITION OF ENTRY J PACKED IN COLUMN K
      TRUE (J, K) = J - PDIAG + K

C///  INITIALIZE

C     PACKED ROW POSITION OF THE DIAGONAL
      PDIAG = LOWER + UPPER + 1

      INFO = 0

C///  TOP OF THE LOOP OVER COLUMNS

      DO 2060 K = 1, N

C///  INITIALIZE THE FILL-IN SPACE

         DO 2010 I = 1, LOWER
            A(I, K) = 0.0
2010     CONTINUE

C///  LOOP OVER THE PREVIOUS COLUMNS

         DO 2030 J = MAX (1, K - LOWER - UPPER), K - 1
            PJK = PACK (PIVOT(J), K)
            JK = PACK (J, K)
            T = A(PJK, K)
            IF (PJK .NE. JK) THEN
               A(PJK, K) = A(JK, K)
               A(JK, K) = T
            END IF

            IF (T .NE. 0.0) THEN
               DO 2020 I = 1, MIN (LOWER, N - J)
                  A(JK + I, K) = T * A(PDIAG + I, J) + A(JK + I, K)
2020           CONTINUE
            END IF
2030     CONTINUE

C///  FIND THE PIVOT

         SAVE = PDIAG
         VALUE = ABS (A(PDIAG, K))
         DO 2040 I = PDIAG + 1, PDIAG + MIN (LOWER, N - K)
            IF (VALUE .LT. ABS (A(I, K))) THEN
               SAVE = I
               VALUE = ABS (A(I, K))
            END IF
2040     CONTINUE
         PIVOT(K) = TRUE (SAVE, K)

C///  INTERCHANGE IF NECESSARY

         IF (SAVE .NE. PDIAG) THEN
            T = A(SAVE, K)
            A(SAVE, K) = A(PDIAG, K)
            A(PDIAG, K) = T
         END IF

C///  SCALE THE LOWER COLUMN

         IF (A(SAVE, K) .NE. 0.0) THEN
            T = - 1.0 / A(PDIAG, K)
            DO 2050 I = PDIAG + 1, PDIAG + MIN (LOWER, N - K)
               A(I, K) = T * A(I, K)
2050        CONTINUE
         ELSE
            INFO = K
         END IF

C///  BOTTOM OF THE LOOP OVER COLUMNS

2060  CONTINUE

C///  THE FINAL COLUMN IS TRIVIAL

      PIVOT(N) = N
      IF (A(PDIAG, N) .EQ. 0.0) INFO = N

C///  EXIT

      RETURN
      END
      SUBROUTINE TWGBSL (ABD, LDA, N, LOWER, UPPER, PIVOT, B)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWGBSL
C
C     SOLVE A SYSTEM OF LINEAR EQUATIONS FOR TWSOLV.  BASED ON _GBSL
C     FROM THE LINPACK LIBRARY.
C
C///////////////////////////////////////////////////////////////////////

C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   ABD, B, T
      INTEGER
     +   J, JDIAG, K, L, LA, LB, LDA, LM, LOWER, N, PIVOT, UPPER
      INTRINSIC
     +   MIN

      DIMENSION
     +   ABD(LDA,*), B(*), PIVOT(*)

      JDIAG = UPPER + LOWER + 1

      IF (0 .LT. LOWER) THEN
         DO 1020 K = 1, N - 1
            L = PIVOT(K)
            T = B(L)
            IF (L .NE. K) THEN
               B(L) = B(K)
               B(K) = T
            END IF

            LM = MIN (LOWER, N - K)
            DO 1010 J = 1, LM
               B(K + J) = T * ABD(JDIAG + J, K) + B(K + J)
1010        CONTINUE
1020     CONTINUE
      END IF

      DO 1040 K = N, 1, - 1
         B(K) = B(K) / ABD (JDIAG, K)
         LM = MIN (K, JDIAG) - 1
         LA = JDIAG - LM - 1
         LB = K - LM - 1
         T = - B(K)
         DO 1030 J = 1, LM
            B(LB + J) = T * ABD(LA + J, K) + B(LB + J)
1030     CONTINUE
1040  CONTINUE

      RETURN
      END
      SUBROUTINE TWGRAB (ERROR, LAST, FIRST, NUMBER)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWGRAB
C
C     RESERVE SPACE IN AN ARRAY.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)

      INTEGER
     +   FIRST, LAST, NUMBER
      INTRINSIC
     +   MAX
      LOGICAL
     +   ERROR

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (0 .LE. LAST)
      IF (ERROR) GO TO 99999

      ERROR = .NOT. (0 .LE. NUMBER)
      IF (ERROR) GO TO 99999

C///  GRAB THE SPACE.

      FIRST = LAST + 1
      LAST = LAST + MAX (1, NUMBER)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWINIT (ERROR, TEXT, FORCE)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWINIT
C
C     INITIALIZE THE CONTROLS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   ID*9
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   RVALUE
      INTEGER
     +   CNTRLS, COUNT, IVALUE, TEXT
      LOGICAL
     +   ERROR, FIRST, FORCE, LVALUE, MESS

      PARAMETER (ID = 'TWINIT:  ')
      PARAMETER (CNTRLS = 22)

      DIMENSION IVALUE(CNTRLS), LVALUE(CNTRLS), RVALUE(CNTRLS)

      COMMON / TWCOMI / IVALUE
      COMMON / TWCOML / LVALUE
      COMMON / TWCOMR / RVALUE

C     THE GNU F77 COMPILER REQUIRES THE SAVE TO PRECEED THE DATA

      SAVE FIRST

      DATA FIRST / .TRUE. /

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      MESS = .FALSE.

      IF (MESS .AND. 0 .LT. TEXT) GO TO 9001

C///  TOP OF THE BLOCK TO SET THE CONTROLS.

      IF (FIRST .OR. FORCE) THEN
         FIRST = .FALSE.

C///  SET THE CONTROLS.

      COUNT = 0

C     ADAPT

      COUNT = COUNT + 1
      LVALUE(COUNT) = .FALSE.

C     LEVELD

      COUNT = COUNT + 1
      IVALUE(COUNT) = 1

C     LEVELM

      COUNT = COUNT + 1
      IVALUE(COUNT) = 1

C     PADD

      COUNT = COUNT + 1
      LVALUE(COUNT) = .FALSE.

C     SSABS

      COUNT = COUNT + 1
      RVALUE(COUNT) = 1.0E-9

C     SSAGE

      COUNT = COUNT + 1
      IVALUE(COUNT) = 10

C     SSREL

      COUNT = COUNT + 1
      RVALUE(COUNT) = 1.0E-6

C     STEADY

      COUNT = COUNT + 1
      LVALUE(COUNT) = .TRUE.

C     STEPS0

      COUNT = COUNT + 1
      IVALUE(COUNT) = 0

C     STEPS1

      COUNT = COUNT + 1
      IVALUE(COUNT) = 200

C     STEPS2

      COUNT = COUNT + 1
      IVALUE(COUNT) = 100

C     STRID0

      COUNT = COUNT + 1
      RVALUE(COUNT) = 1.0E-4

C     TDABS

      COUNT = COUNT + 1
      RVALUE(COUNT) = 1.0E-9

C     TDAGE

      COUNT = COUNT + 1
      IVALUE(COUNT) = 20

C     TDEC

      COUNT = COUNT + 1
      RVALUE(COUNT) = 3.1623

C     TDREL

      COUNT = COUNT + 1
      RVALUE(COUNT) = 1.0E-6

C     TINC

      COUNT = COUNT + 1
      RVALUE(COUNT) = 10.0

C     TMAX

      COUNT = COUNT + 1
      RVALUE(COUNT) = 1.0E-2

C     TMIN

      COUNT = COUNT + 1
      RVALUE(COUNT) = 1.0E-20

C     TOLER0

      COUNT = COUNT + 1
      RVALUE(COUNT) = 1.0E-9

C     TOLER1

      COUNT = COUNT + 1
      RVALUE(COUNT) = 0.2

C     TOLER2

      COUNT = COUNT + 1
      RVALUE(COUNT) = 0.2

C///  BOTTOM OF THE BLOCK TO SET THE CONTROLS.

         ERROR = .NOT. (COUNT .EQ. CNTRLS)
         IF (ERROR) GO TO 9001
      END IF

C///  ERROR MESSAGES.

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, CNTRLS, COUNT
      IF (.NOT. MESS) GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.'
     + //10X, I10, '  CONTROLS'
     +  /10X, I10, '  COUNTED')

C///  EXIT.

      STOP
99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWLAPS (TIMER)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWLAPS
C
C     OBTAIN ELAPSED COMPUTING TIME IN SECONDS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      EXTERNAL TWTIME
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   TEMP, TIMER

      CALL TWTIME (TEMP)
      TIMER = TEMP - TIMER

      RETURN
      END
      SUBROUTINE TWLAST (LENGTH, STRING)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWLAST
C
C     FIND THE LAST NONBLANK CHARACTER IN A STRING.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   STRING*(*)
      INTEGER
     +   J, LENGTH
      INTRINSIC
     +   LEN

      DO 0100 J = LEN (STRING), 1, - 1
         IF (STRING(J : J) .NE. ' ') THEN
            LENGTH = J
            GO TO 0200
         END IF
0100  CONTINUE
      LENGTH = 1
0200  CONTINUE

      RETURN
      END
      SUBROUTINE TWLOGR (STRING, VALUE)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWLOGR
C
C     WRITE A COMMON LOGARITHM TO A CHARACTER STRING.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)

      CHARACTER
     +   STRING*(*)
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   VALUE
      INTEGER
     +   LEN
      INTRINSIC
     +   LEN, LOG10

      IF (6 .LE. LEN (STRING)) THEN
         IF (VALUE .LT. 0.0) THEN
            STRING = ' '
         ELSE IF (VALUE .EQ. 0.0) THEN
            STRING = '  ZERO'
         ELSE
            WRITE (STRING, '(F6.2)') LOG10 (VALUE)
         END IF
      ELSE
         STRING = '******'
      END IF

      RETURN
      END
      SUBROUTINE TWNORM (N, VALUE, X)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWNORM
C
C     COMPUTE THE MAX-NORM OF A VECTOR.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)

C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   VALUE, X
      INTEGER
     +   J, N
      INTRINSIC
     +   ABS, MAX

      DIMENSION X(N)

      VALUE = 0.0
      DO 0100 J = 1, N
         VALUE = MAX (VALUE, ABS (X(J)))
0100  CONTINUE

      RETURN
      END
      SUBROUTINE TWOPNT
     +  (ERROR, TEXT, VERSIO,
     +   ABOVE, ACTIVE, BELOW, BUFFER, COMPS, CONDIT, GROUPA, GROUPB,
     +   ISIZE, IWORK, MARK, NAME, NAMES, PMAX, POINTS, REPORT, RSIZE,
     +   RWORK, SIGNAL, STRIDE, TIME, U, X)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWOPNT
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   COLUMN*80, CTEMP1*80, CTEMP2*80, HEADER*80, ID*9, NAME*(*),
     +   REPORT*(*), SIGNAL*(*), STRING*80, VERSIO*(*), VNMBR*8
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   ABOVE, BELOW, BUFFER, CONDIT, DETAIL, MAXCON, RATIO, RVALUE,
     +   RWORK, SSABS, SSREL, STRID0, STRIDE, TDABS, TDEC, TDREL, TEMP,
     +   TIMER, TINC, TMAX, TMIN, TOLER0, TOLER1, TOLER2, TOTAL, U, X,
     +   YNORM
      EXTERNAL
     +   EVOLVE, REFINE, SEARCH, TWCOPY, TWGRAB, TWLAPS, TWLAST, TWLOGR,
     +   TWNORM, TWSQEZ, TWTIME, TWINIT
      INTEGER
     +   AGE, CNTRLS, COMPS, COUNT, DESIRE, EVENT, GMAX, GRID, GROUPA,
     +   GROUPB, ILAST, ISIZE, IVALUE, IWORK, J, JACOBS, K, LABEL, LEN1,
     +   LEN2, LENGTH, LEVELD, LEVELM, LINES, NAMES, NSTEPS, PADD, PMAX,
     +   POINTS, PSAVE, QABOVE, QBELOW, QBNDS, QDVRG, QENTRY, QEXIT,
     +   QFUNCT, QGRID, QJACOB, QNULL, QOTHER, QRAT1, QRAT2, QREFIN,
     +   QS0, QS1, QSEARC, QSOLVE, QTASK, QTIMST, QTOTAL, QTYPE, QUSAVE,
     +   QV1, QVARY, QVARY1, QVARY2, QVSAVE, QXSAVE, QY0, QY1, RETURN,
     +   RLAST, ROUTE, RSIZE, SIZE, SSAGE, STEP, STEPS, STEPS0, STEPS1,
     +   STEPS2, TDAGE, TEXT, VNMBRS, XREPOR
      INTRINSIC
     +   MAX
      LOGICAL
     +   ACTIVE, ADAPT, ALLOW, ERROR, EXIST, FIRST, FLAG, FOUND, LVALUE,
     +   MARK, MESS, SATISF, STEADY, TIME

      PARAMETER (ID = 'TWOPNT:  ')
      PARAMETER (CNTRLS = 22)
      PARAMETER (GMAX = 100)
      PARAMETER (LINES = 20)
      PARAMETER (VNMBRS = 12)

C     REPORT CODES
      PARAMETER (QNULL = 0, QBNDS = 1, QDVRG = 2)

C     LOCATION OF DATA IN ARRAYS DETAIL, EVENT, TIMER, AND TOTAL.  THE
C     LOCATIONS ARE CHOSEN TO SIMPLIFY WRITE STATEMENTS.  DETAIL USES
C     ONLY 1 : 8, EVENT USES ONLY 5 : 8, TIMER USES 1 : 9, AND TOTAL
C     USES ONLY 2 : 9.  IN ADDITION, 2, 3, 4, 10, AND 11 ARE USED AS
C     MNEMONIC VALUES FOR QTASK.
      PARAMETER
     +  (QGRID  =  1,
     +   QTIMST =  2,
     +   QSEARC =  3,
     +   QREFIN =  4,
     +   QFUNCT =  5,
     +   QJACOB =  6,
     +   QSOLVE =  7,
     +   QOTHER =  8,
     +   QTOTAL =  9,
     +   QENTRY = 10,
     +   QEXIT  = 11)

      DIMENSION
     +   ABOVE(GROUPA + COMPS + GROUPB), ACTIVE(*), BELOW(GROUPA + COMPS
     +   + GROUPB), BUFFER(GROUPA + COMPS * PMAX + GROUPB), COLUMN(3),
     +   DETAIL(GMAX, QTOTAL), EVENT(GMAX, QTOTAL), HEADER(6),
     +   IVALUE(CNTRLS), IWORK(ISIZE), LVALUE(CNTRLS), MARK(*),
     +   NAME(NAMES), RATIO(2), RVALUE(CNTRLS), RWORK(RSIZE),
     +   SIZE(GMAX), TIMER(QTOTAL), TOTAL(QTOTAL), U(GROUPA + COMPS *
     +   PMAX + GROUPB), VNMBR(VNMBRS), X(*)

      COMMON / TWCOMI / IVALUE
      COMMON / TWCOML / LVALUE
      COMMON / TWCOMR / RVALUE

C///  SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.

      SAVE

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  EVERY-TIME INITIALIZATION.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      MESS = .FALSE.

C     TURN OFF ALL REVERSE COMMUNICATION FLAGS.
      TIME = .FALSE.

C///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      IF (SIGNAL .NE. ' ') THEN
         GO TO (9912, 9922, 9932, 9942) ROUTE
         ERROR = .TRUE.
         GO TO 9001
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     ENTRY BLOCK.  INITIALIZE A NEW PROBLEM.
C
C///////////////////////////////////////////////////////////////////////

C///  TURN OFF ALL STATUS REPORTS.

      ERROR = .FALSE.
      REPORT = ' '

C///  WRITE ALL MESSAGES.

      IF (MESS .AND. 0 .LT. TEXT) THEN
         LABEL = 0
         RETURN = 0
         ROUTE = 0

         WRITE (TEXT, 10004) ID, '???'
         WRITE (TEXT, 10020) ID
         WRITE (TEXT, 10017) ID
         WRITE (TEXT, 10014) ID
         STRING = VNMBR(VNMBRS)
         CALL TWLAST (LENGTH, STRING)
         WRITE (TEXT, 10001) ID, 'DOUBLE PRECISION', STRING (1 : LENGTH)
         WRITE (TEXT, 10022) ID
         WRITE (TEXT, 10021) ID
         WRITE (TEXT, 10011) ID, '???', RATIO, TOLER1, TOLER2
         WRITE (TEXT, 10013) ID, '???', RATIO, TOLER1, TOLER2
         WRITE (TEXT, 10012) ID
         WRITE (TEXT, 10002) ID, 'FINAL SOLUTION:'
         WRITE (TEXT, 10002) ID, 'INITIAL GUESS:'
         WRITE (TEXT, 10019) ID
         WRITE (TEXT, 10018) ID
         WRITE (TEXT, 10016) ID
         WRITE (TEXT, 10015) ID
         STRING = VNMBR(VNMBRS)
         CALL TWLAST (LENGTH, STRING)
         WRITE (TEXT, 10001) ID, 'SINGLE PRECISION', STRING (1 : LENGTH)
         WRITE (TEXT, 10002) ID, 'SOLVE THE PROBLEM.'
         WRITE (TEXT, 10010) ID

         GO TO 9001
      END IF

C///  CHECK THE VERSION.

      DATA VNMBR
     +   / '3.18', '3.19', '3.20', '3.21', '3.22', '3.23', '3.24',
     +     '3.25', '3.26', '3.27', '3.28', '3.29' /

      FLAG = .FALSE.
      DO 1010 J = 1, VNMBRS
         FLAG = FLAG .OR.
C*****PRECISION > DOUBLE
     +      VERSIO .EQ. 'DOUBLE PRECISION VERSION ' // VNMBR(J)
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C     +      VERSIO .EQ. 'SINGLE PRECISION VERSION ' // VNMBR(J)
C*****END PRECISION > SINGLE
1010  CONTINUE
      ERROR = .NOT. FLAG
      IF (ERROR) GO TO 9002

C///  SET THE CONTROLS.

C     SUBROUTINE TWINIT (ERROR, TEXT, FORCE)

      CALL TWINIT (ERROR, TEXT, .FALSE.)
      IF (ERROR) GO TO 9003

      COUNT = 0

      COUNT = COUNT + 1
      ADAPT = LVALUE(COUNT)

      COUNT = COUNT + 1
      LEVELD = IVALUE(COUNT)

      COUNT = COUNT + 1
      LEVELM = IVALUE(COUNT)

      COUNT = COUNT + 1
      IF (LVALUE(COUNT)) THEN
         PADD = IVALUE(COUNT)
      ELSE
         PADD = PMAX
      END IF

      COUNT = COUNT + 1
      SSABS = RVALUE(COUNT)

      COUNT = COUNT + 1
      SSAGE = IVALUE(COUNT)

      COUNT = COUNT + 1
      SSREL = RVALUE(COUNT)

      COUNT = COUNT + 1
      STEADY = LVALUE(COUNT)

      COUNT = COUNT + 1
      STEPS0 = IVALUE(COUNT)

      COUNT = COUNT + 1
      STEPS1 = IVALUE(COUNT)

      COUNT = COUNT + 1
      STEPS2 = IVALUE(COUNT)

      COUNT = COUNT + 1
      STRID0 = RVALUE(COUNT)

      COUNT = COUNT + 1
      TDABS = RVALUE(COUNT)

      COUNT = COUNT + 1
      TDAGE = IVALUE(COUNT)

      COUNT = COUNT + 1
      TDEC = RVALUE(COUNT)

      COUNT = COUNT + 1
      TDREL = RVALUE(COUNT)

      COUNT = COUNT + 1
      TINC = RVALUE(COUNT)

      COUNT = COUNT + 1
      TMAX = RVALUE(COUNT)

      COUNT = COUNT + 1
      TMIN = RVALUE(COUNT)

      COUNT = COUNT + 1
      TOLER0 = RVALUE(COUNT)

      COUNT = COUNT + 1
      TOLER1 = RVALUE(COUNT)

      COUNT = COUNT + 1
      TOLER2 = RVALUE(COUNT)

      ERROR = .NOT. (COUNT .EQ. CNTRLS)
      IF (ERROR) GO TO 9004

C///  PRINT THE ENTRY BANNER AT ALL PRINT LEVELS.

      STRING = VNMBR(VNMBRS)
      CALL TWLAST (LENGTH, STRING)
      IF ((0 .LT. LEVELM .OR. MESS) .AND. 0 .LT. TEXT)
     +   WRITE (TEXT, 10001) ID,
C*****PRECISION > DOUBLE
     +   'DOUBLE PRECISION', STRING (1 : LENGTH)
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C     +   'SINGLE PRECISION', STRING (1 : LENGTH)
C*****END PRECISION > SINGLE

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (LEVELD .LE. LEVELM)
      IF (ERROR) GO TO 9005

      ERROR = .NOT. (0 .LE. COMPS .AND. 0 .LE. POINTS .AND.
     +   0 .LE. GROUPA .AND. 0 .LE. GROUPB)
      IF (ERROR) GO TO 9006

      ERROR = .NOT. ((0 .LT. COMPS) .EQV. (0 .LT. POINTS))
      IF (ERROR) GO TO 9007

      ERROR = .NOT. (0 .LT. GROUPA + COMPS * POINTS + GROUPB)
      IF (ERROR) GO TO 9008

      ERROR = .NOT. (NAMES .EQ. 1 .OR.
     +   NAMES .EQ. GROUPA + COMPS + GROUPB)
      IF (ERROR) GO TO 9009

      ERROR = .NOT. (POINTS .LE. PMAX)
      IF (ERROR) GO TO 9010

      COUNT = 0
      DO 1020 J = 1, GROUPA + COMPS + GROUPB
         IF (.NOT. (BELOW(J) .LT. ABOVE(J))) COUNT = COUNT + 1
1020  CONTINUE
      ERROR = COUNT .NE. 0
      IF (ERROR) GO TO 9011

C///  PARTITION THE INTEGER WORK SPACE.

C     SUBROUTINE TWGRAB (ERROR, LAST, FIRST, NUMBER)

      ILAST = 0

C     VARY(PMAX)
      CALL TWGRAB (ERROR, ILAST, QVARY, PMAX)
      IF (ERROR) GO TO 9012

C     VARY1(PMAX)
      CALL TWGRAB (ERROR, ILAST, QVARY1, PMAX)
      IF (ERROR) GO TO 9012

C     VARY2(PMAX)
      CALL TWGRAB (ERROR, ILAST, QVARY2, PMAX)
      IF (ERROR) GO TO 9012

C///  PARTITION THE REAL WORK SPACE.

C     SUBROUTINE TWGRAB (ERROR, LAST, FIRST, NUMBER)

      RLAST = 0

C     ABOVE(GROUPA + COMPS * PMAX + GROUPB)
      CALL TWGRAB (ERROR, RLAST, QABOVE, GROUPA + COMPS * PMAX + GROUPB)
      IF (ERROR) GO TO 9012

C     BELOW(GROUPA + COMPS * PMAX + GROUPB)
      CALL TWGRAB (ERROR, RLAST, QBELOW, GROUPA + COMPS * PMAX + GROUPB)
      IF (ERROR) GO TO 9012

C     RATIO1(PMAX)
      CALL TWGRAB (ERROR, RLAST, QRAT1, PMAX)
      IF (ERROR) GO TO 9012

C     RATIO2(PMAX)
      CALL TWGRAB (ERROR, RLAST, QRAT2, PMAX)
      IF (ERROR) GO TO 9012

C     S0(GROUPA + COMPS * PMAX + GROUPB)
      CALL TWGRAB (ERROR, RLAST, QS0, GROUPA + COMPS * PMAX + GROUPB)
      IF (ERROR) GO TO 9012

C     S1(GROUPA + COMPS * PMAX + GROUPB)
      CALL TWGRAB (ERROR, RLAST, QS1, GROUPA + COMPS * PMAX + GROUPB)
      IF (ERROR) GO TO 9012

C     USAVE(GROUPA + COMPS * PMAX + GROUPB)
      CALL TWGRAB (ERROR, RLAST, QUSAVE, GROUPA + COMPS * PMAX + GROUPB)
      IF (ERROR) GO TO 9012

C     VSAVE(GROUPA + COMPS * PMAX + GROUPB)
      CALL TWGRAB (ERROR, RLAST, QVSAVE, GROUPA + COMPS * PMAX + GROUPB)
      IF (ERROR) GO TO 9012

C     V1(GROUPA + COMPS * PMAX + GROUPB)
      CALL TWGRAB (ERROR, RLAST, QV1, GROUPA + COMPS * PMAX + GROUPB)
      IF (ERROR) GO TO 9012

C     XSAVE(PMAX)
      CALL TWGRAB (ERROR, RLAST, QXSAVE, PMAX)
      IF (ERROR) GO TO 9012

C     Y0(GROUPA + COMPS * PMAX + GROUPB)
      CALL TWGRAB (ERROR, RLAST, QY0, GROUPA + COMPS * PMAX + GROUPB)
      IF (ERROR) GO TO 9012

C     Y1(GROUPA + COMPS * PMAX + GROUPB)
      CALL TWGRAB (ERROR, RLAST, QY1, GROUPA + COMPS * PMAX + GROUPB)
      IF (ERROR) GO TO 9012

C///  CHECK THE WORK SPACES' SIZES.

      ERROR = .NOT. (ILAST .LE. ISIZE .AND. RLAST .LE. RSIZE)
      IF (ERROR) GO TO 9013

C///  ONE-TIME INITIALIZATION.

C     ALLOW FURTHER TIME EVOLUTION
      ALLOW = .TRUE.

C     PRESENT TASK
      QTASK = QENTRY

C     STATISTICS ARRAYS
      DO 1040 K = 1, QTOTAL
         TOTAL(K) = 0.0
         DO 1030 J = 1, GMAX
            DETAIL(J, K) = 0.0
            EVENT(J, K) = 0
1030     CONTINUE
1040  CONTINUE

C     TOTAL TIME STATISTIC
      CALL TWTIME (TIMER(QTOTAL))

C     GRID POINTER AND STATISTICS FOR THE FIRST GRID
      GRID = 1
      SIZE(GRID) = POINTS
      CALL TWTIME (TIMER(QGRID))

C     TIME STEP NUMBER
      STEP = 0

C     SOLUTION FLAG
      FOUND = .TRUE.

C///  EXPAND THE BOUNDS.

      COUNT = 0

      DO 1050 J = 1, GROUPA
         RWORK(QABOVE + COUNT) = ABOVE(J)
         RWORK(QBELOW + COUNT) = BELOW(J)
         COUNT = COUNT + 1
1050  CONTINUE

      DO 1070 K = 1, POINTS
         DO 1060 J = 1, COMPS
            RWORK(QABOVE + COUNT) = ABOVE(GROUPA + J)
            RWORK(QBELOW + COUNT) = BELOW(GROUPA + J)
            COUNT = COUNT + 1
1060     CONTINUE
1070  CONTINUE

      DO 1080 J = 1, GROUPB
         RWORK(QABOVE + COUNT) = ABOVE(GROUPA + COMPS + J)
         RWORK(QBELOW + COUNT) = BELOW(GROUPA + COMPS + J)
         COUNT = COUNT + 1
1080  CONTINUE

C///  SAVE THE INITIAL SOLUTION.

      PSAVE = POINTS
      IF (ADAPT .AND. 0 .LT. POINTS)
     +   CALL TWCOPY (POINTS, X, RWORK(QXSAVE))
      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, U, RWORK(QUSAVE))

C     GO TO 1090 WHEN RETURN = 1
      RETURN = 1
      GO TO 9911
1090  CONTINUE

C///  PRINT LEVELS 11, 21, AND 22.

      IF (0 .LT. LEVELD .AND. 0 .LT. TEXT) THEN
         WRITE (TEXT, 10002) ID, 'INITIAL GUESS:'
C        GO TO 1100 WHEN RETURN = 2
         RETURN = 2
         GO TO 9921
      END IF
1100  CONTINUE

C///  PRINT LEVEL 10 AND 11.

C                  123456789_123456789_123456789_1234
C                  12345678   123456  123456   123456
      HEADER(1) = '            LOG10   LOG10         '
      HEADER(2) = '    TASK   NORM F  COND J   REMARK'

      IF (LEVELM .EQ. 1 .AND. 0 .LT. TEXT) THEN
         IF (0 .LT. LEVELD) WRITE (TEXT, 10002) ID,
     +      'SOLVE THE PROBLEM.'
         WRITE (TEXT, 10003) (HEADER(J), J = 1, 2)
C        GO TO 1110 WHEN LABEL = 1
         LABEL = 1
         GO TO 7010
      END IF
1110  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     DECISION BLOCK.  THE PREVIOUS TASK DETERMINES THE NEXT.
C
C///////////////////////////////////////////////////////////////////////

2010  CONTINUE

C///  ENTRY WAS THE PREVIOUS TASK.

      IF (QTASK .EQ. QENTRY) THEN
         IF (0 .LT. STEPS0) THEN
            QTASK = QTIMST
            DESIRE = STEPS0
         ELSE IF (STEADY) THEN
            QTASK = QSEARC
         ELSE
            ERROR = .TRUE.
            GO TO 9014
         END IF

C///  SEARCH WAS THE PREVIOUS TASK.

      ELSE IF (QTASK .EQ. QSEARC) THEN
         IF (FOUND) THEN
            IF (ADAPT) THEN
               QTASK = QREFIN
            ELSE
               QTASK = QEXIT
               REPORT = ' '
            END IF
         ELSE
            IF (ALLOW .AND. 0 .LT. STEPS1) THEN
               QTASK = QTIMST
               DESIRE = STEPS1
            ELSE
               QTASK = QEXIT
               IF (1 .LT. GRID) THEN
                  REPORT = 'SOME SOLVED'
               ELSE
                  REPORT = 'NOT SOLVED'
               END IF
            END IF
         END IF

C///  REFINE WAS THE PREVIOUS TASK.

      ELSE IF (QTASK .EQ. QREFIN) THEN
         IF (FOUND) THEN
            STEP = 0
            QTASK = QSEARC
            ALLOW = .TRUE.
         ELSE
            QTASK = QEXIT
            IF (SATISF) THEN
               REPORT = ' '
            ELSE
               REPORT = 'NO SPACE'
            END IF
         END IF

C///  EVOLVE WAS THE PREVIOUS TASK.

      ELSE IF (QTASK .EQ. QTIMST) THEN
         IF (FOUND) THEN
            IF (STEADY) THEN
               QTASK = QSEARC
            ELSE
               QTASK = QEXIT
               REPORT = ' '
            END IF
         ELSE
            QTASK = QEXIT
            IF (1 .LT. GRID) THEN
               REPORT = 'SOME SOLVED'
            ELSE
               REPORT = 'NOT SOLVED'
            END IF
         END IF
      END IF

C///  BRANCH TO THE NEXT TASK.

      IF (QTASK .EQ. QEXIT) GO TO 3010
      IF (QTASK .EQ. QSEARC) GO TO 4010
      IF (QTASK .EQ. QREFIN) GO TO 5010
      IF (QTASK .EQ. QTIMST) GO TO 6010
      ERROR = .TRUE.
      GO TO 9015

C///////////////////////////////////////////////////////////////////////
C
C     EXIT BLOCK.
C
C///////////////////////////////////////////////////////////////////////

3010  CONTINUE

C///  COMPLETE STATISTICS FOR THE LAST GRID.

      CALL TWLAPS (TIMER(QGRID))
      IF (GRID .LE. GMAX) THEN
         DETAIL(GRID, QGRID) = TIMER(QGRID)
         DETAIL(GRID, QOTHER)
     +      = DETAIL(GRID, QGRID) - (DETAIL(GRID, QFUNCT)
     +      + DETAIL(GRID, QJACOB) + DETAIL(GRID, QSOLVE))
      END IF

C///  RESTORE THE SOLUTION.

      IF (REPORT .NE. ' ') THEN
C        BE CAREFUL NOT TO ASSIGN A VALUE TO A PARAMETER
         IF (POINTS .NE. PSAVE) POINTS = PSAVE
         IF (ADAPT .AND. 0 .LT. POINTS)
     +      CALL TWCOPY (POINTS, RWORK(QXSAVE), X)
         CALL TWCOPY
     +      (GROUPA + COMPS * POINTS + GROUPB, RWORK(QUSAVE), U)
      END IF

C///  PRINT LEVEL 11 OR 21.

C     SAVE THE STATUS REPORTS DURING REVERSE COMMUNICATION
      STRING = REPORT

      IF (LEVELD .EQ. 1 .AND. 0 .LT. TEXT) THEN
         WRITE (TEXT, 10002) ID, 'FINAL SOLUTION:'
C        GO TO 3020 WHEN RETURN = 3
         RETURN = 3
         GO TO 9921
      END IF
3020  CONTINUE

C     RESTORE THE STATUS REPORTS AFTER REVERSE COMMUNICATION
      REPORT = STRING

C///  COMPLETE THE TOTAL TIME STATISTICS.

      CALL TWLAPS (TIMER(QTOTAL))
      TOTAL(QTOTAL) = TIMER(QTOTAL)
      TOTAL(QOTHER) = TOTAL(QTOTAL)
     +   - (TOTAL(QFUNCT) + TOTAL(QJACOB) + TOTAL(QSOLVE))

C///  TOP OF THE REPORT BLOCK.

      IF (0 .LT. LEVELM .AND. 0 .LT. TEXT) THEN
         IF (0.0 .LT. TOTAL(QTOTAL)) THEN

C///  REPORT TOTAL COMPUTER TIME.

      TEMP = TOTAL(QTOTAL)
      IF (3600.0 .LE. TEMP) THEN
         WRITE (STRING, '(F10.2, A)') TEMP / 3600.0, ' HOURS'
      ELSE IF (60.0 .LE. TEMP) THEN
         WRITE (STRING, '(F10.2, A)') TEMP / 60.0, ' MINUTES'
      ELSE
         WRITE (STRING, '(F10.2, A)') TEMP, ' SECONDS'
      END IF

      CALL TWSQEZ (LENGTH, STRING)
      WRITE (TEXT, 10004) ID, STRING (1 : LENGTH)

C///  REPORT PERCENT OF TOTAL COMPUTER TIME.

      TEMP = 100.0 / TOTAL(QTOTAL)
      IF (ADAPT) THEN

C                  123456789_123456789_123456789_12345678
C                  123456  123456  123456 123456 123456
      HEADER(1) = '                TASK                  '
      HEADER(3) = '  GRID    GRID  --------------------  '
      HEADER(5) = 'POINTS  TOTALS  EVOLVE SEARCH REFINE  '

C                  123456789_123456789_1234567
C                  123456 123456 123456 123456
      HEADER(2) = 'SUBTASK                    '
      HEADER(4) = '---------------------------'
      HEADER(6) = 'EVAL F PREP J  SOLVE  OTHER'

      WRITE (TEXT, 10005) HEADER,
     +   (SIZE(J), (TEMP * DETAIL(J, K), K = 1, 8), J = 1, GRID)
      IF (1 .LT. GRID) WRITE (TEXT, 10006) (TEMP * TOTAL(K), K = 2, 8)
      IF (GMAX .LT. GRID) WRITE (TEXT, 10007)

      ELSE

C                  123456789_123456789_123456789_123456789_123456789_1
C                  123456   123456   123456   123456   123456   123456
      HEADER(1) = 'SUBTASK                             TASK           '
      HEADER(2) = '---------------------------------   ---------------'
      HEADER(3) = 'EVAL F   PREP J    SOLVE    OTHER   EVOLVE   SEARCH'

      WRITE (TEXT, 10008)
     +   (HEADER(J), J = 1, 3), '  % OF TOTAL',
     +   (TEMP * TOTAL(K), K = 5, 8), (TEMP * TOTAL(K), K = 2, 3),
     +   'MEAN SECONDS', (DETAIL(1, K) / EVENT(1, K), K = 5, 7),
     +   '    QUANTITY', (EVENT(1, K), K = 5, 7)

      END IF

C///  REPORT AVERAGE COMPUTER TIME.

C                  123456789_123456789_123456789_1234567
C                  123456   1234567  1234567  1234567
      HEADER(1) = '         AVERAGE SECONDS             '
      HEADER(3) = '  GRID   -------------------------   '
      HEADER(5) = 'POINTS    EVAL F   PREP J    SOLVE   '


C                  123456789_123456789_12345
C                  1234567  1234567  1234567
      HEADER(2) = 'NUMBER OF SUBTASKS       '
      HEADER(4) = '-------------------------'
      HEADER(6) = ' EVAL F   PREP J    SOLVE'

      IF (ADAPT) WRITE (TEXT, 10009) HEADER,
     +   (SIZE(J), (DETAIL(J, K) / EVENT(J, K), K = 5, 7),
     +   (EVENT(J, K), K = 5, 7), J = 1, GRID)

      END IF

C///  REPORT THE COMPLETION STATUS.

      IF (0 .LT. LEVELM) THEN
         IF (REPORT .EQ. ' ') THEN
            WRITE (TEXT, 10010) ID
         ELSE IF (REPORT .EQ. 'NO SPACE') THEN
            WRITE (STRING, '(I10)') POINTS
            CALL TWSQEZ (LENGTH, STRING)
            WRITE (TEXT, 10011)
     +         ID, STRING (1 : LENGTH), RATIO, TOLER1, TOLER2
         ELSE IF (REPORT .EQ. 'NOT SOLVED') THEN
            WRITE (TEXT, 10012) ID
         ELSE IF (REPORT .EQ. 'SOME SOLVED') THEN
            WRITE (STRING, '(I10)') POINTS
            CALL TWSQEZ (LENGTH, STRING)
            WRITE (TEXT, 10013)
     +         ID, STRING (1 : LENGTH), RATIO, TOLER1, TOLER2
         ELSE
            ERROR = .TRUE.
            GO TO 9016
         END IF
      END IF

C///  BOTTOM OF THE REPORT BLOCK.

      END IF

C///  BOTTOM OF THE EXIT BLOCK.

      GO TO 99999

C///////////////////////////////////////////////////////////////////////
C
C     SEARCH BLOCK.
C
C///////////////////////////////////////////////////////////////////////

4010  CONTINUE

C///  INITIALIZE STATISTICS ON ENTRY TO THE SEARCH BLOCK.

      CALL TWTIME (TIMER(QSEARC))
      FIRST = .TRUE.
      JACOBS = 0
      MAXCON = 0.0

C///  PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE SEARCH BLOCK.

      IF (1 .LT. LEVELM) THEN
         IF (0 .LT. TEXT) WRITE (TEXT, 10014) ID
      END IF

C///  PREPARE TO CALL SEARCH.

C     SAVE THE SOLUTION SHOULD THE SEARCH FAIL
      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, U, RWORK(QVSAVE))

      EXIST = .FALSE.

C///  CALL SEARCH.

      AGE = 0
4020  CONTINUE

C     SUBROUTINE SEARCH
C    +  (ERROR, TEXT,
C    +   ABOVE, AGE, BELOW, BUFFER, COMPS, CONDIT, EXIST, GROUPA,
C    +   GROUPB, LEVELD, LEVELM, NAME, NAMES, POINTS, REPORT, S0, S1,
C    +   SIGNAL, STEPS, SUCCES, V0, V1, XXABS, XXAGE, XXREL, Y0, Y0NORM,
C    +   Y1)

      CALL SEARCH
     +  (ERROR, TEXT,
     +   RWORK(QABOVE), AGE, RWORK(QBELOW), BUFFER, COMPS, CONDIT,
     +   EXIST, GROUPA, GROUPB, LEVELD - 1, LEVELM - 1, NAME, NAMES,
     +   POINTS, XREPOR, RWORK(QS0), RWORK(QS1), SIGNAL, NSTEPS, FOUND,
     +   U, RWORK(QV1), SSABS, SSAGE, SSREL, RWORK(QY0), YNORM,
     +   RWORK(QY1))
      IF (ERROR) GO TO 9017

C///  PASS REQUESTS FROM SEARCH TO THE CALLER.

      IF (SIGNAL .NE. ' ') THEN
C        GO TO 4020 WHEN RETURN = 4
         RETURN = 4
         GO TO 9931
      END IF

C///  REACT TO THE COMPLETION OF SEARCH.

      IF (FOUND) THEN
C        SAVE THE LATEST SOLUTION

         PSAVE = POINTS
         IF (ADAPT .AND. 0 .LT. POINTS)
     +      CALL TWCOPY (POINTS, X, RWORK(QXSAVE))
         CALL TWCOPY
     +      (GROUPA + COMPS * POINTS + GROUPB, U, RWORK(QUSAVE))

C        GO TO 4030 WHEN RETURN = 5
         RETURN = 5
         GO TO 9911
      ELSE
C        RESTORE THE SOLUTION
         CALL TWCOPY
     +      (GROUPA + COMPS * POINTS + GROUPB, RWORK(QVSAVE), U)
      END IF
4030  CONTINUE

C///  COMPLETE STATISTICS FOR THE SEARCH BLOCK.

      CALL TWLAPS (TIMER(QSEARC))
      TOTAL(QSEARC) = TOTAL(QSEARC) + TIMER(QSEARC)
      IF (GRID .LE. GMAX)
     +   DETAIL(GRID, QSEARC) = DETAIL(GRID, QSEARC) + TIMER(QSEARC)

C///  PRINT LEVEL 10 OR 11 ON EXIT FROM THE SEARCH BLOCK.

      IF (LEVELM .EQ. 1 .AND. 0 .LT. TEXT) THEN
C        GO TO 4040 WHEN LABEL = 2
         LABEL = 2
         GO TO 7010
      END IF
4040  CONTINUE

C///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE SEARCH BLOCK.

      IF (1 .LT. LEVELM) THEN
         IF (FOUND) THEN
            IF (0 .LT. TEXT) WRITE (TEXT, 10015) ID
         ELSE
            IF (0 .LT. TEXT) WRITE (TEXT, 10016) ID
         END IF
      END IF

C///  BOTTOM OF THE SEARCH BLOCK.

      GO TO 2010

C///////////////////////////////////////////////////////////////////////
C
C     REFINE BLOCK.
C
C///////////////////////////////////////////////////////////////////////

5010  CONTINUE

C///  INITIALIZE STATISTICS ON ENTRY TO THE REFINE BLOCK.

      CALL TWTIME (TIMER(QREFIN))

C///  PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE REFINE BLOCK.

      IF (1 .LT. LEVELM) THEN
         IF (0 .LT. TEXT) WRITE (TEXT, 10017) ID
      END IF

C///  PREPARE TO CALL REFINE.

C     SAVE THE GROUP B VALUES
      DO 5020 J = 1, GROUPB
         RWORK(QVSAVE - 1 + J) = U(GROUPA + COMPS * POINTS + J)
5020  CONTINUE

      EXIST = .FALSE.

C///  CALL REFINE.

5030  CONTINUE

C     SUBROUTINE REFINE
C    +  (ERROR, TEXT,
C    +   ACTIVE, BUFFER, COMPS, LEVELD, LEVELM, MARK, NEWX, PADD, PMAX,
C    +   POINTS, RATIO, RATIO1, RATIO2, SIGNAL, SUCCES, TOLER0, TOLER1,
C    +   TOLER2, U, VARY1, VARY2, WEIGHT, X)

      CALL REFINE
     +  (ERROR, TEXT,
     +   ACTIVE,
     +   BUFFER(GROUPA + 1), COMPS, LEVELD - 1, LEVELM - 1, MARK,
     +   FOUND, PADD, PMAX, POINTS, RATIO, RWORK(QRAT1), RWORK(QRAT2),
     +   SIGNAL, SATISF, TOLER0, TOLER1, TOLER2, U(GROUPA + 1),
     +   IWORK(QVARY1), IWORK(QVARY2), IWORK(QVARY), X)
      IF (ERROR) GO TO 9018

C///  SERVICE REQUESTS FROM REFINE: PASS REQUESTS TO THE CALLER.

      IF (SIGNAL .NE. ' ') THEN
C        INSERT THE GROUP A AND B UNKNOWNS
         DO 5040 J = 1, GROUPA
            BUFFER(J) = U(J)
5040     CONTINUE
         DO 5050 J = 1, GROUPB
            BUFFER(GROUPA + COMPS * POINTS + J) = RWORK(QVSAVE - 1 + J)
5050     CONTINUE

C        GO TO 5030 WHEN RETURN = 6
         RETURN = 6
         GO TO 9931
      END IF

C///  REACT TO THE COMPLETION OF REFINE.

      IF (.NOT. FOUND) GO TO 5110

C        COMPLETE STATISTICS FOR THE OLD GRID
         CALL TWLAPS (TIMER(QGRID))
         IF (GRID .LE. GMAX) THEN
            DETAIL(GRID, QGRID) = TIMER(QGRID)
            DETAIL(GRID, QOTHER)
     +         = DETAIL(GRID, QGRID) - (DETAIL(GRID, QFUNCT)
     +         + DETAIL(GRID, QJACOB) + DETAIL(GRID, QSOLVE))
         END IF

C        INITIALIZE STATISTICS FOR THE NEW GRID
         GRID = GRID + 1
         IF (GRID .LE. GMAX) THEN
            CALL TWTIME (TIMER(QGRID))
            SIZE(GRID) = POINTS
         END IF

C        INSERT THE GROUP B VALUES
         DO 5060 J = 1, GROUPB
            U(GROUPA + COMPS * POINTS + J) = RWORK(QVSAVE - 1 + J)
5060     CONTINUE

C        EXPAND THE BOUNDS
         COUNT = GROUPA

         DO 5080 K = 1, POINTS
            DO 5070 J = 1, COMPS
               RWORK(QABOVE + COUNT) = ABOVE(GROUPA + J)
               RWORK(QBELOW + COUNT) = BELOW(GROUPA + J)
               COUNT = COUNT + 1
5070        CONTINUE
5080     CONTINUE

         DO 5090 J = 1, GROUPB
            RWORK(QABOVE + COUNT) = ABOVE(GROUPA + COMPS + J)
            RWORK(QBELOW + COUNT) = BELOW(GROUPA + COMPS + J)
            COUNT = COUNT + 1
5090     CONTINUE

C        SAVE THE LATEST SOLUTION
C        GO TO 5100 WHEN RETURN = 7
         RETURN = 7
         GO TO 9911
5100     CONTINUE

5110  CONTINUE

C///  COMPLETE STATISTICS FOR THE REFINE BLOCK.

      CALL TWLAPS (TIMER(QREFIN))
      TOTAL(QREFIN) = TOTAL(QREFIN) + TIMER(QREFIN)
      IF (GRID .LE. GMAX)
     +   DETAIL(GRID, QREFIN) = DETAIL(GRID, QREFIN) + TIMER(QREFIN)

C///  PRINT LEVEL 10 OR 11 ON EXIT FROM THE REFINE BLOCK.

      IF (LEVELM .EQ. 1 .AND. 0 .LT. TEXT) THEN
         WRITE (TEXT, '()')
C        GO TO 5120 WHEN LABEL = 3
         LABEL = 3
         GO TO 7010
      END IF
5120  CONTINUE

C///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE REFINE BLOCK.

      IF (1 .LT. LEVELM) THEN
         IF (FOUND) THEN
            IF (0 .LT. TEXT) WRITE (TEXT, 10018) ID
         ELSE
            IF (0 .LT. TEXT) WRITE (TEXT, 10019) ID
         END IF
      END IF

C///  BOTTOM OF THE REFINE BLOCK.

      GO TO 2010

C///////////////////////////////////////////////////////////////////////
C
C     EVOLVE BLOCK.
C
C///////////////////////////////////////////////////////////////////////

6010  CONTINUE

C///  INITIALIZE STATISTICS ON ENTRY TO THE EVOLVE BLOCK.

      CALL TWTIME (TIMER(QTIMST))
      FIRST = .TRUE.
      JACOBS = 0
      MAXCON = 0.0
      STEPS = STEP

C///  PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE EVOLVE BLOCK.

      IF (1 .LT. LEVELM) THEN
         IF (0 .LT. TEXT) WRITE (TEXT, 10020) ID
      END IF

C///  CALL EVOLVE.

6020  CONTINUE

C     SUBROUTINE EVOLVE
C    +  (ERROR, TEXT,
C    +   ABOVE, BELOW, BUFFER, COMPS, CONDIT, DESIRE, GROUPA, GROUPB,
C    +   LEVELD, LEVELM, NAME, NAMES, POINTS, REPORT, S0, S1, SIGNAL,
C    +   STEP, STEPS2, STRID0, STRIDE, SUCCES, TDABS, TDAGE, TDEC,
C    +   TDREL, TIME, TINC, TMAX, TMIN, V0, V1, VSAVE, Y0, Y1, YNORM)

      CALL EVOLVE
     +  (ERROR, TEXT,
     +   RWORK(QABOVE), RWORK(QBELOW), BUFFER, COMPS, CONDIT, DESIRE,
     +   GROUPA, GROUPB, LEVELD - 1, LEVELM - 1, NAME, NAMES, POINTS,
     +   XREPOR, RWORK(QS0), RWORK(QS1), SIGNAL, STEP, STEPS2, STRID0,
     +   STRIDE, FOUND, TDABS, TDAGE, TDEC, TDREL, TIME, TINC, TMAX,
     +   TMIN, U, RWORK(QV1), RWORK(QVSAVE), RWORK(QY0), RWORK(QY1),
     +   YNORM)
      IF (ERROR) GO TO 9019

C///  PASS REQUESTS FROM EVOLVE TO THE CALLER.

      IF (SIGNAL .NE. ' ') THEN
C        GO TO 6020 WHEN RETURN = 8
         RETURN = 8
         GO TO 9931
      END IF

C///  REACT TO THE COMPLETION OF EVOLVE.

      IF (FOUND) THEN
C        SAVE THE LATEST SOLUTION
C        GO TO 6030 WHEN RETURN = 9
         RETURN = 9
         GO TO 9911
      END IF
6030  CONTINUE

C///  ALLOW FURTHER TIME EVOLUTION.

      ALLOW = XREPOR .EQ. QNULL

C///  COMPLETE STATISTICS FOR THE EVOLVE BLOCK.

      CALL TWLAPS (TIMER(QTIMST))
      TOTAL(QTIMST) = TOTAL(QTIMST) + TIMER(QTIMST)
      IF (GRID .LE. GMAX)
     +   DETAIL(GRID, QTIMST) = DETAIL(GRID, QTIMST) + TIMER(QTIMST)
      STEPS = STEP - STEPS

C///  PRINT LEVEL 10 OR 11 ON EXIT FROM THE EVOLVE BLOCK.

      IF (LEVELM .EQ. 1 .AND. 0 .LT. TEXT) THEN
C        GO TO 6040 WHEN LABEL = 4
         LABEL = 4
         GO TO 7010
      END IF
6040  CONTINUE

C///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE EVOLVE BLOCK.

      IF (1 .LT. LEVELM) THEN
         IF (FOUND) THEN
            IF (0 .LT. TEXT) WRITE (TEXT, 10021) ID
         ELSE
            IF (0 .LT. TEXT) WRITE (TEXT, 10022) ID
         END IF
      END IF

C///  BOTTOM OF THE EVOLVE BLOCK.

      GO TO 2010

C///////////////////////////////////////////////////////////////////////
C
C     BLOCK TO PRINT LOG LINES.
C
C///////////////////////////////////////////////////////////////////////

7010  CONTINUE

      DO 7020 J = 1, 3
         COLUMN(J) = ' '
7020  CONTINUE

      STRING = ' '

C     COLUMN 1: NAME OF THE TASK
      IF (QTASK .EQ. QENTRY) COLUMN(1) = '   START'
      IF (QTASK .EQ. QSEARC) COLUMN(1) = '  SEARCH'
      IF (QTASK .EQ. QREFIN) COLUMN(1) = '  REFINE'
      IF (QTASK .EQ. QTIMST) COLUMN(1) = '  EVOLVE'

C     COLUMN 2: NORM OF THE STEADY STATE FUNCTION
      IF (.NOT. FOUND) GO TO 7040
C        GO TO 7030 WHEN RETURN = 10
         RETURN = 10
         GO TO 9941
7030     CONTINUE
         CALL TWNORM (GROUPA + COMPS * POINTS + GROUPB, TEMP, BUFFER)
         CALL TWLOGR (COLUMN(2), TEMP)
7040  CONTINUE

C     COLUMN 3: LARGEST CONDITION NUMBER
      IF (QTASK .EQ. QSEARC .OR. QTASK .EQ. QTIMST) THEN
         IF (MAXCON .NE. 0.0) CALL TWLOGR (COLUMN(3), MAXCON)
      END IF

C     REMARK
      IF (QTASK .EQ. QSEARC) THEN
         IF (XREPOR .EQ. QDVRG) THEN
            STRING = 'DIVERGING'
         ELSE IF (XREPOR .EQ. QNULL) THEN
            IF (NSTEPS .EQ. 1) THEN
               WRITE (STRING, '(I10, A)') NSTEPS, ' SEARCH STEP'
            ELSE
               WRITE (STRING, '(I10, A)') NSTEPS, ' SEARCH STEPS'
            END IF
         ELSE IF (XREPOR .EQ. QBNDS) THEN
            STRING = 'GOING OUT OF BOUNDS'
         ELSE
            STRING = '?'
         END IF
      ELSE IF (QTASK .EQ. QTIMST) THEN
         IF (XREPOR .EQ. QBNDS .OR. XREPOR .EQ. QDVRG .OR.
     +      XREPOR .EQ. QNULL) THEN
            WRITE (STRING, '(I10, A, 1P, E10.1, A)')
     +         STEPS, ' TIME STEPS, ', STRIDE, ' LAST STRIDE'
         ELSE
            STRING = '?'
         END IF
      ELSE IF (QTASK .EQ. QENTRY .AND. ADAPT) THEN
         WRITE (STRING, '(I10, A)') POINTS, ' GRID POINTS'
      ELSE IF (QTASK .EQ. QREFIN) THEN
         IF (FOUND) THEN
            WRITE (STRING, '(F10.2, A, F10.2, A, I10, A)')
     +         RATIO(1), ' AND ', RATIO(2), ' RATIOS, ', POINTS,
     +         ' GRID POINTS'
         ELSE
            WRITE (STRING, '(F10.2, A, F10.2, A)')
     +         RATIO(1), ' AND ', RATIO(2), ' RATIOS'
         END IF
      END IF

      CALL TWSQEZ (LENGTH, STRING)
      IF (0 .LT. TEXT) WRITE (TEXT, 10023) COLUMN, STRING (1 : LENGTH)

      GO TO (1110, 4040, 5120, 6040) LABEL
      ERROR = .TRUE.
      GO TO 9020

C///////////////////////////////////////////////////////////////////////
C
C     REQUEST REVERSE COMMUNICATION.
C
C///////////////////////////////////////////////////////////////////////

C///  SAVE THE SOLUTION.

9911  CONTINUE

      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, U, BUFFER)
      SIGNAL = 'SAVE'
C     GO TO 9912 WHEN ROUTE = 1
      ROUTE = 1
      GO TO 99999
9912  CONTINUE
      SIGNAL = ' '

      GO TO (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030)
     +   RETURN
      ERROR = .TRUE.
      GO TO 9021

C///  PRINT THE LATEST SOLUTION.

9921  CONTINUE

      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, U, BUFFER)
      SIGNAL = 'SHOW'
C     GO TO 9922 WHEN ROUTE = 2
      ROUTE = 2
      GO TO 99999
9922  CONTINUE
      SIGNAL = ' '

      GO TO (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030)
     +   RETURN
      ERROR = .TRUE.
      GO TO 9021

C///  PASS REQUESTS FROM SEARCH, REFINE, OR EVOLVE TO THE CALLER.

9931  CONTINUE

C     IDENTIFY THE REQUEST.  THIS MUST BE SAVED TO GATHER STATISTICS
C     AT REENTRY.  THE REVERSE COMMUNICATION FLAGS WILL NOT BE SAVED
C     BECAUSE THEY ARE CLEARED AT EVERY ENTRY.
      IF (SIGNAL .EQ. 'RESIDUAL') THEN
         QTYPE = QFUNCT
      ELSE IF (SIGNAL .EQ. 'PREPARE') THEN
         QTYPE = QJACOB
      ELSE IF (SIGNAL .EQ. 'SOLVE') THEN
         QTYPE = QSOLVE
      ELSE
         QTYPE = QOTHER
      END IF

C     COUNT THE JACOBIANS
      IF (QTYPE .EQ. QJACOB) JACOBS = JACOBS + 1

      CALL TWTIME (TIMER(QTYPE))

C     GO TO 9932 WHEN ROUTE = 3
      ROUTE = 3
      GO TO 99999
9932  CONTINUE

C     SAVE THE CONDITION NUMBER
      IF (QTYPE .EQ. QJACOB) MAXCON = MAX (MAXCON, CONDIT)

      CALL TWLAPS (TIMER(QTYPE))
      TOTAL(QTYPE) = TOTAL(QTYPE) + TIMER(QTYPE)
      IF (GRID .LE. GMAX) THEN
         DETAIL(GRID, QTYPE) = DETAIL(GRID, QTYPE) + TIMER(QTYPE)
         EVENT(GRID, QTYPE) = EVENT(GRID, QTYPE) + 1
      END IF

      GO TO (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030)
     +   RETURN
      ERROR = .TRUE.
      GO TO 9021

C///  EVALUATE THE STEADY STATE FUNCTION.

9941  CONTINUE
      CALL TWTIME (TIMER(QFUNCT))

      CALL TWCOPY (GROUPA + COMPS * POINTS + GROUPB, U, BUFFER)
      SIGNAL = 'RESIDUAL'
      TIME = .FALSE.
C     GO TO 9942 WHEN ROUTE = 4
      ROUTE = 4
      GO TO 99999
9942  CONTINUE
      SIGNAL = ' '

      CALL TWLAPS (TIMER(QFUNCT))
      TOTAL(QFUNCT) = TOTAL(QFUNCT) + TIMER(QFUNCT)
      IF (GRID .LE. GMAX) THEN
         DETAIL(GRID, QFUNCT) = DETAIL(GRID, QFUNCT) + TIMER(QFUNCT)
         EVENT(GRID, QFUNCT) = EVENT(GRID, QFUNCT) + 1
      END IF

      GO TO (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030)
     +   RETURN
      ERROR = .TRUE.
      GO TO 9021

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +   (/1X, A9, A, ' (TWO POINT BOUNDARY VALUE PROBLEM) SOLVER,'
     +   /10X, 'VERSION ', A,
     +   ' OF APRIL 1998 BY DR. JOSEPH F. GRCAR.')

10002 FORMAT
     +   (/1X, A9, A)

10003 FORMAT
     +   (3(/10X, A35)/)

10004 FORMAT
     +  (/1X, A9, A, ' TOTAL COMPUTER TIME (SEE BREAKDOWN BELOW).')

10005 FORMAT
     +  (/10X, 'PERCENT OF TOTAL COMPUTER TIME FOR VARIOUS TASKS:'
     +   /3(/10X, A38, A27)
     +  //(10X, I6, 2X, F6.1, 1X, 3(1X, F6.1), 1X, 4(1X, F6.1)))

10006 FORMAT
     +  (/12X, 'TASK TOTALS:', 1X, 3(1X, F6.1), 1X, 4(1X, F6.1))

10007 FORMAT
     +  (/10X, 'SOME GRIDS ARE OMITTED, BUT THE TOTALS ARE FOR ALL.')

10008 FORMAT
     +  (3(/24X, A51)
     +  //10X, A12, F8.1, 5F9.1
     +   /10X, A12, F8.3, 2F9.3
     +   /10X, A12, I8, 2I9)

10009 FORMAT
     +  (/10X, 'AVERAGE COMPUTER TIMES FOR, AND NUMBERS OF, SUBTASKS:'
     +   /3(/10X, A37, A25)
     +  //(10X, I6, 3X, F7.3, 2X, F7.3, 2X, F7.3, 1X, 3(2X, I7)))

10010 FORMAT
     +  (/1X, A9, 'SUCCESS.  PROBLEM SOLVED.')

10011 FORMAT
     +  (/1X, A9, 'FAILURE.  A SOLUTION WAS FOUND FOR A GRID WITH ', A
     +  /10X, 'POINTS, BUT ONE OR BOTH RATIOS ARE TOO LARGE.'
C               123456789_  123456789_
     +  //22X, '   RATIO 1     RATIO 2'
     +  //10X, '     FOUND', 2F12.2
     +   /10X, '   DESIRED', 2F12.2
     +  //10X, 'A LARGER GRID COULD NOT BE FORMED.')

10012 FORMAT
     +  (/1X, A9, 'FAILURE.  NO SOLUTION WAS FOUND.')

10013 FORMAT
     +  (/1X, A9, 'FAILURE.  A SOLUTION WAS FOUND FOR A GRID WITH ', A
     +  /10X, 'POINTS, BUT ONE OR BOTH RATIOS ARE TOO LARGE.'
C               123456789_  123456789_
     +  //22X, '   RATIO 1     RATIO 2'
     +  //10X, '     FOUND', 2F12.2
     +   /10X, '   DESIRED', 2F12.2
     +  //10X, 'A SOLUTION COULD NOT BE FOUND FOR A LARGER GRID.')

10014 FORMAT
     +  (/1X, A9, 'CALLING SEARCH TO SOLVE THE STEADY STATE PROBLEM.')

10015 FORMAT
     +  (/1X, A9, 'SEARCH FOUND THE STEADY STATE.')

10016 FORMAT
     +  (/1X, A9, 'SEARCH DID NOT FIND THE STEADY STATE.')

10017 FORMAT
     +  (/1X, A9, 'CALLING REFINE TO PRODUCE A NEW GRID.')

10018 FORMAT
     +  (/1X, A9, 'REFINE SELECTED A NEW GRID.')

10019 FORMAT
     +  (/1X, A9, 'REFINE DID NOT SELECT A NEW GRID.')

10020 FORMAT
     +  (/1X, A9, 'CALLING EVOLVE TO PERFORM TIME EVOLUTION.')

10021 FORMAT
     +  (/1X, A9, 'EVOLVE PERFORMED A TIME EVOLUTION.')

10022 FORMAT
     +  (/1X, A9, 'EVOLVE DID NOT PERFORM A TIME EVOLUTION.')

10023 FORMAT
     +   (10X, A8, 3X, A6, 2X, A6, 3X, A)

80001 FORMAT
     +   ('(', A, ' ', I10, ')')

80002 FORMAT
     +   (10X, 1P, E10.2, 2X, E10.2, 3X, A)

80003 FORMAT
     +   (10X, '  ... MORE')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

C     GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, ROUTE
      IF (.NOT. MESS) GO TO 99999

9002  IF (0 .LT. TEXT) THEN
         CALL TWLAST (LENGTH, VERSIO)
         WRITE (TEXT, 99002) ID, VERSIO (1 : LENGTH), VNMBR(VNMBRS)
         DO 9901 J = VNMBRS - 1, 1, - 1
            WRITE (TEXT, '(10X, A, A)')
C*****PRECISION > DOUBLE
     +         ' CAN REPLACE:  DOUBLE PRECISION VERSION ', VNMBR(J)
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C     +         ' CAN REPLACE:  SINGLE PRECISION VERSION ', VNMBR(J)
C*****END PRECISION > SINGLE
9901     CONTINUE
      END IF
      IF (.NOT. MESS) GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID
      IF (.NOT. MESS) GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, CNTRLS, COUNT
      IF (.NOT. MESS) GO TO 99999

9005  IF (0 .LT. TEXT) WRITE (TEXT, 99005) ID, LEVELD, LEVELM
      IF (.NOT. MESS) GO TO 99999

9006  IF (0 .LT. TEXT) WRITE (TEXT, 99006) ID,
     +   COMPS, POINTS, GROUPA, GROUPB
      IF (.NOT. MESS) GO TO 99999

9007  IF (0 .LT. TEXT) WRITE (TEXT, 99007) ID, COMPS, POINTS
      IF (.NOT. MESS) GO TO 99999

9008  IF (0 .LT. TEXT) WRITE (TEXT, 99008) ID,
     +   COMPS, POINTS, GROUPA, GROUPB, GROUPA + COMPS * POINTS + GROUPB
      IF (.NOT. MESS) GO TO 99999

9009  IF (0 .LT. TEXT) WRITE (TEXT, 99009) ID,
     +   NAMES, COMPS, GROUPA, GROUPB, GROUPA + COMPS + GROUPB
      IF (.NOT. MESS) GO TO 99999

9010  IF (0 .LT. TEXT) WRITE (TEXT, 99010) ID, POINTS, PMAX
      IF (.NOT. MESS) GO TO 99999

9011  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99011) ID,
     +      GROUPA, GROUPB, COMPS, GROUPA + COMPS + GROUPB, COUNT
         COUNT = 0
         DO 8010 J = 1, GROUPA + COMPS + GROUPB
            IF (.NOT. (BELOW(J) .LT. ABOVE(J)) .OR. MESS) THEN
               COUNT = COUNT + 1
               IF (COUNT .LE. LINES) THEN
                  IF (NAMES .EQ. COMPS + GROUPA + GROUPB) THEN
                     CTEMP1 = NAME(J)
                  ELSE
                     CTEMP1 = ' '
                  END IF
                  CALL TWSQEZ (LEN1, CTEMP1)

                  IF (J .LE. GROUPA) THEN
                     WRITE (CTEMP2, 80001) 'A', J
                  ELSE IF (J .LE. GROUPA + COMPS) THEN
                     WRITE (CTEMP2, 80001) 'C', J - GROUPA
                  ELSE
                     WRITE (CTEMP2, 80001) 'B', J - GROUPA - COMPS
                  END IF
                  CALL TWSQEZ (LEN2, CTEMP2)

                  IF (CTEMP1 .EQ. ' ') THEN
                     STRING = CTEMP2
                     LENGTH = LEN2
                  ELSE IF (LEN1 + 2 + LEN2 .LE. 40) THEN
                     STRING = CTEMP1 (1 : LEN1) // '  ' // CTEMP2
                     LENGTH = LEN1 + 2 + LEN2
                  ELSE IF (LEN1 + 1 + LEN2 .LE. 40) THEN
                     STRING = CTEMP1 (1 : LEN1) // ' ' // CTEMP2
                     LENGTH = LEN1 + 1 + LEN2
                  ELSE
                     LEN1 = 40 - LEN2 - 4
                     STRING = CTEMP1 (1 : LEN1) // '... ' // CTEMP2
                     LENGTH = 40
                  END IF

                  WRITE (TEXT, 80002)
     +               BELOW(J), ABOVE(J), STRING (1 : LENGTH)
               END IF
            END IF
8010     CONTINUE
         IF (LINES .LT. COUNT) WRITE (TEXT, 80003)
      END IF
      IF (.NOT. MESS) GO TO 99999

9012  IF (0 .LT. TEXT) WRITE (TEXT, 99012) ID
      IF (.NOT. MESS) GO TO 99999

9013  IF (0 .LT. TEXT) WRITE (TEXT, 99013) ID,
     +   ISIZE, RSIZE, ILAST, RLAST
      IF (.NOT. MESS) GO TO 99999

9014  IF (0 .LT. TEXT) WRITE (TEXT, 99014) ID
      IF (.NOT. MESS) GO TO 99999

9015  IF (0 .LT. TEXT) WRITE (TEXT, 99015) ID
      IF (.NOT. MESS) GO TO 99999

9016  IF (0 .LT. TEXT) WRITE (TEXT, 99016) ID
      IF (.NOT. MESS) GO TO 99999

9017  IF (0 .LT. TEXT) WRITE (TEXT, 99017) ID
      IF (.NOT. MESS) GO TO 99999

9018  IF (0 .LT. TEXT) WRITE (TEXT, 99018) ID
      IF (.NOT. MESS) GO TO 99999

9019  IF (0 .LT. TEXT) WRITE (TEXT, 99019) ID
      IF (.NOT. MESS) GO TO 99999

9020  IF (0 .LT. TEXT) WRITE (TEXT, 99020) ID, LABEL
      IF (.NOT. MESS) GO TO 99999

9021  IF (0 .LT. TEXT) WRITE (TEXT, 99021) ID, RETURN
      IF (.NOT. MESS) GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, I10, '  ROUTE')

99002 FORMAT
     +  (/1X, A9, 'ERROR.  THE CALLING PROGRAM EXPECTS A VERSION OF'
     +  /10X, 'TWOPNT NOT COMPATIBLE WITH THIS VERSION.'
     + //10X, '     EXPECTS:  ', A
C*****PRECISION > DOUBLE
     + //10X, 'THIS VERSION:  DOUBLE PRECISION VERSION ', A)
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C     + //10X, 'THIS VERSION:  SINGLE PRECISION VERSION ', A)
C*****END PRECISION > SINGLE

99003 FORMAT
     +  (/1X, A9, 'ERROR.  TWINIT FAILS.')

99004 FORMAT
     +  (/1X, A9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.'
     + //10X, I10, '  CONTROLS'
     +  /10X, I10, '  COUNTED')

99005 FORMAT
     +  (/1X, A9, 'ERROR.  THE PRINTING LEVELS ARE OUT OF ORDER.'
     +  /10X, 'LEVELD CANNOT EXCEED LEVELM.'
     + //10X, I10, '  LEVELD, FOR SOLUTIONS'
     +  /10X, I10, '  LEVELM, FOR MESSAGES')

99006 FORMAT
     +  (/1X, A9, 'ERROR.  NUMBERS OF ALL TYPES OF UNKNOWNS MUST BE AT'
     +  /10X, 'LEAST ZERO.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS')

99007 FORMAT
     +  (/1X, A9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS')

99008 FORMAT
     +  (/1X, A9, 'ERROR.  TOTAL UNKNOWNS MUST BE POSITIVE.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  TOTAL NUMBER')

99009 FORMAT
     +  (/1X, A9, 'ERROR.  THE NUMBER OF NAMES IS WRONG.'
     + //10X, I10, '  NAMES'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  TOTAL NUMBER')

99010 FORMAT
     +  (/1X, A9, 'ERROR.  THERE ARE TOO MANY POINTS.'
     + //10X, I10, '  POINTS'
     +  /10X, I10, '  PMAX, LIMIT ON POINTS')

99011 FORMAT
     +  (/1X, A9, 'ERROR.  THE LOWER AND UPPER BOUNDS ON SOME UNKNOWNS'
     +  /10X, 'ARE OUT OF ORDER.'
     + //10X, I10, '  GROUP A UNKNOWNS (A)'
     +  /10X, I10, '  GROUP B UNKNOWNS (B)'
     +  /10X, I10, '  COMPONENTS AT POINTS (C)'
     +  /10X, I10, '  TOTAL TYPES OF UNKNOWNS'
     +  /10X, I10, '  NUMBER OF BOUNDS OUT OF ORDER'
C              123456789_  123456789_
     + //10X, '     LOWER       UPPER'
     +  /10X, '     BOUND       BOUND   UNKNOWN'
     +  /)

99012 FORMAT
     +  (/1X, A9, 'ERROR.  TWGRAB FAILS.')

99013 FORMAT
     +  (/1X, A9, 'ERROR.  ONE OR BOTH WORK SPACES ARE TOO SMALL.'
C              123456789_  123456789_
     + //25X, '   INTEGER        REAL'
C              123456789_123
     + //10X, ' PRESENT SIZE', 2I12
     +  /10X, 'REQUIRED SIZE', 2I12)

99014 FORMAT
     +  (/1X, A9, 'ERROR.  NEITHER THE INITIAL TIME EVOLUTION NOR THE'
     +  /10X, 'SEARCH FOR THE STEADY STATE IS ALLOWED.')

99015 FORMAT
     +  (/1X, A9, 'ERROR.  UNKNOWN TASK.')

99016 FORMAT
     +  (/1X, A9, 'ERROR.  UNKNOWN REPORT CODE.')

99017 FORMAT
     +  (/1X, A9, 'ERROR.  SEARCH FAILS.')

99018 FORMAT
     +  (/1X, A9, 'ERROR.  REFINE FAILS.')

99019 FORMAT
     +  (/1X, A9, 'ERROR.  EVOLVE FAILS.')

99020 FORMAT
     +  (/1X, A9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, I10, '  LABEL')

99021 FORMAT
     +  (/1X, A9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, I10, '  RETURN')

C///  EXIT.

      STOP
99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWPREP
     +  (ERROR, TEXT,
     +   A, ASIZE, BUFFER, COMPS, CONDIT, GROUPA, GROUPB, PIVOT, POINTS,
     +   RETURN)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWPREP
C
C     EVALUATE A BLOCK TRIDIAGONAL JACOBIAN MATRIX BY ONE-SIDED FINITE
C     DIFFERENCES AND REVERSE COMMUNICATION, PACK THE MATRIX INTO THE
C     LINPACK BANDED FORM, SCALE THE ROWS, AND FACTOR THE MATRIX USING
C     LINPACK'S SGBCO.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)

      CHARACTER
     +   ID*9, STRING*80
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, ABSOL, BUFFER, CONDIT, DELTA, EPS, RELAT, SUM, TEMP
      EXTERNAL
     +   TWEPS, TWGBCO, TWSQEZ
      INTEGER
     +   ASIZE, BLOCK, BLOCKS, CFIRST, CLAST, COL, COMPS, COUNT, DIAG,
     +   GROUPA, GROUPB, J, LDA, LENGTH, LINES, N, OFFSET, PIVOT,
     +   POINTS, RFIRST, RLAST, ROUTE, ROW, SKIP, TEXT, WIDTH
      INTRINSIC
     +   ABS, INT, MAX, MIN, MOD, SQRT
      LOGICAL
     +   ERROR, FOUND, MESS, RETURN

      PARAMETER (ID = 'TWPREP:  ')
      PARAMETER (LINES = 20)

      DIMENSION
     +   A(ASIZE), PIVOT(GROUPA + COMPS * POINTS + GROUPB),
     +   BUFFER(GROUPA + COMPS * POINTS + GROUPB)

      SAVE

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  EVERY-TIME INITIALIZATION.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      MESS = .FALSE.

C///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      IF (RETURN) THEN
         RETURN = .FALSE.
         GO TO (2030, 3050) ROUTE
         ERROR = .TRUE.
         GO TO 9001
      ENDIF

C///  CHECK THE ARGUMENTS.

      N = GROUPA + COMPS * POINTS + GROUPB
      ERROR = .NOT. (((0 .LT. COMPS) .EQV. (0 .LT. POINTS)) .AND.
     +   0 .LE. COMPS .AND. 0 .LE. POINTS .AND. 0 .LE. GROUPA .AND.
     +   0 .LE. GROUPB .AND. 0 .LT. N)
      IF (ERROR) GO TO 9002

      WIDTH = COMPS + MAX (COMPS, GROUPA, GROUPB) - 1
      ERROR = .NOT. ((3 * WIDTH + 2) * N .LE. ASIZE)
      IF (ERROR) GO TO 9003

C///  WRITE ALL MESSAGES.

      IF (MESS .AND. 0 .LT. TEXT) THEN
         ROUTE = 0
         GO TO 9001
      END IF

C///  FORM MACHINE EPSILON AND THE ABSOLUTE AND RELATIVE PERTURBATIONS.

      CALL TWEPS (EPS)
      ABSOL = SQRT (EPS)
      RELAT = SQRT (EPS)

C///  INITIALIZE COUNTERS AND POINTERS.

C     MAIN DIAGONAL ROW IN THE PACKING WHICH PLACES DIAGONALS IN ROWS

      DIAG = 2 * WIDTH + 1

C     PACKED ROW DIMENSION

      LDA = 3 * WIDTH + 1
      SKIP = 2 * WIDTH + 1

C     BLOCKS AND BLOCK SIZES
C     ARRAY PIVOT HOLDS BLOCK SIZES AND POINTERS TEMPORARILY

      BLOCKS = 0
      IF (0 .LT. GROUPA) THEN
         BLOCKS = BLOCKS + 1
         PIVOT(BLOCKS) = GROUPA
      END IF

      DO 1020 J = 1, POINTS
         BLOCKS = BLOCKS + 1
         PIVOT(BLOCKS) = COMPS
1020  CONTINUE

      IF (0 .LT. GROUPB) THEN
         BLOCKS = BLOCKS + 1
         PIVOT(BLOCKS) = GROUPB
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (2) INITIALIZE THE COLUMNS OF THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

C///  STORE THE EVALUATION VECTOR.

      DO 2010 J = 1, N
         A(J) = BUFFER(J)
2010  CONTINUE

C///  CLEAR THE MATRIX.

       DO 2020 J = N + 1, (3 * WIDTH + 2) * N
          A(J) = 0.0
2020  CONTINUE

C///  EVALUATE THE FUNCTION AT THE UNPERTURBED X.

C     GO TO 2030 WHEN ROUTE = 1
      ROUTE = 1
      RETURN = .TRUE.
      GO TO 99999
2030  CONTINUE

C///  PLACE THE FUNCTION VALUES IN THE MATRIX.

      CLAST = 0
      DO 2060 BLOCK = 1, BLOCKS
         CFIRST = CLAST + 1
         CLAST = CLAST + PIVOT(BLOCK)

         IF (1 .LT. BLOCK) THEN
            RFIRST = CFIRST - PIVOT(BLOCK - 1)
         ELSE
            RFIRST = CFIRST
         END IF

         IF (BLOCK .LT. BLOCKS) THEN
            RLAST = CLAST + PIVOT(BLOCK + 1)
         ELSE
            RLAST = CLAST
         END IF

         DO 2050 COL = CFIRST, CLAST
            OFFSET = N + DIAG - COL + LDA * (COL - 1)
            DO 2040 ROW = RFIRST, RLAST
               A(OFFSET + ROW) = BUFFER(ROW)
2040        CONTINUE
2050     CONTINUE
2060  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (3) FORM THE COLUMNS OF THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP OVER GROUPS OF COLUMNS.

3010  CONTINUE
      FOUND = .FALSE.

C///  RESTORE THE EVALUATION VECTOR.

      DO 3020 J = 1, N
         BUFFER(J) = A(J)
3020  CONTINUE

C///  PERTURB THE VECTOR AT INDEPENDENT POSITIONS.

      BLOCK = 1
      CFIRST = 1
3030  CONTINUE
         IF (0 .LT. PIVOT(BLOCK)) THEN
            FOUND = .TRUE.
            COL = CFIRST - 1 + PIVOT(BLOCK)
            IF (0 .LE. A(COL)) THEN
               DELTA = RELAT * A(COL) + ABSOL
            ELSE
               DELTA = RELAT * A(COL) - ABSOL
            END IF
            BUFFER(COL) = BUFFER(COL) + DELTA

            COUNT = 3
         ELSE
            COUNT = 1
         END IF

         DO 3040 J = 1, COUNT
            IF (BLOCK .EQ. 1 .AND. 0 .LT. GROUPA) THEN
               CFIRST = CFIRST + GROUPA
            ELSE IF (BLOCK .EQ. BLOCKS .AND. 0 .LT. GROUPB) THEN
               CFIRST = CFIRST + GROUPB
            ELSE
               CFIRST = CFIRST + COMPS
            END IF
            BLOCK = BLOCK + 1
3040     CONTINUE

      IF (BLOCK .LE. BLOCKS) GO TO 3030

C///  EXIT OF THE LOOP OVER GROUPS OF COLUMNS.

      IF (.NOT. FOUND) GO TO 3090

C///  EVALUATE THE FUNCTION AT THE PERTURBED VALUES.

C     GO TO 3050 WHEN ROUTE = 2
      ROUTE = 2
      RETURN = .TRUE.
      GO TO 99999
3050  CONTINUE

C///  DIFFERENCE TO FORM THE COLUMNS OF THE JACOBIAN MATRIX.

      BLOCK = 1
      CFIRST = 1
3060  CONTINUE
         IF (0 .LT. PIVOT(BLOCK)) THEN
            COL = CFIRST - 1 + PIVOT(BLOCK)
            PIVOT(BLOCK) = PIVOT(BLOCK) - 1

            IF (0 .LE. A(COL)) THEN
               DELTA = RELAT * A(COL) + ABSOL
            ELSE
               DELTA = RELAT * A(COL) - ABSOL
            END IF
            TEMP = 1.0 / DELTA
            OFFSET = N + DIAG - COL + LDA * (COL - 1)

            IF (BLOCK .EQ. 1 .AND. 0 .LT. GROUPA) THEN
               CLAST = CFIRST + GROUPA - 1
            ELSE IF (BLOCK .EQ. BLOCKS .AND. 0 .LT. GROUPB) THEN
               CLAST = CFIRST + GROUPB - 1
            ELSE
               CLAST = CFIRST + COMPS - 1
            END IF

            IF (1 .LT. BLOCK) THEN
               IF (BLOCK .EQ. 2 .AND. 0 .LT. GROUPA) THEN
                  RFIRST = CFIRST - GROUPA
               ELSE
                  RFIRST = CFIRST - COMPS
               END IF
            ELSE
               RFIRST = CFIRST
            END IF

            IF (BLOCK .LT. BLOCKS) THEN
               IF (BLOCK .EQ. BLOCKS - 1 .AND. 0 .LT. GROUPB) THEN
                  RLAST = CLAST + GROUPB
               ELSE
                  RLAST = CLAST + COMPS
               END IF
            ELSE
               RLAST = CLAST
            END IF

            DO 3070 ROW = RFIRST, RLAST
               A(OFFSET + ROW) = (BUFFER(ROW) - A(OFFSET + ROW)) * TEMP
3070        CONTINUE

            COUNT = 3
         ELSE
            COUNT = 1
         END IF

         DO 3080 J = 1, COUNT
            IF (BLOCK .EQ. 1 .AND. 0 .LT. GROUPA) THEN
               CFIRST = CFIRST + GROUPA
            ELSE IF (BLOCK .EQ. BLOCKS .AND. 0 .LT. GROUPB) THEN
               CFIRST = CFIRST + GROUPB
            ELSE
               CFIRST = CFIRST + COMPS
            END IF
            BLOCK = BLOCK + 1
3080     CONTINUE

      IF (BLOCK .LE. BLOCKS) GO TO 3060

C///  BOTTOM OF THE LOOP OVER GROUPS OF COLUMNS.

      GO TO 3010
3090  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (4) CHECK FOR ZERO COLUMNS.
C
C///////////////////////////////////////////////////////////////////////

      COUNT = 0
      DO 4020 COL = 1, N
         OFFSET = N + DIAG - COL + LDA * (COL - 1)
         SUM = 0.0
         DO 4010 ROW = MAX (COL - WIDTH, 1), MIN (COL + WIDTH, N)
            SUM = SUM + ABS (A(OFFSET + ROW))
4010     CONTINUE
         A(COL) = SUM

         IF (SUM .EQ. 0.0) COUNT = COUNT + 1
4020  CONTINUE

      ERROR = .NOT. (COUNT .EQ. 0)
      IF (ERROR) GO TO 9004

C///////////////////////////////////////////////////////////////////////
C
C     (5) SCALE THE ROWS.
C
C///////////////////////////////////////////////////////////////////////

      COUNT = 0
      DO 5030 ROW = 1, N
         OFFSET = N + DIAG + ROW
         SUM = 0.0
         DO 5010 COL = MAX (ROW - WIDTH, 1), MIN (ROW + WIDTH, N)
            SUM = SUM + ABS (A(OFFSET - COL + LDA * (COL - 1)))
5010     CONTINUE

         IF (SUM .EQ. 0.0) THEN
            COUNT = COUNT + 1
            A(ROW) = SUM
         ELSE
            TEMP = 1.0 / SUM
            A(ROW) = TEMP

            DO 5020 COL = MAX (ROW - WIDTH, 1), MIN (ROW + WIDTH, N)
               A(OFFSET - COL + LDA * (COL - 1))
     +            = A(OFFSET - COL + LDA * (COL - 1)) * TEMP
5020        CONTINUE
         ENDIF
5030  CONTINUE

      ERROR = .NOT. (COUNT .EQ. 0)
      IF (ERROR) GO TO 9005

C///////////////////////////////////////////////////////////////////////
C
C     (6) FACTOR THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

      CALL TWGBCO
     +  (A(N + 1), LDA, N, WIDTH, WIDTH, PIVOT, CONDIT, BUFFER)

      ERROR = CONDIT .EQ. 0.0
      IF (ERROR) GO TO 9006

      CONDIT = 1.0 / CONDIT

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

80001 FORMAT
     +  (10X, A)

80002 FORMAT
     +  (10X, '... MORE')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, ROUTE
      IF (.NOT. MESS) GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID,
     +   COMPS, POINTS, GROUPA, GROUPB, N
      IF (.NOT. MESS) GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID,
     +   COMPS, POINTS, GROUPA, GROUPB, N, WIDTH,
     +   (3 * WIDTH + 2) * N, ASIZE
      IF (.NOT. MESS) GO TO 99999

9004  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99004) ID, COMPS, POINTS, GROUPA, GROUPB,
     +      GROUPA + COMPS * POINTS + GROUPB, COUNT
         COUNT = 0
         DO 8010 J = 1, GROUPA + COMPS * POINTS + GROUPB
            IF (A(J) .EQ. 0.0 .OR. MESS) THEN
               COUNT = COUNT + 1
               IF (COUNT .LE. LINES) THEN
                  IF (J .LE. GROUPA) THEN
                     WRITE (STRING, '(A, I10)') 'GROUP A ', J
                  ELSE IF (J .LE. GROUPA + COMPS * POINTS) THEN
                     WRITE (STRING, '(A, I10, A, I10)')
     +                  ' COMPONENT ', MOD (J - GROUPA - 1, COMPS) + 1,
     +                  ' AT POINT ', INT ((J - GROUPA - 1) / COMPS) + 1
                  ELSE
                     WRITE (STRING, '(A, I10)')
     +                  'GROUP B ', J - GROUPA - COMPS * POINTS
                  END IF
                  CALL TWSQEZ (LENGTH, STRING)
                  WRITE (TEXT, 80001) STRING (1 : LENGTH)
               END IF
            END IF
8010     CONTINUE
         IF (LINES .LT. COUNT) WRITE (TEXT, 80002)
      END IF
      IF (.NOT. MESS) GO TO 99999

9005  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99005) ID, COMPS, POINTS, GROUPA, GROUPB,
     +      GROUPA + COMPS * POINTS + GROUPB, COUNT
         COUNT = 0
         DO 8020 J = 1, GROUPA + COMPS * POINTS + GROUPB
            IF (A(J) .EQ. 0.0 .OR. MESS) THEN
               COUNT = COUNT + 1
               IF (COUNT .LE. LINES) THEN
                  IF (J .LE. GROUPA) THEN
                     WRITE (STRING, '(A, I10)') 'GROUP A ', J
                  ELSE IF (J .LE. GROUPA + COMPS * POINTS) THEN
                     WRITE (STRING, '(A, I10, A, I10)')
     +                  ' COMPONENT ', MOD (J - GROUPA - 1, COMPS) + 1,
     +                  ' AT POINT ', INT ((J - GROUPA - 1) / COMPS) + 1
                  ELSE
                     WRITE (STRING, '(A, I10)')
     +                  'GROUP B ', J - GROUPA - COMPS * POINTS
                  END IF
                  CALL TWSQEZ (LENGTH, STRING)
                  WRITE (TEXT, 80001) STRING (1 : LENGTH)
               END IF
            END IF
8020     CONTINUE
         IF (LINES .LT. COUNT) WRITE (TEXT, 80002)
      END IF
      IF (.NOT. MESS) GO TO 99999

9006  IF (0 .LT. TEXT) WRITE (TEXT, 99006) ID
      IF (.NOT. MESS) GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, I10, '  ROUTE')

99002 FORMAT
     +  (/1X, A9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES'
     +  /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  TOTAL UNKNOWNS')

99003 FORMAT
     +  (/1X, A9, 'ERROR.  THE MATRIX SPACE IS TOO SMALL.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  MATRIX ORDER'
     +  /10X, I10, '  STRICT HALF BANDWIDTH'
     + //10X, I10, '  SPACE REQUIRED'
     +  /10X, I10, '  ASIZE, PROVIDED')

99004 FORMAT
     +  (/1X, A9, 'ERROR.  SOME COLUMNS ARE ZERO.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  TOTAL COLUMNS'
     +  /10X, I10, '  ZERO COLUMNS'
     + //10X, 'UNKNOWNS WITH ZERO COLUMNS:'
     +  /)

99005 FORMAT
     +  (/1X, A9, 'ERROR.  SOME ROWS ARE ZERO.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  TOTAL ROWS'
     +  /10X, I10, '  ZERO ROWS'
     + //10X, 'ZERO ROWS:'
     +  /)

99006 FORMAT
     +  (/1X, A9, 'ERROR.  THE JACOBIAN MATRIX IS SINGULAR.')

C///  EXIT.

      STOP
99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWSETI (ERROR, TEXT, CONTRL, VALUE)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSETI
C
C     SET A CONTROL THAT TAKES AN INTEGER VALUE.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   CONTRL*(*), ID*9, STRING*80
      EXTERNAL
     +   TWINIT, TWLAST
      INTEGER
     +   CNTRLS, COUNT, IVALUE, LENGTH, TEXT, VALUE
      LOGICAL
     +   ERROR, FOUND, MESS, LVALUE

      PARAMETER (ID = 'TWSETI:  ')
      PARAMETER (CNTRLS = 22)

      DIMENSION IVALUE(CNTRLS), LVALUE(CNTRLS)

      COMMON / TWCOMI / IVALUE
      COMMON / TWCOML / LVALUE

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      MESS = .FALSE.

      IF (MESS .AND. 0 .LT. TEXT) GO TO 9001

C///  INITIALIZE THE CONTROLS.

      CALL TWINIT (ERROR, TEXT, .FALSE.)
      IF (ERROR) GO TO 9001

C///  SET THE CONTROLS.

      COUNT = 0
      FOUND = .FALSE.

C     ADAPT

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'ADAPT') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     LEVELD

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'LEVELD') THEN
         FOUND = .TRUE.
         IVALUE(COUNT) = VALUE
      END IF

C     LEVELM

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'LEVELM') THEN
         FOUND = .TRUE.
         IVALUE(COUNT) = VALUE
      END IF

C     PADD

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'PADD') THEN
         FOUND = .TRUE.
         IVALUE(COUNT) = VALUE
         LVALUE(COUNT) = .TRUE.
      END IF

C     SSABS

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'SSABS') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     SSAGE

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'SSAGE') THEN
         FOUND = .TRUE.
         IVALUE(COUNT) = VALUE
      END IF

C     SSREL

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'SSREL') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     STEADY

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEADY') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     STEPS0

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEPS0') THEN
         FOUND = .TRUE.
         IVALUE(COUNT) = VALUE
      END IF

C     STEPS1

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEPS1') THEN
         FOUND = .TRUE.
         IVALUE(COUNT) = VALUE
      END IF

C     STEPS2

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEPS2') THEN
         FOUND = .TRUE.
         IVALUE(COUNT) = VALUE
      END IF

C     STRID0

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STRID0') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TDABS

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDABS') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TDAGE

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDAGE') THEN
         FOUND = .TRUE.
         IVALUE(COUNT) = VALUE
      END IF

C     TDEC

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDEC') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TDREL

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDREL') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TINC

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TINC') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TMAX

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TMAX') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TMIN

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TMIN') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TOLER0

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TOLER0') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TOLER1

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TOLER1') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TOLER2

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TOLER2') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

      ERROR = .NOT. (COUNT .EQ. CNTRLS)
      IF (ERROR) GO TO 9004

      ERROR = .NOT. FOUND
      IF (ERROR) GO TO 9005

C///  ERROR MESSAGES.

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID
      IF (.NOT. MESS) GO TO 99999

9002  IF (0 .LT. TEXT) THEN
         CALL TWLAST (LENGTH, CONTRL)
         IF (LENGTH .LE. 40) THEN
            STRING = CONTRL
         ELSE
            LENGTH = 40
            STRING = CONTRL (1 : 37) // '...'
         END IF
         WRITE (TEXT, 99002) ID, STRING (1 : LENGTH)
      END IF
      IF (.NOT. MESS) GO TO 99999

9003  IF (0 .LT. TEXT) THEN
         CALL TWLAST (LENGTH, CONTRL)
         IF (LENGTH .LE. 40) THEN
            STRING = CONTRL
         ELSE
            LENGTH = 40
            STRING = CONTRL (1 : 37) // '...'
         END IF
         WRITE (TEXT, 99003) ID, STRING (1 : LENGTH)
      END IF
      IF (.NOT. MESS) GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, CNTRLS, COUNT
      IF (.NOT. MESS) GO TO 99999

9005  IF (0 .LT. TEXT) THEN
         CALL TWLAST (LENGTH, CONTRL)
         IF (LENGTH .LE. 40) THEN
            STRING = CONTRL
         ELSE
            LENGTH = 40
            STRING = CONTRL (1 : 37) // '...'
         END IF
         WRITE (TEXT, 99005) ID, STRING (1 : LENGTH)
      END IF
      IF (.NOT. MESS) GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  TWINIT FAILS.')

99002 FORMAT
     +  (/1X, A9, 'ERROR.  THE CONTROL TAKES A LOGICAL VALUE WHICH'
     +  /10X, 'MUST BE SET USING TWSETL.'
     + //10X, '     CONTROL:  ', A)

99003 FORMAT
     +  (/1X, A9, 'ERROR.  THE CONTROL TAKES A REAL VALUE WHICH MUST BE'
     +  /10X, 'SET USING TWSETR.'
     + //10X, '     CONTROL:  ', A)

99004 FORMAT
     +  (/1X, A9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.'
     + //10X, I10, '  CONTROLS'
     +  /10X, I10, '  COUNTED')

99005 FORMAT
     +  (/1X, A9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.'
     + //10X, '     CONTROL:  ', A)

C///  EXIT.

      STOP
99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWSETL (ERROR, TEXT, CONTRL, VALUE)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSETL
C
C     SET A CONTROL THAT TAKES A LOGICAL VALUE.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   CONTRL*(*), ID*9, STRING*80
      EXTERNAL
     +   TWINIT, TWLAST
      INTEGER
     +   CNTRLS, COUNT, LENGTH, TEXT
      LOGICAL
     +   ERROR, FOUND, LVALUE, MESS, VALUE

      PARAMETER (ID = 'TWSETL:  ')
      PARAMETER (CNTRLS = 22)

      DIMENSION LVALUE(CNTRLS)

      COMMON / TWCOML / LVALUE

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      MESS = .FALSE.

      IF (MESS .AND. 0 .LT. TEXT) GO TO 9001

C///  INITIALIZE THE CONTROLS.

      CALL TWINIT (ERROR, TEXT, .FALSE.)
      IF (ERROR) GO TO 9001

C///  SET THE CONTROLS.

      COUNT = 0
      FOUND = .FALSE.

C     ADAPT

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'ADAPT') THEN
         FOUND = .TRUE.
         LVALUE(COUNT)= VALUE
      END IF

C     LEVELD

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'LEVELD') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     LEVELM

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'LEVELM') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     PADD

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'PADD') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     SSABS

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'SSABS') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     SSAGE

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'SSAGE') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     SSREL

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'SSREL') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     STEADY

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEADY') THEN
         FOUND = .TRUE.
         LVALUE(COUNT)= VALUE
      END IF

C     STEPS0

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEPS0') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     STEPS1

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEPS1') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     STEPS2

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEPS2') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     STRID0

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STRID0') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TDABS

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDABS') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TDAGE

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDAGE') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     TDEC

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDEC') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TDREL

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDREL') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TINC

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TINC') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TMAX

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TMAX') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TMIN

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TMIN') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TOLER0

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TOLER0') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TOLER1

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TOLER1') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TOLER2

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TOLER2') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

      ERROR = .NOT. (COUNT .EQ. CNTRLS)
      IF (ERROR) GO TO 9004

      ERROR = .NOT. FOUND
      IF (ERROR) GO TO 9005

C///  ERROR MESSAGES.

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID
      IF (.NOT. MESS) GO TO 99999

9002  IF (0 .LT. TEXT) THEN
         CALL TWLAST (LENGTH, CONTRL)
         IF (LENGTH .LE. 40) THEN
            STRING = CONTRL
         ELSE
            LENGTH = 40
            STRING = CONTRL (1 : 37) // '...'
         END IF
         WRITE (TEXT, 99002) ID, STRING (1 : LENGTH)
      END IF
      IF (.NOT. MESS) GO TO 99999

9003  IF (0 .LT. TEXT) THEN
         CALL TWLAST (LENGTH, CONTRL)
         IF (LENGTH .LE. 40) THEN
            STRING = CONTRL
         ELSE
            LENGTH = 40
            STRING = CONTRL (1 : 37) // '...'
         END IF
         WRITE (TEXT, 99003) ID, STRING (1 : LENGTH)
      END IF
      IF (.NOT. MESS) GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, CNTRLS, COUNT
      IF (.NOT. MESS) GO TO 99999

9005  IF (0 .LT. TEXT) THEN
         CALL TWLAST (LENGTH, CONTRL)
         IF (LENGTH .LE. 40) THEN
            STRING = CONTRL
         ELSE
            LENGTH = 40
            STRING = CONTRL (1 : 37) // '...'
         END IF
         WRITE (TEXT, 99005) ID, STRING (1 : LENGTH)
      END IF
      IF (.NOT. MESS) GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  TWINIT FAILS.')

99002 FORMAT
     +  (/1X, A9, 'ERROR.  THE CONTROL TAKES AN INTEGER VALUE WHICH'
     +  /10X, 'MUST BE SET USING TWSETI.'
     + //10X, '     CONTROL:  ', A)

99003 FORMAT
     +  (/1X, A9, 'ERROR.  THE CONTROL TAKES A REAL VALUE WHICH MUST BE'
     +  /10X, 'SET USING TWSETR.'
     + //10X, '     CONTROL:  ', A)

99004 FORMAT
     +  (/1X, A9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.'
     + //10X, I10, '  CONTROLS'
     +  /10X, I10, '  COUNTED')

99005 FORMAT
     +  (/1X, A9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.'
     + //10X, '     CONTROL:  ', A)

C///  EXIT.

      STOP
99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWSETR (ERROR, TEXT, CONTRL, VALUE)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSETR
C
C     SET A CONTROL THAT TAKES A REAL VALUE.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   CONTRL*(*), ID*9, STRING*80
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   RVALUE, VALUE
      EXTERNAL
     +   TWINIT, TWLAST
      INTEGER
     +   CNTRLS, COUNT, LENGTH, TEXT
      LOGICAL
     +   ERROR, FOUND, MESS

      PARAMETER (ID = 'TWSETR:  ')
      PARAMETER (CNTRLS = 22)

      DIMENSION RVALUE(CNTRLS)

      COMMON / TWCOMR / RVALUE

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      MESS = .FALSE.

      IF (MESS .AND. 0 .LT. TEXT) GO TO 9001

C///  INITIALIZE THE CONTROLS.

      CALL TWINIT (ERROR, TEXT, .FALSE.)
      IF (ERROR) GO TO 9001

C///  SET THE CONTROLS.

      COUNT = 0
      FOUND = .FALSE.

C     ADAPT

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'ADAPT') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     LEVELD

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'LEVELD') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     LEVELM

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'LEVELM') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     PADD

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'PADD') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     SSABS

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'SSABS') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

C     SSAGE

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'SSAGE') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     SSREL

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'SSREL') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

C     STEADY

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEADY') THEN
         ERROR = .TRUE.
         GO TO 9002
      END IF

C     STEPS0

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEPS0') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     STEPS1

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEPS1') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     STEPS2

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STEPS2') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     STRID0

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'STRID0') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

C     TDABS

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDABS') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

C     TDAGE

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDAGE') THEN
         ERROR = .TRUE.
         GO TO 9003
      END IF

C     TDEC

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDEC') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

C     TDREL

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TDREL') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

C     TINC

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TINC') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

C     TMAX

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TMAX') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

C     TMIN

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TMIN') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

C     TOLER0

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TOLER0') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

C     TOLER1

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TOLER1') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

C     TOLER2

      COUNT = COUNT + 1
      IF (CONTRL .EQ. 'TOLER2') THEN
         FOUND = .TRUE.
         RVALUE(COUNT) = VALUE
      END IF

      ERROR = .NOT. (COUNT .EQ. CNTRLS)
      IF (ERROR) GO TO 9004

      ERROR = .NOT. FOUND
      IF (ERROR) GO TO 9005

C///  ERROR MESSAGES.

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID
      IF (.NOT. MESS) GO TO 99999

9002  IF (0 .LT. TEXT) THEN
         CALL TWLAST (LENGTH, CONTRL)
         IF (LENGTH .LE. 40) THEN
            STRING = CONTRL
         ELSE
            LENGTH = 40
            STRING = CONTRL (1 : 37) // '...'
         END IF
         WRITE (TEXT, 99002) ID, STRING (1 : LENGTH)
      END IF
      IF (.NOT. MESS) GO TO 99999

9003  IF (0 .LT. TEXT) THEN
         CALL TWLAST (LENGTH, CONTRL)
         IF (LENGTH .LE. 40) THEN
            STRING = CONTRL
         ELSE
            LENGTH = 40
            STRING = CONTRL (1 : 37) // '...'
         END IF
         WRITE (TEXT, 99003) ID, STRING (1 : LENGTH)
      END IF
      IF (.NOT. MESS) GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, CNTRLS, COUNT
      IF (.NOT. MESS) GO TO 99999

9005  IF (0 .LT. TEXT) THEN
         CALL TWLAST (LENGTH, CONTRL)
         IF (LENGTH .LE. 40) THEN
            STRING = CONTRL
         ELSE
            LENGTH = 40
            STRING = CONTRL (1 : 37) // '...'
         END IF
         WRITE (TEXT, 99005) ID, STRING (1 : LENGTH)
      END IF
      IF (.NOT. MESS) GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  TWINIT FAILS.')

99002 FORMAT
     +  (/1X, A9, 'ERROR.  THE CONTROL TAKES A LOGICAL VALUE WHICH'
     +  /10X, 'MUST BE SET USING TWSETL.'
     + //10X, '     CONTROL:  ', A)

99003 FORMAT
     +  (/1X, A9, 'ERROR.  THE CONTROL TAKES AN INTEGER VALUE WHICH'
     +  /10X, 'MUST BE SET USING TWSETI.'
     + //10X, '     CONTROL:  ', A)

99004 FORMAT
     +  (/1X, A9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.'
     + //10X, I10, '  CONTROLS'
     +  /10X, I10, '  COUNTED')

99005 FORMAT
     +  (/1X, A9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.'
     + //10X, '     CONTROL:  ', A)

C///  EXIT.

      STOP
99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWSHOW
     +  (ERROR, TEXT,
     +   BUFFER, COMPS, GRID, GROUPA, GROUPB, POINTS, X)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSHOW
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)

      CHARACTER
     +   ID*9, STRING*80, TITLE*80
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   BUFFER, X
      EXTERNAL
     +   TWSQEZ
      INTEGER
     +   COLS, COMP, COMPS, COUNT, FIRST, GROUPA, GROUPB, GROUPS, J,
     +   LAST, LENGTH, POINT, POINTS, TEXT
      INTRINSIC
     +   MIN
      LOGICAL
     +   ERROR, GRID, MESS

      PARAMETER (ID = 'TWSHOW:  ')

      DIMENSION
     +   BUFFER(GROUPA + COMPS * POINTS + GROUPB), TITLE(6), X(*)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      MESS = .FALSE.

      IF (MESS .AND. 0 .LT. TEXT) GO TO 9001

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (((0 .LT. COMPS) .EQV. (0 .LT. POINTS)) .AND.
     +   0 .LE. COMPS .AND. 0 .LE. POINTS .AND. 0 .LE. GROUPA .AND.
     +   0 .LE. GROUPB .AND. 0 .LT. GROUPA + COMPS * POINTS + GROUPB)
      IF (ERROR) GO TO 9001

C///  COUNT THE GROUPS.

      GROUPS = 0
      IF (0 .LT. GROUPA) GROUPS = GROUPS + 1
      IF (0 .LT. GROUPB) GROUPS = GROUPS + 1
      IF (0 .LT. COMPS .AND. 0 .LT. POINTS) GROUPS = GROUPS + 1

C///  CHOOSE NUMBER OF DATA COLUMNS.

      IF (GRID) THEN
         COLS = 5
      ELSE
         COLS = 6
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (2) PRINT THE GROUPED DATA.
C
C///////////////////////////////////////////////////////////////////////

      IF (0 .LT. TEXT) THEN

      IF (0 .LT. GROUPA) THEN
         IF (1 .LT. GROUPS) WRITE (TEXT, 10001) 'GROUP A UNKNOWNS'
         WRITE (TEXT, 10002) (J, BUFFER(J), J = 1, GROUPA)
      END IF

      IF (0 .LT. GROUPB) THEN
         IF (1 .LT. GROUPS) WRITE (TEXT, 10001) 'GROUP B UNKNOWNS'
         WRITE (TEXT, 10002)
     +      (J, BUFFER(GROUPA + COMPS * POINTS + J), J = 1, GROUPB)
      END IF

C///////////////////////////////////////////////////////////////////////
C
C     (2) PRINT THE COMPONENTS AT POINTS.
C
C///////////////////////////////////////////////////////////////////////

      IF (0 .LT. COMPS .AND. 0 .LT. POINTS) THEN
         IF (1 .LT. GROUPS) WRITE (TEXT, 10001) 'COMPONENTS AT POINTS'

         DO 2030 FIRST = 1, COMPS, COLS
            COUNT = 0
            LAST = MIN (FIRST + COLS - 1, COMPS)
            DO 2010 COMP = FIRST, LAST
               COUNT = COUNT + 1
               TITLE(COUNT) = ' '
               WRITE (STRING, '(A5, I5)') 'COMP ', COMP
               CALL TWSQEZ (LENGTH, STRING)
               TITLE(COUNT) (11 - LENGTH : 10) = STRING
2010        CONTINUE

            IF (GRID) THEN
               WRITE (TEXT, 10003)
     +            'GRID POINT', (TITLE(J), J = 1, COUNT)
            ELSE
               WRITE (TEXT, 10003) (TITLE(J), J = 1, COUNT)
            END IF

            IF (COUNT .EQ. COLS) THEN
               IF (GRID) THEN
                  WRITE (TEXT, 10004) (POINT, X(POINT),
     +               (BUFFER(GROUPA + COMP + COMPS * (POINT - 1)),
     +               COMP = FIRST, LAST), POINT = 1, POINTS)
               ELSE
                  WRITE (TEXT, 10005) (POINT,
     +               (BUFFER(GROUPA + COMP + COMPS * (POINT - 1)),
     +               COMP = FIRST, LAST), POINT = 1, POINTS)
               END IF
            ELSE
               DO 2020 POINT = 1, POINTS
                  IF (GRID) THEN
                     WRITE (TEXT, 10004) POINT, X(POINT),
     +                  (BUFFER(GROUPA + COMP + COMPS * (POINT - 1)),
     +                  COMP = FIRST, LAST)
                  ELSE
                     WRITE (TEXT, 10005) POINT,
     +                  (BUFFER(GROUPA + COMP + COMPS * (POINT - 1)),
     +                  COMP = FIRST, LAST)
                  END IF
2020           CONTINUE
            END IF
2030     CONTINUE
      END IF

      END IF

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +  (/10X, A)

10002 FORMAT
     + (/(10X, 4(I3, '> ', 1PE10.3)))

10003 FORMAT
     +  (/14X, 6(1X, A10))

10004 FORMAT
     +  (10X, 0P, I3, '>', F11.6, 1P, 5E11.3)

10005 FORMAT
     +  (10X, 0P, I3, '>', 1P, 6E11.3)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID,
     +   COMPS, POINTS, GROUPA, GROUPB, GROUPA + COMPS * POINTS + GROUPB
      IF (.NOT. MESS) GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES'
     +  /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  TOTAL UNKNOWNS')
      IF (.NOT. MESS) GO TO 99999

C///  EXIT.

      STOP
99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWSOLV
     +  (ERROR, TEXT,
     +   A, ASIZE, BUFFER, COMPS, GROUPA, GROUPB, PIVOT, POINTS)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSOLV
C
C     SOLVE A SYSTEM OF LINEAR EQUATIONS USING THE MATRIX PREPARED BY
C     TWPREP.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)

      CHARACTER
     +   ID*9
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   A, BUFFER
      EXTERNAL
     +   TWGBSL
      INTEGER
     +   ASIZE, COMPS, GROUPA, GROUPB, J, N, PIVOT, POINTS, TEXT, WIDTH
      INTRINSIC
     +   MAX
      LOGICAL
     +   ERROR, MESS

      PARAMETER (ID = 'TWSOLV:  ')

      DIMENSION
     +   A(ASIZE), BUFFER(GROUPA + COMPS * POINTS + GROUPB),
     +   PIVOT(GROUPA + COMPS * POINTS + GROUPB)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      MESS = .FALSE.

      IF (MESS .AND. 0 .LT. TEXT) GO TO 9001

C///  CHECK THE ARGUMENTS.

      N = GROUPA + COMPS * POINTS + GROUPB
      ERROR = .NOT. (((0 .LT. COMPS) .EQV. (0 .LT. POINTS)) .AND.
     +   0 .LE. COMPS .AND. 0 .LE. POINTS .AND. 0 .LE. GROUPA .AND.
     +   0 .LE. GROUPB .AND. 0 .LT. N)
      IF (ERROR) GO TO 9001

      WIDTH = COMPS + MAX (COMPS, GROUPA, GROUPB) - 1
      ERROR = .NOT. ((3 * WIDTH + 2) * N .LE. ASIZE)
      IF (ERROR) GO TO 9002

C///////////////////////////////////////////////////////////////////////
C
C     (2) SCALE AND SOLVE THE EQUATIONS.
C
C///////////////////////////////////////////////////////////////////////

      DO 2010 J = 1, N
         BUFFER(J) = BUFFER(J) * A(J)
2010  CONTINUE

      CALL TWGBSL
     +  (A(N + 1), 3 * WIDTH + 1, N, WIDTH, WIDTH, PIVOT, BUFFER)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID,
     +   COMPS, POINTS, GROUPA, GROUPB, N
      IF (.NOT. MESS) GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID,
     +   COMPS, POINTS, GROUPA, GROUPB, N, WIDTH,
     +   (3 * WIDTH + 2) * N, ASIZE
      IF (.NOT. MESS) GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES'
     +  /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  TOTAL UNKNOWNS')

99002 FORMAT
     +  (/1X, A9, 'ERROR.  THE MATRIX SPACE IS TOO SMALL.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS'
     +  /10X, I10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, I10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, I10, '  MATRIX ORDER'
     +  /10X, I10, '  STRICT HALF BANDWIDTH'
     + //10X, I10, '  SPACE EXPECTED'
     +  /10X, I10, '  ASIZE, PROVIDED')

C///  EXIT.

      STOP
99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWSQEZ (LENGTH, STRING)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSQEZ
C
C     SQUEEZE LEADING BLANKS AND MULTIPLE BLANKS FROM A CHARACTER
C     STRING.  RETURN THE LENGTH OF THE SQUEEZED STRING.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   CHAR*1, STRING*(*)
      INTEGER
     +   J, LENGTH
      INTRINSIC
     +   LEN
      LOGICAL
     +   BLANK

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
      SUBROUTINE TWTIME (TIMER)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWTIME
C
C     OBTAIN COMPUTING TIME IN SECONDS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   TIMER

C     THE SYSTEM LIBRARY BASELIB CONTAINS SUBROUTINE IZM00.
C*****COMPUTING TIME > CRI (CRAY) CTSS (LIVERMORE)
C      EXTERNAL IZM00
C      INTEGER FLAG, IZM00, VALUE
C      REAL TEMP
C      DIMENSION VALUE(5)
C
C      FLAG = IZM00 (VALUE)
C      IF (FLAG .EQ. 0) THEN
C         TEMP = REAL (VALUE(1)) / 1.0E6
C      ELSE
C         TEMP = 0.0
C      END IF
C*****END COMPUTING TIME > CRI (CRAY) CTSS (LIVERMORE)

C     THE SYSTEM LIBRARY CFTLIB CONTAINS SUBROUTINE ISECOND.
C*****COMPUTING TIME > CRI (CRAY) CTSS (LOS ALAMOS)
C      EXTERNAL ISECOND
C      INTEGER VALUE
C      REAL TEMP
C
C      CALL ISECOND (VALUE)
C      TEMP = REAL (VALUE) / 1.0E6
C*****END COMPUTING TIME > CRI (CRAY) CTSS (LOS ALAMOS)

C*****COMPUTING TIME > CRI (CRAY) UNICOS
C      REAL SECOND, TEMP
C      TEMP = SECOND ()
C*****END COMPUTING TIME > CRI (CRAY) UNICOS

C     THIS IS A SYSTEM SERVICE CALL FROM FORTRAN.
C*****COMPUTING TIME > DEC (VAX) VMS
C      INTEGER*2 LIST
C      INTEGER*4 ADDRESS, TICS
C      REAL*4 TEMP
C
C      DIMENSION LIST(6)
C
C      EQUIVALENCE (LIST(3), ADDRESS)
C
C      DATA LIST(1), LIST(2), LIST(5), LIST(6) / 4, '0407'X, 0, 0 /
C
C      ADDRESS = %LOC (TICS)
C      CALL SYS$GETJPIW (, , , LIST, , , )
C      TEMP = REAL (TICS) / 1.0E2
C*****END COMPUTING TIME > DEC (VAX) VMS

C*****COMPUTING TIME > generic unix etime
C      EXTERNAL ETIME
C      REAL ETIME, LIST, TEMP
C      DIMENSION LIST(2)
C
C      TEMP = ETIME (LIST)
C*****END COMPUTING TIME > generic unix etime

C*****COMPUTING TIME > IBM (RISC System/6000) AIX
C      INTEGER MCLOCK
C      REAL TEMP
C
C      TEMP = REAL (MCLOCK()) * 0.01
C*****END COMPUTING TIME > IBM (RISC System/6000) AIX

C*****COMPUTING TIME > none
      REAL TEMP

      TEMP = 0.0
C*****END COMPUTING TIME > none

C*****COMPUTING TIME > SUN (SPARCstation) SunOS
C      EXTERNAL ETIME
C      REAL ETIME, LIST, TEMP
C      DIMENSION LIST(2)
C
C      TEMP = ETIME (LIST)
C*****END COMPUTING TIME > SUN (SPARCstation) SunOS

      TIMER = TEMP

      RETURN
      END
