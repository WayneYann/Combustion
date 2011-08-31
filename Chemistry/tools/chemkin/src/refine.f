C  CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      SUBROUTINE REFINE
     +  (ERROR, TEXT,
     +   ACTIVE, BUFFER, COMPS, LEVELD, LEVELM, MARK, NEWX, PADD, PMAX,
     +   POINTS, RATIO, RATIO1, RATIO2, SIGNAL, SUCCES, TOLER0, TOLER1,
     +   TOLER2, U, VARY1, VARY2, WEIGHT, X)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     REFINE
C
C     PERFORM AUTOMATIC GRID SELECTION.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - Z)
      CHARACTER
     +   ID*9, SIGNAL*(*), WORD*80
C*****PRECISION > DOUBLE
      DOUBLE PRECISION
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   BUFFER, DIFFER, LEFT, LENGTH, LOWER, MAXMAG, MEAN, RANGE,
     +   RATIO, RATIO1, RATIO2, RIGHT, TEMP, TEMP1, TEMP2, TOLER0,
     +   TOLER1, TOLER2, U, UPPER, X
      EXTERNAL
     +   TWCOPY
      INTEGER
     +   ACT, COMPS, COUNT, FORMER, ITEMP, J, K, LEAST, LEVELD, LEVELM,
     +   MORE, MOST, NEW, OLD, PADD, PMAX, POINTS, ROUTE, SIGNIF, TEXT,
     +   TOTAL, VARY1, VARY2, WEIGHT
      INTRINSIC
     +   ABS, MAX, MIN
      LOGICAL
     +   ACTIVE, ERROR, MARK, MESS, NEWX, SUCCES

      PARAMETER (ID = 'REFINE:  ')

      DIMENSION
     +   ACTIVE(COMPS), BUFFER(COMPS * PMAX), MARK(PMAX), RATIO(2),
     +   RATIO1(PMAX), RATIO2(PMAX), U(COMPS, PMAX), VARY1(PMAX),
     +   VARY2(PMAX), WEIGHT(PMAX), X(PMAX)

C///  SAVE LOCAL VARIABLES DURING RETURNS FOR REVERSE COMMUNICATION.

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
      NEWX = .FALSE.
      SUCCES = .FALSE.

C///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      IF (SIGNAL .NE. ' ') THEN
         GO TO (4060, 5040) ROUTE
         ERROR = .TRUE.
         GO TO 9001
      END IF

C///  WRITE ALL MESSAGES.

      IF (MESS .AND. 0 .LT. TEXT) THEN
         ROUTE = 0

         WRITE (TEXT, 10001) ID
         WRITE (TEXT, 10002) ID
         WRITE (TEXT, 10004) ID
         WRITE (TEXT, 10005) ID
         WRITE (TEXT, 10010) ID

         GO TO 9001
      END IF

C///  LEVELM PRINTING.

      IF (0 .LT. LEVELM .AND. 0 .LT. TEXT) WRITE (TEXT, 10001) ID

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (1 .LE. COMPS .AND. 2 .LE. POINTS)
      IF (ERROR) GO TO 9002

      ERROR = .NOT. (0 .LE. PADD)
      IF (ERROR) GO TO 9003

      ERROR = .NOT. (POINTS .LE. PMAX)
      IF (ERROR) GO TO 9004

      COUNT = 0
      DO 1010 J = 1, COMPS
         IF (ACTIVE(J)) COUNT = COUNT + 1
1010  CONTINUE
      ERROR = .NOT. (1 .LE. COUNT)
      IF (ERROR) GO TO 9005

      ERROR = .NOT. (0.0 .LE. TOLER0)
      IF (ERROR) GO TO 9006

      ERROR = .NOT. (0.0 .LE. TOLER1 .AND. TOLER1 .LE. 1.0
     +         .AND. 0.0 .LE. TOLER2 .AND. TOLER2 .LE. 1.0)
      IF (ERROR) GO TO 9007

      COUNT = 0
      DO 1020 K = 1, POINTS - 1
         IF (X(K) .LT. X(K + 1)) COUNT = COUNT + 1
1020  CONTINUE
      ERROR = .NOT. (COUNT .EQ. 0 .OR. COUNT .EQ. POINTS - 1)
      IF (ERROR) GO TO 9008

C///////////////////////////////////////////////////////////////////////
C
C     AT EACH INTERVAL, COUNT THE ACTIVE, SIGNIFICANT COMPONENTS THAT
C     VARY TOO GREATLY.
C
C///////////////////////////////////////////////////////////////////////

      ACT = 0
      SIGNIF = 0

      DO 2010 K = 1, POINTS
         MARK(K) = .FALSE.
         RATIO1(K) = 0.0
         RATIO2(K) = 0.0
         VARY1(K) = 0
         VARY2(K) = 0
2010  CONTINUE

C///  TOP OF THE LOOP OVER THE COMPONENTS.

      DO 2060 J = 1, COMPS
         IF (ACTIVE(J)) THEN
            ACT = ACT + 1

C///  FIND THE RANGE AND MAXIMUM MAGNITUDE OF THE COMPONENT.

            LOWER = U(J, 1)
            UPPER = U(J, 1)
            DO 2020 K = 2, POINTS
               LOWER = MIN (LOWER, U(J, K))
               UPPER = MAX (UPPER, U(J, K))
2020        CONTINUE
            RANGE = UPPER - LOWER
            MAXMAG = MAX (ABS (LOWER), ABS (UPPER))

C///  DECIDE WHETHER THE COMPONENT IS SIGNIFICANT.

            IF (ABS (RANGE) .GT. MAX (TOLER0, TOLER0 * MAXMAG)) THEN
               SIGNIF = SIGNIF + 1

C///  AT EACH INTERVAL, SEE WHETHER THE COMPONENT'S CHANGE EXCEEDS SOME
C///  FRACTION OF THE COMPONENT'S GLOBAL CHANGE.

               DO 2030 K = 1, POINTS - 1
                  DIFFER = ABS (U(J, K + 1) - U(J, K))
                  IF (0.0 .LT. RANGE)
     +               RATIO1(K) = MAX (RATIO1(K), DIFFER / RANGE)
                  IF (TOLER1 * RANGE .LT. DIFFER)
     +               VARY1(K) = VARY1(K) + 1
2030           CONTINUE

C///  FIND THE GLOBAL CHANGE OF THE COMPONENT'S DERIVATIVE.

               TEMP = (U(J, 2) - U(J, 1)) / (X(2) - X(1))
               LOWER = TEMP
               UPPER = TEMP
               DO 2040 K = 2, POINTS - 1
                  TEMP = (U(J, K + 1) - U(J, K))
     +                 / (X(K + 1) - X(K))
                  LOWER = MIN (LOWER, TEMP)
                  UPPER = MAX (UPPER, TEMP)
2040           CONTINUE
               RANGE = UPPER - LOWER

C///  AT EACH INTERIOR POINT, SEE WHETHER THE DERIVATIVE'S CHANGE
C///  EXCEEDS SOME FRACTION OF THE DERIVATIVE'S GLOBAL CHANGE.

               RIGHT = (U(J, 2) - U(J, 1)) / (X(2) - X(1))
               DO 2050 K = 2, POINTS - 1
                  LEFT = RIGHT
                  RIGHT = (U(J, K + 1) - U(J, K)) / (X(K + 1) - X(K))
                  DIFFER = ABS (LEFT - RIGHT)
                  IF (0.0 .LT. RANGE)
     +               RATIO2(K) = MAX (RATIO2(K), DIFFER / RANGE)
                  IF (TOLER2 * RANGE .LT. DIFFER)
     +               VARY2(K) = VARY2(K) + 1
2050           CONTINUE

C///  BOTTOM OF THE LOOP OVER THE COMPONENTS.

            END IF
         END IF
2060  CONTINUE

C///  SAVE THE MAXIMUM RATIOS.

      RATIO(1) = 0.0
      DO 2070 K = 1, POINTS - 1
         RATIO(1) = MAX (RATIO(1), RATIO1(K))
2070  CONTINUE

      RATIO(2) = 0.0
      DO 2080 K = 2, POINTS - 1
         RATIO(2) = MAX (RATIO(2), RATIO2(K))
2080  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     SELECT THE INTERVALS TO HALVE.
C
C///////////////////////////////////////////////////////////////////////

C///  WEIGHT THE INTERVALS IN WHICH VARIATIONS THAT ARE TOO LARGE OCCUR.

      MOST = 0
      DO 3010 K = 1, POINTS - 1
         WEIGHT(K) = VARY1(K)
         IF (1 .LT. K) WEIGHT(K) = WEIGHT(K) + VARY2(K)
         IF (K .LT. POINTS - 1) WEIGHT(K) = WEIGHT(K) + VARY2(K + 1)
         IF (0 .LT. WEIGHT(K)) MOST = MOST + 1
3010  CONTINUE

C///  SORT THE WEIGHTS.

      DO 3030 K = 1, POINTS - 1
         DO 3020 J = K + 1, POINTS - 1
            IF (WEIGHT(J) .GT. WEIGHT(K)) THEN
               ITEMP = WEIGHT(J)
               WEIGHT(J) = WEIGHT(K)
               WEIGHT(K) = ITEMP
            END IF
3020     CONTINUE
         IF (WEIGHT(K) .EQ. 0) GO TO 3040
3030  CONTINUE
3040  CONTINUE

C///  FIND THE LEAST WEIGHT OF INTERVALS TO HALVE.

      MORE = MAX (0, MIN (MOST, PADD, PMAX - POINTS))
      IF (0 .LT. MORE) THEN
         LEAST = WEIGHT(MORE)
      ELSE
         LEAST = 1 + WEIGHT(1)
      END IF

C///  RECONSTRUCT THE WEIGHTS.

      DO 3050 K = 1, POINTS - 1
         WEIGHT(K) = VARY1(K)
         IF (1 .LT. K) WEIGHT(K) = WEIGHT(K) + VARY2(K)
         IF (K .LT. POINTS - 1) WEIGHT(K) = WEIGHT(K) + VARY2(K + 1)
3050  CONTINUE

C///  MARK THE INTERVALS TO HALVE.

      COUNT = 0
      DO 3060 K = 1, POINTS - 1
         IF (COUNT .LT. MORE .AND. LEAST .LE. WEIGHT(K)) THEN
            COUNT = COUNT + 1
            MARK(K) = .TRUE.
         END IF
3060  CONTINUE

c HACK  if one point is marked, mark them all
c      if (COUNT .GT. 0) then
c         COUNT = 0
c         DO K = 1, POINTS - 1
c            COUNT = COUNT + 1
c            MARK(K) = .TRUE.
c         ENDDO
c      endif

      MORE = COUNT

C///////////////////////////////////////////////////////////////////////
C
C     HALVE THE INTERVALS, IF ANY.
C
C///////////////////////////////////////////////////////////////////////

C///  FORM THE TOTAL NUMBER OF POINTS IN THE NEW GRID.

      TOTAL = POINTS + MORE
      IF (0 .EQ. MORE) GO TO 4070

C///  TOP OF THE BLOCK TO CREATE THE NEW GRID.  CHECK THE ORDER.

      COUNT = 0
      LENGTH = ABS (X(POINTS) - X(1))
      DO 4010 K = 1, POINTS - 1
         IF (MARK(K)) THEN
            MEAN = 0.5 * (X(K) + X(K + 1))
            IF (.NOT. ((X(K) .LT. MEAN .AND. MEAN .LT. X(K + 1))
     +         .OR. (X(K + 1) .LT. MEAN .AND. MEAN .LT. X(K))))
     +         COUNT = COUNT + 1
         END IF
4010  CONTINUE
      ERROR = .NOT. (COUNT .EQ. 0)
      IF (ERROR) GO TO 9009

C///  ADD THE NEW POINTS.  INTERPOLATE X AND THE BOUNDS.

      NEW = TOTAL
      DO 4040 OLD = POINTS, 2, - 1
         X(NEW) = X(OLD)
         DO 4020 J = 1, COMPS
            U(J, NEW) = U(J, OLD)
4020     CONTINUE
         NEW = NEW - 1

         IF (MARK(OLD - 1)) THEN
            X(NEW) = 0.5 * (X(OLD) + X(OLD - 1))
            DO 4030 J = 1, COMPS
               U(J, NEW) = 0.5 * (U(J, OLD) + U(J, OLD - 1))
4030        CONTINUE
            NEW = NEW - 1
         END IF
4040  CONTINUE

C///  MARK THE NEW POINTS.

      NEW = TOTAL
      DO 4050 OLD = POINTS, 2, - 1
         MARK(NEW) = .FALSE.
         NEW = NEW - 1
         IF (MARK(OLD - 1)) THEN
            MARK(NEW) = .TRUE.
            NEW = NEW - 1
         END IF
4050  CONTINUE
      MARK(NEW) = .FALSE.

C///  UPDATE THE NUMBER OF POINTS.

      FORMER = POINTS
      POINTS = TOTAL

C///  ALLOW THE USER TO UPDATE THE SOLUTION.

      CALL TWCOPY (COMPS * POINTS, U, BUFFER)
      SIGNAL = 'UPDATE'
C     GO TO 4060 WHEN ROUTE = 1
      ROUTE = 1
      GO TO 99999
4060  CONTINUE
      SIGNAL = ' '
      CALL TWCOPY (COMPS * POINTS, BUFFER, U)

C///  BOTTOM OF THE BLOCK TO CREATE A NEW GRID.

4070  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  PRINT.

      IF (0 .LT. LEVELM .AND. 0 .LT. TEXT) THEN
         TEMP1 = RATIO1(1)
         DO 5010 K = 2, FORMER - 1
            TEMP1 = MAX (TEMP1, RATIO1(K))
5010     CONTINUE

         TEMP2 = RATIO2(2)
         DO 5020 K = 3, FORMER - 1
            TEMP2 = MAX (TEMP2, RATIO2(K))
5020     CONTINUE

         IF (SIGNIF .EQ. 0) THEN
            WRITE (TEXT, 10002) ID
         ELSE
            WRITE (TEXT, 10003) TEMP1, TEMP2, TOLER1, TOLER2
            IF (MOST .EQ. 0) THEN
               WRITE (TEXT, 10004) ID
            ELSE IF (MORE .EQ. 0) THEN
               WRITE (TEXT, 10005) ID
            ELSE
               WRITE (TEXT, 10006)

               OLD = 0
               DO 5030 K = 1, POINTS
                  IF (.NOT. MARK(K)) THEN
                     OLD = OLD + 1
c                     dx = 0.0d0
                     IF (1 .LT. K) THEN
c                        dx = x(k)-x(k-1)
                        IF (VARY1(OLD - 1) .NE. 0.0) THEN
                           WRITE (WORD, '(F4.2, I4)') RATIO1(OLD - 1),
     +                        VARY1(OLD - 1)
                        ELSE
                           WRITE (WORD, '(F4.2, I4)') RATIO1(OLD - 1)
                        END IF
                        IF (MARK(K - 1)) THEN
                           WRITE (TEXT, 10007) K - 1, X(K - 1), WORD
                        ELSE
                           WRITE (TEXT, 10008) WORD
                        END IF
                     END IF

                     IF (1 .LT. K .AND. K .LT. POINTS) THEN
                        IF (VARY2(OLD) .NE. 0.0) THEN
                           WRITE (WORD, '(F4.2, I4)') RATIO2(OLD),
     +                        VARY2(OLD)
                        ELSE
                           WRITE (WORD, '(F4.2, I4)') RATIO2(OLD)
                        END IF
                        WRITE (TEXT, 10009) K, X(K), WORD
                     ELSE
                        WRITE (TEXT, 10009) K, X(K)
                     END IF
                  END IF
5030           CONTINUE
            END IF
         END IF

         IF (0 .LT. LEVELD .AND. 0 .LT. MORE) THEN
            WRITE (TEXT, 10010) ID
            CALL TWCOPY (COMPS * POINTS, U, BUFFER)
            SIGNAL = 'SHOW'
C           GO TO 5040 WHEN ROUTE = 2
            ROUTE = 2
            GO TO 99999
         END IF
      END IF
5040  CONTINUE
      SIGNAL = ' '

C///  SET THE COMPLETION STATUS FLAGS.

      NEWX = 0 .LT. MORE
      SUCCES = 0 .EQ. MOST

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 FORMAT
     +  (/1X, A9, 'SELECT A GRID.')

10002 FORMAT
     +  (/1X, A9, 'SUCCESS.  THE GRID IS ADEQUATE BECAUSE ALL ACTIVE'
     +  /10X, 'COMPONENTS ARE INSIGNIFICANT.')

10003 FORMAT
C              123456789-   1234567   1234567
     + (/15X, '             RATIO 1   RATIO 2'
     +  /15X, '             -------   -------'
     +  /15X, '    ACTUAL', 2F10.3
     +  /15X, '   DESIRED', 2F10.3)

10004 FORMAT
     +  (/1X, A9, 'SUCCESS.  THE GRID IS ADEQUATE.')

10005 FORMAT
     +  (/1X, A9, 'FAILURE.  MORE POINTS ARE NEEDED BUT NONE CAN BE'
     +  /10X, 'ADDED.')

10006 FORMAT
     + (/10X, 'THE NEW GRID (* MARKS NEW POINTS):'
C              123456   123456789-123456   12345678   12345678
     + //10X, '                             LARGEST RATIOS AND'
     +  /10X, ' INDEX         GRID POINT      NUMBER TOO LARGE'
     +  /10X, '------   ----------------   -------------------'
     +  /10X, '                             RATIO 1    RATIO 2')

10007 FORMAT
     +  (10X, I6, '*  ', 1P, E16.9, 0P, 3X, A8)

10008 FORMAT
     +  (38X, A8)

10009 FORMAT
     +  (10X, I6, '   ', 1P, E16.9, 0P, 14X, A8)

10010 FORMAT
     +  (/1X, A9, 'THE SOLUTION GUESS FOR THE NEW GRID:')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, ROUTE
      IF (.NOT. MESS) GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, COMPS, POINTS
      IF (.NOT. MESS) GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, PADD
      IF (.NOT. MESS) GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, POINTS, PMAX
      IF (.NOT. MESS) GO TO 99999

9005  IF (0 .LT. TEXT) WRITE (TEXT, 99005) ID
      IF (.NOT. MESS) GO TO 99999

9006  IF (0 .LT. TEXT) WRITE (TEXT, 99006) ID, TOLER0
      IF (.NOT. MESS) GO TO 99999

9007  IF (0 .LT. TEXT) WRITE (TEXT, 99007) ID, TOLER1, TOLER2
      IF (.NOT. MESS) GO TO 99999

9008  IF (0 .LT. TEXT) WRITE (TEXT, 99008) ID
      IF (.NOT. MESS) GO TO 99999

9009  IF (0 .LT. TEXT) WRITE (TEXT, 99009) ID
      IF (.NOT. MESS) GO TO 99999

99001 FORMAT
     +  (/1X, A9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, I10, '  ROUTE')

99002 FORMAT
     +  (/1X, A9, 'ERROR.  THERE MUST BE AT LEAST ONE COMPONENT AND AT'
     +  /10X, 'LEAST TWO POINTS.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS')

99003 FORMAT
     +  (/1X, A9, 'ERROR.  THE LIMIT ON POINTS ADDED TO A GRID MUST BE'
     +  /10X, 'ZERO OR POSITIVE.'
     + //10X, I10, '  PADD, LIMIT ON ADDED POINTS')

99004 FORMAT
     +  (/1X, A9, 'ERROR.  POINTS IS OUT OF RANGE.'
     + //10X, I10, '  POINTS'
     +  /10X, I10, '  PMAX, LIMIT ON POINTS')

99005 FORMAT
     +  (/1X, A9, 'ERROR.  THERE ARE NO ACTIVE COMPONENTS.')

99006 FORMAT
     +  (/1X, A9, 'ERROR.  THE BOUNDS ON MAGNITUDE AND RELATIVE CHANGE'
     +  /10X, 'OF MAGNITUDE FOR INSIGNIFICANT COMPONENTS MUST BE'
     +  /10X, 'POSITIVE.'
     + //10X, 1P, E10.2, '  TOLER0, SIGNIFICANCE LEVEL')

99007 FORMAT
     +  (/1X, A9, 'ERROR.  THE BOUNDS ON RELATIVE CHANGES IN MAGNITUDE'
     +  /10X, 'AND ANGLE MUST LIE BETWEEN 0 AND 1.'
     + //10X, 1P, E10.2, '  TOLER1'
     +  /10X, 1P, E10.2, '  TOLER2')

99008 FORMAT
     +  (/1X, A9, 'ERROR.  THE GRID IS NOT ORDERED.')

99009 FORMAT
     +  (/1X, A9, 'ERROR.  SOME INTERVALS IN THE GRID ARE TOO SHORT.'
     +  /10X, 'THE NEW GRID WOULD NOT BE ORDERED.')

      STOP
99999 CONTINUE
      RETURN
      END
