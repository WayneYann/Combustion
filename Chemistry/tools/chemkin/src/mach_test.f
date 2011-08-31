C     CVS $Revision: 1.1.1.1 $ reposited $Date: 2006/05/26 19:09:33 $

C///////////////////////////////////////////////////////////////////////
C
C     mach_test.f
C
C     test the routines in mach.f
C
C
C     97/10/29  Joseph Grcar
C
C        1. Works on hp, sgi and sun.
C
C///////////////////////////////////////////////////////////////////////

      PROGRAM MAIN

      DOUBLE PRECISION D(5), D1MACH
      INTEGER I(16), I1MACH, J
      REAL R(5), R1MACH

      DO 0100 J = 1, 5
         D(J) = D1MACH(J)
         R(J) = R1MACH(J)
0100  CONTINUE

      DO 0200 J = 1, 16
         I(J) = I1MACH(J)
0200  CONTINUE

      WRITE (6, 10001)
     +   (J, D(J), J = 1, 5), 
     +   (J, R(J), J = 1, 5), 
     +   (J, I(J), J = 1, 16) 

10001 FORMAT
     +  (/1X, 'D1MACH'
     +   /5(/1X, I2, ') ', 1P, D20.10)
     +   //1X, 'R1MACH'
     +   /5(/1X, I2, ') ', 1P, E20.10)
     +   //1X, 'I1MACH'
     +   /16(/1X, I2, ') ', I20))

      END
