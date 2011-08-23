C  CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:34 $
      SUBROUTINE XSGBCO(A, NDIM1, AA, NDIM2, N, ML, MU, IPVT,
     1                  RCOND, DR , DC, Z, ANORM)
C=======================================================================
C                                                         
C                                                          
C                        FILE =  xsgbco.f          
C                                                          
C    -----------------   VERSION = 1.6        
C    |    SCCS  FILE |          
C    |     SUMMARY   |   CURRENT DATE = 8/1/90 at 13:56:11
C    -----------------                           
C                        DATE OF NEWEST DELTA = 8/1/90 at 13:56:00
C                                                          
C           SCCS file name = s.xsgbco.f                   
C           module type    =    
C           q flag         =                                       
C=======================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C    XDGECO factors the system AX = B, using row and
C column scaling.  A is assumed to be in general LINPACK
C band storage format.
C     On return, A will contain the scaled
C matrix.  AA will contain the factorization of the scaled 
C matrix. RCOND will contain an estimate of the inverse of
C the condition number.  DR and DC will contain the row and
C column scaling vectors respectively.  IPVT will contain
C the pivoting information necessary for reconstructing
C the LU decomposition.  ANORM will contain the matrix norm
C (max |a(j)|), where |a(j)| is the L1 norm of the column
C vectors of A) of the new scaled matrix. 
C   Note Z is a work array of length N.
C
      INTEGER N, IPVT(N), NDIM1, NDIM2
      INTEGER ML, MU
      REAL A(NDIM1,N) , AA(NDIM2,N), DR(N), DC(N)
      REAL RCOND, Z(N), ANORM
C*****xlinpk printing
C      REAL ANORMO
C*****END xlinpk printing

      INTEGER I, J, I1, I2, M, JLNGTH, TEXT
      LOGICAL IPRINT
      SAVE TEXT, IPRINT
      REAL SASUM
      EXTERNAL SASUM, SGBCO

      DATA TEXT/6/
      DATA IPRINT /.FALSE./

C*****xlinpk printing
C      IPRINT = .TRUE.
C*****END xlinpk printing

      M = ML + MU + 1

C CHECK DIMENSIONS
      IF (NDIM1 .LT. (MU+2*ML+1)) THEN
        WRITE(TEXT,*)'XSGBCO: NDIM1 IS TOO SMALL', NDIM1,
     1               (MU+2*ML+1)
        STOP
      END IF
      IF (NDIM2 .LT. (MU+2*ML+1)) THEN
        WRITE(TEXT,*)'XSBGCO: NDIM2 IS TOO SMALL',
     1               NDIM2, (MU+2*ML+1)
        STOP
      END IF

C                Scale A and Store it in AA

      DO 5 I = 1 , N
        DR(I) = 0.0
 5    CONTINUE
      DO 10 J = 1 , N
        I1 = MAX0(1,J-MU)
        I2 = MIN0(N,J+ML)
        DO 7 I = I1, I2
          DR(I) =  ABS(A(I-J+M,J)) + DR(I)
 7      CONTINUE
10    CONTINUE
      DO 11 I = 1 , N
        IF (DR(I) .NE. 0.0) THEN
          DR(I) = 1.0 / DR(I)
        ELSE
C*****xlinpk printing
C          WRITE(TEXT,1020) I
C*****END xlinpk printing
          DR(I) = 1.0
        END IF
 11   CONTINUE
      DO 20 J = 1 , N
        I1 = MAX0(1,J-MU)
        I2 = MIN0(N,J+ML)
        DO 15 I = I1, I2
15        AA(I-J+M,J) = A(I-J+M,J)*DR(I)
20    CONTINUE
C
      ANORM = 0.0
C*****xlinpk printing
C      ANORMO = 0.0
C*****END xlinpk printing
      DO 40 J = 1 , N
      I1 = MAX0(1,J-MU)
      I2 = MIN0(N,J+ML)
      JLNGTH = I2 - I1
C*****xlinpk printing
C      ANORMO = MAX( ANORMO,SASUM(JLNGTH,A(I1-J+M,J),1) )
C*****END xlinpk printing
C*****xlinpk column scaling
      DC(J) = 0.0
      DO 30 I = I1, I2
      DC(J) = ABS(AA(I-J+M,J)) + DC(J)
30    CONTINUE
      IF (DC(J) .NE. 0.0) THEN
        DC(J) = 1.0 / DC(J)
      ELSE
        IF (IPRINT) WRITE(TEXT,1030) J
        DC(J) = 1.0
      END IF
C*****END xlinpk column scaling
      DO 35 I = I1, I2
C*****xlinpk column scaling
      AA(I-J+M,J) = AA(I-J+M,J) * DC(J)
C*****END xlinpk column scaling
35    A(I-J+M,J)  = AA(I-J+M,J)
      ANORM = MAX( ANORM,SASUM(JLNGTH,AA(I1-J+M,J),1) )
40    CONTINUE
C
C*****xlinpk printing
C      WRITE(TEXT,1050) ANORM, ANORMO
C*****END xlinpk printing
C
      CALL SGBCO(AA, NDIM2, N, ML, MU, IPVT, RCOND, Z)
C
      RETURN
1020  FORMAT('XSGBCO: ERROR DR(',I3,') IS ZERO')
1030  FORMAT('XSGBCO: ERROR, DC(',I3,') IS ZERO')
1050  FORMAT('XDGBCO:   New matrix norm (after scaling) = ',G10.5,
     1  ' Old matrix norm = ', G10.5)

      END
      SUBROUTINE XSGBSL(A, NDIM1, AA, NDIM2, N, ML, MU, IPVT,
     1                  DR, DC, B,
     1                  RCOND, X, RELERR, INFO, R, ANORM)
C=======================================================================
C                                                         
C                                                          
C                        FILE =  xsgbsl.f          
C                                                          
C    -----------------   VERSION = 1.3        
C    |    SCCS  FILE |          
C    |     SUMMARY   |   CURRENT DATE = 8/1/90 at 14:04:54
C    -----------------                           
C                        DATE OF NEWEST DELTA = 8/1/90 at 14:04:50
C                                                          
C           SCCS file name = s.xsgbsl.f                   
C           module type    =    
C           q flag         =                                       
C=======================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C     This subroutine solves a linear system AX = B, using 
C iterative improvement and row and column scaling to
C minimize round-off error.  It assumes that the matrix
C A (full matrix) has already been scaled and factored.
C The factored matrix is storred in AA.  This subroutine
C needs double the storage of solvers like DGESL, which
C don't perform iterative improvement of the solution.
C However, the extra work involved may not be that large.
C The iterative improvement involves a few n**2 operations
C at most, while the initial factorization involves N**3
C operations.
C
C INPUT
C-------
C
C A(NDIM1,N) = original matrix N by N
C AA(NDIM2,N) = LU decomposition matrix returned by XDGECO
C IPVT(N)  = vector of pivot information returned by XDGECO
C DR(N) = vector of row scalings returned by XDGECO
C DC(N) = vector of column scalings returned by XDGECO
C N = size of the matrix
C B = right hand side.  
C RCOND = condition number returned by XDGECO
C ANORM = Norm of the scaled matrix.
C
C WORK SPACE
C------------
C X(N) , R(N)
C
C OUTPUT
C--------
C
C B = Contains the solution
C INFO = flag to indicate success or failure
C        0: success
C        -1: failure - iterative improvement didn't converge.
C RELERR = relative error in initial solution.
C
C

C                    DECLARATIONS FOR DUMMY VARIABLES
      INTEGER NDIM1, NDIM2, N
      REAL RCOND , X(N) , RELERR, ANORM
      REAL A(NDIM1,N), AA(NDIM2,N), DR(N), DC(N), 
     1                 B(N), R(N)
      INTEGER INFO, ML, MU
      INTEGER IPVT(N)

C                     DECLARARTIONS FOR LOCAL VARIABLES

      REAL XNORM, RNORM, ERREST, RRNORM
      DOUBLE PRECISION SUM
      LOGICAL ITIMP
      INTEGER I, J, ITER, TEXT, M, J1, J2
      REAL T
      SAVE TEXT, ITIMP
C Externals
      REAL SASUM
      EXTERNAL SASUM, SGBSL

      DATA TEXT/6/
      DATA ITIMP /.FALSE./

C*****xlinpk iterative improvement
      ITIMP = .TRUE.
C*****END xlinpk iterative improvement

      INFO = 0
      M = 1 + MU + ML

      IF (ITIMP) THEN

      DO 10 I = 1 , N
      B(I) = B(I) * DR(I)
10    X(I) = B(I)

C Solve AA X = B using iterative improvement 
      CALL SGBSL(AA, NDIM2, N, ML, MU, IPVT, X, 0)
      XNORM = SASUM(N, X, 1)
      RELERR = 0.0
      IF (XNORM .EQ. 0.0) GO TO 60

      DO 50 ITER = 1, 10
C Use DOUBLE PRECISION to form the right hand side
          DO 30  I = 1 , N
            SUM = DBLE(-B(I))
            J1 = MAX0(1,I-ML)
            J2 = MIN0(N,I+MU)
            DO 25 J = J1, J2
25            SUM = SUM + DBLE(A(I-J+M,J)) * DBLE(X(J))
            R(I) = SNGL(SUM)
30        CONTINUE

        RRNORM = SASUM(N, R, 1)
        IF (ITER .EQ. 1) THEN
          ERREST = RRNORM / (ANORM * XNORM * RCOND)
          IF (ERREST .GE. 1.0E-6) THEN
C*****xlinpk printing
C            WRITE(TEXT,1020) ERREST, 1.0E0/RCOND, RRNORM, XNORM
C*****END xlinpk printing
          ELSE
C*****xlinpk printing
C           WRITE(TEXT,*)'XSGBSL: ERREST = ',ERREST
C*****END xlinpk printing
            GO TO 60
          END IF
        END IF
        CALL SGBSL(AA, NDIM2, N, ML, MU, IPVT, R, 0)
        DO 40 I = 1, N
          X(I) = X(I) - R(I)
40      CONTINUE
        RNORM = SASUM(N, R, 1)
        IF (ITER .EQ. 1) THEN
          RELERR = RNORM/XNORM
        END IF
C*****xlinpk printing
C        WRITE(TEXT,1040) ITER, RNORM, RRNORM
C*****END xlinpk printing
C       CHECK FOR CONVERGENCE 
C         (assume converged if good to single precision)
        T = XNORM + RNORM
        IF (T .EQ. XNORM) GO TO 60
50    CONTINUE
      INFO = -1
60    CONTINUE
      DO 200 I = 1 , N
        B(I) = X(I)
C*****xlinpk column scaling
     1              * DC(I)
C*****END xlinpk column scaling
200   CONTINUE
      RETURN

      ELSE

      DO 210 I = 1 , N
        B(I) = B(I) * DR(I)
210   CONTINUE
      CALL SGBSL(AA, NDIM2, N, ML, MU, IPVT, B, 0)
      RELERR = 0.0
C*****xlinpk column scaling
      DO 220 I = 1, N
220     B(I) = B(I) * DC(I)
C*****END xlinpk column scaling
      RETURN

      END IF

1020  FORMAT('XSGBSL: Error estimate indicates that roundoff'
     1             ,' error may be influencing results. Iterative '
     1             ,'improvement will be attempted.'/
     1 '           ERREST = ',G10.5,' COND_NUM = ',
     1             G10.5,' RNORM = ',G10.5,' XNORM = ', G10.5)
1040  FORMAT('XSGBSL:     ITER = ', I3,' DEL_X_NORM = ',
     1        G10.5,' RNORM = ', G10.5)
      END
      SUBROUTINE XSGECO(A, NDIM1, AA, NDIM2, N, IPVT,
     1                  RCOND, DR , DC, Z, ANORM)
C=======================================================================
C                                                         
C                                                          
C                        FILE =  xsgeco.f          
C                                                          
C    -----------------   VERSION = 1.3        
C    |    SCCS  FILE |          
C    |     SUMMARY   |   CURRENT DATE = 8/1/90 at 13:17:18
C    -----------------                           
C                        DATE OF NEWEST DELTA = 8/1/90 at 13:17:13
C                                                          
C           SCCS file name = s.xsgeco.f                   
C           module type    =    
C           q flag         =                                       
C=======================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C    XSGECO factors the system AX = B, using row and
C column scaling, using single precision.
C  On return, A will contain the scaled
C matrix.  AA will contain the factorization of the scaled 
C matrix. RCOND will contain an estimate of the inverse of
C the condition number.  DR and DC will contain the row and
C column scaling vectors respectively.  IPVT will contain
C the pivoting information necessary for reconstructing
C the LU decomposition.  Note Z is a work array.
C   This version also calculates the matrix norm of the 
C scaled matrix, ANORM, and returns it to the calling routine.
C
      INTEGER N, IPVT(N), NDIM1, NDIM2
      REAL A(NDIM1,N) , AA(NDIM2,N), DR(N), DC(N)
      REAL RCOND, Z(N), ANORM
C*****xlinpk printing
C      REAL ANORMO
C*****END xlinpk printing

      REAL SASUM
      EXTERNAL SGECO, SASUM

      INTEGER I, J, TEXT
      LOGICAL IPRINT
      SAVE TEXT, IPRINT
      DATA TEXT/6/
      DATA IPRINT /.FALSE./

C*****xlinpk printing
C      IPRINT = .TRUE.
C*****END xlinpk printing

C Scale A and Store it in AA

      DO 20 I = 1 , N
      DR(I) = 0.0
      DO 10 J = 1 , N
      DR(I) = ABS(A(I,J)) + DR(I)
10    CONTINUE
      IF (DR(I) .NE. 0.0) THEN
        DR(I) = 1.0 / DR(I)
      ELSE
C*****xlinpk printing
C        PRINT *,'XSGECO: ERROR DR(',I,') IS ZERO'
C        PRINT *,'        SETTING DR TO ONE'
C*****END xlinpk printing
        DR(I) = 1.0
      END IF
      DO 15 J = 1 , N
15    AA(I,J) = A(I,J)*DR(I)
20    CONTINUE
C
      ANORM = 0.0
C*****xlinpk printing
C      ANORMO = 0.0
C*****END xlinpk printing
      DO 40 J = 1 , N
C*****xlinpk printing
C      ANORMO = MAX( ANORMO,SASUM(N,A(1,J),1) )
C*****END xlinpk printing
C*****xlinpk column scaling
      DC(J) = 0.0
      DO 30 I = 1, N
      DC(J) = ABS(AA(I,J)) + DC(J)
30    CONTINUE
      IF (DC(J) .NE. 0.0) THEN
        DC(J) = 1.0 / DC(J)
      ELSE
        IF (IPRINT) THEN
          PRINT *,'XSGECO: ERROR, DC(',J,') IS ZERO'
          PRINT *,'        SETTING DC TO ONE'
        END IF
        DC(J) = 1.0
      END IF
C*****END xlinpk column scaling
      DO 35 I = 1 , N
C*****xlinpk column scaling
      AA(I,J) = AA(I,J) * DC(J)
C*****END xlinpk column scaling
35    A(I,J)  = AA(I,J)
      ANORM  = MAX( ANORM,SASUM(N,AA(1,J),1) )
40    CONTINUE

C*****xlinpk printing
C      WRITE(TEXT,*)'XDGECO:   NEW MATRIX NORM = ',ANORM,
C     1  ' OLD MATRIX NORM = ', ANORMO
C*****END xlinpk printing

      CALL SGECO(AA, NDIM2, N, IPVT, RCOND, Z)

      RETURN
      END
      SUBROUTINE XSGESL(A, NDIM1, AA, NDIM2, N, IPVT, DR, DC, B,
     1                  RCOND, X, RELERR, INFO, R, ANORM)
C=======================================================================
C                                                         
C                                                          
C                        FILE =  xsgesl.f          
C                                                          
C    -----------------   VERSION = 1.2        
C    |    SCCS  FILE |          
C    |     SUMMARY   |   CURRENT DATE = 8/1/90 at 13:28:12
C    -----------------                           
C                        DATE OF NEWEST DELTA = 8/1/90 at 13:28:04
C                                                          
C           SCCS file name = s.xsgesl.f                   
C           module type    =    
C           q flag         =                                       
C=======================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C     This subroutine solves a linear system AX = B, using 
C iterative improvement and row and column scaling to
C minimize round-off error.  It assumes that the matrix
C A (full matrix) has already been scaled and factored.
C The factored matrix is storred in AA.  This subroutine
C needs double the storage of solvers like SGESL, which
C don't perform iterative improvement of the solution.
C However, the extra work involved may not be that large.
C The iterative improvement involves a few n**2 operations
C at most, while the initial factorization involves N**3
C operations.
C
C INPUT
C-------
C
C A(NDIM1,N) = original matrix N by N
C AA(NDIM2,N) = LU decomposition matrix returned by XSGECO
C IPVT(N)  = vector of pivot information returned by XSGECO
C DR(N) = vector of row scalings returned by XSGECO
C DC(N) = vector of column scalings returned by XSGECO
C N = size of the matrix
C B = right hand side.  
C RCOND = condition number returned by XSGECO
C ANORM = Norm of the scaled matrix.
C
C WORK SPACE
C------------
C X(N) , R(N)
C
C OUTPUT
C--------
C
C B = Contains the solution
C INFO = flag to indicate success or failure
C        0: success
C        -1: failure - iterative improvement didn't converge.
C RELERR = relative error in initial solution.
C
C

C                    DECLARATIONS FOR DUMMY VARIABLES

      REAL RCOND , X(N) , RELERR, ANORM
      REAL A(NDIM1,N), AA(NDIM2,N), DR(N), DC(N),
     1                 B(N), R(N)
      INTEGER N, NDIM1, NDIM2, INFO, TEXT
      INTEGER IPVT(N)

C                     DECLARARTIONS FOR LOCAL VARIABLES

      REAL XNORM, RNORM, ERREST, RRNORM
      DOUBLE PRECISION SUM
      LOGICAL ITIMP
      SAVE ITIMP
      INTEGER I, J, ITER
      REAL T

C Externals
      REAL SASUM
      EXTERNAL SASUM, SGESL

      DATA TEXT/6/
      DATA ITIMP /.FALSE./

C*****xlinpk iterative improvement
      ITIMP = .TRUE.
C*****END xlinpk iterative improvement

      INFO = 0

      IF (ITIMP) THEN

      DO 10 I = 1 , N
      B(I) = B(I) * DR(I)
10    X(I) = B(I)


C Solve AA X = B using iterative improvement
      CALL SGESL(AA,NDIM2,N,IPVT,X,0)
      XNORM = SASUM(N,X,1)
      RELERR = 0.0
      IF (XNORM .EQ. 0.0) GO TO 60
      DO 50 ITER = 1, 10
        DO 30 I = 1, N
          SUM = DBLE(-B(I))
          DO 25 J = 1 , N
            SUM = SUM + DBLE(A(I,J)) * DBLE(X(J))
 25       CONTINUE
          R(I) = SNGL(SUM)
 30     CONTINUE
        RRNORM = SASUM(N,R,1)
        IF (ITER .EQ. 1) THEN
          ERREST = RRNORM / (ANORM * XNORM * RCOND)
          IF (ERREST .GE. 1.0E-6) THEN
C*****xlinpk printing
C            WRITE(TEXT,1020) ERREST, 
C     1                     1.0D0/RCOND, RRNORM, XNORM
C*****END xlinpk printing
          ELSE
           GO TO 60
          END IF
        END IF
        CALL SGESL(AA, NDIM2, N, IPVT, R, 0)
        DO 40 I = 1, N
          X(I) = X(I) - R(I)
40      CONTINUE
          RNORM = SASUM(N,R,1)
        IF (ITER .EQ. 1) THEN
          RELERR = RNORM/XNORM
        END IF
C*****xlinpk printing
C        WRITE(TEXT,1040) ITER, RNORM, RRNORM
C*****END xlinpk printing
C       CHECK FOR CONVERGENCE
C         (assume converged if good to single precision)
        T = XNORM + RNORM
        IF (T .EQ. XNORM) GO TO 60
50    CONTINUE
      INFO = -1
60    CONTINUE
      DO 200 I = 1 , N
        B(I) = X(I) 
C*****xlinpk column scaling
     1                    * DC(I)
C*****END xlinpk column scaling
200   CONTINUE
      RETURN

      ELSE

      DO 210 I = 1 , N
        B(I) = B(I) * DR(I)
210   CONTINUE
      CALL SGESL(AA,NDIM2,N,IPVT,X,0)
      RELERR = 0.0
C*****xlinpk column scaling
      DO 220 I = 1, N
220     B(I) = B(I) * DC(I)
C*****END xlinpk column scaling
      RETURN

      END IF

1020  FORMAT('XSGESL: Error estimate indicates that roundoff'
     1             ,' error may be influencing results. Iterative '
     1             ,'improvement will be attempted.'/
     1 '      ERREST = ',G10.5,' COND_NUM = ',G10.5,' RNORM = ',G10.5
     1 ,' XNORM = ', G10.5)
1040  FORMAT('XSGESL:     ITER = ',I3,' DEL_X_NORM = ',G10.5
     1      ,' RNORM = ',G10.5) 

      END
      SUBROUTINE XDGBCO(A, NDIM1, AA, NDIM2, N, ML, MU, IPVT,
     1                  RCOND, DR , DC, Z, ANORM)
C=======================================================================
C                                                         
C                                                          
C                        FILE =  xdgbco.f          
C                                                          
C    -----------------   VERSION = 1.8        
C    |    SCCS  FILE |          
C    |     SUMMARY   |   CURRENT DATE = 8/1/90 at 13:47:07
C    -----------------                           
C                        DATE OF NEWEST DELTA = 8/1/90 at 13:47:00
C                                                          
C           SCCS file name = s.xdgbco.f                   
C           module type    =    
C           q flag         =                                       
C=======================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C    XDGBCO factors the system AX = B, using row and
C column scaling.  A is assumed to be in general LINPACK
C band storage format.
C     On return, A will contain the scaled
C matrix.  AA will contain the factorization of the scaled 
C matrix. RCOND will contain an estimate of the inverse of
C the condition number.  DR and DC will contain the row and
C column scaling vectors respectively.  IPVT will contain
C the pivoting information necessary for reconstructing
C the LU decomposition.  The matrix norm of the scaled 
C matrix is returned in ANORM.
C  Note Z is a work array.
C
      INTEGER N, IPVT(N), NDIM1, NDIM2, TEXT
      INTEGER ML, MU
      DOUBLE PRECISION A(NDIM1,N) , AA(NDIM2,N), DR(N), DC(N)
      DOUBLE PRECISION RCOND, Z(N), ANORM
C*****xlinpk printing
C      DOUBLE PRECISION ANORMO
C*****END xlinpk printing

      INTEGER I, J, I1, I2, M, JLNGTH
      LOGICAL IPRINT
      SAVE IPRINT
      DOUBLE PRECISION DASUM
      EXTERNAL DASUM, DGBCO

      DATA TEXT/6/
      DATA IPRINT /.FALSE./

C*****xlinpk printing
C      IPRINT = .TRUE.
C*****END xlinpk printing

      M = ML + MU + 1

C CHECK DIMENSIONS
      IF (NDIM1 .LT. (MU+2*ML+1)) THEN
        WRITE(TEXT,*)'NDIM1 IS TOO SMALL', NDIM1, (MU+2*ML+1)
        STOP
      END IF
      IF (NDIM2 .LT. (MU+2*ML+1)) THEN
        WRITE(TEXT,*)'NDIM2 IS TOO SMALL', NDIM2, (MU+2*ML+1)
        STOP
      END IF

C            Scale A and Store it in AA

      DO 5 I = 1 , N
        DR(I) = 0.0D0
 5    CONTINUE
      DO 10 J = 1 , N
        I1 = MAX0(1,J-MU)
        I2 = MIN0(N,J+ML)
        DO 7 I = I1, I2
          DR(I) = DABS(A(I-J+M,J)) + DR(I)
 7      CONTINUE
10    CONTINUE
      DO 11 I = 1 , N
        IF (DR(I) .GT. 0.0D0) THEN
          DR(I) = 1.0D0 / DR(I)
        ELSE
C*****xlinpk printing
C         PRINT *,'XDGBCO: ERROR DR(',I,') IS ZERO'
C         PRINT *,'        SETTING DR TO ONE'
C*****END xlinpk printing
          DR(I) = 1.0D0
        END IF
 11   CONTINUE
      DO 20 J = 1 , N
        I1 = MAX0(1,J-MU)
        I2 = MIN0(N,J+ML)
        DO 15 I = I1, I2
15        AA(I-J+M,J) = A(I-J+M,J)*DR(I)
20    CONTINUE
C
C*****xlinpk column scaling
      DO 40 J = 1 , N
      DC(J) = 0.0D0
      I1 = MAX0(1,J-MU)
      I2 = MIN0(N,J+ML)
      DO 30 I = I1, I2
      DC(J) = DABS(AA(I-J+M,J)) + DC(J)
30    CONTINUE
      IF (DC(J) .GT. 0.0D0) THEN
        DC(J) = 1.0D0 / DC(J)
      ELSE
        IF (IPRINT) THEN
          PRINT *,'XDGBCO: ERROR, DC(',J,') IS ZERO'
          PRINT *,'        SETTING DC TO ONE'
        END IF
        DC(J) = 1.0D0
      END IF
      DO 35 I = I1, I2
35    AA(I-J+M,J) = AA(I-J+M,J) * DC(J)
40    CONTINUE
C*****END xlinpk column scaling



      ANORM = 0.0D0
C*****xlinpk printing
C      ANORMO = 0.0D0
C*****END xlinpk printing
      DO 110 J = 1, N
      I1 = MAX0(1,J-MU)
      I2 = MIN0(N,J+ML)
      JLNGTH = I2 - I1
C*****xlinpk printing
C      ANORMO = DMAX1( ANORMO,DASUM(JLNGTH,A(I1-J+M,J),1) )
C*****END xlinpk printing
      ANORM = DMAX1( ANORM,DASUM(JLNGTH,AA(I1-J+M,J),1) )
C      WRITE(TEXT,*)'XDGBCO: I = ',J,' DR = ',DR(J)
C     1            ,' DC = ',DC(J)
      DO 105 I = I1, I2
        A(I-J+M,J) = AA(I-J+M,J)
105   CONTINUE
110   CONTINUE
C*****xlinpk printing
C      WRITE(TEXT,*)'XDGBCO:   NEW MATRIX NORM = ',ANORM,
C     1  ' OLD MATRIX NORM = ', ANORMO
C*****END xlinpk printing

      CALL DGBCO(AA, NDIM2, N, ML, MU, IPVT, RCOND, Z)
      RETURN
      END
      SUBROUTINE XDGBSL(A, NDIM1, AA, NDIM2, N, ML, MU, IPVT,
     1                  DR, DC, B,
     1                  RCOND, X, RELERR, INFO, R, ANORM)
C=======================================================================
C                                                         
C                                                          
C                        FILE =  xdgbsl.f          
C                                                          
C    -----------------   VERSION = 1.5        
C    |    SCCS  FILE |          
C    |     SUMMARY   |   CURRENT DATE = 8/1/90 at 14:06:02
C    -----------------                           
C                        DATE OF NEWEST DELTA = 8/1/90 at 14:05:55
C                                                          
C           SCCS file name = s.xdgbsl.f                   
C           module type    =    
C           q flag         =                                       
C=======================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C     This subroutine solves a linear system AX = B, using 
C iterative improvement and row and column scaling to
C minimize round-off error.  It assumes that the matrix
C A (full matrix) has already been scaled and factored.
C The factored matrix is storred in AA.  This subroutine
C needs double the storage of solvers like DGESL, which
C don't perform iterative improvement of the solution.
C However, the extra work involved may not be that large.
C The iterative improvement involves a few n**2 operations
C at most, while the initial factorization involves N**3
C operations.
C
C INPUT
C-------
C
C A(NDIM1,N) = original matrix N by N
C AA(NDIM2,N) = LU decomposition matrix returned by XDGECO
C IPVT(N)  = vector of pivot information returned by XDGECO
C DR(N) = vector of row scalings returned by XDGECO
C DC(N) = vector of column scalings returned by XDGECO
C N = size of the matrix
C B = right hand side.  
C RCOND = condition number returned by XDGECO
C ANORM = Norm of the scaled matrix.
C
C WORK SPACE
C------------
C X(N) , R(N)
C
C OUTPUT
C--------
C
C B = Contains the solution
C INFO = flag to indicate success or failure
C        0: success
C        -1: failure - iterative improvement didn't converge.
C RELERR = relative error in initial solution.
C
C

C                    DECLARATIONS FOR DUMMY VARIABLES

      DOUBLE PRECISION RCOND , X(N) , RELERR, ANORM
      DOUBLE PRECISION A(NDIM1,N), AA(NDIM2,N), DR(N), DC(N), 
     1                 B(N), R(N)
      INTEGER N, NDIM1, NDIM2, INFO, ML, MU
      INTEGER IPVT(N)

C                     DECLARARTIONS FOR LOCAL VARIABLES

      DOUBLE PRECISION XNORM, RNORM, ERREST, RRNORM

      INTEGER I, J, ITER, TEXT, M, I1, I2
      LOGICAL ITIMP
      SAVE TEXT, ITIMP
      DOUBLE PRECISION T

C Externals
      DOUBLE PRECISION DASUM
      EXTERNAL DASUM, DGESL

      DATA TEXT/6/
      DATA ITIMP /.FALSE./

C*****xlinpk iterative improvement
      ITIMP = .TRUE.
C*****END xlinpk iterative improvement
      INFO = 0

      IF (ITIMP) THEN

      M = 1 + MU + ML
      DO 10 I = 1 , N
      B(I) = B(I) * DR(I)
10    X(I) = B(I)

C Solve AA X = B using iterative improvement 
      CALL DGBSL(AA, NDIM2, N, ML, MU, IPVT, X, 0)
      XNORM = DASUM(N, X, 1)
      RELERR = 0.0D0
      IF (XNORM .EQ. 0.0D0) GO TO 60

      DO 50 ITER = 1, 3
C (really should use quad precision through statement 30, here)
        DO 23 I = 1, N
          R(I) = -B(I)
23      CONTINUE
        DO 30 J = 1, N
          I1 = MAX0(1,J-MU)
          I2 = MIN0(N,J+ML)
          DO 25 I = I1, I2
            R(I) = R(I) + A(I-J+M,J) * X(J)
25        CONTINUE
30      CONTINUE
        RRNORM = DASUM(N, R, 1)
        IF (ITER .EQ. 1) THEN
          ERREST = RRNORM / (ANORM * XNORM * RCOND)
          IF (ERREST .GE. 1.0D-3) THEN
C*****xlinpk printing
C            WRITE(TEXT,1020)
C            WRITE(TEXT,1030) ERREST, 1.0D0/RCOND, RRNORM, XNORM
C*****END xlinpk printing
          ELSE
C           WRITE(TEXT,*)'XDGBSL: ERREST = ',ERREST
            GO TO 60
          END IF
        END IF
        CALL DGBSL(AA, NDIM2, N, ML, MU, IPVT, R, 0)
        DO 40 I = 1, N
          X(I) = X(I) - R(I)
40      CONTINUE
        RNORM = DASUM(N, R, 1)
        IF (ITER .EQ. 1) THEN
          RELERR = RNORM/XNORM
        END IF
C*****xlinpk printing
C        WRITE(TEXT,1040) ITER, RNORM, RRNORM
C*****END xlinpk printing
C       CHECK FOR CONVERGENCE 
C         (assume converged if good to single precision)
        T = XNORM + (RNORM*1.0D-3)
        IF (T .EQ. XNORM) GO TO 60
50    CONTINUE
      INFO = -1
60    CONTINUE
      DO 200 I = 1 , N
        B(I) = X(I)
C*****xlinpk column scaling
     1               * DC(I)
C*****END xlinpk column scaling
200   CONTINUE
      RETURN

      ELSE

      DO 210 I = 1 , N
        B(I) = B(I) * DR(I)
210   CONTINUE
      CALL DGBSL(AA, NDIM2, N, ML, MU, IPVT, B, 0)
      RELERR = 0.0D0
C*****xlinpk column scaling
      DO 220 I = 1, N
220     B(I) = B(I) * DC(I)
C*****END xlinpk column scaling
      RETURN

      END IF
1020  FORMAT('XDGBSL: Error estimate indicates that roundoff'
     1 ,' error may be influencing results. Iterative '
     1             ,'improvement will be attempted.')
1030  FORMAT('XDGBSL:      ERREST = ',G10.5,
     1              ' COND_NUM = ',G10.5,' RNORM = ',G10.5,
     1              ' XNORM = ', G10.5)
1040  FORMAT(20X,'XDGBSL:   ITER = ', I2,' DEL_X_NORM = ',E12.5,
     1       ', RNORM = ',E12.5)
      END
      SUBROUTINE XDGECO(A, NDIM1, AA, NDIM2, N, IPVT,
     1                  RCOND, DR , DC, Z, ANORM)
C=======================================================================
C                                                         
C                                                          
C                        FILE =  xdgeco.f          
C                                                          
C    -----------------   VERSION = 1.4        
C    |    SCCS  FILE |          
C    |     SUMMARY   |   CURRENT DATE = 8/1/90 at 13:00:27
C    -----------------                           
C                        DATE OF NEWEST DELTA = 8/1/90 at 13:00:21
C                                                          
C           SCCS file name = s.xdgeco.f                   
C           module type    =    
C           q flag         =                                       
C=======================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C    XDGECO factors the system AX = B, using row and
C column scaling.  On return, A will contain the scaled
C matrix.  AA will contain the factorization of the scaled 
C matrix. RCOND will contain an estimate of the inverse of
C the condition number.  DR and DC will contain the row and
C column scaling vectors respectively.  IPVT will contain
C the pivoting information necessary for reconstructing
C the LU decomposition.  ANORM will contain the matrix
C norm of the scaled matrix.
C
C  Note Z is a work array.
C
      INTEGER N, IPVT(N), NDIM1, NDIM2, TEXT
      DOUBLE PRECISION A(NDIM1,N) , AA(NDIM2,N), DR(N), DC(N)
      DOUBLE PRECISION RCOND, Z(N), ANORM
C*****xlinpk printing
C      DOUBLE PRECISION ANORMO
C*****END xlinpk printing

C         EXTERNALS
      DOUBLE PRECISION DASUM
      EXTERNAL DASUM, DGECO

      INTEGER I, J
      LOGICAL IPRINT
      SAVE IPRINT, TEXT
      DATA TEXT/6/
      DATA IPRINT /.FALSE./

C*****xlinpk printing
C      IPRINT = .TRUE.
C*****END xlinpk printing
C Scale A and Store it in AA

      DO 20 I = 1 , N
      DR(I) = 0.0D0
      DO 10 J = 1 , N
      DR(I) = DABS(A(I,J)) + DR(I)
10    CONTINUE
      IF (DR(I) .GT. 0.0D0) THEN
        DR(I) = 1.0D0 / DR(I)
      ELSE
C*****xlinpk printing
C        WRITE(TEXT,*)'XDGECO: ERROR DR(',I,') IS ZERO'
C        WRITE(TEXT,*)'        SETTING DR TO ONE'
C*****END xlinpk printing
        DR(I) = 1.0D0
      END IF
      DO 15 J = 1 , N
15    AA(I,J) = A(I,J)*DR(I)
20    CONTINUE
C
C*****xlinpk printing
C      ANORMO = 0.0D0
C*****END xlinpk printing
      ANORM  = 0.0D0

      DO 40 J = 1 , N
C*****xlinpk printing
C      ANORMO = DMAX1( ANORMO, DASUM(N,A(1,J),1) ) 
C*****END xlinpk printing
C*****xlinpk column scaling
      DC(J) = 0.0D0
      DO 30 I = 1, N
      DC(J) = DABS(AA(I,J)) + DC(J)
30    CONTINUE
      IF (DC(J) .NE. 0.0D0) THEN
        DC(J) = 1.0D0 / DC(J)
      ELSE
        IF (IPRINT) THEN
          WRITE(TEXT,*)'XDGECO: ERROR, DC(',J,') IS ZERO'
          WRITE(TEXT,*)'        SETTING DC EQUAL TO ONE'
        END IF
        DC(J) = 1.0D0
      END IF
C*****END xlinpk column scaling
      DO 35 I = 1 , N
C*****xlinpk column scaling
      AA(I,J) = AA(I,J) * DC(J)
C*****END xlinpk column scaling
35    A(I,J)  = AA(I,J)
      ANORM  =  DMAX1( ANORM , DASUM(N,AA(1,J),1) )
40    CONTINUE

C*****xlinpk printing
C      WRITE(TEXT,*)'XDGECO:   NEW MATRIX NORM = ',ANORM,
C     1  ' OLD MATRIX NORM = ', ANORMO
C*****END xlinpk printing

      CALL DGECO(AA, NDIM2, N, IPVT, RCOND, Z)

      RETURN
      END
      SUBROUTINE XDGESL(A, NDIM1, AA, NDIM2, N, IPVT, DR, DC, B,
     1                  RCOND, X, RELERR, INFO, R, ANORM)
C=======================================================================
C                                                         
C                                                          
C                        FILE =  xdgesl.f          
C                                                          
C    -----------------   VERSION = 1.5        
C    |    SCCS  FILE |          
C    |     SUMMARY   |   CURRENT DATE = 8/1/90 at 13:09:16
C    -----------------                           
C                        DATE OF NEWEST DELTA = 8/1/90 at 13:09:08
C                                                          
C           SCCS file name = s.xdgesl.f                   
C           module type    =    
C           q flag         =                                       
C=======================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C     This subroutine solves a linear system AX = B, using 
C iterative improvement and row and column scaling to
C minimize round-off error.  It assumes that the matrix
C A (full matrix) has already been scaled and factored.
C The factored matrix is storred in AA.  This subroutine
C needs double the storage of solvers like DGESL, which
C don't perform iterative improvement of the solution.
C However, the extra work involved may not be that large.
C The iterative improvement involves a few n**2 operations
C at most, while the initial factorization involves N**3
C operations.
C
C INPUT
C-------
C
C A(NDIM1,N) = original matrix N by N
C AA(NDIM2,N) = LU decomposition matrix returned by XDGECO
C IPVT(N)  = vector of pivot information returned by XDGECO
C DR(N) = vector of row scalings returned by XDGECO
C DC(N) = vector of column scalings returned by XDGECO
C N = size of the matrix
C B = right hand side.  
C RCOND = condition number returned by XDGECO
C ANORM = Norm of the scaled matrix.
C
C WORK SPACE
C------------
C X(N) , R(N)
C
C OUTPUT
C--------
C
C B = Contains the solution
C INFO = flag to indicate success or failure
C        0: success
C        -1: failure - iterative improvement didn't converge.
C RELERR = relative error in initial solution.
C

C                    DECLARATIONS FOR DUMMY VARIABLES

      DOUBLE PRECISION RCOND , X(N) , RELERR, ANORM
      DOUBLE PRECISION A(NDIM1,N), AA(NDIM2,N), DR(N), DC(N), 
     1                 B(N), R(N)
      INTEGER N, NDIM1, NDIM2, INFO
      INTEGER IPVT(N)

C                     DECLARARTIONS FOR LOCAL VARIABLES

      DOUBLE PRECISION XNORM, RNORM, ERREST, RRNORM

      INTEGER I, J, ITER, TEXT
      LOGICAL ITIMP
      REAL T
      SAVE TEXT, ITIMP

C Externals
      DOUBLE PRECISION DASUM
      EXTERNAL DASUM, DGESL

      DATA TEXT/6/
      DATA ITIMP /.FALSE./

C*****xlinpk iterative improvement
      ITIMP = .TRUE.
C*****END xlinpk iterative improvement

      INFO = 0

      IF (ITIMP) THEN

      DO 10 I = 1 , N
      B(I) = B(I) * DR(I)
10    X(I) = B(I)

C Solve AA X = B using iterative improvement 
      CALL DGESL(AA,NDIM2,N,IPVT,X,0)
      XNORM = DASUM(N,X,1)
      RELERR = 0.0D0
      IF (XNORM .EQ. 0.0D0) GO TO 60
      DO 50 ITER = 1, 3
C (really should use quad precision through statement 30, here)
        DO 30 I = 1, N
          R(I) = -B(I)
          DO 25 J = 1, N
            R(I) = R(I) + A(I,J) * X(J)
25        CONTINUE
30      CONTINUE
        RRNORM = DASUM(N, R, 1)
        IF (ITER .EQ. 1) THEN
          ERREST = RRNORM / (ANORM * XNORM * RCOND)
          IF (ERREST .GE. 1.0D-3) THEN
C*****xlinpk printing
C            WRITE(TEXT,1020)
C            WRITE(TEXT,1030) ERREST, 1.0D0/RCOND, XNORM, RRNORM
C*****END xlinpk printing
          ELSE
            GO TO 60
          END IF
        END IF
        CALL DGESL(AA, NDIM2, N, IPVT, R, 0)
        DO 40 I = 1, N
          X(I) = X(I) - R(I)
40      CONTINUE
        RNORM = DASUM(N, R, 1)
        IF (ITER .EQ. 1) THEN
          RELERR = RNORM/XNORM
        END IF
C*****xlinpk printing
C        WRITE(TEXT,1040) ITER, RNORM, RRNORM
C*****END xlinpk printing
C       CHECK FOR CONVERGENCE 
C         (assume converged if good to single precision)
        T = SNGL(XNORM) + SNGL(RNORM)
        IF (T .EQ. SNGL(XNORM)) GO TO 60
50    CONTINUE
      INFO = -1
60    CONTINUE
      DO 200 I = 1 , N
        B(I) = X(I) 
C*****xlinpk column scaling
     1         * DC(I)
C*****END xlinpk column scaling
200   CONTINUE
      RETURN

      ELSE 

      DO 210 I = 1 , N
        B(I) = B(I) * DR(I)
210   CONTINUE
      CALL DGESL(AA, NDIM2, N, IPVT, B, 0)
      RELERR = 0.0D0
C*****xlinpk column scaling
      DO 220 I = 1, N
220     B(I) = B(I) * DC(I)
C*****END xlinpk column scaling
      RETURN

      END IF

1020  FORMAT('XDGESL: Error estimate indicates that roundoff'
     1             ,' error may be influencing results. Iterative '
     1             ,'improvement will be attempted.')
1030  FORMAT('XDGESL:      ERREST = ',G10.5,' COND_NUM = ',
     1        G10.5,' XNORM = ',G10.5,' RNORM = ', G10.5)
1040  FORMAT('XDGESL:     ITER = ',I3,' DEL_X_NORM = ',G10.5,
     1       ' Rnorm = ',G10.5)
      END
