C  CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      SUBROUTINE SPLIFT (X,Y,YP,YPP,N,W,IERR,ISX,A1,B1,AN,BN)           SPLIFT 2
C                                                                       ADDRESS2
C     SANDIA MATHEMATICAL PROGRAM LIBRARY                               ADDRESS3
C     APPLIED MATHEMATICS DIVISION 2646                                 CHANGE31
C     SANDIA LABORATORIES                                               ADDRESS5
C     ALBUQUERQUE, NEW MEXICO  87185                                    JUN02781
C     CONTROL DATA 6600/7600  VERSION 8.1  AUGUST 1980                  AUG14801
C                   *************************                           AUG22792
C                   *       ISSUED BY       *                           AUG22793
C                   *  SANDIA LABORATORIES, *                           AUG22794
C                   *   A PRIME CONTRACTOR  *                           AUG22795
C                   ********     TO THE     *                           AUG22796
C                          *  UNITED STATES *                           AUG22797
C                          *   DEPARTMENT   *                           AUG22798
C                          *       OF       *                           AUG22799
C                          *     ENERGY     *                           AUG22710
C      *********************  ---NOTICE---  *********************       AUG22711
C      *THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED*       AUG22712
C      *  BY THE UNITED STATES GOVERNMENT.  NEITHER THE UNITED  *       AUG22713
C      *   STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,   *       AUG22714
C      *               NOR ANY OF THEIR EMPLOYEES,              *       AUG22715
C      * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR *       AUG22716
C      * EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR  *       AUG22717
C      * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE  *       AUG22718
C      *          **********    ACCURACY,   **********          *       AUG22719
C      *          *        *  COMPLETENESS  *        *          *       AUG22720
C      *          *        *  OR USEFULNESS *        *          *       AUG22721
C      *          *        *     OF ANY     *        *          *       AUG22722
C      *          *        *  INFORMATION,  *        *          *       AUG22723
C      *          *        *   APPARATUS,   *        *          *       AUG22724
C      *       ****        *     PRODUCT    *        ****       *       AUG22725
C      *       *           *   OR PROCESS   *           *       *       AUG22726
C      *       *           *   DISCLOSED,   *           *       *       AUG22727
C      *       *           *  OR REPRESENTS *           *       *       AUG22728
C      *       *          **    THAT ITS    **          *       *       AUG22729
C      *       *          **  USE WOULD NOT **          *       *       AUG22730
C      *********          **    INFRINGE    **          *********       AUG22731
C                         **    PRIVATELY   **                          AUG22732
C                         **      OWNED     **                          AUG22733
C                         **     RIGHTS.    **                          AUG22734
C                         **                **                          AUG22735
C                         **                **                          AUG22736
C                         **                **                          AUG22737
C                         ********************                          AUG22738
C                                                                       ADDRESS9
C     WRITTEN BY RONDALL E. JONES                                       SPLIFT 4
C                                                                       SPLIFT 5
C     ABSTRACT                                                          SPLIFT 6
C         SPLIFT FITS AN INTERPOLATING CUBIC SPLINE TO THE N DATA POINTSSPLIFT 7
C         GIVEN IN X AND Y AND RETURNS THE FIRST AND SECOND DERIVATIVES SPLIFT 8
C         IN YP AND YPP.  THE RESULTING SPLINE (DEFINED BY X, Y, AND    SPLIFT 9
C         YPP) AND ITS FIRST AND SECOND DERIVATIVES MAY THEN BE         SPLIFT10
C         EVALUATED USING SPLINT.  THE SPLINE MAY BE INTEGRATED USING   SPLIFT11
C         SPLIQ.  FOR A SMOOTHING SPLINE FIT SEE SUBROUTINE SMOO.       SPLIFT12
C                                                                       SPLIFT13
C     DESCRIPTION OF ARGUMENTS                                          SPLIFT14
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST,SPLIFT15
C         E.G.   X(N), Y(N), YP(N), YPP(N), W(3N)                       SPLIFT16
C                                                                       SPLIFT17
C       --INPUT--                                                       SPLIFT18
C                                                                       SPLIFT19
C         X    - ARRAY OF ABSCISSAS OF DATA (IN INCREASING ORDER)       SPLIFT20
C         Y    - ARRAY OF ORDINATES OF DATA                             SPLIFT21
C         N    - THE NUMBER OF DATA POINTS.  THE ARRAYS X, Y, YP, AND   SPLIFT22
C                YPP MUST BE DIMENSIONED AT LEAST N.  (N .GE. 4)        SPLIFT23
C         ISX  - MUST BE ZERO ON THE INITIAL CALL TO SPLIFT.            SPLIFT24
C                IF A SPLINE IS TO BE FITTED TO A SECOND SET OF DATA    SPLIFT25
C                THAT HAS THE SAME SET OF ABSCISSAS AS A PREVIOUS SET,  SPLIFT26
C                AND IF THE CONTENTS OF W HAVE NOT BEEN CHANGED SINCE   SPLIFT27
C                THAT PREVIOUS FIT WAS COMPUTED, THEN ISX MAY BE        SPLIFT28
C                SET TO ONE FOR FASTER EXECUTION.                       SPLIFT29
C         A1,B1,AN,BN - SPECIFY THE END CONDITIONS FOR THE SPLINE WHICH SPLIFT30
C                ARE EXPRESSED AS CONSTRAINTS ON THE SECOND DERIVATIVE  SPLIFT31
C                OF THE SPLINE AT THE END POINTS (SEE YPP).             SPLIFT32
C                THE END CONDITION CONSTRAINTS ARE                      SPLIFT33
C                        YPP(1) = A1*YPP(2) + B1                        SPLIFT34
C                AND                                                    SPLIFT35
C                        YPP(N) = AN*YPP(N-1) + BN                      SPLIFT36
C                WHERE                                                  SPLIFT37
C                        ABS(A1).LT. 1.0  AND  ABS(AN).LT. 1.0.         SPLIFT38
C                                                                       SPLIFT39
C                THE SMOOTHEST SPLINE (I.E., LEAST INTEGRAL OF SQUARE   SPLIFT40
C                OF SECOND DERIVATIVE) IS OBTAINED BY A1=B1=AN=BN=0.    SPLIFT41
C                IN THIS CASE THERE IS AN INFLECTION AT X(1) AND X(N).  SPLIFT42
C                IF THE DATA IS TO BE EXTRAPOLATED (SAY, BY USING SPLINTSPLIFT43
C                TO EVALUATE THE SPLINE OUTSIDE THE RANGE X(1) TO X(N)),SPLIFT44
C                THEN TAKING A1=AN=0.5 AND B1=BN=0 MAY YIELD BETTER     SPLIFT45
C                RESULTS.  IN THIS CASE THERE IS AN INFLECTION          SPLIFT46
C                AT X(1) - (X(2)-X(1)) AND AT X(N) + (X(N)-X(N-1)).     SPLIFT47
C                IN THE MORE GENERAL CASE OF A1=AN=A  AND B1=BN=0,      SPLIFT48
C                THERE IS AN INFLECTION AT X(1) - (X(2)-X(1))*A/(1.0-A) SPLIFT49
C                AND AT X(N) + (X(N)-X(N-1))*A/(1.0-A).                 SPLIFT50
C                                                                       SPLIFT51
C                A SPLINE THAT HAS A GIVEN FIRST DERIVATIVE YP1 AT X(1) SPLIFT52
C                AND YPN AT Y(N) MAY BE DEFINED BY USING THE            SPLIFT53
C                FOLLOWING CONDITIONS.                                  SPLIFT54
C                                                                       SPLIFT55
C                A1=-0.5                                                SPLIFT56
C                                                                       SPLIFT57
C                B1= 3.0*((Y(2)-Y(1))/(X(2)-X(1))-YP1)/(X(2)-X(1))      SPLIFT58
C                                                                       SPLIFT59
C                AN=-0.5                                                SPLIFT60
C                                                                       SPLIFT61
C                BN=-3.0*((Y(N)-Y(N-1))/(X(N)-X(N-1))-YPN)/(X(N)-X(N-1))SPLIFT62
C                                                                       SPLIFT63
C       --OUTPUT--                                                      SPLIFT64
C                                                                       SPLIFT65
C         YP   - ARRAY OF FIRST DERIVATIVES OF SPLINE (AT THE X(I))     SPLIFT66
C         YPP  - ARRAY OF SECOND DERIVATIVES OF SPLINE (AT THE X(I))    SPLIFT67
C         IERR - A STATUS CODE                                          SPLIFT68
C              --NORMAL CODE                                            SPLIFT69
C                 1 MEANS THAT THE REQUESTED SPLINE WAS COMPUTED.       SPLIFT70
C              --ABNORMAL CODES                                         SPLIFT71
C                 2 MEANS THAT N, THE NUMBER OF POINTS, WAS .LT. 4.     SPLIFT72
C                 3 MEANS THE ABSCISSAS WERE NOT STRICTLY INCREASING.   SPLIFT73
C                                                                       SPLIFT74
C       --WORK--                                                        SPLIFT75
C                                                                       SPLIFT76
C         W    - ARRAY OF WORKING STORAGE DIMENSIONED AT LEAST 3N.      SPLIFT77
      DIMENSION X(N),Y(N),YP(N),YPP(N),W(N,3)                           SPLIFT79
C                                                                       SPLIFT80
      IF (N.LT.4) GO TO 200                                             SPLIFT82
      NM1  = N-1                                                        SPLIFT83
      NM2  = N-2                                                        SPLIFT84
      IF (ISX.GT.0) GO TO 40                                            SPLIFT85
      DO 5 I=2,N                                                        SPLIFT86
      IF (X(I)-X(I-1)) 300,300,5                                        SPLIFT87
    5 CONTINUE                                                          SPLIFT88
C                                                                       SPLIFT89
C     DEFINE THE TRIDIAGONAL MATRIX                                     SPLIFT90
C                                                                       SPLIFT91
      W(1,3) = X(2)-X(1)                                                SPLIFT92
      DO 10 I=2,NM1                                                     SPLIFT93
      W(I,2) = W(I-1,3)                                                 SPLIFT94
      W(I,3) = X(I+1)-X(I)                                              SPLIFT95
   10 W(I,1) = 2.0*(W(I,2)+W(I,3))                                      SPLIFT96
      W(1,1) = 4.0                                                      SPLIFT97
      W(1,3) =-4.0*A1                                                   SPLIFT98
      W(N,1) = 4.0                                                      SPLIFT99
      W(N,2) =-4.0*AN                                                   SPLIF100
C                                                                       SPLIF101
C     L U DECOMPOSITION                                                 SPLIF102
C                                                                       SPLIF103
      DO 30 I=2,N                                                       SPLIF104
      W(I-1,3) = W(I-1,3)/W(I-1,1)                                      SPLIF105
   30 W(I,1)   = W(I,1) - W(I,2)*W(I-1,3)                               SPLIF106
C                                                                       SPLIF107
C     DEFINE *CONSTANT* VECTOR                                          SPLIF108
C                                                                       SPLIF109
   40 YPP(1) = 4.0*B1                                                   SPLIF110
      DOLD   = (Y(2)-Y(1))/W(2,2)                                       SPLIF111
      DO 50 I=2,NM2                                                     SPLIF112
      DNEW   = (Y(I+1) - Y(I))/W(I+1,2)                                 SPLIF113
      YPP(I) = 6.0*(DNEW - DOLD)                                        SPLIF114
      YP(I)  = DOLD                                                     SPLIF115
   50 DOLD   = DNEW                                                     SPLIF116
      DNEW   = (Y(N)-Y(N-1))/(X(N)-X(N-1))                              SPLIF117
      YPP(NM1) = 6.0*(DNEW - DOLD)                                      SPLIF118
      YPP(N) = 4.0*BN                                                   SPLIF119
      YP(NM1)= DOLD                                                     SPLIF120
      YP(N)  = DNEW                                                     SPLIF121
C                                                                       SPLIF122
C     FORWARD SUBSTITUTION                                              SPLIF123
C                                                                       SPLIF124
      YPP(1) = YPP(1)/W(1,1)                                            SPLIF125
      DO 60 I=2,N                                                       SPLIF126
   60 YPP(I) = (YPP(I) - W(I,2)*YPP(I-1))/W(I,1)                        SPLIF127
C                                                                       SPLIF128
C     BACKWARD SUBSTITUTION                                             SPLIF129
C                                                                       SPLIF130
      DO 70 J=1,NM1                                                     SPLIF131
      I = N-J                                                           SPLIF132
   70 YPP(I) = YPP(I) - W(I,3)*YPP(I+1)                                 SPLIF133
C                                                                       SPLIF134
C     COMPUTE FIRST DERIVATIVES                                         SPLIF135
C                                                                       SPLIF136
      YP(1)  = (Y(2)-Y(1))/(X(2)-X(1)) - (X(2)-X(1))*(2.0*YPP(1)        SPLIF137
     1         + YPP(2))/6.0                                            SPLIF138
      DO 80 I=2,NM1                                                     SPLIF139
   80 YP(I)  = YP(I) + W(I,2)*(YPP(I-1) + 2.0*YPP(I))/6.0               SPLIF140
      YP(N)  = YP(N) + (X(N)-X(NM1))*(YPP(NM1) + 2.0*YPP(N))/6.0        SPLIF141
C                                                                       SPLIF142
      IERR = 1                                                          SPLIF143
      RETURN                                                            SPLIF144
  200 IERR = 2                                                          SPLIF145
      WRITE (6, *) 'SPLIFT-THERE WERE LESS THAN 4 DATA VALUES.'         SPLIF146
      RETURN                                                            SPLIF147
  300 IERR = 3                                                          SPLIF148
      WRITE (6, *) 'SPLIFT-THE ABSCISSAS WERE NOT STRICTLY INCREASING.' SPLIF149
      RETURN                                                            SPLIF151
      END                                                               SPLIF152
      SUBROUTINE SPLINT (X,Y,YPP,N,XI,YI,YPI,YPPI,NI,KERR)              SPLINT 2
C                                                                       ADDRESS2
C     SANDIA MATHEMATICAL PROGRAM LIBRARY                               ADDRESS3
C     APPLIED MATHEMATICS DIVISION 2646                                 CHANGE31
C     SANDIA LABORATORIES                                               ADDRESS5
C     ALBUQUERQUE, NEW MEXICO  87185                                    JUN02781
C     CONTROL DATA 6600/7600  VERSION 8.1  AUGUST 1980                  AUG14801
C                   *************************                           AUG22792
C                   *       ISSUED BY       *                           AUG22793
C                   *  SANDIA LABORATORIES, *                           AUG22794
C                   *   A PRIME CONTRACTOR  *                           AUG22795
C                   ********     TO THE     *                           AUG22796
C                          *  UNITED STATES *                           AUG22797
C                          *   DEPARTMENT   *                           AUG22798
C                          *       OF       *                           AUG22799
C                          *     ENERGY     *                           AUG22710
C      *********************  ---NOTICE---  *********************       AUG22711
C      *THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED*       AUG22712
C      *  BY THE UNITED STATES GOVERNMENT.  NEITHER THE UNITED  *       AUG22713
C      *   STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,   *       AUG22714
C      *               NOR ANY OF THEIR EMPLOYEES,              *       AUG22715
C      * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR *       AUG22716
C      * EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR  *       AUG22717
C      * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE  *       AUG22718
C      *          **********    ACCURACY,   **********          *       AUG22719
C      *          *        *  COMPLETENESS  *        *          *       AUG22720
C      *          *        *  OR USEFULNESS *        *          *       AUG22721
C      *          *        *     OF ANY     *        *          *       AUG22722
C      *          *        *  INFORMATION,  *        *          *       AUG22723
C      *          *        *   APPARATUS,   *        *          *       AUG22724
C      *       ****        *     PRODUCT    *        ****       *       AUG22725
C      *       *           *   OR PROCESS   *           *       *       AUG22726
C      *       *           *   DISCLOSED,   *           *       *       AUG22727
C      *       *           *  OR REPRESENTS *           *       *       AUG22728
C      *       *          **    THAT ITS    **          *       *       AUG22729
C      *       *          **  USE WOULD NOT **          *       *       AUG22730
C      *********          **    INFRINGE    **          *********       AUG22731
C                         **    PRIVATELY   **                          AUG22732
C                         **      OWNED     **                          AUG22733
C                         **     RIGHTS.    **                          AUG22734
C                         **                **                          AUG22735
C                         **                **                          AUG22736
C                         **                **                          AUG22737
C                         ********************                          AUG22738
C                                                                       ADDRESS9
C     WRITTEN BY RONDALL E. JONES                                       SPLINT 4
C                                                                       SPLINT 5
C     ABSTRACT                                                          SPLINT 6
C                                                                       SPLINT 7
C         SPLINT EVALUATES A CUBIC SPLINE AND ITS FIRST AND SECOND      SPLINT 8
C         DERIVATIVES AT THE ABSCISSAS IN XI.  THE SPLINE (WHICH        SPLINT 9
C         IS DEFINED BY X, Y, AND YPP) MAY HAVE BEEN DETERMINED BY      SPLINT10
C         SPLIFT OR SMOO OR ANY OTHER SPLINE FITTING ROUTINE THAT       SPLINT11
C         PROVIDES SECOND DERIVATIVES.                                  SPLINT12
C                                                                       SPLINT13
C     DESCRIPTION OF ARGUMENTS                                          SPLINT14
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST,SPLINT15
C         E.G.  X(N), Y(N), YPP(N), XI(NI), YI(NI), YPI(NI), YPPI(NI)   SPLINT16
C                                                                       SPLINT17
C       --INPUT--                                                       SPLINT18
C                                                                       SPLINT19
C         X   - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE THESPLINT20
C               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.  SPLINT21
C         Y   - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y ISSPLINT22
C               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.                SPLINT23
C         YPP - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.     SPLINT24
C               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO. SPLINT25
C         N   - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.       SPLINT26
C               THE ARRAYS X, Y, AND YPP MUST BE DIMENSIONED AT LEAST N.SPLINT27
C               N MUST BE GREATER THAN OR EQUAL TO 2.                   SPLINT28
C         XI  - THE ABSCISSA OR ARRAY OF ABSCISSAS (IN ARBITRARY ORDER) SPLINT29
C               AT WHICH THE SPLINE IS TO BE EVALUATED.                 SPLINT30
C               EACH XI(K) THAT LIES BETWEEN X(1) AND X(N) IS A CASE OF SPLINT31
C               INTERPOLATION.  EACH XI(K) THAT DOES NOT LIE BETWEEN    SPLINT32
C               X(1) AND X(N) IS A CASE OF EXTRAPOLATION.  BOTH CASES   SPLINT33
C               ARE ALLOWED.  SEE DESCRIPTION OF KERR.                  SPLINT34
C         NI  - THE NUMBER OF ABSCISSAS AT WHICH THE SPLINE IS TO BE    SPLINT35
C               EVALUATED.  IF NI IS GREATER THAN 1, THEN XI, YI, YPI,  SPLINT36
C               AND YPPI MUST BE ARRAYS DIMENSIONED AT LEAST NI.        SPLINT37
C               NI MUST BE GREATER THAN OR EQUAL TO 1.                  SPLINT38
C                                                                       SPLINT39
C       --OUTPUT--                                                      SPLINT40
C                                                                       SPLINT41
C         YI  - ARRAY OF VALUES OF THE SPLINE (ORDINATES) AT XI.        SPLINT42
C         YPI - ARRAY OF VALUES OF THE FIRST DERIVATIVE OF SPLINE AT XI.SPLINT43
C         YPPI- ARRAY OF VALUES OF SECOND DERIVATIVES OF SPLINE AT XI.  SPLINT44
C         KERR- A STATUS CODE                                           SPLINT45
C             --NORMAL CODES                                            SPLINT46
C                1 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA SPLINT47
C                  IN XI USING ONLY INTERPOLATION.                      SPLINT48
C                2 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA SPLINT49
C                  IN XI, BUT AT LEAST ONE EXTRAPOLATION WAS PERFORMED. SPLINT50
C             -- ABNORMAL CODE                                          SPLINT51
C                3 MEANS THAT THE REQUESTED NUMBER OF EVALUATIONS, NI,  SPLINT52
C                  WAS NOT POSITIVE.                                    SPLINT53
C                                                                       SPLINT55
      DIMENSION X(N),Y(N),YPP(N),XI(NI),YI(NI),YPI(NI),YPPI(NI)         SPLINT56
C                                                                       SPLINT57
C     CHECK INPUT                                                       SPLINT58
C                                                                       SPLINT59
      IF (NI) 1,1,2                                                     SPLINT61
    1 WRITE (6, *) 'SPLINT-THE REQUESTED NUMBER OF INTERPOLATIONS WAS ',SPLINT62
     1'NOT POSITIVE.'                                                     SPLINT63
      KERR = 3                                                          SPLINT64
      RETURN                                                            SPLINT65
    2 KERR = 1                                                          SPLINT66
      NM1= N-1                                                          SPLINT67
C                                                                       SPLINT68
C     K IS INDEX ON VALUE OF XI BEING WORKED ON.  XX IS THAT VALUE.     SPLINT69
C     I IS CURRENT INDEX INTO X ARRAY.                                  SPLINT70
C                                                                       SPLINT71
      K  = 1                                                            SPLINT72
      XX = XI(1)                                                        SPLINT73
      IF (XX.LT.X(1)) GO TO 90                                          SPLINT74
      IF (XX.GT.X(N)) GO TO 80                                          SPLINT75
      IL = 1                                                            SPLINT76
      IR = N                                                            SPLINT77
C                                                                       SPLINT78
C     BISECTION SEARCH                                                  SPLINT79
C                                                                       SPLINT80
   10 I  = (IL+IR)/2                                                    SPLINT81
      IF (I.EQ.IL) GO TO 100                                            SPLINT82
      IF (XX-X(I)) 20,100,30                                            SPLINT83
   20 IR = I                                                            SPLINT84
      GO TO 10                                                          SPLINT85
   30 IL = I                                                            SPLINT86
      GO TO 10                                                          SPLINT87
C                                                                       SPLINT88
C     LINEAR FORWARD SEARCH                                             SPLINT89
C                                                                       SPLINT90
   50 IF (XX-X(I+1)) 100,100,60                                         SPLINT91
   60 IF (I.GE.NM1) GO TO 80                                            SPLINT92
      I  = I+1                                                          SPLINT93
      GO TO 50                                                          SPLINT94
C                                                                       SPLINT95
C     EXTRAPOLATION                                                     SPLINT96
C                                                                       SPLINT97
   80 KERR = 2                                                          SPLINT98
      I  = NM1                                                          SPLINT99
      GO TO 100                                                         SPLIN100
   90 KERR = 2                                                          SPLIN101
      I  = 1                                                            SPLIN102
C                                                                       SPLIN103
C     INTERPOLATION                                                     SPLIN104
C                                                                       SPLIN105
  100 H  = X(I+1) - X(I)                                                SPLIN106
      H2 = H*H                                                          SPLIN107
      XR = (X(I+1)-XX)/H                                                SPLIN108
      XR2= XR*XR                                                        SPLIN109
      XR3= XR*XR2                                                       SPLIN110
      XL = (XX-X(I))/H                                                  SPLIN111
      XL2= XL*XL                                                        SPLIN112
      XL3= XL*XL2                                                       SPLIN113
      YI(K) = Y(I)*XR + Y(I+1)*XL                                       SPLIN114
     1       -H2*(YPP(I)*(XR-XR3) + YPP(I+1)*(XL-XL3))/6.0              SPLIN115
      YPI(K) = (Y(I+1)-Y(I))/H                                          SPLIN116
     1        +H*(YPP(I)*(1.0-3.0*XR2) - YPP(I+1)*(1.0-3.0*XL2))/6.0    SPLIN117
      YPPI(K) = YPP(I)*XR + YPP(I+1)*XL                                 SPLIN118
C                                                                       SPLIN119
C     NEXT POINT                                                        SPLIN120
C                                                                       SPLIN121
      IF (K.GE.NI) RETURN                                               SPLIN122
      K = K+1                                                           SPLIN123
      XX = XI(K)                                                        SPLIN124
      IF (XX.LT.X(1)) GO TO 90                                          SPLIN125
      IF (XX.GT.X(N)) GO TO 80                                          SPLIN126
      IF (XX-XI(K-1)) 110,100,50                                        SPLIN127
  110 IL = 1                                                            SPLIN128
      IR = I+1                                                          SPLIN129
      GO TO 10                                                          SPLIN130
C                                                                       SPLIN131
      END                                                               SPLIN132
