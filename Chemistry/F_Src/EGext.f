C adapted from EGSLTD1
C=======================================================================
      SUBROUTINE EGSPTC1 (T, WEG, IWEG, PTC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        T         temperature
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       partial thermal conductivity
C
C
C     Two CG iterations are performed on the matrix L_[e] in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[e] in order to evaluate D.
C
C     Note
C     ----
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSPTC1 ( NS, T, WEG(IXTR), WEG(IYTR), 
     &        WEG(IEGWT), WEG(IEGRU), WEG(IEGPA), 
     &        PTC, 
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &        IWEG(IEGLIN),
     &        WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &        WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSPTC1 ( NS, TEMPER, XTR, YTR, WT, 
     &           RU, PATMOS, PTC, 
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMLE ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &               LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 2 * NS
C.....
      ITERMX = 2
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG2 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C.....
c wz      ITERMX = 2
c wz      CALL EGSSI2 ( NS, NG, G, DMI, TEMP, XTR, YTR, WT, WW, Y, 
c wz     &             D, ITERMX, AUX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================

C adapted from EGSLTD2
C=======================================================================
      SUBROUTINE EGSPTC2 (T, WEG, IWEG, PTC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        T         temperature
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       thermal conducitvity
C     
C     We form the Choleski decomposition of matrix L_[e] 
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSPTC2 ( NS, T, WEG(IXTR), WEG(IYTR), 
     &        WEG(IEGWT), WEG(IEGRU), WEG(IEGPA), 
     &        PTC, 
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &        IWEG(IEGLIN),
     &        WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &        WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSPTC2 ( NS, TEMPER, XTR, YTR, WT,
     &           RU, PATMOS, PTC,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), YTR(NS)
C-----------------------------------------------------------------------
      CALL EGSEMLE ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &               LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 2 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         G(II1) = G(II1) + YTR(I) * YTR(I)
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            G(IJ1) = G(IJ1) + YTR(I) * YTR(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGSDEC ( NG, G, TEMP, IER )
      CALL DCOPY  ( NG, BETA, 1, AN, 1 )
      CALL EGSSOL ( NG, G, AN )
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      RETURN
      END
C=======================================================================
