C----------------------------------------------------------------------C
C     	
C     26-step reduced mechanism for DME/air
C     (Non-stiff version)
C
C     Developed by:
C
C       Tianfeng Lu & Zhaoyu Luo
C       Department of Mechanical Engineering
C       University of Connecticut
C       191 Audotorium Road U-3139
C       Storrs, CT 06269
C
C       Email: tlu@engr.uconn.edu
C
C     July 18, 2009
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE GETRATES  (P, T, Y, DIFF, DT, ICKWRK, RCKWRK, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*),WDOT(*),ICKWRK(*),RCKWRK(*),DIFF(*)
      DIMENSION RF(175),RB(175),RKLOW(11),C(30)
      DIMENSION XQ(9)
CESR, I have had problems giving this chemical model time steps less than 1ns
CESR  and it does not like a zero timestep, protection against small timesteps has
CESR  has been moved to the input file stiff.in and chemkin_m.f90.
C
      CALL YTCP(P, T, Y, C)
      CALL RATT(T, RF, RB, RKLOW)
      CALL RATX(T, C, RF, RB, RKLOW)
      CALL QSSA(RF, RB, XQ)
Cwqz      CALL STIF(RF, RB, DIFF, DT, C)
      CALL STIF(RF, RB, DIFF, max(1.d-9,DT), C)
      CALL RDOT(RF, RB, WDOT)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
Cwqz      SUBROUTINE GETRATES  (P, T, Y, DIFF, DT, ICKWRK, RCKWRK, WDOT)
      SUBROUTINE DSCKWYP (P, T, Y, ICKWRK, RCKWRK, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
Cwqz      DIMENSION Y(*),WDOT(*),ICKWRK(*),RCKWRK(*),DIFF(*)
      DIMENSION Y(*),WDOT(*),ICKWRK(*),RCKWRK(*)
      DIMENSION RF(175),RB(175),RKLOW(11),C(30)
      DIMENSION XQ(9)
CESR, I have had problems giving this chemical model time steps less than 1ns
CESR  and it does not like a zero timestep, protection against small timesteps has
CESR  has been moved to the input file stiff.in and chemkin_m.f90.
C
      CALL YTCP(P, T, Y, C)
      CALL RATT(T, RF, RB, RKLOW)
      CALL RATX(T, C, RF, RB, RKLOW)
      CALL QSSA(RF, RB, XQ)
Cwqz      CALL STIF(RF, RB, DIFF, DT, C)
      CALL RDOT(RF, RB, WDOT)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE YTCP (P, T, Y, C)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), C(*)
      DATA SMALL/1D-200/
C
      C(1) = Y(1)/1.00796998D0
      C(2) = Y(2)/2.01593995D0
      C(3) = Y(3)/1.50350603D1
      C(4) = Y(4)/1.59994001D1
      C(5) = Y(5)/1.60430303D1
      C(6) = Y(6)/1.70073701D1
      C(7) = Y(7)/1.80153401D1
      C(8) = Y(8)/2.60382407D1
      C(9) = Y(9)/2.80105505D1
      C(10) = Y(10)/2.80541806D1
      C(11) = Y(11)/2.90621506D1
      C(12) = Y(12)/3.00264904D1
      C(13) = Y(13)/3.00701206D1
      C(14) = Y(14)/3.10344604D1
      C(15) = Y(15)/3.19988003D1
      C(16) = Y(16)/3.30067703D1
      C(17) = Y(17)/3.40147402D1
      C(18) = Y(18)/4.40099506D1
      C(19) = Y(19)/4.40535808D1
      C(20) = Y(20)/4.60258906D1
      C(21) = Y(21)/4.60695207D1
      C(22) = Y(22)/5.90450109D1
      C(23) = Y(23)/6.00529809D1
      C(24) = Y(24)/6.20689209D1
      C(25) = Y(25)/7.50444111D1
      C(26) = Y(26)/7.50444111D1
      C(27) = Y(27)/7.7060351D1
      C(28) = Y(28)/9.20517812D1
      C(29) = Y(29)/1.09059151D2
      C(30) = Y(30)/2.80133991D1
C
      SUM = 0D0
      DO K = 1, 30
         SUM = SUM + C(K)
      ENDDO
      SUM = P/(SUM*T*8.314510D7)
C
      DO K = 1, 30
         C(K) = MAX(C(K),SMALL) * SUM
      ENDDO
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RATT (T, RF, RB, RKLOW)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (RU=8.31451D7, SMALL=1.D-200, PATM=1.01325D6)
      DIMENSION RF(*), RB(*), RKLOW(*)
      DIMENSION SMH(38), EG(38)
C
      ALOGT = LOG(T)
      TI = 1D0/T
      TI2 = TI*TI
C
      CALL RDSMH (T, SMH)
      EG(1) = EXP(SMH(1))
      EG(2) = EXP(SMH(2))
      EG(3) = EXP(SMH(3))
      EG(4) = EXP(SMH(4))
      EG(5) = EXP(SMH(5))
      EG(6) = EXP(SMH(6))
      EG(7) = EXP(SMH(7))
      EG(8) = EXP(SMH(8))
      EG(9) = EXP(SMH(9))
      EG(10) = EXP(SMH(10))
      EG(11) = EXP(SMH(11))
      EG(12) = EXP(SMH(12))
      EG(13) = EXP(SMH(13))
      EG(14) = EXP(SMH(14))
      EG(15) = EXP(SMH(15))
      EG(16) = EXP(SMH(16))
      EG(17) = EXP(SMH(17))
      EG(18) = EXP(SMH(18))
      EG(19) = EXP(SMH(19))
      EG(20) = EXP(SMH(20))
      EG(21) = EXP(SMH(21))
      EG(22) = EXP(SMH(22))
      EG(23) = EXP(SMH(23))
      EG(24) = EXP(SMH(24))
      EG(25) = EXP(SMH(25))
      EG(26) = EXP(SMH(26))
      EG(27) = EXP(SMH(27))
      EG(28) = EXP(SMH(28))
      EG(29) = EXP(SMH(29))
      EG(30) = EXP(SMH(30))
      EG(31) = EXP(SMH(31))
      EG(32) = EXP(SMH(32))
      EG(33) = EXP(SMH(33))
      EG(34) = EXP(SMH(34))
      EG(35) = EXP(SMH(35))
      EG(36) = EXP(SMH(36))
      EG(37) = EXP(SMH(37))
      EG(38) = EXP(SMH(38))
      PFAC1 = PATM / (RU*T)
      PFAC2 = PFAC1*PFAC1
      PFAC3 = PFAC2*PFAC1
C
C
      RF(1) = EXP(3.58048786D1 -4.06D-1*ALOGT -8.35289347D3*TI)
      EQK = EG(6)*EG(8)/EG(1)/EG(20)
      RB(1) = RF(1) / MAX(EQK, SMALL)
      RF(2) = EXP(1.08356516D1 +2.67D0*ALOGT -3.16523284D3*TI)
      EQK = EG(1)*EG(8)/EG(2)/EG(6)
      RB(2) = RF(2) / MAX(EQK, SMALL)
      RF(3) = EXP(1.9190789D1 +1.51D0*ALOGT -1.72603317D3*TI)
      EQK = EG(1)*EG(9)/EG(2)/EG(8)
      RB(3) = RF(3) / MAX(EQK, SMALL)
      RF(4) = EXP(1.49040725D1 +2.02D0*ALOGT -6.74310335D3*TI)
      EQK = EG(8)*EG(8)/EG(6)/EG(9)
      RB(4) = RF(4) / MAX(EQK, SMALL)
      RF(5) = EXP(4.52701605D1 -1.4D0*ALOGT -5.25257558D4*TI)
      EQK = EG(1)*EG(1)/EG(2)*PFAC1
      RB(5) = RF(5) / MAX(EQK, SMALL)
      RF(6) = EXP(3.63576645D1 -5D-1*ALOGT)
      EQK = EG(20)/EG(6)/EG(6)/PFAC1
      RB(6) = RF(6) / MAX(EQK, SMALL)
      RF(7) = EXP(4.29970685D1 -1D0*ALOGT)
      EQK = EG(8)/EG(1)/EG(6)/PFAC1
      RB(7) = RF(7) / MAX(EQK, SMALL)
      RF(8) = EXP(5.19918731D1 -2D0*ALOGT)
      EQK = EG(9)/EG(1)/EG(8)/PFAC1
      RB(8) = RF(8) / MAX(EQK, SMALL)
      RF(9) = EXP(2.80196791D1 +6D-1*ALOGT)
      EQK = EG(21)/EG(1)/EG(20)/PFAC1
      RB(9) = RF(9) / MAX(EQK, SMALL)
      RF(10) = EXP(3.04404238D1 -4.14147317D2*TI)
      EQK = EG(2)*EG(20)/EG(1)/EG(21)
      RB(10) = RF(10) / MAX(EQK, SMALL)
      RF(11) = EXP(3.18907389D1 -1.48448917D2*TI)
      EQK = EG(8)*EG(8)/EG(1)/EG(21)
      RB(11) = RF(11) / MAX(EQK, SMALL)
      RF(12) = 3.25D13
      EQK = EG(8)*EG(20)/EG(6)/EG(21)
      RB(12) = RF(12) / MAX(EQK, SMALL)
      RF(13) = EXP(3.09948627D1 +2.50098684D2*TI)
      EQK = EG(9)*EG(20)/EG(8)/EG(21)
      RB(13) = RF(13) / MAX(EQK, SMALL)
      RF(14) = EXP(3.36712758D1 -6.02954211D3*TI)
      EQK = EG(20)*EG(22)/EG(21)/EG(21)
      RB(14) = RF(14) / MAX(EQK, SMALL)
      RF(15) = EXP(2.55908003D1 +8.19890917D2*TI)
      EQK = EG(20)*EG(22)/EG(21)/EG(21)
      RB(15) = RF(15) / MAX(EQK, SMALL)
      RF(16) = EXP(3.33183354D1 -2.43707832D4*TI)
      EQK = EG(8)*EG(8)/EG(22)*PFAC1
      RB(16) = RF(16) / MAX(EQK, SMALL)
      RF(17) = EXP(3.0813233D1 -1.99777017D3*TI)
      EQK = EG(8)*EG(9)/EG(1)/EG(22)
      RB(17) = RF(17) / MAX(EQK, SMALL)
      RF(18) = EXP(3.15063801D1 -4.00057251D3*TI)
      EQK = EG(2)*EG(21)/EG(1)/EG(22)
      RB(18) = RF(18) / MAX(EQK, SMALL)
      RF(19) = EXP(1.60720517D1 +2D0*ALOGT -1.99777017D3*TI)
      EQK = EG(8)*EG(21)/EG(6)/EG(22)
      RB(19) = RF(19) / MAX(EQK, SMALL)
      RF(20) = 1D12
      EQK = EG(9)*EG(21)/EG(8)/EG(22)
      RB(20) = RF(20) / MAX(EQK, SMALL)
      RF(21) = EXP(3.39940492D1 -4.80924169D3*TI)
      EQK = EG(9)*EG(21)/EG(8)/EG(22)
      RB(21) = RF(21) / MAX(EQK, SMALL)
      RF(22) = EXP(2.36136376D1 -1.19966854D3*TI)
      EQK = EG(23)/EG(6)/EG(12)/PFAC1
      RB(22) = RF(22) / MAX(EQK, SMALL)
      RF(23) = EXP(2.85592404D1 -2.4003435D4*TI)
      EQK = EG(6)*EG(23)/EG(12)/EG(20)
      RB(23) = RF(23) / MAX(EQK, SMALL)
      RF(24) = EXP(3.10355463D1 -1.15739834D4*TI)
      EQK = EG(8)*EG(23)/EG(12)/EG(21)
      RB(24) = RF(24) / MAX(EQK, SMALL)
      RF(25) = EXP(1.23144785D1 +1.89D0*ALOGT +5.83077153D2*TI)
      EQK = EG(1)*EG(23)/EG(8)/EG(12)
      RB(25) = RF(25) / MAX(EQK, SMALL)
      RF(26) = EXP(2.68862648D1 +6.59D-1*ALOGT -7.48484471D3*TI)
      EQK = EG(1)*EG(12)/EG(14)*PFAC1
      RB(26) = RF(26) / MAX(EQK, SMALL)
      RF(27) = EXP(2.96565343D1 -2.06318834D2*TI)
      EQK = EG(12)*EG(21)/EG(14)/EG(20)
      RB(27) = RF(27) / MAX(EQK, SMALL)
      RF(28) = 7.23D13
      EQK = EG(2)*EG(12)/EG(1)/EG(14)
      RB(28) = RF(28) / MAX(EQK, SMALL)
      RF(29) = 3.02D13
      EQK = EG(8)*EG(12)/EG(6)/EG(14)
      RB(29) = RF(29) / MAX(EQK, SMALL)
      RF(30) = 3.02D13
      EQK = EG(9)*EG(12)/EG(8)/EG(14)
      RB(30) = RF(30) / MAX(EQK, SMALL)
      RF(31) = 3D13
      EQK = EG(1)*EG(23)/EG(6)/EG(14)
      RB(31) = RF(31) / MAX(EQK, SMALL)
      RF(32) = 3D13
      EQK = EG(1)*EG(8)*EG(23)/EG(14)/EG(21)*PFAC1
      RB(32) = RF(32) / MAX(EQK, SMALL)
      RF(33) = 3D12
      EQK = EG(2)*EG(12)*EG(12)/EG(14)/EG(14)*PFAC1
      RB(33) = RF(33) / MAX(EQK, SMALL)
      RF(34) = 2.65D13
      EQK = EG(7)*EG(12)/EG(5)/EG(14)
      RB(34) = RF(34) / MAX(EQK, SMALL)
      RF(35) = 3D13
      EQK = EG(12)*EG(16)/EG(14)/EG(14)
      RB(35) = RF(35) / MAX(EQK, SMALL)
      RF(36) = EXP(9.09947411D1 -6.3D0*ALOGT -5.02713451D4*TI)
      EQK = EG(1)*EG(14)/EG(16)*PFAC1
      RB(36) = RF(36) / MAX(EQK, SMALL)
      RF(37) = EXP(1.04747731D2 -8D0*ALOGT -4.90686573D4*TI)
      EQK = EG(2)*EG(12)/EG(16)*PFAC1
      RB(37) = RF(37) / MAX(EQK, SMALL)
      RF(38) = EXP(1.78655549D1 +1.9D0*ALOGT -1.38314133D3*TI)
      EQK = EG(2)*EG(14)/EG(1)/EG(16)
      RB(38) = RF(38) / MAX(EQK, SMALL)
      RF(39) = EXP(3.05269331D1 -1.54990734D3*TI)
      EQK = EG(8)*EG(14)/EG(6)/EG(16)
      RB(39) = RF(39) / MAX(EQK, SMALL)
      RF(40) = EXP(2.19558261D1 +1.18D0*ALOGT +2.2493785D2*TI)
      EQK = EG(9)*EG(14)/EG(8)/EG(16)
      RB(40) = RF(40) / MAX(EQK, SMALL)
      RF(41) = EXP(1.40225247D1 +3D0*ALOGT -2.61672667D4*TI)
      EQK = EG(14)*EG(21)/EG(16)/EG(20)
      RB(41) = RF(41) / MAX(EQK, SMALL)
      RF(42) = EXP(1.06237634D1 +2.5D0*ALOGT -5.13784218D3*TI)
      EQK = EG(14)*EG(22)/EG(16)/EG(21)
      RB(42) = RF(42) / MAX(EQK, SMALL)
      RF(43) = EXP(-1.25246264D1 +5.42D0*ALOGT -5.02210234D2*TI)
      EQK = EG(7)*EG(14)/EG(5)/EG(16)
      RB(43) = RF(43) / MAX(EQK, SMALL)
      RF(44) = 8.43D13
      EQK = EG(1)*EG(16)/EG(5)/EG(6)
      RB(44) = RF(44) / MAX(EQK, SMALL)
      RF(45) = EXP(4.21346663D1 -1.57D0*ALOGT -1.47090232D4*TI)
      EQK = EG(6)*EG(19)/EG(5)/EG(20)
      RB(45) = RF(45) / MAX(EQK, SMALL)
      RF(46) = EXP(2.66475216D1 -7.36709201D3*TI)
      EQK = EG(8)*EG(16)/EG(5)/EG(20)
      RB(46) = RF(46) / MAX(EQK, SMALL)
      RF(47) = EXP(2.39054777D1 +7.6D-1*ALOGT +1.16997875D3*TI)
      EQK = EG(8)*EG(19)/EG(5)/EG(21)
      RB(47) = RF(47) / MAX(EQK, SMALL)
      RF(48) = EXP(3.53616352D1 -6.9D-1*ALOGT -8.79924665D1*TI)
      EQK = EG(17)/EG(5)/EG(5)/PFAC1
      RB(48) = RF(48) / MAX(EQK, SMALL)
      RF(49) = EXP(3.70803784D1 -6.3D-1*ALOGT -1.92731984D2*TI)
      EQK = EG(7)/EG(1)/EG(5)/PFAC1
      RB(49) = RF(49) / MAX(EQK, SMALL)
      RF(50) = EXP(1.78173743D1 +1.97D0*ALOGT -5.64105884D3*TI)
      EQK = EG(2)*EG(5)/EG(1)/EG(7)
      RB(50) = RF(50) / MAX(EQK, SMALL)
      RF(51) = EXP(2.87784236D1 +5D-1*ALOGT -5.17809951D3*TI)
      EQK = EG(5)*EG(8)/EG(6)/EG(7)
      RB(51) = RF(51) / MAX(EQK, SMALL)
      RF(52) = EXP(1.55594794D1 +1.96D0*ALOGT -1.32798879D3*TI)
      EQK = EG(5)*EG(9)/EG(7)/EG(8)
      RB(52) = RF(52) / MAX(EQK, SMALL)
      RF(53) = 3.16D12
      EQK = EG(7)*EG(20)/EG(5)/EG(21)
      RB(53) = RF(53) / MAX(EQK, SMALL)
      RF(54) = EXP(2.59217629D1 -9.34976568D3*TI)
      EQK = EG(5)*EG(22)/EG(7)/EG(21)
      RB(54) = RF(54) / MAX(EQK, SMALL)
      RF(55) = EXP(3.22361913D1 -1.26307384D4*TI)
      EQK = EG(1)*EG(16)/EG(18)*PFAC1
      RB(55) = RF(55) / MAX(EQK, SMALL)
      RF(56) = 6D12
      EQK = EG(2)*EG(16)/EG(1)/EG(18)
      RB(56) = RF(56) / MAX(EQK, SMALL)
      RF(57) = 9.635D13
      EQK = EG(5)*EG(8)/EG(1)/EG(18)
      RB(57) = RF(57) / MAX(EQK, SMALL)
      RF(58) = 4.2D13
      EQK = EG(8)*EG(16)/EG(6)/EG(18)
      RB(58) = RF(58) / MAX(EQK, SMALL)
      RF(59) = 2.4D13
      EQK = EG(9)*EG(16)/EG(8)/EG(18)
      RB(59) = RF(59) / MAX(EQK, SMALL)
      RF(60) = EXP(3.3115818D1 -2.52463802D3*TI)
      EQK = EG(16)*EG(21)/EG(18)/EG(20)
      RB(60) = RF(60) / MAX(EQK, SMALL)
      RF(61) = 3.20322444D-4*RF(7)
      EQK = EG(16)*EG(21)/EG(18)/EG(20)
      RB(61) = RF(61) / MAX(EQK, SMALL)
      RF(62) = 1.2D13
      EQK = EG(16)*EG(22)/EG(18)/EG(21)
      RB(62) = RF(62) / MAX(EQK, SMALL)
      RF(63) = 1.5D13
      EQK = EG(16)*EG(16)/EG(14)/EG(18)
      RB(63) = RF(63) / MAX(EQK, SMALL)
      RF(64) = EXP(4.12602021D1 -1.2D0*ALOGT -7.79985835D3*TI)
      EQK = EG(1)*EG(16)/EG(19)*PFAC1
      RB(64) = RF(64) / MAX(EQK, SMALL)
      RF(65) = 3.2D13
      EQK = EG(5)*EG(8)/EG(1)/EG(19)
      RB(65) = RF(65) / MAX(EQK, SMALL)
      RF(66) = 6D12
      EQK = EG(8)*EG(16)/EG(6)/EG(19)
      RB(66) = RF(66) / MAX(EQK, SMALL)
      RF(67) = 1.8D13
      EQK = EG(9)*EG(16)/EG(8)/EG(19)
      RB(67) = RF(67) / MAX(EQK, SMALL)
      RF(68) = EXP(3.21344907D1 -6.02853568D3*TI)
      EQK = EG(16)*EG(21)/EG(19)/EG(20)
      RB(68) = RF(68) / MAX(EQK, SMALL)
      RF(69) = EXP(2.38143083D1 -8.79622735D2*TI)
      EQK = EG(16)*EG(21)/EG(19)/EG(20)
      RB(69) = RF(69) / MAX(EQK, SMALL)
      RF(70) = 3D11
      EQK = EG(16)*EG(22)/EG(19)/EG(21)
      RB(70) = RF(70) / MAX(EQK, SMALL)
      RF(71) = EXP(3.04036098D1 -5.93795668D3*TI)
      EQK = EG(5)*EG(23)/EG(12)/EG(19)
      RB(71) = RF(71) / MAX(EQK, SMALL)
      RF(72) = EXP(2.9238457D1 +1D-1*ALOGT -5.33409668D3*TI)
      EQK = EG(1)*EG(15)/EG(5)/EG(5)
      RB(72) = RF(72) / MAX(EQK, SMALL)
      RF(73) = EXP(1.47156719D1 +2D0*ALOGT -4.16160184D3*TI)
      EQK = EG(5)*EG(5)/EG(3)/EG(7)
      RB(73) = RF(73) / MAX(EQK, SMALL)
      RF(74) = EXP(3.04036098D1 +2.86833501D2*TI)
      EQK = EG(5)*EG(5)/EG(4)/EG(7)
      RB(74) = RF(74) / MAX(EQK, SMALL)
      RF(75) = EXP(1.78408622D1 +1.6D0*ALOGT -2.72743434D3*TI)
      EQK = EG(3)*EG(9)/EG(5)/EG(8)
      RB(75) = RF(75) / MAX(EQK, SMALL)
      RF(76) = 2.501D13
      EQK = EG(4)*EG(9)/EG(5)/EG(8)
      RB(76) = RF(76) / MAX(EQK, SMALL)
      RF(77) = 4D13
      EQK = EG(1)*EG(13)/EG(3)/EG(5)
      RB(77) = RF(77) / MAX(EQK, SMALL)
      RF(78) = 7.5D-1*RF(74)
      EQK = EG(1)*EG(13)/EG(4)/EG(5)
      RB(78) = RF(78) / MAX(EQK, SMALL)
      RF(79) = 1.6D13
      EQK = EG(4)*EG(9)/EG(1)/EG(19)
      RB(79) = RF(79) / MAX(EQK, SMALL)
      RF(80) = EXP(1.85604427D1 +1.9D0*ALOGT -3.78922151D3*TI)
      EQK = EG(2)*EG(15)/EG(1)/EG(17)
      RB(80) = RF(80) / MAX(EQK, SMALL)
      RF(81) = EXP(1.83130955D1 +1.92D0*ALOGT -2.86330284D3*TI)
      EQK = EG(8)*EG(15)/EG(6)/EG(17)
      RB(81) = RF(81) / MAX(EQK, SMALL)
      RF(82) = EXP(1.50796373D1 +2.12D0*ALOGT -4.37798501D2*TI)
      EQK = EG(9)*EG(15)/EG(8)/EG(17)
      RB(82) = RF(82) / MAX(EQK, SMALL)
      RF(83) = EXP(3.13199006D1 -2.56137284D4*TI)
      EQK = EG(15)*EG(21)/EG(17)/EG(20)
      RB(83) = RF(83) / MAX(EQK, SMALL)
      RF(84) = EXP(2.64068456D1 -7.51805701D3*TI)
      EQK = EG(15)*EG(22)/EG(17)/EG(21)
      RB(84) = RF(84) / MAX(EQK, SMALL)
      RF(85) = EXP(1.56303353D1 +1.74D0*ALOGT -5.25861418D3*TI)
      EQK = EG(7)*EG(15)/EG(5)/EG(17)
      RB(85) = RF(85) / MAX(EQK, SMALL)
      RF(86) = EXP(4.07945264D1 -9.9D-1*ALOGT -7.95082335D2*TI)
      EQK = EG(17)/EG(1)/EG(15)/PFAC1
      RB(86) = RF(86) / MAX(EQK, SMALL)
      RF(87) = 2D12
      EQK = EG(2)*EG(13)/EG(1)/EG(15)
      RB(87) = RF(87) / MAX(EQK, SMALL)
      RF(88) = 1.32D14
      EQK = EG(5)*EG(16)/EG(6)/EG(15)
      RB(88) = RF(88) / MAX(EQK, SMALL)
      RF(89) = 2D10
      EQK = EG(13)*EG(21)/EG(15)/EG(20)
      RB(89) = RF(89) / MAX(EQK, SMALL)
      RF(90) = 1.4D12
      EQK = EG(13)*EG(17)/EG(15)/EG(15)
      RB(90) = RF(90) / MAX(EQK, SMALL)
      RF(91) = 1.2D14
      EQK = EG(12)*EG(17)/EG(14)/EG(15)
      RB(91) = RF(91) / MAX(EQK, SMALL)
      RF(92) = 8.02D13
      EQK = EG(1)*EG(24)/EG(6)/EG(15)
      RB(92) = RF(92) / MAX(EQK, SMALL)
      RF(93) = EXP(2.97104627D1 +4.4D-1*ALOGT -4.46705436D4*TI)
      EQK = EG(2)*EG(10)/EG(13)*PFAC1
      RB(93) = RF(93) / MAX(EQK, SMALL)
      RF(94) = EXP(2.77079822D1 +4.54D-1*ALOGT -9.15854335D2*TI)
      EQK = EG(15)/EG(1)/EG(13)/PFAC1
      RB(94) = RF(94) / MAX(EQK, SMALL)
      RF(95) = EXP(2.94360258D1 +2.7D-1*ALOGT -1.40900667D2*TI)
      EQK = EG(13)/EG(1)/EG(11)/PFAC1
      RB(95) = RF(95) / MAX(EQK, SMALL)
      RF(96) = EXP(1.4096923D1 +2.53D0*ALOGT -6.15937201D3*TI)
      EQK = EG(2)*EG(11)/EG(1)/EG(13)
      RB(96) = RF(96) / MAX(EQK, SMALL)
      RF(97) = EXP(1.44032972D1 +2D0*ALOGT -1.25804167D3*TI)
      EQK = EG(9)*EG(11)/EG(8)/EG(13)
      RB(97) = RF(97) / MAX(EQK, SMALL)
      RF(98) = EXP(1.23327053D1 +2D0*ALOGT -4.62959334D3*TI)
      EQK = EG(7)*EG(11)/EG(5)/EG(13)
      RB(98) = RF(98) / MAX(EQK, SMALL)
      RF(99) = EXP(1.67704208D1 +1.83D0*ALOGT -1.10707667D2*TI)
      EQK = EG(5)*EG(14)/EG(6)/EG(13)
      RB(99) = RF(99) / MAX(EQK, SMALL)
      RF(100) = 5D12
      EQK = EG(9)*EG(10)/EG(8)/EG(11)
      RB(100) = RF(100) / MAX(EQK, SMALL)
      RF(101) = EXP(1.65302053D1 +1.91D0*ALOGT -1.88203034D3*TI)
      EQK = EG(8)*EG(11)/EG(6)/EG(13)
      RB(101) = RF(101) / MAX(EQK, SMALL)
      RF(102) = EXP(3.13722558D1 -2.89852801D4*TI)
      EQK = EG(11)*EG(21)/EG(13)/EG(20)
      RB(102) = RF(102) / MAX(EQK, SMALL)
      RF(103) = 9.64D13
      EQK = EG(2)*EG(10)/EG(1)/EG(11)
      RB(103) = RF(103) / MAX(EQK, SMALL)
      RF(104) = EXP(2.32164713D1 +2.99917134D2*TI)
      EQK = EG(13)*EG(21)/EG(11)/EG(22)
      RB(104) = RF(104) / MAX(EQK, SMALL)
      RF(105) = 3.9D11
      EQK = EG(7)*EG(10)/EG(5)/EG(11)
      RB(105) = RF(105) / MAX(EQK, SMALL)
      RF(106) = 9.6D11
      EQK = EG(10)*EG(13)/EG(11)/EG(11)
      RB(106) = RF(106) / MAX(EQK, SMALL)
      RF(107) = EXP(3.83630605D1 -1.39D0*ALOGT -5.10764918D2*TI)
      EQK = EG(14)*EG(16)/EG(11)/EG(20)
      RB(107) = RF(107) / MAX(EQK, SMALL)
      RF(108) = EXP(1.41059389D1 +1.61D0*ALOGT +1.932352D2*TI)
      EQK = EG(10)*EG(21)/EG(11)/EG(20)
      RB(108) = RF(108) / MAX(EQK, SMALL)
      RF(109) = EXP(2.93537877D1 -1.20772D3*TI)
      EQK = EG(11)/EG(1)/EG(10)/PFAC1
      RB(109) = RF(109) / MAX(EQK, SMALL)
      RF(110) = EXP(1.52216075D1 +2D0*ALOGT -9.56111669D2*TI)
      EQK = EG(3)*EG(12)/EG(6)/EG(10)
      RB(110) = RF(110) / MAX(EQK, SMALL)
      RF(111) = EXP(-7.6354939D0 +4D0*ALOGT +1.00643334D3*TI)
      EQK = EG(5)*EG(12)/EG(8)/EG(10)
      RB(111) = RF(111) / MAX(EQK, SMALL)
      RF(112) = EXP(3.77576522D1 -8D-1*ALOGT)
      EQK = EG(5)/EG(1)/EG(3)/PFAC1
      RB(112) = RF(112) / MAX(EQK, SMALL)
      RF(113) = 8D13
      EQK = EG(1)*EG(14)/EG(3)/EG(6)
      RB(113) = RF(113) / MAX(EQK, SMALL)
      RF(114) = 2D13
      EQK = EG(1)*EG(16)/EG(3)/EG(8)
      RB(114) = RF(114) / MAX(EQK, SMALL)
      RF(115) = EXP(1.31223634D1 +2D0*ALOGT -3.63825651D3*TI)
      EQK = EG(1)*EG(5)/EG(2)/EG(3)
      RB(115) = RF(115) / MAX(EQK, SMALL)
      RF(116) = EXP(3.02112379D1 -7.54825001D2*TI)
      EQK = EG(8)*EG(14)/EG(3)/EG(20)
      RB(116) = RF(116) / MAX(EQK, SMALL)
      RF(117) = 2D13
      EQK = EG(8)*EG(16)/EG(3)/EG(21)
      RB(117) = RF(117) / MAX(EQK, SMALL)
      RF(118) = 3.2D13
      EQK = EG(2)*EG(10)/EG(3)/EG(3)
      RB(118) = RF(118) / MAX(EQK, SMALL)
      RF(119) = EXP(2.98282457D1 -3.01930001D2*TI)
      EQK = EG(3)/EG(4)
      RB(119) = RF(119) / MAX(EQK, SMALL)
      RF(120) = 3D13
      EQK = EG(3)/EG(4)
      RB(120) = RF(120) / MAX(EQK, SMALL)
      RF(121) = 9D12
      EQK = EG(3)/EG(4)
      RB(121) = RF(121) / MAX(EQK, SMALL)
      RF(122) = 7D12
      EQK = EG(3)/EG(4)
      RB(122) = RF(122) / MAX(EQK, SMALL)
      RF(123) = 1.5D13
      EQK = EG(2)*EG(12)/EG(4)/EG(6)
      RB(123) = RF(123) / MAX(EQK, SMALL)
      RF(124) = 1.5D13
      EQK = EG(1)*EG(14)/EG(4)/EG(6)
      RB(124) = RF(124) / MAX(EQK, SMALL)
      RF(125) = 3D13
      EQK = EG(1)*EG(16)/EG(4)/EG(8)
      RB(125) = RF(125) / MAX(EQK, SMALL)
      RF(126) = 7D13
      EQK = EG(1)*EG(5)/EG(2)/EG(4)
      RB(126) = RF(126) / MAX(EQK, SMALL)
      RF(127) = 2.8D13
      EQK = EG(1)*EG(8)*EG(12)/EG(4)/EG(20)*PFAC1
      RB(127) = RF(127) / MAX(EQK, SMALL)
      RF(128) = 1.2D13
      EQK = EG(9)*EG(12)/EG(4)/EG(20)
      RB(128) = RF(128) / MAX(EQK, SMALL)
      RF(129) = 1.4D13
      EQK = EG(12)*EG(16)/EG(4)/EG(23)
      RB(129) = RF(129) / MAX(EQK, SMALL)
      RF(130) = EXP(3.64846865D1 -4.10997181D4*TI)
      EQK = EG(5)*EG(14)/EG(24)*PFAC1
      RB(130) = RF(130) / MAX(EQK, SMALL)
      RF(131) = EXP(9.72385961D1 -7.95359D0*ALOGT -4.61986113D4*TI)
      EQK = EG(5)*EG(19)/EG(27)*PFAC1
      RB(131) = RF(131) / MAX(EQK, SMALL)
      RF(132) = EXP(1.57191095D1 +2D0*ALOGT +3.16966115D2*TI)
      EQK = EG(9)*EG(25)/EG(8)/EG(27)
      RB(132) = RF(132) / MAX(EQK, SMALL)
      RF(133) = EXP(1.72066576D1 +2D0*ALOGT -2.02977978D3*TI)
      EQK = EG(2)*EG(25)/EG(1)/EG(27)
      RB(133) = RF(133) / MAX(EQK, SMALL)
      RF(134) = EXP(3.28840189D0 +3.7779D0*ALOGT -4.84663069D3*TI)
      EQK = EG(7)*EG(25)/EG(5)/EG(27)
      RB(134) = RF(134) / MAX(EQK, SMALL)
      RF(135) = EXP(-6.28987058D0 +5.29D0*ALOGT +5.48506168D1*TI)
      EQK = EG(8)*EG(25)/EG(6)/EG(27)
      RB(135) = RF(135) / MAX(EQK, SMALL)
      RF(136) = EXP(3.06267534D1 -8.30307502D3*TI)
      EQK = EG(22)*EG(25)/EG(21)/EG(27)
      RB(136) = RF(136) / MAX(EQK, SMALL)
      RF(137) = EXP(3.13445932D1 -2.25994605D4*TI)
      EQK = EG(21)*EG(25)/EG(20)/EG(27)
      RB(137) = RF(137) / MAX(EQK, SMALL)
      RF(138) = EXP(3.01159278D1 -1.29578292D4*TI)
      EQK = EG(5)*EG(16)/EG(25)*PFAC1
      RB(138) = RF(138) / MAX(EQK, SMALL)
      RF(139) = 2.41D13
      EQK = EG(16)*EG(27)/EG(19)/EG(25)
      RB(139) = RF(139) / MAX(EQK, SMALL)
      RF(140) = EXP(8.61068353D0 +2.8D0*ALOGT -2.94985611D3*TI)
      EQK = EG(14)*EG(27)/EG(16)/EG(25)
      RB(140) = RF(140) / MAX(EQK, SMALL)
      RF(141) = 9D12
      EQK = EG(8)*EG(31)/EG(21)/EG(25)
      RB(141) = RF(141) / MAX(EQK, SMALL)
      RF(142) = EXP(3.7398116D1 -6.6D-1*ALOGT -5.89769934D3*TI)
      EQK = EG(1)*EG(30)/EG(31)*PFAC1
      RB(142) = RF(142) / MAX(EQK, SMALL)
      RF(143) = EXP(2.99336062D1 -2.50098684D4*TI)
      EQK = EG(21)*EG(29)/EG(20)/EG(30)
      RB(143) = RF(143) / MAX(EQK, SMALL)
      RF(144) = EXP(1.69682466D1 +1.61D0*ALOGT +1.76125834D1*TI)
      EQK = EG(9)*EG(29)/EG(8)/EG(30)
      RB(144) = RF(144) / MAX(EQK, SMALL)
      RF(145) = EXP(2.7829872D1 -8.55468335D3*TI)
      EQK = EG(22)*EG(29)/EG(21)/EG(30)
      RB(145) = RF(145) / MAX(EQK, SMALL)
      RF(146) = EXP(1.23673408D1 +2.5D0*ALOGT -1.12217317D3*TI)
      EQK = EG(8)*EG(29)/EG(6)/EG(30)
      RB(146) = RF(146) / MAX(EQK, SMALL)
      RF(147) = EXP(1.53306378D1 +2D0*ALOGT -2.51608334D3*TI)
      EQK = EG(2)*EG(29)/EG(1)/EG(30)
      RB(147) = RF(147) / MAX(EQK, SMALL)
      RF(148) = EXP(-2.8103753D-1 +3.46D0*ALOGT -2.75813056D3*TI)
      EQK = EG(7)*EG(29)/EG(5)/EG(30)
      RB(148) = RF(148) / MAX(EQK, SMALL)
      RF(149) = EXP(2.96393694D1 -1.76D0*ALOGT -8.63016585D3*TI)
      EQK = EG(12)*EG(19)/EG(29)*PFAC1
      RB(149) = RF(149) / MAX(EQK, SMALL)
      RF(150) = EXP(2.80457763D1 -1.78D0*ALOGT -6.95445435D3*TI)
      EQK = EG(5)*EG(23)/EG(29)*PFAC1
      RB(150) = RF(150) / MAX(EQK, SMALL)
      RF(151) = 2D12
      EQK = EG(35)/EG(20)/EG(25)/PFAC1
      RB(151) = RF(151) / MAX(EQK, SMALL)
      RF(152) = EXP(5.3427584D1 -4.5D0*ALOGT)
      EQK = EG(20)*EG(31)*EG(31)/EG(35)/EG(35)*PFAC1
      RB(152) = RF(152) / MAX(EQK, SMALL)
      RF(153) = 4.28553538D-1*RF(152)
      EQK = EG(20)*EG(30)*EG(32)/EG(35)/EG(35)*PFAC1
      RB(153) = RF(153) / MAX(EQK, SMALL)
      RF(154) = EXP(3.68131678D1 -1.1D0*ALOGT -1.0386392D4*TI)
      EQK = EG(16)*EG(19)/EG(31)*PFAC1
      RB(154) = RF(154) / MAX(EQK, SMALL)
      RF(155) = EXP(2.46352888D1 -2.51608334D2*TI)
      EQK = EG(21)*EG(30)/EG(20)/EG(31)
      RB(155) = RF(155) / MAX(EQK, SMALL)
      RF(156) = EXP(2.48176104D1 -1.08191584D4*TI)
      EQK = EG(36)/EG(35)
      RB(156) = RF(156) / MAX(EQK, SMALL)
      RF(157) = EXP(3.03390713D1 -1.03159417D4*TI)
      EQK = EG(8)*EG(16)*EG(16)/EG(36)*PFAC2
      RB(157) = RF(157) / MAX(EQK, SMALL)
      RF(158) = 7D11
      EQK = EG(38)/EG(20)/EG(36)/PFAC1
      RB(158) = RF(158) / MAX(EQK, SMALL)
      RF(159) = EXP(2.44121453D1 -9.30950835D3*TI)
      EQK = EG(8)*EG(37)/EG(38)*PFAC1
      RB(159) = RF(159) / MAX(EQK, SMALL)
      RF(160) = EXP(3.79399738D1 -2.01286667D4*TI)
      EQK = EG(8)*EG(33)/EG(37)*PFAC1
      RB(160) = RF(160) / MAX(EQK, SMALL)
      RF(161) = EXP(2.5328436D1 -7.04503335D3*TI)
      EQK = EG(34)/EG(33)
      RB(161) = RF(161) / MAX(EQK, SMALL)
      RF(162) = EXP(3.76193093D1 -2.69D0*ALOGT -8.65532668D3*TI)
      EQK = EG(12)*EG(28)/EG(34)*PFAC1
      RB(162) = RF(162) / MAX(EQK, SMALL)
      RF(163) = EXP(3.62085565D1 -2.61D0*ALOGT -1.04719389D4*TI)
      EQK = EG(18)*EG(23)/EG(34)*PFAC1
      RB(163) = RF(163) / MAX(EQK, SMALL)
      RF(164) = EXP(3.22361913D1 -7.49792835D3*TI)
      EQK = EG(1)*EG(26)/EG(28)*PFAC1
      RB(164) = RF(164) / MAX(EQK, SMALL)
      RF(165) = EXP(3.60428538D1 -1.11D0*ALOGT)
      EQK = EG(28)/EG(8)/EG(16)/PFAC1
      RB(165) = RF(165) / MAX(EQK, SMALL)
      RF(166) = EXP(3.07665153D1 -2.51608334D4*TI)
      EQK = EG(9)*EG(12)/EG(26)*PFAC1
      RB(166) = RF(166) / MAX(EQK, SMALL)
      RF(167) = EXP(3.72468266D1 -2.86833501D4*TI)
      EQK = EG(2)*EG(23)/EG(26)*PFAC1
      RB(167) = RF(167) / MAX(EQK, SMALL)
      RF(168) = EXP(4.29710651D1 -4.6D-1*ALOGT -5.44983651D4*TI)
      EQK = EG(8)*EG(14)/EG(26)*PFAC1
      RB(168) = RF(168) / MAX(EQK, SMALL)
      RF(169) = EXP(1.47786849D1 +2.06D0*ALOGT -4.60946468D2*TI)
      EQK = EG(1)*EG(9)*EG(23)/EG(8)/EG(26)*PFAC1
      RB(169) = RF(169) / MAX(EQK, SMALL)
      RF(170) = EXP(1.67332813D1 +1.51D0*ALOGT +4.84094434D2*TI)
      EQK = EG(9)*EG(12)/EG(26)*PFAC1
      RB(170) = RF(170) / MAX(EQK, SMALL)
      RF(171) = EXP(1.52600738D1 +2.1D0*ALOGT -2.44965874D3*TI)
      EQK = EG(2)*EG(23)/EG(26)*PFAC1
      RB(171) = RF(171) / MAX(EQK, SMALL)
      RF(172) = EXP(3.17303532D1 -3.5D-1*ALOGT -1.5036114D3*TI)
      EQK = EG(2)*EG(8)*EG(12)/EG(1)/EG(26)*PFAC1
      RB(172) = RF(172) / MAX(EQK, SMALL)
      RF(173) = EXP(-1.47571191D1 +5.8D0*ALOGT -1.10707667D3*TI)
      EQK = EG(7)*EG(8)*EG(12)/EG(5)/EG(26)*PFAC1
      RB(173) = RF(173) / MAX(EQK, SMALL)
      RF(174) = EXP(2.76310211D1 -5.99834268D3*TI)
      EQK = EG(8)*EG(12)*EG(22)/EG(21)/EG(26)*PFAC1
      RB(174) = RF(174) / MAX(EQK, SMALL)
      RF(175) = EXP(4.20175112D1 -1.9D0*ALOGT -1.49706959D3*TI)
      EQK = EG(8)*EG(8)*EG(12)/EG(6)/EG(26)*PFAC1
      RB(175) = RF(175) / MAX(EQK, SMALL)
C
      RKLOW(1) = EXP(4.79026732D1 -1.72D0*ALOGT -2.64088107D2*TI)
      RKLOW(2) = EXP(3.93279334D1 -2.28963584D4*TI)
      RKLOW(3) = EXP(5.57002972D1 -2.79D0*ALOGT -2.10898105D3*TI)
      RKLOW(4) = EXP(7.34663067D1 -3.75D0*ALOGT -4.93957481D2*TI)
      RKLOW(5) = EXP(7.68923562D1 -4.76D0*ALOGT -1.22784867D3*TI)
      RKLOW(6) = EXP(9.50941235D1 -7.08D0*ALOGT -3.36400342D3*TI)
      RKLOW(7) = EXP(1.17075165D2 -9.31D0*ALOGT -5.02512164D4*TI)
      RKLOW(8) = EXP(9.68908955D1 -7.62D0*ALOGT -3.50742017D3*TI)
      RKLOW(9) = EXP(6.9414025D1 -3.86D0*ALOGT -1.67067934D3*TI)
      RKLOW(10) = EXP(9.34384048D1 -7.27D0*ALOGT -3.63322434D3*TI)
      RKLOW(11) = EXP(6.33329483D1 -3.14D0*ALOGT -6.18956501D2*TI)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RDSMH  (T, SMH)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION SMH(*), TN(5)
C
      TLOG = LOG(T)
      TI = 1D0/T
C
      TN(1) = TLOG - 1D0
      TN(2) = T
      TN(3) = TN(2)*T
      TN(4) = TN(3)*T
      TN(5) = TN(4)*T
C H
      IF (T .GT. 1000) THEN
      SMH(1) = -4.601176D-1 -2.547163D4*TI 
     *         +2.5D0*TN(1) 
      ELSE
      SMH(1) = -4.601176D-1 -2.547163D4*TI 
     *         +2.5D0*TN(1) 
      ENDIF
C H2
      IF (T .GT. 1000) THEN
      SMH(2) = -1.35511D0 +8.35034D2*TI 
     *         +2.991423D0*TN(1) +3.500322D-4*TN(2) 
     *         -9.389715D-9*TN(3) -7.69298167D-13*TN(4) 
     *         +7.91376D-17*TN(5) 
      ELSE
      SMH(2) = -3.294094D0 +1.012521D3*TI 
     *         +3.298124D0*TN(1) +4.124721D-4*TN(2) 
     *         -1.35716917D-7*TN(3) -7.896195D-12*TN(4) 
     *         +2.067436D-14*TN(5) 
      ENDIF
C CH2
      IF (T .GT. 1000) THEN
      SMH(3) = 6.17119324D0 -4.6263604D4*TI 
     *         +2.87410113D0*TN(1) +1.82819646D-3*TN(2) 
     *         -2.34824328D-7*TN(3) +2.16816291D-11*TN(4) 
     *         -9.38637835D-16*TN(5) 
      ELSE
      SMH(3) = 1.56253185D0 -4.60040401D4*TI 
     *         +3.76267867D0*TN(1) +4.84436072D-4*TN(2) 
     *         +4.65816402D-7*TN(3) -3.20909294D-10*TN(4) 
     *         +8.43708595D-14*TN(5) 
      ENDIF
C CH2(S)
      IF (T .GT. 1000) THEN
      SMH(4) = 8.62650169D0 -5.09259997D4*TI 
     *         +2.29203842D0*TN(1) +2.32794319D-3*TN(2) 
     *         -3.35319912D-7*TN(3) +3.48255D-11*TN(4) 
     *         -1.69858183D-15*TN(5) 
      ELSE
      SMH(4) = -7.69118967D-1 -5.04968163D4*TI 
     *         +4.19860411D0*TN(1) -1.1833071D-3*TN(2) 
     *         +1.37216037D-6*TN(3) -5.57346651D-10*TN(4) 
     *         +9.71573685D-14*TN(5) 
      ENDIF
C CH3
      IF (T .GT. 1000) THEN
      SMH(5) = 4.7224799D0 -1.6509513D4*TI 
     *         +2.9781206D0*TN(1) +2.898926D-3*TN(2) 
     *         -3.29263333D-7*TN(3) +2.56081583D-11*TN(4) 
     *         -8.958708D-16*TN(5) 
      ELSE
      SMH(5) = 1.6735354D0 -1.6422716D4*TI 
     *         +3.6571797D0*TN(1) +1.06329895D-3*TN(2) 
     *         +9.09731383D-7*TN(3) -5.51508358D-10*TN(4) 
     *         +1.2328537D-13*TN(5) 
      ENDIF
C O
      IF (T .GT. 1000) THEN
      SMH(6) = 4.920308D0 -2.92308D4*TI 
     *         +2.54206D0*TN(1) -1.377531D-5*TN(2) 
     *         -5.17133833D-10*TN(3) +3.79255583D-13*TN(4) 
     *         -2.184026D-17*TN(5) 
      ELSE
      SMH(6) = 2.963995D0 -2.914764D4*TI 
     *         +2.946429D0*TN(1) -8.19083D-4*TN(2) 
     *         +4.03505333D-7*TN(3) -1.3357025D-10*TN(4) 
     *         +1.945348D-14*TN(5) 
      ENDIF
C CH4
      IF (T .GT. 1000) THEN
      SMH(7) = 9.623395D0 +1.008079D4*TI 
     *         +1.683479D0*TN(1) +5.11862D-3*TN(2) 
     *         -6.45854833D-7*TN(3) +5.65465417D-11*TN(4) 
     *         -2.2517115D-15*TN(5) 
      ELSE
      SMH(7) = 1.372219D1 +9.825229D3*TI 
     *         +7.787415D-1*TN(1) +8.73834D-3*TN(2) 
     *         -4.639015D-6*TN(3) +2.54142333D-9*TN(4) 
     *         -6.119655D-13*TN(5) 
      ENDIF
C OH
      IF (T .GT. 1000) THEN
      SMH(8) = 5.70164073D0 -3.68362875D3*TI 
     *         +2.86472886D0*TN(1) +5.2825224D-4*TN(2) 
     *         -4.31804597D-8*TN(3) +2.54348895D-12*TN(4) 
     *         -6.6597938D-17*TN(5) 
      ELSE
      SMH(8) = -6.9043296D-1 -3.34630913D3*TI 
     *         +4.12530561D0*TN(1) -1.6127247D-3*TN(2) 
     *         +1.08794115D-6*TN(3) -4.83211369D-10*TN(4) 
     *         +1.03118689D-13*TN(5) 
      ENDIF
C H2O
      IF (T .GT. 1000) THEN
      SMH(9) = 6.862817D0 +2.989921D4*TI 
     *         +2.672146D0*TN(1) +1.5281465D-3*TN(2) 
     *         -1.45504333D-7*TN(3) +1.00083D-11*TN(4) 
     *         -3.195809D-16*TN(5) 
      ELSE
      SMH(9) = 2.590233D0 +3.020811D4*TI 
     *         +3.386842D0*TN(1) +1.737491D-3*TN(2) 
     *         -1.059116D-6*TN(3) +5.80715083D-10*TN(4) 
     *         -1.253294D-13*TN(5) 
      ENDIF
C C2H2
      IF (T .GT. 1000) THEN
      SMH(10) = -1.23028121D0 -2.59359992D4*TI 
     *         +4.14756964D0*TN(1) +2.98083332D-3*TN(2) 
     *         -3.9549142D-7*TN(3) +3.89510143D-11*TN(4) 
     *         -1.80617607D-15*TN(5) 
      ELSE
      SMH(10) = 1.39397051D1 -2.64289807D4*TI 
     *         +8.08681094D-1*TN(1) +1.16807815D-2*TN(2) 
     *         -5.91953025D-6*TN(3) +2.33460364D-9*TN(4) 
     *         -4.25036487D-13*TN(5) 
      ENDIF
C C2H3
      IF (T .GT. 1000) THEN
      SMH(11) = 7.78732378D0 -3.46128739D4*TI 
     *         +3.016724D0*TN(1) +5.1651146D-3*TN(2) 
     *         -7.80137248D-7*TN(3) +8.480274D-11*TN(4) 
     *         -4.3130352D-15*TN(5) 
      ELSE
      SMH(11) = 8.51054025D0 -3.48598468D4*TI 
     *         +3.21246645D0*TN(1) +7.5739581D-4*TN(2) 
     *         +4.32015687D-6*TN(3) -2.98048206D-9*TN(4) 
     *         +7.35754365D-13*TN(5) 
      ENDIF
C CO
      IF (T .GT. 1000) THEN
      SMH(12) = 6.108218D0 +1.426835D4*TI 
     *         +3.025078D0*TN(1) +7.213445D-4*TN(2) 
     *         -9.38471333D-8*TN(3) +8.488175D-12*TN(4) 
     *         -3.455476D-16*TN(5) 
      ELSE
      SMH(12) = 4.848897D0 +1.431054D4*TI 
     *         +3.262452D0*TN(1) +7.559705D-4*TN(2) 
     *         -6.46959167D-7*TN(3) +4.65162D-10*TN(4) 
     *         -1.2374755D-13*TN(5) 
      ENDIF
C C2H4
      IF (T .GT. 1000) THEN
      SMH(13) = 1.03053693D1 -4.93988614D3*TI 
     *         +2.03611116D0*TN(1) +7.32270755D-3*TN(2) 
     *         -1.11846319D-6*TN(3) +1.22685769D-10*TN(4) 
     *         -6.28530305D-15*TN(5) 
      ELSE
      SMH(13) = 4.09733096D0 -5.08977593D3*TI 
     *         +3.95920148D0*TN(1) -3.78526124D-3*TN(2) 
     *         +9.51650487D-6*TN(3) -5.76323961D-9*TN(4) 
     *         +1.34942187D-12*TN(5) 
      ENDIF
C HCO
      IF (T .GT. 1000) THEN
      SMH(14) = 5.552299D0 -3.916324D3*TI 
     *         +3.557271D0*TN(1) +1.6727865D-3*TN(2) 
     *         -2.22501D-7*TN(3) +2.05881083D-11*TN(4) 
     *         -8.569255D-16*TN(5) 
      ELSE
      SMH(14) = 8.983614D0 -4.159922D3*TI 
     *         +2.89833D0*TN(1) +3.0995735D-3*TN(2) 
     *         -1.60384733D-6*TN(3) +9.081875D-10*TN(4) 
     *         -2.2874425D-13*TN(5) 
      ENDIF
C C2H5
      IF (T .GT. 1000) THEN
      SMH(15) = 8.4602583D-1 -1.2056455D4*TI 
     *         +4.2878814D0*TN(1) +6.2169465D-3*TN(2) 
     *         -7.35651983D-7*TN(3) +5.88784183D-11*TN(4) 
     *         -2.1017568D-15*TN(5) 
      ELSE
      SMH(15) = 4.7100236D0 -1.2841714D4*TI 
     *         +4.305858D0*TN(1) -2.0916819D-3*TN(2) 
     *         +8.284545D-6*TN(3) -4.99215617D-9*TN(4) 
     *         +1.1524239D-12*TN(5) 
      ENDIF
C CH2O
      IF (T .GT. 1200) THEN
      SMH(16) = -5.1213813D0 +1.6230173D4*TI 
     *         +5.1481905D0*TN(1) +1.4339008D-3*TN(2) 
     *         -3.96377217D-8*TN(3) -1.34260858D-11*TN(4) 
     *         +1.42833675D-15*TN(5) 
      ELSE
      SMH(16) = 9.4697599D0 +1.4970793D4*TI 
     *         +2.6962612D0*TN(1) +2.46307115D-3*TN(2) 
     *         +1.38044157D-7*TN(3) -4.58651633D-11*TN(4) 
     *         -1.9805163D-14*TN(5) 
      ENDIF
C C2H6
      IF (T .GT. 1000) THEN
      SMH(17) = -5.239507D0 +1.271779D4*TI 
     *         +4.825938D0*TN(1) +6.920215D-3*TN(2) 
     *         -7.59543167D-7*TN(3) +5.60413917D-11*TN(4) 
     *         -1.7990805D-15*TN(5) 
      ELSE
      SMH(17) = 1.443229D1 +1.123918D4*TI 
     *         +1.462539D0*TN(1) +7.747335D-3*TN(2) 
     *         +9.63417833D-7*TN(3) -1.04819333D-9*TN(4) 
     *         +2.2931335D-13*TN(5) 
      ENDIF
C CH2OH
      IF (T .GT. 750) THEN
      SMH(18) = 5.4281095D0 +3.6664824D3*TI 
     *         +3.7469103D0*TN(1) +4.43230605D-3*TN(2) 
     *         -7.096787D-7*TN(3) +8.4067D-11*TN(4) 
     *         -4.72507805D-15*TN(5) 
      ELSE
      SMH(18) = 2.8351399D0 +3.6040734D3*TI 
     *         +4.6119792D0*TN(1) -1.560188D-3*TN(2) 
     *         +5.92194667D-6*TN(3) -4.11494983D-9*TN(4) 
     *         +1.10136235D-12*TN(5) 
      ENDIF
C CH3O
      IF (T .GT. 1000) THEN
      SMH(19) = 2.929575D0 -1.278325D2*TI 
     *         +3.7708D0*TN(1) +3.9357485D-3*TN(2) 
     *         -4.42730667D-7*TN(3) +3.28702583D-11*TN(4) 
     *         -1.056308D-15*TN(5) 
      ELSE
      SMH(19) = 1.315218D1 -9.786011D2*TI 
     *         +2.106204D0*TN(1) +3.6082975D-3*TN(2) 
     *         +8.89745333D-7*TN(3) -6.14803D-10*TN(4) 
     *         +1.0378055D-13*TN(5) 
      ENDIF
C O2
      IF (T .GT. 1000) THEN
      SMH(20) = 3.189166D0 +1.23393D3*TI 
     *         +3.697578D0*TN(1) +3.0675985D-4*TN(2) 
     *         -2.09807D-8*TN(3) +1.47940083D-12*TN(4) 
     *         -5.682175D-17*TN(5) 
      ELSE
      SMH(20) = 6.034738D0 +1.005249D3*TI 
     *         +3.212936D0*TN(1) +5.63743D-4*TN(2) 
     *         -9.59358333D-8*TN(3) +1.0948975D-10*TN(4) 
     *         -4.384277D-14*TN(5) 
      ENDIF
C HO2
      IF (T .GT. 1000) THEN
      SMH(21) = 3.78510215D0 -1.11856713D2*TI 
     *         +4.0172109D0*TN(1) +1.11991007D-3*TN(2) 
     *         -1.05609692D-7*TN(3) +9.52053083D-12*TN(4) 
     *         -5.39542675D-16*TN(5) 
      ELSE
      SMH(21) = 3.71666245D0 -2.9480804D2*TI 
     *         +4.30179801D0*TN(1) -2.37456026D-3*TN(2) 
     *         +3.52638152D-6*TN(3) -2.02303245D-9*TN(4) 
     *         +4.64612562D-13*TN(5) 
      ENDIF
C H2O2
      IF (T .GT. 1000) THEN
      SMH(22) = 5.01137D-1 +1.800696D4*TI 
     *         +4.573167D0*TN(1) +2.168068D-3*TN(2) 
     *         -2.457815D-7*TN(3) +1.95742D-11*TN(4) 
     *         -7.15827D-16*TN(5) 
      ELSE
      SMH(22) = 6.785363D0 +1.766315D4*TI 
     *         +3.388754D0*TN(1) +3.284613D-3*TN(2) 
     *         -2.47502167D-8*TN(3) -3.85483833D-10*TN(4) 
     *         +1.2357575D-13*TN(5) 
      ENDIF
C CO2
      IF (T .GT. 1000) THEN
      SMH(23) = -9.553959D-1 +4.896696D4*TI 
     *         +4.453623D0*TN(1) +1.5700845D-3*TN(2) 
     *         -2.130685D-7*TN(3) +1.9949975D-11*TN(4) 
     *         -8.345165D-16*TN(5) 
      ELSE
      SMH(23) = 1.018849D1 +4.837314D4*TI 
     *         +2.275725D0*TN(1) +4.961036D-3*TN(2) 
     *         -1.73485167D-6*TN(3) +5.72223917D-10*TN(4) 
     *         -1.05864D-13*TN(5) 
      ENDIF
C CH3HCO
      IF (T .GT. 1000) THEN
      SMH(24) = -3.4807917D0 +2.2593122D4*TI 
     *         +5.4041108D0*TN(1) +5.8615295D-3*TN(2) 
     *         -7.04385617D-7*TN(3) +5.69770425D-11*TN(4) 
     *         -2.04924315D-15*TN(5) 
      ELSE
      SMH(24) = 4.1030159D0 +2.1572878D4*TI 
     *         +4.7294595D0*TN(1) -1.5966429D-3*TN(2) 
     *         +7.92248683D-6*TN(3) -4.78821758D-9*TN(4) 
     *         +1.0965556D-12*TN(5) 
      ENDIF
C CH3OCH2
      IF (T .GT. 1376) THEN
      SMH(25) = -1.78650856D1 +3.41941605D3*TI 
     *         +8.17137842D0*TN(1) +5.50430905D-3*TN(2) 
     *         -6.37253795D-7*TN(3) +4.99697668D-11*TN(4) 
     *         -1.75158756D-15*TN(5) 
      ELSE
      SMH(25) = 1.16066817D1 +1.1884424D3*TI 
     *         +2.91327415D0*TN(1) +1.0168233D-2*TN(2) 
     *         -1.59952057D-6*TN(3) +1.72898771D-10*TN(4) 
     *         -8.5671681D-15*TN(5) 
      ENDIF
C HCOOH
      IF (T .GT. 1376) THEN
      SMH(26) = -1.13104798D1 +4.839954D4*TI 
     *         +6.68733013D0*TN(1) +2.57144684D-3*TN(2) 
     *         -3.03730855D-7*TN(3) +2.41432636D-11*TN(4) 
     *         -8.54460995D-16*TN(5) 
      ELSE
      SMH(26) = 1.72885798D1 +4.64616504D4*TI 
     *         +1.43548185D0*TN(1) +8.1681508D-3*TN(2) 
     *         -1.77095702D-6*TN(3) +2.76777481D-10*TN(4) 
     *         -2.01088051D-14*TN(5) 
      ENDIF
C CH3OCH3
      IF (T .GT. 710) THEN
      SMH(27) = 2.0217436D1 +2.34120975D4*TI 
     *         +8.30815546D-1*TN(1) +1.34586631D-2*TN(2) 
     *         -2.31457962D-6*TN(3) +2.89595899D-10*TN(4) 
     *         -1.70853392D-14*TN(5) 
      ELSE
      SMH(27) = -6.36955496D-1 +2.39755455D4*TI 
     *         +5.68097447D0*TN(1) -2.69717376D-3*TN(2) 
     *         +1.08245458D-5*TN(3) -6.70887765D-9*TN(4) 
     *         +1.63737009D-12*TN(5) 
      ENDIF
C HOCH2O
      IF (T .GT. 1452) THEN
      SMH(28) = -7.29290847D0 +2.47500385D4*TI 
     *         +6.39521515D0*TN(1) +3.71836521D-3*TN(2) 
     *         -4.1737059D-7*TN(3) +3.20733093D-11*TN(4) 
     *         -1.10889344D-15*TN(5) 
      ELSE
      SMH(28) = 6.81381989D0 +2.34414546D4*TI 
     *         +4.11183145D0*TN(1) +3.76925349D-3*TN(2) 
     *         +6.28895617D-7*TN(3) -4.48955004D-10*TN(4) 
     *         +7.28079435D-14*TN(5) 
      ENDIF
C CH3OCO
      IF (T .GT. 1362) THEN
      SMH(29) = -3.27914051D1 +2.466164D4*TI 
     *         +1.308776D1*TN(1) +2.26772475D-3*TN(2) 
     *         -2.75160607D-7*TN(3) +2.22664398D-11*TN(4) 
     *         -7.97884315D-16*TN(5) 
      ELSE
      SMH(29) = 1.66954362D1 +2.14404829D4*TI 
     *         +3.94199159D0*TN(1) +1.21717442D-2*TN(2) 
     *         -2.759926D-6*TN(3) +3.82114509D-10*TN(4) 
     *         -1.65897854D-14*TN(5) 
      ENDIF
C CH3OCHO
      IF (T .GT. 1686) THEN
      SMH(30) = -1.89301478D1 +4.64364769D4*TI 
     *         +8.69123518D0*TN(1) +5.7751561D-3*TN(2) 
     *         -7.1297081D-7*TN(3) +5.85444216D-11*TN(4) 
     *         -2.12166776D-15*TN(5) 
      ELSE
      SMH(30) = 1.25364719D1 +4.41855167D4*TI 
     *         +3.08839783D0*TN(1) +1.01880024D-2*TN(2) 
     *         -1.14129507D-6*TN(3) -6.06821836D-11*TN(4) 
     *         +2.81065108D-14*TN(5) 
      ENDIF
C CH3OCH2O
      IF (T .GT. 2012) THEN
      SMH(31) = -1.75775023D1 +2.13762444D4*TI 
     *         +8.60261845D0*TN(1) +6.78860975D-3*TN(2) 
     *         -8.07769337D-7*TN(3) +6.48138494D-11*TN(4) 
     *         -2.31316812D-15*TN(5) 
      ELSE
      SMH(31) = 1.23680069D1 +1.92377212D4*TI 
     *         +3.25889339D0*TN(1) +1.1107318D-2*TN(2) 
     *         -1.2975939D-6*TN(3) -2.01236798D-11*TN(4) 
     *         +2.25957248D-14*TN(5) 
      ENDIF
C CH3OCH2OH
      IF (T .GT. 2014) THEN
      SMH(32) = -1.80226702D1 +4.76607115D4*TI 
     *         +8.7098157D0*TN(1) +7.6801186D-3*TN(2) 
     *         -9.0167298D-7*TN(3) +7.17144538D-11*TN(4) 
     *         -2.54409876D-15*TN(5) 
      ELSE
      SMH(32) = 1.30511235D1 +4.54488899D4*TI 
     *         +3.15851876D0*TN(1) +1.22162876D-2*TN(2) 
     *         -1.44497464D-6*TN(3) -4.94432773D-12*TN(4) 
     *         +2.18200001D-14*TN(5) 
      ENDIF
C OCH2OCHO
      IF (T .GT. 1475) THEN
      SMH(33) = -3.33691809D1 +4.33647231D4*TI 
     *         +1.20233916D1*TN(1) +4.05631329D-3*TN(2) 
     *         -4.85594103D-7*TN(3) +3.8945032D-11*TN(4) 
     *         -1.38687762D-15*TN(5) 
      ELSE
      SMH(33) = 6.11645828D0 +4.02242792D4*TI 
     *         +5.19690837D0*TN(1) +7.94198615D-3*TN(2) 
     *         +5.89234245D-8*TN(3) -5.08714103D-10*TN(4) 
     *         +9.73309005D-14*TN(5) 
      ENDIF
C HOCH2OCO
      IF (T .GT. 1603) THEN
      SMH(34) = -2.86035265D1 +4.65575743D4*TI 
     *         +1.13737391D1*TN(1) +4.08831949D-3*TN(2) 
     *         -4.86723368D-7*TN(3) +3.88913013D-11*TN(4) 
     *         -1.38138411D-15*TN(5) 
      ELSE
      SMH(34) = 2.54054449D0 +4.39526183D4*TI 
     *         +6.08180801D0*TN(1) +6.43841795D-3*TN(2) 
     *         +3.4069903D-7*TN(3) -5.08462434D-10*TN(4) 
     *         +8.99102795D-14*TN(5) 
      ENDIF
C CH3OCH2O2
      IF (T .GT. 1389) THEN
      SMH(35) = -3.53740145D1 +2.29679238D4*TI 
     *         +1.24249729D1*TN(1) +5.9352993D-3*TN(2) 
     *         -6.7984422D-7*TN(3) +5.29425674D-11*TN(4) 
     *         -1.84713934D-15*TN(5) 
      ELSE
      SMH(35) = 1.91463601D1 +1.9494094D4*TI 
     *         +2.21029612D0*TN(1) +1.84438727D-2*TN(2) 
     *         -4.70935925D-6*TN(3) +9.64421108D-10*TN(4) 
     *         -9.8565235D-14*TN(5) 
      ENDIF
C CH2OCH2O2H
      IF (T .GT. 1393) THEN
      SMH(36) = -4.85706618D1 +1.84114867D4*TI 
     *         +1.51191783D1*TN(1) +4.61859441D-3*TN(2) 
     *         -5.31879175D-7*TN(3) +4.15928898D-11*TN(4) 
     *         -1.45581244D-15*TN(5) 
      ELSE
      SMH(36) = 1.76899251D1 +1.44293306D4*TI 
     *         +2.52895507D0*TN(1) +2.12064145D-2*TN(2) 
     *         -6.22343977D-6*TN(3) +1.38866111D-9*TN(4) 
     *         -1.48221656D-13*TN(5) 
      ENDIF
C HO2CH2OCHO
      IF (T .GT. 1387) THEN
      SMH(37) = -5.38924139D1 +6.23959608D4*TI 
     *         +1.64584298D1*TN(1) +4.26341756D-3*TN(2) 
     *         -5.06855833D-7*TN(3) +4.0466409D-11*TN(4) 
     *         -1.43658167D-15*TN(5) 
      ELSE
      SMH(37) = 1.52521392D1 +5.80629934D4*TI 
     *         +3.47935703D0*TN(1) +2.01476196D-2*TN(2) 
     *         -5.5018216D-6*TN(3) +1.11966764D-9*TN(4) 
     *         -1.0930079D-13*TN(5) 
      ENDIF
C O2CH2OCH2O2H
      IF (T .GT. 1402) THEN
      SMH(38) = -6.51847273D1 +3.79207055D4*TI 
     *         +1.92038046D1*TN(1) +5.21974205D-3*TN(2) 
     *         -6.00971565D-7*TN(3) +4.69827369D-11*TN(4) 
     *         -1.64403607D-15*TN(5) 
      ELSE
      SMH(38) = 2.44215005D1 +3.27628742D4*TI 
     *         +1.99640551D0*TN(1) +2.91613116D-2*TN(2) 
     *         -9.2209963D-6*TN(3) +2.16508783D-9*TN(4) 
     *         -2.38570502D-13*TN(5) 
      ENDIF
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RATX (T, C, RF, RB, RKLOW)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (SMALL = 1D-200)
      DIMENSION C(*), RF(*), RB(*), RKLOW(*)
C
      ALOGT = LOG(T)
      CTOT = 0.0
      DO K = 1, 30
         CTOT = CTOT + C(K)
      ENDDO
C
      RF(1) = RF(1)*C(1)*C(15)
      RB(1) = RB(1)*C(4)*C(6)
      RF(2) = RF(2)*C(2)*C(4)
      RB(2) = RB(2)*C(1)*C(6)
      RF(3) = RF(3)*C(2)*C(6)
      RB(3) = RB(3)*C(1)*C(7)
      RF(4) = RF(4)*C(4)*C(7)
      RB(4) = RB(4)*C(6)*C(6)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(7)+9D-1*C(9)+2.8D0*C(18)
      RF(5) = RF(5)*CTB*C(2)
      RB(5) = RB(5)*CTB*C(1)*C(1)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(7)+9D-1*C(9)+2.8D0*C(18)
      RF(6) = RF(6)*CTB*C(4)*C(4)
      RB(6) = RB(6)*CTB*C(15)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(7)+9D-1*C(9)+2.8D0*C(18)
      RF(7) = RF(7)*CTB*C(1)*C(4)
      RB(7) = RB(7)*CTB*C(6)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(7)+9D-1*C(9)+2.8D0*C(18)
      RF(8) = RF(8)*CTB*C(1)*C(6)
      RB(8) = RB(8)*CTB*C(7)
      CTB = CTOT+C(2)+1D1*C(7)-2.2D-1*C(15)+9D-1*C(9)
     * +2.8D0*C(18)
      PR = RKLOW(1) * CTB / RF(9)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FC = (PRLOG -3.35070291D-1)
     *     /(8.73075717D-1 -0.14D0*(PRLOG -3.35070291D-1))
      FC = EXP(-2.23143551D-1 /(1.0D0 + FC*FC))
      PCOR = FC * PCOR
      RF(9) = RF(9) * PCOR
      RB(9) = RB(9) * PCOR
      RF(9) = RF(9)*C(1)*C(15)
      RB(9) = RB(9)*C(16)
      RF(10) = RF(10)*C(1)*C(16)
      RB(10) = RB(10)*C(2)*C(15)
      RF(11) = RF(11)*C(1)*C(16)
      RB(11) = RB(11)*C(6)*C(6)
      RF(12) = RF(12)*C(4)*C(16)
      RB(12) = RB(12)*C(6)*C(15)
      RF(13) = RF(13)*C(6)*C(16)
      RB(13) = RB(13)*C(7)*C(15)
      RF(14) = RF(14)*C(16)*C(16)
      RB(14) = RB(14)*C(15)*C(17)
      RF(15) = RF(15)*C(16)*C(16)
      RB(15) = RB(15)*C(15)*C(17)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(7)+9D-1*C(9)+2.8D0*C(18)
      PR = RKLOW(2) * CTB / RF(16)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FC = (PRLOG -1.98309903D-1)
     *     /(1.13230809D0 -0.14D0*(PRLOG -1.98309903D-1))
      FC = EXP(-6.93147181D-1 /(1.0D0 + FC*FC))
      PCOR = FC * PCOR
      RF(16) = RF(16) * PCOR
      RB(16) = RB(16) * PCOR
      RF(16) = RF(16)*C(17)
      RB(16) = RB(16)*C(6)*C(6)
      RF(17) = RF(17)*C(1)*C(17)
      RB(17) = RB(17)*C(6)*C(7)
      RF(18) = RF(18)*C(1)*C(17)
      RB(18) = RB(18)*C(2)*C(16)
      RF(19) = RF(19)*C(4)*C(17)
      RB(19) = RB(19)*C(6)*C(16)
      RF(20) = RF(20)*C(6)*C(17)
      RB(20) = RB(20)*C(7)*C(16)
      RF(21) = RF(21)*C(6)*C(17)
      RB(21) = RB(21)*C(7)*C(16)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(7)+9D-1*C(9)+2.8D0*C(18)
      PR = RKLOW(3) * CTB / RF(22)
      PCOR = PR / (1.0 + PR)
      RF(22) = RF(22) * PCOR
      RB(22) = RB(22) * PCOR
      RF(22) = RF(22)*C(4)*C(9)
      RB(22) = RB(22)*C(18)
      RF(23) = RF(23)*C(9)*C(15)
      RB(23) = RB(23)*C(4)*C(18)
      RF(24) = RF(24)*C(9)*C(16)
      RB(24) = RB(24)*C(6)*C(18)
      RF(25) = RF(25)*C(6)*C(9)
      RB(25) = RB(25)*C(1)*C(18)
      CTB = CTOT+1.5D0*C(2)+5D0*C(7)+9D-1*C(9)+2.8D0*C(18)
      RF(26) = RF(26)*CTB
      RB(26) = RB(26)*CTB*C(1)*C(9)
      RF(27) = RF(27)*C(15)
      RB(27) = RB(27)*C(9)*C(16)
      RF(28) = RF(28)*C(1)
      RB(28) = RB(28)*C(2)*C(9)
      RF(29) = RF(29)*C(4)
      RB(29) = RB(29)*C(6)*C(9)
      RF(30) = RF(30)*C(6)
      RB(30) = RB(30)*C(7)*C(9)
      RF(31) = RF(31)*C(4)
      RB(31) = RB(31)*C(1)*C(18)
      RF(32) = RF(32)*C(16)
      RB(32) = RB(32)*C(1)*C(6)*C(18)
      RB(33) = RB(33)*C(2)*C(9)*C(9)
      RF(34) = RF(34)*C(3)
      RB(34) = RB(34)*C(5)*C(9)
      RB(35) = RB(35)*C(9)*C(12)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(7)+9D-1*C(9)+2.8D0*C(18)
      RF(36) = RF(36)*CTB*C(12)
      RB(36) = RB(36)*CTB*C(1)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(7)+9D-1*C(9)+2.8D0*C(18)
      RF(37) = RF(37)*CTB*C(12)
      RB(37) = RB(37)*CTB*C(2)*C(9)
      RF(38) = RF(38)*C(1)*C(12)
      RB(38) = RB(38)*C(2)
      RF(39) = RF(39)*C(4)*C(12)
      RB(39) = RB(39)*C(6)
      RF(40) = RF(40)*C(6)*C(12)
      RB(40) = RB(40)*C(7)
      RF(41) = RF(41)*C(12)*C(15)
      RB(41) = RB(41)*C(16)
      RF(42) = RF(42)*C(12)*C(16)
      RB(42) = RB(42)*C(17)
      RF(43) = RF(43)*C(3)*C(12)
      RB(43) = RB(43)*C(5)
      RF(44) = RF(44)*C(3)*C(4)
      RB(44) = RB(44)*C(1)*C(12)
      RF(45) = RF(45)*C(3)*C(15)
      RB(45) = RB(45)*C(4)*C(14)
      RF(46) = RF(46)*C(3)*C(15)
      RB(46) = RB(46)*C(6)*C(12)
      RF(47) = RF(47)*C(3)*C(16)
      RB(47) = RB(47)*C(6)*C(14)
      CTB = CTOT+4D0*C(7)+C(9)+2D0*C(18)
      PR = RKLOW(4) * CTB / RF(48)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 1D0*EXP(-T/5.7D2)
     *     + EXP(-1D30/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(48) = RF(48) * PCOR
      RB(48) = RB(48) * PCOR
      RF(48) = RF(48)*C(3)*C(3)
      RB(48) = RB(48)*C(13)
      CTB = CTOT+C(2)+5D0*C(7)+C(5)+5D-1*C(9)
     * +C(18)+2D0*C(13)
      PR = RKLOW(5) * CTB / RF(49)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.17D-1*EXP(-T/7.4D1) + 7.83D-1*EXP(-T/2.941D3)
     *     + EXP(-6.964D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(49) = RF(49) * PCOR
      RB(49) = RB(49) * PCOR
      RF(49) = RF(49)*C(1)*C(3)
      RB(49) = RB(49)*C(5)
      RF(50) = RF(50)*C(1)*C(5)
      RB(50) = RB(50)*C(2)*C(3)
      RF(51) = RF(51)*C(4)*C(5)
      RB(51) = RB(51)*C(3)*C(6)
      RF(52) = RF(52)*C(5)*C(6)
      RB(52) = RB(52)*C(3)*C(7)
      RF(53) = RF(53)*C(3)*C(16)
      RB(53) = RB(53)*C(5)*C(15)
      RF(54) = RF(54)*C(5)*C(16)
      RB(54) = RB(54)*C(3)*C(17)
      CTB = CTOT
      RF(55) = RF(55)*CTB
      RB(55) = RB(55)*CTB*C(1)*C(12)
      RF(56) = RF(56)*C(1)
      RB(56) = RB(56)*C(2)*C(12)
      RF(57) = RF(57)*C(1)
      RB(57) = RB(57)*C(3)*C(6)
      RF(58) = RF(58)*C(4)
      RB(58) = RB(58)*C(6)*C(12)
      RF(59) = RF(59)*C(6)
      RB(59) = RB(59)*C(7)*C(12)
      RF(60) = RF(60)*C(15)
      RB(60) = RB(60)*C(12)*C(16)
      RF(61) = RF(61)*C(15)
      RB(61) = RB(61)*C(12)*C(16)
      RF(62) = RF(62)*C(16)
      RB(62) = RB(62)*C(12)*C(17)
      RB(63) = RB(63)*C(12)*C(12)
      CTB = CTOT
      RF(64) = RF(64)*CTB*C(14)
      RB(64) = RB(64)*CTB*C(1)*C(12)
      RF(65) = RF(65)*C(1)*C(14)
      RB(65) = RB(65)*C(3)*C(6)
      RF(66) = RF(66)*C(4)*C(14)
      RB(66) = RB(66)*C(6)*C(12)
      RF(67) = RF(67)*C(6)*C(14)
      RB(67) = RB(67)*C(7)*C(12)
      RF(68) = RF(68)*C(14)*C(15)
      RB(68) = RB(68)*C(12)*C(16)
      RF(69) = RF(69)*C(14)*C(15)
      RB(69) = RB(69)*C(12)*C(16)
      RF(70) = RF(70)*C(14)*C(16)
      RB(70) = RB(70)*C(12)*C(17)
      RF(71) = RF(71)*C(9)*C(14)
      RB(71) = RB(71)*C(3)*C(18)
      RF(72) = RF(72)*C(3)*C(3)
      RB(72) = RB(72)*C(1)*C(11)
      RF(73) = RF(73)*C(5)
      RB(73) = RB(73)*C(3)*C(3)
      RF(74) = RF(74)*C(5)
      RB(74) = RB(74)*C(3)*C(3)
      RF(75) = RF(75)*C(3)*C(6)
      RB(75) = RB(75)*C(7)
      RF(76) = RF(76)*C(3)*C(6)
      RB(76) = RB(76)*C(7)
      RF(77) = RF(77)*C(3)
      RB(77) = RB(77)*C(1)*C(10)
      RF(78) = RF(78)*C(3)
      RB(78) = RB(78)*C(1)*C(10)
      RF(79) = RF(79)*C(1)*C(14)
      RB(79) = RB(79)*C(7)
      RF(80) = RF(80)*C(1)*C(13)
      RB(80) = RB(80)*C(2)*C(11)
      RF(81) = RF(81)*C(4)*C(13)
      RB(81) = RB(81)*C(6)*C(11)
      RF(82) = RF(82)*C(6)*C(13)
      RB(82) = RB(82)*C(7)*C(11)
      RF(83) = RF(83)*C(13)*C(15)
      RB(83) = RB(83)*C(11)*C(16)
      RF(84) = RF(84)*C(13)*C(16)
      RB(84) = RB(84)*C(11)*C(17)
      RF(85) = RF(85)*C(3)*C(13)
      RB(85) = RB(85)*C(5)*C(11)
      CTB = CTOT+C(2)+5D0*C(7)+C(5)+5D-1*C(9)
     * +C(18)+2D0*C(13)
      PR = RKLOW(6) * CTB / RF(86)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 1.578D-1*EXP(-T/1.25D2) + 8.422D-1*EXP(-T/2.219D3)
     *     + EXP(-6.882D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(86) = RF(86) * PCOR
      RB(86) = RB(86) * PCOR
      RF(86) = RF(86)*C(1)*C(11)
      RB(86) = RB(86)*C(13)
      RF(87) = RF(87)*C(1)*C(11)
      RB(87) = RB(87)*C(2)*C(10)
      RF(88) = RF(88)*C(4)*C(11)
      RB(88) = RB(88)*C(3)*C(12)
      RF(89) = RF(89)*C(11)*C(15)
      RB(89) = RB(89)*C(10)*C(16)
      RF(90) = RF(90)*C(11)*C(11)
      RB(90) = RB(90)*C(10)*C(13)
      RF(91) = RF(91)*C(11)
      RB(91) = RB(91)*C(9)*C(13)
      RF(92) = RF(92)*C(4)*C(11)
      RB(92) = RB(92)*C(1)*C(19)
      CTB = CTOT+C(2)+5D0*C(7)+C(5)+5D-1*C(9)
     * +C(18)+2D0*C(13)
      PR = RKLOW(7) * CTB / RF(93)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.655D-1*EXP(-T/1.8D2) + 7.345D-1*EXP(-T/1.035D3)
     *     + EXP(-5.417D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(93) = RF(93) * PCOR
      RB(93) = RB(93) * PCOR
      RF(93) = RF(93)*C(10)
      RB(93) = RB(93)*C(2)*C(8)
      CTB = CTOT+C(2)+5D0*C(7)+C(5)+5D-1*C(9)
     * +C(18)+2D0*C(13)
      PR = RKLOW(8) * CTB / RF(94)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.47D-2*EXP(-T/2.1D2) + 9.753D-1*EXP(-T/9.84D2)
     *     + EXP(-4.374D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(94) = RF(94) * PCOR
      RB(94) = RB(94) * PCOR
      RF(94) = RF(94)*C(1)*C(10)
      RB(94) = RB(94)*C(11)
      CTB = CTOT+C(2)+5D0*C(7)+C(5)+5D-1*C(9)
     * +C(18)+2D0*C(13)
      PR = RKLOW(9) * CTB / RF(95)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.18D-1*EXP(-T/2.075D2) + 7.82D-1*EXP(-T/2.663D3)
     *     + EXP(-6.095D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(95) = RF(95) * PCOR
      RB(95) = RB(95) * PCOR
      RF(95) = RF(95)*C(1)
      RB(95) = RB(95)*C(10)
      RF(96) = RF(96)*C(1)*C(10)
      RB(96) = RB(96)*C(2)
      RF(97) = RF(97)*C(6)*C(10)
      RB(97) = RB(97)*C(7)
      RF(98) = RF(98)*C(3)*C(10)
      RB(98) = RB(98)*C(5)
      RF(99) = RF(99)*C(4)*C(10)
      RB(99) = RB(99)*C(3)
      RF(100) = RF(100)*C(6)
      RB(100) = RB(100)*C(7)*C(8)
      RF(101) = RF(101)*C(4)*C(10)
      RB(101) = RB(101)*C(6)
      RF(102) = RF(102)*C(10)*C(15)
      RB(102) = RB(102)*C(16)
      RF(103) = RF(103)*C(1)
      RB(103) = RB(103)*C(2)*C(8)
      RF(104) = RF(104)*C(17)
      RB(104) = RB(104)*C(10)*C(16)
      RF(105) = RF(105)*C(3)
      RB(105) = RB(105)*C(5)*C(8)
      RB(106) = RB(106)*C(8)*C(10)
      RF(107) = RF(107)*C(15)
      RB(107) = RB(107)*C(12)
      RF(108) = RF(108)*C(15)
      RB(108) = RB(108)*C(8)*C(16)
      CTB = CTOT+C(2)+5D0*C(7)+C(5)+5D-1*C(9)
     * +C(18)+2D0*C(13)
      PR = RKLOW(10) * CTB / RF(109)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.493D-1*EXP(-T/9.85D1) + 7.507D-1*EXP(-T/1.302D3)
     *     + EXP(-4.167D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(109) = RF(109) * PCOR
      RB(109) = RB(109) * PCOR
      RF(109) = RF(109)*C(1)*C(8)
      RF(110) = RF(110)*C(4)*C(8)
      RB(110) = RB(110)*C(9)
      RF(111) = RF(111)*C(6)*C(8)
      RB(111) = RB(111)*C(3)*C(9)
      CTB = CTOT+C(2)+5D0*C(7)+C(5)+5D-1*C(9)
     * +C(18)+2D0*C(13)
      PR = RKLOW(11) * CTB / RF(112)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 3.2D-1*EXP(-T/7.8D1) + 6.8D-1*EXP(-T/1.995D3)
     *     + EXP(-5.59D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(112) = RF(112) * PCOR
      RB(112) = RB(112) * PCOR
      RF(112) = RF(112)*C(1)
      RB(112) = RB(112)*C(3)
      RF(113) = RF(113)*C(4)
      RB(113) = RB(113)*C(1)
      RF(114) = RF(114)*C(6)
      RB(114) = RB(114)*C(1)*C(12)
      RF(115) = RF(115)*C(2)
      RB(115) = RB(115)*C(1)*C(3)
      RF(116) = RF(116)*C(15)
      RB(116) = RB(116)*C(6)
      RF(117) = RF(117)*C(16)
      RB(117) = RB(117)*C(6)*C(12)
      RB(118) = RB(118)*C(2)*C(8)
      CTB = CTOT-C(7)-C(9)-C(18)
      RF(119) = RF(119)*CTB
      RB(119) = RB(119)*CTB
      RF(120) = RF(120)*C(7)
      RB(120) = RB(120)*C(7)
      RF(121) = RF(121)*C(9)
      RB(121) = RB(121)*C(9)
      RF(122) = RF(122)*C(18)
      RB(122) = RB(122)*C(18)
      RF(123) = RF(123)*C(4)
      RB(123) = RB(123)*C(2)*C(9)
      RF(124) = RF(124)*C(4)
      RB(124) = RB(124)*C(1)
      RF(125) = RF(125)*C(6)
      RB(125) = RB(125)*C(1)*C(12)
      RF(126) = RF(126)*C(2)
      RB(126) = RB(126)*C(1)*C(3)
      RF(127) = RF(127)*C(15)
      RB(127) = RB(127)*C(1)*C(6)*C(9)
      RF(128) = RF(128)*C(15)
      RB(128) = RB(128)*C(7)*C(9)
      RF(129) = RF(129)*C(18)
      RB(129) = RB(129)*C(9)*C(12)
      RF(130) = RF(130)*C(19)
      RB(130) = RB(130)*C(3)
      RF(131) = RF(131)*C(21)
      RB(131) = RB(131)*C(3)*C(14)
      RF(132) = RF(132)*C(6)*C(21)
      RB(132) = RB(132)*C(7)
      RF(133) = RF(133)*C(1)*C(21)
      RB(133) = RB(133)*C(2)
      RF(134) = RF(134)*C(3)*C(21)
      RB(134) = RB(134)*C(5)
      RF(135) = RF(135)*C(4)*C(21)
      RB(135) = RB(135)*C(6)
      RF(136) = RF(136)*C(16)*C(21)
      RB(136) = RB(136)*C(17)
      RF(137) = RF(137)*C(15)*C(21)
      RB(137) = RB(137)*C(16)
      RB(138) = RB(138)*C(3)*C(12)
      RF(139) = RF(139)*C(14)
      RB(139) = RB(139)*C(12)*C(21)
      RF(140) = RF(140)*C(12)
      RB(140) = RB(140)*C(21)
      RF(141) = RF(141)*C(16)
      RB(141) = RB(141)*C(6)
      RB(142) = RB(142)*C(1)*C(23)
      RF(143) = RF(143)*C(15)*C(23)
      RB(143) = RB(143)*C(16)*C(22)
      RF(144) = RF(144)*C(6)*C(23)
      RB(144) = RB(144)*C(7)*C(22)
      RF(145) = RF(145)*C(16)*C(23)
      RB(145) = RB(145)*C(17)*C(22)
      RF(146) = RF(146)*C(4)*C(23)
      RB(146) = RB(146)*C(6)*C(22)
      RF(147) = RF(147)*C(1)*C(23)
      RB(147) = RB(147)*C(2)*C(22)
      RF(148) = RF(148)*C(3)*C(23)
      RB(148) = RB(148)*C(5)*C(22)
      RF(149) = RF(149)*C(22)
      RB(149) = RB(149)*C(9)*C(14)
      RF(150) = RF(150)*C(22)
      RB(150) = RB(150)*C(3)*C(18)
      RF(151) = RF(151)*C(15)
      RB(151) = RB(151)*C(27)
      RF(152) = RF(152)*C(27)*C(27)
      RB(152) = RB(152)*C(15)
      RF(153) = RF(153)*C(27)*C(27)
      RB(153) = RB(153)*C(15)*C(23)*C(24)
      RB(154) = RB(154)*C(12)*C(14)
      RF(155) = RF(155)*C(15)
      RB(155) = RB(155)*C(16)*C(23)
      RF(156) = RF(156)*C(27)
      RB(157) = RB(157)*C(6)*C(12)*C(12)
      RF(158) = RF(158)*C(15)
      RB(158) = RB(158)*C(29)
      RF(159) = RF(159)*C(29)
      RB(159) = RB(159)*C(6)*C(28)
      RF(160) = RF(160)*C(28)
      RB(160) = RB(160)*C(6)*C(25)
      RF(161) = RF(161)*C(25)
      RB(161) = RB(161)*C(26)
      RF(162) = RF(162)*C(26)
      RB(162) = RB(162)*C(9)
      RF(163) = RF(163)*C(26)
      RB(163) = RB(163)*C(18)
      RB(164) = RB(164)*C(1)*C(20)
      RF(165) = RF(165)*C(6)*C(12)
      CTB = CTOT
      RF(166) = RF(166)*CTB*C(20)
      RB(166) = RB(166)*CTB*C(7)*C(9)
      CTB = CTOT
      RF(167) = RF(167)*CTB*C(20)
      RB(167) = RB(167)*CTB*C(2)*C(18)
      RF(168) = RF(168)*C(20)
      RB(168) = RB(168)*C(6)
      RF(169) = RF(169)*C(6)*C(20)
      RB(169) = RB(169)*C(1)*C(7)*C(18)
      RF(170) = RF(170)*C(6)*C(20)
      RB(170) = RB(170)*C(6)*C(7)*C(9)
      RF(171) = RF(171)*C(1)*C(20)
      RB(171) = RB(171)*C(1)*C(2)*C(18)
      RF(172) = RF(172)*C(1)*C(20)
      RB(172) = RB(172)*C(2)*C(6)*C(9)
      RF(173) = RF(173)*C(3)*C(20)
      RB(173) = RB(173)*C(5)*C(6)*C(9)
      RF(174) = RF(174)*C(16)*C(20)
      RB(174) = RB(174)*C(6)*C(9)*C(17)
      RF(175) = RF(175)*C(4)*C(20)
      RB(175) = RB(175)*C(6)*C(6)*C(9)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE QSSA(RF, RB, XQ)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RF(*), RB(*), XQ(*)
      PARAMETER (SMALL = 1D-50)
C
      RF(33) = 0.D0
      RF(35) = 0.D0
      RF(63) = 0.D0
      RF(106) = 0.D0
      RF(118) = 0.D0
      RB(152) = 0.D0
C
C     CH2
      DEN = +RF( 73) +RF( 77) +RF(112) +RF(113) +RF(114) 
     *  +RF(115) +RF(116) +RF(117) +RB( 75) +RB(110) +RB(119) 
     *  +RB(120) +RB(121) +RB(122) 
      A1_0 = ( +RB( 73) +RF( 75) +RB( 77) +RF(110) +RB(112) 
     *  +RB(114) +RB(115) +RB(117) +RB(118) +RB(118) )/MAX(DEN,SMALL)
      A1_2 = ( +RF(119) +RF(120) +RF(121) +RF(122) )/MAX(DEN,SMALL)
      A1_4 = ( +RB(113) +RB(116) )/MAX(DEN,SMALL)
C     CH2(S)
      DEN = +RF( 74) +RF( 78) +RF(119) +RF(120) +RF(121) 
     *  +RF(122) +RF(123) +RF(124) +RF(125) +RF(126) +RF(127) 
     *  +RF(128) +RF(129) +RB( 76) +RB( 79) 
      A2_0 = ( +RB( 74) +RF( 76) +RB( 78) +RF( 79) +RB(123) 
     *  +RB(125) +RB(126) +RB(127) +RB(128) +RB(129) )/MAX(DEN,SMALL)
      A2_1 = ( +RB(119) +RB(120) +RB(121) +RB(122) )/MAX(DEN,SMALL)
      A2_4 = ( +RB(124) )/MAX(DEN,SMALL)
C     C2H3
      DEN = +RF( 95) +RF(100) +RF(103) +RF(104) +RF(105) 
     *  +RF(107) +RF(108) +RB( 96) +RB( 97) +RB( 98) +RB(101) 
     *  +RB(102) +RB(109) 
      A3_0 = ( +RB( 95) +RF( 96) +RF( 97) +RF( 98) +RB(100) 
     *  +RF(101) +RF(102) +RB(103) +RB(104) +RB(105) +RB(106) 
     *  +RB(106) +RB(108) +RF(109) )/MAX(DEN,SMALL)
      A3_4 = ( +RB(107) )/MAX(DEN,SMALL)
C     HCO
      DEN = +RF( 26) +RF( 27) +RF( 28) +RF( 29) +RF( 30) 
     *  +RF( 31) +RF( 32) +RF( 34) +RF( 91) +RB( 36) +RB( 38) 
     *  +RB( 39) +RB( 40) +RB( 41) +RB( 42) +RB( 43) +RB( 99) 
     *  +RB(107) +RB(113) +RB(116) +RB(124) +RB(130) +RB(140) 
     *  +RB(168) 
      A4_0 = ( +RB( 26) +RB( 27) +RB( 28) +RB( 29) +RB( 30) 
     *  +RB( 31) +RB( 32) +RB( 33) +RB( 33) +RB( 34) +RB( 35) 
     *  +RB( 35) +RF( 36) +RF( 38) +RF( 39) +RF( 40) +RF( 41) 
     *  +RF( 42) +RF( 43) +RB( 63) +RB( 91) +RF( 99) +RF(130) 
     *  +RF(168) )/MAX(DEN,SMALL)
      A4_1 = ( +RF(113) +RF(116) )/MAX(DEN,SMALL)
      A4_2 = ( +RF(124) )/MAX(DEN,SMALL)
      A4_3 = ( +RF(107) )/MAX(DEN,SMALL)
      A4_6 = ( +RF(140) )/MAX(DEN,SMALL)
C     CH2OH
      DEN = +RF( 55) +RF( 56) +RF( 57) +RF( 58) +RF( 59) 
     *  +RF( 60) +RF( 61) +RF( 62) +RB(163) 
      A5_0 = ( +RB( 55) +RB( 56) +RB( 57) +RB( 58) +RB( 59) 
     *  +RB( 60) +RB( 61) +RB( 62) +RB( 63) +RF(163) )/MAX(DEN,SMALL)
C     CH3OCH2
      DEN = +RF(138) +RF(139) +RF(140) +RF(141) +RF(151) 
     *  +RB(132) +RB(133) +RB(134) +RB(135) +RB(136) +RB(137) 
      A6_0 = ( +RF(132) +RF(133) +RF(134) +RF(135) +RF(136) 
     *  +RF(137) +RB(138) +RB(139) +RB(151) )/MAX(DEN,SMALL)
      A6_4 = ( +RB(140) )/MAX(DEN,SMALL)
      A6_8 = ( +RB(141) )/MAX(DEN,SMALL)
C     HOCH2O
      DEN = +RF(164) +RB(162) +RB(165) 
      A7_0 = ( +RF(162) +RB(164) +RF(165) )/MAX(DEN,SMALL)
C     CH3OCH2O
      DEN = +RF(142) +RF(154) +RF(155) +RB(141) 
      A8_0 = ( +RB(142) +RF(152) +RF(152) +RB(154) +RB(155) )
     *          /MAX(DEN,SMALL)
      A8_6 = ( +RF(141) )/MAX(DEN,SMALL)
C     CH2OCH2O2H
      DEN = +RF(157) +RF(158) +RB(156) 
      A9_0 = ( +RF(156) +RB(157) +RB(158) )/MAX(DEN,SMALL)
C
      A6_0 = A6_0 + A6_8*A8_0
      DEN = 1 -A6_8*A8_6
      A6_0 = A6_0/MAX(DEN,SMALL)
      A6_4 = A6_4/MAX(DEN,SMALL)
      A4_0 = A4_0 + A4_3*A3_0
      DEN = 1 -A4_3*A3_4
      A4_0 = A4_0/MAX(DEN,SMALL)
      A4_1 = A4_1/MAX(DEN,SMALL)
      A4_6 = A4_6/MAX(DEN,SMALL)
      A4_2 = A4_2/MAX(DEN,SMALL)
      A4_0 = A4_0 + A4_2*A2_0
      A4_1 = A4_1 + A4_2*A2_1
      DEN = 1 -A4_2*A2_4
      A4_0 = A4_0/MAX(DEN,SMALL)
      A4_1 = A4_1/MAX(DEN,SMALL)
      A4_6 = A4_6/MAX(DEN,SMALL)
      A1_0 = A1_0 + A1_2*A2_0
      A1_4 = A1_4 + A1_2*A2_4
      DEN = 1 -A1_2*A2_1
      A1_0 = A1_0/MAX(DEN,SMALL)
      A1_4 = A1_4/MAX(DEN,SMALL)
      A4_0 = A4_0 + A4_6*A6_0
      DEN = 1 -A4_6*A6_4
      A4_0 = A4_0/MAX(DEN,SMALL)
      A4_1 = A4_1/MAX(DEN,SMALL)
      A4_0 = A4_0 + A4_1*A1_0
      DEN = 1 -A4_1*A1_4
      A4_0 = A4_0/MAX(DEN,SMALL)
      XQ(4) = A4_0
      XQ(1) = A1_0 +A1_4*XQ(4)
      XQ(6) = A6_0 +A6_4*XQ(4)
      XQ(2) = A2_0 +A2_4*XQ(4) +A2_1*XQ(1)
      XQ(3) = A3_0 +A3_4*XQ(4)
      XQ(8) = A8_0 +A8_6*XQ(6)
      XQ(5) = A5_0
      XQ(7) = A7_0
      XQ(9) = A9_0
C
      RF( 26) = RF( 26)*XQ( 4)
      RF( 27) = RF( 27)*XQ( 4)
      RF( 28) = RF( 28)*XQ( 4)
      RF( 29) = RF( 29)*XQ( 4)
      RF( 30) = RF( 30)*XQ( 4)
      RF( 31) = RF( 31)*XQ( 4)
      RF( 32) = RF( 32)*XQ( 4)
      RF( 34) = RF( 34)*XQ( 4)
      RB( 36) = RB( 36)*XQ( 4)
      RB( 38) = RB( 38)*XQ( 4)
      RB( 39) = RB( 39)*XQ( 4)
      RB( 40) = RB( 40)*XQ( 4)
      RB( 41) = RB( 41)*XQ( 4)
      RB( 42) = RB( 42)*XQ( 4)
      RB( 43) = RB( 43)*XQ( 4)
      RF( 55) = RF( 55)*XQ( 5)
      RF( 56) = RF( 56)*XQ( 5)
      RF( 57) = RF( 57)*XQ( 5)
      RF( 58) = RF( 58)*XQ( 5)
      RF( 59) = RF( 59)*XQ( 5)
      RF( 60) = RF( 60)*XQ( 5)
      RF( 61) = RF( 61)*XQ( 5)
      RF( 62) = RF( 62)*XQ( 5)
      RF( 73) = RF( 73)*XQ( 1)
      RF( 74) = RF( 74)*XQ( 2)
      RB( 75) = RB( 75)*XQ( 1)
      RB( 76) = RB( 76)*XQ( 2)
      RF( 77) = RF( 77)*XQ( 1)
      RF( 78) = RF( 78)*XQ( 2)
      RB( 79) = RB( 79)*XQ( 2)
      RF( 91) = RF( 91)*XQ( 4)
      RF( 95) = RF( 95)*XQ( 3)
      RB( 96) = RB( 96)*XQ( 3)
      RB( 97) = RB( 97)*XQ( 3)
      RB( 98) = RB( 98)*XQ( 3)
      RB( 99) = RB( 99)*XQ( 4)
      RF(100) = RF(100)*XQ( 3)
      RB(101) = RB(101)*XQ( 3)
      RB(102) = RB(102)*XQ( 3)
      RF(103) = RF(103)*XQ( 3)
      RF(104) = RF(104)*XQ( 3)
      RF(105) = RF(105)*XQ( 3)
      RF(107) = RF(107)*XQ( 3)
      RB(107) = RB(107)*XQ( 4)
      RF(108) = RF(108)*XQ( 3)
      RB(109) = RB(109)*XQ( 3)
      RB(110) = RB(110)*XQ( 1)
      RF(112) = RF(112)*XQ( 1)
      RF(113) = RF(113)*XQ( 1)
      RB(113) = RB(113)*XQ( 4)
      RF(114) = RF(114)*XQ( 1)
      RF(115) = RF(115)*XQ( 1)
      RF(116) = RF(116)*XQ( 1)
      RB(116) = RB(116)*XQ( 4)
      RF(117) = RF(117)*XQ( 1)
      RF(119) = RF(119)*XQ( 2)
      RB(119) = RB(119)*XQ( 1)
      RF(120) = RF(120)*XQ( 2)
      RB(120) = RB(120)*XQ( 1)
      RF(121) = RF(121)*XQ( 2)
      RB(121) = RB(121)*XQ( 1)
      RF(122) = RF(122)*XQ( 2)
      RB(122) = RB(122)*XQ( 1)
      RF(123) = RF(123)*XQ( 2)
      RF(124) = RF(124)*XQ( 2)
      RB(124) = RB(124)*XQ( 4)
      RF(125) = RF(125)*XQ( 2)
      RF(126) = RF(126)*XQ( 2)
      RF(127) = RF(127)*XQ( 2)
      RF(128) = RF(128)*XQ( 2)
      RF(129) = RF(129)*XQ( 2)
      RB(130) = RB(130)*XQ( 4)
      RB(132) = RB(132)*XQ( 6)
      RB(133) = RB(133)*XQ( 6)
      RB(134) = RB(134)*XQ( 6)
      RB(135) = RB(135)*XQ( 6)
      RB(136) = RB(136)*XQ( 6)
      RB(137) = RB(137)*XQ( 6)
      RF(138) = RF(138)*XQ( 6)
      RF(139) = RF(139)*XQ( 6)
      RF(140) = RF(140)*XQ( 6)
      RB(140) = RB(140)*XQ( 4)
      RF(141) = RF(141)*XQ( 6)
      RB(141) = RB(141)*XQ( 8)
      RF(142) = RF(142)*XQ( 8)
      RF(151) = RF(151)*XQ( 6)
      RF(154) = RF(154)*XQ( 8)
      RF(155) = RF(155)*XQ( 8)
      RB(156) = RB(156)*XQ( 9)
      RF(157) = RF(157)*XQ( 9)
      RF(158) = RF(158)*XQ( 9)
      RB(162) = RB(162)*XQ( 7)
      RB(163) = RB(163)*XQ( 5)
      RF(164) = RF(164)*XQ( 7)
      RB(165) = RB(165)*XQ( 7)
      RB(168) = RB(168)*XQ( 4)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RDOT(RF, RB, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RF(*), RB(*), WDOT(*)
C
      DO K = 1, 30
         WDOT(K) = 0D0
      ENDDO
C
      ROP = RF(1)-RB(1)
      WDOT(1) = WDOT(1) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(2)-RB(2)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(3)-RB(3)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(4)-RB(4)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP +ROP
      WDOT(7) = WDOT(7) -ROP
      ROP = RF(5)-RB(5)
      WDOT(1) = WDOT(1) +ROP +ROP
      WDOT(2) = WDOT(2) -ROP
      ROP = RF(6)-RB(6)
      WDOT(4) = WDOT(4) -ROP -ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(7)-RB(7)
      WDOT(1) = WDOT(1) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(8)-RB(8)
      WDOT(1) = WDOT(1) -ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(9)-RB(9)
      WDOT(1) = WDOT(1) -ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(10)-RB(10)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(15) = WDOT(15) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(11)-RB(11)
      WDOT(1) = WDOT(1) -ROP
      WDOT(6) = WDOT(6) +ROP +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(12)-RB(12)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(15) = WDOT(15) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(13)-RB(13)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(15) = WDOT(15) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(14)-RB(14)
      WDOT(15) = WDOT(15) +ROP
      WDOT(16) = WDOT(16) -ROP -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(15)-RB(15)
      WDOT(15) = WDOT(15) +ROP
      WDOT(16) = WDOT(16) -ROP -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(16)-RB(16)
      WDOT(6) = WDOT(6) +ROP +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(17)-RB(17)
      WDOT(1) = WDOT(1) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(18)-RB(18)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(19)-RB(19)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(20)-RB(20)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(21)-RB(21)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(22)-RB(22)
      WDOT(4) = WDOT(4) -ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(23)-RB(23)
      WDOT(4) = WDOT(4) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(24)-RB(24)
      WDOT(6) = WDOT(6) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(25)-RB(25)
      WDOT(1) = WDOT(1) +ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(26)-RB(26)
      WDOT(1) = WDOT(1) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(27)-RB(27)
      WDOT(9) = WDOT(9) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(28)-RB(28)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(29)-RB(29)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(30)-RB(30)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(31)-RB(31)
      WDOT(1) = WDOT(1) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(32)-RB(32)
      WDOT(1) = WDOT(1) +ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(33)-RB(33)
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) +ROP +ROP
      ROP = RF(34)-RB(34)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(35)-RB(35)
      WDOT(9) = WDOT(9) +ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(36)-RB(36)
      WDOT(1) = WDOT(1) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(37)-RB(37)
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(38)-RB(38)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(39)-RB(39)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(40)-RB(40)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(41)-RB(41)
      WDOT(12) = WDOT(12) -ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(42)-RB(42)
      WDOT(12) = WDOT(12) -ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(43)-RB(43)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(44)-RB(44)
      WDOT(1) = WDOT(1) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(45)-RB(45)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(14) = WDOT(14) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(46)-RB(46)
      WDOT(3) = WDOT(3) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(47)-RB(47)
      WDOT(3) = WDOT(3) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(14) = WDOT(14) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(48)-RB(48)
      WDOT(3) = WDOT(3) -ROP -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(49)-RB(49)
      WDOT(1) = WDOT(1) -ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      ROP = RF(50)-RB(50)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) +ROP
      WDOT(5) = WDOT(5) -ROP
      ROP = RF(51)-RB(51)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(52)-RB(52)
      WDOT(3) = WDOT(3) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(53)-RB(53)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(15) = WDOT(15) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(54)-RB(54)
      WDOT(3) = WDOT(3) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(55)-RB(55)
      WDOT(1) = WDOT(1) +ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(56)-RB(56)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(57)-RB(57)
      WDOT(1) = WDOT(1) -ROP
      WDOT(3) = WDOT(3) +ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(58)-RB(58)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(59)-RB(59)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(60)-RB(60)
      WDOT(12) = WDOT(12) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(61)-RB(61)
      WDOT(12) = WDOT(12) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(62)-RB(62)
      WDOT(12) = WDOT(12) +ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(63)-RB(63)
      WDOT(12) = WDOT(12) +ROP +ROP
      ROP = RF(64)-RB(64)
      WDOT(1) = WDOT(1) +ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(14) = WDOT(14) -ROP
      ROP = RF(65)-RB(65)
      WDOT(1) = WDOT(1) -ROP
      WDOT(3) = WDOT(3) +ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(14) = WDOT(14) -ROP
      ROP = RF(66)-RB(66)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(14) = WDOT(14) -ROP
      ROP = RF(67)-RB(67)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(14) = WDOT(14) -ROP
      ROP = RF(68)-RB(68)
      WDOT(12) = WDOT(12) +ROP
      WDOT(14) = WDOT(14) -ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(69)-RB(69)
      WDOT(12) = WDOT(12) +ROP
      WDOT(14) = WDOT(14) -ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(70)-RB(70)
      WDOT(12) = WDOT(12) +ROP
      WDOT(14) = WDOT(14) -ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(71)-RB(71)
      WDOT(3) = WDOT(3) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(14) = WDOT(14) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(72)-RB(72)
      WDOT(1) = WDOT(1) +ROP
      WDOT(3) = WDOT(3) -ROP -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(73)-RB(73)
      WDOT(3) = WDOT(3) +ROP +ROP
      WDOT(5) = WDOT(5) -ROP
      ROP = RF(74)-RB(74)
      WDOT(3) = WDOT(3) +ROP +ROP
      WDOT(5) = WDOT(5) -ROP
      ROP = RF(75)-RB(75)
      WDOT(3) = WDOT(3) -ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(76)-RB(76)
      WDOT(3) = WDOT(3) -ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(77)-RB(77)
      WDOT(1) = WDOT(1) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(78)-RB(78)
      WDOT(1) = WDOT(1) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(79)-RB(79)
      WDOT(1) = WDOT(1) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(14) = WDOT(14) -ROP
      ROP = RF(80)-RB(80)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(81)-RB(81)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(82)-RB(82)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(83)-RB(83)
      WDOT(11) = WDOT(11) +ROP
      WDOT(13) = WDOT(13) -ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(84)-RB(84)
      WDOT(11) = WDOT(11) +ROP
      WDOT(13) = WDOT(13) -ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(85)-RB(85)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(86)-RB(86)
      WDOT(1) = WDOT(1) -ROP
      WDOT(11) = WDOT(11) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(87)-RB(87)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(88)-RB(88)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(11) = WDOT(11) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(89)-RB(89)
      WDOT(10) = WDOT(10) +ROP
      WDOT(11) = WDOT(11) -ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(90)-RB(90)
      WDOT(10) = WDOT(10) +ROP
      WDOT(11) = WDOT(11) -ROP -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(91)-RB(91)
      WDOT(9) = WDOT(9) +ROP
      WDOT(11) = WDOT(11) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(92)-RB(92)
      WDOT(1) = WDOT(1) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(11) = WDOT(11) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(93)-RB(93)
      WDOT(2) = WDOT(2) +ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(94)-RB(94)
      WDOT(1) = WDOT(1) -ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(95)-RB(95)
      WDOT(1) = WDOT(1) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(96)-RB(96)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(97)-RB(97)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(98)-RB(98)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(99)-RB(99)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(100)-RB(100)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) +ROP
      ROP = RF(101)-RB(101)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(102)-RB(102)
      WDOT(10) = WDOT(10) -ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(103)-RB(103)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(8) = WDOT(8) +ROP
      ROP = RF(104)-RB(104)
      WDOT(10) = WDOT(10) +ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(105)-RB(105)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(8) = WDOT(8) +ROP
      ROP = RF(106)-RB(106)
      WDOT(8) = WDOT(8) +ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(107)-RB(107)
      WDOT(12) = WDOT(12) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(108)-RB(108)
      WDOT(8) = WDOT(8) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(109)-RB(109)
      WDOT(1) = WDOT(1) -ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(110)-RB(110)
      WDOT(4) = WDOT(4) -ROP
      WDOT(8) = WDOT(8) -ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(111)-RB(111)
      WDOT(3) = WDOT(3) +ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(8) = WDOT(8) -ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(112)-RB(112)
      WDOT(1) = WDOT(1) -ROP
      WDOT(3) = WDOT(3) +ROP
      ROP = RF(113)-RB(113)
      WDOT(1) = WDOT(1) +ROP
      WDOT(4) = WDOT(4) -ROP
      ROP = RF(114)-RB(114)
      WDOT(1) = WDOT(1) +ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(115)-RB(115)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(3) = WDOT(3) +ROP
      ROP = RF(116)-RB(116)
      WDOT(6) = WDOT(6) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(117)-RB(117)
      WDOT(6) = WDOT(6) +ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(118)-RB(118)
      WDOT(2) = WDOT(2) +ROP
      WDOT(8) = WDOT(8) +ROP
      ROP = RF(119)-RB(119)
      ROP = RF(120)-RB(120)
      ROP = RF(121)-RB(121)
      ROP = RF(122)-RB(122)
      ROP = RF(123)-RB(123)
      WDOT(2) = WDOT(2) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(124)-RB(124)
      WDOT(1) = WDOT(1) +ROP
      WDOT(4) = WDOT(4) -ROP
      ROP = RF(125)-RB(125)
      WDOT(1) = WDOT(1) +ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(126)-RB(126)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(3) = WDOT(3) +ROP
      ROP = RF(127)-RB(127)
      WDOT(1) = WDOT(1) +ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(128)-RB(128)
      WDOT(7) = WDOT(7) +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(129)-RB(129)
      WDOT(9) = WDOT(9) +ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(130)-RB(130)
      WDOT(3) = WDOT(3) +ROP
      WDOT(19) = WDOT(19) -ROP
      ROP = RF(131)-RB(131)
      WDOT(3) = WDOT(3) +ROP
      WDOT(14) = WDOT(14) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(132)-RB(132)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(133)-RB(133)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(134)-RB(134)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(135)-RB(135)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(136)-RB(136)
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(137)-RB(137)
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(138)-RB(138)
      WDOT(3) = WDOT(3) +ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(139)-RB(139)
      WDOT(12) = WDOT(12) +ROP
      WDOT(14) = WDOT(14) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(140)-RB(140)
      WDOT(12) = WDOT(12) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(141)-RB(141)
      WDOT(6) = WDOT(6) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(142)-RB(142)
      WDOT(1) = WDOT(1) +ROP
      WDOT(23) = WDOT(23) +ROP
      ROP = RF(143)-RB(143)
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(144)-RB(144)
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(145)-RB(145)
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(146)-RB(146)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(147)-RB(147)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(148)-RB(148)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(149)-RB(149)
      WDOT(9) = WDOT(9) +ROP
      WDOT(14) = WDOT(14) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(150)-RB(150)
      WDOT(3) = WDOT(3) +ROP
      WDOT(18) = WDOT(18) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(151)-RB(151)
      WDOT(15) = WDOT(15) -ROP
      WDOT(27) = WDOT(27) +ROP
      ROP = RF(152)-RB(152)
      WDOT(15) = WDOT(15) +ROP
      WDOT(27) = WDOT(27) -ROP -ROP
      ROP = RF(153)-RB(153)
      WDOT(15) = WDOT(15) +ROP
      WDOT(23) = WDOT(23) +ROP
      WDOT(24) = WDOT(24) +ROP
      WDOT(27) = WDOT(27) -ROP -ROP
      ROP = RF(154)-RB(154)
      WDOT(12) = WDOT(12) +ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(155)-RB(155)
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(23) = WDOT(23) +ROP
      ROP = RF(156)-RB(156)
      WDOT(27) = WDOT(27) -ROP
      ROP = RF(157)-RB(157)
      WDOT(6) = WDOT(6) +ROP
      WDOT(12) = WDOT(12) +ROP +ROP
      ROP = RF(158)-RB(158)
      WDOT(15) = WDOT(15) -ROP
      WDOT(29) = WDOT(29) +ROP
      ROP = RF(159)-RB(159)
      WDOT(6) = WDOT(6) +ROP
      WDOT(28) = WDOT(28) +ROP
      WDOT(29) = WDOT(29) -ROP
      ROP = RF(160)-RB(160)
      WDOT(6) = WDOT(6) +ROP
      WDOT(25) = WDOT(25) +ROP
      WDOT(28) = WDOT(28) -ROP
      ROP = RF(161)-RB(161)
      WDOT(25) = WDOT(25) -ROP
      WDOT(26) = WDOT(26) +ROP
      ROP = RF(162)-RB(162)
      WDOT(9) = WDOT(9) +ROP
      WDOT(26) = WDOT(26) -ROP
      ROP = RF(163)-RB(163)
      WDOT(18) = WDOT(18) +ROP
      WDOT(26) = WDOT(26) -ROP
      ROP = RF(164)-RB(164)
      WDOT(1) = WDOT(1) +ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(165)-RB(165)
      WDOT(6) = WDOT(6) -ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(166)-RB(166)
      WDOT(7) = WDOT(7) +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(167)-RB(167)
      WDOT(2) = WDOT(2) +ROP
      WDOT(18) = WDOT(18) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(168)-RB(168)
      WDOT(6) = WDOT(6) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(169)-RB(169)
      WDOT(1) = WDOT(1) +ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(18) = WDOT(18) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(170)-RB(170)
      WDOT(7) = WDOT(7) +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(171)-RB(171)
      WDOT(2) = WDOT(2) +ROP
      WDOT(18) = WDOT(18) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(172)-RB(172)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(173)-RB(173)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(174)-RB(174)
      WDOT(6) = WDOT(6) +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(175)-RB(175)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(20) = WDOT(20) -ROP
C
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE STIF(RF, RB, DIFF, DT, C)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RF(*), RB(*), DIFF(*), C(*)
CCCC      DATA TC/5.D8/
       TC=2.0D0/DT         !esr24Aug09
C Previously TC was a constant and this sometimes gave problems for dt less than 1ns.
C
C     H
      DDOT=+RF(1)+RB(2)+RB(3)+RB(5)+RB(5)+RF(7)+RF(8)+RF(9)+RF(10)
     *+RF(11)+RF(17)+RF(18)+RB(25)+RB(26)+RF(28)+RB(31)+RB(32)+RB(36)
     *+RF(38)+RB(44)+RF(49)+RF(50)+RB(55)+RF(56)+RF(57)+RB(64)+RF(65)
     *+RB(72)+RB(77)+RB(78)+RF(79)+RF(80)+RF(86)+RF(87)+RB(92)+RF(94)
     *+RF(95)+RF(96)+RF(103)+RF(109)+RF(112)+RB(113)+RB(114)+RB(115)
     *+RB(124)+RB(125)+RB(126)+RB(127)+RF(133)+RB(142)+RF(147)+RB(164)
     *+RB(169)+RF(171)+RB(171)+RF(172)
      TINV=DDOT/C(1)
      IF (TINV .GT. TC) THEN
      CDOT=+RB(1)+RF(2)+RF(3)+RF(5)+RF(5)+RB(7)+RB(8)+RB(9)+RB(10)
     *+RB(11)+RB(17)+RB(18)+RF(25)+RF(26)+RB(28)+RF(31)+RF(32)+RF(36)
     *+RB(38)+RF(44)+RB(49)+RB(50)+RF(55)+RB(56)+RB(57)+RF(64)+RB(65)
     *+RF(72)+RF(77)+RF(78)+RB(79)+RB(80)+RB(86)+RB(87)+RF(92)+RB(94)
     *+RB(95)+RB(96)+RB(103)+RB(109)+RB(112)+RF(113)+RF(114)+RF(115)
     *+RF(124)+RF(125)+RF(126)+RF(127)+RB(133)+RF(142)+RB(147)+RF(164)
     *+RF(169)+RB(171)+RF(171)+RB(172)
      C0=C(1)*(CDOT+DIFF(1)/1.00796998)/DDOT
      C0=C(1)*(CDOT+DIFF(1)/1.00796998+(C(1)-C0)/DT)/DDOT
      R=C0/C(1)
      RF(1)=RF(1)*R
      RB(2)=RB(2)*R
      RB(3)=RB(3)*R
      RB(5)=RB(5)*R
      RB(5)=RB(5)*R
      RF(7)=RF(7)*R
      RF(8)=RF(8)*R
      RF(9)=RF(9)*R
      RF(10)=RF(10)*R
      RF(11)=RF(11)*R
      RF(17)=RF(17)*R
      RF(18)=RF(18)*R
      RB(25)=RB(25)*R
      RB(26)=RB(26)*R
      RF(28)=RF(28)*R
      RB(31)=RB(31)*R
      RB(32)=RB(32)*R
      RB(36)=RB(36)*R
      RF(38)=RF(38)*R
      RB(44)=RB(44)*R
      RF(49)=RF(49)*R
      RF(50)=RF(50)*R
      RB(55)=RB(55)*R
      RF(56)=RF(56)*R
      RF(57)=RF(57)*R
      RB(64)=RB(64)*R
      RF(65)=RF(65)*R
      RB(72)=RB(72)*R
      RB(77)=RB(77)*R
      RB(78)=RB(78)*R
      RF(79)=RF(79)*R
      RF(80)=RF(80)*R
      RF(86)=RF(86)*R
      RF(87)=RF(87)*R
      RB(92)=RB(92)*R
      RF(94)=RF(94)*R
      RF(95)=RF(95)*R
      RF(96)=RF(96)*R
      RF(103)=RF(103)*R
      RF(109)=RF(109)*R
      RF(112)=RF(112)*R
      RB(113)=RB(113)*R
      RB(114)=RB(114)*R
      RB(115)=RB(115)*R
      RB(124)=RB(124)*R
      RB(125)=RB(125)*R
      RB(126)=RB(126)*R
      RB(127)=RB(127)*R
      RF(133)=RF(133)*R
      RB(142)=RB(142)*R
      RF(147)=RF(147)*R
      RB(164)=RB(164)*R
      RB(169)=RB(169)*R
      RF(171)=RF(171)*R
      RB(171)=RB(171)*R
      RF(172)=RF(172)*R
      ENDIF
C
C     H2
      DDOT=+RF(2)+RF(3)+RF(5)+RB(10)+RB(18)+RB(28)+RB(33)+RB(37)+RB(38)
     *+RB(50)+RB(56)+RB(80)+RB(87)+RB(93)+RB(96)+RB(103)+RF(115)+RB(118)
     *+RB(123)+RF(126)+RB(133)+RB(147)+RB(167)+RB(171)+RB(172)
      TINV=DDOT/C(2)
      IF (TINV .GT. TC) THEN
      CDOT=+RB(2)+RB(3)+RB(5)+RF(10)+RF(18)+RF(28)+RF(33)+RF(37)+RF(38)
     *+RF(50)+RF(56)+RF(80)+RF(87)+RF(93)+RF(96)+RF(103)+RB(115)+RF(118)
     *+RF(123)+RB(126)+RF(133)+RF(147)+RF(167)+RF(171)+RF(172)
      C0=C(2)*(CDOT+DIFF(2)/2.01593995)/DDOT
      C0=C(2)*(CDOT+DIFF(2)/2.01593995+(C(2)-C0)/DT)/DDOT
      R=C0/C(2)
      RF(2)=RF(2)*R
      RF(3)=RF(3)*R
      RF(5)=RF(5)*R
      RB(10)=RB(10)*R
      RB(18)=RB(18)*R
      RB(28)=RB(28)*R
      RB(33)=RB(33)*R
      RB(37)=RB(37)*R
      RB(38)=RB(38)*R
      RB(50)=RB(50)*R
      RB(56)=RB(56)*R
      RB(80)=RB(80)*R
      RB(87)=RB(87)*R
      RB(93)=RB(93)*R
      RB(96)=RB(96)*R
      RB(103)=RB(103)*R
      RF(115)=RF(115)*R
      RB(118)=RB(118)*R
      RB(123)=RB(123)*R
      RF(126)=RF(126)*R
      RB(133)=RB(133)*R
      RB(147)=RB(147)*R
      RB(167)=RB(167)*R
      RB(171)=RB(171)*R
      RB(172)=RB(172)*R
      ENDIF
C
C     CH3
      DDOT=+RF(34)+RF(43)+RF(44)+RF(45)+RF(46)+RF(47)+RF(48)+RF(48)
     *+RF(49)+RB(50)+RB(51)+RB(52)+RF(53)+RB(54)+RB(57)+RB(65)+RB(71)
     *+RF(72)+RF(72)+RB(73)+RB(73)+RB(74)+RB(74)+RF(75)+RF(76)+RF(77)
     *+RF(78)+RF(85)+RB(88)+RF(98)+RB(99)+RF(105)+RB(111)+RB(112)
     *+RB(115)+RB(126)+RB(130)+RB(131)+RF(134)+RB(138)+RF(148)+RB(150)
     *+RF(173)
      TINV=DDOT/C(3)
      IF (TINV .GT. TC) THEN
      CDOT=+RB(34)+RB(43)+RB(44)+RB(45)+RB(46)+RB(47)+RB(48)+RB(48)
     *+RB(49)+RF(50)+RF(51)+RF(52)+RB(53)+RF(54)+RF(57)+RF(65)+RF(71)
     *+RB(72)+RB(72)+RF(73)+RF(73)+RF(74)+RF(74)+RB(75)+RB(76)+RB(77)
     *+RB(78)+RB(85)+RF(88)+RB(98)+RF(99)+RB(105)+RF(111)+RF(112)
     *+RF(115)+RF(126)+RF(130)+RF(131)+RB(134)+RF(138)+RB(148)+RF(150)
     *+RB(173)
      C0=C(3)*(CDOT+DIFF(3)/15.03506029)/DDOT
      C0=C(3)*(CDOT+DIFF(3)/15.03506029+(C(3)-C0)/DT)/DDOT
      R=C0/C(3)
      RF(34)=RF(34)*R
      RF(43)=RF(43)*R
      RF(44)=RF(44)*R
      RF(45)=RF(45)*R
      RF(46)=RF(46)*R
      RF(47)=RF(47)*R
      RF(48)=RF(48)*R
      RF(48)=RF(48)*R
      RF(49)=RF(49)*R
      RB(50)=RB(50)*R
      RB(51)=RB(51)*R
      RB(52)=RB(52)*R
      RF(53)=RF(53)*R
      RB(54)=RB(54)*R
      RB(57)=RB(57)*R
      RB(65)=RB(65)*R
      RB(71)=RB(71)*R
      RF(72)=RF(72)*R
      RF(72)=RF(72)*R
      RB(73)=RB(73)*R
      RB(73)=RB(73)*R
      RB(74)=RB(74)*R
      RB(74)=RB(74)*R
      RF(75)=RF(75)*R
      RF(76)=RF(76)*R
      RF(77)=RF(77)*R
      RF(78)=RF(78)*R
      RF(85)=RF(85)*R
      RB(88)=RB(88)*R
      RF(98)=RF(98)*R
      RB(99)=RB(99)*R
      RF(105)=RF(105)*R
      RB(111)=RB(111)*R
      RB(112)=RB(112)*R
      RB(115)=RB(115)*R
      RB(126)=RB(126)*R
      RB(130)=RB(130)*R
      RB(131)=RB(131)*R
      RF(134)=RF(134)*R
      RB(138)=RB(138)*R
      RF(148)=RF(148)*R
      RB(150)=RB(150)*R
      RF(173)=RF(173)*R
      ENDIF
C
C     O
      DDOT=+RB(1)+RF(2)+RF(4)+RF(6)+RF(6)+RF(7)+RF(12)+RF(19)+RF(22)
     *+RB(23)+RF(29)+RF(31)+RF(39)+RF(44)+RB(45)+RF(51)+RF(58)+RF(66)
     *+RF(81)+RF(88)+RF(92)+RF(99)+RF(101)+RF(110)+RF(113)+RF(123)
     *+RF(124)+RF(135)+RF(146)+RF(175)
      TINV=DDOT/C(4)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(1)+RB(2)+RB(4)+RB(6)+RB(6)+RB(7)+RB(12)+RB(19)+RB(22)
     *+RF(23)+RB(29)+RB(31)+RB(39)+RB(44)+RF(45)+RB(51)+RB(58)+RB(66)
     *+RB(81)+RB(88)+RB(92)+RB(99)+RB(101)+RB(110)+RB(113)+RB(123)
     *+RB(124)+RB(135)+RB(146)+RB(175)
      C0=C(4)*(CDOT+DIFF(4)/15.99940014)/DDOT
      C0=C(4)*(CDOT+DIFF(4)/15.99940014+(C(4)-C0)/DT)/DDOT
      R=C0/C(4)
      RB(1)=RB(1)*R
      RF(2)=RF(2)*R
      RF(4)=RF(4)*R
      RF(6)=RF(6)*R
      RF(6)=RF(6)*R
      RF(7)=RF(7)*R
      RF(12)=RF(12)*R
      RF(19)=RF(19)*R
      RF(22)=RF(22)*R
      RB(23)=RB(23)*R
      RF(29)=RF(29)*R
      RF(31)=RF(31)*R
      RF(39)=RF(39)*R
      RF(44)=RF(44)*R
      RB(45)=RB(45)*R
      RF(51)=RF(51)*R
      RF(58)=RF(58)*R
      RF(66)=RF(66)*R
      RF(81)=RF(81)*R
      RF(88)=RF(88)*R
      RF(92)=RF(92)*R
      RF(99)=RF(99)*R
      RF(101)=RF(101)*R
      RF(110)=RF(110)*R
      RF(113)=RF(113)*R
      RF(123)=RF(123)*R
      RF(124)=RF(124)*R
      RF(135)=RF(135)*R
      RF(146)=RF(146)*R
      RF(175)=RF(175)*R
      ENDIF
C
C     CH4
      DDOT=+RB(34)+RB(43)+RB(49)+RF(50)+RF(51)+RF(52)+RB(53)+RF(54)
     *+RF(73)+RF(74)+RB(85)+RB(98)+RB(105)+RB(134)+RB(148)+RB(173)
      TINV=DDOT/C(5)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(34)+RF(43)+RF(49)+RB(50)+RB(51)+RB(52)+RF(53)+RB(54)
     *+RB(73)+RB(74)+RF(85)+RF(98)+RF(105)+RF(134)+RF(148)+RF(173)
      C0=C(5)*(CDOT+DIFF(5)/16.04303026)/DDOT
      C0=C(5)*(CDOT+DIFF(5)/16.04303026+(C(5)-C0)/DT)/DDOT
      R=C0/C(5)
      RB(34)=RB(34)*R
      RB(43)=RB(43)*R
      RB(49)=RB(49)*R
      RF(50)=RF(50)*R
      RF(51)=RF(51)*R
      RF(52)=RF(52)*R
      RB(53)=RB(53)*R
      RF(54)=RF(54)*R
      RF(73)=RF(73)*R
      RF(74)=RF(74)*R
      RB(85)=RB(85)*R
      RB(98)=RB(98)*R
      RB(105)=RB(105)*R
      RB(134)=RB(134)*R
      RB(148)=RB(148)*R
      RB(173)=RB(173)*R
      ENDIF
C
C     OH
      DDOT=+RB(1)+RB(2)+RF(3)+RB(4)+RB(4)+RB(7)+RF(8)+RB(11)+RB(11)
     *+RB(12)+RF(13)+RB(16)+RB(16)+RB(17)+RB(19)+RF(20)+RF(21)+RB(24)
     *+RF(25)+RB(29)+RF(30)+RB(32)+RB(39)+RF(40)+RB(46)+RB(47)+RB(51)
     *+RF(52)+RB(57)+RB(58)+RF(59)+RB(65)+RB(66)+RF(67)+RF(75)+RF(76)
     *+RB(81)+RF(82)+RF(97)+RF(100)+RB(101)+RF(111)+RF(114)+RB(116)
     *+RB(117)+RF(125)+RB(127)+RF(132)+RB(135)+RB(141)+RF(144)+RB(146)
     *+RB(157)+RB(159)+RB(160)+RF(165)+RB(168)+RF(169)+RF(170)+RB(170)
     *+RB(172)+RB(173)+RB(174)+RB(175)+RB(175)
      TINV=DDOT/C(6)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(1)+RF(2)+RB(3)+RF(4)+RF(4)+RF(7)+RB(8)+RF(11)+RF(11)
     *+RF(12)+RB(13)+RF(16)+RF(16)+RF(17)+RF(19)+RB(20)+RB(21)+RF(24)
     *+RB(25)+RF(29)+RB(30)+RF(32)+RF(39)+RB(40)+RF(46)+RF(47)+RF(51)
     *+RB(52)+RF(57)+RF(58)+RB(59)+RF(65)+RF(66)+RB(67)+RB(75)+RB(76)
     *+RF(81)+RB(82)+RB(97)+RB(100)+RF(101)+RB(111)+RB(114)+RF(116)
     *+RF(117)+RB(125)+RF(127)+RB(132)+RF(135)+RF(141)+RB(144)+RF(146)
     *+RF(157)+RF(159)+RF(160)+RB(165)+RF(168)+RB(169)+RB(170)+RF(170)
     *+RF(172)+RF(173)+RF(174)+RF(175)+RF(175)
      C0=C(6)*(CDOT+DIFF(6)/17.00737011)/DDOT
      C0=C(6)*(CDOT+DIFF(6)/17.00737011+(C(6)-C0)/DT)/DDOT
      R=C0/C(6)
      RB(1)=RB(1)*R
      RB(2)=RB(2)*R
      RF(3)=RF(3)*R
      RB(4)=RB(4)*R
      RB(4)=RB(4)*R
      RB(7)=RB(7)*R
      RF(8)=RF(8)*R
      RB(11)=RB(11)*R
      RB(11)=RB(11)*R
      RB(12)=RB(12)*R
      RF(13)=RF(13)*R
      RB(16)=RB(16)*R
      RB(16)=RB(16)*R
      RB(17)=RB(17)*R
      RB(19)=RB(19)*R
      RF(20)=RF(20)*R
      RF(21)=RF(21)*R
      RB(24)=RB(24)*R
      RF(25)=RF(25)*R
      RB(29)=RB(29)*R
      RF(30)=RF(30)*R
      RB(32)=RB(32)*R
      RB(39)=RB(39)*R
      RF(40)=RF(40)*R
      RB(46)=RB(46)*R
      RB(47)=RB(47)*R
      RB(51)=RB(51)*R
      RF(52)=RF(52)*R
      RB(57)=RB(57)*R
      RB(58)=RB(58)*R
      RF(59)=RF(59)*R
      RB(65)=RB(65)*R
      RB(66)=RB(66)*R
      RF(67)=RF(67)*R
      RF(75)=RF(75)*R
      RF(76)=RF(76)*R
      RB(81)=RB(81)*R
      RF(82)=RF(82)*R
      RF(97)=RF(97)*R
      RF(100)=RF(100)*R
      RB(101)=RB(101)*R
      RF(111)=RF(111)*R
      RF(114)=RF(114)*R
      RB(116)=RB(116)*R
      RB(117)=RB(117)*R
      RF(125)=RF(125)*R
      RB(127)=RB(127)*R
      RF(132)=RF(132)*R
      RB(135)=RB(135)*R
      RB(141)=RB(141)*R
      RF(144)=RF(144)*R
      RB(146)=RB(146)*R
      RB(157)=RB(157)*R
      RB(159)=RB(159)*R
      RB(160)=RB(160)*R
      RF(165)=RF(165)*R
      RB(168)=RB(168)*R
      RF(169)=RF(169)*R
      RF(170)=RF(170)*R
      RB(170)=RB(170)*R
      RB(172)=RB(172)*R
      RB(173)=RB(173)*R
      RB(174)=RB(174)*R
      RB(175)=RB(175)*R
      RB(175)=RB(175)*R
      ENDIF
C
C     C2H4
      DDOT=+RB(77)+RB(78)+RB(87)+RB(89)+RB(90)+RF(93)+RF(94)+RB(95)
     *+RF(96)+RF(97)+RF(98)+RF(99)+RF(101)+RF(102)+RB(104)+RB(106)
      TINV=DDOT/C(10)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(77)+RF(78)+RF(87)+RF(89)+RF(90)+RB(93)+RB(94)+RF(95)
     *+RB(96)+RB(97)+RB(98)+RB(99)+RB(101)+RB(102)+RF(104)+RF(106)
      C0=C(10)*(CDOT+DIFF(10)/28.05418062)/DDOT
      C0=C(10)*(CDOT+DIFF(10)/28.05418062+(C(10)-C0)/DT)/DDOT
      R=C0/C(10)
      RB(77)=RB(77)*R
      RB(78)=RB(78)*R
      RB(87)=RB(87)*R
      RB(89)=RB(89)*R
      RB(90)=RB(90)*R
      RF(93)=RF(93)*R
      RF(94)=RF(94)*R
      RB(95)=RB(95)*R
      RF(96)=RF(96)*R
      RF(97)=RF(97)*R
      RF(98)=RF(98)*R
      RF(99)=RF(99)*R
      RF(101)=RF(101)*R
      RF(102)=RF(102)*R
      RB(104)=RB(104)*R
      RB(106)=RB(106)*R
      ENDIF
C
C     C2H5
      DDOT=+RB(72)+RB(80)+RB(81)+RB(82)+RB(83)+RB(84)+RB(85)+RF(86)
     *+RF(87)+RF(88)+RF(89)+RF(90)+RF(90)+RF(91)+RF(92)+RB(94)
      TINV=DDOT/C(11)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(72)+RF(80)+RF(81)+RF(82)+RF(83)+RF(84)+RF(85)+RB(86)
     *+RB(87)+RB(88)+RB(89)+RB(90)+RB(90)+RB(91)+RB(92)+RF(94)
      C0=C(11)*(CDOT+DIFF(11)/29.06215060)/DDOT
      C0=C(11)*(CDOT+DIFF(11)/29.06215060+(C(11)-C0)/DT)/DDOT
      R=C0/C(11)
      RB(72)=RB(72)*R
      RB(80)=RB(80)*R
      RB(81)=RB(81)*R
      RB(82)=RB(82)*R
      RB(83)=RB(83)*R
      RB(84)=RB(84)*R
      RB(85)=RB(85)*R
      RF(86)=RF(86)*R
      RF(87)=RF(87)*R
      RF(88)=RF(88)*R
      RF(89)=RF(89)*R
      RF(90)=RF(90)*R
      RF(90)=RF(90)*R
      RF(91)=RF(91)*R
      RF(92)=RF(92)*R
      RB(94)=RB(94)*R
      ENDIF
C
C     CH2O
      DDOT=+RB(35)+RF(36)+RF(37)+RF(38)+RF(39)+RF(40)+RF(41)+RF(42)
     *+RF(43)+RB(44)+RB(46)+RB(55)+RB(56)+RB(58)+RB(59)+RB(60)+RB(61)
     *+RB(62)+RB(63)+RB(63)+RB(64)+RB(66)+RB(67)+RB(68)+RB(69)+RB(70)
     *+RB(88)+RB(107)+RB(114)+RB(117)+RB(125)+RB(129)+RB(138)+RB(139)
     *+RF(140)+RB(154)+RB(157)+RB(157)+RF(165)
      TINV=DDOT/C(12)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(35)+RB(36)+RB(37)+RB(38)+RB(39)+RB(40)+RB(41)+RB(42)
     *+RB(43)+RF(44)+RF(46)+RF(55)+RF(56)+RF(58)+RF(59)+RF(60)+RF(61)
     *+RF(62)+RF(63)+RF(63)+RF(64)+RF(66)+RF(67)+RF(68)+RF(69)+RF(70)
     *+RF(88)+RF(107)+RF(114)+RF(117)+RF(125)+RF(129)+RF(138)+RF(139)
     *+RB(140)+RF(154)+RF(157)+RF(157)+RB(165)
      C0=C(12)*(CDOT+DIFF(12)/30.02649045)/DDOT
      C0=C(12)*(CDOT+DIFF(12)/30.02649045+(C(12)-C0)/DT)/DDOT
      R=C0/C(12)
      RB(35)=RB(35)*R
      RF(36)=RF(36)*R
      RF(37)=RF(37)*R
      RF(38)=RF(38)*R
      RF(39)=RF(39)*R
      RF(40)=RF(40)*R
      RF(41)=RF(41)*R
      RF(42)=RF(42)*R
      RF(43)=RF(43)*R
      RB(44)=RB(44)*R
      RB(46)=RB(46)*R
      RB(55)=RB(55)*R
      RB(56)=RB(56)*R
      RB(58)=RB(58)*R
      RB(59)=RB(59)*R
      RB(60)=RB(60)*R
      RB(61)=RB(61)*R
      RB(62)=RB(62)*R
      RB(63)=RB(63)*R
      RB(63)=RB(63)*R
      RB(64)=RB(64)*R
      RB(66)=RB(66)*R
      RB(67)=RB(67)*R
      RB(68)=RB(68)*R
      RB(69)=RB(69)*R
      RB(70)=RB(70)*R
      RB(88)=RB(88)*R
      RB(107)=RB(107)*R
      RB(114)=RB(114)*R
      RB(117)=RB(117)*R
      RB(125)=RB(125)*R
      RB(129)=RB(129)*R
      RB(138)=RB(138)*R
      RB(139)=RB(139)*R
      RF(140)=RF(140)*R
      RB(154)=RB(154)*R
      RB(157)=RB(157)*R
      RB(157)=RB(157)*R
      RF(165)=RF(165)*R
      ENDIF
C
C     C2H6
      DDOT=+RB(48)+RF(80)+RF(81)+RF(82)+RF(83)+RF(84)+RF(85)+RB(86)
     *+RB(90)+RB(91)
      TINV=DDOT/C(13)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(48)+RB(80)+RB(81)+RB(82)+RB(83)+RB(84)+RB(85)+RF(86)
     *+RF(90)+RF(91)
      C0=C(13)*(CDOT+DIFF(13)/30.07012057)/DDOT
      C0=C(13)*(CDOT+DIFF(13)/30.07012057+(C(13)-C0)/DT)/DDOT
      R=C0/C(13)
      RB(48)=RB(48)*R
      RF(80)=RF(80)*R
      RF(81)=RF(81)*R
      RF(82)=RF(82)*R
      RF(83)=RF(83)*R
      RF(84)=RF(84)*R
      RF(85)=RF(85)*R
      RB(86)=RB(86)*R
      RB(90)=RB(90)*R
      RB(91)=RB(91)*R
      ENDIF
C
C     CH3O
      DDOT=+RB(45)+RB(47)+RF(64)+RF(65)+RF(66)+RF(67)+RF(68)+RF(69)
     *+RF(70)+RF(71)+RF(79)+RB(131)+RF(139)+RB(149)+RB(154)
      TINV=DDOT/C(14)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(45)+RF(47)+RB(64)+RB(65)+RB(66)+RB(67)+RB(68)+RB(69)
     *+RB(70)+RB(71)+RB(79)+RF(131)+RB(139)+RF(149)+RF(154)
      C0=C(14)*(CDOT+DIFF(14)/31.03446043)/DDOT
      C0=C(14)*(CDOT+DIFF(14)/31.03446043+(C(14)-C0)/DT)/DDOT
      R=C0/C(14)
      RB(45)=RB(45)*R
      RB(47)=RB(47)*R
      RF(64)=RF(64)*R
      RF(65)=RF(65)*R
      RF(66)=RF(66)*R
      RF(67)=RF(67)*R
      RF(68)=RF(68)*R
      RF(69)=RF(69)*R
      RF(70)=RF(70)*R
      RF(71)=RF(71)*R
      RF(79)=RF(79)*R
      RB(131)=RB(131)*R
      RF(139)=RF(139)*R
      RB(149)=RB(149)*R
      RB(154)=RB(154)*R
      ENDIF
C
C     HO2
      DDOT=+RB(9)+RF(10)+RF(11)+RF(12)+RF(13)+RF(14)+RF(14)+RF(15)
     *+RF(15)+RB(18)+RB(19)+RB(20)+RB(21)+RF(24)+RB(27)+RF(32)+RB(41)
     *+RF(42)+RF(47)+RF(53)+RF(54)+RB(60)+RB(61)+RF(62)+RB(68)+RB(69)
     *+RF(70)+RB(83)+RF(84)+RB(89)+RB(102)+RB(104)+RB(108)+RF(117)
     *+RF(136)+RB(137)+RF(141)+RB(143)+RF(145)+RB(155)+RF(174)
      TINV=DDOT/C(16)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(9)+RB(10)+RB(11)+RB(12)+RB(13)+RB(14)+RB(14)+RB(15)
     *+RB(15)+RF(18)+RF(19)+RF(20)+RF(21)+RB(24)+RF(27)+RB(32)+RF(41)
     *+RB(42)+RB(47)+RB(53)+RB(54)+RF(60)+RF(61)+RB(62)+RF(68)+RF(69)
     *+RB(70)+RF(83)+RB(84)+RF(89)+RF(102)+RF(104)+RF(108)+RB(117)
     *+RB(136)+RF(137)+RB(141)+RF(143)+RB(145)+RF(155)+RB(174)
      C0=C(16)*(CDOT+DIFF(16)/33.00677025)/DDOT
      C0=C(16)*(CDOT+DIFF(16)/33.00677025+(C(16)-C0)/DT)/DDOT
      R=C0/C(16)
      RB(9)=RB(9)*R
      RF(10)=RF(10)*R
      RF(11)=RF(11)*R
      RF(12)=RF(12)*R
      RF(13)=RF(13)*R
      RF(14)=RF(14)*R
      RF(14)=RF(14)*R
      RF(15)=RF(15)*R
      RF(15)=RF(15)*R
      RB(18)=RB(18)*R
      RB(19)=RB(19)*R
      RB(20)=RB(20)*R
      RB(21)=RB(21)*R
      RF(24)=RF(24)*R
      RB(27)=RB(27)*R
      RF(32)=RF(32)*R
      RB(41)=RB(41)*R
      RF(42)=RF(42)*R
      RF(47)=RF(47)*R
      RF(53)=RF(53)*R
      RF(54)=RF(54)*R
      RB(60)=RB(60)*R
      RB(61)=RB(61)*R
      RF(62)=RF(62)*R
      RB(68)=RB(68)*R
      RB(69)=RB(69)*R
      RF(70)=RF(70)*R
      RB(83)=RB(83)*R
      RF(84)=RF(84)*R
      RB(89)=RB(89)*R
      RB(102)=RB(102)*R
      RB(104)=RB(104)*R
      RB(108)=RB(108)*R
      RF(117)=RF(117)*R
      RF(136)=RF(136)*R
      RB(137)=RB(137)*R
      RF(141)=RF(141)*R
      RB(143)=RB(143)*R
      RF(145)=RF(145)*R
      RB(155)=RB(155)*R
      RF(174)=RF(174)*R
      ENDIF
C
C     H2O2
      DDOT=+RB(14)+RB(15)+RF(16)+RF(17)+RF(18)+RF(19)+RF(20)+RF(21)
     *+RB(42)+RB(54)+RB(62)+RB(70)+RB(84)+RF(104)+RB(136)+RB(145)
     *+RB(174)
      TINV=DDOT/C(17)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(14)+RF(15)+RB(16)+RB(17)+RB(18)+RB(19)+RB(20)+RB(21)
     *+RF(42)+RF(54)+RF(62)+RF(70)+RF(84)+RB(104)+RF(136)+RF(145)
     *+RF(174)
      C0=C(17)*(CDOT+DIFF(17)/34.01474023)/DDOT
      C0=C(17)*(CDOT+DIFF(17)/34.01474023+(C(17)-C0)/DT)/DDOT
      R=C0/C(17)
      RB(14)=RB(14)*R
      RB(15)=RB(15)*R
      RF(16)=RF(16)*R
      RF(17)=RF(17)*R
      RF(18)=RF(18)*R
      RF(19)=RF(19)*R
      RF(20)=RF(20)*R
      RF(21)=RF(21)*R
      RB(42)=RB(42)*R
      RB(54)=RB(54)*R
      RB(62)=RB(62)*R
      RB(70)=RB(70)*R
      RB(84)=RB(84)*R
      RF(104)=RF(104)*R
      RB(136)=RB(136)*R
      RB(145)=RB(145)*R
      RB(174)=RB(174)*R
      ENDIF
C
C     CH3HCO
      DDOT=+RB(92)+RF(130)
      TINV=DDOT/C(19)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(92)+RB(130)
      C0=C(19)*(CDOT+DIFF(19)/44.05358076)/DDOT
      C0=C(19)*(CDOT+DIFF(19)/44.05358076+(C(19)-C0)/DT)/DDOT
      R=C0/C(19)
      RB(92)=RB(92)*R
      RF(130)=RF(130)*R
      ENDIF
C
C     HCOOH
      DDOT=+RB(164)+RF(166)+RF(167)+RF(168)+RF(169)+RF(170)+RF(171)
     *+RF(172)+RF(173)+RF(174)+RF(175)
      TINV=DDOT/C(20)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(164)+RB(166)+RB(167)+RB(168)+RB(169)+RB(170)+RB(171)
     *+RB(172)+RB(173)+RB(174)+RB(175)
      C0=C(20)*(CDOT+DIFF(20)/46.02589059)/DDOT
      C0=C(20)*(CDOT+DIFF(20)/46.02589059+(C(20)-C0)/DT)/DDOT
      R=C0/C(20)
      RB(164)=RB(164)*R
      RF(166)=RF(166)*R
      RF(167)=RF(167)*R
      RF(168)=RF(168)*R
      RF(169)=RF(169)*R
      RF(170)=RF(170)*R
      RF(171)=RF(171)*R
      RF(172)=RF(172)*R
      RF(173)=RF(173)*R
      RF(174)=RF(174)*R
      RF(175)=RF(175)*R
      ENDIF
C
C     CH3OCH3
      DDOT=+RF(131)+RF(132)+RF(133)+RF(134)+RF(135)+RF(136)+RF(137)
     *+RB(139)+RB(140)
      TINV=DDOT/C(21)
      IF (TINV .GT. TC) THEN
      CDOT=+RB(131)+RB(132)+RB(133)+RB(134)+RB(135)+RB(136)+RB(137)
     *+RF(139)+RF(140)
      C0=C(21)*(CDOT+DIFF(21)/46.06952071)/DDOT
      C0=C(21)*(CDOT+DIFF(21)/46.06952071+(C(21)-C0)/DT)/DDOT
      R=C0/C(21)
      RF(131)=RF(131)*R
      RF(132)=RF(132)*R
      RF(133)=RF(133)*R
      RF(134)=RF(134)*R
      RF(135)=RF(135)*R
      RF(136)=RF(136)*R
      RF(137)=RF(137)*R
      RB(139)=RB(139)*R
      RB(140)=RB(140)*R
      ENDIF
C
C     CH3OCHO
      DDOT=+RB(142)+RF(143)+RF(144)+RF(145)+RF(146)+RF(147)+RF(148)
     *+RB(153)+RB(155)
      TINV=DDOT/C(23)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(142)+RB(143)+RB(144)+RB(145)+RB(146)+RB(147)+RB(148)
     *+RF(153)+RF(155)
      C0=C(23)*(CDOT+DIFF(23)/60.05298090)/DDOT
      C0=C(23)*(CDOT+DIFF(23)/60.05298090+(C(23)-C0)/DT)/DDOT
      R=C0/C(23)
      RB(142)=RB(142)*R
      RF(143)=RF(143)*R
      RF(144)=RF(144)*R
      RF(145)=RF(145)*R
      RF(146)=RF(146)*R
      RF(147)=RF(147)*R
      RF(148)=RF(148)*R
      RB(153)=RB(153)*R
      RB(155)=RB(155)*R
      ENDIF
C
C     OCH2OCHO
      DDOT=+RB(160)+RF(161)
      TINV=DDOT/C(25)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(160)+RB(161)
      C0=C(25)*(CDOT+DIFF(25)/75.04441106)/DDOT
      C0=C(25)*(CDOT+DIFF(25)/75.04441106+(C(25)-C0)/DT)/DDOT
      R=C0/C(25)
      RB(160)=RB(160)*R
      RF(161)=RF(161)*R
      ENDIF
C
C     HOCH2OCO
      DDOT=+RB(161)+RF(162)+RF(163)
      TINV=DDOT/C(26)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(161)+RB(162)+RB(163)
      C0=C(26)*(CDOT+DIFF(26)/75.04441106)/DDOT
      C0=C(26)*(CDOT+DIFF(26)/75.04441106+(C(26)-C0)/DT)/DDOT
      R=C0/C(26)
      RB(161)=RB(161)*R
      RF(162)=RF(162)*R
      RF(163)=RF(163)*R
      ENDIF
C
C     CH3OCH2O2
      DDOT=+RB(151)+RF(152)+RF(152)+RF(153)+RF(153)+RF(156)
      TINV=DDOT/C(27)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(151)+RB(152)+RB(152)+RB(153)+RB(153)+RB(156)
      C0=C(27)*(CDOT+DIFF(27)/77.06035101)/DDOT
      C0=C(27)*(CDOT+DIFF(27)/77.06035101+(C(27)-C0)/DT)/DDOT
      R=C0/C(27)
      RB(151)=RB(151)*R
      RF(152)=RF(152)*R
      RF(152)=RF(152)*R
      RF(153)=RF(153)*R
      RF(153)=RF(153)*R
      RF(156)=RF(156)*R
      ENDIF
C
C     HO2CH2OCHO
      DDOT=+RB(159)+RF(160)
      TINV=DDOT/C(28)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(159)+RB(160)
      C0=C(28)*(CDOT+DIFF(28)/92.05178118)/DDOT
      C0=C(28)*(CDOT+DIFF(28)/92.05178118+(C(28)-C0)/DT)/DDOT
      R=C0/C(28)
      RB(159)=RB(159)*R
      RF(160)=RF(160)*R
      ENDIF
C
C     O2CH2OCH2O2H
      DDOT=+RB(158)+RF(159)
      TINV=DDOT/C(29)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(158)+RB(159)
      C0=C(29)*(CDOT+DIFF(29)/109.05915129)/DDOT
      C0=C(29)*(CDOT+DIFF(29)/109.05915129+(C(29)-C0)/DT)/DDOT
      R=C0/C(29)
      RB(158)=RB(158)*R
      RF(159)=RF(159)*R
      ENDIF
      END
