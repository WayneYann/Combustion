      SUBROUTINE EGINICD (NP, LOUT, IFLAG, ITLS,
     &                    WEG, LWEG, IWEG, LIWEG)
C-----------------------------------------------------------------------
C
C     This subroutine initializes the pointers for the work arrays
C     WEG and IWEG and checks their length.
C     This subroutine should be called by the user once at the
C     beginning of the program.
C
C     Input
C     -----
C        NP        number of nodes
C        LOUT      output file number
C        IFLAG     flag for evaluating parameters and space allocation
C                  (see below)
C        ITLS      flag for space allocation (see below)
C        WEG       double precision work array for EGLIB
C        LWEG      length of WEG declared in main code
C        IWEG      integer work array for EGLIB
C        LIWEG     length of IWEG declared in main code
C
C
C     The value of IFLAG and ITLS depends on the subroutines that
C     will be used as indicated by the following table
C
C
C     Subroutine     ITLS      IFLAG
C
C     EG*D(R)1         1         2
C     EG*D(R)2         1         2
C
C     EG*E1            0         1
C     EG*E2            1         2
C     EG*E3            1         3
C     EG*E4            1         3
C
C     EG*K1            0         4
C     EG*K2            1         4
C     EG*K3            1         5
C     EG*K4            2         4
C     EG*K5            2         5
C     EG*K6            2         5
C
C     EG*L1            0         1
C     EG*L2            1         6
C     EG*L3            1         7
C     EG*L4            2         6
C     EG*L5            2         7
C
C     EG*LC1           1         7
C     EG*LC2           1         7
C     EG*LC3           2         7
C     EG*LC4           2         7
C
C     EG*LTD(R)1       2         7
C     EG*LTD(R)2       2         7
C     EG*LTD(R)3       3         7
C     EG*LTD(R)4       3         7
C     EG*LTD(R)5       3         7
C     EG*LTD(R)6       3         7
C
C     EG*TD(R)1        3         7
C
C     EG*V(R)1         0         2
C
C
C     EGINI should be called with the highest possible values for
C     IFLAG and ITLS as read from the table.
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
C     Check values for IFLAG and ITLS
C-----------------------------------------------------------------------
      IERROR = 0
      IF ( IFLAG .LT. 0 .OR. IFLAG .GT. 7 ) THEN
         WRITE(LOUT,'(1X,''IFLAG should be between 0 and 7'')')
         WRITE(LOUT,'(1X,''value read in EGini is'',2x,I5//)') IFLAG
         IERROR = 1
      ENDIF
      IF ( ITLS .LT. 0 .OR. ITLS .GT. 3 ) THEN
         WRITE(LOUT,'(1X,''ITLS should be between 0 and 3'')')
         WRITE(LOUT,'(1X,''value read in EGini is'',2x,I5//)') ITLS
         IERROR = 1
      ENDIF
      IF ( IERROR .EQ. 1 ) STOP
C-----------------------------------------------------------------------
C     Read the Linkeg file
C-----------------------------------------------------------------------
      LLEG   = 11
C-----------------------------------------------------------------------
c      OPEN (UNIT=LLEG,STATUS='OLD',FORM='UNFORMATTED',FILE='Linkeg')
c        READ (LLEG) NSLK, NO
c      CLOSE(UNIT=LLEG)

      call egtransetKK(NSLK)
      call egtransetNO(NO)

C-----------------------------------------------------------------------
C     Store IFLAG and the number of species in common 'eg.cmn'
C-----------------------------------------------------------------------
      JFLAG = IFLAG
C-----------------------------------------------------------------------
      NS = NSLK
      IF ( NS .LE. 1 ) THEN
         WRITE(LOUT,'(1X,''Error: the number of species must '',
     &                   '' be larger or equal to 2'')')
         STOP
      ENDIF
C-----------------------------------------------------------------------
C     Compute the size of the transport linear system.
C-----------------------------------------------------------------------
      NSS  = ITLS * NS
C-----------------------------------------------------------------------
C     NFIT is the degree for the polynomial fitting Aij -- Cij
C-----------------------------------------------------------------------
      NFIT = 7
      NAIJ = MAX0(NS*NS,NP)
C-----------------------------------------------------------------------
      IEGRU  = 1
      IEGPA  = IEGRU  + 1
      IFITA  = IEGPA  + 1
      IFITB  = IFITA  + NFIT * NS*NS
      IFITC  = IFITB  + NFIT * NS*NS
      IFITA0 = IFITC  + NFIT * NS*NS
      IFITB0 = IFITA0 + NFIT
      IFITC0 = IFITB0 + NFIT
      ICTAIJ = IFITC0 + NFIT
      ICTBIJ = ICTAIJ + NS*NS
      ICTCIJ = ICTBIJ + NS*NS
      IDLT1  = ICTCIJ + NS*NS
      IDLT2  = IDLT1  + NP
      IDLT3  = IDLT2  + NP
      IDLT4  = IDLT3  + NP
      IDLT5  = IDLT4  + NP
      IDLT6  = IDLT5  + NP
      IEGEPS = IDLT6  + NP
      IEGPOL = IEGEPS + NS
      IEGSIG = IEGPOL + NS
      IEGDIP = IEGSIG + NS
      IEPSIJ = IEGDIP + NS
      IEGCFD = IEPSIJ + NS*NS
      IEGCFE = IEGCFD + 4 * NS*NS
      IEGCFL = IEGCFE + 4 * NS
      IEGZRT = IEGCFL + 4 * NS
      IEGWT  = IEGZRT + NS
      IAAA   = IEGWT  + NS
      IBBB   = IAAA   + NP
      IAUX   = IBBB   + NP
      IETA   = IAUX   + NS * NP
      IETALG = IETA   + NS * NP
      IXTR   = IETALG + NS * NP
      IYTR   = IXTR   + NS * NP
      IAIJ   = IYTR   + NS * NP
      IBIJ   = IAIJ   + NAIJ
      ICIJ   = IBIJ   + NAIJ
      IBIN   = ICIJ   + NAIJ
      ICINT  = IBIN   + (NS*(NS+1))/2 * NP
      ICXI   = ICINT  + NS * NP
      IEND   = ICXI   + NS * NP
C.....
      IF ( IFLAG .EQ. 1 ) THEN
         IEND = IBIN
      ELSEIF ( IFLAG .LE. 3 ) THEN
         IEND = ICINT
      ENDIF
C.....
      IDMI   = IEND
      IG     = IDMI   + NS*(ITLS*(ITLS+1))/2 * NP
      IAN    = IG     + (ITLS*NS*(ITLS*NS+1))/2 * NP
      IZN    = IAN    + NSS * NP
      IRN    = IZN    + NSS * NP
      ITEMP  = IRN    + NSS * NP
      IBETA  = ITEMP  + NSS * NP
      INEXT  = IBETA  + NSS * NP - 1
C.....
      IEGLIN = 1
      IINXT  = IEGLIN + NS - 1
C-----------------------------------------------------------------------
      ILOW = 0
      IF ( INEXT .GT. LWEG ) THEN
         WRITE(LOUT,'(//1X,''Error: the length of WEG should be '',
     &                   ''at least'',I12//)') INEXT
         ILOW = 1
      ENDIF
      IF ( IINXT .GT. LIWEG ) THEN
         WRITE(LOUT,'(//1X,''Error: the length of IWEG should be '',
     &                   ''at least'',I12//)') IINXT
         ILOW = 1
      ENDIF
      IF ( ILOW .EQ. 1 ) STOP
c     WRITE(LOUT,'(//1X,''The array WEG requires '',I12,
c    &                '' storage locations'')') INEXT
c     WRITE(LOUT,'(1X,''The array IWEG requires '',I12,
c    &                '' storage locations''//)') IINXT
C-----------------------------------------------------------------------
C     Store the universal gas constant and the atmospheric pressure
C     units: [erg/mol.K] for RU and [dyne/cm^2] for PA
C-----------------------------------------------------------------------
      WEG(IEGRU)  = 8.314D7
      WEG(IEGPA)  = 1.01325D6
C-----------------------------------------------------------------------
C     Read the Linkeg file
C-----------------------------------------------------------------------
      LLEG   = 11
      NONS   = NO*NS
      NONSNS = NO*NS*NS
C-----------------------------------------------------------------------
c
c     Set required data from funcs rather than Linkeg
c
      call egtransetWT(WEG(IEGWT))
      call egtransetEPS(WEG(IEGEPS))
      call egtransetSIG(WEG(IEGSIG))
      call egtransetDIP(WEG(IEGDIP))
      call egtransetPOL(WEG(IEGPOL))
      call egtransetZROT(WEG(IEGZRT))
      call egtransetNLIN(IWEG(IEGLIN))
      call egtransetCOFETA(WEG(IEGCFE))
      call egtransetCOFLAM(WEG(IEGCFL))
      call egtransetCOFD(WEG(IEGCFD))

c      OPEN (UNIT=LLEG,STATUS='OLD',FORM='UNFORMATTED',FILE='Linkeg')
c        READ (LLEG) NSLK, NO, (WEG(IEGWT+K-1), K=1, NS),
c     &              (WEG(IEGEPS+K-1), K=1, NS),
c     &              (WEG(IEGSIG+K-1), K=1, NS),
c     &              (WEG(IEGDIP+K-1), K=1, NS),
c     &              (WEG(IEGPOL+K-1), K=1, NS),
c     &              (WEG(IEGZRT+K-1), K=1, NS),
c     &              (IWEG(IEGLIN+K-1), K=1, NS),
c     &              (WEG(IEGCFE+N-1), N=1, NONS),
c     &              (WEG(IEGCFL+N-1), N=1, NONS),
c     &              (WEG(IEGCFD+N-1), N=1, NONSNS)
c      CLOSE(UNIT=LLEG)
C-----------------------------------------------------------------------
      CALL LEVEPS (NS, WEG(IEGEPS), WEG(IEGSIG), WEG(IEGDIP),
     &             WEG(IEGPOL), WEG(IEPSIJ) )
C-----------------------------------------------------------------------
C     Initialize the coefficients for fitting Aij, Bij and Cij
C-----------------------------------------------------------------------
      WEG(IFITA0    ) =  .1106910525D+01
      WEG(IFITA0 + 1) = -.7065517161D-02
      WEG(IFITA0 + 2) = -.1671975393D-01
      WEG(IFITA0 + 3) =  .1188708609D-01
      WEG(IFITA0 + 4) =  .7569367323D-03
      WEG(IFITA0 + 5) = -.1313998345D-02
      WEG(IFITA0 + 6) =  .1720853282D-03
C.....
      WEG(IFITB0    ) =  .1199673577D+01
      WEG(IFITB0 + 1) = -.1140928763D+00
      WEG(IFITB0 + 2) = -.2147636665D-02
      WEG(IFITB0 + 3) =  .2512965407D-01
      WEG(IFITB0 + 4) = -.3030372973D-02
      WEG(IFITB0 + 5) = -.1445009039D-02
      WEG(IFITB0 + 6) =  .2492954809D-03
C.....
      WEG(IFITC0    ) =  .8386993788D+00
      WEG(IFITC0 + 1) =  .4748325276D-01
      WEG(IFITC0 + 2) =  .3250097527D-01
      WEG(IFITC0 + 3) = -.1625859588D-01
      WEG(IFITC0 + 4) = -.2260153363D-02
      WEG(IFITC0 + 5) =  .1844922811D-02
      WEG(IFITC0 + 6) = -.2115417788D-03
C-----------------------------------------------------------------------
C     Evaluate Aij, Bij and Cij at the reference temperature of 1000K.
C-----------------------------------------------------------------------
      DDD = DLOG(1.0D3)
      DO J = 1, NS
        DO I = 1, NS
         IJ = (J-1) * NS + I-1
         TSLOG = DDD - WEG( IEPSIJ + IJ )
         T1 = TSLOG
         T2 = TSLOG * T1
         T3 = TSLOG * T2
         T4 = TSLOG * T3
         T5 = TSLOG * T4
         T6 = TSLOG * T5
         WEG(ICTAIJ+IJ) = WEG(IFITA0  )    + WEG(IFITA0+1)*T1
     1                    + WEG(IFITA0+2)*T2 + WEG(IFITA0+3)*T3
     2                    + WEG(IFITA0+4)*T4 + WEG(IFITA0+5)*T5
     3                    + WEG(IFITA0+6)*T6
         WEG(ICTBIJ+IJ) = WEG(IFITB0  )    + WEG(IFITB0+1)*T1
     1                    + WEG(IFITB0+2)*T2 + WEG(IFITB0+3)*T3
     2                    + WEG(IFITB0+4)*T4 + WEG(IFITB0+5)*T5
     3                    + WEG(IFITB0+6)*T6
         WEG(ICTCIJ+IJ) = WEG(IFITC0  )    + WEG(IFITC0+1)*T1
     1                    + WEG(IFITC0+2)*T2 + WEG(IFITC0+3)*T3
     2                    + WEG(IFITC0+4)*T4 + WEG(IFITC0+5)*T5
     3                    + WEG(IFITC0+6)*T6
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
C     Evaluate FITA, FITB and FITC
C-----------------------------------------------------------------------
      CALL EGABC ( NS, NFIT, WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &             WEG(IFITA0), WEG(IFITB0), WEG(IFITC0),
     &             WEG(IEPSIJ) )
C-----------------------------------------------------------------------
      RETURN
      END


c$$$      subroutine vodeF_T_RhoY(NEQ, TIME, Z, ZP, RPAR, IPAR)
c$$$      implicit none
c$$$      include 'spec.h'
c$$$      double precision TIME, Z(0:Nspec), ZP(0:Nspec), RPAR(*)
c$$$      integer NEQ, IPAR(*)
c$$$
c$$$      double precision RHO, CPB, SUM, Y(Nspec)
c$$$      double precision HK(Nspec), WDOTK(Nspec), C(Nspec), RWRK
c$$$      integer K, IWRK
c$$$
c$$$      integer NiterMAX, Niter
c$$$      parameter (NiterMAX = 30)
c$$$      double precision res(NiterMAX), errMAX, hmix, T
c$$$
c$$$C     For strang
c$$$C     Variables in Z are:  Z(0)   = T
c$$$C                          Z(K) = Y(K)
c$$$C
c$$$C     For SDC
c$$$C     Variables in Z are:  Z(0)   = rho*h
c$$$C                          Z(K) = rho*Y(K)
c$$$
c$$$
c$$$      if (Pcgs.lt.0.d0) then
c$$$         print *,'vodeF_T_RhoY: Must set Pcgs before calling vode'
c$$$         stop
c$$$      endif
c$$$
c$$$         RHO = 0.d0
c$$$         do K=1,Nspec
c$$$            RHO = RHO + Z(K)
c$$$         enddo
c$$$
c$$$         do K=1,Nspec
c$$$            if (lim_rxns .eq. 0) then
c$$$               C(K) = Z(K)*invmwt(K)
c$$$            else
c$$$               C(K) = MAX(0.d0,Z(K)*invmwt(K))
c$$$            endif
c$$$            Y(K) = Z(K)/RHO
c$$$         enddo
c$$$
c$$$         T = T_INIT
c$$$         hmix = (rhoh_INIT + c_0(0)*
c$$$     &        TIME + c_1(0)*TIME*TIME*0.5d0)/RHO
c$$$         errMax = hmix_TYP*1.e-20
c$$$         call FORT_TfromHYpt(T,hmix,Y,Nspec,errMax,NiterMAX,
c$$$     &                       res,Niter)
c$$$         T_INIT = T
c$$$         if (Niter.lt.0) then
c$$$            print *,'vodeF_T_RhoY: H to T solve failed'
c$$$            print *,'Niter=',Niter
c$$$            stop
c$$$         endif
c$$$
c$$$         call CKWC(T,C,IWRK,RWRK,WDOTK)
c$$$         do k= 1, Nspec
c$$$            ZP(k) = WDOTK(k)*mwt(k) + c_0(k) + c_1(k)*TIME
c$$$         end do
c$$$         ZP(0) = c_0(0) + c_1(0)*TIME
c$$$
c$$$      END
c$$$
c$$$      subroutine vodeJ(NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
c$$$      implicit none
c$$$      integer NEQ, NRPD, ML, MU, IPAR(*)
c$$$      double precision T, Y(NEQ), PD(NRPD,NEQ), RPAR(*)
c$$$      print *,'Should not be in vodeJ'
c$$$      stop
c$$$      end
c$$$
c$$$
c$$$
c$$$
c$$$      subroutine chemsolve(RYnew, Tnew, RYold, Told, FuncCount, dt,
c$$$     &     diag, do_diag, ifail, i)
c$$$      implicit none
c$$$      include 'spec.h'
c$$$
c$$$      double precision YJ_SAVE(80)
c$$$      logical FIRST
c$$$      common /VHACK/ YJ_SAVE, FIRST
c$$$      save   /VHACK/
c$$$
c$$$      integer do_diag, ifail, FuncCount, i
c$$$      double precision RYold(*), RYnew(*), Told, Tnew
c$$$      double precision dt, diag(*)
c$$$
c$$$      integer NEQ, ITOL, IOPT, ITASK, open_vode_failure_file
c$$$      parameter (ITOL=1, IOPT=1, ITASK=1)
c$$$      double precision RTOL, ATOL(Nspec+1), ATOLEPS, TT1, TT2
c$$$C      parameter (RTOL=1.0E-8, ATOLEPS=1.0E-8)
c$$$      parameter (RTOL=1.0D-13, ATOLEPS=1.0D-13)
c$$$      external vodeF_T_RhoY, vodeJ, open_vode_failure_file
c$$$      integer n, MF, ISTATE
c$$$
c$$$      integer nsubchem, nsub, node
c$$$      double precision dtloc, weight, TT1save
c$$$      double precision C(Nspec),Q(Nreac)
c$$$
c$$$c     DVODE workspace requirements
c$$$      integer dvr, dvi
c$$$
c$$$c     METH = 2 (backward diffs)
c$$$c     MAXORD = 5
c$$$c     NEQ = Nspec+1
c$$$c     NYH = NEQ = Nspec+1
c$$$c     MITER = 2
c$$$c     JSV = 1
c$$$c     JSV = SIGN(MF)
c$$$c     MF = JSV*(10*METH + MITER) = 22
c$$$c     LWM = 2*(Nspec+1)**2 + 2    (MITER = 2, MF<0)
c$$$c     lenr = 20 + NYH*(MAXORD + 1) + 3*NEQ + LWM
c$$$c          = 20 + (Nspec+1)*(6) + 3*(Nspec+1) + 2*(Nspec+1)**2 + 2
c$$$c          = 22 + (Nspec+1)*9 + 2*(Nspec+1)**2
c$$$c
c$$$      parameter (dvr = 22 + 9*(maxspec+1) + 2*(maxspec+1)**2)
c$$$      parameter (dvi = 30 + maxspec + 1)
c$$$
c$$$      double precision DVRWRK(dvr)
c$$$      integer DVIWRK(dvi)
c$$$
c$$$      common / VODE_SPACE / dvrwrk, dviwrk
c$$$      save / VODE_SPACE /
c$$$
c$$$      double precision Z(0:Nspec)
c$$$      double precision RPAR, RWRK
c$$$      integer IPAR, IWRK
c$$$
c$$$c     IOPT=1 parameter settings for VODE
c$$$C      DVRWRK(1) = 0.d0
c$$$C      DVRWRK(2) = 0.d0
c$$$C      DVRWRK(3) = 0.d0
c$$$C      DVRWRK(4) = 0.d0
c$$$      DVRWRK(5) = 0.d0
c$$$      DVRWRK(6) = 0.d0
c$$$      DVRWRK(7) = min_vode_timestep
c$$$C      DVRWRK(8) = 0.d0
c$$$C      DVRWRK(9) = 0.d0
c$$$C      DVRWRK(10) = 0.d0
c$$$C      DVIWRK(1) = 0
c$$$C      DVIWRK(2) = 0
c$$$C      DVIWRK(3) = 0
c$$$C      DVIWRK(4) = 0
c$$$      DVIWRK(5) = 0
c$$$      DVIWRK(6) = max_vode_subcycles
c$$$      DVIWRK(7) = 0
c$$$C      DVIWRK(8) = 0
c$$$C      DVIWRK(9) = 0
c$$$C      DVIWRK(10) = 0
c$$$
c$$$      if (do_diag.eq.1) nsubchem = nchemdiag
c$$$
c$$$
c$$$      MF = 22
c$$$      ATOL(1) = ATOLEPS
c$$$      TT1 = 0.d0
c$$$      TT2 = dt
c$$$      if (do_diag.eq.1) then
c$$$         nsub = nsubchem
c$$$         dtloc = dt/nsubchem
c$$$      else
c$$$         nsub = 1
c$$$         dtloc = dt
c$$$      endif
c$$$      ISTATE = 1
c$$$      NEQ = Nspec + 1
c$$$
c$$$
c$$$      Z(0) = Told
c$$$      do n=1,Nspec
c$$$         Z(n) = RYold(n)
c$$$      end do
c$$$
c$$$c     Always form Jacobian to start
c$$$      if (i .eq. 0) then
c$$$         FIRST = .TRUE.
c$$$      else
c$$$         FIRST = .FALSE.
c$$$      end if
c$$$
c$$$      ifail = 0
c$$$      do node = 1,nsub
c$$$         if (node.lt.nsub) then
c$$$            weight = 1.d0
c$$$         else
c$$$            weight = 0.5d0
c$$$         endif
c$$$
c$$$         TT1save = TT1
c$$$         TT2 = TT1 + dtloc
c$$$
c$$$c     HACK
c$$$         FIRST = .TRUE.
c$$$
c$$$         CALL DVODE
c$$$     &        (vodeF_T_RhoY, NEQ, Z(0), TT1, TT2, ITOL, RTOL, ATOL,
c$$$     &        ITASK, ISTATE, IOPT, DVRWRK, dvr, DVIWRK,
c$$$     &        dvi, vodeJ, MF, RPAR, IPAR)
c$$$
c$$$         TT1 = TT2
c$$$
c$$$         FuncCount = DVIWRK(11)
c$$$
c$$$         if (ISTATE.LE.-1  .or.  verbose_vode .eq. 1) then
c$$$            write(6,*) '......dvode done:'
c$$$            write(6,*) ' last successful step size = ',DVRWRK(11)
c$$$            write(6,*) '          next step to try = ',DVRWRK(12)
c$$$            write(6,*) '   integrated time reached = ',DVRWRK(13)
c$$$            write(6,*) '      number of time steps = ',DVIWRK(11)
c$$$            write(6,*) '              number of fs = ',DVIWRK(12)
c$$$            write(6,*) '              number of Js = ',DVIWRK(13)
c$$$            write(6,*) '    method order last used = ',DVIWRK(14)
c$$$            write(6,*) '   method order to be used = ',DVIWRK(15)
c$$$            write(6,*) '            number of LUDs = ',DVIWRK(19)
c$$$            write(6,*) ' number of Newton iterations ',DVIWRK(20)
c$$$            write(6,*) ' number of Newton failures = ',DVIWRK(21)
c$$$            write(6,*) '     comp with largest err = ',DVIWRK(16)
c$$$
c$$$            print *,'Input state:'
c$$$            print *,'T:',Told
c$$$            print *,'RY:',(RYold(n),n=1,Nspec)
c$$$         end if
c$$$
c$$$         Tnew = Z(0)
c$$$         do n=1,Nspec
c$$$            RYnew(n) = Z(n)
c$$$         end do
c$$$
c$$$         if (ISTATE.LE.-1) ifail = 1
c$$$
c$$$      enddo
c$$$      end


      subroutine FORT_TfromHYpt(T,Hin,Y,Nspec,errMax,NiterMAX,res,Niter)
      implicit none
      double precision T,Y(*),H,Hin
      double precision TMIN,TMAX,errMAX
      integer Nspec,NiterMAX,Niter,n,NiterDAMP, Discont_NiterMAX
      parameter (TMIN=250, TMAX=5000, Discont_NiterMAX=100)
      double precision  T0,cp,dH,temp
      double precision res(0:NiterMAX-1),dT, Htarg
      logical out_of_bounds, converged, soln_bad, stalled
      double precision HMIN,cpMIN,HMAX,cpMAX
      integer ihitlo,ihithi,IWRK
      double precision old_T, old_H, Tsec, Hsec, RWRK

      out_of_bounds(temp) = (temp.lt.TMIN-1.d0) .or. (temp.gt.TMAX)

      NiterDAMP = NiterMAX
      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
         T0 = T
      else
         T0 = 0.5d0*(TMIN+TMAX)
         T = T0
      end if
      Niter = 0
      soln_bad = .FALSE.
c     Hin in MKS, convert to CGS
c      Htarg = Hin * 1.d4
      Htarg = Hin
      ihitlo = 0
      ihithi = 0

      CALL CKHBMS(T,Y,IWRK,RWRK,H)

      old_T = T
      old_H = H

      Htarg = Hin

      dH = 2.d0*ABS(H - Htarg)/(1.d0 + ABS(H) + ABS(Htarg))
      res(Niter) = dH
      converged = dH.le.errMAX
      stalled = .false.

      do while ((.not.converged) .and.
     &     (.not.stalled) .and.
     &     (.not.soln_bad))

         CALL CKCPBS(T,Y,IWRK,RWRK,cp)
         dT = (Htarg - H)/cp
         old_T = T
         if ((Niter.le.NiterDAMP).and.(T+dT.ge.TMAX)) then
            T = TMAX
            ihithi = 1
         else if ((Niter.le.NiterDAMP).and.(T+dT.le.TMIN)) then
            T = TMIN
            ihitlo = 1
         else
            T = T + dT
         end if
         soln_bad = out_of_bounds(T)
         if (soln_bad) then
            Niter = -1
            goto 100
         else
            old_H = H
            CALL CKHBMS(T,Y,IWRK,RWRK,H)
            dH = 2.d0*ABS(H - Htarg)/(1.d0 + ABS(H) + ABS(Htarg))
            res(Niter) = min(dH,abs(dT))
            Niter = Niter + 1
         end if
         converged = (dH.le.errMAX)
         stalled = (ABS(dT).le.errMAX)
         if ((ihitlo.eq.1).and.(H.gt.Htarg)) then
            T = TMIN
            CALL CKHBMS(T,Y,IWRK,RWRK,HMIN)
            CALL CKCPBS(T,Y,IWRK,RWRK,cpMIN)
            T=TMIN+(Htarg-HMIN)/cpMIN
            converged = .true.
         endif
         if ((ihithi.eq.1).and.(H.lt.Htarg)) then
            T = TMAX
            CALL CKHBMS(T,Y,IWRK,RWRK,HMAX)
            CALL CKCPBS(T,Y,IWRK,RWRK,cpMAX)
            T=TMAX+(Htarg-HMAX)/cpMAX
            converged = .true.
         endif

c     If the iterations are failing, perhaps it is because the fits are discontinuous
c     The following implements a modified secant method to hone in on a discontinity in h
c     with T.  Discont_NiterMAX is fairly large because this process can be particularly
c     slow to converge if the Htarg value happens to lay between the discontinuous function
c     values.
         if (Niter .ge. NiterMAX) then
            do while (.not. stalled)
               dT = - (H - Htarg) * (old_T - T)/(old_H - H)
               Tsec = T + dT
               soln_bad = out_of_bounds(Tsec)
               if (soln_bad) then
                  Niter = -3
                  goto 100
               endif
               CALL CKHBMS(Tsec,Y,IWRK,RWRK,Hsec)
               if ( (Hsec-Htarg)*(Htarg-H) .gt. 0.d0 ) then
                  old_H = H
                  old_T = T
               endif
               H = Hsec
               T = Tsec
               stalled = (ABS(dT).le.errMAX)
               Niter = Niter + 1
               if (Niter.gt.NiterMAX+Discont_NiterMAX) then
                  Niter = -2
                  goto 100
               endif
            enddo
            converged = .true.
         endif
      end do

      return
c
c     Error condition....dump state and bail out
c
 100  continue

      write(6,997) 'T from (H,Y): failed'
      write(6,997) 'iterations tried = ',Niter
      write(6,998) 'initial T = ',T0
      write(6,998) 'current T = ',T
      write(6,998) 'previous T = ',old_T
      write(6,998) 'target H = ',Htarg
      write(6,998) 'species mass fracs:'
      do n = 1,Nspec
         write(6,998) '  ',Y(n)
      end do
      write(6,998)
      write(6,998) 'residual:'
      do n = 0,NiterMAX-1
         write(6,998) '  ',res(n)
      end do

 997  format(a,3(i4,a))
 998  format(a,d21.12)
      end

      integer function open_vode_failure_file ()
      implicit none
      character*30 name, myproc
      integer lout,i,j,k,idx

c     Hardwire the unit number to 26 for the moment
      lout = 26
c      call bl_pd_myproc(i)
      i = 0
      write(myproc, *) i
      idx = 1
      do j = 1, 30
         if (myproc(j:j) .ne. ' ') then
            idx = j
            goto 1
         end if
      end do
 1    continue
      do k = 30, j+1, -1
         if (myproc(k:k) .ne. ' ') then
            goto 2
         end if
      end do
 2    continue
      write(name, '(2a)') 'vode.failed.', myproc(idx:k)
c      write(name, '(2a)') 'vode.failed.', myproc(idx:30)
      open(unit=lout, file=name, form='formatted', status='replace')
      open_vode_failure_file = lout
      end

