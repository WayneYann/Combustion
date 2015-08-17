      block data chemdat
      include 'spec.h'
      data traninit / -1 /
      data Pcgs / -1 /
      data iH2  / -1 /
      data iO2  / -1 /
      data iN2  / -1 /
      data iCH4 / -1 /
      end

      subroutine compute_diffusion_coefficients(beta, scal)
         implicit none
         
         include 'spec.h'
         
         double precision, intent(out  ) :: beta(-2:nx+1,nscal)
         double precision, intent(in   ) :: scal(-2:nx+1,nscal)
         
         double precision Dt(Nspec), CPMS(Nspec), Y(Nspec)
         double precision Tt, Wavg, rho
         double precision X(Nspec), alpha, l1, l2, cpmix, RWRK
         integer n, i, IWRK

         double precision fourThirds

         fourThirds = 4.d0 / 3.d0

   !     Ensure chem/tran initialized
         if (traninit.lt.0) call initchem()
         
         do i=-2,nx+1
            Tt = MAX(scal(i,Temp),TMIN_TRANS) 
            rho = 0.d0
            do n=1,Nspec
               rho = rho + scal(i,FirstSpec+n-1)
            enddo
            do n=1,Nspec
               Y(n) = scal(i,FirstSpec+n-1) / rho
            enddo
            
            !  given y[species]: mass fractions
            !  returns mean molecular weight (gm/mole)
            CALL CKMMWY(Y,IWRK,RWRK,Wavg)

            !  returns the specific heats at constant pressure
            !  in mass units
            CALL CKCPMS(Tt,IWRK,RWRK,CPMS)

            !  convert y[species] (mass fracs) to x[species] (mole fracs)
            CALL CKYTX(Y,IWRK,RWRK,X)

            ! initialize the thermomolecular parameters that are needed in 
            ! to evaluate the transport linear systems
            CALL EGSPAR(Tt,X,Y,CPMS,EGRWRK,EGIWRK)

            ! compute flux diffusion coefficients
            CALL EGSV1(Pcgs,Tt,Y,Wavg,EGRWRK,Dt)
            
            do n=1,Nspec
               beta(i,FirstSpec+n-1) = rho * Wavg * invmwt(n) * Dt(n)
            end do

            alpha = 1.0D0
            ! compute thermal conductivity
            CALL EGSL1(alpha, Tt, X, EGRWRK, l1)
            alpha = -1.0D0
            !  compute thermal conductivity with a different averating para
            CALL EGSL1(alpha, Tt, X, EGRWRK, l2)
            beta(i,Temp) = .5d0 * (l1 + l2)
            !  Returns the mean specific heat at CP
            CALL CKCPBS(scal(i,Temp),Y,IWRK,RWRK,CPMIX)
            beta(i,RhoH) = beta(i,Temp) / CPMIX
         enddo
         
      end subroutine compute_diffusion_coefficients
      
      subroutine initchem
      implicit none
      include 'spec.h'
      double precision RU, RUC, RWRK, RDUMMY
      integer lout, IWRK, n, len, offset, IDUMMY
      integer kname(maxspnml*maxspec)

      if (traninit.lt.0) then
!         open(unit=51,status='old',form='formatted',file=tranfile)

!     This sets Nspec
         CALL CKINIT()
         CALL CKINDX(IDUMMY,RDUMMY,Nelt,Nspec,Nreac,Nfit)

         lout = 6
         call EGINICD(eg_nodes, lout, eg_IFLAG, eg_ITLS,&
     &        EGRWRK, egr, EGIWRK, egi)

!     Other useful things
         call CKRP(IWRK, RWRK, RU, RUC, P1ATM)
         call CKWT(IWRK, RWRK, mwt)
     
         do n=1,Nspec
            invmwt(n) = 1.d0 / mwt(n)
         end do

         call cksyms(kname,maxspnml)
         do n=1,Nspec
            offset = (n-1)*maxspnml+1
            call convStr(kname(offset),maxspnml,specNames(n),len)
            if (specNames(n).eq.'H2')  iH2=n
            if (specNames(n).eq.'O2')  iO2=n
            if (specNames(n).eq.'N2')  iN2=n
            if (specNames(n).eq.'CH4') iCH4=n
            if (specNames(n).eq.'CH3OCH3') iCH3OCH3=n
            if (specNames(n).eq.'CO2') iCO2=n
            if (specNames(n).eq.'H2O') iH2O=n
         enddo
         traninit = 1
      endif
      end

      subroutine convStr(codedString,maxLen,string,strLen)
      implicit none
      integer codedString(*), maxLen, k, strLen
      character*(*) string
      do k=1,maxLen
         string(k:k) = ' '
         if (codedString(k) .ne. ICHAR(' ')) then
            string(k:k) = CHAR(codedString(k))
            strLen = k
         endif
      enddo
      end

      SUBROUTINE EGINICD (NP, LOUT, IFLAG, ITLS, &
     &                    WEG, LWEG, IWEG, LIWEG)
!-----------------------------------------------------------------------
!
!     This subroutine initializes the pointers for the work arrays
!     WEG and IWEG and checks their length.
!     This subroutine should be called by the user once at the
!     beginning of the program.
!
!     Input
!     -----
!        NP        number of nodes
!        LOUT      output file number
!        IFLAG     flag for evaluating parameters and space allocation
!                  (see below)
!        ITLS      flag for space allocation (see below)
!        WEG       double precision work array for EGLIB
!        LWEG      length of WEG declared in main code
!        IWEG      integer work array for EGLIB
!        LIWEG     length of IWEG declared in main code
!
!        
!     The value of IFLAG and ITLS depends on the subroutines that
!     will be used as indicated by the following table
!
!
!     Subroutine     ITLS      IFLAG
!
!     EG*D(R)1         1         2
!     EG*D(R)2         1         2
!
!     EG*E1            0         1
!     EG*E2            1         2
!     EG*E3            1         3
!     EG*E4            1         3
! 
!     EG*K1            0         4
!     EG*K2            1         4
!     EG*K3            1         5
!     EG*K4            2         4
!     EG*K5            2         5
!     EG*K6            2         5
!  
!     EG*L1            0         1
!     EG*L2            1         6
!     EG*L3            1         7
!     EG*L4            2         6
!     EG*L5            2         7
!  
!     EG*LC1           1         7
!     EG*LC2           1         7
!     EG*LC3           2         7
!     EG*LC4           2         7
!  
!     EG*LTD(R)1       2         7
!     EG*LTD(R)2       2         7
!     EG*LTD(R)3       3         7
!     EG*LTD(R)4       3         7
!     EG*LTD(R)5       3         7
!     EG*LTD(R)6       3         7
!   
!     EG*TD(R)1        3         7
!  
!     EG*V(R)1         0         2
!
!
!     EGINI should be called with the highest possible values for
!     IFLAG and ITLS as read from the table.
!
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
!-----------------------------------------------------------------------
!     Check values for IFLAG and ITLS
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!     Read the Linkeg file
!-----------------------------------------------------------------------
      LLEG   = 11
!-----------------------------------------------------------------------
!      OPEN (UNIT=LLEG,STATUS='OLD',FORM='UNFORMATTED',FILE='Linkeg')
!        READ (LLEG) NSLK, NO
!      CLOSE(UNIT=LLEG)

      call egtransetKK(NSLK)
      call egtransetNO(NO)

!-----------------------------------------------------------------------
!     Store IFLAG and the number of species in common 'eg.cmn'
!-----------------------------------------------------------------------
      JFLAG = IFLAG
!-----------------------------------------------------------------------
      NS = NSLK
      IF ( NS .LE. 1 ) THEN
         WRITE(LOUT,'(1X,''Error: the number of species must '',&
     &                   '' be larger or equal to 2'')')
         STOP
      ENDIF
!-----------------------------------------------------------------------
!     Compute the size of the transport linear system.
!-----------------------------------------------------------------------
      NSS  = ITLS * NS
!-----------------------------------------------------------------------
!     NFIT is the degree for the polynomial fitting Aij -- Cij
!-----------------------------------------------------------------------
      NFIT = 7
      NAIJ = MAX0(NS*NS,NP)
!-----------------------------------------------------------------------
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
!.....
      IF ( IFLAG .EQ. 1 ) THEN
         IEND = IBIN
      ELSEIF ( IFLAG .LE. 3 ) THEN
         IEND = ICINT
      ENDIF
!.....
      IDMI   = IEND
      IG     = IDMI   + NS*(ITLS*(ITLS+1))/2 * NP
      IAN    = IG     + (ITLS*NS*(ITLS*NS+1))/2 * NP
      IZN    = IAN    + NSS * NP
      IRN    = IZN    + NSS * NP
      ITEMP  = IRN    + NSS * NP
      IBETA  = ITEMP  + NSS * NP
      INEXT  = IBETA  + NSS * NP - 1
!.....
      IEGLIN = 1
      IINXT  = IEGLIN + NS - 1
!-----------------------------------------------------------------------
      ILOW = 0
      IF ( INEXT .GT. LWEG ) THEN
         WRITE(LOUT,'(//1X,''Error: the length of WEG should be '',&
     &                   ''at least'',I12//)') INEXT
         ILOW = 1
      ENDIF
      IF ( IINXT .GT. LIWEG ) THEN
         WRITE(LOUT,'(//1X,''Error: the length of IWEG should be '',&
     &                   ''at least'',I12//)') IINXT
         ILOW = 1
      ENDIF
      IF ( ILOW .EQ. 1 ) STOP
!     WRITE(LOUT,'(//1X,''The array WEG requires '',I12,
!    &                '' storage locations'')') INEXT
!     WRITE(LOUT,'(1X,''The array IWEG requires '',I12,
!    &                '' storage locations''//)') IINXT
!-----------------------------------------------------------------------
!     Store the universal gas constant and the atmospheric pressure
!     units: [erg/mol.K] for RU and [dyne/cm^2] for PA
!-----------------------------------------------------------------------
      WEG(IEGRU)  = 8.314D7
      WEG(IEGPA)  = 1.01325D6
!-----------------------------------------------------------------------
!     Read the Linkeg file
!-----------------------------------------------------------------------
      LLEG   = 11
      NONS   = NO*NS
      NONSNS = NO*NS*NS
!-----------------------------------------------------------------------
!
!     Set required data from funcs rather than Linkeg
!
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

!      OPEN (UNIT=LLEG,STATUS='OLD',FORM='UNFORMATTED',FILE='Linkeg')
!        READ (LLEG) NSLK, NO, (WEG(IEGWT+K-1), K=1, NS), 
!     &              (WEG(IEGEPS+K-1), K=1, NS), 
!     &              (WEG(IEGSIG+K-1), K=1, NS), 
!     &              (WEG(IEGDIP+K-1), K=1, NS), 
!     &              (WEG(IEGPOL+K-1), K=1, NS), 
!     &              (WEG(IEGZRT+K-1), K=1, NS), 
!     &              (IWEG(IEGLIN+K-1), K=1, NS), 
!     &              (WEG(IEGCFE+N-1), N=1, NONS), 
!     &              (WEG(IEGCFL+N-1), N=1, NONS), 
!     &              (WEG(IEGCFD+N-1), N=1, NONSNS)
!      CLOSE(UNIT=LLEG)
!-----------------------------------------------------------------------
      CALL LEVEPS (NS, WEG(IEGEPS), WEG(IEGSIG), WEG(IEGDIP), &
     &             WEG(IEGPOL), WEG(IEPSIJ) )
!-----------------------------------------------------------------------
!     Initialize the coefficients for fitting Aij, Bij and Cij
!-----------------------------------------------------------------------
      WEG(IFITA0    ) =  .1106910525D+01
      WEG(IFITA0 + 1) = -.7065517161D-02
      WEG(IFITA0 + 2) = -.1671975393D-01
      WEG(IFITA0 + 3) =  .1188708609D-01
      WEG(IFITA0 + 4) =  .7569367323D-03
      WEG(IFITA0 + 5) = -.1313998345D-02
      WEG(IFITA0 + 6) =  .1720853282D-03
!.....
      WEG(IFITB0    ) =  .1199673577D+01
      WEG(IFITB0 + 1) = -.1140928763D+00
      WEG(IFITB0 + 2) = -.2147636665D-02
      WEG(IFITB0 + 3) =  .2512965407D-01
      WEG(IFITB0 + 4) = -.3030372973D-02
      WEG(IFITB0 + 5) = -.1445009039D-02
      WEG(IFITB0 + 6) =  .2492954809D-03
!.....
      WEG(IFITC0    ) =  .8386993788D+00
      WEG(IFITC0 + 1) =  .4748325276D-01
      WEG(IFITC0 + 2) =  .3250097527D-01
      WEG(IFITC0 + 3) = -.1625859588D-01
      WEG(IFITC0 + 4) = -.2260153363D-02
      WEG(IFITC0 + 5) =  .1844922811D-02
      WEG(IFITC0 + 6) = -.2115417788D-03
!-----------------------------------------------------------------------
!     Evaluate Aij, Bij and Cij at the reference temperature of 1000K.
!-----------------------------------------------------------------------
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
         WEG(ICTAIJ+IJ) = WEG(IFITA0  )    + WEG(IFITA0+1)*T1&
     &                    + WEG(IFITA0+2)*T2 + WEG(IFITA0+3)*T3&
     &                    + WEG(IFITA0+4)*T4 + WEG(IFITA0+5)*T5&
     &                    + WEG(IFITA0+6)*T6 
         WEG(ICTBIJ+IJ) = WEG(IFITB0  )    + WEG(IFITB0+1)*T1&
     &                    + WEG(IFITB0+2)*T2 + WEG(IFITB0+3)*T3&
     &                    + WEG(IFITB0+4)*T4 + WEG(IFITB0+5)*T5&
     &                    + WEG(IFITB0+6)*T6 
         WEG(ICTCIJ+IJ) = WEG(IFITC0  )    + WEG(IFITC0+1)*T1&
     &                    + WEG(IFITC0+2)*T2 + WEG(IFITC0+3)*T3&
     &                    + WEG(IFITC0+4)*T4 + WEG(IFITC0+5)*T5&
     &                    + WEG(IFITC0+6)*T6 
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     Evaluate FITA, FITB and FITC 
!-----------------------------------------------------------------------
      CALL EGABC ( NS, NFIT, WEG(IFITA), WEG(IFITB), WEG(IFITC),&
     &             WEG(IFITA0), WEG(IFITB0), WEG(IFITC0),&
     &             WEG(IEPSIJ) )
!-----------------------------------------------------------------------
      RETURN
      END


      subroutine vodeF_T_RhoY(NEQ, TIME, Z, ZP, RPAR, IPAR)
      implicit none
      include 'spec.h'
      double precision TIME, Z(0:Nspec), ZP(0:Nspec), RPAR(*)
      integer NEQ, IPAR(*)
      
      double precision RHO, Y(Nspec)
      double precision WDOTK(Nspec), C(Nspec), RWRK
      integer K, IWRK

      integer NiterMAX, Niter
      parameter (NiterMAX = 30)
      double precision res(NiterMAX), errMAX, hmix, T

!     For strang
!     Variables in Z are:  Z(0)   = T
!                          Z(K) = Y(K)
!
!     For SDC
!     Variables in Z are:  Z(0)   = rho*h
!                          Z(K) = rho*Y(K)


      if (Pcgs.lt.0.d0) then
         print *,'vodeF_T_RhoY: Must set Pcgs before calling vode'
         stop
      endif

      RHO = 0.d0
      do K=1,Nspec
         RHO = RHO + Z(K)
      enddo
      
      do K=1,Nspec
         if (lim_rxns .eq. 0) then
            C(K) = Z(K)*invmwt(K)
         else
            C(K) = MAX(0.d0,Z(K)*invmwt(K))
         endif
         Y(K) = Z(K)/RHO
      enddo

      T = T_INIT
      hmix = (rhoh_INIT + c_0(0)*&
     &        TIME + c_1(0)*TIME*TIME*0.5d0)/RHO
      errMax = hmix_TYP*1.e-20
      call FORT_TfromHYpt(T,hmix,Y,Nspec,errMax,NiterMAX,&
     &                       res,Niter)
      T_INIT = T
      if (Niter.lt.0) then
         print *,'vodeF_T_RhoY: H to T solve failed'
         print *,'Niter=',Niter
         stop
      endif

      call CKWC(T,C,IWRK,RWRK,WDOTK)
      
      do k= 1, Nspec
         ZP(k) = WDOTK(k)*mwt(k) + c_0(k) + c_1(k)*TIME
      end do
      ZP(0) = c_0(0) + c_1(0)*TIME 

      END

      subroutine vodeJ(NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      implicit none
      integer NEQ, NRPD, ML, MU, IPAR(*)
      double precision T, Y(NEQ), PD(NRPD,NEQ), RPAR(*)
      print *,'Should not be in vodeJ'
      stop
      end


      subroutine chemsolve(RYnew, Tnew, RYold, Told, FuncCount, dt,&
     &     diag, do_diag, ifail, i)
      implicit none
      include 'spec.h'

      double precision YJ_SAVE(80)
      logical FIRST
      common /VHACK/ YJ_SAVE, FIRST
      save   /VHACK/

      integer do_diag, ifail, FuncCount, i
      double precision RYold(*), RYnew(*), Told, Tnew
      double precision dt, diag(*)
   
      integer NEQ, ITOL, IOPT, ITASK, open_vode_failure_file
      parameter (ITOL=1, IOPT=1, ITASK=1)
      double precision RTOL, ATOL(Nspec+1), ATOLEPS, TT1, TT2
      parameter (RTOL=1.0D-13, ATOLEPS=1.0D-13)
      external vodeF_T_RhoY, vodeJ, open_vode_failure_file
      integer n, MF, ISTATE

      integer nsubchem, nsub, node
      double precision dtloc, weight

!     DVODE workspace requirements      
      integer dvr, dvi
!     
      parameter (dvr = 22 + 9*(maxspec+1) + 2*(maxspec+1)**2)
      parameter (dvi = 30 + maxspec + 1)
      
      double precision DVRWRK(dvr)
      integer DVIWRK(dvi)

      common / VODE_SPACE / dvrwrk, dviwrk
      save / VODE_SPACE /

      double precision Z(0:Nspec)
      
      double precision RPAR
      integer IPAR

!     IOPT=1 parameter settings for VODE
!      DVRWRK(1) = 0.d0
!      DVRWRK(2) = 0.d0
!      DVRWRK(3) = 0.d0
!      DVRWRK(4) = 0.d0
      DVRWRK(5) = 0.d0
      DVRWRK(6) = 0.d0
      DVRWRK(7) = min_vode_timestep
!      DVRWRK(8) = 0.d0
!      DVRWRK(9) = 0.d0
!      DVRWRK(10) = 0.d0
!      DVIWRK(1) = 0
!      DVIWRK(2) = 0
!      DVIWRK(3) = 0
!      DVIWRK(4) = 0
!      DVIWRK(5) = 0
      DVIWRK(5) = 2
      DVIWRK(6) = max_vode_subcycles
      DVIWRK(7) = 0
!      DVIWRK(8) = 0
!      DVIWRK(9) = 0
!      DVIWRK(10) = 0

      if (do_diag.eq.1) nsubchem = nchemdiag

      
      MF = 22
!      MF = 10
      
      ATOL(1) = ATOLEPS
      TT1 = 0.d0
      
      if (do_diag.eq.1) then
         nsub = nsubchem
         dtloc = dt/nsubchem
      else
         nsub = 1
         dtloc = dt
      endif
      
      ISTATE = 1
      NEQ = Nspec + 1


      Z(0) = Told
      do n=1,Nspec
         Z(n) = RYold(n)
      end do

      ifail = 0
      do node = 1,nsub
         if (node.lt.nsub) then
            weight = 1.d0
         else
            weight = 0.5d0
         endif
         
         TT2 = TT1 + dtloc

!     HACK
         FIRST = .TRUE.

         IPAR = i

         CALL DVODE&
     &        (vodeF_T_RhoY, NEQ, Z, TT1, TT2, ITOL, RTOL, ATOLEPS,&
     &        ITASK, ISTATE, IOPT, DVRWRK, dvr, DVIWRK,&
     &        dvi, vodeJ, MF, RPAR, IPAR)
         
         TT1 = TT2
         
         FuncCount = DVIWRK(11)

         if (ISTATE.LE.-1  .or.  verbose_vode .eq. 1) then
            write(6,*) '......dvode done:'
            write(6,*) ' last successful step size = ',DVRWRK(11)
            write(6,*) '          next step to try = ',DVRWRK(12)
            write(6,*) '   integrated time reached = ',DVRWRK(13)
            write(6,*) '      number of time steps = ',DVIWRK(11)
            write(6,*) '              number of fs = ',DVIWRK(12)
            write(6,*) '              number of Js = ',DVIWRK(13)
            write(6,*) '    method order last used = ',DVIWRK(14)
            write(6,*) '   method order to be used = ',DVIWRK(15)
            write(6,*) '            number of LUDs = ',DVIWRK(19)
            write(6,*) ' number of Newton iterations ',DVIWRK(20)
            write(6,*) ' number of Newton failures = ',DVIWRK(21)
            write(6,*) '     comp with largest err = ',DVIWRK(16)

            print *,'Input state:'
            print *,'T:',Told
            print *,'RY:',(RYold(n),n=1,Nspec)
         end if

         Tnew = Z(0)
         do n=1,Nspec
            RYnew(n) = Z(n)
         end do

         if (ISTATE.LE.-1) ifail = 1

      enddo
      
      end


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
!     Hin in MKS, convert to CGS
!      Htarg = Hin * 1.d4
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

      do while ((.not.converged) .and. &
     &     (.not.stalled) .and.&
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

!     If the iterations are failing, perhaps it is because the fits are 
!     The following implements a modified secant method to hone in on a 
!     with T.  Discont_NiterMAX is fairly large because this process can
!     slow to converge if the Htarg value happens to lay between the dis
!     values.  
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
!
!     Error condition....dump state and bail out
!
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

!     Hardwire the unit number to 26 for the moment
      lout = 26 
!      call bl_pd_myproc(i)
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
!      write(name, '(2a)') 'vode.failed.', myproc(idx:30)
      open(unit=lout, file=name, form='formatted', status='replace')
      open_vode_failure_file = lout
      end
