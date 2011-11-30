      subroutine settmintrans(TminTRANS)

      use cdwrk_module
      implicit none

      double precision :: TminTRANS
      TMIN_TRANS = TminTRANS

      end subroutine settmintrans

      subroutine setverbosevode()

      use cdwrk_module
      implicit none

      verbose_vode = 1

      end subroutine setverbosevode

      subroutine setvodesubcyc(maxcyc)

      use cdwrk_module
      implicit none

      integer maxcyc

      max_vode_subcycles = maxcyc

      end subroutine setvodesubcyc

      subroutine setspecscaly(name, nlength)

      use cdwrk_module
      implicit none

      integer nlength, name(nlength), i, j, maxlen
      double precision ::  val
      parameter (maxlen=256)
      character filet*(maxlen)
      character*(maxspnml) spname, spinname
!      
!     Convert encoded names to strings, and open file
!
      if (nlength.GT.maxlen) then
         call bl_abort('FORT_SETSPECSCAL: scale file name too long')
      end if
      
      do i = 1, nlength
         filet(i:i) = char(name(i))
      end do
      open(unit=51,status='OLD',form='FORMATTED',  &
           file=filet(1:nlength),err=30)
      
 10   continue 
      read(51,*,end=20) spinname, val
      do j = 1,Nspec
         call get_spec_name(spname,j)
         if (spname .eq. spinname) then
            spec_scalY(j) = ABS(val)
         end if
      end do
      goto 10
 20   close(51)
      goto 40
 30   write(6,*) 'Trouble opening file = ',filet(1:nlength)
      call bl_abort(" ")
 40   continue 
      end subroutine setspecscaly

      subroutine initchem()

      use cdwrk_module
      use conp_module

      logical error
      integer ioproc, myid, ierr, n, RTOT, lout, namlen, i
      integer getckspecname
      parameter (ioproc = 0)
      character*(maxspnml) name
      integer coded(maxspnml)

      call EGINICD(eg_nodes, lout, eg_IFLAG, eg_ITLS,  &
           RWRK(egbr), egr, IWRK(egbi), egi)

!      
!     Set pointers in conp common blocks (used in conpF evaluations)
!
      CALL CKINDX(IWRK(ckbi),RWRK(ckbr),Nelt,Nspec,Nreac,Nfit)
      NEQ   = Nspec + 1
      NP    = dvdbr
      NWT   = NP  + 1
      NZ    = NWT + Nspec
      RTOT  = NZ  + NEQ - 1

      if (RTOT .GT. dvder) then
         write(6,*) 'Memory layout bust, dvdr not big enough'
         write(6,*) RTOT, dvder
         call bl_abort(" ")
      end if
!
!     Set molecular weights
!
      CALL CKWT(IWRK(ckbi), RWRK(ckbr), RWRK(NWT))
!
!     Find N2 in the list.
!
      iN2 = -1
      do n = 1,Nspec
         call get_spec_name(name,n)
         if (name .eq. 'N2' ) iN2 = n
      end do
      if (iN2.eq.-1)
     &     write(6,*) '.....warning: no N2 in chemistry species list'
      end subroutine initchem


      SUBROUTINE eginicd (NP, LOUT, IFLAG, ITLS, WEG, LWEG, IWEG, LIWEG)
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
         WRITE(LOUT,'(1X,''Error: the number of species must '', &
                  '' be larger or equal to 2'')')
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
         WRITE(LOUT,'(//1X,''Error: the length of WEG should be '',
     &                   ''at least'',I12//)') INEXT
         ILOW = 1
      ENDIF
      IF ( IINXT .GT. LIWEG ) THEN
         WRITE(LOUT,'(//1X,''Error: the length of IWEG should be '', &  
                ''at least'',I12//)') IINXT
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
      CALL LEVEPS (NS, WEG(IEGEPS), WEG(IEGSIG), WEG(IEGDIP), 
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
         WEG(ICTAIJ+IJ) = WEG(IFITA0  )    + WEG(IFITA0+1)*T1  &
                          + WEG(IFITA0+2)*T2 + WEG(IFITA0+3)*T3  &
                          + WEG(IFITA0+4)*T4 + WEG(IFITA0+5)*T5  &
                          + WEG(IFITA0+6)*T6 
         WEG(ICTBIJ+IJ) = WEG(IFITB0  )    + WEG(IFITB0+1)*T1  &
                          + WEG(IFITB0+2)*T2 + WEG(IFITB0+3)*T3  &
                          + WEG(IFITB0+4)*T4 + WEG(IFITB0+5)*T5  &
                          + WEG(IFITB0+6)*T6 
         WEG(ICTCIJ+IJ) = WEG(IFITC0  )    + WEG(IFITC0+1)*T1  &
                          + WEG(IFITC0+2)*T2 + WEG(IFITC0+3)*T3  &
                          + WEG(IFITC0+4)*T4 + WEG(IFITC0+5)*T5  &
                          + WEG(IFITC0+6)*T6 
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     Evaluate FITA, FITB and FITC 
!-----------------------------------------------------------------------
      CALL EGABC ( NS, NFIT, WEG(IFITA), WEG(IFITB), WEG(IFITC), & 
           WEG(IFITA0), WEG(IFITB0), WEG(IFITC0), & 
           WEG(IEPSIJ) )
!-----------------------------------------------------------------------
      RETURN
      END subroutine eginicd


      integer function getckmaxnamelen()
      use cdwrk_module
      implicit none
      FORT_GETCKMAXNAMELEN = maxspnml
      end function getchmaxnamelen

      subroutine getckdimparams(imaxreac, imaxspec, imaxelts, & 
                             imaxord, imaxthrdb, imaxtp, imaxsp, & 
                             imaxspnml)
      use cdwrk_module
      implicit none
      integer imaxreac, imaxspec, imaxelts, imaxord
      integer imaxthrdb, imaxtp, imaxsp, imaxspnml
      imaxreac = maxreac
      imaxspec = maxspec
      imaxelts = maxelts
      imaxord = 10
      imaxthrdb = maxthrdb
      imaxtp = maxtp
      imaxsp = maxsp
      imaxspnml = maxspnml
      end subroutine getckdimparams

      subroutine findlhs(reactions, Nreacs, id)
      use cdwrk_module
      implicit none
      integer reactions(*), Nreacs, id
      integer j, n, Ndim, Nids, KI(maxsp), NU(maxsp)
! #ifdef MIKE
!     Ndim = maxsp
!     if ((id.le.0).or.(id.gt.Nspec)) then
!        write(6,*) 'FORT_FINDLHS:  species id out of range: ',id
!        call bl_abort(" ")
!     end if
!     Nreacs = 0
!     do j=1,Nreac
!        CALL CKINU(j, Ndim, IWRK(ckbi), RWRK(ckbr), Nids, KI, NU)
!        do n=1,Nids
!           if ((KI(n).eq.id).and.(NU(n).lt.0)) then
!              Nreacs = Nreacs + 1
!              reactions(Nreacs) = j
!           endif
!        end do
!     end do
! #else
      call bl_abort("FORT_FINDLHS not implemented")
! #endif
      end subroutine findlhs

      subroutine findrhs(reactions, Nreacs, id)
      use cdwrk_module
      implicit none
      integer reactions(*), Nreacs, id
      integer j, n, Ndim, Nids, KI(maxsp), NU(maxsp)
!  #ifdef MIKE
!        Ndim = maxsp
!        if ((id.le.0).or.(id.gt.Nspec)) then
!           write(6,*) 'FORT_FINDRHS:  species id out of range: ',id
!           call bl_abort(" ")
!        end if
!        Nreacs = 0
!        do j=1,Nreac
!           CALL CKINU(j, Ndim, IWRK(ckbi), RWRK(ckbr), Nids, KI, NU)
!           do n=1,Nids
!              if ((KI(n).eq.id).and.(NU(n).gt.0)) then
!                 Nreacs = Nreacs + 1
!                 reactions(Nreacs) = j
!              endif
!           end do
!        end do
!  #else
      call bl_abort("FORT_FINDRHS not implemented")
!  #endif
      end subroutine findrhs

      subroutine setnu(nu,lenNU)
      use cdwrk_module
      implicit none
      integer lenNU
      integer nu(maxreac,maxspec)
      integer i

      if (lenNU .lt. maxreac*maxspec) then
         write(6,*) 'FORT_CKNU:  nu work array too small: '
         call bl_abort(" ")
      endif
      call CKNU(maxreac, IWRK(ckbi), RWRK(ckbr), nu)
      end subroutine setnu

      subroutine ckinu(Nids,KI,lenKI,NU,lenNU,rxnID,nuAll)
      use cdwrk_module
      implicit none
      integer lenKI,lenNU,NDIM1
      integer rxnID, Nids, KI(lenKI), NU(lenNU), nuAll(maxreac,maxspec)
      integer Ndim, k
      Ndim = MIN(lenKI,lenNU)
      if ((rxnID.le.0).or.(rxnID.gt.Nreac)) then
         write(6,*) 'FORT_CKINU:  reaction id out of range: ',rxnID
         call bl_abort(" ")
      end if
      if (Ndim.lt.maxsp) then
         call bl_abort('FORT_CKINU:  KI or NU not long enough')
      end if
      Nids = 0
      do k=1,Nspec
         if (nuAll(rxnID,k).ne.0) then
            Nids = Nids + 1
            KI(Nids) = k
            NU(Nids) = nuAll(rxnID,k)
         endif
      enddo
      end subroutine ckinu

      integer function ckeltxinspy(eltID, spID)
      use cdwrk_module
      implicit none
      integer eltID, spID, i, j
      integer NCF(maxelts,maxspec)
      CALL CKNCF(maxelts, IWRK(ckbi), RWRK(ckbr), NCF)
      ckeltxinspy = NCF(eltID+1,spID+1)
      end function ckeltxinspy

      integer function getcknumspec()
      use cdwrk_module
      implicit none
      getcknumspec = Nspec
      end function getcknumspec

      integer function getcknumelt()
      use cdwrk_module
      implicit none
      getcknumelt = Nelt
      end fucntion getcknumelt

      integer function getcknumreac()
      use cdwrk_module
      implicit none
      getcknumreac = Nreac
      end function getchknumreac

      double precision function runiv()
      use cdwrk_module
      implicit none
      double precision Ruc, Pa
      call CKRP(IWRK(ckbi),RWRK(ckbr),FORT_RUNIV,Ruc,Pa)
!     1 erg/(mole.K) = 1.e-4 J/(kmole.K)
      runiv = runiv*1.d-4
      end runiv

      double precision function p1atmmks()
      use cdwrk_module
      implicit none
      double precision Ru, Ruc, Pa
      call CKRP(IWRK(ckbi),RWRK(ckbr),Ru,Ruc,Pa)
!     1 N/(m.m) = 0.1 dyne/(cm.cm)
      p1atmmks = Pa*1.d-1
      end function p1atmmks

      integer function getckeltname(i, coded)
      use cdwrk_module
      implicit none
      integer i
      integer coded(*)
      integer names(0:maxelts*2)
      integer ctr, nlen
      nlen = 2
      call CKSYME(names,nlen)
      coded(1) = names(2*(i-1)  )
      coded(2) = names(2*(i-1)+1)
      if (coded(2).eq.ICHAR(' ')) then
         getckeltname = 1
      else
         getckeltname = 2
      endif
      end function getckeltname

      integer function getckspecname(i, coded)
      use cdwrk_module
      implicit none
      integer i
      integer coded(*)
      integer names(maxspec*maxspnml)
      integer j, str_len
      call CKSYMS(names, maxspnml)
      do j = 1, maxspnml
         coded(j) = names(maxspnml*(i-1)+j)
      end do
      do j = 1, maxspnml
         if (coded(j).eq.ICHAR(' ')) then
            str_len = j
            exit
         endif 
      end do
      getckspecname = str_len - 1
      end function getckspecname

      integer function cksymr(fortReacIdx, coded)
      use cdwrk_module
      implicit none
      integer fortReacIdx
      integer coded(*)
      character*(72) line 
      integer j, str_len, istr, iend, lout
      logical error
      data error /.false./
!  #ifdef MIKE
!     lout = 6
!     call CKSYMR(fortReacIdx,lout,IWRK(ckbi),RWRK(ckbr),
!    &     CWRK(ckbc),str_len,line,error)
!     if (error) then
!        write(lout,*) 'Could not get reaction name for ',fortReacIdx
!        call bl_abort(" ")
!     end if
!  !      
!  !     Encode the name for transfer to C++
!  !
!     istr = 1
!     do while (line(istr:istr) .EQ. ' ')
!        istr = istr + 1
!     end do
!     do j = 0, str_len-1
!        coded(j+1) = ICHAR(line(istr+j:istr+j))
!     end do
!     FORT_CKSYMR = str_len
!  #else
      cksymr = 0
      call bl_abort("FORT_CKSYMR not implemented")
!  #endif
      end function cksymr

      subroutine get_spec_name(name, j)
      use cdwrk_module
      implicit none
      integer i, j, FORT_GETCKSPECNAME
      integer coded(maxspnml), len
      character*(maxspnml) name
      len = FORT_GETCKSPECNAME(j, coded)
      do i = 1, maxspnml
         name(i:i) = ' '
      end do
      do i = 1, len
         name(i:i) = char(coded(i))
      end do
      end subrouitne get_spec_name

      subroutine get_spec_number(name, j)
      use cdwrk_module
      implicit none
      integer j, n
      character*(*) name
      character*(maxspnml) locName
      
      j = -1
      do n = 1, Nspec
         call get_spec_name(locName, n)
         if (locName .EQ. name) j = n
      end do
      end subroutine get_spec_number

      subroutine getckmwt(mwt)
      use cdwrk_module
      implicit none
      REAL_T mwt(*)
      integer n
!     Result in kg/kmole
      call CKWT(IWRK(ckbi),RWRK(ckbr),mwt)
      end subroutine getckmwt

      subroutine getckawt(awt)
      use cdwrk_module
      implicit none
      REAL_T awt(*)
!  #ifdef MIKE
!        integer n
!  !     Result in kg/kmole
!        call CKAWT(IWRK(ckbi),RWRK(ckbr),awt)
!  #else
      call bl_abort("FORT_GETCKAWT not implemented")
!  #endif
      end subroutine getckawt

      subroutine conpFY(N, TIME, Z, ZP, RPAR, IPAR)
      use cdwrk_module
      use conp_module
      implicit none
!   #include "visc.H"
      REAL_T TIME, Z(NEQ), ZP(NEQ), RPAR(*)
      integer N, IPAR(*)
      
      REAL_T RHO, CPB, SUM, H, WDOT, WT, THFAC
      integer K

      REAL_T CONC(maxspec), WDOTS(maxspec), ENTHALPY(maxspec)
!
!     Variables in Z are:  Z(1)   = T
!                          Z(K+1) = Y(K)
      
      CALL CKRHOY(RPAR(NP),Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),RHO)
      CALL CKCPBS(Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),CPB)
      CALL CKYTCP(RPAR(NP),Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),CONC)
!
!     Get net production rates.  Compute such that production from -ve
!     reactants gives zero contrib.
!      
!      do k=1,Nspec
!         if (CONC(k) .lt. zero) write(6,*) '.....negative C',k
!      end do
      
!  this was in an ifdef
!     do k=1,Nspec
!        CONC(k) = MAX(CONC(k),zero)
!     end do 
      
      CALL CKWC(Z(1), CONC, IPAR(ckbi), RPAR(ckbr), WDOTS)
      CALL CKHMS(Z(1), IPAR(ckbi), RPAR(ckbr), ENTHALPY)
!
!     Form governing equation
!
!     THFAC = one / thickFacCH
      SUM = zero
      DO K = 1, Nspec
         H    = ENTHALPY(K)
         WDOT = WDOTS(K) * THFAC
         WT   = RPAR(NWT+K-1)
         ZP(K+1) = WDOT * WT / RHO
         SUM = SUM + H * WDOT * WT
      END DO
      ZP(1) = -SUM / (RHO*CPB)

!  #if 0
!     print*, 'Z:'
!     do k = 1, Nspec+1
!        write(6,996) Z(K)
!     end do
!     print*, 'ZP:'
!     do k = 1, Nspec+1
!        write(6,996) ZP(K)
!     end do
!     print*, 'WDOT:'
!     do k = 1, Nspec
!        write(6,996) WDOTS(K)
!     end do

!  996   format(e30.22)
!  #endif

      END subroutine conpFY

      subroutine conpJY(N, TN, Y, SAVF, NFE, FTEM, ML,
     &                  MU, PD, NRPD, RPAR, IPAR)
      use cdwrk_module
      use conp_module
      implicit none
      REAL_T SAVF
      REAL_T PD, RPAR(*), TN, Y
      dimension SAVF(*)
      integer N, NRPD, ML, MU, IPAR(*)
      dimension Y(N), PD(NRPD,N)
      REAL_T FTEM
      dimension FTEM(*)
      integer NFE
      call bl_abort("conpJY: SHOULD NOT BE HERE!")
      END subroutine conpJY

      integer function TfromeYpt(T,ein,Y,errMax,NiterMAX,res)
      use cdwrk_module
      implicit none
      integer NiterMAX,Niter,n,NiterDAMP
      REAL_T T,ein,Y(*),errMAX,res(0:NiterMAX-1)
      REAL_T TMIN,TMAX,e

      parameter (TMIN=200, TMAX=5000)
      REAL_T  T0,h,cp,cv,de,temp,RoverWbar,Wbar,RU,RUC,P1ATM
      REAL_T dT, etarg
      logical out_of_bounds, converged, soln_bad, stalled
      REAL_T e300,cv300,e6500,cv6500
      integer ihitlo,ihithi

      out_of_bounds(temp) = (temp.lt.TMIN) .or. (temp.gt.TMAX)

      NiterDAMP = NiterMAX
      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
         T0 = T
      else
         T0 = half*(TMIN+TMAX)
         T = T0
      end if
      Niter = 0
      de = zero
      soln_bad = .FALSE.
      etarg = ein * BL_REAL_E(1.0,4)
      ihitlo = 0
      ihithi = 0

      CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e)

      de = two*ABS(e - etarg)/(one + ABS(e) + ABS(etarg))
      res(Niter) = de
      converged = de.le.errMAX

      do while ((.not.converged) .and. (.not.soln_bad))
         CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv)
         dT = (etarg - e)/cv
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
            TfromeYpt = -1
            goto 100
         else
            CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e)
            de = two*ABS(e - etarg)/(one + ABS(e) + ABS(etarg))
            res(Niter) = de
            Niter = Niter + 1
         end if
         if (Niter .ge. NiterMAX) then
            TfromeYpt = -2
            goto 100
         endif
         converged = (de.le.errMAX) .or. (ABS(dT).le.errMAX)

         if((ihitlo.eq.1).and.(e.gt.etarg))then
            T = 300.d0
            CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e300)
            CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv300)
            T=300.d0+(etarg-e300)/cv300
            converged = .true.
         endif
         if((ihithi.eq.1).and.(e.lt.etarg))then
            T = 6500.d0
            CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e6500)
            CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv6500)
            T=6500.d0+(etarg-e6500)/cv6500
            converged = .true.
         endif

      end do

!     Set max iters taken during this solve, and exit
      TfromeYpt = Niter
      return

!     Error condition....dump state and bail out
 100  continue

      write(6,997) 'T from (e,Y): failed'
      write(6,997) 'iterations tried = ',Niter
      write(6,998) 'initial T = ',T0
      write(6,998) 'current T = ',T
      write(6,998) 'species mass fracs:'
      do n = 1,Nspec
         write(6,998) '  ',Y(n)
      end do
      write(6,998)
      write(6,998) 'residual = e - h + RT/Wbar [cgs]'
      do n = 0,Niter-1
         write(6,998) '  ',res(n)
      end do

 997  format(a,3(i4,a))
 998  format(a,d21.12)

      end function TfromeYpt

      subroutine FORT_TfromHYpt(T,Hin,Y,errMax,NiterMAX,res,Niter)
      use cdwrk_module
      implicit none
      REAL_T T,Y(*),H,Hin
      REAL_T TMIN,TMAX,errMAX
      integer NiterMAX,Niter,n,NiterDAMP, Discont_NiterMAX
      parameter (TMIN=250, TMAX=5000, Discont_NiterMAX=100)
      REAL_T  T0,cp,cv,dH,temp,RoverWbar,Wbar,RU,RUC,P1ATM
      REAL_T res(0:NiterMAX-1),dT, Htarg
      logical out_of_bounds, converged, soln_bad, stalled, discont
      REAL_T HMIN,cpMIN,HMAX,cpMAX
      integer ihitlo,ihithi,j
      REAL_T old_T, old_H, Tsec, Hsec

      out_of_bounds(temp) = (temp.lt.TMIN-one) .or. (temp.gt.TMAX)

      NiterDAMP = NiterMAX
      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
         T0 = T
      else
         T0 = half*(TMIN+TMAX)
         T = T0
      end if
      Niter = 0
      soln_bad = .FALSE.
      Htarg = Hin * 1.d4
      ihitlo = 0
      ihithi = 0

      CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),H)

      old_T = T
      old_H = H

      dH = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
      res(Niter) = dH
      converged = dH.le.errMAX
      stalled = .false.

      do while ((.not.converged) .and.  (.not.stalled) .and.(.not.soln_bad))

         CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cp)
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
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),H)
            dH = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
            res(Niter) = min(dH,abs(dT))
            Niter = Niter + 1
         end if
         converged = (dH.le.errMAX) .or. (ABS(dT).le.errMAX)
         if (Niter .ge. NiterMAX) then
            if(abs(T-1000.d0).le.1.d-3.and. dH.le.1.d-5)then
              converged = .true.
            else  
              write(6,*)" missing my test",T,dH
              Niter = -2
              goto 100
            endif
         endif

         if ((ihitlo.eq.1).and.(H.gt.Htarg)) then
            T = TMIN
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),HMIN)
            CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cpMIN)
            T=TMIN+(Htarg-HMIN)/cpMIN
            converged = .true.
         endif
         if ((ihithi.eq.1).and.(H.lt.Htarg)) then
            T = TMAX
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),HMAX)
            CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cpMAX)
            T=TMAX+(Htarg-HMAX)/cpMAX
            converged = .true.
         endif

!     If the iterations are failing, perhaps it is because the fits are discontinuous
!     The following implements a modified secant method to hone in on a discontinity in h
!     with T.  Discont_NiterMAX is fairly large because this process can be particularly
!     slow to converge if the Htarg value happens to lay between the discontinuous function
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
               CALL CKHBMS(Tsec,Y,IWRK(ckbi),RWRK(ckbr),Hsec)
               if ( (Hsec-Htarg)*(Htarg-H) .gt. 0.d0 ) then
                  old_H = H
                  old_T = T
               endif
               H = Hsec
               T = Tsec
               stalled = (2*ABS(old_T-T)/(old_T+T).le.errMAX)
!               stalled = (ABS(dT).le.errMAX)
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
      end function TfromHYpt
  
      integer function open_vode_failure_file ()
      implicit none
      character*30 name, myproc
      integer lout,i,j,k,idx

!     Hardwire the unit number to 26 for the moment
      lout = 26 
      call bl_pd_myproc(i)
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
      end function open_vode_failure_file
    
      block data tranjunk
      use cdwrk_module
      use conp_module
      
      data verbose_vode / 0 /
      data max_vode_subcycles / 15000 /
      data spec_scalY / maxspec*one /
      data thickFacCH / 1.d0 /
      data nstiff / 1 /

      end block data tranjunk
