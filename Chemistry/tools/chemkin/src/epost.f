C CVS $Revision: 1.1.1.1 $ created $Date: 2006/05/26 19:09:33 $
C///////////////////////////////////////////////////////////////////////
C
C  VERSION 3.6, 97/12/10
C
C  post processor for AURORA - reads user keywords and queries 
C                              solution files for desired data.
C                              creates x-y data files for solution
C                              results obtained by running 
C                              continuations in AURORA, and varying
C                              a single parameter.
C
C See header for subroutine RDKEY for information on user keywords
C
C///////////////////////////////////////////////////////////////////////
C Changes history:
C Changes for 3.6 (97/12/10, E. Meeks)
C 1. Add capability to printout sensitivity coefficients to separate 
C    file for each solution found in save.bin.
C 2. Add logical variable LSOLSN, passed into RDSAVE, to be set when
C    sensitivity data is found and read correctly.  If not set, ignore
C    any requests for sensitivity data printout.
C Changes for 3.5 (97/08/06, E. Meeks)
C 1. Fix problems with multiple PSRs: Only increment NS when IPSR=1 in
C    RDSAVE.
C Changes for 3.4 (97/08/05, E. Meeks)
C 1. Add reading of parameter BIASPW from aurora save file to 
C    correspond to aurora revision 6.9.
C 2. Add print option to write parameter BIASPW values to output.
C 3. Add print option for species outlet flow rates, using keywrd OFRS.
C 4. Add declaration for VERSN as CHARACTER*80 to fix garbage printing
C    in RDSAVE.
C 5. Fix bug#066:  Added logicals LCHEM and LSOL in RDSAVE to check 
C    whether solution is read and gracefully exit if not. 
C Changes for 3.3 (97/07/25, E. Meeks)
C 1. Fix bug #012: remove extra spaces in logical operators on line 689.
C 2. Fix lines #1,1835, and 1932 longer than 72 characters.
C Changes for 3.2 (97/03/01, F. Rupley)
C 1. make new main "driver" program to set up arrays,
C    to satisfy Reaction Design requirement of providing object
C    files instead of source code - move OPENs to driver;
C    RDKEY needs upper-case OPTION string instead of making
C    all arguments upper-case (for the case where species names may
C    be mixed-case)
C 2. fix wrong dimension for ROPS in SUBROUTINE RDSAVE
C 3. do not read SDOTI array in RDSAVE if IISUR(IM).LE.0
C///////////////////////////////////////////////////////////////////////
C
      SUBROUTINE PPSR (LIN, LOUT, LSAVE, LPRNT, MXSPEC, MXREAC, MAXMAT,
     1                 MAXPSR, KBMAX, MXPTS, MXPHSE, MAXION, MAXPAR,
     2                 LREAC, NCMAX,
     2                 REACTS,
     2                 LENICK, ICKWRK, LENRCK, RCKWRK, LENCCK, CCKWRK,
     3                 LENISK, ISKWRK, LENRSK, RSKWRK, LENCSK, CSKWRK,
     4                 LI, I, LR, R, LC, C, LL, L)
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
      DIMENSION ICKWRK(LENICK), RCKWRK(LENRCK),
     1          ISKWRK(LENISK), RSKWRK(LENRSK),
     2          I(LI), R(LR)
      CHARACTER*16 CCKWRK(LENCCK), CSKWRK(LENCSK), C(LC)
      CHARACTER*60 REACTS(LREAC)
      LOGICAL KERR, L(LL)
      DATA KERR/.FALSE./
C
      IKION = 1
      IKKPHA= IKION + MAXION
      IKPSP = IKKPHA+ MXPHSE * MAXMAT
      IKSPC = IKPSP + MAXPAR
      IKRSP = IKSPC + MXSPEC
      IMPSP = IKRSP + MXSPEC
      IKTFL = IMPSP + MAXPAR
      IKCHG = IKTFL + MXSPEC
      IGPHS = IKCHG + MXSPEC
      IKFRST= IGPHS + MXPHSE
      IKLAST= IKFRST+ MXPHSE * MAXMAT
      IKSPN = IKLAST+ MXPHSE * MAXMAT
      IKSPF = IKSPN + MXSPEC
      IKCUR = IKSPF + MXSPEC
      IMISK = IKCUR + MAXION
      IMRSK = IMISK + MAXMAT
      IMCSK = IMRSK + MAXMAT
      IKKSUR= IMCSK + MAXMAT
      IKKBUL= IKKSUR+ MAXMAT
      IKKM  = IKKBUL+ MAXMAT
      INPHA = IKKM  + MAXMAT
      INNSUR= INPHA + MAXMAT
      INFSUR= INNSUR+ MAXMAT
      INLSUR= INFSUR+ MAXMAT
      INNBUL= INLSUR+ MAXMAT
      INFBUL= INNBUL+ MAXMAT
      INLBUL= INFBUL+ MAXMAT
      ISUR  = INLBUL+ MAXMAT
      INSPCH= ISUR  + MAXMAT
      IMAPPH= INSPCH+ MAXMAT
      IKYLD = IMAPPH+ MXPHSE * MAXMAT
      ITOT  = IKYLD + MAXION
      IF (ITOT .GT. LI) THEN
         WRITE (LOUT, *)
     1   'PPSR ERROR: IPAR must be at least ', ITOT
         KERR = .TRUE.
      ENDIF
C
      NAFRC = 1
      NOUTR = NAFRC + MAXMAT
      NS    = NOUTR + NCMAX * MXPTS
      NTOSC = NS    + MXSPEC
      NEQUIV= NTOSC + MXPTS
      NPRES = NEQUIV+ MXPTS * MAXPSR
      NTAU  = NPRES + MXPTS * MAXPSR
      NFLRT = NTAU  + MXPTS * MAXPSR
      NVOL  = NFLRT + MXPTS * MAXPSR
      NAREA = NVOL  + MXPTS * MAXPSR
      NQL   = NAREA + MXPTS * MAXPSR
      NTSURF= NQL   + MXPTS * MAXPSR
      NHTRN = NTSURF+ MXPTS * MAXMAT * MAXPSR
      NTAMB = NHTRN + MXPTS * MAXPSR
      NTIN  = NTAMB + MXPTS * MAXPSR
      NTEIN = NTIN  + MXPTS * MAXPSR
      NXIN  = NTEIN + MXPTS * MAXPSR
      NPOWR = NXIN  + MXPTS * MXSPEC * MAXPSR
      NTGAS = NPOWR + MXPTS * MAXPSR
      NTELEC= NTGAS + MXPTS * MAXPSR
      NXT   = NTELEC+ MXPTS * MAXPSR 
      NXMF  = NXT   + MXSPEC
      NSDEN = NXMF  + MXPTS * MXSPEC * MAXPSR
      NSDEN0= NSDEN + MXPTS * MXPHSE * MAXMAT * MAXPSR
      NWT   = NSDEN0+ MXPHSE * MAXMAT
      NROPG = NWT   + MXSPEC * MAXMAT
      NROPS = NROPG + MXPTS * MXREAC * MXSPEC
      NCONC = NROPS + MXPTS * MXREAC * MXSPEC * MAXMAT
      NSDOTM= NCONC + MXPTS * MXSPEC
      NSDOT = NSDOTM+ MXPTS * MXSPEC * MAXMAT
      NDEN  = NSDOT + MXSPEC
      NGRAT = NDEN  + MXSPEC * MAXMAT
      NTION = NGRAT + MXPTS * MXPHSE * MAXMAT 
      NEION = NTION + MXPTS
      NCURR = NEION + MXPTS * MAXMAT
      NSDOTI= NCURR + MXPTS * MAXION * MAXMAT
      NSITDT= NSDOTI+ MXSPEC * MXREAC * MAXMAT
      NYLD  = NSITDT+ MXPHSE * MAXMAT
      NYTOT = NYLD  + MXPTS * MAXION * KBMAX
      NBHMX = NYTOT + MXPTS * KBMAX
      NTIME = NBHMX + MAXPSR
      NBPWR = NTIME + MXPTS
      NFLROU= NBPWR + MXPTS * MAXMAT
      NFLRS = NFLROU+ MXPTS * MAXPSR
      NSENS = NFLRS + MXPTS * MXSPEC
      NSNSS = NSENS + MXPTS * MXREAC * MXSPEC
      NTOT  = NSNSS + MXPTS * MXREAC * MXSPEC * MAXMAT
      IF (NTOT .GT. LR) THEN
         WRITE (LOUT, *)
     1   'PPSR ERROR: RPAR must be at least ', NTOT
         KERR = .TRUE.
      ENDIF
C
      NCKNAM = 1
      NCPNAM = NCKNAM + MXSPEC * MAXMAT
      NCPSP  = NCPNAM + MXPHSE * MAXMAT
      NCMNAM = NCPSP  + MAXPAR
      NCLABL = NCMNAM + MAXMAT
      NCLABL2= NCLABL + NCMAX
      NCTOT  = NCLABL2+ NCMAX
      IF (NCTOT .GT. LC) THEN
         WRITE (LOUT, *)
     1   'PPSR ERROR: CPAR must be at least ', NCTOT
         KERR = .TRUE.
      ENDIF
C
      NRGAS   = 1
      NRSUR   = NRGAS + MXREAC
      NRTOT   = NRSUR + MXREAC * MAXMAT
      IF (NRTOT .GT. LREAC) THEN
         WRITE (LOUT, *)
     1   'PPSR ERROR: REACTS must be at least ', NRTOT
         KERR = .TRUE.
      ENDIF
C
      LTOT = MXPTS * MAXPSR
      IF (LTOT .GT. LL) THEN
         WRITE (LOUT, *)
     1   'PPSR ERROR: LPAR must be at least ', LTOT
         KERR = .TRUE.
      ENDIF
      IF (KERR) RETURN
C
      CALL POSTPSR (LIN, LOUT, LSAVE, LPRNT, MXSPEC, MXREAC, MAXMAT,
     1              MAXPSR, KBMAX, MXPTS, MXPHSE, MAXION, MAXPAR,
     2              NCMAX,
     2              LENICK, ICKWRK, LENRCK, RCKWRK, LENCCK, CCKWRK,
     3              LENISK, ISKWRK, LENRSK, RSKWRK, LENCSK, CSKWRK,
     4I(IKION),I(IKKPHA),I(IKPSP),I(IKSPC),I(IKRSP),I(IMPSP),I(IKTFL),
     5I(IKCHG),I(IGPHS),I(IKFRST),I(IKLAST),I(IKSPN),I(IKSPF),
     6I(IKCUR), I(IMISK), I(IMRSK), I(IMCSK), I(IKKSUR), I(IKKBUL),
     7I(IKKM), I(INPHA), I(INNSUR), I(INFSUR), I(INLSUR), 
     8I(INNBUL), I(INFBUL), I(INLBUL), I(ISUR), I(INSPCH), I(IMAPPH),
     9I(IKYLD),
     *R(NAFRC), R(NOUTR), R(NS), R(NTOSC), R(NEQUIV),
     1R(NPRES), R(NTAU), R(NFLRT), R(NVOL), R(NAREA),
     2R(NQL), R(NTSURF), R(NHTRN), R(NTAMB), R(NTIN),
     3R(NTEIN), R(NXIN), R(NPOWR), R(NTGAS), R(NTELEC),
     4R(NXT), R(NXMF), R(NSDEN), R(NSDEN0), R(NWT),
     5R(NROPG), R(NROPS), R(NCONC), R(NSDOTM), R(NSDOT),
     6R(NDEN), R(NGRAT), R(NTION), R(NEION), R(NCURR),
     7R(NSDOTI), R(NSITDT), R(NYLD), R(NYTOT), R(NBHMX),
     8R(NTIME), R(NBPWR), R(NFLROU), R(NFLRS), R(NSENS), R(NSNSS),
     9C(NCKNAM), C(NCPNAM), C(NCPSP), C(NCMNAM), C(NCLABL),
     *C(NCLABL2),
     1REACTS(NRGAS), REACTS(NRSUR), L, KERR)
      RETURN
     
      END
C
C----------------------------------------------------------------------C
      SUBROUTINE POSTPSR (LIN, LOUT, LSAVE, LPRNT, KMAX, IMAX, 
     1           MAXMAT, MAXPSR, KBMAX, JMAX, NMAX, KIMAX, MAXPAR,
     2           NCMAX,
     2           LENICK, ICKWRK, LENRCK, RCKWRK, LENCCK, CCKWRK,
     3           LENISK, ISKWRK, LENRSK, RSKWRK, LENCSK, CSKWRK,
     4           KION, KKPHAS, KPSPEC, KSPEC, KRSPEC, MPSPEC,
     5           KTFL, KCHG, KGPHAS, KFIRST, KLAST, KSPN, KSPF, KCURR,
     5           IMISK, IMRSK, IMCSK, KKSURF, KKBULK, KKM, NPHASE,
     7           NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     8           IISUR, NSPHCH, MAPPH, KYLD,
     *           AFRAC, OUTRAY, S, TOSCCM, EQUIV, PRES, TAU, FLRT,
     1           VOL, AREA, QL, TSURF, HTRN, TAMB, TIN, TEIN, XIN,
     2           POWR, TGAS, TELEC, XT, XMF, SDEN, SDEN0, WT,
     3           ROPG, ROPS, CONC, SDOTM, SDOT, DEN, GRATE, TIONP,
     4           EIONSH, CURR, SDOTI, SITDTI, YLD, YLDTOT, BHMXI,
     5           TIME, BIASPW, FLROUT, FLRS, SENS, SENSS,
     6           KSYM, PSYM, PSPEC, MSYM, LABEL, LABEL2, 
     7           IGAS, ISURF,
     8           LHTRN, KERR)
C
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
      PARAMETER(BOLTZ=1.380658D-16, AVOG=6.0221367D23, 
     &          ECHG=1.60217733D-12, EV2K=ECHG/BOLTZ)
      PARAMETER (NT=1, NTE=2, NYS=2, NY1=3, TOWATT=1.E-7,
     1           TOTORR=7.5006D-4 , TOEV=1./EV2K, SMALL=1.E-25)
C
      DIMENSION ICKWRK(LENICK), ISKWRK(LENISK), KION(KIMAX),
     1          KKPHAS(NMAX,MAXMAT), KPSPEC(MAXPAR), KSPEC(KMAX), 
     2          KRSPEC(KMAX), MPSPEC(MAXPAR),
     3          KTFL(KMAX),  KCHG(KMAX), KGPHAS(NMAX), 
     3          KFIRST(NMAX,MAXMAT), KLAST(NMAX,MAXMAT),
     4          KSPN(KMAX), KSPF(KMAX), KCURR(KIMAX), IMISK(MAXMAT),
     5          IMRSK(MAXMAT), IMCSK(MAXMAT), KKSURF(MAXMAT),
     6          KKBULK(MAXMAT), KKM(MAXMAT), NPHASE(MAXMAT), 
     7          NNSURF(MAXMAT), NFSURF(MAXMAT), NLSURF(MAXMAT), 
     8          NNBULK(MAXMAT), NFBULK(MAXMAT), NLBULK(MAXMAT),
     9          IISUR(MAXMAT), NSPHCH(MAXMAT), AFRAC(MAXMAT),
     *          RCKWRK(LENRCK), RSKWRK(LENRSK), OUTRAY(NCMAX,JMAX),
     1          S(KMAX), TEMP(3), TOSCCM(JMAX), EQUIV(JMAX,MAXPSR), 
     3          PRES(JMAX,MAXPSR), TAU(JMAX,MAXPSR), FLRT(JMAX,MAXPSR),
     4          VOL(JMAX,MAXPSR), AREA(JMAX,MAXPSR), QL(JMAX,MAXPSR),
     5          TSURF(JMAX,MAXMAT,MAXPSR), HTRN(JMAX,MAXPSR), 
     6          TAMB(JMAX,MAXPSR), TIN(JMAX,MAXPSR), TEIN(JMAX,MAXPSR),
     7          XIN(JMAX,KMAX,MAXPSR), POWR(JMAX,MAXPSR),
     8          TGAS(JMAX,MAXPSR), TELEC(JMAX,MAXPSR), XT(KMAX),
     9          XMF(JMAX,KMAX,MAXPSR), SDEN(JMAX,NMAX,MAXMAT,MAXPSR),
     *          SDEN0(NMAX,MAXMAT),WT(KMAX,MAXMAT),ROPG(JMAX,IMAX,KMAX),
     1          ROPS(JMAX,IMAX,KMAX,MAXMAT), CONC(KMAX,JMAX), 
     2          SDOTM(JMAX,KMAX,MAXMAT),SDOT(KMAX), DEN(KMAX,MAXMAT), 
     3          GRATE(JMAX,NMAX,MAXMAT), TIONP(JMAX), 
     4          EIONSH(JMAX,MAXMAT),CURR(JMAX,KIMAX,MAXMAT), 
     5          SDOTI(KMAX,IMAX,MAXMAT), SITDTI(NMAX,MAXMAT),
     6          MAPPH(NMAX,MAXMAT),YLD(JMAX,KIMAX,KBMAX), KYLD(KIMAX),
     7          YLDTOT(JMAX,KBMAX), BHMXI(MAXPSR), TIME(JMAX), 
     8          BIASPW(JMAX,MAXMAT),FLROUT(JMAX,MAXPSR),FLRS(KMAX,JMAX),
     9          SENS(JMAX,IMAX,KMAX+2),SENSS(JMAX,IMAX,KMAX+2,MAXMAT)
C
      CHARACTER*16 KSYM(KMAX,MAXMAT), CCKWRK(LENCCK), CSKWRK(LENCSK),
     1             PSYM(NMAX,MAXMAT), CXAXIS, PSPEC(MAXPAR), GRUNIT, 
     2             LABEL(NCMAX), BLNK, LABEL2(NCMAX), MSYM(MAXMAT)
      CHARACTER*60 IGAS(IMAX), FNAME, FILENM, LABL, ISURF(IMAX,MAXMAT),
     1             FTUBE, TUBENM, FSENS, SENSNM
      CHARACTER*3  INUM
      CHARACTER*4  ISNUM
      CHARACTER*6  RLAB
      CHARACTER*10 CHNUM 
      CHARACTER*6  UNITS
C
      LOGICAL LCONT, LHTRN(JMAX,MAXPSR), LGRALL, LCURALL, LCURTOT,
     1        LYLDALL, LYLDTOT, LTUBE, LSENS, KERR, LSOLSN
      EXTERNAL ILASCH
C
      DATA LCONT/.FALSE./, FNAME/' '/, FILENM/' '/
      DATA TUBENM/' '/, LTUBE/.FALSE./, FTUBE/' '/
      DATA LSENS/.FALSE./, SENSNM/' '/, FSENS/' '/, RLAB/' '/
      DATA CHNUM/'0123456789'/
C
      DO 5 K = 1, NMAX
         DO 3 IM = 1, MAXMAT
            PSYM(K,IM) = ' '
 3       CONTINUE
 5    CONTINUE
      DO 10 K = 1, KMAX
         KSYM(K,1) = ' '
 10   CONTINUE
      DO 15 I = 1, 16
         BLNK(I:I) = ' '
 15   CONTINUE
C
C      OPEN (LSAVE,  STATUS='OLD', FORM='UNFORMATTED',
C     1       FILE='save.bin')
C
 20   CONTINUE
      DO 22 L = 1, NCMAX
         LABEL(L) = ' '
         LABEL2(L) = ' '
 22   CONTINUE
      LABL = ' '
C
      FILENM = ' '
      TUBENM = ' '
      SENSNM = ' '
      KERR = .FALSE.
      CALL RDSAVE (LSAVE, LOUT, LENICK, ICKWRK, LENRCK, RCKWRK, LENCCK,
     1             CCKWRK, LENISK, ISKWRK, LENRSK, RSKWRK, LENCSK,
     2             CSKWRK, NT, NTE, NYS, NY1,
     3             KMAX, JMAX, IMAX, II, KSYM, IGAS, ISURF, KKGAS,
     4             KKSURF, KKBULK, IISUR, S, EQUIV, PRES, TAU, FLRT, 
     5             VOL, AREA, QL, TSURF, HTRN, TAMB, LHTRN, TIN, XIN,
     6             TEIN, TGAS, TELEC, POWR, XMF, XT, SDEN, 
     7             MAPPH, KEL, KKION, KION, KKM,
     8             MAXPSR, NMAX, NS, SDEN0, WT, ROPG, ROPS,
     9             PSYM, KKPHAS, NPHASE, NNBULK, NFBULK,
     *             NLBULK, NSPHCH, NUMPSR, KFIRST, KLAST, KCHG,
     1             NNSURF, NFSURF, NLSURF, TIONP, EIONSH, MAXMAT,
     2             NMAT, IMISK, IMRSK, IMCSK, KKTOT, AFRAC, KIMAX,
     3             BHMXI, TIME, MSYM, SDOTI, SITDTI, SDOTM, KERR,
     4             BIASPW, SENS, SENSS, LSOLSN)
      IF (KERR) RETURN
C
      NSDEN1 = NY1 + KKGAS
      NSDENS = NYS + KKGAS
      NPTS = NS
C
      CALL RDKEY (LIN, LOUT, LCONT, FNAME, IPSR, CXAXIS, KXAXIS,
     1            PSPEC, IPSPEC, KPSPEC, KKM, KSPEC, IKSPEC,
     2            KGPHAS, IGPHAS, GRUNIT, LGRALL, KKGAS,
     3            KSYM, PSYM, MSYM, NPHASE, KKPHAS, KRSPEC, IRSPEC,
     4            MPSPEC, MXAXIS, KMAX, NNBULK,
     5            NFBULK, NLBULK, NSPHCH, NUMPSR,
     6            KFIRST, KLAST, MAXPAR, MAPPH, KSPN, IKSPN, KSPF,
     7            IKSPF, NFSURF, 
     8            NLSURF, KEL, KKION, KION, LCURALL, KCURR, ICURR,
     9            LCURTOT, NMAT, IMISK, IMRSK, IMCSK, KKTOT,NMAX,
     *            KKSURF, KKBULK, NNSURF, LYLDALL, KYLD, IYLD, LYLDTOT,
     1            LTUBE, FTUBE, LSENS, FSENS, LSOLSN, KERR)
      IF (KERR) RETURN
C
      WRITE (LOUT, *) ' GAS TEMPERATURE = ', TGAS(1,IPSR)
C
      CALL CKRP (ICKWRK, RCKWRK, RU, RUC, PATM)
      
      DO 46 K = 1, KKGAS
         KTFL(K) = 1
 46   CONTINUE
      IF (KEL.NE.0) THEN
         KTFL(KEL) = 2
         DO 47 KI = 1, KKION
            KTFL(KION(KI)) = 3
 47      CONTINUE
      ENDIF
      CALL CKKTFL (ICKWRK, KTFL)
      DO 48 IMAT = 1, NMAT
         CALL SKKTFL (ISKWRK(IMISK(IMAT)), KTFL)
 48   CONTINUE
C
      DO 60 J = 1, NPTS
         TEMP(1) = TIN(J,IPSR)
         TEMP(2) = TEIN(J,IPSR)
         TEMP(3) = TIONP(J)
C 
C calculate the mass loss to the walls for determining outlet flow
C
         XMSLSS = 0.0
         DO 55 K = 1, KKGAS
            DO 52 IM = 1, NMAT
               XMSLSS = XMSLSS - SDOTM(J,K,IM)*WT(K,1)
 52         CONTINUE
            XT(K) = XIN(J,K,IPSR)
55       CONTINUE
         XMSLSS = XMSLSS * AREA(J,IPSR)
         FLROUT(J,IPSR) = FLRT(J,IPSR) - XMSLSS
C
         CALL CKRHOX (PRES(J,IPSR),TEMP,XT,ICKWRK,
     1                RCKWRK,RHO)
         VOLFIN = FLRT(J,IPSR) / RHO
         TOSCCM(J) =  (298.15/TIN(J,IPSR)) 
     1                      * (PRES(J,IPSR)/PATM) * 60. / RHO
 60   CONTINUE
C
C  check to make sure NCMAX is large enough
C
      IISTOT = 0
      NBSUM = 0
      DO 63 IM = 1, NMAT
         NBSUM = NBSUM + NNBULK(IM)
         IISTOT = IISTOT + IISUR(IM)
 63   CONTINUE
      MAXOUT = IPSPEC+IKSPEC+IKSPN+IKSPF+IGPHAS+NMAT*KKION+ICURR*NMAT
     1         +NBSUM+NBSUM*IYLD+IRSPEC*II+IRSPEC*IISTOT+2
      IF (NCMAX .LT. MAXOUT ) THEN
         WRITE (LOUT, *) ' Error...data output dimension too small,',
     1                   ' NCMAX should be at least ',MAXOUT
         RETURN
      ENDIF
C
C  fill output array, starting with x-axis identification
C
      NCOL = 1
      IF (CXAXIS .EQ. 'SOLUTION_#') THEN
         DO 65 J = 1, NPTS
            OUTRAY(1,J) = FLOAT(J)
 65      CONTINUE
      ELSEIF (CXAXIS .EQ. 'EQUIV_RATIO') THEN
         DO 70 J = 1, NPTS
            OUTRAY(1,J) = EQUIV(J,IPSR)
 70      CONTINUE
      ELSEIF (CXAXIS .EQ. 'PRESSURE(TORR)') THEN
         DO 75 J = 1, NPTS
            OUTRAY(1,J) = PRES(J,IPSR)*TOTORR
 75      CONTINUE
      ELSEIF (CXAXIS .EQ. 'TAU(1/S)') THEN
         DO 80 J = 1, NPTS
            OUTRAY(1,J) = TAU(J,IPSR)
 80      CONTINUE
      ELSEIF (CXAXIS .EQ. 'FLRT(SCCM)') THEN
         DO 85 J = 1, NPTS
            OUTRAY(1,J) = FLRT(J,IPSR)*TOSCCM(J)
 85      CONTINUE
      ELSEIF (CXAXIS .EQ. 'AREA(CM2)') THEN
         DO 90 J = 1, NPTS
            OUTRAY(1,J) = AREA(J,IPSR)
 90      CONTINUE
      ELSEIF (CXAXIS .EQ. 'VOLUME(CM3)') THEN
         DO 95 J = 1, NPTS
            OUTRAY(1,J) = VOL(J,IPSR)
 95      CONTINUE
      ELSEIF (CXAXIS .EQ. 'QL(WATTS)') THEN
         DO 100 J = 1, NPTS
            IF (LHTRN(J,IPSR)) THEN
               OUTRAY(1,J) = HTRN(J,IPSR)*(TGAS(J,IPSR) - TAMB(J,IPSR))
     &                       *TOWATT
            ELSE
               OUTRAY(1,J) = QL(J,IPSR)*TOWATT
            ENDIF
 100     CONTINUE
      ELSEIF (CXAXIS .EQ. 'TSURF') THEN
         DO 105 J = 1, NPTS
            OUTRAY(1,J) = TSURF(J,MXAXIS,IPSR)
 105     CONTINUE
      ELSEIF (CXAXIS .EQ. '1000/TSURF') THEN
         DO 110 J = 1, NPTS
            OUTRAY(1,J) = 1000./TSURF(J,MXAXIS,IPSR)
 110     CONTINUE
      ELSEIF (CXAXIS .EQ. 'HTRN') THEN
         DO 115 J = 1, NPTS
            OUTRAY(1,J) = HTRN(J,IPSR)
 115     CONTINUE
      ELSEIF (CXAXIS .EQ. 'TAMB') THEN
         DO 120 J = 1, NPTS
            OUTRAY(1,J) = TAMB(J,IPSR)
 120     CONTINUE
      ELSEIF (CXAXIS .EQ. 'TIN') THEN
         DO 125 J = 1, NPTS
            OUTRAY(1,J) = TIN(J,IPSR)
 125     CONTINUE
      ELSEIF (CXAXIS .EQ. 'TGAS') THEN
         DO 130 J = 1, NPTS
            OUTRAY(1,J) = TGAS(J,IPSR)
 130     CONTINUE
      ELSEIF (CXAXIS .EQ. '1000/TGAS') THEN
         DO 135 J = 1, NPTS
            OUTRAY(1,J) = 1000./TGAS(J,IPSR)
 135     CONTINUE
      ELSEIF (CXAXIS .EQ. 'TELEC') THEN
         DO 140 J = 1, NPTS
            OUTRAY(1,J) = TELEC(J,IPSR)
 140     CONTINUE
      ELSEIF (CXAXIS .EQ. 'TELEC(EV)') THEN
         DO 145 J = 1, NPTS
            OUTRAY(1,J) = TELEC(J,IPSR)*TOEV
 145     CONTINUE
      ELSEIF (CXAXIS .EQ. 'POWER(WATTS)') THEN
         DO 150 J = 1, NPTS
            OUTRAY(1,J) = POWR(J,IPSR)*TOWATT
 150     CONTINUE
      ELSEIF (CXAXIS .EQ. 'IONE(eV)') THEN
         DO 155 J = 1, NPTS
            OUTRAY(1,J) = EIONSH(J,MXAXIS)*TOEV
 155     CONTINUE
         CXAXIS(9:9) = '_'
C         CXAXIS(10:16) = MSYM(MXAXIS)(1:7)
         CXAXIS(10:10) = 'M'
         CXAXIS(11:11) = CHNUM(MXAXIS+1:MXAXIS+1)
      ELSEIF (CXAXIS .EQ. 'BIAS(W)') THEN
         DO 157 J = 1, NPTS
            OUTRAY(1,J) = BIASPW(J,MXAXIS)/1.e7
 157     CONTINUE
         CXAXIS(8:8) = '_'
C         CXAXIS(9:16) = MSYM(MXAXIS)(1:7)
         CXAXIS(9:9) = 'M'
         CXAXIS(10:10) = CHNUM(MXAXIS+1:MXAXIS+1)
      ELSEIF (CXAXIS .EQ. 'TION') THEN
         DO 165 J = 1, NPTS
            OUTRAY(1,J) = TIONP(J)
 165     CONTINUE
      ELSEIF (CXAXIS .EQ. 'TION(EV)') THEN
         DO 170 J = 1, NPTS
            OUTRAY(1,J) = TIONP(J) *TOEV
 170     CONTINUE
      ELSEIF (CXAXIS .EQ. 'XIN') THEN
         CXAXIS(4:4) = '_'
         CXAXIS(5:16) = KSYM(KXAXIS,1)(1:12)
         DO 180 J = 1, NPTS
            OUTRAY(1,J) = XIN(J,KXAXIS,IPSR)
 180     CONTINUE
      ELSEIF (CXAXIS .EQ. 'PIN') THEN
         CXAXIS(4:4) = '_'
         CXAXIS(5:16) = KSYM(KXAXIS,1)(1:12)
         DO 185 J = 1, NPTS
            OUTRAY(1,J) = PRES(J,IPSR)*XIN(J,KXAXIS,IPSR)*TOTORR
 185     CONTINUE
      ELSEIF (CXAXIS .EQ. 'TIME(SEC)') THEN
         DO 190 J = 1, NPTS
            OUTRAY(1,J) = TIME(J)
 190     CONTINUE
      ELSE
         DO 220 J = 1, NPTS
            OUTRAY(1,J) = FLOAT(J)
 220     CONTINUE
         CXAXIS = 'SOLUTION_#'
      ENDIF
C
C  loop over parameters specified for y-axis data
C
      DO 400 IP = 1, IPSPEC
         IF (PSPEC(IP) .EQ. 'SOLUTION_#') THEN
            DO 265 J = 1, NPTS
               OUTRAY(IP+1,J) = FLOAT(J)
 265        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'EQUIV_RATIO') THEN
            DO 270 J = 1, NPTS
               OUTRAY(IP+1,J) = EQUIV(J,IPSR)
 270        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'PRESSURE(TORR)') THEN
            DO 275 J = 1, NPTS
               OUTRAY(IP+1,J) = PRES(J,IPSR)*TOTORR
 275        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'TAU(1/S)') THEN
            DO 280 J = 1, NPTS
               OUTRAY(IP+1,J) = TAU(J,IPSR)
 280        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'FLRT(SCCM)') THEN
            DO 285 J = 1, NPTS
               OUTRAY(IP+1,J) = FLRT(J,IPSR)*TOSCCM(J)
 285        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'AREA(CM2)') THEN
            DO 290 J = 1, NPTS
               OUTRAY(IP+1,J) = AREA(J,IPSR)
 290        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'VOLUME(CM3)') THEN
            DO 295 J = 1, NPTS
               OUTRAY(IP+1,J) = VOL(J,IPSR)
 295        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'QL(WATTS)') THEN
            DO 300 J = 1, NPTS
               IF (LHTRN(J,IPSR)) THEN
                  OUTRAY(IP+1,J) = HTRN(J,IPSR)*
     &                          (TGAS(J,IPSR) - TAMB(J,IPSR))*TOWATT
               ELSE
                  OUTRAY(IP+1,J) = QL(J,IPSR)*TOWATT
               ENDIF
 300        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'TSURF') THEN
            MP = MPSPEC(IP)
            DO 305 J = 1, NPTS
               OUTRAY(IP+1,J) = TSURF(J,MP,IPSR)
 305        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. '1000/TSURF') THEN
            MP = MPSPEC(IP)
            DO 310 J = 1, NPTS
               OUTRAY(IP+1,J) = 1000./TSURF(J,MP,IPSR)
 310        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'HTRN') THEN
            DO 315 J = 1, NPTS
               OUTRAY(IP+1,J) = HTRN(J,IPSR)
 315        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'TAMB') THEN
            DO 320 J = 1, NPTS 
               OUTRAY(IP+1,J) = TAMB(J,IPSR)
 320        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'TIN') THEN
            DO 325 J = 1, NPTS
               OUTRAY(IP+1,J) = TIN(J,IPSR)
 325        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'TGAS') THEN
            DO 330 J = 1, NPTS
               OUTRAY(IP+1,J) = TGAS(J,IPSR)
 330        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. '1000/TGAS') THEN
            DO 335 J = 1, NPTS
               OUTRAY(IP+1,J) = 1000./TGAS(J,IPSR)
 335        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'TELEC') THEN
            DO 340 J = 1, NPTS
               OUTRAY(IP+1,J) = TELEC(J,IPSR)
 340        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'TELEC(EV)') THEN
            DO 345 J = 1, NPTS
               OUTRAY(IP+1,J) = TELEC(J,IPSR)*TOEV
 345        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'POWER(WATTS)') THEN
            DO 350 J = 1, NPTS
               OUTRAY(IP+1,J) = POWR(J,IPSR)*TOWATT
 350        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'IONE(eV)') THEN
            MP = MPSPEC(IP)
            DO 355 J = 1, NPTS
               OUTRAY(IP+1,J) = EIONSH(J,MP)*TOEV
 355        CONTINUE
            PSPEC(IP)(9:9) = '_'
            PSPEC(IP)(10:10) = 'M'
            PSPEC(IP)(11:11) = CHNUM(MP+1:MP+1)
         ELSEIF (PSPEC(IP) .EQ. 'BIAS(W)') THEN
            MP = MPSPEC(IP)
            DO 357 J = 1, NPTS
               OUTRAY(IP+1,J) = BIASPW(J,MP)/1.E7
 357        CONTINUE
            PSPEC(IP)(8:8) = '_'
            PSPEC(IP)(9:9) = 'M'
            PSPEC(IP)(10:10) = CHNUM(MP+1:MP+1)
         ELSEIF (PSPEC(IP) .EQ. 'TION') THEN
            DO 365 J = 1, NPTS
               OUTRAY(IP+1,J) = TIONP(J)
 365        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'TION(EV)') THEN
            DO 370 J = 1, NPTS
               OUTRAY(IP+1,J) = TIONP(J) * TOEV
 370        CONTINUE
         ELSEIF (PSPEC(IP) .EQ. 'XIN') THEN
            KP = KPSPEC(IP)
            DO 380 J = 1, NPTS
               OUTRAY(IP+1,J) = XIN(J,KP,IPSR)
 380        CONTINUE
            PSPEC(IP)(4:4) = '_'
            PSPEC(IP)(5:16) = KSYM(KP,1)(1:12)
         ELSEIF (PSPEC(IP) .EQ. 'PIN') THEN
            KP = KPSPEC(IP)
            DO 385 J = 1, NPTS
               OUTRAY(IP+1,J) = PRES(J,IPSR)*XIN(J,KP,IPSR)*TOTORR
 385           CONTINUE
            PSPEC(IP)(4:4) = '_'
            PSPEC(IP)(5:16) = KSYM(KP,1)(1:12)
         ELSEIF (PSPEC(IP) .EQ. 'TIME(SEC)') THEN
            DO 390 J = 1, NPTS
               OUTRAY(IP+1,J) = TIME(J)
 390        CONTINUE
         ELSE
            DO 395 J = 1, NPTS
               OUTRAY(IP+1,J) = FLOAT(J)
 395        CONTINUE
            PSPEC(IP) = 'SOLUTION_#'
         ENDIF
 400  CONTINUE
C
      NCOL = NCOL + IPSPEC
      NLAB = 0
C
C loop over species mole fractions as y-axis data
C
      DO 550 N = 1, IKSPEC
         K = KSPEC(N)
         NC = NCOL + N
         L = NLAB + N
         KSTOT = 0
         IMAT = 0
         DO 525 IM = 1, NMAT
            IF (IMAT .EQ. 0) THEN
               KBEG = KKGAS + KSTOT
               KEND = KBEG + KKSURF(IM) + KKBULK(IM)
               IF (K .GT. KBEG .AND. K .LE. KEND) THEN
                  IMAT = IM
                  KM = K - KSTOT
                  DO 520 IPH = 1, NPHASE(IM)
                    IF (KFIRST(IPH,IM) .LE. KM) THEN
                      IF (KLAST(IPH,IM) .GE. KM) THEN
                        IPHASE = IPH
                      END IF
                    END IF
 520              CONTINUE
               ELSEIF (K .LE. KKGAS) THEN
                  IMAT = 1
                  IPHASE = 1
                  KM = K
               ENDIF
               KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
            ENDIF
 525     CONTINUE
         IF (IPHASE .EQ. 1) THEN
            LABL = 'GMF_'
            LNXT = 5
         ELSEIF (IPHASE .LE. NLSURF(IMAT)) THEN
            LABL = 'SSF_'
            LNXT = 5
         ELSE
            LABL = 'ACT_'
            LNXT = 5
         ENDIF
         LT = ILASCH(KSYM(KM,IMAT))
         LABL(LNXT:LNXT+LT-1) = KSYM(KM,IMAT)
         LABEL(L)(1:16) = LABL(1:16)
         DO 530 J = 1, NPTS
            OUTRAY(NC,J) = XMF(J,K,IPSR)
 530     CONTINUE
 550  CONTINUE
      NCOL = NCOL + IKSPEC
      NLAB = NLAB + IKSPEC
C
C convert mole fraction data to number density concentrations
C
      IF (IKSPN .GT. 0) THEN
         DO 560 J = 1, NPTS
           TEMP(1) = TGAS(J,IPSR)
           TEMP(2) = TELEC(J,IPSR)
           TEMP(3) = TIONP(J)
           DO 553 K = 1, KKTOT
              XT(K) = XMF(J,K,IPSR)
 553       CONTINUE
           CALL CKRHOX (PRES(J,IPSR),TEMP,XT,ICKWRK,RCKWRK,RHO)
           CALL CKXTCR (RHO, TEMP(1), XT, ICKWRK, RCKWRK, CONC(1,J))
           DO 558 K = 1, KKTOT
             KSTOT = 0
             IMAT = 0
             DO 557 IM = 1, NMAT
                IF (IMAT .EQ. 0) THEN
                   KBEG = KKGAS + KSTOT
                   KEND = KBEG + KKSURF(IM) + KKBULK(IM)
                   IF (K .GT. KBEG .AND. K .LE. KEND) THEN
                      IMAT = IM
                      DO 556 IPH = 1, NPHASE(IM)
                        IF (KFIRST(IPH,IM) .LE. K) THEN
                          IF (KLAST(IPH,IM) .GE. K) THEN
                            IPHASE = IPH
                          END IF
                        END IF
 556                  CONTINUE
                   ELSEIF (K .LE. KKGAS) THEN
                      IMAT = 1
                      IPHASE = 1
                   ENDIF
                   KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
                ENDIF
 557         CONTINUE
             IF (IPHASE .EQ. 1) THEN
                CONC(K,J) = CONC(K,J) * AVOG
             ELSEIF (IPHASE .LE. NLSURF(IMAT)) THEN
                CONC(K,J) = XT(K) * SDEN(J,IPHASE,IMAT,IPSR) * AVOG
             ELSE
                CONC(K,J) = XT(K)
             ENDIF
 558       CONTINUE
 560     CONTINUE
      ENDIF
C
C loop over species number densities as y-data
C
      DO 580 N = 1, IKSPN
         K = KSPN(N)
         NC = NCOL + N
         L = NLAB + N
         KSTOT = 0
         IMAT = 0 
         DO 573 IM = 1, NMAT
            IF (IMAT .EQ. 0) THEN
               KBEG = KKGAS + KSTOT
               KEND = KBEG + KKSURF(IM) + KKBULK(IM)
               IF ((K .GT. KBEG .AND. K .LE. KEND).OR. K.LE.KKGAS) THEN
                  IMAT = IM
                  KM = K - KSTOT
                  DO 570 IPH = 1, NPHASE(IM)
                     IF (KFIRST(IPH,IM) .LE. KM) THEN
                        IF (KLAST(IPH,IM) .GE. KM) THEN
                           IPHASE = IPH
                        END IF
                     END IF
 570              CONTINUE
               ELSEIF (K .LE. KKGAS) THEN
                  IMAT = 1
                  IPHASE = 1
                  KM = K
               ENDIF
               KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
            ENDIF
 573     CONTINUE
         IF (IPHASE .EQ. 1) THEN
            LABL = 'ND_'
            LNXT = 4
            UNITS = '(/CM3)'
         ELSEIF (IPHASE .LE. NLSURF(IMAT)) THEN
            LABL = 'SSD_'
            LNXT = 5
            UNITS = '(/CM2)'
         ELSE
            LABL = 'ACT_'
            LNXT = 5
            UNITS = ' '
         ENDIF
         LT = ILASCH(KSYM(KM,IMAT))
         LABL(LNXT:LNXT+LT-1) = KSYM(KM,IMAT)
         LNXT = LNXT + LT
         LT = ILASCH(UNITS)
         LABL(LNXT:LNXT+LT-1) = UNITS
         LABEL(L)(1:16) = LABL(1:16)
         DO 575 J = 1, NPTS
            OUTRAY(NC,J) = CONC(K,J)
 575     CONTINUE
 580  CONTINUE
      NCOL = NCOL + IKSPN
      NLAB = NLAB + IKSPN
C
C convert mole fraction data to outlet flow rates of species
C
      IF (IKSPF .GT. 0) THEN
         DO 660 J = 1, NPTS
           TEMP(1) = TGAS(J,IPSR)
           TEMP(2) = TELEC(J,IPSR)
           TEMP(3) = TIONP(J)
           DO 653 K = 1, KKTOT
              XT(K) = XMF(J,K,IPSR)
 653       CONTINUE
           CALL CKRHOX (PRES(J,IPSR),TEMP,XT,ICKWRK,RCKWRK,RHO)
           VOLFOU = FLROUT(J,IPSR) / RHO
           SCCMOU = VOLFOU * (298.15/TGAS(J,IPSR)) 
     1                     * (PRES(J,IPSR)/PATM) * 60.
           DO 658 K = 1, KKGAS
             FLRS(K,J) = SCCMOU*XT(K)
 658       CONTINUE
 660    CONTINUE
      ENDIF
C
C loop over species number densities as y-data
C
      DO 680 N = 1, IKSPF
         K = KSPF(N)
         NC = NCOL + N
         L = NLAB + N
         LABL = 'FR_'
         LNXT = 4
         UNITS = '(SCCM)'
         LT = ILASCH(KSYM(K,1))
         LABL(LNXT:LNXT+LT-1) = KSYM(K,1)
         LNXT = LNXT + LT
         LT = ILASCH(UNITS)
         LABL(LNXT:LNXT+LT-1) = UNITS
         LABEL(L)(1:16) = LABL(1:16)
         DO 675 J = 1, NPTS
            OUTRAY(NC,J) = FLRS(K,J)
 675     CONTINUE
 680  CONTINUE
      NCOL = NCOL + IKSPF
      NLAB = NLAB + IKSPF
C
C
C loop over growth rate data as y-data
C
      IF (LGRALL .OR. IGPHAS .GT. 0) THEN
C
         NDEF = 0
         IF (GRUNIT .EQ. 'ANG/HR') THEN
            GRCONV = 3.6E11
         ELSEIF (GRUNIT .EQ. 'ANG/MIN') THEN
            GRCONV = 6.0E9
         ELSEIF (GRUNIT .EQ. 'ANG/SEC') THEN
            GRCONV = 1.0E8
         ELSEIF (GRUNIT .EQ. 'UM/HR') THEN
            GRCONV = 3.6E7
         ELSEIF (GRUNIT .EQ. 'UM/MIN') THEN
            GRCONV = 6.0E5
         ELSEIF (GRUNIT .EQ. 'UM/SEC') THEN
            GRCONV = 1.0E4
         ELSE
            GRCONV = 6.0E9
            NDEF = 1
         ENDIF
C
C
C       SDOTI AND SITDTI arrays PASSED IN FROM SAVE FILE

         DO 2070 J = 1, NPTS
C
C       Calculate DEN Array 
C
            TEMP(2) = TELEC(J,IPSR)
            TEMP(3) = TIONP(J)
            KSTOT = 0
            DO 2050 IM = 1, NMAT
               DO 800 IPH = 1, NPHASE(IM)
                  SDEN0(IPH,IM) = SDEN(J,IPH,IM,IPSR)
 800           CONTINUE

               DO 853 K = 1, KKTOT
                  XT(K) = XMF(J,K,IPSR)
 853           CONTINUE
               TEMP(1) = TSURF(J,IM,IPSR)
               CALL SKDEN (PRES(J,IPSR), TEMP, XT, 
     1                     SDEN0(1,IM),ISKWRK(IMISK(IM)),
     2                     RSKWRK(IMRSK(IM)), DEN(1,IM))
C
C       Calculate the Growth Rate for bulk phase
C       The basic units are cm/s; the default conversion is to Ang/min
C
               IF (NNBULK(IM).GT.0) THEN
               DO 2000 IPHASE = NFBULK(IM), NLBULK(IM)
                  IF (GRUNIT .EQ. 'MOLE/CM2-S') THEN
                     GRATE(J,IPHASE,IM) = 0.0
                     DO 1600 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                        GRATE(J,IPHASE,IM) = GRATE(J,IPHASE,IM) 
     1                                       + SDOTM(J,KM,IM)
 1600                CONTINUE
                  ELSEIF (GRUNIT .EQ. 'GM/CM2-S') THEN
                     GRATE(J,IPHASE,IM) = 0.0
                     DO 1700 K = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                        GRATE(J,IPHASE,IM) = GRATE(J,IPHASE,IM)
     1                                      + WT(KM,IM)*SDOTM(J,K,IM)
 1700                CONTINUE
                  ELSEIF (GRUNIT .EQ. 'ML/SEC') THEN
                     WTPH = 0.0
                     DENSPH = 0.0
                     GRATE(J,IPHASE,IM) = 0.0
                     DO 1800 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                        K = KM + KSTOT
                        WTPH   = WTPH   + WT(KM,IM)   * XT(K)
                        DENSPH = DENSPH + DEN(KM,IM) * XT(K)
                        GRATE(J,IPHASE,IM) = GRATE(J,IPHASE,IM) +
     1                 WT(KM,IM) / MAX(DEN(KM,IM),SMALL)*SDOTM(J,KM,IM)
 1800                CONTINUE
                     ASPACE = ( WTPH/(DENSPH*AVOG))**(1./3.)
                     GRATE(J,IPHASE,IM) = GRATE(J,IPHASE,IM)/ASPACE
                  ELSE
C     
C    (default growth rate units = Angstrom/minute)
C
                     GRATE(J,IPHASE,IM) = 0.0
                     DO 1900 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                        K = KM + KSTOT
                        GRATE(J,IPHASE,IM) = GRATE(J,IPHASE,IM) +
     1                  WT(KM,IM) / MAX(DEN(KM,IM),SMALL)*SDOTM(J,KM,IM)
 1900                CONTINUE
                     GRATE(J,IPHASE,IM) = GRATE(J,IPHASE,IM)*GRCONV
                     IF (NDEF .EQ. 1) GRUNIT = 'ANG/MIN'
                  ENDIF
 2000          CONTINUE
               ENDIF
               KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
 2050       CONTINUE
 2070    CONTINUE

C
C copy growth rate data to output array
C
         IF (LGRALL) THEN
            DO 2250 IM = 1, NMAT
               IF (NNBULK(IM) .GT. 0) THEN
                  DO 2200 IPHASE = NFBULK(IM), NLBULK(IM)
                     NC = NCOL + 1
                     L = NLAB + 1
                     LABL(1:3) = 'GR_'
                     LT = ILASCH(PSYM(IPHASE,IM))
                     LABL(4:LT+3) = PSYM(IPHASE,IM)
                     LNXT = LT + 4
                     LABL(LNXT:LNXT) = '('
                     LT = ILASCH(GRUNIT) 
                     LABL(LNXT+1:LNXT+LT) = GRUNIT
                     LNXT = LNXT + LT + 1
                     LABL(LNXT:LNXT+1) = ')'
                     LABEL(L)(1:16) = LABL(1:16)
                     DO 2100 J = 1, NPTS
                        OUTRAY(NC,J) = GRATE(J,IPHASE,IM)
 2100                CONTINUE
 2200             CONTINUE
                  NLAB = NLAB + NLBULK(IM) - NFBULK(IM) +1
                  NCOL = NCOL + NLBULK(IM) - NFBULK(IM) +1
               ENDIF
 2250       CONTINUE
         ELSE
            DO 2400 I = 1, IGPHAS
               ISTOT = 0
               IPHASE = 0
               DO 2270 IM = 1, NMAT
                  IF (IPHASE .EQ. 0) THEN
                     IBEG = 1 + ISTOT
                     IEND = NLBULK(IM) + ISTOT
                     IF (KGPHAS(I) .GT. IBEG 
     1                            .AND. KGPHAS(I) .LE. IEND) THEN
                        IMAT = IM
                        IPHASE = KGPHAS(I) - ISTOT
                     ENDIF
                     ISTOT = ISTOT + NNSURF(IM) + NNBULK(IM)
                  ENDIF
 2270          CONTINUE
               NC = NCOL + I
               L = NLAB + I
               LABL(1:3) = 'GR_'
               LT = ILASCH(PSYM(IPHASE,IMAT))
               LABL(4:LT+3) = PSYM(IPHASE,IMAT)
               LNXT = LT + 4
               LABL(LNXT:LNXT) = '('
               LT = ILASCH(GRUNIT) 
               LABL(LNXT+1:LNXT+LT) = GRUNIT
               LNXT = LNXT + LT + 1
               LABL(LNXT:LNXT+1) = ')_'
               LABEL(L)(1:16) = LABL(1:16)
               DO 2300 J = 1, NPTS
                  OUTRAY(NC,J) = GRATE(J,IPHASE,IMAT)
 2300          CONTINUE
 2400       CONTINUE
            NLAB = NLAB + IGPHAS
            NCOL = NCOL + IGPHAS
         ENDIF
      ENDIF
C
C loop over ion current data as y-data
C
      IF (LCURALL .OR. ICURR .GT. 0 .OR. LCURTOT) THEN
C
C
C       Calculate the ion current for each ion using surface rates
C       The basic units are erg/V-cm2-s; convert to A/cm2
C
C
         DO 3200 J = 1, NPTS
            DO 3100 IM = 1, NMAT
               DO 3095 KI = 1, KKION
                  K = KION(KI)
                  CURR(J,KI,IM) 
     1              = -SDOTM(J,K,IM)*ECHG*AVOG*TOWATT*KCHG(K)
     2                 /AFRAC(IM)
 3095          CONTINUE
 3100       CONTINUE
 3200    CONTINUE
C
C copy ion current data to output array
C
         IF (LCURALL) THEN
            DO 3500 KI = 1, KKION
               DO 3450 IM = 1, NMAT
                  K = KION(KI)
                  NC = NCOL + (KI-1)*NMAT + IM
                  L = NLAB + (KI-1)*NMAT + IM
                  LT = ILASCH(KSYM(K,1))
                  LABL = 'J_'
                  LNXT = 3
                  LABL(LNXT:LNXT+LT-1) = KSYM(K,1)
                  LNXT = LNXT+ LT
                  LABL(LNXT:LNXT+1) = '_M'
                  LABL(LNXT+2:LNXT+2)=CHNUM(IM+1:IM+1)
                  LNXT = LNXT+3
                  LABL(LNXT:LNXT) = '('
                  LT = 5
                  LABL(LNXT+1:LNXT+LT) = 'A/CM2'
                  LNXT = LNXT + LT + 1
                  LABL(LNXT:LNXT) = ')'
                  LABEL(L)(1:16) = LABL(1:16)
                  DO 3400 J = 1, NPTS
                     OUTRAY(NC,J) = CURR(J,KI,IM)
 3400             CONTINUE
 3450          CONTINUE
 3500       CONTINUE
            NLAB = NLAB + NMAT*KKION
            NCOL = NCOL + NMAT*KKION
         ELSE
            DO 3700 I = 1, ICURR
               DO 3650 IM = 1, NMAT
                  KI = KCURR(I)
                  NC = NCOL + (I-1)*NMAT + IM
                  L = NLAB + (I-1)*NMAT + IM
                  K = KION(KI)
                  LT = ILASCH(KSYM(K,1))
                  LABL = 'J_'
                  LNXT = 3
                  LABL(LNXT:LNXT+LT-1) = KSYM(K,1)
                  LNXT = LNXT+ LT
                  LABL(LNXT:LNXT+1) = '_M'
                  LABL(LNXT+2:LNXT+2)=CHNUM(IM+1:IM+1)
                  LNXT = LNXT+3
                  LABL(LNXT:LNXT) = '('
                  LT = 5
                  LABL(LNXT+1:LNXT+LT) = 'A/CM2'
                  LNXT = LNXT + LT + 1
                  LABL(LNXT:LNXT) = ')'
                  LABEL(L)(1:16) = LABL(1:16)
                  DO 3600 J = 1, NPTS
                     OUTRAY(NC,J) = CURR(J,KI,IM)
 3600             CONTINUE
 3650          CONTINUE
 3700       CONTINUE
            NLAB = NLAB + ICURR*NMAT
            NCOL = NCOL + ICURR*NMAT
         ENDIF
         IF (LCURTOT) THEN
            DO 3910 IM = 1, NMAT
               NC = NCOL + 1
               L = NLAB + 1
               LABEL(L) = 'JITOT(A/CM2)_M'
               LABEL(L)(15:15) = CHNUM(IM+1:IM+1)
               DO 3900 J = 1, NPTS
                  CURTOT = 0.0
                  DO 3800 KI = 1, KKION
                     K = KION(KI)
                     CURTOT = CURTOT + CURR(J,KI,IM)
 3800             CONTINUE
                  OUTRAY(NC,J) = CURTOT
 3900          CONTINUE
               NLAB = NLAB + 1
               NCOL = NCOL + 1
 3910       CONTINUE
         ENDIF
      ENDIF
C
C loop over ion yield data as y-data
C
      IF (LYLDALL .OR. IYLD .GT. 0 .OR. LYLDTOT) THEN
C
C       Calculate the ion yield for each ion using surface rates
C       The basic units are erg/V-cm2-s; convert to A/cm2
C
         DO 4000 J = 1, NPTS
            DO K = 1, KKTOT
               SDOT(K) = 0.0
            ENDDO
C
            KSTOT = 0
            DO 3930 IM = 1, NMAT
               DO 3920 KM = 1, KKM(IM)
                  IF (KM.GT.KKGAS) THEN
                     K = KM + KSTOT
                  ELSE
                     K = KM
                  ENDIF
                  SDOT(K) = SDOT(K) + SDOTM(J,KM,IM)/AFRAC(IM)
 3920          CONTINUE
               KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
 3930       CONTINUE
C
            KSTOT = 0
            NBETCH = 0
            DO 3990 IM = 1, NMAT
               IF (KKBULK(IM).GT.0) THEN
                  DO 3975 KM = KFIRST(NFBULK(IM),IM),
     2                         KLAST(NLBULK(IM),IM)
                     K = KM + KSTOT
                     IF (SDOT(K).LT.0.0) THEN
                        NBETCH = NBETCH+1
                        SUMFLX = 0.0
                        DO 3970 KI = 1, KKION
                           KG = KION(KI)
                           IF (KCHG(KG).GT.0.AND.SDOT(KG).LT.0.0) THEN
                              YLD(J,KI,NBETCH) = SDOT(K)/SDOT(KG)
                              SUMFLX = SUMFLX + SDOT(KG)
                           ELSE
                              YLD(J,KI,NBETCH) = 0.00
                           ENDIF
 3970                   CONTINUE
                        YLDTOT(J,NBETCH) = SDOT(K)/SUMFLX
                     ENDIF
 3975             CONTINUE
               ENDIF
               KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
 3990       CONTINUE
 4000    CONTINUE
C
C copy ion yield data to output array
C
         IF (LYLDALL) THEN
            NC = NCOL
            L = NLAB
            DO 4120 KI = 1, KKION
               K = KION(KI)
               IF (KCHG(K).GT.0) THEN
                  DO 4110 KB = 1, NBETCH
                     NC = NC+1
                     L = L + 1
                     LT = ILASCH(KSYM(K,1))
                     LABL = 'YLDI_'
                     LNXT = 6
                     LABL(LNXT:LNXT+LT-1) = KSYM(K,1)
                     LNXT = LNXT+ LT
                     LABEL(L)(1:16) = LABL(1:16)
                     DO 4100 J = 1, NPTS
                        OUTRAY(NC,J) = YLD(J,KI,KB)
 4100                CONTINUE
 4110             CONTINUE
                  NLAB = NLAB + NBETCH
                  NCOL = NCOL + NBETCH
               ENDIF
 4120       CONTINUE
         ELSE
            DO 4140 I = 1, IYLD
               KI = KYLD(I)
               NC = NCOL + I
               L = NLAB + I
               K = KION(KI)
               LT = ILASCH(KSYM(K,1))
               LABL = 'YLD_'
               LNXT = 6
               LABL(LNXT:LNXT+LT-1) = KSYM(K,1)
               LNXT = LNXT+ LT
               LABEL(L)(1:16) = LABL(1:16)
               DO 4130 J = 1, NPTS
                  OUTRAY(NC,J) = YLD(J,KI,KB)
 4130          CONTINUE
 4140       CONTINUE
            NLAB = NLAB + IYLD*NBETCH
            NCOL = NCOL + IYLD*NBETCH
         ENDIF
         IF (LYLDTOT) THEN
            DO 4165 KB = 1, NBETCH
               NC = NCOL + 1
               L = NLAB + 1
               LABEL(L) = 'YLDTOT_B'
               LABEL(L)(9:9) = CHNUM(KB+1:KB+1)
               DO 4160 J = 1, NPTS
                  OUTRAY(NC,J) = YLDTOT(J,KB)
 4160          CONTINUE
               NLAB = NLAB + 1
               NCOL = NCOL + 1
 4165       CONTINUE
         ENDIF
      ENDIF
C
C loop over rate-of-production data as y-axis data
C
      KIISTOT = 0
      DO 5000 N = 1, IRSPEC
         K = KRSPEC(N)
         KSTOT = 0
         IMAT = 0
         KIISTOT = 0
         DO 4800 IM = 1, NMAT
            IF (IMAT .EQ. 0) THEN
               KBEG = KKGAS + KSTOT
               KEND = KBEG + KKSURF(IM) + KKBULK(IM)
               IF ((K .GT. KBEG .AND. K .LE. KEND).OR. K.LE.KKGAS) THEN
                  IMAT = IM
                  KM = K - KSTOT
                  DO 4180 IPH = 1, NPHASE(IM)
                     IF (KFIRST(IPH,IM) .LE. KM) THEN
                        IF (KLAST(IPH,IM) .GE. KM) THEN
                           IPHASE = IPH
                        ENDIF
                     ENDIF
 4180             CONTINUE
                  LABL = 'ROP_'
                  LNXT = 5
                  LT = ILASCH(KSYM(KM,IM))
                  LABL(LNXT:LNXT+LT-1) = KSYM(KM,IM)
                  LNXT = LNXT + LT
                  DO 4600 J = 1, NPTS
                     UNITS = '(MOL/CM3-S)'
                     LT = ILASCH(UNITS)
                     LABL(LNXT:LNXT+LT-1) = UNITS
                     DO 4200 I = 1, II
                        NC = NCOL + 1
                        L = NLAB + 1
                        OUTRAY(NC,J) = ROPG(J,I,K)
                        LABEL(L)(1:16) = LABL(1:16)
                        LABEL2(L)(1:16) = IGAS(I)(1:16)
 4200                CONTINUE
                     UNITS = '(MOL/CM2-S)'
                     LT = ILASCH(UNITS)
                     LABL(LNXT:LNXT+LT-1) = UNITS
                     DO 4400 I = 1, IISUR(IM)
                        NC = NCOL + 1
                        L = NLAB + 1
                        OUTRAY(NC,J) = ROPS(J,I,K,IM)
                        LABEL(L)(1:16) = LABL(1:16)
                        LABEL2(L)(1:16) = ISURF(I,IM)(1:16)
 4400                CONTINUE
 4600             CONTINUE
                  KIISTOT = KIISTOT + IISUR(IM)
               ENDIF
               KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
            ENDIF
 4800    CONTINUE
 5000 CONTINUE
      NCOL = NCOL + IRSPEC*(II) + KIISTOT
      NLAB = NLAB + IRSPEC*(II) + KIISTOT
C
C determine number of output files with 8 columns ea.(16 char wide)
C
      MAXCOL = 8
      MAXY = MAXCOL - 1
      NF = NCOL / MAXY
      IF (FLOAT(NCOL)/FLOAT(MAXY) .GT. FLOAT(NF)) NF = NF+1
C
C open files and print output arrays with titles to files
C (parce file names to be 'fname.1', 'fname.2', etc.)
C
      LT = ILASCH(FNAME)
      LDOT = LT+1
      FNAME(LDOT:LDOT) = '.'
      LNUM = LDOT + 1
      IP1 = 1
      IP2 = 0
      IL1 = 1
      IL2 = 0
      I1 = 2
      I2 = MIN(MAXY, NLAB+IPSPEC) +1
      IF (IPSPEC .LE. MAXY) THEN
         IP2 = IPSPEC
         IF (IPSPEC .NE. MAXY) THEN
            IL2 = MIN(NLAB,MAXY-IPSPEC)
         ENDIF
      ELSE
         IP2 = MAXY
         IL2 = 0
      ENDIF
C
C write output x-y array to file; can only make up to 9 files
C (to write more than 9 files, use continuation with new filename)
C
      DO 8000 IFI = 1, MIN(NF,9)
         FILENM(1:LDOT) = FNAME(1:LDOT)
         INUM = CHNUM(IFI+1:IFI+1)
         FILENM(LNUM:LNUM) = INUM(1:1)
         OPEN (LPRNT, STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILENM)
         WRITE (LPRNT,9600) CXAXIS, (PSPEC(IP)(1:15),IP=IP1,IP2),
     1                     (LABEL(IL)(1:15), IL = IL1,IL2)
C         WRITE (LPRNT,9600) BLNK, (BLNK,IP=IP1,IP2),
C     2                     (LABEL2(IL)(1:15), IL = IL1,IL2)
         DO 7800 J = 1, NPTS
            WRITE (LPRNT,9700) OUTRAY(1,J),(OUTRAY(I,J),I=I1,I2)
 7800    CONTINUE
         I1 = I1 + MAXY
         I2 = MIN((IPSPEC+NLAB+1),(I2 + MAXY))
         IP1 = IP1 + MAXY
         IL1 = IL2 + 1
         IF (IP2 .EQ. IFI*MAXY) THEN
            IP2 = MIN(IPSPEC,(IFI+1)*MAXY)
            IF (IP2 .LT. IP1) IP2 = 0
            IF (IPSPEC .LT. (IFI+1)*MAXY) THEN
               IL2 = MIN(NLAB,(IFI+1)*MAXY - IPSPEC)
            ENDIF 
         ELSE
            IP2 = 0
            IL2 = MIN(NLAB,IL2+MAXY)
         ENDIF
         CLOSE (UNIT=LPRNT)
 8000 CONTINUE
C
C
C Write transport-tube inlet file if desired
C
      IF (LTUBE) THEN
C
         LT = ILASCH(FTUBE)
         LDOT = LT+1
         FTUBE(LDOT:LDOT) = '.'
         LNUM = LDOT + 1
         DO 8300 J = 1, NPTS
C
C open file and print composition to file (MAX continuations <100)
C (parce file names to be 'ftube.1', 'ftube.2', etc.)
C
            TUBENM(1:LDOT) = FTUBE(1:LDOT)
            IF (FLOAT(J)/10. .LT. 1.0) THEN
               INUM = CHNUM(J+1:J+1)
               LNDG = 0
            ELSE 
               ITEN = J/10
               INUM(1:1) = CHNUM(ITEN+1:ITEN+1)
               IONE = J - ITEN*10
               INUM(2:2) = CHNUM(IONE+1:IONE+1)
               LNDG = 1
            ENDIF
            TUBENM(LNUM:LNUM+LNDG) = INUM(1:1+LNDG)
            OPEN (LPRNT, STATUS='UNKNOWN',FORM='FORMATTED',FILE=TUBENM)
            WRITE (LPRNT,9750) 'TEMP', TGAS(J,IPSR)
            WRITE (LPRNT,9750) 'PRES', PRES(J,IPSR)*TOTORR/760.
            TEMP(1) = TGAS(J,IPSR)
            TEMP(2) = TELEC(J,IPSR)
            TEMP(3) = TIONP(J)
            DO 8155 K = 1, KKGAS
               XT(K) = XMF(J,K,IPSR)
 8155       CONTINUE
            CALL CKRHOX (PRES(J,IPSR),TEMP,XT,ICKWRK,
     1                   RCKWRK,RHO)
            QCCM = FLRT(J,IPSR) / RHO
            WRITE (LPRNT,9750) 'VDOT', QCCM
            WRITE (LPRNT,9750) 'TEIN', TELEC(J,IPSR)
            WRITE (LPRNT,9750) 'TION', TIONP(J)
            DO 8200 K = 1, KKGAS
               WRITE (LPRNT,9760) 'INIT', KSYM(K,1), XMF(J,K,IPSR)
 8200       CONTINUE
            CLOSE (UNIT=LPRNT)
 8300    CONTINUE
      ENDIF
C
C Write Sensitivity data to file for each solution, if requested
C 
      IF (LSENS) THEN
C
         LT = ILASCH(FSENS)
         LDOT = LT+1
         FSENS(LDOT:LDOT) = '.'
         LNUM = LDOT + 1
         IF (NPTS.GT.99) THEN
            WRITE(LOUT,*) 'Warning...Number of continuations ',
     &                    'in solution file exceeds 99; ',
     &              'Sensitivities will only be written for first 99.'
         ENDIF
         DO 8600 J = 1, MIN(NPTS,99)
C
C open file and print composition to file (MAX continuations <100)
C (parce file names to be 'sens.1', 'sens.2', etc.)
C
            SENSNM(1:LDOT) = FSENS(1:LDOT)
            IF (FLOAT(J)/10. .LT. 1.0) THEN
               INUM = CHNUM(J+1:J+1)
               LNDG = 0
            ELSE 
               ITEN = J/10
               INUM(1:1) = CHNUM(ITEN+1:ITEN+1)
               IONE = J - ITEN*10
               INUM(2:2) = CHNUM(IONE+1:IONE+1)
               LNDG = 1
            ENDIF
            SENSNM(LNUM:LNUM+LNDG) = INUM(1:1+LNDG)
            OPEN (LPRNT, STATUS='UNKNOWN',FORM='FORMATTED',FILE=SENSNM)
C
C hardwire format statement to allow up to 500 columns
C allow ONLY < 1000 reactions.
C
            IIMAX = 999
            WRITE(LPRNT,9765)'REACTION','TEMP', 'TELE', 
     1         (KSYM(K,1), K=1, KKM(1)),(KSYM(K,IM), K=KKGAS+1,KKM(IM))
            RLAB = 'G     '
            SENSNM(LNUM:LNUM+LNDG) = INUM(1:1+LNDG)
            IITOT = II
            DO 8340 IM = 1, NMAT
              IITOT = II + IISUR(IM)
 8340       CONTINUE
            IF (II.GT.IIMAX) WRITE(LOUT,*) 
     1          'Warning...The total number of reactions',
     2          ' exceeds 999; Sensitivities will be ',
     3          ' printed only for the first 999.'
            DO 8400 I = 1, MIN(II,IIMAX)
               IF (FLOAT(I)/10. .LT. 1.0) THEN
                  ISNUM(1:1) = '0'
                  ISNUM(2:2) = '0'
                  ISNUM(3:3) = CHNUM(I+1:I+1)
               ELSEIF (FLOAT(I)/100. .LT. 1.0) THEN
                  ITEN = I/10
                  ISNUM(1:1) = '0'
                  ISNUM(2:2) = CHNUM(ITEN+1:ITEN+1)
                  IONE = I - ITEN*10
                  ISNUM(3:3) = CHNUM(IONE+1:IONE+1)
               ELSE
                  IHUN = I/100
                  ISNUM(1:1) = CHNUM(IHUN+1:IHUN+1)
                  ITEN = (I - IHUN*100)/10
                  ISNUM(2:2) = CHNUM(ITEN+1:ITEN+1)
                  IONE = I - IHUN*100 - ITEN*10
                  ISNUM(3:3) = CHNUM(IONE+1:IONE+1)
               ENDIF
               RLAB(2:4) = ISNUM(1:3)
               WRITE (LPRNT, 9770) RLAB, (SENS(J,I,N), N=1, 2+KKTOT)
 8400       CONTINUE
            RLAB = 'M     '
            DO 8500 IM = 1, MIN(NMAT,9)
               RLAB(2:2) = CHNUM(IM+1:IM+1)
               RLAB(3:3) = 'S'
               IF (IISUR(IM).GT.IIMAX) WRITE(LOUT,*) 
     1            'Warning...The total number of reactions',
     2            ' exceeds 999; Sensitivities will be ',
     3            ' printed only for the first 999.'
               DO 8450 IS = 1, IISUR(IM)
                  IF (FLOAT(IS)/10. .LT. 1.0) THEN
                     ISNUM(1:1) = '0'
                     ISNUM(2:2) = '0'
                     ISNUM(3:3) = CHNUM(IS+1:IS+1)
                  ELSEIF (FLOAT(IS)/100. .LT. 1.0) THEN
                     ITEN = IS/10
                     ISNUM(1:1) = '0'
                     ISNUM(2:2) = CHNUM(ITEN+1:ITEN+1)
                     IONE = IS - ITEN*10
                     ISNUM(3:3) = CHNUM(IONE+1:IONE+1)
                  ELSE
                     IHUN = IS/100
                     ISNUM(1:1) = CHNUM(IHUN+1:IHUN+1)
                     ITEN = (IS - IHUN*100)/10
                     ISNUM(2:2) = CHNUM(ITEN+1:ITEN+1)
                     IONE = IS - IHUN*100 - ITEN*10
                     ISNUM(3:3) = CHNUM(IONE+1:IONE+1)
                  ENDIF
                  RLAB(4:6) = ISNUM(1:3)
                  WRITE (LPRNT, 9770) RLAB,
     1                               (SENSS(J,IS,N,IM), N=1, 2+KKTOT)
 8450          CONTINUE
 8500       CONTINUE
            CLOSE(LPRNT)
 8600    CONTINUE
C
      ENDIF
C
      IF (LCONT) GO TO 20
C
      RETURN
 9600 FORMAT (2X, 8(A14,1X))
 9700 FORMAT (1X, 8(1PE12.5,3X))
 9750 FORMAT (A4, 3x, G12.5)
 9760 FORMAT (A4, 3x, A16, 2x, E12.5)
 9765 FORMAT (A8, 1X, A4, 1X, A4, 1X, 100(A16,1X))
 9770 FORMAT (A6, 2X, 500(1PE10.3,2X))
      END
C----------------------------------------------------------------------C
      SUBROUTINE RDSAVE (LS, LOUT, LENICK, ICKWRK, LENRCK, RCKWRK,
     1                   LENCCK, CCKWRK, LENISK, ISKWRK, LENRSK, RSKWRK,
     2                   LENCSK, CSKWRK, NT, NTE, NYS, NY1,
     3             KMAX, JMAX, IMAX, II, KSYM, IGAS, ISURF, KKGAS,
     4             KKSURF, KKBULK, IISUR, S, EQUIV, PRES, TAU, FLRT, 
     5             VOL, AREA, QL, TSURF, HTRN, TAMB, LHTRN, TIN, XIN,
     6             TEIN, TGAS, TELEC, POWR, XMF, XT, SDEN, 
     5             MAPPH, KEL, KKION, KION, KKM,
     6             MAXPSR, NMAX, NS, SDEN0, WT, ROPG, ROPS,
     7             PSYM, KKPHAS, NPHASE, NNBULK,
     8             NFBULK, NLBULK, NSPHCH, NUMPSR, KFIRST, KLAST, KCHG,
     9             NNSURF, NFSURF, NLSURF, TIONP, EIONSH, MAXMAT, 
     *             NMAT, IMISK, IMRSK, IMCSK, KKTOT, AFRAC, KIMAX,
     1             BHMXI, TIME, MSYM, SDOTI, SITDTI, SDOTM, KERR,BIASPW,
     2             SENS, SENSS, LSOLSN)
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
      DIMENSION ICKWRK(LENICK), ISKWRK(LENISK), MAPPH(NMAX,MAXMAT), 
     1          KKPHAS(NMAX,MAXMAT), KFIRST(NMAX,MAXMAT),
     1          KLAST(NMAX,MAXMAT), KCHG(KMAX), KION(KIMAX), 
     2          IMISK(MAXMAT), IMRSK(MAXMAT), IMCSK(MAXMAT), 
     3          KKSURF(MAXMAT), KKBULK(MAXMAT), 
     3          KKM(MAXMAT), NPHASE(MAXMAT), NNSURF(MAXMAT), 
     4          NFSURF(MAXMAT), NLSURF(MAXMAT), NNBULK(MAXMAT),
     4          NFBULK(MAXMAT), NLBULK(MAXMAT), IISUR(MAXMAT), 
     5          NSPHCH(MAXMAT), S(KMAX), RCKWRK(LENRCK), 
     *          RSKWRK(LENRSK), EQUIV(JMAX,MAXPSR),PRES(JMAX,MAXPSR),
     1          TAU(JMAX,MAXPSR), FLRT(JMAX,MAXPSR), 
     2          VOL(JMAX,MAXPSR), AREA(JMAX,MAXPSR), 
     2          QL(JMAX,MAXPSR), HTRN(JMAX,MAXPSR), 
     3          TAMB(JMAX,MAXPSR), TEIN(JMAX,MAXPSR),
     3          POWR(JMAX,MAXPSR), TGAS(JMAX,MAXPSR), 
     4          TIN(JMAX,MAXPSR), 
     4          XIN(JMAX,KMAX,MAXPSR), XMF(JMAX,KMAX,MAXPSR), 
     5          XT(KMAX),TELEC(JMAX,MAXPSR),
     5          SDEN(JMAX,NMAX,MAXMAT,MAXPSR), 
     6          SDEN0(NMAX,MAXMAT), WT(KMAX,MAXMAT), 
     6          ROPG(JMAX,IMAX,KMAX),
C     7          ROPS(JMAX,IMAX,MAXMAT,*), 
     7          ROPS(JMAX,IMAX,KMAX,MAXMAT),
     8          TIONP(JMAX), EIONSH(JMAX,MAXMAT), AFRAC(MAXMAT), 
     9          BHMXI(MAXPSR),TIME(JMAX), SITDTI(NMAX,MAXMAT),
     9          SDOTI(KMAX,IMAX,MAXMAT),SDOTM(JMAX,KMAX,MAXMAT),
     *          TSURF(JMAX,MAXMAT,MAXPSR), BIASPW(JMAX,MAXMAT),
     1          SENS(JMAX,IMAX,KMAX+2),SENSS(JMAX,IMAX,KMAX+2,MAXMAT)
C
      CHARACTER*18 ICHR
      CHARACTER*60 ISURF(IMAX,MAXMAT), IGAS(*)
      CHARACTER*80 VERSN
      CHARACTER VERS*16, PREC*16, CDUM*16,CCKWRK(LENCCK)*16,
     1          CSKWRK(LENCSK)*16,
     1          KSYM(KMAX,MAXMAT)*16,  
     2          PSYM(NMAX,MAXMAT)*16, MSYM(MAXMAT)*16
      LOGICAL KERR, IERR, LHTRN(JMAX,MAXPSR), MORE, LCHEM, LSOL, LSOLSN
C
      REWIND (LS)
      NS = 0

      KERR = .FALSE.
      IERR = .FALSE.
C
C  Read Header block to see if we have a compatible data file
C
      READ (LS, ERR=8200, END=8100) ICHR
      IF (ICHR .EQ. 'PSR CODE DATA FILE') THEN
        WRITE (LOUT,*) 'RDSAVE: Starting to read PSR CODE DATA FILE'
      ELSE
        WRITE (LOUT,*) 'RDSAVE: First line of binary data file is ',
     1                 'incompatible with expectations'
        KERR = .TRUE.
        RETURN
      END IF
C
C  Check the version number.  
C
      READ (LS, ERR=8500, END=8550) VERSION, VERSN
      WRITE (LOUT,9002)  VERSION, VERSN
C
C  Logicals to check whether solution info exists
C
      LCHEM = .FALSE.
      LSOL  = .FALSE.
      LSOLSN = .FALSE.
C
C-------------------TOP OF THE LOOP-------------------------------------
   10 CONTINUE
      ICHR = ' '
      READ (LS, ERR=8100, END=8100) ICHR
C
      IF (ICHR .EQ. 'CHEMISTRY') THEN
C
C        Initialize Chemkin Library common block
C
         CALL CKPNT (LS, LOUT, NPOINT, VERS, PREC, LENI, LENR, LENC,
     1               IERR)
         KERR = KERR.OR.IERR
C
C*****precision > double
        IF (INDEX(PREC,'DOUB') .LE. 0) THEN
C*****END precision > double
C*****precision > single
C        IF (INDEX(PREC,'SING') .LE. 0) THEN
C*****END precision > single
C
            WRITE (LOUT, 9001) ICHR
            KERR = .TRUE.
            RETURN
        ENDIF
C
C       Initialize Chemkin work arrays
C
        IF (LENI.LE.LENICK .AND. LENR.LE.LENRCK .AND. LENC.LE.LENCCK)
     1     THEN
           READ (LS, ERR=8200, END=8200) (ICKWRK(L), L = 1, LENI)
           READ (LS, ERR=8200, END=8200) (RCKWRK(L), L = 1, LENR)
           READ (LS, ERR=8200, END=8200) (CCKWRK(L), L = 1, LENC)
           CALL CKINDX (ICKWRK, RCKWRK, MM, KK, II, NFIT)
           CALL CKCHRG (ICKWRK, RCKWRK, KCHG)
           CALL CKSYMS (CCKWRK, LOUT, KSYM(1,1), IERR)
           CALL PKINDX (ICKWRK, KEL, KKION)
           KERR = KERR .OR. IERR
           IF (II .LE. IMAX) THEN
              DO 12 I = 1, II
                 CALL CKSYMR (I, LOUT, ICKWRK, RCKWRK, CCKWRK, LT,
     1                        IGAS(I), IERR)
 12           CONTINUE
           ELSE
              WRITE (LOUT, *) ' Error...reaction dimension too small,',
     1                        ' IMAX should be at least ',II
              KERR = .TRUE.
           ENDIF
           IF (KKION .LE. KIMAX) THEN
              CALL CKION (ICKWRK, KION)
           ELSE
              WRITE (LOUT, *) ' Error...ion dimension too small,',
     1                        ' KIMAX should be at least ', KKION
              KERR = .TRUE.
           ENDIF
        ELSE
           READ (LS, ERR=8200, END=8200) (IDUM, L = 1, LENI)
           READ (LS, ERR=8200, END=8200) (RDUM, L = 1, LENR)
           READ (LS, ERR=8200, END=8200) (CDUM, L = 1, LENC)
           WRITE (LOUT, 9000)
     1        ICHR, LENI, LENICK, LENR, LENRCK, LENC, LENCCK
           KERR = .TRUE.
        ENDIF
        NMAT = 0
        LENR = 0
        LENI = 0
        LENC = 0
        IMRSK(1) = 1
        IMISK(1) = 1
        IMCSK(1) = 1
        KKTOT = KK
        LITOT = 0
        LRTOT = 0
        LCTOT = 0
C
 20     CONTINUE
C
        NMAT = NMAT + 1
        IF (NMAT .NE. 1) THEN
           IMISK(NMAT) = IMISK(NMAT-1) + LENI
           IMRSK(NMAT) = IMRSK(NMAT-1) + LENR
           IMCSK(NMAT) = IMCSK(NMAT-1) + LENC
        ENDIF
C
C        Initialize Surface Library common block
C
        CALL SKPNT (LS, LOUT, VERS, PREC, LENI, LENR, LENC,
     1                  IERR)
        LITOT = LITOT + LENI
        LRTOT = LRTOT + LENR
        LCTOT = LCTOT + LENC
        KERR = KERR.OR.IERR
C
C*****precision > double
        IF (INDEX(PREC,'DOUB') .LE. 0) THEN
C*****END precision > double
C*****precision > single
C        IF (INDEX(PREC,'SING') .LE. 0) THEN
C*****END precision > single
C 
           WRITE (LOUT, 9001) ICHR
           KERR = .TRUE.
           RETURN
        ENDIF
C
C      Initialize Surface Chemistry work arrays
C
        IF (LITOT.LE.LENISK .AND. LRTOT.LE.LENRSK 
     1                     .AND. LCTOT.LE.LENCSK) THEN
           READ (LS, ERR=8200, END=8200) (ISKWRK(L), L = IMISK(NMAT), 
     1                                             IMISK(NMAT)+LENI-1)
           READ (LS, ERR=8200, END=8200) (RSKWRK(L), L = IMRSK(NMAT), 
     1                                             IMRSK(NMAT)+LENR-1)
           READ (LS, ERR=8200, END=8200) (CSKWRK(L), L = IMCSK(NMAT), 
     1                                             IMCSK(NMAT)+LENC-1)
           CALL SKINDX (ISKWRK(IMISK(NMAT)), NELEM, KKGAS,KKSURF(NMAT), 
     1                  KKBULK(NMAT), KKM(NMAT),NPHASE(NMAT), 
     2                  NNSURF(NMAT), NFSURF(NMAT), NLSURF(NMAT), 
     3                  NNBULK(NMAT), NFBULK(NMAT),NLBULK(NMAT),
     4                  IISUR(NMAT))
C  Extract the phase pointer vectors
           CALL SKPKK (ISKWRK(IMISK(NMAT)), KKPHAS(1,NMAT), 
     1                 KFIRST(1,NMAT), KLAST(1,NMAT))            
C  Extract the material name 
           CALL SKSYMM (ISKWRK(IMISK(NMAT)), CSKWRK(IMCSK(NMAT)), 
     1                  LOUT, MSYM(NMAT), IERR)
C  Extract the phase names into a vector
           IF (NPHASE(NMAT) .LE. NMAX) THEN
              CALL SKSYMP (ISKWRK(IMISK(NMAT)), CSKWRK(IMCSK(NMAT)), 
     1                     LOUT, PSYM(1,NMAT), IERR)
              KERR = KERR .OR. IERR
           ELSE
              WRITE (LOUT, *) ' Error...phase dimension too small,',
     1                        ' NMAX should be at least ', NPHASE(IM)
              KERR = .TRUE.
              RETURN
           ENDIF
C  Extract the species names into a vector
           KKTOT = KKTOT + KKSURF(NMAT) + KKBULK(NMAT)
           IF (KKTOT .LE. KMAX) THEN
              CALL SKSYMS (ISKWRK(IMISK(NMAT)), CSKWRK(IMCSK(NMAT)), 
     1                     LOUT, KSYM(1,NMAT), IERR)
           ELSE
              WRITE (LOUT, *) ' Error...species dimension too small,',
     1                         ' KMAX should be at least ',KKTOT
              KERR = .TRUE.
           ENDIF
           IF (IISUR(NMAT) .LE. IMAX) THEN
              DO 22 IR = 1, IISUR(NMAT)
                 CALL SKSYMR (IR, LOUT, ISKWRK(IMISK(NMAT)), 
     1                        RSKWRK(IMRSK(NMAT)), CSKWRK(IMCSK(NMAT)), 
     2                        LT,ISURF(IR,NMAT), IERR)
 22           CONTINUE
           ELSE
              WRITE (LOUT, *) ' Error...reaction dimension too small,',
     1                         ' IMAX should be at least ',IISUR(NMAT)
              KERR = .TRUE.
           ENDIF
C   Extract the standard state site densities from the work array
           CALL SKSDEN (ISKWRK(IMISK(NMAT)), RSKWRK(IMRSK(NMAT)), 
     1                   SDEN0(1,NMAT))
C   Extract the molecular weights for the species
           CALL SKWT(ISKWRK(IMISK(NMAT)), RSKWRK(IMRSK(NMAT)), 
     1                WT(1,NMAT))

        ELSE
           READ (LS, ERR=8200, END=8200) (IDUM, L = 1, LENI)
           READ (LS, ERR=8200, END=8200) (RDUM, L = 1, LENR)
           READ (LS, ERR=8200, END=8200) (CDUM, L = 1, LENC)
           WRITE (LOUT, 9200)
     1        ICHR, LITOT, LENISK, LRTOT, LENRSK, LCTOT, LENCSK
           KERR = .TRUE.
        ENDIF
        MORE = .FALSE.
        IF (ISKWRK(IMISK(NMAT) + LENI - 1) .EQ. 1) THEN
           MORE = .TRUE.
        ENDIF
        IF (MORE) GO TO 20
        IF (.NOT.KERR) LCHEM = .TRUE.
C
      ELSEIF (ICHR .EQ. 'SOLUTION') THEN
         IF (KERR) RETURN
C
         READ (LS, ERR=8300, END=8350) NATJ, IPSR, NUMPSR, NMAT,
     1                         (NSPHCH(IM),IDUMIS, IM=1,NMAT)
         IF (IPSR.EQ.1) NS = NS + 1
         IF (NS .GT. JMAX) THEN
            WRITE (LOUT, *) ' Error...dimension too small,',
     1                      ' JMAX should be at least ', NS
            KERR = .TRUE.
            RETURN
         ENDIF
C
C Check bounds on arrays:
C
         IF (NATJ .GT. KMAX) THEN
            WRITE (LOUT, *) ' Error...solution dimension too small,',
     1                      ' KMAX should be at least ',NATJ
            KERR = .TRUE.
            RETURN
         ENDIF
         IF (IPSR .GT. MAXPSR) THEN
            WRITE (LOUT, *) ' Error...solution dimension too small,',
     1                      ' MAXPSR should be at least ',IPSR
            KERR = .TRUE.
            RETURN
         ENDIF
         NSPTOT = 0
         DO 3000 IM = 1, NMAT
            IF (NSPHCH(IM) .GT. NMAX) THEN
               WRITE (LOUT, *) ' Error...solution dimension too small,',
     1                         ' NMAX should be at least ', NSPHCH(IM)
            ENDIF
            NSPTOT = NSPTOT + NSPHCH(IM)
 3000    CONTINUE
C
         READ (LS, ERR=8300, END=8350) EQUIV(NS,IPSR), PRES(NS,IPSR),
     1          TAU(NS,IPSR),  FLRT(NS,IPSR), VOL(NS,IPSR),
     2          AREA(NS,IPSR), QL(NS,IPSR),
     3         (TSURF(NS,IM,IPSR),IM=1,NMAT),
     4          HTRN(NS,IPSR), TAMB(NS,IPSR), LHTRN(NS,IPSR),
     5         (AFRAC(IM),IM=1,NMAT)
         READ (LS, ERR=8300, END=8350) TIN(NS,IPSR),
     1                               (XIN(NS,K,IPSR), K = 1, KKGAS)
         IF (KEL .NE. 0) THEN
           READ (LS, ERR=8300, END=8350) TEIN(NS,IPSR), POWR(NS,IPSR),
     1                             TIONP(NS), BHMXI(IPSR),
     2             (EIONSH(NS,IM),ESHTH, BIASPW(NS,IM),IM=1,NMAT),
     3                              QLEX
         ENDIF
         IF (NSPTOT .GT. 0) THEN
           READ (LS, ERR=8300, END=8350) (S(K), K = 1, NATJ),
     1             ((MAPPH(I,IM), I = 1, MAX(1,NSPHCH(IM))), IM=1,NMAT)
         ELSE
           READ (LS, ERR=8300, END=8350) (S(K), K = 1, NATJ)
         ENDIF
         READ (LS, ERR=8300, END=8350) 
     1    (((SDOTI(K,I,IM),K=1,KKM(IM)),I=1,MAX(1,IISUR(IM))),
     2      IM=1,NMAT), 
     2     ((SITDTI(N,IM), N=1, NPHASE(IM)),IM=1,NMAT)
         READ (LS, ERR=8300, END=8350) SOLTIM
C
         TIME(NS) = SOLTIM
         TGAS(NS,IPSR) = S(NT)
         TELEC(NS,IPSR) = S(NTE)
         CALL CKYTX (S(NY1), ICKWRK, RCKWRK, XT)
         DO 6700 K = 1, KKGAS
            XMF(NS,K,IPSR) = XT(K)
 6700    CONTINUE
         IF (KKTOT .GT. KKGAS) THEN
           DO 6900 K = KKGAS+1,KKTOT
              XMF(NS,K,IPSR) = S(NYS+K)
 6900      CONTINUE
         END IF
         IDTOT = 0
         DO 7250 IM = 1, NMAT
            DO K = 1, KKM(IM)
               SDOTM(NS,K,IM) = 0.0
            ENDDO
            DO I = 1, IISUR(IM)
               DO K = 1, KKM(IM)
                  SDOTM(NS,K,IM) = SDOTM(NS,K,IM) + SDOTI(K,I,IM)
               ENDDO
            ENDDO
            DO 7100 I = 1, NPHASE(IM)
               SDEN(NS,I,IM,IPSR) = SDEN0(I,IM)
 7100       CONTINUE
            IF (NSPHCH(IM) .GT. 0) THEN
              DO 7200 IDM = 1, NSPHCH(IM)
                 I = IDM + IDTOT
                 SDEN(NS,MAPPH(IDM,IM),IM,IPSR) = 
     1                   S(KKTOT+1+I) * SDEN0(MAPPH(IDM,IM),IM)
 7200         CONTINUE
            END IF
            IDTOT = IDTOT + NSPHCH(IM)
 7250    CONTINUE
         IF (.NOT. KERR) LSOL = .TRUE.
C-----------------------------------------------------------------------
      ELSEIF (ICHR .EQ. 'RATE OF PRODUCTION') THEN
         DO 7300 K = 1, KKTOT
             READ (LS, ERR=8370, END=8380) KDUM,
     1      (ROPG(NS,I,K), I=1,II),
     2      ((ROPS(NS,I,K,IM), I=1,MAX(1,IISUR(IM))),IM=1,NMAT)
 7300    CONTINUE
C-----------------------------------------------------------------------
      ELSEIF (ICHR .EQ. 'SENSITIVITY') THEN
         NNBTOT = 0
         DO 7350 IM = 1, NMAT
            NNBTOT = NNBTOT + NNBULK(IM)
 7350    CONTINUE
         DO 7400 I = 1, II
            IF (NNBTOT .GT. 0) THEN
               READ (LS, ERR=8400, END=8450) IS,
     1           ( SENS(NS,I,N), N=1,NATJ ),
     2           (( GRSG,IPHASE=NFBULK(IM),NLBULK(IM) ),IM=1,NMAT)
            ELSE
               READ (LS, ERR=8400, END=8450) IS,
     1           ( SENS(NS,I,N), N=1,NATJ ) 
            ENDIF
 7400    CONTINUE
         DO 7550 IM = 1, NMAT
            DO 7500 I = 1, IISUR(IM)
               IF (NNBULK(IM) .GT. 0) THEN
                  READ (LS, ERR=8400, END=8450) IS,
     1              (SENSS(NS,I,N,IM), N=1,NATJ),
     1              (GRSS, IPHASE=NFBULK(IM),NLBULK(IM))
               ELSE
                  READ (LS, ERR=8400, END=8450) IS,
     1              (SENSS(NS,I,N,IM), N=1,NATJ)
               ENDIF
 7500       CONTINUE
 7550    CONTINUE
         IF (.NOT. KERR) LSOLSN = .TRUE.
      ENDIF
      IF (.NOT. KERR) GO TO 10
      RETURN
C
C----------------RETURN AFTER NORMAL END OF DATA FILE-------------------
 8100 CONTINUE
      IF (.NOT. LCHEM) THEN
         KERR = .TRUE.
         WRITE(LOUT,*)'RDSAVE: Error, did not find chemistry data'
      ELSEIF (.NOT. LSOL) THEN
         KERR = .TRUE.
         WRITE(LOUT,*)'RDSAVE: Error, did not find solution data'
      ELSE
         WRITE(LOUT,*)'RDSAVE: Finished reading psr data file'
      ENDIF
      RETURN
C--------------------STOP AFTER ERROR CONDITION-------------------------
 8200 CONTINUE
      WRITE  (LOUT, 9210) ICHR
 9210 FORMAT (' Error reading ', A, ' data...')
      KERR = .TRUE.
      RETURN
C
 8300 CONTINUE
      WRITE (LOUT, *) ' Error reading solution...'
      WRITE (LOUT, *) ' NATJ=',NATJ,', NSOL=',NS,', IPSR = ', IPSR
      KERR = .TRUE.
      RETURN
 8350 CONTINUE
      WRITE (LOUT, *) ' End of file read while reading the solution'
      WRITE (LOUT, *) ' NATJ=',NATJ,', NSOL=',NS,', IPSR = ', IPSR,
     1                ' NMAT=',NMAT
      KERR = .TRUE.
      RETURN
C
 8370 CONTINUE
      WRITE (LOUT, *) ' Error reading gas-phase rates of progress...'
      KERR = .TRUE.
      RETURN
C
 8380 CONTINUE
      WRITE (LOUT, *) ' End of file read while reading rates of ',
     1                'progress'
      KERR = .TRUE.
      RETURN
C
 8400 CONTINUE
      WRITE (LOUT, *) ' Error reading gas-phase sensitivities...'
      WRITE (LOUT, *) ' NSOL=',NS,', II=',II,', I=',I
      KERR = .TRUE.
      RETURN
C
 8450 CONTINUE
      WRITE (LOUT, *) ' End of file read while reading sensitivities'
      KERR = .TRUE.
      RETURN
C
 8500 CONTINUE
      WRITE (LOUT, *) ' Error: Could not read the first record in the '
     1               ,' binary data file.'
      KERR = .TRUE.
      RETURN
C
 8550 CONTINUE
      WRITE (LOUT, *) ' Error: Could not read the second record in ',
     1                ' binary data file.'
      WRITE (LOUT, *) '        One possible error may be that the ',
     1                ' floating point data type is incompatible.'
      KERR = .TRUE.
      RETURN
C
 9002 FORMAT (' RDSAVE: Version number of data file = ', F5.2,
     1      /,' RDSAVE: ',A80)
 9001 FORMAT (' Error...',A,
     1        ' precision incompatible with post-processor')
 9000 FORMAT (/,' Error...not enough ',A,' work space provided:',
     1        /,'                REQUIRED        PROVIDED',
     2        /,' INTEGER  ', 2I15,
     3        /,' REAL     ', 2I15,
     4        /,' CHARACTER', 2I15,/)
 9200 FORMAT (/,' Error...not enough ',A,
     1          ' work space provided for SURFACE CHEMKIN work space:',
     2        /,'                REQUIRED        PROVIDED',
     3        /,' INTEGER  ', 2I15,
     4        /,' REAL     ', 2I15,
     5        /,' CHARACTER', 2I15,/)
C
      END
C----------------------------------------------------------------------C
      SUBROUTINE RDKEY (LIN, LOUT, LCONT, FNAME, IPSR, CXAXIS, KXAXIS,
     1                  PSPEC, IPSPEC, KPSPEC, KKM, KSPEC, IKSPEC,
     2                  KGPHAS, IGPHAS, GRUNIT, LGRALL, KKGAS,
     3                  KSYM, PSYM, MSYM, NPHASE, KKPHAS, KRSPEC, 
     4                  IRSPEC, MPSPEC, MXAXIS, KMAX, NNBULK,
     5                  NFBULK, NLBULK, NSPHCH, NUMPSR,
     6                  KFIRST, KLAST, MAXPAR, MAPPH, KSPN, IKSPN, 
     7                  KSPF, IKSPF,
     8                  NFSURF, NLSURF, KEL, KKION, KION, LCURALL, 
     9                  KCURR, ICURR, LCURTOT, NMAT, IMISK, IMRSK, 
     *                  IMCSK, KKTOT, NMAX, KKSURF, KKBULK, NNSURF,
     2                  LYLDALL, KYLD, IYLD, LYLDTOT, LTUBE, FTUBE,
     2                  LSENS, FSENS, LSOLSN, KERR)
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C LIST OF KEYWORDS 
C--------------------------------------
C
C END
C       Signals end of the current keyword input
C
C CNTN
C       Signals a new plotfile should be made after the current one
C
C XAXS  [Character expression]  [Species name]*16
C       Specifies the variable for the x axis.  Currently, the
C       following variables may be used as the x-axis
C
C   * THE TOTAL NUMBER OF PARAMETERS IN THE FOLLOWING LIST          *
C   * MUST BE LESS THAN 'MAXPAR' SPECIFIED IN A PARAMETER STATEMENT *
C
C       NS........................Continuation Number (DEFAULT)
C       EQUIV.....................Equivalence Ratio
C       PTOTAL....................Total Reactor Pressure (torr)
C       TAU ......................Residence time (sec)
C       FLRT......................Flow rate (sccm)
C       AREA......................Surface area (cm**2)
C       VOLUME....................Reactor Volume (cm**3)
C       QL........................Rate of Heat loss to the ambient
C       TSURF [Material name]*16..Surface temperature
C       TSURF [Mat name]*16 ARRH..1000 / [Surface temperature]
C       HTRN......................Heat transfer coeff. [erg/cm2-s-K]
C       TAMB......................Ambient temperature
C       TGAS......................Gas Temperature
C       TGAS  ARRH................1000/[Gas Temperature]
C       TIN.......................Inlet temperature
C       XIN [Species name]*16.....Inlet mole fraction of a species
C       PIN [Species name]*16.....Inlet partial pressure of a species
C                                      [PIN = PTOTAL*XIN}
C       TELEC.....................Electron temperature in [K]
C       TELEC EV..................Electron temperature in [eV]
C       POWER.....................Power deposition to electrons [Watts]
C       IONE [Material name]*16...Ion energy after traversing sheath
C                                 [eV] at a specified material
C       BIAS [Material name]*16...Bias power at a specified material [W]
C       TION......................Assumed ion temperature in [K]
C       TION EV...................Assumed ion temperature in [eV]
C       TIME......................Solution time for transient solution
C                                 [seconds]
C       
C
C GRTH [ALL] [Phase name]*16 [Units]
C       PLots the growth rate/etch rate of either all bulk phases
C       or just one bulk phase.  The units may be specified on
C       the command line as well.  Possible Units are:
C          ANG/HR
C          ANG/MIN
C          ANG/SEC
C          UM/HR
C          UM/MIN
C          UM/SEC
C          MOLE/CM**2*SEC
C          GM/CM**2*SEC
C          ML/SEC      monolayer per second
C       If units are not specified angstroms/min are assumed.
C
C MOLS [ALL] [Phase name]*16 [Species name]*16
C       Plot the mole fraction/site fraction
C       of the species specified by the
C       character string.  Up to 10
C       species can be specified on the same keyword line.
C       Example: (mf = mole fraction, sf = site fraction
C          MOLS ALL         Plot all mf in all phases
C          MOLS H20         Plot H2O mf in first phase that it appears
C          MOLS GAS H20     Plot H2O mf in gas phase
C          MOLS BULK1 H20   Plot H20 mf in "BULK1" bulk phase
C          MOLS H2O/SITE2/  Plot H2O sf "SITE2" surface phase
C          MOLS SITE1       Plot all sfs in "SITE1" surface phase
C
C NUMD [ALL] [phase name]*16 [species name]*16 
C       plot number densities of the species [#/cm**3]
C
C OFRS [ALL] [phase name]*16 [species name]*16 
C       plot outlet flow rates of the species [sccm]
C
C RATE [ALL] [Phase name]*16 [Species name]*16 [STREAMS] [TOTAL]
C       Plots the rate of production of the species
C       Units for the plot are in moles/cm**3.
C
C CURR [ALL] [Species name]*16
C      plot ion current densities of the ion species [A/cm2]
C
C YLDI [ALL] [Ion Species name]*16 
C      plot the etch yield of the ion for all bulk species etched 
C
C PRAM  [ALL] [Character expression]  [Species name]*16
C       Plots any of the problem parameters specified in the
C       XAXS keyword description against any other parameter.
C       The arguments to the keyword are the same as the XAXS
C       keyword.
C       If [ALL] is specified, then all parameters are plotted.
C
C SENS [filename]*60
C       List the sensitivity data for all species (columns) vs.
C       all reactions (rows).  Output to separate file for each
C       solution found in the save.bin file.  Must have one or
C       more sensitivity keyword in the aurora input file to work.
C       Default filename is 'sens'.  Filenames are appended with
C       '.#, where # is the aurora solution number'
C
C IPSR  [Integer]
C       Plot the solution from the IPSR'th PSR that is connected in
C       series (default = 1)
C--------------------------------------
C  OUTPUT VARIABLES
C ------------------
C LCONT   Logical variable specifying whether another plot should be
C         made after the current one.
C
C IPSR    Number of the PSR to plot.
C
C LGRALL   Flag for doing the growth rate of all phases plot
C
C IGPHAS   Number of individual phase growth rate plots
C
C KGPHAS(NMAX) Vector of phase numbers for doing growth rate plots
C
C GRUNIT  Character*16 of units for the growth rate plot
C
C LCURALL   Flag for doing the ion currents of all ions plot
C
C LCURTOT   Flag for calculating the total ion current
C
C ICURR  Number of individual ion current plots
C
C KCURR(KIMAX) Vector of ion species numbers for doing current plots
C
C LYLDALL    Flag for doing yield plots for all ions
C
C LYLDTOT   Flag for calculating the total ion yield for etched bulks
C
C IYLD   Number of individual ion-yield plots
C
C KYLD(KIMAX)  Vector of ion species numbers for doing yield plots
C
C IKPHAS Number of phase mole fraction plots to be made
C
C IKSPEC  Number of species mole fraction plots to be made
C
C KSPEC(KMAX)  Integer vector indicating which species to
C              be plotted (mole fractions)
C
C IKSPN   Number of species number densities plots to be made
C
C KSPN(KMAX)   Integer vector indicating which species to 
C              be plotted (number density)
C
C IKSPF   Number of species flow-rate plots to be made
C
C KSPF(KMAX)   Integer vector indicating which species to 
C              be plotted (flow rates)
C
C IRSPEC  Number of rate of production plots to be made
C
C KRSPEC(KMAX)    Integer vector, indicating which species to
C                be plotted (rate of productions). 
C
C IPSPEC    Number of parameter y-variables to be listed
C
C PSPEC(*)  Character*16 array identifying the y-axis parameters
C           (other than mole fractions, or rops)
C KPSPEC(*) Integer variable identifying the species whose inlet mole
C           fraction or partial pressure was identified for the y-axis
C
C MPSPEC(*) Integer variable identifying the material at which the ion
C           energy was identified for the y-axis
C
C CXAXIS    Character*16 identifying x-axis parameter
C
C KXAXIS    Integer variable identifying the species whose inlet mole
C           fraction or partial pressure should be plotted on the
C           x axis.
C MXAXIS    Integer variable identifying the material on which the
C           ion energy should be plotted on the x axis.
C ----------------------
C  NECESSARY CONDITIONS
C ----------------------
C
C  1      One of the following keywords must be found:
C             PRAM, MOLS, NUMD, RATE, CURR, GRTH, YLDI, OFRS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
      PARAMETER (NDIM=10)
C
      CHARACTER KEY*4, LINE*128, KDUM(NDIM)*80,
     1          CXAXIS*16, PSPEC(*)*(*), GRUNIT*16, KSYM(KMAX,*)*(16),
     2          PSYM(NMAX,*)*(16),SUB(NDIM)*20,FNAME*60, FTUBE*60,
     3          MSYM(*)*(16), FSENS*60
      CHARACTER*80 OPTION, CKCHUP
      DIMENSION NRAY(NDIM), KPSPEC(*), KSPEC(*), MPSPEC(*),
     1          KGPHAS(*), KKPHAS(NMAX,*), KRSPEC(*), 
     2          NPHASE(*), KKM(*),
     3          KFIRST(NMAX,*), KLAST(NMAX,*), MAPPH(NMAX,*),
     4          KSPN(*), KCURR(*), KION(*), KYLD(*), KSPF(*),
     5          IMISK(*), IMRSK(*), IMCSK(*), KKSURF(*), KKBULK(*),
     6          NNSURF(*), NNBULK(*), NFSURF(*), NLSURF(*), NFBULK(*),
     7          NLBULK(*), NSPHCH(*), VALUE(5)
C
      LOGICAL LCONT, NEC(1), KERR, IERR, LGRALL, LCURALL, LCURTOT,
     2        LYLDALL, LYLDTOT, LTUBE, LSENS, LSOLSN
      INTEGER CKSLEN, CKLSCH
      EXTERNAL CKSLEN, CKCHUP, CKLSCH
C      
C Initialize variables to defaults:
C
      LGRALL        = .FALSE.
      LCONT         = .FALSE.
      LCURALL       = .FALSE.
      LCURTOT       = .FALSE.
      LYLDALL       = .FALSE.
      LYLDTOT       = .FALSE. 
      LTUBE         = .FALSE.
      LSENS         = .FALSE.
      KERR          = .FALSE.
      IERR          = .FALSE.
      IPSR          = 1
      IGPHAS        = 0
      ICURR         = 0
      IYLD          = 0
      IKSPEC        = 0
      IRSPEC        = 0
      IPSPEC        = 0
      IKSPN         = 0
      IKSPF         = 0
      NEC(1)  = .FALSE.
      CXAXIS = 'NS'
      FNAME = ' '
      FTUBE = ' '
      FSENS = ' '
      DO 5 I = 1, NDIM
        KDUM(I) = ' '
5       SUB(I) = ' '
      DO 7 I = 1, KKTOT
        KSPEC(I)   = 0
        KRSPEC(I) = 0
        KSPN(I) = 0
        KSPF(I) = 0
7     CONTINUE
      DO 8 I = 1, NMAX
        KGPHAS(I) = 0
 8    CONTINUE
      DO 9 I = 1, MAXPAR
        KPSPEC(I) = 0
        MPSPEC(I) = 0
 9    CONTINUE
C
C
   10 CONTINUE
      KEY = ' '
      LINE = ' '
      READ (LIN, '(A)', END=500, ERR=600) LINE
      WRITE (LOUT, '(A)') LINE(1:CKLSCH(LINE))
      LT = CKSLEN(LINE)
      IF (LT .LE. 0) GO TO 10
      IF (LT .LT. LEN(LINE)) LINE(LT+1:) = ' '
C
      KEY = CKCHUP(LINE(1:4), 4)
      LINE(1:4) = ' '
C
C skip line if it is a comment
C
      IF (KEY(1:1) .EQ. ' ' .OR.
     1    KEY(1:1) .EQ. '.' .OR.
     2    KEY(1:1) .EQ. '/' .OR.
     3    KEY(1:1) .EQ. '!'       )      GO TO 10
C
C Strip Line of "!" Comments
C
C continuation marker
C
      IF (KEY .EQ. 'CNTN') THEN
         LCONT = .TRUE.
         RETURN
C
C end marker
C
      ELSEIF (KEY(1:3) .EQ. 'END') THEN
         IF (KERR) THEN
            WRITE (LOUT, *) ' Error in keyword input '
            RETURN
         ENDIF
         LCONT = .FALSE.
         RETURN
C
C plot-file name for x-y data
C
      ELSEIF (KEY .EQ. 'FILE') THEN
         LF = ILASCH(LINE)
         L1=1
         DO 208 LI = 1, LF
           IF (LINE(LI:LI).EQ.' ') L1=LI+1
208      CONTINUE
         FNAME = LINE(L1:LF)
c
c keyword to generate transport tube inlet file and filename
c
      ELSEIF (KEY .EQ. 'TUBE') THEN
         LTUBE = .TRUE.
         LF = ILASCH(LINE)
         L1=1
         DO 209 LI = 1, LF
           IF (LINE(LI:LI).EQ.' ') L1=LI+1
 209     CONTINUE
         FTUBE = LINE(L1:LF)
      ELSEIF (KEY .EQ. 'SENS') THEN
         LSENS = .TRUE.
         IF (.NOT. LSOLSN) THEN
            WRITE(LOUT,*)'RDKEY WARNING: No sensitivity data found',
     &                   '; keyword ignored.'
            LSENS = .FALSE.
         ENDIF
         IF (LSENS) THEN
            LF = ILASCH(LINE)
            L1=1
            DO 210 LI = 1, LF
               IF (LINE(LI:LI).EQ.' ') L1=LI+1
 210        CONTINUE
            FSENS = LINE(L1:LF)
            IF (FSENS .EQ. ' ') FSENS = 'sens'
         ENDIF
C
C x-axis data specification
C####---------------------------------------------------------------XAXS
C
      ELSEIF (KEY .EQ. 'XAXS') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
         IF (NFOUND .GT. 0) THEN
           IF (SUB(1) .EQ. 'NS') THEN
             CXAXIS = 'SOLUTION_#'
           ELSE IF (SUB(1) .EQ. 'EQUIV') THEN
             CXAXIS = 'EQUIV_RATIO'
           ELSE IF (SUB(1) .EQ. 'PTOTAL') THEN
             CXAXIS = 'PRESSURE(TORR)'
           ELSE IF (SUB(1) .EQ. 'VOLUME') THEN
             CXAXIS = 'VOLUME(CM3)'
           ELSE IF (SUB(1) .EQ. 'AREA') THEN
             CXAXIS = 'AREA(CM2)'
           ELSE IF (SUB(1) .EQ. 'FLRT') THEN
             CXAXIS = 'FLRT(SCCM)'
           ELSE IF (SUB(1) .EQ. 'TAU') THEN
             CXAXIS = 'TAU(1/S)'
           ELSE IF (SUB(1) .EQ. 'QL') THEN
             CXAXIS = 'QL(WATTS)'
           ELSE IF (SUB(1) .EQ. 'POWER') THEN
             CXAXIS = 'POWER(WATTS)'
           ELSE IF (SUB(1) .EQ. 'INLOSS') THEN
             CXAXIS = 'INLOSS(EV)'
           ELSE IF (SUB(1) .EQ. 'TIME') THEN
             CXAXIS = 'TIME(SEC)'
           ELSE IF (SUB(1) .EQ. 'TSURF') THEN
             IF (NFOUND .GE. 2) THEN
                LT = ILASCH(SUB(2))
                IF (SUB(2)(1:LT) .EQ. 'ARRH') THEN
                   CXAXIS = '1000/TSURF'
                   IF (NFOUND .GE. 3) THEN
                      LT = ILASCH(SUB(3))
                      CALL CKCRAY(SUB(3)(1:LT), NMAT, MSYM,
     1                            LOUT, NDIM, NRAY, NF, IERR)
                      IF (NF .GT. 0 .AND. .NOT. IERR) THEN
                         MXAXIS = NRAY(1)
                      ELSE
                         WRITE(LOUT,*)
     1                   'RDKEY ERROR: Unknown material name: ', SUB(3)
                         KERR = .TRUE.
                      ENDIF
                   ENDIF
                ELSE
                   CXAXIS = 'TSURF'
                   CALL CKCRAY(SUB(2)(1:LT), NMAT, MSYM,
     1                         LOUT, NDIM, NRAY, NF, IERR)
                   IF (NF .GT. 0 .AND. .NOT. IERR) THEN
                      MXAXIS = NRAY(1)
                   ELSE
                      WRITE(LOUT,*)
     1                  'RDKEY ERROR: Unknown material name: ', SUB(2)
                      KERR = .TRUE.
                   ENDIF
                   IF (NFOUND .GE. 3) THEN
                      LT = ILASCH(SUB(3))
                      IF (SUB(3)(1:LT) .EQ. 'ARRH') THEN
                         CXAXIS = '1000/TSURF'
                      ELSE
                         CXAXIS = 'TSURF'
                      ENDIF
                   ENDIF
                ENDIF
             ELSE
                 CXAXIS = 'TSURF'
                 MXAXIS = 1
             ENDIF
           ELSE IF (SUB(1) .EQ. 'HTRN') THEN
             CXAXIS = 'HTRN'
           ELSE IF (SUB(1) .EQ. 'TAMB') THEN
             CXAXIS = 'TAMB'
           ELSE IF (SUB(1) .EQ. 'TIN') THEN
             CXAXIS = 'TIN'
           ELSE IF (SUB(1) .EQ. 'XIN') THEN
             CXAXIS = 'XIN'
             IF (NFOUND .GE. 2) THEN
               LT = ILASCH(SUB(2))
               CALL CKCRAY(SUB(2)(1:LT), KKGAS, KSYM(1,1),
     1                     LOUT, NDIM, NRAY, NF, IERR)
               IF (NF .LE. 0 .OR. IERR) THEN
                 WRITE(LOUT,*)
     1            'RDKEY ERROR: Unknown gas phase species: ', SUB(2)
                 KERR = .TRUE.
               ELSE
                 KXAXIS = NRAY(1)
               END IF
             ELSE
               WRITE(LOUT,*)
     1          ' RDKEY ERROR: XAXS XIN requires a species name'
               KERR = .TRUE.
             END IF
           ELSE IF (SUB(1) .EQ. 'PIN') THEN
             CXAXIS = 'PIN'
             IF (NFOUND .GE. 2) THEN
               LT = ILASCH(SUB(2))
               CALL CKCRAY(SUB(2), KKGAS, KSYM(1,1),
     1                     LOUT, NDIM, NRAY, NF, IERR)
               IF (NF .LE. 0 .OR. IERR) THEN
                 WRITE(LOUT,*)
     1            'RDKEY ERROR: Unknown gas phase species: ', SUB(2)
                 KERR = .TRUE.
               ELSE
                 KXAXIS = NRAY(1)
               END IF
             ELSE
               WRITE(LOUT,*)
     1          ' RDKEY ERROR: XAXS PIN requires a species name'
               KERR = .TRUE.
             END IF
           ELSE IF (SUB(1) .EQ. 'TGAS') THEN
             IF (NFOUND .GE. 2) THEN
               IF (SUB(2) .EQ. 'ARRH') THEN
                 CXAXIS = '1000/TGAS'
               ELSE
                 CXAXIS = 'TGAS'
               END IF
             ELSE
                 CXAXIS = 'TGAS'
             END IF
C
           ELSE IF (SUB(1) .EQ. 'TELEC') THEN
             IF (NFOUND .GE. 2) THEN
               IF (SUB(2) .EQ. 'EV') THEN
                 CXAXIS = 'TELEC(EV)'
               ELSE
                 CXAXIS = 'TELEC'
               END IF
             ELSE
                 CXAXIS = 'TELEC'
             END IF
           ELSE IF (SUB(1) .EQ. 'TION') THEN
             IF (NFOUND .GE. 2) THEN
               IF (SUB(2) .EQ. 'EV') THEN
                 CXAXIS = 'TION(EV)'
               ELSE
                 CXAXIS = 'TION'
               END IF
             ELSE
                 CXAXIS = 'TION'
             END IF
           ELSE IF (SUB(1) .EQ. 'IONE') THEN
             CXAXIS = 'IONE(eV)'
             IF (NFOUND .GE. 2) THEN
                LT = ILASCH(SUB(2))
                IF (SUB(2)(1:LT) .EQ. 'EV') THEN
                   WRITE(LOUT,*)
     1             'RDKEY WARNING: EV not required after IONE' 
                   MXAXIS = 1
                ELSE
                   CALL CKCRAY(SUB(2)(1:LT), NMAT, MSYM,
     1                         LOUT, NDIM, NRAY, NF, IERR)
                   IF (NF .GT. 0 .AND. .NOT. IERR) THEN
                      MXAXIS = NRAY(1)
                   ELSE
                      WRITE(LOUT,*)
     1                  'RDKEY ERROR: Unknown material name: ', SUB(2)
                      KERR = .TRUE.
                   ENDIF
                ENDIF
             ELSE
                WRITE(LOUT,*)
     1               ' RDKEY WARNING: Assuming MATERIAL1 for IONE'
                MXAXIS = 1
             END IF
C
           ELSE IF (SUB(1) .EQ. 'BIAS') THEN
             CXAXIS = 'BIAS(W)'
             IF (NFOUND .GE. 2) THEN
                LT = ILASCH(SUB(2))
                IF (SUB(2)(1:LT) .EQ. 'W') THEN
                   WRITE(LOUT,*)
     1             'RDKEY WARNING: W not required after BIAS' 
                   MXAXIS = 1
                ELSE
                   CALL CKCRAY(SUB(2)(1:LT), NMAT, MSYM,
     1                         LOUT, NDIM, NRAY, NF, IERR)
                   IF (NF .GT. 0 .AND. .NOT. IERR) THEN
                      MXAXIS = NRAY(1)
                   ELSE
                      WRITE(LOUT,*)
     1                  'RDKEY ERROR: Unknown material name: ', SUB(2)
                      KERR = .TRUE.
                   ENDIF
                ENDIF
             ELSE
                WRITE(LOUT,*)
     1               ' RDKEY WARNING: Assuming MATERIAL1 for BIAS'
                MXAXIS = 1
             END IF
C
           ELSE
             WRITE(LOUT,*)' RDKEY ERROR: Unknown option for XAXS'
             KERR = .TRUE.
           END IF
         ELSE
           CXAXIS = 'SOLUTION_#'
         ENDIF
C
C y-variable specification (other than mole frac, or rop)
C####---------------------------------------------------------------PRAM
      ELSEIF (KEY .EQ. 'PRAM') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, IERR)
 13      CONTINUE
         IF (NFOUND .LE. 0) THEN
           NFOUND = 1
           SUB(1) = 'ALL'
           GO TO 13
         ELSE
           IPSPEC = IPSPEC + 1
           IF (IPSPEC .GT. MAXPAR) THEN
              WRITE(LOUT,*)
     1            'RDKEY ERROR:...parameter list exceeds ',MAXPAR
              KERR = .TRUE.
           ENDIF
           IF (SUB(1) .EQ. 'ALL') THEN
             IPSPEC  = MAXPAR
             PSPEC(IPSPEC) = 'ALL'
           ELSE IF (SUB(1) .EQ. 'NS') THEN
             PSPEC(IPSPEC) = 'SOLUTION_#'
           ELSE IF (SUB(1) .EQ. 'EQUIV') THEN
             PSPEC(IPSPEC) = 'EQUIV_RATIO'
           ELSE IF (SUB(1) .EQ. 'PTOTAL') THEN
             PSPEC(IPSPEC) = 'PRESSURE(TORR)'
           ELSE IF (SUB(1) .EQ. 'VOLUME') THEN
             PSPEC(IPSPEC) = 'VOLUME(CM3)'
           ELSE IF (SUB(1) .EQ. 'AREA') THEN
             PSPEC(IPSPEC) = 'AREA(CM2)'
           ELSE IF (SUB(1) .EQ. 'FLRT') THEN
             PSPEC(IPSPEC) = 'FLRT(SCCM)'
           ELSE IF (SUB(1) .EQ. 'TAU') THEN
             PSPEC(IPSPEC) = 'TAU(1/S)'
           ELSE IF (SUB(1) .EQ. 'QL') THEN
             PSPEC(IPSPEC) = 'QL(WATTS)'
           ELSE IF (SUB(1) .EQ. 'POWER') THEN
             PSPEC(IPSPEC) = 'POWER(WATTS)'
           ELSE IF (SUB(1) .EQ. 'IONE') THEN
             PSPEC(IPSPEC) = 'IONE(eV)'
             IF (NFOUND .GE. 2) THEN
                LT = ILASCH(SUB(2))
                IF (SUB(2)(1:LT) .EQ. 'EV') THEN
                   WRITE(LOUT,*)
     1             'RDKEY WARNING: EV not required after IONE' 
                   MPSPEC(IPSPEC) = 1
                ELSE
                   CALL CKCRAY(SUB(2)(1:LT), NMAT, MSYM,
     1                     LOUT, NDIM, NRAY, NF, IERR)
                   IF (NF .LE. 0 .OR. IERR) THEN
                      WRITE(LOUT,*)
     1                  'RDKEY ERROR: Unknown material name: ', SUB(2)
                      KERR = .TRUE.
                   ELSE
                      MPSPEC(IPSPEC) = NRAY(1)
                   END IF
                ENDIF
             ELSE
                MPSPEC(IPSPEC) = 1
             END IF
           ELSE IF (SUB(1) .EQ. 'BIAS') THEN
             PSPEC(IPSPEC) = 'BIAS(W)'
             IF (NFOUND .GE. 2) THEN
                LT = ILASCH(SUB(2))
                IF (SUB(2)(1:LT) .EQ. 'W') THEN
                   WRITE(LOUT,*)
     1             'RDKEY WARNING: W not required after BIAS' 
                   MPSPEC(IPSPEC) = 1
                ELSE
                   CALL CKCRAY(SUB(2)(1:LT), NMAT, MSYM,
     1                     LOUT, NDIM, NRAY, NF, IERR)
                   IF (NF .LE. 0 .OR. IERR) THEN
                      WRITE(LOUT,*)
     1                  'RDKEY ERROR: Unknown material name: ', SUB(2)
                      KERR = .TRUE.
                   ELSE
                      MPSPEC(IPSPEC) = NRAY(1)
                   END IF
                ENDIF
             ELSE
                MPSPEC(IPSPEC) = 1
             END IF
           ELSE IF (SUB(1) .EQ. 'INLOSS') THEN
             PSPEC(IPSPEC) = 'INLOSS(EV)'
           ELSE IF (SUB(1) .EQ. 'TIME') THEN
             PSPEC(IPSPEC) = 'TIME(SEC)'
           ELSE IF (SUB(1) .EQ. 'TSURF') THEN
             IF (NFOUND .GE. 2) THEN
                LT = ILASCH(SUB(2))
                IF (SUB(2)(1:LT) .EQ. 'ARRH') THEN
                   PSPEC(IPSPEC) = '1000/TSURF'
                   IF (NFOUND .GE. 3) THEN
                      LT = ILASCH(SUB(3))
                      CALL CKCRAY(SUB(3)(1:LT), NMAT, MSYM,
     1                            LOUT, NDIM, NRAY, NF, IERR)
                      IF (NF .GT. 0 .AND. .NOT. IERR) THEN
                         MPSPEC(IPSPEC) = NRAY(1)
                      ELSE
                         WRITE(LOUT,*)
     1                   'RDKEY ERROR: Unknown material name: ', SUB(3)
                         KERR = .TRUE.
                      ENDIF
                   ENDIF
                ELSE
                   PSPEC(IPSPEC) = 'TSURF'
                   CALL CKCRAY(SUB(2)(1:LT), NMAT, MSYM,
     1                         LOUT, NDIM, NRAY, NF, IERR)
                   IF (NF .GT. 0 .AND. .NOT. IERR) THEN
                      MPSPEC(IPSPEC) = NRAY(1)
                   ELSE
                      WRITE(LOUT,*)
     1                  'RDKEY ERROR: Unknown material name: ', SUB(2)
                      KERR = .TRUE.
                   ENDIF
                   IF (NFOUND .GE. 3) THEN
                      LT = ILASCH(SUB(3))
                      IF (SUB(3)(1:LT) .EQ. 'ARRH') THEN
                         PSPEC(IPSPEC) = '1000/TSURF'
                      ELSE
                         PSPEC(IPSPEC) = 'TSURF'
                      ENDIF
                   ENDIF
                ENDIF
             ELSE
                 PSPEC(IPSPEC) = 'TSURF'
                 MPSPEC(IPSPEC) = 1
             ENDIF
           ELSE IF (SUB(1) .EQ. 'HTRN') THEN
             PSPEC(IPSPEC) = 'HTRN'
           ELSE IF (SUB(1) .EQ. 'TAMB') THEN
             PSPEC(IPSPEC) = 'TAMB'
           ELSE IF (SUB(1) .EQ. 'TIN') THEN
             PSPEC(IPSPEC) = 'TIN'
           ELSE IF (SUB(1) .EQ. 'XIN') THEN
             PSPEC(IPSPEC) = 'XIN'
             IF (NFOUND .GE. 2) THEN
               LT = ILASCH(SUB(2))
               CALL CKCRAY(SUB(2)(1:LT), KKGAS, KSYM(1,1),
     1                     LOUT, NDIM, NRAY, NF, IERR)
               IF (NF .LE. 0 .OR. IERR) THEN
                 WRITE(LOUT,*)
     1            'RDKEY ERROR: Unknown gas phase species: ', SUB(2)
                 KERR = .TRUE.
               ELSE
                 KPSPEC(IPSPEC) = NRAY(1)
               END IF
             ELSE
C                (plot last inlet mole fractions)
               KPSPEC(IPSPEC) = KKGAS
             END IF
           ELSE IF (SUB(1) .EQ. 'PIN') THEN
             PSPEC(IPSPEC) = 'PIN'
             IF (NFOUND .GE. 2) THEN
               LT = ILASCH(SUB(2))
               CALL CKCRAY(SUB(2), KKGAS, KSYM(1,1),
     1                     LOUT, NDIM, NRAY, NF, IERR)
               IF (NF .LE. 0 .OR. IERR) THEN
                 WRITE(LOUT,*)
     1            'RDKEY ERROR: Unknown gas phase species: ', SUB(2)
                 KERR = .TRUE.
               ELSE
                 KPSPEC(IPSPEC) = NRAY(1)
               END IF
             ELSE
C                (plot all inlet mole fractions)
               KPSPEC(IPSPEC) = 0
             END IF
           ELSE IF (SUB(1) .EQ. 'TGAS') THEN
             IF (NFOUND .GE. 2) THEN
               IF (SUB(2) .EQ. 'ARRH') THEN
                 PSPEC(IPSPEC) = '1000/TGAS'
               ELSE
                 PSPEC(IPSPEC) = 'TGAS'
               END IF
             ELSE
               PSPEC(IPSPEC) = 'TGAS'
             END IF
C
           ELSE IF (SUB(1) .EQ. 'TELEC') THEN
             IF (NFOUND .GE. 2) THEN
               IF (SUB(2) .EQ. 'EV') THEN
                 PSPEC(IPSPEC) = 'TELEC(EV)'
               ELSE
                 PSPEC(IPSPEC) = 'TELEC'
               END IF
             ELSE
                 PSPEC(IPSPEC) = 'TELEC'
             END IF
           ELSE IF (SUB(1) .EQ. 'TION') THEN
             IF (NFOUND .GE. 2) THEN
               IF (SUB(2) .EQ. 'EV') THEN
                 PSPEC(IPSPEC) = 'TION(EV)'
               ELSE
                 PSPEC(IPSPEC) = 'TION'
               END IF
             ELSE
                 PSPEC(IPSPEC) = 'TION'
             END IF
C
           ELSE
             WRITE(LOUT,*)' RDKEY ERROR: Unknown option for PRAM'
             KERR = .TRUE.
           END IF
         NEC(1) = .TRUE.
         END IF
C
C species mole fractions as y-data
C####---------------------------------------------------------------MOLS
      ELSEIF (KEY .EQ. 'MOLS') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, KDUM, NSPEC, IERR)
         KERR = KERR.OR.IERR
         IF (NSPEC .LE. 0) THEN
           KDUM(1) = 'ALL'
           NSPEC = 1
         ENDIF
         OPTION = CKCHUP(KDUM(1), 80) 
         IF (OPTION .EQ. 'ALL') THEN
           DO 16 K = 1, KKTOT
             IKSPEC = IKSPEC + 1
             KSPEC(IKSPEC) = K
 16        CONTINUE
         ELSE
           IPHASE = 0
           KSTOT = 0
           DO 170 IM = 1, NMAT
              IF (IPHASE .EQ. 0) THEN
                 CALL SKCOMP(KDUM(1)(1:16), PSYM(1,IM)(1:16), 
     1                       NPHASE(IM), IPHDUM, NT)
                 IF (IPHDUM.NE.0) THEN
                    IPHASE = IPHDUM 
                    IMAT = IM
                    KSMAT = KSTOT
                 ENDIF
                 KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
              ENDIF
 170       CONTINUE
           IF (NT .GT. 0) THEN
             OPTION = CKCHUP(KDUM(2), 80)
             IF (NSPEC .EQ. 1 .OR. OPTION .EQ. 'ALL') THEN
               DO 17 KM = KFIRST(IPHASE,IMAT), KLAST(IPHASE,IMAT)
                  IKSPEC = IKSPEC + 1
                  K = KM + KSMAT
                  KSPEC(IKSPEC) = K
 17            CONTINUE
             ELSE
               DO 18 N = 2, MIN(NDIM, NSPEC)
C                Search iphase for matching species names
                 KF = KFIRST(IPHASE,IMAT)
                 CALL SKSNUM (KDUM(N), 0, LOUT, KSYM(KF,IMAT),
     1                       KKPHAS(IPHASE,IMAT), PSYM(IPHASE,IMAT), 1,
     2                       KKPHAS(IPHASE,IMAT), KNUM, NT, NVAL,
     3                       VALUE, IERR)
                 IF (KNUM .GT. 0 .AND. KNUM .LE. KKPHAS(IPHASE,IMAT)) 
     1           THEN
                   KNUM = KNUM + KSMAT
                   IKSPEC = IKSPEC + 1
                   KSPEC(IKSPEC) = KNUM+KFIRST(IPHASE,IMAT)-1
                 ELSE
                   LT = ILASCH(KDUM(N))
                   WRITE(LOUT,*) 'RDKEY ERROR: ',KDUM(N)(1:LT),
     1                           ' is not a species name in the',
     2                           PSYM(IPHASE,IMAT), ' phase'
                   KERR = .TRUE.
                 END IF
 18           CONTINUE
             END IF
           ELSE
C            Search for individual species names
             DO 20 N = 1, MIN(NDIM, NSPEC)
               KNUM = 0
               KSTOT = 0
               DO 19 IM = 1, NMAT
                  IF (KNUM .EQ. 0) THEN
                     CALL SKPCMP (KDUM(N), KSYM(1,IM), KKM(IM), 
     1                            PSYM(1,IM), NPHASE(IM), KKPHAS(1,IM),
     2                            KNUMD, NT)
                     IF (KNUMD .GT. 0) THEN
                        KNUM = KNUMD+KSTOT
                        IMAT = IM
                     ENDIF
                  ENDIF
                  KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
 19            CONTINUE
               IF (KNUM .GT. 0 .AND. KNUM .LE. KKTOT) THEN
                 IKSPEC = IKSPEC + 1
                 KSPEC(IKSPEC) = KNUM
               ELSE
                 LT = ILASCH(KDUM(N))
                 WRITE(LOUT,*) 'RDKEY ERROR: ',KDUM(N)(1:LT),
     1                         ' is not a species name'
                 KERR = .TRUE.
               END IF
   20        CONTINUE
           END IF
         END IF
         NEC(1) = .TRUE.
C
C species number densities as y-data
C####---------------------------------------------------------------NUMD
      ELSEIF (KEY .EQ. 'NUMD') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, KDUM, NSPEC, IERR)
         KERR = KERR.OR.IERR
         IF (NSPEC .LE. 0) THEN
           KDUM(1) = 'ALL'
           NSPEC = 1
         ENDIF
         OPTION = CKCHUP (KDUM(1), 80)
         IF (OPTION .EQ. 'ALL') THEN
           DO 26 K = 1, KKTOT
             IKSPN = IKSPN + 1
             KSPN(IKSPN) = K
 26        CONTINUE
         ELSE
           KSTOT = 0
           IPHASE = 0
           DO 260 IM = 1, NMAT
              IF (IPHASE .EQ. 0) THEN
                 CALL SKCOMP(KDUM(1)(1:16), PSYM(1,IM)(1:16), 
     1                       NPHASE(IM), IPHDUM, NT)
                 IF (IPHDUM .NE. 0) THEN
                    IPHASE = IPHDUM
                    IMAT = IM
                    KSMAT = KSTOT
                 ENDIF
                 KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
              ENDIF
 260       CONTINUE
           IF (NT .GT. 0) THEN
               OPTION = CKCHUP(KDUM(2), 80)
             IF (NSPEC .EQ. 1 .OR. OPTION .EQ. 'ALL') THEN
               DO 27 KM = KFIRST(IPHASE,IMAT), KLAST(IPHASE,IMAT)
                  IKSPN = IKSPN + 1
                  KSPN(IKSPN) = KM + KSMAT
 27            CONTINUE
             ELSE
               DO 28 N = 2, MIN(NDIM, NSPEC)
C                Search iphase for matching species names
                 KF = KFIRST(IPHASE,IMAT)
                 CALL SKSNUM (KDUM(N), 0, LOUT, KSYM(KF,IMAT),
     1                       KKPHAS(IPHASE,IMAT), PSYM(IPHASE,IMAT), 1,
     1                       KKPHAS(IPHASE,IMAT), KNUM, NT, NVAL,
     1                       VALUE, IERR)
                 IF (KNUM .GT. 0 .AND. KNUM .LE. KKPHAS(IPHASE,IMAT)) 
     1           THEN
                   IKSPN = IKSPN + 1
                   KSPN(IKSPN) = KNUM+KFIRST(IPHASE,IMAT)-1 + KSMAT
                 ELSE
                   LT = ILASCH(KDUM(N))
                   WRITE(LOUT,*) 'RDKEY ERROR: ',KDUM(N)(1:LT),
     1                           ' is not a species name in the',
     2                           PSYM(IPHASE,IMAT), ' phase'
                   KERR = .TRUE.
                 END IF
 28            CONTINUE
             END IF
           ELSE
C            Search for individual species names
             DO 30 N = 1, MIN(NDIM, NSPEC)
               KNUM = 0
               KSTOT =0
               DO 301 IM = 1, NMAT
                  IF (KNUM .EQ. 0) THEN
                     CALL SKPCMP (KDUM(N), KSYM(1,IM), KKM(IM), 
     1                            PSYM(1,IM),NPHASE(IM), KKPHAS(1,IM),
     2                            KNUMD, NT)
                     IF (KNUMD .NE. 0) THEN
                        KNUM = KNUMD + KSTOT
                        IMAT = IM
                     ENDIF
                     KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
                  ENDIF
 301           CONTINUE
               IF (KNUM .GT. 0 .AND. KNUM .LE. KKM(IMAT)) THEN
                 IKSPN = IKSPN + 1
                 KSPN(IKSPN) = KNUM 
               ELSE
                 LT = ILASCH(KDUM(N))
                 WRITE(LOUT,*) 'RDKEY ERROR: ',KDUM(N)(1:LT),
     1                         ' is not a species name'
                 KERR = .TRUE.
               END IF
   30        CONTINUE
           END IF
         END IF
         NEC(1) = .TRUE.
C
C species outlet flow rates (sccm) as y-data
C####---------------------------------------------------------------OFRS
      ELSEIF (KEY .EQ. 'OFRS') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, KDUM, NSPEC, IERR)
         KERR = KERR.OR.IERR
         IF (NSPEC .LE. 0) THEN
           KDUM(1) = 'ALL'
           NSPEC = 1
         ENDIF
         OPTION = CKCHUP (KDUM(1), 80)
         IF (OPTION .EQ. 'ALL') THEN
           DO 32 K = 1, KKTOT
             IKSPF = IKSPF + 1
             KSPF(IKSPF) = K
 32       CONTINUE
         ELSE
           KSTOT = 0
           IPHASE = 0
           DO 320 IM = 1, NMAT
              IF (IPHASE .EQ. 0) THEN
                 CALL SKCOMP(KDUM(1)(1:16), PSYM(1,IM)(1:16), 
     1                       NPHASE(IM), IPHDUM, NT)
                 IF (IPHDUM .NE. 0) THEN
                    IPHASE = IPHDUM
                    IMAT = IM
                    KSMAT = KSTOT
                 ENDIF
                 KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
              ENDIF
 320       CONTINUE
           IF (NT .GT. 0) THEN
               OPTION = CKCHUP(KDUM(2), 80)
             IF (NSPEC .EQ. 1 .OR. OPTION .EQ. 'ALL') THEN
               DO 33 KM = KFIRST(IPHASE,IMAT), KLAST(IPHASE,IMAT)
                  IKSPF = IKSPF + 1
                  KSPF(IKSPF) = KM + KSMAT
 33            CONTINUE
             ELSE
               DO 35 N = 2, MIN(NDIM, NSPEC)
C                Search iphase for matching species names
                 KF = KFIRST(IPHASE,IMAT)
                 CALL SKSNUM (KDUM(N), 0, LOUT, KSYM(KF,IMAT),
     1                       KKPHAS(IPHASE,IMAT), PSYM(IPHASE,IMAT), 1,
     1                       KKPHAS(IPHASE,IMAT), KNUM, NT, NVAL,
     1                       VALUE, IERR)
                 IF (KNUM .GT. 0 .AND. KNUM .LE. KKPHAS(IPHASE,IMAT)) 
     1           THEN
                   IKSPF = IKSPF + 1
                   KSPF(IKSPF) = KNUM+KFIRST(IPHASE,IMAT)-1 + KSMAT
                 ELSE
                   LT = ILASCH(KDUM(N))
                   WRITE(LOUT,*) 'RDKEY ERROR: ',KDUM(N)(1:LT),
     1                           ' is not a species name in the',
     2                           PSYM(IPHASE,IMAT), ' phase'
                   KERR = .TRUE.
                 END IF
 35           CONTINUE
             END IF
           ELSE
C            Search for individual species names
             DO 38 N = 1, MIN(NDIM, NSPEC)
               KNUM = 0
               KSTOT =0
               DO 381 IM = 1, NMAT
                  IF (KNUM .EQ. 0) THEN
                     CALL SKPCMP (KDUM(N), KSYM(1,IM), KKM(IM), 
     1                            PSYM(1,IM),NPHASE(IM), KKPHAS(1,IM),
     2                            KNUMD, NT)
                     IF (KNUMD .NE. 0) THEN
                        KNUM = KNUMD + KSTOT
                        IMAT = IM
                     ENDIF
                     KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
                  ENDIF
 381           CONTINUE
               IF (KNUM .GT. 0 .AND. KNUM .LE. KKM(IMAT)) THEN
                 IKSPF = IKSPF + 1
                 KSPF(IKSPF) = KNUM 
               ELSE
                 LT = ILASCH(KDUM(N))
                 WRITE(LOUT,*) 'RDKEY ERROR: ',KDUM(N)(1:LT),
     1                         ' is not a species name'
                 KERR = .TRUE.
               END IF
 38         CONTINUE
          ENDIF
        END IF
        NEC(1) = .TRUE.
C####---------------------------------------------------------------GRTH
      ELSEIF (KEY .EQ. 'GRTH') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, KDUM, NSPEC, IERR)
         KERR = KERR.OR.IERR
         IF (NSPEC .LE. 0) THEN
           KDUM(1) = 'ALL'
           NSPEC = 1
         END IF
         IFOUND = 0
         DO 125 IWORD = 1, NSPEC
            OPTION = CKCHUP(KDUM(IWORD), 80)
            IF (OPTION .EQ. 'ALL') THEN
              LGRALL = .TRUE.
              IFOUND = 1
            ELSE
              ISTOT = 0
              INUM = 0
              DO 40 IM = 1, NMAT
                 IF (INUM .EQ. 0) THEN
                    CALL SKCOMP (KDUM(IWORD), PSYM(1,IM), NPHASE(IM),
     1                           INUMD, NT)
                    IF (INUMD .NE. 0) THEN
                       INUM = INUMD + ISTOT
                       IMAT = IM
                       IFOUND = 1
                    ENDIF
                    ISTOT = ISTOT + NNSURF(IM) + NNBULK(IM)
                 ENDIF
 40           CONTINUE
              IF (NT .GT. 0) THEN
                IFOUND = 1
                IGPHAS = IGPHAS + 1
                KGPHAS(IGPHAS) = INUM 
              ELSE IF (OPTION .EQ. 'ANG/HR') THEN
                GRUNIT = 'ANG/HR'
              ELSE IF (OPTION .EQ. 'ANG/MIN') THEN
                GRUNIT = 'ANG/MIN'
              ELSE IF (OPTION .EQ. 'ANG/SEC') THEN
                GRUNIT = 'ANG/SEC'
              ELSE IF (OPTION .EQ. 'UM/HR') THEN
                GRUNIT = 'UM/HR'
              ELSE IF (OPTION .EQ. 'UM/MIN') THEN
                GRUNIT = 'UM/MIN'
              ELSE IF (OPTION .EQ. 'UM/SEC') THEN
                GRUNIT = 'UM/SEC'
              ELSE IF (OPTION .EQ. 'MOLE/CM**2*SEC') THEN
                GRUNIT = 'MOLE/CM2-S'
              ELSE IF (OPTION .EQ. 'GM/CM**2*SEC') THEN
                GRUNIT = 'G/CM2-S'
              ELSE IF (OPTION .EQ. 'ML/SEC') THEN
                GRUNIT = 'ML/SEC'
              ELSE IF (OPTION .EQ. 'ANG/MIN') THEN
                GRUNIT = 'ANG/MIN'
              ELSE IF (OPTION .EQ. 'CM/SEC') THEN
                GRUNIT = 'CM/SEC'
              ELSE
                IF (IFOUND .LE. 0) THEN
                   WRITE(LOUT,*)'RDKEY ERROR: Phase selected for ',
     1           'growth rate plot, ',KDUM(IWORD),' is not a bulk phase'
                ELSE
                   WRITE(LOUT,*)
     1           'RDKEY ERROR: Unknown option for GRTH keyword: ',
     2            KDUM(IWORD)
                ENDIF
                KERR = .TRUE.
              ENDIF
            END IF
125      CONTINUE
         NEC(1) = .TRUE.
C
C rates of progress as y-axis data
C####---------------------------------------------------------------ROP
      ELSEIF (KEY .EQ. 'RATE') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, KDUM, NSPEC, IERR)
         KERR = KERR.OR.IERR
         IF (NSPEC .LE. 0) THEN
            KDUM(1) = 'ALL'
            NSPEC = 1
         ENDIF
         OPTION = CKCHUP(KDUM(1), 80)
         IF (OPTION .EQ. 'ALL') THEN
            DO 128 K = 1, KKTOT
               IRSPEC = IRSPEC + 1
               KRSPEC(IRSPEC) = K
 128        CONTINUE
         ELSE
            KSTOT = 0
            IPHASE = 0
            DO 129 IM = 1, NMAT
               IF (IPHASE .EQ. 0) THEN
                  CALL SKCOMP (KDUM(1)(1:16), PSYM(1,IM)(1:16), 
     1                         NPHASE(IM), IPHDUM, NT)
                  IF (IPHDUM .NE. 0) THEN
                     IPHASE = IPHDUM
                     IMAT = IM
                     KSMAT = KSTOT
                  ENDIF
                  KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
               ENDIF
 129        CONTINUE
            IF (NT .GT. 0) THEN
               OPTION = CKCHUP(KDUM(2), 80)
               IF (NSPEC .EQ. 1 .OR. OPTION .EQ. 'ALL') THEN
                  DO 130 K = KFIRST(IPHASE,IMAT), KLAST(IPHASE,IMAT)
                     IRSPEC = IRSPEC + 1
                     KRSPEC(IRSPEC) = K + KSMAT
 130              CONTINUE
               ELSE
                  DO 132 N = 2, MIN(NDIM, NSPEC)
                     KF = KFIRST(IPHASE,IMAT)
                     CALL SKSNUM (KDUM(N), 0, LOUT, KSYM(KF,IMAT),
     1                           KKPHAS(IPHASE,IMAT),PSYM(IPHASE,IMAT),
     2                           1, KKPHAS(IPHASE,IMAT), KNUM, NT, 
     3                           NVAL, VALUE, IERR)
                     IF (KNUM.GT.0 .AND. KNUM.LE.KKPHAS(IPHASE,IMAT))
     1               THEN
                        IRSPEC = IRSPEC + 1
                        KRSPEC(IRSPEC) = KNUM + KKPHAS(IPHASE,IMAT)-1
     1                                   + KSMAT
                     ELSE
                        LT = ILASCH(KDUM(N))
                        WRITE(LOUT,*) 'RDKEY ERROR: ', KDUM(N)(1:LT),
     1                                ' is not a species name in the ',
     2                                PSYM(IPHASE,IMAT), ' phase'
                        KERR = .TRUE.
                     ENDIF
 132              CONTINUE
               ENDIF
            ELSE
               DO 134 N = 1, MIN(NDIM, NSPEC)
                  KNUM = 0
                  KSTOT = 0
                  DO 133 IM = 1, NMAT
                     IF (KNUM .EQ. 0) THEN
                        CALL SKPCMP (KDUM(N), KSYM(1,IM), KKM(IM), 
     1                               PSYM(1,IM), NPHASE(IM), 
     2                               KKPHAS(1,IM), KNUMD, NT)
                        IF (KNUMD .NE. 0) THEN
                           KNUM = KNUMD + KSTOT
                           IMAT = IM
                        ENDIF
                        KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
                     ENDIF
 133              CONTINUE
                  IF (KNUM .GT. 0 .AND. KNUM .LE. KKM(IMAT)) THEN
                     IRSPEC = IRSPEC + 1
                     KRSPEC(IRSPEC) = KNUM
                  ELSE
                     LT = ILASCH(KDUM(N))
                     WRITE(LOUT,*) 'RDKEY ERROR: ',KDUM(N)(1:LT),
     1                             ' is not a species name'
                     KERR = .TRUE.
                  ENDIF
 134           CONTINUE
            ENDIF
         ENDIF
         NEC(1) = .TRUE.
C####---------------------------------------------------------------CURR
      ELSEIF (KEY .EQ. 'CURR') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, KDUM, NSPEC, IERR)
         KERR = KERR.OR.IERR
         IF (NSPEC .LE. 0) THEN
           KDUM(1) = 'ALL'
           NSPEC = 1
         END IF
         IFOUND = 0
         DO 138 IWORD = 1, NSPEC
            OPTION = CKCHUP(KDUM(IWORD), 80)
            IF (OPTION .EQ. 'ALL') THEN
              LCURALL = .TRUE.
              IFOUND = 1
            ELSEIF (OPTION .EQ. 'TOTAL') THEN
              LCURTOT = .TRUE.
              IFOUND = 1
            ELSE
              CALL SKCOMP (KDUM(IWORD), KSYM(1,1), KKM(1), KNUM, NT)
              IF (NT .GT. 0) THEN
                DO 137 KI = 1, KKION
                   IF (KNUM.EQ.KION(KI)) THEN
                      IFOUND = 1
                      ICURR = ICURR+1
                      KCURR(ICURR) = KI
                   ENDIF
 137            CONTINUE
              ELSE
                IF (IFOUND .LE. 0) THEN
                   WRITE(LOUT,*)'RDKEY ERROR: species selected for ',
     1           'ion current plot, ',KDUM(IWORD),' is not an ion'
                ELSE
                   WRITE(LOUT,*)
     1           'RDKEY ERROR: Unknown option for CURR keyword: ',
     2            KDUM(IWORD)
                ENDIF
                KERR = .TRUE.
              ENDIF
            END IF
 138     CONTINUE
         NEC(1) = .TRUE.
C####---------------------------------------------------------------YLDI
      ELSEIF (KEY .EQ. 'YLDI') THEN
         CALL CKSUBS (LINE, LOUT, NDIM, KDUM, NSPEC, IERR)
         KERR = KERR.OR.IERR
         IF (NSPEC .LE. 0) THEN
           KDUM(1) = 'ALL'
           NSPEC = 1
         END IF
         IFOUND = 0
         DO 146 IWORD = 1, NSPEC
            OPTION = CKCHUP(KDUM(IWORD), 80)
            IF (OPTION .EQ. 'ALL') THEN
              LYLDALL = .TRUE.
              IFOUND = 1
            ELSEIF (OPTION .EQ. 'TOTAL') THEN
              LYLDTOT = .TRUE.
              IFOUND = 1
            ELSE
              CALL SKCOMP (KDUM(IWORD), KSYM(1,1), KKM(1), KNUM, NT)
              IF (NT .GT. 0) THEN
                DO 144 KI = 1, KKION
                   IF (KNUM.EQ.KION(KI)) THEN
                      IFOUND = 1
                      IYLD = IYLD+1
                      KYLD(IYLD) = KI
                   ENDIF
 144            CONTINUE
              ELSE
                IF (IFOUND .LE. 0) THEN
                   WRITE(LOUT,*)'RDKEY ERROR: species selected for ',
     1           'ion yield plot, ',KDUM(IWORD),' is not an ion'
                ELSE
                   WRITE(LOUT,*)
     1           'RDKEY ERROR: Unknown option for YLDI keyword: ',
     2            KDUM(IWORD)
                ENDIF
                KERR = .TRUE.
              ENDIF
            END IF
 146     CONTINUE
         NEC(1) = .TRUE.
C####---------------------------------------------------------------IPSR
      ELSEIF (KEY .EQ. 'IPSR') THEN
         CALL CKXNUM(LINE, 1, LOUT, NVAL, VALUE(1), IERR)
         IPSR = INT(VALUE(1))
         IF (IPSR .LE. 0 .OR. IPSR .GT. NUMPSR) THEN
          WRITE(LOUT,*)' RDKEY ERROR: IPSR is out of solution bounds:',
     1                   IPSR
           KERR = .TRUE.
         END IF
         KERR = KERR.OR.IERR
      ENDIF
C
      GO TO 10
C####-------------------------------------------------------Return Block
C
C    Return from here
C
 500  CONTINUE
C
C  Check for necessary conditions
C
      IF (.NOT.NEC(1)) THEN
         WRITE(LOUT,*) 'RDKEY ERROR: No plots were specified to be made'
         WRITE(LOUT,*) '             Specify one of MOLS, NUMD, RATE,'
         WRITE(LOUT,*) '             CURR, GRTH, YLDI, or PRAM'
         KERR = .TRUE.
      END IF
C      IF (KERR) STOP
      RETURN
C------------------------------------------------------------Error Block
C
C  Read Errors
C
 600  CONTINUE
      WRITE (LOUT,*) 'RDKEY: ERROR ... error reading keyword file'
      KERR = .TRUE.
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE EIFIND (KK, KNAM, KCHG, ICKWRK, RCKWRK, KKION,
     1                   KION, KEL)
C
C  This routine locates the electron and positive ion species and
C  returns KION(KKION), a pointer to the ions.
C  e.g. if KION(1) = K, the Kth gas species is the first ion encountered
C  KKION  is the total number of positive ions
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ICKWRK(*), RCKWRK(*), KION(*), KCHG(*)
C
      CHARACTER KNAM(*)*(*),  IST*16
C
      CALL CKCHRG (ICKWRK, RCKWRK, KCHG)
      IST = 'E'
      CALL CKCOMP (IST, KNAM, KK, KEL)
      KKION = 0
      DO 100 J = 1, KK
         IF (KCHG(J).NE.0 .AND. J .NE. KEL) THEN
            KKION = KKION + 1
            KION(KKION) = J
         ENDIF
100   CONTINUE
C
      RETURN
      END

