C CVS: $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:32 $
      SUBROUTINE PSRABS
C///////////////////////////////////////////////////////////////////
C
C                 AURORA: PERFECTLY STIRRED REACTOR WITH SURFACE
C                         REACTIONS AND PLASMA CHEMISTRY
C
C     WRITTEN BY:
C
C         ELLEN MEEKS, ROBERT J. KEE AND JOSEPH F. GRCAR
C         THERMAL AND PLASMA PROCESSES DEPARTMENT
C         MAIL STOP 9042
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94551-0969
C     AND
C         HARRY K. MOFFAT
C         SURFACE SCIENCE PROCESSING DIVISION
C         SANDIA NATIONAL LABORATORIES
C         ALBUQUERQUE, NM 87125
C     AND
C         PETER GLARBORG
C         LABORATORY FOR HEATING AND AIR CONDITIONING
C         TECHNICAL UNIVERSITY OF DENMARK
C         2800 LYNGBY
C         DENMARK
C
C//////////////////////////////////////////////////////////////////////
C     VERSION 6.17
C     Changes for 6.17 (E. Meeks, 98/4/15):
C     1. Fixed bug#162:  Replaced call to SKHMS in PSRFUN with a call
C        to SKHML to avoid any chance of divide by zero for surface
C        species that my have zero molecular weights.  
C     2. Modified the enthalpy term by  /WT(K,1) in SUM1 
C        definition in PSRFUN (SUM1 only goes over gas-phase species,
C        for which the molecular weights must be nonzero), so that 
C        term is correct with molar enthalpies.
C     3. Modified the enthalpy term in SUM2 by removing the * WT(K,1),
C        so that term is correct with molar enthalpies.
C     4. Removed * WTS from definition of SUM3 in PSRFUN and removed
C        definitions of WTS, since this is no longer required.
C     5. Removed * WT(KEL,1) from F(NT) correction in the case that
C        KEL is nonzero.
C     6. Corrected comments that give units of HIN as [erg/mole] to
C        reflect actual units of [erg/g].  (HIN is defined through a 
C        call to CKHMS in PSRDRV)
C     Changes for 6.16 (E. Meeks, 98/1/23):
C     1. Fixed bug#145: Added definition of K=KM+KSTOT to correct
C        treatment for multiple materials in loop 3000 in 
C        subroutine PSRSEN.
C     2. Fixed bug#144: Moved 500 continuation statmenet in subroutine
C        PSRTWO to above TWSET  initializations, to correct case
C        when printing level is set to 0.
C     3. Fixed bug#146: Corrected dimensioning of XIN in subroutine
C        PSRKEY to be (KKTOT,MAXPSR), and corrected indexing 
C        throughout this routine.
C     4. Fixed bug#147: overwriting of work space occuring between
C        variables defined in top-level driver and those in PSPNT 
C        routine.  Moved allocation of reals and logicals: 
C        P, PA, V, Q, AREA, T, TIN, FLRT, TAU, HTRN, TAMBNT, POWR, 
C        TEIN, TE, GFAC, SFAC, ESHTH, TSURF, TIONP, LENRGY, LFLRT, 
C        LHTRN, LENRGE, LTION, LELSH, LRFSH, LWALHB, and LTDIFF 
C        down to PSPNT.  Removed allocation of integer space for 
C        NATJ in PSRUN, since NATJ is a scalar.
C     5. Updated comments describing arguments in PSPNT. 
C     Changes for 6.15 (E. Meeks, 97/12/11):
C     1. Fixed bug#136: Added two more BACKSPACE commands before 
C        overwriting recover file (unit=LRECOV) with new solution.
C     2. Updated VERSNN to agree with aurora.f version number.
C     Changes for 6.14 (E. Meeks, 97/12/1):
C     1. Fixed bug#135: removed space in line 3357 in PSRFUN s.t. 
C        line is <80 chars.
C     Changes for 6.13 (E. Meeks, 97/11/07):
C     1. Added pointer HIM in PSPNT, /PSPSPS/, and call to PSRDRV to
C        allow correct definition of species enthalpies for multiple 
C        materials. 
C     2. Passed HIM array into PSRTWO, PSPRNT, PSRFUN, PSSDOT,
C        PSFIXJ, PSRJAC, PSRSEN, and PSRROP. Dimensioned (KKMAX,NMAT).
C     3. Added loop 3180 in PSRFUN to define HIM for all materials.
C     4. In PSSDOT, define HIM coming in as ENRGI and use as
C        ion energy array. Eliminate SCRTCH(1,2) array from call list.
C     5. In PSRFUN, replace use of H array by ROP for scratch array
C        in call to PSSDOT near loop 4050.
C     6. Eliminate SCRTCH(1,2) from call list to PSRFUN for H array
C        which was replaced by HIM array.
C     7. Remove redundant assignment of TEMP(3)=TIONP after loop 2400
C        in PSRFUN.
C     Changes for 6.12 (E. Meeks, 97/10/29):
C     1. Set ABSOL=RELAT=SQRT(D1MACH(4)), rather than 
C        SQRT(2.*D1MACH(4)) to be compatible with old formulation.
C        Result on SGI is now 1.4E-8 (was 2.1E-8).
C     2. Put double-precision change blocks around above, and
C        add single-precision change block using R1MACH(4).
C     Changes for 6.11 (E. Meeks, 97/10/28):
C     1. Fixed bug#102: removed unused variables: SMALL in PSRDRV;
C        ADAP, GRAD, CURV, ADAFLR in PSRTWO; EWALL, TIMRF, VWALL,
C        VWBAR, AJTOT in PSPRNT and PSRSEN; DGESL, SGESL, XDGESL,
C        XSGESL in PSFIXJ; MASLOS, LOWSDN in PSSDOT.
C     2. Fixed bug#021: Removed non-standard f77 code; changed two 
C        occurences of do/enddo to numbered do loops; Replaced names
C        that were longer than 6 chars: LRFSHTH=>LRFTH, 
C        TAMBIENT=>TAMBNT,VERSION=>VERSNN,LASTCSK=>LSTCSK,
C        LASTISK=>LSTISK,LASTRSK=>LSTRSK,DEPRATE=>DEPRAT,
C        SPUTTER=>SPUTTR,TRNFILE=>TRNFIL,DYETEDT=>DYTEDT,
C        LEQCONT=>LEQCNT,LEQPRNT=>LEQPRN,MATETCH=>MTETCH, 
C        CONDBST=>CNDBST.
C     3. Fixed bug#096c: Replaced ABSOL and RELAT definitions to use
C        D1MACH(4) instead of calculated U.
C     Changes for 6.10 (E. Meeks, 97/08/06):
C     1. Expanded options for multiple PSRs.  In PSRKEY allow user to
C        input last numeric value as digit indicating PSR for parameters
C        that may change from one PSR to the next:  ESHTH, EIONSH, 
C        BIASPW, TIONP, LELSH, LRFSHTH, LIONP, and LWALHB.
C     2. Allocate additional workspace and dimension the above variables
C        in top level routine for above variables by *MAXPSR.
C     3. Address action#062 to allow AFRAC to change between PSRs.
C     4. New error statements in PSRKEY when IMAT is not within 1-MAXMAT
C        to avoid illegal subscripting.
C     Changes for 6.9 (E. Meeks, 97/08/05):
C     1. Added new keyword "BPWR" to allow user to specify the BIAS 
C        power applied to a material boundary.  In this case, the ion 
C        energy used in the kinetics will be approximated by the bias 
C        power divided by the total ion current density multiplied by 
C        the surface area of the material. 
C     2. Changed call statements of subroutines to allow BIASPW(NMAT) 
C        and AREA to be passed to PSSDOT wherever it is called.
C     3. Added logic in PSSDOT to determine when BIASPW option is turned
C        on and to calculate the average ion energy in this case.
C     4. Added initialization of BIASPW = 0.0 to PSRKEY and error checks
C        when this option coincides with other methods of specifying ion
C        energy.
C     5. Added writing of BIASPW to LSAVE and LRECOV units for KEL.NE.0,
C        although reading restart is unchanged since this info is not 
C        used on restart.
C     6. Move VERSN and VERSION definition down to PSR where they are
C        declared and read in.
C     7. Get rid of pesky warning messages in PSRSEN when mole fracs
C        are zero; write warning once only.
C     Changes for 6.8 (E. Meeks, 97/07/24):
C     1. Fix comment line (#59) longer than 72 characters. 
C     2. Fix Bug #012: remove extra spaces in logical operators
C        in PSRTWO (2 in loop 710), in PSRROP (in loop 95), in
C        PSRFUN (in loop 1800), in PSRKEY) (after loop 1300), in
C        PSRDRV (2 in loop 7100).
C     Changes for 6.7 (97/07/23):
C     1. Fix Bug #053:  Add call to CKYTX converting S(NY) to X for
C        IPSR+1 just above loop 8200 in PSRDRV.  Corrects inlet mole
C        fractions for multiple PSRs.
C     Changes for 6.6 (97/06/09):
C     1. Fix Bug #039:  problem with flow rate not getting defined on
C        continuations when SCCM keyword is used.
C        a) Added LSCCM to PSRKEY call statement
C        b) Added LSCCM to logical declarations in PSRDRV
C        c) Renamed SCCM to SCCMIN within PSRKEY routine to keep name 
C           same in call and avoid confusion with SCCM used in PSRDRV.
C        d) Added SCCMIN to PSRKEY call statement
C        e) Added initialization of SCCMIN=0.0 in PSRKEY for cold start.
C     Changes for 6.5 (97/05/30):
C     1. Fuxed Bug #008:  In combustion calculations with FUEL and OXID
C        specifications, required user to use all elements in chemkin
C        mechanism in both reactants and products.  Fix allows users to
C        include species/elements in mechanism not used in the problem
C        being solved in AURORA.  The specific changes are as follows:
C        a)  Add pointer IMU in PSPNT for MM long array.
C        b)  Add IMU to common block /PSPSPS/ (two places).
C        c)  Add array I(IMU) => MUSD(MM) to call list of PSRDRV.
C        d)  Add MUSD to call list for PSRKEY; Dimension MUSD(MM).
C        e)  Initialize MUSD to zeros (loop 1920) and add loop 1950 to 
C            determine which elements are actually used, with MMUSED as 
C            counter, fill MUSD array.
C        f)  Change loop 2000 for MU=1 to MMUSED elements.  Change 
C            appropriate indices of matrix elements using MU.
C        g)  Change MM => MMUSED in calls to DGEFA, SGEFA, DGESL, SGESL 
C            and correct error check to be N .NE. MMUSED.
C        h)  Improve error messages written to user.
C     Changes for 6.4 (97/05/30)
C     1. Fixed bug #007:  flow rates not recalculated on continuations 
C        when composition changes and SCCM is used.  Problem was with 
C        LSCCM being reinitialized to .FALSE. on continuation in top of
C        PSRKEY.  Removed this initialization line.
C     Changes for 6.3 (97/05/30)
C     1. Fixed Bug #020:  Missing argument NIK in call to PSRJAC from
C        PSFIXJ.
C     Changes for 6.2 (97/05/01)
C     1. update 97/03/01; make new main "driver" program to set up
C        arrays, to satisfy Reaction Design requirement to provide
C        object files instead of source code,
C     2. do not write SDOTI array if IISUR(IM).LE.0
C/////////////////////////////////////////////////////////////////////
      END
C
      SUBROUTINE PSRUN (LIN, LOUT, LINKCK, LINKSK, LREAD, LSAVE,
     1                  LRECOV, LOUTSH, LENLWK, L, LENIWK, I, LENRWK, 
     2                  R, LENCWK, C, MAXPSR, MAXMAT)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      DIMENSION I(LENIWK), R(LENRWK)
      LOGICAL L(LENLWK), KERR
      CHARACTER C(LENCWK)*16
C
      KERR = .FALSE.
C     reserve some real array space
      NUSED= 0
      IF (LENRWK .LT. NUSED) THEN
         WRITE (LOUT, *) 'PSRUN ERROR: RWORK must be at least ', NUSED
         KERR = .TRUE.
      ENDIF
C     reserve some integer array space
      INSP = 1
      IKKS = INSP + MAXMAT
      IKKB = IKKS+ MAXMAT
      IKK  = IKKB+ MAXMAT
      INPHA= IKK + MAXMAT
      INNS = INPHA+ MAXMAT
      INFS = INNS+ MAXMAT
      INLS = INFS+ MAXMAT
      INNB = INLS+ MAXMAT
      INFB = INNB+ MAXMAT
      INLB = INFB+ MAXMAT
      IIIS = INLB+ MAXMAT
      IMRS = IIIS+ MAXMAT
      IMIS = IMRS+ MAXMAT
      IMCS = IMIS+ MAXMAT
      IMPH = IMCS+ MAXMAT
      IINC = IMPH + MAXMAT
      INIS = IINC + MAXMAT
      INBHM= INIS + MAXMAT 
      IUSED = INBHM+MAXMAT - 1
      IF (LENIWK .LT. IUSED) THEN
         WRITE (LOUT, *) 'PSRUN ERROR: IWORK must be at least ', IUSED
         KERR = .TRUE.
      ENDIF
C     reserve some logical array space
      LUSED = 0
      IF (LENLWK .LT. LUSED) THEN
         WRITE (LOUT, *) 'PSRUN ERROR: LWORK must be at least ', LUSED
         KERR = .TRUE.
      ENDIF
      IF (KERR) RETURN
C
      CALL PSR (LIN, LOUT, LOUTSH, LINKCK, LINKSK, LREAD, LSAVE, 
     1          LRECOV,
     1          LENLWK, L, LENIWK, I, LENRWK, R, LENCWK, C, MAXPSR,
     2          MAXMAT, NUSED,
     6          I(INSP), I(IKKS), I(IKKB), I(IKK),  I(INPHA),
     7          I(INNS), I(INFS), I(INLS), I(INNB), I(INFB),
     8          I(INLB), I(IIIS), I(IMRS), I(IMIS), I(IMCS), 
     9          I(IMPH), I(IINC),  I(INIS), I(INBHM),IUSED, LUSED) 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PSR (LIN, LOUT, LOUTSH, LINKCK, LINKSK, LREAD, LSAVE, 
     1                LRECOV,
     1                LENLWK, L, LENIWK, I, LENRWK, R, LENCWK, C,
     2                MAXPSR, MAXMAT, NUSED, NSPHCH, KKSURF,
     5                KKBULK,KK, NPHASE, NNSURF, NFSURF, NLSURF, NNBULK,
     6                NFBULK, NLBULK, IISUR, IMRSK, IMISK, IMCSK, 
     7                IMPH, INCN, NNISUR, NBHM, IUSED, LUSED)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER(NPTS=1)
      DIMENSION I(LENIWK), R(LENRWK), NSPHCH(MAXMAT), KKSURF(MAXMAT), 
     5          KKBULK(MAXMAT), KK(MAXMAT), NPHASE(MAXMAT), 
     6          NNSURF(MAXMAT), NFSURF(MAXMAT), NLSURF(MAXMAT),  
     7          NNBULK(MAXMAT), NFBULK(MAXMAT), NLBULK(MAXMAT), 
     8          IISUR(MAXMAT), IMRSK(MAXMAT), IMISK(MAXMAT), 
     9          IMCSK(MAXMAT), IMPH(MAXMAT), INCN(MAXMAT),
     1          NNISUR(MAXMAT), NBHM(MAXMAT)
C
      LOGICAL L(LENLWK), KERR
C
      CHARACTER*16 C(LENCWK), VERSN*80 
      CHARACTER*18 DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO
C
      COMMON /PSPSPS/ NCKW, NSKW, NWT, NFU, NOX, NADD, NSTO, NXIN, NYIN,
     1                NHIN, NXST, NSCR, NABV, NBLW, NBUF, NTWP, NS, NSN,
     2                NF, NFN, NSSV, NA, NAA, NDR, NDC, NEQ, NACT, NSDN,
     3                NSD0, NSID, NDEN, NGR, NGRD,NDSK, NSDT, NWDT,
     4                NCDT, NDDT, NROP, NQEK,NQL, NTQL, NREX,NSDI,NSII,
     5                NAFR, NEIO, NBMX, NFLX, NFXC, NSTT, NRPS, NBPW,
     6                NHIM, NTSR, NP, NPA, NV, NQ, NARE, NTMP, NTIN,
     7                NFLR, NTAU, NHTR, NTAM, NPOW, NTEI, NTEL, NGFA,
     8                NSFA, NESH, NTIO,
     9                ICKW, ISKW, INCF, IIP, IKPR,IKSP, IEQ, 
     *                IKKP, IKFT, IKLT, IKST, INST,IKCC, IKCH, IKIO, 
     1                IKTF, IEIM, IEX, ITD, IMU, 
     2                LENR, LFLR, LHTR, LENE, LTIO, LELS, LRFS, LWAL, 
     3                LTDI, LACT, LIS, LIR,  LLES, LETC, ICC,  ICS,  
     4                IKS,  IMM, IPSY, IMSY
C
C                CHARACTER DATA FILE COMMON BLOCK
C
      COMMON /CDATAF/ VERSNN, DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO
C
      WRITE (LOUT, 15)
   15 FORMAT (
     1/' AURORA:  PERFECTLY STIRRED REACTOR CODE',
     2/'             VERSION 6.17, 98/4/15',
C*****precision > double
     3/'       DOUBLE PRECISION')
C*****END precision > double
C*****precision > single
C     3/'       SINGLE PRECISION')
C*****END precision > single
C
C          WRITE HEADER ON DATA FILE
C
      VERSN = '@(#) AURORA.F VERSION 6.17, DATE 98/4/15'
      REWIND (LSAVE)
      WRITE  (LSAVE) DHEAD
      WRITE  (LSAVE) VERSNN, VERSN
C
C          SET UP INTERNAL WORK POINTERS
C
      KERR = .FALSE.
      CALL PSPNT (LINKCK, LINKSK, LOUT, MAXPSR, MM, KK, II, IISUR, 
     1            IISUR0, NSPHCH, KKGAS, KKSURF, KKBULK, NPHASE, NNSURF,
     2            NFSURF, NLSURF, NNBULK, NFBULK, NLBULK, NIK, NATJ,
     3            NPTS, LENIWK, LENRWK, LENCWK, LENTWP, LENIEQ, LENEQ,
     4            LTOT, ITOT, NTOT, ICTOT, NT, NTE, NYS, NY, NSDEN,
     5            NSDENS, LENICK, LENCK, LENISK, LENSK, MAXMAT, IMISK,
     6            IMRSK, IMCSK, NPHMAX, IMPH, INCN, KKTOT, KKMAX, NMAT,
     7            NSPTOT, IISMAX, I, R, C, IUSED, NUSED, LUSED, KERR)
       IF (KERR) RETURN
C
C           CHECK FOR ENOUGH SPACE
C
      WRITE (LOUT, 9200) LENLWK, LTOT, LENIWK, ITOT, LENRWK, NTOT,
     1                   LENCWK, ICTOT
 9200 FORMAT (/,'               WORKING SPACE REQUIREMENTS',
     1        /,'                 PROVIDED        REQUIRED ',
     2        /,' LOGICAL  ', 2I15,
     3        /,' INTEGER  ', 2I15,
     4        /,' REAL     ', 2I15,
     5        /,' CHARACTER', 2I15,/)
C
      IF (LTOT.GT.LENLWK .OR. ITOT.GT.LENIWK .OR. NTOT.GT.LENRWK
     1                   .OR. ICTOT.GT.LENCWK) THEN
         WRITE (LOUT, *) '  FATAL ERROR, NOT ENOUGH WORK SPACE PROVIDED'
         RETURN
      ENDIF
C
      CALL PSRDRV (LIN, LOUT, LOUTSH, LREAD, LSAVE, LRECOV, MAXPSR, MM, 
     *             KK, II, 
     1             IISUR, IISUR0, NSPHCH, KKGAS, KKSURF, KKBULK, NPHASE, 
     2             NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK, NIK, 
     3             NATJ, NPTS, LENTWP, LENICK, LENCK, LENISK, LENSK, 
     4             LENIEQ, LENEQ, C(ICC), C(ICS), C(IKS),C(IMM),C(IPSY), 
     5             C(IMSY), I(ICKW), I(ISKW), I(INCN(1)),I(IMPH(1)), 
     6             I(INCF), I(IIP),  I(IKPR), I(IKSP), I(IEQ), I(IKKP), 
     7             I(IKFT), I(IKLT), I(IKST), I(INST), I(IKCC),I(IKCH), 
     8             I(IKIO), I(IKTF), I(IEIM), I(IEX),  I(ITD), I(IMU),
     9             L(LENR), L(LFLR), L(LHTR), L(LENE), L(LTIO), L(LELS),
     *             L(LRFS), L(LWAL), L(LTDI), L(LACT), L(LIS),  L(LIR),  
     1             L(LLES), L(LETC), 
     2             R(NCKW), R(NSKW), R(NWT),  R(NFU),  R(NOX),  R(NADD), 
     3             R(NSTO), R(NXIN), R(NYIN), R(NHIN), R(NXST), R(NSCR), 
     4             R(NABV), R(NBLW), R(NS),   R(NSN),  R(NF),   R(NFN),  
     5             R(NSSV), R(NA),   R(NAA),  R(NDR),  R(NDC),  R(NBUF), 
     6             R(NTWP), R(NEQ),  R(NACT), R(NSDN), R(NSD0), R(NSID), 
     7             R(NDEN), R(NGR),  R(NGRD), R(NDSK), R(NSDT), R(NWDT), 
     8             R(NCDT), R(NDDT), R(NROP), R(NQEK), R(NQL),  R(NTQL), 
     9             R(NREX), R(NSDI), R(NSII), R(NAFR), R(NEIO), R(NBMX),
     *             R(NFLX), R(NFXC), R(NSTT), R(NRPS), R(NBPW), R(NHIM),
     1             R(NTSR), R(NPOW), R(NP),   R(NPA),  R(NV),   R(NQ), 
     2             R(NHTR), R(NTAM), R(NARE), R(NTMP), R(NTIO), R(NTIN),
     3             R(NFLR), R(NTAU), R(NTEI), R(NTEL), R(NGFA), R(NSFA),
     4             R(NESH),
     5             NT, NTE, NYS, NY, 
     6             NSDEN, NSDENS, IMISK, IMRSK, IMCSK, NPHMAX,
     7             KKTOT, NMAT, NSPTOT, KKMAX, IISMAX, 
     8             NNISUR, NBHM)
      RETURN
      END
      BLOCK DATA
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      CHARACTER*18 DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO
      COMMON /CDATAF/ VERSNN, DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO
      DATA DHEAD/'PSR CODE DATA FILE'/
      DATA ISOLUT/'SOLUTION          '/
      DATA ICHEM/'CHEMISTRY         '/
      DATA IROPRO /'RATE OF PRODUCTION'/
      DATA ISENSI/'SENSITIVITY       '/
      DATA VERSNN/6.17/
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE PSPNT (LINKCK, LINKSK, LOUT, MAXPSR, MM, KK, II, IISUR,
     1                  IISUR0, NSPHCH, KKGAS, KKSURF, KKBULK, NPHASE,
     2                  NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     3                  NIK, NATJ, NPTS, LENIWK, LENRWK, LENCWK, LENTWP,
     4                  LENIEQ, LENEQ, LTOT, ITOT, NTOT, ICTOT, NT, NTE,
     5                  NYS, NY, NSDEN, NSDENS, LENICK, LENRCK, LENISK,
     6                  LENRSK, MAXMAT, IMISK, IMRSK, IMCSK, NPHMAX,
     7                  IMPH, INCN, KKTOT, KKMAX, NMAT, NSPTOT, IISMAX,
     8                  I, R, C, IUSED, NUSED, LUSED, KERR)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C INPUT:
C   LINKCK  - UNIT NUMBER FOR GAS-PHASE LINK FILE
C   LINKSK  - UNIT NUMBER FOR SURFACE LINK FILE
C   LOUT    - UNIT NUMBER FOR STANDARD OUTPUT
C   MAXPSR  - MAXIMUM NUMBER OF PSR'S IN SERIES
C   NPTS    - NUMBER OF TWOPNT POINTS = 1
C   LENIWK  - USER-SUPPLIED LENGTH OF THE I WORK SPACE
C   LENRWK  - USER-SUPPLIED LENGTH OF THE R WORK SPACE
C   LENCWK  - USER-SUPPLIED LENGTH OF THE C WORK SPACE
C   MAXMAT  - THE MAXIMUM NUMBER OF MATERIALS ALLOWED
C   I       - THE INTEGER WORK SPACE
C   R       - THE REAL WORK SPACE
C   C       - THE CHARACTER WORK SPACE
C   LUSED   - THE AMOUNT OF LOGICAL WORK SPACE ALREADY USED IN PSRUN
C   NUSED   - THE AMOUNT OF REAL WORK SPACE ALREADY USED IN PSRUN
C   IUSED   - THE AMOUNT OF INTEGER WORK SPACE ALREADY USED IN PSRUN
C
C OUTPUT:
C   MM      - TOTAL NUMBER OF DIFFERENT ELEMENTS USED IN MECHANISM
C   KK      - TOTAL NUMBER OF SPECIES FOR EACH GAS/MATERIAL MECHANISM
C             DIMENSION AT LEAST KK(MAXMAT)
C   II      - TOTAL NUMBER OF GAS PHASE REACTIONS
C   IISUR   - TOTAL NUMBER OF SURFACE REACTIONS ON EACH MATERIAL
C             DIMENSION AT LEAST IISUR(MAXMAT)
C   IISUR0  - MINIMUM IISUR DIMENSION FOR MULTIDIMENSIONED VARIABLES
C   NSPHCH  - NUMBER OF SURF PHASES WHERE THE SITE DENSITY CAN CHANGE
C             DIMENSION AT LEAST NSPHCH(MAXMAT)
C   KKGAS   - TOTAL NUMBER OF GAS PHASE SPECIES IN THE MECHANISM
C   KKSURF  - TOTAL NUMBER OF SURFACE SPECIES IN THE MECHANISM.
C             DIMENSION AT LEAST KKSUR(MAXMAT)
C   KKBULK  - TOTAL NUMBER OF BULK PHASE SPECIES
C             DIMENSION AT LEAST KKBULK(MAXMAT)
C   NPHASE  - TOTAL NUMBER OF PHASES FOR EACH MATERIAL
C             DIMENSION AT LEAST NPHASE(MAXMAT)
C   NNSURF  - NUMBER OF SURFCE PHASES FOR EACH MATERIAL
C             DIMENSION AT LEAST NNSURF(MAXMAT)
C   NFSURF  - INDEX OF THE FIRST SURFACE PHASE OF EACH MATERIAL
C             DIMENSION AT LEAST NFSURF(MAXMAT)
C   NLSURF  - INDEX OF THE LAST SURFACE PHASE OF EACH MATERIAL
C             DIMENSION AT LEAST NLSURF(MAXMAT)
C   NNBULK  - NUMBER OF BULK PHASES FOR EACH MATERIAL
C             DIMENSION AT LEAST NNBULK(MAXMAT)
C   NFBULK  - INDEX OF FIRST BULK PHASE FOR EACH MATERIAL
C             DIMENSION AT LEAST NFBULK(MAXMAT)
C   NLBULK  - INDEX OF LAST BULK PHASE FOR EACH MATERIAL
C             DIMENSION AT LEAST NLBULK(MAXMAT)
C   NIK     - THE VALUE OF THE MAX NUMBER OF REACTIONS 
C             TIMES THE MAX NUMBER OF SPECIES
C   NATJ    - TOTAL NUMBER OF DEPENDENT VARIABLES = KK + 2 + NSPTOT
C   LENTWP  - LENGTH OF THE TWOPNT WORKSPACE
C   LENIEQ  - LENGTH OF THE INTEGER EQUILIBRIUM-CODE WORKSPACE.
C   LENEQ   - LENGTH OF THE FLOATING-PT EQ-CODE WORKSPACE.
C   LTOT    - TOTAL LOGICAL WORKSPACE ALLOCATED FOR THIS PROBLEM
C   ITOT    - TOTAL INTEGER WORKSPACE ALLOCATED FOR THIS PROBLEM
C   NTOT    - TOTAL REAL WORKSPACE ALLOCATED FOR THIS PROBLEM
C   ICTOT   - TOTAL CHARACTER WORKSPACE ALLOCATED FOR THIS PROBLEM
C   NT      - SOLUTION POINTER TO THE TEMPERATURE VARIABLE
C   NTE     - SOLUTION POINTER TO THE ELECTRON TEMPERATURE VARIABLE
C   NYS     - SOLUTION POINTER TO THE FIRST SPECIES VARIABLE - 1
C   NY      - SOLUTION POINTER TO THE FIRST SPECIES VARIABLE 
C   NSDEN   - SOLUTION POINTER TO THE FIRST VARYING PHASE DENSITY
C   NSDENS  - SOLUTION POINTER TO THE FIRST VARYING PHASE DENSITY - 1
C   LENICK  - LENGTH OF THE CHEMKIN INTEGER WORKSPACE
C   LENRCK  - LENGTH OF THE CHEMKIN REAL WORKSPACE
C   LENISK  - LENGTH OF THE SURFACE CHEMKIN INTEGER WORKSPACE
C   LENRSK  - LENGTH OF THE SURFACE CHEMKIN REAL WORKSPACE
C   IMISK   - POINTER TO THE NTH MATERIAL IN THE ISKWRK ARRAY
C   IMRSK   - POINTER TO THE NTH MATERIAL IN THE RSKWRK ARRAY
C   IMCSK   - POINTER TO THE NTH MATERIAL IN THE CSKWRK ARRAY
C   NPHMAX  - THE MAXIMUM NUMBER OF PHASES ON ANY MATERIAL
C   IMPH    - INDICES OF THE SURFACE PHASES WHICH CAN CHANGE
C             DIMENSION IMPH(MAXMAT).
C   INCN    - INDICES OF THE SURFACE PHASES THAT DO NOT CONSERVE SITES
C             DIMENSION INCN(MAXMAT).
C   KKTOT   - TOTAL NUMBER OF SPECIES OF ALL KINDS
C   KKMAX   - MAXIMUM TOTAL NUMBER OF SPECIES ON ANY MATERIAL.
C   NMAT    - THE NUMBER OF MATERIALS FOUND ON THE SURF LINKING FILE
C   NSPTOT  - THE TOTAL NUMBER OF PHASES W/CHANGING SITE DENSITY
C   IISMAX  - THE MAXIMUM NUMBER OF SURFACE REACTIONS ON ANY MATERIAL
C   KERR    - LOGICAL ERROR FLAG FOR ALLOCATING WORK SPACE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION R(LENRWK), I(LENIWK), KKSURF(MAXMAT), KKBULK(MAXMAT), 
     1          KK(MAXMAT), NPHASE(MAXMAT), NNSURF(MAXMAT), 
     2          NFSURF(MAXMAT), NLSURF(MAXMAT), NNBULK(MAXMAT), 
     3          NFBULK(MAXMAT), NLBULK(MAXMAT), IISUR(MAXMAT),
     4          IMRSK(MAXMAT), IMISK(MAXMAT), IMCSK(MAXMAT),
     5          IMPH(MAXMAT), INCN(MAXMAT), NSPHCH(MAXMAT)
      CHARACTER*16 C(LENCWK)
      LOGICAL KERR
      COMMON /PSPSPS/ NCKW, NSKW, NWT, NFU, NOX, NADD, NSTO, NXIN, NYIN,
     1                NHIN, NXST, NSCR, NABV, NBLW, NBUF, NTWP, NS, NSN,
     2                NF, NFN, NSSV, NA, NAA, NDR, NDC, NEQ, NACT, NSDN,
     3                NSD0, NSID, NDEN, NGR, NGRD,NDSK, NSDT, NWDT,
     4                NCDT, NDDT, NROP, NQEK,NQL, NTQL, NREX,NSDI,NSII,
     5                NAFR, NEIO, NBMX, NFLX, NFXC, NSTT, NRPS, NBPW,
     6                NHIM, NTSR, NP, NPA, NV, NQ, NARE, NTMP, NTIN,
     7                NFLR, NTAU, NHTR, NTAM, NPOW, NTEI, NTEL, NGFA,
     8                NSFA, NESH, NTIO,
     9                ICKW, ISKW, INCF, IIP, IKPR,IKSP, IEQ, 
     *                IKKP, IKFT, IKLT, IKST, INST,IKCC, IKCH, IKIO, 
     1                IKTF, IEIM, IEX, ITD, IMU, 
     2                LENR, LFLR, LHTR, LENE, LTIO, LELS, LRFS, LWAL, 
     3                LTDI, LACT, LIS, LIR,  LLES, LETC, ICC,  ICS,  
     4                IKS,  IMM, IPSY, IMSY
C
      KERR = .FALSE.
C           START TO APPORTION THE WORK SPACE (FOR CHEMKIN AND SURFACE)
C
      CALL CKLEN (LINKCK, LOUT, LENICK, LENRCK, LENCCK, IFLAG)
      IF (IFLAG .NE. 0) THEN
         WRITE (LOUT,9000) IFLAG
C         STOP
         KERR = .TRUE.
         RETURN
      ENDIF
C
C  SET UP COUNTER AND LOOP OVER MULTIPLE MATERIALS IN SK LINK FILE
C
      NMAT = 0
      LENRSK = 0
      LENISK = 0
      LENCSK = 0
      LTOT = 1
C
10    CONTINUE
C
      NMAT = NMAT + 1
      LSTRSK = LENRSK
      LSTISK = LENISK
      LSTCSK = LENCSK
      CALL SKLEN (LINKSK, LOUT, LENISK, LENRSK, LENCSK, IFLAG)
      IF (IFLAG .NE. 0) THEN
         WRITE (LOUT,9005) IFLAG
C         STOP
         KERR = .TRUE.
         RETURN
      ENDIF
C
      IF (NMAT .EQ. 1) THEN
         IMRSK(1) = 1
         IMISK(1) = 1
         IMCSK(1) = 1
C  REAL WORK SPACE ARRAYS
         NCKW = NUSED + 1
         NSKW = NCKW + LENRCK
         NTOT = NSKW + LENRSK
C  INTEGER WORK SPACE ARRAYS
         ICKW = IUSED + 1
         ISKW = ICKW + LENICK
         ITOT = ISKW + LENISK
C  CHARACTER WORK SPACE ARRAYS
         ICC  = 1
         ICS  = ICC + LENCCK
         ICTOT = ICS + LENCSK
C
C              INITIALIZE CHEMKIN AND SURFACE CHEMKIN
C
         IF (ITOT.LT.LENIWK .AND. NTOT.LT.LENRWK .AND.
     1      ICTOT.LT.LENCWK) THEN
            CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT, 
     1                   I(ICKW), R(NCKW), C(ICC), IFLAG)
            CALL CKINDX (I(ICKW), R(NCKW), MM, KKGAS, II, NFIT)
            CALL SKINIT (LENISK, LENRSK, LENCSK, LINKSK, LOUT, 
     1                   I(ISKW), R(NSKW), C(ICS), IFLAG)
            CALL SKINDX (I(ISKW),MM, KKGAS, KKSURF(1), KKBULK(1),  
     1                   KK(1), NPHASE(1), NNSURF(1), NFSURF(1), 
     2                   NLSURF(1), NNBULK(1), NFBULK(1), NLBULK(1), 
     3                   IISUR(1))
         ELSE
            WRITE (LOUT, *)'NOT ENOUGH WORK SPACE ALLOCATED TO START'
            RETURN
         ENDIF
         KKSMAX = KKSURF(1)
         KKBMAX = KKBULK(1)
         KKMAX = KK(1)
         KKTOT = KK(1)
         NPHMAX = NPHASE(1)
         NNSMAX = NNSURF(1)
         NNBMAX = NNBULK(1)
         IISMAX = IISUR(1)
      ELSE
         IMRSK(NMAT) = IMRSK(NMAT-1) + LSTRSK
         IMISK(NMAT) = IMISK(NMAT-1) + LSTISK
         IMCSK(NMAT) = IMCSK(NMAT-1) + LSTCSK
         NTOT = NTOT + LENRSK
         ITOT = ITOT + LENISK
         ICTOT = ICTOT + LENCSK
C
C              INITIALIZE CHEMKIN AND SURFACE CHEMKIN
C
         IF (ITOT.LT.LENIWK .AND. NTOT.LT.LENRWK .AND.
     1      ICTOT.LT.LENCWK) THEN
            CALL SKINIT (LENISK, LENRSK, LENCSK, LINKSK, LOUT, 
     1                   I(ISKW+IMISK(NMAT)-1), R(NSKW+IMRSK(NMAT)-1), 
     2                   C(ICS+IMCSK(NMAT)-1), IFLAG)
            CALL SKINDX (I(ISKW+IMISK(NMAT)-1),MM, KKGAS, KKSURF(NMAT), 
     1                   KKBULK(NMAT), KK(NMAT), NPHASE(NMAT), 
     2                   NNSURF(NMAT), NFSURF(NMAT), NLSURF(NMAT), 
     3                   NNBULK(NMAT), NFBULK(NMAT), NLBULK(NMAT),
     4                   IISUR(NMAT))
         ELSE
            WRITE (LOUT, *)
     1          'NOT ENOUGH WORK SPACE FOR SURFACE MATERIAL ', NMAT
C            STOP
            KERR = .TRUE.
            RETURN
         ENDIF
         KKSMAX = MAX(KKSMAX, KKSURF(NMAT))
         KKBMAX = MAX(KKBMAX, KKBULK(NMAT))
         KKMAX = MAX(KKMAX, KK(NMAT))
         KKTOT = KKTOT + KKSURF(NMAT) + KKBULK(NMAT)
         NPHMAX = MAX(NPHMAX, NPHASE(NMAT))
         NNSMAX = MAX(NNSMAX, NNSURF(NMAT))
         NNBMAX = MAX(NNBMAX, NNBULK(NMAT))
         IISMAX = MAX(IISMAX, IISUR(NMAT))
      ENDIF
C
C   CHECK END OF WORKSPACE FOR FLAG INDICATING ANOTHER MATERIAL
C
      MORE = I(ISKW + IMISK(NMAT) - 1 + LENISK - 1)
      IF (MORE .EQ. 1) THEN
         IF (NMAT+1 .GT. MAXMAT) THEN
            WRITE (LOUT, *)
     1      'STOP...number of materials NMAT exceeds MAXMAT in ',
     2      'calling program...'
            KERR = .TRUE.
            RETURN
         ENDIF
         GO TO 10
      ENDIF
C
C    RETURN THE LENGTH OF THE APPENDED SURFACE LINKING FILES FOR
C    USE IN DIMENSIONING ISKWRK
C
      LENISK = LENISK + IMISK(NMAT)
      LENRSK = LENRSK + IMRSK(NMAT)
      LENCSK = LENCSK + IMCSK(NMAT)
C
C    TAKE CARE OF THE IISUR = 0 CASE
C
      IISUR0 = MAX(1,IISMAX)

C    DETERMINE HOW MANY SURFACE PHASE DENSITIES ARE SUBJECT TO CHANGE
C
      NSPTOT = 0
      DO 20 IMAT = 1, NMAT
         INCN(IMAT) = ITOT
         IMPH(IMAT) = INCN(IMAT) + NPHMAX
         ITOT = IMPH(IMAT) + NPHMAX
         DO 14 NPH = 1, NPHMAX
	    I(INCN(IMAT)-1+NPH) = 0
            I(IMPH(IMAT)-1+NPH) = 0
 14      CONTINUE
         IF (NNSURF(IMAT) .GT. 0) THEN
            CALL SKNCON (I(ISKW+IMISK(IMAT)-1), R(NSKW+IMRSK(IMAT)-1),
     1                   I(INCN(IMAT)))
            NSPHCH(IMAT) = 0
            DO 15 IPHASE = NFSURF(IMAT), NLSURF(IMAT)
               IF (I(INCN(IMAT)-1+IPHASE) .NE. 0 ) THEN
                  NSPHCH(IMAT) = NSPHCH(IMAT) + 1
                  I(IMPH(IMAT)-1+NSPHCH(IMAT)) = IPHASE
               ENDIF
15          CONTINUE
         ELSE
            NSPHCH(IMAT) = 0
         ENDIF
         NSPTOT = NSPTOT + NSPHCH(IMAT)
20    CONTINUE

      NCS = 1
      NCON = 0
      KKBEQ = 0
      NNPEQ = 1
      NNBEQ = 0
C
C.....5/19/96 upgrade to new equilibrium code, per F. Rupley
C
      CALL EQLEN (MM, KKGAS, KKBEQ, KKGAS, NNPEQ, NNBEQ, LENIEQ, 
     1            LENEQ, LENCEQ)
C.....
C
      NATJ = KKTOT + 2 + NSPTOT
C
C          SET THE POINTERS INTO THE SOLUTION VECTOR
C
      NT  = 1
      NTE = 2
      NYS = NTE
      NY  = NTE + 1
      NSDEN = NY + KKTOT
      NSDENS = NSDEN - 1
C
      LENTWP = 7*NATJ + (7*NATJ + 2)
C
C          APPORTION THE REST OF THE FLOATING POINT SPACE
C
      NP   = NTOT
      NPA  = NP   + MAXPSR
      NV   = NPA  + MAXPSR
      NQ   = NV   + MAXPSR
      NARE = NQ   + MAXPSR
      NTMP = NARE + MAXPSR
      NTIN = NTMP + MAXPSR
      NFLR = NTIN + MAXPSR
      NTAU = NFLR + MAXPSR
      NHTR = NTAU + MAXPSR
      NTAM = NHTR + MAXPSR
      NPOW = NTAM + MAXPSR
      NTEI = NPOW + MAXPSR
      NTEL = NTEI + MAXPSR
      NGFA = NTEL + MAXPSR
      NSFA = NGFA + MAXPSR
      NESH = NSFA + MAXPSR
      NTIO = NESH + NMAT * MAXPSR
      NWT  = NTIO + MAXPSR
      NFU  = NWT  + KKMAX*NMAT
      NOX  = NFU  + KKMAX
      NADD = NOX  + KKMAX
      NSTO = NADD + KKMAX
      NXIN = NSTO + KKMAX
      NYIN = NXIN + KKTOT * MAXPSR
      NHIN = NYIN + KKMAX
      NXST = NHIN + KKMAX
      NSCR = NXST + KKTOT
                    NIK = MAX( IISUR0, MAX(NATJ,II), IISUR0*NMAT)
      NABV = NSCR + NIK * MAX(9,NMAT)
      NBLW = NABV + NATJ
      NBUF = NBLW + NATJ
      NTWP = NBUF + NATJ
      NS   = NTWP + LENTWP
      NSN  = NS   + NATJ
      NF   = NSN  + NATJ
      NFN  = NF   + NATJ
      NSSV = NFN  + NATJ
      NA   = NSSV + NATJ * MAXPSR
      NAA  = NA   + NATJ * NATJ
      NDR  = NAA  + NATJ * NATJ
      NDC  = NDR  + NATJ
      NEQ  = NDC  + NATJ
      NACT = NEQ  + LENEQ
      NSDN = NACT + KKMAX * NMAT
      NSD0 = NSDN + NPHMAX * NMAT
      NSID = NSD0 + NPHMAX * NMAT
      NDEN = NSID + NPHMAX * NMAT
      NGR  = NDEN + KKMAX * NMAT
      NGRD = NGR  + NPHMAX * NMAT
      NDSK = NGRD + NPHMAX * NMAT
      NSDT = NDSK + KKTOT*NATJ
      NWDT = NSDT + KKTOT
      NCDT = NWDT + KKTOT
      NDDT = NCDT + KKTOT
      NROP = NDDT + KKTOT
      NQEK = NROP + MAX(II,IISUR0)
      NQL  = NQEK + KKGAS
      NTQL = NQL  + NPTS
      NREX = NTQL + NPTS
      NSDI = NREX + II
      NSII = NSDI + KKMAX * IISUR0 * NMAT
      NAFR = NSII + NPHMAX * NMAT
      NEIO = NAFR + NMAT * MAXPSR
      NBMX = NEIO + NMAT * MAXPSR
      NFLX = NBMX + MAXPSR
      NFXC = NFLX + KKMAX
      NSTT = NFXC + KKMAX
      NRPS = NSTT + KKTOT
      NBPW = NRPS + IISUR0 * NMAT
      NHIM = NBPW + NMAT * MAXPSR
      NTSR = NHIM + NMAT * KKMAX
      NTOT = NTSR + NMAT * MAXPSR
C
C           APPORTION THE REST OF THE INTEGER SPACE
C
      INCF = ITOT
      IIP  = INCF + MM * KKMAX * NMAT
      IKPR = IIP  + NATJ
      IKSP = IKPR + KKMAX
      IEQ  = IKSP + KKMAX
      IKKP = IEQ  + LENIEQ
      IKFT = IKKP + NPHMAX * NMAT
      IKLT = IKFT + NPHMAX * NMAT
      IKST = IKLT  + NPHMAX * NMAT
      INST = IKST  + IISUR0 * KKMAX * NMAT
      IKCC = INST  + IISUR0 * NPHMAX * NMAT
      IKCH = IKCC + KKMAX * NMAT
      IKIO = IKCH + KKMAX
      IKTF = IKIO + KKMAX
      IEIM = IKTF + KKTOT
      IEX  = IEIM + II
      ITD  = IEX  + II
      IMU  = ITD  + II
      ITOT = IMU  + MM
C
C           APPORTION THE LOGICAL SPACE
C
      LENR = LUSED + 1
      LFLR = LENR + MAXPSR
      LHTR = LFLR + MAXPSR
      LENE = LHTR + MAXPSR
      LTIO = LENE + MAXPSR
      LELS = LTIO + MAXPSR
      LRFS = LELS + NMAT * MAXPSR
      LWAL = LRFS + NMAT * MAXPSR
      LTDI = LWAL + NMAT * MAXPSR
      LACT = LTDI + NMAT * MAXPSR
      LIS  = LACT + NATJ
      LIR  = LIS  + NATJ
      LLES = LIR  + NATJ
      LETC = LLES + NPHMAX * NMAT
      LTOT = LETC  + NPHMAX * MAXPSR * NMAT
C
C           APPORTION THE REST OF THE CHARACTER SPACE
C
      IKS  = ICTOT
      IMM  = IKS + KKMAX*NMAT
      IPSY = IMM + MM
      IMSY = IPSY + NPHMAX*NMAT
      ICTOT = IMSY + NMAT
C
      RETURN
9000  FORMAT (1X,'ERROR READING CHEMKIN LINKING FILE, IFLAG = ',I3)
9005  FORMAT (1X,'ERROR READING SURFACE CHEMKIN LINKING FILE, 
     1        IFLAG = ',I3)
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE PSRDRV (LIN, LOUT, LOUTSH, LREAD, LSAVE, LRECOV,
     1                   MAXPSR, MM,
     1                   KK, II, IISUR, IISUR0, NSPHCH, KKGAS, KKSURF, 
     2                   KKBULK, NPHASE, NNSURF, NFSURF, NLSURF, NNBULK,
     3                   NFBULK, NLBULK, NIK, NATJ, NPTS, LENTWP, 
     4                   LENICK, LENCK, LENISK, LENSK, LENIEQ, LENEQ, 
     5                   CCWORK, CSWORK, KSYM, ATOM, PSYM, MSYM, ICKWRK, 
     6                   ISKWRK, NCON, MAPPH, NCF, IPVT, KPROD, KSP, 
     7                   IEQWRK, KKPHAS, KFIRST, KLAST, KSTOIC, NSTOIC,
     8                   KOCC, KCHG, KION, KTFL, IEIMP, IEXC, ITDEP, 
     9                   MUSD, LENRGY, LFLRT, LHTRN, LENRGE, LTION, 
     *                   LELSH, LRFSH, LWALHB, LTDIFF, 
     1                   ACTIVE, ISEN, IROP, LESTIM, ETCH, RCKWRK, 
     2                   RSKWRK, WT, FUEL, OXID, ADD, STOICH, XIN, YIN, 
     3                   HIN, X, SCRTCH, ABOVE, BELOW, S, SN, F, FN, 
     4                   SSAVE, A, AA, DR, DC, BUFFER, TWPWK, EQWRK, 
     5                   ACT, SDEN, SDEN0, SITDOT, DEN, GRATE, GRATED, 
     6                   DSKDPL, SDOT, WDOT, CDOT, DDOT, ROP, QEK, QLSE, 
     7                   TQLSE,REXC, SDOTI, SITDTI, AFRAC, EIONSH, 
     8                   BHMXI, FLXION, FLXCOR, SDOTT, ROPS, BIASPW, 
     9                   HIM, TSURF, POWR, P, PA, V, Q, HTRN, TAMBNT, 
     *                   AREA,T, TIONP, TIN, FLRT, TAU, TEIN, TE, GFAC, 
     1                   SFAC,  ESHTH, 
     2                   NT, NTE, NYS, NY, NSDEN, NSDENS, 
     3                   IMISK, IMRSK, IMCSK, NPHMAX, KKTOT, NMAT, 
     4                   NSPTOT, KKMAX, IISMAX, NNISUR, NBHM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DOUBLE PRECISION MASLOS, INLET
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      REAL MASLOS, INLET
C*****END precision > single
C
      DIMENSION KKSURF(NMAT), KKBULK(NMAT), NPHASE(NMAT), NNSURF(NMAT),
     1          NFSURF(NMAT), NLSURF(NMAT), NNBULK(NMAT), NFBULK(NMAT),
     2          NLBULK(NMAT), NSPHCH(NMAT), IMISK(NMAT),
     3          IMRSK(NMAT), IMCSK(NMAT), IISUR(NMAT), KK(NMAT)
      DIMENSION ICKWRK(LENICK), ISKWRK(LENISK), NCON(NPHMAX,NMAT), 
     1          KPROD(KKMAX), MAPPH(NPHMAX,NMAT), NCF(MM,KKMAX,NMAT),
     1          IPVT(NATJ), KSP(KKMAX), IEQWRK(LENIEQ),
     2          KKPHAS(NPHMAX,NMAT), KFIRST(NPHMAX,NMAT),
     3          KLAST(NPHMAX,NMAT), KSTOIC(IISUR0,KKMAX,*),
     4          NSTOIC(IISUR0,NPHMAX,*), KOCC(KKMAX,NMAT), 
     4          RCKWRK(LENCK), RSKWRK(LENSK), WT(KKMAX,NMAT),
     6          EQWRK(LENEQ), FUEL(KKMAX), OXID(KKMAX), ADD(KKMAX), 
     7          STOICH(KKMAX), XIN(KKTOT,MAXPSR), YIN(KKMAX),
     7          HIN(KKMAX), X(KKTOT), SCRTCH(NIK,9), S(NATJ), SN(NATJ), 
     8          F(NATJ), FN(NATJ), SSAVE(NATJ,MAXPSR), A(NATJ,NATJ),
     8          AA(NATJ,NATJ), DR(NATJ), DC(NATJ), ABOVE(NATJ),
     9          BELOW(NATJ), BUFFER(NATJ), TWPWK(*), ACT(KKMAX,NMAT),
     *          SDEN(NPHMAX,NMAT), SDEN0(NPHMAX,NMAT), 
     1          SITDOT(NPHMAX,NMAT), DEN(KKMAX,NMAT),
     1          GRATE(NPHMAX,NMAT), GRATED(NPHMAX,NMAT),
     2          DSKDPL(KKTOT,NATJ), SDOT(KKTOT), WDOT(KKTOT), 
     3          CDOT(KKTOT), DDOT(KKTOT), P(MAXPSR), PA(MAXPSR),
     3          V(MAXPSR), Q(MAXPSR), AREA(MAXPSR), T(MAXPSR),
     4          TSURF(NMAT,MAXPSR), TIN(MAXPSR), FLRT(MAXPSR), TEMP(3),
     5          TAU(MAXPSR), HTRN(MAXPSR), TAMBNT(MAXPSR), LEVEL(2),
     6          KCHG(KKMAX), KTFL(KKTOT), KION(*), QEK(KKGAS),  
     7          POWR(MAXPSR), TEIN(MAXPSR), TE(MAXPSR), ROP(*),
     7          IEIMP(*), QLSE(NPTS), ITDEP(*), TQLSE(NPTS), IEXC(*),
     8          REXC(*), SDOTI(KKMAX,IISUR0,NMAT), SITDTI(NPHMAX,NMAT),
     9          AFRAC(NMAT,MAXPSR), GFAC(MAXPSR), SFAC(MAXPSR),
     *          EIONSH(NMAT,MAXPSR), BHMXI(MAXPSR), FLXION(*),FLXCOR(*), 
     1          SDOTT(KKTOT), ESHTH(NMAT,MAXPSR), ROPS(IISUR0,NMAT),
     2          NNISUR(*), NBHM(NMAT), MUSD(MM), BIASPW(NMAT,MAXPSR), 
     3          TIONP(MAXPSR), HIM(KKMAX,NMAT)
C
      PARAMETER(BOLTZ=1.380658D-16, AVOG=6.0221367D23, 
     &          ECHG=1.60217733D-12, EV2K=ECHG/BOLTZ)
      CHARACTER*16  KSYM(KKMAX,NMAT), ATOM(MM), CCWORK(*), CSWORK(*), 
     1              PSYM(NPHMAX,NMAT), MSYM(NMAT)
      CHARACTER*18 DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO, ICHR
      CHARACTER    ADDON*14, STRING*120
      CHARACTER*10 CHNUM, TRNFIL
C
      LOGICAL ACTIVE(NATJ), ISEN(KKTOT), IROP(KKTOT), 
     1        LESTIM(NPHMAX,NMAT), ETCH(NPHMAX,NMAT,MAXPSR),
     2        LENRGY(MAXPSR), LTDIFF(NMAT,MAXPSR),
     2        LFLRT(MAXPSR), LHTRN(MAXPSR), LEQUIV, LRSTRT, LCNTUE,
     3        RSTCNT, LSEN, LSENT, LSENG, LROP, LPRTIC, KERR, IERR,
     4        LENRGE(MAXPSR),LNOFT, LCONFN, LDUM, LTION(MAXPSR), 
     5        LELSH(NMAT,MAXPSR), LRFSH(NMAT,MAXPSR), LRFCUR, 
     6        LWALHB(NMAT,MAXPSR), LSTEDY, LTERRN, LSCCM
C
C                CHARACTER DATA FILE COMMON BLOCK
C
      COMMON /CDATAF/ VERSNN, DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO
      DATA LCNTUE/.FALSE./, KERR/.FALSE./
      DATA CHNUM/'0123456789'/
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C         WRITE OUT CHEMISTRY TO DATA FILE
C
      WRITE (LSAVE) ICHEM
      CALL CKSAVE (LOUT, LSAVE, ICKWRK, RCKWRK, CCWORK)
      DO 100 IM = 1, NMAT
         CALL SKSAVE (LOUT, LSAVE, ISKWRK(IMISK(IM)), 
     1                RSKWRK(IMRSK(IM)), CSWORK(IMCSK(IM)))
100   CONTINUE
C
C*****precision > double
      ABSOL = SQRT(D1MACH(4))
      RELAT = SQRT(D1MACH(4))
C*****END precision > double
C*****precision > single
C      ABSOL = SQRT(R1MACH(4))
C      RELAT = SQRT(R1MACH(4))
C*****END precision > single
C
      RSTCNT = .FALSE.
C
C          READ IN VARIOUS CHEMKIN AND SURFACE CHEMKIN ARRAYS
C
      CALL CKRP  (ICKWRK, RCKWRK, RU, RUC, PATM)
      CALL SKSYME (ISKWRK(IMISK(1)), CSWORK(IMCSK(1)), LOUT, ATOM, IERR)
         KERR = KERR.OR.IERR
C
C          LOOP OVER MULTIPLE MATERIALS FOR SURFACE INFORMATION
C
      DO 500 IM = 1, NMAT
         CALL SKSYMM (ISKWRK(IMISK(IM)), CSWORK(IMCSK(IM)), LOUT, 
     1                MSYM(IM), IERR)
         CALL SKSYMS (ISKWRK(IMISK(IM)), CSWORK(IMCSK(IM)), LOUT, 
     1                KSYM(1,IM), IERR)
            KERR = KERR.OR.IERR
         CALL SKSYMP (ISKWRK(IMISK(IM)), CSWORK(IMCSK(IM)), LOUT, 
     1                PSYM(1,IM), IERR)
         KERR = KERR.OR.IERR
         CALL SKWT   (ISKWRK(IMISK(IM)), RSKWRK(IMRSK(IM)), WT(1,IM))
         CALL SKNCF  (MM, ISKWRK(IMISK(IM)), NCF(1,1,IM))
         CALL SKSDEN (ISKWRK(IMISK(IM)), RSKWRK(IMRSK(IM)), SDEN0(1,IM))
         DO 400 IPHASE = 1, NPHASE(IM)
            SDEN(IPHASE,IM) = SDEN0(IPHASE,IM)
400      CONTINUE
         CALL SKPKK (ISKWRK(IMISK(IM)), KKPHAS(1,IM), KFIRST(1,IM), 
     1               KLAST(1,IM))
         CALL SKNU  (IISUR(IM), ISKWRK(IMISK(IM)), RSKWRK(IMRSK(IM)), 
     1               KSTOIC(1,1,IM), NSTOIC(1,1,IM))
         CALL SKCOV (ISKWRK(IMISK(IM)), KOCC(1,IM))
500   CONTINUE
C
C         FIND K FOR ELECTRON AND IONS AND TOTAL NUMBER OF IONS
C
      CALL PKINDX (ICKWRK, KEL, KKION)
      CALL CKION  (ICKWRK, KION)
      CALL CKCHRG (ICKWRK, RCKWRK, KCHG)
C
C         FIND ELECTRON-IMPACT REACTIONS, IF ANY
C         FIND EXCITATION REACTIONS, IF ANY
C
      DO 600 I = 1, II
         IEIMP(I) = 0   
         ITDEP(I) = 0
         IEXC(I)  = 0
600   CONTINUE
      IF (KEL.NE.0) THEN
         CALL CKIEIM (ICKWRK, RCKWRK, IEIMP)
         CALL CKITDE (ICKWRK, RCKWRK, ITDEP)
         CALL CKIEXC (ICKWRK, RCKWRK, IEXC, REXC)
      ENDIF
C
C         SET TEMPERATURE FLAG ARRAY FOR ELECTRONS AND NEUTRALS=IONS
C
      DO 800 K = 1, KKTOT
         KTFL(K) = 1
 800  CONTINUE
      IF (KEL.NE.0) KTFL(KEL) = 2
      DO 810 KI = 1, KKION
         KTFL(KION(KI)) = 3
 810  CONTINUE
C
      CALL CKKTFL (ICKWRK, KTFL)
      DO 820 IM = 1, NMAT
C
C  Initialize species temperature flag in sklib
c
         CALL SKKTFL (ISKWRK(IMISK(IM)), KTFL)
c
c  Determine how many BOHM reactions are included in surf mech
C
         NBHM(IM) = 0
         DO 813 IS = 1, IISUR(IM)
            CALL SKIBHM(IS,ISKWRK(IMISK(IM)),IBHMFL)
            IF (IBHMFL.GT.0.AND.IBHMFL.LE.KKGAS) NBHM(IM)=NBHM(IM)+1
 813     CONTINUE
 820  CONTINUE
C
      ISOL = 0
C
C         RETURN HERE FOR A CONTINUATION PROBLEM
C
1000  CONTINUE
C
      ISOL = ISOL + 1
C
C   REWIND RECOVER FILE AND WRITE ITS HEADER
C
      REWIND (LRECOV)
      WRITE  (LRECOV) DHEAD
      WRITE  (LRECOV) VERSNN
C
C         CALL THE KEYWORD INPUT
C
      CALL PSRKEY (LIN, LOUT, MM, KK, KKGAS, KKSURF, KKBULK, NATJ,
     1             NSPHCH, NPHASE, NUMPSR, MAXPSR, NNSURF, NFSURF,
     2             NLSURF, NNBULK, NFBULK, NLBULK, ICKWRK, RCKWRK,
     3             LENISK, ISKWRK, LENSK, RSKWRK, LENCK, LENICK, KSYM,
     3             ATOM, PSYM, NCF, LENRGY, LFLRT, LEQUIV, LTDIFF,
     4             LRSTRT, LCNTUE, IPRNT, LSEN, LSENT, LSENG, ISEN,
     5             EPSS, EPST, EPSG, ABSOL, RELAT, LROP, IROP, LESTIM,
     6             ETCH, KKPHAS, KFIRST, KLAST, EPSR, TIN, XIN, T, 
     7             TSURF, X, EQUIV, PATM, PA, P, TAU, FLRT, V, Q,
     8             HTRN, TAMBNT, AREA, LHTRN, SDEN, FUEL, OXID, ADD, 
     9             STOICH, KPROD, KSP, ATOL, RTOL, ATIM, RTIM, IRETIR,
     *             NINIT, NJAC, ITJAC, NUMDT, DT1, NUMDT2, DT2, SFLR,
     1             UFAC, DFAC, DTMIN, DTMAX, LENIEQ, LENEQ, IEQWRK,
     2             EQWRK, IPVT, NIK, SCRTCH(1,1), A, NT, NTE, NYS, NY,
     3             NSDEN, NSDENS, WT, KEL, QEK, POWR, LENRGE, TEIN, TE,
     4             EIONSH, LNOFT, NPTS, NQLSE, QLSE, TQLSE, GFAC, SFAC,
     5             TIONP, LCONFN, SPOS, NPHMAX, KKMAX, NMAT, KKTOT,
     6             NSPTOT, AFRAC, MSYM, LTION, LELSH, BHMXI, ESHTH,
     7             LRFSH, LRFCUR, RFFREQ, RFAPAR, RFDPAR, LWALHB,
     8             VISC, THCOND, LSTEDY, ENDTIM, OLDTIM, LTERRN, KERR,
     9             EV2K, MUSD, LSCCM, SCCMIN, BIASPW)
      IF (KERR) RETURN
C
         IF (IPRNT .LT. 10) THEN
            LEVEL(1) = IPRNT
            LEVEL(2) = LEVEL(1)
         ELSE
            LEVEL(1) = IPRNT/10
            LEVEL(2) = IPRNT - 10*LEVEL(1)
            LEVEL(1) = MAX ( LEVEL(1), LEVEL(2) )
         ENDIF
C
C         SET THE SOLUTION BOUNDS
C
         BELOW(NT) = 200.
         ABOVE(NT) = 10000.0E0
         BELOW(NTE) = 200.
         ABOVE(NTE) = 1.E15
         DO 1200 K = 1, KKTOT
            BELOW (NYS+K) = SFLR
            ABOVE (NYS+K) = 1.01
1200     CONTINUE
         IF (NSPTOT .GT.0) THEN
            DO 1300 I = 1, NSPTOT
               BELOW (NSDENS+I) = -1.0E-10
               ABOVE (NSDENS+I) =  1.0E6
1300        CONTINUE
         ENDIF
C
C                            LOOP OVER PSRS IN SERIES
C
      DO 8500 IPSR = 1, NUMPSR
C
         TEMP(1) = TIN(IPSR)
         TEMP(2) = TEIN(IPSR)
         TEMP(3) = TIN(IPSR)
         IF (LENRGY(IPSR)) 
     1       CALL CKHMS (TEMP, ICKWRK, RCKWRK, HIN)
         CALL CKXTY (XIN(1,IPSR), ICKWRK, RCKWRK, YIN)
C
C    INITIAL SOLUTION GUESS FOR A RESTART:
C            - READ NEXT SOLUTION ESTIMATE FROM A DATA FILE
C              (can handle restarts involving multiple
C               psrs in series)
C
         IF (LRSTRT .AND. .NOT. RSTCNT) THEN
1500           CONTINUE
               READ (LREAD, ERR=2300, END=2400) ICHR
               IF (ICHR .EQ. DHEAD) THEN
                  READ (LREAD, ERR=2300)
                  GO TO 1500
               ELSEIF (ICHR .EQ. ICHEM) THEN
                  DO 1800 I = 1, 4
                     READ (LREAD, ERR=2300, END=2400)
1800              CONTINUE
                  DO 1850 IM = 1, NMAT
                     DO 1820 I = 1, 4
                        READ (LREAD, ERR=2300, END=2400)
 1820                CONTINUE
 1850             CONTINUE
                  GO TO 1500
               ELSEIF (ICHR .EQ. ISENSI) THEN
                  IISTOT = 0
                  DO 1900 IM = 1, NMAT
                     IISTOT = IISTOT + NNISUR(IM)
1900              CONTINUE
                  DO 2000 I = 1, (II+IISTOT)
                     READ  (LREAD, ERR=2300, END=2400)
2000              CONTINUE
                  GO TO 1500
               ELSEIF (ICHR .EQ. IROPRO) THEN
                  DO 2100 I = 1, KKTOT
                     READ (LREAD, ERR=2300, END=2400)
2100              CONTINUE
                  GO TO 1500
               ELSEIF (ICHR .NE. ISOLUT) THEN
                  WRITE (LOUT, *)
     1                   'ERROR: CHARACTER HEADER IS NOT KNOWN:',
     2                    ICHR
                  WRITE (LOUT, *)'          ',
     1                   'WILL CONTINUE WITHOUT RESTART FILE'
                  LRSTRT = .FALSE.
                  GO TO 2500
               ENDIF
               GO TO 2450
2300           CONTINUE
               WRITE (LOUT, *)'ERROR WHILE READING RESTART FILE'
               WRITE (LOUT, *)
     1               '       will continue without restart file'
               LRSTRT = .FALSE.
               GO TO 2500
2400           CONTINUE
               WRITE (LOUT, *)'END-OF-FILE REACHED WHILE READING ',
     1                        'THE RESTART FILE'
               WRITE (LOUT, *)
     1               '       will continue without restart file'
               LRSTRT = .FALSE.
               GO TO 2500
2450           CONTINUE
               READ (LREAD, ERR=2300, END=2400) NNNN, NIPSR, NNPSR,
     1                NNMAT, (NNSPHC,NNISUR(IM),IM=1,NNMAT)
               IF (NNNN .NE. NATJ) THEN
                  WRITE (LOUT, *)
     1            ' WARNING ERROR, INCOMPATIBLE RESTART FILE'
                  WRITE (LOUT, * ) 'Files''s natj,', NNNN,
     1            ' is not equal to problem''s natj,', NATJ
                  WRITE (LOUT, *)
     1            '              Will try to use it anyway'
               ENDIF
               IF (NNMAT .NE. NMAT) THEN
                  WRITE (LOUT,*)
     1            ' WARNING ERROR, INCOMPATIBLE RESTART FILE'
                  WRITE (LOUT,*) 'Files''s nmat,', NNMAT,
     1            ' is not equal to problem''s nmat,', NMAT
                  WRITE (LOUT, *)
     1            '              Will try to use it anyway'
               ENDIF
               IF (NIPSR .NE. IPSR) THEN
                  WRITE (LOUT, *)
     1            'RESTART WARNING: Restart file''s ipsr, ',
     2            NIPSR,' is not equal to current ipsr = ', IPSR
               ENDIF
               READ (LREAD, ERR=2300, END=2400) DUM1, DUM2, DUM3, DUM4,
     1                 DUM5,DUM6,DUM7, (DUM8, IM = 1, NNMAT),
     2                 DUM9, DUM10, LDUM,
     3                (DUM1, IM = 1, NNMAT)
               READ (LREAD, ERR=2300, END=2400) DUM1,
     1                 (SCRTCH(M,1), M=1,KKGAS)
               IF (KEL .NE. 0) READ (LREAD, ERR=2300, END=2400) DUM1, 
     1                               DUM2, DUM3, DUM4, DUM5
               READ (LREAD, ERR=2300, END=2400) 
     1                              (S(N), N=1,MIN(NNNN,NATJ))
               READ (LREAD, ERR=2300, END=2400) 
     1            (((DUM1,K=1,KK(IM)),I=1,MAX(1,NNISUR(IM))),
     2             IM=1,NNMAT), 
     2            ((DUM2, N=1, NPHASE(IM)),IM=1,NNMAT)
               READ (LREAD, ERR=2300, END=2400) OLDTIM
               IF (LENRGY(IPSR)) THEN
                 T(IPSR) = S(NT)
               ELSE
                 S(NT) = T(IPSR)
               ENDIF
               IF (KEL .EQ. 0) THEN
                 TE(IPSR) = T(IPSR)
                 TEIN(IPSR) = TIN(IPSR)
               ELSEIF (LENRGE(IPSR)) THEN
                 TE(IPSR) = S(NTE)
               ELSE
                 S(NTE) = TE(IPSR)
               ENDIF
2500           CONTINUE
         ENDIF
         IF (SPOS .GT. 0.0) THEN
            KSTOT = 0
            DO 3300 IM = 1, NMAT
               IF (KKSURF(IM) .GT. 0) THEN
                  DO 3000 KM = KFIRST(NFSURF(IM),IM),
     1                                  KLAST(NLSURF(IM),IM)
                     K = KM + KSTOT
                     S(NYS+K) = MAX(S(NYS+K), SPOS)
3000              CONTINUE
               ENDIF
               IF (KKBULK(IM) .GT. 0) THEN
                  DO 3100 KM = KFIRST(NFBULK(IM),IM), 
     1                                  KLAST(NLBULK(IM),IM)
                     K = KM + KSTOT
                     S(NYS+K) = MAX(S(NYS+K), SPOS)
3100              CONTINUE
               ENDIF
               DO 3200 K = 1, KKGAS
                  IF (KCHG(K).NE.0) THEN
                  ELSE
                     S(NYS+K) = MAX(S(NYS+K), SPOS)
                  ENDIF
3200           CONTINUE
               KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
3300        CONTINUE
         ENDIF

C
C INITIAL CONDITIONS FOR A "COLD" START
C
         IF (.NOT. LRSTRT) THEN
            CALL CKXTY (X, ICKWRK, RCKWRK, S(NY))
            S(NT) = T(IPSR)
            IF (KEL.NE.0) THEN
               S(NTE) = TE(IPSR)
            ELSE
               S(NTE) = T(IPSR)
            ENDIF
            KSTOT = 0
            IDTOT = 0
            DO 3700 IM = 1, NMAT
               IF (NPHASE(IM) .GT. 1) THEN
                  KF = KKGAS + 1
                  KL = KLAST(NPHASE(IM),IM)
                  DO 3500 KM = KF, KL
                     K = KM + KSTOT
                     S(NYS+K) = X(K)
3500              CONTINUE
                  KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
               ENDIF
               IF (NSPHCH(IM) .GT. 0) THEN
                  DO 3600 IDM = 1, NSPHCH(IM)
                     I = IDM + IDTOT
                     S(NSDENS+I) = SDEN(MAPPH(IDM,IM),IM)
     1                              / SDEN0(MAPPH(IDM,IM),IM)
3600              CONTINUE
               ENDIF
3700        CONTINUE
         ENDIF
C
C INITIAL CONDITIONS FOR A CONTINUATION RUN
C      - Change initial guess if user inputs solution estimates
C      - Currently solution estimates can only be given on
C        the first PSR.
C
         IF (RSTCNT) THEN
           IF (LENRGY(IPSR)) THEN
             S(NT) = T(IPSR)
           ELSE
             S(NT) = T(IPSR)
           ENDIF
           IF (KEL.NE.0) THEN
              S(NTE) = TE(IPSR)
           ELSE
              S(NTE) = T(IPSR)
           ENDIF
           IF (LESTIM(1,1) .AND. IPSR .EQ. 1) THEN
             CALL CKXTY (X, ICKWRK, RCKWRK, S(NY))
           ENDIF
           KSTOT = 0
           IDTOT = 0
           DO 4200 IM = 1, NMAT
             IF (NPHASE(IM) .GT. 1 .AND. IPSR .EQ. 1) THEN
               DO 4000 IPHASE = 2, NPHASE(IM)
                 IF (LESTIM(IPHASE,IM)) THEN
                   DO 3900 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                     K = KM + KSTOT
                     S(NYS+K) = X(K)
3900               CONTINUE
                 ENDIF
4000           CONTINUE
             ENDIF
             IF (NSPHCH(IM) .GT. 0) THEN
               DO 4100 IDM = 1, NSPHCH(IM)
                  I = IDM + IDTOT
                  S(NSDENS+I) = SDEN(MAPPH(IDM,IM),IM)
     1                           / SDEN0(MAPPH(IDM,IM),IM)
4100           CONTINUE
             ENDIF
             KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
             IDTOT = IDTOT + NSPHCH(IM)
4200       CONTINUE
         ENDIF
C
C*****precision > double
         CALL DCOPY (NATJ, S, 1, SN, 1)
C*****END precision > double
C*****precision > single
C         CALL SCOPY (NATJ, S, 1, SN, 1)
C*****END precision > single
C
C                PRINT OUT INFORMATION ABOUT THE KINETICS
C
         IF (LEVEL(1) .NE. 0) THEN
            WRITE (LOUT, 9000) KKGAS
            WRITE (LOUT, 9010) II
            WRITE (LOUT, 9007) NMAT
            DO 4900 IM = 1, NMAT
               WRITE (LOUT, 9008) MSYM(IM)
               WRITE (LOUT, 9001) KKSURF(IM), NNSURF(IM)
               IF (NNSURF(IM) .NE. 0) THEN
                  DO 4400 IPHASE = NFSURF(IM), NLSURF(IM)
                     WRITE (LOUT, 9002)PSYM(IPHASE,IM), 
     1                                 KKPHAS(IPHASE,IM)
4400              CONTINUE
               ENDIF
               WRITE (LOUT, 9005) KKBULK(IM), NNBULK(IM)
               IF (NNBULK(IM) .NE. 0) THEN
                  DO 4500 IPHASE = NFBULK(IM), NLBULK(IM)
                     WRITE (LOUT, 9002)PSYM(IPHASE,IM), 
     1                                 KKPHAS(IPHASE,IM)
4500              CONTINUE
               ENDIF
               WRITE (LOUT, 9011) IISUR(IM)
               IF (NSPHCH(IM) .GT. 0) THEN
                  WRITE (LOUT, 9020) NSPHCH(IM),
     1                   (PSYM(MAPPH(I,IM),IM), I = 1, NSPHCH(IM))
                  DO 4700 IPHASE = 1, NSPHCH(IM)
                     WRITE (LOUT, 9030) PSYM(MAPPH(IPHASE,IM),IM),
     1                                  NCON(MAPPH(IPHASE,IM),IM)
                     DO 4600 I = 1, IISUR(IM)
                        IF (NSTOIC(I,MAPPH(IPHASE,IM),IM) .NE. 0) THEN
                          CALL SKSYMR (I, LOUT, ISKWRK(IMISK(IM)), 
     1                            RSKWRK(IMRSK(IM)), CSWORK(IMCSK(IM)), 
     1                            LT, STRING, KERR)
                          WRITE (LOUT, 9040) I, STRING(1:LT)
                        ENDIF
4600                 CONTINUE
4700              CONTINUE
               ENDIF
               DO 4800 IPHASE = 1, NPHASE(IM)
                  IF (ETCH(IPHASE,IM,IPSR)) 
     1                     WRITE (LOUT, 9045) PSYM(IPHASE,IM)
4800           CONTINUE
4900        CONTINUE
            WRITE (LOUT, 9050)
         ENDIF
C
C              CALCULATE AND PRINT OUT INLET CONDITIONS
C
         TEMP(1) = S(NT)
         TEMP(2) = S(NTE)
         TEMP(3) = TIONP(IPSR)
         CALL CKRHOY (P(IPSR), TEMP, S(NY), ICKWRK, RCKWRK, RHO)
         IF (LFLRT(IPSR)) THEN
            TAU(IPSR) = RHO * V(IPSR) / FLRT(IPSR)
            IF (LEVEL(1) .GT. 0) WRITE (LOUT, 9200)
     1          FLRT(IPSR), RHO, V(IPSR), TAU(IPSR)
         ELSE
            FLRT(IPSR) = RHO * V(IPSR) / TAU(IPSR)
            IF (LEVEL(1) .GT. 0)
     1      WRITE (LOUT, 9201) TAU(IPSR), RHO, V(IPSR), FLRT(IPSR)
         ENDIF
         TEMP(1) = TIN(IPSR)
         TEMP(2) = TEIN(IPSR)
         TEMP(3) = TIN(IPSR)
         CALL CKRHOX (P(IPSR), TEMP, XIN(1,IPSR), ICKWRK,
     1                RCKWRK, RHOIN)
         CALL CKMMWX (XIN(1,IPSR), ICKWRK, RCKWRK, WTIN)
         VOLFIN = FLRT(IPSR) / RHOIN
         SCCM = VOLFIN * (298.15/TIN(IPSR)) * (PA(IPSR)/1.0) * 60.
         SLPM = SCCM * 1.0E-3
         EVOLTS = TEIN(IPSR)/EV2K
         IF (LEVEL(1) .GT. 0) THEN
            WRITE (LOUT, 9202) TIN(IPSR), PA(IPSR), RHOIN, WTIN,
     1                         FLRT(IPSR)/WTIN, VOLFIN
            WRITE (LOUT, 9204) SCCM, SLPM
            IF (KEL.NE.0) WRITE (LOUT, 9205) TEIN(IPSR), EVOLTS
            WRITE (LOUT, 9210)
            DO 5100 K = 1, KKGAS
               WRITE (LOUT, 9215)
     1         KSYM(K,1), XIN(K,IPSR),
     2         XIN(K,IPSR)*FLRT(IPSR)/WTIN,
     3         XIN(K,IPSR)*FLRT(IPSR)/WTIN * WT(K,1),
     4         XIN(K,IPSR)*FLRT(IPSR)/RHOIN,
     5         SCCM*XIN(K,IPSR), SLPM*XIN(K,IPSR)
5100        CONTINUE
            WRITE (LOUT, 9220)
         ENDIF
C
C        CALL THE DRIVER TO TWOPNT
C
         CALL PSRTWO (LOUT, LRECOV, IPSR, NUMPSR, NSPHCH, KK, KKGAS, 
     1                KKSURF, KKBULK, NPTS, NATJ, NPHASE, NNSURF, 
     2                NFSURF, NLSURF, NNBULK, NFBULK, NLBULK, KOCC, 
     3                LENRGY(IPSR), LFLRT(IPSR), LEQUIV, LTDIFF(1,IPSR), 
     4                LENTWP, LEVEL, ACTIVE, ETCH(1,1,IPSR), ICKWRK, 
     5                RCKWRK, ISKWRK, RSKWRK, MAPPH, WT, KSYM, PSYM, 
     6                KKPHAS, KFIRST, KLAST, ATIM, RTIM, ATOL, RTOL,
     7                TIN(IPSR), XIN(1,IPSR), YIN, HIN, ABSOL, RELAT, 
     8                ABOVE, BELOW, BUFFER, TWPWK, IPVT, NIK, SCRTCH, 
     9                T(IPSR), TSURF(1,IPSR), ACT, SDEN, SDEN0, SITDOT, 
     *                X, S, SN, F, FN, A, AA, DR, DC, DEN, GRATE, 
     1                IRETIR, NUMDT, NUMDT2, DT1, DT2, UFAC, DFAC,
     2                DTMIN, P(IPSR), PA(IPSR), TAU(IPSR), FLRT(IPSR),
     3                V(IPSR), Q(IPSR), HTRN(IPSR), TAMBNT(IPSR),
     4                LHTRN(IPSR), AREA(IPSR), EQUIV, NT, NTE, NYS, NY,
     5                NSDEN, NSDENS, NINIT, NJAC, ITJAC, DTMAX, KEL,
     6                KKION, KCHG, KION, QEK, POWR(IPSR), LENRGE(IPSR),
     7                TEIN(IPSR), TE(IPSR), ROP, ROPS, IEIMP, 
     8                EIONSH(1,IPSR), II, LNOFT, NQLSE, QLSE, TQLSE, 
     9                GFAC(IPSR), SFAC(IPSR), TIONP(IPSR), LCONFN, IEXC,
     *                REXC, SPOS, SDOTI, SITDTI, IISUR, NMAT, IISUR0, 
     1                KKMAX, NPHMAX, KKTOT, IMISK, IMRSK, IMCSK, 
     2                AFRAC(1,IPSR), MSYM, ITDEP, BHMXI(IPSR), FLXION, 
     3                FLXCOR, SDOTT, ESHTH(1,IPSR), LELSH(1,IPSR), 
     4                LRFSH(1,IPSR), LRFCUR, RFFREQ, RFAPAR, RFDPAR, 
     5                LWALHB(1,IPSR), VISC, THCOND, SDOT, LSTEDY, 
     6                ENDTIM, OLDTIM, LTERRN, ISOL, DEPRAT, SPUTTR, 
     7                NBHM, KERR, LOUTSH, NSPTOT, BOLTZ, AVOG, EV2K, 
     8                BIASPW(1,IPSR), HIM)
         IF (KERR) RETURN
C
C           PRINT FINAL SOLUTION, IF NO PRINTING THROUGHOUT ITERATION
C
         IF (IPRNT .EQ. 0) THEN
            LPRTIC = .TRUE.
            CALL PSPRNT (LOUT, NSPHCH, KK, KKGAS, KKSURF, KKBULK,
     1                   NPTS, NATJ, NPHASE, NNSURF, NFSURF, NLSURF,
     2                   NNBULK, NFBULK, NLBULK, KKPHAS, KFIRST, KLAST,
     3                   KOCC, KSYM, PSYM, ETCH(1,1,IPSR), LENRGY(IPSR),
     4                   LFLRT(IPSR), LTDIFF(1,IPSR), LEQUIV, LPRTIC,
     5                   HIN, YIN, EQUIV, PA(IPSR), P(IPSR), TAU(IPSR),
     6                   FLRT(IPSR), V(IPSR), AREA(IPSR), Q(IPSR),
     7                   HTRN(IPSR), TAMBNT(IPSR), LHTRN(IPSR),
     8                   TIN(IPSR), XIN(1,IPSR), T(IPSR), TSURF(1,IPSR), 
     9                   X, NIK, SCRTCH, SN, S, F, ICKWRK, RCKWRK, 
     *                   ISKWRK, RSKWRK, MAPPH, WT, ACT, SDEN, SDEN0,
     1                   SDOT, DEN, SITDOT, GRATE, NT, NTE, NYS, NY, 
     2                   NSDEN, NSDENS, KEL, KKION, KCHG, KION, QEK,
     3                   POWR(IPSR), LENRGE(IPSR),TEIN(IPSR), TE(IPSR),
     4                   ROP, ROPS, IEIMP, EIONSH(1,IPSR), II, NQLSE, 
     5                   QLSE, TQLSE, GFAC(IPSR), SFAC(IPSR), 
     6                   TIONP(IPSR), LCONFN, IEXC, REXC, SDOTI, SITDTI,
     7                   IISUR, NMAT, IISUR0, KKMAX, NPHMAX, KKTOT, 
     8                   IMISK, IMRSK,IMCSK, AFRAC(1,IPSR), MSYM, ITDEP, 
     9                   BHMXI(IPSR), FLXION, FLXCOR, SDOTT,
     *                   ESHTH(1,IPSR), LELSH(1,IPSR), LRFSH(1,IPSR), 
     1                   LRFCUR, RFFREQ, RFAPAR, RFDPAR, LWALHB(1,IPSR), 
     2                   VISC, THCOND, LSTEDY, ENDTIM, OLDTIM, LTERRN, 
     3                   DEPRAT, SPUTTR, NBHM, LOUTSH, BOLTZ, AVOG, 
     4                   EV2K, BIASPW(1,IPSR), HIM)
         ENDIF
C
C              WRITE TO DATA FILE, LSAVE, WHEN SOLUTION IS COMPLETE
C
         WRITE (LSAVE) ISOLUT
         WRITE (LSAVE) NATJ, IPSR, NUMPSR, NMAT, 
     1                 (NSPHCH(IM), IISUR(IM), IM=1,NMAT)
         WRITE (LSAVE) EQUIV, P(IPSR), TAU(IPSR), FLRT(IPSR),
     1                 V(IPSR), AREA(IPSR), Q(IPSR), 
     2                 (TSURF(IM,IPSR),IM=1,NMAT),
     3                 HTRN(IPSR), TAMBNT(IPSR), LHTRN(IPSR),
     4                 (AFRAC(IM,IPSR),IM=1,NMAT)
         WRITE (LSAVE) TIN(IPSR), (XIN(M,IPSR), M=1,KKGAS)
         IF (NQLSE.GT.0) CALL PLTEMP (NQLSE, S(NTE), TQLSE, QLSE, QLEX)
         IF (KEL.NE.0) WRITE (LSAVE) TEIN(IPSR), POWR(IPSR),TIONP(IPSR), 
     1       BHMXI(IPSR), 
     2    (EIONSH(IM,IPSR), ESHTH(IM,IPSR), BIASPW(IM,IPSR), IM=1,NMAT),
     3       QLEX
         IF (NSPTOT .GT. 0) THEN
            WRITE (LSAVE) (S(N), N=1,NATJ), 
     1                ((MAPPH(I,IM), I=1, MAX(1,NSPHCH(IM))),IM=1,NMAT)
         ELSE
            WRITE (LSAVE) (S(N), N=1,NATJ)
         ENDIF
C
         WRITE (LSAVE) 
     1    (((SDOTI(K,I,IM),K=1,KK(IM)),I=1,MAX(1,IISUR(IM))),IM=1,NMAT), 
     2    ((SITDTI(N,IM), N=1, NPHASE(IM)),IM=1,NMAT)
C
         SOLTIM = OLDTIM + ENDTIM
         WRITE (LSAVE) SOLTIM
C
C              WRITE DETAILED OUTPUT TO ASCII OUTPUT FILE
C
         IF (LEVEL(1) .EQ. 0) GOTO 8000
C
         TEMP(1) = S(NT)
         TEMP(2) = S(NTE)
         TEMP(3) = TIONP(IPSR)
         CALL CKRHOY (P(IPSR), TEMP, S(NY), ICKWRK, RCKWRK, RHO)
         IF (LFLRT(IPSR)) THEN
            TAU(IPSR) = RHO * V(IPSR) / FLRT(IPSR)
         ELSE
            FLRT(IPSR) = RHO * V(IPSR) / TAU(IPSR)
         ENDIF
         CALL CKYTX  (S(NY), ICKWRK, RCKWRK, X)
         CALL CKMMWX (X, ICKWRK, RCKWRK, WTOUT)
         CALL CKCDYP (P(IPSR), TEMP, S(NY), ICKWRK, RCKWRK,
     1                CDOT, DDOT)
         CALL CKWYP (P(IPSR), TEMP, S(NY), ICKWRK, RCKWRK, WDOT)
         DO 5200 K = 1, KKGAS
            WDOT(K) = WDOT(K) * GFAC(IPSR)
 5200    CONTINUE
         IDTOT = 0
         KSTOT = 0
         DO 5700 IM = 1, NMAT
            DO 5300 K = 1 , KKGAS
               ACT(K,IM) = X(K)
5300        CONTINUE
            IF (NNSURF(IM) .GT. 0) THEN
               IF (NSPHCH(IM) .GT. 0) THEN
                  DO 5400 IDM =1, NSPHCH(IM)
                     I = IDM + IDTOT
                     SDEN(MAPPH(IDM,IM),IM) = S(NSDENS+I) 
     1                                         * SDEN0(MAPPH(IDM,IM),IM)
5400              CONTINUE
               ENDIF
               DO 5500 KM = KFIRST(NFSURF(IM),IM), KLAST(NLSURF(IM),IM)
                  K = KM + KSTOT
                  ACT(KM,IM) = S(NYS+K)
                  X(K) = S(NYS+K)
5500           CONTINUE
            ENDIF
            IF (NNBULK(IM) .GT. 0) THEN
               DO 5600 KM = KFIRST(NFBULK(IM),IM),
     1                                KLAST(NLBULK(IM),IM)
                  K = KM + KSTOT
                  X(K) = S(NYS+K)
                  ACT(KM,IM) = X(K)
5600           CONTINUE
            ENDIF
            KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
            IDTOT = IDTOT + NSPHCH(IM)
5700     CONTINUE
C
C  DETERMINE the surface species production rates
C
         CALL PSSDOT (KKION, KEL, KION, KKGAS, KKTOT, KCHG, NMAT, 
     1                NPHASE, IISUR, KFIRST, KLAST, NFSURF, NLSURF, 
     2                NNSURF, KKSURF, NFBULK, NLBULK, NNBULK, KKBULK,
     3                IISUR0, KKMAX, NPHMAX, KK, BHMXI(IPSR), WT,
     4                EIONSH(1,IPSR), NIK, AFRAC, RFFREQ, 
     5                RFAPAR, RFDPAR, LRFSH(1,IPSR), LRFCUR, P, ACT, 
     6                SCRTCH(1,1), TIONP(IPSR), S(NTE), TSURF(1,IPSR), 
     7                SDEN, ISKWRK, RSKWRK, IMISK, IMRSK, FLXCOR, 
     8                FLXION, ESHTH(1,IPSR), LELSH(1,IPSR), SFAC(IPSR), 
     9                SDOTT, ROPS, SDOTI, SITDTI, SITDOT, SDOT, SHLOSS, 
     *                NBHM, LOUTSH, BOLTZ, AVOG, EV2K, BIASPW(1,IPSR), 
     1                AREA(IPSR), HIM(1,1))
C
         MASLOS = 0.0
         DO 6600 K = 1 , KKGAS
            MASLOS = MASLOS - SDOT(K)*WT(K,1)
6600     CONTINUE
         MASLOS = MASLOS*AREA(IPSR)
         FLRTOT = FLRT(IPSR) - MASLOS
C
C                   DETAILED OUTLET STREAM ANALYSIS
C
         IF (LFLRT(IPSR)) THEN
            WRITE (LOUT, 9400) FLRT(IPSR), FLRTOT, MASLOS,
     1                         RHO, V(IPSR), TAU(IPSR)
         ELSE
            WRITE (LOUT, 9401) TAU(IPSR), RHO, V(IPSR),
     1                         FLRT(IPSR), FLRTOT, MASLOS
         ENDIF
         VOLFOT = FLRTOT / RHO
         WRITE (LOUT, 9403) S(NT), PA(IPSR), RHO, WTOUT,
     1                      FLRTOT/WTOUT, VOLFOT
         SCCM = VOLFOT * (298.15/S(NT)) * (PA(IPSR)/1.0) * 60.
         SLPM = SCCM * 1.0E-3
         WRITE (LOUT, 9404) SCCM, SLPM
         IF (KEL.NE.0) THEN
            EVOLTS = S(NTE)/EV2K
            WRITE (LOUT, 9430) S(NTE), EVOLTS
         ENDIF
         WRITE (LOUT, 9410)
         TEMP(1) = S(NT)
         TEMP(2) = S(NTE)
         TEMP(3) = TIONP(IPSR)
         CALL CKRHOX (P(IPSR), TEMP, X, ICKWRK, RCKWRK, RHO)
         CALL CKXTCR (RHO, TEMP, X, ICKWRK, RCKWRK, SCRTCH(1,1))
         DO 6800 K = 1, KKGAS
            XNK = SCRTCH(K,1)*AVOG
            WRITE (LOUT, 9411)
     1      KSYM(K,1), X(K), XNK, X(K)*FLRTOT/WTOUT,
     2               X(K)*FLRTOT/WTOUT*WT(K,1),
     3               X(K)*FLRTOT/RHO, SCCM*X(K)
6800     CONTINUE
         WRITE (LOUT, 9412)
C
C                   DETAILED SPECIES BALANCE
C
         WRITE (LOUT, 9440)
         DO 6900 K = 1, KKGAS
            TNET = XIN(K,IPSR)*FLRT(IPSR)/WTIN
     1              + CDOT(K)*V(IPSR)
     2              + SDOT(K)*AREA(IPSR)
     3              - (X(K)*FLRTOT/WTOUT + DDOT(K)*V(IPSR))
            WRITE (LOUT, 9450)
     1         KSYM(K,1), XIN(K,IPSR)*FLRT(IPSR)/WTIN,
     2                  X(K)*FLRTOT/WTOUT,
     3                  CDOT(K)*V(IPSR), DDOT(K)*V(IPSR),
     4                  SDOT(K)*AREA(IPSR), TNET
6900     CONTINUE
         KSTOT = 0
         DO 7100 IM = 1, NMAT
            IF (KKSURF(IM) .NE. 0) THEN
               DO 7000 KM = KFIRST(NFSURF(IM),IM),KLAST(NLSURF(IM),IM)
                  K = KM + KSTOT
                  TNET =  SDOT(K)*AREA(IPSR)*AFRAC(IM,IPSR)
                  WRITE (LOUT, 9450) KSYM(KM,IM), 0.0, 0.0, 0.0, 0.0,
     1                            SDOT(K)*AREA(IPSR), TNET
7000           CONTINUE
            ENDIF
            IF (KKBULK(IM) .NE. 0) THEN
               DO 7050 KM = KFIRST(NFBULK(IM),IM),KLAST(NLBULK(IM),IM)
                  K = KM + KSTOT
                  TNET =  SDOT(K)*AREA(IPSR)*AFRAC(IM,IPSR)
                  WRITE (LOUT, 9450) KSYM(KM,IM), 0.0, 0.0, 0.0, 0.0,
     1                        SDOT(K)*AREA(IPSR)*AFRAC(IM,IPSR), TNET
7050           CONTINUE
            ENDIF
            KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
7100     CONTINUE
C
C           DETAILED ELEMENT BALANCE
C
         STRING = ' ELEMENT          INLET_FR      OUTLET_FR    '
         LSTRNG = 44
         DO 7350 IM = 1, NMAT
            IF (NNBULK(IM) .GT. 0) THEN
               DO 7300 IPHASE = NFBULK(IM), NLBULK(IM)
               WRITE(ADDON,9470) PSYM(IPHASE,IM)
               STRING = STRING(1:LSTRNG) // ADDON
               LSTRNG = LSTRNG + 14
7300        CONTINUE
         ENDIF
7350     CONTINUE
         STRING = STRING(1:LSTRNG) // ' TOTAL_NET '
         LSTRNG = LSTRNG + 11
         WRITE (LOUT, 9475) STRING
         DO 7800 M = 1, MM
            INLET = 0.0
            OUTLET = 0.0
            DO 7500 K = 1, KKGAS
               OUTLET = OUTLET + X(K)*NCF(M,K,1)
               INLET = INLET + XIN(K,IPSR)*NCF(M,K,1)
7500        CONTINUE
            INLET  = INLET  * FLRT(IPSR)   / WTIN
            OUTLET = OUTLET * FLRTOT / WTOUT
            TNET = INLET - OUTLET
            KSTOT = 0
            NBTOT = 0
            DO 7750 IM = 1, NMAT
               IF (NNBULK(IM) .GT. 0) THEN
                  NBTOT = NBTOT + NNBULK(IM)
                  DO 7700 IPHASE = NFBULK(IM), NLBULK(IM)
                     GRATED(IPHASE,IM) = 0.0
                     DO 7600 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                         K = KM + KSTOT
                         GRATED(IPHASE,IM) = GRATED(IPHASE,IM)
     1                      + SDOT(K)*AREA(IPSR)*AFRAC(IM,IPSR)
     2                      * NCF(M,KM,IM)
7600                 CONTINUE
                     TNET = TNET - GRATED(IPHASE,IM)
7700              CONTINUE
               ENDIF
               KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
7750        CONTINUE
            IF (NBTOT .GT. 0) THEN
               WRITE (LOUT, 9480) ATOM(M), INLET, OUTLET,
     1       ((GRATED(I,IM),I=NFBULK(IM), NLBULK(IM)),IM=1,NMAT), TNET
            ELSE
               WRITE (LOUT, 9480) ATOM(M), INLET, OUTLET, TNET
            ENDIF
7800     CONTINUE
         WRITE (LOUT, 9220)
C
C                 DO SENSITIVITY CALCULATIONS AND ROP CALCS
C
8000      CONTINUE
         IF (LSEN .OR. LTERRN) THEN
C
            CALL PSRSEN (LOUT, LSAVE, LRECOV, KK, II, IISUR, NSPHCH, 
     1                   KKGAS, KKSURF, KKBULK, NPTS, NATJ, NPHASE,  
     2                  NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK, 
     3                  KKPHAS, KFIRST, KLAST, KOCC, KSYM, PSYM, CCWORK, 
     4                   CSWORK, ETCH(1,1,IPSR), LENRGY(IPSR), 
     5                   LFLRT(IPSR), LTDIFF(1,IPSR), ABSOL, RELAT, DT, 
     6                   ABOVE, BELOW, WT, HIN, YIN, X, DEN, TIN(IPSR), 
     7                   P(IPSR), TAU(IPSR), FLRT(IPSR), V(IPSR),
     8                   AREA(IPSR), Q(IPSR), HTRN(IPSR), 
     9                   TAMBNT(IPSR), LHTRN(IPSR), T(IPSR),
     *                   TSURF(1,IPSR), S, SN, F, FN, A, ACT, SDEN,
     1                   SDEN0, SITDOT, SDOT, SCRTCH(1,8), SCRTCH(1,9), 
     2                   GRATE, GRATED, DSKDPL, LSEN, LSENT, LSENG, 
     3                   ISEN, EPSS, EPST, EPSG, ICKWRK, RCKWRK, ISKWRK,
     5                   RSKWRK, MAPPH, IPVT, NIK, SCRTCH, NT, NTE, NYS,
     6                   NY, NSDEN, NSDENS, KEL, KKION, KCHG, KION, 
     7                   QEK, POWR(IPSR), LENRGE(IPSR), TEIN(IPSR),
     8                   TIN(IPSR), ROP, ROPS, IEIMP, EIONSH(1,IPSR), 
     9                   NQLSE, QLSE, TQLSE, GFAC(IPSR), SFAC(IPSR),
     *                   TIONP(IPSR), LCONFN, IEXC, REXC, SDOTI, SITDTI, 
     1                   KKTOT, NMAT, IISUR0, IISMAX, KKMAX, NPHMAX, 
     2                   IMISK, IMRSK, IMCSK, AFRAC(1,IPSR),MSYM, ITDEP, 
     3                   BHMXI(IPSR), FLXION, FLXCOR, SDOTT, 
     4                   ESHTH(1,IPSR), LELSH(1,IPSR), LRFSH(1,IPSR), 
     5                   LRFCUR, RFFREQ, RFAPAR, RFDPAR, LWALHB(1,IPSR), 
     6                   VISC, THCOND, LTERRN, ANISOT, NBHM, KERR, 
     7                   LOUTSH, BOLTZ, AVOG, EV2K, BIASPW(1,IPSR), HIM)
C
            IF (KERR) RETURN
            IF (LEVEL(2) .GT. 0) WRITE (LOUT,'(/A/)')
     1      ' SENSITIVITY CALCULATION COMPLETE'
C
         ENDIF
C
         IF (LROP) THEN
C
            CALL PSRROP (LOUT, LSAVE, LRECOV, KK, II, IISUR, NSPHCH, 
     1                   KKGAS, KKSURF, KKBULK, NPHASE, NNSURF, NFSURF, 
     2                   NLSURF, NNBULK, NFBULK, NLBULK, NATJ, KKPHAS, 
     3                   KFIRST, KLAST, KSYM, PSYM, P(IPSR), 
     4                   TSURF(1,IPSR), V(IPSR), AREA(IPSR), S,
     5                   X, ACT, SDEN, SDEN0, ICKWRK, RCKWRK, CCWORK,
     5                   ISKWRK, RSKWRK, CSWORK, MAPPH, LROP, IROP,
     7                   EPSR, ROP, ROPS, NIK, SCRTCH(1,3),
     8                   SCRTCH(1,4), SCRTCH(1,5), SCRTCH(1,6),
     9                   NT, NTE, NYS, NY, NSDEN, NSDENS, TIONP(IPSR),
     *                   LCONFN, KION, KCHG, KKTOT, NMAT, KKMAX, NPHMAX,
     1                   IMISK, IMRSK, IMCSK, IISUR0, MSYM,
     2                   AFRAC(1,IPSR), BHMXI(IPSR), FLXION, FLXCOR, 
     3                   SFAC(IPSR), SDOTT, SCRTCH(1,7), EIONSH(1,IPSR), 
     4                   WT, KEL, KKION, RFFREQ, RFAPAR, RFDPAR,
     5                   LRFSH(1,IPSR), LRFCUR, ESHTH(1,IPSR), 
     6                   LELSH(1,IPSR), SDOT, SITDOT,
     6                   SDOTI, SITDTI, KERR, LOUTSH, BOLTZ, AVOG, EV2K,
     7                   BIASPW(1,IPSR), HIM)
            IF (KERR) RETURN
C
            IF (LEVEL(2) .GT. 0) WRITE (LOUT,'(/A/)')
     1      ' RATE-OF-PRODUCTION CALCULATION COMPLETE'
C
         ENDIF
C
C              SAVE SOLUTION FOR USE IN CONT. PROB.
         DO 8100 I = 1, NATJ
            SSAVE(I, IPSR) = S(I)
8100      CONTINUE
C                  FIND THE INITIAL CONDITIONS FOR PSRS IN SERIES
         IF (IPSR .LT. NUMPSR) THEN
            TIN(IPSR+1)   = S(NT)
            TEIN(IPSR+1)  = S(NTE)
            FLRT(IPSR+1)  = FLRTOT
            LFLRT(IPSR+1) = .TRUE.
            CALL CKYTX  (S(NY), ICKWRK, RCKWRK, X)
            DO 8200 K = 1, KKGAS
               XIN(K,IPSR+1) = X(K)
8200        CONTINUE
C
C            FOR THE INITIAL GUESS TO THE SOLUTION, USE SAVED SOLUTION
C              VECTOR FOR NEXT PSR, IF AVAILABLE
C              (If it is not available, the current solution, S, for
C               the current PSR is used for the initial guess for the
C               solution to the next PSR)
C
            IF (RSTCNT) THEN
               DO 8300 I = 1, NATJ
                  S(I) = SSAVE(I,IPSR+1)
8300            CONTINUE
            ENDIF
            WRITE (LOUT, 9810) IPSR + 1
         ENDIF

8500  CONTINUE
C
C
C  WRITE DEP RATE AND ETCH RATE DATA TO TERRAIN FILE 'terrain.#' 
C
      IF (LTERRN .AND. ISOL .LE. 9) THEN
         REDEP = 0.5
         TRNFIL = 'terrain.'
         TRNFIL(9:9) = CHNUM(ISOL+1:ISOL+1)
         OPEN (UNIT=18, FILE= TRNFIL, FORM='FORMATTED',
     1         STATUS='UNKNOWN')
         REWIND(18)
         WRITE (18, 9900) DEPRAT, ANISOT, REDEP, SPUTTR
         CLOSE (UNIT=18)
      ENDIF

C            CHECK FOR CONTINUATION
C
      IF (LCNTUE) THEN
C
         WRITE (LOUT,'(/////)')
         DO 8600 L = 1, 5
            WRITE (LOUT, *)
     1    ' ////////////////// CONTINUING TO NEW PROBLEM /////////////'
8600     CONTINUE
         WRITE (LOUT,'(/////)')
C
C            READ IN SAVED SOLUTION FOR PSR #1
         DO 8700 I = 1, NATJ
            S(I) = SSAVE(I,1)
8700     CONTINUE
C
C  Make sure that variables, which may or may not be changed by
C user input in the function, PSRKEY, agree with the calculated
C solution for all PSRs
C - (They must be reread back into the solution vector, S, after
C    the call to PSRKEY)
C
         IDMTOT = 0
         DO 8850 IM = 1, NMAT
            IF (NSPHCH(IM) .GT. 0) THEN
              DO 8800 IDM = 1, NSPHCH(IM)
                 I = IDM + IDMTOT
                 SDEN(MAPPH(IDM,IM),IM) = 
     1                      S(NSDENS+I)*SDEN0(MAPPH(IDM,IM),IM)
8800          CONTINUE
              IDMTOT = IDM + NSPHCH(IM)
            ENDIF
8850     CONTINUE
         DO 8900 IPSR = 1, NUMPSR
           T(IPSR) = SSAVE(NT,IPSR)
           TE(IPSR) = SSAVE(NTE,IPSR)
8900     CONTINUE
C
C  Set flags to indicate that that it is a continuation problem
C and that there exists a good approximation to the solution.
         RSTCNT = .TRUE.
         LRSTRT = .TRUE.
         GO TO 1000
C
      ENDIF
C
      RETURN
C-----------Kinetics printout subsection----------------------------
9000  FORMAT(/120('=')/5X,'KINETICS MECHANISM:'/
     1 8X,'Total number of gas phase species     = ',I4)
9001  FORMAT(
     1 10X,'Total number of surface phase species = ',I4/
     1 10X,'Total number of surface phases = ',I4)
9002  FORMAT(
     1 15X,'Number of species in phase, ',A10,' = ',I4)
9005  FORMAT(
     1 10X,'Total number of bulk phase species = ',I4/
     1 10X,'Total number of bulk phases = ',I4)
9007  FORMAT(8X,'Total number of surface materials = ',I4/1X)
9008  FORMAT(8X,'Material:  ',A10)
9010  FORMAT(
     1 8X,'Total number of gas phase reactions     = ',I4)
9011  FORMAT(
     1 10X,'Total number of surface phase reactions = ',I4)
9020  FORMAT(
     1 10X,'Total number of surface phases whose site densities are ',
     1  'not conserved = ', I4/15X,'their names = ',
     1  8(A10,2X) )
9030  FORMAT(15X,'For surface phase, ', A10,
     1   ', number of non-conserving surface site surface reactions = ',
     1   I4)
9040  FORMAT(20X,'Surface reaction # ',I4,' --> ',A60)
9045  FORMAT(10X,'Bulk phase, ',A10,' is assumed to be undergoing an',
     1    ' etching reaction.'/15X,'Its composition will be fixed at',
     1    ' the initial estimate.')
9050  FORMAT(120('=')/)
C-----------Inlet Conditions printout subsection-------------------
9200  FORMAT(/120('=')/5X,'INLET CONDITIONS:'/
     1 10X,'Inlet mass flow rate = ',G11.3,' gm/sec'/
     1 15X,'(which based on an estimated reactor density = ',
     1 G11.3,' gm/cm**3'/
     1 15X,'        and  on a reactor volume  = ',
     1 G11.3,' cm**3'/
     1 15X,'       produces an estimated residence time) = ',
     1 G11.3,' sec'//)
9201  FORMAT(/120('=')/5X,'INLET CONDITIONS:'/
     1 10X,'Residence time in reactor = ',G11.3,' sec'/
     1 15X,'(which based on an estimated reactor density = ',
     1 G11.3,' gm/cm**3'/
     1 15X,'        and on a reactor volume  = ',
     1 G11.3,' cm**3'/
     1 15X,' produces an estimated inlet mass flow rate) = ',
     1 G11.3,' gm/sec'//)
9202  FORMAT(
     1 10X,'Inlet temperature = ',G11.5,' Kelvin'/
     1 10X,'Inlet pressure (assumed equal to reactor pressure) = ',
     1 G11.3,' atm'/
     1 10X,'Inlet density = ',G11.5,' gm/cm**3'/
     1 10X,'Inlet mean molecular weight = ',G11.5,' gm/mole'/
     1 10X,'Inlet molar flow rate      = ',G11.5,' moles/sec'/
     1 10X,'Inlet volumetric flow rate = ',G11.5,' cm**3/sec'/
     1 10X,'     (based on reactor pressure and inlet temperature)')
9204  FORMAT(
     1 10X,'                           = ',G11.5,' SCCM'/
     1 10X,'                           = ',G11.5,' SLPM'/)
9205  FORMAT(
     1 10X,'Inlet electron temperature = ',G11.5,' Kelvin'/
     1 10X,'                           = ',G11.5,' eV')
9210  FORMAT(5X,'INLET CONDITIONS FOR GAS PHASE MOLECULAR SPECIES:'//
     1 3X,'Species        mole_frac     moles/sec       gm/sec   ',
     1 '    cm**3/sec  ',
     1 '     SCCM            SLPM'/102('-'))
9215  FORMAT(
     1 4X,A10,3X,G11.5,3X,G11.5,3X,G11.5,3X,G11.5,3X,G11.5,3X,G11.5)
9220  FORMAT(/120('=')//)
C--------------------OUTPUT CONDITIONS FORMAT SECTION---------------
9400  FORMAT(/120('=')/5X,'OUTLET CONDITIONS:'/
     1 10X,'Specified inlet mass flow rate = ',G11.3,' gm/sec'/
     1 10X,'Outlet mass flow rate = ',G11.3,' gm/sec'/
     1 10X,'Rate of Mass Loss to the walls = ',G11.3,' gm/sec'/
     1 15X,'(which based on an reactor density = ',
     1 G11.3,' gm/cm**3'/
     1 15X,'        and  on a reactor volume  = ',
     1 G11.3,' cm**3'/
     1 15X,'       produces a residence time) = ',
     1 G11.3,' sec'//)
9401  FORMAT(/120('=')/5X,'OUTLET CONDITIONS:'/
     1 10X,'Specified residence time in reactor = ',G11.3,' sec'/
     1 15X,'(which based on a reactor density = ',
     1 G11.3,' gm/cm**3'/
     1 15X,'        and on a reactor volume  = ',
     1 G11.3,' cm**3'/
     1 15X,' produces an inlet mass flow rate) = ',
     1 G11.3,' gm/sec'/
     1 10X,'Outlet mass flow rate = ',G11.3,' gm/sec'/
     1 10X,'Rate of Mass Loss to the walls = ',G11.3,' gm/sec'//)
9403  FORMAT(
     1 10X,'Outlet and reactor temperature = ',G11.5,' Kelvin'/
     1 10X,'Outlet and reactor pressure = ',
     1 G11.3,' atm'/
     1 10X,'Outlet and reactor density = ',G11.5,' gm/cm**3'/
     1 10X,'Outlet and reactor mean molecular weight = ',
     1     G11.5,' gm/mole'/
     1 10X,'Outlet molar flow rate      = ',G11.5,' moles/sec'/
     1 10X,'Outlet volumetric flow rate = ',G11.5,' cm**3/sec'/
     1 10X,'     (based on reactor pressure and temperature)')
9404  FORMAT(
     1 10X,'                           = ',G11.5,' SCCM'/
     1 10X,'                           = ',G11.5,' SLPM'/)
9410  FORMAT(5X,'OUTLET CONDITIONS FOR GAS PHASE MOLECULAR SPECIES:'//
     1 3X,'Species        mole_frac     #/cm3        moles/sec       ',
     1 'gm/sec       cm**3/sec  ',
     1 '     SCCM'/102('-'))
9411  FORMAT(
     1 4X,A10,3X,G11.5,3X,G11.5,3X,G11.5,3X,G11.5,3X,G11.5,3X,G11.5)
9412  FORMAT(/120('=')//)
9440  FORMAT(/120('=')//25X,'DETAILED SPECIES BALANCE'/
     1       35X,'(all rates are in moles per sec)'//
     1 '   SPECIES           INLET_FR    OUTLET_FR  GAS_PROD_RATE',
     1 '  GAS_DEST_RATE   SURF_NET_PROD    TOTAL_NET'
     1 /130('-'))
9430  FORMAT(
     1 10X,'Outlet and reactor electron temperature = ',G11.5,' Kelvin'/
     1 10X,'                                        = ',G11.5,' eV'/)
9450  FORMAT(1X, A16, 3X, 6(1PG11.3, 3X))
9460  FORMAT(//A/120('-'))
9470  FORMAT(' ',A8,'_GR  ')
9475  FORMAT(//120('=')//25X,'DETAILED ELEMENT BALANCES'/
     1         35X,'(all rates are in moles per sec)'//
     1         A/120('-'))
9480  FORMAT(1X, A14,3X,8(1PG11.3,3X))
C-------------------------------------------------------------------
9810  FORMAT(///////////120('=')/120('=')/40X,
     1       'START PSR CALCULATION FOR NEXT PSR, ',I2/50X,
     1 '(assume inlet is the same as the outlet of the last psr)',
     1 /120('=')/120('=')////)
 9900 FORMAT (1X,'mach.depo hdp name=dsm material=Oxide rate=',f9.5/
     1        '+ anisot=',f5.3,' rdep.ra=',f5.4,' mill.ra=',f9.5/
     2        '+ sput.c1=5.5 sput.c2=-6 sput.c4=1.5')
      END
C
C--------------------------------------------------------------------
C
      SUBROUTINE PSRTWO (LOUT, LRECOV, IPSR, NUMPSR, NSPHCH, KK, KKGAS,
     1                   KKSURF, KKBULK, NPTS, NATJ, NPHASE, NNSURF,
     2                   NFSURF, NLSURF, NNBULK, NFBULK, NLBULK, KOCC,
     3                   LENRGY, LFLRT, LEQUIV, LTDIFF, LENTWP, LEVEL,
     4                   ACTIVE, ETCH, ICKWRK, RCKWRK, ISKWRK, RSKWRK,
     5                   MAPPH, WT, KSYM, PSYM, KKPHAS, KFIRST, KLAST,
     6                   ATIM, RTIM, ATOL, RTOL, TIN, XIN, YIN, HIN,
     7                   ABSOL, RELAT, ABOVE, BELOW, BUFFER, TWPWK,
     8                   IPVT, NIK, SCRTCH, T, TSURF, ACT, SDEN, SDEN0,
     9                   SITDOT, X, S, SN, F, FN, A, AA, DR, DC, DEN,
     *                   GRATE, IRETIR, NUMDT, NUMDT2, DT1, DT2, UFAC,
     1                   DFAC, DTMIN, P, PA, TAU, FLRT, V, Q, HTRN,
     2                   TAMBNT, LHTRN, AREA, EQUIV, NT, NTE, NYS,
     3                   NY, NSDEN, NSDENS, NINIT, NJAC, ITJAC, DTMAX,
     4                   KEL, KKION, KCHG, KION, QEK, POWR, LENRGE,
     5                   TEIN, TE, ROP, ROPS,IEIMP, EIONSH, II, LNOFT,
     6                   NQLSE, QLSE, TQLSE, GFAC, SFAC, TIONP, LCONFN,
     7                   IEXC, REXC, SPOS, SDOTI, SITDTI, IISUR, NMAT,
     8                   IISUR0, KKMAX, NPHMAX, KKTOT, IMISK, IMRSK,
     9                   IMCSK, AFRAC, MSYM, ITDEP, BHMXI, FLXION,
     *                   FLXCOR, SDOTT, ESHTH, LELSH, LRFSH,
     1                   LRFCUR, RFFREQ, RFAPAR, RFDPAR, LWALHB, VISC,
     2                   THCOND, SDOT, LSTEDY, ENDTIM, OLDTIM, LTERRN,
     3                   ISOL, DEPRAT, SPUTTR, NBHM, KERR, LOUTSH,
     4                   NSPTOT, BOLTZ, AVOG, EV2K, BIASPW, HIM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NSPHCH(*), KK(NMAT), KKSURF(NMAT), KKBULK(NMAT), 
     1          NPHASE(NMAT), NNSURF(NMAT), NFSURF(NMAT),
     1          NLSURF(NMAT), NNBULK(NMAT), NFBULK(NMAT),
     2          NLBULK(NMAT), IISUR(NMAT), IMISK(NMAT), IMRSK(NMAT), 
     5          IMCSK(NMAT)
      DIMENSION LEVEL(2), ICKWRK(*), ISKWRK(*), MAPPH(NPHMAX,NMAT),
     1          IPVT(NATJ), KOCC(KKMAX,NMAT), KKPHAS(NPHMAX,NMAT), 
     2          KFIRST(NPHMAX,NMAT), KLAST(NPHMAX,NMAT), RCKWRK(*), 
     3          RSKWRK(*), WT(KKMAX,NMAT), XIN(KKTOT),
     4          YIN(KKMAX), HIN(KKMAX), ABOVE(NATJ), BELOW(NATJ),
     5          BUFFER(NATJ), TWPWK(*), ACT(KKMAX,NMAT), 
     6          SDEN(NPHMAX,NMAT), SDEN0(NPHMAX,NMAT),
     6          SITDOT(NPHMAX,NMAT), SCRTCH(NIK,*), X(KKTOT),
     7          S(NATJ), SN(NATJ), F(NATJ), FN(NATJ), A(NATJ,NATJ),
     9          AA(NATJ,NATJ), DR(NATJ), DC(NATJ), ITDEP(*),
     *          DEN(KKMAX,NMAT), GRATE(NPHMAX,NMAT), KCHG(KKMAX), 
     1          KION(*), QEK(KKGAS), ROP(*), ROPS(IISUR0,NMAT),
     1          IEIMP(*), QLSE(NPTS), TQLSE(NPTS), IEXC(*),
     2          REXC(*),  AFRAC(NMAT), SDOTI(KKMAX,IISUR0,NMAT),
     3          SITDTI(NPHMAX,NMAT), EIONSH(NMAT), FLXION(*),
     4          FLXCOR(*), SDOTT(KKTOT), ESHTH(NMAT), SDOT(KKTOT),
     5          TSURF(NMAT), NBHM(NMAT), BIASPW(NMAT),
     6          HIM(KKMAX,NMAT)
C
      CHARACTER*16 KSYM(KKMAX,NMAT), PSYM(NPHMAX,NMAT), MSYM(NMAT)
      CHARACTER*18 DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO
      CHARACTER VERSIO*80, REPORT*80, SIGNAL*80
C
      LOGICAL ACTIVE(NATJ), MARK(1), LENRGY, LFLRT, LEQUIV, 
     1        LTDIFF(NMAT), ETCH(NPHMAX,NMAT), LHTRN, LPRTIC,ERROR,
     2        KERR, LTIME, ENERGY, LWRITE, XLINPK, FINDJ,
     4        LENRGE, EENRGY, LNOFT, LCONFN, LELSH(NMAT), 
     5        LRFSH(NMAT), LRFCUR, LWALHB(NMAT), LSTEDY, LTERRN

C                CHARACTER DATA FILE COMMON BLOCK
C
      COMMON /CDATAF/ VERSNN, DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO
C
      INTEGER  PMAX, POINTS, CALL, CALLS
      SAVE XLINPK
C
C                             EXTERNALS
C
      EXTERNAL PSRFUN, PSRJAC, PSPRNT
C*****precision > double
      EXTERNAL DCOPY, DGECO, DGESL
C*****END precision > double
C*****xlinpk double precision
      EXTERNAL XDGECO, XDGESL
C*****END xlinpk double precision
C*****precision > single
C      EXTERNAL SCOPY, SGECO, SGESL
C*****END precision > single
C*****xlinpk single precision
C      EXTERNAL XSGESL, XSGECO
C*****END xlinpk single precision
C
      DATA PMAX/1/, POINTS/1/, XLINPK /.FALSE./
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C INPUT-
C   LOUT   - UNIT FOR PRINTED OUTPUT
C   LRECOV - UNIT TO WHICH THE RECOVER FILE IS WRITTEN
C   KKGAS  - NUMBER OF GAS PHASE SPECIES
C   KKSURF - NUMBER OF SURFACE SPECIES
C   KKBULK - NUMBER OF BULK PHASE SPECIES
C   NPHASE - NUMBER OF PHASES
C   KK     - TOTAL NUMBER OF CHEMICAL SPECIES.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES.  NATJ=KK+2.
C   KKPHAS - VECTOR OF THE NUMBER OF SPECIES IN EACH PHASE
C              DIMENSION KKPHAS AT LEAST NPHASE
C   KFIRST - VECTOR OF POINTERS TO THE FIRST SPECIES IN A PHASE
C              DIMENSION KFIRST AT LEAST NPHASE
C   KLAST  - VECTOR OF POINTERS TO THE LAST SPECIES IN EACH PHASE
C              DIMENSION KLAST AT LEAST NPHASE
C   LENRGY - IF LENRGY=.TRUE. THEN THE ENERGY EQUATION IS TO BE
C            SOLVED, OTHERWISE A FIXED TEMPERATURE IS USED.
C   LFLRT  - IF LFLRT=.TRUE. THEN THE MASS FLOW RATE IS SPECIFIED AND
C              THE RESIDENCE TIME HAS TO BE CALCULATED.
C            IF LFLRT=.FALSE. THEN THE RESIDENCE TIME IS SPECIFIED AND
C              THE MASS FLOW RATE HAS TO BE CALCULATED.
C   LEQUIV - IF LEQUIV=.TRUE. THEN THE INLET COMPOSITION IS DEFINED BY
C              THE FUEL EQUIVALENCE RATIO, THE FUEL AND OXIDIZER
C              COMPOSITIONS AND THE PRODUCTS.
C            IF LEQUIV=.FALSE. THEN THE INLET COMPOSITION IS DEFINED
C              DIRECTLY BY THE USER.
C            LEQUIV AND EQUIV ARE USED ONLY FOR PRINTING.  THEREFORE, IF
C            PSRTWO IS CALLED OUTSIDE OF PSR, THEN EQUIV MUST BE GIVEN
C            IF LEQUIV=.TRUE.
C   LHTRN  - IF LHTRN=.TRUE. THEN THE HEAT LOSS FROM THE REACTOR
C            IS GIVEN IN TERMS OF A HEAT TRANSFER COEFFICIENT, HTRN,
C            AND AMBIENT TEMPERATURE, TAMBNT
C   LTDIFF - IF LTDIFF=.TRUE. THEN THE SURFACE TEMPERATURE IS DIFFERENT
C              THAN THE BULK TEMPERATURE FOR A MATERIAL.
C   LEVEL  - LEVEL(1) CONTROLS THE PRINTED OUTPUT FROM TWOPNT AND
C            LEVEL(2) CONTROLS THE PRINTING FROM PSR.  LEVEL(1) MUST
C            ALWAYS BE GREATER THAN OR EQUAL TO LEVEL(2).  LEVEL MAY
C            HAVE VALUES OF 0, 1, 2, OR 3.
C              DIMENSION LEVEL(2).
C   ACTIVE - ACTIVE(*) IS A LOGICAL ARRAY THAT CONTROLS WHICH VARIABLES
C            TWOPNT WILL ATTEMPT ADAPTIVE MESHING.  FOR THE PSR PROBLEM
C            THERE IS NO MESH, THEREFORE ALL ACTIVE(*)=.FALSE.
C              LOGICAL ACTIVE(KK).
C   WT     - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
C              CGS UNITS - GM/MOLE
C              DIMENSION WT(*) AT LEAST KK.
C   KSYM   - CHEMKIN SPECIES NAMES.
C              DIMENSION KSYM AT LEAST KK.
C   PSYM   - SURFACE CHEMKIN PHASE NAMES
C              DIMENSION PSYM AT LEAST NPHASE
C   ATIM   - ABSOLUTE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION
C            AS USED FOR THE TIME STEPS.
C   ATOL   - ABSOLUTE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION.
C   RTIM   - RELATIVE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION
C            AS USED FOR THE TIME STEPS.
C   RTOL   - RELATIVE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION.
C   TIN    - INLET TEMPERATURE.
C              CGS UNITS - K
C   XIN    - ARRAY OF INLET SPECIES MOLE FRACTIONS. USED FOR PRINTING.
C              DIMENSION XIN(*) AT LEAST KK.
C   YIN    - ARRAY OF INLET SPECIES MASS FRACTIONS.
C              DIMENSION YIN(*) AT LEAST KK.
C   HIN    - ARRAY OF SPECIES ENTHALPIES AT INLET TEMPERATURE TIN.
C              CGS UNITS - ERGS/G
C              DIMENSION HIN(*) AT LEAST KK.
C   HTRN   - HEAT TRANSFER COEFFICIENT OF THE REACTOR
C   ABSOL  - ABSOLUTE PERTURBATION FOR COMPUTING JACOBIAN.
C   RELAT  - RELATIVE PERTURBATION FOR COMPUTING JACOBIAN.
C   ABOVE  - ARRAY OF UPPER BOUNDS FOR THE DEPENDENT VARIABLES.  USED
C            BY TWOPNT TO CONTROL THE DAMPING.  ABOVE HAS THE SAME
C            STRUCTURE AS THE SOLUTION VECTOR S(*).  ABOVE(*) SHOULD BE
C            SET TO VALUES THAT ARE ABOVE ACCEPTABLE SOLUTION VALUES.
C              DIMENSION ABOVE(*) AT LEAST NATJ.
C   BELOW  - ARRAY OF LOWER BOUNDS FOR THE DEPENDENT VARIABLES.  USED
C            BY TWOPNT TO CONTROL THE DAMPING.  BELOW HAS THE SAME
C            STRUCTURE AS THE SOLUTION VECTOR S(*).  BELOW(*) SHOULD BE
C            SET TO VALUES THAT ARE BELOW ACCEPTABLE SOLUTION VALUES.
C              DIMENSION BELOW(*) AT LEAST NATJ.
C   T      - THE TEMPERATURE. FOR FIXED-TEMPERATURE PROBLEMS THIS
C            IS THE USER-SPECIFIED TEMPERATURE.
C              CGS UNITS - K
C   TAMBNT - TEMPERATURE OF THE AMBIENT
C   TSURF  - THE SURFACE TEMPERATURE.  USED IF LTDIFF=.TRUE. 
C   S      - DEPENDENT VARIABLE ARRAY.  ON INPUT THE SOLUTION ESTIMATE
C            IS IN S.  THE TEMPERATURE IS STORED IN T=S(NT), AND THE
C            MASS FRACTIONS ARE IN Y(K)=S(NYS+K).
C              DIMENSION S(*) AT LEAST NATJ.
C   NUMDT  - THE NUMBER OF TIME STEPS TO BE TAKEN UPON FAILURE OF A
C            NEWTON ITERATION.  NUMDT IS IN EFFECT WHEN THE TEMPERATURE-
C            FIXED PROBLEM IS BEING SOLVED.
C   NUMDT2 - THE NUMBER OF TIME STEPS TO BE TAKEN UPON FAILURE OF A
C            NEWTON ITERATION.  NUMDT2 IS IN EFFECT WHEN THE ENERGY
C            EQUATION IS BEING SOLVED.
C   DT1    - THE VALUE OF THE TIMESTEP THAT IS TAKEN WHEN THE NEWTON
C            ITERATION FAILS IN THE TEMPERATURE-FIXED PROBLEM.
C              CGS UNITS - SEC
C   DT2    - THE VALUE OF THE TIMESTEP THAT IS TAKEN WHEN THE NEWTON
C            ITERATION FAILS WHEN THE ENERGY EQUATION IS BEING SOLVED.
C              CGS UNITS - SEC
C   P      - THE PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   PA     - THE PRESSURE.
C              CGS UNITS - ATMOSPHERES
C   TAU    - THE NOMINAL RESIDENCE TIME OF THE REACTOR
C              CGS UNITS - SEC
C   FLRT   - THE MASS FLOW RATE.
C              CGS UNITS - GM/SEC
C   V      - THE VOLUME OF THE REACTOR
C              CGS UNITS - CM**3
C   AREA   - SURFACE AREA OF THE REACTOR
C              CGS UNITS - CM**3
C   Q      - THE HEAT LOSS OF THE REACTOR
C              CGS UNITS - ERGS/SEC
C   EQUIV  - FUEL/OXIDIZER EQUIVALENCE RATIO.
C   NT     - POINTER TO THE TEMPERATURE IN THE SOLUTION AND
C            RESIDUAL VECTORS.  S(NT) = TEMPERATURE
C   NYS    - POINTER TO THE SOLUTION ARRAY.  S(NYS+K) IS THE MASS
C            FRACTION OF THE KTH SPECIES.
C   NY     - POINTER TO THE SOLUTION ARRAY.  S(NY) IS THE MASS
C            FRACTION OF THE 1ST SPECIES.  USED IN SUBROUTINE CALLS.
C
C INPUT TO TWOPNT-
C
C   PMAX   - MAX NUMBER OF GRID POINTS = 1
C   POINTS - CURRENT NUMBER OF GRID POINTS = 1
C   TSTEP0 - INITIAL TIME STEP FOR TWOPNT
C
C WORK AND SCRATCH SPACE-
C   RCKWRK - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   BUFFER - WORK SPACE INTO WHICH TWOPNT WRITES DATA AND THROUGH WHICH
C            DATA IS RETURNED TO TWOPNT THROUGH THE REVERSE
C            COMMUNICATION INTERFACE.  SEE TWOPNT DOCUMENTATION.
C              DIMENSION BELOW(*) AT LEAST NATJ.
C   TWPWK  - WORK SPACE FOR TWOPNT.
C              DIMENSION WNEWTN(LENTWP).
C   IPVT   - ARRAY OF PIVOTS USED BY LINPACK LU FACTORIZATION AND SOLVE
C            ROUTINES.
C              DIMENSION IPVT(*) AT LEAST NATJ.
C   SCRTCH - SCRATCH SPACE USED BY THE FUNCTION, JACOBIAN, AND PRINT.
C              DIMENSION SCRTCH(NATJ,9)
C   SN     - DEPENDENT VARIABLE ARRAY AT PREVIOUS TIMESTEP. THE
C            STRUCTURE IS THE SAME AS S.
C              DIMENSION SN(*) AT LEAST NATJ.
C   F      - RESIDUALS OF THE GOVERNING EQUATIONS AS EVALUATED AT
C            S(N).  THE RESIDUAL OF THE KK SPECIES EQUATIONS IS IN
C            F(NYS+K), THE ENERGY EQUATION IN F(NT,J).  IF
C            INERGY=.FALSE. THEN THE ENERGY EQUATION IS REPLACED
C            BY THE GIVEN TEMPERATURE.
C               DIMENSION F(*) AT LEAST NATJ.
C   A      - STORAGE SPACE FOR THE JACOBIAN MATRIX, AND ITS LU FACTORS.
C               DIMENSION A(NATJ,NATJ)
C
C OUTPUT-
C   X      - ARRAY OF MOLE FRACTIONS AT THE SOLUTION.
C              DIMENSION X(*) AT LEAST KK.
C   S      - DEPENDENT VARIABLE ARRAY.  ON OUTPUT S(*) CONTAINS THE
C            SOLUTION.  THE TEMPERATURE IS STORED IN T=S(NT), AND THE
C            MASS FRACTIONS ARE IN Y(K)=S(NYS+K)
C              DIMENSION S(*) AT LEAST NATJ.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  CHANGE BLOCK TO TURN ON SCALING AND ITERATIVE IMPROVEMENT OF
C          LINEAR ALGEGRA PROBLEM
C*****xlinpk double precision
      XLINPK = .TRUE.
C*****END xlinpk double precision
C*****xlinpk single precision
C      XLINPK = .TRUE.
C*****END xlinpk single precision

      DO 20 N = 1, NATJ
         ACTIVE(N) = .FALSE.
   20 CONTINUE
      ERROR = .FALSE.
      KERR = .FALSE.
C
C        SET FLAG TO INDICATE WHETHER YOU HAVE WRITTEN SOLUTION TO
C          LRECOV OR NOT
C
      LWRITE = .FALSE.
C
C        SET NDT AND VDT FOR THE TEMPERATURE FIXED PROBLEM
C
      NDT    = NUMDT
      TSTEP0 = DT1
C
C          DECIDE HOW MANY TIMES TO CALL TWOPNT
C
      ENERGY = LENRGY
      EENRGY = LENRGE
      IF (ENERGY .OR. EENRGY) THEN
         CALLS = 2
      ELSE
         CALLS = 1
      ENDIF
C
C           RETURN TO PRINTING THE INLET CONDITIONS, THEY ARE TURNED OFF
C            AFTER BEING PRINTED ONCE BY PRINT.
C
      LPRTIC = .TRUE.
C
C          TOP OF THE LOOP OVER CALLS TO TWOPNT.
C
      DO 1000 CALL = 1, CALLS
C
         IF (LNOFT .AND. CALLS.EQ.2 .AND. CALL.EQ.1) GOTO 1000
         IF (CALL.EQ.2 .AND. LEVEL(2).GT.0) WRITE (LOUT,'(/A/)')
     1   '   PSRTWO: FINISHED FIXED TEMPERATURE, ADDING ENERGY EQUATION'
C
         IF (ENERGY) LENRGY = (CALL .NE. 1)
         IF (EENRGY) LENRGE = (CALL .NE. 1)
C
         IF (CALL .GT. 1) THEN
            NDT    = NUMDT2
            TSTEP0 = DT2
            IF (.NOT.LSTEDY) THEN
               IF (ENDTIM .GT. 0.0) THEN
                  NINIT = INT(ENDTIM/TSTEP0)
               ENDIF
               ENDTIM = NINIT*TSTEP0
               IRETIR = MAX(IRETIR,NINIT)
            ENDIF
         ENDIF
C
         SIGNAL = ' '
C
         IPASSS = 1
C
C                   SECTION TO PRINT OUT TWOPNT PARAMETERS
C
         IF (LEVEL(1) .EQ. 0) GO TO 500
C
         WRITE (LOUT, 6000)
         WRITE (LOUT, 6001) LEVEL(1), LEVEL(2)
         IF (LEVEL(1) .EQ. 1) THEN
           WRITE (LOUT, 6002)
         ELSE
           WRITE (LOUT, 6003) LEVEL(1)-1
         ENDIF
         WRITE (LOUT, 6004)
         WRITE (LOUT, 6005)'TEMPERATURE', ABOVE(NT), BELOW(NT)
         IF (KEL.NE.0) WRITE (LOUT,6005)'ELECTRON TEMP', ABOVE(NTE),
     1                 BELOW(NTE)
         DO 450 K = 1, KKGAS
            WRITE (LOUT, 6005)KSYM(K,1), ABOVE(NYS+K), BELOW(NYS+K)
450      CONTINUE
         KSTOT = 0
         DO 455 IM = 1, NMAT
C----Check for more than just a gas phase; or KFIRST(2,IM)=0
            IF (NPHASE(IM) .GT. 1) THEN
               DO 453 KM = KFIRST(2,IM), KLAST(NPHASE(IM),IM)
                  K = KM + KSTOT
                  WRITE (LOUT, 6005) KSYM(KM,IM), ABOVE(NYS+K),
     1                               BELOW(NYS+K)                
 453           CONTINUE
               KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
            ENDIF
C-----
455      CONTINUE
         IDTOT = 0
         DO 465 IM = 1, NMAT
            IF (NSPHCH(IM) .GT. 0) THEN
               DO 460 IDM = 1, NSPHCH(IM)
                  I = IDM + IDTOT
                  WRITE (LOUT, 6017) PSYM(MAPPH(IDM,IM),IM),
     1                               ABOVE(NSDENS+I), BELOW(NSDENS+I)
460            CONTINUE
            ENDIF
            IDTOT = IDTOT + NSPHCH(IM)
465      CONTINUE
         WRITE (LOUT,'(//)')
         WRITE (LOUT, 6006) NDT
         WRITE (LOUT, 6007) ATOL
         WRITE (LOUT, 6008) RTOL
         WRITE (LOUT, 6009) ATIM
         WRITE (LOUT, 6010) RTIM
         WRITE (LOUT, 6011) IRETIR
         WRITE (LOUT, 6012) DFAC
         WRITE (LOUT, 6013) UFAC
         WRITE (LOUT, 6014) DTMIN
         WRITE (LOUT, 6015) TSTEP0
         WRITE (LOUT, 6018) ABSOL, RELAT
         WRITE (LOUT, 6016)

500      CONTINUE

         CALL TWSETL (ERROR, LOUT, 'STEADY', LSTEDY)
         KERR = KERR.OR.ERROR

         CALL TWSETL (ERROR, LOUT, 'ADAPT', .FALSE.)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR

         CALL TWSETI (ERROR, LOUT, 'STEPS1', NDT)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR

         CALL TWSETR (ERROR, LOUT, 'SSABS', ATOL)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR
         
         CALL TWSETR (ERROR, LOUT, 'SSREL', RTOL)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR
         
         CALL TWSETR (ERROR, LOUT, 'TDABS', ATIM)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR
         
         CALL TWSETR (ERROR, LOUT, 'TDREL', RTIM)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR
         
         CALL TWSETI (ERROR, LOUT, 'STEPS2', IRETIR)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR

         CALL TWSETR (ERROR, LOUT, 'TDEC', DFAC)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR
         
         CALL TWSETR (ERROR, LOUT, 'TINC', UFAC)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR
         
         CALL TWSETR (ERROR, LOUT, 'TMAX', DTMAX)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR
         
         CALL TWSETR (ERROR, LOUT, 'TMIN', DTMIN)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR
         
         CALL TWSETR (ERROR, LOUT, 'STRID0', TSTEP0)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR

         CALL TWSETI (ERROR, LOUT, 'SSAGE', NJAC)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR

         CALL TWSETI (ERROR, LOUT, 'TDAGE', ITJAC)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR

         CALL TWSETI (ERROR, LOUT, 'STEPS0', NINIT)
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR

         CALL TWSETI (ERROR, LOUT, 'LEVELD', LEVEL(2))
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR

         CALL TWSETI (ERROR, LOUT, 'LEVELM', LEVEL(1))
C         IF (ERROR) STOP
         KERR = KERR.OR.ERROR
C
         IF (KERR) RETURN
C
C*****precision > double
      VERSIO = 'DOUBLE PRECISION VERSION 3.18'
C*****END precision > double
C*****precision > single
C      VERSIO = 'SINGLE PRECISION VERSION 3.18'
C*****END precision > single


C     SUBROUTINE TWOPNT
C    +  (ERROR, TEXT, VERSIO,
C    +   ABOVE, ACTIVE, BELOW, BUFFER, COMPS, CONDIT, GROUPA, GROUPB,
C    +   ISIZE, IWORK, MARK, NAME, NAMES, PMAX, POINTS, REPORT, RSIZE,
C    +   RWORK, SIGNAL, STRIDE, TIME, U, X)

      CALL TWOPNT
     +  (ERROR, LOUT, VERSIO,
     +   ABOVE, ACTIVE, BELOW, BUFFER, NATJ, CONDIT, 0, 0,
     +   3, ITWPWK, MARK, ' ', 1, PMAX, POINTS, REPORT, LENTWP, TWPWK,
     +   SIGNAL, DT, LTIME, S, X)
C      IF (ERROR) STOP 
      IF (ERROR) THEN
         KERR = .TRUE.
         RETURN
      ENDIF

         IF (SIGNAL .EQ. 'RESIDUAL') THEN

               CALL PSRFUN (KK, KKGAS, KKSURF, KKBULK, NSPHCH, NPTS,
     1                      NATJ, NPHASE, NNSURF, NFSURF, NLSURF, 
     2                      NNBULK, NFBULK, NLBULK, KKPHAS, KFIRST, 
     3                      KLAST, KOCC, ETCH, LTIME, LENRGY, LFLRT, 
     4                      LTDIFF, DT, HIN, YIN, TIN, T, TSURF, P, TAU, 
     5                      FLRT, V, AREA, Q, HTRN, TAMBNT, LHTRN,
     6                      ICKWRK, RCKWRK, ISKWRK, RSKWRK, MAPPH, WT,
     6                      SCRTCH(1,1), NIK, SN, BUFFER,
     7                      F, SCRTCH(1,3), SDOT, ACT, SDEN, SDEN0,
     8                      SITDOT, NT, NTE, NYS, NY, NSDEN, NSDENS,
     9                      KEL, KKION, KCHG, KION, QEK, POWR, LENRGE,
     *                      TEIN, TE, ROP, ROPS, IEIMP, EIONSH, II,
     1                      NQLSE, QLSE, TQLSE, GFAC, SFAC, TIONP,
     2                      LCONFN, IEXC, REXC, SDOTI, SITDTI, IISUR,
     3                      KKTOT, NMAT, IISUR0, KKMAX, NPHMAX, IMISK,
     4                      IMRSK, IMCSK, AFRAC, ITDEP, BHMXI, FLXION,
     5                      FLXCOR, SDOTT, SCRTCH(1,8),
     6                      ESHTH, LELSH, LRFSH, LRFCUR, RFFREQ,
     7                      RFAPAR, RFDPAR, LWALHB, VISC, THCOND,
     8                      NBHM, LOUTSH, BOLTZ, AVOG, EV2K, BIASPW,HIM)
C
C
C*****precision > double
               CALL DCOPY (NATJ, F, 1, BUFFER, 1)
C*****END precision > double
C
C*****precision > single
C               CALL SCOPY (NATJ, F, 1, BUFFER, 1)
C*****END precision > single
               GO TO 500
C
            ELSEIF (SIGNAL .EQ. 'PREPARE') THEN

               CALL PSRJAC (LOUT, KK, KKGAS, KKSURF, KKBULK, NSPHCH,
     1                      NPTS,NATJ, NPHASE, NNSURF, NFSURF, NLSURF,
     2                      NNBULK, NFBULK, NLBULK, KKPHAS, KFIRST,
     3                      KLAST, KOCC, ABOVE, BELOW, ETCH, LTIME,
     4                      LENRGY, LFLRT, LTDIFF, ABSOL, RELAT, DT,
     5                      HIN, YIN, TIN, T, TSURF, P, TAU, FLRT, V,
     6                      AREA, Q, HTRN, TAMBNT, LHTRN, ICKWRK,
     7                      RCKWRK, ISKWRK, RSKWRK, MAPPH, WT, NIK,
     8                      SCRTCH, SN, S, F, FN, A, ACT, SDEN, SDEN0,
     9                      SITDOT, NT, NTE, NYS, NY, NSDEN, NSDENS,  
     *                      KEL, KKION, KCHG, KION, QEK, POWR, LENRGE,  
     1                      TEIN, TE, ROP, ROPS, IEIMP, EIONSH, II,  
     2                      NQLSE, QLSE, TQLSE, GFAC, SFAC, TIONP, 
     3                      LCONFN, IEXC, REXC, SDOTI, SITDTI, IISUR, 
     4                      KKTOT, NMAT, IISUR0, KKMAX, NPHMAX, IMISK, 
     5                      IMRSK, IMCSK, AFRAC, ITDEP, BHMXI, FLXION, 
     6                      FLXCOR, SDOTT, ESHTH, LELSH, LRFSH, 
     7                      LRFCUR, RFFREQ, RFAPAR, RFDPAR, LWALHB, 
     8                      VISC, THCOND, SDOT, NBHM, LOUTSH, BOLTZ,
     9                      AVOG, EV2K, BIASPW, HIM)
C
               IF (XLINPK) THEN
C*****xlinpk double precision
                  CALL XDGECO (A, NATJ, AA, NATJ, NATJ, IPVT,
     1                         RCOND, DR , DC, FN, ANORM)
C*****END xlinpk double precision
C*****xlinpk single precision
C                  CALL XSGECO (A, NATJ, AA, NATJ, NATJ, IPVT,
C     1                         RCOND, DR , DC, FN, ANORM)
C*****END xlinpk single precision
               ELSE
C*****precision > double
                  CALL DGECO (A, NATJ, NATJ, IPVT, RCOND, FN)
C*****END precision > double
C*****precision > single
C                  CALL SGECO (A, NATJ, NATJ, IPVT, RCOND, FN)
C*****END precision > single
               ENDIF
C
               IF (RCOND .LE. 0.0E0) THEN
                  WRITE (LOUT, *)
     1            ' FATAL ERROR, SINGULAR JACOBIAN: start diagnostics '
                  FINDJ = .FALSE.
                  CALL PSFIXJ (LOUT, KK, KKGAS, KKSURF, KKBULK, NSPHCH, 
     1                         NPTS, NATJ, NPHASE, KSYM, PSYM, XLINPK, 
     2                         NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2                         NLBULK, KKPHAS, KFIRST, KLAST, KOCC,
     3                         ABOVE, BELOW, ETCH, LTIME, LENRGY, LFLRT,
     5                         LTDIFF, ABSOL, RELAT, RCOND, DT, HIN, 
     6                         YIN, TIN, T, TSURF, P, TAU, FLRT, V, 
     7                         AREA, Q, HTRN, TAMBNT, LHTRN, ICKWRK, 
     8                         RCKWRK, ISKWRK, RSKWRK, MAPPH, WT, 
     9                         NIK,SCRTCH, SN, S, F, FN, A, AA, DR, DC, 
     *                         IPVT, ACT, SDEN, SDEN0, SITDOT, NT, NTE, 
     1                         NYS, NY, NSDEN, NSDENS, FINDJ, KEL, 
     2                         KKION, KCHG, KION, QEK, POWR, LENRGE,
     3                         TEIN, TE, ROP, IEIMP, EIONSH, II, 
     4                         NQLSE, QLSE, TQLSE, GFAC, SFAC, TIONP,
     5                         LCONFN, IEXC, REXC, SDOTI, SITDTI, KKTOT, 
     6                         NMAT, IISUR0, KKMAX, NPHMAX, IMISK,
     7                         IMRSK, IMCSK, IISUR, AFRAC, ITDEP,
     8                         BHMXI, FLXION, FLXCOR, SDOTT, ESHTH,
     9                         LELSH, LRFSH, LRFCUR, RFFREQ, RFAPAR,
     *                         RFDPAR, LWALHB, VISC, THCOND, SDOT, ROPS,
     1                         NBHM, LOUTSH, BOLTZ, AVOG, EV2K, BIASPW,
     2                         HIM)
                  IF (.NOT. FINDJ) THEN
                     WRITE (LOUT, *)
     1               'A NONSINGULAR JACOBIAN COULD NOT BE FOUND'
C                     STOP
                     KERR = .TRUE.
                     RETURN
                  ENDIF
               ENDIF
               CONDIT = 1.0 / RCOND
               GO TO 500

            ELSEIF (SIGNAL .EQ. 'SOLVE') THEN
C
               IF (XLINPK) THEN
C*****xlinpk double precision
                  CALL XDGESL
C*****END xlinpk double precision
C*****xlinpk single precision
C                  CALL XSGESL
C*****END xlinpk single precision
     1            (A, NATJ, AA, NATJ, NATJ, IPVT, DR, DC,
     2            BUFFER, RCOND, FN, RELERR, INFO, SCRTCH, ANORM)
               ELSE
C*****precision > double
                  CALL DGESL (A, NATJ, NATJ, IPVT, BUFFER, 0)
C*****END precision > double
C*****precision > single
C                  CALL SGESL (A, NATJ, NATJ, IPVT, BUFFER, 0)
C*****END precision > single
               ENDIF
               GO TO 500
C
            ELSEIF (SIGNAL .EQ. 'RETAIN') THEN
C
C*****precision > double
               CALL DCOPY (NATJ, BUFFER, 1, SN, 1)
C*****END precision > double
C
C*****precision > single
C               CALL SCOPY (NATJ, BUFFER, 1, SN, 1)
C*****END precision > single
C
               IF (SPOS .GE. 0.0) THEN
                  KSTOT = 0
                  DO 710 IM = 1, NMAT
                     IF (KKSURF(IM) .GT. 0) THEN
                        DO 690 KM = KFIRST(NFSURF(IM),IM),
     1                                      KLAST(NLSURF(IM),IM)
                           K = KM + KSTOT
                           SN(NYS+K) = MAX (SN(NYS+K), SPOS)
 690                    CONTINUE
                     ENDIF
                     IF (KKBULK(IM) .GT. 0) THEN
                        DO 700 KM = KFIRST(NFBULK(IM),IM), 
     1                                      KLAST(NLBULK(IM),IM)
                           K = KM + KSTOT
                           SN(NYS+K) = MAX (SN(NYS+K), SPOS)
 700                    CONTINUE
                     ENDIF
                     KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
710               CONTINUE
                  DO 720 K = 1, KKGAS
                     IF (KCHG(K).NE.0) THEN
                     ELSE
                        SN(NYS+K) = MAX (SN(NYS+K), SPOS)
                     ENDIF
 720              CONTINUE
               ENDIF
               GO TO 500
            ELSEIF (SIGNAL .EQ. 'SHOW') THEN
               CALL PSPRNT (LOUT, NSPHCH, KK, KKGAS, KKSURF, KKBULK,
     1                      NPTS, NATJ, NPHASE, NNSURF, NFSURF, NLSURF,
     2                      NNBULK, NFBULK, NLBULK, KKPHAS, KFIRST,
     3                      KLAST, KOCC, KSYM, PSYM, ETCH, LENRGY,
     4                      LFLRT, LTDIFF, LEQUIV, LPRTIC, HIN, YIN,
     5                      EQUIV, PA, P, TAU, FLRT, V, AREA, Q, HTRN,
     5                      TAMBNT, LHTRN, TIN, XIN, T, TSURF, X, 
     6                      NIK, SCRTCH, SN, S, F, ICKWRK, RCKWRK,
     7                      ISKWRK, RSKWRK, MAPPH, WT, ACT, SDEN, SDEN0,
     8                      SDOT, DEN, SITDOT,  GRATE, NT, NTE, NYS, NY,
     9                      NSDEN, NSDENS, KEL, KKION, KCHG, KION, QEK,
     *                      POWR, LENRGE, TEIN, TE, ROP, ROPS, IEIMP,
     1                      EIONSH, II, NQLSE, QLSE, TQLSE, GFAC, SFAC,
     2                      TIONP, LCONFN, IEXC, REXC, SDOTI, SITDTI,
     3                      IISUR, NMAT, IISUR0, KKMAX, NPHMAX, KKTOT, 
     4                      IMISK, IMRSK, IMCSK, AFRAC, MSYM, ITDEP,
     5                      BHMXI, FLXION, FLXCOR, SDOTT, ESHTH,
     6                      LELSH, LRFSH, LRFCUR, RFFREQ, RFAPAR,
     7                      RFDPAR, LWALHB, VISC, THCOND, LSTEDY,
     8                      ENDTIM, OLDTIM, LTERRN, DEPRAT, SPUTTR,
     9                      NBHM, LOUTSH, BOLTZ, AVOG, EV2K, BIASPW, 
     *                      HIM)

               GO TO 500
C
            ELSEIF (SIGNAL .EQ. 'SAVE') THEN
               IF (LWRITE) THEN
                  BACKSPACE(LRECOV, ERR=2000)
                  BACKSPACE(LRECOV, ERR=2000)
                  BACKSPACE(LRECOV, ERR=2000)
                  BACKSPACE(LRECOV, ERR=2000)
                  IF (KEL.NE.0) BACKSPACE(LRECOV, ERR=2000)
                  BACKSPACE(LRECOV, ERR=2000)
                  BACKSPACE(LRECOV, ERR=2000)
                  BACKSPACE(LRECOV, ERR=2000)
               ELSE
                  LWRITE = .TRUE.
               ENDIF
               WRITE (LRECOV) ISOLUT
               WRITE (LRECOV) NATJ, IPSR, NUMPSR, NMAT,
     1                        (NSPHCH(IM),IISUR(IM),IM=1,NMAT)
               WRITE (LRECOV) EQUIV, P, TAU, FLRT, V, AREA, Q, 
     1                        (TSURF(IM),IM=1,NMAT),
     2                        HTRN, TAMBNT, LHTRN,
     3                        (AFRAC(IM),IM=1,NMAT)
               WRITE (LRECOV) TIN, (XIN(M), M=1,KKGAS)
               IF (NQLSE.NE.0) CALL PLTEMP (NQLSE, S(NTE), TQLSE, QLSE,
     1                                      QLEX)
               IF (KEL.NE.0) WRITE (LRECOV) TEIN, POWR, TIONP, BHMXI,
     1                   (EIONSH(IM),ESHTH(IM),BIASPW(IM),IM=1,NMAT),
     2                   QLEX
               IF (NSPTOT .GT. 0) THEN
                  WRITE (LRECOV) (BUFFER(N), N=1,NATJ),
     1                  ((MAPPH(I,IM),I=1,MAX(1,NSPHCH(IM))),IM=1,NMAT)
               ELSE
                  WRITE (LRECOV) (BUFFER(N), N=1,NATJ)
               ENDIF
C
               WRITE (LRECOV) 
     1          (((SDOTI(K,I,IM),K=1,KK(IM)),I=1,MAX(1,IISUR(IM))),
     2          IM=1,NMAT), 
     2          ((SITDTI(N,IM), N=1, NPHASE(IM)),IM=1,NMAT)
C
               SOLTIM = OLDTIM + ENDTIM
               WRITE (LRECOV) SOLTIM
C
               GO TO 500
C
         ELSE IF (SIGNAL .NE. ' ') THEN
            GO TO 500
         ENDIF

         IF (REPORT .EQ. ' ') THEN
         ELSE IF (REPORT .EQ. 'NONE FOUND') THEN
            WRITE (LOUT, *)
     1      'PSR FATAL ERROR:  Twopnt reports no success ',
     2      'in solving the problem.'
            WRITE (LOUT, *)
     1      '                  PSR gives up after writing ',
     2      'the restart file!'
C            STOP
            KERR = .TRUE.
            RETURN
         ELSE
C            STOP
             RETURN
         END IF
C
1000  CONTINUE
C
      RETURN
2000  CONTINUE
      WRITE (LOUT, *)'PSR ERROR: Trouble backspacing recover file'
C      STOP
      KERR = .TRUE.
      RETURN
C----------------------------------------------------------------
6000  FORMAT (//120('=')/40X,'TWOPNT PARAMETERS:'/)
6001  FORMAT (10X,'PRINTING: Level(1) = ',I2,'  Level(2) = ',I2)
6002  FORMAT (20X,'Only TWOPNT will write messages')
6003  FORMAT (20X,'TWOPNT subroutines ',I2,' below entrance will'
     1          ,' all write to output')
6004  FORMAT (/20X,'SOLUTION BOUNDS FOR VARIABLES:'//
     1        ' Variable      Upper_bound     Lower_bound')
6005  FORMAT (1X, A11, 3X, 1PG11.3, 3X,1PG11.3)
6017  FORMAT (1X,'Surface site density for phase,', A11,
     1        3X, 1PG11.3, 3X,1PG11.3)
6006  FORMAT (10X,'Maximum number of time steps to try if ',
     1         'newton''s method fails = ',i3)
6007  FORMAT (10X,'Absolute bound for the steady-state problem = ',
     1        1PG11.3)
6008  FORMAT (10X,'Relative bound for the steady-state problem = ',
     1        1PG11.3)
6009  FORMAT (10X,'Absolute bound for the time-dependent prob = ',
     1         1PG11.3)
6010  FORMAT (10X,'Relative bound for the time-dependent prob = ',
     1         1PG11.3)
6011  FORMAT (10X,'Number of steps to be taken before increasing ',
     1           'the time step = ',I3)
6012  FORMAT (10X,'Factor by which to decrease time-step when ',
     1        'necessary = ',1PG11.3)
6013  FORMAT (10X,'Factor by which to increase time-step when ',
     1        'necessary = ',1PG11.3)
6014  FORMAT (10X,'Minimum time step allowed = ',1PG11.3)
6015  FORMAT (10X,'Initial time step = ',1PG11.3)
6018  FORMAT (10X,'Absolute delta for numerical differencing = ',
     1           1PG11.3/
     1        10X,'Relative delta for numerical differencing = ',
     1           1PG11.3)
6016  FORMAT (/120('=')//)


      END
C
C------------------------------------------------------------
C
      SUBROUTINE PSRJAC (LOUT, KK, KKGAS, KKSURF, KKBULK, NSPHCH, NPTS,
     1                   NATJ, NPHASE, NNSURF, NFSURF, NLSURF, NNBULK,
     2                   NFBULK, NLBULK, KKPHAS, KFIRST, KLAST, KOCC,
     3                   ABOVE, BELOW, ETCH, LTIME, LENRGY, LFLRT,
     4                   LTDIFF, ABSOL, RELAT, DT, HIN, YIN, TIN, T,
     5                   TSURF, P, TAU, FLRT, V, AREA, Q, HTRN,
     6                   TAMBNT, LHTRN, ICKWRK, RCKWRK, ISKWRK,
     7                   RSKWRK, MAPPH, WT, NIK, SCRTCH, SN, S, F, FN,
     8                   A, ACT, SDEN, SDEN0, SITDOT, NT, NTE, NYS, NY,
     9                   NSDEN, NSDENS, KEL, KKION, KCHG, KION, QEK,
     *                   POWR, LENRGE, TEIN, TE, ROP, ROPS, IEIMP,
     1                   EIONSH, II, NQLSE, QLSE, TQLSE, GFAC, SFAC,
     2                   TIONP, LCONFN, IEXC, REXC, SDOTI, SITDTI,
     3                   IISUR, KKTOT, NMAT, IISUR0, KKMAX, NPHMAX,
     4                   IMISK, IMRSK, IMCSK, AFRAC, ITDEP, BHMXI,
     5                   FLXION, FLXCOR, SDOTT, ESHTH, LELSH, LRFSH,
     6                   LRFCUR, RFFREQ, RFAPAR, RFDPAR, LWALHB, VISC,
     7                   THCOND, SDOT, NBHM, LOUTSH, BOLTZ, AVOG, EV2K,
     8                   BIASPW, HIM)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C 
      DIMENSION KKPHAS(NPHMAX,NMAT), KFIRST(NPHMAX,NMAT), 
     1          KLAST(NPHMAX,NMAT), KOCC(KKMAX,NMAT), ICKWRK(*),
     1          ISKWRK(*), MAPPH(NPHMAX,NMAT), S(NATJ), F(NATJ),
     2          FN(NATJ), A(NATJ,NATJ), KKSURF(NMAT), KKBULK(NMAT),
     3          SN(NATJ), ABOVE(NATJ), BELOW(NATJ), SCRTCH(NIK,*), 
     4          KK(NMAT), WT(KKMAX,NMAT), HIN(KKMAX), YIN(KKMAX),
     5          RCKWRK(*), RSKWRK(*), ACT(KKMAX,NMAT), 
     5          SDEN(NPHMAX,NMAT), SDEN0(NPHMAX,NMAT),
     6          SITDOT(NPHMAX,NMAT), KCHG(KKMAX), KION(*), QEK(KKGAS), 
     7          ROP(*), IEIMP(*), QLSE(NPTS), TQLSE(NPTS), IEXC(*),  
     8          REXC(*), SDOTI(KKMAX,IISUR0,NMAT), SITDTI(NPHMAX,NMAT), 
     8          ROPS(IISUR0,NMAT), NPHASE(NMAT), NNSURF(NMAT),
     9          NFSURF(NMAT), NLSURF(NMAT), NNBULK(NMAT), NFBULK(NMAT),
     *          NLBULK(NMAT), IISUR(NMAT), NSPHCH(NMAT), IMISK(NMAT),
     1          IMRSK(NMAT), IMCSK(NMAT), AFRAC(NMAT), ITDEP(*), 
     1          EIONSH(NMAT), FLXION(*), FLXCOR(*), SDOTT(KKTOT),
     3          ESHTH(NMAT), SDOT(KKTOT), TSURF(NMAT), NBHM(NMAT),
     4          BIASPW(NMAT), HIM(KKMAX,NMAT)
      LOGICAL LENRGY, LFLRT, LTDIFF(NMAT), LTIME, ETCH(NPHMAX,NMAT), 
     1        LHTRN, LENRGE, LCONFN, LELSH(NMAT), LRFSH(NMAT),  
     2        LRFCUR, LWALHB(NMAT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C            ZERO THE MATRIX STORAGE SPACE.
C
C*****precision > double
      CALL DCOPY (NATJ * NATJ, 0.0D0, 0, A, 1)
C*****END precision > double
C
C*****precision > single
C      CALL SCOPY (NATJ * NATJ, 0.0, 0, A, 1)
C*****END precision > single
C
C            CALL THE FUNCTION AT S AND STORE IN FN.
C
      CALL PSRFUN (KK, KKGAS, KKSURF, KKBULK, NSPHCH, NPTS,NATJ, NPHASE,
     1             NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2             KKPHAS, KFIRST, KLAST, KOCC, ETCH, LTIME, LENRGY,
     3             LFLRT, LTDIFF, DT, HIN, YIN, TIN, T, TSURF, P, TAU,
     4             FLRT, V, AREA, Q, HTRN, TAMBNT, LHTRN, ICKWRK,
     4             RCKWRK, ISKWRK, RSKWRK, MAPPH, WT, SCRTCH(1,1), NIK,
     6             SN, S, FN, SCRTCH(1,3), SDOT, ACT, SDEN,
     7             SDEN0, SITDOT, NT, NTE, NYS, NY, NSDEN, NSDENS, KEL,
     8             KKION, KCHG, KION, QEK, POWR, LENRGE, TEIN, TE, ROP,
     9             ROPS, IEIMP, EIONSH, II, NQLSE, QLSE, TQLSE, GFAC,
     *             SFAC, TIONP, LCONFN, IEXC, REXC, SDOTI, SITDTI,
     1             IISUR, KKTOT, NMAT, IISUR0, KKMAX, NPHMAX, IMISK,
     2             IMRSK, IMCSK, AFRAC, ITDEP, BHMXI, FLXION, FLXCOR,
     3             SDOTT, SCRTCH(1,5), ESHTH, LELSH,
     4             LRFSH, LRFCUR, RFFREQ, RFAPAR, RFDPAR, LWALHB,
     5             VISC, THCOND, NBHM, LOUTSH, BOLTZ, AVOG, EV2K,BIASPW,
     6             HIM)
C
C        TOP OF THE LOOPS OVER THE RESIDUE CLASSES AND
C                                              SOLUTION COMPONENTS.
C
      DO 200 M = 1, NATJ
C
C            FOR A GIVEN RESIDUE CLASS AND A GIVEN SOLUTION COMPONENT,
C            PERTRB THE S VECTOR AT POINTS IN THE SAME RESIDUE CLASS.
C
         SAVE = S(M)
         PERTRB = ABS(S(M)) * RELAT + ABSOL
         IF (S(M)+PERTRB .LT. ABOVE(M)) THEN
           S(M) = S(M) + PERTRB
         ELSEIF (S(M)-PERTRB .GT. BELOW(M)) THEN
           PERTRB = -1. * PERTRB
           S(M) = S(M) + PERTRB
         ELSE
           S(M) = S(M) + PERTRB
           WRITE (LOUT, *)'WARNING PSRJAC: variable ',M,'value = ',S(M)
           WRITE (LOUT, *)'                perturbation = ',PERTRB
           WRITE (LOUT, *)'                bounds = ',BELOW(M),
     1                '  ',ABOVE(M)
         ENDIF
C
C             CALL THE FUNCTION AT THE PERTURBED S AND STORE
C
         CALL PSRFUN (KK, KKGAS, KKSURF, KKBULK, NSPHCH, NPTS, NATJ,
     1                NPHASE, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2                NLBULK, KKPHAS, KFIRST, KLAST, KOCC, ETCH, LTIME,
     3                LENRGY, LFLRT, LTDIFF, DT, HIN, YIN, TIN, T,
     4                TSURF, P, TAU, FLRT, V, AREA, Q, HTRN, TAMBNT,
     4                LHTRN, ICKWRK, RCKWRK, ISKWRK, RSKWRK, MAPPH, WT,
     6                SCRTCH(1,1), NIK, SN, S, F,
     7                SCRTCH(1,3), SDOT, ACT, SDEN, SDEN0, SITDOT, NT,
     8                NTE, NYS, NY, NSDEN, NSDENS, KEL, KKION, KCHG,
     9                KION, QEK, POWR, LENRGE, TEIN, TE, ROP, ROPS,
     *                IEIMP, EIONSH, II, NQLSE, QLSE, TQLSE, GFAC, SFAC,
     1                TIONP, LCONFN, IEXC, REXC, SDOTI, SITDTI, IISUR,
     2                KKTOT, NMAT, IISUR0, KKMAX, NPHMAX, IMISK, IMRSK,
     3                IMCSK, AFRAC, ITDEP, BHMXI, FLXION, FLXCOR,
     4                SDOTT, SCRTCH(1,5), ESHTH, LELSH,
     5                LRFSH, LRFCUR, RFFREQ, RFAPAR, RFDPAR, LWALHB,
     6                VISC, THCOND, NBHM, LOUTSH, BOLTZ, AVOG, EV2K,
     7                BIASPW, HIM)
C
C              RESTORE S TO ITS ORIGINAL VALUE.
C
         S(M) = SAVE
C
C              DIFFERENCE TO GET THE COLUMNS OF THE JACOBIAN.
C
         DO 0100 N = 1, NATJ
            A(N, M) = (F(N) - FN(N)) / PERTRB
  100    CONTINUE
C
C          BOTTOM OF THE LOOPS OVER THE RESIDUE CLASSES AND SOLUTION
C          COMPONENTS.
C
  200 CONTINUE
C*****print jacobian
C      PRINT *,'=================================================='
C      PRINT *,'JACOBIAN: '
C      PRINT *
C      DO 210 N = 1,NATJ
C      PRINT *,'ROW # ',N
C      PRINT 276,(A(N,M),M=1,NATJ)
C210   CONTINUE
C      PRINT *,'=================================================='
C276    FORMAT(10(1X,1PG12.4))
C*****END print jacobian
C
      RETURN
      END
C
C--------------------------------------------------------------------
C
      SUBROUTINE PSRFUN (KK, KKGAS, KKSURF, KKBULK, NSPHCH, NPTS, NATJ,
     1                   NPHASE, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2                   NLBULK, KKPHAS, KFIRST, KLAST, KOCC, ETCH,
     3                   LTIME, LENRGY, LFLRT, LTDIFF, DT, HIN, YIN,
     4                   TIN, T, TSURF, P, TAU, FLRT, V, AREA, Q, HTRN,
     5                   TAMBNT, LHTRN, ICKWRK, RCKWRK, ISKWRK,
     6                   RSKWRK, MAPPH, WT, WDOT, NIK, SN, S, F, X,
     7                   SDOT, ACT, SDEN, SDEN0, SITDOT, NT, NTE, NYS,
     8                   NY, NSDEN, NSDENS, KEL, KKION, KCHG, KION,
     9                   QEK, POWR, LENRGE, TEIN, TE, ROP, ROPS, IEIMP,
     *                   EIONSH, II, NQLSE, QLSE, TQLSE, GFAC, SFAC,
     1                   TIONP, LCONFN, IEXC, REXC, SDOTI, SITDTI, 
     2                   IISUR, KKTOT, NMAT, IISUR0, KKMAX, NPHMAX,
     3                   IMISK, IMRSK, IMCSK, AFRAC, ITDEP, BHMXI, 
     4                   FLXION, FLXCOR, SDOTT, CONC, ESHTH,
     5                   LELSH, LRFSH, LRFCUR, RFFREQ, RFAPAR, RFDPAR,
     6                   LWALHB, VISC, THCOND, NBHM, LOUTSH, BOLTZ, 
     7                   AVOG, EV2K, BIASPW, HIM)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DOUBLE PRECISION MASLOS, LOWSDN
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      REAL MASLOS, LOWSDN
C*****END precision > single
C
      DIMENSION KKPHAS(NPHMAX,NMAT), KFIRST(NPHMAX,NMAT), 
     1          KLAST(NPHMAX,NMAT), KOCC(KKMAX,NMAT), ICKWRK(*),
     2          ISKWRK(*), MAPPH(NPHMAX,NMAT), S(NATJ), SN(NATJ),
     2          F(NATJ), WT(KKMAX,NMAT), HIN(KKMAX), YIN(KKMAX),
     3          WDOT(NIK), RCKWRK(*), RSKWRK(*), X(NIK), SDOT(KKTOT),
     4          ACT(KKMAX,NMAT), SDEN(NPHMAX,NMAT), SDEN0(NPHMAX,NMAT),
     5          SITDOT(NPHMAX,NMAT), KCHG(KKMAX), KION(*), TEMP(3), 
     5          QEK(KKGAS), ROP(*), IEIMP(*), QLSE(NPTS), TQLSE(NPTS),
     6          IEXC(*), REXC(*), SDOTI(KKMAX,IISUR0,NMAT),
     7          SITDTI(NPHMAX,NMAT), KKSURF(NMAT), KKBULK(NMAT),
     8          KK(NMAT), NPHASE(NMAT), NNSURF(NMAT), NFSURF(NMAT),
     9          NLSURF(NMAT), NNBULK(NMAT), NFBULK(NMAT), NLBULK(NMAT), 
     *          IISUR(NMAT), NSPHCH(NMAT), IMISK(NMAT), IMRSK(NMAT),
     *          IMCSK(NMAT), AFRAC(NMAT), EIONSH(NMAT), FLXION(*),
     1          FLXCOR(*), SDOTT(KKTOT), CONC(NIK), ESHTH(NMAT),
     2          ITDEP(*), ROPS(IISUR0,NMAT), TSURF(NMAT), NBHM(NMAT),
     3          BIASPW(NMAT), HIM(KKMAX,NMAT)
      LOGICAL LENRGY, LFLRT, LTIME, LTDIFF(NMAT), ETCH(NPHMAX,NMAT), 
     1        LHTRN, LENRGE, LCONFN, LELSH(NMAT), LRFSH(NMAT), 
     1        LRFCUR, LWALHB(NMAT)
      PARAMETER(LOWSDN = 1.0E-10, ZEROE=0.0, ONE=1.0)
      PARAMETER(SMALL=1.E-25)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   DT     - THE TIME STEP THAT IS USED IF LTIME=.TRUE.
C              CGS UNITS - SEC
C   FLRT   - THE MASS FLOW RATE.
C              CGS UNITS - GM/SEC
C   LTIME  - IF LTIME=.TRUE.  A TIME STEP OF DT WILL BE ADDED INTO
C            THE RESIDUAL.
C   HIN    - ARRAY OF SPECIES ENTHALPIES AT INLET TEMPERATURE TIN.
C              CGS UNITS - ERGS/G
C              DIMENSION HIN(*) AT LEAST KK.
C   KK     - NUMBER OF CHEMICAL SPECIES.
C   LENRGY - IF LENRGY=.TRUE. THEN THE ENERGY EQUATION IS TO BE
C            SOLVED, OTHERWISE A FIXED TEMPERATURE IS USED.
C   LFLRT  - IF LFLRT=.TRUE. THEN THE MASS FLOW RATE IS SPECIFIED AND
C              THE RESIDENCE TIME HAS TO BE CALCULATED.
C            IF LFLRT=.FALSE. THEN THE RESIDENCE TIME IS SPECIFIED AND
C              THE MASS FLOW RATE HAS TO BE CALCULATED.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES. NATJ=KK+2.
C   P      - THE PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   Q      - THE HEAT LOSS OF THE REACTOR
C              CGS UNITS - ERGS/SEC
C   S      - DEPENDENT VARIABLE ARRAY. THE TEMPERATURE IS STORED IN
C            T=S(NT), AND THE MASS FRACTIONS ARE IN Y(K)=S(NYS+K)
C              DIMENSION S(*) AT LEAST NATJ.
C   SN     - DEPENDENT VARIABLE ARRAY AT PREVIOUS TIMESTEP.
C   T      - THE TEMPERATURE. FOR FIXED-TEMPERATURE PROBLEMS THIS
C            IS THE USER-SPECIFIED TEMPERATURE.
C              CGS UNITS - K
C   TSURF  - SURFACE TEMPERATURE
C              CGS UNITS - K
C   TAU    - THE NOMINAL RESIDENCE TIME OF THE REACTOR
C              CGS UNITS - SEC
C   TIN    - THE INLET TEMPERATURE.
C              CGS UNITS - K
C   V      - THE VOLUME OF THE REACTOR
C              CGS UNITS - CM**3
C   AREA   - THE AREA OF THE REACTOR
C              CGS UNITS - CM**2
C   WT     - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
C              CGS UNITS - GM/MOLE
C              DIMENSION WT(*) AT LEAST KK.
C   YIN    - ARRAY OF INLET SPECIES MASS FRACTIONS.
C              DIMENSION YIN(*) AT LEAST KK.
C
C WORK AND SCRATCH SPACE-
C   RCKWRK - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   HIM      - ARRAY OF SPECIES ENTHALPIES.
C              CGS UNITS - ERGS/MOLE
C              DIMENSION HIM(*,NMAT) AT LEAST KK.
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   WDOT   - ARRAY OF CHEMICAL PRODUCTION RATES.
C              CGS UNITS - MOLES/(CM**3-SEC)
C              DIMENSION WDOT(*) AT LEAST KK.
C   SDOT   - ARRAY OF CHEMICAL PRODCUTION RATES FROM SURFACE REACTIONS
C              CGS UNITS - MOLES/(CM**2-SEC)
C              DIMENSION SDOT(*) AT LEAST KK.
C   X      - ARRAY OF GAS MOLE FRACTIONS, SURFACE SITE FRACTIONS,
C              AND BULK MOLE FRACTIONS
C              DIMENSION X AT LEAST KK.
C   ACT    - ARRAY OF ACTIVITIES
C              DIMENSION ACT AT LEAST KK.
C
C OUTPUT-
C   F      - RESIDUALS OF THE GOVERNING EQUATIONS AS EVALUATED AT
C            S(N).  THE RESIDUAL OF THE KK SPECIES EQUATIONS IS IN
C            F(NYS+K), THE ENERGY EQUATION IN F(NT,J).  IF
C            INERGY=.FALSE. THEN THE ENERGY EQUATION IS REPLACED
C            BY THE GIVEN TEMPERATURE.
C               DIMENSION F(*) AT LEAST NATJ.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C             FORM THE CHEMICAL RATE TERMS
C
      TEMP(1) = TIN
      TEMP(2) = TEIN
      TEMP(3) = TIN
      CALL CKRHOY (P, TEMP, YIN, ICKWRK, RCKWRK, RHOIN)
      TEMP(1) = S(NT)
      TEMP(2) = S(NTE)
      TEMP(3) = TIONP
      CALL CKWYP (P, TEMP, S(NY), ICKWRK, RCKWRK, WDOT)
      CALL CKRHOY(P, TEMP, S(NY), ICKWRK, RCKWRK, RHO)
      IDTOT = 0
      DO 40 IM = 1, NMAT
C
C            FIND THE SURFACE TEMPERATURE
C
         IF (.NOT. LTDIFF(IM)) TSURF(IM) = S(NT)
         IF (LWALHB(IM)) THEN
C            VISC = 0.00051645
            QOUT = -HTRN*(S(NT)-TAMBNT)
            DEFF = 4*V/AREA
            AXSEC  = 0.25*3.14159*DEFF**2
            LEFF = V/AXSEC
            RENUM = FLRT*DEFF/(AXSEC*VISC)
            PRNUM = 0.75
C            THCOND = 6500.
            REPRLD = (DEFF/LEFF)*RENUM*PRNUM
            XNUSS = 3.66 + (0.0668*REPRLD)/(1.+0.04*REPRLD**0.6667)
            HTRIN = XNUSS*THCOND/DEFF
            TSURF(IM) = S(NT) + QOUT/HTRIN
            IF (TSURF(IM).LE.TAMBNT) TSURF(IM) = TAMBNT
            IF (TSURF(IM).GE.S(NT)) TSURF(IM) = S(NT)
         ENDIF
C
C  calculate the site densities if they change
C
         IF (NNSURF(IM) .GT. 0) THEN
            IF (NSPHCH(IM) .GT. 0) THEN
               DO 30 IDM = 1, NSPHCH(IM)
                  I = IDM + IDTOT
                  SDEN(MAPPH(IDM,IM),IM) = S(NSDENS+I) 
     1                                     * SDEN0(MAPPH(IDM,IM),IM)
30             CONTINUE
            ENDIF
         ENDIF
         IDTOT = IDTOT + NSPHCH(IM)
40    CONTINUE
C
C  FORM X, ACT
C
      CALL CKYTX(S(NY), ICKWRK, RCKWRK, X)
      KSTOT = 0
      NNBTOT = 0
      DO 220 IM = 1, NMAT
         DO 50 K = 1 , KKGAS
            ACT(K,IM) = X(K)
50       CONTINUE
         IF (NNSURF(IM) .GT. 0) THEN
            DO 100 KM = KFIRST(NFSURF(IM),IM), KLAST(NLSURF(IM),IM)
               K = KM + KSTOT
               ACT(KM,IM) = S(NYS+K)
               X(K) = S(NYS+K)
100         CONTINUE
         ENDIF
         IF (NNBULK(IM) .GT. 0) THEN
           DO 200 KM = KFIRST(NFBULK(IM),IM), KLAST(NLBULK(IM),IM)
              K = KM + KSTOT
              X(K) = S(NYS+K)
200        CONTINUE
         ENDIF
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
         NNBTOT = NNBTOT + NNBULK(IM)
220   CONTINUE
C
C  CALL THE USER ROUTINE TO CALCULATION THE BULK PHASE ACTIVITIES
C
      IF (NNBTOT .GT. 0) THEN
         CALL PSRACT (X, KK, KKBULK, NFBULK, NLBULK, NPHASE, KFIRST, 
     1                KLAST, NNBULK, NMAT, NPHMAX, KKSURF, KKMAX, ACT)
      ENDIF
C
C  calculate sdot from surface chemkin, with bohm corrections.
C
      CALL PSSDOT (KKION, KEL, KION, KKGAS, KKTOT, KCHG, NMAT, NPHASE,
     1             IISUR, KFIRST, KLAST, NFSURF, NLSURF, NNSURF, KKSURF,
     2             NFBULK, NLBULK, NNBULK, KKBULK, IISUR0, KKMAX,
     3             NPHMAX, KK, BHMXI, WT, EIONSH, NIK, AFRAC,
     3             RFFREQ, RFAPAR, RFDPAR, LRFSH, LRFCUR, P, ACT,
     4             CONC, TIONP, S(NTE), TSURF, SDEN, ISKWRK, RSKWRK,
     6             IMISK, IMRSK, FLXCOR, FLXION, ESHTH, LELSH, SFAC,
     7             SDOTT, ROPS, SDOTI, SITDTI, SITDOT, SDOT, SHLOSS,
     8             NBHM, LOUTSH, BOLTZ, AVOG, EV2K, BIASPW, AREA, 
     9             HIM(1,1))
C
C  MULTIPLY WDOT BY USER-INPUT FACTORS (SFAC INCLUDED IN PSSDOT)
C
      DO 1200 K = 1, KKGAS
         WDOT(K) = WDOT(K) * GFAC
1200  CONTINUE
C
C  convert between flow rates and residence time
C
      IF (LFLRT) THEN
         TAU = RHO * V / FLRT
      ELSE
         FLRT = RHO * V / TAU
      ENDIF
C
C  calculate net mass loss from gas (AFRAC is included in SDOT here)
C
      MASLOS = ZEROE
      DO 1600 K = 1 , KKGAS
         MASLOS = MASLOS + SDOT(K)*WT(K,1)
1600  CONTINUE
      MASLOS = MASLOS*AREA
C
C  SPECIES CONSERVATION EQUATION; for electron, set species residual to 
C                                 electron energy equation, and v-versa
C
      DUM1 = ONE/TAU
      DUM2 = ONE/RHO
      DUM3 = ONE/(RHO*V)
      DO 1800 K = 1, KKGAS
         IF (KCHG(K).LT.0 .AND. K.NE.KEL .AND. LCONFN) THEN
            F(NYS+K) = YIN(K)*DUM1 + (WT(K,1)*DUM2)*WDOT(K)
     1           + ( WT(K,1)*SDOT(K)*AREA - S(NYS+K)*MASLOS )*DUM3
         ELSEIF (K .EQ. KEL) THEN
            F(NTE) = - (S(NYS+K) - YIN(K))*DUM1 + (WT(K,1)*DUM2)*WDOT(K)
     1           + ( WT(K,1)*SDOT(K)*AREA - S(NYS+K)*MASLOS )*DUM3
         ELSE
            F(NYS+K) = - (S(NYS+K)-YIN(K))*DUM1 + (WT(K,1)*DUM2)*WDOT(K)
     1           + ( WT(K,1)*SDOT(K)*AREA - S(NYS+K)*MASLOS )*DUM3
         ENDIF
1800  CONTINUE
C
C  PLASMA EQUATIONS: calculate collision terms for electron energy
C
      COLLE = 0.0
      COLLI = 0.0
      EPROD = 0.0
      EXCLSS = 0.0
      QLEX = 0.0
      BILOSS = 0.0
C
      IF (KEL .NE. 0) THEN
C
         RU = BOLTZ*AVOG
         CPE = 2.5*RU/WT(KEL,1)
         CVE = 1.5*RU/WT(KEL,1)
C
C           ELASTIC COLLISIONS:
C
         DO 2000 K = 1, KKGAS
            XNK = S(NYS+K)*AVOG*RHO/WT(K,1)
            GE  = SQRT(8.*BOLTZ*S(NTE)*AVOG/(WT(KEL,1)*3.14159))
            VEK = QEK(K)*XNK*GE
            COLLE = COLLE + VEK*WT(KEL,1)/WT(K,1)
2000     CONTINUE
         XNE = S(NYS+KEL)*AVOG*RHO/WT(KEL,1)
         COLLE = COLLE * 3.*V*BOLTZ*XNE*(S(NTE) - S(NT))
C
C           INELASTIC COLLISIONS:
C
         TEMP(1) = S(NT)
         TEMP(2) = S(NTE)
         TEMP(3) = TIONP
         CALL CKQYP (P, TEMP, S(NY), ICKWRK, RCKWRK, ROP)
         CALL CKHML (TEMP, ICKWRK, RCKWRK, HIM(1,1))
         DO 2400 I = 1, II
            IF (IEIMP(I) .GT. 0 .OR. ITDEP(I) .EQ. KEL) THEN
               CALL CKHRX (I, HIM(1,1), ICKWRK, RCKWRK, HRXI)
               COLLI = COLLI + ROP(I) * HRXI * V
               IF (IEXC(I) .GT. 0) THEN
                  ELOSS = REXC(I)*EV2K*RU - HRXI
                  EXCLSS = EXCLSS + ROP(I) * ELOSS * V
               ENDIF
            ENDIF
2400     CONTINUE
         COLLI = COLLI + EXCLSS
         
C
C           LOSS OF POWER COUPLING TO POSITIVE IONS IN SHEATH
C
         SHLOSS = SHLOSS*AVOG*BOLTZ*AREA
C     
C           ELECTRON PRODUCTION (THERMALIZATION) TERM:
C
         EPROD = WDOT(KEL)*V*CPE*WT(KEL,1)*(S(NTE) - S(NT))
C
C           POWER COUPLING TO IONS IN GAS
C
         DO 3000 KI = 1, KKION
            K = KION(KI)
            BILOSS = BILOSS + WDOT(K)*AVOG*V*(TIONP - S(NT))*BOLTZ
3000     CONTINUE
C          
C
C           ELECTRON ENERGY LOSS TO EXCITATION REACTIONS
C           NOT INCLUDED IN THE CHEM MECHANISM BUT INPUT BY USER
C
         IF (NQLSE .GT.0) THEN
            CALL PLTEMP (NQLSE, S(NTE), TQLSE, QLSE, QLEX)
            QLEX = QLEX * EV2K * RU * WDOT(KEL)*V
         ENDIF
C
C  ELECTRON ENERGY EQUATION
C
         IF (LENRGE) THEN
            F(NYS+KEL) = FLRT*YIN(KEL)*CPE*(TEIN-S(NTE))
     2            - (COLLE + COLLI + EPROD + QLEX )
     3            + (POWR - SHLOSS - BILOSS)
         ELSE
           F(NYS+KEL) = S(NTE) - TE
         ENDIF
      ELSE
         F(NTE) = S(NTE) - S(NT)
      ENDIF
C
C  GAS ENERGY EQUATION
C
      IF (LENRGY) THEN
         TEMP(1) = S(NT)
         TEMP(2) = S(NTE)
         TEMP(3) = TIONP
         DO 3180 IM = 1, NMAT
            CALL SKHML (TEMP, ISKWRK(IMISK(IM)), RSKWRK(IMRSK(IM)),
     1                  HIM(1,IM))
 3180    CONTINUE
         CALL CKCPBS (TEMP, S(NY), ICKWRK, RCKWRK, CPB)
         SUM1 = ZEROE
         SUM2 = ZEROE
         SUM3 = ZEROE
         SUM4 = ZEROE
         DO 3200 K = 1, KKGAS
            IF (K.NE.KEL) THEN
               SUM1 = SUM1 + YIN(K) * (HIN(K)-HIM(K,1)/WT(K,1))
               SUM2 = SUM2 + HIM(K,1) * WDOT(K)
            ENDIF
3200     CONTINUE
         DO 3400 K = 1, KKTOT
            IF (K .LE. KKGAS) THEN
               HS = HIM(K,1)
            ELSE
               KSTOT = 0
               DO 3300 IM = 1, NMAT
                  KMBEG = KKGAS + KSTOT
                  KMENDS = KMBEG + KKSURF(IM) 
                  KMENDB = KMBEG + KKSURF(IM) + KKBULK(IM)
                  IF (K .GT. KMBEG .AND. K .LE. KMENDS) THEN
                     KM = K - KSTOT
                     HS = HIM(KM,IM)
                  ELSEIF (K .GT. KMENDS .AND. K .LE. KMENDB) THEN
                     KM = K - KSTOT
                     HS = HIM(KM,IM)*AFRAC(IM)
                  ENDIF
                  KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
3300           CONTINUE
            ENDIF
            IF (KEL.NE.0) THEN
               IF (K.EQ.KEL.AND.SDOT(KEL).GT.0.0) THEN
                  SUM3 = SUM3 + HS*(S(NT)/S(NTE)) * SDOT(K)
               ELSE
                  SUM3 = SUM3 + HS * SDOT(K)
               ENDIF
            ELSE
               SUM3 = SUM3 + HS * SDOT(K)
            ENDIF
3400     CONTINUE
         IF (LHTRN) THEN
           F(NT) = SUM1 / (CPB*TAU) - SUM2 / (RHO*CPB)
     1         - AREA*(SUM3 + HTRN*(S(NT)-TAMBNT)) / (RHO*V*CPB)
         ELSE
           F(NT) = SUM1 / (CPB*TAU) - SUM2 / (RHO*CPB)
     1             - ( SUM3*AREA + Q )/ (RHO*V*CPB)
         ENDIF
C
C Add Plasma corrections to gas energy equation
C
         IF (KEL .NE. 0) THEN
            IF (WDOT(KEL).GT.0.0) THEN
               HEL = HIM(KEL,1)*S(NT)/S(NTE)
            ELSE
               HEL = HIM(KEL,1)
            ENDIF
            F(NT) = F(NT)
     3             - ( WDOT(KEL) * HEL ) /(RHO*CPB)
     6             + (SHLOSS + BILOSS)/(RHO*V*CPB)
     7             + (COLLE + COLLI - EXCLSS + QLEX) /(RHO*V*CPB)
         ENDIF
      ELSE
         F(NT) = S(NT) - T
      ENDIF
C
C SURFACE AND BULK SPECIES EQUATIONS
C
C  LOOP THROUGH MATERIALS 
C
      KSTOT = 0
      IDTOT = 0
      DO 7700 IM = 1, NMAT
         IFLAG = 0
         IF (NNSURF(IM) .GT. 0 .AND. KKSURF(IM) .GT. 0) THEN
            DO 3600 IPHASE = NFSURF(IM), NLSURF(IM)
               IF (SDEN(IPHASE,IM) .LE. ZEROE) THEN
                  X(IPHASE) = SDEN(IPHASE,IM)
                  IFLAG = 1
                  SDEN(IPHASE,IM) = LOWSDN * SDEN0(IPHASE,IM)
               ENDIF
3600        CONTINUE
            IF (IFLAG .EQ. 1) THEN
C
C CASE WHEN SDEN < 0.0: 
C change sden(iphase) to a small value and calculate new sdot and
C sitdot vectors (store them in wdot and rop). 
C Use him(*,1) for ion energies in PSSDOT.
C
               DO 4000 K = 1, KK(IM)
                  WDOT(K) = 0.0
                  SDOTT(K) = 0.0
4000           CONTINUE
               DO 4050 N = 1, NPHASE(IM)
                  ROP(N) = 0.0
4050           CONTINUE
C
               CALL PSSDOT (KKION, KEL, KION, KKGAS, KKTOT, KCHG, NMAT, 
     1                   NPHASE, IISUR, KFIRST, KLAST, NFSURF, NLSURF, 
     2                   NNSURF, KKSURF, NFBULK, NLBULK, NNBULK, KKBULK,
     3                   IISUR0, KKMAX, NPHMAX, KK, BHMXI, WT, EIONSH,
     4                   NIK, AFRAC, RFFREQ, RFAPAR, RFDPAR,
     5                   LRFSH, LRFCUR, P, ACT, CONC, TIONP, S(NTE), 
     6                   TSURF, SDEN, ISKWRK, RSKWRK, IMISK, IMRSK,
     7                   FLXCOR, FLXION, ESHTH, LELSH, SFAC, SDOTT,
     8                   ROPS, SDOTI, SITDTI, ROP, WDOT, SHLOSS, NBHM,
     9                   LOUTSH, BOLTZ, AVOG, EV2K, BIASPW, AREA, 
     *                   HIM(1,1))

               DO 4600 IPHASE = NFSURF(IM), NLSURF(IM)
                  SDEN(IPHASE,IM) = X(IPHASE)
4600           CONTINUE
            ENDIF
C
C SURFACE SPECIES EQUATION:
C
            DO 5600 IPHASE = NFSURF(IM), NLSURF(IM)
               IF (SDEN(IPHASE,IM) .GT. ZEROE) THEN
                  DUM2 = ONE/SDEN(IPHASE,IM)
                  SUM1 = ONE
                  DO 4800 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                     K = KM + KSTOT
                     SUM1 = SUM1 - S(NYS+K)
4800              CONTINUE
                  DO 5000 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                     K = KM + KSTOT
                     F(NYS+K) =   (KOCC(KM,IM)*SDOT(K)
     1                                    - S(NYS+K)*SITDOT(IPHASE,IM))
     2                           * DUM2  + S(NYS+K)*SUM1*DUM1
5000              CONTINUE
               ELSE
                  SDEN(IPHASE,IM) = LOWSDN * SDEN0(IPHASE,IM)
                  DUM2 = ONE/SDEN(IPHASE,IM)
                  SUM1 = ONE
                  DO 5200 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                     K = KM + KSTOT
                     SUM1 = SUM1 - S(NYS+K)
5200              CONTINUE
                  DO 5400 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                     K = KM + KSTOT
                     F(NYS+K) =
     1              (KOCC(KM,IM)*WDOT(KM) - S(NYS+K)*ROP(IPHASE)) * DUM2
     2               + S(NYS+K)*SUM1*DUM1
5400              CONTINUE
                  SDEN(IPHASE,IM) = X(IPHASE)
               ENDIF
5600        CONTINUE
         ENDIF
C
C BULK SPECIES EQUATIONS
C
         IF (NNBULK(IM) .GT. 0 .AND. KKBULK(IM) .GT. 0) THEN
            DO 7400 IPHASE = NFBULK(IM), NLBULK(IM)
               IF (ETCH(IPHASE,IM)) THEN
                  SUM1 = 1.0
                  DO 5800 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                     K = KM + KSTOT
                     SUM1 = SUM1 - S(NYS+K)
5800              CONTINUE
                  DO 6000 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                     K = KM + KSTOT
                     F(NYS+K) = (SN(NYS+K)-S(NYS+K) + S(NYS+K)*SUM1)/TAU
6000              CONTINUE
               ELSE
                  GROWTH = ZEROE
                  DO 6200 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                     K = KM + KSTOT
C*****precision > double
                     GROWTH = GROWTH + DMAX1(ZEROE,SDOT(K))
C*****END precision > double
C*****precision > single
C                     GROWTH = GROWTH + MAX(ZEROE,SDOT(K))
C*****END precision > single
6200              CONTINUE
                  IF (ABS(GROWTH) .GT. ZEROE) THEN
                     THETA =GROWTH*TAU
                     DO 6400  KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                        K = KM + KSTOT
                        IF (SDOT(K) .GE. ZEROE) THEN
                           F(NYS+K) = (SDOT(K) -S(NYS+K)*GROWTH)/THETA
                        ELSE
                           F(NYS+K) = - S(NYS+K)/TAU
                        ENDIF
6400                 CONTINUE
                  ELSE
                     IF (IFLAG .EQ. 1) THEN
                        GROWTH = ZEROE
                        DO 6600 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                           K = KM + KSTOT
C*****precision > double
                           GROWTH = GROWTH 
     1                              + DMAX1(ZEROE,WDOT(KM))
C*****END precision > double
C*****precision > single
C                           GROWTH = GROWTH
C     1                              + MAX(ZEROE,WDOT(KM))
C*****END precision > single
6600                    CONTINUE
                     ENDIF
                     IF (ABS(GROWTH) .GT. ZEROE) THEN
                        THETA =GROWTH*TAU
                        DO 6800 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                           K = KM + KSTOT
                           IF (WDOT(KM) .GE. ZEROE) THEN
                              F(NYS+K) = (WDOT(KM)
     1                                    -S(NYS+K)*GROWTH) /THETA
                           ELSE
                              F(NYS+K) = - S(NYS+K)/TAU
                           ENDIF
6800                    CONTINUE
                     ELSE
C    (give up and revert to alternate equations)
                        IF (KKPHAS(IPHASE,IM) .EQ. 1) THEN
                           F(NYS+KFIRST(IPHASE,IM)+KSTOT) =
     1                     (1.0 - S(NYS+KFIRST(IPHASE,IM)+KSTOT))/TAU
                        ELSE
                           SUM1 = ONE
                           DO 7000 KM = KFIRST(IPHASE,IM),
     1                                               KLAST(IPHASE,IM)
                              K = KM + KSTOT
                              SUM1 = SUM1 - S(NYS+K)
7000                       CONTINUE
                           DO 7200 KM = KFIRST(IPHASE,IM), 
     1                                               KLAST(IPHASE,IM)
                              K = KM + KSTOT
                              F(NYS+K) = (SN(NYS+K)-S(NYS+K)
     1                              + S(NYS+K)*SUM1)/TAU
7200                       CONTINUE
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
7400        CONTINUE
         ENDIF
C
C  SURFACE SITE DENSITY EQUATIONS, IF NECESSARY
C
         IF (NSPHCH(IM) .GT. 0) THEN
            DO 7600 IDM = 1, NSPHCH(IM)
               I = IDM + IDTOT
               F(NSDENS+I) = SITDOT(MAPPH(IDM,IM),IM)
7600        CONTINUE
         ENDIF
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
         IDTOT = IDTOT + NSPHCH(IM)
7700  CONTINUE
C
C  ADD THE TIME STEP, IF NEEDED
C
      IF (LTIME) THEN
C
         DUM1 = ONE/DT
         DO 7800 K = 1, KKTOT
            DYDT = (S(NYS+K) - SN(NYS+K)) * DUM1
            IF (K.EQ.KEL) THEN
               F(NTE) = F(NTE) - DYDT
            ELSE
               F(NYS+K) = F(NYS+K) - DYDT
            ENDIF
7800     CONTINUE
C
         IF (LENRGY) THEN
            DTDT = (S(NT) - SN(NT)) * DUM1
            F(NT) = F(NT) - DTDT 
         ENDIF
C
         IF (LENRGE .AND. KEL .NE. 0) THEN
            DTEDT = (S(NTE) - SN(NTE)) * DUM1 
            DYEDT = (S(NYS+KEL) - SN(NYS+KEL)) * DUM1
            F(NYS+KEL) = F(NYS+KEL) - DTEDT * CVE*V*RHO*S(NYS+KEL)
     1                   + DYEDT * S(NTE)*RHO*V*RU/WT(KEL,1)
            IF (LENRGY) THEN	
               DYTEDT = (S(NTE)*S(NYS+KEL) - SN(NTE)*SN(NYS+KEL))*DUM1
               F(NT) = F(NT) + DYTEDT*(RU/WT(KEL,1))/CPB
            ENDIF
         ENDIF
C
         IDTOT = 0
         DO 8200 IM = 1, NMAT
            IF (NSPHCH(IM) .GT. 0) THEN
               DO 8100 IDM = 1, NSPHCH(IM)
                  I = IDM + IDTOT
                  DYDT = (S(NSDENS+I) - SN(NSDENS+I)) * DUM1
                  F(NSDENS+I) = F(NSDENS+I) - DYDT
8100           CONTINUE
            ENDIF
            IDTOT = IDTOT + NSPHCH(IM)
8200     CONTINUE
C
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE PSPRNT (LOUT, NSPHCH, KK, KKGAS, KKSURF, KKBULK, NPTS,
     1                   NATJ, NPHASE, NNSURF, NFSURF, NLSURF, NNBULK,
     2                   NFBULK, NLBULK, KKPHAS, KFIRST, KLAST, KOCC,
     3                   KSYM, PSYM, ETCH, LENRGY, LFLRT, LTDIFF, 
     4                   LEQUIV, LPRTIC, HIN, YIN, EQUIV, PA, P, TAU,
     5                   FLRT, V, AREA, Q, HTRN, TAMBNT, LHTRN, TIN,
     6                   XIN, T, TSURF, X, NIK, SCRTCH, SN, S, F,
     8                   ICKWRK, RCKWRK, ISKWRK, RSKWRK, MAPPH, WT,
     8                   ACT, SDEN, SDEN0, SDOT, DEN, SITDOT, GRATE,
     9                   NT, NTE,  NYS, NY, NSDEN, NSDENS, KEL, KKION,
     *                   KCHG, KION, QEK, POWR, LENRGE, TEIN, TE, ROP,
     1                   ROPS, IEIMP, EIONSH, II, NQLSE, QLSE, TQLSE,
     2                   GFAC, SFAC, TIONP, LCONFN, IEXC, REXC, SDOTI,
     4                   SITDTI, IISUR, NMAT, IISUR0, KKMAX, NPHMAX, 
     5                   KKTOT, IMISK, IMRSK, IMCSK, AFRAC, MSYM, 
     6                   ITDEP, BHMXI, FLXION, FLXCOR, SDOTT, ESHTH, 
     7                   LELSH, LRFSH, LRFCUR, RFFREQ, RFAPAR, 
     8                   RFDPAR, LWALHB, VISC, THCOND, LSTEDY, ENDTIM, 
     9                   OLDTIM, LTERRN, DEPRAT, SPUTTR, NBHM,
     *                   LOUTSH, BOLTZ, AVOG, EV2K, BIASPW, HIM)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DOUBLE PRECISION LOWSDN
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      REAL LOWSDN
C*****END precision > single
C
      DIMENSION KKPHAS(NPHMAX,NMAT), KFIRST(NPHMAX,NMAT), 
     1          KLAST(NPHMAX,NMAT), KOCC(KKMAX,NMAT), ICKWRK(*),
     1          ISKWRK(*), MAPPH(NPHMAX,NMAT), HIN(KKMAX), YIN(KKMAX),
     2          XIN(KKTOT), X(KKTOT), SCRTCH(NIK,6), SN(NATJ), S(NATJ),
     3          F(NATJ), RCKWRK(*), RSKWRK(*), WT(KKMAX,NMAT),
     4          ACT(KKMAX,NMAT), SDEN(NPHMAX,NMAT), SDOT(KKTOT),
     5          SDEN0(NPHMAX,NMAT), DEN(KKMAX,NMAT),  ITDEP(*), KION(*), 
     6          SITDOT(NPHMAX,NMAT), GRATE(NPHMAX,NMAT), KCHG(KKMAX), 
     7          TEMP(3), QEK(KKGAS), ROP(*), IEIMP(*), QLSE(NPTS),
     7          TQLSE(NPTS), IEXC(*), REXC(*), SDOTI(KKMAX,IISUR0,NMAT),
     8          SITDTI(NPHMAX,NMAT), IMISK(NMAT), IMRSK(NMAT), 
     9          IMCSK(NMAT),KKSURF(NMAT), KKBULK(NMAT), KK(NMAT),
     *          NPHASE(NMAT), NNSURF(NMAT), NFSURF(NMAT), NLSURF(NMAT),
     *          NNBULK(NMAT), NFBULK(NMAT), NLBULK(NMAT), IISUR(NMAT),
     1          NSPHCH(NMAT), AFRAC(NMAT), EIONSH(NMAT), FLXION(*),
     2          FLXCOR(*), SDOTT(KKTOT), ROPS(IISUR0,NMAT),
     3          TSURF(NMAT), NBHM(NMAT), ESHTH(NMAT), BIASPW(NMAT),
     4          HIM(KKMAX,NMAT)
      CHARACTER*16 KSYM(KKMAX,NMAT), PSYM(NPHMAX,NMAT), MSYM(NMAT)
      LOGICAL LENRGY, LFLRT, LTDIFF(NMAT), LEQUIV, LPRTIC, 
     1        ETCH(NPHMAX,NMAT),
     1        LHTRN, LTIME, LENRGE, LCONFN, LELSH(NMAT), LRFSH(NMAT), 
     2        LRFCUR, LWALHB(NMAT), LSTEDY, LTERRN
      PARAMETER (LOWSDN = 1.0E-10, SMALL = 1.E-25)
      PARAMETER(IDIMRF = 101, UMPMIN = 6.E5)
      DIMENSION KYLD(12), PYLD(4)
      DATA CALERG /2.389E-8/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      WRITE (LOUT, 6000)
C
      IF (.NOT. LSTEDY) WRITE (LOUT,6001) OLDTIM+ENDTIM

C  COMPUTE PHASE DENSITIES
      IDTOT = 0
      DO 7 IM = 1, NMAT
         IF (NSPHCH(IM) .GT. 0) THEN
            DO 5 IDM = 1, NSPHCH(IM)
               I = IDM + IDTOT
               SDEN(MAPPH(IDM,IM),IM) = 
     1               S(NSDENS+I) * SDEN0(MAPPH(IDM,IM),IM)
5           CONTINUE
            IDTOT = IDTOT + NSPHCH(IM)
         ENDIF
7     CONTINUE
      TEMP(1) = S(NT)
      TEMP(2) = S(NTE)
      TEMP(3) = TIONP
      CALL CKRHOY (P, TEMP, S(NY), ICKWRK, RCKWRK, RHO)

C  COMPUTE MOLE FRACTION ARRAY
      CALL CKYTX  (S(NY), ICKWRK, RCKWRK, X)
      KSTOT = 0
      DO 17 IM = 1, NMAT
         IF (NNSURF(IM) .GT. 0 ) THEN
            DO 10 KM = KFIRST(NFSURF(IM),IM), KLAST(NLSURF(IM),IM)
               K = KM + KSTOT
               X(K) = S(NYS+K)
10          CONTINUE
         ENDIF
         IF (NNBULK(IM) .GT. 0) THEN
            DO 15 KM = KFIRST(NFBULK(IM),IM), KLAST(NLBULK(IM),IM)
               K = KM + KSTOT
               X(K) = S(NYS+K)
15          CONTINUE
         ENDIF
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
17    CONTINUE
C
C COMPUTE ACTIVITY COEFFICIENT ARRAY
      KSTOT = 0
      DO 19 IM = 1, NMAT
         DO 18 KM = 1, KK(IM)
            IF (KM .GT. KKGAS) THEN
               K = KM + KSTOT
            ELSE
               K = KM
            ENDIF
            ACT(KM,IM) = X(K)
 18      CONTINUE
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
19    CONTINUE
      CALL PSRACT (X, KK, KKBULK, NFBULK, NLBULK, NPHASE, KFIRST,
     1             KLAST, NNBULK, NMAT, NPHMAX, KKSURF, KKMAX, ACT)
C
C CALCULATE SDOT AND SITDOT ARRAYS
C
      CALL PSSDOT (KKION, KEL, KION, KKGAS, KKTOT, KCHG, NMAT, NPHASE,
     1             IISUR, KFIRST, KLAST, NFSURF, NLSURF, NNSURF, KKSURF,
     2             NFBULK, NLBULK, NNBULK, KKBULK, IISUR0, KKMAX,
     3             NPHMAX, KK, BHMXI, WT, EIONSH, NIK,
     3             AFRAC, RFFREQ, RFAPAR, RFDPAR, LRFSH, LRFCUR,
     4             P, ACT, SCRTCH(1,5), TIONP, S(NTE), TSURF, SDEN,
     5             ISKWRK, RSKWRK, IMISK, IMRSK, FLXCOR, FLXION, ESHTH,
     7             LELSH, SFAC, SDOTT, ROPS, SDOTI, SITDTI, SITDOT, 
     8             SDOT, SHLOSS, NBHM, LOUTSH, BOLTZ, AVOG, EV2K,
     9             BIASPW, AREA, HIM(1,1))

      DEPRAT = 0.0
      SPUTTR = 0.0
      KSTOT = 0
      DO 80 IM = 1, NMAT
C
C CALCULATE DEN ARRAY (note only valid for surface and bulk phase)
C
         CALL SKDEN (P, TEMP, ACT(1,IM), SDEN(1,IM), ISKWRK(IMISK(IM)), 
     1               RSKWRK(IMRSK(IM)), DEN(1,IM))
C
C     
C CALCULATE NET DEPOSITION VS. SPUTTR RATES FOR TERRAIN INPUT
C
         DO 74 I = 1, IISUR(IM)
            CALL SKIYLD (I, ISKWRK, RSKWRK, IYLD, IYION, KYLD, PYLD)
            IF (NNBULK(IM) .GT. 0) THEN
               DO 73 KM = KFIRST(NFBULK(IM),IM), KLAST(NLBULK(IM),IM)
                  K = KM + KSTOT
                  IF (LTERRN) THEN
                     IF (IYLD.GT.0 .AND. SDOTI(KM,I,IM) .LT. 0.0) THEN
                        SPUTTR = SPUTTR - UMPMIN *
     1                           WT(KM,IM) / DEN(KM,IM)*SDOTI(KM,I,IM)
                     ELSE
                        DEPRAT = DEPRAT + UMPMIN *
     1                           WT(KM,IM) / DEN(KM,IM)*SDOTI(KM,I,IM)
                     ENDIF
                  ENDIF
 73            CONTINUE
            ENDIF
 74      CONTINUE
C
C CALCULATE THE GROWTH RATE 
C
         IF (NNBULK(IM) .GT. 0) THEN
            DO 75 IPHASE = NFBULK(IM), NLBULK(IM)
               GRATE(IPHASE,IM) = 0.0
               DO 75 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                  K = KM + KSTOT
                  GRATE(IPHASE,IM) = GRATE(IPHASE,IM) +
     1                         WT(KM,IM) / DEN(KM,IM) * SDOT(K)
 75            CONTINUE
         ENDIF
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
 80   CONTINUE
C
      IF (LFLRT) THEN
         TAU = RHO * V / FLRT
      ELSE
         FLRT = RHO * V / TAU
      ENDIF
C
C           PRINT SOLUTION
C
      IF (LEQUIV) WRITE (LOUT, 6020) EQUIV
      WRITE (LOUT, 6030) TAU
      WRITE (LOUT, 6040) FLRT
      WRITE (LOUT, 6050) PA
      WRITE (LOUT, 6055) RHO
      WRITE (LOUT, 6060) V
      WRITE (LOUT, 6065) AREA
      WRITE (LOUT, 6066) AREA/V
C
      IF (LENRGY) THEN
         WRITE (LOUT, 6070) TIN
         WRITE (LOUT, 6080) S(NT)

         IF (LHTRN) THEN
           WRITE (LOUT,6089) AREA, HTRN*CALERG, TAMBNT
         ELSE
           QCAL = Q * CALERG
           WRITE (LOUT, 6090) QCAL
         END IF
      ELSE
         WRITE (LOUT, 6100) T
      ENDIF
      IF (LENRGE .AND. KEL .NE. 0) THEN
         WRITE (LOUT, 6075) TEIN
         WRITE (LOUT, 6085) S(NTE)
         WRITE (LOUT, 6095) POWR
      ELSEIF (KEL .NE. 0) THEN
         WRITE (LOUT, 6098) TE
      ENDIF
      DO 225 IM = 1, NMAT
         IF (LTDIFF(IM)) THEN
            WRITE (LOUT, 6102) MSYM(IM), TSURF(IM)
         ELSE
            WRITE (LOUT, 6103) MSYM(IM)
         ENDIF
 225  CONTINUE
C
      IF (LPRTIC) THEN
         WRITE (LOUT, 6110)
         CALL PSPRT1 (LOUT, KKGAS, KSYM(1,1), XIN)
         LPRTIC = .FALSE.
      ENDIF
C
      WRITE (LOUT, 6120)
      CALL PSPRT1 (LOUT, KKGAS, KSYM(1,1), X)
      KSTOT = 0
      DO 250 IM = 1, NMAT
        IF (NPHASE(IM) .GT. 1) THEN
           IF (NMAT .GT. 1) WRITE (LOUT, 6125) MSYM(IM), AFRAC(IM)*100.
           IF (NNSURF(IM) .GT. 0) THEN
             DO  100 IPHASE = NFSURF(IM), NLSURF(IM)
               IF (ABS(SDEN(IPHASE,IM)) .EQ. 0.0) THEN
                  WRITE (LOUT, 6131) PSYM(IPHASE,IM), 
     1                         LOWSDN*SDEN0(IPHASE,IM),
     2                         SDEN0(IPHASE,IM), SITDOT(IPHASE,IM)
               ELSE
                  WRITE (LOUT, 6130) PSYM(IPHASE,IM), 
     1                         SDEN(IPHASE,IM),  SDEN0(IPHASE,IM), 
     2                         SITDOT(IPHASE,IM)
               ENDIF
               K1 = KFIRST(IPHASE,IM) + KSTOT
               CALL PSPRT1 (LOUT, KKPHAS(IPHASE,IM),
     1              KSYM(KFIRST(IPHASE,IM),IM), X(K1))
100          CONTINUE
           ENDIF
           IF (NNBULK(IM) .GT. 0) THEN
             DO  200 IPHASE = NFBULK(IM), NLBULK(IM)
               DENPH = 0.0
               WTPH  = 0.0
               DO 140 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                  K = KM + KSTOT
                  WTPH  = WTPH  + WT(KM,IM) * X(K)
                  DENPH = DENPH + WT(KM,IM) * X(K) /DEN(KM,IM)
140            CONTINUE
               DENPH = WTPH / DENPH
               WRITE (LOUT, 6140) PSYM(IPHASE,IM), 
     1               GRATE(IPHASE,IM), GRATE(IPHASE,IM)*AREA*DENPH, 
     2               DENPH, WTPH
               IF (GRATE(IPHASE,IM) .LE. 0.0 
     1                              .AND. (.NOT. ETCH(IPHASE,IM)))
     1            WRITE (LOUT, 6143)
               IF (GRATE(IPHASE,IM) .GT. 0.0 .AND. ETCH(IPHASE,IM) )
     1            WRITE (LOUT, 6144)
               WRITE (LOUT, 6146)
               DO 150 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                  K = KM + KSTOT
                  WRITE (LOUT, 6150) KSYM(KM,IM), X(K), ACT(KM,IM),
     1            DEN(KM,IM), SDOT(K), SDOT(K)*WT(KM,IM),
     1            SDOT(K)*WT(KM,IM)/DEN(KM,IM),
     1            SDOT(K)*WT(KM,IM)/DEN(KM,IM)*3600.*1.0E4
150            CONTINUE
200          CONTINUE
           ENDIF
        ENDIF
        KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
250   CONTINUE
C
C compute charge imbalance
C
      IF (KEL.NE.0. .AND. KKION .NE. 0) THEN
         SUMXI = 0.0
         SUMPOS = 0.0
         DO 280 KI = 1, KKION
            K = KION(KI)
            SUMXI = SUMXI + X(K)*FLOAT(KCHG(K))
            IF (KCHG(K).GT.0) SUMPOS = SUMPOS + X(K)
280      CONTINUE
         XIMB = SUMXI - X(KEL)
         PIMB = 100.* XIMB / SUMPOS 
         WRITE (LOUT, 6122) PIMB
C
C write out sheath average ion energy if solved for
C
         DO 283 IM = 1, NMAT
            IF (LRFSH(IM) .OR. BIASPW(IM) .GT. 0.0) THEN
               WRITE (LOUT,6123) MSYM(IM), (EIONSH(IM)/EV2K)
            ENDIF
 283     CONTINUE
      ENDIF
C
C    PRINT OUT RESIDUALS
C
      LTIME = .FALSE.
      CALL PSRFUN (KK, KKGAS, KKSURF, KKBULK, NSPHCH, NPTS, NATJ,
     1             NPHASE, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2             NLBULK, KKPHAS, KFIRST, KLAST, KOCC, ETCH, LTIME,
     3             LENRGY, LFLRT, LTDIFF, DT, HIN, YIN, TIN, T,
     4             TSURF, P, TAU, FLRT, V, AREA, Q, HTRN, TAMBNT,
     5             LHTRN, ICKWRK, RCKWRK, ISKWRK, RSKWRK, MAPPH, WT,
     6             SCRTCH(1,1), NIK, SN, S, F,
     7             SCRTCH(1,3), SDOT, ACT, SDEN, SDEN0, SITDOT, NT,
     8             NTE, NYS, NY, NSDEN, NSDENS, KEL, KKION, KCHG,
     9             KION, QEK, POWR, LENRGE, TEIN, TE, ROP, ROPS,
     *             IEIMP, EIONSH, II, NQLSE, QLSE, TQLSE, GFAC, SFAC,
     1             TIONP, LCONFN, IEXC, REXC, SDOTI, SITDTI, IISUR,
     2             KKTOT, NMAT, IISUR0, KKMAX, NPHMAX, IMISK, IMRSK,
     3             IMCSK, AFRAC, ITDEP, BHMXI, FLXION, FLXCOR, SDOTT,
     4             SCRTCH(1,5), ESHTH, LELSH, LRFSH,
     5             LRFCUR, RFFREQ, RFAPAR, RFDPAR, LWALHB, VISC,
     6             THCOND, NBHM, LOUTSH, BOLTZ, AVOG, EV2K, BIASPW, 
     7             HIM)
      WRITE (LOUT, 7010)
      WRITE (LOUT, 7011) F(NT)
      IF (KEL .NE. 0) WRITE (LOUT, 7020) F(NTE)
      DO 300 K = 1, KKGAS
         WRITE (LOUT, 7022) KSYM(K,1), F(NYS+K)
300   CONTINUE
      KSTOT = 0
      DO 350 IM = 1, NMAT
         IF (NMAT.GT.1 .AND. (NNSURF(IM)+NNBULK(IM)) .GT.0) 
     1      WRITE (LOUT, 7026) MSYM(IM)
         IF (NNSURF(IM) .GT. 0) THEN
            WRITE (LOUT, 7023)
            DO 310 KM = KFIRST(NFSURF(IM),IM), KLAST(NLSURF(IM),IM)
               K = KM + KSTOT
               WRITE (LOUT, 7022) KSYM(KM,IM), F(NYS+K)
310         CONTINUE
         ENDIF
         IF (NNBULK(IM) .GT. 0) THEN
            WRITE (LOUT, 7024)
            DO 320 KM = KFIRST(NFBULK(IM),IM), KLAST(NLBULK(IM),IM)
               K = KM + KSTOT
               WRITE (LOUT, 7022) KSYM(KM,IM), F(NYS+K)
320         CONTINUE
         ENDIF
         IDTOT = 0
         IF (NSPHCH(IM) .GT. 0) THEN
            WRITE (LOUT, 7025)
            DO 330 ID = 1 , NSPHCH(IM)
               I = ID + IDTOT
               IPHASE = MAPPH(ID,IM)
               WRITE (LOUT, 7022) PSYM(IPHASE,IM), F(NSDENS+I)
330         CONTINUE
         ENDIF
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
         IDTOT = IDTOT + NSPHCH(IM)
350   CONTINUE
      WRITE (LOUT, 7000)
C
 6000 FORMAT (//120('=')/30X,
     1 'PSPRNT: Printing of current solution from TWOPNT:',
     1        /120('=')/)
 6001 FORMAT (//3X,'TWOPNT Transient Solution at Time = ',1PE12.4,
     1        ' seconds')
 6020 FORMAT (/15X,'FUEL EQUIVALENCE RATIO ',7X,F11.4)
 6030 FORMAT (/15X,'RESIDENCE TIME ',15X,E11.4,2X,'SEC')
 6040 FORMAT (15X,'MASS FLOW RATE ',15X,E11.4,2X,'GM/SEC')
 6050 FORMAT (15X,'PRESSURE ',21X,G11.4,2X,'ATM')
 6055 FORMAT (15X,'MASS DENSITY ',17X,G11.4,2X,'GM/CM3')
 6060 FORMAT (15X,'VOLUME ',23X,G11.4,2X,'CM3')
 6065 FORMAT (15X,'SURFACE AREA ',17X,G11.4,2X,'CM2')
 6066 FORMAT (15X,'SURFACE TO VOLUME RATIO',7X,G11.4,2X,'CM-1')
 6070 FORMAT (15X,'TEMPERATURE (INLET) ',10X,F11.4,2X,'K')
 6075 FORMAT (15X,'ELECTRON TEMPERATURE (INLET) ',1X,F11.4,2X,'K')
 6080 FORMAT (15X,'TEMPERATURE ',18X,G11.4,2X,'K')
 6085 FORMAT (15X,'ELECTRON TEMPERATURE ',9X,G11.4,2X,'K')
 6089 FORMAT (15X,'HEAT LOSS ', 20X,'(',G9.2,' CM2) (',G11.4,
     1       ' CAL/(CM2*K*SEC)) ( T -',G9.2,' K)')
 6090 FORMAT (15X,'HEAT LOSS ', 20X, G11.4,2X,'CAL/SEC')
 6095 FORMAT (15X,'ELECTRON POWER COUPLING ',6X,G11.4,2X,'ERG/S')
 6098 FORMAT (15X,'ELECTRON TEMPERATURE (FIXED) ',1X,F11.4,2X,'K')
 6100 FORMAT (15X,'TEMPERATURE (FIXED) ',10X,G11.4,2X,'K')
 6102 FORMAT (15X,'SURF TEMP (different than gas) ON MATERIAL ',
     1        A15,' = ', F11.4,2X,'K')
 6103 FORMAT (15X,'SURF TEMP same as gas temp ON MATERIAL ', A15)
 6110 FORMAT (//3X,'INLET  GAS PHASE MOLE FRACTIONS ',/)
 6120 FORMAT (//3X,'EXIT   GAS PHASE MOLE FRACTIONS ',/)
 6122 FORMAT (//15X,'CHARGE IMBALANCE FOR THIS SOLUTION = ',
     1        E12.5,' % OF POSITIVE CHARGE DENSITY'/)
 6123 FORMAT (//15X,'THE AVERAGE ION ENERGY COMPUTED AT MATERIAL ',
     1        A15,' = ', E12.5,' VOLTS',/)
 6125 FORMAT (//3X,'SURFACE AND BULK FRACTIONS FOR MATERIAL, ',A10,
     1        ',',F11.4,' % OF TOTAL SURFACE AREA')
 6130 FORMAT (/21X,'SURFACE SITE FRACTIONS IN SURFACE PHASE, ',
     1        A10/
     1        35X,'Site density = ',g11.4,' mole/cm**2'/
     1        35X,'Standard State Site density = ',g11.4,' mole/cm**2'/
     1        35x,'Rate of change of site density = ',g11.4,
     1             ' mole/(cm**2*sec)'/)
 6131 FORMAT (/21X,'SURFACE SITE FRACTIONS IN SURFACE PHASE, ',
     1        A10/
     1        35X,'Site density = 0.0 mole/cm**2'/
     1        35X,'Site density used in the calculation of surface ',
     1            'site fractions = ',G11.4,' mole/cm**2'/
     1        35X,'Standard State Site density = ',g11.4,' mole/cm**2'/
     1        35x,'Rate of change of site density = ',g11.4,
     1             ' mole/(cm**2*sec)'/)
 6140 FORMAT (/21X,'BULK PHASE MOLE FRACTIONS AND ACTIVITIES IN',
     1       ' BULK PHASE, ', A10/
     1  35X,'Total growth rate of bulk phase = ',G11.4,' cm/sec'/
     1  35X,'                                = ',G11.4,' gm/sec'/
     1  35X,'Density of the bulk phase = ',G11.4,'gm/cm**3'/
     1  35X,'Average molecular weight of bulk phase = ',G11.4,
     1      ' gm/mole')
6143  FORMAT(35X,'WARNING: Growth rate is negative, but an etching',
     1  ' keyword was not input for this phase')
6144  FORMAT(35X,'WARNING: Growth rate is positive, but an etching',
     1  ' reaction was input for this phase')
6146  FORMAT(//
     1 ' Species Name   Mole_frac     Activity      Density      ',
     1 '------------------------Growth Rate-----------------------'
     1 /
     1 '                                          (gm/cm**3)    ',
     1 'mole/(cm**2*sec)  gm/(cm**2*sec))   cm/sec  (microns/hr)'
     1 /)
 6150 FORMAT (3X, A10,'= ', G11.4, 2X, G11.4, 2X, G11.4, 3X, G11.4,
     1        3X, G11.4, 3X, G11.4, 3X, G11.4 )
 7000 FORMAT (/120('=')//)
 7010 FORMAT (//30X,'PRINTOUT OF  RAW RESIDUALS:'/)
 7011 FORMAT (10X,'Temperature     : ',G11.3)
 7020 FORMAT (10X,'Electron temperature     : ',G11.3
     1     /5X,'Gas phase mass fraction residuals:')
 7022 FORMAT (10X,A16,': ',G11.3)
 7023 FORMAT (5X,'Surface phase site fraction residuals:')
 7024 FORMAT (5X,'Bulk phase mole fraction residuals:')
 7025 FORMAT (5X,'Site Density of surface phase residuals:')
 7026 FORMAT (5X,'Material:  ',A)
C
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE PSPRT1 (LOUT, KK, KSYM, X)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION X(KK)
      CHARACTER*16 KSYM(*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 10 K = 1, KK, 3
         WRITE (LOUT, 6010) (KSYM(L), X(L), L=K, MIN(K+2, KK))
   10 CONTINUE
 6010 FORMAT (3X,3(A17,'= ', G11.4, 6X))
C
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE PSRKEY (LIN, LOUT, MM, KK, KKGAS, KKSURF, KKBULK, NATJ,
     1                   NSPHCH, NPHASE, NUMPSR, MAXPSR, NNSURF, NFSURF,
     2                   NLSURF, NNBULK, NFBULK, NLBULK, ICKWRK, RCKWRK,
     4                   LENISK, ISKWRK, LENSK, RSKWRK, LENCK, LENICK,
     4                   KSYM, ATOM, PSYM, NCF, LENRGY, LFLRT, LEQUIV,
     4                   LTDIFF, LRSTRT, LCNTUE, IPRNT, LSEN, LSENT,
     5                   LSENG, ISEN, EPSS, EPST, EPSG, ABSOL, RELAT,
     6                   LROP, IROP, LESTIM, ETCH, KKPHAS, KFIRST, 
     7                   KLAST, EPSR, TIN, XIN, T, TSURF, XEST, EQUIV,
     9                   PATM, PA, P, TAU, FLRT, V, Q, HTRN, TAMBNT,
     *                   AREA, LHTRN, SDEN, FUEL, OXID, ADD, STOICH,
     1                   KPROD, KSP, ATOL, RTOL, ATIM, RTIM, IRETIR,
     2                   NINIT, NJAC, ITJAC, NUMDT, DT, NUMDT2, DT2,
     3                   SFLR, UFAC, DFAC, DTMIN, DTMAX, LENIEQ,
     4                   LENEQ, IEQWRK, EQWRK, IPVT, NIK, RHS, A, NT,
     5                   NTE, NYS, NY, NSDEN, NSDENS, WT, KEL, QEK,
     6                   POWR, LENRGE, TEIN, TE, EIONSH, LNOFT, NPTS,
     8                   NQLSE, QLSE, TQLSE, GFAC, SFAC, TIONP, LCONFN,
     9                   SPOS, NPHMAX, KKMAX, NMAT, KKTOT, NSPTOT,
     *                   AFRAC, MSYM, LTION, LELSH, BHMXI, ESHTH, 
     1                   LRFSH, LRFCUR, RFFREQ, RFAPAR, RFDPAR,
     2                   LWALHB, VISC, THCOND, LSTEDY, ENDTIM, OLDTIM,
     3                   LTERRN, KERR, EV2K, MUSD, LSCCM, SCCMIN,
     4                   BIASPW)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KKSURF(NMAT), KKBULK(NMAT), KK(NMAT), NPHASE(NMAT), 
     1          NNSURF(NMAT), NFSURF(NMAT), NLSURF(NMAT), NNBULK(NMAT),
     2          NFBULK(NMAT), NLBULK(NMAT), NSPHCH(NMAT), TIN(MAXPSR),
     3          T(MAXPSR), TSURF(NMAT,MAXPSR), PA(MAXPSR), P(MAXPSR),
     1          TAU(MAXPSR), FLRT(MAXPSR), V(MAXPSR), Q(MAXPSR),
     2          HTRN(MAXPSR), AREA(MAXPSR), TAMBNT(MAXPSR),
     3          SDEN(NPHMAX,NMAT), FUEL(KKMAX), OXID(KKMAX), ADD(KKMAX),
     4          STOICH(KKMAX), XIN(KKTOT,MAXPSR), XEST(KKTOT), RHS(NIK),
     4          A(NATJ,NATJ), RCKWRK(LENCK), EQWRK(LENEQ),
     5          ISKWRK(LENISK), RSKWRK(LENSK), WT(KKMAX,NMAT),
     6          ICKWRK(LENICK), NCF(MM,KKMAX,NMAT), KPROD(KKMAX),
     7          IPVT(NATJ), KSP(KKMAX), IEQWRK(LENIEQ), TIONP(MAXPSR),
     7          KKPHAS(NPHMAX,NMAT), KFIRST(NPHMAX,NMAT),
     8          KLAST(NPHMAX,NMAT), VALUE(5), TQ(2), PQ(2), HQ(2), 
     9          VQ(2), SQ(2), CQ(2), WM(2), XCON(1), KCON(1), TEMP(3),
     *          QEK(KKGAS), POWR(MAXPSR), TEIN(MAXPSR), TE(MAXPSR),
     1          EIONSH(NMAT,MAXPSR), BIASPW(NMAT,MAXPSR), QLSE(NPTS), 
     2          TQLSE(NPTS), AFRAC(NMAT,MAXPSR), GFAC(MAXPSR), 
     1          SFAC(MAXPSR), BHMXI(MAXPSR), ESHTH(NMAT,MAXPSR),MUSD(MM)
      LOGICAL LESTIM(NPHMAX,NMAT), ETCH(NPHMAX,NMAT,MAXPSR), LRSTRT,
     1        LCNTUE, LFLRT(MAXPSR), LEQUIV, LSEN, LSENT, LSENG, 
     2        LTDIFF(NMAT,MAXPSR), ISEN(KKTOT), LROP, IROP(KKTOT), 
     3        LNOFT, CNTNUD, LENRGY(MAXPSR), IERR, KERR,
     4        LHTRN(MAXPSR), NEC(17), LSCCM, ENRGIN, LABSOL, LEQST,
     5        LENRGE(MAXPSR), LEQPRN, LEQCNT, LCONFN, LTION(MAXPSR), 
     6        LELSH(NMAT,MAXPSR), LRFSH(NMAT,MAXPSR), LRFCUR, 
     6        LWALHB(NMAT,MAXPSR),LSTEDY, LTERRN
C
      CHARACTER*16  KSYM(KKMAX,NMAT), ATOM(MM), PSYM(NPHMAX,NMAT), 
     1              MSYM(NMAT)
      CHARACTER  KEYWRD*4, KSTR*4, LINE*80
C
C equilibrium code variables:
      DIMENSION KKPEQ(1), KFPEQ(1), KLPEQ(1), ACTEQ(1), PDENEQ(1),
     1          VINPUT(6), IRGWRK(1), RRGWRK(1)
      CHARACTER*16 CEQWRK(500)
C
C                                 EXTERNALS
      CHARACTER*4 CKCHUP
      EXTERNAL CKCHUP
      INTEGER PSRPID
      EXTERNAL PSRPID
C                                 DATA STATEMENTS
      DATA CALERG /2.389E-8/, CNTNUD/.FALSE./
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C INPUT-
C
C   MM     - NUMBER OF ELEMENTS.
C   KK     - TOTAL NUMBER OF CHEMICAL SPECIES.
C   KKGAS  - NUMBER OF GAS PHASE SPECIES
C   KKSURF - NUMBER OF SURFACE SPECIES
C   KKBULK - NUMBER OF BULK PHASE SPECIES
C   NPHASE - NUMBER OF PHASES
C   KSYM   - CHEMKIN SPECIES NAMES.
C              DIMENSION KSYM(*) AT LEAST KK.
C   PSYM   - CHEMKIN PHASE NAMES.
C              DIMENSION PSYM(*) AT LEAST NPHASE
C   LENIEQ - LENGTH OF THE INTEGER WORKSPACE FOR THE EQUILIBRIUM
C            SUBROUTINES.
C   LENEQ  - LENGTH OF THE FLOATING POINT WORKSPACE FOR THE EQUILIBRIUM
C            SUBROUTINES.
C   LIN    - UNIT FOR READING KEYWORD INPUT.
C   LOUT   - UNIT FOR PRINTED OUTPUT.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES. NATJ=KK+2.
C   PATM   - PRESSURE OF ONE ATMOSPHERE.
C              CGS UNITS - DYNES/CM**2
C   NT     - POINTER FOR TEMPERATURE UNKNOWN IN THE SOLUTION VECTOR
C   NYS, NY- POINTERS FOR GAS PHASE MOLE FRACTION UNKNOWNS IN THE
C            SOLUTION VECTOR
C
C WORK AND SCRATCH SPACE-
C   A      - JACOBIAN OF STOICHIOMETRIC EQUATION.
C              DIMENSION A(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION
C              AND AT LEAST NATJ FOR THE SECOND.
C   ADD    - ARRAY OF ADDED SPECIES MOLE FRACTIONS.
C              DIMENSION ADD(*) AT LEAST KKGAS.
C   RCKWRK - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   EQWRK  - THE FLOATING POINT WORKSPACE FOR THE EQUILIBRIUM
C            SUBROUTINES.
C              DIMENSION EQWRK(*) AT LEAST LENEQ.
C   FUEL   - ARRAY OF FUEL SPECIES MOLE FRACTIONS.
C              DIMENSION FUEL(*) AT LEAST KKGAS.
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   IEQWRK - THE INTEGER WORKSPACE FOR THE EQUILIBRIUM SUBROUTINES.
C              DIMENSION IEQWRK(*) AT LEAST LENIEQ.
C   IPVT   - ARRAY OF PIVOTS FOR BALANCING THE STOICHIOMETRIC EQUATION.
C              DIMENSION IPVT(*) AT LEAST KK.
C   KPROD  - ARRAY OF PRODUCT SPECIES NUMBERS.
C              DIMENSION KPROD(*) AT LEAST KK.
C   KSP    - ARRAY CONTAINING SPECIES NUMBERS FOR SELECTED SPECIES.
C              DIMENSION KSP(*) AT LEAST KK.
C   OXID   - ARRAY OF OXIDIZER SPECIES MOLE FRACTIONS.
C              DIMENSION OXID(*) AT LEAST KK.
C   RHS    - ARRAY CONTAINING RIGHT-HAND-SIDE OF STOICHIOMETRIC
C            EQUATION.
C              DIMENSION AT LEAST KK.
C   STOICH - ARRAY OF REACTANT SPECIES MOLE FRACTIONS AT STOICHIOMETRIC
C            CONDITIONS.
C              DIMENSION STOICH(*) AT LEAST KKGAS.
C
C OUTPUT-
C   ATIM   - ABSOLUTE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION
C            AS USED FOR THE TIME STEPS.
C   ATOL   - ABSOLUTE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION.
C   DT     - SIZE OF THE TIME STEPS.
C              CGS UNITS - SEC
C   DT2    - SIZE OF THE TIME STEPS AFTER THE ENERGY EQUATION IS
C            INCLUDED.
C              CGS UNITS - SEC
C   EPSR   - TRESHOLD VALUE FOR RATE-OF-PRODUCTION COEFFICIENTS.
C   EPSS   - TRESHOLD VALUE FOR SPECIES SENSITIVITY COEFFICIENTS.
C   EPST   - TRESHOLD VALUE FOR TEMPERATURE SENSITIVITY COEFFICIENTS.
C   EPSG   - TRESHOLD VALUE FOR GROWTH RATE SENSITIVITY COEFFICIENTS.
C   EQUIV  - THE FUEL EQUIVALENCE RATIO.
C   ETCH   - Logical matrix(nphase,nmat,maxpsr) - TRUE if corresponding
C            bulk phase in corresponding reactor is being etched.
C   FLRT   - THE MASS FLOW RATE.
C              CGS UNITS - GM/SEC
C              Dimension FLRT(*) at least MAXPSR in calling program
C   HTRN   - Heat transfer coefficient
C              CGS units - gm/(sec**3*K)
C              Dimension HTRN at least MAXPSR in calling program
C   IPRNT  - FLAG TO SPECIFY THE AMOUNT OF PRINTING.
C              IPRNT = 0, PRINT ONLY THE SOLUTION AFTER CONVERGENCE
C              IPRNT = 1, PRINT THE NORMS AFTER EACH ITERATION.
C              IPRNT = 2, PRINT THE FULL SOLUTION AFTER EACH ITERATION.
C   IROP   - LOGICAL ARRAY. IF IROP(K)=.TRUE. THEN THE RATE-OF-PRODUC-
C            TION COEFFICIENTS OF SPECIES K WILL BE PRINTED OUT.
C              DIMENSION IROP(*) AT LEAST KK.
C   ISEN   - LOGICAL ARRAY. IF ISEN(K)=.TRUE. THEN THE FIRST ORDER
C            SENSITIVITY COEFFICIENTS OF SPECIES K WITH RESPECT TO
C            RATE CONSTANTS WILL BE PRINTED OUT.
C              DIMENSION ISEN(*) AT LEAST KK.
C   LCNTUE - IF LCNTUE=.TRUE. THEN A CONTINUATION PROBLEM WILL FOLLOW.
C            IF LCNTUE=.FLASE. THEN THIS IS THE ONLY PROBLEM FOR THE
C            RUN.
C   LENRGY - IF LENRGY=.TRUE. THEN THE ENERGY EQUATION IS SOLVED.
C            IF LENRGY=.FALSE. THEN A SPECIFIED TEMPERATURE IS USED.
C   LEQUIV - IF LEQUIV=.TRUE. THEN THE INLET COMPOSITION IS DEFINED BY
C              THE FUEL EQUIVALENCE RATIO, THE FUEL AND OXIDIZER
C              COMPOSITIONS AND THE PRODUCTS.
C            IF LEQUIV=.FALSE. THEN THE INLET COMPOSITION IS DEFINED
C              DIRECTLY BY THE USER.
C   LFLRT  - IF LFLRT=.TRUE. THEN THE MASS FLOW RATE IS SPECIFIED AND
C              THE RESIDENCE TIME HAS TO BE CALCULATED.
C            IF LFLRT=.FALSE. THEN THE RESIDENCE TIME IS SPECIFIED AND
C              THE MASS FLOW RATE HAS TO BE CALCULATED.
C   LHTRN  - Logical vector of flags indicating whether heat flux
C            out of the reactor is specified with a heat transfer
C            coefficient or whether with a constant heat flux.
C            Dimension LHTRN(*) at least MAXPSR in the calling program.
C   LTDIFF - IF LTDIFF=.TRUE. THEN A DISTINCT SURFACE TEMPERATURE HAS
C              BEEN SPECIFIED (DISTINCT FROM THE GAS PHASE TEMPERATURE)
C              THIS TEMPERATURE IS USED FOR THE SURFACE PHASE REACTION
C              RATES IRRESPECTIVE OF WHETHER THE GAS PHASE TEMPERATURE
C              IS FIXED OR CALCULATED FROM THE ENERGY EQUATION.
C   LROP   - IF LROP=.TRUE. THEN A RATE-OF-PRODUCTION ANALYSIS
C            IS CARRIED OUT.
C   LRSTRT - IF LRSTRT=.TRUE. THEN START FROM A PREVIOUS PROFILE.
C            IF LRSTRT=.FALSE. THEN START FRESH.
C   LSEN   - IF LSEN=.TRUE. THEN A FIRST ORDER SENSITIVITY ANALYSIS
C            IS CARRIED OUT.
C   LSENT  - IF LSENT=.TRUE. THEN THE TEMPERATURE SENSITIVITY
C            COEFFICIENTS ARE PRINTED OUT.
C   LSENG  - IF LSENG=.TRUE. THEN THE GROWTH RATE SENSITIVITY
C            COEFFICIENTS ARE PRINTED OUT.
C   NUMDT  - NUMBER OF TIME STEPS TO TAKE WHEN DOING A TIME START.
C   NUMDT2 - NUMBER OF TIME STEPS TO TAKE WHEN TIME STEPPING WITH THE
C            ENERGY EQUATION INCLUDED.
C   P      - THE PRESSURE.
C              CGS UNITS - DYNES/CM**2
C            Dimension P(*) at least MAXPSR in the calling program.
C   PA     - THE PRESSURE.
C              UNITS - ATM
C            Dimension PA(*) at least MAXPSR in the calling program.
C   Q      - THE HEAT LOSS OF THE REACTOR
C              CGS UNITS - ERGS/SEC
C            Dimension Q(*) at least MAXPSR in the calling program.
C   RTIM   - RELATIVE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION
C            AS USED FOR THE TIME STEPS.
C   RTOL   - RELATIVE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION.
C   T      - THE USER-SPECIFIED TEMPERATURE.
C              CGS UNITS - K
C            Dimension T(*) at least MAXPSR in the calling program.
C   TAMBNT - Ambient Temperature, for use with the heat transfer
C            coefficient.
C              CGS UNITS - K
C              Dimension TAMBNT(*) at least MAXPSR in the calling
C              program.
C   TAU    - THE NOMINAL RESIDENCE TIME OF THE REACTOR
C              CGS UNITS - SEC
C            Dimension TAU(*) at least MAXPSR in the calling program.
C   TIN    - THE INLET TEMPERATURE.
C              CGS UNITS - K
C            Dimension TIN(*) at least MAXPSR in the calling program.
C   V      - THE VOLUME OF THE REACTOR
C              CGS UNITS - CM**3
C            Dimension V(*) at least MAXPSR in the calling program.
C   AREA   - THE SURFACE AREA OF THE REACTOR
C              CGS UNITS - CM**2 (DEFAULT = 0.0)
C            Dimension AREA(*) at least MAXPSR in the calling program.
C   NUMPSR - NUMBER OF PSRS IN SERIES.  THE PROGRAM ASSUMES THAT
C            THE OUTPUT OF ONE STIRRED TANK IS THE INPUT TO THE
C            NEXT STIRRED TANK.  THE INLET TEMPERATURE OF THE
C            NEXT REACTOR IS THE CALCULATED GAS PHASE TEMPERATURE
C            OF THE CURRENT REACTOR. (DEFAULT = 1)
C   XEST   - ARRAY OF SPECIES MOLE FRACTIONS. STARTING POINT FOR
C            ITERATIONS.  THE FIRST KKGAS ARE THE GAS PHASE MOLE
C            FRACTIONS.  THE NEXT KKSURF ARE THE SURFACE SITE
C            FRACTIONS.  THE LAST KKBULK ARE THE BULK PHASE MOLE
C            FRACTIONS.
C              DIMENSION XEST(*) AT LEAST KK.
C   XIN    - ARRAY OF SPECIES INPUT MOLE FRACTIONS.
C              DIMENSION XIN(*,MAXPSR) AT LEAST KKGAS.
C   AFRAC  - ARRAY OF SURFACE MATERIAL FRACTIONS OF TOTAL AREA
C              DIMENSION AFRAC(NMAT,*) AT LEAST MAXPSR
C   EIONSH - ARRAY OF ION ENERGIES AT SURFACE MATERIALS
C              DIMENSION EIONSH(NMAT,*) AT LEAST MAXPSR
C              CGS UNITS - K
C   BIASPW - ARRAY OF BIAS POWERS APPLIED TO SURFACE MATERIALS
C              CGS UNITS - ERGS/S
C              DIMENSION BIASPW(NMAT,*) AT LEAST MAXPSR
C   LRFSH  - LOGICAL FLAG INDICATING USE OF RF SHEATH MODEL
C             TO DETERMINE ION ENERGIES
C             DIMENSION LRFSH AT LEAST NMAT
C   LRFCUR  - LOGICAL FLAG INDICATING USE OF CURRENT CONTROL FOR
C             RF SHEATH MODEL
C   LWALHB  - LOGICAL FLAG INDICATING HEAT BALANCE TO DETERMINE TSURF
C             DIMENSION LWALHB AT LEAST NMAT
C   RFFREQ  - RF FREQUENCY (HZ) FOR RF SHEATH MODEL
C   RFAPAR  - RF AC PARAMETER (CURRENT OR VOLTAGE) FOR SHEATH MODEL
C   RFDPAR  - RF DC PARAMETER (CURRENT OR VOLTAGE) FOR SHEATH MODEL
C   VISC    - VISCOSITY OF GAS FOR HEAT TRANSFER BALANCE ON WALL
C   THCOND  - THERMAL CONDUCTIVITY OF GAS FOR HEAT TRANSFER BALANCE
C   LSTEDY  - LOGICAL FLAG INDICATING WHETHER OR NOT A STEADY-STATE
C             CALCULATION IS DESIRED IN TWOPNT.
C   ENDTIM  - THE SPECIFIED END TIME FOR A TRANSIENT SOLUTION
C   OLDTIM  - TIME OF SOLUTION INCREMENTED OVER CONTINUATIONS
C   LTERRN  - LOGICAL FLAG TO REQUEST OUTPUT OF FILE FOR TERRAIN
C             TOPOGRAPHICAL SIMULATOR
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
C            INITIALIZE LOCAL VARIABLES
C
      KERR = .FALSE.
      ENRGIN = .FALSE.
      LABSOL = .FALSE.
      NPROD = 0
      MS = 0
      MR = 0
      QDEF = 0.0
      RFVAC = 0.0
      RFIAC = 0.0
      RFVDC = 0.0
C
      DO 100 L = 1, 17
         NEC(L) = .FALSE.
  100 CONTINUE
C
      IF (LCNTUE) THEN
C
C  Continuation Runs:
C      reset some global variables
C      Set local flag, CNTNUD
C
         LCNTUE = .FALSE.
         CNTNUD = .TRUE.
         DO 115 IM = 1, NMAT
            DO 110 I = 1 , NPHASE(IM)
               LESTIM(I,IM) = .FALSE.
110         CONTINUE
115      CONTINUE
         DO 120 K = 1, KKTOT
            XEST(K) = 0.0
 120     CONTINUE
         IF (LRFCUR) THEN
            RFIAC = RFAPAR
         ELSE
            RFVAC = RFAPAR
            RFVDC = RFDPAR
         ENDIF
         OLDTIM = OLDTIM + ENDTIM
      ELSE
C
C  Default values for global variables for
C    normal runs and restart runs.  These values are not needed 
C    for continuation runs.
C
         ATOL   = 1.0E-9
         RTOL   = 1.0E-4
         ATIM   = 1.0E-9
         RTIM   = 1.0E-4
         SFLR =  -1.E-4
         SPOS = -1.0
         NUMDT = 100
         DT    = 1.0E-6
         NUMDT2= 100
         NINIT = 0
         IRETIR = 25
         NJAC   = 20
         ITJAC  = 20
         DT2   = 1.0E-6
         UFAC = 2.0
         DFAC = 2.2
         DTMIN = 1.E-10
         DTMAX = 1.E-04
         ENDTIM = 0.0
         OLDTIM = 0.0
         IPRNT = 1
         LCNTUE = .FALSE.
         LRSTRT = .FALSE.
         LSEN = .FALSE.
         LSENT = .FALSE.
         LSENG = .FALSE.
         LROP = .FALSE.
         LSCCM = .FALSE.
         LNOFT = .FALSE.
         LCONFN = .FALSE.
         LRFCUR = .TRUE.
         LSTEDY = .TRUE.
         LTERRN = .FALSE.
         EPSS = 0.01
         EPST = 1.E-04
         EPSG = 1.E-04
         EPSR = 0.01
         EQUIV = 0.0
         RFFREQ = 13.56E6
         RFAPAR = 0.0
         RFDPAR = 0.0
         VISC   = 0.0
         THCOND = 0.0
         SCCMIN = 0.0
         NUMPSR = 1
         NNU = 0
         NQLSE = 0
         DO 160 N = 1, NPTS
            QLSE(N) = 0.0
            TQLSE(N) = 0.0
 160     CONTINUE
         DO 170 IPSR = 1, MAXPSR
            AREA(IPSR) = 0.0
            Q(IPSR)    = 0.0
            TAU(IPSR)  = 0.0
            FLRT(IPSR) = 0.0
            HTRN(IPSR) = 0.0
            TAMBNT(IPSR) = 298.15
            TIN (IPSR)  = 298.15
            TIONP(IPSR) = 298.15
            LTION(IPSR)  = .FALSE.
            LFLRT(IPSR)  = .TRUE.
            LHTRN(IPSR)  = .FALSE.
            LENRGY(IPSR) = .FALSE.
            LENRGE(IPSR) = .FALSE.
            GFAC(IPSR) = 1.0
            SFAC(IPSR) = 1.0
            BHMXI(IPSR) = 0.0
            DO 168 IM = 1, NMAT
               AFRAC(IM,IPSR) = 1.0
               LTDIFF(IM,IPSR) = .FALSE.
               DO 165 IPHASE = 1, NPHASE(IM)
                 ETCH(IPHASE, IM, IPSR)   = .FALSE.
165            CONTINUE
               LELSH(IM,IPSR) = .FALSE.
               ESHTH(IM,IPSR) = 1.0
               EIONSH(IM,IPSR) = 0.0
               BIASPW(IM,IPSR) = 0.0
               LRFSH(IM,IPSR) = .FALSE.
               LWALHB(IM,IPSR) = .FALSE.
168         CONTINUE
            POWR(IPSR) = 0.0
            TEIN(IPSR) = 298.15
            DO 169 K = 1, KKTOT
               XIN(K,IPSR) = 0.0
 169        CONTINUE
170      CONTINUE
         DO 185 IM = 1, NMAT
            DO 180 IPHASE = 1 , NPHASE(IM)
               LESTIM(IPHASE,IM) = .FALSE.
180         CONTINUE
185      CONTINUE
         DO 190 K = 1, KKGAS
            STOICH(K)=0.0
            FUEL(K)=0.0
            OXID(K)=0.0
            ADD(K)=0.0
            QEK(K)=0.0
190      CONTINUE
         DO 200 K = 1, KKTOT
            XEST(K) = 0.0
            ISEN(K) = .FALSE.
            IROP(K) = .FALSE.
200      CONTINUE
      ENDIF
C
C--------------------------------------------------------------
C
C         READ NEXT INPUT LINE
C
      WRITE (LOUT,'(/A/)') '           KEYWORD INPUT '
C
 1000 CONTINUE
      KEYWRD = ' '
      KSTR = ' '
      LINE = ' '
      READ  (LIN,  7000) KSTR, LINE
      WRITE (LOUT, 8000) KSTR, LINE
      KEYWRD = CKCHUP(KSTR,4)
C
C               IS THIS A KEYWORD COMMENT?
C
      IF (KEYWRD(1:1) .EQ. '.' .OR. KEYWRD(1:1) .EQ. '/'
     1                         .OR. KEYWRD(1:1) .EQ. '!') GO TO 1000
C
C--------------PROBLEM TYPE KEYWORDS--------------------
C
      IF (KEYWRD .EQ. 'TGIV') THEN
C
C         ENERGY EQUATION IS NOT INCLUDED
C
         CALL CKXNUM (LINE, -1, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 0) THEN
            IF (NEC(1)) WRITE (LOUT, *)
     1        'WARNING...BOTH "TGIV" AND "ENRG" GIVEN'
            NEC(1) = .TRUE.
            DO 1010 IPSR = 1, MAXPSR
               LENRGY(IPSR)   = .FALSE.
 1010       CONTINUE
            IERR = .FALSE.
         ELSE
            IPSR = INT(VALUE(1))
            LENRGY(IPSR) = .FALSE.
            IF (IPSR .EQ. 1) NEC(1) = .TRUE.
         ENDIF
         KERR = KERR .OR. IERR
C
C
      ELSEIF (KEYWRD .EQ. 'ENRG') THEN
C
C         ENERGY EQUATION IS INCLUDED
C
         CALL CKXNUM (LINE, -1, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 0) THEN
            IF (NEC(1)) WRITE (LOUT, *)
     1         'WARNING...BOTH "TGIV" AND "ENRG" GIVEN'
            NEC(1)      = .TRUE.
            IF (CNTNUD .AND. (.NOT.LENRGY(1))) ENRGIN = .TRUE.
            DO 1020 IPSR = 1, MAXPSR
               LENRGY(IPSR) = .TRUE.
1020        CONTINUE
            IERR = .FALSE.
         ELSE
            IPSR = INT(VALUE(1))
            LENRGY(IPSR) = .TRUE.
            IF (IPSR .EQ. 1) NEC(1) = .TRUE.
         ENDIF
         KERR = KERR .OR. IERR
C
      ELSEIF (KEYWRD .EQ. 'TEGV') THEN
C
C         ELECTRON ENERGY EQUATION IS NOT INCLUDED
C
         CALL CKXNUM (LINE, -1, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 0) THEN
            IF (NEC(14)) WRITE (LOUT, *)
     1        'WARNING...BOTH "TEGV" AND "ENGE" GIVEN'
            NEC(14) = .TRUE.
            DO 1022 IPSR = 1, MAXPSR
               LENRGE(IPSR)   = .FALSE.
 1022       CONTINUE
            IERR = .FALSE.
         ELSE
            IPSR = INT(VALUE(1))
            LENRGE(IPSR) = .FALSE.
            IF (IPSR .EQ. 1) NEC(14) = .TRUE.
         ENDIF
         KERR = KERR .OR. IERR
C
C
      ELSEIF (KEYWRD .EQ. 'ENGE') THEN
C
C         ELECTRON ENERGY EQUATION IS INCLUDED
C
         CALL CKXNUM (LINE, -1, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 0) THEN
            IF (NEC(14)) WRITE (LOUT, *)
     1         'WARNING...BOTH "TEGV" AND "ENGE" GIVEN'
            NEC(14)      = .TRUE.
            DO 1023 IPSR = 1, MAXPSR
               LENRGE(IPSR) = .TRUE.
 1023       CONTINUE
            IERR = .FALSE.
         ELSE
            IPSR = INT(VALUE(1))
            LENRGE(IPSR) = .TRUE.
            IF (IPSR .EQ. 1) NEC(14) = .TRUE.
         ENDIF
         KERR = KERR .OR. IERR
C
      ELSEIF (KEYWRD .EQ. 'NOFT') THEN
C
C         NO FIXED TEMPERATURE SOLUTION
C
         LNOFT = .TRUE.
C
       ELSEIF (KEYWRD .EQ. 'NSDN') THEN
C
C         DON'T SOLVE SURFACE SITE DENSITY EQUATIONS
C
         NSPTOT = 0
         DO 1024 IM = 1, NMAT
            NSPHCH(IM) = 0
1024     CONTINUE
         NATJ = 2 + KKTOT
C
       ELSEIF (KEYWRD .EQ. 'ETCH') THEN
C
C        SPECIFY THAT A BULK PHASE IS TO BE ETCHED INSTEAD OF GROWN
C
         IPHASE = 0
         MTETCH = 0
         NPH = NPHMAX*NMAT
         CALL CKSNUM (LINE, -1 ,LOUT, PSYM, NPH, 
     1                IPHDUM, NVAL, VALUE, IERR)
         DO 10240 IM = 1, NMAT
            IF (IPHDUM.GT.(NPHMAX*(IM-1)) .AND. 
     1          IPHDUM.LE.(NPHMAX*(IM-1)+NPHASE(IM))) THEN
               IPHDUM = IPHDUM - NPHMAX*(IM-1)
               IF (IPHDUM.GE.NFBULK(IM) .OR. IPHDUM.LE.NLBULK(IM)) THEN
                  IPHASE = IPHDUM
                  MTETCH = IM
                  IVAL = INT(VALUE(1))
               ENDIF
            ENDIF
10240    CONTINUE
         IF (IPHASE .EQ. 0 .AND. NVAL .NE. 0) THEN
           WRITE (LOUT,'(A)')
     1       ' ERROR READING DATA FOR KEYWORD '//KEYWRD
         ELSE
           IF (IPHASE .NE. 0 .AND. NVAL .EQ. 0) THEN
             DO 10250 IPSR = 1, MAXPSR
                ETCH(IPHASE, MTETCH, IPSR) = .TRUE.
10250        CONTINUE
           ELSEIF (IPHASE .EQ. 0 .AND. NVAL .EQ. 0) THEN
             DO 1026 IM = 1, NMAT
                DO 1025 IPSR = 1, MAXPSR
                   IF (NNBULK(IM).GT.0) THEN
                      DO 10255 IPHASE = NFBULK(IM), NLBULK(IM)
                         ETCH(IPHASE, IM, IPSR) = .TRUE.
10255                 CONTINUE
                   ENDIF
 1025           CONTINUE
 1026        CONTINUE
             IERR = .FALSE.
           ELSEIF (NVAL .EQ. 1 .OR. NVAL .EQ. -1) THEN
             IPSR = IVAL
             IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               ETCH(IPHASE, MTETCH, IPSR) = .TRUE.
               IERR = .FALSE.
             ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
             ENDIF
           ENDIF
         ENDIF
         KERR = KERR .OR. IERR
C
C--------------METHOD OPTIONS KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'ATOL') THEN
C
C       ABSOLUTE NEWTON ITERATION CONVERGENCE CRITERIA
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, ATOL, IERR)
         KERR = KERR.OR.IERR
         IF (.NOT. LABSOL) ABSOL = ATOL
C
      ELSEIF (KEYWRD .EQ. 'RTOL') THEN
C
C       RELATIVE NEWTON ITERATION CONVERGENCE CRITERIA
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RTOL, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'ATIM') THEN
C
C       ABSOLUTE NEWTON CONVERGENCE CRITERIA FOR TIMESTEPS
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, ATIM, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'RTIM') THEN
C
C       RELATIVE NEWTON CONVERGENCE CRITERIA FOR TIMESTEPS
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RTIM, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'TIME') THEN
C
C       TIME STEP STARTING PROCEDURE
C
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         NUMDT = INT(VALUE(1))
         DT = VALUE(2)
C
      ELSEIF (KEYWRD .EQ. 'TIM2') THEN
C
C        TIME STEPPING, AFTER ADDING THE ENERGY EQUATION
C
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         NUMDT2 = INT(VALUE(1))
         DT2 = VALUE(2)
C
      ELSEIF (KEYWRD .EQ. 'UFAC') THEN
C
C       TIMESTEP INCREASE WHEN TIMESTEP DOES NOT CHANGE SOLUTION
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, UFAC, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'DFAC') THEN
C
C       TIMESTEP DECREASE WHEN NEWTON FAILS CONVERGENCE ON TIMESTEP
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DFAC, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'DTMN') THEN
C
C       MINIMUM TIMESTEP
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DTMIN, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'DTMX') THEN
C
C       MAXIMUM TIMESTEP
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DTMAX, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'TJAC') THEN
C
C        RETIREMENT AGE OF JACOBIAN DURING TIME STEPPING
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         ITJAC = INT (VALUE(1))
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'ABSL') THEN
C
C      ABSOLUTE PERTURBATION IN FINDING THE NUMERICAL JACOBIAN
C
        CALL CKXNUM (LINE, 1, LOUT, NVAL, ABSOL, IERR)
        LABSOL = .TRUE.
        KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'RELT') THEN
C
C      RELATIVE PERTURBATION IN FINDING THE NUMERICAL JACOBIAN
C
        CALL CKXNUM (LINE, 1, LOUT, NVAL, RELAT, IERR)
        KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'NJAC') THEN
C
C      RETIREMENT AGE OF JACOBIAN DURING STEADY-STATE NEWTON
C
        CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
        NJAC = INT (VALUE(1))
        KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'ISTP') THEN
C
C      NUMBER OF INITIAL TIME STEPS BEFORE NEWTON
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         NINIT = INT (VALUE(1))
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'IRET') THEN
C
C      NUMBER OF TIME STEPS TO TAKE BEFORE CHANGING THE
C        TIME STEP INCREMENT
C
        CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
        IRETIR = INT(VALUE(1))
        KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'STDY') THEN
C
C         TWOPNT WILL SOLVE FOR THE STEADY-STATE SOLUTION (DEFAULT)
C
         LSTEDY = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'TRNS') THEN
C
C         TWOPNT WILL NOT SOLVE FOR THE STEADY-STATE SOLUTION;
C         INSTEAD WILL SOLVE WITH TIMESTEPS TO THE SPECIFIED END TIME
C
         LSTEDY = .FALSE.
         CALL CKXNUM (LINE, -1, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         IF (NVAL.GE.1) ENDTIM = VALUE(1)
C
C--------------REACTOR DEFINITION KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'TEMP') THEN
C
C         TEMPERATURE
C
         NEC(2)  = .TRUE.
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 1030 IPSR = 1, MAXPSR
               T(IPSR) = VALUE(1)
1030        CONTINUE
           IERR = .FALSE.
C
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               T(IPSR) = VALUE(1)
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'ETMP') THEN
C
C         ELECTRON TEMPERATURE
C
         NEC(16)  = .TRUE.
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 1035 IPSR = 1, MAXPSR
               TE(IPSR) = VALUE(1)
 1035       CONTINUE
           IERR = .FALSE.
C
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               TE(IPSR) = VALUE(1)
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'TSRF') THEN
C
C         SURFACE TEMPERATURE FOR SURFACE MATERIAL
C
         CALL CKSNUM(LINE, -2 ,LOUT, MSYM, NMAT, IMAT, NVAL,
     1              VALUE, IERR)
         IF (NVAL.EQ.1) THEN
            IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
               WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
               KERR = .TRUE.
            ELSE
               DO 1036 IPSR = 1, MAXPSR
                  TSURF(IMAT,IPSR) = VALUE(1)
                  LTDIFF(IMAT,IPSR) = .TRUE.
 1036          CONTINUE
            ENDIF
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  TSURF(IMAT,IPSR) = VALUE(1)
                  LTDIFF(IMAT,IPSR) = .TRUE.
               ENDIF
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ELSEIF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'WLHT') THEN
C
C         USE HEAT BALANCE AT WALL MATERIAL TO DETERMINE TSURF
C         ON A SPECIFIED SURFACE MATERIAL;
C         TWO PARAMETERS REQUIRED:  VISCOSITY, THERMAL CONDUCTIVITY
C                                   OF THE GAS (IN ERG-CGS UNITS)
C
         CALL CKSNUM(LINE, -3 ,LOUT, MSYM, NMAT, IMAT, NVAL,
     1              VALUE, IERR)
         IF (NVAL.EQ.2) THEN
            DO 1037 IPSR = 1, MAXPSR
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  LWALHB(IMAT,IPSR) = .TRUE.
                  VISC = VALUE(1)
                  THCOND = VALUE(2)
               ENDIF
 1037       CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 3 .OR. NVAL .EQ. -3) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  AFRAC(IMAT,IPSR) = VALUE(1)
               ENDIF
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ELSEIF (NVAL .LT. 2) THEN
            WRITE (LOUT,'(A)')
     1 ' ERROR...VISCOSITY AND THERMAL CONDUCTIVITY MUST BE SPECIFIED' 
            KERR = .TRUE.
         ELSEIF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'IONE') THEN
C
C          ION ENERGY IN ELECTRON VOLTS FOR SURFACE MATERIAL
C
         CALL CKSNUM(LINE, -2 ,LOUT, MSYM, NMAT, IMAT, NVAL,
     1              VALUE, IERR)
         IF (NVAL.EQ.1) THEN
            DO 1040 IPSR = 1, MAXPSR
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  EIONSH(IMAT,IPSR) = VALUE(1)*EV2K
               ENDIF
 1040       CONTINUE
            IF (LELSH(IMAT,1)) WRITE (LOUT, *)
     1              'WARNING...BOTH IONE AND ELSH SPECIFIED;',
     2              ' ELSH WILL BE USED FOR ELECTRON ENERGY LOSS TERM'
            IF (LRFSH(IMAT,1)) WRITE (LOUT,*)
     1              'WARNING...BOTH IONE AND RFSH SPECIFIED;',
     2              ' IONE WILL BE IGNORED'
            IF (BIASPW(IMAT,1).GT.0.0) WRITE (LOUT,*)
     1              'WARNING...BOTH IONE AND BPWR SPECIFIED;',
     2              ' IONE WILL BE IGNORED'
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  EIONSH(IMAT,IPSR) = VALUE(1)*EV2K
                  IF (LELSH(IMAT,IPSR)) WRITE (LOUT, *)
     1              'WARNING...BOTH IONE AND ELSH SPECIFIED;',
     2              ' ELSH WILL BE USED FOR ELECTRON ENERGY LOSS TERM'
                  IF (LRFSH(IMAT,IPSR)) WRITE (LOUT,*)
     1              'WARNING...BOTH IONE AND RFSH SPECIFIED;',
     2              ' IONE WILL BE IGNORED'
                  IF (BIASPW(IMAT,IPSR).GT.0.0) WRITE (LOUT,*)
     1              'WARNING...BOTH IONE AND BWPR SPECIFIED;',
     2              ' IONE WILL BE IGNORED'
               ENDIF
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ELSEIF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            
         ENDIF
         KERR = KERR.OR.IERR
C
C
      ELSEIF (KEYWRD .EQ. 'BPWR') THEN
C
C          BIAS POWER IN WATTS FOR SURFACE MATERIAL
C          Ion energy will be calculated from Bias Power / Ji*Area
C
         CALL CKSNUM(LINE, -2 ,LOUT, MSYM, NMAT, IMAT, NVAL,
     1              VALUE, IERR)
         IF (NVAL.EQ.1) THEN
            DO 1042 IPSR = 1, MAXPSR
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  BIASPW(IMAT,IPSR) = VALUE(1)*1.E7
               ENDIF
 1042       CONTINUE
            IF (LELSH(IMAT,1)) WRITE (LOUT, *)
     1              'WARNING...BOTH BPWR AND ELSH SPECIFIED;',
     2              ' ELSH WILL BE USED FOR ELECTRON ENERGY LOSS TERM'
            IF (LRFSH(IMAT,1)) WRITE (LOUT,*)
     1              'WARNING...BOTH BPWR AND RFSH SPECIFIED;',
     2              ' BPWR WILL BE IGNORED'
            IF (EIONSH(IMAT,1).NE.0.0) WRITE (LOUT,*)
     1              'WARNING...BOTH BPWR AND IONE SPECIFIED;',
     2              ' IONE WILL BE IGNORED'
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  BIASPW(IMAT,IPSR) = VALUE(1)*1.E7
                  IF (LELSH(IMAT,IPSR)) WRITE (LOUT, *)
     1              'WARNING...BOTH BPWR AND ELSH SPECIFIED;',
     2              ' ELSH WILL BE USED FOR ELECTRON ENERGY LOSS TERM'
                  IF (LRFSH(IMAT,IPSR)) WRITE (LOUT,*)
     1              'WARNING...BOTH BPWR AND RFSH SPECIFIED;',
     2              ' BPWR WILL BE IGNORED'
                  IF (EIONSH(IMAT,IPSR).NE.0.0) WRITE (LOUT,*)
     1              'WARNING...BOTH BPWR AND IONE SPECIFIED;',
     2              ' IONE WILL BE IGNORED'
               ENDIF
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ELSEIF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'ELSH') THEN
C
C         ENERGY LOSS THROUGH SHEATH HEATING OF ELECTRONS AND IONS
C         (VALUE IS MULTIPLIER OF TE), ONLY USED IF IONE NOT USED
C
         CALL CKSNUM(LINE, -2 ,LOUT, MSYM, NMAT, IMAT, NVAL,
     1              VALUE, IERR)
         IF (NVAL.EQ.1) THEN
            DO 1045 IPSR = 1, MAXPSR
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  LELSH(IMAT,IPSR) = .TRUE.
                  ESHTH(IMAT,IPSR) = VALUE(1)
               ENDIF
 1045       CONTINUE
            IF (EIONSH(IMAT,1).NE.0.0) WRITE (LOUT, *)
     1              'WARNING...BOTH IONE AND ELSH SPECIFIED;',
     2              ' ELSH WILL BE USED FOR ELECTRON ENERGY LOSS TERM'
            IF (BIASPW(IMAT,1).NE.0.0) WRITE (LOUT, *)
     1              'WARNING...BOTH BPWR AND ELSH SPECIFIED;',
     2              ' ELSH WILL BE USED FOR ELECTRON ENERGY LOSS TERM'
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  LELSH(IMAT,IPSR) = .TRUE.
                  ESHTH(IMAT,IPSR) = VALUE(1)
                  IF (EIONSH(IMAT,IPSR).NE.0.0) WRITE (LOUT, *)
     1              'WARNING...BOTH IONE AND ELSH SPECIFIED;',
     2              ' ELSH WILL BE USED FOR ELECTRON ENERGY LOSS TERM'
                  IF (BIASPW(IMAT,IPSR).NE.0.0) WRITE (LOUT, *)
     1              'WARNING...BOTH BPWR AND ELSH SPECIFIED;',
     2              ' ELSH WILL BE USED FOR ELECTRON ENERGY LOSS TERM'
               ENDIF
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ELSEIF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'RFSH') THEN
C
C         USE RF SHEATH MODEL TO DETERMINE ION ENERGIES FOR
C         A SPECIFIED SURFACE MATERIAL;
C         INTEGER NUMBER INDICATES VOLTAGE OR CURRENT CONTROL
C
         CALL CKSNUM(LINE, -2 ,LOUT, MSYM, NMAT, IMAT, NVAL,
     1              VALUE, IERR)
         IF (NVAL.LE.1) THEN
            DO 1047 IPSR = 1, MAXPSR
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  LRFSH(IMAT,IPSR) = .TRUE.
                  IF (NVAL.EQ.1) THEN
                     IF (INT(VALUE(1)).EQ.0) THEN
                        LRFCUR = .FALSE.
                     ENDIF
                  ENDIF
               ENDIF
 1047       CONTINUE
            IF (LELSH(IMAT,1)) WRITE (LOUT, *)
     1              'WARNING...BOTH RFSH AND ELSH SPECIFIED;',
     2              ' ELSH WILL BE USED FOR ELECTRON ENERGY LOSS TERM'
            IF (EIONSH(IMAT,1).NE.0.0) WRITE (LOUT,*)
     1              'WARNING...BOTH IONE AND RFSH SPECIFIED;',
     2              ' IONE WILL BE IGNORED'
            IF (BIASPW(IMAT,1).NE.0.0) WRITE (LOUT,*)
     1              'WARNING...BOTH BPWR AND RFSH SPECIFIED;',
     2              ' BPWR WILL BE IGNORED'
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  LRFSH(IMAT,IPSR) = .TRUE.
                  IF (LELSH(IMAT,IPSR)) WRITE (LOUT, *)
     1              'WARNING...BOTH RFSH AND ELSH SPECIFIED;',
     2              ' ELSH WILL BE USED FOR ELECTRON ENERGY LOSS TERM'
                  IF (EIONSH(IMAT,IPSR).NE.0.0) WRITE (LOUT,*)
     1              'WARNING...BOTH IONE AND RFSH SPECIFIED;',
     2              ' IONE WILL BE IGNORED'
                  IF (BIASPW(IMAT,IPSR).NE.0.0) WRITE (LOUT,*)
     1              'WARNING...BOTH BPWR AND RFSH SPECIFIED;',
     2              ' BPWR WILL BE IGNORED'
                  IF (INT(VALUE(1)).EQ.0) THEN
                     LRFCUR = .FALSE.
                  ENDIF
               ENDIF
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ELSEIF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'RFFQ') THEN
C
C         RF FREQUENCY FOR RF SHEATH MODEL
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            RFFREQ = VALUE(1)
            IERR = .FALSE.
         ENDIF
         KERR = KERR .OR. IERR
C
      ELSEIF (KEYWRD .EQ. 'RFVA') THEN
C
C         AC VOLTAGE AMPLITUDE FOR RF SHEATH MODEL (VOLTS)
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            RFVAC = VALUE(1)
            IERR = .FALSE.
         ENDIF
         KERR = KERR .OR. IERR
C
      ELSEIF (KEYWRD .EQ. 'RFVD') THEN
C
C         DC VOLTAGE FOR RF SHEATH MODEL (VOLTS)
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            RFVDC = VALUE(1)
            IERR = .FALSE.
         ENDIF
         KERR = KERR .OR. IERR
C
      ELSEIF (KEYWRD .EQ. 'RFJA' .OR. KEYWRD .EQ. 'RFIA') THEN
C
C         AC CURRENT DENSITY FOR RF SHEATH MODEL (A/M2)
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            RFIAC = VALUE(1)
            IERR = .FALSE.
         ENDIF
         KERR = KERR .OR. IERR
C
C
      ELSEIF (KEYWRD .EQ. 'TION') THEN
C
C         ESTIMATED ION TEMPERATURE (K)
C
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 1048 IPSR = 1, MAXPSR
              TIONP(IPSR) = VALUE(1)
              LTION(IPSR) = .TRUE.
 1048      CONTINUE
           IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               TIONP(IPSR) = VALUE(1)
               LTION(IPSR) = .TRUE.
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'QLSE') THEN
C
C         EXCITATIONAL ELECTRON ENERGY LOSS WITH TEMPERATURE DEPENDENCE
C         (VALUE IS IN eV PER IONIZATION EVENT)
C
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
C
         IF (NQLSE+1 .GT. NPTS) THEN
            WRITE (LOUT, *)
     1      ' ERROR... THIS ARRAY IS ONLY DIMENSIONED FOR ', NPTS,
     2      ' (Te,Qloss-e) PAIRS'
            KERR = .TRUE.
         ELSE
            NQLSE = NQLSE+1
            TQLSE(NQLSE) = VALUE(1)
            QLSE(NQLSE) = VALUE(2)
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'PRES') THEN
C
C         PRESSURE
C
C         IF (NEC(3)) WRITE (LOUT, *)
C     1      'WARNING...BOTH PRES AND PRMT ARE SPECIFIED'
         NEC(3)  = .TRUE.
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 1050 IPSR = 1, MAXPSR
               PA(IPSR) = VALUE(1)
               P(IPSR) = VALUE(1) * PATM
1050        CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               PA(IPSR) = VALUE(1)
               P(IPSR) = VALUE(1)*PATM
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'PRMT') THEN
C
C         PRESSURE in mTorr
C
C         IF (NEC(3)) WRITE (LOUT, *)
C     1      'WARNING...BOTH PRES AND PRMT ARE SPECIFIED'
         NEC(3)  = .TRUE.
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 1052 IPSR = 1, MAXPSR
               PA(IPSR) = VALUE(1) /760.E3
               P(IPSR) = VALUE(1) * PATM /760.E3
1052        CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               PA(IPSR) = VALUE(1) /760.E3
               P(IPSR) = VALUE(1)*PATM /760.E3
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'TAU ') THEN
C
C         RESIDENCE TIME
C
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 1055 IPSR = 1, MAXPSR
               TAU(IPSR) = VALUE(1)
               LSCCM = .FALSE.
               LFLRT(IPSR) = .FALSE.
1055        CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (NEC(4).AND.LFLRT(1)) WRITE (LOUT, *)
     1         'WARNING...BOTH TAU AND FLRT ARE DEFINED'
            LFLRT(IPSR) = .FALSE.
            LSCCM = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               TAU(IPSR) = VALUE(1)
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         NEC(4)  = .TRUE.
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'FLRT') THEN
C
C         MASS FLOW RATE
C
         IF (NEC(4)) WRITE (LOUT, *)
     1      'WARNING...BOTH TAU AND FLRT/SCCM ARE DEFINED'
         NEC(4)  = .TRUE.
         LFLRT(1) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, FLRT(1), IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'SCCM') THEN
C
C         TOTAL FLOW RATE IN TERMS OF SCCM
C
         IF (NEC(4)) WRITE (LOUT, *)
     1      'WARNING...BOTH TAU AND FLRT/SCCM ARE DEFINED'
         NEC(4)  = .TRUE.
         LFLRT(1) = .TRUE.
         LSCCM = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SCCMIN, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'CNFN') THEN
C
C         CONFINE NEGATIVE IONS; NO OUTFLOW OF NEG IONS
C
         LCONFN = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'VOL ') THEN
C
C         VOLUME
C
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            NEC(5)  = .TRUE.
            DO 1060 IPSR = 1, MAXPSR
               V(IPSR) = VALUE(1)
1060        CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .EQ. 1) NEC(5) = .TRUE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               V(IPSR) = VALUE(1)
            ELSE
              WRITE (LOUT, 6002) IPSR, MAXPSR
              IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'AREA') THEN
C
C         Specify the surface area
C
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 1070 IPSR = 1, MAXPSR
              AREA(IPSR) = VALUE(1)
1070        CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               AREA(IPSR) = VALUE(1)
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'QLOS') THEN
C
C         HEAT LOSS
C
         NEC(6) = .TRUE.
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 1080 IPSR = 1, MAXPSR
               Q(IPSR) = VALUE(1) / CALERG
               LHTRN(IPSR) = .FALSE.
1080        CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               Q(IPSR) = VALUE(1) / CALERG
               LHTRN(IPSR) = .FALSE.
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'HTRN') THEN
C
C     HEAT LOSS GIVEN BY HEAT TRANSFER COEFFICIENT AND TAMBNT
C
         NEC(6) = .TRUE.
         CALL CKXNUM (LINE, -3, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            DO 1110 IPSR = 1, MAXPSR
              HTRN(IPSR) = VALUE(1) / CALERG
              TAMBNT(IPSR) = VALUE(2)
              LHTRN(IPSR) = .TRUE.
1110        CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 3 .OR. NVAL .EQ. -3) THEN
            IPSR = INT(VALUE(3))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               HTRN(IPSR) = VALUE(1) / CALERG
               TAMBNT(IPSR) = VALUE(2)
               LHTRN(IPSR) = .TRUE.
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'PWRC') THEN
C
C         POWER COUPLING SOURCE FOR ELECTRON ENERGY IN CALORIES
C
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 1115 IPSR = 1, MAXPSR
               POWR(IPSR) = VALUE(1) / CALERG
 1115       CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               POWR(IPSR) = VALUE(1) / CALERG
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'PWRW') THEN
C
C         POWER COUPLING SOURCE FOR ELECTRON ENERGY IN WATTS
C
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 1120 IPSR = 1, MAXPSR
               POWR(IPSR) = VALUE(1) / 1.E-7
 1120       CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               POWR(IPSR) = VALUE(1) / 1.E-7
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'PWRE') THEN
C
C         POWER COUPLING SOURCE FOR ELECTRON ENERGY IN ERGS
C
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 1125 IPSR = 1, MAXPSR
               POWR(IPSR) = VALUE(1) 
 1125       CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               POWR(IPSR) = VALUE(1) 
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'NPSR') THEN
C
C         NUMBER OF PSR'S IN SERIES
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         NUMPSR = INT(VALUE(1))
         IF (NUMPSR .GT. MAXPSR) THEN
            WRITE (LOUT, *)
     1      'Error...NPSR value exceeds MAXPSR in calling program...'
            IERR = .TRUE.
         ENDIF
         KERR = KERR.OR.IERR
C
C--------------INLET CONDITIONS ------------------------------------
C
      ELSEIF (KEYWRD .EQ. 'TINL') THEN
C
C         INLET TEMPERATURE
C
         NEC(7)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TIN(1), IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'TEIN') THEN
C
C         INLET ELECTRON TEMPERATURE
C
         NEC(17)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TEIN(1), IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'REAC') THEN
C
C         REACTANT COMING INTO FIRST PSR
C
         IF (CNTNUD .AND. (.NOT. NEC(8))) THEN
            DO 300 K = 1, KKGAS
               XIN(K,1) = 0.0
300         CONTINUE
         ENDIF
         NEC(8)  = .TRUE.
         CALL CKSNUM (LINE, 1, LOUT, KSYM(1,1), KKGAS, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ELSE
            LEQUIV = .FALSE.
            XIN(KSPEC,1) = VALUE(1)
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'XSDF') THEN
C
C         DEFAULT VALUE FOR ELECTRON-K COLLISION CROSS-SECTION
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, QDEF, IERR)
         NEC(15) = .TRUE.
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'XSEK') THEN
C
C         ELECTRON-K COLLISION CROSS-SECTION (CM2)
C
         CALL CKSNUM (LINE, 1, LOUT, KSYM(1,1), KKGAS, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR) THEN
            WRITE (LOUT, '(A)')
     1       ' ERROR READING DATA FOR KEYWORD '//KEYWRD
         ELSEIF (KSPEC .NE. KEL) THEN
            QEK(KSPEC) = VALUE(1)
            NNU = NNU + 1
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'EQUI') THEN
C
C         EQUIVALENCE RATIO
C
         NEC(9)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, EQUIV, IERR)
         KERR = KERR.OR.IERR
         LEQUIV = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'FUEL') THEN
C
C         FUEL DEFINITION
C
         NEC(10)  = .TRUE.
         CALL CKSNUM (LINE, 1, LOUT, KSYM(1,1), KKGAS, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ELSE
            FUEL(KSPEC) = VALUE(1)
         ENDIF
         IF (CNTNUD) THEN
            WRITE (LOUT, *) 'ERROR..."FUEL" CANNOT BE REDEFINED'
            KERR = .TRUE.
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'OXID') THEN
C
C         OXIDIZER DEFINITION
C
         NEC(11)  = .TRUE.
         CALL CKSNUM (LINE, 1, LOUT, KSYM(1,1), KKGAS, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ELSE
            OXID(KSPEC) = VALUE(1)
         ENDIF
         IF (CNTNUD) THEN
            WRITE (LOUT, *) 'ERROR..."OXID" CANNOT BE REDEFINED'
            KERR = .TRUE.
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'PROD') THEN
C
C         PRODUCT DEFINITION
C
         NEC(12)  = .TRUE.
         CALL CKCOMP (LINE, KSYM(1,1), KKGAS, KSPEC)
         IF (KSPEC.LE.0 .OR. NPROD+1.GT.KKGAS) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
         ELSE
            NPROD = NPROD + 1
            KPROD(NPROD)=KSPEC
         ENDIF
         IF (CNTNUD) THEN
            WRITE (LOUT, *) 'ERROR..."PROD" CANNOT BE REDEFINED'
C            STOP
            KERR = .TRUE.
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'ADD ') THEN
C
C         ADDED SPECIES
C
         IF (CNTNUD .AND. (.NOT. NEC(13))) THEN
            DO 350 K = 1, KKGAS
               ADD(K) = 0.0
350         CONTINUE
         ENDIF
         NEC(13)  = .TRUE.
         CALL CKSNUM (LINE, 1, LOUT, KSYM(1,1), KKGAS, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ELSE
            ADD(KSPEC) = VALUE(1)
         ENDIF
C
C--------------SOLUTION ESTIMATE--------------------
C
      ELSEIF (KEYWRD .EQ. 'XEST') THEN
C
C              GAS PHASE MOLE FRACTION ESTIMATE
C
         CALL CKSNUM (LINE, 1, LOUT, KSYM(1,1), KKGAS, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ELSE
            LESTIM(1,1) = .TRUE.
            XEST(KSPEC) = VALUE(1)
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'SURF') THEN
C
C             SURFACE SITE FRACTION ESTIMATE
C
         KSTOT = 0
         IFOUND = 0
         DO 320 IM = 1, NMAT
            CALL SKSNUM (LINE, -1, LOUT, KSYM(1,IM), KK(IM), 
     1                   PSYM(1,IM), NPHASE(IM), KKPHAS(1,IM),
     2                   KDUM, NT1, NVAL, VALUE, IERR)
            IF (.NOT. IERR .AND. KDUM .NE. 0) THEN
               IFOUND = 1
               IMAT = IM
               IPHASE = PSRPID (KDUM, KLAST(1,IM), NPHASE(IM), LOUT)
               IF (IPHASE .LE. 0) THEN
                  KERR = .TRUE.
                  RETURN
               ENDIF
               LESTIM(IPHASE,IM) = .TRUE.
               KSPEC = KDUM + KSTOT
               XEST(KSPEC) = VALUE(1)
            ENDIF
            NVSUM = NVSUM + NVAL
            KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
320      CONTINUE
         IF (IFOUND .EQ. 0) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'BULK') THEN
C
C             BULK PHASE MOLE FRACTION ESTIMATE
C
         KSTOT = 0
         IFOUND = 0
         DO 330 IM = 1, NMAT
            CALL SKSNUM (LINE, -1, LOUT, KSYM(1,IM), KK(IM), 
     1                   PSYM(1,IM), NPHASE(IM), KKPHAS(1,IM),
     2                   KDUM, NT1, NVAL, VALUE, IERR)
            IF (.NOT. IERR .AND. KDUM .NE. 0) THEN
               IFOUND = 1
               IMAT = IM
               IPHASE = PSRPID (KDUM, KLAST(1,IM), NPHASE(IM), LOUT)
               IF (IPHASE .LE. 0) THEN
                  KERR = .TRUE.
                  RETURN
               ENDIF
               LESTIM(IPHASE,IM) = .TRUE.
               KSPEC = KDUM + KSTOT
               XEST(KSPEC) = VALUE(1)
            ENDIF
            NVSUM = NVSUM + NVAL
            KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
330      CONTINUE
         IF (IFOUND .EQ. 0) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'SDEN') THEN
C
C            SURFACE SITE DENSITY CHANGES (FROM STANDARD STATE)
C(note this serves to either change the variable or the constant)
C
         IFOUND = 0
         DO 340 IM = 1, NMAT
            CALL CKSNUM(LINE, 1 ,LOUT, PSYM(1,IM), NPHASE(IM), IPHASE, 
     1                  NVAL,VALUE, IERR)
            IF (.NOT. IERR .AND. IPHASE .GT. 0) THEN
               SDEN(IPHASE,IM) = VALUE(1)
               IFOUND = 1
            ENDIF
340      CONTINUE
         IF (IFOUND .EQ. 0) THEN
           WRITE (LOUT,'(A)')
     1       ' ERROR READING DATA FOR KEYWORD '//KEYWRD
           KERR = .TRUE.
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'AFRA') THEN
C
C           FRACTIONAL PART OF AREA FOR SURFACE MATERIAL
C
         CALL CKSNUM(LINE, -2 ,LOUT, MSYM, NMAT, IMAT, NVAL,
     1              VALUE, IERR)
         IF (NVAL.EQ.1) THEN
            DO 345 IPSR = 1, MAXPSR
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  AFRAC(IMAT,IPSR) = VALUE(1)
               ENDIF
 345        CONTINUE
            IERR = .FALSE.
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               IF (IMAT .LT. 0 .OR. IMAT .GT. NMAT) THEN
                  WRITE (LOUT,'(A)')
     1              ' ERROR...VALID SURFACE MATERIAL NAME NOT FOUND'
                  KERR = .TRUE.
               ELSE
                  AFRAC(IMAT,IPSR) = VALUE(1)
               ENDIF
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ELSEIF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'SFAC') THEN
C
C        MULTIPLYING FACTOR FOR SURFACE PRODUCTION RATES
C
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 352 IPSR = 1, MAXPSR
               SFAC(IPSR) = VALUE(1)
 352        CONTINUE
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               SFAC(IPSR) = VALUE(1)
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
C
      ELSEIF (KEYWRD .EQ. 'GFAC') THEN
C
C        MULTIPLYING FACTOR FOR GAS-PHASE PRODUCTION RATES
C
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 355 IPSR = 1, MAXPSR
               GFAC(IPSR) = VALUE(1)
 355        CONTINUE
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               GFAC(IPSR) = VALUE(1)
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'BOHM') THEN
C
C        TRANSPORT-LIMITATION FACTOR FOR BOHM-LIMITED FLUXES
C        FOR BHMXI > 0.0, ESTABLISHES USE OF BOHM CRITERION
C        WITH CORRECTIONS FOR NEGATIVE ION PRESENCE.
C
         CALL CKXNUM (LINE, -2, LOUT, NVAL, VALUE, IERR)
         IF (NVAL .EQ. 1) THEN
            DO 357 IPSR = 1, MAXPSR
               BHMXI(IPSR) = VALUE(1)
 357        CONTINUE
         ELSEIF (NVAL .EQ. 2 .OR. NVAL .EQ. -2) THEN
            IPSR = INT(VALUE(2))
            IERR = .FALSE.
            IF (IPSR .LE. MAXPSR .AND. IPSR .GT. 0) THEN
               BHMXI(IPSR) = VALUE(1)
            ELSE
               WRITE (LOUT, 6002) IPSR, MAXPSR
               IERR = .TRUE.
            ENDIF
         ENDIF
         KERR = KERR.OR.IERR
C
C
C--------------SENSITIVITY KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'ASEN') THEN
C
C        ALL REACTION SENSITIVITY
C
         LSEN = .TRUE.
         LSENT = .TRUE.
         LSENG = .TRUE.
         DO 400 K = 1, KKTOT
            ISEN(K) = .TRUE.
400      CONTINUE
C
      ELSEIF (KEYWRD .EQ. 'SEN ') THEN
C
C         SPECIFIC SPECIES KEYWORDS
C
         LSEN = .TRUE.
         IERR = .FALSE.
         KKALL = KKMAX*NMAT
         CALL CKCRAY (LINE, KKALL, KSYM, LOUT, KKALL,
     1                KSP, NS, IERR)
         IF (.NOT. IERR .AND. NS .GT. 0) THEN
            DO 410 N = 1, NS
               IFOUND = 0
               KSPEC = 0
               IF (KSP(N) .LE. KKGAS) THEN
                  KSPEC = KSP(N)
                  IFOUND = 1
               ELSE
                  KSTOT = 0
                  DO 450 IM = 1, NMAT
                     IF (IFOUND.EQ.0 .AND. 
     1                   KSP(N) .LE. (IM-1)*KKMAX+KK(IM)) THEN
                        KSPEC = KSP(N)-(IM-1)*KKMAX + KSTOT
                        IFOUND = 1
                     ENDIF
                     KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
 450              CONTINUE
               ENDIF
               IF (IFOUND.NE.0.AND.KSPEC.GT.0) THEN
                  ISEN(KSPEC) = .TRUE.
               ELSE
                  IERR = .TRUE.
               ENDIF
 410        CONTINUE
         ENDIF
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'SFLR') THEN
C
C        FLOOR VALUE FOR THE SPECIES BOUNDS
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SFLR, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'SPOS') THEN
C
C        CONVERT NEGATIVE GAS AND SITE SPECIES SOLUTIONS
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SPOS, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'SENT') THEN
C
C         TEMPERATURE SENSITIVITY KEYWORD
C
         LSEN  = .TRUE.
         LSENT = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'SENG') THEN
C
C         GROWTH RATE SENSITIVITY KEYWORD
C
         LSEN  = .TRUE.
         LSENG = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'EPSS') THEN
C
C         SENSITIVITY PRINT OPTION
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, EPSS, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'EPST') THEN
C
C         SENSITIVITY PRINT OPTION (TEMPERATURE)
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, EPST, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'EPSG') THEN
C
C         SENSITIVITY PRINT OPTION (GROWTH RATE)
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, EPSG, IERR)
         KERR = KERR.OR.IERR
C
C--------------RATE-OF-PRODUCTION KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'AROP') THEN
C
C        ALL SPECIES RATE-OF-PRODUCTION
C
         LROP = .TRUE.
         DO 500 K = 1, KKTOT
            IROP(K) = .TRUE.
500      CONTINUE
C
      ELSEIF (KEYWRD .EQ. 'ROP ') THEN
C
C         SPECIFIC SPECIES KEYWORDS
C
         LROP = .TRUE.
         IERR = .FALSE.
         KKALL = KKMAX*NMAT
         CALL CKCRAY (LINE, KKALL, KSYM, LOUT, KKALL,
     1                KSP, NS, IERR)
         IF (.NOT. IERR .AND. NS .GT. 0) THEN
            DO 510 N = 1, NS
               IFOUND = 0
               KSPEC = 0
               IF (KSP(N) .LE. KKGAS) THEN
                  KSPEC = KSP(N)
                  IFOUND = 1
               ELSE
                  KSTOT = 0
                  DO 550 IM = 1, NMAT
                     IF (IFOUND.EQ.0 .AND. 
     1                   KSP(N) .LE. (IM-1)*KKMAX+KK(IM)) THEN
                        KSPEC = KSP(N)-(IM-1)*KKMAX + KSTOT
                        IFOUND = 1
                     ENDIF
                     KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
 550              CONTINUE
               ENDIF
               IF (IFOUND.NE.0.AND.KSPEC.GT.0) THEN
                  IROP(KSPEC) = .TRUE.
               ELSE
                  IERR = .TRUE.
               ENDIF
 510        CONTINUE
         ENDIF
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'EPSR') THEN
C
C         RATE-OF-PRODUCTION PRINT OPTION
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, EPSR, IERR)
         KERR = KERR.OR.IERR
C
C--------------PRINTING AND RESTARTING KEYWORDS-------------------
C
      ELSEIF (KEYWRD .EQ. 'PRNT') THEN
C
C          PRINT CONTROL
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         IPRNT    = INT(VALUE(1))
C
      ELSEIF (KEYWRD .EQ. 'TERR') THEN
C
C         AN OUTPUT FILE FOR THE TERRAIN TOPO SIMULATOR WILL BE WRITTEN
C
         LTERRN = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'RSTR') THEN
C
C          RESTART CHECK
C
         DO 560 IM = 1, NMAT
            DO 560 I = 1 , NPHASE(IM)
               LESTIM(I,IM) = .TRUE.
560      CONTINUE
         LRSTRT = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'CNTN') THEN
C
C          CONTINUATION FLAG
C
         LCNTUE   = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'END ') THEN
C
C        LAST CARD
C
         GO TO 1100
      ELSE
C
C--------------END OF KEYWORDS--------------------
C
C
C        TO GET HERE, AN INVALID KEYWORD WAS READ
C
         WRITE (LOUT, *) ' ERROR...ILLEGAL KEYWORD'
         KERR = .TRUE.
      ENDIF
      GO TO 1000
 1100 CONTINUE
C
C          CHECK FOR NECESSARY INPUT
C
      IF (.NOT. CNTNUD) THEN
         IF (.NOT. NEC(1) ) THEN
            WRITE (LOUT, *)
     1      ' ERROR..."ENRG" OR "TGIV" MUST BE PROVIDED '
            KERR = .TRUE.
C
         ELSEIF (.NOT. NEC(2) ) THEN
            WRITE (LOUT, *)
     1      ' ERROR...NO REACTOR TEMPERATURE IS SPECIFIED OR ESTIMATED'
            KERR = .TRUE.
C
         ELSEIF (.NOT. NEC(3) ) THEN
            WRITE (LOUT, *) ' ERROR...PRESSURE NOT GIVEN'
            KERR = .TRUE.
C
         ELSEIF (.NOT. NEC(4) ) THEN
            WRITE (LOUT, *) ' ERROR..."FLRT" OR "TAU" MUST BE PROVIDED'
            KERR = .TRUE.
C
         ELSEIF (.NOT. NEC(5) ) THEN
            WRITE (LOUT, *)
     1      ' ERROR...VOLUME MUST BE SPECIFIED.'
            KERR = .TRUE.
C
         ENDIF
C
         IF (LENRGY(1) ) THEN
            IF (.NOT. NEC(6) ) THEN
              WRITE (LOUT, *)
     1        'WARNING...HEAT LOSS NEEDS TO BE SPECIFIED ',
     1               'FOR "ENRG" PROBLEMS'
              WRITE (LOUT, *)
     1        '          A HEAT LOSS OF 0.0 IS ASSUMED.'
            END IF
            IF (.NOT. NEC(7) ) THEN
               WRITE (LOUT, *)
     1        ' ERROR...T INLET MUST BE SPECIFIED FOR "ENRG" PROBLEMS'
               KERR = .TRUE.
            ENDIF
         ELSE
            IF (NEC(6) ) WRITE (LOUT, *)
     1       ' WARNING...HEAT LOSS SPECIFIED FOR "TGIV" PROBLEM'
            IF (NEC(7) ) WRITE (LOUT, *)
     1        ' WARNING..."TINL" SPECIFIED FOR "TGIV" PROBLEM'
         ENDIF
C
         IF (NEC(8) ) THEN
            DO 1200 I = 9, 12
               IF (NEC(I)) THEN
                  WRITE (LOUT, *)
     1            'ERROR...MULTIPLE DEFINITION OF INLET COMPOSITION'
                  KERR = .TRUE.
               ENDIF
1200        CONTINUE
         ELSE
            DO 1300 I = 9, 12
               IF (.NOT. NEC(I)) THEN
                  WRITE (LOUT, *)
     1            'ERROR...INLET COMPOSITION NOT DEFINED'
                  KERR = .TRUE.
               ENDIF
1300        CONTINUE
         ENDIF
         IF (KEL .NE. 0) THEN
C
            IF (.NOT. NEC(14)) THEN
               WRITE (LOUT,*) 
     1     ' ERROR..."ENGE" OR "TEGV" MUST BE SPECIFIED WITH ELECTRONS'
               KERR = .TRUE.
            ENDIF
C
            IF (.NOT. NEC(15)) THEN
               IF (NNU .LT. KKGAS-1) THEN
                  WRITE (LOUT,*)
     1     ' ERROR...NEED TO SPECIFY "XSEK" FOR ALL SPECIES OR "XSDF"'
                  KERR = .TRUE.
               ENDIF
            ELSE
               DO 1310 K = 1, KKGAS
                  IF (QEK(K) .LE. 0.0 .AND. K .NE. KEL) THEN
                     QEK(K) = QDEF
                  ENDIF
 1310          CONTINUE
            ENDIF
C
            IF (.NOT. NEC(16)) THEN
               WRITE (LOUT,*)
     1         ' WARNING..."ETMP" NOT SPECIFIED; USING TE=T'
               DO 1320 IPSR = 1, MAXPSR
                  TE(IPSR) = T(IPSR)
 1320          CONTINUE
            ENDIF            
C
            IF (.NOT. NEC(17)) THEN
               WRITE (LOUT,*)
     1         ' WARNING..."TEIN" NOT SPECIFIED; USING TEIN=TIN'
               DO 1330 IPSR = 1, MAXPSR
                  TEIN(IPSR) = TIN(IPSR)
 1330          CONTINUE
            ENDIF            
         ELSE
            DO 1335 IPSR = 1, MAXPSR
               TE(IPSR) = T(IPSR)
               TEIN(IPSR) = TIN(IPSR)
 1335       CONTINUE
         ENDIF
      ELSEIF (ENRGIN) THEN
         IF (.NOT. NEC(6) ) THEN
            WRITE (LOUT, *)
     1     ' ERROR...HEAT LOSS MUST BE SPECIFIED FOR "ENRG" PROBLEMS'
            KERR = .TRUE.
         ELSEIF (.NOT. NEC(7) ) THEN
            WRITE (LOUT, *)
     1        ' ERROR...T INLET MUST BE SPECIFIED FOR "ENRG" PROBLEMS'
            KERR = .TRUE.
         ENDIF
C
      ENDIF
C
C  CHECK FOR CONSISTENT INPUT FOR RF SHEATH MODEL
C
      DO 1342 IPSR = 1, NUMPSR
         DO 1340 IM = 1, NMAT
            IF (LRFSH(IM,IPSR)) THEN
               IF (LRFCUR) THEN
                  IF (RFVAC .NE. 0.0) THEN
                     WRITE (LOUT,*)
     1             ' ERROR...VOLTAGE SPECIFIED WITH RF CURRENT CONTROL'
                     KERR = .TRUE.
                  ELSEIF (RFIAC .EQ. 0.0) THEN
                     WRITE (LOUT,*)
     1          ' ERROR...NO CURRENT AMPLITUDE SPECIFIED FOR RF SHEATH'
                     KERR = .TRUE.
                  ELSE
                     RFAPAR = RFIAC
                     RFDPAR = 0.0
                  ENDIF
               ELSE
                  IF (RFIAC .NE. 0.0) THEN
                     WRITE (LOUT,*)
     1           ' ERROR...CURRENT SPECIFIED WITH RF VOLTAGE CONTROL'
                     KERR = .TRUE.
                  ELSEIF (RFVDC .GE. -ABS(RFVAC)) THEN
                     WRITE (LOUT,*)
     1         ' ERROR...VDC > -ABS(VAC) WILL RESULT IN RF MODEL ERROR'
                     KERR = .TRUE.
                  ELSE
                     RFAPAR = RFVAC
                     RFDPAR = RFVDC
                  ENDIF
               ENDIF
            ENDIF
 1340    CONTINUE
 1342 CONTINUE
C
C          FIND THE INITIAL SURFACE TEMPERATURE, ION TEMPERATURE
C
      DO 1350 IPSR = 1, NUMPSR
         DO 1345 IM = 1, NMAT
            IF (.NOT. LTDIFF(IM,IPSR)) THEN
               TSURF(IM,IPSR) = T(IPSR)
            ENDIF
 1345    CONTINUE
         IF (.NOT. LTION(IPSR)) TIONP(IPSR) = T(IPSR)
1350  CONTINUE
C
C      NORMALIZE THE SOLUTION ESTIMATE
C
C            CHECK FOR BULK PHASES TO BE ETCHED - THEY MUST HAVE
C            INITIAL CONDITIONS SPECIFIED.
C            (currently it only makes sense to check the first psr)
C
      KSTOT = 0
      DO 1570 IM = 1, NMAT
         DO 1550 I = 1 , NPHASE(IM)
            IF (LESTIM(I,IM)) THEN
               SUM = 0.0
               DO 1400 KM = KFIRST(I,IM), KLAST(I,IM)
                  K = KM + KSTOT
                  SUM = SUM + XEST(K)
1400           CONTINUE
               IF (SUM .GT. 0.0) THEN
                  DO 1500 KM = KFIRST(I,IM), KLAST(I,IM)
                     K = KM + KSTOT
                     XEST(K) = XEST(K) / SUM
1500              CONTINUE
                  IF (ABS(SUM-1.0) .GT. 1.E-6) WRITE (LOUT, *)
     1            ' CAUTION...XEST MOLE FRACTIONS SUM TO ', SUM,
     1            ' FOR PHASE', I,' = ',PSYM(I,IM)
               ELSEIF (SUM .EQ. 0.0) THEN
                  DO 1520 KM = KFIRST(I,IM), KLAST(I,IM)
                     K = KM + KSTOT
                     XEST(K) = 1.0/MAX(KKPHAS(I,IM),1)
 1520             CONTINUE
                  WRITE (LOUT, *) 
     1            ' CAUTION...XEST MOLE FRACTIONS SUM TO ', SUM,
     1            ' FOR PHASE', I,' = ',PSYM(I,IM)
               ENDIF
            ELSE
               IF (ETCH(I,IM,1) .AND. .NOT. CNTNUD) THEN
                  IF (KLAST(I,IM) .EQ. KFIRST(I,IM) ) THEN
                     K = KFIRST(I,IM) + KSTOT
                     XEST(K) = 1.0
                  ELSE
                     WRITE (LOUT, *)
     1               'ERROR: ETCH KEYWORD WAS SUPPLIED FOR PHASE',
     1               I, ' = ', PSYM(I,IM),
     1               '.  HOWEVER, NO BULK PHASE COMPOSITIONS WERE GIVEN'
                     KERR = .TRUE.
                  ENDIF
               ENDIF
            ENDIF
1550     CONTINUE
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
1570  CONTINUE
C
C          CHECK THE SUM OF USER SUPPLIED REACTANTS
C
      IF (.NOT. LEQUIV) THEN
         SUM = 0.0
         DO 1600 K = 1, KKGAS
            SUM = SUM + XIN(K,1)
1600     CONTINUE
         DO 1700 K = 1, KKGAS
            XIN(K,1) = XIN(K,1) / SUM
1700     CONTINUE
C
         IF (ABS(SUM-1.0) .GT. 1.0E-6) WRITE (LOUT, *)
     1      ' CAUTION...REACTANT MOLE FRACTIONS SUM TO ', SUM
      ENDIF
C
C           STOP IF ANY ERRORS ENCOUNTERED IN INPUT SO FAR
C
      IF (KERR) THEN
         WRITE (LOUT, *)'PSR FATAL ERROR: error in input file',
     1                ' - stopping'
C         STOP
         RETURN
      ENDIF
C
C         CHECK THE SUMS FOR FUEL, OXIDIZER AND ADDED SPECIES
C
      IF ((.NOT.CNTNUD) .AND. LEQUIV) THEN
         SUMF = 0.0
         SUMO = 0.0
         SUMA = 0.0
         DO 1800 K = 1, KKGAS
            SUMF = SUMF + FUEL(K)
            SUMO = SUMO + OXID(K)
            SUMA = SUMA + ADD(K)
1800     CONTINUE
C
C         NORMALIZE FUEL AND OXIDIZER FRACTIONS
C
         DO 1900 K = 1, KKGAS
            FUEL(K) = FUEL(K) / SUMF
            OXID(K) = OXID(K) / SUMO
1900     CONTINUE
C
         IF (ABS(SUMF-1.0) .GT. 1.E-6) WRITE (LOUT, *)
     1             ' CAUTION...FUEL FRACTIONS SUM TO ', SUMF
         IF (ABS(SUMO-1.0) .GT. 1.E-6) WRITE (LOUT, *)
     1             ' CAUTION...OXIDIZER FRACTIONS SUM TO ',  SUMO
         IF (SUMA .GT. 1.0) THEN
            WRITE (LOUT, *)
     1             ' ERROR...ADDED SPECIES SUM TO ',  SUMA
C            STOP
            KERR = .TRUE.
            RETURN
         ENDIF
      ENDIF
C
C           BALANCE THE EQUATION FOR COMPLETE COMBUSTION
C
      IF (LEQUIV) THEN
         IF (.NOT. CNTNUD) THEN
C
C  Determine how many elements from mechanism are actually used
C  in this problem:
C
            MMUSED = 0
            DO 1920 M = 1, MM
               MUSD(M) = 0
 1920       CONTINUE
C
            DO 1950 M = 1, MM
               IFOUND = 0
               DO 1940 K = 1, KKGAS
                  IF (IFOUND.EQ.0.AND.NCF(M,K,1).NE.0) THEN
                     IF (FUEL(K).NE.0.0 .OR. OXID(K).NE.0.0) THEN
                        MMUSED = MMUSED + 1
                        MUSD(MMUSED) = M
                        IFOUND = 1
                     ENDIF
                  ENDIF
 1940          CONTINUE
 1950       CONTINUE
C
C  Create equilbrium solution matrix only from elements used:
C
            DO 2000 MU = 1, MMUSED
C
               M = MUSD(MU)
               RSUM = 0.0
               ASUM = 0.0
               DO 2100 K = 1, KKGAS
                  RSUM = RSUM + FUEL(K) * NCF(M,K,1)
                  ASUM = ASUM + OXID(K) * NCF(M,K,1)
2100           CONTINUE
               RHS(MU) = RSUM
               A(MU,1) = -ASUM
C
               N = 1
               DO 2300 L = 1, NPROD
                  N = N + 1
                  A(MU,N) = NCF(M,KPROD(L),1)
2300           CONTINUE
2000        CONTINUE
C
C Give user error message if balance matrix mm rows & n cols not equal
C
            IF (N .NE. MMUSED) THEN
               WRITE (LOUT, *)
     1         ' ERROR in PSRKEY:  In equilibrium calculation, ',
     2         ' the number of elements must be the same in reactants',
     3         ' and products.'
C     1         ' THE 
               KERR = .TRUE.
               RETURN
            ENDIF
C
C            CALL LINPAC TO BALANCE THE REACTIONS
C
C*****precision > double
            CALL DGEFA (A, NATJ, MMUSED, IPVT, INFO)
            IF (INFO .NE. 0) THEN
               WRITE (LOUT, *)
     1'ERROR..INCONSISTENT DEFINITION OF FUEL, OXIDIZER, AND PRODUCTS'
C               STOP
               KERR = .TRUE.
               RETURN
            ENDIF
            CALL DGESL (A, NATJ, MMUSED, IPVT, RHS, 0)
C*****END precision > double
C
C*****precision > single
C            CALL SGEFA (A, NATJ, MMUSED, IPVT, INFO)
C            IF (INFO .NE. 0) THEN
C               WRITE (LOUT, *)
C     1'ERROR..INCONSISTENT DEFINITION OF FUEL, OXIDIZER, AND PRODUCTS'
CC               STOP
C               KERR = .TRUE.
C               RETURN
C            ENDIF
C            CALL SGESL (A, NATJ, MMUSED, IPVT, RHS, 0)
C*****END precision > single
C
            SUM = 0.0
            DO 2500 K = 1, KKGAS
               STOICH(K) = FUEL(K) + RHS(1) * OXID(K)
               SUM = SUM + STOICH(K)
2500        CONTINUE
            DO 2600 K = 1, KKGAS
               STOICH(K) = STOICH(K) / SUM
2600        CONTINUE
         ENDIF
C
         SUM = 0.0
         SUMA = 0.0
         DO 2700 K = 1, KKGAS
            IF (OXID(K) .GT. 0.0) THEN
               XIN(K,1) = STOICH(K) / EQUIV
            ELSE
               XIN(K,1) = STOICH(K)
            ENDIF
            SUM = SUM + XIN(K,1)
            SUMA = SUMA + ADD(K)
2700     CONTINUE
C
         DO 2800 K = 1, KKGAS
            XIN(K,1) = XIN(K,1) / SUM
2800     CONTINUE
C
         IF (SUMA .GT. 0.0) THEN
            SUM = 1.0 + SUMA
            DO 3000 K = 1, KKGAS
               XIN(K,1) = XIN(K,1) / SUM + ADD(K)
3000        CONTINUE
         ENDIF
      ENDIF
C
C-------FIND THE INITIAL FLOW RATE IF SCCM KEYWORD IS USED
C
      IF (LFLRT(1)) THEN
         IF (LSCCM) THEN
            TEMP(1) = 298.15
            TEMP(2) = 298.15
            TEMP(3) = 298.15
            CALL CKRHOX (PATM, TEMP, XIN(1,1), ICKWRK, RCKWRK, RHOIN)
            FLRT(1) = SCCMIN * RHOIN / 60.0
         ENDIF
      ENDIF
C
C-------CALCULATION OF EQUILIBRIUM VALUES (FOR INITIAL GUESS)
C
      IF (.NOT. LESTIM(1,1) .AND. .NOT. CNTNUD) THEN
         NOP = 1
         MON = 0
         PQ(1) = P(1)
         TQ(1) = T(1)
C
         NCON = 0
         LEQST = .FALSE.
         LEQPRN = .FALSE.
         LEQCNT = .FALSE.
         KKPEQ(1) = KKGAS
         KFPEQ(1) = 1
         KLPEQ(1) = KKGAS
         NNPEQ = 1
         NNBEQ = 0
         NFBEQ = 0
         NLBEQ = 0
         LENCEQ = 500
         KKSEQ = 0
         KKBEQ = 0
         LSEQ = 0
C        starting temperature (K)
         VINPUT(1) = T(1)
C        starting pressure (atm)
         VINPUT(2) = PA(1)
C
         LSURF = 0
         CALL EQUIL (LOUT, LEQPRN, LSEQ, LEQST, LEQCNT, ICKWRK,
     1               RCKWRK, ISKWRK, RSKWRK, KKPEQ, KFPEQ, KLPEQ,
     2               NNPEQ, NNBEQ, NFBEQ, NLBEQ, ACTEQ, PDENEQ,
     3               LENIEQ, IEQWRK, LENEQ, EQWRK, LENCEQ, CEQWRK,
     4               MM, KKGAS, KKSEQ, KKBEQ, KKGAS, ATOM, KSYM(1,1),
     5               NOP, MON, .FALSE., XIN(1,1), TQ(1), PQ(1), NCON,
     6               KCON, XCON, VINPUT, .FALSE., .FALSE.,
     7               IRGWRK, RRGWRK, LSURF, KERR)
         IF (KERR) RETURN
         CALL EQSOL (KKGAS, EQWRK, XEST, RHS, TQ(2), PQ(2),
     1               HQ(2), VQ(2), SQ(2), WM(2), CQ(2), CDETQ,
     2               KFPEQ, KLPEQ, NNBEQ, NFBEQ, NLBEQ)
C
         WRITE (LOUT, '(/A/)')
     1   ' FIRST SOLUTION ESTIMATE IS EQUILIBRIUM '
C
      ENDIF
C
C    FIND INITIAL GUESSES FOR ALL OTHER PHASES IF NONE WERE PROVIDED
C
      KSTOT = 0
      DO 4200 IM = 1, NMAT
         IF (NPHASE(IM) .GT. 1) THEN
            DO 4000 I = 2, NPHASE(IM)
               IF (.NOT. LESTIM(I,IM) .AND. .NOT. CNTNUD) THEN
                  DO 3990 KM = KFIRST(I,IM), KLAST(I,IM)
                     K = KM + KSTOT
                     XEST(K) = 1.0/KKPHAS(I,IM)
3990              CONTINUE
               ENDIF
4000        CONTINUE
         ENDIF
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
4200  CONTINUE
C
C           MAKE SURE ELECTRON ENERGY LOSS VALUES ARE IN ASCENDING 
C           TEMPERATURE ORDER
C
      DO 4720 N = 2, NQLSE
         IF (TQLSE(N-1) .GE. TQLSE(N)) THEN
            WRITE (LOUT,*)
     1              ' ERROR...E-ENERGY LOSS VALUES ARE OUT OF ORDER'
            KERR = .TRUE.
         ENDIF
 4720 CONTINUE
C
6002  FORMAT ('PSRKEY Error: IPSR, ',I3,' is greater than maxpsr, ',
     1         I3,', or less than or equal to zero')
7000  FORMAT (A4, A)
8000  FORMAT (10X, A4, A76)
C
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE PSRSEN (LOUT, LSAVE, LRECOV, KK, II, IISUR, NSPHCH,
     1                   KKGAS, KKSURF, KKBULK, NPTS, NATJ, NPHASE,  
     2                   NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK, 
     3                   KKPHAS, KFIRST, KLAST, KOCC, KSYM, PSYM, 
     4                   CCWORK, CSWORK, ETCH, LENRGY, LFLRT, LTDIFF, 
     5                   ABSOL, RELAT, DT, ABOVE, BELOW, WT, HIN, YIN,
     6                   X, DEN, TIN, P, TAU, FLRT, V, AREA, Q, HTRN,
     7                   TAMBNT, LHTRN, T, TSURF, S, SN, F, FN, A,
     8                   ACT, SDEN, SDEN0, SITDOT, SDOT, SDOTD, DFDALP,
     9                   GRATE, GRATED, DSKDPL, LSEN, LSENT, LSENG,
     *                   ISEN, EPSS, EPST, EPSG, ICKWRK, RCKWRK, 
     1                   ISKWRK, RSKWRK, MAPPH, IPVT, NIK, SCRTCH, NT,
     3                   NTE, NYS, NY, NSDEN, NSDENS, KEL, KKION, 
     3                   KCHG, KION, QEK, POWR, LENRGE, TEIN, TE, ROP,
     4                   ROPS, IEIMP, EIONSH, NQLSE, QLSE, TQLSE, GFAC,
     5                   SFAC, TIONP, LCONFN, IEXC, REXC, SDOTI, 
     6                   SITDTI, KKTOT, NMAT, IISUR0, IISMAX, KKMAX,
     7                   NPHMAX, IMISK, IMRSK, IMCSK, AFRAC, MSYM,
     8                   ITDEP, BHMXI, FLXION, FLXCOR, SDOTT, ESHTH,
     9                   LELSH, LRFSH, LRFCUR, RFFREQ, RFAPAR,
     *                   RFDPAR, LWALHB, VISC, THCOND, LTERRN, ANISOT,
     1                   NBHM, KERR, LOUTSH, BOLTZ, AVOG, EV2K, BIASPW,
     2                   HIM)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KKPHAS(NPHMAX,NMAT), KFIRST(NPHMAX,NMAT), 
     1          KLAST(NPHMAX,NMAT), KOCC(KKMAX,NMAT), ICKWRK(*),
     1          ISKWRK(*), MAPPH(NPHMAX,NMAT), IPVT(NATJ), 
     2          WT(KKMAX,NMAT), HIN(KKMAX), YIN(KKMAX), X(KKTOT),
     3          DEN(KKMAX,NMAT), SN(NATJ), S(NATJ), ABOVE(NATJ),
     3          BELOW(NATJ), F(NATJ), FN(NATJ), SCRTCH(NIK,6),
     4          RCKWRK(*), RSKWRK(*), A(NATJ,NATJ), ACT(KKMAX,NMAT),
     5          SDEN(NPHMAX,NMAT), SDEN0(NPHMAX,NMAT),
     6          SITDOT(NPHMAX,NMAT), SDOT(KKTOT), SDOTD(NIK),
     7          DFDALP(NIK), GRATE(NPHMAX,NMAT), GRATED(NPHMAX,NMAT), 
     8          DSKDPL(KKTOT,NATJ), KEEP(3), KEEPS(3), KCHG(KKMAX),  
     9          KION(*), QEK(KKGAS), ROP(*), ROPS(IISUR0,NMAT), 
     *          IEIMP(*), TEMP(3), QLSE(NPTS), TQLSE(NPTS), IEXC(*),
     *          REXC(*), KK(NMAT), SITDTI(NPHMAX,NMAT), KKSURF(NMAT),
     1          KKBULK(NMAT), NPHASE(NMAT), NNSURF(NMAT), NFSURF(NMAT),
     2          NLSURF(NMAT), NNBULK(NMAT), NFBULK(NMAT), NLBULK(NMAT),
     4          IISUR(NMAT), NSPHCH(NMAT), IMISK(NMAT), IMRSK(NMAT),
     4          IMCSK(NMAT), AFRAC(NMAT), ITDEP(*),EIONSH(NMAT), 
     5          SDOTT(KKTOT), SDOTI(KKMAX,IISUR0,NMAT), FLXION(*),
     6          FLXCOR(*), ESHTH(NMAT), TSURF(NMAT), KEEPIM(3),
     6          NBHM(NMAT), BIASPW(NMAT), HIM(KKMAX,NMAT)
      CHARACTER*16  KSYM(KKMAX,NMAT), PSYM(NPHMAX,NMAT), CCWORK(*), 
     1              CSWORK(*), MSYM(NMAT)
      CHARACTER*18  DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO
      CHARACTER*80  IHOL
      LOGICAL  LENRGY, LFLRT, LTDIFF(NMAT), LTIME, LSEN, LSENT, LSENG,
     1         ISEN(KKTOT), ETCH(NPHMAX,NMAT), LHTRN, LWRITE, IERR, 
     2         LENRGE, LCONFN, LELSH(NMAT), LRFSH(NMAT), LRFCUR, 
     3         LWALHB(NMAT), LTERRN
      LOGICAL KERR
      PARAMETER (SMALL = 1.E-25)
      PARAMETER(IDIMRF = 101, ZEROE = 0.0, LPMAX=60)
C
C  CHARACTER DATA FILE COMMON BLOCK
C
      COMMON /CDATAF/ VERSNN, DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C INPUT
C
C   ABSOL  - ABSOLUTE PERTURBATION FOR EVALUATING THE JACOBIAN AND DF/DA
C   DT     - SIZE OF THE TIME STEPS.
C              CGS UNITS - SEC
C   EPSS   - THRESHOLD VALUE FOR SPECIES SENSITIVITY COEFFICIENTS.
C   EPST   - THRESHOLD VALUE FOR TEMPERATURE SENSITIVITY COEFFICIENTS.
C   EPSG   - THRESHOLD VALUE FOR GROWTH RATE SENSITIVITY COEFF.
C   FLRT   - THE MASS FLOW RATE.
C              CGS UNITS - GM/SEC
C   LTIME  - LTIME=.FALSE. CAUSES FUNCTION TO SET UP THE ALGEBRAIC
C            EQUATIONS.
C   HIN    - ARRAY OF SPECIES ENTHALPIES AT INLET TEMPERATURE TIN.
C              CGS UNITS - ERGS/G
C              DIMENSION H(*) AT LEAST KK.
C   II     - NUMBER OF ELEMENTARY GAS PHASE CHEMICAL REACTIONS.
C   IISUR  - NUMBER OF ELEMENTARY SURFACE PHASE CHEMICAL REACTIONS.
C   ISEN   - LOGICAL ARRAY. IF ISEN(K)=.TRUE. THEN THE FIRST ORDER
C            SENSITIVITY COEFFICIENTS OF SPECIES K WITH RESPECT TO
C            RATE CONSTANTS WILL BE PRINTED OUT.
C              DIMENSION ISEN(*) AT LEAST KK.
C   KK     - NUMBER OF CHEMICAL SPECIES.
C   KSYM   - CHEMKIN SPECIES NAMES.
C              DIMENSION KSYM(*) AT LEAST KK.
C   PSYM   - SURFACE CHEMKIN PHASE NAMES
C   LENRGY - IF LENRGY=.TRUE. THEN THE ENERGY EQUATION IS SOLVED.
C            IF LENRGY=.FALSE. THEN A SPECIFIED TEMPERATURE IS USED.
C   LFLRT  - IF LFLRT=.TRUE. THEN THE MASS FLOW RATE IS SPECIFIED AND
C              THE RESIDENCE TIME HAS TO BE CALCULATED.
C            IF LFLRT=.FALSE. THEN THE RESIDENCE TIME IS SPECIFIED AND
C              THE MASS FLOW RATE HAS TO BE CALCULATED.
C
C   LOUT   - UNIT FOR PRINTED OUTPUT.
C   LRECOV - UNIT FOR RECOVER FILE
C   LSAVE  - UNIT FOR SAVE FILE.
C   LSEN   - IF LSEN=.TRUE. THEN A FIRST ORDER SENSITIVITY ANALYSIS
C            IS CARRIED OUT.
C   LSENT  - IF LSENT=.TRUE. THEN A FIRST ORDER SENSITIVITIES OF THE
C            TEMPERATURE ARE PRINTED.
C   LSENG  - IF LSENG=.TRUE. THEN A FIRST ORDER SENSITIVITIES OF THE
C            GROWTH RATE WITH RESPECT TO REACTION RATE PREEXPONENTIALS
C            ARE PRINTED.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES. NATJ=KK+2.
C   P      - THE PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   PERTRB - PERTURBATON FOR THE JACOBIAN EVALUATION
C   Q      - THE HEAT LOSS OF THE REACTOR
C              CGS UNITS - ERGS/SEC
C   RELAT  - RELATIVE PERTURBATION FOR EVALUATING THE JACOBIAN AND DF/DA
C   S      - DEPENDENT VARIABLE ARRAY. THE TEMPERATURE IS STORED IN
C            T=S(NT), AND THE MASS FRACTIONS ARE IN Y(K)=S(NYS+K)
C              DIMENSION S(*) AT LEAST NATJ.
C   T      - THE USER-SPECIFIED TEMPERATURE.
C              CGS UNITS - K
C   TAU    - THE NOMINAL RESIDENCE TIME OF THE REACTOR
C              CGS UNITS - SEC
C   TIN    - THE INLET TEMPERATURE.
C              CGS UNITS - K
C   V      - THE VOLUME OF THE REACTOR
C              CGS UNITS - CM**3
C   WT     - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
C              CGS UNITS - GM/MOLE
C              DIMENSION WT(*) AT LEAST KK.
C   YIN    - ARRAY OF SPECIES INPUT MASS FRACTIONS.
C              DIMENSION REAC(*) AT LEAST KK.
C   GRATE  - GROWTH RATE OF BULK PHASES (CM/SEC)
C              DIMENSION GRATE(*) AT LEAST NPHASE.
C   GRATED - SENSITIVITY COEFFICIENT OF GROWTH RATE WITH RESPECT TO
C            PREEXPONENTIAL OF REACTION RATE.
C              DIMENSION GRATED(*) AT LEAST NPHASE.
C   DSKDPL - DERIVATIVE OF THE PRODUCTION RATE OF THE KTH SPECIES WITH
C            RESPECT TO THE LTH SOLUTION UNKNOWN.
C              DIMENSION DSKDPL(*) AT LEAST KK BY NATJ IN CALLING PROG.
C
C
C WORK AND SCRATCH SPACE-
C   A      - THE JACOBIAN OF THE GOVERNING EQUATIONS.
C              DIMENSION A(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION
C              AND AT LEAST NATJ FOR THE SECOND.
C   RCKWRK - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   RSKWRK - FLOATING POINT SURFACE CHEMKIN WORK SPACE.
C   F      - ARRAY OF PERTURBED FUNCTION VALUES.
C               DIMENSION F(*) AT LEAST NATJ.
C   FN     - ARRAY OF THE UNPERTURBED FUNCTION VALUES.
C               DIMENSION F(*) AT LEAST NATJ.
C   DFDALP - PARTIAL DERIVATIVE OF F WRT ALPHA
C   X      - VECTOR OF MOLE FRACTIONS, SITE DENSITIES, AND
C            BULK PHASE MOE FRACTIONS
C              DIMENSION X(*) AT LEAST KK.
C   ACT    - VECTOR OF ACTIVITY COEFFICIENTS FOR ALL PHASES.
C              DIMENSION ACT AT LEAST KK.
C   DEN    - VECTOR OF DENSITIES IN GM AND CM OF SPECIES
C              DIMENSION DEN AT LEAST NPHASE.
C   SDOT   - VECTOR OF PRODUCTION RATES OF KTH SPECIES FROM
C            SURFACE REACTIONS
C              DIMENSION SDOT AT LEAST KK.
C   SDOTD  - VECTOR OF PRODUCTION RATES OF KTH SPECIES FROM
C            SURFACE REACTION - DELTA CHANGE IN LTH SOLUTION VARIABLE
C              DIMENSION SDOTD AT LEAST KK.
C   SAVEP  - STORED PRE-EXPONENTIAL FOR THE ITH REACTION
C   SAVE   - STORED VALUE FOR THE LTH SOLUTION UNKNOWN.
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   ISKWRK - INTEGER SURFACE CHEMKIN WORK SPACE.
C   IPVT   - ARRAY OF PIVOTS.
C              DIMENSION IPVT(*) AT LEAST NATJ.
C
C   SCRTCH - SCRATCH SPACE.
C              DIMENSION(II, 5)
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      LTIME = .FALSE.
      KERR = .FALSE.
C
C          FORM THE JACOBIAN, AND AT THE SAME TIME EVALUATE THE
C               UNPERTURBED FUNCTION FN
C
      CALL PSRJAC (LOUT, KK, KKGAS, KKSURF, KKBULK, NSPHCH, NPTS,NATJ,
     1             NPHASE, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2             NLBULK, KKPHAS, KFIRST, KLAST, KOCC, ABOVE, BELOW,
     3             ETCH, LTIME, LENRGY, LFLRT, LTDIFF, ABSOL, RELAT, 
     4             DT, HIN, YIN, TIN, T, TSURF, P, TAU, FLRT, V, AREA,
     5             Q, HTRN, TAMBNT, LHTRN, ICKWRK, RCKWRK, ISKWRK,
     6             RSKWRK, MAPPH, WT, NIK, SCRTCH, SN, S, F, FN, A,
     7             ACT, SDEN, SDEN0, SITDOT, NT, NTE, NYS, NY, NSDEN, 
     8             NSDENS, KEL, KKION, KCHG, KION, QEK, POWR, LENRGE, 
     9             TEIN, TE, ROP, ROPS, IEIMP, EIONSH, II, NQLSE, 
     *             QLSE, TQLSE, GFAC, SFAC, TIONP, LCONFN, IEXC, REXC, 
     1             SDOTI, SITDTI, IISUR, KKTOT, NMAT, IISUR0, KKMAX,
     2             NPHMAX, IMISK, IMRSK, IMCSK, AFRAC, ITDEP, BHMXI,
     3             FLXION, FLXCOR, SDOTT, ESHTH, LELSH, LRFSH, LRFCUR, 
     4             RFFREQ, RFAPAR, RFDPAR, LWALHB, VISC, THCOND, SDOT,
     5             NBHM, LOUTSH, BOLTZ, AVOG, EV2K, BIASPW, HIM)
C
C          FACTOR THE JACOBIAN
C
C*****precision > double
      CALL DGEFA (A, NATJ, NATJ, IPVT, INFO)
C*****END precision > double
C
C*****precision > single
C      CALL SGEFA (A, NATJ, NATJ, IPVT, INFO)
C*****end precision > single
C
      IF (INFO .NE. 0) THEN
         WRITE (LOUT, *) 'ERROR IN FACTORING THE JACOBIAN, INFO = ',INFO
C         STOP
         KERR = .TRUE.
         RETURN
      ENDIF
C
C        COMPUTE THE BASE GROWTH RATE AND  THE DERIVATIVE OF THE
C        BULK PHASE PRODUCTION RATES WRT SOLUTION VARIABLES, DSKDPL
C
      CALL CKYTX (S(NY), ICKWRK, RCKWRK, X)
      DO 200 K = 1, KKGAS
         SDOT(K) = 0.0
200   CONTINUE
      DO 300 K = KKGAS+1, KKTOT
         X(K) = S(NYS+K)
         SDOT(K) = 0.0
300   CONTINUE
C
C
      CALL PSRACT (X, KK, KKBULK, NFBULK, NLBULK, NPHASE, KFIRST, KLAST, 
     1             NNBULK, NMAT, NPHMAX, KKSURF, KKMAX, ACT)
C
C  FILL IN ACT ARRAY AND CALCULATE THE SDEN ARRAY
C  WRITE WARNING MESSAGES ONCE FOR ZERO MOLE FRACTIONS
C
      KSTOT = 0
      IDTOT = 0
      DO 900 IM = 1, NMAT
         DO 600 K = 1, KKGAS
            ACT(K,IM) = X(K)
            IF (X(K) .EQ. 0.0E0) THEN
               WRITE (LOUT,*)'SPECIES ',KSYM(K,1),
     1                ' HAS ZERO MOLE FRACTION.',
     2                ' ITS SENS. COEFFICIENT CAN''T BE NORMALIZED'
            ENDIF
 600     CONTINUE
         IF (NNBULK(IM) .GT. 0) THEN
            IF (NNSURF(IM) .GT. 0) THEN
               IF (NSPHCH(IM) .GT. 0) THEN
                  DO 700 IDM = 1, NSPHCH(IM)
                     I = IDM + IDTOT
                     SDEN(MAPPH(IDM,IM),IM) 
     1                        = S(NSDENS+I) * SDEN0(MAPPH(IDM,IM),IM)
 700              CONTINUE
               ENDIF
               DO 750 KM = KFIRST(NFSURF(IM),IM), KLAST(NLSURF(IM),IM)
                  K = KM + KSTOT
                  ACT(KM,IM) = S(NYS+K)
                  IF (S(NYS+K) .EQ. 0.0E0) THEN
                     WRITE (LOUT,*)'SPECIES ',KSYM(KM,IM),
     1                  ' HAS ZERO MOLE FRACTION.',
     2                  ' ITS SENS. COEFFICIENT CAN''T BE NORMALIZED'
                  ENDIF
 750           CONTINUE
            ENDIF
         ENDIF
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
 900  CONTINUE
C
C               CALCULATE THE BASE BULK PRODUCTION RATE, SDOT
C
      DO 500 IM = 1, NMAT
         DO 400 N = 1, NPHASE(IM)
            SITDOT(N,IM) = 0.0
400      CONTINUE
500   CONTINUE
C
C  TEMPORARILY USE SCRTCH(K,2) ARRAY FOR ION ENERGIES; REDEFINED BELOW
C
      CALL PSSDOT (KKION, KEL, KION, KKGAS, KKTOT, KCHG, NMAT, NPHASE,
     1             IISUR, KFIRST, KLAST, NFSURF, NLSURF, NNSURF, KKSURF,
     2             NFBULK, NLBULK, NNBULK, KKBULK, IISUR0, KKMAX,
     3             NPHMAX, KK, BHMXI, WT, EIONSH, NIK, 
     3             AFRAC, RFFREQ, RFAPAR, RFDPAR, LRFSH, LRFCUR, P,
     4             ACT, SCRTCH(1,4), TIONP, S(NTE), TSURF, SDEN,
     5             ISKWRK, RSKWRK, IMISK, IMRSK, FLXCOR, FLXION,
     6             ESHTH, LELSH, SFAC, SDOTT, ROPS, SDOTI, SITDTI,
     8             SITDOT, SDOT, SHLOSS, NBHM, LOUTSH, BOLTZ, AVOG, 
     9             EV2K, BIASPW, AREA, HIM(1,1))

C
      KSTOT = 0
      IDTOT = 0
      DO 3600 IM = 1, NMAT
         IF (NNBULK(IM) .GT. 0) THEN
C
C CALCULATE DEN ARRAY (NOTE ONLY VALID FOR SURFACE AND BULK PHASE)
C
            TEMP(1) = S(NT)
            TEMP(2) = S(NTE)
            TEMP(3) = TIONP
            CALL SKDEN (P, TEMP, ACT(1,IM), SDEN(1,IM), 
     1                  ISKWRK(IMISK(IM)), RSKWRK(IMRSK(IM)), DEN(1,IM))
C
C  CALCULATE THE BASE GROWTH RATES, GRATE, IN CM / SEC
C
            DO 2300 IPHASE = NFBULK(IM), NLBULK(IM)
               GRATE(IPHASE,IM) = 0.0
               DO 2200 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                  K = KM + KSTOT
                  GRATE(IPHASE,IM) = GRATE(IPHASE,IM) +
     1                        WT(KM,IM)/DEN(KM,IM) * SDOT(K)
2200           CONTINUE
2300        CONTINUE
C
            DO 3400 L = 1 , NATJ
C               (PERTURB THE LTH SOLUTION VARIABLE)
               SAVE = S(L)
               PERTRB = ABS(S(L)) * RELAT + ABSOL
               S(L) = S(L) + PERTRB
C                   (RECALC SDOT WITH A PERTURBED S(L), SDOTD)
               CALL CKYTX (S(NY), ICKWRK, RCKWRK, X)
               DO 2400 K = 1, KKGAS
                  ACT(K,IM) = X(K)
2400           CONTINUE
               DO 2500 K = 1, KKTOT
                  SDOTD(K) = 0.0
                  SDOTT(K) = 0.0
2500           CONTINUE
               IF (NNSURF(IM) .GT. 0) THEN
                  IF (NSPHCH(IM) .GT. 0) THEN
                     DO 2600 IDM =1, NSPHCH(IM)
                        I = IDM + IDTOT
                        SDEN(MAPPH(IDM,IM),IM) 
     1                       = S(NSDENS+I) * SDEN0(MAPPH(IDM,IM),IM)
2600                 CONTINUE
                  ENDIF
                  DO 2800 KM = KFIRST(NFSURF(IM),IM), 
     1                                          KLAST(NLSURF(IM),IM)
                     K = KM + KSTOT
                     ACT(KM,IM) = S(NYS+K)
                     X(K)   = S(NYS+K)
2800              CONTINUE
               ENDIF
               IF (NNBULK(IM).GT.0) THEN
                  DO 3000 KM = KFIRST(NFBULK(IM),IM),
     1                                 KLAST(NLBULK(IM),IM)
                     K = KM + KSTOT
                     ACT(KM,IM) = S(NYS+K)
                     X(K) = S(NYS+K)
3000              CONTINUE
               ENDIF

               CALL PSRACT (X, KK, KKBULK, NFBULK, NLBULK, NPHASE, 
     1                      KFIRST, KLAST, NNBULK, NMAT, NPHMAX, KKSURF,
     2                      KKMAX, ACT)
C               CALCULATE THE SURFACE TEMPERATURE
               IF (.NOT. LTDIFF(IM)) TSURF(IM) = S(NT)
C
               CALL PSSDOT (KKION, KEL, KION, KKGAS, KKTOT, KCHG, NMAT, 
     1                      NPHASE, IISUR, KFIRST, KLAST, NFSURF,
     2                      NLSURF, NNSURF, KKSURF, NFBULK, NLBULK, 
     3                      NNBULK, KKBULK, IISUR0, KKMAX, NPHMAX, KK, 
     4                      BHMXI, WT, EIONSH, NIK, AFRAC, 
     5                      RFFREQ, RFAPAR, RFDPAR, LRFSH, LRFCUR, P, 
     6                      ACT, SCRTCH(1,4), TIONP, S(NTE), TSURF, 
     7                      SDEN, ISKWRK, RSKWRK, IMISK, IMRSK, FLXCOR, 
     8                      FLXION, ESHTH, LELSH, SFAC, SDOTT, ROPS,
     9                      SDOTI, SITDTI, SITDOT, SDOTD, SHLOSS, NBHM,
     *                      LOUTSH, BOLTZ, AVOG, EV2K, BIASPW, AREA, 
     1                      HIM(1,1))
C
C                 (CALCULATE DERIVATIVE AND STORE IT IN DSKDPL)
C
               DO 3200 K = KFIRST(NFBULK(IM),IM), KLAST(NLBULK(IM),IM)
                  DSKDPL(K,L) = (SDOTD(K)-SDOT(K)) /PERTRB
3200           CONTINUE
C                  (RESTORE SOLUTION VARIABLE BACK TO BASE CONDITION)
               S(L) = SAVE
               IF (L.EQ.NT .AND. (.NOT.LTDIFF(IM))) TSURF(IM) = S(NT)
3400        CONTINUE
         ENDIF
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
         IDTOT = IDTOT + NSPHCH(IM)
3600  CONTINUE

C          COMPUTE THE RAW SENSITIVITY COEFFICIENTS WITH RESPECT TO
C             THE RATE CONSTANTS, D(MASS FRACTION)/D(RATE CONSTANT).
C
      WRITE (LSAVE) ISENSI
C
      WRITE (LOUT, 9001) EPSS, EPST, EPSG
C
      DO 5800 I = 1, II
C
         CALL CKRAEX (I, RCKWRK, SAVEP)
         DP = RELAT*SAVEP + ABSOL
         CALL CKRAEX (-I, RCKWRK, SAVEP+DP )
C
         CALL PSRFUN (KK, KKGAS, KKSURF, KKBULK, NSPHCH, NPTS, NATJ,
     1                NPHASE, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2                NLBULK, KKPHAS, KFIRST, KLAST, KOCC, ETCH, LTIME,
     3                LENRGY, LFLRT, LTDIFF, DT, HIN, YIN, TIN, T,
     4                TSURF, P, TAU, FLRT, V, AREA, Q, HTRN, TAMBNT,
     5                LHTRN, ICKWRK, RCKWRK, ISKWRK, RSKWRK, MAPPH, WT,
     6                SCRTCH(1,1), NIK, SN, S, F,
     7                SCRTCH(1,3), SDOTD, ACT, SDEN, SDEN0, SITDOT,
     8                NT, NTE, NYS, NY, NSDEN, NSDENS, KEL, KKION, 
     9                KCHG, KION, QEK, POWR, LENRGE, TEIN, TE, ROP,
     *                ROPS, IEIMP, EIONSH, II, NQLSE, QLSE, TQLSE,
     1                GFAC, SFAC, TIONP, LCONFN, IEXC, REXC, SDOTI, 
     2                SITDTI, IISUR, KKTOT, NMAT, IISUR0, KKMAX, NPHMAX,
     3                IMISK, IMRSK, IMCSK, AFRAC, ITDEP, BHMXI, FLXION,
     4                FLXCOR, SDOTT, SCRTCH(1,4), ESHTH,
     5                LELSH, LRFSH, LRFCUR, RFFREQ, RFAPAR, RFDPAR,
     6                LWALHB, VISC, THCOND, NBHM, LOUTSH, BOLTZ, AVOG,
     7                EV2K, BIASPW, HIM)
C
         CALL CKRAEX (-I, RCKWRK, SAVEP )
         DO 3800 N = 1, NATJ
            DFDALP(N) =  - (F(N)-FN(N)) / DP
3800     CONTINUE
C
C
C*****precision > double
         CALL DGESL (A, NATJ, NATJ, IPVT, DFDALP, 0)
C*****END precision > double
C
C*****precision > single
C         CALL SGESL (A, NATJ, NATJ, IPVT, DFDALP, 0)
C*****END precision > single
C
C      DFDALP(N) NOW CONTAINS THE RAW SENSIVITY ARRAY DS(N)/DAI
C

C           COMPUTE THE BASE GROWTH RATE AND
C           RAW BULK PHASE PRODUCTION RATES WRT SOLUTION VARIABLES

         KSTOT = 0
         DO 4300 IM = 1, NMAT
            IF (NNBULK(IM) .GT. 0) THEN
               DO 4200 IPHASE = NFBULK(IM), NLBULK(IM)
                  GRATED(IPHASE,IM) = 0.0
                  DO 4100 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                     K = KM + KSTOT
                     GRATED(IPHASE,IM) = GRATED(IPHASE,IM) + WT(KM,IM)
     1                               /DEN(KM,IM) * (SDOTD(K)-SDOT(K))/DP
                     DO 4000 L = 1, NATJ
                        GRATED(IPHASE,IM) = GRATED(IPHASE,IM) +WT(KM,IM)
     1                             /DEN(KM,IM) * DSKDPL(K,L) * DFDALP(L)
4000                 CONTINUE
4100              CONTINUE
4200           CONTINUE
            ENDIF
            KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
4300     CONTINUE
C
C          NORMALIZE THE SENSIVITY COEFFICIENTS
C
         DFDALP(NT) = DFDALP(NT) * SAVEP / S(NT)
         SUM = 0.E0
         DO 4400 L = 1, KKGAS
            SUM = SUM + DFDALP(NYS+L) / WT(L,1)
4400     CONTINUE
         CALL CKMMWY (S(NY), ICKWRK, RCKWRK, WTM)
         DO 4600 K = 1, KKGAS
            IF (S(NYS+K) .NE. 0.0E0) THEN
               DFDALP(NYS+K) = DFDALP(NYS+K) * SAVEP / S(NYS+K)
     1                    - SAVEP * WTM * SUM
            ENDIF
4600     CONTINUE
         IF (KKGAS .LT. KKTOT) THEN
            DO 4800 K = KKGAS+1, KKTOT
               IF (S(NYS+K) .NE. 0.0E0) THEN
                  DFDALP(NYS+K) = DFDALP(NYS+K) * SAVEP / S(NYS+K)
               ENDIF
4800        CONTINUE
         ENDIF
         NNBTOT = 0
         DO 5100 IM = 1, NMAT
            IF (NNBULK(IM) .GT. 0) THEN
               DO 5000 IPHASE = NFBULK(IM), NLBULK(IM)
                  GRATED(IPHASE,IM) = GRATED(IPHASE,IM)
     1                        * SAVEP / MAX(SMALL,ABS(GRATE(IPHASE,IM)))
5000           CONTINUE
            ENDIF
            NNBTOT = NNBTOT + NNBULK(IM)
5100     CONTINUE
C
         IF (NNBTOT .GT. 0) THEN
            WRITE (LSAVE)  I, (DFDALP(N),N=1,NATJ),
     1                     ((GRATED(IPHASE,IM), IPHASE=NFBULK(IM), 
     2                                        NLBULK(IM)),IM = 1, NMAT)
         ELSE
            WRITE (LSAVE)  I, (DFDALP(N),N=1,NATJ)
         ENDIF
C
         LWRITE = .FALSE.
         DO 5200 K = 1, KKTOT
            LWRITE =
     1      ( ((ABS(DFDALP(NYS+K)).GE.EPSS) .AND. ISEN(K)).OR.LWRITE)
            LWRITE =
     1      ( ((ABS(DFDALP(NT)).GE.EPST) .AND. LSENT) .OR. LWRITE )
5200     CONTINUE
         DO 5250 IM = 1, NMAT
            IF (LSENG .AND. NNBULK(IM) .GT. 0) THEN
               DO 5260 IPHASE = NFBULK(IM), NLBULK(IM)
                  LWRITE = 
     1            ( (ABS(GRATED(IPHASE,IM)).GE.EPSG) .OR. LWRITE )
 5260          CONTINUE
            ENDIF
 5250    CONTINUE
C
         IF (LWRITE) THEN
C
C           PRINT THE REACTION HOLLERITH
C
            CALL CKSYMR (I, LOUT, ICKWRK, RCKWRK, CCWORK, LT, IHOL,
     1                   IERR)
            LT = MIN(LT,LPMAX)
            WRITE (LOUT, 9010) I, IHOL(:LT)
C
            KP = 0
            DO 5400 K = 1, KKTOT
               IF ((ABS(DFDALP(NYS+K)).GE.EPSS) .AND. ISEN(K) ) THEN
                  KP = KP + 1
                  KEEP(KP) = K
                  IF (K.LE.KKGAS) THEN 
                     KEEPIM(KP) = 1
                     KEEPS(KP) = K
                  ELSE
                     KSTOT = 0
                     DO 5350 IM = 1, NMAT
                        KMBEG = KKGAS + KSTOT 
                        KMEND = KMBEG + KKSURF(IM) + KKBULK(IM) 
                        IF (K .GT. KMBEG .AND. K .LE. KMEND) THEN
                           KEEPIM(KP) = IM
                           KEEPS(KP) = K - KSTOT
                        ENDIF
                        KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
5350                 CONTINUE
                  ENDIF
                  IF (KP .EQ. 3 ) THEN
                     IF (EPSS .LT. 1.E-03) THEN
                        WRITE (LOUT, 9020)
     1        ( KSYM(KEEPS(L),KEEPIM(L)), DFDALP(NYS+KEEP(L)) , L=1,KP)
                     ELSE
                        WRITE (LOUT, 9025)
     1        ( KSYM(KEEPS(L),KEEPIM(L)), DFDALP(NYS+KEEP(L)) , L=1,KP)
                     ENDIF
                     KP = 0
                  ENDIF
               ENDIF
5400        CONTINUE
            IF (KP .GT. 0 ) THEN
               IF (EPSS .LT. 1.E-03) THEN
                  WRITE (LOUT, 9020)
     1        ( KSYM(KEEPS(L),KEEPIM(L)), DFDALP(NYS+KEEP(L)) , L=1,KP)
               ELSE
                  WRITE (LOUT, 9025)
     1        ( KSYM(KEEPS(L),KEEPIM(L)), DFDALP(NYS+KEEP(L)) , L=1,KP)
               ENDIF
            ENDIF
            IF (LSENT .AND. (ABS(DFDALP(NT)).GE.EPST))
     1         WRITE (LOUT, 9030) DFDALP(NT)
            DO 5700 IM = 1, NMAT
               IF (LSENG .AND. NNBULK(IM) .GT. 0) THEN
                  DO 5600 IPHASE = NFBULK(IM), NLBULK(IM)
                      IF (ABS(GRATED(IPHASE,IM)) .GE. EPSG)
     1                   WRITE (LOUT, 9040) PSYM(IPHASE,IM), 
     2                                      GRATED(IPHASE,IM)
5600              CONTINUE
               ENDIF
5700        CONTINUE
         ENDIF
C
5800  CONTINUE
C
C              DO SENSITIVITY ANALYSIS FOR SURFACE REACTIONS
C 
      IF (IISMAX .EQ. 0) RETURN
      KSTOT = 0
      GDION = 0.0
      GDNEU = 0.0
      DO 8200 IM = 1, NMAT
        IF (IISUR(IM) .GT. 0) THEN
         IF (NMAT .GT. 1) THEN
            WRITE (LOUT, 9002) MSYM(IM)
         ELSE
            WRITE (LOUT, 9003)
         ENDIF
         DO 8000 I = 1, IISUR(IM)
C
            CALL SKRAEX (I, ISKWRK(IMISK(IM)), RSKWRK(IMRSK(IM)), SAVEP)
            DP = RELAT*SAVEP + ABSOL
            CALL SKRAEX (-I, ISKWRK(IMISK(IM)), RSKWRK(IMRSK(IM)), 
     1                   SAVEP+DP )
C
            CALL PSRFUN (KK, KKGAS, KKSURF, KKBULK, NSPHCH, NPTS, NATJ,
     1                NPHASE, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2                NLBULK, KKPHAS, KFIRST, KLAST, KOCC, ETCH, LTIME,
     3                LENRGY, LFLRT, LTDIFF, DT, HIN, YIN, TIN, T,
     4                TSURF, P, TAU, FLRT, V, AREA, Q, HTRN, TAMBNT,
     5                LHTRN, ICKWRK, RCKWRK, ISKWRK, RSKWRK, MAPPH, WT,
     6                SCRTCH(1,1), NIK, SN, S, F,
     7                SCRTCH(1,3), SDOTD, ACT, SDEN, SDEN0, SITDOT,
     8                NT, NTE, NYS, NY, NSDEN, NSDENS, KEL, KKION, 
     9                KCHG, KION, QEK, POWR, LENRGE, TEIN, TE, ROP, 
     *                ROPS, IEIMP, EIONSH, II, NQLSE, QLSE, TQLSE, GFAC,
     1                SFAC, TIONP, LCONFN, IEXC, REXC, SDOTI, SITDTI,
     2                IISUR, KKTOT, NMAT, IISUR0, KKMAX, NPHMAX, IMISK,
     3                IMRSK, IMCSK, AFRAC, ITDEP, BHMXI, FLXION,
     4                FLXCOR, SDOTT, SCRTCH(1,4), ESHTH,
     5                LELSH, LRFSH, LRFCUR, RFFREQ, RFAPAR, RFDPAR,
     6                LWALHB, VISC, THCOND, NBHM, LOUTSH, BOLTZ, AVOG,
     7                EV2K, BIASPW, HIM)
C
            CALL SKRAEX (-I, ISKWRK(IMISK(IM)), RSKWRK(IMRSK(IM)),
     1                   SAVEP )
            DO 6000 N = 1, NATJ
               DFDALP(N) =  - (F(N)-FN(N)) / DP
6000        CONTINUE
C
C
C*****precision > double
            CALL DGESL (A, NATJ, NATJ, IPVT, DFDALP, 0)
C*****END precision > double
C
C*****precision > single
C            CALL SGESL (A, NATJ, NATJ, IPVT, DFDALP, 0)
C*****END precision > single
C
C     DFDALP(N) NOW CONTAINS THE RAW SENSIVITY ARRAY DS(N)/DAi
C

C           COMPUTE THE BASE GROWTH RATE AND
C           RAW BULK PHASE PRODUCTION RATES WRT SOLUTION VARIABLES

            IF (NNBULK(IM) .GT. 0) THEN
               DO 6400 IPHASE = NFBULK(IM), NLBULK(IM)
                  GRATED(IPHASE,IM) = 0.0
                  DO 6300 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                     K = KM + KSTOT
                     GRATED(IPHASE,IM) = GRATED(IPHASE,IM) 
     1                                   + WT(KM,IM)/DEN(KM,IM)
     2                                   * (SDOTD(K)-SDOT(K))/DP
                     DO 6200 L = 1, NATJ
                        GRATED(IPHASE,IM) = GRATED(IPHASE,IM) 
     1                                      + WT(KM,IM)/DEN(KM,IM)
     2                                      * DSKDPL(K,L) * DFDALP(L)
6200                 CONTINUE
6300              CONTINUE
6400           CONTINUE
            ENDIF
C
C          NORMALIZE THE SENSIVITY COEFFICIENTS
C
            DFDALP(NT) = DFDALP(NT) * SAVEP / S(NT)
            SUM = 0.E0
            DO 6600 L = 1, KKGAS
               SUM = SUM + DFDALP(NYS+L) / WT(L,1)
6600        CONTINUE
            CALL CKMMWY (S(NY), ICKWRK, RCKWRK, WTM)
            DO 6800 K = 1, KKGAS
              IF (S(NYS+K) .NE. 0.0E0) THEN
                 DFDALP(NYS+K) = DFDALP(NYS+K) * SAVEP / S(NYS+K)
     1                    - SAVEP * WTM * SUM
              ENDIF
6800        CONTINUE
            IF (KKGAS .LT. KK(IM)) THEN
               DO 7000 KM = KKGAS+1, KK(IM)
                  K = KM + KSTOT
                  IF (S(NYS+K) .NE. 0.0E0) THEN
                     DFDALP(NYS+K) = DFDALP(NYS+K) * SAVEP / S(NYS+K)
                  ENDIF
7000           CONTINUE
            ENDIF
            IF (NNBULK(IM) .GT. 0) THEN
               DO 7200 IPHASE = NFBULK(IM), NLBULK(IM)
                  GRATED(IPHASE,IM) = GRATED(IPHASE,IM)
     1                       * SAVEP / MAX(SMALL,ABS(GRATE(IPHASE,IM)))
7200           CONTINUE
            ENDIF
C
            IF (NNBULK(IM) .GT. 0) THEN
               WRITE (LSAVE)  I, (DFDALP(N),N=1,NATJ),
     1               (GRATED(IPHASE,IM), IPHASE=NFBULK(IM), NLBULK(IM))
            ELSE
               WRITE (LSAVE)  I, (DFDALP(N),N=1,NATJ)
            ENDIF
C
            LWRITE = .FALSE.
            DO 7400 KM = 1, KK(IM)
               K = KM + KSTOT
               LWRITE =
     1         ( ((ABS(DFDALP(NYS+K)).GE.EPSS) .AND. ISEN(K)).OR.LWRITE)
               LWRITE =
     1         ( ((ABS(DFDALP(NT)).GE.EPST) .AND. LSENT) .OR. LWRITE )
7400        CONTINUE
            IF (LSENG .AND. NNBULK(IM) .GT. 0) THEN
               DO 7460 IPHASE = NFBULK(IM), NLBULK(IM)
                  LWRITE = 
     1            ( (ABS(GRATED(IPHASE,IM)).GE.EPSG) .OR. LWRITE )
 7460          CONTINUE
            ENDIF
C
            IF (LWRITE) THEN
C
C       PRINT THE REACTION HOLLERITH
C
               CALL SKSYMR (I, LOUT, ISKWRK(IMISK(IM)), 
     1                      RSKWRK(IMRSK(IM)), CSWORK(IMCSK(IM)), 
     2                      LT, IHOL, IERR)
               LT = MIN(LT,LPMAX)
               WRITE (LOUT, 9015) IM, I, IHOL(:LT)
C
               KP = 0
               DO 7600 KM = 1, KK(IM)
                  K = KM + KSTOT
                  IF ((ABS(DFDALP(NYS+K)).GE.EPSS) .AND. ISEN(K) ) THEN
                     KP = KP + 1
                     KEEP(KP) = K
                     KEEPS(KP) = KM
                     IF (KP .EQ. 3 ) THEN
                        IF (EPSS .LT. 1.E-03) THEN
                           WRITE (LOUT, 9020)
     1              ( KSYM(KEEPS(L),IM), DFDALP(NYS+KEEP(L)) , L=1,KP)
                        ELSE
                           WRITE (LOUT, 9025)
     1              ( KSYM(KEEPS(L),IM), DFDALP(NYS+KEEP(L)) , L=1,KP)
                        ENDIF
                        KP = 0
                     ENDIF
                  ENDIF
7600           CONTINUE
               IF (KP .GT. 0 ) THEN
                  IF (EPSS .LT. 1.E-03) THEN
                     WRITE (LOUT, 9020)
     1              ( KSYM(KEEPS(L),IM), DFDALP(NYS+KEEP(L)) , L=1,KP)
                  ELSE
                     WRITE (LOUT, 9025)
     1              ( KSYM(KEEPS(L),IM), DFDALP(NYS+KEEP(L)) , L=1,KP)
                  ENDIF
               ENDIF
               IF (LSENT .AND. (ABS(DFDALP(NT)).GE.EPST))
     1            WRITE (LOUT, 9030) DFDALP(NT)
               IF (LSENG .AND. NNBULK(IM) .GT. 0) THEN
                  DO 7800 IPHASE = NFBULK(IM), NLBULK(IM)
                     IF (ABS(GRATED(IPHASE,IM)) .GE. EPSG)
     1                  WRITE (LOUT, 9040) PSYM(IPHASE,IM), 
     2                                     GRATED(IPHASE,IM)
7800              CONTINUE
               ENDIF
            ENDIF
C 
C  Calculate the relative anisotropy of the deposition rate by finding
C  which reactions are bohm reactions; for topographical simulation
C
            IF (LTERRN) THEN
               CALL SKIBHM (I, ISKWRK, IBHMFL)
               DO 7810 IPHASE = NFBULK(IM), NLBULK(IM)
                  IF (IBHMFL.GT.0) THEN
                     BMCOR = 1.0
                     IF (BHMXI.GT.0) THEN
                        DO 7805 KI = 1, KKION
                           IF (IBHMFL.EQ.KION(KI)) BMCOR = FLXCOR(KI)
 7805                   CONTINUE
                     ENDIF
                     GDION = GDION + MAX(GRATED(IPHASE,IM),ZEROE)
     1                               /BMCOR
                  ELSE
                     GDNEU = GDNEU + MAX(GRATED(IPHASE,IM),ZEROE)
                  ENDIF
                  ANISOT = GDION / MAX((GDION+GDNEU),SMALL)
 7810          CONTINUE
            ENDIF
C
8000     CONTINUE
        ENDIF
        KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
8200  CONTINUE
C
      RETURN
 9001 FORMAT (//,120('=')//20X,
     1 '    LOGARITHMIC SENSITIVITY COEFFICIENTS FOR',
     1 ' GAS PHASE REACTIONS ',//25X,
     1 'Threshold normalized value for printing of SEN ',
     1  'coefficients = ',G11.5/25X,
     1 'Threshold normalized value for printing of ',
     1  'temperature sensitivity coefficient = ',G11.5/25X,
     1 'Threshold normalized value for printing of ',
     1  'growth rate sensitivity coefficient = ',G11.5//)
 9002 FORMAT (///20X,
     1 '    LOGARITHMIC SENSITIVITY COEFFICIENTS FOR ',
     1 'SURFACE PHASE REACTIONS OF MATERIAL:  ',A//)
 9003 FORMAT (///20X,
     1 '    LOGARITHMIC SENSITIVITY COEFFICIENTS FOR ',
     1 'SURFACE PHASE REACTIONS:'//)
 9010 FORMAT (/'  Gas Phase Reaction', I4, '.  ', A)
 9015 FORMAT (/'  Surface Reaction on Material',I4,':', I4, '.  ', A)
 9020 FORMAT ( 3(5X, A15, '= ', 1PE9.2))
 9025 FORMAT ( 3(5X, A15, '= ', F9.3))
 9030 FORMAT ( 5X,'TEMP      = ',1PE9.2)
 9040 FORMAT ( 5X,'GR(',A10,') = ',1PG11.4)
      END
C
C--------------------------------------------------------------------
C
      SUBROUTINE PSRROP (LOUT, LSAVE, LRECOV, KK, II, IISUR, NSPHCH,
     1                   KKGAS, KKSURF, KKBULK, NPHASE, NNSURF, NFSURF,
     2                   NLSURF, NNBULK, NFBULK, NLBULK, NATJ, KKPHAS,
     3                   KFIRST, KLAST, KSYM, PSYM, P, TSURF, V, AREA,
     4                   S, X, ACT, SDEN, SDEN0, ICKWRK, RCKWRK, 
     5                   CCWORK, ISKWRK, RSKWRK, CSWORK, MAPPH, LROP,
     5                   IROP, EPSR, ROP, ROPS, NIK, CIK, CIKS, CIKN,
     7                   CIKNS, NT, NTE, NYS, NY, NSDEN, NSDENS, TIONP,
     8                   LCONFN, KION, KCHG, KKTOT, NMAT, KKMAX, NPHMAX,
     9                   IMISK, IMRSK, IMCSK, IISUR0, MSYM, AFRAC,BHMXI, 
     *                   FLXION, FLXCOR, SFAC, SDOTT, CONC,  
     1                   EIONSH, WT, KEL, KKION, RFFREQ, RFAPAR, RFDPAR,
     2                   LRFSH, LRFCUR, ESHTH, LELSH, SDOT, SITDOT,
     3                   SDOTI, SITDTI, KERR, LOUTSH, BOLTZ, AVOG, EV2K,
     4                   BIASPW, HIM)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KKPHAS(NPHMAX,NMAT), KFIRST(NPHMAX,NMAT), 
     1          KLAST(NPHMAX,NMAT), ICKWRK(*), ISKWRK(*),
     1          MAPPH(NPHMAX,NMAT), S(NATJ), ROP(*), ROPS(IISUR0,NMAT),
     2          X(KKTOT), ACT(KKMAX,NMAT), SDEN(NPHMAX,NMAT),
     3         SDEN0(NPHMAX,NMAT),CIK(NIK),CIKS(IISUR0,NMAT),CIKN(NIK),
     3          CIKNS(IISUR0,NMAT),RCKWRK(*),RSKWRK(*),TEMP(3),
     4          KCHG(KKMAX), KION(*), KK(NMAT), KKSURF(NMAT), 
     5          KKBULK(NMAT), NPHASE(NMAT), NNSURF(NMAT), NFSURF(NMAT),
     6          NLSURF(NMAT), NNBULK(NMAT), NFBULK(NMAT), NLBULK(NMAT),
     6          IISUR(NMAT), NSPHCH(NMAT), IMISK(NMAT), IMRSK(NMAT),
     7          IMCSK(NMAT), AFRAC(NMAT), CONC(NIK), EIONSH(NMAT),
     9          FLXION(*), FLXCOR(*), WT(KKMAX,NMAT), SDOTT(KKTOT),
     *          SDOT(KKTOT), SITDOT(NPHMAX,NMAT), ESHTH(NMAT),
     1          TSURF(NMAT), SDOTI(KKMAX,IISUR0,NMAT), 
     2          SITDTI(NPHMAX,NMAT), BIASPW(NMAT), HIM(KKMAX,NMAT)
C
      LOGICAL LROP, IROP(KKTOT), IERR, LCONFN, LRFSH(NMAT), 
     1        LRFCUR, LELSH(NMAT), KERR
      CHARACTER*16 KSYM(KKMAX,NMAT), PSYM(NPHMAX,NMAT), CCWORK(*), 
     1             CSWORK(*), MSYM(NMAT)
      CHARACTER*18  DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO
      CHARACTER*80  IHOL
      PARAMETER (SMALL = 1.E-25, LPMAX = 50)
C
C               CHARACTER DATA FILE COMMON BLOCK
C
      COMMON /CDATAF/ VERSNN, DHEAD, ICHEM, ISOLUT, ISENSI, IROPRO
C
C                      EXTERNALS
C
      INTEGER PSRPID
      EXTERNAL PSRPID
C
      DATA CNEPS /1.0E-20/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   EPSR   - TRESHOLD VALUE FOR RATE-OF-PRODUCTION COEFFICIENTS.
C   II     - NUMBER OF ELEMENTARY CHEMICAL REACTIONS.
C   IROP   - LOGICAL ARRAY. IF IROP(K)=.TRUE. THEN THE RATE-OF-PRODUC-
C            TION COEFFICIENTS OF SPECIES K WILL BE PRINTED OUT.
C              DIMENSION IROP(*) AT LEAST KK.
C   KK     - NUMBER OF CHEMICAL SPECIES.
C   KSYM   - CHEMKIN SPECIES NAMES.
C              DIMENSION KSYM(*) AT LEAST KK.
C   LOUT   - UNIT FOR PRINTED OUTPUT.
C   LRECOV - UNIT FOR RECOVER FILE
C   LROP   - IF LROP=.TRUE. THEN A RATE-OF-PRODUCTION ANALYSIS
C            IS CARRIED OUT.
C   LSAVE  - UNIT FOR SAVE FILE.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES. NATJ=KK+2.
C   P      - THE PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   S      - DEPENDENT VARIABLE ARRAY. THE TEMPERATURE IS STORED IN
C            T=S(NT), AND THE MASS FRACTIONS ARE IN Y(K)=S(NYS+K)
C              DIMENSION S(*) AT LEAST NATJ.
C
C WORK AND SCRATCH SPACE-
C   RCKWRK  - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      KERR = .FALSE.
      WRITE (LOUT, 6001) EPSR
      WRITE (LSAVE) IROPRO
C
C         CALCULATE RATE-OF-PRODUCTION IN THE GAS
C
      TEMP(1) = S(NT)
      TEMP(2) = S(NTE)
      TEMP(3) = TIONP
      CALL CKQYP (P, TEMP, S(NY), ICKWRK, RCKWRK, ROP)
C
C         CALCULATE RATE-OF-PRODUCTION AT SURFACES
C
C  DETERMINE BULK ACTIVITIES
C
      CALL PSRACT (X, KK, KKBULK, NFBULK, NLBULK, NPHASE, KFIRST, 
     1             KLAST, NNBULK, NMAT, NPHMAX, KKSURF, KKMAX, ACT)
C
C
C  FILL IN REST OF ACT ARRAY AND CALCULATE SDEN
C
      KSTOT = 0
      IDTOT = 0
      DO 50 IM = 1, NMAT
         IF (NNSURF(IM) .GT. 0) THEN
            IF (NSPHCH(IM) .GT. 0) THEN
               DO 5 IDM = 1, NSPHCH(IM)
                  I = IDM + IDTOT
                  SDEN(MAPPH(IDM,IM),IM) 
     1                      = S(NSDENS+I) * SDEN0(MAPPH(IDM,IM),IM)
5              CONTINUE
            ENDIF
            DO 10 KM = KFIRST(NFSURF(IM),IM), KLAST(NLSURF(IM),IM)
               K = KM + KSTOT
               X(K) = S(NYS+K)
               ACT(KM,IM) = X(K)
10          CONTINUE
         ENDIF
         IF (NNBULK(IM) .GT. 0) THEN
            DO 15 KM = KFIRST(NFBULK(IM),IM), KLAST(NLBULK(IM),IM)
               K = KM + KSTOT
C               X(K) = S(NYS+K)
15          CONTINUE
         ENDIF
         CALL CKYTX (S(NY), ICKWRK, RCKWRK, X)
         DO 20 K = 1, KKGAS
            ACT(K,IM) = X(K)
 20      CONTINUE
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
         IDTOT = IDTOT + NSPHCH(IM)
 50   CONTINUE
C
C
C  TEMPORARILY USE X ARRAY FOR ION ENERGIES
C
      CALL PSSDOT (KKION, KEL, KION, KKGAS, KKTOT, KCHG, NMAT, NPHASE,
     1             IISUR, KFIRST, KLAST, NFSURF, NLSURF, NNSURF, KKSURF,
     2             NFBULK, NLBULK, NNBULK, KKBULK, IISUR0, KKMAX, 
     3             NPHMAX, KK, BHMXI, WT, EIONSH, NIK, AFRAC,
     4             RFFREQ, RFAPAR, RFDPAR, LRFSH, LRFCUR, P, ACT,
     5             CONC, TIONP, S(NTE), TSURF, SDEN, ISKWRK, RSKWRK,
     7             IMISK, IMRSK, FLXCOR, FLXION, ESHTH, LELSH, SFAC,
     8             SDOTT, ROPS, SDOTI, SITDTI, SITDOT, SDOT, SHLOSS,
     9             NBHM, LOUTSH, BOLTZ, AVOG, EV2K, BIASPW, AREA, 
     *             HIM(1,1))
C
      DO 1000 K = 1, KKTOT
         IF (K .LE. KKGAS) THEN
            CALL CKCONT (K, ROP, ICKWRK, RCKWRK, CIK)
         ELSE
            DO 90 L = 1, II
               CIK(L) = 0.0
90          CONTINUE
         ENDIF
         KSTOT = 0
         DO 95 IM = 1, NMAT
            IF (K .LE. KKGAS) THEN
               CALL SKCONT (K, ROPS(1,IM), ISKWRK(IMISK(IM)), 
     1                      RSKWRK(IMRSK(IM)), CIKS(1,IM))
            ELSEIF (K .GT. KKGAS+KSTOT .AND. K .LE. KK(IM)+KSTOT) THEN
               KM = K - KSTOT
               CALL SKCONT (KM, ROPS(1,IM), ISKWRK(IMISK(IM)), 
     1                      RSKWRK(IMRSK(IM)), CIKS(1,IM))
            ELSE
               DO 93 L = 1, IISUR(IM)
                  CIKS(L,IM) = 0.0
 93            CONTINUE
            ENDIF
            KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
95       CONTINUE

         WRITE (LSAVE)  K, (CIK(L),L=1,II),
     1                     ((CIKS(L,IM),L=1,MAX(1,IISUR(IM))), 
     1                      IM=1,NMAT)
C
         IF (IROP(K)) THEN
C
            IF (K .LE. KKGAS ) THEN
               WRITE (LOUT, 6010) K, KSYM(K,1)
            ELSE
               KSTOT = 0
               DO 98 IM = 1, NMAT
                  KMBEG = KKGAS + KSTOT 
                  KMEND = KMBEG + KKSURF(IM) + KKBULK(IM) 
                  IF (K.GT.KMBEG .AND. K.LE.KMEND) THEN
                     KSP = K - KSTOT
                     IPHASE = PSRPID(KSP, KLAST(1,IM), NPHASE(IM), 
     1                               LOUT)
                     IF (IPHASE .LE. 0) THEN
                        KERR = .TRUE.
                        RETURN
                     ENDIF
                     IF (IPHASE .LE. NLSURF(IM)) THEN
                        WRITE (LOUT, 6011) K, KSYM(KSP,IM), 
     1                                     PSYM(IPHASE,IM), MSYM(IM)
                     ELSE
                        WRITE (LOUT, 6012) K, KSYM(KSP,IM), 
     1                                     PSYM(IPHASE,IM), MSYM(IM)
                     ENDIF
                  ENDIF
                  KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
98             CONTINUE
            ENDIF
C
C             NORMALIZATION
C
            CNORM1 = 0.0
            CNORM2 = 0.0
            DO 200 I = 1, II
               IF (CIK(I) .LT. 0.0) CNORM1=CNORM1-CIK(I)
               IF (CIK(I) .GT. 0.0) CNORM2=CNORM2+CIK(I)
               CIKN(I)=0.0
200         CONTINUE
            DO 215 IM = 1, NMAT
               DO 210 I = 1, IISUR(IM)
                  IF (CIKS(I,IM) .LT. 0.0)
     1                CNORM1=CNORM1-CIKS(I,IM)*AREA/V
                  IF (CIKS(I,IM) .GT. 0.0)
     1                CNORM2=CNORM2+CIKS(I,IM)*AREA/V
                  CIKNS(I,IM)=0.0
210            CONTINUE
215         CONTINUE

            DO 300 I = 1, II
               IF (CIK(I) .LT. 0.0) THEN
                  CIKN(I) = CIK(I) / MAX(CNEPS, CNORM1)
               ELSE
                  CIKN(I) = CIK(I) / MAX(CNEPS, CNORM2)
               ENDIF
               IF (ABS(CIKN(I)) .GE. EPSR) THEN
                  CALL CKSYMR (I, LOUT, ICKWRK, RCKWRK, CCWORK, LT,
     1                         IHOL, IERR)
                  LT = MIN(LT,LPMAX)
                  WRITE (LOUT, 6020) I, IHOL(:LT), CIKN(I), CIK(I)
               ENDIF
300         CONTINUE
            DO 320 IM = 1, NMAT
               DO 310 I = 1, IISUR(IM)
                  IF (CIKS(I,IM) .LT. 0.0) THEN
                     CIKNS(I,IM) = CIKS(I,IM)
     1                              / MAX(CNEPS, CNORM1) * AREA / V
                  ELSE
                     CIKNS(I,IM) = CIKS(I,IM)
     1                              / MAX(CNEPS, CNORM2) * AREA / V
                  ENDIF
                  IF (ABS(CIKNS(I,IM)) .GE. EPSR) THEN
                     CALL SKSYMR (I, LOUT, ISKWRK(IMISK(IM)), 
     1                            RSKWRK(IMRSK(IM)), CSWORK(IMCSK(IM)),
     1                            LT, IHOL, IERR)
                     LT = MIN(LT,LPMAX)
                     WRITE (LOUT, 6025) I, IHOL(:LT), CIKNS(I,IM),
     1                               CIKS(I,IM)*AREA/V
                  ENDIF
310            CONTINUE
320         CONTINUE
C
            WRITE (LOUT, 6015) CNORM2, CNORM1, (CNORM2-CNORM1)
C
         ENDIF
1000  CONTINUE
 6001 FORMAT (//120('=')//10X,
     1   'NORMALIZED AND ABSOLUTE RATE-OF-PRODUCTION COEFFICIENTS',
     1  /15X,'(surface reactions are normalized by the surface ',
     1   'area to volume ratio)',
     1  /15X,'Threshold normalized value for printing of ROP ',
     1  'coefficients = ', G10.3//)
 6010 FORMAT (/2X,'Gas Phase Species ', I3, '. ', A10, 58X,
     1       'NORMALIZED', 2X,
     1                                    '(MOLES/CC-SEC)' )
 6011 FORMAT (/2X,'Surface Phase Species ', I3, '. ', A10,
     1       ' In Surface Phase ',A10,' of Material ',A10,':', 2X,
     1       'NORMALIZED', 2X,
     1                                    '(MOLES/CC-SEC)' )
 6012 FORMAT (/2X,'Bulk Phase Species ', I3, '. ', A10,
     1       ' In Bulk Phase ',A10,' of Material ',A10,':', 8X,
     1       'NORMALIZED', 2X,
     1                                    '(moles/cc-sec)' )
 6015 FORMAT (21X, 'Total Rate-of-production  (moles/cc-sec) = ',
     1        1PE11.4,
     1       /21X, 'Total Rate-of-consumption (moles/cc-sec) = ',
     1        1PE11.4,
     1       /21X, 'Net   Rate-of-production  (moles/cc-sec) = ',
     1        1PE11.4 )
 6020 FORMAT (10X,    'Gas Phase Reaction', I4,'. ',
     1           A50, 5X, F9.3, 5X, '(', 1PE11.4, ')')
 6025 FORMAT (6X,'Surface Phase Reaction', I4,'. ',
     1           A50, 5X, F9.3, 5X, '(', 1PE11.4, ')')
C
      RETURN
      END
C
C-------------------------------------------------------------------
C
      INTEGER FUNCTION PSRPID (KSPEC, KLAST, NPHASE, LOUT)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER KSPEC, NPHASE, LOUT, IPHASE, KLAST(NPHASE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PSRPID = 0
      DO 10 IPHASE = 1, NPHASE
         IF (KSPEC .LE. KLAST(IPHASE)) THEN
            PSRPID = IPHASE
            RETURN
         ENDIF
10    CONTINUE
      WRITE (LOUT, *)'PSRPID: Failure, kspec = ',KSPEC
C      STOP
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE PSFIXJ (LOUT, KK, KKGAS, KKSURF, KKBULK, NSPHCH, NPTS,
     1                   NATJ, NPHASE, KSYM, PSYM, XLINPK, NNSURF,
     2                   NFSURF, NLSURF, NNBULK, NFBULK, NLBULK, KKPHAS,
     3                   KFIRST, KLAST, KOCC, ABOVE, BELOW, ETCH, LTIME,
     4                   LENRGY, LFLRT, LTDIFF, ABSOL, RELAT, RCOND, DT,
     5                   HIN, YIN, TIN, T, TSURF, P, TAU, FLRT, V, AREA,
     6                   Q, HTRN, TAMBNT, LHTRN, ICKWRK, RCKWRK,
     7                   ISKWRK, RSKWRK, MAPPH, WT, NIK, SCRTCH, SN, S,
     8                   F, FN, A, AA, DR, DC, IPVT, ACT, SDEN, SDEN0,
     9                   SITDOT, NT, NTE, NYS, NY, NSDEN, NSDENS, FINDJ,
     *                   KEL, KKION, KCHG, KION, QEK, POWR, LENRGE,
     1                   TEIN, TE, ROP, IEIMP, EIONSH, II, NQLSE, QLSE,
     2                   TQLSE, GFAC, SFAC, TIONP, LCONFN, IEXC, REXC,
     3                   SDOTI, SITDTI, KKTOT, NMAT, IISUR0, KKMAX,
     4                   NPHMAX, IMISK, IMRSK, IMCSK, IISUR, AFRAC,
     5                   ITDEP, BHMXI, FLXION, FLXCOR, SDOTT, ESHTH,
     6                   LELSH, LRFSH, LRFCUR, RFFREQ, RFAPAR, 
     7                   RFDPAR, LWALHB, VISC, THCOND, SDOT, ROPS,
     8                   NBHM, LOUTSH, BOLTZ, AVOG, EV2K, BIASPW, HIM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KKPHAS(NPHMAX,NMAT), KFIRST(NPHMAX,NMAT), 
     1          KLAST(NPHMAX,NMAT), KOCC(KKMAX,NMAT), ICKWRK(*),  
     2          ISKWRK(*), MAPPH(NPHMAX,NMAT), IPVT(NATJ), S(NATJ),
     2          F(NATJ), FN(NATJ), A(NATJ,NATJ), AA(NATJ,NATJ),
     3          DR(NATJ), DC(NATJ), SN(NATJ), ABOVE(NATJ), BELOW(NATJ),
     5          SCRTCH(NIK,6), WT(KKMAX,NMAT), HIN(KKMAX), YIN(KKMAX),
     5          RCKWRK(*), RSKWRK(*), ACT(KKMAX,NMAT), 
     6          SDEN(NPHMAX,NMAT), SDEN0(NPHMAX,NMAT), 
     7          SITDOT(NPHMAX,NMAT), CONDR(5), KCHG(KKMAX),
     7          KION(*), QEK(KKGAS), ROP(*), IEIMP(*), QLSE(NPTS),
     8          TQLSE(NPTS), IEXC(*), REXC(*), SITDTI(NPHMAX,NMAT),
     9          KKSURF(NMAT), KKBULK(NMAT), KK(NMAT), NPHASE(NMAT),
     *          NNSURF(NMAT), NFSURF(NMAT), NLSURF(NMAT), NNBULK(NMAT), 
     1          NFBULK(NMAT), NLBULK(NMAT), IISUR(NMAT), NSPHCH(NMAT),
     2          IMISK(NMAT), IMRSK(NMAT), IMCSK(NMAT), AFRAC(NMAT),
     3          ITDEP(*), EIONSH(NMAT), FLXION(*), FLXCOR(*), 
     4          SDOTT(KKTOT), ESHTH(NMAT), SDOT(KKTOT),
     5          SDOTI(KKMAX,IISUR0,NMAT), ROPS(IISUR0,NMAT), 
     6          TSURF(NMAT), NBHM(NMAT), BIASPW(NMAT), HIM(KKMAX,NMAT)
C
      CHARACTER KSYM(KKMAX,NMAT)*16, PSYM(NPHMAX,NMAT)*16
C
      LOGICAL LENRGY, LFLRT, LTDIFF(NMAT), LTIME, FINDJ, XLINPK,
     1        ETCH(NPHMAX,NMAT), LHTRN, LFLAG, LENRGE, LCONFN, 
     2        LELSH(NMAT), LRFSH(NMAT), LRFCUR, LWALHB(NMAT)
C
C               EXTERNALS
C
      EXTERNAL PSRJAC
C*****precision > double
      EXTERNAL DGECO
C*****END precision > double
C*****precision > single
C      EXTERNAL SGECO
C*****END precision > single
C*****xlinpk double precision
      EXTERNAL XDGECO
C*****END xlinpk double precision
C*****xlinpk single precision
C      EXTERNAL XSGECO
C*****END xlinpk single precision
C
C           DATA STATEMENTS
      DATA ZERO /0.0/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PRINT OUT ZERO ROWS
C
      LFLAG = .FALSE.
      WRITE (LOUT, 6000)
      DO 20 I = 1, NATJ
         DO 10 J = 1, NATJ
            IF (ABS(A(I,J)) .NE. ZERO)  GO TO 15
10       CONTINUE
         LFLAG = .TRUE.
         IF (I .LE. 2) THEN
            WRITE (LOUT, 6010) I
         ELSEIF (I .LE. KKGAS+2) THEN
            WRITE (LOUT, 6011) I, KSYM(I-2,1)
         ELSEIF (I .LE. KKTOT+2) THEN
            KSTOT = 0
            IMAT = 0
            DO 12 IM = 1, NMAT
               KMBEG = KKGAS + KSTOT 
               KMEND = KMBEG + KKSURF(IM) + KKBULK(IM)
               IF (I-2 .GT. KMBEG .AND. I-2 .LE. KMEND) THEN
                  IMAT = IM
                  KM = I - 2 - KSTOT
               ENDIF
               KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
12          CONTINUE          
            WRITE (LOUT, 6011) I, KSYM(KM,IMAT)
         ELSE
            IDTOT = 0
            IMAT = 0
            DO 14 IM = 1, NMAT
               IMBEG = KKTOT + IDTOT
               IMEND = IMBEG + NSPHCH(IM)
               IF (I-2 .GT. IMBEG .AND. I-2 .LE. IMEND) THEN
                  IMAT = IM
                  IDM = I - 2 - KKTOT - IDTOT
               ENDIF
14          CONTINUE
            WRITE (LOUT, 6012) I, PSYM(MAPPH(IDM,IMAT),IMAT)
         ENDIF
15    CONTINUE
20    CONTINUE
      IF (.NOT. LFLAG) WRITE (LOUT, *)'        NONE'

C PRINT OUT ZERO COLUMNS

      LFLAG = .FALSE.
      WRITE (LOUT, 6100)
      DO 120 J = 1, NATJ
         DO 110 I = 1, NATJ
            IF (ABS(A(I,J)) .NE. ZERO) GO TO 115
110      CONTINUE
         LFLAG = .TRUE.
         IF (J .EQ. 1) THEN
            WRITE (LOUT, 6010) J
         ELSEIF (J .LE. KKGAS+2) THEN
            WRITE (LOUT, 6011) J, KSYM(J-2,1)
         ELSEIF (J .LE. KKTOT+2) THEN
            KSTOT = 0
            IMAT = 0
            DO 112 IM = 1, NMAT
               KMBEG = KKGAS + KSTOT 
               KMEND = KMBEG + KKSURF(IM) + KKBULK(IM)
               IF (J-2 .GT. KMBEG .AND. J-2 .LE. KMEND) THEN
                  IMAT = IM
                  KM = J - 2 - KSTOT
               ENDIF
               KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
112          CONTINUE          
            WRITE (LOUT, 6011) J, KSYM(KM,IMAT)
         ELSE
            IDTOT = 0
            IMAT = 0
            DO 114 IM = 1, NMAT
               IMBEG = KKTOT + IDTOT
               IMEND = IMBEG + NSPHCH(IM)
               IF (J-2 .GT. IMBEG .AND. J-2 .LE. IMEND) THEN
                  IMAT = IM
                  IDM = J - 2 - KKTOT - IDTOT
               ENDIF
114          CONTINUE
            WRITE (LOUT, 6012) J, PSYM(MAPPH(IDM,IMAT),IMAT)
         ENDIF
115   CONTINUE
120   CONTINUE
      IF (.NOT. LFLAG) WRITE (LOUT, *)'        NONE'

C DO DIAGNOSTICS ON FORMULATION OF NUMERICAL JACOBIAN

      ABSOLB = ABSOL
      RELATB = RELAT
      WRITE (LOUT, 6200) ABSOLB, RELATB, (RELATB*10.**(J-2),J=1,5)
      DO 300 I = 1, 5
         ABSOL = ABSOLB * 20.0**(I-2)
         DO 200 J = 1, 5
            RELAT = RELATB * 10.0**(J-2)
            Call PSRJAC (LOUT, KK, KKGAS, KKSURF, KKBULK, NSPHCH,
     1                   NPTS, NATJ, NPHASE, NNSURF, NFSURF, NLSURF,
     2                   NNBULK, NFBULK, NLBULK, KKPHAS, KFIRST,
     3                   KLAST, KOCC, ABOVE, BELOW, ETCH, LTIME,
     4                   LENRGY, LFLRT, LTDIFF, ABSOL, RELAT, DT,
     5                   HIN, YIN, TIN, T, TSURF, P, TAU, FLRT, V,
     6                   AREA, Q, HTRN, TAMBNT, LHTRN, ICKWRK, RCKWRK,
     7                   ISKWRK, RSKWRK, MAPPH, WT, NIK, SCRTCH, SN, 
     8                   S, F, FN, A, ACT, SDEN, SDEN0, SITDOT, NT, NTE, 
     9                   NYS, NY, NSDEN, NSDENS, KEL, KKION, KCHG, KION, 
     *                   QEK, POWR, LENRGE, TEIN, TE, ROP, ROPS, IEIMP,
     1                   EIONSH, II, NQLSE, QLSE, TQLSE, GFAC, SFAC, 
     2                   TIONP, LCONFN, IEXC, REXC, SDOTI, 
     3                   SITDTI, IISUR, KKTOT, NMAT, IISUR0, KKMAX,
     4                   NPHMAX, IMISK, IMRSK, IMCSK, AFRAC, ITDEP,
     5                   BHMXI, FLXION, FLXCOR, SDOTT, ESHTH, LELSH, 
     6                   LRFSH, LRFCUR, RFFREQ, RFAPAR, RFDPAR,
     7                   LWALHB, VISC, THCOND, SDOT, NBHM, LOUTSH,
     8                   BOLTZ, AVOG, EV2K, BIASPW, HIM)
C
            IF (XLINPK) THEN
C*****xlinpk double precision
               CALL XDGECO (A, NATJ, AA, NATJ, NATJ, IPVT,
     1                      RCOND, DR , DC, FN, ANORM)
C*****END xlinpk double precision
C*****xlinpk single precision
C               CALL XSGECO (A, NATJ, AA, NATJ, NATJ, IPVT,
C     1                      RCOND, DR , DC, FN, ANORM)
C*****END xlinpk single precision
            ELSE
C*****precision > double
               CALL DGECO (A, NATJ, NATJ, IPVT, RCOND, FN)
C*****END precision > double
C*****precision > single
C               CALL SGECO (A, NATJ, NATJ, IPVT, RCOND, FN)
C*****END precision > single
            ENDIF
C
            IF (RCOND .NE. ZERO) THEN
               CONDR(J) = 1.0/RCOND
               IF (.NOT. FINDJ) THEN
                  IBEST = I
                  JBEST = J
                  CNDBST = CONDR(J)
               ELSE
                  IF (I+J .GE. 5) THEN
                     IF ((I+J.LT.IBEST+JBEST) .OR.
     1                   (IBEST+JBEST.LT.5)) THEN
                        IBEST = I
                        JBEST = J
                     ELSEIF (I+J .EQ. IBEST+JBEST) THEN
                        IF (CONDR(J) .LE. CNDBST) THEN
                           IBEST = I
                           JBEST = J
                           CNDBST = CONDR(J)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
               FINDJ = .TRUE.
            ELSE
               CONDR(J) = ZERO
            ENDIF
200      CONTINUE
         WRITE (LOUT, 6210) ABSOL, (CONDR(J),J=1,5)
300   CONTINUE
      IF (FINDJ) THEN
         ABSOL = ABSOLB * 20.0**(IBEST-2)
         RELAT = RELATB * 10.0**(JBEST-2)
         RCOND = 1.0 / CNDBST
         WRITE (LOUT, 6220) ABSOL, RELAT, CNDBST
      ENDIF
      WRITE (LOUT, 7000)
      RETURN
6000  FORMAT(//120('=')/
     1        30X,'SINGULAR JACOBIAN DIAGNOSTICS ROUTINE'/120('-')//
     1      10X,'LIST OF ZERO ROWS:')
6010  FORMAT(5X,'SOLUTION UNKNOWN # ',I4,'    TEMPERATURE')
6011  FORMAT(5X,'SOLUTION UNKNOWN # ',I4,'SPECIES NAME = ',A16)
6012  FORMAT(5X,'SOLUTION UNKNOWN # ',I4,'SURFACE PHASE NAME = ',A16)
6100  FORMAT(//10X,'LIST OF ZERO COLUMNS:')
6200  FORMAT(//10X,'DIAGNOSTICS ON NUMERICAL JACOBIAN:'/
     1          15X,'Print out of the condition number'/
     1          20X,'(A value of zero indicates a singular jacobian)'/
     1          15X,'BASE RELATIVE DELTA = ',G11.3/
     1          15X,'BASE ABSOLUTE DELTA = ',G11.3//
     1  20X,5('  RELAT = ',G10.3)/1X,119('-'))
6210  FORMAT(1X,'ABSOL = ',G10.3,'|',5(3X,G12.5,5X) )
6220  FORMAT(/20X,
     1  'PSFIXJ: RETURNING TO TWOPOINT PROBLEM AFTER HAVING ',
     1  'SUCCESSFULLY FOUND NON-SINGULAR JACOBIAN'/
     1  25X,'New ABSOL = ',G13.4/
     1  25X,'New RELAT = ',G13.4/
     1  25X,'New value of condition number = ',G11.3//)
7000  FORMAT(/120('=')//)

      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PSRACT (X, KK, KKBULK, NFBULK, NLBULK, NPHASE, KFIRST, 
     1                   KLAST, NNBULK, NMAT, NPHMAX, KKSURF, KKMAX, 
     2                   ACT)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Subroutine to calculate the activities of the bulk phase
c  species given the mole fractions, temperature and pressure.
C
C  This can be changed by the user.  What is implemented below
C  is the perfect solution approximation.
C
C INPUT
C
C      X  - GAS PHASE MOLE FRACTIONS, SURFACE SITE FRACTIONS
C              AND BULK PHASE MOLE FRACTIONS
C              ORDERED ACCORDING TO THE SURFACE CHEMKIN ORDER
C     KK  - TOTAL NUMBER OF SPECIES
C  KKBULK - TOTAL NUMBER OF BULK PHASE SPECIES IN ALL BULK PHASES
C  NFBULK - POINTER TO THE FIRST BULK PHASE
C  NLBULK - POINTER TO THE LAST BULK PHASE
C  NPHASE - TOTAL NUMBER OF PHASES
C  KFIRST - VECTOR OF POINTERS TO THE FIRST SPECIES IN EACH PHASE
C  KLAST  - VECTOR OF POINTERS TO THE LAST SPECIES IN EACH PHASE
C
C OUTPUT
C
C   ACT   - VECTOR OF ACTIVITIES - DO NOT TOUCH THE NON-BULK
C                PHASE ENTRIES IN THIS VECTOR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KFIRST(NPHMAX,NMAT), KLAST(NPHMAX,NMAT), 
     1          X(*), ACT(KKMAX,NMAT), KKSURF(NMAT), NNBULK(NMAT),
     2          NFBULK(NMAT), NLBULK(NMAT), KKBULK(NMAT), KK(NMAT), 
     3          NPHASE(NMAT)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      KSTOT = 0
      DO 30 IM = 1, NMAT
         IF (NNBULK(IM) .NE. 0) THEN
            DO 20 IPHASE = NFBULK(IM), NLBULK(IM)
               DO 10 KM = KFIRST(IPHASE,IM), KLAST(IPHASE,IM)
                  K = KM + KSTOT
                  ACT(KM,IM) = X(K)
10            CONTINUE
20          CONTINUE
         ENDIF
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
30    CONTINUE
C
      RETURN
      END
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
      CHARACTER KNAM(*)*16, IST*16
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
C----------------------------------------------------------------------
C
      SUBROUTINE PLTEMP (NPTS, X, XX, TT, T)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION XX(*), TT(*)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     THIS SUBROUTINE USES BISECTION TO LINEARLY INTERPOLATE
C     AN ARRAY OF XX,TT PAIRS.  GIVEN AN XX,TT PAIR THIS ROUTINE
C     RETURNS THE INTERPOLATED VALUE OF THE T AT THE POINT X.
C
C INPUT-
C   NPTS   - NUMBER OF XX,TT PAIRS.
C   X      - LOCATION AT WHICH INTERPOLATED T IS DESIRED.
C   XX     - ARRAY OF X POINTS AT WHICH TT ARE GIVEN.
C   TT     - ARRAY OF T VALUES AT THE XX LOCATIONS.
C
C OUTPUT-
C   T     - INTERPOLATED T AT POINT X
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C             check for x outside (1,npts)
C
      IF (X .LE. XX(2)) THEN
         N = 2
         S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
      ELSEIF (X .GE. XX(NPTS-1)) THEN
         N = NPTS-1
         S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
      ELSE
         NLO = 1
         NHI = NPTS
         S   = 0.0
C
C        bisect interval
C
50       CONTINUE
         N = (NLO+NHI)/2
         IF (X .LT. XX(N)) THEN
            IF (X .LT. XX(N-1)) THEN
               NHI = N
               GO TO 50
            ELSEIF (X .EQ. XX(N-1)) THEN
               N = N-1
            ELSE
               S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
            ENDIF
         ELSEIF (X .GT. XX(N)) THEN
            IF (X .GT. XX(N+1)) THEN
               NLO = N
               GO TO 50
            ELSEIF (X .EQ. XX(N+1)) THEN
               N = N + 1
            ELSE
               S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
            ENDIF
         ENDIF
      ENDIF
C
  100 CONTINUE
      T      = TT(N) + S * (X - XX(N))
      RETURN
C
      END
C
      SUBROUTINE PSRBHM (KKION, KEL, KCHG, KION, KKGAS, C, TI, TE,
     1                   WT, BHMXI, FLXION)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Subroutine to calculate  a correction factor for the Bohm-
c  velocity surface reaction rate-of-production.  When negative ions
C  are present, the positive ion velocity cannot reach the Bohm 
C  velocity and is limited by the thermal speed of negative ions.
C
C INPUT
C  KKION  - TOTAL NUMBER OF IONS
C    KEL  - INDEX OF THE ELECTRON IN THE MOLE FRACTION ARRAY
C   KCHG  - INTEGER ARRAY OF SPECIES ELECTRONIC CHARGE
C           DIMENSION AT LEAST KCHG(KKGAS)
C   KION  - INTEGER ARRAY OF SPECIES INDICES FOR THE KKION IONS
C           DIMENSION AT LEAST KION(KKION)
C  KKGAS  - THE TOTAL NUMBER OF GAS-PHASE SPECIES
C      C  - GAS PHASE MOLAR CONCENTRATIONS
C           DIMENSION AT LEAST C(KKGAS)
C     TI  - ION TEMPERATURE (K)
C     TE  - ELECTRON TEMPERATURE (K)
C     WT  - GAS PHASE MOLECULAR WEIGHTS (G/MOLE)
C  BHMXI  - THE TRANSPORT-LIMITATION MULTIPLIER FOR THE BOHM FLUX
C
C OUTPUT
C
C FLXION  - ION MOLAR FLUX BASED ON BOHM AND NEGATIVE-ION CORRECTION
C           DIMENSION AT LEAST FLXION(KKION)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KCHG(*), KION(*), C(*), FLXION(*), WT(*)

      IF (KKION.EQ.0 .AND. KEL .EQ.0) RETURN
C
      ZEROE = 0.0
      SUM = 0.0
      DO 50 K = 1, KKGAS
         SUM = SUM + MAX(ZEROE,C(K))
 50   CONTINUE
      BHMCOR = 1.0
      SUMNEG = 0.0
      DO 100 KI = 1, KKION
         K = KION(KI)
         IF (KCHG(K).GT.0) THEN
C 
C Ubohm = sqrt (k_boltz * Te / m_i), where m_i = wt(k_ion)/AVOG
C
            FLXION(KI) = C(K) * 9117.76 * SQRT(TE/WT(K)) * BHMXI
         ELSE
            IF (KCHG(K).LT.0) SUMNEG = SUMNEG + MAX(ZEROE,C(K))
            FLXION(KI) = 0.0
         ENDIF
 100  CONTINUE
      XNEG = SUMNEG / SUM
      XEL  = MAX(ZEROE,C(KEL)) / SUM
      DENOM = XNEG*TE+XEL*TI
      IF (DENOM .GT. 0.0 .AND. (XNEG+XEL) .GT. 0.0) THEN
         BHMCOR = SQRT((XNEG+XEL)*TI/ DENOM)
      ELSE
         BHMCOR = 1.0
      ENDIF
C
      DO 200 KI = 1, KKION
         FLXION(KI) = FLXION(KI) * BHMCOR
 200  CONTINUE
C
      RETURN
      END
      SUBROUTINE PSBHM1 (KKION, KEL, KCHG, KION, KKGAS, C, TI, TE,
     1                   FNEG)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Subroutine to calculate  a correction factor for the Bohm-
c  velocity surface reaction rate-of-production.  When negative ions
C  are present, the positive ion velocity cannot reach the Bohm 
C  velocity and is limited by the thermal speed of negative ions.
C
C INPUT
C  KKION  - TOTAL NUMBER OF IONS
C    KEL  - INDEX OF THE ELECTRON IN THE MOLE FRACTION ARRAY
C   KCHG  - INTEGER ARRAY OF SPECIES ELECTRONIC CHARGE
C           DIMENSION AT LEAST KCHG(KKGAS)
C   KION  - INTEGER ARRAY OF SPECIES INDICES FOR THE KKION IONS
C           DIMENSION AT LEAST KION(KKION)
C  KKGAS  - THE TOTAL NUMBER OF GAS-PHASE SPECIES
C      C  - GAS PHASE MOLAR CONCENTRATIONS
C           DIMENSION AT LEAST C(KKGAS)
C     TI  - ION TEMPERATURE (K)
C     TE  - ELECTRON TEMPERATURE (K)
C
C OUTPUT
C
C   FNEG  - THE CORRECTION FACTOR FOR BOHM FLUX WHEN NEGATIVE IONS
C           ARE PRESENT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KCHG(*), KION(*), C(*)

      IF (KKION.EQ.0 .AND. KEL .EQ.0) RETURN
C
      ZEROE = 0.0
      SUM = 0.0
      DO 50 K = 1, KKGAS
         SUM = SUM + MAX(ZEROE,C(K))
 50   CONTINUE
      FNEG = 1.0
      SUMNEG = 0.0
      DO 100 KI = 1, KKION
         K = KION(KI)
         IF (KCHG(K).LT.0) SUMNEG = SUMNEG + MAX(C(K),ZEROE)
 100  CONTINUE
      XNEG = SUMNEG / SUM
      XEL  = MAX(ZEROE,C(KEL)) / SUM
      DENOM = XNEG*TE + XEL*TI
      IF (DENOM .GT. 0.0 .AND. (XNEG+XEL).GT.0) THEN
         FNEG = SQRT((XNEG+XEL)*TI/ DENOM )
      ELSE
         FNEG = 1.0
      ENDIF
C
      RETURN
      END
C
C--------------------------------------------------------------------
C
      SUBROUTINE PSSDOT (KKION, KEL, KION, KKGAS, KKTOT, KCHG, NMAT, 
     1                   NPHASE, IISUR, KFIRST, KLAST, NFSURF, NLSURF, 
     2                   NNSURF, KKSURF, NFBULK, NLBULK, NNBULK, KKBULK,
     3                   IISUR0, KKMAX, NPHMAX, KK, BHMXI, WT, EIONSH,
     4                   NIK, AFRAC, RFFREQ, RFAPAR, RFDPAR, LRFSH,
     5                   LRFCUR, P, ACT, CONC, TIONP, TELEC, TSURF, 
     6                   SDEN, ISKWRK, RSKWRK, IMISK, IMRSK, FLXCOR,
     7                   FLXION, ESHTH, LELSH, SFAC, SDOTT, ROPS,
     8                   SDOTI, SITDTI, SITDOT, SDOT, SHLOSS, NBHM,
     9                   LOUTSH, BOLTZ, AVOG, EV2K, BIASPW, AREA, ENRGI)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KION(*), KCHG(KKMAX), NPHASE(NMAT), IISUR(NMAT), 
     1          KFIRST(NPHMAX,NMAT), KLAST(NPHMAX,NMAT), NFSURF(NMAT), 
     1          NLSURF(NMAT), NNSURF(NMAT), NNBULK(NMAT), 
     2          NFBULK(NMAT), NLBULK(NMAT), KK(NMAT), KKSURF(NMAT),
     3          KKBULK(NMAT), EIONSH(NMAT), AFRAC(NMAT),
     4          WT(KKMAX,NMAT), ACT(KKMAX,NMAT), SDEN(NPHMAX,NMAT), 
     5          CONC(NIK), ISKWRK(*), RSKWRK(*), FLXCOR(*), FLXION(*), 
     6          SDOTT(KKTOT), ROPS(IISUR0,NMAT),
     6          SDOTI(KKMAX,IISUR0,NMAT), SITDTI(NPHMAX,NMAT),
     7          SITDOT(NPHMAX,NMAT), SDOT(KKTOT), NBHM(NMAT),
     8          TEMP(3), IMISK(NMAT), IMRSK(NMAT), ESHTH(NMAT), 
     9          TSURF(NMAT), BIASPW(NMAT), ENRGI(KKMAX)

      LOGICAL LELSH(NMAT), LRFSH(NMAT), LRFCUR
      PARAMETER(SMALL=1.E-25, ZEROE=0.0D0)
      PARAMETER(IDIMRF = 101)
      DIMENSION TIMRF(IDIMRF), VWALL(IDIMRF), VWBAR(IDIMRF), 
     1          EWALL(IDIMRF), AJTOT(IDIMRF)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   KKION  - NUMBER OF IONS IN CHEMKIN MECHANISM.
C   KEL    - INDEX OF THE ELECTRON IN THE SPECIES ARRAY
C   KION   - ARRAY OF ION SPECIES INDICES (NOT INCLUDING THE ELECTRON)
C              DIMENSION KION(*) AT LEAST KKION.
C   KKGAS  - NUMBER OF GAS-PHASE SPECIES IN THE CHEMKIN MECHANISM.
C   KKTOT  - TOTAL NUMBER OF SPECIES IN THE GAS AND SURFACE MECHANISMS.
C   KCHG   - ARRAY OF SPECIES ELECTRONIC CHARGES (INTEGER)
C              DIMENSION KCHG(*) AT LEAST KKGAS.
C   NMAT   - NUMBER OF MATERIALS IN THE SURFACE CHEMKIN MECHANISM.
C   NPHASE - NUMBER OF PHASES FOR EACH MATERIAL IN THE SURFACE MECH.
C              DIMENSION NPHASE(*) AT LEAST NMAT.
C   IISUR  - NUMBER OF SURFACE REACTIONS FOR EACH MATERIAL.
C              DIMENSION IISUR(*) AT LEAST NMAT.
C   KFIRST - INDEX OF THE FIRST SPECIES FOR EACH PHASE IN EACH MATERIAL.
C              DIMENSION KFIRST(NPHMAX,*) AT LEAST NMAT.
C   KLAST  - INDEX OF THE LAST SPECIES FOR EACH PHASE IN EACH MATERIAL.
C              DIMENSION KLAST(NPHMAX,*) AT LEAST NMAT.
C   NFSURF - INDEX OF THE FIRST SURFACE PHASE FOR EACH MATERIAL.
C              DIMENSION NFSURF(*) AT LEAST NMAT.
C   NLSURF - INDEX OF THE LAST SURFACE PHASE FOR EACH MATERIAL.
C              DIMENSION NLSURF(*) AT LEAST NMAT.
C   KKSURF - NUMBER OF SURFACE SPECIES IN EACH MATERIAL.
C              DIMENSION KKSURF(*) AT LEAST NMAT.
C   KKBULK - NUMBER OF BULK SPECIES IN EACH MATERIAL.
C              DIMENSION KKBULK(*) AT LEAST NMAT.
C   NFBULK - INDEX OF THE FIRST BULK PHASE FOR EACH MATERIAL.
C              DIMENSION NFBULK(*) AT LEAST NMAT.
C   NLBULK - INDEX OF THE LAST BULK PHASE FOR EACH MATERIAL.
C              DIMENSION NLBULK(*) AT LEAST NMAT.
C   NNSURF - NUMBER OF SURFACE PHASES IN EACH MATERIAL.
C              DIMENSION NNSURF(*) AT LEAST NMAT.
C   NNBULK - NUMBER OF BULK PHASES IN EACH MATERIAL.
C              DIMENSION NNBULK(*) AT LEAST NMAT.
C   KK     - NUMBER OF SPECIES ASSOCIATED WITH EACH MATERIAL.
C              DIMENSION KK(*) AT LEAST NMAT.
C   KKMAX  - MAXIMUM NUMBER OF SPECIES IN ANY MATERIAL (INCLUDING GAS).
C   NPHMAX - MAXIMUM NUMBER OF PHASES IN ANY MATERIAL (INCLUDING GAS).
C   BHMXI  - THE BOHM VELOCITY CORRECTION FACTOR.
C   WT     - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
C              CGS UNITS - GM/MOLE
C              DIMENSION WT(KK,*) AT LEAST NMAT.
C   EIONSH - THE ION ENERGY AFTER TRAVERSING THE SHEATH 
C              CGS UNITS - K; DIMENSION AT LEAST NMAT.
C   BIASPW - THE BIAS POWER APPLIED TO A MATERIAL BOUNDARY
C              CGS UNITS - erg/s
C              DIMENSION BIASPW(*) AT LEAST NMAT.
C   AFRAC  - THE AREA FRACTION OF EACH MATERIAL IN THE REACTOR.
C              DIMENSION AFRAC(*) AT LEAST NMAT.
C   AREA   - THE TOTAL SURFACE AREA IN THE REACTOR
C              CGS UNITS - CM**2
C   RFFREQ - THE RF FREQUENCY OF THE RF WAFER BIAS, WHEN INCLUDED.
C              CGS UNITS - /SEC
C   RFAPAR - THE AMPLITUDE FACTOR FOR THE RF WAFER BIAS.
C              UNITS - VOLTS OR AMPS DEPENDING ON VOLTAGE/CURR CONTROL.
C   RFDPAR - THE DC PARAMETER FOR THE RF WAFER BIAS.
C              UNITS - VOLTS OR AMPS DEPENDING ON VOLTS/CURR CONTROL.
C   LRFSH  - LOGICAL FLAG FOR SOLVING THE RF SHEATH MODEL.
C              DIMENSION LRFSH(*) AT LEAST NMAT.
C   LRFCUR - LOGICAL FLAG FOR CURRENT CONTROL IN RF SHEATH SIMULATION.
C   P      - THE PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   TIONP  - ION TEMPERATURE BEFORE ENTERING SHEATH
C              CGS UNITS - K
C   TSURF  - SURFACE TEMPERATURE
C              CGS UNITS - K
C   TELEC  - ELECTRON TEMPERATURE
C              CGS UNITS - K
C   SDEN   - PHASE SITE DENSITIES FOR EACH MATERIAL.
C              CGS UNTIS - MOLES/CM**2
C              DIMENSION SDEN(NPHMAX,*) AT LEAST NMAT.
C   IMISK  - INDEX OF THE LOCATION IN THE ISKWRK ARRAY FOR EACH MATERIAL
C              DIMENSION IMISK(*) AT LEAST NMAT.
C   IMRSK  - INDEX OF THE LOCATION IN THE RSKWRK ARRAY FOR EACH MATERIAL
C              DIMENSION IMRSK(*) AT LEAST NMAT.
C   SFAC   - USER MULTIPLYING FACTOR FOR ALL SURFACE CHEMISTRY RATES
C WORK AND SCRATCH SPACE-
C   RSKWRK - REAL SURFACE CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE SURFACE CHEMKIN DOCUMENTATION.
C   ISKWRK - INTEGER SURFACE CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE SURFACE CHEMKIN DOCUMENTATION.
C   ENRGI    - ARRAY OF SPECIES (ION) ENERGIES HITTING A SURFACE.
C              CGS UNITS - ERGS/MOLE
C              DIMENSION ENRGI(*) AT LEAST KK.
C   ACT    - ARRAY OF ACTIVITIES FOR EACH MATERIAL
C              DIMENSION ACT(KKMAX,*) AT LEAST NMAT.
C   CONC   - ARRAY OF SPECIES CONCENTRATIONS.
C              DIMENSION CONC(*) AT LEAST KKMAX.
C   FLXCOR - ARRAY OF FLUX CORRECTION FACTORS FOR ION FLUXES
C              DIMENSION FLXCOR(*) AT LEAST KKION.
C   FLXION - ARRAY OF ION FLUXES
C              DIMENSION FLXION(*) AT LEAST KKION.
C   ESHTH  - ION ENERGY LOSS IN THE SHEATH - FACTOR * ELECTRON ENERGY
C              DIMENSION ESHTH(*) AT LEAST NMAT.
C   LELSH  - LOGICAL FLAG INDICATING WHETHER DEFAULT FORMULATION OR
C            SPECIFIED ION ENERGY LOSS IN SHEATH IS USED.
C   SDOTT  - ARRAY OF TEMPORARY SURFACE PRODUCTION RATES
C              DIMENSION SDOTT(*) AT LEAST KKTOT.
C   ROPS    - ARRAY OF REACTION RATES OF PROGRESS.
C              DIMENSION ROPS(IISUR0,*) AT LEAST NMAT.
C   SDOTI  - ARRAY OF PRODUCTION RATES FOR EACH REACTION.
C              DIMENSION SDOTI(KKMAX,IISUR0,*) AT LEAST NMAT.
C   SITDTI - ARRAY OF PHASE PRODUCTION RATES FOR EACH REACTION.
C              DIMENSION SITDTI(NPHMAX,*) AT LEAST NMAT.
C   NBHM   - ARRAY OF NUMBER OF BOHM REACTIONS ON EACH MATERIAL
C              DIMENSION NBHM(*) AT LEAST NMAT.
C OUTPUT-
C   SDOT   - ARRAY OF CHEMICAL PRODUCTION RATES FROM SURFACE REACTIONS
C              CGS UNITS - MOLES/(CM**2-SEC)
C              DIMENSION SDOT(*) AT LEAST KKTOT.
C   SITDOT - ARRAY OF PRODUCTION RATES OF PHASE SITE DENSITIES
C              CGS UNITS - MOLES/(CM**2-SEC)
C              DIMENSION SITDOT(NPHMAX,*) AT LEAST NMAT.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C initialize ECHG to be molar electron charge:
C
      ECHG=AVOG*BOLTZ*EV2K
C
      DO 500 K = 1, KKTOT
         SDOT(K) = 0.0
500   CONTINUE
C
      DO 520 IM = 1, NMAT
         DO 510 K = 1, KKGAS
            ACT(K,IM) = MAX(ZEROE,ACT(K,IM))
 510     CONTINUE
 520  CONTINUE
C
      TEMP(2) = TELEC
      TEMP(3) = TIONP
      KSTOT = 0
      SHLOSS = 0.0
      DO 1100 IM = 1, NMAT    
         DO 600 N = 1, NPHASE(IM)
            SITDOT(N,IM) = 0.0
600      CONTINUE
         SHLIM = 0.0
C
C  Calculate the gas concentrations for the surface temperature
C
         TEMP(1) = TSURF(IM)
         CALL SKATCZ (P, TEMP, ACT(1,IM), SDEN(1,IM), ISKWRK(IMISK(IM)), 
     1                RSKWRK(IMRSK(IM)),CONC)
C
C  DETERMINE ION FLUX USING BOHM CONDITION AND NEG-ION CORRECTION
C
         FNEG = 1.0
         CURTOT = 0.0
         IF (BHMXI .GT. 0.0 ) THEN
            CALL PSRBHM (KKION, KEL, KCHG, KION, KKGAS, CONC, TIONP, 
     1                   TELEC, WT(1,1), BHMXI, FLXION)
            DO 603 KI = 1, KKION
               CURTOT = CURTOT + FLXION(KI)*ECHG
 603        CONTINUE
         ELSEIF (KEL.NE.0 .AND. NBHM(IM).NE.0) THEN
            CALL PSBHM1 (KKION, KEL, KCHG, KION, KKGAS, CONC, TIONP, 
     1                   TELEC, FNEG)
C
C   for case bias power is specified but ion current not calculated yet,
C   approximate total ion current.
C
            IF (BIASPW(IM) .NE. 0.0 .AND. EIONSH(IM) .EQ. 0.0) THEN
               DO 604 KI = 1, KKION
                  K = KION(KI)
                  VBOHM = FNEG * SQRT(BOLTZ*AVOG*TELEC/WT(K,1))
                  CURTOT = CURTOT + CONC(K)*VBOHM*ECHG
 604           CONTINUE
            ENDIF
         ENDIF
C
C
C  SOLVE FOR THE AVERAGE ION ENERGY USING AN RF SHEATH MODEL IF DESIRED
C
         IF (LRFSH(IM)) THEN
            IF (LRFCUR) THEN
               ICODE = 1
            ELSE
               ICODE = 0
            ENDIF
            TEV = TELEC/EV2K
            IPRNT = 0
C
C   convert to #/m3 for concentrations
C
            DO 605 K = 1, KKGAS
               CONC(K) = CONC(K)*AVOG*1.E6
 605        CONTINUE
            CALL SHEATH (RFFREQ, WT(1,1), TEV, CONC, KCHG, KKION, KION,
     1                   KKGAS, ICODE, RFAPAR, RFDPAR, IPRNT, LOUTSH,
     2                   TIMRF, VWALL, VWBAR, EWALL, AJTOT, VAVEB, IERR)
            EIONSH(IM) = ABS(VAVEB*EV2K)
         ELSEIF (BIASPW(IM).NE. 0.0 ) THEN
            IF (CURTOT .GT. 0.0) THEN
               EIONSH(IM) = BIASPW(IM) * EV2K / (CURTOT*AREA*AFRAC(IM))
            ENDIF
         ENDIF
C
C  DEFINE ION ENERGIES
C  SDOTT IS THE TEMPORARY SDOT ARRAY TO DETERMINE ION FLUX LIMITS
C
         DO 610 K = 1, KKGAS
            ENRGI(K) = 0.0
 610     CONTINUE
         DO 615 KI = 1, KKION
            K = KION(KI)
            ENRGI(K) = EIONSH(IM)
 615     CONTINUE
         DO 620 K = 1, KKTOT
            SDOTT(K) = 0.0
 620     CONTINUE
C     
C  INITIALIZE AUXILIARY PARAMETER INFORMATION FOR RATE CALCULATIONS
C
         CALL SKRPAR(ISKWRK(IMISK(IM)), RSKWRK(IMRSK(IM)), ENRGI)
C
         CALL SKROP (P, TEMP, ACT(1,IM), SDEN(1,IM), 
     1               ISKWRK(IMISK(IM)), RSKWRK(IMRSK(IM)), ROPS(1,IM))
C
         DO 705 I = 1, IISUR(IM)
            CALL SKRATI (I, ROPS(1,IM), ISKWRK(IMISK(IM)), 
     1                   RSKWRK(IMRSK(IM)), 
     1                   SDOTI(1,I,IM), SITDTI(1,IM))
            IF (BHMXI .GT. 0.0 .AND. KKION .GT. 0) THEN
               DO 700 KI = 1, KKION
                  K = KION(KI)
                  SDOTT(K) = SDOTT(K) + SDOTI(K,I,IM)
 700           CONTINUE
            ELSE
               CALL SKIBHM (I, ISKWRK(IMISK(IM)), IBHMFL)
               IF (IBHMFL .GT. 0) THEN
                  DO 703 K = 1, KK(IM)
                     SDOTI(K,I,IM) = SDOTI(K,I,IM)*FNEG
                     ROPS(I,IM) = ROPS(I,IM) * FNEG
 703              CONTINUE
               ENDIF
            ENDIF
 705     CONTINUE
         IF (BHMXI .GT. 0.0 .AND. KKION .GT. 0) THEN
            DO 710 KI = 1, KKION
               K = KION(KI)
               IF (KCHG(K) .GT. 0) THEN
                  FLXCOR(KI) = FLXION(KI)/MAX(ABS(SDOTT(K)),SMALL)
               ELSE
                  FLXCOR(KI) = 1.0
               ENDIF
 710        CONTINUE
         ENDIF
         DO 1000 I = 1, IISUR(IM)
            COR = 1.0
            IF (BHMXI .GT. 0.0 .AND. KKION .GT. 0) THEN
               DO 715 KI = 1, KKION
                  K = KION(KI)
                  IF (ABS(SDOTI(K,I,IM)) .GT. 0.0) THEN
                     COR = COR*FLXCOR(KI)
                  ENDIF
 715           CONTINUE
            ENDIF
            ROPS(I,IM) = ROPS(I,IM) * COR * AFRAC(IM) * SFAC
            DO 730 K = 1, KKGAS
               SDOTI(K,I,IM) = SDOTI(K,I,IM)*COR*AFRAC(IM) * SFAC
               SDOT(K) = SDOT(K) + SDOTI(K,I,IM)
               DO 720 KI = 1, KKION
                  IF (KION(KI) .EQ. K .AND. KCHG(K).GT.0) THEN
                     SHLIM = SHLIM - SDOTI(K,I,IM)
                  ENDIF
 720           CONTINUE
 730        CONTINUE
            IF (NNSURF(IM) .GT. 0) THEN
               DO 750 KM = KFIRST(NFSURF(IM),IM), KLAST(NLSURF(IM),IM)
                  K = KM + KSTOT
                  SDOTI(KM,I,IM) = SDOTI(KM,I,IM)*COR*SFAC
                  SDOT(K) = SDOT(K) + SDOTI(KM,I,IM)
750            CONTINUE
            ENDIF
            IF (NNBULK(IM) .GT. 0) THEN
               DO 770 KM = KFIRST(NFBULK(IM),IM), KLAST(NLBULK(IM),IM)
                  K = KM + KSTOT
                  SDOTI(KM,I,IM) = SDOTI(KM,I,IM)*COR*SFAC
                  SDOT(K) = SDOT(K) + SDOTI(KM,I,IM)
770            CONTINUE
            ENDIF
            DO 800 N = 1, NPHASE(IM)
               SITDTI(N,IM) = SITDTI(N,IM)*COR*SFAC
               SITDOT(N,IM) = SITDOT(N,IM) + SITDTI(N,IM)
800         CONTINUE
C
1000     CONTINUE
C
C Reset ion energy based on surface production rates of ions.
C
         IF (BIASPW(IM).NE. 0.0 ) THEN
            CURTOT = 0.0
            DO 1003 KI = 1, KKION
               K = KION(KI)
               CURTOT = CURTOT - SDOT(K)*ECHG
 1003       CONTINUE
            IF (CURTOT .GT. 0.0) THEN
               EIONSH(IM) = BIASPW(IM) * EV2K / (CURTOT*AREA*AFRAC(IM))
            ENDIF
         ENDIF
         
         KSTOT = KSTOT + KKSURF(IM) + KKBULK(IM)
         IF (LELSH(IM)) THEN
            SHLOSS = SHLOSS + SHLIM*ESHTH(IM)*TELEC
         ELSE
            SHLOSS = SHLOSS + SHLIM*EIONSH(IM)
         ENDIF
 1100 CONTINUE
      RETURN
      END

      subroutine sheath(anurf,amass,etev,ank,kchg,kkion,kion,kkgas,
     1                  icode,xac,xdc,iprnt,lun,time,vw,vwb,ew,aj,vaveb,
     2                  ierr)
       implicit real*8(a-h,o-z)
c This code uses time integration of the sheath model based on
c the approximate first integral of the Poisson equation. Reference:
c Merle E. Riley, "Unified Model of the rf Plasma Sheath," Sandia
c National Laboratories Technical Report SAND95-0775 UC-401, May 1995.
c
c Input:  anurf  - the frequency of the rf field (Hz).
c         amass()- atomic mass of the species (amu).
c                  dimension at least kkgas (number of gas species)
c         etev   - electron temperature (eV).
c         anion0 - electron density at bulk boundary (#/m3).
c         ank()  - gas density for each species (#/m3)
c                  dimension at least kkgas (number of gas species)
c         kchg() - integer charge of the species
c                  dimension at least kkgas (number of gas species)
c         kkgas  - total number of species
c         kkion  - number of positive ions
c         kion() - index array of ions in species array
c                  dimension at least kkion (number of ions + and -)
c         icode  - code for selecting current(1) or voltage(0) control.
c         xac    - ac amplitude of total current or voltage (A/m2 or V).
c         xdc    - dc value of total current or voltage (A/m2 or V).
c         iprnt  - iprnt=1 gives local print, otherwise not.
c         lun    - output unit for writing results
c
c Output: time() - time array used for output arrays
c         vw()   - one rf cycle of wall potential (V)
c         vwb()  - one rf cycle of ion potential (V)
c         ew()   - one rf cycle of field at wall (V/m)
c         aj()   - one rf cycle of total current (A/m2)
c         vaveb  - the time average ion potential (V) at wall
c         ierr   - error status code. 0 is normal, 1 and 2 abnormal.
c
c Calls:  eofv, dedv, dedvb, dvdt, dvbdt
c
      common/phycon/pi,echar,eps0
      common/params/an0,etemp,tave,am2inc,amsq
      dimension time(*),vw(*),vwb(*),ew(*),aj(*)
      dimension amass(*),ank(*), kchg(*), kion(*)
      logical lstor
      ierr=0
c
c Basic physical constants. All quantities are in SI units.
c      data pi,emass,echar,eps0,akbltz,amu/3.1415926535898,9.109534e-31,
c     1 1.6021892e-19,8.85418782e-12,1.380662e-23,1.6605665e-27/
      pi = 3.1415926535898
      emass =9.109534e-31
      echar =1.6021892e-19
      eps0  =8.85418782e-12
      akbltz=1.380662e-23
      amu   =1.6605665e-27
c
c      open (unit=lun,file='last.s',status='unknown')
c
c Some quantities related to input.
c
c The radian frequency.
      omrf=2.0*pi*anurf
c The period.
      rfperd=1.0/anurf
c
c Relocate density and temperature to avoid "common exclusion."
      etemp=etev
c
c Loop over ions to define quantities
c
      an0 = 0.0
      totvi = 0.0
      omegav = 0.0
      ajitot = 0.0
      sumnip = 0.0
      avgmss = 0.0
      zeroe = 0.0
      do 900 ki = 1, kkion
         k = kion(ki)
c
c The electron density from quasineutrality
c
         ani = max(zeroe,ank(k))
         an0 = an0 + kchg(k)*ani
c
c The ion plasma radian frequency.
c
         omegip = dsqrt(ani*echar*echar/eps0/(amass(k)*amu))
         omegav = omegav + omegip*ani
c
c The Bohm velocity and the corresponding ion current density
c
         if (kchg(k).gt.0) then
            vbohm = dsqrt(echar*etemp/(amass(k)*amu))
            aji = echar*kchg(k)*ani*vbohm
            totvi = totvi + vbohm*ani*kchg(k)
            ajitot = ajitot + aji
            avgmss = avgmss + amass(k)*ani
            sumnip = sumnip + ani
         else
            vbohm=0.0
         endif
c
 900  continue
      omegav = omegav/sumnip
      vbavg = totvi/sumnip
      avgmss = (avgmss/sumnip)
c
      if(iprnt.eq.1) write(lun,1000) anurf,avgmss,etemp,an0
 1000 format(' ',//,1x,'rf freq(Hz)',1x,'avg ion wt.',3x,'Te(eV)',7x,
     &        'ne(/m3):'/,1p,4e12.4)
c
c Derived quantities
c
c The electron thermal velocity
      usube=dsqrt(8.0*etemp*echar/pi/emass)
      if (icode.eq.1) xdc=echar*an0*0.9*vbavg
c The ion Mach number squared, initial setting.
      amsq=1.0
c The increment to be used if first integral, E(V), goes imaginary.
      am2inc=0.01
c The Debye length
      debyel=dsqrt(eps0*etemp/echar/an0)
c The floating potential, negative.
      vfloat=-etemp*dlog(0.25*usube/vbavg)
c The damping constant based on radian ion plasma frequency.
      tave=1.0/omegav
c
c Print out the derived quantities.
c
      if(iprnt.eq.1) write(lun,1001) usube,vbavg,debyel,omegav
 1001 format(' ',//,1x,'C_e(m/s)',2x,'Ubohm_avg(m/s)',2x,
     &       'Debye(m)',2x,'avg ion omega:',/,1p,4e12.4)
      if(iprnt.eq.1) write(lun,1002) vfloat,tave,ajitot
 1002 format(' ',//,1x,'Vfloat(V)',3x,'t_damp(s)',2x,
     &       'Ji_tot(A/m2):',/,1p,4e12.4)
c
c Set control current or voltage.
c
      if(icode.eq.1) then
c The amplitude of the ac current in A/m2. The dc current also.
         ajac=xac
         ajdc=xdc
         if(iprnt.eq.1) write(lun,1003) ajac,ajdc
 1003    format(' ',//,' current control, J_ac and J_dc are:',
     &          1p,2e12.4/)
      else
c The sheath voltage, ac and dc
         vac=xac
         vdc=xdc
         if(iprnt.eq.1) write(lun,1004) vac,vdc
 1004    format(' ',//,' voltage control, V_ac and V_dc are:',
     &          1p,2e12.4/)
         if(vdc.ge.-abs(vac)) then
            ierr=1
            write(lun,1008)
 1008       format(//,' ','STOPPING IN SHEATH',/
     &             ' - transiently nonnegative wall potential.')
            return
         endif
      endif
c
c Prepare for the time-dependent integration. Numerics and counting.
c These values are set for adequate accuracy under nominal conditions.
c Any drastic changes may require revisions of step size and logic.
c
c Set the total time of integration in terms of total rf cycles
      ncycle=20
c The time increment fixed as number of steps per rf cycle.
      inccyc=1000
      delt=rfperd/float(inccyc)
c Set print interval in intervals per rf cycle (must divide evenly)
      nptcyc=100
c Total number of steps not counting zero.
      ntotal=ncycle*inccyc
c Steps per print interval, inccyc should be a multiple of nptcyc
      iperp=inccyc/nptcyc
c First point to be stored for return of time-dependent variables.
      istart=(ncycle-1)*inccyc
c
c Initial values of time, potential, average potential.
      t=0.0
      vwall=vfloat
      vwbar=vfloat
c
c Initialize counter for storage of dependent variables.
      istor=0
c Initialize average potential for ions.
      vaveb=0.0
c
      if(iprnt.eq.1) write(lun,1005)
 1005 format(2x,'icount',6x,'t(ns)',7x,'Vwall',7x,'Vbarw',4x,
     &       'Je(A/m2)',2x,'Jd(A/m2)',3x,'Jtot(A/m2)',6x,'M^2')
c
c The integration loop.
c
      do 200 icount=1,ntotal
c Derivatives for current or voltage control
         if(icode.eq.1) then
            ajtot=ajac*sin(omrf*t)+ajdc
            ewall=eofv(vwall,vwbar,ierr)
            aje=-echar*0.25*usube*an0*dexp(vwall/etemp)
            ebyv=dedv(vwall,vwbar,ewall)
            ebyvb=dedvb(vwall,vwbar,ewall)
            vbbyt=dvbdt(vwall,vwbar)
            ajd = ajtot - aje - ajitot
            vbyt=dvdt(ebyv,ajd,ebyvb,vbbyt)
           vwall=vwall+delt*vbyt
         else
            vwall=vac*sin(omrf*t)+vdc
            aje=-echar*0.25*usube*an0*dexp(vwall/etemp)
            ewall=eofv(vwall,vwbar,ierr)
            ebyv=dedv(vwall,vwbar,ewall)
            ebyvb=dedvb(vwall,vwbar,ewall)
            vbyt=vac*omrf*cos(omrf*t)
            vbbyt=dvbdt(vwall,vwbar)
            ajtot=aje + eps0*(vbyt*ebyv+vbbyt*ebyvb) + ajitot
         endif
         if(ierr.ne.0) return
c Increment the average solution and time
         vwbar=vwbar+delt*vbbyt
         t=t+delt
c Prepare storage grid variables and average of ion potential.
         lstor=(icount.ge.istart).and.((icount/iperp)*iperp.eq.icount)
         if(lstor) then
            istor=istor+1
            time(istor)=t
            vw(istor)=vwall
            vwb(istor)=vwbar
            ew(istor)=ewall
            aj(istor)=ajtot
            if(iprnt.eq.1) then
              ajd = ajtot - aje - ajitot
              tns = (icount-istart)*delt*1.0e9
              write(lun,1007) icount,tns,vwall,vwbar,aje,ajd,ajtot,amsq
            endif
 1007       format(i8,2x,e12.4,1x,6e11.3)
            vaveb=vaveb+vwbar
         endif
c
  200 continue
c
c Average the ion potential
      vaveb=vaveb/istor
      if(iprnt.eq.1) write(lun,1006) vaveb
 1006 format(' The average ion potential is ',1p,e12.4)
      close (unit=lun)
      return
      end

      function dedv(vwall,vwbar,ewall)
      implicit real*8(a-h,o-z)
c Evaluates dEwall/dVwall from Vwall, Vwall average, and Ewall.
      common/phycon/pi,echar,eps0
      common/params/an0,etemp,tave,am2inc,amsq
      chi=vwall/etemp
      chibar=vwbar/etemp/amsq
      pref=echar*an0/eps0
      dedv=1.0/ewall*pref*((dsqrt(1.0-2.0*chibar)-1.0)/chibar+dexp(chi))
      return
      end

      function dedvb(vwall,vwbar,ewall)
      implicit real*8(a-h,o-z)
c Evaluates dEwall/dVwbar from Vwall, Vwall average, and Ewall.
      common/phycon/pi,echar,eps0
      common/params/an0,etemp,tave,am2inc,amsq
      chi=vwall/etemp
      chibar=vwbar/etemp/amsq
      pref=echar*an0/eps0/amsq
      root=dsqrt(1.0-2.0*chibar)
      dedvb=1.0/ewall*pref*(1.0-root-chibar/root)
     1      *chi/chibar/(chibar*etemp)
      return
      end

      function eofv(vwall,vwbar,ierr)
      implicit real*8(a-h,o-z)
c Evaluates Ewall from Vwall and Vwall average.
      common/phycon/pi,echar,eps0
      common/params/an0,etemp,tave,am2inc,amsq
      istop=0
      chi=vwall/etemp
    1 chibar=vwbar/etemp/amsq
      istop=istop+1
c Stop possible trapping in infinite loop.
      if(istop.ge.10000) then
        ierr=2
        eofv=0.0
        return
      endif
      pref=echar*an0/eps0*etemp
      arg=2.0*pref*(chi/chibar*(dsqrt(1.0-2.0*chibar)-1.0)
     1         +dexp(chi)-1.0)
c Protect against imaginary field by incrementing ion injection 
c velocity.
      if(arg.le.0.0) then
         amsq=amsq+am2inc
         go to 1
      endif
      eofv=dsqrt(arg+1.0e-12)
      return
      end

      function dvdt(ebyv,ajd,ebyvb,vbbyt)
      implicit real*8(a-h,o-z)
c Evaluates time derivative of wall potential from dEwall/dVwall and
c  the currents. Also needs the Vbar terms.
      common/phycon/pi,echar,eps0
      dvdt=(1.0/eps0*(ajd)-ebyvb*vbbyt)/ebyv
      return
      end

      function dvbdt(vwall,vwbar)
      implicit real*8(a-h,o-z)
c Evaluates time derivative of damped wall potential.
      common/params/an0,etemp,tave,am2inc,amsq
      dvbdt=-(vwbar-vwall)/tave
      return
      end
