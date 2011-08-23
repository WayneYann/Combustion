C    CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
C///////////////////////////////////////////////////////////////////////
C
C            CRESLAF: CHEMICALLY-REACTING SHEAR LAYER FLOW
C                (TWO-DIMENSIONAL BOUNDARY-LAYER MODEL)
C
C     WRITTEN BY:
C         MICHAEL E. COLTRIN
C         SANDIA NATIONAL LABORATORIES
C         SURFACE PROCESSING SCIENCES DIVISION
C         ALBUQUERQUE, NM 87185
C         (505) 844-7843
C       AND
C         ROBERT J. KEE
C         COMPUTATIONAL MECHANICS DIVISION
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94550
C         (415) 294-3272
C
C///////////////////////////////////////////////////////////////////////
C
C     V.5.10 98/03/04 (E. Meeks)
C     1. Action #142: In CRTWO's call to TWOPNT, change argument name 
C        LDIM1 to LDUM1 and declare LDUM1 as logical.
C     V.5.9 97/10/29 (E. Meeks)
C     1. Fix bug#096b: Set ABSOL=RELAT=SQRT(D1MACH(4)), rather than
C        calculating the machine constant.
C     2. Put double-precision change blocks around above, and
C        add single-precision change block using R1MACH(4).
C     3. Add explicit typing of D1MACH as DOUBLE PREC and R1MACH as REAL
C        in CRTWO within appropriate change blocks.
C     4. Remove U, COMP from type declaration statement in CRTWO.
C     V.5.8 97/09/19
C     1. Three places in CRINIT, set S(NU,1)=S(NU,2) for the case
C        of radial coordinates.
C     2. Reformat the creslaf banner statment into separate lines.
C     3. Remove a write to unit 8 that was not needed (in CRRUN).
C     4. Remove un-used variables: PI (CRFUN); SMALL,KPERLN (CRPRNT);
C        VTRACE (CRRDKY); KPERLN,KPRLN, SMALL (CRSPRT).
C     V.5.7 97/08/27
C     1. Fix bug#078: Add test for ICRD.EQ.0 for LONGPT option in 
C        CRPRNT.
C     V.5.6 97/08/27
C     Fix bugs reported in Action #079:
C     1. add KERR as last argument to second CALL CRTWO
C     2. delete XMESH from CALL CRTWO and MESH from SUBROUTINE CRTWO
C        arguments (obsolete TWOPNT requirement)
C     3. delete L(LAC) from CALL CRTWO and ACTIVE from SUBROUTINE CRTWO
C        arguments (obsolete TWOPNT requirement), and delete LAC pointer
C        from COMMON /CRCRCR/ and LAC assignment SUBROUTINE CRPNT
C     4. delete L(LMK) from CALL CRTWO and MARK from SUBROUTINE CRTWO
C        arguments (obsolete TWOPNT requirement), and delete LMK pointer
C        from COMMON /CRCRCR/ and LMK assignment SUBROUTINE CRPNT
C        (3 and 4 to result in LTOT=1 in SUBROUTINE CRPNT)
C     5. delete argument SITDOT from CALL SKRATI
C     6. delete argument IFLAG from CALL SKINDX, and delete error check
C        of IFLAG immediately following the call
C     V.5.5 97/08/18
C     1. Fix bug #075: change error subroutines from single to 
C        double precision for consistent change blocks.
C     2. Change printout to read "Chemkin-III" rather than "II".
C     V.5.4 97/07/23 
C     1. In CRINIT, line 1391, remove space from logical operator.
C        (additional occurrence of Bug #012)  (F. Rupley)
C     V.5.3 97/07/22 (F. Rupley)
C     1. in CRRUN, change DATA INFO/15*0.0/ to DO loop
C        initialization of INFO(I), I=1,15.
C     V.5.2 97/07/15 (F. Rupley)
C     1. remove DATA LIN/5/ from CRRDKY and add to the CALL list
C     2. correct READ of dassl workspace in CRRSTR
C     V.5.1 97/06/13 (M. Coltrin)
C     1. Convert ". AND." to " .AND." (and similarly for "OR", "NE", 
C        and "EQ") to conform to ANSI standard (Chemkin Action #012).
C     2. Removed the variable SUM from CRSFUN, and the variable LDUM1
C        from CRTWO, because they were not ever used.
C     3. Put if-test protection around loops 6280 and 6300 in CRRDKY to
C        account for the case of KKSURF or KKBULK = 0 (Action #042).
C     4. In CRSPRT, two IF tests issued a RETURN if KKSURF or KKBULK
C        were zero. We really didn't want to return, but rather just
C        skip the next block of code. Logic was corrected in this 
C        version.
C     V.5.0 97/04/15
C     1. fix indexing DO 34 loop of SUBROUTINE CRRUN;
C        CHEMKIN Committee bug report #002
C     V.4.9 97/03/21
C     1. update DASERR copy of COMMON /CRCRCR/ to conform to that in
C        other subroutines
C     VERSION 4.8 97/03/01
C     1. move OPEN statements to main "driver" program,
C        to satisfy Reaction Design requirement of providing object
C        files instead of source code.
C     VERSION 4.7, 96/05/24
C     1. initial sccs version
C     2. reverse chronology of change history
C     VERSION 4.6 (4/8/96 F. Rupley)
C     1. incorporate twopnt v.3.22 - call list changed
C     VERSION 4.5 (2/27/95 F. Rupley)
C     1. Change character index "(:" to "(1:"
C     VERSION 4.4 (2/17/95 F. Rupley)
C     1. Error omission of IFLAG in CALL SKINIT corrected.
C     VERSION 4.3 (1/19/95 F. Rupley)
C     1. Add integer error flag variable to CKLEN, CKINIT, MCLEN,
C        MCINIT and SKLEN, SKINIT call lists.
C     VERSION 4.2 (9/29/94 F. Rupley)
C     1. 'IF (KEYWRD' checking modified for keywords less than 
C        4 characters long.
C     VERSION 4.1 (7/25/94 F. Rupley per M. Coltrin)
C     1. Remove C from column1, DO 215 loop, SUBROUTINE CRES
C     VERSION 4.00 (2/2/93) Coltrin, Grcar and Rupley
C     1. Install new twopnt with changes to CRPNT and CRTWO.
C     2. Remove do 2111 loop (some diagnostics print statement)
C     VERSION 3.96 (2/19/92 M. Coltrin)
C     1. Add subroutine CRPROF to read a solution profile.
C     2. Add PROF keyword and associated LPROF logical variable.
C     3. Add unit number LPRO to CRES arguements.
C     4. Change RPAR(NZ1) to RPAR(NZJJ) in first call to CRTWO. 
C     5. Changed logic having to do with LNOTP. 
C     VERSION 3.95 (2/15/92 M. Coltrin)
C     1. Add temperature print-out to CRSPRT
C     VERSION 3.94 (12/9/91 M. Coltrin)
C     1. Convert input temperatures (keywords TSPL, GTMP, STMP) from
C        Centegrade to Kelvin.
C     2. Change keyword input for ICRD from numeric (0, 1, or 2) to
C        characters (PLAN, RAD, or SYMC).
C     3. Change definition of VEL keyword to be the maximum velocity in
C        the channel (it used to be the average velocity). Therefore,
C        in computing an initial parabolic velocity profile we no
C        longer need to multiply the input velocity by 1.5 (planar)
C        or 2.0 (radial).
C     4. Add option for a correction velocity (keyword VCOR). This 
C        involved adding integer variable IVCOR, logical variable LVCOR,
C        a pointer IIVC, and loopS 700 and 800 in CRDIFV, and modifying
C        CRFUN and CRSFUN species equations. 
C     VERSION 3.93 (10/17/91 M. Coltrin)
C     1. Modify CRSFUN to be valid for either upper or lower boundary.
C     2. Call CRTWO for both the upper and lower boundaries in 
C        cartesian coordinates.
C     VERSION 3.92 (10/14/91 M. Coltrin)
C     1. Add keyword GASW for initial guess of gas-phase mole fractions
C        at the wall (for initial TWOPNT problem).
C     2. Add keyword STP0 for initial time-step in TWOPNT problem.
C     3. Scale surface residual by SDEN(N) in CRSFUN.
C     VERSION 3.91 (10/3/91 M. Coltrin)
C     1. When a boundary-layer thickness is given (BLTK), put in a
C        linear temperature gradient between the surface temperature
C        and the gas-temperature over that thickness range.
C     2. Calculate RHO at every node in CRRTOX.
C     VERSION 3.90 (10/2/91 M. Coltrin)
C     1. Add option for a symmetric channel in cartesian coordinates
C        (ICRD=2).
C     VERSION 3.84 (9/28/91 M. Coltrin)
C     1. Add IFIXED to arguements of CRINIT.
C     2. Set T at top wall to TEMP if IFIXED=2 in CRINIT.
C     3. Correct equation for parabolic velocity profile at lower wall
C        in cartesian coordinates if LTHIC is true.
C     VERSION 3.83 (9/28/91 M. Coltrin)
C     1. Comment-out false transient for gas-phase species in CRSFUN.
C     VERSION 3.82 (9/24/91 H. Moffat)
C     1. Added another keyword, 'PRND', to control the printing
C        of diagnostic information from DASSL.
C     2. Protected taking the square root of the radius operation.
C        Added absolute values inside the square root (F. Rupley)
C     VERSION 3.81 (9/16/91 H. Moffat)
C     1. Fixed an error in the pointers for GFAC and SFAC.
C     2. Added another keyword, 'HO', to set the initial time step.
C     3. Added a symmetric temperature profile capability for
C        cartesian coordinates.  The top wall temperature can
C        now be specified to be equal to the bottom wall temperature.
C        The new keyword for this is SYMT.
C     VERSION 3.80 (9/9/91 F. Rupley per Bob Kee)
C     1. Expand LTHIC logic to include planar coordinates.
C     VERSION 3.79 (6/27/91 F. Rupley per Bob Kee)
C     1. Implement boundary layer thickness in CRINIT - new keyword
C        BLTK, logical LTHIC
C     VERSION 3.78 (6/27/91 F. Rupley per Bob Kee)
C     1. Change the pointer for LAC in CRPNT from length KKSURF to
C        KKSURF+KKGAS
C     2. Implement SFAC and GFAC scaling of rates.
C     3. Keyword NOTP to skip calling TWOPNT
C     4. Output in mole fraction rather than partial pressure.
C     VERSION 3.77
C     1. Took out work around for the sun compiler
C     2. Minor changes to write statements so that it will compile
C        on vax the vax ultrix fortran compiler, f77.
C     VERSION 3.76
C     1. Work around for sun compiler bug that causes bombs on long
C        unformatted writes.
C     VERSION 3.75
C     1. Fixed error in the dimensioning of the variable, ACTIVE.  It
C        should be COMPS long instead of KKSURF long.
C     VERSION 3.74
C     1. Add TWOPNT keyword DTMN
C     2. TWOPNT keywords ISTP, IRET, NJAC, TJAC, and DTMX to be added
C        in version 3.7 were actually missing, so they were put-in
C     VERSION 3.73
C     1. Add TWOPNT keywords: TWAB, TWRE, TWTA, TWTR, TWPR, TWST
C     2. Dynamically dimension XLINPK arrays: AA, RTEMP, DR, DC, XTEMP
C     VERSION 3.72
C     1. Only evaluate transport properties and diffusion velocities
C        at one node in CRSFUN
C     VERSION 3.71
C     1. Get rid of variables STIM and NSTM (no longer needed)
C     2. Hold transport properties fixed in evaluation of jacobian
C        for TWOPNT
C     VERSION 3.7
C     1. Do not use linking file information stored in a restart
C        solution, as mechanism may have changed.
C     2. Call list for TWOPNT requires additional input; optional
C        use of new keywords reset default values:
C        'ISTP' n - sets NINIT initial time steps before newton
C                   (default is 0)
C        'IRET' n - set retirement age IRETIR of old time step
C                   (default 50)
C        'NJAC' n - set retirement age NJAC of Jacobian during
C                   steady state newton (default 20)
C        'TJAC' n - set retirement age ITJAC of Jacobian during
C                   time stepping (default 20)
C        'DTMX' x - set maximum time step DTMAX (default 1.0E-4)
C     VERSION 3.62
C       1. Changed the twopoint initial value problem so that the
C          algebraic constraints (the surface site fraction
C          equations) are satisfied at t = 0+
C       2. Added change blocks for linear algegra subroutines
C          that scale the linear algebra problem before
C          factorization and optionally conducts iterative
C          improvement on the linear algebra problem.
C       3. Scale surface species residuals by the site densities
C          in CRFUN
C       4. Add individual time step error controls for each
C          species.
C     VERSION 3.61
C       1. Check that sites, site species, or bulk species exist
C          before filling appropriate arrays.
C       2. SDOT(*) becomes SITDOT(*), WDOT(*) becomes SDOT(*);
C          WDOTI(*) becomes SDOTI(*), and addition array SITDTI(*)
C          is returned from SUBROUTINE SKRATI.
C     VERSION 3.6
C       1. CALL SKSNUM instead of CKSNUM for surface species in
C          CRRDKY
C       2. CALL SKRATI instead of SKRTI
C     VERSION 3.5
C       1. Upgrade to Surface Chemistry Version 3.5
C       2. Add KKBULK and ACT to data file.
C     VERSION 1.9
C       1. Add keyword GRAV and calculate RHOE in CRINIT
C       2. Add gravity term to momentum equation in CRFUN
C     VERSION 1.8
C       1. Additional arguments read from linking file are KERR and
C          IIREV
C       2. VAX open statements changed; will no longer use LINKTP$,
C          LINKCK$, and LINKSK$
C     VERSION 1.7
C       1. Replace initial surface time-integration with TWOPNT
C          iteration
C       2. Get initial temperature from calling CRWALL rather than
C          GTMP to avoid small difference due to single/double
C          precision conversion
C       3. Set XSPSTR and XSPEND for both IDOSP=0 and 1 in CRTMPR
C     VERSION 1.6
C       1. In CRPNT fix read of LINK file from SKLIB.18
C       2. In CRPNT fix read of LINK file from CKLIB.14
C       3. Change calls from SKRAT to SKRAT3
C       4. Change calls from SKRATI to SKRTI3
C       5. Change arguements in calls to SKRNAM
C       6. Add keyword XTMP and make changes to CRTMPR and CRRDKY
C     VERSION 1.5
C       1. Add TSPL as keyword input and necessary code to
C          allow the user to input wall temperatures
C       2. Read LINK file from SKLIB.17, which includes the
C          thermodynamic properties.  This change only requires adding
C          three variables to the READ (LINKSK) XXXXXX  statement.
C       3. Scale residual F(NR,1) for radial coordinates by 1E18
C          to try to avoid negative square root
C       4. Switch order of arguements in calls to SKKWT
C     VERSION 1.4
C       1. Let user have the option of specifying MIX and TDIF
C       2. Add call to MCMCDT in mixture average section of CRTRNP
C     VERSION 1.3
C       1. Add "quad precision"
C       2. Add MORD and KKSURF to restart file
C       3. Add call to CRSAVE at bottom of DASSL integration loop
C     VERSION 1.2
C       1. Add MORD and STIM as keyword inputs
C       2. Take MORD out of subroutine CRSRUN
C       3. Integrate rate equation for KKSURF species in CRSFUN
C       4. Force non-negative surface fractions in CRSRUN
C       5. Change units (scale) surface residuals in CRFUN
C       6. Add print out at bottom of DASSL integration loop
C       7. Took out the non-negative constraint in CRSRUN
C       8. Add option of negative STIM in CRRUN
C     VERSION 1.1
C       1. Replace TWOPNT iteration with time integ. of surface species
C       2. Change initial time derivative of mass-loss equations in
C          loop 300 of CRRUN
C       3. New version of subroutine CRRTOX
C       4. Change CWORK(JCCH) to CWORK in call to CKINIT
C       5. Take P out of the call list for CRTRNP
C       6. Delete 10**30 scaling of F(NR,1) in CRFUN
C     VERSION 1.0
C       1. Add TWOPNT iteration for self-consistent surface fractions
C       2. Correct KKSURF-surface species conservation equation
C       3. Change sign in DMDX for the upper boundary in CRFUN
C       4. Fix errors in momentum, energy, and species eqs. in CRFUN
C       5. Add KKSURF to arguments of CRRSTR
C       6. Add ICRD to arguments of CRSAVE and CRTIM
C       7. Fix limits on DO 200 loop in CRRUN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRES (LIN, LOUT, LINKCK, LINKSK, LINKMC, LREST, LSAVE,
     1                 LDAS, LPRO, LNIPAR, IPAR, LNRPAR, RPAR, LNCWRK, 
     2                 CWORK, LNLWRK, LWORK)
C
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION IPAR(*), RPAR(*)
      CHARACTER CWORK(*)*(*)
      LOGICAL LWORK(*), LTDIF, LNOTP, LMOLF, LTHIC, LVCOR, LUPPER, LPROF
      LOGICAL KERR
C
      COMMON /CRCRCR/ NCK,  NMC,  NSK,  NWT,  NMAS, NSLB, NSRB, NX,
     1                NZ1,  NZJJ, NS,   NSP,  NYAV, NYV,  NHH,  NCON,
     2                NVIS, NCND, NRR,  NY,   NYP,  NROP, NXAV, NXM,
     3                NXP,  ND,   NTDC, NF,   NDEP, NEKR, NS1,  NS2,
     4                NQ1,  NQ2,  NQ12, NWDT, NSDT, NDAS, NCP,  NSPS,
     5                NSPE, NSTS, NDDE, NSDE, NXTW, NFZ,  NZN,  NFZN,
     6                NWDI, NSDI, NGTP, NPR,  NTWP, NABV, NBLW, NBUF,
     7                NA,   NACT, NAA,  NRTP, NDR,  NDC,  NXTP, NASP,
     8                NRSP, NATS, NRTS, NGFC, NSFC, NGASW,ICK,  IMC,
     9                ISK,  IJJ,  INAJ, IKK,  IKKT, IICR, IITD, IIML,
     *                ILRO, ILRI, ILDS, ILOU, IIFX, IPKK, IPKF, IPKL,
     1                IIDS, ISIT, IFS,  ILS,  IBLK, IFB,  ILB,  IKKS,
     2                ITPW, IIPV, IMO,  IIVC, JCCH, JSCH, JKCH, JPCH
      COMMON /LOCS/  NP, NR, NU, NT, NML, NMU, NYS, NY1
      COMMON /CONSTS/RUC, PATM, GRAV, RHOE
C
      DATA JJ/4/
      KERR = .FALSE.
C
C*****precision > quad
C      PI=QATAN2(0.0Q0,-1.0Q0)
C*****END precision > quad
C*****precision > double
      PI=DATAN2(0.0D0,-1.0D0)
C*****END precision > double
C*****precision > single
C      PI=ATAN2(0.,-1.)
C*****END precision > single

C     initialize LOCS common block
      NP = 1
      NR = 2
      NU = 3
      NT = 4
      NML = 5
      NMU = 6
      NYS = 6
      NY1 = 7
C
C     write version number
      WRITE (LOUT,*) '   '
      WRITE (LOUT,*)' CRESLAF: CHEMICALLY-REACTING SHEAR LAYER FLOW '
      WRITE (LOUT,*)'    (TWO-DIMENSIONAL BOUNDARY-LAYER MODEL) '
      WRITE (LOUT,*)'    (CHEMKIN-III Version 5.10, 98/03/04)'
      WRITE (LOUT,*)'   '
C
C*****precision > quad
C      WRITE (LOUT,*)'     QUADRUPLE PRECISION'
C*****END precision > quad
C*****precision > double
      WRITE (LOUT,*)'     DOUBLE PRECISION'
C*****END precision > double
C*****precision > single
C      WRITE (LOUT,*)'     SINGLE PRECISION'
C*****END precision > single
C
C
      IPCALL = 1
      ICRD = 0
      DO 25 I = 1, 2
         IPCALL = I
         CALL CRPNT (LINKCK, LINKMC, LINKSK, LOUT, NATJ, LNIPAR, LNRPAR,
     1               LNCWRK, KKGAS, KKSURF, KKBULK, KKTOT, NPHASE,
     2               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK, II,
     3               IISUR, LNTWOP, IPCALL, JJ, ICRD, LNLWRK, IPAR,
     4               RPAR, CWORK, KERR)
         IF (KERR) RETURN
C
         IF (I .EQ. 1)
     1       CALL CRRDKY (LIN, KKGAS, KKSURF, KKBULK, KKTOT, NFSURF, 
     2                    NLSURF,NFBULK,NLBULK,IPAR(IPKF), IPAR(IPKL),
     3                    IPAR(IPKK), NPHASE, CWORK(JPCH), LOUT,
     4                    RPAR(NYP), PRES, VEL, GTMP, STMP, XEND, JJ,
     5                    ICRD, HITE, DX, IRST, ITDIF, IMULT, IFIXED,
     6                    STCH, ATOL, RTOL, CWORK(JKCH), RPAR(NZ1),
     7                    RPAR(NZJJ), MORD, IDOSP, XTMP, RPAR(NACT),
     8                    GFAC, SFAC, RPAR(NASP), RPAR(NRSP), SSABS,
     9                    SSREL, TDABS, TDREL, IPRNT, LIMIT, NINIT,
     +                    IRETIR, NJAC, ITJAC, DTMAX, DTMIN, LNOTP,
     1                    LMOLF, LTHIC, BLTK, HO, NPRNDS, RPAR(NGASW),
     2                    TSTEP0, IVCOR, LPROF, KERR)
      IF (KERR) RETURN
C
   25 CONTINUE
C
      RPAR(NGFC) = GFAC
      RPAR(NSFC) = SFAC
      RPAR(NPR)  = PRES
      RPAR(NGTP) = GTMP
C
      IF (IRST .EQ. 1) THEN
C
C         restarting the program
          CALL CRRSTR (KKGAS, JJ, NATJ, LREST, LSAVE, LDAS, TIME,
     1                 RPAR(NS), RPAR(NSP), RPAR(NSLB), RPAR(NSRB),
     2                 RPAR(NRR), RPAR(NX), IPAR(IIDS), RPAR(NDAS),
     3                 ICRD, RPAR(NDEP), GTMP, STMP, RPAR(NMAS),
     4                 ITDIF, IMULT, IFIXED, KKSURF, MORD, RPAR(NACT),
     5                 IVCOR)
C
      ENDIF
C
C     put integer flags into IPAR
      IPAR (INAJ)   = NATJ
      IPAR (IJJ)    = JJ
      IPAR (IKK)    = KKGAS
      IPAR (IKKS)   = KKSURF
      IPAR (IKKT)   = KKTOT
      IPAR (ISIT)   = NNSURF
      IPAR (IFS)    = NFSURF
      IPAR (ILS)    = NLSURF
      IPAR (IBLK)   = NNBULK
      IPAR (IFB)    = NFBULK
      IPAR (ILB)    = NLBULK
      IPAR (IICR)   = ICRD
      IPAR (IITD)   = ITDIF
      IPAR (IIVC)   = IVCOR
      IPAR (IIML)   = IMULT
      IPAR (IIFX)   = IFIXED
      IPAR (ILRO)   = LSAVE
      IPAR (ILRI)   = LREST
      IPAR (ILDS)   = LDAS
      IPAR (ILOU)   = LOUT
      IPAR (IMO)    = MORD
C
C     set up splines for susceptor temperature
      CALL CRTMPR (GTMP, STMP, RPAR(NSPS), RPAR(NSPE), RPAR(NSTS),
     1             IDOSP, XTMP)
C
C     IPTGAS points to the location in the solution vector of the
C     first gas-phase node
C
      IF (IRST .EQ. 0) THEN
         IF (ICRD .EQ. 0) THEN
C           planar coordinates
            IPTGAS = KKSURF
         ELSE
C           radial coordinates (or planar with a symmetry axis)
            IPTGAS = 0
         ENDIF
C
C        put the initial conditions into the S array
         CALL CRINIT (KKGAS, RPAR(NWT), IPAR, RPAR, PI, NATJ, ICRD, JJ,
     1                GTMP, PRES, VEL, RPAR(NYP), HITE, STCH,
     2                RPAR(NS+IPTGAS), RPAR(NRR), RPAR(NX), RPAR(NSLB),
     3                RPAR(NSRB), RPAR(NMAS), RPAR(NY), RPAR(NSPS),
     4                RPAR(NSPE), LTHIC, BLTK, IFIXED, RPAR(NGASW))
C
         IF (LPROF) THEN
C           read user-supplied solution profile from unit LPRO
            CALL CRPROF (LPRO, LOUT, RPAR(NS+IPTGAS), RPAR(NRR), 
     1                   RPAR(NX),   RPAR(NACT),      RPAR(NZ1), 
     2                   RPAR(NZJJ), IPAR, RPAR, NATJ, JJ, 
     3                   KKGAS, KKBULK, KKSURF, ICRD,  RPAR(NMAS),
     4                   KERR )
            IF (KERR) RETURN
         ENDIF
C
         IF (ICRD .EQ. 0) THEN
           IPTGAS = KKSURF
         ELSE
           IPTGAS = 0
         ENDIF

         IF (KKSURF.GT.0 .OR. IISUR.GT.0) THEN
             IF (LNOTP) THEN
C               do not solve the TWOPNT problem for the initial guess
C
                IPTTOP = IPTGAS + JJ*NATJ
                DO 210 K = 1, KKSURF
                   RPAR(NS+IPTTOP + K - 1) = RPAR(NZJJ + K - 1)
  210           CONTINUE
                IF ( ICRD .EQ. 0) THEN
                   DO 215 K = 1, KKSURF
                      RPAR(NS+K-1) = RPAR(NZ1 + K -1)
  215              CONTINUE
                ENDIF
C
             ELSE
C               solve TWOPNT problem to get initial guess
C
                IPMAX = 1
                LTDIF = (IPAR(IITD) .EQ. 1)
                LVCOR = (IPAR(IIVC) .EQ. 1)
                LUPPER = .TRUE.
C   
                CALL CRTWO (KKSURF+KKGAS, LOUT, IPMAX, LNTWOP, 
     1           RPAR(NTWP),
     1           RPAR(NABV), RPAR(NBLW), RPAR(NBUF),
     2           RPAR(NXTW), RPAR(NA), IPAR(ITPW),
     3           IPAR, RPAR, IPAR(ISK), RPAR(NSK), PRES,
     4           RPAR(NZJJ), RPAR(NYP), KKGAS, KKSURF, KKBULK,
     5           IPAR(ISIT), IPAR(IBLK), IPAR(IPKK), IPAR(IPKF),
     6           IPAR(IPKL), RPAR(NWDT), RPAR(NSDT), RPAR(NFZ),
     7           RPAR(NZN), RPAR(NFZN), IPAR(IIPV), RPAR(NWT),
     8           RPAR(NDDE), RPAR(NSDE), RPAR(NACT), IPAR(IFS),
     9           IPAR(ILS), IPAR(IFB), IPAR(ILB), CWORK(JKCH),
     *           CWORK(JPCH), PATM, IISUR, RPAR(NROP), CWORK(JSCH),
     1           RPAR(NYAV), RPAR(NS+IPTGAS), JJ, NATJ, LTDIF,
     2           IPAR(IIML), RPAR(NX), RPAR(NMAS), IPAR(IICR),
     3           IPAR(IMC), RPAR(NMC), RPAR(NXAV), RPAR(NHH),
     4           RPAR(NVIS), RPAR(NCND), RPAR(ND), RPAR(NTDC),
     5           RPAR(NYV), RPAR(NAA), RPAR(NRTP), RPAR(NDR), RPAR(NDC),
     6           RPAR(NXTP), SSABS, SSREL, TDABS, TDREL, IPRNT, LIMIT,
     7           NINIT, IRETIR, NJAC, ITJAC, DTMAX, DTMIN, GFAC, SFAC,
     8           TSTEP0, LUPPER, LVCOR, KERR)
                IF (KERR) RETURN
C
C               Store solution from TWOPNT in the solution vector
                IPTTOP = IPTGAS + JJ*NATJ
                DO 225 K = 1, KKSURF
                   RPAR(NS+IPTTOP + K - 1) = RPAR(NXTW + KKGAS + K - 1)
  225           CONTINUE
                DO 230 K = 1 , KKGAS
                   RPAR(NS + IPTTOP + NYS - NATJ + K - 1) =
     1                            RPAR(NXTW + K - 1)
  230           CONTINUE
C
                IF (ICRD .EQ. 0) THEN
C                  planar coordinates
                   LUPPER = .FALSE.
C
C                  use results from TWOPNT problem at the top node
C                  as the starting guess for the bottom node
                   DO 235 K = 1 , KKGAS
                      RPAR(NS + NYS + IPTGAS + K -1) = RPAR(NXTW+K-1)
  235              CONTINUE
                   DO 240 K = 1, KKSURF
                      RPAR(NZ1 + K - 1) = RPAR(NXTW + KKGAS + K - 1)
  240              CONTINUE
                   CALL CRTWO (KKSURF+KKGAS, LOUT, IPMAX, LNTWOP,
     1              RPAR(NTWP),
     1              RPAR(NABV), RPAR(NBLW), RPAR(NBUF),
     2              RPAR(NXTW), RPAR(NA), IPAR(ITPW),
     3              IPAR, RPAR, IPAR(ISK), RPAR(NSK), PRES,
     4              RPAR(NZ1), RPAR(NYP), KKGAS, KKSURF, KKBULK,
     5              IPAR(ISIT), IPAR(IBLK), IPAR(IPKK), IPAR(IPKF),
     6              IPAR(IPKL), RPAR(NWDT), RPAR(NSDT), RPAR(NFZ),
     7              RPAR(NZN), RPAR(NFZN), IPAR(IIPV), RPAR(NWT),
     8              RPAR(NDDE), RPAR(NSDE), RPAR(NACT), IPAR(IFS),
     9              IPAR(ILS), IPAR(IFB), IPAR(ILB), CWORK(JKCH),
     *              CWORK(JPCH), PATM, IISUR, RPAR(NROP), CWORK(JSCH),
     1              RPAR(NYAV), RPAR(NS+IPTGAS), JJ, NATJ, LTDIF,
     2              IPAR(IIML), RPAR(NX), RPAR(NMAS), IPAR(IICR),
     3              IPAR(IMC), RPAR(NMC), RPAR(NXAV), RPAR(NHH),
     4              RPAR(NVIS), RPAR(NCND), RPAR(ND), RPAR(NTDC),
     5              RPAR(NYV), RPAR(NAA), RPAR(NRTP), RPAR(NDR),
     6              RPAR(NDC), RPAR(NXTP), SSABS, SSREL, TDABS, TDREL,
     7              IPRNT, LIMIT, NINIT, IRETIR, NJAC, ITJAC, DTMAX, 
     8              DTMIN, GFAC, SFAC, TSTEP0, LUPPER, LVCOR, KERR)
C
                  DO 245 K = 1, KKSURF
                      RPAR(NS + K - 1) = RPAR(NXTW + KKGAS + K - 1)
  245             CONTINUE
C
                  DO 250 K = 1 , KKGAS
                     RPAR(NS + NYS + IPTGAS + K -1) = RPAR(NXTW + K - 1)
  250             CONTINUE
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C
      CALL CRRUN (IPAR(IKK), IPAR(IKKT), JJ, NATJ, LOUT, ICRD, IRST,
     1            ITDIF, IMULT, TIME, RPAR(NS), RPAR(NSP), RPAR(NF),
     2            DX, XEND, RPAR(NWT), RPAR(NY), RPAR(NYP), RPAR(NXAV),
     3            RPAR(NYAV), RPAR(NYV), RPAR(NX), RPAR(NRR), RPAR(ND),
     4            RPAR(NTDC), IPAR, RPAR, IPAR(IMC), RPAR(NMC),
     5            RPAR(NDEP), IPAR(IPKK), IPAR(IPKF), IPAR(IPKL),
     6            IPAR(IIDS), RPAR(NDAS), LSAVE, LDAS, ATOL, RTOL,
     7            GTMP, STMP, RPAR(NXM), RPAR(NXP), RPAR(NMAS),
     8            RPAR(NCON), IFIXED, IPAR(ISK), RPAR(NSK), CWORK(JKCH),
     9            CWORK(JPCH), RPAR(NWDT), RPAR(NSDT), NNSURF, KKBULK,
     *            RPAR(NDDE), RPAR(NSDE), RPAR(NSPS), RPAR(NSPE),
     1            IPAR(IKKS), RPAR(NWDI), RPAR(NSDI), RPAR(NROP), IISUR,
     2            CWORK(JSCH), MORD, IPAR(IFS), IPAR(ILS), IPAR(IBLK),
     3            IPAR(IFB), IPAR(ILB), RPAR(NACT),
     4            RPAR(NASP), RPAR(NRSP), RPAR(NATS), RPAR(NRTS),
     5            RPAR(NGFC), RPAR(NSFC), LMOLF, HO, NPRNDS, IVCOR)
C
      CALL CRTIM (KKGAS, JJ, NATJ, LSAVE, LDAS, LOUT, IPAR(IIDS),
     1            RPAR(NDAS), IPAR(IKKS), ICRD)
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRDIFV (KKGAS, JJ, NATJ, LTDIF, IMULT, X, S, WT, TMASS,
     1                   ICRD, YAV, XMF, XMFP, D, TDC, IPAR, RPAR, YV,
     2                   LVCOR)
C
C  START PROLOGUE
C
C   KKGAS  - NUMBER OF CHEMICAL SPECIES.
C   JJ     - NUMBER OF MESH POINTS.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES AT EACH MESH POINT, AND
C            THE EXACT FIRST DIMENSION OF S(NATJ,*), SN(NATJ,*), AND
C            F(NATJ,*).  NATJ=KKGAS+2.
C   LTDIF  - IF LTDIF=.TRUE. THEN EVALUATE THERMAL DIFFUSION RATIOS AS
C            WELL AS DIFFUSION COEFFICIENTS.
C              CGS UNITS - GM/(CM**2-SEC)
C   LVCOR  - IF LVCOR=.TRUE. THEN ADD A CORRECTION VELOCITY TO SPECIES
C            DIFFUSION VELOCITIES.
C   IMULT  - IF IMULT = 1, THEN USE MULTICOMPONENT TRANSPORT COEFFICIENT
C            OTHERWISE, USE SIMPLER MIXTURE COEFFICIENTS.
C   X      - THE ARRAY OF MESH POINT LOCATIONS.
C              DIMENSION X(*) AT LEAST JJ
C              CGS UNITS - CM
C   S      - DEPENDENT VARIABLE MATRIX.  THE TEMPERATURES ARE STORED IN
C            T(J)=S(NT,J), THE MASS FRACTIONS ARE IN Y(K,J)=S(NYS+K,J),
C            AND THE FLOW RATES ARE IN FLRT(J)=S(NM,J).
C              DIMENSION S(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION,
C              AND AT LEAST JJ FOR THE SECOND.
C   WT     - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
C              CGS UNITS - GM/MOLE
C              DIMENSION WT(*) AT LEAST KKGAS.
C   TMASS  - THE INITIAL TOTAL MASS FLOW RATE (GM/SEC).
C   ICRD   - THE COORDINATE SYSTEM INDEX, ICRD=0 PLANAR,
C                                         ICRD=1 RADIAL.
C WORK AND SCRATCH SPACE
C   YAV    - ARRAY OF MASS FRACTIONS AT MESH MIDPOINTS.  YAV(K) IS THE
C            MASS FRACTION BETWEEN J AND J+1.
C              DIMENSION YAV(*) AT LEAST KKGAS.
C   XMF    - ARRAY OF MOLE FRACTIONS AT MESH POINT J.
C              DIMENSION XMF(*) AT LEAST KKGAS.
C   XMFP   - ARRAY OF MOLE FRACTIONS AT MESH POINT J+1.
C              DIMENSION XMFP(*) AT LEAST KKGAS.
C   D      - MATRIX OF SPECIES DIFFUSION COEFFICIENTS AT THE MESH
C            MIDPOINTS.  IF LVARTP=.TRUE. THESE ARE COMPUTED EACH
C            TIME THE FUNCTION IS CALLED, OTHERWISE THE STORED VALUES
C            ARE USED.
C              CGS UNITS - CM**2/SEC
C              DIMENSION D(KKGAS,KKGAS,*) EXACTLY KKGAS FOR THE FIRST
C              TWO DIMENSIONS AND AT LEAST JJ FOR THE THIRD.
C   TDC    - MATRIX OF SPECIES THERMAL DIFFUSION RATIOS AT THE MESH
C            MIDPOINTS.  IF LVARTP=.TRUE. THESE ARE COMPUTED EACH
C            TIME THE FUNCTION IS CALLED, OTHERWISE THE STORED VALUES
C            ARE USED.
C              CGS UNITS - NONE
C              DIMENSION TDC(KKGAS,*) EXACTLY KKGAS FOR THE FIRST
C              DIMENSION AND AT LEAST JJ FOR THE SECOND.
C   IPAR    - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   RPAR    - FLOATING POINT CHEMKIN WORK SPACE.
C                                         ICRD=2 PLANAR W/ SYMMETRY AXIS
C
C
C
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C OUTPUT-
C   YV     - MATRIX OF MASS FRACTIONS TIMES DIFFUSION VELOCITIES AT THE
C            MESH MIDPOINTS.  YV(K,J) IS THE FLUX OF KTH SPECIES BETWEEN
C            J AND J+1.
C              DIMENSION YV(KKGAS,*) EXACTLY KKGAS FOR THE FIRST
C              DIMENSION AND AT LEAST JJ FOR THE SECOND.
C
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
C
      DIMENSION S(NATJ,*), X(*), WT(*), YAV(*), XMF(*), XMFP(*),
     1          D(KKGAS,KKGAS,*), TDC(KKGAS,*), YV(KKGAS,*), IPAR(*),
     2          RPAR(*)
C
      LOGICAL LTDIF, LVCOR
C
      DATA EPS/1.0E-30/
C
      CALL CKYTX (S(NY1,1), IPAR, RPAR, XMFP)
C
C     Loop over all mesh points, computing the diffusion velocity 
C     at the midpoints; the indexing is such that YV(K,J) is the
C     diffusion velocity of species K midway between nodes J and J+1
C
      DO 1000 J = 1, JJ-1
         JP1 = J + 1
C
         TMF = TMASS + S(NML,J) + S(NMU,J)
         UAV = 0.5 * (S(NU,J) + S(NU,JP1))
         PAV = 0.5 * (S(NP,J) + S(NP,JP1))
         TAV = 0.5 * (S(NT,J) + S(NT,JP1))
         IF (ICRD .EQ. 1) THEN
C           radial coordinates
            RTICRD = 0.5 * (SQRT(ABS(S(NR,J))) + SQRT(ABS(S(NR,JP1))))
         ELSE
C           planar coordinates
            RTICRD = 1.
         ENDIF
C
         CALL CKAVG (KKGAS, S(NYS+1,J), S(NYS+1,JP1), YAV)
         CALL CKRHOY (PAV, TAV, YAV, IPAR, RPAR, RHOAV)
         CALL CKMMWY (YAV, IPAR, RPAR, WTMAV)
C
         DO 400 K = 1, KKGAS
            XMF(K) = XMFP(K)
400      CONTINUE
         CALL CKYTX (S(NY1,JP1), IPAR, RPAR, XMFP)
C
         TMFAC = ((X(JP1) - X(J))*TMF)
         DO 500 K = 1, KKGAS
C
            IF (IMULT .EQ. 1) THEN
C              multicomponent diffusion coefficients
C
               SUMN = 0.0
               SUMD = 0.0
               DO 450 L=1,KKGAS
                  SUMN = SUMN + WT(L) * D(K,L,J) *
     1                    ((XMFP(L) - XMF(L))) /
     2                    TMFAC
C     2                    ((X(JP1) - X(J))*TMF)
                  SUMD = SUMD + ((XMFP(L) - XMF(L))) /
     1                     TMFAC
C     1                    ((X(JP1) - X(J))*TMF)
450            CONTINUE
               SUMD = SUMD - ((XMFP(K) - XMF(K))) /
     1                  TMFAC
C     1                 ((X(JP1) - X(J))*TMF)
               DKM = (SUMN + EPS) / (WTMAV * (SUMD + EPS))
C
            ELSE
C              "mixture" diffusion coefficients
C
               DKM = D(K,1,J)
            ENDIF
C
            YV(K,J) = - WT(K)/WTMAV * RHOAV * UAV * RTICRD *
     1                  DKM * (XMFP(K) - XMF(K)) /
     2                  TMFAC
C     2                     ((X(JP1) - X(J))*TMF)
500      CONTINUE
C
         IF (LTDIF) THEN
C           add the thermal diffusion, if requested
C
            TMFAC = (S(NT,JP1)-S(NT,J)) / ((X(JP1)-X(J))*TMF)
            DO 600 K = 1, KKGAS
               YV(K,J) = YV(K,J) -
     1                   UAV * RTICRD * TDC(K,J) / TAV *
     2                   TMFAC 
C     2                   (S(NT,JP1)-S(NT,J)) / ((X(JP1)-X(J))*TMF)
600         CONTINUE
         ENDIF
C
         IF (LVCOR) THEN
C           add a correction velocity, if requested
C
            SUM = 0.0
            DO 700 K = 1, KKGAS
               SUM = SUM + YV(K,J)
700         CONTINUE
C
            DO 800 K = 1, KKGAS
               YV(K,J) = YV(K,J) - YAV(K) * SUM
800         CONTINUE
         ENDIF
C
1000  CONTINUE
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRFCN (DIST, S, SP, F, IRES, RPAR, IPAR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION S(*), SP(*), F(*), RPAR(*), IPAR(*)
C
      COMMON /CRCRCR/ NCK,  NMC,  NSK,  NWT,  NMAS, NSLB, NSRB, NX,
     1                NZ1,  NZJJ, NS,   NSP,  NYAV, NYV,  NHH,  NCON,
     2                NVIS, NCND, NRR,  NY,   NYP,  NROP, NXAV, NXM,
     3                NXP,  ND,   NTDC, NF,   NDEP, NEKR, NS1,  NS2,
     4                NQ1,  NQ2,  NQ12, NWDT, NSDT, NDAS, NCP,  NSPS,
     5                NSPE, NSTS, NDDE, NSDE, NXTW, NFZ,  NZN,  NFZN,
     6                NWDI, NSDI, NGTP, NPR,  NTWP, NABV, NBLW, NBUF,
     7                NA,   NACT, NAA,  NRTP, NDR,  NDC,  NXTP, NASP,
     8                NRSP, NATS, NRTS, NGFC, NSFC, NGASW,ICK,  IMC,
     9                ISK,  IJJ,  INAJ, IKK,  IKKT, IICR, IITD, IIML,
     *                ILRO, ILRI, ILDS, ILOU, IIFX, IPKK, IPKF, IPKL,
     1                IIDS, ISIT, IFS,  ILS,  IBLK, IFB,  ILB,  IKKS,
     2                ITPW, IIPV, IMO,  IIVC, JCCH, JSCH, JKCH, JPCH
C
      COMMON /LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
      LOGICAL LVARTP, LTDIF, LVCOR
C
      ITSAFE = 60
      ITLEFT = ITSAFE
C
C*****unicos timelimit
C      CALL TREMAIN (TLEFT)
C      ITLEFT = INT (TLEFT)
C*****END unicos timelimit
C
      IF (ITLEFT .LT. ITSAFE) THEN
         CALL CRTIM (IPAR(IKK), IPAR(IJJ), IPAR(INAJ), IPAR(ILRO),
     1               IPAR(ILDS), IPAR(ILOU), IPAR(IIDS), RPAR(NDAS),
     2               IPAR(IKKS), IPAR(IICR))
         RETURN
      ENDIF
C
C     IPTGAS points to the location in the solution vector of the
C     first gas-phase node
C
      IF (IPAR(IICR) .EQ. 0) THEN
C        planar coordinates
         IPTGAS = IPAR(IKKS) + 1
      ELSE
C        radial coordinates (or planar with a symmetry axis)
         IPTGAS = 1
      ENDIF
      IPTTOP = IPTGAS + IPAR(IJJ)*IPAR(INAJ)
C
      DO 100 J = 1, IPAR(IJJ)
         IF (S(IPTGAS-1 + IPAR(INAJ)*(J-1) + NT).LE.200.)THEN
            IRES = -1
            RETURN
         ENDIF
100   CONTINUE
C
      LVARTP = (IRES .NE. 1)
      LTDIF  = (IPAR(IITD) .EQ. 1)
      LVCOR  = (IPAR(IIVC) .EQ. 1)
C
      CALL CRFUN (IPAR(IKK), IPAR(IJJ), IPAR(INAJ), LTDIF, IPAR(IIFX),
     1            IPAR(IIML), LVARTP, DIST, RPAR(NWT), RPAR(NMAS),
     2            IPAR(IICR), RPAR(NSLB), RPAR(NSRB), RPAR(NX),
     3            S(IPTGAS), S(1), S(IPTTOP), SP(IPTGAS), SP(1),
     4            SP(IPTTOP), RPAR(NYAV), RPAR(NYV), RPAR(NWDT),
     5            RPAR(NSDT), RPAR(NCP),  RPAR(NHH), RPAR(NCON),
     6            RPAR(NVIS), RPAR(NCND), RPAR(ND),  RPAR(NTDC),
     7            IPAR, RPAR, IPAR(IMC),  RPAR(NMC), F(IPTGAS),
     8            F(1), F(IPTTOP),  IPAR(ISK),  RPAR(NSK), RPAR(NSPS),
     9            RPAR(NSPE), RPAR(NSTS), IPAR(ISIT),RPAR(NSDE),
     *            IPAR(IKKS), RPAR(NACT), IPAR(IFS), IPAR(ILS),
     1            IPAR(IPKK), IPAR(IPKF), IPAR(IPKL), RPAR(NGFC),
     2            RPAR(NSFC), LVCOR)
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRFLUX (KKGAS, KKSURF, KKTOT, JJ, NATJ, ICRD, S, WT,
     1                   IPAR, RPAR, CON, DEP, ILOW, Z, ISKWRK,
     2                   RSKWRK, SDOT, SITDOT, SDEN, ACT, NFSURF,
     3                   NLSURF, KFIRST, KLAST, SFAC)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION S(NATJ, *), WT(*), IPAR(*), RPAR(*), DEP(*), CON(*),
     1          Z(*), ISKWRK(*), RSKWRK(*), SDOT(*), SITDOT(*), SDEN(*),
     2          ACT(*), KFIRST(*), KLAST(*)
C
      COMMON /LOCS/    NP, NR, NU, NT, NML, NMU, NYS, NY1
C
C     determine the boundary indices, depending on the 
C     coordinate system
      IF (ILOW .EQ. 0) THEN
C       upper boundary
        N1 = JJ
        N2 = JJ-1
      ELSE
C       lower boundary
        N1 = 1
        N2 = 2
      ENDIF
C
      P = S(NP,1)
      T = S(NT,N1)
      IF (ICRD .EQ. 1) THEN
C        radial coordinates
         R = 0.5 * (SQRT(ABS(S(NR,N1))) + SQRT(ABS(S(NR,N2))))
      ELSE
C        planar coordinates
         R = 1.
      ENDIF
C
      CALL CKYTX (S(NY1,N1), IPAR, RPAR, ACT)
      IF (NFSURF .GT. 0) THEN
        DO 105 K = KFIRST(NFSURF), KLAST(NLSURF)
           ACT(K) = Z(K - KKGAS)
  105   CONTINUE
      ENDIF
C
      CALL SKSDEN (ISKWRK, RSKWRK, SDEN)
      CALL SKRAT  (P, T, ACT, SDEN, ISKWRK, RSKWRK, SDOT, SITDOT)
C
      DO 225 K = 1, KKTOT
         SDOT(K) = SFAC*SDOT(K)
  225 CONTINUE
C
      DO 250 K = 1, KKGAS
         DEP(K) = SDOT(K) * WT(K) * R
  250 CONTINUE
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRFUN (KKGAS, JJ, NATJ, LTDIF, IFIXED, IMULT, LVARTP,
     1                  DIST, WT, TMASS, ICRD, SLB, SRB, X, S, Z1, ZJJ,
     2                  SP, Z1P, ZJJP, YAV, YV, SDOT, SITDOT, CP, H,
     3                  CON, VISC, COND, D, TDC, IPAR, RPAR, IMCWRK,
     4                  RMCWRK, F, FZ1, FZJJ, ISKWRK, RSKWRK, XSPSTR,
     5                  XSPEND, TGTMP, NNSURF, SDEN, KKSURF, ACT,
     6                  NFSURF, NLSURF, KKPHAS, KFIRST, KLAST, GFAC,
     7                  SFAC, LVCOR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/    NP, NR, NU, NT, NML, NMU, NYS, NY1
      COMMON /CONSTS/  RUC, PATM, GRAV, RHOE
C
      DIMENSION S(NATJ,*), Z1(*), ZJJ(*), SP(NATJ,*), Z1P(*), ZJJP(*),
     1          F(NATJ,*), FZ1(*), FZJJ(*), WT(*), X(*), SLB(*), SRB(*),
     2          YAV(*), YV(KKGAS,*), SDOT(*), SITDOT(*), CP(*), H(*),
     3          CON(*), VISC(*), COND(*), D(KKGAS,KKGAS,*),
     4          TDC(KKGAS,*), IPAR(*), RPAR(*), IMCWRK(*), RMCWRK(*),
     5          ISKWRK(*), RSKWRK(*), SDEN(*), ACT(*), KKPHAS(*),
     6          KFIRST(*), KLAST(*)
C
      LOGICAL LTDIF, LVARTP, LVCOR
      EXTERNAL VTRACE

C*****precision > double
      DOUBLE PRECISION ZERO, ONE, HALF, TWO
C*****END precision > double
C*****precision > single
C      REAL ZERO, ONE, HALF, TWO
C*****END precision > single
      SAVE ZERO, ONE, HALF, TWO
C
C*****precision > double
      DATA ZERO/0.0D0/ , ONE/1.0D0/, HALF/0.5D0/, TWO/2.0D0/
C*****END precision > double
C*****precision > single
C      DATA ZERO/0.0/ , ONE/1.0/, HALF/0.5/, TWO/2.0/
C*****END precision > single
C
C     Evaluate and store the transport coefficients
C     (in the following call h(*) and cp(*) are used temporarily
C     for scratch space)
      IF (LVARTP) CALL CRTRNP (KKGAS, JJ, NATJ, LTDIF, IMULT, X, S, YAV,
     1                         CP, IPAR, RPAR, IMCWRK, RMCWRK, VISC,
     2                         COND, D, TDC)
C
C     Evaluate and store the diffusion velocities
      CALL CRDIFV (KKGAS, JJ, NATJ, LTDIF, IMULT, X, S, WT, TMASS, ICRD,
     1             YAV, H, CP, D, TDC, IPAR, RPAR, YV, LVCOR)
C
C-----LOWER BOUNDARY
C
      PAV = HALF * (S(NP,1) + S(NP,2))
      TAV = HALF * (S(NT,1) + S(NT,2))
      UAV = HALF * (S(NU,1) + S(NU,2))
      IF (ICRD .EQ. 1) THEN
C        radial coordinates
         RJ     = ZERO
         R2P    = HALF * (S(NR,1) + S(NR,2))
         RAVP   = HALF * (SQRT(ABS(S(NR,1))) + SQRT(ABS(S(NR,2))))
      ELSE
C        planar coordinates
         RJ     = ONE
         R2P    = ONE
         RAVP   = ONE
      ENDIF
      CALL CKAVG (KKGAS, S(NYS+1,1), S(NYS+1,2), YAV)
      CALL CKRHOY (S(NP,1), S(NT,1), S(NY1,1), IPAR, RPAR, RHO)
      RHOU = RHO * S(NU,1)
      CALL CKRHOY (PAV, TAV, YAV, IPAR, RPAR, RHOAP)
      RHOUAP = RHOAP * UAV
C
C
      IF (ICRD .EQ. 0) THEN
C        planar coordinates (lower boundary is susceptor)
C
         F(NP,1) = S(NP,2) - S(NP,1)
         F(NR,1) = (SLB(NR) - S(NR,1))
         F(NU,1) = SLB(NU) - S(NU,1)
         CALL CRWALL (XSPSTR, XSPEND, DIST, TEMP)
         F(NT,1) = TEMP - S(NT,1)
C
C        total mass loss at lower boundary
         CALL CKYTX (S(NY1,1), IPAR, RPAR, ACT)
         IF (NFSURF .GT. 0) THEN
            DO 115 K = KFIRST(NFSURF), KLAST(NLSURF)
               ACT(K) = Z1(K - KKGAS)
  115       CONTINUE
         ENDIF
C
         CALL SKSDEN (ISKWRK, RSKWRK, SDEN)
         CALL SKRAT  (S(NP,1), S(NT,1), ACT, SDEN, ISKWRK, RSKWRK,
     1                   SDOT, SITDOT)
C
         DO 525 K = 1, KKGAS+KKSURF
            SDOT(K) = SFAC*SDOT(K)
  525    CONTINUE
C
C        calculate convective velocity
         VCON = ZERO
         DO 550 K = 1, KKGAS
            VCON = VCON + SDOT(K) * WT(K) / RHOAP
  550    CONTINUE
         DMDX = RHOAP * VCON
C
         F(NML, 1) = SP(NML, 1) - DMDX
         F(NMU, 1) = S (NMU, 2) - S (NMU, 1)
C
C        boundary conditions for chemical species
         DO 700 K = 1, KKGAS
            N = NYS + K
            F(N,1) = RHOAP * (S(N,1)*VCON + YV(K,1)) -
     1                   SDOT(K)*WT(K)
  700    CONTINUE
C
         IF (.NOT. LVCOR) THEN
C           force the sum of the mass fractions to add to one
            F(NYS+KKGAS,1) = VTRACE (S(NYS+1,1), KKGAS)
         ENDIF
C
C        conservation equations for surface species
         IF (KKSURF .GT. 0) THEN
            DO 775 N = NFSURF, NLSURF
               SUMZ1 = ZERO
               IF (KKPHAS(N) .GT. 1) THEN
                  DO 770 K = KFIRST(N), KLAST(N) - 1
                     FZ1(K - KKGAS) = SDOT(K) / SDEN(N)
                     SUMZ1 = SUMZ1 + Z1(K - KKGAS)
770             CONTINUE
               ENDIF
               FZ1(KLAST(N)-KKGAS) = SUMZ1 + Z1(KLAST(N)-KKGAS) - ONE
  775       CONTINUE
         ENDIF
C
      ELSE
C
C        radial coordinates (or planar with a symmetry axis),
C        lower boundary is centerline
C
         DO 750 N = 1, NATJ
            F(N,1) = S(N,2) - S(N,1)
750      CONTINUE
         F(NR,1)  = (SLB(NR)-S(NR,1)) * 1.0E18
         F(NML,1) = SLB(NML) - S(NML,1)
C
      ENDIF
C
C-----INTERIOR MESH POINTS
C
      DO 2000 J = 2, JJ-1
         JP1 = J + 1
         JM1 = J - 1
C
         TMF = TMASS + S(NML, JP1) + S(NMU, JP1)
         RHOUM = RHOU
         CALL CKRHOY (S(NP,J), S(NT,J), S(NY1,J), IPAR, RPAR, RHO)
         RHOU = RHO * S(NU,J)
C
C        average solution at J + 1/2
         RHOAM = RHOAP
         RHOUAM = RHOUAP
         PAV = HALF * (S(NP,J) + S(NP,JP1))
         TAV = HALF * (S(NT,J) + S(NT,JP1))
         UAV = HALF * (S(NU,J) + S(NU,JP1))
         CALL CKAVG (KKGAS, S(NYS+1,J), S(NYS+1,JP1), YAV)
         CALL CKRHOY (PAV, TAV, YAV, IPAR, RPAR, RHOAP)
         RHOUAP = RHOAP * UAV
C
         RAVM = RAVP
         R2M = R2P
         IF (ICRD .EQ. 1) THEN
C           radial coordinates
            RJ = SQRT(ABS(S(NR, J)))
            R2P = HALF*(S(NR,J)+S(NR,JP1))
            RAVP = HALF * (SQRT(ABS(S(NR,J))) + SQRT(ABS(S(NR,JP1))))
         ELSE
C           planar coordinates
            RJ = ONE
            R2P = ONE
            RAVP = ONE
         ENDIF
C
C        coordinate differences
         DXP = (X(JP1)-X(J))*TMF
         DXM = (X(J)-X(JM1))*TMF
         DXAV = 0.5*(X(JP1)-X(JM1))*TMF
         DXPM = (X(JP1)-X(JM1))*TMF
C
C        thermodynamic properties at J
         CALL CKCPBS (S(NT,J), S(NY1,J), IPAR, RPAR, CPB)
         CALL CKCPMS (S(NT,J), IPAR, RPAR, CP)
         CALL CKHML  (S(NT,J), IPAR, RPAR, H)
C
C        chemistry rates at J
         CALL CKWYP (S(NP,J), S(NT,J), S(NY1,J), IPAR, RPAR, SDOT)
         DO 1040 K = 1, KKGAS
            SDOT(K) = SDOT(K) * GFAC
 1040    CONTINUE
C
         TDOT = ZERO
         DO 1050 K = 1, KKGAS
            TDOT = TDOT + H(K)*SDOT(K)
 1050    CONTINUE
C
C        momentum equation
         T1 = SP(NU,J)
         T2 = SP(NP,J) / (RHO * S(NU,J))
         DUDXP = (S(NU,JP1) - S(NU,J)) / DXP
         DUDXM = (S(NU,J) - S(NU,JM1)) / DXM
         DUDXA = (S(NU,JP1) - S(NU,JM1)) / DXPM
         T3 = (RHOUAP*R2P*VISC(J)*DUDXP -
     1                 RHOUAM*R2M*VISC(JM1)*DUDXM) / DXAV
         T4 = ( -SP(NML,J) + X(J)*(SP(NML,J)+SP(NMU,J)) ) * DUDXA
         T5 = GRAV * (RHOE - RHO) / (RHO * S(NU,J))
         F(NU,J) = (T1 - T4 + T2 - T3 - T5)
C
C        energy equation
         DTDXP = (S(NT,JP1) - S(NT,J))  / DXP
         DTDXM = (S(NT,J) - S(NT,JM1))  / DXM
         DTDXA = (S(NT,JP1) - S(NT,JM1)) / DXPM
         T1 = SP(NT,J)
         T2 = (RHOUAP*R2P*COND(J)*DTDXP -
     1                 RHOUAM*R2M*COND(JM1)*DTDXM) / (DXAV*CPB)
         T3 = TDOT / (RHO * S(NU,J) * CPB)
         SUM = ZERO
         DO 1350 K = 1, KKGAS
            YVAV = HALF * (YV(K,J) + YV(K,JM1))
            SUM = SUM + CP(K)*YVAV
 1350    CONTINUE
         T4 = RHO / CPB * RJ * DTDXA * SUM
         T5 = ( -SP(NML,J) + X(J)*(SP(NML,J)+SP(NMU,J)) ) * DTDXA
         F(NT,J) = (T1 - T5 - T2 + T3 + T4)
C
C        mass loss equations
         F(NML, J) = S(NML, J) - S(NML, JM1)
         F(NMU, J) = S(NMU, J) - S(NMU, JP1)
C
C        pressure equation
         F(NP,J) = (S(NP,JP1) - S(NP,J))
C
C        Y coordinate equation
         IF (ICRD .EQ. 1) THEN
C           radial coordinates
            F(NR,J) = -(S(NR,J)-S(NR,JM1))/DXM + (4.0/(RHOU + RHOUM))
C
         ELSE
C           planar coordinates
            F(NR,J) = - (S(NR,J)-S(NR,JM1))/DXM + TWO/((RHOU + RHOUM))
C
         ENDIF
C
C        species equations
         T5 = ( -SP(NML,J) + X(J)*(SP(NML,J)+SP(NMU,J)) ) 
         DO 1400 K = 1, KKGAS
            N = NYS + K
            T1 = SP(N,J)
            T2 = (RHOAP*RAVP*YV(K,J) - RHOAM*RAVM*YV(K,JM1)) / DXAV
            T3 = SDOT(K)*WT(K) / (RHO * S(NU,J))
            DSDXA = (S(N,JP1) - S(N,JM1)) / DXPM
C            T5 = ( -SP(NML,J) + X(J)*(SP(NML,J)+SP(NMU,J)) ) * DSDXA
C            F(N,J) = (T1 - T5 + T2 - T3)
            F(N,J) = (T1 - T5*DSDXA + T2 - T3)
 1400   CONTINUE
C
        IF (.NOT. LVCOR) THEN
C          force the sum of the mass fractions to add to one
           F(NYS+KKGAS, J) = VTRACE (S(NYS+1,J), KKGAS)
        ENDIF
C
 2000 CONTINUE
C
C-----UPPER BOUNDARY NODE
C
      DXM = DXP
      RHOUM = RHOU
      CALL CKRHOY (S(NP,JJ), S(NT,JJ), S(NY1,JJ), IPAR, RPAR, RHO)
      RHOU = RHO * S(NU,JJ)
C
      JJM1 = JJ - 1
      DO 2050 N = 1, NATJ
         F(N, JJ) = S(N, JJM1) - S(N, JJ)
 2050 CONTINUE
C
      F(NP, JJ) = SRB(NR) - S(NR, JJ)
      F(NU, JJ) = SRB(NU) - S(NU, JJ)
C
      IF (ICRD .EQ. 1) THEN
C         radial coordinates
          CALL CRWALL (XSPSTR, XSPEND, DIST, TEMP)
          F(NT,JJ) = TEMP - S(NT,JJ)
          F(NR,JJ) = - ( S(NR,JJ) - S(NR,JJM1) ) / DXM +
     1              (4.0 / (RHOU + RHOUM) )
C
      ELSE IF (ICRD .EQ. 0) THEN
C         planar coordinates
          IF (IFIXED .EQ. 2) THEN
            CALL CRWALL (XSPSTR, XSPEND, DIST, TEMP)
            F(NT,JJ) = TEMP - S(NT,JJ)
          ELSE IF (IFIXED .EQ. 1) THEN
             F(NT, JJ) = S(NT, JJ) - TGTMP
          ELSE
             F(NT, JJ) = S(NT, JJ) - S(NT, JJM1)
          ENDIF
          F(NR,JJ) = - ( S(NR, JJ) - S(NR, JJM1) ) / DXM +
     1                  (TWO / (RHOU + RHOUM) )
C
      ELSE
C        planar coordinates with a symmetry axis
         CALL CRWALL (XSPSTR, XSPEND, DIST, TEMP)
         F(NT,JJ) = TEMP - S(NT,JJ)
         F(NR,JJ) = - ( S(NR, JJ) - S(NR, JJM1) ) / DXM +
     1                  (TWO / (RHOU + RHOUM) )
      ENDIF
C
      CALL CKYTX (S(NY1,JJ), IPAR, RPAR, ACT)
      IF (KKSURF .GT. 0) THEN
        DO 125 K = KFIRST(NFSURF), KLAST(NLSURF)
           ACT(K) = ZJJ(K - KKGAS)
  125   CONTINUE
      ENDIF
C
      CALL SKSDEN (ISKWRK, RSKWRK, SDEN)
      CALL SKRAT (S(NP,JJ), S(NT,JJ), ACT, SDEN, ISKWRK, RSKWRK,
     1             SDOT, SITDOT)
C
      DO 2075 K = 1, KKGAS+KKSURF
         SDOT(K) = SFAC*SDOT(K)
 2075 CONTINUE
C
C     calculate convective velocity
      VCON = ZERO
      DO 2100 K = 1, KKGAS
         VCON = VCON - SDOT(K) * WT(K) / RHOAP
 2100 CONTINUE
      DMDX = - RAVP * RHOAP * VCON
C
C     boundary conditions on chemical species
      DO 2300 K = 1, KKGAS
         N = NYS + K
         F(N,JJ) = -RHOAP * (S(N,JJ)*VCON + YV(K,JJM1)) -
     1                 SDOT(K) * WT(K)
 2300 CONTINUE
C
      F(NML, JJ) = S (NML, JJ) - S(NML, JJM1)
      F(NMU, JJ) = SP(NMU, JJ) - DMDX
C
      IF (.NOT. LVCOR) THEN
C        force the sum of the mass fractions to add to one.
         F(NYS+KKGAS, JJ) = VTRACE (S(NYS+1,JJ), KKGAS)
      ENDIF
C
      IF (KKSURF .LE. 0) RETURN
C
C     conservation equations for surface species
      DO 2775 N = NFSURF, NLSURF
         SUMZJJ = 0.0
         IF (KKPHAS(N) .GT. 1) THEN
            DO 2770 K = KFIRST(N), KLAST(N) - 1
               FZJJ(K - KKGAS) = SDOT(K) / SDEN(N)
               SUMZJJ = SUMZJJ + ZJJ(K - KKGAS)
 2770       CONTINUE
         ENDIF
         FZJJ(KLAST(N) - KKGAS ) = SUMZJJ + ZJJ(KLAST(N) - KKGAS) - ONE
 2775 CONTINUE
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRINIT (KKGAS, WT, IPAR, RPAR, PI, NATJ, ICRD, JJ,
     1                   GTMP, P, VEL, REAC, HITE, STCH, S, R, X, SLB,
     2                   SRB, TMF, Y, XSPSTR, XSPEND, LTHIC, BLTK,
     3                   IFIXED, GASW)
C
C  START PROLOGUE
C     GIVEN THE INPUT DATA, THIS SUBROUTINE COMPUTES THE MESH AND
C     INSERTS THE INITIAL CONDITIONS INTO THE DEPENDENT VARIABLE
C     ARRAY, S.
C
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C*****precision > quad
C      REAL TSING
C*****END precision > quad
C
      DIMENSION S(NATJ,*), IPAR(*), RPAR(*), WT(*), SLB(*),
     1          SRB(*), R(*), X(*), REAC(*), Y(*), GASW(*)
      LOGICAL LTHIC
C
      COMMON /LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
      COMMON /CONSTS/ RUC, PATM, GRAV, RHOE
C
C     set the cross stream physical coordinate (R)
C     HITE = channel dimension
C     STCH = parameter for skewing the grid
C
C*****precision > quad
C      FAC = HITE/(QEXT(JJ-1))**STCH
C*****END precision > quad
C*****precision > double
      FAC = HITE/(DBLE(JJ-1))**STCH
C*****END precision > double
C*****precision > single
C      FAC = HITE/(REAL(JJ-1))**STCH
C*****END precision > single
C
      IF (ICRD .EQ. 0) THEN
C        planar coordinates
         DO 100 J = 1, JJ
C*****precision > quad
C            R(J) = FAC*(QEXT(J-1))**STCH
C*****END precision > quad
C*****precision > double
            R(J) = FAC*(DBLE(J-1))**STCH
C*****END precision > double
C*****precision > single
C             R(J) = FAC*(REAL(J-1))**STCH
C*****END precision > single
  100    CONTINUE
C
      ELSE
C        radial coordinates (or planar with a symmetry axis)
         DO 120 J = 1, JJ
C*****precision > quad
C            R(J) = HITE - FAC * QEXT(JJ-J)**STCH
C*****END precision > quad
C*****precision > double
            R(J) = HITE - FAC * DBLE(JJ-J)**STCH
C*****END precision > double
C*****precision > single
C            R(J) = HITE - FAC * REAL(JJ-J)**STCH
C*****END precision > single
  120    CONTINUE
C
      ENDIF
C
C     put temperature into S
      DIST = 0.0
      CALL CRWALL (XSPSTR, XSPEND, DIST, TEMP)
      DELTA = TEMP - GTMP
      IF (LTHIC) THEN
         DO 200 J = 1, JJ
            IF ( R(J) .LT. BLTK) THEN
               IF (ICRD .EQ. 0) THEN
C                 planar coordinates, interpolate between the
C                 wall and gas temperatures
                  S(NT,J) = GTMP + (BLTK-R(J))/BLTK * DELTA
               ELSE
C                 radial coordinates (or planar with a symmetry axis)
                  S(NT,J) = GTMP
               ENDIF
C
            ELSE IF ( R(J).GE.BLTK .AND. R(J).LE.(HITE-BLTK)) THEN
               S(NT,J) = GTMP
C
            ELSE IF ( R(J).GT.(HITE-BLTK)) THEN
C              interpolate between wall and gas temperature
               S(NT,J) = GTMP + (R(J)-HITE+BLTK)/BLTK * DELTA
            ENDIF
  200    CONTINUE
      ELSE
         DO 201 J = 1, JJ
            S(NT,J) = GTMP
  201    CONTINUE
      ENDIF
C
      IF (ICRD .EQ. 0) THEN
C        planar coordinates
         S(NT,1)  = TEMP
         IF (IFIXED .EQ. 2) S(NT,JJ) = TEMP
      ELSE
C        radial coordinates (or planar with a symmetry axis)
         S(NT,JJ) = TEMP
      ENDIF
C
C     parabolic velocity profile
      IF (LTHIC) THEN
C
         BLTKSQ = BLTK**2 
         IF (ICRD .EQ. 0) THEN
C           planar coordinates
            DO 310 J = 1, JJ
               IF (R(J) .LT. BLTK) THEN
                  S(NU,J) = -VEL * (R(J)-BLTK)**2 / BLTKSQ + VEL
               ELSEIF (R(J).GE.BLTK .AND. R(J).LE.(HITE-BLTK)) THEN
                  S(NU,J) = VEL
               ELSEIF (R(J) .GT. (HITE-BLTK)) THEN
                  S(NU,J) = -VEL * (R(J)-HITE+BLTK)**2 / BLTKSQ + VEL
               ENDIF
  310       CONTINUE
         ELSE
C           radial coordinates, or planar with a symmetry axis
            DO 340 J = 1, JJ
               IF (R(J) .LT. (HITE-BLTK)) THEN
                  S(NU,J) = VEL
               ELSE
                  S(NU,J) = -VEL * (R(J)-HITE+BLTK)**2 / BLTKSQ + VEL
               ENDIF
  340       CONTINUE
            S(NU,1) = S(NU,2)
         ENDIF
C
      ELSE
         IF (ICRD .EQ. 0) THEN
C           planar coordinates
            A = -4.0 * VEL/HITE**2
            DO 300 J = 1, JJ
               S(NU,J) = A * (R(J)-0.5*HITE)**2 + VEL
  300       CONTINUE
C
         ELSEIF (ICRD .EQ. 1) THEN
C           radial coordinates
            A = -4.0 * VEL/HITE**2
            DO 320 J = 1, JJ
               S(NU,J) = 0.25 * A * R(J)**2 + VEL
  320       CONTINUE
            S(NU,1) = S(NU,2)
C
         ELSE
C           planar with a symmetry axis
            A = -VEL/HITE**2
            DO 360 J = 1, JJ
               S(NU,J) = A * R(J)**2 + VEL
  360       CONTINUE
            S(NU,1) = S(NU,2)
         ENDIF
      ENDIF    
C
C     put P and R into S vector
      IF (ICRD .EQ. 1) THEN
C        radial coordinates
         DO 380 J = 1, JJ
            S(NP,J) = P
            S(NR,J) = R(J)**2
  380    CONTINUE
      ELSE
         DO 400 J = 1, JJ
            S(NP,J) = P
            S(NR,J) = R(J)
  400    CONTINUE
      ENDIF
C
C     convert REAC to mass fraction, and put in S
      CALL CKXTY (REAC, IPAR, RPAR, Y)
      CALL CKXTY (GASW, IPAR, RPAR, S(NY1,JJ))
      DO 500 J = 1, JJ-1
         CALL CKXTY (REAC, IPAR, RPAR, S(NY1,J))
         S(NYS+KKGAS, J) = VTRACE (S(NYS+1, J), KKGAS-1)
  500 CONTINUE
C
C     convert physical dimension (R) to streamfunction (X)
      CALL CRRTOX (ICRD, NATJ, JJ, R, S, IPAR, RPAR, X)
C
C     compute total mass flux, and store in S(NM,J)
      TMF = X(JJ)
      IF (ICRD .EQ. 1) TMF = TMF * 2.0 * PI
C
      DO 550 J = 1, JJ
         S(NML, J) = 0.
         S(NMU, J) = 0.
C        normalize streamfunction
         X(J) = X(J) / TMF
  550 CONTINUE
C
C     set up boundary conditions
      DO 650 N = 1, NATJ
         SLB(N) = S(N,1)
         SRB(N) = S(N,JJ)
  650 CONTINUE
C
C     calculate reference density for buoyancy term
      IF (GRAV .NE. 0.0) CALL CKRHOY (P, TEMP, Y, IPAR, RPAR, RHOE)
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRPNT (LINKCK, LINKMC, LINKSK, LOUT, NATJ, LNIPAR,
     1                  LNRPAR, LNCWRK, KKGAS, KKSURF, KKBULK, KKTOT,
     2                  NPHASE, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     3                  NLBULK, II, IISUR, LNTWOP, IPCALL, JJ, ICRD,
     4                  LNLWRK, IPAR, RPAR, CWORK, KERR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION IPAR(*), RPAR(*)
      CHARACTER CWORK(*)*(*)
      LOGICAL IERR, KERR
      INTEGER COMPS
      INTEGER NOLD, NASPO, NRSPO
      SAVE    NOLD, NASPO, NRSPO
C
      COMMON /CRCRCR/ NCK,  NMC,  NSK,  NWT,  NMAS, NSLB, NSRB, NX,
     1                NZ1,  NZJJ, NS,   NSP,  NYAV, NYV,  NHH,  NCON,
     2                NVIS, NCND, NRR,  NY,   NYP,  NROP, NXAV, NXM,
     3                NXP,  ND,   NTDC, NF,   NDEP, NEKR, NS1,  NS2,
     4                NQ1,  NQ2,  NQ12, NWDT, NSDT, NDAS, NCP,  NSPS,
     5                NSPE, NSTS, NDDE, NSDE, NXTW, NFZ,  NZN,  NFZN,
     6                NWDI, NSDI, NGTP, NPR,  NTWP, NABV, NBLW, NBUF,
     7                NA,   NACT, NAA,  NRTP, NDR,  NDC,  NXTP, NASP,
     8                NRSP, NATS, NRTS, NGFC, NSFC, NGASW,ICK,  IMC,
     9                ISK,  IJJ,  INAJ, IKK,  IKKT, IICR, IITD, IIML,
     *                ILRO, ILRI, ILDS, ILOU, IIFX, IPKK, IPKF, IPKL,
     1                IIDS, ISIT, IFS,  ILS,  IBLK, IFB,  ILB,  IKKS,
     2                ITPW, IIPV, IMO,  IIVC, JCCH, JSCH, JKCH, JPCH 

      COMMON /CONSTS/RUC, PATM, GRAV, RHOE
C
      KERR = .FALSE.
      REWIND (LINKSK)
      CALL CKLEN (LINKCK, LOUT, LENICK, LENRCK, LENCCK, IFLAG1)
      CALL MCLEN (LINKMC, LOUT, LENIMC, LENRMC, IFLAG2)
      CALL SKLEN (LINKSK, LOUT, LENISK, LENRSK, LENCSK, IFLAG3)
      IF (IFLAG1.GT.0 .OR. IFLAG2.GT.0 .OR. IFLAG3.GT.0) THEN
         KERR = .TRUE.
         RETURN
      ENDIF
C
C     real chemkin work space
      NCK   = 1
C     real transport work space
      NMC   = NCK   + LENRCK
C     real surface work space
      NSK   = NMC   + LENRMC
      NTOT  = NSK   + LENRSK
C
C     integer chemkin work space
      ICK    = 1
C     integer transport work space
      IMC    = ICK    + LENICK
C     integer surface work space
      ISK    = IMC    + LENIMC
      ITOT   = ISK    + LENISK
C
C     chemkin character work space
      JCCH  = 1
C     surf. char. work space
      JSCH  = JCCH + LENCCK
      JTOT  = JSCH + LENCSK
C
      IF (LNIPAR.GT.ITOT .AND. LNRPAR.GT.NTOT .AND. LNCWRK.GT.JTOT)
     1    THEN
         CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT, IPAR, RPAR,
     1                CWORK, IFLAG)
         KERR = KERR.OR. (IFLAG .GT. 0) 
C
         CALL CKINDX (IPAR, RPAR, MM, KKGAS, II, NFIT)
         CALL MCINIT (LINKMC, LOUT, LENIMC, LENRMC, IPAR(IMC),
     1                RPAR(NMC), IFLAG)
         KERR = KERR.OR. (IFLAG .GT. 0)
C
         CALL SKINIT (LENISK, LENRSK, LENCSK, LINKSK, LOUT, IPAR(ISK),
     1                RPAR(NSK), CWORK(JSCH), IFLAG)
         KERR = KERR.OR. (IFLAG .GT. 0)
C
         IF (KERR) RETURN
         CALL SKINDX (IPAR(ISK), MM, KKGAS, KKSURF, KKBULK, KKTOT, 
     1                NPHASE, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2                NLBULK, IISUR)
      ENDIF
C
      NATJ = KKGAS + 6
C
C     TWOPNT space
C     COMPS is the number of unknowns in the initialization problem
      COMPS = KKGAS + KKSURF
      LNTWOP = 3 + 9 * COMPS
C
C     APPORTION THE REAL SPACE
      NWT   = NTOT
      NMAS  = NWT   + KKTOT
      NYP   = NMAS  + 1
      NSPS  = NYP   + KKGAS
      NGASW = NSPS  + 1
      NSPE  = NGASW + KKGAS
      NSTS  = NSPE  + 1
      NSLB  = NSTS  + 1
      NZ1   = NSLB  + NATJ
      NZJJ  = NZ1   + COMPS
      NSRB  = NZJJ  + KKSURF
      NX    = NSRB  + NATJ
      NS    = NX    + JJ
      NSP   = NS    + NATJ*JJ + 2*KKSURF
      NYAV  = NSP   + NATJ*JJ + 2*KKSURF
      NYV   = NYAV  + KKGAS
      NCP   = NYV   + KKGAS*JJ
      NHH   = NCP   + KKGAS
      NCON  = NHH   + KKGAS
      NVIS  = NCON  + KKGAS
      NCND =  NVIS  + JJ
      NRR   = NCND  + JJ
      NY    = NRR   + JJ
      NXAV  = NY    + KKGAS
      NXM   = NXAV  + KKGAS
      NXP   = NXM   + KKGAS
      ND    = NXP   + KKGAS
      NTDC  = ND    + KKGAS*KKGAS*JJ
      NF    = NTDC  + KKGAS*JJ
      NDEP  = NF    + NATJ*JJ + 2*KKSURF
      NEKR  = NDEP  + 2*KKGAS
      NS1   = NEKR  + KKGAS**2
      NS2   = NS1   + KKGAS**2
      NQ1   = NS2   + KKGAS**2
      NQ2   = NQ1   + KKGAS**2
      NQ12  = NQ2   + KKGAS**2
      NWDT  = NQ12  + KKGAS**2
      NSDT  = NWDT  + 2*KKTOT
      NDDE  = NSDT  + NPHASE
      NSDE  = NDDE  + KKTOT
C
      NXTW  = NSDE  + NPHASE
      NFZ   = NXTW + COMPS
      NZN   = NFZ  + COMPS
      NWDI  = NZN  + COMPS
      NSDI  = NWDI + KKTOT
      NSMN  = NSDI + NPHASE
      NSMP  = NSMN + KKTOT
      NROP  = NSMP + KKTOT
      NFZN  = NROP + IISUR
      NGTP  = NFZN + COMPS
      NPR   = NGTP + 1
      NTWP  = NPR  + 1
      NABV  = NTWP + LNTWOP
      NBLW  = NABV + COMPS
      NBUF  = NBLW + COMPS
      NA    = NBUF + COMPS
      NAA   = NA   + COMPS*COMPS
      NRTP  = NAA  + COMPS*COMPS
      NDR   = NRTP + COMPS
      NDC   = NDR  + COMPS
      NXTP  = NDC  + COMPS
C
C     Save value of pointers from last call to CRPNT.  They contain
C     information input from CRRDKY that needs to be moved.
      NOLD  = NACT
      NASPO = NASP
      NRSPO = NRSP
C
      NACT  = NXTP + COMPS
      NASP  = NACT + KKTOT
      NRSP  = NASP + KKTOT
      NATS  = NRSP + KKTOT
      NRTS  = NATS + NATJ*JJ + 2*KKSURF
      NDAS  = NRTS + NATJ*JJ + 2*KKSURF
C
      IF (ICRD .EQ. 1) THEN
C        radial coordinates (or planar with a symmetry axis)
         NEQ = NATJ * JJ + KKSURF
      ELSE
C        planar coordinates
         NEQ = NATJ * JJ + 2*KKSURF
      ENDIF
      ML = MAX0(2*NATJ-1, NATJ+KKSURF-1)
      NGFC = NDAS + 40 + 9*NEQ + (3*ML+1)*NEQ+2*(NEQ/(2*ML+1)+1)
      NSFC = NGFC + 1
      NTOT = NSFC + 1
C
C     INTEGER SPACE
      IPKK   = ITOT
      IPKF   = IPKK   + NPHASE
      IPKL   = IPKF   + NPHASE
      IFS    = IPKL   + NPHASE
      ILS    = IFS    + 1
      IBLK   = ILS    + 1
      IFB    = IBLK   + 1
      ILB    = IFB    + 1
      IIDS   = ILB    + 1
      IJJ    = IIDS   + 20 + NATJ*JJ + 2*KKSURF
      ITPW   = IJJ    + 1
      IIPV   = ITPW   + 3 
      INAJ   = IIPV   + COMPS
      IKK    = INAJ   + 1
      IKKT   = IKK    + 1
      IICR   = IKKT   + 1
      IITD   = IICR   + 1
      IIVC   = IITD   + 1
      IIML   = IIVC   + 1
      IIFX   = IIML   + 1
      ILRO   = IIFX   + 1
      ILRI   = ILRO   + 1
      ILDS   = ILRI   + 1
      ILOU   = ILDS   + 1
      ISIT   = ILOU   + 1
      IKKS   = ISIT   + 1
      IMO    = IKKS   + 1
      ITOT   = IMO    + 1 - 1
C
C     APPORTION THE CHARACTER*16 SPACE
C     species names
      JKCH  = JTOT
C     phase names
      JPCH  = JKCH + KKTOT
      JTOT  = JPCH + NPHASE - 1
C
C     APPORTION THE LOGICAL SPACE
      LTOT = 1

C     On second call to CRPNT, move information to the correct position
      IF (IPCALL .GT. 1) THEN
C
         WRITE (LOUT, 7010) LNIPAR, ITOT, LNRPAR, NTOT, LNCWRK, JTOT,
     1                      LNLWRK, LTOT
C
         DO 30 K = 1, KKTOT
            RPAR(NACT + K - 1) = RPAR(NOLD + K - 1)
            RPAR(NASP + K - 1) = RPAR(NASPO + K - 1)
            RPAR(NRSP + K - 1) = RPAR(NRSPO + K - 1)
   30    CONTINUE
      ENDIF
C
C     check for enough space
      IF (ITOT.GT.LNIPAR .OR. NTOT.GT.LNRPAR .OR. JTOT.GT.LNCWRK
     1                   .OR. LTOT.GT.LNLWRK) THEN
         WRITE (LOUT, 7010) LNIPAR, ITOT, LNRPAR, NTOT, LNCWRK, JTOT,
     1                      LNLWRK, LTOT
         WRITE (LOUT, *) '  FATAL ERROR, NOT ENOUGH WORK SPACE PROVIDED'
         KERR = .TRUE.
         RETURN
      ENDIF
C
      CALL CKRP   (IPAR(ICK), RPAR(NCK), RU, RUC, PATM)
      CALL SKPKK  (IPAR(ISK), IPAR(IPKK), IPAR(IPKF), IPAR(IPKL))
      CALL SKSDEN (IPAR(ISK), RPAR(NSK), RPAR(NSDE))
      CALL SKWT   (IPAR(ISK), RPAR(NSK), RPAR(NWT))
      CALL SKSYMS (IPAR(ISK), CWORK(JSCH), LOUT, CWORK(JKCH), IERR)
      CALL SKSYMP (IPAR(ISK), CWORK(JSCH), LOUT, CWORK(JPCH), IERR)
C
      RETURN
C
C     FORMAT statements
7010  FORMAT (/,'                WORKING SPACE REQUIREMENTS',
     1        /,'                 PROVIDED        REQUIRED ',
     2        /,' INTEGER  ' , 2I15,
     3        /,' REAL     ' , 2I15,
     4        /,' CHARACTER' , 2I15,
     5        /,' LOGICAL  ' , 2I15,/)
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRPRNT (KKGAS, JJ, NATJ, ICRD, TIME, LOUT, S, Y, YP,
     1                   R, DEP, IDAS, RDAS, IPAR, RPAR, TMASS, Z1,
     2                   ZJJ, WT, ISKWRK, RSKWRK, CSKWRK, KNAM, PNAM,
     3                   NNSURF, KKBULK, SDOT, SITDOT, SDEN, DEN,
     4                   SDOTI, SITDTI, ROP, KKTOT, IISUR, KKSURF, ACT,
     5                   NFBULK, NLBULK, NFSURF, NLSURF, KKPHAS,
     6                   KFIRST, KLAST, LMOLF)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION S(NATJ,*), Z1(*), ZJJ(*), IPAR(*), RPAR(*), IDAS(*),
     1          RDAS(*), R(*), Y(*), YP(*), DEP(*), ISKWRK(*),
     2          RSKWRK(*), SDOT(*), SITDOT(*), DEN(*), ROP(*), WT(*),
     3          SDOTI(*), SITDTI(*), SDEN(*), ACT(*), KKPHAS(*),
     4          KFIRST(*), KLAST(*)
C
      CHARACTER*48 RNAM
      CHARACTER*80 ISTR
      CHARACTER*(*) CSKWRK(*), KNAM(*), PNAM(*)
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      LOGICAL LONGPT, IERR, LMOLF
C
      COMMON /LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
      COMMON /CONSTS/ RUC, PATM, GRAV, RHOE
C
      DATA DEPPRT/60.0E4/, KPRLN/5/,
     1     LT/48/, LONGPT/.FALSE./
C
      KPLN = 10
      K2 = KPLN - 3
      LS = K2 + 1
      K2 = MIN (K2, KKGAS)
C
      WRITE (LOUT, 7000)
C
      TMF = TMASS + S(NML, 1) + S(NMU, 1)
      PTORR = S(NP,1) / PATM * 760.
C
C     extract the physical coordinate (R)
      IF (ICRD .EQ. 1) THEN
C        radial coordinates
         DO 150 J = 1, JJ
           R(J) = SQRT(ABS(S(NR,J)))
  150    CONTINUE
      ELSE
C        planar coordinates
         DO 200 J = 1, JJ
            R(J) = S(NR,J)
  200    CONTINUE
      ENDIF
C
C     get statics from DASSL
      HLAST  = RDAS(7)
      NORDER = IDAS(8)
      NSTEP  = IDAS(11)
      NRES   = IDAS(12)
      NNJAC  = IDAS(13)
      NERR   = IDAS(14)
      NCON   = IDAS(15)
C
      WRITE (LOUT, 7010)  TIME, HLAST, NSTEP
      WRITE (LOUT, 7015)  NORDER, NRES, NNJAC, NERR, NCON
      WRITE (LOUT, 7080)  TMF, PTORR
C
      WRITE (LOUT, 7070) (KNAM(K), K=1,K2)
      WRITE (LOUT, 7140)
C
C     print deposition rate at lower wall only for planar
C     coordinates
      IF (ICRD .EQ. 0) WRITE (LOUT, 7120) (DEP(K),K=1,K2)
C
      WRITE (LOUT, 7150) (DEP(K), K = KKGAS+1, KKGAS+K2)
      WRITE (LOUT, 7140)
C
      DO 350 I = 1, JJ
         J = JJ - I + 1
         OMS = 0.
         DO 250 K = 1, KKGAS
            Y(K) = S(NYS+K,J)
         OMS = OMS + Y(K)
  250    CONTINUE
         OMS = 1. - OMS
C
C        compute partial pressures in torr for print
         CALL CKYTX (Y, IPAR, RPAR, YP)
         DO 300 K = 1, K2
            Y(K) = YP(K) * PTORR
  300    CONTINUE
C
         IF (LMOLF) THEN
            WRITE (LOUT, 7020) J, R(J), S(NU,J), S(NT,J), OMS,
     1                         (YP(K),K=1,K2)
         ELSE
            WRITE (LOUT, 7020) J, R(J), S(NU,J), S(NT,J), OMS,
     1                         (Y(K),K=1,K2)
         ENDIF
C
  350 CONTINUE
C
      IF (K2 .NE. KKGAS) THEN
C
         DO 500 L = LS, KKGAS, KPLN
            K1 = L
            K2 = L + KPLN - 1
            IF (K2 .GT. KKGAS) K2 = KKGAS
C
            WRITE (LOUT, 7060) (KNAM(K), K=K1,K2)
            WRITE (LOUT, 7140)
C
C           print deposition rate at lower wall only for planar
C           coordinates
            IF (ICRD .EQ. 0) WRITE (LOUT, 7130) (DEP(K), K=K1,K2)
C
            WRITE (LOUT, 7160) (DEP(K), K = KKGAS + K1, KKGAS + K2)
            WRITE (LOUT, 7140)
C
            DO 500 I = 1, JJ
               J = JJ - I + 1
               DO 400 K = 1, KKGAS
                  Y(K) = S(NYS+K,J)
  400          CONTINUE
               CALL CKYTX (Y, IPAR, RPAR, YP)
               DO 450 K = K1, K2
                  Y(K) = YP(K) * PTORR
  450          CONTINUE
               IF (LMOLF) THEN
                  WRITE (LOUT, 7020) J, R(J), (YP(K), K = K1, K2)
               ELSE
                  WRITE (LOUT, 7020) J, R(J), (Y(K), K = K1, K2)
               ENDIF
  500    CONTINUE
      ENDIF
C
  525 CONTINUE
C
      IF (NNSURF .GT. 0) THEN
         WRITE (LOUT, '(/1X,A/)') 'SURFACE SPECIES SITE FRACTIONS'
C
C        surface species
         DO 1500 N = NFSURF, NLSURF
            ILS = CKLSCH(PNAM(N))
            ISTR = ' '
            ISTR = PNAM(N)(1:ILS)
            ILS = N - 1
            IF (ICRD .EQ. 0) THEN
C              planar coordinates
               WRITE (LOUT, 7045) ILS, ISTR
            ELSE
C              radial coordinates (or planar with a symmetry axis)
               WRITE (LOUT, 7050) ILS, ISTR
            ENDIF
            IF (KKPHAS(N) .GT. 0) THEN
               DO 1520 K = KFIRST(N), KLAST(N)
                  IF (ICRD .EQ. 0) THEN
C                    planar coordinates
                     WRITE (LOUT, 7055) KNAM(K), Z1(K-KKGAS),
     1                                  ZJJ(K-KKGAS)
                  ELSE
C                    radial coordinates (or planar with a symmetry axis)
                     WRITE (LOUT, 7055) KNAM(K), ZJJ(K-KKGAS)
                  ENDIF
1520           CONTINUE
            ELSE
               WRITE (LOUT, *) ' -----NO SPECIES ON THIS SITE-----'
            ENDIF
1500     CONTINUE
      ENDIF
C
      CALL SKWT (ISKWRK, RSKWRK, WT)
C
      CALL CKYTX (Y, IPAR, RPAR, ACT)
      IF (KKSURF .GT. 0) THEN
        DO 1510 K = KFIRST(NFSURF), KLAST(NLSURF)
           ACT(K) = Z1(K-KKGAS)
 1510   CONTINUE
      ENDIF
      CALL SKSDEN (ISKWRK, RSKWRK, SDEN)
      CALL SKDEN  (S(NP,1), S(NT,1), ACT, SDEN, ISKWRK, RSKWRK, DEN)
C
      IF (ICRD.EQ.0 .AND. LONGPT) THEN
         IF (IISUR .GT. 0) THEN
            CALL SKROP (S(NP,1), S(NT,1), ACT, SDEN, ISKWRK, RSKWRK,
     1                     ROP)
            WRITE (LOUT,'(2(/1X,A),/)')
     1'CONTRIBUTION OF SURFACE REACTIONS TO SPECIES DESTRUCTION RATES',
     1'-- POSITIVE AND NEGATIVE RATES ARE EACH NORMALIZED TO ONE'
C
C           surface reaction rates
            DO 1300 L = 1, KKTOT, KPRLN
C
               K1 = L
               K2 = L+KPRLN-1
               K2 = MIN (K2, KKTOT)
               WRITE (LOUT, 7270) (KNAM(K), K=K1,K2)
C
               DO 1200 I = 1, IISUR
                  RNAM = '  '
                  CALL SKSYMR (I, LOUT, ISKWRK, RSKWRK, CSKWRK, LT,
     1                            RNAM, IERR)
                  CALL SKRATI (I, ROP, ISKWRK, RSKWRK, SDOTI, 
     1                         SITDTI)
                  WRITE (LOUT, 7280) I, RNAM, (SDOTI(K), K=K1,K2)
1200           CONTINUE
C
               WRITE (LOUT, 7290) (SDOT(K)*WT(K), K=K1,K2)
1300        CONTINUE
         ENDIF
      ENDIF
C
      IF (KKBULK .LE. 0) RETURN
C
C     print total deposition ratesfor all deposit species
      WRITE (LOUT, '(/1X,A/)') 'DEPOSITION RATE (MICRON/MIN)'
C
      CALL SKDEN (S(NP,1), S(NT,1), ACT, SDEN, ISKWRK, RSKWRK, DEN)
      DO 2000 K = KFIRST(NFBULK), KLAST(NLBULK)
         ILS = CKLSCH(KNAM(K))
         ISTR = ' '
         ISTR = KNAM(K)(1:ILS)
         IF (ICRD .EQ. 0) THEN
C           planar coordinates
            DEPLOW = SDOT(K)      *WT(K)/DEN(K) * DEPPRT
            DEPHI  = SDOT(K+KKTOT)*WT(K)/DEN(K) * DEPPRT
            ILS = K - KKGAS - KKSURF
            WRITE (LOUT, 7100)ILS,ISTR,DEPLOW,DEPHI
         ELSE
C           radial coordinates
            DEPHI  = SDOT(K+KKTOT)*WT(K)/DEN(K) * DEPPRT
            ILS = K - KKGAS - KKSURF
            WRITE (LOUT, 7100) ILS, ISTR, DEPHI
         ENDIF
2000  CONTINUE
C
      RETURN
C
C     FORMAT statements
 7000 FORMAT (1H1)
 7010 FORMAT (10X,  'DISTANCE=',  8X, E10.3, 7X, 'LAST STEP SIZE=',
     1         7X, E10.3, 7X, 'NUMBER OF STEPS=', 5X, I6)
 7015 FORMAT (10X, 'ORDER OF INTEGRATION=', I6,
     1         7X,  'NUMBER OF FUNCTION CALLS =', I6,
     2         7X,  'JACOBIAN EVALUATIONS=', I6,
     3        /10X, 'ERROR TEST FAILURES =', I6,
     4         7X,  'CONVERGENCE TEST FAILURES=', I6)
 7020 FORMAT (I4, 12E11.3)
 7045 FORMAT (' SITE ',I2,': ', A16, 5X,
     1        ' LOWER WALL            UPPER WALL')
 7050 FORMAT (' SITE ',I2,': ', A16, 5X, ' OUTER WALL')
 7055 FORMAT (4X,  A12, 14X, 2(1PE11.3, 11X))
 7060 FORMAT (/10X, 'Y(CM)', 5X, 10(A10,1X))
 7070 FORMAT (/10X, 'Y(CM)', 8X, 'U', 8X, 'T', 8X, 'OMS', 7X,
     1           8(A10,1X))
 7080 FORMAT (10X, 'MASS IN THE FLOW=', E10.3, 7X,
     1         'P(TORR)=', 14X, E10.3, 7X)
 7100 FORMAT (I4, 2X, A16, 8X, 2(1PE11.3, 11X))
 7120 FORMAT (21X,'DEPOSITION RATE--LOWER WALL', 7E11.3)
 7130 FORMAT (5X, 'LOWER WALL', 10E11.3)
 7140 FORMAT (10X,'  ')
 7150 FORMAT (21X,'DEPOSITION RATE--UPPER WALL', 7E11.3)
 7160 FORMAT (5X, 'UPPER WALL', 10E11.3)
 7270 FORMAT (/57X, 5(A10, 1X))
 7280 FORMAT (I4, 2X, A48, 5(1PE11.3))
 7290 FORMAT (' TOTAL (GM/CM**2-SEC)', 33X, 5(1PE11.3))
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRPROF (LPRO, LOUT, S, R, X, ACT, ZLOW, ZHIGH, 
     1          IPAR, RPAR, NATJ, JJ, KK, KKBULK, KKSURF, ICRD, 
     2          TMF, KERR)
C
C  START PROLOGUE
C        READ A SOLUTION PROFILE FROM LOGICAL UNIT LPRO
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION S(NATJ,*), X(*), ACT(*), ZLOW(*), ZHIGH(*), 
     1          IPAR(*), RPAR(*), R(*)
      LOGICAL KERR
      COMMON/LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
C
C      OPEN ( UNIT=LPRO, STATUS='OLD', FORM='UNFORMATTED',
C     1       FILE='cres.pro', ERR=300)
C
      KERR = .FALSE.
      READ (LPRO) ICRDIN, JJIN, NATJIN, KKIN, KKSIN, KKBIN
C
      IF (ICRD .NE. ICRDIN) THEN
         WRITE ( LOUT, *) 
     1         ' ERROR IN CRPROF: ICRD DOES NOT MATCH VALUE READ'
        KERR = .TRUE.
      ENDIF
C
      IF (JJ .NE. JJIN) THEN
         WRITE ( LOUT, *) 
     1         ' ERROR IN CRPROF: JJ DOES NOT MATCH VALUE READ'
         KERR = .TRUE.
      ENDIF
C
      IF (NATJ .NE. NATJIN) THEN
         WRITE ( LOUT, *) 
     1         ' ERROR IN CRPROF: NATJ DOES NOT MATCH VALUE READ'
        KERR = .TRUE.
      ENDIF
C
      IF (KK .NE. KKIN) THEN
         WRITE ( LOUT, *) 
     1         ' ERROR IN CRPROF: KK DOES NOT MATCH VALUE READ'
         KERR = .TRUE.
      ENDIF
C
      IF (KKSURF .NE. KKSIN) THEN
         WRITE ( LOUT, *) 
     1         ' ERROR IN CRPROF: KKSURF DOES NOT MATCH VALUE READ'
         KERR = .TRUE.
      ENDIF
C
      IF (KKBULK .NE. KKBIN) THEN
         WRITE ( LOUT, *) 
     1         ' ERROR IN CRPROF: KKBULK DOES NOT MATCH VALUE READ'
         KERR = .TRUE.
      ENDIF
C
      IF (KERR) RETURN
        IF (ICRD .EQ. 0) THEN
            READ (LPRO,ERR=400,END=400) (ZLOW(K),K=1,KKSURF),
     1                   ((S(N,J),N=1,NATJ),J=1,JJ),
     2                   (ZHIGH(K),K=1,KKSURF),  
     3                   (ACT(K),K=KK+KKSURF+1,KK+KKSURF+KKBULK)
        ELSE
            READ (LPRO,ERR=400,END=400) ((S(N,J),N=1,NATJ),J=1,JJ),
     1                   (ZHIGH(K),K=1,KKSURF),
     2                   (ACT(K),K=KK+KKSURF+1,KK+KKSURF+KKBULK)
        ENDIF
C
      IF (ICRD .EQ. 1) THEN
C        radial coordinates
         C = 4.0
      ELSE
C        planar coordinates
         C = 2.0
      ENDIF
C
      X(1) = 0.
      R(1) = 0.
      CALL CKRHOY (S(NP,1), S(NT,1), S(NY1,1), IPAR, RPAR, RHO)
C
      DO 100 J = 2, JJ
         JM1 = J - 1
         RHOM = RHO
         CALL CKRHOY (S(NP,J), S(NT,J), S(NY1,J), IPAR, RPAR, RHO)
         X(J) = X(JM1) + (S(NR,J)-S(NR,JM1)) *
     1                   (S(NU,J)*RHO+S(NU,JM1)*RHOM) / C
         IF ( ICRD .EQ. 1) THEN
            R(J) = SQRT(S(NR,J))
         ELSE
            R(J) = S(NR,J)
         ENDIF
  100 CONTINUE
C
      TMF = X(JJ)
      IF (ICRD .EQ. 1) TMF = TMF * 2.0 * ACOS(-1.0)
C
C     normalize streamfunction
      DO 200 J = 1, JJ
         X(J) = X(J) / TMF
  200 CONTINUE
C
C
      RETURN
C
  300 CONTINUE
      WRITE ( LOUT, *) ' ERROR IN CRPROF: COULD NOT OPEN UNIT LPRO'
      KERR = .TRUE.
      RETURN
  400 CONTINUE
      WRITE (LOUT, *) ' ERROR IN CRPROF: ERROR READING RECORD',
     1              ' OR END OF FILE'
      KERR = .TRUE.
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRRDKY (LIN, KKGAS, KKSURF, KKBULK, KKTOT, NFSURF,
     1                   NLSURF, NFBULK, NLBULK, KFIRST, KLAST,
     2                   KKPHAS, NPHASE, PNAM, LOUT, REAC, PRES, VEL,
     3                   GTMP, STMP, XEND, NPTS, ICRD, HITE, DX, IRST,
     4                   ITDIF, IMULT, IFIXED, STCH, ATOL, RTOL, KNAM,
     5                   Z1, ZJJ, MORD, IDOSP, XTMP, ACT, GFAC,
     6                   SFAC, ATOLSP, RTOLSP, SSABS, SSREL, TDABS,
     7                   TDREL, IPRNT, LIMIT, NINIT, IRETIR, NJAC,
     8                   ITJAC, DTMAX, DTMIN, LNOTP, LMOLF, LTHIC,
     9                   BLTK, HO, NPRNDS, GASW, TSTEP0, IVCOR, LPROF,
     *                   KERR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
      DIMENSION REAC(*), VALUE(5), Z1(*), ZJJ(*), ACT(*),
     1          KFIRST(*), KLAST(*), KKPHAS(*),
     1          ATOLSP(KKTOT), RTOLSP(KKTOT), GASW(*)
      LOGICAL NEC(10), KERR, IERR, LNOTP, LMOLF, LTHIC, LPROF
C
      CHARACTER LINE*80, KEY*4, KNAM(*)*(*), PNAM(*)*(*), 
     1          SUB*80
C
      REAL SPLX(50), SPLY(50), SPLYP(50), SPLYPP(50), SPLW(150), TMPSTR,
     1     TMPEND
C
      COMMON /SPLN/    SPLX, SPLY, SPLYP, SPLYPP,
     1                 SPLW, ISPLN, TMPSTR, TMPEND
      COMMON /CONSTS/ RUC, PATM, GRAV, RHOE
C
      DATA (NEC(I),I=1,9)/9*.FALSE./, NMAX/50/
C
C     initialize variables
      LNOTP = .FALSE.
      LMOLF = .TRUE.
      LTHIC = .FALSE.
C
      DO 10 K = 1, KKGAS
         REAC(K) = 0.
         GASW(K) = 0.
10    CONTINUE
      DO 20 K = 1, KKSURF
         Z1(K)  = 0.0
         ZJJ(K) = 0.0
   20 CONTINUE
      DO 30 K = 1, KKTOT
         ACT(K) = 0.0
   30 CONTINUE
C
      WRITE (LOUT, 9000)
C
      GFAC = 1.0
      SFAC = 1.0
      KERR  = .FALSE.
      LPROF = .FALSE.
      NREAC = 0
      NGWALL = 0
      ATOL  = 1.0E-8
      RTOL  = 1.0E-4
      MORD  = 5
      HO    = 0.0
      DX    = .5
      IRST  = 0
      ITDIF = 0
      IMULT = 0
      IVCOR = 0
      STCH  = 1.0
      IFIXED= 0
      ISPLN = 0
      IDOSP = 0
      XTMP  = 0.5
      GRAV  = 0.0
C     TWOPNT parameters
      SSABS = 1.0E-13
      SSREL = 1.0E-14
      TDABS = 1.0E-12
      TDREL = 1.0E-4
      TSTEP0 = 1.0E-6
      LIMIT = 100
      IPRNT = 22
      NINIT = 0
      NPRNDS = 0
      IRETIR = 50
      NJAC = 20
      ITJAC = 20
      DTMAX = 1.0E-4
      DTMIN = 1.0E-10
      BLTK  = 0.0
C
      DO 35 I = 1, KKTOT
        ATOLSP(I) = ATOL
35      RTOLSP(I) = RTOL
C
C
90    CONTINUE
C     read input line
      KEY = ' '
      LINE   = ' '
      IERR   = .FALSE.
      READ  (LIN,  7000, END=6000) KEY, LINE
      CALL CKDTAB (LINE)
C
C     is this a keyword comment?
      IF (KEY(1:1) .EQ. '.' .OR. KEY(1:1) .EQ. ' ' .OR.
     1    KEY(1:1) .EQ. '!') GO TO 90
      WRITE (LOUT, 8000) KEY, LINE
      IF (KEY(1:1) .EQ. '/') GO TO 90
C
C-----METHOD OPTIONS KEYWORDS----------------------------
C
      IF (KEY .EQ. 'ATOL') THEN
C        absolute error tolerance
C
         ATOLO = ATOL
         CALL CKXNUM (LINE, 1, LOUT, NVAL, ATOL, IERR)
         DO 41 I = 1, KKTOT
           IF (ATOLSP(I) .EQ. ATOLO) ATOLSP(I) = ATOL
41       CONTINUE
C
      ELSEIF (KEY .EQ. 'RTOL') THEN
C        relative rror tolerance
C
         RTOLO = RTOL
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RTOL, IERR)
         DO 42 I = 1, KKTOT
           IF (RTOLSP(I) .EQ. RTOLO) RTOLSP(I) = RTOL
42       CONTINUE
C
      ELSEIF (KEY .EQ. 'IERW') THEN
C        individual error tolerances
C
         CALL SKSNUM(LINE, 2, LOUT, KNAM, KKTOT, PNAM, NPHASE,
     1         KKPHAS, KSP, NKF, NVAL, VALUE, IERR)
         IF (.NOT. IERR) THEN
           ATOLSP(KSP) = VALUE(1)
           RTOLSP(KSP) = VALUE(2)
         ENDIF
C
      ELSEIF (KEY(1:2) .EQ. 'HO') THEN
C        value of the initial time step
         CALL CKXNUM (LINE, 1, LOUT, NVAL, HO, IERR)
C
      ELSEIF (KEY .EQ. 'NPTS') THEN
C        number of grid points
         NEC(6)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         NPTS   = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'MORD') THEN
C        maximum order of integration
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         MORD   = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'STCH') THEN
C        coordinate stretching parameter
         CALL CKXNUM (LINE, 1, LOUT, NVAL, STCH, IERR)
C
      ELSEIF (KEY .EQ. 'XEND') THEN
C        end of channel
         NEC(5)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XEND, IERR)
C
      ELSEIF (KEY .EQ. 'PROF') THEN
C        read solution profile from unit LPRO
         LPROF  = .TRUE.
C
C-----TWOPNT KEYWORDS---------------------------------
C
      ELSEIF (KEY .EQ. 'NOTP') THEN
C        do not call twopnt
         LNOTP = .TRUE.
C
      ELSEIF (KEY .EQ. 'TWAB') THEN
C        absolute error tolerance
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SSABS, IERR)
C
      ELSEIF (KEY .EQ. 'TWRE') THEN
C        relative error tolerance
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SSREL, IERR)
C
      ELSEIF (KEY .EQ. 'TWTA') THEN
C        time-step absolute error tolerance
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TDABS, IERR)
C
      ELSEIF (KEY .EQ. 'TWTR') THEN
C        time-step relative error tolerance
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TDREL, IERR)
C
      ELSEIF (KEY .EQ. 'TWPR') THEN
C        print flag
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IPRNT = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'TWST') THEN
C        number of time steps before trying another Newton step
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         LIMIT = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'ISTP') THEN
C        number of initial time steps before the first Newton step
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         NINIT = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'IRET') THEN
C        retirement age of old time step
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IRETIR = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'NJAC') THEN
C        retirement age of Jacobian during steady-state Newton
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         NJAC = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'TJAC') THEN
C        retirement age of Jacobian during time-stepping
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         ITJAC = INT(VALUE(1))
C
      ELSEIF (KEY .EQ. 'DTMX') THEN
C        maximum time-step
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DTMAX, IERR)
C
      ELSEIF (KEY .EQ. 'DTMN') THEN
C        maximum time-step
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DTMIN, IERR)
C
      ELSEIF (KEY .EQ. 'STP0') THEN
C        initial time-step
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TSTEP0, IERR)
C
      ELSEIF (KEY .EQ. 'MOLF') THEN
C        print mol;e fractions
         LMOLF = .TRUE.
C
      ELSEIF (KEY .EQ. 'PARP') THEN
C        print partial pressures
         LMOLF = .FALSE.
C
C-----REACTOR DEFINITION KEYWORDS-------------------------
C
      ELSEIF (KEY .EQ. 'PRES') THEN
C        pressure
         NEC(2)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, PRES, IERR)
         PRES = PRES*PATM
C
      ELSEIF (KEY .EQ. 'GTMP') THEN
C        gas temperature (Kelvin)
         NEC(4)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, GTMP, IERR)
C
      ELSEIF (KEY .EQ. 'STMP') THEN
C        susceptor temperature (Kelvin)
         NEC(9)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, STMP, IERR)
C
      ELSEIF (KEY .EQ. 'XTMP') THEN
C         distance for temperature rise (cm)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XTMP, IERR)
C
      ELSEIF (KEY .EQ. 'TSPL') THEN
C        read specified temperature profile (X,T) pairs
C
         IDOSP = 1
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         IF (ISPLN+1 .GT. NMAX) THEN
               WRITE (LOUT, *)
     1          ' ERROR... THE PROBLEM IS ONLY DIMENSIONED FOR ', NMAX,
     2          ' (X,T) PAIRS'
               IERR = .TRUE.
         ELSE
            ISPLN = ISPLN + 1
C*****precision > quad
C            SPLX(ISPLN) = SNGLQ(VALUE(1))
C            SPLY(ISPLN) = SNGLQ(VALUE(2))
C*****END precision > quad
C*****precision > double
            SPLX(ISPLN) = SNGL(VALUE(1))
            SPLY(ISPLN) = SNGL(VALUE(2))
C*****END precision > double
C*****precision > single
C            SPLX(ISPLN) = VALUE(1)
C            SPLY(ISPLN) = VALUE(2)
C*****END precision > single
C
            IF (ISPLN .EQ. 1) THEN
               GTMP = VALUE(2)
            ELSE
               STMP = VALUE(2)
            ENDIF
         ENDIF
C
      ELSEIF (KEY .EQ. 'FIXT') THEN
C        fix the top wall temperature at the inlet gas temperature
         IFIXED = 1
C
      ELSEIF (KEY .EQ. 'SYMT') THEN
C        set the top wall temperature equal to the bottom wall
C        temperature
         IFIXED = 2
C
      ELSEIF (KEY(1:3) .EQ. 'VEL') THEN
C        velocity
         NEC(3)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VEL, IERR)
C
      ELSEIF (KEY .EQ. 'ICRD') THEN
C        coordinate system, PLAN = planar
C                           RAD  = radial coordinates
C                           SYMC = planar channel w/ symmetry axis
C
         NEC(7)  = .TRUE.
         CALL CKSUBS (LINE, LOUT, 1, SUB, NFOUND, IERR)
         IF (NFOUND .LE. 0) THEN
           IERR = .TRUE.
           WRITE(LOUT,*) 
     1         ' RDKEY ERROR: No argument given for ICRD keyword'
         ELSE
           IF (NFOUND .GT. 1) THEN
             WRITE(LOUT,*)
     1         'RDKEY WARNING: More than one argument on ICRD keyword ',
     2         'will be ignored.'
           ENDIF
           IF (SUB .EQ. 'PLAN') THEN
             ICRD = 0
           ELSE IF (SUB(1:3) .EQ. 'RAD') THEN
             ICRD = 1
           ELSE IF (SUB .EQ. 'SYMC') THEN
             ICRD = 2
           ELSE
             WRITE(LOUT,*)
     1         'RDKEY ERROR: Unknown parameter for ICRD keyword'
             IERR = .TRUE.
           ENDIF
         ENDIF
C
      ELSEIF (KEY .EQ. 'GRAV') THEN
C        acceleration of gravity: > 0, flow is up,
C                                 < 0, flow is down,
C                                 = 0, flow is horizontal
         CALL CKXNUM (LINE, 1, LOUT, NVAL, GRAV, IERR)
C
      ELSEIF (KEY .EQ. 'HITE') THEN
C        channel dimension
         NEC(8)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, HITE, IERR)
C
      ELSEIF (KEY .EQ. 'BLTK') THEN
C        boundary layer thickness
         LTHIC = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, BLTK, IERR)
C
      ELSEIF (KEY .EQ. 'REAC') THEN
C        reactant
C
         NEC(1) = .TRUE.
         CALL CKSNUM (LINE, 1, LOUT, KNAM, KKGAS, KSP, NVAL, VALUE,
     1                IERR)
         IERR = IERR.OR. NREAC+1.GT.KKGAS
         IF (.NOT. IERR) THEN
            NREAC       = NREAC+1
            REAC(KSP) = VALUE(1)
         ENDIF
C
      ELSEIF (KEY .EQ. 'GASW') THEN
C        gas mole fractions at the wall
C
         CALL CKSNUM (LINE, 1, LOUT, KNAM, KKGAS, KSP, NVAL, VALUE,
     1                IERR)
         IERR = IERR.OR. NGWALL+1.GT.KKGAS
         IF (.NOT. IERR) THEN
            NGWALL       = NGWALL+1
            GASW(KSP) = VALUE(1)
         ENDIF
C
      ELSEIF (KEY .EQ. 'SURF') THEN
C        surface species initial fractions
C
         IF (KKSURF .LE. 0) THEN
            IERR = .TRUE.
            WRITE (LOUT, *) ' Error...no site-phase species exist'
         ELSE
           CALL SKSNUM (LINE, 1, LOUT, KNAM, KKTOT, PNAM, NPHASE,
     1                  KKPHAS, KSP, NKF, NVAL, VALUE, IERR)
           IF (IERR) THEN
           ELSEIF (KSP.LT.KFIRST(NFSURF).OR.KSP.GT.KLAST(NLSURF)) THEN
              WRITE (LOUT, *)
     1        ' Error...SURF must be site-phase species '
              IERR = .TRUE.
           ELSE
              IF (NKF .GT. 1) WRITE (LOUT, *)
     1        ' Warning...non-unique species name given '
              Z1(KSP - KKGAS)   = VALUE(1)
              ZJJ(KSP- KKGAS)  = VALUE(1)
           ENDIF
         ENDIF
C
      ELSEIF (KEY(1:3) .EQ. 'ACT') THEN
C        bulk species activities
C
         IF (KKBULK .LE. 0) THEN
            IERR = .TRUE.
            WRITE (LOUT, *) ' Error...no bulk-phase species exist'
         ELSE
            CALL SKSNUM (LINE, 1, LOUT, KNAM, KKTOT, PNAM, NPHASE,
     1                   KKPHAS, KSP, NKF, NVAL, VALUE, IERR)
            IF (IERR) THEN
            ELSEIF (KSP.LT.KFIRST(NFBULK).OR.KSP.GT.KLAST(NLBULK)) THEN
               WRITE (LOUT, *)
     1         ' Error...ACT must be bulk-phase species'
               IERR = .TRUE.
            ELSE
               ACT(KSP) = VALUE(1)
            ENDIF
         ENDIF
C
      ELSEIF (KEY .EQ. 'GFAC') THEN
C        factor for gas-phase rate constants
         CALL CKXNUM (LINE, 1, LOUT, NVAL, GFAC, IERR)
C
      ELSEIF (KEY .EQ. 'SFAC') THEN
C        factor for surface-phase rate constants
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SFAC, IERR)
C
C-----TRANSPORT OPTIONS KEYWORDS---------------------------
C
      ELSEIF (KEY .EQ. 'MULT') THEN
C        multicomponent diffusion included
         IMULT   = 1
C
      ELSEIF (KEY(1:3) .EQ. 'MIX') THEN
C        mixture-averaged diffusion used
         IMULT   = 0
C
      ELSEIF (KEY .EQ. 'TDIF') THEN
C        thermal diffusion included
         ITDIF    = 1
C
      ELSEIF (KEY .EQ. 'VCOR') THEN
C        correction velocity included
         IVCOR    = 1
C
C-----PRINTING AND RESTARTING KEYWORDS-----------------------
C
      ELSEIF (KEY(1:2) .EQ. 'DX') THEN
C        output interval
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DX, IERR)
C
      ELSEIF (KEY .EQ. 'PRND') THEN
C        print diagnostic information from CRESLAF
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE(1), IERR)
         NPRNDS = INT(VALUE(1))
         IF (NPRNDS .GT. 2) NPRNDS = 2
         IF (NPRNDS .LT. 0) NPRNDS = 0
C
      ELSEIF (KEY .EQ. 'IRST') THEN
C        restart flag, 1=restart, 0=new problem
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IRST    = INT(VALUE(1))
C
      ELSEIF (KEY(1:3) .EQ. 'END') THEN
C        end of problem input
         KERR = KERR.OR.IERR
         GO TO 6000
C
      ELSE
C        to get here, an invalid keyword was read
         WRITE (LOUT,*) ' ERROR...ILLEGAL KEYWORD'
         KERR = .TRUE.
C
      ENDIF
C
C     go back up and read the next line
      KERR = KERR.OR.IERR
      GO TO 90
C
6000  CONTINUE
C
C     check the reactant and surface sums
C
      IF ( NGWALL .EQ. 0) THEN
         DO 6050 K = 1, KKGAS
            GASW(K) = REAC(K)
 6050    CONTINUE
      ENDIF
C
      SUMR       = 0.
      SUMG       = 0.
      DO 6100 K = 1, KKGAS
         SUMR = SUMR + REAC(K)
         SUMG = SUMG + GASW(K)
6100  CONTINUE
C
C     normalize reactant fractions
      DO 6200 K = 1, KKGAS
         REAC(K) = REAC(K)/SUMR
         GASW(K) = GASW(K)/SUMG
6200  CONTINUE
      IF (ABS(SUMR-1.0) .GT. 1.E-6) THEN
         WRITE (LOUT, *)
     1                ' CAUTION...REACTANT FRACTIONS SUM TO ', SUMR
      ENDIF
C
      IF (ABS(SUMG-1.0) .GT. 1.E-6) THEN
         WRITE (LOUT, *)
     1                ' CAUTION...WALL MASS FRACTIONS SUM TO ', SUMG
      ENDIF
C
      IF (KKSURF.GT.0) THEN
         DO 6280 N = NFSURF, NLSURF
            SUMS = 0.
            DO 6120 K = KFIRST(N), KLAST(N)
               SUMS = SUMS + Z1(K-KKGAS)
6120        CONTINUE
C           NORMALIZE REACTANT FRACTIONS
            DO 6250 K = KFIRST(N), KLAST(N)
               Z1(K-KKGAS)  = Z1(K-KKGAS) /SUMS
               ZJJ(K-KKGAS) = ZJJ(K-KKGAS)/SUMS
6250        CONTINUE
            IF (ABS(SUMS-1.0) .GT. 1.E-6) THEN
               WRITE (LOUT, *)
     1                   ' CAUTION...SURFACE FRACTIONS SUM TO ', SUMS
            ENDIF
 6280    CONTINUE
      ENDIF
C
      IF (KKBULK.GT.0) THEN
C        normalize bulk activities
         DO 6300 N = NFBULK, NLBULK
            SUMA = 0.0
            DO 6350 K =  KFIRST(N), KLAST(N)
               SUMA = SUMA + ACT(K)
 6350       CONTINUE
C
            IF (SUMA .GT. 0) THEN
               DO 6360 K = KFIRST(N), KLAST(N)
                  ACT(K) = ACT(K) / SUMA
 6360          CONTINUE
            ENDIF
            IF (ABS(SUMA-1.0) .GT. 1.E-6) WRITE (LOUT,'(1X,A,E13.5)')
     1      ' CAUTION...BULK ACTIVITIES SUM TO ',SUMA
 6300    CONTINUE
      ENDIF
C
      IF (IDOSP .EQ. 1) THEN
         TMPSTR = SPLY(1)
         TMPEND = SPLY(ISPLN)
      ENDIF
C
C     check for necessary input
      IF (.NOT. NEC(1)) THEN
         WRITE (LOUT, *)
     1                ' ERROR...NO REACTANTS GIVEN'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(2)) THEN
         WRITE (LOUT, *) ' ERROR...PRESSURE NOT GIVEN'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(3)) THEN
         WRITE (LOUT, *) ' ERROR...VELOCITY NOT GIVEN'
         KERR = .TRUE.
      ENDIF
C
      IF (IDOSP .EQ. 0) THEN
         IF (.NOT. NEC(4))THEN
            WRITE (LOUT, *) ' ERROR...GAS TEMPERTURE NOT GIVEN'
            KERR = .TRUE.
         ENDIF
         IF (.NOT. NEC(9)) THEN
            WRITE (LOUT, *) ' ERROR...SUSCEPTOR TEMPERATURE NOT GIVEN'
            KERR = .TRUE.
         ENDIF
      ENDIF
C
      IF (.NOT. NEC(5)) THEN
         WRITE (LOUT, *) ' ERROR..."XEND" NOT GIVEN '
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(6)) THEN
         WRITE (LOUT, *) ' ERROR...NUMBER OF GRID POINTS NOT GIVEN'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(7)) THEN
         WRITE (LOUT, *) ' ERROR..."ICRD" NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(8)) THEN
         WRITE (LOUT, *) ' ERROR...CHANNEL DIMENSION NOT GIVEN'
         KERR = .TRUE.
      ENDIF
C
C     check for valid value of MORD
      IF (MORD.LT.1  .OR. MORD .GT. 5) THEN
         WRITE (LOUT, *) ' CAUTION...MUST HAVE 1.GE.MORD.LE.5 '
         WRITE (LOUT, *) '           CONTINUING WITH MORD SET TO 5'
         MORD = 5
      ENDIF
C
C     write out error weights to be used
      WRITE(LOUT,9500) ATOL, RTOL
      DO 34 I = 1 , KKTOT
34      WRITE(LOUT,9600) KNAM(I),ATOLSP(I),RTOLSP(I)
C
      RETURN
C
C     FORMATS
7000  FORMAT (A4, A)
8000  FORMAT (10X, A4, A76)
9000  FORMAT (////10X, ' KEYWORD INPUT', /)
9500  FORMAT (/' ERROR TOLERANCES TO BE USED:'/
     1 ,10X,'(default values for all other species: ATOL = '
     1 ,G10.4,'RTOL = ',G10.4//5X,'SPECIES      ATOL       '
     1 ,' RTOL')
9600  FORMAT (2X,A16,4X,G10.4,4X,G10.4)
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRRSTR (KKGAS, JJ, NATJ, LREST, LSAVE, LDAS, TIME,
     1                   S, SP, SLB, SRB, R, X, IDAS, RDAS, ICRD,
     2                   DEP, GTMP, STMP, TMASS, ITDIF, IMULT, IFIXED,
     3                   KKSURF, MORD, ACT, IVCOR)
C
C  START PROLOGUE
C        RESTART THE CODE. READ PREVIOUS SOLUTION FROM UNIT LREST,
C        COPY IT TO UNIT LSAVE, AND READ DASSL WORK SPACE.
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION S(*), SP(*), SLB(*), SRB(*), X(*), IDAS(*), RDAS(*),
     1          DEP(*), R(*), ACT(*)
C
      COMMON/LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
C
C     the following parameters must match the previous run or else
C     they will be overwritten.
C
      READ (LREST) ICRD, ITDIF, IMULT, IVCOR, IFIXED, KKGAS, KKSURF, 
     1            KKBULK, JJ, NATJ, GTMP, STMP, TMASS, MORD
      WRITE (LSAVE) ICRD, ITDIF, IMULT, IVCOR, IFIXED, KKGAS, KKSURF,
     2            KKBULK, JJ, NATJ, GTMP, STMP, TMASS, MORD
C
      IF (ICRD .EQ. 0) THEN
C        planar coordinates
         NEQ = NATJ * JJ + 2*KKSURF
      ELSE
C        radial coordinates (or planar with a symmetry axis)
         NEQ = NATJ * JJ + KKSURF
      ENDIF
      ML = MAX0(2*NATJ-1, NATJ+KKSURF-1)
C              
C     read solution
 1000 CONTINUE
      READ (LREST) IFLAG
C
      IF (IFLAG .NE. -1) THEN
C
         READ (LREST) TIME, (S(N), N = 1, NEQ), (SP(N), N = 1, NEQ),
     1               (X(J), J = 1, JJ), (DEP(K), K = 1,2*KKGAS),
     2               (ACT(K), K=KKGAS+KKSURF+1,KKGAS+KKSURF+KKBULK)
C
         WRITE (LSAVE) IFLAG
         WRITE (LSAVE) TIME, (S(N), N = 1, NEQ), (SP(N), N = 1, NEQ),
     1                 (X(J), J = 1, JJ), (DEP(K), K = 1,2*KKGAS),
     2                 (ACT(K), K=KKGAS+KKSURF+1,KKGAS+KKSURF+KKBULK)
         GO TO 1000
      ENDIF
C
C     have read complete solution; now read DASSL workspace
      LIDWRK = 20 + NEQ
      LDWRK = 40 + 9*NEQ + (3*ML+1)*NEQ + 2*(NEQ/(2*ML+1)+1)
C
      READ (LREST) (IDAS(I), I = 1, LIDWRK)
      READ (LREST) (RDAS(I), I = 1, LDWRK)
C
C     save most recent DASSL work space on unit LDAS
      WRITE (LDAS) (IDAS(I), I = 1, LIDWRK)
      WRITE (LDAS) (RDAS(I), I = 1, LDWRK)
      REWIND LDAS
C
      IF (ICRD .EQ. 0) THEN
C        planar coordinates
         DO 2050 J = 1, JJ
            IOFF = KKSURF
            IRJ  = IOFF + (J-1)*NATJ + NR
            R(J) = S(IRJ)
 2050    CONTINUE
      ELSE IF (ICRD .EQ. 1) THEN
C        radial coordinates
         DO 2075 J = 1, JJ
            IOFF = 0
            IRJ  = IOFF + (J-1)*NATJ + NR
            R(J) = SQRT(ABS(S(IRJ)))
 2075    CONTINUE
C
      ELSE
C        planar with a symmetry axis
         DO 2080 J = 1, JJ
            IOFF = 0
            IRJ  = IOFF + (J-1)*NATJ + NR
            R(J) = S(IRJ)
 2080    CONTINUE
      ENDIF
C
C     set up boundary conditions
      DO 2100 N = 1, NATJ
         ILB    = IOFF + N
         SLB(N) = S(ILB)
         IRB    = IOFF + N + NATJ*(JJ-1)
         SRB(N) = S(IRB)
 2100 CONTINUE
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRRTOX (ICRD, NATJ, JJ, R, S, IPAR, RPAR, X)
C
C  START PROLOGUE
C     CONVERT NORMAL COORDINATE (R) TO STREAMFUNCTION (X)
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION S(NATJ,*), R(*), X(*), IPAR(*), RPAR(*)
C
      COMMON /LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
C
      IF (ICRD .EQ. 1) THEN
C        radial coordinates
         C = 4.0
      ELSE
C        planar coordinates
         C = 2.0
      ENDIF
C
      X(1) = 0.
      CALL CKRHOY (S(NP,1), S(NT,1), S(NY1,1), IPAR, RPAR, RHO)
C
      DO 50 J = 2, JJ
         JM1 = J - 1
         RHOM = RHO
         CALL CKRHOY (S(NP,J), S(NT,J), S(NY1,J), IPAR, RPAR, RHO)
         X(J) = X(JM1) + (S(NR,J)-S(NR,JM1)) *
     1                   (S(NU,J)*RHO+S(NU,JM1)*RHOM) / C
 50   CONTINUE
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRRUN (KKGAS, KKTOT, JJ, NATJ, LOUT, ICRD, IRST, ITDIF,
     1                  IMULT, TIME, S, SP, F, DX, XEND, WT, Y, YP,
     2                  XAV, YAV, YV, X, R, D, TDC, IPAR, RPAR,
     3                  IMCWRK, RMCWRK, DEP, KKPHAS, KFIRST, KLAST,
     4                  IDAS, RDAS, LSAVE, LDAS, ATOL, RTOL, GTMP,
     5                  STMP, XM, XP, TMASS, CON, IFIXED, ISKWRK,
     6                  RSKWRK, KNAM, PNAM, SDOT, SITDOT, NNSURF,
     7                  KKBULK, DEN, SDEN, XSPSTR, XSPEND, KKSURF,
     8                  SDOTI, SITDTI, ROP, IISUR, CSKWRK, MORD, NFSURF,
     9                  NLSURF, NNBULK, NFBULK, NLBULK, ACT,
     *                  ATOLSP, RTOLSP, ATOLY, RTOLY, GFAC, SFAC, LMOLF,
     1                  HO, NPRNDS, IVCOR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) KNAM(*), PNAM(*), CSKWRK(*)
C
      DIMENSION S(*), SP(*), F(*), IPAR(*), RPAR(*), X(*), R(*),
     1          IDAS(*), RDAS(*), DEP(*), XM(*), XP(*), IMCWRK(*),
     2          RMCWRK(*), WT(*), XAV(*), YAV(*), YV(KKGAS,*),
     3          TDC(KKGAS,*), INFO(15), Y(*), YP(*), CON(*), ISKWRK(*),
     4          RSKWRK(*), SDOT(*), SITDOT(*), DEN(*), SDEN(*),
     1          SDOTI(*), SITDTI(*),
     5          ROP(*), D(KKGAS,KKGAS,*), ACT(*), KKPHAS(*), KFIRST(*),
     6          KLAST(*),
     7          ATOLSP(KKTOT), RTOLSP(KKTOT), ATOLY(*), RTOLY(*)
C
      LOGICAL LMOLF
      COMMON /LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
C
      EXTERNAL CRFCN
C
      DO 5 I = 1, 15
         INFO(I) = 0
    5 CONTINUE
C
C     set up integration control parameters
      IF (ICRD .EQ. 0) THEN
C        planar coordinates
         NEQ = NATJ * JJ + 2*KKSURF
         IOFF = KKSURF
      ELSE
C        radial coordinates (or planar with a symmetry axis)
         NEQ = NATJ * JJ + KKSURF
         IOFF = 0
      ENDIF
      ML = MAX0(2*NATJ-1, NATJ+KKSURF-1)
C
C     set vectors for ATOL and RTOL
      INFO(2) = 1
      DO 31 I = 1 , NEQ
         ATOLY(I) = ATOL
         RTOLY(I) = RTOL
   31 CONTINUE
      IF (ICRD .EQ. 0) THEN
C        planar coordinates
         IF (NFSURF .GT. 0) THEN
            DO 32 I = KFIRST(NFSURF),KLAST(NLSURF)
               ATOLY(I-KKGAS) = ATOLSP(I)
               RTOLY(I-KKGAS) = RTOLSP(I)
   32       CONTINUE
         ENDIF
      ENDIF
C
      DO 33 J = 1 , JJ
         JSKIP = IOFF + (J-1)*NATJ + NYS
         DO 33 I = 1 , KKGAS
            ATOLY(I+JSKIP) = ATOLSP(I)
            RTOLY(I+JSKIP) = RTOLSP(I)
   33 CONTINUE
C
      IF (NFSURF .GT. 0) THEN
C         JSKIP = KKGAS + IOFF + NATJ*JJ
         JSKIP = -KKGAS + IOFF + NATJ*JJ
         DO 34 I = KFIRST(NFSURF), KLAST(NLSURF)
C            ATOLY(I-JSKIP) = ATOLSP(I)
C            RTOLY(I-JSKIP) = RTOLSP(I)
            ATOLY(I+JSKIP) = ATOLSP(I)
            RTOLY(I+JSKIP) = RTOLSP(I)
   34    CONTINUE
      ENDIF

C     for printout only at tout
      INFO(3) = 0
C      for printout every step
C      INFO(3) = 1
      INFO(6)=1
C     INFO(10)=1
      IDAS(1) = ML
      IDAS(2) = ML
      INFO(9) = 1
      IDAS(3) = MORD

C     Set the initial time step as per user's request
      IF (HO .GT. 0.0) THEN
        INFO(8) = 1
        RDAS(3) = HO
      ENDIF

C     information printed out at every error
      INFO(12) = NPRNDS
C
      LIDWRK = 20 + NEQ
      LDWRK = 40+9*NEQ+(3*ML+1)*NEQ+2*(NEQ/(2*ML+1)+1)
C
C     test for restart
      IF (IRST .EQ. 0) THEN
C        new problem; re-initialize some arrays for DASSL
         RDAS(7)  = 0.0
         IDAS(8)  = 0
         IDAS(11) = 0
         IDAS(12) = 0
         IDAS(13) = 0
         IDAS(14) = 0
         IDAS(15) = 0
         WRITE (LSAVE) ICRD, ITDIF, IMULT, IVCOR, IFIXED, KKGAS, 
     1                 KKSURF,KKBULK, JJ, NATJ, GTMP, STMP, 
     2                 TMASS, MORD
         TIME = 0.
C
C        set the marching derivatives (SP) at time=0.
         DO 100 N = 1, NEQ
            SP(N) = 0.
  100    CONTINUE
C
         IRES = 0
         CALL CRFCN (TIME, S, SP, F, IRES, RPAR, IPAR)
C
         JJM1 = JJ-1
C        set SP for nodes 2 through JJ-1
         DO 200 N = IOFF+NATJ+1, IOFF+NATJ*JJM1
            SP(N) = -F(N)
  200    CONTINUE
C
         DO 300 J = 2, JJM1
            IJPT = IOFF + NATJ*(J-1)
            SP(NP+IJPT)  = 0.0
            SP(NR+IJPT)  = 0.0
            SP(NML+IJPT) = -F(IOFF+NML)
            SP(NMU+IJPT) = -F(IOFF+(JJM1)*NATJ+NMU)
  300    CONTINUE
C
         SP(IOFF+NML)             = - F(IOFF+NML)
         SP(IOFF+(JJM1)*NATJ+NML) = - F(IOFF+NML)
         SP(IOFF+NMU)             = - F(IOFF+(JJM1)*NATJ+NMU)
         SP(IOFF+(JJM1)*NATJ+NMU) = - F(IOFF+(JJM1)*NATJ+NMU)
C
      ELSE
C        restarting problem
         INFO(1) = 1
      ENDIF
      XOUT = TIME
C
C-----INTEGRATION LOOP
C
 1000 CONTINUE
C
      XOUT = XOUT + DX
C
      IF (ICRD .EQ. 0) THEN
C        compute deposition rates at lower boundary only for planar
C        coordinates
         ILOW = 1
         IOFF = KKSURF
C
         CALL CRFLUX (KKGAS, KKSURF, KKTOT, JJ, NATJ, ICRD, S(IOFF+1),
     1                WT, IPAR, RPAR, CON, DEP, ILOW, S(1), ISKWRK,
     2                RSKWRK, SDOT, SITDOT, SDEN, ACT, NFSURF, NLSURF,
     3                KFIRST, KLAST, SFAC)
      ELSE
C        radial coordinates (or planar with a symmetry axis)
         IOFF = 0
      ENDIF
C
C     compute deposition rates at upper boundary
      ILOW = 0
C
      CALL CRFLUX (KKGAS, KKSURF, KKTOT, JJ, NATJ, ICRD, S(IOFF+1),
     1             WT, IPAR, RPAR, CON, DEP(KKGAS+1), ILOW,
     2             S(IOFF+NATJ*JJ+1), ISKWRK, RSKWRK, SDOT(KKTOT+1),
     3             SITDOT, SDEN, ACT, NFSURF, NLSURF, KFIRST, KLAST,
     4             SFAC)
C
      CALL CRPRNT (KKGAS, JJ, NATJ, ICRD, TIME, LOUT, S(IOFF+1), Y, YP,
     1             R, DEP, IDAS, RDAS, IPAR, RPAR, TMASS, S(1),
     2             S(IOFF+NATJ*JJ+1), WT, ISKWRK, RSKWRK, CSKWRK, KNAM,
     3             PNAM, NNSURF, KKBULK, SDOT, SITDOT, SDEN, DEN, SDOTI,
     4             SITDTI, ROP, KKTOT, IISUR, KKSURF, ACT, NFBULK,
     5             NLBULK, NFSURF, NLSURF, KKPHAS, KFIRST, KLAST, LMOLF)
C
      CALL CRSAVE (KKGAS, KKSURF, KKBULK, JJ, TIME, NATJ, S, SP, X,
     2             DEP, IDAS, RDAS, LSAVE, LDAS, ICRD, MORD, ACT)
C
  500 CONTINUE
C
C*****precision > single
C      CALL SDASSL(CRFCN, NEQ, TIME, S, SP, XOUT, INFO, RTOLY, ATOLY,
C     1           IDID, RDAS, LDWRK, IDAS, LIDWRK, RPAR, IPAR, JAC)
C*****END precision > single
C*****precision > double
      CALL DDASSL(CRFCN, NEQ, TIME, S, SP, XOUT, INFO, RTOLY, ATOLY,
     1            IDID, RDAS, LDWRK, IDAS, LIDWRK, RPAR, IPAR, JAC)
C*****END precision > double
C
       IF (IDID .EQ. -1) THEN
           WRITE (LOUT, 9000) IDID
           INFO(1)=1
           GO TO 500
       ENDIF
C
      IF (IDID .LT. 0) THEN
C        trouble with DASSL; print statistics, then quit
         WRITE (LOUT, 9010) IDID
C
C        get statistics from DASSL
         HLAST  = RDAS(7)
         NORDER = IDAS(8)
         NSTEP  = IDAS(11)
         NRES   = IDAS(12)
         NNJAC   = IDAS(13)
         NERR   = IDAS(14)
         NCON   = IDAS(15)
C
         WRITE (LOUT, 9015)
         WRITE (LOUT, 9020)  TIME, HLAST, NSTEP
         WRITE (LOUT, 9030)  NORDER, NRES, NNJAC, NERR, NCON
C
         RETURN
      ENDIF
C
      IF (TIME .LT. XEND) GO TO 1000
C
C-----FINISHED INTEGRATION -------------
C
      IF (ICRD .EQ. 0) THEN
C        compute deposition rates at lower boundary only for
C        planar coordinates
C
         ILOW = 1
         IOFF = KKSURF
C
         CALL CRFLUX (KKGAS, KKSURF, KKTOT, JJ, NATJ, ICRD, S(IOFF+1),
     1                WT, IPAR, RPAR, CON, DEP, ILOW, S(1), ISKWRK,
     2                RSKWRK, SDOT, SITDOT, SDEN, ACT, NFSURF, NLSURF,
     3                KFIRST, KLAST, SFAC)
      ELSE
C        radial coordinates (or planar with a symmetry axis)
         IOFF = 0
      ENDIF
C
C     compute deposition rates at upper boundary
      ILOW = 0
C
      CALL CRFLUX (KKGAS, KKSURF, KKTOT, JJ, NATJ, ICRD, S(IOFF+1), WT,
     1             IPAR, RPAR, CON, DEP(KKGAS+1), ILOW,
     2             S(IOFF+NATJ*JJ+1), ISKWRK, RSKWRK, SDOT(KKTOT+1),
     3             SITDOT, SDEN, ACT, NFSURF, NLSURF, KFIRST, KLAST,
     4             SFAC)
C
      CALL CRPRNT (KKGAS, JJ, NATJ, ICRD, TIME, LOUT, S(IOFF+1), Y, YP,
     1             R, DEP, IDAS, RDAS, IPAR, RPAR, TMASS, S(1),
     2             S(IOFF+NATJ*JJ+1), WT, ISKWRK, RSKWRK, CSKWRK, KNAM,
     3             PNAM, NNSURF, KKBULK, SDOT, SITDOT, SDEN, DEN, SDOTI,
     4             SITDTI, ROP, KKTOT, IISUR, KKSURF, ACT, NFBULK,
     5             NLBULK, NFSURF, NLSURF, KKPHAS, KFIRST, KLAST, LMOLF)
C
      CALL CRSAVE (KKGAS, KKSURF, KKBULK, JJ, TIME, NATJ, S, SP, X,
     1             DEP, IDAS, RDAS, LSAVE, LDAS, ICRD, MORD, ACT)
C
C     RETURN
C
C     FORMAT statements
 9000 FORMAT (/10X, 'A LOT OF WORK IN DASSL, IDID=', I6)
 9010 FORMAT (///10X, 'INTEGRATION STOPPED, IDID=', I6)
 9015 FORMAT (1H1)
 9020 FORMAT (10X,  'DISTANCE=',  8X, E10.3, 7X, 'LAST STEP SIZE=',
     1         7X, E10.3, 7X, 'NUMBER OF STEPS=', 5X, I6)
 9030 FORMAT (10X, 'ORDER OF INTEGRATION=', I6,
     1         7X,  'NUMBER OF FUNCTION CALLS =', I6,
     2         7X,  'JACOBIAN EVALUATIONS=', I6,
     3        /10X, 'ERROR TEST FAILURES =', I6,
     4         7X,  'CONVERGENCE TEST FAILURES=', I6)
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRSAVE (KKGAS, KKSURF, KKBULK, JJ, TIME, NATJ,
     1                   S, SP, X, DEP, IDAS, RDAS, LSAVE, LDAS,
     2                   ICRD, MORD, ACT)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION S(*), SP(*), X(*), DEP(*), IDAS(*), RDAS(*), ACT(*)
C
      IF (ICRD .EQ. 0) THEN
C        planar coordinates
         NEQ = NATJ * JJ + 2*KKSURF
      ELSE
C        radial coordinates (or planar with a symmetry axis)
         NEQ = NATJ * JJ + KKSURF
      ENDIF
      ML = MAX0(2*NATJ-1, NATJ+KKSURF-1)
C
      IFLAG = 1
C
      WRITE (LSAVE) IFLAG
      WRITE (LSAVE) TIME, (S(N), N = 1, NEQ), (SP(N), N = 1, NEQ),
     1              (X(J), J = 1, JJ), (DEP(K), K = 1,2*KKGAS),
     2              (ACT(K), K=KKGAS+KKSURF+1,KKGAS+KKSURF+KKBULK)
C
      LIDWRK = 20 + NEQ
      LDWRK = 40 + 9*NEQ + (3*ML+1)*NEQ + 2*(NEQ/(2*ML+1)+1)
C
      WRITE (LDAS)(IDAS(I), I = 1, LIDWRK)
      WRITE (LDAS)(RDAS(I), I = 1, LDWRK)
      REWIND LDAS
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRSFUN (KKGAS, KKSURF, KKBULK, BUFFER, ISKWRK,
     1                   RSKWRK, ICKWRK, RCKWRK, T, P, SDOT, SITDOT,
     2                   NNSURF, FTIME, LVARTP, FZ, ZN, SDEN, ACT,
     3                   NFSURF, NLSURF, KKPHAS, KFIRST, KLAST, TSTEP1,
     4                   WT, YAV, S, JJ, NATJ, LTDIF, IMULT, X,
     5                   TMASS, ICRD, IMCWRK, RMCWRK, XAV, XMF,
     6                   VISC, COND, D, TDC, YV, GFAC, SFAC, LUPPER,
     7                   LVCOR)
C  
C  START PROLOGUE
C  END PROLOGUE
C
C     dummy variables
      INTEGER KKGAS
      INTEGER KKSURF
      INTEGER KKBULK
      INTEGER NNSURF
C*****precision > double
C Current value of the solution variables. The first KKGAS are
C      the mass fractions of the gas phase species at the top
C      surface.  The next KKSURF are the values for the site
C      fractions at the upper surface
      DOUBLE PRECISION BUFFER(*)
C*****END precision > double
C Work space for surface chemkin
      INTEGER ISKWRK(*)
C*****precision > double
      DOUBLE PRECISION RSKWRK(*)
C*****END precision > double
C Work array for chemkin
      INTEGER ICKWRK(*)
C*****precision > double
      DOUBLE PRECISION RCKWRK(*)
C Temperature at the surface
      DOUBLE PRECISION T
C Pressure of the gas in cgs units
      DOUBLE PRECISION P
C Rate of production of species from surface reactions
      DOUBLE PRECISION SDOT(*)
C Rate of production of surface phases
      DOUBLE PRECISION SITDOT(*)
C Density of surface phases
      DOUBLE PRECISION SDEN(*)
C Vector of activities
      DOUBLE PRECISION ACT(*)
C*****END precision > double
C Pointers to the first and last surface phase, respectively
      INTEGER NFSURF , NLSURF
C Index arrays for phases
      INTEGER KKPHAS(*) , KFIRST(*) , KLAST(*)
C*****precision > double
C Vector of the residuals that will be returned
      DOUBLE PRECISION FZ(*)
C Solution at previous time step
      DOUBLE PRECISION ZN(*)
C Current time step
      DOUBLE PRECISION TSTEP1
C Molecular weight of the gas phase species
      DOUBLE PRECISION WT(KKGAS)
C Vector of average mass fractions at the half node
      DOUBLE PRECISION YAV(KKGAS)
C Solution array
      DOUBLE PRECISION S(NATJ,JJ)
C*****END precision > double
C Number of nodes in the full problem
      INTEGER JJ
C Total number of nodes plus 6
      INTEGER NATJ
C Flag indicating whether we are considering the upper or the lower wall
      LOGICAL LUPPER
C Flag to indicate if thermal diffusion is turned on
      LOGICAL LTDIF
C Flag to indicate if a correction velocity is used
      LOGICAL LVCOR
C Flag to indicate that the residual for the transient problem
C          should be evaluated
      LOGICAL FTIME
C Flag to determine whether or not to evaluate transport coefficients
      LOGICAL LVARTP
C Flag to indicate whether multicomponent diffusion formulation is used
      INTEGER IMULT
C The coordinate system index
      INTEGER ICRD
C Work space for transport package
      INTEGER IMCWRK(*)
C*****precision > double
      DOUBLE PRECISION RMCWRK(*)
C Vector of mesh point locations
      DOUBLE PRECISION X(JJ)
C The initial total mass flow rate
      DOUBLE PRECISION TMASS
C Work space for subroutines
      DOUBLE PRECISION XAV(KKGAS) , XMF(KKGAS)
C Vector of viscosities at midponts
      DOUBLE PRECISION VISC(JJ)
C Vector of conductivities
      DOUBLE PRECISION COND(JJ)
C Array of diffusion coefficients
      DOUBLE PRECISION D(KKGAS,KKGAS,JJ)
C Array of thermal diffusion coefficients
      DOUBLE PRECISION TDC(KKGAS,JJ)
C Array of diffusion velocities
      DOUBLE PRECISION YV(KKGAS,JJ)
C
      DOUBLE PRECISION GFAC, SFAC
C*****END precision > double
C*****precision > single
C      REAL BUFFER(*), RCKWRK(*), RSKWRK(*), T, P, SDOT(*),
C     1     SITDOT(*), SDEN(*), ACT(*), GFAC, SFAC
C      REAL RMCWRK(*),XAV(KKGAS), XMF(KKGAS),VISC(JJ),COND(JJ),
C     1     D(KKGAS,KKGAS,JJ), TDC(KKGAS,JJ), YV(KKGAS,JJ), X(JJ),
C     2     TMASS, FZ(*), ZN(*), TSTEP1, WT(KKGAS), YAV(KKGAS),
C     3     S(NATJ,JJ)
C*****END precision > single
C
C     local variables
C
C*****precision > double
      DOUBLE PRECISION SUMZ
C Densities of the gas at the cell boundary and the wall, respectivley
      DOUBLE PRECISION RHOAP
C      DOUBLE PRECISION RHO
C Convective velocity at the surface (actually at JJ-1/2)
      DOUBLE PRECISION VCON
C Pressure and temperature at JJ-1/2
      DOUBLE PRECISION PAV , TAV
      DOUBLE PRECISION DZDT
C Constants
      DOUBLE PRECISION ZERO, ONE
C
      DOUBLE PRECISION VTRACE
C*****END precision > double
C*****precision > single
C      REAL SUMZ, RHOAP, VCON, PAV, TAV, DZDT, ZERO, ONE, VTRACE
C*****END precision > single
      SAVE ZERO, ONE
C Index for species
      INTEGER K
C Loop over phases
      INTEGER N
C
C     LOCS COMMON BLOCK
      INTEGER       NP, NR, NU, NT, NML, NMU, NYS, NY1
      COMMON /LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
C
C     externals
      EXTERNAL CRTRNP , CRDIFV
      EXTERNAL CKYTX , CKRHOY
      EXTERNAL SKSDEN , SKRAT
      EXTERNAL VTRACE
C
      DATA ZERO/0.0/ , ONE/1.0/
C
C-----------------------------------------------------------------------
C     store current solution in array, S
      IF (LUPPER) THEN
         JEND = JJ
      ELSE
         JEND = 1
      ENDIF
      DO 10 K = 1, KKGAS
         S(NYS + K, JEND) = BUFFER(K)
   10 CONTINUE
C
C     chemical production rates on the surface
      CALL CKYTX (BUFFER , ICKWRK , RCKWRK , ACT)
      IF (NFSURF .GT. 0) THEN
        DO 105 K = KFIRST(NFSURF), KLAST(NLSURF)
           ACT(K) = BUFFER(K)
  105   CONTINUE
      ENDIF
C
      CALL SKSDEN (ISKWRK, RSKWRK, SDEN)
      CALL SKRAT (P, T, ACT, SDEN, ISKWRK, RSKWRK, SDOT, SITDOT)
C
      DO 110 K = 1, KKGAS+KKSURF
         SDOT(K) = SFAC*SDOT(K)
  110 CONTINUE
C
      IF (LUPPER) THEN
         JJM1 = JJ - 1
      ELSE
         JJM1 = 1
      ENDIF
C
C     evaluate and store the transport coefficients
      IF(LVARTP .AND. .NOT.FTIME)
     1       CALL CRTRNP (KKGAS, 2, NATJ, LTDIF, IMULT, X(JJM1),
     2                    S(1,JJM1), YAV, XAV, ICKWRK, RCKWRK,
     3                    IMCWRK, RMCWRK, VISC(JJM1), COND(JJM1),
     4                    D(1,1,JJM1), TDC(1,JJM1))
C
C     evaluate and store the diffusion velocities
      CALL CRDIFV (KKGAS, 2, NATJ, LTDIF, IMULT, X(JJM1), S(1,JJM1),
     1             WT, TMASS, ICRD, YAV, XMF, XAV, D(1,1,JJM1),
     2             TDC(1,JJM1), ICKWRK, RCKWRK, YV(1,JJM1), LVCOR)
C
C     calculate some averages at the half-node position
      PAV = 0.5 * (S(NP,JJM1+1) + S(NP,JJM1))
      TAV = 0.5 * (S(NT,JJM1+1) + S(NT,JJM1))
      CALL CKAVG (KKGAS, S(NYS+1,JJM1+1), S(NYS+1,JJM1), YAV)
      CALL CKRHOY (PAV, TAV, YAV, ICKWRK, RCKWRK, RHOAP)
C
C     calculate convective velocity, 
C     boundary conditions on chemical species
C
      IF (LUPPER) THEN
         VCON = 0.0
         DO 2100 K = 1, KKGAS
            VCON = VCON - SDOT(K) * WT(K) / RHOAP
 2100    CONTINUE
         JJM1 = JJ - 1
         DO 2300 K = 1, KKGAS
            FZ(K) = -RHOAP * (BUFFER(K)*VCON + YV(K,JJM1)) -
     1               SDOT(K) * WT(K)
 2300    CONTINUE
      ELSE
C
         VCON = 0.0
         DO 2101 K = 1, KKGAS
            VCON = VCON + SDOT(K) * WT(K) / RHOAP
 2101    CONTINUE
         DO 2301 K = 1, KKGAS
            FZ(K) = RHOAP * (BUFFER(K)*VCON + YV(K,1)) -
     1              SDOT(K) * WT(K)
 2301    CONTINUE
      ENDIF
C
      IF (.NOT. LVCOR) THEN
C        force the sum of the mass fractions to add to one.
         FZ(KKGAS) = - VTRACE(BUFFER, KKGAS)
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF (KKSURF .GT. 0) THEN
         DO 55 N = NFSURF, NLSURF
            SUMZ = ZERO
            IF (KKPHAS(N) .GT. 1) THEN
               DO 50 K = KFIRST(N), KLAST(N) - 1
                  FZ(K) = SDOT(K) / SDEN(N)
                  SUMZ = SUMZ + BUFFER(K)
   50          CONTINUE
            ENDIF
            FZ(K) = SUMZ + BUFFER(K) - ONE
   55    CONTINUE
      ENDIF
C
      IF (.NOT.FTIME .OR. KKSURF.LE.1) RETURN
C
C     add the time step, if needed
      DO 2055 N = NFSURF, NLSURF
         IF (KKPHAS(N) .GT. 1) THEN
            DO 2050 K = KFIRST(N), KLAST(N) - 1
               DZDT = (BUFFER(K) - ZN(K)) / TSTEP1
               FZ(K) = DZDT - FZ(K)
 2050       CONTINUE
         ENDIF
 2055 CONTINUE
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRSJAC (KKGAS, KKSURF, KKBULK, COMPS, LTIME, ABSOL,
     1                   RELAT, TSURF, P, ISKWRK, RSKWRK, IPAR, RPAR,
     2                   ZN, Z, FZ, FZN, A, Y, SDOT, SITDOT, NNSURF,
     3                   SDEN, ACT, NFSURF, NLSURF, KKPHAS, KFIRST,
     4                   KLAST, TSTEP1, WT, YAV, S, JJ, NATJ, LTDIF,
     5                   IMULT, X, TMASS, ICRD, IMCWRK, RMCWRK, XAV,
     6                   XMF, VISC, COND, D, TDC, YV, GFAC, SFAC,
     7                   LUPPER, LVCOR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INTEGER COMPS
C*****precision > double
C Molecular weight of the gas phase species
      DOUBLE PRECISION WT(KKGAS)
C Vector of averages mass fractions at JJ-1/2
      DOUBLE PRECISION YAV(KKGAS)
C Solution array
      DOUBLE PRECISION S(NATJ,JJ)
C*****END precision > double
C Number of nodes in the full problem
      INTEGER JJ
C Total number of nodes plus 6
      INTEGER NATJ
C Flag to indicate if thermal diffusion is turned on
      LOGICAL LTDIF
C Flag to indicate if a correction velocity is used
      LOGICAL LVCOR
C Flag to indicate whether multicomponent diffusion formulation is used
      INTEGER IMULT
C Flag to determine whether to evaluate transport coefficients
      LOGICAL LVARTP
C The coordinate system index
      INTEGER ICRD
C Work space for transport package
      INTEGER IMCWRK(*)
C*****precision > double
      DOUBLE PRECISION RMCWRK(*)
C Vector of mesh point locations
      DOUBLE PRECISION X(JJ)
C The initial total mass flow rate
      DOUBLE PRECISION TMASS
C Work space for subroutines
      DOUBLE PRECISION XAV(KKGAS) , XMF(KKGAS)
C Vector of viscosities at midponts
      DOUBLE PRECISION VISC(JJ)
C Vector of conductivities
      DOUBLE PRECISION COND(JJ)
C Array of diffusion coefficients
      DOUBLE PRECISION D(KKGAS,KKGAS,JJ)
C Array of thermal diffusion coefficients
      DOUBLE PRECISION TDC(KKGAS,JJ)
C Array of diffusion velocities
      DOUBLE PRECISION YV(KKGAS,JJ)
C*****END precision > double
C*****precision > single
C      REAL RMCWRK(*),XAV(KKGAS),XMF(KKGAS),VISC(JJ),COND(JJ),
C     1     D(KKGAS,KKGAS,JJ),TDC(KKGAS,JJ),YV(KKGAS,JJ), X(JJ),
C     2     TMASS, WT(KKGAS),YAV(KKGAS), S(NATJ,JJ)
C*****END precision > single
C
      DIMENSION Z(*), FZ(*), FZN(*), A(COMPS, *), ZN(*), SDOT(*),
     1          SITDOT(*), SDEN(*), Y(*), ISKWRK(*), RSKWRK(*), IPAR(*),
     2          RPAR(*), ACT(*), KKPHAS(*), KFIRST(*), KLAST(*)
      LOGICAL LTIME, LUPPER
C
C     zero the matrix storage space.
C*****precision > double
      CALL DCOPY (COMPS * COMPS, 0.0D0, 0, A, 1)
C*****END precision > double
C
C*****precision > single
C      CALL SCOPY (COMPS * COMPS, 0.0, 0, A, 1)
C*****END precision > single
C
C     call the function at S and store in FZN.
      LVARTP = .TRUE.
      CALL CRSFUN (KKGAS, KKSURF, KKBULK, Z, ISKWRK, RSKWRK, IPAR,
     1             RPAR, TSURF, P, SDOT, SITDOT, NNSURF, LTIME, LVARTP,
     2             FZN, ZN, SDEN, ACT, NFSURF, NLSURF, KKPHAS,
     3             KFIRST, KLAST, TSTEP1, WT, YAV, S, JJ, NATJ,
     4             LTDIF, IMULT, X, TMASS, ICRD, IMCWRK, RMCWRK, XAV,
     5             XMF, VISC, COND, D, TDC, YV, GFAC, SFAC, LUPPER,
     6             LVCOR)
      LVARTP = .FALSE.
C
C     top of the loops over the residue classes and solution
C     components
      DO 0200 M = 1, COMPS
C        for a given residue class and a given solution componnt,
C        perturb the S vector at points in the same residue class.
C
         SAVE = Z(M)
         PERTRB = ABS(Z(M)) * RELAT + ABSOL
         Z(M) = Z(M) + PERTRB
C
C        call the function at the perturbed S and store the result
C        in FZ.
         CALL CRSFUN (KKGAS, KKSURF, KKBULK, Z, ISKWRK, RSKWRK, IPAR,
     1                RPAR, TSURF, P, SDOT, SITDOT, NNSURF, LTIME,
     2                LVARTP, FZ, ZN, SDEN, ACT, NFSURF, NLSURF,
     3                KKPHAS, KFIRST, KLAST, TSTEP1, WT, YAV, S, JJ,
     4                NATJ, LTDIF, IMULT, X, TMASS, ICRD, IMCWRK,
     5                RMCWRK, XAV, XMF, VISC, COND, D, TDC, YV, GFAC,
     6                SFAC, LUPPER, LVCOR)
C
C        restore S to its original value.
         Z(M) = SAVE
C
C        difference to get the columns of the Jacobian
         DO 0100 N = 1, COMPS
            A(N, M) = (FZ(N) - FZN(N)) / PERTRB
  100    CONTINUE
C
C     bottom of the loops over the residue classes and solution
C     components
  200 CONTINUE
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRSPRT (KNAM, PNAM, ICKWRK, RCKWRK, TEXT, BUFFER,
     1                   ISKWRK, RSKWRK, P, TSURF, KKGAS, KKSURF,NNSURF,
     2                   NNBULK, NFSURF, NLSURF, NFBULK, NLBULK, KKPHAS,
     3                   KFIRST, KLAST, SDOT, SITDOT, SDEN, ACT, WT,
     4                   DEN, Y, PATM, IISUR, ROP, CSKWRK, YAV, S, JJ,
     5                   NATJ, LTDIF, IMULT, X, TMASS, ICRD, IMCWRK,
     6                   RMCWRK, XAV, XMF, VISC, COND, D, TDC, YV, FZ,
     7                   SFAC, LUPPER, LVCOR)
C
C  START PROLOGUE
C  END PROLOGUE
C
      INTEGER TEXT, KKGAS, KKSURF, NNSURF, NNBULK, NFSURF, NLSURF,
     1        NFBULK, NLBULK, IISUR
C*****precision > double
C Pressure of the gas in cgs units
      DOUBLE PRECISION P
C Value of 1 atm in the current units
      DOUBLE PRECISION PATM
C Temperature of the surface
      DOUBLE PRECISION TSURF
C Molecular weight of the gas phase species
      DOUBLE PRECISION WT(*)
C Vector of averages mass fractions at JJ-1/2
      DOUBLE PRECISION YAV(KKGAS)
C Solution array
      DOUBLE PRECISION S(NATJ,JJ)
C*****END precision > double
C Number of nodes in the full problem
      INTEGER JJ
C Total number of nodes plus 6
      INTEGER NATJ
C Flag to indicate if thermal diffusion is turned on
      LOGICAL LTDIF
C Flag to indicate if a correction velocity is used
      LOGICAL LVCOR
C Flag to indicate whether we are solving for the upper or lower wall
      LOGICAL LUPPER
C Flag to indicate whether multicomponent diffusion formulation is used
      INTEGER IMULT
C The coordinate system index
      INTEGER ICRD
C Work space for transport package
      INTEGER IMCWRK(*)
C*****precision > double
      DOUBLE PRECISION RMCWRK(*)
C Vector of mesh point locations
      DOUBLE PRECISION X(JJ)
C The initial total mass flow rate
      DOUBLE PRECISION TMASS
C Work space for subroutines
      DOUBLE PRECISION XAV(KKGAS) , XMF(KKGAS)
C Vector of viscosities at midponts
      DOUBLE PRECISION VISC(JJ)
C Vector of conductivities
      DOUBLE PRECISION COND(JJ)
C Array of diffusion coefficients
      DOUBLE PRECISION D(KKGAS,KKGAS,JJ)
C Array of thermal diffusion coefficients
      DOUBLE PRECISION TDC(KKGAS,JJ)
C Array of diffusion velocities
      DOUBLE PRECISION YV(KKGAS,JJ)
C Current value of the gas phase residuals
      DOUBLE PRECISION FZ(*)
C
      DOUBLE PRECISION SFAC
C*****END precision > double
C*****precision > single
C      REAL RMCWRK(*), X(JJ), TMASS, XAV(KKGAS), XMF(KKGAS),
C     1     VISC(JJ), COND(JJ), D(KKGAS,KKGAS,JJ), TDC(KKGAS,JJ),
C     2     YV(KKGAS,JJ), WT(*), YAV(KKGAS), S(NATJ,JJ),
C     3     FZ(*), P, TSURF, PATM, SFAC
C*****END precision > single

      CHARACTER*(*) CSKWRK(*), KNAM(*), PNAM(*)
      CHARACTER*48 RNAM
      CHARACTER*80 ISTR
      LOGICAL IERR
C*****precision > double
      DOUBLE PRECISION BUFFER(*), Y(*), SDOT(*), SITDOT(*),
     1                 SDEN(*), DEN(*), ROP(*)
      DOUBLE PRECISION RSKWRK(*), ACT(*), RCKWRK(*)
C*****END precision > double
C*****precision > single
C      REAL BUFFER(*), Y(*), SDOT(*), SITDOT(*),
C     1                 SDEN(*), DEN(*), ROP(*)
C      REAL RSKWRK(*), ACT(*), RCKWRK(*)
C*****END precision > single
      INTEGER ISKWRK(*), KKPHAS(*), KFIRST(*), KLAST(*), ICKWRK(*)
C
C     LOCS COMMON BLOCK
      INTEGER       NP, NR, NU, NT, NML, NMU, NYS, NY1
      COMMON /LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
C
C     local variables
      INTEGER K , N, I, ILS, LT
C*****precision > double
      DOUBLE PRECISION PAV, TAV, RHOAP, VCON, DRATE
C*****END precision > double
C*****precision > single
C      REAL PAV, TAV, RHOAP, VCON, DRATE
C*****END precision > single

C     externals
      INTEGER  CKLSCH
      EXTERNAL CKLSCH
C
      DATA LT/48/
C
C     Store current gas species solution in array, S
      IF ( LUPPER) THEN
         JWALL = JJ
         JGAS  = JJ-1
      ELSE
         JWALL = 1
         JGAS  = 2
      ENDIF
C
      DO 10 K = 1 , KKGAS
         S(NYS + K, JWALL) = BUFFER(K)
   10 CONTINUE

      CALL CKYTX (BUFFER, ICKWRK, RCKWRK, ACT)
      IF (NFSURF .GT. 0) THEN
        DO 105 K = KFIRST(NFSURF), KLAST(NLSURF)
           ACT(K) = BUFFER(K)
  105   CONTINUE
      ENDIF
C
      CALL SKSDEN (ISKWRK, RSKWRK, SDEN)
      CALL SKRAT  (P, TSURF, ACT, SDEN, ISKWRK, RSKWRK, SDOT, SITDOT)
C
      DO 110 K = 1, KKGAS+KKSURF
         SDOT(K) = SFAC*SDOT(K)
  110 CONTINUE
C
C     Evaluate and store the transport coefficients
      CALL CRTRNP (KKGAS, JJ, NATJ, LTDIF, IMULT, X, S, YAV,
     1                         XAV, ICKWRK, RCKWRK, IMCWRK, RMCWRK,
     2                         VISC, COND, D, TDC)
C
C     Evaluate and store the diffusion velocities
      CALL CRDIFV (KKGAS, JJ, NATJ, LTDIF, IMULT, X, S, WT, TMASS, ICRD,
     1             YAV, XMF, XAV, D, TDC, ICKWRK, RCKWRK, YV, LVCOR)
C
C-----UPPER BOUNDARY NODE
C
C     Calculate some averages at the half-node position
      PAV = 0.5 * (S(NP,JWALL) + S(NP,JGAS))
      TAV = 0.5 * (S(NT,JWALL) + S(NT,JGAS))
      CALL CKAVG (KKGAS, S(NYS+1,JWALL), S(NYS+1,JGAS), YAV)
      CALL CKRHOY (PAV, TAV, YAV, ICKWRK, RCKWRK, RHOAP)
C
      WRITE (TEXT,'(/1X,A/)') 'GAS PHASE SPECIES MASS FRACTIONS'
      WRITE(TEXT,290)
C
C     calculate convective velocity, boundary conditions on
C     chemical species
C
      IF (LUPPER) THEN
         VCON = 0.0
         DO 2100 K = 1, KKGAS
            VCON = VCON - SDOT(K) * WT(K) / RHOAP
 2100    CONTINUE
         JJM1 = JJ - 1
         DO 300 K = 1, KKGAS
            N = NYS + K
            FZ(K) = -RHOAP * (BUFFER(K)*VCON + YV(K,JJM1)) -
     1               SDOT(K) * WT(K)
            WRITE(TEXT,295) KNAM(K), S(N,JWALL) , S(N,JGAS),
     1                      SDOT(K) , FZ(K), -RHOAP*BUFFER(K)*VCON,
     2                      -RHOAP*YV(K,JJM1), SDOT(K)*WT(K)
  300    CONTINUE
      ELSE
C
         VCON = 0.0
         DO 2101 K = 1, KKGAS
            VCON = VCON + SDOT(K) * WT(K) / RHOAP
 2101    CONTINUE
         DO 301 K = 1, KKGAS
            N = NYS + K
            FZ(K) = RHOAP * (BUFFER(K)*VCON + YV(K,1)) -
     1              SDOT(K) * WT(K)
            WRITE(TEXT,295) KNAM(K), S(N,JWALL) , S(N,JGAS),
     1                      SDOT(K) , FZ(K), RHOAP*BUFFER(K)*VCON,
     2                      RHOAP*YV(K,1), SDOT(K)*WT(K)
  301    CONTINUE
      ENDIF
C
      WRITE ( TEXT, 310) VCON
      WRITE ( TEXT, 320) S(NT,JWALL) 
      WRITE ( TEXT, 330) S(NT,JGAS)
C
      IF (NNSURF .GT. 0) THEN
C
      WRITE (TEXT,'(/1X,A/)') 'SURFACE SPECIES SITE FRACTIONS'
C        surface species
         DO 500 N = NFSURF, NLSURF
            ILS = CKLSCH(PNAM(N))
            ISTR = ' '
            ISTR = PNAM(N)(1:ILS)
            ILS = N-1
            WRITE (TEXT, 7050) ILS, ISTR
            IF (KKPHAS(N) .LE. 0) THEN
               WRITE (TEXT, *) ' -----NO SPECIES ON THIS SITE-----'
            ELSE
               DO 450 K = KFIRST(N), KLAST(N)
                  WRITE (TEXT, 7060) KNAM(K), BUFFER(K), SDOT(K)
  450          CONTINUE
            ENDIF
500      CONTINUE
      ENDIF
C
      IF (NNBULK .GT. 0) THEN
C
C        print total deposition rates for all deposit species
         WRITE (TEXT,'(/1X,A,/)') 'DEPOSITION RATE (MICRONS/MIN)'
C
         CALL SKDEN (P, TSURF, ACT, SDEN, ISKWRK, RSKWRK, DEN)
         DO 2000 K = KFIRST(NFBULK), KLAST(NLBULK)
            ILS = CKLSCH(KNAM(K))
            ISTR = ' '
            ISTR = KNAM(K)(1:ILS)
            DRATE = SDOT(K)*WT(K)/DEN(K) * 60.0E4
            ILS = K-KKGAS
            WRITE (TEXT, 7100) ILS, ISTR, DRATE
 2000    CONTINUE
      ENDIF
C
      IF (IISUR .LE. 0) RETURN
      WRITE (TEXT,'(/1X,A)') 'RATE OF PROGRESS OF REACTIONS'
C
      DO 2200 I = 1, IISUR
         CALL SKROP  (P, TSURF, ACT, SDEN, ISKWRK, RSKWRK, ROP)
         CALL SKSYMR (I, TEXT, ISKWRK, RSKWRK, CSKWRK, LT, RNAM, IERR)
         WRITE (TEXT, 7200) I, RNAM, ROP(I)
 2200 CONTINUE
C
C     RETURN
C
290   FORMAT (1h ,'Species_Name    Mass_frac(Wall) Mass_frac(Gas)    '
     1       ,'   Wdot       Residual  Bulk_Vel_Contrb'
     1       ,'    Dif_Vel_Cont  Surface_Mass_Rate'
     1 /120('='))
295   FORMAT (1H , A16,3X,G11.5,4X,G11.5,4X,G11.5,4X,G11.5,4X,G11.5,4X
     1        ,G11.5,4X,G11.5)
310   FORMAT (120('=')//10X,'Total Velocity into wall'
     1       ,' (positive indicates a net flux into wall) = '
     1        ,G11.5,' cm/sec'/)
320   FORMAT ( ' Temperature at the wall = ', F10.3)
330   FORMAT ( ' Temperature in the gas  = ', F10.3)
7050  FORMAT (' SITE ',I2,': ', A16, 20X, ' SITDOT')
7060  FORMAT (4X,  A12, 2(1PE11.3, 14X))
7100  FORMAT (I4, 2X, A16, 1PE11.3)
7150  FORMAT (4X,  A12, 2(1PE11.3, 4X))
7200  FORMAT (I4, 2X, A48, 5(1PE11.3))
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRTIM (KKGAS, JJ, NATJ, LSAVE, LDAS, LOUT, IDAS, RDAS,
     1                  KKSURF, ICRD)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION IDAS(*), RDAS(*)
C
C     get statistics from DASSL
      TIME   = RDAS(4)
      HLAST  = RDAS(7)
      NORDER = IDAS(8)
      NSTEP  = IDAS(11)
      NRES   = IDAS(12)
      NNJAC  = IDAS(13)
      NERR   = IDAS(14)
      NCON   = IDAS(15)
C
      WRITE (LOUT, 9010)
      WRITE (LOUT, 9020)  TIME, HLAST, NSTEP
      WRITE (LOUT, 9030)  NORDER, NRES, NNJAC, NERR, NCON
C
      IF (ICRD .EQ. 0) THEN
C        planar coordinates
         NEQ = NATJ * JJ + 2*KKSURF
      ELSE
C        radial coordinates (or planar with a symmetry axis)
         NEQ = NATJ * JJ + KKSURF
      ENDIF
C
      ML = MAX0(2*NATJ-1, NATJ+KKSURF-1)
      LIDWRK = 20 + NEQ
      LDWRK  = 40 + 9*NEQ + (3*ML+1)*NEQ + 2*(NEQ/(2*ML+1)+1)
C
C     read DASSL work space at the last printout and save it
      READ (LDAS) (IDAS(I), I = 1, LIDWRK)
      READ (LDAS) (RDAS(I), I = 1, LDWRK)
      IFLAG = -1
C
      WRITE (LSAVE) IFLAG
      WRITE (LSAVE) (IDAS(I), I = 1,LIDWRK)
      WRITE (LSAVE) (RDAS(I), I = 1,LDWRK)
      WRITE (LSAVE)
C
      RETURN
C
C     FORMAT statements
 9000 FORMAT (' QUITTING EXECUTION BECAUSE OF TIME LIMIT')
 9010 FORMAT (1H1)
 9020 FORMAT (10X, 'DISTANCE=', 8X, E10.3, 7X, 'LAST STEP SIZE=',
     1         7X, E10.3, 7X, 'NUMBER OF STEPS=', 5X, I6)
 9030 FORMAT (10X, 'ORDER OF INTEGRATION=', I6,
     1         7X,  'NUMBER OF FUNCTION CALLS =', I6,
     2         7X,  'JACOBIAN EVALUATIONS=', I6,
     3        /10X, 'ERROR TEST FAILURES =', I6,
     4         7X,  'CONVERGENCE TEST FAILURES=', I6)
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRTMPR (GTMP, STMP, XSPSTR, XSPEND, TGTMP, IDOSP, XTMP)
C
C  START PROLOGUE
C   This subroutine sets the spline temperature points (SPLX, SPLY) on
C the susceptor.  The splines are later used by the subroutine, CRWALL,
C to calculate the temperature at a given axial position.
C
C  a) The variable IDOSP is first checked to see if the user input
C     spline temperature values.  If this is the case, then these
C     are honored.
C  b) If not then the values of STMP, the susceptor temperature,
C     GTMP, the starting temperature, and XTMP, the distance over
C     which the temperature ramp takes place are utilized to
C     create a spline fit to a gradual temperature rize on the
C     susceptor.  At X=0.0, the susceptor starts at a temperature
C     of GTMP.  At X>=XTMP, the susceptor attains a temperature of
C     STMP.
C  c) If XTMP is set to be equal or less than zero, then the
C     susceptor is set to be equal to STMP at X=0.0, i.e. no
C     ramping of the temperature profile is done.  This may
C     have repercussions on the problem startup.
C
C   This subroutine does not deal with the upper wall's temperature
C for planar coordinates.
C
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      REAL SPLX(50), SPLY(50), SPLYP(50), SPLYPP(50), SPLW(150), TMPSTR,
     1     TMPEND, A1, AN, B1, BN, DELTT, SXTMP
C
      COMMON /SPLN/    SPLX, SPLY, SPLYP, SPLYPP,
     1                 SPLW, ISPLN, TMPSTR, TMPEND
C
      IF (IDOSP .EQ. 0) THEN
C*****precision > quad
C         TMPEND = SNGLQ(STMP)
C         TMPSTR = SNGLQ(GTMP)
C         SXTMP  = SNGLQ(XTMP)
C*****END precision > quad
C*****precision > double
         TMPEND = SNGL(STMP)
         TMPSTR = SNGL(GTMP)
         SXTMP  = SNGL(XTMP)
C*****END precision > double
C*****precision > single
C         TMPEND = STMP
C         TMPSTR = GTMP
C         SXTMP  = XTMP
C*****END precision > single
         IF (SXTMP .LE. 0.0) THEN
           TMPSTR = TMPEND
           SXTMP  = 0.5
         ENDIF
         DELTT = TMPEND - TMPSTR
         SPLX(1) = 0.
         SPLY(1) = TMPSTR
         SPLX(2) = SPLX(1) + 0.000001
         SPLY(2) = TMPSTR
         SPLX(3) = SPLX(1) + 0.5*SXTMP
         SPLY(3) = TMPSTR + 0.5*DELTT
         SPLX(4) = SPLX(1) + SXTMP
         SPLY(4) = TMPEND
         SPLX(5) = SPLX(4) + 0.000001
         SPLY(5) = TMPEND
         ISPLN = 5
      ENDIF
C
      XSPSTR = SPLX(1)
      XSPEND = SPLX(ISPLN)
      TGTMP = GTMP
C
C     set parameters for spline routine
      ISX = 0
      A1 = -0.5
      AN = -0.5
      B1 = 3.0*((SPLY(2)-SPLY(1))/(SPLX(2)-SPLX(1)))/(SPLX(2)-SPLX(1))
      BN = 3.0*((SPLY(ISPLN)-SPLY(ISPLN-1)) /
     1        (SPLX(ISPLN)-SPLX(ISPLN-1)))/(SPLX(ISPLN)-SPLX(ISPLN-1))
C
C     call spline routine
      CALL SPLIFT (SPLX, SPLY, SPLYP, SPLYPP, ISPLN, SPLW, ISPERR,
     1             ISX, A1, B1, AN, BN)
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRTRNP (KKGAS, JJ, NATJ, LTDIF, IMULT, X, S, YAV, XAV,
     1                   IPAR, RPAR, IMCWRK, RMCWRK, VISC, COND, D,
     2                   TDC)
C
C  START PROLOGUE
C
C   KKGAS  - NUMBER OF CHEMICAL SPECIES.
C   JJ     - NUMBER OF MESH POINTS.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES AT EACH MESH POINT.  ALSO,
C            NATJ MUST BE THE EXACT FIRST DIMENSION OF S(NATJ,*).
C   LTDIF  - IF LTDIF=.TRUE. THEN EVALUATE THERMAL DIFFUSION RATIOS AS
C            WELL AS DIFFUSION COEFFICIENTS.
C   IMULT  - IF IMULT = 1, THEN USE MULTICOMPONENT TRANSPORT COEFFICIENT
C            OTHERWISE, USE SIMPLER MIXTURE COEFFICIENTS.
C   P      - PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   X      - THE ARRAY OF MESH POINT LOCATIONS.
C              DIMENSION X(*) AT LEAST JJ
C              CGS UNITS - CM
C   S      - DEPENDENT VARIABLE MATRIX.  THE TEMPERATURES ARE STORED IN
C            T(J)=S(NT,J), THE MASS FRACTIONS ARE IN Y(K,J)=S(NYS+K,J),
C            AND THE FLOW RATES ARE IN FLRT(J)=S(NM,J).
C              DIMENSION S(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION,
C              AND AT LEAST JJ FOR THE SECOND.
C
C WORK AND SCRATCH SPACE
C   YAV    - ARRAY OF MASS FRACTIONS AT MESH MIDPOINTS.
C              DIMENSION Y(*) AT LEAST KKGAS.
C   XAV    - ARRAY OF MOLE FRACTIONS AT MESH MIDPOINTS.
C              DIMENSION X(*) AT LEAST KKGAS.
C   IPAR   - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   RPAR   - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   IMCWRK - INTEGER TRANSPORT PROPERTY WORK SPACE.
C              DIMENSIONING - SEE TRANSPORT DOCUMENTATION.
C   RMCWRK - FLOATING POINT TRANSPORT PROPERTY WORK SPACE.
C
C              DIMENSIONING - SEE TRANSPORT DOCUMENTATION.
C OUTPUT-
C   COND   - ARRAY OF CONDUCTIVITIES AT THE MESH MID-POINTS.
C   VISC   - ARRAY OF VISCOSITIES AT THE MESH MID-POINTS.
C   D      - MATRIX OF MULTICOMPONENT DIFFUSION COEFFICIENTS AT THE MESH
C            MID-POINTS.
C              DIMENSION D(KKGAS,KKGAS,*) EXACTLY KKGAS FOR THE FIRST
C              TWO DIMENSIONS, AND AT LEAST JJ FOR THE THIRD.
C   TDC    - MATRIX OF THERMAL DIFFUSION COEFFICIENTS AT THE MESH
C              MID-POINTS.
C              DIMENSION TDC(KKGAS,*) EXACTLY KKGAS FOR THE FIRST
C              DIMENSION, AND AT LEAST JJ FOR THE SECOND.
C
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
      DIMENSION S(NATJ,*), YAV(*), XAV(*), IPAR(*), RPAR(*), IMCWRK(*),
     1          RMCWRK(*), VISC(*), COND(*), D(KKGAS,KKGAS,*),
     2          TDC(KKGAS,*)
      LOGICAL LTDIF
C
      IF (IMULT .EQ. 1) THEN
         DO 200 J = 1, JJ-1
            JP1 = J + 1
            TAV = 0.5 * (S(NT,J) + S(NT,JP1))
            PAV = 0.5 * (S(NP,J) + S(NP,JP1))
            CALL CKYTX (S(NY1,J),   IPAR, RPAR, YAV)
            CALL CKYTX (S(NY1,JP1), IPAR, RPAR, XAV)
            DO 100 K = 1, KKGAS
               N = NYS + K
               XAV(K) = 0.5 * (YAV(K) + XAV(K))
               YAV(K) = 0.5 * (S(N,J) + S(N,JP1))
  100       CONTINUE
            CALL MCMDIF (PAV, TAV, XAV, KKGAS, IMCWRK, RMCWRK, D(1,1,J))
            CALL MCAVIS (TAV, XAV, RMCWRK, VISC(J))
            CALL MCMCDT (PAV, TAV, XAV, IMCWRK, RMCWRK, IPAR, RPAR,
     1                   TDC(1,J), COND(J))
  200    CONTINUE
         RETURN
      ENDIF
C
      DO 300 J = 1, JJ-1
         JP1 = J + 1
         TAV = 0.5 * (S(NT,J) + S(NT,JP1))
         PAV = 0.5 * (S(NP,J) + S(NP,JP1))
         CALL CKYTX (S(NY1,J),   IPAR, RPAR, YAV)
         CALL CKYTX (S(NY1,JP1), IPAR, RPAR, XAV)
         DO 250 K = 1, KKGAS
            N = NYS + K
            XAV(K) = 0.5 * (YAV(K) + XAV(K))
            YAV(K) = 0.5 * (S(N,J) + S(N,JP1))
250      CONTINUE
         CALL MCADIF (PAV, TAV, XAV, RMCWRK, D(1,1,J))
         CALL MCAVIS (TAV, XAV, RMCWRK, VISC(J))
         IF (LTDIF) CALL MCMCDT (PAV, TAV, XAV, IMCWRK, RMCWRK,
     1                   IPAR, RPAR, TDC(1,J), COND(J))
         CALL MCACON (TAV, XAV, RMCWRK, COND(J))
300   CONTINUE
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRTWO (COMPS, TEXT, PMAX, RSIZE, RTPWRK, ABOVE, 
     1                  BELOW, BUFFER, X, A, ITPWRK, IPAR,
     2                  RPAR, ISKWRK, RSKWRK, P, Z, Y, KKGAS, KKSURF,
     3                  KKBULK, NNSURF, NNBULK, KKPHAS, KFIRST, KLAST,
     4                  SDOT, SITDOT, FZ, ZN, FZN, IPVT, WT, DEN, SDEN,
     5                  ACT, NFSURF, NLSURF, NFBULK, NLBULK, KNAM, PNAM,
     6                  PATM, IISUR, ROP, CSKWRK, YAV, S, JJ, NATJ,
     7                  LTDIF, IMULT, GRID, TMASS, ICRD, IMCWRK, RMCWRK,
     8                  XAV, XMF, VISC, COND, D, TDC, YV, AA, RTEMP, DR,
     9                  DC, XTEMP, SSABS, SSREL, TDABS, TDREL, IPRNT,
     *                  LIMIT, NINIT, IRETIR, NJAC, ITJAC, DTMAX, DTMIN,
     2                  GFAC, SFAC, TSTEP0, LUPPER, LVCOR, KERR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C Dummy Variables

      IMPLICIT INTEGER (Q)
      INTEGER   KKGAS, KKSURF, KKBULK, NNSURF, NNBULK, NFSURF, NLSURF,
     1          NFBULK, NLBULK, IISUR, IPRNT, LIMIT, COMPS, PMAX
      INTEGER   ISKWRK(*), ITPWRK(*), IPVT(*),
     1          KKPHAS(*), KFIRST(*), KLAST(*), IPAR(*)
      CHARACTER*(*) CSKWRK(*), KNAM(*), PNAM(*)
      CHARACTER REPORT*80, SIGNAL*80, VERSIO*80
      LOGICAL LDUM2, KERR
C*****precision > double
C Vector of averages mass fractions at JJ-1/2
      DOUBLE PRECISION YAV(KKGAS)
C Solution array
      DOUBLE PRECISION S(NATJ,JJ)
C*****END precision > double
C Number of nodes in the full problem
      INTEGER JJ
C Total number of nodes plus 6
      INTEGER NATJ
C Flag to indicate if thermal diffusion is turned on
      LOGICAL LTDIF
C Flag to indicate if a correction velocity is used
      LOGICAL LVCOR
C Flag to indicate whether multicomponent diffusion formulation is used
      INTEGER IMULT
C The coordinate system index
      INTEGER ICRD
C Work space for transport package
      INTEGER IMCWRK(*)
C*****precision > double
      DOUBLE PRECISION RMCWRK(*)
C Vector of mesh point locations
      DOUBLE PRECISION GRID(JJ)
C The initial total mass flow rate
      DOUBLE PRECISION TMASS
C Work space for subroutines
      DOUBLE PRECISION XAV(KKGAS) , XMF(KKGAS)
C Vector of viscosities at midponts
      DOUBLE PRECISION VISC(JJ)
C Vector of conductivities
      DOUBLE PRECISION COND(JJ)
C Array of diffusion coefficients
      DOUBLE PRECISION D(KKGAS,KKGAS,JJ)
C Array of thermal diffusion coefficients
      DOUBLE PRECISION TDC(KKGAS,JJ)
C Array of diffusion velocities
      DOUBLE PRECISION YV(KKGAS,JJ)
C*****END precision > double
C*****precision > single
C      REAL RMCWRK(*), GRID(JJ), TMASS, XAV(KKGAS), XMF(KKGAS),
C     1     VISC(JJ), COND(JJ), D(KKGAS,KKGAS,JJ), TDC(KKGAS,JJ),
C     2     YV(KKGAS,JJ), YAV(KKGAS), S(NATJ,JJ)
C*****END precision > single

C*****precision > double
      DOUBLE PRECISION
C*****END precision > double
C*****precision > single
C      REAL
C*****END precision > single
     1          X(*),  RSKWRK(*), Z(*), Y(*), SDOT(*),
     1          SITDOT(*), WT(*), DEN(*), ROP(*), FZ(*), ZN(*), FZN(*),
     2          ABOVE(*), BELOW(*), BUFFER(*), RTPWRK(*), A(COMPS, *),
     3          ACT(*), SDEN(*), PATM, P, RPAR(*), RDUMMY,
     4          GFAC, SFAC, SSABS, SSREL, TDABS, TDREL

C  XLINPK Variables

C*****precision > double
      DOUBLE PRECISION
C*****END precision > double
C*****precision > single
C      REAL
C*****END precision > single
     1       AA(COMPS,COMPS), RTEMP(COMPS), DR(COMPS), DC(COMPS),
     1       XTEMP(COMPS)
C*****xlinpk: precision > double
C      INTEGER INFO
C      DOUBLE PRECISION   ANORM, RELERR
C      EXTERNAL XDGECO, XDGESL
C*****END xlinpk: precision > double
C*****xlinpk: precision > single
C      INTEGER INFO
C      REAL     ANORM, RELERR
C      EXTERNAL XSGECO, XSGESL
C*****END xlinpk: precision > single

C     Local Variables
      INTEGER TEXT, RSIZE, K, NINIT, IRETIR, NJAC, ITJAC
      LOGICAL ERROR, TIME, XLINPK, LVARTP
      LOGICAL LUPPER
C
C*****precision > double
      DOUBLE PRECISION D1MACH
      DOUBLE PRECISION
C*****END precision > double
C*****precision > single
C      REAL R1MACH
C      REAL
C*****END precision > single
     1            TSURF, ABSOL, RELAT, SFLR,
     2            TSTEP0, TSTEP1, DTMAX, DTMIN,
     3            CONDIT, RCOND
      CHARACTER CNTRL*80
      INTEGER CNTRLS, IVALUE
      LOGICAL LVALUE, LDUM1
      DOUBLE PRECISION RVALUE
C
C     LOCS COMMON BLOCK
      INTEGER       NP, NR, NU, NT, NML, NMU, NYS, NY1
      COMMON /LOCS/ NP, NR, NU, NT, NML, NMU, NYS, NY1
C
C     pointers to TWOPNT's control values
      PARAMETER
     +  (QADAPT = 01,  QLEVD = 02,  QLEVM = 03,  QPADD = 04,
     +   QSSABS = 05, QSSAGE = 06, QSSREL = 07,  QSTDY = 08,
     +   QSTPS0 = 09, QSTPS1 = 10, QSTPS2 = 11, QSTRID = 12,
     +   QTDABS = 13, QTDAGE = 14,  QTDEC = 15, QTDREL = 16,
     +    QTINC = 17,  QTMAX = 18,  QTMIN = 19,  QTOL0 = 20,
     +    QTOL1 = 21,  QTOL2 = 22)
C
      PARAMETER (CNTRLS = 22)
C
      DIMENSION 
     +   CNTRL(CNTRLS), IVALUE(CNTRLS), LVALUE(CNTRLS), RVALUE(CNTRLS)
C
C-----------------------------------------------------------------------
C     
      XLINPK = .FALSE.
      KERR = .FALSE.
C*****xlinpk: precision > single
C      XLINPK = .TRUE.
C*****END xlinpk: precision > single
C*****xlinpk: precision > double
C      XLINPK = .TRUE.
C*****END xlinpk: precision > double
C
C
      CNTRL(QADAPT) = 'ADAPT'
      LVALUE(QADAPT) = .FALSE.
      CALL TWSETL (ERROR, LOUT, CNTRL(QADAPT), LVALUE(QADAPT))
      KERR = KERR.OR.ERROR
C
      IF (IPRNT .LT. 10) THEN
         CNTRL(QLEVD) = 'LEVELD'
         IVALUE(QLEVD) = MIN (4, IPRNT)
C
         CNTRL(QLEVM) = 'LEVELM'
         IVALUE(QLEVM) = MIN (4, IPRNT)
      ELSE
         CNTRL(QLEVM) = 'LEVELM'
         IVALUE(QLEVM) = MAX (IPRNT/10, MOD (IPRNT, 10))
C
         CNTRL(QLEVD) = 'LEVELD'
         IVALUE(QLEVD) = MOD (IPRNT, 10)
      ENDIF
      CALL TWSETI (ERROR, LOUT, CNTRL(QLEVM), IVALUE(QLEVM))
      KERR = KERR.OR.ERROR
      CALL TWSETI (ERROR, LOUT, CNTRL(QLEVD), IVALUE(QLEVD))
      KERR = KERR.OR.ERROR
C
      CNTRL(QPADD) = 'PADD'
      IVALUE(QPADD) = 0
      CALL TWSETI (ERROR, LOUT, CNTRL(QPADD), IVALUE(QPADD))
      KERR = KERR.OR.ERROR
C
      CNTRL(QSSABS) = 'SSABS'
      RVALUE(QSSABS) = SSABS
      CALL TWSETR (ERROR, LOUT, CNTRL(QSSABS), RVALUE(QSSABS))
      KERR = KERR.OR.ERROR
C
      CNTRL(QSSAGE) = 'SSAGE'
      IVALUE(QSSAGE) = NJAC
      CALL TWSETI (ERROR, LOUT, CNTRL(QSSAGE), IVALUE(QSSAGE))
      KERR = KERR.OR.ERROR
C
      CNTRL(QSSREL) = 'SSREL'
      RVALUE(QSSREL) = SSREL
      CALL TWSETR (ERROR, LOUT, CNTRL(QSSREL), RVALUE(QSSREL))
      KERR = KERR.OR.ERROR
C
      CNTRL(QSTDY) = 'STEADY'
      LVALUE(QSTDY) = .TRUE.
      CALL TWSETL (ERROR, LOUT, CNTRL(QSTDY), LVALUE(QSTDY))
      KERR = KERR.OR.ERROR
C
      CNTRL(QSTPS0) = 'STEPS0'
      IVALUE(QSTPS0) = NINIT
      CALL TWSETI (ERROR, LOUT, CNTRL(QSTPS0), IVALUE(QSTPS0))
      KERR = KERR.OR.ERROR
C
      CNTRL(QSTPS1) = 'STEPS1'
      IVALUE(QSTPS1) = LIMIT
      CALL TWSETI (ERROR, LOUT, CNTRL(QSTPS1), IVALUE(QSTPS1))
      KERR = KERR.OR.ERROR
C
      CNTRL(QSTPS2) = 'STEPS2'
      IVALUE(QSTPS2) = IRETIR
      CALL TWSETI (ERROR, LOUT, CNTRL(QSTPS2), IVALUE(QSTPS2))
      KERR = KERR.OR.ERROR
C
      CNTRL(QSTRID) = 'STRID0'
      RVALUE(QSTRID) = TSTEP0
      CALL TWSETR (ERROR, LOUT, CNTRL(QSTRID), RVALUE(QSTRID))
      KERR = KERR.OR.ERROR
C
      CNTRL(QTDABS) = 'TDABS'
      RVALUE(QTDABS) = TDABS
      CALL TWSETR (ERROR, LOUT, CNTRL(QTDABS), RVALUE(QTDABS))
      KERR = KERR.OR.ERROR
C
      CNTRL(QTDAGE) = 'TDAGE'
      IVALUE(QTDAGE) = ITJAC
      CALL TWSETI (ERROR, LOUT, CNTRL(QTDAGE), IVALUE(QTDAGE))
      KERR = KERR.OR.ERROR
C
      CNTRL(QTDEC) = 'TDEC'
      RVALUE(QTDEC) = 2.0
      CALL TWSETR (ERROR, LOUT, CNTRL(QTDEC), RVALUE(QTDEC))
      KERR = KERR.OR.ERROR
C
      CNTRL(QTDREL) = 'TDREL'
      RVALUE(QTDREL) = TDREL
      CALL TWSETR (ERROR, LOUT, CNTRL(QTDREL), RVALUE(QTDREL))
      KERR = KERR.OR.ERROR
C
      CNTRL(QTINC) = 'TINC'
      RVALUE(QTINC) = 2.2
      CALL TWSETR (ERROR, LOUT, CNTRL(QTINC), RVALUE(QTINC))
      KERR = KERR.OR.ERROR
C
      CNTRL(QTMAX) = 'TMAX'
      RVALUE(QTMAX) = DTMAX
      CALL TWSETR (ERROR, LOUT, CNTRL(QTMAX), RVALUE(QTMAX))
      KERR = KERR.OR.ERROR
C
      CNTRL(QTMIN) = 'TMIN'
      RVALUE(QTMIN) = DTMIN
      CALL TWSETR (ERROR, LOUT, CNTRL(QTMIN), RVALUE(QTMIN))
      KERR = KERR.OR.ERROR
C
      CNTRL(QTOL0) = 'TOLER0'
      RVALUE(QTOL0) = 1.0E-8
      CALL TWSETR (ERROR, LOUT, CNTRL(QTOL0), RVALUE(QTOL0))
      KERR = KERR.OR.ERROR
C
      CNTRL(QTOL1) = 'TOLER1'
      RVALUE(QTOL1) = 0.5
      CALL TWSETR (ERROR, LOUT, CNTRL(QTOL1), RVALUE(QTOL1))
      KERR = KERR.OR.ERROR
C
      CNTRL(QTOL2) = 'TOLER2'
      RVALUE(QTOL2) = 0.5
      CALL TWSETR (ERROR, LOUT, CNTRL(QTOL2), RVALUE(QTOL2))
      KERR = KERR.OR.ERROR
C
      IF (KERR) RETURN
C*****precision > double
      ABSOL = SQRT(D1MACH(4))
      RELAT = SQRT(D1MACH(4))
C*****END precision > double
C*****precision > single
C      ABSOL = SQRT(R1MACH(4))
C      RELAT = SQRT(R1MACH(4))
C*****END precision > single
      SFLR = -1.0E-5
C
C     set the bound on the TWOPNT solution
      DO 10 K = 1, COMPS
         ABOVE(K) = 1.1
         BELOW(K) = SFLR
   10 CONTINUE
C
C     set the initial guess on the TWOPNT solution
      IF (LUPPER) THEN
         JEND = JJ
      ELSE
         JEND = 1
      ENDIF
C
      DO 17 K = 1, KKGAS
         X(K) = S(NYS+K,JEND)
   17 CONTINUE
C
C     surface temperature
      TSURF = S(NT, JEND)
C
      IF (KKSURF .GT. 0) THEN
        DO 16 K = KFIRST(NFSURF),KLAST(NLSURF)
           X(K) = Z(K-KKGAS)
   16   CONTINUE
      ENDIF
C
C///  CHOOSE THE TWOPNT VERSION.
C
C*****precision > double
      VERSIO = 'DOUBLE PRECISION VERSION 3.22'
C*****END precision > double
C*****precision > single
C      VERSIO = 'SINGLE PRECISION VERSION 3.22'
C*****END precision > single
C
C///  TOP OF THE LOOP OVER CALLS TO TWOPNT.
C

      SIGNAL = ' '
500   CONTINUE
      CALL TWOPNT (ERROR, TEXT, VERSIO, ABOVE, LDUM1, BELOW,
     1             BUFFER, COMPS, CONDIT, 0, 0, 3, ITPWRK,
     2             LDUM2, KNAM, COMPS, 1, 1, REPORT, RSIZE,
     3             RTPWRK, SIGNAL, TSTEP1, TIME, X, RDUMMY)
      IF (ERROR) THEN
         KERR = .TRUE.
         RETURN
      ENDIF

C///  TOP OF THE BLOCK TO SERVICE REQUESTS FROM TWOPNT.

      IF (SIGNAL .NE. ' ') THEN
C
         IF (SIGNAL .EQ. 'RESIDUAL') THEN
C           evaluate the residual
C
            LVARTP = .TRUE.
            CALL CRSFUN (KKGAS, KKSURF, KKBULK, BUFFER, ISKWRK,
     1                   RSKWRK, IPAR, RPAR, TSURF, P, SDOT, SITDOT,
     2                   NNSURF, TIME, LVARTP, FZ, ZN, SDEN, ACT,
     3                   NFSURF, NLSURF, KKPHAS, KFIRST, KLAST, TSTEP1,
     4                   WT, YAV, S, JJ, NATJ, LTDIF, IMULT, GRID,
     5                   TMASS, ICRD, IMCWRK, RMCWRK, XAV, XMF,
     6                   VISC, COND, D, TDC, YV, GFAC, SFAC, LUPPER,
     7                   LVCOR)
C
C*****precision > double
           CALL DCOPY (COMPS, FZ, 1, BUFFER, 1)
C*****END precision > double
C
C*****precision > single
C           CALL SCOPY (COMPS, FZ, 1, BUFFER, 1)
C*****END precision > single
C
        ELSEIF (SIGNAL .EQ. 'PREPARE') THEN
C
           CALL CRSJAC (KKGAS, KKSURF, KKBULK, COMPS, TIME, ABSOL,
     1                  RELAT, TSURF, P, ISKWRK, RSKWRK, IPAR, RPAR,
     2                  ZN, BUFFER, FZ, FZN, A, Y, SDOT, SITDOT, NNSURF,
     3                  SDEN, ACT, NFSURF, NLSURF, KKPHAS, KFIRST,
     4                  KLAST, TSTEP1, WT, YAV, S, JJ, NATJ, LTDIF,
     5                  IMULT, GRID, TMASS, ICRD, IMCWRK, RMCWRK, XAV,
     6                  XMF, VISC, COND, D, TDC, YV, GFAC, SFAC, LUPPER,
     7                  LVCOR)
C
           IF (.NOT. XLINPK) THEN
C*****precision > double
             CALL DGECO (A, COMPS, COMPS, IPVT, RCOND, FZN)
C*****END precision > double
C*****precision > single
C             CALL SGECO (A, COMPS, COMPS, IPVT, RCOND, FZN)
C*****END precision > single
           ELSE
C*****xlinpk: precision > double
C            CALL XDGECO(A, COMPS, AA, COMPS, COMPS, IPVT, RCOND,
C     1                     DR, DC, FZN, ANORM)
C*****END xlinpk: precision > double
C*****xlinpk: precision > single
C            CALL XSGECO(A, COMPS, AA, COMPS, COMPS, IPVT, RCOND,
C     1                     DR, DC, FZN, ANORM)
C*****END xlinpk: precision > single
           ENDIF
C
            IF (RCOND .LE. 0.0E0) THEN
               WRITE (TEXT,*) ' FATAL ERROR, SINGULAR JACOBIAN '
               KERR = .TRUE.
               RETURN
            ENDIF
            CONDIT = 1.0 / RCOND
C
         ELSEIF (SIGNAL .EQ. 'SOLVE') THEN
            IF (.NOT. XLINPK) THEN
C*****precision > double
              CALL DGESL (A, COMPS, COMPS, IPVT, BUFFER, 0)
C*****END precision > double
C*****precision > single
C            CALL SGESL (A, COMPS, COMPS, IPVT, BUFFER, 0)
C*****END precision > single
            ELSE
C*****xlinpk: precision > double
C             CALL XDGESL (A, COMPS, AA, COMPS, COMPS, IPVT, DR, DC,
C     1         BUFFER, RCOND, XTEMP, RELERR, INFO, RTEMP, ANORM)
C*****END xlinpk: precision > double
C*****xlinpk: precision > single
C             CALL XSGESL (A, COMPS, AA, COMPS, COMPS, IPVT, DR, DC,
C     1         BUFFER, RCOND, XTEMP, RELERR, INFO, RTEMP, ANORM)
C*****END xlinpk: precision > single
            ENDIF
C
        ELSEIF (SIGNAL .EQ. 'RETAIN') THEN
C
C*****precision > double
          CALL DCOPY (COMPS, BUFFER, 1, ZN, 1)
C*****END precision > double
C
C*****precision > single
C          CALL SCOPY (COMPS, BUFFER, 1, ZN, 1)
C*****END precision > single
C
        ELSEIF (SIGNAL .EQ. 'SHOW') THEN
           CALL CRSPRT (KNAM, PNAM, IPAR, RPAR, TEXT, BUFFER, ISKWRK,
     1                  RSKWRK, P, TSURF, KKGAS, KKSURF,NNSURF, NNBULK,
     2                  NFSURF, NLSURF, NFBULK, NLBULK, KKPHAS, KFIRST,
     3                  KLAST, SDOT, SITDOT, SDEN, ACT, WT, DEN, Y,
     4                  PATM, IISUR, ROP, CSKWRK, YAV, S, JJ, NATJ,
     5                  LTDIF, IMULT, GRID, TMASS, ICRD, IMCWRK,
     6                  RMCWRK, XAV, XMF, VISC, COND, D, TDC, YV, FZ,
     7                  SFAC, LUPPER, LVCOR)
        ENDIF
C
        GO TO 500
C
      ENDIF
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRWALL (XMIN, XMAX, X, TEMP)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > quad
C        IMPLICIT REAL*16 (A-H, O-Z), INTEGER (I-N)
C*****END precision > quad
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      REAL SPLX(50), SPLY(50), SPLYP(50), SPLYPP(50), SPLW(150), TMPSTR,
     1     TMPEND, YPI, YPPI, XSNGL, TSNGL
C
      COMMON /SPLN/    SPLX, SPLY, SPLYP, SPLYPP,
     1                 SPLW, ISPLN, TMPSTR, TMPEND
C
C*****precision > quad
C      IF (X .LE. XMIN) THEN
C         TEMP = QEXT(TMPSTR)
C      ELSEIF (X .GE. XMAX) THEN
C         TEMP = QEXT(TMPEND)
C      ELSE
C         IPPTS = 1
C         XSNGL = SNGLQ(X)
C         CALL SPLINT (SPLX, SPLY, SPLYPP, ISPLN, XSNGL, TSNGL,
C     1                YPI, YPPI, IPPTS, KERR)
C         TEMP = QEXT(TSNGL)
C      ENDIF
C*****END precision > quad
C*****precision > double
      IF (X .LE. XMIN) THEN
         TEMP = DBLE(TMPSTR)
      ELSEIF (X .GE. XMAX) THEN
         TEMP = DBLE(TMPEND)
      ELSE
         IPPTS = 1
         XSNGL = SNGL(X)
         CALL SPLINT (SPLX, SPLY, SPLYPP, ISPLN, XSNGL, TSNGL,
     1                YPI, YPPI, IPPTS, KERR)
         TEMP = DBLE(TSNGL)
      ENDIF
C*****END precision > double
C*****precision > single
C      IF (X .LE. XMIN) THEN
C         TEMP = TMPSTR
C      ELSEIF (X .GE. XMAX) THEN
C         TEMP = TMPEND
C      ELSE
C         IPPTS = 1
C         XSNGL = X
C         CALL SPLINT (SPLX, SPLY, SPLYPP, ISPLN, XSNGL, TSNGL,
C     1                YPI, YPPI, IPPTS, KERR)
C         TEMP = TSNGL
C      ENDIF
C*****END precision > single
C
      RETURN
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION VTRACE (S, KK)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION VTRACE (S, KK)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C  START PROLOGUE
C  In the case of LVCOR=.FALSE., use 'trace' approximation,
C  lumping all transport errors into final species
C
C  S(*) - real array, species at solution gridpoint
C  KK   - integer scalar, species count
C
C  END PROLOGUE
C
      DIMENSION S(KK)
      SUMYK = 0.0
      DO 10 K = 1, KK
         SUMYK = SUMYK + S(K)
  10  CONTINUE
      VTRACE = 1.0 - SUMYK
      RETURN
      END
C=======================================================================
C
C
C                         FILE =  daserr.f
C
C    -----------------    VERSION = 1.6
C    |    SCCS  FILE |
C    |     SUMMARY   |    CURRENT DATE = 9/17/91 at 17:18:40
C    -----------------
C                         DATE OF NEWEST DELTA = 9/17/91 at 17:18:39
C
C          SCCS file name = SCCS/s.daserr.f
C          module type    =
C          q flag         =
C=======================================================================
      SUBROUTINE DASERR (T,NEQ,H,E,NFUNCT,Y,YPRIME,WT,DELTA,IPAR,RPAR)
C-----------------------------------------------------------------------
C Outer subroutine to handle error output from dassl
C
C               Nfunct = 1 : Convergence failure
C                        2 : Time step error control failure
C                        3 : Successful completion of step in sdaini
C                            (Yprime is the calculated initial deriv.)
C                        4 : Convergence failure, due to nonnegativity
C                            constraint being violated.  E is replaced
C                            by deltaA, the negative part of the
C                            solution.
C                        5 : Other failures
C
C  Dummy Variables
      INTEGER NFUNCT, NEQ, IPAR(*)
C*****precision > double
      DOUBLE PRECISION
C*****END precision > double
C*****precision > single
C      REAL
C*****END precision > single
     1         T, H, E(NEQ), Y(NEQ), YPRIME(NEQ), WT(NEQ), DELTA(NEQ),
     2         RPAR(*)

C Creslaf Pointer Common Block
      INTEGER         NCK,  NMC,  NSK,  NWT,  NMAS, NSLB, NSRB, NX,
     1                NZ1,  NZJJ, NS,   NSP,  NYAV, NYV,  NHH,  NCON,
     2                NVIS, NCND, NRR,  NY,   NYP,  NROP, NXAV, NXM,
     3                NXP,  ND,   NTDC, NF,   NDEP, NEKR, NS1,  NS2,
     4                NQ1,  NQ2,  NQ12, NWDT, NSDT, NDAS, NCP,  NSPS,
     5                NSPE, NSTS, NDDE, NSDE, NXTW, NFZ,  NZN,  NFZN,
     6                NWDI, NSDI, NGTP, NPR,  NTWP, NABV, NBLW, NBUF,
     7                NA,   NACT, NAA,  NRTP, NDR,  NDC,  NXTP, NASP,
     8                NRSP, NATS, NRTS, NGFC, NSFC, NGASW,ICK,  IMC,
     9                ISK,  IJJ,  INAJ, IKK,  IKKT, IICR, IITD, IIML,
     *                ILRO, ILRI, ILDS, ILOU, IIFX, IPKK, IPKF, IPKL,
     1                IIDS, ISIT, IFS,  ILS,  IBLK, IFB,  ILB,  IKKS,
     2                ITPW, IIPV, IMO,  IIVC, JCCH, JSCH, JKCH, JPCH
      COMMON /CRCRCR/ NCK,  NMC,  NSK,  NWT,  NMAS, NSLB, NSRB, NX,
     1                NZ1,  NZJJ, NS,   NSP,  NYAV, NYV,  NHH,  NCON,
     2                NVIS, NCND, NRR,  NY,   NYP,  NROP, NXAV, NXM,
     3                NXP,  ND,   NTDC, NF,   NDEP, NEKR, NS1,  NS2,
     4                NQ1,  NQ2,  NQ12, NWDT, NSDT, NDAS, NCP,  NSPS,
     5                NSPE, NSTS, NDDE, NSDE, NXTW, NFZ,  NZN,  NFZN,
     6                NWDI, NSDI, NGTP, NPR,  NTWP, NABV, NBLW, NBUF,
     7                NA,   NACT, NAA,  NRTP, NDR,  NDC,  NXTP, NASP,
     8                NRSP, NATS, NRTS, NGFC, NSFC, NGASW,ICK,  IMC,
     9                ISK,  IJJ,  INAJ, IKK,  IKKT, IICR, IITD, IIML,
     *                ILRO, ILRI, ILDS, ILOU, IIFX, IPKK, IPKF, IPKL,
     1                IIDS, ISIT, IFS,  ILS,  IBLK, IFB,  ILB,  IKKS,
     2                ITPW, IIPV, IMO,  IIVC, JCCH, JSCH, JKCH, JPCH
C
C Common Blocks from the Creslaf_main
      INTEGER LNCWRK
      PARAMETER (LNCWRK=450)
      CHARACTER CWORK*16
      COMMON /CCWORK/ CWORK(LNCWRK)

C Right Now, I am only interested in convergence failures.
      IF (NFUNCT .NE. 1) RETURN

      CALL DASER1( T, NEQ, H, E, NFUNCT, Y, YPRIME ,WT, DELTA,
     1    IPAR(IJJ), IPAR(INAJ), IPAR(IKK), IPAR(IKKS), IPAR(IKKT),
     2    IPAR(IICR), IPAR(ICK), RPAR(NCK), CWORK(JCCH),
     3    IPAR(ISK), RPAR(NSK), CWORK(JSCH), CWORK(JKCH)
     2           )
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE DASER1( T, NEQ, H, E, NFUNCT, Y, YPRIME, WT, DELTA,
     1                 JJ, NATJ, KKGAS, KKSURF, KKTOT,
     2                 IICD, ICKWRK, RCKWRK, CCKWRK,
     3                 ISKWRK, RSKWRK, CSKWRK, KNAM
     4                 )
C-----------------------------------------------------------------------
C Inner subroutine to handle error output from dassl
C
C               Nfunct = 1 : Convergence failure
C                        2 : Time step error control failure
C                        3 : Successful completion of step in sdaini
C                            (Yprime is the calculated initial deriv.)
C                        4 : Convergence failure, due to nonnegativity
C                            constraint being violated.  E is replaced
C                            by deltaA, the negative part of the
C                            solution.
C                        5 : Other failures
C
C
C OUTPUT FLAG TO INDICATE REASON FOR FAILURE
      INTEGER NFUNCT
C NUMBER OF EQUATIONS
      INTEGER NEQ
C*****precision > double
C TIME
      DOUBLE PRECISION T
C current step-size
      DOUBLE PRECISION H
C Error output vector - Difference between the guess and the actual
C                       value of y at the current time.
      DOUBLE PRECISION E(NEQ)
C Current solution vector
      DOUBLE PRECISION Y(NEQ)
C Current solution time derivative
      DOUBLE PRECISION YPRIME(NEQ)
C weighting vector used in the error calculations
      DOUBLE PRECISION WT(NEQ)
C Current value of delta X - For Nfunct = 1, this is the change in
C                            the solution vector from the last
C                            modified newton step in dassl.
      DOUBLE PRECISION DELTA(NEQ)
C*****END precision > double
C*****precision > single
C      REAL T, H, E(NEQ), Y(NEQ), YPRIME(NEQ), WT(NEQ), DELTA(NEQ)
C*****END precision > single

      INTEGER JJ, NATJ, KKGAS, KKSURF, KKTOT, IICD, ICKWRK, ISKWRK
C*****precision > double
      DOUBLE PRECISION
C*****END precision > double
C*****precision > single
C      REAL
C*****END precision > single
     1          RCKWRK(*), RSKWRK(*)
      CHARACTER*16 CCKWRK(*), CSKWRK(*), KNAM(KKTOT)

C Local Variables
      INTEGER I, IPTGAS, IPTTOP, NODE, LT, ID
      CHARACTER STRING_ID*30

C External Functions
      INTEGER  CKLSCH
      EXTERNAL CKLSCH

C-----------------------------------------------------------------------
C          NFUNCT = 1 Convergence Failures.
C-----------------------------------------------------------------------
      IF (NFUNCT .EQ. 1) THEN
        IF (IICD .EQ. 0) THEN
C              PLANAR COORDINATES
          IPTGAS = 1 + KKSURF
          IPTTOP = 1 + KKSURF + JJ*NATJ
        ELSE
C              CYLINDRICAL COORDINATES
          IPTGAS = 1
          IPTTOP = 1 + JJ*NATJ
        END IF
        WRITE (6,30) T,  IPTTOP + KKSURF - 1, NATJ, JJ, KKGAS, KKSURF
        DO 10 I = 1 , NEQ
C*****precision > double
        IF (DABS(DELTA(I)/WT(I)) .GT. 0.2D0) THEN
C*****END precision > double
C*****precision > single
C        IF (ABS(DELTA(I)/WT(I)) .GT. 0.20) THEN
C*****END precision > single
C
C         (find out what I is)
          STRING_ID = ' '
          IF (I .LT. (IPTGAS-1)) THEN
            NODE = 1
            LT = CKLSCH(KNAM(I+KKGAS))
            STRING_ID = KNAM(I+KKGAS)(1:LT) // ' Site Fraction'
          ELSE IF (I .GT. (IPTGAS-1+JJ*NATJ)) THEN
            NODE = JJ
            LT = CKLSCH(KNAM(I-IPTTOP+1+KKGAS))
            STRING_ID = KNAM(I-IPTTOP+1+KKGAS)(1:LT) // ' Site Fraction'
          ELSE
            NODE = (I - IPTGAS) / NATJ + 1
            ID = I - IPTGAS + 1 - (NODE-1)*NATJ
            IF (ID .EQ. 1) THEN
              STRING_ID = 'Pressure Unknown'
            ELSE IF (ID .EQ. 2) THEN
              STRING_ID = 'Position Unknown'
            ELSE IF (ID .EQ. 3) THEN
              STRING_ID = 'Axial Velocity'
            ELSE IF (ID .EQ. 4) THEN
              STRING_ID = 'Temperature Unknown'
            ELSE IF (ID .EQ. 5) THEN
              STRING_ID = 'Mass Loss - Lower Boundary'
            ELSE IF (ID .EQ. 6) THEN
              STRING_ID = 'Mass Loss - Upper Boundary'
            ELSE
              LT = CKLSCH(KNAM(ID-6))
              STRING_ID = KNAM(ID-6)(1:LT) // ' Mass Fraction'
            END IF
          END IF
          WRITE(6,31) I, STRING_ID, NODE, DELTA(I), Y(I),
     1                Y(I)-E(I), YPRIME(I), WT(I), DELTA(I)/WT(I)
        END IF
10      CONTINUE
        WRITE(6,35)
      ELSE
C-----------------------------------------------------------------------
C          Other NFUNCT
C-----------------------------------------------------------------------
      END IF
      RETURN

30    FORMAT(/130('=')/
     1 10X,'DASERR: Convergence Failure Occurred, Time = ',G11.4,
     1 ': Guilty Unknowns Printout:'/
     1 20X,'    Total number of unknowns = ',I4/
     1 20X,'    Number of unknowns per node = ',I4/
     1 20X,'    Number of grid points = ',I4/
     1 20X,'    Number of gas-phase species = ',I4/
     1 20X,'    Number of surface-phase species = ',I4/
     1 130('=')/
     1 ' Var_Number         Identity            Node    ',
     1 ' Delta_Y       Y_new       ',
     1 ' Y_Guess      Yprime_new  Weight_Fac   Rel_Delta_Y'
     1 /130('-'))
31    FORMAT(I7,3X,A30, I4, 2X, G11.5,3X,G11.5,3X,G11.5,3X,G11.5,
     1       3X,G11.5,3X,G11.5)
35    FORMAT(130('-')/)
      END


