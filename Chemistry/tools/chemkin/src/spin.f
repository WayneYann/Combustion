C     CVS: $Revision: 1.1.1.1 $ deposited $Date: 2006/05/26 19:09:33 $
C///////////////////////////////////////////////////////////////////////
C
C     SPIN
C     ONE DIMENSIONAL LAMINAR SPINNING-DISK REACTOR CODE
C
C///////////////////////////////////////////////////////////////////////
C
C
C     WRITTEN BY:  ROBERT J. KEE - AND - GREGORY H. EVANS
C                  COMPUTATIONAL MECHANICS DEPARTMENT
C                  SANDIA NATIONAL LABORATORIES
C                  LIVERMORE, CA 94551-0969 USA
C                  (510) 294-3272 / (510) 294-2795
C
C         AND BY:  MICHAEL E. COLTRIN
C                  CHEMICAL PROCESSING SCIENCES DEPARTMENT
C                  SANDIA NATIONAL LABORATORIES
C                  ALBUQUERQUE, NM 87185-5800 USA
C                  (505) 844-7843
C
C     CHANGED BY:  ELLEN MEEKS - AND - FRAN RUPLEY
C                  COMPUTATIONAL MECHANICS DEPARTMENT
C                  SANDIA NATIONAL LABORATORIES
C                  LIVERMORE, CA 94551-0969 USA
C                  (510) 294-3669 / (510) 294-3657
C
C         AND BY:  JOSEPH F. GRCAR
C                  COMPUTATIONAL MECHANICS DEPARTMENT
C                  SANDIA NATIONAL LABORATORIES
C                  LIVERMORE, CA 94551-0969 USA
C                  (510) 294-2662
C
C         AND BY:  HARRY K. MOFFAT
C                  CHEMICAL PROCESSING SCIENCES DEPARTMENT
C                  SANDIA NATIONAL LABORATORIES
C                  ALBUQUERQUE, NM 87185-5800 USA
C                  (505) 844-6912
C
C///////////////////////////////////////////////////////////////////////
C
C     V 5.27 98/05/15 (F. Rupley)
C     1. Fix bug 176, line 2438 CALL TWSTRT had transposed arguments
C        CNTRLS, NMAX instead of NMAX, CNTRLS.
C     2. Fix bug 175, in SPINDR and SPRDKY, implement use of variable
C        NTOT to allow user to limit the number of gridpoints to be
C        less than the parameterized NMAX.
C     V 5.26 98/04/20 (E. Meeks)
C     1. Fix bug#0162: Replace call to SKHMS with call to SKHML in
C        subroutine SPPRNT, to avoid divide by zero in case some
C        surface species have 0 molecular weight.  Remove multiplication
C        by WT(K) for H(K) in two lines in loop 1840.
C     V 5.25 98/03/19 (J. Grcar)
C     1. Fix bug 159 regarding variable SOLTYP in subroutine SPREAD.
C     V 5.24 98/03/19 (J. Grcar)
C     1. Fix bug 158 by adding factor 1.283153 to lines 3711 and 3791.
C     2. Use change to remove trailing blanks from several lines.
C     V.5.23 97/10/29 (E. Meeks)
C     1. Fix bug#096b: Set ABSOL=RELAT=SQRT(D1MACH(4)), rather than
C        calculating the machine constant.
C     2. Put double-precision change blocks around above, and
C        add single-precision change block using R1MACH(4).
C     V.5.22 97/09/18 (M.E. Coltrin)
C     1. Fix bug#086: In SPFUN and SPNRFN (two occurences each),
C        remove special treatment of boundary nodes when defining
C        the velocity scaling factor, VELOC; comment out IF blocks
C        that set VELOC=-S(NU,1) or VELOC=-S(NU,JJ).
C     2. Fix bug#087:  Remove unused variables: SMALL in SPPRNT;
C        STEADY in SPREAD. (E. Meeks)
C     V.5.21 97/07/29 (F. Rupley)
C     1. Fix bug#072a: Add dimensioning for CP and H in SPINDR
C     2. Fix bug#072b: Add dimensioning for TEMISG and EMISG in SPFUN
C     V.5.20 97/07/23 (F. Rupley)
C     1. Additional fixes for bug#012: removed spaces from logical
C        operators.  Found 10 occurrences of . GT. and 5 of . NE.
C     V.5.19 97/06/12 (M. Coltrin)
C     1.  Added NNSURF and NNBULK to arguments of SPSTRT, put and
C         if-test protection around loops 1191 and 1192 to test for
C         the case of no surface or no bulk species.
C     2.  Protect against the case of no bulk species for loop 120
C         in SPRDKY.
C     V.5.18 97/06/11 (M. Coltrin)
C     1.  Changed all occurances of ". AND" to ".AND" to conform to
C         ANSI standard.
C     2.  Moved DATA statements after all variable type declarations to
C         conform to ANSI standard.
C     V.5.17 97/05/29 (M. Coltrin)
C     1.  New solution-file format, v 1.00
C     2.  Add subroutine SPSOLN to store the solution, to avoid
C         duplicating code when saving the "recover" and "solution"
C         files
C     3.  Moved call to SPSAVE from SPPNT1 to the beginning of
C         SPINDR; this made change number 4 easier.
C     4.  Added IVERSN and VERNUM to arguments of SPSAVE, and
C         saved them to the solution/recover file ahead of the
C         Chemkin-type records.
C     5.  Changed the local variable LSAVE to IUNIT in SPSAVE,
C         to make it more generic.
C     V.5.16 97/03/26
C     1. add KERR as last argument to CALL SPSENS
C     2. initialize KERR=.FALSE. in SUBROUTINE SPSENS
C     CHANGES FOR VERSION 5.15 (9/5/96 H. Moffat)
C     1.  Changed misspelled EMP to EMG in SPNRFN.
C     2.  Duplicate array declarations in several places.
C     3.  SPNORM -> Add check for division by zero.
C     4.  SPSENS -> Protect against division by zero.
C     CHANGES FOR VERSION 5.14 (4/8/96 F. Rupley)
C     1.  incorporate new twopnt v.3.22 - modify CALL TWOPNT list,
C         remake SUBROUTINE TWINIT into TWSTRT for initialization
C     2.  use some JP1/JM1 to save using J+1/J-1 in calculations,
C         and in sub loops over KK or NATJ
C     3.  incorporate N=NYS+K to save using NYS+K often in calculations
C     4.  bugfix for RADIATION LOSS at JJ
C     5.  CALL CKKFRT need P as first argument?
C     6.  DEC-machine does not allow format specification continued to
C         next line (i.e. WRITE (LOUT, '(1X,A.............,
C        $      E13.5)')
C     7.  DO 2150 loop does same as DO 2200 so deleted
C     8.  added START/END PROLOGUE for comments and put any existing
C         comments inside
C     CHANGES FOR VERSION 5.13 (2/29/96 H. Moffat)
C      1. Fixed U velocity printout -> needed a weighted average
C         of U_J and U_J-1 to calculate U(X_J) for uneven grids.
C     CHANGES FOR VERSION 5.12 (2/13/96 H. Moffat)
C      1. Fixed the KDEX pointer logic for surface and bulk species
C         in SPFUN.
C      2. Initialized the variable ZERO in SPFLJJ.
C      3. Made flux printouts at top a debugging change block.
C      4. Took out the denomenators in the total mass balance columns
C         of the utilization tables printout. Reformatted that section.
C      5. Took out ad hoc flux printout that was used for debugging
C      6. Started reverse chronological order for change comment
C         statements.
C
C   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C     CHANGES FROM VERSION 1.0
C       1. Add logical variable LPRNT to SPPRNT
C       2. Add sensitivity w.r.t. surface reaction rates
C       3. After call to SPRXSE, change LASEN to .FALSE.
C       4. Move print-out for "ISENSI" above loop 1000 in SPRXSE
C       5. Make VAX open statment "readonly" on old files
C       6. Divide transient term in mixture continuity equation by
C          RHO in subroutines SPFUN and SPNRFN
C       7. Take RHO outside the parenthesis in the transient terms in
C          the radial momentum, circumferential momentum and species
C          continuity equations in subroutines SPFUN and SPNRFN
C       8. Change signs in the radial momentum, circumferential
C          momentum and energy equations in subroutines SPFUN
C          and SPNRFN
C     CHANGES FROM VERSION 1.1
C       1. Add logical variable LCHEM and keywords NOCH and CHEM.
C       2. Let the user specify TDIF without having to use MULT
C     CHANGES FROM VERSION 1.2
C       1. Read LINK file from SKLIB.17, which includes the
C          thermodynamic properties.  This change only requires adding
C          three variables to the READ (LINKSK) XXXXXX  statement.
C       2. Add the micron/min deposition rate print statement.
C       3. Include the possibility to use the correction velocity
C          formalism instead of always lumping the transport errors
C          into the last species.  New keywords VCOR and TRCE added.
C     CHANGES FROM VERSION 1.3
C       1. Fix error below do 35 loop in SPFUN.  F(NYS+K,J)..F(NYS+K,1)
C     CHANGES FROM VERSION 1.3a
C       1. Read array lengths from linkfiles
C       2. Changed REAL*8 to DOUBLE PRECISION
C       3. Normalize surface species sensitivity coefficients in SPRXSE
C       4. Fix another error below do 35 loop in SPFUN.
C          S(NYS+K,J)..S(NYS+K,1)
C       5. Fix arguments in call to SKRNAM
C       6. Allow upper/lower case input
C     CHANGES FROM VERSION 1.4
C       1. Surface linkfile has one additional parameter to
C          indicate the number of reactions with reverse Arrhenius
C          parameters defined.
C       2. Allow mixed case species names
C     CHANGES FROM VERSION 1.5
C       1. Chemkin and Surface linkfiles have additional variables
C          to read, KERR and MAXTB
C     CHANGES FROM VERSION 1.6
C       1. Add "unix" change blocks
C       2. Get rid of non-standard comment statements
C     CHANGES FOR VERSION 2.5
C       1. Incorporate mixture and pure species
C       2. Change arrays dimensioned (NMIX) or (KKSOL) to dimension (1)
C       3. Fix mistakes in SPRXSE
C          a. Loops 210 and 1210 were normalizing the sensitivity
C             coefficients for the surface species. However, the next
C             loops needed to use the raw coefficients. Therefore,
C             moved these loops below the calculation of DRDAI and
C             renumbered them 260 and 1260.
C          b. Loops 225 and 1225 used SN(L) when L>KK when they should
C             have used SN(L-KK) in order to point at the surface
C             sensitivities.
C     CHANGES FOR VERSION 3.0
C       1. Adapted to SKLIB.30 and SKINTERP.30 for 'phases':
C          new pointers, new surface linkfile, eliminates
C          calls to SKLIB subroutines which return different arrays
C          for gas, surface, and solid species.
C     CHANGES FOR VERSION 3.2
C       1. Fixed indices on loop 105 in SPFUN.
C       2. Changed NFSUR to NFSURF in call to SPFUN in SPINDR
C       3. Limits on loop 115 in SPRXSE were incorrect
C       4. Changed upper limit on loop 215 in SPRXSE from KKTOT
C          to KKGAS+KKSURF+KKBULK
C       5. Added KKBULK to arguments of SPRXSE
C       6. Changed KPHASE to KKPHAS in SPINDR
C       7. In SPPNT changed length allocated to NSDEN from NSURF to
C          NPHASE
C     CHANGES FOR VERSION 3.3
C       1. Variable NIISTK for sticking reactions added to surface
C          linkfile
C     CHANGES FOR VERSION 3.4
C       1. Rate subroutines now require addition input SDEN(*);
C          SKROP additional output SDOT(*) so additional pointer NSDT
C       2. Correct calculations which use Avogadro's number in
C          conjunction with SDEN(*), since SDEN(*) is now moles/cm**2
C          in interpreter and linking fil.
C       3. Add UNIX change blocks (OPEN statements, 'tremain' instead
C          of ctss 'timeleft' function
C       4. Add error check for presence of keywords PROD and REAC
C       5. Remove driver.
C       6. Rate factor keywords GFAC and SFAC.
C     CHANGES FOR VERSION 3.5
C       1. Add integer array KDEX for the species of maximum
C          concentration in a phase, and integer array JDEX of length
C          NMAX
C       2. Add logical LREORD input to SPFUN from calling subroutines
C          SPBJAC, SPRXSE
C       3. It is no longer necessary to read the surface linking
C          file directly; instead, call SKLEN to find lengths for the
C          work arrays, and SKINDX to find integer pointer data.
C       4. Call list has changed for SKINDX and SKDEN; SKINDX no longer
C          returns KFIRST, KLAST, and KKPHAS arrays, so SKPKK must
C          be called to return the arrays; SKDEN requires SDEN input
C          so SKSDEN must be called in SPPRNT.  SKINDX no longer
C          required ISKWRK for input.
C       5. Changed SPPNT to initialize Chemkin, Transport, and Surface
C          work arrays and get integer constants to pass to SPINDR.
C       6. Call SKNCON in SPPNT to check for presence of non-conserved
C          sites reactions.
C       7. MM as argument in SPPNT and SPINDR, add pointers INCF
C          for elemental composition, and JECH for element names,
C          required by subroutine SPFPNT
C       8. Keyword input "REOR" for choice of reordering or not
C          in SPFUN
C       9. Change values of ABOVE and ABOVEG from 200.0 to 1000.0
C     CHANGES FOR VERSION 3.6
C       1. Change CALL CKSNUM to CALL SKSNUM in SPRDKY.
C       2. CALL SKRATI instead of SKRTI.
C       3. Keyword NMAX allows user to limit number of grid points
C     CHANGES FOR VERSION 3.61
C       1. On DO loops which involve surface species, check to make
C          sure there are surface species
C       2. On DO loops which involve bulk species, check to make sure
C          there are bulk species
C       3. Above 2 changes are rescinded because we have changed the
C          interpreter such that a site or bulk without species is
C          an error.
C       4. SDOT(*) becomes SITDOT(*) and WDOT(*) becomes SDOT(*);
C          WDOTI(*) becomes SDOTI(*), and additional array SITDTI(*)
C          is returned from SUBROUTINE SKRATI, so need new pointer
C          NSDI.
C       5. Store cklink, mclink, and sklink information, and ickwrk,
C          rckwrk, cckwrk, imcwrk, rmcwrk, iskwrk, rskwrk, and cskwrk
C          arrays in savefile
C     CHANGES FOR VERSION 3.62
C       1. Call list for SKDRDA and SKDRDC.
C     CHANGES FOR VERSION 3.63
C       1. Ellen Meeks' version of 3.62 with new keywords and variables
C          QDOT, XSRC, WSRC
C     CHANGES FOR VERSION 3.71 MADE BY JOSEPH GRCAR.
C       1. Replace subroutines SPBJAC and SPNRJC by SPJACO.  The matrix
C          is evaluated by reverse communication and is row scaled.
C       2. Remove arrays PERTRB and SSAVE (pointers NDS and NSSV) from
C          the workspace.  Add array SCALE (pointer NSCL).
C       3. Remove variables NCOL and MBAND.
C       4. Correct some indentation spacing (2 where there should be 3
C          spaces) in subroutines SPINDR and SPNRDR.
C       5. Add arguments LDA and SCALE to SPNRDR.
C       6. (F. Rupley per M. Coltrin) Additional pointer ICON for
C          array of length NPHASE, due to change in SKNCON.
C     CHANGES FOR VERSION 3.72 (F. Rupley per R. Kee)
C       1. New keyword NONR and logical variable LNONR to skip the
C          non-reacting problem
C     CHANGES FOR VERSION 3.73 (F. Rupley per R. Kee)
C       1. New keyword SPOS and variable SPOS to convert negative
C          gas and site species to a small positive value.
C       2. Add binary linkfile data to RESTORE file, as in SAVE
C          file
C      CHANGES FOR VERSION 3.74 (F. Rupley, E. Meeks per R. Kee)
C       1. Add keywords TWAL, EMIS, POWR to allow radiation energy
C          balance for determination of disk temperature.
C      CHANGES FOR VERSION 3.75 (F. Rupley per R. Kee)
C       1. Keyword RADB and calculation of initial guess for
C          disk temperature.
C      CHANGES FOR VERSION 3.76 (F. Rupley per J. Grcar)
C       1. Call list for TWOPNT requires additional input; optional
C          use of new keywords reset default values:
C          'ISTP' n - sets NINIT initial time steps before newton
C                     (default is 0)
C          'IRET' n - set retirement age IRETIR of old time step
C                     (default 50)
C          'NJAC' n - set retirement age NJAC of Jacobian during
C                     steady state newton (default 20)
C          'TJAC' n - set retirement age ITJAC of Jacobian during
C                     time stepping (default 20)
C          'DTMX' x - set maximum time step DTMAX (default 1.0E-4)
C      CHANGES FOR VERSION 3.77 (F. Rupley)
C       1. If no bulk species, do not calculate production rate
C          sensitivities DRDAI in Subroutine SPRXSE; this also
C          affects the reading of restart solutions, as there are
C          fewer records to read
C      CHANGES FOR VERSION 3.78 (F. Rupley)
C       1. In reading solution file for a restart, do not store
C          old linkfile arrays, just read over these records
C          and use new linkfiles.
C      CHANGES FOR VERSION 3.79 (F. Rupley per R. J. Kee)
C      1.  Add reordering for FZ array in SPFUN (DO 2025 loop)
C      CHANGES FOR VERSION 3.80
C      1.  Take-out calculation of temperature dependent transport
C          properties from SPNRFN
C      CHANGES FOR VERSION 3.81 (D. S. Dandy)
C      1.  In SPRDKY, normalize surface site fractions such that
C          they sum to unity in each phase.
C      2.  In SPFUN, normalize surface residuals FZ(K) by site
C          density of each phase, SDEN(N).
C      CHANGES FOR VERSION 3.82 (F. RUPLEY PER ELLEN MEEKS)
C      1.  Change VSUR2 (VSUR**2) to VSUR in SPFUN.
C      CHANGES FOR VERSION 3.83 (F. RUPLEY PER BOB KEE)
C      1.  Add re-gridding option
C      CHANGES FOR VERSION 3.84 (F. RUPLEY PER BOB KEE)
C      1.  In SPRXSE, changed CALL CKRAEX to CALL CKRDEX, and
C          changed CALL SKRAEX to SKRDEX, for
C          perturbation factors used in computing sensitivities.
C      2.  Add Subroutine SPHXSE to generate sensitivity to
C          heat of formation for species.
C      CHANGES FOR VERSION 3.85 (F. Rupley per M. Coltrin)
C      1.  Allow regridding immediately after reading a restart
C          file.
C      CHANGES FOR VERSION 3.86 (F. Rupley per R. Kee)
C      1.  Changes in SPHXSE
C      CHANGES FOR VERSION 4.01 (F. Rupley per J. Grcar)
C      1.  Much restructuring due to new TWOPNT.300
C      CHANGES FOR VERSION 4.02 (F. Rupley per R. Kee & M. Coltrin)
C      1.  Reordering for grid node #1 in SPFUN (M. Coltrin)
C      2.  Injection phase keywords INJM,INJX,INJW,INJS,INJT and
C          variables FINJ, XINJ, WINJ, YINJ(K) for K=1,KKGAS,TINJ
C      CHANGES FOR VERSION 4.03 (F. Rupley/E. Meeks)
C      1.  Initialize variable XSRC to 0.0 in SPRDKY, correct NTOT
C          in SPPNT.
C      CHANGES FOR VERSION 4.04 (F. Rupley/E. Meeks)
C      1.  Correct bug in SPRDKY regarding TINJ initialization.
C      CHANGES FOR VERSION 4.05 (8/28/91 J. Grcar)
C      1.  Correct logic for LNONR
C      CHANGES FOR VERSION 4.06
C      1.  CALL TWPREP and TWSOLV instead of SPJACO. (9/5/91 J. Grcar)
C      2.  Correct injection energy equation, requires additional
C          pointer NINH for the array HINJ(K) (9/5/91 F. Rupley per Kee)
C      CHANGES FOR VERSION 4.07
C      1.  Change some twopnt controls temporarily in SPNRDR.
C      CHANGES FOR VERSION 4.08 (9/17/91 F. Rupley)
C      1.  Further correction of LNONR logic in SPNRDR.
C      CHANGES FOR VERSION 4.09 (1/15/92 D. S. Dandy & E. Meeks)
C      1.  Allow for radiation exchange between substrate and disk
C          of temperature TRAD and emmisivity ERAD. Still allow for
C          radiation between substrate and black body of temperature
C          TWAL. Keyword RADB must be used. New keywords in SPRDKY are
C            ERAD - Emissivity of radiation disk.
C            TRAD - Temperature of radiation disk (K).
C            RRAD - Ratio of radiation disk radius to distance between
C                   radiation disk and substrate (R_1 / a).
C            RDSK - Ratio of substrate radius to distance between
C                   substrate and radiation disk (R_2 / a).
C      2.  Allow for conduction heat loss through substrate. Keyword
C          RADB must be used. New keywords in SPRDKY are
C            CDCT - Indicates inclusion of substrate conduction.
C            COND - Thermal conductivity of substrate (W/cm K). If not
C                   given, molydbenum is assumed.
C            CNDT - Temperature of boundary at far side of substrate,
C                   that is, the side away from the growing surface (K).
C            CNDX - Substrate thickness, or conduction path length (cm).
C      3.  Allow for gas swirl at inlet. New keyword in SPRDKY is
C            OINL - Swirl rate at inlet (rpm).
C      4.  Modify SPRDKY to set intial site fractions to same value
C          in each phase in the event that keyword SURF not present.
C          Previously, intial site fractions were assumed to be zero
C          if SURF not used.
C      5.  Include temperature-dependent emissivities for gas radiation
C          to disk and inlet surfaces.  Additional Keywords are:
C            EMSG - gas emissivity; should be followed by a temperature
C                   and emissivity value (should give at least 2 sets).
C      6.  Fixed printout of heat fluxes to correct W/m2 units and
C          made sure heat fluxes are printed with no surface or bulk
C          sites.
C      CHANGES for V.4.10 (2/5/92 F. Rupley per J. Grcar)
C      1.  Character array TPNAM is new input to Twopnt.308 for solution
C          (groupA + comps + groupB = surface species names +
C          NYS solution names + gas species names)
C      CHANGES FOR V.4.11 (3/2/92 F. Rupley per M. Coltrin)
C      1.  Add H to SPPRNT call list, print heat release at the surface
C      CHANGES FOR V.4.20 (4/9/92 M. Coltrin)
C      1.  Add conservation equations for bulk species, surface phases.
C          The solution data structure changes with this version.
C      2.  Added logic for reordering, etc. for gas phase species at
C          node number one, in SPFUN.
C      3.  Made numerous cosmetic changes, i.e. changing line numbers
C          throughout the code.
C      4.  New keywords SDEN and ETCH.
C      5.  A number of new pointers for the work arrays (changes common
C          block SPSPSP).
C      6.  New keyword NSDN and associated logical variable LSDEN.
C      7.  Replaced KKSURF in all calls to TW* with a more general
C          LGRA (the length of "group A").
C      8.  Bulk phase activities are computed solely in the subroutine
C          SPACT.  The user can change this subroutine to incorporate
C          different models for the bulk phase activities.  By default
C          a perfect solution model is assumed in SPACT.
C      CHANGES FOR V.4.21 (4/28/92 F. Rupley)
C      1.  In CALL SPACT and SUBROUTINE SPACT variable list,
C          changed KK to KKTOT.
C      2.  In CALL SPFUN and SUBROUTINE SPFUN variable list, add
C          NPHASE as last argument.
C      3.  In SUBROUTINE SPFUN variable list, ad KKPHAS after SDEN0
C          (was in CALL SPFUN list, but not in SUBROUTINE SPFUN)
C      4.  In CALL SPRXSE, CALL SPHXSE, SUBROUTINE SPRXSE and
C          SUBROUTINE SPHXSE variable list, need NPHASE (for CALL
C          SPACT)
C      CHANGES FOR V.4.21B (7/1/92 F. Rupley per R. Kee)
C      1.  LOGICAL LCOMP to flag use of inlet mass fractions
C          of (default) inlet flux fractions
C      2.  SUBROUTINE SPPRNT now prints utilization indices
C          (required additional arrays and pointers)
C      CHANGES FOR V.4.21C (7/12/92 E. Meeks)
C      1.  Change SPFUN to use CKLIB version 3.7 routines CKKFRT, CKWYPK
C          in place of CKWYP when possible.
C      2.  Add 2nd dimension to IVALUE, RVALUE, and LVALUE (TWOPNT
C          arrays) to include all call parameters (like PREMIX.28)
C      3.  Change SPFUN to evaluate reaction rates only when T or X
C          is perturbed during Jacobian preparation.
C      CHANGES FOR V.4.21G (9/30/92 F. Rupley per M. Coltrin)
C      1.  Change in logic (IF (LENRGY .AND. LTDIF)...) in SPTRNP
C          saves much time.
C      2.  Logical LPUTIL (keyword PUTL), default=.FALSE., skip
C          printing utiliztion tables unless PUTIL used.
C      3.  SPSAVE (and CALL CKSAVE, MCSAVE, SKSAVE) no longer need
C          LINKCK, LINKMC, and LINKSK arguments
C      CHANGES FOR V.4.21H (10/6/92  M. Coltrin)
C      1.  Add about twenty variables to the end of the call list
C          to SPINDR in preparation for making the SPIN and ALE
C          pointers the same
C      2.  Add twenty new pointers into SPSPSP and to SPPNT
C      CHANGES FOR V.4.21I (10/8/92  M. Coltrin)
C      1.  Break SPPNT into two subroutines and SPSPSP into two
C          common blocks in preparation for merging with ALE
C      CHANGES FOR V.4.21J (10/9/92  M. Coltrin)
C      1.  Add LCALSP and SP to arguments of SPFUN
C      2.  Change logical array ETCH to integer variable IETCH
C      3.  Add arguments RU, DEP, DEPP, and FDEP to SPFUN
C      4.  Add total amount deposited to solution structure; increases
C          the number of solution components by KKBULK
C      CHANGES FOR V.4.21K (10/30/92 E. Meeks)
C      1.  Add SUBROUTINE SPFLAM to perform grid transformation for
C          flame restarts
C      2.  Add logical LFLAME and new solution pointer NUIN, used
C          only when LFLAME is .TRUE.
C      3.  Add keyword FLAM for which FLMX and FLMT are specified (the
C          position and value of the fixed temp. for flame restrt).
C          Also added keyword FLWD which returns FLMW - the flame width
C      4.  NO change to SAVE file structure
C      5.  Changed SPFUN to fix T at JFIXT when LFLAME is .TRUE.
C      6.  Renamed pointers NH=>NU, NF=>NV, NG=>NW
C      CHANGES FOR V.4.21N (4/1/93 F. Rupley and E. Meeks)
C      1.  Fixed bugs in sensitivity for LFLAM implementation
C      2.  Fixed bug in call to SPFLAM:  S(IPGAS), SN(IPGAS)
C      CHANGES FOR V.4.30 (2/9/94 M. Coltrin)
C      1.  Re-do sensitivity coefficients
C      2.  Replace the trivial equations for the second set of bulk
C          species residuals with equations to make the deposition
C          rate a dependent variable in SPFUN.
C      3.  Every call to SPFUN used a pointer of IPSDEN-1 for the
C          site densities at the old time.  This should have been
C          IPSDEN.  Error is corrected in this version.
C      4.  Change the normalization in SPNORM to use the absolute
C          value of the third argument.
C      5.  Get rid of variables: ICHR in SPINDR, and LADAPT, LSTEAD
C          in SPNRDR (unused).
C      CHANGES FOR V. 4.31 (2/18/94 M. Coltrin)
C      1.  Fix error in reading an old solution file. Got rid of two
C          lines that doubled the variable NREAD.
C      CHANGES FOR V. 4.32 (4/18/94 M. Coltrin)
C      1.  change call to SKNCON to be compatible with new version of
C          surface library (with real coefficients).
C      CHANGES FOR V. 4.40 (6/17/94 E. Meeks)
C      1.  Merged v. 4.21N with M. Coltrin's v. 4.32
C      2.  Added keyword AINL for V gradient at inlet
C      CHANGES FOR V. 4.41 (7/14/94 F. Rupley)
C      1.  Add ISKWRK to some CALL SK routines.
C      CHANGES FOR V. 4.42 (9/2/94 M. Coltrin)
C      1.  Remove arguments TINFTY, TDISK, RHOINF, AOMEGA from SPDIFV
C          (unused)
C      2.  Remove arguments TINFTY, TDISK from SPTRNP (unused)
C      CHANGES FOR V. 4.43 (9/30/94 H. Moffat)
C      1.  New discretization strategy - SPFUN and SPNRFN rewritten
C      2.  Axial velocity taken out of grid refinement strategy.
C      3.  Added in a change option for mixture-averaged diffusion
C          for small components, even when multicomponent diffusion
C          is used for large components (previous treatment was
C          numerically unstable for small components).
C      4.  Added in a change option for complete evaluation of the
C          jacobian and residual, including variation in transport
C          properties.
C      5.  Switched to solving the enthalpy formulation of the
C          energy conservation equation, so that a macroscopic
C          balance expression for the energy could hold rigorously.
C      6.  Put in a false transient for the Lambda unknown
C          (this is currently a change option).
C      7.  Took out the time derivative in the total continuity
C          equation (anelastic approximation). 6&7 produce much
C          better small delta_t behavior for the code.
C          (this is currently a change option)
C      8.  Modified Correction velocity case so that Sum_Yk condition
C          is applied at the reacting surface, anyway (change option).
C      9.  Fixed regrid algorithm, so that second node is not done
C          in a special fashion.
C      CHANGES FOR V. 5.00 (12/1/94 H. Moffat)
C      1.  SPTRNP - Unused parameters and calculations taken out.
C      2.  Got rid of some tab delimited formating - solaris is
C          finicky about seek ability on stdin and stdout
C      3.  Jumped version numbers from 4.x to 5.x to indicate
C          major changes that were made in SPFUN.
C      CHANGES FOR V. 5.01 (1/20/95 M. Coltrin)
C      1.  Add integer error flag to CKLEN,CKINIT,MCLEN,MCINIT,
C          SKLEN,SKINIT call lists
C      2.  Correction by E. Meeks for SPFLAM regridding.
C      3.  Change several CHARACTER C(*)*(*) declarations to
C          CHARACTER*(*) C(*) (for example).
C      4.  Include modified SUBROUTINE REFINE, and UPKINT, PCKINT
C          to over-ride current REFINE subroutine in the TWOPNT
C          library.
C      CHANGES FOR V. 5.02 (1/31/95 F. Rupley)
C      1.  Make all PRECISION > change blocks lower-case.
C      CHANGES FOR V. 5.03 (4/17/95 M. Coltrin)
C      1.  Add many parameters to first line of solution record:
C          POWR, CNDT, IVCOR, ICOMP, GFAC, SFAC, ICHEM, ICNDCT,
C          CNDFAC, IMULTI, ITDIF, IUINF, IENRGY, IRADB, SINJ,
C          XINJ, WINJ, FINJ, TINJ, QDOT, XSRC, WSRC, TWAL, EMIS,
C          ERAD, RDSK, RRAD, TRAD, NMAX, TEMISG, EMISG
C      2.  Add keyword SPRT to turn-off the long printout.
C      3.  Delete transient term from F(NT,1), F(NT,JJ), F(NV,1),
C          F(NV,JJ), F(NW,1), F(NW,JJ) for LCALSP=.FALSE. in SPFUN.
C      4.  Change "GT" to SPOS.GE.0 is SPINDR.
C      5.  Add logic to initialize values of SINJ in SPRDKY.
C      6.  Sometimes the integer flag for multi-component
C          transport was IMULT and in other places IMULTI.
C          Changed all occurrances to the latter.
C      7.  Remove transient term from bulk species residual
C          for LCALSP=.FALSE. in SPFUN
C      CHANGES FOR V. 5.04 (6/28/95 M. Coltrin per T. Badgwell, Rice U.)
C      1.  Remove ERAD from calculation of RADIN in SPFUN and SPPRNT.
C      2.  Correct units in format statement 7090
C      3.  The terms involving RRAD and RDSK in the calculation of
C          XTEMP in SPRDKY were switched.
C      CHANGES FOR V. 5.05 (7/8/95 H.K. Moffat)
C      1.  Changed top velocity specification to a total mass flux
C          specification in SPFUN
C      2.  Simplified algorithm for adding time transients in SPFUN.
C      3.  Simplified algorithm for bulk mole fraction equations in
C          SPFUN.
C      4.  Added int routine to calculate max mole fraction in a phase.
C          SPKMAX is called from SPFUN.
C      CHANGES FOR V. 5.06 (7/26/95 M. Coltrin)
C      1.  In SPFUN, the residual for the surface sites summing
C          to 1 was scaled by the variable TMP in version 5.05.
C          I changed this scaling to DXINV.
C      2.  In SPFUN, forced the residual for the JDEX(J) species
C          to be "1 - sum of the mass fractions" to be applied always.
C          This considerably simplified the logic for the species
C          equations.
C      3.  Added JDEX to the arguments of SPSTRT and initialized the
C          array to KKGAS.
C      4.  Added the integer flag IWANT to solution file.
C      5.  Move IJDX and IPMX pointers from SPONLY common block to
C          SPSPSP common block, since they will now be used in ALE;
C          move their initialization from subroutine SPPNT1 to SPPNT.
C      6.  Change the residual for F(NYS+K,JJ) to depend on the
C          desired mass fraction instead of one depending on the
C          time-derivative (for the case of LCOMP and LCALSP).
C      7.  Removed all uses of LWANT from SPFUN; it is now the
C          only option for gas, surface or bulk species.
C      8.  Changed the variable LREORD to LKMAX; changed LWANT to
C          LREORD; changed IWANT to IREORD.
C      CHANGES FOR V. 5.07 (7/30/95 M. Coltrin)
C      1.  Add keyword GDOT, for the destruction rate of gas-phase
C          species at the surface.
C      2.  Add array pointer NGDT.
C      3.  Add logical variable LGDOT.
C      4.  Add integer flag IGDOT to first record of solution.
C      5.  Add GDOT logic to SPFUN.
C      6.  Renumber loops in SPFUN.
C      7.  Add LGDOT and GDOT to arguments of SPPRNT.
C      CHANGES FOR V. 5.08 (8/4/95 M. Coltrin)
C      1.  In SPFUN everytime there was a statment like IF(LSDEN),
C          it should have read IF(LSDEN .AND. NNSURF.GT.0).
C          There were 3 occurrances of this error.
C      2.  Made the same change several other places in the code.
C      CHANGES FOR V. 5.09 (8/30/95 M. Coltrin)
C      1.  Fix logic error in SPPRNT. It would not do the LPRNT
C          print-out for gas species for the special case of KKSURF=0.
C      2.  Change some format statement in SPPRNT to give more
C          digits in the LPRT printouts.
C      3.  Add EPS and UINF to the solution file.
C      CHANGES FOR V. 5.10 (11/9/95 M. Coltrin)
C      1.  Drop "P" from the arguments of CKKFRT to be compatible
C      CHANGES FOR V. 5.11 (1/29/96 M. Coltrin)
C      1.  Fix two bugs relating to KDEX near loops 700 and 900
C          in SPFUN.
C///////////////////////////////////////////////////////////////////////
C
      SUBROUTINE SPIN (NMAX, LIN, LOUT, LINKCK, LINKSK, LINKMC, LRIN,
     $                 LROUT, LRCRVR, LFLUX, LENLWK, L, LENIWK, I,
     $                 LENRWK, R, LENCWK, C)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION I(LENIWK), R(LENRWK)
      LOGICAL L(LENLWK), KERR
      CHARACTER*(*) C(LENCWK)
C
      COMMON /SPSPSP/ NCKW, NMCW, NSKW, NEPS, NWT,  NRE,  NX,   NVIS,
     1                NCON, NTGV, NXGV, ND,   NDTC, NYV,  NS,   NACT,
     2                NFF,  NFN,  NYY,  NXAV, NDKJ, NYAV, NHHH, NCPP,
     3                NXMF, NXMP, NDDN, NSDE, NROP, NTEM, NEMG, NBLI,
     4                NSD0, NSDI, NFLX, NIFX, NRFX, NDRA, NDRC, NCDT,
     5                NDDT, NCMD, NAWT, NKA6, NWDI, NINJ, NINH, NRKF,
     6                NRKR, NRSA, NUSA, NCDF, NCNT, NPRE, NTDS, NTIN,
     7                NOME, NUIN, NRHO, NGFA, NSFA, NSWL, NQDT, NXSR,
     8                NWSR, NTWL, NEMS, NPOW, NTRD, NERD, NVFC, NRU,
     9                NINX, NINW, NITJ, NINM, NSP,  NSN,  NDP,  NDPP,
     +                NFDP, NGDT, IPST, IPND, IPKK, INCF, ICON, ICKW,
     1                IMCW, ISKW, IKR,  IETC, IPMX, IJDX, JCCH, JWCH,
     2                JKCH, JPCH, JECH
C
      COMMON /SPONLY/ NABV, NBLW, NBUF, NTWP, NA,   NPRD, NINT, NWDT,
     1                NSDT, NBUL, NFB,  NFS,  NSDN, NZ,   NZN,  NFZ,
     2                NTFL, ITWP, IKI,  IKP,  IIP,  LAC,  LMK,  JTWP
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      COMMON /SPCON / MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
C
C     write version number
      WRITE  (LOUT, 15)
   15 FORMAT (
     1/' SPIN:  1-D steady rotating disk/stagnation-point reactor',
     2/'        solving dimensional equations',
     3/'        CHEMKIN-III Version 5.27 May 15, 1998')
C
C*****precision > double
      WRITE (LOUT, '(A)') '        DOUBLE PRECISION'
C*****END precision > double
C*****precision > single
C      WRITE (LOUT, '(A)') '        SINGLE PRECISION'
C*****END precision > single
C
C     SET UP INTERNAL WORK POINTERS
C
C     get array pointers that are used only by SPIN (and not by ALE)
      KERR = .FALSE.
      CALL SPPNT1 (LINKCK, LINKMC, LINKSK, NMAX, LOUT, LENIWK,
     1            LENRWK, LENCWK, LTOT, ITOT, NTOT, JTOT, MAXNEQ,
     2            I, R, C, KERR)
      IF (KERR) RETURN
C     get the array pointers that are common to both SPIN and ALE
      CALL SPPNT ( NTOT, ITOT, JTOT, NMAX, MAXNEQ)
C
C     check for enough space
      WRITE (LOUT, 7000) LENLWK, LTOT, LENIWK, ITOT, LENRWK, NTOT,
     1                   LENCWK, JTOT
7000  FORMAT (/,'                 WORKING SPACE REQUIREMENTS',
     1        /,'                  PROVIDED        REQUIRED ',
     2        /,' LOGICAL  ' , 2I15,
     3        /,' INTEGER  ' , 2I15,
     4        /,' REAL     ' , 2I15,
     5        /,' CHARACTER' , 2I15,/)
C
      IF (LTOT.GT.LENLWK .OR. ITOT.GT.LENIWK .OR. NTOT.GT.LENRWK
     1                   .OR. JTOT.GT.LENCWK) THEN
         WRITE (LOUT, *) '  FATAL ERROR, NOT ENOUGH WORK SPACE PROVIDED'
         RETURN
      ENDIF
C
      CALL SPINDR (LIN,LOUT,LINKCK,LINKMC,LINKSK,LRIN,LROUT,LRCRVR,
     1    LFLUX,NMAX,R(NA),R(NWT),R(NXGV),R(NTGV),R(NEPS),R(NX),
     2    R(NBUF),R(NS),R(NACT),R(NFF),R(NSN),R(NFN),R(NZ),R(NZN),
     3    R(NFZ),I(ICKW),R(NCKW),I(IMCW),R(NMCW),I(ISKW),R(NSKW),
     4    C(JCCH),C(JWCH),R(NROP),I(IPST),I(IPND),I(IPKK),I(IPMX),
     5    I(IJDX),R(NVIS),R(NCON),R(ND),R(NDKJ),R(NDTC),R(NWDI),
     6    R(NSDI),R(NXAV),R(NYAV),R(NYV),R(NCPP),R(NHHH),R(NXMF),
     7    R(NXMP),R(NYY),C(JKCH),C(JPCH),C(JECH),C(JTWP),R(NDDN),
     8    R(NSDE),I(IKR),I(IKI),I(IKP),I(IIP),I(INCF),L(LAC),L(LMK),
     9    R(NRE),R(NINT),R(NPRD),R(NABV),R(NBLW),R(NTWP),I(ITWP),
     +    R(NFLX),R(NTFL),R(NIFX),R(NRFX),R(NDRA),R(NDRC),R(NWDT),
     1    R(NCDT),R(NDDT),R(NCMD),R(NAWT),R(NSDT),R(NKA6),R(NINJ),
     2    R(NINH),I(ICON),R(NTEM),R(NEMG),I(IETC), R(NBUL),R(NBLI),
     3    R(NSDN),R(NFB), R(NFS), R(NSD0), R(NUSA),R(NRSA),R(NRKF),
     4    R(NRKR),R(NCDF),R(NCNT),R(NPRE),R(NTDS),R(NTIN),R(NOME),
     5    R(NUIN),R(NRHO),R(NGFA),R(NSFA),R(NSWL),R(NQDT),R(NXSR),
     6    R(NWSR),R(NTWL),R(NEMS),R(NPOW),R(NTRD),R(NERD),R(NVFC),
     7    R(NRU),R(NINM),R(NINX),R(NINW),R(NITJ),R(NSP),R(NDP),
     8    R(NDPP),R(NFDP), R(NGDT) )
C     end of SUBROUTINE SPIN
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SPACT (ISKWRK, SKWRK, T, P, XGAS, XSURF, XBULK,
     1                  KFIRST, KLAST, ACT)
C
C  START PROLOGUE
C  Subroutine to calculate the activities of the bulk phase
c  species given the mole fractions, temperature and pressure.
C
C  This can be changed by the user.  What is implemented below
C  is the perfect solution approximation.
C
C INPUT
C
C ISKWRK  - SURFACE CHEMKIN INTEGER WORK ARRAY
C  SKWRK  - SURFACE CHEMKIN WORK ARRAY
C      T  - TEMPERATURE OF THE BULK PHASE
C             CGS UNITS - K
C      P  - PRESSURE OF THE BULK PHASE
C              CGS UNITS - ERGS/CM**2
C   XGAS  - GAS PHASE MOLE FRACTIONS
C  XSURF  - SURFACE SITE FRACTIONS
C  XBULK  - BULK PHASE MOLE FRACTIONS
C  KKTOT  - TOTAL NUMBER OF SPECIES
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
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SPCON / MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
C
      DIMENSION XGAS(KKGAS), XSURF(KKSURF), XBULK(KKBULK), ACT(KKTOT),
     1          SKWRK(*), KFIRST(NPHASE), KLAST(NPHASE), ISKWRK(*)
C
      IF (KKBULK .GT. 0) THEN
         DO 20 K = 1, KKBULK
            ACT(K) = XBULK(K)
   20    CONTINUE
      ENDIF
C
C     end of SUBROUTINE SPACT
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE SPDIFV (KKGAS, JJ, NATJ, LTDIF, LVCOR, X, S, WT, YAV,
     1                   XMF, XMFP, P,
     2                   D, DT, ICKWRK, RCKWRK, YV, DKJ, LMULTI)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
C
      DIMENSION X(JJ), S(NATJ,JJ), WT(KKGAS), YAV(KKGAS),
     1          XMF(KKGAS), XMFP(KKGAS), D(KKGAS,*), DT(KKGAS,*),
     2          DKJ(KKGAS, KKGAS, *), YV(KKGAS,*), ICKWRK(*),
     3          RCKWRK(*)
C
      LOGICAL LTDIF, LVCOR, LMULTI
C
C*****OPTIONS 3) Mixture-averaged diffusion for small components
      XABS = 1.0E-11
C*****END OPTIONS 3) Mixture-averaged diffusion for small components
C
      CALL CKYTX (S(NY1,1), ICKWRK, RCKWRK, XMFP)
C
C     Loop over all mesh points, computing the diffusion velocity at
C     the midpoints.  The indexing is such that YV(K,J) is the
C     diffusion velocity of species K midway between nodes J and J+1
C
      DO 1000 J = 1, JJ-1
C
         JP1 = J + 1
         TAV  = 0.5 * (S(NT,JP1) + S(NT,J))
         TDIF = S(NT,JP1) - S(NT,J)
         XDIF = X(JP1) - X(J)
         DO 200 K = 1, KKGAS
            N = NYS + K
            YAV(K) = 0.5 * (S(N,J) + S(N,JP1))
            XMF(K) = XMFP(K)
  200    CONTINUE
C
         CALL CKMMWY (YAV, ICKWRK, RCKWRK, WTM)
C
C        Find the mole fractions at J+1
         CALL CKYTX (S(NY1,JP1), ICKWRK, RCKWRK, XMFP)
C
         IF (LMULTI) THEN
            WTM2 = WTM**2
            DO 400 K = 1, KKGAS
               YV(K,J) = 0.0
               DO 300 L = 1, KKGAS
                  YV(K,J) = YV(K,J) + WT(L) * DKJ(K,L,J) *
     $                     (XMFP(L)-XMF(L))
  300          CONTINUE
               YV(K,J) = YV(K,J) * WT(K) / (WTM2 * XDIF)
C
C*****OPTIONS 3) Mixture-averaged diffusion for small components
C
               XMAX = MAX (XMFP(K), XMF(K))
               IF (XMAX .LT. XABS) THEN
                  DD = 0.0
                  DO 350 L = 1, KKGAS
                     DD = DD + WT(L) * DKJ(K,L,J) * XMFP(L)
  350             CONTINUE
                  DD = DD / (1.0 - XMFP(K))
                  YVMIX = - DD * WT(K) * (XMFP(K)-XMF(K))
     $                          / ( XDIF * WTM2)
                  ZEROE = 0.0
                  THETA1 = MAX(ZEROE, ((XMAX-0.1*XABS)/(0.9*XABS)))
                  YV(K,J) = YVMIX + THETA1*(YV(K,J)-YVMIX)
               ENDIF
C
C*****END OPTIONS 3) Mixture-averaged diffusion for small components
C
  400       CONTINUE
C
         ELSE
C
C           use mixture-averaged form for Fickian diffusion
            DO 500 K = 1, KKGAS
               YV(K,J) = - D(K,J) * (WT(K)/WTM) * (XMFP(K)-XMF(K))/XDIF
500         CONTINUE
         ENDIF
C
         IF (LTDIF) THEN
C           add thermal diffusion, if requested
C
            CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOAV)
            TRAV = TDIF / (TAV * RHOAV * XDIF)
            DO 600 K = 1, KKGAS
               YV(K,J) = YV(K,J) - DT(K,J) * TRAV
600         CONTINUE
C
         ENDIF
C
         IF (LVCOR) THEN
C           compute and add the correction velocity
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
C     end of SUBROUTINE SPDIFV
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE SPFLAM (FLMX, FLMT, FLMW, SNEW, SOLD, NVAR, JJNEW,
     1                   JJOLD, XNEW, XOLD, NT, NMAX, JFIXT, IERR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION SNEW(NVAR,*), SOLD(NVAR,*), XNEW(*), XOLD(*)
      LOGICAL  LBURN
      DATA J4ADD/9/, JWDTH/16/, DFRAC/.002/
      IERR = 0
C
C     Determine if previous solution is burner-stabilized:
      GRADJJ = (SOLD(NT,JJOLD-1) - SOLD(NT,JJOLD)) / SOLD(NT,JJOLD)
      LBURN = (GRADJJ .GE. 0.01)
C
C     Find TMAX location and old XFIX; adjust TFIX accordingly
      IFIX = 0
      TMAX = SOLD(NT,JJOLD)
      DO 100 J = JJOLD-1, 2, -1
         IF (SOLD(NT,J) .GT. TMAX) THEN
            TMAX = SOLD(NT,J)
            JTMAX = J
         ENDIF
         IF (IFIX .EQ. 1) GOTO 100
         IF ((SOLD(NT,J-1).GT.FLMT) .AND. (SOLD(NT,J+1).LT.FLMT)) THEN
            IFIX = 1
            FXOLD = XOLD(J)
            FTOLD = SOLD(NT,J)
            JFIXT = J
         ENDIF
100   CONTINUE
      FLMT = FTOLD
C
C     Check to make sure there's a flame
      IF (TMAX.LE.SOLD(NT,1) .OR. TMAX.LE.SOLD(NT,JJOLD) .OR.
     1        IFIX.EQ.0) THEN
         IERR = 1
         RETURN
      ENDIF
C
C     Define transformation parameters
      IF (LBURN) THEN
         J3 = JJOLD
         J4 = MIN(JJOLD+J4ADD, NMAX)
         JWDTH = JJOLD - JFIXT
         J2 = JFIXT - JWDTH
         IF (JJOLD .GE. NMAX) IERR=2
      ELSE
         J4 = JJOLD
         J3 = JFIXT + JWDTH
         J2 = JFIXT - JWDTH
         IF (J3 .GE. NMAX) IERR=3
         IF (J3 .GE. J4) J4 = MIN(J3+J4ADD, NMAX)
      ENDIF
      IF (IERR.NE.0) RETURN
C
      J1 = 0
C
      DO 200 J = 2, J2-1
         IF (SOLD(NT,J+1).GT.(.75*TMAX)
     1         .AND. SOLD(NT,J-1).LT.(.75*TMAX)) J1 = J
200   CONTINUE
      IF (J1 .EQ. 0) J1 = J2 - 1
      DELTA = DFRAC * (XOLD(JJOLD) - XOLD(1))
      H1 = XOLD(J1)
      H2 = FLMX - FLMW
      DELTA = MIN(DELTA,.25*FLMW)
      IF (H2 .LE. H1) H2 = H1 + DELTA
      H3 = FLMX + FLMW
      H4 = XOLD(JJOLD)
      IF (H3 .GE. H4) H3 = H4 - DELTA
      S1 = H1 / FLOAT(J1)
      S2 = (H3 - H2) / FLOAT(J3 - J2)
      A1 = 3.0*(H2-H1)/FLOAT(J2-J1)**2 - (2.0*S1+S2)/FLOAT(J2-J1)
      A2 = (S1+S2)/FLOAT(J2-J1)**2 - 2.0*(H2-H1)/FLOAT(J2-J1)**3
      A3 = (H4-H3)/FLOAT(J4-J3)**2 - S2/FLOAT(J4-J3)
C
C     Perform transformation of X-coordinate
      XNEW(1) = XOLD(1)
      DO 300 J = 2, J1
         XNEW(J) = S1*FLOAT(J)
         DO 295 N = 1, NVAR
            SNEW(N,J) = SOLD(N,J)
295      CONTINUE
300   CONTINUE
C
      DO 325 J = J1+1, J2
         XNEW(J) = H1 + S1*FLOAT(J-J1) + A1*(FLOAT(J-J1))**2
     1             + A2*(FLOAT(J-J1))**3
         DO 315 N = 1, NVAR
            SNEW(N,J) = SOLD(N,J)
315      CONTINUE
325   CONTINUE
C
      DO 350 J = J2+1, J3
         XNEW(J) = H2 + S2*FLOAT(J-J2)
         DO 335 N = 1, NVAR
            IF (J .LT. JJOLD) THEN
               SNEW(N,J) = SOLD(N,J)
            ELSE
               SNEW(J,J) = SOLD(N,JJOLD)
            ENDIF
335      CONTINUE
350   CONTINUE
C
      DO 375 J = J3+1, J4-1
         XNEW(J) = H3 + S2*FLOAT(J-J3) + A3*FLOAT(J-J3)**2
         DO 365 N = 1, NVAR
            IF (J.LT.JJOLD) THEN
               SNEW(N,J) = SOLD(N,J)
            ELSE
               SNEW(N,J) = SOLD(N,JJOLD)
            ENDIF
365      CONTINUE
375   CONTINUE
C
      JJNEW = J4
      XNEW(JJNEW) = XOLD(JJOLD)
      DO 380 N = 1, NVAR
         SNEW(N,JJNEW) = SOLD(N,JJOLD)
380   CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE SPFUN (LOUT, JJ, KFIRST, KLAST, LENRGY, LTDIF, LVCOR,
     1   LVARTP, FTIME, LMULTI, WT, SITSOL, DT, NTEMP, XGIVEN, TGIVEN,
     2   EPS, P, TDISK, TINFTY, AOMEGA, ASWIRL, RHOINF, UINF, LUINF,
     3   LCHEM, X, Z, SN, ACT, S, F, ZN, FZ, ICKWRK, RCKWRK, IMCWRK,
     4   RMCWRK, ISKWRK, RSKWRK, SDOT, GFAC, SFAC, SITDOT, VIS, COND, D,
     5   DKJ, DTCOEF, XAV, YAV, YV, CP, H, XMF, XMFP, LKMAX, KDEX,
     6   JDEX, QDOT, XSRC, WSRC, TWAL, EMIS, POWR, TRAD, ERAD, VFAC,
     7   LRADB, FINJ, XINJ, WINJ, TINJ, YINJ, HINJ, LCNDCT, CNDFAC,
     8   CNDT, NEMSG, TEMISG, EMISG, IETCH, BUL, BULIN, BULN, SITSLN,
     9   FB, FS, NIICON, SDEN0, KKPHAS, LSDEN, SDEN, LCOMP, USAVE,
     +   RSAVE, RKFT, RKRT, IRATE, LCALSP, ZP, BACTP, SSDENP, SP, RU,
     1   DEP, DEPP, FDEP, LFLAME, FLMT, JFIXT, ASPREAD, LGDOT, GDOT)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (PI=3.141592654D0, SIGMA = 5.67D-5, ZERO=0.0D0,
     $           ONE=1.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (PI=3.141592654, SIGMA = 5.67E-5, ZERO=0.0, ONE=1.0)
C*****END precision > single
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      COMMON /SPCON / MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
C
C     Integer arrays
      DIMENSION ICKWRK(*), IMCWRK(*), ISKWRK(*), KFIRST(NPHASE),
     1          KLAST(NPHASE), KKPHAS(NPHASE),
     2          KDEX(*), JDEX(*), NIICON(*), IETCH(*)
C     Real arrays
      DIMENSION WT(*), EPS(*), XGIVEN(*), TGIVEN(*), X(*), S(NATJ,*),
     1          SN(NATJ,*), ACT(*), F(NATJ,*), Z(*), ZN(*), FZ(*),
     2          RCKWRK(*), RMCWRK(*), RSKWRK(*), SDOT(*), SITDOT(*),
     3          XMF(*), XMFP(*), XAV(*), YAV(*),
     4          YV(KKGAS,*), CP(*), H(*), VIS(*), COND(*), D(KKGAS,*),
     5          DKJ(KKGAS,KKGAS,*), DTCOEF(KKGAS,*), SDEN(*), YINJ(*),
     7          HINJ(*), BUL(*), BULIN(*), BULN(*), SITSLN(*), FB(*),
     8          FS(*), SDEN0(*), SITSOL(*), USAVE(NATJ,*),
     9          RSAVE(KKGAS,*), RKFT(II,*), RKRT(II,*), SP(NATJ,*),
     +          ZP(*), SSDENP(*), BACTP(*), DEP(*), DEPP(*),
     1          FDEP(*), GDOT(*), TEMISG(*), EMISG(*)
C
      LOGICAL LENRGY, LTDIF, LVCOR, LVARTP, FTIME, LMULTI, LUINF,
     1        LKMAX, LCHEM, LRADB, LCNDCT, LSDEN,
     2        LCOMP, LCALSP, LFLAME, LGDOT
C
C     This variable specifies the heat flux at the top of the domain
C     rather than the temperature, when LCOMP is turned off.
C     (More consistent with the flux bc's in general)
      LOGICAL   LHEATFLX
      SAVE      LHEATFLX
C
      INTEGER   SPKMAX
      EXTERNAL  SPKMAX
C
C*****OPTIONS: Use Heat Flux BC at top of domain
C      DATA    LHEATFLX /.TRUE./
C*****END OPTIONS: Use Heat Flux BC at top of domain
C*****OPTIONS: Don't Use Heat Flux BC at top of domain
      DATA    LHEATFLX /.FALSE./
C*****END OPTIONS: Don't Use Heat Flux BC at top of domain
C
C======================================================================
C              PREPROCESSING STEPS
C======================================================================
C
C     Find a characteristic velocity for relaxation of
C     Dirchelet boundary conditions at the disk surface
CC      IF (S(NU,1) .LT. ZERO) THEN
CC        VELOC = - S(NU,1)
CC      ELSEIF (AOMEGA .GT. ZERO .OR. LUINF) THEN
      IF (AOMEGA .GT. ZERO .OR. LUINF) THEN
        VELOC = 0.5 * (ABS((X(JJ)-X(1)) * AOMEGA) + ABS(UINF))
      ELSE
        VELOC = 1000.
      ENDIF
C
C     Evaluate and store the transport coefficients
      IF (LVARTP) THEN
         IF(.NOT.LCALSP .OR. .NOT.FTIME)
     $   CALL SPTRNP (KKGAS, JJ, NATJ, LENRGY, LMULTI, LTDIF,
     $                P, X, S, YAV, ICKWRK, RCKWRK, IMCWRK, RMCWRK,
     $                XAV, COND, VIS, D, DTCOEF, DKJ)
      ENDIF
C
C     Evaluate and store the diffusion velocities
      CALL SPDIFV (KKGAS, JJ, NATJ, LTDIF, LVCOR, X, S, WT, YAV, XMF,
     $             XMFP, P, D, DTCOEF,
     $             ICKWRK, RCKWRK, YV, DKJ, LMULTI)
C
C-----------------------------------------------------------------------
C     Steady-State Equations at Node J = 1, the disk surface
C-----------------------------------------------------------------------
C
C     Calculate Mixture Densities (RHOP) at J = 1 + 1/2
      TAV = 0.5*(S(NT,1)+S(NT,2))
      DO 50 K = 1, KKGAS
         N = NYS + K
         YAV(K) = 0.5 * (S(N,1)+S(N,2))
   50 CONTINUE
      CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOP)
C
C     calculate Density (RHOJ) at the disk surface (J=1)
      CALL CKRHOY (P, S(NT,1), S(NY1,1), ICKWRK, RCKWRK, RHOJ)
C
C     calculate Density (RHOJ) at Node J+1, J = 2
      CALL CKRHOY (P, S(NT,2), S(NY1,2), ICKWRK, RCKWRK, RHOJP1)
C
C     form the differences to be used at the disk
      DXINV = 2.0 / (X(2) - X(1))
C
C     Determine the amount of injectant for current node
      FDOT = ZERO
      IF (FINJ .GT. ZERO) THEN
        XCUT = 2 * WINJ
        IF (ABS (X(1)-XINJ) .LE. XCUT) THEN
           FO = 1.283153 * FINJ * SQRT (3.0/PI) / WINJ
           FDOT = FO * EXP (-3.0*(X(1)-XINJ)**2 / WINJ**2)
        ENDIF
      ENDIF
C
      IF (LCHEM .AND. II.GT.0) THEN
C         Form the chemical rate term, SDOT, in an efficient
C         manner at the disk surface
C
          CALL SPCHEM (1, JJ, IRATE, S, P, ICKWRK, RCKWRK, SDOT,
     $                 USAVE, RKFT, RKRT, RSAVE, GFAC)
       ELSE
          DO 100 K = 1, KKGAS
            SDOT(K) = ZERO
 100      CONTINUE
       ENDIF
C
C     BC on the total continuity condition eqn at disk surface
      F(NU,1) = RHOP * S(NU,1) * DXINV + 2*RHOJ*S(NV,1) - FDOT
C
C     Boundary Condition on the radial velocity at disk surface
      F(NV,1) = RHOP * VELOC * DXINV * S(NV,1)
C
C     BC on circumferential velocity at disk surface
      F(NW,1) = RHOP * VELOC * DXINV * (S(NW,1) - AOMEGA)
C
C     Zero-devirative condition for LAMBDA at disk surface
      F(NL,1) =  VELOC * DXINV * (S(NL,1) - S(NL,2))
C
C     Determine the specific heat at the surface
      CALL CKCPBS (S(NT,1), S(NY1,1), ICKWRK, RCKWRK, CPB)
C
C     Start off the energy calculations
C
      F(NT,1) = RHOP*CPB*VELOC*DXINV * (S(NT,1) - TDISK)
      IF (LENRGY .OR. LRADB) THEN
C
C        Calculate the mean specific enthalpy of the gas
C        at the current node, J = 1
         CALL CKHBMS (S(NT,1), S(NY1,1), ICKWRK, RCKWRK, HBMSJ)
C
C        Calculate the mixture specific enthalpy at the J+1 node
         CALL CKHBMS (S(NT,2), S(NY1,2), ICKWRK, RCKWRK, HBMSJP1)
C
C        Determine the specific enthalpy of each species in
C        the mixture at J+1/2
         CALL CKHMS (TAV, ICKWRK, RCKWRK, H)
C
C        Determine the energy transport by interdiffusion
C        of chemical species at J + 1/2
         EDIFFP  = ZERO
         DO 150 K = 1, KKGAS
            EDIFFP = EDIFFP + YV(K,1) * H(K)
 150     CONTINUE
         EDIFFP = EDIFFP * RHOP
C
C        Equation for the temperature unknown at the disk surface
C
         IF (LENRGY .AND. LRADB) THEN
C           Calculate the energy source term
C
            SRCJ = ZERO
            IF (QDOT .GT. ZERO) THEN
               XCUT = 2.*WSRC
               IF (ABS(X(1) - XSRC) .LE. XCUT) THEN
                  Q0   = 1.283153 * QDOT * SQRT(3./PI) / WSRC
                  SRCJ = Q0 * EXP(-3. * (X(1)-XSRC)**2 / WSRC**2)
               ENDIF
            ENDIF
C
C           Calculate the enthalpy source term due to injector flow
            TIDOT = ZERO
            IF (FDOT .GT. ZERO) THEN
               DO 200 K = 1, KKGAS
                  TIDOT = TIDOT + YINJ(K) * HINJ(K)
 200           CONTINUE
               TIDOT = FDOT * (TIDOT - HBMSJ)
            ENDIF
C
C           Calculate the gas-phase radiation loss term
            RADJ = ZERO
            IF (NEMSG .GT. 0) THEN
               EMG = SPTEMP (NEMSG, S(NT,1), TEMISG, EMISG)
               RADJ =  EMG * SIGMA * (S(NT,1)**4 - TDISK**4)
     $               + EMG * SIGMA * (S(NT,1)**4 - TINFTY**4)
            ENDIF
C
C           Calculate the surface radiation terms
            RADIN = VFAC*SIGMA*TRAD**4
     $              + (ONE - VFAC)*SIGMA*TWAL**4
            RADOUT = EMIS*SIGMA*S(NT,1)**4 + (ONE - EMIS)*RADIN
C
            F(NT,1) = - COND(1) * (S(NT,2) - S(NT,1)) / (X(2) - X(1))
     $                                 * DXINV
     $                + DXINV * EDIFFP
     $                + RHOP*S(NU,1)*DXINV* (HBMSJP1 - HBMSJ)
     $                + DXINV * (RADOUT - RADIN - POWR)
     $                + (RADJ - TIDOT - SRCJ)
C
            IF (LCNDCT)
C              add in loss through conduction to the back side
     1         F(NT,1) = F(NT,1) + CNDFAC*DXINV*(S(NT,1) - CNDT)
         ENDIF
      ENDIF
C
C     Boundary condition on the species equations at the disk
C     surface.  The mass flux of all species needs to be
C     specified, relative to stationary coordinates.
C
      IF (LKMAX) JDEX(1) = SPKMAX (KKGAS, S(NY1,1))
      SUMYK = - ONE
      DO 250 K = 1, KKGAS
         N = NYS + K
         SUMYK = SUMYK + S(N,1)
         F(N,1) = (RHOP*YV(K,1) + RHOP*S(NU,1)*S(N,2))
     $                    * DXINV    - SDOT(K)*WT(K)
     $                    - FDOT * (YINJ(K) - S(N, 1))
  250 CONTINUE
      F(NYS+JDEX(1),1) = RHOP*VELOC*DXINV*SUMYK
C
C     BC for the INLET VELOCITY FOR FIXED FLAME POSITION
C     at the surface
C
      IF (LFLAME)
     1   F(NUL,1) = RHOP * VELOC * DXINV * (S(NUL,1) - S(NUL,2))
C
C --------------------------------------------------------------------
C  ADD IN THE SURFACE CONTRIBUTIONS TO STEADY STATE RESIDUALS AT J = 1
C --------------------------------------------------------------------
C
C Big IF BLOCK for IISUR>0
C    - If IISUR=0, then can't have any surface or bulk species or
C      deposition rates or surface site densities
C
      IF (IISUR .GT. 0) THEN
C
C        Call the SPIN surface chemistry routine that calculates
C          - SDOT
C          - SITDOT
C        by filling up the activity coefficient vector, ACT,
C        and then calling SKRAT
C
         CALL SPSCHM (ICKWRK, RCKWRK, ISKWRK, RSKWRK, S, ACT, Z,
     $                KFIRST, KLAST, P,
     $                BUL, LSDEN, SDEN, SITSOL, SDOT, SITDOT, SFAC)
C
C        Determine the gas velocity at the surface
         VSUR = ZERO
         DO 300 K = 1, KKGAS
            VSUR = VSUR + WT(K)*SDOT(K)
  300    CONTINUE
C
C        Surf Rxn contribution to BC for the continuity equation
C        at the surface
         F(NU,1) = F(NU,1) - VSUR*DXINV
C
C        Determine the contribution to the energy equation
C        due to surface production terms
C
         IF (LENRGY .AND. LRADB) THEN
C
C           Determine the partial molar enthalpies of all species
C           (gas, surface, and bulk) at the surface
C
            CALL SKHML (S(NT,1), ISKWRK, RSKWRK, H)
C
            TDOT = ZERO
            IF (NNSURF .GT. 0 .AND. KKSURF .GT. 0) THEN
               DO 350 K = KFIRST(NFSURF), KLAST(NLSURF)
                  TDOT = TDOT + SDOT(K)*H(K)
 350           CONTINUE
            ENDIF
            IF (NNBULK .GT. 0 .AND. KKBULK .GT. 0) THEN
               DO 400 K = KFIRST(NFBULK), KLAST(NLBULK)
                  TDOT = TDOT + SDOT(K)*H(K)
 400           CONTINUE
            ENDIF
C
            F(NT,1) = F(NT,1) + DXINV*(TDOT + VSUR * HBMSJ)
C
C           Add in contributions to the total increase in surface
C           sites
C
            IF (LSDEN .AND. NNSURF.GT.0) THEN
               DO 500 N = NFSURF, NLSURF
               IF (NIICON(N) .NE. 0 .AND. SITDOT(N) .NE. ZERO) THEN
C
C                 Calculate the molar enthalpy of current phase
                  HBMSJ = ZERO
                  DO 450 K = KFIRST(NFSURF), KLAST(NLSURF)
                    HBMSJ = HBMSJ + H(K)*ACT(K)
 450              CONTINUE
C
                  F(NT,1) = F(NT,1) + DXINV * SITDOT(N) * HBMSJ
               ENDIF
 500         CONTINUE
           ENDIF
         ENDIF
C
C        Determine the contribution to species continuity
C        equation at the disk surface from surface reactions
         DO 550 K = 1, KKGAS
            N = NYS + K
            F(N,1) = F(N,1) - SDOT(K)*WT(K)*DXINV
  550    CONTINUE
C        THE LAST TERM SHOULDN'T BE USED FOR THE SPECIAL EQUATION
         F(NYS+JDEX(1),1) = F(NYS+JDEX(1),1) +
     $                         SDOT(JDEX(1))*WT(JDEX(1))*DXINV
C
C        Equation for the Surface site fractions
C
         IF (KKSURF .GT. 0) THEN
C
C           Find the offset into the solution vector
C           FZ(1) refers to the first surface site fraction
C           unknown. However, Surface Chemkin refers to this
C           species as KFIRST(NFSURF).
C
            K_OFST = KFIRST(NFSURF) - 1
C
C           Find the special equations if requested
C           for the surface species unknowns. KDEX() has the same
C           values as the solution vectors FZ() and FB().
C           Therefore, FZ(KDEX(N)) points to the special
C           species in the phase, N.
C
            IF (LKMAX) THEN
               DO 600 N = NFSURF, NLSURF
                 INDEX = KFIRST(N) - KFIRST(NFSURF)
                 KDEX(N) = INDEX +
     $                     SPKMAX (KLAST(N)-KFIRST(N)+1, Z(INDEX+1))
  600          CONTINUE
            ENDIF
C
            DO 700 N = NFSURF, NLSURF
               TMP  = DXINV / SDEN(N)
               SUMZ = - ONE
               DO 650 K = KFIRST(N), KLAST(N)
                  INDEX = K - K_OFST
                  SUMZ = SUMZ + Z(INDEX)
                  FZ(INDEX) = - SDOT(K) * TMP
  650          CONTINUE
               FZ(KDEX(N)) = DXINV * SUMZ
  700       CONTINUE
         ENDIF
C
C --------------------------------------------------------------------
C    Bulk Mole Fraction and Deposition Rate Equations
C --------------------------------------------------------------------
C
C        Bulk Mole Fraction Equation
C
         IF (NNBULK .GT. 0 .AND. KKBULK .GT. 0) THEN
C
C           Find the offset into the solution vector
C           FB(1) refers to the first bulk mole fraction
C           unknown. However, Surface Chemkin refers to this
C           species as KFIRST(NFBULK).
C
            K_OFST = KFIRST(NFBULK) - 1
C
C           Do all processing as a loop over bulk phases
C
            DO 900 IPHASE = NFBULK, NLBULK
C
C              Find the special equation if requested
C              for the bulk mole fraction unknowns. KDEX() has the
C              same values as the solution vectors FZ() and FB().
C              Therefore, FZ(KDEX(IPHASE)) points to the special
C              species in the phase, IPHASE.
C
               IF (LKMAX) THEN
                  INDEX = KFIRST(IPHASE) - KFIRST(NFBULK)
                  KDEX(IPHASE) = INDEX +
     $              SPKMAX(KLAST(IPHASE)-KFIRST(IPHASE)+1,BUL(INDEX+1))
               ENDIF
C
C              The user expects this phase to be growing
C              The growth rate will be calculated as the sum
C              of all positive rates of production.
C
               IF (IETCH(IPHASE) .NE. 1) THEN
                  GROWTH = ZERO
                  DO 750 K = KFIRST(IPHASE), KLAST(IPHASE)
                     GROWTH = GROWTH + MAX(ZERO, SDOT(K))
  750             CONTINUE
               ENDIF
C
C              The form of the bulk mole fraction equations depend
C              on whether the phase is being etched or not
C
               IF (IETCH(IPHASE).EQ.1 .OR. GROWTH .LE. ZERO) THEN
C
C                 The user has specified that this bulk phase is to
C                 be etched - set the mole fraction equal to the
C                 initial mole fractions - BULIN()
C                   -or-
C                 User specified that bulk phase is to be grown,
C                 but all production rates are le zero.
C
                  SUMB = - ONE
                  DO 800 K = KFIRST(IPHASE), KLAST(IPHASE)
                     INDEX = K - K_OFST
                     SUMB = SUMB + BUL(INDEX)
                     FB(INDEX) = BUL(INDEX) - BULIN(INDEX)
  800             CONTINUE
                  FB(KDEX(IPHASE)) = SUMB
               ELSE
C
C                 At least one of the species has a positive
C                 production rate and the user expects phase
C                 to be growing
C
C                 Set species with positive production rates
C                 Eqn. 18 of Surface PSR manual.
C                    - its species mole fraction will equal the
C                      ratio of the individual production rate
C                      to the total phase production rate
C                 Species which aren't growing have their
C                 mole fractions set to zero.
C
                  SUMB = - ONE
                  DO 850  K = KFIRST(IPHASE), KLAST(IPHASE)
                     INDEX = K - K_OFST
                     SUMB = SUMB + BUL(INDEX)
                     IF (SDOT(K) .GE. ZERO) THEN
                        FB(INDEX) = BUL(INDEX) - SDOT(K)/GROWTH
                     ELSE
                        FB(INDEX) = BUL(INDEX)
                     ENDIF
  850             CONTINUE
                  FB(KDEX(IPHASE)) = SUMB
               ENDIF
  900       CONTINUE
C
C Bulk Deposition Rate Equations for all the phases at once
C    - DEP(K) is the unknown -> it is set equal to SDOT.
C
            DO 950 K = 1, KKBULK
               FDEP(K) = DEP(K) - SDOT(K+KKGAS+KKSURF)
  950       CONTINUE
C
         ENDIF
C
C        Equations for specification of the surface site
C        densities, when they are part of the solution vector
C           - At steady state, the production of each surface site
C             must be equal to zero.
C
         IF (LSDEN .AND. NNSURF.GT.0) THEN
            DO 1000 N = NFSURF, NLSURF
               IF (NIICON(N) .NE. 0) THEN
                  FS(N-1) = - SITDOT(N)
               ELSE
                  FS(N-1) = SITSOL(N-1) - SDEN0(N)
               ENDIF
 1000       CONTINUE
         ENDIF
C
C        End of BIG loop over IISUR > 0
C
      ENDIF
C
      IF (LGDOT) THEN
C        Handle the case that the destruction rates for the gas-
C        phase species at the surface are given by the user via
C        the GDOT keyword; in this case, KKSURF, KKBULK and IISUR
C        are all zero and the above big IF block will not have been
C        executed.
C
C        Determine the gas velocity at the surface
         VSUR = ZERO
         DO 1010 K = 1, KKGAS
            VSUR = VSUR + WT(K)*GDOT(K)
 1010    CONTINUE
C
C        Surf Rxn contribution to BC for the continuity equation
C        at the surface
         F(NU,1) = F(NU,1) - VSUR*DXINV
C
         IF (LENRGY.AND. LRADB)
     1      F(NT,1) = F(NT,1) + DXINV * VSUR * HBMSJ
C
C        Determine the contribution to species continuity
C        equation at the disk surface from surface reactions
C
         DO 1020 K = 1, KKGAS
            N = NYS + K
            F(N,1) = F(N,1) - GDOT(K)*WT(K)*DXINV
 1020    CONTINUE
C
C        the last term shouldn't be used for the special equation
         F(NYS+JDEX(1),1) = F(NYS+JDEX(1),1) +
     $                         GDOT(JDEX(1))*WT(JDEX(1))*DXINV
      ENDIF
C
C --------------------------------------------------------------------
C                   INTERIOR MESH POINTS
C --------------------------------------------------------------------
C
      VMAX = ZERO
C
      DO 2000 J = 2, JJ-1
C
C        Densities at J+1/2 and J-1/2
         JP1 = J + 1
         JM1 = J - 1
         TAV = 0.5*(S(NT,J)+S(NT,JP1))
         DO 1050 K = 1, KKGAS
            N = NYS + K
            YAV(K) = 0.5 * (S(N,J)+S(N,JP1))
 1050    CONTINUE
C
         RHOM = RHOP
         CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOP)
C
C        Densities at nodes J, and J + 1
         RHOJ   = RHOJP1
         CALL CKRHOY (P, S(NT,JP1), S(NY1,JP1), ICKWRK, RCKWRK, RHOJP1)
C
         IF (LCHEM .AND. II.GT.0) THEN
C
C          Form the chemical rate term, SDOT, in an efficient manner
           CALL SPCHEM (J, JJ, IRATE, S, P, ICKWRK, RCKWRK, SDOT,
     $                  USAVE, RKFT, RKRT, RSAVE, GFAC)
         ELSE
            DO 1100 K = 1, KKGAS
               SDOT(K) = ZERO
 1100       CONTINUE
         ENDIF
C
C        form the mesh differences
         DXP =        (X(JP1) - X(J)  )
         DXM =        (X(J)   - X(JM1))
         DXAV = 0.5 * (X(JP1) - X(JM1))
         DXPM =       (X(JP1) - X(JM1))
         DXINV = ONE / DXAV
C
C        form the coefficients for central differences
         CENDFM = - DXP / (DXM*DXPM)
         CENDFC =   (DXP-DXM) / (DXP*DXM)
         CENDFP =   DXM / (DXP*DXPM)
C
C        determine the amount of injectant for the current node
         FDOT = ZERO
         IF (FINJ .GT. ZERO) THEN
           XCUT = 2 * WINJ
           IF (ABS (X(J)-XINJ) .LE. XCUT) THEN
             FO = 1.283153 * FINJ * SQRT (3.0/PI) / WINJ
             FDOT = FO * EXP (-3.0*(X(J)-XINJ)**2 / WINJ**2)
           ENDIF
         ENDIF
C
C        form the common temporary variable for the convection term
         RHOJUP1 = RHOP * S(NU, J) * DXINV
C
C        mixture continuity equation
         F(NU,J) = RHOJUP1 - RHOM*S(NU,JM1)*DXINV
     $            + 2 * RHOJ * S(NV,J) - FDOT
C
C        mixture radial momentum equation
         VMAX = MAX (VMAX, ABS(S(NV,J)))
C
         F(NV,J) = (VIS(JM1)*(S(NV,J)   - S(NV,JM1))/DXM -
     $              VIS(J)  *(S(NV,JP1) - S(NV,J)  )/DXP  ) * DXINV
     $             + RHOJUP1 * (S(NV,JP1) - S(NV,J))
     $             + RHOJ * (S(NV,J)**2 - S(NW,J)**2)
     $             + S(NL,J)
C
C        mixture circumferential momentum equation
         F(NW,J) = (VIS(JM1)*(S(NW,J)   - S(NW,JM1))/DXM -
     $              VIS(J)  *(S(NW,JP1) - S(NW,J))  /DXP  ) * DXINV
     $              + RHOJUP1 * (S(NW,JP1) - S(NW,J))
     $              + 2 * RHOJ * S(NV,J) * S(NW,J)
C
C        zero-derivative condition for LAMBDA
         F(NL,J) =  VELOC * DXINV * (S(NL,J) - S(NL,JP1))
C
C        Energy Equation
C
         IF (LENRGY) THEN
C
C          Calculate the mean specific enthalpy
C          at the current node, J
           HBMSJ = HBMSJP1
C
C          Calculate the mixture specific enthalpy at the J+1 node
           CALL CKHBMS (S(NT,JP1), S(NY1,JP1), ICKWRK, RCKWRK, HBMSJP1)
C
C          Determine the specific enthalpy of each species in
C          the mixture at J+1/2
           CALL CKHMS (TAV, ICKWRK, RCKWRK, H)
C
C          Determine the energy transport by interdiffusion
C          of chemical species at J - 1/2
           EDIFFM = EDIFFP
C
C          Determine the energy transport by interdiffusion
C          of chemical species at J + 1/2
           EDIFFP  = ZERO
           DO 1150 K = 1, KKGAS
             EDIFFP = EDIFFP + YV(K,J) * H(K)
 1150      CONTINUE
           EDIFFP = EDIFFP * RHOP
C
C          Calculate the energy source term
           SRCJ = ZERO
           IF (QDOT .GT. ZERO) THEN
              XCUT = 2.*WSRC
              IF (ABS(X(J) - XSRC) .LE. XCUT) THEN
                Q0   = 1.283153 * QDOT * SQRT(3./PI) / WSRC
                SRCJ = Q0 * EXP(-3. * (X(J)-XSRC)**2 / WSRC**2)
              ENDIF
           ENDIF
C
C          Calculate the enthalpy source term due to injection
           TIDOT = ZERO
           IF (FDOT .GT. ZERO) THEN
              DO 1200 K = 1, KKGAS
                 TIDOT = TIDOT + YINJ(K) * HINJ(K)
 1200         CONTINUE
              TIDOT = FDOT * (TIDOT - HBMSJ)
           ENDIF
C
C           Calculate the radiation loss term
            RADJ = ZERO
            IF (NEMSG .GT. 0) THEN
               EMG = SPTEMP (NEMSG, S(NT,J), TEMISG, EMISG)
               RADJ =  EMG * SIGMA * (S(NT,J)**4 - TDISK**4)
     $               + EMG * SIGMA * (S(NT,J)**4 - TINFTY**4)
            ENDIF
C
C           Actually form the energy equation residual
            F(NT,J) = (COND(JM1)*(S(NT,J)   - S(NT,JM1))/DXM -
     $                 COND(J)  *(S(NT,JP1) - S(NT,J))  /DXP   )
     $                                 * DXINV
     $            + DXINV  * (EDIFFP - EDIFFM)
     $            + RHOJUP1* (HBMSJP1 - HBMSJ)
     $            + (RADJ - SRCJ - TIDOT)
C
         ELSE
C
C           note: Use the specific heat at node 1, if not
C                 solving for the energy equation
C
            F(NT,J) = RHOP * CPB * VELOC * DXINV * (S(NT,J) -
     1                SPTEMP (NTEMP, X(J), XGIVEN, TGIVEN))
C
         ENDIF
C
C        Species Conservation Equations
         IF (LKMAX) JDEX(J) = SPKMAX (KKGAS, S(NY1,J))
         SUMYK = - ONE
         DO 1250 K = 1, KKGAS
            N = NYS + K
            SUMYK = SUMYK + S(N,J)
            F(N,J) = (RHOP*YV(K,J) - RHOM*YV(K,JM1))*DXINV
     $                + RHOJUP1 * (S(N,JP1) - S(N,J))
     $                - WT(K) * SDOT(K)
     $                - FDOT * (YINJ(K) - S(N, J))
 1250    CONTINUE
         F(NYS+JDEX(J),J) = RHOP*VELOC*DXINV*SUMYK
C
C        Inlet velocity for a fixed flame position problem
         IF (LFLAME) THEN
            IF (J.GT.JFIXT) THEN
               F(NUL,J) = VELOC * DXINV * (S(NUL,J) - S(NUL,JM1))
            ELSEIF (J.EQ.JFIXT) THEN
               F(NUL,J) = VELOC * DXINV * (S(NT,J) - FLMT)
            ELSE
               F(NUL,J) = VELOC * DXINV * (S(NUL,J) - S(NUL,JP1))
            ENDIF
         ENDIF
C
 2000 CONTINUE
C
C-------------------------------------------------------------------
C           BOUNDARY CONDITIONS FAR FROM THE DISK, J = JJ
C-------------------------------------------------------------------
C
C     Density at JJ-1/2
      RHOM   = RHOP
C
C     Density at the top of the domain is equal to RHO_JJ
      CALL CKRHOY (P, S(NT,JJ), S(NY1,JJ), ICKWRK, RCKWRK, RHOJ)
C
C     Density at inlet conditions = RHOP
C         - may not be equal to RHOJ, due to Damkwertz bc causing
C           discontinuities in T and Y_k's at the inlet.
      CALL CKRHOY (P, TINFTY, EPS(1), ICKWRK, RCKWRK, RHOP)
C
C     Calculate del_x's for the top control volume
      DXM    = DXP
      DXINV  = 2.0 / (X(JJ) - X(JJ-1))
C
C     Find a characteristic positive velocity for relaxation of
C     Dirchelet boundary conditions at JJ
C      IF (S(NU,JJ) .LT. ZERO) THEN
C        VELOC = - S(NU,JJ)
C      ELSEIF (AOMEGA .GT. ZERO .OR. LUINF) THEN
      IF (AOMEGA .GT. ZERO .OR. LUINF) THEN
        VELOC = 0.5 * (ABS(X(JJ)-X(1)) * ABS(AOMEGA) + ABS(UINF))
      ELSE
        VELOC = 1000.
      ENDIF
C
      IF (LCHEM .AND. II.GT.0) THEN
C
C         Form the chemical rate term, SDOT, in an efficient
C         manner, at the top of the domain
         CALL SPCHEM (JJ, JJ, IRATE, S, P, ICKWRK, RCKWRK, SDOT,
     $                USAVE, RKFT, RKRT, RSAVE, GFAC)
       ELSE
         DO 2050 K = 1, KKGAS
           SDOT(K) = ZERO
 2050    CONTINUE
       ENDIF
C
C      Determine the amount of injectant at the top
       FDOT = ZERO
       IF (FINJ .GT. ZERO) THEN
         XCUT = 2 * WINJ
         IF (ABS (X(JJ)-XINJ) .LE. XCUT) THEN
           FO = 1.283153 * FINJ * SQRT (3.0/PI) / WINJ
           FDOT = FO * EXP (-3.0*(X(JJ)-XINJ)**2 / WINJ**2)
         ENDIF
       ENDIF
C
C     Common temporary variable for convection
      RHOJUP1 = RHOJ*S(NU,JJ)*DXINV
C
C     BC on Total continuity equation at the top
      F(NU,JJ) = RHOJUP1 - RHOM*S(NU,JJ-1)*DXINV
     $           + 2.0*RHOJ*S(NV,JJ) - FDOT
C
C     radial velocity is zero at infinity
      VMAX = MAX (VMAX, ABS(S(NV,JJ)))
C
      F(NV,JJ) = RHOJ * VELOC * DXINV * (S(NV,JJ) - ASPREAD)
C
C     circumferential velocity = swirl velocity at infinity
      F(NW,JJ) = RHOJ * VELOC * DXINV * (S(NW,JJ) - ASWIRL)
C
C     BC on the temperature at the top of the domain
C     - this may also be a flux condition, if LHEATFLX is true
C
      IF (LENRGY) THEN
C
C       determine the specific heat of the mixture at the top
        CALL CKCPBS (S(NT,JJ), S(NY1,JJ), ICKWRK, RCKWRK, CPB)
C
C       When LHEAFLX is true, the temperature at the top of
C       the domain is calculated from a balance equation.
        IF (LHEATFLX) THEN
C
C          calculate the mean specific enthalpy
C          at the current node, J = JJ
           HBMSJ = HBMSJP1
C
C          calculate the mixture specific enthalpy at the J+1=top
C          node.  This is the enthalpy of the injection mixture
C          evaluated at the temperature, TINFTY.
           CALL CKHBMS (TINFTY, EPS(1), ICKWRK, RCKWRK, HBMSJP1)
C
C          determine the energy transport by interdiffusion
C          of chemical species at J - 1/2
           EDIFFM = EDIFFP
C
C          calculate the energy source term at the top
           SRCJ = ZERO
           IF (QDOT .GT. ZERO) THEN
             XCUT = 2.*WSRC
             IF (ABS(X(JJ) - XSRC) .LE. XCUT) THEN
                Q0   = 1.283153 * QDOT * SQRT(3./PI) / WSRC
                SRCJ = Q0 * EXP(-3. * (X(JJ)-XSRC)**2 / WSRC**2)
             ENDIF
           ENDIF
C
C          calculate the injection of energy from injector at top
           TIDOT = ZERO
           IF (FDOT .GT. ZERO) THEN
             DO 2100 K = 1, KKGAS
               TIDOT = TIDOT + YINJ(K) * HINJ(K)
 2100        CONTINUE
             TIDOT = FDOT * (TIDOT - HBMSJ)
           ENDIF
C
C          calculate the gas-phase radiation loss term
           RADJ = ZERO
           IF (NEMSG .GT. 0) THEN
             EMG = SPTEMP (NEMSG, S(NT,JJ), TEMISG, EMISG)
             RADJ =  EMG * SIGMA * (S(NT,JJ)**4 - TDISK**4)
     $             + EMG * SIGMA * (S(NT,JJ)**4 - TINFTY**4)
           ENDIF
C
C          actually form the energy equation residual
           F(NT,JJ) = + COND(JJ-1) * DXINV *
     $                 (S(NT,JJ) - S(NT,JJ-1))/(X(JJ) - X(JJ-1))
     $                - DXINV * EDIFFM
     $                + RHOJUP1*(HBMSJP1 - HBMSJ)
     $                + (RADJ - SRCJ - TIDOT)
C
         ELSE
C          just specify the temperature at the top of the domain
           F(NT,JJ) = RHOJ*CPB*VELOC * DXINV * (S(NT,JJ) - TINFTY)
         ENDIF
      ELSE
         F(NT,JJ) = RHOJ * CPB * VELOC * DXINV * (S(NT,JJ) -
     1              SPTEMP(NTEMP, X(JJ), XGIVEN, TGIVEN))
      ENDIF
C
C     the LAMBDA equation determines the specified total
C     input mass flux
C
      IF (LUINF) THEN
         IF (LFLAME) THEN
            F(NL,JJ) = VMAX**2 * (RHOJ*S(NU,JJ) - RHOP*S(NUL,JJ)) /
     $                            ABS(X(JJ) - X(1))
         ELSE
            F(NL,JJ) = VMAX**2 * (RHOJ*S(NU,JJ) + RHOP*UINF) /
     $                            ABS(X(JJ) - X(1))
         ENDIF
C
      ELSE
C        input Velocity wasn't specified.  Therefore, set Lambda
C        to zero - the ideal, rotating disk solution
         F(NL,JJ) = VELOC * S(NL,JJ) / ABS(X(JJ) - X(1))
C
      ENDIF
C
      IF (LFLAME)
C       inlet velocity eigenvalue for a fixed flame problem
     1  F(NUL,JJ) =  RHOJ*VELOC * DXINV * (S(NUL,JJ) - S(NUL,JJ-1))
C
C     species continuity equation at the top of the domain
      IF (LCOMP) THEN
C        inlet mass fractions
         SUMYK = - ONE
         DO 2200 K = 1, KKGAS
            N = NYS + K
            SUMYK = SUMYK + S(N,JJ)
            F(N,JJ) = RHOJ*VELOC*DXINV*(S(N,JJ) - EPS(K))
 2200    CONTINUE
C
      ELSE
C        Specify the inlet mass flux of each species
         SUMYK = - ONE
         JJM1 = JJ - 1
         DO 2250 K = 1, KKGAS
            N = NYS + K
            SUMYK = SUMYK + S(N,JJ)
            F(N,JJ) = (RHOJ*S(NU,JJ)*(EPS(K) - S(N,JJ))
     $                      - RHOM * YV(K,JJM1)) * DXINV
     $                      - WT(K) * SDOT(K)
     $                      - FDOT * (YINJ(K) - S(N, JJ))
 2250    CONTINUE
      ENDIF
C
      IF (LKMAX) JDEX(JJ) = SPKMAX (KKGAS, S(NY1,JJ))
      F(NYS+JDEX(JJ),JJ) = RHOJ*VELOC*DXINV*SUMYK

C======================================================================
C           ADD THE TIME STEP, IF NEEDED,
C           IF NOT, RETURN
C======================================================================
C
      IF (.NOT. FTIME) RETURN
C
      TMP = ONE / DT
C     add transient term to the surface site fractions
C     conservation
C
      IF (KKSURF .GT. 0) THEN
C        precalculate the surface site fraction time derivative for
C        the SPIN case
C
         IF (LCALSP) THEN
            DO 2300 K = KFIRST(NFSURF), KLAST(NLSURF)
               INDEX = K - KFIRST(NFSURF) + 1
               ZP(INDEX) = (Z(INDEX) - ZN(INDEX)) * TMP
 2300       CONTINUE
         ENDIF
C
C        don't add a transient term to the special equation
         K_OFST = KFIRST(NFSURF) - 1
         DO 2400 N = NFSURF, NLSURF
            IF (KLAST(N) .GT. KFIRST(N)) THEN
               DO 2350 K = KFIRST(N), KLAST(N)
                  INDEX = K - K_OFST
                  IF (INDEX .NE. KDEX(N)) THEN
                     FZ(INDEX) = FZ(INDEX) + ZP(INDEX)
                  ENDIF
 2350          CONTINUE
            ENDIF
 2400    CONTINUE
      ENDIF
C
C Add a transient term to the bulk mole fraction equations and bulk
C deposition rate equations
C      - only for SPIN, ALE treats these as algebraic variables
C      - Don't add a transient to the special equation
C
      IF (LCALSP .AND. (KKBULK .GT. 0)) THEN
         K_OFST = KFIRST(NFBULK) - 1
         DO 2500 N = NFBULK, NLBULK
            DO 2450 K = KFIRST(N), KLAST(N)
               INDEX = K - K_OFST
               IF (INDEX .NE. KDEX(N)) THEN
                  FB(INDEX) = FB(INDEX) +
     $                        (BUL(INDEX) - BULN(INDEX))/DT
               ENDIF
 2450       CONTINUE
 2500    CONTINUE
      ENDIF
C
C Add transient term to deposition rate equation for the ALE case
C only
C   - ALE code uses a totally different equation which represents
C     the net deposited amount of bulk species K from the start of
C     the calculation.
C
      IF (.NOT. LCALSP) THEN
         DO 2550 K = 1, KKBULK
            FDEP(K) = DEPP(K) - SDOT(K+KKGAS+KKSURF)
 2550    CONTINUE
      ENDIF
C
C     Add a transient term to the surface site density equation
C
      IF (LSDEN.AND. NNSURF.GT.0) THEN
         DO 2600 N = NFSURF, NLSURF
            IF ( NIICON(N) .NE. 0) THEN
               IF (LCALSP) THEN
                  DSDT = (SITSLN(N-1)-SITSOL(N-1)) / DT
               ELSE
                  DSDT = SSDENP(N-1)
               ENDIF
               FS(N-1) = FS(N-1) + DSDT
            ENDIF
 2600    CONTINUE
      ENDIF
C
C-------------------------------------------------------------------
C  Transient terms for gas-phase grid points variables
C-------------------------------------------------------------------
C
C     Form the time derivative for the SPIN case, where one only has
C     an old soln.
C
      IF (LCALSP) THEN
        TMP = ONE / DT
        DO 2700 J = 1, JJ
          SP(NV,J)  = (S(NV,J) - SN(NV,J)) * TMP
          SP(NW,J)  = (S(NW,J) - SN(NW,J)) * TMP
          SP(NT,J)  = (S(NT,J) - SN(NT,J)) * TMP
          DO 2700 K = 1, KKGAS
            N = NYS + K
            SP(N,J) = (S(N,J)-SN(N,J)) * TMP
 2700   CONTINUE
      ENDIF
C
C     Loop over all grid points
C
      DO 3000 J = 1, JJ
C
C        density at the current grid point
         CALL CKRHOY (P, S(NT,J), S(NY1,J), ICKWRK, RCKWRK, RHOJ)
C
C        add the transient term to the species continuity residual
         IF (J.EQ.JJ .AND. LCOMP) THEN
C
C           for the top node and for a composition boundary
C           condition we don't want to add a transient term,
C           so do nothing
C
         ELSE
            DO 2750 K = 1, KKGAS
               N = NYS + K
               F(N,J) = F(N,J) + RHOJ * SP(N,J)
 2750       CONTINUE
C
C           A transient term shouldn't be used for the special equn.
            N = NYS + JDEX(J)
            F(N,J) = F(N,J) - RHOJ*SP(N,J)
         ENDIF
C
C         Add transient term to the mixture continuity equation
C
C
C*****TIME_DERIV - No Anelastic approximation
C
C         IF (LCALSP) THEN
C           CALL CKRHOY (P, SN(NT,J), SN(NY1,J), ICKWRK, RCKWRK, RHOJN)
C           RHODOT = (RHOJ    - RHOJN)     / DT
C         ELSE
C           CALL CKMMWY (S(NY1,J), ICKWRK, RCKWRK, WTM)
C           RHODOT = ZERO
C           DO 2800 K = 1, KKGAS
C             RHODOT = RHODOT - SP(NYS+K,J) / WT(K)
C 2800      CONTINUE
C           RHODOT = RHODOT * WTM**2 * P / (RU * S(NT,J))
C     $                 - P * WTM / (RU * S(NT,J)**2) * SP(NT,J)
C         ENDIF
C         F(NU,J) = F(NU,J) + RHODOT
C
C*****END TIME_DERIV - No Anelastic approximation
C
C
C        Add transients to the radial and circumferential velocity eqns
C
         IF (J.NE.1  .AND. J.NE.JJ) THEN
            F(NV,J) = F(NV,J) + RHOJ * SP(NV,J)
            F(NW,J) = F(NW,J) + RHOJ * SP(NW,J)
         ENDIF
C
C        Add a transient term to the Lambda equation
C             In ALE, there is no false transient for LAMBDA.
C             The top mass flux is specified as a function of time.
C
C*****TIME_DERIV - False Transient for Lambda
CC
C         IF (LCALSP) THEN
C           F(NL,J) = F(NL,J) + (S(NL,J) - SN(NL,J)) / DT
C         ENDIF
CC
C*****END TIME_DERIV - False Transient for Lambda
C
C        Add a transient term to the energy equation
C
         IF (LENRGY) THEN
C            Precalculate HDOT, when solving for the energy
C
           IF (LCALSP) THEN
             CALL CKHBMS (SN(NT,J), SN(NY1,J), ICKWRK, RCKWRK, HBMSJP1)
             CALL CKHBMS ( S(NT,J),  S(NY1,J), ICKWRK, RCKWRK, HBMSJ)
             HDOT = (HBMSJ - HBMSJP1) / DT
           ELSE
             CALL CKCPBS (S(NT,J), S(NY1,J), ICKWRK, RCKWRK, CPB)
             CALL CKHMS  (S(NT,J), ICKWRK, RCKWRK, H)
             HDOT = CPB * SP(NT,J)
             DO 2850 K = 1, KKGAS
               HDOT = HDOT + H(K)*SP(NYS+K,J)
 2850        CONTINUE
           ENDIF
C
C          special situations for boundaries and SPIN
C
           IF (J .EQ. JJ) THEN
             IF (LHEATFLX) THEN
               F(NT,JJ) = F(NT,JJ) + RHOJ * HDOT
             ELSE
               IF (LCALSP) THEN
                 CALL CKCPBS(S(NT,JJ), S(NY1,JJ), ICKWRK, RCKWRK, CPB)
                 F(NT,JJ) = F(NT,JJ) + RHOJ * CPB * SP(NT,J)
               ENDIF
             ENDIF
           ELSEIF (J .EQ. 1) THEN
             IF (LCALSP) THEN
               F(NT,J) = F(NT,J) + RHOJ*10000.*CPB*SP(NT,J)
             ENDIF
           ELSE
             F(NT,J) = F(NT,J) + RHOJ * HDOT
           ENDIF
         ELSE
C
C          Section for fixed temperature situations
           IF (J .EQ. 1) THEN
             CALL CKCPBS (S(NT,J), S(NY1,J), ICKWRK, RCKWRK, CPB)
           ENDIF
           F(NT,J) = F(NT,J) + RHOJ * CPB * SP(NT,J)
         ENDIF
C
 3000 CONTINUE
C
C     end of SUBROUTINE SPFUN
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SPCHEM
     $      (J, JJ, IRATE, S, P, ICKWRK, RCKWRK, SDOT, USAVE,
     $       RKFT, RKRT, RSAVE, GFAC)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INTEGER     J, JJ, IRATE, ICKWRK(*)
*
C*****precision > double
      DOUBLE PRECISION
C*****END precision > double
C*****precision > single
C      REAL
C*****END precision > single
     $            S(NATJ,JJ), RCKWRK(*), SDOT(KKGAS),
     $            USAVE(NATJ,JJ), RKFT(II,JJ), RKRT(II,JJ),
     $            RSAVE(KKGAS,JJ)
C
C ==============================================================
C
C SPCHEM: Form the homogeneous chemical rate terms
C
C     J = current node
C     IRATE = 1 - do a normal (expensive call to rate routines)
C             2 - Calculate and save the production rates
C             3 - Compare production rates against saved soln.
C                 Use as much of the saved production rates
C                 as possible.
C ==============================================================
C
      INTEGER       NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      INTEGER         MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
      COMMON /SPCON / MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
C
      INTEGER    I, K
      LOGICAL    MATCH, TMATCH
C
C             ========================================
C
      IF (IRATE .EQ. 1) THEN
C        normal call to rate routine - very expensive
C
        CALL CKWYP (P, S(NT,J), S(NY1,J), ICKWRK, RCKWRK, SDOT)
C
      ELSEIF (IRATE .EQ. 2) THEN
C        save all of the solution components at this node
C
        DO 305 I = 1, NATJ
          USAVE(I,J) = S(I,J)
305     CONTINUE
C
C       get the forward and reverse rate constants
        CALL CKKFRT (P, S(NT,J), ICKWRK, RCKWRK, RKFT(1,J),
     $               RKRT(1,J) )
C
C       calculate molar production rates
        CALL CKWYPK (P, S(NT,J), S(NY1,J), RKFT(1,J), RKRT(1,J),
     $               ICKWRK, RCKWRK, SDOT)
C
C       save the production rates for all species at this node
        DO 307 K = 1, KKGAS
          RSAVE(K,J) = SDOT(K)
307     CONTINUE
C
      ELSE
C
C       compare current solution against the saved solution
        TMATCH = S(NT,J) .EQ. USAVE(NT,J)
        IF (TMATCH) THEN
          MATCH = .TRUE.
          DO 310 K = 1, KKGAS
            N = NYS + K
            MATCH = MATCH .AND. (S(N,J) .EQ. USAVE(N,J))
310       CONTINUE
        ELSE
          MATCH = .FALSE.
        ENDIF
C
C       if match is true after the two tests, then neither the
C       temperature nor the species concentrations have changed,
C       so don't need to recalculate the production rates
C
        IF (MATCH) THEN
          DO 312 K = 1, KKGAS
            SDOT(K) = RSAVE(K,J)
312       CONTINUE
C
        ELSEIF (TMATCH) THEN
C         temperature has not changed, so don't need to recalculate
C         rate constants to get new production rates
C
          CALL CKWYPK (P, S(NT,J), S(NY1,J), RKFT(1,J),
     $                 RKRT(1,J), ICKWRK, RCKWRK, SDOT)
C
        ELSE
C          temperature has changed, so calculate new rate constants
C          and production rates (expensive)
C
           CALL CKWYP (P, S(NT,J), S(NY1,J), ICKWRK, RCKWRK, SDOT)
C
        ENDIF
C
      ENDIF
C
C     damp the rates, no matter how they were calculated
      IF (GFAC .NE. 1.0) THEN
        DO 315 K = 1, KKGAS
          SDOT(K) = SDOT(K)*GFAC
  315   CONTINUE
      ENDIF
C
C     end of SUBROUTINE SPCHEM
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SPSCHM (ICKWRK, RCKWRK, ISKWRK, RSKWRK, S, ACT, Z,
     $                   KFIRST, KLAST, P, BUL, LSDEN, SDEN, SITSOL,
     $                   SDOT, SITDOT, SFAC)
C
C  START PROLOGUE
C  END PROLOGUE
c
C*****precision > double
      DOUBLE PRECISION        ZERO,         ONE
      PARAMETER              (ZERO = 0.0D0, ONE = 1.0D0)
C*****END precision > double
C*****precision > single
C      REAL                    ZERO,         ONE
C      PARAMETER              (ZERO = 0.0E0, ONE = 1.0E0)
C*****END precision > single
C
C   ==================================================================
C
C     This subroutine calculates the chemistry
C     production rates on the surface.  There is some preparation
C     involved.  Therefore, it's a separate routine.
C   ==================================================================
C
C
      INTEGER       NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      INTEGER        MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
      COMMON /SPCON / MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
C
      INTEGER        ICKWRK(*), ISKWRK(*), KFIRST(NPHASE),
     $               KLAST(NPHASE)
      LOGICAL        LSDEN
*
C*****precision > double
      DOUBLE PRECISION
C*****END precision > double
C*****precision > single
C      REAL
C*****END precision > single
     $               RCKWRK(*), RSKWRK(*), S(NATJ,*), ACT(KKTOT),
     $               Z(KKTOT), BUL(*), SDEN(NPHASE), SITSOL(NPHASE),
     $               SDOT(KKTOT), SITDOT(NPHASE), P, SFAC
C
      INTEGER        K, INDEX, N
C
C      -----------------------
      IF (IISUR .LE. 0) THEN
         DO 100 K = 1, KKTOT
            SDOT(K) = ZERO
  100    CONTINUE
         DO 102 N = 2, NNSURF+1
           SITDOT(N) = ZERO
  102    CONTINUE
         RETURN
      ENDIF
C
C     store mole fractions in ACT(*) - surface is at J = 1
      CALL CKYTX (S(NY1,1), ICKWRK, RCKWRK, ACT)
C
C     store site fractions in the next KKSURF places in ACT
      IF (KKSURF .GT. 0) THEN
         DO 50 K = KFIRST(NFSURF), KLAST(NLSURF)
            INDEX = K - KFIRST(NFSURF) + 1
            ACT(K) = Z(INDEX)
   50    CONTINUE
      ENDIF
C
C     Calculate and store the bulk activities in the next
C     place (note BUL are the bulk mole fractions, and for
C     ideal solid solution, ACT(*) is equal to BUL(*))
C
      IF (KKBULK .GT. 0) THEN
        CALL SPACT (ISKWRK, RSKWRK, S(NT,1), P, ACT(1), Z(1), BUL,
     $              KFIRST, KLAST, ACT(KFIRST(NFBULK)) )
      ENDIF
C
C     Store the current site densities for the surface phases
C     when and if they are part of the solution vector
C
      IF (LSDEN .AND. NNSURF.GT.0) THEN
         DO 65 N = 1, NNSURF
            SDEN(N+1) = SITSOL(N)
   65    CONTINUE
      ENDIF
C
C     Call the surface chemkin routine that calculates the
C     production rate of species, SDOT, and the production
C     rate of sites, SITDOT.
C
      CALL SKRAT (P, S(NT,1), ACT, SDEN, ISKWRK, RSKWRK, SDOT,
     $               SITDOT)
C
C     damp the production terms, if called for
      IF (SFAC .NE. ONE) THEN
        DO 70 K = 1, KKTOT
           SDOT(K) = SFAC*SDOT(K)
   70   CONTINUE
        DO 72 N = 2, NNSURF+1
           SITDOT(N) = SFAC*SITDOT(N)
   72   CONTINUE
      ENDIF
C
C     end of SUBROUTINE SPSCHM
      RETURN
      END
C
C----------------------------------------------------------------------
C
      INTEGER FUNCTION SPKMAX (NUM, VALUE)
C
C  START PROLOGUE
C        Returns the index of the maximum value of a vector
C        (NUM >= 1)
C  END PROLOGUE
C
      INTEGER                  NUM
C*****precision > double
      DOUBLE PRECISION         VALUE(NUM)
C*****END precision > double
C*****precision > single
C      REAL                    VALUE(NUM)
C*****END precision > single
C
      INTEGER                  K, KMAX
C
      KMAX = NUM
      DO 10 K = 1, NUM - 1
        IF (VALUE(K) .GT. VALUE(KMAX)) KMAX = K
10    CONTINUE
      SPKMAX = KMAX
C
C     end of FUNCTION SPKMAX
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SPINDR (LIN, LOUT, LINKCK, LINKMC, LINKSK, LRIN, LROUT,
     1    LRCRVR, LFLUX, NMAX, A, WT, XGIVEN, TGIVEN, EPS, X, BUFFER,
     2    S, ACT, F, SN, FN, Z, ZN, FZ, ICKWRK, RCKWRK, IMCWRK, RMCWRK,
     3    ISKWRK, RSKWRK, CCKWRK, CSKWRK, ROP, KFIRST, KLAST, KKPHAS,
     4    KDEX, JDEX, VIS, COND, D, DKJ, DTCOEF, SDOTI, SITDTI, XAV,
     5    YAV, YV, CP, H, XMF, XMFP, Y, KNAM, PNAM, ENAM, TPNAM, DEN,
     6    SDEN, KR, KI, KP, IP, KNCF, ACTIVE, MARK, REAC, XINTM, PROD,
     7    ABOVE, BELOW, RTWWRK, ITWWRK, FLUX, TFLUXA, REAFLX, RADFLX,
     8    DRDAI, DRDCL, SDOT, CDOT, DDOT, WDOT, AWT, SITDOT, A6, YINJ,
     9    HINJ, NIICON, TEMISG, EMISG, IETCH, BUL, BULIN,SITSLN, FB, FS,
     +    SDEN0, USAVE, RSAVE, RKFT, RKRT, CNDFAC, CNDT, P, TDISK,
     1    TINFTY, AOMEGA, UINF, RHOINF, GFAC, SFAC, ASWIRL, QDOT, XSRC,
     2    WSRC, TWAL, EMIS, POWR, TRAD, ERAD, VFAC, RU, FINJ, XINJ,
     3    WINJ, TINJ, SP, DEP, DEPP, FDEP, GDOT)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INTEGER CNTRLS, CALL, CALLS
      PARAMETER (CNTRLS = 22)
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      COMMON /SPCON / MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
C     Integer arrays
      DIMENSION ICKWRK(*), RCKWRK(*), IMCWRK(*), RMCWRK(*),
     1          ISKWRK(*), RSKWRK(*), IVALUE(CNTRLS,3), KR(*), KI(*),
     2          KP(*), IP(*), KFIRST(*), KLAST(*), KKPHAS(*),
     3          KDEX(*), JDEX(*), KNCF(MM,KKTOT), NIICON(*),
     4          ITWWRK(*), IETCH(*)
C     Real arrays
      DIMENSION RVALUE(CNTRLS,3), EPS(*), REAC(*), XINTM(*), PROD(*),
     1          FLUX(*), TFLUXA(*), YAV(*), XAV(*), XMF(*), XMFP(*),
     2          Y(*), WT(*), SDOT(*), SITDOT(*), SDOTI(*), SITDTI(*),
     3          DRDAI(*), DRDCL(*), ACT(*), X(*), VIS(*), COND(*),
     4          D(KKGAS,*), DKJ(KKGAS,KKGAS,*), DTCOEF(KKGAS,*),
     5          YV(KKGAS,*), XGIVEN(*), TGIVEN(*), S(*), SN(*), F(*),
     6          FN(*), ABOVE(*), BELOW(*), BUFFER(*), A(*), Z(*),
     8          ZN(*), FZ(*), ROP(*), SDEN(*), RTWWRK(LENRTW),
     9          A6(*), YINJ(*), HINJ(*), TEMISG(*), EMISG(*),
     *          BUL(*), BULIN(*), SITSLN(*), FB(*), FS(*), SDEN0(*),
     1          CDOT(*), DDOT(*), WDOT(*), AWT(*), REAFLX(*),
     2          RADFLX(*), USAVE(*), RSAVE(KKGAS,*), RKFT(II,*),
     4          RKRT(II,*), SP(*), DEP(*), DEPP(*), FDEP(*), GDOT(*),
     5          CP(KKGAS), H(KKTOT)
C
      CHARACTER*(*) CSKWRK(*), CCKWRK(*), KNAM(*), PNAM(*), ENAM(*),
     1              TPNAM(*)
      CHARACTER*16 STEADY, ISENSI, IHSENS, IVERSN, VERNUM,
     1             SIGNAL, REPORT
C
      CHARACTER VERSIO*80
      LOGICAL LTIME2, LENRGY, LTDIF, LVCOR, LUMESH, LRSTRT, LCNTUE,
     1        LASEN, LVARTP, RSTCNT, ERROR, FTIME, FIRST, ACTIVE(*),
     2        MARK(*), ENERGY, LUSTGV, LMULTI, LUINF, LCHEM, LREORD,
     3        LKMAX, IERR, LFIRST, LPRNT, FLAG, LNONR, LRADB, LREGRD,
     4        LHSEN, LVALUE(CNTRLS,3), LCNDCT, LSDEN, LTOFF,
     5        LPUTIL, LCALSP, LFLAME, LOLDSP, LCOMP, LGDOT, KERR
C
      DATA STEADY/'STEADY SOLN SPIN'/, ISENSI/'SENSITIVITY     '/,
     1     IHSENS/'HSENSITIVITY    '/,
     2     IVERSN/'VERSION         '/, VERNUM/'1.00            '/
C
C*****precision > double
      DATA VERSIO/'DOUBLE PRECISION VERSION 3.22'/
C*****END precision > double
C*****precision > single
C      DATA VERSIO/'SINGLE PRECISION VERSION 3.22'/
C*****END precision > single
C
      DATA LCNTUE/.FALSE./, NSOL/0/, LCALSP/.TRUE./, ERROR/.FALSE./
C
      KERR = .FALSE.
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
      CALL SKNCON (ISKWRK, RSKWRK, NIICON)
C
      CALL SKRP   (ISKWRK, RSKWRK, RU, RUC, PATM)
      CALL SKPKK  (ISKWRK, KKPHAS, KFIRST, KLAST)
      CALL SKSYMS (ISKWRK, CSKWRK, LOUT, KNAM, IERR)
      CALL SKSYMP (ISKWRK, CSKWRK, LOUT, PNAM, IERR)
      CALL SKWT   (ISKWRK, RSKWRK, WT)
      CALL SKNCF  (MM, ISKWRK, KNCF)
      CALL SKSYME (ISKWRK, CSKWRK, LOUT, ENAM, IERR)
      CALL CKAWT  (ICKWRK, RCKWRK, AWT)
      CALL SKNCF  (MM, ISKWRK, KNCF)
C
      REWIND LROUT
      CALL SPSAVE (ICKWRK, RCKWRK, CCKWRK,
     1             IMCWRK, RMCWRK, ISKWRK, RSKWRK,
     2             CSKWRK, LOUT, LROUT, IVERSN, VERNUM)
C
  100 CONTINUE
C
C           CALL THE KEYWORD INPUT
C
      LFIRST = .TRUE.
C
      CALL SPRDKY (PNAM, KFIRST, KLAST, KKPHAS, NMAX, LIN, LOUT, PATM,
     1             LTIME2, LUSTGV, LENRGY, RSTCNT, LCHEM, LTDIF, LVCOR,
     2             LUMESH, LRSTRT, LCNTUE, LMULTI, LUINF, LASEN, LPRNT,
     3             LNONR, LRADB, SFLR, P, NPTS, NTOT, NADP, X, NREAC,
     4             NINTM, NPROD, REAC, XINTM, PROD, KR, KI, KP, Z,
     5             TINFTY, TDISK, AOMEGA, UINF, ASWIRL, XSTR, XCEN,
     6             XEND, WMIX, NTEMP, XGIVEN, TGIVEN, N1CALL, MFILE,
     7             KNAM, ACT, GFAC, SFAC, LREORD, QDOT, XSRC, WSRC,
     8             SPOS, LREGRD, JJREGD, PCTADP, RATGTC, LHSEN, IVALUE,
     9             RVALUE, LVALUE, FINJ, XINJ, WINJ, BUFFER,
     *             TINJ, HINJ, TWAL, EMIS, POWR, TRAD, ERAD, VFAC,
     1             LCNDCT, CNDFAC, CNDT, NEMSG, TEMISG, EMISG, IETCH,
     2             SDEN, SDEN0, ISKWRK, RSKWRK, BUL, BULIN, LSDEN,
     3             LCOMP, CNTRLS, LTOFF, LPUTIL, FLMX, FLMT, LFLAME,
     4             FLMW, LOLDSP, ASPREAD, GDOT, LGDOT, KERR)
      IF (KERR) RETURN
C
C
      IF (LVCOR) THEN
         IVCOR = 1
      ELSE
         IVCOR = 0
      ENDIF
C
      IF (LCOMP) THEN
         ICOMP = 1
      ELSE
         ICOMP = 0
      ENDIF
C
      IF (LCHEM) THEN
         ICHEM = 1
      ELSE
         ICHEM = 0
      ENDIF
C
      IF (LCNDCT)THEN
         ICNDCT= 1
      ELSE
         ICNDCT= 0
      ENDIF
C
      IF (LMULTI)THEN
         IMULTI= 1
      ELSE
         IMULTI= 0
      ENDIF
C
      IF (LTDIF) THEN
         ITDIF = 1
      ELSE
         ITDIF = 0
      ENDIF
C
      IF (LUINF) THEN
         IUINF = 1
      ELSE
         IUINF = 0
      ENDIF
C
      IF (LENRGY)THEN
         IENRGY= 1
      ELSE
         IENRGY= 0
      ENDIF
C
      IF (LRADB) THEN
         IRADB = 1
      ELSE
         IRADB = 0
      ENDIF
C
      IF (LREORD) THEN
         IREORD = 1
      ELSE
         IREORD = 0
      ENDIF
C
      IF (LGDOT) THEN
         IGDOT = 1
      ELSE
         IGDOT = 0
      ENDIF
C
C     NSSOLV is the number of surface site density equations
C     that will be solved
      IF (LSDEN) THEN
         NSSOLV = NNSURF
      ELSE
         NSSOLV = 0
      ENDIF
C
C     pointers to locations in solution
      IPSURF = 1
      IPBULK = IPSURF + KKSURF
      IPSDEN = IPBULK + KKBULK*2
      IPGAS  = IPSDEN + NSSOLV
C
C     TWOPNT character space
      IF (KKSURF .GT. 0) THEN
         DO 20 K = 1, KKSURF
            TPNAM(IPSURF+K-1) = KNAM(KKGAS+K)
   20    CONTINUE
      ENDIF
      IF (KKBULK .GT. 0) THEN
         DO 30 K = 1, KKBULK
            TPNAM(IPBULK+K-1)        = KNAM(KKGAS+KKSURF+K)
            TPNAM(IPBULK+KKBULK+K-1) = KNAM(KKGAS+KKSURF+K)
   30    CONTINUE
      ENDIF
      IF (LSDEN .AND. NNSURF.GT.0) THEN
         DO 40 N = 1, NNSURF
            TPNAM(IPSDEN+N-1) = PNAM(N+1)
   40    CONTINUE
      ENDIF
      DO 50 K = 1, KKGAS
         TPNAM(IPGAS + K + NYS - 1) = KNAM(K)
   50 CONTINUE
      TPNAM(IPGAS - 1 + NU) = 'VEL (AXIAL)'
      TPNAM(IPGAS - 1 + NV) = 'VEL (RADIAL)'
      TPNAM(IPGAS - 1 + NW) = 'VEL (CIRC)'
      TPNAM(IPGAS - 1 + NT) = 'TEMPERATURE'
      TPNAM(IPGAS - 1 + NL) = 'LAMBDA'
      ITPNAM = KKSURF + 2*KKBULK + NSSOLV + NYS + KKGAS
C
      IF (LFLAME .AND..NOT. (RSTCNT .OR. LRSTRT)) LFLAME = .FALSE.
      IF (LFLAME) THEN
         ITPNAM = ITPNAM + 1
         TPNAM(IPGAS - 1 + NUL) = 'UINLET'
         NATJ = NATJF
      ELSE
         NATJ = NATJF - 1
      ENDIF
C
      JJ = NPTS
      NBEQ = KKSURF + 2*KKBULK + NSSOLV
      NEQ = NATJ*JJ + NBEQ
C
      DO 110 K = 1, KKGAS
         YINJ(K) = 0.0
         HINJ(K) = 0.0
  110 CONTINUE
C
      DO 115 K = 1, KKBULK
         DEP(K) = 0.0
  115 CONTINUE
C
      IF (FINJ .GT. 0.0) THEN
         CALL CKXTY (BUFFER, ICKWRK, RCKWRK, YINJ)
         IF (TINJ .GT. 0.0) CALL CKHMS (TINJ, ICKWRK, RCKWRK, HINJ)
      ENDIF
C
C     SET THE VARIABLES ON WHICH TO ADAPT
C
C*****REFINE: Use the axial Velocity in grid refinement
      ACTIVE(NU) = .TRUE.
C*****END REFINE: Use the axial Velocity in grid refinement
C*****REFINE: Don't the axial Velocity in grid refinement
C      ACTIVE(NU) = .FALSE.
C*****END REFINE: Don't the axial Velocity in grid refinement
C
      ACTIVE(NV) = .TRUE.
      ACTIVE(NW) = .TRUE.
      ACTIVE(NT) =  LENRGY .AND. .NOT. LTOFF
      DO 120 K = 1, KKGAS
         ACTIVE(NYS+K) = .TRUE.
  120 CONTINUE
      IF (LFLAME) ACTIVE(NUL) = .FALSE.
C
      IF (LRSTRT) THEN
         IF (.NOT. RSTCNT) THEN
C           this is a restart
C
            NATJ = NATJF - 1
            CALL SPREAD (LRIN, LOUT, MFILE, NATJ, JJ, KKTOT, KKSURF,
     1                   KKBULK, NSSOLV, II, IISUR, X, S, LOLDSP,
     2                   KERR)
            IF (KERR) RETURN
C
         ELSEIF (NATJ .EQ. NATJF) THEN
            NEQ = NATJ*JJ + NBEQ
            CALL TWCOPY (NEQ, S, SN)
            NATJ = NATJF - 1
            DO 124 J = 1, JJ
               NKJ = IPGAS-1 + (J-1)*NATJF + 1
               NKJ2 = IPGAS-1 + (J-1)*NATJ + 1
               CALL TWCOPY (NATJ, S(NKJ2), SN(NKJ))
124         CONTINUE
            NEQ = NATJ*JJ + NBEQ
         ENDIF
C
         CALL SPRSTR (KKGAS, NTOT, NATJ, JJ, LOUT, LUMESH, LUSTGV,
     1                NREAC, NINTM, NPROD, REAC, XINTM, PROD, TINFTY,
     2                TDISK, AOMEGA, ASWIRL, P, KR, KI, KP, NTEMP,
     3                XGIVEN, TGIVEN, XSTR, XCEN, XEND, WMIX, ICKWRK,
     4                RCKWRK, Y, COND, EPS, X, S(IPGAS), IMCWRK,
     5                RMCWRK, RHOINF, ASPREAD, KERR)
         IF (KERR) RETURN
C
         IF (LREGRD .AND. (RSTCNT.OR.LRSTRT) .AND. JJ.GT.JJREGD) THEN
            WRITE (LOUT, *) ' REGRIDDING TO ', JJREGD, ' POINTS'
C
            CALL REGRID (SN(IPGAS), S(IPGAS), NATJ, JJREGD, JJ,
     $                   VIS(1), X(1), COND, PCTADP, RATGTC, NT,
     1                   KERR)
            IF (KERR) RETURN
C
            JJ = JJREGD
            DO 280 J = 2, JJ
               X(J) = VIS(J)
               NKJ = IPGAS-1 + (J-1)*NATJ + 1
               CALL TWCOPY (NL, SN(NKJ), S(NKJ))
               CALL TWCOPY (KKGAS, SN(NKJ+NYS), S(NKJ+NYS))
  280       CONTINUE
C
            IF (SPOS .GE. 0.0) THEN
               IF (KKSURF .GT. 0) CALL SPPOS (S(IPSURF), KKSURF, SPOS)
               IF (KKBULK .GT. 0) CALL SPPOS (S(IPBULK), KKBULK, SPOS)
               IF (LSDEN .AND. NNSURF.GT.0)
     1                            CALL SPPOS (S(IPSDEN), NNSURF, SPOS)
C
               DO 330 J = 1, JJ
                  NKJ = IPGAS-1 + (J-1)*NATJ + NYS + 1
                  CALL SPPOS (S(NKJ), KKGAS, SPOS)
  330          CONTINUE
            ENDIF
C
            NTEMP = JJ
            DO 340 N = 1, NTEMP
               XGIVEN(N) = X(N)
               TGIVEN(N) = S(IPGAS-1 + (N-1)*NATJ + NT)
  340       CONTINUE
         ENDIF
C
         IF (LFLAME) THEN
            IF (FLMW .LT. 0.0) THEN
              IFIX = 0
              DO 345 J = JJ-1, 2, -1
                 IF (IFIX .EQ. 1) GOTO 345
                 NKT1 = IPGAS-1 + (J-2)*NATJ + NT
                 NKT2 = IPGAS-1 + J*NATJ + NT
                 NKT = IPGAS-1 + (J-1)*NATJ + NT
                 IF ((S(NKT1).GT.FLMT)
     1                  .AND. (S(NKT2).LT.FLMT)) THEN
                   IFIX = 1
                   FXOLD = X(J)
                   FTOLD = S(NKT)
                   JFIXT = J
                 ENDIF
345           CONTINUE
              FLMT = FTOLD
              DO 356 J = 2, JJ
                 NKJ = IPGAS-1 + (J-1)*NATJ + 1
                 CALL TWCOPY (NATJ, SN(NKJ), S(NKJ))
C
                 NKJ2 = IPGAS-1 + (J-1)*NATJF + 1
                 CALL TWCOPY (NL, S(NKJ2), SN(NKJ))
                 CALL TWCOPY (KKGAS, S(NKJ2+NYS), SN(NKJ+NYS))
C
                 NKJ2 = IPGAS-1 + (J-1)*NATJF + NUL
                 S(NKJ2) = - UINF
356           CONTINUE
              S(IPGAS-1 + NUL) = - UINF
            ELSE
              CALL SPFLAM (FLMX, FLMT, FLMW, SN(IPGAS), S(IPGAS), NATJ,
     1                     JJNEW, JJ, VIS, X, NT, NTOT, JFIXT, IFLERR)
              IF (IFLERR.EQ.1) THEN
                 WRITE (LOUT, '(/1X, A/)')
     1            'FATAL ERROR IN FLAME RESTART: NO FLAME FOUND '
                 RETURN
              ELSEIF (IFLERR.GE.2) THEN
                 WRITE (LOUT,'(/1X, A/)')
     1            'FATAL ERROR IN FLAME RESTART: NMAX TOO SMALL'
                 RETURN
              ENDIF
              JJ = JJNEW
              DO 150 J = 2, JJ
                 X(J) = VIS(J)
                 NKJ = IPGAS-1 + (J-1)*NATJ + 1
                 NKJ2 = IPGAS-1 + (J-1)*NATJF + 1
                 CALL TWCOPY (NL, S(NKJ2), SN(NKJ))
                 CALL TWCOPY (KKGAS, S(NKJ2+NYS), SN(NKJ+NYS))
C
                 NKJ2 = IPGAS-1 + (J-1)*NATJF + NUL
                 S(NKJ2) = - UINF
150           CONTINUE
              S(IPGAS-1 + NUL) = - UINF
            ENDIF
            NATJ = NATJF
            XFIXT = X(JFIXT)
         ENDIF
      ELSE
C
        NATJ = NATJF - 1
C
        IF (.NOT. LNONR) WRITE (LOUT, '(/A,/A/)')
     1  ' ///////////////////////////////////',
     2  ' SPINDR: SOLVING THE NONREACTING PROBLEM'
C
         NSHT = 5
         CALL = 1
         CALL TWSTRT (CALL, LOUT, NTOT, CNTRLS, LVALUE, IVALUE,
     1                RVALUE, ERROR)
         CALL SPNRDR (LOUT, KKGAS, JJ, NSHT, P, WT, REAC, DT, X,
     1                XCEN, WMIX, LENRGY, LUINF, NTEMP, XGIVEN,
     2                TGIVEN, TDISK, TINFTY, AOMEGA, ASWIRL, RHOINF,
     3                UINF, SN, S, VIS, COND, ICKWRK, RCKWRK, IMCWRK,
     4                RMCWRK, F, FN, IASIZE, A, XMF, XSTR, XEND,
     5                LUMESH, ABOVE, ACTIVE, MARK, BELOW, BUFFER,
     6                NTOT, LENITW, ITWWRK, LENRTW, RTWWRK, NSHT,
     7                TPNAM(IPGAS), IP, QDOT, XSRC, WSRC, LNONR,
     8                NEMSG, TEMISG, EMISG, ASPREAD, LCOMP)
C
         IF (LNONR) THEN
            WRITE (LOUT, '(/A,/A/)')
     1      ' ////////////////////////////////////',
     2      ' SPINDR: SKIPPING THE NONREACTING PROBLEM'
         ELSE
            WRITE (LOUT, '(/A,/A/)')
     1      ' ///////////////////////////////////',
     2      ' SPINDR: FINISHED SOLVING THE NONREACTING PROBLEM'
         ENDIF
C
C        SET THE STARTING PROFILES
C        (in the following, COND is used as length njj scratch space)
C
         CALL SPSTRT (KKGAS, NTOT, NATJ, JJ, NSHT, LOUT, LUMESH, NREAC,
     1                NINTM, NPROD, REAC, XINTM, PROD, TINFTY, TDISK,
     2                AOMEGA, ASWIRL, P, KR, KI, KP, NTEMP, XGIVEN,
     3                TGIVEN, XSTR, XCEN, XEND, WMIX, ICKWRK, RCKWRK,
     4                YAV, COND, EPS, X, SN, S(IPGAS), IMCWRK,
     5                RMCWRK, RHOINF, ASPREAD, JDEX, KDEX, KFIRST,
     6                KLAST, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     7                NLBULK)
C
         IF (KKSURF .GT. 0) CALL TWCOPY (KKSURF, Z, S(IPSURF))
         IF (KKBULK .GT. 0) CALL TWCOPY (KKBULK, BUL, S(IPBULK))
         IF (KKBULK .GT. 0) CALL TWCOPY (KKBULK, DEP, S(IPBULK+KKBULK))
         IF (NSSOLV .GT. 0) CALL TWCOPY (NSSOLV, SDEN(2), S(IPSDEN))
C
      ENDIF
C
C///  DECIDE HOW MANY TIMES TO CALL TWOPNT.
C
      ENERGY = LENRGY
      IF (ENERGY) THEN
         CALLS = 3
      ELSE
         CALLS = 2
      ENDIF
C
C     set the solution bounds
      CALL SPSETB (NATJ, NTOT, KKGAS, KKSURF, KKBULK, NNSURF, SFLR,
     1             ABOVE(IPGAS), BELOW(IPGAS), ABOVE(IPSURF),
     2             BELOW(IPSURF), IPSURF, IPBULK, IPSDEN, IPGAS, LSDEN,
     3             LFLAME)
C
C///  TOP OF THE LOOP OVER CALLS TO TWOPNT.
C
      DO 2000 CALL = N1CALL, CALLS
C
         IF (CALL .EQ. 2 .AND. LRSTRT .AND. ENERGY) GO TO 2000
         IF (CALL .EQ. 3) WRITE (LOUT, '(/A,/A)')
     1   ' ///////////////////////////////////',
     2   ' SPINDR: INCLUDING THE ENERGY EQUATION'
C
         IF (CALL .EQ. 3) THEN
            LENRGY = .TRUE.
         ELSE
            LENRGY = .FALSE.
         ENDIF
C
         ACTIVE(NT) = LENRGY .AND. .NOT. LTOFF
         SIGNAL = ' '
C
  350    CONTINUE
C
         LGRA = KKSURF + 2*KKBULK + NSSOLV
C
         CALL TWSTRT (CALL, LOUT, NTOT, CNTRLS, LVALUE, IVALUE,
     1                RVALUE, ERROR)
         CALL TWOPNT (ERROR, LOUT, VERSIO, ABOVE, ACTIVE, BELOW, BUFFER,
     1                NATJ, CONDIT, LGRA, 0, LENITW, ITWWRK, MARK,
     3                TPNAM, ITPNAM, NTOT, JJ, REPORT, LENRTW,
     4                RTWWRK, SIGNAL, DT, FTIME, S, X)
C
         NBEQ = KKSURF + 2*KKBULK + NSSOLV
         NEQ = NATJ*JJ + NBEQ
C
         IF (ERROR) THEN
            RETURN
C
         ELSEIF (SIGNAL .EQ. 'RESIDUAL') THEN
C
C*****OPTIONS - LVARTP always positive
            LVARTP = .TRUE.
C*****END OPTIONS - LVARTP always positive
C*****OPTIONS - LVARTP positive only on jacobian reevaluations
C            LVARTP = .FALSE.
C*****END OPTIONS - LVARTP positive only on jacobian reevaluations
C
            LKMAX = .FALSE.
            IRATE = 1
            CALL SPFUN (LOUT, JJ, KFIRST, KLAST, LENRGY, LTDIF, LVCOR,
     1                  LVARTP, FTIME, LMULTI, WT, BUFFER(IPSDEN), DT,
     2                  NTEMP, XGIVEN, TGIVEN, EPS, P, TDISK, TINFTY,
     3                  AOMEGA, ASWIRL, RHOINF, UINF, LUINF, LCHEM, X,
     4                  BUFFER(IPSURF), SN(IPGAS), ACT, BUFFER(IPGAS),
     5                  F, SN(IPSURF), FZ, ICKWRK, RCKWRK, IMCWRK,
     6                  RMCWRK, ISKWRK, RSKWRK, SDOT, GFAC, SFAC,
     7                  SITDOT, VIS, COND, D, DKJ, DTCOEF, XAV, YAV,
     8                  YV, CP, H, XMF, XMFP, LKMAX, KDEX,
     9                  JDEX, QDOT, XSRC, WSRC, TWAL, EMIS, POWR,
     *                  TRAD, ERAD, VFAC, LRADB, FINJ, XINJ, WINJ, TINJ,
     1                  YINJ, HINJ, LCNDCT, CNDFAC, CNDT, NEMSG, TEMISG,
     2                  EMISG,IETCH, BUFFER(IPBULK), BULIN, SN(IPBULK),
     3                  SN(IPSDEN), FB, FS, NIICON, SDEN0, KKPHAS,
     4                  LSDEN, SDEN, LCOMP, USAVE, RSAVE, RKFT, RKRT,
     5                  IRATE, LCALSP, SP(IPSURF), SP(IPBULK),
     6                  SP(IPSDEN), SP(IPGAS),RU,BUFFER(IPBULK+KKBULK),
     7                  SP(IPBULK+KKBULK), FDEP, LFLAME, FLMT, JFIXT,
     8                  ASPREAD, LGDOT, GDOT)
C
            IF (KKSURF .GT. 0) CALL TWCOPY (KKSURF, FZ, BUFFER(IPSURF))
            IF (KKBULK .GT. 0) CALL TWCOPY (KKBULK, FB, BUFFER(IPBULK))
            IF (KKBULK .GT. 0) CALL TWCOPY (KKBULK, FDEP,
     1                                      BUFFER(IPBULK+KKBULK))
            IF (NSSOLV .GT. 0) CALL TWCOPY (NSSOLV, FS, BUFFER(IPSDEN))
            CALL TWCOPY (NATJ*JJ, F, BUFFER(IPGAS))
C
         ELSEIF (SIGNAL .EQ. 'PREPARE') THEN
C
            FIRST = .TRUE.
            FLAG = .FALSE.
            IRATE = 2
  375       CONTINUE
C
            CALL TWPREP
     +        (ERROR, LOUT,
     +         A, IASIZE, BUFFER, NATJ, CONDIT, LGRA, 0, IP, JJ, FLAG)
C
            IF (ERROR) THEN
               WRITE (LOUT, '(/1X, A /)')
     +            'FATAL ERROR, COULD NOT PREPARE JACOBIAN MATRIX.'
               RETURN
            ENDIF
C
            IF (FLAG) THEN
               CALL SPFUN (LOUT, JJ, KFIRST, KLAST, LENRGY, LTDIF,
     1              LVCOR, FIRST, FTIME, LMULTI, WT, BUFFER(IPSDEN),
     2              DT, NTEMP, XGIVEN, TGIVEN, EPS, P, TDISK, TINFTY,
     3              AOMEGA, ASWIRL, RHOINF, UINF, LUINF, LCHEM, X,
     4              BUFFER(IPSURF), SN(IPGAS), ACT, BUFFER(IPGAS),
     5              F(IPGAS), SN(IPSURF), F(IPSURF), ICKWRK, RCKWRK,
     6              IMCWRK, RMCWRK, ISKWRK, RSKWRK, SDOT, GFAC, SFAC,
     7              SITDOT, VIS, COND, D, DKJ, DTCOEF, XAV, YAV, YV,
     8              CP, H, XMF, XMFP, FIRST, KDEX, JDEX, QDOT,
     9              XSRC, WSRC, TWAL, EMIS, POWR, TRAD, ERAD, VFAC,
     *              LRADB, FINJ, XINJ, WINJ, TINJ, YINJ, HINJ, LCNDCT,
     1              CNDFAC, CNDT, NEMSG, TEMISG, EMISG, IETCH,
     2              BUFFER(IPBULK), BULIN, SN(IPBULK), SN(IPSDEN),
     3              F(IPBULK), F(IPSDEN), NIICON, SDEN0, KKPHAS,
     4              LSDEN, SDEN, LCOMP, USAVE, RSAVE, RKFT, RKRT,
     5              IRATE, LCALSP, SP(IPSURF), SP(IPBULK), SP(IPSDEN),
     6              SP(IPGAS),RU, BUFFER(IPBULK+KKBULK),
     7              SP(IPBULK+KKBULK), F(IPBULK+KKBULK), LFLAME, FLMT,
     8              JFIXT, ASPREAD, LGDOT, GDOT)
C
               CALL TWCOPY (NEQ, F, BUFFER)
               IRATE = 3
               FIRST = .FALSE.
               GO TO 375
            ENDIF
C
         ELSEIF (SIGNAL .EQ. 'SOLVE') THEN
C
            CALL TWSOLV (ERROR, LOUT, A, IASIZE, BUFFER, NATJ, LGRA,
     1                   0, IP, JJ)
            IF (ERROR) THEN
               WRITE (LOUT, '(/1X, A /)')
     +            'FATAL ERROR, COULD NOT SOLVE LINEAR EQUATIONS.'
               RETURN
            ENDIF
C
         ELSEIF (SIGNAL .EQ. 'RETAIN') THEN
C
            CALL TWCOPY (NEQ, BUFFER, SN)
C
            IF (SPOS .GE. 0.0) THEN
               IF (KKSURF .GT. 0) CALL SPPOS (SN(IPSURF), KKSURF, SPOS)
               IF (KKBULK .GT. 0) CALL SPPOS (SN(IPBULK), KKBULK, SPOS)
               IF (LSDEN .AND. NNSURF.GT.0)
     1                            CALL SPPOS (SN(IPSDEN), NNSURF, SPOS)
               DO 430 J = 1, JJ
                  NKJ = IPGAS-1 + (J-1)*NATJ + NYS + 1
                  CALL SPPOS (SN(NKJ), KKGAS, SPOS)
  430          CONTINUE
            ENDIF
C
         ELSEIF (SIGNAL .EQ. 'SHOW') THEN
            LVARTP = .TRUE.
            LKMAX = LFIRST
            IRATE = 1
            CALL SPFUN (LOUT, JJ, KFIRST, KLAST, LENRGY, LTDIF, LVCOR,
     1                  LVARTP, FTIME, LMULTI, WT, BUFFER(IPSDEN), DT,
     2                  NTEMP, XGIVEN, TGIVEN, EPS, P, TDISK, TINFTY,
     3                  AOMEGA, ASWIRL, RHOINF, UINF, LUINF, LCHEM, X,
     4                  BUFFER(IPSURF), SN(IPGAS), ACT, BUFFER(IPGAS),
     5                  F(IPGAS), SN(IPSURF), F(IPSURF), ICKWRK, RCKWRK,
     6                  IMCWRK,RMCWRK, ISKWRK, RSKWRK, SDOT, GFAC, SFAC,
     7                  SITDOT, VIS, COND, D, DKJ, DTCOEF, XAV, YAV,
     8                  YV, CP, H, XMF, XMFP, LKMAX, KDEX,
     9                  JDEX, QDOT, XSRC, WSRC, TWAL, EMIS, POWR, TRAD,
     *                  ERAD, VFAC, LRADB, FINJ, XINJ, WINJ, TINJ, YINJ,
     1                  HINJ, LCNDCT, CNDFAC, CNDT, NEMSG, TEMISG,
     2                  EMISG,IETCH, BUFFER(IPBULK), BULIN, SN(IPBULK),
     3                  SN(IPSDEN), F(IPBULK), F(IPSDEN), NIICON, SDEN0,
     4                  KKPHAS, LSDEN, SDEN, LCOMP, USAVE, RSAVE, RKFT,
     5                  RKRT, IRATE, LCALSP, SP(IPSURF), SP(IPBULK),
     6                  SP(IPSDEN),SP(IPGAS),RU,BUFFER(IPBULK+KKBULK),
     7                  SP(IPBULK+KKBULK), F(IPBULK+KKBULK), LFLAME,
     8                  FLMT, JFIXT, ASPREAD, LGDOT, GDOT)
C
            LFIRST = .FALSE.
C
            CALL SPPRNT (LOUT, JJ, KFIRST, KLAST, P, X, BUFFER(IPSURF),
     1                   BUFFER(IPGAS), BUFFER(IPBULK), BUFFER(IPSDEN),
     $                   ACT, Y, XMF, FLUX, ROP, SDOT,
     2                   SITDOT, SDOTI, SITDTI, KNAM, PNAM, SDEN, DEN,
     3                   WT, H, ICKWRK, RCKWRK, ISKWRK, RSKWRK, CSKWRK,
     4                   RMCWRK, LPRNT, ENAM, KNCF, AWT, CDOT, DDOT,
     5                   WDOT, YV, REAFLX, RADFLX, LCOMP, LPUTIL, LVCOR,
     6                   EPS, LREORD, LRADB, COND, TINFTY, VFAC, ERAD,
     7                   TRAD, TWAL, EMIS, POWR, LSDEN, SFAC, LCNDCT,
     8                   CNDFAC, CNDT, LGDOT, GDOT)
C
         ELSEIF (SIGNAL .EQ. 'SAVE') THEN
            REWIND LRCRVR
            CALL SPSAVE (ICKWRK, RCKWRK, CCKWRK,
     1                   IMCWRK, RMCWRK, ISKWRK, RSKWRK,
     2                   CSKWRK, LOUT, LRCRVR, IVERSN, VERNUM)
            NATJ = NATJF - 1
            NEQ = NATJ * JJ + NBEQ
            TIME = 0.0
            CALL SPSOLN (LRCRVR, STEADY, JJ,     NEQ,    NATJ,
     1                   KKGAS,  KKSURF, KKBULK, NSSOLV, IMULTI,
     2                   ITDIF,  IVCOR,  ICOMP,  ICHEM,  IGDOT,
     3                   IREORD, ICNDCT, IUINF,  IENRGY, IRADB,
     4                   P,      UINF,   GFAC,   SFAC,   CNDFAC,
     5                   POWR,   CNDT,   XINJ,   WINJ,   FINJ,
     6                   TINJ,   QDOT,   XSRC,   WSRC,   TWAL,
     7                   EMIS,   ERAD,   VFAC,   TRAD,   NEMSG,
     8                   EPS,    GDOT,   YINJ,   TEMISG, EMISG,
     9                   TIME,   LFLAME, IPGAS,  NATJF,  NBEQ,
     +                   X,      BUFFER )
C
         ELSEIF (SIGNAL .EQ. 'UPDATE') THEN
            IF (.NOT. LENRGY) THEN
              DO 440 J = 1, JJ
                 BUFFER(IPGAS-1 + (J-1)*NATJ + NT) =
     1                SPTEMP (NTEMP, X(J), XGIVEN, TGIVEN)
  440         CONTINUE
            ENDIF
            IF (LFLAME) THEN
               DO 445 J = 1, JJ
                  IF (X(J) .EQ. XFIXT) JFIXT = J
445            CONTINUE
            ENDIF
C
         ELSEIF (REPORT.NE.' ' .AND. SIGNAL.EQ.' ') THEN
            IF (REPORT .EQ. 'NO SPACE') WRITE (LOUT, *)
     1      '     TWOPNT requires more mesh points, but NMAX too small'
C
         ENDIF
C
         IF (SIGNAL .NE. ' ') GO TO 350
C
 2000 CONTINUE
C
C              WRITE TO LROUT WHEN SOLUTION IS COMPLETE
C
      NATJ = NATJF - 1
      NEQ = NATJ*JJ + NBEQ
      TIME = 0.0
      CALL SPSOLN       (LROUT,  STEADY, JJ,     NEQ,    NATJ,
     1                   KKGAS,  KKSURF, KKBULK, NSSOLV, IMULTI,
     2                   ITDIF,  IVCOR,  ICOMP,  ICHEM,  IGDOT,
     3                   IREORD, ICNDCT, IUINF,  IENRGY, IRADB,
     4                   P,      UINF,   GFAC,   SFAC,   CNDFAC,
     5                   POWR,   CNDT,   XINJ,   WINJ,   FINJ,
     6                   TINJ,   QDOT,   XSRC,   WSRC,   TWAL,
     7                   EMIS,   ERAD,   VFAC,   TRAD,   NEMSG,
     8                   EPS,    GDOT,   YINJ,   TEMISG, EMISG,
     9                   TIME,   LFLAME, IPGAS,  NATJF,  NBEQ,
     +                   X,      S )
C
      IF (LASEN .OR. LHSEN) THEN
         IF (LFLAME) THEN
            NEQ = NATJ*JJ + NBEQ
            CALL TWCOPY (NEQ, S, BUFFER)
            NATJ = NATJF - 1
            DO 2015 J = 1, JJ
               DO 2010 N = 1, NATJ
                  NKJ = IPGAS-1 + (J-1)*NATJF + N
                  NKJ2 = IPGAS-1 + (J-1)*NATJ + N
                  S(NKJ2) = BUFFER(NKJ)
2010           CONTINUE
2015        CONTINUE
            NEQ = NATJ*JJ + NBEQ
         ENDIF
      ENDIF
C
      IF (LASEN) THEN
C
         CALL SPSENS (ISENSI, NEQ, LROUT, LRCRVR, LOUT, A,
     1                JJ, KFIRST, KLAST, LENRGY, LTDIF, LVCOR, LVARTP,
     2                FTIME, LMULTI, WT, SDEN, DT, NTEMP, XGIVEN,
     3                TGIVEN, EPS, P, TDISK, TINFTY, AOMEGA, ASWIRL,
     4                RHOINF, UINF, LUINF, LCHEM, X, S, ACT, F, SN, FN,
     5                ICKWRK, RCKWRK, IMCWRK, RMCWRK, ISKWRK, RSKWRK,
     6                VIS, COND, D, DKJ, DTCOEF, XAV, YAV, YV, CP, H,
     7                XMF, XMFP, IP, DRDAI, DRDCL, SDOT, GFAC, SFAC,
     8                SITDOT, RU, LKMAX, KDEX, JDEX, QDOT, XSRC,
     9                WSRC, TWAL, EMIS, POWR, TRAD, ERAD, VFAC, LRADB,
     *                A6, FINJ, XINJ, WINJ, TINJ, YINJ, HINJ, BUFFER,
     1                LCNDCT, CNDFAC, CNDT, NEMSG, TEMISG, EMISG,IETCH,
     2                BULIN, NIICON, SDEN0, KKPHAS, LSDEN, NSSOLV, LGRA,
     3                LCOMP, USAVE, RSAVE, RKFT, RKRT, LCALSP, SP,
     4                DEP, DEPP, FDEP, LFLAME, FLMT, JFIXT, ASPREAD,
     5                LGDOT, GDOT, KERR)
         IF (KERR) RETURN
C
C
         WRITE (LOUT, '(/A/)') ' RATE SENSITIVITY CALCULATION COMPLETE'
         LASEN = .FALSE.
C
      ENDIF
C
      IF (LHSEN) THEN
C
         CALL SPSENS (IHSENS, NEQ, LROUT, LRCRVR, LOUT, A,
     1                JJ, KFIRST, KLAST, LENRGY, LTDIF, LVCOR, LVARTP,
     2                FTIME, LMULTI, WT, SDEN, DT, NTEMP, XGIVEN,
     3                TGIVEN, EPS, P, TDISK, TINFTY, AOMEGA, ASWIRL,
     4                RHOINF, UINF, LUINF, LCHEM, X, S, ACT, F, SN, FN,
     5                ICKWRK, RCKWRK, IMCWRK, RMCWRK, ISKWRK, RSKWRK,
     6                VIS, COND, D, DKJ, DTCOEF, XAV, YAV, YV, CP, H,
     7                XMF, XMFP, IP, DRDAI, DRDCL, SDOT, GFAC, SFAC,
     8                SITDOT, RU, LKMAX, KDEX, JDEX, QDOT, XSRC,
     9                WSRC, TWAL, EMIS, POWR, TRAD, ERAD, VFAC, LRADB,
     *                A6, FINJ, XINJ, WINJ, TINJ, YINJ, HINJ, BUFFER,
     1                LCNDCT, CNDFAC, CNDT, NEMSG, TEMISG, EMISG,IETCH,
     2                BULIN, NIICON, SDEN0, KKPHAS, LSDEN, NSSOLV, LGRA,
     3                LCOMP, USAVE, RSAVE, RKFT, RKRT, LCALSP, SP,
     4                DEP, DEPP, FDEP, LFLAME, FLMT, JFIXT, ASPREAD,
     5                LGDOT, GDOT, KERR)
C
         WRITE (LOUT, '(/A/)') ' H SENSITIVITY CALCULATION COMPLETE'
         LHSEN = .FALSE.
C
      ENDIF
C
      NSOL = NSOL + 1
      CALL SPFPNT (LFLUX, NSOL, ENAM, MM, KNAM, KKGAS, FLUX,
     1             KNCF, WT, TFLUXA)
C
C     check for continuation
      IF (.NOT. LCNTUE) RETURN
C
      WRITE (LOUT, '(/////)')
      DO 2020 L = 1, 5
         WRITE (LOUT, *)
     1    ' ////////////////// CONTINUING TO NEW PROBLEM /////////////'
 2020 CONTINUE
      WRITE (LOUT, '(/////)')
C
      RSTCNT = .TRUE.
      LRSTRT = .TRUE.
      NPTS = JJ
      GO TO 100
C
C     end of SUBROUTINE SPINDR
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE SPLNWM (WMIX, XCEN, XNODE, XRE, XPD)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      IF (XNODE .LE. (XCEN-WMIX/2.)) THEN
         XPD = 0.0
      ELSE
         IF (XNODE .LT. (XCEN+WMIX/2.)) THEN
            XPD = (1.0/WMIX)*(XNODE-XCEN) + 0.5
         ELSE
            XPD = 1.0
         ENDIF
      ENDIF
C
      XRE = 1.0 - XPD
C
C     end of SUBROUTINE SPLNWM
      RETURN
      END
C
C---------------------------------------------------------------------
      SUBROUTINE SPNORM (S1, S2, X, NN)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION S1(NN), S2(NN)
C
      DO 100 N = 1, NN
         IF (S2(N) .NE. 0.0) THEN
           S1(N) = S1(N) * X / ABS(S2(N))
         ELSE
           S1(N) = X
         ENDIF
  100 CONTINUE
C
C     end of SUBROUTINE SPNORM
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE SPNRDR (LOUT, KKGAS, JJ, NATJ,  P, WT, REAC, DT, X,
     1                   XCEN, WMIX, LENRGY, LUINF, NTEMP, XGIVEN,
     2                   TGIVEN, TDISK, TINFTY, AOMEGA, ASWIRL, RHOINF,
     3                   UINF, SN, S, VIS, COND, ICKWRK, RCKWRK,
     4                   IMCWRK, RMCWRK, F, FN, IASIZE, A, EPS, XSTR,
     5                   XEND, LUMESH, ABOVE, ACTIVE, MARK, BELOW,
     6                   BUFFER, NMAX, LENITW, ITWWRK, LENRTW, RTWWRK,
     7                   ITPNAM, TPNAM, IP, QDOT, XSRC, WSRC, LNONR,
     8                   NEMSG, TEMISG, EMISG, ASPREAD, LCOMP)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
C
C     Integer arrays
      DIMENSION ICKWRK(*), IMCWRK(*), IP(*), ITWWRK(*)
C     Real arrays
      DIMENSION RCKWRK(*), RMCWRK(*), REAC(*), WT(*), EPS(*),
     1          X(*), VIS(*), COND(*), S(NATJ,*), SN(NATJ,JJ),
     2          F(*), XGIVEN(*), TGIVEN(*), FN(NATJ,*),
     3          TEMISG(*), EMISG(*), ABOVE(NATJ), BELOW(NATJ),
     4          BUFFER(*), RTWWRK(*), A(IASIZE)
C
      CHARACTER*(*) TPNAM(NATJ)
      CHARACTER VERSIO*80, REPORT*16, SIGNAL*16
C
      LOGICAL LUMESH, LENRGY, ERROR, REENTR, FTIME, FIRST,
     1        ACTIVE(*), MARK(*), LUINF, FLAG, LNONR, LCOMP
C
C*****precision > double
      DATA VERSIO/'DOUBLE PRECISION VERSION 3.22'/
C*****END precision > double
C*****precision > single
C      DATA VERSIO/'SINGLE PRECISION VERSION 3.22'/
C*****END precision > single
C
C     SET THE SOLUTION BOUNDS
C
C     axial velocity
      BELOW(NU) = -1.E30
      ABOVE(NU) =  1.E30
C
C     radial velocity
      BELOW(NV) = -1.E30
      ABOVE(NV) = 1.E30
C
C     circumferential velocity
      BELOW(NW) = -10.0
      ABOVE(NW) = 1000.0
C
C     temperature
      BELOW(NT) = 200.0
      ABOVE(NT) = 5000.0
C
C     lambda
      BELOW(NL) = -1.E30
      ABOVE(NL) =  1.E30
C
C     set the initial profiles
      IF (LUMESH) THEN
C         set uniform X mesh coordinates
C
         DX = (XEND-XSTR) / FLOAT(JJ-1)
         DO 100 J = 1,JJ
            X(J) = XSTR + DX*FLOAT(J-1)
  100    CONTINUE
      ENDIF
C
C     reference properties at inlet
      CALL CKXTY (REAC, ICKWRK, RCKWRK, EPS)
C
C     transport properties at infinity
      CALL CKRHOY (P, TINFTY, EPS, ICKWRK, RCKWRK, RHOINF)
      CALL MCAVIS (TINFTY, REAC, RMCWRK, XMUINF)
      XNUINF = XMUINF / RHOINF
C
C     INITIAL GUESSES FOR THE NONREACTING CASE
C
C     set temperature, F, G, H, and lambda profiles
      TSLOPE = - (TDISK - TINFTY) / (X(JJ) - X(1))
C
      IF (LUINF) THEN
C        user has specified an input axial velocity at X=L
C
         XJJ2 = X(JJ)**2
         IF (AOMEGA .EQ. 0.0) THEN
C           pure stagnation-point flow
C
            DO 200 J = 1, JJ
               S(NU,J) = -UINF * X(J) / X(JJ)
               IF (X(J) .LE. X(JJ)/3.0) THEN
                  S(NV,J) = 2.4 * UINF * X(J) / XJJ2
     1                       + ASPREAD * X(J) / X(JJ)
               ELSE
                  S(NV,J) = 1.2 * UINF * (X(JJ) - X(J)) / XJJ2
     1                       + ASPREAD * X(J) / X(JJ)
               ENDIF
               S(NW,J) = ASWIRL * X(J) / X(JJ)
               S(NL,J) = -1.5  * RHOINF * (ASWIRL + UINF / X(JJ))**2
               S(NT,J) = TSLOPE*(X(J) - X(1)) + TDISK
  200       CONTINUE
C
         ELSE
C           the disk is rotating; determine whether it is being
C           forced or starved.
C
C           calculate two Reynolds numbers
            REW = XJJ2 * AOMEGA / XNUINF
            REU = UINF * X(JJ) / XNUINF
            UIDEAL = 0.8838 * XNUINF * SQRT(REW + REU) / X(JJ)
C
            ASLOPE = (AOMEGA - ASWIRL) / (X(1) - X(JJ))
C
            IF (UINF .GT. UIDEAL) THEN
C              forced flow greater than required by the infinite
C              rotating disk alone.
C
               DO 250 J = 1, JJ
                  S(NU,J) = -UINF * X(J) / X(JJ)
                  IF (X(J) .LE. X(JJ)/5.0) THEN
                     S(NV,J) = AOMEGA * X(J) / X(JJ)
     1                       + ASPREAD * X(J) / X(JJ)
                  ELSE
                     S(NV,J) = AOMEGA * (X(JJ) - X(J)) / (4.0 * X(JJ))
     1                       + ASPREAD * X(J) / X(JJ)
                  ENDIF
C                  IF (X(J) .LE. X(JJ)/3.0) THEN
C                     S(NW,J) = AOMEGA * (1.0 - 3.0*X(J)/X(JJ))
C                  ELSE
C                     S(NW,J) = 0.0
C                  ENDIF
                  S(NW,J) = AOMEGA + ASLOPE * (X(J) - X(1))
                  S(NL,J) = -0.15  * RHOINF
     +                 * (AOMEGA + ASWIRL + UINF/X(JJ))**2
                  S(NT,J) = TSLOPE*(X(J) - X(1)) + TDISK
  250          CONTINUE
C
            ELSE
C              forced flow less than required by infinite disk alone.
C
               DO 300 J = 1, JJ
                  IF (X(J) .LE. X(JJ)/3.0) THEN
                     S(NU,J) = -2.4 * XNUINF * X(J) * SQRT(REW + REU) /
     1                          XJJ2
C                     S(NW,J) = AOMEGA * (1.0 - 3.0*X(J)/X(JJ))
                  ELSE
                     S(NU,J) = -UINF + (X(JJ)-X(J)) * ( 1.5*UINF/X(JJ) -
     1                       1.2 * XNUINF * SQRT(REW + REU) / XJJ2)
C                     S(NW,J) = 0.0
                  ENDIF
                  IF (X(J) .LE. X(JJ)/5.0) THEN
                     S(NV,J) = AOMEGA * X(J) / X(JJ)
     1                       + ASPREAD * X(J) / X(JJ)
                  ELSEIF (X(J) .LE. X(JJ)/2.0) THEN
                     S(NV,J) = 2.0 * AOMEGA * (0.5*X(JJ) - X(J)) /
     1                         (3.0 * X(JJ)) + ASPREAD * X(J) / X(JJ)
                  ELSEIF (X(J) .LE. 0.75*X(JJ)) THEN
                     S(NV,J) = -2.0 * AOMEGA * (X(J) - 0.5*X(JJ)) /
     1                         (10.0 * X(JJ)) + ASPREAD * X(J) / X(JJ)
                  ELSE
                     S(NV,J) = -2.0 * AOMEGA * ( X(JJ) - X(J)) /
     1                         (10.0 * X(JJ))
                  ENDIF
                  S(NW,J) = AOMEGA + ASLOPE * (X(J) - X(1))
                  S(NL,J) = 0.15  * RHOINF
     +                 * (AOMEGA + ASWIRL + UINF/X(JJ))**2
                  S(NT,J) = TSLOPE*X(J) + TDISK
  300          CONTINUE
            ENDIF
         ENDIF
C
      ELSE
C        LUINF was false; ideal rotating disk profiles.
C
         DO 400 J = 1, JJ
            S(NU,J) =  -SQRT(AOMEGA)*X(J)/X(JJ)
            IF (X(J)/X(JJ) .LT. 0.2) THEN
               S(NV,J) = 0.5 * AOMEGA * X(J)/(0.2*X(JJ))
            ELSE
               S(NV,J) = 5./8.*AOMEGA * (X(JJ)-X(J)) / X(JJ)
            ENDIF
            S(NW,J) = AOMEGA*(1.-X(J)/X(JJ))
            S(NL,J) = 0.0
            S(NT,J) = TSLOPE*X(J) + TDISK
  400    CONTINUE
      ENDIF
C
      IF (LNONR) THEN
         DO 425 J = 1, JJ
            S(NT, J) = SPTEMP(NTEMP, X(J), XGIVEN, TGIVEN)
  425    CONTINUE
C
      ELSE
C
         REENTR=.FALSE.
         FIRST = .TRUE.
C
C///  TOP OF THE LOOP OVER CALLS TO TWOPNT.
         SIGNAL = ' '
  500    CONTINUE
C
         CALL TWOPNT (ERROR, LOUT, VERSIO, ABOVE, ACTIVE, BELOW,
     1                BUFFER, NATJ, CONDIT, 0, 0, LENITW, ITWWRK,
     2                MARK, TPNAM, ITPNAM, NMAX, JJ, REPORT,
     3                LENRTW, RTWWRK, SIGNAL, DT, FTIME, S, X)
C
         IF (ERROR) THEN
            RETURN
C
         ELSEIF (SIGNAL .EQ. 'RESIDUAL') THEN
C
            CALL SPNRFN (LOUT, JJ, NATJ, LENRGY, FTIME, LUINF, P, WT,
     1                   EPS, DT, X, NTEMP, XGIVEN, TGIVEN, TDISK,
     2                   TINFTY, AOMEGA, ASWIRL, RHOINF, UINF, SN,
     3                   BUFFER, VIS, COND, ICKWRK, RCKWRK, IMCWRK,
     4                   RMCWRK, F, REAC, QDOT, XSRC, WSRC, FIRST,
     5                   NEMSG, TEMISG, EMISG, ASPREAD)
            FIRST = .FALSE.
            CALL TWCOPY (NATJ*JJ, F, BUFFER)
C
          ELSEIF (SIGNAL .EQ. 'PREPARE') THEN
C
             FLAG = .FALSE.
 1010        CONTINUE
C
            CALL TWPREP (ERROR, LOUT, A, IASIZE, BUFFER, NATJ,
     1                   CONDIT, 0, 0, IP, JJ, FLAG)
C
            IF (ERROR) THEN
               WRITE (LOUT, '(/1X, A /)')
     +            'FATAL ERROR, COULD NOT PREPARE JACOBIAN MATRIX.'
               RETURN
            ENDIF
C
            IF (FLAG) THEN
               CALL SPNRFN (LOUT, JJ, NATJ, LENRGY, FTIME, LUINF,
     1                      P, WT, EPS, DT, X, NTEMP, XGIVEN,
     2                      TGIVEN, TDISK, TINFTY, AOMEGA, ASWIRL,
     3                      RHOINF, UINF, SN, BUFFER, VIS, COND,
     4                      ICKWRK, RCKWRK, IMCWRK, RMCWRK, F, REAC,
     5                      QDOT, XSRC, WSRC, FIRST, NEMSG,
     6                      TEMISG, EMISG, ASPREAD)
C
               CALL TWCOPY (NATJ*JJ, F, BUFFER)
               FIRST = .FALSE.
               GO TO 1010
            ENDIF
C
         ELSEIF (SIGNAL .EQ. 'SOLVE') THEN
C
            CALL TWSOLV (ERROR, LOUT, A, IASIZE, BUFFER, NATJ,
     1                   0, 0, IP, JJ)
C
            IF (ERROR) THEN
               WRITE (LOUT, '(/1X, A /)')
     +            'FATAL ERROR, COULD NOT SOLVE LINEAR EQUATIONS.'
               RETURN
            ENDIF
C
         ELSEIF (SIGNAL .EQ. 'RETAIN') THEN
            CALL TWCOPY (NATJ*JJ, BUFFER, SN)
C
         ELSEIF (SIGNAL .EQ. 'SHOW') THEN
            CALL SPNRPR (LOUT, JJ, NATJ, X, BUFFER, ICKWRK, RCKWRK)
C
         ELSEIF (REPORT.NE.' ' .AND. SIGNAL.EQ.' ') THEN
            IF (REPORT .EQ. 'NO SPACE') WRITE (LOUT, *)
     1   '     TWOPNT requires more mesh points, but NMAX too small'
         ENDIF
C
         IF (SIGNAL .NE. ' ') GO TO 500
C
C        END OF CODE FOR CALLING TWOPNT
      ENDIF
C
C     copy S into SN
      CALL TWCOPY (NATJ*JJ, S, SN)
C
      IF (.NOT. LENRGY) RETURN
C
C     copy the nonreacting temperatures into TGIVEN
      NTEMP = JJ
      DO 700 J = 1, JJ
         XGIVEN(J) = X(J)
         TGIVEN(J) = S(NT,J)
700   CONTINUE
C
C     end of SUBROUTINE SPNRDR
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE SPNRFN (LOUT, JJ, NATJ, LENRGY, FTIME, LUINF, P, WT,
     1                   EPS, DT, X, NTEMP, XGIVEN, TGIVEN, TDISK,
     2                   TINFTY, AOMEGA, ASWIRL, RHOINF, UINF, SN, S,
     3                   VIS, COND, ICKWRK, RCKWRK, IMCWRK, RMCWRK, F,
     4                   REAC, QDOT, XSRC, WSRC, FIRST, NEMSG,
     5                   TEMISG, EMISG, ASPREAD)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (PI=3.141592654D0, SIGMA = 5.67D-5, ZERO=0.0D0,
     $           ONE=1.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (PI=3.141592654, SIGMA = 5.67E-5, ZERO=0.0,
C     $           ONE=1.0)
C*****END precision > single
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
C     Integer arrays
      DIMENSION ICKWRK(*), IMCWRK(*)
C     Real arrays
      DIMENSION WT(*), EPS(*), X(*), S(NATJ,*), SN(NATJ,*),
     1          XGIVEN(*), TGIVEN(*), VIS(*), COND(*),
     2          RCKWRK(*), RMCWRK(*), F(NATJ,*), REAC(*),
     3          TEMISG(*), EMISG(*)
C
      LOGICAL FTIME, LENRGY, LUINF, FIRST
C
      LOGICAL LHEATFLX
      SAVE    LHEATFLX
C
C*****OPTIONS: Use Heat Flux BC at top of domain
C      DATA    LHEATFLX /.TRUE./
C*****END OPTIONS: Use Heat Flux BC at top of domain
C*****OPTIONS: Don't Use Heat Flux BC at top of domain
      DATA    LHEATFLX /.FALSE./
C*****END OPTIONS: Don't Use Heat Flux BC at top of domain
C
C
      IF (FIRST) THEN
         DO 100 J = 1, JJ-1
            TAV = 0.5*(S(NT,J+1)+S(NT,J))
            CALL MCAVIS (TAV, REAC, RMCWRK, VIS(J))
            CALL MCACON (TAV, REAC, RMCWRK, COND(J))
  100    CONTINUE
      ENDIF
C
C     calculate a characteristic velocity for satisfaction of
C     Dirichelet conditions
CC      IF (S(NU,1) .LT. ZERO) THEN
CC        VELOC = - S(NU,1)
CC      ELSEIF (ABS(AOMEGA) .GT. ZERO) THEN
      IF (ABS(AOMEGA) .GT. ZERO) THEN
        VELOC = ABS(X(JJ)-X(1)) * ABS(AOMEGA)
      ELSEIF (UINF .GT. ZERO) THEN
        VELOC = ABS(UINF)
      ELSE
        VELOC = ONE
      ENDIF
C
C     mixture density (RHOP) at J = 1+1/2
      TAV = 0.5*(S(NT,1)+S(NT,2))
      CALL CKRHOY (P, TAV, EPS, ICKWRK, RCKWRK, RHOP)
C
C     mixture density at the disk surface
      CALL CKRHOY (P, S(NT,1), EPS, ICKWRK, RCKWRK, RHOJ)
C
C     mixture density at node J+1 = 2
      CALL CKRHOY (P, S(NT,2), EPS, ICKWRK, RCKWRK, RHOJP1)
C
C     specific heat of mixture at J
      CALL CKCPBS (S(NT,1), EPS, ICKWRK, RCKWRK, CPB)
C
C     DISK SURFACE
C
C     zero slip conditions on the radial and circumferential velocity
C     components
      DXINV = 2.0 / (X(2) - X(1))
      F(NV,1) = RHOP * VELOC * DXINV *  S(NV,1)
      F(NW,1) = RHOP * VELOC * DXINV * (S(NW,1) - AOMEGA)
      F(NU,1) = RHOP * VELOC * DXINV *
     $                    (S(NU,1) * DXINV + 2.0 * S(NV,1))
C
C     temperature of disk
      F(NT,1) = RHOP * CPB*VELOC * DXINV * (S(NT,1) - TDISK)
C
C     zero-derivative condition for lambda
      F(NL,1) = VELOC * DXINV * (S(NL,1) - S(NL,2))
C
C-----INTERIOR MESH POINTS
C
      DO 1000 J = 2, JJ-1
C
         JP1 = J + 1
         JM1 = J - 1
         TAV = 0.5*(S(NT,J)+S(NT,JP1))
         RHOM = RHOP
         CALL CKRHOY (P, TAV, EPS, ICKWRK, RCKWRK, RHOP)
C
C        densities at nodes J and J+1
         RHOJ   = RHOJP1
         CALL CKRHOY (P, S(NT,JP1), EPS, ICKWRK, RCKWRK, RHOJP1)
C
C        form the mesh differences
         DXP =        (X(JP1) - X(J)  )
         DXM =        (X(J)   - X(JM1))
         DXAV = 0.5 * (X(JP1) - X(JM1))
         DXPM =       (X(JP1) - X(JM1))
         DXINV = ONE / DXAV
C
C        form the coefficients for central differences
         CENDFM = - DXP / (DXM*DXPM)
         CENDFC =   (DXP-DXM) / (DXP*DXM)
         CENDFP =   DXM / (DXP*DXPM)
C
C        mixture continuity equation (changed to conservation form)
         F(NU,J) = (RHOP*S(NU,J) - RHOM*S(NU,JM1)) * DXINV
     $            + 2.0*RHOJ*S(NV,J)
C
C        mixture radial momentum equation
         F(NV,J) = (VIS(JM1)*(S(NV,J)   - S(NV,JM1))/DXM -
     $              VIS(J)  *(S(NV,JP1) - S(NV,J)  )/DXP) *DXINV
     $           + RHOP*S(NU,J)*(S(NV,JP1) - S(NV,J)) *DXINV
     $           + RHOJ * (S(NV,J)**2 - S(NW,J)**2)
     $           + S(NL,J)
C
C        mixture circumferential momentum equation
         F(NW,J) = (VIS(JM1)*(S(NW,J)   - S(NW,JM1))/DXM -
     $              VIS(J)  *(S(NW,JP1) - S(NW,J)  )/DXP  )*DXINV
     $           + RHOP*S(NU,J)*DXINV*(S(NW,JP1) - S(NW,J))
     $           + 2.0*RHOJ*S(NV,J)*S(NW,J)
C
C        energy equation
         IF (LENRGY) THEN
C
C           determine the specific heat of the mixture at J
            CALL CKCPBS (S(NT,J), EPS, ICKWRK, RCKWRK, CPB)
C
C           calculate energy source term
            SRCJ = ZERO
            IF (QDOT .GT. ZERO) THEN
               XCUT = 2.*WSRC
               IF (ABS(X(J) - XSRC).LE.XCUT) THEN
                  Q0 = 1.283153 * QDOT * SQRT (3 / PI) / WSRC
                  SRCJ = Q0 * EXP(-3.*(X(J)-XSRC)**2./WSRC**2.)
               ENDIF
            ENDIF
C
C           calculate radiation loss term
            RADJ = ZERO
            IF (NEMSG .GT. 0) THEN
               EMG  = SPTEMP (NEMSG, S(NT,J), TEMISG, EMISG)
               RADJ =  EMG * SIGMA * (S(NT,J)**4 - TDISK**4)
     1               + EMG * SIGMA * (S(NT,J)**4 - TINFTY**4)
            ENDIF
C
            F(NT,J) = (COND(JM1)*(S(NT,J) - S(NT,JM1))/DXM -
     $                 COND(J)*(S(NT,JP1) - S(NT,J)) / DXP   )
     $                          * DXINV
     $             + RHOP*S(NU,J)*DXINV*CPB*(S(NT,JP1) - S(NT,J))
     $             - SRCJ + RADJ
C
         ELSE
C
            F(NT,J) = RHOP*VELOC*DXINV*CPB*(S(NT,J) -
     1                SPTEMP (NTEMP, X(J), XGIVEN, TGIVEN))
C
         ENDIF
C
C        zero-derivative condition for lambda
         F(NL,J) = VELOC * DXINV * (S(NL,J) - S(NL,JP1))
C
1000  CONTINUE
C
C-----BOUNDARY CONDITIONS FAR FROM THE DISK
C
      RHOJ   = RHOJP1
      RHOM   = RHOP
      RHOP   = RHOJ
C
      DXM    = DXP
      DXINV  = 2.0 / (X(JJ) - X(JJ-1))
C
C     calculate a characteristic velocity for satisfaction of
C     Dirichelet conditions
C      IF (S(NU,JJ) .LT. ZERO) THEN
C        VELOC = - S(NU,JJ)
C      ELSEIF (AOMEGA .GT. ZERO) THEN
      IF (AOMEGA .GT. ZERO) THEN
        VELOC = ABS((X(JJ)-X(1)) * AOMEGA)
      ELSEIF (UINF .GT. ZERO) THEN
        VELOC = ABS(UINF)
      ELSE
        VELOC = ONE
      ENDIF
C
C
C     MIXTURE CONTINUITY EQUATION
C     - stays intact up to last node
C       U velocity at top is set by lambda equation
C      (changed to conservation form)
C
      F(NU,JJ) = (RHOP*S(NU,JJ) - RHOM*S(NU,JJ-1)) * DXINV
     $           + 2.0 * RHOJ * S(NV,JJ)
C
C     radial velocity is zero at infinity
      F(NV,JJ) = RHOP * VELOC * DXINV * (S(NV,JJ) - ASPREAD)
C
C     circumferential velocity = swirl velocity at infinity
      F(NW,JJ) = RHOP * VELOC * DXINV * (S(NW,JJ) - ASWIRL)
C
C     energy equation at the top of the domain
      IF (LENRGY) THEN
C
C        determine the specific heat of the mixture at the top
         CALL CKCPBS (S(NT,JJ), EPS, ICKWRK, RCKWRK, CPB)
C
         IF (LHEATFLX) THEN
C           calculate energy source term
C
            IF (QDOT .GT. ZERO) THEN
               XCUT = 2.*WSRC
               IF (ABS(X(JJ) - XSRC) .LE. XCUT) THEN
                  Q0 = 1.283153 * QDOT * SQRT (3 / PI) / WSRC
                  SRCJ = Q0 * EXP(-3.*(X(JJ)-XSRC)**2./WSRC**2.)
               ELSE
                  SRCJ = ZERO
               ENDIF
            ELSE
               SRCJ = ZERO
            ENDIF
C
C           calculate ratiation loss term
            RADJ = ZERO
            IF (NEMSG .GT. 0) THEN
               EMG = SPTEMP (NEMSG, S(NT,JJ), TEMISG, EMISG)
               RADJ = EMG * SIGMA * (S(NT,JJ)**4 - TDISK**4)
     $              + EMG * SIGMA * (S(NT,JJ)**4 - TINFTY**4)
            ENDIF
C
             F(NT,JJ) = (COND(JJ-1)*(S(NT,JJ)-S(NT,JJ-1))/DXM)
     $                        * DXINV
     $           + RHOP*S(NU,JJ)*DXINV*CPB*(TINFTY - S(NT,JJ))
     $           - SRCJ + RADJ
C
         ELSE
           F(NT,JJ) = RHOP * CPB * VELOC * DXINV
     $                     * (S(NT,JJ) - TINFTY)
         ENDIF
C
      ELSE
         F(NT,JJ) = RHOP * CPB * VELOC * DXINV * (S(NT,JJ) -
     1              SPTEMP(NTEMP, X(JJ), XGIVEN, TGIVEN))
      ENDIF
C
      IF (LUINF) THEN
C        the lambda equation determines the specified input velocity
         F(NL,JJ) = - RHOP * VELOC * DXINV * (S(NU,JJ) + UINF)
      ELSE
C        input velocity was not specified, so lambda equation is
C        trivially zero (idel-disk axial velocity)
         F(NL,JJ) = VELOC * DXINV * S(NL,JJ)
      ENDIF
C
C-----ADD THE TIME STEP, IF NEEDED
C
      IF (.NOT. FTIME) RETURN
C
      DO 2500 J = 1, JJ
C
C        add transient term to the species continuity residual;
C        first get mixture density
         CALL CKRHOY (P, S(NT,J), EPS, ICKWRK, RCKWRK, RHOJ)
C
C        obtain the mixture density at the old time step
         CALL CKRHOY (P, SN(NT,J), EPS, ICKWRK, RCKWRK, RHOJN)
C
C        add transient term to the mixture continuity residual
         F(NU,J) = F(NU,J) + (RHOJ - RHOJN) / DT
C
C        add transient term to the radial momentum residual
         F(NV,J) = F(NV,J) + RHOJ * (S(NV,J) - SN(NV,J)) / DT
C
C        add transient term to the circumferential momentum residual
         F(NW,J) = F(NW,J) + RHOJ * (S(NW,J) - SN(NW,J)) / DT
C
C        add transient term to the energy equation residual;
C        determine the specific heat of the misture
         IF (LENRGY) CALL CKCPBS (S(NT,J), EPS, ICKWRK, RCKWRK, CPB)
         F(NT,J) = F(NT,J) + RHOJ * CPB * (S(NT,J) - SN(NT,J)) / DT
C
C        add a transient term to the lambda equation
         IF (.NOT. LUINF) F(NL,J) = F(NL,J) + (S(NL,J) - SN(NL,J)) / DT
C
2500  CONTINUE
C
C     end of SUBROUTINE SPNRFN
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE SPNRPR (LOUT, JJ, NATJ, X, S, ICKWRK, RCKWRK)
C
C  START PROLOGUE
C
C INPUT-
C   LOUT   - UNIT NUMBER FOR OUTPUT.
C   JJ     - NUMBER OF MESH POINTS.
C   X      - THE ARRAY OF MESH POINT LOCATIONS.
C              DIMENSION X(*) AT LEAST JJ
C              CGS UNITS - CM
C   S      - DEPENDENT VARIABLE MATRIX.  THE TEMPERATURES ARE STORED IN
C            T(J)=S(NT,J)*SF(NT), THE MASS FLOW RATES ARE STORED IN
C            FLRT(J)=S(NU,J), AND THE MASS FRACTIONS  ARE IN
C            Y(K,J)=S(NYS+K,J)*SF(NYS+K).
C              DIMENSION S(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION,
C              AND AT LEAST JJ FOR THE SECOND.
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   RCKWRK - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
C
      DIMENSION X(JJ), S(NATJ,JJ), ICKWRK(*), RCKWRK(*)
C
C     print the first line
      WRITE (LOUT, 9000)
      DO 100 J = JJ, 1, -1
         WRITE (LOUT, 9010) J, X(J), S(NT,J), S(NU,J), S(NV,J), S(NW,J)
  100 CONTINUE
C
      RETURN
C
 9000 FORMAT (/8X, 'X(CM)' , 4X, 'T(K)', 7X, 'U(C/S)', 4X, 'V/R(1/S)',
     1         3X, 'W/R(1/S)')
 9010 FORMAT (I4, 2X, F7.4, 10(1PE11.3))
C
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE SPPNT1 (LINKCK, LINKMC, LINKSK, NMAX, LOUT,
     1                   LENIWK, LENRWK, LENCWK, LTOT, ITOT, NTOT,
     2                   JTOT, MAXNEQ, I, R, C, KERR)
C
C  START PROLOGUE
C          SET UP THE POINTERS THAT ARE USED ONLY IN SPIN
C          (AND NOT IN ALE)
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION I(*), R(*)
      CHARACTER*(*) C(*)
      LOGICAL KERR
C
      COMMON /SPSPSP/ NCKW, NMCW, NSKW, NEPS, NWT,  NRE,  NX,   NVIS,
     1                NCON, NTGV, NXGV, ND,   NDTC, NYV,  NS,   NACT,
     2                NFF,  NFN,  NYY,  NXAV, NDKJ, NYAV, NHHH, NCPP,
     3                NXMF, NXMP, NDDN, NSDE, NROP, NTEM, NEMG, NBLI,
     4                NSD0, NSDI, NFLX, NIFX, NRFX, NDRA, NDRC, NCDT,
     5                NDDT, NCMD, NAWT, NKA6, NWDI, NINJ, NINH, NRKF,
     6                NRKR, NRSA, NUSA, NCDF, NCNT, NPRE, NTDS, NTIN,
     7                NOME, NUIN, NRHO, NGFA, NSFA, NSWL, NQDT, NXSR,
     8                NWSR, NTWL, NEMS, NPOW, NTRD, NERD, NVFC, NRU,
     9                NINX, NINW, NITJ, NINM, NSP,  NSN,  NDP,  NDPP,
     +                NFDP, NGDT, IPST, IPND, IPKK, INCF, ICON, ICKW,
     1                IMCW, ISKW, IKR,  IETC, IPMX, IJDX, JCCH, JWCH,
     2                JKCH, JPCH, JECH
C
      COMMON /SPONLY/ NABV, NBLW, NBUF, NTWP, NA,   NPRD, NINT, NWDT,
     1                NSDT, NBUL, NFB,  NFS,  NSDN, NZ,   NZN,  NFZ,
     2                NTFL, ITWP, IKI,  IKP,  IIP,  LAC,  LMK,  JTWP
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      COMMON /SPCON / MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
C
      KERR = .FALSE.
      CALL CKLEN (LINKCK, LOUT, LENICK, LENRCK, LENCCK, IFLAG1)
      CALL MCLEN (LINKMC, LOUT, LENIMC, LENRMC, IFLAG2)
      CALL SKLEN (LINKSK, LOUT, LENISK, LENRSK, LENCSK, IFLAG3)
      IF (IFLAG1.GT.0 .OR. IFLAG2.GT.0 .OR. IFLAG3.GT.0) THEN
         KERR = .TRUE.
         RETURN
      ENDIF
C
C! real chemkin work space
         NCKW = 1
C! real transport work space
         NMCW = NCKW + LENRCK
C! real surface work space
         NSKW = NMCW + LENRMC
         NTOT = NSKW + LENRSK
C
C! integer chemkin space
         ICKW = 1
C! integer transport space
         IMCW = ICKW + LENICK
C! integer surface space
         ISKW = IMCW + LENIMC
         ITOT = ISKW + LENISK
C
C! chem. char. work space
         JCCH  = 1
C! surf. char. work space
         JWCH  = JCCH + LENCCK
         JTOT  = JWCH + LENCSK
C
      IF (ITOT.LT.LENIWK .AND. NTOT.LT.LENRWK .AND. JTOT.LT.LENCWK)
     1    THEN
C
         CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT, I, R, C,
     1                IFLAG)
         KERR = KERR.OR. (IFLAG .GT. 0)
         CALL MCINIT (LINKMC, LOUT, LENIMC, LENRMC, I(IMCW), R(NMCW),
     1                IFLAG)
         KERR = KERR.OR. (IFLAG .GT. 0)
         CALL SKINIT (LENISK, LENRSK, LENCSK, LINKSK, LOUT,
     1                I(ISKW), R(NSKW), C(JWCH), IFLAG)
         KERR = KERR.OR. (IFLAG .GT. 0)
         IF (KERR) RETURN
C
         CALL CKINDX (I, R, MM, KKGAS, II, NFIT)
         CALL SKINDX (I(ISKW), MM, KKGAS, KKSURF, KKBULK, KKTOT,
     1                NPHASE, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2                NLBULK, IISUR)
         CALL SKMXTP (I(ISKW), MAXTP)
         NTR = MAXTP - 1
C
      ENDIF
C
C     set the pointers into the solution vector
      NU   = 1
      NV   = 2
      NW   = 3
      NT   = 4
      NL   = 5
      NYS  = 5
      NY1  = NYS + 1
      NUL = NY1 + KKGAS
      NATJF = NUL
C
C     Jacobian information
      NBEQ = KKSURF + 2*KKBULK + NNSURF
      MAXNEQ = NMAX*NATJF + NBEQ
      IASIZE = (3 * (NATJF + MAX (NATJF, NBEQ) - 1) + 2)
     1                * MAXNEQ
C
C     TWOPNT working space
      LENRTW = 3*NMAX + 9*(NBEQ + NATJF*NMAX)
      LENITW = 3*NMAX
C
C     APPORTION THE FLOATING POINT SPACE
C
C! above bounds for TWOPNT
      NABV = NTOT
C! below bounds for TWOPNT
      NBLW = NABV + NATJF + NBEQ
C! TWOPNT buffer space
      NBUF = NBLW + NATJF + NBEQ
C! TWOPNT working space
      NTWP = NBUF + MAXNEQ
C! surf. spec. solution vector
      NZ   = NTWP + LENRTW
C! surf. spec. at last time
      NZN  = NZ   + KKSURF
C! surf. spec. residual vector
      NFZ  = NZN  + KKSURF
C! Jacobian
      NA   = NFZ  + KKSURF
C! product species mole fractions
      NPRD = NA + IASIZE
C! intermediate species mole fra.
      NINT = NPRD + KKGAS
C! sitdot
      NSDT = NINT + KKGAS
C! total flux over the elements
      NTFL = NSDT + NPHASE
C! production rates
      NWDT = NTFL + MM
C! bulk species mole fractions
      NBUL = NWDT + KKTOT
C! residual for bulk species equations
      NFB  = NBUL + KKBULK
C! "old values" of surface site densities
      NSDN = NFB  + KKBULK
C! residuals for surface site density equations
      NFS  = NSDN + NPHASE
      NTOT = NFS + NPHASE
C
C     APPORTION THE INTEGER SPACE
C
C! intermediate species indicies
      IKI  = ITOT
C! product species indicies
      IKP  = IKI  + KKGAS
C! linpack pivots
      IIP  = IKP  + KKGAS
C  twopnt working space
      ITWP = IIP + MAXNEQ
C
      ITOT = ITWP + LENITW
C
C     APPORTION THE CHARACTER*16 SPACE
C
C! twopnt working space
      JTWP  = JTOT
C!
      JTOT  = JTWP + NATJF + NBEQ
C
C     APPORTION THE LOGICAL SPACE
C
C! adaptive comp. flags
      LAC  = 1
C! new mesh point markers
      LMK  = LAC + NATJF
C
      LTOT = LMK + NMAX - 1
C
C     end of SUBROUTINE SPPNT1
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE SPPNT ( NTOT, ITOT, JTOT, NMAX, MAXNEQ)
C
C  START PROLOGUE
C         SET UP THE ARRAY POINTERS THAT ARE COMMON TO SPIN AND ALE.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SPSPSP/ NCKW, NMCW, NSKW, NEPS, NWT,  NRE,  NX,   NVIS,
     1                NCON, NTGV, NXGV, ND,   NDTC, NYV,  NS,   NACT,
     2                NFF,  NFN,  NYY,  NXAV, NDKJ, NYAV, NHHH, NCPP,
     3                NXMF, NXMP, NDDN, NSDE, NROP, NTEM, NEMG, NBLI,
     4                NSD0, NSDI, NFLX, NIFX, NRFX, NDRA, NDRC, NCDT,
     5                NDDT, NCMD, NAWT, NKA6, NWDI, NINJ, NINH, NRKF,
     6                NRKR, NRSA, NUSA, NCDF, NCNT, NPRE, NTDS, NTIN,
     7                NOME, NUIN, NRHO, NGFA, NSFA, NSWL, NQDT, NXSR,
     8                NWSR, NTWL, NEMS, NPOW, NTRD, NERD, NVFC, NRU,
     9                NINX, NINW, NITJ, NINM, NSP,  NSN,  NDP,  NDPP,
     +                NFDP, NGDT, IPST, IPND, IPKK, INCF, ICON, ICKW,
     1                IMCW, ISKW, IKR,  IETC, IPMX, IJDX, JCCH, JWCH,
     2                JKCH, JPCH, JECH
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      COMMON /SPCON / MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
C
C
C        APPORTION THE BALANCE OF THE FLOATING POINT SPACE
C
C! EPS
      NEPS = NTOT
C! species molecular weights WT
      NWT  = NEPS + KKGAS
C! reactant species mole fract.
      NRE  = NWT  + KKTOT
C! mesh points, X(j)
      NX   = NRE  + KKGAS
C! viscosities at X(j)
      NVIS = NX   + NMAX
C! conductivities at X(j)
      NCON = NVIS + NMAX
C! specified temperatures
      NTGV = NCON + NMAX
C! x locations for spec. temp.
      NXGV = NTGV + NMAX
C! temperatures for specified gas emissivities
      NTEM = NXGV + NMAX
C! emissivities for specified temperatures
      NEMG = NTEM + NMAX
C! diffusion coefficients
      ND   = NEMG + NMAX
C! thermal diffusion coef.
      NDTC = ND   + KKGAS*NMAX
C! YV diffusion velocities
      NYV  = NDTC + KKGAS*NMAX
C! solution vector
      NS   = NYV  + KKGAS*NMAX
C! solution time derivative
      NSP  = NS   + MAXNEQ
C! solution vector at last time
      NSN  = NSP  + MAXNEQ
C! species activities (mole fractions, site fractions, bulk activities)
      NACT = NSN  + MAXNEQ
C! residual vector
      NFF  = NACT + KKTOT
C! res. vector for Jac. perturb.
      NFN  = NFF  + MAXNEQ
C! mass fractions
      NYY  = NFN  + MAXNEQ
C! average mole fractions.
      NXAV = NYY  + KKGAS
C! multicomp. diff. coef.
      NDKJ = NXAV + KKGAS
C! average mass fractions
      NYAV = NDKJ + KKGAS*KKGAS*NMAX
C! species enthalpies
      NHHH = NYAV + KKGAS
C! species specific heats
      NCPP = NHHH + KKTOT
C! species mole fractions
      NXMF = NCPP + KKGAS
C! species mole fractions
      NXMP = NXMF + KKGAS
C! species densities
      NDDN = NXMP + KKGAS
C! site densities
      NSDE = NDDN + KKTOT
C! rate of progress for sur. rea.
      NROP = NSDE + NPHASE
C! sdot for reactions
      NWDI = NROP + IISUR
C! sitdot for reactions
      NSDI = NWDI + KKTOT
C! flux of gas-phase species
      NFLX = NSDI + NPHASE
C! inlet flux over the reagents
      NIFX = NFLX + KKGAS
C! radial flux over the reagents
      NRFX = NIFX + KKGAS
C! creation rates for gas-phase species
      NCDT = NRFX + KKGAS
C! destruction rates for gas-phase species
      NDDT = NCDT + KKGAS
C! creation - destruction (wdot)
      NCMD = NDDT + KKGAS
C! atomic weights
      NAWT = NCMD + KKGAS
C! surface destruction rates for gas-phase species (gdot)
      NGDT = NAWT + MM
C! sens w.r.t. pre-exp constant
      NDRA = NGDT + KKGAS
C! sens w.r.t. concentration
      NDRC = NDRA + KKTOT
C! thermodynamic coefficients a6
      NKA6 = NDRC + KKTOT
C! injected species
      NINJ = NKA6 + NTR
C! H of injected species
      NINH = NINJ + KKGAS
C! initial values of bulk species mole fractions
      NBLI = NINH + KKGAS
C! standard-state site densities for surface phase
      NSD0 = NBLI + KKBULK
C! saved forward rates
      NRKF = NSD0 + NPHASE
C! saved reverse rates
      NRKR = NRKF + II * NMAX
C! saved solution in FUN
      NUSA = NRKR + II * NMAX
C! saved production rates
      NRSA = NUSA + NATJF * NMAX
C! total amount of bulk species grown
      NDP  = NRSA + KKTOT * NMAX
C! time derivative of total amount of bulk species grown
      NDPP = NDP  + KKBULK
C! residual of total amount of bulk species grown
      NFDP = NDPP + KKBULK
C! pressure
      NPRE = NFDP + KKBULK
C! disk temperature
      NTDS = NPRE + 1
C! temperature at infinity
      NTIN = NTDS + 1
C! spin rate
      NOME = NTIN + 1
C! velocity at infinity
      NUIN = NOME + 1
C! density at infinity
      NRHO = NUIN + 1
C! multiplicative constant for gas reaction rates
      NGFA = NRHO + 1
C! multiplicative constant for surface reaction rates
      NSFA = NGFA + 1
C! inlet swirl rate
      NSWL = NSFA + 1
C! energy source term
      NQDT = NSWL + 1
C! x coordinate for center of energy source
      NXSR = NQDT + 1
C! base width of gaussian source term
      NWSR = NXSR + 1
C! ambient temperature for radiation balance
      NTWL = NWSR + 1
C! emissivity of the substrate
      NEMS = NTWL + 1
C! power heating the disk
      NPOW = NEMS + 1
C! temperature of upper radiating disk
      NTRD = NPOW + 1
C! emissivity of a radiating disk above the substrate
      NERD = NTRD + 1
C! vfac term used in radiation balance
      NVFC = NERD + 1
C! cndfac term used in substrate energy balance
      NCDF = NVFC + 1
C! temperature on far side of substrate
      NCNT = NCDF + 1
C! gas constant
      NRU  = NCNT + 1
C! flow rate of injection
      NINM = NRU + 1
C! position of injector
      NINX = NINM + 1
C! width of injector
      NINW = NINX + 1
C! temperature of injector
      NITJ = NINW + 1
      NTOT = NITJ + 1 - 1
C
C     APPORTION THE BALANCE OF THE INTEGER SPACE
C
C! starting species of phases
      IPST = ITOT
C! ending species of phases
      IPND = IPST + NPHASE
C! total species in a phase
      IPKK = IPND + NPHASE
C! reactant species indicies
      IKR  = IPKK + NPHASE
C! etching flag
      IETC = IKR  + KKGAS
C! elemental composition of species
      INCF = IETC + NPHASE
C! conservation of sites array
      ICON = INCF + MM*KKTOT
C! species in a phase with maximum concentration
      IPMX = ICON + NPHASE
C! jindex array over JJ
      IJDX = IPMX + NPHASE
      ITOT = IJDX + NMAX -1
C
C     APPORTION THE BALANCE OF THE CHARACTER*16 SPACE
C
C! species names
      JKCH  = JTOT
C! phase names
      JPCH  = JKCH + KKTOT
C! element names
      JECH  = JPCH + NPHASE
      JTOT  = JECH + MM - 1
C
C     end of SUBROUTINE SPPNT
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE SPPRNT (LOUT, JJ, KFIRST, KLAST, P, X, Z, S, BUL,
     1   SITSOL, ACT, Y, XMF, FLUX, ROP, SDOT, SITDOT, SDOTI, SITDTI,
     2   KNAM, PNAM, SDEN, DEN, WT, H, ICKWRK, RCKWRK, ISKWRK, RSKWRK,
     3   CSKWRK, RMCWRK, LPRNT, ENAM, KNCF, AWT, CDOT, DDOT, WDOT, YV,
     4   REAFLX, RADFLX, LCOMP, LPUTIL, LVCOR, EPS, LREORD, LRADB,
     5   COND, TINFTY, VFAC, ERAD, TRAD, TWAL, EMIS, POWR, LSDEN,
     6   SFAC, LCNDCT, CNDFAC, CNDT, LGDOT, GDOT)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (PI=3.141592654D0, SIGMA = 5.67D-5, ZERO=0.0D0,
     $           ONE=1.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (PI=3.141592654E0, SIGMA = 5.67E-5, ZERO=0.0E0,
C     $           ONE=1.0E0)
C*****END precision > single
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      COMMON /SPCON / MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
C     Integer arrays
      DIMENSION ICKWRK(*), ISKWRK(*), KNCF(MM,KKTOT), KFIRST(NPHASE),
     1          KLAST(NPHASE)
C     Real arrays
      DIMENSION X(JJ), S(NATJ,JJ), Z(*), Y(*), XMF(*), RCKWRK(*),
     1          RSKWRK(*), FLUX(*), SDOTI(*), SITDTI(*), SDOT(*),
     3          SITDOT(*), SDEN(*), DEN(*), ROP(*), WT(*), ACT(*),
     4          H(*), RMCWRK(*), AWT(MM), CDOT(*), DDOT(*), WDOT(*),
     5          YV(KKGAS,*), REAFLX(*), RADFLX(*), EPS(KKGAS),
     6          COND(JJ), BUL(KKBULK), SITSOL(NPHASE), GDOT(*)
C
      DIMENSION YAV(200)
C
      CHARACTER*(*) CSKWRK(*), KNAM(KKTOT), PNAM(NPHASE), ENAM(MM)
      CHARACTER*48 RNAM
C
      LOGICAL LPRNT, IERR, LCOMP, LPUTIL, LVCOR, LREORD, LRADB,
     $        LHEATFLX, LSDEN, LCNDCT, LYCHECK, LGDOT
      SAVE    LHEATFLX, LYCHECK
C
      DATA KPERLN /6/, KPRLN/5/, LT/48/
C
C*****OPTIONS: Use Heat Flux BC at top of domain
C      DATA    LHEATFLX /.TRUE./
C*****END OPTIONS: Use Heat Flux BC at top of domain
C*****OPTIONS: Don't Use Heat Flux BC at top of domain
      DATA    LHEATFLX /.FALSE./
C*****END OPTIONS: Don't Use Heat Flux BC at top of domain
      DATA     LYCHECK /.FALSE./
C
C*****OPTIONS: Check Sum_Yk condition
      LYCHECK = .TRUE.
C*****END OPTIONS: Check Sum_Yk condition
C
C     Do some initial Calculations, using the solution stored
C     in S, Z, and BUL.
C        - Calculate and fill up ACT(*)
C        - Calculate SDOT and SITDOT
C
C     calculate and store surface flux of gas-phase species
      IF (LGDOT) THEN
         DO 10 K = 1, KKGAS
            SDOT(K) = GDOT(K)
            FLUX(K) = SDOT(K)
   10    CONTINUE
      ELSE
         CALL SPSCHM (ICKWRK, RCKWRK, ISKWRK, RSKWRK, S, ACT, Z,
     $             KFIRST, KLAST, P, BUL, LSDEN, SDEN, SITSOL,
     $             SDOT, SITDOT, SFAC)
         DO 20 K = 1, KKGAS
            FLUX(K) = SDOT(K) * WT(K)
   20    CONTINUE
      ENDIF
C
C     calculate rate of progress variables
      IF (IISUR .GT. 0)
     $   CALL SKROP (P, S(NT,1), ACT, SDEN, ISKWRK, RSKWRK, ROP)
C
C     correct ROP variable for an overall reduction in surface
C     chemistry rate
      DO 30 I = 1, IISUR
         ROP(I) = ROP(I)*SFAC
   30 CONTINUE
C
C     calculate surface velocity
      VSUR = ZERO
      DO 50 K = 1, KKGAS
        VSUR = VSUR + SDOT(K) * WT(K)
  50  CONTINUE
      CALL CKRHOY (P, S(NT,1), S(NY1,1), ICKWRK, RCKWRK, RHO)
      VSUR = VSUR / RHO
C
C     print the first line
      CALL CKRP (ICKWRK, RCKWRK, RU, RUC, PA)
      IF (LYCHECK) THEN
        WRITE (LOUT, 7040)
      ELSE
        WRITE (LOUT, 7041)
      ENDIF
C
      DO 100 J = JJ, 1, -1
         CALL CKRHOY (P, S(NT,J), S(NY1,J), ICKWRK, RCKWRK, RHO)
         IF (LYCHECK) THEN
           SOUND = 0.0
           DO 90 K = 1, KKGAS
             SOUND = SOUND + S(NYS+K,J)
  90       CONTINUE
         ELSE
            CALL CKMMWY (S(NY1,J), ICKWRK, RCKWRK, WTM)
            CALL CKCVBS (S(NT,J), S(NY1,J), ICKWRK, RCKWRK, CVBMS)
            CALL CKCPBS (S(NT,J), S(NY1,J), ICKWRK, RCKWRK, CPBMS)
            GAMM = CPBMS/CVBMS
            SOUND = SQRT(GAMM*(RU/WTM)*S(NT,J))
         ENDIF
         IF (J .EQ. JJ) THEN
           WRITE (LOUT, 7010) J, X(J), S(NT,J), S(NU,J), S(NV,J),
     $                        S(NW,J), RHO, S(NL,J), SOUND
         ELSEIF (J .EQ. 1) THEN
           WRITE (LOUT, 7010) J, X(J), S(NT,J), VSUR, S(NV,J),
     $                        S(NW,J), RHO, S(NL,J), SOUND
         ELSE
           WRITE (LOUT, 7010) J, X(J), S(NT,J),
     $                        (S(NU,J-1)*(X(J+1)- X(J)  )+
     $                         S(NU,J)  *(X(J)  - X(J-1)))/
     $                        (X(J+1) - X(J-1)),
     $                        S(NV,J), S(NW,J), RHO, S(NL,J), SOUND
         ENDIF
  100 CONTINUE
C
C     write out gas-phase mass fractions and their fluxes to the
C     surface
      DO 200 L = 1, KKGAS, KPERLN
         K2   = MIN (L+KPERLN-1,KKGAS)
         WRITE (LOUT, 7030) (KNAM(K), K = L, K2)
C
         DO 160 J = JJ, 1, -1
            DO 120 K = 1, KKGAS
               Y(K) = S(NYS+K, J)
  120       CONTINUE
            CALL CKYTX (Y, ICKWRK, RCKWRK, XMF)
            DO 140 K = 1, KKGAS
               Y(K) = XMF(K)
  140       CONTINUE
            WRITE (LOUT, 7010) J, X(J), (Y(K), K = L, K2)
  160     CONTINUE
C
        WRITE (LOUT, 7020) (FLUX(K), K = L, K2)
  200 CONTINUE
C
      IF (KKSURF .GT. 0) THEN
C        write out surface species site fractions
         WRITE (LOUT, '(/A/)') ' SURFACE SPECIES SITE FRACTIONS'
C
C        surface species
         DO 220 N = NFSURF, NLSURF
            WRITE (LOUT, 7050) (N-1), PNAM(N), SDEN(N)
            WRITE (LOUT, 7060)
     1       (KNAM(K), Z(K-KKGAS), K=KFIRST(N),KLAST(N))
  220    CONTINUE
      ENDIF
C
      IF (IISUR.GT.0 .AND. (LPRNT)) THEN
C
            WRITE (LOUT, '(/A)')
     1' CONTRIBUTION OF SURFACE REACTIONS TO SPECIES DESTRUCTION RATES'
C
C        surface reaction rates
         DO 260 L = 1, KKGAS, KPRLN
            K2 = MIN (L+KPRLN-1,KKGAS)
            WRITE (LOUT, 7070) (KNAM(K), K = L, K2)
C
            DO 240 I = 1, IISUR
               RNAM = '  '
               CALL SKSYMR (I, LOUT, ISKWRK, RSKWRK, CSKWRK, LT,
     1                      RNAM,IERR)
               CALL SKRATI (I, ROP, ISKWRK, RSKWRK, SDOTI, SITDTI)
               WRITE (LOUT, 7080) I, RNAM, (SDOTI(K), K = L, K2)
  240       CONTINUE
            WRITE (LOUT, 7090) (SDOT(K), K = L, K2)
  260    CONTINUE
C
         IF (KKSURF .GT. 0) THEN
            DO 300 L = KFIRST(NFSURF), KLAST(NLSURF), KPRLN
               K2 = MIN (L+KPRLN-1,KLAST(NLSURF))
               WRITE (LOUT, 7070) (KNAM(K), K = L, K2)
C
               DO 280 I = 1, IISUR
                  RNAM = '  '
                  CALL SKSYMR (I, LOUT, ISKWRK, RSKWRK, CSKWRK, LT,
     1                      RNAM,IERR)
                  CALL SKRATI (I, ROP, ISKWRK, RSKWRK, SDOTI, SITDTI)
                  WRITE (LOUT, 7080) I, RNAM, (SDOTI(K), K = L, K2)
  280          CONTINUE
               WRITE (LOUT, 7090) (SDOT(K), K = L, K2)
  300       CONTINUE
         ENDIF
      ENDIF
C
      IF (KKBULK .GT. 0) THEN
C       write out the rates of production of bulk species from each
C       gas-phase reaction
C
        IF (LPRNT) THEN
           DO 340 L = KFIRST(NFBULK), KLAST(NLBULK), KPRLN
              K2 = MIN (L+KPRLN-1,KKTOT)
              WRITE (LOUT, 7070) (KNAM(K), K = L, K2)
C
              DO 320 I = 1, IISUR
                 RNAM = '  '
                 CALL SKSYMR (I, LOUT, ISKWRK, RSKWRK, CSKWRK, LT,
     $                        RNAM, IERR)
                 CALL SKRATI (I, ROP, ISKWRK, RSKWRK, SDOTI, SITDTI)
                 WRITE (LOUT, 7080) I, RNAM, (SDOTI(K), K = L, K2)
  320         CONTINUE
              WRITE (LOUT, 7090) (SDOT(K), K = L, K2)
  340      CONTINUE
        ENDIF
C
C       print total deposition rates for all solid species
        WRITE (LOUT, '(/A)')
     $' DEPOSITION RATE       (CM/SEC) (MICRON/MIN)      MOLE FRACTION'
C
        CALL SKDEN  (P, S(NT,1), ACT, SDEN, ISKWRK, RSKWRK, DEN)
C
        DO 380 N = NFBULK, NLBULK
           WRITE (LOUT, 7055) (N-MAX(1,NLSURF)),PNAM(N)
           DO 380 K = KFIRST(N), KLAST(N)
              DRAT = SDOT(K)*WT(K)/DEN(K)
              WRITE (LOUT, 7100) KNAM(K), DRAT, DRAT*60.0E4, ACT(K)
  380   CONTINUE
      ENDIF
C
C     print heat release due to surface reaction
      CALL SKHML (S(NT,1), ISKWRK, RSKWRK, H)
      QCHEM = ZERO
      DO 390 K = 1, KKTOT
         QCHEM = QCHEM + SDOT(K)*H(K)
  390 CONTINUE
      QCHEMW = QCHEM * 1.0E-7
      WRITE (LOUT, 7114) QCHEM, QCHEMW
C
C     print heat flux loss to burner and to substrate
      Q1   = COND(1)   *(S(NT,2)  - S(NT,1))/(X(2) - X(1))
      QJJ  = COND(JJ-1)*(S(NT,JJ) - S(NT,JJ-1))/(X(JJ) - X(JJ-1))
      Q1W  = Q1*1.0E-7
      QJJW = QJJ *1.0E-7
      WRITE (LOUT, 7110) QJJ, QJJW
      WRITE (LOUT, 7112) Q1, Q1W
C
C-------------------------------------------------------------------
C  Print out Utilization tables and check global mass balances
C-------------------------------------------------------------------
C
      IF (.NOT. LPUTIL) RETURN
C
C     zero counters
      TOTINL = ZERO
      TOTCRE = ZERO
      TOTDES = ZERO
      TOTRAD = ZERO
      TOTSRF = ZERO
C
C     print headers for utilization tables
C
      WRITE (LOUT, 453) 'UTILIZATION TABLES'
 453  FORMAT(//105('-')/T20,A/105('-')/)
      WRITE (LOUT,454)
 454  FORMAT(/'Species Utilization Table (all units are in gm/cm2Sec)'
     $       /105('-'))
      WRITE (LOUT, 455)
     $     '   Species      Inlet      Creation  Destruction',
     $     'Radial   Surface     Total    Utilizatn   Consumption'
 455  FORMAT(1X, A, 3X, A)
      WRITE (LOUT, 456) 'Flux', 'Rate', 'Rate', 'Flux', 'Flux',
     $                  'Mass_Balnc', 'Fraction', 'Fraction'
 456  FORMAT(T18,A, T29,A, T40,A, T53,A, T62,A,
     $       T71,A, T86, A, T97, A/105('-'))
C
C     find the density at the top node
      CALL CKRHOY (P, S(NT,JJ), S(NY1,JJ), ICKWRK, RCKWRK, RHOJJ)
C
C     find the average density RHOM between JJ and JJ-1
      JJM1 = JJ - 1
      DO 590 K = 1, KKGAS
         N = NYS + K
         YAV(K) = 0.5 * (S(N,JJ) + S(N,JJM1))
590   CONTINUE
      TAV = 0.5 * (S(NT,JJ) + S(NT,JJM1))
      CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOM)
C
C      Find the fix-up diffusion velocity at JJ-1/2
C        - When a special equation is used, where the
C          sum of mass fractions equals one equation is
C          employed, the diffusion velocity for that
C          species must be calculated from the other
C          diffusion velocities.
C
      IF (.NOT. LVCOR) THEN
         YVSUM = ZERO
         DO 554 K = 1, KKGAS
           YVSUM = YVSUM + YV(K,JJM1)
554      CONTINUE
         KINDEX = KKGAS
         IF (LREORD) THEN
           DO 555 K = 1, KKGAS-1
             IF (S(NYS+K,JJM1) .GT. S(NYS+KINDEX,JJM1)) KINDEX = K
  555      CONTINUE
         ENDIF
         YVFIX = YV(KINDEX,JJM1) - YVSUM
      ENDIF
C
C----------------------------------------------------------------------
C           Loop over each species, calculating a mass balance
C           just for that species
C----------------------------------------------------------------------
C
      DO  1000 K = 1, KKGAS
C
C        zero counters
         UTIL    = ZERO
         CONS    = ZERO
C
C        contributions from J=1, the susceptor
         SURFLX = - (SDOT(K) * WT(K))
C
         XDIFF = (X(2) - X(1)) * 0.5
         CALL CKRHOY (P, S(NT,1), S(NY1,1), ICKWRK, RCKWRK, RHOJ)
         CALL CKCDYP (P, S(NT,1), S(NY1,1), ICKWRK, RCKWRK, CDOT,
     1                DDOT)
         CREFLX = CDOT(K) * WT(K) * XDIFF
         DESFLX = DDOT(K) * WT(K) * XDIFF
         RADFLX(K) = 2.0 * RHOJ * S(NYS+K,1) * S(NV,1) * XDIFF
C
C        contributions from the interior nodes
         DO 950 J = 2, JJ-1
C
            CALL CKCDYP (P, S(NT,J), S(NY1,J), ICKWRK, RCKWRK, CDOT,
     1                   DDOT)
            XDIFF = (X(J+1) - X(J-1))*0.5
            CALL CKRHOY (P, S(NT,J), S(NY1,J), ICKWRK, RCKWRK, RHOJ)
C
C           radial flux
            RADFLX(K) = RADFLX(K) + 2.*RHOJ*S(NYS+K,J)*S(NV,J)*XDIFF
C
C           creation flux
            CREFLX = CREFLX + CDOT(K) * WT(K) * XDIFF
C
C           destruction flux
            DESFLX = DESFLX + DDOT(K) * WT(K) * XDIFF
C
  950    CONTINUE
C
C        contributions from J=JJ, the top of the domain
         JJM1 = JJ - 1
         XDIFF = (X(JJ) - X(JJM1)) * 0.5
         CALL CKCDYP (P, S(NT,JJ), S(NY1,JJ), ICKWRK, RCKWRK, CDOT,
     1                DDOT)
C
         IF (LCOMP) THEN
C           base calculation - but watch out for non-correction velocity
C           cases!
C
            IF ((.NOT. LVCOR) .AND. (KINDEX .EQ. K)) THEN
              REAFLX(K) = - RHOM * S(NYS+K,JJ) * S(NU,JJM1)
     $                    - RHOM * YVFIX
            ELSE
              REAFLX(K) = - RHOM * S(NYS+K,JJ) * S(NU,JJM1)
     $                    - RHOM * YV(K,JJM1)
            ENDIF
C
         ELSE
C
C           creation flux
            CREFLX = CREFLX + CDOT(K) * WT(K) * XDIFF
C           destruction flux
            DESFLX = DESFLX + DDOT(K) * WT(K) * XDIFF
            REAFLX(K) = - RHOJJ * EPS(K) * S(NU,JJ)
         ENDIF
C
         RADFLX(K) = RADFLX(K) + 2.*RHOJJ*S(NYS+K,JJ)*S(NV,JJ)
     $                           *XDIFF
C
C        Figure out the fractional amount of reactant that is
C        utilized in surface reactions.
C        If a species is not a reactant, use this to find out
C        the fractional amount of that species created in
C        gas phase reactions that is utilized in surface reactions.
C
         IF (REAFLX(K) .GT. ABS(RHOJJ * 1.0E-10 * S(NU,JJ))) THEN
           UTIL = SURFLX / REAFLX(K)
         ELSEIF (CREFLX .GT. ZERO) THEN
           UTIL = SURFLX / CREFLX
         ELSE
           UTIL = ZERO
         ENDIF
C
C        calculate the fractional amount of a reactant
C        consumed either by gas-phase or surface-phase reactions
         IF (REAFLX(K) .GT. ZERO) THEN
            CONS = (DESFLX - CREFLX + SURFLX) / REAFLX(K)
         ELSE
            CONS = ZERO
         ENDIF
C
C        mass balance
         RIGHT = REAFLX(K) + CREFLX
         RLEFT = DESFLX + RADFLX(K) + SURFLX
         BAL   = (RIGHT - RLEFT)
C
         WRITE (LOUT, '(A16, 8(1PG11.4))')
     $                  KNAM(K), REAFLX(K), CREFLX, DESFLX,
     $                  RADFLX(K), SURFLX, BAL, UTIL, CONS
C
         TOTINL = TOTINL + REAFLX(K)
         TOTCRE = TOTCRE + CREFLX
         TOTDES = TOTDES + DESFLX
         TOTRAD = TOTRAD + RADFLX(K)
         TOTSRF = TOTSRF + SURFLX
C
1000  CONTINUE
C
C     find radial flux of bulk flow as check
      CALL CKRHOY (P, S(NT,1), S(NY1,1), ICKWRK, RCKWRK, RHOJ)
      BLKOUT = RHOJ * S(NV,1) * (X(2) - X(1))
      DO 1010 J = 2, JJ-1
        CALL CKRHOY (P, S(NT,J), S(NY1,J), ICKWRK, RCKWRK, RHOJ)
        BLKOUT = BLKOUT + RHOJ*S(NV,J) * (X(J+1)-X(J-1))
1010  CONTINUE
      CALL CKRHOY (P, S(NT,JJ), S(NY1,JJ), ICKWRK, RCKWRK, RHOJ)
      BLKOUT = BLKOUT + RHOJ*S(NV,JJ) * (X(JJ)-X(JJ-1))
C
C     find the bulk flow at the inlet
      BLKINL = - RHOJJ * S(NU,JJ)
C
C     find the bulk flow to the surface, J=1
      DO 1990 K = 1, KKGAS
         N = NYS + K
         YAV(K) = 0.5 * (S(N,2) + S(N,1))
1990  CONTINUE
      TAV = 0.5 * (S(NT,2) + S(NT,1))
      CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOJ)
      BLKSRF = - RHOJ * S(NU,1)
C
      TOTR = TOTINL + TOTCRE
      TOTL = TOTDES + TOTRAD + TOTSRF
      TOTBAL = (TOTR - TOTL)
C
      WRITE (LOUT, 1354) 'TOTAL',
     $                  TOTINL,TOTCRE,TOTDES,TOTRAD,TOTSRF,TOTBAL
1354  FORMAT(105('-')/
     $          A16,1PG11.4, 1PG11.4, 1PG11.4,
     $                   1PG11.4, 1PG11.4, 1PG11.4)
      WRITE (LOUT, '(A16,1PG11.4, 22X, 1PG11.4, 1PG11.4, 1PG11.4)')
     $                'CHECK', BLKINL, BLKOUT, BLKSRF,
     $                (BLKINL - BLKOUT - BLKSRF)
C
C-----------------------------------------------------------------------
C     Calculate Balances for the Elements - obtained from the
C     balances for the species.
C-----------------------------------------------------------------------
C
      WRITE (LOUT,1454)
1454  FORMAT(/'Element Utilization Table (all units are in gm/cm2Sec)'
     $       /105('-'))
      WRITE (LOUT,'(1X,A,T14,A,T28,A,T43,A,T58,A,T73,A)')
     1      'Atom','Inlet','Radial','Surface','Total',
     2      'Incorporation'
      WRITE (LOUT, 1456) 'Flux', 'Flux', 'Flux', 'Mass_Balnc',
     $      'Fraction'
1456  FORMAT(T14,A, T28,A, T43,A, T56,A, T76, A /105('-'))
C
      DO 1500 M = 1, MM
C
C        element flux at inlet
         AINFLX = ZERO
C        bulk flux at surface
         BLKFLX = ZERO
C        radial flux
         OUTFLO = ZERO
C
         DO 1050 K = 1, KKGAS
            ARAT = KNCF(M,K) * AWT(M) / WT(K)
            AINFLX = AINFLX + REAFLX(K) * ARAT
            OUTFLO = OUTFLO + RADFLX(K) * ARAT
 1050    CONTINUE
C
         DO 1080 K = KKGAS + KKSURF + 1, KKGAS + KKSURF + KKBULK
            BLKFLX = BLKFLX + SDOT(K) * AWT(M) * KNCF(M,K)
 1080    CONTINUE
C
         IF (AINFLX .GT. ZERO) THEN
           EINC = BLKFLX / AINFLX
         ELSE
           EINC = ZERO
         ENDIF
         BAL = (AINFLX - BLKFLX - OUTFLO)
         WRITE (LOUT,
     1'(1X,A9,1PE13.5,T25,1PE13.5,T40,1PE13.5,T55,1PE13.5,T70,1PE13.5)')
     2         ENAM(M), AINFLX, OUTFLO, BLKFLX, BAL, EINC
 1500 CONTINUE
      WRITE (LOUT,1501)
 1501 FORMAT(105('-')/)
C
C -----------------------------------------------------------------
C           Energy Equation Balance Section
C -----------------------------------------------------------------
C
C     Base Normalization for the enthalpy flow terms (ergs/gm)
      IF (LCOMP) THEN
        CALL CKHBMS (TINFTY,       EPS, ICKWRK, RCKWRK, HBASE)
        CALL CKCPBS (TINFTY,       EPS, ICKWRK, RCKWRK, CPBASE)
      ELSE
        CALL CKHBMS (TINFTY, S(NY1,JJ), ICKWRK, RCKWRK, HBASE)
        CALL CKCPBS (TINFTY, S(NY1,JJ), ICKWRK, RCKWRK, CPBASE)
      ENDIF
      HBASE = HBASE - CPBASE*TINFTY
C
C     contributions from J = 1
      XDIFF = 0.5 * (X(2) - X(1))
      CALL CKRHOY (P, S(NT,1), S(NY1,1), ICKWRK, RCKWRK, RHOJ)
C
      IF (LRADB) THEN
         ECONDB   = ZERO
         EDIFFB   = ZERO
         EDIFFB_N = ZERO
         ECONVB   = ZERO
         ECONVB_N = ZERO
C
C        calculate the surface radiation terms
         RADIN = VFAC*SIGMA*TRAD**4
     $           + (ONE - VFAC)*SIGMA*TWAL**4
         RADOUT = EMIS*SIGMA*S(NT, 1)**4 + (ONE - EMIS)*RADIN
         RADNET = RADIN - RADOUT
C
C        calculate the steady state bulk phase enthalpy creation term
         EBULK   = ZERO
         EBULK_N = ZERO
         IF (KKBULK .GT. 0 .AND. NNBULK .GT. 0) THEN
           CALL SKHML (S(NT,1), ISKWRK, RSKWRK, H)
           DO 1840 K = KFIRST(NFBULK) , KLAST(NLBULK)
             EBULK   = EBULK   + SDOT(K) * H(K)
             EBULK_N = EBULK_N + SDOT(K) * (H(K)-HBASE* WT(K))
 1840      CONTINUE
         ENDIF
C
C        calculate the conductive losses to the backside
         IF (LCNDCT) THEN
           EHTRAN = - CNDFAC * (S(NT,1) - CNDT)
         ELSE
           EHTRAN = ZERO
         ENDIF
      ELSE
C
        CALL CKHBMS (S(NT,2), S(NY1,2), ICKWRK, RCKWRK, HBMSJP1)
C
        TAV = 0.5 * (S(NT,2) + S(NT,1))
        DO 1850 K = 1, KKGAS
          N = NYS + K
          YAV(K) = 0.5 * (S(N,2) + S(N,1))
 1850   CONTINUE
        CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOP)
        CALL CKHMS (TAV, ICKWRK, RCKWRK, H)
C
        ECONVB   = + RHOP * S(NU,1) *  HBMSJP1
        ECONVB_N = + RHOP * S(NU,1) * (HBMSJP1 - HBASE)
        ECONDB = - COND(1) * (S(NT,2) - S(NT,1)) / (X(2) - X(1))
        EDIFFB   = ZERO
C
        DO 1900 K = 1, KKGAS
          EDIFFB   = EDIFFB   + RHOP * YV(K,1) *  H(K)
 1900   CONTINUE
      ENDIF
      CALL CKHBMS (S(NT,1), S(NY1,1), ICKWRK, RCKWRK, HBMSJ)
      ECONVR   = 2.*RHOJ*S(NV,1)*XDIFF* HBMSJ
      ECONVR_N = 2.*RHOJ*S(NV,1)*XDIFF*(HBMSJ-HBASE)
C
      DO 2000 J = 2, JJ-1
C        contributions from interior nodes
C
         CALL CKRHOY (P, S(NT,J), S(NY1,J), ICKWRK, RCKWRK, RHOJ)
         RHODIF = RHOJ*S(NV,J)*(X(J+1) - X(J-1))
         CALL CKHBMS (S(NT,J), S(NY1,J), ICKWRK, RCKWRK, HBMSJ)
         ECONVR   = ECONVR   + RHODIF* HBMSJ
         ECONVR_N = ECONVR_N + RHODIF*(HBMSJ-HBASE)
C
 2000 CONTINUE
C
C     contributions from the top of the domain, J=JJ
      JJM1 = J - 1
      CALL CKRHOY (P, S(NT,JJ), S(NY1,JJ), ICKWRK, RCKWRK, RHOJ)
      RHODIF = RHOJ*S(NV,JJ)*(X(JJ) - X(JJM1))
      TAV = 0.5 * (S(NT,JJ) + S(NT,JJM1))
C
      DO 2110 K = 1, KKGAS
        N = NYS + K
        YAV(K) = 0.5 * (S(N,JJ) + S(N,JJM1))
 2110 CONTINUE
C
      CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOM)
      CALL CKHMS (TAV, ICKWRK, RCKWRK, H)
      CALL CKHBMS (S(NT,JJ), S(NY1,JJ), ICKWRK, RCKWRK, HBMSJ)
      ECONVR   = ECONVR   + RHODIF* HBMSJ
      ECONVR_N = ECONVR_N + RHODIF*(HBMSJ-HBASE)
C
C     convective flux of enthalpy into domain from the top
      IF (.NOT. LHEATFLX) THEN
c       convective flux of enthalpy into domain from the top
        ECONVT   = - RHOJ * S(NU,JJ) *  HBMSJ
        ECONVT_N = - RHOJ * S(NU,JJ) * (HBMSJ-HBASE)
C       diffusive flux of enthalpy into domain from the top from
C       thermal conduction
         JM1 = JJ - 1
         ECONDT = + COND(JJM1) * (S(NT,JJ)-S(NT,JJM1)) /
     1              (X(JJ) - X(JJM1))
C       diffusive flux of enthalpy into domain from the top from
C       species diffusion
        EDIFFT = ZERO
        DO 2400 K = 1, KKGAS
           EDIFFT = EDIFFT - RHOM * YV(K,JJM1) * H(K)
 2400   CONTINUE
      ELSE
        CALL CKHBMS (TINFTY, EPS, ICKWRK, RCKWRK, HBMSJ)
        ECONVT   = - RHOJ * S(NU,JJ) *  HBMSJ
        ECONVT_N = - RHOJ * S(NU,JJ) * (HBMSJ-HBASE)
        ECONDT = ZERO
        EDIFFT = ZERO
      ENDIF
C
C     print out the energy balance
      WRITE (LOUT, 8000)
      WRITE (LOUT, 8002) ECONVT, ECONVT_N
      WRITE (LOUT, 8004) ECONDT, ECONDT
      WRITE (LOUT, 8006) EDIFFT, EDIFFT
      TOTALT   = ECONDT + ECONVT   + EDIFFT
      TOTALT_N = ECONDT + ECONVT_N + EDIFFT
      WRITE (LOUT, 8008) TOTALT, TOTALT_N
      WRITE (LOUT, 8010) ECONVR, ECONVR_N
      WRITE (LOUT, 8012) ECONVR, ECONVR_N
      IF (LRADB) THEN
        WRITE (LOUT, 8022) POWR,   POWR
        WRITE (LOUT, 8024) RADNET, RADNET
        WRITE (LOUT, 8025) EHTRAN, EHTRAN
        TOTALB   = POWR + RADNET + EHTRAN
        TOTALB_N = POWR + RADNET + EHTRAN
        WRITE (LOUT, 8020) TOTALB, TOTALB_N
         WRITE (LOUT, 8026) EBULK,  EBULK_N
         WRITE (LOUT, 8028) EBULK,  EBULK_N
      ELSE
        WRITE (LOUT, 8014) ECONVB, ECONVB_N
        WRITE (LOUT, 8016) ECONDB, ECONDB
        WRITE (LOUT, 8018) EDIFFB, EDIFFB
        TOTALB   = ECONDB + ECONVB   + EDIFFB
        TOTALB_N = ECONDB + ECONVB_N + EDIFFB
        WRITE (LOUT, 8020) TOTALB, TOTALB_N
         EBULK = ZERO
         EBULK_N = ZERO
      ENDIF
C
      ETOT   = - TOTALB   + ECONVR   - TOTALT   + EBULK
      ETOT_N = - TOTALB_N + ECONVR_N - TOTALT_N + EBULK_N
C
      WRITE (LOUT, 8030) ETOT, ETOT_N
C
      EMAX = MAX (ABS(ECONVT), ABS(ECONDT), ABS(EDIFFT),
     1            ABS(ECONVR), ABS(ECONVB), ABS(ECONDB),
     2            ABS(EDIFFB))
      EMAX_N = MAX (ABS(ECONVT_N), ABS(ECONDT), ABS(EDIFFT),
     1              ABS(ECONVR_N), ABS(ECONVB_N), ABS(ECONDB),
     2              ABS(EDIFFB))
C
      WRITE (LOUT, 8040) ETOT/EMAX, ETOT_N/EMAX_N
C
7010  FORMAT (I4, 2X, 1PG12.4,   10(1PE11.3))
7020  FORMAT (' FLX(G/CM2 S)', 5X , 10(1PE11.3))
7030  FORMAT (/9X, 'X(CM)' , 7X, 10(A10, 1X))
7040  FORMAT (/9X, 'X(CM)' , 7X, 'T(K)', 7X, 'U(C/S)', 4X, 'V/R(1/S)',
     1         3X, 'W/R(1/S)', 2X, 'RHO(G/CC)',2X, 'LAMBDA',6X,'Sum_Yk')
7041  FORMAT (/9X, 'X(CM)' , 7X, 'T(K)', 7X, 'U(C/S)', 4X, 'V/R(1/S)',
     1         3X, 'W/R(1/S)', 2X, 'RHO(G/CC)',2X, 'LAMBDA',6X,'A(C/S)')
7050  FORMAT (' SITE ',I2,': ', A16, ',  SITE DENSITY ', 1PE15.6)
7055  FORMAT (' BULK ',I2,': ', A16)
7060  FORMAT (4X,  A12, 1PE11.3)
7070  FORMAT (/57X, 5(A10, 2X))
7080  FORMAT (I4, 2X, A48, 5(1PE12.4))
7090  FORMAT (' TOTAL (MOL/CM**2-SEC)', 32X, 5(1PE12.4))
7100  FORMAT (4X, A16, 1PE11.3, 1PE11.3, 7X, 1PE11.3)
7110  FORMAT (' HEAT FLUX TO INLET:     ',5X, 1PE11.3,' ERG/CM2-S',
     1        '  (',1PE11.3,' W/CM2)')
7112  FORMAT (' HEAT FLUX TO SUBSTRATE: ', 5X, 1PE11.3,' ERG/CM2-S',
     1       '  (',1PE11.3,' W/CM2)')
7114  FORMAT (/,' HEAT FROM SURFACE REACTIONS: ', 1PE11.3,' ERG/CM2-S',
     1          '  (',1PE11.3,' W/CM2)')
C
8000  FORMAT (/,' ENTHALPY BALANCES: (all units are erg/cm**2sec)'//
     $        36X,' UNNORMALIZED       NORMALIZED'/
     $        36X,'-------------------------------')
8002  FORMAT (1X,'TOP:  Convective Flux In         = ',
     $                                     1PE12.5,5X,1PE12.5)
8004  FORMAT (1X,'      Conductive Flux In         = ',
     $                                     1PE12.5,5X,1PE12.5)
8006  FORMAT (1X,'      SpecDiffEnrgy Flux In      = ',
     $                                     1PE12.5,5X,1PE12.5)
8008  FORMAT (1X,'TOTAL TOP                        = ',
     $                                     1PE12.5,5X,1PE12.5/)
8010  FORMAT (1X,'SIDE: Radial Convective Flux out = ',
     $                                     1PE12.5,5X,1PE12.5)
8012  FORMAT (1X,'TOTAL SIDE                       = ',
     $                                     1PE12.5,5X,1PE12.5/)
8014  FORMAT (1X,'BOT:  Convective Flux Into Domain= ',
     $                                     1PE12.5,5X,1PE12.5)
8016  FORMAT (1X,'      Conductive Flux Into Domain= ',
     $                                     1PE12.5,5X,1PE12.5)
8018  FORMAT (1X,'      SpecDiffEnrgy Flux Into    = ',
     $                                     1PE12.5,5X,1PE12.5)
*
8022  FORMAT (1X,'      Power Input into disk      = ',
     $                                     1PE12.5,5X,1PE12.5)
8024  FORMAT (1X,'      Radiation Input into disk  = ',
     $                                     1PE12.5,5X,1PE12.5)
8025  FORMAT (1X,'      HeatTran into disk from Bck= ',
     $                                     1PE12.5,5X,1PE12.5)
8020  FORMAT (1X,'TOTAL BOTTOM                     = ',
     $                                     1PE12.5,5X,1PE12.5/)
8026  FORMAT (1X,'RETAIN:Bulk Enthalpy CreationRate= ',
     $                                     1PE12.5,5X,1PE12.5)
8028  FORMAT (1X,'TOTAL RETAIN RATE                = ',
     $                                     1PE12.5,5X,1PE12.5/)
8030  FORMAT (1X,'NET TOTAL                        = ',
     $                                     1PE12.5,5X,1PE12.5)
8040  FORMAT (1X,'NET Fractional Energy Balance    = ',
     $                                     1PE12.5,5X,1PE12.5/)
C
C     end of SUBROUTINE SPPRNT
      RETURN
      END
C--------------------------------------------------------------------
C
      SUBROUTINE SPRDKY (PNAM, KFIRST, KLAST, KKPHAS, NMAX, LIN, LOUT,
     1     PATM, LTIME2, LUSTGV, LENRGY, RSTCNT, LCHEM, LTDIF, LVCOR,
     2     LUMESH, LRSTRT, LCNTUE, LMULTI, LUINF, LASEN,  LPRNT, LNONR,
     3     LRADB, SFLR, P, NPTS, NTOT, NADP, X, NREAC, NINTM, NPROD,
     4     REAC, XINTM, PROD, KR, KI, KP, Z, TINFTY, TDISK, AOMEGA,
     5     UINF, ASWIRL, XSTR, XCEN, XEND, WMIX, NTEMP, XGIVEN, TGIVEN,
     6     N1CALL, MFILE, KNAM, ACT, GFAC, SFAC, LREORD, QDOT, XSRC,
     7     WSRC,SPOS, LREGRD, JJREGD, PCTADP, RATGTC, LHSEN, IVALUE,
     8     RVALUE, LVALUE, FINJ, XINJ, WINJ, SINJ, TINJ, HINJ,
     9     TWAL, EMIS, POWR, TRAD, ERAD, VFAC, LCNDCT, CNDFAC, CNDT,
     +     NEMS, TEMISG, EMISG, IETCH, SDEN, SDEN0, ISKWRK, RSKWRK,
     1     BUL, BULIN, LSDEN, LCOMP, CNTRLS, LTOFF, LPUTIL, FLMX, FLMT,
     2     LFLAME, FLMW, LOLDSP, ASPREAD, GDOT, LGDOT, KERR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INTEGER
     +   QADAPT, QLEVD, QLEVM, QPADD, QSSABS, QSSAGE, QSSREL, QSTDY,
     +   QSTPS0, QSTPS1, QSTPS2, QSTRID, QTDAGE, QTDABS, QTDEC, QTDREL,
     +   QTINC, QTMAX, QTMIN, QTOL0, QTOL1, QTOL2, CNTRLS
      COMMON /QPTR/
     +   QADAPT, QLEVD, QLEVM, QPADD, QSSABS, QSSAGE, QSSREL, QSTDY,
     +   QSTPS0, QSTPS1, QSTPS2, QSTRID, QTDAGE, QTDABS, QTDEC, QTDREL,
     +   QTINC, QTMAX, QTMIN, QTOL0, QTOL1, QTOL2
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      COMMON /SPCON / MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
C     Integer arrays
      INTEGER   ISKWRK(*), KR(*), KI(*), KP(*), KFIRST(NPHASE),
     1          KLAST(NPHASE), KKPHAS(NPHASE), IVALUE(CNTRLS,*),
     2          IETCH(*)
C     Real arrays
      DIMENSION RSKWRK(*), REAC(*), XINTM(*), PROD(*),
     1          XGIVEN(*), TGIVEN(*), X(*),
     2          Z(*), VALUE(5), ACT(*),
     3          RVALUE(CNTRLS,*),
     4          SINJ(*), TEMISG(*), EMISG(*), SDEN(*), SDEN0(*),
     5          BUL(*), BULIN(*),  GDOT(*)
C
      CHARACTER*(*) KNAM(KKTOT), PNAM(NPHASE)
      CHARACTER LINE*80, KEYWRD*4, CKCHUP*4
      EXTERNAL CKCHUP
C
      LOGICAL LTIME2, LENRGY, LUINF, RSTCNT, LTDIF, LVCOR, LUMESH,
     1        LRSTRT, LCNTUE, LUSTGV, LASEN, LMULTI, LSTAG, IERR,
     2        KERR, NEC(10), NOPT(5), CNTNUD, LFIRST, LCHEM, LREORD,
     3        LPRNT, LNONR, LRADB, LREGRD, LHSEN, LVALUE(CNTRLS,*),
     4        LCNDCT, LSDEN, LCOMP, LTOFF, LPUTIL, LFLAME, LOLDSP,
     5        LFRSTI, LFRSTG, LGDOT
C
      SAVE    NEC, NOPT, CNTNUD, LSTAG
      DATA NEC/10*.FALSE./, NOPT/5*.FALSE./, CNTNUD/.FALSE./
C
      KERR = .FALSE.
C
C     INITIALIZE VARIABLES
C
      IF (LCNTUE) THEN
C
         LCNTUE = .FALSE.
         CNTNUD = .TRUE.
         LFIRST = .TRUE.
         LFRSTI = .TRUE.
         LFRSTG = .TRUE.
         NP = 0
         NEMSG = 0
C
C        need only Q pointers for continuation input values
         CALL TWSTRT (-1, LOUT, NMAX, CNTRLS, LVALUE, IVALUE,
     1                RVALUE, KERR)
      ELSE
C
C        need Q pointers and initial default values
         CALL TWSTRT (0, LOUT, NMAX, CNTRLS, LVALUE, IVALUE,
     1                RVALUE, KERR)
C
         DO 10 K = 1, KKGAS
            REAC(K)  = 0.
            XINTM(K) = 0.
            PROD(K)  = 0.
            SINJ(K)  = 0.
   10    CONTINUE
C
         IF ( KKSURF .GT. 0) THEN
            DO 15 K = 1, KKSURF
               Z(K) = 0.0
   15       CONTINUE
         ENDIF
C
         DO 20 K = 1, KKTOT
            ACT(K) = 0.0
   20    CONTINUE
C
         IF ( KKBULK .GT. 0) THEN
            DO 25 K = 1, KKBULK
               BUL(K)   = 0.0
               BULIN(K) = 0.0
   25       CONTINUE
         ENDIF
C
         DO 30 N = 1, NPHASE
            IETCH(N) = 0
   30    CONTINUE
C
         CALL SKSDEN (ISKWRK, RSKWRK, SDEN)
         CALL SKSDEN (ISKWRK, RSKWRK, SDEN0)
C
         FINJ  = 0.0
         TINJ  = 0.0
         WINJ  = 0.0
         XINJ  = 0.0
C
         NREAC = 0
         NINTM = 0
         NPROD = 0
         NPTS  = 6
         NTOT  = NMAX
         SFLR  = -1.E-4
         SPOS  = -1.0
         N1CALL = 2
         MFILE = 1
         NTEMP = 0
         NEMS  = 0
         NEMSG  = 0
         AOMEGA = 0.0
         ASWIRL = 0.0
         ASPREAD = 0.0
         NP = 0
         WNDFAC = 1.0
         XSTR   = 0.0
         GFAC = 1.0
         SFAC = 1.0
         QDOT = 0.0
         WSRC = 0.0
         XSRC = 0.0
         TDISK = 1.0
         TWAL = 500.0
         EMIS = 0.85
         ERAD = 0.85
         TRAD = 1000.0
         VFAC = 0.0
         POWR = 0.0
         PCTADP = 0.75
         JJREGD = 40
         RATGTC = 1.0
         COND = 1.38
         CNDX = 0.0
         CNDT = 300.0
         FLMX = 0.0
         FLMT = 1.0
         FLMW = 0.025
C
         LFIRST = .TRUE.
         LFRSTI = .TRUE.
         LFRSTG = .TRUE.
         LGDOT  = .FALSE.
         LUMESH = .TRUE.
         LUSTGV = .FALSE.
         LCNTUE = .FALSE.
         LRSTRT = .FALSE.
         LTDIF  = .FALSE.
         LVCOR  = .FALSE.
         LTIME2 = .FALSE.
         LENRGY = .FALSE.
         LASEN  = .FALSE.
         LHSEN  = .FALSE.
         LMULTI = .FALSE.
         LUINF  = .FALSE.
         UINF   = 0.0
         LSTAG  = .FALSE.
         LCHEM  = .TRUE.
         LREORD = .FALSE.
         LPRNT  = .FALSE.
         LNONR  = .FALSE.
         LRADB  = .FALSE.
         LREGRD = .FALSE.
         LCNDCT = .FALSE.
         LSDEN  = .TRUE.
         LCOMP  = .FALSE.
         LTOFF  = .FALSE.
         LPUTIL = .FALSE.
         LFLAME = .FALSE.
         LOLDSP = .FALSE.
       ENDIF
C
C--------------------------------------------------------------
C
      WRITE (LOUT, '(/A/)') '           KEYWORD INPUT '
C
C     Read next input line
90    CONTINUE
      KEYWRD = ' '
      LINE = ' '
      IERR = .FALSE.
      READ  (LIN,  '(A)', END=500, ERR=500) LINE
      WRITE (LOUT, '(10X,A)') LINE
      KEYWRD = CKCHUP(LINE(1:4), 4)
      LINE(1:4) = ' '
C
C     Is this a keyword comment?
C
      IF (KEYWRD(1:1) .EQ. '.' .OR. KEYWRD(1:1) .EQ. '/' .OR.
     1    KEYWRD(1:1) .EQ. '!') GO TO 90
C
C
C-----PROBLEM TYPE KEYWORDS--------------------
C
      IF (KEYWRD .EQ. 'TGIV') THEN
C        Energy equation is not included
         NEC(1)   = .TRUE.
         LENRGY   = .FALSE.
         LVALUE(QADAPT,2) = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'ENRG') THEN
C        Energy equation is included
         NEC(1) = .TRUE.
         LENRGY = .TRUE.
         LVALUE(QADAPT,3) = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'STAG') THEN
C        Stagnation-point flow
         NEC(6) = .TRUE.
         LSTAG  = .TRUE.
C
C-----METHOD OPTIONS KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'ATOL') THEN
C        Absolute Newton iteration convergence criteria
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RVALUE(QSSABS,1), IERR)
         RVALUE(QSSABS,2) = RVALUE(QSSABS,1)
         RVALUE(QSSABS,3) = RVALUE(QSSABS,1)
C
      ELSEIF (KEYWRD .EQ. 'NJAC') THEN
C        Retirement age of Jacobian during steady-state Newton
        CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
        IVALUE(QSSAGE,1) =  INT (VALUE(1))
        IVALUE(QSSAGE,2) = IVALUE(QSSAGE,1)
        IVALUE(QSSAGE,3) = IVALUE(QSSAGE,1)
C
      ELSEIF (KEYWRD .EQ. 'RTOL') THEN
C        Relative Newton iteration convergence criteria
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RVALUE(QSSREL,1), IERR)
         RVALUE(QSSREL,2) = RVALUE(QSSREL,1)
         RVALUE(QSSREL,3) = RVALUE(QSSREL,1)
C
      ELSEIF (KEYWRD .EQ. 'ATIM') THEN
C        Absolute Newton convergence criteria for timesteps
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RVALUE(QTDABS,1), IERR)
         RVALUE(QTDABS,2) = RVALUE(QTDABS,1)
         RVALUE(QTDABS,3) = RVALUE(QTDABS,1)
C
      ELSEIF (KEYWRD .EQ. 'RTIM') THEN
C        Relative Newton convergence criteria for timesteps
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RVALUE(QTDREL,1), IERR)
         RVALUE(QTDREL,2) = RVALUE(QTDREL,1)
         RVALUE(QTDREL,3) = RVALUE(QTDREL,1)
C
      ELSEIF (KEYWRD .EQ. 'TIME') THEN
C        Time step starting procedure
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         IVALUE(QSTPS1,1) = INT(VALUE(1))
         IVALUE(QSTPS1,2) = IVALUE(QSTPS1,1)
         RVALUE(QSTRID,1) = VALUE(2)
         RVALUE(QSTRID,2) = RVALUE(QSTRID,1)
C
      ELSEIF (KEYWRD .EQ. 'ISTP') THEN
C        Number of initial time steps before Newton
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IVALUE(QSTPS0,2) = ABS (INT(VALUE(1)))
         IVALUE(QSTPS0,3) = IVALUE(QSTPS0,2)
         IF (VALUE(1) .LT. 0) LVALUE(QSTDY,2) = .FALSE.
         LVALUE(QSTDY,3) = LVALUE(QSTDY,2)
C
      ELSEIF (KEYWRD .EQ. 'IRET') THEN
C        Retirement age of old time step (default-50)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IVALUE(QSTPS2,2) = INT(VALUE(1))
         IVALUE(QSTPS2,3) = IVALUE(QSTPS2,2)
C
      ELSEIF (KEYWRD .EQ. 'TIM2') THEN
C        Time stepping, after adding the energy equation
         LTIME2    = .TRUE.
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         IVALUE(QSTPS1,3) = INT(VALUE(1))
         RVALUE(QSTRID,3) = VALUE(2)
C
      ELSEIF (KEYWRD .EQ. 'UFAC') THEN
C        Timestep increase when timestep does not change solution
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RVALUE(QTINC,1), IERR)
         RVALUE(QTINC,2) = RVALUE(QTINC,1)
         RVALUE(QTINC,3) = RVALUE(QTINC,1)
C
      ELSEIF (KEYWRD .EQ. 'DFAC') THEN
C       Timestep decrease when Newton fails convergence on timestep
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RVALUE(QTDEC,1), IERR)
         RVALUE(QTDEC,2) = RVALUE(QTDEC,1)
         RVALUE(QTDEC,3) = RVALUE(QTDEC,1)
C
      ELSEIF (KEYWRD .EQ. 'DTMN') THEN
C        minimum time step
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RVALUE(QTMIN,1), IERR)
         RVALUE(QTMIN,2) = RVALUE(QTMIN,1)
         RVALUE(QTMIN,3) = RVALUE(QTMIN,1)
C
      ELSEIF (KEYWRD .EQ. 'DTMX') THEN
C        Maximum time step
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RVALUE(QTMAX,1), IERR)
         RVALUE(QTMAX,2) = RVALUE(QTMAX,1)
         RVALUE(QTMAX,3) = RVALUE(QTMAX,1)
C
      ELSEIF (KEYWRD .EQ. 'TJAC') THEN
C        Retirement age of Jacobian during time stepping
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IVALUE(QTDAGE,1) = INT (VALUE(1))
         IVALUE(QTDAGE,2) = IVALUE(QTDAGE,1)
         IVALUE(QTDAGE,3) = IVALUE(QTDAGE,1)
C
      ELSEIF (KEYWRD .EQ. 'NADP') THEN
C        Number of points that can be added during an adaption
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IVALUE(QPADD,1) = INT(VALUE(1))
         IVALUE(QPADD,2) = IVALUE(QPADD,1)
         IVALUE(QPADD,3) = IVALUE(QPADD,1)
C
      ELSEIF (KEYWRD .EQ. 'TOFF') THEN
C        Turn temperature inactive for regridding criteria
         LTOFF = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'NOFT') THEN
C        Do not do the fixed temperature problem
         N1CALL = 3
C
      ELSEIF (KEYWRD .EQ. 'SFLR') THEN
C        Floor value for the species bounds
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SFLR, IERR)
C
      ELSEIF (KEYWRD .EQ. 'SPOS') THEN
C        convert negative gas and site species solutions
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SPOS, IERR)
C
      ELSEIF (KEYWRD .EQ. 'NPTS') THEN
C        number of initial mesh points (this is overwritten GRID input)
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         NPTS   = INT(VALUE(1))
      ELSEIF (KEYWRD .EQ. 'NMAX') THEN
C        Maximum number of mesh points
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IF (VALUE(1) .GT. NMAX) THEN
            WRITE (LOUT, *) ' NMAX IGNORED, MUST BE < OR = ',NMAX
         ELSEIF (VALUE(1) .LT. NP) THEN
            WRITE (LOUT, *) ' NMAX IGNORED, ALREADY HAVE ', NP,' POINTS'
         ELSE
            NTOT = INT(VALUE(1))
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'GRID') THEN
C        initial mesh
C
         LUMESH    = .FALSE.
         IERR = (NP+1 .GT. NMAX)
         IF (.NOT. IERR) THEN
            NP = NP + 1
            CALL CKXNUM (LINE, 1, LOUT, NVAL, X(NP), IERR)
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'GRAD') THEN
C        gradient mesh adaption parameter
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RVALUE(QTOL1,1), IERR)
         RVALUE(QTOL1,2) = RVALUE(QTOL1,1)
         RVALUE(QTOL1,3) = RVALUE(QTOL1,1)
C
      ELSEIF (KEYWRD .EQ. 'CURV') THEN
C        curvature mesh adaption parameter
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RVALUE(QTOL2,1), IERR)
         RVALUE(QTOL2,2) = RVALUE(QTOL2,1)
         RVALUE(QTOL2,3) = RVALUE(QTOL2,2)
C
      ELSEIF (KEYWRD .EQ. 'XSTR') THEN
C        point for left boundary condition
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XSTR, IERR)
C
      ELSEIF (KEYWRD .EQ. 'XCEN') THEN
C        center of mixing region
         NOPT(2) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XCEN, IERR)
C
      ELSEIF (KEYWRD .EQ. 'XEND') THEN
C         distance at which end boundary condition is applied
         NEC(2)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XEND, IERR)
C
      ELSEIF (KEYWRD .EQ. 'WMIX') THEN
C        width of mixing zone
         NOPT(1) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, WMIX, IERR)
C
C-----REACTOR DEFINITION KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'PRES') THEN
C         pressure
         NEC(3)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, P, IERR)
         P = P*PATM
C
      ELSEIF (KEYWRD .EQ. 'TINF') THEN
C        temperature at infinity
         NEC(4)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TINFTY, IERR)
C
      ELSEIF (KEYWRD .EQ. 'TDSK') THEN
C        rotating disk temperature (Kelvin)
         NEC(5)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TDISK, IERR)
C
      ELSEIF (KEYWRD .EQ. 'OMEG') THEN
C        disk spin rate
         NEC(6)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, AOMEGA, IERR)
         AOMEGA = AOMEGA * 2.0 * 3.1415926535 / 60.0
C
      ELSEIF (KEYWRD .EQ. 'AINL') THEN
C        inlet swirl rate
         CALL CKXNUM (LINE, 1, LOUT, NVAL, ASPREAD, IERR)
C
      ELSEIF (KEYWRD .EQ. 'OINL') THEN
C        inlet swirl rate
         CALL CKXNUM (LINE, 1, LOUT, NVAL, ASWIRL, IERR)
         ASWIRL = ASWIRL * 2.0 * 3.1415926535 / 60.0
C
      ELSEIF (KEYWRD .EQ. 'UINF') THEN
C        axial velocity at infinity
         CALL CKXNUM (LINE, 1, LOUT, NVAL, UINF, IERR)
         LUINF = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'TEMP') THEN
C        read specified temperature profile (X,T) pairs
C
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         IF (NTEMP+1 .GT. NMAX) THEN
            WRITE (LOUT, *)
     1      ' ERROR... THE PROBLEM IS ONLY DIMENSIONED FOR ', NMAX,
     2      ' (X,T) PAIRS'
            IERR = .TRUE.
         ELSE
            NTEMP = NTEMP+1
            XGIVEN(NTEMP) = VALUE(1)
            TGIVEN(NTEMP) = VALUE(2)
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'USTG') THEN
C        on a restart use given temperature profile, not the one on
C        the restart file
         LUSTGV = .TRUE.
C
C
C-----RADIATION BC OPTIONS KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'RADB') THEN
C        calculate initial TDISK
         LRADB = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'TWAL') THEN
C        ambient wall temperature for radiation balance
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TWAL, IERR)
C
      ELSEIF (KEYWRD .EQ. 'EMIS') THEN
C        emissivity of the substrate
         CALL CKXNUM (LINE, 1, LOUT, NVAL, EMIS, IERR)
C
      ELSEIF (KEYWRD .EQ. 'ERAD') THEN
C        emissivity of a radiating disk above the substrate
         CALL CKXNUM (LINE, 1, LOUT, NVAL, ERAD, IERR)
C
      ELSEIF (KEYWRD .EQ. 'RDSK') THEN
C        ratio of substrate radius to sep. distance to upper disk
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RDSK, IERR)
C
      ELSEIF (KEYWRD .EQ. 'RRAD') THEN
C        ratio of upper disk radius to sep. distance to lower disk
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RRAD, IERR)
C
      ELSEIF (KEYWRD .EQ. 'TRAD') THEN
C        temperature of upper ratiating disk
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TRAD, IERR)
C
      ELSEIF (KEYWRD .EQ. 'POWR') THEN
C        power heating the disk (ergs/cm**2-sec)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, POWR, IERR)
C
      ELSEIF (KEYWRD .EQ. 'EMSG') THEN
C        gas emissivity with temperature dependence
C
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         IF (NEMSG+1 .GT. NMAX) THEN
            WRITE (LOUT, *)
     1      ' ERROR... THE PROBLEM IS ONLY DIMENSIONED FOR ', NMAX,
     2      ' (T,emiss) PAIRS'
            IERR = .TRUE.
         ELSE
            NEMSG = NEMSG+1
            TEMISG(NEMSG) = VALUE(1)
            EMISG(NEMSG) = VALUE(2)
         ENDIF
C
C
C-----ENERGY SOURCE TERM KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'QDOT') THEN
C        energy source term
         CALL CKXNUM (LINE, 1, LOUT, NVAL, QDOT, IERR)
         NEC(9) = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'XSRC') THEN
C        X coordinate for center of energy source
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XSRC, IERR)
         NEC(10) = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'WSRC') THEN
C        base width in X of Gaussian source term
         CALL CKXNUM (LINE, 1, LOUT, NVAL, WSRC, IERR)
C
C-----CHEMISTRY KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'NOCH') THEN
C        do not solve gas-phase chemistry
         LCHEM = .FALSE.
C
      ELSEIF (KEYWRD .EQ. 'GFAC') THEN
C        factor for gas-phase rate constants
         CALL CKXNUM (LINE, 1, LOUT, NVAL, GFAC, IERR)
C
      ELSEIF (KEYWRD .EQ. 'SFAC') THEN
C        factor for surface-phase rate constants
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SFAC, IERR)
C
      ELSEIF (KEYWRD .EQ. 'CHEM') THEN
C        do solve gas-phase chemistry
         LCHEM = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'REAC') THEN
C        reactant
C
         IF (LFIRST) THEN
            LFIRST = .FALSE.
            NREAC = 0
            DO 100 K = 1, KKGAS
               REAC(K) = 0.
  100       CONTINUE
         ENDIF
         CALL CKSNUM (LINE, 1, LOUT, KNAM, KKGAS, KSP, NVAL, VALUE,
     1                IERR)
         IERR = IERR.OR. NREAC+1.GT.KKGAS
         IF (.NOT. IERR) THEN
            NREAC       = NREAC+1
            KR(NREAC)   = KSP
            REAC(KSP)   = VALUE(1)
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'COMP') THEN
C        inlet mass fraction
         LCOMP = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'FLUX') THEN
C        inlet mass flux fraction
         LCOMP = .FALSE.
C
      ELSEIF (KEYWRD .EQ. 'INTM') THEN
C        intermediate
C
         CALL CKSNUM (LINE, 1, LOUT, KNAM, KKGAS, KSP, NVAL, VALUE,
     1                IERR)
         IERR = IERR.OR. NINTM+1.GT.KKGAS
         IF (.NOT. IERR) THEN
            NINTM       = NINTM + 1
            KI(NINTM)   = KSP
            XINTM(KSP)  = VALUE(1)
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'PROD') THEN
C        product
C
         CALL CKSNUM (LINE, 1, LOUT, KNAM, KKGAS, KSP, NVAL, VALUE,
     1                IERR)
         IERR = IERR.OR. NPROD+1.GT.KKGAS
         IF (.NOT. IERR) THEN
            NPROD = NPROD + 1
            KP(NPROD) = KSP
            PROD(KSP) = VALUE(1)
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'SURF') THEN
C        surface species initial guess
C
         IF (KKSURF .LE. 0) THEN
            IERR = .TRUE.
            WRITE (LOUT, *) ' Error...no site-phase species exist '
         ELSE
            CALL SKSNUM (LINE, 1, LOUT, KNAM, KKTOT, PNAM, NPHASE,
     1                   KKPHAS, KSP, NKF, NVAL, VALUE, IERR)
            IF (IERR) THEN
            ELSEIF (KSP.LT.KFIRST(NFSURF).OR.KSP.GT.KLAST(NLSURF)) THEN
               WRITE (LOUT, *)
     1         ' Error...SURF must be site-phase species '
               IERR = .TRUE.
            ELSE
               IF (NKF .GT. 1) WRITE (LOUT, *)
     1         ' Warning...non-unique species name given '
               Z(KSP-KKGAS) = VALUE(1)
            ENDIF
         ENDIF
C
      ELSEIF (KEYWRD(:3) .EQ. 'ACT') THEN
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
               WRITE (LOUT, *) ' Error...ACT must be bulk-phase species'
               IERR = .TRUE.
            ELSE
               IF (NKF .GT. 1) WRITE (LOUT, *)
     1         ' Warning...non-unique species name given '
               ACT(KSP) = VALUE(1)
            ENDIF
        ENDIF
C
      ELSEIF (KEYWRD .EQ. 'SDEN') THEN
C        surface site densities
C
         IF (NNSURF .LE. 0) THEN
            IERR = .TRUE.
            WRITE (LOUT, *) ' Error...no surface phases exist'
         ELSE
            CALL CKSNUM (LINE, 1, LOUT, PNAM, NPHASE, IPHASE,
     1                   NVAL, VALUE, IERR)
           IF (IPHASE .LT. NFSURF .OR. IPHASE .GT. NLSURF) IERR = .TRUE.
           IF (IERR) THEN
               WRITE (LOUT,'(A)')
     1                  ' ERROR READING DATA FOR KEYWORD '//KEYWRD
           ELSE
              SDEN (IPHASE) = VALUE(1)
              SDEN0(IPHASE) = VALUE(1)
           ENDIF
        ENDIF
C
      ELSEIF (KEYWRD .EQ. 'NSDN') THEN
C        do not solve an equation for surface site densities, no matter
C        what the mechanism
         LSDEN = .FALSE.
C
      ELSEIF (KEYWRD .EQ. 'GDOT') THEN
C        destruction rate of gas-phase species due to surface reactions
C
         IF (LFRSTG) THEN
            LFRSTG = .FALSE.
            LGDOT  = .TRUE.
            NGDOT  = 0
            DO 1050 K = 1, KKGAS
               GDOT(K) = 0.
 1050       CONTINUE
         ENDIF
         CALL CKSNUM (LINE, 1, LOUT, KNAM, KKGAS, KSP, NVAL, VALUE,
     1                IERR)
         IERR = IERR.OR. NGDOT+1.GT.KKGAS
         IF (.NOT. IERR) THEN
            NGDOT = NGDOT + 1
            GDOT(KSP) = VALUE(1)
         ENDIF
C
       ELSEIF (KEYWRD .EQ. 'ETCH') THEN
C         specify that a bulk phase is to be etched instead of grown
C
         CALL CKCRAY (LINE, NPHASE, PNAM, LOUT, 5, IVAL, NVAL, IERR)
C         IF (NVAL .GT. 5) KERR = .TRUE.
         IF (NVAL .GT. 0) THEN
           IF (NVAL .GT. 5) KERR = .TRUE.
           CALL CKSNUM (LINE, -1 ,LOUT, PNAM, NPHASE, IPHASE, NVAL,
     1                  VALUE, IERR)
           IF (IPHASE .LT. NFBULK .OR. IPHASE .GT. NLBULK) IERR = .TRUE.
           IF (IERR) THEN
               WRITE (LOUT,'(A)')
     1                  ' ERROR READING DATA FOR KEYWORD '//KEYWRD
           ELSE
              IETCH(IPHASE) = 1
           ENDIF
         ELSE
           IF (NNBULK.GT.0) THEN
              DO 120 IPHASE = NFBULK, NLBULK
                IETCH(IPHASE) = 1
  120         CONTINUE
           ELSE
              IERR = .TRUE.
              WRITE (LOUT,'(A)')
     1      ' ETCH KEYWORD USED WHEN MECHANISM CONTAINS NO BULK PHASES'
           ENDIF
         ENDIF
C
C-----INJECTION KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'INJS') THEN
C        injected species
C
         IF (LFRSTI) THEN
            LFRSTI = .FALSE.
            NINJS = 0
            DO 1200 K = 1, KKGAS
               SINJ(K) = 0.
1200        CONTINUE
         ENDIF
         CALL CKSNUM (LINE, 1, LOUT, KNAM, KKGAS, KSP, NVAL, VALUE,
     1                IERR)
         IERR = IERR.OR. NINJS+1.GT.KKGAS
         IF (.NOT. IERR) THEN
            NINJS     = NINJS+1
            SINJ(KSP) = VALUE(1)
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'INJX') THEN
C        location of height of injection
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XINJ, IERR)
C
      ELSEIF (KEYWRD .EQ. 'INJW') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, WINJ, IERR)
C
      ELSEIF (KEYWRD .EQ. 'INJM') THEN
C        flow rate of injection
         CALL CKXNUM (LINE, 1, LOUT, NVAL, FINJ, IERR)
C
      ELSEIF (KEYWRD .EQ. 'INJT') THEN
C        temperature of injection
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TINJ, IERR)
C
C-----FLAME RESTART KEYWORDS---------------------------
C
      ELSEIF (KEYWRD .EQ. 'FLAM') THEN
C        position and value of fixed temperature in flame
         LFLAME = .TRUE.
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         FLMX = VALUE(1)
         FLMT = VALUE(2)
C
      ELSEIF (KEYWRD .EQ. 'FLWD') THEN
C        flame width (optional input)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, FLMW, IERR)
C
C-----SUBSTRATE CONDUCTION KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'CDCT') THEN
C        include heat conduction through substrate
         LCNDCT   = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'COND') THEN
C        thermal conductivity of substrate (w/cm K)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, COND, IERR)
C
      ELSEIF (KEYWRD .EQ. 'CNDX') THEN
C        conduction path length, or substrate thickness (cm)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, CNDX, IERR)
C
      ELSEIF (KEYWRD .EQ. 'CNDT') THEN
C        temperature on far side of substrate (K)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, CNDT, IERR)
C
C-----TRANSPORT OPTIONS KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'MULT') THEN
C        multicomponent diffusion included
         LMULTI   = .TRUE.
C
      ELSEIF (KEYWRD(:3) .EQ. 'MIX') THEN
C        mixture-averaged diffusion used
         LMULTI   = .FALSE.
C
      ELSEIF (KEYWRD .EQ. 'TDIF') THEN
C        thermal diffusion included
         LTDIF    = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'VCOR') THEN
C        use correction velocity formalism
         LVCOR    = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'TRCE') THEN
C        use "trace" approximation, lump all the transport errors
C        into the last species
         LVCOR    = .FALSE.
C
      ELSEIF (KEYWRD .EQ. 'NONR') THEN
C        do not solve the non-reacting problem
         LNONR = .TRUE.
C
C-----SENSITIVITY KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'ASEN') THEN
C        all reaction snsitivity
         LASEN = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'HSEN') THEN
C        H298 sensitivity
         LHSEN = .TRUE.
C
C-----REORDER
C
      ELSEIF (KEYWRD .EQ. 'REOR') THEN
          LREORD = .TRUE.
C
C-----PRINTING AND RESTARTING KEYWORDS--------------------
C
      ELSEIF (KEYWRD .EQ. 'PRNT') THEN
C        print control
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IPRNT = INT(VALUE(1))
         DO 400 IC = 1,3
            IF (IPRNT .LT. 10) THEN
               IVALUE(QLEVM,IC) = MIN (3, IPRNT)
               IVALUE(QLEVD,IC) = IVALUE(QLEVM,IC)
            ELSE
               IVALUE(QLEVM,IC) = IPRNT/10
               IVALUE(QLEVD,IC) = IPRNT - 10*IVALUE(QLEVM,IC)
               IVALUE(QLEVM,IC) =
     1                     MAX (IVALUE(QLEVM,IC), IVALUE(QLEVD,IC))
            ENDIF
400      CONTINUE
C
      ELSEIF (KEYWRD .EQ. 'LPRT') THEN
         LPRNT = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'SPRT') THEN
         LPRNT = .FALSE.
C
      ELSEIF (KEYWRD .EQ. 'PUTL') THEN
         LPUTIL = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'RSTR') THEN
C        restart check
         LRSTRT = .TRUE.
         RSTCNT = .FALSE.
C
      ELSEIF (KEYWRD .EQ. 'FILE') THEN
C        restart files
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         MFILE    = INT(VALUE(1))
C
      ELSEIF (KEYWRD .EQ. 'OLDS') THEN
C        restart file is from old SPIN, before ALE vers
         LOLDSP = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'CNTN') THEN
C        continuation flag
         LCNTUE   = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'JJRG') THEN
C        number of mesh points in the regrid
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         JJREGD = INT(VALUE(1))
         LREGRD = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'PCAD') THEN
C        percentage of regrid points dedicated to adaption
         CALL CKXNUM (LINE, 1, LOUT, NVAL, PCTADP, IERR)
C
      ELSEIF (KEYWRD .EQ. 'RGTC') THEN
C        ratio of gradient regrid points to curvature points
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RATGTC, IERR)
C
      ELSEIF (KEYWRD(:3) .EQ. 'END') THEN
C        last line
         KERR = KERR.OR.IERR
         GO TO 500
C
      ELSE
C
C        to get here, an invalid keyword was reac
         WRITE (LOUT, *) ' ERROR...ILLEGAL KEYWORD'
         KERR = .TRUE.
      ENDIF
C
C     bo gack up and read the next line
      KERR = KERR.OR.IERR
      GO TO 90
C
C-----END OF KEYWORDS--------------------
C
C     CHECK THE REACTANT, PRODUCT, INJECTION SUMS
C
  500 CONTINUE
C
      SUMR       = 0.
      SUMP       = 0.
      SUMI       = 0.
      DO 510 K = 1, KKGAS
         SUMR = SUMR+REAC(K)
         SUMP = SUMP+PROD(K)
         SUMI = SUMI+SINJ(K)
  510 CONTINUE
C
C     normalize reactant and product fractions
      IF (SUMR .GT. 0.0) THEN
         DO 520 K = 1, KKGAS
            REAC(K) = REAC(K)/SUMR
  520    CONTINUE
      ELSE
         WRITE (LOUT, *) ' ERROR..."REAC" NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF (SUMP .GT. 0.0) THEN
         DO 530 K = 1, KKGAS
            PROD(K) = PROD(K)/SUMP
  530    CONTINUE
      ELSE
         WRITE (LOUT, *) ' WARNING..."PROD" NOT SPECIFIED'
         IF (SUMI .GT. 0.0) THEN
           NPROD = KKGAS
           DO 532 K = 1, KKGAS
             KP(K) = K
             PROD(K) = SINJ(K) / SUMI
  532      CONTINUE
         ELSE
           NPROD = KKGAS
           DO 534 K = 1, KKGAS
             KP(K) = K
             PROD(K) = (9.*REAC(K) + 1.0/FLOAT(KKGAS)) / 10.
  534      CONTINUE
         ENDIF
      ENDIF
C
      IF (SUMI .GT. 0.0) THEN
         DO 540 K = 1, KKGAS
            SINJ(K) = SINJ(K)/SUMI
  540    CONTINUE
      ENDIF
C
      IF (ABS(SUMR-1.0) .GT. 1.E-6)
     1 WRITE (LOUT, *) ' CAUTION...REACTANT FRACTIONS SUM TO ', SUMR
C
      IF ((.NOT. CNTNUD) .AND. ABS(SUMP-1.0).GT.1.E-6)
     1 WRITE (LOUT, *) ' CAUTION...PRODUCT FRACTIONS SUM TO ',  SUMP
C
C     normalize the surface site fraction sums
      IF (KKSURF .GT. 0) THEN
C
         DO 600 N = NFSURF, NLSURF
            SUMZ = 0.0
            DO 550 K = KFIRST(N), KLAST(N)
               SUMZ = SUMZ + Z(K-KKGAS)
  550       CONTINUE
C
            IF (SUMZ .GT. 0.0) THEN
C
               DO 560 K = KFIRST(N), KLAST(N)
                  Z(K-KKGAS) = Z(K-KKGAS) / SUMZ
  560          CONTINUE
C
               IF (ABS(SUMZ-1.0) .GT. 1.E-6)
     1             WRITE (LOUT, 8305) SUMZ, PNAM(N)
C
            ELSE
C
               SPECFR = 1.0 / FLOAT( KLAST(N) - KFIRST(N) + 1 )
               DO 570 K = KFIRST(N), KLAST(N)
                  Z(K-KKGAS) = SPECFR
  570          CONTINUE
               WRITE (LOUT, 8306) PNAM(N)
C
            ENDIF
C
  600    CONTINUE
C
      ENDIF
C
      IF (KKBULK .GT. 0) THEN
C         normalize bulk activities
C
         DO 650 N = NFBULK, NLBULK
            SUMA = 0.0
            DO 620 K = KFIRST(N), KLAST(N)
               SUMA = SUMA + ACT(K)
  620       CONTINUE
C
            IF (SUMA .GT. 0) THEN
               DO 630 K = KFIRST(N), KLAST(N)
                  ACT(K) = ACT(K) / SUMA
                  BUL  (K-KKGAS-KKSURF) = ACT(K)
                  BULIN(K-KKGAS-KKSURF) = ACT(K)
  630          CONTINUE
            ENDIF
C
            IF (ABS(SUMA-1.0) .GT. 1.E-6) WRITE (LOUT, '(1X,A,E13.5)')
     1      ' CAUTION...BULK ACTIVITIES SUM TO ',SUMA
  650    CONTINUE
      ENDIF
C
C     CHECK FOR NECESSARY INPUT
C
      IF (.NOT. NEC(1)) THEN
         WRITE (LOUT, *)
     1                ' ERROR...EITHER "TGIV" OR "ENRG" NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(2)) THEN
         WRITE (LOUT, *) ' ERROR..."XEND" NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(3)) THEN
         WRITE (LOUT, *) ' ERROR...PRESSURE NOT GIVEN'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(4))THEN
         WRITE (LOUT, *) ' ERROR..."TINF" NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF ((.NOT.NEC(5)) .AND. (KKBULK+KKSURF .EQ. 0)) THEN
         WRITE (LOUT, *)
     1   ' ERROR..."TDSK" NOT SPECIFIED AND NO SURFACE SPECIES PRESENT'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(6)) THEN
         WRITE (LOUT, *)
     1      ' ERROR...EITHER "STAG" NOT SPECIFIED OR "OMEG" NOT GIVEN'
         KERR = .TRUE.
      ENDIF
C
      IF (NEC(9)) THEN
         IF (.NOT.NEC(10)) THEN
            WRITE (LOUT,*)
     1       'ERROR..."QDOT" WAS SPECIFIED BUT "XSRC" NOT GIVEN'
            KERR = .TRUE.
         ENDIF
      ELSEIF (NEC(10)) THEN
         IF (.NOT.NEC(9)) THEN
            WRITE (LOUT,*)
     1       'ERROR..."XSRC" WAS SPECIFIED WITH NO "QDOT"'
            KERR = .TRUE.
         ENDIF
      ENDIF
C
C     Check to make sure that the temperatures have been given
C
      IF (NTEMP .EQ. 0) THEN
            WRITE (LOUT, *)
     1              ' WARNING...TEMPERATURE PROFILE WASN''T GIVEN'
         NTEMP = 2
         XGIVEN(1) = XSTR
         TGIVEN(1) = TINFTY
         XGIVEN(2) = XEND
         TGIVEN(2) = TDISK
      ENDIF
C
C     Make sure the (X,T) pairs are in order
C
      DO 700 N = 2, NTEMP
         IF (XGIVEN(N-1) .GE. XGIVEN(N)) THEN
            WRITE (LOUT, *)
     1              ' ERROR...SPECIFIED TEMPERATURES ARE OUT OF ORDER'
            KERR = .TRUE.
         ENDIF
  700 CONTINUE
C
C     Make sure specified emissivities are in ascending T order
C
      DO 720 N = 2, NEMSG
         IF (TEMISG(N-1) .GE. TEMISG(N)) THEN
            WRITE (LOUT,*)
     1              ' ERROR...SPECIFIED EMISSIVITIES ARE OUT OF ORDER'
            KERR = .TRUE.
         ENDIF
  720 CONTINUE
      IF (NEMSG .NE. 0) NEMS = NEMSG
C
C     Make sure the initial grid points are in order
C
      IF (.NOT.CNTNUD .AND. .NOT.LUMESH) THEN
         NPTS = NP
         DO 740 N = 2, NPTS
            IF (X(N-1) .GE. X(N)) THEN
              WRITE (LOUT, *)
     1                 ' ERROR...INITIAL GRID IS OUT OF ORDER'
              KERR = .TRUE.
            ENDIF
  740    CONTINUE
      ENDIF
C
C     Make sure the given temperatures span the XEND-XSTR domain
C
      IF (.NOT.LRSTRT .OR. .NOT.CNTNUD .OR. LUSTGV) THEN
         IF (XGIVEN(1).GT.XSTR .OR. XGIVEN(NTEMP).LT.XEND) THEN
            WRITE (LOUT, *)
     1      ' ERROR...GIVEN TEMPERATURE PROFILE DOES NOT SPAN XEND-XSTR'
            KERR = .TRUE.
         ENDIF
      ENDIF
C
C     Check for consistency if have radiating disk above substrate
C
      IF ( LRADB ) THEN
         IF (RRAD*RDSK .EQ. 0.0 .AND. RRAD+RDSK .NE. 0.0 ) THEN
            WRITE (LOUT, *)
     +          ' ERROR...TO HAVE A RADIATING DISK ABOVE THE SUBSTRATE,'
            WRITE (LOUT, *)
     +          '        BOTH RDSK AND RRAD MUST BE SPECIFIED (> 0).'
            KERR = .TRUE.
         ENDIF
C
         IF (RRAD*RDSK .GT. 0.0 ) THEN
            XTEMP = 1.0 + (1.0 + RRAD**2)/RDSK**2
            VFAC = 0.5 * (XTEMP
     +           - SQRT( XTEMP**2 - 4.0*(RDSK/RRAD)**2 ) )
         ENDIF
      ENDIF
C
C     Make sure have substrate thickness if including conduction
C
      IF ( LCNDCT ) THEN
         IF ( CNDX .EQ. 0.0 ) THEN
            WRITE (LOUT, *)
     +          ' ERROR...SUBSTRATE THICKNESS NOT SPECIFIED.'
            KERR = .TRUE.
         ELSE
            CNDFAC = 1.0E07 * COND / CNDX
         ENDIF
      ENDIF
C
C     Set optional input if needed
C
      IF (.NOT. NOPT(1)) WMIX = (XEND-XSTR)*0.50
      IF (.NOT. NOPT(2)) XCEN = (XEND-XSTR)*0.35
      IF (.NOT. LTIME2) THEN
         IVALUE(QSTPS1,3) = IVALUE(QSTPS1,2)
         RVALUE(QSTRID,3) = RVALUE(QSTRID,2)
      ENDIF
C
C     Check for consistence of "STAG" and AOMEGA
C
      IF (LSTAG) THEN
C        MUST INPUT UINF
         IF (.NOT. LUINF) THEN
            WRITE (LOUT, *)
     1         ' ERROR...MUST SUPPLY "UINF" FOR STAGNATION-POINT FLOW'
            KERR = .TRUE.
         ENDIF
C        CANNOT HAVE NON-ZERO AOMEGA FOR STAGNATION-POINT FLOW
         IF (AOMEGA .NE. 0.0) THEN
            WRITE (LOUT, *)
     1         ' ERROR...NON-ZERO OMEG INCONSISTENT WITH "STAG" OPTION'
            KERR = .TRUE.
         ENDIF
      ELSE
         IF (AOMEGA.EQ.0.0 .AND. .NOT.LUINF) THEN
C           MUST INPUT UINF WHEN DISK IS NOT SPINNING
               WRITE (LOUT, *)
     1             ' ERROR...MUST SUPPLY "UINF" WHEN OMEG=0 '
               KERR = .TRUE.
         ENDIF
      ENDIF
C
C     Check that no surface mechanism was used if GTOT was read
C
      IF (LGDOT .AND. KKSURF.NE.0 .AND. KKBULK.NE.0) THEN
C           WE HAVE DEFINED THIS TO BE AN ERROR
            WRITE (LOUT, *)
     1      ' SURFACE CHEMKIN REACTION MECHANISM MUST BE EMPTY '
            WRITE (LOUT, *)
     1      ' WHEN GDOT KEYWORD IS USED FOR SURFACE BOUNDARY CONDITION '
            KERR = .TRUE.
      ENDIF
C
C     FORMATS
 7000 FORMAT (A4, A)
 8000 FORMAT (10X, A4, A)
 8305 FORMAT (1X,' CAUTION...SITE FRACTIONS SUM TO', F11.8,
     1        ' FOR PHASE ', A16)
 8306 FORMAT (1X,
     1' CAUTION...KEYWORD SURF NOT PRESENT. ASSIGNING EQUAL SITE'/1X,
     2           '           FRACTIONS FOR PHASE ', A16)
C
C     end of SUBROUTINE SPRDKY
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE SPREAD (LREST, LOUT, MFILE, COMPS, POINTS,
     1                   KKTOT, KKSURF, KKBULK, NSSOLV, II, IISUR,
     2                   X, S, LOLDSP, KERR)
C
C  START PROLOGUE
C
C///  READ THE RESTART DATA FILE FOR STARTING PROFILES.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INTEGER COMPS, POINTS
      DIMENSION X(*), S(*)
      CHARACTER*16 ICHR, ICKLNK, IMCLNK, ISKLNK, ISENSI,
     1             IHSENS, IVERSN, VERNUM, SOLTYP
      LOGICAL LOLDSP, KERR
      DATA ISENSI/'SENSITIVITY     '/,
     1     ICKLNK/'CKLINK          '/, IMCLNK/'MCLINK          '/,
     2     ISKLNK/'SKLINK          '/, IHSENS/'HSENSITIVITY    '/,
     3     IVERSN/'VERSION         '/, VERNUM/'1.00            '/,
     4     SOLTYP/'STEADY SOLN SPIN'/
C
      KERR = .FALSE.
      REWIND (LREST)
      NREST = 0
  140 CONTINUE
      READ (LREST, END=250, ERR=250) ICHR
C
      IF (ICHR .EQ. ICKLNK) THEN
         DO 150 L = 1, 4
            READ (LREST, END=250, ERR=250)
  150    CONTINUE
C
      ELSEIF (ICHR .EQ. IMCLNK) THEN
         DO 160 L = 1, 3
            READ (LREST, END=250, ERR=250)
  160    CONTINUE
C
      ELSEIF (ICHR .EQ. ISKLNK) THEN
         DO 170 L = 1, 4
            READ (LREST, END=250, ERR=250)
  170    CONTINUE
C
      ELSEIF (ICHR .EQ. IVERSN) THEN
C        READ THE VERSION NUMBER OF THIS SOLUTION FILE
         READ (LREST, END=250, ERR=250) ICHR
         IF (ICHR .NE. VERNUM) THEN
            WRITE (LOUT, *) 'WRONG VERSION NUMBER ON SOLUTION FILE'
            WRITE (LOUT, *) 'VERSION READ WAS:    ', ICHR
            WRITE (LOUT, *) 'VERSION REQUIRED IS: ', VERNUM
            STOP
         ENDIF
C
      ELSEIF (ICHR .EQ. SOLTYP) THEN
C        HEADER INFO (SHOULD MATCH "STEADY" CHARACTER VARIABLE CONTENTS)
         READ (LREST, END=250, ERR=250) POINTS, NEQ, NNNN, KKGDM, KKSDM,
     1                                  KKBDM, NNSDM
C        SKIPPING OVER THE NEXT 4 RECORDS
         READ (LREST, END=250, ERR=250)
         READ (LREST, END=250, ERR=250)
         READ (LREST, END=250, ERR=250)
         READ (LREST, END=250, ERR=250)
C
         READ (LREST, END=250, ERR=250) (X(J), J=1,POINTS)
         READ (LREST, END=250, ERR=250) (S(N), N=1,NEQ)
         IF (NNNN .NE.COMPS   .OR. KKSDM.NE.KKSURF .OR.
     1       KKBDM.NE.KKBULK .OR. NNSDM.NE.NSSOLV) THEN
            WRITE (LOUT, *) ' FATAL ERROR, INCOMPATIBLE RESTART FILE'
            KERR = .TRUE.
            RETURN
         ENDIF
         IF (LOLDSP) THEN
            NEQOLD = NEQ
            NEQ = NEQOLD + KKBULK
            NB1 = KKSURF + KKBULK
            NB2 = KKSURF + 2*KKBULK
            DO 175 N = NEQ, NB2+1, -1
               S(N) = S(N-KKBULK)
175         CONTINUE
            DO 177 N = NB2, NB1+1
               S(N) = 0.0
177         CONTINUE
         ENDIF
         NREST = NREST + 1
         IF (NREST .EQ. MFILE) GO TO 250
C
      ELSEIF (ICHR .EQ. ISENSI) THEN
         NREAD = II + IISUR
         DO 180 I = 1, NREAD
            READ (LREST, END=250, ERR=250)
  180    CONTINUE
C
      ELSEIF (ICHR .EQ. IHSENS) THEN
         NREAD = KKTOT
         DO 220 K = 1, NREAD
            READ (LREST, END=250, ERR=250)
  220    CONTINUE
C
      ELSE
         WRITE (LOUT, *)
     1      'FATAL ERROR, NOT A SOLUTION ON RESTART FILE'
         KERR = .TRUE.
         RETURN
      ENDIF
C
      GO TO 140
C
  250 CONTINUE
      IF (NREST .NE. MFILE) THEN
         WRITE (LOUT, *) 'NEQ=',NEQ,', NATJ=',NNNN,', JJ=',POINTS
         WRITE (LOUT, *) ' Error reading solution file...'
         KERR = .TRUE.
         RETURN
      ENDIF
C
C     end of SUBROUTINE SPREAD
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE SPRSTR (KKGAS, NMAX, NATJ, JJ, LOUT, LUMESH, LUSTGV,
     1                   NREAC, NINTM, NPROD, REAC, XINTM, PROD,
     2                   TINFTY, TDISK, AOMEGA, ASWIRL, P, KR, KI, KP,
     3                   NTEMP, XGIVEN, TGIVEN, XSTR, XCEN, XEND, WMIX,
     4                   ICKWRK, RCKWRK, Y, SI, EPS, X, S, IMCWRK,
     5                   RMCWRK, RHOINF, ASPREAD, KERR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
C     Integer arrays
      DIMENSION KR(NREAC), KI(NINTM), KP(NPROD), ICKWRK(*),
     1          IMCWRK(*)
C     Real arrays
      DIMENSION REAC(*), XINTM(*), PROD(*), RCKWRK(*), RMCWRK(*),
     2          Y(*), SI(*), EPS(*), XGIVEN(*), TGIVEN(*), X(*),
     3          S(NATJ,*)
C
      LOGICAL LUMESH, LUSTGV, KERR
C
      KERR = .FALSE.
C     initialize mass flux fractions to zero
      DO 100 K = 1, KKGAS
         EPS(K) = 0.0
100   CONTINUE
C
C     shift the solution to the right if a new XSTR is < X(1)
      IF (XSTR .LT. X(1)) THEN
         JJ = JJ + 1
         IF (JJ .GT. NMAX) THEN
            WRITE (LOUT, *) ' ERROR...NEW XSTR NEEDS TOO MANY POINTS'
            KERR = .TRUE.
            RETURN
         ENDIF
         DO 225 I = 2, JJ
            J = JJ + 2 - I
            X(J) = X(J-1)
            DO 200 N = 1, NATJ
               S(N,J) = S(N,J-1)
  200       CONTINUE
  225    CONTINUE
         X(1) = XSTR
         DO 300 N = 1, NATJ
            S(N,1) = S(N,2)
  300    CONTINUE
      ENDIF
C
C     if a new XEND > X(JJ) then add a point at JJ+1,
C     else reduce JJ and set X(JJ) = XEND
C
      IF (XEND .GT. (X(JJ)+1.E-4)) THEN
         JJ = JJ + 1
         IF (JJ .GT. NMAX) THEN
            WRITE (LOUT, *) ' ERROR...NEW XEND NEEDS TOO MANY POINTS'
            KERR = .TRUE.
            RETURN
         ENDIF
         X(JJ) = XEND
         DO 400 N = 1, NATJ
            S(N,JJ) = S(N,JJ-1)
  400    CONTINUE
      ENDIF
C
C     unless 'USTG' is specified, use the restart temperature profiles
C     for the 'given' temperatures
      IF (.NOT. LUSTGV) THEN
         NTEMP = JJ
         DO 500 N = 1, NTEMP
            XGIVEN(N) = X(N)
            TGIVEN(N) = S(NT,N)
  500    CONTINUE
      ENDIF
C
C     set the mass flux fraction boundary conditions
      CALL CKXTY (REAC, ICKWRK, RCKWRK, EPS)
C
C     compute the reference properties at the inlet JJ
      CALL CKRHOY (P, TINFTY, EPS, ICKWRK, RCKWRK, RHOINF)
C
C     end of SUBROUTINE SPRSTR
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE SPSAVE (ICKWRK, RCKWRK, CCKWRK,
     1                   IMCWRK, RMCWRK,  ISKWRK, RSKWRK,
     2                   CSKWRK, LOUT, IUNIT, IVERSN, VERNUM)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ICKWRK(*), RCKWRK(*), IMCWRK(*), RMCWRK(*),
     1          ISKWRK(*), RSKWRK(*)
      CHARACTER*(*) CCKWRK(*), CSKWRK(*)
      CHARACTER*16 ILINK, IVERSN, VERNUM
C
C     HEADER INFORMATION
      WRITE (IUNIT) IVERSN
C
C     VERSION NUMBER OF THIS SOLUTION FORMAT
      WRITE (IUNIT) VERNUM
C
      ILINK = 'CKLINK          '
      WRITE (IUNIT) ILINK
      CALL CKSAVE (LOUT, IUNIT, ICKWRK, RCKWRK, CCKWRK)
C
      ILINK = 'MCLINK          '
      WRITE (IUNIT) ILINK
      CALL MCSAVE (LOUT, IUNIT, IMCWRK, RMCWRK)
C
      ILINK = 'SKLINK          '
      WRITE (IUNIT) ILINK
      CALL SKSAVE (LOUT, IUNIT, ISKWRK, RSKWRK, CSKWRK)
C
C     end of SUBROUTINE SPSAVE
      RETURN
      END
C
C---------------------------------------------------------------
      SUBROUTINE SPSENS (ISENS, NEQ, LROUT, LRCRVR, LOUT, A, JJ,
     1   KFIRST, KLAST, LENRGY, LTDIF, LVCOR, LVARTP, FTIME, LMULTI,
     2   WT, SDEN, DT, NTEMP, XGIVEN, TGIVEN, EPS, P, TDISK, TINFTY,
     3   AOMEGA, ASWIRL, RHOINF, UINF, LUINF, LCHEM, X, S, ACT, F, SN,
     4   FN, ICKWRK, RCKWRK, IMCWRK, RMCWRK, ISKWRK, RSKWRK, VIS, COND,
     5   D, DKJ, DTCOEF, XAV, YAV, YV, CP, H, XMF, XMFP, IP, DRDAI,
     6   DRDCL, SDOT, GFAC, SFAC, SITDOT, RU, LKMAX, KDEX, JDEX, QDOT,
     7   XSRC, WSRC, TWAL, EMIS, POWR, TRAD, ERAD, VFAC, LRADB, A6,
     8   FINJ, XINJ, WINJ, TINJ, YINJ, HINJ, BUFFER, LCNDCT, CNDFAC,
     9   CNDT, NEMSG, TEMISG, EMISG,IETCH, BULIN, NIICON, SDEN0,
     +   KKPHAS, LSDEN, NSSOLV, LGRA, LCOMP, USAVE, RSAVE, RKFT, RKRT,
     1   LCALSP, SP, DEP, DEPP, FDEP, LFLAME, FLMT, JFIXT, ASPREAD,
     2   LGDOT, GDOT, KERR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /SPCON / MM, KKGAS, NATJ, KKSURF, KKBULK, KKTOT, NPHASE,
     1               NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, NLBULK,
     2               II, IISUR, LENRTW, LENITW, IASIZE, IPSURF,
     3               IPBULK, IPSDEN, IPGAS, NATJF
C     Integer arrays
      DIMENSION IP(*), ICKWRK(*), IMCWRK(*), ISKWRK(*),
     1          KFIRST(NPHASE), KLAST(NPHASE), KDEX(*), JDEX(*),
     2          NIICON(*), IETCH(*)
      DIMENSION WT(*), EPS(*), XGIVEN(*), TGIVEN(*), X(*), SN(*),
     1          S(*), F(*), FN(*), YV(KKGAS,*), VIS(*), COND(*),
     2          D(KKGAS,*), DTCOEF(KKGAS,*), RCKWRK(*), RMCWRK(*),
     4          RSKWRK(*), A(IASIZE), XAV(*), YAV(*),
     5          DKJ(KKGAS,KKGAS,*), CP(*), H(*), XMF(*), XMFP(*),
     6          SDEN(*), SDOT(*), SITDOT(*), DRDAI(*), DRDCL(*),
     8          ACT(*), A6(*), YINJ(*), HINJ(*), BUFFER(*),
     9          TEMISG(*), EMISG(*), BULIN(*), SDEN0(*),
     *          KKPHAS(*), USAVE(NATJ,*), RSAVE(KKGAS,*),
     1          RKFT(II,*), RKRT(II,*), SP(*), DEP(*),
     2          DEPP(*), FDEP(*), GDOT(*)
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
C
      LOGICAL LENRGY, LTDIF, LVCOR, LVARTP, FTIME, LMULTI, LUINF,
     1        LKMAX, LCHEM, ERROR, FIRST, FLAG, LRADB,
     2        LCNDCT, LSDEN, LCOMP, LCALSP, LFLAME, LSENG, LGDOT
      LOGICAL KERR
C
      CHARACTER*16 ISENS
C
      KERR = .FALSE.
C     compute the unit roundoff and the relative and absolute
C     perturbations for the sensitivity calculation
      U = 1.0
   10 CONTINUE
      U = U*0.5
      COMP = 1.0 + U
      IF (COMP .NE. 1.0) GO TO 10
      ABSOL = SQRT(2.0*U)
      RELAT = SQRT(2.0*U)
C
C     need the gas constant for the heat-of-formation
C     sensitivity normalization
      CALL SKRP   (ISKWRK, RSKWRK, RU, RUC, PATM)
C
C
      WRITE (LROUT)  ISENS
      IF (ISENS .EQ. 'HSENSITIVITY') THEN
C
C        doing KKTOT heat-of-formation sensitivities
         ITOT = KKTOT
      ELSEIF (ISENS .EQ. 'SENSITIVITY') THEN
C
C        doing II gas-phase reaction rate constants first,
C        then  IISUR surface reactions
         LSENG = .TRUE.
         ITOT = II
      ENDIF
C
C     FORM THE JACOBIAN, AND AT THE SAME TIME EVALUATE THE
C     UNPERTURBED FUNCTION FN
C
C     put a copy of the solution into BUFFER work array
      CALL TWCOPY (NEQ, S, BUFFER)
      ERROR = .FALSE.
C
C     flag indicates first time through loop which forms the
C     Jacobian
      FIRST = .TRUE.
C
C     specify that SPFUN will not add the transient term
      FTIME = .FALSE.
C
C     first time through the Jacobian loop, evaluate the
C     transport properties and check for species reordering
      LVARTP = .TRUE.
      LKMAX = .TRUE.
C
C     initialize instruction flag used by TWPREP
      FLAG = .FALSE.
C
C     first time through the Jacobian loop, set the IRATE flag
C     for SPFUN to 2, to store the forward and reverse rate constants
C     so that we can use IRATE=3 on subsequent calls
      IRATE = 2
C
 1500 CONTINUE
C
C     this loop is forming the Jacobian
      CALL TWPREP
     +  (ERROR, LOUT,
     +   A, IASIZE, BUFFER, NATJ, CONDIT, LGRA, 0, IP, JJ, FLAG)
C
      IF (ERROR) THEN
         WRITE (LOUT, '(/1X, A /)')
     +      'FATAL ERROR IN SPSENS, JACOBIAN EVALUATION FAILS'
         KERR = .TRUE.
         RETURN
      ENDIF
C
      IF (FLAG) THEN
C
C        evaluate the function for forming the Jacobian; reverse
C        communication from TWPREP is setting the variable FLAG
C        telling us to do this; upon input to SPFUN, BUFFER contains
C        the (perhaps perturbed) solution vector at which TWPREP
C        wants the residual evaluated.
C
         CALL SPFUN (LOUT, JJ, KFIRST, KLAST, LENRGY, LTDIF, LVCOR,
     1               LVARTP, FTIME, LMULTI, WT, BUFFER(IPSDEN), DT,
     2               NTEMP, XGIVEN, TGIVEN, EPS, P, TDISK, TINFTY,
     3               AOMEGA, ASWIRL, RHOINF, UINF, LUINF, LCHEM, X,
     4               BUFFER(IPSURF), SN(IPGAS), ACT, BUFFER(IPGAS),
     5               F(IPGAS), SN(IPSURF), F(IPSURF), ICKWRK, RCKWRK,
     6               IMCWRK, RMCWRK, ISKWRK, RSKWRK, SDOT, GFAC, SFAC,
     7               SITDOT, VIS, COND, D, DKJ, DTCOEF, XAV, YAV, YV,
     8               CP, H, XMF, XMFP, LKMAX, KDEX, JDEX, QDOT,
     9               XSRC, WSRC, TWAL, EMIS, POWR, TRAD, ERAD, VFAC,
     *               LRADB, FINJ, XINJ, WINJ, TINJ, YINJ, HINJ,
     1               LCNDCT, CNDFAC, CNDT, NEMSG, TEMISG, EMISG,IETCH,
     2               BUFFER(IPBULK), BULIN, SN(IPBULK), SN(IPSDEN),
     3               F(IPBULK), F(IPSDEN), NIICON, SDEN0, KKPHAS,
     4               LSDEN, SDEN, LCOMP, USAVE, RSAVE, RKFT, RKRT,
     5               IRATE, LCALSP, SP(IPSURF), SP(IPBULK),
     6               SP(IPSDEN), SP(IPGAS),RU,BUFFER(IPBULK+KKBULK),
     7               SP(IPBULK+KKBULK), F(IPBULK+KKBULK), LFLAME,
     8               FLMT, JFIXT, ASPREAD, LGDOT, GDOT)
C
C        the residual evaluated for this solution vector is put into
C        BUFFER to pass back to TWPREP
         CALL TWCOPY (NEQ, F, BUFFER)
C
         IF (FIRST) THEN
C
C           save the first (unperturbed) function value in FN
            CALL TWCOPY (NEQ, F, FN)
C
C           after the first time through, don't need to evaluate
C           transport properties
C*****OPTIONS - LVARTP always positive
            LVARTP = .TRUE.
C*****END OPTIONS - LVARTP always positive
C*****OPTIONS - LVARTP positive only on jacobian reevaluations
C            LVARTP = .FALSE.
C*****END OPTIONS - LVARTP positive only on jacobian reevaluations
C
C           after the first time through, don't need species reorder
            LKMAX = .TRUE.
C
C           after the first time through, switch to the most-efficient
C           computation of reaction rates
            IRATE = 3
C
C           switch the flag indicating that this is not the first time
C           through the Jacobian loop
            FIRST = .FALSE.
         ENDIF
C
         GO TO 1500
      ENDIF
C
C     the Jacobian has been formed;
C     compute the raw sensitivity coefficicients with respect to the
C     enthalpy, D(mass fraction)/D(enthalpy)
   75 CONTINUE
C
C     THIS IS THE TOP OF A LOOP FOR EVALUATING SENSITIVITY COEFFICIENTS.
C     IT IS EXECUTED TWICE FOR REACTION-RATE SENSITIVITIES
C     (ONCE FOR GAS REACTIONS AND ONCE FOR SURFACE REACTIONS).
C     THE SAME LOOP ALSO DOES THE HEAT-OF-FORMATION SENSITIVITIES
C     DEPENDING ON THE ISENS FLAG.
C
C     evaluate reaction rates in SPFUN in the normal way
      IRATE = 1
C
C     no species reorder
      LKMAX = .FALSE.
C
      DO 1000 IND = 1, ITOT
C
         IF (ISENS .EQ. 'HSENSITIVITY') THEN
C
C           we are doing heat-of-formation (HOF) sensitivities;
C           get the 6th thermodynamic polynomial coefficient for
C           this species (in every temperature regime).
C           The HOF is directly proportional to this coefficient.
            CALL SKRHEX (IND, ISKWRK, RSKWRK, A6)
C
C           find an amount to perturb this coefficient by
            DP    = RELAT*A6(1) + ABSOL
C
C           perturb the H.O.F. coefficient for each temperature range
            DO 50 L = 1, NTR
               A6(L) = A6(L) + DP
   50       CONTINUE
C
C           put perturbed A6 coefficients into surface CHEMKIN
            CALL SKRHEX (-IND, ISKWRK, RSKWRK, A6)
C
C           if gas-phase, put perturbed A6 coefficient into CHEMKIN
            IF (IND .LE. KKGAS) CALL CKRHEX (-IND, RCKWRK, A6)
         ELSEIF (ISENS .EQ. 'SENSITIVITY') THEN
            IF (LSENG) THEN
C
C              gas-phase reaction sensitivities; get pre-exponential
               CALL CKRDEX (IND, RCKWRK, SAVEP)
C
C              find an amount to perturb the pre-exponential
               DP = RELAT*SAVEP + ABSOL
C
C              put perturbed pre-exponential back into CHEMKIN workspace
               CALL CKRDEX (-IND, RCKWRK, SAVEP+DP)
            ELSE
C
C              surface reaction sensitivities; get pre-exponential
               CALL SKRDEX (IND, ISKWRK, RSKWRK, SAVEP)
C
C              find an amount to perturb the pre-exponential
               DP = RELAT*SAVEP + ABSOL
C
C              put perturbed pre-exponential back into SURFACE workspace
               CALL SKRDEX (-IND, ISKWRK, RSKWRK, SAVEP+DP)
            ENDIF
         ENDIF
C
C        get the residuals using perturbed parameter values
         CALL SPFUN (LOUT, JJ, KFIRST, KLAST, LENRGY, LTDIF, LVCOR,
     1               LVARTP, FTIME, LMULTI, WT, S(IPSDEN), DT,
     2               NTEMP, XGIVEN, TGIVEN, EPS, P, TDISK, TINFTY,
     3               AOMEGA, ASWIRL, RHOINF, UINF, LUINF, LCHEM, X,
     4               S(IPSURF), SN(IPGAS), ACT, S(IPGAS), F(IPGAS),
     5               SN(IPSURF), F(IPSURF), ICKWRK, RCKWRK, IMCWRK,
     6               RMCWRK, ISKWRK, RSKWRK, SDOT, GFAC, SFAC,
     7               SITDOT, VIS, COND, D, DKJ, DTCOEF, XAV, YAV, YV,
     8               CP, H, XMF, XMFP, LKMAX, KDEX, JDEX,
     9               QDOT, XSRC, WSRC, TWAL, EMIS, POWR, TRAD, ERAD,
     *               VFAC, LRADB, FINJ, XINJ, WINJ, TINJ, YINJ, HINJ,
     1               LCNDCT, CNDFAC, CNDT, NEMSG, TEMISG, EMISG,IETCH,
     2               S(IPBULK), BULIN, SN(IPBULK), SN(IPSDEN),
     3               F(IPBULK), F(IPSDEN), NIICON, SDEN0, KKPHAS,
     4               LSDEN, SDEN, LCOMP, USAVE, RSAVE, RKFT, RKRT,
     5               IRATE, LCALSP, SP(IPSURF), SP(IPBULK),
     6               SP(IPSDEN), SP(IPGAS),RU, S(IPBULK+KKBULK),
     7               SP(IPBULK+KKBULK), F(IPBULK+KKBULK), LFLAME,
     8               FLMT, JFIXT, ASPREAD, LGDOT, GDOT)
C
         IF (ISENS .EQ. 'HSENSITIVITY') THEN
C
C           return H.O.F. coefficient to original value
            DO 60 L = 1, NTR
               A6(L) = A6(L) - DP
   60       CONTINUE
            CALL SKRHEX (-IND, ISKWRK, RSKWRK, A6)
            IF (IND .LE. KKGAS) CALL CKRHEX (-IND, RCKWRK, A6)
C
C           normalize H.O.F. sensitivity by 1 Kcal/mole
            SAVEP = 1.E3 / RUC
C
         ELSEIF (ISENS .EQ. 'SENSITIVITY') THEN
            IF (LSENG) THEN
C              restore gas-phase reaction pre-exponential
C
               CALL CKRDEX (-IND, RCKWRK, SAVEP)
            ELSE
C              restore surface reaction pre-exponential
C
               CALL SKRDEX (-IND, ISKWRK, RSKWRK, SAVEP)
            ENDIF
         ENDIF
C
C        calculate the partial of F (the residual) with respect to
C        the parameter, store it's negative in SN
         DO 100 N = 1, NEQ
            SN(N) =  - (F(N)-FN(N)) / DP
  100    CONTINUE
C
C        SOLVE THE MATRIX EQUATION
C        - D F / D ALPHA =  SUM_N ( D F / D PHI_N * D PHI_N/ D ALPHA)
C                SN      =   JACOBIAN * SENSITIVITY_COEFF
C
         CALL TWSOLV (ERROR, LOUT, A, IASIZE, SN, NATJ, LGRA, 0, IP, JJ)
C
C        upon exit from TWSOLV, SN contains the sensitiity coefficient
         IF (ERROR) THEN
            WRITE (LOUT, '(/1X, A /)')
     +         'FATAL ERROR, COULD NOT SOLVE LINEAR EQUATIONS.'
            KERR = .TRUE.
            RETURN
         ENDIF
C
C        normalize the sensitivity coefficients SN,  containing the raw
C        sensitivity matrix DS(n)/DA(ind)
         DO 120 J = 1, JJ
            JSKIP = IPGAS-1 + (J-1)*NATJ
C
C           temperature
            NTJ = JSKIP + NT
            SN(NTJ) = SN(NTJ) * SAVEP / S(NTJ)
C
C           axial velocity
            NHJ = JSKIP + NU
            IF (S(NHJ) .NE. 0.0) SN(NHJ) = SN(NHJ) * SAVEP / S(NHJ)
C
C           radial velocity
            NHV = JSKIP + NV
            IF (S(NHV) .NE. 0.0) SN(NHV)=SN(NHV)*SAVEP/ ABS(S(NHV))
C
C           circumferential velocity
            NHW = JSKIP + NW
            IF (S(NHW) .NE. 0.0) SN(NHW)=SN(NHW)*SAVEP/ ABS(S(NHW))
  120    CONTINUE
C
C        normalize info for surface site fractions
         IF (KKSURF .GT. 0)
     1      CALL SPNORM (SN(IPSURF), S(IPSURF), SAVEP, KKSURF)
C
C        normalize info for bulk species activities
         IF (KKBULK .GT. 0)
     1      CALL SPNORM (SN(IPBULK), S(IPBULK), SAVEP, KKBULK)
C
C        normalize info for surface site densities
         IF (LSDEN .AND. NNSURF.GT.0)
     1      CALL SPNORM (SN(IPSDEN), S(IPSDEN), SAVEP, NNSURF)
C
C        normalize and convert mass-fraction sensitivities to
C        mole fraction sensitivities
         DO 400 J = 1, JJ-1
            SUM = 0.0E0
            JSKIP = IPGAS-1 + (J-1)*NATJ
            NKJ = JSKIP + NYS
            DO 360 M = 1, KKGAS
               SUM = SUM + SN(NKJ+M) / WT(M)
  360       CONTINUE
            NYJ = JSKIP + NY1
            CALL CKMMWY (S(NYJ), ICKWRK, RCKWRK, WTM)
C
            WTMSUM = WTM * SUM
            DO 380 M = 1, KKGAS
               IF (S(NKJ+M) .EQ. 0.0) THEN
                 SN(NKJ+M) = SAVEP
               ELSE
                 SN(NKJ+M) = SAVEP * (SN(NKJ+M) / S(NKJ+M) - WTMSUM)
               ENDIF
  380       CONTINUE
  400    CONTINUE
C
C        normalize production rate sensitivities
         IF (KKBULK .GT. 0)
     1      CALL SPNORM(SN(IPBULK+KKBULK),S(IPBULK+KKBULK),SAVEP,KKBULK)
C
C        store the normalized sensitivity coefficients
         WRITE (LROUT)  IND, (SN(N),N=1,NEQ)
C
 1000 CONTINUE
C
      IF (ISENS .EQ. 'SENSITIVITY') THEN
C        two passes may be needed; once for gas-reactions,
C        once for surface reactions
C
         IF (LSENG .AND. IISUR.GT.0) THEN
C           LSENG is TRUE so have made only one pass through,
C           IISUR is > 0  so there are surface reactions;
C           need to go back for surface reaction rate constant
C           sensitivities
            ITOT = IISUR
C
C           set LSENG to turn off gas-phase rate constant sensitivity
            LSENG = .FALSE.
C
C           go to start of sensitivity loop
            GO TO 75
         ENDIF
      ENDIF
C
C     end of SUBROUTINE SPSENS
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE SPSETB (NATJ, NMAX, KKGAS, KKSURF, KKBULK, NNSURF,
     1                   SFLR, ABOVEG,
     1                   BELOWG, ABOVES, BELOWS, IPSURF, IPBULK,
     2                   IPSDEN, IPGAS, LSDEN, LFLAME)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      LOGICAL LSDEN, LFLAME
C
      DIMENSION ABOVEG(*), BELOWG(*), ABOVES(*), BELOWS(*)
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
C
      IF (KKSURF .GT. 0) THEN
         DO 50 K = 1, KKSURF
            BELOWS(IPSURF+K-1) = SFLR
            ABOVES(IPSURF+K-1) = 1.5
   50    CONTINUE
      ENDIF
C
      IF (KKBULK .GT. 0) THEN
         DO 60 K = 1, KKBULK
            BELOWS(IPBULK+K-1) = SFLR
            ABOVES(IPBULK+K-1) = 1.5
            BELOWS(IPBULK+KKBULK+K-1) = -1.0
            ABOVES(IPBULK+KKBULK+K-1) =  1.0
   60    CONTINUE
      ENDIF
C
      IF (LSDEN .AND. NNSURF.GT.0) THEN
         DO 70 N = 1, NNSURF
            BELOWS(IPSDEN+N-1) = SFLR
            ABOVES(IPSDEN+N-1) = 1.5
   70    CONTINUE
      ENDIF
C
C     axial velocity
      BELOWG(NU) = -1.E30
      ABOVEG(NU) =  1.E30
C
C     radial velocity
      BELOWG(NV) = -1.E30
      ABOVEG(NV) = 1.E30
C
C     circumferential velocity
      BELOWG(NW) = -10.0
      ABOVEG(NW) = 1000.0
C
C     temperature
      BELOWG(NT) = 200.0
      ABOVEG(NT) = 5000.0
C
C     lambda
      BELOWG(NL) = -1.E30
      ABOVEG(NL) =  1.E30
C
C     species
      DO 100 K = 1, KKGAS
         N = NYS + K
         BELOWG(N) = SFLR
         ABOVEG(N) = 10.
  100 CONTINUE
C
      IF (LFLAME) THEN
C        inlet velocity (flame restart)
C
         BELOWG(NUL) = -1.E30
         ABOVEG(NUL) = 1.E30
      ENDIF
C
C     end of SUBROUTINE SPSETB
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE SPSTRT (KKGAS, NMAX, NATJ, JJ, NSHT, LOUT, LUMESH,
     1                   NREAC, NINTM, NPROD, REAC, XINTM, PROD,
     2                   TINFTY, TDISK, AOMEGA, ASWIRL, P, KR, KI, KP,
     3                   NTEMP, XGIVEN, TGIVEN, XSTR, XCEN, XEND, WMIX,
     4                   ICKWRK, RCKWRK, Y, SI, EPS, X, SN, S, IMCWRK,
     5                   RMCWRK, RHOINF, ASPREAD, JDEX, KDEX, KFIRST,
     6                   KLAST, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     7                   NLBULK)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
C     Integer arrays
      DIMENSION KR(NREAC), KI(NINTM), KP(NPROD), ICKWRK(*),
     1          IMCWRK(*), JDEX(*), KDEX(*), KFIRST(*), KLAST(*)
C     Real arrays
      DIMENSION REAC(*), XINTM(*), PROD(*), RCKWRK(*), RMCWRK(*),
     2          Y(*), SI(*), EPS(*), XGIVEN(*), TGIVEN(*), X(*),
     3          S(NATJ,*), SN(NSHT,*)
C
      LOGICAL LUMESH
C
C     initialize mass flux fractions and mass fractions to zero
      DO 100 K = 1, KKGAS
         EPS(K) = 0.0
         N = NYS + K
         DO 100 J = 1, NMAX
            S(N, J) = 0.0
  100 CONTINUE
C
      DO 110 J = 1, NMAX
         JDEX(J) = KKGAS
  110 CONTINUE
C
      IF (NNSURF.GT.0) THEN
         DO 1191 N = NFSURF, NLSURF
           KDEX(N)=KLAST(N) - KFIRST(NFSURF) + 1
 1191    CONTINUE
      ENDIF
C
      IF (NNBULK.GT.0) THEN
         DO 1192 N = NFBULK, NLBULK
            KDEX(N)=KLAST(N) - KFIRST(NFBULK) + 1
 1192    CONTINUE
      ENDIF
C
C     set intermediate gaussians
      DO 140 N = 1,NINTM
         GBAS = 0.0
         GM   = XINTM(KI(N))
         GMIX = 0.15*GM + GBAS
         W    = -LOG((GMIX-GBAS)/GM)/(WMIX/2.)**2
C
         DO 120 J = 1,JJ
            S(NYS+KI(N), J) = GM*EXP(-W*(X(J)-XCEN)**2) + GBAS
  120    CONTINUE
  140 CONTINUE
C
C     sum the intermediates at each J
      DO 180 J = 1, JJ
         SI(J) = 0.0
         DO 175 N = 1, NINTM
            SI(J) = SI(J) + S(NYS+KI(N), J)
  175    CONTINUE
  180 CONTINUE
C
C     set starting species profiles
      DO 260 J = 1, JJ
         CALL SPLNWM (WMIX, XCEN, X(J), XPD, XRE)
         FAC = 1.0-SI(J)
         DO 200 N = 1, NREAC
            S(NYS+KR(N), J) = (XPD*PROD(KR(N)) + XRE*REAC(KR(N))) * FAC
  200    CONTINUE
C
         DO 240 N = 1, NPROD
            DO 220 L = 1, NREAC
               IF (KP(N) .EQ. KR(L)) GO TO 240
  220       CONTINUE
            S(NYS+KP(N), J) = (XPD*PROD(KP(N)) + XRE*REAC(KP(N))) * FAC
  240    CONTINUE
  260 CONTINUE
C
C     set mass flux fraction boundary conditions
      CALL CKXTY (REAC, ICKWRK, RCKWRK, EPS)
C
C     convert starting estimates to mass fraction, if needed
      DO 300 J = 1, JJ
         CALL CKXTY (S(NY1, J), ICKWRK, RCKWRK, Y)
         DO 280 K = 1,KKGAS
            S(NYS+K, J) = Y(K)
  280    CONTINUE
  300 CONTINUE
C
C     compute the reference properties at the inlet JJ
      CALL CKRHOY (P, TINFTY, EPS, ICKWRK, RCKWRK, RHOINF)
C
C     set the temperature profile
      IF (TDISK - TINFTY .GT. 1.0) THEN
         DO 320 J = 1, JJ
            S(NT,J) = SPTEMP(NTEMP, X(J), XGIVEN, TGIVEN)
  320    CONTINUE
      ELSE
         DO 340 J = 1, JJ
            S(NT,J) = TDISK
  340    CONTINUE
      ENDIF
C
C     set the temperature, F, G, L, and H profiles
      DO 360 J = 1, JJ
         S(NT,J) =  SN(NT,J)
         S(NV,J) =  SN(NV,J)
         S(NW,J) =  SN(NW,J)
         S(NU,J) =  SN(NU,J)
         S(NL,J) =  SN(NL,J)
  360 CONTINUE
C
C     end of SUBROUTINE SPSTRT
      RETURN
      END
C
      SUBROUTINE SPPOS (S, IDIM, SPOS)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER(I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER(I-N)
C*****END precision > single
C
      DIMENSION S(IDIM)
      DO 10 I = 1, IDIM
         S(I) = MAX (S(I), SPOS)
   10 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
C
C*****precision > double
      DOUBLE PRECISION FUNCTION SPTEMP (NPTS, X, XX, TT)
      DOUBLE PRECISION   XX(NPTS), TT(NPTS), X, S
C*****END precision > double
C*****precision > single
C      REAL FUNCTION SPTEMP (NPTS, X, XX, TT)
C      REAL              XX(NPTS), TT(NPTS), X, S
C*****END precision > single
C
C  START PROLOGUE
C
C     THIS FUNCTION USES BISECTION TO LINEARLY INTERPOLATE
C     AN ARRAY OF XX,TT PAIRS.  GIVEN AN XX,TT PAIR THIS ROUTINE
C     RETURNS THE INTERPOLATED VALUE OF THE T AT THE POINT X.
C
C      It uses linear interpolation from the boundary (XX,TT)
C      pair to interpolate outside the XX(1) to XX(NPTS) region.
C
C INPUT-
C   NPTS   - NUMBER OF XX,TT PAIRS.
C   X      - LOCATION AT WHICH INTERPOLATED T IS DESIRED.
C   XX     - ARRAY OF X POINTS AT WHICH TT ARE GIVEN.
C   TT     - ARRAY OF T VALUES AT THE XX LOCATIONS.
C
C OUTPUT-
C   SPTEMP - INTERPOLATED T AT POINT X
C
C  END PROLOGUE
C
      INTEGER            NPTS, N, NLO, NHI
C
C Check arguments
C
      IF (NPTS .LE. 1) THEN
         IF (NPTS .EQ. 1) THEN
           SPTEMP = TT(1)
         ELSE
           SPTEMP = 0.0
         ENDIF
         RETURN
      ENDIF
C
      S = 0.0
C
C Check for x outside (1,npts)
C
      IF (X .LE. XX(2)) THEN
         N = 2
         S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
      ELSEIF (X .GE. XX(NPTS-1)) THEN
         N = NPTS-1
         S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
      ELSE
C
C Do a binary search to find the correct interval
C
         NLO = 1
         NHI = NPTS
50       CONTINUE
         N = (NLO+NHI)/2
         IF (X .LT. XX(N)) THEN
            IF (X .LT. XX(N-1)) THEN
               NHI = N
               GO TO 50
            ELSEIF (X .EQ. XX(N-1)) THEN
               SPTEMP = TT(N-1)
               RETURN
            ELSE
               S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
            ENDIF
         ELSEIF (X .GT. XX(N)) THEN
            IF (X .GT. XX(N+1)) THEN
               NLO = N
               GO TO 50
            ELSEIF (X .EQ. XX(N+1)) THEN
               SPTEMP = TT(N+1)
               RETURN
            ELSE
               S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
            ENDIF
         ENDIF
      ENDIF
C
      SPTEMP = TT(N) + S * (X - XX(N))
C
C     end of FUNCTION SPTEMP
      RETURN
      END
C
C--------------------------------------------------------------------
C
      SUBROUTINE SPTRNP (KKGAS, JJ, NATJ, LENRGY, LMULTI, LTDIF,
     $                   P, X, S, YAV, ICKWRK, RCKWRK, IMCWRK, RMCWRK,
     $                   XAV, COND, VIS, D, DT, DKJ)
C
C  START PROLOGUE
C  END PROLOGUE
C
C     dummy parameters
      INTEGER     KKGAS, JJ, NATJ, ICKWRK(*), IMCWRK(*)
      LOGICAL     LENRGY, LTDIF, LMULTI
C*****precision > double
      DOUBLE PRECISION
C*****END precision > double
C*****precision > single
C      REAL
C*****END precision > single
     $            P, X(KKGAS), S(NATJ,JJ), YAV(KKGAS), XAV(KKGAS),
     $            RCKWRK(*), RMCWRK(*), COND(JJ), VIS(JJ),
     $            DKJ(KKGAS,KKGAS,JJ), D(KKGAS,JJ), DT(KKGAS,JJ)
C
      INTEGER       NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
C
C     local variables
      INTEGER     J, K, JP1, N
C*****precision > double
      DOUBLE PRECISION         TAV
C*****END precision > double
C*****precision > single
C      REAL                    TAV
C*****END precision > single
C
      IF (LMULTI) THEN
C        Multicomponent formulas
C
         DO 200 J = 1, JJ-1
            JP1 = J + 1
            TAV = 0.5 * (S(NT,J) + S(NT,JP1))
C
            DO 100 K = 1, KKGAS
               N = NYS + K
               YAV(K) = 0.5 * (S(N,J) + S(N,JP1))
100         CONTINUE
            CALL CKYTX  (YAV, ICKWRK, RCKWRK, XAV)
C
C           cetermine the multicomponent diffusion coefficients at J
            CALL MCMDIF (P, TAV, XAV, KKGAS, IMCWRK, RMCWRK, DKJ(1,1,J))
C
C           determine the mixture conductivity and thermal
C           diffusion coefficients at J
            IF (LENRGY .OR. LTDIF)
     $         CALL MCMCDT (P, TAV, XAV, IMCWRK, RMCWRK, ICKWRK,
     $                      RCKWRK, DT(1,J), COND(J))
C
C           determine the mixture viscosity at J
            CALL MCAVIS (TAV, XAV, RMCWRK, VIS(J))
200      CONTINUE
         RETURN
      ENDIF
C
C     Mixture-averaged formulas
      DO 400 J = 1, JJ-1
         JP1 = J + 1
C
         TAV = 0.5 * (S(NT,J) + S(NT,JP1))
         DO 300 K = 1, KKGAS
            N = NYS + K
            YAV(K) = 0.5 * (S(N,J) + S(N,JP1))
300      CONTINUE
         CALL CKYTX  (YAV, ICKWRK, RCKWRK, XAV)
C
C        determine the mixture averaged diffusion coefficients at J
         CALL MCADIF (P, TAV, XAV, RMCWRK, D(1,J))
C
         IF (LTDIF) THEN
C           determine the thermal conductivity and thermal diffusion
C           coefficients at J
            CALL MCMCDT (P, TAV, XAV, IMCWRK, RMCWRK,
     $                   ICKWRK, RCKWRK, DT(1,J), COND(J))
C
        ELSEIF (LENRGY) THEN
C          determine the mixture conductivitity at J from an
C          approximate expression
C
           CALL MCACON (TAV, XAV, RMCWRK, COND(J))
        ENDIF
C
C        determine the mixture viscosity at J
         CALL MCAVIS (TAV, XAV, RMCWRK, VIS(J))
400   CONTINUE
C
C     end of SUBROUTINE SPTRNP
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SPFPNT (LFLUX, NSOL, ENAM, MM, KNAM, KKGAS, FLUX,
     1                   KNCF, WT, TFLUXA)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION FLUX(*), KNCF(MM,*), WT(*), TFLUXA(*)
      CHARACTER*(*) ENAM(*), KNAM(*)
C
      WRITE (LFLUX, 9000) NSOL, (ENAM(M), M = 1, MM)
 9000 FORMAT (/1X,'SOLUTION: ', I3,
     1        /1X,'Species',T15,'Flux (g/cm**2-sec)',5X,
     2 '* #atoms/wt.="atom-moles"',
     3    /,T32,15(3X,A3,6X))
C
      TFLUX = 0.0
      DO 100 M = 1, MM
         TFLUXA(M) = 0.0
  100 CONTINUE
C
      DO 140 K = 1, KKGAS
         WRITE (LFLUX, 9010) KNAM(K), FLUX(K),
     1         (FLOAT(KNCF(M,K))*FLUX(K)/WT(K), M=1,MM)
         TFLUX = TFLUX + FLUX(K)
         DO 120 M = 1, MM
            TFLUXA(M) = TFLUXA(M) + FLOAT(KNCF(M,K))*FLUX(K)/WT(K)
  120    CONTINUE
  140 CONTINUE
C
      WRITE (LFLUX, 9010) 'Total:', TFLUX, (TFLUXA(M), M = 1, MM)
 9010 FORMAT (1X, A12, T16, 15(1PE12.3))
C
C     end of SUBROUTINE SPFPNT
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE REGRID (SNEW, SOLD, NVAR, JJNEW, JJOLD, XNEW, XOLD,
     1                   WORK, P, R, N , KERR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
        DOUBLE PRECISION         ZERO,         ONE
        PARAMETER               (ZERO = 0.0D0, ONE = 1.0D0)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C        REAL                     ZERO,       ONE
C        PARAMETER               (ZERO = 0.0, ONE = 1.0)
C*****END precision > single
C
      DIMENSION SNEW(NVAR,JJNEW), SOLD(NVAR,JJOLD), XNEW(JJNEW),
     $          XOLD(JJOLD), WORK(JJOLD)
      LOGICAL KERR
      KERR = .FALSE.
C
C     interpolate the solution onto an equidistributing mesh
C
C       compute weight coefficients - Right now, this error
C       weighting function just covers the temperature unknown.
C       However, in the future, it should cover all variables in
C       an appropriately weighted fashion that reflects the error
C       norms in the problem.
C
      R0 = ONE - P
      R1 = P * R / (R + ONE)
      R2 = P - R1
C
C     Find the total magnitude of the change in values
      TV1 = ZERO
      DO 10 I = 2, JJOLD
        TV1 = TV1 + ABS( SOLD(N,I) - SOLD(N,I-1) )
10    CONTINUE
C
C     Find the total magnitude of the change in slopes
      TV2  = ZERO
      DXP1 = XOLD(2)   - XOLD(1)
      SXP1 = SOLD(N,2) - SOLD(N,1)
      DO 20 I = 2, JJOLD-1
         DX      = DXP1
         DXP1    = XOLD(I+1)   - XOLD(I)
         SX      = SXP1
         SXP1    = SOLD(N,I+1) - SOLD(N,I)
         TV2 = TV2 + ABS( SXP1/DXP1 - SX/DX )
20    CONTINUE
C
C     Find the Normalization constants
      B1 = R1 * (XOLD(JJOLD) - XOLD(1)) / ((ONE - P) * TV1)
      B2 = R2 * (XOLD(JJOLD) - XOLD(1)) / ((ONE - P) * TV2)
C
C     Compute the error weighting function
      WORK(1) = ZERO
      DXP1    = XOLD(2)   - XOLD(1)
      SXP1    = SOLD(N,2) - SOLD(N,1)
      DO 50 I = 2, JJOLD - 1
         DX      = DXP1
         DXP1    = XOLD(I+1)   - XOLD(I)
         SX      = SXP1
         SXP1    = SOLD(N,I+1) - SOLD(N,I)
         WORK(I) = WORK(I-1) + DX + B1*ABS( SX )
     $                            + B2*ABS( SXP1/DXP1 - SX/DX )
   50 CONTINUE
      DX = DXP1
      SX = SXP1
      WORK(JJOLD) = WORK(JJOLD-1) + DX + B1*ABS( SX )
C
C     Renormalize the error function to 1
      DX = ONE / WORK(JJOLD)
      DO 60 I = 2, JJOLD
         WORK(I) = WORK(I) * DX
   60 CONTINUE
C
C     Interpolate onto uniform eta grid to find new x
      XNEW(1)     = XOLD(1)
      XNEW(JJNEW) = XOLD(JJOLD)
      ISTART  = 2
      DETA = ONE / (JJNEW-1)
      ETAJ = ZERO
      DO 80 J = 2, JJNEW-1
         ETAJ = ETAJ + DETA
         DO 70 I = ISTART, JJOLD
            IF (ETAJ .LE. WORK(I)) THEN
               DEL     = (ETAJ - WORK(I-1)) / (WORK(I) - WORK(I-1))
               XNEW(J) = XOLD(I-1) +(XOLD(I) - XOLD(I-1)) * DEL
               GO TO 80
            ELSE
               ISTART = I
            ENDIF
70       CONTINUE
         WRITE (6,*) ' REGRID ERROR: *** VALUE OF ETA NOT FOUND ***'
         KERR = .TRUE.
         RETURN
80    CONTINUE
C
C     interpolate solution...
      CALL INTPL8 (SOLD, SNEW, XOLD, XNEW, NVAR, JJOLD, JJNEW)
C
C     end of SUBROUTINE REGRID
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE INTPL8 (F1, F2, X1, X2, MVAR, N1, N2)
C
C  START PROLOGUE
C       interpolate to get f2(x2) from f1(x1)
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION F1(MVAR,*), F2(MVAR,*), X1(*), X2(*)
C
      DO 100 J = 2, N2-1
         XVAL = X2(J)
         DO 50 I = 2, N1
            IF(XVAL.LE.X1(I)) THEN
               DEL = (XVAL-X1(I-1))/(X1(I)-X1(I-1))
               DO 40 L = 1, MVAR
                  F2(L,J)=F1(L,I-1)+(F1(L,I)-F1(L,I-1))*DEL
40             CONTINUE
               GOTO 80
            ENDIF
50       CONTINUE
         WRITE(7,*) ' *** STOP...INTERPOLATION ERROR.'
80       CONTINUE
100   CONTINUE
C
C     endpoints..
      DO 110 L = 1, MVAR
         F2(L,1)  = F1(L,1)
         F2(L,N2) = F1(L,N1)
110   CONTINUE
C
C     end of SUBROUTINE INTPL8
      RETURN
      END
C
C----------------------------------------------------------------------

      SUBROUTINE TWSTRT (NCALL, LOUT, PMAX, CNTRLS, LVALUE, IVALUE,
     1                   RVALUE, ERROR)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INTEGER CNTRLS, QADAPT, QLEVD, QLEVM, QPADD, QSSABS, QSSAGE,
     1        QSSREL, QSTDY, QSTPS0, QSTPS1, QSTPS2, QSTRID,
     2        QTDABS, QTDAGE, QTDEC, QTDREL, QTINC, QTMAX, QTMIN,
     3        QTOL0, QTOL1, QTOL2, PMAX
      COMMON /QPTR/
     +   QADAPT, QLEVD, QLEVM, QPADD, QSSABS, QSSAGE, QSSREL, QSTDY,
     +   QSTPS0, QSTPS1, QSTPS2, QSTRID, QTDAGE, QTDABS, QTDEC, QTDREL,
     +   QTINC, QTMAX, QTMIN, QTOL0, QTOL1, QTOL2
C
      DIMENSION IVALUE(CNTRLS, 3), RVALUE(CNTRLS, 3)
      LOGICAL LVALUE(CNTRLS, 3), ERROR
C
C///  SET DEFAULE VALUES FOR TWOPNT'S CONTROLS, FOR THREE
C///  DIFFERENT PROBLEMS:  (1)
C                          (2)
C                          (3)
C
      QADAPT =  1
      QLEVD  =  2
      QLEVM  =  3
      QPADD  =  4
      QSSABS =  5
      QSSAGE =  6
      QSSREL =  7
      QSTDY  =  8
      QSTPS0 =  9
      QSTPS1 = 10
      QSTPS2 = 11
      QSTRID = 12
      QTDABS = 13
      QTDAGE = 14
      QTDEC  = 15
      QTDREL = 16
      QTINC  = 17
      QTMAX  = 18
      QTMIN  = 19
      QTOL0  = 20
      QTOL1  = 21
      QTOL2  = 22
      IF (NCALL .LT. 0) RETURN
C
      IF (NCALL .EQ. 0) THEN
         DO 10 N = 1, CNTRLS
            LVALUE(N, 1) = .FALSE.
            IVALUE(N, 1) = 0
            RVALUE(N, 1) = 0.0
   10    CONTINUE
C
         LVALUE(QADAPT,1) = .FALSE.
         IVALUE(QLEVD,1) = 1
         IVALUE(QLEVM,1) = 1
         IVALUE(QPADD,1) = PMAX
         RVALUE(QSSABS,1) = 1.0E-9
         IVALUE(QSSAGE,1) = 20
         RVALUE(QSSREL,1) = 1.0E-4
         LVALUE(QSTDY,1) = .TRUE.
         IVALUE(QSTPS1,1) = 100
         IVALUE(QSTPS2,1) = 50
         RVALUE(QSTRID,1) = 1.0E-6
         RVALUE(QTDABS,1) = 1.0E-9
         IVALUE(QTDAGE,1) = 20
         RVALUE(QTDEC,1) = 2.2
         RVALUE(QTDREL,1) = 1.0E-4
         RVALUE(QTINC,1) = 2.0
         RVALUE(QTMAX,1) = 1.0E-4
         RVALUE(QTMIN,1) = 1.0E-10
         RVALUE(QTOL0,1) = 1.0E-8
         RVALUE(QTOL1,1) = 0.1
         RVALUE(QTOL2,1) = 0.5
         DO 0050 N = 1, CNTRLS
            DO 0040 L = 2, 3
               LVALUE (N, L) = LVALUE (N, 1)
               RVALUE (N, L) = RVALUE (N, 1)
               IVALUE (N, L) = IVALUE (N, 1)
0040        CONTINUE
0050     CONTINUE
         RETURN
      ENDIF
C
      CALL TWSETL (ERROR, LOUT, 'ADAPT', LVALUE(QADAPT,NCALL))
      CALL TWSETI (ERROR, LOUT, 'LEVELD', IVALUE(QLEVD,NCALL))
      CALL TWSETI (ERROR, LOUT, 'LEVELM', IVALUE(QLEVM,NCALL))
      CALL TWSETI (ERROR, LOUT, 'PADD',   IVALUE(QPADD,NCALL))
      CALL TWSETR (ERROR, LOUT, 'SSABS',  RVALUE(QSSABS,NCALL))
      CALL TWSETI (ERROR, LOUT, 'SSAGE',  IVALUE(QSSAGE,NCALL))
      CALL TWSETR (ERROR, LOUT, 'SSREL',  RVALUE(QSSREL,NCALL))
      CALL TWSETL (ERROR, LOUT, 'STEADY', LVALUE(QSTDY,NCALL))
      CALL TWSETI (ERROR, LOUT, 'STEPS0', IVALUE(QSTPS0,NCALL))
      CALL TWSETI (ERROR, LOUT, 'STEPS1', IVALUE(QSTPS1,NCALL))
      CALL TWSETI (ERROR, LOUT, 'STEPS2', IVALUE(QSTPS2,NCALL))
      CALL TWSETR (ERROR, LOUT, 'STRID0', RVALUE(QSTRID,NCALL))
      CALL TWSETR (ERROR, LOUT, 'TDABS',  RVALUE(QTDABS,NCALL))
      CALL TWSETI (ERROR, LOUT, 'TDAGE',  IVALUE(QTDAGE,NCALL))
      CALL TWSETR (ERROR, LOUT, 'TDEC',   RVALUE(QTDEC,NCALL))
      CALL TWSETR (ERROR, LOUT, 'TDREL',  RVALUE(QTDREL,NCALL))
      CALL TWSETR (ERROR, LOUT, 'TINC',   RVALUE(QTINC,NCALL))
      CALL TWSETR (ERROR, LOUT, 'TMAX',   RVALUE(QTMAX,NCALL))
      CALL TWSETR (ERROR, LOUT, 'TMIN',   RVALUE(QTMIN,NCALL))
      CALL TWSETR (ERROR, LOUT, 'TOLER0', RVALUE(QTOL0,NCALL))
      CALL TWSETR (ERROR, LOUT, 'TOLER1', RVALUE(QTOL1,NCALL))
      CALL TWSETR (ERROR, LOUT, 'TOLER2', RVALUE(QTOL2,NCALL))
C
C     end of SUBROUTINE TWSTRT
      RETURN
      END
C
      SUBROUTINE REFINE
     +  (ERROR, TEXT,
     +   ACTIVE, BUFFER, COMPS, LEVELD, LEVELM, MARK, NEWX, PADD, PMAX,
     +   POINTS, RATIO, RATIO1, RATIO2, SIGNAL, SUCCES, TOLER0, TOLER1,
     +   TOLER2, U, VARY1, VARY2, WEIGHT, X)
C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     REFINE
C
C     PERFORM AUTOMATIC GRID SELECTION.
C
C///////////////////////////////////////////////////////////////////////
      IMPLICIT COMPLEX (A - Z)
      CHARACTER ID*9, SIGNAL*(*), WORD*80
      EXTERNAL TWCOPY
      INTEGER
     +   ACT, COMPS, COUNT, FORMER, ITEMP, J, K, LEAST, LEVELD, LEVELM,
     +   MORE, MOST, NEW, OLD, PADD, PMAX, POINTS, ROUTE, SIGNIF, TEXT,
     +   TOTAL, VARY1, VARY2, WEIGHT, IBASIS, I
      LOGICAL ACTIVE, ERROR, MARK, MESS, NEWX, SUCCES
C*****precision > double
      DOUBLE PRECISION
C*****END precision > double
C*****precision > single
C      REAL
C*****END precision > single
     +   BUFFER, DIFFER, LEFT, LENGTH, LOWER, MEAN, RANGE, RATIO,
     +   RATIO1, RATIO2, RIGHT, SCALE, TEMP, TEMP1, TEMP2, TOLER0,
     +   TOLER1, TOLER2, U, UPPER, X
      PARAMETER (ID = 'REFINE:  ')
C*****LIST MESSAGES > NO
      PARAMETER (MESS = .FALSE.)
C*****END LIST MESSAGES > NO
C*****LIST MESSAGES > YES
C      PARAMETER (MESS = .TRUE.)
C*****END LIST MESSAGES > YES
      DIMENSION
     +   ACTIVE(COMPS), BUFFER(COMPS * PMAX), MARK(PMAX), RATIO(2),
     +   RATIO1(PMAX), RATIO2(PMAX), U(COMPS, PMAX), VARY1(PMAX),
     +   VARY2(PMAX), WEIGHT(PMAX), X(PMAX)
*
C*****REFINE: New Weighting + Diagnostics
      LOGICAL JSPECIAL
      INTEGER    ISCALE
      PARAMETER (ISCALE = 100)
C*****END REFINE: New Weighting + Diagnostics
*
C///  SAVE LOCAL VARIABLES DURING RETURNS FOR REVERSE COMMUNICATION.
      SAVE
C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////
C///  EVERY-TIME INITIALIZATION.
C     TURN OFF ALL COMPLETION STATUS FLAGS.
      ERROR = .FALSE.
      NEWX = .FALSE.
      SUCCES = .FALSE.
C     Calculate the packing constant
      IBASIS = COMPS + 1
C///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.
      IF (SIGNAL .NE. ' ') GO TO ROUTE
C///  WRITE ALL MESSAGES.
      IF (MESS .AND. 0 .LT. TEXT) THEN
         WRITE (TEXT, 10005) ID
         WRITE (TEXT, 10001) ID
         WRITE (TEXT, 10004) ID
         WRITE (TEXT, 10002) ID
         WRITE (TEXT, 10010) ID
         GO TO 9001
      ENDIF
C///  LEVELM PRINTING.
      IF (0 .LT. LEVELM .AND. 0 .LT. TEXT) WRITE (TEXT, 10001) ID
C///  CHECK THE ARGUMENTS.
      ERROR = .NOT. (1 .LE. COMPS .AND. 2 .LE. POINTS)
      IF (ERROR) GO TO 9001
      ERROR = .NOT. (0 .LE. PADD)
      IF (ERROR) GO TO 9002
      ERROR = .NOT. (POINTS .LE. PMAX)
      IF (ERROR) GO TO 9003
      COUNT = 0
      DO 1010 J = 1, COMPS
         IF (ACTIVE(J)) COUNT = COUNT + 1
1010  CONTINUE
      ERROR = .NOT. (1 .LE. COUNT)
      IF (ERROR) GO TO 9004
      ERROR = .NOT. (0.0 .LE. TOLER0)
      IF (ERROR) GO TO 9005
      ERROR = .NOT. (0.0 .LE. TOLER1 .AND. TOLER1 .LE. 1.0
     +         .AND. 0.0 .LE. TOLER2 .AND. TOLER2 .LE. 1.0)
      IF (ERROR) GO TO 9006
      COUNT = 0
      DO 1020 K = 1, POINTS - 1
         IF (X(K) .LT. X(K + 1)) COUNT = COUNT + 1
1020  CONTINUE
      ERROR = .NOT. (COUNT .EQ. 0 .OR. COUNT .EQ. POINTS - 1)
      IF (ERROR) GO TO 9007
C///////////////////////////////////////////////////////////////////////
C
C     AT EACH INTERVAL, COUNT THE ACTIVE, SIGNIFICANT COMPONENTS THAT
C     VARY TOO GREATLY.
C
C///////////////////////////////////////////////////////////////////////
      ACT = 0
      SIGNIF = 0
      DO 2010 K = 1, POINTS
         MARK(K) = .FALSE.
         RATIO1(K) = 0.0
         RATIO2(K) = 0.0
         VARY1(K) = 0
         VARY2(K) = 0
2010  CONTINUE
C///  TOP OF THE LOOP OVER THE COMPONENTS.
      DO 2060 J = COMPS, 1, -1
*
C*****REFINE: New Weighting + Diagnostics
         JSPECIAL = .FALSE.
C*****END REFINE: New Weighting + Diagnostics
*
C*****REFINE: J=1 is a mid-node variable
         IF (J .EQ. 1) JSPECIAL = .TRUE.
C*****END REFINE: J=1 is a mid-node variable
*
         IF (ACTIVE(J)) THEN
            ACT = ACT + 1
C///  FIND THE RANGE AND MAXIMUM MAGNITUDE OF THE COMPONENT.
            LOWER = U(J, 1)
            UPPER = U(J, 1)
            DO 2020 K = 2, POINTS
               LOWER = MIN (LOWER, U(J, K))
               UPPER = MAX (UPPER, U(J, K))
2020        CONTINUE
            RANGE = UPPER - LOWER
C///  DECIDE WHETHER THE COMPONENT IS SIGNIFICANT.
            TEMP = MAX (ABS (LOWER), ABS (UPPER))
            IF (ABS (RANGE) .GT. MAX (TOLER0, TOLER0 * TEMP)) THEN
               SIGNIF = SIGNIF + 1
C///  AT EACH INTERVAL, SEE WHETHER THE COMPONENT'S CHANGE EXCEEDS SOME
C///  FRACTION OF THE COMPONENT'S GLOBAL CHANGE.
               DO 2030 K = 1, POINTS - 1
                  DIFFER = ABS (U(J, K + 1) - U(J, K))
                  IF (0.0 .LT. RANGE)
     +               RATIO1(K) = MAX (RATIO1(K), DIFFER / RANGE)
                  IF (TOLER1 * RANGE .LT. DIFFER) THEN
*
C*****REFINE: J=1 is a mid-node variable
*
*                       These variables at K lie midway between
*                       X(K) and X(K+1).  Therefore, a decision must
*                       be made as to which interval to halve.
*
                     IF (JSPECIAL .AND. (K .LT. (POINTS-1))) THEN
                         IF (ABS(X(K+2)-X(K+1)) .GT.
     $                       ABS(X(K+1)-X(K))       ) THEN
                           VARY1(K+1) = VARY1(K+1) + 1
                           CALL PCKINT (VARY1(K+1), J, IBASIS)
                         ELSE
                           VARY1(K) = VARY1(K) + 1
                           CALL PCKINT (VARY1(K), J, IBASIS)
                         ENDIF
                     ELSE
                       VARY1(K) = VARY1(K) + 1
                       CALL PCKINT (VARY1(K), J, IBASIS)
                     ENDIF
C*****END REFINE: J=1 is a mid-node variable
*
C*****REFINE: Old
C                     VARY1(K) = VARY1(K) + 1
C                     CALL PCKINT (VARY1(K), J, IBASIS)
C*****END REFINE: Old
*
                  ENDIF
2030           CONTINUE
C///  FIND THE GLOBAL CHANGE OF THE COMPONENT'S ARCTANGENT.
               SCALE = ABS ((X(POINTS) - X(1)) / RANGE)
               TEMP = ATAN (SCALE * (U(J, 2) - U(J, 1)) / (X(2) - X(1)))
               LOWER = TEMP
               UPPER = TEMP
               DO 2040 K = 2, POINTS - 1
                  TEMP = ATAN (SCALE * (U(J, K + 1) - U(J, K))
     +               / (X(K + 1) - X(K)))
                  LOWER = MIN (LOWER, TEMP)
                  UPPER = MAX (UPPER, TEMP)
2040           CONTINUE
               RANGE = UPPER - LOWER
C///  AT EACH INTERIOR POINT, SEE WHETHER THE DERIVATIVE'S CHANGE
C///  EXCEEDS SOME FRACTION OF THE DERIVATIVE'S GLOBAL CHANGE.
*
C*****REFINE: J=1 is a mid-node variable
               IF (JSPECIAL) THEN
                  RIGHT = ATAN (SCALE * (U(J, 2) - U(J, 1)) * 2.0
     +            / (X(3) - X(1)))
               ELSE
                  RIGHT = ATAN (SCALE * (U(J, 2) - U(J, 1))
     +            / (X(2) - X(1)))
               ENDIF
C*****END REFINE: J=1 is a mid-node variable
*
C*****REFINE: Old
C               RIGHT = ATAN (SCALE * (U(J, 2) - U(J, 1))
C     +            / (X(2) - X(1)))
C*****END REFINE: Old
*
C*****REFINE: New Weighting + Diagnostics
               RATIO2(1)      = 0.0
               RATIO2(POINTS) = 0.0
C*****END REFINE: New Weighting + Diagnostics
*
               DO 2050 K = 2, POINTS - 1
                  LEFT = RIGHT
*
C*****REFINE: J=1 is a mid-node variable
                  IF (JSPECIAL) THEN
                    IF (K .LT. (POINTS-1)) THEN
                      RIGHT = ATAN (SCALE * (U(J, K + 1) - U(J, K))
     +                              * 2.0 / (X(K + 2) - X(K)))
                    ELSE
                      RIGHT = ATAN (SCALE * (U(J, K + 1) - U(J, K))
     +                              * 2.0 / (X(K + 1) - X(K)))
                    ENDIF
                  ELSE
                    RIGHT = ATAN (SCALE * (U(J, K + 1) - U(J, K))
     +                            / (X(K + 1) - X(K)))
                  ENDIF
C*****END REFINE: J=1 is a mid-node variable
*
C*****REFINE: Old
C                  RIGHT = ATAN (SCALE * (U(J, K + 1) - U(J, K))
C     +               / (X(K + 1) - X(K)))
C*****END REFINE: Old
*
                  DIFFER = ABS (LEFT - RIGHT)
                  IF (0.0 .LT. RANGE)
     +               RATIO2(K) = MAX (RATIO2(K), DIFFER / RANGE)
*
C*****REFINE: Old
C                  IF (TOLER2 * RANGE .LT. DIFFER)
C     +               VARY2(K) = VARY2(K) + 1
C*****END REFINE: Old
*
C*****REFINE: New Weighting + Diagnostics
                   IF (TOLER2 * RANGE .LT. DIFFER) THEN
                      VARY2(K) = VARY2(K) + 1
                      CALL PCKINT (VARY2(K), J, IBASIS)
                   ENDIF
C*****END REFINE: New Weighting + Diagnostics
*
2050           CONTINUE
C///  BOTTOM OF THE LOOP OVER THE COMPONENTS.
            ENDIF
         ENDIF
2060  CONTINUE
C///  SAVE THE MAXIMUM RATIOS.
      RATIO(1) = 0.0
      DO 2070 K = 1, POINTS - 1
         RATIO(1) = MAX (RATIO(1), RATIO1(K))
2070  CONTINUE
      RATIO(2) = 0.0
      DO 2080 K = 2, POINTS - 1
         RATIO(2) = MAX (RATIO(2), RATIO2(K))
2080  CONTINUE
C///////////////////////////////////////////////////////////////////////
C
C     SELECT THE INTERVALS TO HALVE.
C
C///////////////////////////////////////////////////////////////////////
C///  WEIGHT THE INTERVALS IN WHICH VARIATIONS THAT ARE TOO LARGE OCCUR.
*
C*****REFINE: Old
C      MOST = 0
C*****END REFINE: Old
*
C*****REFINE: New Weighting + Diagnostics
      DO 3009 K = 1, POINTS - 1
        WEIGHT(K) = 0
 3009 CONTINUE
C*****END REFINE: New Weighting + Diagnostics
*
      DO 3010 K = 1, POINTS - 1
*
C*****REFINE: Old
C         WEIGHT(K) = VARY1(K)
C         IF (1 .LT. K) WEIGHT(K) = WEIGHT(K) + VARY2(K)
C         IF (K .LT. POINTS - 1) WEIGHT(K) = WEIGHT(K) + VARY2(K + 1)
C         IF (0 .LT. WEIGHT(K)) MOST = MOST + 1
C*****END REFINE: Old
*
C*****REFINE: New Weighting + Diagnostics
         ITMP = MAX (1, INT (RATIO1(K) / TOLER1))
         ITMPX = INT (  ABS(X(K+1) -X(K)) * POINTS * ISCALE
     $                / ABS(X(POINTS)-X(1)))
         WEIGHT(K) = WEIGHT(K) +
     $               (ITMP*ISCALE + ITMPX) * MOD(VARY1(K),IBASIS)
*
*          Decide where to put the weighting factor for VARY2
*          depending on which of the two intervals is larger.
*
         IF (K .GT. 1) THEN
           IF (ABS(X(K)-X(K-1)) .GT. ABS(X(K+1)-X(K))) THEN
             ITMPX = INT (  ABS(X(K) -X(K-1)) * POINTS * ISCALE
     $                / ABS(X(POINTS)-X(1)))
             ITMP = MAX (1, INT ((RATIO2(K)+RATIO2(K-1)) / TOLER2))
             WEIGHT(K-1) = WEIGHT(K-1) +
     $               (ITMP*ISCALE + ITMPX)* MOD(VARY2(K),IBASIS)
           ELSE
             ITMP = MAX (1, INT ((RATIO2(K)+RATIO2(K+1)) / TOLER2))
             WEIGHT(K) = WEIGHT(K) +
     $               (ITMP*ISCALE + ITMPX)* MOD(VARY2(K),IBASIS)
           ENDIF
         ENDIF
*
C*****END REFINE: New Weighting + Diagnostics
*
3010  CONTINUE
*
C*****REFINE: New Weighting + Diagnostics
      MOST = 0
      DO 3011 K = 1, POINTS-1
        IF (0 .LT. WEIGHT(K)) MOST = MOST + 1
3011  CONTINUE
C*****END REFINE: New Weighting + Diagnostics
*
C///  SORT THE WEIGHTS.
      DO 3030 K = 1, POINTS - 1
         DO 3020 J = K + 1, POINTS - 1
            IF (WEIGHT(J) .GT. WEIGHT(K)) THEN
               ITEMP = WEIGHT(J)
               WEIGHT(J) = WEIGHT(K)
               WEIGHT(K) = ITEMP
            ENDIF
3020     CONTINUE
         IF (WEIGHT(K) .EQ. 0) GO TO 3040
3030  CONTINUE
3040  CONTINUE
C///  FIND THE LEAST WEIGHT OF INTERVALS TO HALVE.
      MORE = MAX (0, MIN (MOST, PADD, PMAX - POINTS))
      IF (0 .LT. MORE) THEN
         LEAST = WEIGHT(MORE)
      ELSE
         LEAST = 1 + WEIGHT(1)
      ENDIF
C///  RECONSTRUCT THE WEIGHTS.
*
C*****REFINE: New Weighting + Diagnostics
      DO 3049 K = 1, POINTS - 1
        WEIGHT(K) = 0
 3049 CONTINUE
C*****END REFINE: New Weighting + Diagnostics
*
      DO 3050 K = 1, POINTS - 1
*
C*****REFINE: Old
C         WEIGHT(K) = VARY1(K)
C         IF (1 .LT. K) WEIGHT(K) = WEIGHT(K) + VARY2(K)
C         IF (K .LT. POINTS - 1) WEIGHT(K) = WEIGHT(K) + VARY2(K + 1)
C*****END REFINE: Old
*
C*****REFINE: New Weighting + Diagnostics
*
         ITMP = MAX (1, INT (RATIO1(K) / TOLER1))
         ITMPX = INT (  ABS(X(K+1) -X(K)) * POINTS * ISCALE
     $                / ABS(X(POINTS)-X(1)))
         WEIGHT(K) =  (ITMP*ISCALE + ITMPX) * MOD(VARY1(K),IBASIS)
*
*          Decide where to put the weighting factor for VARY2
*          depending on which of the two intervals is larger.
*
         IF (K .GT. 1) THEN
           IF (ABS(X(K)-X(K-1)) .GT. ABS(X(K+1)-X(K))) THEN
             ITMPX = INT (  ABS(X(K) -X(K-1)) * POINTS * ISCALE
     $                / ABS(X(POINTS)-X(1)))
             ITMP = MAX (1, INT ((RATIO2(K)+RATIO2(K-1)) / TOLER2))
             WEIGHT(K-1) = WEIGHT(K-1) +
     $               (ITMP*ISCALE + ITMPX)* MOD(VARY2(K),IBASIS)
           ELSE
             ITMP = MAX (1, INT ((RATIO2(K)+RATIO2(K+1)) / TOLER2))
             WEIGHT(K) = WEIGHT(K) +
     $               (ITMP*ISCALE + ITMPX)* MOD(VARY2(K),IBASIS)
           ENDIF
         ENDIF
*
C*****END REFINE: New Weighting + Diagnostics
*
3050  CONTINUE
C///  MARK THE INTERVALS TO HALVE.
      COUNT = 0
      DO 3060 K = 1, POINTS - 1
         IF (COUNT .LT. MORE .AND. LEAST .LE. WEIGHT(K)) THEN
            COUNT = COUNT + 1
            MARK(K) = .TRUE.
         ENDIF
3060  CONTINUE
      MORE = COUNT
C///////////////////////////////////////////////////////////////////////
C
C     HALVE THE INTERVALS, IF ANY.
C
C///////////////////////////////////////////////////////////////////////
C///  FORM THE TOTAL NUMBER OF POINTS IN THE NEW GRID.
      TOTAL = POINTS + MORE
      IF (0 .EQ. MORE) GO TO 4070
C///  TOP OF THE BLOCK TO CREATE THE NEW GRID.  CHECK THE ORDER.
      COUNT = 0
      LENGTH = ABS (X(POINTS) - X(1))
      DO 4010 K = 1, POINTS - 1
         IF (MARK(K)) THEN
            MEAN = 0.5 * (X(K) + X(K + 1))
            IF (.NOT. ((X(K) .LT. MEAN .AND. MEAN .LT. X(K + 1))
     +         .OR. (X(K + 1) .LT. MEAN .AND. MEAN .LT. X(K))))
     +         COUNT = COUNT + 1
         ENDIF
4010  CONTINUE
      ERROR = .NOT. (COUNT .EQ. 0)
      IF (ERROR) GO TO 9008
C///  ADD THE NEW POINTS.  INTERPOLATE X AND THE BOUNDS.
      NEW = TOTAL
      DO 4040 OLD = POINTS, 2, - 1
         X(NEW) = X(OLD)
         DO 4020 J = 1, COMPS
            U(J, NEW) = U(J, OLD)
4020     CONTINUE
         NEW = NEW - 1
         IF (MARK(OLD - 1)) THEN
            X(NEW) = 0.5 * (X(OLD) + X(OLD - 1))
            DO 4030 J = 1, COMPS
               U(J, NEW) = 0.5 * (U(J, OLD) + U(J, OLD - 1))
4030        CONTINUE
            NEW = NEW - 1
         ENDIF
4040  CONTINUE
C///  MARK THE NEW POINTS.
      NEW = TOTAL
      DO 4050 OLD = POINTS, 2, - 1
         MARK(NEW) = .FALSE.
         NEW = NEW - 1
         IF (MARK(OLD - 1)) THEN
            MARK(NEW) = .TRUE.
            NEW = NEW - 1
         ENDIF
4050  CONTINUE
      MARK(NEW) = .FALSE.
C///  UPDATE THE NUMBER OF POINTS.
      FORMER = POINTS
      POINTS = TOTAL
C///  ALLOW THE USER TO UPDATE THE SOLUTION.
      CALL TWCOPY (COMPS * POINTS, U, BUFFER)
      SIGNAL = 'UPDATE'
      ASSIGN 4060 TO ROUTE
      GO TO 99999
4060  CONTINUE
      SIGNAL = ' '
      CALL TWCOPY (COMPS * POINTS, BUFFER, U)
C///  BOTTOM OF THE BLOCK TO CREATE A NEW GRID.
4070  CONTINUE
C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////
C///  PRINT.
      IF (0 .LT. LEVELM .AND. 0 .LT. TEXT) THEN
         TEMP1 = RATIO1(1)
         DO 5010 K = 2, FORMER - 1
            TEMP1 = MAX (TEMP1, RATIO1(K))
5010     CONTINUE
         TEMP2 = RATIO2(2)
         DO 5020 K = 3, FORMER - 1
            TEMP2 = MAX (TEMP2, RATIO2(K))
5020     CONTINUE
         IF (SIGNIF .EQ. 0) THEN
            WRITE (TEXT, 10002) ID
         ELSE
            WRITE (TEXT, 10003) TEMP1, TEMP2, TOLER1, TOLER2
            IF (MOST .EQ. 0) THEN
               WRITE (TEXT, 10004) ID
            ELSEIF (MORE .EQ. 0) THEN
               WRITE (TEXT, 10005) ID
            ELSE
               WRITE (TEXT, 10006)
               OLD = 0
               DO 5030 K = 1, POINTS
                  WORD = ' '
                  IF (.NOT. MARK(K)) THEN
                     OLD = OLD + 1
                     IF (1 .LT. K) THEN
                        IF (VARY1(OLD - 1) .NE. 0) THEN
*
C*****REFINE: Old
C                           WRITE (WORD, '(F4.2, I4)') RATIO1(OLD-1),
C     +                                                VARY1(OLD-1)
C                           I = 8
C*****END REFINE: Old
*
C*****REFINE: New Weighting + Diagnostics
*
                           WRITE (WORD, '(F4.2, I4)') RATIO1(OLD - 1),
     +                        MOD(VARY1(OLD - 1),IBASIS)
                           I = 23
                           WRITE (WORD(I-2:I),'('' ('')')
 5024                      CONTINUE
                           CALL UPKINT(VARY1(OLD-1),J,IBASIS)
                           IF (J .NE. 0) THEN
                             IF (I .GT. 23) THEN
                                WRITE (WORD(I:I+2),'('', '')')
                                I = I + 2
                             ENDIF
                             IF (J .LT. 10) THEN
                               WRITE (WORD(I:I+1), '(I1)') J
                               I = I + 1
                             ELSEIF (J .LT. 100) THEN
                               WRITE (WORD(I:I+2), '(I2)') J
                               I = I + 2
                             ELSE
                               WRITE (WORD(I:I+3), '(I3)') J
                               I = I + 3
                             ENDIF
                             IF (I .LT. 50)   GO TO 5024
                           ENDIF
                           WRITE (WORD(I:I+1),'('')'')')
                           I = I + 1
*
C*****END REFINE: New Weighting + Diagnostics
*
                        ELSE
                           WRITE (WORD, '(F4.2, I4)') RATIO1(OLD - 1)
                        ENDIF
                        IF (MARK(K - 1)) THEN
                           WRITE (TEXT, 10007) K - 1, X(K - 1),
     +                           WORD(1:I)
                        ELSE
                           WRITE (TEXT, 10008) WORD(1:I)
                        ENDIF
                     ENDIF
*
C*****REFINE: Old
C                    IF (1 .LT. K .AND. K .LT. POINTS) THEN
C*****END REFINE: Old
*
                        IF (VARY2(OLD) .NE. 0) THEN
*
C*****REFINE: Old
C                           WRITE (WORD, '(F4.2, I4)') RATIO2(OLD),
C     +                                                VARY2(OLD)
C                           I = 8
C*****END REFINE: Old
*
C*****REFINE: New Weighting + Diagnostics
*
                           WRITE (WORD, '(F4.2, I4)') RATIO2(OLD),
     +                                 MOD(VARY2(OLD),IBASIS)
                           I = 12
                           WRITE (WORD(I-1:I),'(''('')')
 5026                      CONTINUE
                           CALL UPKINT (VARY2(OLD), J, IBASIS)
                           IF (J .NE. 0) THEN
                             IF (I .GT. 12) THEN
                                WRITE (WORD(I:I+2),'('', '')')
                                I = I + 2
                             ENDIF
                             IF (J .LT. 10) THEN
                               WRITE (WORD(I:I+1), '(I1)') J
                               I = I + 1
                             ELSEIF (J .LT. 100) THEN
                               WRITE (WORD(I:I+2), '(I2)') J
                               I = I + 2
                             ELSE
                               WRITE (WORD(I:I+3), '(I3)') J
                               I = I + 3
                             ENDIF
                             IF (I .LT. 50) GO TO 5026
                           ENDIF
                           WRITE (WORD(I:I+1),'('')'')')
                           I = I + 1
*
C*****END REFINE: New Weighting + Diagnostics
*
                        ELSE
                           WRITE (WORD, '(F4.2, I4)') RATIO2(OLD)
                        ENDIF
                        WRITE (TEXT, 10009) K, X(K), WORD(1:I)
*
C*****REFINE: Old
C                     ELSE
C                        WRITE (TEXT, 10009) K, X(K)
C                     ENDIF
C*****END REFINE: Old
*
                  ENDIF
5030           CONTINUE
            ENDIF
         ENDIF
         IF (0 .LT. LEVELD .AND. 0 .LT. MORE) THEN
            WRITE (TEXT, 10010) ID
            CALL TWCOPY (COMPS * POINTS, U, BUFFER)
            SIGNAL = 'SHOW'
            ASSIGN 5040 TO ROUTE
            GO TO 99999
         ENDIF
      ENDIF
5040  CONTINUE
      SIGNAL = ' '
C///  SET THE COMPLETION STATUS FLAGS.
      NEWX = 0 .LT. MORE
      SUCCES = 0 .EQ. MOST
C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////
10001 FORMAT
     +  (/1X, A9, 'SELECT A GRID.')
10002 FORMAT
     +  (/1X, A9, 'SUCCESS.  THE GRID IS ADEQUATE BECAUSE ALL ACTIVE'
     +  /10X, 'COMPONENTS ARE INSIGNIFICANT.')
10003 FORMAT
C              123456789-   1234567   1234567
     + (/15X, '             RATIO 1   RATIO 2'
     +  /15X, '             -------   -------'
     +  /15X, '    ACTUAL', 2F10.3
     +  /15X, '   DESIRED', 2F10.3)
10004 FORMAT
     +  (/1X, A9, 'SUCCESS.  THE GRID IS ADEQUATE.')
10005 FORMAT
     +  (/1X, A9, 'FAILURE.  MORE POINTS ARE NEEDED BUT NONE CAN BE'
     +  /10X, 'ADDED.')
10006 FORMAT
     + (/10X, 'THE NEW GRID (* MARKS NEW POINTS):'
C              123456   123456789-123456   12345678   12345678
     + //10X, '                             LARGEST RATIOS AND',
     +    3X, '(SOME OF THE UNREFINED)'
     +  /10X, ' INDEX         GRID POINT      NUMBER TOO LARGE',
     +    3X, '( SOLUTION COMPONENTS )'
     +  /10X, '------   ----------------   -------------------',
     +    3X, '-----------------------'
     +  /10X, '                             RATIO 1    RATIO 2')
10007 FORMAT
     +  (10X, I6, '*  ', 1P, E16.9, 0P, 3X, A)
10008 FORMAT
     +  (38X, A)
10009 FORMAT
     +  (10X, I6, '   ', 1P, E16.9, 0P, 14X, A)
10010 FORMAT
     +  (/1X, A9, 'THE SOLUTION GUESS FOR THE NEW GRID:')
C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////
      GO TO 99999
9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, COMPS, POINTS
      IF (.NOT. MESS) GO TO 99999
9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, PADD
      IF (.NOT. MESS) GO TO 99999
9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, POINTS, PMAX
      IF (.NOT. MESS) GO TO 99999
9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID
      IF (.NOT. MESS) GO TO 99999
9005  IF (0 .LT. TEXT) WRITE (TEXT, 99005) ID, TOLER0
      IF (.NOT. MESS) GO TO 99999
9006  IF (0 .LT. TEXT) WRITE (TEXT, 99006) ID, TOLER1, TOLER2
      IF (.NOT. MESS) GO TO 99999
9007  IF (0 .LT. TEXT) WRITE (TEXT, 99007) ID
      IF (.NOT. MESS) GO TO 99999
9008  IF (0 .LT. TEXT) WRITE (TEXT, 99008) ID
      IF (.NOT. MESS) GO TO 99999
99001 FORMAT
     +  (/1X, A9, 'ERROR.  THERE MUST BE AT LEAST ONE COMPONENT AND AT'
     +  /10X, 'LEAST TWO POINTS.'
     + //10X, I10, '  COMPS, COMPONENTS'
     +  /10X, I10, '  POINTS')
99002 FORMAT
     +  (/1X, A9, 'ERROR.  THE LIMIT ON POINTS ADDED TO A GRID MUST BE'
     +  /10X, 'ZERO OR POSITIVE.'
     + //10X, I10, '  PADD, LIMIT ON ADDED POINTS')
99003 FORMAT
     +  (/1X, A9, 'ERROR.  POINTS IS OUT OF RANGE.'
     + //10X, I10, '  POINTS'
     +  /10X, I10, '  PMAX, LIMIT ON POINTS')
99004 FORMAT
     +  (/1X, A9, 'ERROR.  THERE ARE NO ACTIVE COMPONENTS.')
99005 FORMAT
     +  (/1X, A9, 'ERROR.  THE BOUNDS ON MAGNITUDE AND RELATIVE CHANGE'
     +  /10X, 'OF MAGNITUDE FOR INSIGNIFICANT COMPONENTS MUST BE'
     +  /10X, 'POSITIVE.'
     + //10X, 1P, E10.2, '  TOLER0, SIGNIFICANCE LEVEL')
99006 FORMAT
     +  (/1X, A9, 'ERROR.  THE BOUNDS ON RELATIVE CHANGES IN MAGNITUDE'
     +  /10X, 'AND ANGLE MUST LIE BETWEEN 0 AND 1.'
     + //10X, 1P, E10.2, '  TOLER1'
     +  /10X, 1P, E10.2, '  TOLER2')
99007 FORMAT
     +  (/1X, A9, 'ERROR.  THE GRID IS NOT ORDERED.')
99008 FORMAT
     +  (/1X, A9, 'ERROR.  SOME INTERVALS IN THE GRID ARE TOO SHORT.'
     +  /10X, 'THE NEW GRID WOULD NOT BE ORDERED.')
C
C     end of SUBROUTINE REFINE
      STOP
C
99999 CONTINUE
      RETURN
      END
*
C*****REFINE: New Weighting + Diagnostics
*
      SUBROUTINE UPKINT (IVALUE, ICOMP, IBASIS)
*
      INTEGER            IVALUE, ICOMP, IBASIS
C
C            2**32 - assumed to be the max int value allowed
C
      INTEGER    IMAX
      PARAMETER (IMAX = 2147000000)
      INTEGER    IDIVISOR, JBASIS, ITEMP, IMAX1
C
C Unpacks a small value from an integer
C   If there is nothing to unpack - it returns ICOMP=0
C   and IVALUE unchanged.
C
      IF (IBASIS .GT. 2) THEN
        ITEMP  = IBASIS
        JBASIS = IBASIS
      ELSE
        ITEMP  = 3
        JBASIS = 3
      ENDIF
      IDIVISOR = 1
      IMAX1 = IMAX / IBASIS
10    CONTINUE
      IF (IVALUE/ITEMP .GT. 0) THEN
         IDIVISOR = ITEMP
         IF (ITEMP .LT. IMAX1) THEN
           ITEMP = ITEMP * JBASIS
           GO TO 10
         ENDIF
      ENDIF
      IF (IDIVISOR .GT. 1) THEN
        ICOMP = IVALUE / IDIVISOR
        IVALUE = IVALUE - ICOMP * IDIVISOR
      ELSE
        ICOMP = 0
      ENDIF
C
C     end of SUBROUTINE UPKINT
      RETURN
      END
      SUBROUTINE PCKINT (IVALUE, ICOMP, IBASIS)
*
      INTEGER            IVALUE, ICOMP, IBASIS
C
C            2**31 - assumed to be the max int value allowed
C
      INTEGER    IMAX
      PARAMETER (IMAX = 2147000000)
      INTEGER    JBASIS, IMAX1, IDIVISOR
C
C Packs an integer with a small value, if there is room
C
      IF (IBASIS .GT. 2) THEN
        JBASIS = IBASIS
      ELSE
        JBASIS = 3
      ENDIF
      IDIVISOR = JBASIS
      IMAX1 = IMAX / JBASIS
10    CONTINUE
      IF ((IVALUE/IDIVISOR) .EQ. 0) THEN
          IVALUE = IVALUE + ICOMP*IDIVISOR
          RETURN
      ENDIF
      IDIVISOR = IDIVISOR * JBASIS
      IF (IDIVISOR .LT. IMAX1) GO TO 10
C
C Integer is maxed out - return without an error message
C
C     end of SUBROUTINE PCKINT
      RETURN
      END
C*****END REFINE: New Weighting + Diagnostics
C
C---------------------------------------------------------------------
C
      SUBROUTINE SPFLJJ (KKGAS, NATJ, JJ, ICKWRK, RCKWRK, P, S,
     1                   X, YAV, WT, SDOT, FINJ, WINJ, XINJ, YINJ,
     2                   YV, FLUX)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (PI=3.141592654D0, ZERO=0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (PI=3.141592654, ZERO=0.0)
C*****END precision > single
C
      COMMON /LOCS/ NU, NV, NW, NT, NL, NYS, NY1, NTR, NUL
      DIMENSION ICKWRK(*), RCKWRK(*), S(NATJ,*),
     1          YV(KKGAS,*), FLUX(*), X(*), WT(*), SDOT(*),
     2          YINJ(*), YAV(*)
C
      JJ1 = JJ - 1
      CALL CKRHOY (P, S(NT,JJ), S(NY1,JJ), ICKWRK, RCKWRK, RHOJJ)
      DELTAX = (X(JJ) - X(JJ1)) / 2.0
      UJJ = S(NU,JJ)
C
      DO 100 K = 1, KKGAS
         N = NYS + K
         YAV(K) = 0.5 * (S(N,JJ) + S(N,JJ1))
  100 CONTINUE
      TAV = 0.5 * (S(NT,JJ) + S(NT,JJ1))
      CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOM)
      CALL CKWYP  (P, S(NT,JJ), S(NY1,JJ), ICKWRK, RCKWRK, SDOT)
C
C     Determine the amount of injectant at the top
      FDOT = ZERO
      IF (FINJ .GT. ZERO) THEN
         XCUT = 2 * WINJ
         IF (ABS (X(JJ)-XINJ) .LE. XCUT) THEN
            FO = 1.283153 * FINJ * SQRT (3.0/PI) / WINJ
            FDOT = FO * EXP (-3.0*(X(JJ)-XINJ)**2 / WINJ**2)
         ENDIF
      ENDIF
C
      JJM1 = JJ - 1
      DO 1100 K = 1, KKGAS
         N = NYS + K
         FLUX(K) = RHOJJ*S(NU,JJ)*S(N,JJ) +
     1             RHOM*YV(K,JJM1) +
     2             DELTAX*WT(K)*SDOT(K) +
     3             DELTAX*FDOT*(YINJ(K)-S(N,JJ))
 1100 CONTINUE
C
C     end of SUBROUTINE SPFLJJ
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SPSOLN (IUNIT,  SOLTYP, JJ,     NEQ,    NATJ,
     1                   KKGAS,  KKSURF, KKBULK, NSSOLV, IMULTI,
     2                   ITDIF,  IVCOR,  ICOMP,  ICHEM,  IGDOT,
     3                   IREORD, ICNDCT, IUINF,  IENRGY, IRADB,
     4                   P,      UINF,   GFAC,   SFAC,   CNDFAC,
     5                   POWR,   CNDT,   XINJ,   WINJ,   FINJ,
     6                   TINJ,   QDOT,   XSRC,   WSRC,   TWAL,
     7                   EMIS,   ERAD,   VFAC,   TRAD,   NEMSG,
     8                   EPS,    GDOT,   YINJ,   TEMISG, EMISG,
     9                   TIME,   LFLAME, IPGAS,  NATJF,  NBEQ,
     +                   X,      S )
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION EPS(*), GDOT(*), YINJ(*), TEMISG(*), EMISG(*),
     1          X(*), S(*)
C
      CHARACTER*16 SOLTYP
C
      LOGICAL LFLAME
C
C     TELL WHAT KIND OF SOLUTION THIS IS
      WRITE (IUNIT) SOLTYP
C
C     PROBLEM SIZE INFORMATION
      WRITE (IUNIT) JJ, NEQ, NATJ, KKGAS, KKSURF, KKBULK, NSSOLV
C
C     OPTIONS USED IN RUNNING THE CALCULATION
      WRITE (IUNIT) IMULTI, ITDIF,  IVCOR, ICOMP,  ICHEM, IGDOT,
     1              IREORD, ICNDCT, IUINF, IENRGY, IRADB
C
C     PARAMETERS IN THE CALCULATION
      WRITE (IUNIT) P,      UINF,   GFAC,  SFAC,   CNDFAC,
     1              POWR,   CNDT,   XINJ,  WINJ,   FINJ,
     2              TINJ,   QDOT,   XSRC,  WSRC,   TWAL,
     3              EMIS,   ERAD,   VFAC,  TRAD,   NEMSG
C
C     ARRAYS SPECIFYING MORE INPUT PARAMETERS
      WRITE (IUNIT) (EPS(K),   K=1,KKGAS), (GDOT(K), K=1,KKGAS),
     1              (YINJ(K),  K=1,KKGAS),
     2              (TEMISG(N),N=1,NEMSG), (EMISG(N),N=1,NEMSG)
C
C     TIME (WILL HAVE MEANING ONLY FOR TRANSIENT SOLUTION)
      WRITE (IUNIT) TIME
C
C     LOCATION OF GRID POINTS
      WRITE (IUNIT) (X(J), J=1,JJ)
C
C     THE SOLUTION ITSELF
      IF (LFLAME) THEN
         WRITE (IUNIT) (S(N), N=1,NBEQ),
     1                 ((S(IPGAS-1+(J-1)*NATJF+N),N=1,NATJ),J=1,JJ)
         NATJ = NATJF
         NEQ = NATJ*JJ + NBEQ
      ELSE
         WRITE (IUNIT) (S(N), N=1,NEQ)
      ENDIF
C
      RETURN
      END
