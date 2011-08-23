C    CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
C
      SUBROUTINE SKABS
C
C     SUBROUTIKNE SKABS contains the introduction, comments, and
C     corrections for the CHEMKIN_III surface kinetics library of
C     subroutines.
C
C///////////////////////////////////////////////////////////////////
C
C            SKLIB: SURFACE KINETICS SUBROUTINE LIBRARY
C                   VERSION 7.9, 98/04/08
C
C     WRITTEN BY:
C         FRAN M. RUPLEY AND
C         ROBERT J. KEE
C         COMPUTATIONAL MECHANICS DEPARTMENT
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94550
C         (415) 294-3272
C       AND
C         MICHAEL E. COLTRIN
C         SANDIA NATIONAL LABORATORIES
C         CHEMICAL PROCESSING SCIENCES DEPARTMENT
C         ALBUQUERQUE, NM 87185
C         (505) 844-7843
C
C       Copyright 1990, Sandia Corporation.
C       The U.S. Goverment retains a limited license in this
C       software.
C
C       The U.S. Government retains, in this software, a paid-up,
C       nonexclusive, irrevocable worldwide license to reproduce,
C       prepare derivative works, perform publicly and display
C       publicly by or for the Government, including the right to
C       distribute to other Government contractors.
C
C       Neither the United States, the U.S. Dept. of Energy, nor
C       any of their employees, makes any warranty, express or
C       implied, or assumes any legal liability or responsibility
C       for the accuracy, completeness, or usefulness of any
C       information, apparatus, product, or process disclosed, or
C       represents that its use would not infringe privately owned
C       rights.
C
C/////////////////////////////////////////////////////////////////////
C      V. 7.9 98/04/08 E. Meeks
C      1. Action#0167: Separate out temperature-dependent rate
C         calculations from composition-dependent part.
C         a) Add subroutine SKRRPT and SKRRPX that, combined, perform
C            the same set of calculations as SKRROP.  SKRRPT only
C            does the temperature-dependent part, calculating RKFT and
C            RKRT; SKRRPX takes the RKFT and RKRT values and returns 
C            the RKR and RKF values.
C         b) In SKRROP, replace contents with call to SKRRPT, followed
C            by call to SKRRPX.
C         d) Add new user-callable subroutines SKKFRT, SKRATK, and
C            SKROPK.  The SKKFRT routine calls the SKRRPT subroutine,
C            returning RKFT and RKRT arrays.  The SKROPK routine
C            takes the RKFT and RKRT arrays as input, and returns the
C            net rates of progress of the reactions (similar to SKROP),
C            through a call to SKRRPX.  The SKRATK routine returns
C            species rates of production (similar to SKRAT) through a
C            call to SKRRPX.
C      V. 7.8 98/04/07 E. Meeks
C      1. Fix bug#165: In SKIRNU, replace variable name in SUBROUTINE
C         statement, NIRNU, with NIIRNU.  Also replace in comments of
C         the subroutine's abstract.
C      V. 7.7 97/11/07 E. Meeks
C      1. Fix bug#126: in SKRROP, correct treatment of site-specific
C         sticking reactions that use FORD change of order; was 
C         dividing by the total site density raised to the 
C         stoichiometric coefficient, rather than the revised order 
C         for the site species.
C      2. Fix bug#127: correct treatment of site-specific BOHM 
C         reactions that use FORD change of order in SKRROP, as above.
C      3. Correct indexing of RNU in NBOHM reaction loop in SKRROP.
c      V. 7.6 97/10/28 H. Moffat
c      1. Fix bug#031a: rewrite SKSYMR to correctly treat real
c         coefficients in reaction strings.
C      V. 7.5 97/08/28 F. Rupley
C      1. Fix bug#005b: in SKRROP ensure non-negative concentrations 
C         used for rate calculations with real stoichiometric 
C         coefficients (DO 350), and with real changed orders (DO 450)
C      V. 7.4 97/07/23 F. Rupley
C      1. Fix bug #012: remove spaces from logical operator . EQ. on
C         line 7376 in SKRROP. 
C      V. 7.3 97/06/10 (H. Moffat / E. Meeks)
C      1. Fixed Bug #032: Corrected indexing of real coeffient 
C         reaction numbers in subroutine SKCONT, loop 300;  
C         Added definition of I_IRNU=ISKWRK(IiIRNU),
C         and used this pointer in definition I=ISKWRK(I_IRNU+N).
C      V. 7.2 97/05/30 (H. Moffat)
C      1. Fix Bug #015: Corrected indexing on NUNK in SKSDOT;
C         was NUNK(xx,L) corrected to be NUNK(xx,I).  Affects real 
C         coefficient rate calculations that use SKSDOT.
C      2. Added two explanatory comments in the SKSDOT routine 
C      V. 7.1 97/04/15
C      1. in SKFLGS, ISKWRK(IiMOTZ) becomes ISKWRK(IiMSTK)
C         (CHEMKIN Committee bug report #004)
C      2. in several subroutines SAVE needs to precede DATA
C         (CHEMKIN Committee bug report #006)
C      V. 7.0 97/03/01
C      1. fix case of NTR=1, and other bugfixes;
C         "final" version for Reaction Design requirements.
C      V. 6.9 (97/01/21 F. Rupley)
C      1. in SKDSDC1 need loop over NIICOV reactions to modify
C         coverage reaction rates, due to the fact that reaction I
C         may appear more than once in the NIICOV reactions with
C         different coverage species and/or parameters.
C      V. 6.8 (96/12/11 F. Rupley)
C      1. several bugfixes
C      V. 6.7 (96/11/11 F. Rupley)
C      1. new time-saving logic in thermodynamic and rate routines
C      V. 6.6 (96/09/04 F. Rupley per H. Moffat)
C      1. new SKDSDC and SKSDOT from H. Moffat
C      V. 6.5 (96/08/19 F. Rupley, per M. Coltrin)
C      1. in SKRROP set default EQKC(I)=1.0, then check IF (SUMSMH.NE.0
C         to set EQKC(I)
C      V. 6.4 (96/08/05 F. Rupley)
C      1. minor edit in SKATCZ
C      V. 6.3 (96/05/24 F. Rupley)
C      1. SKIYLD argument NYDIM deleted.
C      2. SKIBHM argument IBMFL now reflects K, the ion species index
C         of a BOHM reaction.
C      V. 6.2 initial sccs version
C      VERSION 6.1 (4/29/96 F. Rupley per E. Meeks/M. Coltrin)
C      CAUTION...THIS WORK IS STILL IN PROGRESS!
C      1. in SKRROP rate RKFT for sticking reaction is MIN(RKFT,1.0)
C      CHANGES FOR VERSION 6.0
C      (Initial CHEMKIN-III version, F. Rupley, March 1, 1996)
C      1.  add workspace, pointers, and code for yield-modification
C          of stoichiometric coefficients.
C      2.  add, modify comments sections throughout in preparation
C          for updated documentation.
C      3.  add ascii linkfile options.
C      4.  separate linkfile version vs. code versions.
C      5.  PROLOGUEs
C      6.  more checking of change-order species to ensure required
C          coef=1.0 for certain species in some rate options
C      7.  MOTZ array of flags over the sticking reactions
C      8.  ETH energy threshold checking
C      9.  add yield-modify logic for SKRATI
C
C      CHANGES FOR VERSION 5.5 (E. Meeks)
C      1.  Added subroutine SKRPAR for auxiliary rate-parameter
C          information. The routine is used to initialize this
C          information by the user, rather than adding arguments
C          to SKRROP.  Contains ion energies for ion-energy dependent
C          reactions.
C      2.  Added pointer IrENGI to skstrt.h and SKINIT; required
C          expanding LENRSK and LENISK in skinterp.f
C      3.  Added argument RSKWRK(ISKWRK(IrENGI)) to SKRROP calls.
C      4.  Changed energy-dependence calculation in SKRROP, so that
C          reaction rate only changes if T_ion > E_threshold.
C      CHANGES FOR VERSION 5.4 (H. K. Moffat)
C      1.  Bug fixes to jacobian routines, due to SKRROP parameter
C          list additions.
C      2.  Fixed a bug in SKINIT - IiNEDP was spelled wrong
C      3.  Changed SKATCZ - It created jacobian entries.
C          This routine may need to be changed again.
C      4.  Some efficieny changes to SKRROP.
C      CHANGES FOR VERSION 5.3 (6/2/95 M. Coltrin)
C      1.  Add SKDSDC for the partial derivative of production rates
C          w.r.t. concentrations. (H. Moffat)
C      2.  Remove IBT array (Bohm temperatures); remove IiTBHM pointer
C          everywhere it occurred.
C      3.  Add pointers: IiIONS, IiKION, IiKTFL, IiNEDP, IiKEDP,
C          IiIEDP, IiELEC, IrPEDP.
C      4.  Read KION array from linking file.
C      5.  Read IEDP, KEDEP, PEDEP arrays from linking file.
C      6.  Add NEDPAR, NKKION, KEL, NIIEDP to list of variables
C          read from linking file.
C      7.  Add the following arguments to SKRROP: NIIEDP, IEDP,
C          KEDEP, PEDEP, KTFL, NEDPAR, KEL
C      8.  Make temperature an array in every subroutine that uses it:
C          SKAML, SKAMS, SKATCZ, SKCPML, SKCPMS, SKCPOR, SKDEN, SKDRDA,
C          SKDRDC, SKEQ, SKGML, SKGMS, SKHML, SKHMS, SKHORT, SKRAT,
C          SKROP, SKRROP, SKSMH, SKSML, SKSMS, SKSOR, SKUML, SKUMS.
C      9.  Added SUBROUTINE SKIENR to return an integer flag if the
C          input reaction number has an ion-energy dependence.
C      10. Add SUBROUTINE SKKTFL to allow the user to set the location
C          in the temperature array to be used for each gas-phase
C          species.
C      11. Reset NLIST to 1, because the linking file has changed, and
C          made the only item in LIST '4.7' to correspond to the new
C          interpreter version (these changes are in SKLEN).
C      12. Add SUBROUTINE SKKION to return the species index for the
C          the electron, the number of positive ions in the gas-phase
C          and a array containing their species numbers.
C      13. Removed all PSK... subroutines (related to plasmas) because
C          they are now obsolete.
C      CHANGES FOR VERSION 5.2 (2/27/95 F. Rupley)
C      1.  Change character index "(:" to "(1:"
C      2.  Correct SKRROP to get rid of requirement to have specified
C          sites (per E. Meeks)
C      CHANGES FOR VERSION 5.0 (1/19/95 F. Rupley per M. Coltrin)
C      1.  Add integer error flag to SKLEN, SKINIT call lists)
C      CHANGES FOR VERSION 4.9 (12/15/94 F. Rupley per E. Meeks)
C      1.  Correct value of TEMP in PSKHML.
C      CHANGES FOR VERSION 4.8 (9/30/94 F. Rupley per G. Evans)
C      1.  Correct SKRROP (DO 70); variable SDOT should be SDTOT
C      CHANGES FOR VERSION 4.7 (9/7/94 F. Rupley per M. Coltrin)
C      1.  Loop 10 in SKEQ contained a wrong pointer (IrEQ should have
C          been IRIT3).  This subroutine has not returned the
C          correct value since version 4.2.
C      CHANGES FOR VERSION 4.52 (8/26 F. Rupley)
C      1.  Correct above value for RUC (RSKWRK(ISKWRK(IrRUC)))
C          in SKINIT.
C      CHANGES FOR VERSION 4.51 (8/10/94 H. Moffat)
C      1.  Changed gas constant to conform to 1986 CODATA
C          recommendations. RUC, obtained from RU by dividing
C          by 4.184, is now compatible with RU up to DP machine
C          precision.
C      CHANGES FOR VERSION 4.5 (8/6/94 F. Rupley)
C      1.  Assign local variables before CALL SKRROP in order to
C          shorten call list.
C      CHANGES FOR VERSION 4.4 (8/4/94 F. Rupley)
C      1.  Correct READ statement in SKINIT for 'IF (NIICOV'.
C      2.  Correct 'IND =' index in SKABE
C      3.  Use REWIND in SKINIT if last material in linking file
C      CHANGES FOR VERSION 4.3 (8/2/94 F. Rupley)
C      1.  Bugfixes regarding using old constants vs. newer
C          ISKWRK pointers
C      2.  Add pointers to SKSTRT in order to store work array lengths
C      3.  Accept linking file V.4.3.
C      CHANGES FOR VERSION 4.2 (7/27/94 F. Rupley per E. Meeks)
C      1.  Correction NIBOHM > NIIBHM; accept linking file 4.2.
C      2.  Correct DO 70 in SKRROP for BOHM rates.
C      CHANGES FOR VERSION 4.14 (7/14/94 F.Rupley)
C      1.  Changes to allow multiple materials include:
C          pointers are now stored in ISKWRK
C          most subroutines now require ISKWRK input
C          backspaces instead of rewind in SKLEN
C          etc.
C      CHANGES FOR VERSION 4.13 (6/28/94 E. Meeks)
C      1.  Add subroutine SKIBHM and correct Bohm reaction rate formula
C      2.  Correct common block /SKSTRT/ to include IiIBHM
C      3.  Added option for correction to sticking parameter; Option
C          is turned off in the interpreter input using REACTION line
C          keyword MWOFF and passed through linking file to SKRROP.
C      CHANGES FOR VERSION 4.12c (6/3/94 F. Rupley)
C      1.  Accept linking file 4.08c (bugfixes per H. Moffat)
C      CHANGES FOR VERSION 4.12b (5/20/94 F. Rupley per E. Meeks)
C      1.  Incorporate plasma options (linking file 4.08b)
C      CHANGES FOR VERSION 4.12 (4/28/94 F. Rupley, per M. Coltrin)
C      1.  New Subroutines SKINU, SKIRNU, SKIORD for real
C          stoichiometric coefficients and change of order.
C      CHANGES FOR VERSION 4.11 (4/19/94 F. Rupley)
C      1.  accept linking file 4.07 (correction in SKBAL, SKRBAL)
C      CHANGES FOR VERSION 4.10 (4/14/94 F. Rupley, per E. Meeks)
C      1.  use INCLUDE 'skstrt.h' instead of having COMMON /SKSTRT/ in
C          every subroutine
C      CHANGES FOR VERSION 4.09 (3/15/94 F. Rupley)
C      1.  DOS/PC compatibility effort includes adding file names to
C          OPEN statements, removing unused variables in CALL lists,
C          unusued but possibly initialized variables.
C       CHANGES FOR VERSION 4.08 (1/27/94 F. Rupley per R. Kee)
C       1. Real stoichometric coefficients added; NRNU,  IRNU,
C          RNU (additional pointers)
C       2. Integer phase (site) balance INCF becomes real RNCF;
C          require RSKWRK be added to calls to SKNCON, SKNU
C       CHANGES FOR VERSION 4.07 (10/1/92 F. Rupley per M. Coltrin)
C       1. COMMON /SKCONS/ VERS, PREC, KERR, LENI, LENR, LENC
C          eliminates the need for LINSK argument to SKSAVE
C       CHANGES FOR VERSION 4.06 (4/13/92 F. Rupley)
C       1. Accept Interpreter binary file V.4.04. (correction to
C          Subroutine SKDUP)
C       CHANGES FOR VERSION 4.04 (8/28/91 F. Rupley per J. Grcar)
C       1. Instead of using NSPAR+1 in the calls to SKROP for the
C          dimension of SPAR, pass NSPAR, then dimension
C          SPAR(NSPAR+1,*) in SKROP, to avoid confusion as to the
C          value of NSPAR.
C       CHANGES FOR VERSION 4.03 (8/1/91 F. Rupley per M. Coltrin)
C       1. New SKSTRT pointer IrEQ added to record read in SKPNT and
C          added to record written in SKSAVE.
C       CHANGES FOR VERSION 4.02 (7/17/91 F. Rupley)
C       1. Accept Interpreter binary file V.4.02.
C       CHANGES FOR VERSION 4.01 (7/11/91 F. Rupley)
C       1. Accept Interpreter binary file V.4.01.
C       CHANGES FOR VERSION 4.0 (6/5/91 M. Coltrin)
C       1. Add new subroutines SKFLGS, SKIREV, SKNUF (H. Moffat)
C       2. Change units on surface-coverage modification of the the
C          rate of progress (concentration units --> site fractions)
C       3. Fix error in calculation of equilibrium constant. Added
C          array of length IISUR to hold correction factor for each
C          reaction. Added pointer IrEQ to common block SKSTRT.
C          Added EQFAC to arguments of SKRROP.
C       4. In SKATCZ and SKDEN we were dividing a real number by
C          ISKWRK(IiNSCV). Put in change blocks to first convert
C          the number to floating point (either REAL or DBLE).
C       CHANGES FOR VERSION 3.87 (5/23/91 H. Moffat)
C       1. In SKRROP, check to see if sticking coefficient is greater
C          than one. If it is, set it to one.
C       CHANGES FOR VERSION 3.86 (5/9/91 F. Rupley)
C       1. In SKRROP, perform the "d" (PAR(4,I)) perturbation
C          before the checking for change of sign of the rates instead
C          of after.
C       CHANGES FOR VERSION 3.85 (4/29/91 F. Rupley per M. Coltrin)
C       1. Correct indexing in SKRDEX
C       CHANGES FOR VERSION 3.84 (4/2/91 F. Rupley)
C       1. Add Subroutine SKRHEX to get/put thermodynamic polynomial
C          coefficient a6, to enable applications codes to perturb
C          heat of formation for species.
C       2. Add Subroutine SKMAXTP to find the maximum number of
C          temperatures used in fitting thermodynamic data for the
C          species.
C       3. Change Subroutine SKABE to return either Arhennius or
C          sticking coefficients for the reactions, and an array
C          of integer flags to indicate the type of the coefficients.
C      CHANGES FOR VERSION 3.83 (2/18/91 F. Rupley, per R. Kee)
C      1.  Modify equation using sticking parameters for rates in
C          SKRROP, per Peter Glarborg.
C       2. Add a fourth parameter to the array of Arhennius
C          coefficients for the IISUR reactions;
C          increase the value of NSPAR in COMMON /SKSTRT/ by one.
C          (this also increases the length of the array of reverse
C          Arhennius parameters);
C          initialize the value of the fourth parameter to 1.0 in
C          SKINIT;
C          use this value as a "perturbation factor" for the forward
C          rates in SKRROP;
C          add SUBROUTINE SKRDEX to allow applications codes to change
C          the perturbation factor RD(I) in sensitivity calculations.
C       3. Accept binary file V.3.78 (LENRSK was increased by NIISUR
C          + NIIREV to reflect above changes in RSKWRK array.
C      CHANGES FOR VERSION 3.82 (F. Rupley)
C      1.  Accept surface binary file V.3.77 (correction)
C      CHANGES FOR VERSION 3.81 (F. Rupley)
C      1.  Accept surface binary file V.3.76 (correction)
C      CHANGES FOR VERSION 3.80 (F. Rupley, per M. Coltrin 1/11/91)
C      1.  Add COMMON MACHN and set values of BIG, SMALL, EXPARG in
C          Subroutine SKPNT.
C      CHANGES FOR VERSION 3.79 (F. Rupley, per M. Coltrin 1/10/91)
C      1.  Expand logic in SKCOMP to ignore leading/trailing spaces
C          in strings being compared.
C      CHANGES FOR VERSION 3.78 (F. Rupley, per M. Coltrin 1/3/91)
C      1.  Correct bug in Subroutine SKNU
C      CHANGES FOR VERSION 3.77
C      1.  Accept Surface binary file 3.75 (correction)
C      CHANGES FOR VERSION 3.76
C      1.  Accept Surface binary file 3.74 (correction made to
C          interpreter)
C      CHANGES FOR VERSION 3.75
C      1.  Accept Surface binary file 3.73 (correction made to
C          interpreter)
C      CHANGES FOR VERSION 3.74
C      1.  Correction to SKRROP
C      CHANGES FOR VERSION 3.73
C      1.  Accept Surface binary file 3.72 (correction made to
C          interpreter)
C      CHANGES FOR VERSION 3.72
C      1.  In SKRROP, if SDTOT=0 then RKF(I)=0.
C      CHANGES FOR VERSION 3.71
C      1.  Accept Surface binary file 3.71.
C     CHANGES FOR VERSION 3.7
C      1.  SKNCON requires ISKWRK input, and returns array of
C          length NPHASE
C      2.  SKNU returns array KSTOIC, stoichiometric coefficients
C          for the surface reactions, and array NSTOIC, the net change
C          in phases for the surface reactions.
C      3.  SKSDEN return array SDEN changed to SDEN0
C      4.  Correct indexing for ISKWRK(IiNCF) in SKINIT, SKRAT and
C          SKRATI
C      5.  Initialize NSIG=0 in Subroutine SKRROP after the start of
C          the DO 80 loop, not before it.
C     CHANGES FOR VERSION 3.64
C      1.  Accept Surface binary file 3.64
C     CHANGES FOR VERSION 3.63
C      1.  Accept binary file versions 3.62 and 3.63
C      2.  Modify SKRROP according to new Equation 35, need
C          NFSUR, NLSUR, KFIRST, KLAST, SDEN, and RNCF in call list.
C     CHANGES FOR VERSION 3.62
C      1.  New versions of SKDRDA and SKDRDC by M. Coltrin.
C      2.  New SKDRDA reflects reversal of change #10 above.
C      3.  SKDRDA call list changed from (IR, ROP, ISKWRK, ...) to
C          (IR, P, T, ACT, SDEN, ISKWRK, ...)
C      CHANGES TO VERSION 3.61
C      1.  Correct SKABE conversion of sticking coefficients
C      2.  Change SKDEN to return densities of surface species
C          in gm/cm**2 instead of moles/cm**2
C      3.  Check that there site-phase species, bulk-phase species,
C          sites, or site phases before filling appropriate arrays.
C      4.  Correct SUMS in SKABE and SKRROP to include only
C          stoichiometric coefficients of Surface Species, and
C          check that SUMS is greater than 0.0
C      5.  SIGN function in SKRROP to preserve sign of rates
C      6.  Add VERS 3.61 to list of correct versions in SKLEN
C      7.  Correct DO 50 loop in SKICOV
C      8.  DO 115 in SKABE  only if NNSUR.GT.0
C          DO  60 in SKDRDC  "       "
C          DO  60 in SKEQ    "       "
C          DO  20 in SKRAT   "       "
C             (and combine two loops over NNSUR into one loop)
C          DO  60 in SKROP   "       "
C      9.  In SKSDEN, if NNSUR.LE.0 RETURN
C     10.  In SKDRDA, CALL SKABE instead of SKRAEX to get RA for
C          reaction IR.
C     11.  In SKDRDC, ROPNU is a sum if species present as reactant
C          and as product
C     12.  In SKRAT, replace the variable SDOT with SITDOT,
C          then replace WDOT with SDOT
C     13.  In SKRATI, replace the variable WDOTI with SDOTI, and
C          add the array SITDTI (phase change due to a reaction)
C     14.  Add pointer IrPT1 to COMMON/SKSTRT for real dummy work
C          space of length NPHASE, in order to accept the array
C          SITDTI(*) where SUBROUTINE SKDRDA has a CALL SKRATI
C     15.  In SKSYMR, include slash-delimited phase name in a
C          reaction string where a species name is not unique
C     16.  New subroutine SKPNT, SKSAVE to read, write binary
C          file information and work arrays.
C     17.  SKSYMR corrected to start the reaction string in the
C          first character space.
C     CHANGES TO VERSION 3.6
C      1.  Bring up to date with manual changes
C      2.  SKCOMP has additional argument NT to indicate how many
C          occurrences of a character string occurs in an array
C          of character strings.
C      3.  New Subroutine SKPCMP to find and species index number
C          for a SPECNAME/PHASENAME/ character string.
C      4.  Subroutine SKABE must convert sticking
C          coefficients to Arrhenius coefficients
C      5.  Subroutine SKSNUM
C      6.  Subroutine SKISTK
C      7.  User subroutine SKRTI renamed to SKRATI, and internal
C          subroutine SKRATI becomes SKRTI
C      8.  Eliminate SKRTI, since same logic is in SKRATI; the
C          only subroutine that called SKRTI is SKDRDA; SKDRDA
C          now calls SKRATI instead.
C      9.  Correct DO 20 loop in SKRAT to go from NFSUR to NLSUR
C          instead of 1 to NNSUR
C      CHANGES FOR VERSION 3.5
C      1.  Additional pointer IiKCOV and modified read for the species
C          numbers of the surface species associated with the coverage
C          parameters in a surface reaction; additional subroutine
C          SKICOV provides information for surface reaction I as to
C          it total number of coverage species, their species numbers,
C          and their coverage parameter numbers.
C      2.  Modified SKRAT to correct SDOT loop
C      3.  Modified SKDEN to correct units and calculation of site
C          species densities
C      4.  Modified SKRROP for coverage factors in calculations of
C          forward rates.
C      5.  Subroutine SKLEN reads first record of surface binary
C          file in order to provide lengths required for work arrays.
C      6.  Removed arrays KFIRST, KLAST, and KKPHAS from argument
C          list of Subroutine SKINDX.
C      7.  Change Subroutine SKPKK to return arrays KFIRST, KLAST,
C          and KKPHAS for all the phases, instead of just for
C          a specified phase.
C      8.  Additional input SDEN required for Subroutine SKDEN.
C      9.  Added integer constant NIICON to first record of binary
C          file, to COMMON/SKSTRT/, and to Subroutine ???,
C           for the total number of reactions which do not conserve
C           sites.
C     10.  No longer need argument ISKWRK for SKINDX.
C     11.  New Subroutine SKCONT
C     12.  First record of binary file contains character strings
C          VERS and PREC to indication binary file version and
C          precision, and the logical flag KERR to indicate an
C          error in the binary file; KERR previously was the
C          first item in the binary file
C     13.  Correct SKRAT DO 50, SKRATI DO 100, and SKRDWC DO 400 to
C          make sure NUNK(N,I) .ne. 0.
C     14.  Delete SKRPAR, since somewhat duplicates SKABE.
C     15.  Change SKKEL to SKNCF to parallel Chemkin.
C     16.  Change SKHRT to SKHORT     "        "    .
C     17.  CHANGE SKTHM TO SKATHM     "
C     18.  Delete SKRDWC and include it's calculation in SKDRDC,
C          since SKDRDC is the only subroutine that calls SKRDWC.
C     19.  COMMON /MACH/ renamed to COMMON /MACHN/ to avoid conflict
C          with CKLIB
C     20.  Additional factor in forward rate calculation
C     CHANGES FOR VERSION 3.4
C       1. Remove pointer and read statement for RSKWRK(IrSTK); real
C          sticking coefficients are now stored in the Arrhenius
C          parameter work space.
C       2. Site densities now have units of moles/cm**2 instead of
C          # densities.
C       3. SKRROP modified (DO 120 loop) for 'sticking' coefficients
C          of reaction I:
C          If equivalent first Arrhenius coefficient is
C             PAR1 = PAR(1,I) * 3637.601 / SQRT(WT1),
C             where WT1 is the molecular wt. of the gas-phase reactant,
C          and equivalent second Arrhenius coefficient is
C             PAR2 = 0.5 + PAR(2,I),
C          then the forward rate expression
C             RKF(I) = PAR(1,I) * EXP(PAR(2,I)*ALOGT - PAR(3,I)/T)
C          becomes
C             RKF(I) = PAR1 * EXP(PAR2*ALOGT - PAR(3,I)/T) / (SDTOT**M)
C          where SDTOT is the total of the site densities (moles/cm**2),
C          and M is the sum of the stoichiometric coefficients of the
C          site-phase (surface) reactants.
C          This affects all subroutines that CALL SKRROP, since SDTOT,
C          WT(*), NIISTK, IISTK(*), NKKGAS and NKKSUR must be known.
C      4.  Subroutine SKATCZ and all subroutines which call it now
C          require additional input argument SDEN, the site densities.
C      5.  Subroutine SKRATS is renamed SKRAT to replace previous SKRAT
C          and returns additional array SDOT.
C      6.  Correct SKATCZ: divide by coverage parameter rather than
C          multiply
C     CHANGES FOR VERSION 3.3
C       1. CZ(K) for site species multiplied by their site coverage
C                parameter
C       2. Add sticking coefficient, and 'NONCON' options to allow
C          non-balancing site information for reactions;
C          RNCF(N,I), N=1,NPHASE is site information for the species
C          in reaction I
C     CHANGES FOR VERSION 3.2
C       1. Phase pointer for SDEN was incorrect in loop 210 of SKDEN
C          and loop 130 of SKATCZ
C     CHANGES FOR VERSION 3.1
C       1. Mixture/pure species eliminated in favor of bulk species
C          only; binary file and common skstrt restructured
C     CHANGES FOR VERSION 3.1
C       1. vax/cray change blocks for machine constants changed to
C          smallexp/bigexp blocks
C     CHANGES FOR VERSION 3.0
C       1. Implement phases, bulk species
C     CHANGES FOR VERSION 2.5
C       1. Implement new mixture and pure solids.
C     CHANGES FROM VERSION 2.0
C       1. Get rid of non-standard comments
C     CHANGES FROM VERSION 1.9
C       1. Character manipulation subroutines have additional
C          error handling and arguments
C       2. Change IckNAME to IckNAM in SKDNAM
C     CHANGES FROM VERSION 1.8
C       1. Changed DDOT/SDOT to DO loops
C       2. Allow reversible Arhennius coefficients, with changes
C          to /SKSTRT/ of NIIREV, the number of reactions with
C          reverse Arrhenius parameters defined, ISKWRK pointer IiIREV
C          for the array of reaction numbers, RSKWRK pointer IrRPAR
C          for the array of reverse parameters, and additional
C          reaction workspace pointer IrIT3
C       3. New SKRROP to allow reverse Arrhenius parameters.
C     CHANGES FROM VERSION 1.7
C       1. Added PROLOGUE modifications
C       2. Added subroutine SKNU, SKSTOI and SKCOVI
C       3. Eliminate SKSEL
C       4. Switch old SKRAT  to SKRAT3 and add SKRAT,
C          switch old SKRATI to SKRTI3 and add SKRTI,
C          switch old SKRTI  to SKRATI
C       5. Added subroutines SKEQYP and SKYTCZ
C       6. Correct call list for CALL SKSTOI in subroutine SKNU
C       7. Added DO 20 loop to subroutine SKRROP
C       8. Pointer in 2nd call to DCOPY in SKYTCZ
C       9. Eliminated SKWT, since essentially a duplicate of SKKWT
C      10. Added missing DIMENSION ISKWRK(*) to subroutine SKSITP
C     CHANGES FROM VERSION 1.6:
C       1. Added thermodynamic properties subroutines -
C             SKTHM3  SKHML3  SKHMS3 SKHRT3
C             SKAML   SKAML3  SKAMS  SKAMS3
C             SKCPML  SKCPL3  SKCPMS SKCPS3  SKCPOR
C             SKGML   SKGML3  SKGMS  SKGMS3
C             SKSML   SKSML3  SKSMS  SKSMS3  SKSMH  SKSMH3 SKSOR
C             SKUML   SKUML3  SKUMS  SKUMS3
C     CHANGES FROM VERSION 1.5:
C       1. Correct call lists in some subroutines using NKK...'s as
C          an argument instead of defined in COMMON /SKSTRT/
C       2. Add double precision change blocks with DBLE(I) to
C          parallel REAL(I) used in calculations
C       3. Replace Subroutine SKRTK
C     CHANGES FROM VERSION 1.4:
C       1. Added thermodynamic properties subroutines SKTHM, SKHML,
C          SKHMS, SKHRT
C       2. Changed from SCOPY for integer arrays to DO loops
C       3. Changed order of call list for SKWT, SKKNAM
C     CHANGES FROM VERSION 1.3:
C       1. Added subroutines SKDRDA, SKDRDC, and SKRDWC
C     CHANGES FROM VERSION 1.2:
C       1. Replaced "REAL*8" statements with "DOUBLE PRECISION"
C     CHANGES FROM VERSION 1.1:
C       1. Changed from individual site densities to SUMDEN in
C          Subroutines SKRROP and SKRAT
C     CHANGES FROM VERSION 1.0:
C       1. Added subroutine SKRAEX
C/////////////////////////////////////////////////////////////////////
C
C  START PROLOGUE
C  SUBROUTINE SKABS
C
C  Work arrays ISKWRK, RSKWRK, and CSKWRK contain information about the
C  elements, species and reaction in the mechanism; they also contain
C  some work space needed for internal manipulations.  A user wishing
C  to modify a subroutine or to write new routines will probably want
C  to use the work arrays directly.  The pointers described below are
C  starting addresses for information stored in the work arrays, and
C  are found in the labeled common block COMMON /SKSTRT/, declared by
C  the use of the include file skstrt.h.
C
C  It should be noted that storage in work array ISKWRK serves
C  several purposes, for instance,
C  1. ISKWRK(N) may be a constant in regards to a surface described
C     by the mechanism, and is subject to change if there are multiple
C     surfaces,
C  2. ISKWRK(N) may be a pointer, where L=ISKWRK(N) is used to locate
C     locate information specific to a surface, either in
C     ISKWRK, RSKWRK, or CSKWRK, then
C  3. ISKWRK(L=ISKWRK(N)) may the first value of an integer array
C     specific to a surface, RSKWRK(L=ISKWRK(N)) the first of a real
C     array, and CSKWRK(L=ISKWRK(N)) the first of a character string
C     array.
C
C  COMMON /SKSTRT/
C
C  Integer constants
C
C 1   MAXSPR, NELEM, NKKGAS, NSPAR, NSCOV, NEDPAR, NYPAR, MAXORD,
C 2   MAXTP, NCP,    NCP1,  NCP2,   NCP2T,
C
C  ISKWRK pointers to integer variables
C
C 3   IiLENI, IiLENR, IiLENC, IiKSUR, IiKBLK, IiKTOT, IiNPHA, IiFSUR,
C 4   IiLSUR, IiNSUR, IiFBLK, IiLBLK, IiNBLK, IiNIIS, IiNCOV, IiNREV,
C 5   IiNSTK, IiNCON, IiNBHM, IiNRNU, IiNORD, IiNEDP, IiELEC, IiNYLD,
C
C  ISKWRK pointers to integer arrays
C
C 6   IiPKST, IiPKND, IiPTOT, IiKPHS, IiKCHG, IiKCMP, IiNSCV, IiKNT,
C 7   IiNRPP, IiNREA, IiNUNK, IiNU,   IiNSUM, IiICOV, IiKCOV, IiIREV,
C 8   IiISTK, IiMSTK, IiIBHM, IiKBHM, IiIRNU, IiIORD, IiKORD, IiIONS,
C 9   IiKION, IiKTFL, IiNEDP, IiIEDP, IiKEDP, IiIYLD, IiYION, IiKYLD,
C
C  ISKWRK pointers to real variables
C
C *   IrSKMN, IrPATM, IrRU,   IrRUC,
C
C  ISKWRK pointers to real arrays
C
C 1   IrSDEN, IrKTMP, IrKTHM, IrKDEN, IrAWT,  IrKWT,  IrPAR, IrKCOV,
C 2   IrRPAR, IrEQ,   IrRNU,  IrNCF,  IrKORD, IrKFT,  IrKRT, IrKT1,
C 3   IrKT2,  IrPT1,  IrIT1,  IrIT2,  IrIT3,  IrPEDP, IrENGI,
C 4   IrPYLD, IrYNCF,
C
C  ISKWRK pointers to character string arrays
C
C 4   IcENAM, IcKNAM, IcMNAM, IcPNAM
C
C  END include file for sklib.f
C
C  INTEGER CONSTANTS:
C
C  MAXSPR, Maximum number of species in any surface reaction.
C          Unless changed in the interpreter MAXSPR=12.
C  NELEM,  Total count, elements in problem.
C  NKKGAS, Total count, gas-phase species in problem.
C  NSPAR,  Number of parameters  required in the rate expression
C          for reactions;  in the current formulation NSPAR=3,
C          however, a 4th parameter can be used for purposes of scaling.
C  NSCOV,  Number of parameters required in the rate expression
C          for a coverage reaction;  NSCOV=3.
C  NEDPAR, Number of parameters required in the rate expression
C          for ion-energy dependence reactions;  NEDPAR=3.
C  NYPAR,  Number of parameters required in the yield-modified
C          reactions; NYPAR=4.
C  MAXORD, Maximum number of change-orders allowed in a reaction.
C  MAXTP,  Maximum number of temperatures allowed in fits of
C          thermodynamic properties for any species; MAXTP=3.
C  NCP,    Number of polynomial coefficients to fits of CP/R for
C          a species; NCP=5.
C  NCP1,   NCP+1.
C  NCP2,   NCP+2.
C  NCP2T,  Total number of thermodynamic fit coefficients for
C          species; NCP2T = (MAXTP-1)*NCP2 = 14.
C
C  ISKWRK POINTERS TO INTEGER VARIABLES:
C
C  IiLENI, ISKWRK(IiLENI) is the total length of ISKWRK required.
C  IiLENR, ISKWRK(IiLENR) is the total length of RSKWRK required.
C  IiLENC, ISKWRK(IiLENC) is the total length of CSKWRK required.
C  IiKSUR, ISKWRK(IiKSUR) is the total surface species count.
C  IiKBLK, ISKWRK(IiKTOT) is the total species count (gas+surface+bulk).
C  IiNPHA, ISKWRK(IiNPHA) is the total phase count.
C  IiFSUR, ISKWRK(IiFSUR) is the phase index of the first site.
C  IiLSUR, ISKWRK(IiLSUR) is the phase index of the last site.
C  IiNSUR, ISKWRK(IiNSUR) is the total surface phage count.
C  IiFBLK, ISKWRK(IiFBLK) is the phase index of the first bulk.
C  IiLBLK, ISKWRK(IiLBLK) is the phase index of the last bulk.
C  IiNBLK, ISKWRK(IiNBLK) is the total bulk phase count.
C  IiNIIS, ISKWRK(IiNIIS) is the total surface reaction count.
C  IiNCOV, ISKWRK(IiNCOV) is the total coverage reaction count.
C  IiNREV, ISKWRK(IiNREV) is the total count of surface reactions which
C          use explicit reverse parameters.
C  IiNSTK, ISKWRK(IiNSTK) is the total count of sticking surface
C          reactions.
C  IiNCON, ISKWRK(IiNCON) is the total count of surface reactions which
C          do not conserve sites.
C  IiNBHM, ISKWRK(IiNBHM) is the total count of Bohm surface reactions.
C  IiNRNU, ISKWRK(IiNRNU) is the total count of surface reactions with
C          real stoichiometry coefficients.
C  IiNORD, ISKWRK(IiNORD) is the total count of surface reactions with
C          changed-order species.
C  IiNEDP, ISKWRK(IiNEDP) is the total count of surface reactions with
C          ion-energy dependence.
C  IiELEC, ISKWRK(IiELEC) is the location in a species array of the
C          electron species.
C  IiNYLD, ISKWRK(IiNYLD) is the total count of surface reactions with
C          yield-modified species.
C
C  ISKWRK POINTERS TO THE START OF ISKWRK ARRAY WORKSPACE:
C
C  IiPKST, ISKWRK(I = ISKWRK(IiPKST)) starts an array of species indices
C          for the first species of the phases;
C          ISKWRK(I + N - 1) is the index of the first species of
C          phase N.
C  IiPKND, ISKWRK(I = ISKWRK(IiPKND)) starts an array of species indices
C          for the last species of the phases;
C          ISKWRK(I + N - 1) is the indes of the final species of
C          phase N.
C  IiPTOT, ISKWRK(I = ISKWRK(IiPTOT)) starts an array of total counts
C          of species in the phases;
C          ISKWRK(I + N - 1) is the total species count of phase N.
C  IiKPHS, ISKWRK(I = ISKWRK(IiKPHS)) starts an array of physical
C          phases for the species;
C          ISKWRK(I + K - 1)  = -1, species K is a solid,
C                             =  0, speciies K is a gas,
C                             = +1, species K is a liquid.
C  IiKCHG, ISKWRK(I = ISKWRK(IiKCHG)) starts an array of electronic
C          charges for the species;
C          ISKWRK(I + K - 1) is the charge of species K,
C          for example, a value of -2 indicates two excess electrons.
C  IiKCMP, ISKWRK(I = ISKWRK(IiKCMP)) starts a matrix of elemental
C          composition for the species;
C          ISKWRK(I + (K-1)*NELEM + M - 1) is the quantity of
C          element M in species K.
C  IiNSCV, ISKWRK(I = ISKWRK(IiNSCV)) starts an array of site coverage
C          for the species;
C          ISKWRK(I + K - 1) is the site coverage of species K.
C  IiKNT,  ISKWRK(I = ISKWRK(IiKNT)) starts an array of the total
C          number of temperatures dividing the ranges of thermodynamic
C          fits of the species.
C          ISKWRK(I + K - 1) is the number of dividing temperatures
C          for thermodynamic fits for species K;
C  IiNRPP, ISKWRK(I = ISKWRK(IiNRPP)) starts an array of the total
C          of participant species for the surface reactions,
C          and indidates the reversibility of the reactions;
C          ISKWRK(I + IS - 1) = +N, reaction IS has N participant
C                                   species (reactants + products),
C                                   and reaction IS is reversible,
C                             = -N, N participant species in an
C                                   irreversible reaction.
C  IiNREA, ISKWRK(I = ISKWRK(IiNREA)) starts an array of the total
C          count of reactants only for the surface reactions.
C          ISKWRK(I + N - 1) is the reactant count for the Nth
C          surface reaction.
C  IiNUNK, ISKWRK(I = ISKWRK(IiNUNK)) starts a matrix of indices
C          for the species in the surface reactions;
C          ISKWRK(I + (N-1)*MAXSPR + L - 1) is the species index for
C          Lth species in the Nth surface reaction.
C  IiNU,   ISKWRK(I = ISKWRK(IiNU)) starts a matrix of stoichiometric
C          coefficients for the species in the surface reactions;
C          ISKWRK(I + (N-1)*MAXSPR + L - 1) is the stoichiometric
C          of the Lth species in the Nth surface reaction.
C  IiNSUM, ISKWRK(I = ISKWRK(IiNSUM)) starts an array containing sums
C          of the stoichiometric coefficients of the (gas-phase only)
C          species in the surface reactions;
C          ISKWRK(I + N - 1) is the sum of gas-phase species
C          stoichiometric coefficients for the Nth surface reaction.
C  IiICOV, ISKWRK(I = ISKWRK(IiCOV)) starts an array of reaction
C          indices for those with coverage parameters;
C          ISKWRK(I + N - 1) is the reaction index of the Nth reaction
C          with coverage parameters.
C  IiKCOV, ISKWRK(I = ISKWRK(IiKCOV)) starts an array of species
C          indices for the coverage-dependendent species in
C          the surface reactions with coverage parameters;
C          ISKWRK(I + N - 1) is the coverage-dependent species for the
C          Nth reaction with coverage parameters.
C  IiIREV, ISKWRK(I = ISKWRK(IiIREV)) starts an array of reaction
C          indices for those with explicit reverse Arrhenius parameters;
C          ISKWRK(I + N - 1) is the reaction index of the Nth reaction
C          with explicit reverse parameters.
C  IiISTK, ISKWRK(I = ISKWRK(IiISTK)) starts an array of reaction
C          indices for those with sticking coefficients;
C          ISKWRK(I + N - 1) is the reaction index of the Nth sticking
C          coefficient reaction.
C  IiMSTK, ISKWRK(I = ISKWRK(IiMSTK)) starts an array of 0/1 flags
C          for Motz-Wise rate correction for the surface reactions
C          with sticking coefficients;
C          ISKWRK(I + N - 1) = 1, use Motz-Wise rate correction for
C          the Nth sticking reaction.
C  IiIBHM, ISKWRK(I = ISKWRK(IiIBHM)) starts an array of reaction
C          indices for the Bohm surface reactions;
C          ISKWRK(I + N - 1) is the reaction index of the Nth Bohm
C          reaction.
C  IiKBHM, ISKWRK(I = ISKWRK(IiKBHM)) starts an array of species
C          indices used in the Bohm formulation;
C          ISKWRK(I + N - 1) is the index of the species used for the
C          Nth Bohm reaction.
C  IiIRNU, ISKWRK(I = ISKWRK(IiIRNU)) starts an array of reaction
C          indices for those with real stoichiometric coefficients;
C          ISKWRK(I + N - 1) is the reaction index of the Nth real
C          stoichiometry reaction.
C  IiIORD, ISKWRK(I = ISKWRK(IiIORD)) starts an array of reaction
C          indices for those with changed-order species;
C          ISKWRK(I + N - 1) is the reaction index of the
C          Nth change-order reaction.
C  IiKORD, ISKWRK(I = ISKWRK(IiKORD)) starts a matrix of changed-
C          species indices for the change-order reactions;
C          ISKWRK(I + (N-1)*MAXORD + L - 1) is the index of the Lth
C          changed-order species in the Nth change-order reaction.
C  IiIONS, ISKWRK(I = ISKWRK(IiIONS)) starts an array of species
C          indices for the ionic species.
C          ISKWRK(I + N - 1) is the species index of the Nth ionic
C          species.
C  IiKTFL, ISKWRK(I = ISKWRK(IiKTFL)) starts an array of indices
C          into a temperature array for the species;
C          ISKWRK(I + K - 1) is the temperature index for species K.
C  IiIEDP, ISKWRK(I = ISKWRK(IiIEDP)) starts an array of reaction
C          indices for the surface reactions with ion-energy dependence.
C          ISKWRK(I + N - 1) is the species index of the Nth ion-
C          -energy-dependent reaction.
C  IiKEDP, ISKWRK(I = ISKWRK(IiKEDP)) starts an array of species
C          indices for the energy-dependent ions in the
C          ion-energy-dependence reactions;
C          ISKWRK(I + N - 1) is the species index of the ion for the Nth
C          energy-dependent reaction.
C  IiIYLD, ISKWRK(I = ISKWRK(IiIYLD)) starts an array of reaction
C          indices for those with modified yield;
C          ISKWRK(I + N - 1) is the reaction index of the Nth yield-
C          modified reaction.
C  IiYION, ISKWRK(I = ISKWRK(IiYION)) starts an array of species
C          indices for the ion in a yield-modified reaction;
C          ISKWRK(I + N - 1) is the species index of the ion for the
C          Nth yield-modified reaction.
C  IiKYLD, ISKWRK(I = ISKWRK(IiKYLD)) starts a matrix of yield-
C          modification flags for species in yield-modified reactions;
C          ISKWRK(I + (N-1)*MAXSPR + L - 1)
C          = 1, the Lth species of the Nth yield-modify reaction is to
C               be modified,
C          = 0, the species is not modified.
C
C  ISKWRK POINTERS TO RSKWRK REAL VARIABLES:
C
C  IrSKMN, RSKWRK(I = ISKWRK(IrSKMN)) is the minimum difference
C          allowed for conservation of mass and site.
C  IrPATM, RSKWRK(I = ISKWRK(IrPATM)) is the pressure of one standard
C          atmosphere (dynes/cm**2).
C  IrRU,   RSKWRK(I = ISKWRK(IrRU)) is the universal gas constant
C          (ergs/mole-K).
C  IrRUC,  RSKWRK(I = ISKWRK(IrRUC)) is the universal gas constant
C          (cal/mole-K).
C
C  ISKWRK POINTERS TO START OF RSKWRK ARRAY WORKSPACE:
C
C  IrSDEN, RSKWRK(I = ISKWRK(IrSDEN)) starts an array of phase
C          densities;
C          RSKWRK(I + N - 1) is the density of the Nth phase.
C  IrKTMP, RSKWRK(I = ISKWRK(IrKTMP)) starts a matrix of the
C          dividing temperatures in the thermodynamic fits for
C          the species;
C          RSKWRK(I + (K-1)*MAXTP + N - 1) is the Nth temperature
C          dividing the ranges of coefficients for species K.
C  IrKTHM, RSKWRK(I = ISKWRK(IrKTHM)) starts a three-dimensional
C          array of coefficients for the fits to thermodynamic
C          properties for the species;
C          RSKWRK(I + (L-1)*NCP2 + (K-1)*NCP2T + N - 1) is the Nth
C          polynomial coefficient A(N,L,K) for species K, in the Lth
C          temperature range.
C  IrKDEN, RSKWRK(I = ISKWRK(IrKDEN)) starts an array of species
C          densities;
C          RSKWRK(I + K - 1) is the density of the species K
C          (gm/cm**3 for gas or bulk species, gm/cm**2 for surface
C          species).
C  IrAWT,  RSKWRK(I = ISKWRK(IrAWT)) starts an array of atomic weights;
C          RSKWRK(I + M - 1) is the atomic weight of element M.
C  IrKWT,  RSKWRK(I = ISKWRK(IrKWT)) starts an array of molecular
C          weights;
C          RSKWRK(I + K - 1) is the molecular weight of species K.
C  IrPAR,  RSKWRK(I = ISKWRK(IrPAR)) starts a matrix of Arrhenius
C          parameters for the surface reactions;
C          RSKWRK(I + (N-1)*(NSPAR+1) + L - 1) is, if
C          L=1, the pre-exponential factor (mole-cm-sec-K)
C          L=2, the temperature exponent
C          L=3, the activation energy (K)
C          L=4 is used as a scalar, in sensitivity analysis,
C          for the Nth surface reaction.
C  IrCOV,  RSKWRK(I = ISKWRK(IrKCOV)) starts a matrix of coverage
C          parameters for the coverage surface reactions;
C          RSKWRK(I + (N-1)*NSCOV + L - 1) is the Lth coverage
C          parameter for the Nth coverage reaction.
C  IrRPAR, RSKWRK(I = ISKWRK(IrRPAR)) starts a matrix of reverse
C          Arrhenius parameters for surface reactions which give them
C          explicitly;
C          RSKWRK(I + (N-1)*NSPAR + N - 1), for N=1,3, the
C          reverse parameters for the Nth reverse-parameter reaction,
C          and for N=4, a scaling factor.
C  IrEQ,   RSKWRK(I = ISKWRK(IrEQ)) starts an array of scalars for
C          the surface reaction equilibrium constants;
C          RSKWRK(I + N - 1) is the scalar for the Nth surface reaction.
C  IrRNU,  RSKWRK(I = ISKWRK(IrRNU)) starts a matrix of stoichiometric
C          coefficients for the surface reactions with real
C          coefficients;
C          RSKWRK (I + (N-1)*MAXSPR + L - 1) is the stoichiometric
C          coefficient for the Lth species in the Nth real-
C          stoichiometry reaction.
C  IrNCF,  RSKWRK(I = ISKWRK(IrNCF)) starts a matrix of net site-changes
C          due to the surface reactions;
C          RSKWRK(I + (N-1)*NNPHAS + L - 1) is the net change in sites
C          for phase L due to the Nth surface reaction.
C  IrKORD, RSKWRK(I = ISKWRK(IrKORD)) starts a matrix of species
C          orders for the surface reactions with species change-orders;
C          RSKWRK(I + (N-1)*MAXORD + L - 1) is the order for the Lth
C          change-order species in the Nth change-order reaction.
C  IrKFT,  RSKWRK(I = ISKWRK(IrKFT)) starts an array of the temperature-
C          dependent portion of forward reaction rates for the surface
C          reactions;
C          RSKWRK(I + N - 1) is the forward temperature-dependent rate
C          for the Nth surface reaction.
C  IrKRT,  RSKWRK(I = ISKWRK(IrKRT)) starts an array of the temperature-
C          dependent portion of reverse reaction rates for the surface
C          reactions;
C          RSKWRK(I + N - 1) is the reverse temperature-dependent rate
C          for the Nth surface reaction.
C  IrKT1,  RSKWRK(I = ISKWRK(IrKT1)) starts species scratch space.
C  IrKT2,  RSKWRK(I = ISKWRK(IrKT2)) starts species scratch space.
C  IrPT1,  RSKWRK(I = ISKWRK(IrPT1)) starts phase scratch space.
C  IrPt2,  RSKWRK(I = ISKWRK(IrPT2)) starts phase scratch space.
C  IrIT1,  RSKWRK(I = ISKWRK(IrIT1)) starts reaction scratch space.
C  IrIT2,  RSKWRK(I = ISKWRK(IrIT2)) starts reaction scratch space.
C  IrIT3,  RSKWRK(I = ISKWRK(IrIT3)) starts reaction scratch space.
C  IrPEDP, RSKWRK(I = ISKWRK(IrPEDP)) starts a matrix of parameters
C          for surface reactions with ion-energy dependence;
C          RSKWRK(I + (N-1)*NEDPAR + L - 1) is the Lth parameter for
C          the Nth ion-energy-dependent reaction.
C  IrENGI, RSKWRK(I = ISKWRK(IrENGI)) starts an array of ion energies
C          for the gas-phase species only;
C          RSKWRK(I + K - 1) is the ion energy of species K.
C  IrPYLD, RSKWRK(I = ISKWRK(IrPYLD)) starts a matrix of yield
C          parameters for the surface reactions with yield modification;
C          RSKWRK(I + (N-1)*NYPAR + L - 1) is the Lth parameter for
C          yield-modification reaction N.
C  IrYNCF, RSKWRK(I = ISKWRK(IrYNCF)) starts a matrix of net site
C          changes due to yield changes in the yield-modified reactions;
C          RSKWRK(I + (N-1)*NNPHAS + L - 1) is the net change in sites
C          for phase L due to the Nth yield-modify reaction.
C
C  ISKWRK POINTERS TO CHARACTER WORKSPACE (for one surface):
C
C  IcENAM, CSKWRK(I = ISKWRK(IcENAM)) starts an array of element names;
C          CSKWRK(I + M - 1) is the name of element M.
C  IcKNAM, CSKWRK(I = ISKWRK(IcKNAM)) starts an array of species names;
C          CSKWRK(I + K - 1) is the name of species K.
C  IcPNAM, CSKWRK(I = ISKWRK(IcPNAM)) starts an array of phase names;
C          CSKWRK(I + N - 1) is the name of phase N.
C
C  STORING DATA INTO THE ARRAYS is usually accomplished by a
C  CALL SKINIT, which reads a linkfile generated by the surface
C  mechanism interpreter;
C  the linkfile consists of the following records:
C
C  Linkfile information:
C  1. FILVER, character*16, the linkfile format version
C  2. PRVERS,  character*16, the interpreter program
C  3. PREC,    character*16, the machine precision of the linkfile data
C              (SINGLE, DOUBLE)
C  4. KERR,    logical to indicate whether or not an error was found by
C              the interpreter program
C
C  Parameters and constants:
C  5. LENISK, LENRSK, LENCSK, minimum lengths required to store linkfile
C     data into the integer, real, and character workspace arrays
C  6. MAXSPR, MAXTP, NCP, NSPAR, NSCOV,
C     NEDPAR, NYPAR, MAXORD
C  7. NELEM, NKKGAS, NKKSUR, NKKBLK, NKKTOT,
C     NPHASE, NFSUR, NLSUR, NNSUR, NFBLK,
C     NLBLK, NNBLK, NIISUR, NIICOV, NIIREV,
C     NIISTK, NIICON, NIIBHM, NIIRNU, NIIORD,
C     NIIEDP, NIIYLD, NKKION, KEL, MORE
C  8. SKMIN
C
C  Surface data:
C  9. CSKWRK(ISKWRK(IcMNAM))                            material name
C
C  Element data:
C 10. ( CSKWRK(ISKWRK(IcENAM) + M - 1), M = 1, NELEM)   names
C 11. ( RSKWRK(ISKWRK(IrAWT)  + M - 1), M = 1, NELEM)   weight
C
C  Species data:
C 12. ( CSKWRK(ISKWRK(IcKNAM) + K - 1), K = 1, NKKTOT)  names
C 13. ( RSKWRK(ISKWRK(IrKWT)  + K - 1), K = 1, NKKTOT)  weight
C 14. ( (ISKWRK(ISKWRK(IiKCMP)+(K - 1)*NELEM + M - 1),
C        M = 1, NELEM), K = 1, NKKTOT)                  composition
C 15. ( ISKWRK(ISKWRK(IiKCHG) + K - 1), K = 1, NKKTOT)  charge
C 16. ( ISKWRK(ISKWRK(IiKNT)  + K - 1), K = 1, NKKTOT)  #fit temp's
C 17. ( ISKWRK(ISKWRK(IiKPHS) + K - 1), K = 1, NKKTOT)  phase
C 18. ( ISKWRK(ISKWRK(IiNSCV) + K - 1), K = 1, NKKTOT)  coverage
C 19. ( RSKWRK(ISKWRK(IrKDEN) + K - 1), K = 1, NKKTOT)  density
C 20. ( (RSKWRK(ISKWRK(IrKTMP)+(K - 1)*MAXTP + L - 1),  fit temperatures
C        L = 1, MAXTP), K = 1, NKKTOT)
C 21. ((( RSKWRK(ISKWRK(IrKTHM)+(L-1)*NCP2+(K-1)*NCP2T + N-1),
C         N = 1, NCP2), L = 1, NTR), K = 1, NKKTOT)     thermo coeff'nts
C
C  Ion data (if NKKION > 0):
C 22. NKKION
C 23. ( ISKWRK(ISKWRK(IiIONS) + K - 1), K = 1, NKKION)  species indices
C
C  Phase data:
C 24. ( CSKWRK(ISKWRK(IcPNAM) + N - 1), N = 1, NPHASE)  names
C 25. ( ISKWRK(ISKWRK(IiPKST) + N - 1), N = 1, NPHASE)  starting species
C 26. ( ISKWRK(ISKWRK(IiPKND) + N - 1), N = 1, NPHASE)  ending species
C 27. ( ISKWRK(ISKWRK(IiPTOT) + N - 1), N = 1, NPHASE)  species count
C 28. ( RSKWRK(ISKWRK(IrSDEN) + N - 1), N = 1, NPHASE)  density
C
C  Reaction data (if NIISUR > 0):
C 29. ( ISKWRK(ISKWRK(IiNRPP) + I - 1), I = 1, NIISUR)  species count
C 30. ( ISKWRK(ISKWRK(IiNREA) + I - 1), I = 1, NIISUR)  reactant count
C 31. ( (ISKWRK(ISKWRK(IiNU) + (I-1)*MAXSPR + N - 1),   stoichiometry
C        ISKWRK(ISKWRK(IiNUNK)+(I-1)*MAXSPR + N - 1),   species indices
C       N = 1, MAXSPR), I = 1, NIISUR)
C 32. ( ISKWRK(ISKWRK(IiNSUM) + I - 1), I = 1, NIISUR)  stoich. sum
C 33. ( (RSKWRK(ISKWRK(IrPAR)+(I-1)*(NSPAR+1)+N-1),
C       N = 1, NSPAR), I = 1, NIISUR)                   Arrh. coeff'nts
C 34. ( RSKWRK(ISKWRK(IrEQ) + I - 1), I = 1, NIISUR)    equil. factor
C 35. ( (RSKWRK(ISKWRK(IrNCF) + (I-1)*NPHASE+N-1),      phase balance
C       N = 1, NPHASE), I = 1, NIISUR)
C
C  Coverage reaction data (if NIICOV > 0):
C 36. NIICOV, NSCOV
C 37. ( ISKWRK(ISKWRK(IiICOV) + N - 1), N = 1, NIICOV)  reaction indices
C 38. ( ISKWRK(ISKWRK(IiKCOV) + N - 1), N = 1, NIICOV)  species indices
C 39. ( (RSKWRK(ISKWRK(IrKCOV)+(N-1)*NSCOV+L-1),        parameters
C       L = 1, NSCOV), N = 1, NIICOV)
C
C  Reverse reaction data (if NIIREV > 0):
C 40. NIIREV
C 41. ( ISKWRK(ISKWRK(IiIREV) + N - 1), N = 1, NIIREV)  reaction indices
C 42. ( (RSKWRK(ISKWRK(IrPAR)+(N-1)*(NSPAR_1)+L-1),     rev. parameters
C       L = 1, NSPAR), N = 1, NIIREV)
C
C  Sticking reaction data (if NIISTK > 0):
C 43. NIISTK
C 44. ( ISKWRK(ISKWRK(IiISTK) + N - 1), N = 1, NIISTK)  reaction indices
C 45. ( ISKWRK(ISKWRK(IiMSTK) + N - 1), N = 1, NIISTK)  Motz-wise flag
C
C  Bohm reaction data (if NIIBHM > 0):
C 46. NIIBHM
C 47. ( ISKWRK(ISKWRK(IiBHM) + N - 1),  N = 1, NIIBHM)  reaction indices
C 48. ( ISKWRK(ISKWRK(IiKBHM)+ N - 1),  N = 1, NIIBHM)
C
C Real stoichiometry data (if NIIRNU > 0):
C 49. NIIRNU
C 50. ( ISKWRK(ISKWRK(IiRNU) + N - 1), N = 1, NIIRNU)   reaction indices
C 51. ( (RSKWRK(ISKWRK(IrRNU) + (N-1)*MAXSPR + L - 1),  real coeff'nts
C       L = 1, MAXSPR), N = 1, NIIRNU)
C
C Change-order data (if NIIORD > 0):
C 52. NIIORD, MAXORD
C 53. ( ISKWRK(ISKWRK(IiIORD)+ N - 1), N = 1, NIIORD)  reaction indices
C 54. ( (ISKWRK(ISKWRK(IiKORD)+(N-1)*MAXORD + L - 1),  species indices
C       L = 1, MAXORD), N = 1, NIIORD)
C 55. ( (RSKWRK(ISKWRK(IrKORD)+(N-1)*MAXORD + L - 1),  order values
C       L = 1, MAXORD), N = 1, NIIORD)
C
C Temperature-dependent reaction data (if NIIEDP > 0):
C 56. NIIEDP, NEDPAR
C 57. ( ISKWRK(ISKWRK(IiIEDP)+ N - 1), N = 1, NIIEDP)  reaction indices
C 58. ( ISKWRK(ISKWRK(IiKEDP)+ N - 1), N = 1, NIIEDP)  species indices
C 59. ( (RSKWRK(ISKWRK(IrPEDP)+(L-1)*NEDPAR + L - 1),  parameters
C       L = 1, NEDPAR), I = 1, NIIEDP)
C
C Yield-modify data (if NIIYLD > 0):
C 60. NIIYLD, NYPAR
C 61. ( ISKWRK(ISKWRK(IiIYLD)+ N - 1), N = 1, NIIYLD)  reaction indices
C 62. ( ISKWRK(ISKWRK(IiYION)+ N - 1), N = 1, NIIYLD)  ion indices
C 63. ( (ISKWRK(ISKWRK(IiKYLD)+(N-1)*MAXSPR + L - 1),  yield flags
C       L = 1, MAXSPR), N = 1, NIIYLD)
C 64. ( (RSKWRK(ISKWRK(IrPYLD)+(N-1)*NYPAR + L - 1),   yield parameters
C       L = 1, NYPAR), N = 1, NIIYLD)
C 65. ( (RSKWRY(ISKWRY(IrYNCF)+(N-1)*NPHASE + L - 1),  phase balance
C       L = 1, NPHASE), N = 1, NIIYLD)
C
C  END PROLOGUE
C
C     End of SUBROUTINE SKABS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKABE  (ISKWRK, RSKWRK, RA, RB, RE, ISTFL)
C
C  START PROLOGUE
C
C  SUBROUTINE SKABE  (ISKWRK, RSKWRK, RA, RB, RE, ISTFL)
C  Returns the Arrhenius coefficients or the sticking coefficients
C  of the surface reactions, and integer flags to indicate the type
C  of the coefficients.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real workspace array;    dimension at least LENRSK.
C
C  OUTPUT
C  RA(*)     - Real array, pre-exponential constants for reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, mole-cm-sec-K
C  RB(*)     - Real array, temperature dependence exponents for
C              reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C  RE(*)     - Real array, activation energies for reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, K
C  ISTFL(*)  - Integer array, sticking reaction information;
C              dimension at least IISUR, the total surface reaction
C              count.
C              =1, a reaction uses sticking coefficients.
C              =0, a rection does not.
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
      INCLUDE 'skstrt.h'
      DIMENSION RA(*), RB(*), RE(*), ISKWRK(*), RSKWRK(*), ISTFL(*)
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN

      I_PAR = ISKWRK(IrPAR)
      DO 100 I = 1, ISKWRK(IiNIIS)
         RA(I) = RSKWRK(I_PAR)
         RB(I) = RSKWRK(I_PAR+1)
         RE(I) = RSKWRK(I_PAR+2)
         I_PAR = I_PAR + (NSPAR+1)
         ISTFL(I) = 0
  100 CONTINUE
C
      IF (ISKWRK(IiNSTK) .LE. 0) RETURN
C
      I_ISTK = ISKWRK(IiISTK) - 1
      DO 150 N = 1, ISKWRK(IiNSTK)
         I     = ISKWRK(I_ISTK + N)
         ISTFL(I) = 1
  150 CONTINUE
C
C     end of SUBROUTINE SKABE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKAML  (T, ISKWRK, RSKWRK, AML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKAML  (T, ISKWRK, RSKWRK, AML)
C  Returns the standard state Helmholtz free energies in molar units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  AML(*)    - Real array, standard state Helmholtz free energies
C              for species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), AML(*), T(*)
C
      I_SML = ISKWRK(IrKT1)
      I_HML = ISKWRK(IrKT2)
      CALL SKSML (T, ISKWRK, RSKWRK, RSKWRK(I_SML))
      CALL SKHML (T, ISKWRK, RSKWRK, RSKWRK(I_HML))
C
      RU = RSKWRK(ISKWRK(IrRU))
      I_KTFL = ISKWRK(IiKTFL)
C
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 100 K = 0, NKM1
         IF (K .LT. NKKGAS) THEN
            TK = T(ISKWRK(I_KTFL+K))
         ELSE
            TK = T(1)
         ENDIF
         AML(K+1) = RSKWRK(I_HML+K) - TK * (RU + RSKWRK(I_SML+K))
100   CONTINUE
C
C     end of SUBROUTINE SKAML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKAMS  (T, ISKWRK, RSKWRK, AMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKAMS  (T, ISKWRK, RSKWRK, AMS)
C  Returns an the standard state Helmholtz free energies in mas units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  AMS(*)    - Real array, standard state Helmholtz free energies
C              for species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), AMS(*), T(*)
C
      I_SMS = ISKWRK(IrKT1)
      I_HMS = ISKWRK(IrKT2)
      CALL SKSMS (T, ISKWRK, RSKWRK, RSKWRK(I_SMS))
      CALL SKHMS (T, ISKWRK, RSKWRK, RSKWRK(I_HMS))
C
      RU = RSKWRK(ISKWRK(IrRU))
      I_KWT = ISKWRK(IrKWT)
      I_KTFL = ISKWRK(IiKTFL)
C
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 100 K = 0, NKM1
         IF (K .LT. NKKGAS) THEN
            TK = T(ISKWRK(I_KTFL+K))
         ELSE
            TK = T(1)
         ENDIF
         AMS(K+1) = RSKWRK(I_HMS + K) - TK *
     1            (RU / RSKWRK(I_KWT + K) + RSKWRK(I_SMS + K))
100   CONTINUE
C
C     end of SUBROUTINE SKAMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, CZ)
C
C  START PROLOGUE
C
C  SUBROUTINE SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, CZ)
C  Returns the concentrations of the species, given the pressure,
C  temperature and activities.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  CZ(*)     - Real array, gas-phase and surface species concentrations,
C              and bulk species activities;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS gas-phase concentrations are moles/cm**3,
C              the next KKSURF site concentrations are moles/cm**2, and
C              the final KKBULK entries are bulk species activities.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), CZ(*), T(*)
C
C       Compute the molar concentrations from the mole fractions.
C       Temperature is computed from an average temperature
C       formulation!??
C       Sum of X_k term in the numerator was introduced so that no
C       extra Jacobian entries would be created by the average
C       temperature term.
C
      SUMXT = 0.0
      PRUT  = 0.0
      RU = RSKWRK(ISKWRK(IrRU))
      I_KTFL = ISKWRK(IiKTFL) - 1
      DO 50 K = 1, NKKGAS
         SUMXT = SUMXT + ACT(K)*T(ISKWRK(I_KTFL + K))
         PRUT  = PRUT  + ACT(K)
   50 CONTINUE
      PRUT = PRUT * P / (RU * SUMXT)
      DO 100 K = 1, NKKGAS
         CZ(K) = ACT(K) * PRUT
  100 CONTINUE
C
      IF (ISKWRK(IiKSUR) .GT. 0) THEN
C        Surface concentrations are calculated from the site
C        site density (which is a variable)
C        divided by the number of sites per species.
C
         KCOV = ISKWRK(IiNSCV) - 1
         IKST = ISKWRK(IiPKST) - 1
         IKND = ISKWRK(IiPKND) - 1
         DO 210 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            KFIRST = ISKWRK(IKST + N)
            KLAST  = ISKWRK(IKND + N)
            DO 205 K = KFIRST, KLAST
              CZ(K) = ACT(K) * SDEN(N) / ISKWRK(KCOV + K)
  205       CONTINUE
  210    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiKBLK) .GT. 0) THEN
C        Bulk molefractions; the "concentration vector"
C        for these unknowns is just the bulk mole fraction itself
C
         KFIRST = ISKWRK(ISKWRK(IiPKST) + ISKWRK(IiFBLK) - 1)
         KLAST  = ISKWRK(ISKWRK(IiPKND) + ISKWRK(IiLBLK) - 1)
         DO 220 K = KFIRST, KLAST
            CZ(K) = ACT(K)
  220    CONTINUE
      ENDIF
C
C     end of SUBROUTINE SKATCX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKATHM (MDIM, NDIM1, NDIM2, ISKWRK, RSKWRK, NT, TMP,
     1                   A)
C
C  START PROLOGUE
C
C  SUBROUTINE SKATHM (MDIM, NDIM1, NDIM2, ISKWRK, RSKWRK, NT, TMP,
C                     A)
C  Returns the polynomial coefficients of the fits for
C  thermodynamic properties of all of the species.
C
C  INPUT
C  MDIM      - Integer scalar, first dimension of an array of
C              temperatures used in thermodynamic fits for species;
C              MDIM must be at least MAXTP, the maximum number of
C              temperatures used to fit the thermodynamics.
C  NDIM1     - Integer scalar, first dimension of A, the three-
C              dimensional array of thermodynamic fit coefficients;
C              NDIM1 must be at least NPCP2, the total number of
C              coefficients for one temperature range.
C  NDIM2     - Integer scalar, second dimension of A; NDIM2 must be
C              at least MAXTP-1, the total number of temperature ranges.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  Where NT(K) is the number of temperatures used in fitting the
C  thermodynamic properties of species K, TMP(N) is the Nth
C  temperature, NT(K)-1 is the number of temperature ranges for
C  which the polynomial coefficients are valid, then
C  A (L, N, K) is the Lth polynomial coefficient, for the Nth
C  temperature range, and the Kth species; i.e.,
C
C          | <   N = 1  >. <N=2> .                    .<  N = NT - 1>
C    P  E  |    .        .       .                    .        .
C    O  X  |    .        .       .                    .        .
C    L  P  |    .        .       .                    .        .
C    Y  R  |    .        .       .                    .        .
C    N  E  |    .        .       .                    .        .
C    O  S  |    .        .       .                    .        .
C    M  S  |    .        .       .                    .        .
C    I  I  |    .        .       .                    .        .
C    A  O  |    .        .       .                    .        .
C    L  N  |____.________._______.____________________.________.______
C             TMP(1)  TMP(2)  TMP(3) . .  .  .  . TMP(NT-1)  TMP(NT)
C
C  NT(*)     - Integer array, total number of temperatures used in
C              fitting coefficients of thermodynamic properties for
C              the species;
C              dimension at least KKTOT, the total species count.
C  TMP(*,*)  - Real matrix, temperatures for dividing the
C              thermodynamic fits for species; dimension at least
C              MAXTP for the first, and at least KKTOT for the second,
C              the total species count.
C                 cgs units, K
C  A(*,*,*)  - Real three-dimensioned array of fit coefficients to the
C              thermodynamic data for species;
C              dimension exactly NPCP2 for the first, exactly MAXTP-1
C              for the second, and at least KKTOT for the third, the
C              total species count.
C              The indicies in  A(N,L,K) mean-
C              N = 1,NN represent polynomial coefficients in CP/R
C                CP/R(K)=A(1,L,K) + A(2,L,K)*T + A(3,L,K)*T**2 + ...
C              N = NN+1 is for the formation enthalpies, i.e.,
C                HO/R = A(NN+1,L,K)
C              N = NN+2 is for the formation entropies, i.e.,
C                SO/R = A(NN+2,L,K)
C              L = 1 is for temperature <= TMP(2,K)
C              L = 2 is for TMP(2,K) < temperature <= TMP(3)
C                :
C              L = (NTMP-1) is for TMP(NTMP-1) <= temperature;
C              K  is  the  species index
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION NT(*),TMP(MDIM,*),A(NDIM1,NDIM2,*),ISKWRK(*),
     1          RSKWRK(*)
C
      DO 100 K = 1, ISKWRK(IiKTOT)
         NT(K) = ISKWRK(ISKWRK(IiKNT) + K - 1)
  100 CONTINUE
C
      IF (MDIM .LT. MAXTP) RETURN
      I_KTMP = ISKWRK(IrKTMP)
      DO 140 L = 1, MAXTP
         DO 130 K = 1, ISKWRK(IiKTOT)
            TMP(L,K) = RSKWRK(I_KTMP)
            I_KTMP = I_KTMP + 1
  130    CONTINUE
  140 CONTINUE
C
      I_NA1 = ISKWRK(IrKTHM)
      DO 170 K = 1, ISKWRK(IiKTOT)
         DO 160 L = 1, MAXTP-1
             DO 150 M = 1, NCP2
                A(M, L, K) = RSKWRK(I_NA1)
                I_NA1 = I_NA1 + 1
  150        CONTINUE
  160    CONTINUE
  170 CONTINUE
C
C     end of SUBROUTINE SKATHM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCHRG (ISKWRK, RSKWRK, KCHARG)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCHRG (ISKWRK, RSKWRK, KCHARG)
C  Returns an array containing electronic charges of the species.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  KCHARG(*) - Integer array, electronic charges of the species;
C              dimension at least KKTOT, the total species count.
C              KCHARG(K)=-2 indicates that the species K has two excess
C              electrons.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), KCHARG(*)
C
      I_KCHG = ISKWRK(IiKCHG) - 1
      DO 100 K = 1, ISKWRK(IiKTOT)
         KCHARG(K) = ISKWRK(I_KCHG + K)
  100 CONTINUE
C
C     end of SUBROUTINE SKCHRG
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCOMP (ISTR, IRAY, NN, IND, NT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCOMP (ISTR, IRAY, NN, IND, NT)
C  Search for the occurrence of character string ISTR, in the NN
C  character strings of array IRAY;
C  IND is the first location in IRAY of ISTR if found, or 0 if not
C  found, and NT is the total number of times it occurs.
C
C     Consider the following example,
C        IRAY = {"BOOK","BLUE","BEAR","BOOK"}
C        NN=4.
C
C     If ISTR="BLUE" then IND=2 and NT=1;
C     if ISTR="RED"  then IND=0 and NT=0; and
C     if ISTR="BOOK",then IND=1 and NT=2.
C
C  INPUT
C  ISTR      - Character string.
C  IRAY(*)   - Character string array.
C  NN        - Integer scalar, length of IRAY(*).
C
C  OUTPUT
C  IND       - Integer scalar, location in IRAY of the character string
C              ISTR, or 0 if ISTR does not appear in IRAY.
C  NT        - Integer scalar, total number of times ISTR occurs
C              in IRAY.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) ISTR,IRAY(NN)
C
      CALL CKNCMP (ISTR, IRAY, NN, IND, NT)
C
C     end of SUBROUTINE SKCOMP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCONT (KSPEC, ROP, ISKWRK, RSKWRK, CIK)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCONT (KSPEC, ROP, ISKWRK, RSKWRK, CIK)
C  Returns the contributions of the surface reactions to the molar
C  production rate of species KSPEC.
C
C  INPUT
C  KSPEC     - Integer scalar, species index.
C  ROP(*)    - Real array, rates of progress for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec)
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  CIK(*)    - Real array, contributions of the surface reactions to the
C              production rate of species KSPEC;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, mole/(cm**2*sec)
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
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), ROP(*), CIK(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      NIISUR = ISKWRK(IiNIIS)
      IF (NIISUR .LE. 0) RETURN
C
      I_NK = ISKWRK(IiNUNK)
      I_NU = ISKWRK(IiNU)
C
      DO 200 I = 1, NIISUR
         CIK(I) = 0.0
C
         DO 190 N = 0, MAXSPR - 1
            IF (ISKWRK(I_NK + N) .EQ. KSPEC)
     1         CIK(I) = CIK(I) + ROP(I)*ISKWRK(I_NU + N)
  190       CONTINUE
            I_NU = I_NU + MAXSPR
            I_NK = I_NK + MAXSPR
200   CONTINUE
C
      IF (ISKWRK(IiNRNU) .GT. 0) THEN
C        RNC's above were 0
         I_NU  = ISKWRK(IrRNU)
         I_IRNU = ISKWRK(IiIRNU)
         DO 300 N = 0, ISKWRK(IiNRNU) - 1
            I = ISKWRK(I_IRNU + N)
            I_NK = ISKWRK(IiNUNK) + MAXSPR*(I-1)
            DO 290 L = 0, MAXSPR - 1
               IF (ISKWRK(I_NK + L) .EQ. KSPEC)
     1            CIK(I) = CIK(I) + ROP(I)*RSKWRK(I_NU + L)
  290       CONTINUE
            I_NU = I_NU + MAXSPR
  300    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiNYLD) .LE. 0) RETURN
C
      I_IYLD = ISKWRK(IiIYLD)
      DO 400 L = 0, ISKWRK(IiNYLD) - 1
         I     = ISKWRK(I_IYLD + L)
         ISRNU = CKLKUP (I, ISKWRK(ISKWRK(IiIRNU)), ISKWRK(IiNRNU))
         IF (ISRNU .GT. 0) THEN
            I_NU = ISKWRK(IrRNU) + (ISRNU-1)*MAXSPR
         ELSE
            I_NU = ISKWRK(IiNU) + (I-1)*MAXSPR
         ENDIF
         I_NK = ISKWRK(IiNUNK) + MAXSPR*(I-1)
C
         DO 390 N = 0, MAXSPR - 1
            K = ISKWRK(I_NK + N)
C           is this an active species
            IF (K .EQ. KSPEC) THEN
C              does this species require yield-modify?
               KYLD = ISKWRK(ISKWRK(IiKYLD) + L*MAXSPR + N)
C
               IF (KYLD .GT. 0) THEN
C                 does this species have integer or real coefficient?
                  IF (ISRNU .GT. 0) THEN
                     COEF = RSKWRK(I_NU + N)
                  ELSE
                     COEF = ISKWRK(I_NU + N)
                  ENDIF
C
C                 actual species index of the ion
                  KION  = ISKWRK( ISKWRK(IiYION) + L)
C
C                 location of first yield parameter
                  IPYLD = ISKWRK(IrPYLD) + L * NYPAR
                  EI    = RSKWRK( ISKWRK(IrENGI) + KION - 1 )
                  ETH   = RSKWRK( IPYLD + 1)
                  IF (EI .GE. ETH) THEN
                     ASCAL = RSKWRK( IPYLD )
                     A     = RSKWRK( IPYLD + 2)
                     B     = RSKWRK( IPYLD + 3)
                     SCALE = ASCAL * (EI**A - ETH**A) ** B
                  ELSE
                     SCALE = 0.0
                  ENDIF
C                 need to delete previous ROP*COEF term,
C                 and modify before adding back in
                  CIK(I) = CIK(I) + COEF*(SCALE - 1.0) * ROP(I)
               ENDIF
            ENDIF
  390    CONTINUE
  400 CONTINUE
C
C     end of SUBROUTINE SKCONT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCOV  (ISKWRK, KOCC)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCOV  (ISKWRK, KOCC)
C  Returns an array of site occupancy numbers for the species.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C
C  OUTPUT
C  KOCC(*)   - Integer array, site occupancy numbers for the species;
C              dimension at least KKTOT, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), KOCC(*)
C
      I_KCOV = ISKWRK(IiNSCV) - 1
      DO 50 K = 1, ISKWRK(IiKTOT)
         KOCC(K) = ISKWRK(I_KCOV + K)
   50 CONTINUE
C
C     end of SUBROUTINE SKCOV
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCPML (T, ISKWRK, RSKWRK, CPML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCPML (T, ISKWRK, RSKWRK, CPML)
C  Returns an array of the specific heats at constant pressure
C  in molar units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  CPML(*)   - Real array, specific heats at constant pressure for the
C              species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, ergs/(mole*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), CPML(*), T(*)
      SAVE TN1, TN2, TN3, TN4
      DATA TN1/0.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
      ENDIF
      I_KTFL = ISKWRK(IiKTFL)
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 250 K = 0, NKM1
         IF (K.LT.NKKGAS .AND. ISKWRK(I_KTFL+K).GT.1) THEN
            TK1 = T(ISKWRK(I_KTFL+K))
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
         ELSE
            TK1 = TN1
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ISKWRK(ISKWRK(IiKNT)+K) - 1
C        location of FIRST set of thermodynamic coefficients
         NA1 = ISKWRK(IrKTHM) + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = ISKWRK(IrKTMP) + K*MAXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RSKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         CPML(K+1) = RSKWRK(ISKWRK(IrRU)) * (RSKWRK(NA1)
     1             + RSKWRK(NA1+1) * TK1 + RSKWRK(NA1+2) * TK2
     2             + RSKWRK(NA1+3) * TK3 + RSKWRK(NA1+4) * TK4)
C
250   CONTINUE
C
C     end of SUBROUTINE SKCPML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCPMS (T, ISKWRK, RSKWRK, CPMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCPMS (T, ISKWRK, RSKWRK, CPMS)
C  Returns an array of the specific heats at constant pressure
C  in mass units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C   OUTPUT
C   CPMS(*)  - Real array, specific heats at constant pressure for the
C              species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, ergs/(gm*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), CPMS(*), T(*)
      SAVE TN1, TN2, TN3, TN4
      DATA TN1/0.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
      ENDIF
      RU = RSKWRK(ISKWRK(IrRU))
      I_KTFL = ISKWRK(IiKTFL)
      I_KWT  = ISKWRK(IrKWT)
C
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 250 K = 0, NKM1
         IF (K.LT.NKKGAS .AND. ISKWRK(I_KTFL+K).GT.1) THEN
            TK1 = T(ISKWRK(I_KTFL+K))
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
         ELSE
            TK1 = TN1
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ISKWRK(ISKWRK(IiKNT) + K) - 1
C        location of FIRST set of thermodynamic coefficients
         NA1 = ISKWRK(IrKTHM) + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = ISKWRK(IrKTMP) + K*MAXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RSKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         CPMS(K+1) = RU * (RSKWRK(NA1) + RSKWRK(NA1+1) * TK1
     1             + RSKWRK(NA1+2) * TK2 + RSKWRK(NA1+3) * TK3
     2             + RSKWRK(NA1+4) * TK4) / RSKWRK(I_KWT + K)
250   CONTINUE
C
C     end of SUBROUTINE SKCPMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCPOR (T, ISKWRK, RSKWRK, CPOR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCPOR (T, ISKWRK, RSKWRK, CPOR)
C  Returns an array of the nondimensional specific heats at constant
C  pressure.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  CPOR(*)   - Real array, nondimensional specific heats at constant
C              pressure for the species;
C              dimension at least KKTOT, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), CPOR(*), T(*)
      SAVE TN1, TN2, TN3, TN4
      DATA TN1/0.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
      ENDIF
      I_KTFL = ISKWRK(IiKTFL)
C
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 250 K = 0, NKM1
         IF (K.LT.NKKGAS .AND. ISKWRK(I_KTFL+K).GT.1) THEN
            TK1 = T(ISKWRK(I_KTFL+K))
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
         ELSE
            TK1 = TN1
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ISKWRK(ISKWRK(IiKNT) + K) - 1
C        location of FIRST set of thermodynamic coefficients
         NA1 = ISKWRK(IrKTHM) + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = ISKWRK(IrKTMP) + K*MAXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RSKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         CPOR(K+1) = RSKWRK(NA1)
     1             + RSKWRK(NA1+1) * TK1 + RSKWRK(NA1+2) * TK2
     2             + RSKWRK(NA1+3) * TK3 + RSKWRK(NA1+4) * TK4
250   CONTINUE
C
C     end of SUBROUTINE SKCPOR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCZTA (T, CZ, SDEN, ISKWRK, RSKWRK, ACT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCZTA (T, CZ, SDEN, ISKWRK, RSKWRK, ACT)
C  Returns the activities of the species, given the pressure,
C  temperature and concentrations.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  CZ(*)     - Real array, gas-phase and surface species concentrations,
C              and bulk species activities;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS gas-phase concentrations are moles/cm**3,
C              the next KKSURF site concentrations are moles/cm**2, and
C              the final KKBULK entries are bulk species activities.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), CZ(*), T(*)
C
      SUMCZ = 0.0
      DO 50 K = 1, NKKGAS
         SUMCZ = SUMCZ + CZ(K)
   50 CONTINUE
      DO 100 K = 1, NKKGAS
         ACT(K) = CZ(K) / SUMCZ
  100 CONTINUE
C
      IF (ISKWRK(IiKSUR) .GT. 0) THEN
C
C        Surface concentrations are calculated from the site
C        density (which is a variable)
C        divided by the number of sites per species.
         KCOV = ISKWRK(IiNSCV)
         NKF  = ISKWRK(IiPKST) - 1
         NKL  = ISKWRK(IiPKND) - 1
         DO 210 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            KFIRST = ISKWRK(NKF + N)
            KLAST  = ISKWRK(NKL + N)
            DO 205 K = KFIRST, KLAST
              ACT(K) = CZ(K) * ISKWRK(KCOV + K - 1) / SDEN(N)
  205       CONTINUE
  210    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiKBLK) .GT. 0) THEN
C
C       Bulk molefractions - The "concentration vector"
C       for these unknowns is just the bulk mole fraction itself
         KFIRST = ISKWRK(ISKWRK(IiPKST) + ISKWRK(IiFBLK) - 1)
         KLAST  = ISKWRK(ISKWRK(IiPKND) + ISKWRK(IiLBLK) - 1)
         DO 220 K = KFIRST, KLAST
            ACT(K) = CZ(K)
  220    CONTINUE
      ENDIF
C
C     end of SUBROUTINE SKCTZA
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKDEN  (P, T, ACT, SDEN, ISKWRK, RSKWRK, DEN)
C
C  START PROLOGUE
C
C  SUBROUTINE SKDEN  (P, T, ACT, SDEN, ISKWRK, RSKWRK, DEN)
C  Returns a real array of species densities.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  DEN(*)    - Real array, densities for the species;
C              dimension at least KKTOT, the total species count.
C                 cgm units, gm/cm**3 for gas-phase species
C                            bm/cm**2 for surface species
C                            gm/cm**3 for bulk species
C
C              NOTE:  mass densities are not required to be input to
C                     the Interpreter for bulk-phase species.
C                     If they are input, they are returned by this
C                     subroutine.  If not, DEN = -1.0 for the bulk
C                     species
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), DEN(*), T(*)
C
      SUMXT = 0.0
      I_KTFL = ISKWRK(IiKTFL) - 1
      DO 50 K = 1, NKKGAS
         SUMXT = SUMXT + ACT(K)*T(ISKWRK(I_KTFL+K))
   50 CONTINUE
C
      PRUT = P/(RSKWRK(ISKWRK(IrRU))*SUMXT)
C
      I_KWT = ISKWRK(IrKWT) - 1
      DO 120 K = 1, NKKGAS
         DEN(K) = ACT(K) * RSKWRK(I_KWT + K) * PRUT
  120 CONTINUE
C
      IF (ISKWRK(IiKSUR) .GT. 0) THEN
         NKF = ISKWRK(IiPKST) - 1
         NKL = ISKWRK(IiPKND) - 1
         NKV = ISKWRK(IiNSCV) - 1
         DO 130 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            KFIRST = ISKWRK(NKF + N)
            KLAST  = ISKWRK(NKL + N)
            DO 125 K = KFIRST, KLAST
               RCOV = ISKWRK(NKV + K)
               DEN(K) =  ACT(K)*SDEN(N)*RSKWRK(I_KWT + K) / RCOV
  125       CONTINUE
  130    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiKBLK) .GT. 0) THEN
         I_KDEN = ISKWRK(IrKDEN) - 1
         KFIRST = ISKWRK(ISKWRK(IiPKST) + ISKWRK(IiFBLK) - 1)
         KLAST  = ISKWRK(ISKWRK(IiPKND) + ISKWRK(IiLBLK) - 1)
         DO 140 K = KFIRST, KLAST
            DEN(K) = RSKWRK(I_KDEN + K)
  140    CONTINUE
      ENDIF
C
C     end of SUBROUTINE SKDEN
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKDRDA (IR, P, T, ACT, SDEN, ISKWRK, RSKWRK, DKDAI)
C
C  START PROLOGUE
C
C  SUBROUTINE SKDRDA (IR, P, T, ACT, SDEN, ISKWRK, RSKWRK, DKDAI)
C  Returns the partial of the rates of production of the species with
C  respect to the pre-exponential constant of surface reaction IR.
C
C  INPUT
C  IR        - Integer scalar, surface reaction index.
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  DKDAI(*)  - Real array, partials of the partial of production rates
C              of the species with respect to the pre-exponential
C              constant for surface reaction IR;
C              dimension at least KKTOT, the total species count.
C                 cgs units, moles/(cm**2*sec) / (units of A)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), DKDAI(*), T(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      SDTOT = ZERO
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
         IF (ISKWRK(IiNYLD) .GT. 0) CALL SKEQY (ISKWRK, RSKWRK)
      ENDIF
C
C     initialize some local pointers in COMM~ON /SKPTR/
      CALL SKLOC (ISKWRK)
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, RSKWRK(I_CZ))
      CALL SKRROP (ISKWRK, RSKWRK, RSKWRK(I_SMH), T, RSKWRK(I_CZ),
     1             ACT, RSKWRK(I_KWT), ISKWRK(I_NREA), ISKWRK(I_NRPP),
     2             ISKWRK(I_NU), ISKWRK(I_NK), ISKWRK(I_IRNU),
     3             RSKWRK(I_RNU), ISKWRK(I_NUSUM), RSKWRK(I_PAR),
     4             RSKWRK(I_RPAR), ISKWRK(I_IREV), ISKWRK(I_ICOV),
     5             ISKWRK(I_KCOV), RSKWRK(I_CPAR), ISKWRK(I_ISTK),
     6             ISKWRK(I_MSTK), ISKWRK(I_IBHM), ISKWRK(I_KBHM),
     7             ISKWRK(I_IEDP), ISKWRK(I_KEDP), RSKWRK(I_PEDP),
     8             ISKWRK(I_KTFL), SDTOT, ISKWRK(I_IORD),
     9             ISKWRK(I_KORD), RSKWRK(I_ORD), ISKWRK(I_IYLD),
     *             ISKWRK(I_YION), ISKWRK(I_KYLD), RSKWRK(I_PYLD),
     1             RSKWRK(I_RKFT), RSKWRK(I_RKRT), RSKWRK(I_RKF),
     2             RSKWRK(I_RKR), RSKWRK(I_EQKC), RSKWRK(I_EQFAC),
     3             RSKWRK(I_ENRG))
C
C     PROCESS REACTION IR
C
C     reverse Arrhenius parameters specified?
      ISREV = CKLKUP (IR, ISKWRK(I_IREV), ISKWRK(IiNREV))
C     real stoichiometry?
      ISRNU = CKLKUP (IR, ISKWRK(I_IRNU), ISKWRK(IiNRNU))
C     #-modified stoichiometry?
      ISYLD = CKLKUP (IR, ISKWRK(I_IYLD), ISKWRK(IiNYLD))
C
      RKF = RSKWRK(I_RKF + IR - 1)
      RKR = RSKWRK(I_RKR + IR - 1)
      CALL SKRAEX (IR, ISKWRK, RSKWRK, RA)
C
      DO 50 K = 1, ISKWRK(IiKTOT)
         DKDAI(K) = ZERO
   50 CONTINUE
      I_NK = I_NK + (IR-1)*MAXSPR - 1
      IF (ISRNU .GT. 0) THEN
         I_NU = I_RNU + (ISRNU-1)*MAXSPR - 1
      ELSE
         I_NU = I_NU + (IR-1)*MAXSPR - 1
      ENDIF
C
      DO 80 N = 1, MAXSPR
C
C        is this an active species?
         K = ISKWRK(I_NK + N)
         IF (K .LE. 0) GO TO 80
C
         IF (ISRNU .GT. 0) THEN
            COEF = RSKWRK(I_NU + N)
         ELSE
            COEF = ISKWRK(I_NU + N)
         ENDIF
C
C        is this a yield-modify reaction
         IF (ISYLD .GT. 0) THEN
C
C           does the species require modification?
            KYLD = ISKWRK(I_KYLD + (ISYLD-1)*MAXSPR + N - 1)
            IF (KYLD .GT. 0) THEN
C
C              species index for the ion
               KION  = ISKWRK( I_YION + ISYLD - 1 )
C
C              location of first yield parameter
               IPYLD = I_PYLD + NYPAR*(ISYLD-1)
               EI    = RSKWRK( I_ENRG + KION - 1 )
               ETH   = RSKWRK( IPYLD + 1)
               IF (EI .GE. ETH) THEN
                  ASCAL = RSKWRK( IPYLD)
                  A     = RSKWRK( IPYLD + 2)
                  B     = RSKWRK( IPYLD + 3)
                  COEF = COEF * ASCAL * (EI**A - ETH**A) ** B
               ELSE
                  COEF = ZERO
               ENDIF
            ENDIF
         ENDIF
C
         IF (COEF .NE. ZERO) THEN
            DKDAI(K) = DKDAI(K) + COEF*RKF / RA
C           if reverse parameters were not given, IREV=0
            IF (ISREV .EQ. 0) DKDAI(K) = DKDAI(K) - COEF*RKR / RA
         ENDIF
   80 CONTINUE
C
C     end of SUBROUTINE SKDRDA
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKDRDC (KSPEC, P, T, ACT, SDEN, ISKWRK, RSKWRK, DKDC)
C
C  START PROLOGUE
C
C  SUBROUTINE SKDRDC (KSPEC, P, T, ACT, SDEN, ISKWRK, RSKWRK, DKDC)
C  Returns the partial derivative of the production rates of the
C  species with respect to the concentration of species KSPEC.
C
C  INPUT
C  KSPEC     - Integer scalar, species index
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  DKDC(*)   - Real array, partial of the production rates of the
C              species with respect to the concentration of species
C              KSPEC;
C              dimension at least KKTOT, the total species count.
C                 cgs units, moles/(cm**2*sec) / (units of KSPEC)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0)
C*****END precision > single
C
      INCLUDE          'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ISKWRK(*), RSKWRK(*), T(*), ACT(*), SDEN(*), DKDC(*)
C
C     Zero out the result vector
C
      DO 50 K = 1, ISKWRK(IiKTOT)
         DKDC(K) = ZERO
   50 CONTINUE
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      SDTOT  = ZERO
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
         IF (ISKWRK(IiNYLD) .GT. 0) CALL SKEQY (ISKWRK, RSKWRK)
      ENDIF
C
C     initialize some local pointers in COMMON /SKPTR/
      CALL SKLOC (ISKWRK)
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, RSKWRK(I_CZ))
      CALL SKRROP (ISKWRK, RSKWRK, RSKWRK(I_SMH), T, RSKWRK(I_CZ),
     1             ACT, RSKWRK(I_KWT), ISKWRK(I_NREA), ISKWRK(I_NRPP),
     2             ISKWRK(I_NU), ISKWRK(I_NK), ISKWRK(I_IRNU),
     3             RSKWRK(I_RNU), ISKWRK(I_NUSUM), RSKWRK(I_PAR),
     4             RSKWRK(I_RPAR), ISKWRK(I_IREV), ISKWRK(I_ICOV),
     5             ISKWRK(I_KCOV), RSKWRK(I_CPAR), ISKWRK(I_ISTK),
     6             ISKWRK(I_MSTK), ISKWRK(I_IBHM), ISKWRK(I_KBHM),
     7             ISKWRK(I_IEDP), ISKWRK(I_KEDP), RSKWRK(I_PEDP),
     8             ISKWRK(I_KTFL), SDTOT, ISKWRK(I_IORD),
     9             ISKWRK(I_KORD), RSKWRK(I_ORD), ISKWRK(I_IYLD),
     *             ISKWRK(I_YION), ISKWRK(I_KYLD), RSKWRK(I_PYLD),
     1             RSKWRK(I_RKFT), RSKWRK(I_RKRT), RSKWRK(I_RKF),
     2             RSKWRK(I_RKR), RSKWRK(I_EQKC), RSKWRK(I_EQFAC),
     3             RSKWRK(I_ENRG))
C
C          Call a subroutine to calculate the partial derivatives
      CALL SKDRDC1 (KSPEC, P, T, ACT, SDEN, ISKWRK, RSKWRK, DKDC,
     1              RSKWRK(I_CZ), ISKWRK(I_NU), ISKWRK(I_NK),
     2              ISKWRK(I_IREV), ISKWRK(I_IRNU), RSKWRK(I_RNU),
     3              ISKWRK(I_ICOV), ISKWRK(I_KCOV), RSKWRK(I_CPAR),
     4              ISKWRK(I_IBHM), ISKWRK(I_KBHM), ISKWRK(I_IORD),
     5              ISKWRK(I_KORD), RSKWRK(I_ORD), ISKWRK(I_IYLD),
     6              RSKWRK(I_RKF), RSKWRK(I_RKR))
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKDRDC1(KSPEC, P, T, ACT, SDEN, ISKWRK, RSKWRK, DKDC,
     1                   CZ, NU, NUNK, IREV, IRNU, RNU, ICOV, KCOV,
     2                   CPAR, IBOHM, IBK, IORD, KORD, RORD, IYLD,
     3                   RKF, RKR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKDRDC1(KSPEC, P, T, ACT, SDEN, ISKWRK, RSKWRK, DKDC,
C                     CZ, NU, NUNK, IREV, IRNU, RNU, ICOV, KCOV,
C                     CPAR, IBOHM, IBK, IORD, KORD, RORD, RKF, RKR)
C
C  This is an internal routine that does most of the work for the user-
C  callable routine, SKDRDC.
C  Returns the partial derivative of the production rates of the
C  species with respect to the concentration of species KSPEC.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TEN=10.0D0, SMALLC=1.0D-20)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0, ONE=1.0, TEN=10.0, SMALLC=1.0E-20)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), NUNK(MAXSPR,*), NU(MAXSPR,*), IREV(*),
     1          IRNU(*), ICOV(*), KCOV(*), IBOHM(*), IBK(*), IORD(*),
     2          KORD(MAXORD,*), IYLD(*), T(*), ACT(*), SDEN(*),
     3          RSKWRK(*), DKDC(*), CZ(*), RNU(MAXSPR,*),
     4          CPAR(NSCOV,*), RORD(MAXORD,*), RKF(*), RKR(*)
C
      INTEGER          CKLKUP
      EXTERNAL         CKLKUP
C
C    Obtain local values of integer numbers in work arrays
C
      NIISUR = ISKWRK(IiNIIS)
      NIIRNU = ISKWRK(IiNRNU)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
      NIIYLD = ISKWRK(IiNYLD)
C
C     This is the concentration of the species in question
C     - protect against division by zero
C
      IF (CZ(KSPEC) .LE. ZERO) THEN
         CZK = SMALLC
      ELSE
         CZK = CZ(KSPEC)
      ENDIF
C
C     Process each reaction in a big outer loop
C
      DO 1000 I = 1, NIISUR
C
C        index for supplemental reverse rate data, if given
         ISREV = CKLKUP(I, IREV, NIIREV)
C        index for supplemental real stoichiometry, if given
         ISRNU = CKLKUP(I, IRNU, NIIRNU)
C        index for supplemental species change-orders, if given
         ISORD = CKLKUP(I, IORD, NIIORD)
C        index for supplemental Bohm-type data, if given
         ISBHM = CKLKUP(I, IBOHM, NIIBHM)
C        index for supplemental yield-modify data, if given
         ISYLD = CKLKUP(I, IYLD, NIIYLD)
C
C          Initialize the forward and reverse partial deriv factors
C
         FFACT = ZERO
         RFACT = ZERO
C
C          Main Loop for finding the factors
C
         IF (ISORD .EQ. 0) THEN
C         (artifice to start processing at the 1st and 7th nunk value)
          DO 350 NSTART = 1, 7, 6
           DO 300 N = NSTART, NSTART + 5
             NK = NUNK(N,I)
             IF (NK .EQ. 0) GO TO 325
             IF (NK .EQ. KSPEC) THEN
C
C                Check to see if this rxn has non-integer stoichiometry
C
               IF (ISRNU .GT. 0) THEN
                 RORDER = RNU(N,ISRNU)
               ELSE
                 RORDER = NU(N,I)
               END IF
C
C                Set the multiplicative factors for forward or reverse
C                directions, depending on whether the species appears
C                in the forward or reverse rate expression.
C                Irreversible reactions are handled by RKR_I=0.0.
C
               IF (RORDER .LT. ZERO) THEN
                  FFACT = -RORDER / CZK
               ELSE
                  RFACT =  RORDER / CZK
               ENDIF
             ENDIF
  300      CONTINUE
  325      CONTINUE
  350     CONTINUE
         ELSE
C
C           Handle reactions with different order vs stoichiometries
C           - They are packed tight, so can bail out if NK=0
C
          DO 360 N = 1, MAXORD
            NK = KORD(N,ISORD)
            IF (ABS(NK) .EQ. KSPEC) THEN
              IF (NK .LT. 0) THEN
                FFACT = RORD(N,ISORD) / CZK
              ELSE
                RFACT = RORD(N,ISORD) / CZK
              ENDIF
            ELSE IF (NK .EQ. 0) THEN
              GOTO 361
            ENDIF
  360     CONTINUE
  361     CONTINUE
         ENDIF
C
C          Special section for BOHM reactions
C            Restrictions observed in SKRROP:
C            - always integer stoichiometry
C            - reverse reaction rate is zero.
C
         IF (ISBHM .GT. 0) THEN
           IF (IBK(N) .EQ. KSPEC) THEN
             FFACT = ONE / CZK
           ELSE IF (KSPEC .GT. NKKGAS) THEN
             DO 365 N = 1, 6
               NK = NUNK(N,I)
               IF (NK .NE. IBK(N) .AND. NK .EQ. KSPEC) THEN
                 FFACT = -NU(N,I) / CZK
               ELSE IF (NK .EQ. 0) THEN
                 GO TO 366
               ENDIF
  365        CONTINUE
  366        CONTINUE
           ENDIF
         ENDIF
C
C
C          Fix up the factors for coverage dependent reaction
C          rate constants
         IF (NIICOV .GT. 0) THEN
            DO 370 N = 1, NIICOV
               IF (ICOV(N).EQ.I .AND. KCOV(N).EQ.KSPEC) THEN
C
C               Calc the coverage modification for the forward
C               direction
C                - ACT/CZK is the (sites per surface spec)
C                  divided by (density of sites for that phase)
C                   -- look at routine skatcz.
C
                 CMOD = (ACT(KSPEC)
     $                 * (CPAR(1,N)*LOG(TEN) - CPAR(3,N)/T(1))
     $                 + CPAR(2,N)) / CZK
C
C               Add the modification into the main factor
C
                 FFACT = FFACT + CMOD
C
C               Add in the coverage modifications to the reverse
C               direction, if reverse arrhenius parameters were
C               not given for the reaction.
C
                 IF (ISREV .EQ. 0) RFACT = RFACT + CMOD
              ENDIF
  370       CONTINUE
         ENDIF
C
C          Calculate the partial derivative of the reaction rate wrt
C          the concentration of the KSPEC species, TFACT
C
         TFACT = RKF(I) * FFACT  -  RKR(I) * RFACT
C
C          Now combine the value of TFACT with the stoichiometric
C          coefficients of the reaction to yield the partial derivatives
C          of source terms wrt the concentration of the KSPEC species
C
         DO 500 NSTART = 1, 7, 6
           DO 400 N = NSTART, NSTART + 5
C
C             Bail out of loop if no more stoich coeff. to process
C
             NK = NUNK(N,I)
             IF (NK .EQ. 0) GO TO 450
C
C              The regular stoichiometric coefficient can be
C              overriden in real var. stoich's are being used
C
             IF (ISRNU .GT. 0) THEN
               RSTOICH = RNU(N,ISRNU)
             ELSE
               RSTOICH = NU(N,I)
             ENDIF
C
C              Add in the contribution from this reaction to the
C              NKth source term
C
             DKDC(NK) = DKDC(NK) + RSTOICH * TFACT
  400      CONTINUE
  450      CONTINUE
  500    CONTINUE
C
 1000 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKDSDC (P, T, X, ACT, SDEN, ISKWRK, RSKWRK, DSDC,
     $                   KKTOT, SDOT, SITDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKDSDC (P, T, X, ACT, SDEN, ISKWRK, RSKWRK, DSDC, KKTOT,
C                     SDOT, SITDOT)
C  Returns the partial derivative of the production rates of the
C  species with respect to the concentration of each species.
C  It also returns the matching production rates.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fraction (or its equivalent) of the
C              species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS X are mole fractions,
C              the next KKSURF X are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  X(K)*SITE_DENSITY / # sites per species),
C              the next KKBULK X are bulk species mole fractions.
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C  KKTOT     - Integer scalar, total species count.
C
C  OUTPUT
C  DSDC(*,*) - Real matrix, the partial derivatives of the production
C              rates of the species with respect to the concentration
C              of species KSPEC;
C              dimension at least KKTOT, the total species count, for
C              both the first and second dimensions.
C                 cgs units, moles/(cm**2*sec) / (units of KSPEC)
C  SDOT(*)   - Real array, production rates of the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, moles/(cm**2*sec)
C              for 1,KKGAS, the production rates of gas-phase species,
C              for KKGAS+1,KKGAS+KKSUR, the production rates of surface
C              species,
C              for KKGAS+KKSUR+1,KKTOT, the production rate of bulk
C              species.
C  SITDOT(*) - Real array, production rates of the surface phases;
C              dimension at least NPHASE, the total phase count, but
C              subroutine only calculates entries for site phases.
C                 cgs units, moles/(cm**2*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0, SMALLC=1.0E-30)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0, SMALLC=1.0E-20)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ISKWRK(*), RSKWRK(*), T(*), X(*), ACT(*), SDEN(*),
     1          DSDC(KKTOT,KKTOT), SDOT(KKTOT), SITDOT(*)
C
C      Stop if KKTOT is wrong
C
      IF (KKTOT .NE. ISKWRK(IiKTOT)) STOP
C
C      Zero out the jacobian and SDOT result vector
C
      DO 50 J = 1, KKTOT
        SDOT(J) = ZERO
        DO 40 K = 1, KKTOT
          DSDC(K,J) = ZERO
   40   CONTINUE
   50 CONTINUE
C
C     Calculate the sum of all the site densities
C     - Zero the SITDOT() array at the same time.
C
      SDTOT  = ZERO
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SITDOT(N) = ZERO
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
         IF (ISKWRK(IiNYLD) .GT. 0) CALL SKEQY (ISKWRK, RSKWRK)
      ENDIF
C
C       Quick return if there are no surface reactions
C       - All of the output arrays have been zeroed
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
C     get local array pointers
      CALL SKLOC (ISKWRK)
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, RSKWRK(I_CZ))
C
C       The concentrations of the species
C       must be protected against division by zero
C       This version allows negative concentrations, because some
C       activities may be negative. - CZ comes from internal storage
C       so change values is OK.
C         Note - It is necessary to protect against a zero mole or
C                site fraction,  before the call to SKRROP.
C                Because, some of jacobian entries would be zero,
C                and this will cause a singular jacobian to occur.
C
      NKM1 = KKTOT - 1
      DO 70 K = 0, NKM1
        IF (RSKWRK(I_CZ+K) .EQ. ZERO) RSKWRK(I_CZ+K) = SMALLC
   70 CONTINUE
      CALL SKRROP (ISKWRK, RSKWRK, RSKWRK(I_SMH), T, RSKWRK(I_CZ),
     1             ACT, RSKWRK(I_KWT), ISKWRK(I_NREA), ISKWRK(I_NRPP),
     2             ISKWRK(I_NU), ISKWRK(I_NK), ISKWRK(I_IRNU),
     3             RSKWRK(I_RNU), ISKWRK(I_NUSUM), RSKWRK(I_PAR),
     4             RSKWRK(I_RPAR), ISKWRK(I_IREV), ISKWRK(I_ICOV),
     5             ISKWRK(I_KCOV), RSKWRK(I_CPAR), ISKWRK(I_ISTK),
     6             ISKWRK(I_MSTK), ISKWRK(I_IBHM), ISKWRK(I_KBHM),
     7             ISKWRK(I_IEDP), ISKWRK(I_KEDP), RSKWRK(I_PEDP),
     8             ISKWRK(I_KTFL), SDTOT, ISKWRK(I_IORD),
     9             ISKWRK(I_KORD), RSKWRK(I_ORD), ISKWRK(I_IYLD),
     *             ISKWRK(I_YION), ISKWRK(I_KYLD), RSKWRK(I_PYLD),
     1             RSKWRK(I_RKFT), RSKWRK(I_RKRT), RSKWRK(I_RKF),
     2             RSKWRK(I_RKR), RSKWRK(I_EQKC), RSKWRK(I_EQFAC),
     3             RSKWRK(I_ENRG))
C
C     Call a subroutine to calculate the partial derivatives
C     of the production rate, SDOT, with respect to the
C     concentrations
C
      I_DRDK = ISKWRK(IrKT2)
      CALL SKDSDC1 (P, T, ACT, SDEN, ISKWRK, RSKWRK, DSDC, KKTOT,
     1              RSKWRK(I_CZ), ISKWRK(I_NU), ISKWRK(I_NK),
     2              ISKWRK(I_IREV), ISKWRK(I_IRNU), RSKWRK(I_RNU),
     3              ISKWRK(I_ICOV), ISKWRK(I_KCOV), RSKWRK(I_CPAR),
     4              ISKWRK(I_IBHM), ISKWRK(I_KBHM), ISKWRK(I_IORD),
     5              ISKWRK(I_KORD), RSKWRK(I_ORD), RSKWRK(I_DRDK),
     6              RSKWRK(I_RKF), RSKWRK(I_RKR))
C
C     Now, calculate SDOT and SITDOT normally.
      I_NCF = ISKWRK(IrNCF)
      I_YNCF= ISKWRK(IrYNCF)
      CALL SKSDOT (ISKWRK(IiNIIS), RSKWRK(I_RKF), RSKWRK(I_RKR),
     1             ISKWRK, ISKWRK(I_NU), ISKWRK(I_NK),
     2             ISKWRK(I_IRNU), RSKWRK(I_RNU), RSKWRK(I_NCF),
     3             ISKWRK(I_IYLD), ISKWRK(I_KYLD), RSKWRK(I_PYLD),
     4             ISKWRK(I_YION), RSKWRK(I_ENRG), ISKWRK(IiNPHA),
     5             RSKWRK(I_YNCF), SDOT, SITDOT)
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
       SUBROUTINE SKDSDC1(P, T, ACT, SDEN, ISKWRK, RSKWRK, DSDC, KKTOT,
     1                    CZ, NU, NUNK, IREV, IRNU, RNU, ICOV, KCOV,
     2                    CPAR, IBOHM, IBK, IORD, KORD, RORD, DRopDK,
     3                    RKF, RKR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKDSDC1(P, T, ACT, SDEN, ISKWRK, RSKWRK, DSDC, KKTOT,
C                     CZ, NU, NUNK, IREV, IRNU, RNU, ICOV, KCOV,
C                     CPAR, IBOHM, IBK, IORD, KORD, RORD, DRopDK,
C                     RKF, RKR)
C
C  This is an internal routine that does most of the work for the user-
C  callable routines, SKDSDC and SKDSDX.
C  Returns the partial derivative of the production rates of the
C  species with respect to either the concentration activitity of
C  each species or the activity.
C  The difference between the DSDC and DSDX case is that in the
C  later case CZ is replaced with ACT in the calling routine, SKDSDX.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C  KKTOT     - Integer scalar, total species count.
C  CZ(*)     - Real array, gas-phase and surface species concentrations,
C              and bulk species activities;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS gas-phase concentrations are moles/cm**3,
C              the next KKSURF site concentrations are moles/cm**2,
C              and
C              the final KKBULK entries are bulk species activities.
C  NU(*,*)   - Integer matrix, stoichiometric coefficients for
C              species in surface reactions;
C              dimension at least MAXSPR for the first and at least
C              IISUR for the second.
C              NU(N,IR) is the stoichiometric coefficient of the Nth
C              species in reaction IR, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(*,*) - Integer matrix, indices of species in reactions;
C              dimension at least MAXSP for the first, and at least
C              IISUR for the second.
C              NUNK(N,IR) is the species index for the Nth species in
C              reaction IR.
C  IREV(*)   - Integer array, reaction indices for reactions with
C              explicit reverse parameters;
C              dimension at least NIIREV, the total count of reactions
C              with explicit reverse parameters.
C  IRNU(*)   - Integer array, reaction indices for reactions with real
C              stoichiometric coefficients;
C              dimension at least NIIRNU, the total count of reactions
C              with real stoichiometry.
C              Note - the NIIRNU reactions with have real stoichiometric
C              coefficients have their NU rows identically zero.
C  RNU(*,*)  - Real matrix, stoichiometric coefficients for reactions
C              with real stoichiometry;
C              dimension at least MAXSPR for the first, the maximum
C              number of species in a reaction, and at least NIIRNU for
C              the second, the total count of reactions with real
C              stoichiometry.
C  ICOV(*)   - Integer array, reaction indices for the NCOV reactions;
C              dimension at least NCOV.
C  KCOV(*)   - Integer array, coverage species indices for the NCOV
C              reactions;
C              dimension at least NCOV.
C  CPAR(*,*) - Real matrix, coverage parameters for the NCOV reactions;
C              dimension at least NSCOV for the first, the number of
C              coverage parameters allowed, and at least NCOV for the
C              second, the total coverage reaction count.
C  IBOHM(*)  - Integer array, reaction indices for the Bohm reactions;
C              dimension at least NBOHM, the total Bohm reaction count.
C  IBK(*)    - Integer array, species indices for the Bohm reactions;
C              dimension at least NBOHM, the total Bohm reaction count.
C  IORD(*)   - Integer array, reaction indices for the NORD reactions;
C              dimension at least NORD.
C              IORD(N) is the index of the Nth change-order reaction.
C  KORD(*,*) - Integer matrix, species indices for the order changes in
C              the NORD reactions; dimension at least MXORD for the
C              first and at least NORD for the second.
C              KORD(L,N) is the species index for the Lth order change
C              in the Nth change-order reaction.
C              KORD < 0 indicates change in forward order;
C              KORD > 0 indicates change in reverse order.
C  RORD(*,*) - Real matrix, order values for the NORD reactions;
C              dimension at least MXORD for the first and at least NORD
C              for the second.
C              RORD(L,N) is the order for the Lth order change in the
C              Nth change-order reaction.
C  DRopDK    - Derivative of the rate of progress of the current rxn
C              wrt the kth species (KKTOT long)
C  RKF(*)    - Real array, forward rates of progress for the surface
C              reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec)
C  RKR(*)    - Real array, reverse rates of progress for the surface
C              reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec)
C  OUTPUT
C  DSDC(*,*) - Real matrix, partial derivatives of the production
C              rates of the species with respect to the concentration
C              of species KSPEC;
C              dimension at least KKTOT, the total species count,
C              for both the first and second dimensions.
C                 cgs units, moles/(cm**2*sec) / (units of KSPEC)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TEN=10.0D0, SMALLC=1.0D-20)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0, ONE=1.0, TEN=10.0, SMALLC=1.0E-20)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), NUNK(MAXSPR,*), NU(MAXSPR,*),
     1          IREV(*), IRNU(*), ICOV(*), KCOV(*), IBOHM(*), IBK(*),
     2          IORD(*), KORD(MAXORD,*), T(*), ACT(KKTOT), SDEN(*),
     3          DSDC(KKTOT,KKTOT), CZ(KKTOT), RNU(MAXSPR,*),
     4          CPAR(NSCOV,*), RORD(MAXORD,*), DRopDK(KKTOT), RKF(*),
     5          RKR(*)
C
      EXTERNAL CKLKUP
      INTEGER  CKLKUP
C
C       Obtain local values of integer numbers in work arrays
C
      NIISUR = ISKWRK(IiNIIS)
      NIIRNU = ISKWRK(IiNRNU)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
      IF (NIICOV .GT. 0) TENLOG = LOG(TEN)
C
C        Zero the rate of progress DK vector
C
      DO 50 K = 1, KKTOT
        DRopDK(K) = ZERO
 50   CONTINUE
C
C  Process each reaction in the outer loop
C
      DO 1000 I = 1, NIISUR
C
C           Store the local forward and reverse reaction rates
C
C        index for supplemental reverse rate data, if given
         ISREV = CKLKUP(I, IREV, NIIREV)
C        index for supplemental real stoichiometry, if given
         ISRNU = CKLKUP(I, IRNU, NIIRNU)
C        index for supplemental species change-orders, if given
         ISORD = CKLKUP(I, IORD, NIIORD)
C        index for supplemental Bohm-type data, if given
         ISBHM = CKLKUP(I, IBOHM, NIIBHM)
C
C          Main Loop for finding the factors
C
         IF (ISORD .EQ. 0) THEN
C         (artifice to start processing at the 1st and 7th nunk value)
          DO 350 NSTART = 1, 7, 6
           DO 300 N = NSTART, NSTART + 5
             NK = NUNK(N,I)
             IF (NK .EQ. 0) GO TO 325
C
C              Check to see if this rxn has non-integer stoichiometry
C
             IF (ISRNU .GT. 0) THEN
               RORDER = RNU(N,ISRNU)
             ELSE
               RORDER = NU(N,I)
             END IF
C
C              Set the multiplicative factors for forward or reverse
C              directions, depending on whether the species appears
C              in the forward or reverse rate expression.
C              Irreversible reactions are handled by RKR_I=0.0.
C
             IF (RORDER .LT. ZERO) THEN
                DRopDK(NK) = DRopDK(NK) - RKF(I)*RORDER/CZ(NK)
             ELSE
                DRopDK(NK) = DRopDK(NK) - RKR(I)*RORDER/CZ(NK)
             ENDIF
  300      CONTINUE
  325      CONTINUE
  350     CONTINUE
         ELSE
C
C           Handle reactions with different order vs stoichiometries
C           - They are packed tight, so can bail out if NK=0
C
          DO 360 N = 1, MAXORD
            NK = KORD(N,ISORD)
            IF (NK .EQ. 0) GOTO 361
            IF (NK .LT. 0) THEN
C                    NK is a reactant
              NK = ABS(NK)
              DRopDK(NK) = DRopDK(NK) + RKF(I)*RORD(N,ISORD)/CZ(NK)
            ELSE
              DRopDK(NK) = DRopDK(NK) - RKR(I)*RORD(N,ISORD)/CZ(NK)
            ENDIF
  360     CONTINUE
  361     CONTINUE
         ENDIF
C
C          Special section for BOHM reactions
C            Restrictions observed in SKRROP:
C            - always integer stoichiometry
C            - reverse reaction rate is zero.
C
         IF (ISBHM .GT. 0) THEN
           K = IBK(N)
           DRopDK(K) = DRopDK(K) + RKF(I)/CZ(K)
           DO 365 N = 1, 6
             NK = NUNK(N,I)
             IF (NK .EQ. 0) GOTO 366
             IF (NK .GT. NKKGAS) THEN
               DRopDK(NK) = DRopDK(NK) - RKF(I)*NU(N,I)/CZ(NK)
             ENDIF
  365      CONTINUE
  366      CONTINUE
         ENDIF
C
C          Fix up the factors for coverage dependent reaction
C          rate constants
C
          IF (NIICOV .GT. 0) THEN
              DO 370 N = 1, NIICOV
                 IF (ICOV(N) .EQ. I) THEN
C                  coverage dependence species index
                   NK = KCOV(N)
C
C                  Calc the coverage modification for the forward
C                  direction
C                - ACT/CZ is the (sites per surface spec)
C                  divided by (density of sites for that phase)
C                   -- look at routine skatcz.
C
                   CMOD = (ACT(NK)
     $                 * (CPAR(1,N)*TENLOG - CPAR(3,N)/T(1))
     $                 + CPAR(2,N)) / CZ(NK)
C
C                  Add the modification into the main factor
C
                   DRopDK(NK) = DRopDK(NK) + RKF(I)*CMOD
C
C                  Add in the coverage modifications to the reverse
C                  direction, if reverse arrhenius parameters were
C                  not given for the reaction.
C
                 IF (ISREV .EQ. 0) DRopDK(NK) = DRopDK(NK) - RKR(I)*CMOD
               ENDIF
  370       CONTINUE
         ENDIF
C
C        Obtain the final partial derivatives by first summing over
C        the column entry.  Check to see if it is zero first, since
C        most entries will be zero.
C
         DO 600 K = 1, KKTOT
           IF (DRopDK(K) .NE. ZERO) THEN
C
C          Now combine the value of TFACT with the stoichiometric
C          coefficients of the reaction to yield the partial
C          derivatives of the source terms wrt the concentration
C          of the KSPEC species
C
             DO 400 N = 1, 6
C
C               Bail out of loop if there are no more stoich coeff.
C               to process
C
               NK = NUNK(N,I)
               IF (NK .EQ. 0) GO TO 450
C
C              The regular stoichiometric coefficient can be
C              overriden in real var. stoich's are being used
C
               IF (ISRNU .GT. 0) THEN
                 RSTOICH = RNU(N,ISRNU)
               ELSE
                 RSTOICH = NU(N,I)
               ENDIF
C
C              Add in the contribution from this reaction to the
C              NKth source term
C
               DSDC(NK,K) = DSDC(NK,K) + RSTOICH * DRopDK(K)
  400        CONTINUE
  450        CONTINUE
             DO 500 N = 7, 12
C
C                Bail out of loop if there are no more stoich coeff.
C                to process
C
               NK = NUNK(N,I)
               IF (NK .EQ. 0) GO TO 550
C
C                The regular stoichiometric coefficient can be
C                overriden in real var. stoich's are being used
C
               IF (ISRNU .GT. 0) THEN
                 RSTOICH = RNU(N,ISRNU)
               ELSE
                 RSTOICH = NU(N,I)
               ENDIF
C
C                Add in the contribution from this reaction to the
C                NKth source term
C
               DSDC(NK,K) = DSDC(NK,K) + RSTOICH * DRopDK(K)
  500        CONTINUE
  550        CONTINUE
C
C              Initialize the non-zero entry for the next reaction
C
             DRopDK(K) = ZERO
           ENDIF
  600    CONTINUE
C
 1000 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKDSDX (P, T, X, ACT, SDEN, ISKWRK, RSKWRK, DSDX,
     $                   KKTOT, SDOT, SITDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKDSDX (P, T, X, ACT, SDEN, ISKWRK, RSKWRK, DSDX, KKTOT,
C                     SDOT, SITDOT)
C  Returns the partial derivative of the production rates of the
C  species with respect to the activity for each species.
C  It also returns the matching production rates.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fraction (or its equivalent) of the
C              species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS X are mole fractions,
C              the next KKSURF X are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  X(K)*SITE_DENSITY / # sites per species),
C              the next KKBULK X are bulk species mole fractions.
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C  KKTOT     - Integer scalar, total species count.
C
C  OUTPUT
C  DSDX(*,*) - Real matrix, partial derivatives of the production rates
C              of the species with respect to the activity of species
C              KSPEC;
C              dimension at least KKTOT, the total species count, for
C              both the first and second dimensions.
C                 cgs units, moles/(cm**2*sec) / (units of KSPEC)
C  SDOT(*)   - Real array, production rates of the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, moles/(cm**2*sec)
C              SDOT(K) is
C              for 1,KKGAS, the production rate of gas-phase species,
C              for KKGAS+1,KKGAS+KKSUR, the production rate of surface
C              species,
C              for KKGAS+KKSUR+1,KKTOT, the production rate of bulk
C              species.
C  SITDOT(*) - Real array, production rates of the surface phases;
C              dimension at least NPHASE, the total phase count, but
C              subroutine only calculates entries for site phases.
C                 cgs units, moles/(cm**2*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0, SMALLC=1.0D-30)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0, SMALLC=1.0E-30)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ISKWRK(*), RSKWRK(*), T(*), X(*), ACT(*), SDEN(*),
     1          DSDX(KKTOT,KKTOT), SDOT(KKTOT), SITDOT(*)
C
C     Stop if KKTOT is wrong
C
      IF (KKTOT .NE. ISKWRK(IiKTOT)) STOP
C
C      Zero out the jacobian and SDOT result vector
C
      DO 50 J = 1, KKTOT
         SDOT(J) = ZERO
         DO 40 K = 1, KKTOT
           DSDX(K,J) = ZERO
   40    CONTINUE
   50 CONTINUE
C
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 20 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SITDOT(N) = ZERO
   20    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      SDTOT  = ZERO
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
         IF (ISKWRK(IiNYLD) .GT. 0) CALL SKEQY (ISKWRK, RSKWRK)
      ENDIF
C
C     initialize some local pointers in COMMON /SKPTR/
      CALL SKLOC (ISKWRK)
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, RSKWRK(I_CZ))
C
C     The concentrations of the species must be protected against
C     against division by zero
C     This version allows negative concentrations, because some
C     activities may be negative. - CZ comes from internal storage
C     so change values is OK.
C         Note - It is necessary to protect against a zero mole or
C                site fraction,  before the call to SKRROP.
C                Because, some of jacobian entries would be zero,
C                and this will cause a singular jacobian to occur.
C
      NKM1 = KKTOT - 1
      DO 70 K = 0, NKM1
        IF (RSKWRK(I_CZ+K) .EQ. ZERO) RSKWRK(I_CZ+K) = SMALLC
   70 CONTINUE
C
C     Calculate the forward and reverse rates for the surface mech
      CALL SKRROP (ISKWRK, RSKWRK, RSKWRK(I_SMH), T, RSKWRK(I_CZ),
     1             ACT, RSKWRK(I_KWT), ISKWRK(I_NREA), ISKWRK(I_NRPP),
     2             ISKWRK(I_NU), ISKWRK(I_NK), ISKWRK(I_IRNU),
     3             RSKWRK(I_RNU), ISKWRK(I_NUSUM), RSKWRK(I_PAR),
     4             RSKWRK(I_RPAR), ISKWRK(I_IREV), ISKWRK(I_ICOV),
     5             ISKWRK(I_KCOV), RSKWRK(I_CPAR), ISKWRK(I_ISTK),
     6             ISKWRK(I_MSTK), ISKWRK(I_IBHM), ISKWRK(I_KBHM),
     7             ISKWRK(I_IEDP), ISKWRK(I_KEDP), RSKWRK(I_PEDP),
     8             ISKWRK(I_KTFL), SDTOT, ISKWRK(I_IORD),
     9             ISKWRK(I_KORD), RSKWRK(I_ORD), ISKWRK(I_IYLD),
     *             ISKWRK(I_YION), ISKWRK(I_KYLD), RSKWRK(I_PYLD),
     1             RSKWRK(I_RKFT), RSKWRK(I_RKRT), RSKWRK(I_RKF),
     2             RSKWRK(I_RKR), RSKWRK(I_EQKC), RSKWRK(I_EQFAC),
     3             RSKWRK(I_ENRG))
C
C     Change the concentration vector that we will feed into
C     SKDSDC1 into an activity vector.  Again, protect against
C     division by zero.
C
      DO 170 K = 1, KKTOT
         IF (ACT(K) .EQ. ZERO) THEN
            RSKWRK(I_CZ + K - 1) = SMALLC
         ELSE
            RSKWRK(I_CZ + K - 1) = ACT(K)
         ENDIF
  170 CONTINUE
C
C     Call a subroutine to calculate the partial derivatives
C     of the production rate, SDOT, with respect to the activity.
C         - We are making use here of the simple transformation
C           between the concentration activity vector and the
C           activity vector.
C         - Replace the ACT vector as well, with the non-zero
C           CZ vector -> necessary for coverage dependence
C           section.
C
C
      I_DRDK = ISKWRK(IrKT2)
      CALL SKDSDC1 (P, T, RSKWRK(I_CZ), SDEN, ISKWRK, RSKWRK, DSDX,
     1              KKTOT, RSKWRK(I_CZ), ISKWRK(I_NU), ISKWRK(I_NK),
     2              ISKWRK(I_IREV), ISKWRK(I_IRNU), RSKWRK(I_RNU),
     3              ISKWRK(I_ICOV), ISKWRK(I_KCOV), RSKWRK(I_CPAR),
     4              ISKWRK(I_IBHM), ISKWRK(I_KBHM), ISKWRK(I_IORD),
     5              ISKWRK(I_KORD), RSKWRK(I_ORD), RSKWRK(I_DRDK),
     6              RSKWRK(I_RKF), RSKWRK(I_RKR))
C
C      Now, calculate SDOT and SITDOT normally.
      I_NCF = ISKWRK(IrNCF)
      I_YNCF = ISKWRK(IrYNCF)
      CALL SKSDOT (ISKWRK(IiNIIS), RSKWRK(I_RKF), RSKWRK(I_RKR),
     1             ISKWRK, ISKWRK(I_NU), ISKWRK(I_NK), ISKWRK(I_IRNU),
     2             RSKWRK(I_RNU), RSKWRK(I_NCF), ISKWRK(I_IYLD),
     3             ISKWRK(I_KYLD), RSKWRK(I_PYLD), ISKWRK(I_YION),
     4             RSKWRK(I_ENRG), ISKWRK(IiNPHA), RSKWRK(I_YNCF),
     5             SDOT, SITDOT)
C
      RETURN
      END
C
      SUBROUTINE SKEQ   (P, T, ACT, SDEN, ISKWRK, RSKWRK, EQKC)
C
C  START PROLOGUE
C
C  SUBROUTINE SKEQ    (P, T, ACT, SDEN, ISKWRK, RSKWRK, EQKC)
C  Returns the equilibrium constants for the surface reactions given
C  pressure, temperature, species activities, and the site densities.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  EQKC(*)   - Real array, equilibrium constants in concentration units
C              for the reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, depends on reaction (moles, cm)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO = 0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ISKWRK(*), RSKWRK(*), ACT(*), SDEN(*), EQKC(*), T(*)
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      SDTOT  = ZERO
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
         IF (ISKWRK(IiNYLD) .GT. 0) CALL SKEQY (ISKWRK, RSKWRK)
      ENDIF
C
C     initialize some local pointers in COMMON /SKPTR/
      CALL SKLOC (ISKWRK)
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, RSKWRK(I_CZ))
      CALL SKRROP (ISKWRK, RSKWRK, RSKWRK(I_SMH), T, RSKWRK(I_CZ),
     1             ACT, RSKWRK(I_KWT), ISKWRK(I_NREA), ISKWRK(I_NRPP),
     2             ISKWRK(I_NU), ISKWRK(I_NK), ISKWRK(I_IRNU),
     3             RSKWRK(I_RNU), ISKWRK(I_NUSUM), RSKWRK(I_PAR),
     4             RSKWRK(I_RPAR), ISKWRK(I_IREV), ISKWRK(I_ICOV),
     5             ISKWRK(I_KCOV), RSKWRK(I_CPAR), ISKWRK(I_ISTK),
     6             ISKWRK(I_MSTK), ISKWRK(I_IBHM), ISKWRK(I_KBHM),
     7             ISKWRK(I_IEDP), ISKWRK(I_KEDP), RSKWRK(I_PEDP),
     8             ISKWRK(I_KTFL), SDTOT, ISKWRK(I_IORD),
     9             ISKWRK(I_KORD), RSKWRK(I_ORD), ISKWRK(I_IYLD),
     *             ISKWRK(I_YION), ISKWRK(I_KYLD), RSKWRK(I_PYLD),
     1             RSKWRK(I_RKFT), RSKWRK(I_RKRT), RSKWRK(I_RKF),
     2             RSKWRK(I_RKR), EQKC, RSKWRK(I_EQFAC),
     3             RSKWRK(I_ENRG))
C
C     end of SUBROUTINE SKEQ
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKEQY (ISKWRK, RSKWRK)
C
C  START PROLOGUE
C
C  SUBROUTINE SKEQY (ISKWRK, RSKWRK)
C  Not normally called by a user, this subroutine resets the multiplica-
C  tive factor of the equilibrium constant due to the yield-modify
C  option in conjunction with ion energy, therefore is called before
C  SKRROP which uses the constant.
C
C  END PROLOGUE

C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      NIIYLD = ISKWRK(IiNYLD)
      NNSUR  = ISKWRK(IiNSUR)
C
      IF (NIIYLD .LE. 0 .OR. NNSUR.LE.0) RETURN
      NIIRNU = ISKWRK(IiNRNU)
      NFSUR  = ISKWRK(IiFSUR)
      NLSUR  = ISKWRK(IiLSUR)
      SKMIN  = RSKWRK(IrSKMN)
C
      DO 250 L = 1, NIIYLD
C        loop over yield-modify reactions
         I = ISKWRK (ISKWRK(IiIYLD) + L - 1)
C        locator for species indices
         INK = ISKWRK(IiNUNK) + (I-1)*MAXSPR - 1
C        locator for reactions integeger stoichiometry
         I_NU = ISKWRK(IiNU)   + (I-1)*MAXSPR - 1
C        locator for the equilibrium factor
         IEQ = ISKWRK(IrEQ)   + (I-1)
C
C        actual species index for the ion
         IKION = ISKWRK (ISKWRK(IiYION) + L - 1)
C        locator for reaction species yield-modify flags
         IKYLD = ISKWRK(IiKYLD) + MAXSPR*(L-1)
C        locator for yield-modify parameters
         IPYLD = ISKWRK(IrPYLD) + NYPAR *(L-1)
C        does this reaction have real coefficients
         ISRNU  = CKLKUP (I, ISKWRK(ISKWRK(IiIRNU)), NIIRNU)
C
         DO 240 N = NFSUR, NLSUR
C           loop over sites
C
C           first site species index
            KST = ISKWRK( ISKWRK(IiPKST) + N - 1)
C           last site species index
            KND = ISKWRK( ISKWRK(IiPKND) + N - 1)
C           site density
            SDEN = RSKWRK( ISKWRK(IrSDEN) + N - 1)
C
            RNUSUM = 0.0
C           loop over species in sites
            DO 230 K = KST, KND
C              site coverage for this species
               RKCOV = ISKWRK( ISKWRK(IiNSCV) + K - 1)
C
               DO 220 M = 1, MAXSPR
C                 is this species in the reaction
                  KSPEC = ISKWRK(INK + M)
                  IF (K .NE. KSPEC) GO TO 220
C
C                 stoichiometric coefficient
                  IF (ISRNU .GT. 0) THEN
                     STOICH =
     1               RSKWRK( ISKWRK(IrRNU) + (ISRNU-1)*MAXSPR + M - 1)
                  ELSE
                     STOICH = ISKWRK(I_NU + M)
                  ENDIf
C
C                 yield-modify flag
                  IF (ISKWRK(IKYLD + M - 1) .GT. 0) THEN
                     EI = RSKWRK(ISKWRK(IrENGI) + IKION - 1)
                     ETH = RSKWRK( IPYLD + 1)
                     IF (EI .LT. ETH) THEN
                        STOICH = 0.0
                     ELSE
                        ASCAL = RSKWRK( IPYLD)
                        A = RSKWRK( IPYLD + 2)
                        B = RSKWRK( IPYLD + 3)
                        STOICH = STOICH * ASCAL * (EI**A - ETH**A)**B
                     ENDIF
                  ENDIF
                  IF (ABS(STOICH) .GT. SKMIN) THEN
                     RSKWRK(IEQ) = RSKWRK(IEQ) * RKCOV**(-STOICH)
                     RNUSUM = RNUSUM + STOICH
                  ENDIF
  220          CONTINUE
  230       CONTINUE
C
            RSKWRK(IEQ) = RSKWRK(IEQ) * SDEN**RNUSUM
  240    CONTINUE
  250 CONTINUE
C
C     end of SUBROUTINE SKEQY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKFLGS (IR, ISKWRK, NRPP, IREV, ISTFL, ICOV, IMOTZ,
     1                   IEDP, IBHM, IORD, IYLD)
C
C  START PROLOGUE
C
C  SUBROUTINE SKFLGS (IR, ISKWRK, NRPP, IREV, ISTFL, ICOV, IMOTZ,
C                     IEDP, IBHM, IORD, IYLD)
C  Returns several integer flags describing surface reaction IR.
C
C
C  INPUT
C  IR        - Integer scalar, surface reaction index.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C
C  OUTPUT
C  NRPP      - Integer scalar, number of species (reactants+products)
C              for surface reaction IR, combined with reversibility
C              flag.
C              NRPP > 0, NRPP species, reversible surface reaction,
C                  < 0, ABS(NRPP) species, irreversible reaction.
C  IREV      - Integer scalar, flag for explicit reverse Arrhenius
C              parameters.
C              =1, reaction has explicit reverse Arrhenius parameters
C              =0, no (may or may not be reversible, see NRPP).
C  ISTFL     - Integer scalar, flag for sticking coefficients;
C              =1, reaction does not use sticking coefficients
C              =0, no
C  IMOTZ     - Integer scalar, flag for Motz-Wise correction of
C              sticking coefficients;
C              =1, sticking reaction with Motz-Wise correction
C              =0, no (may or may not be sticking reaction, see ISTFL)
C  ICOV      - Integer scalar, flag to indidicate that reaction has
C              coverage dependence;
C              =1, reaction has coverage dependence
C              =0, no.
C  IEDP      - Integer scalar, flag for energy-dependence;
C              =1, reaction is energy-dependent,
C              =0, no.
C  IBHM      - Integer scalar, flag for Bohm correction;
C              =1, Bohm reaction,
C              =0, no
C  IORD      - Integer scalar, flag for species order change;
C              =1, reaction has species order change,
C              =0, no
C  IYLD      - Integer scalar, flag for yield-modification;
C              =1, yield-modification in reaction;
C              =0, no
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      NRPP = 0
      IREV = 0
      ISTFL = 0
      IMOTZ = 0
      ICOV = 0
      IEDP = 0
      IBHM = 0
      IORD = 0
      IYLD = 0
      IF (ISKWRK(IiNIIS) .LE. 0)   RETURN
C
      NRPP = ISKWRK(ISKWRK(IiNRPP) + IR - 1)
      IF (ISKWRK(IiNREV) .GT. 0) IREV = MIN (1,
     1       CKLKUP (IR, ISKWRK(ISKWRK(IiIREV)), ISKWRK(IiNREV)) )
      IF (ISKWRK(IiNSTK) .GT. 0) ISSTK =
     1       CKLKUP(IR, ISKWRK(ISKWRK(IiISTK)), ISKWRK(IiNSTK))
      IF (ISSTK .GT. 0) THEN
         ISTFL = 1
         IMOTZ = ISKWRK(ISKWRK(IiMSTK) + ISSTK - 1)
      ENDIF
C
      IF (ISKWRK(IiNCOV) .GT. 0) ICOV = MIN (1,
     1       CKLKUP (IR, ISKWRK(ISKWRK(IiICOV)), ISKWRK(IiNCOV)) )
      IF (ISKWRK(IiNEDP) .GT. 0) IEDP = MIN (1,
     1       CKLKUP (IR, ISKWRK(ISKWRK(IiIEDP)), ISKWRK(IiNEDP)) )
      IF (ISKWRK(IiNBHM) .GT. 0) IBHM = MIN (1,
     1       CKLKUP (IR, ISKWRK(ISKWRK(IiIBHM)), ISKWRK(IiNBHM)) )
      IF (ISKWRK(IiNORD) .GT. 0) IORD = MIN (1,
     1       CKLKUP (IR, ISKWRK(ISKWRK(IiIORD)), ISKWRK(IiNORD)) )
      IF (ISKWRK(IiNYLD) .GT. 0) IYLD = MIN (1,
     1       CKLKUP (IR, ISKWRK(ISKWRK(IiIYLD)), ISKWRK(IiNYLD)) )
C
C     end of SUBROUTINE SKFLGS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKGML  (T, ISKWRK, RSKWRK, GML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKGML  (T, ISKWRK, RSKWRK, GML)
C  Returns an array of the standard state Gibbs free energies
C  in molar units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  GML(*)    - Real array, standard state Gibbs free energies
C              for the species;
C              dimension KKTOT, the total species count.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), GML(*), T(*)
C
      IHML = ISKWRK(IrKT1)
      ISML = ISKWRK(IrKT2)
      CALL SKHML (T, ISKWRK, RSKWRK, RSKWRK(IHML))
      CALL SKSML (T, ISKWRK, RSKWRK, RSKWRK(ISML))
C
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 100 K = 0, NKM1
         IF (K .LT. NKKGAS) THEN
            TK = T(ISKWRK(ISKWRK(IiKTFL) + K))
         ELSE
            TK = T(1)
         ENDIF
         GML(K+1) = RSKWRK(IHML + K) - TK*RSKWRK(ISML + K)
100   CONTINUE
C
C     end of SUBROUTINE SKGML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKGMS  (T, ISKWRK, RSKWRK, GMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKGMS  (T, ISKWRK, RSKWRK, GMS)
C  Returns an array of the standard state Gibbs free energies
C  in mass units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  GMS(*)    - Real array, standard state Gibbs free energies
C              for the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), GMS(*), T(*)
C
      I_HMS = ISKWRK(IrKT1)
      I_SMS = ISKWRK(IrKT2)
      CALL SKHMS (T, ISKWRK, RSKWRK, RSKWRK(I_HMS))
      CALL SKSMS (T, ISKWRK, RSKWRK, RSKWRK(I_SMS))
C
      I_KTFL = ISKWRK(IiKTFL)
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 100 K = 0, NKM1
         IF (K .LT. NKKGAS) THEN
            TK = T(ISKWRK(I_KTFL + K))
         ELSE
            TK = T(1)
         ENDIF
         GMS(K+1) = RSKWRK(I_HMS + K) - TK*RSKWRK(I_SMS + K)
100   CONTINUE
C
C     end of SUBROUTINE SKGMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKHML  (T, ISKWRK, RSKWRK, HML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKHML  (T, ISKWRK, RSKWRK, HML)
C  Returns an array of the enthalpies in molar units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  HML(*)    - Real array, enthalpies for the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), HML(*), T(*)
      SAVE TN1, TN2, TN3, TN4, TN5
      DATA TN1/0.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN5 = TN1*TN4
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
         TN5 = TN5 / 5
      ENDIF
      I_KTFL = ISKWRK(IiKTFL)
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 250 K = 0, NKM1
         IF (K.LT.NKKGAS .AND. ISKWRK(I_KTFL+K).GT.1) THEN
            TK1 = T(ISKWRK(I_KTFL+K))
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK5 = TK1*TK4
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
            TK5 = TK5 / 5
         ELSE
            TK1 = TN1
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
            TK5 = TN5
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ISKWRK(ISKWRK(IiKNT)+K) - 1
C        location of FIRST set of thermodynamic coefficients
         NA1 = ISKWRK(IrKTHM) + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = ISKWRK(IrKTMP) + K*MAXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RSKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         HML(K+1) = RSKWRK(ISKWRK(IrRU)) * (RSKWRK(NA1)*TK1
     1            + RSKWRK(NA1+1) * TK2 + RSKWRK(NA1+2) * TK3
     2            + RSKWRK(NA1+3) * TK4 + RSKWRK(NA1+4) * TK5
     3            + RSKWRK(NA1 + NCP1 - 1))
250   CONTINUE
C
C     end of SUBROUTINE SKHML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKHMS  (T, ISKWRK, RSKWRK, HMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKHMS  (T, ISKWRK, RSKWRK, HMS)
C  Returns an array of the enthalpies in mass units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  HMS(*)    - Real array, enthalpies for the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), HMS(*), T(*)
      SAVE TN1, TN2, TN3, TN4, TN5
      DATA TN1/0.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN5 = TN1*TN4
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
         TN5 = TN5 / 5
      ENDIF
      I_KTFL = ISKWRK(IiKTFL)
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 250 K = 0, NKM1
         IF (K.LT.NKKGAS .AND. ISKWRK(I_KTFL+K).GT.1) THEN
            TK1 = T(ISKWRK(I_KTFL+K))
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK5 = TK1*TK4
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
            TK5 = TK5 / 5
         ELSE
            TK1 = TN1
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
            TK5 = TN5
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ISKWRK(ISKWRK(IiKNT)+K) - 1
C        location of FIRST set of thermodynamic coefficients
         NA1 = ISKWRK(IrKTHM) + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = ISKWRK(IrKTMP) + K*MAXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RSKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         HMS(K+1) = RSKWRK(ISKWRK(IrRU)) * (RSKWRK(NA1) * TK1
     1            + RSKWRK(NA1+1) * TK2 + RSKWRK(NA1+2) * TK3
     2            + RSKWRK(NA1+3) * TK4 + RSKWRK(NA1+4) * TK5
     3            + RSKWRK(NA1 + NCP1 - 1)) / RSKWRK(ISKWRK(IrKWT)+K)
250   CONTINUE
C
C     end of SUBROUTINE SKHMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKHORT (T, ISKWRK, RSKWRK, HORT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKHORT (T, ISKWRK, RSKWRK, HORT)
C  Returns an array of the nondimensional enthalpies.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  HORT(*)   - Real array, nondimensional enthalpies for the species;
C              dimension at least KKTOT, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), HORT(*), T(*)
      SAVE TN1, TNHALF, TN2, TN3, TN4
      DATA TN1/0.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNHALF = TN1 / 2
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2/3
         TN3 = TN3/4
         TN4 = TN4/5
      ENDIF
      I_KTFL = ISKWRK(IiKTFL)
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 250 K = 0, NKM1
         IF (K.LT.NKKGAS .AND. ISKWRK(I_KTFL+K).GT.1) THEN
            TK1 = T(ISKWRK(I_KTFL+K))
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TKHALF = TK1 / 2
            TK2 = TK2/3
            TK3 = TK3/4
            TK4 = TK4/5
         ELSE
            TK1 = TN1
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
            TKHALF = TNHALF
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ISKWRK(ISKWRK(IiKNT) + K) - 1
C        location of FIRST set of thermodynamic coefficients
         NA1 = ISKWRK(IrKTHM) + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = ISKWRK(IrKTMP) + K*MAXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RSKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         HORT(K+1) = RSKWRK(NA1)
     1             + RSKWRK(NA1+1) * TKHALF  + RSKWRK(NA1+2) * TK2
     2             + RSKWRK(NA1+3) * TK3 + RSKWRK(NA1+4) * TK4
     3             + RSKWRK(NA1 + NCP1 - 1) / TK1
250   CONTINUE
C
C     end of SUBROUTINE SKHORT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKIBHM (IR, ISKWRK, IBMFL)
C
C  START PROLOGUE
C
C  SUBROUTINE SKIBHM (IR, ISKWRK, IBMFL)
C  Returns an integer flag to indicate whether reaction IR uses
C  BOHM sticking coefficients.
C
C  INPUT
C  IR        - Integer scalar, surface reaction index.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C
C  OUTPUT
C  IBMFL     - Integer scalar, flag for Bohm reactions;
C              0, reaction IR does not use BOHM sticking coefficients
C              K, reaction IR does use BOHM sticking coefficients,
C                 and K is the index of the BOHM-correction ion.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C

      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      ISBHM = CKLKUP (IR, ISKWRK(ISKWRK(IiIBHM)), ISKWRK(IiNBHM))
      IF (ISBHM .GT. 0) THEN
         IBMFL = ISKWRK(ISKWRK(IiKBHM) + ISBHM - 1)
      ELSE
         IBMFL = 0
      ENDIF
C
C     end of SUBROUTINE SKIBHM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKICOV (IR, NDIM, ISKWRK, RSKWRK, NCOVI, KCOVI, CPARI)
C
C  START PROLOGUE
C
C  SUBROUTINE SKICOV (IR, NDIM, ISKWRK, RSKWRK, NCOVI, KCOVI, CPARI)
C  Returns the coverage species index numbers and their coverage
C  parameters for reaction IR.
C
C  INPUT
C  IR        - Integer scalar, surface reaction index.
C  NDIM      - Integer scalar, first dimension of array CPAR, the
C              coverage parameters; NDIM must be at least NSCOV,
C              the total number of coverage parameters.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  NCOVI     - Integer scalar, total number of species that modify the
C              rate of reaction IR through coverage dependence.
C  KCOVI(*)  - Integer array, species indices for the NCOVI species that
C              modify the rate of a coverage dependence reaction;
C              dimension at least KKTOT, the total species count.
C  CPARI(*,*)- Real matrix, coverage parameters for the coverage species
C              of reaction IR;
C              dimension at least NSCOV for the first, the number of
C              coverage parameters required, and at least KKTOT for the
C              second, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), KCOVI(*), CPARI(NDIM,*)
C
      NCOVI = 0
      DO 10 K = 1, ISKWRK(IiKTOT)
         KCOVI(K) = 0
         DO 05 N = 1, NSCOV
            CPARI(N, K) = 0.0
   05    CONTINUE
   10 CONTINUE

      DO 50 N = 0, ISKWRK(IiNCOV) - 1
         IF (IR .EQ. ISKWRK(ISKWRK(IiICOV) + N)) THEN
            NCOVI = NCOVI + 1
            KCOVI(NCOVI) = ISKWRK(ISKWRK(IiKCOV) + N)
            DO 25 J = 1, NSCOV
               CPARI(J, NCOVI) =
     1            RSKWRK(ISKWRK(IrKCOV) + NSCOV*N + J - 1)
   25       CONTINUE
         ENDIF
   50 CONTINUE
C
C     end of SUBROUTINE SKICOV
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKIENR (IR, ISKWRK, RSKWRK, IENRFL, IEION, PEDEP)
C
C  START PROLOGUE
C
C  SUBROUTINE SKIENR (IR, ISKWRK, SKWRK, IENRFL, IEION, PEDEP)
C  Returns an integer flag to indicate if reaction IR is ion-energy-
C  dependent, and if so, formulation-specific parameters.
C
C  INPUT
C  IR        - Integer scalar, reaction index;
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  IENRFL    - Integer scalar,
C              0, reaction IR does not have an ion-energy dependence
C              1, reaction IR does have an ion-energy dependence
C  IEION     - Integer scalar, species index of the ion on which
C              reaction is dependent;
C  PEDEP(*)  - Real array, supplemental parameters for an
C              ion-energy-dependent reaction rate formulation;
C              dimension at least NEDPAR, the number of supplemental
C              rate parameters required.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), PEDEP(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      ISIEDP = CKLKUP (IR, ISKWRK(ISKWRK(IiIEDP)), ISKWRK(IiNEDP))
      IF (ISIEDP .GT. 0) THEN
         IENRFL = 1
         IEION = ISKWRK(ISKWRK(IiKEDP + ISIEDP - 1))
         I_PAR = ISKWRK(IrPEDP) + (ISIEDP-1)*NEDPAR - 1
         DO 5 L = 1, NEDPAR
            PEDEP(L) = RSKWRK(I_PAR + L)
    5    CONTINUE
      ELSE
         IENRFL = 0
         IEION = 0
      ENDIF
C
C     end of SUBROUTINE SKIENR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKINDX (ISKWRK, NELM, KKGAS, KKSUR, KKBULK, KKTOT,
     1                   NNPHAS, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2                   NLBULK, IISUR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKINDX (ISKWRK, NELM, KKGAS, KKSUR, KKBULK, KKTOT,
C                     NNPHAS, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
C                     NLBULK, IISUR)
C  Returns a group of indices defining the size of the surface
C  reaction mechanism.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C
C  OUTPU
C  NELM      - Integer scalar, total element count.
C  KKGAS     - Integer scalar, total gas-phase species count.
C  KKSUR     - Integer scalar, total surface species count.
C  KKBULK    - Integer scalar, total bulk species count.
C  KKTOT     - Integer scalar, total species count (KKGAS+KKSUR+KKBULK).
C  NNPHAS    - Integer scalar, total phase count (gas + sites + bulks).
C  NNSURF    - Integer scalar, total surface phase count.
C  NFSURF    - Integer scalar, phase index of the first surface phase.
C  NLSURF    - Integer scalar, phase index of the last surface phase.
C  NNBULK    - Integer scalar, total bulk phase count.
C  NFBULK    - Integer scalar, phase index of the first bulk phase.
C  NLBULK    - Integer scalar, phase index of the last bulk phase.
C  IISUR     - Integer scalar, total surface reaction count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*)
C
      NELM  = NELEM
      KKGAS = NKKGAS
      KKSUR = ISKWRK(IiKSUR)
      KKBULK= ISKWRK(IiKBLK)
      KKTOT = ISKWRK(IiKTOT)
      NNPHAS= ISKWRK(IiNPHA)
      NNSURF= ISKWRK(IiNSUR)
      NFSURF= ISKWRK(IiFSUR)
      NLSURF= ISKWRK(IiLSUR)
      NNBULK= ISKWRK(IiNBLK)
      NFBULK= ISKWRK(IiFBLK)
      NLBULK= ISKWRK(IiLBLK)
      IISUR = ISKWRK(IiNIIS)
C
C     end of SUBROUTINE SKINDX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKINIT (LENISK, LENRSK, LENCSK, LINSK, LOUT,
     1                   ISKWRK, RSKWRK, CSKWRK, IFLAG)
C
C  START PROLOGUE
C
C  SUBROUTINE SKINIT (LENISK, LENRSK, LENCSK, LINSK, LOUT,
C                     ISKWRK, RSKWRK, CSKWRK, IFLAG)
C  Reads the surface linkfile and creates internal work arrays ISKWRK,
C  RSKWRK, and CSKWRK.  SKINIT must be called before any other Surface
C  Chemkin subroutine can be used, as the work arrays must be available
C  as their input.
C
C  INPUT
C  LENISK     - Integer scalar, length of the integer array ISKWRK.
C  LENRSK     - Integer scalar, length of the real array RSKWRK.
C  LENCSK     - Integer scalar, length of the character string array
C               CSKWRK.
C  LINSK      - Integer scalar, linkfile input file unit number.
C  LOUT       - Integer scalar, formatted output file unit number.
C
C  OUTPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C  CSKWRK(*) - Character string workspace array; dimension at least
C              LENCSK.
C  IFLAG     - Integer scalar to indicate successful reading of
C              linkfile; IFLAG>0 is an error type.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (RU=8.314510D7, RUC=RU/4.184D7, PA=1.01325D6)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C      PARAMETER (RU=8.314510E7, RUC=RU/4.184E7, PA=1.01325E6)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /SKCONS/ PREC, FILVER, KERR
      COMMON /MACHN/ SMALL, BIG, EXPARG
C
      DIMENSION ISKWRK(*), RSKWRK(*)
      CHARACTER*16 CSKWRK(*), FILVER, PRVERS, PREC, IFMT, RFMT, CFMT,
     1             SKVERS, SKPREC, SKDATE
      PARAMETER (IFMT='(10I12)', CFMT='(8A16)', RFMT='(1P,5E24.16)')
      LOGICAL IOK,ROK,COK,KERR,LBIN
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C*****precision > double
      SKPREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C      SKPREC = 'SINGLE'
C*****END precision > single
C
C----------------------------------------------------------------------C
C     Data about the machine dependent constants is carried in         C
C     COMMON/MACHN/SMALL,BIG,EXPARG
C
C      THIS STATEMENT WILL NOT COMPILE, MACHINE-DEPENDENT CONSTANTS
C
C*****exponent range > +/-30
C      SMALL = 1.0E-30
C      BIG   = 1.0E+30
C*****END exponent range > +/-30
C*****exponent range > +/-300
      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
C*****END exponent range > +/-300
C
      EXPARG = LOG(BIG)
C
C----------------------------------------------------------------------C
C
      SKVERS = '7.9'
      SKDATE = '98/04/08'
      WRITE (LOUT, '(/A, /1X,A, A, A, A, /A, /A, //)')
     1   ' SKLIB: CHEMKIN-III SURFACE KINETICS LIBRARY,',
     2   SKPREC(1:CKLSCH(SKPREC)), ' PRECISION Vers. ',
     3   SKVERS(1:CKLSCH(SKVERS)+1), SKDATE,
     4   ' Copyright 1995, Sandia Corporation.',
     5' The U.S. Government retains a limited license in this software.'
C
      CALL SKLEN (LINSK, LOUT, LENI, LENR, LENC, IFLAG)
      IF (IFLAG .GT. 0) RETURN
C
      IOK = (LENISK .GE. LENI)
      ROK = (LENRSK .GE. LENR)
      COK = (LENCSK .GE. LENC)
      IF (.NOT.IOK) WRITE (LOUT,300) LENI
      IF (.NOT.ROK) WRITE (LOUT,350) LENR
      IF (.NOT.COK) WRITE (LOUT,375) LENC
      IF (.NOT.IOK .OR. .NOT.ROK .OR. .NOT.COK) THEN
         IFLAG = 20
         RETURN
      ENDIF
C
      MORE = 0
C
C*****linkfile (surface) > binary
C      LBIN = .TRUE.
C*****END linkfile (surface) > binary
C*****linkfile (surface) > ascii
      LBIN = .FALSE.
C*****END linkfile (surface) > ascii
C
      IF (LBIN) THEN
         NREC = 1
         READ (LINSK, ERR=100) FILVER
         NREC = 2
         READ (LINSK, ERR=100) PRVERS
         NREC = 3
         READ (LINSK, ERR=100) PREC
         NREC = 4
         READ (LINSK, ERR=100) KERR
         NREC = 5
         READ (LINSK, ERR=100) LENI, LENR, LENC
         NREC = 6
         READ (LINSK, ERR=100) MAXSPR, MAXTP, NCP, NSPAR, NSCOV,
     1                         NEDPAR, NYPAR, MAXORD
         NREC = 7
         READ (LINSK, ERR=100) NELEM, NKKGAS, NKKSUR, NKKBLK, NKKTOT,
     1                         NPHASE, NFSUR, NLSUR, NNSUR, NFBLK,
     2                         NLBLK, NNBLK, NIISUR, NIICOV, NIIREV,
     3                         NIISTK, NIICON, NIIBHM, NIIRNU, NIIORD,
     4                         NIIEDP, NIIYLD, NKKION, KEL, MORE
         NREC = 8
         READ (LINSK, ERR=100) SKMIN
      ELSE
         NREC = 1
         READ (LINSK, '(A16)', ERR=100) FILVER
         NREC = 2
         READ (LINSK, '(A16)', ERR=100) PRVERS
         NREC = 3
         READ (LINSK, '(A16)', ERR=100) PREC
         NREC = 4
         READ (LINSK, '(L8)', ERR=100) KERR
         NREC = 5
         READ (LINSK, '(3I12)', ERR=100) LENI, LENR, LENC
         NREC = 6
         READ (LINSK, IFMT, ERR=100) MAXSPR, MAXTP, NCP, NSPAR, NSCOV,
     1                               NEDPAR, NYPAR, MAXORD
         NREC = 7
         READ (LINSK, IFMT, ERR=100)
     1         NELEM, NKKGAS, NKKSUR, NKKBLK, NKKTOT,
     1         NPHASE, NFSUR, NLSUR, NNSUR, NFBLK,
     2         NLBLK, NNBLK, NIISUR, NIICOV, NIIREV,
     3         NIISTK, NIICON, NIIBHM, NIIRNU, NIIORD,
     4         NIIEDP, NIIYLD, NKKION, KEL, MORE
         NREC = 8
         READ (LINSK, RFMT, ERR=100) SKMIN
      ENDIF
C
      NCP1  = NCP+1
      NCP2  = NCP+2
      NCP2T = NCP2*(MAXTP-1)
C
C     ISKWRK pointers to integer variables, array locations
C
      IiLENI = 1
      IiLENR = IiLENI + 1
      IiLENC = IiLENR + 1
      IiKSUR = IiLENC + 1
      IiKBLK = IiKSUR + 1
      IiKTOT = IiKBLK + 1
      IiNPHA = IiKTOT + 1
      IiFSUR = IiNPHA + 1
      IiLSUR = IiFSUR + 1
      IiNSUR = IiLSUR + 1
      IiFBLK = IiNSUR + 1
      IiLBLK = IiFBLK + 1
      IiNBLK = IiLBLK + 1
      IiNIIS = IiNBLK + 1
      IiNCOV = IiNIIS + 1
      IiNREV = IiNCOV + 1
      IiNSTK = IiNREV + 1
      IiNCON = IiNSTK + 1
      IiNBHM = IiNCON + 1
      IiNRNU = IiNBHM + 1
      IiNORD = IiNRNU + 1
      IiKION = IiNORD + 1
      IiKTFL = IiKION + 1
      IiNEDP = IiKTFL + 1
      IiKEDP = IiNEDP + 1
      IiIEDP = IiKEDP + 1
      IiNYLD = IiIEDP + 1
      IiIYLD = IiNYLD + 1
      IiYION = IiIYLD + 1
      IiKYLD = IiYION + 1
      IiIONS = IiKYLD + 1
      IiELEC = IiIONS + 1
      IiPKST = IiELEC + 1
      IiPKND = IiPKST + 1
      IiPTOT = IiPKND + 1
      IiKPHS = IiPTOT + 1
      IiKCHG = IiKPHS + 1
      IiKCMP = IiKCHG + 1
      IiNSCV = IiKCMP + 1
      IiKNT  = IiNSCV + 1
      IiNRPP = IiKNT  + 1
      IiNREA = IiNRPP + 1
      IiNUNK = IiNREA + 1
      IiNU   = IiNUNK + 1
      IiNSUM = IiNU   + 1
      IiICOV = IiNSUM + 1
      IiKCOV = IiICOV + 1
      IiIREV = IiKCOV + 1
      IiISTK = IiIREV + 1
      IiMSTK = IiISTK + 1
      IiIBHM = IiMSTK + 1
      IiKBHM = IiIBHM + 1
      IiIRNU = IiKBHM + 1
      IiIORD = IiIRNU + 1
      IiKORD = IiIORD + 1
C
C     ISKWRK pointers for real variables, arrays
C
      IrSKMN = IiKORD + 1
      IrPATM = IrSKMN + 1
      IrRU   = IrPATM + 1
      IrRUC  = IrRU   + 1
      IrSDEN = IrRUC  + 1
      IrKTMP = IrSDEN + 1
      IrKTHM = IrKTMP + 1
      IrKDEN = IrKTHM + 1
      IrAWT  = IrKDEN + 1
      IrKWT  = IrAWT  + 1
      IrPAR  = IrKWT  + 1
      IrKCOV = IrPAR  + 1
      IrRPAR = IrKCOV + 1
      IrEQ   = IrRPAR + 1
      IrRNU  = IrEQ   + 1
      IrNCF  = IrRNU  + 1
      IrKORD = IrNCF  + 1
      IrPEDP = IrKORD + 1
      IrENGI = IrPEDP + 1
      IrPYLD = IrENGI + 1
      IrYNCF = IrPYLD + 1
      IrKFT  = IrYNCF + 1
      IrKRT  = IrKFT  + 1
      IrKT1  = IrKRT  + 1
      IrKT2  = IrKT1  + 1
      IrPT1  = IrKT2  + 1
      IrIT1  = IrPT1  + 1
      IrIT2  = IrIT1  + 1
      IrIT3  = IrIT2  + 1
C
C     ISKWRK pointers to character variables, arrays
C
      IcENAM = IrIT3  + 1
      IcKNAM = IcENAM + 1
      IcMNAM = IcKNAM + 1
      IcPNAM = IcMNAM + 1
C
      NPNTS = IcPNAM
C
      ISKWRK(IiLENI) = LENI
      ISKWRK(IiLENR) = LENR
      ISKWRK(IiLENC) = LENC
C! number of site species
      ISKWRK(IiKSUR) = NKKSUR
C! number of bulk species
      ISKWRK(IiKBLK) = NKKBLK
C! total number of species
      ISKWRK(IiKTOT) = NKKTOT
C! number of phases
      ISKWRK(IiNPHA) = NPHASE
C! first site phase
      ISKWRK(IiFSUR) = NFSUR
C! last site phase
      ISKWRK(IiLSUR) = NLSUR
C! number of site phases
      ISKWRK(IiNSUR) = NNSUR
C! first bulk phase
      ISKWRK(IiFBLK) = NFBLK
C! last bulk phase
      ISKWRK(IiLBLK) = NLBLK
C! number of bulk phases
      ISKWRK(IiNBLK) = NNBLK
C! number of surface reactions
      ISKWRK(IiNIIS) = NIISUR
C! number of reactions with coverage parameters
      ISKWRK(IiNCOV) = NIICOV
C! number of reactions with reverse parameters
      ISKWRK(IiNREV) = NIIREV
C! number of reactions with sticking coefficients
      ISKWRK(IiNSTK) = NIISTK
C! number of reactions with non-conservation of sites
      ISKWRK(IiNCON) = NIICON
C! number of reactions with Bohm parameters
      ISKWRK(IiNBHM) = NIIBHM
C! number of reactions with real stoichiometry
      ISKWRK(IiNRNU) = NIIRNU
C! number of reactions with change-of-order
      ISKWRK(IiNORD) = NIIORD
C! number of reactions with ion-energy dependence
      ISKWRK(IiNEDP) = NIIEDP
C! number of reactions with #-modified stoichiometry
      ISKWRK(IiNYLD) = NIIYLD
C! number gas-phase ions
      ISKWRK(IiKION) = NKKION
C! species number for the electron
      ISKWRK(IiELEC) = KEL
C! starting species location of phases
      ISKWRK(IiPKST) = NPNTS + 1
C! ending species location of phases
      ISKWRK(IiPKND) = ISKWRK(IiPKST) + NPHASE
C! number of species each phase
      ISKWRK(IiPTOT) = ISKWRK(IiPKND) + NPHASE
C! species phases
      ISKWRK(IiKPHS) = ISKWRK(IiPTOT) + NPHASE
C! species charges
      ISKWRK(IiKCHG) = ISKWRK(IiKPHS) + NKKTOT
C! species elemental composition
      ISKWRK(IiKCMP) = ISKWRK(IiKCHG) + NKKTOT
C! site coverage of species
      ISKWRK(IiNSCV) = ISKWRK(IiKCMP) + NKKTOT*NELEM
C! # of temperatures for fit
      ISKWRK(IiKNT)  = ISKWRK(IiNSCV) + NKKTOT
C! num of species in reactions
      ISKWRK(IiNRPP) = ISKWRK(IiKNT) + NKKTOT
C! num of reactants in reactions
      ISKWRK(IiNREA) = ISKWRK(IiNRPP) + NIISUR
C! spec. numbers in reactions
      ISKWRK(IiNUNK) = ISKWRK(IiNREA) + NIISUR
C! species coefficients
      ISKWRK(IiNU)   = ISKWRK(IiNUNK) + NIISUR*MAXSPR
C! total gas-phase coefficients
      ISKWRK(IiNSUM) = ISKWRK(IiNU) + NIISUR*MAXSPR
C! reaction num. for coverage parameters
      ISKWRK(IiICOV)  = ISKWRK(IiNSUM) + NIISUR
C! species num. for coverage parameters
      ISKWRK(IiKCOV)  = ISKWRK(IiICOV)  + NIICOV
C!  " with reverse parameters
      ISKWRK(IiIREV)  = ISKWRK(IiKCOV)  + NIICOV
C!  " with sticking parameters
      ISKWRK(IiISTK)  = ISKWRK(IiIREV)  + NIIREV
C!  Motz-Wise correction flag for sticking
      ISKWRK(IiMSTK)  = ISKWRK(IiISTK) + NIISTK
C!  reactions with Bohm velocity conditions
      ISKWRK(IiIBHM)  = ISKWRK(IiMSTK)  + NIISTK
C! species numbers used in Bohm velocity formulation
      ISKWRK(IiKBHM)   = ISKWRK(IiIBHM)  + NIIBHM
C! species with real stoichometric coeff's
      ISKWRK(IiIRNU)  = ISKWRK(IiKBHM)  + NIIBHM
C! species with change of order
      ISKWRK(IiIORD)  = ISKWRK(IiIRNU) + NIIRNU
C!  change of order species
      ISKWRK(IiKORD)  = ISKWRK(IiIORD) + NIIORD
C!  ionic species indices for reactions with energy dependence
      ISKWRK(IiKEDP)  = ISKWRK(IiKORD) + NIIORD*MAXORD
C!  reaction numbers for reactions with energy dependence
      ISKWRK(IiIEDP)  = ISKWRK(IiKEDP) + NIIEDP
C!  reaction numbers for the #-modified reactions
      ISKWRK(IiIYLD)  = ISKWRK(IiIEDP) + NIIEDP
C!  the gas-phase ion for the #-modified reactions
      ISKWRK(IiYION)  = ISKWRK(IiIYLD) + NIIYLD
C!  integer flags for the #-modified species
      ISKWRK(IiKYLD)  = ISKWRK(IiYION) + NIIYLD
C!  indices of the gas-phase ions
      ISKWRK(IiIONS)  = ISKWRK(IiKYLD) + NIIYLD*MAXSPR
C!  array indicating which temperature (in the temperature array)
C   that each gas-phase species will use
      ISKWRK(IiKTFL)  = ISKWRK(IiIONS) + NKKION
C! total size req. for ISKWRK
      IiTOT   = ISKWRK(IiKTFL) + NKKGAS
      ISKWRK(IiTOT) = MORE
C
C! minimum difference allowed for balancing stoichiometry,
C  sites, etc.
      ISKWRK(IrSKMN) = 1
      RSKWRK(ISKWRK(IrSKMN)) = SKMIN
C! pressure 1 atm
      ISKWRK(IrPATM) = ISKWRK(IrSKMN) + 1
      RSKWRK(ISKWRK(IrPATM)) = PA
C
      ISKWRK(IrRU)   = ISKWRK(IrPATM) + 1
      RSKWRK(ISKWRK(IrRU)) = RU
C
      ISKWRK(IrRUC)  = ISKWRK(IrRU) + 1
      RSKWRK(ISKWRK(IrRUC)) = RU / 4.184E7
C! densities for phases
      ISKWRK(IrSDEN) = ISKWRK(IrRUC) + 1
C! temperatures for fit
      ISKWRK(IrKTMP)  = ISKWRK(IrSDEN) + NPHASE
C! thermodynamic properties
      ISKWRK(IrKTHM)  = ISKWRK(IrKTMP) + MAXTP *NKKTOT
C! densities for the species
      ISKWRK(IrKDEN)  = ISKWRK(IrKTHM) + NCP2T*NKKTOT
C! atomic weights of elements
      ISKWRK(IrAWT)   = ISKWRK(IrKDEN) + NKKTOT
C! molecular weights of species
      ISKWRK(IrKWT)   = ISKWRK(IrAWT)  + NELEM
C! Arrhenius parameters
      ISKWRK(IrPAR)   = ISKWRK(IrKWT)  + NKKTOT
C! coverage parameters
      ISKWRK(IrKCOV)   = ISKWRK(IrPAR)  + NIISUR*(NSPAR+1)
C! reverse parameters
      ISKWRK(IrRPAR)   = ISKWRK(IrKCOV)  + NIICOV*NSCOV
C! multiplicative constant in surface equilibrium constants
      ISKWRK(IrEQ)    = ISKWRK(IrRPAR) + NIIREV*(NSPAR+1)
C! real stoichometric coefficients
      ISKWRK(IrRNU)   = ISKWRK(IrEQ) + NIISUR
C! real phase (site) balances
      ISKWRK(IrNCF)   = ISKWRK(IrRNU) + NIIRNU*MAXSPR
C! yield-modify parameters
      ISKWRK(IrPYLD)  = ISKWRK(IrNCF) + NIISUR*NPHASE
C! real phase (site) balances for the #-modified sites
      ISKWRK(IrYNCF)  = ISKWRK(IrPYLD) + NIIYLD*NYPAR
C! species order coefficients
      ISKWRK(IrKORD)   = ISKWRK(IrYNCF) + NIIYLD*NPHASE
C! parameters for ion-energy dependence reactions
      ISKWRK(IrPEDP)  = ISKWRK(IrKORD) + NIIORD*MAXORD
C! ion energy array used internally
      ISKWRK(IrENGI)  = ISKWRK(IrPEDP) + NIIEDP*NEDPAR
C
      ISKWRK(IrKFT)   = ISKWRK(IrENGI) + NKKGAS
C!
      ISKWRK(IrKRT)   = ISKWRK(IrKFT) + NIISUR
C! dummy species work space
      ISKWRK(IrKT1)   = ISKWRK(IrKRT) + NIISUR
C! dummy species work space
      ISKWRK(IrKT2)   = ISKWRK(IrKT1)  + NKKTOT
C! dummy phase work space
      ISKWRK(IrPT1)   = ISKWRK(IrKT2)  + NKKTOT
C! dummy reaction work space
      ISKWRK(IrIT1)   = ISKWRK(IrPT1)  + NPHASE
C! dummy reaction work space
      ISKWRK(IrIT2)   = ISKWRK(IrIT1)  + NIISUR
C! dummy reaction work space
      ISKWRK(IrIT3)   = ISKWRK(IrIT2)  + NIISUR
C! total size req. for RSKWRK
      IrTOT   = ISKWRK(IrIT3)  + NIISUR-1
      NTOT = 4 + NPHASE + MAXTP*NKKTOT + NCP2T*NKKTOT
     1         + NKKTOT + NELEM + NKKTOT + NIISUR*(NSPAR+1)
     1         + NIICOV*NSCOV + NIIREV*(NSPAR+1) + NIISUR
     1         + NIIRNU*MAXSPR + NIISUR*NPHASE + NIIYLD*NYPAR
     2         + NIIYLD*NPHASE + NIIORD*MAXORD + NIIEDP*NEDPAR
     3         + NKKGAS + NIISUR + NIISUR + NKKTOT + NKKTOT
     4         + NPHASE + NIISUR + NIISUR + NIISUR
C
C! element names
      ISKWRK(IcENAM) = 1
C! species names
      ISKWRK(IcKNAM) = ISKWRK(IcENAM) + NELEM
C! material name
      ISKWRK(IcMNAM) = ISKWRK(IcKNAM) + NKKTOT
C! phase names
      ISKWRK(IcPNAM) = ISKWRK(IcMNAM) + 1
C! total size req. for CSKWRK
      IcTOT  = ISKWRK(IcPNAM) + NPHASE - 1
C
C     RSKWRK(ISKWRK(IrSKMN)) = SKMIN
C! dynes per sq. cm.
      RSKWRK(ISKWRK(IrPATM)) = PA
C! calories/ mole K
      RSKWRK(ISKWRK(IrRUC))  = RUC
C! ergs/ mole K
      RSKWRK(ISKWRK(IrRU))   = RU
C
      IF (LBIN) THEN
         NREC = 9
         READ (LINSK, ERR=100) CSKWRK(ISKWRK(IcMNAM))
         NREC = 10
         READ (LINSK, ERR=100)
     1   ( CSKWRK(ISKWRK(IcENAM) + M - 1), M = 1, NELEM)
         NREC = 11
         READ (LINSK, ERR=100)
     1   ( RSKWRK(ISKWRK(IrAWT) + M - 1), M = 1, NELEM)
         NREC = 12
         READ (LINSK, ERR=100)
     1   ( CSKWRK(ISKWRK(IcKNAM) + K - 1), K = 1, NKKTOT)
         NREC = 13
         READ (LINSK, ERR=100)
     1   ( RSKWRK(ISKWRK(IrKWT) + K - 1), K = 1, NKKTOT)
         NREC = 14
         READ (LINSK, ERR=100)
     1   ( (ISKWRK(ISKWRK(IiKCMP)+(K-1)*NELEM+M-1), M=1,NELEM),
     2       K=1,NKKTOT)
         NREC = 15
         READ (LINSK, ERR=100)
     1   ( ISKWRK(ISKWRK(IiKCHG) + K - 1), K = 1, NKKTOT)
         NREC = 16
         READ (LINSK, ERR=100)
     1   ( ISKWRK(ISKWRK(IiKNT) + K - 1), K = 1, NKKTOT)
         NREC = 17
         READ (LINSK, ERR=100)
     1   ( ISKWRK(ISKWRK(IiKPHS) + K - 1), K = 1, NKKTOT)
         NREC = 18
         READ (LINSK, ERR=100)
     1   ( ISKWRK(ISKWRK(IiNSCV) + K - 1), K = 1, NKKTOT)
         NREC = 19
         READ (LINSK, ERR=100)
     1   ( RSKWRK(ISKWRK(IrKDEN) + K - 1), K = 1, NKKTOT)
         NREC = 20
         READ (LINSK, ERR=100)
     1   ( (RSKWRK(ISKWRK(IrKTMP)+(K-1)*MAXTP + L - 1), L=1,MAXTP),
     2      K = 1, NKKTOT)
         NREC = 21
         READ (LINSK, ERR=100)
     1   (( (RSKWRK(ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T + N - 1),
     2      N=1,NCP2), L=1,(MAXTP-1)), K = 1, NKKTOT)
C
         IF (NKKION .GT. 0) THEN
            NREC = 22
            READ (LINSK, ERR=100)
            NREC = 23
            READ (LINSK, ERR=100)
     1      ( ISKWRK(ISKWRK(IiIONS) + K - 1), K = 1, NKKION)
         ENDIF
C
         NREC = 24
         READ (LINSK, ERR=100)
     1   ( CSKWRK(ISKWRK(IcPNAM) + N - 1), N = 1, NPHASE)
         NREC = 25
         READ (LINSK, ERR=100)
     1   ( ISKWRK(ISKWRK(IiPKST) + N - 1), N = 1, NPHASE)
         NREC = 26
         READ (LINSK, ERR=100)
     1   ( ISKWRK(ISKWRK(IiPKND) + N - 1), N = 1, NPHASE)
         NREC = 27
         READ (LINSK, ERR=100)
     1   ( ISKWRK(ISKWRK(IiPTOT) + N - 1), N = 1, NPHASE)
         NREC = 28
         READ (LINSK, ERR=100)
     1   ( RSKWRK(ISKWRK(IrSDEN) + N - 1), N = 1, NPHASE)
C
      ELSE
         NREC = 9
         READ (LINSK, CFMT, ERR=100) CSKWRK(ISKWRK(IcMNAM))
         NREC = 10
         READ (LINSK, CFMT, ERR=100)
     1   ( CSKWRK(ISKWRK(IcENAM) + M - 1), M = 1, NELEM)
         NREC = 11
         READ (LINSK, RFMT, ERR=100)
     1   ( RSKWRK(ISKWRK(IrAWT) + M - 1), M = 1, NELEM)
         NREC = 12
         READ (LINSK, CFMT, ERR=100)
     1   ( CSKWRK(ISKWRK(IcKNAM) + K - 1), K = 1, NKKTOT)
         NREC = 13
         READ (LINSK, RFMT, ERR=100)
     1   ( RSKWRK(ISKWRK(IrKWT) + K - 1), K = 1, NKKTOT)
         NREC = 14
         READ (LINSK, IFMT, ERR=100)
     1   ( (ISKWRK(ISKWRK(IiKCMP)+(K-1)*NELEM+M-1), M=1,NELEM),
     2      K=1,NKKTOT)
         NREC = 15
         READ (LINSK, IFMT, ERR=100)
     1   ( ISKWRK(ISKWRK(IiKCHG) + K - 1), K = 1, NKKTOT)
         NREC = 16
         READ (LINSK, IFMT, ERR=100)
     1   ( ISKWRK(ISKWRK(IiKNT) + K - 1), K = 1, NKKTOT)
         NREC = 17
         READ (LINSK, IFMT, ERR=100)
     1   ( ISKWRK(ISKWRK(IiKPHS) + K - 1), K = 1, NKKTOT)
         NREC = 18
         READ (LINSK, IFMT, ERR=100)
     1   ( ISKWRK(ISKWRK(IiNSCV) + K - 1), K = 1, NKKTOT)
         NREC = 19
         READ (LINSK, RFMT, ERR=100)
     1   ( RSKWRK(ISKWRK(IrKDEN) + K - 1), K = 1, NKKTOT)
         NREC = 20
         READ (LINSK, RFMT, ERR=100)
     1   ( (RSKWRK(ISKWRK(IrKTMP)+(K-1)*MAXTP + L - 1), L=1,MAXTP),
     2      K = 1, NKKTOT)
         NREC = 21
         READ (LINSK, RFMT, ERR=100)
     1   (( (RSKWRK(ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T + N - 1),
     2      N=1,NCP2), L=1,(MAXTP-1)), K = 1, NKKTOT)
C
         IF (NKKION .GT. 0) THEN
            NREC = 22
            READ (LINSK, IFMT, ERR=100)
            NREC = 23
            READ (LINSK, IFMT, ERR=100)
     1      ( ISKWRK(ISKWRK(IiIONS) + K - 1), K = 1, NKKION)
         ENDIF
C
         NREC = 24
         READ (LINSK, CFMT, ERR=100)
     1   ( CSKWRK(ISKWRK(IcPNAM) + N - 1), N = 1, NPHASE)
         NREC = 25
         READ (LINSK, IFMT, ERR=100)
     1   ( ISKWRK(ISKWRK(IiPKST) + N - 1), N = 1, NPHASE)
         NREC = 26
         READ (LINSK, IFMT, ERR=100)
     1   ( ISKWRK(ISKWRK(IiPKND) + N - 1), N = 1, NPHASE)
         NREC = 27
         READ (LINSK, IFMT, ERR=100)
     1   ( ISKWRK(ISKWRK(IiPTOT) + N - 1), N = 1, NPHASE)
         NREC = 28
         READ (LINSK, RFMT, ERR=100)
     1   ( RSKWRK(ISKWRK(IrSDEN) + N - 1), N = 1, NPHASE)
      ENDIF
C
      I_KTFL = ISKWRK(IiKTFL)
      NKM1 = NKKGAS - 1
      DO 4 K = 0, NKM1
C
C        ASSIGN THE DEFAULT VALUE OF 1 TO EACH ELEMENT OF KTFL FOR
C        THE GAS-PHASE SPECIES. THUS, UNLESS OTHERWISE SPECIFIED
C        BY SETTING THE VALUES IN A CALL TO SKKTFL, EACH GAS-PHASE
C        SPECIES WILL USE THE FIRST ELEMENT IN THE TEMPERATURE
C        ARRAY, I.E., THE THERMAL TEMPERATURE
C
         ISKWRK(I_KTFL + K) = 1
 4    CONTINUE
c
c        ASSIGN THE DEFAULT VALUE OF 0.0 FOR THE ION ENERGIES.
C        THESE WILL ONLY BE SET DIFFERENTLY IF SKRPAR IS CALLED.
C
      I_ENGR = ISKWRK(IrENGI)
      NKM1 = NKKGAS - 1
      DO 6 K = 0, NKM1
         RSKWRK(I_ENGR + K) = 0.0
 6    CONTINUE
C
      IF (NIISUR .LE. 0) THEN
         IF (MORE .LE. 0) REWIND (LINSK)
         RETURN
      ENDIF
C
      NRPP = ISKWRK(IiNRPP)
      NREA = ISKWRK(IiNREA)
      NIEQ = ISKWRK(IrEQ)
      IF (LBIN) THEN
         NREC = 29
         READ (LINSK, ERR=100) (ISKWRK(NRPP + I - 1), I=1,NIISUR)
         NREC = 30
         READ (LINSK, ERR=100) (ISKWRK(NREA + I - 1), I=1,NIISUR)
         NREC = 31
         READ (LINSK, ERR=100)
     1   ( (ISKWRK(ISKWRK(IiNU)+(I-1)*MAXSPR+N-1),
     2      ISKWRK(ISKWRK(IiNUNK)+(I-1)*MAXSPR+N-1), N=1,MAXSPR),
     3     I = 1, NIISUR)
         NREC = 32
         READ (LINSK, ERR=100)
     1   ( ISKWRK(ISKWRK(IiNSUM) + I - 1), I = 1, NIISUR)
         NREC = 33
         READ (LINSK, ERR=100)
     1   ( (RSKWRK(ISKWRK(IrPAR)+(I-1)*(NSPAR+1)+N-1), N=1,NSPAR),
     2      I = 1, NIISUR)
         NREC = 34
         READ (LINSK, ERR=100) (RSKWRK(NIEQ + I - 1), I=1,NIISUR)
         NREC = 35
         READ (LINSK, ERR=100)
     1   ( (RSKWRK(ISKWRK(IrNCF) + (I-1)*NPHASE+N-1), N=1,NPHASE),
     2      I = 1, NIISUR)
      ELSE
         NREC = 29
         READ (LINSK, IFMT, ERR=100) (ISKWRK(NRPP+I-1), I=1,NIISUR)
         NREC = 30
         READ (LINSK, IFMT, ERR=100) (ISKWRK(NREA+I-1), I=1,NIISUR)
         NREC = 31
         READ (LINSK, IFMT, ERR=100)
     1   ( (ISKWRK(ISKWRK(IiNU)+(I-1)*MAXSPR+N-1),
     2      ISKWRK(ISKWRK(IiNUNK)+(I-1)*MAXSPR+N-1), N=1,MAXSPR),
     3     I = 1, NIISUR)
         NREC = 32
         READ (LINSK, IFMT, ERR=100)
     1   ( ISKWRK(ISKWRK(IiNSUM) + I - 1), I = 1, NIISUR)
         NREC = 33
         READ (LINSK, RFMT, ERR=100)
     1   ( (RSKWRK(ISKWRK(IrPAR)+(I-1)*(NSPAR+1)+N-1), N=1,NSPAR),
     2      I = 1, NIISUR)
         NREC = 34
         READ (LINSK, RFMT, ERR=100) (RSKWRK(NIEQ+I-1), I=1,NIISUR)
         NREC = 35
         READ (LINSK, RFMT, ERR=100)
     1   ( (RSKWRK(ISKWRK(IrNCF) + (I-1)*NPHASE+N-1), N=1,NPHASE),
     2      I = 1, NIISUR)
      ENDIF
C
      IF (NIICOV .GT. 0) THEN
         INUM = ISKWRK(IiICOV)
         KNUM = ISKWRK(IiKCOV)
         IF (LBIN) THEN
            NREC = 36
            READ (LINSK, ERR=100)
            NREC = 37
            READ (LINSK, ERR=100) (ISKWRK(INUM + N - 1), N=1,NIICOV)
            NREC = 38
            READ (LINSK, ERR=100) (ISKWRK(KNUM + N - 1), N=1,NIICOV)
            NREC = 39
            READ (LINSK, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrKCOV)+(N-1)*NSCOV+L-1),L=1,NSCOV),
     2         N=1,NIICOV)
         ELSE
            NREC = 36
            READ (LINSK, IFMT, ERR=100)
            NREC = 37
            READ (LINSK, IFMT, ERR=100) (ISKWRK(INUM+N-1),N=1,NIICOV)
            NREC = 38
            READ (LINSK, IFMT, ERR=100) (ISKWRK(KNUM+N-1),N=1,NIICOV)
            NREC = 39
            READ (LINSK, RFMT, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrKCOV)+(N-1)*NSCOV+L-1),L=1,NSCOV),
     2         N=1,NIICOV)
         ENDIF
      ENDIF
C
      IF (NIIREV .GT. 0) THEN
         INUM = ISKWRK(IiIREV)
         IF (LBIN) THEN
            NREC = 40
            READ (LINSK, ERR=100)
            NREC = 41
            READ (LINSK, ERR=100) (ISKWRK(INUM + N - 1), N=1,NIIREV)
            NREC = 42
            READ (LINSK, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrRPAR)+(N-1)*(NSPAR+1)+L-1),L=1,NSPAR),
     2         N=1,NIIREV)
         ELSE
            NREC = 40
            READ (LINSK, IFMT, ERR=100)
            NREC = 41
            READ (LINSK, IFMT, ERR=100) (ISKWRK(INUM+N-1), N=1,NIIREV)
            NREC = 42
            READ (LINSK, RFMT, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrRPAR)+(N-1)*(NSPAR+1)+L-1),L=1,NSPAR),
     2         N=1,NIIREV)
         ENDIF
      ENDIF
C
      IF (NIISTK .GT. 0) THEN
         INUM = ISKWRK(IiISTK)
         MNUM = ISKWRK(IiMSTK)
         IF (LBIN) THEN
            NREC = 43
            READ (LINSK, ERR=100)
            NREC = 44
            READ (LINSK, ERR=100) (ISKWRK(INUM + N - 1), N=1,NIISTK)
            NREC = 45
            READ (LINSK, ERR=100) (ISKWRK(MNUM + N - 1), N=1,NIISTK)
         ELSE
            NREC = 43
            READ (LINSK, IFMT, ERR=100)
            NREC = 44
            READ (LINSK, IFMT, ERR=100) (ISKWRK(INUM+N-1), N=1,NIISTK)
            NREC = 45
            READ (LINSK, IFMT, ERR=100) (ISKWRK(MNUM+N-1), N=1,NIISTK)
         ENDIF
      ENDIF
C
      IF (NIIBHM .GT. 0) THEN
         INUM = ISKWRK(IiIBHM)
         KNUM = ISKWRK(IiKBHM)
         IF (LBIN) THEN
            NREC = 46
            READ (LINSK, ERR=100)
            NREC = 47
            READ (LINSK, ERR=100) (ISKWRK(INUM + N - 1), N=1,NIIBHM)
            NREC = 48
            READ (LINSK, ERR=100) (ISKWRK(KNUM + N - 1), N=1,NIIBHM)
         ELSE
            NREC = 46
            READ (LINSK, IFMT, ERR=100)
            NREC = 47
            READ (LINSK, IFMT, ERR=100) (ISKWRK(INUM+N-1), N=1,NIIBHM)
            NREC = 48
            READ (LINSK, IFMT, ERR=100) (ISKWRK(KNUM+N-1), N=1,NIIBHM)
         ENDIF
      ENDIF
C
      IF (NIIRNU .GT. 0) THEN
         INUM = ISKWRK(IiIRNU)
         IF (LBIN) THEN
            NREC = 49
            READ (LINSK, ERR=100)
            NREC = 50
            READ (LINSK, ERR=100) (ISKWRK(INUM+N-1), N=1,NIIRNU)
            NREC = 51
            READ (LINSK, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrRNU) + (N-1)*MAXSPR + L-1),L=1,MAXSPR),
     2         N = 1, NIIRNU)
         ELSE
            NREC = 49
            READ (LINSK, IFMT, ERR=100)
            NREC = 50
            READ (LINSK, IFMT, ERR=100) (ISKWRK(INUM+N-1), N=1,NIIRNU)
            NREC = 51
            READ (LINSK, RFMT, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrRNU) + (N-1)*MAXSPR + L-1),L=1,MAXSPR),
     2         N = 1, NIIRNU)
         ENDIF
      ENDIF
C
      IF (NIIORD .GT. 0) THEN
         INUM = ISKWRK(IiIORD)
         IF (LBIN) THEN
            NREC = 52
            READ (LINSK, ERR=100)
            NREC = 53
            READ (LINSK, ERR=100) (ISKWRK(INUM + N - 1), N=1,NIIORD)
            NREC = 54
            READ (LINSK, ERR=100)
     1      ( (ISKWRK(ISKWRK(IiKORD) + (N-1)*MAXORD + L - 1),
     2          L=1,MAXORD), N = 1, NIIORD)
            NREC = 55
            READ (LINSK, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrKORD) + (N-1)*MAXORD + L - 1),
     2         L=1,MAXORD), N = 1, NIIORD)
         ELSE
            NREC = 52
            READ (LINSK, IFMT, ERR=100)
            NREC = 53
            READ (LINSK, IFMT, ERR=100) (ISKWRK(INUM+N-1), N=1,NIIORD)
            NREC = 54
            READ (LINSK, IFMT, ERR=100)
     1      ( (ISKWRK(ISKWRK(IiKORD) + (N-1)*MAXORD + L - 1),
     2          L=1,MAXORD), N = 1, NIIORD)
            NREC = 55
            READ (LINSK, RFMT, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrKORD) + (N-1)*MAXORD + L - 1),
     2         L=1,MAXORD), N = 1, NIIORD)
         ENDIF
      ENDIF
C
      IF (NIIEDP .GT. 0) THEN
         INUM = ISKWRK(IiIEDP)
         KNUM = ISKWRK(IiKEDP)
         IF (LBIN) THEN
            NREC = 56
            READ (LINSK, ERR=100)
            NREC = 57
            READ (LINSK, ERR=100) (ISKWRK(INUM + N - 1), N=1,NIIEDP)
            NREC = 58
            READ (LINSK, ERR=100) (ISKWRK(KNUM + N - 1), N=1,NIIEDP)
            NREC = 59
            READ (LINSK, ERR=100)
     1     ( (RSKWRK(ISKWRK(IrPEDP) + (L-1)*NEDPAR + L-1), L=1,NEDPAR),
     2        I=1,NIIEDP)
         ELSE
            NREC = 56
            READ (LINSK, IFMT, ERR=100)
            NREC = 57
            READ (LINSK, IFMT, ERR=100) (ISKWRK(INUM+N-1), N=1,NIIEDP)
            NREC = 58
            READ (LINSK, IFMT, ERR=100) (ISKWRK(KNUM+N-1), N=1,NIIEDP)
            NREC = 59
            READ (LINSK, RFMT, ERR=100)
     1     ( (RSKWRK(ISKWRK(IrPEDP) + (N-1)*NEDPAR + L-1), L=1,NEDPAR),
     2        N=1,NIIEDP)
         ENDIF
      ENDIF
C
      IF (NIIYLD .GT. 0) THEN
         INUM = ISKWRK(IiIYLD)
         IION = ISKWRK(IiYION)
         IF (LBIN) THEN
            NREC = 60
            READ (LINSK, ERR=100)
            NREC = 61
            READ (LINSK, ERR=100) (ISKWRK(INUM + N - 1), N=1,NIIYLD)
            NREC = 62
            READ (LINSK, ERR=100) (ISKWRK(IION + N - 1), N=1,NIIYLD)
            NREC = 63
            READ (LINSK, ERR=100)
     1      ( (ISKWRK(ISKWRK(IiKYLD) + (N-1)*MAXSPR + L-1), L=1,MAXSPR),
     2        N = 1, NIIYLD)
            NREC = 64
            READ (LINSK, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrPYLD) + (N-1)*NYPAR + L-1), L=1,NYPAR),
     2        N = 1, NIIYLD)
            NREC = 65
            READ (LINSK, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrYNCF) + (N-1)*NPHASE + L-1), L=1,NPHASE),
     2         N = 1, NIIYLD)
         ELSE
            NREC = 60
            READ (LINSK, IFMT, ERR=100)
            NREC = 61
            READ (LINSK, IFMT, ERR=100) (ISKWRK(INUM+N-1),N=1,NIIYLD)
            NREC = 62
            READ (LINSK, IFMT, ERR=100) (ISKWRK(IION+N-1),N=1,NIIYLD)
            NREC = 63
            READ (LINSK, IFMT, ERR=100)
     1      ( (ISKWRK(ISKWRK(IiKYLD) + (N-1)*MAXSPR + L-1), L=1,MAXSPR),
     2        N = 1, NIIYLD)
            NREC = 64
            READ (LINSK, RFMT, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrPYLD) + (N-1)*NYPAR + L-1), L=1,NYPAR),
     2        N = 1, NIIYLD)
            NREC = 65
            READ (LINSK, RFMT, ERR=100)
     1      ( (RSKWRK(ISKWRK(IrYNCF) + (N-1)*NPHASE + L-1), L=1,NPHASE),
     2         N = 1, NIIYLD)
         ENDIF
      ENDIF
C
      DO 10 I = 0, NIISUR - 1
C        PERTURBATION FACTOR
         I_PAR = ISKWRK(IrPAR) + I*(NSPAR+1) + NSPAR
         RSKWRK(I_PAR) = 1.0
   10 CONTINUE
C
      IF (MORE .LE. 0) REWIND (LINSK)
      RETURN
C
  100 CONTINUE
      IFLAG = NREC
      WRITE (LOUT, '(/A,/A,I5)')
     1   ' Error reading surface linkfile,',
     2   ' SUBROUTINE SKINIT record index #', IFLAG
      IF (MORE .LE. 0) REWIND (LINSK)
C
  300 FORMAT (10X,'ISKWRK MUST BE DIMENSIONED AT LEAST ',I5)
  350 FORMAT (10X,'RSKWRK MUST BE DIMENSIONED AT LEAST ',I5)
  375 FORMAT (10X,'CSKWRK MUST BE DIMENSIONED AT LEAST ',I5)
C
C     end of SUBROUTINE SKINIT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKINU   (IR, NDIM, ISKWRK, RSKWRK, NSPEC, KI, NU)
C
C  START PROLOGUE
C
C  SUBROUTINE SKINU   (IR, NDIM, ISKWRK, RSKWRK, NSPEC, KI, NU)
C  Returns the number of species in a surface reaction, and the
C  species indices and stoichiometric coefficients.
C
C  INPUT
C  IR        - Integer scalar, index number of a surface reaction;
C              IR must be greater than 0 and less than or equal to
C              IISUR, the total surface reaction count.
C  NDIM      - Integer scalar, dimension of the arrays KI and NU;
C              NDIM must be at least MAXSPR, the total number of
C              species allowed in a surface reaction.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  NSPEC     - Integer scalar, the number of species (reactants +
C              products) in surface reaction IR.
C  KI(*)     - Integer array, species indices for the species in surface
C              reaction IR;
C              dimension at least MAXSPR, the total number of species
C              allowed in a surface reaction.
C  NU(*)     - Integer array, stoichiometric coefficients of the
C              species in surface reaction IR;
C              dimension at least MAXSPR, the total number of species
C              allowed in a surface reaction.
C              NU is negative if the Nth species is a reactant;
C              NU is positive if the Nth species is a product.
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
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), KI(*), NU(*)
C
      NSPEC = 0
      DO 50 N = 1, NDIM
         KI(N) = 0
         NU(N) = 0
   50 CONTINUE
C
      IF (ISKWRK(IiNIIS).LE.0 .OR. IR.LE.0 .OR.
     1    IR.GT.ISKWRK(IiNIIS).OR. MAXSPR.GT.NDIM) RETURN
C
      I_NK = ISKWRK(IiNUNK) + (IR-1)*MAXSPR - 1
      I_NU = ISKWRK(IiNU)   + (IR-1)*MAXSPR - 1
      DO 200 N = 1, MAXSPR
         K = ISKWRK(I_NK + N)
         IF (K .NE. 0) THEN
            NSPEC = NSPEC + 1
            KI(NSPEC) = K
            NU(NSPEC) = ISKWRK(I_NU + N)
         ENDIF
200   CONTINUE
C
C     end of SUBROUTINE SKINU
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKIORD (IDIM, KDIM, ISKWRK, RSKWRK, NFORD, IFORD, FORD,
     1                   NRORD, IRORD, RORD)
C
C  START PROLOGUE
C
C  SUBROUTINE SKIORD (IDIM, KDIM, ISKWRK, RSKWRK, NFORD, IFORD, FORD,
C                     NRORD, IRORD, RORD)
C  Returns the number and indices of surface reactions with modified
C  species orders, and the order values for the species in the
C  surface mechanism.
C
C  INPUT
C  IDIM      - Integer scalar, dimension of arrays IFORD and IRORD;
C              IDIM must be at least NORD, the total number of
C              surface reactions with modified species orders.
C  KDIM      - Integer scalar, first dimension of the arrays FORD and
C              RORD;
C              KDIM must be at least NKK, the total species count.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  NFORD     - Integer scalar, total number of surface reactions with
C              modified forward species orders.
C  IFORD(*)  - Integer array, indices of surface reactions with modified
C              forward species orders; dimension at least NFORD.
C  FORD(*,*) - Real matrix, the modified forward species orders for the
C              NFORD surface reactions;
C              dimension at least KKTOT, the total species count, for
C              the first, and at least NFORD for the second.
C              FORD(K,N) is the forward order of species K for the Nth
C              surface change-order reaction.
C  NRORD     - Integer scalar, total number of surface reactions with
C              modified reverse species orders.
C  IRORD(*)  - Integer array, indices of surface reactions with modified
C              reverse species orders; dimension at least NRORD.
C  RORD(*,*) - Real matrix, the modified reverse species orders for the
C              NRORD surface reactions;
C              dimension at least KKTOT for the first, the total species
C              count, and at least NRORD for the second.
C              RORD(K,N) is the reverse order of species K for the Nth
C              surface change-order reaction.
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
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ISKWRK(*), RSKWRK(*), IFORD(*), FORD(KDIM,*),
     1          IRORD(*), RORD(KDIM,*)
      LOGICAL LFORD, LRORD
C
      NFORD = 0
      NRORD = 0
      DO 100 N = 1, IDIM
         IFORD(N) = 0
         IRORD(N) = 0
         DO 50 K = 1, KDIM
            FORD(K,N) = 0.0
            RORD(K,N) = 0.0
   50    CONTINUE
  100 CONTINUE
C
      IF (IDIM.LE.0 .OR. IDIM.LT.ISKWRK(IiNORD) .OR.
     1    KDIM.LT.ISKWRK(IiKTOT)) RETURN
C
C     initialize some local pointers in COMMON /SKPTR/
      CALL SKLOC (ISKWRK)
      DO 200 N = 1, ISKWRK(IiNORD)
         LFORD = .FALSE.
         LRORD = .FALSE.
         I  = ISKWRK(I_IORD + N - 1)
         DO 150 K = 1, MAXORD
            KSPEC = ISKWRK(I_KORD + N - 1)
            IF (KSPEC .EQ. 0) THEN
               I_KORD = I_KORD + MAXORD
               I_ORD  = I_ORD + MAXORD
               GO TO 200
            ELSEIF (KSPEC .LT. 0) THEN
               IF (.NOT. LFORD) THEN
                  LFORD = .TRUE.
                  NFORD = NFORD + 1
                  IFORD(N) = I
               ENDIF
               FORD(KSPEC,NFORD) = RSKWRK(I_ORD + N - 1)
            ENDIF
  150    CONTINUE
  200 CONTINUE
C
C     end of SUBROUTINE SKIORD
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKIREV (IR, ISKWRK, RSKWRK, IREV, RAR, RBR, RER)
C
C  START PROLOGUE
C
C  SUBROUTINE SKIREV (IR, ISKWRK, RSKWRK, IREV, RAR, RBR, RER)
C  Returns an integer flag to indicate whether reaction IR has an
C  explicitly assigned reverse rate constant.  It also returns the
C  reverse Arrhenius expression values for surface reaction IR,
C  if it was explicitly assigned in the Surface Chemkin interpreter.
C  If reverse Arrhenius values were not explicitly assigned,
C  RAR, RBR and RER will be zero.
C
C  INPUT
C  IR        - Integer scalar, surface reaction index.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  IREV      - Integer scalar,
C              1, reaction IR has explicit reverse rate parameters
C              0, no.
C  RAR       - Real scalar, explicit pre-exponential constants
C              for reaction IR.
C                 cgs units, mole-cm-sec-K
C  RBR       - Real scalar, explicit temperature dependence exponents
C              for reaction IR.
C  RER       - Real scalar, explicit activation energy for reaction IR.
C                 cgs units, Kelvins
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      ISIREV = CKLKUP (IR, ISKWRK(ISKWRK(IiIREV)), ISKWRK(IiNREV))
      IF (ISIREV .GT. 0) THEN
         IREV = 1
         I_PAR = ISKWRK(IrRPAR) + (ISIREV-1)*(NSPAR+1)
         RAR = RSKWRK(I_PAR)
         RBR = RSKWRK(I_PAR + 1)
         RER = RSKWRK(I_PAR + 2)
      ELSE
         RAR = 0.0
         RBR = 0.0
         RER = 0.0
         IREV = 0
      ENDIF
C
C     end of SUBROUTINE SKIREV
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKIRNU (IDIM, NDIM, ISKWRK, RSKWRK, NIIRNU, IRNU,
     1                   NSPEC, KI, RNU)
C
C  START PROLOGUE
C
C  SUBROUTINE SKIRNU (IDIM, NDIM, ISKWRK, RSKWRK, NIIRNU, IRNU, 
C                     NSPEC, KI, RNU)
C  Returns the number and indices of surface reactions with real
C  stoichiometric coefficients, number of species in the reactions,
C  and the species indices and coefficients;
C
C  INPUT
C  IDIM      - Integer scalar, dimension of the arrays IRNU and NSPEC,
C              and the second dimension of matrices KI and RNU;
C              IDIM must be at least NIIRNU, the number of surface
C              reactions with real stoichiometric coefficients.
C  NDIM      - Integer scalar, first dimension of matrices KI and RNU;
C              NDIM must be at least MAXSPR, the maximum number of
C              species allowed in a surface reaction.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  NIIRNU    - Integer scalar, total number of surface reactions with
C              real stoichiometric coefficients.
C  IRNU(*)   - Integer array, indices of surface reactions with real
C              stoichiometric coefficients; dimension at least NIIRNU.
C  NSPEC(*)  - Integer array, total number of species in a surface
C              reaction;
C              dimension at least NIIRNU.
C  KI(*,*)   - Integer matrix, species indices for species in a surface
C              reaction;
C              dimension at least MAXSPR for the first, and at least
C              NIIRNU for the second.
C              KI(M,N) is the species index of the Mth species in the
C              Nth real coefficient surface reaction.
C  RNU(*,*)  - Real matrix, stoichiometric coefficients for species
C              in the NIIRNU reactions; dimension at least MAXSPR for
C              the first, and at least NIIRNU for the second.
C              RNU(M,N) is the stoichiometric coefficient of the Mth
C              species in the Nth real coefficient surface reaction, and
C              RNU(M,*) < 0 if the Mth species is a reactant;
C              RNU(M,*) > 0 if the Mth species is a product.
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
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), IRNU(*), NSPEC(*),
     1          KI(NDIM,*), RNU(NDIM,*)
C
      DO 100 N = 1, IDIM
         NSPEC(N) = 0
         IRNU(N) = 0
         DO 50 M = 1, NDIM
            KI(M,N) = 0
            RNU(M,N) = 0.0
   50    CONTINUE
  100 CONTINUE
C
      NIIRNU = ISKWRK(IiNRNU)
      IF (NIIRNU.LE.0 .OR. IDIM.LT.ISKWRK(IiNRNU) .OR.
     1    NDIM.LT.MAXSPR) RETURN
C
      DO 200 N = 1, NIIRNU
         I = ISKWRK(ISKWRK(IiIRNU) + N - 1)
         IRNU(N) = I
         I_NK = ISKWRK(IiNUNK) + (I-1)*MAXSPR
         I_NU = ISKWRK(IrRNU) + (N-1) *MAXSPR
         DO 150 M = 0, MAXSPR - 1
            IF (ISKWRK(I_NK + M) .NE. 0) THEN
               NSPEC(N) = NSPEC(N) + 1
               KI(NSPEC(N),N) = ISKWRK(I_NK + M)
               RNU(NSPEC(N),N) = RSKWRK(I_NU + M)
            ENDIF
  150    CONTINUE
  200 CONTINUE
C
C     end of SUBROUTINE SKIRNU
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKISTK (IR, ISKWRK, ISTFL)
C
C  START PROLOGUE
C
C  SUBROUTINE SKISTK (IR, ISKWRK, ISTFL)
C  Returns an integer flag to indicate whether reaction IR uses
C  sticking coefficients.
C
C  INPUT
C  IR        - Integer scalar, index of a surface reaction.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C
C  OUTPUT
C  ISTFL     - Integer scalar,
C              0, reaction IR does not use sticking coefficients
C              1, reaction IR does use sticking coefficients
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      ISTFL = CKLKUP (IR, ISKWRK(ISKWRK(IiISTK)), ISKWRK(IiNSTK))
      ISTFL = MIN (ISTFL, 1)
C
C     end of SUBROUTINE SKISTK
      RETURN
      END
C                                                                      C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKIYLD (IR, ISKWRK, RSKWRK, IYLD, IYION, KYLD, PYLD)
C
C  START PROLOGUE
C
C  SUBROUTINE SKIYLD (IR, ISKWRK, RSKWRK, IYLD, IYION, KYLD, PYLD)
C  Returns an integer flag to indicate whether reaction IR has yield-
C  modified species, the species index of its ion, yield-modify flags
C  for its reactants and products, and parameters for the yield
C  expression.
C
C  INPUT
C  IR        - Integer scalar, surface reaction index.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  IYLD      - Integer scalar, flag for yield-modify reactions;
C              1, reaction IR uses yield-modification
C              0, no
C  IYION     - Integer scalar, species index of the ion in a yield-
C              modify reaction.
C  KYLD(*)   - Integer array, yield flags for the species in a yield-
C              modify reaction;
C              dimension at least MAXSPR, the maximum number of species
C              allowed in a surface reaction.
C              1, species is yield-modified
C              0, no
C  PYLD(*)   - Real array, parameters for the yield-expression in
C              a yield-modify reaction;
C              dimension at least NYPAR, the number of parameters
C              required.
C              If IYLD=1, and KYLD of the Nth species in the reaction
C              is 1, the stoichiometric coefficient NU of the species is
C              scaled by the results of the expression
C                 PYLD(1) * [Ei**PYLD(3) - PYLD(2)**PYLD(3)] **PYLD(4)
C              where Ei is the ion energy of species IYION.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), KYLD(*), PYLD(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      IYLD = 0
      NIIYLD = ISKWRK(IiNYLD)
      IF (NIIYLD .LE. 0) RETURN
      N = CKLKUP (IR, ISKWRK(ISKWRK(IiIYLD)), NIIYLD)
      IF (N .GT. 0) THEN
         IYLD = 1
         IYION = ISKWRK(ISKWRK(IiYION) + N - 1)
         DO 50 K = 1, MAXSPR
            KYLD(K) = ISKWRK( ISKWRK(IiKYLD) + (N-1)*MAXSPR + K - 1)
   50    CONTINUE
         DO 60 L = 1, NYPAR
            PYLD(L) = RSKWRK(ISKWRK(IrPYLD + L - 1))
   60    CONTINUE
      ENDIF
C
C     end of SUBROUTINE SKIYLD
      RETURN
      END
C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKKFRT (P, T, ISKWRK, RSKWRK, RKFT, RKRT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKKFRT (P, T, ISKWRK, RSKWRK, RKFT, RKRT)
C  Returns the temperature-dependent forward and reverse reaction
C  rate coefficients for reactions given pressure and temperature.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  RKFT(*)   - Real array, forward reaction rates for reactions;
C              dimension at least IISUR, the total reaction count.
C                 cgs units, depends on the reaction
C  RKRT(*)   - Real array, reverse reaction rates for reactions;
C              dimension at least IISUR, the total reaction count.
C                 cgs units, depends on the reaction
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ISKWRK(*), RSKWRK(*), T(*), RKFT(*), RKRT(*)
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
C     initialize some local pointers in COMMON /SKPTR/
      CALL SKLOC (ISKWRK)
      CALL SKRRPT (ISKWRK, RSKWRK, RSKWRK(I_SMH), T, RSKWRK(I_KWT),
     1             ISKWRK(I_ NREA), ISKWRK(I_NRPP), ISKWRK(I_NU), 
     2             ISKWRK(I_NK), ISKWRK(I_IRNU), RSKWRK(I_RNU), 
     3             ISKWRK(I_NUSUM), RSKWRK(I_PAR), RSKWRK(I_RPAR),
     4             ISKWRK(I_IREV), ISKWRK(I_ISTK), ISKWRK(I_MSTK), 
     5             ISKWRK(I_IBHM), ISKWRK(I_KBHM), ISKWRK(I_IEDP),
     6             ISKWRK(I_KEDP), RSKWRK(I_PEDP), ISKWRK(I_KTFL), 
     7             ISKWRK(I_IORD), ISKWRK(I_KORD), RSKWRK(I_ORD),
     8             ISKWRK(I_IYLD), ISKWRK(I_YION), ISKWRK(I_KYLD), 
     9             RSKWRK(I_PYLD), RKFT, RKRT, RSKWRK(I_EQKC),
     *             RSKWRK(I_EQFAC), RSKWRK(I_ENRG))
C
C     end of SUBROUTINE SKKFRT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKKION (ISKWRK, KEL, KKION, KION)
C
C  START PROLOGUE
C
C  SUBROUTINE SKKION (ISKWRK, KEL, KKION, KION)
C  Returns the species number of the electron, the number of positive
C  ions in the gas phase, and an array of species number for each
C  positive ion
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C
C  OUTPUT
C  KEL       - Integer scalar, species index of the electron species.
C  KKION     - Integer scalar, total gas-phase positive ion count.
C  KION      - Integer array, species indices for the gas-phase positive
C              ions;
C              dimension at least NKKGAS, the gas-phase species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), KION(*)
C
      KEL = ISKWRK(IiELEC)
      KKION = ISKWRK(IiKION)
      IF (KKION.EQ.0) RETURN
C
      DO 100 K = 1, KKION
         KION(K) = ISKWRK(ISKWRK(IiIONS) + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE SKKION
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKKTFL (ISKWRK, KTFL)
C
C  START PROLOGUE
C
C  SUBROUTINE SKKTFL (ISKWRK, KTFL)
C
C  Allows the user to assign a location in the temperature array
C  to use for the gas-phase species.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  KTFL(*)   - Integer array, ndices into the temperature(s) for
C              gas-phase species;
C              dimension at least KKGAS, the total gas-phase species
C              count.
C              Default value stored in ISKWRK is set to 1 in SKINIT.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), KTFL(*)
C
      I_KTFL = ISKWRK(IiKTFL) - 1
      DO 100 K = 1, NKKGAS
         ISKWRK(I_KTFL+K) = KTFL(K)
  100 CONTINUE
C
C     end of SUBROUTINE SKKTFL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKLEN  (LINSK, LOUT, LENI, LENR, LENC, IFLAG)
C
C  START PROLOGUE
C
C  SUBROUTINE SKLEN  (LINSK, LOUT, LENI, LENR, LENC, IFLAG)
C  Reads the first record of the linkfile to return the lengths
C  required for the integer, real, and character work arrays.
C
C  INPUT
C  LINSK     - Integer scalar, input unit number assigned to linkfile.
C  LOUT      - Integer scalar, formatted output unit file number.
C
C  OUTPUT
C  LENI      - Integer scalar, dimension required for integer work
C              array, ISKWRK.
C  LENR      - Integer scalar, dimension required for real work
C              array, RSKWRK.
C  LENC      - Integer scalar, dimension required for character work
C              array, CSKWRK.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (NLIST=1)
      LOGICAL KERR, VOK, POK, LBIN
      CHARACTER*16 LIST(NLIST), FILVER, PREC, PRVERS, IFMT, CFMT,
     1             RFMT, LFMT
      PARAMETER (IFMT='(10I12)', CFMT='(8A16)', RFMT='(1P,5E24.16)',
     1           LFMT='(L8)')
C
      COMMON /SKCONS/ PREC, FILVER, KERR
      DATA LIST(1) /'1.0'/
C
      FILVER = ' '
      PRVERS  = ' '
      PREC = ' '
      LENI = 0
      LENR = 0
      LENC = 0
C
      IFLAG = 0
      KERR = .FALSE.
      MORE = 0
C
C*****linkfile (surface) > binary
C      LBIN = .TRUE.
C*****END linkfile (surface) > binary
C*****linkfile (surface) > ascii
      LBIN = .FALSE.
C*****END linkfile (surface) > ascii
C
      NREC = 1
      IF (LBIN) THEN
         READ (LINSK, ERR=100) FILVER
      ELSE
         READ (LINSK, CFMT, ERR=100) FILVER
      ENDIF
      VOK = .FALSE.
      DO 5 N = 1, NLIST
         IF (FILVER .EQ. LIST(N)) VOK = .TRUE.
    5 CONTINUE
C
      IF (LBIN) THEN
         NREC = 2
         READ (LINSK, ERR=100) PRVERS
         NREC = 3
         READ (LINSK, ERR=100) PREC
         NREC = 4
         READ (LINSK, ERR=100) KERR
      ELSE
         NREC = 2
         READ (LINSK, CFMT, ERR=100) PRVERS
         NREC = 3
         READ (LINSK, CFMT, ERR=100) PREC
         NREC = 4
         READ (LINSK, LFMT, ERR=100) KERR
      ENDIF
C
      POK = .FALSE.
C*****precision > double
      IF (INDEX(PREC, 'DOUB') .GT. 0) POK = .TRUE.
C*****END precision > double
C
C*****precision > single
C      IF (INDEX(PREC, 'SING') .GT. 0) POK = .TRUE.
C*****END precision > single
C
      IF (KERR .OR. (.NOT.POK) .OR. (.NOT.VOK)) THEN
         IF (KERR) THEN
            WRITE (LOUT,'(/A,/A)')
     1      ' There is an error in the surface linkfile...',
     2      ' Check SURFACE INTERPRETER output for error conditions.'
         ENDIF
C
         IF (.NOT. VOK) WRITE (LOUT, '(/A)')
     1   ' Surface linkfile is incompatible with Surface Library.'
C
         IF (.NOT. POK) THEN
            WRITE (LOUT,'(/A,A)')
     1      ' Precision of surface linkfile does not agree with',
     2      ' precision of surface library'
         ENDIF
         IFLAG = 20
         DO 15 N = 1, NREC
            BACKSPACE (LINSK)
   15    CONTINUE
         RETURN
      ENDIF
C
      NREC = 5
      IF (LBIN) THEN
         READ (LINSK, ERR=100) LENI, LENR, LENC
      ELSE
         READ (LINSK, IFMT, ERR=100) LENI, LENR, LENC
      ENDIF
C
      DO 10 N = 1, NREC
         BACKSPACE (LINSK)
   10 CONTINUE
      RETURN
C
  100 CONTINUE
      IFLAG = NREC
      WRITE (LOUT,'(/A,/A,I5)')
     1   ' Error reading surface linkfile,',
     2   ' SUBROUTINE SKLEN record index #',IFLAG
      DO 20 N = 1, NREC
         BACKSPACE (LINSK)
   20 CONTINUE
C
C     end of SUBROUTINE SKLEN
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKLOC (ISKWRK)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C  START PROLOGUE
C
C  SUBROUTINE SKLOC
C  Initialize some local pointers into the arrays and store them
C  in COMMON /SKPTR/
C
C  END PROLOGUE
C
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ISKWRK(*)
C
C     convert some ISKWRK pointers to local indices:
      I_CZ  = ISKWRK(IrKT1)
      I_SMH = ISKWRK(IrKT2)
      I_KWT = ISKWRK(IrKWT)
      I_NREA = ISKWRK(IiNREA)
      I_NRPP = ISKWRK(IiNRPP)
      I_NU   = ISKWRK(IiNU)
      I_NK   = ISKWRK(IiNUNK)
      I_IRNU = ISKWRK(IiIRNU)
      I_RNU  = ISKWRK(IrRNU)
      I_NUSUM = ISKWRK(IiNSUM)
      I_PAR  = ISKWRK(IrPAR)
      I_RPAR = ISKWRK(IrRPAR)
      I_IREV = ISKWRK(IiIREV)
      I_ICOV = ISKWRK(IiICOV)
      I_KCOV = ISKWRK(IiKCOV)
      I_CPAR = ISKWRK(IrKCOV)
      I_ISTK = ISKWRK(IiISTK)
      I_MSTK = ISKWRK(IiMSTK)
      I_IBHM = ISKWRK(IiIBHM)
      I_KBHM = ISKWRK(IiKBHM)
      I_IEDP = ISKWRK(IiIEDP)
      I_KEDP = ISKWRK(IiKEDP)
      I_PEDP = ISKWRK(IrPEDP)
      I_KTFL = ISKWRK(IiKTFL)
      I_IORD = ISKWRK(IiIORD)
      I_KORD = ISKWRK(IiKORD)
      I_ORD = ISKWRK(IrKORD)
      I_IYLD = ISKWRK(IiIYLD)
      I_YION = ISKWRK(IiYION)
      I_KYLD = ISKWRK(IiKYLD)
      I_PYLD = ISKWRK(IrPYLD)
      I_RKFT = ISKWRK(IrKFT)
      I_RKRT = ISKWRK(IrKRT)
      I_RKF  = ISKWRK(IrIT1)
      I_RKR  = ISKWRK(IrIT2)
      I_EQKC = ISKWRK(IrIT3)
      I_EQFAC= ISKWRK(IrEQ)
      I_ENRG = ISKWRK(IrENGI)
C
C     end of SUBROUTINE SKLOC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKMXTP (ISKWRK, MXTP)
C
C  START PROLOGUE
C
C  SUBROUTINE SKMXTP (ISKWRK, MXTP)
C  Returns the maximum number of temperatures used in
C  fitting the thermodynamic properties of the species.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C
C  OUTPUT
C  MXTP      - Integer scalar, maximum number of temperatures used in
C              fitting the thermodynamic properties of the species.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*)
C
      MXTP = MAXTP
C
C     end of SUBROUTINE SKMXTP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNCF  (NELDIM, ISKWRK, NEL)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNCF (NELDIM, ISKWRK, NEL)
C  Returns the elemental composition of the species.
C
C  INPUT
C  NELDIM    - Integer scalar, first dimension of the matrix NEL;
C              must be at least NELEM, the total element count.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C
C  OUTPUT
C  NEL(*,*)  - Integer matrix, elemental compositions of the species;
C              dimension at least NELEM for the first, the total
C              element count, and at least KKTOT for the second, the
C              total species count.
C              NEL(M,K) is the quantity of element M in species K.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), NEL(NELDIM,*)
C
      IF (NELDIM .LT. NELEM) RETURN
      DO 100 M = 1, NELEM
         DO 50 N = 1, ISKWRK(IiKTOT)
            NEL(M,N) = ISKWRK(ISKWRK(IiKCMP) + (N-1)*NELEM + M - 1)
   50    CONTINUE
  100 CONTINUE
C
C     end of SUBROUTINE SKNCF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNCON (ISKWRK, RSKWRK, NCON)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNCON (ISKWRK, RSKWRK, NCON)
C  Returns the total number of surface reactions which do not conserve
C  sites of the phases.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C
C  OUTPUT
C  NCON(*)   - Integer array, count of surface reactions which do not
C              conserve sites in the phases;
C              dimension at least NPHASE, the total phase count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), NCON(*)
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      NIISUR = ISKWRK(IiNIIS)
      NPHASE = ISKWRK(IiNPHA)
      DO 50 N = 1, NPHASE
         NCON(N) = 0
         DO 25 I = 1, NIISUR
            RNCF = RSKWRK (ISKWRK(IrNCF) + (I-1)*NPHASE + N - 1)
            IF (ABS(RNCF) .GT. RSKWRK(ISKWRK(IrSKMN)))
     1         NCON(N) = NCON(N) + 1
   25    CONTINUE
   50 CONTINUE
C
C     end of SUBROUTINE SKNCON
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNU   (IDIM, ISKWRK, RSKWRK, KSTOIC, NSTOIC)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNU   (IDIM, ISKWRK, RSKWRK, KSTOIC, NSTOIC)
C  Returns the stoichiometric coefficients of the species and the net
C  change in phases for all of the surface reactions in a mechanism.
C
C  INPUT
C  IDIM      - Integer scalar, first dimension of the array NSTOIC;
C              must be at least IISUR, the total surface reaction count.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  KSTOIC(*,*)-Integer matrix, stoichiometric coefficients for the
C              species in the surface reactions;
C              the first dimension must be at least IISUR, the total
C              surface reaction count, and at least KKTOT for the
C              second, the total species count.
C  NSTOIC(*,*)-Integer matrix, net change of the phases for the surface
C              reactions;
C              the first dimension must be at least IISUR, the total
C              surface reaction count, and at least NPHASE for the
C              second the total phase count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), KSTOIC(IDIM,*), NSTOIC(IDIM,*)
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
C             PROCESS EACH REACTION
C
      NIISUR = ISKWRK(IiNIIS)
      DO 10 I = 1, NIISUR
         DO 05 K = 1, ISKWRK(IiKTOT)
            KSTOIC(I,K) = 0
   05    CONTINUE
   10 CONTINUE
C
      DO 100 N = 1, MAXSPR
         DO 50 I = 1, NIISUR
            K = ISKWRK(ISKWRK(IiNUNK) + (I-1)*MAXSPR + N - 1)
            NU = ISKWRK(ISKWRK(IiNU)  + (I-1)*MAXSPR + N - 1)
            IF (K .GT. 0) KSTOIC(I,K) = KSTOIC(I,K) + NU
   50    CONTINUE
  100 CONTINUE
C
      NPHASE = ISKWRK(IiNPHA)
      IF (NPHASE .LE. 0) RETURN
      DO 200 N = 1, NPHASE
         DO 150 I = 1, NIISUR
            NSTOIC(I,N) = RSKWRK(ISKWRK(IrNCF) + (N-1)*NIISUR + I - 1)
     1                  + RSKWRK(ISKWRK(IrSKMN))
  150    CONTINUE
200   CONTINUE
C
C     end of SUBROUTINE SKNU
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNUF   (IDIM, ISKWRK, KSTOIF)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNUF   (IDIM, ISKWRK, KSTOIF)
C  Returns the stoichiometric coefficients of the species
C  for all reactants in all surface reactions in a mechanism.
C  (note - reactants only! - they will all be negative)
C
C  INPUT
C  IDIM      - Integer scalar, first dimension of the array NSTOIC;
C              must be at least ISUR, the total surface reaction count.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C
C  OUTPUT
C  KSTOIF(*,*)-Integer matrix,  stoichiometric coefficients for the
C              reactants in the surface reactions;
C              dimension at least IISUR for the first, the total surface
C              reaction  ount, and at least KKTOT for the second, the
C              total species count.
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
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), KSTOIF(IDIM,*)
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
C             PROCESS EACH REACTION
C
      DO 10 I = 1, ISKWRK(IiNIIS)
         DO 05 K = 1, ISKWRK(IiKTOT)
            KSTOIF(I,K) = 0
   05    CONTINUE
   10 CONTINUE
C
      I_NK = ISKWRK(IiNUNK)
      I_NU = ISKWRK(IiNU)
      DO 20 N = 0, MAXSPR/2 - 1
         DO 15 I = 0, ISKWRK(IiNIIS) - 1
            K  = ISKWRK(I_NK + I*MAXSPR + N)
            IF (K .GT. 0) THEN
               NU = ISKWRK(I_NU  + I*MAXSPR + N)
               KSTOIF(I+1,K) = KSTOIF(I+1,K) + NU
            ENDIF
   15    CONTINUE
   20 CONTINUE
C
C     end of SUBROUTINE SKNUF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKPCMP (ISTR, IRAY, NN, SETS, NSETS, ISET, IND, NT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKPCMP (ISTR, IRAY, NN, SETS, NSETS, ISET, IND, NT)
C  This subroutine can do everything that the subroutine SKCOMP can do,
C  and additionally, has the capabilities of separating the elements of
C  IRAY into categories and then search IRAY by element and category.
C  The categories that each element of IRAY will be assigned to are
C  specified by the input character string vector SETS of vector length
C  NSETS.  Elements of each category in IRAY must be grouped congrously.
C  The number of elements in each category within IRAY is specified by
C  the input integer vector ISET.  To search for the existence of an
C  element within acategory ISTR may additionally be composed of two
C  substrings, ISTR="ELEMENT_NAME/CATEGORY_NAME/", where CATEGORY_NAME
C  is one of the categories specified in SETS.  In this case, IND will
C  return the first position in IRAY where ELEMENT_NAME occurred within
C  the category CATEGORY_NAME.  NT will return the total number of
C  times ELEMENT_NAME occurred within the category CATEGORY_NAME.
C  If ELEMENT_NAME is not found within the specified category, IND and
C  NT are returned with a value of zero.  If no category is specified
C  within ISTR, IND and NT return with the same values as they would
C  from subroutine SKCOMP.
C
C  Consider the following example,
C        IRAY = {"RED", "BLUE", "JADE", "RUBY", "TOPAZ", "JADE"}
C        NN = 6
C        SETS = {"COLORS", "STONES"},
C        NSETS = 2
C        ISET = {4, 2}.
C  This assumes that the elements of IRAY were grouped into two
C  sets, consisting of 4 and 2 elements, respectively, and the
C  following names
C        "COLORS"  = {"RED", "BLUE", "JADE", "RUBY"}, and
C        "STONES"  = {"TOPAZ", "JADE"}.
C
C  If ISTR="BLUE" then IND=2 and NT=1;
C  if ISTR="PINK" then IND=0 and NT=0; and
C  if ISTR="JADE",then IND=3 and NT=2.
C
C  If ISTR="BLUE/COLORS/" then IND=2 and NT=1;
C  if ISTR="BLUE/STONES/" then IND=0 and NT=0;
C  if ISTR="JADE/GEMS/"   then IND=0 and NT=0; and
C  if ISTR="JADE/STONES/",then IND=6 and NT=1.
C
C  INPUT
C  ISTR        - Character string, which may or may not end with a
C                slash-delimited substring.
C  IRAY(*)     - Character string array;
C                dimension at least NN.
C  NN          - Integer scalar, number of entries in IRAY(*).
C  SETS(*)     - Character string array, cross-reference set to relate
C                with elements of IRAY;
C                dimension at least NSETS.
C  NSETS       - Integer scalar, number of entries in SETS(*)
C  ISET(*)     - Integer array, total number of entries in a subset of
C                IRAY;  dimension at least NSETS.
C
C  OUTPUT
C  IND         - Integer scalar, index of ISTR in IRAY(*).
C                If ISTR is not in IRAY(*), IND = 0.
C                If the slash-delimited substring of ISTR is not
C                in SETS(*), IND = 0.
C                If the slash-delimited substring of ISTR is in
C                SETS(N), but the substring before the slash is
C                not a member of the subset associated with SETS(N),
C                IND = 0, whether or not the substring is in IRAY(*).
C  NT          - Integer scalar, total occurrence of ISTR in IRAY(*),
C                or total number of times ISTR occurs in a subset
C                of IRAY(*).
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISET(*)
      CHARACTER*(*) ISTR, IRAY(*), SETS(*)
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IND = 0
      NT  = 0
C
      IF (ISTR .EQ. ' ') RETURN
C
      I2 = CKLSCH(ISTR)
      IF (ISTR(I2:I2) .NE. '/') THEN
         CALL SKCOMP (ISTR, IRAY, NN, IND, NT)
         RETURN
      ENDIF
C
      I1 = 0
      DO 10 L = I2-1, 1, -1
         IF (ISTR(L:L).EQ.'/' .AND. I1.EQ.0) I1 = L
   10 CONTINUE
      IF (I1 .LE. 0) RETURN
C
      K = 0
      DO 50 N = 1, NSETS
         DO 25 J = 1, ISET(N)
            K = K + 1
            IF (ISTR(I1+1:I2-1) .EQ. SETS(N) .AND.
     1          ISTR(1:I1-1) .EQ. IRAY(K) .AND. IND.EQ.0) THEN
                IND = K
                NT = NT + 1
            ENDIF
   25    CONTINUE
   50 CONTINUE
C
C     end of SUBROUTINE SKPCMP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKPKK  (ISKWRK, KKPHAS, KFIRST, KLAST)
C
C  START PROLOGUE
C
C  SUBROUTINE SKPKK  (ISKWRK, KKPHAS, KFIRST, KLAST)
C  Returns arrays of species pointers for the phases.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C
C  OUTPUT
C  KKPHAS(*) - Integer array, the total species counts for phases;
C              dimension at least NPHASE, the total phase count.
C  KFIRST(*) - Integer array, species indices for the first species of
C              the phases;
C              dimension at least NPHASE, the total phase count.
C  KLAST(*)  - Integer array, species indices for the last species of
C              the phases;
C              dimension at least NPHASE, the total phase count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), KFIRST(*), KLAST(*), KKPHAS(*)
C
      IF (ISKWRK(IiNPHA) .LE. 0) RETURN
C
      DO 50 N = 1, ISKWRK(IiNPHA)
         KKPHAS(N) = ISKWRK(ISKWRK(IiPTOT) + N - 1)
         KFIRST(N) = ISKWRK(ISKWRK(IiPKST) + N - 1)
         KLAST(N)  = ISKWRK(ISKWRK(IiPKND) + N - 1)
   50 CONTINUE
C
C     end of SUBROUTINE SKPKK
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKPNT  (LSAVE, LOUT, V, P, LENI, LENR, LENC, IERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKPNT  (LSAVE, LOUT, VERS, PREC, LENI, LENR, LENC, IERR)
C  Reads from a file information about a Surface Chemkin linkfile,
C  pointers for the Surface Chemkin Library, and returns lengths
C  of work arrays.
C
C  INPUT
C  LSAVE     - Integer scalar, input unit for binary data file.
C  LOUT      - Integer scalar, formatted output file unit number.
C
C  OUTPUT
C  VERS      - Real scalar, version number of the Surface Chemkin
C              linkfile.
C  PREC      - Character string, machine precision of the linkfile.
C  LENI      - Integer scalar, length required for integer work array.
C  LENR      - Integer scalar, length required for real work array.
C  LENC      - Integer scalar, length required for character work array.
C  KERR      - Logical, error flag.
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
      INCLUDE 'skstrt.h'
      COMMON /MACHN/ SMALL, BIG, EXPARG
      COMMON /SKCONS/ PREC, FILVER, KERR
C
      LOGICAL KERR, IERR
      CHARACTER*16 PREC, FILVER, V, P
C----------------------------------------------------------------------C
C     Data about the machine dependent constants is carried in         C
C     COMMON/MACHN/SMALL,BIG,EXPARG
C
C      THIS STATEMENT WILL NOT COMPILE, MACHINE-DEPENDENT CONSTANTS
C
C*****exponent range > +/-30
C      SMALL = 1.0E-30
C      BIG   = 1.0E+30
C*****END exponent range > +/-30
C*****exponent range > +/-300
      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
C*****END exponent range > +/-300
C
      EXPARG = LOG(BIG)
C
C----------------------------------------------------------------------C
C
      KERR = .FALSE.
      IERR = KERR
C
      READ (LSAVE, ERR=100) FILVER, PREC, LENI, LENR, LENC,
C
C     Include file for CHEMKIN-III sklib.f, dated:  March 1, 1966
C
C     Integer constants
C
     1   MAXSPR, NELEM, NKKGAS, NSPAR, NSCOV, NEDPAR, NYPAR, MAXORD,
     2   MAXTP, NCP,    NCP1,  NCP2,   NCP2T,
C
C     ISKWRK pointers to integer variables
C
     3   IiLENI, IiLENR, IiLENC, IiKSUR, IiKBLK, IiKTOT, IiNPHA, IiFSUR,
     4   IiLSUR, IiNSUR, IiFBLK, IiLBLK, IiNBLK, IiNIIS, IiNCOV, IiNREV,
     5   IiNSTK, IiNCON, IiNBHM, IiNRNU, IiNORD, IiMOTZ, IiELEC, IiNYLD,
C
C     ISKWRK pointers to integer arrays
C
     6   IiPKST, IiPKND, IiPTOT, IiKPHS, IiKCHG, IiKCMP, IiNSCV, IiKNT,
     7   IiNRPP, IiNREA, IiNUNK, IiNU,   IiNSUM, IiICOV, IiKCOV, IiIREV,
     8   IiISTK, IiMSTK, IiIBHM, IiKBHM, IiIRNU, IiIORD, IiKORD, IiIONS,
     9   IiKION, IiKTFL, IiNEDP, IiIEDP, IiKEDP, IiIYLD, IiYION, IiKYLD,
C
C     ISKWRK pointers to real variables
C
     *   IrSKMN, IrPATM, IrRU,   IrRUC,
C
C     ISKWRK pointers to real arrays
C
     1   IrSDEN, IrKTMP, IrKTHM, IrKDEN, IrAWT,  IrKWT,  IrPAR, IrKCOV,
     2   IrRPAR, IrEQ,   IrRNU,  IrNCF,  IrKORD, IrKFT,  IrKRT, IrKT1,
     3   IrKT2,  IrPT1,  IrIT1,  IrIT2,  IrIT3,  IrPEDP, IrENGI,
     4   IrPYLD, IrYNCF,
C
C     ISKWRK pointers to character arrays
C
     4   IcENAM, IcKNAM, IcMNAM, IcPNAM
C
C     END include file for sklib.f
C
      V  = FILVER
      P  = PREC
      IERR = KERR
      RETURN
C
  100 CONTINUE
      WRITE (LOUT, *) ' Error reading Surface binary file data...'
      KERR   = .TRUE.
      IERR   = KERR
      NPOINT = 0
      FILVER   = ' '
      V      = FILVER
      PREC   = ' '
      P      = PREC
C
C     end of SUBROUTINE SKPNT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRAEX (IR, ISKWRK, RSKWRK, RA)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRAEX (IR, ISKWRK, RSKWRK, RA)
C
C  Returns the Pre-exponential rate constant
C  (or sticking coefficient) of the IRth reaction, or changes its
C  value, depending on the sign of IR.
C
C  INPUT
C  IR        - Integer scalar, reaction index;
C              IR> 0 gets RA(I) from RSKWRK
C              IR< 0 puts RA(I) into RSKWRK
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  If IR< 0:
C  RA        - Real scalar, pre-exponential or sticking coefficient for
C              reaction IR.
C                 cgs units, mole-cm-sec-K for pre-exponential,
C                            none for sticking coefficients
C
C  OUTPUT
C  If IR> 0:
C  RA        - Real scalar, pre-exponential or sticking coefficient for
C              reaction IR.
C                 cgs units, mole-cm-sec-K for pre-exponential,
C                            none for sticking coefficients
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*)
C
      NI = ISKWRK(IrPAR) + (ABS(IR)-1)*(NSPAR+1)
      IF (IR .GT. 0) THEN
         RA = RSKWRK(NI)
      ELSE
         RSKWRK(NI) = RA
      ENDIF
C
C     end of SUBROUTINE SKRAEX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRAT  (P, T, ACT, SDEN, ISKWRK, RSKWRK, SDOT, SITDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRAT  (P, T, ACT, SDEN, ISKWRK, RSKWRK, SDOT, SITDOT)
C  Returns production rates for the species and sites.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  SDOT(*)   - Real array, production rates of the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, moles/(cm**2*sec)
C              for 1,KKGAS, the production rate of gas-phase species,
C              for KKGAS+1,KKGAS+KKSUR, the production rate of surface
C              species,
C              for KKGAS+KKSUR+1,KKTOT, the production rate of bulk
C              species.
C  SITDOT(*) - Real array, production rates of the surface phases;
C              dimension at least NPHASE, the total phase count, but
C              subroutine only calculates entries for site phases.
C                 cgs units, moles/(cm**2*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ISKWRK(*), RSKWRK(*), ACT(*), T(*), SDOT(*), SITDOT(*),
     1          SDEN(*)
C
      DO 10 K = 1, ISKWRK(IiKTOT)
         SDOT(K) = ZERO
   10 CONTINUE
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 20 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SITDOT(N) = ZERO
   20    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      SDTOT  = ZERO
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
         IF (ISKWRK(IiNYLD) .GT. 0) CALL SKEQY (ISKWRK, RSKWRK)
      ENDIF
C
C     initialize some local pointers in COMMON /SKPTR/
      CALL SKLOC (ISKWRK)
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, RSKWRK(I_CZ))
      CALL SKRROP (ISKWRK, RSKWRK, RSKWRK(I_SMH), T, RSKWRK(I_CZ),
     1             ACT, RSKWRK(I_KWT), ISKWRK(I_NREA), ISKWRK(I_NRPP),
     2             ISKWRK(I_NU), ISKWRK(I_NK), ISKWRK(I_IRNU),
     3             RSKWRK(I_RNU), ISKWRK(I_NUSUM), RSKWRK(I_PAR),
     4             RSKWRK(I_RPAR), ISKWRK(I_IREV), ISKWRK(I_ICOV),
     5             ISKWRK(I_KCOV), RSKWRK(I_CPAR), ISKWRK(I_ISTK),
     6             ISKWRK(I_MSTK), ISKWRK(I_IBHM), ISKWRK(I_KBHM),
     7             ISKWRK(I_IEDP), ISKWRK(I_KEDP), RSKWRK(I_PEDP),
     8             ISKWRK(I_KTFL), SDTOT, ISKWRK(I_IORD),
     9             ISKWRK(I_KORD), RSKWRK(I_ORD), ISKWRK(I_IYLD),
     *             ISKWRK(I_YION), ISKWRK(I_KYLD), RSKWRK(I_PYLD),
     1             RSKWRK(I_RKFT), RSKWRK(I_RKRT), RSKWRK(I_RKF),
     2             RSKWRK(I_RKR), RSKWRK(I_EQKC), RSKWRK(I_EQFAC),
     3             RSKWRK(I_ENRG))
C
      I_NCF  = ISKWRK(IrNCF)
      I_YNCF = ISKWRK(IrYNCF)
      CALL SKSDOT(ISKWRK(IiNIIS), RSKWRK(I_RKF), RSKWRK(I_RKR), ISKWRK,
     1            ISKWRK(I_NU), ISKWRK(I_NK), ISKWRK(I_IRNU),
     2            RSKWRK(I_RNU), RSKWRK(I_NCF), ISKWRK(I_IYLD),
     3            ISKWRK(I_KYLD), RSKWRK(I_PYLD), ISKWRK(I_YION),
     4            RSKWRK(I_ENRG), ISKWRK(IiNPHA), RSKWRK(I_YNCF), SDOT,
     5            SITDOT)
C
C     end of SUBROUTINE SKRAT
      RETURN
      END
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRATI (IR, ROP, ISKWRK, RSKWRK, SDOTI, SITDTI)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRATI (IR, ROP, ISKWRK, RSKWRK, SDOTI, SITDTI)
C  Returns rates of production of the species by surface reaction IR.
C
C  INPUT
C  IR        - Integer scalar, reaction index;
C  ROP(*)    - Real array, rates of progress for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec).
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  SDOTI(*)  - Real array, production rates of the species by reaction
C              IR;
C              dimension at least KKTOT, the total species count.
C                 cgs units, moles/(cm**2*sec)
C              for 1,KKGAS, the production rate of gas-phase species,
C              for KKGAS+1,KKGAS+KKSUR, the production rate of surface
C              species,
C              for KKGAS+KKSUR+1,KKTOT, the production rate of bulk
C              species.
C SITDTI(*)  - Real array, production rates of the surface phases due to
C              reaction IR;
C              dimension at least NPHASE, the total phase count, but
C              subroutine calculates entries only for site phases.
C                 cgs units, moles/(cm**2*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ROP(*), ISKWRK(*), RSKWRK(*), SDOTI(*), SITDTI(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      DO 50 K = 1, ISKWRK(IiKTOT)
         SDOTI(K) = 0.0
   50 CONTINUE
C
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SITDTI(N) = 0.0
   60    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
C     initialize some local pointers in COMMON /SKPTR/
      CALL SKLOC (ISKWRK)
C
C     PROCESS REACTION IR
C
      ISRNU = CKLKUP (IR, ISKWRK(I_IRNU), ISKWRK(IiNRNU))
      ISYLD = CKLKUP (IR, ISKWRK(I_IYLD), ISKWRK(IiNYLD))
C
      I_NK = I_NK + (IR-1)*MAXSPR - 1
      IF (ISRNU .GT. 0) THEN
         I_NU   = I_RNU   + (ISRNU-1)*MAXSPR - 1
      ELSE
         I_NU   = I_NU  + (IR-1)*MAXSPR - 1
      ENDIF
C
      DO 100 N = 1, MAXSPR
        K = ISKWRK(I_NK + N)
C
C       is this an active species
        IF (K .LE. 0) GO TO 100
C
C       does the species have integer, or real coefficient
        IF (ISRNU .NE. 0) THEN
           RNU = RSKWRK(I_NU + N)
        ELSE
           RNU = ISKWRK(I_NU + N)
        ENDIF
C
        IF (ISYLD .GT. 0) THEN
C          this is a yield-modify reaction
C
           KYLD = ISKWRK(I_KYLD + (ISYLD-1)*MAXSPR + N - 1)
C
C          is this a yield-modify species
           IF (KYLD .GT. 0) THEN
C
C             ion species index
              KION = ISKWRK(I_YION + ISYLD - 1)
C
C             locator for yield parameters
              IPYLD = I_PYLD + NYPAR*(ISYLD-1)
              EI    = RSKWRK(I_ENRG + KION - 1)
              ETH   = RSKWRK(IPYLD + 1)
              IF (EI .GE. ETH) THEN
                 ASCAL = RSKWRK( IPYLD)
                 A     = RSKWRK( IPYLD + 2)
                 B     = RSKWRK( IPYLD + 3)
                 SCALE = ASCAL * (EI**A - ETH**A) **B
              ELSE
                 SCALE = 0.0
              ENDIF
              RNU = RNU * SCALE
           ENDIF
        ENDIF
        SDOTI(K) = SDOTI(K) + ROP(IR)*RNU
  100 CONTINUE
C
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
C
         NPHASE = ISKWRK(IiNPHA)
         I_NCF = ISKWRK(IrNCF) + (IR-1)*NPHASE - 1
         DO 125 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SITDTI(N) = SITDTI(N) + ROP(IR)*RSKWRK(I_NCF + N)
  125    CONTINUE
C
         IF (ISYLD .GT. 0) THEN
C
            I_NCF = ISKWRK(IrYNCF) + (ISYLD-1)*ISKWRK(IiNPHA) - 1
            DO 150 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
               SITDTI(N) = SITDTI(N) + ROP(IR)*RSKWRK(I_NCF + N)
  150       CONTINUE
         ENDIF
      ENDIF
C
C     end of SUBROUTINE SKRATI
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRATK  (P, T, ACT, SDEN, RKFT, RKRT, ISKWRK, RSKWRK,
     1                    SDOT, SITDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRATK  (P, T, ACT, SDEN, RKFT, RKRT, ISKWRK, RSKWRK, 
C                      SDOT, SITDOT)
C  Returns production rates for the species and sites, using the
C  pre-calculated temperature-dependent rate coefficients RKFT, RKRT.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  RKFT(*)   - Real array, forward reaction rates for reactions;
C              dimension at least IISUR, the total reaction count.
C                 cgs units, depends on the reaction
C  RKRT(*)   - Real array, reverse reaction rates for reactions;
C              dimension at least IISUR, the total reaction count.
C                 cgs units, depends on the reaction
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  SDOT(*)   - Real array, production rates of the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, moles/(cm**2*sec)
C              for 1,KKGAS, the production rate of gas-phase species,
C              for KKGAS+1,KKGAS+KKSUR, the production rate of surface
C              species,
C              for KKGAS+KKSUR+1,KKTOT, the production rate of bulk
C              species.
C  SITDOT(*) - Real array, production rates of the surface phases;
C              dimension at least NPHASE, the total phase count, but
C              subroutine only calculates entries for site phases.
C                 cgs units, moles/(cm**2*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ISKWRK(*), RSKWRK(*), ACT(*), T(*), SDOT(*), SITDOT(*),
     1          SDEN(*)
C
      DO 10 K = 1, ISKWRK(IiKTOT)
         SDOT(K) = ZERO
   10 CONTINUE
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 20 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SITDOT(N) = ZERO
   20    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      SDTOT  = ZERO
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
         IF (ISKWRK(IiNYLD) .GT. 0) CALL SKEQY (ISKWRK, RSKWRK)
      ENDIF
C
C     initialize some local pointers in COMMON /SKPTR/
      CALL SKLOC (ISKWRK)
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, RSKWRK(I_CZ))
      CALL SKRRPX (ISKWRK, RSKWRK, T, RSKWRK(I_CZ), ACT, ISKWRK(I_NREA),
     1             ISKWRK(I_NRPP), ISKWRK(I_NU), ISKWRK(I_NK), 
     2             ISKWRK(I_IRNU), RSKWRK(I_RNU), RSKWRK(I_PAR), 
     2             ISKWRK(I_IREV), ISKWRK(I_ICOV), ISKWRK(I_KCOV),
     3             RSKWRK(I_CPAR), ISKWRK(I_ISTK), ISKWRK(I_IBHM), 
     4             ISKWRK(I_KBHM), ISKWRK(I_KTFL), SDTOT, 
     5             ISKWRK(I_IORD), ISKWRK(I_KORD), RSKWRK(I_ORD), 
     6             RKFT, RKRT, RSKWRK(I_RKF), RSKWRK(I_RKR))
C
      I_NCF  = ISKWRK(IrNCF)
      I_YNCF = ISKWRK(IrYNCF)
      CALL SKSDOT(ISKWRK(IiNIIS), RSKWRK(I_RKF), RSKWRK(I_RKR), ISKWRK,
     1            ISKWRK(I_NU), ISKWRK(I_NK), ISKWRK(I_IRNU),
     2            RSKWRK(I_RNU), RSKWRK(I_NCF), ISKWRK(I_IYLD),
     3            ISKWRK(I_KYLD), RSKWRK(I_PYLD), ISKWRK(I_YION),
     4            RSKWRK(I_ENRG), ISKWRK(IiNPHA), RSKWRK(I_YNCF), SDOT,
     5            SITDOT)
C
C     end of SUBROUTINE SKRATK
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRDEX (IR, ISKWRK, RSKWRK, RD)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRDEX (IR, ISKWRK, RSKWRK, RD)
C
C  Returns the perturbation factor of the IRth reaction,
C  or changes its value, depending on the sign of IR.
C
C  INPUT
C  IR        - Integer scalar, reaction index;
C              IR> 0 gets RD(I) from RSKWRK
C              IR< 0 puts RD(I) into RSKWRK
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  If IR< 0:
C  RD        - Real scalar, perturbation factor for reaction IR.
C
C  OUTPUT
C  If IR> 0:
C  RD        - Real scalar, perturbation factor for reaction IR.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*)
C
      NI = ISKWRK(IrPAR) + (ABS(IR)-1)*(NSPAR+1) + NSPAR
      IF (IR .GT. 0) THEN
         RD = RSKWRK(NI)
      ELSE
         RSKWRK(NI) = RD
      ENDIF
C
C     end of SUBROUTINE SKRDEX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRHEX (K, ISKWRK, RSKWRK, A6)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRHEX (K, ISKWRK, RSKWRK, A6)
C
C  Returns an array of the sixth thermodynamic polynomial
C  coefficients for a species, or changes their value,
C  depending on the sign of K.
C
C  INPUT
C  K         - Integer scalar, species index;
C              K > 0 gets A6(*) from RSKWRK
C              K < 0 puts A6(*) into RSKWRK
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  If K < 0:
C  A6(*)     - Integer array, the 6th thermodynamic polynomial
C              coefficients for species K, over the number of
C              temperature ranges used in fitting thermodynamic
C              properties;
C              dimension at least MAXTP-1, where MAXTP is the
C              maximum number of temperatures used in fitting the
C              thermodynamic properties of the species.
C
C  OUTPUT
C  If K > 0:
C  A6(*)     - Integer array, the 6th thermodynamic polynomial
C              coefficients for species K, over the number of
C              temperature ranges used in fitting thermodynamic
C              properties;
C              dimension at least MAXTP-1, where MAXTP is the
C              maximum number of temperatures used in fitting the
C              thermodynamic properties of the species.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), A6(*)
C
      DO 100 L = 1, MAXTP-1
         NA6 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (ABS(K)-1)*NCP2T + NCP
         IF (K .GT. 0) THEN
            A6(L) = RSKWRK(NA6)
         ELSE
            RSKWRK(NA6) = A6(L)
         ENDIF
  100 CONTINUE
C
C     end of SUBROUTINE SKRHEX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKROP  (P, T, ACT, SDEN, ISKWRK, RSKWRK, ROP)
C
C  START PROLOGUE
C
C  SUBROUTINE SKROP  (P, T, ACT, SDEN, ISKWRK, RSKWRK, ROP)
C  Returns rates of progress for the surface reactions.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  ROP(*)    - Real array, rates of progress for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec).
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), ROP(*), T(*)
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      SDTOT  = ZERO
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
         IF (ISKWRK(IiNYLD) .GT. 0) CALL SKEQY (ISKWRK, RSKWRK)
      ENDIF
C
C     initialize some local pointers in COMMON /SKPTR/
      CALL SKLOC(ISKWRK)
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, RSKWRK(I_CZ))
      CALL SKRROP (ISKWRK, RSKWRK, RSKWRK(I_SMH), T, RSKWRK(I_CZ),
     1             ACT, RSKWRK(I_KWT), ISKWRK(I_NREA), ISKWRK(I_NRPP),
     2             ISKWRK(I_NU), ISKWRK(I_NK), ISKWRK(I_IRNU),
     3             RSKWRK(I_RNU), ISKWRK(I_NUSUM), RSKWRK(I_PAR),
     4             RSKWRK(I_RPAR), ISKWRK(I_IREV), ISKWRK(I_ICOV),
     5             ISKWRK(I_KCOV), RSKWRK(I_CPAR), ISKWRK(I_ISTK),
     6             ISKWRK(I_MSTK), ISKWRK(I_IBHM), ISKWRK(I_KBHM),
     7             ISKWRK(I_IEDP), ISKWRK(I_KEDP), RSKWRK(I_PEDP),
     8             ISKWRK(I_KTFL), SDTOT, ISKWRK(I_IORD),
     9             ISKWRK(I_KORD), RSKWRK(I_ORD), ISKWRK(I_IYLD),
     *             ISKWRK(I_YION), ISKWRK(I_KYLD), RSKWRK(I_PYLD),
     1             RSKWRK(I_RKFT), RSKWRK(I_RKRT), RSKWRK(I_RKF),
     2             RSKWRK(I_RKR), RSKWRK(I_EQKC), RSKWRK(I_EQFAC),
     3             RSKWRK(I_ENRG))
C
      DO 100 I = 0, ISKWRK(IiNIIS) - 1
         ROP(I+1) = RSKWRK(I_RKF + I) - RSKWRK(I_RKR + I)
  100 CONTINUE
C
C     end of SUBROUTINE SKROP
      RETURN
      END
C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKROPK (P, T, ACT, SDEN, RKFT, RKRT, ISKWRK, RSKWRK,
     1                   ROP)
C
C  START PROLOGUE
C
C  SUBROUTINE SKROPK (P, T, ACT, SDEN, RKFT, RKRT, ISKWRK, RSKWRK,
C                     ROP)
C  Returns rates of progress for the surface reactions, given the
C  temperature-dependent rate coefficients RKFT, and RKRT.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  SDEN(*)   - Real array, site densities for the site types;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  RKFT(*)   - Real array, forward temperature-dependent rate
C              coefficient for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, depends on reaction order
C  RKRT(*)   - Real array, reverse temperature-dependent rate
C              coefficient for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, depends on reaction order
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  ROP(*)    - Real array, rates of progress for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec).
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /SKPTR/ I_CZ, I_SMH, I_KWT, I_NREA, I_NRPP, I_NU, I_NK,
     1               I_IRNU, I_RNU, I_NUSUM, I_PAR, I_RPAR, I_IREV,
     2               I_ICOV, I_KCOV, I_CPAR, I_ISTK, I_MSTK, I_IBHM,
     3               I_KBHM, I_IEDP, I_KEDP, I_PEDP, I_KTFL, I_IORD,
     4               I_KORD, I_ORD, I_IYLD, I_YION, I_KYLD, I_PYLD,
     5               I_RKFT, I_RKRT, I_RKF, I_RKR, I_EQKC, I_EQFAC,
     6               I_ENRG
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), ROP(*), T(*)
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      SDTOT  = ZERO
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
         IF (ISKWRK(IiNYLD) .GT. 0) CALL SKEQY (ISKWRK, RSKWRK)
      ENDIF
C
C     initialize some local pointers in COMMON /SKPTR/
      CALL SKLOC(ISKWRK)
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, RSKWRK(I_CZ))
      CALL SKRRPX (ISKWRK, RSKWRK, T, RSKWRK(I_CZ), ACT, ISKWRK(I_NREA),
     1             ISKWRK(I_NRPP), ISKWRK(I_NU), ISKWRK(I_NK), 
     2             ISKWRK(I_IRNU), RSKWRK(I_RNU), RSKWRK(I_PAR), 
     2             ISKWRK(I_IREV), ISKWRK(I_ICOV), ISKWRK(I_KCOV),
     3             RSKWRK(I_CPAR), ISKWRK(I_ISTK), ISKWRK(I_IBHM), 
     4             ISKWRK(I_KBHM), ISKWRK(I_KTFL), SDTOT, 
     5             ISKWRK(I_IORD), ISKWRK(I_KORD), RSKWRK(I_ORD), 
     6             RKFT, RKRT, RSKWRK(I_RKF), RSKWRK(I_RKR))
C
      DO 100 I = 0, ISKWRK(IiNIIS) - 1
         ROP(I+1) = RSKWRK(I_RKF + I) - RSKWRK(I_RKR + I)
  100 CONTINUE
C
C     end of SUBROUTINE SKROPK
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRP   (ISKWRK, RSKWRK, RU, RUC, PATM)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRP   (ISKWRK, RSKWRK, RU, RUC, PATM)
C  Returns universal gas constants and the pressure of one standard
C  atmosphere.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  RU        - Real scalar, universal gas constant.
C                 cgs units, 8.314510E7 ergs/(mole*K)
C  RUC       - Real scalar, universal gas constant used only in
C              conjuction with activation energy.
C                 preferred units, RU / 4.184 cal/(mole*K)
C  PA        - Real scalar, pressure of one standard atmosphere.
C                 cgs units, 1.01325E6 dynes/cm**2
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
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*)
C
      RU   = RSKWRK(ISKWRK(IrRU))
      RUC  = RSKWRK(ISKWRK(IrRUC))
      PATM = RSKWRK(ISKWRK(IrPATM))
C
C     end of SUBROUTINE SKRP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRPAR (ISKWRK, RSKWRK, ENRGI)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRPAR (ISKWRK, RSKWRK, ENRGI)
C
C  Allows the user to input auxiliary reaction-rate parameters for
C  special types of reactions.  The first parameter is the species (ion)
C  directed energy for ion-energy-dependent reactions.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C  ENRGI(*)  - Real array, species ion energies used in the NIIEDP
C              reactions;
C              dimension at least KKGAS, the total gas-phase species
C              count.
C              Default value stored in RSKWRK is set to 0.0 in SKINIT.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), ENRGI(*)
C
      I_ENGR = ISKWRK(IrENGI) - 1
      DO 100 K = 1, NKKGAS
         RSKWRK(I_ENGR+K) = ENRGI(K)
  100 CONTINUE
C
C     end of SUBROUTINE SKRPAR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRROP (ISKWRK, RSKWRK, SMH, T, CZ, ACT, WT, NREAC,
     1                   NRPP, NU, NUNK, IRNU, RNU, NUSUMK, PAR, RPAR,
     2                   IREV, ICOV, KCOV, CPAR, ISTK, MOTZ, IBOHM, IBK,
     3                   IEDP, KEDEP, PEDEP, KTFL, SDTOT, IORD, KORD,
     4                   RORD, IYLD, IKION, KYLD, PYLD, RKFT, RKRT, RKF,
     5                   RKR, EQKC, EQFAC, ENRGI)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRROP (ISKWRK, RSKWRK, SMH, T, CZ, ACT, WT, NREAC,
C                     NRPP, NU, NUNK, IRNU, RNU, NUSUMK, PAR, RPAR,
C                     IREV, ICOV, KCOV, CPAR, ISTK, MOTZ, IBOHM, IBK,
C                     IEDP, KEDEP, PEDEP, KTFL, SDTOT, IORD, KORD,
C                     RORD, IYLD, IKION, KYLD, PYLD, RKFT, RKRT, RKF,
C                     RKR, EQKC, EQFAC, ENRGI)
C  Returns forward and reverse rates of progress and equilibrium
C  constants for the surface reactions.
C  It is not normally called by the user application code.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C  SMH(*)    - Real array, entropy minus enthalpy for species,
C              SMH(K) = S(K)/R - H(K)/RT;
C              dimension at least KKTOT, the total species count.
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  CZ(*)     - Real array, gas-phase and surface species concentrations,
C              and bulk species activities;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS gas-phase concentrations are moles/cm**3,
C              the next KKSURF site concentrations are moles/cm**2,
C              and
C              the final KKBULK entries are bulk species activities.
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  WT(*)     - Real array, molecular weights of the species;
C              dimension at least KKTOT, the total species count.
C  NREAC(*)  - Integer array, reactant counts for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C  NRPP(*)   - Integer array, number of species (reactants+products) for
C              the surface reactions, combined with reversibility flag;
C              dimension at least IISUR, the total surface reaction
C              count.
C              NRPP > 0, NRPP species, reversible surface reaction,
C                   < 0, ABS(NRPP) species, irreversible reaction.
C  NU(*,*)   - Integer matrix, stoichiometric coefficients for
C              species in surface reactions;
C              dimension at least MAXSPR for the first and at least
C              IISUR for the second.
C              NU(N,IR) is the stoichiometric coefficient of the Nth
C              species in reaction IR, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(*,*) - Integer matrix, indices of species in reactions;
C              dimension at least MAXSP for the first, and at least
C              IISUR for the second.
C              NUNK(N,IR) is the species index for the Nth species in
C              reaction IR.
C  NUSUMK(*) - Integer array, sum of coefficients of the gas-phase
C              species in the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C  PAR(*,*)  - Real matrix, Arrhenius coefficients for reactions;
C              dimension at least NSPAR for the first, and at least
C              IISUR for the second.  For any surface reaction IR,
C              PAR(1,IR) is the pre-exponential constant
C                 cgs units, mole-cm-sec-K
C              PAR(2,IR) is the temperature dependent exponent
C                 cgs units, none
C              PAR(3,IR) is the activation energy
C                 cgs units, K
C              PAR(4,IR) is used as a perturbation factor in
C              sensitivity analyses.
C  RPAR(*,*) - Real matrix, reverse Arrhenius rate coefficients for
C              the surface reactions with explicit reverse coefficients;
C              dimension at least NSPAR for the first, and at least NREV
C              for the second, the total count of reactions with reverse
C              coefficients.
C  IREV(*)   - Integer array, reaction indices for the NREV reactions;
C              dimension at least NREV.
C  ICOV(*)   - Integer array, reaction indices for the NCOV reactions;
C              dimension at least NCOV.
C  KCOV(*)   - Integer array, coverage species indices for the NCOV
C              reactions;
C              dimension at least NCOV.
C  CPAR(*,*) - Real matrix, coverage parameters for the NCOV reactions;
C              dimension at least NSCOV for the first, the number of
C              coverage parameters allowed, and at least NCOV for the
C              second, the total coverage reaction count.
C  ISTK(*)   - Integer array, reaction indices for the NSTK reactions;
C              dimension at least NSTK, the total sticking reaction
C              count.
C  MOTZ(*)  - Integer array, 0/1 flag for Motz-Wise correction of
C              sticking reaction rate;
C              dimension at least NSTK, the total sticking reaction
C              count.
C  IBOHM(*)  - Integer array, reaction indices for the Bohmreactions;
C              dimension at least NBOHM, the total Bohm reaction count.
C  IBK(*)    - Integer array, species indices for the Bohm reactions;
C              dimension at least NBOHM, the total Bohm reaction count.
C  IEDP(*)   - Integer array, reaction indices for the ion-energy-
C              dependent reactions; dimension at least NIIEDP.
C  KEDEP(*)  - Integer array, species indices for the ion-energy-
C              dependent reactions; dimension at least NIIEDP.
C  PEDEP(*,*)- Real matrix, parameters for the ion-energy-dependent
C              reactions;
C              dimension at least NEDPAR for the first, the number of
C              parameters required, and at least NIIEDP for the second.
C  KTFL(*)   - Integer array, indices into the temperature(s) for
C              gas-phase species;
C              dimension at least KKGAS, the total gas-phase species
C              count.
C              Default value stored in ISKWRK is set to 1 in SKINIT.
C  SDTOT     - Real scalar, the sum of the densities of the phases.
C  IORD(*)   - Integer array, reaction indices for the NORD reactions;
C              dimension at least NORD.
C              IORD(N) is the index of the Nth change-order reaction.
C  KORD(*,*) - Integer matrix, species indices for the order changes in
C              the NORD reactions; dimension at least MXORD for the
C              first and at least NORD for the second.
C              KORD(L,N) is the species index for the Lth order change
C              in the Nth change-order reaction.
C              KORD < 0 indicates change in forward order;
C              KORD > 0 indicates change in reverse order.
C  RORD(*,*) - Real matrix, order values for the NORD reactions;
C              dimension at least MXORD for the first and at least NORD
C              for the second.
C              RORD(L,N) is the order for the Lth order change in the
C              Nth change-order reaction.
C  IYLD(*)   - Integer array, reaction indices for yield-modify
C              reactions;
C              dimension at least NIIYLD.
C  IKION(*)  - Integer array, ion species indices of yield-modify
C              reactions;
C              dimension at least NIIYLD.
C  KYLD(*,*) - Integer matrix, yield flags for the species in the yield-
C              modify reactions;
C              dimension at least MAXSPR for the first, and at least
C              NIIYLD for the second.
C  PYLD(*,*) - Real matrix, the parameters for the yield-modify
C              reactions;
C              dimension at least NYPAR for the first, and at least
C              NIIYLD for the second.
C  EQFAC(*)  - Real array, factors for equilibrium constants for surface
C              reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C  ENRGI(*)  - Real array, species ion energies used in the NIIEDP
C              reactions;
C              dimension at least KKGAS, the total gas-phase species
C              count.
C
C  OUTPUT
C  RKF(*)    - Real array, forward rates of progress for the surface
C              reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec)
C  RKR(*)    - Real array, reverse rates of progress for the surface
C              reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec)
C  EQKC(*)   - Real array, equilibrium constants in concentration
C              units for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C      PARAMETER (ZERO = 0.0, ONE = 1.0)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /MACHN/ SMALL,BIG,EXPARG
      DIMENSION ISKWRK(*), RSKWRK(*), SMH(*), CZ(*), ACT(*), WT(*),
     1          NREAC(*), NRPP(*), NU(MAXSPR,*), NUNK(MAXSPR,*),
     2          NUSUMK(*), PAR(NSPAR+1,*), RPAR(NSPAR+1,*), IREV(*),
     3          ICOV(*), KCOV(*), CPAR(NSCOV,*), ISTK(*),
     4          RKRT(*), RKFT(*), RKF(*), RKR(*), EQKC(*), EQFAC(*),
     5          IRNU(*), RNU(MAXSPR,*), IORD(*), KORD(MAXORD,*),
     6          RORD(MAXORD,*), T(*), IBOHM(*), IBK(*), IEDP(*),
     7          KEDEP(*), PEDEP(NEDPAR,*), KTFL(*), ENRGI(*),
     8          IYLD(*), IKION(*), KYLD(MAXSPR,*), PYLD(NYPAR,*),
     9          MOTZ(*)
C
      CALL SKRRPT (ISKWRK, RSKWRK, SMH, T, WT, NREAC,
     1             NRPP, NU, NUNK, IRNU, RNU, NUSUMK, PAR, RPAR,
     2             IREV, ISTK, MOTZ, IBOHM, IBK,
     3             IEDP, KEDEP, PEDEP, KTFL, IORD, KORD,
     4             RORD, IYLD, IKION, KYLD, PYLD, RKFT, RKRT, 
     5             EQKC, EQFAC, ENRGI)
C
      CALL SKRRPX (ISKWRK, RSKWRK, T, CZ, ACT, NREAC,
     1             NRPP, NU, NUNK, IRNU, RNU, PAR, 
     2             IREV, ICOV, KCOV, CPAR, ISTK, IBOHM, IBK,
     3             KTFL, SDTOT, IORD, KORD,
     4             RORD, RKFT, RKRT, RKF, RKR)
C
C     end of SUBROUTINE SKRROP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRRPT (ISKWRK, RSKWRK, SMH, T, WT, NREAC,
     1                   NRPP, NU, NUNK, IRNU, RNU, NUSUMK, PAR, RPAR,
     2                   IREV, ISTK, MOTZ, IBOHM, IBK,
     3                   IEDP, KEDEP, PEDEP, KTFL, IORD, KORD,
     4                   RORD, IYLD, IKION, KYLD, PYLD, RKFT, RKRT, 
     5                   EQKC, EQFAC, ENRGI)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRRPT (ISKWRK, RSKWRK, SMH, T, WT, NREAC,
C                     NRPP, NU, NUNK, IRNU, RNU, NUSUMK, PAR, RPAR,
C                     IREV, ISTK, MOTZ, IBOHM, IBK,
C                     IEDP, KEDEP, PEDEP, KTFL, IORD, KORD,
C                     RORD, IYLD, IKION, KYLD, PYLD, RKFT, RKRT, 
C                     EQKC, EQFAC, ENRGI)
C  Returns forward and reverse rates of progress and equilibrium
C  constants for the surface reactions.
C  It is not normally called by the user application code.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C  SMH(*)    - Real array, entropy minus enthalpy for species,
C              SMH(K) = S(K)/R - H(K)/RT;
C              dimension at least KKTOT, the total species count.
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  WT(*)     - Real array, molecular weights of the species;
C              dimension at least KKTOT, the total species count.
C  NREAC(*)  - Integer array, reactant counts for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C  NRPP(*)   - Integer array, number of species (reactants+products) for
C              the surface reactions, combined with reversibility flag;
C              dimension at least IISUR, the total surface reaction
C              count.
C              NRPP > 0, NRPP species, reversible surface reaction,
C                   < 0, ABS(NRPP) species, irreversible reaction.
C  NU(*,*)   - Integer matrix, stoichiometric coefficients for
C              species in surface reactions;
C              dimension at least MAXSPR for the first and at least
C              IISUR for the second.
C              NU(N,IR) is the stoichiometric coefficient of the Nth
C              species in reaction IR, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(*,*) - Integer matrix, indices of species in reactions;
C              dimension at least MAXSP for the first, and at least
C              IISUR for the second.
C              NUNK(N,IR) is the species index for the Nth species in
C              reaction IR.
C  NUSUMK(*) - Integer array, sum of coefficients of the gas-phase
C              species in the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C  PAR(*,*)  - Real matrix, Arrhenius coefficients for reactions;
C              dimension at least NSPAR for the first, and at least
C              IISUR for the second.  For any surface reaction IR,
C              PAR(1,IR) is the pre-exponential constant
C                 cgs units, mole-cm-sec-K
C              PAR(2,IR) is the temperature dependent exponent
C                 cgs units, none
C              PAR(3,IR) is the activation energy
C                 cgs units, K
C              PAR(4,IR) is used as a perturbation factor in
C              sensitivity analyses.
C  RPAR(*,*) - Real matrix, reverse Arrhenius rate coefficients for
C              the surface reactions with explicit reverse coefficients;
C              dimension at least NSPAR for the first, and at least NREV
C              for the second, the total count of reactions with reverse
C              coefficients.
C  IREV(*)   - Integer array, reaction indices for the NREV reactions;
C              dimension at least NREV.
C  ISTK(*)   - Integer array, reaction indices for the NSTK reactions;
C              dimension at least NSTK, the total sticking reaction
C              count.
C  MOTZ(*)  - Integer array, 0/1 flag for Motz-Wise correction of
C              sticking reaction rate;
C              dimension at least NSTK, the total sticking reaction
C              count.
C  IBOHM(*)  - Integer array, reaction indices for the Bohmreactions;
C              dimension at least NBOHM, the total Bohm reaction count.
C  IBK(*)    - Integer array, species indices for the Bohm reactions;
C              dimension at least NBOHM, the total Bohm reaction count.
C  IEDP(*)   - Integer array, reaction indices for the ion-energy-
C              dependent reactions; dimension at least NIIEDP.
C  KEDEP(*)  - Integer array, species indices for the ion-energy-
C              dependent reactions; dimension at least NIIEDP.
C  PEDEP(*,*)- Real matrix, parameters for the ion-energy-dependent
C              reactions;
C              dimension at least NEDPAR for the first, the number of
C              parameters required, and at least NIIEDP for the second.
C  KTFL(*)   - Integer array, indices into the temperature(s) for
C              gas-phase species;
C              dimension at least KKGAS, the total gas-phase species
C              count.
C              Default value stored in ISKWRK is set to 1 in SKINIT.
C  IORD(*)   - Integer array, reaction indices for the NORD reactions;
C              dimension at least NORD.
C              IORD(N) is the index of the Nth change-order reaction.
C  KORD(*,*) - Integer matrix, species indices for the order changes in
C              the NORD reactions; dimension at least MXORD for the
C              first and at least NORD for the second.
C              KORD(L,N) is the species index for the Lth order change
C              in the Nth change-order reaction.
C              KORD < 0 indicates change in forward order;
C              KORD > 0 indicates change in reverse order.
C  RORD(*,*) - Real matrix, order values for the NORD reactions;
C              dimension at least MXORD for the first and at least NORD
C              for the second.
C              RORD(L,N) is the order for the Lth order change in the
C              Nth change-order reaction.
C  IYLD(*)   - Integer array, reaction indices for yield-modify
C              reactions;
C              dimension at least NIIYLD.
C  IKION(*)  - Integer array, ion species indices of yield-modify
C              reactions;
C              dimension at least NIIYLD.
C  KYLD(*,*) - Integer matrix, yield flags for the species in the yield-
C              modify reactions;
C              dimension at least MAXSPR for the first, and at least
C              NIIYLD for the second.
C  PYLD(*,*) - Real matrix, the parameters for the yield-modify
C              reactions;
C              dimension at least NYPAR for the first, and at least
C              NIIYLD for the second.
C  EQFAC(*)  - Real array, factors for equilibrium constants for surface
C              reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C  ENRGI(*)  - Real array, species ion energies used in the NIIEDP
C              reactions;
C              dimension at least KKGAS, the total gas-phase species
C              count.
C
C  OUTPUT
C  RKFT(*)   - Real array, forward temperature-dependent rate
C              coefficient for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, depends on reaction order
C  RKRT(*)   - Real array, reverse temperature-dependent rate
C              coefficient for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, depends on reaction order
C  EQKC(*)   - Real array, equilibrium constants in concentration
C              units for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C      PARAMETER (ZERO = 0.0, ONE = 1.0)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /MACHN/ SMALL,BIG,EXPARG
      DIMENSION ISKWRK(*), RSKWRK(*), SMH(*), WT(*),
     1          NREAC(*), NRPP(*), NU(MAXSPR,*), NUNK(MAXSPR,*),
     2          NUSUMK(*), PAR(NSPAR+1,*), RPAR(NSPAR+1,*), IREV(*),
     3          ISTK(*),
     4          RKRT(*), RKFT(*), EQKC(*), EQFAC(*),
     5          IRNU(*), RNU(MAXSPR,*), IORD(*), KORD(MAXORD,*),
     6          RORD(MAXORD,*), T(*), IBOHM(*), IBK(*), IEDP(*),
     7          KEDEP(*), PEDEP(NEDPAR,*), KTFL(*), ENRGI(*),
     8          IYLD(*), IKION(*), KYLD(MAXSPR,*), PYLD(NYPAR,*),
     9          MOTZ(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      IISUR = ISKWRK(IiNIIS)
      IF (IISUR .LE. 0) RETURN
C
      ALOGT = LOG(T(1))
      DO 20 I = 1, IISUR
C        Forward rate constant for reaction I
C        ("standard" Arrhenius expression)
         RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*ALOGT - PAR(3,I)/T(1))
   20 CONTINUE
C
      NSTK = ISKWRK(IiNSTK)
      NRNU = ISKWRK(IiNRNU)
      NORD = ISKWRK(IiNORD)
      KKGAS = NKKGAS
      KKSUR = ISKWRK(IiKSUR)
      SKMIN = RSKWRK(ISKWRK(IrSKMN))
C
      DO 50 N = 1, NSTK
C        sticking coefficient reaction
         I = ISTK(N)
C        real coefficient reaction?
         ISRNU = CKLKUP (I, IRNU, NRNU)
C
C        find a gas-phase species with order approx. 1.0;
C        interpreter error-checkint has ensured that it is there
         KGAS = 0
C        change-order reaction?
         ISORD = CKLKUP(I, IORD, NORD)
         IF (ISORD .GT. 0) THEN
            L = ISORD
            DO 52 J = 1, MAXORD
C              latest species orders
               K = KORD(J,L)
               IF (K.LT.0 .AND. ABS(K).LE.KKGAS .AND.
     1             RORD(J,L).LE.1.0+SKMIN .AND.
     2             RORD(J,L).GE.1.0-SKMIN) KGAS = ABS(K)
   52       CONTINUE
         ENDIF
C
         IF (KGAS .EQ. 0) THEN
            DO 54 L = 1, NREAC(I)
               K = NUNK(L,I)
               IF (K.GT.0 .AND. K.LE.KKGAS) THEN
                  IF (ISRNU .GT. 0) THEN
                     IF (ABS(RNU(L,ISRNU)).LE.1.0+SKMIN .AND.
     1                   ABS(RNU(L,ISRNU)).GE.1.0-SKMIN) KGAS = K
                  ELSE
                     IF (ABS(NU(L,I)) .EQ. 1) KGAS = K
                  ENDIF
               ENDIF
   54       CONTINUE
         ENDIF
C
C        temperature for gas may have been set by application code;
C        default is T(1)
C
         TEMP = T(KTFL(KGAS))
         RKFT(I) = PAR(1,I) * EXP(PAR(2,I) * LOG(TEMP) - PAR(3,I)/TEMP)
         RKFT(I) = MIN(RKFT(I), ONE)
C
         IF (MOTZ(N) .EQ. 0) THEN
            RKFT(I) = RKFT(I) * 3637.6011 * SQRT(TEMP / WT(KGAS))
         ELSE
            RKFT(I) = RKFT(I) * 3637.6011 * SQRT(TEMP / WT(KGAS))
     1                / (1.0 - RKFT(I) / 2.0)
         ENDIF
C
   50 CONTINUE
C
      NBOHM = ISKWRK(IiNBHM)
      KEL   = ISKWRK(IiELEC)
      DO 70 N = 1, NBOHM
C        BOHM reactions
         I = IBOHM(N)
         K = IBK(N)
         TEMP = T(KTFL(KEL))
C        Bohm velocity:  sqrt(kTe/mi) x gamma
         RKFT(I) = PAR(1,I) * 9117.76 * SQRT(TEMP/WT(K))
   70 CONTINUE
C
      DO 80 N = 1, ISKWRK(IiNEDP)
C        Reactions with ion-energy dependence
         I = IEDP(N)
C        directed energy of the ion in Kelvin
         TION = ENRGI(KEDEP(N))
C        energy threshold for this reaction
         ETH = PEDEP(1,N)
C
         IF (TION .GE. ETH) THEN
C           apply the ion-energy dependence of the reaction rate;
C           first exponent,
            A   = PEDEP(2,N)
C           second exponent,
            B   = PEDEP(3,N)
            RKFT(I) = RKFT(I) * (TION**A - ETH**A)**B
         ELSE
            RKFT(I) = 0.0
         ENDIF
   80 CONTINUE
C
      CALL SKSMH (T(1), ISKWRK, RSKWRK, SMH)
      PFAC = RSKWRK(ISKWRK(IrPATM)) / (RSKWRK(ISKWRK(IrRU))*T(1))
      DO 100 I = 1, IISUR
C
C        Gibbs free energy and equilibrium constant for every
C        reaction
        EQKC(I) = 1.0
C
        SUMSMH = NU(1,I)*SMH(NUNK(1,I)) + NU(7,I)*SMH(NUNK(7,I))
        IF (NUNK(2,I) .NE. 0) THEN
          SUMSMH = SUMSMH + NU(2,I)*SMH(NUNK(2,I))
          IF (NUNK(3,I) .NE. 0) THEN
            SUMSMH = SUMSMH + NU(3,I)*SMH(NUNK(3,I))
            IF (NUNK(4,I) .NE. 0) THEN
              SUMSMH = SUMSMH + NU(4,I)*SMH(NUNK(4,I))
              IF (NUNK(5,I) .NE. 0) THEN
                SUMSMH = SUMSMH + NU(5,I)*SMH(NUNK(5,I))
                IF (NUNK(6,I) .NE. 0) THEN
                  SUMSMH = SUMSMH + NU(6,I)*SMH(NUNK(6,I))
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        IF (NUNK(8,I) .NE. 0) THEN
           SUMSMH = SUMSMH + NU(8,I)*SMH(NUNK(8,I))
           IF (NUNK(9,I) .NE. 0) THEN
              SUMSMH = SUMSMH + NU(9,I)*SMH(NUNK(9,I))
              IF (NUNK(10,I) .NE. 0) THEN
                 SUMSMH = SUMSMH + NU(10,I)*SMH(NUNK(10,I))
                 IF (NUNK(11,I) .NE. 0) THEN
                    SUMSMH = SUMSMH + NU(11,I)*SMH(NUNK(11,I))
                    IF (NUNK(12,I) .NE. 0) THEN
                       SUMSMH = SUMSMH + NU(12,I)*SMH(NUNK(12,I))
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
C
C        Equilibrium constant; EQKC(I)=1.0 if SUMSMH=0.0
C
         IF (SUMSMH .NE. 0.0)
     1   EQKC(I) = EXP(MIN(SUMSMH,EXPARG)) * PFAC**NUSUMK(I) * EQFAC(I)
C
  100 CONTINUE
C
      DO 125 N = 1, NRNU
C
C        Gibbs free energy and equilibrium constant for those reactions
C        with real stoichiometric coefficients
C
         I = IRNU(N)
         EQKC(I) = 1.0
         SUMSMH = 0.0
         RNUSUM = 0.0
C
         DO 110 L = 1, NREAC(I)
            NK = NUNK(L,I)
            SUMSMH = SUMSMH + RNU(L,N) * SMH(NK)
            IF (NK .LE. KKGAS) RNUSUM = RNUSUM + RNU(L,N)
  110    CONTINUE
C
         DO 120 L = 7, 12
            NK = NUNK(L,I)
            IF (NK .EQ. 0) GO TO 121
            SUMSMH = SUMSMH + RNU(L,N) * SMH(NK)
            IF (NK .LE. KKGAS) RNUSUM = RNUSUM + RNU(L,N)
  120    CONTINUE
  121    CONTINUE
C
         PFRNU = EXP (RNUSUM * LOG(PFAC))
C        Equilibrium constant = 1.0 if SUMSMH=0.0
         IF (SUMSMH .NE. 0.0)
     1   EQKC(I) = EXP(MIN(SUMSMH,EXPARG)) * PFRNU * EQFAC(I)
  125 CONTINUE
C
      DO 150 I = 1, IISUR
C
C        Default reverse reaction rate constant for every reaction;
C        0.0 for irreversible reactions, else RKFT / MAX(EQKC,SMALL)
C        = RKFT / MAX(EQK,SMALL)
C
         IF (NRPP(I) .GT. 0) THEN
           RKRT(I) = RKFT(I) / MAX(EQKC(I),SMALL)
         ELSE
           RKRT(I) = ZERO
         ENDIF
  150 CONTINUE
C
      NREV = ISKWRK(IiNREV)
      DO 200 N = 1, NREV
C
C        Reverse reaction rate constant if explicit reverse rate
C        parameters given
C
         I = IREV(N)
         RKRT(I) = RPAR(1,N) * EXP(RPAR(2,N)*ALOGT-RPAR(3,N)/T(1))
         IF (RKRT(I) .NE. ZERO) EQKC(I) = RKFT(I) / RKRT(I)
  200 CONTINUE
C
C
C     end of SUBROUTINE SKRRPT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRRPX (ISKWRK, RSKWRK, T, CZ, ACT, NREAC,
     1                   NRPP, NU, NUNK, IRNU, RNU, PAR, 
     2                   IREV, ICOV, KCOV, CPAR, ISTK, IBOHM, IBK,
     3                   KTFL, SDTOT, IORD, KORD,
     4                   RORD, RKFT, RKRT, RKF, RKR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRRPX (ISKWRK, RSKWRK, T, CZ, ACT, NREAC,
C                     NRPP, NU, NUNK, IRNU, RNU, PAR, 
C                     IREV, ICOV, KCOV, CPAR, ISTK, IBOHM, IBK,
C                     KTFL, SDTOT, IORD, KORD,
C                     RORD, RKFT, RKRT, RKF, RKR)
C  Returns forward and reverse rates of progress and equilibrium
C  constants for the surface reactions.
C  It is not normally called by the user application code.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  CZ(*)     - Real array, gas-phase and surface species concentrations,
C              and bulk species activities;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS gas-phase concentrations are moles/cm**3,
C              the next KKSURF site concentrations are moles/cm**2,
C              and
C              the final KKBULK entries are bulk species activities.
C  ACT(*)    - Real array, activities of the species;
C              dimension at least KKTOT, the total species count.
C              The first KKGAS activities are mole fractions,
C              the next KKSURF activities are site fractions
C                 (species density normalized by the site density;
C                  surface concentration in moles/cm**2 is
C                  ACT(K)*SITE_DENSITY / # sites per species), and
C              the next KKBULK activities for bulk phase species
C              should be from 0 to 1, and should sum to 1 for each
C              phase.
C  NREAC(*)  - Integer array, reactant counts for the surface reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C  NRPP(*)   - Integer array, number of species (reactants+products) for
C              the surface reactions, combined with reversibility flag;
C              dimension at least IISUR, the total surface reaction
C              count.
C              NRPP > 0, NRPP species, reversible surface reaction,
C                   < 0, ABS(NRPP) species, irreversible reaction.
C  NU(*,*)   - Integer matrix, stoichiometric coefficients for
C              species in surface reactions;
C              dimension at least MAXSPR for the first and at least
C              IISUR for the second.
C              NU(N,IR) is the stoichiometric coefficient of the Nth
C              species in reaction IR, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(*,*) - Integer matrix, indices of species in reactions;
C              dimension at least MAXSP for the first, and at least
C              IISUR for the second.
C              NUNK(N,IR) is the species index for the Nth species in
C              reaction IR.
C  PAR(*,*)  - Real matrix, Arrhenius coefficients for reactions;
C              dimension at least NSPAR for the first, and at least
C              IISUR for the second.  For any surface reaction IR,
C              PAR(1,IR) is the pre-exponential constant
C                 cgs units, mole-cm-sec-K
C              PAR(2,IR) is the temperature dependent exponent
C                 cgs units, none
C              PAR(3,IR) is the activation energy
C                 cgs units, K
C              PAR(4,IR) is used as a perturbation factor in
C              sensitivity analyses.
C  IREV(*)   - Integer array, reaction indices for the NREV reactions;
C              dimension at least NREV.
C  ICOV(*)   - Integer array, reaction indices for the NCOV reactions;
C              dimension at least NCOV.
C  KCOV(*)   - Integer array, coverage species indices for the NCOV
C              reactions;
C              dimension at least NCOV.
C  CPAR(*,*) - Real matrix, coverage parameters for the NCOV reactions;
C              dimension at least NSCOV for the first, the number of
C              coverage parameters allowed, and at least NCOV for the
C              second, the total coverage reaction count.
C  ISTK(*)   - Integer array, reaction indices for the NSTK reactions;
C              dimension at least NSTK, the total sticking reaction
C              count.
C  IBOHM(*)  - Integer array, reaction indices for the Bohmreactions;
C              dimension at least NBOHM, the total Bohm reaction count.
C  IBK(*)    - Integer array, species indices for the Bohm reactions;
C              dimension at least NBOHM, the total Bohm reaction count.
C  KTFL(*)   - Integer array, indices into the temperature(s) for
C              gas-phase species;
C              dimension at least KKGAS, the total gas-phase species
C              count.
C              Default value stored in ISKWRK is set to 1 in SKINIT.
C  SDTOT     - Real scalar, the sum of the densities of the phases.
C  IORD(*)   - Integer array, reaction indices for the NORD reactions;
C              dimension at least NORD.
C              IORD(N) is the index of the Nth change-order reaction.
C  KORD(*,*) - Integer matrix, species indices for the order changes in
C              the NORD reactions; dimension at least MXORD for the
C              first and at least NORD for the second.
C              KORD(L,N) is the species index for the Lth order change
C              in the Nth change-order reaction.
C              KORD < 0 indicates change in forward order;
C              KORD > 0 indicates change in reverse order.
C  RORD(*,*) - Real matrix, order values for the NORD reactions;
C              dimension at least MXORD for the first and at least NORD
C              for the second.
C              RORD(L,N) is the order for the Lth order change in the
C              Nth change-order reaction.
C
C  OUTPUT
C  RKF(*)    - Real array, forward rates of progress for the surface
C              reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec)
C  RKR(*)    - Real array, reverse rates of progress for the surface
C              reactions;
C              dimension at least IISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C      PARAMETER (ZERO = 0.0, ONE = 1.0)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      COMMON /MACHN/ SMALL,BIG,EXPARG
      DIMENSION ISKWRK(*), RSKWRK(*), CZ(*), ACT(*), 
     1          NREAC(*), NRPP(*), NU(MAXSPR,*), NUNK(MAXSPR,*),
     2          PAR(NSPAR+1,*), IREV(*),
     3          ICOV(*), KCOV(*), CPAR(NSCOV,*), ISTK(*),
     4          RKRT(*), RKFT(*), RKF(*), RKR(*), 
     5          IRNU(*), RNU(MAXSPR,*), IORD(*), KORD(MAXORD,*),
     6          RORD(MAXORD,*), T(*), IBOHM(*), IBK(*), 
     7          KTFL(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      IISUR = ISKWRK(IiNIIS)
      IF (IISUR .LE. 0) RETURN
C
      NSTK = ISKWRK(IiNSTK)
      NRNU = ISKWRK(IiNRNU)
      NORD = ISKWRK(IiNORD)
      KKGAS = NKKGAS
      KKSUR = ISKWRK(IiKSUR)
      SKMIN = RSKWRK(ISKWRK(IrSKMN))
C
      DO 50 N = 1, NSTK
C        sticking coefficient reaction
         I = ISTK(N)
C        real coefficient reaction?
         ISRNU = CKLKUP (I, IRNU, NRNU)
C        change-order reaction?
         ISORD = CKLKUP(I, IORD, NORD)
C
         IF (SDTOT .GT. 0.0) THEN
C
            SDTR = 1.0
            IF (ISORD .GT. 0) THEN
C              change-of-order reaction
               RSUM = 0.0
               DO 46 L = 1, MAXORD
                  K = KORD(L, ISORD)
                  IF (K .LT. 0) THEN
                     IF (ABS(K).GT.KKGAS .AND. ABS(K).LE.KKGAS+KKSUR)
     1                  RSUM = RSUM + ABS(RORD(L,ISORD))
                  ENDIF
 46            CONTINUE
               IF (RSUM .GT. SKMIN) THEN
                  SDTR = SDTOT**RSUM
               ENDIF
            ELSEIF (ISRNU .GT. 0) THEN
C              real stoichiometry coefficients
               RSUM = 0.0
               DO 42 L = 1, NREAC(I)
C                 include site species
                  K = NUNK(L,I)
                  IF (K.GT.KKGAS .AND. K.LE.KKGAS+KKSUR)
     1               RSUM = RSUM + ABS(RNU(L,ISRNU))
   42          CONTINUE
               IF (RSUM .GT. SKMIN) THEN
                  SDTR = SDTOT**RSUM
               ENDIF
            ELSE
C              integer stoichiometry coefficients
               NSUM = 0
               DO 25 L = 1, NREAC(I)
C                 include site species
                  K = NUNK(L,I)
                  IF (K.GT.KKGAS .AND. K.LE.KKGAS+KKSUR)
     1               NSUM = NSUM + ABS(NU(L,I))
   25          CONTINUE
               IF (NSUM .GT. 0) SDTR = SDTOT**NSUM
            ENDIF
            RKFT(I) = RKFT(I) / SDTR
            RKRT(I) = RKRT(I) / SDTR
         ENDIF
   50 CONTINUE
C
      NCOV = ISKWRK(IiNCOV)
      DO 60 N = 1, NCOV
C        coverage reactions
         I = ICOV(N)
         COEF = ACT(KCOV(N))
         IF (COEF .GT. ZERO) THEN
            COVFAC =  10.0**(CPAR(1,N) * COEF)
     1                * COEF**CPAR(2,N)
     2                * EXP(-CPAR(3, N) * COEF / T(1))
            RKFT(I) = RKFT(I) * COVFAC
            RKRT(I) = RKRT(I) * COVFAC
         ENDIF
   60 CONTINUE
C
      NBOHM = ISKWRK(IiNBHM)
      KEL   = ISKWRK(IiELEC)
      DO 70 N = 1, NBOHM
C        BOHM reactions
         I = IBOHM(N)
C
         IF (SDTOT .EQ. 0.0 .OR. KKSUR.EQ.0) GO TO 70
C
         ISORD = CKLKUP(I, IORD, NORD)
         ISRNU = CKLKUP (I, IRNU, NRNU)
         SDTR = 1.0
         IF (ISORD .GT. 0) THEN
C            change-of-order reaction
            RSUM = 0.0
            DO 62 L = 1, MAXORD
               K = KORD(L, ISORD)
C              Find reactant site species only
               IF (K .LT. 0) THEN
                  IF (ABS(K).GT.KKGAS .AND. ABS(K).LE.KKGAS+KKSUR)
     1                 RSUM = RSUM + ABS(RORD(L,ISORD))
               ENDIF
 62         CONTINUE
            IF (RSUM .GT. SKMIN) THEN
               SDTR = SDTOT**RSUM
            ENDIF
         ELSEIF (ISRNU .GT. 0) THEN
C           real stoichiometry coefficients
            RSUM = 0.0
            DO 63 M = 1, NREAC(I)
C              Find site species only
               NK = NUNK(M,I)
               IF (NK.GT.KKGAS .AND. NK.LE.KKGAS+KKSUR)
     1             RSUM = RSUM + ABS(RNU(M,ISRNU))
   63       CONTINUE
            IF (RSUM .GT. SKMIN) THEN
               SDTR = SDTOT**RSUM
            ENDIF
         ELSE
C           integer stoichiometry coefficients
            NSUM = 0
            DO 65 L = 1, NREAC(I)
C              Find site species only
               NK = NUNK(L,I)
               IF (NK.GT.KKGAS .AND. NK.LE.KKGAS+KKSUR)
     1            NSUM = NSUM + ABS(NU(L,I))
   65       CONTINUE
            IF (NSUM .GT. 0) SDTR = SDTOT**NSUM
         ENDIF
         RKFT(I) = RKFT(I) / SDTR
         RKRT(I) = ZERO
   70 CONTINUE
C
      DO 250 I = 1, IISUR
C
C        all integer stoichiometric coefficients are zero if this
C        reaction has real coefficients
         IF (NU(1,I) .EQ. 0) GO TO 250
C
C        Forward rates from rate constant and reactant concentrations
         RKF(I) = RKFT(I) * CZ(NUNK(1,I))**IABS(NU(1,I)) * PAR(4,I)
         IF (NUNK(2,I) .NE. 0) THEN
            RKF(I) = RKF(I) * CZ(NUNK(2,I))**IABS(NU(2,I))
            IF (NUNK(3,I) .NE. 0) THEN
               RKF(I) = RKF(I) * CZ(NUNK(3,I))**IABS(NU(3,I))
               IF (NUNK(4,I) .NE. 0) THEN
                  RKF(I) = RKF(I) * CZ(NUNK(4,I))**IABS(NU(4,I))
                  IF (NUNK(5,I) .NE. 0) THEN
                     RKF(I) = RKF(I) * CZ(NUNK(5,I))**IABS(NU(5,I))
                     IF (NUNK(6,I) .NE. 0) THEN
                        RKF(I) = RKF(I) * CZ(NUNK(6,I))**IABS(NU(6,I))
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
         IF (RKRT(I) .EQ. ZERO) THEN
C           Reverse reaction rate constant is zero for irreversible
C           reaction (most often case)
            RKR(I) = ZERO
            GO TO 250
         ENDIF
C
C        Reverse rates from rate constant and product concentrations
         RKR(I) = RKRT(I) * CZ(NUNK(7,I))**NU(7,I) * PAR(4,I)
         IF (NUNK(8,I) .NE. 0) THEN
            RKR(I) = RKR(I) * CZ(NUNK(8,I))**NU(8,I)
            IF (NUNK(9,I) .NE. 0) THEN
               RKR(I) = RKR(I) * CZ(NUNK(9,I))**NU(9,I)
               IF (NUNK(10,I) .NE. 0) THEN
                  RKR(I) = RKR(I) * CZ(NUNK(10,I))**NU(10,I)
                  IF (NUNK(11,I) .NE. 0) THEN
                     RKR(I) = RKR(I) * CZ(NUNK(11,I))**NU(11,I)
                     IF (NUNK(12,I) .NE. 0) THEN
                        RKR(I) = RKR(I) * CZ(NUNK(12,I))**NU(12,I)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  250 CONTINUE
C
      DO 350 N = 1, NRNU
C
C        Rates for the reactions with real stoichiometric coefficients
C        (the integer coefficients were zero)
         I = IRNU(N)
C
         RKF(I) = RKFT(I) * MAX(ZERO,CZ(NUNK(1,I)))**ABS(RNU(1,N)) 
     1            * PAR(4,I)
         IF (NUNK(2,I) .NE. 0) THEN
            RKF(I) = RKF(I) * MAX(ZERO,CZ(NUNK(2,I)))**ABS(RNU(2,N))
            IF (NUNK(3,I) .NE. 0) THEN
               RKF(I) = RKF(I) * MAX(ZERO,CZ(NUNK(3,I)))**ABS(RNU(3,N))
               IF (NUNK(4,I) .NE. 0) THEN
                  RKF(I) = RKF(I) 
     1                     * MAX(ZERO,CZ(NUNK(4,I)))**ABS(RNU(4,N))
                  IF (NUNK(5,I) .NE. 0) THEN
                     RKF(I) = RKF(I) *
     1                        MAX(ZERO,CZ(NUNK(5,I)))**ABS(RNU(5,N))
                     IF (NUNK(6,I) .NE. 0) THEN
                        RKF(I) = RKF(I) *
     1                           MAX(ZERO,CZ(NUNK(6,I)))**ABS(RNU(6,N))
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
         IF (RKRT(I) .EQ. ZERO) THEN
C           irreversible reaction
            RKR(I) = ZERO
            GO TO 350
         ENDIF
C
          RKR(I) = RKRT(I) * MAX(ZERO,CZ(NUNK(7,I)))**RNU(7,N) 
     1                     * PAR(4,I)
          IF (NUNK(8,I) .NE. 0) THEN
             RKR(I) = RKR(I) * MAX(ZERO,CZ(NUNK(8,I)))**RNU(8,N)
            IF (NUNK(9,I) .NE. 0) THEN
               RKR(I) = RKR(I) * MAX(ZERO,CZ(NUNK(9,I)))**RNU(9,N)
               IF (NUNK(10,I) .NE. 0) THEN
                  RKR(I) = RKR(I) * MAX(ZERO,CZ(NUNK(10,I)))**RNU(10,N)
                  IF (NUNK(11,I) .NE. 0) THEN
                     RKR(I) = RKR(I) * 
     1                        MAX(ZERO,CZ(NUNK(11,I)))**RNU(11,N)
                     IF (NUNK(12,I) .NE. 0) THEN
                        RKR(I) = RKR(I) * 
     1                           MAX(ZERO,CZ(NUNK(12,I)))**RNU(12,N)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  350 CONTINUE
C
      DO 450 N = 1, NORD
C
C        Recalculate rates for the change-order reactions;
C        may have different orders than their stoichiometry
C
         I = IORD(N)
         RKF(I) = RKFT(I)
         RKR(I) = RKRT(I)
C
         DO 440 L = 1, MAXORD
            K = KORD(L,N)
            ORD = RORD(L,N)
            IF (K .LT. 0) THEN
               RKF(I) = RKF(I) * MAX(ZERO,CZ(-K))**ORD
            ELSEIF (K .GT. 0) THEN
               RKR(I) = RKR(I) * MAX(ZERO,CZ(K))**ORD
            ENDIF
  440    CONTINUE
  450 CONTINUE
C
      NBOHM = ISKWRK(IiNBHM)
      DO 550 N = 1, NBOHM
C
C        Recalculation for Bohm reactions
C        only the ion IBK concentration used in the rate
         I = IBOHM(N)
         RKR(I) = ZERO
         RKF(I) = RKFT(I) * CZ(IBK(N))
C
         ISORD = CKLKUP(I, IORD, NORD)
         ISRNU = CKLKUP(I, IRNU, NRNU)
         IF (ISORD .GT. 0) THEN
C           Find site or bulk species
            DO 480 L = 1, MAXORD
               NK = KORD(L,ISORD)
               ORD = RORD(L,ISORD)
               IF (NK .LT. 0 .AND. ABS(NK) .GT. KKGAS) THEN
                  RKF(I) = RKF(I) * MAX(ZERO,CZ(-NK))**ORD
               ENDIF
 480        CONTINUE
         ELSEIF (ISRNU .GT. 0) THEN
            DO 485 NS = 1, NREAC(I)
C              find site or bulk species
               NK = NUNK(M,I)
               IF (NK .GT. KKGAS) THEN
                  RKF(I) = RKF(I) * MAX(ZERO,CZ(NK))**ABS(RNU(NS,ISRNU))
               ENDIF
 485        CONTINUE
         ELSE
            DO 500 NS = 1, NREAC(I)
C              Find site or bulk species
               NK = NUNK(NS,I)
               IF (NK.GT.KKGAS)
     1              RKF(I) = RKF(I) * CZ(NK)**IABS(NU(NS,I))
 500        CONTINUE
         ENDIF
  550 CONTINUE
C
C     end of SUBROUTINE SKRRPX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSAVE (LOUT, LSAVE, ISKWRK, RSKWRK, CSKWRK)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSAVE (LOUT, LSAVE, ISKWRK, RSKWRK, CSKWRK)
C  Writes to a binary file information about a Surface Chemkin
C  linkfile, pointers for the Surface Chemkin Library, and
C  Surface Chemkin work arrays.
C
C  INPUT
C  LOUT      - Integer scalar, formatted output file unit number.
C  LSAVE     - Integer scalar, unformatted output file unit number.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C  CSKWRK(*) - Character string workspace array; dimension at least
C              LENCSK.
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
      INCLUDE 'skstrt.h'
      COMMON /SKCONS/ PREC, FILVER, KERR
C
      DIMENSION ISKWRK(*), RSKWRK(*)
      CHARACTER*(*) CSKWRK(*)
      CHARACTER*16 FILVER, PREC
      LOGICAL KERR
C
      LENI = ISKWRK(IiLENI)
      LENR = ISKWRK(IiLENR)
      LENC = ISKWRK(IiLENC)
      WRITE (LSAVE, ERR=999) FILVER, PREC, LENI, LENR, LENC,
C
C     Include file for CHEMKIN-III sklib.f, dated:  March 1, 1966
C
C     Integer constants
C
     1   MAXSPR, NELEM, NKKGAS, NSPAR, NSCOV, NEDPAR, NYPAR, MAXORD,
     2   MAXTP, NCP,    NCP1,  NCP2,   NCP2T,
C
C     ISKWRK pointers to integer variables
C
     3   IiLENI, IiLENR, IiLENC, IiKSUR, IiKBLK, IiKTOT, IiNPHA, IiFSUR,
     4   IiLSUR, IiNSUR, IiFBLK, IiLBLK, IiNBLK, IiNIIS, IiNCOV, IiNREV,
     5   IiNSTK, IiNCON, IiNBHM, IiNRNU, IiNORD, IiMOTZ, IiELEC, IiNYLD,
C
C     ISKWRK pointers to integer arrays
C
     6   IiPKST, IiPKND, IiPTOT, IiKPHS, IiKCHG, IiKCMP, IiNSCV, IiKNT,
     7   IiNRPP, IiNREA, IiNUNK, IiNU,   IiNSUM, IiICOV, IiKCOV, IiIREV,
     8   IiISTK, IiMSTK, IiIBHM, IiKBHM, IiIRNU, IiIORD, IiKORD, IiIONS,
     9   IiKION, IiKTFL, IiNEDP, IiIEDP, IiKEDP, IiIYLD, IiYION, IiKYLD,
C
C     ISKWRK pointers to real variables
C
     *   IrSKMN, IrPATM, IrRU,   IrRUC,
C
C     ISKWRK pointers to real arrays
C
     1   IrSDEN, IrKTMP, IrKTHM, IrKDEN, IrAWT,  IrKWT,  IrPAR, IrKCOV,
     2   IrRPAR, IrEQ,   IrRNU,  IrNCF,  IrKORD, IrKFT,  IrKRT, IrKT1,
     3   IrKT2,  IrPT1,  IrIT1,  IrIT2,  IrIT3,  IrPEDP, IrENGI,
     4   IrPYLD, IrYNCF,
C
C     ISKWRK pointers to character arrays
C
     4   IcENAM, IcKNAM, IcMNAM, IcPNAM
C
C     END include file for sklib.f
C
      WRITE (LSAVE, ERR=999) (ISKWRK(L), L = 1, LENI)
      WRITE (LSAVE, ERR=999) (RSKWRK(L), L = 1, LENR)
      WRITE (LSAVE, ERR=999) (CSKWRK(L), L = 1, LENC)
C
      GO TO 1000
  999 CONTINUE
      WRITE (LOUT, *)
     1 ' Error writing Surface binary file information...'
C
C     end of SUBROUTINE SKSAVE
 1000 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSDEN (ISKWRK, RSKWRK, SDEN0)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSDEN (ISKWRK, RSKWRK, SDEN0)
C  Returns a real array of standard-state phase densities as given
C  on input to the interpreter.
C
C  INPUT
C  RSKWRK(*) - Real workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  SDEN0(*)  - Real array; standard-state densities for the
C              site types, AS READ BY THE INTERPRETER;
C              dimension at least NPHASE, the total phase count,
C              but the subroutine only uses site phase entries,
C              NFSURF <= N <= NLSURF.
C                 cgs units, moles/cm**2.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), SDEN0(*)
C
      IF (ISKWRK(IiNSUR) .LE. 0) RETURN
C
      DO 100 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
         SDEN0(N) = RSKWRK(ISKWRK(IrSDEN) + N - 1)
100   CONTINUE
C
C     end of SUBROUTINE SKSDEN
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C                                                                      C
      SUBROUTINE SKSDOT(NIISUR, RKF, RKR, ISKWRK, NU, NUNK, IRNU, RNU,
     1                  RNCF, IYLD, KYLD, PYLD, KION, ENRGI, NPHASE,
     2                  YNCF, SDOT, SITDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSDOT(NIISUR, RKF, RKR, ISKWRK, NU, NUNK, IRNU, RNU,
C                    RNCF, IYLD, KYLD, PYLD, KION, ENRGI, NPHASE,
C                    YNCF, SDOT, SITDOT)
C
C  Returns the production rates of the species, SDOT, and sites,
C  SITDOT, given the rates of progess, RKF, RKR, of the forward
C  and reverse reaction surface reaction rates.
C  This subroutine is not normally called at the user level.
C
C
C  INPUT
C  NIISUR    - Integer scalar, the total surface reaction count,
C              and the first dimension of the matrix RNCF.
C  RKF(I)    - Real array, forward rates of progress for the surface
C              reactions;
C              dimension at least NIISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec)
C  RKR(I)    - Real array, reverse rates of progress for the surface
C              reactions;
C              dimension at least NIISUR, the total surface reaction
C              count.
C                 cgs units, moles/(cm**2*sec)
C  ISKWRK(*) - Integer workspace array; dimension at least LENIXK.
C  NU(N,I)   - Integer matrix, stoichiometric coefficients for
C              species, NUNK(N,I), in surface reaction, I;
C              dimension at least MAXSPR for the first and at least
C              IISUR for the second.
C              NU(N,IR) is the stoichiometric coefficient of the Nth
C              species in reaction IR, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(N,I) - Integer matrix, index, N, of species in reaction, I;
C              dimension at least MAXSP for the first, and at least
C              IISUR for the second.
C              NUNK(N,IR) is the species index for the Nth species in
C              reaction IR.
C  IRNU(L)   - Integer array, reaction indices for the surface reactions
C              with real stoichiometric coefficients;
C              dimension at least NRNU, the total count of real-
C              stoichiometry surface reactions.
C              Note:  these reactions have their NU rows exactly zero.
C  RNU(N,L)  - Real matrix, stoichiometric coefficients for species,
C              NUNK(N,I), in the real-stoichiometry surface reaction,
C              IRNU(L)
C              Dimension at least MAXSPR for the first, the maximum
C              number of species allowed in a reaction, and at least
C              NRNU for the second, the total count of real-
C              stoichiometry surface reactions.
C  RNCF(N,I) - Real matrix, the change in site concentrations for
C              surface phase N due to surface reaction I;
C              dimension at least NIISUR for the first dimension,
C              the total surface reaction count, and at least NPHASE
C              for the second, the total phase count.
C  IYLD(*)   -
C  KYLD(*,*) -
C  PYLD(*)   -
C  KION(*)   -
C  ENRGI(*)  -
C  NPHASE    - Number of phases in the current material type
C  YNCF(*,*) -
C
C  OUTPUT
C  SDOT(*)   - Real array, production rates of the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, moles/(cm**2*sec)
C              for 1,KKGAS, the production rate of gas-phase species,
C              for KKGAS+1,KKGAS+KKSUR, the production rate of surface
C              species,
C              for KKGAS+KKSUR+1,KKTOT, the production rate of bulk
C              species.
C  SITDOT(*) - Real array, production rates of the surface phases;
C              dimension at least NPHASE, the total phase count, but
C              subroutine only calculates entries for site phases.
C                 cgs units, moles/(cm**2*sec)
C
C  Note:       SDOT and SITDOT are expected to be zero on input to this
C              routine.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), NUNK(MAXSPR,*), NU(MAXSPR,*), IRNU(*),
     1          IYLD(*), KION(*), KYLD(MAXSPR,*), RKF(NIISUR),
     2          RKR(NIISUR), RNCF(NPHASE,*), RNU(MAXSPR,*),
     3          PYLD(NYPAR,*), ENRGI(*), YNCF(NPHASE,*), SDOT(*),
     4          SITDOT(*)
C
      EXTERNAL CKLKUP
      INTEGER  CKLKUP
C
C       Loop over surface reactions gathering production rates
C        (Largest cost in algorithm so optimize this loop)
C
      DO 100 I = 1, NIISUR
C
C           NU(1,I) is only equal to zero if reaction I has 
C           real stoichiometric coefficients. It will be taken
C           care of in the following loop.
C
         IF (NU(1,I) .EQ. 0) GO TO 100
C
C        Store the rate of progress for the reaction
         ROP = RKF(I) - RKR(I)
C
C        This algorithm depends on the fact that forward
C        species indices start at 1 and reverse species indices
C        start at 7, and that there is at least one reactant
C        and one product so NUNK(1,I) and NU(1,I), NUNK(7,I) and
C        NU(7,I) are nonzero.
C
         K = NUNK(1,I)
         SDOT(K) = SDOT(K) + ROP * NU(1,I)
         K = NUNK(2,I)
         IF (K .GT. 0) THEN
            SDOT(K) = SDOT(K) + ROP * NU(2,I)
            K = NUNK(3,I)
            IF (K .GT. 0) THEN
               SDOT(K) = SDOT(K) + ROP * NU(3,I)
               K = NUNK(4,I)
               IF (K .GT. 0) THEN
                  SDOT(K) = SDOT(K) + ROP * NU(4,I)
                   K = NUNK(5,I)
                   IF (K .GT. 0) THEN
                      SDOT(K) = SDOT(K) + ROP * NU(5,I)
                      K = NUNK(6,I)
                      IF (K .GT. 0) SDOT(K) = SDOT(K) + ROP*NU(6,I)
                   ENDIF
               ENDIF
            ENDIF
         ENDIF
C
         K = NUNK(7,I)
         SDOT(K) = SDOT(K) + ROP * NU(7,I)
         K = NUNK(8,I)
         IF (K .GT. 0) THEN
           SDOT(K) = SDOT(K) + ROP * NU(8,I)
           K = NUNK(9,I)
           IF (K .GT. 0) THEN
              SDOT(K) = SDOT(K) + ROP * NU(9,I)
              K = NUNK(10,I)
              IF (K .GT. 0) THEN
                 SDOT(K) = SDOT(K) + ROP * NU(10,I)
                 K = NUNK(11,I)
                 IF (K .GT. 0) THEN
                    SDOT(K) = SDOT(K) + ROP * NU(11,I)
                    K = NUNK(12,I)
                    IF (K .GT. 0) SDOT(K) = SDOT(K) + ROP*NU(12,I)
                 ENDIF
              ENDIF
           ENDIF
         ENDIF
  100 CONTINUE
C
C       Calculate the rate of change of the total site density for all
C       surface phases.
C
      NKKSUR = ISKWRK(IiKSUR)
      IF (NKKSUR .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SITDOT(N) = ZERO
            DO 55 I = 1, NIISUR
               SITDOT(N) = SITDOT(N) + RNCF(N,I)*(RKF(I)-RKR(I))
   55       CONTINUE
   60    CONTINUE
      ENDIF
C
C       Now fix up the results for those reactions which have
C       real stoichiometric coefficients. The reason why we
C       are not double counting reactions with real coefficients
C       between loops 100 and 200 is because the NU(N,I)'s are set to
C       zero for reactions with real stoichiometric coefficients.
C
      NIIRNU = ISKWRK(IiNRNU)
      DO 200 L = 1, NIIRNU
         I = IRNU(L)
         ROP = RKF(I) - RKR(I)
         K = NUNK(1,I)
         SDOT(K) = SDOT(K) + ROP * RNU(1,L)
         K = NUNK(2,I)
         IF (K .GT. 0) THEN
            SDOT(K) = SDOT(K) + ROP * RNU(2,L)
            K = NUNK(3,I)
            IF (K .GT. 0) THEN
               SDOT(K) = SDOT(K) + ROP * RNU(3,L)
               K = NUNK(4,I)
               IF (K .GT. 0) THEN
                  SDOT(K) = SDOT(K) + ROP * RNU(4,L)
                  K = NUNK(5,I)
                  IF (K .GT. 0) THEN
                     SDOT(K) = SDOT(K) + ROP * RNU(5,L)
                     K = NUNK(6,I)
                     IF (K .GT. 0) SDOT(K) = SDOT(K) + ROP*RNU(6,L)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
         K = NUNK(7,I)
         SDOT(K) = SDOT(K) + ROP * RNU(7,L)
         K = NUNK(8,I)
         IF (K .GT. 0) THEN
            SDOT(K) = SDOT(K) + ROP * RNU(8,L)
            K = NUNK(9,I)
            IF (K .GT. 0) THEN
               SDOT(K) = SDOT(K) + ROP * RNU(9,L)
               K = NUNK(10,I)
               IF (K .GT. 0) THEN
                  SDOT(K) = SDOT(K) + ROP * RNU(10,L)
                  K = NUNK(11,I)
                  IF (K .GT. 0) THEN
                     SDOT(K) = SDOT(K) + ROP * RNU(11,L)
                     K = NUNK(12,I)
                     IF (K .GT. 0) SDOT(K) = SDOT(K) + ROP * RNU(12,L)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
 200  CONTINUE
C
      NIIYLD = ISKWRK(IiNYLD)
      IF (NIIYLD .LE. 0) RETURN
C
C  Variable yield surface reactions require special treatment
C
      DO 250 L = 1, NIIYLD
C
C          Determine the rxn which has the variable yield and
C          store its rate of progress
C
         I = IYLD(L)
         ROP = RKF(I) - RKR(I)
C
C          Determine the multiplicative scaling factor to
C          multiply all stoichiometric coefficients with
C          for the reaction
C             - mult ROP into scale factor as well
C
         EI    = ENRGI(KION(L))
         ETH   = PYLD(2,L)
         IF (EI .LT. ETH) THEN
           SCALE = - ROP
         ELSE
           ASCAL = PYLD(1,L)
           A     = PYLD(3,L)
           B     = PYLD(4,L)
           SCALE = ROP * (ASCAL * (EI**A - ETH**A) ** B - 1.0)
         ENDIF
C
C          Look to see if this rxn, I, is also a real coefficient
C          reactions -> this could be simplified by ALWAYS
C          assuming that these reactions are real coefficients!
C
         ISRNU = CKLKUP(I, IRNU, NIIRNU)
C
         DO 230 N = 1, MAXSPR
            K = NUNK(N,I)
            IF (K.GT.0 .AND. KYLD(N,L).GT.0) THEN
               IF (ISRNU .GT. 0) THEN
                  COEF = RNU(N,ISRNU)
               ELSE
                  COEF = NU(N,I)
               ENDIF
               SDOT(K) = SDOT(K) + COEF * SCALE
            ENDIF
 230     CONTINUE
  250 CONTINUE
C
      IF (NKKSUR .GT. 0) THEN
         DO 240 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            DO 235 L = 1, NIIYLD
               I = IYLD(L)
               SITDOT(N) = SITDOT(N) + YNCF(N,L) * (RKF(I) - RKR(I))
  235       CONTINUE
  240    CONTINUE
      ENDIF
C
C     end of SUBROUTINE SKSDOT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSMH  (T, ISKWRK, RSKWRK, SMH)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSMH  (T, ISKWRK, RSKWRK, SMH)
C  Returns the array of dimensionless entropies minus enthalpies
C  for the species.  It is normally not called directly by the user.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  SMH(*)    - Real array, dimensionless entropies minus enthalpies
C              for the species;
C              dimension at least KKTOT, the total species count.
C              SMH(K) = S(K)/R - H(K)/RT.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), SMH(*), T(*)
      SAVE TN1, TNHALF, TNLOG, TN2, TN3, TN4
      DATA TN1/0.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNHALF = TN1/2
         TNLOG = LOG(TN1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2/6
         TN3 = TN3/12
         TN4 = TN4/20
      ENDIF
      I_KTFL = ISKWRK(IiKTFL)
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 250 K = 0, NKM1
         IF (K.LT.NKKGAS .AND. ISKWRK(I_KTFL+K).GT.1) THEN
            TK1 = T(ISKWRK(I_KTFL+K))
            TKLOG = LOG(TK1)
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TKHALF = TK1 / 2
            TK2 = TK2 / 6
            TK3 = TK3 / 12
            TK4 = TK4 / 20
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TKHALF= TNHALF
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ISKWRK(ISKWRK(IiKNT) + K) - 1
C        location of FIRST set of thermodynamic coefficients
         NA1 = ISKWRK(IrKTHM) + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = ISKWRK(IrKTMP) + K*MAXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RSKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SMH(K+1) = RSKWRK(NA1) * (TKLOG - 1.0)
     1            + RSKWRK(NA1+1) * TKHALF + RSKWRK(NA1+2) * TK2
     2            + RSKWRK(NA1+3) * TK3    + RSKWRK(NA1+4) * TK4
     3            - RSKWRK(NA1 + NCP1 - 1) / TK1
     4            + RSKWRK(NA1 + NCP2 - 1)
 250  CONTINUE
C
C     end of SUBROUTINE SKSMH
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSML  (T, ISKWRK, RSKWRK, SML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSML  (T, ISKWRK, RSKWRK, SML)
C  Returns an array of the standard state entropies in molar units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  SML(*)    - Real array, standard state entropies for the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, ergs/(mole*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), SML(*), T(*)
      SAVE TN1, TNLOG, TN2, TN3, TN4
      DATA TN1/0.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNLOG = LOG(TN1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2/2
         TN3 = TN3/3
         TN4 = TN4/4
      ENDIF
      I_KTFL = ISKWRK(IiKTFL)
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 250 K = 0, NKM1
         IF (K.LT.NKKGAS .AND. ISKWRK(I_KTFL+K).GT.1) THEN
            TK1 = T(ISKWRK(I_KTFL+K))
            TKLOG = LOG(TK1)
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TK2=TN2
            TK3=TN3
            TK4=TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ISKWRK(ISKWRK(IiKNT)+K) - 1
C        location of FIRST set of thermodynamic coefficients
         NA1 = ISKWRK(IrKTHM) + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = ISKWRK(IrKTMP) + K*MAXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RSKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SML(K+1) = RSKWRK(ISKWRK(IrRU)) * (RSKWRK(NA1) * TKLOG
     1            + RSKWRK(NA1+1) * TK1 + RSKWRK(NA1+2) * TK2
     2            + RSKWRK(NA1+3) * TK3 + RSKWRK(NA1+4) * TK4
     3            + RSKWRK(NA1 + NCP2 - 1))
250   CONTINUE
C
C     end of SUBROUTINE SKSML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSMS  (T, ISKWRK, RSKWRK, SMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSMS  (T, ISKWRK, RSKWRK, SMS)
C  Returns an array of the standard state entropies in mass units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  SMS(*)    - Real array, standard state entropies for the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, ergs/(gm*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), SMS(*), T(*)
      SAVE TN1, TNLOG, TN2, TN3, TN4
      DATA TN1/0.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNLOG = LOG(TN1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
      ENDIF
      I_KTFL = ISKWRK(IiKTFL)
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 250 K = 0, NKM1
         IF (K.LT.NKKGAS .AND. ISKWRK(I_KTFL+K).GT.1) THEN
            TK1 = T(ISKWRK(I_KTFL+K))
            TKLOG = LOG(TK1)
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TK2=TN2
            TK3=TN3
            TK4=TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ISKWRK(ISKWRK(IiKNT)+K) - 1
C        location of FIRST set of thermodynamic coefficients
         NA1 = ISKWRK(IrKTHM) + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = ISKWRK(IrKTMP) + K*MAXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RSKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SMS(K+1) = RSKWRK(ISKWRK(IrRU)) * (RSKWRK(NA1) * TKLOG
     1            + RSKWRK(NA1+1) * TK1 + RSKWRK(NA1+2) * TK2
     2            + RSKWRK(NA1+3) * TK3 + RSKWRK(NA1+4) * TK4
     3            + RSKWRK(NA1 + NCP2 - 1)) / RSKWRK(ISKWRK(IrKWT)+K)
250   CONTINUE
C
C     end of SUBROUTINE SKSMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSNUM (LINE, NEXP, LOUT, KNAM, KKTOT, PNAM, NNPHAS,
     1                   KKPHAS, KNUM, NT, NVAL, RVAL, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSNUM (LINE, NEXP, LOUT, KNAM, KKTOT, PNAM, NNPHAS,
C                     KKPHAS, KNUM, NT, NVAL, RVAL, KERR)
C  This subroutine is used to read a format-free input line of
C  combined alphanumerical data.  It can be used to parse an input
C  character string, LINE, which may be composed of several blank-
C  delimited substrings.  This subroutine assumes that the first
C  substring in LINE is the name of a species in the Surface Chemkin
C  mechanism.  If the species name is not unique within the Surface
C  Chemkin mechanism, the phase of the species should be input
C  immediately after the species name, delimited by slashes.
C  Upon return from the subroutine, KNUM returns the index position
C  of the species within the Surface Chemkin linkfile.  If the
C  species name is not unique, KNUM returns the first position and
C  NT returns the number of the times the species occurs within the
C  linkfile.  If the species name is not found, or there is a
C  syntax error, on return, KNUM=0, NT=0, and KERR=.TRUE.
C  The substrings in LINE following the first are expected to
C  represent numbers.  They are converted into floating point
C  values and stored in the output vector, RVAL(*).  Upon input,
C  NEXP is equal to the number of values expected to be found.
C  If NEXP numbers are not found, KERR will be set to .TRUE. on
C  return from the subroutine.
C
C     Example input:
C             LINE     = GA(S)/BULK1/ 1.2
C             NEXP     = 1, the number of values expected
C             LOUT     = 6, a logical unit number on which to write
C                        diagnostic messages
C             KNAM(*)  = Array of character species names
C             KKTOT    = Total number of species
C             PNAM(*)  = Array of character phase names
C             NNPHAS   = Total number of phases
C             KKPHAS(*)= Index array of the number of species in the
C                        phases
C     Output:
C             KNUM     = The index number of the species which
C                        has the name "GA(S)" and resides in phase
C                        "BULK1"
C             NT       = 1, if there is only one species GA(S)
C                        in phase BULK1
C             NVAL     = 1, the number of values found in LINE
C                        following the species name
C             RVAL(1)  = 1.200E+00, the substring converted to a
C                        real number
C             KERR     = .FALSE.
C
C  INPUT
C  LINE      - Character string; length depends on calling routine.
C  NEXP      - Integer scalar, number of values to be found in LINE.
C              If NEXP < 0, then IABS(NEXP) values are expected, but
C              it is not an error condition if less values are found.
C  LOUT      - Integer scalar, formatted output file unit number.
C  KNAM(*)   - Character string array, species names;
C              dimension at least KKTOT, the total species count.
C  KKTOT     - Integer scalar, the total species count.
C  PNAM(*)   - Character string array, phase names;
C              dimension at least NNPHAS, the total phase count.
C  NNPHAS    - Integer scalar, the total phase count.
C  KKPHAS(*) - Integer array, total species counts for the phases;
C              dimension at least NNPHAS, the total phase count.
C
C  OUTPUT
C  KNUM      - Integer scalar, species index if the species name appears
C              in LINE.
C  NT        - Integer scalar, number of times the species name occurs
C              in the linkfile.
C  NVAL      - Integer scalar, number of value character strings found
C              in LINE.
C  RVAL(*)   - Real array, real values for their character strings in
C              LINE;
C              dimension at least NEXP, the number of values expected.
C  KERR      - Logical, syntax or dimensioning error flag;
C
C  END PROLOG
C     A !' will comment out a line, or remainder of the line.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION RVAL(*), KKPHAS(*)
      CHARACTER*(*) KNAM(*), PNAM(*), LINE
      CHARACTER*80 ISTR
      LOGICAL KERR
      INTEGER CKFRCH, CKLSCH, CKSLEN
      EXTERNAL CKFRCH, CKLSCH, CKSLEN
C
      KNUM = 0
      NVAL = 0
      KERR = .FALSE.
      NT   = 0
C
      ILEN = MIN (CKSLEN(LINE), CKLSCH(LINE))
      IF (ILEN .LE. 0) RETURN
C
      IF (INDEX(LINE,'/') .GT. 0) THEN
         I1 = INDEX (LINE, '/')
         I2 = I1 + INDEX(LINE(I1+1:), '/')
      ELSE
         I2 = CKFRCH(LINE) + INDEX(LINE(CKFRCH(LINE):),' ') - 2
      ENDIF
C
      I1 = CKFRCH(LINE)
      IF (I2 .GE. I1) CALL SKPCMP (LINE(I1:I2), KNAM, KKTOT, PNAM,
     1                             NNPHAS, KKPHAS, KNUM, NT)
C
      IF (KNUM .EQ. 0) THEN
         ISTR = ' '
         ISTR = LINE(1:ILEN)
         IF (NEXP .GT. 0) THEN
            WRITE (LOUT, '(A)')
     1      ' Error in SKSNUM...'//ISTR(1:ILEN)//' not found...'
            KERR = .TRUE.
         ENDIF
      ENDIF
      IF (NEXP .NE. 0)
     1   CALL CKXNUM (LINE(I2+1:), NEXP, LOUT, NVAL, RVAL, KERR)
C
C     end of SUBROUTINE SKSNUM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSOR  (T, ISKWRK, RSKWRK, SOR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSOR  (T, ISKWRK, RSKWRK, SOR)
C  Returns an array of the nondimensional entropies.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  SOR(*)    - Real array, nondimensional entropies for the species;
C              dimension at least KKTOT, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), SOR(*), T(*)
      SAVE TN1, TNLOG, TN2, TN3, TN4
      DATA TN1/0.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNLOG = LOG(TN1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
      ENDIF
      I_KTFL = ISKWRK(IiKTFL)
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 250 K = 0, NKM1
         IF (K.LT.NKKGAS .AND. ISKWRK(I_KTFL+K).GT.1) THEN
            TK1 = T(ISKWRK(I_KTFL+K))
            TKLOG = LOG(TK1)
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ISKWRK(ISKWRK(IiKNT)+K) - 1
C        location of FIRST set of thermodynamic coefficients
         NA1 = ISKWRK(IrKTHM) + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = ISKWRK(IrKTMP) + K*MAXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RSKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SOR(K+1) = RSKWRK(NA1) * TKLOG
     1            + RSKWRK(NA1+1) * TK1 + RSKWRK(NA1+2) * TK2
     2            + RSKWRK(NA1+3) * TK3 + RSKWRK(NA1+4) * TK4
     3            + RSKWRK(NA1 + NCP2 - 1)
250   CONTINUE
C
C     end of SUBROUTINE SKSOR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSYME (ISKWRK, CSKWRK, LOUT, ENAM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSYME (ISKWRK, CSKWRK, LOUT, ENAM, KERR)
C  Returns a character string array of element names.
C
C  INPUT
C  CSKWRK(*) - Character string workspace array; dimension at least
C              LENCSK.
C  LOUT      - Integer scalar, formatted output file unit number.
C
C  OUTPUT
C  ENAM(*)   - Character string array, element names;
C              dimension at least NELEM, the total element count.
C  KERR      - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*)
      CHARACTER*(*) ENAM(*), CSKWRK(*)
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      KERR = .FALSE.
      ILEN = LEN(ENAM(1))
      I_ENAM = ISKWRK(IcENAM) - 1
C
      DO 100 N = 1, NELEM
         ENAM(N) = ' '
         LT = CKLSCH(CSKWRK(I_ENAM + N))
         IF (LT .GT. ILEN) THEN
            KERR = .TRUE.
         ELSE
            ENAM(N) = CSKWRK(I_ENAM + N)
         ENDIF
  100 CONTINUE
C
C     end of SUBROUTINE SKSYME
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSYMM (ISKWRK, CSKWRK, LOUT, MATNAM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSYMM (ISKWRK, CSKWRK, LOUT, MATNAM, KERR)
C  Returns the character string name of a material.
C
C  INPUT
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  CSKWRK(*) - Character string workspace array; dimension at least
C              LENCSK.
C  LOUT      - Integer scalar, formatted output file unit number.
C
C  OUTPUT
C  MATNAM    - Character string, material name.
C  KERR      - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*)
      CHARACTER*(*) CSKWRK(*), MATNAM
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      KERR = .FALSE.
      ILEN = LEN(MATNAM)
      MATNAM = ' '
      I_MNAM = ISKWRK(IcMNAM)
      LT = CKLSCH(CSKWRK(I_MNAM))
      IF (LT .GT. ILEN) THEN
         KERR = .TRUE.
         WRITE (LOUT,*)
     1   ' Error in SKSYMM...character string length too small '
      ELSE
         MATNAM = CSKWRK(I_MNAM)
      ENDIF
C
C     end of SUBROUTINE SKSYMM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSYMP (ISKWRK, CSKWRK, LOUT, PNAM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSYMP (ISKWRK, CSKWRK, LOUT, PNAM, KERR)
C  Returns a character string array of phase names.
C
C  INPUT
C  CSKWRK(*) - Character string workspace array; dimension at least
C              LENCSK.
C  LOUT      - Integer scalar, formatted output file unit number.
C
C  OUTPUT
C  PNAM(*)   - Character string array, phase names;
C              dimension at least NNPHAS, the total phase count.
C  KERR      - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*)
      CHARACTER*(*) PNAM(*), CSKWRK(*)
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      KERR = .FALSE.
      ILEN = LEN(PNAM(1))
      I_PNAM = ISKWRK(IcPNAM) - 1
      DO 100 N = 1, ISKWRK(IiNPHA)
         PNAM(N) = ' '
         LT = CKLSCH(CSKWRK(I_PNAM + N))
         IF (LT .GT. ILEN) THEN
            KERR = .TRUE.
            WRITE (LOUT,*)
     1      ' Error in SKSYMP...character string length too small '
         ELSE
            PNAM(N) = CSKWRK(I_PNAM + N)
         ENDIF
  100 CONTINUE
C
C     end of SUBROUTINE SKSYM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSYMR (IR, LOUT, ISKWRK, RSKWRK, CSKWRK, LT, RNAM,
     1                   KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSYMR (IR, LOUT, ISKWRK, RSKWRK, CSKWRK, LT, RNAM,
C                     KERR)
C  Returns the character string representation of reaction IR.
C
C  INPUT
C  IR        - Integer scalar, reaction index.
C  LOUT      - Integer scalar, formatted output file unit number.
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C  CSKWRK(*) - Character string workspace array; dimension at least
C              LENCSK.
C
C  OUTPUT
C  LT        - Integer scalar, nunber of non-blank characters in the
C              reaction string.
C  RNAM      - Character string, representation of reaction.
C  KERR      - Logical, character length error flag.
C
C  END PROLOGUE
C
      INCLUDE 'skstrt.h'
      INTEGER       IR, ISKWRK(*), LT, J, NP, ILEN, ISRNU, 
     $              I_RNU, I_NU, I_NK, NU, I_KNAM, ISIGN, I_PNAM,
     $              LOUT, KFIRST, KLAST, KNUM, LK, LL, N, K, NSTP,
     $              NS, NSTRT, NK

      CHARACTER RNAM*(*), CSKWRK(*)*(*), IDUM*16, ISTR*255
      LOGICAL   KERR
C*****precision > double
      DOUBLE PRECISION RSKWRK(*), RNU
C*****END precision > double
C*****precision > single
C       REAL RSKWRK(*), RNU
C*****END precision > single
C
      INTEGER  CKLSCH, CKLKUP
      EXTERNAL CKLSCH, CKLKUP, CKR2CH, CKI2CH
C
      RNAM = ' '
      ISTR = ' '
      ILEN = LEN(ISTR)
      IF (IR .LE. 0 .OR. IR .GT. ISKWRK(IiNIIS)) THEN
        LT = 0
        KERR = .TRUE.
        RETURN
      ENDIF
C
C       Check to see whether the IRth reaction has real stoichiometric
C       coefficients
C
      ISRNU = CKLKUP(IR, ISKWRK(ISKWRK(IiIRNU)), ISKWRK(IiNRNU))
C
C       Calculate the beginning of the stoich coefficients in the
C       workspace vector for the IRth reaction
C 
      IF (ISRNU .GT. 0)
     1    I_RNU = ISKWRK(IrRNU) + (ISRNU-1)*MAXSPR - 1

      I_NU   = ISKWRK(IiNU) + (IR-1)*MAXSPR - 1
      I_NK   = ISKWRK(IiNUNK) + (IR-1)*MAXSPR - 1
      I_KNAM = ISKWRK(IcKNAM)
      I_PNAM = ISKWRK(IcPNAM)
C
C       Loop over the reactants and then the products
C
      DO 100 J = 1, 2
         NS = 0
         NSTRT = (J-1) * MAXSPR/2 + 1
         NSTP  = NSTRT + MAXSPR/2 - 1
C
C          Calculate the correct multiplying sign for products and
C          reactants. All stoich coefficients should turn out positive
C          except for weird cases.
C
         IF (J .EQ. 1) THEN
           ISIGN = -1
         ELSE
           ISIGN = 1
         ENDIF
C
C         Loop over the species on the reactant or product side
C
         DO 50 N = NSTRT, NSTP
            NS = NS + 1
            NU = ISIGN * ISKWRK(I_NU + N)
            K  = ISKWRK(I_NK + N)
C
C Get the stoichiometric coefficient and store it in the character 
C array, IDUM
C
            IDUM = ' '
            IF (ISRNU .EQ. 0) THEN
              IF (NU .EQ. 0) GOTO 50
              IF (NU .GT. 1) THEN
C                 convert NU to character string
                CALL CKI2CH (NU, IDUM, LK, KERR)
                IF (KERR) GOTO 500
              ENDIF
            ELSE
              RNU = ISIGN * RSKWRK(I_RNU + N)
              IF (RNU .EQ. 0.0) GOTO 50
              NU = INT(RNU)
              IF ((RNU-NU) .EQ. 0.0) THEN
                IF (NU .NE. 1) THEN
                  CALL CKI2CH(NU, IDUM, LK, KERR)
                  IF (KERR) GO TO 500
                ENDIF
              ELSE
C                  convert RNU to character string
                CALL CKR2CH(RNU, IDUM, LK, KERR)
                IF (KERR) GOTO 500
              ENDIF
            ENDIF
C
C Add the Plus symbol between species to the string
C
            IF (NS .GT. 1) THEN
              LL = CKLSCH(ISTR)
              IF (LL+1 .GT. ILEN) GOTO 400
              ISTR(LL+1:LL+1) = '+'
            ENDIF
C
C Add the stoichiometric coefficient to the string
C           
            IF (IDUM .NE. ' ') THEN
              LL = CKLSCH(ISTR)
              IF (LL+LK .GT. ILEN) GO TO 400
              ISTR(LL+1:LL+LK) = IDUM(1:LK)
            ENDIF
C
C Add the name of the species to the string
C
            LK = CKLSCH(CSKWRK(I_KNAM + K - 1))
            LL = CKLSCH(ISTR)
            IF (LL+LK .GT. ILEN) GOTO 400
            ISTR(LL+1:LL+LK) = CSKWRK(I_KNAM+K-1)(1:LK)
C
C Does the species name occur more than once?
C    If it does, we have to add the string /PSYM/ to the name of the
C    species.
C
            CALL SKCOMP(CSKWRK(I_KNAM+K-1), CSKWRK(I_KNAM),
     $                  ISKWRK(IiKTOT), KNUM, NK)
            IF (NK .GT. 1) THEN
C
C             Judge what phase the species occurs by checking 
C             the boundary species numbers.
C             When the phase has been found, add the phase name
C             to the string.
C
              DO 40 NP = 1, ISKWRK(IiNPHA)
                KFIRST = ISKWRK(ISKWRK(IiPKST) + NP - 1)
                KLAST  = ISKWRK(ISKWRK(IiPKND) + NP - 1)
                IF (K.GE.KFIRST .AND. K.LE.KLAST) THEN
                  LK = CKLSCH(CSKWRK(I_PNAM + NP - 1))
                  LL = CKLSCH(ISTR)
                  IF (LL+LK+2 .GT. ILEN) GOTO 400
                  ISTR(LL+1:LL+LK+2) =
     1              '/' // CSKWRK(I_PNAM + NP - 1)(1:LK) // '/'
                ENDIF
   40         CONTINUE
            ENDIF
   50    CONTINUE
C
C Inbetween the reactants and products, we need to add the
C the reaction directional symbol
C
         IF (J .EQ. 1) THEN
           LL = CKLSCH(ISTR)
           IF (ISKWRK(ISKWRK(IiNRPP)+IR-1) .LT. 0) THEN
             IF (LL+2 .GT. ILEN) GOTO 400
             ISTR(LL+1:LL+2) = '=>'
           ELSE
             IF (LL+3 .GT. ILEN) GOTO 400
             ISTR(LL+1:LL+3) = '<=>'
           ENDIF
         ENDIF
  100 CONTINUE
C
C Normal return from the program -> copy into return string.
C
      KERR = .FALSE.
      LT = CKLSCH(ISTR)
      IF (LT .GT. LEN(RNAM)) KERR = .TRUE.
      DO 160 N = 1, MIN(LT, LEN(RNAM))
        RNAM(N:N) = ISTR(N:N)
  160 CONTINUE
      RETURN
C
  400 CONTINUE
      WRITE (LOUT, 401) IR
  401 FORMAT (
     $ ' Error in SKSYMR... internal char string too small, IR = ', I3)
      LT = 0
      KERR = .TRUE.
      RETURN
C
  500 CONTINUE
      WRITE (LOUT, 501) IR
  501 FORMAT (' Error in SKSYMR... syntax error, IR = ', I3)
      LT = 0
      KERR = .TRUE.
      RETURN
C
C     end of SUBROUTINE SKSYMR
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSYMS (ISKWRK, CSKWRK, LOUT, KNAM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSYMS (ISKWRK, CSKWRK, LOUT, KNAM, KERR)
C  Returns a character array of species names.
C
C  INPUT
C  CSKWRK(*) - Character string workspace array; dimension at least
C              LENCSK.
C  LOUT      - Integer scalar, formatted output file unit number.
C
C  OUTPUT
C  KNAM(*)   - Character string array, species names;
C              dimension at least KKTOT, the total species count.
C  KERR      - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*)
      CHARACTER*(*) CSKWRK(*), KNAM(*)
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      KERR = .FALSE.
      ILEN = LEN(KNAM(1))
      I_KNAM = ISKWRK(IcKNAM) - 1
C
      DO 100 N = 1, ISKWRK(IiKTOT)
         I_KNAM = I_KNAM + 1
         KNAM(N) = ' '
         LT = CKLSCH(CSKWRK(I_KNAM))
         IF (LT .GT. ILEN) THEN
            KERR = .TRUE.
            WRITE (LOUT,500) CSKWRK(I_KNAM)
         ELSE
            KNAM(N) = CSKWRK(I_KNAM)
         ENDIF
  100 CONTINUE
C
  500 FORMAT (' Error in SKSYMS...character string length too small ')
C
C     end of SUBROUTINE SKSYMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKUML  (T, ISKWRK, RSKWRK, UML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKUML  (T, ISKWRK, RSKWRK, UML)
C  Returns an array of the internal energies in molar units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  UML(*)    - Real array, internal energies of the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), UML(*), T(*)
C
      CALL SKHML (T, ISKWRK, RSKWRK, UML)
      RU = RSKWRK(ISKWRK(IrRU))
      I_KTFL = ISKWRK(IiKTFL) - 1
      DO 100 K = 1, NKKGAS
         UML(K) = UML(K) - RU * T(ISKWRK(I_KTFL + K))
  100 CONTINUE
      DO 200 K = NKKGAS+1, ISKWRK(IiKTOT)
         UML(K) = UML(K) - RU * T(1)
  200 CONTINUE
C
C     end of SUBROUTINE SKUML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKUMS  (T, ISKWRK, RSKWRK, UMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKUMS  (T, ISKWRK, RSKWRK, UMS)
C  Returns an array of the internal energies in mass units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ISKWRK(*) - Integer workspace array; dimension at least LENISK.
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  UMS(*)    - Real array, internal energies of the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), UMS(*), T(*)
C
      CALL SKHMS (T, ISKWRK, RSKWRK, UMS)
      RU = RSKWRK(ISKWRK(IrRU))
      I_KWT = ISKWRK(IrKWT)
      I_KTFL = ISKWRK(IiKTFL)
C
      NKM1 = ISKWRK(IiKTOT) - 1
      DO 100 K = 0, NKM1
         IF (K .LT. NKKGAS) THEN
            TK = T(ISKWRK(I_KTFL + K))
         ELSE
            TK = T(1)
         ENDIF
         UMS(K+1) = UMS(K+1) - RU * TK / RSKWRK(I_KWT + K)
100   CONTINUE
C
C     end of SUBROUTINE SKUMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKWT   (ISKWRK, RSKWRK, WT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKWT   (ISKWRK, RSKWRK, WT)
C  Returns the molecular weights of the species.
C
C  INPUT
C  RSKWRK(*) - Real    workspace array; dimension at least LENRSK.
C
C  OUTPUT
C  WT(*)     - Real array, molecular masses for the species;
C              dimension at least KKTOT, the total species count.
C                 cgs units, gm/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.h'
      DIMENSION ISKWRK(*), RSKWRK(*), WT(*)
C
      I_KWT = ISKWRK(IrKWT) - 1
      DO 100 N = 1, ISKWRK(IiKTOT)
         WT(N) = RSKWRK(I_KWT + N)
  100 CONTINUE
C
C     end of SUBROUTINE SKWT
      RETURN
      END
