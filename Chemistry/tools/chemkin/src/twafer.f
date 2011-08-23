C     CVS $Revision: 1.1.1.1 $ reposited $Date: 2006/05/26 19:09:34 $

C///////////////////////////////////////////////////////////////////////
C
C     TWAFER-II
C
C     VERSION 1.48 OF JANUARY 1998
C
C     STEADY-STATE OR TRANSIENT SOLUTION TO CONDUCTIVE AND RADIATIVE
C     HEAT TRANSPORT IN A MULTIPLE-WAFER INTUBE LPCVD FURNANCE.
C
C     WRITTEN BY:
C
C        DR. WILLIAM G. HOUF            DR. JOSEPH F. GRCAR
C        DIVISION 8345                  DIVISION 8345
C        SANDIA NATIONAL LABORATORIES   SANDIA NATIONAL LABORATORIES
C        LIVERMORE, CA 94551-0969 USA   LIVERMORE, CA 94551-0969 USA
C        (510) 294-3184                 (510) 294-2662
C
C///////////////////////////////////////////////////////////////////////
C
C        will@california.sandia.gov     sepp@california.sandia.gov
C                                       na.grcar@na-net.ornl.gov
C
C///////////////////////////////////////////////////////////////////////
C
C     DOCUMENTATION:
C
C     W. G. Houf and J. F. Grcar, "TWAFER:  A Model for Determining
C     Temperatures in Batch LPCVD Reactors," Sandia National
C     Laboratories Report SANDXX-XXXX, Livermore, California.
C
C///////////////////////////////////////////////////////////////////////
C
C     USER LICENSE AGREEMENT:
C
C     TWAFER SOFTWARE
C     COPYRIGHT 1994 SANDIA CORPORATION
C     ALL RIGHTS RESERVED
C     EXCEPT AS SET FORTH IN THIS END USER LICENSE AGREEMENT
C
C     THIS Software was produced by Sandia Corporation under its
C     Contract No. DE_AC04-94AL85000 with the United States Department
C     of Energy for the management and operation of Sandia National
C     Laboratories, Livermore, California 94551-0969.  If you retain
C     possession of this SOFTWARE, you are agreeing to the following
C     terms and conditions of the license.  If you do not agree to
C     these terms and conditions of license you must destroy all forms
C     of this Software in your possession.
C
C     1.    Grants
C
C     1.1   Subject to the terms and conditions of this Agreement,
C           you are granted a royalty-free, nontransferable,
C           nonexclusive right and license, without the right to
C           sublicense or to sell or offer for sale in any form to:
C           use, modify, make copies, make derivative works or
C           compilations of the Software.
C
C     1.2   You agree to use this Software in the United States in
C           accordance with 35 USC 204.
C
C     1.3   If you desire a license right to sell: offer for sale;
C           display; or perform any portion of this Software or any
C           derivations or compilations thereof, you can obtain a
C           royalty bearing commercial license to do so by contacting
C           Gib Marguth at 510-294-2114.
C
C     1.4   You agree to provide Sandia with any enhancements,
C           including complete translations, made to the Software.
C           Such enhancements shall be useable royalty-free by Sandia
C           at Sandia National Laboratories, and may be transferred by
C           Sandia to other users of the Software.
C
C     1.5   You agree to use the Software at a single computer Site.
C           A separate License is required for each Site on which the
C           Software will be used.
C
C     1.6   You agree that a suitable protection plan will be in place
C           before the Software is installed on a computer system at
C           the agreed upon Site.  This plan will include the name of
C           the responsible individual at the facility who is charged
C           with enforcing its provisions and acting as the point of
C           contact with Sandia.
C
C     2.    Term of Agreement
C
C     2.1   The term of this End User Software License Agreement shall
C           continue for as long as you comply with all of the terms
C           and conditions of this License.
C
C     3.    Early Termination
C
C     3.1   If you fail to comply with any provision of this license,
C           your license will automatically terminate with your failure
C           to comply.
C
C     3.2   Upon termination, you shall, within thirty (30) days,
C           destroy all physical forms of the Software including any
C           and all derivatives or compilations thereof.
C
C     3.3   Sandia reserves all rights or remedies which Sandia may
C           have at law or in equity against those who do not accept
C           the rights or licenses provided herein and those who
C           violate the terms or conditions of this license.
C
C     4.    Disclaimer
C
C     4.1   NOTICE:  The Government is granted for itself and others
C           acting in its behalf a paid-up, nonexclusive, irrevocable
C           worldwide license in this data to reproduce, prepare
C           derivative works, and perform publicly and display publicly.
C           Beginning five (5) years after (date permission to assert
C           copyright was obtained) the Government is granted for itself
C           and others acting in its behalf a paid-up, nonexclusive,
C           irrevocable worldwide license in this data to reproduce,
C           prepare derivative works, distribute copies to the public
C           and display publicly, and to permit others to do so.
C           NEITHER SANDIA, THE UNITED STATES NOR THE UNITED STATES
C           DEPARTMENT OF ENERGY, NOR ANY OF THEIR EMPLOYEES MAKES ANY
C           WARRANTY, EXPRESS OR IMPLIED, INCLUDING ANY WARRANTY OF
C           MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR
C           ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE
C           ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION,
C           APPARATUS, PRODUCT, OR PROCESS OR REPRESENTS THAT ITS USE
C           WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, OR ASSUMES ANY
C           LIABILITY FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES RESULTING
C           FROM ITS USE BY ANYONE.
C
C     4.2   The U.S. Government is neither a party to nor assumes any
C           liability for activities of Sandia Corporation in
C           connection with this license.
C
C     5.    Indemnity
C
C     5.1   With regard to this Software, you indemnify Sandia
C           Corporation, the U.S. Government and all persons acting on
C           their behalf from and against all damages, costs, and
C           expenses, including attorney's fees, arising from personal
C           injury or property damage occurring as a result of the
C           making, using or selling of a product, process or service
C           by or on behalf of you, your assignees or licensees which
C           is derived from work performed under this agreement.
C
C     6.    Export Control Notice
C
C     6.1   THE LICENSED MATERIALS MAY BE EXPORT CONTROLLED AND
C           DISCLOSURE TO ANY NON U.S. CITIZEN WITHOUT OBTAINING A
C           LICENSE FROM THE U.S. GOVERNMENT MAY RESULT IN CRIMINAL
C           PROSECUTION OF THE DISCLOSER UNDER THE U.S. EXPORT CONTROL
C           LAWS OR REGULATIONS.
C
C     7.    Software Liaisons
C
C     7.1   All communication between the Licensee and Sandia shall be
C           with the following:
C
C           Dr. William G. Houf            Dr. Joseph F. Grcar
C           Sandia National Laboratories   Sandia National Laboratories
C           Mail Stop 9042                 Mail Stop 9042
C           Livermore, CA 94551-0969       Livermore, CA 94551-0969
C
C///////////////////////////////////////////////////////////////////////
C
C     SUBROUTINES:
C
C     TWAFER       MAIN SUBROUTINE
C       TWAF0        GRID, SHAPE FACTORS, AND INITIAL GUESS
C       . CONEND       FACTORS FOR THE REACTOR ENDS
C       . . CONSES       FACTOR BETWEEN A TUDE WALL AND ITSELF
C       . . CONSW        FACTOR BETWEEN A TUBE WALL AND TUBE END
C       . . CONWW        FACTOR BETWEEN FACING RINGS
C       . CONWAF       FACTORS BETWEEN TWO WAFERS
C       . . CONSS        FACTOR BETWEEN A TUBE WALL AND ITSELF
C       . . SCONSW       FACTOR BETWEEN A TUBE WALL AND A WAFER
C       . . SCONWW       FACTOR BETWEEN TWO WAFER REGIONS
C       . INTERP       LINEAR INTERPOLATION
C       TWAF1        EVALUATE THE RESIDUAL
C       . POW          DISTRIBUTE POWER FROM HEATERS TO WALL CELLS
C       TWAF2        PRINT THE SOLUTION
C       TWAF3        WRITE A PLOT2D FILE
C       TWAF4        WRITE A TEMPERATURE FILE FOR OVEND
C       TWAF5        WRITE THE MATRIX DATA FILE FOR SUNIL
C       . TWAF51       WRITE THE MATRIX ENTRIES.
C       TWAF7        PRINT THE SOLUTION
C
C///////////////////////////////////////////////////////////////////////
C
C     CHANGES WITH VERSION 1.12:
C
C     (1) RENAME SUBROUTINE TWAF11 TO TWAF0.
C     (2) REMOVE UNUSED ARGUMENT TEND (AND SIMILARLY NAMED ARGUMENTS)
C         FROM SUBROUTINES TWAF1 AND TWAF2.
C     (3) SHORTEN ARRAY TWALL TO TWAL IN SUBROUTINE TWAF1, AND SO ON,
C         TO ALLOW CONSISTENT NAMING: TWAL, TWAL0, TWAL4, ETC.
C     (4) ADD KEYWORD "WRITE THE COEFFICIENT MATRIX" AND SUBROUTINES
C         TWAF5 AND TWAF51 TO DO SAME.
C     (5) ADD NEW UNIT NUMBER TO MAIN PROGRAM AND RENAME UNIT NUMBERS
C         SO DATA1, DATA2, AND DATA3 ARE UNITS FOR THE THREE DATA FILES.
C     (6) CORRECT TWAF3.
C
C     CHANGES WITH VERSION 1.13:
C
C     (1) HEATER POWERS NOW SPECIFIED ON INPUT.
C
C     CHANGES WITH VERSION 1.14:
C
C     (1) CORRECT TWAF5.
C     (2) BREAK OVERLONG DATA STATEMENT IN TWAFER.
C
C     CHANGES WITH VERSION 1.15:
C
C     (1) CORRECT THE LABELING OF EXIT REGION RESIDUALS AND VARIABLES
C         IN TWAF5.
C
C     CHANGES WITH VERSION 1.16:
C
C     (1) UPDATED TO MODEL GAP BETWEEN END WAFERS AND REACTOR WALL AS
C         A RERADIATING SURFACE.
C     (2) ADDED SUBOUTINE TO FIND POWER DISSIPATED IN EACH WALL CELL.
C     (3) POWER AND PWALL ARE NOW IN WATTS AND NOT POWER PER UNIT
C         VOLUME.
C     (4) CORRECTED ERROR IN UNITS FOR THE HEAT TRANSFER COEFFICIENT H.
C         CHANGED FROM 4.55 TO 4.55E-04.
C     (5) CHANGED SUM IN DO LOOP 3020 TO INCLUDE WALL ELEMENT
C         SELF-CONFIGURATION FACTOR IN SUMMATION.  CHANGED CODE SO THAT
C         SELF-CONFIGURATION FACTOR IS NOT INCLUDE IN THE WALL RADIOSITY
C         RESIDUAL BUT IS INCLUDED IN THE WALL IRRADIATION (GWALL).
C
C     CHANGES WITH VERSION 1.17:
C
C     (1) CORRECTED INDEX OF CP AND RHO USED FOR INSULATION ON INLET
C         AND EXIT.
C     (2) CHANGE GRID TO REMOVE DUPLICATE POINTS WHERE ENTRANCE AND
C         EXIT SECTIONS CONNECT TO WAFER SECTION.
C
C     CHANGES WITH VERSION 1.18:
C
C     (1) ADDED TWO SEMITRANSPARENT QUARTZ TUBES BETWEEN HEATERS AND
C         WAFERS.
C     (2) UPDATED MATRIX OUTPUT FOR ADDED NEW VARIABLES.
C
C     CHANGES WITH VERSION 1.19:
C
C     (1) SHAPE FACTORS BETWEEN SEMITRANSPARENT QUARTZ TUBES TREATED
C         AS INFINITE CONCENTRIC CYLINDERS.
C     (2) SHAPE FACTORS AND RADIATION EXCHANGE BETWEEN QUARTZ TUBES
C         TREATED BY ASSUMING INFINITELY THIN CYLINDERS WITH RADII
C         EQUAL TO THE INNER RADII OF THE QUARTZ TUBES.
C
C     CHANGES WITH VERSION 1.20:
C     (1) ADDED MODEL FOR THERMOCOUPLES BETWEEN QUARTZ WALLS. (TWAF6.F)
C
C     CHANGES WITH VERSION 1.21:
C     (1) SEVERAL BUG FIXES.
C     (2) ADDITION OF NEW KEYWORDS FOR VARIOUS MATERIAL PROPERTIES IN
C         THE REACTOR MODEL.
C     (3) CORRECTIONS IN THE MATRIX WRITING ROUTINES, ISSUED AS 120B,
C         ARE INCORPORATED HERE.
C     (4) SUBROUTINE TWAF6 FOR THERMOCOUPLE TEMPERATURES INCORPORATED
C         INTO TWAF1.  TWAF6 DISCARDED.
C     (5) MATRIX ROUTINES ALTERED TO WRITE THERMOCOUPLE MATRIX TOO.
C
C     CHANGES WITH VERSION 1.23:
C     (1) MERGE VERSION 1.21 (WITH CORRECT MATRIX PRINTING) WITH VERSION
C         1.22 (WITH STAGGERED WAFER AND TUBE GRIDS).
C
C     CHANGES WITH VERSION 1.24:
C     (1) CORRECTED RADIATION FROM WAFER EDGE FOR STAGGERED GRID.
C     (2) DELETE DUMMY EXTERNAL POINTS FROM THE MATRIX.
C     (3) REMOVE CHANGE BLOCKS FOR OPENING THE MATRIX FILE.
C
C     CHANGES WITH VERSION 1.25:
C     (1) REVERSE CHANGE (2) OF VERSION 1.24.
C     (2) CHANGE TWAF1 TO PRODUCE FOURTH POWER OF THERMOCOUPLE VALUES.
C
C     CHANGES WITH VERSION 1.26:
C     (1) ADDED DIFFERENT HEAT TRANSFER COEFFICIENTS AND EMISSIVITIES
C         FOR THE ENTRANCE, MAIN, AND EXIT SECTIONS OF THE REACTOR
C         EXTERIOR SKIN.
C
C     CHANGES WITH VERSION 1.27:
C     (1) ADDED THERMOCOUPLE MODEL COMMUNICATION WITH ENTRANCE AND
C         EXIT SECTIONS.
C
C     CHANGES WITH VERSION 1.28:
C     (1) CORRECTED TRANSIENT TERMS FOR HEATER AND INSULATION OVER
C         WAFER REGION.
C
C     CHANGES WITH VERSION 1.29:
C     (1) ALLOW HEAT LOSS FROM OUTER CYLINDER OF END INSULATION.
C     (2) FIRST GRID POINT ON ENTRANCE AND EXIT MOVED TO CENTER OF CELL.
C     (3) ADDED HEAT LOSS FROM ENDS OF QUARTZ TUBES, HEATERS, AND
C         HEATER INSULATION BLANKET.
C
C     CHANGES WITH VERSION 1.30:
C     (1) ENCAPSULATED END REGIONS IN SUBROUTINE EBAL1.
C
C     CHANGE WITH VERSION 1.31:
C     (1) ADDED POWER DISSIPATION TO ENTRANCE AND EXIT CAPS.
C     (2) RENAMED ARRAY RWORK IN SUBROUTINE TWAFER TO R SO AS TO
C         SHORTEN LONG ARGUMENT LISTS TO 19 CONTINUATION LINES.
C     (3) MAKE SUBROUTINE TWAFER INTO A FORTRAN MAIN PROGRAM, THUS
C         DISCARDING TWAFER_MAIN.
C
C     CHANGES WITH VERSION 1.32:
C     (1) CHANGED SUBROUTINE TWAF4 TO OUTPUT A TEMPERATURE FILE
C         COMPATIBLE WITH OVEND_241.
C
C     CHANGES FROM VERSION 1.33:
C     BEGINNING TO MAKE CHANGES TO CONVERT TO TIME-ACCURATE CODE.
C     (1) REMOVED ENTRANCE AND EXIT RADIOSITIES FROM TWOPNT UNKNOWNS.
C     (2) ALL VARIABLES FOR SUB. TWAF1 (EXCEPT LOGICALS) PACKED IN
C         R(*), AND IWORK(*).
C     (3) VOL*RHO*CP TERM REMOVED FROM TRANSIENT TERM OF RESIDUAL AND
C         MOVED TO OTHER SIDE OF EQUATION AS DIVISOR.
C
C     CHANGES FROM VERSION 1.34:
C     (1) DASSL TIME INTEGRATION ADDED FOR TRANSIENT SOLUTION.
C     (2) KEYWORDS ADDED SO STEADY-STATE SOLUTION CAN BE FOUND USING
C         TWOPNT OR TRANSIENT SOLUTION WITH DASSL.
C     (3) TWAF2 AND TWAF3 CHANGED TO ALLOW OUTPUT AND PLOTTING FILES
C         TO BE WRITTEN FOR EITHER STEADY-STATE OR TRANSIENT CASE.
C     (4) T.DAT FILE ADDED FOR OUTPUTTING WAFER CENTERLINE, WAFER EDGE,
C         QUARTZ WALLS, HEATER, AND INSULATION TEMPERATURES FOR THE
C         WAFER AT THE CENTER OF THE LOAD.
C
C     CHANGES FROM VERSION 1.35:
C     (1) MEMORY UTILIZATION IMPROVED.
C     (2) WORK ARRAYS FOR DASSL COMBINED INTO R AND IWORK.
C
C     CHANGES WITH VERSION 1.36:
C     (1) CORRECTED BUG IN IPOINT COMMON BLOCK.
C
C     CHANGE WITH VERSION 1.37:
C     (1) RING BOAT MODEL ADDED (WAFER GRID CHANGED SO THAT NEAREST
C         CELL FACE IS MOVED TO LIE ON RING EDGE).
C     (2) ADDITIONAL VARIABLES OPTIONS ADDED TO RESIDUAL (NOT AVAILABLE
C         THROUGH INPUT FILE).
C
C     CHANGE WITH VERSION 1.38:
C     (1) SAME AS TWAFER1.37 BUT LICENSE AGREEMENT ADDED
C         THE BELOW CHANGES.
C     (2) SOME BUG FIXES WITH READING INPUT FILE.
C     (3) CHANGED NAME TO TWAFER II.
C
C     CHANGE WITH VERSION 1.39:
C     (1) FIRST VERSION WITH AIR COOLING BETWEEN QUARTZ JARS
C         AND AIR COOLING BETWEEN OUTER QUARTZ JAR AND HEATER.
C         PRESENTLY ONLY AIR FLOW FROM BOTTOM UP IS ALLOWED.
C         USER INTERFACE FOR COOLING NOT INCORPORATED IN KEYWORD
C         INPUT AS YET.
C
C     CHANGES WITH VERSION 1.41:
C     (1) AIR COOLING FROM BOTTOM UP OR TOP DOWN ALLOWED THROUGH
C         EITHER CHANNEL.  USER INTERFACE FOR COOLING NOT
C         INCORPORATED IN KEYWORD INPUT AS YET.
C     (2) SOME BUG FIXES ON AIR-COOLING MODEL.
C     (3) T.DAT - FILE MODIFIED FOR MORE OUTPUT.
C
C     CHANGES WITH VERSION 1.42:
C     (1) APPROXIMATE PEDESTAL MODEL - CONTROLLED BY LOGICAL VARIABLE
C         PED IN SUBROUTINE TWAFER.F.  IF INVOKED FIRST WAFER IS
C         TEATED AS QUARTZ DISC AND REST OF PEDESTAL PLACED IN INNER
C         QUARTZ TUBE OF ENTRANCE REGION.
C
C     CHANGES WITH VERSION 1.43:
C     (1) ADDED MODEL FOR RADIAL INJECTION BETWEEN OUTER QUARTZ TUBES
C         AND HEATERS BY ADDING COLD GAS TO AXIAL FLOW.  NEW SUBROUTINE
C         NUSSL.F IS USED THAT SMOOTHES LAMINAR AND TURBULENT HEAT
C         TRANSFER CORRELATIONS TOGETHER INTO ONE CURVE.
C
C     CHANGES WITH VERSION 1.44:
C     (1) REYNOLDS NUMBER FOR HEAT TRANSFER CORRELATIONS FOR COOLING
C         FLOW NOW COMPUTED BASED ON LOCAL REYNOLDS NUMBER AND
C         INDIVIDUAL CELL MASS FLOWRATES.
C
C     CHANGES WITH VERSION 1.45:
C     (1) SEVERAL BUG FIXES FOUND BY SEMITHERM, DEALING WITH CALL LISTS
C         OF SUBROUTINE TWAF1.
C
C     CHANGES WITH VERSION 1.46:
C     (1) INSTALL NEW TWOPNT, VERSION 3.27.
C
C     CHANGES WITH VERSION 1.47:
C     (1) SOME BUG FIXES FOR DO LOOP INDICIES BEYOND ARRAY BOUNDS.
C     (2) INITIALIZED RADIOSITIES IN ENTRANCE AND EXIT FOR INITIAL
C         PRINTOUT.
C
C     CHANGES WITH VERSION 1.48:
C     (1) ADDED ERF AND DERF OBTAINED FROM SLATEC VIA NETLIB.
C
C///////////////////////////////////////////////////////////////////////

      SUBROUTINE CONCC
     +  (ERROR, TEXT,
     +   F12, F21, R1, R2, Z1H, Z1L, Z2H, Z2L)

C///////////////////////////////////////////////////////////////////////
C
C     CONCC
C
C///  COMPUTES THE CONFIGUATION FACTOR BETWEEN RINGS ON CONCENTRIC
C     CYLINDERS.
C
C     AREA1 = AREA OF ANNULAR RING ON OUTSIDE OF THE INNER CYLINDER.
C     AREA2 = AREA OF ANNULAR RING ON THE INSIDE OF THE OUTER CYLINDER.
C
C     FORMULA FROM REA REFERENCE.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     INPUT:
C
C     R1      INPUT REAL - RADIUS OF THE INNER CYLINDER.
C
C     R2      INPUT REAL - RADIUS OF THE OUTER CYLINDER.
C
C     Z1H     INPUT REAL - Z COORDINATE OF UPPER EDGE OF RING ON
C                          INNER CYLINDER
C
C     Z1L     INPUT REAL - Z COORDINATE OF LOWER EDGE OF RING ON
C                          INNER CYLINDER
C
C     Z2H     INPUT REAL - Z COORDINATE OF UPPER EDGE OF RING ON
C                          OUTER CYLINDER
C
C     Z2L     INPUT REAL - Z COORDINATE OF LOWER EDGE OF RING ON
C                          OUTER CYLINDER
C
C     OUTPUT:
C
C     F12     OUTPUT REAL - SHAPE FACTOR FROM THE RING ON THE INNER
C                           CYLINDER TO THE RING ON THE OUTER CYLINDER.
C
C     F21     OUTPUT REAL - SHAPE FACTOR FROM THE RING ON THE OUTER
C                           CYLINDER TO THE RING ON THE INNER CYLINDER.
C
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER TEXT
      LOGICAL ERROR, CASE1, CAS2R, CAS2L, CASE3, CAS4R, CAS4L
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   A1, A2, DUM, F12, F21, FA, FB, FC, FD, PI, R1, R2, Z1H, Z1L,
     +   Z2H, Z2L

      PARAMETER (ID = ' CONCC:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE
C
C///////////////////////////////////////////////////////////////////////

C///  CHECKS THE ARGUMENTS.

      ERROR = (R1 .GE. R2)
      IF (ERROR) GO TO 9001

      ERROR = ((Z1H .LT. 0.0).OR.(Z1L .LT. 0.0))
      ERROR = ((Z2H .LT. 0.0).OR.(Z2L .LT. 0.0))
      IF (ERROR) GO TO 9002

      ERROR = ((R1 .LE. 0.0) .OR. (R2 .LE. 0.0))
      IF (ERROR) GO TO 9003

      ERROR = (Z1H .LT. Z1L)
      IF (ERROR) GO TO 9004

      ERROR = (Z2H .LT. Z2L)
      IF (ERROR) GO TO 9005

      CASE1 = .FALSE.
      CAS2R = .FALSE.
      CAS2L = .FALSE.
      CASE3 = .FALSE.
      CAS4L = .FALSE.
      CAS4R = .FALSE.

C///////////////////////////////////////////////////////////////////////
C
C     (2) FIND THE CONFIGURATION FACTORS.
C
C///////////////////////////////////////////////////////////////////////

C///  DETERMINE THE ORIENTATION OF THE RINGS.

      CASE1 = ((Z1L .LE. Z2L) .AND. (Z1H .GE. Z2H))
      CAS2L = ((Z1L .LE. Z2L) .AND. ((Z2L .LE. Z1H).AND.(Z1H .LE. Z2H)))
      CAS2R = ((Z1H .GE. Z2H) .AND. ((Z2L .LE. Z1L).AND.(Z1L .LE. Z2H)))
      CASE3 = (((Z2L .LE. Z1L).AND.(Z1L .LE. Z2H)) .AND.
     +         ((Z2L .LE. Z1H).AND.(Z1H .LE. Z2H)))
      CAS4L = ((Z1L .LE. Z2L) .AND. (Z1H .LE. Z2L))
      CAS4R = ((Z1L .GE. Z2H) .AND. (Z1H .GE. Z2H))

C///  FIND COMPONET CONFIGURATION FACTORS FOR CYLINDER TO DISC

      IF (CASE1) THEN

C        SUBROUTINE CONCD
C       +  (ERROR, TEXT, F12, F21, L1, R1, R2)

         CALL CONCD(ERROR, TEXT, FA, DUM, (Z2L - Z1L), R1, R2)
         CALL CONCD(ERROR, TEXT, FB, DUM, (Z1H - Z2H), R1, R2)
         CALL CONCD(ERROR, TEXT, FC, DUM, (Z2H - Z1L), R1, R2)
         CALL CONCD(ERROR, TEXT, FD, DUM, (Z1H - Z2L), R1, R2)

         F12 = ((Z2H - Z2L)+ FA*(Z2L - Z1L) + FB*(Z1H - Z2H) -
     +          FC*(Z2H - Z1L) - FD*(Z1H - Z2L))/(Z1H - Z1L)

      ELSEIF (CAS2L) THEN

C        SUBROUTINE CONCD
C       +  (ERROR, TEXT, F12, F21, L1, R1, R2)

         CALL CONCD(ERROR, TEXT, FA, DUM, (Z2L - Z1L), R1, R2)
         CALL CONCD(ERROR, TEXT, FB, DUM, (Z2H - Z1H), R1, R2)
         CALL CONCD(ERROR, TEXT, FC, DUM, (Z1H - Z2L), R1, R2)
         CALL CONCD(ERROR, TEXT, FD, DUM, (Z2H - Z1L), R1, R2)

         F12 = ((Z1H - Z2L) + FA*(Z2L - Z1L) + FB*(Z2H - Z1H) -
     +          FC*(Z1H - Z2L) - FD*(Z2H - Z1L))/(Z1H - Z1L)

      ELSEIF (CAS2R) THEN

C        SUBROUTINE CONCD
C       +  (ERROR, TEXT, F12, F21, L1, R1, R2)

         CALL CONCD(ERROR, TEXT, FA, DUM, (Z1H - Z2H), R1, R2)
         CALL CONCD(ERROR, TEXT, FB, DUM, (Z1L - Z2L), R1, R2)
         CALL CONCD(ERROR, TEXT, FC, DUM, (Z2H - Z1L), R1, R2)
         CALL CONCD(ERROR, TEXT, FD, DUM, (Z1H - Z2L), R1, R2)

         F12 = ((Z2H - Z1L) + FA*(Z1H - Z2H) + FB*(Z1L - Z2L) -
     +          FC*(Z2H - Z1L) - FD*(Z1H - Z2L))/(Z1H - Z1L)

      ELSEIF (CASE3) THEN

C        SUBROUTINE CONCD
C       +  (ERROR, TEXT, F12, F21, L1, R1, R2)

         CALL CONCD(ERROR, TEXT, FA, DUM, (Z1L - Z2L), R1, R2)
         CALL CONCD(ERROR, TEXT, FB, DUM, (Z2H - Z1H), R1, R2)
         CALL CONCD(ERROR, TEXT, FC, DUM, (Z1H - Z2L), R1, R2)
         CALL CONCD(ERROR, TEXT, FD, DUM, (Z2H - Z1L), R1, R2)

         F12 = 1.0 + (FA*(Z1L - Z2L) + FB*(Z2H - Z1H) -
     +         FC*(Z1H - Z2L) - FD*(Z2H - Z1L))/(Z1H - Z1L)

      ELSEIF (CAS4L) THEN

C        SUBROUTINE CONCD
C       +  (ERROR, TEXT, F12, F21, L1, R1, R2)

         CALL CONCD(ERROR, TEXT, FA, DUM, (Z2L - Z1L), R1, R2)
         CALL CONCD(ERROR, TEXT, FB, DUM, (Z2H - Z1H), R1, R2)
         CALL CONCD(ERROR, TEXT, FC, DUM, (Z2L - Z1H), R1, R2)
         CALL CONCD(ERROR, TEXT, FD, DUM, (Z2H - Z1L), R1, R2)

         F12 = (FA*(Z2L - Z1L) + FB*(Z2H - Z1H) -
     +          FC*(Z2L - Z1H) - FD*(Z2H - Z1L))/(Z1H - Z1L)

      ELSEIF (CAS4R) THEN

C        SUBROUTINE CONCD
C       +  (ERROR, TEXT, F12, F21, L1, R1, R2)

         CALL CONCD(ERROR, TEXT, FA, DUM, (Z1H - Z2H), R1, R2)
         CALL CONCD(ERROR, TEXT, FB, DUM, (Z1L - Z2L), R1, R2)
         CALL CONCD(ERROR, TEXT, FC, DUM, (Z1L - Z2H), R1, R2)
         CALL CONCD(ERROR, TEXT, FD, DUM, (Z1H - Z2L), R1, R2)

         F12 = (FA*(Z1H - Z2H) + FB*(Z1L - Z2L) -
     +          FC*(Z1L - Z2H) - FD*(Z1H - Z2L))/(Z1H - Z1L)

      ELSE
      ENDIF

C///  AREAS

      A1 = 2.*R1*PI*(Z1H - Z1L)
      A2 = 2.*R2*PI*(Z2H - Z2L)

C///  FINAL CONFIGURATION FACTORS

      F21 = A1/A2*F12

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, R1, R2
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, Z1H, Z1L, Z2H, Z2L
      GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, R1, R2
      GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, Z1H, Z1L
      GO TO 99999

9005  IF (0 .LT. TEXT) WRITE (TEXT, 99005) ID, Z2H, Z2L
      GO TO 99999

99001 FORMAT
     +   (/1X, A, 'ERROR. INNER CYLINDER RADIUS IS GREATER THAN',
     +    /10X, 'OUTER CYLINDER RADIUS.',
     +   //10X, 'R1 = ',E10.3,
     +    /10X, 'R2 = ',E10.3)

99002 FORMAT
     +   (/1X, A, 'ERROR. Z COORDINATES OF RINGS MUST BE POSITIVE',
     +   //10X, 'Z1H = ',E10.3,
     +   //10X, 'Z1L = ',E10.3,
     +   //10X, 'Z2H = ',E10.3,
     +   //10X, 'Z2L = ',E10.3)

99003 FORMAT
     +   (/1X, A, 'ERROR. ONE OR MORE OF THE CYLINDER RADII'
     +    /10X, 'ARE LESS THAN OR EQUAL TO ZERO.',
     +   //10X, 'R1 = ',E10.3,
     +    /10X, 'R2 = ',E10.3)

99004 FORMAT
     +  (/1X, A, 'ERROR. THE COORDINATE OF THE LOWER EDGE OF THE INNER'
     +    /10X, 'CYLINDER RING IS LESS THAN THE UPPER COORDINATE.',
     +   //10X, 'Z1H = ',E10.3,
     +    /10X, 'Z1L = ',E10.3)

99005 FORMAT
     +  (/1X, A, 'ERROR. THE COORDINATE OF THE LOWER EDGE OF THE OUTER'
     +    /10X, 'CYLINDER RING IS LESS THAN THE UPPER COORDINATE.',
     +   //10X, 'Z2H = ',E10.3,
     +    /10X, 'Z2L = ',E10.3)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE CONCD
     +  (ERROR, TEXT,
     +   F12, F21, L1, R1, R2)

C///////////////////////////////////////////////////////////////////////
C
C     CONCD
C
C///  COMPUTES THE CONFIGUATION FACTOR BETWEEN THE OUTER SURFACE OF A
C     CYLINDER TO THE ANNULAR DISK AT THE END OF THE CYLINDER.
C
C     AREA1 = AREA OF THE CYLINDER.
C     AREA2 = ANNULAR RING ATTACHED TO BASE OF THE CYLINDER.
C
C     FORMULA C-72 FROM "A CATALOG OF RADIATION CONFIGURATION FACTORS",
C     J.R. HOWELL, (EXPRESSION CORRECTED ACCORDING TO ORIGONAL
C     REFERENCE FROM REA.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     L1      INPUT REAL - LENGTH OF THE CYLINDER.
C
C     R1      INPUT REAL - RADIUS OF THE CYLINDER.
C
C     R2      INPUT REAL - RADIUS OF THE DISK.
C
C     F12     OUTPUT REAL - SHAPE FACTOR FROM THE CYLINDER TO
C                           ANNULAR DISC AT THE BASE OF CYLINDER.
C
C     F21     OUTPUT REAL - SHAPE FACTOR FROM THE ANNULAR DISC AT THE
C                           BASE OF CYLINDER TO THE CYLINDER.
C
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   TEXT

      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   A, A1, A2, B, F12, F21, L, L1, PI, R, R1, R2

      PARAMETER (ID = 'CONCD:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE
C
C///////////////////////////////////////////////////////////////////////

C///  CHECKS THE ARGUMENTS.

      ERROR = (R1 .GE. R2)
      IF (ERROR) GO TO 9001

      ERROR = (L1 .LT. 0.0)
      IF (ERROR) GO TO 9002

      ERROR = ((R1 .LE. 0.0) .OR. (R2 .LE. 0.0))
      IF (ERROR) GO TO 9003

C///////////////////////////////////////////////////////////////////////
C
C     (2) FIND THE CONFIGURATION FACTORS.
C
C///////////////////////////////////////////////////////////////////////

      IF (L1 .EQ. 0.0) THEN
         F12 = 0.0
         F21 = 0.0
         RETURN
      ELSE
      ENDIF

      R = R1/R2
      L = L1/R2
      A = L*L + R*R - 1.0
      B = L*L - R*R + 1.0

C///  CONFIGURATION FACTOR FROM THE CYLINDER TO THE DISC
C     FROM HOWELL J.R. HOWELL (CORRECTED FROM REA REF.)

C*****DOUBLE PRECISION
      F12 =
     + B/(8.0*R*L) + (DACOS(A/B) - DSQRT((L**2 + R**2 + 1.0)**2/R**2 -
     + 4.0)/(2.0*L)*DACOS(A*R/B) - A/(2.0*R*L)*DASIN(R))/(2.0*PI)
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      F12 =
C     + B/(8.0*R*L) + (ACOS(A/B) - SQRT((L**2 + R**2 + 1.0)**2/R**2 -
C     + 4.0)/(2.0*L)*ACOS(A*R/B) - A/(2.0*R*L)*ASIN(R))/(2.0*PI)
C*****END SINGLE PRECISION

C///  AREAS

      A1 = 2.*R1*PI*(L1)
      A2 = PI*(R2**2 - R1**2)

C///  FINAL CONFIGURATION FACTORS

      F21 = A1/A2*F12

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, R1, R2

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, L

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, R1, R2

      GO TO 99999

99001 FORMAT
     +   (/1X, A, 'ERROR. CYLINDER RADIUS IS GREATER THAN',
     +    /10X, 'DISC RADIUS.',
     +   //10X, 'R1 = ',E10.3,
     +    /10X, 'R2 = ',E10.3)
99002 FORMAT
     +   (/1X, A, 'ERROR. LENGTH OF CYLINDER LESS THAN ZERO',
     +   //10X, 'L = ',E10.3)

99003 FORMAT
     +   (/1X, A, 'ERROR. RADIUS OF CYLINDER OR RADIUS OF DISC'
     +    /10X, 'IS LESS THAN OR EQUAL TO ZERO.',
     +   //10X, 'R1 = ',E10.3,
     +    /10X, 'R2 = ',E10.3)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE CONEND
     +  (ERROR, TEXT,
     +   FEE, FEG, FESE, FEW, FGE, FGSE, FSEE, FSEG, FSES, FSEW, FWE,
     +   FWSE, IOPT, PDIM, REFLEC, RTUBE, RPNTS, RWALL, ZPNTS, ZWALL)

C///////////////////////////////////////////////////////////////////////
C
C     CONEND
C
C///  COMPUTES THE CONFIGUATION FACTORS FOR THE END GEOMETRY OF THE
C     REACTOR ASSUMING THE WAFER IS A SPECULAR REFLECTOR AND ALL OTHER
C     SURFACES ARE DIFFUSE. ASSUMES ZWALL(1) IS AT THE FACE OF THE
C     REACTOR END WALL AND ZWALL(ZPNTS+1) IS AT THE FACE OF THE FIRST
C     WAFER.
C
C     CHANGES:
C     (1) CORRECTED INDEXING ERROR IN CALL TO SUBROUTINE CONWW.
C     (2) ADDED SHAPE FACTORS TO GAP BETWEEN FIRST WAFER EDGE AND
C         REACTOR WALL.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     IOPT    INPUT INTEGER - IOPT = 2 THEN RADIUS = (RADIUS OF THE
C             WAFERS) ELSE RADIUS = (RADIUS OF THE REACTOR)
C
C     REFLEC  INPUT REAL - REFLECTIVITY OF WAFERS.
C
C     RTUBE   INPUT REAL - RADIUS OF THE REACTOR TUBE.
C
C     RPNTS   INPUT INTEGER - NUMBER OF RADIAL GRID POINTS ON A WAFER.
C
C     RWALLS  INPUT REAL DIMENSIONED (RPNTS + 1) - RADIAL LOCATIONS OF
C             CELL WALLS ON WAFER.
C
C     SPACE   INPUT REAL - SPACING BETWEEN WAFERS.
C
C     ZPNTS   INPUT REAL - NUMBER OF POINTS ALONG THE WALL BETWEEN
C             REACTOR END WALL AND FIRST WAFER.
C
C     ZWALL   INPUT REAL DIMENSIONED (ZPNTS+1) - AXIAL LOCATIONS OF
C             CELL WALLS ALONG REACTOR
C                          ZWALL(1) = 0.0,
C                          ZWALL(ZPNTS+1) = LOCATION OF FIRST WAFER.
C
C     FEE     OUTPUT REAL - SHAPE FACTOR FROM REACTOR END TO ITSELF.
C
C     FEG     OUTPUT REAL - SHAPE FACTOR FROM REACTOR END TO THE GAP
C             BETWEEN THE FIRST WAFER AND THE WALL.
C
C     FESE    OUTPUT REAL DIMENSIONED (ZPNTS) - SHAPE FACTORS FROM
C             REACTOR END WALL TO ELEMENTS ON REACTOR SIDE WALL.
C
C     FEW     OUTPUT REAL DIMENSIONED (RPNTS) - SHAPE FACTORS FROM
C             RING ELEMENTS ON FIRST WAFER TO REACTOR END WALL.
C
C     FGE     OUTPUT REAL - SHAPE FACTOR FROM THE GAP BETWEEN THE
C             FIRST WAFER AND THE WALL TO THE REACTOR END WALL.
C
C     FGSE    OUTPUT REAL DIMENSIONED (ZPNTS) - SHAPE FACTORS FROM
C             GAP BETWEEN THE FIRST WAFER AND WALL TO RING ELEMENTS
C             ON REACTOR SIDE WALL.
C
C     FSEE    OUTPUT REAL DIMENSIONED (ZPNTS) - SHAPE FACTORS FROM
C             ELEMENTS ON REACTOR SIDE WALL TO REACTOR END WALL.
C
C     FSEG    OUTPUT REAL DIMENSIONED (ZPNTS) - SHAPE FACTORS FROM
C             RING ELEMENTS ON REACTOR SIDE WALL TO GAP BETWEEN THE
C             FIRST WAFER AND WALL.
C
C     FSES    OUTPUT REAL DIMENSIONED (PDIM, ZPNTS) - SHAPE FACTORS
C             BETWEEN ELEMENTS ON REACTOR SIDE WALL IN THE ENTRANCE
C             REGION.
C
C     FSEW    OUTPUT REAL DIMENSIONED (PDIM, RPNTS) - SHAPE FACTORS
C             FROM RING ELEMENTS ON REACTOR SIDE WALL ELEMENTS TO
C             RING ELEMENTS ON FIRST WAFER.
C
C     FWE     OUTPUT REAL DIMENSIONED (RPNTS) - SHAPE FACTORS FROM
C             REACTOR END WALL TO A RING ELEMENTS ON THE FIRST WAFER.
C
C     FWSE    OUTPUT REAL DIMENSIONED (RPNTS, ZPNTS) - SHAPE FACTORS
C             FROM RING ELEMENTS ON FIRST WAFER TO REACTOR SIDE WALL
C             ELEMENTS.
C
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   RPNTS, I, IOPT, J, PDIM, TEXT, ZPNTS

      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   DUMMY, FDIR1, FDIR2, FEE, FEG, FESE, FEW, FGSE, FGE, FREFL1,
     +   FREFL2, FSEE, FSEG, FSES, FSEW, FWE, FWSE, PI, RADIUS, REFLEC,
     +   RTUBE, RWALL, ZERO, ZLEN, ZLEN2, ZWALL

      PARAMETER (ID = 'CONEND:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)
      PARAMETER (ZERO = 0.0D0)

      DIMENSION
     +   FESE(ZPNTS), FEW(RPNTS), FGSE(ZPNTS), FSEE(ZPNTS), FSEG(ZPNTS),
     +   FSES(PDIM, ZPNTS), FSEW(PDIM, RPNTS), FWE(RPNTS),
     +   FWSE(RPNTS, ZPNTS), RWALL(RPNTS+1), ZWALL(ZPNTS+1)

C///////////////////////////////////////////////////////////////////////
C
C     (1)  CONFIGURATION FACTORS
C
C///////////////////////////////////////////////////////////////////////

      RADIUS = RTUBE
      IF (IOPT. EQ. 2) RADIUS = RWALL(RPNTS+1)

      ZLEN = ZWALL(ZPNTS+1) - ZWALL(1)
      ZLEN2 = 2.0*ZLEN

      DO 200 J = 1, RPNTS
         DO 100 I = 1, ZPNTS

C///        WAFER TO ENTRANCE SIDE WALL

C           SUBROUTINE CONSW
C       +     (ERROR, TEXT,
C       +      F12, F21, LLOWER, LUPPER, RADIUS, RIN, ROUT)

            CALL CONSW
     +        (ERROR, TEXT,
     +         FSEW(I, J), FWSE(J, I), (ZWALL(ZPNTS+1) - ZWALL(I+1)),
     +         (ZWALL(ZPNTS+1) - ZWALL(I)), RADIUS, RWALL(J),
     +         RWALL(J+1))
100      CONTINUE

C///     WAFER TO END AND END TO WAFER

C        SUBROUTINE CONWW
C     +    (ERROR, TEXT,
C     +     A, F12, F21, RIN1, RIN2, ROUT1, ROUT2)

         CALL CONWW
     +     (ERROR, TEXT,
     +      ZWALL(ZPNTS+1), FEW(J), FWE(J), ZERO, RWALL(J),
     +      RADIUS, RWALL(J+1))

200   CONTINUE

C///  CONFIGURATION FACTORS BETWEEN SIDEWALL ELEMENTS AND GAP BETWEEN
C     THE FIRST WAFER AND REACTOR WALL

      DO 250 I = 1, ZPNTS

C         SUBROUTINE CONSW
C       +   (ERROR, TEXT,
C       +    F12, F21, LLOWER, LUPPER, RADIUS, RIN, ROUT)

          CALL CONSW
     +      (ERROR, TEXT,
     +       FSEG(I), FGSE(I), (ZWALL(ZPNTS+1) - ZWALL(I+1)),
     +       (ZWALL(ZPNTS+1) - ZWALL(I)), RADIUS, RWALL(RPNTS+1),
     +       RADIUS)
250   CONTINUE

C///  CONFIGURATION FACTOR BETWEEN REACTOR END WALL AND GAP BETWEEN
C     THE FIRST WAFER AND REACTOR WALL

C        SUBROUTINE CONWW
C     +    (ERROR, TEXT,
C     +     A, F12, F21, RIN1, RIN2, ROUT1, ROUT2)

         CALL CONWW
     +     (ERROR, TEXT,
     +      ZWALL(ZPNTS+1), FEG, FGE, ZERO, RWALL(RPNTS+1),
     +      RADIUS, RADIUS)

C///  CONFIGURATION FACTOR BETWEEN REACTOR END WALL AND SIDE ELEMENTS
C     ON REACTOR ENTRANCE

      DO 300 I = 1, ZPNTS

C        SUBROUTINE CONSW
C       +  (ERROR, TEXT,
C       +   F12, F21, LLOWER, LUPPER, RADIUS, RIN, ROUT)

C        DIRECT VIEW

         CALL CONSW
     +     (ERROR, TEXT,
     +      FDIR2, FDIR1, ZWALL(I), ZWALL(I+1), RADIUS,
     +      ZERO, RADIUS)

C        REFLECTED VIEW

         CALL CONSW
     +     (ERROR, TEXT,
     +      FREFL2, FREFL1, (ZLEN2 - ZWALL(I+1)), (ZLEN2 - ZWALL(I)),
     +      RADIUS, ZERO, RADIUS)

         FESE(I) = FDIR1 + REFLEC*FREFL1
         FSEE(I) = FDIR2 + REFLEC*FREFL2

300   CONTINUE

C///  CONFIGURATIONS FACTORS BETWEEN ELEMENTS ON THE REACTOR SIDE WALL
C     IN THE REACTOR ENTRANCE.

C     SUBROUTINE CONSES
C    +  (ERROR, TEXT,
C    +   F, PDIM, RADIUS, REFLEC, ZPNTS, ZWALL)

      CALL CONSES
     +  (ERROR, TEXT,
     +   FSES, PDIM, RADIUS, REFLEC, ZPNTS, ZWALL)

C///  SELF CONFIGURATION FACTOR FOR REACTOR END WALL

C     SUBROUTINE CONWW
C    +  (ERROR, TEXT,
C    +   A, F12, F21, RIN1, RIN2, ROUT1, RO
      CALL CONWW
     +  (ERROR, TEXT,
     +   ZLEN2, FREFL1, DUMMY, ZERO, ZERO, RADIUS,
     +   RADIUS)

      FEE = REFLEC*FREFL1

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE CONSES
     +  (ERROR, TEXT,
     +   F, PDIM, RADIUS, REFLEC, ZPNTS, ZWALL)

C///////////////////////////////////////////////////////////////////////
C
C     CONSES
C
C///  COMPUTES THE CONFIGUATION FACTORS FROM THE RING AREAS ON THE
C     INTERIOR SIDE OF THE REACTOR ENTRANCE BETWEEN THE END WALL AND THE
C     FIRST WAFER. THE FIRST RING IS NEXT TO THE END WALL AND RING ZPNTS
C     IS NEXT TO THE FIRST WAFER.  ASSUMES ZWALL(ZPNTS+1) LIES AT THE
C     Z LOCATION OF THE FACE OF THE FIRST WAFER.
C
C     FORMULA C-82 FROM "A CATALOG OF RADIATION CONFIGURATION
C     FACTORS",J.R. HOWELL, IS USED IN CONJUNCTION WITH CONFIGURATION
C     FACTOR ALGEBRA.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     RADIUS  INPUT REAL - RADIUS OF CYLINDER.
C
C     REFLEC  INPUT REAL - REFLECTIVITY OF BOUNDING WAFER
C
C     ZPNTS   INPUT REAL - NUMBER OF POINTS ALONG THE WALL BETWEEN
C             REACTOR END WALL AND FIRST WAFER.
C
C     ZWALL   INPUT REAL DIMENSIONED (ZPNTS+1) - AXIAL LOCATIONS OF
C             CELL WALLS ALONG REACTOR
C                          ZWALL(1) = 0.0,
C                          ZWALL(ZPNTS+1) = LOCATION OF FIRST WAFER.
C
C     F       OUTPUT REAL DIMENSIONED (PDIM, ZPNTS) - F(I,J) IS THE
C             SHAPE FACTOR FROM RING I TO RING J ON THE REACTOR WALL.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   I, J, NUMTOT, PDIM, TEXT, ZPNTS

      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   F, L1, L2, L3, PI, RADIUS, REFLEC, ZWALL

      PARAMETER (ID = 'CONSES:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)
      PARAMETER (NUMTOT = 150)

      DIMENSION F(PDIM, ZPNTS), ZWALL(ZPNTS+1)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE
C
C///////////////////////////////////////////////////////////////////////

C///  CHECKS THE ARGUMENTS.

      ERROR = (RADIUS .LE. 0.0)
      IF (ERROR) GO TO 9001

C///////////////////////////////////////////////////////////////////////
C
C     (2) FIND THE CONFIGURATION FACTORS.
C
C///////////////////////////////////////////////////////////////////////

      DO 200 I = 1, ZPNTS
         DO 200 J = I, ZPNTS

C           DIRECT EXCHANGE COMPONENT

            IF (I .EQ. J) THEN
               F(I,J) = 1. + (ZWALL(J+1) - ZWALL(J))/(2.*RADIUS)
     +                  - SQRT(1. + ((ZWALL(J+1)
     +                  - ZWALL(J))**2)/(4.*RADIUS**2))
            ELSE
               L1 = (ZWALL(J+1) - ZWALL(J))/RADIUS
               L2 = (ZWALL(J+1) - ZWALL(I+1))/RADIUS
               L3 = (ZWALL(J+1) - ZWALL(I))/RADIUS
               F(I,J) = 1./(4.*(L3 - L2))*(2.*L1*(L3 - L2) +
     +                  (L3 - L1)*SQRT((L3 - L1)**2 + 4.) -
     +                  (L2 - L1)*SQRT((L2 - L1)**2 + 4.) -
     +                  L3*SQRT(L3**2 + 4.) +
     +                  L2*SQRT(L2**2 + 4.))
            ENDIF

C           COMPUTE AND ADD THE REFLECTED COMPONENT TO THE DIRECT

            L1 = (ZWALL(J+1) - ZWALL(J))/RADIUS
            L2 = ((2.*ZWALL(ZPNTS+1) - ZWALL(I+1)) - ZWALL(J))/RADIUS
            L3 = ((2.*ZWALL(ZPNTS+1) - ZWALL(I)) - ZWALL(J))/RADIUS
            F(I,J) = F(I,J) +
     +               REFLEC*(1./(4.*(L3 - L2))*(2.*L1*(L3 - L2) +
     +               (L3 - L1)*SQRT((L3 - L1)**2 + 4.) -
     +               (L2 - L1)*SQRT((L2 - L1)**2 + 4.) -
     +               L3*SQRT(L3**2 + 4.) +
     +               L2*SQRT(L2**2 + 4.)))

C           APPLY RECIPROCITY TO FIND F(J,I) FROM F(I,J)

            F(J,I) = (ZWALL(I+1) - ZWALL(I))/(ZWALL(J+1) -
     +               ZWALL(J))*F(I,J)

200   CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, RADIUS

      GO TO 99999

99001 FORMAT
     +   (/1X, A, 'ERROR. RADIUS IS LESS THAN OR EQUAL TO ZERO.',
     +   //10X, 'RADIUS = ',E10.3)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE CONSS
     +  (ERROR, TEXT,
     +   F11, RADIUS, REFLEC, SPACE)

C///////////////////////////////////////////////////////////////////////
C
C     CONSS
C
C///  COMPUTES THE CONFIGUATION FACTOR FROM A RING AREA ON THE INTERIOR
C     SIDE OF A RIGHT CIRCULAR CYLINDER TO ITSELF INCLUDING SPECULAR
C     REFLECTIONS WITH BOUNDING WAFERS. THE CONFIGURATION FACTOR FOR
C     DIRECT EXCHANGE CAN BE FOUND BY SETTING REFLEC = 0.
C
C     FORMULAS C-74 AND C-82 FROM "A CATALOG OF RADIATION CONFIGURATION
C     FACTORS",J.R. HOWELL, ARE USED IN CONJUNCTION WITH CONFIGURATION
C     FACTOR ALGEBRA.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     SPACE   INPUT REAL - LENGTH OF RING ELEMENT ON CYLINDER WALL
C
C     RADIUS  INPUT REAL - RADIUS OF CYLINDER.
C
C     REFLEC  INPUT REAL - REFLECTIVITY OF BOUNDING WAFERS
C
C     F11     OUTPUT REAL - SELF CONFIGURATION FRACTOR FOR RING
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   K, NUMTOT, TEXT

      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   F11, L1, L2, L3, PI, RADIUS, REFLEC, SPACE, SUM

      PARAMETER (ID = 'CONSW:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)
      PARAMETER (NUMTOT = 150)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE
C
C///////////////////////////////////////////////////////////////////////

C///  CHECKS THE ARGUMENTS.

      ERROR = (RADIUS .LE. 0.0)
      IF (ERROR) GO TO 9001

C///////////////////////////////////////////////////////////////////////
C
C     (2) FIND THE CONFIGURATION FACTORS.
C
C///////////////////////////////////////////////////////////////////////

C///  DIRECT EXCHANGE FACTOR

      SUM = 1. + SPACE/(2.*RADIUS) - SQRT(1. + SPACE**2/(4.*RADIUS**2))

C///  LOOP THROUGH AND SUM EXCHANGE WITH REFLECTED IMAGES

      DO 100 K = 1, NUMTOT
          L1 = SPACE/RADIUS
          L2 = SPACE*FLOAT(K)/RADIUS
          L3 = SPACE*FLOAT(K+1)/RADIUS
          SUM = SUM + 2.*(REFLEC**K)/(4.*(L3 - L2))*(2.*L1*(L3 - L2) +
     +          (L3 - L1)*SQRT((L3 - L1)**2 + 4.) -
     +          (L2 - L1)*SQRT((L2 - L1)**2 + 4.) -
     +           L3*SQRT(L3**2 + 4.) +
     +           L2*SQRT(L2**2 + 4.))
100   CONTINUE

C///  FINAL CONFIGURATION FACTOR

      F11 = SUM


C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, RADIUS

      GO TO 99999

99001 FORMAT
     +   (/1X, A, 'ERROR. RADIUS IS LESS THAN OR EQUAL TO ZERO.',
     +   //10X, 'RADIUS = ',E10.3)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE CONSW
     +  (ERROR, TEXT,
     +   F12, F21, LLOWER, LUPPER, RADIUS, RIN, ROUT)

C///////////////////////////////////////////////////////////////////////
C
C     CONSW
C
C///  COMPUTES THE CONFIGUATION FACTOR BETWEEN A FINITE RING AREA ON
C     THE INTERIOR SIDE OF A RIGHT CIRCULAR CYLINDER TO AN ANNULAR RING
C     AREA IN THE BASE OF THE CYLINDER.  IN THE LIMIT AS RIN = 0,
C     THE CONFIGURATION FACTOR FROM THE RING ELEMENT ON THE CYLINDER
C     SIDE TO A DISK IN THE BASE OF THE CYLINDER IS COMPUTED.
C
C     AREA1 = RING AREA ON INTERIOR OF THE CYLINDER.
C     AREA2 = ANNULAR RING IN THE BASE OF THE CYLINDER.
C
C     FORMULA C-77 FROM "A CATALOG OF RADIATION CONFIGURATION FACTORS",
C     J.R. HOWELL, IS USED IN CONJUNCTION WITH CONFIGURATION FACTOR
C     ALGEBRA.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     LLOWER  INPUT REAL - DISTANCE FROM CYLINDER BASE TO LOWER EDGE OF
C                          RING ON INTERIOR OF CYLINDER.
C
C     LUPPER  INPUT REAL - DISTANCE FROM CYLINDER BASE TO UPPER EDGE OF
C                          RING ON INTERIOR OF CYLINDER.
C
C     RADIUS  INPUT REAL - RADIUS OF CYLINDER.
C
C     RIN     INPUT REAL - INNER RADIUS OF ANNULAR RING IN BASE OF
C                          CYLINDER.
C
C     ROUT    INPUT REAL - OUTER RADIUS OF ANNULAR RING IN BASE OF
C                          CYLINDER.
C
C     F12     OUTPUT REAL - SHAPE FACTOR FROM RING ON CYLINDER SIDE TO
C                           ANNULAR RING IN BASE OF CYLINDER.
C
C     F21     OUTPUT REAL - SHAPE FACTOR FROM ANNULAR RING IN BASE OF
C                           CYLINDER TO RING ON CYLINDER SIDE.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   TEXT

      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   A1, A2, F12, F21, FINNER, FOUTER, LLOWER, LUPPER, PI, RADIUS,
     +   RIN, ROUT

      PARAMETER (ID = 'CONSW:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE
C
C///////////////////////////////////////////////////////////////////////

C///  CHECKS THE ARGUMENTS.

      ERROR = (RIN .GE. ROUT)
      IF (ERROR) GO TO 9001

      ERROR = (LLOWER .GE. LUPPER)
      IF (ERROR) GO TO 9002

      ERROR = ((RIN .GE. RADIUS) .OR. (ROUT .GT. RADIUS))
      IF (ERROR) GO TO 9003

C///////////////////////////////////////////////////////////////////////
C
C     (2) FIND THE CONFIGURATION FACTORS.
C
C///////////////////////////////////////////////////////////////////////

      IF ((ROUT .EQ. 0.0) .OR. (LUPPER .EQ. LLOWER)) THEN
         F12 = 0.0
         F21 = 0.0
         RETURN
      ELSE
      ENDIF

C///  CONFIGURATION FACTOR TO THE OUTER DISK AND INNER DISK

         FOUTER = 1./(4.*RADIUS/ROUT*(LUPPER/ROUT - LLOWER/ROUT))*
     +            (((LLOWER/ROUT)**2 - (LUPPER/ROUT)**2) -
     +            SQRT(((LLOWER/ROUT)**2 +
     +            (RADIUS/ROUT)**2 + 1.)**2 - 4.*(RADIUS/ROUT)**2) +
     +            SQRT(((LUPPER/ROUT)**2 +
     +            (RADIUS/ROUT)**2 + 1.)**2 - 4.*(RADIUS/ROUT)**2))

      IF (RIN .EQ. 0.0) THEN
         FINNER = 0.0
      ELSE
         FINNER = 1./(4.*RADIUS/RIN*(LUPPER/RIN - LLOWER/RIN))*
     +            (((LLOWER/RIN)**2 - (LUPPER/RIN)**2) -
     +            SQRT(((LLOWER/RIN)**2 +
     +            (RADIUS/RIN)**2 + 1.)**2 - 4.*(RADIUS/RIN)**2) +
     +            SQRT(((LUPPER/RIN)**2 +
     +            (RADIUS/RIN)**2 + 1.)**2 - 4.*(RADIUS/RIN)**2))
      ENDIF

C///  AREAS

      A1 = 2.*RADIUS*PI*(LUPPER - LLOWER)
      A2 = PI*(ROUT**2 - RIN**2)

C///  FINAL CONFIGURATION FACTORS

      F12 = FOUTER - FINNER
      F21 = A1/A2*F12

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, RIN, ROUT

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, LLOWER, LUPPER

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, RADIUS, RIN, ROUT

      GO TO 99999

99001 FORMAT
     +   (/1X, A, 'ERROR. INNER RADIUS IS GREATER THAN OR EQUAL',
     +    /10X, 'TO OUTER RADIUS.',
     +   //10X, 'RIN = ',E10.3,
     +    /10X, 'ROUT = ',E10.3)
99002 FORMAT
     +   (/1X, A, 'ERROR. DISTANCE TO LOWER EDGE OF RING ON CYLINDER',
     +    /10X, 'IS GREATER THAN OR EQUAL TO DISTANCE TO UPPER EDGE.',
     +   //10X, 'LLOWER = ',E10.3,
     +    /10X, 'LUPPER = ',E10.3)

99003 FORMAT
     +   (/1X, A, 'ERROR. INNER OR OUTER RADIUS OF ANNULAR RING IS'
     +    /10X, 'GREATER THAN OR EQUAL THE RADIUS OF THE CYLINDER.',
     +   //10X, 'RADIUS = ',E10.3,
     +    /10X, 'RIN = ',E10.3,
     +    /10X, 'ROUT = ',E10.3)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE CONWAF
     +  (ERROR, TEXT,
     +   FSS, FSW, FWS, FWW, REFLEC, RPNTS, RWALL, SPACE)

C///////////////////////////////////////////////////////////////////////
C
C     CONWAF
C
C///  COMPUTES THE CONFIGUATION FACTORS BETWEEN TWO WAFERS, THE WAFERS
C     AND THE SIDE WALL, AND THE SELF CONFIGURATION FACTOR FROM THE SIDE
C     WALL AND ITSELF. SPECULAR REFLECTION BETWEEN THE WAFERS IS
C     CONSIDERED.  DIFFUSE EXCAHNGE FACTORS CAN BE CALCULATED BY SETTING
C     REFLEC = 0.0.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     REFLEC  INPUT REAL - REFLECTIVITY OF WAFERS.
C
C     RPNTS   INPUT REAL - NUMBER OF RADIAL GRID POINTS ON A WAFER.
C
C     RWALLS  INPUT REAL DIMENSIONED (RPNTS+1) - RADIAL LOCATIONS OF
C             CELL WALLS ON WAFER.
C
C     SPACE   INPUT REAL - SPACING BETWEEN WAFERS.
C
C     FSS     OUTPUT REAL - SELF SHAPE FACTOR FOR SIDE WALL.
C
C     FSW     OUTPUT REAL DIMENSIONED (2*RPNTS) - SHAPE FACTORS FROM
C             SIDE WALL TO WAFER ELEMENTS.
C
C     FWS     OUTPUT REAL DIMENSIONED (2*RPNTS) - SHAPE FACTORS FROM
C             WAFER ELEMENTS TO SIDE WALL.
C
C     FWW     OUTPUT REAL DIMENSIONED (2*RPNTS, 2*RPNTS) - SHAPE
C             FACTORS BETWEEN WAFER ELEMENTS.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   RPNTS, I, J, TEXT

      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   FSS, FSW, FWS, FWW, PI, REFLEC, RWALL, SPACE

      PARAMETER (ID = 'CONWAF:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)

      DIMENSION  FSW(2*RPNTS), FWS(2*RPNTS), FWW(2*RPNTS, 2*RPNTS),
     +           RWALL(RPNTS+1)

C///////////////////////////////////////////////////////////////////////
C
C     (1)  CONFIGURATION FACTORS
C
C///////////////////////////////////////////////////////////////////////

      DO 200 I = 1, RPNTS
         DO 100 J = 1, RPNTS

C///     WAFER TO WAFER

C           SUBROUTINE SCONWW
C          +  (ERROR, TEXT,
C          +   ASTAR, F12, F21, REFLEC, RIN1, RIN2, ROUT1, ROUT2, SAME)

C///        WITHIN THE SAME WAFER

            CALL SCONWW
     +        (ERROR, TEXT,
     +         SPACE, FWW(I,J), FWW(J,I), REFLEC, RWALL(I), RWALL(J),
     +         RWALL(I+1), RWALL(J+1), .TRUE.)

            FWW(RPNTS + I, RPNTS + J) = FWW(I, J)

C///        BETWEEN ADJACENT WAFERS

            CALL SCONWW
     +        (ERROR, TEXT,
     +         SPACE, FWW(I, RPNTS + J), FWW(RPNTS + J, I), REFLEC,
     +         RWALL(I), RWALL(J), RWALL(I+1), RWALL(J+1), .FALSE.)

100      CONTINUE

C///     SIDE WALL TO WAFER

C        SUBROUTINE SCONSW
C       +  (ERROR, TEXT,
C       +   F12, F21, RADIUS, REFLEC, RIN, ROUT, SPACE)

         CALL SCONSW
     +     (ERROR, TEXT,
     +      FSW(I), FWS(I), RWALL(RPNTS+1), REFLEC, RWALL(I),
     +      RWALL(I+1), SPACE)

         FSW(RPNTS + I) = FSW(I)
         FWS(RPNTS + I) = FWS(I)

200   CONTINUE

C///  SIDE WALL SELF CONFIGURATION FACTOR

C     SUBROUTINE CONSS
C    +  (ERROR, TEXT,
C    +   F11, RADIUS, REFLEC, SPACE)

      CALL CONSS
     +  (ERROR, TEXT,
     +   FSS, RWALL(RPNTS + 1), REFLEC, SPACE)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE CONWW
     +  (ERROR, TEXT,
     +   A, F12, F21, RIN1, RIN2, ROUT1, ROUT2)

C///////////////////////////////////////////////////////////////////////
C
C     CONWW
C
C///  COMPUTES THE CONFIGUATION FACTOR BETWEEN TWO PARALLEL, DIRECTLY
C     APPOSED PLANE RING AREAS. IN THE LIMIT AS RIN1 AND/OR RIN2 GOES
C     TO ZERO, THE CONFIGURATION FACTOR BETWEEN A RING AND A DISK OR
C     TWO DISKS IS COMPUTED.
C
C
C     CONFIGURATION FACTOR CALCULATED USING SHAPE FACTOR ALGEBRA AND
C     REFERENCE C-35 FROM "A CATALOG OF RADIATION CONFIGURATION
C     FACTORS", J.R. HOWELL.
C
C     WRITTEN BY: W.G. HOUF
C                 SANDIA NATIONAL LABORATORIES
C                 LIVERMORE, CA
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     A       INPUT REAL - ON AXIS DISTANCE BETWEEN TWO RINGS.
C
C     RIN1    INPUT REAL - INNER RADIUS OF RING 1.
C
C     RIN2    INPUT REAL - INNER RADIOUS OF RING 2.
C
C     ROUT1   INPUT REAL - OUTER RADIUS OF RING 1.
C
C     ROUT2   INPUT REAL - OUTER RADIUS OF RING 2.
C
C     F12     OUTPUT REAL - SHAPE FACTOR FROM RING 1 TO RING 2
C
C     F21     OUTPUT REAL - SHAPE FACTOR FROM RING 2 TO RING 1
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   TEXT

      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   A, ARING1, ARING2, F12, F21, FL1L2, FL1S2, FS1L2, FS1S2,
     +   PI, RIN1, RIN2, ROUT1, ROUT2, SMALL

      PARAMETER (ID = 'CONWW:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)
      PARAMETER (SMALL = 1.0E-30)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////


C///  CHECKS

      ERROR = (RIN1 .GE. ROUT1)
      IF (ERROR) GO TO 9001

      ERROR = (RIN2 .GE. ROUT2)
      IF (ERROR) GO TO 9002

      ERROR = (A .LE. 0.0)
      IF (ERROR) GO TO 9003

C///////////////////////////////////////////////////////////////////////
C
C     (2) FIND THE CONFIGURATION FACTORS
C
C///////////////////////////////////////////////////////////////////////

      IF ((ROUT1 .EQ. 0.0) .OR. (ROUT2 .EQ. 0.0)) THEN
         F12 = 0.0
         F21 = 0.0
         RETURN
      ELSE
      ENDIF

C///  CONFIGURATION FACTOR FROM LARGE DISK1 TO LARGE DISK2 (FL1L2) AND
C     LARGE DISK1 TO SMALL DISK2 (FL1S1)

         FL1L2 = 0.5*(1. + (1. + (ROUT2/A)**2)/(ROUT1/A)**2 -
     +           SQRT((1. + (1. + (ROUT2/A)**2)/(ROUT1/A)**2)**2 -
     +           4.*(ROUT2/ROUT1)**2))
         FL1S2 = 0.5*(1. + (1. + (RIN2/A)**2)/(ROUT1/A)**2 -
     +           SQRT((1. + (1. + (RIN2/A)**2)/(ROUT1/A)**2)**2 -
     +           4.*(RIN2/ROUT1)**2))


C///  CONFIGURATION FACTOR FROM SMALL DISK1 TO LARGE DISK2 (FS1L1) AND
C     SMALL DISK1 TO SMALL DISK2 (FS1S1)

      IF (RIN1 .LE. SMALL) THEN
         FS1L2 = 0.0
         FS1S2 = 0.0
      ELSE
         FS1L2 = 0.5*(1. + (1. + (ROUT2/A)**2)/(RIN1/A)**2 -
     +           SQRT((1. + (1. + (ROUT2/A)**2)/(RIN1/A)**2)**2 -
     +           4.*(ROUT2/RIN1)**2))
         FS1S2 = 0.5*(1. + (1. + (RIN2/A)**2)/(RIN1/A)**2 -
     +           SQRT((1. + (1. + (RIN2/A)**2)/(RIN1/A)**2)**2 -
     +           4.*(RIN2/RIN1)**2))
      ENDIF

C///  AREAS

      ARING1 = PI*(ROUT1**2 - RIN1**2)
      ARING2 = PI*(ROUT2**2 - RIN2**2)

C///  FINAL CONFIGURATION FACTORS

      F12 = (PI*RIN1**2 + ARING1)/ARING1*(FL1L2 - FL1S2) -
     +      (PI*RIN1**2)/ARING1*(FS1L2 - FS1S2)

      F21 = ARING1/ARING2*F12

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, RIN1, ROUT1

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, RIN2, ROUT2

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, A

      GO TO 99999

99001 FORMAT
     +   (/1X, A, 'ERROR. RIN1 GREATER THAN OR EQUAL TO ROUT1.'
     +   //10X, 'RIN1 = ',E10.3,
     +    /10X, 'ROUT1 = ',E10.3)

99002 FORMAT
     +   (/1X, A, 'ERROR. RIN2 GREATER THAN OR EQUAL TO ROUT2.'
     +   //10X, 'RIN2 = ',E10.3,
     +    /10X, 'ROUT2 = ',E10.3)

99003 FORMAT
     +   (/1X, A, 'ERROR. SPACING BETWEEN DISKS IS ZERO OR LESS.'
     +   //10X, 'A = ', E10.3)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE EBAL1
     +  (ERROR, TEXT,
     +   AEX, AIN, AREA, CFLOW, COMPS, COND, CONDUC, CP, DZ, DZEND,
     +   DZENDI, EMIS, EMISS, EMISSE, EMISSI, EMISSW, F12, F21, F22,
     +   F2H, FEE, FEG, FENDS, FESE, FEW, FGE, FGSE, FH2, FHH, FSEE,
     +   FSEG, FSES, FSEW, FWE, FWSE, GBOT, GROUPA, GROUPB, GTOP,
     +   HCOEF, IBEG, IEND, IPVT3, IPVT4, J1, JBOT, JEND, JEX,
     +   JIN, JTOP, JWALL, LOS1, LOS2, LOS3, MAXHT, MFLOW, MINJE,
     +   MWZ, NUMEX, NUMHT, NUMT, NUMTC, NUMWF, PCAP, PDIM, PEX,
     +   PIN, PMAX, POWERH, PWALL, QWAF, R, REFL, RHO, RIN, RINI,
     +   RINJEC,
     +   RINO, RINQ, ROUTQ, RPNTS, RTC, RTUBE, RTUBEI, RTUBEO, SUMB,
     +   SUMT, TAMB, TAMB4, TC, TCAP, TCAPS, TEMMAX, TEND, TEND0,
     +   TEND4, TENDI, TENDI0, TENDI4, TEX, TEX0, TEX4, TFIX, THICKI,
     +   TIME, TRAN, TSIDE, TSTEP, TTGAS, TURB, TWAF, TWAF0, TWAF4,
     +   TX, TX0, TX4, WR, WZ, Z, ZBEG, ZDAT, ZEND, ZLEN, ZTC, ZWAF)

C///////////////////////////////////////////////////////////////////////
C
C     EBAL1
C
C///  EVALUATE THE RESIDUAL OF THE ENERGY EQUATIONS FOR THE ENTRANCE
C     AND EXIT.
C
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     INPUT:
C     (SEE SUBROUTINE TWAF1)
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   AEXDIM, AINDIM, COMPS, COUNT, FLOW, GROUPA, GROUPB, I,
     +   IBEG, IEND, IOPT, IPVT3, IPVT4, J, K, KRON, LOC, MAXHT,
     +   NUMEX, NUMHT, NUMT, NUMTC, NUMWF, PDIM, PEX, PIN, PMAX, RPNTS,
     +   TEMMAX, TEXT, WAFER, WALL, ZPNTS
      LOGICAL CFLOW, ERROR, PCON, PED, RINJEC, TCAPS, TFIX, TIME, TURB
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   ACOND, AIN, AINN, AEX, AOUT, AREA, AREAF, COND,
     +   CONDE, CONDUC, CONVN, CONVS, CP, CPAIR, DH, DUMMY, DZ,
     +   DZEND, DZENDI, EEAST, EGAP, EINJ, EWEST, EMIS, EMISS,
     +   EMISSE, EMISSI, EMISSW, F12, F21, F22, F2H, FACTE, FBOT,
     +   FEE, FEG, FENDS, FESE, FEW, FGE, FGSE, FH2, FHH, FSEE, FSEG,
     +   FSES, FSEW, FTOP, FWE, FWSE, GBOT, GEND, GGAP, GTOP, GWALL,
     +   HCOEF, HCOF, HINN, HOUT, J1, JBOT, JEND, JEX, JIN, JTOP,
     +   JWALL, KAIR, LOS1, LOS2, LOS3, MFLOW, MFLOWE, MFLOWW, MINJE,
     +   MWZ, MUAIR, NU, PCAP,
     +   PEDRI, PEDRO,
     +   PI, POWERH, PRAIR, PWALL, QCONV, QEAST,
     +   QEDGE, QFLUX, QNORTH, QSOUTH, QWAF, QWEST, R, RCPV,
     +   RE, REFL,
     +   REYN, RHO, RHOAIR, RIN, RINI,
     +   RINO, RINQ, ROUTQ, RTC, RTUBE, RTUBEI, RTUBEO, SIG,
     +   STEADY, SUMB, SUMSE, SUMSE1, SUMT, SUMW, SUMWE, TAIR,
     +   TAIRIN, TAMB, TAMB4, TC, TCAP, TEAST1,
     +   TEND, TEND0, TEND4, TENDI, TENDI0, TENDI4, TEX, TEX0, TEX4,
     +   TFILM, THICKI, TTHIC,
     +   TRAN, TREF, TSIDE, TSTEP, TTGAS, TWAF, TWAF0,
     +   TWAF4,
     +   TWEST1, TX, TX0, TX4, VOL, WR, WZ, Z, ZBEG, ZDAT, ZEND, ZLEN,
     +   ZLOWER, ZTC, ZUPPER, ZWAF,
     +   TMASS

      PARAMETER (ID = ' EBAL1:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)
      PARAMETER (QINNER = 0, QOUTER = 1)
      PARAMETER (SIG = 5.6696E-12)

      DIMENSION
     +   AEX(5*PEX+2, 5*PEX+2), AIN(5*PIN+2, 5*PIN+2), AINN(2),
     +   AOUT(2), AREA(2), AREAF(2), CFLOW(2), CONDUC(8), CONVN(2),
     +   CONVS(2), CP(8), DH(2), DZEND(2),
     +   DZENDI(2), EGAP(2), EMIS(2), EMISSI(5), FEE(2), FEG(2),
     +   FENDS(PDIM+1, NUMEX+1, 2), FESE(PDIM, 2), FEW(RPNTS, 2),
     +   FGE(2), FGSE(PDIM, 2), FSEE(PDIM, 2), FSEG(PDIM, 2),
     +   FSES(PDIM, PDIM, 2), FSEW(PDIM, RPNTS, 2), FWE(RPNTS, 2),
     +   FWSE(RPNTS, PDIM, 2), GBOT(NUMTC), GGAP(2), GTOP(NUMTC),
     +   HCOEF(5), HINN(2), HOUT(2), IBEG(3), IEND(3),
     +   IPVT3(5*PIN + 2),
     +   IPVT4(5*PEX + 2), J1(4, PDIM, 2), JEX(5*PEX+2), JIN(5*PIN+2),
     +   JBOT(NUMTC), JEND(2), JTOP(NUMTC), JWALL(PDIM, 2),
     +   LOS1(NUMWF), LOS2(PDIM, 2), LOS3(2), MFLOW(2),
     +   MINJE(PDIM, 2), MWZ(PDIM+1, 2), NU(2),
     +   PCAP(2), POWERH(MAXHT)

      DIMENSION
     +   PWALL(PDIM, 2), QEDGE(4, 2), QWAF(RPNTS, 2), R(RPNTS), RE(2),
     +   REFL(2),
     +   RHO(8), RINQ(2), ROUTQ(2), SUMB(NUMTC), SUMT(NUMTC),
     +   TAIRIN(2), TC(NUMTC),
     +   TCAP(2), TCAPS(2), TEND(2), TEND0(2), TEND4(2), TENDI(2),
     +   TENDI0(2), TENDI4(2), TEX(NUMWF, NUMEX), TEX0(NUMWF, NUMEX),
     +   TEX4(NUMWF, NUMEX), TRAN(2), TSIDE(PDIM, 2), TTGAS(2),
     +   TWAF(NUMWF+1, RPNTS), TWAF0(NUMWF+1, RPNTS),
     +   TWAF4(NUMWF+1, RPNTS), TX(PDIM, NUMEX, 2),
     +   TX0(PDIM, NUMEX, 2), TX4(PDIM, NUMEX, 2), TURB(2), VOL(2),
     +   WR(RPNTS+1), WZ(PDIM+1, 2), Z(PDIM,2),
     +   ZBEG(MAXHT), ZDAT(TEMMAX + 2), ZEND(MAXHT), ZTC(NUMTC),
     +   ZWAF(NUMWF)

      COMMON /PED1/ PEDRI, PEDRO, TTHIC
      COMMON /PED2/ PCON, PED

C///////////////////////////////////////////////////////////////////////
C
C     (1) ENERGY BALANCE TERMS FOR ENTRANCE AND EXIT REGION.
C
C///////////////////////////////////////////////////////////////////////


C///  DEFINE STATEMENT FUNCTION FOR KRONECKER DELTA

      KRON(I, J) = MAX(-1, -ABS(I - J)) + 1

C///  DEFINE STATEMENT FUNCTION FOR VISCOSITY OF AIR (G/(CM SEC)).
C     TAIR - TEMPERATURE OF AIR IN K.

      MUAIR(TAIR) = 1.716E-04*383.66/(TAIR + 110.56)*(TAIR/273.1)**1.5

C///  DEFINE STATEMENT FUNCTION FOR CONDUCTIVITY OF AIR (W/(CM K)).

      KAIR(TAIR) = 2.414E-04*467.5/(TAIR + 194.4)*(TAIR/273.1)**1.5

C///  DEFINE STATEMENT FUNCTION FOR DENSITY OF AIR AT 1 ATM (G/CM**3).

      RHOAIR(TAIR) = 1.0133E02/287./TAIR

C///  SET THE SPECIFIC HEAT OF AIR (J/G/K)

      CPAIR = 1.189

C///  PRANDTL NUMBER OF AIR

      PRAIR = 0.719

C///  SET GAP TREATMENT OPTION

      IOPT = 1

C*****PRINT SEPARATION BARS
C      WRITE(TEXT, *) ''
C      WRITE(TEXT, *) '--------- ENTRANCE AND EXIT ------------------'
C      WRITE(TEXT, *) ' '
C*****END PRINT SEPARATION BARS

C-----NEW ENTRANCE AND EXIT MATRICES

      LOC = 1
      WAFER = 1
      ZPNTS = PIN
      AINDIM = 5*PIN+2

      COUNT = 0
      DO 200 I = 1, ZPNTS
         SUMW = 0.0
         DO 100 J = 1, RPNTS
            SUMW = SUMW + EMISS*FSEW(I, J, LOC)*TWAF4(WAFER, J)
100      CONTINUE

         JIN(COUNT + 1) = EMISSW*SIG*TX4(I,3,LOC)
         JIN(COUNT + 2) = EMIS(2)*SIG*TX4(I,2,LOC)
         JIN(COUNT + 3) = EMIS(2)*SIG*TX4(I,2,LOC)
         JIN(COUNT + 4) = EMIS(1)*SIG*TX4(I,1,LOC) + TRAN(1)*SIG*SUMW
         JIN(COUNT + 5) = EMIS(1)*SIG*TX4(I,1,LOC) + REFL(1)*SIG*SUMW
         COUNT = COUNT + 5
200   CONTINUE
      SUMWE = 0.0
      DO 300 J = 1, RPNTS
         SUMWE = SUMWE + EMISS*FEW(J, LOC)*TWAF4(WAFER, J)
300   CONTINUE
      JIN(5*PIN+1) = EMISSE*SIG*TEND4(LOC) + (1. - EMISSE)*SIG*SUMWE
      JIN(5*PIN+2) = 0.0

C*****DOUBLE PRECISION
      CALL DGESL(AIN, AINDIM, AINDIM, IPVT3, JIN, 0)
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      CALL SGESL(AIN, AINDIM, AINDIM, IPVT3, JIN, 0)
C*****END SINGLE PRECISION

C     -- EXIT SECTION --

      LOC = 2
      WAFER = NUMWF
      ZPNTS = PEX
      AEXDIM = 5*PEX+2

      COUNT = 0
      DO 500 I = 1, ZPNTS
         SUMW = 0.0
         DO 400 J = 1, RPNTS
            SUMW = SUMW + EMISS*FSEW(I, J, LOC)*TWAF4(WAFER, J)
400      CONTINUE

         JEX(COUNT + 1) = EMISSW*SIG*TX4(I,3,LOC)
         JEX(COUNT + 2) = EMIS(2)*SIG*TX4(I,2,LOC)
         JEX(COUNT + 3) = EMIS(2)*SIG*TX4(I,2,LOC)
         JEX(COUNT + 4) = EMIS(1)*SIG*TX4(I,1,LOC) + TRAN(1)*SIG*SUMW
         JEX(COUNT + 5) = EMIS(1)*SIG*TX4(I,1,LOC) + REFL(1)*SIG*SUMW
         COUNT = COUNT + 5
500   CONTINUE
      SUMWE = 0.0
      DO 600 J = 1, RPNTS
         SUMWE = SUMWE + EMISS*FEW(J, LOC)*TWAF4(WAFER, J)
600   CONTINUE
      JEX(5*PEX+1) = EMISSE*SIG*TEND4(LOC) + (1. - EMISSE)*SIG*SUMWE
      JEX(5*PEX+2) = 0.0

C*****DOUBLE PRECISION
      CALL DGESL(AEX, AEXDIM, AEXDIM, IPVT4, JEX, 0)
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      CALL SGESL(AEX, AEXDIM, AEXDIM, IPVT4, JEX, 0)
C*****END SINGLE PRECISION

C     PACK JIN AND JEX IN J1(4, PDIM, 2) ARRAY

      DO 900 LOC = 1, 2
      IF (LOC.EQ.1) THEN
         COUNT = 0
         DO 700 K = 1, PIN
            J1(1, K, LOC) = JIN(COUNT + 1)
            J1(2, K, LOC) = JIN(COUNT + 2)
            J1(3, K, LOC) = JIN(COUNT + 3)
            J1(4, K, LOC) = JIN(COUNT + 4)
            JWALL(K, LOC) = JIN(COUNT + 5)
            COUNT = COUNT + 5
700      CONTINUE
         JEND(LOC) = JIN(5*PIN + 1)
         EGAP(LOC) = JIN(5*PIN + 2)
      ELSE
         COUNT = 0
         DO 800 K = 1, PEX
            J1(1, K, LOC) = JEX(COUNT + 1)
            J1(2, K, LOC) = JEX(COUNT + 2)
            J1(3, K, LOC) = JEX(COUNT + 3)
            J1(4, K, LOC) = JEX(COUNT + 4)
            JWALL(K, LOC) = JEX(COUNT + 5)
            COUNT = COUNT + 5
800      CONTINUE
         JEND(LOC) = JEX(5*PEX + 1)
         EGAP(LOC) = JEX(5*PEX + 2)
      ENDIF
900   CONTINUE

C///  DEFINE AREAS FOR FLOW AND OTHER PARAMETERS FOR
C     CONVECTIVE HEAT TRANSFER COOLING.

C     FLOW AREAS FOR ANNULLI.
      AREAF(1) = PI*(RINQ(2)**2 - ROUTQ(1)**2)
      AREAF(2) = PI*(RTUBEI**2 - ROUTQ(2)**2)

C     INSIDE AND OUTSIDE AREAS AND VOLUME TERMS OF ANNULLI.
      AOUT(1)  = 2.*PI*RINQ(2)
      AINN(1)  = 2.*PI*ROUTQ(1)
      AOUT(2)  = 2.*PI*RTUBEI
      AINN(2)  = 2.*PI*ROUTQ(2)
      VOL(1)   = PI*(RINQ(2)**2 - ROUTQ(1)**2)
      VOL(2)   = PI*(RTUBEI**2 - ROUTQ(2)**2)

C     HYDRAULIC DIAMETERS.
      DH(1) = 2.*(RINQ(2) - ROUTQ(1))
      DH(2) = 2.*(RTUBEI - ROUTQ(2))

C     SET REFERENCE TEMPERATURE FOR AIR VISCOSITY TO 1000K.
      TREF = 1000.0

C     REYNOLDS NUMBER CHECKED AT TURB(1) AND TURB(2) SET IN
C     TWAF1.

C-----END NEW ENTRANCE AND EXIT MATRICES

C///  TOP OF THE LOOP OVER THE ENTRANCE OR EXIT SECTION.

      DO 1120 LOC = 1, 2
         IF (LOC .EQ. 1) THEN
            ZPNTS = PIN
            WAFER = 1
         ELSE
            ZPNTS = PEX
            WAFER = NUMWF
         END IF

C///  IRRADIAION ON THE GAP BETWEEN THE EDGE OF THE END WAFER AND THE
C///  REACTOR WALL.

      GGAP(LOC) = JEND(LOC)*FGE(LOC)
      DO 1010 I = 1, ZPNTS
         GGAP(LOC) = GGAP(LOC) + JWALL(I, LOC)*FGSE(I, LOC)
1010  CONTINUE

      IF (IOPT .EQ. 1) THEN
         EGAP(LOC) = GGAP(LOC)
      ELSE
         EGAP(LOC) = 0.0
      END IF

C///  TOP OF THE LOOP OVER THE POINTS IN THE ENTRANCE OR EXIT.

      DO 1060 I = 1, ZPNTS
         SUMSE = 0.0
         DO 1020 K = 1, ZPNTS
            SUMSE = SUMSE + JWALL(K, LOC)*FSES(I, K, LOC)
1020     CONTINUE

         SUMW = 0.0
         DO 1030 J = 1, RPNTS
            SUMW = SUMW + EMISS*FSEW(I, J, LOC)*TWAF4(WAFER, J)
1030     CONTINUE

C*****PRINT QUARTZ RADIOSITIES
C      WRITE(TEXT, *) '--QUARTZ RADIOSITES--'
C      WRITE(TEXT, *) 'LOC = ', LOC
C      WRITE(TEXT, *) 'I = ',I,' J1(1) = ',J1(1, I, LOC),
C     +               ' J1(2) = ',J1(2, I, LOC),
C     +               ' J1(3) = ',J1(3, I, LOC),
C     +               ' J1(4) = ',J1(4, I, LOC)
C      WRITE(TEXT, *) 'JEND = ',JEND(LOC),' EGAP = ',EGAP(LOC)
C*****END PRINT QUARTZ RADIOSITIES

C--------------------THERMOCOUPLE SECTION --------------------------
C///  ACCUMULATE IRRADIANCE ON EACH THERMOCOUPLE OVER ENTRANCE
C     AND EXIT.

      DO 1040 K = 1, NUMTC

C        SUBROUTINE CONCC
C     +    (ERROR, TEXT,
C     +     F12, F21, R1, R2, Z1H, Z1L, Z2H, Z2L)

         IF (LOC .EQ. 1) THEN
            ZUPPER = WZ(I+1, LOC)
            ZLOWER = WZ(I, LOC)
         ELSE
            ZUPPER = ZLEN - WZ(I, LOC)
            ZLOWER = ZLEN - WZ(I+1, LOC)
         ENDIF

         CALL CONCC(ERROR, TEXT, FTOP, DUMMY, RTC, RINQ(2), ZTC(K)+DZ,
     +             ZTC(K), ZUPPER, ZLOWER)
         IF (ERROR) GO TO 9103

         SUMT(K) = SUMT(K) + FTOP

         CALL CONCC(ERROR, TEXT, DUMMY, FBOT, RINQ(1), RTC,
     +             ZUPPER, ZLOWER, ZTC(K)+DZ, ZTC(K))
         IF (ERROR) GO TO 9103

         SUMB(K) = SUMB(K) + FBOT

         GTOP(K) = GTOP(K) + J1(3, I, LOC)*FTOP
         GBOT(K) = GBOT(K) + J1(4, I, LOC)*FBOT

         IF((ZUPPER.LT.ZTC(K)).AND.
     +      (ZTC(K).LE.ZLOWER)) THEN
            JTOP(K) = J1(3, I, LOC)
            JBOT(K) = J1(4, I, LOC)
         ENDIF

C*****WRITE THERMOCOUPLE SHAPEFACTORS
C      WRITE(TEXT, *) ' THERMOCOUPLE SHAPE FACTORS'
C      WRITE(TEXT, *) ' OVER ENTRANCE AND EXIT'
C      WRITE(TEXT, *) 'LOC = ', LOC
C      WRITE(TEXT, *) 'I = ', I
C      WRITE(TEXT, *) 'TC NUMBER = ', K
C      WRITE(TEXT, *) 'FTOP = ', FTOP
C      WRITE(TEXT, *) 'FBOT = ', FBOT
C      WRITE(TEXT, *) 'FBOT = ', FBOT
C*****END WRITE THERMOCOUPLE SHAPEFACTORS

1040  CONTINUE
C-------------------------------------------------------------------
C*****OLD FENDS TERMS
C      FENDS(I, 1, LOC)
C     +      = JWALL(I, LOC)
C     +      - (REFL(1)*(SUMSE - JWALL(I, LOC)*FSES(I, I, LOC))
C     +      + TRAN(1)*J1(3)*F12
C     +      + REFL(1)*JEND(LOC)*FSEE(I, LOC)
C     +      + REFL(1)*SIG*SUMW
C     +      + REFL(1)*EMISS*EGAP(LOC)*FSEG(I, LOC)
C     +      + EMIS(1)*SIG*TX4(I, 1, LOC))
C     +      / (1. - REFL(1)*FSES(I, I, LOC))
C*****END OLD FENDS TERMS

C*****NEW FENDS TERMS
              SUMSE1 = 0.0
              DO 1045 K = 1, ZPNTS
              SUMSE1 = SUMSE1 +
     +                (KRON(I,K) - REFL(1)*FSES(I,K,LOC))*JWALL(K,LOC)
1045          CONTINUE
C       FENDS(I, 1, LOC)
C     +      = SUMSE1 +
C     +      - TRAN(1)*J1(3)*F12
C     +      - REFL(1)*JEND(LOC)*FSEE(I, LOC)
C     +      - REFL(1)*SIG*SUMW
C     +      - REFL(1)*EMISS*EGAP(LOC)*FSEG(I, LOC)
C     +      - EMIS(1)*SIG*TX4(I, 1, LOC)
C*****END NEW FENDS TERMS

C      EXPRESSION FOR JWALL BEFORE QUARTZ WALLS WERE ADDED.
C      FENDS(I, 1, LOC)
C     +      = JWALL(I, LOC)
C     +      - (EMISSW * SIG * TWAL4(I, LOC)
C     +      + (1.0 - EMISSW) * JEND(LOC) * FSEE(I, LOC)
C     +      + (1.0 - EMISSW) * (SUMSE - JWALL(I, LOC)*FSES(I, I, LOC))
C     +      + (1.0 - EMISSW) * SIG * SUMW
C     +      + (1.0 - EMISSW) * EMISS * EGAP(LOC) * FSEG(I, LOC))
C     +      / (1.0 - (1.0 - EMISSW) * FSES(I, I, LOC))

C///  ENERGY BALANCE FOR THE QUARTZ WALLS 1 AND 2.

      DO 1050 WALL = 1, 2
         IF (WALL. EQ. 1) THEN
            GWALL = JEND(LOC)*FSEE(I, LOC) + SUMSE + SIG*SUMW +
     +              EMISS*EGAP(LOC)*FSEG(I, LOC)
            QNORTH = J1(4, I, LOC) - J1(3, I, LOC)*F12
            QSOUTH = GWALL - JWALL(I, LOC)
C*****FIX HEAT FLUX ON OUTSIDE OF WALL 1 IN ENTRANCE AND EXIT
C            QNORTH = -200.0/(AREA(WALL)*(WZ(I+1, LOC) - WZ(I, LOC)))
C*****END FIX HEAT FLUX ON OUTSIDE OF WALL 1 IN ENTRANCE AND EXIT
         ELSE IF (WALL .EQ. 2) THEN
            QNORTH = J1(2, I, LOC) - J1(1, I, LOC)*F2H
            QSOUTH = J1(4, I, LOC)*F21 +
     +               J1(3, I, LOC)*F22 - J1(3, I, LOC)
         END IF

         IF (I .EQ. ZPNTS) THEN
            IF (LOC. EQ. 1) QEAST = -CONDUC(4)*(TEX(WAFER, WALL) -
     +           TX(I, WALL, LOC))/(ZWAF(WAFER) - Z(I, LOC))
            IF (LOC. EQ. 2) QEAST = -CONDUC(4)*(TEX(WAFER, WALL) -
     +           TX(I, WALL, LOC))/((ZLEN - ZWAF(WAFER-1)) - Z(I, LOC))
C*****INSULATE QUARTZ WALLS ENTRANCE AND EXIT FROM MIDSECTION
C            QEAST = 0.0
C*****END INSULATE QUARTZ WALLS ENTRANCE AND EXIT FROM MIDSECTION
         ELSE
            QEAST = -CONDUC(4)*(TX(I+1, WALL, LOC) -
     +              TX(I, WALL, LOC))/(Z(I+1, LOC) - Z(I, LOC))
         END IF

         IF (I .EQ. 1) THEN
            FACTE = 0.5*Z(1, LOC)/(DZENDI(LOC) + DZEND(LOC) +
     +              0.5*Z(1, LOC))
            CONDE = 1.0/((1. - FACTE)/CONDUC(6+LOC) +
     +              FACTE/CONDUC(4))
            QEDGE(WALL, LOC) = -CONDE*(TX(1, WALL, LOC) -
     +                         TENDI(LOC))/(0.5*(DZENDI(LOC) +
     +                         DZEND(LOC) + Z(1, LOC)))
            QWEST = QEDGE(WALL, LOC)
C*****ZERO END HEAT FLUXES FOR QUARTZ WALL 1 AND 2
C      QEDGE(WALL, LOC) = 0.0
C      QWEST = 0.0
C*****END ZERO END HEAT FLUXES FOR QUARTZ WALL 1 AND 2
         ELSE
            QWEST = -CONDUC(4)*(TX(I, WALL, LOC) -
     +              TX(I-1, WALL, LOC))/(Z(I, LOC) - Z(I-1, LOC))
         END IF

C*****ZERO SIDE HEAT FLUXES FOR QUARTZ WALL 1 AND 2
C      QEAST = 0.0
C      QWEST = 0.0
C*****END ZERO SIDE HEAT FLUXES FOR QUARTZ WALL 1 AND 2

C*****INSULATE THE INSIDE OF QUARTZ WALL 1
C      IF (WALL .EQ. 1) QSOUTH = 0.0
C*****END INSULATE THE INSIDE OF QUARTZ WALL 1

C*****INSULATE THE OUTSIDE OF QUARTZ WALL 1
C      IF (WALL .EQ. 1) QNORTH = 0.0
C*****END INSULATE THE OUTSIDE OF QUARTZ WALL 1

C*****INSULATE THE INSIDE OF QUARTZ WALL 2
C      IF (WALL .EQ. 2) QSOUTH = 0.0
C*****END INSULATE THE INSIDE OF QUARTZ WALL 2

C*****WRITE RADIOSITIES
C      WRITE(TEXT, *) ' '
C      WRITE(TEXT, *) 'LOC = ', LOC, ' WALL = ', WALL, ' I = ', I
C      WRITE(TEXT, *) 'J1(1) = ', J1(1, I, LOC)
C      WRITE(TEXT, *) 'J1(2) = ', J1(2, I, LOC)
C      WRITE(TEXT, *) 'J1(3) = ', J1(3, I, LOC)
C      WRITE(TEXT, *) 'J1(4) = ', J1(4, I, LOC)
C      WRITE(TEXT, *) 'JWALL(I, LOC) = ', JWALL(LOC, WALL)
C*****END WRITE RADIOSITIES

C///  COOLING GAS CONTRIBUTION

      IF (CFLOW(1)) THEN
          TFILM = (TX(I, 4+1, LOC) + TX(I, 1, LOC))/2.
          REYN = ABS(MFLOW(1))/AREAF(1)*DH(1)/MUAIR(TFILM)
          ACOND = KAIR(TFILM)
C         SUBROUTINE NUSSL(ERROR, TEXT,
C     +                    DIA, H, COND, RE, PR, NU, TURB)
          CALL NUSSL(ERROR, TEXT,
     +               DH(1), HCOF, ACOND, REYN, PRAIR, NU(1), TURB(1))

C*****FIX CONVECTIVE HEAT TRANSFER COEFFICIENT
C          HCOF = 25.0E-04
C*****END FIX CONVECTIVE HEAT TRANSFER COEFFICIENT

          CONVN(1) = HCOF*(TX(I, 4+1, LOC) -
     +               TX(I, 2, LOC))*AOUT(1)*(WZ(I+1, LOC)
     +               - WZ(I, LOC))

          CONVS(1) = HCOF*(TX(I, 1, LOC) -
     +               TX(I, 4+1, LOC))*AINN(1)*(WZ(I+1, LOC)
     +               - WZ(I, LOC))
      ELSE
          CONVN(1) = 0.0
          CONVS(1) = 0.0
      ENDIF

      IF (CFLOW(2)) THEN
          TFILM = (TX(I, 4+2, LOC) + TX(I, 2, LOC))/2.
          IF ((LOC.EQ.1) .AND. RINJEC) THEN
             REYN = ABS(MWZ(I,LOC))/AREAF(2)*DH(2)/MUAIR(TFILM)
          ELSEIF ((LOC.EQ.2) .AND. RINJEC) THEN
             REYN = ABS(MWZ(I+1,LOC))/AREAF(2)*DH(2)/MUAIR(TFILM)
          ELSE
             REYN = ABS(MFLOW(2))/AREAF(2)*DH(2)/MUAIR(TFILM)
          ENDIF
          ACOND = KAIR(TFILM)
C         SUBROUTINE NUSSL(ERROR, TEXT,
C     +                    DIA, H, COND, RE, PR, NU, TURB)
          CALL NUSSL(ERROR, TEXT,
     +               DH(2), HCOF, ACOND, REYN, PRAIR, NU(2), TURB(2))

C*****FIX CONVECTIVE HEAT TRANSFER COEFFICIENT
C          HCOF = 25.0E-04
C*****END FIX CONVECTIVE HEAT TRANSFER COEFFICIENT

C*****PRINT HEAT TRANSFER CORRELATION INFORMATION FOR OUTER FLOW
C          WRITE(TEXT, *)'--HEAT TRANS. COEF. FOR OUTER FLOW (EBAL1)--'
C          WRITE(TEXT, *)'LOC = ', LOC
C          WRITE(TEXT, *)'I = ', I
C          WRITE(TEXT, *)'CFLOW(2) = ', CFLOW(2)
C          WRITE(TEXT, *)'TURB(2) = ',TURB(2)
C          WRITE(TEXT, *)'MWZ(I, LOC) = ', MWZ(I, LOC)
C          WRITE(TEXT, *)'MWZ(I+1, LOC) = ', MWZ(I+1, LOC)
C          WRITE(TEXT, *)'REYN = ',REYN
C          WRITE(TEXT, *)'DH(2) = ',DH(2)
C          WRITE(TEXT, *)'PRAIR = ',PRAIR
C          WRITE(TEXT, *)'ACOND = ',ACOND
C          WRITE(TEXT, *)'NU(2) = ',NU(2)
C          WRITE(TEXT, *)'HCOF = ',HCOF
C*****END PRINT HEAT TRANSFER CORRELATION INFORMATION FOR OUTER FLOW

          CONVN(2) = HCOF*(TX(I, 4+2, LOC) -
     +               TX(I, 3, LOC))*AOUT(2)*(WZ(I+1, LOC)
     +               - WZ(I, LOC))
          CONVS(2) = HCOF*(TX(I, 2, LOC) -
     +               TX(I, 4+2, LOC))*AINN(2)*(WZ(I+1, LOC)
     +               - WZ(I, LOC))
      ELSE
          CONVN(2) = 0.0
          CONVS(2) = 0.0
      ENDIF

C*****FIX CONVECTIVE COOLING HEAT TRANSFER FLUXES
C
C      QFLUX = 100.0
C      IF (CFLOW(1)) THEN
C          CONVN(1) = QFLUX*AOUT(1)*(WZ(I+1, LOC) - WZ(I, LOC))
C          CONVS(1) = QFLUX*AINN(1)*(WZ(I+1, LOC) - WZ(I, LOC))
C      ELSE
C          CONVN(1) = 0.0
C          CONVS(1) = 0.0
C      ENDIF
C
C      IF (CFLOW(2)) THEN
C          CONVN(2) = QFLUX*AOUT(2)*(WZ(I+1, LOC) - WZ(I, LOC))
C          CONVS(2) = QFLUX*AINN(2)*(WZ(I+1, LOC) - WZ(I, LOC))
C      ELSE
C          CONVN(2) = 0.0
C          CONVS(2) = 0.0
C      ENDIF
C
C*****END FIX CONVECTIVE COOLING HEAT TRANSFER FLUXES

      IF (WALL .EQ. 1) QCONV = CONVS(1)
      IF (WALL .EQ. 2) QCONV = CONVS(2) - CONVN(1)

C*****PRINT QUARTZ WALL ENERGY BALANCES IN ENTRANCE AND EXIT
C      WRITE(TEXT, *) 'WALL = ', WALL, ' LOC = ', LOC
C      WRITE(TEXT, *) 'CELL = ', I
C      WRITE(TEXT, *) 'NORTH TERM = ',
C     +   AREA(WALL)*(WZ(I+1, LOC) - WZ(I, LOC))*QNORTH
C      WRITE(TEXT, *) 'SOUTH TERM = ',
C     +   AREA(WALL)*(WZ(I+1, LOC) - WZ(I, LOC))*QSOUTH
C      WRITE(TEXT, *) 'EASTERN TERM = ',
C     +   PI*(ROUTQ(WALL)**2 - RINQ(WALL)**2)*QEAST
C      WRITE(TEXT, *) 'WESTERN TERM = ',
C     +   PI*(ROUTQ(WALL)**2 - RINQ(WALL)**2)*QWEST
C      WRITE(TEXT, *) 'QCONV = ', QCONV
C      WRITE(TEXT, *) 'CONVN(1) = ',CONVN(1), ' CONVS(1) = ',CONVS(1)
C      WRITE(TEXT, *) 'CONVN(2) = ',CONVN(2), ' CONVS(2) = ',CONVS(2)
C*****END PRINT QUARTZ WALL ENERGY BALANCES IN ENTRANCE AND EXIT

         FENDS(I, WALL, LOC) =
     +      AREA(WALL)*(WZ(I+1, LOC) - WZ(I, LOC))*QNORTH -
     +      AREA(WALL)*(WZ(I+1, LOC) - WZ(I, LOC))*QSOUTH +
     +      PI*(ROUTQ(WALL)**2 - RINQ(WALL)**2)*(QEAST - QWEST) +
     +      QCONV

         IF ((PED .AND. (LOC.EQ.1)) .AND. ((WALL.EQ.1).AND.PCON)) THEN
            IF (I .EQ. ZPNTS) THEN
                FENDS(I, WALL, LOC) = FENDS(I, WALL, LOC) +
     +          PI*(PEDRO**2 - PEDRI**2)*(-QWEST)
            ELSE
                FENDS(I, WALL, LOC) = FENDS(I, WALL, LOC) +
     +          PI*(PEDRO**2 - PEDRI**2)*(QEAST - QWEST)
            ENDIF
         ELSE
         ENDIF

C*****SET CONSTANT TEMPERATURE FOR QUARTZ WALL 1
C      IF (WALL .EQ. 1) FENDS(I, WALL, LOC) = TX(I, WALL, LOC) - 800.0
C*****END SET CONSTANT TEMPERATURE FOR QUARTZ WALL 1

C*****SET CONSTANT TEMPERATURE FOR QUARTZ WALL 2
C      IF (WALL .EQ. 2) FENDS(I, WALL, LOC) = TX(I, WALL, LOC) - 800.0
C*****END SET CONSTANT TEMPERATURE FOR QUARTZ WALL 2

C        DIVIDED VOL*RHO*CP FOR THE ELEMENT

         RCPV = PI*(ROUTQ(WALL)**2 -
     +   RINQ(WALL)**2)*(WZ(I+1, LOC) - WZ(I, LOC))*RHO(4)*CP(4)

         IF ((PED .AND. (LOC.EQ.1)) .AND. (WALL.EQ.1)) THEN
             RCPV = RCPV + PI*(PEDRO**2 - PEDRI**2)*(WZ(I+1, LOC) -
     +       WZ(I, LOC))*RHO(4)*CP(4)
         ELSE
         ENDIF

         FENDS(I, WALL, LOC) = FENDS(I, WALL, LOC)/RCPV

         STEADY = FENDS(I, WALL, LOC)

         IF (TIME) THEN
            FENDS(I, WALL, LOC) = FENDS(I, WALL, LOC) +
     +                            TX0(I, WALL, LOC)
         END IF

C*****FIX TEMPERATURES OF ENDS OF QUARTZ TUBES AT CAP TEMPERATURES
C      IF (I .EQ. 1) THEN
C         FENDS(1, WALL, LOC) = TX(I, WALL, LOC) - TEND(LOC)
C      ELSE
C      ENDIF
C*****END FIX TEMPERATURES OF ENDS OF QUARTZ TUBES AT CAP TEMPERATURES

C*****PRINT STEADY-STATE AND TRANSIENT RESIDUAL
C      WRITE(TEXT, *) ' '
C      WRITE(TEXT, *) 'TIME = ', TIME
C      WRITE(TEXT, *) 'LOC = ', LOC, ' WALL = ', WALL, ' I = ', I
C      WRITE(TEXT, *) 'STEADY RESID. = ', STEADY
C      WRITE(TEXT, *) 'TRANSIENT RESID. = ', FENDS(I, WALL, LOC)
C      WRITE(TEXT, *) 'TRANSIENT TERM = ',
C     +   PI*(ROUTQ(WALL)**2 - RINQ(WALL)**2)*
C     +   (WZ(I+1, LOC) - WZ(I, LOC))*RHO(4)*CP(4)*TX0(I, WALL, LOC)
C*****END PRINT STEADY-STATE AND TRANSIENT RESIDUAL

1050  CONTINUE

C///  ENERGY BALANCE FOR THE WALL HEATER ELEMENTS

      QNORTH = -COND*(TX(I, 4, LOC) - TX(I, 3, LOC))/(RIN - RTUBE)
      QSOUTH = J1(2, I, LOC)*FH2 + J1(1, I, LOC)*FHH - J1(1, I, LOC)

C*****INSULATE BOTTOM OF HEATERS
C      QSOUTH = 0.0
C*****END INSULATE BOTTOM OF HEATERS

      IF (I .EQ. 1) THEN
         FACTE = 0.5*Z(1, LOC)/(DZENDI(LOC) + DZEND(LOC) +
     +           0.5*Z(1, LOC))
         CONDE = 1.0/((1. - FACTE)/CONDUC(6+LOC) + FACTE/CONDUC(2))
         QEDGE(3, LOC) = -CONDE*(TX(1, 3, LOC) -
     +                   TENDI(LOC))/(0.5*(DZENDI(LOC) +
     +                   DZEND(LOC) + Z(1, LOC)))
         QWEST = QEDGE(3, LOC)
C*****ZERO EDGE HEAT FLUX FOR HEATER CYLINDER
C      QEDGE(3, LOC) = 0.0
C      QWEST = 0.0
C*****END ZERO EDGE HEAT FLUX FOR HEATER CYLINDER
      ELSE
         QWEST = -CONDUC(2)*(TX(I, 3, LOC) -
     +            TX(I-1, 3, LOC))/(Z(I, LOC) - Z(I-1, LOC))
      END IF

      IF (I .EQ. ZPNTS) THEN
         IF (LOC. EQ. 1) QEAST = -CONDUC(2)*(TEX(WAFER, 3) -
     +           TX(I, 3, LOC))/(ZWAF(WAFER) - Z(I, LOC))
         IF (LOC. EQ. 2) QEAST = -CONDUC(2)*(TEX(WAFER, 3) -
     +           TX(I, 3, LOC))/((ZLEN - ZWAF(WAFER-1)) - Z(I, LOC))
C*****INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
C         QEAST = 0.0
C*****END INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
      ELSE
         QEAST = -CONDUC(2)*(TX(I+1, 3, LOC) -
     +           TX(I, 3, LOC))/(Z(I+1, LOC) - Z(I, LOC))
      END IF

C*****ZERO SIDE HEAT FLUX IN HEATERS
C      QEAST = 0.0
C      QWEST = 0.0
C*****END ZERO SIDE HEAT FLUX IN HEATERS

C*****PRINT HEATER ENERGY BALANCES
C      WRITE(TEXT, *) 'HEATER ENERGY BALANCES'
C      WRITE(TEXT, 566) LOC, I,
C     +            2.*PI*RTUBEO*(WZ(I+1, LOC) - WZ(I, LOC))*QNORTH,
C     +            2.*PI*RTUBEI*(WZ(I+1, LOC) - WZ(I, LOC))*QSOUTH,
C     +            PI*(RTUBEO**2 - RTUBEI**2)*QEAST,
C     +            PI*(RTUBEO**2 - RTUBEI**2)*QWEST,
C     +            PWALL(I, LOC)
C      WRITE(TEXT, *) 'CONVN(2) = ', CONVN(2)
C      WRITE(TEXT, 567) QNORTH, COND,
C     +                 WZ(I+1, LOC), WZ(I, LOC), RTUBEO
C
C  566 FORMAT('WALL',1X,'LOC =',I2,' I = ',I2,5(E11.4))
C  567 FORMAT(2X,'QNORTH = ',E11.4,1X,'COND = ',E11.4/
C     +       2X,'WZ(I+1, LOC) = ',E11.4,1X,'WZ(I, LOC) = ',E11.4,
C     +       1X,'RTUBEO = ',E11.4)
C*****END PRINT HEATER ENERGY BALANCES

      FENDS(I, 3, LOC) =
     +   2.*PI*RTUBEO*(WZ(I+1, LOC) - WZ(I, LOC))*QNORTH -
     +   2.*PI*RTUBEI*(WZ(I+1, LOC) - WZ(I, LOC))*QSOUTH +
     +   PI*(RTUBEO**2 - RTUBEI**2)*(QEAST - QWEST) -
     +   PWALL(I, LOC) - CONVN(2)

C     DIVIDED VOL*RHO*CP FOR THE ELEMENT

      FENDS(I, 3, LOC) = FENDS(I, 3, LOC)/(PI * (RTUBEO**2
     +   - RTUBEI**2) * (WZ(I + 1, LOC) - WZ(I, LOC))*RHO(2)*CP(2))

C*****SET CONSTANT TEMPERATURE FOR HEATERS
C      IF (LOC .EQ. 1) THEN
C         FENDS(I, 3, LOC) = TX(I, 3, LOC) - 800.0
C      ELSE
C         FENDS(I, 3, LOC) = TX(I, 3, LOC) - 800.0
C      END IF
C*****END SET CONSTANT TEMPERATURE FOR HEATERS

C///  TRANSIENT TERM.

      IF (TIME) THEN
         FENDS(I, 3, LOC) = FENDS(I, 3, LOC) + TX0(I, 3, LOC)
      END IF

C///  ENERGY BALANCE FOR THE INSULATION ELEMENT.

      QSOUTH = QNORTH
      QNORTH = HCOEF(LOC+3)*(TX(I, 4, LOC) - TAMB) +
     +         EMISSI(LOC+3)*SIG*(TX4(I, 4, LOC) - TAMB4)

      IF (I .EQ. 1) THEN
         FACTE = 0.5*Z(1, LOC)/(DZENDI(LOC) + DZEND(LOC) +
     +           0.5*Z(1, LOC))
         CONDE = 1.0/((1. - FACTE)/CONDUC(6+LOC) + FACTE/CONDUC(3))
         QEDGE(4, LOC) = -CONDE*(TX(1, 4, LOC) -
     +                   TENDI(LOC))/(0.5*(DZENDI(LOC) +
     +                   DZEND(LOC) + Z(1, LOC)))
         QWEST = QEDGE(4, LOC)
C*****ZERO EDGE HEAT FLUX FOR HEATER INSULATION BLANKET
C      QWEST = 0.0
C      QEDGE(4, LOC) = 0.0
C*****END ZERO EDGE HEAT FLUX FOR HEATER INSULATION BLANKET
      ELSE
         QWEST = -CONDUC(3)*(TX(I, 4, LOC) -
     +            TX(I-1, 4, LOC))/(Z(I, LOC) - Z(I-1, LOC))
      END IF

      IF (I .EQ. ZPNTS) THEN
         IF (LOC. EQ. 1) QEAST = -CONDUC(3)*(TEX(WAFER, 4) -
     +           TX(I, 4, LOC))/(ZWAF(WAFER) - Z(I, LOC))
         IF (LOC. EQ. 2) QEAST = -CONDUC(3)*(TEX(WAFER, 4) -
     +           TX(I, 4, LOC))/((ZLEN - ZWAF(WAFER-1)) - Z(I, LOC))
C*****INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
C         QEAST = 0.0
C*****END INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
      ELSE
         QEAST = -CONDUC(3)*(TX(I+1, 4, LOC) -
     +            TX(I, 4, LOC))/(Z(I+1, LOC) - Z(I, LOC))
      END IF

C*****ZERO SIDE HEAT FLUXES IN INSULATION
C      QEAST = 0.0
C      QWEST = 0.0
C*****END ZERO SIDE HEAT FLUXES IN INSULATION

C*****PRINT INSULATION ENERGY BALANCES
C      WRITE(TEXT, *) 'INSULATION ENERGY BALANCES'
C      WRITE(TEXT, 577) LOC, I,
C     +            2.*PI*RINO*(WZ(I+1, LOC) - WZ(I, LOC))*QNORTH,
C     +            2.*PI*RINI*(WZ(I+1, LOC) - WZ(I, LOC))*QSOUTH,
C     +            PI*(RINO**2 - RINI**2)*QEAST,
C     +            PI*(RINO**2 - RINI**2)*QWEST
C  577 FORMAT('INSU',1X,'LOC =',I2,' I = ',I2,4(E11.4))
C*****END PRINT INSULATION ENERGY BALANCES

      FENDS(I, 4, LOC) =
     +   2.*PI*RINO*(WZ(I+1, LOC) - WZ(I, LOC))*QNORTH -
     +   2.*PI*RINI*(WZ(I+1, LOC) - WZ(I, LOC))*QSOUTH +
     +   PI*(RINO**2 - RINI**2)*(QEAST - QWEST)

C     DIVIDED VOL*RHO*CP FOR THE ELEMENT

      FENDS(I, 4, LOC) =
     +   FENDS(I, 4, LOC)/(PI*(RINO**2 - RINI**2)*(WZ(I+1, LOC) -
     +   WZ(I, LOC))*RHO(3)*CP(3))

      LOS2(I, LOC) = 2.*PI*RINO*(WZ(I+1, LOC) - WZ(I, LOC))*QNORTH

C*****SET CONSTANT INSULATION TEMPERATURE
C      IF (LOC .EQ. 1) THEN
C         FENDS(I, 4, LOC) = TX(I, 4, LOC) - 800.0
C      ELSE
C         FENDS(I, 4, LOC) = TX(I, 4, LOC) - 800.0
C      END IF
C*****END SET CONSTANT INSULATION TEMPERATURE

C///  TRANSIENT TERM.

      IF (TIME) THEN
         FENDS(I, 4, LOC) =
     +      FENDS(I, 4, LOC) + TX0(I, 4, LOC)
      END IF

      IF (NUMEX .EQ. 4) GO TO 1060

C-----------------------AIR ENERGY BALANCE----------------------

C///  ENERGY BALANCE FOR AIR BETWEEN HEATERS AND OUTER QUARTZ JAR
C     AND AIR BETWEEN TWO QUARTZ JARS.

C*****PRINT GAS ENERGY BALANCE TERMS
C
C      WRITE(TEXT, *) ' '
C      WRITE(TEXT, *) ' GAS COOLING ENERGY BALANCE TERMS'
C
C*****END PRINT GAS ENERGY BALANCE TERMS

      DO 1055 FLOW = 1, 2

         IF (LOC .EQ. 1) THEN
            IF (I .EQ. 1) THEN
                IF(MFLOW(FLOW).GE. 0.0) TWEST1 = TTGAS(FLOW)
                TEAST1 = TX(I+1, 4+FLOW, LOC)
            ELSEIF (I .EQ. ZPNTS) THEN
                TWEST1 = TX(I-1, 4+FLOW, LOC)
                TEAST1 = TEX(1, 4+FLOW)
            ELSE
                TWEST1 = TX(I-1, 4+FLOW, LOC)
                TEAST1 = TX(I+1, 4+FLOW, LOC)
            ENDIF
         ELSE
         ENDIF

         IF (LOC .EQ. 2) THEN
            IF (I .EQ. 1) THEN
                TWEST1 = TX(I+1, 4+FLOW, LOC)
                IF(MFLOW(FLOW) .LT. 0.0) TEAST1 = TTGAS(FLOW)
            ELSEIF (I .EQ. ZPNTS) THEN
                TWEST1 = TEX(WAFER, 4+FLOW)
                TEAST1 = TX(I-1, 4+FLOW, LOC)
            ELSE
                TWEST1 = TX(I+1, 4+FLOW, LOC)
                TEAST1 = TX(I-1, 4+FLOW, LOC)
            ENDIF
         ELSE
         ENDIF

         IF (FLOW .EQ. 1) THEN
            MFLOWW = MFLOW(1)
            MFLOWE = MFLOW(1)
            EINJ = 0.0
         ELSEIF ((FLOW .EQ. 2).AND.(LOC .EQ. 1)) THEN
            IF (RINJEC) THEN
               MFLOWW = MWZ(I, LOC)
               MFLOWE = MWZ(I+1, LOC)
               EINJ = MINJE(I, LOC)*CPAIR*TTGAS(2)
            ELSE
               MFLOWW = MFLOW(2)
               MFLOWE = MFLOW(2)
               EINJ = 0.0
            ENDIF
         ELSEIF ((FLOW .EQ. 2).AND.(LOC .EQ. 2)) THEN
            IF (RINJEC) THEN
               MFLOWW = MWZ(I+1, LOC)
               MFLOWE = MWZ(I, LOC)
               EINJ = MINJE(I, LOC)*CPAIR*TTGAS(2)
            ELSE
               MFLOWW = MFLOW(2)
               MFLOWE = MFLOW(2)
               EINJ = 0.0
            ENDIF
         ELSE
         ENDIF

         EEAST = TX(I, 4+FLOW, LOC)*MAX(MFLOWE, 0.0)*CPAIR -
     +           TEAST1*MAX(-MFLOWE, 0.0)*CPAIR
         EWEST = TWEST1*MAX(MFLOWW, 0.0)*CPAIR -
     +           TX(I, 4+FLOW, LOC)*MAX(-MFLOWW, 0.0)*CPAIR

         FENDS(I, 4+FLOW, LOC) =
     +      CONVN(FLOW) - CONVS(FLOW) + EEAST - EWEST - EINJ

C*****PRINT GAS ENERGY BALANCE TERMS
C
C         WRITE(TEXT, *) ' FLOW = ', FLOW, ' LOC = ', LOC,' I = ',I
C         WRITE(TEXT, *) ' MFLOWW = ', MFLOWW, ' MFLOWE = ', MFLOWE
C         WRITE(TEXT, *) ' MINJE = ', MINJE(I, LOC)
C         WRITE(TEXT, *) ' EINJ = ', EINJ
C         WRITE(TEXT, *) ' EWEST = ',EWEST,' EEAST = ', EEAST
C         WRITE(TEXT, *) ' CONVN - CONVS + EEAST - EWEST - EINJ = ',
C     +                    CONVN(FLOW)-CONVS(FLOW)+EEAST-EWEST-EINJ
C         WRITE(TEXT, *) ' CONVN = ',CONVN(FLOW),' CONVS = ',CONVS(FLOW)
C         WRITE(TEXT, *) ' CONVN - CONVS = ', CONVN(FLOW)-CONVS(FLOW)
C         WRITE(TEXT, *) ' TEAST1 = ',TEAST1, ' TWEST1 = ', TWEST1
C         WRITE(TEXT, *) ' MFLOW = ',MFLOW(FLOW),' CPAIR = ',CPAIR
C         WRITE(TEXT, *) ' TGAS1 = ',TX(I,4+1,LOC),
C     +                  ' TGAS2 = ',TX(I,4+2,LOC)
C         WRITE(TEXT, *) ' TQUAR1 = ',TX(I,1,LOC),
C     +                  ' TQUAR2 = ',TX(I,2,LOC)
C         WRITE(TEXT, *) ' THEAT = ',TX(I,3,LOC)
C
C*****END PRINT GAS ENERGY BALANCE TERMS

C*****INCLUDE ONLY ADVECTIVE TERMS IN GAS ENERGY BALANCE
C
C         FENDS(I, 4+FLOW, LOC) = EEAST - EWEST
C
C*****END INCLUDE ONLY ADVECTIVE TERMS IN GAS ENERGY BALANCE

C*****CONVECTION IN ENTRANCE
C
C      IF (LOC.EQ.1)
C     +   FENDS(I, 4+FLOW, LOC) =  EEAST - EWEST
C     +                            + CONVN(FLOW) - CONVS(FLOW)
C
C*****END CONVECTION IN ENTRANCE

C*****CONVECTION IN EXIT
C
C      IF (LOC.EQ.2)
C     +   FENDS(I, 4+FLOW, LOC) =  EEAST - EWEST
C     +                            + CONVN(FLOW) - CONVS(FLOW)
C
C*****END CONVECTION IN EXIT

C*****FIX GAS TEMPERATURES IN INLET
C
C         IF (LOC .EQ. 1) FENDS(I, 4+FLOW, LOC) =
C     +                   TX(I, 4+FLOW, LOC) - 400.0
C
C*****END FIX GAS TEMPERATURES IN INLET

C*****FIX GAS TEMPERATURES IN OUTLET
C
C         IF (LOC .EQ. 2)  FENDS(I, 4+FLOW, LOC) =
C     +                    TX(I, 4+FLOW, LOC) - 500.0
C
C*****END FIX GAS TEMPERATURES IN OUTLET

         TMASS = (VOL(FLOW)*(WZ(I+1, LOC) -
     +   WZ(I, LOC))*RHOAIR(TX(I, 4+FLOW, LOC))*CPAIR)

C        DIVIDED VOL*RHO*CP FOR THE ELEMENT.

         FENDS(I, 4+FLOW, LOC) =
     +   FENDS(I, 4+FLOW, LOC)/(VOL(FLOW)*(WZ(I+1, LOC) -
     +   WZ(I, LOC))*RHOAIR(TX(I, 4+FLOW, LOC))*CPAIR)

C*****SET THE GAS TEMPERATURE DERIVATIVE TO A CONSTANT.
C
C         FENDS(I, 4+FLOW, LOC) = 0.0
C
C*****END SET THE GAS TEMPERATURE DERIVATIVE TO A CONSTANT.

C///  TRANSIENT TERM.

      IF (TIME) THEN
         FENDS(I, 4+FLOW, LOC) =
     +      FENDS(I, 4+FLOW, LOC) + TX0(I, 4+FLOW, LOC)
      END IF

1055  CONTINUE

C-------------------------------------------------------------

C///  BOTTOM OF THE LOOP OVER THE POINTS IN THE ENTRANCE OR EXIT.

1060  CONTINUE

C///  RADIOSITY FROM REACTOR END WALL.

      SUMSE = 0.0
      DO 1070 I = 1, ZPNTS
         SUMSE = SUMSE + JWALL(I, LOC)*FESE(I, LOC)
1070  CONTINUE

C*****PRINT RADIOSITY INFORMATION IN ENTRANCE AND EXIT
C      WRITE(TEXT, *) '-----------------------------------'
C*****END PRINT RADIOSITY INFORMATION IN ENTRANCE AND EXIT
      SUMWE = 0.0
      DO 1080 J = 1, RPNTS
         SUMWE = SUMWE + EMISS*FEW(J, LOC)*TWAF4(WAFER, J)
C*****PRINT RADIOSITY INFORMATION IN ENTRANCE AND EXIT
C      WRITE(TEXT, *) ''
C      WRITE(TEXT, *) ' LOC = ', LOC, ' WAFER = ', WAFER, ' J = ', J
C      WRITE(TEXT, *) ' FEW(J, LOC) = ', FEW(J, LOC)
C      WRITE(TEXT, *) ' TWAF4(WAFER, J) = ', TWAF4(WAFER, J)
C*****END PRINT RADIOSITY INFORMATION IN ENTRANCE AND EXIT
1080  CONTINUE

C///  RADIOSITY EQUATION

C*****PRINT RADIOSITY INFORMATION IN ENTRANCE AND EXIT
C      WRITE(TEXT, *) ' '
C      WRITE(TEXT, *) ' LOC = ', LOC, ' ZPNTS = ', ZPNTS
C      WRITE(TEXT, *) ' EMISSE = ',EMISSE, ' EMISS = ', EMISS
C      WRITE(TEXT, *) ' SUMWE = ',SUMWE, ' SUMSE = ', SUMSE
C      WRITE(TEXT, *) ' EGAP(LOC) = ', EGAP(LOC),
C     +               ' JEND(LOC) = ', JEND(LOC)
C      WRITE(TEXT, *) ' FEG(LOC) = ', FEG(LOC),
C     +               ' FEE(LOC) = ', FEG(LOC)
C      WRITE(TEXT, *) ' TEND4(LOC) = ',TEND4(LOC), ' SIG = ', SIG
C*****END PRINT RADIOSITY INFORMATION IN ENTRANCE AND EXIT

C      FENDS(ZPNTS + 1, 1, LOC) =
C     +   JEND(LOC) - (EMISSE*SIG*TEND4(LOC) +
C     +   (1. - EMISSE)*(SUMSE + SIG*SUMWE +
C     +   EMISS*EGAP(LOC)*FEG(LOC)))/(1. - (1. - EMISSE)*FEE(LOC))

C///  ENERGY BALANCE ON REACTOR END WALL ELEMENTS.

      GEND = JEND(LOC) * FEE(LOC) + SUMSE + SIG * SUMWE
     +   + EMISS * EGAP(LOC) * FEG(LOC)

      FACTE = DZEND(LOC)/2./((DZEND(LOC)+ DZENDI(LOC))/2.)
      CONDE = 1.0/((1. - FACTE)/CONDUC(6+LOC) + FACTE/CONDUC(4+LOC))

      QEAST = JEND(LOC) - GEND
      QWEST = -CONDE*(TEND(LOC) - TENDI(LOC))/((DZEND(LOC) +
     +        DZENDI(LOC))/2.)

C*****ZERO HEAT LOSS FROM REACTOR ENDS
C      QWEST = 0.0
C*****END ZERO HEAT LOSS FROM REACTOR ENDS

C*****WRITE ENERGY BALANCE FOR REACTOR ENDS
C
C      WRITE(TEXT, *) 'END CAP ENERGY BALANCE'
C      WRITE(TEXT, *) 'LOC = ',LOC
C      WRITE(TEXT, *) 'GEND = ', GEND, ' JEND = ', JEND(LOC)
C      WRITE(TEXT, *) 'QEAST = ',QEAST, ' QWEST = ', QWEST
C      WRITE(TEXT, *) 'EASTERN TERM = ', PI*RINQ(1)**2*QEAST
C      WRITE(TEXT, *) 'WESTERN TERM = ', PI*RINQ(1)**2*QWEST
C      WRITE(TEXT, *) 'SOURCE TERM = ', PCAP(LOC)
C
CC      WRITE(TEXT, *) 'EASTERN TERM = ', PI*RTUBEI**2*QEAST
CC      WRITE(TEXT, *) 'WESTERN TERM = ', PI*RTUBEI**2*QWEST
C
C*****END WRITE ENERGY BALANCE FOR REACTOR ENDS

C      FENDS(ZPNTS + 1, 1, LOC) = PI*RTUBEI**2*(QEAST - QWEST) -
C                                 PCAP(LOC)

      FENDS(ZPNTS + 1, 1, LOC) = PI*RINQ(1)**2*(QEAST - QWEST) -
     +                           PCAP(LOC)

C     DIVIDED VOL*RHO*CP FOR THE ELEMENT

      FENDS(ZPNTS + 1, 1, LOC) =
     +   FENDS(ZPNTS+1, 1, LOC)/
     +  (PI*RINQ(1)**2*DZEND(LOC)*RHO(4+LOC)*CP(4+LOC))

C*****SET CONSTANT END TEMPERATURE
C      FENDS(ZPNTS + 1, 1, LOC) = TEND(LOC) - 800.0
C*****END SET CONSTANT END TEMPERATURE

C///  TRANSIENT TERM.

      IF (TIME) THEN
         FENDS(ZPNTS + 1, 1, LOC) =
     +      FENDS(ZPNTS+1, 1, LOC) + TEND0(LOC)
      END IF

C///  SET CONSTANT END TEMPERATURE

      IF (TCAPS(LOC)) FENDS(ZPNTS + 1, 1, LOC) = TEND(LOC) - TCAP(LOC)

C///  ENERGY BALANCE ON REACTOR END WALL INSULATION ELEMENTS.

      QEAST = QWEST
      QWEST = HCOEF(LOC)*(TAMB - TENDI(LOC)) +
     +        EMISSI(LOC)*SIG*(TAMB4 - TENDI4(LOC))

C      FENDS(ZPNTS + 1, 2, LOC) = PI*RTUBEI**2*(QEAST - QWEST)
C      FENDS(ZPNTS + 1, 2, LOC) = PI*RINQ(1)**2*(QEAST - QWEST)

      FENDS(ZPNTS + 1, 2, LOC) = PI*RINQ(1)**2*QEAST -
     +   PI*RINO**2*QWEST -
     +   (DZEND(LOC) + DZENDI(LOC))*2.0*PI*RINO*QWEST +
     +   PI*(ROUTQ(1)**2 - RINQ(1)**2)*QEDGE(1, LOC) +
     +   PI*(ROUTQ(2)**2 - RINQ(2)**2)*QEDGE(2, LOC) +
     +   PI*(RTUBEO**2 - RTUBEI**2)*QEDGE(3, LOC) +
     +   PI*(RINO**2 - RINI**2)*QEDGE(4, LOC)

      IF (PCON .AND. (LOC.EQ.1)) FENDS(ZPNTS + 1, 2, LOC) =
     +   FENDS(ZPNTS + 1, 2, LOC) +
     +   PI*(PEDRO**2 - PEDRI**2)*QEDGE(1, LOC)

C     DIVIDED VOL*RHO*CP FOR THE ELEMENT

      FENDS(ZPNTS + 1, 2, LOC) =
     +   FENDS(ZPNTS + 1, 2, LOC)
     +   /(PI*RINO**2*DZENDI(LOC)*RHO(6+LOC)*CP(6+LOC))

      LOS3(LOC) = PI*RINO**2*QWEST +
     +           (DZEND(LOC) + DZENDI(LOC))*2.0*PI*RINO*QWEST

C*****SET CONSTANT INSULATION TEMPERATURE
C      FENDS(ZPNTS + 1, 2, LOC) = TENDI(LOC) - 800.0
C*****END SET CONSTANT INSULATION TEMPERATURE

C///  TRANSIENT TERM.

      IF (TIME) THEN
         FENDS(ZPNTS + 1, 2, LOC) =
     +      FENDS(ZPNTS + 1, 2, LOC) + TENDI0(LOC)
      END IF

C///  HEAT FLUX ON THE WEST FACE OF THE FIRST WAFER.

      DO 1100 J = 1, RPNTS
         SUMSE = 0.0
         DO 1090 I = 1, ZPNTS
            SUMSE = SUMSE + FWSE(J, I, LOC)*JWALL(I, LOC)
1090     CONTINUE
         QWAF(J, LOC) = EMISS*(JEND(LOC)*FWE(J, LOC) + SUMSE)
     +                  - EMISS*SIG*TWAF4(WAFER, J)
C*****PRINT HEAT FLUX ON END WAFERS (ENTRANCE AND EXIT CALC.)
C         WRITE(TEXT, *) 'LOC = ', LOC, ' CELL = ', J
C         WRITE(TEXT, *) 'END WAFER HEAT FLUX = ',
C     +                   (EMISS*(JEND(LOC)*FWE(J, LOC) + SUMSE)
C     +                  - EMISS*SIG*TWAF4(WAFER, J))
C*****END PRINT HEAT FLUX ON END WAFERS (ENTRANCE AND EXIT CALC.)
1100  CONTINUE

C///  SET FIXED TEMPERATURES IN INLET AND EXIT SECTIONS.

      IF (TFIX) THEN
         DO 1110 I = IBEG(LOC), IEND(LOC)
            FENDS(I, 3, LOC) = TX(I, 3, LOC) - TSIDE(I, LOC)
1110     CONTINUE
      END IF

C///  BOTTOM OF THE LOOP OVER THE ENTRANCE OR EXIT SECTION.

1120  CONTINUE


C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID

99103 FORMAT
     +   (/1X, A, 'ERROR. ERROR IN CALCULATION OF THERMOCOUPLE'
     +   /10X, 'SHAPEFACTORS.')

C///  EXIT.

99999 CONTINUE

      RETURN
      END
      SUBROUTINE INTERP (NPTS, X, XX, TT, T)
C
C*****DOUBLE PRECISION
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END DOUBLE PRECISION
C
C*****SINGLE PRECISION
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END SINGLE PRECISION
C
      DIMENSION XX(NPTS), TT(NPTS)
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
      END
      SUBROUTINE MASSD
     +  (ERROR, TEXT,
     +   MAXIN, MFLOW, MINJ, MINJE, MWZ, MWZW, NUMIN, NUMWF,
     +   PDIM, PEX, PIN, VINJ, WZ,
     +   WZWAF, ZBEGI, ZENDI, ZLEN)

C///////////////////////////////////////////////////////////////////////
C
C     MASSD
C
C///  COMPUTES THE MASS FLOW DISTRIBUTION FOR COOLING AIR INJECTION
C     BETWEEN QUARTZ WALLS AND OUTER QUARTZ WALL AND HEATERS.
C
C     FOR ANY INJECTOR "I", THE MASS INJECTION VINJ(I) IS ASSUMED
C     TO TAKE PLACE UNIFORMLY BETWEEN ZBEGI(I) AND ZENDI(I) (UNITS
C     ARE CFM CUBIC FEET/MIN FOR VINJ(I) AND CM FOR ZBEGI AND ZENDI).
C
C     THE ROUTINE RETURNS MINJE(I, LOC) AND MINJ(I) INJECTED IN EACH
C     CELL IN THE ENTRANCE AND EXIT REGION AND THE REGION OVER THE
C     THE WAFERS.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO
C             AND NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     MAXIN   INPUT INTEGER - MAXIMUM NUMBER OF INJECTORS.
C
C     MFLOW   INPUT REAL DIMENSIONED (2) - FLOWRATES IN COOLING
C             CHANNELS (G/SEC).
C
C     MWZ     OUTPUT REAL DIMENSIONED (PDIM+1, 2) - FLOWRATE AT
C             CELL WALLS IN THE ENTRANCE AND EXIT REGION.
C
C     MWZW    OUTPUT REAL DIMENSIONED (NUMWF+1) - FLOWRATE AT
C             CELL WALLS OVER THE WAFER REGION.
C
C     NUMIN   INPUT INTEGER - NUMBER OF INJECTORS.
C
C     PEX     INPUT INTEGER - NUMBER OF POINTS IN THE EXIT.
C
C     PIN     INPUT INTEGER - NUMBER OF POINTS IN THE INLET.
C
C     PDIM    INPUT INTEGER - MAXIMUM DIMENSION OF SOME ARRAYS.
C
C     MINJ    OUTPUT REAL DIMENSIONED (NUMWF) - WALL MASS FLOW AT GRID
C             POINTS ABOVE THE WAFERS (G/SEC).
C
C     VINJ    INPUT REAL DIMENSIONED (MAXIN) - VINJ(I) IS THE TOTAL
C             VOLUMETRIC FLOWRATE (CFM) FROM INJECTOR "I".
C
C     MINJE   OUTPUT REAL DIMENSIONED (PDIM, 2) - MINJE(I, 1) IS THE
C             INJECTED MASS FLOW IN THE WALL REGION IN THE ENTRANCE.
C             MINJE(I, 2) IS THE MASS FLOW IN THE EXIT REGION (G/SEC).
C
C     WZ      INPUT REAL DIMENSION (PDIM+1, 2) - LOCATION OF CELL
C             WALLS.
C
C     WZWAF   INPUT REAL DIMENSION (NUMWF+1) - LOCATION OF CELL
C             WALLS.
C
C     ZBEGI   INPUT REAL DIMENSIONED (MAXIN) - ZBEGI(I) IS THE
C             BEGINNING Z COORDINATE FOR INJECTOR "I".
C
C     ZENDI   INPUT REAL DIMENSIONED (MAXIN) - ZENDI IS THE END Z
C             COORDINATE FOR INJECTOR "I".
C
C     ZLEN    INPUT REAL - TOTAL LENGTH OF THE REACTOR (CM).
C
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   CELL, CELLS, I, INNER, K, LOC, MAXIN, NUMIN,
     +   NUMWF, OUTER,
     +   PDIM, PEX, PIN, QINNER, QOUTER, TEXT, WALLS

      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   FRAC, PART1, PART2, PI, SUM1, SUM2, MFLOW, MINJ, MINJE, MWZ,
     +   MWZW, VINJ, WZ, WZINN, WZOUT,
     +   WZWAF, ZBEGI, ZENDI, ZLEN

      PARAMETER (ID = 'POW:  ')
      PARAMETER (QINNER = 0, QOUTER = 1)
      DIMENSION
     +   MFLOW(2), MINJ(NUMWF), MINJE(PDIM, 2), MWZ(PDIM+1, 2),
     +   MWZW(NUMWF+1), VINJ(MAXIN),
     +   WZ(PDIM+1, 2), WZWAF(NUMWF+1), ZBEGI(MAXIN), ZENDI(MAXIN)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE
C
C///////////////////////////////////////////////////////////////////////

      DATA PI /3.1415926/

C///  ZERO MASS INJECTED FOR ALL CELLS.

      DO 1080 LOC = 1, 2
         DO 1070 K = 1, PDIM
            MINJE(K, LOC) = 0.0
 1070    CONTINUE
 1080 CONTINUE
      DO 1090 K = 1, NUMWF
         MINJ(K) = 0.0
 1090 CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (2) FIND THE MASS INJECTED TO EACH CELL IN THE ENTRANCE AND EXIT.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP FOR THE ENTRANCE OR EXIT.

      DO 3095 LOC = 1, 2

C///  TOP OF THE LOOP FOR CELLS.

      IF (LOC .EQ. 1) THEN
         CELLS = PIN
         WALLS = CELLS + 1
      ELSE
         CELLS = PEX
         WALLS = CELLS + 1
      ENDIF

      DO 3090 CELL = 1, CELLS
         INNER = CELL + QINNER
         OUTER = CELL + QOUTER

C///  TOP OF THE LOOP FOR THE INJECTORS.

      DO 3070 I = 1, NUMIN

         IF (LOC .EQ. 1) THEN
             WZINN = WZ(INNER, LOC)
             WZOUT = WZ(OUTER, LOC)
         ELSE
C            COORDINATES ARE REVERSED IN THE EXIT
             WZOUT = ZLEN - WZ(INNER, LOC)
             WZINN = ZLEN - WZ(OUTER, LOC)
         ENDIF

         IF ((WZOUT .GE. ZBEGI(I)) .AND. (WZINN .LE. ZENDI(I))) THEN
             PART1 = WZINN - ZBEGI(I)
             PART2 = ZENDI(I) - WZOUT
             IF (PART1 .LE. 0.0) PART1 = 0.0
             IF (PART2 .LE. 0.0) PART2 = 0.0
             FRAC = (ZENDI(I) - ZBEGI(I) - PART1 - PART2)
     +              /(ZENDI(I) - ZBEGI(I))
         ELSE
             FRAC = 0.0
         ENDIF

C///  ACCUMULATE MASS FOR EACH CELL AND CONVERT FROM
C     (FT**3/MIN TO G/SEC).

         MINJE(CELL, LOC) = FRAC*472.01/773.2*VINJ(I) + MINJE(CELL, LOC)

C///  BOTTOM OF THE LOOP FOR INJECTORS.

 3070 CONTINUE

C///  BOTTOM OF THE LOOP FOR CELLS.

 3090 CONTINUE

C///  BOTTOM OF THE LOOP OVER ENTRANCE AND EXIT.

 3095 CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (3) FIND THE MASS INJECTED TO EACH CELL IN THE WAFER REGION.
C
C///////////////////////////////////////////////////////////////////////

      CELLS = NUMWF - 1

C///  TOP OF THE LOOP FOR CELLS.

      DO 4090 CELL = 1, CELLS
         INNER = CELL + QINNER
         OUTER = CELL + QOUTER

C///  TOP OF THE LOOP FOR HEATERS.

      DO 4070 I = 1, NUMIN

         WZINN = WZWAF(INNER)
         WZOUT = WZWAF(OUTER)

         IF ((WZOUT .GE. ZBEGI(I)) .AND. (WZINN .LE. ZENDI(I))) THEN
             PART1 = WZINN - ZBEGI(I)
             PART2 = ZENDI(I) - WZOUT
             IF (PART1 .LE. 0.0) PART1 = 0.0
             IF (PART2 .LE. 0.0) PART2 = 0.0
             FRAC = (ZENDI(I) - ZBEGI(I) - PART1 - PART2)
     +              /(ZENDI(I) - ZBEGI(I))
         ELSE
             FRAC = 0.0
         ENDIF

C///  ACCUMULATE MASS FOR EACH FOR EACH CELL AND CONVERT FROM
C     (FT**3/MIN TO G/SEC).

         MINJ(CELL) = FRAC*472.01/773.2*VINJ(I) + MINJ(CELL)

C///  BOTTOM OF THE LOOP FOR HEATERS.

 4070 CONTINUE

C///  BOTTOM OF THE LOOP FOR CELLS.

 4090 CONTINUE

      SUM1 = 0.0
      SUM2 = 0.0

      DO 5020 LOC = 1, 2

C///  TOP OF THE LOOP FOR CELLS.

      IF (LOC .EQ. 1) THEN
         CELLS = PIN
         WALLS = CELLS + 1
      ELSE
         CELLS = PEX
         WALLS = CELLS + 1
      ENDIF

      DO 5010 CELL = 1, CELLS
         SUM1 = SUM1 + MINJE(CELL, LOC)
5010  CONTINUE
5020  CONTINUE

      DO 5030 CELL = 1, NUMWF
         SUM2 = SUM2 + MINJ(CELL)
5030  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (4) DETERMINE MASS FLOWRATES AT CELL WALLS.
C
C///////////////////////////////////////////////////////////////////////

C     ZERO MASS FLOWRATES AT CELL WALLS.

      DO 6020 LOC = 1, 2
         DO 6010 K = 1, PDIM+1
            MWZ(K, LOC) = 0.0
 6010    CONTINUE
 6020 CONTINUE
      DO 6030 K = 1, NUMWF+1
         MWZW(K) = 0.0
 6030 CONTINUE

C     FIND THE FLOWRATES AT THE CELL WALLS.

      IF (MFLOW(2) .GE. 0.0) THEN

C        ENTRANCE SECTION
         LOC = 1
         WALLS = PIN + 1
         MWZ(1 ,LOC) = MFLOW(2)
         DO 6040 I = 2, WALLS
            MWZ(I, LOC) = MWZ(I-1, LOC) + MINJE(I-1, LOC)
6040     CONTINUE

C        WAFER SECTION
         MWZW(1) = MWZ(PIN+1, 1)
         WALLS = NUMWF
         DO 6050 I = 2, WALLS
            MWZW(I) = MWZW(I-1) + MINJ(I-1)
6050     CONTINUE

C        EXIT SECTION
         LOC = 2
         MWZ(PEX+1, LOC) = MWZW(NUMWF)
         DO 6060 I = PEX, 1, -1
            MWZ(I, LOC) = MWZ(I+1, LOC) + MINJE(I, LOC)
6060     CONTINUE
      ELSE
      ENDIF

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE NUSSL
     +  (ERROR, TEXT,
     +   DIA, H, COND, RE, PR, NU, TURB)

C///////////////////////////////////////////////////////////////////////
C
C     NUSSL
C
C///  RETURNS NUSSLET NUMBER FOR FLOW BETWEEN CONCENTRIC CLYLINDERS.
C
C     THE ROUTINE RETURNS NU GIVEN THE REYNOLDS NUMBER, PRANDTL NUMBER,
C     AND TURB (PARAMETERS THAT INDICATES LAMINAR OR TURBULENT FLOW).
C     ALSO COMPUTES THE HEAT TRANSFER COEFFICIENT (H) GIVEN THE
C     HYDRAULIC DIAMETER (DIA) AND THERMAL CONDUCTIVITY (COND).
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO
C             AND NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     COND    INPUT REAL - THERMAL CONDUCTIVITY OF COOLING GAS.
C                          (W/CM/K)
C
C     DIA     INPUT REAL - HYDRAULIC DIAMETER OF ANNULUS (CM).
C
C     H       OUTPUT REAL - CONVECTIVE HEAT TRANSFER COEFFICIENT
C                           (W/CM**2/K).
C
C     NU      OUTPUT REAL - NUSSLET NUMBER.
C
C     PR      INPUT REAL - PRANDTL NUMBER.
C
C     RE      INPUT REAL - REYNOLDS NUMBER.
C
C     TURB    INPUT LOGICAL
C             IF TURB = .TRUE.  - USE TURBULENT FLOW CORRELATION.
C             IF TURB = .FALSE. - USE CORRELATION COVERING BOTH LAMINAR
C                                 AND TURBULENT REGIME.
C
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   TEXT

      LOGICAL ERROR, TURB
C*****DOUBLE PRECISION
      DOUBLE PRECISION
     +   DERF,
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C     +   ERF,
C*****END SINGLE PRECISION
     +   COND, DIA, F, H, NU, NU1, NU2, PI, PR, RE, X

      PARAMETER (ID = 'NUSSL:  ')

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE
C
C///////////////////////////////////////////////////////////////////////

      DATA PI /3.1415926/

C///////////////////////////////////////////////////////////////////////
C
C     (2) DETERMINE THE NUSSLET NUMBER.
C
C///////////////////////////////////////////////////////////////////////

C     LAMINAR FLOW CORRELATION (FULLY DEVELOPED FLOW).
      NU1 = 4.86

C     TURBULENT FLOW CORRELATION.
      NU2 = 0.023*RE**(4./5.)*PR**0.4

C///////////////////////////////////////////////////////////////////////
C
C     (3) COMBINED NUSSLET NUMBER EXPRESSION FOR LAMINAR AND TURB.
C
C///////////////////////////////////////////////////////////////////////

      X = 0.0005*(RE-6000.0)

C*****DOUBLE PRECISION
      F = DERF(X)
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      F = ERF(X)
C*****END SINGLE PRECISION

      IF (TURB) THEN
          NU = NU2
      ELSE
          NU = 0.5*(1.-F)*NU1 + 0.5*(1.+F)*NU2
      ENDIF

      H = NU * COND/DIA

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE POW
     +  (ERROR, TEXT,
     +   MAXHT, NUMHT, NUMWF, PDIM, PEX, PIN, POWER, POWERH, PWALL, WZ,
     +   WZWAF, ZBEG, ZEND, ZLEN)

C///////////////////////////////////////////////////////////////////////
C
C     POW
C
C///  COMPUTES THE POWER IN EACH WALL CELL FROM THE HEATERS.
C
C     FOR ANY HEATER "I", THE POWER IS ASSUMED TO BE UNIFORM
C     BETWEEN ZBEG(I) AND ZEND(I) WITH A TOTAL POWER OF
C     POWERH(I) IN WATTS.
C
C     THE ROUTINE RETURNS PWALL(I, LOC) AND POWER(I) THE POWER (IN
C     WATTS) DISSIPATED IN THE WALL REGIONS OF THE ENTRANCE AND EXIT
C     AND ABOVE THE WAFERS.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO
C             AND NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     MAXHT   INPUT INTEGER - MAXIMUM NUMBER OF HEATERS.
C
C     NUMHT   INPUT INTEGER - NUMBER OF HEATERS.
C
C     PEX     INPUT INTEGER - NUMBER OF POINTS IN THE EXIT.
C
C     PIN     INPUT INTEGER - NUMBER OF POINTS IN THE INLET.
C
C     PDIM    INPUT INTEGER - MAXIMUM DIMENSION OF SOME ARRAYS.
C
C     POWER   OUTPUT REAL DIMENSIONED (NUMWF) - WALL POWERS AT GRID
C             POINTS ABOVE THE WAFERS (WATTS).
C
C     POWERH  INPUT REAL DIMENSIONED (MAXHT) - POWERH(I) IS THE TOTAL
C             POWER (WATTS) FROM HEATER "I".
C
C     PWALL   OUTPUT REAL DIMENSIONED (PDIM, 2) - PWALL(I, 1) IS THE
C             POWER DISSIPATED IN THE WALL REGION IN THE ENTRANCE.
C             PWALL(I, 2) IS THE POWER IN THE EXIT REGION.
C
C     WZ      INPUT REAL DIMENSION (PDIM+1, 2) - LOCATION OF CELL
C             WALLS.
C
C     WZWAF   INPUT REAL DIMENSION (NUMWF+1) - LOCATION OF CELL
C             WALLS.
C
C     ZBEG    INPUT REAL DIMENSIONED (MAXHT) - ZBEG(I) IS THE
C             BEGINNING Z COORDINATE FOR HEATER "I".
C
C     ZEND    INPUT REAL DIMENSIONED (MAXHT) - ZEND(I) IS THE ENDING
C             Z COORDINATE FOR HEATER "I".
C
C     ZLEN    INPUT REAL - TOTAL LENGTH OF THE REACTOR (CM).
C
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   CELL, CELLS, I, INNER, K, LOC, MAXHT, NUMHT, NUMWF, OUTER,
     +   PDIM, PEX, PIN, QINNER, QOUTER, TEXT, WALLS

      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   FRAC, PART1, PART2, PI, POWER, POWERH, PWALL, WZ, WZINN, WZOUT,
     +   WZWAF, ZBEG, ZEND, ZLEN

      PARAMETER (ID = 'POW:  ')
      PARAMETER (QINNER = 0, QOUTER = 1)
      DIMENSION
     +   POWER(NUMWF), POWERH(MAXHT), PWALL(PDIM, 2), WZ(PDIM+1, 2),
     +   WZWAF(NUMWF+1), ZBEG(MAXHT), ZEND(MAXHT)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE
C
C///////////////////////////////////////////////////////////////////////

      DATA PI /3.1415926/

C///  ZERO POWERS FOR ALL CELLS.

      DO 1080 LOC = 1, 2
         DO 1070 K = 1, PDIM
            PWALL(K, LOC) = 0.0
 1070    CONTINUE
 1080 CONTINUE
      DO 1090 K = 1, NUMWF
         POWER(K) = 0.0
 1090 CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (2) FIND THE POWER TO EACH CELL IN THE ENTRANCE AND EXIT.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP FOR THE ENTRANCE OR EXIT.

      DO 3095 LOC = 1, 2

C///  TOP OF THE LOOP FOR CELLS.

      IF (LOC .EQ. 1) THEN
         CELLS = PIN
         WALLS = CELLS + 1
      ELSE
         CELLS = PEX
         WALLS = CELLS + 1
      ENDIF

      DO 3090 CELL = 1, CELLS
         INNER = CELL + QINNER
         OUTER = CELL + QOUTER

C///  TOP OF THE LOOP FOR THE HEATERS.

      DO 3070 I = 1, NUMHT

         IF (LOC .EQ. 1) THEN
             WZINN = WZ(INNER, LOC)
             WZOUT = WZ(OUTER, LOC)
         ELSE
C            COORDINATES ARE REVERSED IN THE EXIT
             WZOUT = ZLEN - WZ(INNER, LOC)
             WZINN = ZLEN - WZ(OUTER, LOC)
         ENDIF

         IF ((WZOUT .GE. ZBEG(I)) .AND. (WZINN .LE. ZEND(I))) THEN
             PART1 = WZINN - ZBEG(I)
             PART2 = ZEND(I) - WZOUT
             IF (PART1 .LE. 0.0) PART1 = 0.0
             IF (PART2 .LE. 0.0) PART2 = 0.0
             FRAC = (ZEND(I) - ZBEG(I) - PART1 - PART2)
     +              /(ZEND(I) - ZBEG(I))
         ELSE
             FRAC = 0.0
         ENDIF

C///  ACCUMULATE POWER FOR EACH CELL.

         PWALL(CELL, LOC) = FRAC*POWERH(I) + PWALL(CELL, LOC)

C///  BOTTOM OF THE LOOP FOR POWER.

 3070 CONTINUE

C///  BOTTOM OF THE LOOP FOR CELLS.

 3090 CONTINUE

C///  BOTTOM OF THE LOOP OVER ENTRANCE AND EXIT.

 3095 CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (3) FIND THE POWER TO EACH CELL IN THE WAFER REGION.
C
C///////////////////////////////////////////////////////////////////////

      CELLS = NUMWF

C///  TOP OF THE LOOP FOR CELLS.

      DO 4090 CELL = 1, CELLS
         INNER = CELL + QINNER
         OUTER = CELL + QOUTER

C///  TOP OF THE LOOP FOR HEATERS.

      DO 4070 I = 1, NUMHT

         WZINN = WZWAF(INNER)
         WZOUT = WZWAF(OUTER)

         IF ((WZOUT .GE. ZBEG(I)) .AND. (WZINN .LE. ZEND(I))) THEN
             PART1 = WZINN - ZBEG(I)
             PART2 = ZEND(I) - WZOUT
             IF (PART1 .LE. 0.0) PART1 = 0.0
             IF (PART2 .LE. 0.0) PART2 = 0.0
             FRAC = (ZEND(I) - ZBEG(I) - PART1 - PART2)
     +              /(ZEND(I) - ZBEG(I))
         ELSE
             FRAC = 0.0
         ENDIF

C///  ACCUMULATE POWER FOR EACH CELL.

         POWER(CELL) = FRAC*POWERH(I) + POWER(CELL)

C///  BOTTOM OF THE LOOP FOR HEATERS.

 4070 CONTINUE

C///  BOTTOM OF THE LOOP FOR CELLS.

 4090 CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE RES
     +  (T, S, SP, F, IRES, RPAR, IPAR)

C///////////////////////////////////////////////////////////////////////
C
C     RES
C
C///  RESIDUAL SUBROUTINE FOR DASSL CALLING FORMAT.
C
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO
C             AND NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     INPUT:
C
C     F       REAL DIMENSION (GROUPA + COMPS*PMAX + GROUPB) -
C             RESIDUAL ARRAY.
C
C     IPAR    INTEGER DIMENSION (ISIZE1).
C
C     IRES    INTEGER - INTEGER FLAG EQUAL TO 0 ON INPUT.
C
C     RPAR    REAL DIMENSION (RSIZE1).
C
C     S       REAL DIMENSION (GROUPA + COMPS*PMAX + GROUPB + NADD) -
C             CONTAINS SOLUTION COMPONENTS AT T.
C
C     SP      REAL DIMENSION (GROUPA + COMPS*PMAX + GROUPB + NADD) -
C             CONTAINS DERIVATIVE OF SOLUTION COMPONENTS AT T.
C
C     T       REAL - CURRENT VALUE OF INDEPENDENT VARIABLE (TIME).
C
C     INTERNAL:
C
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   CHOICE, COUNT, I, IPAR, IRES, IWAF, MAXIN, MID, NADD, NUMEX,
     +   NUMHT, NUMIN, NUMWF, OFFSET, PEX, PIN,
     +   RNUM, RPNTS, WNUM, TEXT
      LOGICAL ERROR, CFLOW, RBOAT, RINJEC, SHAPE, TCAPS, TFIX, TIME

C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +    F, RPAR, S, SP, T

C------------------------------------------------------------
C      DON'T FORGET ABOUT THIS INFORMATION IN INPUT TO TWAF1.
CCCC---LOGIALS BFOUND(QTEMP), SHAPE, CHOICE ------------
C------------------------------------------------------------

      PARAMETER (ID = ' RES:  ')

      DIMENSION
     +   F(*), IPAR(*), RPAR(*), S(*), SP(*)

      DIMENSION CFLOW(2), TCAPS(2)

C///  POINTERS FOR RPAR AND IPAR FOR DASSL

      COMMON /RPOINT/
     +   QA2,    QAEX,   QAIN,   QAREAC, QTAMB,  QCOND,
     +   QCONW,  QCP,    QDZEND, QDZNDI, QEMIS,  QEMISS,
     +   QEMISE, QEMISI, QEMISW, QF,     QFEE,
     +   QFEG,   QFENDS, QFENG,  QFESE,  QFEW,
     +   QFGE,   QFGSE,  QFSEE,  QFSEG,  QFSES,
     +   QFSEW,  QFSS,   QFSW,   QFWE,   QFWS,   QFWSE,
     +   QFWW,   QHCOEF, QJ1,    QJ2,    QJEND,  QJEX,
     +   QJIN,   QJWALL, QLOS1,  QLOS2,  QLOS3,  QMFLO,
     +   QMINJ,  QMINJE, QMWZ,   QMWZW,
     +   QPCAP,  QPOWER, QPOWRH, QPWALL, QQWAF,
     +   QR,     QRCPV,  QREFL,  QRHO,   QRING,  QRINQ,
     +   QROUTQ,
     +   QRTUBI, QRTUBO, QBUF,   QX0,    QSPACE,
     +   QTCAP,  QTDAT,  QTEAST, QTEX,
     +   QTEX0,  QTEX4,  QTHICK, QTHCKI,
     +   QTMID,  QTRAN,  QTSIDE, QTSTEP, QTTGAS, QTWAF,
     +   QTWAF0, QTWAF4, QTWEST, QTX,    QTX0,
     +   QTX4,   QWR,    QWZ,    QWZWAF, QZ,
     +   QZBEG,  QZDAT,  QZEND,  QZLEN,  QZWAF,
     +   QGBOT,  QGTOP,  QJBOT,  QJTOP,  QSUMB,
     +   QSUMT,  QTCT,   QVINJ,  QZBEGI, QZENDI, QZTC

      COMMON /IPOINT/
     +   QA2DIM, QCHOIC, QCOMPS, QGROPA,
     +   QGROPB, QIBEG,  QIEND,  QIPV2,  QIPV3,  QIPV4,
     +   QMAXHT, QMAXIN, QNADD,  QNUMEX, QNUMHT, QNUMIN,
     +   QNUMT,  QNUMTC,
     +   QNUMWF, QPDIM,  QPEX,   QPIN,   QPMAX,
     +   QRPNTS, QSHAPE, QTEXT,  QTEMAX, QTFIX,
     +   QTIME

      COMMON /LOGIC/
     +   ERROR, CFLOW, RBOAT, RINJEC, SHAPE, TCAPS, TFIX, TIME

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  DEFINE STATEMENT FUNCTION FOR RADIAL POINT ON A WAFER FOR
C     LOCATION IN S AND SP ARRAYS.

C     WNUM = WAFER NUMBER,  RNUM = RADIAL POINT NUMBER

      IWAF(WNUM, RNUM) = IPAR(QPIN)*IPAR(QNUMEX) + 2
     +                 + (IPAR(QRPNTS)+IPAR(QNUMEX))*(WNUM - 1) + RNUM

C     OFFSET FOR ADDITIONAL VARIABLES

      OFFSET = IPAR(QGROPA) + IPAR(QCOMPS)*IPAR(QPMAX) + IPAR(QGROPB)

C     SET SOME INDICIES
      MAXIN = IPAR(QMAXIN)
      NADD = IPAR(QNADD)
      NUMEX = IPAR(QNUMEX)
      NUMHT = IPAR(QNUMHT)
      NUMIN = IPAR(QNUMIN)
      NUMWF = IPAR(QNUMWF)
      RPNTS = IPAR(QRPNTS)
      PEX = IPAR(QPEX)
      PIN = IPAR(QPIN)

C///  SET POWERS

      IF (NADD.NE.0) THEN

      DO 100 I = 1, NUMHT
         RPAR(QPOWRH + I - 1) = S(OFFSET + I)
100   CONTINUE

      RPAR(QPCAP) = S(OFFSET + NUMHT + 1)
      RPAR(QPCAP + 1) = S(OFFSET + NUMHT + 2)

      COUNT = 0

      ELSE
      ENDIF

C///////////////////////////////////////////////////////////////////////
C
C     (2) CALL THE RESIDUAL.
C
C///////////////////////////////////////////////////////////////////////

C      SUBROUTINE TWAF1
C     +  (ERROR, TEXT,
C     +   A2, A2DIM, AEX, AIN, AREAC, BTEMP, CFLOW, CHOICE,
C     +   COMPS, CONDUC,
C     +   CONW, CP, DZEND, DZENDI, EMIS, EMISS, EMISSE, EMISSI, EMISSW,
C     +   F, FEE, FEG, FENDS, FENG, FESE, FEW, FGE, FGSE, FSEE, FSEG,
C     +   FSES, FSEW, FSS,
C     +   FSW, FWE, FWS, FWSE, FWW, GROUPA, GROUPB, HCOEF, IBEG, IEND,
C     +   IPVT2, IPVT3, IPVT4, J1, J2, JEND, JEX, JIN, JWALL, LOS1,
C     +   LOS2, LOS3, MAXHT, MAXIN, MFLOW, MINJ, MINJE, MWZ, MWZW,
C     +   NUMEX, NUMHT, NUMIN, NUMT, NUMWF, PCAP, PDIM, PEX, PIN,
C     +   PMAX, POWER, POWERH, PWALL, QWAF, R, RBOAT, RCPV, REFL,
C     +   RHO, RING, RINJEC, RINQ, ROUTQ,
C     +   RPNTS, RTUBEI, RTUBEO,
C     +   S, S0, SHAPE, SPACE, TCAP, TCAPS, TDAT, TEAST, TEMMAX, TEX,
C     +   TEX0, TEX4, TFIX, THICK, THICKI, TIME, TMID, TRAN, TSIDE,
C     +   TSTEP, TTGAS, TWAF, TWAF0, TWAF4, TWEST, TX, TX0, TX4,
C     +   VINJ, WR, WZ, WZWAF,
C     +   Z, ZBEG, ZBEGI, ZDAT, ZEND, ZENDI, ZLEN, ZWAF,
C     +   GBOT, GTOP, JBOT, JTOP, NUMTC, SUMB, SUMT, TC, ZTC)

      CHOICE = 4
      TIME = .TRUE.
      SHAPE = .TRUE.
      CALL TWAF1
     +  (ERROR, IPAR(QTEXT),
     +   RPAR(QA2), IPAR(QA2DIM), RPAR(QAEX), RPAR(QAIN),
     +   RPAR(QAREAC), RPAR(QTAMB), CFLOW,
     +   CHOICE, IPAR(QCOMPS), RPAR(QCOND), RPAR(QCONW), RPAR(QCP),
     +   RPAR(QDZEND), RPAR(QDZNDI), RPAR(QEMIS), RPAR(QEMISS),
     +   RPAR(QEMISE), RPAR(QEMISI), RPAR(QEMISW), F,
     +   RPAR(QFEE), RPAR(QFEG), RPAR(QFENDS), RPAR(QFENG),
     +   RPAR(QFESE), RPAR(QFEW), RPAR(QFGE), RPAR(QFGSE), RPAR(QFSEE),
     +   RPAR(QFSEG), RPAR(QFSES), RPAR(QFSEW), RPAR(QFSS),
     +   RPAR(QFSW), RPAR(QFWE), RPAR(QFWS), RPAR(QFWSE), RPAR(QFWW),
     +   IPAR(QGROPA), IPAR(QGROPB), RPAR(QHCOEF), IPAR(QIBEG),
     +   IPAR(QIEND), IPAR(QIPV2), IPAR(QIPV3), IPAR(QIPV4),
     +   RPAR(QJ1), RPAR(QJ2), RPAR(QJEND), RPAR(QJEX), RPAR(QJIN),
     +   RPAR(QJWALL), RPAR(QLOS1), RPAR(QLOS2), RPAR(QLOS3),
     +   IPAR(QMAXHT), IPAR(QMAXIN), RPAR(QMFLO),
     +   RPAR(QMINJ), RPAR(QMINJE), RPAR(QMWZ), RPAR(QMWZW),
     +   IPAR(QNUMEX), IPAR(QNUMHT), IPAR(QNUMIN), IPAR(QNUMT),
     +   IPAR(QNUMWF),
     +   RPAR(QPCAP), IPAR(QPDIM), IPAR(QPEX), IPAR(QPIN),
     +   IPAR(QPMAX), RPAR(QPOWER), RPAR(QPOWRH), RPAR(QPWALL),
     +   RPAR(QQWAF), RPAR(QR), RBOAT, RPAR(QRCPV), RPAR(QREFL),
     +   RPAR(QRHO), RPAR(QRING), RINJEC, RPAR(QRINQ), RPAR(QROUTQ),
     +   IPAR(QRPNTS), RPAR(QRTUBI), RPAR(QRTUBO),
     +   S, SP, SHAPE, RPAR(QSPACE), RPAR(QTCAP), TCAPS,
     +   RPAR(QTDAT), RPAR(QTEAST), IPAR(QTEMAX),
     +   RPAR(QTEX), RPAR(QTEX0), RPAR(QTEX4), TFIX,
     +   RPAR(QTHICK), RPAR(QTHCKI), TIME, RPAR(QTMID), RPAR(QTRAN),
     +   RPAR(QTSIDE), RPAR(QTSTEP), RPAR(QTTGAS), RPAR(QTWAF),
     +   RPAR(QTWAF0), RPAR(QTWAF4), RPAR(QTWEST), RPAR(QTX),
     +   RPAR(QTX0), RPAR(QTX4), RPAR(QVINJ),
     +   RPAR(QWR), RPAR(QWZ), RPAR(QWZWAF), RPAR(QZ),
     +   RPAR(QZBEG), RPAR(QZBEGI), RPAR(QZDAT), RPAR(QZEND),
     +   RPAR(QZENDI), RPAR(QZLEN),
     +   RPAR(QZWAF), RPAR(QGBOT), RPAR(QGTOP), RPAR(QJBOT),
     +   RPAR(QJTOP), IPAR(QNUMTC), RPAR(QSUMB), RPAR(QSUMT),
     +   RPAR(QTCT), RPAR(QZTC))

      IF (ERROR) GO TO 9101

      MID = (NUMWF + 1)/2

C///  SET THE ADDITIONAL VARIABLES IF NECESSARY

      IF (NADD .NE. 0) THEN

      DO 200 I = 1, NADD
         F(OFFSET + I) = SP(OFFSET + I) - 1.0*(1. - EXP(-10.*T))
200   CONTINUE

      I = IWAF(MID, RPNTS+3)

C       F(OFFSET + 1) = SP(IWAF(1, RPNTS))
C     +                 - 100.0/60.0*(1. - EXP(-10.*T))
C       F(OFFSET + 2) = SP(IWAF(MID, RPNTS))
C     +                 - 100.0/60.*(1. - EXP(-10.*T))
C       F(OFFSET + 3) = SP(IWAF(1, RPNTS))
C     +                 - 100.0/60.0*(1. - EXP(-10.*T))
C       F(OFFSET + 4) = SP(IWAF(NUMWF, RPNTS))
C     +                 - 100.0/60.0*(1. - EXP(-10.*T))

      ELSE
      ENDIF

C///////////////////////////////////////////////////////////////////////
C
C     (3) EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////


C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID
      GO TO 99999

99101 FORMAT
     +   (/1X, A, 'ERROR. PROBLEM WITH SUBROUTINE RES.')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE SCONSW
     +  (ERROR, TEXT,
     +   F12, F21, RADIUS, REFLEC, RIN, ROUT, SPACE)

C///////////////////////////////////////////////////////////////////////
C
C     SCONSW
C
C///  COMPUTES THE CONFIGUATION FACTOR BETWEEN THE FINITE RING AREA
C     BOUNDED BY TWO WAFERS ON THE INTERIOR SIDE OF A RIGHT CIRCULAR
C     CYLINDER TO AN ANNULAR RING ON ONE OF THE WAFERS ACCOUNTING FOR
C     MULTIPLE REFLECTIONS BY THE WAFERS.  IN THE LIMIT AS RIN = 0,
C     THE CONFIGURATION FACTOR FROM THE RING ELEMENT ON THE CYLINDER
C     SIDE TO A DISK OF RADIUS ROUT IS COMPUTED.
C
C     AREA1 = RING AREA ON INTERIOR OF THE CYLINDER.
C     AREA2 = ANNULAR RING OR DISK ON WAFER.
C
C     FORMULA C-77 FROM "A CATALOG OF RADIATION CONFIGURATION FACTORS",
C     J.R. HOWELL, IS USED IN CONJUNCTION WITH CONFIGURATION FACTOR
C     ALGEBRA.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     RADIUS  INPUT REAL - RADIUS OF CYLINDER.
C
C     REFLEC  INPUT REAL - REFLECTIVITY OF WAFERS.
C
C     RIN     INPUT REAL - INNER RADIUS OF ANNULAR RING IN BASE OF
C                          CYLINDER.
C
C     ROUT    INPUT REAL - OUTER RADIUS OF ANNULAR RING IN BASE OF
C                          CYLINDER.
C
C     SPACE   INPUT REAL - SPACE BETWEEN WAFERS
C
C     F12     OUTPUT REAL - SHAPE FACTOR FROM RING ON CYLINDER SIDE TO
C                           ANNULAR RING IN BASE OF CYLINDER.
C
C     F21     OUTPUT REAL - SHAPE FACTOR FROM ANNULAR RING IN BASE OF
C                           CYLINDER TO RING ON CYLINDER SIDE.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   K, NUMTOT, TEXT

      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   A1, A2, F12, F21, FINNER, FOUTER, LLOW, LUPP, PI,
     +   RADIUS, REFLEC, RIN, ROUT, SPACE, SUM

      PARAMETER (ID = 'SCONSW:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)
      PARAMETER (NUMTOT = 150)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE
C
C///////////////////////////////////////////////////////////////////////

C///  CHECKS THE ARGUMENTS.

      ERROR = (RIN .GE. ROUT)
      IF (ERROR) GO TO 9001

      ERROR = (SPACE .LE. 0.0)
      IF (ERROR) GO TO 9002

      ERROR = ((RIN .GE. RADIUS) .OR. (ROUT .GT. RADIUS))
      IF (ERROR) GO TO 9003

C///////////////////////////////////////////////////////////////////////
C
C     (2) FIND THE CONFIGURATION FACTORS.
C
C///////////////////////////////////////////////////////////////////////

      IF ((ROUT .EQ. 0.0) .OR. (SPACE .EQ. 0.0)) THEN
         F12 = 0.0
         F21 = 0.0
         RETURN
      ELSE
      ENDIF

      SUM = 0.0

C///  AREAS

      A1 = 2.*RADIUS*PI*SPACE
      A2 = PI*(ROUT**2 - RIN**2)

      DO 100 K = 0, NUMTOT
         LLOW = FLOAT(K)*SPACE
         LUPP = FLOAT(K+1)*SPACE

C///  CONFIGURATION FACTOR TO THE OUTER DISK AND INNER DISK

         FOUTER = 1./(4.*RADIUS/ROUT*(LUPP/ROUT - LLOW/ROUT))*
     +            (((LLOW/ROUT)**2 - (LUPP/ROUT)**2) -
     +            SQRT(((LLOW/ROUT)**2 +
     +            (RADIUS/ROUT)**2 + 1.)**2 - 4.*(RADIUS/ROUT)**2) +
     +            SQRT(((LUPP/ROUT)**2 +
     +            (RADIUS/ROUT)**2 + 1.)**2 - 4.*(RADIUS/ROUT)**2))

         IF (RIN .EQ. 0.0) THEN
            FINNER = 0.0
         ELSE
            FINNER = 1./(4.*RADIUS/RIN*(LUPP/RIN - LLOW/RIN))*
     +               (((LLOW/RIN)**2 - (LUPP/RIN)**2) -
     +               SQRT(((LLOW/RIN)**2 +
     +               (RADIUS/RIN)**2 + 1.)**2 - 4.*(RADIUS/RIN)**2) +
     +               SQRT(((LUPP/RIN)**2 +
     +               (RADIUS/RIN)**2 + 1.)**2 - 4.*(RADIUS/RIN)**2))
         ENDIF

         SUM = SUM + (REFLEC**K)*(FOUTER - FINNER)

100   CONTINUE

C///  FINAL CONFIGURATION FACTORS

      F12 = SUM
      F21 = A1/A2*F12

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, RIN, ROUT

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, SPACE

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, RADIUS, RIN, ROUT

      GO TO 99999

99001 FORMAT
     +   (/1X, A, 'ERROR. INNER RADIUS IS GREATER THAN OR EQUAL',
     +    /10X, 'TO OUTER RADIUS.',
     +   //10X, 'RIN = ',E10.3,
     +    /10X, 'ROUT = ',E10.3)
99002 FORMAT
     +   (/1X, A, 'ERROR. SPACE BETWEEN WAFERS IS ZERO OR LESS',
     +   //10X, 'SPACE = ',E10.3)
99003 FORMAT
     +   (/1X, A, 'ERROR. INNER OR OUTER RADIUS OF ANNULAR RING IS'
     +    /10X, 'GREATER THAN OR EQUAL THE RADIUS OF THE CYLINDER.',
     +   //10X, 'RADIUS = ',E10.3,
     +    /10X, 'RIN = ',E10.3,
     +    /10X, 'ROUT = ',E10.3)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE SCONWW
     +  (ERROR, TEXT,
     +   ASTAR, F12, F21, REFLEC, RIN1, RIN2, ROUT1, ROUT2, SAME)

C///////////////////////////////////////////////////////////////////////
C
C     SCONWW
C
C///  COMPUTES THE CONFIGUATION FACTOR BETWEEN TWO PARALLEL, DIRECTLY
C     APPOSED PLANE RING AREAS WITH SPECULAR REFLECTION. IN THE LIMIT
C     AS RIN1 AND/OR RIN2 GOES TO ZERO, THE CONFIGURATION FACTOR
C     BETWEEN A RING AND A DISK OR TWO DISKS IS COMPUTED. DIFFUSE
C     EXCHANGE FACTOR CAN BE CALCULATED BY SETTING REFLEC = 0.0.
C
C
C     CONFIGURATION FACTORS CALCULATED USING SHAPE FACTOR ALGEBRA AND
C     REFERENCE C-35 FROM "A CATALOG OF RADIATION CONFIGURATION
C     FACTORS", J.R. HOWELL.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     ASTAR   INPUT REAL - ON AXIS DISTANCE BETWEEN TWO RINGS.
C
C     REFLEC  INPUT REAL - REFLECTIVITY OF WAFERS
C
C     RIN1    INPUT REAL - INNER RADIUS OF RING 1.
C
C     RIN2    INPUT REAL - INNER RADIOUS OF RING 2.
C
C     ROUT1   INPUT REAL - OUTER RADIUS OF RING 1.
C
C     ROUT2   INPUT REAL - OUTER RADIUS OF RING 2.
C
C     SAME    INPUT LOGICAL - IF SAME = .TRUE. THEN EXCHANGE FACTOR IS
C                             FOR RING ELEMENTS WITHIN THE SAME WAFER
C                             IF SAME = .FALSE THEN EXCHANGE FACTOR IS
C                             FOR RING ELEMENTS IN ADJACENT WAFERS
C
C     F12     OUTPUT REAL - SHAPE FACTOR FROM RING 1 TO RING 2
C
C     F21     OUTPUT REAL - SHAPE FACTOR FROM RING 2 TO RING 1
C
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   K, KSTART, NUMTOT, TEXT

      LOGICAL ERROR, SAME
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   A, ARING1, ARING2, ASTAR, F12, F21, FL1L2, FL1S2, FS1L2,
     +   FS1S2, PI, REFLEC, RIN1, RIN2, ROUT1, ROUT2, SUM

      PARAMETER (ID = 'CONWW:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)
      PARAMETER (NUMTOT = 150)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////


C///  CHECKS

      ERROR = (RIN1 .GE. ROUT1)
      IF (ERROR) GO TO 9001

      ERROR = (RIN2 .GE. ROUT2)
      IF (ERROR) GO TO 9002

      ERROR = (ASTAR .LE. 0.0)
      IF (ERROR) GO TO 9003

C///  CHECK TO SEE IF CALCULATION IS FOR RINGS WITHIN THE SAME WAFER
C     OR FOR RINGS WITHIN ADJACENT WAFERS

      IF (SAME) THEN
         KSTART = 1
      ELSE
         KSTART = 0
      ENDIF

C///////////////////////////////////////////////////////////////////////
C
C     (2) FIND THE CONFIGURATION FACTORS
C
C///////////////////////////////////////////////////////////////////////

      IF ((ROUT1 .EQ. 0.0) .OR. (ROUT2 .EQ. 0.0)) THEN
         F12 = 0.0
         F21 = 0.0
         RETURN
      ELSE
      ENDIF

      SUM = 0.0

C///  AREAS

      ARING1 = PI*(ROUT1**2 - RIN1**2)
      ARING2 = PI*(ROUT2**2 - RIN2**2)

C///  TOP OF THE LOOP OVER REFLECTIONS

      DO 100 K = KSTART, NUMTOT, 2

         A = ASTAR*FLOAT(K+1)

C///  CONFIGURATION FACTOR FROM LARGE DISK1 TO LARGE DISK2 (FL1L2) AND
C     LARGE DISK1 TO SMALL DISK2 (FL1S1)

         FL1L2 = 0.5*(1. + (1. + (ROUT2/A)**2)/(ROUT1/A)**2 -
     +           SQRT((1. + (1. + (ROUT2/A)**2)/(ROUT1/A)**2)**2 -
     +           4.*(ROUT2/ROUT1)**2))
         FL1S2 = 0.5*(1. + (1. + (RIN2/A)**2)/(ROUT1/A)**2 -
     +           SQRT((1. + (1. + (RIN2/A)**2)/(ROUT1/A)**2)**2 -
     +           4.*(RIN2/ROUT1)**2))


C///  CONFIGURATION FACTOR FROM SMALL DISK1 TO LARGE DISK2 (FS1L1) AND
C     SMALL DISK1 TO SMALL DISK2 (FS1S1)

         IF (RIN1 .EQ. 0.0) THEN
            FS1L2 = 0.0
            FS1S2 = 0.0
         ELSE
            FS1L2 = 0.5*(1. + (1. + (ROUT2/A)**2)/(RIN1/A)**2 -
     +              SQRT((1. + (1. + (ROUT2/A)**2)/(RIN1/A)**2)**2 -
     +              4.*(ROUT2/RIN1)**2))
            FS1S2 = 0.5*(1. + (1. + (RIN2/A)**2)/(RIN1/A)**2 -
     +              SQRT((1. + (1. + (RIN2/A)**2)/(RIN1/A)**2)**2 -
     +              4.*(RIN2/RIN1)**2))
         ENDIF

C///  FINAL CONFIGURATION FACTORS

         SUM = SUM + (REFLEC**K)*((PI*RIN1**2 + ARING1)/ARING1*(FL1L2 -
     +         FL1S2) - (PI*RIN1**2)/ARING1*(FS1L2 - FS1S2))

100   CONTINUE

      F12 = SUM
      F21 = ARING1/ARING2*F12

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, RIN1, ROUT1

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, RIN2, ROUT2

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, ASTAR

      GO TO 99999

99001 FORMAT
     +   (/1X, A, 'ERROR. RIN1 GREATER THAN OR EQUAL TO ROUT1.'
     +   //10X, 'RIN1 = ',E10.3,
     +    /10X, 'ROUT1 = ',E10.3)

99002 FORMAT
     +   (/1X, A, 'ERROR. RIN2 GREATER THAN OR EQUAL TO ROUT2.'
     +   //10X, 'RIN2 = ',E10.3,
     +    /10X, 'ROUT2 = ',E10.3)

99003 FORMAT
     +   (/1X, A, 'ERROR. SPACING BETWEEN DISKS IS ZERO OR LESS.'
     +   //10X, 'ASTAR = ', E10.3)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWAF0
     +  (ERROR, TEXT,
     +   A2, A2DIM, AEX, AIN, CFLOW, COMPS, EMIS, EMISS,
     +   EMISSE, EMISSW, FEE, FEG, FESE,
     +   FEW, FGE, FGSE, FSEE, FSEG, FSES, FSEW, FSS, FSW, FWE, FWS,
     +   FWSE, FWW, GROUPA, GROUPB, IBEG, IEND, IPVT2, IPVT3,
     +   IPVT4, JEND, JWALL, NUMEX, NUMT, NUMWF, PDIM, PEX, PIN,
     +   PMAX, R, RADIUS, RBOAT, REFL, RING, RINQ, ROUTQ, RPNTS,
     +   RTUBEI,
     +   S, S0, SHAPE, SPACE, TDAT, TEMMAX, TFIXED,
     +   TMID, TRAN,
     +   TSIDE, TTEND, TTGAS, TTIN, TTQ, TTUBE, TTWAF,
     +   WORK2, WORK3, WORK4, WR, WZ,
     +   WZWAF, Z, ZDAT, ZFIRST, ZLAST, ZLEN, ZWAF)

C///////////////////////////////////////////////////////////////////////
C
C     TWAF0
C
C///  PREPARE TO SOLVE THE THERMAL PROBLEM.
C     COMPUTE GRID, SHAPE FACTORS, AND AN INITIAL GUESS FOR THE
C     SOLUTION ARRAY.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     INPUT:
C
C     COMPS   INPUT INTEGER - NUMBER OF VARIABLES AT EACH CELL.
C
C     CFLOW   INPUT LOGICAL DIMENSION (2) - LOGICAL VARIABLE TO
C             IDENTIFY IF GAS FLOW IS ACTIVE.
C             CFLOW(1) = .TRUE. - COOLING FLOW BETWEEN QUARTZ TUBES.
C             CLFOW(2) = .TRUE. - COOLING FLOW BETWEEN HEATER AND
C                                 OUTER QUARTZ TUBE.
C
C     EMIS    INPUT REAL DIMENSIONED (2) - EMIS(1) AND EMIS(2) ARE
C             THE EMITTANCES OF QUARTZ TUBES 1 AND 2.
C
C     EMISS   INPUT REAL - EMISSIVITY OF WAFERS.
C
C     EMISSE  INPUT REAL - EMISSIVITY OF REACTOR END WALL.
C
C     EMISSW  INPUT REAL - EMISSIVITY OF REACTOR HEATERS.
C
C     NUMT    INPUT INTERER - NUMBER OF TEMPERATURE DATA POINTS ALONG
C             WALL.
C
C     NUMWF   INPUT INTEGER - NUMBER OF SIMULATED WAFERS.
C
C     PDIM    INPUT INTEGER - MAXIMUM DIMENSION OF SOME ARRAYS.
C
C     PIN     INPUT INTEGER - NUMBER OF POINTS BETWEEN REACTOR END WALL
C             AND FIRST WAFER IN THE ENTRANCE REGION.
C
C     PEX     INPUT INTEGER - NUMBER OF POINTS BETWEEN REACTOR END WALL
C             AND LAST WAFER IN THE EXIT REGION.
C
C     PMAX    INPUT INTEGER - NUMBER OF POSITIONS AT WHICH SOLUTION IS
C             COMPUTED.
C
C     R       INPUT REAL DIMENSIONED (RPNTS) - RADIAL GRID POINTS IN
C             WAFERS.
C
C     RADIUS  INPUT REAL - RADIUS OF THE WAFERS.
C
C     RBOAT   INPUT LOGICAL - IF .TRUE. INCLUDE RING BOAT.
C
C     REFL    INPUT REAL DIMENSIONED (2) - REFL(1) AND REFL(2) ARE THE
C             REFLECTANCES OF QUARTZ TUBES 1 AND 2.
C
C     RING    INPUT REAL DIMENSIONED (5) - RING BOAT PARAMETERS.
C             RING(1) - RING BOAT INSIDE RADIUS.
C             RING(2) - RING BOAT THICKNESS.
C             RING(3) - DENISITY OF RING BOAT MATERIAL.
C             RING(4) - SPECIFIC HEAT OF RING BOAT MATERIAL.
C             RING(5) - THERMAL CONDUCTIVITY OF RING BOAT MATERIAL.
C
C     RINQ    INPUT REAL DIMENSIONED (2) - INNER RADII OF QUARTZ TUBES
C             1 AND 2.
C
C     ROUTQ   INPUT REAL DIMENSIONED (2) - OUTER RADII OF QUARTZ TUBES
C             1 AND 2.
C
C     RPNTS   INPUT INTEGER - NUMBER OF RADIAL POINTS IN A WAFER.
C
C     RTUBEI  INPUT REAL - INNER RADIUS OF THE REACTOR TUBE.
C
C     SHAPE   INPUT LOGICAL - IF .TRUE. PREFORM A SHAPE FACTOR
C             CALCULATION BEFORE COMPUTING THE RESIDUAL.
C
C     SPACE   INPUT REAL - SPACING BETWEEN WAFERS
C
C     TDAT    INPUT REAL DIMENSIONED (NUMT) - TEMPERATURE DATA POINTS
C             ALONG REACTOR WALL AT ZDAT(I) LOCATIONS.
C
C     TEND    REAL DIMENSIONED (2) - TEND(1) AND TEND(2) ARE THE
C             RESPECTIVE TEMPERATURES OF THE END WALLS IN THE ENTRANCE
C             AND EXIT OF THE REACTOR.
C
C     TRAN    REAL DIMENSIONED (2) - TRAN(1) AND TRAN(2) ARE THE
C             TRANSMITTANCES OF QUARTZ TUBES 1 AND 2.
C
C     TTIN    INPUT REAL - INITIAL INSULATION TEMPERATURE.
C
C     TTEND   INPUT REAL - INITIAL END TEMPERATURES OF TUBE.
C
C     TTGAS   INPUT REAL DIMENSIONED (2) - INTIAL TEMPERATURE OF GAS.
C             TTGAS(1) - INITIAL GAS TEMPERATURE BETWEEN HEATER
C                        AND OUTER QUARTZ TUBE.
C             TTGAS(2) - INITIAL GAS TEMPERATURE BETWEEN QUARTZ TUBES.
C
C     TTUBE   INPUT REAL - INITIAL TUBE (HEATER) TEMPERATURE.
C
C     TTQ     INPUT REAL - INITIAL TEMPERATURE OF QUARTZ WALLS.
C
C     TTWAF   INPUT REAL - INITIAL WAFER TEMPERATURE.
C
C     Z       INPUT REAL DIMENSIONED (PDIM, 2) - Z(I, 1) ARE THE AXIAL
C             GRID POINTS BETWEEN THE REACTOR END WALL AND FIRST WAFER.
C
C     ZDAT    INPUT REAL DIMENSIONED (NUMT) - AXIAL LOCATIONS OF
C             TEMPERATURE DATA POINTS ALONG REACTOR WALL.
C
C     ZFIRST  INPUT REAL - POSITION OF FIRST WAFER.
C
C     ZLAST   INPUT REAL - POSITION OF LAST WAFER.
C
C     ZLEN    INPUT REAL - LENGTH OF THE REACTOR.
C
C     ZWAF    INPUT REAL DIMEMSIONED (NUMWF)  - ZWAF(I) IS THE AXIAL
C             LOCATION OF SIMULATED WAFER I.
C
C     OUTPUT:
C
C     FEE     REAL DIMENSIONED (2) - SHAPE FACTOR FROM REACTOR END WALL
C             TO ITSELF.
C
C     FEG     REAL DIMENSIONED (2) - FEG(1) IS THE SHAPE FACTOR FROM
C             REACTOR END WALL TO THE GAP BETWEEN THE FIRST WAFER AND
C             THE WALL.
C
C     FESE    REAL DIMENSIONED (PDIM, 2) - FESE(I, 1) ARE THE SHAPE
C             FACTORS FROM THE REACTOR END WALL TO ELEMENTS I ON THE
C             WALL OF THE REACTOR ENTRANCE.
C
C     FEW     REAL DIMENSIONED (RPNTS, 2) - FEW(I, 1) ARE THE SHAPE
C             FACTORS FROM THE RING ELEMENTS I ON THE FIRST WAFER TO
C             THE REACTOR END WALL.
C
C     FGE     REAL DIMENSIONED (2) - FGE(1) IS THE SHAPE FACTOR FROM
C             THE GAP BETWEEN THE FIRST WAFER AND THE WALL TO THE END
C             WALL OF THE REACTOR.
C
C     FGSE    REAL DIMENSIONED (PDIM, 2) - FGSE(I, 1) ARE THE SHAPE
C             FACTORS FROM THE GAP BETWEEN THE FIRST WAFER AND THE WALL
C             TO THE ELEMENTS I ON THE SIDE WALL IN THE ENTRANCE REGION.
C
C     FSEE    REAL DIMENSIONED (PDIM, 2) - SHAPE FACTORS FROM
C             ELEMENTS ON REACTOR SIDE WALL TO REACTOR END WALL.
C
C     FSES    REAL DIMENSIONED (PDIM, PDIM, 2) - SHAPE FACTORS
C             BETWEEN ELEMENTS ON REACTOR SIDE WALL IN THE ENTRANCE
C             REGION.
C
C     FSEG    REAL DIMENSIONED (PDIM, 2) - FSEG(I, 1) ARE THE SHAPE
C             FACTORS FROM THE ELEMENTS I ON THE ENTRANCE WALL OF THE
C             REACTOR TO THE GAP BETWEEN THE FIRST WAFER AND THE WALL.
C
C     FSEW    REAL DIMENSIONED (PDIM, RPNTS, 2) - FSEW(I,J) ARE THE
C             SHAPE FACTORS FROM ELEMENTS I ON THE WALL OF THE
C             REACTOR ENTRANCE TO RING ELEMENTS J ON THE FIRST WAFER.
C
C     FSS     REAL - CONFIGURATION FACTOR FROM THE REACTOR SIDE
C             WALL AREA BETWEEN TWO WAFERS TO ITSELF.
C
C     FSW     REAL DIMENSIONED (2*RPNTS) - FSW(I) ARE THE SHAPE
C             FACTORS FROM THE REACTOR WALL AREA BETWEEN TWO WAFERS TO
C             THE RING ELEMENTS ON THE WAFER SURFACES.
C
C     FWE     REAL DIMENSIONED (RPNTS, 2) - FWE(I, 1) ARE THE SHAPE
C             FACTORS FROM THE RING ELEMENTS I ON THE FIRST WAFER TO
C             THE REACTOR END WALL.
C
C     FWS     REAL DIMENSIONED (2*RPNTS) - FWS(I) ARE THE SHAPE
C             FACTORS FROM THE RING ELEMENTS ON THE WAFER SURFACES TO
C             THE REACTOR WALL AREA BETWEEN TWO WAFERS.
C
C     FWSE    REAL DIMENSIONED (RPNTS, PDIM, 2) - FWSE(I, J, 1) ARE
C             THE SHAPE FACTORS FROM RING ELEMENTS I ON THE FIRST
C             WAFER TO ELEMENTS J ON THE WALL OF THE REACTOR ENTRANCE.
C
C     FWW     REAL DIMENSIONED (2*RPNTS, 2*RPNTS) - FWW(I,J) ARE
C             THE SHAPE FACTORS FROM WAFER RING ELEMENT I TO WAFER RING
C             ELEMENT J. ELEMENTS 1 TO RPNTS ARE ON THE LEFT WAFER FROM
C             THE CENTERLINE TO THE EDGE. ELEMENTS RPNTS+1 TO 2*RPNTS
C             ARE ON THE RIGHT WAFER FROM THE CENTERLINE TO THE EDGE.
C
C     JEND    REAL DIMENSIONED (2) - FEND(1) AND FEND(2) ARE THE
C             RESPECTIVE RADIOSITIES FROM THE REACTOR ENDS AT THE
C             ENTRANCE AND EXIT.
C
C     JWALL   REAL DIMENSIONED (PDIM, 2) - JWALL(I, 1) ARE VALUES OF
C             THE SIDE-WALL RADIOSITIES IN THE ENTRANCE REGION.
C             JWALL(I, 2) ARE VALUES OF THE SIDE-WALL RADIOSITIES
C             IN THE EXIT REGION. JWALL(1,1) AND JWALL(1,2) ARE NEXT
C             TO THE END WALLS.
C
C     S       OUTPUT REAL DIMENSIONED (GROUPA + COMPS*PMAX + GROUPB) -
C             INITIAL GUESS FOR SOLUTION.
C
C     S0      OUTPUT REAL DIMENSIONED (GROUPA + COMPS*PMAX + GROUPB) -
C             INTIAL GUESS FOR SOLUTION AT OLD TIME.
C
C     SIG     REAL - STEFAN BOLTZMAN CONSTANT (W/CM**2/K**4).
C
C     WR      REAL DIMENSIONED (RPNTS+1) - CELL WALL LOCATIONS
C             IN WAFERS.
C
C     WZ      REAL DIMENSIONED (PDIM+1, 2) - WZ(I, 1) ARE THE CELL WALL
C             LOCATIONS BETWEEN REACTOR END WALL AND FIRST WAFER.
C
C     WZWAF   REAL DIMENSIONED (NUMWF+1) - WZWAF(I) ARE THE CELL WALL
C             LOCATIONS IN THE WAFER LOAD.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID * 9
      INTEGER
     +   A2DIM, AINDIM, AEXDIM, CELL, CELLS, COMPS, COUNT,
     +   GROUPA, GROUPB,
     +   IBEG, IEND, IOPT, IPOINT,
     +   IPVT2, IPVT3, IPVT4, I, J, K, KRON, LOC, NUMEX,
     +   NUMT, NUMWF, PDIM, PEX, PIN, PMAX, RPNTS, RWALLS, TEMMAX,
     +   TEXT, WALL, ZWALLS
      LOGICAL CFLOW, ERROR, RBOAT, SHAPE, TFIXED
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   A1P, A1M, A2, A2P, A2M,  AEX, AH, AIN, AREA, AREA1, AWAF,
     +   EMIS, EMISS,
     +   EMISSE, EMISSW, F12W, AREAW,
     +   F12, F21, F21W, F22, F22W, F2H, FH2, FHH,
     +   FEE, FEG,
     +   FESE, FEW, FGE, FGSE, FSEE, FSEG, FSES, FSEW, FSS, FSW, FWE,
     +   FWS, FWSE, FWW, JEND, JWALL, MIN, OFFSET, PI, R,
     +   RADIUS, REFL, REFLEC, RING,
     +   RCOND1, RCOND2, RINQ, ROUTQ, RTUBEI, S, S0, SIG, SPACE,
     +   SUM1, SUM2, TDAT, TEND, TMID, TRAN, TSIDE, TTEND, TTGAS,
     +   TTIN, TTQ, TTUBE, TTWAF, WORK2, WORK3, WORK4, WR,
     +   WZ, WZWAF, Z, ZDAT, ZFIRST, ZLAST, ZLEN, ZTEST, ZWAF

      PARAMETER (ID = ' TWAF0:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)
      PARAMETER (QINNER = 0, QOUTER = 1)
      PARAMETER (SIG = 5.6696E-12)

      DIMENSION
     +   A2(A2DIM, A2DIM),
     +   AEX(5*PEX+2, 5*PEX+2), AIN(5*PIN+2, 5*PIN+2),
     +   AREA(2), AREAW(2), CFLOW(2),
     +   EMIS(2), FEE(2), FEG(2),
     +   FESE(PDIM, 2), FEW(RPNTS, 2), FGE(2),
     +   FGSE(PDIM, 2), FSEE(PDIM, 2), FSEG(PDIM, 2),
     +   FSES(PDIM, PDIM, 2), FSEW(PDIM, RPNTS, 2), FSW(2*RPNTS),
     +   FWE(RPNTS, 2), FWS(2*RPNTS), FWSE(RPNTS, PDIM, 2),
     +   FWW(2*RPNTS, 2*RPNTS), IBEG(3), IEND(3),
     +   IPVT2(A2DIM), IPVT3(5*PIN+2), IPVT4(5*PEX+2), JEND(2),
     +   JWALL(PDIM, 2), OFFSET(2), REFL(2), R(RPNTS), RING(5),
     +   RINQ(2), ROUTQ(2),
     +   S(GROUPA + COMPS*PMAX + GROUPB),
     +   S0(GROUPA + COMPS*PMAX + GROUPB), TDAT(TEMMAX + 2),
     +   TEND(2), TMID(NUMWF), TRAN(2), TSIDE(PDIM, 2), TTGAS(2),
     +   WORK2(A2DIM), WORK3(5*PIN+2), WORK4(5*PEX+2),
     +   WR(RPNTS+1), WZ(PDIM+1, 2),
     +   WZWAF(NUMWF+1), Z(PDIM, 2), ZDAT(TEMMAX + 2), ZWAF(NUMWF)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  DEFINE STATEMENT FUNCTION FOR KRONECKER DELTA

      KRON(I, J) = MAX(-1, -ABS(I - J)) + 1

      I = KRON(3, 5)
      J = KRON(3, 3)

      AINDIM = 5*PIN+2
      AEXDIM = 5*PEX+2

C///  R GRID ACROSS WAFERS

      ERROR = .NOT. (2 .LE. RPNTS)
      IF (ERROR) GO TO 9405

      DO 1100 J = 1, RPNTS
         R(J) =  RADIUS * (REAL (J - 1) / REAL (RPNTS - 1))
1100  CONTINUE

      RWALLS = RPNTS + 1

C///  LOCATION OF CELL WALLS IN WAFERS.

      DO 1200 WALL = 1, RWALLS
         IF (WALL .EQ. 1) THEN
            CELL = WALL - QINNER
            WR(WALL) = R(CELL)
         ELSE IF (WALL .EQ. RWALLS) THEN
            CELL = WALL - QOUTER
            WR(WALL) = R(CELL)
         ELSE
            CELL = WALL - QOUTER
            WR(WALL) = 0.5 * (R(CELL) + R(CELL + 1))
         ENDIF
1200  CONTINUE

C///  IF RING BOAT IS PRESENT MODIFY LOCATION OF NEAREST CELL WALL
C     TO MATCH INNER RADIUS OF RING BOAT.

      MIN = R(RPNTS)
      IPOINT = 2

      IF (RBOAT) THEN
         DO 1250 WALL = 2, RWALLS-1
            IF (ABS(WR(WALL) - RING(1)) .LE. MIN) THEN
               IPOINT = WALL
               MIN = ABS(WR(WALL) - RING(1))
            ELSE
            ENDIF
1250     CONTINUE
         IF(R(1) .EQ. RING(1)) IPOINT = 1

C        MOVE THE CLOSEST CELL BOUNDARY TO THE INNER RADIUS OF THE
C        RING BOAT

         WR(IPOINT) = RING(1)

C        REPOSITION WAFER POINTS AT THE CENTER OF CELLS.

         DO 1275 J = 2, RPNTS-1
            R(J) = (WR(J+1) + WR(J))/2.
1275     CONTINUE
      ELSE
      ENDIF

C///  CELL WALLS AND POINTS FOR THE ENTRANCE REGION

      ERROR = .NOT. (2 .LE. PIN)
      IF (ERROR) GO TO 9402

      LOC = 1
      ZWALLS = PIN + 1

      DO 1300 WALL = 1, ZWALLS
         WZ(WALL, LOC) = REAL(WALL - 1)*ZFIRST/PIN
1300  CONTINUE

      DO 1350 CELL = 1, PIN
         WALL = CELL + 1
         Z(CELL, LOC) = (WZ(WALL, LOC) + WZ(WALL-1, LOC))/2.
1350  CONTINUE

C///  CELL WALLS AND POINTS FOR THE EXIT REGION

      ERROR = .NOT. (2 .LE. PEX)
      IF (ERROR) GO TO 9403

      LOC = 2
      ZWALLS = PEX + 1

      DO 1400 WALL = 1, ZWALLS
         WZ(WALL, LOC) = REAL(WALL - 1)*(ZLEN - ZLAST)/PEX
1400  CONTINUE

      DO 1450 CELL = 1, PEX
         WALL = CELL + 1
         Z(CELL, LOC) = (WZ(WALL, LOC) + WZ(WALL-1, LOC))/2.
1450  CONTINUE

      ERROR = .NOT. (2 .LE. NUMWF)
      IF (ERROR) GO TO 9404

C///  Z GRID AND CELL WALL IN THE WAFER LOAD

      WZWAF(1) = ZFIRST
      DO 1500 J = 2, NUMWF
         WZWAF(J) = WZWAF(J-1) + (ZLAST - ZFIRST)/REAL(NUMWF - 1)
1500  CONTINUE

      DO 1550 J = 1, NUMWF-1
         ZWAF(J) = 0.5*(WZWAF(J+1) + WZWAF(J))
1550  CONTINUE

      TEND(1) = TTEND
      TEND(2) = TTEND

C///  DEFINE AREAS FOR QUARTZ WALLS AND HEATER

      A1M = 2.*PI*RINQ(1)
      A1P = 2.*PI*ROUTQ(1)
      A2M = 2.*PI*RINQ(2)
      A2P = 2.*PI*ROUTQ(2)
      AH  = 2.*PI*RTUBEI
      AWAF = 2.*PI*WR(RPNTS+1)

      AREA(1) = A1M
      AREA(2) = A2M
      AREAW(1) = AWAF
      AREAW(2) = A2M

C///  ZERO AIN AND AEX MATRICES

      DO 1600 I = 1, 5*PIN+2
         DO 1600 J = 1, 5*PIN+2
            AIN(I, J) = 0.0
1600  CONTINUE

      DO 1650 I = 1, 5*PEX+2
         DO 1650 J = 1, 5*PEX+2
            AEX(I, J) = 0.0
1650  CONTINUE


C///////////////////////////////////////////////////////////////////////
C
C     (2) SHAPE FACTORS.
C
C///////////////////////////////////////////////////////////////////////

C///  CALCULATE SHAPEFACTORS FOR QUARTZ WALLS AND HEATER

C      FHH = 1.0 - A2P/AH
C      FH2 = A2P/AH
C      F2H = 1.0
C      F22 = 1.0 - A1P/A2M
C      F21E = 0.5*A1P/A2M
C      F21W = 0.5*A1P/A2M
C      F21 = A1P/A2M
C      F1E2 = 1.0
C      F1W2 = 1.0
C      F12 = 1.0


C///  SHAPE FACTORS FOR ENTRANCE AND EXIT REGION

      FHH = 1.0 - AREA(2)/AH
      FH2 = AREA(2)/AH
      F2H = 1.0
      F22 = 1.0 - AREA(1)/AREA(2)
      F21 = AREA(1)/AREA(2)
      F12 = 1.0

C///  MODIFICATIONS FOR SHAPE FACTORS OVER WAFER REGION

      F12W = 1.0
      F21W = AREAW(1)/AREAW(2)
      F22W = 1.0 - AREAW(1)/AREAW(2)

C*****USE SAME WALL SHAPE FACTORS OVER WAFERS AND ENTRANCE/EXIT
C
C      F12W = F12
C      F21W = F21
C      F22W = F22
C
C*****END USE SAME WALL SHAPE FACTORS OVER WAFERS AND ENTRANCE/EXIT

      IF (SHAPE) THEN

      REFLEC = 1.0 - EMISS

C///  SHAPE FACTORS BETWEEN TWO WAFERS AND SIDE WALL.

C     SUBROUTINE CONWAF
C    +  (ERROR, TEXT,
C    +   FSS, FSW, FWS, FWW, REFLEC, RPNTS, RWALL, SPACE)

      CALL CONWAF
     +  (ERROR, TEXT,
     +   FSS, FSW, FWS, FWW, REFLEC, RPNTS, WR, SPACE)

C///  ENCLOSURE CHECK FOR SHAPE FACTORS BETWEEN WAFERS AND SIDE WALL.

C*****PRINT SHAPE FACTOR CHECK FOR WAFER TO WAFER AND SIDEWALL
C
C      SUM2 = 0.0
C
C      DO 5400 I = 1, 2*RPNTS
C         SUM1 = 0.0
C         DO 5300 J = 1, 2*RPNTS
C            SUM1 = SUM1 + (1.0 - REFLEC)*FWW(I,J)
C5300     CONTINUE
C         IF (I .LE. RPNTS) THEN
C             AREA1 = PI*(WR(I+1)**2 - WR(I)**2)
C         ELSE
C             AREA1 = PI*(WR(I+1-RPNTS)**2 - WR(I-RPNTS)**2)
C         ENDIF
C         SUM1 = 2.0*PI*R(RPNTS)*SPACE/AREA1*FSW(I) + SUM1
C         SUM2 = SUM2 + (1.0 - REFLEC)*FSW(I)
C         WRITE(TEXT, 5004) I, SUM1
C5400  CONTINUE
C         SUM2 = SUM2 + FSS
C      WRITE (TEXT, 5005) SUM2
C
C5004  FORMAT(/10X,'CLOSURE FOR RING ',I2,' = ',E14.7)
C
C5005  FORMAT(/10X,'CLOSURE FOR SIDE WALL = ',E14.7)
C
C*****END PRINT SHAPE FACTOR CHECK FOR WAFER TO WAFER AND SIDEWALL

C///  SHAPE FACTORS BETWEEN ENTRANCE REGION AND FIRST WAFER.

      IOPT = 1

C     SUBROUTINE CONEND
C    +  (ERROR, TEXT,
C    +   FEE, FEG, FESE, FEW, FGE, FGSE, FSEE, FSEG, FSES, FSEW, FWE,
C    +   FWSE, IOPT, PDIM, REFLEC, RTUBEI, RPNTS, RWALL, ZPNTS, ZWALL)

C     CHANGE INNER RADIUS OF WALL FROM RTUBEI TO RINQ(1) FOR ADDITION
C     OF QUARTZ WALLS

      CALL CONEND
     +  (ERROR, TEXT,
     +   FEE(1), FEG(1), FESE(1, 1), FEW(1, 1), FGE(1), FGSE(1, 1),
     +   FSEE(1, 1), FSEG(1, 1), FSES(1, 1, 1), FSEW(1, 1, 1),
     +   FWE(1, 1), FWSE(1, 1, 1), IOPT, PDIM, REFLEC,
     +   RINQ(1), RPNTS, WR, PIN, WZ(1, 1))

C///  SHAPE FACTORS BETWEEN EXIT REGION AND LAST WAFER.

C     SUBROUTINE CONEND
C    +  (ERROR, TEXT,
C    +   FEE, FEG, FESE, FEW, FGE, FGSE, FSEE, FSEG, FSES, FSEW, FWE,
C    +   FWSE, IOPT, PDIM, REFLEC, RTUBEI, RPNTS, RWALL, ZPNTS, ZWALL)

      CALL CONEND
     +  (ERROR, TEXT,
     +   FEE(2), FEG(2), FESE(1, 2), FEW(1, 2), FGE(2), FGSE(1, 2),
     +   FSEE(1, 2), FSEG(1, 2), FSES(1, 1, 2), FSEW(1, 1, 2),
     +   FWE(1, 2), FWSE(1, 1, 2), IOPT, PDIM, REFLEC,
     +   RINQ(1), RPNTS, WR, PEX, WZ(1, 2))

C///  TURN OFF SHAPE FACTOR CALCULATIONS FOR NEXT EVALUATION.

      SHAPE = .FALSE.

      ELSE
      ENDIF

C///////////////////////////////////////////////////////////////////////
C
C     (3) DEFINE [A1] MATRIX AND FACTOR IT TO FIND RADIOSITES FOR THE
C         ENTRANCE AND EXIT SECTION.
C
C///////////////////////////////////////////////////////////////////////

C      OLD A1 MATRIX AT EVERY GRID POINT.

C      A1(1,1) = 1.0 - (1.0 - EMISSW)*FHH
C      A1(1,2) = (EMISSW - 1.0)*FH2
C      A1(1,3) = 0.0
C      A1(1,4) = 0.0
C
C      A1(2,1) = -REFL(2)*F2H
C      A1(2,2) = 1.0
C      A1(2,3) = -TRAN(2)*F22
C      A1(2,4) = -TRAN(2)*F21
C
C      A1(3,1) = -TRAN(2)*F2H
C      A1(3,2) = 0.0
C      A1(3,3) = 1.0 - REFL(2)*F22
C      A1(3,4) = -REFL(2)*F21
C
C      A1(4,1) = 0.0
C      A1(4,2) = 0.0
C      A1(4,3) = -REFL(1)*F12
C      A1(4,4) = 1.0

C /// ENTRANCE AND EXIT MATRICES

      LOC = 1
      COUNT = 0
      DO 1750 I = 1, PIN
         AIN(1 + COUNT, 1 + COUNT) = 1.0 - (1.0 - EMISSW)*FHH
         AIN(1 + COUNT, 2 + COUNT) = (EMISSW - 1.0)*FH2
         AIN(1 + COUNT, 3 + COUNT) = 0.0
         AIN(1 + COUNT, 4 + COUNT) = 0.0
         AIN(1 + COUNT, 5 + COUNT) = 0.0

         AIN(2 + COUNT,1 + COUNT) = -REFL(2)*F2H
         AIN(2 + COUNT,2 + COUNT) = 1.0
         AIN(2 + COUNT,3 + COUNT) = -TRAN(2)*F22
         AIN(2 + COUNT,4 + COUNT) = -TRAN(2)*F21
         AIN(2 + COUNT,5 + COUNT) = 0.0

         AIN(3 + COUNT,1 + COUNT) = -TRAN(2)*F2H
         AIN(3 + COUNT,2 + COUNT) = 0.0
         AIN(3 + COUNT,3 + COUNT) = 1.0 - REFL(2)*F22
         AIN(3 + COUNT,4 + COUNT) = -REFL(2)*F21
         AIN(3 + COUNT,5 + COUNT) = 0.0

         AIN(4 + COUNT,1 + COUNT) = 0.0
         AIN(4 + COUNT,2 + COUNT) = 0.0
         AIN(4 + COUNT,3 + COUNT) = -REFL(1)*F12
         AIN(4 + COUNT,4 + COUNT) = 1.0
         DO 1710 J = 1, PIN
            AIN(4 + COUNT,5*J) = -TRAN(1)*FSES(I,J,LOC)
1710     CONTINUE
         AIN(4 + COUNT,5*PIN+1) = -TRAN(1)*FSEE(I,LOC)
         AIN(4 + COUNT,5*PIN+2) = -TRAN(1)*EMISS*FSEG(I,LOC)

         AIN(5 + COUNT,1 + COUNT) = 0.0
         AIN(5 + COUNT,2 + COUNT) = 0.0
         AIN(5 + COUNT,3 + COUNT) = -TRAN(1)*F12
         AIN(5 + COUNT,4 + COUNT) = 0.0
         DO 1720 J = 1, PIN
            AIN(5 + COUNT,5*J) = KRON(I,J) - REFL(1)*FSES(I,J,LOC)
1720     CONTINUE
         AIN(5 + COUNT,5*PIN+1) = -REFL(1)*FSEE(I,LOC)
         AIN(5 + COUNT,5*PIN+2) = -REFL(1)*EMISS*FSEG(I,LOC)

         COUNT = COUNT + 5
1750  CONTINUE

C     END REGION ON INSIDE OF INNER QUARTZ TUBE

      DO 1760 J = 1, PIN
         AIN(5*PIN + 1, 5*J) = - (1.0 - EMISSE)*FESE(J,LOC)
         AIN(5*PIN + 2, 5*J) = FGSE(J,LOC)
1760  CONTINUE

      AIN(5*PIN + 1, 5*PIN + 1) = (1.0 - (1.0 - EMISSE)*FEE(LOC))
      AIN(5*PIN + 1, 5*PIN + 2) = -(1.0 - EMISSE)*EMISS*FEG(LOC)

      AIN(5*PIN + 2, 5*PIN + 1) = FGE(LOC)
C      AIN(5*PIN + 2, 5*PIN + 2) = -EMISS
      AIN(5*PIN + 2, 5*PIN + 2) = -1.0

C*****WRITE AIN MATRIX
C      DO 1765 M=1, 5*PIN+2
C         WRITE(TEXT, 99100) (AIN(M,L), L=1,5*PIN+2)
C1765  CONTINUE
C99100 FORMAT(17E9.2)
C*****END WRITE AIN MATRIX

C*****DOUBLE PRECISION

      CALL DGECO(AIN, AINDIM, AINDIM, IPVT3, RCOND1, WORK3)

C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      CALL SGECO(AIN, AINDIM, AINDIM, IPVT3, RCOND1, WORK3)
C*****END SINGLE PRECISION

      ERROR = (1.0 + RCOND1 .EQ. 1.0)
      IF (ERROR) GO TO 9003

C     -- EXIT SECTION --

      LOC = 2
      COUNT = 0
      DO 1850 I = 1, PEX
         AEX(1 + COUNT, 1 + COUNT) = 1.0 - (1.0 - EMISSW)*FHH
         AEX(1 + COUNT, 2 + COUNT) = (EMISSW - 1.0)*FH2
         AEX(1 + COUNT, 3 + COUNT) = 0.0
         AEX(1 + COUNT, 4 + COUNT) = 0.0
         AEX(1 + COUNT, 5 + COUNT) = 0.0

         AEX(2 + COUNT,1 + COUNT) = -REFL(2)*F2H
         AEX(2 + COUNT,2 + COUNT) = 1.0
         AEX(2 + COUNT,3 + COUNT) = -TRAN(2)*F22
         AEX(2 + COUNT,4 + COUNT) = -TRAN(2)*F21
         AEX(2 + COUNT,5 + COUNT) = 0.0

         AEX(3 + COUNT,1 + COUNT) = -TRAN(2)*F2H
         AEX(3 + COUNT,2 + COUNT) = 0.0
         AEX(3 + COUNT,3 + COUNT) = 1.0 - REFL(2)*F22
         AEX(3 + COUNT,4 + COUNT) = -REFL(2)*F21
         AEX(3 + COUNT,5 + COUNT) = 0.0

         AEX(4 + COUNT,1 + COUNT) = 0.0
         AEX(4 + COUNT,2 + COUNT) = 0.0
         AEX(4 + COUNT,3 + COUNT) = -REFL(1)*F12
         AEX(4 + COUNT,4 + COUNT) = 1.0
         DO 1810 J = 1, PEX
            AEX(4 + COUNT,5*J) = -TRAN(1)*FSES(I,J,LOC)
1810     CONTINUE
         AEX(4 + COUNT,5*PEX+1) = -TRAN(1)*FSEE(I,LOC)
         AEX(4 + COUNT,5*PEX+2) = -TRAN(1)*EMISS*FSEG(I,LOC)

         AEX(5 + COUNT,1 + COUNT) = 0.0
         AEX(5 + COUNT,2 + COUNT) = 0.0
         AEX(5 + COUNT,3 + COUNT) = -TRAN(1)*F12
         AEX(5 + COUNT,4 + COUNT) = 0.0
         DO 1820 J = 1, PEX
            AEX(5 + COUNT,5*J) = KRON(I,J) - REFL(1)*FSES(I,J,LOC)
1820     CONTINUE
         AEX(5 + COUNT,5*PEX+1) = -REFL(1)*FSEE(I,LOC)
         AEX(5 + COUNT,5*PEX+2) = -REFL(1)*EMISS*FSEG(I,LOC)

         COUNT = COUNT + 5
1850  CONTINUE

C     END REGION ON INSIDE OF INNER QUARTZ TUBE

      DO 1860 J = 1, PEX
         AEX(5*PEX + 1, 5*J) = - (1.0 - EMISSE)*FESE(J,LOC)
         AEX(5*PEX + 2, 5*J) = FGSE(J,LOC)
1860  CONTINUE

      AEX(5*PEX + 1, 5*PEX + 1) = (1.0 - (1.0 - EMISSE)*FEE(LOC))
      AEX(5*PEX + 1, 5*PEX + 2) = -(1.0 - EMISSE)*EMISS*FEG(LOC)

      AEX(5*PEX + 2, 5*PEX + 1) = FGE(LOC)
C      AEX(5*PEX + 2, 5*PEX + 2) = -EMISS
      AEX(5*PEX + 2, 5*PEX + 2) = -1.0

C*****WRITE AEX MATRIX
C      DO 1865 M=1, 5*PEX+2
C         WRITE(TEXT, 99200) (AEX(M,L), L=1,5*PEX+2)
C1865  CONTINUE
C99200 FORMAT(17E9.2)
C*****END WRITE AEX MATRIX

C*****DOUBLE PRECISION

      CALL DGECO(AEX, AEXDIM, AEXDIM, IPVT4, RCOND1, WORK3)

C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      CALL SGECO(AEX, AEXDIM, AEXDIM, IPVT4, RCOND1, WORK3)
C*****END SINGLE PRECISION

      ERROR = (1.0 + RCOND1 .EQ. 1.0)
      IF (ERROR) GO TO 9004

C///////////////////////////////////////////////////////////////////////
C
C     (4) DEFINE [A2] MATRIX AND FACTOR IT TO FIND RADIOSITES FOR TUBE
C         WALLS OF THE WAFER REGION.
C
C///////////////////////////////////////////////////////////////////////

      A2(1,1) = 1.0 - (1.0 - EMISSW)*FHH
      A2(1,2) = (EMISSW - 1.0)*FH2
      A2(1,3) = 0.0
      A2(1,4) = 0.0

      A2(2,1) = -REFL(2)*F2H
      A2(2,2) = 1.0
      A2(2,3) = -TRAN(2)*F22W
      A2(2,4) = -TRAN(2)*F21W

      A2(3,1) = -TRAN(2)*F2H
      A2(3,2) = 0.0
      A2(3,3) = 1.0 - REFL(2)*F22W
      A2(3,4) = -REFL(2)*F21W

      A2(4,1) = 0.0
      A2(4,2) = 0.0
      A2(4,3) = -(TRAN(1)**2*FSS*F12W/(1.-REFL(1)*FSS) + REFL(1)*F12W)
      A2(4,4) = 1.0

C*****DOUBLE PRECISION

      CALL DGECO(A2, A2DIM, A2DIM, IPVT2, RCOND2, WORK2)

C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      CALL SGECO(A2, A2DIM, A2DIM, IPVT2, RCOND2, WORK2)
C*****END SINGLE PRECISION

      ERROR = (1.0 + RCOND2 .EQ. 1.0)
      IF (ERROR) GO TO 9002

C///////////////////////////////////////////////////////////////////////
C
C     (5) ESTIMATE JWALL(I) AND JEND FOR THE ENTRANCE AND EXIT REGIONS.
C
C///////////////////////////////////////////////////////////////////////

      DO 2050 I = 1, PIN
         JWALL(I, 1) = 0.0
2050  CONTINUE
      JEND(1) = 0.0

      DO 2060 I = 1, PEX
         JWALL(PEX - I + 1, 2) = 0.0
2060  CONTINUE
      JEND(2) = 0.0

C///  STORE ESTIMATES IN S AND S0 ARRAYS

      COUNT = 0
      DO 2070 I = 1, PIN
         DO 2065 K = 1, 2
            COUNT = COUNT + 1
            S(COUNT) = TTQ
            S0(COUNT) = TTQ
2065     CONTINUE
         COUNT = COUNT + 1
         S(COUNT)  = TTUBE
         S0(COUNT) = TTUBE
         COUNT = COUNT + 1
         S(COUNT)  = TTIN
         S0(COUNT) = TTIN
         IF (CFLOW(1) .OR. CFLOW(2)) THEN
            COUNT = COUNT + 1
            S(COUNT)  = TTGAS(1)
            S0(COUNT) = TTGAS(1)
            COUNT = COUNT + 1
            S(COUNT)  = TTGAS(2)
            S0(COUNT) = TTGAS(2)
         ELSE
         ENDIF
2070  CONTINUE
      COUNT = COUNT + 1
      S(COUNT) = TTEND
      S0(COUNT) = TTEND
      COUNT = COUNT + 1
      S(COUNT) = TTIN
      S0(COUNT) = TTIN
      DO 2080 J = 1, NUMWF
         DO 2075 K = 1, RPNTS
            COUNT = COUNT + 1
            S(COUNT)  = TTWAF
            S0(COUNT) = TTWAF
2075     CONTINUE
         DO 2080 I = 1, NUMEX
            IF (I .EQ. 1) THEN
               COUNT = COUNT + 1
               S(COUNT)  = TTQ
               S0(COUNT) = TTQ
            ELSE IF (I .EQ. 2) THEN
               COUNT = COUNT + 1
               S(COUNT)  = TTQ
               S0(COUNT) = TTQ
            ELSE IF (I .EQ. 3) THEN
               COUNT = COUNT + 1
               S(COUNT)  = TTUBE
               S0(COUNT) = TTUBE
            ELSE IF (I .EQ. 4) THEN
               COUNT = COUNT + 1
               S(COUNT)  = TTIN
               S0(COUNT) = TTIN
            ELSE IF (I .EQ. 5) THEN
               COUNT = COUNT + 1
               S(COUNT)  = TTGAS(1)
               S0(COUNT) = TTGAS(1)
            ELSE IF (I .EQ. 6) THEN
               COUNT = COUNT + 1
               S(COUNT)  = TTGAS(2)
               S0(COUNT) = TTGAS(2)
            ENDIF
2080  CONTINUE
      DO 2090 I = 1, PEX
         DO 2085 K = 1, 2
            COUNT = COUNT + 1
            S(COUNT) = TTQ
            S0(COUNT) = TTQ
2085     CONTINUE
         COUNT = COUNT + 1
         S(COUNT) = TTUBE
         S0(COUNT) = TTUBE
         COUNT = COUNT + 1
         S(COUNT) = TTIN
         S0(COUNT) = TTIN
         IF (CFLOW(1) .OR. CFLOW(2)) THEN
            COUNT = COUNT + 1
            S(COUNT)  = TTGAS(1)
            S0(COUNT) = TTGAS(1)
            COUNT = COUNT + 1
            S(COUNT)  = TTGAS(2)
            S0(COUNT) = TTGAS(2)
         ELSE
         ENDIF
2090  CONTINUE
      COUNT = COUNT + 1
      S(COUNT) = TTEND
      S0(COUNT) = TTEND
      COUNT = COUNT + 1
      S(COUNT) = TTIN
      S0(COUNT) = TTIN

C///////////////////////////////////////////////////////////////////////
C
C     (6) INTERPOLATE HEATER TEMPERATURES.
C
C///////////////////////////////////////////////////////////////////////

      IF (.NOT. TFIXED) GO TO 99999

C///  FIND THE MAX. AND MIN INDICES OF GRID POINTS BETWEEN ZDAT(1) AND
C     ZDAT(NUMT).

      OFFSET(1) = 0.0
      OFFSET(2) = ZLEN
      IBEG(1) = PIN
      IBEG(2) = PEX
      IBEG(3) = NUMWF
      IEND(1) = 0
      IEND(2) = 0
      IEND(3) = 0

      DO 3030 LOC = 1, 2
         IF (LOC .EQ. 1) CELLS = PIN
         IF (LOC. EQ. 2) CELLS = PEX
         DO 3010 CELL = 1, CELLS
            ZTEST = ABS(OFFSET(LOC) - Z(CELL, LOC))
            IF ((ZDAT(1).LE.ZTEST) .AND. (ZTEST.LE.ZDAT(NUMT))) THEN
               IBEG(LOC) = MIN0(CELL, IBEG(LOC))
               IEND(LOC) = MAX0(CELL, IEND(LOC))
            ELSE
            ENDIF
3010     CONTINUE
         DO 3020 I = IBEG(LOC), IEND(LOC)
            ZTEST = ABS(OFFSET(LOC) - Z(I, LOC))
            CALL INTERP (NUMT, ZTEST, ZDAT, TDAT, TSIDE(I, LOC))
3020     CONTINUE
3030  CONTINUE

      DO 3040 CELL = 1, NUMWF - 1
         ZTEST = ZWAF(CELL)
         IF ((ZDAT(1).LE.ZTEST) .AND. (ZTEST.LE.ZDAT(NUMT))) THEN
            IBEG(3) = MIN0(CELL, IBEG(3))
            IEND(3) = MAX0(CELL, IEND(3))
         ELSE
         ENDIF
3040  CONTINUE
      DO 3050 I = IBEG(3), IEND(3)
         CALL INTERP (NUMT, ZWAF(I), ZDAT, TDAT, TMID(I))
3050  CONTINUE


C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, RCOND1
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, RCOND2
      GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, RCOND1
      GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, RCOND1
      GO TO 99999

9402  IF (0 .LT. TEXT) WRITE (TEXT, 99402) ID, PIN
      GO TO 99999

9403  IF (0 .LT. TEXT) WRITE (TEXT, 99403) ID, PEX
      GO TO 99999

9404  IF (0 .LT. TEXT) WRITE (TEXT, 99404) ID, NUMWF
      GO TO 99999

9405  IF (0 .LT. TEXT) WRITE (TEXT, 99405) ID, RPNTS
      GO TO 99999

99001 FORMAT
     +   (/1X, A9, 'ERROR.  THE A1 MATRIX FOR WALL RADIOSITES OVER'
     +   /10X, 'THE ENTRANCE AND EXIT REGION IS SINGULAR TO MACHINE'
     +   /10X, 'ACCURACY.'
     +  //10X, 'A1 MATRIX CONDITION NUMBER = ', E14.7)

99002 FORMAT
     +   (/1X, A9, 'ERROR.  THE A2 MATRIX FOR WALL RADIOSITES OVER'
     +   /10X, 'THE WAFER REGION IS SINGULAR TO MACHINE ACCURACY.'
     +  //10X, 'A2 MATRIX CONDITION NUMBER = ', E14.7)

99003 FORMAT
     +   (/1X, A9, 'ERROR.  THE AIN MATRIX FOR WALL RADIOSITES OVER'
     +   /10X, 'THE ENTRANCE REGION IS SINGULAR TO MACHINE'
     +   /10X, 'ACCURACY.'
     +  //10X, 'AIN MATRIX CONDITION NUMBER = ', E14.7)

99004 FORMAT
     +   (/1X, A9, 'ERROR.  THE AEX MATRIX FOR WALL RADIOSITES OVER'
     +   /10X, 'THE EXIT REGION IS SINGULAR TO MACHINE'
     +   /10X, 'ACCURACY.'
     +  //10X, 'AEX MATRIX CONDITION NUMBER = ', E14.7)

99402 FORMAT
     +   (/1X, A9, 'ERROR.  THE ENTRANCE GRID MUST HAVE AT LEAST'
     +   /10X, 'TWO POINTS.'
     +  //10X, I10, '  POINTS')

99403 FORMAT
     +   (/1X, A9, 'ERROR.  THE EXIT GRID MUST HAVE AT LEAST'
     +   /10X, 'TWO POINTS.'
     +  //10X, I10, '  POINTS')

99404 FORMAT
     +   (/1X, A9, 'ERROR.  THE WAFER LOAD Z GRID MUST HAVE AT LEAST'
     +   /10X, 'TWO POINTS.'
     +  //10X, I10, '  POINTS')

99405 FORMAT
     +   (/1X, A9, 'ERROR.  THE WAFER LOAD R GRID MUST HAVE AT LEAST'
     +   /10X, 'TWO POINTS.'
     +  //10X, I10, '  POINTS')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWAF1
     +  (ERROR, TEXT,
     +   A2, A2DIM, AEX, AIN, AREAC, BTEMP, CFLOW, CHOICE,
     +   COMPS, CONDUC,
     +   CONW, CP, DZEND, DZENDI, EMIS, EMISS, EMISSE, EMISSI, EMISSW,
     +   F, FEE, FEG, FENDS, FENG, FESE, FEW, FGE, FGSE, FSEE, FSEG,
     +   FSES, FSEW, FSS,
     +   FSW, FWE, FWS, FWSE, FWW, GROUPA, GROUPB, HCOEF, IBEG, IEND,
     +   IPVT2, IPVT3, IPVT4, J1, J2, JEND, JEX, JIN, JWALL, LOS1,
     +   LOS2, LOS3, MAXHT, MAXIN, MFLOW, MINJ, MINJE, MWZ, MWZW,
     +   NUMEX, NUMHT, NUMIN, NUMT, NUMWF, PCAP, PDIM, PEX, PIN,
     +   PMAX, POWER, POWERH, PWALL, QWAF, R, RBOAT, RCPV, REFL,
     +   RHO, RING, RINJEC, RINQ, ROUTQ,
     +   RPNTS, RTUBEI, RTUBEO,
     +   S, S0, SHAPE, SPACE, TCAP, TCAPS, TDAT, TEAST, TEMMAX, TEX,
     +   TEX0, TEX4, TFIX, THICK, THICKI, TIME, TMID, TRAN, TSIDE,
     +   TSTEP, TTGAS, TWAF, TWAF0, TWAF4, TWEST, TX, TX0, TX4,
     +   VINJ, WR, WZ, WZWAF,
     +   Z, ZBEG, ZBEGI, ZDAT, ZEND, ZENDI, ZLEN, ZWAF,
     +   GBOT, GTOP, JBOT, JTOP, NUMTC, SUMB, SUMT, TC, ZTC)

C///////////////////////////////////////////////////////////////////////
C
C     TWAF1
C
C///  EVALUATE THE RESIDUAL OF THE WAFER THERMAL ENERGY EQUATION.
C
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     INPUT:
C
C     CFLOW   INPUT LOGICAL DIMENSIONED (2) - CONTROLS GAS COOLING.
C             CFLOW(1) = .TRUE. - INCLUDE GAS COOLING BETWEEN QUARTZ
C                                 TUBE.
C             CFLOW(2) = .TRUE. - INCLUDE GAS COOLING BETWEEN OUTER
C                                 QUARTZ TUBE AND HEATER.
C
C     CHOICE  INPUT INTEGER - UNPACKING VARIABLE.
C             CHOICE = 0 - NORMAL MODE FOR TWOPNT.
C             CHOICE = 1,2,3 - MODES USED TO OUTPUT MATRIX FOR MATRIXX
C                              DURING TWOPNT MODE.
C             CHOICE = 4 - MODE USED FOR DASSL, ASSUMES S0 ARRAY
C                          CONTAINS YPRIME FROM DASSL.  IF
C                          TIME = .TRUE. THEN S0 (DASSL YP) IS ADDED
C                          TO THE RESIDUAL.
C
C     CONDUC  INPUT REAL DIMENSION (8) - THERMAL CONDUCTIVITY ARRAY.
C             CONDUC(1) - WAFER MATERIAL.
C             CONDUC(2) - HEATERS.
C             CONDUC(3) - INSULATION
C             CONDUC(4) - QUARTZ.
C             CONDUC(5) - ENTRANCE CAP.
C             CONDUC(6) - EXIT CAP.
C             CONDUC(7) - ENTRANCE CAP INSULATION.
C             CONDUC(8) - EXIT CAP INSULATION.
C
C     CP      INPUT REAL DIMENSION (8) - SPECIFIC HEAT ARRAY.
C             CP(1) - WAFER MATERIAL.
C             CP(2) - HEATER.
C             CP(3) - INSULATION
C             CP(4) - QUARTZ.
C             CP(5) - ENTRANCE CAP.
C             CP(6) - EXIT CAP.
C             CP(7) - ENTRANCE CAP INSULATION.
C             CP(8) - EXIT CAP INSULATION.
C
C     CPAIR   INPUT REAL - SPECIFIC HEAT OF AIR FOR COOLING.
C
C     DZEND   INPUT REAL DIMENSIONED 2 - DZEND(1) IS THE THICKNESS OF
C             OF CAP ON THE REACTOR ENTRANCE END.  DZEND(2) IS THE
C             THICKNESS OF THE CAP ON THE EXIT.
C
C     DZENDI  INPUT REAL DIMENSIONED 2 - DZENDI(1) IS THE THICKNESS OF
C             INSULATION ON THE ENTRANCE CAP.  DZENDI(2) IS THE
C             THICKNESS OF INSULATION OF THE EXIT CAP.
C
C     EMIS    INPUT REAL DIMENSIONED (2) - EMIS(1) AND EMIS(2) ARE
C             THE EMITTANCES OF QUARTZ TUBES 1 AND 2.
C
C     EMISS   INPUT REAL - EMISSIVITY OF WAFERS.
C
C     EMISSE  INPUT REAL - EMISSIVITY OF REACTOR END WALL.
C
C     EMISSI  INPUT REAL - EMISSIVITY OF THE REACTOR OUTER SKIN.
C             EMISSI(1) - OUTER SKIN OF ENTRANCE CAP.
C             EMISSI(2) - OUTER SKIN OF EXIT CAP.
C             EMISSI(3) - OUTER SKIN OF MAIN REACTOR BODY (OVER WAFERS).
C             EMISSI(4) - OUTER SKIN OF REACTOR ENTRANCE SECTION.
C             EMISSI(5) - OUTER SKIN OF REACTOR EXIT SECTION.
C
C
C     EMISSW  INPUT REAL - EMISSIVITY OF REACTOR HEATERS.
C
C     HCOEF   INPUT REAL DIMENSION (3) - EXTERNAL CONVECTIVE HEAT
C             TRANSFER COEFFICIENTS (W/cm**2/K).
C             HCOEF(1) - OUTER SKIN OF ENTRANCE CAP.
C             HCOEF(2) - OUTER SKIN OF EXIT CAP.
C             HCOEF(3) - OUTER SKIN OF MAIN REACTOR BODY (OVER WAFERS).
C             HCOEF(4) - OUTER SKIN OF ENTRANCE SECTION REACTOR BODY.
C             HCOEF(5) - OUTER SKIN OF EXIT SECTION REACTOR BODY.
C
C     MAXHT   INPUT INTEGER - MAXIMUM NUMBER OF HEATER ELEMENTS.
C
C     MFLOW   INPUT REAL DIMENSION (2) - COOLING GAS MASS FLOW RATE
C             (G/SEC).
C
C     NUMHT   INPUT INTEGER - NUMBER OF HEATER ELEMENTS.
C
C     NUMT    INPUT INTERER - NUMBER OF TEMPERATURE DATA POINTS ALONG
C             WALL.
C
C     NUMTC   INPUT INTEGER - NUMBER OF THERMOCOUPLES.
C
C     NUMWF   INPUT INTEGER - NUMBER OF SIMULATED WAFERS.
C
C     NUMEX   INPUT INTEGER - NUMBER OF EXTERNAL RADIAL GRID POINTS.
C
C     PCAP    INPUT REAL DIMENSION (2) - INPUT POWER FOR ENTRANCE AND
C             EXIT CAPS.
C
C     PDIM    INPUT INTERGER - MAXIMUM DIMENSION OF SOME ARRAYS.
C
C     PIN     INPUT INTEGER - NUMBER OF POINTS BETWEEN REACTOR END WALL
C             AND FIRST WAFER IN THE ENTRANCE REGION.
C
C     PEX     INPUT INTEGER - NUMBER OF POINTS BEWTEEN REACTOR END WALL
C             AND LAST WAFER IN THE EXIT REGION.
C
C     POWER   INPUT REAL DIMENSION (NUMWF) - WALL POWERS AT GRID POINTS
C             ABOVE THE WAFERS. POWER(I) IS THE POWER (WATTS)
C             DISSIPATED IN THE WALL ELEMENT AT WAFER I.
C
C     POWERH  INPUT REAL DIMENSION (MAXHT) - POWER(I) IS THE POWER IN
C             WATTS TO HEATER I.
C
C     PMAX    INPUT INTEGER - NUMBER OF POSITIONS AT WHICH SOLUTION IS
C             COMPUTED.
C
C     PWALL   INPUT REAL DIMENSION (PDIM, 2) - PWALL(I, 1) IS THE POWER
C             (WATTS) DISSIPATED IN WALL ELEMENT "I" IN THE ENTRANCE.
C
C     R       INPUT REAL DIMENSIONED (RPNTS) - RADIAL GRID POINTS IN
C             WAFERS.
C
C     RBOAT   INPUT LOGICAL - IF .TRUE. INCLUDE RING BOAT IN
C             CALCULATION.
C
C     REFL    INPUT REAL DIMENSIONED (2) - REFL(1) AND REFL(2) ARE THE
C             REFLECTANCES OF QUARTZ TUBES 1 AND 2.
C
C     RHO     INPUT REAL DIMENSIONED (8) - DENSITY ARRAY.
C             RHO(1) - WAFER MATERIAL.
C             RHO(2) - HEATERS.
C             RHO(3) - DENSITY INSULATION.
C             RHO(4) - DENSITY OF QUARTZ.
C             RHO(5) - ENTRANCE CAP MATERIAL.
C             RHO(6) - EXIT CAP MATERIAL.
C             RHO(7) - INSULATION ON ENTRANCE CAP.
C             RHO(8) - INSULATION ON EXIT CAP.
C
C     RING    INPUT REAL DIMENSIONED (5) - RING BOAT PARAMETERS.
C             RING(1) - RING BOAT INSIDE RADIUS.
C             RING(2) - RING BOAT THICKNESS.
C             RING(3) - DENISITY OF RING BOAT MATERIAL.
C             RING(4) - SPECIFIC HEAT OF RING BOAT MATERIAL.
C             RING(5) - THERMAL CONDUCTIVITY OF RING BOAT MATERIAL.
C
C     RINJEC  INPUT LOGICAL - IF .TRUE. INCLUDE RADIAL INJECTED MASS
C             FLOWS IN COOLING GAS ENERGY BALANCE.  RADIAL INJECTION
C             ONLY ALLOWED FOR OUTER FLOW.
C
C     RINQ    INPUT REAL DIMENSIONED (2) - RINQ(1) AND RINQ(2) ARE THE
C             INNER RADII OF QUARTZ TUBE 1 AND 2, RESPECTIVELY.
C
C     ROUTQ   INPUT REAL DIMENSIONED (2) - ROUTQ(1) AND ROUTQ(2) ARE
C             THE OUTER RADII OF QUARTZ TUBE 1 AND 2, RESPECTIVELY.
C
C     RTUBEI  INPUT REAL - INNER RADIUS OF THE HEATER.
C
C     RTUBEO  INPUT REAL - OUTER RADIUS OF THE HEATER.
C
C     RPNTS   INPUT INTEGER - NUMBER OF RADIAL POINTS IN A WAFER.
C
C     SHAPE   INPUT LOGICAL - IF .TRUE. PERFORM A SHAPE FACTOR
C             CALCULATION BEFORE COMPUTING THE RESIDUAL.
C
C     SPACE   INPUT REAL - SPACING BETWEEN WAFERS.
C
C     TCAPS   INPUT LOCICAL DIMENSIONED 2
C             TCAPS(1) = .TRUE. USE TCAP(1) AS THE ENTRANCE CAP TEMP.
C             TCAPS(2) = .TRUE. USE TCAP(2) AS THE EXIT CAP TEMP.
C
C     TDAT    INPUT REAL DIMENSIONED (TEMMAX + 2) - TEMPERATURE DATA
C             POINTS ALONG REACTOR WALL AT ZDAT(I) LOCATIONS.
C
C     TIME    INPUT LOGICAL - IF TIME = .TRUE. THEN TRANSIENT TERM
C             IS ADDED TO RESIDUAL.
C
C     TRAN    REAL DIMENSIONED (2) - TRAN(1) AND TRAN(2) ARE THE
C             TRANSMITTANCES OF QUARTZ TUBES 1 AND 2.
C
C     TWALL   INPUT REAL DIMENSIONED (PDIM, 2) - TWALL(I, 1) IS THE
C             REACTOR WALL TEMPERATURE AT Z(I, 1) IN THE ENTRANCE
C             REGION WALL.
C
C     THICK   INPUT REAL - THICKNESS OF THE WAFER.
C
C     THICKI  INPUT REAL - THICKNESS OF THE INSULATION.
C
C     WZ      REAL DIMENSIONED (PDIM+1, 2) - WZ(I, 1) ARE THE CELL WALL
C             LOCATIONS BETWEEN REACTOR END WALL AND FIRST WAFER.
C
C     WZWAF   REAL DIMENSIONED (NUMWF+1) - WZWAF(I) ARE THE CELL WALL
C             LOCATIONS IN THE WAFER LOAD.
C
C     Z       INPUT REAL DIMENSIONED (PDIM, 2) - AXIAL GRID POINTS
C             BETWEEN THE REACTOR END WALLS AND FIRST WAFER OR LAST
C             WAFERS.
C
C     ZBEG    INPUT REAL DIMENSIONED (MAXHT) - ZBEG(I) IS BEGINNING
C             Z COORDINATE FOR HEATER "I".
C
C     ZDAT    INPUT REAL DIMENSIONED (TEMMAX + 2) - AXIAL LOCATIONS OF
C             TEMPERATURE DATA POINTS ALONG REACTOR WALL.
C
C     ZEND    INPUT REAL DIMENSIONED (MAXHT) - ZEND(I) IS BEGINNING
C             Z COORDINATE FOR HEATER "I".
C
C     ZLEN    INPUT REAL - LENGTH OF THE REACTOR (CM).
C
C     ZTC     REAL DIMENSIONED (NUMTC) - AXIAL LOCATION OF
C             THERMOCOUPLES.
C
C     ZWAF    INPUT REAL DIMEMSIONED (NUMWF)  - ZWAF(I) IS THE AXIAL
C             LOCATION OF SIMULATED WAFER I.
C
C     INTERNAL:
C
C     AEX     REAL DIMENSIONED (5*PEX+2, 5*PEX+2) - AEX IS THE
C             FACTORED MATRIX FOR RADIOSITIES IN THE EXIT.
C
C     AIN     REAL DIMENSIONED (5*PIN+2, 5*PIN+2) - AIN IS THE
C             FACTORED MATRIX FOR RADIOSITIES IN THE ENTRANCE.
C
C     AREAC   REAL DIMENSIONED (RPNTS+1) - CROSS SECTIONAL AREA AT
C             CELL FACES IN A WAFER.
C
C     A2      REAL DIMENSIONED (A2DIM, A2DIM) - A2 IS THE FACTORED
C             MATRIX FOR RADIOSITIES OVER THE WAFER REGION.
C
C     CONVN   REAL DIMENSIONED (2) - CONVECTIVE AIR COOLING TERMS.
C             CONVN(FLOW) IS THE NORTHERN CONVECTIVE ENERGY TERM
C             FOR THE AIR FOR FLOW = 1 OR 2.
C
C     CONVS   REAL DIMENSIONED (2) - CONVECTIVE AIR COOLING TERMS.
C             CONVN(FLOW) IS THE SOUTHERM CONVECTIVE ENERGY TERM
C             FOR THE AIR FOR FLOW = 1 OR 2.
C
C     CONW    CONW(RPNTS+1) - CONDUCTIVITY AT WAFER CELL INTERFACES.
C
C     EGAP    REAL DIMENSIONED (2) - EGAP(1) IS THE BLACKBODY EMISSIVE
C             POWER FROM THE GAP REGION BETWEEN THE FIRST WAFER AND THE
C             REACTOR WALL.
C
C     FEE     REAL DIMENSIONED (2) - SHAPE FACTOR FROM REACTOR END WALL
C             TO ITSELF.
C
C     FEG     REAL DIMENSIONED (2) - FEG(1) IS THE SHAPE FACTOR FROM
C             REACTOR END WALL TO THE GAP BETWEEN THE FIRST WAFER AND
C             THE WALL.
C
C     FESE    REAL DIMENSIONED (PDIM, 2) - FESE(I) ARE THE SHAPE
C             FACTORS FROM THE REACTOR END WALL TO ELEMENTS I ON THE
C             WALL OF THE REACTOR ENTRANCE.
C
C     FEW     REAL DIMENSIONED (RPNTS, 2) - FEW(I) ARE THE SHAPE
C             FACTORS FROM THE THE REACTOR END WALL TO RING ELEMENTS I
C             ON THE FIRST WAFER.
C
C     FGE     REAL DIMENSIONED (2) - FGE(1) IS THE SHAPE FACTOR FROM
C             THE GAP BETWEEN THE FIRST WAFER AND THE WALL TO THE END
C             WALL OF THE REACTOR.
C
C     FGSE    REAL DIMENSIONED (PDIM, 2) - FGSE(I, 1) ARE THE SHAPE
C             FACTORS FROM THE GAP BETWEEN THE FIRST WAFER AND THE WALL
C             TO THE ELEMENTS I ON THE SIDE WALL IN THE ENTRANCE REGION.
C
C     FSEE    REAL DIMENSIONED (PDIM, 2) - SHAPE FACTORS FROM
C             ELEMENTS ON REACTOR SIDE WALL TO REACTOR END WALL.
C
C     FSEG    REAL DIMENSIONED (PDIM, 2) - FSEG(I, 1) ARE THE SHAPE
C             FACTORS FROM THE ELEMENTS I ON THE ENTRANCE WALL OF THE
C             REACTOR TO THE GAP BETWEEN THE FIRST WAFER AND THE WALL.
C
C     FSES    REAL DIMENSIONED (PDIM, PDIM, 2) - SHAPE FACTORS
C             BETWEEN ELEMENTS ON REACTOR SIDE WALL IN THE ENTRANCE
C             REGION.
C
C     FSEW    REAL DIMENSIONED (PDIM, PDIM, 2) - FSEW(I,J) ARE THE
C             SHAPE FACTORS FROM ELEMENTS I ON THE WALL OF THE
C             REACTOR ENTRANCE TO RING ELEMENTS J ON THE FIRST WAFER.
C
C     FSS     REAL - CONFIGURATION FACTOR FROM THE REACTOR SIDE
C             WALL AREA BETWEEN TWO WAFERS TO ITSELF.
C
C     FSW     REAL DIMENSIONED (2*RPNTS) - FSW(I) ARE THE SHAPE
C             FACTORS FROM THE REACTOR WALL AREA BETWEEN TWO WAFERS TO
C             THE RING ELEMENTS ON THE WAFER SURFACES.
C
C     FWE     REAL DIMENSIONED (RPNTS, 2) - FWE(I, 1) ARE THE SHAPE
C             FACTORS FROM THE RING ELEMENTS I ON THE FIRST WAFER TO
C             THE REACTOR END WALL.
C
C     FWS     REAL DIMENSIONED (2*RPNTS) - FWS(I) ARE THE SHAPE
C             FACTORS FROM THE RING ELEMENTS ON THE WAFER SURFACES TO
C             THE REACTOR WALL AREA BETWEEN TWO WAFERS.
C
C     FWSE    REAL DIMENSIONED (RPNTS, PDIM, 2) - FWSE(I,J,1) ARE
C             THE SHAPE FACTORS FROM RING ELEMENTS I ON THE FIRST
C             WAFER TO ELEMENTS J ON THE WALL OF THE REACTOR ENTRANCE.
C
C     FWW     REAL DIMENSIONED (2*RPNTS, 2*RPNTS) - FWW(I,J) ARE
C             THE SHAPE FACTORS FROM WAFER RING ELEMENT I TO WAFER RING
C             ELEMENT J. ELEMENTS 1 TO RPNTS ARE ON THE LEFT WAFER FROM
C             THE CENTERLINE TO THE EDGE. ELEMENTS RPNTS+1 TO 2*RPNTS
C             ARE ON THE RIGHT WAFER FROM THE CENTERLINE TO THE EDGE.
C
C     GBOT    REAL DIMENSIONED (NUMTC) - IRRADIATION ON BOTTOM OF TC.
C
C     GTOP    REAL DIMENSIONED (NUMTC) - IRRADIATION ON TOP OF TC.
C
C     IOPT    INTEGER - CONTROLS TREATMENT OF THE GAP BETWEEN THE EDGE
C             OF THE BOUNDING WAFER AND THE REACTOR WALL IN THE ENTRANCE
C             AND EXIT REGIONS OF THE REACTOR.
C             IF IOPT = 1 THEN THE GAP REGION IS TREATED AS A
C             RERADIATING SURFACE OTHERWISE RADIATION IS LOST THROUGH
C             THE GAP.
C
C     MWZW    REAL DIMENSIONED (NUMWF+1) - MASS FLOWRATES AT CELL
C             WALLS, FOR THE FLOW CHANNEL BETWEEN THE OUTER QUARTZ
C             JAR AND THE HEATERS. COMPUTED BY SUBROUTINE MASSD.
C
C     J1      REAL DIMENSIONED (4, PDIM, 2) - WALL RADIOSITIES
C             IN ENTRANCE OR EXIT.
C             J1(1,I,LOC) - J, HEAT (AT Z INDEX I, IF LOC = 1 IN
C                                    ENTRANCE)
C             J1(2,I,LOC) - J2, PLUS
C             J1(3,I,LOC) - J2, MINUS
C             J1(4,I,LOC) - J1, PLUS
C
C     J2      REAL DIMENSIONED (A2DIM, NUMWF) - WALL RADIOSITIES OVER
C             THE WAFER REGION.
C             J2(1,WAFER) - J, HEAT
C             J2(2,WAFER) - J2, PLUS
C             J2(3,WAFER) - J2, MINUS
C             J2(4,WAFER) - J1, EAST, PLUS
C
C     JBOT    REAL DIMENSIONED (NUMTC) - RADIOSITY ON BOTTOM OF TC.
C
C     JTOP    REAL DIMENSIONED (NUMTC) - RADIOSITY ON TOP OF TC.
C
C     JWALL   REAL DIMENSIONED (PDIM, 2) - JWALL(I, 1) ARE VALUES OF
C             THE SIDE-WALL RADIOSITIES IN THE ENTRANCE
C             REGION.  JWALL(I, 2) ARE  RADIOSITIES OF THE SIDE-WALL
C             IN THE EXIT REGION.
C
C     QEDGE   REAL DIMENSIONED (4, 2) - HEAT FLUX FROM CYLINDER ENDS.
C             QEDGE(1, LOC) - QUARTZ WALL 1 (INNER MOST).
C             QEDGE(2, LOC) - QUARTZ WALL 2.
C             QEDGE(3, LOC) - HEATERS.
C             QEDGE(4, LOC) - HEATER INSULATION BLANKET.
C
C     RCPV    REAL DIMENSIONED (RPNTS) - PRODUCT OF RHO, CP, AND
C             CELL VOLUME FOR EACH CELL IN A WAFERS.
C
C     RINI    REAL - INNER RADIUS OF INSULATION ON HEATERS.
C
C     RINO    REAL - OUTER RADIUS OF INSULATION ON HEATERS.
C
C     SIG     REAL - STEFAN BOLTZMAN CONSTANT (W/CM**2/K**4).
C
C     SUMB    REAL DIMENSIONED (NUMTC)
C
C     SUMT    REAL DIMENSIONED (NUMTC)
C
C     TEAST   REAL DIMENSIONED (NUMWF) - TEAST(I) IS THE EASTERN REACTOR
C             WALL TEMPERATURE FOR WAFER I.
C
C     TTGAS   REAL DIMENSIONED (2) - INTIAL COOLING GAS
C             TEMPERATURES (K).
C
C     TWEST   REAL DIMENSIONED (NUMWF) - TWEST(I) IS THE WESTERN REACTOR
C             WALL TEMPERATURE FOR WAFER I.
C
C     WR      REAL DIMENSIONED (RPNTS+1) - CELL WALL LOCATIONS
C             IN WAFERS.
C
C     WZ      REAL DIMENSIONED (PDIM+1, 2) -  WZ(I, 1) ARE THE CELL WALL
C             LOCATIONS BETWEEN THE REACTOR END WALL AND FIRST WAFER.
C
C     OUTPUT:
C
C     F       OUTPUT REAL DIMENSIONED (GROUPA + COMPS*PMAX + GROUPB) -
C             PACKED RESIDUAL.
C
C     FENG    OUTPUT REAL DIMESIONED (NUMWF, RPNTS+NUMEX) - FENG(I, J)
C             IS THE ENERGY BALANCE FOR CELL J OF WAFER I.
C             FENG(NUMWF, RPNTS+1) = ENERGY BALANCE FOR QUARTZ WALL 1.
C             FENG(NUMWF, RPNTS+2) = ENERGY BALANCE FOR QUARTZ WALL 2.
C             FENG(NUMWF, RPNTS+3) = ENERGY BALANCE FOR HEATER.
C             FENG(NUMWF, RPNTS+4) = ENERGY BALANCE FOR INSULATION.
C
C     FENDS   OUTPUT REAL DIMENSIONED (PDIM + 1, NUMEX + 1, 2) -
C             FENDS(I, 1, 1) IS THE RADIOSITY ENERGY BALANCE FOR
C             ELEMENT I ON THE INSIDE OF THE REACTOR QUARTZ
C             ENTRANCE WALL.
C             FENDS(I, 1, 2) SAME AS FENDS(I, 1, 1) BUT FOR THE
C             EXIT QUARTZ WALL.
C             FENDS(I, 2, 1) IS THE ENERGY BALANCE FOR THE QUARTZ
C             1 TEMPERATURE.
C             FENDS(I, 3, 1) IS THE ENERGY BALANCE FOR THE QUARTZ
C             WALL 2 TEMPERATURE.
C             FENDS(I, 4, 1) IS THE ENERGY BALANCE FOR THE HEATER
C             TEMPERATURE.
C             FENDS(I, 5, 1) IS THE ENERGY BALANCE FOR THE INSULATION
C             TEMPERATURE.
C
C     LOS1    OUTPUT REAL DIMENSIONED (NUMWF) - EXTERIOR HEAT LOSS FROM
C             ELEMENTS OVER WAFERS.
C
C     LOS2    OUTPUT REAL DIMENSIONED (PDIM, 2) - EXTERIOR HEAT LOSS
C             FROM ENTRANCE AND EXIT ELEMENTS.
C
C     LOS3    OUTPUT REAL DIMENSIONED (2) - EXTERIOR HEAT LOSS FROM
C             END CAPS INSULATION (INCLUDING ENDS AND SIDES).
C
C     S       OUTPUT REAL DIMENSIONED (GROUPA + COMPS*PMAX + GROUPB) -
C             TWAF(I, J) PACKED.
C
C     TC      REAL DIMENSIONED (NUMTC) - THERMOCOUPLE TEMPERATURES
C             TO THE FOURTH POWER.
C
C     TEX     OUTPUT REAL DIMENSIONED (NUMWF, NUMEX) -
C             TEX(I, J) IS THE TEMPERATURE OF AN EXTERNAL COMPONENT AT
C             THE AXIAL LOCATION OF WAFER I AND RADIAL LOCATION J.
C             TEX(I, 1) - QUARTZ WALL 1.
C             TEX(I, 2) - QUARTZ WALL 2.
C             TEX(I, 3) - HEATER.
C             TEX(I, 4) - INSULATION.
C
C     TWAF    OUTPUT REAL DIMENSIONED (NUMWF+1, RPNTS) - TWAF(I, J) IS
C             THE TEMPERATURE OF WAFER I AT RADIAL POINT J.
C
C     TX      OUTPUT REAL DIMENSIONED (PDIM, NUMEX, 2) -
C             TX(I, J, LOC) IS THE TEMPERATURE OF AN EXTERNAL COMPONENT
C             J AT GRID LOCATION I IN THE ENTRANCE OR EXIT SECTION
C             (LOC = 1 OR 2) OF THE REACTOR.
C             TX(I, 1, LOC) - QUARTZ WALL 1.
C             TX(I, 2, LOC) - QUARTZ WALL 2.
C             TX(I, 3, LOC) - HEATER.
C             TX(I, 4, LOC) - INSULATION.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   A2DIM, CELL, CHOICE, COMPS, COUNT, FLOW, GROUPA, GROUPB, I,
     +   IBEG, IEND, IOPT, IPVT2, IPVT3, IPVT4, J, K, LOC,
     +   MAXHT, MAXIN, NUMEX, NUMHT, NUMIN, NUMT, NUMTC, NUMWF, PDIM,
     +   PEX, PIN, PMAX,
     +   RPNTS, TEMMAX, TEXT, TOTAL, WAFER, WALL
      LOGICAL CORR, ERROR, FLAG, FLAG0, FLAG00,
     +        CFLOW, FLAG4, PCON,
     +        PED, RBOAT, RINJEC, SHAPE, TCAPS, TFIX, TIME, TURB
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   A1M, A1P, A2, A2M, A2P, ABSBOT, ABSTOP, ACOND, AEX, AIN,
     +   AH, AINN, AOUT, AREA, AREAC, AREAF, AREAW, AWAF, BTEMP,
     +   COND, CONDUC, CONVN, CONVS, CONW, CP, CPAIR, DH, DUMMY,
     +   DZ, DZEND, DZENDI,
     +   EEAST, EINJ, EMIS, EMISS, EMISSE, EMISSI, EMISSW, EWEST,
     +   F, F12,
     +   F12W, F21, F21W, F22, F22W, F2H, FACTOR, FBOT, FEE, FEG,
     +   FENDS, FENG, FESE, FEW, FGE,
     +   FGSE, FH2, FHH, FSEE, FSEG, FSES, FSEW, FSS, FSW, FTOP, FWE,
     +   FWS, FWSE, FWW, GBOT, GEAST, GTOP, GWEST, HCOEF, HCOF,
     +   HINN, HOUT,
     +   J1, J2, JBOT, JEAST, JEND, JEX, JIN, JTOP, JWALL, JWEST,
     +   KAIR, LOS1, LOS2, LOS3, MFLOW, MFLOWE, MFLOWW, MINJ, MINJE,
     +   MWZ, MWZW, NU, MUAIR, PCAP,
     +   PEDRI, PEDRO,
     +   PI, POWER, POWERH, PRAIR,
     +   PWALL, QCONV, QEAST, QFLUX,
     +   QNORTH, QSOUTH, QTIP, QTIP2, QWAF, QWEST, R, RE, RCPV,
     +   REFL, REYN, RHO, RHOAIR, RIN, RING, RINI, RINO, RINQ, ROUTQ,
     +   RTC, RTUBE, RTUBEI, RTUBEO, S, S0,
     +   SAV1, SAV2, SAV3,
     +   SCALE, SCALEE, SCALEW,
     +   SIG, SPACE, SUMB, SUME, SUMT, SUMW, SUMWE, SUMWW, TAIR, TAMB,
     +   TAMB4, TC, TCAP, TDAT, TEAST, TEAST1, TEND, TEND0, TEND4,
     +   TENDI, TENDI0, TENDI4, TEX, TEX0,
     +   TEX4, TFILM, TTHIC,
     +   THICK, THICKI, TMID, TRAN, TREF, TSIDE,
     +   TSTEP, TTGAS, TWAF, TWAF0, TWAF4, TWEST, TWEST1, TX, TX0,
     +   TX4, VINJ, VOL, WR, WZ, WZWAF, Z, ZBEG, ZBEGI, ZDAT, ZEND,
     +   ZENDI, ZLEN, ZTC, ZWAF, TMASS

      PARAMETER (ID = ' TWAF1:  ')
      PARAMETER (PI = 3.141592653589793238462643383279502884197169399D0)
      PARAMETER (QINNER = 0, QOUTER = 1)
      PARAMETER (SIG = 5.6696E-12)

      DIMENSION
     +   A2(A2DIM, A2DIM), AEX(5*PEX+2, 5*PEX+2),
     +   AIN(5*PIN+2, 5*PIN+2), AINN(2), AOUT(2), AREA(2),
     +   AREAC(RPNTS+1), AREAF(2), AREAW(2), CFLOW(2),
     +   CONDUC(8), CONVN(2), CONVS(2), CONW(RPNTS+1), CP(8), DH(2),
     +   DZEND(2), DZENDI(2), EMIS(2), EMISSI(5),
     +   F(GROUPA + COMPS*PMAX + GROUPB), FEE(2), FEG(2),
     +   FENDS(PDIM+1, NUMEX+1, 2), FENG(NUMWF, RPNTS+NUMEX),
     +   FESE(PDIM, 2), FEW(RPNTS, 2), FGE(2), FGSE(PDIM, 2),
     +   FSEE(PDIM, 2), FSEG(PDIM, 2), FSES(PDIM, PDIM, 2),
     +   FSEW(PDIM, RPNTS, 2), FSW(2*RPNTS), FWE(RPNTS, 2),
     +   FWS(2*RPNTS), FWSE(RPNTS, PDIM, 2), FWW(2*RPNTS, 2*RPNTS),
     +   GBOT(NUMTC), GTOP(NUMTC), HCOEF(5), HINN(2), HOUT(2),
     +   IBEG(3), IEND(3),
     +   IPVT2(A2DIM), IPVT3(5*PIN + 2), IPVT4(5*PEX + 2),
     +   J1(4, PDIM, 2), J2(A2DIM, NUMWF)

      DIMENSION
     +   JBOT(NUMTC), JEND(2), JEX(5*PEX + 2),
     +   JIN(5*PIN + 2), JTOP(NUMTC), JWALL(PDIM, 2), LOS1(NUMWF),
     +   LOS2(PDIM, 2), LOS3(2), MFLOW(2), MINJ(NUMWF), MINJE(PDIM, 2),
     +   MWZ(PDIM+1, 2), MWZW(NUMWF+1), NU(2), PCAP(2),
     +   POWER(NUMWF), POWERH(MAXHT),
     +   PWALL(PDIM, 2), QWAF(RPNTS, 2), R(RPNTS), RCPV(RPNTS),
     +   RE(2), REFL(2), RHO(8),
     +   RING(5), RINQ(2), ROUTQ(2), S(GROUPA + COMPS*PMAX + GROUPB),
     +   S0(GROUPA + COMPS*PMAX + GROUPB), SUMB(NUMTC), SUMT(NUMTC),
     +   TC(NUMTC), TCAP(2), TCAPS(2), TDAT(TEMMAX + 2), TEAST(NUMWF),
     +   TEND(2), TEND0(2), TEND4(2), TENDI(2), TENDI0(2), TENDI4(2),
     +   TEX(NUMWF, NUMEX), TEX0(NUMWF, NUMEX), TEX4(NUMWF, NUMEX),
     +   TMID(NUMWF), TRAN(2), TSIDE(PDIM, 2), TTGAS(2),
     +   TWAF(NUMWF+1, RPNTS),
     +   TWAF0(NUMWF+1, RPNTS), TWAF4(NUMWF+1, RPNTS), TWEST(NUMWF)

      DIMENSION
     +   TX(PDIM, NUMEX, 2), TX0(PDIM, NUMEX, 2), TX4(PDIM, NUMEX, 2),
     +   VINJ(MAXIN), VOL(2), WR(RPNTS+1), WZ(PDIM+1, 2),
     +   WZWAF(NUMWF+1), Z(PDIM,2),
     +   ZBEG(MAXHT), ZBEGI(MAXIN), ZDAT(TEMMAX + 2),
     +   ZEND(MAXHT), ZENDI(MAXIN), ZTC(NUMTC), ZWAF(NUMWF)

      DIMENSION TURB(2)

      DIMENSION SAV1(30), SAV2(30), SAV3(30)

      COMMON /PED1/ PEDRI, PEDRO, TTHIC
      COMMON /PED2/ PCON, PED

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  DEFINE STATEMENT FUNCTION FOR VISCOSITY OF AIR (G/(CM SEC)).
C     TAIR - TEMPERATURE OF AIR IN K.

      MUAIR(TAIR) = 1.716E-04*383.66/(TAIR + 110.56)*(TAIR/273.1)**1.5

C///  DEFINE STATEMENT FUNCTION FOR CONDUCTIVITY OF AIR (W/(CM K)).

      KAIR(TAIR) = 2.414E-04*467.5/(TAIR + 194.4)*(TAIR/273.1)**1.5

C///  DEFINE STATEMENT FUNCTION FOR DENSITY OF AIR AT 1 ATM (G/CM**3).

      RHOAIR(TAIR) = 1.0133E02/287./TAIR

C///  SET THE SPECIFIC HEAT OF AIR (J/G/K)

      CPAIR = 1.189

C///  PRANDTL NUMBER OF AIR

      PRAIR = 0.719

C///  SET GAP TREATMENT OPTION

      IOPT = 1

C///  DEFINE AREAS FOR QUARTZ WALLS AND HEATER

      A1M  = 2.*PI*RINQ(1)
      A1P  = 2.*PI*ROUTQ(1)
      A2M  = 2.*PI*RINQ(2)
      A2P  = 2.*PI*ROUTQ(2)
      AH   = 2.*PI*RTUBEI
      AWAF = 2.*PI*WR(RPNTS+1)

      AREA(1) = A1M
      AREA(2) = A2M

C     OVER THE WAFERS

      AREAW(1) = AWAF
      AREAW(2) = A2M

C      AREAW(1) = A1M
C      AREAW(2) = A2M

C///  CALCULATE SHAPEFACTORS FOR QUARTZ WALLS AND HEATER

C      FHH = 1.0 - A2P/AH
C      FH2 = A2P/AH
C      F2H = 1.0
C      F22 = 1.0 - A1P/A2M
C      F21E = 0.5*A1P/A2M
C      F21W = 0.5*A1P/A2M
C      F21 = A1P/A2M
C      F1E2 = 1.0
C      F1W2 = 1.0
C      F12 = 1.0

C///  MEAN SHAPE FACTORS

      FHH = 1.0 - AREA(2)/AH
      FH2 = AREA(2)/AH
      F2H = 1.0
      F22 = 1.0 - AREA(1)/AREA(2)
      F21 = AREA(1)/AREA(2)
      F12 = 1.0

C///  MODIFICATIONS FOR SHAPE FACTORS OVER WAFER REGION

      F12W = 1.0
      F21W = AREAW(1)/AREAW(2)
      F22W = 1.0 - AREAW(1)/AREAW(2)

C*****USE SAME WALL SHAPE FACTORS OVER WAFERS AND ENTRANCE/EXIT
C
C      F12W = F12
C      F21W = F21
C      F22W = F22
C
C*****END USE SAME WALL SHAPE FACTORS OVER WAFERS AND ENTRANCE/EXIT

C///  SET POWER, RADIOSITY AND TEMPERATURE VALUES TO ZERO.

      DO 1010 K = 1, NUMWF
         POWER(K) = 0.0
1010  CONTINUE

      TAMB = 0.0
      TAMB4 = 0.0

      DO 1020 LOC = 1, 2
         JEND(LOC) = 0.0
         TEND(LOC) = 0.0
         TEND4(LOC) = 0.0
         TEND0(LOC) = 0.0
         TENDI(LOC) = 0.0
         TENDI4(LOC) = 0.0
         TENDI0(LOC) = 0.0
         DO 1020 K = 1, PDIM
            JWALL(K, LOC) = 0.0
            PWALL(K, LOC) = 0.0
            DO 1020 I = 1, NUMEX
               TX(K, I, LOC) = 0.0
               TX4(K, I, LOC) = 0.0
               TX0(K, I, LOC) = 0.0
1020  CONTINUE

      DO 1040 J = 1, NUMWF
         DO 1030 K = 1, RPNTS
            TWAF(J, K) = 0.0
            TWAF4(J, K) = 0.0
            TWAF0(J, K) = 0.0
1030     CONTINUE
         DO 1040 I = 1, NUMEX
            TEX(J, I) = 0.0
            TEX4(J, I) = 0.0
            TEX0(J, I) = 0.0
1040  CONTINUE

C///  DECIDE WHICH TEMPERATURE VALUES TO UNPACK.

      IF (CHOICE .EQ. 0) THEN
         FLAG = .TRUE.
         FLAG0 = TIME
         FLAG00 = .FALSE.
         FLAG4 = .TRUE.
      ELSEIF (CHOICE .EQ. 4) THEN
         FLAG = .TRUE.
         FLAG0 = .FALSE.
         FLAG00 = .TRUE.
         FLAG4 = .TRUE.
      ELSE
         FLAG = CHOICE .EQ. 1
         FLAG0 = CHOICE .EQ. 2
         FLAG00 = .FALSE.
         FLAG4 = CHOICE .EQ. 3
      END IF

C///  UNPACK THE BOUNDARY VALUES.

      IF (FLAG) TAMB = BTEMP
      IF (FLAG4) TAMB4 = BTEMP**4

C///  GET THE POWER DISSIPATED IN EACH CELL OF THE WALL.

C*****PRINT HEATER POWERS
C      WRITE(TEXT, *) 'HEATER POWERS'
C      WRITE(TEXT, *) ('POWERH(',I,') = ',POWERH(I), I=1, NUMHT)
C*****END PRINT HEATER POWERS

C      SUBROUTINE POW
C     +  (ERROR, TEXT,
C     +   MAXHT, NUMHT, NUMWF, PDIM, PEX, PIN, POWER, POWERH, PWALL,
C     +   WZ, WZWAF, ZBEG, ZEND, ZLEN)

      CALL POW
     +  (ERROR, TEXT,
     +   MAXHT, NUMHT, NUMWF, PDIM, PEX, PIN, POWER, POWERH, PWALL,
     +   WZ, WZWAF, ZBEG, ZEND, ZLEN)

C///  GET COOLING GAS MASS FLOWRATE INJECTED IN EACH CELL AND
C     FLOWRATES AT CELL BOUNDARIES FOR OUTER FLOW.

C      SUBROUTINE MASSD
C     +  (ERROR, TEXT,
C     +   MAXIN, MFLOW, MINJ, MINJE, MWZ, MWZW, NUMIN, NUMWF, PDIM,
C     +   PEX, PIN, VINJ, WZ,
C     +   WZWAF, ZBEGI, ZENDI, ZLEN)

      CALL MASSD
     +  (ERROR, TEXT,
     +   MAXIN, MFLOW, MINJ, MINJE, MWZ, MWZW, NUMIN, NUMWF,
     +   PDIM, PEX, PIN, VINJ, WZ,
     +   WZWAF, ZBEGI, ZENDI, ZLEN)

C///  UNPACK THE DEPENDENT VARIABLES.

      COUNT = 0
      DO 1050 K = 1, PIN
C         COUNT = COUNT + 1
C         JWALL(K, 1) = S(COUNT)
         DO 1050 I = 1, NUMEX
            COUNT = COUNT + 1
            IF (FLAG) TX(K, I, 1) = S(COUNT)
            IF (FLAG0) TX0(K, I, 1) = (S(COUNT) - S0(COUNT)) / TSTEP
            IF (FLAG00) TX0(K, I, 1) = S0(COUNT)
            IF (FLAG4) TX4(K, I, 1) = S(COUNT)**4
1050  CONTINUE

C      COUNT = COUNT + 1
C      JEND(1) = S(COUNT)
      COUNT = COUNT + 1
      IF (FLAG) TEND(1) = S(COUNT)
      IF (FLAG0) TEND0(1) = (S(COUNT) - S0(COUNT)) / TSTEP
      IF (FLAG00) TEND0(1) = S0(COUNT)
      IF (FLAG4) TEND4(1) = S(COUNT)**4
      COUNT = COUNT + 1
      IF (FLAG) TENDI(1) = S(COUNT)
      IF (FLAG0) TENDI0(1) = (S(COUNT) - S0(COUNT)) / TSTEP
      IF (FLAG00) TENDI0(1) = S0(COUNT)
      IF (FLAG4) TENDI4(1) = S(COUNT)**4

      DO 1070 J = 1, NUMWF
         DO 1060 K = 1, RPNTS
            COUNT = COUNT + 1
            IF (FLAG) TWAF(J, K) = S(COUNT)
            IF (FLAG0) TWAF0(J, K) = (S(COUNT) - S0(COUNT)) / TSTEP
            IF (FLAG00) TWAF0(J, K) = S0(COUNT)
            IF (FLAG4) TWAF4(J, K) = S(COUNT)**4
1060     CONTINUE
         DO 1070 I = 1, NUMEX
            COUNT = COUNT + 1
            IF (FLAG) TEX(J, I) = S(COUNT)
            IF (FLAG0) TEX0(J, I) = (S(COUNT) - S0(COUNT)) / TSTEP
            IF (FLAG00) TEX0(J, I) = S0(COUNT)
            IF (FLAG4) TEX4(J, I) = S(COUNT)**4
1070  CONTINUE

      DO 1080 K = PEX, 1, - 1
C         COUNT = COUNT + 1
C         JWALL(K, 2) = S(COUNT)
         DO 1080 I = 1, NUMEX
            COUNT = COUNT + 1
            IF (FLAG) TX(K, I, 2) = S(COUNT)
            IF (FLAG0) TX0(K, I, 2) = (S(COUNT) - S0(COUNT)) / TSTEP
            IF (FLAG00) TX0(K, I, 2) = S0(COUNT)
            IF (FLAG4) TX4(K, I, 2) = S(COUNT)**4
1080  CONTINUE

C      COUNT = COUNT + 1
C      JEND(2) = S(COUNT)
      COUNT = COUNT + 1
      IF (FLAG) TEND(2) = S(COUNT)
      IF (FLAG0) TEND0(2) = (S(COUNT) - S0(COUNT)) / TSTEP
      IF (FLAG00) TEND0(2) = S0(COUNT)
      IF (FLAG4) TEND4(2) = S(COUNT)**4
      COUNT = COUNT + 1
      IF (FLAG) TENDI(2) = S(COUNT)
      IF (FLAG0) TENDI0(2) = (S(COUNT) - S0(COUNT)) / TSTEP
      IF (FLAG00) TENDI0(2) = S0(COUNT)
      IF (FLAG4) TENDI4(2) = S(COUNT)**4

C///  SET SOME USEFUL PARAMETERS

      RINO = RTUBEO + THICKI
      RINI = RTUBEO
      RIN = (RINO + RINI)/2.
      RTUBE = (RTUBEO + RTUBEI)/2.
      FACTOR = (RIN - RTUBEO)/(RIN - RTUBE)
      COND = 1./((1.-FACTOR)/CONDUC(2) + FACTOR/CONDUC(3))

C///  CHECK THERMOCOUPLE LOCATIONS.

      ERROR = .NOT. (NUMTC.GT.0)
      IF (ERROR) GO TO 9101

      DO 1090 I = 1, NUMTC
         ERROR = .NOT. ((0.0.LE.ZTC(I)).AND.
     +                 (ZTC(I).LE.ZLEN))
         IF (ERROR) GO TO 9102
1090  CONTINUE

C///  APPLY CORRECTION TO SHAPEFACTORS

      CORR = .TRUE.

C///  THERMOCOUPLE WIDTH AND RADIUS

      DZ = 0.2
      RTC = 0.5*(ROUTQ(1) + RINQ(2))
      RTC = ROUTQ(1) + DZ

C///  ZERO THERMOCOUPLE IRRADIANCE

      DO 1100 I = 1, NUMTC
         GBOT(I) = 0.0
         GTOP(I) = 0.0
         SUMB(I) = 0.0
         SUMT(I) = 0.0
         JTOP(I) = 0.0
         JBOT(I) = 0.0
1100  CONTINUE

C///  INFORMATION NEEDED FOR GAS-FLOW COOLING MODELS

C     DEFINE AREAS FOR FLOW AND CONVECTIVE HEAT TRANSFER.

C     FLOW AREAS FOR ANNULLI.
      AREAF(1) = PI*(RINQ(2)**2 - ROUTQ(1)**2)
      AREAF(2) = PI*(RTUBEI**2 - ROUTQ(2)**2)

C     INSIDE AND OUTSIDE AREAS AND VOLUME TERMS OF ANNULLI.
      AOUT(1)  = 2.*PI*RINQ(2)
      AINN(1)  = 2.*PI*ROUTQ(1)
      AOUT(2)  = 2.*PI*RTUBEI
      AINN(2)  = 2.*PI*ROUTQ(2)
      VOL(1)   = PI*(RINQ(2)**2 - ROUTQ(1)**2)
      VOL(2)   = PI*(RTUBEI**2 - ROUTQ(2)**2)

C     HYDRAULIC DIAMETERS.
      DH(1) = 2.*(RINQ(2) - ROUTQ(1))
      DH(2) = 2.*(RTUBEI - ROUTQ(2))

C     SET REFERENCE TEMPERATURE FOR AIR VISCOSITY TO 1000K.
      TREF = 1000.0

C     REYNOLDS NUMBERS.
      RE(1) = ABS(MFLOW(1))/AREAF(1)*DH(1)/MUAIR(TREF)
      RE(2) = ABS(MFLOW(2))/AREAF(2)*DH(2)/MUAIR(TREF)

C     LET SUBROUTINE NUSSL USE THE COMBINED CORRELATION FOR
C     LAMINAR AND TURBULENT FLOW.
      TURB(1) = .FALSE.
      TURB(2) = .FALSE.

C///////////////////////////////////////////////////////////////////////
C
C     (3) ENERGY BALANCE TERMS FOR ENTRANCE AND EXIT REGION.
C
C///////////////////////////////////////////////////////////////////////

C      SUBROUTINE EBAL1
C     +  (ERROR, TEXT,
C     +   AEX, AIN, AREA, CFLOW, COMPS, COND, CONDUC, CP, DZ, DZEND,
C     +   DZENDI, EMIS, EMISS, EMISSE, EMISSI, EMISSW, F12, F21, F22,
C     +   F2H, FEE, FEG, FENDS, FESE, FEW, FGE, FGSE, FH2, FHH, FSEE,
C     +   FSEG, FSES, FSEW, FWE, FWSE, GBOT, GROUPA, GROUPB, GTOP,
C     +   HCOEF, IBEG, IEND, IPVT3, IPVT4, J1, JBOT, JEND, JEX,
C     +   JIN, JTOP, JWALL, LOS1, LOS2, LOS3, MAXHT, MFLOW, MINJE,
C     +   MWZ, NUMEX, NUMHT, NUMT, NUMTC, NUMWF, PCAP, PDIM, PEX,
C     +   PIN, PMAX, POWERH, PWALL, QWAF, R, REFL, RHO, RIN, RINI,
C     +   RINJEC,
C     +   RINO, RINQ, ROUTQ, RPNTS, RTC, RTUBE, RTUBEI, RTUBEO, SUMB,
C     +   SUMT, TAMB, TAMB4, TC, TCAP, TCAPS, TEMMAX, TEND, TEND0,
C     +   TEND4, TENDI, TENDI0, TENDI4, TEX, TEX0, TEX4, TFIX, THICKI,
C     +   TIME, TRAN, TSIDE, TSTEP, TTGAS, TURB, TWAF, TWAF0, TWAF4,
C     +   TX, TX0, TX4, WR, WZ, Z, ZBEG, ZDAT, ZEND, ZLEN, ZTC, ZWAF)

      CALL EBAL1
     +  (ERROR, TEXT,
     +   AEX, AIN, AREA, CFLOW, COMPS, COND, CONDUC, CP, DZ, DZEND,
     +   DZENDI,
     +   EMIS, EMISS, EMISSE, EMISSI, EMISSW, F12, F21, F22, F2H, FEE,
     +   FEG, FENDS, FESE, FEW, FGE, FGSE, FH2, FHH, FSEE, FSEG, FSES,
     +   FSEW, FWE, FWSE, GBOT, GROUPA, GROUPB, GTOP, HCOEF, IBEG,
     +   IEND, IPVT3, IPVT4, J1, JBOT, JEND, JEX, JIN, JTOP,
     +   JWALL, LOS1, LOS2, LOS3, MAXHT, MFLOW, MINJE, MWZ,
     +   NUMEX, NUMHT, NUMT, NUMTC, NUMWF, PCAP, PDIM, PEX, PIN, PMAX,
     +   POWERH, PWALL, QWAF, R, REFL, RHO, RIN, RINI, RINJEC,
     +   RINO, RINQ,
     +   ROUTQ, RPNTS, RTC, RTUBE, RTUBEI, RTUBEO, SUMB, SUMT, TAMB,
     +   TAMB4, TC, TCAP, TCAPS, TEMMAX, TEND, TEND0, TEND4, TENDI,
     +   TENDI0, TENDI4, TEX, TEX0, TEX4, TFIX, THICKI, TIME, TRAN,
     +   TSIDE, TSTEP, TTGAS, TURB, TWAF, TWAF0, TWAF4, TX, TX0,
     +   TX4, WR, WZ, Z, ZBEG, ZDAT, ZEND, ZLEN, ZTC, ZWAF)

C///////////////////////////////////////////////////////////////////////
C
C     (5) ENERGY BALANCE FOR WAFERS.
C
C///////////////////////////////////////////////////////////////////////

C*****PRINT SEPARATION BARS
C      WRITE(TEXT, *) ''
C      WRITE(TEXT, *) '--------- WAFER SECTION ------------------'
C      WRITE(TEXT, *) ' '
C*****END PRINT SEPARATION BARS

C///  COMPUTE WAFER CELL AREAS AND RHO*CP*VOL PRODUCTS.

      IF (RBOAT) THEN

C     RING BOAT SECTION
         DO 5001 WALL = 1, RPNTS+1

C        WALL AREAS
         IF (WR(WALL) .GT. RING(1)) THEN
            AREAC(WALL) = 2.*PI*(THICK+RING(2))*WR(WALL)
         ELSE
            AREAC(WALL) = 2.*PI*THICK*WR(WALL)
         ENDIF

C        CONDUCTIVITY AT WALL
         IF (WR(WALL) .GE. RING(1)) THEN
            CONW(WALL) = CONDUC(1)
         ELSE
            CONW(WALL) = (THICK*CONDUC(1) + RING(2)*RING(5))/(THICK
     +                   + RING(2))
         ENDIF
5001     CONTINUE

C        RHO*CP*VOL TERMS
         DO 5002 CELL = 1, RPNTS
         IF (WR(CELL) .GE. RING(1)) THEN
            RCPV(CELL) = PI*(WR(CELL+1)**2 -
     +                WR(CELL)**2)*(THICK*RHO(1)*CP(1)
     +                + RING(2)*RING(3)*RING(3))
         ELSE
            RCPV(CELL) = PI*(WR(CELL+1)**2
     +                   - WR(CELL)**2)*THICK*RHO(1)*CP(1)
         ENDIF
5002     CONTINUE

      ELSE

C     WITHOUT RING BOAT SECTION
         DO 5003 J = 1, RPNTS
         AREAC(J) = 2.*PI*THICK*WR(J)
         CONW(J) = CONDUC(1)
         RCPV(J) = PI*(WR(J+1)**2 - WR(J)**2)*THICK*RHO(1)*CP(1)
5003     CONTINUE
         AREAC(RPNTS+1) = 2.*PI*THICK*WR(RPNTS+1)
         CONW(RPNTS+1) = CONDUC(1)
      ENDIF

C///  SAVE VALUES OF AREAC, CONW, AND RCPV

      DO 5004 J = 1, RPNTS
         SAV1(J) = AREAC(J)
         SAV2(J) = CONW(J)
         SAV3(J) = RCPV(J)
5004  CONTINUE
      SAV1(RPNTS+1) = AREAC(RPNTS+1)
      SAV2(RPNTS+1) = CONW(RPNTS+1)

C///  TOP OF THE LOOP OVER WAFERS.

      DO 5020 WAFER = 1, NUMWF-1

      SCALE = 1.0

C*****USE WALL RADIOSITIES SCALE FACTORS
C
C      SCALE = (2.0*PI*WR(RPNTS+1)*SPACE)/((WZWAF(WAFER+1) -
C     +            WZWAF(WAFER))*AREAW(1))
C
C*****END USE WALL RADIOSITIES SCALE FACTORS

C///  REACTOR WALL RADIOSITIES.

      SUMWE = 0.0
      SUMWW = 0.0
      DO 5010 J = 1, RPNTS
         IF (WAFER .NE. NUMWF)
     +      SUMWE = SUMWE + EMISS*TWAF4(WAFER, J)*FSW(J)
     +            + EMISS*TWAF4(WAFER+1, J)*FSW(J+RPNTS)
         IF (WAFER .NE. 1)
     +      SUMWW = SUMWW + EMISS*TWAF4(WAFER, J)*FSW(J)
     +            + EMISS*TWAF4(WAFER-1, J)*FSW(J+RPNTS)
5010  CONTINUE

C///  FIND THE RADIOSITIES ON THE HEATER AND QUARTZ WALLS FOR
C     AN ELEMENT OVER THE WAFERS.

      J2(1, WAFER) = EMISSW*SIG*TEX4(WAFER, 3)
      J2(2, WAFER) = EMIS(2)*SIG*TEX4(WAFER, 2)
      J2(3, WAFER) = EMIS(2)*SIG*TEX4(WAFER, 2)
      J2(4, WAFER) = EMIS(1)*SIG*TEX4(WAFER, 1) +
     +        SCALE*TRAN(1)*SIG*SUMWE +
     +        TRAN(1)*FSS/(1.- REFL(1)*FSS)*(SCALE*REFL(1)*SIG*SUMWE +
     +        EMIS(1)*SIG*TEX4(WAFER, 1))

C*****WRITE RADIOSITIES OVER WAFERS
C      WRITE(TEXT, *) ' '
C      WRITE(TEXT, *) 'WAFER = ' , WAFER
C      WRITE(TEXT, *) 'J2 FOR RHS'
C      WRITE(TEXT, *) 'J2(1) = ', J2(1, WAFER)
C      WRITE(TEXT, *) 'J2(2) = ', J2(2, WAFER)
C      WRITE(TEXT, *) 'J2(3) = ', J2(3, WAFER)
C      WRITE(TEXT, *) 'J2(4) = ', J2(4, WAFER)
C*****END WRITE RADIOSITIES OVER WAFERS

C*****DOUBLE PRECISION
      CALL DGESL(A2, A2DIM, A2DIM, IPVT2, J2(1, WAFER), 0)
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      CALL SGESL(A2, A2DIM, A2DIM, IPVT2, J2(1, WAFER), 0)
C*****END SINGLE PRECISION

C*****WRITE RADIOSITIES OVER WAFERS
C      WRITE(TEXT, *) ' '
C      WRITE(TEXT, *) 'J2(1) = ', J2(1, WAFER)
C      WRITE(TEXT, *) 'J2(2) = ', J2(2, WAFER)
C      WRITE(TEXT, *) 'J2(3) = ', J2(3, WAFER)
C      WRITE(TEXT, *) 'J2(4) = ', J2(4, WAFER)
C      WRITE(TEXT, *) 'A1M = ', A1M, ' A1P = ', A1P
C      WRITE(TEXT, *) 'A2M = ', A2M, ' A2M = ', A2M
C      WRITE(TEXT, *) 'AH = ', AH
C*****END WRITE RADIOSITIES OVER WAFERS

C///  ACCUMULATE IRRADIANCE ON EACH THERMOCOUPLE OVER WAFERS

      DO 5040 I = 1, NUMTC

C        SUBROUTINE CONCC
C     +    (ERROR, TEXT,
C     +     F12, F21, R1, R2, Z1H, Z1L, Z2H, Z2L)

         CALL CONCC(ERROR, TEXT, FTOP, DUMMY, RTC, RINQ(2), ZTC(I)+DZ,
     +             ZTC(I), WZWAF(WAFER+1), WZWAF(WAFER))
         IF (ERROR) GO TO 9103

         SUMT(I) = SUMT(I) + FTOP

         CALL CONCC(ERROR, TEXT, DUMMY, FBOT, WR(RPNTS+1), RTC,
     +             WZWAF(WAFER+1), WZWAF(WAFER), ZTC(I)+DZ, ZTC(I))
         IF (ERROR) GO TO 9103

         SUMB(I) = SUMB(I) + FBOT

         GTOP(I) = GTOP(I) + J2(3, WAFER)*FTOP
         GBOT(I) = GBOT(I) + J2(4, WAFER)*FBOT

         IF((WZWAF(WAFER).LT.ZTC(I)).AND.
     +      (ZTC(I).LE.WZWAF(WAFER+1))) THEN
            JTOP(I) = J2(3, WAFER)
            JBOT(I) = J2(4, WAFER)
         ENDIF

C*****WRITE THERMOCOUPLE SHAPEFACTORS
C      WRITE(TEXT, *) ' '
C      WRITE(TEXT, *) 'WAFER = ', WAFER
C      WRITE(TEXT, *) 'TC NUMBER = ', I
C      WRITE(TEXT, *) 'FTOP = ', FTOP
C      WRITE(TEXT, *) 'FBOT = ', FBOT
C      WRITE(TEXT, *) 'FBOT = ', FBOT
C*****END WRITE THERMOCOUPLE SHAPEFACTORS

5040  CONTINUE

C///  BOTTOM OF THE LOOP OVER WAFERS.

5020  CONTINUE

C----------------------------------------------------------------------

C///  TOP OF THE LOOP OVER WAFERS.

      DO 5080 WAFER = 1, NUMWF

C///  TREAT THE FIRST WAFER AS A QUARTZ DISC

      IF (PED .AND. (WAFER.EQ.1)) THEN

         DO 5025 J = 1, RPNTS
            AREAC(J) = 2.*PI*TTHIC*WR(J)
            CONW(J) = CONDUC(4)
            RCPV(J) = PI*(WR(J+1)**2 - WR(J)**2)*TTHIC*RHO(4)*CP(4)
5025     CONTINUE
         AREAC(RPNTS+1) = 2.*PI*TTHIC*WR(RPNTS+1)
         CONW(RPNTS+1) = CONDUC(4)

      ELSE

         DO 5026 I = 1, RPNTS
            AREAC(I) = SAV1(I)
            CONW(I) = SAV2(I)
            RCPV(I) = SAV3(I)
5026     CONTINUE
         AREAC(RPNTS+1) = SAV1(RPNTS+1)
         CONW(RPNTS+1) = SAV2(RPNTS+1)

      ENDIF

      SCALEE = 1.0
      SCALEW = 1.0

C*****USE WALL RADIOSITIES SCALE FACTORS
C
C      IF (WAFER .NE. NUMWF) THEN
C         SCALEE = ((WZWAF(WAFER+1) -
C     +            WZWAF(WAFER))*AREAW(1))/(2.0*PI*WR(RPNTS+1)*SPACE)
C      ELSE
C         SCALEE = 0.0
C      END IF
C
C      IF (WAFER .NE. 1) THEN
C         SCALEW = ((WZWAF(WAFER) -
C     +            WZWAF(WAFER-1))*AREAW(1))/(2.0*PI*WR(RPNTS+1)*SPACE)
C      ELSE
C         SCALEW = 0.0
C      END IF
C
C*****END USE WALL RADIOSITIES SCALE FACTORS

C///  REACTOR WALL RADIOSITIES.

      SUMWE = 0.0
      SUMWW = 0.0
      DO 5030 J = 1, RPNTS
         IF (WAFER .NE. NUMWF)
     +      SUMWE = SUMWE + EMISS*TWAF4(WAFER, J)*FSW(J)
     +            + EMISS*TWAF4(WAFER+1, J)*FSW(J+RPNTS)
         IF (WAFER .NE. 1)
     +      SUMWW = SUMWW + EMISS*TWAF4(WAFER, J)*FSW(J)
     +            + EMISS*TWAF4(WAFER-1, J)*FSW(J+RPNTS)
5030  CONTINUE

      IF (WAFER .NE. NUMWF) THEN
         JEAST = (EMIS(1) * SIG * TEX4(WAFER, 1)
     +         + REFL(1)/SCALEE * SIG * SUMWE)/ (1. - REFL(1) * FSS)
     +         + TRAN(1)*F12W*J2(3, WAFER)/ (1. - REFL(1) * FSS)
      ELSE
         JEAST = 0.0
      END IF

      IF (WAFER .NE. 1) THEN
         JWEST = (EMIS(1) * SIG * TEX4(WAFER-1, 1)
     +         + REFL(1)/SCALEW * SIG * SUMWW) / (1. - REFL(1) * FSS)
     +         + TRAN(1)*F12W*J2(3, WAFER-1)/ (1. - REFL(1) * FSS)
      ELSE
         JWEST = 0.0
      END IF

C     JEAST AND JWEST EXPRESSIONS BEFORE QUARTZ WALLS WERE ADDED
C
C         JEAST
C     +      = (EMISSW * SIG * TEX4(WAFER, 1)
C     +      + (1. - EMISSW) * SIG * SUMWE) / (1. - (1. - EMISSW) * FSS)
C         JWEST
C     +      = (EMISSW * SIG * TEX4(WAFER, 1)
C     +      + (1. - EMISSW) * SIG * SUMWW) / (1. - (1. - EMISSW) * FSS)

C///  TOP OF THE LOOP OVER CELLS IN THE WAFER.

      DO 5060 CELL = 1, RPNTS

C///  SUMMATIONS.

      SUME = 0.0
      SUMW = 0.0
      DO 5050 J = 1, RPNTS
         IF (WAFER .NE. NUMWF)
     +      SUME = SUME + EMISS*TWAF4(WAFER, J)*FWW(CELL, J)
     +           + EMISS*TWAF4(WAFER+1, J)*FWW(CELL, J+RPNTS)
         IF (WAFER .NE. 1)
     +      SUMW = SUMW + EMISS*TWAF4(WAFER, J)*FWW(CELL, J)
     +           + EMISS*TWAF4(WAFER-1, J)*FWW(CELL, J+RPNTS)
5050  CONTINUE

C///  IRRADIATION ON WAFER CELL SIDES.

      GEAST = JEAST*SCALEE*FWS(CELL) + SIG*SUME
      GWEST = JWEST*SCALEW*FWS(CELL) + SIG*SUMW

      QEAST = EMISS*SIG*TWAF4(WAFER, CELL) - EMISS*GEAST
      QWEST = EMISS*GWEST - EMISS*SIG*TWAF4(WAFER, CELL)

      IF (WAFER .EQ. 1) QWEST = QWAF(CELL, 1)

      IF (WAFER .EQ. NUMWF) QEAST = -QWAF(CELL, 2)

C*****SPECIFY HEAT FLUX ON INSIDE OF END WAFERS
C
C      IF (WAFER .EQ. 1)
C     +    QEAST = -100.0/(PI*(WR(CELL+1)**2 - WR(CELL)**2))
C      IF (WAFER .EQ. NUMWF)
C     +    QWEST = 100.0/(PI*(WR(CELL+1)**2 - WR(CELL)**2))
C
C*****END SPECIFY HEAT FLUX ON INSIDE OF END WAFERS

C*****SPECIFY HEAT FLUX ON THE OUTSIDE OF WAFERS 2 AND 3
C
C      IF (WAFER .EQ. 2)
C     +    QWEST = -50.0/(PI*(WR(CELL+1)**2 - WR(CELL)**2))
C      IF (WAFER .EQ. 3)
C     +    QEAST = 40.0/(PI*(WR(CELL+1)**2 - WR(CELL)**2))
C
C*****END SPECIFY HEAT FLUX ON THE OUTSIDE OF WAFERS 2 AND 3

C///  HEAT CONDUCTION FLUXES AT WAFER CELL BOUNDARIES.

      IF (CELL .EQ. RPNTS) THEN
C         QNORTH = SIG*(TWAF4(WAFER, CELL) - TEX4(WAFER, 1))
C     +      / (1./EMISS + 1./EMISSW - 1.)
          QNORTH = EMISS*SIG*TWAF4(WAFER, CELL) - EMISS*JEAST
          IF(WAFER.EQ.NUMWF) QNORTH = EMISS*SIG*TWAF4(WAFER, CELL)
     +                                -EMISS*JWEST
         QTIP = AREAC(RPNTS+1)*QNORTH
         IF (WAFER .EQ. NUMWF-1) THEN
            QTIP2 =
     +      AREAC(RPNTS+1)*(EMISS*SIG*TWAF4(WAFER+1, CELL)
     +      - EMISS*JEAST)
         ELSE
         ENDIF
C*****INSULATE TOPS OF WAFERS
C         QNORTH = 0.0
C         QTIP = 0.0
C         QTIP2 = 0.0
C*****END INSULATE TOPS OF WAFERS
      ELSE
         QNORTH = -CONW(CELL+1)*(TWAF(WAFER, CELL+1) -
     +            TWAF(WAFER, CELL))/(R(CELL+1) - R(CELL))
      END IF

      IF (CELL .EQ. 1) THEN
         QSOUTH = 0.0
      ELSE
         QSOUTH = -CONW(CELL)*(TWAF(WAFER, CELL) -
     +            TWAF(WAFER, CELL-1))/(R(CELL) - R(CELL-1))
      END IF

C*****INSULATE OUTSIDE AND TOP OF END WAFERS
C      IF (WAFER .EQ. 1) QWEST = 0.0
C      IF (WAFER .EQ. NUMWF) QEAST = 0.0
C      IF ((WAFER .EQ. 1) .OR. (WAFER. EQ. NUMWF)) QNORTH = 0.0
C*****END INSULATE OUTSIDE AND TOP OF END WAFERS

C*****INSULATE INSIDE AND TOP OF END WAFERS
C      IF (WAFER .EQ. 1) QEAST = 0.0
C      IF (WAFER .EQ. NUMWF) QWEST = 0.0
C      IF (((WAFER.EQ.1).OR.(WAFER.EQ.NUMWF)) .AND. (CELL.EQ.RPNTS))
C     +   QNORTH = 0.0
C*****END INSULATE INSIDE AND TOP OF END WAFERS

C*****WRITE WAFER ENERGY BALANCES
C      WRITE(TEXT, *) 'WAFER ENERGY BALANCE'
C      WRITE(TEXT, *) 'WAFER = ', WAFER, ' CELL = ', CELL
C      WRITE(TEXT, *) 'QNORTH = ', QNORTH, ' QSOUTH = ', QSOUTH
C      WRITE(TEXT, *) 'QEAST = ', QEAST, ' QWEST = ', QWEST
C      WRITE(TEXT, *) 'NORTH TERM = ', AREAC(CELL+1)*QNORTH
C      WRITE(TEXT, *) 'SOUTH TERM = ', AREAC(CELL)*QSOUTH
C      WRITE(TEXT, *) 'EAST TERM = ', PI*(WR(CELL+1)**2 -
C     +                WR(CELL)**2)*QEAST
C      WRITE(TEXT, *) 'WEST TERM = ', PI*(WR(CELL+1)**2 -
C     +                WR(CELL)**2)*QWEST
C      WRITE(TEXT, *) 'RHO*CP*VOL TERM = ',RCPV(CELL)
C      IF (CELL .EQ. RPNTS) WRITE(TEXT, *) 'QTIP = ', QTIP
C*****END WRITE WAFER ENERGY BALANCES

C///  WAFER CELL ENERGY BALANCES.

      FENG(WAFER, CELL) =
     +   AREAC(CELL+1)*QNORTH - AREAC(CELL)*QSOUTH +
     +   PI*(WR(CELL+1)**2 - WR(CELL)**2)*(QEAST - QWEST)

C     DIVIDED VOL*RHO*CP FOR THE ELEMENT

      FENG(WAFER, CELL) =
     +   FENG(WAFER, CELL)/RCPV(CELL)

C///  TRANSIENT TERM.

      IF (TIME) THEN
         FENG(WAFER, CELL) =
     +      FENG(WAFER, CELL) + TWAF0(WAFER, CELL)
      END IF

C*****SET CONSTANT WAFER TEMPERATURE
C      FENG(WAFER, CELL) = TWAF(WAFER, CELL) - 800.0
C*****END SET CONSTANT WAFER TEMPERATURE

C*****FIX THE TEMPERATURE OF THE END WAFERS
C      IF (WAFER .EQ. 1) FENG(WAFER, CELL) = TWAF(WAFER, CELL) - 800.0
C      IF (WAFER .EQ. NUMWF)
C     +   FENG(WAFER, CELL) = TWAF(WAFER, CELL) - 800.0
C*****END FIX THE TEMPERATURE OF THE END WAFERS

C///  BOTTOM OF THE LOOP OVER CELLS.

5060  CONTINUE

      IF (WAFER .EQ. NUMWF) GO TO 5090

C///  ENERGY BALANCE ON EXTERNAL COMPONENTS

C///  TOP OF THE LOOP FOR ENERGY BALANCE OVER THE QUARTZ WALLS.

      DO 5070 WALL = 1, 2

C///  COMPUTE IRRADIATION ON THE FIRST QUARTZ WALL AND OTHER
C///  NEEDED TERMS.

      IF (WALL. EQ. 1) THEN
         GEAST = JEAST*FSS + SIG*SUMWE
         GWEST = JWEST*FSS + SIG*SUMWW
         QNORTH = J2(4, WAFER) - J2(3, WAFER)*F12W
         QSOUTH = GEAST - JEAST
      ELSE
         QNORTH = J2(2, WAFER) - J2(1, WAFER)*F2H
         QSOUTH = J2(4, WAFER)*F21W + J2(3, WAFER)*F22W
     +            - J2(3, WAFER)
      END IF

      IF (WAFER .EQ. 1) THEN
         QWEST = -CONDUC(4)*(TEX(WAFER, WALL) -
     +            TX(PIN, WALL, 1))/(ZWAF(WAFER) - Z(PIN, 1))
C*****INSULATE QUARTZ WALLS ENTRANCE AND EXIT FROM MIDSECTION
C         QWEST = 0.0
C*****END INSULATE QUARTZ WALLS ENTRANCE AND EXIT FROM MIDSECTION
      ELSE
         QWEST = -CONDUC(4)*(TEX(WAFER, WALL) -
     +            TEX(WAFER-1, WALL))/(ZWAF(WAFER) - ZWAF(WAFER-1))
      END IF

      IF (WAFER .EQ. (NUMWF-1)) THEN
         QEAST = -CONDUC(4)*(TEX(WAFER, WALL) -
     +            TX(PEX, WALL, 2))/((ZLEN - ZWAF(WAFER)) - Z(PEX, 2))
C        CHANGE SIGN BECAUSE OF CHANGE IN COORDINATE SYSTEM
         QEAST = -QEAST
C*****INSULATE QUARTZ WALLS ENTRANCE AND EXIT FROM MIDSECTION
C         QEAST = 0.0
C*****END INSULATE QUARTZ WALLS ENTRANCE AND EXIT FROM MIDSECTION
      ELSE
         QEAST = -CONDUC(4)*(TEX(WAFER+1, WALL) -
     +            TEX(WAFER, WALL))/(ZWAF(WAFER+1) - ZWAF(WAFER))
      END IF

C*****ZERO SIDE HEAT FLUXES FOR QUARTZ WALL 1 AND 2
C      QEAST = 0.0
C      QWEST = 0.0
C*****END ZERO SIDE HEAT FLUXES FOR QUARTZ WALL 1 AND 2

C*****INSULATE THE INSIDE OF QUARTZ WALL 1
C      IF (WALL .EQ. 1) THEN
C         QSOUTH = 0.0
C      END IF
C*****END INSULATE THE INSIDE OF QUARTZ WALL 1

C*****INSULATE THE OUTSIDE OF QUARTZ WALL 1 OVER WAFERS
C      IF (WALL .EQ. 1) THEN
C         QNORTH = 0.0
C      END IF
C*****END INSULATE THE OUTSIDE OF QUARTZ WALL 1 OVER WAFERS

C*****FIX HEAT FLUX ON OUTSIDE OF WALL 1 OVER WAFERS
C      QNORTH = -200.0/(AREAW(WALL)*(WZWAF(WAFER+1) - WZWAF(WAFER)))
C*****END FIX HEAT FLUX ON OUTSIDE OF WALL 1 OVER WAFERS

C*****INSULATE THE INSIDE OF QUARTZ WALL 2
C      IF (WALL .EQ. 2) THEN
C         QSOUTH = 0.0
C      END IF
C*****END INSULATE THE INSIDE OF QUARTZ WALL 2

C///  COOLING GAS CONTRIBUTION

      IF (CFLOW(1)) THEN
          TFILM = (TEX(WAFER, 4+1) + TEX(WAFER, 1))/2.
          REYN = ABS(MFLOW(1))/AREAF(1)*DH(1)/MUAIR(TFILM)
          ACOND = KAIR(TFILM)
C         SUBROUTINE NUSSL(ERROR, TEXT,
C     +                    DIA, H, COND, RE, PR, NU, TURB)
          CALL NUSSL(ERROR, TEXT,
     +              DH(1), HCOF, ACOND, REYN, PRAIR, NU(1), TURB(1))

C*****FIX CONVECTIVE HEAT TRANSFER COEFFICIENT
C          HCOF = 25.0E-04
C*****END FIX CONVECTIVE HEAT TRANSFER COEFFICIENT

          CONVN(1) = HCOF*(TEX(WAFER, 4+1) -
     +              TEX(WAFER, 2))*AOUT(1)*(WZWAF(WAFER+1)
     +              - WZWAF(WAFER))

          CONVS(1) = HCOF*(TEX(WAFER, 1) -
     +              TEX(WAFER, 4+1))*AINN(1)*(WZWAF(WAFER+1)
     +              - WZWAF(WAFER))
      ELSE
          CONVN(1) = 0.0
          CONVS(1) = 0.0
      ENDIF

      IF (CFLOW(2)) THEN
          TFILM = (TEX(WAFER, 4+2) + TEX(WAFER, 2))/2.
          REYN = ABS(MWZW(WAFER))/AREAF(2)*DH(2)/MUAIR(TFILM)
          ACOND = KAIR(TFILM)
C         SUBROUTINE NUSSL(ERROR, TEXT,
C     +                    DIA, H, COND, RE, PR, NU, TURB)
          CALL NUSSL(ERROR, TEXT,
     +              DH(2), HCOF, ACOND, REYN, PRAIR, NU(2), TURB(2))

C*****PRINT HEAT TRANSFER CORRELATION INFORMATION FOR OUTER FLOW
C          WRITE(TEXT, *)'--HEAT TRAN. COEF. FOR OUTER FLOW (TWAF1)--'
C          WRITE(TEXT, *)'WAFER = ', WAFER
C          WRITE(TEXT, *)'CFLOW(2) = ', CFLOW(2)
C          WRITE(TEXT, *)'TURB(2) = ',TURB(2)
C          WRITE(TEXT, *)'MWZW(WAFER) = ', MWZW(WAFER)
C          WRITE(TEXT, *)'REYN = ',REYN
C          WRITE(TEXT, *)'DH(2) = ',DH(2)
C          WRITE(TEXT, *)'PRAIR = ',PRAIR
C          WRITE(TEXT, *)'ACOND = ',ACOND
C          WRITE(TEXT, *)'NU(2) = ',NU(2)
C          WRITE(TEXT, *)'HCOF = ',HCOF
C*****END PRINT HEAT TRANSFER CORRELATION INFORMATION FOR OUTER FLOW

C*****FIX CONVECTIVE HEAT TRANSFER COEFFICIENT
C          HCOF = 25.0E-04
C*****END FIX CONVECTIVE HEAT TRANSFER COEFFICIENT

          CONVN(2) = HCOF*(TEX(WAFER, 4+2) -
     +              TEX(WAFER, 3))*AOUT(2)*(WZWAF(WAFER+1)
     +              - WZWAF(WAFER))
          CONVS(2) = HCOF*(TEX(WAFER, 2) -
     +              TEX(WAFER, 4+2))*AINN(2)*(WZWAF(WAFER+1)
     +              - WZWAF(WAFER))
      ELSE
          CONVN(2) = 0.0
          CONVS(2) = 0.0
      ENDIF

C*****FIX CONVECTIVE COOLING HEAT TRANSFER FLUXES
C
C      QFLUX = 100.0
C      IF (CFLOW(1)) THEN
C          CONVN(1) = QFLUX*AOUT(1)*(WZWAF(WAFER+1) - WZWAF(WAFER))
C          CONVS(1) = QFLUX*AINN(1)*(WZWAF(WAFER+1) - WZWAF(WAFER))
C      ELSE
C          CONVN(1) = 0.0
C          CONVS(1) = 0.0
C      ENDIF
C
C      IF (CFLOW(2)) THEN
C          CONVN(2) = QFLUX*AOUT(2)*(WZWAF(WAFER+1) - WZWAF(WAFER))
C          CONVS(2) = QFLUX*AINN(2)*(WZWAF(WAFER+1) - WZWAF(WAFER))
C      ELSE
C          CONVN(2) = 0.0
C          CONVS(2) = 0.0
C      ENDIF
C
C*****END FIX CONVECTIVE COOLING HEAT TRANSFER FLUXES

      IF (WALL .EQ. 1) QCONV = CONVS(1)
      IF (WALL .EQ. 2) QCONV = CONVS(2) - CONVN(1)

C*****PRINT QUARTZ WALL ENERGY BALANCES OVER WAFERS
CC      WRITE(TEXT, *) '------- QUARTZ TUBE ENERGY BALANCE ---------'
C      WRITE(TEXT, *) 'WALL = ', WALL, ' WAFER = ', WAFER
C      WRITE(TEXT, *) 'GEAST = ', GEAST, ' GWEST = ', GWEST
C      WRITE(TEXT, *) 'JEAST = ', JEAST, ' JWEST = ', JWEST
C      WRITE(TEXT, *) 'QNORTH = ', QNORTH
C      WRITE(TEXT, *) 'QSOUTH = ', QSOUTH
C      WRITE(TEXT, *) 'QNORTH TERM = ',
C     +            AREAW(WALL)*(WZWAF(WAFER+1) - WZWAF(WAFER))*QNORTH
C      WRITE(TEXT, *) 'QSOUTH TERM = ',
C     +            AREAW(WALL)*(WZWAF(WAFER+1) - WZWAF(WAFER))*QSOUTH
C      WRITE(TEXT, *) 'QEAST TERM = ',
C     +            PI*(ROUTQ(WALL)**2 - RINQ(WALL)**2)*QEAST
C      WRITE(TEXT, *) 'QWEST TERM = ',
C     +            PI*(ROUTQ(WALL)**2 - RINQ(WALL)**2)*QWEST
C      WRITE(TEXT, *) 'QCONV = ', QCONV
C      WRITE(TEXT, *) 'CONVN(1) = ',CONVN(1), ' CONVS(1) = ',CONVS(1)
C      WRITE(TEXT, *) 'CONVN(2) = ',CONVN(2), ' CONVS(2) = ',CONVS(2)
CC      WRITE(TEXT, *) '--------------------------------------------'
C*****END PRINT QUARTZ WALL ENERGY BALANCES OVER WAFERS


      FENG(WAFER, RPNTS + WALL) =
     +   AREAW(WALL)*(WZWAF(WAFER+1) - WZWAF(WAFER))*QNORTH -
     +   AREAW(WALL)*(WZWAF(WAFER+1) - WZWAF(WAFER))*QSOUTH +
     +   PI*(ROUTQ(WALL)**2 - RINQ(WALL)**2)*(QEAST - QWEST) + QCONV

C*****USE AREA BETWEEN WAFER FOR QUARTZ WALL 1
C      FENG(WAFER, RPNTS + WALL) =
C     +   AREAW(WALL)*(WZWAF(WAFER+1) - WZWAF(WAFER))*QNORTH -
C     +   AREAW(WALL)*SPACE*QSOUTH +
C     +   PI*(ROUTQ(WALL)**2 - RINQ(WALL)**2)*(QEAST - QWEST)
C*****END USE AREA BETWEEN WAFER FOR QUARTZ WALL 1

      IF (WALL .EQ. 1) THEN
         FENG(WAFER, RPNTS + WALL) = FENG(WAFER, RPNTS + WALL) - QTIP
         IF (WAFER .EQ. NUMWF-1) FENG(WAFER, RPNTS + WALL) =
     +                           FENG(WAFER, RPNTS + WALL) - QTIP2
      ELSE
      ENDIF

C     DIVIDED VOL*RHO*CP FOR THE ELEMENT.

      FENG(WAFER, RPNTS + WALL) =
     +   FENG(WAFER, RPNTS + WALL)/(PI*(ROUTQ(WALL)**2 -
     +   RINQ(WALL)**2)*(WZWAF(WAFER+1) -
     +   WZWAF(WAFER))*RHO(4)*CP(4))

C///  TRANSIENT TERM.

      IF (TIME) THEN
         FENG(WAFER, RPNTS + WALL) =
     +      FENG(WAFER, RPNTS + WALL) + TEX0(WAFER, WALL)
      END IF

C*****SET CONSTANT TEMPERATURE FOR QUARTZ WALL 1
C      IF (WALL .EQ. 1)
C     +   FENG(WAFER, RPNTS + WALL) = TEX(WAFER, WALL) - 800.0
C*****END SET CONSTANT TEMPERATURE FOR QUARTZ WALL 1

C*****SET CONSTANT TEMPERATURE FOR QUARTZ WALL 2
C      IF (WALL .EQ. 2)
C     +   FENG(WAFER, RPNTS + WALL) = TEX(WAFER, WALL) - 800.0
C*****END SET CONSTANT TEMPERATURE FOR QUARTZ WALL 2

C///  BOTTOM OF THE LOOP FOR ENERGY BALANCE OVER THE QUARTZ WALLS.

5070  CONTINUE

C///  ENERGY BALANCE FOR THE HEATING ELEMENT

      QNORTH = -COND*(TEX(WAFER, 4) - TEX(WAFER, 3))/(RIN - RTUBE)
      QSOUTH = J2(2, WAFER)*FH2 + J2(1, WAFER)*FHH - J2(1, WAFER)

C*****INSULATE BOTTOM OF HEATERS
C      QSOUTH = 0.0
C*****END INSULATE BOTTOM OF HEATERS

C*****INSULATE TOP OF HEATERS
C      QNORTH = 0.0
C*****END INSULATE TOP OF HEATERS

      IF (WAFER .EQ. 1) THEN
         QWEST = -CONDUC(2)*(TEX(WAFER, 3) -
     +            TX(PIN, 3, 1))/(ZWAF(WAFER) - Z(PIN, 1))
C*****INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
C         QWEST = 0.0
C*****END INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
      ELSE
         QWEST = -CONDUC(2)*(TEX(WAFER, 3) -
     +            TEX(WAFER-1, 3))/(ZWAF(WAFER) - ZWAF(WAFER-1))
      END IF

      IF (WAFER .EQ. (NUMWF-1)) THEN
         QEAST = -CONDUC(2)*(TEX(WAFER, 3) -
     +            TX(PEX, 3, 2))/((ZLEN - ZWAF(WAFER)) - Z(PEX, 2))
C        CHANGE SIGN BECAUSE OF CHANGE IN COORDINATE SYSTEM
         QEAST = - QEAST
C*****INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
C         QEAST = 0.0
C*****END INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
      ELSE
         QEAST = -CONDUC(2)*(TEX(WAFER+1, 3) -
     +            TEX(WAFER, 3))/(ZWAF(WAFER+1) - ZWAF(WAFER))
      END IF

C*****ZERO SIDE HEAT FLUX IN HEATERS
C      QEAST = 0.0
C      QWEST = 0.0
C*****END ZERO SIDE HEAT FLUX IN HEATERS

C*****PRINT HEATER ENERGY BALANCES
C      WRITE(TEXT, *)
C      WRITE(TEXT, *) '----------- HEATER ENERGY BALANCE ----------'
C      WRITE(TEXT, 568) WAFER,
C     +         2.*PI*RTUBEO*(WZWAF(WAFER+1) - WZWAF(WAFER))*QNORTH,
C     +         2.*PI*RTUBEI*(WZWAF(WAFER+1) - WZWAF(WAFER))*QSOUTH,
C     +         PI*(RTUBEO**2 - RTUBEI**2)*QEAST,
C     +         PI*(RTUBEO**2 - RTUBEI**2)*QWEST,
C     +         POWER(WAFER)
C      WRITE(TEXT, *) 'CONVN(2) = ', CONVN(2)
C      WRITE(TEXT, 569) QNORTH, COND,
C     +                 WZWAF(WAFER+1), WZWAF(WAFER), RTUBEO
C      WRITE(TEXT, *) '--------------------------------------------'
C  568 FORMAT('WAFER = ',I2,5(E11.4))
C  569 FORMAT(2X,'QNORTH = ',E11.4,1X,'COND = ',E11.4/
C     +       2X,'WZWAF(WAFER+1) = ',E11.4,1X,'WZWAF(WAFER) = ',
C     +       E11.4,1X,'RTUBEO = ',E11.4)
C*****END PRINT HEATER ENERGY BALANCES

      FENG(WAFER, RPNTS + 3) =
     +   2.*PI*RTUBEO*(WZWAF(WAFER+1) - WZWAF(WAFER))*QNORTH -
     +   2.*PI*RTUBEI*(WZWAF(WAFER+1) - WZWAF(WAFER))*QSOUTH +
     +   PI*(RTUBEO**2 - RTUBEI**2)*(QEAST - QWEST) -
     +   POWER(WAFER) - CONVN(2)

C*****SET CONSTANT TEMPERATURE FOR HEATERS
C      FENG(WAFER, RPNTS + 3) = TEX(WAFER, 3) - 800.0
C*****END SET CONSTANT TEMPERATURE FOR HEATERS

C     DIVIDED VOL*RHO*CP FOR THE ELEMENT.

      FENG(WAFER, RPNTS + 3) =
     +   FENG(WAFER, RPNTS + 3)/(PI*(RTUBEO**2 -
     +   RTUBEI**2)*(WZWAF(WAFER+1) - WZWAF(WAFER))*RHO(2)*CP(2))

C///  TRANSIENT TERM.

      IF (TIME) THEN
         FENG(WAFER, RPNTS + 3) =
     +      FENG(WAFER, RPNTS + 3) + TEX0(WAFER, 3)
      END IF

C///  ENERGY BALANCE FOR THE INSULATION ELEMENT.

      QSOUTH = QNORTH
      QNORTH = HCOEF(3)*(TEX(WAFER, 4) - TAMB) +
     +         EMISSI(3)*SIG*(TEX4(WAFER, 4) - TAMB4)

      IF (WAFER .EQ. 1) THEN
         QWEST = -CONDUC(3)*(TEX(WAFER, 4) -
     +            TX(PIN, 4, 1))/(ZWAF(WAFER) - Z(PIN, 1))
C*****INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
C         QWEST = 0.0
C*****END INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
      ELSE
         QWEST = -CONDUC(3)*(TEX(WAFER, 4) -
     +         TEX(WAFER-1, 4))/(ZWAF(WAFER) - ZWAF(WAFER-1))
      END IF

      IF (WAFER .EQ. (NUMWF-1)) THEN
         QEAST = -CONDUC(3)*(TEX(WAFER, 4) -
     +            TX(PEX, 4, 2))/((ZLEN - ZWAF(WAFER)) - Z(PEX, 2))
C        CHANGE SIGN BECUASE OF CHANGE IN COORDINATE SYSTEM
         QEAST = - QEAST
C*****INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
C         QEAST = 0.0
C*****END INSULATE MIDSECTION OF HEATERS AND INSULATION FROM ENDS
      ELSE
         QEAST = -CONDUC(3)*(TEX(WAFER+1, 4) -
     +            TEX(WAFER, 4))/(ZWAF(WAFER+1) - ZWAF(WAFER))
      END IF

C*****ZERO SIDE HEAT FLUXES IN INSULATION
C      QEAST = 0.0
C      QWEST = 0.0
C*****END ZERO SIDE HEAT FLUXES IN INSULATION

C*****ZERO REACTOR EXTERIOR HEAT LOSS OVER WAFER SECTION
C      QNORTH = 0.0
C*****END ZERO REACTOR EXTERIOR HEAT LOSS OVER WAFER SECTION

C*****PRINT INSULATION ENERGY BALANCES
C      WRITE(TEXT, *) '-------- INSULATION ENERGY BALANCE ---------'
C      WRITE(TEXT, 557) WAFER,
C     +             2.*PI*RINO*(WZWAF(WAFER+1) - WZWAF(WAFER))*QNORTH,
C     +             2.*PI*RINI*(WZWAF(WAFER+1) - WZWAF(WAFER))*QSOUTH,
C     +             PI*(RINO**2 - RINI**2)*QEAST,
C     +             PI*(RINO**2 - RINI**2)*QWEST
C      WRITE(TEXT, *) '--------------------------------------------'
C  557 FORMAT('INSU',1X,'W =',I2,4(E11.4))
C*****END PRINT INSULATION ENERGY BALANCES

      FENG(WAFER, RPNTS + 4) =
     +   2.*PI*RINO*(WZWAF(WAFER+1) - WZWAF(WAFER))*QNORTH -
     +   2.*PI*RINI*(WZWAF(WAFER+1) - WZWAF(WAFER))*QSOUTH +
     +   PI*(RINO**2 - RINI**2)*(QEAST - QWEST)

C     DIVIDED VOL*RHO*CP FOR THE ELEMENT.

      FENG(WAFER, RPNTS + 4) =
     +   FENG(WAFER, RPNTS + 4)/(PI*(RINO**2 -
     +   RINI**2)*(WZWAF(WAFER+1) -
     +   WZWAF(WAFER))*RHO(3)*CP(3))

      LOS1(WAFER) = 2.*PI*RINO*(WZWAF(WAFER+1) - WZWAF(WAFER))*QNORTH

C*****SET LOS1 AND LOS2
C      LOS1(WAFER) = 3.00
C*****END SET LOS1 AND LOS2

C*****SET CONSTANT INSULATION TEMPERATURE
C      FENG(WAFER, RPNTS + 4) = TEX(WAFER, 4) - 800.0
C*****END SET CONSTANT INSULATION TEMPERATURE

C///  TRANSIENT TERM.

      IF (TIME) THEN
         FENG(WAFER, RPNTS + 4) =
     +      FENG(WAFER, RPNTS + 4) + TEX0(WAFER, 4)
      END IF

C-----------------------AIR ENERGY BALANCE----------------------

      IF (NUMEX .EQ. 4) GO TO 5080

C///  ENERGY BALANCE FOR AIR BETWEEN HEATERS AND OUTER QUARTZ JAR
C     AND AIR BETWEEN TWO QUARTZ JARS.

C*****PRINT GAS ENERGY BALANCE TERMS
C
C      WRITE(TEXT, *) ' '
C      WRITE(TEXT, *) ' GAS COOLING ENERGY BALANCE TERMS'
C
C*****END PRINT GAS ENERGY BALANCE TERMS

      DO 5071 FLOW = 1, 2

C        CHECK FOR ENDS OF WAFER REGION.
         IF (WAFER .EQ. 1) THEN
            TWEST1 = TX(PIN, 4+FLOW, 1)
            TEAST1 = TEX(WAFER+1, 4+FLOW)
         ELSEIF (WAFER .GE. (NUMWF-1)) THEN
            TWEST1 = TEX(WAFER-1, 4+FLOW)
            TEAST1 = TX(PEX, 4+FLOW, 2)
         ELSE
            TWEST1 = TEX(WAFER-1, 4+FLOW)
            TEAST1 = TEX(WAFER+1, 4+FLOW)
         ENDIF

         IF (FLOW .EQ. 1) THEN
            MFLOWW = MFLOW(1)
            MFLOWE = MFLOW(1)
            EINJ = 0.0
         ELSEIF (FLOW .EQ. 2) THEN
            IF (RINJEC) THEN
               MFLOWW = MWZW(WAFER)
               MFLOWE = MWZW(WAFER+1)
               EINJ = MINJ(WAFER)*CPAIR*TTGAS(2)
            ELSE
               MFLOWW = MFLOW(2)
               MFLOWE = MFLOW(2)
               EINJ = 0.0
            ENDIF
         ELSE
         ENDIF

         EEAST = TEX(WAFER, 4+FLOW)*MAX(MFLOWE, 0.0)*CPAIR -
     +           TEAST1*MAX(-MFLOWE, 0.0)*CPAIR
         EWEST = TWEST1*MAX(MFLOWW, 0.0)*CPAIR -
     +           TEX(WAFER, 4+FLOW)*MAX(-MFLOWW, 0.0)*CPAIR

         FENG(WAFER, RPNTS+4+FLOW) = CONVN(FLOW) - CONVS(FLOW)
     +   + EEAST - EWEST - EINJ

C*****INCLUDE ONLY ADVECTIVE TERMS IN GAS ENERGY BALANCE
C
C        FENG(WAFER, RPNTS+4+FLOW) = EEAST - EWEST
C
C*****END INCLUDE ONLY ADVECTIVE TERMS IN GAS ENERGY BALANCE

C*****INCLUDE CONVS(1) ONLY FOR FLOW 1
C
C       IF ((FLOW.EQ.1).AND.(WAFER.EQ.2)) FENG(WAFER, RPNTS+4+FLOW) =
C     + EEAST - EWEST - CONVS(FLOW)
C
C*****END INCLUDE CONVS(1) ONLY FOR FLOW 1

C*****FIX GAS TEMPERATURES
C
C      FENG(WAFER, RPNTS+4+FLOW) = TEX(WAFER, 4+FLOW) - 450.0
C
C*****END FIX GAS TEMPERATURES

C*****PRINT GAS ENERGY BALANCE TERMS
C
C         WRITE(TEXT, *) '------- COOLING GAS ENERGY BALANCE ---------'
C         WRITE(TEXT, *) ' FLOW = ', FLOW,' WAFER = ', WAFER
C         WRITE(TEXT, *) ' MFLOWW = ', MFLOWW, ' MFLOWE = ', MFLOWE
C         WRITE(TEXT, *) ' MINJ = ', MINJ(WAFER)
C         WRITE(TEXT, *) ' EWEST = ',EWEST,' EEAST = ', EEAST
C         WRITE(TEXT, *) ' EINJ = ', EINJ
C         WRITE(TEXT, *) ' CONVN + CONVS + EEAST - EWEST - EINJ = ',
C     +                    CONVN(FLOW)-CONVS(FLOW)+EEAST-EWEST-EINJ
C         WRITE(TEXT, *) ' CONVN = ',CONVN(FLOW),' CONVS = ',CONVS(FLOW)
C         WRITE(TEXT, *) ' CONVN - CONVS = ', CONVN(FLOW)-CONVS(FLOW)
C         WRITE(TEXT, *) ' MFLOW = ',MFLOW(FLOW),' CPAIR = ',CPAIR
C         WRITE(TEXT, *) ' TEAST1 = ',TEAST1, ' TWEST1 = ',TWEST1
C         WRITE(TEXT, *) ' TGAS1 = ',TEX(WAFER,4+1),
C     +                  ' TGAS2 = ',TEX(WAFER,4+2)
C         IF((WAFER.GT.1).AND.(WAFER.LT.NUMWF))
C     +   WRITE(TEXT, *) ' TGAS1W = ',TEX(WAFER-1,4+1),
C     +                  ' TGAS1E = ',TEX(WAFER+1,4+1)
C         IF((WAFER.GT.1).AND.(WAFER.LT.NUMWF))
C     +   WRITE(TEXT, *) ' TGAS2W = ',TEX(WAFER-1,4+2),
C     +                  ' TGAS2E = ',TEX(WAFER+1,4+2)
C         WRITE(TEXT, *) ' TQUAR1 = ',TEX(WAFER,1),
C     +                  ' TQUAR2 = ',TEX(WAFER,2)
C         WRITE(TEXT, *) ' THEAT = ',TEX(WAFER,3)
C         WRITE(TEXT, *) '--------------------------------------------'
C
C*****END PRINT GAS ENERGY BALANCE TERMS

         TMASS = (VOL(FLOW)*(WZWAF(WAFER+1) -
     +   WZWAF(WAFER))*RHOAIR(TEX(WAFER, 4+FLOW))*CPAIR)

C        DIVIDED VOL*RHO*CP FOR THE ELEMENT.

         FENG(WAFER, RPNTS+4+FLOW) =
     +   FENG(WAFER, RPNTS+4+FLOW)/(VOL(FLOW)*(WZWAF(WAFER+1) -
     +   WZWAF(WAFER))*RHOAIR(TEX(WAFER, 4+FLOW))*CPAIR)

C*****SET THE GAS TEMPERATURE DERIVATIVE TO A CONSTANT.
C
C         FENG(WAFER, RPNTS+4+FLOW) = 0.0
C
C*****END SET THE GAS TEMPERATURE DERIVATIVE TO A CONSTANT.

C///  TRANSIENT TERM.

      IF (TIME) THEN
         FENG(WAFER, RPNTS+4+FLOW) =
     +      FENG(WAFER, RPNTS+4+FLOW) + TEX0(WAFER, 4+FLOW)
      END IF

5071  CONTINUE

C///  BOTTOM OF THE LOOP OVER WAFERS.

5080  CONTINUE

5090  CONTINUE

      IF (TFIX) THEN
         DO 5100 I = IBEG(3), IEND(3)
            FENG(I, RPNTS + 3) = TEX(I, 3) - TMID(I)
5100     CONTINUE
      END IF

C///  COPY EXTERNAL TEMPERATURES OVER WAFERS FROM NUMWF-1 TO NUMWF
C///  TO MAKE MATRIX SYMMETRIC

C*****COPY RESIDUAL FROM (NUMWF-1) TO NUMWF
C
C      DO 5110 I = 1, NUMEX
C         FENG(NUMWF, RPNTS+I) = FENG(NUMWF-1, RPNTS+I)
C         IF (TIME) THEN
C            FENG(WAFER, RPNTS+I) = FENG(NUMWF - 1, RPNTS+I)
C     +                            + TEX0(NUMWF, I)
C         END IF
C5110  CONTINUE
C
C*****END COPY RESIDUAL FROM (NUMWF-1) TO NUMWF

C*****COPY TEMPERATURES FROM (NUMWF-1) TO NUMWF
      DO 5120 I = 1, NUMEX
         FENG(NUMWF, RPNTS+I) = TEX(NUMWF, I) - TEX(NUMWF-1, I)
5120  CONTINUE
C*****END COPY TEMPERATURES FROM (NUMWF-1) TO NUMWF

C///////////////////////////////////////////////////////////////////////
C
C     (7) EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  THERMOCOUPLE ENERGY BALANCES.

      ABSTOP = 1.0
      ABSBOT = 1.0

      DO 7010 I = 1, NUMTC
C*****WRITE THERMOCOUPLE CORRECTION TERMS
C      WRITE(TEXT, *) ' '
C      WRITE(TEXT, *) 'TC NUMBER = ', I
C      WRITE(TEXT, *) 'SUMT = ', SUMT(I), ' SUMB = ', SUMB(I)
C      WRITE(TEXT, *) 'JTOP = ', JTOP(I), ' JBOT = ', JBOT(I)
C      WRITE(TEXT, *) 'UNCORR. GTOP = ', GTOP(I)
C      WRITE(TEXT, *) 'UNCORR. GBOT = ', GBOT(I)
C*****END WRITE THERMOCOUPLE CORRECTION TERMS
         IF (CORR) THEN
            GTOP(I) = GTOP(I) + (1.0 - SUMT(I))*JTOP(I)
            GBOT(I) = GBOT(I) + (1.0 - SUMB(I))*JBOT(I)
C*****WRITE THERMOCOUPLE CORRECTION TERMS
C      WRITE(TEXT, *) 'GTOP = ', GTOP(I)
C      WRITE(TEXT, *) 'GBOT = ', GBOT(I)
C*****END WRITE THERMOCOUPLE CORRECTION TERMS
         ENDIF
C         TC(I) = (0.5*(GTOP(I) + GBOT(I))/SIG)
         TC(I) = (ABSTOP*GTOP(I) + ABSBOT*GBOT(I))/(SIG*(ABSTOP +
     +           ABSBOT))
7010  CONTINUE

C///  PACK.

      COUNT = 0
      DO 7020 I = 1, PIN
CC       DO 7020 J = 1, NUMEX+1
         DO 7020 J = 1, NUMEX
            COUNT = COUNT + 1
            F(COUNT) = FENDS(I, J, 1)
7020  CONTINUE

CC    DO 7030 J = 1, 3
      DO 7030 J = 1, 2
         COUNT = COUNT + 1
         F(COUNT) = FENDS(PIN + 1, J, 1)
7030  CONTINUE

      DO 7050 J = 1, NUMWF
         DO 7040 K = 1, RPNTS
            COUNT = COUNT + 1
            F(COUNT) = FENG(J, K)
7040     CONTINUE

         DO 7050 I = 1, NUMEX
            COUNT = COUNT + 1
            F(COUNT) = FENG(J, RPNTS + I)
7050  CONTINUE

      DO 7060 I = PEX, 1, - 1
CC       DO 7060 J = 1, NUMEX+1
         DO 7060 J = 1, NUMEX
            COUNT = COUNT + 1
            F(COUNT) = FENDS(I, J, 2)
7060  CONTINUE

CC    DO 7070 J = 1, 3
      DO 7070 J = 1, 2
         COUNT = COUNT + 1
         F(COUNT) = FENDS(PEX + 1, J, 2)
7070  CONTINUE

      TOTAL = COUNT

C*****PRINT RESIDUAL
C      WRITE(TEXT, *) 'TIME = ', TIME
C      DO 7080 I = 1, TOTAL
C         WRITE(TEXT, *) 'F(',I,') = ',F(I)
C7080  CONTINUE
C*****END PRINT RESIDUAL


C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID, NUMTC
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID,
     +   ZLEN, WZWAF(1), WZWAF(NUMWF+1), ZTC(I), I
      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID

99101 FORMAT
     +   (/1X, A, 'ERROR. NUMBER OF THERMOCOUPLES MUST BE'
     +   /10X, 'GREATER THAN ZERO.'
     +  //10X, 'NUMTC = ',I5)

99102 FORMAT
     +   (/1X, A, 'ERROR. THERMOCOUPLES LOCATIONS IN ERROR.'
     +  //10X, 'REACTOR LENGTH = ', F7.2,
     +  //10X, 'POSITION OF FIRST WAFER = ', F7.2,
     +   /10X, 'POSITION OF LAST WAFER = ', F7.2,
     +   /10X, 'THERMOCOUPLE POSITION = ', F7.2,
     +   /10X, 'THERMOCOUPLE NUMBER = ', I7)

99103 FORMAT
     +   (/1X, A, 'ERROR. ERROR IN CALCULATION OF THERMOCOUPLE'
     +   /10X, 'SHAPEFACTORS.')

C///  EXIT.

99999 CONTINUE

      RETURN
      END
      SUBROUTINE TWAF2
     +  (ERROR, TEXT,
     +   CFLOW, COLMAX, COMPS, DATA4, GROUPA, GROUPB, HEADER, JEND,
     +   JEX, JIN,
     +   JWALL, KOUNT, LOS1, LOS2, LOS3, NADD,
     +   NUMEX, NUMTC, NUMWF, PDIM, PEX, PIN, PMAX, R, RINQ, ROUTQ,
     +   RPNTS, RTUBEI, RTUBEO, SUCCES, TC, TEX, THICKI, TIME, TIME2,
     +   TSTEP, TWAF, TX, VAL2, VAL3, VALUE, WZWAF, X, Z, ZWAF)

C///////////////////////////////////////////////////////////////////////
C
C     TWAF2
C
C     PRINT THE SOLUTION.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9, HEADER*(*), STRING*80
      INTEGER
     +   COL, COLMAX, COLS, COMPS, COUNT, DATA4, FIRST, GROUPA, GROUPB,
     +   I, J, K, KOUNT, LAST, LOC, MID, NADD, NUMEX, NUMTC, NUMWF,
     +   PDIM, PEX,
     +   PIN, PMAX, POINT, POINTS, RPNTS, TEXT, WIDTH
      LOGICAL BLANK, ERROR, CFLOW, SUCCES, TIME
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   DT, JEND, JEX, JIN, JWALL, LOS1, LOS2, LOS3, R, RINQ,
     +   ROUTQ, RTUBEI, RTUBEO,
     +   SUM1, SUM2, SUMT, TC, TEND, TENDI, TEX, THICKI, TIME2, TSTEP,
     +   TWAF, TWAF1, TX, VAL2, VAL3, VALUE, WZWAF, X, Z, ZWAF

      PARAMETER (ID = 'TWAF2:  ')
      PARAMETER (WIDTH = 12)

      DIMENSION
     +   CFLOW(2), HEADER(5, COLMAX), JEND(2), JEX(5*PEX + 2),
     +   JIN(5*PIN + 2),
     +   JWALL(PDIM, 2), LOS1(NUMWF),
     +   LOS2(PDIM, 2), LOS3(2), R(RPNTS), RINQ(2), ROUTQ(2), SUM2(2),
     +   TC(NUMTC), TEND(2), TENDI(2), TEX(NUMWF, NUMEX),
     +   TWAF(NUMWF+1, RPNTS), TX(PDIM, NUMEX, 2), VAL2(PDIM, COLMAX),
     +   VAL3(PMAX), VALUE(RPNTS+NUMEX, COLMAX), WZWAF(NUMWF + 1),
     +   X(GROUPA + COMPS*PMAX + GROUPB + NADD), Z(PDIM, 2),
     +   ZWAF(NUMWF)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  TOSHIBA/TEL YIELD STRESS DELTA T (T AND DELTA IN K).
C     TAIR - TEMPERATURE OF AIR IN K.

      DT(TWAF1) = 1.72*EXP(3500.0/TWAF1)

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (1 .LE. COLMAX)
      IF (ERROR) GO TO 9102

C///  UNPACK THE SOLUTION.

      COUNT = 0
      DO 1010 I = 1, PIN
CC         COUNT = COUNT + 1
CC         JWALL(I, 1) = X(COUNT)
      DO 1010 K = 1, NUMEX
         COUNT = COUNT + 1
         TX(I, K, 1) = X(COUNT)
1010  CONTINUE

CC      COUNT = COUNT + 1
CC      JEND(1) = X(COUNT)
      COUNT = COUNT + 1
      TEND(1) = X(COUNT)
      COUNT = COUNT + 1
      TENDI(1) = X(COUNT)

      DO 1020 J = 1, NUMWF
         DO 1015 K = 1, RPNTS
            COUNT = COUNT + 1
            TWAF(J, K) = X(COUNT)
1015     CONTINUE
         DO 1020 I = 1, NUMEX
            COUNT = COUNT + 1
            TEX(J, I) = X(COUNT)
1020  CONTINUE

      DO 1030 I = 1, PEX
CC         COUNT = COUNT + 1
CC         JWALL(PEX - I + 1, 2) = X(COUNT)
         DO 1030 K = 1, NUMEX
         COUNT = COUNT + 1
         TX(PEX - I + 1, K, 2) = X(COUNT)
1030  CONTINUE

CC      COUNT = COUNT + 1
CC      JEND(2) = X(COUNT)
      COUNT = COUNT + 1
      TEND(2) = X(COUNT)
      COUNT = COUNT + 1
      TENDI(2) = X(COUNT)

C///  UNPACK ADD ON VARIABLES

      DO 1040 I = 1, NADD
         COUNT = COUNT + 1
         VAL3(I) = X(COUNT)
1040  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     TOP OF THE BLOCK TO PREPARE ENTRANCE RADIOSITY VALUES.
C
C///////////////////////////////////////////////////////////////////////

      COL = 0
      POINTS = PIN

C///////////////////////////////////////////////////////////////////////
C
C     LOCATIONS OF RADIOSITY AND TEMPERATURES IN THE ENTRANCE.
C
C///////////////////////////////////////////////////////////////////////

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      HEADER(1, COL) = '       AXIAL'
      HEADER(2, COL) = '    POSITION'
      HEADER(3, COL) = '(CENTIMETERS'
      HEADER(4, COL) = '        FROM'
      HEADER(5, COL) = '   ENTRANCE)'

      DO 2010 POINT = 1, POINTS
         VAL2(POINT, COL) = Z(POINT, 1)
2010  CONTINUE

      END IF

C///////////////////////////////////////////////////////////////////////
C
C     RADIOSITIES AND TEMPERATURES IN THE ENTRANCE.
C
C///////////////////////////////////////////////////////////////////////

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      DO 3010 POINT = 1, POINTS
         VAL2(POINT, COL) = JWALL(POINT, 1)
3010  CONTINUE

      HEADER(1, COL) = '            '
      HEADER(2, COL) = '            '
      HEADER(3, COL) = '            '
      HEADER(4, COL) = '        WALL'
      HEADER(5, COL) = ' RADIOSITIES'

      END IF

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      DO 3015 POINT = 1, POINTS
         VAL2(POINT, COL) = TX(POINT, 1, 1)
3015  CONTINUE

      HEADER(1, COL) = '            '
      HEADER(2, COL) = '            '
      HEADER(3, COL) = '            '
      HEADER(4, COL) = ' QUARTZ WALL'
      HEADER(5, COL) = 'TEMPERATURES'

      END IF

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      DO 3020 POINT = 1, POINTS
         VAL2(POINT, COL) = TX(POINT, 2, 1)
3020  CONTINUE

      HEADER(1, COL) = '            '
      HEADER(2, COL) = '            '
      HEADER(3, COL) = '            '
      HEADER(4, COL) = ' QUARTZ WALL'
      HEADER(5, COL) = 'TEMPERATURES'

      END IF

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      DO 3025 POINT = 1, POINTS
         VAL2(POINT, COL) = TX(POINT, 3, 1)
3025  CONTINUE

      HEADER(1, COL) = '            '
      HEADER(2, COL) = '            '
      HEADER(3, COL) = '            '
      HEADER(4, COL) = '      HEATER'
      HEADER(5, COL) = 'TEMPERATURES'

      END IF

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      DO 3030 POINT = 1, POINTS
         VAL2(POINT, COL) = TX(POINT, 4, 1)
3030  CONTINUE

      HEADER(1, COL) = '            '
      HEADER(2, COL) = '            '
      HEADER(3, COL) = '        WALL'
      HEADER(4, COL) = '  INSULATION'
      HEADER(5, COL) = 'TEMPERATURES'

      END IF

      DO 3040 J = 1, 2
         IF (.NOT.CFLOW(J)) GO TO 3040
         COL = COL + 1
         IF (COL .LE. COLMAX) THEN

         DO 3035 POINT = 1, POINTS
            VAL2(POINT, COL) = TX(POINT, 4+J, 1)
3035     CONTINUE

         IF (CFLOW(1)) THEN
            HEADER(1, COL) = '            '
            HEADER(2, COL) = ' COOLING-GAS'
            HEADER(3, COL) = 'TEMPERATURES'
            HEADER(4, COL) = '     BETWEEN'
            HEADER(5, COL) = 'QUARTZ WALLS'
         ELSEIF (CFLOW(2)) THEN
            HEADER(1, COL) = ' COOLING-GAS'
            HEADER(2, COL) = 'TEMPERATURES'
            HEADER(3, COL) = '     BETWEEN'
            HEADER(4, COL) = '  HEATER AND'
            HEADER(5, COL) = ' QUARTZ WALL'
         ELSE
         END IF

         END IF
3040  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     BOTTOM OF THE BLOCK TO PREPARE THE COLUMNS FOR THE CELL VALUES.
C
C///////////////////////////////////////////////////////////////////////

      COLS = COL
      ERROR = .NOT. (COLS .LE. COLMAX)
      IF (ERROR) GO TO 9103

C///////////////////////////////////////////////////////////////////////
C
C     WRITE THE TIME IF A TRANSIENT RUN.
C
C///////////////////////////////////////////////////////////////////////

      MID = (NUMWF+1)/2

      IF (TIME) THEN
         WRITE(TEXT, 99109)
         WRITE(TEXT, 99108) TIME2, TIME2/60.
C         IF (KOUNT.EQ.0) WRITE(DATA4, 99112)
C99112    FORMAT(7X,'T(MIN)', ' T(SEC)',
C     +        ' TC',' TE', ' TQ1', ' TQ2', ' TH', ' TIN',
C     +        ' DT', ' ADT',
C     +        ' TEND1', ' TEND2')
C         WRITE(DATA4, 99107) TIME2/60., TIME2,
C     +        TWAF(MID, 1), TWAF(MID, RPNTS),
C     +        TEX(MID, 1), TEX(MID, 2), TEX(MID, 3), TEX(MID, 4),
C     +        ABS(TWAF(MID, 1) - TWAF(MID, RPNTS)), DT(TWAF(MID,1)),
C     +        TEND(1), TEND(2)
C         KOUNT = KOUNT + 1
      ELSE
      ENDIF

C///////////////////////////////////////////////////////////////////////
C
C     PRINT THE COLUMNS FOR THE RADIOSITY AND TEMPERATURE VALUES.
C
C///////////////////////////////////////////////////////////////////////

C*****WRITE ENTRANCE RADIOSITIES
C      WRITE (TEXT, *) '-- ENTRANCE RADIOSITIES --'
C      WRITE (TEXT, 99108) (JIN(K), K = 1, 5*PIN+2)
C      WRITE (TEXT, *) '-- EXIT RADIOSITIES --'
C      WRITE (TEXT, 99108) (JEX(K), K = 1, 5*PEX+2)
C99108 FORMAT(6(2X, F9.5))
C*****END WRITE ENTRANCE RADIOSITIES

      DO 4070 FIRST = 2, COLS, 4
         LAST = MIN (FIRST + 3, COLS)

         WRITE (TEXT, '()')
         WRITE (STRING, '(A, I2, A)')
     +      '(10X, A12,', LAST - FIRST + 1, '(2X, A12))'

         DO 4060 J = 1, 5
            BLANK = HEADER(J, 1) .EQ. ' '
            DO 4050 K = FIRST, LAST
               BLANK = BLANK .AND. (HEADER(J, K) .EQ. ' ')
4050        CONTINUE
            IF (.NOT. BLANK) WRITE (TEXT, STRING)
     +         HEADER(J, 1), (HEADER(J, K), K = FIRST, LAST)
4060     CONTINUE

         WRITE (STRING, '(A, I2, A)')
     +      '(/ 1P, (1X, I7, 2X, E12.4,',
     +      LAST - FIRST + 1, '(2X, E12.4)))'
         WRITE (TEXT, STRING)
     +      (J, VAL2(J, 1), (VAL2(J, K), K = FIRST, LAST),
     +      J = 1, POINTS)
4070  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     PRINT THE RADIOSITY AND TEMPERATURE FROM THE END.
C
C///////////////////////////////////////////////////////////////////////

      WRITE (TEXT, 99104) JEND(1)
      WRITE (TEXT, 99105) TEND(1)
      WRITE (TEXT, 99106) TENDI(1)

C///////////////////////////////////////////////////////////////////////
C
C     TOP OF THE BLOCK TO PREPARE FOR EXIT RADIOSITY VALUES.
C
C///////////////////////////////////////////////////////////////////////

      COL = 0
      POINTS = PEX

C///////////////////////////////////////////////////////////////////////
C
C     LOCATIONS OF RADIOSITIES AND TEMPERATURES IN THE EXIT.
C
C///////////////////////////////////////////////////////////////////////

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      HEADER(1, COL) = '       AXIAL'
      HEADER(2, COL) = '    POSITION'
      HEADER(3, COL) = '(CENTIMETERS'
      HEADER(4, COL) = '        FROM'
      HEADER(5, COL) = '       EXIT)'

      DO 5010 POINT = 1, POINTS
         VAL2(POINT, COL) = Z(POINT, 2)
5010  CONTINUE

      END IF

C///////////////////////////////////////////////////////////////////////
C
C     RADIOSITIES AND TEMPERATURES IN THE EXIT.
C
C///////////////////////////////////////////////////////////////////////

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      DO 6010 POINT = 1, POINTS
         VAL2(POINT, COL) = JWALL(POINT, 2)
6010  CONTINUE

      HEADER(1, COL) = '            '
      HEADER(2, COL) = '            '
      HEADER(3, COL) = '            '
      HEADER(4, COL) = '        WALL'
      HEADER(5, COL) = ' RADIOSITIES'

      END IF

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      DO 6015 POINT = 1, POINTS
         VAL2(POINT, COL) = TX(POINT, 1, 2)
6015  CONTINUE

      HEADER(1, COL) = '            '
      HEADER(2, COL) = '            '
      HEADER(3, COL) = '            '
      HEADER(4, COL) = ' QUARTZ WALL'
      HEADER(5, COL) = 'TEMPERATURES'

      END IF

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      DO 6020 POINT = 1, POINTS
         VAL2(POINT, COL) = TX(POINT, 2, 2)
6020  CONTINUE

      HEADER(1, COL) = '            '
      HEADER(2, COL) = '            '
      HEADER(3, COL) = '            '
      HEADER(4, COL) = ' QUARTZ WALL'
      HEADER(5, COL) = 'TEMPERATURES'

      END IF

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      DO 6025 POINT = 1, POINTS
         VAL2(POINT, COL) = TX(POINT, 3, 2)
6025  CONTINUE

      HEADER(1, COL) = '            '
      HEADER(2, COL) = '            '
      HEADER(3, COL) = '            '
      HEADER(4, COL) = '      HEATER'
      HEADER(5, COL) = 'TEMPERATURES'

      END IF

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      DO 6030 POINT = 1, POINTS
         VAL2(POINT, COL) = TX(POINT, 4, 2)
6030  CONTINUE

      HEADER(1, COL) = '            '
      HEADER(2, COL) = '            '
      HEADER(3, COL) = '        WALL'
      HEADER(4, COL) = '  INSULATION'
      HEADER(5, COL) = 'TEMPERATURES'

      END IF

      DO 6040 J = 1, 2
         IF (.NOT.CFLOW(J)) GO TO 6040
         COL = COL + 1
         IF (COL .LE. COLMAX) THEN

         DO 6035 POINT = 1, POINTS
            VAL2(POINT, COL) = TX(POINT, 4+J, 2)
6035     CONTINUE

         IF (CFLOW(1)) THEN
            HEADER(1, COL) = '            '
            HEADER(2, COL) = ' COOLING-GAS'
            HEADER(3, COL) = 'TEMPERATURES'
            HEADER(4, COL) = '     BETWEEN'
            HEADER(5, COL) = 'QUARTZ WALLS'
         ELSEIF (CFLOW(2)) THEN
            HEADER(1, COL) = ' COOLING-GAS'
            HEADER(2, COL) = 'TEMPERATURES'
            HEADER(3, COL) = '     BETWEEN'
            HEADER(4, COL) = '  HEATER AND'
            HEADER(5, COL) = ' QUARTZ WALL'
         ELSE
         END IF

         END IF
6040  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     BOTTOM OF THE BLOCK TO PREPARE THE COLUMNS FOR THE CELL VALUES.
C
C///////////////////////////////////////////////////////////////////////

      COLS = COL
      ERROR = .NOT. (COLS .LE. COLMAX)
      IF (ERROR) GO TO 9103

C///////////////////////////////////////////////////////////////////////
C
C     PRINT THE COLUMNS FOR THE RADIOSITY AND TEMPERATURE VALUES.
C
C///////////////////////////////////////////////////////////////////////

      DO 7070 FIRST = 2, COLS, 4
         LAST = MIN (FIRST + 3, COLS)

         WRITE (TEXT, '()')
         WRITE (STRING, '(A, I2, A)')
     +      '(10X, A12,', LAST - FIRST + 1, '(2X, A12))'

         DO 7060 J = 1, 5
            BLANK = HEADER(J, 1) .EQ. ' '
            DO 7050 K = FIRST, LAST
               BLANK = BLANK .AND. (HEADER(J, K) .EQ. ' ')
7050        CONTINUE
            IF (.NOT. BLANK) WRITE (TEXT, STRING)
     +         HEADER(J, 1), (HEADER(J, K), K = FIRST, LAST)
7060     CONTINUE

         WRITE (STRING, '(A, I2, A)')
     +      '(/ 1P, (1X, I7, 2X, E12.4,',
     +      LAST - FIRST + 1, '(2X, E12.4)))'
         WRITE (TEXT, STRING)
     +      (J, VAL2(J, 1), (VAL2(J, K), K = FIRST, LAST),
     +      J = 1, POINTS)
7070  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     PRINT THE RADIOSITY AND TEMPERATURE FROM THE END.
C
C///////////////////////////////////////////////////////////////////////

      WRITE (TEXT, 99104) JEND(2)
      WRITE (TEXT, 99105) TEND(2)
      WRITE (TEXT, 99106) TENDI(2)

C///////////////////////////////////////////////////////////////////////
C
C     TOP OF THE BLOCK TO PREPARE THE COLUMNS FOR WAFER TEMPERATURES.
C
C///////////////////////////////////////////////////////////////////////

      COL = 0
      POINTS = RPNTS + NUMEX

C///////////////////////////////////////////////////////////////////////
C
C     CELL CENTERS.
C
C///////////////////////////////////////////////////////////////////////

      COL = COL + 1
      IF (COL .LE. COLMAX) THEN

      HEADER(1, COL) = '      RADIAL'
      HEADER(2, COL) = '    POSITION'
      HEADER(3, COL) = '(CENTIMETERS'
      HEADER(4, COL) = '  FROM WAFER'
      HEADER(5, COL) = ' CENTERLINE)'

      DO 8010 POINT = 1, RPNTS
            VALUE(POINT, COL) = R(POINT)
8010  CONTINUE

      VALUE(RPNTS+1, COL) = (RINQ(1) + ROUTQ(1))/2.
      VALUE(RPNTS+2, COL) = (RINQ(2) + ROUTQ(2))/2.
      VALUE(RPNTS+3, COL) = (RTUBEO + RTUBEI)/2.
      VALUE(RPNTS+4, COL) = RTUBEO + THICKI/2.

      END IF

C///////////////////////////////////////////////////////////////////////
C
C     WAFER TEMPERATURES.
C
C///////////////////////////////////////////////////////////////////////

      IF (COL + NUMWF .LE. COLMAX) THEN

      DO 9010 POINT = 1, POINTS
         DO 9010 J = 1, NUMWF
            IF (POINT .LE. RPNTS) THEN
               VALUE(POINT, COL + J) = TWAF(J, POINT)
            ELSE
               VALUE(POINT, COL + J) = TEX(J, POINT - RPNTS)
            ENDIF
9010  CONTINUE

      DO 9030 J = 1, NUMWF
         HEADER(1, COL + J) = '       WAFER'
         HEADER(2, COL + J) = ' TEMPERATURE'
         HEADER(3, COL + J) = '          AT'
         WRITE(STRING, '(F12.2)') WZWAF(J)
         CALL RIGHTJ (STRING, HEADER(4, COL + J), WIDTH)
         HEADER(5, COL + J) = ' CENTIMETERS'

9030  CONTINUE

      END IF
      COL = COL + NUMWF

C///////////////////////////////////////////////////////////////////////
C
C     BOTTOM OF THE BLOCK TO PREPARE THE COLUMNS FOR THE CELL VALUES.
C
C///////////////////////////////////////////////////////////////////////

      COLS = COL
      ERROR = .NOT. (COLS .LE. COLMAX)
      IF (ERROR) GO TO 9103

C///////////////////////////////////////////////////////////////////////
C
C     PRINT THE COLUMNS FOR THE CELL VALUES.
C
C///////////////////////////////////////////////////////////////////////

      DO 10030 FIRST = 2, COLS, 4
         LAST = MIN (FIRST + 3, COLS)

         WRITE (TEXT, '()')
         WRITE (STRING, '(A, I2, A)')
     +      '(10X, A12,', LAST - FIRST + 1, '(2X, A12))'

         DO 10020 J = 1, 5
            BLANK = HEADER(J, 1) .EQ. ' '
            DO 10010 K = FIRST, LAST
               BLANK = BLANK .AND. (HEADER(J, K) .EQ. ' ')
10010       CONTINUE
            IF (.NOT. BLANK) WRITE (TEXT, STRING)
     +         HEADER(J, 1), (HEADER(J, K), K = FIRST, LAST)
10020    CONTINUE

         WRITE (STRING, '(A, I2, A)')
     +      '(/ 1P, (1X, I7, 2X, E12.4,',
     +      LAST - FIRST + 1, '(2X, E12.4)))'
         WRITE (TEXT, STRING)
     +      (J, VALUE(J, 1), (VALUE(J, K), K = FIRST, LAST),
     +      J = 1, POINTS)
10030 CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     PRINT ADDITIONAL VARIABLES.
C
C///////////////////////////////////////////////////////////////////////

      IF (NADD .NE. 0) THEN
         WRITE(TEXT, 99110)
         DO 9040 I = 1, NADD
            WRITE (TEXT, 99111) I, VAL3(I)
9040     CONTINUE
      ELSE
      ENDIF

C     PRINT EXTERIOR HEAT LOSS TERMS

C*****PRINT EXTERIOR HEAT LOSS TERMS
C      SUM1 = 0.0
C      SUM2(1) = 0.0
C      SUM2(2) = 0.0
C      SUMT = 0.0
C      WRITE(TEXT, *) ' '
C      WRITE (TEXT, *) 'ENTRANCE AND EXIT TERMS'
C      DO 10055 LOC = 1, 2
C         IF (LOC .EQ. 1) POINTS = PIN
C         IF (LOC .EQ. 2) POINTS = PEX
C         WRITE (TEXT, *) 'LOC = ', LOC, ' POINTS = ', POINTS
C         DO 10050 POINT = 1, POINTS
C            SUM2(LOC) = SUM2(LOC) + LOS2(POINT, LOC)
C            WRITE(TEXT, *) POINT, LOS2(POINT, LOC)
C10050 CONTINUE
C         WRITE(TEXT, *) 'SUM = ', SUM2(LOC)
C10055 CONTINUE
C      WRITE(TEXT, *) ' '
C      WRITE (TEXT, *) 'OVER THE WAFER LOAD'
C      DO 10060 POINT = 1, NUMWF
C         SUM1 = SUM1 + LOS1(POINT)
C         WRITE (TEXT, *) POINT, LOS1(POINT)
C10060 CONTINUE
C      WRITE (TEXT, *) 'SUM OVER WAFER LOAD = ', SUM1
C      WRITE (TEXT, *) ' '
C      WRITE (TEXT, *) 'ENTRANCE AND EXIT CAP'
C      WRITE (TEXT, *) '1', LOS3(1)
C      WRITE (TEXT, *) '2', LOS3(2)
C      WRITE (TEXT, *) ' '
C      SUMT = SUM1 + SUM2(1) + SUM2(2)
C      WRITE (TEXT, *) 'SUM OVER ENTRANCE, EXIT, AND WAFERS = ', SUMT
C      SUMT = SUMT - LOS3(1) - LOS3(2)
C      WRITE (TEXT, *) 'SUM TOTAL = ', SUMT
C*****END PRINT EXTERIOR HEAT LOSS TERMS

C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////


C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, COLMAX
      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID, COLS, COLMAX
      GO TO 99999

99102 FORMAT
     +   (/1X, A9, 'ERROR.  THE SPACE FOR COLUMNS MUST BE POSITIVE.'
     +  //10X, I10, '  COLMAX')

99103 FORMAT
     +   (/1X, A9, 'ERROR.  THE SPACE FOR COLUMNS IS TOO SMALL.'
     +  //10X, I10, '  COLS'
     +   /10X, I10, '  COLMAX')

99104 FORMAT
     +   (/12X,'JEND  = ',2X, 1PE12.4)

99105 FORMAT
     +   (12X, 'TEND  = ',2X, 1PE12.4)

99106 FORMAT
     +   (12X, 'TENDI = ',2X, 1PE12.4//)

99107 FORMAT
     +   (2X, F9.3, 2X, F9.3, 10(2X, F9.3))

99108 FORMAT
     +   (/12X,'TIME  = ',2X,1PE12.4,1X,'SEC',2X,'(',1PE12.4,1X,'MIN)')

99109 FORMAT
     +   (/7X,71('*'))

99110 FORMAT(/12X,'ADDITIONAL VARIABLES')
99111 FORMAT(/12X,'VARIABLE(',I2,') = ',1PE12.4)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWAF3
     +  (ERROR, TEXT,
     +   COMPS, COORD, FLAG, GROUPA, GROUPB, NUMEX, NUMWF,
     +   PIN, PMAX,
     +   R, REXT, RPNTS, TEMPER, TEX, TIME, TIME2, TWAF, UNIT, WZWAF,
     +   X, ZFIRST, ZLAST, ZWAF)

C///////////////////////////////////////////////////////////////////////
C
C     TWAF3
C
C     WRITE A PLOT2D FILE FOR THE SOLUTION.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9, LEGE*80, NAME*80, STRING*80, SUBT*80,
     +   TITL*80, XLAB*80, YLAB*80, ZLAB*80
      INTEGER
     +   COUNT, COMPS, CURVES, GROUPA, GROUPB, I, J, K, LENGTH, NUMBER,
     +   NUMEX, NUMWF, OFFSET, PIN, PMAX, RPNTS, STATUS, TEXT, UNIT,
     +   XINC,
     +   XSIZE, YINC, YSIZE
      LOGICAL ERROR, CFLOW, FLAG, TIME
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   COORD, R, REXT, TEMPER, TEX, TIME2, TWAF, WZWAF, X, ZFIRST,
     +   ZLAST, ZWAF

      PARAMETER (ID = 'TWAF3:  ')

      DIMENSION
     +   CFLOW(2), COORD(NUMWF), R(RPNTS), REXT(NUMEX),
     +   TEMPER(RPNTS+NUMEX, NUMWF),
     +   TEX(NUMWF, NUMEX), TWAF(NUMWF+1, RPNTS), WZWAF(NUMWF+1),
     +   X(GROUPA + COMPS*PMAX + GROUPB),
     +   ZWAF(NUMWF)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////


C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (1 .LE. NUMWF .AND. 0 .LE. RPNTS)
      IF (ERROR) GO TO 9102

C///  UNPACK THE SOLUTION.

CC      OFFSET = (NUMEX+1)*PIN + 3
      OFFSET = (NUMEX)*PIN + 2

      COUNT = OFFSET
      DO 1015 J = 1, NUMWF
         DO 1010 K = 1, RPNTS
            COUNT = COUNT + 1
            TWAF(J, K) = X(COUNT)
            TEMPER(K, J) = TWAF(J, K)
1010     CONTINUE
         DO 1015 I = 1, NUMEX
            COUNT = COUNT + 1
            TEX(J, I) = X(COUNT)
            TEMPER(I+RPNTS, J) = TEX(J, I)
1015  CONTINUE

      IF (FLAG) THEN

C///  WRITE PLOT2D HEADER.

         REWIND UNIT

         STRING = 'PLOT2D VERSION 1.09'
         CALL EXTENT (LENGTH, STRING)
         WRITE (UNIT, '(1X, A)', ERR = 9103, IOSTAT = STATUS)
     +         STRING (1 : LENGTH)
         FLAG = .FALSE.

      ELSE
      ENDIF

C///  OFFSET COORDINATES FROM FIRST WAFER.

      DO 1020 J = 1, NUMWF
         COORD(J) = WZWAF(J) - ZFIRST
1020  CONTINUE

C///  INITIALIZE THE COUNT OF CURVES.

      CURVES = 0

C///////////////////////////////////////////////////////////////////////
C
C     WAFER TEMPERATURES.
C
C///////////////////////////////////////////////////////////////////////

C     SUBROUTINE WRT2D1
C    +  (ERROR, TEXT,
C    +   DATA, LEGE, NAME, SUBT, TITL, XCOORD, XINC, XLAB, XSIZE,
C    +   YCOORD, YINC, YLAB, YSIZE, ZCOORD, ZLAB)

      CURVES = CURVES + 1

      NAME = 'WAFER TEMPERATURES'
      TITL = NAME
      SUBT = NAME
      XINC = 1
      YINC = RPNTS + NUMEX
      XLAB = 'RADIAL POSITION (CM)'
      YLAB = 'AXIAL POSITION (CM)'
      ZLAB = 'WAFER TEMPERATURE (K) '
      XSIZE = RPNTS
      YSIZE = NUMWF

C///  IF TRANSIENT SIMULATION SET LEGEND TO TIME.

      IF (TIME) THEN
         WRITE (LEGE, '(A, 1X, F10.2)') 'TIME = ',TIME2
         CALL SQUEEZ (LENGTH, LEGE)
      ELSE
         LEGE = NAME
      ENDIF

C*****DOUBLE PRECISION
         CALL WRT2D2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C         CALL WRT2D1
C*****END SINGLE PRECISION
     +     (ERROR, TEXT,
     +      UNIT, LEGE, NAME, SUBT, TITL, R, XINC, XLAB, XSIZE,
     +      COORD, YINC, YLAB, YSIZE, TEMPER, ZLAB)
      IF (ERROR) GO TO 9201

C///////////////////////////////////////////////////////////////////////
C
C     EXTERNAL TEMPERATURES.
C
C///////////////////////////////////////////////////////////////////////

C///  DON'T WRITE OUT GAS TEMPERATURES

      NUMBER = 4

C///  OFFSET COORDINATES FROM FIRST WAFER.

      DO 2010 J = 1, NUMWF-1
         COORD(J) = ZWAF(J) - ZFIRST
2010  CONTINUE

      DO 2020 J = 1, NUMWF-1
         DO 2020 K = 1, NUMBER
            TEMPER(K, J) = TEX(J, K)
2020  CONTINUE

C     SUBROUTINE WRT2D1
C    +  (ERROR, TEXT,
C    +   DATA, LEGE, NAME, SUBT, TITL, XCOORD, XINC, XLAB, XSIZE,
C    +   YCOORD, YINC, YLAB, YSIZE, ZCOORD, ZLAB)

      CURVES = CURVES + 1

      NAME = 'EXTERNAL TEMPERATURES'
      TITL = NAME
      SUBT = NAME
      XINC = 1
      YINC = RPNTS + NUMEX
      XLAB = 'RADIAL POSITION (CM)'
      YLAB = 'AXIAL POSITION (CM)'
      ZLAB = 'TEMPERATURE (K) '
      XSIZE = NUMBER
      YSIZE = NUMWF-1

C///  IF TRANSIENT SIMULATION SET LEGEND TO TIME.

      IF (TIME) THEN
         WRITE (LEGE, '(A, 1X, F10.2)') 'TIME = ',TIME2
         CALL SQUEEZ (LENGTH, LEGE)
      ELSE
         LEGE = NAME
      ENDIF

C*****DOUBLE PRECISION
         CALL WRT2D2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C         CALL WRT2D1
C*****END SINGLE PRECISION
     +     (ERROR, TEXT,
     +      UNIT, LEGE, NAME, SUBT, TITL, REXT, XINC, XLAB, XSIZE,
     +      COORD, YINC, YLAB, YSIZE, TEMPER, ZLAB)
      IF (ERROR) GO TO 9201

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, NUMWF, RPNTS
      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID, UNIT, STATUS
      GO TO 99999

9201  IF (0 .LT. TEXT) WRITE (TEXT, 99201) ID, CURVES
      GO TO 99999

99102 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +  //10X, I10, '  NUMWF'
     +   /10X, I10, '  RPNTS')

99103 FORMAT
     +   (/1X, A9, 'ERROR.  FORMATTED WRITE FILE FAILS.'
     +  //10X, I10, '  UNIT NUMBER'
     +   /10X, I10, '  I/O STATUS')

99201 FORMAT
     +   (/1X, A9, 'ERROR.  WRT2D1 OR WRT2D2 FAILS TO WRITE A PLOT2D'
     +   /10X, 'RECORD.'
     +   //10X, I10, '  CURVE')

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWAF4
     +  (ERROR, TEXT,
     +   NUMEX, NUMWF, PDIM, PEX, PIN, R, RPNTS,
     +   TEX, TOUT, TWAF, TX, UNIT, WZWAF, Z, ZFIRST,
     +   ZLAST, ZLEN, ZOUT, ZWAF)

C///////////////////////////////////////////////////////////////////////
C
C     TWAF4
C
C     WRITE A TEMPERATURE FILE TO BE READ BY OVEND.
C
C///////////////////////////////////////////////////////////////////////
C
C     DESCRIPTION OF THE SUBROUTINE ARGUMENTS:
C
C     ERROR   OUTPUT LOGICAL - ERROR FLAG.  IF TRUE, THEN AN ERROR
C             BLOCKS EXECUTION.  ERROR MESSAGES APPEAR IN THE OUTPUT
C             TEXT FILE.
C
C     TEXT    INPUT INTEGER - UNIT NUMBER FOR AN OUTPUT FILE.  ZERO AND
C             NEGATIVE VALUES FOR "TEXT" SUSPEND OUTPUT.
C
C     INTERNAL:
C
C     TOUT    REAL DIMENSIONED (2*PDIM + NUMWF + 2) - ARRAY FOR STORING
C             TEMPERATURES FOR OUTPUT.
C
C     ZOUT    REAL DIMENSIONED (2*PDIM + NUMWF + 2) - ARRAY FOR STORING
C             AXIAL LOCATIONS FOR OUTPUT.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   COUNT, I, J, NUMEX, NUMWF,
     +   PIN, PDIM, PEX, RPNTS, STATUS, TEXT, UNIT
      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   R, TEX, TOUT, TWAF, TX, WZWAF, Z, ZFIRST, ZLAST,
     +   ZLEN, ZOUT, ZWAF

      PARAMETER (ID = 'TWAF4:  ')

      DIMENSION
     +   R(RPNTS), TEX(NUMWF, NUMEX), TOUT(2*PDIM + NUMWF + 2),
     +   TWAF(NUMWF+1, RPNTS), TX(PDIM, NUMEX, 2), WZWAF(NUMWF+1),
     +   Z(PDIM, 2), ZOUT(2*PDIM + NUMWF+ 2), ZWAF(NUMWF)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////


C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (1 .LE. NUMWF .AND. 0 .LE. RPNTS)
      IF (ERROR) GO TO 9102

C///  FILL ZOUT AND TOUT ARRAYS

      COUNT = 0

C///  THE FIRST POINT IS Z = 0.0 (USE TEMPERATURE AT NEAREST NODE)

      COUNT = COUNT + 1
      ZOUT(COUNT) = 0.0
      TOUT(COUNT) = TX(1, 1, 1)

      DO 100 I = 1, PIN
         COUNT = COUNT + 1
         ZOUT(COUNT) = Z(I, 1)
         TOUT(COUNT) = TX(I, 1, 1)
100   CONTINUE
      DO 200 I = 1, NUMWF
         COUNT = COUNT + 1
         ZOUT(COUNT) = WZWAF(I)
         TOUT(COUNT) = TEX(I, 1)
200   CONTINUE
      DO 300 I = 1, PEX
         COUNT = COUNT + 1
         ZOUT(COUNT) = ZLEN - Z(PEX - I + 1, 2)
         TOUT(COUNT) = TX(PEX - I + 1, 1, 2)
300   CONTINUE

C///  THE LAST POINT IS REACTOR LENGTH (USE TEMPERATURE AT NEAREST NODE)

      COUNT = COUNT + 1
      ZOUT(COUNT) = ZLEN
      TOUT(COUNT) = TX(1, 1, 2)

C///////////////////////////////////////////////////////////////////////
C
C     REACTOR WALL TEMPERATURES.
C
C///////////////////////////////////////////////////////////////////////

       WRITE (UNIT, 99201, ERR = 9103, IOSTAT = STATUS)
     +       COUNT

       WRITE (UNIT, 99202, ERR = 9103, IOSTAT = STATUS)
     +       (ZOUT(I), TOUT(I), I = 1, COUNT)

C///////////////////////////////////////////////////////////////////////
C
C     WAFER TEMPERATURES.
C
C///////////////////////////////////////////////////////////////////////


      WRITE (UNIT, 99301, ERR = 9103, IOSTAT = STATUS)
     +      NUMWF, RPNTS

      WRITE (UNIT, 99302, ERR = 9103, IOSTAT = STATUS)
     +      (WZWAF(I), I = 1, NUMWF)

      WRITE (UNIT, 99302, ERR = 9103, IOSTAT = STATUS)
     +      (R(I), I = 1, RPNTS)

      WRITE (UNIT, 99302, ERR = 9103, IOSTAT = STATUS)
     +      ((TWAF(I, J), J = 1, RPNTS), I = 1, NUMWF)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID, NUMWF, RPNTS
      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID, UNIT, STATUS
      GO TO 99999

99102 FORMAT
     +   (/1X, A9, 'ERROR.  SOME DIMENSIONS ARE NOT POSITIVE.'
     +  //10X, I10, '  NUMWF'
     +   /10X, I10, '  RPNTS')

99103 FORMAT
     +   (/1X, A9, 'ERROR.  FORMATTED WRITE FILE FAILS.'
     +  //10X, I10, '  UNIT NUMBER'
     +   /10X, I10, '  I/O STATUS')

99201 FORMAT(I3)
99202 FORMAT(E14.7, 1X, E14.7)

99301 FORMAT(2I3)
99302 FORMAT(E14.7)

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWAF7
     +  (ERROR, TEXT,
     +   LOS1, LOS2, LOS3, MAXHT, NUMHT, NUMWF, PCAP, PDIM, PEX, PIN,
     +   POWERH)

C///////////////////////////////////////////////////////////////////////
C
C     TWAF7
C
C     ENERGY BALANCE CHECK.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   I, LOC, MAXHT, NUMHT, NUMWF, PDIM, PEX, PIN, POINT, POINTS,
     +   TEXT
      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   LOS1, LOS2, LOS3, PCAP, POWERH, SUM1, SUM2, SUM3, SUMT

      PARAMETER (ID = 'TWAF7:  ')

      DIMENSION
     +   LOS1(NUMWF), LOS2(PDIM, 2), LOS3(2), PCAP(2), POWERH(MAXHT),
     +   SUM2(2)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C     SUM HEATER POWERS.

      SUM1 = 0.0
      SUM2(1) = 0.0
      SUM2(2) = 0.0
      SUM3 = 0.0
      DO 1010 I = 1, NUMHT
         SUM1 = SUM1 + POWERH(I)
1010  CONTINUE

C     ADD POWER IN END CAP REGIONS.

      SUM1 = SUM1 + PCAP(1) + PCAP(2)

C     SUM EXTERIOR LOSSES OVER ENTRANCE AND EXIT.

      DO 1030 LOC = 1, 2
         IF(LOC .EQ. 1) POINTS = PIN
         IF(LOC .EQ. 2) POINTS = PEX
         DO 1020 POINT = 1, POINTS
            SUM2(LOC) = SUM2(LOC) + LOS2(POINT, LOC)
1020     CONTINUE
1030  CONTINUE

C     SUM EXTERIOR LOSSES OVER WAFERS

      DO 1040 I = 1, NUMWF-1
         SUM3 = SUM3 + LOS1(I)
1040  CONTINUE

C     ADD TOGETHER LOSSES INCLUDING LOSSES FROM ENDS.

      SUMT = SUM2(1) + SUM2(2) + SUM3 - LOS3(1) - LOS3(2)

C     PERCENTAGE ERROR.

      WRITE(TEXT, *) ' '
      WRITE(TEXT, *) 'SUM OF HEATER POWERS = ', SUM1
      WRITE(TEXT, *) ' '
      WRITE(TEXT, *) 'TOTAL HEAT LOSS FROM REACTOR EXTERIOR = ', SUMT
      WRITE(TEXT, *) ' '

C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////


C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      SUBROUTINE TWAF8
     +  (ERROR, TEXT,
     +   COMPS, DATA8, IRW, MAXHT, NADD, NEQ, NUMHT, PCAP, POWERH, X)

C///////////////////////////////////////////////////////////////////////
C
C     TWAF8
C
C     SUBROUTINE FOR READING AND WRITING RESTART FILE.
C
C     IRW = 0 WRITE VARIABLES.
C     IRW = 1 READ VARIABLES.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER ID*9
      INTEGER
     +   DATA8,
     +   I, IRW, MAXHT, NADD, NEQ, NUMHT, TEXT
      LOGICAL ERROR
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   PCAP, POWERH, X

      PARAMETER (ID = 'TWAF8:  ')

      DIMENSION
     +   PCAP(2), POWERH(MAXHT), X(NEQ)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///////////////////////////////////////////////////////////////////////
C
C     (1) READ OR WRITE VARIABLES.
C
C///////////////////////////////////////////////////////////////////////

      IF (IRW .EQ. 0) THEN
         WRITE (DATA8, 10001) (X(I), I = 1, NEQ)
         WRITE (DATA8, 10001) (POWERH(I), I = 1, NUMHT)
         WRITE (DATA8, 10001) (PCAP(I), I = 1, 2)
      ELSEIF (IRW .EQ. 1) THEN
         READ (DATA8, 10001) (X(I), I = 1, NEQ)
         READ (DATA8, 10001) (POWERH(I), I = 1, NUMHT)
         READ (DATA8, 10001) (PCAP(I), I = 1, 2)
      ELSE
      ENDIF
10001 FORMAT (1X, E20.13)

C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////


C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

C///  EXIT.

99999 CONTINUE
      RETURN
      END
      PROGRAM TWAFER

C///////////////////////////////////////////////////////////////////////
C
C     TWAFER
C
C     MODEL THE RADIATION EXCHANGE AND HEAT TRANSFER IN A HIGH
C     TEMPERATURE MULTIPLE-WAFERS-IN-TUBE LPCVD FURNANCE.
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER
     +   BANNER*80, BLANK*1, BNAME*80, CKEY*32, CWORK*16,
     +   ID*9, KEY*80, LINE*82, REPORT*80, SIGNAL*16, STRING*80,
     +   TYPE*80, VERSIO*80, VMSS*80, WORD*80
      INTEGER
     +   A2DIM, AMAX, BLOCKS, CHOICE, CLAST, CMAX,
     +   COLMAX, COMPS, COUNT, CSIZE, DATA1, DATA2, DATA3, DATA4,
     +   GROUPA,
     +   GROUPB, I, IDID, IKEY, ILAST, IMAX, IN, INFO,
     +   IRES, ISIZE,
     +   IWORK, J, K, KMAX, KOUNT, LENGTH, LEVEL,
     +   LIW, LLAST, LMAX, LRW, LSIZE, MAXHT, MAXIN,
     +   MAXORD, ML, MU, NADD, NG, NEQ, NUMBER, NUMEX, NUMHT,
     +   NUMIN, NUMT, NUMTC, NUMWF, OFFSET, PDIM, PEX, PIN, PMAX,
     +   POINTS, PRINT, RLAST, RMAX, RPNTS,
     +   RSIZE, SCRIPT, SSAGE,
     +   STEPS0, STEPS1, STEPS2, T, TDAGE, TEMMAX, TEXT, VMSSS, W,
     +   XISIZE, XRSIZE
      LOGICAL
     +   ADAPT, CFLOW, BFOUND, BNEED, ERROR, FLAG, FLAG2,
     +   FLAG3, FOUND, LKEY,
     +   LWORK, MATRIX, NEED, PCON, PED, RBOAT, RINJEC, SHAPE,
     +   SUCCES, TCAPS, TFIX, TIME
C*****DOUBLE PRECISION
      DOUBLE PRECISION
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      REAL
C*****END SINGLE PRECISION
     +   ATOL, CONDIT, CONDUC, CP, DTIME, DZEND, DZENDI,
     +   EMIS, EMISS, EMISSE,
     +   EMISSI, EMISSW, HCOEF, MFLOW, PCAP,
     +   PEDRI, PEDRO,
     +   R, RADIUS, REFL, RHO, RINQ,
     +   RKEY, ROUTQ, RTOL,
     +   RTUBEI, RTUBEO, SIG, SPACE, SSABS, SSREL, TAMB,
     +   TCAP, TDABS, TDEC, TDREL, THICK, THICKI, THIGH, TIME1,
     +   TIME2, TINC,
     +   TLOW, TMAX, TMIN, TOLER0, TOLER1, TOLER2, TOUT, TRAN, TSTEP,
     +   TSTEP0, TSTOP,
     +   TTEND, TTGAS, TTIN, TTHIC,
     +   TTQ, TTUBE, TTWAF, VFLOW, VINJ,
     +   ZBEGI, ZENDI, ZFIRST, ZLAST, ZLEN

      PARAMETER (ID = 'TWAFER:  ')
      PARAMETER
     +  (CSIZE = 1 000, ISIZE = 3 000,
     +   LSIZE = 1 000, RSIZE = 5 000 000)
      PARAMETER
     +  (DATA1 = 11, DATA2 = 12, DATA3 = 13, SCRIPT = 14, TEXT = 15,
     +   DATA4 = 16)
      PARAMETER (BLANK = ' ')
      PARAMETER (PRINT = 1)
      PARAMETER (KMAX = 60)
      PARAMETER (MAXHT = 10, NUMTC = 5, TEMMAX = 10)
      PARAMETER (MAXIN = 20)
      PARAMETER (W = 1, T = 2, IN = 3)
      PARAMETER
     +   (QGRID = 1, QREAC = 2, QTEMP = 3, QTWO = 4, QINIT = 5, QHT = 6,
     +    QTC = 7, QENDT = 8, QENDP = 9, QDAS = 10, QSS = 11,
     +    QTRANS =12, QRBOAT=13, BLOCKS = 13)
      PARAMETER (VMSSS = 5)
      PARAMETER (SIG = 5.6696E-12)

      DIMENSION
     +   BANNER(2, 10), BFOUND(BLOCKS), BNAME(BLOCKS), BNEED(BLOCKS),
     +   CFLOW(2), CKEY(KMAX), CONDUC(8), CP(8),
     +   CWORK(CSIZE), DZEND(2), DZENDI(2), EMIS(2), EMISSI(5),
     +   FOUND(KMAX), HCOEF(5),
     +   IKEY(KMAX), INFO(15), IWORK(ISIZE),
     +   KEY(KMAX, BLOCKS),
     +   LEVEL(2), LKEY(KMAX), LWORK(LSIZE),
     +   MFLOW(2), NEED(KMAX, BLOCKS), PCAP(2), VFLOW(2), R(RSIZE),
     +   REFL(2), RHO(8),
     +   RINQ(2), RKEY(KMAX), ROUTQ(2),
     +   TCAP(2), TCAPS(2), TRAN(2), TTGAS(2), TYPE(KMAX, BLOCKS),
     +   VINJ(MAXIN), VMSS(VMSSS), ZBEGI(MAXIN), ZENDI(MAXIN)

C///  POINTERS FOR RPAR AND IPAR FOR DASSL

      COMMON /RPOINT/
     +   QA2,    QAEX,   QAIN,   QAREAC, QTAMB,  QCOND,
     +   QCONW,  QCP,    QDZEND, QDZNDI, QEMIS,  QEMISS,
     +   QEMISE, QEMISI, QEMISW, QF,     QFEE,
     +   QFEG,   QFENDS, QFENG,  QFESE,  QFEW,
     +   QFGE,   QFGSE,  QFSEE,  QFSEG,  QFSES,
     +   QFSEW,  QFSS,   QFSW,   QFWE,   QFWS,   QFWSE,
     +   QFWW,   QHCOEF, QJ1,    QJ2,    QJEND,  QJEX,
     +   QJIN,   QJWALL, QLOS1,  QLOS2,  QLOS3,  QMFLO,
     +   QMINJ,  QMINJE, QMWZ,   QMWZW,
     +   QPCAP,  QPOWER, QPOWRH, QPWALL, QQWAF,
     +   QR,     QRCPV,  QREFL,  QRHO,   QRING,  QRINQ,
     +   QROUTQ,
     +   QRTUBI, QRTUBO, QBUF,   QX0,    QSPACE,
     +   QTCAP,  QTDAT,  QTEAST, QTEX,
     +   QTEX0,  QTEX4,  QTHICK, QTHCKI,
     +   QTMID,  QTRAN,  QTSIDE, QTSTEP, QTTGAS, QTWAF,
     +   QTWAF0, QTWAF4, QTWEST, QTX,    QTX0,
     +   QTX4,   QWR,    QWZ,    QWZWAF, QZ,
     +   QZBEG,  QZDAT,  QZEND,  QZLEN,  QZWAF,
     +   QGBOT,  QGTOP,  QJBOT,  QJTOP,  QSUMB,
     +   QSUMT,  QTCT,   QVINJ,  QZBEGI, QZENDI, QZTC

      COMMON /IPOINT/
     +   QA2DIM, QCHOIC, QCOMPS, QGROPA,
     +   QGROPB, QIBEG,  QIEND,  QIPV2,  QIPV3,  QIPV4,
     +   QMAXHT, QMAXIN, QNADD,  QNUMEX, QNUMHT, QNUMIN,
     +   QNUMT,  QNUMTC,
     +   QNUMWF, QPDIM,  QPEX,   QPIN,   QPMAX,
     +   QRPNTS, QSHAPE, QTEXT,  QTEMAX, QTFIX,
     +   QTIME

      COMMON /LOGIC/
     +   ERROR, CFLOW, RBOAT, RINJEC, SHAPE, TCAPS, TFIX, TIME

      COMMON /PED1/ PEDRI, PEDRO, TTHIC
      COMMON /PED2/ PCON, PED

      EXTERNAL RES
C*****DOUBLE PRECISION DASSRT
C      EXTERNAL ROOT
C*****END DOUBLE PRECISION DASSRT

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  OPEN UNICOS FILES.

C*****UNICOS: OPEN STATEMENTS
C      OPEN (FILE = 'twafer.out', FORM = 'FORMATTED',
C     +   STATUS = 'UNKNOWN', UNIT = TEXT)
C      OPEN (FILE = 'plot2d_twafer.dat', FORM = 'FORMATTED',
C     +   STATUS = 'UNKNOWN', UNIT = DATA1)
C      OPEN (FILE = 'temp.dat', FORM = 'FORMATTED',
C     +   STATUS = 'UNKNOWN', UNIT = DATA2)
C      OPEN (FILE = 'twafer.in', FORM = 'FORMATTED',
C     +   STATUS = 'OLD', UNIT = SCRIPT)
C*****END UNICOS: OPEN STATEMENTS

C///  OPEN VAX FILES.

C*****VAX: OPEN STATEMENTS
C      OPEN (FILE = 'PLOT2D_TWAFER.DAT',
C     +   FORM = 'FORMATTED', STATUS = 'NEW', UNIT = DATA1)
C      OPEN (FILE = 'TEMP.DAT',
C     +   FORM = 'FORMATTED', STATUS = 'NEW', UNIT = DATA2)
C      OPEN (FILE = 'TWAFER.IN',
C     +   FORM = 'FORMATTED', STATUS = 'OLD', UNIT = SCRIPT)
C*****END VAX: OPEN STATEMENTS

C*****VAX: OUTPUT TO SYS$OUTPUT
C      OPEN (FILE = 'SYS$OUTPUT',
C     +   FORM = 'FORMATTED', STATUS = 'OLD', UNIT = TEXT)
C*****END VAX: OUTPUT TO SYS$OUTPUT

C*****VAX: OUTPUT TO TWAFER.OUT
C      OPEN (FILE = 'TWAFER.OUT',
C     +   FORM = 'FORMATTED', STATUS = 'NEW', UNIT = TEXT)
C*****END VAX: OUTPUT TO TWAFER.OUT

C///  OPEN SUN FILES.

C*****SUN: OPEN STATEMENTS
      OPEN (FILE = 'plot2d_twafer.dat', FORM = 'FORMATTED',
     +   STATUS = 'UNKNOWN', UNIT = DATA1)
      OPEN (FILE = 'temp.dat', FORM = 'FORMATTED', STATUS = 'UNKNOWN',
     +   UNIT = DATA2)
      OPEN (FILE = 'twafer.in', FORM = 'FORMATTED', STATUS = 'OLD',
     +   UNIT = SCRIPT)
      OPEN (FILE = 'twafer.out', FORM = 'FORMATTED',
     +   STATUS = 'UNKNOWN', UNIT = TEXT)
      OPEN (FILE = 't.dat', FORM = 'FORMATTED', STATUS = 'UNKNOWN',
     +   UNIT = DATA4)
C*****END SUN: OPEN STATEMENTS

C///  PRINT.

      IF (0 .LT. TEXT) WRITE (TEXT, 10001)
     +   ID, 'VERSION 1.48 OF JANUARY 1998', BANNER

C///  CHECK THE ARGUMENTS.

      ERROR = .NOT. (1 .LE. CSIZE)
      IF (ERROR) GO TO 9001

      ERROR = .NOT. (1 .LE. ISIZE)
      IF (ERROR) GO TO 9002

      ERROR = .NOT. (1 .LE. LSIZE)
      IF (ERROR) GO TO 9003

      ERROR = .NOT. (1 .LE. RSIZE)
      IF (ERROR) GO TO 9004

C///  INITIALIZE THE STACK POINTERS.

      CLAST = 0
      CMAX = 0
      ILAST = 0
      IMAX = 0
      LLAST = 0
      LMAX = 0
      RLAST = 0
      RMAX = 0

      TIME2 = 0.0
      FLAG2 = .FALSE.

C     SET FLAG TO PRINT PLOT2D HEADER IN PLOT FILE.

      FLAG3 = .TRUE.
      KOUNT = 0

C///////////////////////////////////////////////////////////////////////
C
C     (2) SET DEFAULTS.
C
C     MANY OF THESE DEFAULTS MAY NOT BE ENFORCED IN THE SENSE THAT THE
C     INPUT SCRIPT MUST BE REQUIRED TO CONTAIN OVERIDING VALUES.
C
C///////////////////////////////////////////////////////////////////////

C///  INCLUDE PEDESTAL APPROXIMATION

      PED = .FALSE.
      PCON = .FALSE.

C     PEDESTAL MODEL
C     PED  = .TRUE. - TURN PEDESTAL MODEL ON.
C     PCON = .TRUE. - INCLUDE CONDUCTION IN THE CYLINDRICAL SLEEVE
C                     DEFINED BY PEDRI AND PEDRO, THE ADDED THERMAL
C                     MASS IS ALWAYS INCLUDED.
C     TTHIC - THICKNESS OF TOP QUARTZ DISC.
C     PEDRI - INNER RADIUS OF CYLINDER.
C     PEDR0 - OUTER RADIUS OF CYLINDER.

C*****USE PEDESTAL ONE APPROXIMATION
C      PED = .TRUE.
C      PCON = .TRUE.
C      TTHIC = 1.7078
C      PEDRI = 8.8
C      PEDRO = 9.2
C*****END USE PEDESTAL ONE APPROXIMATION

C*****USE PEDESTAL TWO APPROXIMATION
      PED = .FALSE.
      PCON = .FALSE.
      TTHIC = 3.04
      PEDRI = 10.91
      PEDRI = 11.8
C*****END USE PEDESTAL TWO APPROXIMATION

      IF (.NOT.PED) PCON = .FALSE.

C///  TURN AIR COOLING ON OR OF AND SET INITIAL TEMPERATURES.

      CFLOW(1) = .FALSE.
      CFLOW(2) = .FALSE.

C///  INLET COOLING GAS TEMPERATURES (K).

      TTGAS(1) = 300.0
      TTGAS(2) = 300.0

C///  AIR FLOW RATES IN (FT**3/MIN AT STP).

      VFLOW(1) = 5.0
CC      VFLOW(2) = 0.01
      VFLOW(2) = 0.10

C///  SPECIFY NUMBER OF RADIAL INJECTORS AND FLOWRATES (FT**3/MIN AT
C     STP) AND AXIAL LOCATIONS.

C     RADIAL INJECTION MODEL (BETWEEN OUTER TUBE AND HEATERS)
C     MODEL ASSUMES INJECTION UNIFORMLY BETWEE ZBEG(I) AND ZEND(I).
C     FLOW ONLY ALLOWED IN THE UPWARD DIRECTION AND A SMALL AMOUNT
C     OF AXIAL FLOW (SAY 0.5) MUST BE SET FOR THE CODE TO WORK
C     PROPERLY (CFLOW(2) MUST EQUAL .TRUE.).
C
C     RINJEC = .TRUE. - TURN RADIAL INJECTION MODEL ON.
C     ZBEG(I) - BEGINNING AXIAL LOCATION OF INJECTOR I.
C     ZEND(I) - END LOCATION OF INJECTOR I.
C     VINJ(I) - AIR FLOW INJECTION RATE FOR INJECTOR I (FT**3/MIN).
C     NUMIN - NUMBER OF INECTORS (MAXIMUM OF 10).

      RINJEC = .FALSE.

      NUMIN = 1

      ZBEGI(1) = 1.0
      ZENDI(1) = 95.0
CC      VINJ(1) = 300.0
      VINJ(1) = 0.0

C     CONVERT FROM (FT**3/MIN TO G/SEC).

      MFLOW(1) = VFLOW(1)*472.01/773.2
      MFLOW(2) = VFLOW(2)*472.01/773.2

C///  SPECIFY THE DISKS.

C     RADIUS OF THE DISKS
C     SPACING BETWEEN THE DISKS

      RADIUS = 4.0
      SPACE = 1.0

C///  SPECIFY THE GRID.

C     ADAPT THE GRID
C     MAXIMUM POINTS
C     MAXIMUM POINTS TO ADD AT ONCE
C     INITIAL POINTS

      ADAPT = .FALSE.

C///  SPECIFY THE OVEN.

      ZLEN = 50.0
      ZFIRST = 0.0
      ZLAST = ZLEN

C///  SPECIFY THE PRINT LEVEL.

C     INFORMATION DEPTH
C     SOLUTION DEPTH

      LEVEL(1) = 1
      LEVEL(2) = 1

C///  SPECIFY TWOPNT CONTROL DEFAULTS.

C     INITIAL STEPS
C     MAXIMUM SIZE
C     NEWTON ABSOLUTE ACCURACY
C     NEWTON RELATIVE ACCURACY
C     NEWTON STEPS BETWEEN JACOBIANS FOR TIME DEPENDENT PROBLEM
C     NEWTON STEPS BETWEEN JACOBIANS FOR STEADY STATE
C     SIZE
C     SIZE DECREASE FACTOR
C     SIZE INCREASE FACTOR
C     SMOOTHING STEPS
C     STEPS BEFORE SIZE INCREASES

      STEPS0 = 5
      TMAX = 100.0
      TMIN = 1.0E-10
      TDABS = 1.0E-9
      TDREL = 1.0E-6
      TDAGE = 20
      SSAGE = 20
      TSTEP0 = 1.0
      TDEC = 4.0
      TINC = 4.0
      STEPS1 = 50
      STEPS2 = 200

C     ABSOLUTE BOUND FOR STEADY STATE CONVERGENCE
C     RELATIVE BOUND FOR STEADY STATE CONVERGENCE
C     STEADY STATE RETIREMENT AGE FOR JACOBIAN

      SSABS = 1.0E-09
      SSREL = 1.0E-06
      SSAGE = 20

C     SOLUTION PRINTING DEPTH
C     INFORMATION PRINTING DEPTH

      LEVEL(1) = 0
      LEVEL(2) = 0

C     RELATIVE CHANGE IN SLOPE
C     RELATIVE CHANGE IN VALUE
C     THRESHOLD FOR ABSOLUTE AND RELATIVE CHANGE IN VALUE

      TOLER2 = 0.5
      TOLER1 = 0.5
      TOLER0 = 1.0E-9

C///  WRITE THE COEFFICIENT MATRIX.

      MATRIX = .FALSE.

C///  SPECIFY DASSL CONTROL DEFAULTS.

C     ABSOLUTE ERROR = ?! (real)
C     RELATIVE ERROR = ?! (real)
C     INTEGRATION STOP TIME = ?! (real)
C     TIME INTERVAL FOR OUTPUT = ?! (real)
C     MAXIMUM ORDER OF TIME INTEGRATION = ?! (integer.le.5)

      ATOL = 1.0E-08
      RTOL = 1.0E-06
      TSTOP = 1800.0
      DTIME = 60.0
      MAXORD = 5

C///////////////////////////////////////////////////////////////////////
C
C     (3) READ THE SCRIPT FILE.
C
C///////////////////////////////////////////////////////////////////////

C///  RESERVE SPACE FOR THE HEATER BEGINNING AND ENDING
C     Z LOCATIONS AND THE HEATER POWERS

C     ZBEG(MAXHT)
C     REAL INJECT

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QZBEG', MAXHT, QZBEG)
      IF (ERROR) GO TO 9101

C     ZEND(MAXHT)
C     REAL INJECT

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QZEND', MAXHT, QZEND)
      IF (ERROR) GO TO 9101

C     POWERH(MAXHT)
C     REAL INJECT

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QPOWRH', MAXHT, QPOWRH)
      IF (ERROR) GO TO 9101

C///  RESERVE SPACE FOR REACTOR WALL TEMPERATURE DATA
C     Z LOCATIONS AND THE TEMPERATURES

C     ZDAT(TEMMAX + 2)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QZDAT', TEMMAX + 2, QZDAT)
      IF (ERROR) GO TO 9101

C     TDAT(TEMMAX + 2)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTDAT', TEMMAX + 2, QTDAT)
      IF (ERROR) GO TO 9101

C///  RESERVE SPACE FOR REACTOR THERMOCOUPLE CALCULATIONS.

C     GBOT(NUMTC)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QGBOT', NUMTC, QGBOT)
      IF (ERROR) GO TO 9101

C     GTOP(NUMTC)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QGTOP', NUMTC, QGTOP)
      IF (ERROR) GO TO 9101

C     JBOT(NUMTC)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QJBOT', NUMTC, QJBOT)
      IF (ERROR) GO TO 9101

C     JTOP(NUMTC)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QJTOP', NUMTC, QJTOP)
      IF (ERROR) GO TO 9101

C     RING(4)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QRING', 5, QRING)
      IF (ERROR) GO TO 9101

C     SUMB(NUMTC)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QSUMB', NUMTC, QSUMB)
      IF (ERROR) GO TO 9101

C     SUMT(NUMTC)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QSUMT', NUMTC, QSUMT)
      IF (ERROR) GO TO 9101

C     TC(NUMTC)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTCT', NUMTC, QTCT)
      IF (ERROR) GO TO 9101

C     ZTC(NUMTC)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QZTC', NUMTC, QZTC)
      IF (ERROR) GO TO 9101

C///  PRINT.

      IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 10002) ID, 'READING THE SCRIPT FILE.'
         WRITE (TEXT, '()')
      END IF

C///  INITIALIZE THE LINE COUNTER AND THE BLOCK FOUND FLAGS.

      NUMBER = 0

      DO 3010 J = 1, BLOCKS
         BFOUND(J) = .FALSE.
3010  CONTINUE

C///  CHECK THE VERSION NUMBER.
C 1>  TWAFER VERSION ?.??

      CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT,
     +   'TWAFER VERSION')
      IF (ERROR) GO TO 9301

      CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)
      IF (ERROR) GO TO 9302

      DATA VMSS
     +   /'1.44', '1.45', '1.46', '1.47', '1.48' /

      FLAG = .FALSE.
      DO 3020 J = 1, VMSSS
         FLAG = FLAG .OR. WORD .EQ. VMSS(J)
3020  CONTINUE
      ERROR = .NOT. FLAG
      IF (ERROR) GO TO 9303

C///  ZERO HEATER POWERS

      DO 3025 J = 1, MAXHT
         R(QPOWRH + J - 1) = 0.0
3025  CONTINUE

      PCAP(1) = 0.0
      PCAP(2) = 0.0

C///  TOP OF THE LOOP THROUGH THE SCRIPT FILE.

3030  CONTINUE

      CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)
      IF (ERROR) GO TO 9302

C///  TOP OF THE SPECIFICATION BLOCKS.

      IF (WORD .EQ. 'SPECIFY') THEN

      CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, 'THE')
      IF (ERROR) GO TO 9301

      CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)
      IF (ERROR) GO TO 9302

C///  SPECIFY THE GRID.
C 1>  SPECIFY THE GRID
C 2>     NUMBER OF POINTS IN THE ENTRANCE REGION = ?!  (integer)
C 2>     NUMBER OF POINTS IN THE EXIT REGION = ?! (integer)
C 2>     NUMBER OF WAFERS = ?!  (integer)
C 2>     NUMBER OF RADIAL POINTS ACROSS WAFERS = ?! (integer)
C 2>     END

      DATA
     +   BNAME(QGRID), BNEED(QGRID),
     +   (TYPE(J, QGRID), NEED(J, QGRID), KEY(J, QGRID), J = 1, 4)
     + / 'SPECIFY THE GRID', .TRUE.,
     +   'I', .TRUE., 'NUMBER OF POINTS IN THE ENTRANCE REGION =',
     +   'I', .TRUE., 'NUMBER OF POINTS IN THE EXIT REGION =',
     +   'I', .TRUE., 'NUMBER OF WAFERS =',
     +   'I', .TRUE., 'NUMBER OF RADIAL POINTS ACROSS WAFERS ='/

      IF (WORD .EQ. 'GRID') THEN
         Q = QGRID
         ERROR = BFOUND(Q)
         IF (ERROR) GO TO 9304
         BFOUND(Q) = .TRUE.

C     SUBROUTINE READV
C    +  (ERROR, TEXT,
C    +   CKEY, ESIGN, FOUND, KEY, KEYS, IKEY, LINE, LKEY, NEED,
C    +   NUMBER, PRINT, RKEY, SCRIPT, TYPE)

C*****DOUBLE PRECISION
         CALL READV2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C         CALL READV1
C*****END SINGLE PRECISION
     +     (ERROR, TEXT,
     +      CKEY, .FALSE., FOUND, KEY(1, Q), 4, IKEY, LINE, LKEY,
     +      NEED(1, Q), NUMBER, PRINT, RKEY, SCRIPT, TYPE(1, Q))
         IF (ERROR) GO TO 9305

         PIN = IKEY(1)
         PEX = IKEY(2)
         NUMWF = IKEY(3)
         RPNTS = IKEY(4)

C///  SPECIFY THE REACTOR.
C 1>  SPECIFY THE REACTOR
C 2>     HEATER INNER RADIUS = ?! (cm)
C 2>     HEATER OUTER RADIUS = ?! (cm)
C 2>     TUBE LENGTH = ?! (cm)
C 2>     INSULATION THICKNESS ON HEATERS = ?! (cm)
C 2>     LOCATION OF FIRST WAFER = ?! (cm)
C 2>     LOCATION OF LAST WAFER = ?! (cm)
C 2>     WAFER RADIUS = ?! (cm)
C 2>     WAFER THICKNESS = ?! (cm)
C 2>     EMISS. OF WAFERS = ?! (unitless)
C 2>     EMISS. OF REACTOR HEATERS = ?! (unitless)
C 2>     EMISS. OF REACTOR ENTRANCE AND EXIT CAPS = ?! (unitless)
C 2>     SPECIFIC HEAT OF WAFER MATERIAL = ?! (J/g/K)
C 2>     DENSITY OF WAFER MATERIAL = ?! (g/cm**3)
C 2>     THERMAL CONDUCTIVITY OF WAFER MATERIAL = ?! (W/cm/K)
C 2>     SPECIFIC HEAT OF HEATER MATERIAL = ?! (J/g/K)
C 2>     DENSITY OF HEATER MATERIAL = ?! (g/cm**3)
C 2>     THERMAL CONDUCTIVITY OF HEATER MATERIAL = ?! (W/cm/K)
C 2>     SPECIFIC HEAT OF INSULATION MATERIAL = ?! (J/g/K)
C 2>     DENSITY OF INSULATION MATERIAL = ?! (g/cm**3)
C 2>     THERMAL CONDUCTIVITY OF INSULATION MATERIAL = ?! (W/cm/K)
C 2>     THICKNESS ENTRANCE CAP = ?! (cm)
C 2>     THICKNESS INSULATION ON ENTRANCE CAP = ?! (cm)
C 2>     THICKNESS EXIT CAP = ?! (cm)
C 2>     THICKNESS INSULATION ON EXIT CAP = ?! (cm)
C 2>     DENSITY OF QUARTZ TUBE MATERIAL = ?! (g/cm**3)
C 2>     THERMAL CONDUCTIVITY OF QUARTZ TUBE MATERIAL = ?! (W/cm/K)
C 2>     SPECIFIC HEAT OF QUARTZ TUBE MATERIAL = ?! (J/g/k)
C 2>     INNER RADIUS OF FIRST QUARTZ TUBE = ?! (cm)
C 2>     OUTER RADIUS OF FIRST QUARTZ TUBE = ?! (cm)
C 2>     TRANSMITTANCE OF FIRST QUARTZ TUBE = ?! (unitless, T+R+A=1)
C 2>     REFLECTANCE OF FIRST QUARTZ TUBE = ?! (unitless, T+R+A=1)
C 2>     ABSORPTANCE OF FIRST QUARTZ TUBE = ?! (unitless, T+R+A=1)
C 2>     INNER RADIUS OF SECOND QUARTZ TUBE = ?! (cm)
C 2>     OUTER RADIUS OF SECOND QUARTZ TUBE = ?! (cm)
C 2>     TRANSMITTANCE OF SECOND QUARTZ TUBE = ?! (unitless, T+R+A=1)
C 2>     REFLECTANCE OF SECOND QUARTZ TUBE = ?! (unitless, T+R+A=1)
C 2>     ABSORPTANCE OF SECOND QUARTZ TUBE = ?! (unitless, T+R+A=1)
C 2>     EMISS. OF MAIN SECTION OUTER SKIN = ?! (unitless)
C 2>     EMISS. OF ENTRANCE SECTION OUTER SKIN = ?! (unitless)
C 2>     EMISS. OF EXIT SECTION OUTER SKIN = ?! (unitless)
C 2>     MAIN SECTION CONVECTIVE HEAT LOSS COEF. = ?! (W/cm**2/K)
C 2>     ENTRANCE CAP CONVECTIVE HEAT LOSS COEF. = ?! (W/cm**2/K)
C 2>     EXIT CAP CONVECTIVE HEAT LOSS COEF. = ?! (W/cm**2/K)
C 2>     ENTRANCE SECTION CONVECTIVE HEAT LOSS COEF. = ?! (W/cm**2/K)
C 2>     EXIT SECTION CONVECTIVE HEAT LOSS COEF. = ?! (W/cm**2/K)
C 2>     DENSITY ENTRANCE CAP MATERIAL = ?! (g/cm**3)
C 2>     DENSITY OF EXIT CAP MATERIAL = ?! (g/cm**3)
C 2>     DENSITY OF INSULATION ON ENTRANCE CAP = ?! (g/cm**3)
C 2>     DENSITY OF INSULATION ON EXIT CAP = ?! (g/cm**3)
C 2>     THERMAL CONDUCTIVITY OF ENTRANCE CAP = ?! (W/cm/K)
C 2>     THERMAL CONDUCTIVITY OF EXIT CAP = ?! (W/cm/K)
C 2>     THERMAL CONDUCTIVITY OF ENTRANCE CAP INSULATION = ?! (W/cm/K)
C 2>     THERMAL CONDUCTIVITY OF EXIT CAP INSULATION = ?! (W/cm/K)
C 2>     SPECIFIC HEAT OF ENTRANCE CAP = ?! (J/g/K)
C 2>     SPECIFIC HEAT OF EXIT CAP = ?! (J/g/K)
C 2>     SPECIFIC HEAT OF ENTRANCE CAP INSULATION = ?! (J/g/K)
C 2>     SPECIFIC HEAT OF EXIT CAP INSULATION = ?! (J/g/K)
C 2>     EMISS. OF ENTRANCE CAP OUTER SKIN = ?! (unitless)
C 2>     EMISS. OF EXIT CAP OUTER SKIN = ?! (unitless)
C 2>     END

      DATA
     +   BNAME(QREAC), BNEED(QREAC),
     +   (TYPE(J, QREAC), NEED(J, QREAC), KEY(J, QREAC), J = 1, 10)
     + / 'SPECIFY THE OVEN', .TRUE.,
     +   'R', .TRUE.,  'HEATER INNER RADIUS =',
     +   'R', .TRUE.,  'HEATER OUTER RADIUS =',
     +   'R', .TRUE.,  'TUBE LENGTH =',
     +   'R', .TRUE.,  'INSULATION THICKNESS ON HEATERS =',
     +   'R', .TRUE.,  'LOCATION OF FIRST WAFER =',
     +   'R', .TRUE.,  'LOCATION OF LAST WAFER =',
     +   'R', .TRUE.,  'WAFER RADIUS =',
     +   'R', .TRUE.,  'WAFER THICKNESS =',
     +   'R', .TRUE.,  'EMISS. OF WAFERS =',
     +   'R', .TRUE.,  'EMISS. OF REACTOR HEATERS =' /

      DATA
     +   (TYPE(J, QREAC), NEED(J, QREAC), KEY(J, QREAC), J = 11, 20)
     + / 'R', .TRUE.,  'EMISS. OF REACTOR ENTRANCE AND EXIT CAPS =',
     +   'R', .TRUE.,  'SPECIFIC HEAT OF WAFER MATERIAL =',
     +   'R', .TRUE.,  'DENSITY OF WAFER MATERIAL =',
     +   'R', .TRUE.,  'THERMAL CONDUCTIVITY OF WAFER MATERIAL =',
     +   'R', .TRUE.,  'SPECIFIC HEAT OF HEATER MATERIAL =',
     +   'R', .TRUE.,  'DENSITY OF HEATER MATERIAL =',
     +   'R', .TRUE.,  'THERMAL CONDUCTIVITY OF HEATER MATERIAL =',
     +   'R', .TRUE.,  'SPECIFIC HEAT OF INSULATION MATERIAL =',
     +   'R', .TRUE.,  'DENSITY OF INSULATION MATERIAL =',
     +   'R', .TRUE.,  'THERMAL CONDUCTIVITY OF INSULATION MATERIAL =' /

      DATA
     +   (TYPE(J, QREAC), NEED(J, QREAC), KEY(J, QREAC), J = 21, 30)
     + / 'R', .TRUE.,  'THICKNESS ENTRANCE CAP =',
     +   'R', .TRUE.,  'THICKNESS INSULATION ON ENTRANCE CAP =',
     +   'R', .TRUE.,  'THICKNESS EXIT CAP =',
     +   'R', .TRUE.,  'THICKNESS INSULATION ON EXIT CAP =',
     +   'R', .TRUE.,  'EMISS. OF MAIN SECTION OUTER SKIN =',
     +   'R', .TRUE.,  'DENSITY OF QUARTZ TUBE MATERIAL =',
     +   'R', .TRUE.,  'THERMAL CONDUCTIVITY OF QUARTZ TUBE MATERIAL =',
     +   'R', .TRUE.,  'SPECIFIC HEAT OF QUARTZ TUBE MATERIAL =',
     +   'R', .TRUE.,  'INNER RADIUS OF FIRST QUARTZ TUBE =',
     +   'R', .TRUE.,  'OUTER RADIUS OF FIRST QUARTZ TUBE =' /

      DATA
     +   (TYPE(J, QREAC), NEED(J, QREAC), KEY(J, QREAC), J = 31, 40)
     + / 'R', .TRUE.,  'TRANSMITTANCE OF FIRST QUARTZ TUBE =',
     +   'R', .TRUE.,  'REFLECTANCE OF FIRST QUARTZ TUBE =',
     +   'R', .TRUE.,  'ABSORPTANCE OF FIRST QUARTZ TUBE =',
     +   'R', .TRUE.,  'INNER RADIUS OF SECOND QUARTZ TUBE =',
     +   'R', .TRUE.,  'OUTER RADIUS OF SECOND QUARTZ TUBE =',
     +   'R', .TRUE.,  'TRANSMITTANCE OF SECOND QUARTZ TUBE =',
     +   'R', .TRUE.,  'REFLECTANCE OF SECOND QUARTZ TUBE =',
     +   'R', .TRUE.,  'ABSORPTANCE OF SECOND QUARTZ TUBE =',
     +   'R', .TRUE.,  'MAIN SECTION CONVECTIVE HEAT LOSS COEF. =',
     +   'R', .TRUE.,  'ENTRANCE CAP CONVECTIVE HEAT LOSS COEF. =' /

      DATA
     +   (TYPE(J, QREAC), NEED(J, QREAC), KEY(J, QREAC), J = 41, 50)
     + / 'R', .TRUE.,  'EXIT CAP CONVECTIVE HEAT LOSS COEF. =',
     +   'R', .TRUE.,  'DENSITY ENTRANCE CAP MATERIAL =',
     +   'R', .TRUE.,  'DENSITY OF EXIT CAP MATERIAL =',
     +   'R', .TRUE.,  'DENSITY OF INSULATION ON ENTRANCE CAP =',
     +   'R', .TRUE.,  'DENSITY OF INSULATION ON EXIT CAP =',
     +   'R', .TRUE.,  'THERMAL CONDUCTIVITY OF ENTRANCE CAP =',
     +   'R', .TRUE.,  'THERMAL CONDUCTIVITY OF EXIT CAP =',
     +   'R', .TRUE.,
     +   'THERMAL CONDUCTIVITY OF ENTRANCE CAP INSULATION =',
     +   'R', .TRUE.,  'THERMAL CONDUCTIVITY OF EXIT CAP INSULATION =',
     +   'R', .TRUE.,  'SPECIFIC HEAT OF ENTRANCE CAP =' /

      DATA
     +   (TYPE(J, QREAC), NEED(J, QREAC), KEY(J, QREAC), J = 51, 59)
     + / 'R', .TRUE.,  'SPECIFIC HEAT OF EXIT CAP =',
     +   'R', .TRUE.,  'SPECIFIC HEAT OF ENTRANCE CAP INSULATION =',
     +   'R', .TRUE.,  'SPECIFIC HEAT OF EXIT CAP INSULATION =',
     +   'R', .TRUE.,  'EMISS. OF ENTRANCE CAP OUTER SKIN =',
     +   'R', .TRUE.,  'EMISS. OF EXIT CAP OUTER SKIN =',
     +   'R', .TRUE.,  'EMISS. OF ENTRANCE SECTION OUTER SKIN =',
     +   'R', .TRUE.,  'EMISS. OF EXIT SECTION OUTER SKIN =',
     +   'R', .TRUE.,  'ENTRANCE SECTION CONVECTIVE HEAT LOSS COEF. =',
     +   'R', .TRUE.,  'EXIT SECTION CONVECTIVE HEAT LOSS COEF. ='/

      ELSE IF (WORD .EQ. 'REACTOR') THEN
         Q = QREAC
         ERROR = BFOUND(Q)
         IF (ERROR) GO TO 9304
         BFOUND(Q) = .TRUE.

C     SUBROUTINE READV
C    +  (ERROR, TEXT,
C    +   CKEY, ESIGN, FOUND, KEY, KEYS, IKEY, LINE, LKEY, NEED,
C    +   NUMBER, PRINT, RKEY, SCRIPT, TYPE)

C*****DOUBLE PRECISION
         CALL READV2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C         CALL READV1
C*****END SINGLE PRECISION
     +     (ERROR, TEXT,
     +      CKEY, .FALSE., FOUND, KEY(1, Q), 59, IKEY, LINE, LKEY,
     +      NEED(1, Q), NUMBER, PRINT, RKEY, SCRIPT, TYPE(1, Q))
         IF (ERROR) GO TO 9305

         RTUBEI = RKEY(1)
         RTUBEO = RKEY(2)
         ZLEN   = RKEY(3)
         THICKI = RKEY(4)
         ZFIRST = RKEY(5)
         ZLAST  = RKEY(6)
         RADIUS = RKEY(7)
         THICK  = RKEY(8)
         EMISS  = RKEY(9)
         EMISSW = RKEY(10)
         EMISSE = RKEY(11)
         CP(1)     = RKEY(12)
         RHO(1)    = RKEY(13)
         CONDUC(1) = RKEY(14)
         CP(2)     = RKEY(15)
         RHO(2)    = RKEY(16)
         CONDUC(2) = RKEY(17)
         CP(3)     = RKEY(18)
         RHO(3)    = RKEY(19)
         CONDUC(3) = RKEY(20)
         DZEND(1)  = RKEY(21)
         DZENDI(1) = RKEY(22)
         DZEND(2)  = RKEY(23)
         DZENDI(2) = RKEY(24)
         EMISSI(3) = RKEY(25)
         RHO(4) = RKEY(26)
         CONDUC(4) = RKEY(27)
         CP(4) = RKEY(28)
         RINQ(1)  = RKEY(29)
         ROUTQ(1) = RKEY(30)
         TRAN(1) = RKEY(31)
         REFL(1) = RKEY(32)
         EMIS(1) = RKEY(33)
         RINQ(2)  = RKEY(34)
         ROUTQ(2) = RKEY(35)
         TRAN(2) = RKEY(36)
         REFL(2) = RKEY(37)
         EMIS(2) = RKEY(38)
         HCOEF(3) = RKEY(39)
         HCOEF(1) = RKEY(40)
         HCOEF(2) = RKEY(41)
         RHO(5) = RKEY(42)
         RHO(6) = RKEY(43)
         RHO(7) = RKEY(44)
         RHO(8) = RKEY(45)
         CONDUC(5) = RKEY(46)
         CONDUC(6) = RKEY(47)
         CONDUC(7) = RKEY(48)
         CONDUC(8) = RKEY(49)
         CP(5) = RKEY(50)
         CP(6) = RKEY(51)
         CP(7) = RKEY(52)
         CP(8) = RKEY(53)
         EMISSI(1) = RKEY(54)
         EMISSI(2) = RKEY(55)
         EMISSI(4) = RKEY(56)
         EMISSI(5) = RKEY(57)
         HCOEF(4) = RKEY(58)
         HCOEF(5) = RKEY(59)

         ERROR = (RTUBEO .LE. RTUBEI)
         IF (ERROR) GO TO 9317
         ERROR = (ROUTQ(1) .LE. RINQ(1))
         IF (ERROR) GO TO 9318
         ERROR = (ROUTQ(2) .LE. RINQ(2))
         IF (ERROR) GO TO 9319

         ERROR = .NOT.((ROUTQ(2).LT.RTUBEO).AND.(ROUTQ(2).LT.RTUBEI))
         IF (ERROR) GO TO 9320
         ERROR = .NOT.((ROUTQ(1).LT.ROUTQ(2)).AND.(ROUTQ(1).LT.RINQ(2)))
         IF (ERROR) GO TO 9321
         ERROR = .NOT.(RADIUS .LE. RINQ(1))
         IF (ERROR) GO TO 9322

         ERROR = .NOT.(ABS(EMIS(1)+REFL(1)+TRAN(1)-1.0).LT.1.0E-08)
         IF (ERROR) GO TO 9323
         ERROR = .NOT.(ABS(EMIS(2)+REFL(2)+TRAN(2)-1.0).LT.1.0E-08)
         IF (ERROR) GO TO 9324

C///  SPECIFY THE HEATER TEMPERATURES.
C 1>  SPECIFY THE HEATER TEMPERATURES.
C 2>
C 3>     Z(1) = ?! (real) !CM
C 4>     T(1) = ?! (real) !K
C 5>     ..... etc.   .....up to a maximum of 10 points
C 6>     END

      DATA
     +   BNAME(QTEMP), BNEED(QTEMP),
     +   (TYPE(J, QTEMP), NEED(J, QTEMP), KEY(J, QTEMP), J=1,
     +   (TEMMAX*2))
     + / 'SPECIFY THE HEATER TEMPERATURES', .FALSE.,
     +   'R', .FALSE., 'Z(1) =',  'R', .FALSE., 'T(1) =',
     +   'R', .FALSE., 'Z(2) =',  'R', .FALSE., 'T(2) =',
     +   'R', .FALSE., 'Z(3) =',  'R', .FALSE., 'T(3) =',
     +   'R', .FALSE., 'Z(4) =',  'R', .FALSE., 'T(4) =',
     +   'R', .FALSE., 'Z(5) =',  'R', .FALSE., 'T(5) =',
     +   'R', .FALSE., 'Z(6) =',  'R', .FALSE., 'T(6) =',
     +   'R', .FALSE., 'Z(7) =',  'R', .FALSE., 'T(7) =',
     +   'R', .FALSE., 'Z(8) =',  'R', .FALSE., 'T(8) =',
     +   'R', .FALSE., 'Z(9) =',  'R', .FALSE., 'T(9) =',
     +   'R', .FALSE., 'Z(10) =', 'R', .FALSE., 'T(10) ='/

C///  SPECIFY THE HEATER POWERS
C 1>  SPECIFY THE HEATER POWERS
C 2>     ZBEG(1)  = ?! (real) !cm
C 3>     ZEND(1)  = ?! (real) !cm
C 4>     POWER(1) = ?! (real) !WATTS
C 5>     ..... etc.   .....up to a maximum of 10 heater elements
C 6>     END

      DATA
     +   BNAME(QHT), BNEED(QHT)
     + / 'SPECIFY THE HEATER POWERS', .FALSE. /

      DATA
     +   (TYPE(J, QHT), NEED(J, QHT), KEY(J, QHT), J=1, 3*5)
     + / 'R', .FALSE., 'ZBEG(1) =', 'R', .FALSE., 'ZEND(1) =',
     +   'R', .FALSE., 'POWER(1) =',
     +   'R', .FALSE., 'ZBEG(2) =', 'R', .FALSE., 'ZEND(2) =',
     +   'R', .FALSE., 'POWER(2) =',
     +   'R', .FALSE., 'ZBEG(3) =', 'R', .FALSE., 'ZEND(3) =',
     +   'R', .FALSE., 'POWER(3) =',
     +   'R', .FALSE., 'ZBEG(4) =', 'R', .FALSE., 'ZEND(4) =',
     +   'R', .FALSE., 'POWER(4) =',
     +   'R', .FALSE., 'ZBEG(5) =', 'R', .FALSE., 'ZEND(5) =',
     +   'R', .FALSE., 'POWER(5) =' /

      DATA
     +   (TYPE(J, QHT), NEED(J, QHT), KEY(J, QHT), J=3*5 + 1, 3*MAXHT)
     + / 'R', .FALSE., 'ZBEG(6) =', 'R', .FALSE., 'ZEND(6) =',
     +   'R', .FALSE., 'POWER(6) =',
     +   'R', .FALSE., 'ZBEG(7) =', 'R', .FALSE., 'ZEND(7) =',
     +   'R', .FALSE., 'POWER(7) =',
     +   'R', .FALSE., 'ZBEG(8) =', 'R', .FALSE., 'ZEND(8) =',
     +   'R', .FALSE., 'POWER(8) =',
     +   'R', .FALSE., 'ZBEG(9) =', 'R', .FALSE., 'ZEND(9) =',
     +   'R', .FALSE., 'POWER(9) =',
     +   'R', .FALSE., 'ZBEG(10) =', 'R', .FALSE., 'ZEND(10) =',
     +   'R', .FALSE., 'POWER(10) =' /

      ELSE IF (WORD .EQ. 'HEATER') THEN
         CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)
         IF (ERROR) GO TO 9302

         IF (WORD .EQ. 'POWERS') THEN
            Q = QHT
            ERROR = BFOUND(Q)
            IF (ERROR) GO TO 9304
            BFOUND(Q) = .TRUE.
            ERROR = BFOUND(QTEMP)
            IF (ERROR) GO TO 9306

C     SUBROUTINE READV
C    +  (ERROR, TEXT,
C    +   CVALUE, ESIGN, FOUND, KEY, KEYS, IVALUE, LINE, LVALUE,
C    +   NEED, NUMBER, PRINT, RVALUE, SCRIPT, TYPE)

C*****DOUBLE PRECISION
            CALL READV2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C           CALL READV1
C*****END SINGLE PRECISION
     +        (ERROR, TEXT,
     +         CKEY, .FALSE., FOUND, KEY(1,Q), 3*MAXHT,
     +         IKEY, LINE, LKEY, NEED(1,Q), NUMBER,
     +         PRINT, RKEY, SCRIPT, TYPE(1,Q))
            IF (ERROR) GO TO 9305
            NUMHT = 0
            DO 3081 J = 1, 3*MAXHT
               IF(FOUND(J)) NUMHT = NUMHT + 1
3081        CONTINUE
            ERROR = (MOD(NUMHT, 3) .GT. 0)
            IF (ERROR) GO TO 9314
            NUMHT = NUMHT/3

C///  CHECK THAT ALL VALUES ARE DEFINED FOR EACH HEATER

            DO 3085 J = 1, NUMHT
               ERROR = .NOT. (FOUND(J).AND.FOUND(J+1).AND.FOUND(J+2))
               R(QZBEG + J - 1) = RKEY((J-1)*3 + 1)
               R(QZEND + J - 1) = RKEY((J-1)*3 + 2)
               R(QPOWRH + J - 1) = RKEY((J-1)*3 + 3)
3085        CONTINUE
            IF (ERROR) GO TO 9314

         ELSE IF (WORD .EQ. 'TEMPERATURES') THEN
            Q = QTEMP
            ERROR = BFOUND(Q)
            IF (ERROR) GO TO 9304
            BFOUND(Q) = .TRUE.
            ERROR = BFOUND(QHT)
            IF (ERROR) GO TO 9306

C///  CALL READV.

C     SUBROUTINE READV
C    +  (ERROR, TEXT,
C    +   CKEY, ESIGN, FOUND, KEY, KEYS, IKEY, LINE, LKEY,
C    +   NEED, NUMBER, PRINT, RKEY, SCRIPT, TYPE)

C*****DOUBLE PRECISION
            CALL READV2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C           CALL READV1
C*****END SINGLE PRECISION
     +        (ERROR, TEXT,
     +         CKEY, .FALSE., FOUND, KEY(1,Q), 2*TEMMAX,
     +         IKEY, LINE, LKEY, NEED(1,Q), NUMBER,
     +         PRINT, RKEY, SCRIPT, TYPE(1,Q))
            IF (ERROR) GO TO 9305

            NUMT = 0
            DO 3086 J = 1, 2*TEMMAX
               IF(FOUND(J+1)) NUMT = NUMT + 1
3086        CONTINUE

            ERROR = (MOD(NUMT, 2) .GT. 0)
            IF (ERROR) GO TO 9316
            NUMT = NUMT/2

            COUNT = 1
            DO 3087 J = 1, NUMT
               ERROR = .NOT. (FOUND(J) .AND. FOUND(J + 1))
               R(QZDAT + J - 1) = RKEY(COUNT)
               COUNT = COUNT + 1
               R(QTDAT + J - 1) = RKEY(COUNT)
               COUNT = COUNT + 1
 3087       CONTINUE
            IF(ERROR) GO TO 9316
         ELSE
            ERROR = .TRUE.
            GO TO 9602
         ENDIF

C-----------------------------------------------------------------------
C///  SPECIFY THE RING BOAT.
C 1>  SPECIFY THE RING BOAT.
C 2>     RING INSIDE RADIUS = ?! (real)
C 2>     RING THICKNESS = ?! (real)
C 2>     RING MATERIAL DENSITY = ?! (real)
C 2>     RING MATERIAL SPECIFIC HEAT = ?! (real)
C 2>     RING MATERIAL THERMAL CONDUCTIVITY ?! (real)
C 2>     END

      DATA
     +   BNAME(QRBOAT), BNEED(QRBOAT),
     +   (TYPE(J, QRBOAT), NEED(J, QRBOAT), KEY(J, QRBOAT), J = 1, 5)
     + / 'SPECIFY THE RING BOAT', .FALSE.,
     +   'R', .FALSE., 'RING INSIDE RADIUS =',
     +   'R', .FALSE., 'RING THICKNESS =',
     +   'R', .FALSE., 'RING MATERIAL DENSITY =',
     +   'R', .FALSE., 'RING MATERIAL SPECIFIC HEAT =',
     +   'R', .FALSE., 'RING MATERIAL THERMAL CONDUCTIVITY ='/
      ELSE IF (WORD .EQ. 'RING') THEN
         CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT,
     +      'BOAT')
         IF (ERROR) GO TO 9301

         Q = QRBOAT
         ERROR = BFOUND(Q)
         IF (ERROR) GO TO 9304
         BFOUND(Q) = .TRUE.

C     SUBROUTINE READV
C    +  (ERROR, TEXT,
C    +   CKEY, ESIGN, FOUND, KEY, KEYS, IKEY, LINE, LKEY, NEED,
C    +   NUMBER, PRINT, RKEY, SCRIPT, TYPE)

C*****DOUBLE PRECISION
         CALL READV2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C         CALL READV1
C*****END SINGLE PRECISION
     +     (ERROR, TEXT,
     +      CKEY, .FALSE., FOUND, KEY(1, Q), 5, IKEY, LINE, LKEY,
     +      NEED(1, Q), NUMBER, PRINT, RKEY, SCRIPT, TYPE(1, Q))
         IF (ERROR) GO TO 9305

         RBOAT = BFOUND(Q)
         IF (FOUND(1)) R(QRING) = RKEY(1)
         IF (FOUND(2)) R(QRING+1) = RKEY(2)
         IF (FOUND(3)) R(QRING+2) = RKEY(3)
         IF (FOUND(4)) R(QRING+3) = RKEY(4)
         IF (FOUND(4)) R(QRING+4) = RKEY(5)
C-----------------------------------------------------------------------

C///  SPECIFY THE TWOPNT CONTROLS.
C 1>  SPECIFY THE TWOPNT CONTROLS
C 2>     SOLULTION PRINTING DEPTH = ?! (integer)
C 2>     INFORMATION PRINTING DEPTH = ?! (integer)
C 2>     STEADY STATE ABSOLUTE CONVERGENCE TEST = ?! (real)
C 2>     STEADY STATE RELATIVE CONVERGENCE TEST = ?! (real)
C 2>     STEADY STATE JACOBIAN RETIREMENT AGE =?! (real)
C 2>     INITIAL NUMBER OF TIME STEPS = ?! (integer)
C 2>     MAXIMUM TIME STEP SIZE = ?!  (real)
C 2>     MINIMUM TIME STEP SIZE = ?!  (real)
C 2>     INITIAL TIME STEP SIZE = ?!  (real)
C 2>     DIVISOR FOR DECREASE OF TIME STEP SIZE = ?!  (real)
C 2>     MULTIPLIER FOR INCREASE OF TIME STEP SIZE = ?!  (real)
C 2>     TIME STEPS BEFORE STEP SIZE INCREASES = ?! (integer)
C 2>     TIME STEPS AFTER CONVERGENCE FAILURE = ?!  (integer)
C 2>     TIME DEPENDENT ABSOLUTE CONVERGENCE TEST = ?!  (real)
C 2>     TIME DEPENDENT RELATIVE CONVERGENCE TEST = ?!  (real)
C 2>     TIME DEPENDENT JACOBIAN RETIREMENT AGE = ?! (integer)
C 2>     END

      DATA
     +   BNAME(QTWO), BNEED(QTWO),
     +   (TYPE(J, QTWO), NEED(J, QTWO), KEY(J, QTWO), J = 1, 16)
     + / 'SPECIFY THE TWOPNT CONTROLS', .FALSE.,
     +   'I', .FALSE., 'SOLUTION PRINTING DEPTH =',
     +   'I', .FALSE., 'INFORMATION PRINTING DEPTH =',
     +   'R', .FALSE., 'STEADY STATE ABSOLUTE CONVERGENCE TEST =',
     +   'R', .FALSE., 'STEADY STATE RELATIVE CONVERGENCE TEST =',
     +   'I', .FALSE., 'STEADY STATE JACOBIAN RETIREMENT AGE =',
     +   'I', .FALSE., 'INITIAL NUMBER OF TIME STEPS =',
     +   'R', .FALSE., 'MAXIMUM TIME STEP SIZE =',
     +   'R', .FALSE., 'MINIMUM TIME STEP SIZE =',
     +   'R', .FALSE., 'INITIAL TIME STEP SIZE =',
     +   'R', .FALSE., 'DIVISOR FOR DECREASE OF TIME STEP SIZE =',
     +   'R', .FALSE., 'MULTIPLIER FOR INCREASE OF TIME STEP SIZE =',
     +   'I', .FALSE., 'TIME STEPS BEFORE STEP SIZE INCREASES =',
     +   'I', .FALSE., 'TIME STEPS AFTER CONVERGENCE FAILURE =',
     +   'R', .FALSE., 'TIME DEPENDENT ABSOLUTE CONVERGENCE TEST =',
     +   'R', .FALSE., 'TIME DEPENDENT RELATIVE CONVERGENCE TEST =',
     +   'I', .FALSE., 'TIME DEPENDENT JACOBIAN RETIREMENT AGE =' /
      ELSE IF (WORD .EQ. 'TWOPNT') THEN
         CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT,
     +      'CONTROLS')
         IF (ERROR) GO TO 9301

         Q = QTWO
         ERROR = BFOUND(Q)
         IF (ERROR) GO TO 9304
         BFOUND(Q) = .TRUE.

C     SUBROUTINE READV
C    +  (ERROR, TEXT,
C    +   CKEY, ESIGN, FOUND, KEY, KEYS, IKEY, LINE, LKEY, NEED,
C    +   NUMBER, PRINT, RKEY, SCRIPT, TYPE)

C*****DOUBLE PRECISION
         CALL READV2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C         CALL READV1
C*****END SINGLE PRECISION
     +     (ERROR, TEXT,
     +      CKEY, .FALSE., FOUND, KEY(1, Q), 16, IKEY, LINE, LKEY,
     +      NEED(1, Q), NUMBER, PRINT, RKEY, SCRIPT, TYPE(1, Q))
         IF (ERROR) GO TO 9305

         IF (FOUND(1)) LEVEL(1) = IKEY(1)
         IF (FOUND(2)) LEVEL(2) = IKEY(2)
         IF (FOUND(3)) SSABS = RKEY(3)
         IF (FOUND(4)) SSREL = RKEY(4)
         IF (FOUND(5)) SSAGE = IKEY(5)
         IF (FOUND(6)) STEPS0 = IKEY(6)
         IF (FOUND(7)) TMAX = RKEY(7)
         IF (FOUND(8)) TMIN = RKEY(8)
         IF (FOUND(9)) TSTEP0 = RKEY(9)
         IF (FOUND(10)) TDEC = RKEY(10)
         IF (FOUND(11)) TINC = RKEY(11)
         IF (FOUND(12)) STEPS2 = IKEY(12)
         IF (FOUND(13)) STEPS1 = IKEY(13)
         IF (FOUND(14)) TDABS = RKEY(14)
         IF (FOUND(15)) TDREL = RKEY(15)
         IF (FOUND(16)) TDAGE = IKEY(16)

C///  SPECIFY THE DASSL CONTROLS.
C 1>  SPECIFY THE DASSL CONTROLS.
C 2>     ABSOLUTE ERROR = ?! (real)
C 2>     RELATIVE ERROR = ?! (real)
C 2>     INTEGRATION STOP TIME = ?! (real)
C 2>     TIME INTERVAL FOR OUTPUT = ?! (real)
C 2>     MAXIMUM ORDER OF TIME INTEGRATION = ?! (integer.le.5)
C 2>     END

      DATA
     +   BNAME(QDAS), BNEED(QDAS),
     +   (TYPE(J, QDAS), NEED(J, QDAS), KEY(J, QDAS), J = 1, 5)
     + / 'SPECIFY THE TWOPNT CONTROLS', .FALSE.,
     +   'R', .FALSE., 'ABSOLUTE ERROR =',
     +   'R', .FALSE., 'RELATIVE ERROR =',
     +   'R', .FALSE., 'INTEGRATION STOP TIME =',
     +   'R', .FALSE., 'TIME INTERVAL FOR OUTPUT =',
     +   'I', .FALSE., 'MAXIMUM ORDER OF TIME INTEGRATION =' /
      ELSE IF (WORD .EQ. 'DASSL') THEN
         CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT,
     +      'CONTROLS')
         IF (ERROR) GO TO 9301

         Q = QDAS
         ERROR = BFOUND(Q)
         IF (ERROR) GO TO 9304
         BFOUND(Q) = .TRUE.

C     SUBROUTINE READV
C    +  (ERROR, TEXT,
C    +   CKEY, ESIGN, FOUND, KEY, KEYS, IKEY, LINE, LKEY, NEED,
C    +   NUMBER, PRINT, RKEY, SCRIPT, TYPE)

C*****DOUBLE PRECISION
         CALL READV2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C         CALL READV1
C*****END SINGLE PRECISION
     +     (ERROR, TEXT,
     +      CKEY, .FALSE., FOUND, KEY(1, Q), 5, IKEY, LINE, LKEY,
     +      NEED(1, Q), NUMBER, PRINT, RKEY, SCRIPT, TYPE(1, Q))
         IF (ERROR) GO TO 9305

         IF (FOUND(1)) ATOL = RKEY(1)
         IF (FOUND(2)) RTOL = RKEY(2)
         IF (FOUND(3)) TSTOP = RKEY(3)
         IF (FOUND(4)) DTIME = RKEY(4)
         IF (FOUND(5)) MAXORD = IKEY(1)

C///  SPECIFY THE INITIAL CONDITIONS.
C 1>  SPECIFY THE INITIAL CONDITIONS
C 2>       TEMPERATURE OF WAFERS = ?!  (real)
C 2>       TEMPERATURE OF HEATERS = ?!  (real)
C 2>       TEMPERATURE OF END CAPS = ?!  (real)
C 2>       TEMPERATURE OF INSULATION = ?!  (real)
C 2>       AMBIENT TEMPERATURE = ?!  (real)
C 2>       TEMPERATURE OF QUARTZ = ?! (real)
C 2>  END

      DATA
     +   BNAME(QINIT), BNEED(QINIT),
     +   (TYPE(J, QINIT), NEED(J, QINIT), KEY(J, QINIT), J = 1, 6)
     + / 'SPECIFY THE INITIAL CONDITIONS', .TRUE.,
     +   'R', .TRUE., 'TEMPERATURE OF WAFERS =',
     +   'R', .TRUE., 'TEMPERATURE OF HEATERS =',
     +   'R', .TRUE., 'TEMPERATURE OF END CAPS =',
     +   'R', .TRUE., 'TEMPERATURE OF INSULATION =',
     +   'R', .TRUE., 'AMBIENT TEMPERATURE =',
     +   'R', .TRUE., 'TEMPERATURE OF QUARTZ ='/

      ELSE IF (WORD .EQ. 'INITIAL') THEN
         CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT,
     +      'CONDITIONS')
         IF (ERROR) GO TO 9301

         Q = QINIT
         ERROR = BFOUND(Q)
         IF (ERROR) GO TO 9304
         BFOUND(Q) = .TRUE.

C     SUBROUTINE READV
C    +  (ERROR, TEXT,
C    +   CKEY, ESIGN, FOUND, KEY, KEYS, IKEY, LINE, LKEY, NEED,
C    +   NUMBER, PRINT, RKEY, SCRIPT, TYPE)

C*****DOUBLE PRECISION
         CALL READV2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C         CALL READV1
C*****END SINGLE PRECISION
     +     (ERROR, TEXT,
     +      CKEY, .FALSE., FOUND, KEY(1, Q), 6, IKEY, LINE, LKEY,
     +      NEED(1, Q), NUMBER, PRINT, RKEY, SCRIPT, TYPE(1, Q))
         IF (ERROR) GO TO 9305

         IF (FOUND(1)) TTWAF = RKEY(1)
         IF (FOUND(2)) TTUBE = RKEY(2)
         IF (FOUND(3)) TTEND = RKEY(3)
         IF (FOUND(4)) TTIN = RKEY(4)
         IF (FOUND(5)) TAMB = RKEY(5)
         IF (FOUND(6)) TTQ = RKEY(6)

C///  SPECIFY THE CAP TEMPERATURES.
C 1>  SPECIFY THE CAP TEMPERATURES.
C 2>       TEMPERATURE OF ENTRANCE CAP = ?!  (real)
C 2>       TEMPERATURE OF EXIT CAP = ?!  (real)
C 2>  END

      DATA
     +   BNAME(QENDT), BNEED(QENDT),
     +   (TYPE(J, QENDT), NEED(J, QENDT), KEY(J, QENDT), J = 1, 2)
     + / 'SPECIFY THE CAP TEMPERATURES', .FALSE.,
     +   'R', .FALSE., 'TEMPERATURE OF ENTRANCE CAP =',
     +   'R', .FALSE., 'TEMPERATURE OF EXIT CAP ='/

C///  SPECIFY THE CAP POWERS.
C 1>  SPECIFY THE CAP POWERS.
C 2>       POWER TO ENTRANCE CAP = ?!  (real)
C 2>       POWER TO EXIT CAP = ?!  (real)
C 2>  END

      DATA
     +   BNAME(QENDP), BNEED(QENDP),
     +   (TYPE(J, QENDP), NEED(J, QENDP), KEY(J, QENDP), J = 1, 2)
     + / 'SPECIFY THE CAP POWERS', .FALSE.,
     +   'R', .FALSE., 'POWER TO ENTRANCE CAP =',
     +   'R', .FALSE., 'POWER TO EXIT CAP ='/

      ELSE IF (WORD .EQ. 'CAP') THEN
         CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)
         IF (ERROR) GO TO 9301

         IF (WORD .EQ. 'TEMPERATURES') THEN
            Q = QENDT
            ERROR = BFOUND(Q)
            IF (ERROR) GO TO 9304
            BFOUND(Q) = .TRUE.
            ERROR = BFOUND(QENDP)
            IF (ERROR) GO TO 9307

C     SUBROUTINE READV
C    +  (ERROR, TEXT,
C    +   CKEY, ESIGN, FOUND, KEY, KEYS, IKEY, LINE, LKEY, NEED,
C    +   NUMBER, PRINT, RKEY, SCRIPT, TYPE)

C*****DOUBLE PRECISION
         CALL READV2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C         CALL READV1
C*****END SINGLE PRECISION
     +     (ERROR, TEXT,
     +      CKEY, .FALSE., FOUND, KEY(1, Q), 2, IKEY, LINE, LKEY,
     +      NEED(1, Q), NUMBER, PRINT, RKEY, SCRIPT, TYPE(1, Q))
         IF (ERROR) GO TO 9305

         IF (FOUND(1)) THEN
            TCAP(1) = RKEY(1)
            TCAPS(1) = .TRUE.
         ELSE
         ENDIF
         IF (FOUND(2)) THEN
            TCAP(2) = RKEY(2)
            TCAPS(2) = .TRUE.
         ELSE
         ENDIF

         ELSE IF (WORD .EQ. 'POWERS') THEN
            Q = QENDP
            ERROR = BFOUND(Q)
            IF (ERROR) GO TO 9304
            BFOUND(Q) = .TRUE.
            ERROR = BFOUND(QENDT)
            IF (ERROR) GO TO 9307

C     SUBROUTINE READV
C    +  (ERROR, TEXT,
C    +   CKEY, ESIGN, FOUND, KEY, KEYS, IKEY, LINE, LKEY, NEED,
C    +   NUMBER, PRINT, RKEY, SCRIPT, TYPE)

C*****DOUBLE PRECISION
         CALL READV2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C         CALL READV1
C*****END SINGLE PRECISION
     +     (ERROR, TEXT,
     +      CKEY, .FALSE., FOUND, KEY(1, Q), 2, IKEY, LINE, LKEY,
     +      NEED(1, Q), NUMBER, PRINT, RKEY, SCRIPT, TYPE(1, Q))
         IF (ERROR) GO TO 9305

         IF (FOUND(1)) THEN
            PCAP(1) = RKEY(1)
            TCAPS(1) = .FALSE.
         ELSE
         ENDIF
         IF (FOUND(2)) THEN
            PCAP(2) = RKEY(2)
            TCAPS(2) = .FALSE.
         ELSE
         ENDIF
      ELSE
         ERROR = .TRUE.
         GO TO 9603
      ENDIF

C///  SPECIFY THE THERMOCOUPLE LOCATIONS.
C 1>  SPECIFY THE THERMOCOUPLE LOCATIONS
C 2>       Z(1) = ?!  (real)
C 2>       Z(2) = ?!  (real)
C 2>       Z(3) = ?!  (real)
C 2>       Z(4) = ?!  (real)
C 2>       Z(5) = ?!  (real)
C 2>  END

      DATA
     +   BNAME(QTC), BNEED(QTC),
     +   (TYPE(J, QTC), NEED(J, QTC), KEY(J, QTC), J = 1, 5)
     + / 'SPECIFY THE THERMOCOUPLE LOCATIONS', .FALSE.,
     +   'R', .TRUE., 'Z(1) =',
     +   'R', .TRUE., 'Z(2) =',
     +   'R', .TRUE., 'Z(3) =',
     +   'R', .TRUE., 'Z(4) =',
     +   'R', .TRUE., 'Z(5) ='/

      ELSE IF (WORD .EQ. 'THERMOCOUPLE') THEN
         CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT,
     +      'LOCATIONS')
         IF (ERROR) GO TO 9301

         Q = QTC
         ERROR = BFOUND(Q)
         IF (ERROR) GO TO 9304
         BFOUND(Q) = .TRUE.

C     SUBROUTINE READV
C    +  (ERROR, TEXT,
C    +   CKEY, ESIGN, FOUND, KEY, KEYS, IKEY, LINE, LKEY, NEED,
C    +   NUMBER, PRINT, RKEY, SCRIPT, TYPE)

C*****DOUBLE PRECISION
         CALL READV2
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C         CALL READV1
C*****END SINGLE PRECISION
     +     (ERROR, TEXT,
     +      CKEY, .FALSE., FOUND, KEY(1, Q), 5, IKEY, LINE, LKEY,
     +      NEED(1, Q), NUMBER, PRINT, RKEY, SCRIPT, TYPE(1, Q))
         IF (ERROR) GO TO 9305

       DO 3040 J = 1, NUMTC
         IF (FOUND(J)) R(QZTC + J - 1) = RKEY(J)
3040   CONTINUE

C///  BOTTOM OF THE SPECIFICATION BLOCKS.

      END IF

C///  WRITE THE COEFFICIENT MATRIX.
C 1>  WRITE THE COEFFICIENT MATRIX

      ELSE IF (WORD .EQ. 'WRITE') THEN

      CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT,
     +   'THE COEFFICIENT MATRIX')
      IF (ERROR) GO TO 9301

      MATRIX = .TRUE.

C-----------------------------------------------

C///  FIND THE STEADY-STATE SOLUTION.
C 1>  FIND THE STEADY-STATE SOLUTION

      DATA
     +   BNAME(QSS), BNEED(QSS)
     + / 'FIND THE STEADY-STATE SOLUTION', .FALSE./

      ELSE IF (WORD .EQ. 'FIND') THEN

      CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, 'THE')
      IF (ERROR) GO TO 9301

      CALL READW (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT, WORD)
      IF (ERROR) GO TO 9302

      IF (WORD .EQ. 'STEADY-STATE') THEN

         CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT,
     +   'SOLUTION')
         IF (ERROR) GO TO 9301

         Q = QSS
         ERROR = BFOUND(Q)
         IF (ERROR) GO TO 9304
         BFOUND(Q) = .TRUE.

C///  FIND THE TRANSIENT SOLUTION.
C 1>  FIND THE TRANSIENT SOLUTION

      DATA
     +   BNAME(QTRANS), BNEED(QTRANS)
     + / 'FIND THE TRANSIENT SOLUTION', .FALSE./

      ELSE IF (WORD .EQ. 'TRANSIENT') THEN

         CALL VERIFY (ERROR, TEXT, LINE, NUMBER, PRINT, SCRIPT,
     +   'SOLUTION')
         IF (ERROR) GO TO 9301

         Q = QTRANS
         ERROR = BFOUND(Q)
         IF (ERROR) GO TO 9304
         BFOUND(Q) = .TRUE.

      END IF

C-----------------------------------------------

C///  BOTTOM OF THE LOOP THROUGH THE SCRIPT FILE.
C 1>  END

      ELSE IF (WORD .EQ. 'END') THEN
         GO TO 3050
      ELSE
         ERROR = .TRUE.
         GO TO 9312
      END IF
      GO TO 3030
3050  CONTINUE

C///  CHECK THAT ALL SPECIFICATION BLOCKS HAVE BEEN FOUND.

      DO 3060 J = 1, BLOCKS
         ERROR = BNEED(J) .AND. .NOT. BFOUND(J)
         IF (ERROR) GO TO 9313

C///  CHECK FOR ERRORS IN BLOCK SPECIFICATIONS.

      ERROR = (BFOUND(QSS) .AND. BFOUND(QTRANS))
         IF (ERROR) GO TO 9314

      ERROR = (BFOUND(QENDT) .AND. BFOUND(QENDP))
         IF (ERROR) GO TO 9326

3060  CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (4) PARTITION AND INITIALIZE THE WORK SPACE.
C
C///////////////////////////////////////////////////////////////////////

C///  DEFINE DIMENSIONAL PARAMETERS.

C     RES:  NUMBER OF ADDITIONAL UNKNOWNS FOR TRANSIENT SOLUTION.

      NADD = 0

C     TWAF1: NUMEX, NUMBER OF EXTERNAL TEMPERATURE POINTS (WALL,
C            WALL INSULATION, AIR FLOW).

      NUMEX = 4

      IF (CFLOW(1).OR.CFLOW(2)) NUMEX = NUMEX + 2

C     TWAF0: PDIM, MAXIMUM DIMENSION OF ENTRANCE OR EXIT REGION

      PDIM = MAX0(PIN, PEX)

C     TWOPNT: PMAX, TOTAL NUMBER OF COMPUTATIONAL CELLS

      PMAX = NUMWF

C     POINTS = PMAX
      POINTS = NUMWF

C     SPACING BETWEEN WAFERS

      IF (RBOAT) THEN
C          SPACE = (ZLAST - ZFIRST - REAL(NUMWF -
C     +            1)*(THICK+R(QRING+1)))/REAL(NUMWF - 1)
C
C     USE OLD SPACING FOR NOW
C
          SPACE = (ZLAST - ZFIRST)/REAL(NUMWF - 1)
      ELSE
          SPACE = (ZLAST - ZFIRST - REAL(NUMWF -
     +            1)*THICK)/REAL(NUMWF - 1)
          SPACE = (ZLAST - ZFIRST)/REAL(NUMWF - 1)
      ENDIF

C     TWOPNT: NUMBER OF ELEMENTS IN GROUPA AND GROUPB

C      GROUPA = PIN*(NUMEX + 1) + 3
C      GROUPB = PEX*(NUMEX + 1) + 3

      GROUPA = PIN*(NUMEX) + 2
      GROUPB = PEX*(NUMEX) + 2

C     TWOPNT: COMPS, SOLUTION COMPONENTS PER MESH POINT
C     (RADIAL POINTS IN A WAFER ARE TREATED AS COMPS AT
C     EACH AXIAL POINT PLUS ANY POINTS IN THE WALL AND INSULATION)

      COMPS = RPNTS + NUMEX

C     TWOPNT: ISIZE, SIZE OF THE INTEGER WORKSPACE

      XISIZE = 3*PMAX

C     TWOPNT: RSIZE, SIZE OF THE REAL WORKSPACE

      XRSIZE = 3*PMAX+9*(GROUPA + COMPS*PMAX + GROUPB)

C     TWPREP, TWSOLV: AMAX, DIMENSION OF PACKED MATRIX

      AMAX = (3*(COMPS + AMAX0(COMPS, GROUPA, GROUPB) - 1)
     +        + 2)*(GROUPA + COMPS*PMAX + GROUPB)

C     TWAF2: COLMAX, COLUMNS OF OUTPUT

      COLMAX = MAX(NUMWF+1, 8)

C     TWAF0: SIZE OF THE A1 AND A2 MATRIX FOR RADIOSITIES

      A2DIM = 4

C///  PARTITION THE WORK SPACE.

C     ACTIVE(COMPS)
C     LOGICAL ACTIVE

      CALL RESERV (ERROR, TEXT, .FALSE., LLAST, LMAX, LSIZE,
     +   'QACTIV', COMPS, QACTIV)
      IF (ERROR) GO TO 9101

C     AEX(5*PEX+2, 5*PEX+2)
C     REAL TWAF0, TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QAEX', (5*PEX+2)*(5*PEX+2), QAEX)
      IF (ERROR) GO TO 9101

C     AIN(5*PIN+2, 5*PIN+2)
C     REAL TWAF0, TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QAIN', (5*PIN+2)*(5*PIN+2), QAIN)
      IF (ERROR) GO TO 9101

C     AREAC(RPNTS+1)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QAREAC', RPNTS+1, QAREAC)
      IF (ERROR) GO TO 9101

C     A2(A2DIM, A2DIM)
C     REAL TWAF0, TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QA2', A2DIM*A2DIM, QA2)
      IF (ERROR) GO TO 9101

C     A2DIM
C     INTEGER TWAF0, TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QA2DIM', 1, QA2DIM)
      IF (ERROR) GO TO 9101

C     BUFFER(GROUPA + COMPS*PMAX + GROUPB + NADD)
C     REAL TWOPNT, TWAF1, TWPREP, TWSOLV

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QBUF', GROUPA + COMPS*PMAX + GROUPB + NADD, QBUF)
      IF (ERROR) GO TO 9101

C     CFLOW(2)
C     LOGICAL CFLOW

      CALL RESERV (ERROR, TEXT, .FALSE., LLAST, LMAX, LSIZE,
     +   'QCFLO', 2, QCFLO)
      IF (ERROR) GO TO 9101

C     CHOICE
C     INTEGER TWAFO, TWAF1, TWAF5

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QCHOIC', 1, QCHOIC)
      IF (ERROR) GO TO 9101

C     COMPS
C     INTEGER TWAFO, TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QCOMPS', 1, QCOMPS)
      IF (ERROR) GO TO 9101

C     CONDUC(8)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QCOND', 8, QCOND)
      IF (ERROR) GO TO 9101

C     CONW(RPNTS+1)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QCONW', RPNTS+1, QCONW)
      IF (ERROR) GO TO 9101

C     COORD(NUMWF)
C     REAL TWAF3

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QCOOR', NUMWF, QCOOR)
      IF (ERROR) GO TO 9101

C     CP(8)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QCP', 8, QCP)
      IF (ERROR) GO TO 9101

C     DZEND(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QDZEND', 2, QDZEND)
      IF (ERROR) GO TO 9101

C     DZENDI(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QDZNDI', 2, QDZNDI)
      IF (ERROR) GO TO 9101

C     EMIS(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QEMIS', 2, QEMIS)
      IF (ERROR) GO TO 9101

C     EMISS
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QEMISS', 1, QEMISS)
      IF (ERROR) GO TO 9101

C     EMISSE
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QEMISE', 1, QEMISE)
      IF (ERROR) GO TO 9101

C     EMISSI(5)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QEMISI', 5, QEMISI)
      IF (ERROR) GO TO 9101

C     EMISSW
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QEMISW', 1, QEMISW)
      IF (ERROR) GO TO 9101

C     F(GROUPA + COMPS*PMAX + GROUPB + NADD)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QF', GROUPA + COMPS*PMAX + GROUPB + NADD, QF)
      IF (ERROR) GO TO 9101

C     FEE(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFEE', 2, QFEE)
      IF (ERROR) GO TO 9101

C     FEG(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFEG', 2, QFEG)
      IF (ERROR) GO TO 9101

C     FENDS(PDIM + 1, NUMEX + 1, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFENDS', (PDIM + 1)*(NUMEX + 1)*2, QFENDS)
      IF (ERROR) GO TO 9101

C     FENG(NUMWF*(RPNTS+NUMEX))
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFENG', NUMWF*(RPNTS+NUMEX), QFENG)
      IF (ERROR) GO TO 9101

C     FESE(PDIM, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFESE', 2*PDIM, QFESE)
      IF (ERROR) GO TO 9101

C     FEW(RPNTS, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFEW', 2*RPNTS, QFEW)
      IF (ERROR) GO TO 9101

C     FGE(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFGE', 2, QFGE)
      IF (ERROR) GO TO 9101

C     FGSE(PDIM, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFGSE', 2*PDIM, QFGSE)
      IF (ERROR) GO TO 9101

C     FSEE(PDIM, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFSEE', 2*PDIM, QFSEE)
      IF (ERROR) GO TO 9101

C     FSEG(PDIM, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFSEG', 2*PDIM, QFSEG)
      IF (ERROR) GO TO 9101

C     FSES(PDIM, PDIM, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFSES', 2*PDIM*PDIM, QFSES)
      IF (ERROR) GO TO 9101

C     FSEW(PDIM, RPNTS, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFSEW', 2*PDIM*RPNTS, QFSEW)
      IF (ERROR) GO TO 9101

C     FSS
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFSS', 1, QFSS)
      IF (ERROR) GO TO 9101

C     FSW(2*RPNTS)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFSW', 2*RPNTS, QFSW)
      IF (ERROR) GO TO 9101

C     FWE(RPNTS, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFWE', 2*RPNTS, QFWE)
      IF (ERROR) GO TO 9101

C     FWS(2*RPNTS)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFWS', 2*RPNTS, QFWS)
      IF (ERROR) GO TO 9101

C     FWSE(RPNTS, PDIM, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFWSE', 2*RPNTS*PDIM, QFWSE)
      IF (ERROR) GO TO 9101

C     FWW(2*RPNTS, 2*RPNTS)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QFWW', 4*RPNTS*RPNTS, QFWW)
      IF (ERROR) GO TO 9101

C     GROUPA
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QGROPA', 1, QGROPA)
      IF (ERROR) GO TO 9101

C     GROUPB
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QGROPB', 1, QGROPB)
      IF (ERROR) GO TO 9101

C     HEADER(5, COLMAX)
C     CHARACTER TWAF2

      CALL RESERV (ERROR, TEXT, .FALSE., CLAST, CMAX, CSIZE,
     +   'QHEAD', 5*COLMAX, QHEAD)
      IF (ERROR) GO TO 9101

C     HCOEF(5)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QHCOEF', 5, QHCOEF)
      IF (ERROR) GO TO 9101

C     IBEG(3)
C     INTEGER TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QIBEG', 3, QIBEG)
      IF (ERROR) GO TO 9101

C     IEND(3)
C     INTEGER TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QIEND', 3, QIEND)
      IF (ERROR) GO TO 9101

C     IPVT2(A2DIM)
C     INTEGER TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QIPV2', A2DIM, QIPV2)
      IF (ERROR) GO TO 9101

C     IPVT3(5*PIN+2)
C     INTEGER TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QIPV3', 5*PIN+2, QIPV3)
      IF (ERROR) GO TO 9101

C     IPVT4(5*PEX+2)
C     INTEGER TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QIPV4', 5*PEX+2, QIPV4)
      IF (ERROR) GO TO 9101

C     IWORK(XISIZE)
C     INTEGER TWOPNT

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QIWORK', XISIZE, QIWORK)
      IF (ERROR) GO TO 9101

C     J1(4, PDIM, 2)
C     REAL TWAF1, TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QJ1', 4*PDIM*2, QJ1)
      IF (ERROR) GO TO 9101

C     J2(A2DIM, NUMWF)
C     REAL TWAF1, TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QJ2', A2DIM*NUMWF, QJ2)
      IF (ERROR) GO TO 9101

C     JEND(2)
C     REAL TWAF1, TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QJEND', 2, QJEND)
      IF (ERROR) GO TO 9101

C     JEX(5*PEX + 2)
C     REAL TWAF1, TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QJEX', 5*PEX + 2, QJEX)
      IF (ERROR) GO TO 9101

C     JIN(5*PIN + 2)
C     REAL TWAF1, TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QJIN', 5*PIN + 2, QJIN)
      IF (ERROR) GO TO 9101

C     JWALL(PDIM, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QJWALL', PDIM*2, QJWALL)
      IF (ERROR) GO TO 9101

C     LOS1(NUMWF)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QLOS1', NUMWF, QLOS1)
      IF (ERROR) GO TO 9101

C     LOS2(PDIM, 2)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QLOS2', PDIM*2, QLOS2)
      IF (ERROR) GO TO 9101

C     LOS3(2)
C     REAL WAFER

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QLOS3', 2, QLOS3)
      IF (ERROR) GO TO 9101

C     MARK(PMAX)
C     LOGICAL TWOPNT

      CALL RESERV (ERROR, TEXT, .FALSE., LLAST, LMAX, LSIZE,
     +   'QMARK', PMAX, QMARK)
      IF (ERROR) GO TO 9101

C     MAXHT
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QMAXHT', 1, QMAXHT)
      IF (ERROR) GO TO 9101

C     MAXIN
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QMAXIN', 1, QMAXIN)
      IF (ERROR) GO TO 9101

C     MFLOW(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QMFLO', 2, QMFLO)
      IF (ERROR) GO TO 9101

C     MINJ(NUMWF)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QMINJ', NUMWF, QMINJ)
      IF (ERROR) GO TO 9101

C     MINJE(NUMWF)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QMINJE', NUMWF, QMINJE)
      IF (ERROR) GO TO 9101

C     MWZ(PDIM+1, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QMWZ', (PDIM+1)*2, QMWZ)
      IF (ERROR) GO TO 9101

C     MWZW(MUMWF+1)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QMWZW', NUMWF+1, QMWZW)
      IF (ERROR) GO TO 9101

C     NADD
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QNADD', 1, QNADD)
      IF (ERROR) GO TO 9101

C     NUMEX
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QNUMEX', 1, QNUMEX)
      IF (ERROR) GO TO 9101

C     NUMHT
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QNUMHT', 1, QNUMHT)
      IF (ERROR) GO TO 9101

C     NUMIN
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QNUMIN', 1, QNUMIN)
      IF (ERROR) GO TO 9101

C     NUMT
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QNUMT', 1, QNUMT)
      IF (ERROR) GO TO 9101

C     NUMTC
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QNUMTC', 1, QNUMTC)
      IF (ERROR) GO TO 9101

C     NUMWF
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QNUMWF', 1, QNUMWF)
      IF (ERROR) GO TO 9101

C     PCAP(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QPCAP', 2, QPCAP)
      IF (ERROR) GO TO 9101

C     PDIM
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QPDIM', 1, QPDIM)
      IF (ERROR) GO TO 9101

C     PEX
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QPEX', 1, QPEX)
      IF (ERROR) GO TO 9101

C     PIN
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QPIN', 1, QPIN)
      IF (ERROR) GO TO 9101

C     PIVOT(GROUPA + COMPS*PMAX + GROUPB)
C     INTEGER TWPREP TWSOLV

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QPIVOT', GROUPA + COMPS*PMAX + GROUPB, QPIVOT)
      IF (ERROR) GO TO 9101

C     PMAX
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QPMAX', 1, QPMAX)
      IF (ERROR) GO TO 9101

C     POWER(NUMWF)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QPOWER', NUMWF, QPOWER)
      IF (ERROR) GO TO 9101

C     PWALL(PDIM, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QPWALL', PDIM*2, QPWALL)
      IF (ERROR) GO TO 9101

C     QWAF(RPNTS, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QQWAF', RPNTS*2, QQWAF)
      IF (ERROR) GO TO 9101

C     R(RPNTS)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QR', RPNTS, QR)
      IF (ERROR) GO TO 9101

C     REFL(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QREFL', 2, QREFL)
      IF (ERROR) GO TO 9101

C     REXT(NUMEX)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QREXT', NUMEX, QREXT)
      IF (ERROR) GO TO 9101

C     RCPV(RPNTS)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QRCPV', RPNTS, QRCPV)
      IF (ERROR) GO TO 9101

C     RHO(8)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QRHO', 8, QRHO)
      IF (ERROR) GO TO 9101

C     RINQ(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QRINQ', 2, QRINQ)
      IF (ERROR) GO TO 9101

C     RINJEC
C     LOGICAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., LLAST, LMAX, LSIZE,
     +   'QRINJE', 2, QRINJE)
      IF (ERROR) GO TO 9101

C     ROUTQ(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QROUTQ', 2, QROUTQ)
      IF (ERROR) GO TO 9101

C     RPNTS
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QRPNTS', 1, QRPNTS)
      IF (ERROR) GO TO 9101

C     RTUBEI
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QRTUBI', 1, QRTUBI)
      IF (ERROR) GO TO 9101

C     RTUBEO
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QRTUBO', 1, QRTUBO)
      IF (ERROR) GO TO 9101

C     R(XRSIZE)
C     REAL TWOPNT

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QRWORK', XRSIZE, QRWORK)
      IF (ERROR) GO TO 9101

C     SPACE
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QSPACE', 1, QSPACE)
      IF (ERROR) GO TO 9101

C     TAMB
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTAMB', 1, QTAMB)
      IF (ERROR) GO TO 9101

C     TCAP(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTCAP', 2, QTCAP)
      IF (ERROR) GO TO 9101

C     TCAPS(2)
C     LOGICAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., LLAST, LMAX, LSIZE,
     +   'QTCAPS', 2, QTCAPS)
      IF (ERROR) GO TO 9101

C     TEAST(NUMWF)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTEAST', NUMWF, QTEAST)
      IF (ERROR) GO TO 9101

C     TEMPER(RPNTS+NUMEX, NUMWF)
C     REAL TWAF3

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTEMPE', NUMWF*(RPNTS+NUMEX), QTEMPE)
      IF (ERROR) GO TO 9101

C     TEMMAX
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QTEMAX', 1, QTEMAX)
      IF (ERROR) GO TO 9101

C     THICK
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTHICK', 1, QTHICK)
      IF (ERROR) GO TO 9101

C     THICKI
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTHCKI', 1, QTHCKI)
      IF (ERROR) GO TO 9101

C     TIME
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QTIME', 1, QTIME)
      IF (ERROR) GO TO 9101

C     TOUT(2*PDIM + NUMWF + 2)
C     REAL TWAF4

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTOUT2', 2*PDIM + NUMWF + 2, QTOUT2)
      IF (ERROR) GO TO 9101

C     TRAN(2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTRAN', 2, QTRAN)
      IF (ERROR) GO TO 9101

C     TWAF(NUMWF+1*RPNTS)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTWAF', (NUMWF+1)*RPNTS, QTWAF)
      IF (ERROR) GO TO 9101

C     TWAF0(NUWWF+1*RPNTS)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTWAF0', (NUMWF+1)*RPNTS, QTWAF0)
      IF (ERROR) GO TO 9101

C     TWAF4(NUMWF+1*RPNTS)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTWAF4', (NUMWF+1)*RPNTS, QTWAF4)
      IF (ERROR) GO TO 9101

C     TEX(NUMWF, NUMEX)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTEX', NUMWF*NUMEX, QTEX)
      IF (ERROR) GO TO 9101

C     TEX0(NUWWF, NUMEX)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTEX0', NUMWF*NUMEX, QTEX0)
      IF (ERROR) GO TO 9101

C     TEX4(NUMWF, NUMEX)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTEX4', NUMWF*NUMEX, QTEX4)
      IF (ERROR) GO TO 9101

C     TEXT
C     INTEGER TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QTEXT', 1, QTEXT)
      IF (ERROR) GO TO 9101

C     TMID(NUMWF)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTMID', NUMWF, QTMID)
      IF (ERROR) GO TO 9101

C     TSIDE(PDIM, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTSIDE', 2*PDIM, QTSIDE)
      IF (ERROR) GO TO 9101

C     TSTEP
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTSTEP', 1, QTSTEP)
      IF (ERROR) GO TO 9101

C     TTGAS
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTTGAS', 2, QTTGAS)
      IF (ERROR) GO TO 9101

C     TWEST(NUMWF)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTWEST', NUMWF, QTWEST)
      IF (ERROR) GO TO 9101

C     TX(PDIM, NUMEX, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTX', PDIM*NUMEX*2, QTX)
      IF (ERROR) GO TO 9101

C     TX0(PDIM, NUMEX, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTX0', PDIM*NUMEX*2, QTX0)
      IF (ERROR) GO TO 9101

C     TX4(PDIM, NUMEX, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QTX4', PDIM*NUMEX*2, QTX4)
      IF (ERROR) GO TO 9101

C     VALUE(RPNTS+NUMEX, COLMAX)
C     REAL TWAF2

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QVALUE', (RPNTS+NUMEX)*COLMAX, QVALUE)
      IF (ERROR) GO TO 9101

C     VAL2(PDIM, COLMAX)
C     REAL VALUE

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QVAL2', PDIM*COLMAX, QVAL2)
      IF (ERROR) GO TO 9101

C     VAL3(PMAX)
C     REAL VALUE

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QVAL3', PMAX, QVAL3)
      IF (ERROR) GO TO 9101

C     VINJ(MAXIN)
C     REAL VALUE

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QVINJ', MAXIN, QVINJ)
      IF (ERROR) GO TO 9101

C     WORK2(A2DIM)
C     REAL TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QWOR2', A2DIM, QWOR2)
      IF (ERROR) GO TO 9101

C     WORK3(5*PIN+2)
C     REAL TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QWOR3', 5*PIN+2, QWOR3)
      IF (ERROR) GO TO 9101

C     WORK4(5*PEX+2)
C     REAL TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QWOR4', 5*PEX+2, QWOR4)
      IF (ERROR) GO TO 9101

C     WR(RPNTS + 1)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QWR', RPNTS + 1, QWR)
      IF (ERROR) GO TO 9101

C     WZ(PDIM + 1, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QWZ', 2*(PDIM + 1), QWZ)
      IF (ERROR) GO TO 9101

C     WZWAF(NUMWF + 1)
C     REAL TWAF1, TWAF0

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QWZWAF', NUMWF+1, QWZWAF)
      IF (ERROR) GO TO 9101

C     X(GROUPA + COMPS*PMAX + GROUPB + NADD)
C     REAL TWAF1 TWOPNT

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QX', GROUPA + COMPS*PMAX + GROUPB + NADD, QX)
      IF (ERROR) GO TO 9101

C     X0(GROUPA + COMPS*PMAX + GROUPB + NADD)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QX0', GROUPA + COMPS*PMAX + GROUPB + NADD, QX0)
      IF (ERROR) GO TO 9101

C     Z(PDIM, 2)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QZ', 2*PDIM, QZ)
      IF (ERROR) GO TO 9101

C     ZBEGI(MAXIN)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QZBEGI', MAXIN, QZBEGI)
      IF (ERROR) GO TO 9101

C     ZENDI(MAXIN)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QZENDI', MAXIN, QZENDI)
      IF (ERROR) GO TO 9101

C     ZLEN
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QZLEN', 1, QZLEN)
      IF (ERROR) GO TO 9101

C     ZOUT(2*PDIM + NUMWF + 2)
C     REAL TWAF4

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QZOUT', 2*PDIM + NUMWF + 2, QZOUT)
      IF (ERROR) GO TO 9101

C     ZWAF(NUMWF)
C     REAL TWAF1

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QZWAF', NUMWF, QZWAF)
      IF (ERROR) GO TO 9101

      ERROR = .NOT. (CLAST .LE. CSIZE .AND. ILAST .LE. ISIZE .AND.
     +   LLAST .LE. LSIZE .AND. RLAST .LE. RSIZE)
      IF (ERROR) GO TO 9401

C///  SET FIXED HEATER TEMPERATURE FLAG.

      TFIX = BFOUND(QTEMP)

C///  TWOPNT: ACTIVE

      LWORK(QACTIV + 0) = .FALSE.

C///  CALCULATE SHAPE FACTORS, GRID, AND INITIAL GUESS.

      SHAPE = .TRUE.

C///  STORE COOLING FLAGS.

      LWORK(QCFLO) = CFLOW(1)
      LWORK(QCFLO+1) = CFLOW(2)

      LWORK(QRINJE) = RINJEC

C///  STORE SOME IN PARAMETERS IN IPAR AND RPAR.

      IWORK(QA2DIM) = A2DIM
      R(QTAMB) = TAMB
      IWORK(QCOMPS) = COMPS
      DO 100 I = 1, 8
         R(QCOND + I-1) = CONDUC(I)
         R(QCP + I-1) = CP(I)
         R(QRHO + I-1) = RHO(I)
100   CONTINUE
      DO 200 I = 1, 2
         R(QDZEND + I-1) = DZEND(I)
         R(QDZNDI + I-1) = DZENDI(I)
         R(QEMIS + I-1) = EMIS(I)
         R(QMFLO + I-1) = MFLOW(I)
         R(QPCAP + I-1) = PCAP(I)
         R(QREFL + I-1) = REFL(I)
         R(QRINQ + I-1) = RINQ(I)
         R(QROUTQ + I-1) = ROUTQ(I)
         R(QTCAP + I-1) = TCAP(I)
         R(QTRAN + I-1) = TRAN(I)
         R(QTTGAS + I -1) = TTGAS(I)
         LWORK(QTCAPS + I-1) = TCAPS(I)
200   CONTINUE
      DO 300 I = 1, 5
         R(QEMISI + I-1) = EMISSI(I)
         R(QHCOEF + I-1) = HCOEF(I)
300   CONTINUE
      DO 400 I = 1, MAXIN
         R(QVINJ + I-1) = VINJ(I)
         R(QZBEGI + I-1) = ZBEGI(I)
         R(QZENDI + I-1) = ZENDI(I)
400   CONTINUE
      R(QEMISS) = EMISS
      R(QEMISE) = EMISSE
      R(QEMISW) = EMISSW

      IWORK(QGROPA) = GROUPA
      IWORK(QGROPB) = GROUPB
      IWORK(QMAXHT) = MAXHT
      IWORK(QMAXIN) = MAXIN
      IWORK(QNADD) = NADD
      IWORK(QNUMEX) = NUMEX
      IWORK(QNUMHT) = NUMHT
      IWORK(QNUMIN) = NUMIN
      IWORK(QNUMT) = NUMT
      IWORK(QNUMTC) = NUMTC
      IWORK(QNUMWF) = NUMWF
      IWORK(QPDIM) = PDIM
      IWORK(QPEX) = PEX
      IWORK(QPIN) = PIN
      IWORK(QPMAX) = PMAX
      IWORK(QRPNTS) = RPNTS
      IWORK(QTEXT) = TEXT
      R(QRTUBI) = RTUBEI
      R(QRTUBO) = RTUBEO
      R(QSPACE) = SPACE
      R(QTHICK) = THICK
      R(QTHCKI) = THICKI
      IWORK(QTEMAX) = TEMMAX
      R(QTSTEP) = TSTEP
      R(QZLEN) = ZLEN

C     SUBROUTINE TWAF0
C    +  (ERROR, TEXT,
C    +   A2, A2DIM, AEX, AIN, CFLOW, COMPS, EMIS, EMISS,
C    +   EMISSE, EMISSW,
C    +   FEE, FEG, FESE, FEW, FGE, FGSE, FSEE, FSEG, FSES, FSEW,
C    +   FSS, FSW, FWE, FWS, FWSE, FWW, GROUPA, GROUPB, IBEG, IEND,
C    +   IPVT2, IPVT3, IPVT4, JEND, JWALL, NUMEX, NUMT,
C    +   NUMWF, PDIM, PEX, PIN, PMAX, R, RADIUS, RBOAT, REFL,
C    +   RING, RINQ, ROUTQ, RPNTS, RTUBEI, S, S0, SHAPE,
C    +   SPACE, TDAT, TEMMAX, TFIX, TMID, TRAN, TSIDE, TTEND,
C    +   TTGAS, TTIN, TTQ, TTUBE, TTWAF, WORK2, WORK3, WORK4,
C    +   WR, WZ, WZWAF, Z, ZDAT, ZFIRST, ZLAST, ZLEN, ZWAF)

      CALL TWAF0
     +  (ERROR, TEXT,
     +   R(QA2), A2DIM, R(QAEX), R(QAIN), CFLOW, COMPS, EMIS,
     +   EMISS, EMISSE, EMISSW, R(QFEE), R(QFEG), R(QFESE), R(QFEW),
     +   R(QFGE), R(QFGSE),
     +   R(QFSEE), R(QFSEG), R(QFSES), R(QFSEW), R(QFSS), R(QFSW),
     +   R(QFWE),
     +   R(QFWS), R(QFWSE), R(QFWW), GROUPA, GROUPB, IWORK(QIBEG),
     +   IWORK(QIEND), IWORK(QIPV2), IWORK(QIPV3),
     +   IWORK(QIPV4), R(QJEND), R(QJWALL), NUMEX, NUMT, NUMWF,
     +   PDIM, PEX, PIN, PMAX, R(QR), RADIUS, RBOAT, REFL, R(QRING),
     +   RINQ, ROUTQ, RPNTS, RTUBEI, R(QX), R(QX0), SHAPE, SPACE,
     +   R(QTDAT), TEMMAX, TFIX, R(QTMID), TRAN, R(QTSIDE),
     +   TTEND, TTGAS, TTIN, TTQ, TTUBE, TTWAF, R(QWOR2), R(QWOR3),
     +   R(QWOR4), R(QWR), R(QWZ), R(QWZWAF), R(QZ), R(QZDAT), ZFIRST,
     +   ZLAST, ZLEN, R(QZWAF))
      IF (ERROR) GO TO 9406

C///  SET RADIAL POSITIONS OF QUARTZ TUBES, HEATERS, AND INSULATION.

      R(QREXT) = (RINQ(1) + ROUTQ(1))/2.
      R(QREXT+1) = (RINQ(2) + ROUTQ(2))/2.
      R(QREXT+2) = (RTUBEO + RTUBEI)/2.
      R(QREXT+3) = RTUBEO + THICKI/2.

C     GAS FLOWS
C      R(QREXT+4) = (ROUTQ(1) + RINQ(2))/2.
C      R(QREXT+5) = (ROUTQ(2) + RTUBEI)/2.

C     SET FAKE RADIAL POSITION FOR GAS FLOWS
      R(QREXT+4) = R(QREXT+3) + 0.5
      R(QREXT+5) = R(QREXT+4) + 0.5


C///  CHECK TRANSIENT OR STEADY-STATE SOLUTION.

      IF (BFOUND(QTRANS)) GO TO 77777

      ERROR = .NOT. (BFOUND(QSS))
      IF (ERROR) GO TO 9104

C///////////////////////////////////////////////////////////////////////
C
C     (5) STEADY-STATE SECTION (USING TWOPNT).
C
C///////////////////////////////////////////////////////////////////////

C     A(AMAX)
C     REAL TWPREP, TWSOLV

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QA', AMAX, QA)
      IF (ERROR) GO TO 9101

C     ABOVE(GROUPA + COMPS + GROUPB)
C     REAL TWOPNT

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QABOVE', GROUPA + COMPS + GROUPB, QABOVE)
      IF (ERROR) GO TO 9101

C     BELOW(GROUPA + COMPS + GROUPB)
C     REAL TWOPNT

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QBELOW', GROUPA + COMPS + GROUPB, QBELOW)
      IF (ERROR) GO TO 9101

C///  SET ABOVE AND BELOW

      OFFSET = NUMWF*(RPNTS + NUMEX) + GROUPA
      THIGH = 4000.0
      TLOW  = 0.0

      DO 5005 J = 0, GROUPA + COMPS + GROUPB -1
         R(QABOVE + J) = THIGH
         R(QBELOW + J) = TLOW
5005  CONTINUE

C///  CALL TWOPNT.

      IF (0 .LT. TEXT) WRITE (TEXT, 10002) ID, 'CALLING TWOPNT.'

C///  CHOOSE THE VERSION.

C*****DOUBLE PRECISION
      VERSIO = 'DOUBLE PRECISION VERSION 3.27'
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      VERSIO = 'SINGLE PRECISION VERSION 3.27'
C*****END SINGLE PRECISION

C///  SET THE CONTROLS.

      CALL TWINIT (ERROR, TEXT, .FALSE.)
      IF (ERROR) GO TO 9507

      CALL TWSETL (ERROR, TEXT, 'ADAPT', ADAPT)
      IF (ERROR) GO TO 9508

      CALL TWSETI (ERROR, TEXT, 'LEVELD', LEVEL(1))
      IF (ERROR) GO TO 9508

      CALL TWSETI (ERROR, TEXT, 'LEVELM', LEVEL(2))
      IF (ERROR) GO TO 9508

      CALL TWSETI (ERROR, TEXT, 'PADD', 0)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'SSABS', SSABS)
      IF (ERROR) GO TO 9508

      CALL TWSETI (ERROR, TEXT, 'SSAGE', SSAGE)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'SSREL', SSREL)
      IF (ERROR) GO TO 9508

      CALL TWSETL (ERROR, TEXT, 'STEADY', .TRUE.)
      IF (ERROR) GO TO 9508

      CALL TWSETI (ERROR, TEXT, 'STEPS0', STEPS0)
      IF (ERROR) GO TO 9508

      CALL TWSETI (ERROR, TEXT, 'STEPS1', STEPS1)
      IF (ERROR) GO TO 9508

      CALL TWSETI (ERROR, TEXT, 'STEPS2', STEPS2)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'STRID0', TSTEP0)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'TDABS', TDABS)
      IF (ERROR) GO TO 9508

      CALL TWSETI (ERROR, TEXT, 'TDAGE', TDAGE)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'TDEC', TDEC)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'TDREL', TDREL)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'TINC', TINC)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'TMAX', TMAX)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'TMIN', TMIN)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'TOLER0', TOLER0)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'TOLER1', TOLER1)
      IF (ERROR) GO TO 9508

      CALL TWSETR (ERROR, TEXT, 'TOLER2', TOLER2)
      IF (ERROR) GO TO 9508

C///  CALL TWOPNT.

C     SUBROUTINE TWOPNT
C    +  (ERROR, TEXT, VERSIO,
C    +   ABOVE, ACTIVE, BELOW, BUFFER, COMPS, CONDIT, GROUPA, GROUPB,
C    +   ISIZE, IWORK, MARK, NAME, NAMES, PMAX, POINTS, REPORT, RSIZE,
C    +   RWORK, SIGNAL, STRIDE, TIME, U, X)

      SIGNAL = ' '
5010  CONTINUE
      CALL TWOPNT
     +  (ERROR, TEXT, VERSIO,
     +   R(QABOVE), LWORK(QACTIV), R(QBELOW), R(QBUF),
     +   COMPS, CONDIT, GROUPA, GROUPB, XISIZE, IWORK(QIWORK),
     +   LWORK(QMARK), BLANK, 1, PMAX, POINTS, REPORT, XRSIZE,
     +   R(QRWORK), SIGNAL, TSTEP, TIME, R(QX), R(QZ))
      IF (ERROR) GO TO 9501

C///  TOP OF THE BLOCK TO SERVICE REQUESTS FROM TWOPNT.

      IF (SIGNAL .NE. ' ') THEN

C///  EVALUATE THE RESIDUAL.

      IF (SIGNAL .EQ. 'RESIDUAL' ) THEN

C      SUBROUTINE TWAF1
C     +  (ERROR, TEXT,
C     +   A2, A2DIM, AEX, AIN, AREAC, BTEMP, CFLOW, CHOICE,
C     +   COMPS, CONDUC,
C     +   CONW, CP, DZEND, DZENDI, EMIS, EMISS, EMISSE, EMISSI, EMISSW,
C     +   F, FEE, FEG, FENDS, FENG, FESE, FEW, FGE, FGSE, FSEE, FSEG,
C     +   FSES, FSEW, FSS,
C     +   FSW, FWE, FWS, FWSE, FWW, GROUPA, GROUPB, HCOEF, IBEG, IEND,
C     +   IPVT2, IPVT3, IPVT4, J1, J2, JEND, JEX, JIN, JWALL, LOS1,
C     +   LOS2, LOS3, MAXHT, MAXIN, MFLOW, MINJ, MINJE, MWZ, MWZW,
C     +   NUMEX, NUMHT, NUMIN, NUMT, NUMWF, PCAP, PDIM, PEX, PIN,
C     +   PMAX, POWER, POWERH, PWALL, QWAF, R, RBOAT, RCPV, REFL,
C     +   RHO, RING, RINJEC, RINQ, ROUTQ,
C     +   RPNTS, RTUBEI, RTUBEO,
C     +   S, S0, SHAPE, SPACE, TCAP, TCAPS, TDAT, TEAST, TEMMAX, TEX,
C     +   TEX0, TEX4, TFIX, THICK, THICKI, TIME, TMID, TRAN, TSIDE,
C     +   TSTEP, TTGAS, TWAF, TWAF0, TWAF4, TWEST, TX, TX0, TX4,
C     +   VINJ, WR, WZ, WZWAF,
C     +   Z, ZBEG, ZBEGI, ZDAT, ZEND, ZENDI, ZLEN, ZWAF,
C     +   GBOT, GTOP, JBOT, JTOP, NUMTC, SUMB, SUMT, TC, ZTC)

C     THE ARGUMENT LIST FOR TWAF1 IS THE SAME FOR ALL CALLS IN TWAFER,
C     BUT THE LIST IS DIFFERENT FOR THE CALLS IN TWAF5.  THE ARGUMENT
C     LIST FOR TWAF5 (EXCEPT FOR THE LAST LINE) IS THE SAME AS THE
C     ARGUMENT LIST FOR TWAF1.

      CHOICE = 0
      CALL TWAF1
     +  (ERROR, TEXT,
     +   R(QA2), IWORK(QA2DIM), R(QAEX), R(QAIN), R(QAREAC),
     +   R(QTAMB), LWORK(QCFLO), CHOICE,
     +   IWORK(QCOMPS), R(QCOND), R(QCONW), R(QCP),
     +   R(QDZEND), R(QDZNDI), R(QEMIS), R(QEMISS),
     +   R(QEMISE), R(QEMISI), R(QEMISW), R(QF),
     +   R(QFEE), R(QFEG), R(QFENDS), R(QFENG), R(QFESE), R(QFEW),
     +   R(QFGE), R(QFGSE), R(QFSEE), R(QFSEG), R(QFSES), R(QFSEW),
     +   R(QFSS), R(QFSW), R(QFWE), R(QFWS), R(QFWSE), R(QFWW),
     +   IWORK(QGROPA), IWORK(QGROPB), R(QHCOEF), IWORK(QIBEG),
     +   IWORK(QIEND), IWORK(QIPV2), IWORK(QIPV3), IWORK(QIPV4),
     +   R(QJ1), R(QJ2), R(QJEND), R(QJEX), R(QJIN),
     +   R(QJWALL), R(QLOS1), R(QLOS2), R(QLOS3), IWORK(QMAXHT),
     +   IWORK(QMAXIN), R(QMFLO), R(QMINJ), R(QMINJE), R(QMWZ),
     +   R(QMWZW),
     +   IWORK(QNUMEX), IWORK(QNUMHT), IWORK(QNUMIN), IWORK(QNUMT),
     +   IWORK(QNUMWF), R(QPCAP), IWORK(QPDIM), IWORK(QPEX),
     +   IWORK(QPIN), IWORK(QPMAX),
     +   R(QPOWER), R(QPOWRH), R(QPWALL), R(QQWAF), R(QR), RBOAT,
     +   R(QRCPV), R(QREFL), R(QRHO), R(QRING), LWORK(QRINJE),
     +   R(QRINQ), R(QROUTQ),
     +   IWORK(QRPNTS), R(QRTUBI), R(QRTUBO), R(QBUF), R(QX0), SHAPE,
     +   R(QSPACE), R(QTCAP), LWORK(QTCAPS), R(QTDAT), R(QTEAST),
     +   IWORK(QTEMAX), R(QTEX),
     +   R(QTEX0), R(QTEX4), TFIX, R(QTHICK), R(QTHCKI), TIME,
     +   R(QTMID), R(QTRAN), R(QTSIDE), R(QTSTEP), R(QTTGAS),
     +   R(QTWAF), R(QTWAF0),
     +   R(QTWAF4), R(QTWEST), R(QTX), R(QTX0), R(QTX4), R(QVINJ),
     +   R(QWR), R(QWZ),
     +   R(QWZWAF), R(QZ), R(QZBEG), R(QZBEGI), R(QZDAT), R(QZEND),
     +   R(QZENDI),  R(QZLEN),
     +   R(QZWAF), R(QGBOT), R(QGTOP), R(QJBOT), R(QJTOP),
     +   IWORK(QNUMTC), R(QSUMB), R(QSUMT), R(QTCT), R(QZTC))

      IF (ERROR) GO TO 9502

      DO 5020 J = 0, GROUPA + COMPS * POINTS + GROUPB - 1
         R(QBUF + J) = R(QF + J)
5020  CONTINUE

C///  EVALUATE AND FACTOR THE JACOBIAN.

      ELSE IF (SIGNAL .EQ. 'PREPARE' ) THEN

      FLAG = .FALSE.
5030  CONTINUE

C     SUBROUTINE TWPREP
C    +  (ERROR, TEXT,
C    +   A, ASIZE, BUFFER, COMPS, CONDIT, GROUPA, GROUPB, PIVOT, POINTS,
C    +   RETURN)

      CALL TWPREP
     +  (ERROR, TEXT,
     +   R(QA), AMAX, R(QBUF), COMPS, CONDIT, GROUPA, GROUPB,
     +   IWORK(QPIVOT), POINTS, FLAG)
      IF (ERROR) GO TO 9503

      IF (FLAG) THEN

C      SUBROUTINE TWAF1
C    +  (ERROR, TEXT,
C    +   A2, A2DIM, AEX, AIN, AREAC, BTEMP, CFLOW, CHOICE, COMPS,
C    +   CONDUC, CONW, CP, DZEND,
C    +   DZENDI, EMIS, EMISS, EMISSE, EMISSI, EMISSW, F, FEE, FEG,
C    +   FENDS, FENG, FESE, FEW, FGE, FGSE, FSEE, FSEG, FSES, FSEW,
C    +   FSS, FSW, FWE, FWS, FWSE, FWW, GROUPA, GROUPB, HCOEF, IBEG,
C    +   IEND, IPVT2, IPVT3, IPVT4, J1, J2, JEND, JEX, JIN, JWALL,
C    +   LOS1, LOS2, LOS3, MAXHT, MAXIN, MFLOW, MINJ, MINJE, MWZ,
C    +   MWZW, NUMEX, NUMHT, NUMIN, NUMT, NUMWF,
C    +   PCAP, PDIM, PEX, PIN, PMAX, POWER, POWERH,
C    +   PWALL, QWAF, R, RBOAT, RCPV, REFL, RHO, RING, RINJEC,
C    +   RINQ, ROUTQ, RPNTS, RTUBEI, RTUBEO,
C    +   S, S0, SHAPE, SPACE, TCAP, TCAPS, TDAT, TEAST, TEMMAX, TEX,
C    +   TEX0, TEX4, TFIX, THICK, THICKI, TIME, TMID, TRAN, TSIDE,
C    +   TSTEP, TTGAS, TWAF, TWAF0, TWAF4, TWEST, TX, TX0, TX4,
C    +   VINJ, WR, WZ,
C    +   WZWAF, Z, ZBEG, ZBEGI, ZDAT, ZEND, ZENDI, ZLEN, ZWAF,
C    +   GBOT, GTOP, JBOT, JTOP, NUMTC, SUMB, SUMT, TC, ZTC)

C     THE ARGUMENT LIST FOR TWAF1 IS THE SAME FOR ALL CALLS IN TWAFER,
C     BUT THE LIST IS DIFFERENT FOR THE CALLS IN TWAF5.  THE ARGUMENT
C     LIST FOR TWAF5 (EXCEPT FOR THE LAST LINE) IS THE SAME AS THE
C     ARGUMENT LIST FOR TWAF1.

      CHOICE = 0
      CALL TWAF1
     +  (ERROR, TEXT,
     +   R(QA2), IWORK(QA2DIM), R(QAEX), R(QAIN), R(QAREAC),
     +   R(QTAMB), LWORK(QCFLO), CHOICE,
     +   IWORK(QCOMPS), R(QCOND), R(QCONW), R(QCP),
     +   R(QDZEND), R(QDZNDI), R(QEMIS), R(QEMISS),
     +   R(QEMISE), R(QEMISI), R(QEMISW), R(QF),
     +   R(QFEE), R(QFEG), R(QFENDS), R(QFENG), R(QFESE), R(QFEW),
     +   R(QFGE), R(QFGSE), R(QFSEE), R(QFSEG), R(QFSES), R(QFSEW),
     +   R(QFSS), R(QFSW), R(QFWE), R(QFWS), R(QFWSE), R(QFWW),
     +   IWORK(QGROPA), IWORK(QGROPB), R(QHCOEF), IWORK(QIBEG),
     +   IWORK(QIEND), IWORK(QIPV2), IWORK(QIPV3), IWORK(QIPV4),
     +   R(QJ1), R(QJ2), R(QJEND), R(QJEX), R(QJIN),
     +   R(QJWALL), R(QLOS1), R(QLOS2), R(QLOS3), IWORK(QMAXHT),
     +   IWORK(QMAXIN), R(QMFLO), R(QMINJ), R(QMINJE), R(QMWZ),
     +   R(QMWZW),
     +   IWORK(QNUMEX), IWORK(QNUMHT), IWORK(QNUMIN), IWORK(QNUMT),
     +   IWORK(QNUMWF), R(QPCAP), IWORK(QPDIM), IWORK(QPEX),
     +   IWORK(QPIN), IWORK(QPMAX),
     +   R(QPOWER), R(QPOWRH), R(QPWALL), R(QQWAF), R(QR), RBOAT,
     +   R(QRCPV), R(QREFL), R(QRHO), R(QRING), LWORK(QRINJE),
     +   R(QRINQ), R(QROUTQ),
     +   IWORK(QRPNTS), R(QRTUBI), R(QRTUBO), R(QBUF), R(QX0), SHAPE,
     +   R(QSPACE), R(QTCAP), LWORK(QTCAPS), R(QTDAT), R(QTEAST),
     +   IWORK(QTEMAX), R(QTEX),
     +   R(QTEX0), R(QTEX4), TFIX, R(QTHICK), R(QTHCKI), TIME,
     +   R(QTMID), R(QTRAN), R(QTSIDE), R(QTSTEP), R(QTTGAS),
     +   R(QTWAF), R(QTWAF0),
     +   R(QTWAF4), R(QTWEST), R(QTX), R(QTX0), R(QTX4), R(QVINJ),
     +   R(QWR), R(QWZ),
     +   R(QWZWAF), R(QZ), R(QZBEG), R(QZBEGI), R(QZDAT), R(QZEND),
     +   R(QZENDI),  R(QZLEN),
     +   R(QZWAF), R(QGBOT), R(QGTOP), R(QJBOT), R(QJTOP),
     +   IWORK(QNUMTC), R(QSUMB), R(QSUMT), R(QTCT), R(QZTC))

      IF (ERROR) GO TO 9502

      DO 5025 J = 0, GROUPA + COMPS * POINTS + GROUPB - 1
         R(QBUF + J) = R(QF + J)
5025  CONTINUE

      GO TO 5030
      END IF

C///  SAVE THE SOLUTION FOR RESTARTING.

      ELSE IF (SIGNAL .EQ. 'SAVE' ) THEN

C///  SHOW THE SOLUTION.

      ELSE IF (SIGNAL .EQ. 'SHOW' ) THEN

C     SUBROUTINE TWAF2
C    +  (ERROR, TEXT,
C    +   CFLOW, COLMAX, COMPS, DATA4, GROUPA, GROUPB, HEADER,
C    +   JEND, JEX,
C    +   JIN, JWALL, KOUNT, LOS1, LOS2, LOS3, NADD,
C    +   NUMEX, NUMTC, NUMWF, PDIM, PEX, PIN, PMAX, R, RINQ, ROUTQ,
C    +   RPNTS, RTUBEI, RTUBEO, SUCCES, TC, TEX, THICKI, TIME, TIME2,
C    +   TSTEP, TWAF, TX, VAL2, VAL3, VALUE, WZWAF, X, Z, ZWAF)

      CALL TWAF2
     +  (ERROR, TEXT,
     +   CFLOW, COLMAX, COMPS, DATA4, GROUPA, GROUPB, CWORK(QHEAD),
     +   R(QJEND),
     +   R(QJEX), R(QJIN), R(QJWALL), KOUNT,
     +   R(QLOS1), R(QLOS2), R(QLOS3),
     +   NADD, NUMEX, NUMTC, NUMWF, PDIM, PEX,
     +   PIN, PMAX, R(QR), RINQ, ROUTQ, RPNTS, RTUBEI, RTUBEO, SUCCES,
     +   R(QTCT), R(QTEX), THICKI, TIME, TIME2, TSTEP, R(QTWAF),
     +   R(QTX), R(QVAL2), R(QVAL3), R(QVALUE), R(QWZWAF), R(QX),
     +   R(QZ), R(QZWAF))

      IF (ERROR) GO TO 9505

      IF (0 .LT. TEXT) WRITE (TEXT, 10003)

C///  SOLVE THE LINEAR EQUATIONS.

      ELSE IF (SIGNAL .EQ. 'SOLVE' ) THEN

C     SUBROUTIINE TWSOLV
C    +  (ERROR, TEXT,
C    +   A, ASIZE, BUFFER, COMPS, GROUPA, GROUPB, PIVOT, POINTS)

      CALL TWSOLV
     +  (ERROR, TEXT,
     +   R(QA), AMAX, R(QBUF), COMPS, GROUPA, GROUPB, IWORK(QPIVOT),
     +   POINTS)

      IF (ERROR) GO TO 9504

C///  STORE THE SOLUTION FOR TIME INTEGRATION.

      ELSE IF (SIGNAL .EQ. 'RETAIN' ) THEN

      DO 5040 J = 0, GROUPA + COMPS * POINTS  + GROUPB - 1
         R(QX0 + J) = R(QBUF + J)
5040  CONTINUE

C///  UPDATE INTERPOLATED VALUES.

      ELSE IF (SIGNAL .EQ. 'UPDATE' ) THEN

C///  BOTTOM OF THE BLOCK TO SERVICE REQUESTS FROM TWOPNT.

         END IF
         GO TO 5010
      END IF

C///  TWOPNT FINISHES.

      ERROR = .NOT. (SIGNAL .EQ. ' ' )
      IF (ERROR) GO TO 9506

C///////////////////////////////////////////////////////////////////////
C
C     (6) EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  WRITE THERMCOUPLE TEMPERATURES

C     COPY THE SOLUTION TO BUFFER SO THAT TWAF1'S ARGUMENTS ARE THE
C     SAME AS IN EARLIER CALLS
      DO 6010 J = 0, GROUPA + COMPS * POINTS + GROUPB - 1
         R(QBUF + J) = R(QX + J)
6010  CONTINUE

C      SUBROUTINE TWAF1
C     +  (ERROR, TEXT,
C     +   A2, A2DIM, AEX, AIN, AREAC, BTEMP, CFLOW, CHOICE,
C     +   COMPS, CONDUC,
C     +   CONW, CP, DZEND, DZENDI, EMIS, EMISS, EMISSE, EMISSI, EMISSW,
C     +   F, FEE, FEG, FENDS, FENG, FESE, FEW, FGE, FGSE, FSEE, FSEG,
C     +   FSES, FSEW, FSS,
C     +   FSW, FWE, FWS, FWSE, FWW, GROUPA, GROUPB, HCOEF, IBEG, IEND,
C     +   IPVT2, IPVT3, IPVT4, J1, J2, JEND, JEX, JIN, JWALL, LOS1,
C     +   LOS2, LOS3, MAXHT, MAXIN, MFLOW, MINJ, MINJE, MWZ, MWZW,
C     +   NUMEX, NUMHT, NUMIN, NUMT, NUMWF, PCAP, PDIM, PEX, PIN,
C     +   PMAX, POWER, POWERH, PWALL, QWAF, R, RBOAT, RCPV, REFL,
C     +   RHO, RING, RINJEC, RINQ, ROUTQ,
C     +   RPNTS, RTUBEI, RTUBEO,
C     +   S, S0, SHAPE, SPACE, TCAP, TCAPS, TDAT, TEAST, TEMMAX, TEX,
C     +   TEX0, TEX4, TFIX, THICK, THICKI, TIME, TMID, TRAN, TSIDE,
C     +   TSTEP, TTGAS, TWAF, TWAF0, TWAF4, TWEST, TX, TX0, TX4,
C     +   VINJ, WR, WZ, WZWAF,
C     +   Z, ZBEG, ZBEGI, ZDAT, ZEND, ZENDI, ZLEN, ZWAF,
C     +   GBOT, GTOP, JBOT, JTOP, NUMTC, SUMB, SUMT, TC, ZTC)

C     THE ARGUMENT LIST FOR TWAF1 IS THE SAME FOR ALL CALLS IN TWAFER,
C     BUT THE LIST IS DIFFERENT FOR THE CALLS IN TWAF5.  THE ARGUMENT
C     LIST FOR TWAF5 (EXCEPT FOR THE LAST LINE) IS THE SAME AS THE
C     ARGUMENT LIST FOR TWAF1.

      CHOICE = 0
      CALL TWAF1
     +  (ERROR, TEXT,
     +   R(QA2), A2DIM, R(QAEX), R(QAIN), R(QAREAC), TAMB,
     +   LWORK(QCFLO), CHOICE, COMPS, CONDUC, R(QCONW), CP,
     +   DZEND, DZENDI, EMIS, EMISS, EMISSE, EMISSI, EMISSW, R(QF),
     +   R(QFEE), R(QFEG), R(QFENDS), R(QFENG), R(QFESE), R(QFEW),
     +   R(QFGE), R(QFGSE), R(QFSEE), R(QFSEG), R(QFSES), R(QFSEW),
     +   R(QFSS), R(QFSW), R(QFWE), R(QFWS), R(QFWSE), R(QFWW),
     +   GROUPA, GROUPB, HCOEF, IWORK(QIBEG), IWORK(QIEND),
     +   IWORK(QIPV2), IWORK(QIPV3), IWORK(QIPV4), R(QJ1),
     +   R(QJ2), R(QJEND), R(QJEX), R(QJIN),
     +   R(QJWALL), R(QLOS1), R(QLOS2), R(QLOS3), MAXHT,
     +   IWORK(QMAXIN), R(QMFLO), R(QMINJ), R(QMINJE), R(QMWZ),
     +   R(QMWZW), NUMEX, NUMHT, IWORK(QNUMIN), NUMT, NUMWF,
     +   PCAP, PDIM, PEX, PIN, PMAX,
     +   R(QPOWER), R(QPOWRH), R(QPWALL), R(QQWAF), R(QR), RBOAT,
     +   R(QRCPV), REFL, RHO, R(QRING), LWORK(QRINJE),
     +   RINQ, ROUTQ, RPNTS, RTUBEI, RTUBEO, R(QBUF), R(QX0), SHAPE,
     +   SPACE, TCAP, TCAPS, R(QTDAT), R(QTEAST), TEMMAX, R(QTEX),
     +   R(QTEX0), R(QTEX4), TFIX, THICK, THICKI, TIME,
     +   R(QTMID), TRAN, R(QTSIDE), TSTEP, R(QTTGAS),
     +   R(QTWAF), R(QTWAF0), R(QTWAF4), R(QTWEST), R(QTX),
     +   R(QTX0), R(QTX4), R(QVINJ), R(QWR), R(QWZ), R(QWZWAF), R(QZ),
     +   R(QZBEG), R(QZBEGI), R(QZDAT), R(QZEND), R(QZENDI), ZLEN,
     +   R(QZWAF), R(QGBOT), R(QGTOP), R(QJBOT), R(QJTOP), NUMTC,
     +   R(QSUMB), R(QSUMT), R(QTCT), R(QZTC))
      IF (ERROR) GO TO 9502

      WRITE (TEXT, 10004)
     +   ID, (R(QZTC + I), R(QTCT + I)**0.25, I = 0, NUMTC - 1)

C///  WRITE GLOBAL ENERGY BALANCE

C      SUBROUTINE TWAF7
C     +  (ERROR, TEXT,
C     +   LOS1, LOS2, LOS3, MAXHT, NUMHT, NUMWF, PCAP, PDIM, PEX, PIN,
C     +   POWERH)

      IF(BFOUND(QHT)) THEN

      CALL TWAF7
     +  (ERROR, TEXT,
     +   R(QLOS1), R(QLOS2), R(QLOS3), MAXHT, NUMHT, NUMWF, PCAP, PDIM,
     +   PEX, PIN, R(QPOWRH))

      ELSE
      ENDIF

C///  WRITE THE PLOT2D DATA FILE.

      IF (0 .LT. DATA1) THEN

      IF (0 .LT. TEXT) WRITE (TEXT, 10002)
     +   ID, 'WRITING THE PLOT2D DATA FILE.'

C     SUBROUTINE TWAF3
C    +  (ERROR, TEXT,
C    +   COMPS, COORD, FLAG, GROUPA, GROUPB, NUMEX, NUMWF,
C    +   PIN, PMAX, R, REXT, RPNTS, TEMPER, TEX, TIME, TIME2,
C    +   TWAF, UNIT, WZWAF, X, ZFIRST, ZLAST, ZWAF)

      CALL TWAF3
     +  (ERROR, TEXT,
     +   COMPS, R(QCOOR), FLAG3, GROUPA, GROUPB, NUMEX, NUMWF,
     +   PIN, PMAX, R(QR), R(QREXT), RPNTS, R(QTEMPE), R(QTEX),
     +  .FALSE., TIME2, R(QTWAF), DATA1, R(QWZWAF), R(QX),
     +   ZFIRST, ZLAST, R(QZWAF))

      END IF

C///  WRITE THE TEMPERATURE DATA FILE.

      IF (0 .LT. DATA2) THEN

      WRITE (TEXT, 10002) ID, 'WRITING THE TEMPERATURE DATA FILE.'

C      SUBROUTINE TWAF4
C     +  (ERROR, TEXT,
C     +   NUMEX, NUMWF, PDIM, PEX, PIN, R, RPNTS,
C     +   TEX, TOUT, TWAF, TX, UNIT, WZWAF, Z, ZFIRST,
C     +   ZLAST, ZLEN, ZOUT, ZWAF)

      CALL TWAF4
     +  (ERROR, TEXT,
     +   NUMEX, NUMWF, PDIM, PEX, PIN, R(QR), RPNTS, R(QTEX),
     +   R(QTOUT2), R(QTWAF), R(QTX), DATA2, R(QWZWAF), R(QZ),
     +   ZFIRST, ZLAST, ZLEN, R(QZOUT), R(QZWAF))

      ENDIF

C*****TWAF5: SUBROUTINE
C
C      IF (MATRIX) THEN
C
C      IF (0 .LT. DATA2) WRITE (TEXT, 10002)
C     +   ID, 'WRITING THE COEFFICIENT MATRIX.'
C
CC     SUBROUTINE TWAF5
CC    +  (ERROR, TEXT,
CC    +   A2, A2DIM, AEX, AIN, BTEMP, CHOICE, COMPS, CONDUC, CP, DZEND,
CC    +   DZENDI, EMIS, EMISS, EMISSE, EMISSI, EMISSW, F, FEE, FEG,
CC    +   FENDS, FENG, FESE, FEW, FGE, FGSE, FSEE, FSEG, FSES, FSEW,
CC    +   FSS,
CC    +   FSW, FWE, FWS, FWSE, FWW, GROUPA, GROUPB, HCOEF, IBEG, IEND,
CC    +   IPVT2, IPVT3, IPVT4, J1, J2, JEND, JEX, JIN, JWALL, LOS1,
CC    +   LOS2, LOS3, MAXHT, NUMEX,
CC    +   NUMHT, NUMT, NUMWF, PCAP, PDIM, PEX, PIN, PMAX, POWER, POWERH,
CC    +   PWALL, QWAF, R, REFL, RHO, RINQ, ROUTQ, RPNTS, RTUBEI, RTUBEO,
CC    +   S, S0, SHAPE, SPACE, TCAP, TCAPS, TDAT, TEAST, TEMMAX, TEX,
CC    +   TEX0, TEX4, TFIX, THICK, THICKI, TIME, TMID, TRAN, TSIDE,
CC    +   TSTEP, TWAF, TWAF0, TWAF4, TWEST, TX, TX0, TX4, WR, WZ, WZWAF,
CC    +   Z, ZBEG, ZDAT, ZEND, ZLEN, ZWAF,
CC    +   GBOT, GTOP, JBOT, JTOP, NUMTC, SUMB, SUMT, TC, ZTC,
CC    +   DATA3, POINTS)
C
CC     THE ARGUMENT LIST FOR TWAF1 IS THE SAME FOR ALL CALLS IN TWAFER,
CC     BUT THE LIST IS DIFFERENT FOR THE CALLS IN TWAF5.  THE ARGUMENT
CC     LIST FOR TWAF5 (EXCEPT FOR THE LAST LINE) IS THE SAME AS THE
CC     ARGUMENT LIST FOR TWAF1.
C
C      CALL TWAF5
C     +  (ERROR, TEXT,
C     +   R(QA2), A2DIM, AEX, AIN, TAMB, CHOICE, COMPS, CONDUC, CP,
C     +   DZEND, DZENDI, EMIS, EMISS, EMISSE, EMISSI, EMISSW, R(QF),
C     +   R(QFEE), R(QFEG), R(QFENDS), R(QFENG), R(QFESE), R(QFEW),
C     +   R(QFGE), R(QFGSE), R(QFSEE), R(QFSEG), R(QFSES), R(QFSEW),
C     +   R(QFSS), R(QFSW), R(QFWE), R(QFWS), R(QFWSE), R(QFWW), GROUPA,
C     +   GROUPB, HCOEF, IWORK(QIBEG), IWORK(QIEND), IWORK(QIPV2),
C     +   IWORK(QIPV3), IWORK(QIPV4), R(QJ1), R(QJ2), R(QJEND),
C     +   R(QJEX), R(QJIN), R(QJWALL), R(QLOS1), R(QLOS2), R(QLOS3),
C     +   MAXHT,
C     +   NUMEX, NUMHT, NUMT, NUMWF, PCAP, PDIM, PEX, PIN, PMAX,
C     +   R(QPOWER), R(QPOWRH), R(QPWALL), R(QQWAF), R(QR), REFL, RHO,
C     +   RINQ, ROUTQ, RPNTS, RTUBEI, RTUBEO, R(QBUF), R(QX0), SHAPE,
C     +   SPACE, TCAP, TCAPS, R(QTDAT), R(QTEAST), TEMMAX, R(QTEX),
C     +   R(QTEX0), R(QTEX4), TFIX, THICK, THICKI, TIME,
C     +   R(QTMID), TRAN, R(QTSIDE), TSTEP, R(QTWAF), R(QTWAF0),
C     +   R(QTWAF4), R(QTWEST), R(QTX), R(QTX0), R(QTX4), R(QWR),
C     +   R(QWZ),
C     +   R(QWZWAF), R(QZ), R(QZBEG), R(QZDAT), R(QZEND), ZLEN,
C     +   R(QZWAF),
C     +   R(QGBOT), R(QGTOP), R(QJBOT), R(QJTOP), NUMTC, R(QSUMB),
C     +   R(QSUMT), R(QTCT), R(QZTC),
C     +   DATA3, POINTS)
C      END IF
C*****END TWAF5: SUBROUTINE

C///  WRITE MEMORY USE.

      IF (0 .LT. TEXT) WRITE (TEXT, 10005) ID,
     +   CMAX, IMAX, LMAX, RMAX,
     +   CSIZE - CMAX, ISIZE - IMAX, LSIZE - LMAX, RSIZE - RMAX,
     +   CSIZE, ISIZE, LSIZE, RSIZE

      GO TO 99999

77777 CONTINUE

C///////////////////////////////////////////////////////////////////////
C
C     (7) TRANSIENT SECTION OF THE CODE (USING DASSL).
C
C///////////////////////////////////////////////////////////////////////

C///  SPECIFY DASSL CONTROL DEFAULTS.

      IRES = 0

C     NUMBER OF EQUATIONS
C     UPPER AND LOWER BAND WIDTH
C     NUMBER OF CONSTRAINT EQUATIONS (NOT ACTIVE).

      NEQ = GROUPA + COMPS*NUMWF + GROUPB + NADD

      MU = COMPS + MAX(GROUPA, GROUPB, COMPS) - 1
      ML = MU
      NG = 0

C     SET LRW AND LIW THE LENGTH OF RWORK AND IWORK

C*****SET DASSL STORAGE AND FLAGS FOR FULLY DENSE MATRIX
C      LRW = 50 + (MAXORD+4)*NEQ + NEQ**2
C*****END SET DASSL STORAGE AND FLAGS FOR FULLY DENSE MATRIX

C*****SET DASSL STORAGE AND FLAGS FOR BANDED JACOBIAN
      LRW = 50 + (MAXORD + 4)*NEQ + (2*ML*MU+1)*NEQ
     +      + 2*(NEQ/(ML+MU+1)+1)
C*****END SET DASSL STORAGE AND FLAGS FOR BANDED JACOBIAN
      LIW = 20 + NEQ

C///  CREATE SOME STORAGE SPACE

C     IWOK(LIW)
C     INTEGER DASSL

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QIWOK', LIW, QIWOK)
      IF (ERROR) GO TO 9101

C     JROOT(NG)
C     INTEGER DASSL

      CALL RESERV (ERROR, TEXT, .FALSE., ILAST, IMAX, ISIZE,
     +   'QJROOT', NG+1, QJROOT)
      IF (ERROR) GO TO 9101

C     RWOK(LRW)
C     REAL DASSL

      CALL RESERV (ERROR, TEXT, .FALSE., RLAST, RMAX, RSIZE,
     +   'QRWOK', LRW, QRWOK)
      IF (ERROR) GO TO 9101

      ERROR = .NOT. (CLAST .LE. CSIZE .AND. ILAST .LE. ISIZE .AND.
     +   LLAST .LE. LSIZE .AND. RLAST .LE. RSIZE)
      IF (ERROR) GO TO 9401

C///  STEP-UP SOME PARAMETERS FOR THE INTEGRATION LOOP.

      DATA INFO /15*0/

C*****PRINTOUT AT EVERY STEP
C      INFO(3) = 1
C*****END PRINTOUT AT EVERY STEP
C*****PRINTOUT ONLY AT TOUT
      INFO(3) = 0
C*****END PRINTOUT ONLY AT TOUT

      INFO(6) = 1

      IWORK(QIWOK) = ML
      IWORK(QIWOK + 1) = MU
      IWORK(QIWOK + 2) = MAXORD

C*****SET DASSL STORAGE AND FLAGS FOR FULLY DENSE MATRIX
C      INFO(6) = 0
C*****END SET DASSL STORAGE AND FLAGS FOR FULLY DENSE MATRIX

      IF (MAXORD .NE. 5) INFO(9) = 1

      TIME1 = 0.0D16
      TOUT = 0.0D16

C///  SET ADD-ON VARIABLES IF REQUIRED.

      IF (NADD .NE. 0) THEN
          OFFSET = GROUPA + COMPS*NUMWF + GROUPB
          DO 610 I = 1, NADD
             R(QX + OFFSET + I - 1) = 0.0
             R(QX0 + OFFSET + I - 1) = 0.0
610       CONTINUE
      ELSE
      ENDIF

C///  FIND THE INITIAL YPRIME BY CALLING RESIDUAL
C     WITH YPRIME SET EQUAL TO ZERO.

      DO 700 I = 1, NEQ
         R(QBUF+I-1) = 0.0
700   CONTINUE

C      SUBROUTINE RES
C     +  (T, S, SP, F, IRES, RPAR, IPAR)

      CALL RES(TIME1, R(QX), R(QBUF), R(QF), IRES, R, IWORK)

      DO 710 I = 1, NEQ
         R(QX0+I-1) = -R(QF+I-1)
C
C        CHANGE DUE TO ERROR W.G.H. 4/27/94
C
C         R(QX0+I-1) = -R(QBUF+I-1)
710   CONTINUE

C///  PRINT THE INITIAL SOLUTION.

C      SUBROUTINE TWAF2
C     +  (ERROR, TEXT,
C     +   CFLOW, COLMAX, COMPS, DATA4, GROUPA, GROUPB, HEADER, JEND,
C     +   JEX, JIN, JWALL, KOUNT, LOS1, LOS2, LOS3, NADD,
C     +   NUMEX, NUMTC, NUMWF, PDIM, PEX, PIN, PMAX, R, RINQ, ROUTQ,
C     +   RPNTS, RTUBEI, RTUBEO, SUCCES, TC, TEX, THICKI, TIME, TIME2,
C     +   TSTEP, TWAF, TX, VAL2, VAL3, VALUE, WZWAF, X, Z, ZWAF)

      CALL TWAF2
     +  (ERROR, TEXT,
     +   CFLOW, COLMAX, COMPS, DATA4, GROUPA, GROUPB, CWORK(QHEAD),
     +   R(QJEND), R(QJEX), R(QJIN), R(QJWALL), KOUNT,
     +   R(QLOS1), R(QLOS2), R(QLOS3),
     +   NADD, NUMEX, NUMTC, NUMWF, PDIM, PEX,
     +   PIN, PMAX, R(QR), RINQ, ROUTQ, RPNTS, RTUBEI, RTUBEO, SUCCES,
     +   R(QTCT), R(QTEX), THICKI, TIME, TIME1, TSTEP, R(QTWAF), R(QTX),
     +   R(QVAL2), R(QVAL3), R(QVALUE), R(QWZWAF), R(QX), R(QZ),
     +   R(QZWAF))

C///  PRINT THE INITIAL SOLUTION ON THE PLOT FILE.

C     SUBROUTINE TWAF3
C    +  (ERROR, TEXT,
C    +   COMPS, COORD, FLAG, GROUPA, GROUPB, NUMEX, NUMWF, PIN,
C    +   PMAX, R, REXT, RPNTS, TEMPER, TEX, TIME, TIME2, TWAF, UNIT,
C    +   WZWAF, X, ZFIRST, ZLAST, ZWAF)

C*****PLOT FILE PRINTING ON
C      CALL TWAF3
C     +  (ERROR, TEXT,
C     +   COMPS, R(QCOOR), FLAG3, GROUPA, GROUPB, NUMEX, NUMWF, PIN,
C     +   PMAX, R(QR), R(QREXT), RPNTS, R(QTEMPE), R(QTEX), .TRUE.,
C     +   TIME1, R(QTWAF), DATA1, R(QWZWAF), R(QX), ZFIRST, ZLAST,
C     +   R(QZWAF))
C*****END PLOT FILE PRINTING ON

C///  TOP OF THE INTEGRATION LOOP

790   CONTINUE

      TOUT = TOUT + DTIME

C      SUBROUTINE DDASSL (RES,NEQ,T,Y,YPRIME,TOUT,
C     *  INFO,RTOL,ATOL,IDID,
C     *  RWORK,LRW,IWORK,LIW,RPAR,IPAR,
C     *  JAC)

C      SUBROUTINE DDASRT (RES,NEQ,T,Y,YPRIME,TOUT,
C     *  INFO,RTOL,ATOL,IDID,RWORK,LRW,IWORK,LIW,RPAR,IPAR,JAC,
C     *  G,NG,JROOT)

720   CONTINUE

C*****DOUBLE PRECISION DASSL
      CALL DDASSL(RES, NEQ, TIME1, R(QX), R(QX0), TOUT,
     +            INFO, RTOL, ATOL, IDID, R(QRWOK), LRW,
     +            IWORK(QIWOK), LIW, R, IWORK, JAC)
C*****END DOUBLE PRECISION DASSL

C*****DOUBLE PRECISION DASSRT
C      CALL DDASRT(RES, NEQ, TIME1, R(QX), R(QX0), TOUT,
C     +            INFO, RTOL, ATOL, IDID, R(QRWOK), LRW,
C     +            IWORK(QIWOK), LIW, R, IWORK, JAC,
C     +            ROOT, NG, IWOK(QJROOT))
C*****END DOUBLE PRECISION DASSRT

C*****SINGLE PRECISION
C      CALL SDASSL(RES, NEQ, TIME1, R(QX), R(QX0), TOUT,
C     +            INFO, RTOL, ATOL, IDID, R(QRWOK), LRW,
C     +            IWORK(QIWOK), LIW, RPAR, IPAR, JAC)
C*****END SINGLE PRECISION

      IF (IDID .EQ. -1) THEN
          WRITE(TEXT, 10006) IDID
          INFO (1) = 1
          GO TO 720
      ENDIF

      ERROR = (IDID .LT. -1)
      IF (ERROR) GO TO 9007

      IF (IDID .EQ. 1) THEN
          WRITE(TEXT, 10007) IDID, TIME1, TOUT
          GO TO 720
      ELSEIF ((IDID.EQ.2) .OR. (IDID.EQ.3) .OR. (IDID.EQ.4)) THEN

C///  PRINT RESULTS.

C      SUBROUTINE TWAF2
C     +  (ERROR, TEXT,
C     +   CFLOW, COLMAX, COMPS, DATA4, GROUPA, GROUPB, HEADER, JEND,
C     +   JEX, JIN, JWALL, KOUNT, LOS1, LOS2, LOS3, NADD,
C     +   NUMEX, NUMTC, NUMWF, PDIM, PEX, PIN, PMAX, R, RINQ, ROUTQ,
C     +   RPNTS, RTUBEI, RTUBEO, SUCCES, TC, TEX, THICKI, TIME, TIME2,
C     +   TSTEP, TWAF, TX, VAL2, VAL3, VALUE, WZWAF, X, Z, ZWAF)

      CALL TWAF2
     +  (ERROR, TEXT,
     +   CFLOW, COLMAX, COMPS, DATA4, GROUPA, GROUPB, CWORK(QHEAD),
     +   R(QJEND), R(QJEX), R(QJIN), R(QJWALL), KOUNT,
     +   R(QLOS1), R(QLOS2), R(QLOS3),
     +   NADD, NUMEX, NUMTC, NUMWF, PDIM, PEX, PIN, PMAX, R(QR), RINQ,
     +   ROUTQ, RPNTS, RTUBEI, RTUBEO, SUCCES, R(QTCT), R(QTEX),
     +   THICKI, TIME, TOUT, TSTEP, R(QTWAF), R(QTX), R(QVAL2),
     +   R(QVAL3), R(QVALUE), R(QWZWAF), R(QX), R(QZ), R(QZWAF))

C///  PLOT RESULTS.

C     SUBROUTINE TWAF3
C    +  (ERROR, TEXT,
C    +   COMPS, COORD, FLAG, GROUPA, GROUPB, NUMEX, NUMWF, PIN,
C    +   PMAX, R, REXT, RPNTS, TEMPER, TEX, TIME, TIME2, TWAF,
C    +   UNIT, WZWAF, X, ZFIRST, ZLAST, ZWAF)

      CALL TWAF3
     +  (ERROR, TEXT,
     +   COMPS, R(QCOOR), FLAG3, GROUPA, GROUPB, NUMEX, NUMWF,
     +   PIN, PMAX, R(QR), R(QREXT), RPNTS, R(QTEMPE), R(QTEX),
     +  .TRUE., TOUT, R(QTWAF), DATA1, R(QWZWAF), R(QX), ZFIRST,
     +   ZLAST, R(QZWAF))

         IF (IDID .EQ. 4) THEN
            INFO(1) = 1
            GO TO 720
         ELSE
         ENDIF

      ELSE
      ENDIF

      IF (TIME1 .LT. TSTOP) GO TO 790

C///  WRITE THE TEMPERATURE DATA FILE AT THE END OF THE TRANSIENT.

      IF (0 .LT. DATA2) THEN

      WRITE (TEXT, 10002) ID, 'WRITING THE TEMPERATURE DATA FILE.'

C      SUBROUTINE TWAF4
C     +  (ERROR, TEXT,
C     +   NUMEX, NUMWF, PDIM, PEX, PIN, R, RPNTS,
C     +   TEX, TOUT, TWAF, TX, UNIT, WZWAF, Z, ZFIRST,
C     +   ZLAST, ZLEN, ZOUT, ZWAF)

      CALL TWAF4
     +  (ERROR, TEXT,
     +   NUMEX, NUMWF, PDIM, PEX, PIN, R(QR), RPNTS, R(QTEX),
     +   R(QTOUT2), R(QTWAF), R(QTX), DATA2, R(QWZWAF), R(QZ),
     +   ZFIRST, ZLAST, ZLEN, R(QZOUT), R(QZWAF))

      ENDIF

C///  WRITE MEMORY USE.

      IF (0 .LT. TEXT) WRITE (TEXT, 10005) ID,
     +   CMAX, IMAX, LMAX, RMAX,
     +   CSIZE - CMAX, ISIZE - IMAX, LSIZE - LMAX, RSIZE - RMAX,
     +   CSIZE, ISIZE, LSIZE, RSIZE

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      DATA (BANNER(1, K), K=1, 10)
     +   /
     +   'OOOOOOOOOOO OO       OO   OOOOOOO   OOOOOOOOO  OOOOOOOOO',
     +   'OOOOOOOOOOO OO       OO  OOOOOOOOO  OOOOOOOOO  OOOOOOOOO',
     +   '    OOO     OO       OO  OO     OO  OO         OO       ',
     +   '    OOO     OO       OO  OO     OO  OO         OO       ',
     +   '    OOO     OO   O   OO  OOOOOOOOO  OOOOOOO    OOOOOOO  ',
     +   '    OOO     OO  OOO  OO  OOOOOOOOO  OOOOOOO    OOOOOOO  ',
     +   '    OOO     OO OO OO OO  OO     OO  OO         OO       ',
     +   '    OOO     OOOO   OOOO  OO     OO  OO         OO       ',
     +   '    OOO     OOO     OOO  OO     OO  OO         OOOOOOOOO',
     +   '    OOO     OO       OO  OO     OO  OO         OOOOOOOOO'
C         123456789_123456789_123456789_123456789_123456789_123456789_
     +   /
      DATA (BANNER(2, K), K = 1, 10)
     +   /
     +   'OOOOOOOO',
     +   'OOOOOOOOO',
     +   'OO     OO',
     +   'OO     OO',
     +   'OOOOOOOOO',
     +   'OOOOOOOO',
     +   'OO   OO',
     +   'OO    OO',
     +   'OO     OO',
     +   'OO     OO'
     +   /
C         123456789_

10001 FORMAT
     + (/10X, 35('/ ')
     + //1X, A9, A
     +///6X, A56, 2X, A10,
     +  /6X, A56, 2X, A10,
     +  /6X, A56, 2X, A10,
     +  /6X, A56, 2X, A10,
     +  /6X, A56, 2X, A10,
     +  /6X, A56, 2X, A10,
     +  /6X, A56, 2X, A10,
     +  /6X, A56, 2X, A10,
     +  /6X, A56, 2X, A10,
     +  /6X, A56, 2X, A10,
     +///10X, '         WAFER TEMPERATURE DETERMINATION IN A HIGH'
     +  /10X, '            TEMPERATURE MULTIPLE-WAFERS-IN-TUBE'
     +  /10X, '                      LPCVD REACTOR'
     + //10X, '          WILLIAM G. HOUF          JOSEPH F. GRCAR'
     +  /10X, '          DIVISION 8345            DIVISION 8345'
     +  /10X, '          (510) 294-3184           (510) 294-2662'
     + //10X,'SANDIA NATIONAL LABORATORY, LIVERMORE, CA 94551-0969 USA')

10002 FORMAT
     +  (/10X, 35('/ ')
     +   //1X, A9, A)

10003 FORMAT
     +  (/10X, 35('/ '))

10004 FORMAT
     +  (/10X, 35('/ ')
     +   //1X, A9, 'THERMOCOUPLE TEMPERATURES'
     +  //15X, 'Z   TEMPERATURE'
     + //(10X, F7.1, 2X, F7.1))

10005 FORMAT
     +  (/10X, 35('/ ')
     +   //1X, A9, 'THE WORK SPACE REQUIREMENTS ARE AS FOLLOWS.'
C               123456789  123456789  123456789  123456789  123456789
     +  //10X, '           CHARACTER    INTEGER    LOGICAL       REAL'
     +  //10X, '     USED', 4(2X, I9)
     +   /10X, '   EXCESS', 4(2X, I9)
     +   /10X, '    TOTAL', 4(2X, I9))

60001 FORMAT (/1X, A)

60002 FORMAT (/1X, I5)

60003 FORMAT (1X, I5, 2X, A)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID, CSIZE
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID, ISIZE
      GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID, LSIZE
      GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID, RSIZE
      GO TO 99999

9007  IF (0 .LT. TEXT) WRITE (TEXT, 99007) ID, IDID,
     +   TIME1, R(QRWOK+6), IWORK(QIWOK+7),
     +   IWORK(QIWOK+10), IWORK(QIWOK+11), IWORK(QIWOK+12),
     +   IWORK(QIWOK+13), IWORK(QIWOK+14)
      GO TO 99999

9101  IF (0 .LT. TEXT) WRITE (TEXT, 99101) ID
      GO TO 99999

9102  IF (0 .LT. TEXT) WRITE (TEXT, 99102) ID,
     +   CSIZE, ISIZE, LSIZE, RSIZE, CLAST, ILAST, LLAST, RLAST
      GO TO 99999

9103  IF (0 .LT. TEXT) WRITE (TEXT, 99103) ID
      GO TO 99999

9104  IF (0 .LT. TEXT) WRITE (TEXT, 99104) ID
      GO TO 99999

9301  IF (0 .LT. TEXT) WRITE (TEXT, 99301) ID
      GO TO 99999

9302  IF (0 .LT. TEXT) WRITE (TEXT, 99302) ID
      GO TO 99999

9303  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, WORD)
         WRITE (TEXT, 99303) ID, WORD (1 : LENGTH)
         DO 9010 J = 1, VMSSS
            STRING = VMSS(J)
            CALL EXTENT (LENGTH, STRING)
            WRITE (TEXT, '(10X, A, A)')
     +         ' ACCEPTED:  ', STRING (1 : LENGTH)
9010     CONTINUE
      END IF
      GO TO 99999

9304  IF (0 .LT. TEXT) THEN
         WORD = BNAME(Q)
         CALL EXTENT (LENGTH, WORD)
         IF (58 .LT. LENGTH) WORD (54 : ) = ' ...'
         CALL EXTENT (LENGTH, WORD)
         WRITE (TEXT, 99304) ID, NUMBER, WORD (1 : LENGTH)
      END IF
      GO TO 99999

9305  IF (0 .LT. TEXT) WRITE (TEXT, 99305) ID
      GO TO 99999
9306  IF (0 .LT. TEXT) WRITE (TEXT, 99306) ID
      GO TO 99999
9307  IF (0 .LT. TEXT) WRITE (TEXT, 99305) ID
      GO TO 99999

9312  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, WORD)
         WRITE (TEXT, 99312) ID, NUMBER, WORD (1 : LENGTH)
      END IF
      GO TO 99999

9313  IF (0 .LT. TEXT) THEN
         WRITE (TEXT, 99313) ID
         DO 9030 J = 1, BLOCKS
            IF (BNEED(J) .AND. .NOT. BFOUND(J)) THEN
               WORD = BNAME(J)
               CALL EXTENT (LENGTH, WORD)
               IF (58 .LT. LENGTH) WORD (54 : ) = ' ...'
               CALL EXTENT (LENGTH, WORD)
               WRITE (TEXT, '(10X, A, A)')
     +            '  KEYWORD:  ', WORD (1 : LENGTH)
            END IF
9030     CONTINUE
      END IF
      GO TO 99999

9314  IF (0 .LT. TEXT) WRITE (TEXT, 99314) ID
      GO TO 99999

9315  IF (0 .LT. TEXT) WRITE (TEXT, 99315) ID
      GO TO 99999

9316  IF (0 .LT. TEXT) WRITE (TEXT, 99316) ID
      GO TO 99999

9317  IF (0 .LT. TEXT) WRITE (TEXT, 99317) ID, RTUBEO, RTUBEI
      GO TO 99999

9318  IF (0 .LT. TEXT) WRITE (TEXT, 99318) ID, ROUTQ(1), RINQ(1)
      GO TO 99999

9319  IF (0 .LT. TEXT) WRITE (TEXT, 99319) ID, ROUTQ(2), RINQ(2)
      GO TO 99999

9320  IF (0 .LT. TEXT) WRITE (TEXT, 99320) ID, RINQ(2), ROUTQ(2),
     +                 RTUBEI, RTUBEO
      GO TO 99999

9321  IF (0 .LT. TEXT) WRITE (TEXT, 99321) ID, RINQ(1), ROUTQ(1),
     +                 RINQ(2), ROUTQ(2)
      GO TO 99999

9322  IF (0 .LT. TEXT) WRITE (TEXT, 99321) ID, RADIUS, RINQ(1)
      GO TO 99999

9323  IF (0 .LT. TEXT) WRITE (TEXT, 99323) ID, EMIS(1), REFL(1), TRAN(1)
      GO TO 99999

9324  IF (0 .LT. TEXT) WRITE (TEXT, 99324) ID, EMIS(2), REFL(2), TRAN(2)
      GO TO 99999

9325  IF (0 .LT. TEXT) WRITE (TEXT, 99325) ID, BNAME(QSS), BNAME(QTRANS)
      GO TO 99999

9326  IF (0 .LT. TEXT) WRITE (TEXT, 99326) ID
      GO TO 99999

9401  IF (0 .LT. TEXT) WRITE (TEXT, 99401) ID,
     +   CSIZE, ISIZE, LSIZE, RSIZE, CLAST, ILAST, LLAST, RLAST
      GO TO 99999

9406  IF (0 .LT. TEXT) WRITE (TEXT, 99406) ID
      GO TO 99999

9501  IF (0 .LT. TEXT) WRITE (TEXT, 99501) ID
      GO TO 99999

9502  IF (0 .LT. TEXT) WRITE (TEXT, 99502) ID
      GO TO 99999

9503  IF (0 .LT. TEXT) WRITE (TEXT, 99503) ID
      GO TO 99999

9504  IF (0 .LT. TEXT) WRITE (TEXT, 99504) ID
      GO TO 99999

9505  IF (0 .LT. TEXT) WRITE (TEXT, 99505) ID
      GO TO 99999

9506  IF (0 .LT. TEXT) THEN
         CALL EXTENT (LENGTH, REPORT)
         WRITE (TEXT, 99506) ID, REPORT(1 : LENGTH)
      END IF
      GO TO 99999

9507  IF (0 .LT. TEXT) WRITE (TEXT, 99507) ID
      GO TO 99999

9508  IF (0 .LT. TEXT) WRITE (TEXT, 99508) ID
      GO TO 99999

9601  IF (0 .LT. TEXT) WRITE (TEXT, 99601) ID
      GO TO 99999

9602  IF (0 .LT. TEXT) WRITE (TEXT, 99602) ID
      GO TO 99999

9603  IF (0 .LT. TEXT) WRITE (TEXT, 99603) ID
      GO TO 99999

10006 FORMAT
     +  (/10X, 'A LOT OF WORK EXPENDED IN DASSL, IDID = ', I6)

10007 FORMAT
     +  (/10X, 'STEP TAKEN IN INTERMEDIATE OUTPUT MODE, IDID = ', I6,
     +   /10X, 'THE CODE HAS NOT YET REACHED TOUT',
     +   /10X, 'TIME1 = ',E14.7,1X,'TOUT = ',E14.7)

99001 FORMAT
     +   (/1X, A9, 'ERROR.  THE SIZE OF THE CHARACTER WORK SPACE MUST'
     +   /10X, 'BE POSITIVE.'
     +   //10X, I10, '  SIZE')

99002 FORMAT
     +   (/1X, A9, 'ERROR.  THE SIZE OF THE INTEGER WORK SPACE MUST BE'
     +   /10X, 'POSITIVE.'
     +   //10X, I10, '  SIZE')

99003 FORMAT
     +   (/1X, A9, 'ERROR.  THE SIZE OF THE LOGICAL WORK SPACE MUST BE'
     +   /10X, 'POSITIVE.'
     +   //10X, I10, '  SIZE')

99004 FORMAT
     +   (/1X, A9, 'ERROR.  THE SIZE OF THE REAL WORK SPACE MUST BE'
     +   /10X, 'POSITIVE.'
     +   //10X, I10, '  SIZE')

99007 FORMAT
     +   (/1X, A9, 'ERROR.  TROUBLE WITH DASSL -- IDID = ', I6,
     +    /10X, 'CURRENT STATISTICS FROM DASSL.'
     +    /10X, 'TIME = ', E14.7,
     +    /10X, 'LAST STEP SIZE = ', E14.7,
     +    /10X, 'ORDER OF THE METHOD = ', I6,
     +    /10X, 'NUMBER OF STEPS TAKEN = ', I6,
     +    /10X, 'NUMBER OF CALL TO RES = ', I6,
     +    /10X, 'NUMBER OF JACOBIAN EVALUATIONS = ', I6,
     +    /10X, 'NUMBER OF ERROR TEST FAILURES = ', I6,
     +    /10X, 'NUMBER OF CONVERGENCE TEST FAILURES = ', I6)

99101 FORMAT
     +   (/1X, A9, 'ERROR.  RESERV FAILS.')

99102 FORMAT
     +   (/1X, A9, 'ERROR.  ONE OR MORE WORKSPACES ARE TOO SMALL TO'
     +   /10X, 'INITIALIZE THE KINETICS AND TRANPORT LIBRARIES.'
C               1234567890  1234567890  1234567890  1234567890
     +  //25X, ' CHARACTER     INTEGER     LOGICAL        REAL'
C               1234567890123
     +  //10X, ' PRESENT SIZE', 4I12
     +   /10X, 'REQUIRED SIZE', 4I12
     +  //10X, 'MORE SPACE MAY BE NEEDED LATER.')

99103 FORMAT
     +   (/1X, A9, 'ERROR.  SKSYMS FAILS.')

99104 FORMAT
     +   (/1X, A9, 'ERROR.  NEITHER TRANSIENT OR STEADY-STATE'
     +   /10X, 'SOLUTION WAS SPECIFIED.')

99301 FORMAT
     +   (/1X, A9, 'ERROR.  VERIFY FAILS.')

99302 FORMAT
     +   (/1X, A9, 'ERROR.  READW FAILS.')

99303 FORMAT
     +   (/1X, A9, 'ERROR.  THE INPUT SCRIPT IS INTENDED FOR A VERSION'
     +   /10X, 'OF TWAFER THAT MAY BE INCOMPATIBLE WITH THE PRESENT'
     +   /10X, 'VERSION.  CHANGE THE VERSION NUMBER IN THE SCRIPT FILE'
     +   /10X, 'TO ONE THAT IS ACCEPTED.'
     +  //10X, ' INTENDED:  ', A
     +   /)

99304 FORMAT
     +   (/1X, A9, 'ERROR.  THE INPUT SCRIPT REPEATS A KEYWORD.'
     +  //10X, I10, '  LINE NUMBER'
     +  //10X, '  KEYWORD:  ', A)

99305 FORMAT
     +   (/1X, A9, 'ERROR.  READV FAILS.')

99306 FORMAT
     +   (/1X, A9, 'ERROR.  BOTH HEATER POWERS AND HEATER'
     +    /10X, 'TEMPERATURES CANNOT BE SPECIFIED.')

99307 FORMAT
     +   (/1X, A9, 'ERROR.  BOTH END CAP POWERS AND END CAP'
     +    /10X, 'TEMPERATURES CANNOT BE SPECIFIED.')

99312 FORMAT
     +   (/1X, A9, 'ERROR.  A KEYWORD IS MISPLACED OR UNKNOWN.'
     +   //10X, I10, '  LINE NUMBER'
     +   //10X, '  KEYWORD:  ', A)

99313 FORMAT
     +   (/1X, A9, 'ERROR.  ONE OR MORE KEYWORDS ARE MISSING.'
     +   /)

99314 FORMAT
     +   (/1X, A9, 'ERROR. ONE OR MORE HEATER PARAMETERS ARE'
     +   /10X, 'NOT SPECIFIED OR ARE SPECIFIED INCORRECTLY'
     +  //10X, 'REQUIRED PARAMATERS FOR EACH HEATER  - I :',
     +   /10X, 'ZBEG(I) = ',
     +   /10X, 'ZEND(I) = ',
     +   /10X, 'POWER(I) = ')

99315 FORMAT
     +   (/1X, A9,'ERROR. THE KEYWORDS'
     +   /10X, 'TBEG = '
     +   /10X, 'TEND = '
     +   /10X, 'ARE REQUIRED')

99316 FORMAT
     +   (/1X, A9, 'ERROR. ONE OR MORE TEMPERATURE PARAMETERS ARE'
     +   /10X, 'NOT SPECIFIED OR ARE SPECIFIED INCORRECTLY'
     +  //10X, 'REQUIRED PARAMATERS FOR EACH TEMPERATURE - I :',
     +   /10X, 'Z(I) = ',
     +   /10X, 'T(I) = ')

99317 FORMAT
     +   (/1X, A9, 'ERROR. THE INNER RADIUS OF THE HEATERS IS',
     +   /10X, 'LARGER THAN THE OUTER RADIUS.'
     +   /10X, 'OUTER RADIUS = ',E14.7,
     +   /10X, 'INNER RADIUS = ',E14.7)

99318 FORMAT
     +   (/1X, A9, 'ERROR. THE INNER RADIUS OF THE FIRST QUARTZ',
     +   /10X, 'WALL IS LARGER THAN THE OUTER RADIUS.'
     +   /10X, 'OUTER RADIUS = ',E14.7,
     +   /10X, 'INNER RADIUS = ',E14.7)

99319 FORMAT
     +   (/1X, A9, 'ERROR. THE INNER RADIUS OF THE SECOND QUARTZ',
     +   /10X, 'WALL IS LARGER THAN THE OUTER RADIUS.'
     +   /10X, 'OUTER RADIUS = ',E14.7,
     +   /10X, 'INNER RADIUS = ',E14.7)

99320 FORMAT
     +   (/1X, A9, 'ERROR. THE RADII OF THE SECOND QUARTZ',
     +   /10X, 'WALL ARE GREATER THAN THE RADII OF THE HEATERS.'
     +   /10X, 'QUARTZ WALL INNER RADIUS = ', E14.7,
     +   /10X, 'QUARTZ WALL OUTER RADIUS = ', E14.7,
     +   /10X, 'HEATER INNER RADIUS = ',E14.7,
     +   /10X, 'HEATER OUTER RADIUS = ',E14.7)

99321 FORMAT
     +   (/1X, A9, 'ERROR. THE RADII OF THE FIRST QUARTZ WALL ARE',
     +   /10X, 'GREATER THAN THE RADII OF THE SECOND QUARTZ WALL.'
     +   /10X, 'QUARTZ WALL 1 INNER RADIUS = ', E14.7,
     +   /10X, 'QUARTZ WALL 1 OUTER RADIUS = ', E14.7,
     +   /10X, 'QUARTZ WALL 2 INNER RADIUS = ', E14.7,
     +   /10X, 'QUARTZ WALL 2 OUTER RADIUS = ', E14.7)

99322 FORMAT
     +   (/1X, A9, 'ERROR. THE RADII OF THE WAFERS ARE GREATER',
     +   /10X, 'THAN THE INNER RADII OF THE FIRST QUARTZ WALL'
     +   /10X, 'WAFER RADIUS = ',E14.7,
     +   /10X, 'QUARTZ WALL 1 INNER RADIUS = ',E14.7)

99323 FORMAT
     +   (/1X, A9, 'ERROR. THE SUM OF THE ABSORPTANCE, REFLECTANCE,',
     +   /10X, 'AND TRANSMITTANCE DOES NOT EQUAL ONE FOR THE FIRST'
     +   /10X, 'QUARTZ WALL.',
     +   /10X, 'ABSORPTANCE = ', E14.7,
     +   /10X, 'REFLECTANCE = ', E14.7,
     +   /10X, 'TRANSMITTANCE = ', E14.7)

99324 FORMAT
     +   (/1X, A9, 'ERROR. THE SUM OF THE ABSORPTANCE, REFLECTANCE,',
     +   /10X, 'AND TRANSMITTANCE DOES NOT EQUAL ONE FOR THE SECOND'
     +   /10X, 'QUARTZ WALL.',
     +   /10X, 'ABSORPTANCE = ',  E14.7,
     +   /10X, 'REFLECTANCE = ',  E14.7,
     +   /10X, 'TRANSMITTANCE = ', E14.7)

99325 FORMAT
     +   (/1X, A9, 'ERROR.  KEYWORDS FOR BOTH THE STEADY-STATE AND'
     +   /10X, 'TRANSIENT SOLUTIONS HAVE BEEN CHOSEN:  ONLY'
     +   /10X, 'ONE KEYWORD OR THE OTHER IS ALLOWED.'
     +   //10X,'STEADY-STATE KEYWORD:  ', A,
     +    /10X,'TRANSIENT KEYWORD:  ', A)

99326 FORMAT
     +   (/1X, A9, 'ERROR.  KEYWORDS FOR BOTH THE CAP POWERS AND'
     +   /10X, 'CAP TEMPERATURES HAVE BEEN SPECIFIED:  ONLY'
     +   /10X, 'ONE KEYWORD OR THE OTHER IS ALLOWED.')

99401 FORMAT
     +   (/1X, A9, 'ERROR.  ONE OR MORE WORKSPACES ARE TOO SMALL.'
C               1234567890  1234567890  1234567890  1234567890
     +  //25X, ' CHARACTER     INTEGER     LOGICAL        REAL'
C               123456789012345
     +  //10X, ' PRESENT SIZE', 4I12
     +   /10X, 'REQUIRED SIZE', 4I12
     +  //10X, 'NO MORE SPACE WILL BE NEEDED LATER.')

99406 FORMAT
     +   (/1X, A9, 'ERROR.  TWAF0 FAILS.')

99501 FORMAT
     +   (/1X, A9, 'ERROR.  TWOPNT FAILS.')

99502 FORMAT
     +   (/1X, A9, 'ERROR.  TWAF1 FAILS.')

99503 FORMAT
     +   (/1X, A9, 'ERROR.  TWPREP FAILS.')

99504 FORMAT
     +   (/1X, A9, 'ERROR.  TWSOLV FAILS.')

99505 FORMAT
     +   (/1X, A9, 'ERROR.  TWAF2 FAILS.')

99506 FORMAT
     +   (/1X, A9, 'ERROR.  TWOPNT FAILS.'
     +   //10X, I10, '  REPORT')

99507 FORMAT
     +   (/1X, A9, 'ERROR.  TWINIT FAILS.')

99508 FORMAT
     +   (/1X, A9, 'ERROR.  TWSETI, TWSETL OR TWSETR FAILS.')

99601 FORMAT
     +   (/1X, A9, 'ERROR.  TWAF5 FAILS.')

99602 FORMAT
     +   (/1X, A9, 'ERROR.  ERROR IN SPECIFICATION OF HEATERS.')

99603 FORMAT
     +   (/1X, A9, 'ERROR.  ERROR IN SPECIFICATION OF END CAP POWERS'
     +    /10X, 'OR TEMPERATURES.')

C///  EXIT.

99999 CONTINUE
      END
*DECK ERF
      FUNCTION ERF (X)
C***BEGIN PROLOGUE  ERF
C***PURPOSE  Compute the error function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C8A, L5A1E
C***TYPE      SINGLE PRECISION (ERF-S, DERF-D)
C***KEYWORDS  ERF, ERROR FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C ERF(X) calculates the single precision error function for
C single precision argument X.
C
C Series for ERF        on the interval  0.          to  1.00000D+00
C                                        with weighted error   7.10E-18
C                                         log weighted error  17.15
C                               significant figures required  16.31
C                                    decimal places required  17.71
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, ERFC, INITS, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900727  Added EXTERNAL statement.  (WRB)
C   920618  Removed space from variable name.  (RWC, WRB)
C***END PROLOGUE  ERF
      DIMENSION ERFCS(13)
      LOGICAL FIRST
      EXTERNAL ERFC
      SAVE ERFCS, SQRTPI, NTERF, XBIG, SQEPS, FIRST
      DATA ERFCS( 1) /   -.0490461212 34691808E0 /
      DATA ERFCS( 2) /   -.1422612051 0371364E0 /
      DATA ERFCS( 3) /    .0100355821 87599796E0 /
      DATA ERFCS( 4) /   -.0005768764 69976748E0 /
      DATA ERFCS( 5) /    .0000274199 31252196E0 /
      DATA ERFCS( 6) /   -.0000011043 17550734E0 /
      DATA ERFCS( 7) /    .0000000384 88755420E0 /
      DATA ERFCS( 8) /   -.0000000011 80858253E0 /
      DATA ERFCS( 9) /    .0000000000 32334215E0 /
      DATA ERFCS(10) /   -.0000000000 00799101E0 /
      DATA ERFCS(11) /    .0000000000 00017990E0 /
      DATA ERFCS(12) /   -.0000000000 00000371E0 /
      DATA ERFCS(13) /    .0000000000 00000007E0 /
      DATA SQRTPI /1.772453850 9055160E0/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  ERF
      IF (FIRST) THEN
         NTERF = INITS (ERFCS, 13, 0.1*R1MACH(3))
         XBIG = SQRT(-LOG(SQRTPI*R1MACH(3)))
         SQEPS = SQRT(2.0*R1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.1.) GO TO 20
C
C ERF(X) = 1. - ERFC(X) FOR -1. .LE. X .LE. 1.
C
      IF (Y.LE.SQEPS) ERF = 2.0*X/SQRTPI
      IF (Y.GT.SQEPS) ERF = X*(1.0 + CSEVL(2.*X**2-1., ERFCS, NTERF))
      RETURN
C
C ERF(X) = 1. - ERFC(X) FOR  ABS(X) .GT. 1.
C
 20   IF (Y.LE.XBIG) ERF = SIGN (1.0-ERFC(Y), X)
      IF (Y.GT.XBIG) ERF = SIGN (1.0, X)
C
      RETURN
      END
*DECK CSEVL
      FUNCTION CSEVL (X, CS, N)
C***BEGIN PROLOGUE  CSEVL
C***PURPOSE  Evaluate a Chebyshev series.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C3A2
C***TYPE      SINGLE PRECISION (CSEVL-S, DCSEVL-D)
C***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C  Evaluate the N-term Chebyshev series CS at X.  Adapted from
C  a method presented in the paper by Broucke referenced below.
C
C       Input Arguments --
C  X    value at which the series is to be evaluated.
C  CS   array of N terms of a Chebyshev series.  In evaluating
C       CS, only half the first coefficient is summed.
C  N    number of terms in array CS.
C
C***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
C                 Chebyshev series, Algorithm 446, Communications of
C                 the A.C.M. 16, (1973) pp. 254-256.
C               L. Fox and I. B. Parker, Chebyshev Polynomials in
C                 Numerical Analysis, Oxford University Press, 1968,
C                 page 56.
C***ROUTINES CALLED  R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900329  Prologued revised extensively and code rewritten to allow
C           X to be slightly outside interval (-1,+1).  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CSEVL
      REAL B0, B1, B2, CS(*), ONEPL, TWOX, X
      LOGICAL FIRST
      SAVE FIRST, ONEPL
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  CSEVL
      IF (FIRST) ONEPL = 1.0E0 + R1MACH(4)
      FIRST = .FALSE.
      IF (N .LT. 1) CALL XERMSG ('SLATEC', 'CSEVL',
     +   'NUMBER OF TERMS .LE. 0', 2, 2)
      IF (N .GT. 1000) CALL XERMSG ('SLATEC', 'CSEVL',
     +   'NUMBER OF TERMS .GT. 1000', 3, 2)
      IF (ABS(X) .GT. ONEPL) CALL XERMSG ('SLATEC', 'CSEVL',
     +   'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
C
      B1 = 0.0E0
      B0 = 0.0E0
      TWOX = 2.0*X
      DO 10 I = 1,N
         B2 = B1
         B1 = B0
         NI = N + 1 - I
         B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
C
      CSEVL = 0.5E0*(B0-B2)
C
      RETURN
      END
*DECK INITS
      FUNCTION INITS (OS, NOS, ETA)
C***BEGIN PROLOGUE  INITS
C***PURPOSE  Determine the number of terms needed in an orthogonal
C            polynomial series so that it meets a specified accuracy.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C3A2
C***TYPE      SINGLE PRECISION (INITS-S, INITDS-D)
C***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
C             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C  Initialize the orthogonal series, represented by the array OS, so
C  that INITS is the number of terms needed to insure the error is no
C  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
C  machine precision.
C
C             Input Arguments --
C   OS     single precision array of NOS coefficients in an orthogonal
C          series.
C   NOS    number of coefficients in OS.
C   ETA    single precision scalar containing requested accuracy of
C          series.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891115  Modified error message.  (WRB)
C   891115  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  INITS
      REAL OS(*)
C***FIRST EXECUTABLE STATEMENT  INITS
      IF (NOS .LT. 1) CALL XERMSG ('SLATEC', 'INITS',
     +   'Number of coefficients is less than 1', 2, 1)
C
      ERR = 0.
      DO 10 II = 1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(OS(I))
        IF (ERR.GT.ETA) GO TO 20
   10 CONTINUE
C
   20 IF (I .EQ. NOS) CALL XERMSG ('SLATEC', 'INITS',
     +   'Chebyshev series too short for specified accuracy', 1, 1)
      INITS = I
C
      RETURN
      END
*DECK DERF
      DOUBLE PRECISION FUNCTION DERF (X)
C***BEGIN PROLOGUE  DERF
C***PURPOSE  Compute the error function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C8A, L5A1E
C***TYPE      DOUBLE PRECISION (ERF-S, DERF-D)
C***KEYWORDS  ERF, ERROR FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DERF(X) calculates the double precision error function for double
C precision argument X.
C
C Series for ERF        on the interval  0.          to  1.00000E+00
C                                        with weighted error   1.28E-32
C                                         log weighted error  31.89
C                               significant figures required  31.05
C                                    decimal places required  32.55
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, DERFC, INITDS
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900727  Added EXTERNAL statement.  (WRB)
C   920618  Removed space from variable name.  (RWC, WRB)
C***END PROLOGUE  DERF
      DOUBLE PRECISION X, ERFCS(21), SQEPS, SQRTPI, XBIG, Y, D1MACH,
     1  DCSEVL, DERFC
      LOGICAL FIRST
      EXTERNAL DERFC
      SAVE ERFCS, SQRTPI, NTERF, XBIG, SQEPS, FIRST
      DATA ERFCS(  1) / -.4904612123 4691808039 9845440333 76 D-1     /
      DATA ERFCS(  2) / -.1422612051 0371364237 8247418996 31 D+0     /
      DATA ERFCS(  3) / +.1003558218 7599795575 7546767129 33 D-1     /
      DATA ERFCS(  4) / -.5768764699 7674847650 8270255091 67 D-3     /
      DATA ERFCS(  5) / +.2741993125 2196061034 4221607914 71 D-4     /
      DATA ERFCS(  6) / -.1104317550 7344507604 1353812959 05 D-5     /
      DATA ERFCS(  7) / +.3848875542 0345036949 9613114981 74 D-7     /
      DATA ERFCS(  8) / -.1180858253 3875466969 6317518015 81 D-8     /
      DATA ERFCS(  9) / +.3233421582 6050909646 4029309533 54 D-10    /
      DATA ERFCS( 10) / -.7991015947 0045487581 6073747085 95 D-12    /
      DATA ERFCS( 11) / +.1799072511 3961455611 9672454866 34 D-13    /
      DATA ERFCS( 12) / -.3718635487 8186926382 3168282094 93 D-15    /
      DATA ERFCS( 13) / +.7103599003 7142529711 6899083946 66 D-17    /
      DATA ERFCS( 14) / -.1261245511 9155225832 4954248533 33 D-18    /
      DATA ERFCS( 15) / +.2091640694 1769294369 1705002666 66 D-20    /
      DATA ERFCS( 16) / -.3253973102 9314072982 3641600000 00 D-22    /
      DATA ERFCS( 17) / +.4766867209 7976748332 3733333333 33 D-24    /
      DATA ERFCS( 18) / -.6598012078 2851343155 1999999999 99 D-26    /
      DATA ERFCS( 19) / +.8655011469 9637626197 3333333333 33 D-28    /
      DATA ERFCS( 20) / -.1078892517 7498064213 3333333333 33 D-29    /
      DATA ERFCS( 21) / +.1281188399 3017002666 6666666666 66 D-31    /
      DATA SQRTPI / 1.772453850 9055160272 9816748334 115D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DERF
      IF (FIRST) THEN
         NTERF = INITDS (ERFCS, 21, 0.1*REAL(D1MACH(3)))
         XBIG = SQRT(-LOG(SQRTPI*D1MACH(3)))
         SQEPS = SQRT(2.0D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.1.D0) GO TO 20
C
C ERF(X) = 1.0 - ERFC(X)  FOR  -1.0 .LE. X .LE. 1.0
C
      IF (Y.LE.SQEPS) DERF = 2.0D0*X*X/SQRTPI
      IF (Y.GT.SQEPS) DERF = X*(1.0D0 + DCSEVL (2.D0*X*X-1.D0,
     1  ERFCS, NTERF))
      RETURN
C
C ERF(X) = 1.0 - ERFC(X) FOR ABS(X) .GT. 1.0
C
 20   IF (Y.LE.XBIG) DERF = SIGN (1.0D0-DERFC(Y), X)
      IF (Y.GT.XBIG) DERF = SIGN (1.0D0, X)
C
      RETURN
      END
*DECK DCSEVL
      DOUBLE PRECISION FUNCTION DCSEVL (X, CS, N)
C***BEGIN PROLOGUE  DCSEVL
C***PURPOSE  Evaluate a Chebyshev series.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C3A2
C***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
C***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C  Evaluate the N-term Chebyshev series CS at X.  Adapted from
C  a method presented in the paper by Broucke referenced below.
C
C       Input Arguments --
C  X    value at which the series is to be evaluated.
C  CS   array of N terms of a Chebyshev series.  In evaluating
C       CS, only half the first coefficient is summed.
C  N    number of terms in array CS.
C
C***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
C                 Chebyshev series, Algorithm 446, Communications of
C                 the A.C.M. 16, (1973) pp. 254-256.
C               L. Fox and I. B. Parker, Chebyshev Polynomials in
C                 Numerical Analysis, Oxford University Press, 1968,
C                 page 56.
C***ROUTINES CALLED  D1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900329  Prologued revised extensively and code rewritten to allow
C           X to be slightly outside interval (-1,+1).  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DCSEVL
      DOUBLE PRECISION B0, B1, B2, CS(*), ONEPL, TWOX, X, D1MACH
      LOGICAL FIRST
      SAVE FIRST, ONEPL
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DCSEVL
      IF (FIRST) ONEPL = 1.0D0 + D1MACH(4)
      FIRST = .FALSE.
      IF (N .LT. 1) CALL XERMSG ('SLATEC', 'DCSEVL',
     +   'NUMBER OF TERMS .LE. 0', 2, 2)
      IF (N .GT. 1000) CALL XERMSG ('SLATEC', 'DCSEVL',
     +   'NUMBER OF TERMS .GT. 1000', 3, 2)
      IF (ABS(X) .GT. ONEPL) CALL XERMSG ('SLATEC', 'DCSEVL',
     +   'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
C
      B1 = 0.0D0
      B0 = 0.0D0
      TWOX = 2.0D0*X
      DO 10 I = 1,N
         B2 = B1
         B1 = B0
         NI = N + 1 - I
         B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
C
      DCSEVL = 0.5D0*(B0-B2)
C
      RETURN
      END
*DECK INITDS
      FUNCTION INITDS (OS, NOS, ETA)
C***BEGIN PROLOGUE  INITDS
C***PURPOSE  Determine the number of terms needed in an orthogonal
C            polynomial series so that it meets a specified accuracy.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C3A2
C***TYPE      DOUBLE PRECISION (INITS-S, INITDS-D)
C***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
C             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C  Initialize the orthogonal series, represented by the array OS, so
C  that INITDS is the number of terms needed to insure the error is no
C  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
C  machine precision.
C
C             Input Arguments --
C   OS     double precision array of NOS coefficients in an orthogonal
C          series.
C   NOS    number of coefficients in OS.
C   ETA    single precision scalar containing requested accuracy of
C          series.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891115  Modified error message.  (WRB)
C   891115  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  INITDS
      DOUBLE PRECISION OS(*)
C***FIRST EXECUTABLE STATEMENT  INITDS
      IF (NOS .LT. 1) CALL XERMSG ('SLATEC', 'INITDS',
     +   'Number of coefficients is less than 1', 2, 1)
C
      ERR = 0.
      DO 10 II = 1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(REAL(OS(I)))
        IF (ERR.GT.ETA) GO TO 20
   10 CONTINUE
C
   20 IF (I .EQ. NOS) CALL XERMSG ('SLATEC', 'INITDS',
     +   'Chebyshev series too short for specified accuracy', 1, 1)
      INITDS = I
C
      RETURN
      END
