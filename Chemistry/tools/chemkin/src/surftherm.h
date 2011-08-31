C       CVS Revision:$Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
*@(#)==============================================
*@(#) SCCS File Name surftherm.h
*@(#) SCCS Version = 1.2
*@(#) SCCS Modification date = 05/06/96
*@(#)==============================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C       Globally Defined Variables for the Surftherm Program
C
C----------------------------------------------------------------------C
C          SURFTHERM.H
C                    Version 1.10
C
C          WRITTEN BY H. K. Moffat, 1126
C
C----------------------------------------------------------------------C
C
* Parameter specifying the maximum number of temperatures in a table
      INTEGER     MAXTEMP
      PARAMETER  (MAXTEMP=50)
C----------------------------------------------------------------------C
C
C      Key Variables NEEDED FOR ARRAY DIMENSIONING
C
C Common block for Storage of the size of the chemistry problem
C
      INTEGER         LENRCK, LENICK, LENCCK, LENRSK, LENISK, LENCSK,
     1                LENIMC, LENRMC,
     1                NELEM,  KKGAS,  KKSURF, KKBULK, KKTOT,
     2                NNPHAS, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     3                NLBULK, NFIT,   IIGAS,  IISUR,  NFAR

      COMMON /SF_CHM/ LENRCK, LENICK, LENCCK, LENRSK, LENISK, LENCSK,
     1                LENIMC, LENRMC,
     1                NELEM,  KKGAS,  KKSURF, KKBULK, KKTOT,
     2                NNPHAS, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     3                NLBULK, NFIT,   IIGAS,  IISUR,  NFAR
C
C Common Block for Global Integer data
C
*
*  ID_CARR     Identity of the carrier gas species
*              - defaults to the last species in the mechanism
*  ID_MAJOR    Identity of the marjor reactant in the gas
*              - defaults to the first species in the mechanism
*
*  TTABLE_NUM  Number of temperatures in tables which raster over
*              the temperatures
*              - default (TTABLE_MAX - TTABLE_MIN)/TTABLE_DELTAT + 3
*  PTABLE_NUM  Number of pressures in tables which raster over
*              a list of pressures
*
*
      INTEGER ID_CARR, ID_MAJOR, TTABLE_NUM, PTABLE_NUM
*
      COMMON /SF_INT/ 
     $        ID_CARR, ID_MAJOR, TTABLE_NUM, PTABLE_NUM
C
C Common Block For Key real*8 data
C
*
* P_BATH     Pressure of the bath gas (units of torr)
*            - User-definable through keywords
*            - default = 760. torr
* P_CGS      Pressure of the bath gas (cgs units)
*            - Calculated from P_BATH
*
* T_BATH     Temperature of the bath gas (Kelvin)
*            - User-definable through keywords
*            - default = 298.15 Kelvin
*
* TTABLE_MIN Minimum temperature used in temperature tables
*            - User-defined through keywords
*            - default set to 300 K
* TTABLE_MAX Maximum temperature used in temperature tables
*            - User-defined through keywords
*            - default set to 1500 K
* TTABLE_DELTAT Delta temperature difference between rows
*            in temperature tables
*            - User-defined through keywords
*            - default set to 100 K
* PTABLE_MIN Minimum pressure to use in pressure tables (torr)
*            - User-defined through keywords
*            - default = 1 torr
* PTABLE_MAX Maximum pressure to use in pressure tables (torr)
*            - User-defined through keywords
*            - default = 1000 torr
* PTABLE_NUM Number of rows of items in the pressure tables
*            - User-defined through keywords
*            - default = 14
*              log scale is employed between MIN and MAX.
*
* LENGTH_SCALE Length scale to use in non-dimensional numbers
*              (cm)
*            - User-defined through keywords
*            - default = 1 cm
*
* RKCAL      Value of the gas constant in kcal/moleK
*            - Calculated from RUC
*            - value = 1.987E-3
*
* RU         Gas constant in CGS units
*            - Returned from Chemkin
*
* RUC        Gas constant in cals/moleK
*            - Returned from chemkin
*
* PATM      Conversion of pressure from atm to cgs units
*            - Returned from chemkin
*
* CTOT_BATH Concentration of the gas at the bath gas conditions
*           ( mole/cm**3)
C
C*****B-2) precision > double
      DOUBLE PRECISION
C*****END B-2) precision > double
C*****B-1) precision > single
C      REAL
C*****END B-1) precision > single
*
     1                 P_BATH, P_CGS, T_BATH, TTABLE_MIN,
     2                 TTABLE_MAX, TTABLE_DELTAT, PTABLE_MIN,
     3                 PTABLE_MAX, LENGTH_SCALE, RKCAL,
     4                 RU, RUC, PATM, CTOT_BATH
*
      COMMON /SF_REAL/ P_BATH, P_CGS, T_BATH, TTABLE_MIN,
     2                 TTABLE_MAX, TTABLE_DELTAT, PTABLE_MIN,
     3                 PTABLE_MAX, LENGTH_SCALE, RKCAL,
     4                 RU, RUC, PATM, CTOT_BATH
C
C Machine Common block from surface Chemkin
C
C*****B-2) precision > double
      DOUBLE PRECISION
C*****END B-2) precision > double
C*****B-1) precision > single
C      REAL
C*****END B-1) precision > single
     1                SMALL, BIG, EXPARG
      COMMON /MACHN / SMALL, BIG, EXPARG
