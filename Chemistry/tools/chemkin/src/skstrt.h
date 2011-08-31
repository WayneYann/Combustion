C       CVS Revision:$Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $

C	this is CHEMKIN-III file skstrt.h V.3.0 January 1997;
C	it contains pointers into the surface kinetics
C       data storage arrays.
C
C
C     Include file for CHEMKIN-III sklib.f, dated:  March 1, 1966
C
      INTEGER
     1   MAXSPR, NELEM, NKKGAS, NSPAR, NSCOV, NEDPAR, NYPAR, MAXORD, 
     2   MAXTP, NCP,    NCP1,  NCP2,   NCP2T,
     3   IiLENI, IiLENR, IiLENC, IiKSUR, IiKBLK, IiKTOT, IiNPHA, IiFSUR,
     4   IiLSUR, IiNSUR, IiFBLK, IiLBLK, IiNBLK, IiNIIS, IiNCOV, IiNREV,
     5   IiNSTK, IiNCON, IiNBHM, IiNRNU, IiNORD, IiMOTZ, IiELEC, IiNYLD,
     6   IiPKST, IiPKND, IiPTOT, IiKPHS, IiKCHG, IiKCMP, IiNSCV, IiKNT,
     7   IiNRPP, IiNREA, IiNUNK, IiNU,   IiNSUM, IiICOV, IiKCOV, IiIREV,
     8   IiISTK, IiMSTK, IiIBHM, IiKBHM, IiIRNU, IiIORD, IiKORD, IiIONS, 
     9   IiKION, IiKTFL, IiNEDP, IiIEDP, IiKEDP, IiIYLD, IiYION, IiKYLD,
     *   IrSKMN, IrPATM, IrRU,   IrRUC,
     1   IrSDEN, IrKTMP, IrKTHM, IrKDEN, IrAWT,  IrKWT,  IrPAR, IrKCOV,
     2   IrRPAR, IrEQ,   IrRNU,  IrNCF,  IrKORD, IrKFT,  IrKRT, IrKT1,
     3   IrKT2,  IrPT1,  IrIT1,  IrIT2,  IrIT3,  IrPEDP, IrENGI,
     4   IrPYLD, IrYNCF,
     4   IcENAM, IcKNAM, IcMNAM, IcPNAM

      COMMON /SKSTRT/
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
