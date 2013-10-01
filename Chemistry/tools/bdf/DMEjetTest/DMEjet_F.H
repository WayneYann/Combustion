#ifndef _DMEJET_F_H_
#define _DMEJET_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(CHEMSOLVE, chemsolve)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(Uin),
     BL_FORT_FAB_ARG(Uou),
     const Real& dt, const int& verbose, const int& use_vode);

#endif
