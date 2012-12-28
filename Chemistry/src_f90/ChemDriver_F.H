#ifndef _Chemdriver_F_H_
#define _Chemdrivererive_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(CA_DERLAPVAR,ca_derlapvar)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

BL_FORT_PROC_DECL(CA_UMDRV,ca_umdrv)
    (const int    lo[], const int    hi[],
     const BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(stateout),
     const BL_FORT_FAB_ARG(ugdnvx),
#if (BL_SPACEDIM >= 2)
     const BL_FORT_FAB_ARG(ugdnvy),
#endif
#if (BL_SPACEDIM == 3)
     const BL_FORT_FAB_ARG(ugdnvz),
#endif
     const BL_FORT_FAB_ARG(src),
     const BL_FORT_FAB_ARG(grav),
     const Real dx[], const Real* dt,
     D_DECL(const BL_FORT_FAB_ARG(xflux),
            const BL_FORT_FAB_ARG(yflux),
            const BL_FORT_FAB_ARG(zflux)),
     D_DECL(const BL_FORT_FAB_ARG(xarea),
            const BL_FORT_FAB_ARG(yarea),
            const BL_FORT_FAB_ARG(zarea)),
#if (BL_SPACEDIM < 3)
     const BL_FORT_FAB_ARG(dloga),
#endif
     const BL_FORT_FAB_ARG(volume),
     const Real* cflLoc,
     const int& verbose);


#endif
