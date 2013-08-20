#ifndef _ChemDriver_F_H_
#define _ChemDriver_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(CD_INITCHEM, cd_initchem)
   (int& len_en, int& len_sn, int& ne, int& ns);
BL_FORT_PROC_DECL(CD_CLOSECHEM, cd_closechem)();

BL_FORT_PROC_DECL(CD_GETELEMNAME, cd_getelemname)
   (const int& ielem, int& lname, int* elem_name);
BL_FORT_PROC_DECL(CD_GETSPECNAME, cd_getspecname)
   (const int& ispec, int& lname, int* spec_name);

BL_FORT_PROC_DECL(CD_INITVODE, cd_initvode)
   (const int& neq, const int& itol, const Real& rtol, const Real& atol,
    const int& order, const int& maxstep, const int& use_ajac, const int& save_ajac,
    const int& always_new_j, const int& stiff, const int& verbose);
BL_FORT_PROC_DECL(CD_CLOSEVODE, cd_closevode)();

BL_FORT_PROC_DECL(CD_INITEGLIB, cd_initeglib)
   (const int& use_bulk_visc);
BL_FORT_PROC_DECL(CD_CLOSEEGLIB, cd_closeeglib)();

#endif
