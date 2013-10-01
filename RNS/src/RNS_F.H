// This file needs cleanup

#ifndef _RNS_F_H_
#define _RNS_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(RNS_DUDT,rns_dudt)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(U),
     BL_FORT_FAB_ARG(dUdt),
     D_DECL(BL_FORT_FAB_ARG(xflux),
	    BL_FORT_FAB_ARG(yflux),
	    BL_FORT_FAB_ARG(zflux)),
     const Real dx[]);

BL_FORT_PROC_DECL(RNS_ADVCHEM, rns_advchem)
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(U),
     const Real& dt);

BL_FORT_PROC_DECL(RNS_COMPUTE_TEMP,rns_compute_temp)
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(U));

BL_FORT_PROC_DECL(RNS_ENFORCE_CONSISTENT_Y,rns_enforce_consistent_y)
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(U));

BL_FORT_PROC_DECL(FILCC,filcc)
                (const Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                 const int * domlo, const int * domhi,
                 const Real * dx_crse, const Real * xlo,
                 const int * bc);

BL_FORT_PROC_DECL(GET_METHOD_PARAMS,get_method_params)
    (const int* ngrow, const int* nriemann);

BL_FORT_PROC_DECL(SET_METHOD_PARAMS,set_method_params)
    (const int& dm, const int& Density, const int& Xmom, 
     const int& Eden, const int& Temp, const int& FirstSpec, 
     const int& NUM_STATE, const int& NumSpec,
     const Real& small_dens, const Real& small_temp, const Real& small_pres,
     const Real& gamma, const Real& gravity, const Real& Treference,
     const int& riemann, const Real& difmag, const int* blocksize,
     const int& use_vode, const int& do_cc_burning);

BL_FORT_PROC_DECL(SET_PROBLEM_PARAMS,set_problem_params)
    (const int& dm,
     const int* physbc_lo, const int* physbc_hi,
     const Real* phys_lo, const Real* phys_hi,
     const int& Outflow_value, 
     const int& Symmetry_value, const int& coord_type);

BL_FORT_PROC_DECL(RNS_INITDATA,rns_initdata)
    (const int& level, const Real& time, 
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG(state),
     const Real dx[], const Real xlo[], const Real xhi[]);

BL_FORT_PROC_DECL(RNS_DENERROR,rns_denerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(den),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(RNS_VELERROR,rns_velerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(vel),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(RNS_VORTERROR,rns_vorterror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(vel),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(RNS_TEMPERROR,rns_temperror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(temp),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(RNS_PRESSERROR,rns_presserror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(press),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(RNS_TRACERROR,rns_tracerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(press),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(RNS_GRPFILL,rns_grpfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(RNS_REGFILL,rns_regfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(RNS_DENFILL,rns_denfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(RNS_MXFILL,rns_mxfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(RNS_MYFILL,rns_myfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(RNS_MZFILL,rns_mzfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(RNS_ESTDT,rns_estdt)
    (const BL_FORT_FAB_ARG(state),
     const int lo[], const int hi[],
     const Real dx[], Real* dt);

BL_FORT_PROC_DECL(RNS_AVGDOWN,rns_avgdown)
    (BL_FORT_FAB_ARG(crse_fab), const int& nc,
     const BL_FORT_FAB_ARG(crse_vol),
     const BL_FORT_FAB_ARG(fine_fab),
     const BL_FORT_FAB_ARG(fine_vol),
     const int ovlo[], const int ovhi[], const int rat[]);


/* problem-specific stuff goes here */

BL_FORT_PROC_DECL(PROBLEM_CHECKPOINT,problem_checkpoint)(int * int_dir_name, int * len);

BL_FORT_PROC_DECL(PROBLEM_RESTART,problem_restart)(int * int_dir_name, int * len);

/* for plotfiles */

BL_FORT_PROC_DECL(RNS_CTOPRIM, rns_ctoprim)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(cons),
     BL_FORT_FAB_ARG(prim));

BL_FORT_PROC_DECL(RNS_MAKEPLOTVAR, rns_makeplotvar)
    (const int lo[], const int hi[], const Real dx[],
     const BL_FORT_FAB_ARG(prim),
     BL_FORT_FAB_ARG(plot),
     const int& npv, const int& icomp_magvel, const int& icomp_Mach, 
     const int& icomp_divu, const int& icomp_magvort, 
     const int& icomp_X, const int& icomp_omegadot, const int& icomp_dYdt, 
     const int& icomp_heatRelease, const int& icomp_fuelConsumption,
     const int& fuelID);

#endif
