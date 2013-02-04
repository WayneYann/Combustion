#ifndef _CNSReact_F_H_
#define _CNSReact_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(FILCC,filcc)
                (const Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                 const int * domlo, const int * domhi,
                 const Real * dx_crse, const Real * xlo,
                 const int * bc);

BL_FORT_PROC_DECL(GET_METHOD_PARAMS,get_method_params)
    (const int* nHYP);

BL_FORT_PROC_DECL(SET_METHOD_PARAMS,set_method_params)
    (const int& dm,
     const int& Density, const int& Xmom, 
     const int& Eden,    const int& Eint, 
     const int& Temp     , const int& FirstAdv, 
     const int& FirstSpec, const int& NumAdv, 
     const Real& small_dens, const Real& small_temp,
     const Real& small_pres, const int& ppm_type, 
     const int& normalize_species);

BL_FORT_PROC_DECL(SET_PROBLEM_PARAMS,set_problem_params)
    (const int& dm,
     const int* physbc_lo, const int* physbc_hi,
     const Real* phys_lo, const Real* phys_hi,
     const int& Outflow_value, 
     const int& Symmetry_value, const int& coord_type,
     const Real& gravx, const Real& gravy, const Real& gravz);

BL_FORT_PROC_DECL(CNS_INITDATA,cns_initdata)
    (const int& level, const Real& time, 
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG(state),
     const Real dx[], const Real xlo[], const Real xhi[]);

BL_FORT_PROC_DECL(CNS_ENFORCE_CONSISTENT_E,cns_enforce_consistent_e)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG(state));

BL_FORT_PROC_DECL(CNS_DENERROR,cns_denerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(den),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CNS_VELERROR,cns_velerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(vel),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CNS_VORTERROR,cns_vorterror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(vel),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CNS_TEMPERROR,cns_temperror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(temp),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CNS_PRESSERROR,cns_presserror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(press),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CNS_MFRACERROR,cns_mfracerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(press),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CNS_ENTERROR,cns_enterror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(den),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CNS_STATE_ERROR,cns_state_error)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(state),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CNS_HYPFILL,cns_hypfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CNS_DENFILL,cns_denfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CNS_ESTDT,cns_estdt)
    (const BL_FORT_FAB_ARG(state),
     const int lo[], const int hi[],
     const Real dx[], Real* dt);

BL_FORT_PROC_DECL(CNS_UMDRV,cns_umdrv)
    (const int    lo[], const int    hi[],
     BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(stateout),
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

BL_FORT_PROC_DECL(CA_CORRGSRC,ca_corrgsrc)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(grav_src_old),
     const BL_FORT_FAB_ARG(grad_phi_cc),
     const BL_FORT_FAB_ARG(S_old),
     BL_FORT_FAB_ARG(S_new),
     const Real& dt);

BL_FORT_PROC_DECL(CNS_RESET_INTERNAL_ENERGY,cns_reset_internal_energy)
    (BL_FORT_FAB_ARG(S_new), 
     const int lo[], const int hi[], 
     const int& verbose);

BL_FORT_PROC_DECL(CA_SYNCGSRC,ca_syncgsrc)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(gPhi),
     const BL_FORT_FAB_ARG(gdPhi),
     const BL_FORT_FAB_ARG(S),
     const BL_FORT_FAB_ARG(dS),
     BL_FORT_FAB_ARG(src),
     const Real& dt);

#if (BL_SPACEDIM == 3)
BL_FORT_PROC_DECL(CNS_SUMMASS,cns_summass)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real dx[], Real* mass);
#else
BL_FORT_PROC_DECL(CNS_SUMMASS,cns_summass)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real dx[], Real* mass,
     const Real* r, const int& rlo, const int& rhi);
#endif

BL_FORT_PROC_DECL(CNS_ENFORCE_NONNEGATIVE_SPECIES,cns_enforce_nonnegative_species)
    (BL_FORT_FAB_ARG(S_new), const int lo[], const int hi[]);

BL_FORT_PROC_DECL(CNS_AVGDOWN,cns_avgdown)
    (BL_FORT_FAB_ARG(crse_fab), const int& nc,
     const BL_FORT_FAB_ARG(crse_vol),
     const BL_FORT_FAB_ARG(fine_fab),
     const BL_FORT_FAB_ARG(fine_vol),
     const int ovlo[], const int ovhi[], const int rat[]);

BL_FORT_PROC_DECL(CA_REACT_STATE,ca_react_state)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG(Sold),
     const BL_FORT_FAB_ARG(Snew),
     const BL_FORT_FAB_ARG(ReactionTerms),
     const Real& time, const Real& dt_react);

BL_FORT_PROC_DECL(CNS_COMPUTE_TEMP,cns_compute_temp)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state));

BL_FORT_PROC_DECL(CA_FILL_TEMP_COND,ca_fill_temp_cond)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     D_DECL(const BL_FORT_FAB_ARG(xcoeffs),
            const BL_FORT_FAB_ARG(ycoeffs),
            const BL_FORT_FAB_ARG(zcoeffs)),
     const Real* dx);

BL_FORT_PROC_DECL(CA_SET_SPECIAL_TAGGING_FLAG,ca_set_special_tagging_flag)
    (const Real& max_den, const int* flag);

BL_FORT_PROC_DECL(CA_GET_RHOE,ca_get_rhoe)
  (const int lo[], const int hi[],
   const BL_FORT_FAB_ARG(rhoe),
   const BL_FORT_FAB_ARG(temp),
   const BL_FORT_FAB_ARG(ye),
   const BL_FORT_FAB_ARG(state));

#if (BL_SPACEDIM > 1)
BL_FORT_PROC_DECL(CA_COMPUTE_AVGSTATE,ca_compute_avgstate)
    (const int lo[], const int hi[], 
     const Real* dx, const Real* dr, const int* nc,
     const BL_FORT_FAB_ARG(S  ), const Real* avgden,
     const BL_FORT_FAB_ARG(Vol), const Real* avgvol,
     const Real* problo, const int* numpts_1d);
#endif

BL_FORT_PROC_DECL(CA_APPLY_METRIC,ca_apply_metric)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     D_DECL(const BL_FORT_FAB_ARG(xedge),
            const BL_FORT_FAB_ARG(yedge),
            const BL_FORT_FAB_ARG(zedge)), 
            const Real* dx, const int* coord_type);

BL_FORT_PROC_DECL(CA_UNWEIGHT_CC,ca_unweight_cc)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     const Real* dx, const int* coord_type);

#endif
