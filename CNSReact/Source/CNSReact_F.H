#ifndef _CNSReact_F_H_
#define _CNSReact_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(FILCC,filcc)
                (const Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                 const int * domlo, const int * domhi,
                 const Real * dx_crse, const Real * xlo,
                 const int * bc);

BL_FORT_PROC_DECL(CA_NETWORK_INIT,ca_network_init)();

BL_FORT_PROC_DECL(GET_NUM_SPEC,get_num_spec)(int* nspec);

BL_FORT_PROC_DECL(GET_NUM_AUX,get_num_aux)(int* naux);

BL_FORT_PROC_DECL(GET_SPEC_NAMES,get_spec_names)(int* spec_names, int* ispec, int* len);

BL_FORT_PROC_DECL(GET_AUX_NAMES,get_aux_names)(int* aux_names, int* iaux, int* len);

BL_FORT_PROC_DECL(GET_METHOD_PARAMS,get_method_params)(int* HYP_GROW);

BL_FORT_PROC_DECL(SET_METHOD_PARAMS,set_method_params)
    (const int& dm,
     const int& Density, const int& Xmom, 
     const int& Eden,    const int& Eint, 
     const int& Temp     , const int& FirstAdv, 
     const int& FirstSpec, const int& FirstAux, 
     const int& NumAdv, 
     const Real& small_dens, const Real& small_temp,
     const Real& small_pres, const int& allow_negative_energy,
     const int& ppm_type, const Real& gamma_in,
     const int& normalize_species);

BL_FORT_PROC_DECL(SET_PROBLEM_PARAMS,set_problem_params)
    (const int& dm,
     const int* physbc_lo, const int* physbc_hi,
     const int& Outflow_value, 
     const int& Symmetry_value, const int& coord_type);

BL_FORT_PROC_DECL(CA_INITDATA,ca_initdata)
    (const int& level, const Real& time, 
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG(state),
     const Real dx[], const Real xlo[], const Real xhi[]);

BL_FORT_PROC_DECL(CA_ENFORCE_CONSISTENT_E,ca_enforce_consistent_e)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG(state));

BL_FORT_PROC_DECL(CA_DENERROR,ca_denerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(den),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_VELERROR,ca_velerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(vel),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_VORTERROR,ca_vorterror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(vel),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_TEMPERROR,ca_temperror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(temp),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_PRESSERROR,ca_presserror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(press),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_MFRACERROR,ca_mfracerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(press),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_RADERROR,ca_raderror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(den),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_ENTERROR,ca_enterror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(den),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_STATE_ERROR,ca_state_error)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(state),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_HYPFILL,ca_hypfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_DENFILL,ca_denfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_XMOMFILL,ca_xmomfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_YMOMFILL,ca_ymomfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_ZMOMFILL,ca_zmomfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_ZEROFILL,ca_zerofill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_DEFAULT_FILL,ca_default_fill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

#ifdef REACTIONS
BL_FORT_PROC_DECL(CA_REACTFILL,ca_reactfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);
#endif

BL_FORT_PROC_DECL(CA_ESTDT,ca_estdt)
    (const BL_FORT_FAB_ARG(state),
     const int lo[], const int hi[],
     const Real dx[], Real* dt);

BL_FORT_PROC_DECL(CA_ESTDT_BURNING,ca_estdt_burning)
    (const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(reaction_terms),
     const int lo[], const int hi[],
     Real* dt_hydro, Real* dt_nuc, 
     int* initial_step);

BL_FORT_PROC_DECL(CA_UMDRV,ca_umdrv)
    (const int* is_finest_level,
     const Real* time, 
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
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

BL_FORT_PROC_DECL(CA_CORRGSRC,ca_corrgsrc)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(grav_src_old),
     const BL_FORT_FAB_ARG(grad_phi_cc),
     const BL_FORT_FAB_ARG(S_old),
     BL_FORT_FAB_ARG(S_new),
     const Real& dt);

BL_FORT_PROC_DECL(CA_RESET_INTERNAL_ENERGY,ca_reset_internal_energy)
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
BL_FORT_PROC_DECL(CA_SUMMASS,ca_summass)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real dx[], Real* mass);
BL_FORT_PROC_DECL(CA_SUMSQUARED,ca_sumsquared)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real dx[], Real* mass);
BL_FORT_PROC_DECL(CA_SUMLOCMASS,ca_sumlocmass)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real* problo, const Real dx[], Real* mass, const int& idir);
#else
BL_FORT_PROC_DECL(CA_SUMMASS,ca_summass)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real dx[], Real* mass,
     const Real* r, const int& rlo, const int& rhi);
BL_FORT_PROC_DECL(CA_SUMSQUARED,ca_sumsquared)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real dx[], Real* mass,
     const Real* r, const int& rlo, const int& rhi);
BL_FORT_PROC_DECL(CA_SUMLOCMASS,ca_sumlocmass)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real* problo, const Real dx[], Real* mass, 
     const Real* r, const int& rlo, const int& rhi,
     const int& idir);
#endif

BL_FORT_PROC_DECL(CA_ENFORCE_NONNEGATIVE_SPECIES,ca_enforce_nonnegative_species)
    (BL_FORT_FAB_ARG(S_new), const int lo[], const int hi[]);

BL_FORT_PROC_DECL(CA_AVGDOWN,ca_avgdown)
    (BL_FORT_FAB_ARG(crse_fab), const int& nc,
     const BL_FORT_FAB_ARG(crse_vol),
     const BL_FORT_FAB_ARG(fine_fab),
     const BL_FORT_FAB_ARG(fine_vol),
     const int ovlo[], const int ovhi[], const int rat[]);

BL_FORT_PROC_DECL(CA_EXT_SRC,ca_ext_src)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(old_state),
     BL_FORT_FAB_ARG(new_state),
     BL_FORT_FAB_ARG(ext_src),
     const Real* dx, const Real* time, const Real* dt);

BL_FORT_PROC_DECL(CA_REACT_STATE,ca_react_state)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG(Sold),
     const BL_FORT_FAB_ARG(Snew),
     const BL_FORT_FAB_ARG(ReactionTerms),
     const Real& time, const Real& dt_react);

BL_FORT_PROC_DECL(CA_COMPUTE_TEMP,ca_compute_temp)
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
