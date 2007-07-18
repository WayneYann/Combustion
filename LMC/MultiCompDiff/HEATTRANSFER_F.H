
#ifndef _HEATTRANSFER_H_
#define _HEATTRANSFER_H_  $Id: HEATTRANSFER_F.H,v 1.8 2007-07-18 20:19:38 sepp Exp $

#ifdef BL_LANG_FORT
#    define FORT_SET_SCAL_NUMB             setscnum
#    define FORT_GETPAMB                   getpamb
#    define FORT_COMPUTE_RHODGRADHDOTGRADY computerhodgradhdotgrady
#    define FORT_FLOOR_AND_SCALE_SPEC      floorandscalespec 
#    define FORT_DIVUFILL                  divufill
#    define FORT_DSDTFILL                  dsdtfill
#    define FORT_DQRADFILL                 dqradfill
#    define FORT_YDOTFILL                  ydotfill
#    define FORT_EST_DIVU_DT               estdivudt
#    define FORT_CHECK_DIVU_DT             checkdivudt
#    define FORT_INIT_VISC_DATA            initviscdata
#    define FORT_SET_HT_VISC_COMMON        setvisccm
#    define FORT_VELVISC                   velvisc
#    define FORT_TEMPVISC                  tempvisc
#    define FORT_SPECVISC                  specvisc
#    define FORT_SPECTEMPVISC              spectempvisc
#    define FORT_FABMINMAX                 fabminmax
#    define FORT_REPAIR_FLUX               rparflux
#    define FORT_RECOMP_UPDATE             recompupdt
#    define FORT_INCRWEXTFLXDIV            incrwextflxdiv
#    define FORT_SET_COMMON                setcommon
#    define FORT_COMPUTE_UGRADP            computeugradp
#    define FORT_RK_MIXTURE_AVERAGED       rkmixtureaveraged
#    define FORT_RK_MULTICOMPONENT         rkmulticomponent
#    define FORT_RK_STEP_SELECTION         rkstepselection
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_SET_SCAL_NUMB             SETSCNUM
#    define FORT_GETPAMB                   GETPAMB
#    define FORT_COMPUTE_RHODGRADHDOTGRADY COMPUTERHODGRADHDOTGRADY
#    define FORT_FLOOR_AND_SCALE_SPEC      FLOORANDSCALESPEC 
#    define FORT_DIVUFILL                  DIVUFILL
#    define FORT_DSDTFILL                  DSDTFILL
#    define FORT_DQRADFILL                 DQRADFILL
#    define FORT_YDOTFILL                  YDOTFILL
#    define FORT_EST_DIVU_DT               ESTDIVUDT
#    define FORT_CHECK_DIVU_DT             CHECKDIVUDT
#    define FORT_INIT_VISC_DATA            INITVISCDATA
#    define FORT_SET_HT_VISC_COMMON        SETVISCCM
#    define FORT_VELVISC                   VELVISC
#    define FORT_TEMPVISC                  TEMPVISC
#    define FORT_SPECVISC                  SPECVISC
#    define FORT_SPECTEMPVISC              SPECTEMPVISC
#    define FORT_FABMINMAX                 FABMINMAX
#    define FORT_REPAIR_FLUX               RPARFLUX
#    define FORT_RECOMP_UPDATE             RECOMPUPDT
#    define FORT_INCRWEXTFLXDIV            INCRWEXTFLXDIV
#    define FORT_SET_COMMON                SETCOMMON
#    define FORT_COMPUTE_UGRADP            COMPUTEUGRADP
#    define FORT_RK_MIXTURE_AVERAGED       RKMIXTUREAVERAGED
#    define FORT_RK_MULTICOMPONENT         RKMULTICOMPONENT
#    define FORT_RK_STEP_SELECTION         RKSTEPSELECTION
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_SET_SCAL_NUMB             setscnum
#    define FORT_GETPAMB                   getpamb
#    define FORT_COMPUTE_RHODGRADHDOTGRADY computerhodgradhdotgrady
#    define FORT_FLOOR_AND_SCALE_SPEC      floorandscalespec
#    define FORT_DIVUFILL                  divufill
#    define FORT_DSDTFILL                  dsdtfill
#    define FORT_DQRADFILL                 dqradfill
#    define FORT_YDOTFILL                  ydotfill
#    define FORT_EST_DIVU_DT               estdivudt
#    define FORT_CHECK_DIVU_DT             checkdivudt
#    define FORT_INIT_VISC_DATA            initviscdata
#    define FORT_SET_HT_VISC_COMMON        setvisccm
#    define FORT_VELVISC                   velvisc
#    define FORT_TEMPVISC                  tempvisc
#    define FORT_SPECVISC                  specvisc
#    define FORT_SPECTEMPVISC              spectempvisc
#    define FORT_FABMINMAX                 fabminmax
#    define FORT_REPAIR_FLUX               rparflux
#    define FORT_RECOMP_UPDATE             recompupdt
#    define FORT_INCRWEXTFLXDIV            incrwextflxdiv
#    define FORT_SET_COMMON                setcommon
#    define FORT_COMPUTE_UGRADP            computeugradp
#    define FORT_RK_MIXTURE_AVERAGED       rkmixtureaveraged
#    define FORT_RK_MULTICOMPONENT         rkmulticomponent
#    define FORT_RK_STEP_SELECTION         rkstepselection
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_SET_SCAL_NUMB             setscnum_
#    define FORT_GETPAMB                   getpamb_
#    define FORT_COMPUTE_RHODGRADHDOTGRADY computerhodgradhdotgrady_
#    define FORT_FLOOR_AND_SCALE_SPEC      floorandscalespec_
#    define FORT_DIVUFILL                  divufill_
#    define FORT_DSDTFILL                  dsdtfill_
#    define FORT_DQRADFILL                 dqradfill_
#    define FORT_YDOTFILL                  ydotfill_
#    define FORT_EST_DIVU_DT               estdivudt_
#    define FORT_CHECK_DIVU_DT             checkdivudt_
#    define FORT_INIT_VISC_DATA            initviscdata_
#    define FORT_SET_HT_VISC_COMMON        setvisccm_
#    define FORT_VELVISC                   velvisc_
#    define FORT_TEMPVISC                  tempvisc_
#    define FORT_SPECVISC                  specvisc_
#    define FORT_SPECTEMPVISC              spectempvisc_
#    define FORT_FABMINMAX                 fabminmax_
#    define FORT_REPAIR_FLUX               rparflux_
#    define FORT_RECOMP_UPDATE             recompupdt_
#    define FORT_INCRWEXTFLXDIV            incrwextflxdiv_
#    define FORT_SET_COMMON                setcommon_
#    define FORT_COMPUTE_UGRADP            computeugradp_
#    define FORT_RK_MIXTURE_AVERAGED       rkmixtureaveraged_
#    define FORT_RK_MULTICOMPONENT         rkmulticomponent_
#    define FORT_RK_STEP_SELECTION         rkstepselection_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_SET_SCAL_NUMB(const int* Density,
			    const int* Temp,
			    const int* Trac,
			    const int* RhoH,
			    const int* FirstSpec,
			    const int* LastSpec);
    
    void FORT_GETPAMB(Real* pamb, Real* dpdt_factor);

    void FORT_COMPUTE_RHODGRADHDOTGRADY(const Real* dx,
              const int* lo, const int* hi,
              ARLIM_P(speclo),ARLIM_P(spechi), const Real* specdat,
              ARLIM_P(hlo),ARLIM_P(hhi), const Real* hdat,
              ARLIM_P(betaxlo),ARLIM_P(betaxhi), const Real* betax,
              ARLIM_P(betaylo),ARLIM_P(betayhi), const Real* betay,
#if (BL_SPACEDIM == 3) 
              ARLIM_P(betazlo),ARLIM_P(betazhi), const Real* betaz,
#endif            
              ARLIM_P(ddrhslo),ARLIM_P(ddrhshi), Real* ddrhsdat);

    void FORT_FLOOR_AND_SCALE_SPEC(const int& nspecies, 
                    const int* lo, const int* hi, ARLIM_P(s_lo),
                    ARLIM_P(s_hi), Real* sum, Real* species); 
 
    void FORT_DIVUFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
			 const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

    void FORT_DSDTFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

    void FORT_DQRADFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

    void FORT_YDOTFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

    void FORT_EST_DIVU_DT(const int& flag,  const Real* divu_dt_factor,
                        const Real* dx, 
                        const Real* sdat,ARLIM_P(slo),ARLIM_P(shi),
                        const Real* dsdtdat,
                        const Real* rhodat,ARLIM_P(rholo),ARLIM_P(rhohi),
                        const Real* vel,ARLIM_P(ulo),ARLIM_P(uhi),
                        const Real* volume, ARLIM_P (v_lo), ARLIM_P (v_hi), 
                        const Real* areax,
                        ARLIM_P(areax_lo),ARLIM_P(areax_hi),
                        const Real* areay,
                        ARLIM_P(areay_lo),ARLIM_P(areay_hi),
#if (BL_SPACEDIM == 3) 
                        const Real* areaz,
                        ARLIM_P(areaz_lo),ARLIM_P(areaz_hi),
#endif            
                        const int* lo,const int* hi,
                        Real* dt, const Real* min_rho_divu_ceiling);

    void FORT_CHECK_DIVU_DT(const int& flag, const Real* divu_dt_factor,
                        const Real* dx, 
                        const Real* sdat, ARLIM_P(slo), ARLIM_P(shi),
                        const Real* dsdtdat,
                        const Real* rhodat,ARLIM_P(rholo),ARLIM_P(rhohi),
                        const Real* vel,ARLIM_P(ulo),ARLIM_P(uhi),
                        const Real* volume, ARLIM_P (v_lo), ARLIM_P (v_hi), 
                        const Real* areax,
                        ARLIM_P(areax_lo),ARLIM_P(areax_hi),
                        const Real* areay,
                        ARLIM_P(areay_lo),ARLIM_P(areay_hi),
#if (BL_SPACEDIM == 3) 
                        const Real* areaz,
                        ARLIM_P(areaz_lo),ARLIM_P(areaz_hi),
#endif            
                        const int* lo,const int* hi,
                        const Real* dt, const Real* min_rho_divu_ceiling);

    void FORT_INIT_VISC_DATA();
    void FORT_SET_HT_VISC_COMMON(const int* muIsVar,     const Real* muVal,
				 const int* lambdaIsVar, const Real* lambdaVal,
				 const int* rhoDIsVar,   const Real* rhoDVal,
				 const Real* Pr, const Real* Sc, const int* unityLe);
      
    void FORT_VELVISC(const int* lo, const int* hi,
		      ARLIM_P(tlo), ARLIM_P(thi), const Real* tempdat,
		      ARLIM_P(slo), ARLIM_P(shi), const Real* specdat,
		      ARLIM_P(blo), ARLIM_P(bhi), Real* betadat);

    void FORT_TEMPVISC(const int* lo, const int* hi,
		       ARLIM_P(tlo), ARLIM_P(thi), const Real* temp,
		       ARLIM_P(slo), ARLIM_P(shi), const Real* spec,
		       ARLIM_P(blo), ARLIM_P(bhi), Real* betadat);
    
    void FORT_SPECVISC(const int* lo, const int* hi,
		       ARLIM_P(tlo), ARLIM_P(thi), const Real* tempdat,
		       ARLIM_P(slo), ARLIM_P(shi), const Real* specdat,
		       ARLIM_P(blo), ARLIM_P(bhi), Real* betadat,
		       const int* nCompD, const Real* P1ATM_MKS);

    void FORT_SPECTEMPVISC(const int* lo, const int* hi,
                           ARLIM_P(tlo), ARLIM_P(thi), const Real* tempdat,
                           ARLIM_P(slo), ARLIM_P(shi), const Real* specdat,
                           ARLIM_P(blo), ARLIM_P(bhi), Real* betadat,
                           const int* nCompD, const Real* P1ATM_MKS,
                           const int* do_temp, const int* do_VelVisc);

    void FORT_FABMINMAX(const int* lo, const int* hi,
			Real* fabdat, ARLIM_P(flo), ARLIM_P(fhi),
			const Real* min, const Real* max, const int* nc);
    
    void FORT_REPAIR_FLUX(const int* lo, const int* hi,
                          Real* flux,  ARLIM_P(flo), ARLIM_P(fhi),
                          const Real* state, ARLIM_P(slo), ARLIM_P(shi),
                          const int* nc, const int* dir);

    void FORT_RECOMP_UPDATE(const int* lo, const int* hi,
			    const Real* updat, ARLIM_P(uplo), ARLIM_P(uphi),
			    const Real* xflux, ARLIM_P(fxlo), ARLIM_P(fxhi),
			    const Real* yflux, ARLIM_P(fylo), ARLIM_P(fyhi),
#if BL_SPACEDIM == 3
			    const Real* zflux, ARLIM_P(fzlo), ARLIM_P(fzhi),
#endif
			    const Real* vol,   ARLIM_P(vlo),  ARLIM_P(vhi),
			    const int* nc);
    
    void FORT_INCRWEXTFLXDIV(const int* lo, const int* hi,
                             const Real* xflux, ARLIM_P(fxlo),  ARLIM_P(fxhi),
                             const Real* yflux, ARLIM_P(fylo),  ARLIM_P(fyhi),
#if BL_SPACEDIM == 3
                             const Real* zflux, ARLIM_P(fzlo),  ARLIM_P(fzhi),
#endif
                             const Real* sto,     ARLIM_P(stolo),    ARLIM_P(stohi),
                             Real*       stn,     ARLIM_P(stnlo),    ARLIM_P(stnhi),
                             const Real* vol,     ARLIM_P(vollo),    ARLIM_P(volhi),
                             const int* nc, const Real* dt);

    void FORT_SET_COMMON (const Real* time, const int* iteration);
    void FORT_COMPUTE_UGRADP(const Real* p, ARLIM_P(p_lo),  ARLIM_P(p_hi),
                             const Real* ugradp, ARLIM_P(up_lo),  ARLIM_P(up_hi),
                             const Real* umac, ARLIM_P(umac_lo),  ARLIM_P(umac_hi),
                             const Real* vmac, ARLIM_P(vmac_lo),  ARLIM_P(vmac_hi),
#if BL_SPACEDIM == 3
                             const Real* wmac, ARLIM_P(wmac_lo),  ARLIM_P(wmac_hi),
#endif
                             const int* lo, const int* hi, const Real* dx);

    void FORT_RK_MIXTURE_AVERAGED (
	const int* domain_lo, const int* domain_hi,
	const int* lo, const int* hi,
	Real* areax, ARLIM_P(areax_lo), ARLIM_P(areax_hi),
	Real* areay, ARLIM_P(areay_lo), ARLIM_P(areay_hi),
	const int* bc,
	const Real* dt,
	const Real* dx,
	int* index_of_firstY,
	int* index_of_lastY,
	int* index_of_rho,
	int* index_of_rhoH,
	int* index_of_T,
	Real* maximum_error,
	int* maximum_iterations,
	int* ncomps,
	int* nspcecies,
	Real* state, ARLIM_P(state_lo), ARLIM_P(state_hi),
	Real* update_for_H, ARLIM_P(update_for_H_lo), ARLIM_P(update_for_H_hi),
	Real* update_for_Y, ARLIM_P(update_for_Y_lo), ARLIM_P(update_for_Y_hi),
	Real* volume, ARLIM_P(volume_lo), ARLIM_P(volume_hi),
	Real* xflux_for_H, ARLIM_P(xflux_for_H_lo), ARLIM_P(xflux_for_H_hi),
	Real* xflux_for_Y, ARLIM_P(xflux_for_Y_lo), ARLIM_P(xflux_for_Y_hi),
	Real* yflux_for_Y, ARLIM_P(yflux_for_Y_lo), ARLIM_P(yflux_for_Y_hi),
	Real* yflux_for_H, ARLIM_P(yflux_for_H_lo), ARLIM_P(yflux_for_H_hi)
	);

    void FORT_RK_MULTICOMPONENT (
	const int* domain_lo, const int* domain_hi,
	const int* lo, const int* hi,
	Real* areax, ARLIM_P(areax_lo), ARLIM_P(areax_hi),
	Real* areay, ARLIM_P(areay_lo), ARLIM_P(areay_hi),
	const int* bc,
	const Real* dt,
	const Real* dx,
	int* index_of_firstY,
	int* index_of_lastY,
	int* index_of_rho,
	int* index_of_rhoH,
	int* index_of_T,
	Real* maximum_error,
	int* maximum_iterations,
	int* ncomps,
	int* nspcecies,
	Real* state, ARLIM_P(state_lo), ARLIM_P(state_hi),
	Real* update_for_H, ARLIM_P(update_for_H_lo), ARLIM_P(update_for_H_hi),
	Real* update_for_Y, ARLIM_P(update_for_Y_lo), ARLIM_P(update_for_Y_hi),
	Real* volume, ARLIM_P(volume_lo), ARLIM_P(volume_hi),
	Real* xflux_for_H, ARLIM_P(xflux_for_H_lo), ARLIM_P(xflux_for_H_hi),
	Real* xflux_for_Y, ARLIM_P(xflux_for_Y_lo), ARLIM_P(xflux_for_Y_hi),
	Real* yflux_for_Y, ARLIM_P(yflux_for_Y_lo), ARLIM_P(yflux_for_Y_hi),
	Real* yflux_for_H, ARLIM_P(yflux_for_H_lo), ARLIM_P(yflux_for_H_hi)
	);

    void FORT_RK_STEP_SELECTION (
	const int* lo, const int* hi,
	const Real* dx,
	int* index_of_firstY,
	int* index_of_lastY,
	int* index_of_rho,
	int* index_of_rhoH,
	int* index_of_T,
	Real* maximum_error,
	int* maximum_iterations,
	int* ncomps,
	int* nspcecies,
	Real* smallest_time_step,
	Real* state, ARLIM_P(state_lo), ARLIM_P(state_hi)
	);
}
#endif

#endif /*_HEATTRANSFER_H_*/
