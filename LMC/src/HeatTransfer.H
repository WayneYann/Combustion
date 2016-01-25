
#ifndef _HeatTransfer_H_
#define _HeatTransfer_H_

#include <AuxBoundaryData.H>
#include <NavierStokesBase.H>
#include <ChemDriver.H>

#include <list>
#include <map>
#include <utility>

//
// Note: define TEMPERATURE if you want a variable T in the 
//       State_Type part of the state
//       whose component index is Temp and 
//       whose evoution equation is 
//       \pd (rho T)/\pd t + diver (\rho U T) = \diver k/c_p grad T
//       define RADIATION only if TEMPERATURE is also defined and the evolution equation
//       for T is
//       \pd (rho T)/\pd t + diver (\rho U T) = \diver k/c_p grad T - 1/c_p diver q_rad
//
//       Note that component Temp is T, not rho T. This was done so
//       that we could use the existing diffusion operator and 
//       multigrid code.
//


class HeatTransfer
    :
    public NavierStokesBase
{

public:

    HeatTransfer ();
    HeatTransfer (Amr&            papa,
                  int             lev,
                  const Geometry& level_geom,
                  const BoxArray& bl,
                  Real            time);

    virtual ~HeatTransfer () BL_OVERRIDE;

    virtual void restart (Amr&          papa,
                          std::istream& is,
                          bool          bReadSpecial = false) BL_OVERRIDE;

    Real predict_velocity (Real dt, Real& comp_cfl);

    void checkTimeStep (Real dt);

    virtual void initData () BL_OVERRIDE;

    void initDataOtherTypes ();
    //
    // Reset time levels for the initial iterations.
    //
    void resetState (Real time,
		     Real dt_old,
		     Real dt_new);
    //
    // Set time levels of state data.
    //
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new) BL_OVERRIDE;
    //
    // Init data on this level from another NavierStokes (during regrid).
    //
    virtual void init (AmrLevel& old) BL_OVERRIDE;
    //
    // Init data on this level after regridding if old level
    // did not exist previously.
    //
    virtual void init () BL_OVERRIDE;
    //
    // Setup for a level timestep.
    //
    virtual void advance_setup (Real time,
				Real dt,
				int  iteration,
				int  ncycle) BL_OVERRIDE;
    //
    // Advance grids at this level in time.
    //
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) BL_OVERRIDE;
    //
    // Worker functions used by level advance function.
    //
    virtual void create_mac_rhs (MultiFab& mac_rhs, int nGrow, Real time, Real dt) BL_OVERRIDE;

    virtual void getViscTerms (MultiFab& visc_terms,
			       int       src_comp, 
			       int       num_comp,
			       Real      time) BL_OVERRIDE;

    void getTempViscTerms (MultiFab& visc_terms,
                           int       src_comp, 
                           Real      time);

    void getRhoHViscTerms (MultiFab& visc_terms,
                           int       src_comp, 
                           Real      time);

    void diffuse_scalar_setup (Real       dt,
			       int        sigma,
			       int&       rho_flag, 
			       MultiFab*& delta_rhs,
			       MultiFab*& alpha, 
			       FluxBoxes& fb_betan,
			       FluxBoxes& fb_betanp1);

    void diffuse_rhoh_setup (Real       time,
                             Real       dt,
                             MultiFab*& delta_rhs);

    void diffuse_spec_setup (int        istate,
                             Real       time,
                             Real       dt,
                             MultiFab*& delta_rhs);

    void diffuse_temp_setup (Real       time,
                             Real       dt,
                             MultiFab*& delta_rhs,
                             MultiFab*& alpha);

    void momentum_advection (Real dt,
			     bool do_adv_reflux);

    void scalar_advection (Real dt,
			   int  first_scalar,
			   int  last_scalar,
			   bool do_adv_reflux);
    
    void scalar_update (Real dt,
                        int  first_scalar,
                        int  last_scalar,
                        int  do_viscreflux);

    void differential_spec_diffusion_update (Real dt,
					     int  corrector);
    
    void adjust_spec_diffusion_update (MultiFab&              S_new,
				       const MultiFab*        S_old,
				       int                    sCompS,
				       Real                   dt,
				       Real                   time,
				       const Array<int>&      rho_flag,
				       const MultiFab&        rho_half,
				       int                    dataComp,
				       const MultiFab*        delta_rhs, 
				       const MultiFab*        alpha, 
				       const MultiFab* const* betanp1);
    
    void scalar_diffusion_update (Real dt,
				  int  first_scalar, 
				  int  last_scalar,
				  int  do_viscreflux);

    virtual void velocity_diffusion_update (Real dt) BL_OVERRIDE;

    void diffuse_velocity_setup (Real        dt,
				 MultiFab*&  delta_rhs,
				 FluxBoxes&  fb_betan, 
				 FluxBoxes&  fb_betanp1);
    //
    // Compute timesteps.
    //
    virtual void post_timestep (int iteration) BL_OVERRIDE;

    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep (Real cumtime) BL_OVERRIDE;

    void post_init_press (Real&        dt_init,
			  Array<int>&  nc_save,
			  Array<Real>& dt_save);

    virtual void post_regrid (int lbase,
                              int new_finest) BL_OVERRIDE;

    virtual void post_restart () BL_OVERRIDE;

    virtual void post_init (Real stop_time) BL_OVERRIDE;

    void init_once ();
    //
    // Static member functions.
    //
    static ChemDriver& getChemSolve () { return *chemSolve; }

    //
    // Returns a MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the MultiFab is built on the appropriately
    // grown BoxArray.
    //
    virtual MultiFab* derive (const std::string& name,
                              Real               time,
                              int                ngrow) BL_OVERRIDE;
    //
    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    //
    virtual void derive (const std::string& name,
                         Real               time,
                         MultiFab&          mf,
                         int                dcomp) BL_OVERRIDE;
    //
    // Define data descriptors.
    //
    static void variableSetUp ();

    static void ydotSetUp ();
    //
    // Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    
protected:

    //
    // Timestep estimation functions follow.
    //
    virtual Real estTimeStep () BL_OVERRIDE;

    void define_data ();

    //
    // Protected member functions.
    //
    virtual void sum_integrated_quantities () BL_OVERRIDE;

    virtual void mac_sync () BL_OVERRIDE;
    //
    // Crse/fine fixup functions.
    //
    virtual void reflux () BL_OVERRIDE;
    virtual void avgDown () BL_OVERRIDE; // Average down for all the state types.
    //
    // Note: these two functions must be supplied in NS_setup.C
    //       If divu (dsdt) is not included in the state, then
    //       they can be no-op functions
    //       See example in the _tube_flow problem dependent directory
    //
    virtual void calc_divu (Real      time,
			    Real      dt,
			    MultiFab& divu) BL_OVERRIDE;

    virtual void calc_dsdt (Real      time,
			    Real      dt,
			    MultiFab& dsdt) BL_OVERRIDE;

    void calc_dpdt (Real      time,
		    Real      dt,
		    MultiFab& dpdt,
		    MultiFab* u_mac);

    void compute_vel_visc (Real      time,
			   MultiFab* beta);

    void setThermoPress(Real time);

    void compute_rhoRT (const MultiFab& S,
                        MultiFab&       P,
                        int             pComp,
                        const MultiFab* T = 0);

    void temperature_stats (MultiFab& S);

    void compute_rhohmix (Real       time,
                          MultiFab& rhohmix);

    void compute_cp (Real      time,
                     MultiFab& cp);

    void compute_cp (const FArrayBox& temp, 
                     const FArrayBox& species,
                     FArrayBox&       cp);

    void compute_h (Real      time,
                    MultiFab& cp);

    void compute_rhoDgradYgradH (Real      time,
                                 MultiFab& rdgydgh);

    void set_rho_to_species_sum (MultiFab& S,
                                 int       strtcomp, 
                                 int       nghost_in = 0,
                                 int       minzero = 0);

    void set_rho_to_species_sum (MultiFab& S_in,
                                 int       s_in_start,
                                 MultiFab& S_out,
                                 int       s_out_start, 
                                 int       nghost_in = 0,
                                 int       minzero = 0);

    void scale_species (MultiFab& S,
                        int       strtcomp,
                        int       minzero = 0);

    void rhoh_update (Real time,
                      Real dt,
                      int  corrector);

    void tracer_update (Real dt,
                        int  corrector);

    void temp_update (Real dt,
                      int  corrector);

    void spec_update (Real time,
                      Real dt,
                      int  corrector);
    //
    // Protected static member functions.
    //
    static void Initialize ();
    static void Finalize ();

    static void RhoH_to_Temp (MultiFab& S, int nGrow = 0);

    static void RhoH_to_Temp (MultiFab& S,
                              MultiFab& temp,
                              int nGrow = 0,
                              int dominmax = 1);

    static void compute_cp_and_hmix (const MultiFab& S,
                                     MultiFab&       cp, 
                                     MultiFab&       hmix,
                                     MultiFab*       temp,
                                     int             nGrow,
                                     int             calchmix,
                                     int             floor_spec);
    //
    // advance functions
    //
    void compute_edge_states (Real dt,
			      std::vector<int>* state_comps_to_compute = 0);
    void differential_spec_diffuse_sync (Real dt);
    void compute_differential_diffusion_terms (MultiFab& visc_terms,
					       int       sComp,
					       Real      time);
    void adjust_spec_diffusion_fluxes (Real time);

    enum YdotAction { HT_EstimateYdotNew, HT_ImproveYdotOld, HT_LeaveYdotAlone };

    void strang_chem (MultiFab&  state,
		      Real       dt,
		      YdotAction Ydot_action,
		      int        ngrow = 0);

    virtual void calcViscosity (const Real time,
				const Real dt,
				const int  iteration,
				const int  ncycle) BL_OVERRIDE;

    virtual void calcDiffusivity (const Real time) BL_OVERRIDE;

    void calcDiffusivity (const Real time,
                          bool       doCalcVisc);

    virtual void getViscosity (MultiFab* viscosity[BL_SPACEDIM],
			       const Real time) BL_OVERRIDE;

    void getDiffusivity (MultiFab* diffusivity[BL_SPACEDIM],
			 const Real time,
			 const int state_comp,
			 const int dst_comp,
			 const int num_comp);
    
    void zeroBoundaryVisc (MultiFab*  beta[BL_SPACEDIM],
                           const Real time,
                           const int  state_comp,
                           const int  dst_comp,
                           const int  ncomp) const;

    void compute_OT_radloss (Real      time,
                             int       nGrow,
                             MultiFab& radloss);

    void reset_rho_in_rho_states (const MultiFab& rho,
                                  Real            time,
                                  const int       sComp,
                                  const int       nComp);

    virtual void set_preferred_boundary_values (MultiFab& S,
                                                int       state_index,
                                                int       src_comp,
                                                int       dst_comp,
                                                int       num_comp,
                                                Real      time) const BL_OVERRIDE;

    void set_overdetermined_boundary_cells (Real time);
    //
    // Set the variables that are put in the plotfile....
    //
    virtual void setPlotVariables () BL_OVERRIDE;
    //
    // Write plot file stuff to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                VisMF::How         how) BL_OVERRIDE;
    //
    // Functions for interpolating from cell centers to cell edges
    //
    enum FPLoc { HT_Edge = 0, HT_Center };

    static FPLoc fpi_phys_loc (int p_bc);

    static void center_to_edge_fancy (const FArrayBox& cfab,
                                      FArrayBox&       efab,
                                      const Box&       ccBox,
                                      int              sComp,
                                      int              dComp,
                                      int              nComp,
                                      const Box&       domain,
                                      const FPLoc&     bc_lo,
                                      const FPLoc&     bc_hi);

    void add_heat_sources(MultiFab& sum,
                          int       dComp,
                          Real      time,
                          int       nGrow,
                          Real      scale);

    enum Solver_Status {HT_InProgress, HT_Stalled, HT_Solved};

    void set_htt_hmixTYP ();

    void set_typical_values(bool restart);
    void reset_typical_values(const MultiFab& S);
    
    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os,
                             VisMF::How         how = VisMF::OneFilePerCPU,
                             bool               dump_old = true) BL_OVERRIDE;

#ifdef PARTICLES
    virtual int timestamp_num_extras () BL_OVERRIDE;
    virtual void timestamp_add_extras (int lev, Real time, MultiFab& mf) BL_OVERRIDE;
#endif

    DistributionMapping getFuncCountDM (const BoxArray& bxba, int ngrow = 0);

private:

    HeatTransfer& getLevel (int lev)
    {
        return *(HeatTransfer*) &parent->getLevel(lev);
    }

#ifdef PARTICLES
    static void read_particle_params ();
#endif

    MultiFab**  EdgeState;
    MultiFab**  SpecDiffusionFluxn;
    MultiFab**  SpecDiffusionFluxnp1;

    enum CurrentFluxType {HT_None, HT_Diffusion, HT_SyncDiffusion, HT_ExplicitDiffusion};

    Array<CurrentFluxType> spec_diffusion_flux_computed;
    AuxBoundaryData        aux_boundary_data_old;
    AuxBoundaryData        aux_boundary_data_new;
    bool                   FillPatchedOldState_ok;
    bool                   FillPatchedNewState_ok;

    static bool                     plot_reactions;
    static bool                     plot_consumption;
    static bool                     plot_heat_release;
    std::map<std::string,MultiFab*> auxDiag;
    static std::map<std::string,Array<std::string> > auxDiag_names;

    //
    // Protected static data.
    //
    static int         num_divu_iters;
    static int         init_once_done;
    static int         RhoH;
    static int         do_OT_radiation;
    static int         do_heat_sink;
    static int         unity_Le;
    static int         FuncCount_Type;
    static int         Ydot_Type;
    static int         divu_ceiling;
    static Real        min_rho_divu_ceiling;
    static Real        divu_dt_factor;
    static Real        htt_tempmin;
    static Real        htt_tempmax;
    static Real        htt_hmixTYP;
    static int         Trac;
    static int         have_trac;
    static int         RhoRT;
    static int         have_rhort;
    static int         first_spec;
    static int         last_spec;
    static int         nspecies;
    static int         floor_species;
    static ChemDriver* chemSolve;
    static Real        rgas;
    static Real        prandtl;
    static Real        schmidt;
    static Real        constant_mu_val;
    static Real        constant_rhoD_val;
    static Real        constant_lambda_val;
    static int         do_diffuse_sync;
    static int         do_reflux_visc;
    static int         dpdt_option;
    static int         do_set_rho_to_species_sum;
    static int         siegel_test;
    static int         zeroBndryVisc;
    static int         do_add_nonunityLe_corr_to_rhoh_adv_flux;
    static int         do_check_divudt;
    static int         hack_nochem;
    static int         hack_nospecdiff;
    static int         hack_noavgdivu;
    static Real        trac_diff_coef;
    static int         use_tranlib;
    static Real        P1atm_MKS;
    static std::string turbFile;
    static std::string fuelName;
    static Array<std::string> consumptionName;
    static int         calcDiffusivity_count;
    static Real        new_T_threshold;
    static Array<Real> typical_values;
    static int         reset_typical_vals_int;
    static std::map<std::string,Real> typical_values_FileVals;
};

#endif /*_HeatTransfer_H_*/
