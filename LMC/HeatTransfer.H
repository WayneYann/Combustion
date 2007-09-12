
#ifndef _HeatTransfer_H_
#define _HeatTransfer_H_

#include <AuxBoundaryData.H>
#include <NavierStokes.H>
#include <ChemDriver.H>
#include <DDOp.H>

#include <list>
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
    public NavierStokes
{

public:

    HeatTransfer ();
    HeatTransfer (Amr&            papa,
                  int             lev,
                  const Geometry& level_geom,
                  const BoxArray& bl,
                  Real            time);

    virtual ~HeatTransfer ();

    virtual void restart (Amr&          papa,
                          std::istream& is,
                          bool          bReadSpecial = false);
    //
    // Timestep estimation functions follow.
    //
    virtual Real estTimeStep ();

    void checkTimeStep (Real dt);

    virtual void initData ();

    virtual void initDataOtherTypes ();
    //
    // Reset time levels for the initial iterations.
    //
    virtual void resetState (Real time,
                             Real dt_old,
                             Real dt_new);
    //
    // Set time levels of state data.
    //
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new);
    //
    // Init data on this level from another NavierStokes (during regrid).
    //
    virtual void init (AmrLevel& old);
    //
    // Init data on this level after regridding if old level
    // did not exist previously.
    //
    virtual void init ();
    //
    // Setup for a level timestep.
    //
    virtual void advance_setup (Real time,
                                Real dt,
                                int  iteration,
                                int  ncycle);
    //
    // Advance grids at this level in time.
    //
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle);
    //
    // Worker functions used by level advance function.
    //
    virtual MultiFab* create_mac_rhs (Real time, Real dt);

    virtual void getViscTerms (MultiFab& visc_terms,
                               int       src_comp, 
                               int       num_comp,
                               Real      time);

    void getTempViscTerms (MultiFab& visc_terms,
                           int       src_comp, 
                           Real      time);

    void getRhoHViscTerms (MultiFab& visc_terms,
                           int       src_comp, 
                           Real      time);

    virtual void diffuse_scalar_setup (Real        dt,
                                       int         sigma,
                                       int*        rho_flag, 
                                       MultiFab*&  delta_rhs,
                                       MultiFab*&  alpha, 
                                       MultiFab**& betan,
                                       MultiFab**& betanp1);

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

    void diffuse_cleanup (MultiFab*&  delta_rhs, 
                          MultiFab**& betan,
                          MultiFab**& betanp1,
                          MultiFab*&  alpha);

    void diffuse_cleanup (MultiFab*&  delta_rhs, 
                          MultiFab**& betan,
                          MultiFab**& betanp1);

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
				       const MultiFab*        rho_half,
				       int                    dataComp,
				       const MultiFab*        delta_rhs, 
				       const MultiFab*        alpha, 
				       const MultiFab* const* betanp1);
    
    virtual void scalar_diffusion_update (Real dt,
                                          int  first_scalar, 
                                          int  last_scalar,
                                          int  do_viscreflux);

    virtual void velocity_diffusion_update (Real dt);

    virtual void diffuse_velocity_setup (Real        dt,
                                         MultiFab*&  delta_rhs,
                                         MultiFab**& betan, 
                                         MultiFab**& betanp1);

    //
    // Runge-Kutta diffusion functions used by level advance function.
    //
    void rk_diffusion_operator (const Real time,
				const Real scale_factor,
				MultiFab *& div_of_flux_for_H,
				MultiFab *& div_of_flux_for_Y,
				MultiFab **& flux_for_H,
				MultiFab **& flux_for_Y);

    Real rk_step_selection (const Real time);

    //
    // Make rho at time n.
    //
    virtual void make_rho_prev_time ();
    //
    // Make rho at time n+1.
    //
    virtual void make_rho_curr_time ();
    //
    // Compute timesteps.
    //
    virtual void post_timestep (int iteration);

    virtual void post_init_press (Real&        dt_init,
                                  Array<int>&  nc_save,
                                  Array<Real>& dt_save);

    virtual void post_regrid (int lbase,
                              int new_finest);

    virtual void post_restart ();

    virtual void post_init (Real stop_time);

    void init_once ();
    //
    // Static member functions.
    //
    static ChemDriver& getChemSolve () { return *chemSolve; }
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
    // Protected member functions.
    //
    virtual void sum_integrated_quantities ();

    HeatTransfer& getLevel (int lev)
    {
        return *(HeatTransfer*) &parent->getLevel(lev);
    }

    virtual void mac_sync ();
    //
    // Crse/fine fixup functions.
    //
    virtual void reflux ();
    virtual void avgDown (); // Average down for all the state types.
    //
    // Note: these two functions must be supplied in NS_setup.C
    //       If divu (dsdt) is not included in the state, then
    //       they can be no-op functions
    //       See example in the _tube_flow problem dependent directory
    //
    virtual void calc_divu (Real      time,
                            Real      dt,
                            MultiFab& divu);

    virtual void calc_dsdt (Real      time,
                            Real      dt,
                            MultiFab& dsdt);

    virtual void calc_dpdt (Real      time,
                            Real      dt,
                            MultiFab& dpdt,
                            MultiFab* u_mac);

    virtual void compute_vel_visc (Real      time,
                                   MultiFab* beta);

    void setThermoPress(Real time);

    void compute_rhoRT (const MultiFab& S,
                        MultiFab&       P,
                        int             pComp,
                        const MultiFab* T = 0);

    void compute_rhoh_new (Real time,
                           Real dt);

    virtual void temperature_stats (MultiFab& S);

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
    static void read_params ();

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
    virtual void compute_edge_states (Real dt,
                                      std::vector<bool>* state_comps_to_compute = 0);
    virtual void differential_spec_diffuse_sync (Real dt);
    virtual void compute_differential_diffusion_terms (MultiFab& visc_terms,
						       int       sComp,
						       Real      time);

    enum YdotAction { HT_EstimateYdotNew, HT_ImproveYdotOld, HT_LeaveYdotAlone };

    virtual void strang_chem (MultiFab&  state,
			      Real       dt,
			      YdotAction Ydot_action,
                              int        ngrow = 0);

    virtual void calcViscosity (const Real time,
                                const Real dt,
                                const int  iteration,
                                const int  ncycle);

    virtual void calcDiffusivity (const Real time,
                                  const Real dt,
                                  const int  iteration,
                                  const int  ncycle,
                                  const int  src_comp = 0,
                                  const int  num_comp = 1);

    void calcDiffusivity (const Real time,
                          const Real dt,
                          const int  iteration,
                          const int  ncycle,
                          const int  src_comp,
                          const int  num_comp,
                          bool       doCalcVisc);

    virtual void getViscosity (MultiFab* viscosity[BL_SPACEDIM],
                               const Real time);

    virtual void getDiffusivity (MultiFab* diffusivity[BL_SPACEDIM],
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
                                                Real      time) const;

    void set_overdetermined_boundary_cells (Real time);
    //
    // Set the variables that are put in the plotfile....
    //
    virtual void setPlotVariables ();
    //
    // Write plot file stuff to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                VisMF::How         how);
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

    //
    // Functions related to full-transport model
    //
    void fill_mcdd_boundary_data(Real time);

    void compute_mcdd_visc_terms(MultiFab&           vtermsY,
                                 int                 dCompY,
                                 MultiFab&           vtermsH,
                                 int                 dCompH,
                                 Real                time,
                                 int                 nGrow,
                                 DDOp::DD_ApForTorRH whichApp,
                                 PArray<MultiFab>*   flux = 0,
                                 int                 dCompFY = 0,
                                 int                 dCompFH = 0);

    void compute_mcdd_relax_parameter(MultiFab&       lambda,
                                      int             dCompY,
                                      int             dCompH,
                                      const MultiFab& Y,
                                      int             sCompY,
                                      const MultiFab& T,
                                      int             sCompT,
                                      const MultiFab& Rho,
                                      int             sCompR,
                                      Real            dt,
                                      int             mgLevel = 0);

    void mcdd_average_state(const MultiFab& SF,
                            int             sCompY,
                            int             sCompH,
                            const MultiFab& rhoF,
                            int             sCompFR,
                            MultiFab&       SC,
                            int             dCompY,
                            int             dCompH,
                            const MultiFab& rhoC,
                            int             sCompCR) const;

    void mcdd_interp_state(const MultiFab& SC,
                           int             sCompY,
                           int             sCompH,
                           const MultiFab& rhoC,
                           int             sCompCR,
                           MultiFab&       SF,
                           int             dCompY,
                           int             dCompH,
                           const MultiFab& rhoF,
                           int             sCompFR) const;

    void mcdd_v_cycle(MultiFab&         S,
                      int               sCompYS,
                      int               sCompHS,
                      const MultiFab&   Rhs,
                      int               sCompYR,
                      int               sCompHR,
                      MultiFab&         T,
                      int               sCompT,
                      const MultiFab&   Rho,
                      int               sCompR,
                      PArray<MultiFab>& flux,
                      int               dCompFY,
                      int               dCompFH,
                      Real              time,
                      Real              dt,
                      int               level = 0);

    void mcdd_update(Real time,
                     Real dt);

    void mcdd_apply(MultiFab& AH, int dCompH, MultiFab& AY, int dCompY,
                    const MultiFab& H, int sCompH, const MultiFab& Y, int sCompY,
                    MultiFab& T, int sCompT, const MultiFab& Rho, int sCompR,
                    PArray<MultiFab>& fluxH, int dCompFH,
                    PArray<MultiFab>& fluxY, int dCompFY,
                    Real time, Real dt, int level);

    void mcdd_residual(MultiFab& ResH, int dCompH, MultiFab& ResY, int dCompY,
                       const MultiFab& H, int sCompH, const MultiFab& Y, int sCompY,
                       MultiFab& T, int sCompT, const MultiFab& Rho, int sCompR,
                       const MultiFab& RhsH, int sCompRH,
                       const MultiFab& RhsY, int sCompRY,
                       PArray<MultiFab>& fluxH, int dCompFH,
                       PArray<MultiFab>& fluxY, int dCompFY,
                       Real time, Real dt, int level = 0);

    void mcdd_diffuse_sync(Real dt);
    
protected:

    DistributionMapping getFuncCountDM (const BoxArray& bxba);
    //
    // Protected data.
    //
    MultiFab*   Spec_adv_sync;    // Species advective sync update storage
    MultiFab*   Spec_visc_sync;   // Species diffusive sync update storage
    MultiFab*   RhoH_adv_sync;    // RhoH advective sync update storage
    MultiFab*   RhoH_visc_sync;   // RhoH diffusive sync update storage
    MultiFab**  EdgeState;
    MultiFab**  SpecDiffusionFluxn;
    MultiFab**  SpecDiffusionFluxnp1;

    enum CurrentFluxType {HT_None, HT_Diffusion, HT_SyncDiffusion, HT_ExplicitDiffusion};

    Array<CurrentFluxType> spec_diffusion_flux_computed;
    AuxBoundaryData        aux_boundary_data_old;
    AuxBoundaryData        aux_boundary_data_new;
    bool                   FillPatchedOldState_ok;
    PArray<MultiFab>       derive_mfs;

    static bool               plot_auxDiags;
    MultiFab*                 auxDiag;
    static Array<std::string> auxDiag_names;

    DDOp MCDDOp;
    //
    // Protected static data.
    //
    static int         num_divu_iters;
    static int         init_once_done;
    static int         RhoH;
    static int         do_OT_radiation;
    static int         do_heat_sink;
    static ChemDriver::Chem_Evolve chem_integrator;
    static int         unity_Le;
    static int         FuncCount_Type;
    static int         Ydot_Type;
    static int         divu_ceiling;
    static Real        min_rho_divu_ceiling;
    static Real        divu_dt_factor;
    static Real        htt_tempmin;
    static Real        htt_tempmax;
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
    static Real        P1atm_MKS;
    static std::string turbFile;
    static std::string fuelName;
    static std::string consumptionName;
    static int         do_mcdd;
    static int         mcdd_NitersMAX;
    static Real        mcdd_relaxFactor;
    static Real        mcdd_relaxFactor1;
    static Real        mcdd_rtol;
    static int         mcdd_presmooth;
    static int         mcdd_postsmooth;
    static Real        mcdd_cfRelaxFactor;
    static bool        do_rk_diffusion;
    static bool        rk_mixture_averaged;
    static Real        rk_time_step_multiplier;
    static int         calcDiffusivity_count;
};

#endif /*_HeatTransfer_H_*/
