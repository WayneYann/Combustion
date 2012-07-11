
#ifndef _HeatTransfer_H_
#define _HeatTransfer_H_

#include <AuxBoundaryData.H>
#include <NavierStokes.H>
#include <ChemDriver.H>
#include <DDOp.H>

#ifdef PARTICLES
#include <Particles.H>
//
// The "BL_SPACEDIM" bytes in m_data[] are to hold temporary positions and velocities in AdvectWithUmac().
//
typedef ParticleContainer<2*BL_SPACEDIM> HTParticleContainer;
#endif

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

#if 1
    virtual Real predict_velocity (Real  dt,
                                   Real& comp_cfl);
#endif
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
    // Advance grids at this level in time.
    //
    Real advance_sdc (Real time,
		      Real dt,
		      int  iteration,
		      int  ncycle);
    //
    // Worker functions used by level advance function.
    //
    virtual void create_mac_rhs (MultiFab& rhs, Real time, Real dt, int nGrow = 0);

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

    void diffuse_cleanup (MultiFab*&  delta_rhs, 
                          MultiFab**& betan,
                          MultiFab**& betanp1,
                          MultiFab*&  alpha);

    void diffuse_cleanup (MultiFab*&  delta_rhs, 
                          MultiFab**& betan,
                          MultiFab**& betanp1);

    virtual void scalar_advection_update (Real dt,
                                          int  first_scalar,
                                          int  last_scalar);

    void compute_scalar_advection_fluxes_and_divergence (MultiFab& Force,
                                                         Real      dt);

    void differential_diffusion_update (MultiFab& Force,
                                        int       FComp,
                                        Real      theta,
                                        MultiFab& D,
                                        int       DComp,
                                        MultiFab& DD,
                                        Real      theta_enthalpy);

    void flux_divergence (MultiFab&        fdiv,
                          int              fdivComp,
                          const MultiFab* const* extensive_fluxes,
                          int              fluxComp,
                          int              nComp,
                          Real             scale) const;

    void compute_differential_diffusion_fluxes (const Real& time);

    void set_differential_diffusion_fluxes (const Real& time);

    void adjust_spec_diffusion_fluxes (Real                   time,
                                       const MultiFab* const* beta,
                                       bool                   grow_cells_already_filled);

    void compute_enthalpy_fluxes (Real                   time,
                                  const MultiFab* const* beta,
                                  bool                   grow_cells_already_filled);

    virtual void velocity_diffusion_update (Real dt);

    virtual void diffuse_velocity_setup (Real        dt,
                                         MultiFab*&  delta_rhs,
                                         MultiFab**& betan, 
                                         MultiFab**& betanp1);


    enum HowToFillGrow {HT_ZERO_GROW_CELLS, HT_EXTRAP_GROW_CELLS, HT_NUM_GROW_OPTIONS};

    void compute_instantaneous_reaction_rates (MultiFab&       R,
					       const MultiFab& S,
					       int             nGrow = 0,
					       HowToFillGrow   how = HT_ZERO_GROW_CELLS);

    //
    // Compute timesteps.
    //
    virtual void post_timestep (int iteration);

    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep (Real cumtime);

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

#ifdef PARTICLES
    static HTParticleContainer* theHTPC ();

    MultiFab* ParticleDerive (const std::string& name, Real time, int ngrow);
#endif
    //
    // Returns a MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the MultiFab is built on the appropriately
    // grown BoxArray.
    //
    MultiFab* derive (const std::string& name,
                      Real               time,
                      int                ngrow);
    //
    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    //
    void derive (const std::string& name,
                 Real               time,
                 MultiFab&          mf,
                 int                dcomp);
    //
    // Define data descriptors.
    //
    static void variableSetUp ();

    static void rhoydotSetUp ();
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

    //
    // Protected static member functions.
    //
    static void Initialize ();
    static void Finalize ();

    static void RhoH_to_Temp (MultiFab& S,
                              int       nGrow = 0,
                              int       dominmax = false);

    static int RhoH_to_Temp (FArrayBox& S,
                             const Box& box,
                             int        dominmax = false);

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
                                      std::vector<int>* state_comps_to_compute = 0);
    virtual void differential_spec_diffuse_sync (Real dt);
    virtual void compute_differential_diffusion_terms (MultiFab& D,
                                                       MultiFab& DD,
						       Real      time);

    enum YdotAction { HT_EstimateYdotNew, HT_ImproveYdotOld, HT_LeaveYdotAlone };

    virtual void strang_chem (MultiFab&  state,
			      Real       dt,
			      YdotAction Ydot_action,
                              int        ngrow = 0);

    virtual void advance_chemistry (MultiFab&       mf_old,
                                    MultiFab&       mf_new,
                                    Real            dt,
                                    const MultiFab& Force,
                                    int             nCompF);

    virtual void calcViscosity (const Real time);

    virtual void calcDiffusivity (const Real time);

    void calcDiffusivity (const Real time,
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

    void compute_mcdd_visc_terms(MultiFab&           vtermsYH,
                                 Real                time,
                                 int                 nGrow,
                                 DDOp::DD_ApForTorRH whichApp,
                                 PArray<MultiFab>*   flux = 0);

    enum Solver_Status {HT_InProgress, HT_Stalled, HT_Solved};
    struct MCDD_MGParams {
        MCDD_MGParams() {}
        MCDD_MGParams(int nComp_)
            {
                nComp = nComp_;
                maxRes.resize(nComp);
                maxRes_initial.resize(nComp,0);
                maxCor.resize(nComp);
                status = HT_InProgress;
            }
        MCDD_MGParams(const MCDD_MGParams& rhs)
            {
                nComp = rhs.nComp;
                maxRes.resize(nComp);
                maxRes_initial.resize(nComp,0);
                maxCor.resize(nComp);
                for (int i=0; i<nComp; ++i) {
                    maxRes[i] = rhs.maxRes[i];
                    maxCor[i] = rhs.maxCor[i];
                    maxRes_initial[i] = rhs.maxRes_initial[i];
                }
                maxRes_norm = rhs.maxRes_norm;
                maxCor_norm = rhs.maxCor_norm;
                res_nu1_rtol = rhs.res_nu1_rtol;
                res_nu2_rtol = rhs.res_nu2_rtol;
                res_redux_tol = rhs.res_redux_tol;
                res_abs_tol = rhs.res_abs_tol;
                stalled_tol = rhs.stalled_tol;
                gamma = rhs.gamma;
                num_coarser = rhs.num_coarser;
                nu1 = rhs.nu1;
                nu2 = rhs.nu2;
                mg_level = rhs.mg_level;
                status = rhs.status;
                iter = rhs.iter;
            }
        Array<Real> maxRes;
        Array<Real> maxRes_initial;
        Array<Real> maxCor;
        Real        maxRes_norm;
        Real        maxCor_norm;
        Real        res_nu1_rtol;
        Real        res_nu2_rtol;
        Real        res_redux_tol;
        Real        res_abs_tol;
        Real        stalled_tol;
        int         iter;
        int         nComp;
        Solver_Status status;
        int         num_coarser;
        int         gamma;
        int         nu1;
        int         nu2;
        int         mg_level;
    };

    void mcdd_fas_cycle(MCDD_MGParams&      status,
                        MultiFab&           S,
                        const MultiFab&     Rhs,
                        const MultiFab&     rhoCpInv,
                        PArray<MultiFab>&   flux,
                        Real                time,
                        Real                dt,
                        DDOp::DD_ApForTorRH whichApp);
    
    void mcdd_update(Real time,
                     Real dt);

    void mcdd_diffuse_sync(Real dt);

    void set_htt_hmixTYP ();

    void set_typical_values(bool restart);
    
    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os,
                             VisMF::How         how = VisMF::OneFilePerCPU,
                             bool               dump_old = true);

protected:

    DistributionMapping getFuncCountDM (const BoxArray& bxba, int ngrow = 0);

    MultiFab**  EdgeState;
    MultiFab**  SpecDiffusionFluxn;
    MultiFab**  SpecDiffusionFluxnp1;
    MultiFab    sumSpecFluxDotGradHn;
    MultiFab    sumSpecFluxDotGradHnp1;

    AuxBoundaryData        aux_boundary_data_old;
    AuxBoundaryData        aux_boundary_data_new;
    bool                   FillPatchedOldState_ok;

    static bool                     plot_reactions;
    static bool                     plot_consumption;
    static bool                     plot_heat_release;
    std::map<std::string,MultiFab*> auxDiag;
    static std::map<std::string,Array<std::string> > auxDiag_names;

    DDOp MCDDOp;

    static int         num_divu_iters;
    static int         init_once_done;
    static int         RhoH;
    static int         do_OT_radiation;
    static int         do_heat_sink;
    static int         unity_Le;
    static int         FuncCount_Type;
    static int         RhoYdot_Type;
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
    static int         zeroBndryVisc;
    static int         do_add_nonunityLe_corr_to_rhoh_adv_flux;
    static int         do_check_divudt;
    static int         hack_nochem;
    static int         hack_nospecdiff;
    static int         hack_nomcddsync;
    static int         hack_noavgdivu;
    static Real        trac_diff_coef;
    static Real        P1atm_MKS;
    static std::string turbFile;
    static std::string fuelName;
    static Array<std::string> consumptionName;
    static int         do_mcdd;
    static std::string mcdd_transport_model;
    static int         mcdd_NitersMAX;
    static Real        mcdd_relax_factor_T;
    static Real        mcdd_relax_factor_Y;
    static int         mcdd_mgLevelsMAX;
    static Array<int>  mcdd_nu1;
    static Array<int>  mcdd_nu2;
    static int         mcdd_nub;
    static int         mcdd_numcycles;
    static Real        mcdd_res_nu1_rtol;
    static Real        mcdd_res_nu2_rtol;
    static Real        mcdd_res_redux_tol;
    static Real        mcdd_res_abs_tol;
    static Real        mcdd_stalled_tol;
    static Real        mcdd_advance_temp;
    static int         mcdd_verbose;
    static int         calcDiffusivity_count;
    static Real        new_T_threshold;
    static Array<Real> typical_values;

    static int sdc_iterMAX;
    static int nGrowAdvForcing;
    MultiFab   Forcing;
    MultiFab   Dn;
    MultiFab   Dnp1;
    MultiFab   Dhat;
    MultiFab   DDn;
    MultiFab   DDnp1;
};

#endif /*_HeatTransfer_H_*/
