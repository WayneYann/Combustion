
#ifndef _HeatTransfer_H_
#define _HeatTransfer_H_

#include <AuxBoundaryData.H>
#include <NavierStokes.H>
#include <ChemDriver.H>

#ifdef PARTICLES
#include <AmrParticles.H>
//
// The "BL_SPACEDIM" bytes in m_data[] are to hold temporary positions and velocities in AdvectWithUmac().
//
typedef AmrParticleContainer<2*BL_SPACEDIM> HTParticleContainer;
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
    // Worker functions used by level advance function.
    //
    virtual void create_mac_rhs (MultiFab& rhs, int nGrow, Real time, Real dt);

    virtual void getViscTerms (MultiFab& visc_terms,
                               int       src_comp, 
                               int       num_comp,
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

    void compute_scalar_advection_fluxes_and_divergence (const MultiFab& Force,
                                                         const MultiFab& divu,
                                                         Real            dt);

    void differential_diffusion_update (MultiFab& Force,
                                        int       FComp,
                                        MultiFab& D,
                                        int       DComp,
                                        MultiFab& DD);

    void flux_divergence (MultiFab&        fdiv,
                          int              fdivComp,
                          const MultiFab* const* extensive_fluxes,
                          int              fluxComp,
                          int              nComp,
                          Real             scale) const;

    void compute_differential_diffusion_fluxes (const Real& time,
					        const Real& dt);

    void adjust_spec_diffusion_fluxes (Real time);

    void compute_enthalpy_fluxes (Real                   time,
                                  const MultiFab* const* beta);

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

    void ParticleDerive (const std::string& name, Real time, MultiFab& mf, int dcomp);
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
                            MultiFab& divu,
			    bool      is_divu_iter = false);

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

    virtual void temperature_stats (MultiFab& S);

    void compute_rhohmix (Real       time,
                          MultiFab& rhohmix);

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
    //
    // advance functions
    //
    virtual void set_reasonable_grow_cells_for_R(Real time);
    virtual void differential_spec_diffuse_sync (Real dt,
					         bool Wbar_corrector);

    virtual void compute_differential_diffusion_terms (MultiFab& D,
                                                       MultiFab& DD,
#ifdef USE_WBAR
						       MultiFab& DWbar,
#endif
						       Real      time,
						       Real      dt);

#ifdef USE_WBAR
    void compute_Wbar_fluxes(Real time,
                             Real inc);
#endif

    enum YdotAction { HT_EstimateYdotNew, HT_ImproveYdotOld, HT_LeaveYdotAlone };

    virtual void advance_chemistry (MultiFab&       mf_old,
                                    MultiFab&       mf_new,
                                    Real            dt,
                                    const MultiFab& Force,
                                    int             nCompF,
                                    bool            use_stiff_solver = true);

    virtual void calcViscosity (const Real time,
                                const Real dt,
                                const int  iteration,
                                const int  ncycle);

    virtual void calcDiffusivity (const Real time);

#ifdef USE_WBAR
    void calcDiffusivity_Wbar (const Real time);
#endif

    virtual void getViscosity (MultiFab* viscosity[BL_SPACEDIM],
                               const Real time);

    virtual void getDiffusivity (MultiFab* diffusivity[BL_SPACEDIM],
                                 const Real time,
                                 const int state_comp,
                                 const int dst_comp,
                                 const int num_comp);

#ifdef USE_WBAR
    void getDiffusivity_Wbar (MultiFab* diffusivity[BL_SPACEDIM],
			      const Real time);
#endif

    void zeroBoundaryVisc (MultiFab*  beta[BL_SPACEDIM],
                           const Real time,
                           const int  state_comp,
                           const int  dst_comp,
                           const int  ncomp) const;

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

    DistributionMapping getFuncCountDM (const BoxArray& bxba, int ngrow);

    enum Solver_Status {HT_InProgress, HT_Stalled, HT_Solved};

    void set_htt_hmixTYP ();

    void set_typical_values(bool restart);
    void reset_typical_values(const MultiFab& S);
    
    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os,
                             VisMF::How         how = VisMF::OneFilePerCPU,
                             bool               dump_old = true);

protected:

    MultiFab**  EdgeState;
    MultiFab**  EdgeFlux;
    MultiFab**  SpecDiffusionFluxn;
    MultiFab**  SpecDiffusionFluxnp1;
#ifdef USE_WBAR
    MultiFab**  SpecDiffusionFluxWbar;
#endif
    MultiFab    sumSpecFluxDotGradHn;
    MultiFab    sumSpecFluxDotGradHnp1;

#ifdef USE_WBAR
    MultiFab *diffWbar_cc;
#endif

    static bool                     plot_reactions;
    static bool                     plot_consumption;
    static bool                     plot_heat_release;
    std::map<std::string,MultiFab*> auxDiag;
    static std::map<std::string,Array<std::string> > auxDiag_names;

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
    static bool        avg_down_chem;
    static int         reset_typical_vals_int;
    static std::map<std::string,Real> typical_values_FileVals;
    static int         sdc_iterMAX;
    
    static int nGrowAdvForcing;
    bool       updateFluxReg;
    bool       is_predictor;
};

#endif /*_HeatTransfer_H_*/
