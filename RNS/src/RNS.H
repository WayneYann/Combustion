#ifndef _RNS_H_
#define _RNS_H_

#define BEGIN_EXTERN_C extern "C" {
#define END_EXTERN_C   }

#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>
#include <ChemDriver.H>
#include <BL_CXX11.H>

#ifdef USE_SDCLIB
class SDCAmr;
#include <MLSDCAmr.H>
#endif

#include <iostream>

using std::istream;
using std::ostream;

enum StateType { State_Type = 0,
                 NUM_STATE_TYPE };

class RNS
    :
#ifdef USE_SDCLIB
    public MLSDCAmrLevel
#else
    public AmrLevel
#endif
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    RNS ();
    //
    //The basic constructor.
    //
    RNS (Amr&            papa,
	 int             lev,
	 const Geometry& level_geom,
	 const BoxArray& bl,
	 Real            time);
    //
    //The destructor.
    //
    virtual ~RNS () BL_OVERRIDE;
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (Amr&     papa,
			  istream& is,
			  bool bReadSpecial = false) BL_OVERRIDE;
    //
    //Call AmrLevel::checkPoint
    //
    virtual void checkPoint(const std::string& dir,
			    std::ostream&      os,
			    VisMF::How         how,
			    bool               dump_old) BL_OVERRIDE;

    /*A string written as the first item in writePlotFile() at
      level zero. It is so we can distinguish between different
      types of plot files. For RNS it has the form: RNS-Vnnn
    */
    virtual std::string thePlotFileType () const BL_OVERRIDE;

    virtual void setPlotVariables () BL_OVERRIDE;
    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
				ostream&       os,
				VisMF::How     how) BL_OVERRIDE;
    //
    //Define data descriptors.
    //
    static void variableSetUp ();
    //
    //Define tagging functions.
    //
    static void ErrorSetUp ();
    //
    //Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    //Initialize grid data at problem start-up.
    //
    virtual void initData () BL_OVERRIDE;

    //
    //Get problem metrics.
    //
    MultiFab* Area ();
    MultiFab& Area (int dir);
    //
    //The volume of the multifab.
    //
    MultiFab& Volume ();
    //
    //Initialize data on this level from another RNS (during regrid).
    //
    virtual void init (AmrLevel& old) BL_OVERRIDE;
    //
    //Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () BL_OVERRIDE;
    //
    //Proceed with next timestep?
    //
    virtual int okToContinue () BL_OVERRIDE;
    //
    //Advance grids at this level in time.
    //
    virtual Real advance (Real time,
			  Real dt,
			  int  iteration,
			  int  ncycle) BL_OVERRIDE;
    //
    // do chemistry
    // When it's called, 2 ghost cells of U should have already been filled properly.
    //
    void advance_chemistry(MultiFab& U, Real dt);
    void advance_chemistry(MultiFab& U, const MultiFab& Uguess, Real dt);
    //
    // do advection and diffusion
    // When it's called, Unew including ghost cells should have the state at t=time.
    //
    void advance_AD(MultiFab& Unew, Real time, Real dt, int iteration, int ncycle);
    //
    // fill boundary
    //
    void fill_boundary(MultiFab& U, Real time, int type,
		       bool isCorrection=false, bool isFEval=false);
    // choices for "type" in fill_boundary:
    static const int no_fill               = 0;
    static const int use_FillPatchIterator = 1;
    static const int use_FillCoarsePatch   = 2;
    static const int use_FillBoundary      = 3;
    static const int set_PhysBoundary      = 4;
    //
    // special fill boundary for RK4
    //
    void fill_rk_boundary(MultiFab& U, Real time, Real dt, int stage, int iteration, int ncycle);
    //
    // Given U, compute Uprime (i.e., dUdt) for AD
    //
    void dUdt_AD(MultiFab& U, MultiFab& Uprime, Real time, int fill_boundary_type,
		 FluxRegister* fine=0, FluxRegister* current=0, Real dt=0.0,
		 bool partialUpdate = false);
    //
    // Given U, compute Uprime (i.e., dUdt) for chemistry
    // When it's called, 2 ghost cells of U should have already been filled properly.
    //
    void dUdt_chemistry(const MultiFab& U, MultiFab& Uprime, bool partialUpdate = false);
    //
    // Compute U1 = U2 + c Uprime.
    //
    void update_rk(MultiFab& U1, const MultiFab& U2, Real c, const MultiFab& Uprime);
    //
    // Compute U = a Ua + b Ub + c Uprime.
    //
    void update_rk(MultiFab& U, Real a, const MultiFab& Ua, Real b, const MultiFab& Ub,
		   Real c, const MultiFab& Uprime);
    //
    // renormalize mass fractions and compute temperature
    //
    void post_update(MultiFab& U);
    //
    //Estimate time step.
    //
    Real estTimeStep (Real dt_old);
    //
    //Compute initial time step.
    //
    Real initialTimeStep ();
    //
    //Compute initial `dt'.
    //
    virtual void computeInitialDt (int                   finest_level,
				   int                   sub_cycle,
				   Array<int>&           n_cycle,
				   const Array<IntVect>& ref_ratio,
				   Array<Real>&          dt_level,
				   Real                  stop_time) BL_OVERRIDE;
    //
    //Compute new `dt'.
    //
    virtual void computeNewDt (int                   finest_level,
			       int                   sub_cycle,
			       Array<int>&           n_cycle,
			       const Array<IntVect>& ref_ratio,
			       Array<Real>&          dt_min,
			       Array<Real>&          dt_level,
			       Real                  stop_time,
			       int                   post_regrid_flag) BL_OVERRIDE;
    //
    //Do work after timestep().
    //
    virtual void post_timestep (int iteration) BL_OVERRIDE;
    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep (Real cumtime) BL_OVERRIDE;
    //
    //Do work after regrid().
    //
    virtual void post_regrid (int lbase,
			      int new_finest) BL_OVERRIDE;
    //
    //Do work after a restart().
    //
    virtual void post_restart () BL_OVERRIDE;
    //
    //Do work after init().
    //
    virtual void post_init (Real stop_time) BL_OVERRIDE;
    //
    //Error estimation for regridding.
    //
    virtual void errorEst (TagBoxArray& tb,
			   int          clearval,
			   int          tagval,
			   Real         time,
			   int          n_error_buf = 0,
			   int ngrow = 0) BL_OVERRIDE;

    // Returns a MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the MultiFab is built on the appropriately
    // grown BoxArray.
    virtual MultiFab* derive (const std::string& name,
			      Real               time,
			      int                ngrow) BL_OVERRIDE;

    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    virtual void derive (const std::string& name,
			 Real               time,
			 MultiFab&          mf,
			 int                dcomp) BL_OVERRIDE;

    static int numGrow();

    static int NUM_STATE;
    static int Density, Xmom, Ymom, Zmom, Eden, Temp;
    static int FirstSpec, LastSpec, NumSpec;
    static int verbose;
    static int do_chemistry;

    // species information
    static std::string  fuelName     ;
    static int          fuelID       ;
    static std::string  oxidizerName ;
    static int          oxidizerID   ;
    static std::string  productName  ;
    static int          productID    ;
    static std::string  flameTracName;
    static int          flameTracID  ;

    // how to compute dUdt in f2comp; "simple" means setting dUdt = \Delta U / \Delta t
    static int f2comp_simple_dUdt;
    static int f2comp_nbdf;  // only use bdf/vode for the first ? times on each node
    std::vector<int> f2comp_timer; // how many times f2comp has been called on each sdc node
    void reset_f2comp_timer(int nnodes) {
	f2comp_timer.clear();
	f2comp_timer.resize(nnodes-1,0); // last node doesn't call f2comp
    }

    void avgDown ();
    void avgDown (int state_indx);

    static class ChemDriver* chemSolve;

    // which grids touch fine?
    std::vector<int> touchFine;
    void buildTouchFine();
    void clearTouchFine() { touchFine.clear(); }

    MultiFab* chemstatus;
    void zeroChemStatus();

    static int check_imex_order(int ho_imex);

#ifdef USE_SDCLIB
    void restrict(MLSDCAmrEncap& F, MLSDCAmrEncap& G, Real t, MLSDCAmrLevel& levelG);
    void interpolate(MLSDCAmrEncap& F, MLSDCAmrEncap& G, Real t, bool correction, bool evaluation, MLSDCAmrLevel& levelG);
#endif

protected:

    static void read_params ();

    RNS& getLevel (int lev);

    FluxRegister& getFluxReg ();
    FluxRegister& getFluxReg (int lev);

    MultiFab& getRKk (int i);
    MultiFab& getRKk (int lev, int i);

    FluxRegister& getFluxRegRK ();
    FluxRegister& getFluxRegRK (int lev);

    void reflux ();

    void buildMetrics ();

    void sum_conserved_variables ();
    void volWgtSumCons(Array<Real>& s);

    void sum_chemstatus();

    //
    // The data.
    //
    MultiFab             volume;
    MultiFab             area[BL_SPACEDIM];
    FluxRegister*        flux_reg;

    // For RK
    MultiFab*        RK_k;
    FluxRegister*    flux_reg_RK;

    //
    // Static data members.
    //
    static bool      dump_old;
    static Real      cfl;
    static Real      diff_cfl;
    static Real      init_shrink;
    static Real      change_max;
    static int       do_reflux;
    static Array<int> lo_bc;
    static Array<int> hi_bc;
    static BCRec     phys_bc;
    static int       NUM_GROW;

    static ErrorList err_list;
    static int       allow_untagging;
    static int       do_density_ref;
    static int       do_temperature_ref;
    static int       do_pressure_ref;
    static int       do_velocity_ref;
    static int       do_vorticity_ref;
    static int       do_flametrac_ref;

    static Real      small_dens;
    static Real      small_temp;
    static Real      small_pres;
    static Real      gamma;

    static int       gravity_dir; // 1: x, 2: y, 3: z
    static Real      gravity;
    static Real      Treference;

    static int       RK_order;

    enum RiemannType { HLL = 0,
		       JBB,
		       HLLC,
		       NUM_RIEMANN_TYPE };

    static RiemannType Riemann;
    static Real        difmag; // for JBB & HLLC Riemann solvers
    static Real        HLL_factor;  // mixing some HLL into JBB or HLLC

    // what to store in plotfiles
    static int         plot_cons;
    static int         plot_prim;
    static int         plot_magvel;
    static int         plot_Mach;
    static int         plot_divu;
    static int         plot_magvort;
    static int         plot_X;
    static int         plot_omegadot;
    static int         plot_dYdt;
    static int         plot_heatRelease;
    static int         plot_fuelConsumption;
    static int         plot_primplus;

    static int icomp_cons, icomp_prim, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort,
	icomp_X, icomp_omegadot, icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption;
    static std::vector<std::string> plot_names;

    static std::string job_name;

    static std::vector<int> blocksize;

    static int         do_quartic_interp;

    static int         do_weno;
    static int         do_mdcd_weno;
    static int         weno_p;
    static Real        weno_eps;
    static Real        weno_gauss_phi;

    static int         use_vode;
    static int         new_J_cell;
    static int         chem_do_weno;

    enum ChemSolverType { CC_BURNING = 0, // 0: burn at cell centers
			  GAUSS_BURNING,  // 1: burn at Gauss points using BDF/VODE
			  SPLIT_BURNING,  // 2: burn at Gauss points by splitting
			  BECC_BURNING,   // 3: burn at cell center using backward Euler
			  BEGP_BURNING,   // 4: burn at Gauss points using backward Euler
			  NUM_CHEMSOLVER_TYPE };
    // For CC, GAUSS, and SPLIT, mlsdc.ho_imex should be 1,
    // whereas for BECC and BEGP, mlsdc.ho_imex could be either 0 or 1.
    static ChemSolverType chem_solver;

    // for keeping track of the amount of CPU time used -- this will persist
    // after restarts
    static Real        previousCPUTimeUsed;
    static Real        startCPUTime;

    Real getCPUTime();

    int num_ad_evals;
    int num_chem_evals;
};

//
// Inlines.
//

inline
int
RNS::numGrow()
{
    return NUM_GROW;
}

inline
MultiFab*
RNS::Area ()
{
    return area;
}

inline
MultiFab&
RNS::Area (int dir)
{
    return area[dir];
}

inline
MultiFab&
RNS::Volume ()
{
    return volume;
}

inline
RNS&
RNS::getLevel (int lev)
{
    return *(RNS *) &parent->getLevel(lev);
}

inline
FluxRegister&
RNS::getFluxReg ()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
FluxRegister&
RNS::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}

inline
MultiFab&
RNS::getRKk (int i)
{
    return RK_k[i];
}

inline
MultiFab&
RNS::getRKk (int lev, int i)
{
    return getLevel(lev).getRKk(i);
}

inline
FluxRegister&
RNS::getFluxRegRK ()
{
    BL_ASSERT(flux_reg_RK);
    return *flux_reg_RK;
}

inline
FluxRegister&
RNS::getFluxRegRK (int lev)
{
    return getLevel(lev).getFluxRegRK();
}


#endif /*_RNS_H_*/
