#ifndef _RNS_H_
#define _RNS_H_

#define BEGIN_EXTERN_C extern "C" {
#define END_EXTERN_C   }

#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>
#include <ChemDriver.H>

#ifdef USE_SDCLIB
class SDCAmr;
#endif

#include <iostream>

using std::istream;
using std::ostream;

enum StateType { State_Type = 0,
                 NUM_STATE_TYPE };

class RNS
    :
    public AmrLevel
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
    virtual ~RNS ();
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (Amr&     papa,
			  istream& is,
			  bool bReadSpecial = false);
    //
    //Call AmrLevel::checkPoint
    //
    virtual void checkPoint(const std::string& dir,
			    std::ostream&      os,
			    VisMF::How         how,
			    bool               dump_old);

    /*A string written as the first item in writePlotFile() at
      level zero. It is so we can distinguish between different
      types of plot files. For RNS it has the form: RNS-Vnnn
    */
    virtual std::string thePlotFileType () const;

    virtual void setPlotVariables ();
    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
				ostream&       os,
				VisMF::How     how);
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
    virtual void initData ();

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
    virtual void init (AmrLevel& old);
    //
    //Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init ();
    //
    //Proceed with next timestep?
    //
    virtual int okToContinue ();
    //
    //Advance grids at this level in time.
    //
    virtual Real advance (Real time,
			  Real dt,
			  int  iteration,
			  int  ncycle);
    //
    // do chemistry
    // When it's called, 2 ghost cells of U should have already been filled properly.
    //
    void advance_chemistry(MultiFab& U, Real dt);
    //
    // do advection and diffusion
    // When it's called, Unew including ghost cells should have the state at t=time.
    //
    void advance_AD(MultiFab& Unew, Real time, Real dt);
    //
    // fill boundary
    //
    void fill_boundary(MultiFab& U, Real time, int type);
    // choices for "type" in fill_boundary:
    static const int no_fill               = 0;
    static const int use_FillPatchIterator = 1;
    static const int use_FillCoarsePatch   = 2;
    static const int use_FillBoundary      = 3;
    static const int set_PhysBoundary      = 4;
    //
    // Given U, compute Uprime (i.e., dUdt) for AD
    //
    void dUdt_AD(MultiFab& U, MultiFab& Uprime, Real time, int fill_boundary_type,
	      FluxRegister* fine=0, FluxRegister* current=0, Real dt=0.0);
    //
    // Given U, compute Uprime (i.e., dUdt) for chemistry
    // When it's called, 2 ghost cells of U should have already been filled properly.
    //
    void dUdt_chemistry(const MultiFab& U, MultiFab& Uprime);
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
				   Real                  stop_time);
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
			       int                   post_regrid_flag);
    //
    //Do work after timestep().
    //
    virtual void post_timestep (int iteration);
    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep (Real cumtime);
    //
    //Do work after regrid().
    //
    virtual void post_regrid (int lbase,
			      int new_finest);
    //
    //Do work after a restart().
    //
    virtual void post_restart ();
    //
    //Do work after init().
    //
    virtual void post_init (Real stop_time);
    //
    //Error estimation for regridding.
    //
    virtual void errorEst (TagBoxArray& tb,
			   int          clearval,
			   int          tagval,
			   Real         time,
			   int          n_error_buf = 0, int ngrow = 0);

    // Returns a MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the MultiFab is built on the appropriately
    // grown BoxArray.
    MultiFab* derive (const std::string& name,
		      Real               time,
		      int                ngrow);

    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    void derive (const std::string& name,
		 Real               time,
		 MultiFab&          mf,
		 int                dcomp);

    static int numGrow();

    static int NUM_STATE;
    static int Density, Xmom, Ymom, Zmom, Eden, Temp;
    static int FirstSpec, LastSpec, NumSpec;
    static int verbose;

    void avgDown ();
    void avgDown (int state_indx);
    void avgDown (MultiFab& S_crse, const MultiFab& S_fine);

    static class ChemDriver* chemSolve;

#ifdef USE_SDCLIB
    SDCAmr* getSDCAmr() { return (SDCAmr*) parent; }
#endif


protected:

    static void read_params ();

    RNS& getLevel (int lev);

    FluxRegister& getFluxReg ();
    FluxRegister& getFluxReg (int lev);

    void reflux ();


    void buildMetrics ();

    void sum_conserved_variables ();
    void volWgtSumCons(Array<Real>& s);

    //
    // The data.
    //
    MultiFab             volume;
    MultiFab             area[BL_SPACEDIM];
    FluxRegister*        flux_reg;
    //
    // Static data members.
    //
    static bool      dump_old;
    static Real      cfl;
    static Real      init_shrink;
    static Real      change_max;
    static int       do_reflux;
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

    static Real      gravity;
    static Real      Treference;

    static int       RK_order;

    enum RiemannType { HLL = 0,
		       JBB,
		       HLLC,
		       NUM_RIEMANN_TYPE };

    static RiemannType Riemann;
    static Real        difmag; // for JBB & HLLC Riemann solvers

    // species information
    static std::string  fuelName     ;
    static int          fuelID       ;
    static std::string  oxidizerName ;
    static int          oxidizerID   ;
    static std::string  productName  ;
    static int          productID    ;
    static std::string  flameTracName;
    static int          flameTracID  ;

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

    static int         do_component_weno;

    static int         use_vode;
    static int         do_cc_burning;
    static int         split_burning;

    // for keeping track of the amount of CPU time used -- this will persist
    // after restarts
    static Real        previousCPUTimeUsed;
    static Real        startCPUTime;

    Real getCPUTime();
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

#endif /*_RNS_H_*/
