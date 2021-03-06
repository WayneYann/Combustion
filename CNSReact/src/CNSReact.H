#ifndef _CNSReact_H_
#define _CNSReact_H_

#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>
#include <ChemDriver.H>

#include <iostream>

using std::istream;
using std::ostream;

enum StateType { State_Type = 0,
                 NUM_STATE_TYPE };

class CNSReact
    :
    public AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    CNSReact ();
    //
    //The basic constructor.
    //
    CNSReact (Amr&            papa,
	      int             lev,
	      const Geometry& level_geom,
	      const BoxArray& bl,
	      Real            time);
    //
    //The destructor.
    //
    virtual ~CNSReact ();
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
               types of plot files. For CNSReact it has the form: CNSReact-Vnnn
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
    //Set time levels of state data.
    //
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new);
    //
    //Initialize data on this level from another CNSReact (during regrid).
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

    Real advance_hydro (Real time,
                        Real dt,
                        int  iteration,
                        int  ncycle);
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
    //Allocate data at old time.
    //
    virtual void allocOldData ();
    //
    //Remove data at old time.
    //
    virtual void removeOldData ();
    //
    //Print information about energy budget.
    //
    void do_energy_diagnostics();
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

    void reset_internal_energy (MultiFab& State);

    void computeTemp (MultiFab& State);

    Real volWgtSum (const std::string& name, Real time);

    void set_special_tagging_flag (Real time);

    static int       NUM_STATE;
    static int       Density, Xmom, Ymom, Zmom, Eden, Eint, Temp;
    static int       FirstAdv,  LastAdv,  NumAdv;
    static int       FirstSpec, LastSpec, NumSpec;

protected:

    static void read_params ();

    CNSReact& getLevel (int lev);

    FluxRegister& getFluxReg ();
    FluxRegister& getFluxReg (int lev);

    void reflux ();

    void enforce_nonnegative_species (MultiFab& S_new);

    void enforce_consistent_e (MultiFab& S);

    void avgDown ();

    void avgDown (int state_indx);

    void buildMetrics ();

    // integrate derived quantities over domain

    Real sumDerive (const std::string& name, Real time);

    virtual void sum_integrated_quantities ();

    //
    // Interpolate cell-centered cync correction from coarse to fine.
    //
    enum SyncInterpType
    {
        PC_T,
        CellCons_T,
        CellConsLin_T,
        CellConsProt_T
    };

    void SyncInterp (MultiFab&      CrseSync,
                     int            c_lev,
                     MultiFab&      FineSync,
                     int            f_lev,
                     IntVect&       ratio,
                     int            src_comp,
                     int            dest_comp,
                     int            num_comp,
                     int            increment,
                     Real           dt_clev,
                     int**          bc_orig_qty,
                     SyncInterpType which_interp = CellCons_T,
                     int            state_comp   = -1);

    //
    // The data.
    //
    MultiFab             volume;
    MultiFab             area[BL_SPACEDIM];
    MultiFab             dLogArea[1];
    Array< Array<Real> > radius;
    FluxRegister*        flux_reg;
    //
    // Static data members.
    //
    static bool      dump_old;
    static int       radius_grow;
    static int       verbose;
    static Real      cfl;
    static Real      init_shrink;
    static Real      change_max;
    static int       do_reflux;
    static ErrorList err_list;
    static BCRec     phys_bc;
    static int       NUM_GROW;

    static int       allow_untagging;
    static int       normalize_species;
    static int       do_special_tagging;
    static int       ppm_type;

    static Real      small_dens;
    static Real      small_temp;
    static Real      small_pres;

    static Real      gravx;
    static Real      gravy;
    static Real      gravz;

    bool             FillPatchedOldState_ok;

    static class ChemDriver* chemSolve;

};    

//
// Inlines.
//

inline
int
CNSReact::numGrow()
{
    return NUM_GROW;
}

inline
MultiFab*
CNSReact::Area ()
{
    return area;
}

inline
MultiFab&
CNSReact::Area (int dir)
{
    return area[dir];
}

inline
MultiFab&
CNSReact::Volume ()
{
    return volume;
}

inline
CNSReact&
CNSReact::getLevel (int lev)
{
    return *(CNSReact *) &parent->getLevel(lev);
}

inline
FluxRegister&
CNSReact::getFluxReg ()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
FluxRegister&
CNSReact::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}

#endif /*_CNSReact_H_*/
