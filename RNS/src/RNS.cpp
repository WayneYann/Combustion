#include <winstd.H>

#ifndef WIN32
#include <unistd.h>
#endif

#include <iomanip>

#include <algorithm>
#include <cstdio>
#include <vector>
#include <iostream>
#include <string>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <Utility.H>
#include <CONSTANTS.H>
#include <Derive_F.H>
#include <VisMF.H>
#include <TagBox.H>
#include <ParmParse.H>

#include <ChemDriver.H>

#include <RNS.H>
#include <RNS_F.H>

#ifndef NDEBUG
static int  sum_int = 1;
#else
static int  sum_int = -1;
#endif
static Real fixed_dt     = -1.0;
static Real initial_dt   = -1.0;
static Real dt_cutoff    = 0.0;

bool         RNS::dump_old      = false;

int          RNS::verbose       = 0;
#if (BL_SPACEDIM == 1) 
Real         RNS::cfl           = 0.8;
Real         RNS::diff_cfl      = 0.8;
#elif (BL_SPACEDIM == 2)
Real         RNS::cfl           = 0.4;
Real         RNS::diff_cfl      = 0.4;
#else
Real         RNS::cfl           = 0.27;
Real         RNS::diff_cfl      = 0.27;
#endif
Real         RNS::init_shrink   = 1.0;
Real         RNS::change_max    = 1.1;
BCRec        RNS::phys_bc;
int          RNS::NUM_STATE     = -1;
int          RNS::do_reflux     = 1;
int          RNS::NUM_GROW      = -1;

int          RNS::Density       = -1;
int          RNS::Xmom          = -1;
int          RNS::Ymom          = -1;
int          RNS::Zmom          = -1;
int          RNS::Eden          = -1;
int          RNS::Temp          = -1;

ChemDriver*  RNS::chemSolve     = 0;
int          RNS::NumSpec       = 1;
int          RNS::FirstSpec     = -1;
int          RNS::LastSpec      = -1;

Real         RNS::small_dens    = -1.e200;
Real         RNS::small_temp    = -1.e200;
Real         RNS::small_pres    = -1.e200;
Real         RNS::gamma         = 1.4;

Real         RNS::gravity       = 0.0;
Real         RNS::Treference    = 298.0; 

int          RNS::RK_order      = 2;

RNS::RiemannType RNS::Riemann   = RNS::HLL;
Real             RNS::difmag    = -1.0;  // for JBB & HLLC Riemann solvers

std::string  RNS::fuelName           = "";
int          RNS::fuelID             = -1;
std::string  RNS::oxidizerName       = "O2";
int          RNS::oxidizerID         = -1;
std::string  RNS::productName        = "";
int          RNS::productID          = -1;
std::string  RNS::flameTracName      = "H";
int          RNS::flameTracID        = -1;

ErrorList    RNS::err_list;
int          RNS::allow_untagging    = 0;
int          RNS::do_density_ref     = 0;
int          RNS::do_temperature_ref = 0;
int          RNS::do_pressure_ref    = 0;
int          RNS::do_velocity_ref    = 0;
int          RNS::do_vorticity_ref   = 0;
int          RNS::do_flametrac_ref   = 0;

int          RNS::plot_cons            = 0;
int          RNS::plot_prim            = 1;
int          RNS::plot_magvel          = 1;
int          RNS::plot_Mach            = 1;
int          RNS::plot_divu            = 1;
int          RNS::plot_magvort         = 1;
int          RNS::plot_X               = 0;
int          RNS::plot_omegadot        = 0;
int          RNS::plot_dYdt            = 0;
int          RNS::plot_heatRelease     = 1;
int          RNS::plot_fuelConsumption = 1;
int          RNS::plot_primplus        = 1;

int          RNS::icomp_cons            = -1;
int          RNS::icomp_prim            = -1; 
int          RNS::icomp_magvel		= -1;
int          RNS::icomp_Mach		= -1;
int          RNS::icomp_divu		= -1;
int          RNS::icomp_magvort		= -1;
int          RNS::icomp_X		= -1;
int          RNS::icomp_omegadot	= -1;
int          RNS::icomp_dYdt		= -1;
int          RNS::icomp_heatRelease	= -1;
int          RNS::icomp_fuelConsumption = -1;
std::vector<std::string> RNS::plot_names;

std::string  RNS::job_name = "";

#ifdef _OPENMP
std::vector<int> RNS::blocksize(BL_SPACEDIM, 8);
#else
std::vector<int> RNS::blocksize(BL_SPACEDIM, 2048);
#endif

int          RNS::do_quartic_interp   = 1;

int          RNS::do_weno             = 1;
int          RNS::do_quadrature_weno  = 0;
int          RNS::do_component_weno   = 0;

int          RNS::do_chemistry        = 1;
int          RNS::use_vode            = 0;
int          RNS::new_J_cell          = 0; // new Jacobian for each cell?
#ifdef USE_SDCLIB
RNS::ChemSolverType RNS::chem_solver  = RNS::BE_BURNING;
#else
RNS::ChemSolverType RNS::chem_solver  = RNS::GAUSS_BURNING;
#endif
int          RNS::f2comp_simple_dUdt  = 1; // set dUdt = \Delta U / \Delta t in f2comp
int          RNS::f2comp_niter_good   = 1; // starting from iter ?, U is considered good.

// this will be reset upon restart
Real         RNS::previousCPUTimeUsed = 0.0;
Real         RNS::startCPUTime = 0.0;

void
RNS::variableCleanUp () 
{
    desc_lst.clear();
    delete chemSolve;
    chemSolve = 0;
}

void
RNS::read_params ()
{
    static bool done = false;

    if (done) return;
    
    done = true;
    
    ParmParse pp("rns");   
    
    pp.query("v",verbose);
    pp.get("init_shrink",init_shrink);
    pp.get("cfl",cfl);
    pp.query("diff_cfl",diff_cfl);
    pp.query("change_max",change_max);
    pp.query("fixed_dt",fixed_dt);
    pp.query("initial_dt",initial_dt);
    pp.query("sum_int",sum_int);
    pp.query("do_reflux",do_reflux);
    do_reflux = (do_reflux ? 1 : 0);
    pp.get("dt_cutoff",dt_cutoff);
    
    pp.query("dump_old",dump_old);
    
    pp.query("small_dens",small_dens);
    pp.query("small_temp",small_temp);
    pp.query("small_pres",small_pres);
    pp.query("gamma",gamma);

    pp.query("gravity", gravity);
    pp.query("Treference",Treference);
    
    pp.query("RK_order",RK_order);

    {
	int riemann;
	if (pp.query("Riemann", riemann)) {
	    Riemann = static_cast<RiemannType>(riemann);
	}
    }

    // some Riemann solvers need artificial viscosity to suppress odd-even decouping
    if (Riemann == JBB)
    {
	difmag = 0.1;
    }
    else if (Riemann == HLLC)
    {
	difmag = 0.1;
    }
    else
    {
	difmag = -1.0;
    }
    pp.query("difmag", difmag);

    // Get boundary conditions
    Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	phys_bc.setLo(i,lo_bc[i]);
	phys_bc.setHi(i,hi_bc[i]);
    }

    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (Geometry::isAnyPeriodic())
    {
	//
	// Do idiot check.  Periodic means interior in those directions.
	//
	for (int dir = 0; dir<BL_SPACEDIM; dir++)
        {
	    if (Geometry::isPeriodic(dir))
            {
		if (lo_bc[dir] != Interior)
                {
		    std::cerr << "RNS::read_params:periodic in direction "
			      << dir
			      << " but low BC is not Interior\n";
		    BoxLib::Error();
                }
		if (hi_bc[dir] != Interior)
                {
		    std::cerr << "RNS::read_params:periodic in direction "
			      << dir
			      << " but high BC is not Interior\n";
		    BoxLib::Error();
                }
            }
        }
    }
    else
    {
	//
	// Do idiot check.  If not periodic, should be no interior.
	//
	for (int dir=0; dir<BL_SPACEDIM; dir++)
        {
	    if (lo_bc[dir] == Interior)
            {
		std::cerr << "RNS::read_params:interior bc in direction "
			  << dir
			  << " but not periodic\n";
		BoxLib::Error();
            }
	    if (hi_bc[dir] == Interior)
            {
		std::cerr << "RNS::read_params:interior bc in direction "
			  << dir
			  << " but not periodic\n";
		BoxLib::Error();
            }
        }
    }

    pp.query("fuelName"     , fuelName);
    pp.query("oxidizerName" , oxidizerName);
    pp.query("productName"  , productName);
    pp.query("flameTracName", flameTracName);
    
    pp.query("allow_untagging"   , allow_untagging);
    pp.query("do_density_ref"    , do_density_ref);
    pp.query("do_temperature_ref", do_temperature_ref);
    pp.query("do_pressure_ref"   , do_pressure_ref);
    pp.query("do_velocity_ref"   , do_velocity_ref);
    pp.query("do_vorticity_ref"  , do_vorticity_ref);
    pp.query("do_flametrac_ref"  , do_flametrac_ref);

    pp.query("plot_cons"           , plot_cons           );
    pp.query("plot_prim"           , plot_prim           );
    pp.query("plot_magvel"         , plot_magvel         );
    pp.query("plot_Mach"           , plot_Mach           );
    pp.query("plot_divu"           , plot_divu           );
    pp.query("plot_magvort"        , plot_magvort        );
    pp.query("plot_X"              , plot_X              );
    pp.query("plot_omegadot"       , plot_omegadot       );
    pp.query("plot_dYdt"           , plot_dYdt           );
    pp.query("plot_heatRelease"    , plot_heatRelease    );
    pp.query("plot_fuelConsumption", plot_fuelConsumption);
    pp.query("plot_primplus"       , plot_primplus);

    pp.query("job_name",job_name);  

    pp.queryarr("blocksize", blocksize);

    pp.query("do_quartic_interp", do_quartic_interp);

    pp.query("do_weno", do_weno);
    pp.query("do_quadrature_weno", do_quadrature_weno);
    pp.query("do_component_weno", do_component_weno);
    if (ChemDriver::isNull())
    {
	do_weno = 1;  // must do weno
    }

    pp.query("do_chemistry", do_chemistry);
    pp.query("use_vode", use_vode);
    pp.query("new_J_cell", new_J_cell);
    {
	int chem_solver_i;
	if (pp.query("chem_solver", chem_solver_i)) {
	    chem_solver = static_cast<ChemSolverType>(chem_solver_i);
	}
    }
    if (chem_solver == RNS::BE_BURNING) {
	f2comp_simple_dUdt = 1;
	f2comp_niter_good = 1;
    }
    else {
	f2comp_simple_dUdt = 0;
	f2comp_niter_good = 100000;
    }
    pp.query("f2comp_simple_dUdt", f2comp_simple_dUdt);
    pp.query("f2comp_niter_good" , f2comp_niter_good);
}

RNS::RNS ()
{
    flux_reg = 0;
}

RNS::RNS (Amr&            papa,
	  int             lev,
	  const Geometry& level_geom,
	  const BoxArray& bl,
	  Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,time) 
{
    buildMetrics();
    
    flux_reg = 0;
    if (level > 0 && do_reflux) 
    {
	flux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
    }
}

RNS::~RNS () 
{
    delete flux_reg;
}

void
RNS::buildMetrics ()
{
    if ( Geometry::IsSPHERICAL() || Geometry::IsRZ() ) 
	BoxLib::Abort("We don't support curvilinear coordinate systems.");

    //
    // Build volume and face area arrays.
    // volume is not PArrayManaged, must manually delete.
    //
    volume.clear();
    //
    // area is not PArrayManaged, must manually delete.
    //
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	area[dir].clear();
    }
    
    geom.GetVolume(volume,grids,0);
    
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	geom.GetFaceArea(area[dir],grids,dir,0);
    }
}

void
RNS::initData ()
{
    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    int ns          = NUM_STATE;
    const Real* dx  = geom.CellSize();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time   = state[State_Type].curTime();
    
    // make sure dx = dy = dz -- that's all we guarantee to support
    const Real SMALL = 1.e-13;
#if (BL_SPACEDIM == 2)
    if (fabs(dx[0] - dx[1]) > SMALL*dx[0])
    {
	BoxLib::Abort("We don't support dx != dy");
    }
#elif (BL_SPACEDIM == 3)
    if ( (fabs(dx[0] - dx[1]) > SMALL*dx[0]) || (fabs(dx[0] - dx[2]) > SMALL*dx[0]) )
    {
	BoxLib::Abort("We don't support dx != dy != dz");
    }
#endif
    
    if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "Initializing the data at level " << level << std::endl;

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
	RealBox    gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
	const Box& box     = mfi.validbox();
	const int* lo      = box.loVect();
	const int* hi      = box.hiVect();
	
	BL_FORT_PROC_CALL(RNS_INITDATA,rns_initdata)
	    (level, cur_time, lo, hi, ns,
	     BL_TO_FORTRAN(S_new[mfi]), dx,
	     gridloc.lo(), gridloc.hi());
    }
    
    if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "Done initializing the level " << level << " data " << std::endl;
}

void
RNS::init (AmrLevel &old)
{
    RNS* oldlev = (RNS*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[State_Type].curTime();
    Real prev_time = oldlev->state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);
    
    MultiFab& S_new = get_new_data(State_Type);
    
    for (FillPatchIterator fpi(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
	 fpi.isValid();
	 ++fpi)
    {
	S_new[fpi].copy(fpi());
    }
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
RNS::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();
    
    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);
    
    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NUM_STATE);
}

Real
RNS::initialTimeStep ()
{
    Real dummy_dt = 0.0;
    return (initial_dt > 0.0) ? initial_dt : init_shrink*estTimeStep(dummy_dt);
}

Real
RNS::estTimeStep (Real dt_old)
{
    if (fixed_dt > 0.0)
	return fixed_dt;
    
    // This is just a dummy value to start with 
    Real estdt  = 1.0e+20;
    
    const Real* dx    = geom.CellSize();
    const MultiFab& stateMF = get_new_data(State_Type);
    
    for (MFIter mfi(stateMF); mfi.isValid(); ++mfi)
    {
	const Box& box = mfi.validbox();
	Real dt = estdt;
	BL_FORT_PROC_CALL(RNS_ESTDT,rns_estdt)
	    (BL_TO_FORTRAN(stateMF[mfi]),
	     box.loVect(),box.hiVect(),dx,&dt);
	
	estdt = std::min(estdt,dt);
    }
    estdt *= cfl;

#ifdef NULLCHEMISTRY

    ParallelDescriptor::ReduceRealMin(estdt);

    if (verbose && ParallelDescriptor::IOProcessor())
	cout << "At level " << level << ":  estdt = " << estdt << std::endl;

    return estdt;
    
#else

    Real estdt_diff = 1.e+20;
    for (MFIter mfi(stateMF); mfi.isValid(); ++mfi)
    {
	const Box& box = mfi.validbox();
	Real dt = estdt_diff;
	BL_FORT_PROC_CALL(RNS_ESTDT_DIFF,rns_estdt_diff)
	    (BL_TO_FORTRAN(stateMF[mfi]),
	     box.loVect(),box.hiVect(),dx,&dt);
	
	estdt_diff = std::min(estdt_diff,dt);
    }
    estdt_diff *= diff_cfl;

    Real estdts[2];
    estdts[0] = estdt;
    estdts[1] = estdt_diff;

    ParallelDescriptor::ReduceRealMin(estdts, 2);

    Real estdt_min = std::min(estdts[0], estdts[1]);

    if (verbose && ParallelDescriptor::IOProcessor()){
	cout << "At level " << level << ":  estdt = " 
	     << estdts[0] << " (advection),  " << estdts[1] << " (diffusion)" << std::endl;
    }

    return estdt_min;
    
#endif
}

void
RNS::computeNewDt (int                   finest_level,
		   int                   sub_cycle,
		   Array<int>&           n_cycle,
		   const Array<IntVect>& ref_ratio,
		   Array<Real>&          dt_min,
		   Array<Real>&          dt_level,
		   Real                  stop_time,
		   int                   post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
	return;
    
    int i;
    
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
	RNS& adv_level = getLevel(i);
	dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }
    
    if (fixed_dt <= 0.0)
    {
	if (post_regrid_flag == 1) 
	{
	    //
	    // Limit dt's by pre-regrid dt
	    //
	    for (i = 0; i <= finest_level; i++)
	    {
		dt_min[i] = std::min(dt_min[i],dt_level[i]);
	    }
	} 
	else 
	{
	    //
	    // Limit dt's by change_max * old dt
	    //
	    for (i = 0; i <= finest_level; i++)
	    {
		if (verbose && ParallelDescriptor::IOProcessor())
		{
		    if (dt_min[i] > change_max*dt_level[i])
		    {
			cout << "RNS::computeNewDt : limiting dt at level " << i << std::endl;
			cout << " ... new dt computed: " << dt_min[i] << std::endl;
			cout << " ... but limiting to: " << change_max*dt_level[i] <<
			    " = " << change_max << " * " << dt_level[i] << std::endl;
		    }
		}
		dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
	    }
	} 
    }
    
    //
    // Find the minimum over all levels
    //
    for (i = 0; i <= finest_level; i++)
    {
	n_factor *= n_cycle[i];
	dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }
    
    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) 
    {
	if ((cur_time + dt_0) > (stop_time - eps))
	{
	    dt_0 = stop_time - cur_time;
	}
    }
    
    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
	n_factor *= n_cycle[i];
	dt_level[i] = dt_0/n_factor;
    }
}

void
RNS::computeInitialDt (int                   finest_level,
		       int                   sub_cycle,
		       Array<int>&           n_cycle,
		       const Array<IntVect>& ref_ratio,
		       Array<Real>&          dt_level,
		       Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
	return;
    
    int i;
    
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
	dt_level[i] = getLevel(i).initialTimeStep();
	n_factor   *= n_cycle[i];
	dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }
    
    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) 
    {
	if ((cur_time + dt_0) > (stop_time - eps))
	{
	    dt_0 = stop_time - cur_time;
	}
    }
    
    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
	n_factor *= n_cycle[i];
	dt_level[i] = dt_0/n_factor;
    }
}

void
RNS::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();
    
    if (do_reflux && level < finest_level) 
    {
	reflux();

	avgDown();

	MultiFab& S_new_crse = get_new_data(State_Type);
	post_update(S_new_crse);
    }
    else if (level < finest_level) 
    {
	avgDown();
    }
}

void
RNS::post_restart ()
{
}

void
RNS::postCoarseTimeStep (Real cumtime)
{
    if (sum_int <= 0) return;

    int nstep = parent->levelSteps(0);
    if (nstep % sum_int == 0) {
	sum_conserved_variables();
    }
}

void
RNS::post_regrid (int lbase,
		  int new_finest)
{
}

void
RNS::post_init (Real stop_time)
{
    if (level > 0)
	return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--) 
    {
	getLevel(k).avgDown();
    }
}

int
RNS::okToContinue ()
{
    if (level > 0)
	return 1;
    
    int test = (parent->dtLevel(0) < dt_cutoff) ? 0 : 1;
    return test;
}

void
RNS::reflux ()
{
    BL_ASSERT(level<parent->finestLevel());
  
    const Real strt = ParallelDescriptor::second();
    
    getFluxReg(level+1).Reflux(get_new_data(State_Type),volume,1.0,0,0,NUM_STATE,geom);
    
    if (verbose)
    {
	const int IOProc = ParallelDescriptor::IOProcessorNumber();
	Real      end    = ParallelDescriptor::second() - strt;
	
	ParallelDescriptor::ReduceRealMax(end,IOProc);
      
	if (ParallelDescriptor::IOProcessor())
	    std::cout << "RNS::reflux() at level " << level << " : time = " << end << std::endl;
    }
}

void
RNS::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(State_Type);
}

void
RNS::avgDown (MultiFab& S_crse, const MultiFab& S_fine)
{
    if (level == parent->finestLevel()) return;
    
    RNS&       fine_lev = getLevel(level+1);
    MultiFab&  fvolume  = fine_lev.volume;
    const int  ncomp    = S_fine.nComp();
    
    BL_ASSERT(S_crse.boxArray() == volume.boxArray());
    BL_ASSERT(fvolume.boxArray() == S_fine.boxArray());
    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_S_fine_BA(S_fine.boxArray().size());
    
    for (int i = 0; i < S_fine.boxArray().size(); ++i)
    {
	crse_S_fine_BA.set(i,BoxLib::coarsen(S_fine.boxArray()[i],fine_ratio));
    }
    
    MultiFab crse_S_fine(crse_S_fine_BA,ncomp,0);
    MultiFab crse_fvolume(crse_S_fine_BA,1,0);
    
    crse_fvolume.copy(volume);
    
    for (MFIter mfi(S_fine); mfi.isValid(); ++mfi)
    {
	const int        i        = mfi.index();
	const Box&       ovlp     = crse_S_fine_BA[i];
	FArrayBox&       crse_fab = crse_S_fine[i];
	const FArrayBox& crse_vol = crse_fvolume[i];
	const FArrayBox& fine_fab = S_fine[i];
	const FArrayBox& fine_vol = fvolume[i];
	
	BL_FORT_PROC_CALL(RNS_AVGDOWN,rns_avgdown)
	    (BL_TO_FORTRAN(crse_fab), ncomp,
	     BL_TO_FORTRAN(crse_vol),
	     BL_TO_FORTRAN(fine_fab),
	     BL_TO_FORTRAN(fine_vol),
	     ovlp.loVect(),ovlp.hiVect(),fine_ratio.getVect());
    }
    
    S_crse.copy(crse_S_fine);
}


void
RNS::avgDown (int state_indx)
{
    RNS&       fine_lev = getLevel(level+1);
    MultiFab&  S_crse   = get_new_data(state_indx);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    avgDown(S_crse, S_fine);
}

void
RNS::errorEst (TagBoxArray& tags,
	       int          clearval,
	       int          tagval,
	       Real         time,
	       int          n_error_buf,
	       int          ngrow)
{
    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();
    Array<int>  itags;
    
    for (int j = 0; j < err_list.size(); j++)
    {
	MultiFab* mf = derive(err_list[j].name(), time, err_list[j].nGrow());
	
	BL_ASSERT(!(mf == 0));
	
	for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
        {
	    int         idx     = mfi.index();
	    RealBox     gridloc = RealBox(grids[idx],geom.CellSize(),geom.ProbLo());
	    itags               = tags[idx].tags();
	    int*        tptr    = itags.dataPtr();
	    const int*  tlo     = tags[idx].box().loVect();
	    const int*  thi     = tags[idx].box().hiVect();
	    const int*  lo      = mfi.validbox().loVect();
	    const int*  hi      = mfi.validbox().hiVect();
	    const Real* xlo     = gridloc.lo();
	    Real*       dat     = (*mf)[mfi].dataPtr();
	    const int*  dlo     = (*mf)[mfi].box().loVect();
	    const int*  dhi     = (*mf)[mfi].box().hiVect();
	    const int   ncomp   = (*mf)[mfi].nComp();
	    
	    err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
				  &clearval, dat, ARLIM(dlo), ARLIM(dhi),
				  lo,hi, &ncomp, domain_lo, domain_hi,
				  dx, xlo, prob_lo, &time, &level);
	    //
	    // Don't forget to set the tags in the TagBox.
	    //
	    if (allow_untagging == 1) 
            {
		tags[idx].tags_and_untags(itags);
            } 
	    else 
	    {
		tags[idx].tags(itags);
	    }
        }
	
	delete mf;
    }
}

MultiFab*
RNS::derive (const std::string& name,
	     Real           time,
	     int            ngrow)
{
    return AmrLevel::derive(name,time,ngrow);
}

void
RNS::derive (const std::string& name,
	     Real           time,
	     MultiFab&      mf,
	     int            dcomp)
{
    AmrLevel::derive(name,time,mf,dcomp);
}

Real
RNS::getCPUTime()
{
    int numCores = ParallelDescriptor::NProcs();
#ifdef _OPENMP
    numCores *= omp_get_max_threads();
#endif    

    Real T = numCores*(ParallelDescriptor::second() - startCPUTime) + 
	previousCPUTimeUsed;
    
    return T;
}
